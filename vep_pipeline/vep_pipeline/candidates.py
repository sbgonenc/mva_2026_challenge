"""Read the merged VEP VCF and turn genotype calls into findings.

The VCF is streamed once: `bcftools view -s <sample> -i 'GT="alt"' [-f PASS,.] | bcftools query` with the CSQ string
per variant. Only variants with usable evidence (under the full set and the evidence-removed sets below, which cover
the sensitivity scenarios and any configuration that only removes evidence) and popmax <= af.max_recessive are kept. Expect tens of minutes on a WGS VCF; run it
through `prioritise --slurm`, never on the login node.
"""
from __future__ import annotations
import itertools
import re
import shlex
import subprocess
import tempfile
from dataclasses import dataclass, replace

from .checks import X_NONPAR
from .config import ConfigError
from .evidence import (LEVELS, MISSENSE_TOOLS, PANEL_CLASSES, PANEL_REASON, Evidence, gene_evidence, lower,
                       variant_level)
from .jobs import JobError, _key, _tool

FORMAT_TAGS = ("GT", "AD", "DP", "GQ", "PS", "PGT", "PID")
_CSQ = re.compile(r'##INFO=<ID=CSQ,.*Format:\s*([^"]+)"')


@dataclass(frozen=True)
class Call:
    chrom: str
    pos: int
    ref: str
    alt: str
    filt: str
    orig: str                 # INFO/EXO_ORIG or '.'
    gt: str
    ad: str
    dp: str
    gq: str
    ps: str
    pgt: str
    pid: str

    @property
    def key(self) -> str:
        return _key(self.chrom, self.pos, self.ref, self.alt)


@dataclass(frozen=True)
class Candidate:
    call: Call
    ev: Evidence


def header_fields(vcf) -> tuple:
    """(CSQ field names, INFO ids, FORMAT ids) from the VCF header."""
    text = _tool(["bcftools", "view", "-h", str(vcf)])
    m = _CSQ.search(text)
    if not m:
        raise ConfigError(f"{vcf}: no CSQ header (is this the VEP output?)")
    info = set(re.findall(r"##INFO=<ID=([^,>]+)", text))
    fmt = set(re.findall(r"##FORMAT=<ID=([^,>]+)", text))
    return tuple(m.group(1).strip().split("|")), info, fmt


def query_format(info_ids, format_ids) -> str:
    """bcftools query -f string; absent INFO/FORMAT tags become a literal '.' so bcftools never errors."""
    orig = "%INFO/EXO_ORIG" if "EXO_ORIG" in info_ids else "."
    sample = "\t".join(f"%{t}" if t in format_ids else "." for t in FORMAT_TAGS)
    return f"%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t{orig}[\t{sample}]\t%INFO/CSQ\n"


def stream_cmd(vcf, sample, fmt, pass_only) -> str:
    q = shlex.quote
    keep = " -f PASS,." if pass_only else ""
    return (f"bcftools view -s {q(sample)} -i 'GT=\"alt\"'{keep} -Ou {q(str(vcf))} | "
            f"bcftools query -f {q(fmt)} -")


def _lines(cmd: str):
    with tempfile.TemporaryFile("w+") as err:
        p = subprocess.Popen(["bash", "-o", "pipefail", "-c", cmd], stdout=subprocess.PIPE, stderr=err, text=True)
        for line in p.stdout:
            yield line.rstrip("\n")
        rc = p.wait()
        if rc:
            err.seek(0)
            raise JobError(f"VCF stream failed (exit {rc}): {err.read().strip()[-300:]}")


def parse_line(line: str, csq_fields) -> tuple:
    f = line.split("\t")
    call = Call(f[0], int(f[1]), f[2], f[3], f[4], f[5], *f[6:6 + len(FORMAT_TAGS)])
    csq = f[6 + len(FORMAT_TAGS)]
    entries = [{k: v for k, v in zip(csq_fields, e.split("|")) if v not in ("", ".", "-")}
               for e in csq.split(",") if e and e != "."]
    return call, entries


def _bare(chrom: str) -> str:
    return chrom[3:] if chrom.startswith("chr") else chrom


def excluded_region(call: Call, regions) -> str | None:
    """Name of the first configured region containing the call, else None (contig names with or without 'chr')."""
    c = _bare(call.chrom)
    return next((r["name"] for r in regions
                 if _bare(r["chrom"]) == c and r["start"] <= call.pos <= r["end"]), None)


def deep_panel_genes(settings) -> frozenset:
    """Symbols whose variants get the panel-candidate search: the panel genes when deep_second_hit is on, else none."""
    gp = settings["gene_panel"]
    return frozenset(gp["genes"]) if gp["deep_second_hit"] else frozenset()


def collect(vcf, sample, settings, genes, remm=None, excluded=None) -> list:
    """Candidates (one per variant and gene) worth keeping in memory; see the module docstring. `remm` maps a
    variant key to its Genomiser ReMM score; it is attached before the evidence pre-filter so ReMM-only variants stay.
    Records inside settings["exclude_regions"] are skipped; when `excluded` is a list, (key, region) is appended
    for each."""
    csq_fields, info, fmt = header_fields(vcf)
    for need in ("SYMBOL", "Consequence"):
        if need not in csq_fields:
            raise ConfigError(f"{vcf}: CSQ has no {need} field")
    pops = tuple(settings["af"]["populations"])
    gate = None if "population" in settings["exclude_evidence_classes"] else settings["af"]["max_recessive"]
    # Every evidence set x both missense tools: a scenario or setting that removes evidence or switches the tool must
    # never need a variant this pre-filter dropped.
    evidence_sets = (frozenset(), frozenset({"known_variant"}),
                     frozenset({"known_variant", "clinically_trained_predictor"}))
    # A SpliceAI score <= 0.1 (BP4) vetoes a splice-region panel candidate; without unsupervised_predictor the term
    # alone counts, so panel genes are also tested with that class removed.
    panel_sets = evidence_sets + tuple(x | {"unsupervised_predictor"} for x in evidence_sets)
    panel = deep_panel_genes(settings)
    out = []
    cmd = stream_cmd(vcf, sample, query_format(info, fmt), settings["qc"]["pass_only"])
    for line in _lines(cmd):
        call, entries = parse_line(line, csq_fields)
        if call.chrom in ("MT", "chrM"):
            continue                                       # heteroplasmy is not assessed here
        region = excluded_region(call, settings["exclude_regions"])
        if region:
            if excluded is not None:
                excluded.append((call.key, region))
            continue
        for ev in gene_evidence(entries, pops).values():
            if remm:
                ev = replace(ev, remm=remm.get(call.key))
            if gate is not None and ev.popmax > gate:
                continue
            deep = ev.symbol in panel
            if any(variant_level(ev, genes.get(ev.symbol), x, tool, deep)[0]
                   for x in (panel_sets if deep else evidence_sets) for tool in MISSENSE_TOOLS):
                out.append(Candidate(call, ev))
    return out


# ----------------------------------------------------------------------------- genotypes, alleles, findings

def _alleles_of(gt: str) -> list:
    return re.split(r"[/|]", gt) if gt not in ("", ".") else []


def zygosity(call: Call, sex: str) -> str | None:
    a = _alleles_of(call.gt)
    if not a or "." in a or all(x == "0" for x in a):
        return None
    if len(a) == 1:
        return "hemi"
    hom = all(x != "0" for x in a)
    nonpar_x = call.chrom in ("X", "chrX") and X_NONPAR[0] <= call.pos <= X_NONPAR[1]
    if hom and nonpar_x and sex == "MALE":
        return "hemi"
    return "hom" if hom else "het"


def _int(v):
    try:
        return int(v)
    except ValueError:
        return None


def vaf(call: Call) -> float | None:
    parts = [_int(x) for x in call.ad.split(",")[:2]]
    if len(parts) < 2 or None in parts or sum(parts) == 0:
        return None
    return parts[1] / sum(parts)


def qc_ok(call: Call, qc: dict) -> bool:
    if qc["pass_only"] and call.filt not in ("PASS", "."):
        return False
    gq, dp = _int(call.gq), _int(call.dp)
    return (gq is None or gq >= qc["min_gq"]) and (dp is None or dp >= qc["min_dp"])


def orig_record(call: Call):
    """(chrom, pos, ref, alts, alt_index) from EXO_ORIG='CHROM|POS|REF|ALT|IDX' (prepare-patient-vcf.sh sets it on
    every record, changed or not; Task 1); None without the tag."""
    if call.orig in ("", "."):
        return None
    parts = call.orig.split("|")
    if len(parts) != 5:
        raise ConfigError(f"{call.key}: unexpected EXO_ORIG {call.orig!r} (expected CHROM|POS|REF|ALT|IDX)")
    chrom, pos, ref, alts, idx = parts
    return chrom, int(pos), ref, alts, int(idx)


@dataclass(frozen=True)
class Allele:
    pieces: tuple             # Candidates of one original allele (an atomised MNV has several), position order
    level: str
    second_hit_only: bool
    reasons: tuple
    zyg: str
    low_vaf: bool

    @property
    def first(self) -> Candidate:
        return self.pieces[0]

    @property
    def panel_only(self) -> bool:
        """True when the only evidence is a gene-panel candidate (never a single finding, never a pair of two)."""
        return bool(self.reasons) and all(r.startswith(PANEL_REASON) for r in self.reasons)

    @property
    def panel_rank(self) -> int:
        """Index in PANEL_CLASSES of a panel-only allele (lower is a better candidate), else -1."""
        if not self.panel_only:
            return -1
        return min(PANEL_CLASSES.index(r[len(PANEL_REASON) + 2:-1]) for r in self.reasons)

    @property
    def popmax(self) -> float:
        return max(p.ev.popmax for p in self.pieces)

    @property
    def orig(self):
        return orig_record(self.first.call)

    @property
    def key(self) -> str:
        return self.first.call.key


def build_alleles(cands, genes, settings, sex) -> dict:
    """{symbol: [Allele]} after QC. Atomised pieces of one original allele (same EXO_ORIG record and ALT index) are
    merged: they sit on one haplotype and must never form a compound-het pair."""
    exclude = frozenset(settings["exclude_evidence_classes"])
    panel = deep_panel_genes(settings)
    groups = {}
    for c in cands:
        if qc_ok(c.call, settings["qc"]):
            o = orig_record(c.call)
            groups.setdefault((c.ev.symbol, o if o else c.call.key), []).append(c)
    out = {}
    for (symbol, _), pieces in sorted(groups.items(), key=lambda kv: (kv[1][0].call.chrom, kv[1][0].call.pos,
                                                                     kv[1][0].call.alt, kv[0][0])):
        pieces.sort(key=lambda p: p.call.pos)
        zyg = zygosity(pieces[0].call, sex)
        levels = [variant_level(p.ev, genes.get(symbol), exclude, settings["missense_tool"],
                               symbol in panel) for p in pieces]
        usable = [x for x in levels if x[0]]
        if zyg is None or not usable:
            continue
        level = max((x[0] for x in usable), key=LEVELS.index)
        v = vaf(pieces[0].call)
        low = zyg == "het" and v is not None and v < settings["qc"]["low_vaf"]
        out.setdefault(symbol, []).append(Allele(
            tuple(pieces), lower(level) if low else level, all(x[1] for x in usable),
            tuple(sorted({r for x in usable for r in x[2]})), zyg, low))
    return out


def _hap(gt: str):
    """Haplotype index (0/1) carrying the ALT in a phased het GT, else None."""
    a = gt.split("|")
    if len(a) != 2 or "|" not in gt or a[0] == a[1]:
        return None
    return 0 if a[0] != "0" else 1


PHASE_SOURCES = ("multiallelic", "ps", "pid", "whatshap")


def _decide(source, a: Allele, b: Allele, readback):
    """'cis' / 'trans' when one phase source decides the pair, else None."""
    ca, cb = a.first.call, b.first.call
    if source == "multiallelic":
        oa, ob = a.orig, b.orig
        return "trans" if oa and ob and oa[:4] == ob[:4] and oa[4] != ob[4] else None   # ALTs of one record
    if source == "whatshap":
        if not readback or a.key not in readback or b.key not in readback:
            return None
        (ga, sa), (gb, sb) = readback[a.key], readback[b.key]
    elif source == "ps":
        ga, gb, sa, sb = ca.gt, cb.gt, ca.ps, cb.ps
    elif source == "pid":
        ga, gb, sa, sb = ca.pgt, cb.pgt, ca.pid, cb.pid
    else:
        return None                                      # unknown source name: no edge
    if sa in ("", ".") or sa != sb:
        return None
    ha, hb = _hap(ga), _hap(gb)
    return None if ha is None or hb is None else "cis" if ha == hb else "trans"


class PhaseGraph:
    """Union-find with parity over the het alleles of one gene: an edge says two alleles are on the same (cis) or
    opposite (trans) haplotypes. Relations follow transitively; a component with a contradicting edge is conflicted."""

    def __init__(self, alleles, readback=None, sources=PHASE_SOURCES):
        self._ix = {id(a): i for i, a in enumerate(alleles)}
        n = len(alleles)
        self._up = list(range(n))                       # parent; parity[i] = parity of i relative to its parent
        self._par = [0] * n
        self._direct = {}                               # (i, j) -> sources of every edge between i and j
        tree = []                                       # (i, j) spanning-tree edges
        bad = []
        for (i, a), (j, b) in itertools.combinations(enumerate(alleles), 2):
            for source in sources:
                state = _decide(source, a, b, readback)
                if state is None:
                    continue
                want = int(state == "trans")
                self._direct.setdefault((i, j), set()).add(source)
                (ri, pi), (rj, pj) = self._find(i), self._find(j)
                if ri != rj:
                    self._up[ri], self._par[ri] = rj, pi ^ pj ^ want
                    tree.append((i, j))
                elif pi ^ pj != want:
                    bad.append(i)
        self._members = list(alleles)
        self._conflicted = {self._find(i)[0] for i in bad}
        self._adj = {}
        for i, j in tree:
            srcs = self._direct[(i, j)]
            self._adj.setdefault(i, []).append((j, srcs))
            self._adj.setdefault(j, []).append((i, srcs))

    def _find(self, i):
        par = 0
        while self._up[i] != i:
            par ^= self._par[i]
            i = self._up[i]
        return i, par

    def _node(self, a):
        return self._ix.get(id(a))

    def _conflict(self, i):
        return self._find(i)[0] in self._conflicted

    def _path_sources(self, i, j) -> tuple:
        seen, stack = {i: ()}, [i]
        while stack:
            u = stack.pop()
            if u == j:
                return tuple(sorted(set(seen[u])))
            for v, srcs in self._adj.get(u, ()):
                if v not in seen:
                    seen[v] = seen[u] + tuple(srcs)
                    stack.append(v)
        return ()

    def relation(self, a: Allele, b: Allele) -> tuple:
        """(state, sources): 'cis' / 'trans' / 'unknown' and the phase sources on the path between a and b."""
        i, j = self._node(a), self._node(b)
        if i is None or j is None:
            return "unknown", ()
        (ri, pi), (rj, pj) = self._find(i), self._find(j)
        if ri != rj:
            return "unknown", ()
        if ri in self._conflicted:
            return "unknown", ("conflict",)
        direct = self._direct.get((min(i, j), max(i, j)))
        return ("trans" if pi ^ pj else "cis"), tuple(sorted(direct)) if direct else self._path_sources(i, j)

    def sides(self, a: Allele) -> tuple:
        """(cis mates, trans mates) of a within its component; none when it is unlinked or conflicted."""
        i = self._node(a)
        if i is None or self._conflict(i):
            return (), ()
        root, par = self._find(i)
        cis, trans = [], []
        for j, m in enumerate(self._members):
            if j != i and self._find(j)[0] == root:
                (trans if self._find(j)[1] ^ par else cis).append(m)
        return tuple(cis), tuple(trans)


def phase_graph(alleles, readback=None, sources=PHASE_SOURCES) -> PhaseGraph:
    """Phase graph over het alleles; `readback` maps an allele key to (gt, ps) from read-backed phasing."""
    return PhaseGraph(list(alleles), readback, sources)


def phase(a: Allele, b: Allele) -> str:
    """'trans', 'cis' or 'unknown' for two het alleles of one gene."""
    return phase_graph([a, b]).relation(a, b)[0]


@dataclass(frozen=True)
class PairHyp:
    a: Allele
    b: Allele
    phase: str                # 'trans' | 'unknown'
    sources: tuple            # phase sources behind it, or ('conflict',)
    companions: tuple = ()    # cis mates folded into a and b


def _rep_key(a: Allele):
    return -LEVELS.index(a.level), a.second_hit_only, a.popmax, a.first.call.pos, a.first.call.alt


def all_pairs(hets: list, graph: PhaseGraph) -> list:
    """Every pair hypothesis, best first (weaker level, stronger level, proven trans, better panel class, rarity). Alleles proven cis
    form one haplotype, represented by its strongest member and paired as one; never cis, never one variant twice,
    at least one allele not second-hit-only."""
    haps = {}
    for a in hets:
        cis, _ = graph.sides(a)
        mates = [m for m in hets if m is a or any(m is c for c in cis)]
        haps.setdefault(frozenset(id(m) for m in mates), mates)
    reps = [(min(m, key=_rep_key), m) for m in haps.values()]
    reps.sort(key=lambda r: _rep_key(r[0])[3:])
    out = []
    for (a, ma), (b, mb) in itertools.combinations(reps, 2):
        if a.key == b.key or (a.second_hit_only and b.second_hit_only):
            continue
        ph, srcs = graph.relation(a, b)
        if ph == "cis":
            continue
        comp = tuple(m for m in ma + mb if m is not a and m is not b)
        lo, hi = sorted((LEVELS.index(a.level), LEVELS.index(b.level)))
        out.append(((lo, hi, ph == "trans", -max(a.panel_rank, b.panel_rank), -(a.popmax + b.popmax)), PairHyp(a, b, ph, srcs, comp)))
    out.sort(key=lambda x: x[0], reverse=True)
    return [h for _, h in out]


def best_pair(hets: list):
    """((a, b), phase) for the best non-cis pair (see all_pairs), or None if no pair qualifies."""
    pairs = all_pairs(hets, phase_graph(hets))
    return ((pairs[0].a, pairs[0].b), pairs[0].phase) if pairs else None


@dataclass(frozen=True)
class Finding:
    gene: str
    alleles: tuple            # 1 Allele, or 2 for a compound-het pair (position order)
    hypothesis: str           # 'dominant' | 'recessive'
    fit: str                  # 'full' | 'single_hit' | 'unknown_moi' | 'panel_partner'
    phase: str                # 'n/a' for singles; 'trans' | 'unknown' for pairs (cis pairs are never built)
    level: str                # weakest allele level
    phase_source: str = ""    # pairs: phase sources joined with '+', or 'conflict'
    companions: tuple = ()    # pairs: cis mates folded into the two alleles
    pair_rank: int = 0        # pairs: index in the gene's all_pairs order (0 = best); singles: 0

    @property
    def distance(self):
        if len(self.alleles) < 2:
            return None
        return abs(self.alleles[1].first.call.pos - self.alleles[0].first.call.pos)


def pairs_tried(info) -> bool:
    """True when compound-het pairs are tried in a gene: recessive (AR/XR) or unknown MOI."""
    modes = info.modes if info else frozenset()
    return not modes or bool(modes & {"AR", "XR"})


def findings_for_gene(symbol, alleles, info, settings, readback=None) -> list:
    """Single, homozygous/hemizygous and pair findings for one gene. A het alone is a dominant hypothesis (dominant
    AF cut-off); hom/hemi and pairs are recessive hypotheses (recessive cut-off). Unknown MOI tries both."""
    af = settings["af"]
    gate_off = "population" in settings["exclude_evidence_classes"]
    dom = lambda a: gate_off or a.popmax <= af["max_dominant"]
    rec = lambda a: gate_off or a.popmax <= af["max_recessive"]
    modes = info.modes if info else frozenset()
    if modes == {"MT"}:
        return []
    out = []
    for a in alleles:
        if a.second_hit_only:
            continue
        if a.zyg in ("hom", "hemi") and rec(a):
            hyp = "recessive" if (modes & {"AR", "XR"} or not modes) else "dominant"
            out.append(Finding(symbol, (a,), hyp, "full" if modes else "unknown_moi", "n/a", a.level))
        elif a.zyg == "het":
            if (modes & {"AD", "XD"} or not modes) and dom(a):
                out.append(Finding(symbol, (a,), "dominant", "full" if modes else "unknown_moi", "n/a", a.level))
            if modes & {"AR", "XR"} and rec(a):
                out.append(Finding(symbol, (a,), "recessive", "single_hit", "n/a", a.level))
    if pairs_tried(info):
        hets = [a for a in alleles if a.zyg == "het"]
        graph = phase_graph(hets, readback, settings.get("phase_sources", PHASE_SOURCES))
        for rank, h in enumerate(all_pairs([a for a in hets if rec(a)], graph)):
            a, b = sorted((h.a, h.b), key=lambda x: (x.first.call.pos, x.first.call.alt))
            fit = "panel_partner" if a.panel_only or b.panel_only else "full" if modes else "unknown_moi"
            out.append(Finding(symbol, (a, b), "recessive", fit, h.phase,
                               min(a.level, b.level, key=LEVELS.index), "+".join(h.sources), h.companions, rank))
    return out
