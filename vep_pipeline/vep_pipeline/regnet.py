"""Rank rare non-coding variants that fall in scDoRI gene-linked peaks, using S-LDSC TF enrichment.

This module holds the peak index, the hit capture hook for `candidates.collect` and the pure scoring helpers; see
the README section "Regulatory network (`regnet`)". Stdlib only.
"""
from __future__ import annotations
import csv
import random
import statistics
from bisect import bisect_right
from dataclasses import dataclass, replace
from math import comb
from pathlib import Path

from .candidates import Candidate, _bare, qc_ok
from .config import ConfigError
from .evidence import (CADD_NONCODING, INFRAME, NULL_SO, PROTEIN_ALTERING, REMM_SECOND_HIT, SPLICE_REGION, Evidence,
                       _max, clinvar_class, gene_evidence)
from .exomiser import gene_pheno, pheno_strength
from .regnet_build import bh

__all__ = ["CODING", "HitCollector", "PeakIndex", "RegHit", "RegRow", "Resource", "SldscRow", "add_second_hits", "bh",
           "hypergeom_sf", "load_resource", "reg_tier", "regulon_flags", "score", "sldsc_flags", "trait_null", "write_regulatory",
           "write_summary"]

CODING = (PROTEIN_ALTERING | NULL_SO | INFRAME | SPLICE_REGION |
          {"synonymous_variant", "stop_retained_variant", "start_retained_variant",
           "incomplete_terminal_codon_variant", "coding_sequence_variant"})


@dataclass(frozen=True)
class RegHit:
    call: object        # candidates.Call
    entries: tuple      # CSQ dicts of this record (kept for Evidence later)
    peaks: tuple        # overlapping linked peak names
    popmax: float
    cadd: float | None


class PeakIndex:
    """Overlap lookup over {peak: (chrom, start, end)}; start and end are inclusive, contigs match with or without 'chr'."""

    def __init__(self, peaks: dict):
        by_contig: dict = {}
        for name, (chrom, start, end) in peaks.items():
            by_contig.setdefault(_bare(chrom), []).append((start, end, name))
        self._starts, self._peaks, self._maxw = {}, {}, {}
        for contig, rows in by_contig.items():
            rows.sort()
            self._peaks[contig] = rows
            self._starts[contig] = [r[0] for r in rows]
            self._maxw[contig] = max(r[1] - r[0] for r in rows)

    def find(self, chrom: str, start: int, end: int) -> list[str]:
        """Names of peaks overlapping [start, end], sorted by peak start."""
        contig = _bare(chrom)
        rows = self._peaks.get(contig)
        if not rows:
            return []
        out = []
        i = bisect_right(self._starts[contig], end) - 1
        while i >= 0 and rows[i][0] >= start - self._maxw[contig]:
            if rows[i][1] >= start:
                out.append(rows[i][2])
            i -= 1
        return out[::-1]


def _terms(entries):
    return {t for e in entries for t in e.get("Consequence", "").split("&") if t}


class HitCollector:
    """Hook for `candidates.collect(hook=...)`: keeps rare non-coding records that overlap a linked peak."""

    def __init__(self, index: PeakIndex, pops: tuple, gate: float | None):
        self.index, self.gate, self.hits = index, gate, []
        self.af_fields = tuple(f"gnomAD{s}_{p}_AF" for s in ("e", "g") for p in pops)   # as evidence.gene_evidence

    def __call__(self, call, entries) -> None:
        if not entries or call.alt == "*":            # no CSQ entry, or a spanning-deletion allele
            return
        end = call.pos + max(len(call.ref), 1) - 1
        peaks = self.index.find(call.chrom, call.pos, end)
        if not peaks or _terms(entries) & CODING:
            return
        popmax = _max(entries, self.af_fields) or 0.0
        if self.gate is not None and popmax > self.gate:
            return
        self.hits.append(RegHit(call, tuple(entries), tuple(peaks), popmax, _max(entries, ("CADD_PHRED",))))


def hypergeom_sf(k: int, N: int, K: int, n: int) -> float:
    """P(X >= k) for X ~ Hypergeometric(population N, K successes, n draws)."""
    lo, hi = max(k, 0), min(K, n)
    if lo > hi:
        return 0.0
    return sum(comb(K, i) * comb(N - K, n - i) for i in range(lo, hi + 1)) / comb(N, n)


def reg_tier(pheno: str, v: bool, s: bool) -> str:
    """R1..R5 from the gene's phenotype strength and the variant (v) and S-LDSC TF (s) support flags."""
    if pheno == "strong":
        return "R1" if v and s else "R2" if v or s else "R3"
    if pheno == "moderate":
        return "R2" if v and s else "R3" if v or s else "R4"
    return "R5"


# ----------------------------------------------------------------------------- compiled resource

BUILD_FILES = ("linked_peaks.tsv", "peak_tfs.tsv", "regulons.tsv", "gene_universe.txt", "sldsc.tsv", "traits.tsv")


@dataclass(frozen=True)
class SldscRow:
    trait: str
    tf: str
    coef: float         # S-LDSC tau
    fdr: float          # BH within the trait, across TFs


@dataclass(frozen=True)
class Resource:
    index: PeakIndex
    links: dict         # peak -> tuple of linked genes
    peak_tfs: dict      # peak -> tuple of TFs with an eGRN edge to it (any cell type)
    regulons: dict      # tf -> frozenset of target genes (union over cell types)
    universe: frozenset # scDoRI-modelled genes
    sldsc: list         # SldscRow
    traits: dict        # trait -> label ('' when unlabelled)


def _tsv(path: Path) -> list:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_resource(build_dir: Path) -> Resource:
    """Read the TSVs `regnet-build` wrote into build_dir."""
    build_dir = Path(build_dir)
    for f in BUILD_FILES:
        if not (build_dir / f).is_file():
            raise ConfigError(f"{build_dir / f}: regnet build file missing; run `bin/vep-pipeline regnet-build`")
    links, coords = {}, {}
    for r in _tsv(build_dir / "linked_peaks.tsv"):
        links.setdefault(r["peak"], []).append(r["gene"])
        coords[r["peak"]] = (r["chrom"], int(r["start"]), int(r["end"]))
    peak_tfs, regulons = {}, {}
    for r in _tsv(build_dir / "peak_tfs.tsv"):
        peak_tfs.setdefault(r["peak"], []).append(r["tf"])
    for r in _tsv(build_dir / "regulons.tsv"):
        regulons.setdefault(r["tf"], set()).add(r["gene"])
    universe = frozenset(g for g in (build_dir / "gene_universe.txt").read_text().split() if g)
    sldsc = [SldscRow(r["trait"], r["tf"], float(r["coef"]), float(r["fdr"])) for r in _tsv(build_dir / "sldsc.tsv")]
    traits = {r["trait"]: r["label"] for r in _tsv(build_dir / "traits.tsv")}
    return Resource(PeakIndex(coords), {p: tuple(g) for p, g in links.items()},
                    {p: tuple(t) for p, t in peak_tfs.items()}, {t: frozenset(g) for t, g in regulons.items()},
                    universe, sldsc, traits)


# ----------------------------------------------------------------------------- TF flags and the trait null

def sldsc_flags(sldsc_rows, traits, fdr: float) -> dict:
    """{tf: [(trait, fdr), ...] sorted by trait} for rows of a declared trait with tau > 0 and FDR < the cut."""
    want = set(traits)
    out = {}
    for r in sldsc_rows:
        if r.trait in want and r.coef > 0 and r.fdr < fdr:
            out.setdefault(r.tf, []).append((r.trait, r.fdr))
    return {tf: sorted(v) for tf, v in out.items()}


def regulon_flags(regulons, universe, pheno_genes, fdr: float) -> dict:
    """{tf: BH q} for TFs whose regulon (within the universe) is enriched for phenotype genes (hypergeometric, BH
    across the TFs with >= 1 target in the universe), q < the cut. Empty without phenotype genes."""
    hits = set(pheno_genes) & universe
    if not hits:
        return {}
    tested, pvals = [], []
    for tf in sorted(regulons):
        targets = regulons[tf] & universe
        if targets:
            tested.append(tf)
            pvals.append(hypergeom_sf(len(targets & hits), len(universe), len(hits), len(targets)))
    return {tf: q for tf, q in zip(tested, bh(pvals)) if q < fdr}


def trait_null(sldsc_rows, all_traits, k: int, fdr: float, draws: int, seed: int) -> list[int]:
    """Number of distinct S-LDSC-flagged TFs for `draws` random sets of k traits from all_traits (seeded)."""
    by_trait = {}
    for tf, flags in sldsc_flags(sldsc_rows, all_traits, fdr).items():
        for trait, _q in flags:
            by_trait.setdefault(trait, set()).add(tf)
    rng, pool = random.Random(seed), sorted(all_traits)
    out = []
    for _ in range(draws):
        out.append(len(set().union(*(by_trait.get(t, ()) for t in rng.sample(pool, k)))))
    return out


# ----------------------------------------------------------------------------- scoring

TIERS = ("R1", "R2", "R3", "R4", "R5")


@dataclass(frozen=True)
class RegRow:
    hit: RegHit
    gene: str
    peak: str
    pheno: str
    panel: str
    v: bool             # variant signal: CADD >= CADD_NONCODING or ReMM >= REMM_SECOND_HIT (classes not excluded)
    tier: str
    tfs: tuple          # every TF on the peak
    sldsc_tfs: tuple    # ((tf, ((trait, fdr), ...)), ...) for the S-LDSC-flagged TFs on the peak
    regulon_tfs: tuple  # ((tf, q), ...) for the regulon-flagged TFs on the peak
    remm: float | None
    second_hit: bool


def _variant_signal(cadd, remm, exclude) -> bool:
    return ((cadd is not None and cadd >= CADD_NONCODING and "unsupervised_predictor" not in exclude)
            or (remm is not None and remm >= REMM_SECOND_HIT and "clinically_trained_predictor" not in exclude))


def _dedupe_key(r: RegRow) -> tuple:
    """Best row of one (variant, gene) over the variant's peaks."""
    return (r.tier, -len(r.sldsc_tfs), -len(r.regulon_tfs), r.peak)


def _row_key(r: RegRow) -> tuple:
    return (r.tier, -len(r.sldsc_tfs), -len(r.regulon_tfs), -(r.remm or 0), -(r.hit.cadd or 0), r.hit.call.key, r.gene)


def score(hits, res: Resource, pe, remm, settings, rn) -> tuple[list[RegRow], dict]:
    """One row per (QC-passing hit, linked gene), best tier over the variant's peaks; rows sorted; plus the summary."""
    remm = remm or {}
    exclude = frozenset(settings["exclude_evidence_classes"])
    pheno_of = {g: gene_pheno(pe, g, settings) for g in res.universe}
    pheno_genes = {g for g, (p, _panel) in pheno_of.items() if p in ("strong", "moderate")}
    reg_flags = regulon_flags(res.regulons, res.universe, pheno_genes, rn["fdr"])
    traits = list(dict.fromkeys(rn["traits"]))
    sl_flags = sldsc_flags(res.sldsc, traits, rn["fdr"])
    best, n_qc = {}, 0
    for hit in hits:
        if not qc_ok(hit.call, settings["qc"]):
            continue
        n_qc += 1
        r = remm.get(hit.call.key)
        v = _variant_signal(hit.cadd, r, exclude)
        for peak in hit.peaks:
            tfs = res.peak_tfs.get(peak, ())
            s_tfs = tuple((tf, tuple(sl_flags[tf])) for tf in tfs if tf in sl_flags)
            r_tfs = tuple((tf, reg_flags[tf]) for tf in tfs if tf in reg_flags)
            for gene in res.links.get(peak, ()):
                pheno, panel = pheno_of.get(gene) or gene_pheno(pe, gene, settings)
                row = RegRow(hit, gene, peak, pheno, panel, v, reg_tier(pheno, v, bool(s_tfs)), tfs, s_tfs, r_tfs,
                             r, pheno == "strong" and v)
                key = (hit.call.key, gene)
                old = best.get(key)
                if old is None or _dedupe_key(row) < _dedupe_key(old):
                    best[key] = row
    rows = sorted(best.values(), key=_row_key)
    null = trait_null(res.sldsc, list(res.traits), len(traits), rn["fdr"], rn["null_draws"], rn["seed"])
    observed = len(sl_flags)
    summary = {"hits": len(hits), "hits_qc_pass": n_qc, "rows": len(rows)}
    summary.update({f"rows_{t}": sum(r.tier == t for r in rows) for t in TIERS})
    summary["second_hit_rows"] = sum(r.second_hit for r in rows)
    linked = {r.gene for r in rows}
    summary["linked_genes_with_hits"] = len(linked)
    # Exomiser's own ranking (before the panel floor): an unranked gene is 'unknown' and so R5; a symbol mismatch
    # between scDoRI and Exomiser shows up here.
    summary["linked_genes_unranked"] = sum(pheno_strength(pe, g, settings["pheno_rank"]) == "unknown" for g in linked)
    summary["declared_traits"] = [f"{t}:{res.traits.get(t, '')}" for t in traits]
    for t in traits:
        summary[f"flagged_tfs_{t}"] = sorted(tf for tf, fl in sl_flags.items() if any(x == t for x, _q in fl))
    summary.update({"observed_flagged_tfs": observed, "null_draws": len(null),
                    "null_median": statistics.median(null) if null else None,
                    "null_p": (1 + sum(x >= observed for x in null)) / (len(null) + 1),
                    "pheno_genes_in_universe": len(pheno_genes),
                    "regulon_flagged_tfs": [f"{tf}({q:.4g})" for tf, q in sorted(reg_flags.items())]})
    return rows, summary


# ----------------------------------------------------------------------------- second-hit hook into the ranking

def _grn(row: RegRow) -> str:
    tfs = ",".join(tf for tf, _flags in row.sldsc_tfs)
    return f"peak {row.peak} -> {row.gene}; {row.tier}" + (f"; S-LDSC TFs {tfs}" if tfs else "")


def _grn_classes(row: RegRow) -> frozenset:
    """Predictor classes behind the row's V, from its raw CADD and ReMM (exclusions ignored: variant_level re-gates)."""
    out = set()
    if row.hit.cadd is not None and row.hit.cadd >= CADD_NONCODING:
        out.add("unsupervised_predictor")
    if row.remm is not None and row.remm >= REMM_SECOND_HIT:
        out.add("clinically_trained_predictor")
    return frozenset(out)


def _link_evidence(row: RegRow, pops: tuple, remm: float | None) -> Evidence:
    """Evidence of a hit for its linked gene. With a CSQ entry for the gene: gene_evidence, plus the variant's ReMM.
    Link-only (the gene comes from the peak link alone): a record that contributes nothing but the grn reason (set by
    the caller) and the ClinVar B/LB veto, so the hook stays supporting and second-hit-only, and excluding
    regulatory_network removes it entirely. Popmax, consequence and the peak name (as hgvs) are kept for reporting."""
    hit, gene = row.hit, row.gene
    if any(e.get("SYMBOL") == gene for e in hit.entries):
        return replace(gene_evidence(list(hit.entries), pops)[gene], remm=remm)
    cv = next((e for e in hit.entries if e.get("ClinVar_CLNSIG")), {})
    clin, stars = clinvar_class(cv.get("ClinVar_CLNSIG", ""), cv.get("ClinVar_CLNREVSTAT", ""))
    clin, stars = (clin, stars) if clin == "blb" else (None, 0)
    return Evidence(symbol=gene, consequence="&".join(sorted(_terms(hit.entries))), hgvs=row.peak, popmax=hit.popmax,
                    clinvar=clin, stars=stars, revel=None, am=None, splice=None, cadd=None, nmd_escape=False,
                    loeuf=None, phaplo=None, noncoding_signal=False, remm=None)


def add_second_hits(cands, rows, remm, pops) -> list:
    """A new candidate list with each second_hit row's variant as a candidate of its linked gene, carrying Evidence.grn
    and Evidence.grn_classes (a supporting, second-hit-only signal in evidence.variant_level, dropped when every
    predictor class behind V is excluded). An existing (variant, gene) candidate is replaced in place, keeping its
    Evidence and only setting grn, grn_classes (and remm when it had none), so build_alleles never sees a duplicate
    piece; new candidates are appended in row order (see _link_evidence). The first row of a (variant, gene) wins."""
    remm = remm or {}
    out = list(cands)
    at = {}
    for i, c in enumerate(out):
        at.setdefault((c.call.key, c.ev.symbol), i)
    done = set()
    for row in rows:
        key = (row.hit.call.key, row.gene)
        if not row.second_hit or key in done:
            continue
        done.add(key)
        r = remm.get(row.hit.call.key)
        if key in at:
            c = out[at[key]]
            out[at[key]] = replace(c, ev=replace(c.ev, grn=_grn(row), grn_classes=_grn_classes(row),
                                                 remm=c.ev.remm if c.ev.remm is not None else r))
        else:
            out.append(Candidate(row.hit.call, replace(_link_evidence(row, pops, r), grn=_grn(row),
                                                       grn_classes=_grn_classes(row))))
    return out


# ----------------------------------------------------------------------------- outputs

REG_HEADER = ("chrom", "pos", "ref", "alt", "gene", "tier", "pheno", "panel", "second_hit", "peak", "popmax", "cadd",
              "remm", "v", "n_tfs", "sldsc_tfs", "regulon_tfs", "consequence")


def _chrom(c, prefix) -> str:
    # Same as prioritise._chrom (not imported: prioritise imports this module).
    c = c[3:] if c.startswith("chr") else c
    return prefix + ("M" if c == "MT" and prefix == "chr" else c)


def _num(x) -> str:
    return "" if x is None else f"{x:.4g}"


def _yn(b) -> str:
    return "yes" if b else "no"


def write_regulatory(rows, path, prefix) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(REG_HEADER)
        for r in rows:
            c = r.hit.call
            w.writerow([_chrom(c.chrom, prefix), c.pos, c.ref, c.alt, r.gene, r.tier, r.pheno, r.panel,
                        _yn(r.second_hit), r.peak, _num(r.hit.popmax), _num(r.hit.cadd), _num(r.remm), _yn(r.v),
                        len(r.tfs),
                        ";".join(f"{tf}({','.join(f'{t}:{q:.4g}' for t, q in fl)})" for tf, fl in r.sldsc_tfs),
                        ";".join(f"{tf}({q:.4g})" for tf, q in r.regulon_tfs),
                        "&".join(sorted(_terms(r.hit.entries)))])


def _value(v) -> str:
    if isinstance(v, (list, tuple)):
        return ";".join(map(str, v))
    if isinstance(v, float):
        return f"{v:.4g}"
    return "" if v is None else str(v)


def write_summary(summary, path) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["key", "value"])
        w.writerows([k, _value(v)] for k, v in summary.items())
