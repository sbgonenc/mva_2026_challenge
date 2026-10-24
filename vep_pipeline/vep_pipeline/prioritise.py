"""Rank findings into the top-N submission CSV.

Primary epcr = the tier's fixed value (prioritise.tiers), times p_trans_unphased for a pair not proven trans. Each gene
keeps its best finding and, at a phenotype >= alternative_pairs.min_phenotype, up to max_per_gene - 1 further pairs;
a gene's kept pairs exclude each other, so they share one budget (share_budget). The reported primary rows are then
scaled together so their sum is at most primary_sum_cap. Secondary findings (ACMG SF v3.3 genes without a phenotype
match) get the confidence that the variant is a reportable P/LP. All of these are heuristic numbers, not calibrated
probabilities (see README).
"""
from __future__ import annotations
import copy
import csv
import itertools
import shlex
from dataclasses import dataclass, replace
from pathlib import Path

import yaml

from .candidates import build_alleles, collect, findings_for_gene, pairs_tried, phase_graph
from .config import ConfigError
from .evidence import LEVELS, NULL_SO
from . import regnet
from .exomiser import gene_pheno, load_pheno, load_remm
from .genedb import g2p_path, load_g2p, load_sf
from .jobs import PKG_ROOT, _sbatch_header, final_vcf_path, work_dir
from .phasing import load_readback, readback_path
from .rank_settings import PHENOS, exomiser_settings, phasing_settings, regnet_settings

CSV_HEADER = ("proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2", "alt_2", "epcr",
              "finding_type", "notes")
GRADE_ORDER = ("null_lof", "clinvar_1star", "clinvar_2star")
TRUNCATING = NULL_SO - {"transcript_ablation"}


@dataclass(frozen=True)
class Scored:
    finding: object           # candidates.Finding
    finding_type: str         # 'primary' | 'secondary'
    tier: str                 # T1..T5, or 'SF'
    epcr: float
    pheno: str
    note: str = ""
    panel: str = ""           # gene-panel name for a panel gene; '<name> (phenotype floor)' when the floor raised pheno
    tier_epcr: float | None = None   # the tier value before the phase discount (secondary: the epcr)
    phase_factor: float = 1.0        # p_trans_unphased for a pair not proven trans, else 1


def tier_for(level, fit, pheno, tiers) -> tuple:
    for t in tiers:
        for rule in t["rules"]:
            if all(r in ("*", v) for r, v in zip(rule, (level, fit, pheno))):
                return t["name"], t["epcr"]
    raise AssertionError("unreachable: the last tier has a catch-all rule")


def phase_factor_for(f, p_trans) -> float:
    """p_trans for a pair not proven trans, else 1."""
    return p_trans if len(f.alleles) == 2 and f.phase != "trans" else 1.0


def primary_epcr(f, tiers, pheno, p_trans) -> tuple:
    name, e = tier_for(f.level, f.fit, pheno, tiers)
    return name, e * phase_factor_for(f, p_trans)


def _sf_grade(a, exclude) -> str | None:
    if "known_variant" not in exclude:
        stars = max((p.ev.stars for p in a.pieces if p.ev.clinvar == "plp"), default=None)
        if stars is not None and stars >= 2:
            return "clinvar_2star"
        if stars == 1:
            return "clinvar_1star"
    return "null_lof" if "null:strong" in a.reasons else None


def secondary_for(f, sf, exclude) -> str | None:
    """Grade of `f` as an ACMG SF secondary finding under the gene's reporting rule, or None if not reportable."""
    grades = [_sf_grade(a, exclude) for a in f.alleles]
    if None in grades:
        return None
    a = f.alleles[0]
    two = len(f.alleles) == 2 or a.zyg in ("hom", "hemi")
    ok = {"all_plp": True, "xl_any": True,
          "truncating_only": len(f.alleles) == 1 and bool(set(a.first.ev.consequence.split("&")) & TRUNCATING),
          "biallelic": two, "xl_biallelic": two,
          "hfe_c282y_hom": a.zyg == "hom" and "p.Cys282Tyr" in a.first.ev.hgvs}[sf.rule]
    return min(grades, key=GRADE_ORDER.index) if ok else None


def score_all(findings, pe, settings, sf_genes) -> tuple:
    tiers, sec = settings["tiers"], settings["secondary"]
    exclude = frozenset(settings["exclude_evidence_classes"])
    primary, secondary = [], []
    for f in findings:
        pheno, panel = gene_pheno(pe, f.gene, settings)
        sfg = sf_genes.get(f.gene) if sec["enabled"] else None
        if sfg and pheno not in ("strong", "moderate"):
            grade = secondary_for(f, sfg, exclude)
            if grade:
                note = (f"Incidental finding unrelated to primary phenotype; ACMG SF v3.3 {f.gene} "
                        f"({sfg.inheritance}, {sfg.rule}); {grade}"
                        + ("; phenotype relevance not assessed (Exomiser off)" if pheno == "none" else ""))
                e = sec["epcr"][grade]
                secondary.append(Scored(f, "secondary", "SF", e, pheno, note, panel, tier_epcr=e))
                continue
        name, e = tier_for(f.level, f.fit, pheno, tiers)
        k = phase_factor_for(f, settings["p_trans_unphased"])
        primary.append(Scored(f, "primary", name, e * k, pheno, panel=panel, tier_epcr=e, phase_factor=k))
    return primary, secondary


def _sort_key(s, pe):
    """epcr, then calibrated before panel_partner, Exomiser score, weakest level, the gene's all_pairs order (stronger
    allele, proven trans, panel class, rarity; Model 2's best_pair), popmax and position. A higher-epcr pair still wins
    over an earlier one in all_pairs (a proven-trans pair over a stronger unphased one; user ruling 2026-10-20)."""
    f = s.finding
    c = f.alleles[0].first.call
    combined = pe.combined.get(f.gene, 0.0) if pe else 0.0
    return (-s.epcr, f.fit == "panel_partner", -combined, -LEVELS.index(f.level), f.pair_rank,
            max(a.popmax for a in f.alleles), c.chrom, c.pos, c.ref, c.alt, f.gene)


def _keys(s) -> frozenset:
    return frozenset(a.key for a in s.finding.alleles)


def per_gene(rows, pe, settings) -> list:
    """Best finding per gene, sorted; a finding whose allele set was already kept (overlapping genes) is skipped. A
    gene whose phenotype is >= alternative_pairs.min_phenotype also keeps up to max_per_gene - 1 further pairs."""
    ap = settings["alternative_pairs"]
    k, weakest = ap["max_per_gene"], PHENOS.index(ap["min_phenotype"])
    count, variants, out = {}, set(), []
    for s in sorted(rows, key=lambda s: _sort_key(s, pe)):
        n = count.get(s.finding.gene, 0)
        further = n < k and len(s.finding.alleles) == 2 and PHENOS.index(s.pheno) <= weakest
        if _keys(s) in variants or (n and not further):
            continue
        count[s.finding.gene] = n + 1
        variants.add(_keys(s))
        out.append(s)
    return out


def best_per_gene(rows, pe) -> list:
    """per_gene with one finding per gene (Model 2)."""
    return per_gene(rows, pe, {"alternative_pairs": {"max_per_gene": 1, "min_phenotype": PHENOS[0]}})


def share_budget(rows, pe) -> list:
    """A gene's kept pairs exclude each other: with >= 2 of them, they share B = their largest tier value in proportion
    to their tier values, each then times its phase factor. Single rows (nested hypotheses) keep their epcr. Sorted."""
    pairs = {}
    for s in rows:
        if len(s.finding.alleles) == 2:
            pairs.setdefault(s.finding.gene, []).append(s)
    out = []
    for s in rows:
        group = pairs.get(s.finding.gene, ()) if len(s.finding.alleles) == 2 else ()
        if len(group) >= 2:
            total, b = sum(r.tier_epcr for r in group), max(r.tier_epcr for r in group)
            s = replace(s, epcr=b * s.tier_epcr / total * s.phase_factor if total else 0.0)
        out.append(s)
    return sorted(out, key=lambda s: _sort_key(s, pe))


def cap_primary(rows, cap) -> list:
    total = sum(s.epcr for s in rows)
    return [replace(s, epcr=s.epcr * cap / total) for s in rows] if total > cap else list(rows)


def select(primary, secondary, settings) -> list:
    n, k_max = settings["top_n"], settings["secondary"]["max_rows"]
    prim = primary[: n - min(k_max, len(secondary))]
    sec = secondary[: min(k_max, n - len(prim))]
    return cap_primary(prim, settings["primary_sum_cap"]) + sec


def _chrom(c, prefix) -> str:
    c = c[3:] if c.startswith("chr") else c
    return prefix + ("M" if c == "MT" and prefix == "chr" else c)


def changed(a) -> bool:
    """True when normalisation altered the record (atomised, split, trimmed or left-aligned)."""
    o, c = a.orig, a.first.call
    if not o:
        return False
    alt = o[3].split(",")[o[4] - 1]
    return len(a.pieces) > 1 or "," in o[3] or (o[0], o[1], o[2], alt) != (c.chrom, c.pos, c.ref, c.alt)


def coords(a, prefix) -> tuple:
    """(chrom, pos, ref, alt) of the original VCF record when normalisation changed it, else the normalised key."""
    o = a.orig
    if o:
        chrom, pos, ref, alts, idx = o
        return _chrom(chrom, prefix), pos, ref, alts.split(",")[idx - 1]
    c = a.first.call
    return _chrom(c.chrom, prefix), c.pos, c.ref, c.alt


def notes_for(s, pe) -> str:
    f = s.finding
    parts = ([s.note] if s.note else []) + [s.tier, f.gene, f"{f.hypothesis}/{f.fit}", f"phenotype={s.pheno}"]
    if s.panel:
        parts.append(f"panel={s.panel}")
    if pe and f.gene in pe.rank:
        parts.append(f"pheno_rank={pe.rank[f.gene]}")
    if pe and f.gene in pe.exomiser_rank:
        parts.append(f"exomiser_rank={pe.exomiser_rank[f.gene]}")
    for i, a in enumerate(f.alleles, 1):
        ev = a.first.ev
        parts.append(f"v{i}: {ev.hgvs or ev.consequence} {a.zyg} level={a.level} [{', '.join(a.reasons)}] "
                     f"popmax={a.popmax:.2g}" + (" low-VAF" if a.low_vaf else "")
                     + (f" norm={'+'.join(p.call.key for p in a.pieces)}" if changed(a) else ""))
    if len(f.alleles) == 2:
        src = f" ({f.phase_source})" if f.phase_source else ""
        parts.append(f"phase={f.phase}{src}; distance={f.distance}bp")
    return "; ".join(parts)


def csv_row(s, proband, prefix, pe) -> list:
    a = s.finding.alleles
    c1 = coords(a[0], prefix)
    c2 = coords(a[1], prefix) if len(a) == 2 else ("", "", "", "")
    return [proband, *map(str, c1), *map(str, c2), f"{s.epcr:.3f}", s.finding_type, notes_for(s, pe)]


@dataclass(frozen=True)
class Ranking:
    """The outcome of one settings scenario."""
    alleles: dict             # {symbol: [Allele]} after QC
    every: list               # every scored finding, primary then secondary
    kept: list                # per gene: the best finding and any alternative pairs, budget-shared; primary, secondary
    chosen: list              # the CSV rows (capped)


def rank(cands, genes, sf_genes, pe, settings, sex, readback=None) -> Ranking:
    """Score every finding, keep the rows per gene, share pair budgets and select the CSV rows for one scenario."""
    by_gene = build_alleles(cands, genes, settings, sex)
    findings = [f for g in sorted(by_gene) for f in findings_for_gene(g, by_gene[g], genes.get(g), settings, readback)]
    primary, secondary = score_all(findings, pe, settings, sf_genes)
    prim, sec = (share_budget(per_gene(rows, pe, settings), pe) for rows in (primary, secondary))
    return Ranking(by_gene, primary + secondary, prim + sec, select(prim, sec, settings))

# ----------------------------------------------------------------------------- outputs

CAND_COLUMNS = ("gene", "finding_type", "tier", "epcr_before_cap", "phenotype", "hypothesis", "fit", "level", "phase",
                "variant_1", "variant_2", "notes")

# Declared before the first real run; every scenario is reported, none is chosen afterwards.
# p_trans_unphased defaults to 0.5. Guo 2024 Nat Genet (doi:10.1038/s41588-023-01608-3): 86% of co-occurring rare
# pLoF pairs in gnomAD (general population) are in cis, whereas ~96% of phaseable pairs in their diagnosed patients
# were in trans (a selected cohort). No paper derives P(trans) for an affected singleton, so 0.5 is a neutral
# default and 0.15 / 0.9 bracket it.
# The last three isolate the gene panel, the alternative pair rows and read-backed phase. REGNET_SCENARIO (the
# regulatory-network second hit removed) is reported in addition only when regnet.enabled and regnet.second_hit, so a
# run without the hook keeps the same sensitivity table.
SCENARIOS = (
    ("base", {}),
    ("no_known_variant", {"exclude_evidence_classes": ["known_variant"]}),
    ("no_known_no_clinical_predictor", {"exclude_evidence_classes": ["known_variant", "clinically_trained_predictor"]}),
    ("pheno_rank_10_100", {"pheno_rank": {"strong": 10, "moderate": 100}}),
    ("pheno_rank_50_500", {"pheno_rank": {"strong": 50, "moderate": 500}}),
    ("p_trans_0.15", {"p_trans_unphased": 0.15}),
    ("p_trans_0.9", {"p_trans_unphased": 0.9}),
    ("no_gene_panel", {"gene_panel": {"name": "", "genes": [], "phenotype_floor": None, "deep_second_hit": False}}),
    ("one_finding_per_gene", {"alternative_pairs": {"max_per_gene": 1, "min_phenotype": "strong"}}),
    ("no_read_backed_phase", {"phase_sources": ["multiallelic", "ps", "pid"]}),
)
REGNET_SCENARIO = ("no_regulatory_network", {"exclude_evidence_classes": ["regulatory_network"]})


def write_csv(rows, path, proband, prefix, pe):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(CSV_HEADER)
        w.writerows(csv_row(s, proband, prefix, pe) for s in rows)


def write_excluded(excluded, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["variant", "region"])
        w.writerows(sorted(set(excluded)))


def write_candidates(scored, path, pe):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(CAND_COLUMNS)
        for s in sorted(scored, key=lambda s: _sort_key(s, pe)):
            f = s.finding
            keys = [a.key for a in f.alleles] + [""]
            w.writerow([f.gene, s.finding_type, s.tier, f"{s.epcr:.3f}", s.pheno, f.hypothesis, f.fit, f.level,
                        f.phase, keys[0], keys[1], notes_for(s, pe)])


CH_COLUMNS = ("gene", "variant_1", "variant_2", "hgvs_1", "hgvs_2", "level_1", "level_2", "reasons_1", "reasons_2",
              "phase", "phase_sources", "distance", "collapsed_into", "fit", "tier", "epcr_before_cap", "reportable",
              "csv_row")


def _scored_pairs(r) -> dict:
    """{gene: {allele-key set: scored pair finding}}: the kept (budget-shared) row when kept, else the best one."""
    out = {}
    for s in r.every:
        pairs = out.setdefault(s.finding.gene, {})
        if len(s.finding.alleles) == 2 and (_keys(s) not in pairs or s.epcr > pairs[_keys(s)].epcr):
            pairs[_keys(s)] = s
    for s in r.kept:
        if len(s.finding.alleles) == 2:
            out[s.finding.gene][_keys(s)] = s
    return out


def _pair_status(x, y, ph, scored) -> tuple:
    """(reportable, collapsed_into) of two het alleles of one gene, given the gene's scored pairs."""
    keys = frozenset((x.key, y.key))
    if len(keys) == 1:
        return "no (same variant)", ""
    if ph == "cis":
        return "no (cis)", ""
    if x.second_hit_only and y.second_hit_only:
        return "no (two second-hit-only alleles)", ""
    if keys in scored:
        return "yes", ""
    rep = next((s.finding for s in scored.values()
                if keys <= {m.key for m in s.finding.alleles + s.finding.companions}), None)
    return ("no (collapsed)", "|".join(a.key for a in rep.alleles)) if rep else ("no (not paired)", "")


def _gene_pair_rows(g, alleles, settings, readback, scored, in_csv) -> list:
    """compound_het_rows for one gene; `scored` holds the gene's scored pairs (see _scored_pairs)."""
    hets = [a for a in alleles if a.zyg == "het"]
    graph = phase_graph(hets, readback, settings["phase_sources"])
    rows = []
    for x, y in itertools.combinations(sorted(hets, key=lambda a: (a.first.call.pos, a.first.call.alt)), 2):
        ph, srcs = graph.relation(x, y)
        status, into = _pair_status(x, y, ph, scored)
        s = scored[frozenset((x.key, y.key))] if status == "yes" else None
        rows.append([g, x.key, y.key, x.first.ev.hgvs, y.first.ev.hgvs, x.level, y.level, ", ".join(x.reasons),
                     ", ".join(y.reasons), ph, "+".join(srcs), abs(y.first.call.pos - x.first.call.pos), into,
                     *((s.finding.fit, s.tier, f"{s.epcr:.3f}", status, in_csv.get(_keys(s), "")) if s
                       else ("", "", "", status, ""))])
    return rows


def compound_het_rows(r, genes, settings, readback=None) -> list:
    """One row (CH_COLUMNS) per pair of het alleles in each recessive or unknown-MOI gene, reportable or not; gene,
    then position order. Phase comes from the same phase graph as the ranking; csv_row is the 1-based CSV row."""
    scored, in_csv = _scored_pairs(r), {}
    for i, s in enumerate(r.chosen, 1):
        in_csv.setdefault(_keys(s), i)
    return [row for g in sorted(r.alleles) if pairs_tried(genes.get(g))
            for row in _gene_pair_rows(g, r.alleles[g], settings, readback, scored.get(g, {}), in_csv)]


def write_compound_het(r, genes, settings, path, readback=None):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(CH_COLUMNS)
        w.writerows(compound_het_rows(r, genes, settings, readback))


def scenario_settings(settings, change) -> dict:
    s = copy.deepcopy(settings)
    for k, v in change.items():
        s[k] = sorted(set(s[k]) | set(v)) if k == "exclude_evidence_classes" else v
    return s


def write_sensitivity(cands, genes, sf_genes, pe, settings, sex, path, readback=None, scenarios=SCENARIOS):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["scenario", "row", "gene", "finding_type", "tier", "epcr", "variant_1", "variant_2"])
        for name, change in scenarios:
            chosen = rank(cands, genes, sf_genes, pe, scenario_settings(settings, change), sex, readback).chosen
            for i, s in enumerate(chosen, 1):
                keys = [a.key for a in s.finding.alleles] + [""]
                w.writerow([name, i, s.finding.gene, s.finding_type, s.tier, f"{s.epcr:.3f}", keys[0], keys[1]])


def run(cfg, settings, sensitivity=False) -> tuple:
    """Scan the VEP VCF once, rank, and write the CSV, the candidates, excluded and compound-het TSVs, (optionally)
    the sensitivity table and, with regnet.enabled, the regulatory TSV and its summary; with regnet.second_hit too, its
    second_hit rows join the candidates (regnet.add_second_hits) before ranking."""
    ex = exomiser_settings(cfg)
    ph = phasing_settings(cfg)
    vcf = final_vcf_path(cfg)
    if not vcf.is_file():
        raise ConfigError(f"{vcf} not found: run `vep-pipeline run` first")
    genes = load_g2p(g2p_path(cfg, settings), settings["g2p"]["confidence"])
    sf_genes = load_sf(settings["secondary"]["genes"])
    pe = load_pheno(cfg, ex) if ex["enabled"] else None
    remm = load_remm(cfg, ex) if ex["enabled"] else {}
    readback = None
    if ph["enabled"]:
        side = readback_path(cfg, ph)
        if not side.is_file():
            raise ConfigError(f"{side} not found: run `vep-pipeline phase --submit` first")
        readback = load_readback(side)
    rn = regnet_settings(cfg)
    res = regnet.load_resource(Path(rn["build_dir"])) if rn["enabled"] else None
    gate = None if "population" in settings["exclude_evidence_classes"] else settings["af"]["max_recessive"]
    hook = regnet.HitCollector(res.index, tuple(settings["af"]["populations"]), gate) if res else None
    out = Path(settings["output"])
    out.parent.mkdir(parents=True, exist_ok=True)
    excluded = []
    cands = collect(vcf, cfg["sample"], settings, genes, remm, excluded, hook=hook)
    if res:
        rows, summary = regnet.score(hook.hits, res, pe, remm, settings, rn)
        regnet.write_regulatory(rows, out.with_suffix(".regulatory.tsv"), settings["chrom_prefix"])
        regnet.write_summary(summary, out.with_suffix(".regulatory.summary.tsv"))
        if rn["second_hit"]:
            cands = regnet.add_second_hits(cands, rows, remm, tuple(settings["af"]["populations"]))
    r = rank(cands, genes, sf_genes, pe, settings, ex["sex"], readback)
    write_csv(r.chosen, out, cfg["sample"], settings["chrom_prefix"], pe)
    cand_path = out.with_suffix(".candidates.tsv")
    write_candidates(r.kept, cand_path, pe)
    write_excluded(excluded, out.with_suffix(".excluded.tsv"))
    write_compound_het(r, genes, settings, out.with_suffix(".compound_het.tsv"), readback)
    if sensitivity:
        scenarios = SCENARIOS + (REGNET_SCENARIO,) if res and rn["second_hit"] else SCENARIOS
        write_sensitivity(cands, genes, sf_genes, pe, settings, ex["sex"], out.with_suffix(".sensitivity.tsv"),
                          readback, scenarios)
    return out, cand_path


def render_slurm(cfg, cli_args) -> Path:
    """work/prioritise.sbatch running this command on a compute node with a frozen copy of the resolved config."""
    wd = work_dir(cfg)
    wd.mkdir(parents=True, exist_ok=True)
    Path(cfg["slurm"]["log_dir"]).mkdir(parents=True, exist_ok=True)
    cfg_path = wd / "prioritise.config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    q = shlex.quote
    script = wd / "prioritise.sbatch"
    script.write_text(_sbatch_header(cfg, "prioritise", "prioritise_%j.out", 1)
                      + f"export PYTHONPATH={q(str(PKG_ROOT))}${{PYTHONPATH:+:$PYTHONPATH}}\n"
                      + f"micromamba run -n {q(cfg['tools']['env'])} python -m vep_pipeline.cli prioritise "
                      + f"-c {q(str(cfg_path))} {' '.join(q(a) for a in cli_args)}\n")
    return script
