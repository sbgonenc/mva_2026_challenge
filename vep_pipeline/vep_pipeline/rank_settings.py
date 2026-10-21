"""Defaults and validation for the optional `exomiser:`, `phasing:` and `prioritise:` config sections.

Both sections may be absent from config.yaml; every key has a default here. CLI overrides (--top-n, -o) are applied
before the output-path default is computed, so the default file name `{prefix}.top{N}.csv` follows --top-n.
"""
from __future__ import annotations
import copy
import re
from pathlib import Path

from .candidates import PHASE_SOURCES
from .config import ConfigError
from .evidence import LEVELS, MISSENSE_TOOLS

PKG_DIR = Path(__file__).resolve().parent
FITS = ("full", "single_hit", "unknown_moi", "panel_partner")
PANEL_FLOORS = (None, "strong", "moderate")
PHENOS = ("strong", "moderate", "weak", "unknown", "none")
EXCLUDABLE = ("known_variant", "clinically_trained_predictor", "unsupervised_predictor", "population",
              "gene_prior_dominant")
SEXES = ("MALE", "FEMALE", "UNKNOWN_SEX", "OTHER_SEX")
GNOMAD_GROUPS = ("AFR", "AMR", "ASJ", "EAS", "FIN", "MID", "NFE", "SAS", "AMI", "REMAINING")  # gnomAD v4
_HPO = re.compile(r"^HP:\d{7}$")

# rule = [variant level, MOI fit, phenotype]; '*' matches anything; the first matching tier wins.
DEFAULT_TIERS = [
    {"name": "T1", "epcr": 0.90, "rules": [["strong", "full", "strong"]]},
    {"name": "T2", "epcr": 0.60, "rules": [["strong", "full", "moderate"], ["moderate", "full", "strong"],
                                           ["strong", "unknown_moi", "strong"]]},
    {"name": "T3", "epcr": 0.30, "rules": [["moderate", "full", "moderate"], ["supporting", "full", "strong"],
                                           ["moderate", "unknown_moi", "strong"]]},
    {"name": "T4", "epcr": 0.10, "rules": [["*", "full", "moderate"], ["*", "unknown_moi", "moderate"],
                                           ["strong", "single_hit", "strong"], ["strong", "full", "unknown"],
                                           ["strong", "full", "none"], ["*", "panel_partner", "strong"]]},
    {"name": "T5", "epcr": 0.02, "rules": [["*", "*", "*"]]},
]

EXOMISER_DEFAULTS = {
    "enabled": True,
    "hpo": [],
    "sex": "UNKNOWN_SEX",
    "phenopacket": None,
    "root": "/fast/AG_Bunina/Berk/tools/exomiser",
    "mem": "48G",
    "time": "24:00:00",
    "results_dir": None,                       # default {output.dir}/exomiser
}

PHASING_DEFAULTS = {
    "enabled": False,
    "bam": None,
    "reference": None,
    "env": None,                               # micromamba env with whatshap and samtools
    "results_dir": None,                       # default {output.dir}/phasing
    "pad": 1000,                               # bp around each candidate variant fetched from the BAM
    "indels": True,
    "mem": "16G",
    "time": "08:00:00",
}

DEFAULT_EXCLUDE_REGIONS = [{"name": "MHC", "chrom": "6", "start": 28510120, "end": 33480577}]  # GRC MHC, GRCh38

PRIORITISE_DEFAULTS = {
    "top_n": 10,
    "output": None,                            # default {output.dir}/{output.prefix}.top{top_n}.csv
    "chrom_prefix": "chr",
    "exclude_evidence_classes": [],
    "missense_tool": "alphamissense",
    "exclude_regions": DEFAULT_EXCLUDE_REGIONS,
    "qc": {"pass_only": True, "min_gq": 20, "min_dp": 10, "low_vaf": 0.30},
    "af": {"populations": ["AFR", "AMR", "EAS", "NFE", "SAS"], "max_dominant": 0.0001, "max_recessive": 0.01},
    "g2p": {"file": None, "confidence": ["definitive", "strong", "moderate"]},
    "pheno_rank": {"strong": 20, "moderate": 200},
    "phase_sources": list(PHASE_SOURCES),
    "gene_panel": {"name": "", "genes": [], "phenotype_floor": None, "deep_second_hit": False},   # off: Model 2
    "alternative_pairs": {"max_per_gene": 1, "min_phenotype": "strong"},                          # 1: Model 2
    "p_trans_unphased": 0.5,
    "primary_sum_cap": 0.95,
    "tiers": DEFAULT_TIERS,
    "secondary": {"enabled": True, "max_rows": 2, "genes": str(PKG_DIR / "acmg_sf_v3.3.tsv"),
                  "epcr": {"clinvar_2star": 0.30, "clinvar_1star": 0.20, "null_lof": 0.10}},
}


def _merge(default, user):
    if isinstance(default, dict):
        if user is None:
            return copy.deepcopy(default)
        if not isinstance(user, dict):
            raise ConfigError(f"expected a mapping, got {user!r}")
        unknown = set(user) - set(default)
        if unknown:
            raise ConfigError(f"unknown key(s) {sorted(unknown)}; expected some of {sorted(default)}")
        return {k: _merge(v, user[k]) if k in user else copy.deepcopy(v) for k, v in default.items()}
    return copy.deepcopy(default if user is None else user)


def _is_num(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _is_int(v) -> bool:
    return isinstance(v, int) and not isinstance(v, bool)


def _number(value, where, lo=0.0, hi=1.0, lo_open=False):
    """Real number (not bool) in [lo, hi], or (lo, hi] when lo_open."""
    if not _is_num(value):
        raise ConfigError(f"{where} must be a number, got {value!r}")
    in_range = value > lo if lo_open else value >= lo
    if not (in_range and value <= hi):
        interval = f"{'(' if lo_open else '['}{lo}, {hi}]"
        raise ConfigError(f"{where} must be a number in {interval}, got {value!r}")
    return value


def _bool(value, where):
    if not isinstance(value, bool):
        raise ConfigError(f"{where} must be true or false, got {value!r}")
    return value


def _section(cfg, name, defaults) -> dict:
    raw = cfg.get(name)
    if isinstance(raw, bool):
        if "enabled" in defaults:          # exomiser: false / true, as for Plugin in config._normalise
            raw = {"enabled": raw}
        elif raw:
            raw = None
        else:
            raise ConfigError(f"{name}: false is not valid here; remove the section to use the defaults")
    try:
        return _merge(defaults, raw)
    except ConfigError as e:
        raise ConfigError(f"{name}: {e}") from None


def exomiser_settings(cfg: dict) -> dict:
    ex = _section(cfg, "exomiser", EXOMISER_DEFAULTS)
    _bool(ex["enabled"], "exomiser.enabled")
    if ex["sex"] not in SEXES:
        raise ConfigError(f"exomiser.sex must be one of {SEXES}, got {ex['sex']!r}")
    if not isinstance(ex["hpo"], list):
        raise ConfigError(f"exomiser.hpo must be a list of HPO ids, got {ex['hpo']!r}")
    bad = [t for t in ex["hpo"] if not _HPO.match(str(t))]
    if bad:
        raise ConfigError(f"exomiser.hpo: not HPO ids (HP:0000000): {bad}")
    if ex["enabled"] and not ex["hpo"] and not ex["phenopacket"]:
        raise ConfigError("exomiser.enabled needs exomiser.hpo (your HPO terms) or exomiser.phenopacket")
    ex["results_dir"] = ex["results_dir"] or str(Path(cfg["output"]["dir"]) / "exomiser")
    return ex


def _alignment_index(bam: Path) -> Path:
    """The existing index of a BAM/CRAM (`<file>.bai|.crai` or `<stem>.bai|.crai`), else ConfigError."""
    ext = {".bam": ".bai", ".cram": ".crai"}.get(bam.suffix)
    if ext is None:
        raise ConfigError(f"phasing.bam must be a BAM or CRAM file, got {str(bam)!r}")
    found = [p for p in (Path(f"{bam}{ext}"), bam.with_suffix(ext)) if p.is_file()]
    if not found:
        raise ConfigError(f"phasing.bam: index {bam}{ext} (or {bam.with_suffix(ext).name}) not found; "
                          f"run `samtools index`")
    return found[0]


def phasing_settings(cfg: dict) -> dict:
    ph = _section(cfg, "phasing", PHASING_DEFAULTS)
    _bool(ph["enabled"], "phasing.enabled")
    _bool(ph["indels"], "phasing.indels")
    if not _is_int(ph["pad"]) or ph["pad"] < 1:
        raise ConfigError(f"phasing.pad must be a positive integer, got {ph['pad']!r}")
    for k in ("mem", "time"):
        if not isinstance(ph[k], str):
            raise ConfigError(f"phasing.{k} must be a string, got {ph[k]!r}")
    ph["results_dir"] = ph["results_dir"] or str(Path(cfg["output"]["dir"]) / "phasing")
    if ph["enabled"]:
        for k in ("bam", "reference", "env"):
            if not ph[k] or not isinstance(ph[k], str):
                raise ConfigError(f"phasing.{k} must be set when phasing.enabled, got {ph[k]!r}")
        for k in ("bam", "reference"):
            if not Path(ph[k]).is_file():
                raise ConfigError(f"phasing.{k}: {ph[k]} not found")
        _alignment_index(Path(ph["bam"]))
        if not Path(f"{ph['reference']}.fai").is_file():
            raise ConfigError(f"phasing.reference: index {ph['reference']}.fai not found; run `samtools faidx`")
    return ph


def _check_regions(regions):
    if not isinstance(regions, list):
        raise ConfigError(f"prioritise.exclude_regions must be a list, got {regions!r}")
    for i, r in enumerate(regions):
        if not isinstance(r, dict) or set(r) != {"name", "chrom", "start", "end"} \
                or not isinstance(r["name"], str) or not isinstance(r["chrom"], str) \
                or not _is_int(r["start"]) or not _is_int(r["end"]) or not 1 <= r["start"] <= r["end"]:
            raise ConfigError(f"prioritise.exclude_regions[{i}]: expected {{name, chrom (string), start <= end "
                              f"(positive integers)}}, got {r!r}")


def _check_tiers(tiers):
    allowed = (LEVELS, FITS, PHENOS)
    if not isinstance(tiers, list) or not tiers:
        raise ConfigError("prioritise.tiers must be a non-empty list of tiers")
    for i, t in enumerate(tiers):
        if not isinstance(t, dict) or set(t) != {"name", "epcr", "rules"}:
            raise ConfigError(f"prioritise.tiers[{i}]: bad tier {t!r} (keys name, epcr, rules)")
        if not isinstance(t["name"], str):
            raise ConfigError(f"prioritise.tiers[{i}]: name must be a string, got {t['name']!r}")
        _number(t["epcr"], f"prioritise.tiers {t['name']}.epcr")
        if not isinstance(t["rules"], list):
            raise ConfigError(f"prioritise.tiers {t['name']}: rules must be a list, got {t['rules']!r}")
        for rule in t["rules"]:
            if not isinstance(rule, list) or len(rule) != 3 or \
                    any(v != "*" and v not in ok for v, ok in zip(rule, allowed)):
                raise ConfigError(f"prioritise.tiers {t['name']}: bad rule {rule!r}; expected "
                                  f"[level {LEVELS}, fit {FITS}, phenotype {PHENOS}] or '*'")
    if ["*", "*", "*"] not in tiers[-1]["rules"]:
        raise ConfigError("prioritise.tiers: the last tier must contain the catch-all rule ['*', '*', '*']")


def _check_gene_panel(gp):
    genes = gp["genes"]
    if not isinstance(gp["name"], str):
        raise ConfigError(f"prioritise.gene_panel.name must be a string, got {gp['name']!r}")
    if not isinstance(genes, list) or not all(isinstance(g, str) for g in genes):
        raise ConfigError(f"prioritise.gene_panel.genes must be a list of gene symbols, got {genes!r}")
    if gp["phenotype_floor"] not in PANEL_FLOORS:
        raise ConfigError(f"prioritise.gene_panel.phenotype_floor must be one of {PANEL_FLOORS}, "
                          f"got {gp['phenotype_floor']!r}")
    _bool(gp["deep_second_hit"], "prioritise.gene_panel.deep_second_hit")
    if not genes and (gp["deep_second_hit"] or gp["phenotype_floor"]):
        raise ConfigError("prioritise.gene_panel: deep_second_hit and phenotype_floor need at least one gene")


def _check_alternative_pairs(ap):
    if not _is_int(ap["max_per_gene"]) or ap["max_per_gene"] < 1:
        raise ConfigError(f"prioritise.alternative_pairs.max_per_gene must be a positive integer, "
                          f"got {ap['max_per_gene']!r}")
    if ap["min_phenotype"] not in PHENOS:
        raise ConfigError(f"prioritise.alternative_pairs.min_phenotype must be one of {PHENOS}, "
                          f"got {ap['min_phenotype']!r}")


def prioritise_settings(cfg: dict, top_n: int | None = None, output: str | None = None) -> dict:
    s = _section(cfg, "prioritise", PRIORITISE_DEFAULTS)
    if top_n is not None:
        s["top_n"] = top_n
    if output is not None:
        s["output"] = output
    if not _is_int(s["top_n"]) or s["top_n"] < 1:
        raise ConfigError(f"prioritise.top_n must be a positive integer, got {s['top_n']!r}")
    sec = s["secondary"]
    _bool(sec["enabled"], "prioritise.secondary.enabled")
    if not _is_int(sec["max_rows"]) or not 0 <= sec["max_rows"] <= s["top_n"]:
        raise ConfigError(f"prioritise.secondary.max_rows must be an integer between 0 and "
                          f"top_n={s['top_n']}, got {sec['max_rows']!r}")
    for k, v in sec["epcr"].items():
        _number(v, f"prioritise.secondary.epcr.{k}")
    bad = [c for c in s["exclude_evidence_classes"] if c not in EXCLUDABLE]
    if bad:
        raise ConfigError(f"prioritise.exclude_evidence_classes: {bad} not in {EXCLUDABLE}")
    if s["missense_tool"] not in MISSENSE_TOOLS:
        raise ConfigError(f"prioritise.missense_tool must be one of {MISSENSE_TOOLS}, got {s['missense_tool']!r}")
    ps = s["phase_sources"]
    if not isinstance(ps, list) or any(x not in PHASE_SOURCES for x in ps):
        raise ConfigError(f"prioritise.phase_sources must be a list of {PHASE_SOURCES}, got {ps!r}")
    _number(s["p_trans_unphased"], "prioritise.p_trans_unphased", lo_open=True)
    _number(s["primary_sum_cap"], "prioritise.primary_sum_cap", lo_open=True)
    pops = s["af"]["populations"]
    if not isinstance(pops, list) or not pops:
        raise ConfigError(f"prioritise.af.populations must be a non-empty list, got {pops!r}")
    bad = [p for p in pops if p not in GNOMAD_GROUPS]
    if bad:
        raise ConfigError(f"prioritise.af.populations: {bad} not gnomAD groups {GNOMAD_GROUPS}")
    for k in ("max_dominant", "max_recessive"):
        _number(s["af"][k], f"prioritise.af.{k}", lo_open=True)
    pr = s["pheno_rank"]
    for k in ("strong", "moderate"):
        if not _is_int(pr[k]) or pr[k] < 1:
            raise ConfigError(f"prioritise.pheno_rank.{k} must be a positive integer, got {pr[k]!r}")
    if pr["strong"] > pr["moderate"]:
        raise ConfigError(f"prioritise.pheno_rank: strong ({pr['strong']}) must be <= moderate ({pr['moderate']})")
    _check_regions(s["exclude_regions"])
    _check_gene_panel(s["gene_panel"])
    _check_alternative_pairs(s["alternative_pairs"])
    _check_tiers(s["tiers"])
    s["output"] = s["output"] or str(Path(cfg["output"]["dir"]) / f"{cfg['output']['prefix']}.top{s['top_n']}.csv")
    return s
