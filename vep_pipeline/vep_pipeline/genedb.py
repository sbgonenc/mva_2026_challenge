"""Gene-level knowledge read from files: G2P (allelic requirement, molecular mechanism) and the ACMG SF list.

G2P is read from the downloaded all-panels CSV rather than the G2P plugin's per-variant `G2P_gene_req`, which is
only set for variants that pass the plugin's own frequency and consequence filters.
"""
from __future__ import annotations
import csv
from dataclasses import dataclass
from pathlib import Path

from .config import ConfigError

G2P_COLUMNS = ("gene symbol", "allelic requirement", "confidence", "molecular mechanism")
SF_RULES = ("all_plp", "truncating_only", "biallelic", "hfe_c282y_hom", "xl_any", "xl_biallelic")


@dataclass(frozen=True)
class GeneInfo:
    modes: frozenset          # subset of {"AD", "AR", "XD", "XR", "MT"}
    mechanisms: frozenset     # lower-case G2P molecular mechanism strings


def moi_codes(allelic_requirement: str) -> set:
    """G2P allelic requirement -> MOI codes; several requirements may be ';'- or ','-joined."""
    out = set()
    for part in allelic_requirement.lower().replace(",", ";").split(";"):
        p = part.strip()
        if not p:
            continue
        if "mitochondrial" in p:
            out.add("MT")
        elif "x_hemizygous" in p or p == "hemizygous":
            out.add("XR")
        elif "x_heterozygous" in p:
            out.add("XD")
        elif p == "monoallelic_x":                         # X-linked, unspecified: both hypotheses
            out.update({"XD", "XR"})
        elif "biallelic" in p:
            out.add("AR")
        elif "monoallelic" in p:
            out.add("AD")
    return out


def load_g2p(path, confidences) -> dict:
    """{gene symbol: GeneInfo} from the G2P all-panels CSV, keeping rows whose confidence is in `confidences`."""
    path = Path(path)
    if not path.is_file():
        raise ConfigError(f"G2P file not found: {path} (run `vep-pipeline download --only G2P`)")
    keep = {c.lower() for c in confidences}
    modes, mechs = {}, {}
    with open(path, newline="") as fh:
        rows = csv.DictReader(fh)
        cols = {c.lower().strip(): c for c in rows.fieldnames or ()}
        missing = [c for c in G2P_COLUMNS if c not in cols]
        if missing:
            raise ConfigError(f"{path}: G2P columns {missing} missing; header has {sorted(cols)}")
        for r in rows:
            if r[cols["confidence"]].strip().lower() not in keep:
                continue
            gene = r[cols["gene symbol"]].strip()
            modes.setdefault(gene, set()).update(moi_codes(r[cols["allelic requirement"]]))
            mech = r[cols["molecular mechanism"]].strip().lower()
            mechs.setdefault(gene, set()).update({mech} if mech else set())
    return {g: GeneInfo(frozenset(modes[g]), frozenset(mechs.get(g, ()))) for g in modes}


def g2p_path(cfg, settings) -> Path:
    return Path(settings["g2p"]["file"] or Path(cfg["data"]["plugin_data_dir"]) / "G2P" / "G2P_all.csv")


@dataclass(frozen=True)
class SfGene:
    gene: str
    inheritance: str          # AD, AR, XL or SD as printed in the SF table
    rule: str                 # one of SF_RULES


def load_sf(path) -> dict:
    """{gene: SfGene} from the shipped ACMG SF table (gene, inheritance, rule; '#' lines are comments)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            gene, inheritance, rule = line.rstrip("\n").split("\t")
            if gene == "gene":
                continue
            if rule not in SF_RULES:
                raise ConfigError(f"{path}: unknown rule {rule!r} for {gene}")
            out[gene] = SfGene(gene, inheritance, rule)
    return out
