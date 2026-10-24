"""Per-variant, per-gene evidence from VEP CSQ entries, and the variant evidence level.

Thresholds (verified against the source tables, 2026-10-17):
  REVEL PP3 supporting [0.644,0.773), moderate [0.773,0.932), strong >=0.932 (Pejaver 2022 Table 2,
    doi:10.1016/j.ajhg.2022.10.013).
  AlphaMissense PP3 supporting [0.792,0.906), moderate [0.906,0.990), strong >=0.990 (Bergquist 2025 Table 1,
    doi:10.1016/j.gim.2025.101402). One missense tool per variant (Pejaver 2022): `missense_tool` (default
    alphamissense) decides, and the other tool is used only when the chosen one has no score or its evidence class is
    excluded. Model 2 used revel; Model 3 (decided after Model 2's score) uses alphamissense.
  SpliceAI / SpliceAI2 delta >=0.2, outside the +-1,2 dinucleotides, supporting weight only (Walker 2023,
    doi:10.1016/j.ajhg.2023.06.002). SpliceAI2 has no ClinGen calibration; the SpliceAI v1 cut-off is reused.
  CADD PHRED >=25.3 supporting (Pejaver 2022; a missense calibration), only when no other tool scored the variant.
    Standalone only for protein-altering consequences (PROTEIN_ALTERING). For every other consequence it is a
    second-hit-only signal: missense predictors "cannot be applied to variants in non-coding regions" and such
    variants are interpreted only in phenotype-matched genes, as an in-trans second hit (Ellingford 2022 Genome Med,
    doi:10.1186/s13073-022-01073-3), and the missense cut-off does not transfer: the promoter calibration gives
    CADD PP3 at >= 10 (Villani 2024 AJHG, doi:10.1016/j.ajhg.2024.05.002).
  Genomiser ReMM >= 0.963 (REMM_SECOND_HIT) on a non-protein-altering variant is a second-hit-only signal. ReMM is
    from Smedley 2016 AJHG (doi:10.1016/j.ajhg.2016.07.005); the cut-off is the one Xi 2025 Brief Bioinform
    (doi:10.1093/bib/bbaf475) used as a filter on 985 WGS (best as a filter, not as evidence). ReMM is trained on
    curated pathogenic regulatory variants, so it is excluded with the clinically_trained_predictor class.
  A rare non-coding variant in an scDoRI peak linked to a strong-phenotype gene, with CADD >= 20 or ReMM >= 0.963
    (regnet second_hit rows), is a second-hit-only signal for that gene (Evidence.grn); it is excluded with the
    regulatory_network class, and also when every predictor class behind V (Evidence.grn_classes) is excluded.
  Null variants follow the PVS1 tree in outline (Abou Tayoun 2018, doi:10.1002/humu.23626): start_lost and
    NMD-escaping are moderate; a G2P mechanism other than loss of function gives supporting; a dominant or
    unknown-MOI gene without a G2P LoF mechanism needs LOEUF < 0.6 or pHaplo >= 0.86 (Collins 2022,
    doi:10.1016/j.cell.2022.06.036) for strong.
Evidence levels do not add up: each evidence type is scored once and the best level wins.
Canonical splice terms are ignored when the HGVS-shifted indel lies wholly in one intron > 2 bp from the exon
(splice_shifted_out): a tandem-repeat indel at a splice site is not a PVS1 null.
"""
from __future__ import annotations
import re
from dataclasses import dataclass

LEVELS = ("supporting", "moderate", "strong")
REVEL_PP3 = ((0.932, "strong"), (0.773, "moderate"), (0.644, "supporting"))
AM_PP3 = ((0.990, "strong"), (0.906, "moderate"), (0.792, "supporting"))
MISSENSE_TOOLS = ("alphamissense", "revel")
SPLICE_PP3 = 0.2
SPLICE_CANDIDATE = 0.1    # Walker 2023 (doi:10.1016/j.ajhg.2023.06.002): <= 0.1 is BP4; a candidate needs > 0.1
CADD_PP3 = 25.3
CADD_NONCODING = 20.0
REMM_SECOND_HIT = 0.963   # Xi 2025 Brief Bioinform, doi:10.1093/bib/bbaf475 (ReMM filter cut-off on 985 WGS)
PROTEIN_ALTERING = {"missense_variant", "protein_altering_variant", "inframe_insertion", "inframe_deletion",
                    "stop_lost", "start_lost"}
LOEUF_MAX, PHAPLO_MIN = 0.6, 0.86
LOF_MECH, UNDETERMINED = "loss of function", "undetermined"
NULL_SO = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant",
           "transcript_ablation"}
CANONICAL_SPLICE = {"splice_acceptor_variant", "splice_donor_variant"}
INFRAME = {"inframe_insertion", "inframe_deletion"}
SPLICE_REGION = {"splice_region_variant", "splice_donor_5th_base_variant", "splice_donor_region_variant",
                 "splice_polypyrimidine_tract_variant"}
UTR5, UTR3 = {"5_prime_UTR_variant"}, {"3_prime_UTR_variant"}
PANEL_CANDIDATE_TERMS = PROTEIN_ALTERING | INFRAME | SPLICE_REGION | UTR5 | UTR3
PANEL_CLASSES = ("splice", "protein_altering", "5utr", "3utr")   # priority order; the first match names the class
PANEL_REASON = "panel candidate"
SPLICE_FIELDS = ("SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL",
                 "SpliceAI2_spliceai2_summary_score")
_STARS = (("practice_guideline", 4), ("reviewed_by_expert_panel", 3),
          ("criteria_provided,_multiple_submitters,_no_conflicts", 2), ("criteria_provided,_single_submitter", 1))
_PLP = {"pathogenic", "likely_pathogenic"}
_BLB = {"benign", "likely_benign"}


@dataclass(frozen=True)
class Evidence:
    symbol: str
    consequence: str          # '&'-joined SO terms of the representative transcript
    hgvs: str                 # HGVSp, else HGVSc, else ''
    popmax: float             # max AF over the configured gnomAD populations; 0.0 when absent
    clinvar: str | None       # 'plp', 'blb' or None
    stars: int
    revel: float | None
    am: float | None
    splice: float | None      # max SpliceAI / SpliceAI2 delta score
    cadd: float | None
    nmd_escape: bool
    loeuf: float | None
    phaplo: float | None
    noncoding_signal: bool    # 5UTR_annotation set, or ENCODE cCRE overlap with CADD >= 20: second-hit candidate only
    remm: float | None = None  # Genomiser ReMM score, attached from the Exomiser variant TSV by candidates.collect
    hgvsc: str = ""            # HGVSc of the representative transcript (3'-shifted by VEP)
    grn: str = ""              # regulatory-network link (regnet.add_second_hits): a second-hit-only signal
    grn_classes: frozenset = frozenset()   # predictor classes behind the link's V (CADD: unsupervised_predictor,
                                           # ReMM: clinically_trained_predictor); grn needs one not excluded


def _f(v):
    """Float value of a CSQ field; '&'-joined values give their max; blanks give None."""
    if v in (None, "", ".", "-"):
        return None
    vals = []
    for part in str(v).split("&"):
        try:
            vals.append(float(part))
        except ValueError:
            pass
    return max(vals) if vals else None


def _max(entries, fields):
    vals = [x for e in entries for f in fields if (x := _f(e.get(f))) is not None]
    return max(vals) if vals else None


def _representative(entries):
    """MANE Select, else canonical, else PICK, else the lowest Feature id: one transcript per gene."""
    return min(entries, key=lambda e: (not e.get("MANE_SELECT"), e.get("CANONICAL") != "YES", e.get("PICK") != "1",
                                       e.get("Feature", "")))


_HGVS_INTRONIC = re.compile(r"^(?P<b1>[-*]?\d+)(?P<o1>[+-]\d+)(?:_(?P<b2>[-*]?\d+)(?P<o2>[+-]\d+))?(?:del|dup|ins|inv)")


def splice_shifted_out(hgvsc: str) -> bool:
    """True when the HGVSc (3'-shifted) puts an indel wholly inside one intron, more than 2 bp from the exon.

    VEP's consequence uses the left-aligned VCF position; for an indel in a tandem repeat spanning a splice site the
    shifted, equivalent allele can leave the +-1,2 dinucleotide intact (e.g. a 45-bp repeat-unit deletion called
    splice_donor_variant whose HGVSc is c.147+77_147+121del)."""
    if ":c." not in hgvsc:
        return False
    m = _HGVS_INTRONIC.match(hgvsc.split(":c.", 1)[1])
    if not m:
        return False
    offsets = [int(m["o1"])] + ([int(m["o2"])] if m["o2"] else [])
    same_intron = m["b2"] is None or (m["b2"] == m["b1"] and (offsets[0] > 0) == (offsets[1] > 0))
    return same_intron and min(abs(o) for o in offsets) > 2


def effective_terms(ev) -> set:
    """SO terms of the representative transcript, without canonical splice terms that the HGVS shift moves out."""
    terms = set(ev.consequence.split("&"))
    if terms & CANONICAL_SPLICE and splice_shifted_out(ev.hgvsc):
        terms -= CANONICAL_SPLICE
    return terms


def clinvar_class(sig: str, revstat: str) -> tuple:
    """('plp'|'blb'|None, review stars). Conflicting or mixed classifications give None. VEP may write ',' as '&'."""
    sig, revstat = (sig or "").lower(), (revstat or "").lower().replace("&", ",")
    if not sig or "conflicting" in sig or "conflicting" in revstat:
        return None, 0
    tokens = set(re.split(r"[/&,|]", sig))
    stars = next((n for key, n in _STARS if key in revstat), 0)
    if tokens & _PLP and not tokens & _BLB:
        return "plp", stars
    if tokens & _BLB and not tokens & _PLP:
        return "blb", stars
    return None, stars


def gene_evidence(entries: list, pops: tuple) -> dict:
    """{SYMBOL: Evidence} for one variant; allele-level scores are the max over that gene's CSQ entries."""
    af_fields = tuple(f"gnomAD{s}_{p}_AF" for s in ("e", "g") for p in pops)
    by_gene = {}
    for e in entries:
        if e.get("SYMBOL"):
            by_gene.setdefault(e["SYMBOL"], []).append(e)
    out = {}
    for symbol, es in by_gene.items():
        rep = _representative(es)
        cv = next((e for e in es if e.get("ClinVar_CLNSIG")), {})
        clin, stars = clinvar_class(cv.get("ClinVar_CLNSIG", ""), cv.get("ClinVar_CLNREVSTAT", ""))
        cadd = _max(es, ("CADD_PHRED",))
        out[symbol] = Evidence(
            symbol=symbol, consequence=rep.get("Consequence", ""),
            hgvs=rep.get("HGVSp") or rep.get("HGVSc") or "", hgvsc=rep.get("HGVSc", ""),
            popmax=_max(es, af_fields) or 0.0, clinvar=clin, stars=stars,
            revel=_max(es, ("REVEL",)), am=_max(es, ("am_pathogenicity",)), splice=_max(es, SPLICE_FIELDS),
            cadd=cadd, nmd_escape=rep.get("NMD") == "NMD_escaping_variant",
            loeuf=_max(es, ("LOEUF",)), phaplo=_max(es, ("pHaplo",)),
            noncoding_signal=bool(rep.get("5UTR_annotation"))
            or (bool(rep.get("ENCODE_cCRE")) and (cadd or 0.0) >= CADD_NONCODING))
    return out


def lower(level: str) -> str:
    """One level down; 'supporting' stays 'supporting'."""
    return LEVELS[max(LEVELS.index(level) - 1, 0)]


def _from(thresholds, score):
    return next((lvl for cut, lvl in thresholds if score >= cut), None)


def null_level(ev: Evidence, gene, exclude=frozenset()) -> str | None:
    terms = effective_terms(ev)
    if not terms & (NULL_SO | {"start_lost"}):
        return None
    if not terms & NULL_SO or ev.nmd_escape:
        return "moderate"                                    # start_lost, or escapes NMD
    mechs = gene.mechanisms if gene else frozenset()
    lof = LOF_MECH in mechs
    if not lof and any(m not in (LOF_MECH, UNDETERMINED) for m in mechs):
        return "supporting"                                  # G2P says the mechanism is not loss of function
    if lof or (gene is not None and gene.modes and gene.modes <= {"AR", "XR"}):
        return "strong"
    if "gene_prior_dominant" in exclude:
        return "moderate"
    constrained = ((ev.loeuf is not None and ev.loeuf < LOEUF_MAX)
                   or (ev.phaplo is not None and ev.phaplo >= PHAPLO_MIN))
    return "strong" if constrained else "moderate"


def _missense(ev: Evidence, exclude, tool: str) -> dict:
    """{reason: level} from one missense tool: `tool` first; the other only when the first has no usable score."""
    am = ev.am if "unsupervised_predictor" not in exclude else None
    revel = ev.revel if "clinically_trained_predictor" not in exclude else None
    order = [("AlphaMissense", am, AM_PP3), ("REVEL", revel, REVEL_PP3)]
    if tool == "revel":
        order.reverse()
    for name, score, cuts in order:
        if score is not None:
            lvl = _from(cuts, score)
            return {f"{name} {score:.3f}": lvl} if lvl else {}
    return {}


def panel_class(ev: Evidence, terms: set, unsupervised: bool) -> str | None:
    """Panel candidate class of a variant (see PANEL_CLASSES), or None when it is not a candidate.

    `splice`: a usable SpliceAI score (present, unsupervised_predictor not excluded) > SPLICE_CANDIDATE, or a
    splice-region term when there is no usable score. A usable score <= SPLICE_CANDIDATE is BP4 and vetoes the term
    (Walker 2023; user ruling 2026-10-20); the variant may still be protein_altering or UTR."""
    score = ev.splice if unsupervised else None
    splice = score > SPLICE_CANDIDATE if score is not None else bool(terms & SPLICE_REGION)
    if splice:
        return "splice"
    if terms & (PROTEIN_ALTERING | INFRAME):
        return "protein_altering"
    return "5utr" if terms & UTR5 else "3utr" if terms & UTR3 else None


def variant_level(ev: Evidence, gene, exclude=frozenset(), missense_tool="alphamissense", panel=False) -> tuple:
    """(level or None, second_hit_only, reasons) for one variant in one gene.

    None means no usable evidence, or ClinVar benign/likely benign with >= 1 star. With `panel` (the gene is in the
    configured gene panel), a variant with no other evidence that is a splice, protein-altering or UTR candidate is a
    supporting second-hit-only "panel candidate"."""
    known = "known_variant" not in exclude
    unsupervised = "unsupervised_predictor" not in exclude
    if known and ev.clinvar == "blb" and ev.stars >= 1:
        return None, False, [f"ClinVar B/LB {ev.stars}*"]
    terms = effective_terms(ev)
    found = {}
    if known and ev.clinvar == "plp" and ev.stars >= 1:
        found[f"ClinVar P/LP {ev.stars}*"] = "strong"
    nul = null_level(ev, gene, exclude)
    if nul:
        found[f"null:{nul}"] = nul
    if "missense_variant" in terms:
        found.update(_missense(ev, exclude, missense_tool))
    if ev.splice is not None and ev.splice >= SPLICE_PP3 and not terms & CANONICAL_SPLICE and unsupervised:
        found[f"SpliceAI {ev.splice:.2f}"] = "supporting"
    if terms & INFRAME:
        found["in-frame indel"] = "supporting"
    protein_altering = bool(terms & PROTEIN_ALTERING)
    cadd_only = (ev.cadd is not None and ev.cadd >= CADD_PP3 and unsupervised
                 and ev.revel is None and ev.am is None and ev.splice is None)
    if cadd_only and protein_altering:
        found[f"CADD {ev.cadd:.1f}"] = "supporting"
    if found:
        return max(found.values(), key=LEVELS.index), False, sorted(found)
    second = []                                              # second-hit-only signals; they may combine
    if ev.noncoding_signal and unsupervised:
        second.append("non-coding signal (5'UTR uORF, or cCRE with CADD >= 20)")
    if cadd_only and not protein_altering:
        second.append(f"non-coding CADD {ev.cadd:.1f} (second hit only)")
    if (ev.remm is not None and ev.remm >= REMM_SECOND_HIT and not protein_altering
            and "clinically_trained_predictor" not in exclude):
        second.append(f"ReMM {ev.remm:.3f} (second hit only)")
    if ev.grn and "regulatory_network" not in exclude and ev.grn_classes - exclude:
        second.append(f"regulatory network: {ev.grn} (second hit only)")
    if second:
        return "supporting", True, second
    cls = panel_class(ev, terms, unsupervised) if panel else None
    if cls:
        return "supporting", True, [f"{PANEL_REASON} ({cls})"]
    return None, False, []
