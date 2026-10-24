from dataclasses import replace

from prioritise_helpers import make_ev

from vep_pipeline.evidence import clinvar_class, gene_evidence, lower, splice_shifted_out, variant_level
from vep_pipeline.genedb import GeneInfo

AR_LOF = GeneInfo(frozenset({"AR"}), frozenset({"loss of function"}))
AD_GOF = GeneInfo(frozenset({"AD"}), frozenset({"gain of function"}))
AD_UNK = GeneInfo(frozenset({"AD"}), frozenset())
POPS = ("AFR", "AMR", "EAS", "NFE", "SAS")


def test_clinvar_class():
    assert clinvar_class("Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts") == ("plp", 2)
    assert clinvar_class("Pathogenic/Likely_pathogenic", "criteria_provided&_single_submitter") == ("plp", 1)
    assert clinvar_class("Benign/Likely_benign", "reviewed_by_expert_panel") == ("blb", 3)
    assert clinvar_class("Conflicting_classifications_of_pathogenicity", "x") == (None, 0)
    assert clinvar_class("Uncertain_significance", "criteria_provided,_single_submitter") == (None, 1)
    assert clinvar_class("", "") == (None, 0)


def test_gene_evidence_uses_mane_and_core_populations():
    entries = [
        {"SYMBOL": "G", "Feature": "ENST2", "Consequence": "intron_variant", "REVEL": "0.2",
         "gnomADg_FIN_AF": "0.2", "gnomADg_NFE_AF": "0.001"},
        {"SYMBOL": "G", "Feature": "ENST1", "MANE_SELECT": "NM_1.1", "Consequence": "missense_variant",
         "HGVSp": "ENSP1:p.Arg1Cys", "REVEL": "0.95&0.90"},
        {"SYMBOL": "H", "Feature": "ENST9", "Consequence": "upstream_gene_variant"},
    ]
    out = gene_evidence(entries, POPS)
    g = out["G"]
    assert g.consequence == "missense_variant" and g.hgvs == "ENSP1:p.Arg1Cys"
    assert g.revel == 0.95 and g.popmax == 0.001          # FIN (bottlenecked) is ignored
    assert set(out) == {"G", "H"}


def test_revel_levels_and_one_tool_rule():
    assert variant_level(make_ev(revel=0.95), None)[0] == "strong"
    assert variant_level(make_ev(revel=0.80), None)[0] == "moderate"
    assert variant_level(make_ev(revel=0.70), None)[0] == "supporting"
    assert variant_level(make_ev(revel=0.50, am=0.99), None, missense_tool="revel")[0] is None   # REVEL decides
    assert variant_level(make_ev(am=0.95), None)[0] == "moderate"                # AM only when REVEL is blank
    assert variant_level(make_ev(revel=0.50, am=0.99), None, frozenset({"clinically_trained_predictor"}))[0] == "strong"


def test_null_variants_follow_pvs1_outline():
    stop = dict(consequence="stop_gained")
    assert variant_level(make_ev(**stop), AR_LOF)[:3:2] == ("strong", ["null:strong"])
    assert variant_level(make_ev(nmd_escape=True, **stop), AR_LOF)[0] == "moderate"
    assert variant_level(make_ev(**stop), AD_GOF)[0] == "supporting"
    non_lof = GeneInfo(frozenset({"AD"}), frozenset({"undetermined non-loss-of-function"}))
    assert variant_level(make_ev(loeuf=0.1, **stop), non_lof)[0] == "supporting"
    assert variant_level(make_ev(loeuf=0.3, **stop), AD_UNK)[0] == "strong"
    assert variant_level(make_ev(loeuf=0.9, **stop), AD_UNK)[0] == "moderate"
    assert variant_level(make_ev(loeuf=0.3, **stop), AD_UNK, frozenset({"gene_prior_dominant"}))[0] == "moderate"
    assert variant_level(make_ev(consequence="start_lost"), AR_LOF)[0] == "moderate"


def test_splice_and_cadd():
    assert variant_level(make_ev(consequence="intron_variant", splice=0.3), None)[0] == "supporting"
    assert variant_level(make_ev(consequence="intron_variant", splice=0.15), None)[0] is None
    assert variant_level(make_ev(consequence="intron_variant", cadd=26.0), None)[0] == "supporting"
    assert variant_level(make_ev(consequence="intron_variant", cadd=26.0, splice=0.05), None)[0] is None


def test_noncoding_cadd_is_second_hit_only_but_missense_cadd_stands_alone():
    intron = variant_level(make_ev(consequence="intron_variant", cadd=26.0), None)
    assert intron[:2] == ("supporting", True) and intron[2] == ["non-coding CADD 26.0 (second hit only)"]
    assert variant_level(make_ev(cadd=26.0), None)[:2] == ("supporting", False)          # missense, no REVEL/AM
    assert variant_level(make_ev(consequence="intron_variant", cadd=26.0), None,
                         frozenset({"unsupervised_predictor"}))[0] is None
    assert variant_level(make_ev(consequence="synonymous_variant", cadd=26.0), None)[:2] == ("supporting", True)
    assert variant_level(make_ev(consequence="inframe_deletion", cadd=26.0), None)[1] is False


def test_remm_is_a_second_hit_signal_for_non_protein_altering_variants():
    hit = variant_level(make_ev(consequence="intron_variant", remm=0.97), None)
    assert hit[:2] == ("supporting", True) and hit[2] == ["ReMM 0.970 (second hit only)"]
    assert variant_level(make_ev(consequence="intron_variant", remm=0.90), None)[0] is None
    assert variant_level(make_ev(consequence="intron_variant", remm=0.97), None,
                         frozenset({"clinically_trained_predictor"}))[0] is None
    assert variant_level(make_ev(remm=0.97), None)[0] is None                            # missense: not a ReMM target
    both = variant_level(make_ev(consequence="intron_variant", cadd=26.0, remm=0.97), None)
    assert both[:2] == ("supporting", True) and len(both[2]) == 2


UNSUP, CLIN = "unsupervised_predictor", "clinically_trained_predictor"


def test_grn_is_supporting_second_hit_only():
    ev = make_ev(consequence="intergenic_variant", grn="peak chr1:100-600 -> GENEA; R1", grn_classes=frozenset({UNSUP}))
    assert variant_level(ev, None) == ("supporting", True,
        ["regulatory network: peak chr1:100-600 -> GENEA; R1 (second hit only)"])
    assert variant_level(ev, None, frozenset({"regulatory_network"}))[0] is None


def test_grn_combines_with_other_second_hit_signals():
    ev = make_ev(consequence="intron_variant", remm=0.97, grn="peak p -> G; R2", grn_classes=frozenset({CLIN}))
    assert variant_level(ev, None)[2] == ["ReMM 0.970 (second hit only)",
                                          "regulatory network: peak p -> G; R2 (second hit only)"]
    assert variant_level(ev, None, frozenset({"regulatory_network"}))[2] == ["ReMM 0.970 (second hit only)"]


def test_grn_with_clinvar_blb_still_vetoed():
    ev = make_ev(consequence="intergenic_variant", clinvar="blb", stars=1, grn="peak p -> G; R1",
                 grn_classes=frozenset({UNSUP}))
    assert variant_level(ev, None) == (None, False, ["ClinVar B/LB 1*"])
    assert variant_level(ev, None, frozenset({"known_variant"}))[:2] == ("supporting", True)   # veto lifted
    assert variant_level(replace(ev, stars=0), None)[:2] == ("supporting", True)              # 0 stars: no veto


def test_grn_is_regated_by_the_predictor_classes_behind_v():
    remm_only = make_ev(consequence="intergenic_variant", grn="peak p -> G; R1", grn_classes=frozenset({CLIN}))
    assert variant_level(remm_only, None)[:2] == ("supporting", True)
    assert variant_level(remm_only, None, frozenset({UNSUP}))[:2] == ("supporting", True)
    assert variant_level(remm_only, None, frozenset({CLIN}))[0] is None
    cadd_only = replace(remm_only, grn_classes=frozenset({UNSUP}))
    assert variant_level(cadd_only, None, frozenset({CLIN}))[:2] == ("supporting", True)
    assert variant_level(cadd_only, None, frozenset({UNSUP}))[0] is None
    both = replace(remm_only, grn_classes=frozenset({UNSUP, CLIN}))
    assert variant_level(both, None, frozenset({UNSUP}))[:2] == ("supporting", True)
    assert variant_level(both, None, frozenset({CLIN}))[:2] == ("supporting", True)
    assert variant_level(both, None, frozenset({UNSUP, CLIN}))[0] is None
    assert variant_level(replace(remm_only, grn_classes=frozenset()), None)[0] is None     # no class behind V


def test_clinvar_evidence_and_ablation():
    plp = make_ev(revel=0.1, clinvar="plp", stars=1)
    assert variant_level(plp, None)[0] == "strong"
    assert variant_level(plp, None, frozenset({"known_variant"}))[0] is None
    assert variant_level(make_ev(revel=0.99, clinvar="blb", stars=2), None)[0] is None


def test_noncoding_signal_is_second_hit_only():
    assert variant_level(make_ev(consequence="5_prime_UTR_variant", noncoding_signal=True), None)[:2] == \
        ("supporting", True)


def test_lower():
    assert lower("strong") == "moderate" and lower("supporting") == "supporting"


def test_alphamissense_is_the_default_missense_tool():
    assert variant_level(make_ev(revel=0.472, am=0.923), None)[:3:2] == ("moderate", ["AlphaMissense 0.923"])
    assert variant_level(make_ev(revel=0.95, am=0.50), None)[0] is None          # AM decides; REVEL not scanned
    assert variant_level(make_ev(revel=0.95), None)[0] == "strong"               # REVEL only when AM has no score
    assert variant_level(make_ev(revel=0.95, am=0.99), None,
                         frozenset({"unsupervised_predictor"}))[0] == "strong"   # AM excluded -> REVEL
    assert variant_level(make_ev(revel=0.472, am=0.923), None, missense_tool="revel")[0] is None


def test_splice_call_shifted_out_by_hgvs():
    assert splice_shifted_out("ENST1.1:c.100+77_100+121del")
    assert splice_shifted_out("ENST1.1:c.-103+77_-103+121del")
    assert splice_shifted_out("ENST1.1:c.148-30_148-10del")
    assert not splice_shifted_out("ENST1.1:c.147+1_147+4del")     # touches +1,2
    assert not splice_shifted_out("ENST1.1:c.147_147+5del")       # includes the last exon base
    assert not splice_shifted_out("ENST1.1:c.148-2del")
    assert not splice_shifted_out("ENST1.1:c.147+3_148-3del")     # spans an exon
    assert not splice_shifted_out("ENSP1.1:p.Arg100Ter")
    assert not splice_shifted_out("")


def test_repeat_indel_at_a_splice_site_is_not_a_null():
    ev = make_ev(consequence="splice_donor_variant&splice_donor_5th_base_variant&coding_sequence_variant&intron_variant",
                 hgvsc="ENST1.1:c.100+77_100+121del", splice=0.81)
    assert variant_level(ev, AR_LOF)[:3:2] == ("supporting", ["SpliceAI 0.81"])
    assert variant_level(make_ev(consequence="splice_donor_variant", hgvsc="ENST1.1:c.147+1del"), AR_LOF)[0] == "strong"


def test_gene_evidence_keeps_the_representative_hgvsc():
    out = gene_evidence([{"SYMBOL": "G", "Feature": "E1", "MANE_SELECT": "NM_1", "Consequence": "intron_variant",
                          "HGVSc": "E1:c.100+77_100+121del"}], POPS)
    assert out["G"].hgvsc == "E1:c.100+77_100+121del"


def test_panel_candidate_needs_the_panel_flag_and_no_calibrated_evidence():
    mis = make_ev(am=0.5)
    assert variant_level(mis, None)[0] is None                                       # outside the panel
    assert variant_level(mis, None, panel=True) == ("supporting", True, ["panel candidate (protein_altering)"])
    assert variant_level(make_ev(am=0.95), None, panel=True)[1:] == (False, ["AlphaMissense 0.950"])   # calibrated wins


def test_panel_candidate_classes_and_exclusions():
    def cls(**kw):
        return variant_level(make_ev(**kw), None, panel=True)
    intron = dict(consequence="intron_variant")
    assert cls(splice=0.15, **intron) == ("supporting", True, ["panel candidate (splice)"])
    assert cls(consequence="splice_region_variant&intron_variant")[2] == ["panel candidate (splice)"]
    assert cls(consequence="splice_polypyrimidine_tract_variant&intron_variant")[2] == ["panel candidate (splice)"]
    assert cls(consequence="missense_variant", am=0.5, splice=0.15)[2] == ["panel candidate (splice)"]   # splice first
    assert cls(consequence="5_prime_UTR_variant")[2] == ["panel candidate (5utr)"]
    assert cls(consequence="3_prime_UTR_variant")[2] == ["panel candidate (3utr)"]
    assert cls(consequence="synonymous_variant")[0] is None
    assert cls(consequence="synonymous_variant", splice=0.05)[0] is None             # SpliceAI <= 0.1: benign range
    assert cls(splice=0.1, **intron)[0] is None                                       # 0.1 itself is BP4
    assert cls(splice=0.11, **intron)[2] == ["panel candidate (splice)"]              # > SPLICE_CANDIDATE
    for term in ("upstream_gene_variant", "downstream_gene_variant", "intron_variant"):
        assert cls(consequence=term)[0] is None
    assert variant_level(make_ev(splice=0.15, **intron), None, frozenset({"unsupervised_predictor"}),
                         panel=True)[0] is None                                       # SpliceAI path honours exclusion
    assert variant_level(make_ev(consequence="splice_region_variant", splice=0.15), None,
                         frozenset({"unsupervised_predictor"}), panel=True)[2] == ["panel candidate (splice)"]


def test_spliceai_bp4_vetoes_a_splice_region_term():
    """Walker 2023: SpliceAI <= 0.1 is BP4. A usable score decides the splice class (> 0.1); the splice-region term
    alone counts only when there is no usable score (none, or unsupervised_predictor excluded). User ruling
    2026-10-20 (final review)."""
    def cls(exclude=frozenset(), **kw):
        return variant_level(make_ev(**kw), None, exclude, panel=True)[2]
    region = dict(consequence="splice_region_variant&intron_variant")
    assert cls(splice=0.05, **region) == []                                             # BP4: no candidate
    assert cls(**region) == ["panel candidate (splice)"]                                # no SpliceAI score
    assert cls(splice=0.1, **region) == []                                              # 0.1 is in the BP4 range
    assert cls(splice=0.15, **region) == ["panel candidate (splice)"]
    assert cls(consequence="missense_variant&splice_region_variant", splice=0.02) == [
        "panel candidate (protein_altering)"]                                           # falls through to its codon
    assert cls(frozenset({"unsupervised_predictor"}), splice=0.05, **region) == ["panel candidate (splice)"]
    assert cls(consequence="splice_region_variant&3_prime_UTR_variant", splice=0.05) == ["panel candidate (3utr)"]


def test_repeat_indel_with_a_low_spliceai_is_no_panel_candidate():
    ev = make_ev(consequence="splice_donor_variant&splice_donor_5th_base_variant&coding_sequence_variant&intron_variant",
                 hgvsc="ENST1.1:c.100+77_100+121del", splice=0.02)
    assert variant_level(ev, AR_LOF, panel=True) == (None, False, [])
    assert variant_level(make_ev(consequence=ev.consequence, hgvsc=ev.hgvsc), AR_LOF, panel=True)[2] == [
        "panel candidate (splice)"]                                                     # no SpliceAI: the term counts


def test_panel_candidate_yields_to_clinvar_benign_and_other_signals():
    benign = make_ev(am=0.5, clinvar="blb", stars=1)
    assert variant_level(benign, None, panel=True) == (None, False, ["ClinVar B/LB 1*"])
    assert variant_level(make_ev(am=0.5, clinvar="blb", stars=0), None, panel=True)[2] == [
        "panel candidate (protein_altering)"]
    remm = make_ev(consequence="intron_variant", remm=0.97)
    assert variant_level(remm, None, panel=True)[2] == ["ReMM 0.970 (second hit only)"]    # a second-hit signal wins
