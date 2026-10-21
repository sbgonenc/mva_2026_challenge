from dataclasses import replace

import pytest
from prioritise_helpers import make_call, make_cand

from vep_pipeline import prioritise as pr
from vep_pipeline.candidates import build_alleles, findings_for_gene
from vep_pipeline.exomiser import PhenoEvidence
from vep_pipeline.genedb import GeneInfo, SfGene
from vep_pipeline.rank_settings import DEFAULT_TIERS, PRIORITISE_DEFAULTS

S = {**PRIORITISE_DEFAULTS, "exclude_evidence_classes": []}
AR = GeneInfo(frozenset({"AR"}), frozenset({"loss of function"}))
AD = GeneInfo(frozenset({"AD"}), frozenset({"loss of function"}))


@pytest.mark.parametrize("args, tier", [
    (("strong", "full", "strong"), ("T1", 0.90)), (("strong", "unknown_moi", "strong"), ("T2", 0.60)),
    (("supporting", "full", "strong"), ("T3", 0.30)), (("strong", "single_hit", "strong"), ("T4", 0.10)),
    (("strong", "full", "none"), ("T4", 0.10)), (("supporting", "single_hit", "weak"), ("T5", 0.02)),
])
def test_tier_for(args, tier):
    assert pr.tier_for(*args, DEFAULT_TIERS) == tier


def _findings(cands, info, sex="MALE"):
    al = build_alleles(cands, {"G": info} if info else {}, S, sex)
    return findings_for_gene("G", al["G"], info, S)


def test_unphased_pair_is_discounted_trans_pair_is_not():
    pair = _findings([make_cand(1, revel=0.95), make_cand(9, revel=0.95)], AR)[-1]
    assert pr.primary_epcr(pair, DEFAULT_TIERS, "strong", 0.5) == ("T1", 0.45)
    t2 = _findings([make_cand(8000, call=make_call(8000, alt="C", orig="1|8000|A|C,G|1"), revel=0.95),
                    make_cand(8000, call=make_call(8000, alt="G", orig="1|8000|A|C,G|2"), revel=0.95)], AR)[-1]
    assert t2.phase == "trans" and pr.primary_epcr(t2, DEFAULT_TIERS, "strong", 0.5) == ("T1", 0.90)


def test_cap_primary_scales_without_mutating():
    rows = [pr.Scored(None, "primary", "T1", 0.9, "strong"), pr.Scored(None, "primary", "T2", 0.6, "strong")]
    capped = pr.cap_primary(rows, 0.95)
    assert sum(s.epcr for s in capped) == pytest.approx(0.95) and rows[0].epcr == 0.9
    assert capped[0].epcr > capped[1].epcr


def test_select_puts_at_most_two_secondaries_last():
    prim = [pr.Scored(None, "primary", "T5", 0.02, "weak") for _ in range(12)]
    sec = [pr.Scored(None, "secondary", "SF", 0.3, "weak") for _ in range(3)]
    chosen = pr.select(prim, sec, {**S, "primary_sum_cap": 0.95})
    assert [s.finding_type for s in chosen] == ["primary"] * 8 + ["secondary"] * 2


def test_secondary_rules():
    tp53 = {"TP53": SfGene("TP53", "AD", "all_plp")}
    f = _sf_finding("TP53", AD, clinvar="plp", stars=2)
    pe = PhenoEvidence({"TP53": 500}, {}, {})
    prim, sec = pr.score_all([f], pe, S, tp53)
    assert not prim and sec[0].epcr == 0.30 and sec[0].note.startswith("Incidental finding unrelated")
    pe_match = PhenoEvidence({"TP53": 3}, {}, {})
    prim, sec = pr.score_all([f], pe_match, S, tp53)
    assert prim and not sec                                    # phenotype match: scored as primary
    ttn = SfGene("TTN", "AD", "truncating_only")
    assert pr.secondary_for(_sf_finding("TTN", AD, clinvar="plp", stars=2), ttn, frozenset()) is None
    hfe = SfGene("HFE", "AR", "hfe_c282y_hom")
    hom = _sf_finding("HFE", AR, gt="1/1", clinvar="plp", stars=3, hgvs="ENSP:p.Cys282Tyr")
    assert pr.secondary_for(hom, hfe, frozenset()) == "clinvar_2star"
    assert pr.secondary_for(_sf_finding("TP53", AD, clinvar="plp", stars=2), SfGene("TP53", "AD", "all_plp"),
                            frozenset({"known_variant"})) is None


def _sf_finding(gene, info, gt="0/1", **ev):
    c = make_cand(1, symbol=gene, call=make_call(1, gt=gt), revel=0.1, **ev)
    al = build_alleles([c], {gene: info}, S, "MALE")[gene]
    return findings_for_gene(gene, al, info, S)[0]


def test_coords_use_the_original_record():
    al = build_alleles([make_cand(5000, call=make_call(5000, orig="1|5000|AC|GT|1"), revel=0.95),
                        make_cand(5001, call=make_call(5001, orig="1|5000|AC|GT|1"), revel=0.95)], {}, S, "MALE")["G"]
    assert pr.coords(al[0], "chr") == ("chr1", 5000, "AC", "GT")
    multi = build_alleles([make_cand(8, call=make_call(8, alt="G", orig="1|8|A|C,G|2"), revel=0.95)], {}, S,
                          "MALE")["G"]
    assert pr.coords(multi[0], "chr") == ("chr1", 8, "A", "G")
    assert pr.changed(al[0]) and pr.changed(multi[0])
    same = build_alleles([make_cand(7, call=make_call(7, orig="1|7|A|G|1"), revel=0.95)], {}, S, "MALE")["G"]
    assert not pr.changed(same[0])                                  # EXO_ORIG is on every record
    mt = build_alleles([make_cand(3, call=make_call(3, chrom="MT"), revel=0.95)], {}, S, "MALE")["G"]
    assert pr.coords(mt[0], "chr")[0] == "chrM" and pr.coords(mt[0], "")[0] == "MT"


def test_csv_row_single_has_empty_slot_two():
    f = _findings([make_cand(1, revel=0.95)], AD)[0]
    row = pr.csv_row(pr.Scored(f, "primary", "T1", 0.9, "strong"), "S1", "chr", None)
    assert len(row) == len(pr.CSV_HEADER) and row[5:9] == ["", "", "", ""] and row[9] == "0.900"


def test_best_per_gene_and_variant_dedupe():
    a = pr.Scored(_findings([make_cand(1, revel=0.95)], AD)[0], "primary", "T2", 0.6, "moderate")
    b = pr.Scored(_findings([make_cand(1, revel=0.95)], AD)[0], "primary", "T1", 0.9, "strong")
    assert pr.best_per_gene([a, b], None) == [b]


# ----------------------------------------------------------------------------- gene panel

def _panel_settings(floor="strong", deep=True, genes=("G",)):
    return {**S, "gene_panel": {"name": "P", "genes": list(genes), "phenotype_floor": floor, "deep_second_hit": deep}}


def test_panel_partner_pair_is_tier_t4_and_halved_when_unphased():
    ps = _panel_settings()
    al = build_alleles([make_cand(1, consequence="stop_gained"), make_cand(9, am=0.5)], {"G": AR}, ps, "MALE")["G"]
    pair = [f for f in findings_for_gene("G", al, AR, ps) if len(al) == 2 and len(f.alleles) == 2][0]
    assert pair.fit == "panel_partner"
    assert pr.primary_epcr(pair, DEFAULT_TIERS, "strong", 0.5) == ("T4", 0.05)
    assert pr.primary_epcr(pair, DEFAULT_TIERS, "strong", 1.0) == ("T4", 0.10)
    for weaker in ("moderate", "weak", "unknown", "none"):                                    # only a strong phenotype
        assert pr.tier_for("strong", "panel_partner", weaker, DEFAULT_TIERS) == ("T5", 0.02)


def test_panel_partner_pair_never_outranks_a_calibrated_pair_at_a_weak_phenotype():
    ps = _panel_settings(floor=None)
    cands = [make_cand(1, consequence="stop_gained"), make_cand(9, am=0.80, popmax=0.001),
             make_cand(20, consequence="3_prime_UTR_variant")]
    al = build_alleles(cands, {"G": AR}, ps, "MALE")["G"]
    pairs = [f for f in findings_for_gene("G", al, AR, ps) if len(f.alleles) == 2]
    assert {f.fit for f in pairs} == {"full", "panel_partner"}
    pe = PhenoEvidence({"G": 500}, {}, {})
    primary, _ = pr.score_all(pairs, pe, ps, {})
    assert {s.pheno for s in primary} == {"weak"} and {s.tier for s in primary} == {"T5"}
    assert len({s.epcr for s in primary}) == 1                                   # a tie on epcr, level and Exomiser
    ordered = sorted(primary, key=lambda s: pr._sort_key(s, pe))
    assert [s.finding.fit for s in ordered] == ["full", "panel_partner", "panel_partner"]
    assert [s.finding.fit for s in pr.best_per_gene(primary, pe)] == ["full"]

def _single(gene="G"):
    c = make_cand(1, symbol=gene, revel=0.95)
    al = build_alleles([c], {gene: AD}, S, "MALE")[gene]
    return findings_for_gene(gene, al, AD, S)[0]


def test_phenotype_floor_raises_a_panel_gene_only():
    pe = PhenoEvidence({"G": 500, "H": 500}, {}, {})
    (row,), _ = pr.score_all([_single("G")], pe, _panel_settings(), {})
    assert row.pheno == "strong" and row.tier == "T1" and row.panel == "P (phenotype floor)"
    (outside,), _ = pr.score_all([_single("H")], pe, _panel_settings(), {})
    assert outside.pheno == "weak" and outside.panel == ""
    (off,), _ = pr.score_all([_single("G")], pe, _panel_settings(floor=None), {})
    assert off.pheno == "weak" and off.panel == "P"                                  # panel gene, no floor
    (default,), _ = pr.score_all([_single("G")], pe, S, {})
    assert default.pheno == "weak" and default.panel == ""


def test_phenotype_floor_never_lowers():
    pe = PhenoEvidence({"G": 3}, {}, {})
    (row,), _ = pr.score_all([_single("G")], pe, _panel_settings(floor="moderate"), {})
    assert row.pheno == "strong" and row.panel == "P"                                # already stronger: not raised
    pe_far = PhenoEvidence({"G": 500}, {}, {})
    (row,), _ = pr.score_all([_single("G")], pe_far, _panel_settings(floor="moderate"), {})
    assert row.pheno == "moderate" and row.panel == "P (phenotype floor)"
    (row,), _ = pr.score_all([_single("G")], None, _panel_settings(), {})            # Exomiser off: 'none' -> floor
    assert row.pheno == "strong"


def test_notes_name_the_panel():
    pe = PhenoEvidence({"G": 500}, {}, {})
    (row,), _ = pr.score_all([_single("G")], pe, _panel_settings(), {})
    assert "phenotype=strong" in pr.notes_for(row, pe) and "panel=P (phenotype floor)" in pr.notes_for(row, pe)
    (plain,), _ = pr.score_all([_single("G")], pe, S, {})
    assert "panel=" not in pr.notes_for(plain, pe)


# ----------------------------------------------------------------------------- alternative pairs, shared budget

def _alt_settings(k=3, min_pheno="strong", floor="strong"):
    return {**_panel_settings(floor=floor), "alternative_pairs": {"max_per_gene": k, "min_phenotype": min_pheno}}


def _three_pairs(settings, pe):
    """Scored primary rows of an AR panel gene with a stop (1), a moderate missense (5) and a 5'UTR panel candidate
    (20): pairs (1,5) full T2, (1,20) and (5,20) panel_partner T4; singles 1 (T4) and 5 (T5)."""
    cands = [make_cand(1, consequence="stop_gained"), make_cand(5, am=0.95),
             make_cand(20, consequence="5_prime_UTR_variant")]
    al = build_alleles(cands, {"G": AR}, settings, "MALE")["G"]
    primary, _ = pr.score_all(findings_for_gene("G", al, AR, settings), pe, settings, {})
    return primary


def _pos(rows):
    return [[a.first.call.pos for a in s.finding.alleles] for s in rows]


def test_score_all_keeps_the_tier_value_and_the_phase_factor():
    pe = PhenoEvidence({"G": 500}, {}, {})
    rows = _three_pairs(_alt_settings(), pe)
    by_pos = {tuple(p): s for p, s in zip(_pos(rows), rows)}
    pair, single = by_pos[(1, 5)], by_pos[(1,)]
    assert (pair.tier, pair.tier_epcr, pair.phase_factor, pair.epcr) == ("T2", 0.60, 0.5, 0.30)
    assert (single.tier, single.tier_epcr, single.phase_factor, single.epcr) == ("T4", 0.10, 1.0, 0.10)


def test_secondary_rows_keep_their_epcr_as_the_tier_value():
    tp53 = {"TP53": SfGene("TP53", "AD", "all_plp")}
    _, (sec,) = pr.score_all([_sf_finding("TP53", AD, clinvar="plp", stars=2)], PhenoEvidence({"TP53": 500}, {}, {}),
                             S, tp53)
    assert (sec.tier_epcr, sec.phase_factor) == (0.30, 1.0)


def test_per_gene_keeps_further_pairs_of_a_strong_gene():
    pe = PhenoEvidence({"G": 500}, {}, {})
    rows = _three_pairs(_alt_settings(), pe)                     # floor: strong
    assert _pos(pr.per_gene(rows, pe, _alt_settings(k=3))) == [[1, 5], [1, 20], [5, 20]]
    assert _pos(pr.per_gene(rows, pe, _alt_settings(k=2))) == [[1, 5], [1, 20]]
    assert _pos(pr.per_gene(rows, pe, _alt_settings(k=1))) == [[1, 5]]
    assert _pos(pr.per_gene(rows, pe, S)) == [[1, 5]]                     # the default is one per gene


def test_a_weak_phenotype_gene_gets_no_alternatives():
    s = _alt_settings(floor=None)
    pe = PhenoEvidence({"G": 500}, {}, {})
    rows = _three_pairs(s, pe)
    assert {r.pheno for r in rows} == {"weak"}
    assert _pos(pr.per_gene(rows, pe, s)) == [[1]]                         # the T5 single (0.02) beats 0.02 x 0.5
    weak = _alt_settings(k=3, min_pheno="weak", floor=None)
    assert _pos(pr.per_gene(rows, pe, weak)) == [[1], [1, 5], [1, 20]]      # pairs added after a single best row


def test_per_gene_skips_an_allele_set_already_kept():
    pe = PhenoEvidence({"G": 500}, {}, {})
    rows = _three_pairs(_alt_settings(), pe)
    pairs_in_h = [replace(s, finding=replace(s.finding, gene="H")) for s in rows if len(s.finding.alleles) == 2]
    kept = pr.per_gene(rows + pairs_in_h, pe, _alt_settings(k=3))                # an overlapping gene
    assert _pos(kept) == [[1, 5], [1, 20], [5, 20]] and {s.finding.gene for s in kept} == {"G"}
    assert pr.per_gene(rows + rows, pe, _alt_settings(k=3)) == kept               # duplicate pairs in one gene


def test_share_budget_splits_the_largest_tier_value_over_the_pairs_of_a_gene():
    pe = PhenoEvidence({"G": 500}, {}, {})
    kept = pr.per_gene(_three_pairs(_alt_settings(), pe), pe, _alt_settings(k=3))
    shared = pr.share_budget(kept, pe)
    assert [s.finding for s in shared] == [s.finding for s in kept]             # the order is kept
    tier_shares = [s.epcr / s.phase_factor for s in shared]
    assert sum(tier_shares) == pytest.approx(0.60)                              # B = the largest tier value (T2)
    assert tier_shares == pytest.approx([0.45, 0.075, 0.075])
    assert [s.epcr for s in shared] == pytest.approx([0.225, 0.0375, 0.0375])   # x p_trans 0.5 (unphased)
    assert [s.tier_epcr for s in shared] == [0.60, 0.10, 0.10]


def test_share_budget_example_t2_and_one_panel_partner():
    pe = PhenoEvidence({"G": 500}, {}, {})
    kept = pr.per_gene(_three_pairs(_alt_settings(), pe), pe, _alt_settings(k=2))
    assert [s.epcr for s in pr.share_budget(kept, pe)] == pytest.approx([0.257, 0.043], abs=5e-4)


def test_share_budget_applies_each_phase_factor_and_leaves_singles_and_lone_pairs():
    pe = PhenoEvidence({"G": 500}, {}, {})
    rows = _three_pairs(_alt_settings(), pe)
    kept = pr.per_gene(rows, pe, _alt_settings(k=2))
    trans = [kept[0], replace(kept[1], phase_factor=1.0)]                        # a proven-trans alternative
    assert [s.epcr for s in pr.share_budget(trans, pe)] == pytest.approx([0.6 * 6 / 7 * 0.5, 0.6 / 7])
    single = next(s for s in rows if _pos([s]) == [[1]])
    assert pr.share_budget([kept[0], single], pe) == [kept[0], single]          # one pair, one single: unchanged
    assert pr.share_budget(pr.best_per_gene(rows, pe), pe) == pr.best_per_gene(rows, pe)


# ----------------------------------------------------------------------------- compound_het.tsv

CH_HEADER = ("gene", "variant_1", "variant_2", "hgvs_1", "hgvs_2", "level_1", "level_2", "reasons_1", "reasons_2",
             "phase", "phase_sources", "distance", "collapsed_into", "fit", "tier", "epcr_before_cap", "reportable",
             "csv_row")


def test_compound_het_header_is_exact():
    assert pr.CH_COLUMNS == CH_HEADER


def _ch_cands():
    """Stop (1), moderate missense (5) and splice-region candidate (30) phased cis by PS, a 5'UTR candidate (20)."""
    return [make_cand(1, consequence="stop_gained", hgvs="p.Leu1*"),
            make_cand(5, call=make_call(5, gt="0|1", ps="5"), am=0.95, hgvs="p.Asn5Lys"),
            make_cand(20, consequence="5_prime_UTR_variant"),
            make_cand(30, call=make_call(30, gt="0|1", ps="5"), consequence="splice_region_variant")]


def test_compound_het_rows_list_every_het_pair_with_its_status():
    s = _alt_settings(k=3)
    pe = PhenoEvidence({"G": 500}, {}, {})
    r = pr.rank(_ch_cands(), {"G": AR}, {}, pe, s, "MALE")
    assert [[a.first.call.pos for a in c.finding.alleles] for c in r.chosen] == [[1, 5], [1, 20], [5, 20]]
    rows = {(x[1].split(":")[1], x[2].split(":")[1]): dict(zip(CH_HEADER, x))
            for x in pr.compound_het_rows(r, {"G": AR}, s)}
    assert list(rows) == [("1", "5"), ("1", "20"), ("1", "30"), ("5", "20"), ("5", "30"), ("20", "30")]
    best = rows[("1", "5")]
    assert (best["hgvs_1"], best["hgvs_2"], best["level_1"], best["level_2"]) == ("p.Leu1*", "p.Asn5Lys",
                                                                                "strong", "moderate")
    assert (best["fit"], best["tier"], best["epcr_before_cap"], best["reportable"], best["csv_row"]) == \
        ("full", "T2", "0.225", "yes", 1)
    assert (best["phase"], best["phase_sources"], best["distance"], best["reasons_1"]) == ("unknown", "", 4,
                                                                                          "null:strong")
    assert [rows[k]["csv_row"] for k in (("1", "20"), ("5", "20"))] == [2, 3]
    assert rows[("1", "20")]["fit"] == "panel_partner" and rows[("1", "20")]["reasons_2"] == "panel candidate (5utr)"
    cis = rows[("5", "30")]
    assert (cis["phase"], cis["phase_sources"], cis["reportable"]) == ("cis", "ps", "no (cis)")
    assert (cis["fit"], cis["tier"], cis["epcr_before_cap"], cis["csv_row"]) == ("", "", "", "")
    folded = rows[("1", "30")]
    assert folded["reportable"] == "no (collapsed)" and folded["collapsed_into"] == "1:1:A>G|1:5:A>G"
    assert rows[("20", "30")]["reportable"] == "no (two second-hit-only alleles)"
    assert all(rows[k]["collapsed_into"] == "" for k in rows if k != ("1", "30"))


def test_compound_het_rows_flag_one_variant_twice_and_skip_dominant_genes():
    twice = [make_cand(40, call=make_call(40, orig="1|40|A|G|1"), am=0.95),
             make_cand(40, call=make_call(40, orig="1|40|A|G,C|1"), am=0.95)]
    r = pr.rank(twice, {"G": AR}, {}, None, S, "MALE")
    (row,) = pr.compound_het_rows(r, {"G": AR}, S)
    assert row[1] == row[2] == "1:40:A>G" and row[-2] == "no (same variant)" and row[-1] == ""
    r = pr.rank(_ch_cands(), {"G": AD}, {}, None, S, "MALE")
    assert pr.compound_het_rows(r, {"G": AD}, S) == []
    r = pr.rank(_ch_cands(), {}, {}, None, _alt_settings(), "MALE")              # unknown MOI: pairs are listed
    assert len(pr.compound_het_rows(r, {}, _alt_settings())) == 6


def test_rank_k1_reproduces_one_row_per_gene():
    pe = PhenoEvidence({"G": 500}, {}, {})
    r = pr.rank(_ch_cands(), {"G": AR}, {}, pe, _alt_settings(k=1), "MALE")
    assert _pos(r.kept) == _pos(r.chosen) == [[1, 5]] and r.chosen[0].epcr == pytest.approx(0.30)


# ----------------------------------------------------------------------------- pair choice: the all_pairs order

def test_one_pair_per_gene_takes_the_stronger_second_allele_as_model_2_did():
    """Two moderate missense (1, 2) and a null (3, popmax 1e-3): (1, 2) and (1, 3) tie on epcr (T2 0.60 x 0.5) and
    on the weakest level. Model 2's best_pair took (1, 3): the stronger allele outranks rarity."""
    cands = [make_cand(1, am=0.95), make_cand(2, am=0.95), make_cand(3, consequence="stop_gained", popmax=0.001)]
    r = pr.rank(cands, {"G": AR}, {}, PhenoEvidence({"G": 3}, {}, {}), S, "MALE")
    assert _pos(r.chosen) == [[1, 3]]
    assert (r.chosen[0].tier, r.chosen[0].finding.level) == ("T2", "moderate")
    assert r.chosen[0].epcr == pytest.approx(0.30)


def test_alternative_pairs_follow_the_panel_class_order_not_rarity():
    """A null (100), a calibrated missense (10000) and three panel candidates: splice region (5000, popmax 2e-4),
    AlphaMissense 0.30 missense (7000, 1e-4) and a 3'UTR (20000, popmax 0). The panel pairs tie on epcr and level;
    the all_pairs order (better panel class first) picks (100, 5000) and (100, 7000), not the rarest 3'UTR."""
    cands = [make_cand(100, consequence="stop_gained"),
             make_cand(5000, consequence="splice_region_variant&intron_variant", popmax=0.0002),
             make_cand(7000, am=0.30, popmax=0.0001),
             make_cand(10000, am=0.95),
             make_cand(20000, consequence="3_prime_UTR_variant")]
    r = pr.rank(cands, {"G": AR}, {}, PhenoEvidence({"G": 500}, {}, {}), _alt_settings(k=3), "MALE")
    assert _pos(r.kept) == [[100, 10000], [100, 5000], [100, 7000]]
    assert [s.finding.fit for s in r.kept] == ["full", "panel_partner", "panel_partner"]


def test_a_proven_trans_pair_outranks_a_stronger_unphased_pair():
    """User ruling 2026-10-20, the one deviation from Model 2: null + null with unknown phase (T1 0.90 x 0.5 = 0.45)
    comes first in all_pairs, but null + moderate missense proven trans (T2 0.60) has the higher epcr and is the
    gene's row."""
    cands = [make_cand(1, consequence="stop_gained", call=make_call(1, gt="0|1", ps="1")),
             make_cand(5, call=make_call(5, gt="1|0", ps="1"), am=0.95),
             make_cand(9, consequence="stop_gained")]
    al = build_alleles(cands, {"G": AR}, S, "MALE")["G"]
    pairs = [f for f in findings_for_gene("G", al, AR, S) if len(f.alleles) == 2]
    ranked = {tuple(a.first.call.pos for a in f.alleles): (f.pair_rank, f.phase) for f in pairs}
    assert ranked[(1, 9)] == (0, "unknown") and ranked[(1, 5)] == (1, "trans")
    r = pr.rank(cands, {"G": AR}, {}, PhenoEvidence({"G": 3}, {}, {}), S, "MALE")
    assert _pos(r.chosen) == [[1, 5]]
    assert (r.chosen[0].tier, r.chosen[0].finding.phase, r.chosen[0].epcr) == ("T2", "trans", pytest.approx(0.60))
