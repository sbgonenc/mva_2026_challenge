import pytest
from pathlib import Path
from prioritise_helpers import write_vcf

from vep_pipeline.candidates import collect, excluded_region, header_fields, query_format
from vep_pipeline.config import ConfigError
from vep_pipeline.rank_settings import PRIORITISE_DEFAULTS

CSQ = ["Allele", "Consequence", "SYMBOL", "Feature", "MANE_SELECT", "REVEL", "gnomADg_NFE_AF", "gnomADg_FIN_AF"]
SETTINGS = {**PRIORITISE_DEFAULTS, "exclude_evidence_classes": []}


def _mis(symbol, revel, **af):
    return {"Consequence": "missense_variant", "SYMBOL": symbol, "Feature": f"ENST_{symbol}", "MANE_SELECT": "NM_1",
            "REVEL": str(revel), **af}


def test_collect_keeps_rare_variants_with_evidence(tmp_path):
    vcf = write_vcf(tmp_path / "v.vcf.gz", [
        {"pos": 100, "csq": [_mis("A1", 0.95), {"SYMBOL": "A2", "Consequence": "intron_variant", "Feature": "E2"}]},
        {"pos": 200, "csq": [_mis("B1", 0.95, gnomADg_NFE_AF="0.2")]},           # common: dropped
        {"pos": 300, "csq": [_mis("C1", 0.95, gnomADg_FIN_AF="0.3")]},           # FIN only: kept
        {"pos": 400, "csq": [{"SYMBOL": "D1", "Consequence": "intron_variant", "Feature": "E4"}]},  # no evidence
        {"pos": 500, "filter": "LowQual", "csq": [_mis("E1", 0.95)]},            # FILTER: dropped by view -f
        {"pos": 600, "GT": "0/0", "csq": [_mis("F1", 0.95)]},                    # no alt allele
        {"pos": 700, "orig": "1|700|AC|GT|1", "csq": [_mis("G1", 0.95)]},
    ], CSQ)
    cands = collect(vcf, "S1", SETTINGS, {})
    assert [(c.call.pos, c.ev.symbol) for c in cands] == [(100, "A1"), (300, "C1"), (700, "G1")]
    assert cands[2].call.orig == "1|700|AC|GT|1" and cands[0].call.ps == "."


def test_collect_keeps_a_variant_whose_only_evidence_is_remm(tmp_path):
    intron = {"SYMBOL": "R1", "Consequence": "intron_variant", "Feature": "E1"}
    vcf = write_vcf(tmp_path / "r.vcf.gz", [
        {"pos": 100, "csq": [intron]},
        {"pos": 200, "csq": [intron]},
    ], CSQ)
    assert collect(vcf, "S1", SETTINGS, {}) == []
    cands = collect(vcf, "S1", SETTINGS, {}, remm={"1:100:A>G": 0.97, "1:200:A>G": 0.5})
    assert [(c.call.pos, c.ev.remm) for c in cands] == [(100, 0.97)]


def test_query_format_replaces_absent_tags():
    fmt = query_format({"CSQ"}, {"GT", "AD"})
    assert "%INFO/EXO_ORIG" not in fmt and "%PS" not in fmt and "%GT" in fmt and fmt.endswith("%INFO/CSQ\n")


def test_header_without_csq_is_an_error(tmp_path):
    plain = tmp_path / "n.vcf"
    plain.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    with pytest.raises(ConfigError, match="no CSQ header"):
        header_fields(plain)


def test_real_csq_header_has_the_fields_we_read():
    line = (Path(__file__).parent / "data" / "prioritise" / "csq_header.txt").read_text()
    fields = set(line.split("Format:")[1].split('"')[0].strip().split("|"))
    assert {"SYMBOL", "Consequence", "Feature", "MANE_SELECT", "CANONICAL", "HGVSp", "REVEL", "am_pathogenicity",
            "ClinVar_CLNSIG", "ClinVar_CLNREVSTAT", "NMD", "LOEUF", "gnomADg_NFE_AF"} <= fields


def test_collect_keeps_variant_scored_only_when_revel_is_excluded(tmp_path):
    fields = CSQ + ["am_pathogenicity"]
    vcf = write_vcf(tmp_path / "am.vcf.gz", [
        {"pos": 100, "csq": [_mis("H1", 0.5, am_pathogenicity="0.95")]},
    ], fields)
    cands = collect(vcf, "S1", SETTINGS, {})
    assert [(c.call.pos, c.ev.symbol) for c in cands] == [(100, "H1")]


from prioritise_helpers import make_call, make_cand

from vep_pipeline.candidates import (all_pairs, best_pair, build_alleles, findings_for_gene, orig_record, phase,
                                     phase_graph, qc_ok, vaf, zygosity)
from vep_pipeline.genedb import GeneInfo

AR = GeneInfo(frozenset({"AR"}), frozenset({"loss of function"}))
AD = GeneInfo(frozenset({"AD"}), frozenset({"loss of function"}))
QC = PRIORITISE_DEFAULTS["qc"]


def test_zygosity():
    assert zygosity(make_call(1, "0/1"), "MALE") == "het"
    assert zygosity(make_call(1, "1/0"), "MALE") == "het"                          # split multiallelic sibling
    assert zygosity(make_call(1, "1/1"), "MALE") == "hom"
    assert zygosity(make_call(1, "1"), "UNKNOWN_SEX") == "hemi"
    assert zygosity(make_call(5_000_000, "1/1", chrom="X"), "MALE") == "hemi"
    assert zygosity(make_call(100_000, "1/1", chrom="X"), "MALE") == "hom"          # PAR1
    assert zygosity(make_call(1, "./."), "MALE") is None


def test_vaf_and_qc():
    assert vaf(make_call(1, ad="10,2")) == pytest.approx(2 / 12)
    assert qc_ok(make_call(1), QC)
    assert not qc_ok(make_call(1, filt="LowQual"), QC)
    assert not qc_ok(make_call(1, gq="10"), QC)
    assert qc_ok(make_call(1, gq=".", dp="."), QC)


def test_orig_record():
    assert orig_record(make_call(1)) is None
    assert orig_record(make_call(1, orig="1|8000|A|C,G|2")) == ("1", 8000, "A", "C,G", 2)
    with pytest.raises(ConfigError, match="EXO_ORIG"):
        orig_record(make_call(1, orig="1:8000"))


def test_atomised_mnv_pieces_become_one_allele():
    o = "1|5000|AC|GT|1"
    cands = [make_cand(5000, call=make_call(5000, orig=o), revel=0.95),
             make_cand(5001, call=make_call(5001, orig=o), revel=0.95)]
    al = build_alleles(cands, {"G": AR}, SETTINGS, "MALE")["G"]
    assert len(al) == 1 and len(al[0].pieces) == 2 and al[0].level == "strong"
    assert findings_for_gene("G", al, AR, SETTINGS)[0].fit == "single_hit"          # never a pair with itself


def test_low_vaf_het_is_kept_and_lowered():
    al = build_alleles([make_cand(1, call=make_call(1, ad="18,2"), revel=0.95)], {}, SETTINGS, "MALE")["G"]
    assert al[0].low_vaf and al[0].level == "moderate"


def _alleles(*cands):
    return build_alleles(list(cands), {"G": AR}, SETTINGS, "MALE")["G"]


def test_phase():
    a, b = _alleles(make_cand(8000, call=make_call(8000, alt="C", orig="1|8000|A|C,G|1"), revel=0.95),
                    make_cand(8000, call=make_call(8000, alt="G", orig="1|8000|A|C,G|2"), revel=0.95))
    assert phase(a, b) == "trans"                                                    # two ALTs of one record
    a, b = _alleles(make_cand(1, call=make_call(1, gt="0|1", ps="7"), revel=0.95),
                    make_cand(9, call=make_call(9, gt="1|0", ps="7"), revel=0.95))
    assert phase(a, b) == "trans"
    a, b = _alleles(make_cand(1, call=make_call(1, gt="0|1", ps="7"), revel=0.95),
                    make_cand(9, call=make_call(9, gt="0|1", ps="7"), revel=0.95))
    assert phase(a, b) == "cis"
    a, b = _alleles(make_cand(1, call=make_call(1, gt="0/1", pgt="0|1", pid="p"), revel=0.95),
                    make_cand(9, call=make_call(9, gt="0/1", pgt="1|0", pid="p"), revel=0.95))
    assert phase(a, b) == "trans"
    a, b = _alleles(make_cand(1, revel=0.95), make_cand(9, revel=0.95))
    assert phase(a, b) == "unknown"


def test_best_pair_skips_cis_and_needs_one_real_hit():
    a, b, c = _alleles(make_cand(1, call=make_call(1, gt="0|1", ps="7"), revel=0.95),
                       make_cand(9, call=make_call(9, gt="0|1", ps="7"), revel=0.95),
                       make_cand(20, revel=0.70))
    (x, y), ph = best_pair([a, b, c])
    assert {x.first.call.pos, y.first.call.pos} == {1, 20} and ph == "unknown"
    s1, s2 = _alleles(make_cand(1, consequence="5_prime_UTR_variant", noncoding_signal=True),
                      make_cand(9, consequence="5_prime_UTR_variant", noncoding_signal=True))
    assert best_pair([s1, s2]) is None


def test_findings_by_inheritance():
    two = _alleles(make_cand(1, revel=0.95), make_cand(9, consequence="stop_gained"))
    fs = findings_for_gene("G", two, AR, SETTINGS)
    assert [(f.fit, len(f.alleles)) for f in fs] == [("single_hit", 1), ("single_hit", 1), ("full", 2)]
    assert fs[-1].level == "strong" and fs[-1].phase == "unknown" and fs[-1].distance == 8
    hom = build_alleles([make_cand(1, call=make_call(1, gt="1/1"), revel=0.95)], {}, SETTINGS, "MALE")["G"]
    assert findings_for_gene("G", hom, AR, SETTINGS)[0].fit == "full"
    het = build_alleles([make_cand(1, revel=0.95)], {}, SETTINGS, "MALE")["G"]
    assert findings_for_gene("G", het, AD, SETTINGS)[0].fit == "full"
    assert findings_for_gene("G", het, None, SETTINGS)[0].fit == "unknown_moi"
    rare_not_ultra = build_alleles([make_cand(1, revel=0.95, popmax=0.001)], {}, SETTINGS, "MALE")["G"]
    assert findings_for_gene("G", rare_not_ultra, AD, SETTINGS) == []               # 1e-3 > dominant cut-off 1e-4


def test_best_pair_never_pairs_an_allele_with_its_own_normalised_key():
    # two original records that normalise to the same variant 1:100:CA>C
    a, b = _alleles(make_cand(100, call=make_call(100, ref="CA", alt="C", orig="1|100|CAA|CA|1"), revel=0.95),
                    make_cand(100, call=make_call(100, ref="CA", alt="C", orig="1|100|CA|C|1"), revel=0.95))
    assert a.key == b.key == "1:100:CA>C"
    assert best_pair([a, b]) is None


def _split(pos, alt, idx, **kw):
    return make_cand(pos, call=make_call(pos, alt=alt, orig=f"1|8000|A|C,G|{idx}", **kw), revel=0.95)


def test_phase_graph_is_transitive_across_sources():
    a, b, c = _alleles(_split(8000, "C", 1), _split(8000, "G", 2, gt="0|1", ps="7"),
                       make_cand(9000, call=make_call(9000, gt="0|1", ps="7"), revel=0.95))
    g = phase_graph([a, b, c])
    assert g.relation(a, b) == ("trans", ("multiallelic",))
    assert g.relation(b, c) == ("cis", ("ps",))
    assert g.relation(a, c) == ("trans", ("multiallelic", "ps"))
    assert g.sides(a) == ((), (b, c))
    assert g.sides(b) == ((c,), (a,))


def test_relation_lists_every_source_of_a_direct_link_whatever_else_is_in_the_graph():
    a, b, c = _alleles(make_cand(1, call=make_call(1, gt="0|1", ps="7"), revel=0.95),
                       make_cand(9, call=make_call(9, gt="0|1", ps="7", pgt="0|1", pid="p"), revel=0.95),
                       make_cand(20, call=make_call(20, gt="0|1", ps="7", pgt="0|1", pid="p"), revel=0.95))
    assert phase_graph([b, c]).relation(b, c) == ("cis", ("pid", "ps"))
    assert phase_graph([a, b, c]).relation(b, c) == ("cis", ("pid", "ps"))
    assert phase_graph([a, b, c]).relation(a, c) == ("cis", ("ps",))
    assert phase_graph([a, b], sources=("bogus",)).relation(a, b) == ("unknown", ())      # unknown names add no edge


def test_phase_graph_conflict_marks_the_whole_component():
    a, b, c = _alleles(make_cand(1, call=make_call(1, gt="0|1", ps="7", pgt="0|1", pid="p"), revel=0.95),
                       make_cand(9, call=make_call(9, gt="0|1", ps="7", pgt="1|0", pid="p"), revel=0.95),
                       make_cand(20, call=make_call(20, gt="0|1", ps="7"), revel=0.95))
    g = phase_graph([a, b, c])
    assert g.relation(a, b) == ("unknown", ("conflict",))
    assert g.relation(a, c) == ("unknown", ("conflict",))                            # whole component
    assert g.sides(c) == ((), ())
    assert phase_graph([a, b], sources=("ps",)).relation(a, b) == ("cis", ("ps",))   # a disabled source adds no edge
    assert all_pairs([a, b, c], g)[0].sources == ("conflict",)


def test_phase_graph_unlinked_and_readback():
    a, b, c = _alleles(make_cand(1, revel=0.95), make_cand(9, revel=0.95), make_cand(20, revel=0.95))
    assert phase_graph([a, b]).relation(a, b) == ("unknown", ())
    rb = {a.key: ("0|1", "3"), b.key: ("1|0", "3"), c.key: ("0|1", "4")}
    g = phase_graph([a, b, c], readback=rb)
    assert g.relation(a, b) == ("trans", ("whatshap",))
    assert g.relation(a, c) == ("unknown", ())                                       # different PS
    assert phase_graph([a, b], readback=rb, sources=("ps",)).relation(a, b) == ("unknown", ())


def test_all_pairs_collapse_a_cis_haplotype_onto_its_strongest_member():
    a, b, c = _alleles(make_cand(1, consequence="stop_gained", call=make_call(1, gt="0|1", ps="7")),
                       make_cand(9, call=make_call(9, gt="0|1", ps="7"), revel=0.80),
                       make_cand(20, revel=0.70))
    assert (a.level, b.level, c.level) == ("strong", "moderate", "supporting")
    pairs = all_pairs([a, b, c], phase_graph([a, b, c]))
    assert [(p.a, p.b, p.phase, p.companions) for p in pairs] == [(a, c, "unknown", (b,))]


def test_all_pairs_lists_every_pair_best_first_and_proven_trans_wins_ties():
    a, b, c = _alleles(_split(8000, "C", 1), _split(8000, "G", 2),
                       make_cand(9000, revel=0.95))
    pairs = all_pairs([a, b, c], phase_graph([a, b, c]))
    assert [(p.a, p.b, p.phase) for p in pairs] == [(a, b, "trans"), (a, c, "unknown"), (b, c, "unknown")]
    assert pairs[0].sources == ("multiallelic",) and pairs[0].companions == ()


def test_a_gene_whose_only_variants_are_in_cis_has_singles_only():
    a, b = _alleles(make_cand(1, call=make_call(1, gt="0|1", ps="7"), revel=0.95),
                    make_cand(9, call=make_call(9, gt="0|1", ps="7"), revel=0.95))
    assert all_pairs([a, b], phase_graph([a, b])) == [] and best_pair([a, b]) is None
    fs = findings_for_gene("G", [a, b], AR, SETTINGS)
    assert [(f.fit, len(f.alleles)) for f in fs] == [("single_hit", 1), ("single_hit", 1)]


def test_findings_carry_phase_source_and_companions():
    a, b, c = _alleles(make_cand(1, consequence="stop_gained", call=make_call(1, gt="0|1", ps="7")),
                       make_cand(9, call=make_call(9, gt="0|1", ps="7"), revel=0.80),
                       make_cand(20, call=make_call(20, gt="1|0", ps="7"), revel=0.70))
    pair = [f for f in findings_for_gene("G", [a, b, c], AR, SETTINGS) if len(f.alleles) == 2]
    assert len(pair) == 1 and pair[0].phase == "trans" and pair[0].phase_source == "ps"
    assert pair[0].companions == (b,)
    off = [f for f in findings_for_gene("G", [a, b, c], AR, {**SETTINGS, "phase_sources": ()}) if len(f.alleles) == 2]
    assert len(off) == 3 and {f.phase for f in off} == {"unknown"} and {f.phase_source for f in off} == {""}


def test_pair_findings_carry_their_all_pairs_rank():
    a, b, c = _alleles(make_cand(1, am=0.95), make_cand(2, am=0.95),
                       make_cand(3, consequence="stop_gained", popmax=0.001))
    order = [(h.a.first.call.pos, h.b.first.call.pos) for h in all_pairs([a, b, c], phase_graph([a, b, c]))]
    assert order == [(1, 3), (2, 3), (1, 2)]                       # stronger second allele before rarity
    fs = findings_for_gene("G", [a, b, c], AR, SETTINGS)
    pairs = [f for f in fs if len(f.alleles) == 2]
    assert [(tuple(x.first.call.pos for x in f.alleles), f.pair_rank) for f in pairs] == [
        ((1, 3), 0), ((2, 3), 1), ((1, 2), 2)]
    assert {f.pair_rank for f in fs if len(f.alleles) == 1} == {0}


def test_build_alleles_uses_the_configured_missense_tool():
    c = make_cand(1, revel=0.472, am=0.923)
    assert build_alleles([c], {}, SETTINGS, "MALE")["G"][0].level == "moderate"            # default: AlphaMissense
    assert build_alleles([c], {}, {**SETTINGS, "missense_tool": "revel"}, "MALE") == {}    # REVEL 0.472: no evidence


def test_excluded_region_matches_with_or_without_chr():
    regions = [{"name": "MHC", "chrom": "6", "start": 100, "end": 200}]
    assert excluded_region(make_call(150, chrom="6"), regions) == "MHC"
    assert excluded_region(make_call(150, chrom="chr6"), regions) == "MHC"
    assert excluded_region(make_call(201, chrom="6"), regions) is None
    assert excluded_region(make_call(150, chrom="7"), regions) is None


def test_collect_skips_and_reports_excluded_regions(tmp_path):
    vcf = write_vcf(tmp_path / "v.vcf.gz", [{"pos": 100, "csq": [_mis("A1", 0.95)]},
                                            {"pos": 300, "csq": [_mis("C1", 0.95)]}], CSQ)
    settings = {**SETTINGS, "exclude_regions": [{"name": "TEST", "chrom": "1", "start": 50, "end": 150}]}
    excluded = []
    cands = collect(vcf, "S1", settings, {}, excluded=excluded)
    assert [c.ev.symbol for c in cands] == ["C1"] and excluded == [("1:100:A>G", "TEST")]


def test_whatshap_readback_phases_a_pair_and_is_switchable():
    a, b = _alleles(make_cand(1, consequence="stop_gained"), make_cand(9, revel=0.95))
    rb = {a.key: ("0|1", "3"), b.key: ("1|0", "3")}
    pair = lambda s, **kw: [f for f in findings_for_gene("G", [a, b], AR, s, **kw) if len(f.alleles) == 2]
    on = pair(SETTINGS, readback=rb)
    assert on[0].phase == "trans" and on[0].phase_source == "whatshap"
    off = pair({**SETTINGS, "phase_sources": ["multiallelic", "ps", "pid"]}, readback=rb)
    assert off[0].phase == "unknown" and off[0].phase_source == ""


# ----------------------------------------------------------------------------- gene panel (deep second-hit search)

PANEL = {**SETTINGS, "gene_panel": {"name": "P", "genes": ["G"], "phenotype_floor": None, "deep_second_hit": True}}


def _panel_alleles(*cands, settings=PANEL):
    return build_alleles(list(cands), {"G": AR}, settings, "MALE")["G"]


def test_collect_keeps_panel_candidates_only_in_panel_genes_with_deep_search(tmp_path):
    vcf = write_vcf(tmp_path / "p.vcf.gz", [
        {"pos": 100, "csq": [_mis("G", 0.2)]},                   # no calibrated evidence
        {"pos": 200, "csq": [_mis("H", 0.2)]},
        {"pos": 300, "csq": [{"Consequence": "synonymous_variant", "SYMBOL": "G", "Feature": "E3"}]},
    ], CSQ)
    assert collect(vcf, "S1", SETTINGS, {}) == []
    assert [c.call.pos for c in collect(vcf, "S1", PANEL, {})] == [100]
    shallow = {**PANEL, "gene_panel": {**PANEL["gene_panel"], "deep_second_hit": False}}
    assert collect(vcf, "S1", shallow, {}) == []


def test_panel_only_allele_is_second_hit_only_and_never_a_single():
    null, cand = _panel_alleles(make_cand(1, consequence="stop_gained"), make_cand(9, am=0.5))
    assert not null.panel_only and cand.panel_only and cand.second_hit_only and cand.level == "supporting"
    assert cand.reasons == ("panel candidate (protein_altering)",)
    fs = findings_for_gene("G", [null, cand], AR, PANEL)
    assert [(f.fit, len(f.alleles)) for f in fs] == [("single_hit", 1), ("panel_partner", 2)]
    assert build_alleles([make_cand(9, am=0.5)], {"G": AR}, SETTINGS, "MALE") == {}       # panel off: dropped
    other = build_alleles([make_cand(9, symbol="H", am=0.5)], {"H": AR}, PANEL, "MALE")
    assert other == {}                                                                      # not a panel gene


def test_panel_flag_does_not_make_a_calibrated_allele_panel_only():
    (a,) = _panel_alleles(make_cand(1, am=0.95))
    assert not a.panel_only and not a.second_hit_only


def test_two_panel_candidates_never_pair():
    a, b = _panel_alleles(make_cand(1, am=0.5), make_cand(9, am=0.5))
    assert all_pairs([a, b], phase_graph([a, b])) == []
    assert findings_for_gene("G", [a, b], AR, PANEL) == []


def test_a_pair_with_a_panel_partner_is_panel_partner_whatever_the_moi():
    null, cand = _panel_alleles(make_cand(1, consequence="stop_gained"), make_cand(9, am=0.5))
    assert [f.fit for f in findings_for_gene("G", [null, cand], None, PANEL) if len(f.alleles) == 2] == [
        "panel_partner"]
    real = _alleles(make_cand(1, revel=0.95), make_cand(9, consequence="stop_gained"))
    assert [f.fit for f in findings_for_gene("G", real, AR, PANEL) if len(f.alleles) == 2] == ["full"]


def test_all_pairs_prefers_the_lower_panel_class():
    null, utr, splice, mis = _panel_alleles(
        make_cand(1, consequence="stop_gained"),
        make_cand(5, consequence="3_prime_UTR_variant"),
        make_cand(9, consequence="intron_variant", splice=0.15),
        make_cand(20, am=0.5))
    pairs = all_pairs([null, utr, splice, mis], phase_graph([null, utr, splice, mis]))
    assert [p.b for p in pairs] == [splice, mis, utr]
    assert all(p.a is null for p in pairs)


def test_collect_keeps_a_bp4_splice_region_variant_for_a_setting_without_spliceai(tmp_path):
    """SpliceAI <= 0.1 (BP4) vetoes a splice-region term; with unsupervised_predictor excluded the term alone makes a
    panel candidate, so the pre-filter must keep the variant."""
    vcf = write_vcf(tmp_path / "s.vcf.gz", [
        {"pos": 100, "csq": [{"Consequence": "splice_region_variant&intron_variant", "SYMBOL": "G", "Feature": "E1",
                              "SpliceAI_pred_DS_AG": "0.05"}]}], CSQ + ["SpliceAI_pred_DS_AG"])
    assert collect(vcf, "S1", SETTINGS, {}) == []                                      # not a panel gene
    cands = collect(vcf, "S1", PANEL, {})
    assert [c.call.pos for c in cands] == [100]
    assert build_alleles(cands, {"G": AR}, PANEL, "MALE") == {}                        # BP4: no candidate
    no_spliceai = {**PANEL, "exclude_evidence_classes": ["unsupervised_predictor"]}
    (a,) = build_alleles(collect(vcf, "S1", no_spliceai, {}), {"G": AR}, no_spliceai, "MALE")["G"]
    assert a.reasons == ("panel candidate (splice)",)
