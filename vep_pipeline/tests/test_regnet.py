import pytest
from prioritise_helpers import make_call, write_vcf

from vep_pipeline.candidates import collect
from vep_pipeline.rank_settings import PRIORITISE_DEFAULTS
from vep_pipeline.regnet import HitCollector, PeakIndex, bh, hypergeom_sf, reg_tier

POPS = ("NFE", "FIN")
INDEX = {"chr1:100-600": ("chr1", 100, 600), "chr1:1000-1500": ("chr1", 1000, 1500)}


def test_peak_index_inclusive_bounds_and_chr_prefix():
    ix = PeakIndex(INDEX)
    assert ix.find("1", 100, 100) == ["chr1:100-600"]
    assert ix.find("chr1", 600, 600) == ["chr1:100-600"]
    assert ix.find("1", 601, 999) == []
    assert ix.find("1", 590, 1010) == ["chr1:100-600", "chr1:1000-1500"]   # deletion spanning both
    assert ix.find("7", 100, 100) == []                                    # unknown contig


def test_hypergeom_matches_known_value():
    # N=20, K=5 successes, n=4 draws, P(X>=2) = 1 - [C(15,4) + 5*C(15,3)]/C(20,4) = 1 - (1365+2275)/4845
    assert abs(hypergeom_sf(2, 20, 5, 4) - (1 - 3640 / 4845)) < 1e-12
    assert hypergeom_sf(6, 20, 5, 4) == 0.0


def test_bh_monotone_and_capped():
    assert bh([0.01, 0.04, 0.03]) == pytest.approx([0.03, 0.04, 0.04])


@pytest.mark.parametrize("pheno,v,s,tier", [("strong", 1, 1, "R1"), ("strong", 1, 0, "R2"), ("strong", 0, 1, "R2"),
    ("strong", 0, 0, "R3"), ("moderate", 1, 1, "R2"), ("moderate", 0, 1, "R3"), ("moderate", 0, 0, "R4"),
    ("weak", 1, 1, "R5"), ("unknown", 1, 1, "R5"), ("none", 1, 1, "R5")])
def test_reg_tier_table(pheno, v, s, tier):
    assert reg_tier(pheno, bool(v), bool(s)) == tier


def test_collector_skips_coding_common_and_unlinked():
    hc = HitCollector(PeakIndex(INDEX), POPS, 0.01)
    inter = {"Consequence": "intergenic_variant", "CADD_PHRED": "12.5"}
    hc(make_call(200), [{"Consequence": "missense_variant", "SYMBOL": "G1"}])               # coding
    hc(make_call(200), [{"Consequence": "intron_variant&synonymous_variant", "SYMBOL": "G1"}])  # coding term in '&'
    hc(make_call(210), [{**inter, "gnomADg_NFE_AF": "0.02"}])                               # common
    hc(make_call(800), [inter])                                                              # outside every peak
    assert hc.hits == []
    hc(make_call(220), [inter, {"Consequence": "intron_variant", "SYMBOL": "G2", "CADD_PHRED": "20.1",
                                "gnomADe_FIN_AF": "0.004"}])
    (hit,) = hc.hits
    assert hit.call.pos == 220 and hit.peaks == ("chr1:100-600",)
    assert hit.cadd == pytest.approx(20.1) and hit.popmax == pytest.approx(0.004)
    assert isinstance(hit.entries, tuple) and len(hit.entries) == 2


def test_collector_uses_ref_span_and_no_gate():
    hc = HitCollector(PeakIndex(INDEX), POPS, None)
    hc(make_call(595, ref="A" * 10, alt="A"), [{"Consequence": "intron_variant"}])   # deletion reaching peak end
    hc(make_call(598, ref="A", alt="G"), [{"Consequence": "intron_variant", "gnomADg_NFE_AF": "0.4"}])
    assert [h.call.pos for h in hc.hits] == [595, 598] and hc.hits[0].cadd is None
    assert hc.hits[1].popmax == pytest.approx(0.4)


def test_collector_skips_records_without_csq_and_spanning_deletions():
    hc = HitCollector(PeakIndex(INDEX), POPS, None)
    hc(make_call(200), [])                                                  # no CSQ entry
    hc(make_call(210, alt="*"), [{"Consequence": "intron_variant"}])        # spanning deletion allele
    assert hc.hits == []
    hc(make_call(220), [{"Consequence": "intron_variant"}])
    assert [h.call.pos for h in hc.hits] == [220]


def test_collect_calls_hook_after_excluded_region(tmp_path):
    csq = ["Allele", "Consequence", "SYMBOL", "Feature"]
    ent = {"Consequence": "intron_variant", "SYMBOL": "A1", "Feature": "E1"}
    vcf = write_vcf(tmp_path / "v.vcf.gz", [{"pos": 100, "csq": [ent]}, {"pos": 200, "csq": [ent]},
                                            {"pos": 300, "csq": [ent]}], csq)
    settings = {**PRIORITISE_DEFAULTS, "exclude_evidence_classes": [],
                "exclude_regions": [{"name": "MHC", "chrom": "1", "start": 150, "end": 250}]}
    seen = []
    collect(vcf, "S1", settings, {}, hook=lambda call, entries: seen.append((call.pos, len(entries))))
    assert seen == [(100, 1), (300, 1)]


# ----------------------------------------------------------------------------- Task 5: scoring and outputs

from test_regnet_build import make_scdori, make_sldsc

from vep_pipeline.config import ConfigError
from vep_pipeline.exomiser import PhenoEvidence, gene_pheno
from vep_pipeline.regnet import (RegHit, Resource, SldscRow, load_resource, regulon_flags, score, sldsc_flags,
                                 trait_null, write_regulatory, write_summary)
from vep_pipeline.regnet_build import build

PEAK_A, PEAK_B = "chr1:100-600", "chr1:1000-1500"
RN = {"traits": ["T1", "T2"], "fdr": 0.10, "second_hit": True, "null_draws": 30, "seed": 7}


def _pe(ranks):
    return PhenoEvidence(rank=dict(ranks), combined={}, exomiser_rank={})


def _settings(**kw):
    return {**PRIORITISE_DEFAULTS, "exclude_evidence_classes": [], "pheno_rank": {"strong": 2, "moderate": 4}, **kw}


def _resource(peak_tfs=None, sldsc=None, regulons=None, links=None, universe=None, traits=None):
    links = links or {PEAK_A: ("GENEA", "GENE_B"), PEAK_B: ("GENEA",)}
    return Resource(index=PeakIndex({PEAK_A: ("chr1", 100, 600), PEAK_B: ("chr1", 1000, 1500)}),
                    links=links,
                    peak_tfs=peak_tfs if peak_tfs is not None else {PEAK_A: ("TF1", "TF2"), PEAK_B: ("TF2",)},
                    regulons=regulons or {}, universe=frozenset(universe or ("GENEA", "GENE_B")),
                    sldsc=sldsc if sldsc is not None else [SldscRow("T1", "TF1", 1e-9, 0.004)],
                    traits=traits or {"T1": "trait one", "T2": "trait two", "T3": ""})


def _hit(pos, peaks=(PEAK_A,), cadd=None, call=None, **ckw):
    return RegHit(call or make_call(pos, **ckw), ({"Consequence": "intron_variant", "SYMBOL": "X"},), tuple(peaks),
                  0.0, cadd)


def test_gene_pheno_applies_panel_floor():
    pe = _pe({"G1": 1, "G2": 3, "G3": 9})
    s = _settings()
    assert gene_pheno(pe, "G1", s) == ("strong", "")
    assert gene_pheno(pe, "G3", s) == ("weak", "")
    assert gene_pheno(pe, "GX", s) == ("unknown", "")
    assert gene_pheno(None, "G1", s) == ("none", "")
    s = _settings(gene_panel={"name": "P", "genes": ["G2", "G3", "GX"], "phenotype_floor": "moderate",
                              "deep_second_hit": False})
    assert gene_pheno(pe, "G2", s) == ("moderate", "P")                     # already at the floor
    assert gene_pheno(pe, "G3", s) == ("moderate", "P (phenotype floor)")
    assert gene_pheno(pe, "GX", s) == ("moderate", "P (phenotype floor)")
    s["gene_panel"] = {"name": "", "genes": ["G3"], "phenotype_floor": None, "deep_second_hit": False}
    assert gene_pheno(pe, "G3", s) == ("weak", "gene_panel")


def test_sldsc_flags_requires_positive_tau_and_fdr():
    rows = [SldscRow("T1", "TF1", 1e-9, 0.004), SldscRow("T2", "TF1", 2e-9, 0.05),
            SldscRow("T1", "TF2", -1e-9, 0.001),          # tau <= 0: never
            SldscRow("T1", "TF3", 0.0, 0.001),
            SldscRow("T1", "TF4", 1e-9, 0.10),            # fdr not below the cut
            SldscRow("T3", "TF5", 1e-9, 0.001)]           # undeclared trait
    assert sldsc_flags(rows, ["T2", "T1"], 0.10) == {"TF1": [("T1", 0.004), ("T2", 0.05)]}
    assert sldsc_flags(rows, [], 0.10) == {}


def test_regulon_flags_hypergeom_bh():
    universe = frozenset(f"G{i}" for i in range(20))
    pheno = {"G0", "G1", "G2", "G3", "G4"}
    regulons = {"A": frozenset({"G0", "G1", "G2", "G3", "OUTSIDE"}), "B": frozenset({"G10", "G11", "G12", "G13"}),
                "C": frozenset({"OUTSIDE"})}                       # no gene in the universe: not tested
    got = regulon_flags(regulons, universe, pheno, 0.10)
    assert set(got) == {"A"}
    assert got["A"] == pytest.approx(min(1.0, 2 * 5 / 4845))      # BH over the two tested TFs
    assert regulon_flags(regulons, universe, set(), 0.10) == {}


def test_score_rows_one_per_linked_gene_and_order():
    pe = _pe({"GENEA": 1})
    res = _resource()
    hits = [_hit(200, cadd=25.0),                                   # PEAK_A: GENEA (strong) and GENE_B (unknown)
            _hit(300, cadd=None),                                   # no V
            _hit(400, cadd=30.0, gq="5"),                           # QC fails: dropped
            _hit(1200, peaks=(PEAK_B,), cadd=22.0),                 # PEAK_B: TF2 only, not S-LDSC flagged
            _hit(250, cadd=21.0)]                                   # same tier as 200, lower CADD but ReMM
    rows, summary = score(hits, res, pe, {"1:250:A>G": 0.5}, _settings(), RN)
    got = [(r.hit.call.pos, r.gene, r.tier) for r in rows]
    # R1: ReMM before CADD; R2: one S-LDSC TF before none; R5: ReMM, then CADD, then none
    assert got == [(250, "GENEA", "R1"), (200, "GENEA", "R1"), (300, "GENEA", "R2"), (1200, "GENEA", "R2"),
                   (250, "GENE_B", "R5"), (200, "GENE_B", "R5"), (300, "GENE_B", "R5")]
    r = rows[1]
    assert r.peak == PEAK_A and r.tfs == ("TF1", "TF2") and r.v and r.second_hit and r.pheno == "strong"
    assert [tf for tf, _ in r.sldsc_tfs] == ["TF1"] and r.sldsc_tfs[0][1] == (("T1", 0.004),)
    assert rows[0].remm == 0.5 and rows[1].remm is None
    assert summary["hits"] == 5 and summary["rows"] == 7 and summary["rows_R1"] == 2 and summary["rows_R5"] == 3
    assert summary["second_hit_rows"] == 3


def test_score_dedupes_variant_gene_over_peaks_keeping_best_tier():
    pe = _pe({"GENEA": 1})
    res = _resource()
    hits = [_hit(590, peaks=(PEAK_A, PEAK_B), cadd=None, ref="A" * 500, alt="A")]   # spans both peaks
    rows, _ = score(hits, res, pe, {}, _settings(), RN)
    assert [(r.gene, r.peak, r.tier) for r in rows] == [("GENEA", PEAK_A, "R2"), ("GENE_B", PEAK_A, "R5")]


def test_score_dedupe_breaks_ties_on_regulon_tfs_before_peak():
    universe = ("GENEA", "GENE_B", *(f"G{i}" for i in range(19)))
    res = _resource(peak_tfs={PEAK_A: ("TF1",), PEAK_B: ("TF2",)}, sldsc=[], links={PEAK_A: ("GENEA",), PEAK_B: ("GENEA",)},
                    regulons={"TF2": frozenset({"GENEA"})}, universe=universe)
    hits = [_hit(590, peaks=(PEAK_A, PEAK_B), cadd=None, ref="A" * 500, alt="A")]   # same tier on both peaks
    (r,) = score(hits, res, _pe({"GENEA": 1}), {}, _settings(), RN)[0]
    assert (r.tier, r.peak, [tf for tf, _q in r.regulon_tfs]) == ("R3", PEAK_B, ["TF2"])


def test_v_respects_excluded_classes():
    pe = _pe({"GENEA": 1})
    res = _resource(links={PEAK_A: ("GENEA",)})
    hit = _hit(200, cadd=30.0)
    remm = {"1:200:A>G": 0.99}
    (r,) = score([hit], res, pe, {}, _settings(exclude_evidence_classes=["unsupervised_predictor"]), RN)[0]
    assert not r.v and r.tier == "R2"
    (r,) = score([hit], res, pe, remm, _settings(exclude_evidence_classes=["unsupervised_predictor"]), RN)[0]
    assert r.v and r.tier == "R1"
    (r,) = score([hit], res, pe, remm,
                 _settings(exclude_evidence_classes=["unsupervised_predictor", "clinically_trained_predictor"]),
                 RN)[0]
    assert not r.v
    (r,) = score([_hit(200, cadd=19.9)], res, pe, {"1:200:A>G": 0.962}, _settings(), RN)[0]
    assert not r.v


def test_second_hit_flag_only_strong_with_v():
    pe = _pe({"GENEA": 3, "GENE_B": 1})                              # GENEA moderate, GENE_B strong
    res = _resource()
    rows, _ = score([_hit(200, cadd=25.0), _hit(300)], res, pe, None, _settings(), RN)
    flags = {(r.hit.call.pos, r.gene): r.second_hit for r in rows}
    assert flags == {(200, "GENEA"): False, (200, "GENE_B"): True, (300, "GENEA"): False, (300, "GENE_B"): False}


def test_summary_counts_linked_genes_exomiser_does_not_rank():
    hits = [_hit(200, cadd=25.0), _hit(1200, peaks=(PEAK_B,))]
    _, s = score(hits, _resource(), _pe({"GENEA": 1}), {}, _settings(), RN)
    assert s["linked_genes_with_hits"] == 2 and s["linked_genes_unranked"] == 1     # GENE_B is not in the ranking
    floor = _settings(gene_panel={"name": "P", "genes": ["GENE_B"], "phenotype_floor": "strong",
                                  "deep_second_hit": False})
    rows, s = score(hits, _resource(), _pe({"GENEA": 1}), {}, floor, RN)
    assert {r.pheno for r in rows if r.gene == "GENE_B"} == {"strong"}             # panel floor raised it ...
    assert s["linked_genes_unranked"] == 1                                         # ... but Exomiser did not rank it
    _, s = score(hits, _resource(), None, {}, _settings(), RN)
    assert s["linked_genes_with_hits"] == 2 and s["linked_genes_unranked"] == 0     # Exomiser off: nothing to miss
    _, s = score([], _resource(), _pe({}), {}, _settings(), RN)
    assert s["linked_genes_with_hits"] == 0 and s["linked_genes_unranked"] == 0


def test_score_without_exomiser_is_all_r5():
    rows, summary = score([_hit(200, cadd=25.0)], _resource(), None, {}, _settings(), RN)
    assert {r.tier for r in rows} == {"R5"} and summary["pheno_genes_in_universe"] == 0


def test_trait_null_reproducible_with_seed():
    rows = [SldscRow(f"T{i}", f"TF{i % 3}", 1e-9, 0.01 if i % 2 else 0.5) for i in range(10)]
    traits = [f"T{i}" for i in range(10)]
    a = trait_null(rows, traits, 2, 0.10, 25, 11)
    assert a == trait_null(rows, list(reversed(traits)), 2, 0.10, 25, 11) and len(a) == 25
    assert set(a) <= {0, 1, 2}


def test_summary_null_and_traits():
    res = _resource(sldsc=[SldscRow("T1", "TF1", 1e-9, 0.004), SldscRow("T2", "TF2", 1e-9, 0.5),
                           SldscRow("T3", "TF2", 1e-9, 0.001)])
    _, s = score([], res, _pe({"GENEA": 1}), {}, _settings(), RN)
    assert s["declared_traits"] == ["T1:trait one", "T2:trait two"]
    assert s["flagged_tfs_T1"] == ["TF1"] and s["flagged_tfs_T2"] == []
    assert s["observed_flagged_tfs"] == 1 and s["null_draws"] == 30
    # every 2-trait draw from {T1, T2, T3} flags >= 1 TF, so the null never falls below the observed value
    assert s["null_p"] == pytest.approx(1.0) and s["pheno_genes_in_universe"] == 1


def test_write_regulatory_columns(tmp_path):
    pe = _pe({"GENEA": 1})
    res = _resource(regulons={"TF2": frozenset({"GENEA"})})
    rows, summary = score([_hit(200, cadd=25.0)], res, pe, {"1:200:A>G": 0.97}, _settings(), RN)
    path = tmp_path / "x.regulatory.tsv"
    write_regulatory(rows, path, "chr")
    lines = path.read_text().splitlines()
    assert lines[0].split("\t") == ("chrom pos ref alt gene tier pheno panel second_hit peak popmax cadd remm v n_tfs "
                                    "sldsc_tfs regulon_tfs consequence").split()
    rec = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    assert rec["chrom"] == "chr1" and rec["pos"] == "200" and rec["gene"] == "GENEA" and rec["tier"] == "R1"
    assert rec["sldsc_tfs"] == "TF1(T1:0.004)" and rec["second_hit"] == "yes" and rec["v"] == "yes"
    assert rec["cadd"] == "25" and rec["remm"] == "0.97" and rec["n_tfs"] == "2" and rec["panel"] == ""
    assert rec["consequence"] == "intron_variant" and rec["regulon_tfs"] == ""
    write_regulatory(rows, path, "")
    assert path.read_text().splitlines()[1].startswith("1\t200\t")
    sp = tmp_path / "x.summary.tsv"
    write_summary(summary, sp)
    kv = dict(l.split("\t") for l in sp.read_text().splitlines())
    assert kv["declared_traits"] == "T1:trait one;T2:trait two" and kv["rows_R1"] == "1"


def test_sldsc_cell_joins_traits_and_tfs(tmp_path):
    res = _resource(sldsc=[SldscRow("T1", "TF1", 1e-9, 0.004), SldscRow("T2", "TF1", 1e-9, 0.02),
                           SldscRow("T2", "TF2", 1e-9, 0.05)],
                    regulons={"TF1": frozenset({"GENEA"})},                  # p = 1/21, the only TF tested
                    universe=("GENEA", "GENE_B", *(f"G{i}" for i in range(19))))
    rows, _ = score([_hit(200)], res, _pe({"GENEA": 1}), {}, _settings(), RN)
    path = tmp_path / "y.tsv"
    write_regulatory(rows, path, "chr")
    rec = dict(zip(*(l.split("\t") for l in path.read_text().splitlines()[:2])))
    assert rec["sldsc_tfs"] == "TF1(T1:0.004,T2:0.02);TF2(T2:0.05)" and rec["cadd"] == ""
    assert rec["regulon_tfs"] == "TF1(0.04762)"


def test_load_resource_from_build(tmp_path):
    make_scdori(tmp_path / "scdori")
    build(tmp_path / "scdori", [make_sldsc(tmp_path, "srcA")], "M1", tmp_path / "b", None)
    res = load_resource(tmp_path / "b")
    assert res.links == {PEAK_A: ("GENEA",), PEAK_B: ("GENEA", "GENE_B")}
    assert res.peak_tfs == {PEAK_A: ("TF1",)}
    assert res.regulons == {"TF1": frozenset({"GENEA"})}
    assert res.universe == frozenset({"GENEA", "GENE_B", "GENEC"})
    assert set(res.traits) == {"T1", "T2"}
    t1 = {r.tf: r for r in res.sldsc if r.trait == "T1"}
    assert t1["TF1"].coef == pytest.approx(1e-9) and t1["TF1"].fdr == pytest.approx(2e-3)
    assert res.index.find("1", 1100, 1100) == [PEAK_B]
    (tmp_path / "b" / "regulons.tsv").unlink()
    with pytest.raises(ConfigError, match="regulons.tsv"):
        load_resource(tmp_path / "b")


# ----------------------------------------------------------------------------- Task 6: second-hit hook

from dataclasses import replace

from prioritise_helpers import make_cand

from vep_pipeline.candidates import build_alleles
from vep_pipeline.evidence import REMM_SECOND_HIT, variant_level
from vep_pipeline.regnet import RegRow, add_second_hits

S_TF1 = (("TF1", (("T1", 0.004),)),)


def _row(hit, gene="GENEA", peak=PEAK_A, tier="R1", sldsc_tfs=S_TF1, second_hit=True, remm=None):
    return RegRow(hit, gene, peak, "strong", "", True, tier, ("TF1", "TF2"), sldsc_tfs, (), remm, second_hit)


def test_add_second_hits_builds_evidence_for_unsymbolled_variant():
    entries = ({"Consequence": "intergenic_variant"},
               {"Consequence": "regulatory_region_variant", "ClinVar_CLNSIG": "Likely_benign",
                "ClinVar_CLNREVSTAT": "criteria_provided&_single_submitter"})
    hit = RegHit(make_call(200), entries, (PEAK_A,), 0.003, 22.5)
    (c,) = add_second_hits([], [_row(hit)], {"1:200:A>G": 0.5}, POPS)
    assert c.call is hit.call
    ev = c.ev
    assert (ev.symbol, ev.consequence, ev.hgvs, ev.hgvsc) == ("GENEA", "intergenic_variant&regulatory_region_variant",
                                                             PEAK_A, "")
    assert ev.popmax == 0.003 and ev.cadd is None and ev.remm is None             # link-only: no CADD, no ReMM
    assert (ev.clinvar, ev.stars) == ("blb", 1)
    assert (ev.revel, ev.am, ev.splice, ev.loeuf, ev.phaplo, ev.nmd_escape, ev.noncoding_signal) == (
        None, None, None, None, None, False, False)
    assert ev.grn == "peak chr1:100-600 -> GENEA; R1; S-LDSC TFs TF1"
    assert ev.grn_classes == frozenset({"unsupervised_predictor"})                # CADD 22.5; the row has no ReMM
    (c,) = add_second_hits([], [_row(hit, tier="R2", sldsc_tfs=())], None, POPS)
    assert c.ev.grn == "peak chr1:100-600 -> GENEA; R2" and c.ev.remm is None and c.ev.clinvar == "blb"


def _link_only(clnsig="", revstat="", cadd=None, remm=None):
    entry = {"Consequence": "intergenic_variant", "CADD_PHRED": str(cadd) if cadd else ""}
    if clnsig:
        entry.update(ClinVar_CLNSIG=clnsig, ClinVar_CLNREVSTAT=revstat)
    hit = RegHit(make_call(200), (entry,), (PEAK_A,), 0.0, cadd)
    (c,) = add_second_hits([], [_row(hit, tier="R2", sldsc_tfs=(), remm=remm)], {"1:200:A>G": remm} if remm else {},
                           POPS)
    return c.ev


GRN_REASON = ["regulatory network: peak chr1:100-600 -> GENEA; R2 (second hit only)"]


def test_link_only_cadd_and_remm_do_not_survive_the_network_exclusion():
    ev = _link_only(cadd=30.0, remm=0.99)
    assert variant_level(ev, None) == ("supporting", True, GRN_REASON)
    assert variant_level(ev, None, frozenset({"regulatory_network"}))[0] is None


UNSUP, CLIN = "unsupervised_predictor", "clinically_trained_predictor"


def test_link_only_remm_only_v_is_dropped_when_remm_is_excluded():
    ev = _link_only(cadd=15.0, remm=0.99)                                       # CADD below 20: V rests on ReMM
    assert ev.grn_classes == frozenset({CLIN})
    assert variant_level(ev, None) == ("supporting", True, GRN_REASON)
    assert variant_level(ev, None, frozenset({UNSUP})) == ("supporting", True, GRN_REASON)
    assert variant_level(ev, None, frozenset({CLIN}))[0] is None


def test_link_only_cadd_only_v_is_dropped_when_cadd_is_excluded():
    ev = _link_only(cadd=25.0, remm=0.5)
    assert ev.grn_classes == frozenset({UNSUP})
    assert variant_level(ev, None, frozenset({CLIN})) == ("supporting", True, GRN_REASON)
    assert variant_level(ev, None, frozenset({UNSUP}))[0] is None


def test_link_only_with_both_signals_survives_one_exclusion():
    ev = _link_only(cadd=20.0, remm=REMM_SECOND_HIT)                             # both exactly at the cut-off
    assert ev.grn_classes == frozenset({UNSUP, CLIN})
    for one in (UNSUP, CLIN):
        assert variant_level(ev, None, frozenset({one})) == ("supporting", True, GRN_REASON)
    assert variant_level(ev, None, frozenset({UNSUP, CLIN}))[0] is None


def test_link_only_clinvar_plp_stays_second_hit_only():
    ev = _link_only("Pathogenic", "criteria_provided&_multiple_submitters&_no_conflicts", cadd=25.0)
    assert (ev.clinvar, ev.stars) == (None, 0)
    assert variant_level(ev, None) == ("supporting", True, GRN_REASON)


def test_link_only_clinvar_blb_vetoes():
    ev = _link_only("Likely_benign", "criteria_provided&_single_submitter")
    assert variant_level(ev, None)[0] is None


def test_gene_entry_path_keeps_its_own_evidence_and_remm():
    entries = ({"Consequence": "intron_variant", "SYMBOL": "GENEA", "Feature": "E1", "CADD_PHRED": "30"},)
    hit = RegHit(make_call(200), entries, (PEAK_A,), 0.0, 30.0)
    (c,) = add_second_hits([], [_row(hit, remm=0.99)], {"1:200:A>G": 0.99}, POPS)
    assert c.ev.cadd == 30.0 and c.ev.remm == 0.99 and c.ev.grn_classes == frozenset({UNSUP, CLIN})


def test_add_second_hits_uses_gene_entries_when_the_variant_has_them():
    entries = ({"Consequence": "intron_variant", "SYMBOL": "GENEA", "Feature": "E1", "HGVSc": "E1:c.10+50A>G",
                "CADD_PHRED": "21", "gnomADg_NFE_AF": "0.002"},
               {"Consequence": "upstream_gene_variant", "SYMBOL": "OTHER", "Feature": "E2", "CADD_PHRED": "21"})
    hit = RegHit(make_call(200), entries, (PEAK_A,), 0.002, 21.0)
    (c,) = add_second_hits([], [_row(hit, sldsc_tfs=S_TF1 + (("TF2", (("T2", 0.05),)),))], {}, POPS)
    assert (c.ev.symbol, c.ev.consequence, c.ev.hgvs) == ("GENEA", "intron_variant", "E1:c.10+50A>G")
    assert c.ev.popmax == pytest.approx(0.002) and c.ev.cadd == 21.0
    assert c.ev.grn == "peak chr1:100-600 -> GENEA; R1; S-LDSC TFs TF1,TF2"


def test_add_second_hits_replaces_existing_candidate_no_duplicate_piece():
    old = make_cand(200, "GENEA", consequence="intron_variant", cadd=26.0)
    other = make_cand(200, "OTHER", consequence="intron_variant", cadd=26.0)
    kept = make_cand(300, "GENEA", consequence="intron_variant", remm=0.99)
    cands = [old, other, kept]
    rows = [_row(_hit(200), remm=0.97), _row(_hit(300, cadd=21.0), tier="R2", sldsc_tfs=(), remm=0.10)]
    out = add_second_hits(cands, rows, {"1:200:A>G": 0.97, "1:300:A>G": 0.10}, POPS)
    assert cands == [old, other, kept] and old.ev.grn == ""                     # input untouched
    assert [c.call.key for c in out] == ["1:200:A>G", "1:200:A>G", "1:300:A>G"]
    assert out[0].ev == replace(old.ev, remm=0.97, grn="peak chr1:100-600 -> GENEA; R1; S-LDSC TFs TF1",
                                grn_classes=frozenset({CLIN}))                  # from the row: hit CADD None, ReMM 0.97
    assert out[1] is other
    assert out[2].ev == replace(kept.ev, grn="peak chr1:100-600 -> GENEA; R2",   # its own ReMM is kept
                                grn_classes=frozenset({UNSUP}))                  # the row's CADD 21, ReMM 0.10
    alleles = build_alleles(out, {}, _settings(), "MALE")["GENEA"]
    assert [len(a.pieces) for a in alleles] == [1, 1]
    assert "regulatory network: peak chr1:100-600 -> GENEA; R1; S-LDSC TFs TF1 (second hit only)" in alleles[0].reasons


def test_only_second_hit_rows_added():
    old = make_cand(200, "GENEA", consequence="intron_variant", cadd=26.0)
    rows = [_row(_hit(200), second_hit=False), _row(_hit(300), second_hit=False),
            _row(_hit(400), gene="GENE_B"), _row(_hit(500)),
            _row(_hit(400), gene="GENE_B", tier="R3", sldsc_tfs=())]           # same (key, gene) again: first wins
    out = add_second_hits([old], rows, {}, POPS)
    assert out[0] is old
    assert [(c.call.pos, c.ev.symbol, c.ev.grn) for c in out[1:]] == [
        (400, "GENE_B", "peak chr1:100-600 -> GENE_B; R1; S-LDSC TFs TF1"),
        (500, "GENEA", "peak chr1:100-600 -> GENEA; R1; S-LDSC TFs TF1")]
    assert add_second_hits([], [], {}, POPS) == []
