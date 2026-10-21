import pytest

from vep_pipeline.config import ConfigError
from vep_pipeline.rank_settings import DEFAULT_TIERS, exomiser_settings, phasing_settings, prioritise_settings


def test_defaults_when_sections_absent(minimal_cfg):
    s = prioritise_settings(minimal_cfg)
    assert s["top_n"] == 10
    assert s["output"].endswith("/out/S1.top10.csv")
    assert s["tiers"] == DEFAULT_TIERS
    assert s["secondary"]["max_rows"] == 2
    assert s["af"]["populations"] == ["AFR", "AMR", "EAS", "NFE", "SAS"]


def test_top_n_override_sets_the_default_filename(minimal_cfg):
    assert prioritise_settings(minimal_cfg, top_n=5)["output"].endswith("S1.top5.csv")


def test_explicit_output_wins(minimal_cfg):
    assert prioritise_settings(minimal_cfg, output="/x/y.csv")["output"] == "/x/y.csv"


def test_unknown_key_is_an_error(minimal_cfg):
    minimal_cfg["prioritise"] = {"qc": {"min_gq": 20, "min_qual": 5}}
    with pytest.raises(ConfigError, match="min_qual"):
        prioritise_settings(minimal_cfg)


def test_tier_rules_are_validated(minimal_cfg):
    minimal_cfg["prioritise"] = {"tiers": [{"name": "T1", "epcr": 0.9, "rules": [["strong", "full", "great"]]},
                                           {"name": "T5", "epcr": 0.0, "rules": [["*", "*", "*"]]}]}
    with pytest.raises(ConfigError, match="bad rule"):
        prioritise_settings(minimal_cfg)


def test_last_tier_needs_catch_all(minimal_cfg):
    minimal_cfg["prioritise"] = {"tiers": [{"name": "T1", "epcr": 0.9, "rules": [["strong", "full", "strong"]]}]}
    with pytest.raises(ConfigError, match="catch-all"):
        prioritise_settings(minimal_cfg)


def test_exclude_classes_are_validated(minimal_cfg):
    minimal_cfg["prioritise"] = {"exclude_evidence_classes": ["regulatory_overlap"]}
    with pytest.raises(ConfigError, match="regulatory_overlap"):
        prioritise_settings(minimal_cfg)


def test_secondary_rows_cannot_exceed_top_n(minimal_cfg):
    minimal_cfg["prioritise"] = {"top_n": 1, "secondary": {"max_rows": 2}}
    with pytest.raises(ConfigError, match="max_rows"):
        prioritise_settings(minimal_cfg)


def test_exomiser_enabled_needs_hpo_or_phenopacket(minimal_cfg):
    with pytest.raises(ConfigError, match="exomiser.hpo"):
        exomiser_settings(minimal_cfg)


def test_exomiser_hpo_ids_are_checked(minimal_cfg):
    minimal_cfg["exomiser"] = {"hpo": ["HP:0001250", "seizure"]}
    with pytest.raises(ConfigError, match="seizure"):
        exomiser_settings(minimal_cfg)


def test_exomiser_disabled_needs_nothing(minimal_cfg):
    minimal_cfg["exomiser"] = {"enabled": False}
    assert exomiser_settings(minimal_cfg)["sex"] == "UNKNOWN_SEX"


def test_exomiser_results_dir_default(minimal_cfg):
    minimal_cfg["exomiser"] = {"hpo": ["HP:0001250"], "sex": "MALE"}
    ex = exomiser_settings(minimal_cfg)
    assert ex["results_dir"].endswith("/out/exomiser") and ex["sex"] == "MALE"


# ---- fix round 1: type and range validation -------------------------------------------------

CATCH_ALL = {"name": "T5", "epcr": 0.02, "rules": [["*", "*", "*"]]}


def _tier(**kw):
    t = {"name": "T1", "epcr": 0.9, "rules": [["strong", "full", "strong"]]}
    t.update(kw)
    return t


def test_literal_defaults(minimal_cfg):
    s = prioritise_settings(minimal_cfg)
    assert s["p_trans_unphased"] == 0.5
    assert s["primary_sum_cap"] == 0.95
    assert s["pheno_rank"] == {"strong": 20, "moderate": 200}
    assert s["af"]["max_dominant"] == 0.0001 and s["af"]["max_recessive"] == 0.01
    assert s["secondary"]["epcr"] == {"clinvar_2star": 0.30, "clinvar_1star": 0.20, "null_lof": 0.10}
    assert [t["epcr"] for t in s["tiers"]] == [0.90, 0.60, 0.30, 0.10, 0.02]


@pytest.mark.parametrize("bad", [True, "0.9", None, 1.5, -0.1])
def test_tier_epcr_must_be_a_number_in_0_1(minimal_cfg, bad):
    minimal_cfg["prioritise"] = {"tiers": [_tier(epcr=bad), CATCH_ALL]}
    with pytest.raises(ConfigError, match="T1"):
        prioritise_settings(minimal_cfg)


@pytest.mark.parametrize("bad", [None, "strong", [["strong", "full"]], [["strong", "full", "strong", "x"]], ["strong"]])
def test_tier_rules_must_be_three_item_lists(minimal_cfg, bad):
    minimal_cfg["prioritise"] = {"tiers": [_tier(rules=bad), CATCH_ALL]}
    with pytest.raises(ConfigError, match="T1"):
        prioritise_settings(minimal_cfg)


@pytest.mark.parametrize("bad", [True, "5", 2.0, 0])
def test_top_n_must_be_a_positive_int(minimal_cfg, bad):
    minimal_cfg["prioritise"] = {"top_n": bad, "secondary": {"max_rows": 0}}
    with pytest.raises(ConfigError, match="top_n"):
        prioritise_settings(minimal_cfg)


@pytest.mark.parametrize("bad", [True, "2", 2.0, -1])
def test_max_rows_must_be_an_int_in_range(minimal_cfg, bad):
    minimal_cfg["prioritise"] = {"secondary": {"max_rows": bad}}
    with pytest.raises(ConfigError, match="max_rows"):
        prioritise_settings(minimal_cfg)


def test_max_rows_error_shows_the_values(minimal_cfg):
    minimal_cfg["prioritise"] = {"top_n": 1, "secondary": {"max_rows": 2}}
    with pytest.raises(ConfigError, match="top_n=1, got 2"):
        prioritise_settings(minimal_cfg)


@pytest.mark.parametrize("key", ["p_trans_unphased", "primary_sum_cap"])
@pytest.mark.parametrize("bad", [True, "0.5", 0, 0.0, 1.01])
def test_fractions_are_numbers_in_unit_interval(minimal_cfg, key, bad):
    minimal_cfg["prioritise"] = {key: bad}
    with pytest.raises(ConfigError, match=key):
        prioritise_settings(minimal_cfg)


@pytest.mark.parametrize("bad", [True, "0.3", 1.2, -0.1])
def test_secondary_epcr_values_are_numbers_in_0_1(minimal_cfg, bad):
    minimal_cfg["prioritise"] = {"secondary": {"epcr": {"null_lof": bad}}}
    with pytest.raises(ConfigError, match="null_lof"):
        prioritise_settings(minimal_cfg)


@pytest.mark.parametrize("bad", [[], "AFR", ["AFR", "XYZ"], [None]])
def test_af_populations_must_be_gnomad_groups(minimal_cfg, bad):
    minimal_cfg["prioritise"] = {"af": {"populations": bad}}
    with pytest.raises(ConfigError, match="populations"):
        prioritise_settings(minimal_cfg)


def test_all_gnomad_v4_groups_are_accepted(minimal_cfg):
    groups = ["AFR", "AMR", "ASJ", "EAS", "FIN", "MID", "NFE", "SAS", "AMI", "REMAINING"]
    minimal_cfg["prioritise"] = {"af": {"populations": groups}}
    assert prioritise_settings(minimal_cfg)["af"]["populations"] == groups


@pytest.mark.parametrize("key", ["max_dominant", "max_recessive"])
@pytest.mark.parametrize("bad", [True, "0.1", 0, 1.5])
def test_af_cutoffs_are_numbers_in_unit_interval(minimal_cfg, key, bad):
    minimal_cfg["prioritise"] = {"af": {key: bad}}
    with pytest.raises(ConfigError, match=key):
        prioritise_settings(minimal_cfg)


@pytest.mark.parametrize("bad", [{"strong": 0}, {"strong": True}, {"moderate": "200"}, {"strong": 300, "moderate": 200}])
def test_pheno_rank_must_be_positive_ints_ordered(minimal_cfg, bad):
    minimal_cfg["prioritise"] = {"pheno_rank": bad}
    with pytest.raises(ConfigError, match="pheno_rank"):
        prioritise_settings(minimal_cfg)


def test_exomiser_hpo_must_be_a_list(minimal_cfg):
    minimal_cfg["exomiser"] = {"hpo": "HP:0001250"}
    with pytest.raises(ConfigError, match="exomiser.hpo"):
        exomiser_settings(minimal_cfg)


def test_exomiser_false_disables_it(minimal_cfg):
    minimal_cfg["exomiser"] = False
    ex = exomiser_settings(minimal_cfg)
    assert ex["enabled"] is False and ex["sex"] == "UNKNOWN_SEX"


def test_exomiser_true_means_defaults(minimal_cfg):
    minimal_cfg["exomiser"] = True
    with pytest.raises(ConfigError, match="exomiser.hpo"):
        exomiser_settings(minimal_cfg)


def test_exomiser_enabled_must_be_a_bool(minimal_cfg):
    minimal_cfg["exomiser"] = {"enabled": "no", "hpo": ["HP:0001250"]}
    with pytest.raises(ConfigError, match="exomiser.enabled"):
        exomiser_settings(minimal_cfg)


def test_prioritise_false_is_rejected(minimal_cfg):
    minimal_cfg["prioritise"] = False
    with pytest.raises(ConfigError, match="prioritise"):
        prioritise_settings(minimal_cfg)


def test_secondary_enabled_must_be_a_bool(minimal_cfg):
    minimal_cfg["prioritise"] = {"secondary": {"enabled": "yes"}}
    with pytest.raises(ConfigError, match="secondary.enabled"):
        prioritise_settings(minimal_cfg)


def test_missense_tool_default_and_validation(minimal_cfg):
    assert prioritise_settings(minimal_cfg)["missense_tool"] == "alphamissense"
    minimal_cfg["prioritise"] = {"missense_tool": "cadd"}
    with pytest.raises(ConfigError, match="missense_tool"):
        prioritise_settings(minimal_cfg)


def test_exclude_regions_default_and_validation(minimal_cfg):
    assert prioritise_settings(minimal_cfg)["exclude_regions"] == \
        [{"name": "MHC", "chrom": "6", "start": 28510120, "end": 33480577}]
    minimal_cfg["prioritise"] = {"exclude_regions": [{"name": "X", "chrom": "6", "start": 10, "end": 5}]}
    with pytest.raises(ConfigError, match="exclude_regions"):
        prioritise_settings(minimal_cfg)
    minimal_cfg["prioritise"] = {"exclude_regions": []}
    assert prioritise_settings(minimal_cfg)["exclude_regions"] == []


def _phasing_files(tmp_path, bam="s.bam", index="s.bam.bai", fai=True):
    (tmp_path / bam).write_text("")
    if index:
        (tmp_path / index).write_text("")
    (tmp_path / "ref.fa").write_text("")
    if fai:
        (tmp_path / "ref.fa.fai").write_text("")
    return {"enabled": True, "bam": str(tmp_path / bam), "reference": str(tmp_path / "ref.fa"), "env": "whatshap"}


def test_phasing_defaults(minimal_cfg):
    ph = phasing_settings(minimal_cfg)
    assert ph == {"enabled": False, "bam": None, "reference": None, "env": None,
                  "results_dir": str(minimal_cfg["output"]["dir"]) + "/phasing", "pad": 1000, "indels": True,
                  "mem": "16G", "time": "08:00:00"}


def test_phasing_disabled_checks_nothing(minimal_cfg):
    minimal_cfg["phasing"] = {"bam": "/no/such.bam"}
    assert phasing_settings(minimal_cfg)["enabled"] is False


@pytest.mark.parametrize("index", ["s.bam.bai", "s.bai"])
def test_phasing_enabled_accepts_either_bam_index(minimal_cfg, tmp_path, index):
    minimal_cfg["phasing"] = _phasing_files(tmp_path, index=index)
    assert phasing_settings(minimal_cfg)["enabled"] is True


def test_phasing_enabled_accepts_a_cram_index(minimal_cfg, tmp_path):
    minimal_cfg["phasing"] = _phasing_files(tmp_path, bam="s.cram", index="s.cram.crai")
    assert phasing_settings(minimal_cfg)["bam"].endswith("s.cram")


@pytest.mark.parametrize("change, match", [
    ({"bam": None}, "phasing.bam"),
    ({"reference": None}, "phasing.reference"),
    ({"env": None}, "phasing.env"),
    ({"bam": "/no/such.bam"}, "not found"),
    ({"reference": "/no/such.fa"}, "not found"),
    ({"pad": 0}, "pad"),
    ({"indels": "yes"}, "indels"),
    ({"enabled": "yes"}, "enabled"),
    ({"mem": 16}, "mem"),
])
def test_phasing_enabled_validation(minimal_cfg, tmp_path, change, match):
    minimal_cfg["phasing"] = {**_phasing_files(tmp_path), **change}
    with pytest.raises(ConfigError, match=match):
        phasing_settings(minimal_cfg)


def test_phasing_enabled_needs_the_bam_index(minimal_cfg, tmp_path):
    minimal_cfg["phasing"] = _phasing_files(tmp_path, index=None)
    with pytest.raises(ConfigError, match="index"):
        phasing_settings(minimal_cfg)


def test_phasing_enabled_needs_the_reference_index(minimal_cfg, tmp_path):
    minimal_cfg["phasing"] = _phasing_files(tmp_path, fai=False)
    with pytest.raises(ConfigError, match=r"\.fai"):
        phasing_settings(minimal_cfg)


def test_phasing_rejects_an_unknown_alignment_suffix(minimal_cfg, tmp_path):
    minimal_cfg["phasing"] = _phasing_files(tmp_path, bam="s.sam", index=None)
    with pytest.raises(ConfigError, match="BAM or CRAM"):
        phasing_settings(minimal_cfg)


def test_phase_sources_default_and_validation(minimal_cfg):
    assert prioritise_settings(minimal_cfg)["phase_sources"] == ["multiallelic", "ps", "pid", "whatshap"]
    minimal_cfg["prioritise"] = {"phase_sources": []}
    assert prioritise_settings(minimal_cfg)["phase_sources"] == []
    minimal_cfg["prioritise"] = {"phase_sources": ["ps", "bogus"]}
    with pytest.raises(ConfigError, match="phase_sources"):
        prioritise_settings(minimal_cfg)
    minimal_cfg["prioritise"] = {"phase_sources": "ps"}
    with pytest.raises(ConfigError, match="phase_sources"):
        prioritise_settings(minimal_cfg)


def test_gene_panel_defaults_are_off(minimal_cfg):
    assert prioritise_settings(minimal_cfg)["gene_panel"] == {
        "name": "", "genes": [], "phenotype_floor": None, "deep_second_hit": False}


def test_gene_panel_is_accepted(minimal_cfg):
    panel = {"name": "P", "genes": ["A", "B"], "phenotype_floor": "strong", "deep_second_hit": True}
    minimal_cfg["prioritise"] = {"gene_panel": panel}
    assert prioritise_settings(minimal_cfg)["gene_panel"] == panel


@pytest.mark.parametrize("panel", [
    {"name": 1},
    {"genes": "A"},
    {"genes": ["A", 2]},
    {"phenotype_floor": "weak"},
    {"phenotype_floor": "strong"},                            # a floor needs genes
    {"deep_second_hit": True},                                # so does the deep search
    {"deep_second_hit": "yes", "genes": ["A"]},
])
def test_gene_panel_is_validated(minimal_cfg, panel):
    minimal_cfg["prioritise"] = {"gene_panel": panel}
    with pytest.raises(ConfigError, match="prioritise.gene_panel"):
        prioritise_settings(minimal_cfg)


def test_gene_panel_unknown_key_is_an_error(minimal_cfg):
    minimal_cfg["prioritise"] = {"gene_panel": {"genez": ["A"]}}
    with pytest.raises(ConfigError, match="genez"):
        prioritise_settings(minimal_cfg)


def test_panel_partner_is_a_known_fit_in_t4():
    from vep_pipeline.rank_settings import FITS
    assert "panel_partner" in FITS
    t4 = next(t for t in DEFAULT_TIERS if t["name"] == "T4")
    assert ["*", "panel_partner", "strong"] in t4["rules"]


def test_alternative_pairs_default_is_one_per_gene(minimal_cfg):
    assert prioritise_settings(minimal_cfg)["alternative_pairs"] == {"max_per_gene": 1, "min_phenotype": "strong"}
    minimal_cfg["prioritise"] = {"alternative_pairs": {"max_per_gene": 3, "min_phenotype": "moderate"}}
    assert prioritise_settings(minimal_cfg)["alternative_pairs"] == {"max_per_gene": 3, "min_phenotype": "moderate"}


@pytest.mark.parametrize("ap", [{"max_per_gene": 0}, {"max_per_gene": 2.0}, {"max_per_gene": True},
                                {"min_phenotype": "high"}, {"min_phenotype": 1}])
def test_alternative_pairs_are_validated(minimal_cfg, ap):
    minimal_cfg["prioritise"] = {"alternative_pairs": ap}
    with pytest.raises(ConfigError, match="prioritise.alternative_pairs"):
        prioritise_settings(minimal_cfg)
