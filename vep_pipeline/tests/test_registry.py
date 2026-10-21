import pytest
from pathlib import Path
from vep_pipeline.registry import (load_registry, enabled_plugins, write_fields_tsv, RUNNABLE,
                                   fill_prep, Download, FileSpec, PluginSpec)
from vep_pipeline.config import ConfigError

REG = Path(__file__).parents[1] / "plugins.yaml"
PLUGDIR = Path("/fast/AG_Bunina/Berk/projects/mva_challenge_26/scripts/VEP_plugins")

def test_registry_covers_every_plugin_module():
    mods = {s.module for s in load_registry(REG).values()}
    assert {p.stem for p in PLUGDIR.glob("*.pm")} <= mods

def test_defaults_runnable_and_have_fields():
    for s in load_registry(REG).values():
        if s.default:
            assert s.category in RUNNABLE and s.output_fields and s.evidence_class, s.name

def test_output_fields_declared_in_plugin_source():
    for s in load_registry(REG).values():
        if s.default and s.category != "existing":
            src = (PLUGDIR / f"{s.module}.pm").read_text()
            for f in s.output_fields:
                stem = (f.split("_DS_")[0] if f.startswith("SpliceAI_pred")
                        else f[len("SpliceAI2_"):] if f.startswith("SpliceAI2_") else f)   # SpliceAI2.pm adds $PREFIX
                assert stem in src, (s.name, f)

def test_render_public_plugin(minimal_cfg):
    minimal_cfg["plugins"] = {"AlphaMissense": {"enabled": True}}
    (rp,) = enabled_plugins(minimal_cfg, load_registry(REG))
    exp = Path(minimal_cfg["data"]["plugin_data_dir"]) / "AlphaMissense/AlphaMissense_hg38.tsv.gz"
    assert rp.render("21", Path("/w")) == f"AlphaMissense,file={exp}"

def test_g2p_per_chunk_paths(minimal_cfg):
    minimal_cfg["plugins"] = {"G2P": {"enabled": True}}
    (rp,) = enabled_plugins(minimal_cfg, load_registry(REG))
    a = rp.render("1", Path("/w/1")); b = rp.render("2", Path("/w/2"))
    assert "log_dir=/w/1/g2p_log" in a and "g2p.1.txt" in a and "g2p.2.html" in b
    assert f"txt_report={minimal_cfg['output']['dir']}/reports/g2p.1.txt" in a

def test_user_path_overrides(minimal_cfg):
    minimal_cfg["plugins"] = {"CADD": {"enabled": True, "snv": "/x/s.tsv.gz", "indels": "/x/i.tsv.gz"}}
    (rp,) = enabled_plugins(minimal_cfg, load_registry(REG))
    assert rp.render("1", Path("/w")) == "CADD,snv=/x/s.tsv.gz,indels=/x/i.tsv.gz"

def test_manual_without_file_raises(minimal_cfg):
    minimal_cfg["plugins"] = {"dbNSFP": {"enabled": True}}
    with pytest.raises(ConfigError, match="dbnsfp.org"):
        enabled_plugins(minimal_cfg, load_registry(REG))

def test_unsupported_refused(minimal_cfg):
    minimal_cfg["plugins"] = {"MTR": {"enabled": True}}
    with pytest.raises(ConfigError, match="grch37_only"):
        enabled_plugins(minimal_cfg, load_registry(REG))

def test_conflict(minimal_cfg):
    minimal_cfg["plugins"] = {"SpliceAI": {"enabled": True},
                              "SpliceAI_Illumina": {"enabled": True, "snv": "/a", "indel": "/b"}}
    with pytest.raises(ConfigError, match="conflict"):
        enabled_plugins(minimal_cfg, load_registry(REG))

def test_escape_hatch(minimal_cfg):
    minimal_cfg["plugins"] = {"Blosum62": {"enabled": True, "args": "Blosum62"}}
    (rp,) = enabled_plugins(minimal_cfg, load_registry(REG))
    assert rp.render("1", Path("/w")) == "Blosum62"

def test_unknown_plugin(minimal_cfg):
    minimal_cfg["plugins"] = {"Nope": {"enabled": True}}
    with pytest.raises(ConfigError, match="unknown plugin"):
        enabled_plugins(minimal_cfg, load_registry(REG))

def test_fields_tsv(minimal_cfg, tmp_path):
    minimal_cfg["plugins"] = {"REVEL": {"enabled": True}}
    write_fields_tsv(enabled_plugins(minimal_cfg, load_registry(REG)), tmp_path / "f.tsv")
    rows = dict(l.split("\t")[:2] for l in (tmp_path / "f.tsv").read_text().splitlines()[1:])
    assert rows["REVEL"] == "clinically_trained_predictor" and rows["CLIN_SIG"] == "known_variant"

# ---- additional coverage

def test_mini_registry_fixture_is_specs(mini_registry, minimal_cfg):
    assert set(mini_registry) == {"Toy", "Nodata"}
    minimal_cfg["plugins"] = {"Toy": True, "Nodata": False}
    (rp,) = enabled_plugins(minimal_cfg, mini_registry)
    assert rp.spec.name == "Toy"
    assert rp.render("1", Path("/w")).startswith("Toy,file=")

def test_loader_maps_as_and_defaults(mini_registry):
    f = mini_registry["Toy"].files[0]
    assert isinstance(f, FileSpec) and f.key == "file" and f.index == "tbi"
    assert f.downloads == (Download("file:///dev/null", "Toy/toy.tsv.gz", None),)

def test_curated_downloads_exact():
    reg = load_registry(REG)
    (d,) = reg["REVEL"].files[0].downloads
    assert d.as_ == "REVEL/revel.zip" and d.size == 667102707
    assert reg["REVEL"].files[0].dest == "REVEL/new_tabbed_revel_grch38.tsv.gz"
    assert len(reg["MechPredict"].files[0].downloads) == 3
    sai = reg["SpliceAI"].files[0]
    assert sai.downloads[0].size == 28643031420
    assert sai.index_url == sai.downloads[0].url + ".tbi"
    assert sai.checksums_url.endswith("release-116/variation/SpliceAI/CHECKSUMS")
    assert reg["G2P"].files[0].validate["min_lines"] == 1000

def test_controls_follow_schema():
    for s in load_registry(REG).values():
        c = s.control
        if not c:
            continue
        assert c["source"] in {"data", "vcf", "clinvar", "header"}, s.name
        if c["source"] == "data":
            assert len(c["cols"]) == 4 and all(isinstance(i, int) and i >= 1 for i in c["cols"])
        if c["source"] == "clinvar":
            assert isinstance(c["mc"], str) and c["mc"]

def test_bad_control_rejected(tmp_path):
    p = tmp_path / "r.yaml"
    p.write_text("X: {module: X, category: public, control: {source: data}}\n")
    with pytest.raises(ConfigError, match="control"):
        load_registry(p)

def test_bad_category_rejected(tmp_path):
    p = tmp_path / "r.yaml"
    p.write_text("X: {module: X, category: bogus}\n")
    with pytest.raises(ConfigError, match="category"):
        load_registry(p)

def test_unknown_registry_key_rejected(tmp_path):
    p = tmp_path / "r.yaml"
    p.write_text("X: {module: X, category: public, argz: foo}\n")
    with pytest.raises(ConfigError, match="argz"):
        load_registry(p)

def test_nonrunnable_needs_allow_unsupported(minimal_cfg):
    reg = load_registry(REG)
    minimal_cfg["plugins"] = {"MTR": {"enabled": True, "args": "MTR,file=/x"}}
    with pytest.raises(ConfigError, match="allow_unsupported"):
        enabled_plugins(minimal_cfg, reg)
    minimal_cfg["plugins"] = {"MTR": {"enabled": True, "args": "MTR,file=/x", "allow_unsupported": True}}
    (rp,) = enabled_plugins(minimal_cfg, reg)
    assert rp.render("1", Path("/w")) == "MTR,file=/x"

def test_public_without_recipe_needs_args(minimal_cfg):
    minimal_cfg["plugins"] = {"GWAS": {"enabled": True}}
    with pytest.raises(ConfigError, match="args"):
        enabled_plugins(minimal_cfg, load_registry(REG))

def test_unknown_entry_key_rejected(minimal_cfg):
    minimal_cfg["plugins"] = {"AlphaMissense": {"enabled": True, "bogus": 1}}
    with pytest.raises(ConfigError, match="bogus"):
        enabled_plugins(minimal_cfg, load_registry(REG))

def test_disabled_and_null_values(minimal_cfg):
    minimal_cfg["plugins"] = {"AlphaMissense": False,
                              "SpliceAI_Illumina": {"enabled": False, "snv": None, "indel": None}}
    assert enabled_plugins(minimal_cfg, load_registry(REG)) == []

def test_manual_null_file_raises(minimal_cfg):
    minimal_cfg["plugins"] = {"dbNSFP": {"enabled": True, "file": None}}
    with pytest.raises(ConfigError, match="dbnsfp.org"):
        enabled_plugins(minimal_cfg, load_registry(REG))

def test_dbnsfp_columns_and_fields(minimal_cfg, tmp_path):
    minimal_cfg["plugins"] = {"dbNSFP": {"enabled": True, "file": "/d/dbNSFP5.gz",
                                         "columns": ["VEST4_score", "ESM1b_score"]}}
    plugins = enabled_plugins(minimal_cfg, load_registry(REG))
    assert plugins[0].render("1", Path("/w")) == "dbNSFP,/d/dbNSFP5.gz,VEST4_score,ESM1b_score"
    write_fields_tsv(plugins, tmp_path / "f.tsv")
    rows = dict(l.split("\t")[:2] for l in (tmp_path / "f.tsv").read_text().splitlines()[1:])
    assert rows["VEST4_score"] == "clinically_trained_predictor"
    assert rows["ESM1b_score"] == "unsupervised_predictor"
    assert "BayesDel_noAF_score" not in rows

def test_spliceai_split_output_and_paths(minimal_cfg):
    minimal_cfg["plugins"] = {"SpliceAI": True}
    (rp,) = enabled_plugins(minimal_cfg, load_registry(REG))
    out = rp.render("1", Path("/w"))
    assert out.startswith("SpliceAI,snv=") and out.endswith(",split_output=1")
    # SpliceAI.pm add_file()s `indel` unconditionally: snv alone dies with "No file specified" (Task 4, empirical)
    assert f",indel={rp.paths['indel']}," in out
    assert rp.paths["indel"] == Path(minimal_cfg["data"]["plugin_data_dir"]) / "SpliceAI/empty_indel.vcf.gz"

def test_spliceai_empty_indel_is_generated_not_downloaded():
    snv, indel = load_registry(REG)["SpliceAI"].files
    assert (snv.key, indel.key) == ("snv", "indel")
    assert indel.downloads == () and indel.index == "tbi" and indel.manual is None and indel.prep

def test_fields_tsv_header_and_sources(minimal_cfg, tmp_path):
    minimal_cfg["plugins"] = {"REVEL": True, "SpliceRegion": True}
    write_fields_tsv(enabled_plugins(minimal_cfg, load_registry(REG)), tmp_path / "f.tsv")
    lines = (tmp_path / "f.tsv").read_text().splitlines()
    assert lines[0] == "field\tevidence_class\tsource"
    rows = {l.split("\t")[0]: l.split("\t") for l in lines[1:]}
    assert rows["REVEL"][2] == "REVEL" and rows["SIFT"][1] == "clinically_trained_predictor"
    assert rows["ENCODE_cCRE"][1] == "regulatory_overlap" and rows["GENE_PHENO"][1] == "gene_phenotype"
    assert len(rows) == len(lines) - 1

def test_conflicts_with_enforced_in_mini(tmp_path, minimal_cfg):
    p = tmp_path / "r.yaml"
    p.write_text("A: {module: A, category: no_data, args: A, conflicts_with: [B]}\n"
                 "B: {module: B, category: no_data, args: B}\n")
    minimal_cfg["plugins"] = {"A": True, "B": True}
    with pytest.raises(ConfigError, match="conflict"):
        enabled_plugins(minimal_cfg, load_registry(p))

def test_fill_prep_only_known_placeholders():
    cmd = "cd {dir} && awk '{print $1}' {path} {empty_indel}"
    out = fill_prep(cmd, Path("/d/x/f.gz"), empty_indel="/e.vcf.gz")
    assert out == "cd /d/x && awk '{print $1}' /d/x/f.gz /e.vcf.gz"

def test_mechpredict_prep_keeps_awk_braces():
    (prep,) = load_registry(REG)["MechPredict"].files[0].prep
    out = fill_prep(prep, Path("/p/MechPredict/MechPredict_input.tsv"))
    assert "{print $1" in out and "cd /p/MechPredict" in out and "{dir}" not in out

def test_list_plugins_handler(minimal_cfg, capsys):
    from vep_pipeline.cli import HANDLERS
    minimal_cfg["plugins"] = {"AlphaMissense": {"enabled": True}, "Nope": {"enabled": True}}
    assert HANDLERS["list-plugins"](minimal_cfg, None) == 0
    out = capsys.readouterr().out
    rows = {l.split()[0]: l for l in out.splitlines()[1:]}
    assert len(rows) >= 83
    assert " yes " in rows["AlphaMissense"] and "public" in rows["AlphaMissense"]
    assert "grch37_only" in rows["MTR"] and "Paralogues" in rows
    assert "unknown plugin" in rows["Nope"]

def test_spliceai2_renders_user_directories(minimal_cfg):
    minimal_cfg["plugins"] = {"SpliceAI2": {"enabled": True, "snv": "/d/snv", "indel": "/d/indel"}}
    (rp,) = enabled_plugins(minimal_cfg, load_registry(REG))
    assert rp.render("1", Path("/w")) == "SpliceAI2,snv=/d/snv,indel=/d/indel"
    assert {fs.key: fs.glob for fs in rp.spec.files} == {"snv": "chr*.tsv.gz", "indel": "chr*.tsv.gz"}

def test_spliceai2_without_directories_quotes_instructions(minimal_cfg):
    minimal_cfg["plugins"] = {"SpliceAI2": {"enabled": True}}
    with pytest.raises(ConfigError, match="Hugging Face"):
        enabled_plugins(minimal_cfg, load_registry(REG))

def test_spliceai2_fields_are_the_plugin_defaults():
    src = (PLUGDIR / "SpliceAI2.pm").read_text()
    for f in load_registry(REG)["SpliceAI2"].output_fields:
        assert f.startswith("SpliceAI2_") and f[len("SpliceAI2_"):] in src, f

def test_glob_entry_must_be_user_supplied(tmp_path):
    import yaml
    p = tmp_path / "r.yaml"
    p.write_text(yaml.safe_dump({"T": {"module": "T", "category": "manual", "files": {"d": {
        "dest": "T/d", "glob": "*.gz", "manual": "get it", "downloads": [{"url": "file:///dev/null"}]}}}}))
    with pytest.raises(ConfigError, match="glob"):
        load_registry(p)

def test_filespec_members(tmp_path):
    for n in ("chr2.tsv.gz", "chr10.tsv.gz", "chr2.tsv.gz.tbi", "notes.txt"):
        (tmp_path / n).write_text("x")
    assert [p.name for p in FileSpec(key="d", dest="d", glob="chr*.tsv.gz").members(tmp_path)] == \
        ["chr10.tsv.gz", "chr2.tsv.gz"]
    assert FileSpec(key="f", dest="f").members(tmp_path / "notes.txt") == [tmp_path / "notes.txt"]
