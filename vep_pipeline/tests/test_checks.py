from vep_pipeline import checks


def test_unknown_ini_key(tmp_path):
    p = tmp_path / "params.txt"; p.write_text("offline\ncache\nplugin\n")
    assert checks.validate_ini_keys("offline 1\nbogus 3\nplugin X\n", p) == ["bogus"]


def test_ini_inline_comment_rejected(tmp_path):
    p = tmp_path / "params.txt"; p.write_text("plugin\n")
    assert checks.validate_ini_keys("plugin X,a=1 # note\n", p) == ["inline-comment:plugin"]


# ----------------------------------------------------------------------------- Task 4 preflight checks

import subprocess
from pathlib import Path

import pytest
import yaml

from vep_pipeline import cli, download as dl
from vep_pipeline.registry import PluginSpec, ResolvedPlugin, enabled_plugins, load_registry

REG = Path(__file__).resolve().parents[1] / "plugins.yaml"
POS = 25585733


def _bgz_tabix(path, *tabix_args):
    subprocess.run(["bgzip", "-f", str(path)], check=True)
    subprocess.run(["tabix", "-f", *tabix_args, f"{path}.gz"], check=True)
    return Path(f"{path}.gz")


def test_contig_length_mismatch_is_error(tmp_path, vcf_factory, fasta_factory):
    vcf = vcf_factory(contigs={"1": 249250621}); fa = fasta_factory({"1": 248956422})
    r = checks.check_vcf_vs_fasta(vcf, fa, env="vep116")
    assert not r.ok and "assembly" in r.message


def test_chr_prefix_mapped(tmp_path, vcf_factory, fasta_factory):
    vcf = vcf_factory(contigs={"chr1": 1000}, records=[("chr1", 10, "A", "G")])
    fa = fasta_factory({"1": 1000}, seq={"1": "A" * 1000})
    assert checks.check_vcf_vs_fasta(vcf, fa, env="vep116").ok


def test_chrM_mapped_to_MT(tmp_path, vcf_factory, fasta_factory):
    vcf = vcf_factory(contigs={"chrM": 100}, records=[("chrM", 10, "A", "G")])
    fa = fasta_factory({"MT": 100}, seq={"MT": "A" * 100})
    assert checks.check_vcf_vs_fasta(vcf, fa, env="vep116").ok


def test_ref_mismatch_rate(tmp_path, vcf_factory, fasta_factory):
    recs = [("1", i, "A" if i > 3 else "C", "G") for i in range(1, 11)]
    vcf = vcf_factory(contigs={"1": 100}, records=recs); fa = fasta_factory({"1": 100}, seq={"1": "A" * 100})
    r = checks.check_vcf_vs_fasta(vcf, fa, env="vep116")
    assert not r.ok and "30.0%" in r.message and r.level == "error"


def test_ref_mismatch_below_threshold_passes(vcf_factory, fasta_factory):
    recs = [("1", i, "A" if i > 1 else "C", "G") for i in range(1, 101)]      # 1 %
    vcf = vcf_factory(contigs={"1": 200}, records=recs); fa = fasta_factory({"1": 200}, seq={"1": "A" * 200})
    assert checks.check_vcf_vs_fasta(vcf, fa, env="vep116").ok


def test_vcf_vs_fasta_never_raises(tmp_path, fasta_factory):
    r = checks.check_vcf_vs_fasta(tmp_path / "nope.vcf.gz", fasta_factory({"1": 10}), env="vep116")
    assert isinstance(r, checks.CheckResult) and not r.ok and r.level == "error"


# --- sample (f)

def test_sample_present(vcf_factory):
    vcf = vcf_factory(contigs={"1": 100}, samples=("S1",))
    assert checks.check_sample(vcf, "S1", env="vep116").ok


def test_sample_missing_is_error(vcf_factory):
    r = checks.check_sample(vcf_factory(contigs={"1": 100}, samples=("A", "B")), "S1", env="vep116")
    assert not r.ok and r.level == "error" and "S1" in r.message


def test_multi_sample_warns(vcf_factory):
    r = checks.check_sample(vcf_factory(contigs={"1": 100}, samples=("S1", "B")), "S1", env="vep116")
    assert not r.ok and r.level == "warn" and "subset" in r.message


# --- sex (g)

X_LEN, Y_LEN = 156040895, 57227415


def test_sex_xx(vcf_factory):
    recs = [("X", 5_000_000 + i, "A", "G", "0/1") for i in range(40)]
    r = checks.check_sex(vcf_factory({"X": X_LEN, "Y": Y_LEN}, recs), "S1", env="vep116")
    assert r.level == "info" and r.ok and "XX" in r.message


def test_sex_xy(vcf_factory):
    recs = [("X", 5_000_000 + i, "A", "G", "1/1") for i in range(40)]
    recs += [("Y", 3_000_000 + i, "A", "G", "1/1") for i in range(400)]
    r = checks.check_sex(vcf_factory({"X": X_LEN, "Y": Y_LEN}, recs), "S1", env="vep116")
    assert "XY" in r.message and "XX" not in r.message


def test_sex_par_excluded(vcf_factory):
    recs = [("X", 1_000_000 + i, "A", "G", "0/1") for i in range(40)]            # PAR1 only
    r = checks.check_sex(vcf_factory({"X": X_LEN, "Y": Y_LEN}, recs), "S1", env="vep116")
    assert "unclear" in r.message


# --- env leaks (h), disk (k)

def test_env_leaks_warn(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    for k in [k for k in __import__("os").environ if k.startswith("VEP_")]:
        monkeypatch.delenv(k)
    monkeypatch.setenv("VEP_DIR_PLUGINS", "/p")
    assert checks.check_env_leaks().ok
    monkeypatch.setenv("VEP_CACHE", "/x")
    r = checks.check_env_leaks()
    assert not r.ok and r.level == "warn" and "VEP_CACHE" in r.message and "VEP_DIR_PLUGINS" not in r.message
    monkeypatch.delenv("VEP_CACHE")
    (tmp_path / ".vep").mkdir(); (tmp_path / ".vep" / "vep.ini").write_text("offline 1\n")
    r = checks.check_env_leaks()
    assert not r.ok and "vep.ini" in r.message


def test_disk_space(testcache_cfg, vcf_factory, monkeypatch):
    testcache_cfg["input"]["vcf"] = str(vcf_factory({"1": 100}))
    assert checks.check_disk(testcache_cfg).ok
    monkeypatch.setattr(checks, "_free_bytes", lambda p: 1)
    r = checks.check_disk(testcache_cfg)
    assert not r.ok and r.level == "warn"


# --- files (a)

def test_files_check(testcache_cfg, mini_registry, tmp_path):
    testcache_cfg["plugins"] = {"Toy": {"enabled": True}}
    plugins = enabled_plugins(testcache_cfg, mini_registry)
    toy = plugins[0].paths["file"]
    r = checks.check_files(testcache_cfg, plugins)
    assert not r.ok and r.level == "error" and "Toy" in r.message
    toy.parent.mkdir(parents=True)
    src = toy.with_suffix("")
    src.write_text("1\t10\tA\tG\n")
    subprocess.run(["bgzip", "-f", str(src)], check=True)       # toy.tsv.gz, no index yet
    r = checks.check_files(testcache_cfg, plugins)
    assert not r.ok and ".tbi" in r.message
    subprocess.run(["tabix", "-f", "-s1", "-b2", "-e2", str(toy)], check=True)
    assert checks.check_files(testcache_cfg, plugins).ok


def test_files_check_flags_empty_file_and_missing_cache(testcache_cfg, mini_registry):
    testcache_cfg["plugins"] = {"Toy": {"enabled": True}}
    plugins = enabled_plugins(testcache_cfg, mini_registry)
    plugins[0].paths["file"].parent.mkdir(parents=True)
    plugins[0].paths["file"].write_text("")
    Path(f"{plugins[0].paths['file']}.tbi").write_text("x")
    assert "empty" in checks.check_files(testcache_cfg, plugins).message
    testcache_cfg["data"]["cache_dir"] = str(Path(testcache_cfg["data"]["root"]) / "nocache")
    assert "info.txt" in checks.check_files(testcache_cfg, plugins).message


def test_files_check_custom_track(testcache_cfg, tmp_path):
    testcache_cfg["custom"] = {"Trk": {"enabled": True, "file": str(tmp_path / "t.vcf.gz"), "format": "vcf"}}
    r = checks.check_files(testcache_cfg, [])
    assert not r.ok and "Trk" in r.message


# --- ini keys (b)

def test_ini_keys_check_passes_for_generated_ini(testcache_cfg, mini_registry, tmp_path):
    testcache_cfg["plugins"] = {"Toy": {"enabled": True}}
    plugins = enabled_plugins(testcache_cfg, mini_registry)
    r = checks.check_ini_keys(testcache_cfg, plugins, tmp_path / "w")
    assert r.ok, r.message


def test_ini_keys_check_detects_stale_params_file(testcache_cfg, tmp_path):
    names = [l for l in checks.PARAMS_FILE.read_text().splitlines() if l.strip()]
    stale = tmp_path / "params.txt"
    stale.write_text("\n".join(names + ["bogus_option"]) + "\n")
    r = checks.check_ini_keys(testcache_cfg, [], tmp_path / "w", params_file=stale)
    assert not r.ok and "bogus_option" in r.message


def test_ini_keys_check_reports_unknown_key(testcache_cfg, tmp_path):
    testcache_cfg["run"]["extra_flags"] = ["not_a_vep_flag"]
    r = checks.check_ini_keys(testcache_cfg, [], tmp_path / "w")
    assert not r.ok and "not_a_vep_flag" in r.message


def test_committed_params_file_matches_vep():
    assert checks.regenerate_params("vep116") == checks.committed_params()


# --- dbNSFP (j)

def test_dbnsfp_columns(testcache_cfg, tmp_path):
    src = tmp_path / "dbnsfp.tsv"
    src.write_text("#chr\tpos(1-based)\tref\talt\tVEST4_score\tESM1b_score\n1\t10\tA\tG\t0.1\t-1\n")
    gz = _bgz_tabix(src, "-s1", "-b2", "-e2")
    testcache_cfg["plugins"] = {"dbNSFP": {"enabled": True, "file": str(gz), "columns": ["VEST4_score", "ESM1b_score"]}}
    plugins = enabled_plugins(testcache_cfg, load_registry(REG))
    assert checks.check_dbnsfp(plugins, env="vep116").ok
    testcache_cfg["plugins"]["dbNSFP"]["columns"] = ["VEST4_score", "Nope_score"]
    r = checks.check_dbnsfp(enabled_plugins(testcache_cfg, load_registry(REG)), env="vep116")
    assert not r.ok and "Nope_score" in r.message and r.level == "error"


def test_dbnsfp_not_enabled_is_ok():
    r = checks.check_dbnsfp([], env="vep116")
    assert r.ok and "not enabled" in r.message


# --- prepare script (i)

def _have_env(name):
    out = subprocess.run(["micromamba", "env", "list"], capture_output=True, text=True).stdout
    return any(l.split()[:1] == [name] for l in out.splitlines())


needs_exomiser = pytest.mark.skipif(not _have_env("exomiser"), reason="exomiser env not installed")
PREPARE = Path(__file__).resolve().parents[2] / "exomiser" / "bin" / "prepare-patient-vcf.sh"


@needs_exomiser
def test_prepare_script_usage(testcache_cfg):
    testcache_cfg["input"]["normalize"] = True
    testcache_cfg["tools"]["prepare_vcf"] = str(PREPARE)
    assert checks.check_prepare_script(testcache_cfg).ok


@needs_exomiser
def test_prepare_script_broken(testcache_cfg, tmp_path):
    bad = tmp_path / "bad.sh"; bad.write_text("echo hello\n")
    testcache_cfg["input"]["normalize"] = True
    testcache_cfg["tools"]["prepare_vcf"] = str(bad)
    r = checks.check_prepare_script(testcache_cfg)
    assert not r.ok and r.level == "error"


def test_prepare_script_skipped_without_normalize(testcache_cfg):
    testcache_cfg["input"]["normalize"] = False
    testcache_cfg["tools"]["prepare_vcf"] = "/nonexistent.sh"
    r = checks.check_prepare_script(testcache_cfg)
    assert r.ok and "normalize" in r.message


# --- plugin load (c) on the VEP test cache

def _cache_plugin(name="SpliceRegion", fields=("SpliceRegion",)):
    return ResolvedPlugin(PluginSpec(name=name, module=name, category="no_data", output_fields=tuple(fields)),
                          {}, {}, name)


def test_plugin_load_failure_detected(tmp_path, testcache_cfg):
    bad = tmp_path / "plugins/Bad.pm"; bad.parent.mkdir(); bad.write_text("package Bad; die 'boom';\n1;\n")
    testcache_cfg["tools"]["plugins_dir"] = str(bad.parent)
    rp = ResolvedPlugin(PluginSpec(name="Bad", module="Bad", category="no_data", output_fields=("X",)), {}, {}, "Bad")
    r = checks.plugin_load_test(testcache_cfg, rp, tmp_path / "w")
    assert not r.ok and "Bad" in r.message


def test_plugin_load_ok(tmp_path, testcache_cfg):
    r = checks.plugin_load_test(testcache_cfg, _cache_plugin(), tmp_path / "w")
    assert r.ok and r.level == "error", r.message


def test_plugin_load_missing_output_field(tmp_path, testcache_cfg):
    r = checks.plugin_load_test(testcache_cfg, _cache_plugin(fields=("SpliceRegion", "NOT_A_FIELD")), tmp_path / "w")
    assert not r.ok and "NOT_A_FIELD" in r.message


def test_plugin_load_keeps_outputs_in_workdir(tmp_path, testcache_cfg):
    checks.plugin_load_test(testcache_cfg, _cache_plugin(), tmp_path / "w")
    assert (tmp_path / "w").is_dir() and any((tmp_path / "w").iterdir())
    assert not (Path(testcache_cfg["data"]["cache_dir"]) / "homo_sapiens" / "84_GRCh38" / "stats.html").exists()


def test_plugin_load_ignores_inherited_vep_env(tmp_path, testcache_cfg, monkeypatch):
    monkeypatch.setenv("VEP_CACHE", "1")
    monkeypatch.setenv("VEP_DIR_CACHE", str(tmp_path / "bogus"))
    assert checks.plugin_load_test(testcache_cfg, _cache_plugin(), tmp_path / "w").ok


def test_plugin_load_does_not_use_enabled_custom_tracks(tmp_path, testcache_cfg):
    testcache_cfg["custom"] = {"Gone": {"enabled": True, "file": str(tmp_path / "gone.vcf.gz"), "format": "vcf"}}
    assert checks.plugin_load_test(testcache_cfg, _cache_plugin(), tmp_path / "w").ok


# --- SpliceAI (d): empirical result, synthetic SpliceAI-format file in the test-cache region

@pytest.fixture
def spliceai_snv(tmp_path):
    src = tmp_path / "spliceai_snv.vcf"
    src.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3.1 variant annotation. Format: '
        'ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">\n'
        "##contig=<ID=21,length=46709983>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"21\t{POS}\t.\tC\tT\t.\t.\tSpliceAI=T|MRPL39|0.01|0.02|0.03|0.04|1|2|3|4\n")
    return _bgz_tabix(src, "-p", "vcf")


def _spliceai(cfg, snv):
    cfg["plugins"] = {"SpliceAI": {"enabled": True, "snv": str(snv)}}
    (rp,) = enabled_plugins(cfg, load_registry(REG))
    return rp


def test_spliceai_without_indel_fails_to_load(tmp_path, testcache_cfg, spliceai_snv):
    rp = _spliceai(testcache_cfg, spliceai_snv)
    assert checks.needs_empty_indel(testcache_cfg, rp, tmp_path / "w") is True


def test_spliceai_registry_form_loads_with_empty_indel(tmp_path, testcache_cfg, spliceai_snv):
    rp = _spliceai(testcache_cfg, spliceai_snv)
    assert "indel=" in rp.arg_template and rp.arg_template.startswith("SpliceAI,snv=")
    checks.make_empty_indel_vcf(rp.paths["indel"])
    r = checks.plugin_load_test(testcache_cfg, rp, tmp_path / "w")
    assert r.ok, r.message


def test_make_empty_indel_vcf_is_header_only_and_indexed(tmp_path):
    p = tmp_path / "d" / "empty.vcf.gz"
    checks.make_empty_indel_vcf(p)
    assert Path(f"{p}.tbi").is_file()
    text = subprocess.run(["bgzip", "-dc", str(p)], capture_output=True, text=True, check=True).stdout
    assert text.startswith("##fileformat=VCF") and all(l.startswith("#") for l in text.splitlines())


def test_registry_indel_prep_matches_make_empty_indel_vcf(testcache_cfg, spliceai_snv, tmp_path):
    rp = _spliceai(testcache_cfg, spliceai_snv)
    fs = next(f for f in rp.spec.files if f.key == "indel")
    assert fs.downloads == () and fs.dest == "SpliceAI/empty_indel.vcf.gz" and fs.index == "tbi"
    out = tmp_path / "prep" / "x.vcf.gz"
    out.parent.mkdir()
    dl.run_prep(fs.prep, out, "vep116")
    ref = tmp_path / "ref.vcf.gz"
    checks.make_empty_indel_vcf(ref)
    dc = lambda p: subprocess.run(["bgzip", "-dc", str(p)], capture_output=True, text=True, check=True).stdout
    assert dc(out) == dc(ref) and Path(f"{out}.tbi").is_file()


def test_spliceai_check_skipped_when_indel_in_args(tmp_path, testcache_cfg, spliceai_snv):
    rp = _spliceai(testcache_cfg, spliceai_snv)
    r = checks.check_spliceai_form(testcache_cfg, [rp], tmp_path / "w")
    assert r.ok and "indel" in r.message


# --- run_all_checks + cli

def _real_ref():
    return subprocess.run(["samtools", "faidx", str(Path(__file__).resolve().parents[2] / "ensembl-vep/t/testdata/cache/"
                           "homo_sapiens/84_GRCh38/test.fa"), f"21:{POS}-{POS}"], capture_output=True, text=True,
                          check=True).stdout.split("\n")[1]


@pytest.fixture
def full_cfg(testcache_cfg, vcf_factory):
    testcache_cfg["input"]["vcf"] = str(vcf_factory({"21": 46709983}, [("21", POS, _real_ref(), "T")]))
    testcache_cfg["plugins"] = {"SpliceRegion": True}
    return testcache_cfg


def test_run_all_checks_returns_results_and_never_raises(full_cfg):
    results = checks.run_all_checks(full_cfg, enabled_plugins(full_cfg, load_registry(REG)))
    assert results and all(isinstance(r, checks.CheckResult) for r in results)
    names = {r.name for r in results}
    assert {"files", "ini_keys", "plugin:SpliceRegion", "vcf_vs_fasta", "sample", "sex", "env", "prepare_script",
            "disk"} <= names
    assert all(r.level in ("error", "warn", "info") for r in results)
    assert not [r for r in results if r.level == "error" and not r.ok], results


def test_run_all_checks_survives_a_broken_config(full_cfg):
    full_cfg["input"]["vcf"] = "/nonexistent/in.vcf.gz"
    results = checks.run_all_checks(full_cfg, enabled_plugins(full_cfg, load_registry(REG)))
    assert any(r.name == "sample" and not r.ok for r in results)


def _write_cfg(cfg, tmp_path):
    cfg = {k: v for k, v in cfg.items()}
    p = tmp_path / "cfg.yaml"
    p.write_text(yaml.safe_dump(cfg))
    return str(p)


def test_cli_check_requires_input_vcf(testcache_cfg, tmp_path, capsys):
    testcache_cfg["input"]["vcf"] = None
    assert cli.main(["check", "-c", _write_cfg(testcache_cfg, tmp_path)]) == 2
    assert "input.vcf" in capsys.readouterr().err


def test_cli_check_exit_codes(full_cfg, tmp_path, capsys, monkeypatch):
    cfg_path = _write_cfg(full_cfg, tmp_path)
    monkeypatch.setattr(cli, "DEFAULT_REGISTRY", REG)
    code = cli.main(["check", "-c", cfg_path])
    out = capsys.readouterr().out
    assert "plugin:SpliceRegion" in out and "FAIL" not in out
    assert code == 0, out
    full_cfg["sample"] = "NOT_THERE"
    assert cli.main(["check", "-c", _write_cfg(full_cfg, tmp_path)]) == 1
    assert "FAIL" in capsys.readouterr().out


def test_files_check_glob_directory(testcache_cfg, tmp_path):
    reg = tmp_path / "reg_glob.yaml"
    reg.write_text(yaml.safe_dump({"Dir": {"module": "Dir", "category": "manual", "args": "Dir,snv={snv}",
                                           "files": {"snv": {"dest": "Dir/snv", "glob": "chr*.tsv.gz",
                                                             "index": "tbi", "manual": "Get it."}}}}))
    d = tmp_path / "snv"
    testcache_cfg["plugins"] = {"Dir": {"enabled": True, "snv": str(d)}}
    plugins = enabled_plugins(testcache_cfg, load_registry(reg))
    assert "missing directory" in checks.check_files(testcache_cfg, plugins).message
    d.mkdir()
    assert "no files matching" in checks.check_files(testcache_cfg, plugins).message
    for c in ("chr1", "chr2"):
        (d / f"{c}.tsv").write_text(f"chrom\tpos\tref\talt\n{c}\t10\tA\tG\n")
        _bgz_tabix(d / f"{c}.tsv", "-s1", "-b2", "-e2", "-S1")
    Path(f"{d}/chr2.tsv.gz.tbi").unlink()
    r = checks.check_files(testcache_cfg, plugins)
    assert not r.ok and "chr2.tsv.gz" in r.message and "chr1.tsv.gz" not in r.message
    subprocess.run(["tabix", "-f", "-s1", "-b2", "-e2", "-S1", str(d / "chr2.tsv.gz")], check=True)
    assert checks.check_files(testcache_cfg, plugins).ok
