import os, re, subprocess
from pathlib import Path

import pytest, yaml
from vep_pipeline.config import load_config, resolve_placeholders, ConfigError

def write(tmp_path, data):
    p = tmp_path / "config.yaml"; p.write_text(yaml.safe_dump(data)); return p

def test_placeholders_resolve_nested():
    out = resolve_placeholders({"data": {"root": "/d", "cache_dir": "{data.root}/cache"}, "x": "{data.cache_dir}/v"})
    assert out["data"]["cache_dir"] == "/d/cache" and out["x"] == "/d/cache/v"

def test_unknown_placeholder_raises():
    with pytest.raises(ConfigError, match="unknown placeholder"):
        resolve_placeholders({"a": "{nope.key}"})

def test_cycle_raises():
    with pytest.raises(ConfigError, match="cycle"):
        resolve_placeholders({"a": "{b}", "b": "{a}"})

def test_plugin_bool_normalised(minimal_cfg, tmp_path):
    minimal_cfg["plugins"] = {"NMD": True, "REVEL": False}
    cfg = load_config(write(tmp_path, minimal_cfg))
    assert cfg["plugins"] == {"NMD": {"enabled": True}, "REVEL": {"enabled": False}}

@pytest.mark.parametrize("key", ["sample", "input.vcf", "input.reference_fasta", "data.root"])
def test_required_keys(minimal_cfg, tmp_path, key):
    node = minimal_cfg; *parents, leaf = key.split(".")
    for p in parents: node = node[p]
    node[leaf] = None
    with pytest.raises(ConfigError, match=key):
        load_config(write(tmp_path, minimal_cfg))

def test_assembly_must_be_grch38(minimal_cfg, tmp_path):
    minimal_cfg["assembly"] = "GRCh37"
    with pytest.raises(ConfigError, match="GRCh38"):
        load_config(write(tmp_path, minimal_cfg))

def test_whitespace_in_path_rejected(minimal_cfg, tmp_path):
    minimal_cfg["output"]["dir"] = str(tmp_path / "has space")
    with pytest.raises(ConfigError, match="whitespace"):
        load_config(write(tmp_path, minimal_cfg))

WRAPPER = Path(__file__).resolve().parents[1] / "bin" / "vep-pipeline"


def test_wrapper_resolves_relative_config_from_the_callers_cwd(minimal_cfg, tmp_path):
    """I4: bin/vep-pipeline used to cd into the package dir, so `-c sub/my.yaml` resolved there, not in the caller's
    directory (and a bare `-c config.yaml` silently loaded the shipped config)."""
    minimal_cfg["plugins"] = {"OnlyInTmpConfig": True}
    cwd = tmp_path / "elsewhere"
    (cwd / "sub").mkdir(parents=True)
    (cwd / "sub" / "my.yaml").write_text(yaml.safe_dump(minimal_cfg))
    r = subprocess.run(["bash", str(WRAPPER), "list-plugins", "-c", "sub/my.yaml"], cwd=cwd,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert re.search(r"^OnlyInTmpConfig\s.*unknown plugin$", r.stdout, re.M)


def test_wrapper_default_config_is_the_shipped_one_from_any_cwd(tmp_path):
    from vep_pipeline import cli
    shipped = WRAPPER.parents[1] / "config.yaml"
    assert cli.DEFAULT_CONFIG.is_absolute() and cli.DEFAULT_CONFIG == shipped
    r = subprocess.run(["bash", str(WRAPPER), "list-plugins", "--help"], cwd=tmp_path, capture_output=True, text=True,
                       env={**os.environ, "COLUMNS": "400"})
    assert r.returncode == 0, r.stderr
    assert f"(default: {shipped})" in r.stdout


def test_shipped_config_parses_when_vcf_set(tmp_path):
    from pathlib import Path
    raw = yaml.safe_load((Path(__file__).parents[1] / "config.yaml").read_text())
    raw["input"]["vcf"] = "/x.vcf.gz"
    assert load_config(write(tmp_path, raw))["cache_version"] == 116
