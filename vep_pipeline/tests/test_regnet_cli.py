import yaml

from test_regnet_build import make_scdori, make_sldsc, tsv
from vep_pipeline import cli


def test_regnet_build_command_writes_the_resource(minimal_cfg, tmp_path, capsys):
    scdori = tmp_path / "scdori"
    make_scdori(scdori)
    a = make_sldsc(tmp_path, "srcA")
    minimal_cfg["regnet"] = {"scdori_dir": str(scdori), "sldsc_dirs": [str(a)], "sldsc_model": "M1"}
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(yaml.safe_dump(minimal_cfg))
    assert cli.main(["regnet-build", "-c", str(cfg)]) == 0
    out = tmp_path / "out" / "regnet"
    assert {r["trait"] for r in tsv(out / "traits.tsv")} == {"T1", "T2"}
    assert f"regnet resource in {out}: " in capsys.readouterr().out


def test_regnet_build_without_inputs_is_a_config_error(minimal_cfg, tmp_path, capsys):
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(yaml.safe_dump(minimal_cfg))
    assert cli.main(["regnet-build", "-c", str(cfg)]) == 2
    assert "regnet.scdori_dir" in capsys.readouterr().err


def test_prioritise_slurm_checks_regnet_before_rendering(minimal_cfg, tmp_path, capsys, monkeypatch):
    rendered = []
    monkeypatch.setattr(cli.pr, "render_slurm", lambda *a: rendered.append(a))
    minimal_cfg["regnet"] = {"enabled": True, "build_dir": str(tmp_path / "missing"), "traits": ["T1"]}
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(yaml.safe_dump(minimal_cfg))
    assert cli.main(["prioritise", "-c", str(cfg), "--slurm"]) == 2
    assert "regnet-build" in capsys.readouterr().err and rendered == []
