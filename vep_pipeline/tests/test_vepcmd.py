import subprocess
from pathlib import Path

import pytest

from vep_pipeline import checks, vepcmd
from vep_pipeline.registry import enabled_plugins, load_registry

ROOT = Path(__file__).resolve().parents[1]
PARAMS = ROOT / "vep_pipeline" / "vep_params.txt"
CLINVAR_FIELDS = "CLNSIG%CLNREVSTAT%CLNDN%CLNSIGCONF%ALLELEID"


@pytest.fixture
def cfg(minimal_cfg):
    minimal_cfg["plugins"] = {"AlphaMissense": {"enabled": True}, "G2P": {"enabled": True},
                              "NMD": {"enabled": True}}
    minimal_cfg["custom"] = {
        "ClinVar": {"enabled": True, "file": "/d/clinvar.vcf.gz", "format": "vcf", "type": "exact",
                    "fields": ["CLNSIG", "CLNREVSTAT", "CLNDN", "CLNSIGCONF", "ALLELEID"]},
        "ENCODE_cCRE": {"enabled": True, "source": "/d/c.bed", "format": "bed", "type": "overlap"},
    }
    return minimal_cfg


@pytest.fixture
def plugins(cfg):
    return enabled_plugins(cfg, load_registry(ROOT / "plugins.yaml"))


@pytest.fixture
def custom_paths():
    return {"ClinVar": Path("/d/clinvar.vcf.gz"), "ENCODE_cCRE": Path("/d/ccre.bed.gz")}


def test_ini_core_lines(cfg, plugins, custom_paths, tmp_path):
    ini = vepcmd.build_ini(cfg, plugins, custom_paths, "7", tmp_path / "wd")
    lines = ini.splitlines()
    assert "safe 1" in lines and "offline 1" in lines and "cache_version 116" in lines
    assert "assembly GRCh38" in lines and "species homo_sapiens" in lines
    assert f"fasta {cfg['input']['reference_fasta']}" in lines
    assert f"dir_plugins {cfg['tools']['plugins_dir']}" in lines
    assert f"dir_cache {cfg['data']['cache_dir']}" in lines
    assert f"skipped_variants_file {tmp_path / 'wd'}/skipped.txt" in lines
    assert f"warning_file {tmp_path / 'wd'}/warnings.txt" in lines
    assert f"stats_file {tmp_path / 'wd'}/stats.html" in lines
    assert "fork 1" in lines and "buffer_size 100" in lines
    assert not any(l.startswith("database") for l in lines)
    assert "#" not in ini


def test_ini_plugin_and_custom_order(cfg, plugins, custom_paths, tmp_path):
    wd = tmp_path / "wd"
    ini = vepcmd.build_ini(cfg, plugins, custom_paths, "7", wd)
    pl = [l for l in ini.splitlines() if l.startswith("plugin ")]
    assert pl == ["plugin " + p.render("7", wd) for p in plugins]
    assert [p.spec.name for p in plugins] == ["AlphaMissense", "G2P", "NMD"]
    cl = [l for l in ini.splitlines() if l.startswith("custom ")]
    assert [l.split("short_name=")[1].split(",")[0] for l in cl] == ["ClinVar", "ENCODE_cCRE"]
    last_plugin = max(i for i, l in enumerate(ini.splitlines()) if l.startswith("plugin "))
    first_custom = min(i for i, l in enumerate(ini.splitlines()) if l.startswith("custom "))
    assert last_plugin < first_custom


def test_ini_keys_valid(cfg, plugins, custom_paths, tmp_path):
    ini = vepcmd.build_ini(cfg, plugins, custom_paths, "1", tmp_path)
    assert checks.validate_ini_keys(ini, PARAMS) == []


def test_ini_skips_disabled_custom(cfg, plugins, custom_paths, tmp_path):
    cfg["custom"]["ENCODE_cCRE"]["enabled"] = False
    ini = vepcmd.build_ini(cfg, plugins, custom_paths, "1", tmp_path)
    assert "ENCODE_cCRE" not in ini and "ClinVar" in ini


def test_ini_cache_version_from_cfg(cfg, plugins, custom_paths, tmp_path):
    cfg["cache_version"] = 84
    assert "cache_version 84" in vepcmd.build_ini(cfg, plugins, custom_paths, "1", tmp_path).splitlines()


def test_ini_extra_flags(cfg, plugins, custom_paths, tmp_path):
    cfg["run"]["extra_flags"] = ["--regulatory", "--distance=500", "--bogus_flag"]
    ini = vepcmd.build_ini(cfg, plugins, custom_paths, "1", tmp_path)
    lines = ini.splitlines()
    assert "regulatory 1" in lines and "distance 500" in lines and "bogus_flag 1" in lines
    assert checks.validate_ini_keys(ini, PARAMS) == ["bogus_flag"]


def test_g2p_line_per_chunk(cfg, plugins, custom_paths, tmp_path):
    wd = tmp_path / "wd"
    ini = vepcmd.build_ini(cfg, plugins, custom_paths, "7", wd)
    (g2p,) = [l for l in ini.splitlines() if l.startswith("plugin G2P")]
    assert "g2p.7.txt" in g2p and f"{wd}/g2p_log" in g2p
    assert f"{cfg['output']['dir']}/reports/g2p.7.txt" in g2p


def test_custom_arg_clinvar():
    spec = {"enabled": True, "format": "vcf", "type": "exact",
            "fields": ["CLNSIG", "CLNREVSTAT", "CLNDN", "CLNSIGCONF", "ALLELEID"]}
    assert vepcmd.custom_arg("ClinVar", spec, Path("/d/cv.vcf.gz")) == (
        "file=/d/cv.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=" + CLINVAR_FIELDS)


def test_custom_arg_ccre():
    spec = {"enabled": True, "format": "bed", "type": "overlap"}
    assert vepcmd.custom_arg("ENCODE_cCRE", spec, Path("/d/c.bed.gz")) == (
        "file=/d/c.bed.gz,short_name=ENCODE_cCRE,format=bed,type=overlap")


# ----------------------------------------------------------------------------- contig chunks

def _make_vcf(tmp_path, header_contigs, records):
    lines = ["##fileformat=VCFv4.2"] + [f"##contig=<ID={c},length=1000000>" for c in header_contigs]
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    lines += [f"{c}\t{p}\t.\tA\tG\t.\t.\t." for c, p in records]
    plain = tmp_path / "f.vcf"
    plain.write_text("\n".join(lines) + "\n")
    subprocess.run(f"bgzip -f {plain} && tabix -f -p vcf {plain}.gz", shell=True, check=True)
    return Path(f"{plain}.gz")


def test_contig_chunks_primary_and_other(tmp_path):
    vcf = _make_vcf(tmp_path, ["1", "2", "3", "KI270706.1", "GL000194.1"],
                    [("1", 10), ("2", 5), ("2", 9), ("KI270706.1", 4)])
    assert vepcmd.contig_chunks(vcf, "vep116") == [("1", ["1"]), ("2", ["2"]), ("other", ["KI270706.1"])]


def test_contig_chunks_header_order_chr_prefix_and_mt(tmp_path):
    vcf = _make_vcf(tmp_path, ["chr2", "chr1", "chrM", "chrUn_x", "chrX"],
                    [("chr1", 1), ("chr2", 1), ("chrX", 1), ("chrM", 1), ("chrUn_x", 1)])
    assert vepcmd.contig_chunks(vcf, "vep116") == [
        ("chr2", ["chr2"]), ("chr1", ["chr1"]), ("chrM", ["chrM"]), ("chrX", ["chrX"]),
        ("other", ["chrUn_x"])]


def test_contig_chunks_no_other_when_only_primary(tmp_path):
    vcf = _make_vcf(tmp_path, ["MT", "Y"], [("MT", 1)])
    assert vepcmd.contig_chunks(vcf, "vep116") == [("MT", ["MT"])]


def test_contig_chunks_requires_index(tmp_path):
    p = tmp_path / "x.vcf.gz"; p.write_bytes(b"")
    with pytest.raises(Exception):
        vepcmd.contig_chunks(p, "vep116")


# ----------------------------------------------------------------------------- chunk script

def test_chunk_script(cfg, tmp_path):
    wd = tmp_path / "wd"
    s = vepcmd.chunk_script(cfg, "other", wd / "regions.bed", wd / "vep.ini", tmp_path / "in.vcf.gz", wd)
    p = tmp_path / "s.sh"; p.write_text(s)
    assert subprocess.run(["bash", "-n", str(p)]).returncode == 0
    assert s.startswith("#!/usr/bin/env bash\nset -euo pipefail\n")
    assert "export LC_ALL=C PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0" in s
    assert "grep -v '^VEP_DIR_PLUGINS$'" in s
    assert f"cd {wd}" in s and f"mkdir -p {cfg['output']['dir']}/reports" in s
    assert f"rm -rf {wd}/g2p_log" in s
    assert f"bcftools view -R {wd}/regions.bed -Oz -o in.vcf.gz {tmp_path}/in.vcf.gz" in s
    assert f"vep --config {wd}/vep.ini -i in.vcf.gz -o out.vcf.gz" in s
    assert "tabix -f -p vcf out.vcf.gz" in s and "exit 3" in s
    assert "micromamba" not in s


def test_chunk_script_removes_stale_state_before_vep(cfg, tmp_path):
    wd = tmp_path / "wd"
    s = vepcmd.chunk_script(cfg, "1", wd / "regions.bed", wd / "vep.ini", tmp_path / "in.vcf.gz", wd)
    rm = "rm -f out.vcf.gz out.vcf.gz.tbi reconcile.tsv in.vcf.gz in.vcf.gz.tbi warnings.txt skipped.txt\n"
    assert rm in s
    assert s.index(rm) < s.index("bcftools view") < s.index("vep --config")


STALE = ("out.vcf.gz", "out.vcf.gz.tbi", "reconcile.tsv", "in.vcf.gz", "in.vcf.gz.tbi", "warnings.txt",
         "skipped.txt")


def test_chunk_script_failed_rerun_leaves_no_stale_outputs(cfg, tmp_path):
    """I2(b): a rerun that dies before VEP finishes must not leave an older out.vcf.gz and reconcile.tsv
    (unexplained=0) behind, or `run --merge-only` would merge a chunk annotated from an older input."""
    wd = tmp_path / "wd"; wd.mkdir()
    bindir = tmp_path / "bin"; bindir.mkdir()
    (bindir / "bcftools").write_text("#!/usr/bin/env bash\nexit 1\n")
    (bindir / "bcftools").chmod(0o755)
    for name in STALE:
        (wd / name).write_text("old\n")
    script = tmp_path / "s.sh"
    script.write_text(vepcmd.chunk_script(cfg, "1", wd / "regions.bed", wd / "vep.ini", tmp_path / "in.vcf.gz", wd))
    r = subprocess.run(["bash", str(script)], env={"PATH": f"{bindir}:/usr/bin:/bin", "HOME": str(tmp_path)},
                       capture_output=True, text=True)
    assert r.returncode == 1
    assert [n for n in STALE if (wd / n).exists()] == []


def test_chunk_script_per_run_tmpdir_is_created_and_removed(cfg, tmp_path):
    """S5: each run gets its own `mktemp -d` TMPDIR (under $TMPDIR, else the chunk workdir), removed on exit."""
    wd = tmp_path / "wd"; wd.mkdir()
    bindir = tmp_path / "bin"; bindir.mkdir()
    for tool, body in (("bcftools", "exit 0"), ("tabix", "exit 0"),
                       ("vep", 'echo "$TMPDIR" > seen_tmp; [[ -d "$TMPDIR" ]]')):
        (bindir / tool).write_text(f"#!/usr/bin/env bash\n{body}\n")
        (bindir / tool).chmod(0o755)
    script = tmp_path / "s.sh"
    script.write_text(vepcmd.chunk_script(cfg, "1", wd / "regions.bed", wd / "vep.ini", tmp_path / "in.vcf.gz", wd))
    scratch = tmp_path / "scratch"; scratch.mkdir()
    base = {"PATH": f"{bindir}:/usr/bin:/bin", "HOME": str(tmp_path)}
    for env, parent in (({**base, "TMPDIR": str(scratch)}, scratch), (base, wd)):
        r = subprocess.run(["bash", str(script)], env=env, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        seen = Path((wd / "seen_tmp").read_text().strip())
        assert seen.parent == parent and seen.name.startswith("vep_tmp.") and not seen.exists()
    assert list(scratch.iterdir()) == []


def test_chunk_script_stale_warnings_do_not_fail_clean_rerun(cfg, tmp_path):
    """VEP only truncates warnings.txt when it emits a warning, so the script must clear an old one."""
    wd = tmp_path / "wd"; wd.mkdir()
    bindir = tmp_path / "bin"; bindir.mkdir()
    for tool in ("bcftools", "tabix", "vep"):      # stubs: succeed, write no warnings
        (bindir / tool).write_text("#!/usr/bin/env bash\nexit 0\n")
        (bindir / tool).chmod(0o755)
    (wd / "warnings.txt").write_text("Failed to open plugin file\n")
    (wd / "skipped.txt").write_text("old\n")
    script = tmp_path / "s.sh"
    script.write_text(vepcmd.chunk_script(cfg, "1", wd / "regions.bed", wd / "vep.ini", tmp_path / "in.vcf.gz", wd))
    env = {"PATH": f"{bindir}:/usr/bin:/bin", "HOME": str(tmp_path)}
    r = subprocess.run(["bash", str(script)], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert not (wd / "warnings.txt").exists() and not (wd / "skipped.txt").exists()
    # and a warning emitted by the run itself still fails the chunk
    (bindir / "vep").write_text("#!/usr/bin/env bash\necho 'went wrong' > warnings.txt\n")
    r = subprocess.run(["bash", str(script)], env=env, capture_output=True, text=True)
    assert r.returncode == 3
