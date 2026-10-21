"""End-to-end run of the real VEP against the cache that ships with ensembl-vep (no downloads, ~minutes)."""
import shutil
import subprocess
from pathlib import Path

import pytest

from vep_pipeline import jobs
from vep_pipeline.registry import enabled_plugins, load_registry

from conftest import TESTCACHE, TESTCACHE_FASTA, TESTDATA

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.integration


def _sh(*cmd):
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def _sample_vcf(tmp_path) -> Path:
    """test.vcf with its sample renamed S1 and every genotype 0/1, bgzipped and indexed."""
    out = tmp_path / "in" / "snv.vcf"
    out.parent.mkdir()
    lines = []
    for line in (TESTDATA / "input" / "test.vcf").read_text().splitlines():
        f = line.split("\t")
        if line.startswith("#CHROM"):
            lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
            f[9] = "S1"
        elif not line.startswith("#"):
            f[9] = "0/1"
        lines.append("\t".join(f))
        if len(f) > 1 and f[1] == "25592911":                     # a record VEP skips: a <DEL> with neither END nor SVLEN
            lines.append("\t".join(["21", "25593000", "skipme", "C", "<DEL>", ".", ".", "SVTYPE=DEL", "GT", "0/1"]))
    out.write_text("\n".join(lines) + "\n")
    _sh("bgzip", "-f", str(out))
    _sh("tabix", "-f", "-p", "vcf", f"{out}.gz")
    return Path(f"{out}.gz")


@pytest.fixture
def run_cfg(testcache_cfg, tmp_path):
    cfg = testcache_cfg
    fasta = tmp_path / "ref" / "test.fa"          # VEP may write a .fai next to the FASTA: keep it out of the repo
    fasta.parent.mkdir()
    shutil.copy(TESTCACHE_FASTA, fasta)
    shutil.copy(f"{TESTCACHE_FASTA}.fai", f"{fasta}.fai")
    cfg["input"]["reference_fasta"] = str(fasta)
    cfg["input"]["vcf"] = str(_sample_vcf(tmp_path))
    cfg["input"]["normalize"] = False
    cfg["plugins"] = {"NMD": {"enabled": True}, "SpliceRegion": {"enabled": True}}
    cfg["custom"] = {}
    cfg["run"] = {"fork": 1, "buffer_size": 50, "extra_flags": []}
    cfg["tools"]["plugins_dir"] = str(Path(cfg["tools"]["plugins_dir"]))
    cfg["output"]["dir"] = str(tmp_path / "out")
    return cfg


def _csq_column(vcf):
    r = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\t%INFO/CSQ\n", str(vcf)],
                       capture_output=True, text=True, check=True)
    return r.stdout.splitlines()


def _run(cfg):
    plugins = enabled_plugins(cfg, load_registry(ROOT / "plugins.yaml"))
    return jobs.run_local(jobs.render_run(cfg, plugins), parallel=2)


def test_snv_run_end_to_end_and_deterministic(run_cfg, tmp_path):
    out_dir = Path(run_cfg["output"]["dir"])
    final = _run(run_cfg)
    assert final == out_dir / "S1.vep.vcf.gz" and final.is_file() and Path(f"{final}.tbi").is_file()
    assert (out_dir / "fields.tsv").is_file()

    hdr = subprocess.run(["bcftools", "view", "-h", str(final)], capture_output=True, text=True).stdout
    fmt = next(l for l in hdr.splitlines() if "ID=CSQ" in l).split("Format: ")[1].rstrip('">').split("|")
    assert "NMD" in fmt and "SpliceRegion" in fmt

    rec = (out_dir / "work" / "21" / "reconcile.tsv").read_text().splitlines()
    cols = dict(zip(rec[0].split("\t"), map(int, rec[1].split("\t"))))
    print("reconcile.tsv chunk 21:", cols)
    assert cols["unexplained"] == 0 and cols["out"] > 0 and cols["in"] == cols["out"] + cols["skipped"]
    assert cols["skipped"] == 1, "the incomplete <DEL> must be reported by VEP and accounted for by line number"

    first = _csq_column(final)
    assert any("missense_variant" in line for line in first)
    shutil.rmtree(out_dir)
    assert _csq_column(_run(run_cfg)) == first


def test_slurm_scripts_run_end_to_end(run_cfg, tmp_path):
    """Execute the generated array task and merge job by hand (no sbatch): they must produce the same final VCF."""
    import os
    import yaml
    from vep_pipeline import cli
    cfg_file = tmp_path / "cfg.yaml"
    cfg_file.write_text(yaml.safe_dump(run_cfg))
    assert cli.main(["run", "-c", str(cfg_file), "--slurm"]) == 0
    work = Path(run_cfg["output"]["dir"]) / "work"
    env = {**os.environ, "SLURM_ARRAY_TASK_ID": "1", "SLURM_JOB_ID": "77", "SLURM_TMPDIR": str(tmp_path / "scratch")}
    (tmp_path / "scratch").mkdir()
    r = subprocess.run(["bash", str(work / "array.sbatch")], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert (work / "21" / "reconcile.tsv").read_text().splitlines()[1].split("\t")[-1] == "0"
    assert not (tmp_path / "scratch" / "vep.77.1").exists()          # per-task TMPDIR is cleaned up
    r = subprocess.run(["bash", str(work / "merge.sbatch")], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    final = Path(run_cfg["output"]["dir"]) / "S1.vep.vcf.gz"
    assert final.is_file() and Path(f"{final}.tbi").is_file()
    assert cli.main(["reconcile", "-c", str(cfg_file), "--chunk", "21"]) == 0
    assert cli.main(["reconcile", "-c", str(cfg_file), "--chunk", "nope"]) == 2


def test_sv_run(run_cfg, tmp_path):
    sv = tmp_path / "in" / "sv.vcf"
    sv.write_text("\n".join([
        "##fileformat=VCFv4.2",
        "##contig=<ID=21,length=46709983>",
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
        '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="SV length">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
        "21\t25585000\tsv1\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=25600000;SVLEN=-15000\tGT\t0/1",
        "21\t25640000\tsv2\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=25700000;SVLEN=-60000\tGT\t0/1",
    ]) + "\n")
    _sh("bgzip", "-f", str(sv))
    _sh("tabix", "-f", "-p", "vcf", f"{sv}.gz")
    run_cfg["sv"] = {"enabled": True, "vcf": f"{sv}.gz"}
    script = jobs.render_sv_run(run_cfg)
    assert script is not None and script.is_file()
    out = jobs.run_sv(run_cfg)
    out_dir = Path(run_cfg["output"]["dir"])
    assert out == out_dir / "S1.sv.vep.vcf.gz" and Path(f"{out}.tbi").is_file()
    n = subprocess.run(["bcftools", "view", "-H", str(out)], capture_output=True, text=True).stdout.splitlines()
    assert len(n) == 2
    assert "everything 1" in (out_dir / "work" / "sv" / "vep.ini").read_text().splitlines()   # regulatory stays on
    assert not any(l.startswith("plugin ") for l in (out_dir / "work" / "sv" / "vep.ini").read_text().splitlines())
    hdr = subprocess.run(["bcftools", "view", "-h", str(out)], capture_output=True, text=True).stdout
    assert "ID=CSQ" in hdr
    assert all(c for c in _csq_column(out))                       # every SV record carries a CSQ
    rec = (out_dir / "work" / "sv" / "reconcile.tsv").read_text().splitlines()
    assert dict(zip(rec[0].split("\t"), rec[1].split("\t"))) == {
        "in": "2", "out": "2", "skipped": "0", "no_alt_gt": "0", "unexplained": "0"}
