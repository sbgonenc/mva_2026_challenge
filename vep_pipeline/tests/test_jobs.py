import json
import subprocess
from pathlib import Path

import pytest
import yaml

from vep_pipeline import cli, jobs
from vep_pipeline.registry import enabled_plugins, load_registry

ROOT = Path(__file__).resolve().parents[1]


def _sh(*cmd):
    subprocess.run(cmd, check=True, capture_output=True, text=True)


@pytest.fixture
def registry():
    return load_registry(ROOT / "plugins.yaml")


@pytest.fixture
def cfg(minimal_cfg, vcf_factory):
    minimal_cfg["input"]["vcf"] = str(vcf_factory(
        {"1": 1000, "2": 1000, "3": 1000, "KI270": 500},
        [("1", 10, "A", "C"), ("2", 20, "A", "C"), ("3", 30, "A", "C"), ("KI270", 5, "A", "C")],
        samples=("S1", "S2")))
    minimal_cfg["plugins"] = {"NMD": {"enabled": True}}
    return minimal_cfg


@pytest.fixture
def plan(cfg, registry):
    return jobs.render_run(cfg, enabled_plugins(cfg, registry))


def _vcf(path, records, samples=("S1",), header_extra=()):
    lines = ["##fileformat=VCFv4.2", "##contig=<ID=21,length=1000>",
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', *header_extra,
             "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples])]
    for pos, ref, alt, *gt in records:
        lines.append("\t".join(["21", str(pos), ".", ref, alt, ".", ".", ".", "GT", *(gt or ["0/1"] * len(samples))]))
    Path(path).write_text("\n".join(lines) + "\n")
    return Path(path)


# ----------------------------------------------------------------------------- prepare_input

def test_prepare_input_subsets_to_sample(cfg):
    out = jobs.prepare_input(cfg)
    assert out == Path(cfg["output"]["dir"]) / "S1.input.vcf.gz"
    assert Path(f"{out}.tbi").exists()
    names = subprocess.run(["bcftools", "query", "-l", str(out)], capture_output=True, text=True).stdout.split()
    assert names == ["S1"]


def test_prepare_input_normalize_calls_prepare_script(cfg, monkeypatch):
    cfg["input"]["normalize"] = True
    cfg["run"]["fork"] = 3
    calls = []

    def fake_run(cmd, **kw):
        calls.append(list(map(str, cmd)))
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(jobs.subprocess, "run", fake_run)
    out = jobs.prepare_input(cfg)
    out_dir = Path(cfg["output"]["dir"])
    norm = f"{out_dir}/S1.norm.vcf.gz"
    first = calls[0]
    assert first[:5] == ["micromamba", "run", "-n", "exomiser", "bash"]
    assert first[5] == cfg["tools"]["prepare_vcf"]
    assert first[first.index("-i") + 1] == cfg["input"]["vcf"]
    assert first[first.index("-r") + 1] == cfg["input"]["reference_fasta"]
    assert first[first.index("-o") + 1] == norm
    assert first[first.index("-s") + 1] == "S1" and first[first.index("-c") + 1] == "3"
    assert "-q" not in first
    subset = calls[1]
    assert subset[:4] == ["micromamba", "run", "-n", "vep116"] and "bcftools" in subset
    assert norm in subset and f"{out_dir}/S1.input.vcf.gz" in subset
    assert any("tabix" in c for c in calls[2])
    assert out == out_dir / "S1.input.vcf.gz"


@pytest.fixture
def norm_cfg(cfg, monkeypatch, tmp_path, fasta_factory):
    """cfg with input.normalize on; `prepare-patient-vcf.sh` faked (writes the .norm.vcf.gz and its .tbi, records the
    call; fails while `fail` is non-empty), subsetting stubbed out. Returns (cfg, calls, fail)."""
    cfg["input"]["normalize"] = True
    cfg["input"]["reference_fasta"] = str(fasta_factory({"1": 100}))
    script = tmp_path / "prepare.sh"
    script.write_text("#!/usr/bin/env bash\necho v1\n")
    cfg["tools"]["prepare_vcf"] = str(script)
    calls, fail = [], []

    def fake_run(cmd, env_name, **kw):
        cmd = list(map(str, cmd))
        calls.append(cmd)
        out = Path(cmd[cmd.index("-o") + 1])
        out.write_bytes(b"partial")
        if fail:
            return subprocess.CompletedProcess(cmd, 1, "", "died")
        Path(f"{out}.tbi").write_bytes(b"tbi")
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(jobs, "_run", fake_run)
    monkeypatch.setattr(jobs, "_subset", lambda cfg, src, dest: dest)
    return cfg, calls, fail


def test_normalised_vcf_reused_with_matching_fingerprint(norm_cfg):
    cfg, calls, _fail = norm_cfg
    jobs.prepare_input(cfg)
    jobs.prepare_input(cfg)
    assert len(calls) == 1
    norm = Path(cfg["output"]["dir"]) / "S1.norm.vcf.gz"
    fp = json.loads(Path(f"{norm}.fingerprint.json").read_text())
    assert fp["source"] == str(Path(cfg["input"]["vcf"]).resolve()) and fp["sample"] == "S1"
    assert set(fp) == {"source", "source_size", "source_mtime_ns", "fasta", "fasta_size", "prepare_vcf",
                       "prepare_vcf_sha256", "sample"}


def _older_copy(cfg, tmp_path):
    """Point input.vcf at a different file whose mtime is older than the existing normalised VCF."""
    import os, shutil
    other = tmp_path / "other_input.vcf.gz"
    shutil.copyfile(cfg["input"]["vcf"], other)
    os.utime(other, (1_000_000_000, 1_000_000_000))
    cfg["input"]["vcf"] = str(other)


@pytest.mark.parametrize("change", ["older_input_path", "fasta", "prepare_script", "sample", "no_fingerprint"])
def test_normalised_vcf_redone_when_fingerprint_differs(norm_cfg, tmp_path, fasta_factory, change):
    """I3: `_fresh(norm, src)` reused {prefix}.norm.vcf.gz whenever it was newer than input.vcf, ignoring an input path
    change to an older file, a different FASTA, a changed prepare script and a different sample."""
    cfg, calls, _fail = norm_cfg
    jobs.prepare_input(cfg)
    norm = Path(cfg["output"]["dir"]) / "S1.norm.vcf.gz"
    if change == "older_input_path":
        _older_copy(cfg, tmp_path)
    elif change == "fasta":
        cfg["input"]["reference_fasta"] = str(fasta_factory({"1": 200}))
    elif change == "prepare_script":
        Path(cfg["tools"]["prepare_vcf"]).write_text("#!/usr/bin/env bash\necho v2\n")
    elif change == "sample":
        cfg["sample"] = "S2"
    else:
        Path(f"{norm}.fingerprint.json").unlink()            # a normalised VCF from before fingerprints existed
    jobs.prepare_input(cfg)
    assert len(calls) == 2


def test_failed_normalisation_leaves_no_reusable_fingerprint(norm_cfg):
    cfg, calls, fail = norm_cfg
    jobs.prepare_input(cfg)
    cfg["sample"] = "S2"
    fail.append(True)
    with pytest.raises(jobs.JobError, match="prepare-patient-vcf.sh failed"):
        jobs.prepare_input(cfg)
    norm = Path(cfg["output"]["dir"]) / "S1.norm.vcf.gz"
    assert not Path(f"{norm}.fingerprint.json").exists()
    fail.clear()
    cfg["sample"] = "S1"                                       # back to the first fingerprint: still redone
    jobs.prepare_input(cfg)
    assert len(calls) == 3


# ----------------------------------------------------------------------------- render_run

def test_render_run_writes_chunks(plan, cfg):
    wd = Path(cfg["output"]["dir"]) / "work"
    assert plan.workdir == wd
    assert [c.name for c in plan.chunks] == ["1", "2", "3", "other"]
    assert (wd / "contigs.txt").read_text().split() == ["1", "2", "3", "other"]
    for c in plan.chunks:
        assert c.workdir == wd / c.name
        assert (c.workdir / "vep.ini").is_file() and c.script == c.workdir / "run.sh"
        assert f"bash" in c.script.read_text() and "vep --config" in c.script.read_text()
        assert "plugin NMD" in (c.workdir / "vep.ini").read_text()
    assert plan.final_vcf == Path(cfg["output"]["dir"]) / "S1.vep.vcf.gz"
    assert plan.merge_script.is_file() and "--merge-only" in plan.merge_script.read_text()
    assert not plan.individual_mode and plan.sv is None
    # the input handed to the chunks has one sample only
    assert plan.input_vcf == Path(cfg["output"]["dir"]) / "S1.input.vcf.gz"


def test_render_run_adds_sv_job_after_chunks(cfg, registry, vcf_factory, tmp_path):
    cfg["sv"] = {"enabled": True, "vcf": str(vcf_factory({"1": 1000}, [("1", 10, "N", "<DEL>")], samples=("S1",)))}
    plan = jobs.render_run(cfg, enabled_plugins(cfg, registry))
    assert plan.sv is not None and plan.sv.kind == "sv" and plan.sv_final == Path(cfg["output"]["dir"]) / "S1.sv.vep.vcf.gz"
    assert [c.name for c in plan.chunks] == ["1", "2", "3", "other"]          # sv is not merged into the SNV VCF
    assert (plan.workdir / "contigs.txt").read_text().split() == ["1", "2", "3", "other", "sv"]
    array, _ = jobs.render_slurm(plan, cfg)
    assert "--array=1-5" in array.read_text()
    again = jobs.plan_from_disk(cfg, enabled_plugins(cfg, registry))
    assert [j.name for j in again.jobs] == ["1", "2", "3", "other", "sv"]


def test_render_run_g2p_forces_individual(cfg, registry):
    cfg["plugins"] = {"G2P": {"enabled": True}}
    assert jobs.render_run(cfg, enabled_plugins(cfg, registry)).individual_mode


def test_render_sv_run_ini_and_script(cfg, registry, vcf_factory, tmp_path):
    from vep_pipeline import checks
    assert jobs.render_sv_run(cfg) is None
    cfg["sv"] = {"enabled": True, "vcf": str(vcf_factory({"1": 1000}, [("1", 10, "N", "<DEL>")], samples=("S1", "S2")))}
    cfg["custom"] = {"Trk": {"enabled": True, "file": "/x.vcf.gz", "format": "vcf"}}
    script = jobs.render_sv_run(cfg)
    ini = (script.parent / "vep.ini").read_text()
    lines = ini.splitlines()
    assert script.parent == Path(cfg["output"]["dir"]) / "work" / "sv"
    assert "max_sv_size 100000000" in lines and "everything 1" in lines
    assert not any(l.startswith("plugin ") for l in lines)
    assert not any(l.startswith("custom ") for l in lines) and "safe 1" in lines and "offline 1" in lines
    assert checks.validate_ini_keys(ini, ROOT / "vep_pipeline" / "vep_params.txt") == []
    names = subprocess.run(["bcftools", "query", "-l", str(Path(cfg["output"]["dir"]) / "S1.sv.input.vcf.gz")],
                           capture_output=True, text=True).stdout.split()
    assert names == ["S1"]
    cfg["sv"]["vcf"] = None
    with pytest.raises(Exception, match="sv.vcf"):
        jobs.render_sv_run(cfg)


# ----------------------------------------------------------------------------- slurm

def test_render_slurm_three_chunks(cfg, registry, vcf_factory):
    cfg["input"]["vcf"] = str(vcf_factory({"1": 100, "2": 100, "3": 100},
                                          [("1", 1, "A", "C"), ("2", 1, "A", "C"), ("3", 1, "A", "C")]))
    cfg["slurm"].update(partition="short", cpus=6, mem="12G", time="02:00:00")
    plan = jobs.render_run(cfg, enabled_plugins(cfg, registry))
    array, merge = jobs.render_slurm(plan, cfg)
    text = array.read_text()
    assert "#SBATCH --array=1-3" in text and "TMPDIR=${SLURM_TMPDIR" in text
    assert "#SBATCH --partition=short" in text and "#SBATCH --cpus-per-task=6" in text
    assert "#SBATCH --mem=12G" in text and "#SBATCH --time=02:00:00" in text
    assert "#SBATCH --export=ALL" in text
    log_dir = cfg["slurm"]["log_dir"]
    assert f"#SBATCH --output={log_dir}/vep_%A_%a.out" in text and Path(log_dir).is_dir()
    assert "SLURM_ARRAY_TASK_ID" in text and "contigs.txt" in text
    assert "micromamba run -n vep116 bash" in text and "run.sh" in text
    assert "vep_pipeline.cli reconcile --chunk" in text
    assert "#SBATCH" in merge.read_text() and str(plan.merge_script) in merge.read_text()
    assert subprocess.run(["bash", "-n", str(array)]).returncode == 0
    assert subprocess.run(["bash", "-n", str(merge)]).returncode == 0


def test_submit_slurm_chains_merge_with_dependency(tmp_path, monkeypatch):
    calls = []

    def fake_run(cmd, **kw):
        calls.append(list(cmd))
        return subprocess.CompletedProcess(cmd, 0, "123\n" if len(calls) == 1 else "456\n", "")

    monkeypatch.setattr(jobs.subprocess, "run", fake_run)
    a, m = tmp_path / "a.sbatch", tmp_path / "m.sbatch"
    assert jobs.submit_slurm(a, m) == ("123", "456")
    assert calls[0][0] == "sbatch" and str(a) in calls[0] and not any("dependency" in c for c in calls[0])
    assert "--dependency=afterok:123" in calls[1] and str(m) in calls[1]


def test_submit_slurm_failure_is_loud(tmp_path, monkeypatch):
    monkeypatch.setattr(jobs.subprocess, "run",
                        lambda cmd, **kw: subprocess.CompletedProcess(cmd, 1, "", "invalid partition"))
    with pytest.raises(jobs.JobError, match="invalid partition"):
        jobs.submit_slurm(tmp_path / "a", tmp_path / "m")


# ----------------------------------------------------------------------------- reconcile

SKIPPED = "[VEP skipped the following variants from in.vcf.gz]\nLine {n:<5}\t21 {pos} . A T\tdeletion looks incomplete\n"


def test_reconcile_skipped_and_hom_ref_ok(tmp_path):
    inp = _vcf(tmp_path / "in.vcf", [(100, "A", "T", "0/1"), (200, "A", "T", "0/1"), (300, "A", "T", "0/0")])
    out = _vcf(tmp_path / "out.vcf", [(100, "A", "T", "0/1")])
    skipped = tmp_path / "skipped.txt"
    skipped.write_text(SKIPPED.format(n=2, pos=200))
    rep = jobs.reconcile_chunk(inp, out, skipped, "S1", individual_mode=True)
    assert (rep.n_in, rep.n_out, rep.skipped, rep.no_alt_gt, rep.unexplained) == (3, 1, 1, 1, 0)
    rows = (tmp_path / "reconcile.tsv").read_text().splitlines()
    assert rows[0].split("\t") == ["in", "out", "skipped", "no_alt_gt", "unexplained"]
    assert rows[1].split("\t") == ["3", "1", "1", "1", "0"]


def test_reconcile_unexplained_key_raises(tmp_path):
    inp = _vcf(tmp_path / "in.vcf", [(100, "A", "T"), (200, "A", "T"), (300, "A", "T")])
    out = _vcf(tmp_path / "out.vcf", [(100, "A", "T")])
    skipped = tmp_path / "skipped.txt"
    skipped.write_text(SKIPPED.format(n=2, pos=200))
    with pytest.raises(jobs.ReconcileError, match="21:300:A>T"):
        jobs.reconcile_chunk(inp, out, skipped, "S1", individual_mode=True)
    assert (tmp_path / "reconcile.tsv").read_text().splitlines()[1].split("\t")[-1] == "1"


def test_reconcile_hom_ref_only_explained_in_individual_mode(tmp_path):
    inp = _vcf(tmp_path / "in.vcf", [(100, "A", "T", "0/1"), (300, "A", "T", "0/0")])
    out = _vcf(tmp_path / "out.vcf", [(100, "A", "T", "0/1")])
    with pytest.raises(jobs.ReconcileError, match="21:300:A>T"):
        jobs.reconcile_chunk(inp, out, tmp_path / "skipped.txt", "S1", individual_mode=False)


def test_reconcile_output_must_be_subset(tmp_path):
    inp = _vcf(tmp_path / "in.vcf", [(100, "A", "T")])
    out = _vcf(tmp_path / "out.vcf", [(100, "A", "T"), (150, "G", "C")])
    with pytest.raises(jobs.ReconcileError, match="21:150:G>C"):
        jobs.reconcile_chunk(inp, out, tmp_path / "skipped.txt", "S1", individual_mode=False)


def test_reconcile_lists_at_most_20_examples(tmp_path):
    inp = _vcf(tmp_path / "in.vcf", [(i, "A", "T") for i in range(1, 41)])
    out = _vcf(tmp_path / "out.vcf", [])
    with pytest.raises(jobs.ReconcileError) as e:
        jobs.reconcile_chunk(inp, out, tmp_path / "skipped.txt", "S1", individual_mode=False)
    assert "21:20:A>T" in str(e.value) and "21:21:A>T" not in str(e.value) and "40" in str(e.value)


def test_reconcile_sv_counts_records(tmp_path):
    inp = _vcf(tmp_path / "in.vcf", [(100, "A", "<DEL>"), (200, "A", "<DEL>")])
    jobs.reconcile_sv(inp, _vcf(tmp_path / "out.vcf", [(100, "A", "<DEL>"), (200, "A", "<DEL>")]),
                      tmp_path / "skipped.txt")
    assert (tmp_path / "reconcile.tsv").read_text().splitlines()[1].split("\t")[-1] == "0"
    with pytest.raises(jobs.ReconcileError):
        jobs.reconcile_sv(inp, _vcf(tmp_path / "o2.vcf", [(100, "A", "<DEL>")]), tmp_path / "skipped.txt")


# ----------------------------------------------------------------------------- CSQ headers and merge

def _csq(path, fields, records=()):
    hdr = ('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
           f'Format: {"|".join(fields)}">')
    return _vcf(path, records, header_extra=[hdr])


def test_csq_headers_identical_ok(tmp_path):
    a = _csq(tmp_path / "a.vcf", ["Allele", "Consequence", "NMD"])
    b = _csq(tmp_path / "b.vcf", ["Allele", "Consequence", "NMD"])
    jobs.assert_csq_headers_identical([a, b], ["NMD", "Consequence"])


def test_csq_header_missing_field_raises(tmp_path):
    a = _csq(tmp_path / "a.vcf", ["Allele", "Consequence", "NMD"])
    b = _csq(tmp_path / "b.vcf", ["Allele", "Consequence"])
    with pytest.raises(jobs.JobError, match="differ"):
        jobs.assert_csq_headers_identical([a, b], [])
    with pytest.raises(jobs.JobError, match="NMD"):
        jobs.assert_csq_headers_identical([b], ["NMD"])


def test_csq_header_missing_raises(tmp_path):
    with pytest.raises(jobs.JobError, match="CSQ"):
        jobs.assert_csq_headers_identical([_vcf(tmp_path / "a.vcf", [])], [])


def test_csq_prefix_pattern_for_custom_tracks(tmp_path):
    a = _csq(tmp_path / "a.vcf", ["Allele", "ClinVar_CLNSIG"])
    jobs.assert_csq_headers_identical([a], ["ClinVar*"])
    with pytest.raises(jobs.JobError, match="Other"):
        jobs.assert_csq_headers_identical([a], ["Other*"])


def test_merge_chunks_keeps_chunk_order(tmp_path, vcf_factory):
    a = vcf_factory({"1": 100, "2": 100}, [("1", 5, "A", "C")])
    b = vcf_factory({"1": 100, "2": 100}, [("2", 7, "G", "T")])
    out = jobs.merge_chunks([a, b], tmp_path / "m" / "final.vcf.gz")
    recs = subprocess.run(["bcftools", "query", "-f", "%CHROM:%POS\n", str(out)], capture_output=True, text=True).stdout
    assert recs.split() == ["1:5", "2:7"] and Path(f"{out}.tbi").exists()


# ----------------------------------------------------------------------------- chunk inputs

STUB_VEP = """#!/usr/bin/env bash
while [[ $# -gt 0 ]]; do case "$1" in -i) i=$2; shift;; -o) o=$2; shift;; esac; shift; done
cp "$i" "$o"
"""


def _run_chunk_with_stub_vep(job, tmp_path):
    """Run a rendered run.sh with the real bcftools/tabix and a `vep` that copies its input to its output."""
    import os
    bindir = tmp_path / "stubbin"
    bindir.mkdir(exist_ok=True)
    (bindir / "vep").write_text(STUB_VEP)
    (bindir / "vep").chmod(0o755)
    r = subprocess.run(["bash", str(job.script)], capture_output=True, text=True,
                       env={**os.environ, "PATH": f"{bindir}:{os.environ['PATH']}"})
    assert r.returncode == 0, r.stderr
    return r


def test_chunk_regions_survive_contig_names_with_colon_and_star(minimal_cfg, registry, vcf_factory, tmp_path):
    """I2(c): `bcftools view -r <contigs>` parses `HLA-A*01:01:01:01` as a region string (and the shell sees a glob);
    the chunk must select whole contigs from a BED file instead."""
    hla = "HLA-A*01:01:01:01"
    minimal_cfg["input"]["vcf"] = str(vcf_factory({"1": 1000, hla: 3000, "KI270": 500},
                                                  [("1", 10, "A", "C"), (hla, 2999, "A", "C"), (hla, 5, "A", "G"),
                                                   ("KI270", 5, "A", "C")]))
    minimal_cfg["plugins"] = {}
    plan = jobs.render_run(minimal_cfg, enabled_plugins(minimal_cfg, registry))
    other = next(j for j in plan.chunks if j.name == "other")
    assert (other.workdir / "regions.bed").read_text() == f"{hla}\t0\t3000\nKI270\t0\t500\n"
    _run_chunk_with_stub_vep(other, tmp_path)
    got = subprocess.run(["bcftools", "query", "-f", "%CHROM:%POS\n", str(other.in_vcf)],
                         capture_output=True, text=True, check=True).stdout.split()
    assert got == [f"{hla}:5", f"{hla}:2999", "KI270:5"]


def test_chunk_regions_length_falls_back_to_fasta_fai(minimal_cfg, registry, tmp_path, fasta_factory):
    """A contig without a length in the VCF header takes it from the reference .fai; with neither, render fails."""
    plain = tmp_path / "nolen.vcf"
    plain.write_text("##fileformat=VCFv4.2\n##contig=<ID=1>\n"
                     '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1\n")
    _sh("bgzip", "-f", str(plain))
    _sh("tabix", "-f", "-p", "vcf", f"{plain}.gz")
    minimal_cfg["input"]["vcf"] = f"{plain}.gz"
    minimal_cfg["plugins"] = {}
    with pytest.raises(jobs.JobError, match=r"no length for contig\(s\) \['1'\]"):
        jobs.render_run(minimal_cfg, [])
    minimal_cfg["input"]["reference_fasta"] = str(fasta_factory({"1": 700}))
    plan = jobs.render_run(minimal_cfg, [])
    assert (plan.chunks[0].workdir / "regions.bed").read_text() == "1\t0\t700\n"


def _reconciled(job, n_in):
    """A chunk that 'finished' cleanly: an out.vcf.gz and a reconcile.tsv with in == out, unexplained 0."""
    job.workdir.mkdir(parents=True, exist_ok=True)
    job.out_vcf.write_bytes(b"")
    (job.workdir / "reconcile.tsv").write_text(f"in\tout\tskipped\tno_alt_gt\tunexplained\n{n_in}\t{n_in}\t0\t0\t0\n")


def test_finalize_run_detects_records_lost_before_the_chunks(plan):
    """I2(a): every chunk reconciles cleanly against its own in.vcf.gz, yet a record never reached any chunk."""
    for job, n in zip(plan.chunks, (1, 1, 1, 0)):            # the input holds 4 records: 1, 2, 3, KI270
        _reconciled(job, n)
    with pytest.raises(jobs.JobError, match=r"chunk inputs hold 3 record\(s\), but .*S1.input.vcf.gz has 4"):
        jobs.finalize_run(plan)
    assert not plan.final_vcf.exists()


def test_finalize_run_detects_records_lost_before_the_sv_job(cfg, registry, vcf_factory):
    cfg["sv"] = {"enabled": True, "vcf": str(vcf_factory({"1": 1000}, [("1", 10, "N", "<DEL>"), ("1", 50, "N", "<DUP>")]))}
    plan = jobs.render_run(cfg, enabled_plugins(cfg, registry))
    assert plan.sv_input == Path(cfg["output"]["dir"]) / "S1.sv.input.vcf.gz"
    for job in plan.chunks:
        _reconciled(job, 1)
    _reconciled(plan.sv, 1)
    with pytest.raises(jobs.JobError, match=r"chunk inputs hold 1 record\(s\), but .*S1.sv.input.vcf.gz has 2"):
        jobs.finalize_run(plan)
    assert not plan.final_vcf.exists()
    assert jobs.plan_from_disk(cfg, enabled_plugins(cfg, registry)).sv_input == plan.sv_input


# ----------------------------------------------------------------------------- run_local

def _fake_plan(tmp_path, scripts):
    chunks = []
    for name, body in scripts.items():
        wd = tmp_path / "work" / name
        wd.mkdir(parents=True)
        (wd / "run.sh").write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + body)
        chunks.append(jobs.ChunkJob(name=name, contigs=(name,), workdir=wd, ini=wd / "vep.ini", script=wd / "run.sh",
                                    in_vcf=wd / "in.vcf.gz", out_vcf=wd / "out.vcf.gz"))
    return jobs.RunPlan(workdir=tmp_path / "work", chunks=chunks, merge_script=tmp_path / "work" / "merge.sh",
                        final_vcf=tmp_path / "final.vcf.gz", env="vep116", sample="S1")


def test_run_local_failed_chunk_stops_loudly(tmp_path):
    plan = _fake_plan(tmp_path, {"1": "exit 0\n", "2": "echo boom >&2; exit 3\n"})
    with pytest.raises(jobs.JobError, match=r"chunk 2.*exit 3") as e:
        jobs.run_local(plan, parallel=2)
    assert "boom" in str(e.value)
    assert not plan.final_vcf.exists()


# ----------------------------------------------------------------------------- CLI

@pytest.fixture
def cfg_path(cfg, tmp_path):
    path = tmp_path / "cfg.yaml"
    path.write_text(yaml.safe_dump(cfg))
    return path


def test_cli_render_writes_scripts(cfg_path, cfg, capsys):
    assert cli.main(["render", "-c", str(cfg_path)]) == 0
    wd = Path(cfg["output"]["dir"]) / "work"
    assert (wd / "1" / "vep.ini").is_file() and (wd / "other" / "run.sh").is_file()
    assert (wd / "contigs.txt").is_file()


def test_cli_run_slurm_prints_commands_without_submitting(cfg_path, cfg, capsys, monkeypatch):
    real = subprocess.run

    def guard(cmd, *a, **kw):
        assert cmd[0] != "sbatch", "sbatch must not run without --submit"
        return real(cmd, *a, **kw)

    monkeypatch.setattr(jobs.subprocess, "run", guard)
    assert cli.main(["run", "-c", str(cfg_path), "--slurm"]) == 0
    out = capsys.readouterr().out
    assert "sbatch" in out and "--dependency=afterok:" in out and "array.sbatch" in out and "merge.sbatch" in out


def test_cli_run_slurm_submit_runs_sbatch(cfg_path, capsys, monkeypatch):
    real = subprocess.run
    calls = []

    def fake(cmd, *a, **kw):
        if cmd[0] == "sbatch":
            calls.append(list(cmd))
            return subprocess.CompletedProcess(cmd, 0, f"{700 + len(calls)}\n", "")
        return real(cmd, *a, **kw)

    monkeypatch.setattr(jobs.subprocess, "run", fake)
    assert cli.main(["run", "-c", str(cfg_path), "--slurm", "--submit"]) == 0
    assert len(calls) == 2 and "--dependency=afterok:701" in calls[1]
    assert "701" in capsys.readouterr().out


def test_cli_run_needs_a_mode(cfg_path, capsys):
    assert cli.main(["run", "-c", str(cfg_path)]) == 2
    assert "--local" in capsys.readouterr().err
