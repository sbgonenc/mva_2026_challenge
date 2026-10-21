"""Run VEP: prepare the input, render per-chunk jobs, run them locally or on SLURM, reconcile, merge.

The functions that call bcftools/tabix directly (reconcile, merge, header checks) expect to run inside the
tools env (`bin/vep-pipeline` and the generated SLURM scripts guarantee that); everything that needs a *different*
env (chunk scripts, the Exomiser-side normalisation) goes through `micromamba run -n <env>`.

Layout under `{output.dir}`: `{prefix}.input.vcf.gz` (one sample), `work/contigs.txt` (chunk names, one per line;
the SLURM array task N runs line N), `work/<chunk>/` (vep.ini, run.sh, in.vcf.gz, out.vcf.gz, skipped.txt,
reconcile.tsv, run.log), `work/sv/` for the optional SV run, `{prefix}.vep.vcf.gz`, `fields.tsv`.
"""
from __future__ import annotations
import copy
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, replace
from pathlib import Path

import yaml

from . import download as dl
from .checks import _run, _run_ok
from .config import ConfigError
from .registry import _plugin_fields, report_dir, write_fields_tsv
from .vepcmd import build_ini, chunk_script, contig_chunks, header_contigs, regions_bed

PKG_ROOT = Path(__file__).resolve().parents[1]
SV_CHUNK = "sv"
SV_MAX_SIZE = 100_000_000
EXAMPLES = 20
_SKIPPED_LINE = re.compile(r"^Line\s+(\d+)\b")
_FORMAT = re.compile(r'Format:\s*([^"]+)"')


class JobError(RuntimeError):
    """A failed chunk, merge or submission; always stops the run."""


class ReconcileError(JobError):
    """Input and output records of a chunk do not add up."""


@dataclass(frozen=True)
class ChunkJob:
    name: str
    contigs: tuple
    workdir: Path
    ini: Path
    script: Path
    in_vcf: Path
    out_vcf: Path
    kind: str = "snv"                       # "snv" | "sv"

    @property
    def skipped_txt(self) -> Path:
        return self.workdir / "skipped.txt"

    @property
    def regions(self) -> Path:
        return self.workdir / "regions.bed"


@dataclass(frozen=True)
class RunPlan:
    workdir: Path
    chunks: list
    merge_script: Path
    final_vcf: Path
    env: str = "vep116"
    sample: str = ""
    individual_mode: bool = False
    expected_fields: tuple = ()
    plugins: tuple = ()
    input_vcf: Path | None = None
    cfg_path: Path | None = None
    sv: ChunkJob | None = None
    sv_final: Path | None = None
    sv_input: Path | None = None

    @property
    def jobs(self) -> list:
        return list(self.chunks) + ([self.sv] if self.sv else [])


@dataclass(frozen=True)
class ReconcileReport:
    n_in: int
    n_out: int
    skipped: int
    no_alt_gt: int
    unexplained: int
    examples: tuple = ()


# ----------------------------------------------------------------------------- small helpers

def _out_dir(cfg) -> Path:
    return Path(cfg["output"]["dir"])


def _prefix(cfg) -> str:
    return cfg["output"]["prefix"]


def work_dir(cfg) -> Path:
    return _out_dir(cfg) / "work"


def final_vcf_path(cfg) -> Path:
    return _out_dir(cfg) / f"{_prefix(cfg)}.vep.vcf.gz"


def sv_final_path(cfg) -> Path:
    return _out_dir(cfg) / f"{_prefix(cfg)}.sv.vep.vcf.gz"


def sv_input_path(cfg) -> Path:
    return _out_dir(cfg) / f"{_prefix(cfg)}.sv.input.vcf.gz"


def is_individual_mode(plugins) -> bool:
    """G2P forces VEP's `individual` output (variants where the sample carries an alt allele)."""
    return any(p.spec.module == "G2P" for p in plugins)


def _tool(cmd, env=None) -> str:
    """stdout of a tool; `env` None runs it directly (we are already inside the tools env)."""
    cmd = list(map(str, cmd))
    r = subprocess.run(cmd if env is None else ["micromamba", "run", "-n", env, *cmd],
                       capture_output=True, text=True)
    if r.returncode:
        raise JobError(f"{' '.join(cmd[:3])} failed (exit {r.returncode}): {r.stderr.strip()[-500:]}")
    return r.stdout


def _tabix(vcf, env=None):
    _tool(["tabix", "-f", "-p", "vcf", vcf], env)


def _count(vcf, env=None) -> int:
    return sum(1 for _ in _tool(["bcftools", "query", "-f", "x\n", vcf], env).splitlines())


def _key(chrom, pos, ref, alt) -> str:
    return f"{chrom}:{pos}:{ref}>{alt}"


# ----------------------------------------------------------------------------- input

def _fresh(path: Path) -> bool:
    """A finished bgzipped VCF: non-empty, with a .tbi at least as new as the file."""
    idx = Path(f"{path}.tbi")
    return (path.is_file() and idx.is_file() and path.stat().st_size > 0
            and idx.stat().st_mtime >= path.stat().st_mtime)


def _sha256(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def _norm_fingerprint(cfg, src: Path) -> dict:
    """Everything the normalised VCF depends on. `-c` (bcftools threads) does not change the output, so it is not part
    of it; the FASTA is identified by path and size (hashing 3 GB on every run is too slow). A missing file gives
    None and is left to prepare-patient-vcf.sh to report."""
    fasta, script = Path(cfg["input"]["reference_fasta"]), Path(cfg["tools"]["prepare_vcf"])
    st = src.stat() if src.is_file() else None
    return {"source": os.path.realpath(src), "source_size": st and st.st_size, "source_mtime_ns": st and st.st_mtime_ns,
            "fasta": os.path.realpath(fasta), "fasta_size": fasta.stat().st_size if fasta.is_file() else None,
            "prepare_vcf": os.path.realpath(script), "prepare_vcf_sha256": _sha256(script),
            "sample": cfg["sample"]}


def _read_json(path: Path):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


def _subset(cfg, src: Path, dest: Path) -> Path:
    """One-sample bgzipped, tabix-indexed copy of `src` (VEP never sees a multi-sample VCF)."""
    env, sample = cfg["tools"]["env"], cfg["sample"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    _run_ok(["bcftools", "view", "-s", sample, "-Oz", "-o", dest, src], env)
    _run_ok(["tabix", "-f", "-p", "vcf", dest], env)
    return dest


def prepare_input(cfg) -> Path:
    """Optionally normalise (Exomiser's prepare-patient-vcf.sh), then always subset to cfg['sample'].

    `{prefix}.norm.vcf.gz` is reused only when `{prefix}.norm.vcf.gz.fingerprint.json` (written after a successful
    normalisation) matches `_norm_fingerprint` exactly; otherwise it is redone."""
    out, prefix = _out_dir(cfg), _prefix(cfg)
    src = Path(cfg["input"]["vcf"])
    if cfg["input"].get("normalize"):
        out.mkdir(parents=True, exist_ok=True)
        norm = out / f"{prefix}.norm.vcf.gz"
        stamp = Path(f"{norm}.fingerprint.json")
        want = _norm_fingerprint(cfg, src)
        if not (_fresh(norm) and _read_json(stamp) == want):
            stamp.unlink(missing_ok=True)        # an interrupted run must never leave a matching fingerprint
            r = _run(["bash", cfg["tools"]["prepare_vcf"], "-i", src, "-r", cfg["input"]["reference_fasta"],
                      "-o", norm, "-s", cfg["sample"], "-c", (cfg.get("run") or {}).get("fork", 1)],
                     cfg["tools"]["prepare_env"])
            if r.returncode:
                raise JobError(f"prepare-patient-vcf.sh failed (exit {r.returncode}): {(r.stderr or r.stdout).strip()[-800:]}")
            stamp.write_text(json.dumps(want, indent=1, sort_keys=True) + "\n")
        src = norm
    return _subset(cfg, src, out / f"{prefix}.input.vcf.gz")


# ----------------------------------------------------------------------------- rendering

def _job(cfg, name, contigs, kind="snv") -> ChunkJob:
    wd = work_dir(cfg) / name
    return ChunkJob(name=name, contigs=tuple(contigs), workdir=wd, ini=wd / "vep.ini", script=wd / "run.sh",
                    in_vcf=wd / "in.vcf.gz", out_vcf=wd / "out.vcf.gz", kind=kind)


def _expected_fields(cfg, plugins) -> tuple:
    fields = [f for rp in plugins for f, _cls in _plugin_fields(rp)]
    fields += [f"{name}*" for name, spec in (cfg.get("custom") or {}).items() if spec.get("enabled")]
    return tuple(dict.fromkeys(fields))


def _contig_lengths(cfg, vcf, contigs) -> dict:
    """{contig: length} for `contigs`, from the VCF header, else from the reference FASTA .fai."""
    lengths = {c: n for c, n in header_contigs(vcf, cfg["tools"]["env"]).items() if n}
    fai = Path(f"{cfg['input']['reference_fasta']}.fai")
    if any(c not in lengths for c in contigs) and fai.is_file():
        for line in fai.read_text().splitlines():
            name, length = (line.split("\t") + [""])[:2]
            if name in contigs and name not in lengths:
                lengths[name] = int(length)
    missing = [c for c in contigs if c not in lengths]
    if missing:
        raise JobError(f"{vcf}: no length for contig(s) {missing} in its header or in {fai}; "
                       f"the chunk regions need whole-contig lengths")
    return lengths


def _write_job(cfg, job: ChunkJob, plugins, custom_paths, vcf, lengths, ini_cfg=None):
    # {report_dir} comes from the cfg rendered here, not from the one the plugins were resolved against (smoke
    # reuses the production plugins), so the ini's report paths and the script's `mkdir -p` are the same dir.
    plugins = [replace(p, report_dir=report_dir(cfg)) for p in plugins]
    job.workdir.mkdir(parents=True, exist_ok=True)
    ini = build_ini(ini_cfg or cfg, plugins, custom_paths, job.name, job.workdir)
    job.ini.write_text(ini)
    job.regions.write_text(regions_bed(job.contigs, lengths))
    job.script.write_text(chunk_script(cfg, job.name, job.regions, job.ini, vcf, job.workdir))


def _write_config(cfg, wd: Path) -> Path:
    path = wd / "config.resolved.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return path


def _merge_sh(cfg, wd, cfg_path) -> Path:
    path = wd / "merge.sh"
    path.write_text(f"""#!/usr/bin/env bash
# Merge the reconciled chunks into {final_vcf_path(cfg)}; run inside the {cfg['tools']['env']} env.
set -euo pipefail
export PYTHONPATH={shlex.quote(str(PKG_ROOT))}${{PYTHONPATH:+:$PYTHONPATH}}
python -m vep_pipeline.cli run --merge-only -c {shlex.quote(str(cfg_path))}
""")
    return path


def _plan(cfg, plugins, names, chunks_contigs, input_vcf, cfg_path, merge_script) -> RunPlan:
    sv = _job(cfg, SV_CHUNK, (), "sv") if SV_CHUNK in names else None
    chunks = [_job(cfg, n, chunks_contigs.get(n, ())) for n in names if n != SV_CHUNK]
    return RunPlan(workdir=work_dir(cfg), chunks=chunks, merge_script=merge_script, final_vcf=final_vcf_path(cfg),
                   env=cfg["tools"]["env"], sample=cfg["sample"], individual_mode=is_individual_mode(plugins),
                   expected_fields=_expected_fields(cfg, plugins), plugins=tuple(plugins), input_vcf=input_vcf,
                   cfg_path=cfg_path, sv=sv, sv_final=sv_final_path(cfg) if sv else None,
                   sv_input=sv_input_path(cfg) if sv else None)


def render_run(cfg, plugins, input_vcf=None) -> RunPlan:
    """Prepare the input (unless given) and write vep.ini + run.sh for every chunk, contigs.txt and merge.sh."""
    input_vcf = Path(input_vcf) if input_vcf else prepare_input(cfg)
    chunks = contig_chunks(input_vcf, cfg["tools"]["env"])
    if not chunks:
        raise JobError(f"{input_vcf} has no variant records")
    lengths = _contig_lengths(cfg, input_vcf, [c for _n, cs in chunks for c in cs])
    wd = work_dir(cfg)
    wd.mkdir(parents=True, exist_ok=True)
    custom_paths = dl.prepare_custom(cfg)
    cfg_path = _write_config(cfg, wd)
    sv_script = render_sv_run(cfg)
    names = [n for n, _ in chunks] + ([SV_CHUNK] if sv_script else [])
    plan = _plan(cfg, plugins, names, dict(chunks), input_vcf, cfg_path, _merge_sh(cfg, wd, cfg_path))
    for job in plan.chunks:
        _write_job(cfg, job, plugins, custom_paths, input_vcf, lengths)
    (wd / "contigs.txt").write_text("".join(f"{n}\n" for n in names))
    return plan


def plan_from_disk(cfg, plugins) -> RunPlan:
    """The plan of an already rendered run (what a SLURM task or the merge job sees); runs no tool."""
    wd = work_dir(cfg)
    names = (wd / "contigs.txt").read_text().split()
    if not names:
        raise JobError(f"{wd / 'contigs.txt'} is empty: run `render` first")
    cfg_path = wd / "config.resolved.yaml"
    return _plan(cfg, plugins, names, {}, _out_dir(cfg) / f"{_prefix(cfg)}.input.vcf.gz", cfg_path, wd / "merge.sh")


def render_sv_run(cfg) -> Path | None:
    """Subset the SV VCF to the sample and write its ini + script (one job, no chunking); None if sv is off.

    Returns the script path. The ini is the SNV ini (same core options, `everything 1` included) plus
    max_sv_size, with no plugins and no custom tracks."""
    sv = cfg.get("sv") or {}
    if not sv.get("enabled"):
        return None
    if not sv.get("vcf"):
        raise ConfigError("sv.enabled is true but sv.vcf is not set")
    src = Path(sv["vcf"])
    if not src.is_file():
        raise ConfigError(f"sv.vcf: file not found: {src}")
    in_vcf = _subset(cfg, src, sv_input_path(cfg))
    contigs = [c for _n, cs in contig_chunks(in_vcf, cfg["tools"]["env"]) for c in cs]
    if not contigs:
        raise JobError(f"{in_vcf} has no variant records")
    job = _job(cfg, SV_CHUNK, contigs, "sv")
    sv_cfg = copy.deepcopy(cfg)
    sv_cfg["custom"] = {}
    sv_cfg["run"] = {**(cfg.get("run") or {}),
                     "extra_flags": [*((cfg.get("run") or {}).get("extra_flags") or []),
                                     f"max_sv_size={SV_MAX_SIZE}"]}
    # No plugins on SVs: DosageSensitivity crashes on SVs with regulatory on (StructuralVariationOverlapAllele has no transcript); user chose regulatory over dosage.
    _write_job(cfg, job, [], {}, in_vcf, _contig_lengths(cfg, in_vcf, contigs), ini_cfg=sv_cfg)
    return job.script


# ----------------------------------------------------------------------------- reconcile

def _skipped_lines(skipped_txt) -> set:
    """Record numbers (1-based among the VCF records) VEP listed in skipped_variants_file.

    VEP truncates the echoed line to 40 characters, so the `Line N` counter is the only exact reference."""
    p = Path(skipped_txt)
    if not p.is_file():
        return set()
    return {int(m.group(1)) for line in p.read_text().splitlines() if (m := _SKIPPED_LINE.match(line))}


def _has_alt(gt: str) -> bool:
    return any(a not in ("0", ".", "") for a in re.split(r"[/|]", gt))


def _write_reconcile(out_vcf, rep: ReconcileReport):
    Path(out_vcf).parent.joinpath("reconcile.tsv").write_text(
        "in\tout\tskipped\tno_alt_gt\tunexplained\n"
        f"{rep.n_in}\t{rep.n_out}\t{rep.skipped}\t{rep.no_alt_gt}\t{rep.unexplained}\n")


def _raise_unexplained(rep: ReconcileReport, in_vcf, out_vcf):
    raise ReconcileError(
        f"{out_vcf}: {rep.unexplained} record(s) of {in_vcf} are unexplained (in={rep.n_in} out={rep.n_out} "
        f"skipped={rep.skipped} no_alt_gt={rep.no_alt_gt}); first {min(EXAMPLES, rep.unexplained)}: "
        f"{', '.join(rep.examples)}")


def reconcile_chunk(in_vcf, out_vcf, skipped_txt, sample, individual_mode: bool, env=None) -> ReconcileReport:
    """Every (CHROM,POS,REF,ALT) of the input must be in the output, skipped by VEP, or (individual mode only)
    have no alt genotype for `sample`; the output must not hold anything the input lacks. Writes reconcile.tsv
    next to `out_vcf` and raises ReconcileError (naming up to 20 keys) when anything is unexplained."""
    rows = [l.split("\t") for l in _tool(
        ["bcftools", "query", "-s", sample, "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", in_vcf], env).splitlines()]
    out_keys = Counter(_key(*l.split("\t")) for l in _tool(
        ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", out_vcf], env).splitlines())
    n_out = sum(out_keys.values())
    skipped_at = _skipped_lines(skipped_txt)
    skipped = no_alt = 0
    unexplained = []
    for i, (chrom, pos, ref, alt, *gts) in enumerate(rows, 1):
        key = _key(chrom, pos, ref, alt)
        if out_keys[key] > 0:
            out_keys[key] -= 1
        elif i in skipped_at:
            skipped += 1
        elif individual_mode and not any(_has_alt(g) for g in gts):
            no_alt += 1
        else:
            unexplained.append(key)
    unexplained += [f"{k} (only in output)" for k, n in out_keys.items() for _ in range(n)]
    rep = ReconcileReport(len(rows), n_out, skipped, no_alt, len(unexplained), tuple(unexplained[:EXAMPLES]))
    _write_reconcile(out_vcf, rep)
    if rep.unexplained:
        _raise_unexplained(rep, in_vcf, out_vcf)
    return rep


def reconcile_sv(in_vcf, out_vcf, skipped_txt, env=None) -> ReconcileReport:
    """SV run: counts records only (in = out + skipped)."""
    n_in, n_out = _count(in_vcf, env), _count(out_vcf, env)
    skipped = len(_skipped_lines(skipped_txt))
    gap = abs(n_in - skipped - n_out)
    rep = ReconcileReport(n_in, n_out, skipped, 0, gap, (f"{gap} record(s) differ",) if gap else ())
    _write_reconcile(out_vcf, rep)
    if gap:
        _raise_unexplained(rep, in_vcf, out_vcf)
    return rep


def reconcile_job(job: ChunkJob, plan: RunPlan) -> ReconcileReport:
    if job.kind == "sv":
        return reconcile_sv(job.in_vcf, job.out_vcf, job.skipped_txt)
    return reconcile_chunk(job.in_vcf, job.out_vcf, job.skipped_txt, plan.sample, plan.individual_mode)


def _reconciled(job: ChunkJob) -> dict:
    """The counts of a chunk's reconcile.tsv ({"in": n, "out": n, ..., "unexplained": n})."""
    path = job.workdir / "reconcile.tsv"
    if not path.is_file():
        raise JobError(f"chunk {job.name}: no reconcile.tsv ({path}); the chunk did not finish or was not reconciled")
    head, row = [l.split("\t") for l in path.read_text().splitlines()[:2]]
    return {k: int(v) for k, v in zip(head, row)}


def _assert_inputs_covered(jobs, source) -> None:
    """The chunk inputs (reconcile.tsv `in`) must add up to the records of the VCF they were cut from: reconcile
    only compares each chunk with its own in.vcf.gz, so a record the contig selection dropped is caught only here."""
    if source is None:
        raise JobError("the run plan names no input VCF to check the chunk inputs against")
    n_chunks, n_source = sum(_reconciled(j)["in"] for j in jobs), _count(source)
    if n_chunks != n_source:
        raise JobError(f"chunk inputs hold {n_chunks} record(s), but {source} has {n_source}: records were lost or "
                       f"duplicated between the input and the chunks; refusing to merge")


# ----------------------------------------------------------------------------- merge

def _csq_header(vcf, env=None) -> str:
    for line in _tool(["bcftools", "view", "-h", vcf], env).splitlines():
        if line.startswith("##INFO=<ID=CSQ,"):
            return line
    raise JobError(f"{vcf}: no CSQ INFO header (VEP wrote no annotation header)")


def assert_csq_headers_identical(vcfs, expected_fields, env=None) -> None:
    """The CSQ header line must be byte-identical in every VCF and its `Format:` must hold every expected field.

    An expected field ending in `*` is a prefix (custom tracks write `<name>_<field>`)."""
    headers = {str(v): _csq_header(v, env) for v in vcfs}
    ref_vcf, ref = next(iter(headers.items()))
    for vcf, hdr in headers.items():
        if hdr != ref:
            raise JobError(f"CSQ header of {vcf} differs from {ref_vcf}: chunks were annotated with different "
                           f"plugin/field sets, refusing to merge")
    m = _FORMAT.search(ref)
    if not m:
        raise JobError(f"{ref_vcf}: CSQ header has no 'Format:' list")
    present = m.group(1).strip().split("|")
    missing = [f for f in expected_fields
               if not (any(p.startswith(f[:-1]) for p in present) if f.endswith("*") else f in present)]
    if missing:
        raise JobError(f"CSQ Format lacks expected field(s) {missing} (plugin failed to load or wrote no header); "
                       f"present: {present}")


def merge_chunks(vcfs, out, env=None) -> Path:
    """`bcftools concat` in the order given, then tabix; written atomically."""
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    part = out.with_name(out.name + ".part.vcf.gz")
    _tool(["bcftools", "concat", "-Oz", "-o", part, *vcfs], env)
    os.replace(part, out)
    _tabix(out, env)
    return out


def _publish(src: Path, dest: Path, env=None):
    part = dest.with_name(dest.name + ".part.vcf.gz")
    shutil.copyfile(src, part)
    os.replace(part, dest)
    _tabix(dest, env)


def _finish_sv(plan: RunPlan) -> Path:
    assert_csq_headers_identical([plan.sv.out_vcf], ())
    _publish(plan.sv.out_vcf, plan.sv_final)
    return plan.sv_final


def finalize_run(plan: RunPlan) -> Path:
    """Check every chunk was reconciled cleanly and the chunk inputs add up to the input VCF (SV job: to the SV input),
    verify the CSQ headers, concat, index, write fields.tsv."""
    for job in plan.jobs:
        if _reconciled(job)["unexplained"]:
            raise JobError(f"chunk {job.name}: reconcile.tsv reports unexplained records")
        if not job.out_vcf.is_file():
            raise JobError(f"chunk {job.name}: missing {job.out_vcf}")
    _assert_inputs_covered(plan.chunks, plan.input_vcf)
    if plan.sv:
        _assert_inputs_covered([plan.sv], plan.sv_input)
    vcfs = [c.out_vcf for c in plan.chunks]
    assert_csq_headers_identical(vcfs, plan.expected_fields)
    final = merge_chunks(vcfs, plan.final_vcf)
    n_final, n_chunks = _count(final), sum(_count(v) for v in vcfs)
    if n_final != n_chunks:
        raise JobError(f"{final}: {n_final} records, but the chunks hold {n_chunks}")
    write_fields_tsv(plan.plugins, final.parent / "fields.tsv")
    if plan.sv:
        _finish_sv(plan)
    return final


# ----------------------------------------------------------------------------- local run

def _run_job(plan: RunPlan, job: ChunkJob):
    log = job.workdir / "run.log"
    with open(log, "w") as fh:
        r = subprocess.run(["micromamba", "run", "-n", plan.env, "bash", str(job.script)], stdout=fh,
                           stderr=subprocess.STDOUT, env=dl.https_env())
    if r.returncode:
        tail = "\n".join(log.read_text(errors="replace").splitlines()[-15:])
        raise JobError(f"chunk {job.name} failed (exit {r.returncode}); log {log}:\n{tail}")
    reconcile_job(job, plan)


def run_local(plan: RunPlan, parallel: int) -> Path:
    """Run every chunk (and the SV job) with `parallel` workers, reconcile each, then merge. Any failure raises
    after the running chunks have finished; nothing is merged."""
    errors = {}
    with ThreadPoolExecutor(max_workers=max(1, int(parallel))) as pool:
        futs = {pool.submit(_run_job, plan, j): j for j in plan.jobs}
        for fut in as_completed(futs):
            try:
                fut.result()
            except Exception as e:                      # noqa: BLE001 - collected and re-raised below
                errors[futs[fut].name] = e
                for other in futs:
                    other.cancel()
    if errors:
        raise JobError("; ".join(f"{e}" for _n, e in sorted(errors.items())))
    return finalize_run(plan)


def run_sv(cfg) -> Path | None:
    """Render, run, reconcile and publish the SV job alone (None if sv is disabled)."""
    script = render_sv_run(cfg)
    if script is None:
        return None
    job = _job(cfg, SV_CHUNK, (), "sv")
    plan = RunPlan(workdir=work_dir(cfg), chunks=[], merge_script=work_dir(cfg) / "merge.sh",
                   final_vcf=final_vcf_path(cfg), env=cfg["tools"]["env"], sample=cfg["sample"], sv=job,
                   sv_final=sv_final_path(cfg), sv_input=sv_input_path(cfg))
    _run_job(plan, job)
    _assert_inputs_covered([job], plan.sv_input)
    return _finish_sv(plan)


# ----------------------------------------------------------------------------- SLURM

def _sbatch_header(cfg, name, log_pattern, cpus, extra=()) -> str:
    s = cfg["slurm"]
    lines = [f"#SBATCH --job-name={name}", f"#SBATCH --partition={s['partition']}", f"#SBATCH --time={s['time']}",
             f"#SBATCH --mem={s['mem']}", f"#SBATCH --cpus-per-task={cpus}",
             f"#SBATCH --output={s['log_dir']}/{log_pattern}", "#SBATCH --export=ALL", *extra]
    return "#!/usr/bin/env bash\n" + "\n".join(lines) + "\nset -euo pipefail\n"


def render_slurm(plan: RunPlan, cfg) -> tuple[Path, Path]:
    """array.sbatch (task N runs the chunk on line N of contigs.txt, then reconciles it) and merge.sbatch."""
    wd, env, q = plan.workdir, plan.env, shlex.quote
    Path(cfg["slurm"]["log_dir"]).mkdir(parents=True, exist_ok=True)
    cfg_path = plan.cfg_path or wd / "config.resolved.yaml"
    n = len(plan.jobs)
    pypath = f"export PYTHONPATH={q(str(PKG_ROOT))}${{PYTHONPATH:+:$PYTHONPATH}}\n"
    array = wd / "array.sbatch"
    array.write_text(
        _sbatch_header(cfg, "vep", "vep_%A_%a.out", cfg["slurm"]["cpus"], [f"#SBATCH --array=1-{n}"])
        + pypath
        + f"""CHUNK=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {q(str(wd / 'contigs.txt'))})
[[ -n "$CHUNK" ]] || {{ echo "no chunk on line $SLURM_ARRAY_TASK_ID of contigs.txt" >&2; exit 1; }}
TMPDIR=${{SLURM_TMPDIR:-/tmp}}/vep.$SLURM_JOB_ID.$SLURM_ARRAY_TASK_ID
export TMPDIR
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
micromamba run -n {env} bash {q(str(wd))}/"$CHUNK"/run.sh
micromamba run -n {env} python -m vep_pipeline.cli reconcile --chunk "$CHUNK" -c {q(str(cfg_path))}
""")
    merge = wd / "merge.sbatch"
    merge.write_text(_sbatch_header(cfg, "vep_merge", "vep_merge_%j.out", 1)
                     + f"micromamba run -n {env} bash {q(str(plan.merge_script))}\n")
    return array, merge


def submit_slurm(array_script, merge_script) -> tuple[str, str]:
    """sbatch the array, then the merge job with --dependency=afterok:<array id>; returns both job ids."""
    def sbatch(*args):
        r = subprocess.run(["sbatch", "--parsable", *args], capture_output=True, text=True)
        if r.returncode:
            raise JobError(f"sbatch {' '.join(args)} failed (exit {r.returncode}): {r.stderr.strip()}")
        return r.stdout.strip().split(";")[0]
    job_id = sbatch(str(array_script))
    merge_id = sbatch(f"--dependency=afterok:{job_id}", "--kill-on-invalid-dep=yes", str(merge_script))
    return job_id, merge_id


def submit_one(script) -> str:
    """sbatch --parsable <script>; returns the job id."""
    r = subprocess.run(["sbatch", "--parsable", str(script)], capture_output=True, text=True)
    if r.returncode:
        raise JobError(f"sbatch {script} failed (exit {r.returncode}): {r.stderr.strip()}")
    return r.stdout.strip().split(";")[0]
