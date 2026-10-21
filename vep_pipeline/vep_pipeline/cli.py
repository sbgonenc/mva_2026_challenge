from __future__ import annotations
import argparse, sys
from pathlib import Path

from . import checks as chk, download as dl, exomiser as exo, jobs, phasing, prioritise as pr, rank_settings, smoke as smk
from .config import ConfigError, load_config, _normalise
from .registry import RUNNABLE, enabled_plugins, load_registry, unsupported_reason

DEFAULT_CONFIG = Path(__file__).resolve().parents[1] / "config.yaml"
DEFAULT_REGISTRY = Path(__file__).resolve().parents[1] / "plugins.yaml"

COMMANDS = {
    "list-plugins": "list every plugin in the registry and whether it is enabled",
    "download": "download cache, plugin and custom data under data.root",
    "check": "verify data, tools and the input VCF",
    "smoke": "annotate a small region and report plugin coverage",
    "render": "prepare the input and write ini files and chunk scripts without running VEP",
    "run": "annotate the full VCF: --local N, --slurm [--submit], or --merge-only",
    "reconcile": "check one chunk's output against its input (the SLURM array task runs this)",
    "exomiser": "write the phenopacket and the Exomiser SLURM job (phenotype-only + exome + genome); --submit runs it",
    "prioritise": "rank findings and write the top-N CSV (run it with --slurm on a WGS VCF)",
    "phase": "WhatsHap read-backed phasing of candidate genes: writes the SLURM job (--submit runs it)",
}

# Later tasks register handlers here: name -> callable(cfg, args) -> int
HANDLERS: dict = {}

def list_plugins(cfg, args) -> int:
    """Print every registry plugin with its category, default flag, config state and reason."""
    registry = load_registry(DEFAULT_REGISTRY)
    section = _normalise(cfg.get("plugins"))
    rows = [("plugin", "category", "default", "enabled", "reason")]
    for name, spec in registry.items():
        reason = unsupported_reason(spec) if spec.category not in RUNNABLE else (spec.reason or "")
        rows.append((name, spec.category, "yes" if spec.default else "no",
                     "yes" if section.get(name, {}).get("enabled") else "no", reason))
    for name, entry in section.items():
        if name not in registry:
            rows.append((name, "?", "no", "yes" if entry["enabled"] else "no", "unknown plugin"))
    widths = [max(len(r[i]) for r in rows) for i in range(4)]
    for r in rows:
        print("  ".join(c.ljust(w) for c, w in zip(r[:4], widths)) + "  " + r[4])
    return 0

HANDLERS["list-plugins"] = list_plugins

def download(cfg, args) -> int:
    """Print the download plan; unless --dry-run, fetch plugin data, prepare custom tracks, then the VEP cache."""
    plugins = enabled_plugins(cfg, load_registry(DEFAULT_REGISTRY))
    only = args.only
    if only:
        unknown = [n for n in only if n not in {rp.spec.name for rp in plugins} | {"cache", "custom"}]
        if unknown:
            raise ConfigError(f"--only: {unknown} is not an enabled plugin, 'cache' or 'custom'")
        plugins = [rp for rp in plugins if rp.spec.name in only]
    want = lambda part: not only or part in only
    items = dl.plan_downloads(cfg, plugins, include_skipped=True)
    do_cache = want("cache") and not args.skip_cache
    if do_cache:
        items.append(dl.cache_item(cfg))
    dl.print_plan(cfg, items)
    if args.dry_run:
        return 0
    dl.download_all(cfg, plugins)
    if want("custom"):
        for name, path in dl.prepare_custom(cfg).items():
            print(f"custom {name}: {path}")
    if do_cache:
        dl.download_cache(cfg)
    return 0

HANDLERS["download"] = download

def _status(r) -> str:
    return "ok" if r.ok else {"error": "FAIL", "warn": "WARN"}.get(r.level, "info")


def check(cfg, args) -> int:
    """Run every preflight check, print a table, exit 1 if any error-level check failed."""
    if not (cfg.get("input") or {}).get("vcf"):
        raise ConfigError("input.vcf is not set: `check` needs the proband VCF (set input.vcf in the config)")
    results = chk.run_all_checks(cfg, enabled_plugins(cfg, load_registry(DEFAULT_REGISTRY)))
    rows = [("status", "level", "check", "message")] + [
        (_status(r), r.level, r.name, " ".join(r.message.split())) for r in results]
    widths = [max(len(r[i]) for r in rows) for i in range(3)]
    for r in rows:
        print("  ".join(c.ljust(w) for c, w in zip(r[:3], widths)) + "  " + r[3])
    failed = [r.name for r in results if not r.ok and r.level == "error"]
    print(f"\n{len(results)} checks, {len(failed)} failed" + (f": {', '.join(failed)}" if failed else ""))
    return 1 if failed else 0

HANDLERS["check"] = check

def _plugins(cfg):
    return enabled_plugins(cfg, load_registry(DEFAULT_REGISTRY))


def render(cfg, args) -> int:
    """Subset (and, if input.normalize, normalise) the input, then write every chunk's vep.ini and run.sh."""
    plan = jobs.render_run(cfg, _plugins(cfg))
    print(f"input    {plan.input_vcf}")
    for job in plan.jobs:
        print(f"chunk    {job.name:<8} {job.script}")
    print(f"merge    {plan.merge_script}")
    print(f"final    {plan.final_vcf}")
    return 0

HANDLERS["render"] = render

def run(cfg, args) -> int:
    """--merge-only: merge reconciled chunks. --local N: render, run N chunks at a time, merge.
    --slurm: render and print the sbatch commands; with --submit also run them."""
    plugins = _plugins(cfg)
    if args.merge_only:
        print(jobs.finalize_run(jobs.plan_from_disk(cfg, plugins)))
        return 0
    if args.local is None and not args.slurm:
        raise ConfigError("run needs a mode: --local N, --slurm [--submit] or --merge-only")
    if args.submit and not args.slurm:
        raise ConfigError("--submit only applies to --slurm")
    plan = jobs.render_run(cfg, plugins)
    if args.local is not None:
        print(jobs.run_local(plan, args.local))
        return 0
    array, merge = jobs.render_slurm(plan, cfg)
    if not args.submit:
        print(f"sbatch --parsable {array}")
        print(f"sbatch --parsable --dependency=afterok:<ARRAY_JOB_ID> --kill-on-invalid-dep=yes {merge}")
        print("(not submitted: add --submit)")
        return 0
    job_id, merge_id = jobs.submit_slurm(array, merge)
    print(f"submitted array job {job_id} ({len(plan.jobs)} tasks) and merge job {merge_id} "
          f"(afterok:{job_id}); final output {plan.final_vcf}")
    return 0

HANDLERS["run"] = run

def reconcile(cfg, args) -> int:
    """Reconcile one rendered chunk (by name, as listed in work/contigs.txt); exit 1 on unexplained records."""
    plan = jobs.plan_from_disk(cfg, _plugins(cfg))
    job = next((j for j in plan.jobs if j.name == args.chunk), None)
    if job is None:
        raise ConfigError(f"unknown chunk {args.chunk!r}; chunks: {[j.name for j in plan.jobs]}")
    r = jobs.reconcile_job(job, plan)
    print(f"chunk {job.name}: in={r.n_in} out={r.n_out} skipped={r.skipped} no_alt_gt={r.no_alt_gt} "
          f"unexplained={r.unexplained}")
    return 0

HANDLERS["reconcile"] = reconcile

def smoke(cfg, args) -> int:
    """Annotate the smoke locus through the normal pipeline, print the per-plugin report, exit 1 on any FAIL."""
    rows = smk.run_smoke(cfg, _plugins(cfg), seed=args.seed, parallel=(cfg.get("run") or {}).get("fork", 1))
    table = [smk.REPORT_COLUMNS] + [tuple("" if r[c] is None else str(r[c]) for c in smk.REPORT_COLUMNS)
                                    for r in rows]
    widths = [max(len(r[i]) for r in table) for i in range(len(smk.REPORT_COLUMNS) - 1)]
    for r in table:
        print("  ".join(c.ljust(w) for c, w in zip(r[:-1], widths)) + "  " + r[-1])
    named = lambda status: [r["plugin"] for r in rows if r["status"] == status]
    count = lambda names, what: f"{len(names)} {what}" + (f" ({', '.join(names)})" if names else "")
    failed, unchecked = named("FAIL"), named("unchecked")
    print(f"\n{len(rows)} plugins/tracks, {count(failed, 'failed')}, {count(unchecked, 'unchecked')}"
          f"; report {smk.smoke_dir(cfg) / 'smoke.tsv'}")
    return 1 if failed else 0

HANDLERS["smoke"] = smoke

def exomiser_cmd(cfg, args) -> int:
    """Render the Exomiser job; with --submit, sbatch it."""
    script = exo.render_exomiser(cfg)
    if not args.submit:
        print(f"sbatch --parsable {script}")
        print("(not submitted: add --submit)")
        return 0
    print(f"submitted Exomiser job {jobs.submit_one(script)}; results under {script.parent}")
    return 0

HANDLERS["exomiser"] = exomiser_cmd

def prioritise_cmd(cfg, args) -> int:
    """Rank findings and write the CSV here, or (--slurm) write/submit a job that does it."""
    settings = rank_settings.prioritise_settings(cfg, top_n=args.top_n, output=args.output)
    if args.submit and not args.slurm:
        raise ConfigError("--submit only applies to --slurm")
    if args.slurm:
        extra = ([f"--top-n={args.top_n}"] if args.top_n else []) + ([f"--output={args.output}"] if args.output else []) \
            + (["--sensitivity"] if args.sensitivity else [])
        script = pr.render_slurm(cfg, extra)
        if not args.submit:
            print(f"sbatch --parsable {script}")
            print("(not submitted: add --submit)")
            return 0
        print(f"submitted prioritise job {jobs.submit_one(script)}; output {settings['output']}")
        return 0
    out, cand = pr.run(cfg, settings, sensitivity=args.sensitivity)
    print(f"top {settings['top_n']}: {out}\ncandidates: {cand}")
    return 0

HANDLERS["prioritise"] = prioritise_cmd

def phase_cmd(cfg, args) -> int:
    """Render the phasing job (--submit: sbatch it), or run one of its steps here (--prepare, --finish)."""
    ph = rank_settings.phasing_settings(cfg)
    if not ph["enabled"]:
        raise ConfigError("phasing.enabled is false: set phasing.{enabled, bam, reference, env} to use `phase`")
    if args.prepare:
        phasing.prepare(cfg, ph)
        return 0
    if args.finish:
        phasing.finish(cfg, ph)
        return 0
    script = phasing.render_phase(cfg, ph)
    if not args.submit:
        print(f"sbatch --parsable {script}")
        print("(not submitted: add --submit)")
        return 0
    print(f"submitted phasing job {jobs.submit_one(script)}; sidecar {phasing.readback_path(cfg, ph)}")
    return 0

HANDLERS["phase"] = phase_cmd

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="vep-pipeline", description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    for name, help_ in COMMANDS.items():
        sp = sub.add_parser(name, help=help_)
        sp.add_argument("-c", "--config", default=str(DEFAULT_CONFIG),
                        help="config YAML (default: %(default)s)")
        if name == "download":
            sp.add_argument("--dry-run", action="store_true", help="print the plan and stop")
            sp.add_argument("--only", nargs="+", metavar="NAME",
                            help="only these enabled plugins, and/or 'cache' / 'custom'")
            sp.add_argument("--skip-cache", action="store_true", help="do not download the VEP cache")
        if name == "run":
            sp.add_argument("--local", type=int, metavar="N", help="run locally, N chunks in parallel")
            sp.add_argument("--slurm", action="store_true", help="render SLURM scripts and print the sbatch commands")
            sp.add_argument("--submit", action="store_true", help="with --slurm: actually run sbatch")
            sp.add_argument("--merge-only", action="store_true",
                            help="merge already reconciled chunks (the SLURM merge job runs this)")
        if name == "smoke":
            sp.add_argument("--seed", type=int, default=smk.DEFAULT_SEED, help="control sampling seed (default: %(default)s)")
        if name == "exomiser":
            sp.add_argument("--submit", action="store_true", help="run sbatch (default: print the command)")
        if name == "prioritise":
            sp.add_argument("-o", "--output", help="CSV path (default: {output.dir}/{prefix}.top{N}.csv)")
            sp.add_argument("--top-n", type=int, help="rows in the CSV (default: prioritise.top_n = 10)")
            sp.add_argument("--sensitivity", action="store_true", help="also write the pre-declared sensitivity table")
            sp.add_argument("--slurm", action="store_true", help="write a SLURM job instead of running here")
            sp.add_argument("--submit", action="store_true", help="with --slurm: run sbatch")
        if name == "phase":
            mode = sp.add_mutually_exclusive_group()
            mode.add_argument("--submit", action="store_true", help="run sbatch (default: print the command)")
            mode.add_argument("--prepare", action="store_true", help="job step: regions, contig maps, het-call subset")
            mode.add_argument("--finish", action="store_true", help="job step: write the read-back sidecar")
        if name == "reconcile":
            sp.add_argument("--chunk", required=True, help="chunk name from work/contigs.txt")
    return p

def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        cfg = load_config(args.config)
    except ConfigError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    handler = HANDLERS.get(args.command)
    if handler is None:
        print(f"ERROR: '{args.command}' is not implemented yet", file=sys.stderr)
        return 1
    try:
        return handler(cfg, args)
    except ConfigError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    except jobs.JobError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())
