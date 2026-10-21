"""Build the per-chunk VEP ini file and the bash script that runs one chunk."""
from __future__ import annotations
import shlex
import subprocess
from pathlib import Path

from .registry import report_dir

PRIMARY = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
OTHER = "other"

# Plugin errors VEP only logs as warnings; the script turns them into a failure.
WARNING_PATTERN = "Failed to|went wrong|did not return"
# Removed first thing by every chunk run: a rerun that dies early must not leave an older out.vcf.gz plus a clean
# reconcile.tsv that `run --merge-only` would merge.
STALE_FILES = ("out.vcf.gz", "out.vcf.gz.tbi", "reconcile.tsv", "in.vcf.gz", "in.vcf.gz.tbi", "warnings.txt",
               "skipped.txt")


def _is_primary(contig: str) -> bool:
    name = contig[3:] if contig.startswith("chr") else contig
    return name in PRIMARY or contig == "chrM"


def _flag_line(flag: str) -> str:
    key, sep, value = str(flag).lstrip("-").partition("=")
    return f"{key} {value}" if sep else f"{key} 1"


def custom_arg(name, spec, path) -> str:
    """The value of one `custom` ini line. VCF tracks use coords=0 (original VCF coordinates)."""
    fmt = spec.get("format", "vcf")
    parts = [f"file={path}", f"short_name={name}", f"format={fmt}", f"type={spec.get('type', 'exact')}"]
    if fmt == "vcf":
        parts.append(f"coords={spec.get('coords', 0)}")
    fields = spec.get("fields")
    if fields:
        parts.append("fields=" + "%".join(fields))
    return ",".join(parts)


def build_ini(cfg, plugins, custom_paths, chunk: str, workdir) -> str:
    """The vep.ini of one chunk: no comments, one `plugin` line per plugin, one `custom` line per enabled track."""
    wd = Path(workdir)
    run = cfg.get("run") or {}
    lines = [
        "offline 1",
        "cache 1",
        "safe 1",
        f"dir_cache {cfg['data']['cache_dir']}",
        f"cache_version {cfg['cache_version']}",
        "assembly GRCh38",
        "species homo_sapiens",
        f"fasta {cfg['input']['reference_fasta']}",
        f"dir_plugins {cfg['tools']['plugins_dir']}",
        "everything 1",
        "hgvsg 1",
        "flag_pick 1",
        "vcf 1",
        "compress_output bgzip",
        "force_overwrite 1",
        f"fork {run.get('fork', 1)}",
        f"buffer_size {run.get('buffer_size', 5000)}",
        f"skipped_variants_file {wd}/skipped.txt",
        f"warning_file {wd}/warnings.txt",
        f"stats_file {wd}/stats.html",
    ]
    lines += [_flag_line(f) for f in run.get("extra_flags") or []]
    lines += [f"plugin {p.render(chunk, wd)}" for p in plugins]
    for name, spec in (cfg.get("custom") or {}).items():
        if spec.get("enabled"):
            lines.append(f"custom {custom_arg(name, spec, custom_paths[name])}")
    return "\n".join(lines) + "\n"


def _micromamba(env, args):
    r = subprocess.run(["micromamba", "run", "-n", env, *args], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"{' '.join(args)} failed (exit {r.returncode}): {r.stderr.strip()}")
    return r.stdout


def header_contigs(vcf, env) -> dict:
    """{contig: length or None} from the ##contig lines of the VCF header, in header order."""
    out = {}
    for line in _micromamba(env, ["bcftools", "view", "-h", str(vcf)]).splitlines():
        if line.startswith("##contig=<"):
            f = dict(kv.split("=", 1) for kv in line[len("##contig=<"):].rstrip(">").split(",") if "=" in kv)
            if "ID" in f:
                out[f["ID"]] = int(f["length"]) if f.get("length", "").isdigit() else None
    return out


def regions_bed(contigs, lengths: dict) -> str:
    """BED text (0-based, half-open) selecting each contig whole. `bcftools view -R x.bed` reads it as columns, so
    contig names with ':' or '*' (HLA alts) are never parsed as region strings or shell globs."""
    return "".join(f"{c}\t0\t{lengths[c]}\n" for c in contigs)


def contig_chunks(vcf, env) -> list[tuple[str, list[str]]]:
    """[(chunk, [contigs])] in VCF header order: one chunk per primary contig, then one ("other", [...]).

    Contigs without records are dropped. Needs a .tbi/.csi index (`bcftools index -s`).
    """
    counts = {}
    for line in _micromamba(env, ["bcftools", "index", "-s", str(vcf)]).splitlines():
        name, _length, n = (line.split("\t") + ["0", "0"])[:3]
        if name and int(n) > 0:
            counts[name] = int(n)
    order = list(header_contigs(vcf, env))
    order += [c for c in counts if c not in order]        # indexed but not declared in the header
    present = [c for c in order if c in counts]
    chunks = [(c, [c]) for c in present if _is_primary(c)]
    other = [c for c in present if not _is_primary(c)]
    if other:
        chunks.append((OTHER, other))
    return chunks


def chunk_script(cfg, chunk, regions, ini, vcf, wd) -> str:
    """Bash for one chunk, run by Task 6 inside `micromamba run -n <env> bash <script>` (so no micromamba here).

    `regions` is the chunk's BED file (`regions_bed`) of whole contigs."""
    q = lambda p: shlex.quote(str(p))
    return f"""#!/usr/bin/env bash
set -euo pipefail
# chunk {chunk}
export LC_ALL=C PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0
for v in $(env | grep -o '^VEP_[A-Z_]*' | grep -v '^VEP_DIR_PLUGINS$' || true); do unset "$v"; done
cd {q(wd)}
rm -rf {q(Path(wd) / 'g2p_log')}
rm -f {' '.join(STALE_FILES)}
TMPDIR=$(mktemp -d "${{TMPDIR:-$PWD}}/vep_tmp.XXXXXX")
export TMPDIR
trap 'rm -rf "$TMPDIR"' EXIT
mkdir -p {q(report_dir(cfg))}
# bcftools (-Oz) and VEP (compress_output bgzip) write BGZF themselves, so there is no separate `bgzip -f` step.
bcftools view -R {q(regions)} -Oz -o in.vcf.gz {q(vcf)}
tabix -f -p vcf in.vcf.gz
vep --config {q(ini)} -i in.vcf.gz -o out.vcf.gz
tabix -f -p vcf out.vcf.gz
if [[ -s warnings.txt ]] && grep -Eq '{WARNING_PATTERN}' warnings.txt; then
  echo "plugin error in warnings.txt" >&2; exit 3
fi
"""
