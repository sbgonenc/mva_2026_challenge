"""Preflight checks: each returns a CheckResult and never raises.

Convention: `ok` False with level "error" fails `cli check`; level "warn" is shown but does not fail it;
level "info" is informational. Every external command goes through `_run` (micromamba env + CA bundle).
"""
from __future__ import annotations
import dataclasses
import functools
import gzip
import os
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from . import download as dl
from .registry import _plugin_fields
from .vepcmd import WARNING_PATTERN, build_ini

PARAMS_FILE = Path(__file__).with_name("vep_params.txt")

# Test variant: 21:25585733 C>T lies in the test cache (21:25000001-26000000) and in MRPL39 of the real cache.
TEST_CHROM, TEST_POS = "21", 25585733
TEST_VCF = ("##fileformat=VCFv4.2\n"
            "##contig=<ID=21,length=46709983>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
            "21\t25585733\t.\tC\tT\t.\t.\t.\tGT\t0/1\n")
LOAD_TIMEOUT = 900

MISMATCH_LIMIT = 0.02            # (e) REF mismatch rate above which the VCF is not on this FASTA
REF_SAMPLE = 2000                # (e) records compared
# (g) GRCh38 non-PAR chrX / chrY (PAR1 1-2781479, PAR2 155701383-end on X; PAR1 1-2781479, PAR2 56887903-end on Y)
X_NONPAR = (2_781_480, 155_701_382)
Y_NONPAR = (2_781_480, 56_887_902)
X_HET_FEMALE, X_HET_MALE = 0.20, 0.05
Y_CALLS_MALE, Y_CALLS_FEMALE = 200, 50


def validate_ini_keys(ini_text: str, params_file) -> list[str]:
    """Keys of `ini_text` that VEP would not accept, in ini order.

    `params_file` lists one valid VEP option per line (the set is compared, never the order).
    VEP skips only lines that start with `#` (Config.pm read_config_from_file); a `#` after content becomes part
    of the value, so such a line is reported as `inline-comment:<key>`.
    """
    valid = {l.strip() for l in Path(params_file).read_text().splitlines() if l.strip()}
    bad = []
    for line in ini_text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split(None, 1)
        key = parts[0].lstrip("-")
        if len(parts) > 1 and "#" in parts[1]:
            bad.append(f"inline-comment:{key}")
        elif key not in valid:
            bad.append(key)
    return bad


# ----------------------------------------------------------------------------- result type and helpers

@dataclass(frozen=True)
class CheckResult:
    name: str
    ok: bool
    level: str          # error | warn | info
    message: str


def _guard(name: str, level: str = "error"):
    """Decorator: an exception inside a check becomes a failed CheckResult carrying the error text."""
    def deco(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            try:
                return fn(*args, **kwargs)
            except Exception as e:            # noqa: BLE001 - a check must never raise
                return CheckResult(name, False, level, f"{type(e).__name__}: {e}")
        return wrapper
    return deco


def _child_env(extra=None, strip_vep=False) -> dict:
    env = dl.https_env()
    if strip_vep:
        for k in [k for k in env if k.startswith("VEP_") and k != "VEP_DIR_PLUGINS"]:
            del env[k]
    env.update(extra or {})
    return env


def _run(cmd, env_name, *, extra_env=None, strip_vep=False, timeout=None) -> subprocess.CompletedProcess:
    """`micromamba run -n env_name cmd`, with the CA bundle fixed (https_env); never raises on a non-zero exit."""
    return subprocess.run(["micromamba", "run", "-n", env_name, *map(str, cmd)], capture_output=True, text=True,
                          env=_child_env(extra_env, strip_vep), timeout=timeout)


def _run_ok(cmd, env_name, **kw) -> str:
    r = _run(cmd, env_name, **kw)
    if r.returncode:
        raise RuntimeError(f"{' '.join(map(str, cmd[:3]))} failed (exit {r.returncode}): {r.stderr.strip()[-300:]}")
    return r.stdout


def _stream(cmd, env_name, limit=None):
    """Yield stdout lines of `cmd`; stops (and terminates the command) after `limit` lines.

    A command that fails before the stream is exhausted raises RuntimeError with its stderr."""
    with tempfile.TemporaryFile("w+") as err:
        p = subprocess.Popen(["micromamba", "run", "-n", env_name, *map(str, cmd)], stdout=subprocess.PIPE,
                             stderr=err, text=True, env=_child_env())
        n, cut = 0, False
        try:
            for line in p.stdout:
                yield line.rstrip("\n")
                n += 1
                if limit and n >= limit:
                    cut = True
                    break
        finally:
            if cut or p.poll() is None:
                p.terminate()
            p.stdout.close()
            rc = p.wait()
        if not cut and rc:
            err.seek(0)
            raise RuntimeError(f"{' '.join(map(str, cmd[:3]))} failed (exit {rc}): {err.read().strip()[-300:]}")


def _vcf_contigs(vcf, env) -> dict:
    """{contig: length or None} from the VCF header, in header order."""
    out = {}
    for line in _run_ok(["bcftools", "view", "-h", vcf], env).splitlines():
        if line.startswith("##contig=<"):
            f = dict(kv.split("=", 1) for kv in line[len("##contig=<"):].rstrip(">").split(",") if "=" in kv)
            if "ID" in f:
                out[f["ID"]] = int(f["length"]) if f.get("length", "").isdigit() else None
    return out


def _norm_contig(c: str) -> str:
    c = c[3:] if c.startswith("chr") else c
    return "MT" if c == "M" else c


def _custom_paths(cfg) -> dict:
    """{name: path} of the enabled custom tracks where prepare_custom puts them, without touching the disk."""
    custom_dir = Path(cfg["data"]["custom_data_dir"])
    out = {}
    for name, e in (cfg.get("custom") or {}).items():
        if not e.get("enabled"):
            continue
        if name == "ENCODE_cCRE":
            out[name] = custom_dir / dl.CCRE_OUT
            continue
        path = Path(str(e.get("file") or custom_dir / f"{name}.unset"))
        if e.get("format") == "vcf" and not dl._has_index(path):
            path = custom_dir / path.name                        # the symlink download/prepare_custom indexes
        out[name] = path
    return out


# ----------------------------------------------------------------------------- (a) files

def _file_problems(label, path: Path, indexed: bool) -> list[str]:
    if not path.is_file():
        return [f"{label}: missing: {path}"]
    if path.stat().st_size == 0:
        return [f"{label}: empty: {path}"]
    if indexed and not dl._has_index(path):
        return [f"{label}: no .tbi/.csi index next to {path}"]
    return []


def _dir_problems(label, path: Path, fs) -> list[str]:
    """`_file_problems` for every file of a `glob` entry (a directory of per-chromosome files)."""
    if not path.is_dir():
        return [f"{label}: missing directory: {path}"]
    files = fs.members(path)
    if not files:
        return [f"{label}: no files matching {fs.glob} in {path}"]
    return [p for f in files for p in _file_problems(label, f, bool(fs.index))]


@_guard("files")
def check_files(cfg, plugins) -> CheckResult:
    problems = []
    for rp in plugins:
        for fs in rp.spec.files:
            path = rp.paths.get(fs.key)
            if path is not None and fs.glob:
                problems += _dir_problems(f"{rp.spec.name}.{fs.key}", Path(path), fs)
            elif path is not None:
                problems += _file_problems(f"{rp.spec.name}.{fs.key}", Path(path), bool(fs.index))
    for name, e in (cfg.get("custom") or {}).items():
        if e.get("enabled") and not (e.get("file") or e.get("source")):
            problems.append(f"custom {name}: neither `file` nor `source` is set")
    for name, path in _custom_paths(cfg).items():
        fmt = (cfg["custom"][name].get("format") or "vcf")
        problems += _file_problems(f"custom {name}", path, fmt in ("vcf", "bed", "gff", "gtf"))
    fasta = Path(cfg["input"]["reference_fasta"])
    problems += _file_problems("reference FASTA", fasta, False)
    info = dl._cache_info(cfg)
    if not info.is_file():
        problems.append(f"VEP cache: info.txt missing: {info}")
    if problems:
        return CheckResult("files", False, "error", "; ".join(problems))
    return CheckResult("files", True, "error", "all plugin, custom, FASTA and cache files present")


# ----------------------------------------------------------------------------- (b) ini keys

_PARAMS_PERL = ("for (@Bio::EnsEMBL::VEP::Config::VEP_PARAMS){my($n)=/^([^=!+:]+)/; "
                'print "$_\\n" for split /\\|/, $n}')


def committed_params(params_file=PARAMS_FILE) -> set:
    return {l.strip() for l in Path(params_file).read_text().splitlines() if l.strip()}


def regenerate_params(env) -> set:
    """Option names VEP accepts, read from @Bio::EnsEMBL::VEP::Config::VEP_PARAMS (as a set: file order is locale-dependent)."""
    out = _run_ok(["perl", "-MBio::EnsEMBL::VEP::Config", "-e", _PARAMS_PERL], env)
    return {l.strip() for l in out.splitlines() if l.strip()}


@_guard("ini_keys")
def check_ini_keys(cfg, plugins, workdir, params_file=PARAMS_FILE) -> CheckResult:
    ini = build_ini(cfg, plugins, _custom_paths(cfg), "1", workdir)
    bad = validate_ini_keys(ini, params_file)
    fresh, committed = regenerate_params(cfg["tools"]["env"]), committed_params(params_file)
    problems = []
    if bad:
        problems.append("keys VEP would reject: " + ", ".join(bad))
    if fresh != committed:
        problems.append(f"{Path(params_file).name} is out of date: only in VEP {sorted(fresh - committed)}, "
                        f"only in file {sorted(committed - fresh)}")
    if problems:
        return CheckResult("ini_keys", False, "error", "; ".join(problems))
    return CheckResult("ini_keys", True, "error", f"generated ini uses only valid keys; {len(fresh)} VEP options match")


# ----------------------------------------------------------------------------- (c) plugin load, (d) SpliceAI

def _csq_fields(vcf_gz) -> list[str]:
    with gzip.open(vcf_gz, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("##INFO=<ID=CSQ"):
                m = re.search(r'Format: ([^"]+)"', line)
                if m:
                    return m.group(1).split("|")
    return []


def _vep_msg(r) -> str:
    text = r.stderr.strip() or r.stdout.strip()
    msgs = [l.strip() for l in text.splitlines() if l.startswith("MSG:")]
    tail = msgs[-1] if msgs else " | ".join(text.splitlines()[-3:])
    return " ".join(tail.split())[:400]


def _load_test(cfg, plugin, workdir, template, name) -> CheckResult:
    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    (wd / "reports").mkdir(exist_ok=True)
    rp = dataclasses.replace(plugin, report_dir=wd / "reports",
                             **({"arg_template": template} if template else {}))
    one = {**cfg, "custom": {}, "run": {**(cfg.get("run") or {}), "fork": 1, "buffer_size": 100, "extra_flags": []}}
    (wd / "in.vcf").write_text(TEST_VCF.format(sample=cfg["sample"]))
    (wd / "vep.ini").write_text(build_ini(one, [rp], {}, TEST_CHROM, wd))
    for stale in ("warnings.txt", "skipped.txt", "out.vcf.gz"):
        (wd / stale).unlink(missing_ok=True)
    r = _run(["vep", "--config", wd / "vep.ini", "-i", wd / "in.vcf", "-o", wd / "out.vcf.gz"], cfg["tools"]["env"],
             extra_env={"PERL_HASH_SEED": "0", "PERL_PERTURB_KEYS": "0"}, strip_vep=True, timeout=LOAD_TIMEOUT)
    if r.returncode:
        return CheckResult(name, False, "error", f"VEP exited {r.returncode}: {_vep_msg(r)}")
    warnings = wd / "warnings.txt"
    hits = [l.strip() for l in (warnings.read_text().splitlines() if warnings.is_file() else [])
            if re.search(WARNING_PATTERN, l)]
    if hits:
        return CheckResult(name, False, "error", f"plugin warnings in {warnings}: {hits[0][:300]}")
    present = set(_csq_fields(wd / "out.vcf.gz"))
    missing = [f for f, _cls in _plugin_fields(rp) if f not in present]
    if missing:
        return CheckResult(name, False, "error", f"output field(s) {missing} not in the CSQ Format of the test run")
    return CheckResult(name, True, "error", f"loaded and annotated the test variant ({TEST_CHROM}:{TEST_POS})")


def plugin_load_test(cfg, plugin, workdir, template=None) -> CheckResult:
    """Run VEP for real on the 1-record test VCF with only this plugin (`safe 1`, the configured cache and FASTA).

    Passes if VEP exits 0, warnings.txt has no WARNING_PATTERN line and every output field is in the CSQ Format.
    `template` replaces the plugin's argument string (used by `needs_empty_indel`)."""
    name = f"plugin:{plugin.spec.name}"
    try:
        return _load_test(cfg, plugin, workdir, template, name)
    except subprocess.TimeoutExpired:
        return CheckResult(name, False, "error", f"timed out after {LOAD_TIMEOUT}s")
    except Exception as e:                    # noqa: BLE001
        return CheckResult(name, False, "error", f"{plugin.spec.name}: {type(e).__name__}: {e}")


_EMPTY_INDEL_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def make_empty_indel_vcf(path, env="vep116") -> Path:
    """Header-only VCF, bgzipped to `path` with a .tbi. SpliceAI.pm new() calls add_file(indel) unconditionally."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.parent / f"{path.name}.part.vcf"
    part.write_text(_EMPTY_INDEL_HEADER)
    _run_ok(["bgzip", "-f", part], env)
    Path(f"{part}.gz").replace(path)
    _run_ok(["tabix", "-f", "-p", "vcf", path], env)
    return path


def needs_empty_indel(cfg, spliceai_plugin, workdir) -> bool:
    """True when the plugin does not load with `snv=` alone (so an empty `indel=` file must be passed).

    Verified: SpliceAI 116 dies with 'ERROR: No file specified', which is why the registry always passes indel=."""
    snv_only = re.sub(r",indel=[^,]*", "", spliceai_plugin.arg_template)
    return not plugin_load_test(cfg, spliceai_plugin, Path(workdir) / "snv_only", template=snv_only).ok


@_guard("spliceai_form")
def check_spliceai_form(cfg, plugins, workdir) -> CheckResult:
    spl = [rp for rp in plugins if rp.spec.module == "SpliceAI"]
    if not spl:
        return CheckResult("spliceai_form", True, "error", "SpliceAI not enabled")
    for rp in spl:
        if "indel=" not in rp.arg_template and needs_empty_indel(cfg, rp, workdir):
            return CheckResult("spliceai_form", False, "error",
                               f"{rp.spec.name}: snv-only arguments fail to load; add indel=<file> (an empty "
                               f"header-only VCF works)")
    return CheckResult("spliceai_form", True, "error", "SpliceAI arguments include indel= (snv-only dies with "
                       "'No file specified')" if all("indel=" in rp.arg_template for rp in spl)
                       else "snv-only SpliceAI arguments load")


# ----------------------------------------------------------------------------- (e) VCF vs FASTA

def fasta_seqs(fasta, regions, env) -> list[str]:
    """Upper-case sequence of each `contig:start-end` region (1-based, inclusive), in order (`samtools faidx -r`)."""
    with tempfile.TemporaryDirectory(prefix="vep_check_") as tmp:
        path = Path(tmp) / "regions.txt"
        path.write_text("".join(f"{r}\n" for r in regions))
        out = _run_ok(["samtools", "faidx", "-r", path, fasta], env)
    seqs, cur = [], None
    for line in out.splitlines():
        if line.startswith(">"):
            cur = []
            seqs.append(cur)
        elif cur is not None:
            cur.append(line)
    if len(seqs) != len(regions):
        raise RuntimeError(f"samtools faidx returned {len(seqs)} sequences for {len(regions)} regions")
    return ["".join(s).upper() for s in seqs]


@_guard("vcf_vs_fasta")
def check_vcf_vs_fasta(vcf, fasta, env, n=REF_SAMPLE, limit=MISMATCH_LIMIT) -> CheckResult:
    fai = {}
    for line in Path(f"{fasta}.fai").read_text().splitlines():
        name, length = line.split("\t")[:2]
        fai[_norm_contig(name)] = (name, int(length))
    contigs = _vcf_contigs(vcf, env)
    shared = {c: _norm_contig(c) for c in contigs if _norm_contig(c) in fai}
    if contigs and not shared:
        return CheckResult("vcf_vs_fasta", False, "error", "no VCF header contig matches a FASTA contig "
                           "(after chrN->N, chrM->MT): different assembly or contig naming")
    wrong = [f"{c}: VCF {contigs[c]} vs FASTA {fai[k][1]}" for c, k in shared.items()
             if contigs[c] is not None and contigs[c] != fai[k][1]]
    if wrong:
        return CheckResult("vcf_vs_fasta", False, "error", "contig length mismatch, VCF is on a different "
                           "assembly than the FASTA: " + "; ".join(wrong[:5]))
    recs = []
    for line in _stream(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\n", vcf], env, limit=n):
        chrom, pos, ref = line.split("\t")
        recs.append((chrom, int(pos), ref))
    if not recs:
        return CheckResult("vcf_vs_fasta", True, "error", "contig lengths agree; VCF has no records to compare")
    bad, queries = 0, []
    for chrom, pos, ref in recs:
        entry = fai.get(_norm_contig(chrom))
        if entry is None or pos + len(ref) - 1 > entry[1]:
            bad += 1
        else:
            queries.append((f"{entry[0]}:{pos}-{pos + len(ref) - 1}", ref))
    if queries:
        seqs = fasta_seqs(fasta, [q for q, _ in queries], env)
        bad += sum(s != ref.upper() for s, (_q, ref) in zip(seqs, queries))
    rate = bad / len(recs)
    msg = f"REF mismatch {rate:.1%} ({bad}/{len(recs)} of the first {n} records)"
    if rate > limit:
        return CheckResult("vcf_vs_fasta", False, "error", msg + f" exceeds {limit:.0%}: the VCF does not match "
                           "this FASTA (wrong assembly?)")
    return CheckResult("vcf_vs_fasta", True, "error", msg + "; contig lengths agree")


# ----------------------------------------------------------------------------- (f) sample, (g) sex

@_guard("sample")
def check_sample(vcf, sample, env) -> CheckResult:
    samples = _run_ok(["bcftools", "query", "-l", vcf], env).split()
    if sample not in samples:
        return CheckResult("sample", False, "error", f"sample {sample!r} not in the VCF (samples: {samples})")
    if len(samples) > 1:
        return CheckResult("sample", False, "warn", f"VCF has {len(samples)} samples; it will be subset to "
                           f"{sample!r} before VEP runs (VEP never sees several samples)")
    return CheckResult("sample", True, "error", f"sample {sample!r} present")


def _calls(vcf, sample, env, contigs, names, span):
    """GT strings (called, at least one alt allele) of `sample` in the non-PAR span of the first matching contig."""
    name = next((c for c in contigs if c in names), None)
    if name is None:
        return
    for gt in _stream(["bcftools", "query", "-r", f"{name}:{span[0]}-{span[1]}", "-s", sample, "-f", "[%GT]\n", vcf], env):
        alleles = re.split(r"[/|]", gt)
        if "." not in alleles and any(a != "0" for a in alleles):
            yield alleles


@_guard("sex", "info")
def check_sex(vcf, sample, env) -> CheckResult:
    contigs = list(_vcf_contigs(vcf, env))
    xs = list(_calls(vcf, sample, env, contigs, ("X", "chrX"), X_NONPAR))
    n_y = sum(1 for _ in _calls(vcf, sample, env, contigs, ("Y", "chrY"), Y_NONPAR))
    het = sum(len(set(a)) > 1 for a in xs) / len(xs) if xs else None
    if het is not None and het >= X_HET_FEMALE and n_y < Y_CALLS_MALE:
        verdict = "likely XX"
    elif het is not None and het <= X_HET_MALE and n_y >= Y_CALLS_MALE:
        verdict = "likely XY"
    else:
        verdict = "unclear"
    het_txt = "n/a" if het is None else f"{het:.2f}"
    return CheckResult("sex", True, "info", f"{verdict} (chrX non-PAR: {len(xs)} calls, het fraction {het_txt}; "
                       f"chrY non-PAR calls: {n_y}); compare with the phenopacket sex")


# ----------------------------------------------------------------------------- (h)-(k)

_free_bytes = dl._free_bytes


@_guard("env", "warn")
def check_env_leaks() -> CheckResult:
    leaks = sorted(k for k in os.environ if k.startswith("VEP_") and k != "VEP_DIR_PLUGINS")
    ini = Path.home() / ".vep" / "vep.ini"
    problems = []
    if leaks:
        problems.append(f"VEP_* variables set (VEP may read them): {', '.join(leaks)}")
    if ini.is_file():
        problems.append(f"{ini} exists and can change VEP defaults")
    if problems:
        return CheckResult("env", False, "warn", "; ".join(problems))
    return CheckResult("env", True, "warn", "no VEP_* variables (besides VEP_DIR_PLUGINS) and no ~/.vep/vep.ini")


@_guard("prepare_script")
def check_prepare_script(cfg) -> CheckResult:
    if not cfg["input"].get("normalize"):
        return CheckResult("prepare_script", True, "error", "skipped: input.normalize is false")
    env = cfg["tools"].get("prepare_env", "exomiser")
    r = _run(["bash", cfg["tools"]["prepare_vcf"], "-h"], env)
    if r.returncode == 1 and "Usage:" in r.stdout + r.stderr:
        return CheckResult("prepare_script", True, "error", f"{cfg['tools']['prepare_vcf']} -h prints Usage: (env {env})")
    return CheckResult("prepare_script", False, "error",
                       f"`bash {cfg['tools']['prepare_vcf']} -h` in env {env}: exit {r.returncode}, expected exit 1 "
                       f"with 'Usage:': {(r.stderr or r.stdout).strip()[:200]}")


@_guard("dbnsfp")
def check_dbnsfp(plugins, env) -> CheckResult:
    rp = next((p for p in plugins if p.spec.name == "dbNSFP"), None)
    if rp is None:
        return CheckResult("dbnsfp", True, "error", "dbNSFP not enabled")
    path, columns = rp.paths.get("file"), rp.params.get("columns")
    if path is None or not isinstance(columns, list):
        return CheckResult("dbnsfp", True, "error", "dbNSFP configured through `args`; columns not checked")
    header = [l for l in _run_ok(["tabix", "-H", path], env).splitlines() if l.startswith("#")]
    have = {c.lstrip("#") for c in header[-1].split("\t")} if header else set()
    missing = [c for c in columns if c not in have]
    if missing:
        return CheckResult("dbnsfp", False, "error", f"column(s) {missing} not in the header of {path}")
    return CheckResult("dbnsfp", True, "error", f"all {len(columns)} requested columns are in {Path(path).name}")


@_guard("disk", "warn")
def check_disk(cfg) -> CheckResult:
    need = 3 * Path(cfg["input"]["vcf"]).stat().st_size
    free = _free_bytes(Path(cfg["output"]["dir"]))
    if free < need:
        return CheckResult("disk", False, "warn", f"{free} bytes free under {cfg['output']['dir']}, "
                           f"3x the input VCF is {need}")
    return CheckResult("disk", True, "warn", f"{free} bytes free, 3x the input VCF is {need}")


# ----------------------------------------------------------------------------- all checks

def run_all_checks(cfg, plugins) -> list[CheckResult]:
    """Every check, in table order. Never raises: a broken check is reported as failed."""
    def attempt(name, level, fn):
        try:
            return fn()
        except Exception as e:                # noqa: BLE001
            return CheckResult(name, False, level, f"{type(e).__name__}: {e}")

    env = cfg["tools"]["env"]
    vcf = (cfg.get("input") or {}).get("vcf")
    results = []
    with tempfile.TemporaryDirectory(prefix="vep_check_") as tmp:
        wd = Path(tmp)
        results.append(attempt("files", "error", lambda: check_files(cfg, plugins)))
        results.append(attempt("ini_keys", "error", lambda: check_ini_keys(cfg, plugins, wd / "ini")))
        seen = set()
        for rp in plugins:
            if rp.spec.name not in seen:
                seen.add(rp.spec.name)
                results.append(attempt(f"plugin:{rp.spec.name}", "error",
                                       lambda rp=rp: plugin_load_test(cfg, rp, wd / "plugins" / rp.spec.name)))
        results.append(attempt("spliceai_form", "error", lambda: check_spliceai_form(cfg, plugins, wd / "spliceai")))
        vcf_checks = [
            ("vcf_vs_fasta", "error", lambda: check_vcf_vs_fasta(vcf, cfg["input"]["reference_fasta"], env)),
            ("sample", "error", lambda: check_sample(vcf, cfg["sample"], env)),
            ("sex", "info", lambda: check_sex(vcf, cfg["sample"], env)),
        ]
        for name, level, fn in vcf_checks:
            results.append(attempt(name, level, fn) if vcf
                           else CheckResult(name, False, level, "input.vcf is not set"))
        results.append(attempt("env", "warn", check_env_leaks))
        results.append(attempt("prepare_script", "error", lambda: check_prepare_script(cfg)))
        results.append(attempt("dbnsfp", "error", lambda: check_dbnsfp(plugins, env)))
        results.append(attempt("disk", "warn", lambda: check_disk(cfg)) if vcf
                       else CheckResult("disk", False, "warn", "input.vcf is not set"))
    return results
