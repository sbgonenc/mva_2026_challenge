"""Positive-control smoke test: annotate a small locus and check that every plugin put values on its controls.

Controls (all heterozygous, GT 0/1, normalised with `bcftools norm -f <fasta> -m -any`):
  * up to CONTROLS_PER_FILE records drawn from each tabix-backed plugin's *own* data file inside the region
    (registry `control: {source: data|vcf}`),
  * ClinVar records in the region, CLINVAR_PER_STRATUM per consequence stratum (INFO/MC). Any clinical
    significance and review status: a control needs its consequence class, not its pathogenicity.

Each control record carries its provenance in the VCF ID column (`ctl_<Plugin>`, `cv_<stratum>`, `;`-joined), so
`control_keys_from_vcf` rebuilds the control sets from the VCF alone. The smoke run itself is the normal chunk
pipeline (`jobs.render_run` + `jobs.run_local`) on that VCF. `ref_mismatches` finds controls whose REF is not the
FASTA base (`bcftools norm -c w` only warns about them), which flags plugin data on another assembly.
"""
from __future__ import annotations
import copy
import gzip
import re
from pathlib import Path

from .checks import MISMATCH_LIMIT, _run_ok, fasta_seqs
from .config import ConfigError
from .jobs import _key, render_run, run_local
from .registry import _plugin_fields

DEFAULT_REGION = "17:43044295-43125483"        # BRCA1, GRCh38
DEFAULT_SEED = 7
CONTROLS_PER_FILE = 25
CLINVAR_PER_STRATUM = 10
CLINVAR_KEY = "ClinVar"                         # control_keys entry: every ClinVar control (coverage universe)
CLINVAR_SIG = "ClinVar_CLNSIG"                  # custom-track field every ClinVar control must carry
CCRE = "ENCODE_cCRE"
REPORT_COLUMNS = ("plugin", "status", "header_ok", "controls_n", "controls_annotated", "eligible_n",
                  "annotated_n", "coverage", "message")

# stratum -> substrings of INFO/MC that put a ClinVar record in it
STRATA = {
    "missense": ("missense_variant",),
    "nonsense_frameshift": ("nonsense", "frameshift_variant"),
    "splice_region": ("splice",),
    "5_prime_UTR": ("5_prime_UTR_variant",),
    "intron": ("intron_variant",),
}
MAX_ALLELE_LEN = 50                             # longer REF/ALT (ClinVar's kb-scale deletions) are not smoke controls
_ACGT = re.compile(r"^[ACGT]+$")
_FORMAT = re.compile(r'Format:\s*([^"]+)"')
_VCF_COLS = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")


# ----------------------------------------------------------------------------- naming and sampling

def parse_region(region: str) -> tuple[str, int, int]:
    m = re.fullmatch(r"([^:\s]+):(\d+)-(\d+)", str(region).strip())
    if not m or int(m.group(2)) > int(m.group(3)):
        raise ConfigError(f"smoke.region must look like 17:43044295-43125483, got {region!r}")
    return m.group(1), int(m.group(2)), int(m.group(3))


def _ensembl(contig: str) -> str:
    """chr17 -> 17, chrM -> MT (the reference FASTA's naming)."""
    if contig == "chrM":
        return "MT"
    return contig[3:] if contig.startswith("chr") else contig


def _contig_names(contig: str) -> list:
    """The contig as given first, then the other naming convention."""
    base = _ensembl(contig)
    names = [contig]
    for alt in (base, "chrM" if base == "MT" else f"chr{base}"):
        if alt not in names:
            names.append(alt)
    return names


def _stride_pick(rows: list, n: int, seed: int) -> list:
    """`n` rows spread evenly over the (sorted) candidates: fixed stride, fixed offset from the seed."""
    if len(rows) <= n:
        return list(rows)
    stride = len(rows) // n
    return rows[seed % stride::stride][:n]


def _alts(ref: str, alt_field: str) -> list:
    """Plain ACGT (REF, ALT) pairs of a data row; anything else (symbolic, `-`, `.`, lowercase, longer than
    MAX_ALLELE_LEN) is skipped."""
    if not _ACGT.match(ref) or len(ref) > MAX_ALLELE_LEN:
        return []
    return [(ref, a) for a in alt_field.split(",") if _ACGT.match(a) and len(a) <= MAX_ALLELE_LEN]


def _tabix_lines(path, contig, start, end, env) -> list:
    """Data lines of `path` in the region; tries the contig as given, then with `chr` added/removed."""
    last = None
    for name in _contig_names(contig):
        try:
            out = _run_ok(["tabix", path, f"{name}:{start}-{end}"], env)
        except RuntimeError as e:
            last = e
            continue
        lines = [l for l in out.splitlines() if l and not l.startswith("#")]
        if lines:
            return lines
        last = None
    if last is not None:
        raise last
    return []


def _indexed(path: Path) -> bool:
    return path.is_file() and any(Path(f"{path}{ext}").is_file() for ext in (".tbi", ".csi"))


# ----------------------------------------------------------------------------- control selection

def _plugin_controls(rp, contig, start, end, seed, env) -> list:
    """(chrom, pos, ref, alt) controls for a data/vcf-source plugin, from every tabix-indexed file it has (each file
    of a `glob` directory counts as one)."""
    ctl = rp.spec.control
    cols = [1, 2, 4, 5] if ctl["source"] == "vcf" else ctl["cols"]
    specs = {fs.key: fs for fs in rp.spec.files}
    files = [f for key, path in rp.paths.items() for f in (specs[key].members(path) if key in specs else [path])]
    picked = []
    for path in files:
        if not _indexed(Path(path)):
            continue
        rows = set()
        for line in _tabix_lines(str(path), contig, start, end, env):
            f = line.split("\t")
            if len(f) < max(cols):
                continue
            chrom, pos, ref, alt = (f[c - 1] for c in cols)
            if not pos.isdigit() or not start <= int(pos) <= end:
                continue
            for r, a in _alts(ref, alt):
                rows.add((_ensembl(chrom), int(pos), r, a))
        picked += _stride_pick(sorted(rows), CONTROLS_PER_FILE, seed)
    return picked


def _info(field: str) -> dict:
    return dict(kv.split("=", 1) for kv in field.split(";") if "=" in kv)


def _clinvar_candidates(path, contig, start, end, env) -> list:
    """Sorted (chrom, pos, ref, alt, MC) of the ClinVar records in the region, whatever their significance or review
    status (the strata only need the consequence)."""
    rows = set()
    for line in _tabix_lines(str(path), contig, start, end, env):
        f = line.split("\t")
        if len(f) < 8 or not f[1].isdigit() or not start <= int(f[1]) <= end:
            continue
        info = _info(f[7])
        for ref, alt in _alts(f[3], f[4]):
            rows.add((_ensembl(f[0]), int(f[1]), ref, alt, info.get("MC", "")))
    return sorted(rows)


def _clinvar_file(custom) -> str | None:
    """The ClinVar custom track's file if the track is enabled and has one, else None."""
    cv = (custom or {}).get("ClinVar") or {}
    return cv.get("file") if cv.get("enabled") else None


def collect_controls(cfg, plugins, region=DEFAULT_REGION, seed=DEFAULT_SEED) -> dict:
    """{(chrom, pos, ref, alt): set(tags)} before normalisation (chrom in FASTA naming)."""
    contig, start, end = parse_region(region)
    env = cfg["tools"]["env"]
    controls: dict = {}
    for rp in plugins:
        if rp.spec.control.get("source") in ("data", "vcf"):
            for rec in _plugin_controls(rp, contig, start, end, seed, env):
                controls.setdefault(rec, set()).add(f"ctl_{rp.spec.name}")
    clinvar = _clinvar_file(cfg.get("custom"))
    if clinvar:
        cands = _clinvar_candidates(clinvar, contig, start, end, env)
        for stratum, terms in STRATA.items():
            hits = [c for c in cands if any(t in c[4] for t in terms)]
            for rec in _stride_pick(hits, CLINVAR_PER_STRATUM, seed):
                controls.setdefault(rec[:4], set()).add(f"cv_{stratum}")
    return controls


# ----------------------------------------------------------------------------- the smoke VCF

def smoke_dir(cfg) -> Path:
    return Path(cfg["output"]["dir"]) / "smoke"


def _fai_lengths(fasta) -> dict:
    fai = Path(f"{fasta}.fai")
    if not fai.is_file():
        raise ConfigError(f"input.reference_fasta has no .fai index: {fai}")
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in fai.read_text().splitlines() if l}


def _raw_vcf(cfg, controls, contigs: dict) -> str:
    lines = ["##fileformat=VCFv4.2", *(f"##contig=<ID={c},length={n}>" for c, n in contigs.items()),
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             "\t".join([*_VCF_COLS, cfg["sample"]])]
    for (chrom, pos, ref, alt), tags in sorted(controls.items(), key=lambda kv: (list(contigs).index(kv[0][0]),
                                                                                  *kv[0][1:])):
        lines.append("\t".join([chrom, str(pos), ";".join(sorted(tags)), ref, alt, ".", ".", ".", "GT", "0/1"]))
    return "\n".join(lines) + "\n"


def _normalise(cfg, raw: Path, norm: Path):
    _run_ok(["bcftools", "norm", "-f", cfg["input"]["reference_fasta"], "-m", "-any", "-c", "w", "-Ov",
             "-o", norm, raw], cfg["tools"]["env"])


def _canonical(norm_text: str, order: dict) -> str:
    """Drop bcftools' command/version header lines, sort, and merge records that normalised to the same allele."""
    head, merged = [], {}
    for line in norm_text.splitlines():
        if line.startswith("#"):
            if not line.startswith("##bcftools_"):
                head.append(line)
            continue
        f = line.split("\t")
        key = (f[0], int(f[1]), f[3], f[4])
        if key in merged:
            merged[key][2] = ";".join(sorted(set(merged[key][2].split(";")) | set(f[2].split(";"))))
        else:
            merged[key] = f
    body = ["\t".join(f) for _k, f in sorted(merged.items(), key=lambda kv: (order[kv[0][0]], *kv[0][1:]))]
    return "\n".join(head + body) + "\n"


def _build(cfg, plugins, region, seed) -> Path:
    contig, _s, _e = parse_region(region)
    lengths = _fai_lengths(cfg["input"]["reference_fasta"])
    contig = _ensembl(contig)
    if contig not in lengths:
        raise ConfigError(f"smoke.region contig {contig!r} is not in the reference FASTA index")
    controls = collect_controls(cfg, plugins, region, seed)
    contigs = {contig: lengths[contig]}
    out_dir = smoke_dir(cfg)
    out_dir.mkdir(parents=True, exist_ok=True)
    raw, norm = out_dir / "smoke.controls.raw.vcf", out_dir / "smoke.controls.norm.vcf"
    raw.write_text(_raw_vcf(cfg, controls, contigs))
    _normalise(cfg, raw, norm)
    final = out_dir / "smoke.controls.vcf"
    final.write_text(_canonical(norm.read_text(), {contig: 0}))
    raw.unlink()
    norm.unlink()
    env = cfg["tools"]["env"]
    _run_ok(["bgzip", "-f", final], env)
    _run_ok(["tabix", "-f", "-p", "vcf", f"{final}.gz"], env)
    return Path(f"{final}.gz")


def build_smoke_vcf(cfg, plugins, region=None, seed=DEFAULT_SEED) -> Path:
    """Write `{output.dir}/smoke/smoke.controls.vcf.gz` (bgzipped, tabix-indexed, byte-stable for a given input).

    `region` defaults to cfg['smoke']['region'], then to BRCA1."""
    region = region or (cfg.get("smoke") or {}).get("region") or DEFAULT_REGION
    return _build(cfg, plugins, region, seed)


def ref_mismatches(cfg, vcf) -> set:
    """Keys (jobs._key) of the smoke VCF records whose REF is not the reference FASTA sequence at that position."""
    recs = [l.split("\t") for l in _read_lines(vcf) if l and not l.startswith("#")]
    if not recs:
        return set()
    seqs = fasta_seqs(cfg["input"]["reference_fasta"], [f"{f[0]}:{f[1]}-{int(f[1]) + len(f[3]) - 1}" for f in recs],
                      cfg["tools"]["env"])
    return {_key(f[0], f[1], f[3], f[4]) for f, seq in zip(recs, seqs) if seq != f[3].upper()}


def _stratum_matches(mc: str, stratum: str) -> bool:
    terms = (stratum, *STRATA[stratum])
    return any(mc in t or t in mc for t in terms)


def control_keys_from_vcf(vcf, plugins) -> dict:
    """{plugin: {key}} plus {"ClinVar": {key}} read back from the ID tags of a smoke VCF.

    Key format is jobs._key: CHROM:POS:REF>ALT. A clinvar-source plugin's controls are the ClinVar controls of
    the strata matching its `control.mc`."""
    tags: dict = {}
    for line in _read_lines(vcf):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        tags[_key(f[0], f[1], f[3], f[4])] = set(f[2].split(";"))
    keys = {CLINVAR_KEY: {k for k, t in tags.items() if any(x.startswith("cv_") for x in t)}}
    for rp in plugins:
        name, ctl = rp.spec.name, rp.spec.control
        if ctl.get("source") == "clinvar":
            strata = [s for s in STRATA if _stratum_matches(ctl["mc"], s)]
            keys[name] = {k for k, t in tags.items() if any(f"cv_{s}" in t for s in strata)}
        else:
            keys[name] = {k for k, t in tags.items() if f"ctl_{name}" in t}
    return keys


# ----------------------------------------------------------------------------- report

def _read_lines(path) -> list:
    path = Path(path)
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    opener = gzip.open if gz else open
    with opener(path, "rt") as fh:
        return fh.read().splitlines()


def _parse_vep_vcf(vcf):
    """(CSQ field names, {key: [ {field: value} per CSQ entry ]}); empty and `-` values are dropped."""
    fmt, records = None, {}
    for line in _read_lines(vcf):
        if line.startswith("##INFO=<ID=CSQ,"):
            m = _FORMAT.search(line)
            fmt = m.group(1).strip().split("|") if m else None
            continue
        if line.startswith("#") or not line:
            continue
        f = line.split("\t")
        csq = next((kv[4:] for kv in f[7].split(";") if kv.startswith("CSQ=")), "")
        entries = [{k: v for k, v in zip(fmt or (), e.split("|")) if v not in ("", "-")}
                   for e in csq.split(",") if e]
        records.setdefault(_key(f[0], f[1], f[3], f[4]), []).extend(entries)
    return fmt, records


def _has_value(entries, fields) -> bool:
    return any(e.get(f) for e in entries for f in fields)


def _eligible(entries, eligible_set) -> list:
    if not eligible_set:
        return list(entries)
    return [e for e in entries if eligible_set & set(e.get("Consequence", "").split("&"))]


def _row(plugin, status, header_ok, controls_n, controls_annotated, message, eligible_n=None, annotated_n=None):
    return {"plugin": plugin, "status": status, "header_ok": header_ok, "controls_n": controls_n,
            "controls_annotated": controls_annotated, "eligible_n": eligible_n, "annotated_n": annotated_n,
            "coverage": round(annotated_n / eligible_n, 3) if eligible_n else None, "message": message}


def _track_rows(records, present, control_keys, custom) -> list:
    """`custom:ClinVar` (every ClinVar control must carry ClinVar_CLNSIG: a free positive control of the track) and
    `custom:ENCODE_cCRE` (info: smoke records overlapping a cCRE), for the tracks that are enabled."""
    rows = []
    if _clinvar_file(custom):
        keys = sorted(control_keys.get(CLINVAR_KEY, ()))
        n = sum(1 for k in keys if _has_value(records.get(k, ()), [CLINVAR_SIG]))
        ok = CLINVAR_SIG in present
        if not ok:
            status, msg = "FAIL", f"{CLINVAR_SIG} is not in the CSQ header (custom.ClinVar.fields must include CLNSIG)"
        elif not keys:
            status, msg = "FAIL", "no ClinVar record in the smoke region to test the track on: choose another smoke.region"
        elif not n:
            status, msg = "FAIL", f"none of the {len(keys)} ClinVar controls got a {CLINVAR_SIG} value"
        else:
            status, msg = "OK", f"{n}/{len(keys)} ClinVar controls carry {CLINVAR_SIG}"
        rows.append(_row("custom:ClinVar", status, ok, len(keys), n, msg))
    if ((custom or {}).get(CCRE) or {}).get("enabled"):
        n = sum(1 for entries in records.values() if _has_value(entries, [CCRE]))
        ok = CCRE in present
        rows.append(_row(f"custom:{CCRE}", "info" if ok else "FAIL", ok, len(records), n,
                         f"{n}/{len(records)} smoke records overlap an ENCODE cCRE" if ok
                         else f"{CCRE} is not in the CSQ header"))
    return rows


def smoke_report(vep_vcf, plugins, control_keys: dict, custom=None, ref_mismatch=frozenset()) -> list:
    """One row per plugin, then one per enabled custom track (`custom` = cfg['custom']); see REPORT_COLUMNS.

    FAIL: an output field is missing from the CSQ header; a data/vcf control source with no control records
    (the locus has no data for the plugin); a plugin that annotates none of its controls; more than MISMATCH_LIMIT of
    a data/vcf plugin's controls in `ref_mismatch` (REF not the FASTA base: possible wrong assembly in plugin data);
    a clinvar-source plugin while the ClinVar track is off (`custom` given); see `_track_rows` for the tracks.
    OK: header present and (header-only control, or >=1 control annotated). `unchecked`: no control defined, or a
    ClinVar stratum without any control in the region. Coverage is annotated/eligible over the ClinVar controls
    (control_keys['ClinVar']; every record of the VCF when that entry is absent) whose consequence is eligible."""
    fmt, records = _parse_vep_vcf(vep_vcf)
    present = set(fmt or ())
    universe = control_keys.get(CLINVAR_KEY)
    universe = list(records) if universe is None else [k for k in sorted(universe) if k in records]
    rows = []
    for rp in plugins:
        spec = rp.spec
        fields = [f for f, _cls in _plugin_fields(rp)]
        missing = [f for f in fields if f not in present]
        header_ok = not missing
        ctl_keys = sorted(control_keys.get(spec.name, ()))
        annotated = sum(1 for k in ctl_keys if _has_value(records.get(k, ()), fields))
        eligible_set = set(spec.eligible_consequences)
        elig = [k for k in universe if _eligible(records[k], eligible_set)]
        ann = sum(1 for k in elig if _has_value(_eligible(records[k], eligible_set), fields))
        source = spec.control.get("source")
        status, msgs = "OK", []
        if missing:
            status = "FAIL"
            msgs.append(f"field(s) {missing} missing from the CSQ header")
        if not source:
            status = status if status == "FAIL" else "unchecked"
            msgs.append("no control defined in the registry")
        elif source in ("data", "vcf"):
            if not ctl_keys:
                status = "FAIL"
                msgs.append("the control locus has no data for this plugin: choose another smoke.region "
                            "or check the plugin data file")
            elif not annotated:
                status = "FAIL"
                msgs.append(f"none of the {len(ctl_keys)} control records from the plugin's own data was annotated")
            else:
                msgs.append(f"{annotated}/{len(ctl_keys)} controls from the plugin's own data annotated")
            bad = sum(1 for k in ctl_keys if k in ref_mismatch)
            if ctl_keys and bad / len(ctl_keys) > MISMATCH_LIMIT:
                status = "FAIL"
                msgs.append(f"{bad}/{len(ctl_keys)} controls have a REF that is not the reference FASTA base: "
                            f"possible wrong assembly in plugin data")
        elif source == "clinvar":
            if custom is not None and not _clinvar_file(custom):
                status = "FAIL"
                msgs.append("the ClinVar custom track is disabled or has no file, so this plugin's ClinVar controls "
                            "cannot be built: enable custom.ClinVar")
            elif not ctl_keys:
                status = status if status == "FAIL" else "unchecked"
                msgs.append(f"no ClinVar control matching {spec.control.get('mc')!r} in the smoke region")
            elif not annotated:
                status = "FAIL"
                msgs.append(f"none of the {len(ctl_keys)} ClinVar {spec.control.get('mc')} controls was annotated")
            else:
                msgs.append(f"{annotated}/{len(ctl_keys)} ClinVar {spec.control.get('mc')} controls annotated")
        elif source == "header" and header_ok:
            msgs.append("header present")
        rows.append(_row(spec.name, status, header_ok, len(ctl_keys), annotated, "; ".join(msgs), len(elig), ann))
    return rows + _track_rows(records, present, control_keys, custom)


def write_report(rows, path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    def cell(v):
        return "" if v is None else ("yes" if v is True else "no" if v is False else str(v))
    path.write_text("\t".join(REPORT_COLUMNS) + "\n"
                    + "".join("\t".join(cell(r[c]) for c in REPORT_COLUMNS) + "\n" for r in rows))
    return path


# ----------------------------------------------------------------------------- the run

def smoke_cfg(cfg, smoke_vcf) -> dict:
    """A copy of cfg that runs the normal pipeline on the smoke VCF under `{output.dir}/smoke` (no SV run)."""
    s = copy.deepcopy(cfg)
    s["input"] = {**s["input"], "vcf": str(smoke_vcf), "normalize": False}
    s["output"] = {**s["output"], "dir": str(smoke_dir(cfg))}
    s["sv"] = {"enabled": False}
    return s


def run_smoke(cfg, plugins, region=None, seed=DEFAULT_SEED, parallel=1) -> list:
    """Build the controls, run them through the normal pipeline, write and return `smoke/smoke.tsv` rows."""
    vcf = build_smoke_vcf(cfg, plugins, region, seed)
    scfg = smoke_cfg(cfg, vcf)
    final = run_local(render_run(scfg, plugins), parallel)
    rows = smoke_report(final, plugins, control_keys_from_vcf(vcf, plugins), custom=cfg.get("custom") or {},
                        ref_mismatch=ref_mismatches(cfg, vcf))
    write_report(rows, smoke_dir(cfg) / "smoke.tsv")
    return rows
