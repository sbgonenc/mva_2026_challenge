"""Downloads: atomic, resumable, size- and checksum-verified fetches of plugin data, the VEP cache and custom tracks."""
from __future__ import annotations
import hashlib
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from .config import ConfigError
from .registry import FileSpec, fill_prep

CA = "/etc/ssl/certs/ca-certificates.crt"
GIB = 2 ** 30

ENSEMBL_FTP = "https://ftp.ensembl.org/pub"
# Compressed size of the indexed cache tarball; for any other release/assembly the size comes from HEAD.
KNOWN_CACHE_SIZES = {(116, "GRCh38"): 27644657162}
CACHE_FREE_FACTOR = 2.2          # tarball + extracted tree must fit at the same time

CCRE_OUT = "GRCh38-cCREs.class.bed.gz"
CCRE_CMD = ('awk \'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,$6":"$5}}\' {src} | LC_ALL=C sort -k1,1 -k2,2n '
            '| bgzip -f > {part} && mv -f {part} {out} && tabix -f -p bed {out}')


# ----------------------------------------------------------------------------- primitives

def https_env() -> dict:
    """os.environ with a usable CA bundle: the login env points SSL_CERT_FILE at a missing file (curl error 77)."""
    env = dict(os.environ)
    if not Path(env.get("SSL_CERT_FILE", "/nonexistent")).is_file():
        env["SSL_CERT_FILE"] = CA
    env["CURL_CA_BUNDLE"] = env["SSL_CERT_FILE"]
    return env


def fetch(url, dest, size=None, pget=None, big=2 ** 31, verify=None):
    """Download url to dest through dest.part, verify the size when known and run `verify(part)` (checksums), then
    rename (atomic, resumable). The final name therefore only ever holds a verified file; a part that fails
    `verify` is deleted."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    if size is not None and part.exists() and part.stat().st_size > size:
        part.unlink()                                  # cannot be a prefix of the wanted file
    if size is None and part.exists():
        part.unlink()                                  # no size to tell a prefix from a stale file; curl -C - would
                                                       # append to it, or fail (error 33) where ranges are unsupported
    if size is not None and part.exists() and part.stat().st_size == size:
        pass                                           # an earlier run finished the transfer but not the rename
    elif size and size > big and pget and Path(pget).exists() and url.startswith("http"):
        subprocess.run(["bash", str(pget), url, str(part), "8", "512"], check=True, env=https_env())
    else:
        subprocess.run(["curl", "-fL", "--retry", "5", "--retry-delay", "5", "-C", "-", "-o", str(part), url],
                       check=True, env=https_env())
    if size is not None and part.stat().st_size != size:
        raise RuntimeError(f"size mismatch for {url}: {part.stat().st_size} != {size} (partial kept: {part})")
    if verify is not None:
        try:
            verify(part)
        except RuntimeError:
            part.unlink()                              # complete but wrong: resuming it would never fix it
            raise
    part.rename(dest)


def content_length(url) -> int:
    """Content-Length of url from a HEAD request (the last one after redirects)."""
    out = subprocess.run(["curl", "-sIL", "--max-time", "60", url], capture_output=True, text=True,
                         env=https_env()).stdout
    lengths = [int(l.split(":", 1)[1]) for l in out.splitlines() if l.lower().startswith("content-length:")]
    if not lengths or lengths[-1] <= 0:
        raise RuntimeError(f"could not determine the size of {url}")
    return lengths[-1]


def read_text(url) -> str:
    return subprocess.run(["curl", "-fsSL", "--retry", "3", url], capture_output=True, text=True,
                          check=True, env=https_env()).stdout


def bsd_sum(path):
    """BSD checksum and 1 KiB block count, as printed by `sum -r` (tests only; verification shells out)."""
    s = n = 0
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            n += len(chunk)
            for b in chunk:
                s = ((s >> 1) + ((s & 1) << 15) + b) & 0xFFFF
    return s, (n + 1023) // 1024


def verify_checksums(path, checksums_text, name=None):
    """Compare `sum -r path` with the `<sum> <blocks> <name>` line in an Ensembl CHECKSUMS text; `name` defaults to
    path.name (pass the final name when checking a `.part`)."""
    path = Path(path)
    name = name or path.name
    want = next((l.split()[:2] for l in checksums_text.splitlines() if l.split()[-1:] == [name]), None)
    if want is None:
        raise RuntimeError(f"checksum entry for {name} not found")
    got = subprocess.run(["sum", "-r", str(path)], capture_output=True, text=True, check=True).stdout.split()[:2]
    if [int(x) for x in got] != [int(x) for x in want]:
        raise RuntimeError(f"checksum mismatch for {name}: {got} != {want}")


def file_info(path):
    """(bytes, md5) of a file."""
    md5 = hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 22):
            md5.update(chunk)
    return Path(path).stat().st_size, md5.hexdigest()


def append_manifest(root, name, url, path, info=None):
    """Append `name url bytes md5 date` to <root>/manifest.tsv (header written on creation)."""
    manifest = Path(root) / "manifest.tsv"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    nbytes, md5 = info or file_info(path)
    new = not manifest.exists()
    with open(manifest, "a") as fh:
        if new:
            fh.write("name\turl\tbytes\tmd5\tdate\n")
        fh.write(f"{name}\t{url}\t{nbytes}\t{md5}\t{datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}\n")


def run_prep(cmds, path, env, empty_indel=None):
    """Run each prep command in `micromamba run -n env bash -euo pipefail -c`, {path}/{dir}/{empty_indel} filled in."""
    for cmd in cmds:
        filled = fill_prep(cmd, Path(path), empty_indel=empty_indel)
        r = subprocess.run(["micromamba", "run", "-n", env, "bash", "-euo", "pipefail", "-c", filled])
        if r.returncode:
            raise RuntimeError(f"prep command failed (exit {r.returncode}) for {path}: {filled}")


def validate_file(path, rules: dict):
    """Content checks for downloads that can fail with HTTP 200 (error pages).

    Rules: min_lines; header_contains (alias first_line_contains): every term in the first line, case-insensitive;
    not_gzip (alias plain_text): the file must not start with the gzip magic bytes.
    """
    path = Path(path)
    known = {"min_lines", "header_contains", "first_line_contains", "not_gzip", "plain_text"}
    if set(rules) - known:
        raise RuntimeError(f"unknown validate rule(s) {sorted(set(rules) - known)} for {path}")
    if rules.get("not_gzip") or rules.get("plain_text"):
        with open(path, "rb") as fh:
            if fh.read(2) == b"\x1f\x8b":
                raise RuntimeError(f"{path} is gzip-compressed but a plain text file was expected (not_gzip)")
    terms = rules.get("header_contains") or rules.get("first_line_contains") or []
    with open(path, "rb") as fh:
        first = fh.readline().decode("utf-8", "replace").lower()
        missing = [t for t in terms if t.lower() not in first]
        if missing:
            raise RuntimeError(f"first line of {path} lacks {missing} (header_contains): {first[:120]!r}")
        if "min_lines" in rules:
            n = 1 + sum(chunk.count(b"\n") for chunk in iter(lambda: fh.read(1 << 22), b""))
            if n < rules["min_lines"]:
                raise RuntimeError(f"{path} has {n} lines, fewer than min_lines={rules['min_lines']}")


# ----------------------------------------------------------------------------- planning

@dataclass
class DownloadItem:
    name: str
    final: Path
    action: str                      # fetch | skip | missing
    nbytes: int | None = None        # bytes still to transfer; None when a size is unknown
    reason: str = ""
    manual: str | None = None
    plugin: str = ""
    spec: FileSpec | None = None


def _nonempty(p: Path) -> bool:
    return p.is_file() and p.stat().st_size > 0


def _sidecar(final: Path, index) -> Path | None:
    return Path(f"{final}.{index}") if index else None


def _complete(final: Path, fs) -> bool:
    """Every data file behind `final` (one file, or the `glob` matches of a directory) is non-empty and indexed."""
    files = fs.members(final)
    return bool(files) and all(_nonempty(f) and (not fs.index or _sidecar(f, fs.index).exists()) for f in files)


def _user_path_problem(final: Path, fs) -> str:
    if not fs.glob:
        return ("user-supplied file does not exist" if not _nonempty(final)
                else f"user-supplied file has no index {_sidecar(final, fs.index).name} (run tabix)")
    files = fs.members(final)
    if not files:
        return f"user-supplied directory has no files matching {fs.glob}"
    bad = [f.name for f in files if not (_nonempty(f) and (not fs.index or _sidecar(f, fs.index).exists()))]
    return f"{len(bad)} of {len(files)} files are empty or have no .{fs.index} index (run tabix): {', '.join(bad[:5])}"


def plan_downloads(cfg, plugins, include_skipped=False) -> list:
    """One DownloadItem per data file of the enabled plugins that still has to be fetched (or is missing)."""
    data_dir = Path(cfg["data"]["plugin_data_dir"])
    items = []
    for rp in plugins:
        spec = rp.spec
        for fs in spec.files:
            final = rp.paths.get(fs.key)
            if final is None:                            # `args:` escape hatch: the user owns every file
                continue
            name = spec.name if len(spec.files) == 1 else f"{spec.name}.{fs.key}"
            item = DownloadItem(name, final, "fetch", plugin=spec.name, spec=fs, manual=fs.manual)
            if _complete(final, fs):
                item.action = "skip"
            elif final != data_dir / fs.dest:           # user-set path: verified, never written
                item.action = "missing"
                item.reason = _user_path_problem(final, fs)
            elif not fs.downloads and not fs.prep:
                item.action = "missing"
                item.reason = "no download recipe in the registry"
            else:
                todo = [d for d in fs.downloads if not _nonempty(data_dir / d.as_)]
                item.nbytes = None if any(d.size is None for d in todo) else sum(d.size for d in todo)
            if item.action != "skip" or include_skipped:
                items.append(item)
    return items


def _free_bytes(path: Path) -> int:
    p = Path(path)
    while not p.exists() and p != p.parent:
        p = p.parent
    return shutil.disk_usage(p).free


def print_plan(cfg, items, out=None):
    """Table of name, bytes, destination and action, then total bytes still to fetch and free disk space."""
    out = out or sys.stdout
    rows = [("name", "bytes", "action", "destination")]
    total, unknown = 0, 0
    for it in items:
        if it.action == "fetch":
            total += it.nbytes or 0
            unknown += it.nbytes is None
        shown = "-" if it.action == "skip" else ("?" if it.nbytes is None else str(it.nbytes))
        action = {"fetch": "fetch", "skip": "skip", "missing": "MISSING"}[it.action]
        rows.append((it.name, shown, action, str(it.final)))
    widths = [max(len(r[i]) for r in rows) for i in range(3)]
    for r in rows:
        print("  ".join(c.ljust(w) for c, w in zip(r[:3], widths)) + "  " + r[3], file=out)
    for it in items:
        if it.action == "missing":
            print(f"\n{it.name}: {it.reason}: {it.final}", file=out)
        if it.action != "skip" and it.manual:
            print(f"{it.name}: " + " ".join(it.manual.split()), file=out)
    note = f" (+ {unknown} item(s) of unknown size)" if unknown else ""
    print(f"\ntotal to fetch: {total} bytes ({total / GIB:.2f} GiB){note}", file=out)
    print(f"free disk space under {cfg['data']['root']}: {_free_bytes(Path(cfg['data']['root']))} bytes "
          f"({_free_bytes(Path(cfg['data']['root'])) / GIB:.1f} GiB)", file=out)


# ----------------------------------------------------------------------------- execution

def _fetch_item(cfg, item: DownloadItem):
    fs, data_dir = item.spec, Path(cfg["data"]["plugin_data_dir"])
    pget, root = cfg["tools"].get("pget"), cfg["data"]["root"]
    checksums = read_text(fs.checksums_url) if fs.checksums_url else None
    rows = []                                            # manifest rows, recorded once the item is complete
    for d in fs.downloads:
        dest = data_dir / d.as_
        if _nonempty(dest):
            continue                                     # left by an interrupted earlier run (verified before rename)
        verify = None if checksums is None else (lambda part, name=dest.name: verify_checksums(part, checksums, name))
        fetch(d.url, dest, d.size, pget, verify=verify)
        rows.append((d.url, dest, file_info(dest)))
    idx = _sidecar(item.final, fs.index)
    if fs.index_url and idx is not None:
        fetch(fs.index_url, idx)
        rows.append((fs.index_url, idx, file_info(idx)))
    if fs.prep:
        run_prep(fs.prep, item.final, cfg["tools"]["env"])
    if not _nonempty(item.final):
        raise RuntimeError(f"{item.name}: {item.final} is missing or empty after download/prep")
    if idx is not None and not idx.exists():
        raise RuntimeError(f"{item.name}: index {idx} was not created (registry needs index_url or a prep command)")
    if fs.validate:
        try:
            validate_file(item.final, fs.validate)
        except RuntimeError:
            item.final.unlink()
            raise
    for url, path, info in rows:
        append_manifest(root, item.name, url, path, info)


def download_all(cfg, plugins, dry_run=False) -> list:
    """Fetch every missing plugin data file; returns the items acted on. dry_run prints the plan and stops."""
    items = plan_downloads(cfg, plugins, include_skipped=dry_run)
    if dry_run:
        print_plan(cfg, items)
        return items
    missing = [i for i in items if i.action == "missing"]
    if missing:
        raise ConfigError("data files that cannot be downloaded: " + "; ".join(
            f"{i.name}: {i.reason}: {i.final}" + (f" -- {' '.join(i.manual.split())}" if i.manual else "")
            for i in missing))
    for n, item in enumerate(items, 1):
        print(f"[{n}/{len(items)}] {item.name} -> {item.final}", flush=True)
        _fetch_item(cfg, item)
    return items


# ----------------------------------------------------------------------------- VEP cache

def _cache_tag(cfg):
    return f"{cfg['cache_version']}_{cfg['assembly']}"


def _cache_base(cfg):
    return f"{ENSEMBL_FTP}/release-{cfg['cache_version']}/variation/indexed_vep_cache"


def cache_name(cfg):
    return f"homo_sapiens_vep_{_cache_tag(cfg)}.tar.gz"


def cache_url(cfg):
    return f"{_cache_base(cfg)}/{cache_name(cfg)}"


def cache_checksums_url(cfg):
    return f"{_cache_base(cfg)}/CHECKSUMS"


def cache_size(cfg):
    """Tarball bytes: the known value for this release, else None (download_cache then asks the server)."""
    return KNOWN_CACHE_SIZES.get((cfg["cache_version"], cfg["assembly"]))


def _cache_info(cfg) -> Path:
    return Path(cfg["data"]["cache_dir"]) / "homo_sapiens" / _cache_tag(cfg) / "info.txt"


def cache_ready(cfg) -> bool:
    info = _cache_info(cfg)
    return info.is_file() and "var_type\ttabix" in info.read_text()


def cache_item(cfg) -> DownloadItem:
    return DownloadItem("VEP cache", Path(cfg["data"]["cache_dir"]), "skip" if cache_ready(cfg) else "fetch",
                        nbytes=cache_size(cfg))


def download_cache(cfg, url=None, checksums_url=None, size=None) -> Path:
    """Fetch, verify and unpack the tabix-indexed VEP cache; returns <cache_dir>/homo_sapiens/<ver>_<assembly>.

    url/checksums_url/size default to the Ensembl release derived from cache_version/assembly; the size comes from
    KNOWN_CACHE_SIZES, else from a HEAD request. They are parameters so tests can use a synthetic file:// tarball.
    """
    cache_dir = Path(cfg["data"]["cache_dir"])
    if cache_ready(cfg):
        return _cache_info(cfg).parent
    url = url or cache_url(cfg)
    checksums_url = checksums_url or cache_checksums_url(cfg)
    size = size or cache_size(cfg) or content_length(url)
    cache_dir.mkdir(parents=True, exist_ok=True)
    tarball = cache_dir / ".dl" / cache_name(cfg)
    reuse = tarball.is_file() and tarball.stat().st_size == size    # complete tarball from an interrupted run
    need = CACHE_FREE_FACTOR * size - (size if reuse else 0)
    free = _free_bytes(cache_dir)
    if free < need:
        raise RuntimeError(f"not enough free disk space under {cache_dir}: {free} bytes free, need {int(need)}")
    checksums = read_text(checksums_url)
    if reuse:                       # left by an interrupted extraction, or by an older run that renamed before verifying
        try:
            verify_checksums(tarball, checksums)
        except RuntimeError:
            tarball.unlink()                                        # else the next run would reuse a corrupt file
            raise
    else:                                                           # verified as .part, before the rename
        fetch(url, tarball, size, cfg["tools"].get("pget"),
              verify=lambda part: verify_checksums(part, checksums, tarball.name))
    # Extract into a staging dir and move into place only when complete: an interrupted tar must never leave a
    # partial homo_sapiens/<ver>_<assembly> that cache_ready() would accept.
    stage = cache_dir / ".extract"
    shutil.rmtree(stage, ignore_errors=True)
    stage.mkdir()
    subprocess.run(["tar", "-xzf", str(tarball), "-C", str(stage)], check=True)
    staged = stage / "homo_sapiens" / _cache_tag(cfg)
    info = staged / "info.txt"
    if not (info.is_file() and "var_type\ttabix" in info.read_text()):
        raise RuntimeError(f"{info} is missing or does not contain 'var_type\\ttabix': "
                           f"the extracted cache is not the tabix-indexed one")
    target = _cache_info(cfg).parent
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.rmtree(target, ignore_errors=True)                       # leftover of a pre-staging partial extraction
    staged.rename(target)
    shutil.rmtree(stage, ignore_errors=True)
    append_manifest(cfg["data"]["root"], "VEP cache", url, tarball)
    tarball.unlink()
    return _cache_info(cfg).parent


# ----------------------------------------------------------------------------- custom tracks

def _require(path, what):
    if not path or not Path(path).is_file():
        raise ConfigError(f"{what}: file not found: {path}")
    return Path(path)


def _bash_in_env(cfg, script):
    r = subprocess.run(["micromamba", "run", "-n", cfg["tools"]["env"], "bash", "-euo", "pipefail", "-c", script])
    if r.returncode:
        raise RuntimeError(f"custom track command failed (exit {r.returncode}): {script}")


def _has_index(path: Path) -> bool:
    return any(Path(f"{path}.{e}").exists() for e in ("tbi", "csi"))


def prepare_custom(cfg) -> dict:
    """Prepare the enabled custom tracks under data.custom_data_dir; returns {name: path of the indexed file}.

    ENCODE_cCRE: `source` BED -> bgzipped, sorted chrom/start/class:accession BED with a tabix index.
    Any other entry needs `file`; a VCF without an index is symlinked into custom_data_dir and the symlink is
    indexed, so nothing is ever written next to the user's file.
    """
    custom_dir = Path(cfg["data"]["custom_data_dir"])
    out = {}
    for name, entry in (cfg.get("custom") or {}).items():
        if not entry.get("enabled"):
            continue
        if name == "ENCODE_cCRE":
            src = _require(entry.get("source"), "custom.ENCODE_cCRE.source")
            dest = custom_dir / CCRE_OUT
            if not (dest.exists() and Path(f"{dest}.tbi").exists() and dest.stat().st_mtime >= src.stat().st_mtime):
                custom_dir.mkdir(parents=True, exist_ok=True)
                _bash_in_env(cfg, CCRE_CMD.format(src=shlex.quote(str(src)), out=shlex.quote(str(dest)),
                                                  part=shlex.quote(f"{dest}.part")))
            out[name] = dest
            continue
        path = _require(entry.get("file"), f"custom.{name}.file")
        if entry.get("format") == "vcf" and not _has_index(path):
            custom_dir.mkdir(parents=True, exist_ok=True)
            link = custom_dir / path.name
            if (link.is_symlink() or link.exists()) and link.resolve() != path.resolve():
                link.unlink()                            # left by an earlier config: points at another file
                for ext in ("tbi", "csi"):
                    Path(f"{link}.{ext}").unlink(missing_ok=True)
            if not (link.is_symlink() or link.exists()):
                link.symlink_to(path.resolve())
            if not _has_index(link):
                _bash_in_env(cfg, f"tabix -f -p vcf {shlex.quote(str(link))}")
            path = link
        out[name] = path
    return out
