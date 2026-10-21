import gzip
import random
import subprocess
import tarfile
import zipfile
from argparse import Namespace
from pathlib import Path

import pytest
import yaml

from vep_pipeline import download as dl
from vep_pipeline.config import ConfigError
from vep_pipeline.registry import enabled_plugins, load_registry

PLUGINS_YAML = Path(__file__).resolve().parents[1] / "plugins.yaml"


def test_fetch_atomic(tmp_path):
    src = tmp_path / "s"; src.write_bytes(b"x" * 1000); dest = tmp_path / "o/d"
    dl.fetch(src.as_uri(), dest, size=1000)
    assert dest.read_bytes() == b"x" * 1000 and not (tmp_path / "o/d.part").exists()


def test_fetch_size_mismatch(tmp_path):
    src = tmp_path / "s"; src.write_bytes(b"x" * 10)
    with pytest.raises(RuntimeError, match="size"):
        dl.fetch(src.as_uri(), tmp_path / "d", size=11)
    assert not (tmp_path / "d").exists()


def test_fetch_discards_oversized_partial(tmp_path):
    src = tmp_path / "s"; src.write_bytes(b"x" * 10)
    (tmp_path / "d.part").write_bytes(b"junk" * 100)
    dl.fetch(src.as_uri(), tmp_path / "d", size=10)
    assert (tmp_path / "d").read_bytes() == b"x" * 10


def test_fetch_unknown_size_discards_stale_partial(tmp_path):
    """S4(b): without a size a leftover .part cannot be told from a resumable prefix; `curl -C -` would append to
    it (or fail with error 33 on servers without ranges, e.g. the G2P API and OSF)."""
    src = tmp_path / "s"; src.write_bytes(b"fresh content " * 50)
    (tmp_path / "d.part").write_bytes(b"stale")
    dl.fetch(src.as_uri(), tmp_path / "d")
    assert (tmp_path / "d").read_bytes() == src.read_bytes() and not (tmp_path / "d.part").exists()


def test_fetch_verifies_the_part_before_renaming(tmp_path):
    src = tmp_path / "s"; src.write_bytes(b"x" * 10)
    seen = []
    dl.fetch(src.as_uri(), tmp_path / "d", size=10, verify=lambda p: seen.append((p.name, (tmp_path / "d").exists())))
    assert seen == [("d.part", False)] and (tmp_path / "d").read_bytes() == b"x" * 10

    def bad(p):
        raise RuntimeError("checksum mismatch")
    with pytest.raises(RuntimeError, match="checksum"):
        dl.fetch(src.as_uri(), tmp_path / "e", size=10, verify=bad)
    assert not (tmp_path / "e").exists() and not (tmp_path / "e.part").exists()


def test_bsd_sum_matches_coreutils(tmp_path):
    f = tmp_path / "a"; f.write_bytes(bytes(range(256)) * 50)
    out = subprocess.run(["sum", "-r", str(f)], capture_output=True, text=True).stdout.split()
    assert dl.bsd_sum(f) == (int(out[0]), int(out[1]))


def test_verify_checksums_mismatch(tmp_path):
    f = tmp_path / "x.gz"; f.write_bytes(b"abc")
    with pytest.raises(RuntimeError, match="checksum"):
        dl.verify_checksums(f, "00001 1 x.gz\n")


def test_verify_checksums_ok_and_missing_entry(tmp_path):
    f = tmp_path / "x.gz"; f.write_bytes(b"abc")
    s, k = dl.bsd_sum(f)
    dl.verify_checksums(f, f"1 1 other.gz\n{s:05d} {k} x.gz\n")
    with pytest.raises(RuntimeError, match="not found"):
        dl.verify_checksums(f, "1 1 other.gz\n")


def test_validate_g2p_rules(tmp_path):
    f = tmp_path / "g.csv"; f.write_text("g2p id,gene symbol,allelic requirement\n" + "r\n" * 5)
    with pytest.raises(RuntimeError, match="min_lines"):
        dl.validate_file(f, {"min_lines": 1000, "header_contains": ["allelic requirement", "gene symbol"], "not_gzip": True})


def test_validate_registry_rule_names(tmp_path):
    """plugins.yaml spells the rules plain_text / first_line_contains."""
    f = tmp_path / "g.csv"; f.write_text("G2P ID,Gene Symbol,Allelic Requirement\n" + "r\n" * 5)
    dl.validate_file(f, {"plain_text": True, "first_line_contains": ["allelic requirement"], "min_lines": 5})
    with pytest.raises(RuntimeError, match="first line"):
        dl.validate_file(f, {"first_line_contains": ["no such column"]})
    with pytest.raises(RuntimeError, match="unknown"):
        dl.validate_file(f, {"bogus": 1})


def test_validate_rejects_gzip_and_html(tmp_path):
    g = tmp_path / "g"; g.write_bytes(b"\x1f\x8b\x08\x00rest")
    with pytest.raises(RuntimeError, match="gzip"):
        dl.validate_file(g, {"not_gzip": True})
    h = tmp_path / "h"; h.write_text("<html>error</html>\n")
    with pytest.raises(RuntimeError, match="first line"):
        dl.validate_file(h, {"header_contains": ["gene symbol"]})


def test_plan_skips_complete_items(tmp_path, minimal_cfg, mini_registry):
    minimal_cfg["plugins"] = {"Toy": {"enabled": True}}
    plugins = enabled_plugins(minimal_cfg, mini_registry)
    assert len(dl.plan_downloads(minimal_cfg, plugins)) == 1
    final = Path(minimal_cfg["data"]["plugin_data_dir"]) / "Toy/toy.tsv.gz"
    final.parent.mkdir(parents=True); final.write_bytes(b"z"); Path(str(final) + ".tbi").write_bytes(b"i")
    assert dl.plan_downloads(minimal_cfg, plugins) == []
    (item,) = dl.plan_downloads(minimal_cfg, plugins, include_skipped=True)
    assert item.action == "skip" and item.final == final


def test_plan_missing_index_is_not_complete(tmp_path, minimal_cfg, mini_registry):
    minimal_cfg["plugins"] = {"Toy": {"enabled": True}}
    plugins = enabled_plugins(minimal_cfg, mini_registry)
    final = Path(minimal_cfg["data"]["plugin_data_dir"]) / "Toy/toy.tsv.gz"
    final.parent.mkdir(parents=True); final.write_bytes(b"z")
    (item,) = dl.plan_downloads(minimal_cfg, plugins)
    assert item.action == "fetch"


def test_https_env_fix(monkeypatch):
    monkeypatch.setenv("SSL_CERT_FILE", "/nonexistent.crt")
    env = dl.https_env()
    assert env["SSL_CERT_FILE"] == "/etc/ssl/certs/ca-certificates.crt" == env["CURL_CA_BUNDLE"]


def test_prepare_custom_ccre(tmp_path, minimal_cfg):
    bed = tmp_path / "c.bed"; bed.write_text("chr2\t5\t9\tD1\tE1\tPLS\nchr1\t1\t4\tD2\tE2\tpELS\n")
    minimal_cfg["custom"] = {"ENCODE_cCRE": {"enabled": True, "source": str(bed), "format": "bed", "type": "overlap"}}
    out = dl.prepare_custom(minimal_cfg)["ENCODE_cCRE"]
    first = subprocess.run(["bash", "-c", f"zcat {out} | head -1"], capture_output=True, text=True).stdout
    assert first.rstrip("\n").split("\t") == ["chr1", "1", "4", "pELS:E2"] and Path(f"{out}.tbi").exists()


def _bgzip_vcf(path: Path, pos=5):
    path.with_suffix("").write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                                    f"1\t{pos}\t.\tA\tC\t.\t.\t.\n")
    subprocess.run(["micromamba", "run", "-n", "vep116", "bgzip", "-f", str(path.with_suffix(""))], check=True)


def test_prepare_custom_clinvar_indexed_used_in_place(tmp_path, minimal_cfg):
    vcf = tmp_path / "refs/clinvar.vcf.gz"; vcf.parent.mkdir()
    _bgzip_vcf(vcf)
    Path(f"{vcf}.tbi").write_bytes(b"i")
    minimal_cfg["custom"] = {"ClinVar": {"enabled": True, "file": str(vcf), "format": "vcf", "type": "exact"}}
    assert dl.prepare_custom(minimal_cfg)["ClinVar"] == vcf
    custom = Path(minimal_cfg["data"]["custom_data_dir"])
    assert not custom.exists() or not any(custom.iterdir())  # nothing linked or indexed


def test_prepare_custom_clinvar_unindexed_symlinked_and_indexed(tmp_path, minimal_cfg):
    vcf = tmp_path / "refs/clinvar.vcf.gz"; vcf.parent.mkdir()
    _bgzip_vcf(vcf)
    minimal_cfg["custom"] = {"ClinVar": {"enabled": True, "file": str(vcf), "format": "vcf", "type": "exact"}}
    out = dl.prepare_custom(minimal_cfg)["ClinVar"]
    assert out.parent == Path(minimal_cfg["data"]["custom_data_dir"])
    assert out.is_symlink() and out.resolve() == vcf.resolve()
    assert Path(f"{out}.tbi").exists() and not Path(f"{vcf}.tbi").exists()


def test_prepare_custom_clinvar_symlink_follows_the_configured_file(tmp_path, minimal_cfg):
    """S4(c): a link left by an earlier config (same file name, other file) is replaced and re-indexed."""
    old = tmp_path / "old/clinvar.vcf.gz"; old.parent.mkdir()
    _bgzip_vcf(old)
    new = tmp_path / "new/clinvar.vcf.gz"; new.parent.mkdir()
    _bgzip_vcf(new, pos=7)
    minimal_cfg["custom"] = {"ClinVar": {"enabled": True, "file": str(old), "format": "vcf", "type": "exact"}}
    link = dl.prepare_custom(minimal_cfg)["ClinVar"]
    assert link.resolve() == old.resolve()
    minimal_cfg["custom"]["ClinVar"]["file"] = str(new)
    assert dl.prepare_custom(minimal_cfg)["ClinVar"] == link
    assert link.is_symlink() and link.resolve() == new.resolve()
    got = subprocess.run(["tabix", str(link), "1"], capture_output=True, text=True, check=True).stdout
    assert [l.split("\t")[1] for l in got.splitlines()] == ["7"]            # the index is the new file's


def test_prepare_custom_missing_input(minimal_cfg, tmp_path):
    minimal_cfg["custom"] = {"ENCODE_cCRE": {"enabled": True, "source": str(tmp_path / "nope.bed")}}
    with pytest.raises(ConfigError, match="nope.bed"):
        dl.prepare_custom(minimal_cfg)


def test_prepare_custom_skips_disabled(minimal_cfg):
    minimal_cfg["custom"] = {"ENCODE_cCRE": {"enabled": False, "source": "/nonexistent"}}
    assert dl.prepare_custom(minimal_cfg) == {}


def test_run_prep_fills_placeholders_and_fails_loudly(tmp_path):
    final = tmp_path / "d/f.gz"; final.parent.mkdir()
    dl.run_prep(["echo {path} > {dir}/p.txt", "awk 'BEGIN{print 1}' >> {dir}/p.txt"], final, "vep116")
    assert (final.parent / "p.txt").read_text().split() == [str(final), "1"]
    with pytest.raises(RuntimeError, match="prep"):
        dl.run_prep(["false", "echo never > {dir}/n.txt"], final, "vep116")
    assert not (final.parent / "n.txt").exists()


def test_revel_prep_from_plugins_yaml_survives_sigpipe(tmp_path):
    """The real REVEL prep recipe on a zip big enough that `head -n1` closes the pipe while unzip is still writing.

    Regression: under `bash -euo pipefail` unzip dies of SIGPIPE (141) in the header-extraction pipeline.
    """
    spec = load_registry(PLUGINS_YAML)["REVEL"].files[0]
    final = tmp_path / spec.dest
    final.parent.mkdir(parents=True)
    header = "chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid"
    rng = random.Random(7)
    kept, n_dot = [], 0
    lines = [header]
    for i in range(200_000):
        chrom = rng.choice(["1", "2", "10", "X"])
        if i % 50_000 == 7:  # a few unmapped rows: grch38_pos "." must be dropped
            lines.append(f"{chrom},{i + 1},.,A,C,M,V,0.5,ENST{i:011d}")
            n_dot += 1
            continue
        pos = rng.randint(1, 10**8)
        lines.append(f"{chrom},{i + 1},{pos},A,C,M,V,0.123,ENST{i:011d}")
        kept.append((chrom, pos))
    assert n_dot == 4
    csv = ("\n".join(lines) + "\n").encode()
    assert len(csv) > 4_000_000  # >> the 64 KiB pipe buffer
    with zipfile.ZipFile(final.parent / "revel.zip", "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("revel_with_transcript_ids", csv)

    dl.run_prep(list(spec.prep), final, "vep116")

    assert final.exists() and Path(f"{final}.tbi").exists()
    assert not (final.parent / "hdr.tsv").exists() and not (final.parent / "revel.zip").exists()
    with gzip.open(final, "rt") as fh:
        out = fh.read().splitlines()
    assert out[0].startswith("#chr\thg19_pos\tgrch38_pos\t")
    rows = [r.split("\t") for r in out[1:]]
    assert len(rows) == len(kept) and all(r[2] != "." for r in rows)
    keys = [(r[0], int(r[2])) for r in rows]
    assert keys == sorted(keys)  # chr (bytewise, LC_ALL=C) then grch38_pos numeric
    assert sorted(keys) == sorted(kept)


def test_revel_prep_missing_zip_member_still_fails(tmp_path):
    """Disabling pipefail for the header pipeline must not mask a zip without the expected member."""
    spec = load_registry(PLUGINS_YAML)["REVEL"].files[0]
    final = tmp_path / spec.dest
    final.parent.mkdir(parents=True)
    with zipfile.ZipFile(final.parent / "revel.zip", "w") as z:
        z.writestr("something_else", "a,b\n")
    with pytest.raises(RuntimeError, match="prep command failed"):
        dl.run_prep(list(spec.prep), final, "vep116")
    assert not final.exists()


def test_append_manifest(tmp_path):
    f = tmp_path / "f"; f.write_bytes(b"abc")
    dl.append_manifest(tmp_path, "X", "https://u/f", f)
    dl.append_manifest(tmp_path, "Y", "https://u/g", f)
    rows = [l.split("\t") for l in (tmp_path / "manifest.tsv").read_text().splitlines()]
    assert rows[0] == ["name", "url", "bytes", "md5", "date"]
    assert rows[1][:4] == ["X", "https://u/f", "3", "900150983cd24fb0d6963f7d28e17f72"]
    assert len(rows) == 3


# ------------------------------------------------------------------ download_all (file:// only)

def _registry(tmp_path, entry):
    p = tmp_path / "reg.yaml"
    p.write_text(yaml.safe_dump({"Toy": {"module": "Toy", "category": "public", "args": "Toy,file={file}",
                                         "output_fields": ["T"], "evidence_class": "consequence",
                                         "files": {"file": entry}}}))
    return load_registry(p)


def _entry(tmp_path, **kw):
    src = tmp_path / "src.tsv"; src.write_text("a\tb\n" * 20)
    sum_ = dl.bsd_sum(src)
    ck = tmp_path / "CHECKSUMS"; ck.write_text(f"{sum_[0]:05d} {sum_[1]} src.tsv\n")
    entry = {"dest": "T/src.tsv", "downloads": [{"url": src.as_uri(), "as": "T/src.tsv", "size": src.stat().st_size}],
             "index": "tbi", "checksums_url": ck.as_uri()}
    entry.update(kw)
    return entry


def _plugins(cfg, reg):
    cfg["plugins"] = {"Toy": {"enabled": True}}
    return enabled_plugins(cfg, reg)


def test_download_all_fetch_checksum_prep_manifest(tmp_path, minimal_cfg):
    idx = tmp_path / "src.tsv.tbi"; idx.write_bytes(b"IDX")
    reg = _registry(tmp_path, _entry(tmp_path, index_url=idx.as_uri(), prep=["echo prepped > {dir}/prepped.txt"]))
    plugins = _plugins(minimal_cfg, reg)
    done = dl.download_all(minimal_cfg, plugins)
    root = Path(minimal_cfg["data"]["root"]); pdir = Path(minimal_cfg["data"]["plugin_data_dir"])
    assert [i.name for i in done] == ["Toy"]
    assert (pdir / "T/src.tsv").read_text() == "a\tb\n" * 20
    assert (pdir / "T/src.tsv.tbi").read_bytes() == b"IDX"
    assert (pdir / "T/prepped.txt").read_text() == "prepped\n"
    assert not list(pdir.rglob("*.part"))
    rows = [l.split("\t") for l in (root / "manifest.tsv").read_text().splitlines()[1:]]
    assert [r[0] for r in rows] == ["Toy", "Toy"] and rows[0][1].startswith("file://")
    assert dl.plan_downloads(minimal_cfg, plugins) == []
    assert dl.download_all(minimal_cfg, plugins) == []  # nothing left to do


def test_download_all_prep_creates_index(tmp_path, minimal_cfg):
    reg = _registry(tmp_path, _entry(tmp_path, prep=["touch {path}.tbi"]))
    dl.download_all(minimal_cfg, _plugins(minimal_cfg, reg))
    assert Path(minimal_cfg["data"]["plugin_data_dir"], "T/src.tsv.tbi").exists()


def test_download_all_index_never_created_is_an_error(tmp_path, minimal_cfg):
    reg = _registry(tmp_path, _entry(tmp_path))
    with pytest.raises(RuntimeError, match="index"):
        dl.download_all(minimal_cfg, _plugins(minimal_cfg, reg))


def test_download_all_bad_checksum_removes_file(tmp_path, minimal_cfg):
    entry = _entry(tmp_path, prep=["touch {path}.tbi"])
    Path(entry["checksums_url"][len("file://"):]).write_text("00001 1 src.tsv\n")
    reg = _registry(tmp_path, entry)
    plugins = _plugins(minimal_cfg, reg)
    with pytest.raises(RuntimeError, match="checksum"):
        dl.download_all(minimal_cfg, plugins)
    assert len(dl.plan_downloads(minimal_cfg, plugins)) == 1  # not mistaken for complete


def test_download_all_checksum_runs_before_the_final_name_exists(tmp_path, minimal_cfg, monkeypatch):
    """S4(a): verification interrupted after the rename left an unverified final file that later runs skipped."""
    plugins = _plugins(minimal_cfg, _registry(tmp_path, _entry(tmp_path, prep=["touch {path}.tbi"])))
    final = Path(minimal_cfg["data"]["plugin_data_dir"]) / "T/src.tsv"
    seen = []

    def interrupted(path, *a, **kw):
        seen.append((Path(path).name, final.exists()))
        raise KeyboardInterrupt
    monkeypatch.setattr(dl, "verify_checksums", interrupted)
    with pytest.raises(KeyboardInterrupt):
        dl.download_all(minimal_cfg, plugins)
    assert seen == [("src.tsv.part", False)] and not final.exists()
    assert [i.action for i in dl.plan_downloads(minimal_cfg, plugins)] == ["fetch"]
    monkeypatch.undo()
    dl.download_all(minimal_cfg, plugins)                     # the complete .part is verified, then renamed
    assert final.read_text() == "a\tb\n" * 20 and not final.with_name("src.tsv.part").exists()


def test_download_all_failed_validation_removes_file(tmp_path, minimal_cfg):
    reg = _registry(tmp_path, _entry(tmp_path, index=None, checksums_url=None, validate={"min_lines": 1000}))
    plugins = _plugins(minimal_cfg, reg)
    with pytest.raises(RuntimeError, match="min_lines"):
        dl.download_all(minimal_cfg, plugins)
    assert not Path(minimal_cfg["data"]["plugin_data_dir"], "T/src.tsv").exists()
    assert len(dl.plan_downloads(minimal_cfg, plugins)) == 1


def test_dry_run_changes_nothing(tmp_path, minimal_cfg, capsys):
    reg = _registry(tmp_path, _entry(tmp_path, prep=["touch {path}.tbi"]))
    dl.download_all(minimal_cfg, _plugins(minimal_cfg, reg), dry_run=True)
    assert not Path(minimal_cfg["data"]["root"]).exists()
    assert "Toy" in capsys.readouterr().out


def test_user_path_is_never_downloaded_only_verified(tmp_path, minimal_cfg):
    reg = _registry(tmp_path, _entry(tmp_path))
    own = tmp_path / "own/file.tsv.gz"
    minimal_cfg["plugins"] = {"Toy": {"enabled": True, "file": str(own)}}
    plugins = enabled_plugins(minimal_cfg, reg)
    (item,) = dl.plan_downloads(minimal_cfg, plugins)
    assert item.action == "missing" and item.final == own
    with pytest.raises(ConfigError, match="own/file.tsv.gz"):
        dl.download_all(minimal_cfg, plugins)
    assert not own.parent.exists()
    own.parent.mkdir(); own.write_bytes(b"z"); Path(f"{own}.tbi").write_bytes(b"i")
    assert dl.plan_downloads(minimal_cfg, plugins) == []


def test_plan_totals_and_manual_text(tmp_path, minimal_cfg, capsys):
    reg = _registry(tmp_path, _entry(tmp_path, manual="Get it from the portal."))
    own = tmp_path / "nowhere.tsv.gz"
    minimal_cfg["plugins"] = {"Toy": {"enabled": True, "file": str(own)}}
    items = dl.plan_downloads(minimal_cfg, enabled_plugins(minimal_cfg, reg), include_skipped=True)
    dl.print_plan(minimal_cfg, items)
    out = capsys.readouterr().out
    assert "Get it from the portal." in out and str(own) in out and "MISSING" in out
    assert "free" in out.lower() and "total" in out.lower()


# ------------------------------------------------------------------ VEP cache

def _make_cache_tarball(tmp_path, info="var_type\ttabix\n", version="116_GRCh38"):
    work = tmp_path / "tarsrc"; (work / f"homo_sapiens/{version}").mkdir(parents=True)
    (work / f"homo_sapiens/{version}/info.txt").write_text(info)
    tarball = tmp_path / "srv/homo_sapiens_vep_116_GRCh38.tar.gz"; tarball.parent.mkdir()
    with tarfile.open(tarball, "w:gz") as t:
        t.add(work / "homo_sapiens", arcname="homo_sapiens")
    s, k = subprocess.run(["sum", "-r", str(tarball)], capture_output=True, text=True).stdout.split()[:2]
    (tmp_path / "srv/CHECKSUMS").write_text(f"{s} {k} {tarball.name}\n")
    return tarball, tmp_path / "srv/CHECKSUMS"


def test_cache_constants():
    cfg = {"cache_version": 116, "assembly": "GRCh38"}
    assert dl.cache_url(cfg) == ("https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/"
                                 "homo_sapiens_vep_116_GRCh38.tar.gz")
    assert dl.cache_checksums_url(cfg).endswith("release-116/variation/indexed_vep_cache/CHECKSUMS")
    assert dl.cache_size(cfg) == 27644657162


def test_download_cache_full_flow(tmp_path, minimal_cfg):
    tarball, ck = _make_cache_tarball(tmp_path)
    cache = Path(minimal_cfg["data"]["cache_dir"])
    out = dl.download_cache(minimal_cfg, url=tarball.as_uri(), checksums_url=ck.as_uri(), size=tarball.stat().st_size)
    assert out == cache / "homo_sapiens/116_GRCh38"
    assert "var_type\ttabix" in (out / "info.txt").read_text()
    assert not (cache / ".dl" / tarball.name).exists()
    assert "homo_sapiens_vep_116_GRCh38" in (Path(minimal_cfg["data"]["root"]) / "manifest.tsv").read_text()
    # second call is a no-op: the source may vanish
    tarball.unlink()
    assert dl.download_cache(minimal_cfg, url=tarball.as_uri(), checksums_url=ck.as_uri(), size=1) == out


def test_download_cache_interrupted_extraction_is_not_ready_and_resumes(tmp_path, minimal_cfg, monkeypatch):
    tarball, ck = _make_cache_tarball(tmp_path)
    cache = Path(minimal_cfg["data"]["cache_dir"])
    real_run = subprocess.run
    def killed_tar(cmd, *a, **kw):
        if cmd[0] != "tar":
            return real_run(cmd, *a, **kw)
        real_run(cmd, *a, **kw)                         # tar got as far as writing info.txt ...
        raise subprocess.CalledProcessError(1, cmd)     # ... and then died
    monkeypatch.setattr(dl.subprocess, "run", killed_tar)
    with pytest.raises(subprocess.CalledProcessError):
        dl.download_cache(minimal_cfg, url=tarball.as_uri(), checksums_url=ck.as_uri(), size=tarball.stat().st_size)
    assert not dl.cache_ready(minimal_cfg)
    assert not (cache / "homo_sapiens/116_GRCh38").exists()
    assert dl.cache_item(minimal_cfg).action == "fetch"
    assert (cache / ".dl" / tarball.name).is_file()     # kept for the retry
    monkeypatch.undo()
    tarball.unlink()                                    # retry must not need the source again
    out = dl.download_cache(minimal_cfg, url="file:///nonexistent", checksums_url=ck.as_uri(),
                            size=(cache / ".dl" / tarball.name).stat().st_size)
    assert dl.cache_ready(minimal_cfg) and out.is_dir()
    assert not (cache / ".extract").exists() and not (cache / ".dl" / tarball.name).exists()


def test_download_cache_replaces_partial_leftover_without_staging(tmp_path, minimal_cfg):
    """A truncated tree left by an older, non-staged extraction is replaced, not trusted."""
    tarball, ck = _make_cache_tarball(tmp_path)
    leftover = Path(minimal_cfg["data"]["cache_dir"]) / "homo_sapiens/116_GRCh38"
    leftover.mkdir(parents=True); (leftover / "stale.txt").write_text("x")
    out = dl.download_cache(minimal_cfg, url=tarball.as_uri(), checksums_url=ck.as_uri(), size=tarball.stat().st_size)
    assert (out / "info.txt").exists() and not (out / "stale.txt").exists()


def test_download_cache_reuses_complete_tarball_without_fetching(tmp_path, minimal_cfg, monkeypatch):
    tarball, ck = _make_cache_tarball(tmp_path)
    kept = Path(minimal_cfg["data"]["cache_dir"]) / ".dl" / tarball.name
    kept.parent.mkdir(parents=True); kept.write_bytes(tarball.read_bytes())
    def no_fetch(*a, **kw):
        raise AssertionError("fetch must not be called when the tarball is already complete")
    monkeypatch.setattr(dl, "fetch", no_fetch)
    out = dl.download_cache(minimal_cfg, url="file:///nonexistent", checksums_url=ck.as_uri(),
                            size=tarball.stat().st_size)
    assert (out / "info.txt").exists() and not kept.exists()


def test_download_cache_wrong_size_tarball_is_refetched(tmp_path, minimal_cfg):
    tarball, ck = _make_cache_tarball(tmp_path)
    kept = Path(minimal_cfg["data"]["cache_dir"]) / ".dl" / tarball.name
    kept.parent.mkdir(parents=True); kept.write_bytes(b"truncated")
    out = dl.download_cache(minimal_cfg, url=tarball.as_uri(), checksums_url=ck.as_uri(), size=tarball.stat().st_size)
    assert (out / "info.txt").exists()


def test_download_cache_checksum_mismatch_deletes_tarball(tmp_path, minimal_cfg):
    tarball, ck = _make_cache_tarball(tmp_path)
    ck.write_text(f"00001 1 {tarball.name}\n")
    kept = Path(minimal_cfg["data"]["cache_dir"]) / ".dl" / tarball.name
    kept.parent.mkdir(parents=True); kept.write_bytes(tarball.read_bytes())     # reuse path
    with pytest.raises(RuntimeError, match="checksum"):
        dl.download_cache(minimal_cfg, url="file:///nonexistent", checksums_url=ck.as_uri(), size=tarball.stat().st_size)
    assert not kept.exists()
    with pytest.raises(RuntimeError, match="checksum"):                          # fresh-fetch path
        dl.download_cache(minimal_cfg, url=tarball.as_uri(), checksums_url=ck.as_uri(), size=tarball.stat().st_size)
    assert not kept.exists() and not kept.with_name(kept.name + ".part").exists()


def test_download_cache_checksum_runs_before_the_tarball_is_renamed(tmp_path, minimal_cfg, monkeypatch):
    """S4(a): an unverified complete tarball under the final name would be reused as complete."""
    tarball, ck = _make_cache_tarball(tmp_path)
    kept = Path(minimal_cfg["data"]["cache_dir"]) / ".dl" / tarball.name
    seen = []

    def interrupted(path, *a, **kw):
        seen.append((Path(path).name, kept.exists()))
        raise KeyboardInterrupt
    monkeypatch.setattr(dl, "verify_checksums", interrupted)
    with pytest.raises(KeyboardInterrupt):
        dl.download_cache(minimal_cfg, url=tarball.as_uri(), checksums_url=ck.as_uri(), size=tarball.stat().st_size)
    assert seen == [(tarball.name + ".part", False)] and not kept.exists()


def test_download_cache_checksum_mismatch_stops_before_extract(tmp_path, minimal_cfg):
    tarball, ck = _make_cache_tarball(tmp_path)
    ck.write_text(f"00001 1 {tarball.name}\n")
    with pytest.raises(RuntimeError, match="checksum"):
        dl.download_cache(minimal_cfg, url=tarball.as_uri(), checksums_url=ck.as_uri(), size=tarball.stat().st_size)
    assert not (Path(minimal_cfg["data"]["cache_dir"]) / "homo_sapiens").exists()


def test_download_cache_requires_tabix_indexed_cache(tmp_path, minimal_cfg):
    tarball, ck = _make_cache_tarball(tmp_path, info="var_type\tsqlite\n")
    with pytest.raises(RuntimeError, match="tabix"):
        dl.download_cache(minimal_cfg, url=tarball.as_uri(), checksums_url=ck.as_uri(), size=tarball.stat().st_size)


def test_download_cache_disk_space_guard(tmp_path, minimal_cfg):
    tarball, ck = _make_cache_tarball(tmp_path)
    with pytest.raises(RuntimeError, match="free"):
        dl.download_cache(minimal_cfg, url=tarball.as_uri(), checksums_url=ck.as_uri(), size=10**18)
    assert not (Path(minimal_cfg["data"]["cache_dir"]) / ".dl").exists() or \
        not any((Path(minimal_cfg["data"]["cache_dir"]) / ".dl").iterdir())


def test_fetch_uses_pget_only_for_big_https(tmp_path, monkeypatch):
    calls = []
    def fake_run(cmd, check, env):
        calls.append((cmd, env["SSL_CERT_FILE"]))
        Path(cmd[cmd.index("-o") + 1] if "-o" in cmd else cmd[3]).write_bytes(b"x" * 5)
    monkeypatch.setattr(dl.subprocess, "run", fake_run)
    pget = tmp_path / "pget.sh"; pget.write_text("")
    dl.fetch("https://h/big", tmp_path / "big", size=5, pget=pget, big=4)
    dl.fetch("https://h/small", tmp_path / "small", size=5, pget=pget, big=10)
    dl.fetch("file:///f", tmp_path / "local", size=5, pget=pget, big=4)
    assert calls[0][0] == ["bash", str(pget), "https://h/big", str(tmp_path / "big.part"), "8", "512"]
    assert calls[1][0][0] == "curl" and calls[2][0][0] == "curl"
    assert all(Path(c[1]).is_file() for c in calls)  # https_env gave every call a real CA bundle


def test_content_length_head(tmp_path):
    f = tmp_path / "f"; f.write_bytes(b"x" * 123)
    assert dl.content_length(f.as_uri()) == 123


# ------------------------------------------------------------------ CLI

def test_cli_download_dry_run(minimal_cfg, capsys):
    from vep_pipeline import cli
    minimal_cfg["plugins"] = {"LOEUF": {"enabled": True}}
    rc = cli.HANDLERS["download"](minimal_cfg, Namespace(dry_run=True, only=None, skip_cache=False))
    out = capsys.readouterr().out
    assert rc == 0
    assert "LOEUF" in out and "354938316" in out and "loeuf_dataset_grch38.tsv.gz" in out
    assert "fetch" in out and "27644657162" in out and "total" in out.lower() and "free" in out.lower()
    assert not Path(minimal_cfg["data"]["root"]).exists()


def test_cli_download_only_and_skip_cache(minimal_cfg, capsys):
    from vep_pipeline import cli
    minimal_cfg["plugins"] = {"LOEUF": True, "MaveDB": True}
    cli.HANDLERS["download"](minimal_cfg, Namespace(dry_run=True, only=["LOEUF"], skip_cache=True))
    out = capsys.readouterr().out
    assert "LOEUF" in out and "MaveDB" not in out and "27644657162" not in out


def test_cli_download_unknown_only_name(minimal_cfg):
    from vep_pipeline import cli
    with pytest.raises(ConfigError, match="Nope"):
        cli.HANDLERS["download"](minimal_cfg, Namespace(dry_run=True, only=["Nope"], skip_cache=True))


def test_cli_parser_download_flags():
    from vep_pipeline.cli import build_parser
    a = build_parser().parse_args(["download", "--dry-run", "--only", "A", "B", "--skip-cache"])
    assert a.dry_run and a.only == ["A", "B"] and a.skip_cache


def _glob_registry(tmp_path):
    p = tmp_path / "reg_glob.yaml"
    p.write_text(yaml.safe_dump({"Dir": {"module": "Dir", "category": "manual", "args": "Dir,snv={snv}",
                                         "output_fields": ["D"], "evidence_class": "consequence",
                                         "files": {"snv": {"dest": "Dir/snv", "glob": "chr*.tsv.gz",
                                                           "index": "tbi", "manual": "Get it."}}}}))
    return load_registry(p)


def test_plan_glob_directory_needs_every_file_indexed(tmp_path, minimal_cfg):
    d = tmp_path / "snv"
    minimal_cfg["plugins"] = {"Dir": {"enabled": True, "snv": str(d)}}
    plugins = enabled_plugins(minimal_cfg, _glob_registry(tmp_path))
    (item,) = dl.plan_downloads(minimal_cfg, plugins)
    assert item.action == "missing" and "no files matching chr*.tsv.gz" in item.reason
    d.mkdir()
    for c in ("chr1", "chr2"):
        (d / f"{c}.tsv.gz").write_bytes(b"z")
    (d / "chr1.tsv.gz.tbi").write_bytes(b"i")
    (item,) = dl.plan_downloads(minimal_cfg, plugins)
    assert item.action == "missing" and "1 of 2" in item.reason and "chr2.tsv.gz" in item.reason
    (d / "chr2.tsv.gz.tbi").write_bytes(b"i")
    assert dl.plan_downloads(minimal_cfg, plugins) == []
    assert not (Path(minimal_cfg["data"]["plugin_data_dir"]) / "Dir").exists()
