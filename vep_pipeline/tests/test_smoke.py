import argparse
import subprocess
from pathlib import Path

import pytest

from vep_pipeline import cli, smoke
from vep_pipeline.registry import PluginSpec, ResolvedPlugin

SEQ = "ACGT" * 250                       # contig 17, 1000 bp; base at 1-based p is SEQ[p-1]
REGION = "17:100-900"


def _sh(*cmd):
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def _rp(name, fields, control, eligible=(), paths=None):
    spec = PluginSpec(name=name, module=name, category="public", output_fields=tuple(fields),
                      control=control, eligible_consequences=tuple(eligible))
    return ResolvedPlugin(spec, paths or {}, {}, name)


def _alt(p):
    return "A" if SEQ[p - 1] != "A" else "C"


def _bgzip_tabix(path, preset_args):
    _sh("bgzip", "-f", str(path))
    _sh("tabix", "-f", *preset_args, f"{path}.gz")
    return Path(f"{path}.gz")


# ----------------------------------------------------------------------------- smoke_report

CSQ_HEADER = ('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
              'Format: Allele|Consequence|A_val|B_val">')


def _vep_vcf(tmp_path, records, header=CSQ_HEADER):
    lines = ["##fileformat=VCFv4.2", "##contig=<ID=17,length=1000>", header,
             "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])]
    for pos, ref, alt, csq in records:
        lines.append("\t".join(["17", str(pos), ".", ref, alt, ".", ".", f"CSQ={csq}", "GT", "0/1"]))
    p = tmp_path / "vep.vcf"
    p.write_text("\n".join(lines) + "\n")
    return p


@pytest.fixture
def vep_vcf(tmp_path):
    # r1 missense with A_val; r2 missense, empty A_val and `-` B_val; r3 synonymous with A_val (not eligible)
    return _vep_vcf(tmp_path, [
        (100, "A", "C", "C|missense_variant|0.9|,C|synonymous_variant||"),
        (200, "A", "C", "C|missense_variant|-|-"),
        (300, "A", "C", "C|synonymous_variant|0.1|"),
    ])


K = {p: f"17:{p}:A>C" for p in (100, 200, 300)}


def test_smoke_report_statuses_and_coverage(vep_vcf):
    plugins = [
        _rp("A", ["A_val"], {"source": "data", "cols": [1, 2, 3, 4]}, eligible=["missense_variant"]),
        _rp("B", ["B_val"], {"source": "data", "cols": [1, 2, 3, 4]}, eligible=["missense_variant"]),
        _rp("C", ["C_val"], {"source": "data", "cols": [1, 2, 3, 4]}),
    ]
    keys = {"A": {K[100]}, "B": {K[200]}, "C": {K[100]}, "ClinVar": set(K.values())}
    rows = {r["plugin"]: r for r in smoke.smoke_report(vep_vcf, plugins, keys)}
    assert list(rows["A"]) == list(smoke.REPORT_COLUMNS)
    a, b, c = rows["A"], rows["B"], rows["C"]
    assert (a["status"], a["header_ok"], a["controls_n"], a["controls_annotated"]) == ("OK", True, 1, 1)
    # eligible = the two missense records; only r1 has a value on a missense entry (r3 is not eligible)
    assert (a["eligible_n"], a["annotated_n"], a["coverage"]) == (2, 1, 0.5)
    assert (b["status"], b["header_ok"], b["controls_annotated"]) == ("FAIL", True, 0)
    assert "none of the 1" in b["message"]
    assert (c["status"], c["header_ok"]) == ("FAIL", False)
    assert "C_val" in c["message"]


def test_smoke_report_eligible_empty_means_all_records(vep_vcf):
    rows = smoke.smoke_report(vep_vcf, [_rp("A", ["A_val"], {"source": "data", "cols": [1, 2, 3, 4]})],
                              {"A": {K[100]}, "ClinVar": set(K.values())})
    assert (rows[0]["eligible_n"], rows[0]["annotated_n"], rows[0]["coverage"]) == (3, 2, 0.667)


def test_smoke_report_universe_is_the_clinvar_controls(vep_vcf):
    rows = smoke.smoke_report(vep_vcf, [_rp("A", ["A_val"], {"source": "header"})], {"ClinVar": {K[300]}})
    assert (rows[0]["eligible_n"], rows[0]["annotated_n"]) == (1, 1)


def test_smoke_report_header_only_clinvar_unchecked_and_no_data(vep_vcf):
    plugins = [
        _rp("A", ["A_val"], {"source": "header"}),
        _rp("B", ["B_val"], {"source": "clinvar", "mc": "missense_variant"}),
        _rp("Nc", ["A_val"], {}),
        _rp("D", ["A_val"], {"source": "vcf"}),
        _rp("E", ["A_val"], {"source": "clinvar", "mc": "missense_variant"}),
    ]
    rows = {r["plugin"]: r for r in smoke.smoke_report(vep_vcf, plugins, {"B": set(), "E": {K[100]}})}
    assert rows["A"]["status"] == "OK" and rows["A"]["controls_n"] == 0
    assert rows["B"]["status"] == "unchecked" and "no ClinVar control" in rows["B"]["message"]
    assert rows["Nc"]["status"] == "unchecked"
    assert rows["D"]["status"] == "FAIL" and "no data for this plugin" in rows["D"]["message"]
    assert "smoke.region" in rows["D"]["message"]
    assert rows["E"]["status"] == "OK"


CV_ON = {"ClinVar": {"enabled": True, "file": "/cv.vcf.gz"}}


def test_smoke_report_clinvar_source_plugin_fails_without_clinvar_track(vep_vcf):
    """S1(b): with the ClinVar track off a clinvar-source plugin has no controls; that is a FAIL, not `unchecked`."""
    plugins = [_rp("B", ["B_val"], {"source": "clinvar", "mc": "missense_variant"})]
    for custom in ({}, {"ClinVar": {"enabled": False, "file": "/cv.vcf.gz"}}, {"ClinVar": {"enabled": True}}):
        (row,) = smoke.smoke_report(vep_vcf, plugins, {"B": set()}, custom=custom)
        assert row["status"] == "FAIL" and "custom.ClinVar" in row["message"], custom
    row, track = smoke.smoke_report(vep_vcf, plugins, {"B": set(), "ClinVar": set()}, custom=CV_ON)
    assert row["status"] == "unchecked" and track["plugin"] == "custom:ClinVar"


def test_smoke_report_fails_plugin_whose_controls_disagree_with_fasta(vep_vcf):
    """S2: more than 2% of a data plugin's controls with a REF the FASTA does not have -> wrong-assembly FAIL."""
    plugin = _rp("A", ["A_val"], {"source": "data", "cols": [1, 2, 3, 4]})
    ctl = {K[100]} | {f"17:{p}:A>C" for p in range(1000, 1059)}                  # 60 controls, K[100] annotated
    others = sorted(ctl - {K[100]})
    (ok,) = smoke.smoke_report(vep_vcf, [plugin], {"A": ctl}, ref_mismatch={others[0]})          # 1/60 = 1.7%
    assert ok["status"] == "OK"
    (bad,) = smoke.smoke_report(vep_vcf, [plugin], {"A": ctl}, ref_mismatch=set(others[:2]))     # 2/60 = 3.3%
    assert bad["status"] == "FAIL" and "possible wrong assembly in plugin data" in bad["message"]
    assert "2/60" in bad["message"]


TRACK_HEADER = ('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
                'Format: Allele|Consequence|ClinVar_CLNSIG|ENCODE_cCRE">')


def test_smoke_report_custom_track_rows(tmp_path):
    """S3: ClinVar controls must get a ClinVar_CLNSIG value (FAIL if none does); the cCRE overlap count is info."""
    vcf = _vep_vcf(tmp_path, [(100, "A", "C", "C|missense_variant|Pathogenic|"),
                              (200, "A", "C", "C|intron_variant||pELS:EH38E1"),
                              (300, "A", "C", "C|intron_variant||")], header=TRACK_HEADER)
    custom = {**CV_ON, "ENCODE_cCRE": {"enabled": True, "source": "/c.bed"}}
    cv, cc = smoke.smoke_report(vcf, [], {"ClinVar": {K[100], K[200]}}, custom=custom)
    assert list(cv) == list(smoke.REPORT_COLUMNS)
    assert (cv["plugin"], cv["status"], cv["header_ok"], cv["controls_n"], cv["controls_annotated"]) == (
        "custom:ClinVar", "OK", True, 2, 1)
    assert (cc["plugin"], cc["status"], cc["controls_n"], cc["controls_annotated"]) == ("custom:ENCODE_cCRE", "info", 3, 1)
    assert "1/3" in cc["message"]
    (cv,) = smoke.smoke_report(vcf, [], {"ClinVar": {K[200], K[300]}}, custom=CV_ON)
    assert cv["status"] == "FAIL" and "ClinVar_CLNSIG" in cv["message"]
    (cv,) = smoke.smoke_report(vcf, [], {"ClinVar": set()}, custom=CV_ON)
    assert cv["status"] == "FAIL" and "smoke.region" in cv["message"]
    (cv,) = smoke.smoke_report(_vep_vcf(tmp_path, [(100, "A", "C", "C|missense_variant|0.9|")]), [],
                               {"ClinVar": {K[100]}}, custom=CV_ON)
    assert cv["status"] == "FAIL" and cv["header_ok"] is False
    assert smoke.smoke_report(vcf, [], {"ClinVar": {K[100]}}, custom={}) == []


def test_smoke_report_missing_header_fails_even_without_control(vep_vcf):
    rows = smoke.smoke_report(vep_vcf, [_rp("N", ["NOPE"], {"source": "header"}), _rp("U", ["NOPE"], {})], {})
    assert [r["status"] for r in rows] == ["FAIL", "FAIL"]


def test_smoke_report_reads_bgzipped_vep_vcf(vep_vcf):
    gz = _bgzip_tabix(vep_vcf, ["-p", "vcf"])
    rows = smoke.smoke_report(gz, [_rp("A", ["A_val"], {"source": "header"})], {})
    assert rows[0]["status"] == "OK" and rows[0]["header_ok"] is True


def test_write_report_tsv(tmp_path):
    rows = [{"plugin": "A", "status": "OK", "header_ok": True, "controls_n": 1, "controls_annotated": 1,
             "eligible_n": 0, "annotated_n": 0, "coverage": None, "message": "m"}]
    path = smoke.write_report(rows, tmp_path / "sub" / "smoke.tsv")
    head, row = path.read_text().splitlines()
    assert head.split("\t") == ["plugin", "status", "header_ok", "controls_n", "controls_annotated",
                                "eligible_n", "annotated_n", "coverage", "message"]
    assert row.split("\t") == ["A", "OK", "yes", "1", "1", "0", "0", "", "m"]


# ----------------------------------------------------------------------------- build_smoke_vcf

CLNSIG_OK = "Pathogenic"
TWO = "criteria_provided,_multiple_submitters,_no_conflicts"


@pytest.fixture
def smoke_setup(minimal_cfg, fasta_factory, tmp_path):
    cfg = minimal_cfg
    cfg["input"]["reference_fasta"] = str(fasta_factory({"17": 1000, "MT": 100}, seq={"17": SEQ, "MT": "A" * 100}))
    # plugin data: chr-prefixed contig (the registry's data files differ from the FASTA naming), 60 SNVs in the
    # region + decoys: outside the region, a multi-allelic row, a non-ACGT row
    rows = [("chr17", p, SEQ[p - 1], _alt(p), "0.5") for p in range(110, 110 + 60 * 5, 5)]
    rows += [("chr17", 50, "A", "C", "x"), ("chr17", 950, "A", "C", "x"),
             ("chr17", 301, SEQ[300], "T,G", "x"), ("chr17", 303, "A", "-", "x")]
    plugin = _data_plugin(tmp_path, "Toy", rows)

    def cv(pos, mc, sig=CLNSIG_OK, rev=TWO, alt=None):
        alt = alt or _alt(pos)
        return ("17", pos, SEQ[pos - 1], alt, f"CLNSIG={sig};CLNREVSTAT={rev};MC={mc}")
    mis = [cv(p, "SO:0001583|missense_variant") for p in range(120, 120 + 14 * 7, 7)]       # 14 -> pick 10
    non = [cv(p, "SO:0001587|nonsense,SO:0001583|missense_variant") for p in (400, 410)]   # also missense
    utr = [cv(500, "SO:0001623|5_prime_UTR_variant")]
    # any significance / review status is a control: only the consequence (MC) matters
    anysig = [cv(600, "SO:0001583|missense_variant", sig="Benign"),
              cv(610, "SO:0001583|missense_variant", rev="criteria_provided,_single_submitter"),
              cv(620, "SO:0001583|missense_variant", sig="Likely_pathogenic",
                 rev="no_assertion_criteria_provided"),
              cv(640, "SO:0001630|splice_region_variant", sig="Uncertain_significance",
                 rev="criteria_provided,_single_submitter")]
    decoys = [cv(1, "SO:0001583|missense_variant"),                                      # outside the region
              ("17", 660, SEQ[659], _alt(660), f"CLNSIG={CLNSIG_OK};CLNREVSTAT={TWO}"),  # no MC: no stratum
              ("17", 700, SEQ[699:760], SEQ[699], f"CLNSIG={CLNSIG_OK};CLNREVSTAT={TWO};MC=SO:0001589|frameshift_variant")]
    cvrows = sorted(mis + non + utr + anysig + decoys, key=lambda r: r[1])
    clin = tmp_path / "clinvar.vcf"
    clin.write_text("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    + "".join(f"{c}\t{p}\t.\t{r}\t{a}\t.\t.\t{i}\n" for c, p, r, a, i in cvrows))
    clin_gz = _bgzip_tabix(clin, ["-p", "vcf"])
    cfg["custom"] = {"ClinVar": {"enabled": True, "file": str(clin_gz), "format": "vcf", "type": "exact"}}
    return cfg, [plugin]


def _data_plugin(tmp_path, name, rows):
    """A data-source plugin (output field NAME) backed by a tabix-indexed CHROM POS REF ALT SCORE file."""
    body = sorted(((*r, "0.5")[:5] for r in rows), key=lambda r: r[1])
    data = tmp_path / f"{name.lower()}.tsv"
    data.write_text("#CHROM\tPOS\tREF\tALT\tSCORE\n" + "".join("\t".join(map(str, r)) + "\n" for r in body))
    data_gz = _bgzip_tabix(data, ["-s", "1", "-b", "2", "-e", "2", "-c", "#"])
    return _rp(name, [name.upper()], {"source": "data", "cols": [1, 2, 3, 4]}, paths={"file": data_gz})


def _body(vcf):
    out = subprocess.run(["bcftools", "view", "-H", str(vcf)], capture_output=True, text=True, check=True)
    return [l.split("\t") for l in out.stdout.splitlines()]


def test_build_smoke_vcf_is_byte_identical(smoke_setup):
    cfg, plugins = smoke_setup
    first = smoke.build_smoke_vcf(cfg, plugins, region=REGION, seed=7)
    data = first.read_bytes()
    second = smoke.build_smoke_vcf(cfg, plugins, region=REGION, seed=7)
    assert second == first and second.read_bytes() == data
    assert first == Path(cfg["output"]["dir"]) / "smoke" / "smoke.controls.vcf.gz"
    assert Path(f"{first}.tbi").exists()
    assert smoke.build_smoke_vcf(cfg, plugins, region=REGION, seed=8).read_bytes() != data
    assert not list((first.parent).glob("*.raw.vcf")) and not list(first.parent.glob("*.norm.vcf"))


def test_build_smoke_vcf_content(smoke_setup):
    cfg, plugins = smoke_setup
    vcf = smoke.build_smoke_vcf(cfg, plugins, region=REGION, seed=7)
    header = "\n".join(l for l in smoke._read_lines(vcf) if l.startswith("#"))       # the file's own header
    assert "##bcftools" not in header and "##contig=<ID=17,length=1000>" in header
    assert "##contig=<ID=MT" not in header
    assert header.splitlines()[-1].split("\t")[-1] == "S1"
    rows = _body(vcf)
    assert all(r[0] == "17" and r[9] == "0/1" for r in rows)            # FASTA naming, GT 0/1
    assert [int(r[1]) for r in rows] == sorted(int(r[1]) for r in rows)
    pos = {int(r[1]) for r in rows}
    assert not pos & {50, 950, 1, 303, 660, 700}                 # outside region / non-ACGT / no MC / 60-bp REF
    keys = smoke.control_keys_from_vcf(vcf, plugins)
    assert len(keys["Toy"]) == smoke.CONTROLS_PER_FILE
    assert all(int(k.split(":")[1]) in range(110, 410) for k in keys["Toy"])
    # strata: 19 missense records -> 10; the two nonsense records; the 5'UTR record; the VUS splice record
    cv = {r[2] for r in rows if "cv_" in r[2]}
    tags = [t for r in rows for t in r[2].split(";") if t.startswith("cv_")]
    assert tags.count("cv_missense") == 10 and tags.count("cv_nonsense_frameshift") == 2
    assert tags.count("cv_5_prime_UTR") == 1 and tags.count("cv_splice_region") == 1 and "cv_intron" not in tags
    assert 640 in pos
    assert len(keys["ClinVar"]) == len([r for r in rows if "cv_" in r[2]]) and cv
    plain = [_rp("M", ["X"], {"source": "clinvar", "mc": "missense_variant"}),
             _rp("U", ["X"], {"source": "clinvar", "mc": "5_prime_UTR"}),
             _rp("S", ["X"], {"source": "clinvar", "mc": "splice_region"}),
             _rp("H", ["X"], {"source": "header"})]
    k2 = smoke.control_keys_from_vcf(vcf, plain)
    assert len(k2["M"]) >= 10 and len(k2["U"]) == 1 and len(k2["S"]) == 1 and k2["H"] == set()


def test_clinvar_candidates_ignore_significance_and_review_status(smoke_setup):
    """S1(c): ClinVar smoke controls need the consequence class (MC), not pathogenicity or review stars."""
    cfg, _plugins = smoke_setup
    cands = smoke._clinvar_candidates(cfg["custom"]["ClinVar"]["file"], "17", 100, 900, "vep116")
    assert {600, 610, 620, 640, 660} <= {c[1] for c in cands}
    assert not hasattr(smoke, "PLP") and not hasattr(smoke, "TWO_STARS")


def test_ref_mismatches_flags_controls_off_the_fasta(smoke_setup, tmp_path):
    """S2: `bcftools norm -c w` keeps records whose REF is not the FASTA base; ref_mismatches must find them."""
    cfg, plugins = smoke_setup
    bad = _data_plugin(tmp_path, "Bad", [("chr17", p, _alt(p), SEQ[p - 1]) for p in range(112, 400, 5)])
    vcf = smoke.build_smoke_vcf(cfg, plugins + [bad], region=REGION, seed=7)
    keys = smoke.control_keys_from_vcf(vcf, plugins + [bad])
    mism = smoke.ref_mismatches(cfg, vcf)
    assert keys["Bad"] and keys["Bad"] == mism
    assert not keys["Toy"] & mism and not keys["ClinVar"] & mism


def test_build_smoke_vcf_plugin_data_with_no_record_in_region(smoke_setup):
    cfg, plugins = smoke_setup
    vcf = smoke.build_smoke_vcf(cfg, plugins, region="17:960-990", seed=7)
    assert smoke.control_keys_from_vcf(vcf, plugins)["Toy"] == set()


def test_build_smoke_vcf_region_default_from_cfg_and_validation(smoke_setup):
    cfg, plugins = smoke_setup
    cfg["smoke"]["region"] = REGION
    assert smoke.build_smoke_vcf(cfg, plugins).read_bytes() == smoke.build_smoke_vcf(cfg, plugins, REGION).read_bytes()
    with pytest.raises(smoke.ConfigError):
        smoke.build_smoke_vcf(cfg, plugins, region="17-1-2")
    with pytest.raises(smoke.ConfigError):
        smoke.build_smoke_vcf(cfg, plugins, region="18:1-100")


def test_stride_pick_is_spread_deterministic_and_bounded():
    rows = list(range(100))
    assert smoke._stride_pick(rows, 25, 7) == smoke._stride_pick(rows, 25, 7)
    assert len(smoke._stride_pick(rows, 25, 7)) == 25 and smoke._stride_pick(rows, 25, 7)[0] == 7 % 4
    assert smoke._stride_pick(rows[:5], 25, 7) == rows[:5]


def test_contig_naming_candidates():
    assert smoke._contig_names("17") == ["17", "chr17"]
    assert smoke._contig_names("chr17") == ["chr17", "17"]
    assert smoke._contig_names("MT") == ["MT", "chrM"] and smoke._contig_names("chrM") == ["chrM", "MT"]
    assert smoke._ensembl("chrM") == "MT" and smoke._ensembl("chr17") == "17" and smoke._ensembl("MT") == "MT"


# ----------------------------------------------------------------------------- run_smoke and the CLI

def test_run_smoke_uses_normal_pipeline_on_derived_cfg(smoke_setup, tmp_path, monkeypatch):
    cfg, plugins = smoke_setup
    seen = {}
    final = _bgzip_tabix(_vep_vcf(tmp_path, [(100, "A", "C", "C|missense_variant|0.9|")]), ["-p", "vcf"])

    def fake_render(scfg, pl):
        seen["cfg"], seen["plugins"] = scfg, pl
        return "PLAN"

    def fake_run_local(plan, parallel):
        seen["plan"], seen["parallel"] = plan, parallel
        return final

    monkeypatch.setattr(smoke, "render_run", fake_render)
    monkeypatch.setattr(smoke, "run_local", fake_run_local)
    rows = smoke.run_smoke(cfg, plugins, region=REGION, parallel=2)
    scfg = seen["cfg"]
    smoke_dir = Path(cfg["output"]["dir"]) / "smoke"
    assert scfg["input"]["vcf"] == str(smoke_dir / "smoke.controls.vcf.gz") and scfg["input"]["normalize"] is False
    assert scfg["output"]["dir"] == str(smoke_dir) and scfg["sv"] == {"enabled": False}
    assert scfg["sample"] == cfg["sample"] and cfg["input"]["normalize"] is False      # original untouched
    assert cfg["output"]["dir"] != scfg["output"]["dir"]
    assert (seen["plan"], seen["parallel"], seen["plugins"]) == ("PLAN", 2, plugins)
    assert [r["plugin"] for r in rows] == ["Toy", "custom:ClinVar"]
    assert rows[0]["status"] == "FAIL"                                                  # TOY not in the CSQ header
    assert (smoke_dir / "smoke.tsv").read_text().splitlines()[0].startswith("plugin\tstatus\theader_ok")


def test_smoke_render_keeps_g2p_reports_under_smoke_dir(minimal_cfg, vcf_factory):
    """C1: plugins resolved against the production cfg must not write G2P reports into the production reports dir;
    the ini's report paths and the script's `mkdir -p` must name the same smoke directory."""
    from vep_pipeline import jobs
    from vep_pipeline.registry import enabled_plugins, load_registry
    cfg = minimal_cfg
    cfg["plugins"] = {"G2P": {"enabled": True}}
    plugins = enabled_plugins(cfg, load_registry(Path(__file__).resolve().parents[1] / "plugins.yaml"))
    plan = jobs.render_run(smoke.smoke_cfg(cfg, vcf_factory({"17": 1000}, [("17", 100, "A", "C")])), plugins)
    smoke_out = Path(cfg["output"]["dir"]) / "smoke"
    reports = smoke_out / "reports"
    (job,) = plan.chunks
    ini = job.ini.read_text()
    g2p = next(l for l in ini.splitlines() if l.startswith("plugin G2P,"))
    opts = dict(kv.split("=", 1) for kv in g2p.split(",")[1:])
    assert opts["txt_report"] == f"{reports}/g2p.17.txt"
    assert opts["html_report"] == f"{reports}/g2p.17.html"
    assert opts["log_dir"] == f"{smoke_out}/work/17/g2p_log"
    assert f"mkdir -p {reports}\n" in job.script.read_text()
    assert f"{cfg['output']['dir']}/reports" not in ini + job.script.read_text()


def _fake_vep(smoke_vcf, out, fields, values):
    """A VEP-like VCF of every smoke record with one CSQ entry `values(record fields)` in the order of `fields`."""
    lines = []
    for line in smoke._read_lines(smoke_vcf):
        if line.startswith("#CHROM"):
            lines.append('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl '
                         f'VEP. Format: {"|".join(fields)}">')
        if line.startswith("#"):
            lines.append(line)
            continue
        f = line.split("\t")
        f[7] = "CSQ=" + "|".join(values(f))
        lines.append("\t".join(f))
    out.write_text("\n".join(lines) + "\n")
    return out


def test_run_smoke_wires_ref_mismatch_clinvar_track_and_track_rows(smoke_setup, tmp_path, monkeypatch):
    """run_smoke hands the report the REF mismatches (S2), the custom tracks (S1b) and gets the track rows (S3)."""
    cfg, plugins = smoke_setup
    bad = _data_plugin(tmp_path, "Bad", [("chr17", p, _alt(p), SEQ[p - 1]) for p in range(112, 400, 5)])
    cvp = _rp("Cvp", ["TOY"], {"source": "clinvar", "mc": "missense_variant"})
    fields = ["Allele", "Consequence", "TOY", "BAD", "ClinVar_CLNSIG"]
    monkeypatch.setattr(smoke, "render_run", lambda scfg, pl: scfg)
    monkeypatch.setattr(smoke, "run_local", lambda scfg, parallel: _fake_vep(
        Path(scfg["input"]["vcf"]), tmp_path / "vep.vcf", fields,
        lambda f: [f[4], "missense_variant", "0.5", "0.5", "Pathogenic"]))
    rows = {r["plugin"]: r for r in smoke.run_smoke(cfg, plugins + [bad, cvp], region=REGION)}
    assert rows["Toy"]["status"] == "OK" and rows["Cvp"]["status"] == "OK"
    assert rows["Bad"]["status"] == "FAIL" and "possible wrong assembly in plugin data" in rows["Bad"]["message"]
    assert rows["custom:ClinVar"]["status"] == "OK" and rows["custom:ClinVar"]["controls_annotated"] > 0
    cfg["custom"]["ClinVar"]["enabled"] = False
    rows = {r["plugin"]: r for r in smoke.run_smoke(cfg, plugins + [cvp], region=REGION)}
    assert rows["Cvp"]["status"] == "FAIL" and "custom.ClinVar" in rows["Cvp"]["message"]
    assert "custom:ClinVar" not in rows


def test_cli_smoke_prints_and_exits_1_on_fail(minimal_cfg, monkeypatch, capsys):
    row = {"plugin": "Toy", "status": "FAIL", "header_ok": True, "controls_n": 3, "controls_annotated": 0,
           "eligible_n": 0, "annotated_n": 0, "coverage": None, "message": "none annotated"}
    unchecked = {**row, "plugin": "Nc", "status": "unchecked", "message": "no control defined in the registry"}
    monkeypatch.setattr(cli, "_plugins", lambda cfg: [])
    monkeypatch.setattr(cli.smk, "run_smoke", lambda *a, **k: [row, unchecked])
    assert cli.HANDLERS["smoke"](minimal_cfg, argparse.Namespace(seed=7)) == 1
    out = capsys.readouterr().out
    assert "Toy" in out and "FAIL" in out and "1 failed (Toy)" in out
    assert "1 unchecked (Nc)" in out                                      # S1(a): unchecked rows are counted
    monkeypatch.setattr(cli.smk, "run_smoke", lambda *a, **k: [{**row, "status": "OK"}])
    assert cli.HANDLERS["smoke"](minimal_cfg, argparse.Namespace(seed=7)) == 0


def test_cli_smoke_subcommand_parses_seed():
    assert cli.build_parser().parse_args(["smoke", "--seed", "3"]).seed == 3


def test_plugin_controls_read_every_file_of_a_glob_directory(minimal_cfg, tmp_path):
    from vep_pipeline.registry import FileSpec
    d = tmp_path / "snv"
    d.mkdir()
    for contig, pos in (("chr17", 200), ("chr1", 200)):            # per-chromosome files, header line skipped (-S 1)
        data = d / f"{contig}.tsv"
        data.write_text(f"chrom\tpos\tref\talt\tstrand\n{contig}\t{pos}\t{SEQ[pos - 1]}\t{_alt(pos)}\t+\n")
        _bgzip_tabix(data, ["-s", "1", "-b", "2", "-e", "2", "-S", "1"])
    (d / "notes.tsv").write_text("ignored\n")
    spec = PluginSpec(name="Dir", module="Dir", category="manual", output_fields=("DIR",),
                      files=(FileSpec(key="snv", dest="Dir/snv", index="tbi", glob="chr*.tsv.gz"),),
                      control={"source": "data", "cols": [1, 2, 3, 4]})
    rp = ResolvedPlugin(spec, {"snv": d}, {}, "Dir")
    assert smoke._plugin_controls(rp, "17", 100, 900, 7, minimal_cfg["tools"]["env"]) == [("17", 200, SEQ[199], _alt(200))]
