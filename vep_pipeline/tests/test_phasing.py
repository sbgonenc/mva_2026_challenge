import pytest

from vep_pipeline.config import ConfigError
from vep_pipeline.phasing import load_readback, readback_path
from vep_pipeline.rank_settings import phasing_settings


def test_load_readback_parses_and_skips_unphased(tmp_path):
    p = tmp_path / "S1.readback.tsv"
    p.write_text("key\tgt\tps\n1:10:A>G\t0|1\t10\n1:20:C>T\t1|0\t10\n1:30:A>C\t0/1\t.\n1:40:A>T\t0|1\t.\n"
                 "1:50:A>T\t0|1\t\n")
    assert load_readback(p) == {"1:10:A>G": ("0|1", "10"), "1:20:C>T": ("1|0", "10")}


def test_load_readback_needs_the_exact_header(tmp_path):
    p = tmp_path / "S1.readback.tsv"
    p.write_text("variant\tgt\tps\n1:10:A>G\t0|1\t10\n")
    with pytest.raises(ConfigError, match="key.*gt.*ps"):
        load_readback(p)


@pytest.mark.parametrize("row", ["1:20:C>T\t1|", "1:20:C>T", "1:20:C>T\t1|0\t10\t5", ""])
def test_load_readback_rejects_a_row_without_three_fields(tmp_path, row):
    p = tmp_path / "S1.readback.tsv"
    p.write_text(f"key\tgt\tps\n1:10:A>G\t0|1\t10\n{row}\n")                  # e.g. a truncated (partial) sidecar
    with pytest.raises(ConfigError, match=r"S1\.readback\.tsv line 3: .*3 tab-separated fields"):
        load_readback(p)


def test_readback_path_follows_results_dir(minimal_cfg):
    ph = phasing_settings(minimal_cfg)
    assert str(readback_path(minimal_cfg, ph)).endswith("/out/phasing/S1.readback.tsv")


# ----------------------------------------------------------------------------- phase command

import subprocess
from pathlib import Path

import yaml
from prioritise_helpers import make_call, make_cand, write_vcf as _write_vcf

from vep_pipeline import cli, jobs, phasing
from vep_pipeline.candidates import Allele
from vep_pipeline.phasing import contig_map, phase_regions, write_regions_bed

CSQ = ["Allele", "Consequence", "SYMBOL", "Feature", "MANE_SELECT"]


def write_vcf(path, records, csq_fields, **kw):
    return _write_vcf(path, [{"csq": [], **r} for r in records], csq_fields, **kw)


def allele(pos, zyg="het", chrom="1", ref="A", alt="G", symbol="G", pieces=None):
    pcs = pieces or [make_cand(pos, symbol=symbol, call=make_call(pos, chrom=chrom, ref=ref, alt=alt))]
    return Allele(tuple(pcs), "strong", False, (), zyg, False)


def test_phase_regions_skips_genes_with_one_het_and_applies_pad():
    by_gene = {"ONE": [allele(5000), allele(6000, zyg="hom")],
               "TWO": [allele(10000), allele(10500, ref="ACG")]}
    assert phase_regions(by_gene, 100) == [("1", 9900, 10602)]


def test_phase_regions_clamps_at_contig_start():
    assert phase_regions({"G": [allele(10), allele(50)]}, 100) == [("1", 1, 150)]


def test_phase_regions_merges_overlapping_and_touching_per_contig():
    by_gene = {"A": [allele(1000), allele(2000)],
               "B": [allele(2100), allele(3000)],             # overlaps A after the pad
               "C": [allele(3201), allele(4000)],             # touches B (3100 -> 3101..): merged
               "D": [allele(9000, chrom="2"), allele(9500, chrom="2")],
               "E": [allele(5000), allele(5100)]}             # far from the first block: its own region
    assert phase_regions(by_gene, 100) == [("1", 900, 4100), ("1", 4900, 5200), ("2", 8900, 9600)]


def test_phase_regions_spans_a_multi_piece_allele():
    mnv = allele(0, pieces=[make_cand(2000, call=make_call(2000)), make_cand(2005, call=make_call(2005, ref="AT"))])
    assert phase_regions({"G": [allele(1000), mnv]}, 50) == [("1", 950, 2056)]


def test_write_regions_bed_is_zero_based_half_open(tmp_path):
    p = tmp_path / "regions.bed"
    write_regions_bed([("1", 900, 1100), ("2", 1, 10)], p)
    assert p.read_text() == "1\t899\t1100\n2\t0\t10\n"


@pytest.mark.parametrize("vcf, bam, expected", [
    (["1", "2"], ["1", "2", "chr1"], {"1": "1", "2": "2"}),                      # identity wins
    (["1", "X"], ["chr1", "chrX"], {"1": "chr1", "X": "chrX"}),                  # add chr
    (["chr1", "chr2"], ["1", "2"], {"chr1": "1", "chr2": "2"}),                  # strip chr
    (["MT"], ["chrM"], {"MT": "chrM"}),
    (["chrM"], ["MT"], {"chrM": "MT"}),
])
def test_contig_map(vcf, bam, expected):
    assert contig_map(bam, vcf) == expected


def test_contig_map_names_the_unmappable_contig():
    with pytest.raises(ConfigError, match="'7'"):
        contig_map(["chr1"], ["1", "7"])


def phasing_cfg(cfg, tmp_path, **extra):
    bam, ref = tmp_path / "s.bam", tmp_path / "ref.fa"
    for f in (bam, tmp_path / "s.bam.bai", ref, tmp_path / "ref.fa.fai"):
        f.write_text("")
    cfg["phasing"] = {"enabled": True, "bam": str(bam), "reference": str(ref), "env": "wh_env", **extra}
    return cfg


def test_render_phase_writes_the_job_script_and_frozen_config(minimal_cfg, tmp_path):
    cfg = phasing_cfg(minimal_cfg, tmp_path, mem="7G", time="03:00:00")
    cfg["output"]["dir"] = str(tmp_path / "odd dir")             # a space: every path must be quoted
    ph = phasing_settings(cfg)
    script = phasing.render_phase(cfg, ph)
    res = Path(ph["results_dir"])
    text = script.read_text()
    assert script == res / "phase.sbatch"
    assert yaml.safe_load((res / "phase.config.yaml").read_text())["sample"] == "S1"
    assert "#SBATCH --mem=7G" in text and "#SBATCH --time=03:00:00" in text
    assert "#SBATCH --cpus-per-task=4" in text and "phase_%j.out" in text
    assert "set -euo pipefail" in text and "export PYTHONPATH=" in text
    q = lambda p: "'" + str(p) + "'"
    assert f"micromamba run -n wh_env samtools view -H {ph['bam']} > {q(res / 'bam_header.sam')}" in text
    assert f"python -m vep_pipeline.cli phase --prepare -c {q(res / 'phase.config.yaml')}" in text
    assert "micromamba run -n vep116 python -m vep_pipeline.cli phase --prepare" in text
    assert (f"micromamba run -n wh_env whatshap phase --reference {ph['reference']} --ignore-read-groups --indels "
            f"-o {q(res / 'phased.vcf.gz')} {q(res / 'subset.vcf.gz')} {ph['bam']}") in text
    assert f"micromamba run -n vep116 bcftools index -t {q(res / 'phased.vcf.gz')}" in text
    assert f"python -m vep_pipeline.cli phase --finish -c {q(res / 'phase.config.yaml')}" in text
    assert text.index("samtools view") < text.index("--prepare") < text.index("whatshap phase") \
        < text.index("bcftools index") < text.index("--finish")
    assert f"[[ -s {q(res / 'regions.bed')} ]]" in text           # nothing to phase: skip WhatsHap
    rm = (f"rm -f -- {q(res / 'S1.readback.tsv')} {q(res / 'phased.vcf.gz')} "
          f"{q(res / 'phased.vcf.gz.tbi')}\n")                       # an earlier run's outputs, before any step
    assert rm in text and text.index("export PYTHONPATH=") < text.index(rm) < text.index("samtools view")


def test_render_phase_omits_indels_when_disabled(minimal_cfg, tmp_path):
    cfg = phasing_cfg(minimal_cfg, tmp_path, indels=False)
    text = phasing.render_phase(cfg, phasing_settings(cfg)).read_text()
    assert "--indels" not in text and "--ignore-read-groups" in text


def test_finish_writes_the_sidecar_with_contigs_mapped_back(minimal_cfg, tmp_path, capsys):
    cfg = phasing_cfg(minimal_cfg, tmp_path)
    ph = phasing_settings(cfg)
    res = Path(ph["results_dir"])
    res.mkdir(parents=True)
    write_vcf(res / "phased.vcf.gz", [
        {"chrom": "chr1", "pos": 100, "GT": "0|1", "PS": "100"},
        {"chrom": "chr1", "pos": 200, "ref": "AC", "alt": "A", "GT": "1|0", "PS": "100"},
        {"chrom": "chr1", "pos": 300, "GT": "0/1", "PS": "."},                 # unphased: not in the sidecar
        {"chrom": "chr1", "pos": 400, "GT": "0|1", "PS": "400"},
        {"chrom": "chrX", "pos": 500, "GT": "1|0", "PS": "500"},
    ], CSQ, fmt=("GT", "PS"))
    (res / "contigs_from_bam.txt").write_text("chr1\t1\nchrX\tX\n")
    (res / "region_alleles.tsv").write_text("G\t1:100:A>G\nG\t1:200:AC>A\nH\t1:200:AC>A\nH\t1:400:A>G\nH\t1:300:A>G\n")
    assert phasing.finish(cfg, ph) == readback_path(cfg, ph)
    assert readback_path(cfg, ph).read_text() == ("key\tgt\tps\n1:100:A>G\t0|1\t100\n1:200:AC>A\t1|0\t100\n"
                                                  "1:400:A>G\t0|1\t400\nX:500:A>G\t1|0\t500\n")
    assert load_readback(readback_path(cfg, ph))["X:500:A>G"] == ("1|0", "500")
    out = capsys.readouterr().out
    assert "4 phased records" in out and "3 blocks" in out
    assert "1 candidate pair" in out                    # G: 100 + 200 share block 100; H: 200/400 do not, 300 is unphased


def test_finish_never_leaves_a_partial_sidecar(minimal_cfg, tmp_path, monkeypatch):
    cfg = phasing_cfg(minimal_cfg, tmp_path)
    ph = phasing_settings(cfg)
    res = Path(ph["results_dir"])
    res.mkdir(parents=True)
    write_vcf(res / "phased.vcf.gz", [{"chrom": "chr1", "pos": 100, "GT": "0|1", "PS": "100"},
                                      {"chrom": "chr1", "pos": 200, "GT": "1|0", "PS": "100"}], CSQ, fmt=("GT", "PS"))
    (res / "contigs_from_bam.txt").write_text("chr1\t1\n")
    (res / "region_alleles.tsv").write_text("G\t1:100:A>G\nG\t1:200:A>G\n")
    real = Path.write_text

    def disk_full(self, text, *a, **kw):                     # half the text reaches the disk, then the write fails
        real(self, text[:len(text) // 2], *a, **kw)
        raise OSError("No space left on device")
    monkeypatch.setattr(Path, "write_text", disk_full)
    with pytest.raises(OSError, match="No space"):
        phasing.finish(cfg, ph)
    assert not readback_path(cfg, ph).exists()
    assert sorted(p.name for p in res.iterdir()) == ["contigs_from_bam.txt", "phased.vcf.gz", "region_alleles.tsv"]
    monkeypatch.undo()
    phasing.finish(cfg, ph)
    assert load_readback(readback_path(cfg, ph)) == {"1:100:A>G": ("0|1", "100"), "1:200:A>G": ("1|0", "100")}
    assert sorted(p.name for p in res.iterdir()) == ["S1.readback.tsv", "contigs_from_bam.txt", "phased.vcf.gz",
                                                    "region_alleles.tsv"]


def test_finish_needs_the_phased_vcf(minimal_cfg, tmp_path):
    cfg = phasing_cfg(minimal_cfg, tmp_path)
    with pytest.raises(ConfigError, match="phased.vcf.gz"):
        phasing.finish(cfg, phasing_settings(cfg))


def prepare_setup(cfg, tmp_path, vep_records, input_records, bam_contigs=("chr1", "chrX")):
    """Fake VEP output + normalised input VCF, g2p table, BAM header; prioritisation reads them as in production."""
    cfg = phasing_cfg(cfg, tmp_path)
    cfg["exomiser"] = {"enabled": False}
    g2p = tmp_path / "g2p.csv"
    g2p.write_text('"gene symbol","allelic requirement","confidence","molecular mechanism"\n'
                   "GENEAR,biallelic_autosomal,definitive,loss of function\n")
    cfg["prioritise"] = {"g2p": {"file": str(g2p)}}
    out = Path(cfg["output"]["dir"])
    out.mkdir(parents=True, exist_ok=True)
    write_vcf(out / "S1.vep.vcf.gz", vep_records, CSQ)
    inp = write_vcf(out / "S1.input.vcf.gz", input_records, CSQ)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(inp)], check=True)
    res = Path(phasing_settings(cfg)["results_dir"])
    res.mkdir(parents=True, exist_ok=True)
    (res / "bam_header.sam").write_text("@HD\tVN:1.6\n" + "".join(f"@SQ\tSN:{c}\tLN:1000000\n" for c in bam_contigs))
    return cfg, res


def ge(symbol="GENEAR", **kw):
    return {"SYMBOL": symbol, "Consequence": "stop_gained", "Feature": f"ENST_{symbol}", "MANE_SELECT": "NM_1", **kw}


def test_prepare_writes_regions_subset_and_maps(minimal_cfg, tmp_path, capsys):
    vep = [{"pos": 1000, "csq": [ge()]}, {"pos": 5000, "csq": [ge()]}, {"pos": 20000, "csq": [ge("OTHER")]}]
    inp = [{"pos": 800, "GT": "0/1"}, {"pos": 1000, "GT": "0/1"}, {"pos": 3000, "GT": "0/1"},
           {"pos": 3100, "GT": "1/1"}, {"pos": 5000, "GT": "0/1"}, {"pos": 5050, "GT": "0/1", "filter": "LowQual"},
           {"pos": 20000, "GT": "0/1"}]
    cfg, res = prepare_setup(minimal_cfg, tmp_path, vep, inp)
    phasing.prepare(cfg, phasing_settings(cfg))
    assert (res / "regions.bed").read_text() == "1\t0\t6000\n"                  # 1000-1000 .. 5000+1000, 1-based start 1
    assert (res / "region_alleles.tsv").read_text() == "GENEAR\t1:1000:A>G\nGENEAR\t1:5000:A>G\n"
    assert (res / "contigs_to_bam.txt").read_text() == "1\tchr1\n"
    assert (res / "contigs_from_bam.txt").read_text() == "chr1\t1\n"
    body = subprocess.run(["bcftools", "query", "-f", "%CHROM:%POS\n", str(res / "subset.vcf.gz")], check=True,
                          capture_output=True, text=True).stdout.split()
    assert body == ["chr1:800", "chr1:1000", "chr1:3000", "chr1:5000"]         # het + PASS only, in the regions
    assert (res / "subset.vcf.gz.tbi").is_file()
    assert "1 region" in capsys.readouterr().out


def test_prepare_without_candidate_pairs_leaves_an_empty_sidecar(minimal_cfg, tmp_path, capsys):
    cfg, res = prepare_setup(minimal_cfg, tmp_path, [{"pos": 1000, "csq": [ge()]}], [{"pos": 1000, "GT": "0/1"}])
    (res / "regions.bed").write_text("1\t0\t5\n")                                # stale file of an earlier run
    ph = phasing_settings(cfg)
    phasing.prepare(cfg, ph)
    assert (res / "regions.bed").read_text() == ""
    assert readback_path(cfg, ph).read_text() == "key\tgt\tps\n"
    assert "nothing to phase" in capsys.readouterr().out


STALE = ("S1.readback.tsv", "phased.vcf.gz", "phased.vcf.gz.tbi")


def _stale(res):
    """Outputs of an earlier phase job: a later prioritise must never read them as this run's."""
    (res / "S1.readback.tsv").write_text("key\tgt\tps\n1:1000:A>G\t0|1\t7\n1:5000:A>G\t1|0\t7\n")
    for name in STALE[1:]:
        (res / name).write_text("old")


def test_prepare_removes_the_outputs_of_an_earlier_run(minimal_cfg, tmp_path):
    vep = [{"pos": 1000, "csq": [ge()]}, {"pos": 5000, "csq": [ge()]}]
    cfg, res = prepare_setup(minimal_cfg, tmp_path, vep, [{"pos": 1000, "GT": "0/1"}, {"pos": 5000, "GT": "0/1"}])
    _stale(res)
    phasing.prepare(cfg, phasing_settings(cfg))                                   # regions to phase: no sidecar yet
    assert not any((res / name).exists() for name in STALE)
    assert (res / "subset.vcf.gz").is_file()


@pytest.mark.parametrize("failure", ["bam_header", "contig"])
def test_a_failed_prepare_leaves_no_stale_sidecar(minimal_cfg, tmp_path, failure):
    vep = [{"pos": 1000, "csq": [ge()]}, {"pos": 5000, "csq": [ge()]}]
    cfg, res = prepare_setup(minimal_cfg, tmp_path, vep, [{"pos": 1000}, {"pos": 5000}],
                             bam_contigs=("chr2",) if failure == "contig" else ("chr1",))
    if failure == "bam_header":
        (res / "bam_header.sam").unlink()
    _stale(res)
    with pytest.raises(ConfigError):
        phasing.prepare(cfg, phasing_settings(cfg))
    assert not any((res / name).exists() for name in STALE)


def test_prepare_fails_on_a_contig_the_bam_lacks(minimal_cfg, tmp_path):
    vep = [{"pos": 1000, "csq": [ge()]}, {"pos": 5000, "csq": [ge()]}]
    cfg, _ = prepare_setup(minimal_cfg, tmp_path, vep, [{"pos": 1000}, {"pos": 5000}], bam_contigs=("chr2",))
    with pytest.raises(ConfigError, match="'1'"):
        phasing.prepare(cfg, phasing_settings(cfg))


def test_prepare_needs_the_bam_header(minimal_cfg, tmp_path):
    vep = [{"pos": 1000, "csq": [ge()]}, {"pos": 5000, "csq": [ge()]}]
    cfg, res = prepare_setup(minimal_cfg, tmp_path, vep, [{"pos": 1000}, {"pos": 5000}])
    (res / "bam_header.sam").unlink()
    with pytest.raises(ConfigError, match="bam_header.sam"):
        phasing.prepare(cfg, phasing_settings(cfg))


# ----------------------------------------------------------------------------- CLI

def write_cfg(cfg, tmp_path):
    p = tmp_path / "cfg.yaml"
    p.write_text(yaml.safe_dump(cfg))
    return str(p)


def test_cli_phase_prints_the_sbatch_command(minimal_cfg, tmp_path, capsys, monkeypatch):
    cfg = phasing_cfg(minimal_cfg, tmp_path)
    monkeypatch.setattr(jobs, "submit_one", lambda s: pytest.fail("must not submit"))
    assert cli.main(["phase", "-c", write_cfg(cfg, tmp_path)]) == 0
    out = capsys.readouterr().out
    assert "sbatch --parsable " in out and "phase.sbatch" in out and "not submitted: add --submit" in out


def test_cli_phase_submit(minimal_cfg, tmp_path, capsys, monkeypatch):
    cfg = phasing_cfg(minimal_cfg, tmp_path)
    monkeypatch.setattr(jobs, "submit_one", lambda s: "4242")
    assert cli.main(["phase", "-c", write_cfg(cfg, tmp_path), "--submit"]) == 0
    assert "4242" in capsys.readouterr().out


@pytest.mark.parametrize("flags", [[], ["--submit"], ["--prepare"], ["--finish"]])
def test_cli_phase_needs_phasing_enabled(minimal_cfg, tmp_path, capsys, flags):
    assert cli.main(["phase", "-c", write_cfg(minimal_cfg, tmp_path), *flags]) == 2
    assert "phasing.enabled" in capsys.readouterr().err


def test_cli_phase_modes_are_mutually_exclusive(minimal_cfg, tmp_path):
    with pytest.raises(SystemExit):
        cli.main(["phase", "-c", write_cfg(minimal_cfg, tmp_path), "--submit", "--finish"])


def test_cli_phase_finish_dispatches(minimal_cfg, tmp_path, capsys):
    cfg = phasing_cfg(minimal_cfg, tmp_path)
    assert cli.main(["phase", "-c", write_cfg(cfg, tmp_path), "--finish"]) == 2     # no phased.vcf.gz yet
    assert "phased.vcf.gz" in capsys.readouterr().err
