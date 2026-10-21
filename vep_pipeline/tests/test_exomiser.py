import argparse
import re

import pytest
import yaml
from prioritise_helpers import write_genes_tsv, write_variants_tsv

from vep_pipeline import cli
from vep_pipeline.config import ConfigError
from vep_pipeline.exomiser import load_pheno, load_remm, phenopacket_text, pheno_strength, render_exomiser
from vep_pipeline.rank_settings import exomiser_settings


@pytest.fixture
def exo_cfg(minimal_cfg, tmp_path):
    minimal_cfg["exomiser"] = {"hpo": ["HP:0001250", "HP:0001263"], "sex": "MALE"}
    (tmp_path / "out").mkdir()
    (tmp_path / "out" / "S1.input.vcf.gz").write_bytes(b"x")
    return minimal_cfg


def test_phenopacket_text():
    doc = yaml.safe_load(phenopacket_text("S1", ["HP:0001250"], "MALE", "2026-10-17T00:00:00Z"))
    assert doc["subject"] == {"id": "S1", "sex": "MALE"}
    assert doc["phenotypicFeatures"] == [{"type": {"id": "HP:0001250"}}]
    assert doc["metaData"]["phenopacketSchemaVersion"] == 1.0


def test_render_exomiser_writes_phenopacket_and_job(exo_cfg, tmp_path):
    script = render_exomiser(exo_cfg)
    text = script.read_text()
    res = tmp_path / "out" / "exomiser"
    assert "--preset phenotype_only" in text and "#SBATCH --mem=48G" in text
    assert f"run-exomiser.sh S1 {tmp_path}/out/S1.input.vcf.gz {res}/S1.phenopacket.yml {res}" in text
    assert not re.search(r"^set -\w*[eu]", text, re.M)   # exomiser-env.sh must not be sourced under -e/-u
    assert "exec " not in text                        # exec would bypass the EXIT trap that removes $TMPDIR
    assert re.search(r"run-exomiser\.sh .*; exit \$\?\n$", text)
    assert yaml.safe_load((res / "S1.phenopacket.yml").read_text())["subject"]["id"] == "S1"


def test_phenotype_pass_uses_an_empty_grch38_vcf(exo_cfg, tmp_path):
    # Without a VCF, Exomiser 15.1 defaults the sample assembly to hg19 and aborts (job 3713229); `--assembly`
    # alone is a CLI usage error. A header-only VCF with the proband column makes the assembly explicit.
    text = render_exomiser(exo_cfg).read_text()
    res = tmp_path / "out" / "exomiser"
    empty = (res / "phenotype_only.vcf").read_text().splitlines()
    assert empty[0] == "##fileformat=VCFv4.2" and empty[-1].split("\t")[-2:] == ["FORMAT", "S1"]
    assert not [l for l in empty if not l.startswith("#")]           # no records
    pheno_cmd = text.split("run-exomiser.sh")[0]
    assert f"--preset phenotype_only --sample {res}/S1.phenopacket.yml --vcf {res}/phenotype_only.vcf " \
           "--assembly GRCh38" in pheno_cmd


def test_render_needs_the_normalised_input(minimal_cfg):
    minimal_cfg["exomiser"] = {"hpo": ["HP:0001250"]}
    with pytest.raises(ConfigError, match="input.vcf.gz"):
        render_exomiser(minimal_cfg)


def test_given_phenopacket_must_match_sample(exo_cfg, tmp_path):
    p = tmp_path / "p.yml"
    p.write_text("subject: {id: OTHER}\n")
    exo_cfg["exomiser"]["phenopacket"] = str(p)
    with pytest.raises(ConfigError, match="subject.id"):
        render_exomiser(exo_cfg)


def test_load_pheno_ranks_with_ties(exo_cfg, tmp_path):
    res = tmp_path / "out" / "exomiser"
    write_genes_tsv(res / "phenotype" / "S1.genes.tsv",
                    [(1, "A", "ANY", 0.9, 0.9), (2, "B", "ANY", 0.8, 0.8), (3, "C", "ANY", 0.8, 0.8),
                     (4, "D", "ANY", 0.5, 0.5)])
    write_genes_tsv(res / "exome" / "S1.genes.tsv", [(3, "A", "AR", 0.7, 0.9), (5, "A", "AD", 0.4, 0.6)])
    write_genes_tsv(res / "genome" / "S1.genes.tsv", [(1, "A", "AR", 0.95, 0.9)])
    pe = load_pheno(exo_cfg, exomiser_settings(exo_cfg))
    assert pe.rank == {"A": 1, "B": 2, "C": 2, "D": 4}
    assert pe.combined["A"] == 0.95 and pe.exomiser_rank["A"] == 1
    cuts = {"strong": 1, "moderate": 2}
    assert [pheno_strength(pe, g, cuts) for g in ("A", "B", "D", "Z")] == ["strong", "moderate", "weak", "unknown"]
    assert pheno_strength(None, "A", cuts) == "none"


def test_missing_results_say_how_to_get_them(exo_cfg):
    with pytest.raises(ConfigError, match="exomiser --submit"):
        load_pheno(exo_cfg, exomiser_settings(exo_cfg))


def test_cli_prints_sbatch_without_submit(exo_cfg, capsys):
    assert cli.exomiser_cmd(exo_cfg, argparse.Namespace(submit=False)) == 0
    assert "sbatch --parsable" in capsys.readouterr().out


def test_load_remm_parses_all_path_and_takes_the_max(exo_cfg, tmp_path):
    ex = exomiser_settings(exo_cfg)
    write_variants_tsv(tmp_path / "out" / "exomiser" / "genome" / "S1.variants.tsv", [
        ("1", 100, "CA", "C", "CADD=0.99241424,REMM=0.993508"),
        ("1", 100, "CA", "C", "CADD=0.5,REMM=0.970000", "AR"),          # same variant, another mode of inheritance
        ("1", 200, "A", "G", "CADD=0.9"),                               # no REMM token
        ("1", 300, "T", "C", "."),
        ("X", 400, "G", "A", ""),
        ("X", 500, "G", "A", "REMM=0.5"),
    ])
    assert load_remm(exo_cfg, ex) == {"1:100:CA>C": 0.993508, "X:500:G>A": 0.5}


def test_load_remm_missing_file_and_columns(exo_cfg, tmp_path):
    ex = exomiser_settings(exo_cfg)
    with pytest.raises(ConfigError, match="exomiser --submit"):
        load_remm(exo_cfg, ex)
    f = tmp_path / "out" / "exomiser" / "genome" / "S1.variants.tsv"
    f.parent.mkdir(parents=True)
    f.write_text("CONTIG\tSTART\tREF\tALT\n1\t1\tA\tG\n")
    with pytest.raises(ConfigError, match="ALL_PATH"):
        load_remm(exo_cfg, ex)
