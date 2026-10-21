"""Exomiser for `prioritise`: phenopacket, one SLURM job (a phenotype-only pass, then the site's two-pass
run-exomiser.sh), and reading gene results back.

The job runs the install at `exomiser.root` (default /fast/AG_Bunina/Berk/tools/exomiser); its bin/run-exomiser.sh
reads the analysis YAMLs from that tree. The phenotype-only pass scores every gene independently of variant filters,
so a gene the genome pass drops (phenotype < 0.501, frequency, inheritance) still gets a phenotype rank.
"""
from __future__ import annotations
import csv
import datetime
import re
import shlex
from dataclasses import dataclass
from pathlib import Path

import yaml

from .config import ConfigError
from .jobs import _key, _out_dir, _prefix
from .rank_settings import exomiser_settings

# The phenotype-only pass needs an explicit GRCh38 assembly: without a VCF, Exomiser 15.1 defaults the sample to hg19
# and aborts with UnsupportedGenomeAssemblyException (job 3713229; also the site's gate-1 test), and `--assembly` is
# only accepted together with `--vcf`. A header-only VCF carrying the proband column satisfies both; the pass then
# scores every gene (19,193 rows on data release 2602).
HPO_VERSION = "2026-01-08"            # hp.obo of Exomiser data release 2602 (exomiser/mva_patient.phenopacket.yaml)
GENE_COLS = ("#RANK", "GENE_SYMBOL", "EXOMISER_GENE_PHENO_SCORE", "EXOMISER_GENE_COMBINED_SCORE")


def phenopacket_text(sample, hpo, sex, created) -> str:
    doc = {"id": f"{sample}_vep_pipeline",
           "subject": {"id": sample, "sex": sex},
           "phenotypicFeatures": [{"type": {"id": t}} for t in hpo],
           "metaData": {"created": created, "createdBy": "vep_pipeline",
                        "resources": [{"id": "hp", "name": "human phenotype ontology",
                                       "url": "http://purl.obolibrary.org/obo/hp.owl", "version": HPO_VERSION,
                                       "namespacePrefix": "HP", "iriPrefix": "http://purl.obolibrary.org/obo/HP_"}],
                        "phenopacketSchemaVersion": 1.0}}
    return yaml.safe_dump(doc, sort_keys=False)


def empty_vcf_text(sample) -> str:
    """Header-only VCF with the proband column, for the phenotype-only pass (see the note above)."""
    return "##fileformat=VCFv4.2\n" + "\t".join(
        ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample]) + "\n"


def job_script(cfg, ex, vcf: Path, ppkt: Path, empty_vcf: Path) -> str:
    q, s, sample = shlex.quote, cfg["slurm"], cfg["sample"]
    res = Path(ex["results_dir"])
    return f"""#!/usr/bin/env bash
#SBATCH --job-name=exomiser
#SBATCH --partition={s['partition']}
#SBATCH --time={ex['time']}
#SBATCH --mem={ex['mem']}
#SBATCH --cpus-per-task=4
#SBATCH --output={s['log_dir']}/exomiser_%j.out
#SBATCH --export=ALL
# No errexit/nounset: exomiser-env.sh must be sourced without them (see its header).
set -o pipefail
export TMPDIR="${{SLURM_TMPDIR:-/tmp}}/exomiser.${{SLURM_JOB_ID:-$$}}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
source {q(ex['root'] + '/bin/exomiser-env.sh')} || exit 1
mkdir -p {q(str(res / 'phenotype'))}
java -Xmx"$(exo_heap_gib exome)g" $(exo_jvm_flags {q(sample + '.phenotype')}) -jar "$EXO_JAR" analyse \\
  --preset phenotype_only --sample {q(str(ppkt))} --vcf {q(str(empty_vcf))} --assembly GRCh38 \\
  --output-directory {q(str(res / 'phenotype'))} --output-filename {q(sample)} --output-format TSV_GENE || exit 1
{q(ex['root'] + '/bin/run-exomiser.sh')} {q(sample)} {q(str(vcf))} {q(str(ppkt))} {q(str(res))}; exit $?
"""


def render_exomiser(cfg) -> Path:
    """Write the phenopacket (unless exomiser.phenopacket is given), the header-only VCF for the phenotype-only pass
    and {results_dir}/exomiser.sbatch."""
    ex = exomiser_settings(cfg)
    if not ex["enabled"]:
        raise ConfigError("exomiser.enabled is false: nothing to run")
    vcf = _out_dir(cfg) / f"{_prefix(cfg)}.input.vcf.gz"
    if not vcf.is_file():
        raise ConfigError(f"{vcf} not found: run `vep-pipeline render` (or `run`) first; Exomiser uses the same "
                          "normalised input as VEP")
    res = Path(ex["results_dir"])
    res.mkdir(parents=True, exist_ok=True)
    Path(cfg["slurm"]["log_dir"]).mkdir(parents=True, exist_ok=True)
    if ex["phenopacket"]:
        ppkt = Path(ex["phenopacket"])
        subject = ((yaml.safe_load(ppkt.read_text()) or {}).get("subject") or {}).get("id")
        if subject != cfg["sample"]:
            raise ConfigError(f"{ppkt}: subject.id {subject!r} != sample {cfg['sample']!r}")
    else:
        ppkt = res / f"{cfg['sample']}.phenopacket.yml"
        created = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        ppkt.write_text(phenopacket_text(cfg["sample"], ex["hpo"], ex["sex"], created))
    empty_vcf = res / "phenotype_only.vcf"
    empty_vcf.write_text(empty_vcf_text(cfg["sample"]))
    script = res / "exomiser.sbatch"
    script.write_text(job_script(cfg, ex, vcf, ppkt, empty_vcf))
    return script


@dataclass(frozen=True)
class PhenoEvidence:
    rank: dict                # symbol -> phenotype-only rank (1 = best; ties share a rank)
    combined: dict            # symbol -> max EXOMISER_GENE_COMBINED_SCORE over the exome and genome passes
    exomiser_rank: dict       # symbol -> best #RANK over the exome and genome passes


def _read_genes(path: Path) -> list:
    if not path.is_file():
        raise ConfigError(f"Exomiser result missing: {path} (run `vep-pipeline exomiser --submit` and wait for "
                          "the job, or set exomiser.enabled: false)")
    with open(path, newline="") as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in GENE_COLS if c not in (rows.fieldnames or ())]
        if missing:
            raise ConfigError(f"{path}: columns {missing} missing")
        return list(rows)


def _competition_ranks(scores: dict) -> dict:
    ranks, prev, rank = {}, None, 0
    for i, (g, s) in enumerate(sorted(scores.items(), key=lambda kv: (-kv[1], kv[0])), 1):
        if s != prev:
            rank, prev = i, s
        ranks[g] = rank
    return ranks


def load_pheno(cfg, ex) -> PhenoEvidence:
    res, sample = Path(ex["results_dir"]), cfg["sample"]
    pheno = {}
    for r in _read_genes(res / "phenotype" / f"{sample}.genes.tsv"):
        g = r["GENE_SYMBOL"]
        pheno[g] = max(pheno.get(g, 0.0), float(r["EXOMISER_GENE_PHENO_SCORE"]))
    combined, ex_rank = {}, {}
    for sub in ("exome", "genome"):
        for r in _read_genes(res / sub / f"{sample}.genes.tsv"):
            g = r["GENE_SYMBOL"]
            combined[g] = max(combined.get(g, 0.0), float(r["EXOMISER_GENE_COMBINED_SCORE"]))
            ex_rank[g] = min(ex_rank.get(g, 10 ** 9), int(r["#RANK"]))
    return PhenoEvidence(_competition_ranks(pheno), combined, ex_rank)


VARIANT_COLS = ("CONTIG", "START", "REF", "ALT", "ALL_PATH")
_REMM = re.compile(r"(?:^|,)REMM=([0-9.eE+-]+)")


def load_remm(cfg, ex) -> dict:
    """{variant key: max Genomiser ReMM} from {results_dir}/genome/{sample}.variants.tsv (Exomiser TSV_VARIANT).

    ALL_PATH looks like `CADD=0.99,REMM=0.99` or is empty/'.'. A variant has one row per mode of inheritance, so
    the max is taken. Only variants the Genomiser pass kept appear here (its own filters apply). Smedley 2016 AJHG,
    doi:10.1016/j.ajhg.2016.07.005."""
    path = Path(ex["results_dir"]) / "genome" / f"{cfg['sample']}.variants.tsv"
    if not path.is_file():
        raise ConfigError(f"Exomiser result missing: {path} (run `vep-pipeline exomiser --submit` and wait for "
                          "the job, or set exomiser.enabled: false)")
    out = {}
    with open(path, newline="") as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in VARIANT_COLS if c not in (rows.fieldnames or ())]
        if missing:
            raise ConfigError(f"{path}: columns {missing} missing")
        for r in rows:
            m = _REMM.search(r["ALL_PATH"] or "")
            if m:
                k = _key(r["CONTIG"], r["START"], r["REF"], r["ALT"])
                out[k] = max(out.get(k, 0.0), float(m.group(1)))
    return out


def pheno_strength(pe, symbol, cuts) -> str:
    """'none' (Exomiser off), 'unknown' (gene not scored), else 'strong' / 'moderate' / 'weak' by rank."""
    if pe is None:
        return "none"
    r = pe.rank.get(symbol)
    if r is None:
        return "unknown"
    return "strong" if r <= cuts["strong"] else "moderate" if r <= cuts["moderate"] else "weak"
