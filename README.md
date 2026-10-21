# MVA Hackathon 2026, Track 1: variant prioritisation pipeline (sbgonenc96)

Code for ranking candidate causal variants in one WGS proband with Mosaic Variegated Aneuploidy.

## Contents

| Path | What it is |
|---|---|
| `vep_pipeline/` | The main pipeline (Models 2 and 3), a Python CLI on SLURM: Ensembl VEP 116 annotation, Exomiser/Genomiser phenotype ranking, WhatsHap read-backed phasing, and a rule-based top-10 ranking with compound-heterozygous pairs. See [vep_pipeline/README.md](vep_pipeline/README.md). |
| `exomiser/` | Install and run scripts for Exomiser 15.1.0 + Genomiser (data release 2602): `bin/` (environment, data fetch, VCF preparation, two-pass run), `etc/` (application properties, Genomiser analysis). See [exomiser/README.md](exomiser/README.md). |
| `vep_plugins_local/` | The SpliceAI2 VEP plugin the pipeline uses, and its `plugin_config.txt` entry. See [vep_plugins_local/README.md](vep_plugins_local/README.md). |
| `mva_exomiser_postprocess.py` | marimo notebook that turned Exomiser/Genomiser variant tables into the challenge CSV (Model 1). Input paths and the sample ID are placeholders. |
| `annotate_vcf.sh` | Helper that annotates a VCF against a GFF and a known-variants VCF (`-h` for usage). |

## Not included

- **Patient data and results.** These are never committed: the proband's VCF/BAM, phenotype (HPO terms), run outputs, run reports and the challenge methods forms. `vep_pipeline/config.yaml` holds placeholders for the sample, input VCF, HPO terms and sex.
- **Third-party code.** Fetch these separately at the versions used:
  - [Ensembl VEP](https://github.com/Ensembl/ensembl-vep) `release/116` (116.3, commit f479099e);
  - [VEP_plugins](https://github.com/Ensembl/VEP_plugins) `release/116` (commit 59b34e5), plus `vep_plugins_local/SpliceAI2.pm`;
  - [Illumina SpliceAI2](https://github.com/Illumina/SpliceAI2) (v2.0-2-g1dcd2a1; precomputed scores v2.0, non-commercial terms);
  - the Exomiser 15.1.0 CLI release and its data (see `exomiser/README.md`).

## Tests

`cd vep_pipeline && micromamba run -n vep116 pytest tests -m "not integration"`. Some tests use VEP's own test cache and plugins, so they expect the `ensembl-vep/` and `VEP_plugins/` checkouts next to `vep_pipeline/` (as listed above).
