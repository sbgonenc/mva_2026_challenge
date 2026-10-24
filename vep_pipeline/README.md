# vep_pipeline

Config-driven wrapper around Ensembl VEP 116 (offline, GRCh38) for one WGS proband. It annotates every variant with consequence, population, predictor, splicing and gene-level evidence, so that candidates, compound hets included, can be ranked against the patient's phenotypes. It does no allele-frequency or consequence filtering.

One YAML file (`config.yaml`) sets the VCF, the reference, the download location and the plugin set. The CLI downloads what is missing, checks it, smoke-tests each plugin on known positive controls, and runs VEP as one job per chromosome.

The design document, implementation plans and run reports are not in this public repository: they contain patient-level results.

## Quick start

```bash
cd /fast/AG_Bunina/Berk/projects/mva_challenge_26/scripts/vep_pipeline
# 1. edit config.yaml: input.vcf (required, see below), sample, output.dir, plugins
bin/vep-pipeline download --dry-run      # plan, sizes, what is reused or manual
bin/vep-pipeline download                # login node; resumable
bin/vep-pipeline check                   # files, ini keys, plugin load tests, VCF vs FASTA, sample, sex
bin/vep-pipeline smoke                   # positive-control run on BRCA1
bin/vep-pipeline run --slurm             # renders scripts, prints the sbatch commands
bin/vep-pipeline run --slurm --submit    # submits array job + dependent merge job
```

`bin/vep-pipeline` runs `python -m vep_pipeline.cli` inside the micromamba env `vep116`. Every subcommand takes `-c/--config` (default: `config.yaml` next to the package). The wrapper runs from your current directory, so a relative `-c` path (and any relative path inside the config) resolves from where you call it. Without `-c` it loads the shipped `config.yaml`, whatever your current directory.

### Rules that apply to every command

- **`input.vcf` is required by every command**, including `download`, `list-plugins` and `render`. Config loading fails without it. Set a path first; a placeholder is fine for `download`, which does not read the VCF. `check` additionally needs a real file for its VCF checks.
- **`render` and `run --slurm` normalise the VCF on the calling node**, even without `--submit`. Normalisation (`exomiser/bin/prepare-patient-vcf.sh`, in the `exomiser` env) plus the one-sample subset run right there, before any SLURM job exists. An existing normalised file `{output.dir}/{prefix}.norm.vcf.gz` is reused only when its fingerprint `{prefix}.norm.vcf.gz.fingerprint.json` matches exactly: input VCF real path, size and modification time; FASTA real path and size; `prepare-patient-vcf.sh` real path and SHA-256; and `sample`. Any difference re-normalises. The thread count (`-c`) does not change the output and is not part of the fingerprint. If login-node work is restricted, run these two commands from an interactive compute session (`srun --pty ...`).
- Paths in the config must not contain whitespace. `{a.b}` refers to another config key.
- `render` and `run` write under `output.dir`, and `smoke` under `output.dir/smoke`. All three also prepare the enabled custom tracks, so they can write the derived cCRE BED and the ClinVar symlink and index into `data.custom_data_dir`. `download` writes only under `data.root`. Nothing is written to `references/` or `exomiser/data/`.

## Commands

| Command | What it does |
|---|---|
| `list-plugins` | Table of all 83 registry plugins: category, default, enabled, reason when not enableable. |
| `download [--dry-run] [--only NAME...] [--skip-cache]` | Prints the plan (name, bytes, fetch/skip/MISSING, destination, manual instructions, total, free disk), then fetches. `--only` takes enabled plugin names plus `cache` and/or `custom`. |
| `check` | Preflight table (status, level, check, message); exits 1 if any error-level check fails. |
| `smoke` | Positive-control run through the normal chunk pipeline; writes `{output.dir}/smoke/smoke.tsv`; exits 1 on any FAIL. `--seed N` changes the control sampling (default 7). |
| `render` | Normalises and subsets the input, writes `vep.ini` and `run.sh` per chunk, `contigs.txt`, `merge.sh`; runs no VEP. |
| `run --local N` | Renders, runs N chunks in parallel on this machine, reconciles each, merges. |
| `run --slurm [--submit]` | Renders and writes `array.sbatch` and `merge.sbatch`. Without `--submit` it prints the two `sbatch` commands and stops. With `--submit` it submits the array, then the merge job with `--dependency=afterok:<array id> --kill-on-invalid-dep=yes`. |
| `run --merge-only` | Merges already reconciled chunks (this is what the SLURM merge job runs). |
| `reconcile --chunk NAME` | Reconciles one rendered chunk (the array task runs this after VEP). Chunk names are listed in `{output.dir}/work/contigs.txt`. |
| `exomiser [--submit]` | Writes the phenopacket from `exomiser.hpo` and the Exomiser SLURM job (phenotype-only, exome and genome passes). Without `--submit` it prints the `sbatch` command. Results go under `{output.dir}/exomiser`. |
| `phase [--submit]` | With `phasing.enabled`: writes the WhatsHap read-backed phasing job for the candidate genes (`{phasing.results_dir}/phase.sbatch`). Without `--submit` it prints the `sbatch` command. The job writes `{sample}.readback.tsv`, which `prioritise` reads; run it before `prioritise` and wait until it has finished. `--prepare` and `--finish` are the job's own steps. See "Read-backed phasing". |
| `regnet-build` | With the `regnet:` block (inputs only; `enabled` is not needed): compiles the scDoRI regulatory network and the S-LDSC results into small TSVs under `{regnet.build_dir}` (default `{output.dir}/regnet`) for `prioritise`. Run it inside an allocation. See "Regulatory network (`regnet`)". |
| `prioritise [--slurm [--submit]] [--top-n N] [-o CSV] [--sensitivity]` | Ranks findings in the final VEP VCF with the Exomiser phenotype evidence and writes the top-N CSV, `.candidates.tsv`, `.compound_het.tsv`, `.excluded.tsv`, (with `--sensitivity`) `.sensitivity.tsv` and (with `regnet.enabled`) `.regulatory.tsv` and `.regulatory.summary.tsv`. Use `--slurm` for a WGS VCF; see "Prioritising the top 10". |

`run` needs exactly one of `--local N`, `--slurm` or `--merge-only`; `--submit` only applies to `--slurm`.

### What `check` verifies

| Check | Level |
|---|---|
| `files`: every enabled plugin and custom file exists, is non-empty and indexed; cache `info.txt` exists | error |
| `ini_keys`: generated ini keys against the committed `vep_params.txt`; no inline comments | error |
| `plugin:<name>`: real VEP run on a 1-record VCF with `safe 1` and only that plugin; every declared output field must appear in the CSQ `Format:` | error |
| `spliceai_form`: SpliceAI arguments include `indel=` (see SpliceAI below) | error |
| `vcf_vs_fasta`: contig lengths against the `.fai` (after `chr` mapping), REF of the first 2,000 records against the FASTA, error above 2% mismatch | error |
| `sample`: `sample` is a column of the VCF | error |
| `sex`: X-het fraction and Y calls give "likely XX / XY / unclear" (informational only) | info |
| `env`: `VEP_*` variables other than `VEP_DIR_PLUGINS`, `~/.vep/vep.ini` | warn |
| `prepare_script`: `prepare-patient-vcf.sh -h` prints `Usage:` | error |
| `dbnsfp`: requested columns exist in the file header | error |
| `disk`: free space under `output.dir` at least 3x the input VCF | warn |

### Smoke test

`smoke` builds a VCF of heterozygous control records in `smoke.region` (default BRCA1, `17:43044295-43125483`): up to 25 records drawn from each tabix-backed plugin's own data file, plus ClinVar records, 10 per consequence stratum (missense, nonsense/frameshift, splice region, 5' UTR, intron, from `INFO/MC`). The ClinVar controls can have any clinical significance or review status, because a control only needs its consequence class. It then runs the normal pipeline on them. `smoke/smoke.tsv` has one row per enabled plugin, then one per enabled custom track (`custom:ClinVar`, `custom:ENCODE_cCRE`):

| Column | Meaning |
|---|---|
| `plugin` | Plugin name, or `custom:<track>`. |
| `status` | `OK`, `FAIL`, `unchecked` (no control defined, or no ClinVar control of that stratum in the region) or `info` (`custom:ENCODE_cCRE`). |
| `header_ok` | All declared output fields are in the CSQ header (`yes`/`no`). |
| `controls_n`, `controls_annotated` | Control records for this plugin, and how many received a value. |
| `eligible_n`, `annotated_n`, `coverage` | Over the ClinVar controls whose consequence the plugin can annotate: how many, how many got a value, and the ratio. Informational. |
| `message` | Reason for FAIL/unchecked. |

FAIL means one of:

- a field is missing from the header;
- the locus has no data for a data-file plugin, or a plugin annotated none of its controls;
- more than 2% of a data-file plugin's controls have a REF that is not the FASTA base (message: possible wrong assembly in plugin data). With 25 controls a single mismatch fails. `bcftools norm -c w` only warns about such records, so they would otherwise pass;
- a plugin whose controls come from ClinVar is enabled while the ClinVar custom track is disabled or has no file;
- `custom:ClinVar`: none of the ClinVar controls received a `ClinVar_CLNSIG` value, or there is no ClinVar record in the region.

`custom:ENCODE_cCRE` reports how many smoke records overlap a cCRE (informational). Plugins that are sparse by design (NMD, G2P) are checked for header presence, or against ClinVar controls of a matching consequence, instead of random coverage. The summary line counts and names the failed and the `unchecked` rows: an `unchecked` plugin passed without being tested.

## Data

### Reused, not downloaded

Paths relative to `/fast/AG_Bunina/Berk/projects/mva_challenge_26/`. These are read in place and never modified.

| Item | Path | Note |
|---|---|---|
| Reference FASTA | `references/Homo_sapiens.GRCh38.dna.primary_assembly.fa` (+`.fai`; 194 contigs; Ensembl names `1..22,X,Y,MT`) | identical copy in `exomiser/data/` |
| CADD v1.7 GRCh38 | `exomiser/data/whole_genome_SNVs.tsv.gz` (82 GB), `exomiser/data/gnomad.genomes.r4.0.indel.tsv.gz` (1.2 GB), both `.tbi` | header `##CADD GRCh38-v1.7`, which the plugin requires |
| ClinVar | `references/clinvar.vcf.gz` + `.tbi`, fileDate 2026-08-22 | custom track, `type=exact` |
| ENCODE cCREs | `references/GRCh38-cCREs.bed` (plain text, `chr` names) | a derived, indexed copy goes under `data.root/custom` (`GRCh38-cCREs.class.bed.gz`, name column `class:accession`) |
| Normaliser | `exomiser/bin/prepare-patient-vcf.sh` | contig-length gate, `chr` renaming, aborts above 2% REF mismatch, `-m -any`, `--atomize`, drops `*` alleles; does not subset samples; `-q` is opt-in and never used here, because it would drop mosaic hets below 15% VAF |
| Parallel downloader | `exomiser/bin/pget.sh <url> <out> [nconn] [chunk_MB]` | resumable, checks size |
| SLURM conventions | `exomiser/bin/exomiser.sbatch` | `partition=normal`, `--export=ALL`, `TMPDIR` on `SLURM_TMPDIR` |

### Downloaded by `download`

Ensembl files are pinned to `release-116`. URLs and sizes were checked by HEAD request on 2026-10-11. All land under `data.root` (`/fast/AG_Bunina/Berk/projects/mva_challenge_26/data/vep`).

| Item | URL | Bytes |
|---|---|---|
| VEP cache | `https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz` (CHECKSUMS `56036 26996736`) | 27,644,657,162 |
| SpliceAI, Ensembl, MANE v1.4, masked, SNV (only when `SpliceAI` is enabled; off by default, SpliceAI2 is the default) | `https://ftp.ensembl.org/pub/release-116/variation/SpliceAI/spliceai_scores.masked.snv.ensembl_mane_v1.4.grch38.vcf.gz` (+`.tbi`; CHECKSUMS `21141 27971711`) | 28,643,031,420 |
| popEVE | `https://data.evemodel.org/popeve/v1.1/downloads/grch38_popEVE_ukbb_20250715.vcf.gz` | 1,476,124,568 |
| SpliceVault | `https://ftp.ensembl.org/pub/release-116/variation/SpliceVault/SpliceVault_data_GRCh38.tsv.gz` (+`.tbi`) | 883,907,365 |
| REVEL v1.3 | `https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip` | 667,102,707 |
| AlphaMissense | `https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz` | 642,961,469 |
| MaveDB | `https://ftp.ensembl.org/pub/release-116/variation/MaveDB/MaveDB_variants.tsv.gz` (+`.tbi`) | 560,599,773 |
| LOEUF | `https://ftp.ensembl.org/pub/release-116/variation/LOEUF/loeuf_dataset_grch38.tsv.gz` (+`.tbi`) | 354,938,316 |
| PhenotypeOrthologous | `https://ftp.ensembl.org/pub/release-116/variation/PhenotypeOrthologous/PhenotypesOrthologous_homo_sapiens_112_GRCh38.gff3.gz` (+`.tbi`) | 3,424,946 |
| G2P (all panels; plain CSV, which the plugin requires) | `https://www.ebi.ac.uk/gene2phenotype/api/panel/all/download` (GET works; HEAD returns 405) | about 2.6 MB |
| UTRAnnotator | `https://raw.githubusercontent.com/Ensembl/UTRannotator/master/uORF_5UTR_GRCh38_PUBLIC.txt` | 2,376,724 |
| MechPredict | `https://osf.io/download/h45ns/` (`pgof_svm_poly_2023-07-25.tsv`), `https://osf.io/download/xfy38/` (`pdn_svm_poly_2023-07-25.tsv`), `https://osf.io/download/dj4qg/` (`plof_svm_poly_2023-07-28.tsv`) | about 0.7 MB each |
| DosageSensitivity | `https://zenodo.org/record/6347673/files/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz` | 386,018 |

- **Total:** about 61 GB with the Ensembl SpliceAI file, about 32 GB without it (the default; SpliceAI2 is a manual download). Peak disk use is about 90 GB while the cache tarball is extracted; `download` needs at least 2.2x the tarball size free, and deletes the tarball after verification and extraction.
- Downloads are atomic (`.part`, then rename), resumable and size-checked, verified against Ensembl `CHECKSUMS` where available, and recorded in `{data.root}/manifest.tsv` (`name url bytes md5 date`). The checksum is verified on the `.part` before the rename, so a file under its final name is always verified. Re-running `download` skips items that are complete. A download of unknown size (G2P API, OSF) restarts from zero instead of resuming a leftover `.part`.
- An unindexed custom VCF (ClinVar) is symlinked into `data.custom_data_dir` and indexed there; if `file` changes, the old link and its index are replaced.
- HTTPS calls set `SSL_CERT_FILE` and `CURL_CA_BUNDLE` to `/etc/ssl/certs/ca-certificates.crt` (the login environment points at a missing bundle and curl fails with error 77).
- `download` also derives the REVEL, MechPredict and cCRE files and the empty SpliceAI indel VCF. The cache `info.txt` must contain `var_type tabix`.
- Not downloaded: dbNSFP (needs registration, see below), SpliceAI_Illumina (manual BaseSpace download) and SpliceAI2 (manual Hugging Face download). `download --dry-run` lists them with instructions when enabled without a file.

### SpliceAI and the indel file

The SpliceAI plugin calls `add_file($indel)` even when no indel file is given, so `snv=` alone dies with "No file specified". Ensembl publishes SNV scores only, so the pipeline passes an empty, bgzipped, tabix-indexed, header-only VCF as `indel=`. `download` generates it at `{data.plugin_data_dir}/SpliceAI/empty_indel.vcf.gz`; `check` (`spliceai_form`) confirms the argument form.

## Output

Under `output.dir` (`/fast/AG_Bunina/Berk/projects/mva_challenge_26/results/vep`), with `prefix` = `{sample}`:

| Path | Content |
|---|---|
| `{prefix}.norm.vcf.gz` | Normalised VCF from `prepare-patient-vcf.sh` (only if `input.normalize`), with `{prefix}.norm.vcf.gz.fingerprint.json` (see the reuse rule above). |
| `{prefix}.input.vcf.gz` | Normalised VCF subset to `sample`, indexed. VEP never sees multi-sample input. |
| `work/contigs.txt` | One chunk name per line; SLURM array task N runs line N. Chunks are whole chromosomes; all non-primary contigs go into one `other` chunk. |
| `work/<chunk>/` | `vep.ini`, `run.sh`, `regions.bed` (the chunk's whole contigs), `in.vcf.gz`, `out.vcf.gz`, `skipped.txt`, `warnings.txt`, `stats.html`, `reconcile.tsv`, `run.log` (local runs). |
| `work/config.resolved.yaml`, `work/array.sbatch`, `work/merge.sh`, `work/merge.sbatch` | Resolved config and job scripts. |
| `{prefix}.vep.vcf.gz` (+`.tbi`) | Final merged VEP VCF with the CSQ field. |
| `fields.tsv` | CSQ field to evidence class (below). |
| `reports/` | G2P txt/html reports, one per chunk. |
| `logs/` | SLURM logs (`slurm.log_dir`): `vep_<arrayjob>_<task>.out`, `vep_merge_<job>.out`. |
| `smoke/` | Smoke VCF, its own `work/` and `reports/`, and `smoke.tsv`. |

Safety checks before and during a run:

- VEP runs with `safe 1`, so a plugin that fails to load aborts the chunk. `warnings.txt` is scanned for `Failed to|went wrong|did not return`, and the chunk script exits 3 on a match.
- Each chunk selects its whole contigs with `bcftools view -R work/<chunk>/regions.bed` (lengths from the VCF header, else the FASTA `.fai`), so contig names containing `:` or `*` (HLA alts) are never parsed as region strings. Each run first deletes the chunk's previous `in.vcf.gz`, `out.vcf.gz`, their indexes, `reconcile.tsv`, `warnings.txt` and `skipped.txt`, and uses its own `mktemp -d` TMPDIR, removed on exit.
- Every chunk is reconciled per variant key `(CHROM,POS,REF,ALT)`: each input record must appear in the output, be listed in VEP's skipped-variants file, or (G2P forces `individual` mode) have no alt genotype for the sample. The result is `reconcile.tsv` (`in out skipped no_alt_gt unexplained`). Any `unexplained` above 0 fails the chunk.
- Before `bcftools concat`, the `in` counts of all chunks must add up to the records of `{prefix}.input.vcf.gz` (SV job: `{prefix}.sv.input.vcf.gz`), so a record lost between the input and the chunks stops the merge. The `##INFO=<ID=CSQ` header line must be byte-identical across chunks and contain every expected plugin field. The merged record count must equal the chunk total.
- Generated scripts set `PERL_HASH_SEED=0`, `PERL_PERTURB_KEYS=0`, `LC_ALL=C` and unset every `VEP_*` variable except `VEP_DIR_PLUGINS`.

To redo one failed chunk, run what the array task runs, then the merge, with the resolved config the run was rendered with (not `config.yaml`, which may have changed since):

```bash
micromamba run -n vep116 bash <output.dir>/work/<chunk>/run.sh
bin/vep-pipeline reconcile --chunk <chunk> -c <output.dir>/work/config.resolved.yaml
bin/vep-pipeline run --merge-only -c <output.dir>/work/config.resolved.yaml
```

`run.sh` deletes the chunk's previous outputs and `reconcile.tsv` before it starts, so a rerun that fails cannot leave an older chunk for `--merge-only` to merge. `--merge-only` refuses to merge until every chunk has a clean `reconcile.tsv` and the chunk inputs add up to the input VCF.

### Optional SV run

For structural variants (for example calls from a later SV caller):

```yaml
sv:
  enabled: true
  vcf: /path/to/sv_calls.vcf.gz
```

- `render` / `run` subset the SV VCF to `sample` and add one extra job named `sv` (one more array task), with `max_sv_size=100000000`. Output: `{prefix}.sv.vep.vcf.gz`.
- The SV run uses the same ini as the SNV run, with full `--everything` including regulatory features, and **no plugins and no custom tracks**. DosageSensitivity crashes on SVs when regulatory is on (the SV overlap object has no transcript); the user chose regulatory annotation over dosage scores.
- SV reconciliation counts records only (`in = out + skipped`).
- Off by default; `enabled: true` without `vcf` is a config error.

## Enabling other plugins

`bin/vep-pipeline list-plugins` shows every registry entry. Categories: `no_data`, `public`, `existing`, `manual`, `user_file` can be enabled; `filter`, `grch37_only`, `db_online`, `broken` cannot.

```yaml
plugins:
  MaxEntScan: {enabled: true, args: "MaxEntScan,/path/to/maxentscan/fordownload"}   # registry has no recipe
  Blosum62: true                                                                    # no_data plugin, curated recipe
  SomeGRCh37Plugin: {enabled: true, args: "Foo,file=/x", allow_unsupported: true}   # non-runnable category
```

- **`args:`** is the escape hatch. It replaces the whole `--plugin` argument string, and you supply and index every file yourself (nothing is downloaded for it). `{chunk}`, `{workdir}` and `{report_dir}` are filled per chunk.
- **`allow_unsupported: true`** is needed in addition to `args:` for a plugin in a non-runnable category (for example GRCh37-only or database-backed ones). The registry says why each is blocked (`list-plugins` last column). Enabling one is at your own risk: it will probably fail, and `check` runs VEP on it with `safe 1`.
- Plugins with a `manual` category (dbNSFP, SpliceAI_Illumina, PrimateAI, ...) need a file path in the config; the error quotes the download instructions. dbNSFP needs registration at https://www.dbnsfp.org/download (institutional email, CC BY-NC-ND), a GRCh38 VEP-ready file with `.tbi` whose filename contains the version, and `columns:` (default `BayesDel_noAF_score`, `VEST4_score`, `ESM1b_score`). BayesDel and VEST4 are only available on GRCh38 through dbNSFP.
- SpliceAI_Illumina (manual BaseSpace download of the masked SNV and indel files) conflicts with SpliceAI; enable only one.
- SpliceAI2 (`SpliceAI2.pm`, added locally to VEP_plugins) reads Illumina's precomputed SpliceAI2 scores (GENCODE v48 protein-coding genes; non-commercial terms on Hugging Face). `snv` and `indel` are the directories of per-chromosome files (`chr*.tsv.gz`, a `glob` entry in the registry); every file needs a `.tbi` (`tabix -s 1 -b 2 -e 2 -S 1`). Scores are matched to each transcript's strand. SpliceAI2 is the default splice predictor (SpliceAI is off by default); it can run next to SpliceAI: the fields are `SpliceAI2_spliceai2_summary_score` (author cutoffs 0.1 / 0.25 / 0.5 for high recall / balanced / high precision) and the delta scores and distances of the strongest donor/acceptor gain/loss.
- Plugins without a registry recipe have no declared `output_fields`, so `fields.tsv` and the header checks do not cover their columns. Classify them yourself downstream.
- Plugin and custom-track paths can be overridden in the config (`snv`, `indels`, `file`, ...), which makes the file user-owned: it is verified, never written.

## `fields.tsv` and evidence-class ablation

`fields.tsv` (`field`, `evidence_class`, `source`) maps every CSQ field of this run to the kind of evidence it carries, so that a ranking can be rerun with a class removed. Entries ending in `*` are prefixes (`ClinVar_*`, `gnomADe_*`, `gnomADg_*`).

| Evidence class | Fields (source) |
|---|---|
| `consequence` | Consequence, IMPACT, SYMBOL, Gene, Feature, BIOTYPE, HGVSc, HGVSp, HGVSg, CANONICAL, MANE_SELECT, PICK (vep); SpliceRegion, NMD, Existing_uORFs, 5UTR_annotation (plugins) |
| `known_variant` | Existing_variation, CLIN_SIG, SOMATIC, PHENO, PUBMED, `ClinVar_*` |
| `population` | `gnomADe_*`, `gnomADg_*`, AF, MAX_AF, MAX_AF_POPS |
| `clinically_trained_predictor` | SIFT, PolyPhen, REVEL (and dbNSFP BayesDel, VEST4 if enabled) |
| `unsupervised_predictor` | am_pathogenicity, am_class, popEVE_SCORE, popEVE_gene, popEVE_gap_frequency, SpliceAI2_* (summary score, delta scores and distances), SpliceAI_pred_DS_AG/AL/DG/DL if SpliceAI is enabled, CADD_PHRED, CADD_RAW (and dbNSFP ESM1b if enabled) |
| `experimental_rnaseq` | SpliceVault_* |
| `experimental_assay` | MaveDB_score, MaveDB_urn |
| `gene_prior_dominant` | LOEUF, pHaplo, pTriplo, MechPredict_pDN/pGOF/pLOF/prediction |
| `gene_phenotype` | G2P_flag, G2P_complete, G2P_gene_req, PhenotypeOrthologous_Mouse_phenotype, PhenotypeOrthologous_Rat_phenotype, GENE_PHENO |
| `regulatory_overlap` | ENCODE_cCRE |

Ablation example: rank once with everything, once without `known_variant` and `clinically_trained_predictor`. A candidate that only ranks high with those classes is relying on curated knowledge or on training data that already contains it. Keep one extra run without `population` if you want to see what the ranking looks like without frequency filtering.

## Joining with Exomiser

The input is normalised by the same script as the Exomiser input (left-aligned, split, atomised), and VEP keeps the input coordinates, so both outputs carry identical keys. Join the VEP VCF (`CHROM`, `POS`, `REF`, `ALT`) to the Exomiser variants table on:

```
(CONTIG, START, REF, ALT)   <->   (CHROM, POS, REF, ALT)
```

Contig names follow the reference FASTA (`1..22,X,Y,MT`) on both sides. Each variant has one CSQ entry per transcript: select `PICK=1` for one row per variant, or keep all and group by `Gene`. Use the same `sample` as the Exomiser phenopacket `subject.id`.

## Prioritising the top 10

`prioritise` ranks the findings in the final VEP VCF and writes a CSV of the top N rows (default 10). It combines the VEP annotation, the G2P allelic requirement, the ACMG SF gene list, optional read-backed phase from the proband's BAM (`phase`) and, as the only source of phenotype evidence, Exomiser.

### Workflow

1. `bin/vep-pipeline run ...` produces `{prefix}.vep.vcf.gz` (VEP).
2. `bin/vep-pipeline exomiser --submit` writes the phenopacket (`exomiser.hpo` is required) and runs Exomiser as a SLURM job: a phenotype-only pass, then an exome and a genome pass. Results go under `{output.dir}/exomiser`.
3. Only with `phasing.enabled` (Model 3): when both are done, `bin/vep-pipeline phase --submit` runs WhatsHap on the candidate genes and writes the read-back sidecar `{phasing.results_dir}/{sample}.readback.tsv` (see "Read-backed phasing"). Wait until this job has **finished**, not only started. The job deletes the previous sidecar as its first step: a `prioritise` job that starts while the phase job is still pending reads the previous run's sidecar (if there is one), and one that starts while it runs finds no sidecar and stops.
4. Only with `regnet.enabled` (Model 4): `bin/vep-pipeline regnet-build` once, inside an allocation, before `prioritise` (see "Regulatory network (`regnet`)"). It does not read the sample's data, and the result is reused by every later run.
5. Then: `bin/vep-pipeline prioritise --slurm --submit --sensitivity`. The VCF scan runs as a SLURM job, never on the login node. Without `--submit` the command writes `work/prioritise.sbatch` (with a frozen copy of the resolved config, `work/prioritise.config.yaml`; both overwrite an earlier run's files) and prints the `sbatch` command. With `phasing.enabled` and no sidecar, `prioritise` stops with "run `vep-pipeline phase --submit` first". `--slurm` does not check the sidecar when it writes the job, so in that case the job itself fails.

There is no SLURM client in the code-server shell, so `--submit` (here, for `exomiser` and for `phase`) is run from a cluster login shell. The code only renders the scripts. `prioritise` without `--slurm` runs in the calling shell, which is fine for a small VCF or a test, not for a WGS VCF.

Flags: `-o/--output` (CSV path), `--top-n N`, `--sensitivity`, `--slurm`, `--submit` (only with `--slurm`).

### Outputs

Default path `{output.dir}/{prefix}.top{N}.csv`; the other files replace the `.csv` suffix. The read-back sidecar is the exception: `phase` writes it to `phasing.results_dir`.

| File | Content |
|---|---|
| `*.csv` | Header `proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type,notes`. One row per finding: primary rows first (by epcr), then at most `secondary.max_rows` secondary rows. A gene can fill up to `alternative_pairs.max_per_gene` rows (its best finding plus alternative pairs). Slot 2 is empty unless the finding is a pair (two alleles in one gene). `proband_id` = `sample`. Coordinates are the original VCF record (from `EXO_ORIG`) when normalisation changed it, else the normalised key; contigs get `chrom_prefix` (`MT` becomes `chrM`). `notes` holds tier, gene, hypothesis and fit, phenotype strength (and `panel=<name>` for a panel gene) and ranks, and for each allele the HGVS, zygosity, evidence level and reasons, popmax, low-VAF flag, and for a pair the phase with its sources and the distance. |
| `*.candidates.tsv` | The kept rows per gene, chosen separately for primary and secondary findings: the best finding and any alternative pairs, after budget sharing but before selection and the cap (`gene, finding_type, tier, epcr_before_cap, phenotype, hypothesis, fit, level, phase, variant_1, variant_2, notes`). Not every scored finding: it is for auditing why a gene did or did not make the list. |
| `*.compound_het.tsv` | Every pair of het candidate alleles (after QC) in every recessive or unknown-MOI gene, cis and second-hit-only pairs included, in gene then position order. See below. |
| `*.excluded.tsv` | `variant, region`: every record with an alt genotype (PASS) inside an `exclude_regions` interval, sorted and de-duplicated. None of them formed a finding. |
| `*.sensitivity.tsv` | With `--sensitivity`: the selected rows under every pre-declared scenario (see "Sensitivity scenarios"). |
| `{phasing.results_dir}/{sample}.readback.tsv` | Written by `phase`, read by `prioritise` (see "Read-backed phasing"). |
| `*.regulatory.tsv`, `*.regulatory.summary.tsv` | Only with `regnet.enabled`: the ranked rare non-coding variants in gene-linked peaks, and a key/value summary of the run (see "Regulatory-network outputs"). |

`*.compound_het.tsv` columns: `gene, variant_1, variant_2, hgvs_1, hgvs_2, level_1, level_2, reasons_1, reasons_2, phase, phase_sources, distance, collapsed_into, fit, tier, epcr_before_cap, reportable, csv_row`. The phase comes from the same phase graph as the ranking.

- `reportable`:
  - `yes`: the pair was scored as a finding. `fit`, `tier` and `epcr_before_cap` are filled from the kept (budget-shared) row when the pair was kept, else from its best-scoring finding. `yes` does not mean the pair is in the CSV; `csv_row` says that.
  - `no (cis)`, `no (two second-hit-only alleles)`, `no (same variant)`;
  - `no (collapsed)`: one allele is a cis mate folded into a representative pair, named in `collapsed_into` (`key|key`);
  - `no (not paired)`: no pair finding covers the two alleles.
- `csv_row`: the 1-based row in the CSV, blank when the pair is not in it.

### Evidence levels

One predictor per evidence type (Pejaver 2022: avoid scanning multiple tools).

- **Missense:** the tool in `missense_tool` (default `alphamissense`) decides. The other tool, REVEL, is used only when the chosen one has no score or its evidence class is excluded. A score below the supporting cut-off is final: the other tool is not consulted. `missense_tool: revel` restores Model 2's rule (REVEL first, AlphaMissense only when REVEL is blank or excluded).
- **Splicing:** SpliceAI / SpliceAI2, at supporting weight only.
- **CADD:** used only when nothing else scored the variant, and stands alone only for protein-altering consequences (see below).

Thresholds, checked against the source tables on 2026-10-17:

| Tool | PP3 supporting | PP3 moderate | PP3 strong | Source |
|---|---|---|---|---|
| REVEL | [0.644, 0.773) | [0.773, 0.932) | >= 0.932 | Pejaver 2022 Table 2, doi:10.1016/j.ajhg.2022.10.013 |
| AlphaMissense | [0.792, 0.906) | [0.906, 0.990) | >= 0.990 | Bergquist 2025 Table 1, doi:10.1016/j.gim.2025.101402 |
| CADD PHRED | [25.3, 28.1) | >= 28.1 (literature only; not used) | - | Pejaver 2022 Table 2 (missense calibration) |
| SpliceAI delta | >= 0.2 (applied at supporting weight) | - | - | Walker 2023, doi:10.1016/j.ajhg.2023.06.002 |

The table is the literature reference. The pipeline applies CADD at supporting weight only (PHRED >= 25.3, never moderate), and only when REVEL, AlphaMissense and SpliceAI / SpliceAI2 are all blank for the variant.

Non-coding signals are second-hit-only. A second-hit-only signal never ranks a gene alone, but can complete a recessive pair (see "epcr is heuristic"). They are:

- **Non-coding CADD.** CADD >= 25.3 is a standalone supporting hit only when the representative consequence is protein-altering (missense, protein_altering, inframe insertion or deletion, stop_lost, start_lost). For every other consequence (intronic, UTR, up/downstream, synonymous, regulatory, non-coding transcript) it is second-hit-only. Basis: missense tools "cannot be applied to variants in non-coding regions", and non-coding variants are interpreted only in phenotype-matched genes, as an in-trans second hit (Ellingford 2022 Genome Med, doi:10.1186/s13073-022-01073-3). The missense cut-off also does not transfer: the promoter calibration gives CADD PP3 at >= 10, far from 25.3 (Villani 2024 AJHG, doi:10.1016/j.ajhg.2024.05.002).
- **Genomiser ReMM >= 0.963** on a non-protein-altering variant (Smedley 2016 AJHG, doi:10.1016/j.ajhg.2016.07.005). The cut-off is the one Xi 2025 used as a filter on 985 WGS (Brief Bioinform, doi:10.1093/bib/bbaf475); they found ReMM best as a filter, not as evidence, hence second-hit-only. ReMM is read from `exomiser/genome/<sample>.variants.tsv` (`ALL_PATH`, max over modes of inheritance). It exists only for variants the Genomiser pass kept (its own filters apply), and is absent when Exomiser is off. ReMM is trained on curated pathogenic regulatory variants, so it is removed with the `clinically_trained_predictor` class.
- A 5' UTR uORF annotation, or an ENCODE cCRE overlap with CADD >= 20.
- In a gene-panel gene only, a "panel candidate" from the deep second-hit search (see "Gene panel").

Null variants follow the PVS1 decision tree in outline (Abou Tayoun 2018, doi:10.1002/humu.23626):

- start_lost or NMD-escaping: moderate.
- G2P molecular mechanism other than loss of function: supporting.
- In a dominant or unknown-MOI gene without a G2P loss-of-function mechanism, `strong` needs LOEUF < 0.6 or pHaplo >= 0.86 (Collins 2022, doi:10.1016/j.cell.2022.06.036).

Two rules decide what is not evidence at all:

- **Canonical splice calls moved by the HGVS shift (`splice_shifted_out`).** VEP calls the consequence at the left-aligned VCF position. For an indel in a tandem repeat that spans a splice site, the 3'-shifted (equivalent) allele can leave the +-1,2 dinucleotide intact. When the HGVSc (`c.`) puts the indel wholly inside one intron, more than 2 bp from the exon, the terms `splice_donor_variant` and `splice_acceptor_variant` are dropped before the evidence rules run, so the variant is not a PVS1 null. Example: a repeat-unit deletion that VEP calls `splice_donor_variant` at an exon's donor site, but whose 3'-shifted HGVSc (for example c.100+77_100+121del) lies wholly inside the intron. The rule is always on (no setting). Only `c.` HGVS is checked; a non-coding transcript (`n.`) keeps VEP's call.
- **Excluded regions (`exclude_regions`).** Records inside a listed interval are skipped before any evidence is computed: they never form a primary or secondary finding, in any scenario, and are listed in `*.excluded.tsv`. Default: the GRC-defined MHC, chr6:28,510,120-33,480,577 (GRCh38), where mapping and reference-bias artefacts formed findings in Model 2. Contig names match with or without `chr`. `exclude_regions: []` is Model 2 behaviour.

Population frequency is the popmax over gnomAD AFR, AMR, EAS, NFE and SAS (exomes and genomes), the same set Exomiser uses. The dominant (1e-4) or recessive (1e-2) cut-off is applied once the inheritance mode is known. Low-VAF hets are kept, flagged and lowered one level; nothing is dropped for VAF.

### Tier table

The default table (`prioritise.tiers`; first match wins; `*` matches anything). A rule is `[variant level, MOI fit, phenotype]`.

| tier | epcr | rules |
|---|---|---|
| T1 | 0.90 | `[strong, full, strong]` |
| T2 | 0.60 | `[strong, full, moderate]`, `[moderate, full, strong]`, `[strong, unknown_moi, strong]` |
| T3 | 0.30 | `[moderate, full, moderate]`, `[supporting, full, strong]`, `[moderate, unknown_moi, strong]` |
| T4 | 0.10 | `[*, full, moderate]`, `[*, unknown_moi, moderate]`, `[strong, single_hit, strong]`, `[strong, full, unknown]`, `[strong, full, none]`, `[*, panel_partner, strong]` |
| T5 | 0.02 | `[*, *, *]` |

`single_hit` is one het in a recessive gene with no second allele found. It stays below an unphased pair with a weak second allele, because 0.30 x 0.5 = 0.15 > 0.10. `unknown_moi` is a gene without a usable G2P allelic requirement. `panel_partner` is a pair in which one allele is only a gene-panel candidate (see "Gene panel"); it reaches T4 only at a strong phenotype and falls to T5 otherwise.

### epcr is heuristic

`epcr` is not a calibrated probability.

- A primary row gets the tier's fixed value.
- A pair not proven trans is multiplied by `p_trans_unphased` (0.5). Basis: in gnomAD (general population), 86% of co-occurring rare pLoF pairs are in cis, whereas ~96% of phaseable pairs in diagnosed patients were in trans, a selected cohort (Guo 2024 Nat Genet, doi:10.1038/s41588-023-01608-3). No paper derives P(trans) for an affected singleton, so 0.5 is a neutral default; the sensitivity scenarios use 0.15 and 0.9 to bracket it.
- Second-hit-only signals (5' UTR uORF, cCRE with CADD >= 20, non-coding CADD >= 25.3, ReMM >= 0.963, panel candidates) never rank alone, but they do enter epcr when they complete a recessive pair with another qualifying allele: the pair is scored at its weaker level and discounted by `p_trans_unphased` unless proven trans. The partner does not have to be a null.
- The alternative pairs of one gene share one budget (see "Alternative pairs and the shared budget").
- The primary rows are then scaled together so their sum is at most `primary_sum_cap` (0.95), because only 4.9% of diagnosed patients had two or more molecular diagnoses (Posey 2017, doi:10.1056/NEJMoa1516767). The unscaled value is in `*.candidates.tsv`.
- A secondary row gets the confidence that the variant is a reportable P/LP: ClinVar >= 2 stars 0.30, 1 star 0.20, null variant in a loss-of-function gene without ClinVar 0.10. Secondary findings are ACMG SF v3.3 (84 genes; Lee 2025, doi:10.1016/j.gim.2025.101454) without a strong or moderate phenotype match (a gene with one is scored as primary), and their notes start with "Incidental finding unrelated to primary phenotype".
- Exomiser scores are ranks, not probabilities. They are used only to rank genes, never as a number in epcr.

### Phenotype

- Phenotype evidence comes only from Exomiser. The strength is the gene's rank in a **phenotype-only** pass: strong <= 20, moderate <= 200 (`prioritise.pheno_rank`; fixed before the first real run, do not tune on the result). The phenotype-only pass is used because the exome and genome passes drop genes without a passing variant, which would make the rank depend on the variants.
- A gene absent from that output is `unknown`; a gene present but ranked beyond the moderate cut-off is `weak`. With `exomiser.enabled: false` the phenotype is `none` and no row can rank above T4.
- A gene-panel gene can have its phenotype raised by `gene_panel.phenotype_floor` (see "Gene panel").
- TP53, WT1 and RB1 can be primary findings when the phenotype supports them; without a phenotype match they are secondary findings (they are in ACMG SF).

### Phase of pairs

Phase comes from a phase graph over the het alleles of each gene. Its sources are `prioritise.phase_sources` (default: all four):

| Source | Decides a pair when |
|---|---|
| `multiallelic` | the two alleles are ALTs of one original multiallelic record (same `EXO_ORIG` record, different ALT index): trans |
| `ps` | both GTs are phased (`\|`) with the same `FORMAT/PS`: cis if the ALTs are on the same haplotype, else trans |
| `pid` | both `FORMAT/PGT` are phased with the same `FORMAT/PID` (GATK HaplotypeCaller): as for `ps` |
| `whatshap` | both alleles are in the read-back sidecar with the same PS (see "Read-backed phasing"): as for `ps` |

- Every source that decides a pair adds an edge, and relations follow transitively. Example: A and B are split from one multiallelic record (trans), and B and C share a PS block on the same haplotype (cis), so A and C are trans with sources `multiallelic+ps`.
- An edge that contradicts what the graph already implies (for example cis by PS and trans by PID) marks the whole connected component as conflicted. Every pair in it is then `unknown` with source `conflict`.
- **Cis haplotype collapse.** Alleles proven cis form one haplotype. Only its strongest member (level, then not second-hit-only, then lower popmax, then position) is paired; its cis mates are folded into it (`collapsed_into` in `*.compound_het.tsv`; the CSV notes do not list them). A cis pair is never reported, so a gene whose only two alleles are in cis gives single findings only.
- The remaining pairs are hypotheses: never one variant twice, never two second-hit-only alleles. A gene's pairs are ordered best first by the weaker allele's level, the stronger allele's level, proven trans, the better panel class, then the lower popmax sum. This is Model 2's `best_pair` order plus the panel class.
- A pair that is not proven trans (`unknown`, including `conflict`) is multiplied by `p_trans_unphased`. The notes and `*.candidates.tsv` show `phase=<state> (<sources>)`.

Model 2 decided phase pair by pair (multiallelic, then PS, then PID; the first source that decided won), without transitivity, conflict detection or cis collapse.

### Read-backed phasing (`phase`)

`phase` writes one SLURM job, `{phasing.results_dir}/phase.sbatch`, with a frozen copy of the resolved config (`phase.config.yaml`). `--submit` submits it; without it the command prints the `sbatch` command. `--prepare` and `--finish` are steps of the job. WhatsHap (read-based phasing; Patterson 2015, doi:10.1089/cmb.2014.0157) and samtools run in the micromamba env `phasing.env`; the Python steps and bcftools run in `tools.env`.

| `phasing` key | Default | Meaning |
|---|---|---|
| `enabled` | false (Model 2) | true: `phase` can run, and `prioritise` requires the sidecar |
| `bam` | none | the proband's BAM or CRAM; its `.bai` / `.crai` must exist |
| `reference` | none | the FASTA the BAM was aligned to; it needs a `.fai`, and a `.gzi` when bgzipped |
| `env` | none | micromamba env with `whatshap` and `samtools` |
| `results_dir` | `{output.dir}/phasing` | job script, work files and the sidecar |
| `pad` | 1000 | bp added on each side of a gene's span of candidate alleles |
| `indels` | true | pass WhatsHap's legacy `--indels` flag (see below) |
| `mem`, `time` | `16G`, `08:00:00` | job resources (4 CPUs) |

With `enabled`, the settings check that `bam`, `reference` and `env` are set, that both files exist, and that the BAM/CRAM index and the FASTA `.fai` exist; otherwise they raise a config error. The `.gzi` is not checked.

The job runs these steps in order:

1. It deletes the previous run's sidecar, `phased.vcf.gz` and its index (`--prepare` deletes them again). A job that fails later therefore leaves no sidecar for `prioritise` to read.
2. `samtools view -H` reads the BAM's contig names.
3. `phase --prepare` scans the VEP VCF with the prioritise settings, as `prioritise` does (with Exomiser enabled it needs the Genomiser output for ReMM). It writes `regions.bed`: for every gene with at least two het candidate alleles that pass QC, the span from the first to the last allele, plus `pad` on each side, with overlaps merged. It maps the VCF contigs to the BAM's (same name, `chr` added or stripped, `MT` to `chrM`; a contig that cannot be mapped is an error). It then extracts **all** het PASS calls of the proband in those regions from `{prefix}.input.vcf.gz`, because the intermediate het SNPs are what chain reads into a block, and renames their contigs to the BAM's (`subset.vcf.gz`). When no gene has two het alleles it writes an empty sidecar, and the job skips the next two steps.
4. `whatshap phase --reference <reference> --ignore-read-groups [--indels] -o phased.vcf.gz subset.vcf.gz <bam>`. The `.gz` output name makes WhatsHap write BGZF; `bcftools index -t` then indexes it.
5. `phase --finish` maps the contigs back, writes the sidecar atomically (a temporary file moved into place) and prints the counts: phased records, blocks, and candidate pairs whose two alleles share a block.

`indels`: WhatsHap 2 and later phases indels by default (only `--only-snvs` restricts it to SNVs) and has no `--indels` option, so passing it fails. `indels: true` passes the legacy flag (WhatsHap 1.x); `false` omits it. This is a tool-compatibility setting, decided before the run; Model 3 runs WhatsHap 2.8 with `indels: false`.

The sidecar `{phasing.results_dir}/{sample}.readback.tsv` is tab-separated with the header `key`, `gt`, `ps`. It has one row per phased record (GT with `|`, PS set). `key` is the normalised VEP key `chrom:pos:ref>alt` with the VCF's contig names. `prioritise` loads it as the `whatshap` phase source: it skips unphased rows, and a row without exactly three fields is an error ("re-run `vep-pipeline phase --submit`").

The regions follow the base settings only. The sensitivity scenarios reuse the same sidecar, and `no_read_backed_phase` drops the `whatshap` source.

### Gene panel (`gene_panel`)

Off by default (`name: ""`, `genes: []`, `phenotype_floor: null`, `deep_second_hit: false`; Model 2 behaviour). `phenotype_floor` and `deep_second_hit` need at least one gene. No gene name is in the code. Model 3 sets the MVA panel, BUB1B, CEP57 and TRIP13 (Hanks 2004, doi:10.1038/ng1449; Snape 2011, doi:10.1038/ng.822; Yost 2017, doi:10.1038/ng.3883), with `phenotype_floor: strong` and `deep_second_hit: true`. The panel was chosen after Model 2's score (see "Model history").

- **Phenotype floor.** A panel gene's phenotype is raised to at least `phenotype_floor` (`strong` or `moderate`); a stronger Exomiser match is kept. The notes add `panel=<name>`, plus `(phenotype floor)` when the floor raised it. The raised phenotype is used everywhere: tiers, secondary-finding eligibility and `alternative_pairs.min_phenotype`.
- **Deep second-hit search** (`deep_second_hit`). In a panel gene, a variant with no calibrated evidence and no other second-hit-only signal becomes a supporting, second-hit-only "panel candidate (<class>)" when it falls in one of these classes. They are checked in priority order, and the first match names the class:
  1. `splice`: a usable SpliceAI / SpliceAI2 score **> 0.1** (usable: present, and `unsupervised_predictor` not excluded); or, when there is no usable score, a splice-region term (`splice_region_variant`, `splice_donor_5th_base_variant`, `splice_donor_region_variant`, `splice_polypyrimidine_tract_variant`). A usable score <= 0.1 is in the SpliceAI BP4 range (Walker 2023, doi:10.1016/j.ajhg.2023.06.002) and vetoes the term, so the threshold is strictly above 0.1. Scores in (0.1, 0.2) are "not benign" without being PP3 evidence. A vetoed variant can still fall into a later class: an exonic splice-region missense with SpliceAI 0.02 is `protein_altering`, while an intronic splice-region variant with SpliceAI 0.05 is no candidate.
  2. `protein_altering`: missense, protein_altering, in-frame insertion or deletion, stop_lost, start_lost. There is no BP4 filter here: the AlphaMissense BP4 cut-offs are not applied, so a missense variant that AlphaMissense scores as benign is still a candidate.
  3. `5utr`, then `3utr`.

  ClinVar benign/likely benign with >= 1 star still excludes the variant, and the frequency gate is unchanged (popmax <= `af.max_recessive`). Upstream/downstream, plain intronic and synonymous variants without a splice signal are not candidates (no cis-regulatory search; see "Not covered"). Basis: second-hit search in phenotype-matched recessive genes (Ellingford 2022), and MVA second alleles that are often hypomorphic or non-coding (Matsuura 2006, doi:10.1002/ajmg.a.31069; Rio Frio 2010, doi:10.1056/NEJMoa1006565; Sieben 2020, doi:10.1172/JCI126863).
- **`panel_partner`.** A panel candidate never forms a single finding, not even as a homozygote, and two panel candidates never pair. A pair with one panel candidate gets the fit `panel_partner`: tier T4 (0.10) only at a strong phenotype (`[*, panel_partner, strong]`), T5 (0.02) otherwise, times `p_trans_unphased` unless proven trans. At equal epcr a `panel_partner` row sorts after the other rows (across genes too), and among one gene's pairs the better panel class ranks first. In a recessive gene at a strong phenotype (always the case for a panel gene under `phenotype_floor: strong`), every calibrated pair is at least T3 (0.30, or 0.15 unphased), so a `panel_partner` pair (at most 0.10) never outranks a calibrated pair of the same gene, with or without a shared budget. Outside that setting it can: without a floor, at a weak or unknown phenotype both are T5 and a proven-trans `panel_partner` pair (0.02) outranks an unphased calibrated pair (0.01); and in an unknown-MOI panel gene a supporting calibrated pair is T5 while the `panel_partner` pair is T4.

### Alternative pairs and the shared budget (`alternative_pairs`)

Default `{max_per_gene: 1, min_phenotype: strong}`: one finding per gene, Model 2 behaviour. Model 3 sets `max_per_gene: 3`.

- **Selection.** Each gene keeps its best finding, as before; a finding whose allele set was already kept for another (overlapping) gene is skipped. If the gene's phenotype (after the floor) is at least `min_phenotype`, it also keeps up to `max_per_gene - 1` further **pair** findings, in sort order, each with a new allele set. Alternatives are separate rows and can enter the CSV.
- **Pair choice.** Rows sort by epcr, then calibrated before `panel_partner`, the Exomiser combined score, the weaker allele's level, the gene's pair order (`all_pairs`, see "Phase of pairs"), popmax and position. Within a gene the pair choice therefore follows the `all_pairs` order as a tie-break after the level. One deviation from Model 2 (user ruling 2026-10-20): when a stronger unphased pair and a weaker proven-trans pair compete, the pair with the higher epcr is reported. Example: null + null with unknown phase (T1, 0.90 x 0.5 = 0.45) loses to null + moderate missense proven trans (T2, 0.60); Model 2's `best_pair` took the null + null pair.
- **Shared budget.** A gene's kept pairs exclude each other (only one can be the patient's genotype), so two or more of them share one budget. With B = the largest tier value among them (before the phase discount), each gets epcr_i = B x tier_i / sum(tier) x phase factor_i, where the phase factor is `p_trans_unphased` for a pair not proven trans, else 1. Single-hit rows are nested hypotheses, not exclusive ones, and keep their epcr. The rows are then re-sorted and the primary cap applies. Worked example, a BUB1B-like gene with one panel partner: null + moderate missense (T2, 0.60) and null + panel candidate (T4, 0.10) give 0.60 x 0.60 / 0.70 = 0.514 and 0.60 x 0.10 / 0.70 = 0.086, both x 0.5 unphased: **0.257 and 0.043**, still above T5 noise (0.02).
- **Secondary rows follow the same rule** (user ruling 2026-10-20). By definition they have a weak, unknown or none phenotype, so with `min_phenotype: strong` or `moderate` they never get alternatives. A `min_phenotype` of `weak`, `unknown` or `none` lets an ACMG SF gene at or above it keep alternative pairs, and their shared budget then splits the SF confidence (0.30 / 0.20 / 0.10) between them.

### Regulatory network (`regnet`)

Off by default (`enabled: false`; Model 2 and 3 behaviour): with it off, `prioritise` writes byte-identical outputs and `regnet-build` is not needed. Use it to look for a rare non-coding allele that VEP does not link to a gene (no `SYMBOL`, so it is dropped today) but that lies in an enhancer peak linked to a phenotype gene in a single-cell regulatory network (scDoRI), where the transcription factors (TFs) binding the peak are enriched for heritability of a pre-declared GWAS trait (S-LDSC). The result is a separate report (`*.regulatory.tsv`), and optionally a supporting, second-hit-only allele in the recessive-pair ranking (see "Regulatory-network second hit"). Standalone non-coding rows never enter the top-N CSV.

The idea follows sc-linker (Jagadeesh 2022: single-cell programs, enhancer-gene links and S-LDSC) and the convergence of common and Mendelian disease genes on phenotype-matched traits (Freund 2018). It is a prioritisation heuristic, not a calibrated pathogenicity score; see "Regulatory-network limitations".

**`regnet-build`** is proband-independent. It streams the scDoRI eGRN edge lists once (about 112 M lines, one file per cell type), keeps only edges on gene-linked peaks, and compiles everything into `{regnet.build_dir}` (default `{output.dir}/regnet`). On the atlas it was developed with (about 38 k linked peaks, 134 cell types) it takes about 3 minutes and 5.5 GB of memory; run it inside an allocation.

| Input | Read from | Used for |
|---|---|---|
| scDoRI `gene_peak_interactions/<GENE>_enhancers_peak_links_.tsv` | `regnet.scdori_dir` (the scDoRI postprocess directory) | peak to gene links: rows are peaks named `chrN:start-end`, columns are cell types; the gene comes from the file name (it may contain `_`). A peak listed in a gene's file is linked to it. |
| scDoRI `eGRN_act_peak_gene/<cell type>_eGRN_edge_list.tsv` | same | TF to peak edges (`source target weight`), kept only for linked peaks; the TFs of a peak are those with an edge in any cell type |
| scDoRI `GRN_act_TF_gene/<cell type>_GRN_edge_list.tsv`, `<cell type>_GRN.tsv` | same | TF to gene regulons (union over cell types) and the gene universe (the row names of the `_GRN.tsv` matrices) |
| S-LDSC `h2_cts/<model>/<trait>.cell_type_results.txt`, `full/passing.txt` | each directory of `regnet.sldsc_dirs` | per-TF heritability coefficient (tau) and one-sided p-value for each trait that passes QC; a `passing.txt` line with a comma is an rg pair and is ignored |
| trait labels | `regnet.trait_labels` | trait id to label (public ontology labels), for `traits.tsv` and the summary |

S-LDSC directories are read in the listed order and a trait is taken from the **first** directory that has a result file for it. Within each trait the p-values are Benjamini-Hochberg (BH) adjusted across the TFs. `regnet.sldsc_model` selects the TF annotation (`M1` reads the `TF_<name>` rows, `M1bin` the `TFbin_<name>` rows).

| Output in `{regnet.build_dir}` | Columns / content |
|---|---|
| `linked_peaks.tsv` | `peak chrom start end gene link_max` (`link_max` = maximum over cell types; the sign is not used) |
| `peak_tfs.tsv` | `peak tf n_celltypes max_weight celltypes` |
| `regulons.tsv` | `tf gene n_celltypes` |
| `gene_universe.txt` | one scDoRI gene symbol per line |
| `sldsc.tsv` | `trait source tf coef se p fdr` (`source` = the S-LDSC directory name) |
| `traits.tsv` | `trait source label n_tf` |
| `build.json` | input paths, file counts, row counts, regulon size summary, UTC time |

`regnet` block (validated by `prioritise`, also at `--slurm` render time, and by `regnet-build`; an unknown key is an error):

| Key | Default | Meaning |
|---|---|---|
| `enabled` | false | true: `prioritise` loads the resource, writes the regulatory TSVs and (with `second_hit`) adds the second hits. Needs `traits` and a built `build_dir` that contains them. |
| `scdori_dir` | none | scDoRI postprocess directory; `regnet-build` only |
| `sldsc_dirs` | `[]` | S-LDSC result directories, in precedence order (first directory wins per trait); `regnet-build` only |
| `sldsc_model` | `M1` | `M1` or `M1bin` (the TF annotation read from `h2_cts`) |
| `build_dir` | `{output.dir}/regnet` | where `regnet-build` writes and `prioritise` reads |
| `trait_labels` | `vep_pipeline/data/sldsc_trait_labels.tsv` | `trait<TAB>label` TSV joined into `traits.tsv`; null or absent gives the packaged file (labels cannot be switched off), an empty string is an error |
| `traits` | `[]` | the **pre-declared** S-LDSC trait ids used for ranking. They are matched to the HPO terms by hand, from the label table, before any run and before looking at any variant; none is added afterwards. Required when `enabled`, and each must be in `traits.tsv`. No match for a term is a valid outcome. |
| `fdr` | 0.10 | cut-off for an S-LDSC TF flag and a regulon flag; a number in (0, 1] |
| `second_hit` | true | true: eligible rows become second-hit-only alleles in the pair ranking; false: only the report is written |
| `null_draws` | 1000 | number of random trait sets in the null (positive integer) |
| `seed` | 2026 | seed of the null (non-negative integer); the same seed gives the same draws |

#### Using your own single-cell GRN and S-LDSC results

`regnet-build` finds its inputs by directory and file-name pattern only, so the output of your own scDoRI run (or of another single-cell GRN tool, reshaped to this layout) and your own S-LDSC runs can be plugged in through `regnet.scdori_dir` and `regnet.sldsc_dirs`:

```
<scdori_dir>/
  gene_peak_interactions/<GENE>_enhancers_peak_links_.tsv  # header; then per peak: chrN:start-end, one link score per cell type
  eGRN_act_peak_gene/<cell type>_eGRN_edge_list.tsv        # header; then exactly 3 columns: TF, peak, weight
  GRN_act_TF_gene/<cell type>_GRN_edge_list.tsv            # header; then exactly 3 columns: TF, gene, weight
  GRN_act_TF_gene/<cell type>_GRN.tsv                      # header; first column = gene symbol (other columns ignored)
<sldsc_dir>/                                               # one or more, listed in regnet.sldsc_dirs
  full/passing.txt                                         # trait ids that passed your QC, one per line
  h2_cts/M1/<trait>.cell_type_results.txt                  # ldsc.py --h2-cts output (h2_cts/M1bin/ for sldsc_model: M1bin)
```

- **Files are tab-separated** and each has one header row, which is skipped. Columns are read by position, except the S-LDSC files (by name).
- **Peaks** are named `chrN:start-end` (with `chr`; N is 1-22, X, Y or M), and the eGRN `peak` column must use exactly the names of the links files. A malformed peak name in a links file stops the build. Blank or `nan` link scores are ignored, and the score itself is not a filter: every peak listed in `<GENE>_..._.tsv` counts as linked to `<GENE>`.
- **Genes** are taken from the links file name with the suffix `_enhancers_peak_links_.tsv` stripped (a symbol may contain `_`). Use the gene symbols Exomiser reports; a symbol Exomiser does not rank becomes `unknown` (R5).
- **Cell types** are the file-name prefixes of the edge lists and need not match across subdirectories; a TF counts on a peak if it has an edge in any cell type.
- **S-LDSC, one annotation per TF.** Run `ldsc.py --h2-cts` with an `.ldcts` file whose first column (the `Name` column of the results) is `TF_<TF>` for `M1` (or `TFbin_<TF>` for `M1bin`), where `<TF>` is spelled exactly as in the eGRN TF column; rows with another prefix (cell-type or baseline annotations) are ignored. The results need the standard columns `Name`, `Coefficient`, `Coefficient_std_error` and `Coefficient_P_value`. The trait id is the file-name stem. Only traits listed in `full/passing.txt` are read (a line containing a comma is an rg pair and is skipped; a listed trait without a result file is skipped).
- **Trait labels.** `regnet.trait_labels` is a TSV with the header `trait<TAB>label`. The packaged `vep_pipeline/data/sldsc_trait_labels.tsv` holds public ontology labels (EFO, MONDO, HP ids written with `_`); for other trait ids pass your own file. A trait without a label still works (blank label); labels are only for reading `traits.tsv` and the summary.

Then run `regnet-build`, read `{build_dir}/traits.tsv` to choose the traits that match the phenotype (by label, before any run and before looking at variants), list them in `regnet.traits` and set `regnet.enabled: true`.

### Regulatory-network scoring

During the VCF scan (`candidates.collect`, no second pass) a hook captures a record when it is:

- **not mitochondrial and not in an excluded region** (`exclude_regions`), and its alt allele is not `*`;
- **rare**: popmax over `af.populations` (gnomAD exomes and genomes) at most `af.max_recessive`. The gate is off when the `population` class is excluded;
- **non-coding**: no CSQ entry of the record has a coding term (null, protein-altering, in-frame, splice-region, synonymous, stop/start-retained, incomplete terminal codon or coding-sequence variant). A record with no CSQ entry is not captured. VEP's `SYMBOL` is not needed;
- **in a gene-linked peak**: the record's REF span overlaps a peak of `linked_peaks.tsv` (inclusive bounds, `start <= pos <= end`; contig names match with or without `chr`).

With `qc.pass_only` true (the default) non-PASS records are never read, so they are never captured. Calls that fail the GQ or DP cut-off are dropped at scoring (`hits` vs `hits_qc_pass`). Each remaining call gives one row per linked gene of each overlapping peak, de-duplicated per (variant, gene) to its best peak (tier, then flagged-TF counts, then peak name).

For a row, with G the linked gene:

- **pheno** is G's phenotype strength as in the ranking (`pheno_rank` cut-offs on the Exomiser phenotype-only rank, with the gene-panel floor). A gene Exomiser did not rank is `unknown`; with Exomiser off the phenotype is `none` (every row is R5).
- **V** (variant signal) is CADD PHRED >= 20.0 (the maximum over the record's CSQ entries; the same constant as the cCRE signal) or ReMM >= 0.963 (the second-hit cut-off above). CADD is ignored when `unsupervised_predictor` is excluded, and ReMM when `clinically_trained_predictor` is excluded.
- **S** (S-LDSC support): at least one TF on the peak is S-LDSC-flagged. A TF is flagged when, for **any** declared trait, its tau is positive and its BH q-value within that trait (across TFs) is below `fdr`.

| pheno of G | V and S | V or S (one) | neither |
|---|---|---|---|
| strong | R1 | R2 | R3 |
| moderate | R2 | R3 | R4 |
| weak, unknown, none | R5 | R5 | R5 |

Rows are sorted by tier, then by the number of S-LDSC-flagged TFs on the peak (more first), the number of regulon-flagged TFs, ReMM, CADD (all descending), then variant key and gene. S only raises a tier: it never creates evidence on its own, because S-LDSC measures where common-variant heritability concentrates.

- **Regulon enrichment (supporting; used only to sort).** The phenotype genes are the scDoRI-universe genes with a strong or moderate pheno. For each TF with at least one regulon gene in the universe, a one-sided hypergeometric test asks whether its regulon (union of TF-to-gene edges over cell types, restricted to the universe) holds more phenotype genes than expected. The p-values are BH-adjusted across TFs; a TF with q below `fdr` is regulon-flagged. It does not change the tier.
- **Random-trait null.** "Any of k declared traits" inflates the chance that some TF is flagged, so the summary reports a null: `null_draws` random sets of k traits (k = number of declared traits), drawn from every trait in `traits.tsv` with the seed above, each counted by the number of distinct S-LDSC-flagged TFs. `null_p = (1 + #(null >= observed)) / (null_draws + 1)`, where `observed` is the flagged-TF count of the declared traits. A small `null_p` means the declared traits flag more TFs than arbitrary trait sets of the same size; a large one means the flagged TFs are what any traits would give, and S then carries little information. The null tests the trait choice, not any variant.

### Regulatory-network second hit

With `second_hit: true`, a row is **eligible** (`second_hit=yes` in the report) when its pheno is `strong` (the panel floor counts) and V holds. That is R1, or R2 through V. The variant then joins the candidates of gene G with a new reason, `regulatory network: peak <peak> -> <G>; <tier>[; S-LDSC TFs <list>] (second hit only)`, at level `supporting` and second-hit-only: it never ranks a gene alone, and can only complete a recessive pair as its weaker allele, like ReMM and non-coding CADD (see "epcr is heuristic"). Rules:

- **Candidate with a CSQ entry for G** (VEP already linked the variant to G): it has its normal gene evidence and the network reason is added only when the allele has no standalone evidence: an allele that already has one (for example ClinVar P/LP, or SpliceAI >= 0.2) keeps that level and shows no network reason. An existing (variant, G) candidate is replaced in place, so a pair never contains the same variant twice.
- **Link-only candidate** (the gene comes from the peak link alone): the evidence carries only the network reason and the ClinVar benign / likely-benign veto. Popmax, consequence and the peak name (in place of an HGVS) are kept for reporting. No CADD, REVEL, AlphaMissense, splice or panel signal is attached, so the variant cannot be anything but supporting and second-hit-only.
- **Exclusions.** The `regulatory_network` evidence class removes the reason entirely. A row also needs one predictor class behind its V that is not excluded: CADD maps to `unsupervised_predictor`, ReMM to `clinically_trained_predictor`, and the reason is dropped when every class behind V is excluded.
- **Bound.** A pair with a strong gene and a network allele is `[supporting, full, strong]`, tier T3 (0.30), times `p_trans_unphased`: at most 0.15 while the phase is unproven, so it cannot displace a T1 or T2 pair and can appear as an alternative pair. If a VCF `PS` / `PID` proves the pair trans it reaches 0.30. Two second-hit-only alleles never pair.
- **Scenario.** `no_regulatory_network` (below) removes the hook, and is reported only when `enabled` and `second_hit` are both true, so other runs keep the same sensitivity table.
- **Eligibility is fixed at base settings.** `second_hit` uses the base phenotype cut-offs and panel floor; scenarios that change `pheno_rank` or the panel do not recompute it (the class exclusions are the only re-gating).
- **Phase.** The read-back step (`phase --prepare`) scans the VCF without the hook, so a network allele is not read-back phased and its pairs stay at the unphased factor unless the VCF phases them (`PS` / `PID`).

### Regulatory-network outputs

`*.regulatory.tsv` has one row per (variant, linked gene), in the sort order above, with `chrom_prefix` applied to the contig. Columns:

| Column | Meaning |
|---|---|
| `chrom`, `pos`, `ref`, `alt` | the normalised call (the CSV instead uses the original VCF record when normalisation changed it); only the `chrom_prefix` is shared |
| `gene` | the linked gene G |
| `tier` | R1 to R5 |
| `pheno`, `panel` | G's phenotype strength, and `panel=<name>` (with `(phenotype floor)` when the floor raised it) |
| `second_hit` | `yes` when the row is **eligible** (strong pheno and V), not necessarily used: pairs form only in recessive or unknown-MOI genes, and a hom/hemi call, a call with a missing allele, or a gene without a partner allele never produces a pair. The column and `second_hit_rows` are computed whatever `regnet.second_hit` says; that key only controls whether the hook runs |
| `peak` | the peak that gave the best tier |
| `popmax`, `cadd`, `remm` | popmax over `af.populations`; CADD PHRED (maximum over CSQ entries); ReMM from Genomiser. Blank when absent. |
| `v` | `yes` when V holds |
| `n_tfs` | TFs on the peak |
| `sldsc_tfs` | `TF(trait:q,...);...`: the flagged TFs on the peak with their trait and q-value |
| `regulon_tfs` | `TF(q);...`: the regulon-flagged TFs on the peak |
| `consequence` | the record's consequence terms, `&`-joined |

`*.regulatory.summary.tsv` is a two-column `key`, `value` table (lists `;`-joined): `hits` (captured records), `hits_qc_pass`, `rows`, `rows_R1` to `rows_R5`, `second_hit_rows`, `linked_genes_with_hits`, `linked_genes_unranked` (hit genes Exomiser's phenotype-only pass does not rank, before the panel floor; they are `unknown` and so R5, which is where a gene-symbol mismatch between scDoRI and Exomiser shows up), `declared_traits` (`id:label`), `flagged_tfs_<trait>`, `observed_flagged_tfs`, `null_draws`, `null_median`, `null_p`, `pheno_genes_in_universe`, `regulon_flagged_tfs` (`TF(q)`).

### Regulatory-network limitations

- **Coverage.** scDoRI models a fixed gene, TF and peak universe (about 4,000 genes, 500 TFs, 79 k peaks; S-LDSC covers 112 of the TFs), and only a subset of genes has a linked peak. Absence of a link is not evidence against a gene. A gene whose symbol differs between scDoRI and Exomiser is `unknown`, so its rows are R5; the summary counts them (`linked_genes_unranked`).
- **No known distal allele.** A published distal regulatory allele (an intergenic variant 44 kb upstream of a mosaic variegated aneuploidy gene; Ochiai 2014) lies outside every scDoRI peak, so this method would not have found it. The cis window is the peak, not a distance.
- **Common to rare transfer.** S-LDSC measures where common-variant heritability concentrates. Using it to up-weight a rare variant rests on the convergence of common and Mendelian disease genes (Freund 2018), not on direct evidence, which is why S only raises a tier.
- **Missing tissues.** The atlas has no kidney, placenta or gonad, so phenotypes whose biology sits in an absent cell type are not covered.
- **Trait choice and multiplicity.** A TF is flagged when it passes in any declared trait. BH runs within a trait across TFs, not across traits, so more declared traits give more flags; the random-trait null reports this, and it only calibrates the count of flagged TFs. A phenotype with no matching trait gets no S at all.
- **No motif-disruption estimate.** Whether the variant changes a TF motif is not scored. Every rare variant in a peak is otherwise alike, so V (CADD or ReMM) is required for the upper tiers and the second hit.
- **Permissive regulons.** A regulon is the union of TF-to-gene edges over all cell types, so regulons are large and the enrichment is used only to sort within a tier.
- **Peak coordinates.** The 0- or 1-based convention is unverified, so overlap is inclusive and errs by at most 1 bp toward inclusion. The link sign is ignored: any listed peak counts as linked.
- **Eligibility and phase.** Second-hit eligibility uses the base phenotype cut-offs in scenarios, and network pairs are not read-back phased (see above).
- **Post-hoc design.** The module was designed after Model 2's score and Model 3's result (see "Model history"). Thresholds reuse existing constants (CADD 20, ReMM 0.963, `pheno_rank`) or are fixed in this section, and the config is frozen before the run.

Citations: S-LDSC, Finucane 2015 (doi:10.1038/ng.3404); cell-type S-LDSC with a one-sided coefficient test and FDR across annotations, Finucane 2018 (doi:10.1038/s41588-018-0081-4); sc-linker, Jagadeesh 2022 (doi:10.1038/s41588-022-01187-9); Mendelian genes near GWAS loci of phenotype-matched traits, Freund 2018 (doi:10.1016/j.ajhg.2018.08.017); non-coding variants as an in-trans second hit in phenotype-matched genes, Ellingford 2022 (doi:10.1186/s13073-022-01073-3); a distal regulatory MVA allele, Ochiai 2014 (doi:10.1073/pnas.1317008111).

### Sensitivity scenarios

`--sensitivity` writes `*.sensitivity.tsv` (`scenario, row, gene, finding_type, tier, epcr, variant_1, variant_2`): the selected rows under ten scenarios (eleven with the regulatory-network hook), each of which replaces only its own setting. Every scenario is reported; none is chosen afterwards.

| Scenario | Change | Declared |
|---|---|---|
| `base` | none | before Model 2's run |
| `no_known_variant` | `exclude_evidence_classes` + `known_variant` | before Model 2's run |
| `no_known_no_clinical_predictor` | + `known_variant`, `clinically_trained_predictor` | before Model 2's run |
| `pheno_rank_10_100`, `pheno_rank_50_500` | `pheno_rank` strong/moderate cut-offs 10/100 and 50/500 (default 20/200) | before Model 2's run |
| `p_trans_0.15`, `p_trans_0.9` | `p_trans_unphased` (brackets the 0.5 default) | before Model 2's run |
| `no_gene_panel` | panel off (no floor, no deep search) | for Model 3, before its run |
| `one_finding_per_gene` | `alternative_pairs: {max_per_gene: 1, min_phenotype: strong}` | for Model 3, before its run |
| `no_read_backed_phase` | `phase_sources: [multiallelic, ps, pid]` | for Model 3, before its run |
| `no_regulatory_network` | `exclude_evidence_classes: [regulatory_network]`; reported only with `regnet.enabled` and `regnet.second_hit` | for Model 4, before its run |

- With AlphaMissense as the missense tool, `no_known_no_clinical_predictor` is close to `no_known_variant`: removing `clinically_trained_predictor` now only removes the REVEL fallback and ReMM. The informative ablation of the missense tool is `unsupervised_predictor`, which is not a declared scenario (set it in `exclude_evidence_classes` for a separate run).
- `no_gene_panel`, `one_finding_per_gene` and `no_read_backed_phase` equal `base` when their feature is off (no panel; `max_per_gene: 1`; phasing disabled, so there is no `whatshap` source). `no_regulatory_network` equals `base` when no network row was used in a pair.

### Not covered

- No cis-regulatory window. A second allele far upstream, like the intergenic variant 44 kb upstream of BUB1B in Ochiai 2014 (doi:10.1073/pnas.1317008111), cannot be found: VEP links up- and downstream variants to a gene only within 5 kb, entries without a gene symbol are dropped, and the deep search excludes up/downstream and plain intronic variants (user decision 2026-10-20: phase only).
- The optional `regnet` module adds a search bounded by scDoRI peaks only (see "Regulatory-network limitations"); it does not change the bullet above.
- SVs and CNVs, so an SV second allele cannot be found either. A single hit in a recessive gene may have an undetected second allele; that is why `single_hit` is its own category.
- Read-backed phase blocks from short reads are short. WhatsHap links two alleles only through overlapping reads and intermediate het calls, so a pair several kb apart may stay `unknown`. No setting is tuned to change that.
- No statistical or hybrid phasing (SHAPEIT, Beagle with a population panel; Choi 2018, doi:10.1371/journal.pgen.1007308). A rare variant absent from the reference panel gets an arbitrary statistical phase, so it is deferred.
- A homozygous deep-search candidate is never reported: panel candidates are second-hit-only and never form a single or homozygous finding.
- Missense panel candidates have no BP4 filter (the AlphaMissense BP4 cut-offs are not applied).
- REVEL is now used only as a fallback, when AlphaMissense has no score.
- Repeat indels are handled only at canonical splice sites (`splice_shifted_out`). Repeat indels elsewhere (for example frameshifts in homopolymers) are not specially treated.
- Mitochondrial variants are skipped.
- SpliceAI2 has no ClinGen calibration; the SpliceAI v1 cut-off (0.2, supporting) is reused.
- The non-coding second-hit signals and the panel candidates are not calibrated for pathogenicity. They never rank a gene alone and count only as the weaker allele of a recessive pair, discounted by `p_trans_unphased` (see "epcr is heuristic").
- Secondary findings cover rare variants only (popmax <= 1e-2, and <= 1e-4 for a dominant het). Common founder P/LP alleles such as an HFE p.C282Y homozygote are not reported, so the `hfe_c282y_hom` rule rarely fires (user decision 2026-10-18).
- The same variant can appear once as a primary and once as a secondary row in rare cases, for example a recessive SF gene pair that is secondary while one of its alleles is also a primary single hit in the same gene.
- Secondary rows are not kept below the primary rows in epcr. The challenge text says secondary findings do not affect the automated score, but the Space's public scorer (`evaluation.py`, read 2026-10-20) sorts every row by epcr and ignores `finding_type`. An SF row at 0.30 or 0.20 can therefore tie with or outrank a primary pair, and it adds a false positive to F-max. Models 2 and 3 have no secondary rows.
- Hom and hemi calls are always frequency-gated at the recessive cut-off, also in a dominant gene.

### Model history

Each submitted model is a code version plus `config.yaml`. The patient-specific run settings (sample, input VCF, HPO terms, sex) are not in this repository; `config.yaml` holds placeholders for them.

- **Model 2** used REVEL as the single missense tool, no excluded regions, no gene panel, one finding per gene and no read-backed phasing. With today's code the closest settings are `missense_tool: revel`, `exclude_regions: []`, the gene panel off, `alternative_pairs.max_per_gene: 1` and `phasing.enabled: false`. The output can still differ in the one pair-choice edge case (a weaker proven-trans pair with a higher epcr now beats a stronger unphased pair), and in what has no setting: the splice-shift rule, and the phase graph (transitive relations, conflicts, cis collapse) in genes with three or more phase-linked alleles.
- **Model 3** = the settings in `config.yaml`: AlphaMissense as the missense tool, the splice-shift rule, the MHC exclusion, the phase graph with cis collapse, read-backed phasing (`phase`), the MVA gene panel (phenotype floor and deep second-hit search), and alternative pairs with a shared budget, with the three new sensitivity scenarios. Its settings were fixed before its phase and prioritise runs; no value was changed after seeing the output.
- The MVA gene panel, the deep second-hit search and the alternative pairs were added after Model 2's challenge score, so Model 3 is not an independent evaluation.
- **Model 4 (regulatory network)** = Model 3 plus the `regnet` module (see "Regulatory network (`regnet`)"): scDoRI enhancer-peak links, TF regulons and S-LDSC trait enrichment, a regulatory report and a supporting, second-hit-only allele in the pair ranking, with the `no_regulatory_network` scenario. It was designed after the Model 2 and Model 3 results, so it is not an independent evaluation either. Its trait list and every other setting were fixed before its run; no value was changed after seeing the output. The module is off by default, so Model 2 and 3 settings reproduce their outputs unchanged.

## Downstream interpretation caveats

**Missing evidence types**

- **SVs/CNVs, repeat expansions, SMN1 copy number, mtDNA heteroplasmy and UPD/imprinting are not covered by the SNV/indel run.** The optional SV run gives VEP consequences (including regulatory overlaps) for supplied calls, with no dosage-sensitivity scores: DosageSensitivity crashes on SVs when regulatory is on. In TAR the null allele is usually a 1q21.1 deletion. No CNV in this output does not mean no CNV.

**Gaps and blanks**

- **A blank is NA, not benign.** SpliceAI2 (the default) covers SNVs and indels in GENCODE v48 protein-coding genes; its summary-score cut-offs (0.1 / 0.25 / 0.5) are the authors', with no ClinGen calibration. SpliceAI (Ensembl, off by default) covers SNVs on MANE transcripts only, with masked scores, matched by symbol; 50 bp is its delta window. AlphaMissense, REVEL, popEVE, MechPredict and MaveDB cover missense only. For intronic, UTR or indel variants in candidate genes, consider running SpliceAI or Pangolin locally, and see the Genomiser ReMM second-hit signal above.

**Interpreting predictor scores**

- **The predictors are correlated.** ClinGen applies one calibrated predictor for PP3/BP4 ([Pejaver 2022](https://doi.org/10.1016/j.ajhg.2022.10.013)), using thresholds derived on disease genes. Genome-wide use inflates the yield ([Stenton 2024](https://doi.org/10.1016/j.gim.2024.101213)).
- **Known-variant leakage:** ClinVar, CLIN_SIG/PUBMED, REVEL's training data and Exomiser's PP5 all carry curated knowledge. Use `fields.tsv` to rank with and without `known_variant` / `clinically_trained_predictor`. Consider a ClinVar release from before the challenge.
- **popEVE matches on allele only:** require `popEVE_gene == SYMBOL`. It is comparable across genes but is not a PP3 strength, and its calibration population is UK Biobank (mostly European).
- **LOEUF, pHaplo and MechPredict assume a dominant model.** Recessive genes aren't LoF-depleted, so these are not evidence for biallelic candidates. LOEUF `match_by=gene` reports the most constrained transcript.
- **MaveDB scores are raw per score set**, with no shared scale or direction. Keep `MaveDB_urn` and use only calibrated sets.

**G2P limitations** (the plugin is kept by user decision)

- `G2P_complete` depends on how variants fall into forks.
- It counts coding consequence types only, so it can't pair a null with a UTR or intronic hypomorph.
- The default `af_biallelic=0.005` hides alleles like the ~3% TAR SNP.
- The default confidence levels exclude "limited".
- Use it as a gene-level label, not as a pair finder.

**Compound hets in an unphased singleton**

- A pair is only a candidate, and about 50% of random pairs are in cis.
- Build pairs downstream with a phase state: trans-certain (split 1/2 sites, opposite `PS`/`PGT`), cis-certain, or unknown. `prioritise` does this with its phase graph (see "Phase of pairs").
- Record the distance between the pair. Read-backed phasing (WhatsHap) is feasible only for close pairs if a BAM is available; `phase` runs it on the candidate genes.
- Estimate the background count of qualifying pairs, ideally on GIAB HG002 run through the same pipeline.
- Keep a "single hit in a biallelic gene" category (cf. [Matsuura 2006](https://doi.org/10.1002/ajmg.a.31069)).

**Sample and data quality**

- **QC:** the original QUAL/FILTER/FORMAT fields (AD, DP, GQ, PS) are kept. The ranking step must use them, because artefacts look "novel".
- **Mosaicism:** `prepare-patient-vcf.sh -q` is deliberately not used, since it would drop hets below 15% VAF. Bulk WGS has little power for random mosaic aneuploidy, so a karyotype or PCS assay is the decisive test.
- **Sample source:** the sex estimate is informational, so confirm sex. Confirm the DNA is germline; tumour DNA would add somatic variants and 11p15 LOH.
- **Allele frequencies:** `MAX_AF` includes small and bottlenecked populations. gnomAD nhomalt matters for recessive models; record the proband's ancestry.
- **Version skew:** the cache is Ensembl 116, while PhenotypeOrthologous was built on 112, SpliceAI on MANE v1.4 and LOEUF on gnomAD v4. The smoke report quantifies coverage.

## Scientific references

- No allele-frequency or consequence pre-filter at annotation time: filtering later can be undone, and a pre-filter could hide a low-frequency regulatory allele in trans with a rare null (TAR syndrome, [Albers 2012](https://doi.org/10.1038/ng.1083)), a second allele with reduced expression and no coding mutation (PCS/MVA, [Matsuura 2006](https://doi.org/10.1002/ajmg.a.31069)), or hypomorphic BUBR1 alleles ([Sieben 2020](https://doi.org/10.1172/JCI126863)). These precedents support a second-hit search within candidate genes; they do not make genome-wide non-coding ranking powerful.
- REVEL is ClinGen-calibrated up to PP3_Strong for missense ([Pejaver 2022](https://doi.org/10.1016/j.ajhg.2022.10.013)); so is AlphaMissense ([Bergquist 2025](https://doi.org/10.1016/j.gim.2025.101402)). popEVE: [Orenbuch 2025](https://doi.org/10.1038/s41588-025-02400-1).
- SpliceAI use follows the ClinGen splicing recommendations ([Walker 2023](https://doi.org/10.1016/j.ajhg.2023.06.002)). G2P allelic requirements: [Thormann 2019](https://doi.org/10.1038/s41467-019-10016-3). MechPredict: [Badonyi & Marsh 2024](https://doi.org/10.1371/journal.pone.0307312).
- MVA gene panel: BUB1B ([Hanks 2004](https://doi.org/10.1038/ng1449)), CEP57 ([Snape 2011](https://doi.org/10.1038/ng.822)), TRIP13 ([Yost 2017](https://doi.org/10.1038/ng.3883)). Non-coding or hypomorphic second BUB1B alleles: [Rio Frio 2010](https://doi.org/10.1056/NEJMoa1006565), [Ochiai 2014](https://doi.org/10.1073/pnas.1317008111), [Sieben 2020](https://doi.org/10.1172/JCI126863). Read-backed phasing: WhatsHap ([Patterson 2015](https://doi.org/10.1089/cmb.2014.0157)).
- Regulatory network: S-LDSC ([Finucane 2015](https://doi.org/10.1038/ng.3404)); cell-type S-LDSC ([Finucane 2018](https://doi.org/10.1038/s41588-018-0081-4)); sc-linker ([Jagadeesh 2022](https://doi.org/10.1038/s41588-022-01187-9)); common and Mendelian disease genes converge on phenotype-matched traits ([Freund 2018](https://doi.org/10.1016/j.ajhg.2018.08.017)); non-coding second hits in phenotype-matched genes ([Ellingford 2022](https://doi.org/10.1186/s13073-022-01073-3)); a distal regulatory allele upstream of an MVA gene ([Ochiai 2014](https://doi.org/10.1073/pnas.1317008111)).

## Tests

```bash
cd /fast/AG_Bunina/Berk/projects/mva_challenge_26/scripts/vep_pipeline
micromamba run -n vep116 pytest tests -q                 # unit tests and the integration test
micromamba run -n vep116 pytest tests -q -m "not integration"
```

The integration test runs real VEP on the cache bundled with `ensembl-vep` (SNV and SV).
