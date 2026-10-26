# Exomiser 15.1.0 + Genomiser install

GRCh38 only, data release **2602**. Verified by known-answer tests, not by
absence-of-error. Both acceptance gates re-run and passing on 2602.

## Quick start

```bash
source bin/exomiser-env.sh          # micromamba env `exomiser` (OpenJDK 21) + PATH + SSL fix

# 1. normalise the patient VCF (REQUIRED — see "Why normalisation is not optional")
bin/prepare-patient-vcf.sh \
    -i patient.vcf.gz \
    -r data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -o patient.norm.vcf.gz -s PATIENT_01 -c 8

# 2. two-pass analysis
bin/run-exomiser.sh PATIENT_01 patient.norm.vcf.gz patient.phenopacket.yml results/
# -> results/exome/PATIENT_01.{html,genes.tsv,variants.tsv,parquet}
# -> results/genome/PATIENT_01.{html,genes.tsv,variants.tsv,parquet}
```

The phenopacket's `subject.id` **must equal** the VCF sample column name.
Template: `test/genome-smoke.phenopacket.yml`.

## Layout

| Path | What |
|---|---|
| `exomiser-cli-15.1.0/` | distribution (thin jar + `lib/` + `examples/`) |
| `data/2602_hg38/`, `data/2602_phenotype/` | Exomiser reference data, **in use** |
| `data/2512_hg38/`, `data/2512_phenotype/` | previous release, superseded — 55 GB, safe to delete |
| `data/remm/ReMM.v0.4.hg38.tsv.gz` | ReMM — this is what makes it Genomiser |
| `data/whole_genome_SNVs.tsv.gz`, `data/gnomad.genomes.r4.0.indel.tsv.gz` | CADD v1.7 GRCh38 (pre-existing; see `data/filelist`) |
| `etc/application.properties` | site config, absolute paths |
| `etc/analysis-genome-cadd.yml` | the only custom analysis YAML (genome preset + CADD) |
| `etc/analyse-help.txt` | captured `analyse --help` for this exact build |
| `bin/` | env helper, runner, VCF prep, parallel downloader |
| `test/` | hg38 known-answer fixtures + recorded ground truth |
| `data/Homo_sapiens.GRCh38.dna.primary_assembly.fa` | (Ensembl, **unprefixed** contigs) |

## Verified acceptance gates

| Gate | Claim | Result |
|---|---|---|
| 2 — exome | `FGFR2` rank 1 on ClinVar `10:121496701 T>G` + Pfeiffer HPO | ✅ rank 1, combined 0.9452, ACMG **PATHOGENIC** (`PP3_Strong,PP5_Strong,PP4_Moderate`), ACMG disease `OMIM:101600`, ClinVar 13294 ★★ |
| 2 — CADD isolation | CADD configured but absent from exome scoring | ✅ `ALL_PATH=REVEL=0.965,MVP=0.9518,ALPHA_MISSENSE=0.9918` — no CADD |
| 3 — Genomiser | Exomiser's `REMM` **equals** an independent `tabix` read | ✅ both `0.993508` at `chr13:28659088` (POMP 5'UTR, KLICK) |
| 3 — contig naming | `13` vs `chr13` in the VCF must not change results | ✅ byte-identical; Exomiser normalises internally |
| 3 — non-coding evidence | CADD contributes in the genome pass | ✅ `ALL_PATH=CADD=0.99241424,REMM=0.993508` |
| 4 — VCF prep, matching naming | `exome-smoke.vcf.gz` (contig `10`) normalises | ✅ 0/1 REF mismatch, 1 variant out |
| 4 — VCF prep, `chr` naming | `chr13` VCF vs unprefixed reference | ✅ auto-renamed to `13`, normalises, `EXO_ORIG` provenance kept |
| 4 — round trip | prepped VCF reproduces the Gate 3 answer | ✅ identical `CADD=0.99241424,REMM=0.993508` |
| 4 — **wrong-assembly gate** | real hg19 `Pfeiffer.vcf.gz` vs GRCh38 must abort | ✅ 24/37 contig lengths disagree (`1: 249250621 vs 248956422`), rc=1, **no output written** |
| 4 — scaffold tolerance | 14 hg19-era scaffolds (`GL000192.1`…) | ✅ dropped with a note, not a fatal error |

Re-run any gate from `test/`; ground truth is in `test/remm-ground-truth.tsv`
and `test/remm-at-klick.tsv`.

## Why two passes

Cooperstein et al. 2025 (Genome Med, doi:10.1186/s13073-025-01546-1) found Genomiser
alone ranks ~30% fewer *coding* diagnostic variants in the top 30 and runs ~3× longer,
while Exomiser alone cannot recover a compound heterozygote made of one coding plus one
non-coding variant. So `run-exomiser.sh` runs both and you read both.

Demonstrated here on the POMP non-coding variant:

| pass | ALL_PATH |
|---|---|
| exome | *(empty — no evidence at all)* |
| genome | `CADD=0.99241424,REMM=0.993508` |

## Why CADD is in the genome pass only

Exomiser normalises CADD PHRED so 10 → 0.9 and takes the **maximum** across
`pathogenicitySources`, so CADD saturates most exonic variants at ≥0.98 and destroys
discrimination (Cooperstein et al. 2025). But the non-coding gate in
`pathogenicityFilter` is `(CADD ≥ 15) OR (SpliceAI > 0.1) OR (ReMM > 0.914)`, and ReMM's
training set contains **no splice variants at all** — on non-coding ClinVar, CADD v1.6
scored AUPRC 0.160 vs ReMM's 0.035 (Schubach et al. 2023,
doi:10.1093/gigascience/giad024). So CADD covers a real gap there and nowhere else.

## Why normalisation is not optional

Exomiser reads frequency and pathogenicity from its own MVStore keyed on an **exact
match** of `(contig, pos, ref, alt)` against gnomAD-derived, **left-aligned** keys. A
differently-represented indel returns empty allele properties — no frequency, no
pathogenicity, **no warning**. Wrong science, clean exit code.

`bin/prepare-patient-vcf.sh` is derived from `bin/annotate_vcf.sh` but deliberately
**drops `bcftools norm -g/--gff-annot`**, whose documented behaviour is *"Follow HGVS
3'rule and right-align variants in transcripts"* — the opposite of what Exomiser
matches against. It also drops `bcftools csq` (Exomiser annotates with Jannovar) and
`bcftools annotate` (Exomiser reads its own ClinVar/frequency DBs), and changes
`--check-ref ws` → `w` plus a hard abort, because `s` silently rewrites mismatching REF
and would launder a wrong-assembly VCF into a plausible-looking one.

## Gotchas found the hard way

1. **`SSL_CERT_FILE` in the login environment points at `/etc/ssl/certs/ca-bundle.crt`,
   which does not exist on this host.** All HTTPS dies with `curl: (77)`. The real
   bundle is `/etc/ssl/certs/ca-certificates.crt`. `bin/exomiser-env.sh` repairs it.
   *This affects any tool you run here that does HTTPS, not just Exomiser.*
2. **Every score path in `application.properties` is opened EAGERLY at startup**, even
   by analyses that never use it. A missing ReMM file kills a `phenotype_only` run.
3. **`--assembly` is mandatory whenever `--vcf` is given.** Omit it and the sample
   defaults to hg19 → `UnsupportedGenomeAssemblyException`. Accepts `GRCh38` or `hg38`.
4. **The preset is `phenotype_only`** (underscore), not `phenotype-only` as docs say.
5. **15.x requires the `analyse` subcommand.** `java -jar … --analysis x.yml` (the
   13.x/14.x form in most tutorials) fails.
6. **Every shipped example is hg19**, including
   `NA19722_601952_…_5UTR_38.yml` — the `_38` labels the variant, not the assembly.
   There is no official hg38 test data anywhere; `test/` holds ours.
7. **ReMM v0.4 hg38 contigs are `chr`-prefixed** (`chr1`…`chrY`), contrary to the docs.
8. **Two different files are both named `ReMM.v0.4.hg38.tsv.gz.`** Zenodo record
   6138102 is a buggy build; 6576087 is correct. Only the MD5 tells them apart —
   ours is verified `2f8c8fcc75ea29f6798621d0446fe194`.
9. **`bgzip` will not clobber an existing `.gz` and still exits 0**, leaving the old
   content in place. With a fixed temp dir and no `set -e` that is cross-patient
   contamination; `prepare-patient-vcf.sh` uses a fresh `mktemp -d` per run.
10. **Exomiser 15.0.0 recalibrated the combined score.** Any cutoff inherited from a
    13.x/14.x pipeline is invalid.
11. The genome DB is `2512_hg38_genome.mv.db` — **MVStore, opened read-only with a
    shared lock**, so many concurrent jobs can share `data/` safely. No `.lock.db`.

## Data release history

Now on **2602** (Feb 2026). Both gates were re-run after the switch; ranks and
pathogenicity evidence are identical to 2512. The only movement is the exome
phenotype score (0.8127 → 0.8194, combined 0.9452 → 0.9486), consistent with 2602
being primarily a phenotype release (OMIM 2026-02, HPO 2026-01, ClinVar 2026-02).

`etc/application.properties.2512.bak` is the pre-switch config. `data/2512_*` is
superseded and can be removed to reclaim ~55 GB.

## Input requirements

- **GRCh38** only. `prepare-patient-vcf.sh` aborts on a contig-length mismatch against
  the reference `.fai`, on a REF-mismatch rate above 2%, and if `bcftools norm` cannot
  resolve the VCF against the reference at all.
- **Contig naming is handled for you.** `-C auto` (the default) rewrites the VCF to the
  reference's own style via `bcftools --rename-chrs`. Our reference is Ensembl, so a
  GATK/UCSC `chr13` VCF is rewritten to `13` *before* normalisation — which is
  mandatory, because `norm` resolves contigs by exact name and otherwise dies with
  `The sequence "chr13" was not found`.
- **Non-primary contigs are dropped, not fatal.** `chrUn_*`, `*_alt`, `*_random`,
  `HLA-*` and decoys absent from the primary assembly are removed with a note.
  A missing *primary* chromosome is still a hard error.
- Sample column name **must equal** the phenopacket `subject.id`.
- Multi-sample needs `--ped`, and the **proband must be the first sample column**
  (a known bug otherwise yields zero gene-variant scores for comp-het-compatible genes).
- No gVCFs — genotype them first; `<NON_REF>` blocks match nothing.

## Not done yet

- Trio / multi-sample + PED path untested.
- The `-q` quality prefilter (VAF 15–85%, GQ ≥ 20) is implemented but has not been
  exercised on a VCF that actually carries AD/DP/GQ.
- No SLURM wrapper (deferred).
