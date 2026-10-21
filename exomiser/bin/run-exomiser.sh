#!/usr/bin/env bash
# Two-pass Exomiser/Genomiser runner.
#
#   run-exomiser.sh <sample-id> <vcf.gz> <phenopacket.yml> <outdir> [assembly]
#
# Two passes because Cooperstein et al. 2025 (Genome Med) found Genomiser
# alone ranks ~30% fewer coding diagnostic variants in the top 30 and takes
# ~3x longer, while Exomiser alone cannot recover compound heterozygotes made
# of one coding plus one non-coding variant. Run both; read both.
#
# Deliberately no `set -e`: the conditional-append idiom used for optional
# flags returns non-zero when the condition is false, which would kill the
# job under errexit.
set -uo pipefail

source /fast/AG_Bunina/Berk/tools/exomiser/bin/exomiser-env.sh || {
  echo "ERROR: environment setup failed; not running Exomiser" >&2; exit 1; }

SAMPLE="${1:?sample id required}"
VCF="${2:?bgzipped vcf required}"
PPKT="${3:?phenopacket required}"
OUTDIR="${4:?output directory required}"
ASSEMBLY="${5:-GRCh38}"

for f in "$VCF" "$PPKT" "$EXO_JAR" "$EXO_CONF"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
done

# Every score path named in application.properties is opened EAGERLY at
# startup, even by analyses that never use it. A missing ReMM/CADD file is a
# hard startup failure, so check up front with a clear message.
while IFS= read -r p; do
  [[ -e "$p" ]] || { echo "ERROR: configured score file absent: $p" >&2; exit 1; }
done < <(grep -oP '(?<=^exomiser\.hg38\.(remm-path|cadd-snv-path|cadd-in-del-path)=).*' "$EXO_CONF")

# 15.0.0 refuses to write results into the installation tree, and slurmd/
# Exomiser will not mkdir for you.
mkdir -p "$OUTDIR/exome" "$OUTDIR/genome" || exit 1

FORMATS=HTML,TSV_VARIANT,TSV_GENE,PARQUET

# Size the heaps to the actual allocation. Hardcoding -Xmx12g got a 4.96M-variant
# WGS genome pass SIGKILLed on a session that did not have 12 GiB to give.
HEAP_EXOME="$(exo_heap_gib exome)g"
HEAP_GENOME="$(exo_heap_gib genome)g"
exo_mem_report
echo "heaps: exome=-Xmx$HEAP_EXOME  genome=-Xmx$HEAP_GENOME"
echo "variants in input: $(bcftools index -n "$VCF" 2>/dev/null || echo unknown)"

# ── Pass 1: coding candidates ────────────────────────────────────────────
# Stock shipped preset. CADD is configured but NOT listed in its
# pathogenicitySources, so it contributes nothing here -- verified: ALL_PATH
# shows only REVEL/MVP/ALPHA_MISSENSE.
echo "=== pass 1/2  exome    -> $OUTDIR/exome"
t0=$SECONDS
java -Xmx"$HEAP_EXOME" $(exo_jvm_flags "$SAMPLE.exome") -jar "$EXO_JAR" analyse \
  --analysis "$EXO_EXAMPLES/preset-exome-analysis.yml" \
  --sample "$PPKT" --vcf "$VCF" --assembly "$ASSEMBLY" \
  --output-directory "$OUTDIR/exome" --output-filename "$SAMPLE" \
  --output-format "$FORMATS"
rc1=$?
echo "=== pass 1 rc=$rc1 in $((SECONDS-t0))s"

# ── Pass 2: non-coding (Genomiser) ───────────────────────────────────────
# Site YAML adds REMM + CADD. PASS_ONLY keeps heap near 4 GB on a real
# genome instead of ~12 GB.
echo "=== pass 2/2  genomiser -> $OUTDIR/genome"
t0=$SECONDS
java -Xmx"$HEAP_GENOME" -XX:+UseG1GC $(exo_jvm_flags "$SAMPLE.genome") -jar "$EXO_JAR" analyse \
  --analysis "$EXO_ROOT/etc/analysis-genome-cadd.yml" \
  --sample "$PPKT" --vcf "$VCF" --assembly "$ASSEMBLY" \
  --output-directory "$OUTDIR/genome" --output-filename "$SAMPLE" \
  --output-format "$FORMATS"
rc2=$?
echo "=== pass 2 rc=$rc2 in $((SECONDS-t0))s"
if (( rc2 == 137 )); then
  echo "" >&2
  echo "rc=137 is SIGKILL: the genome pass was killed by a memory limit, not by a" >&2
  echo "Java heap error (that would print java.lang.OutOfMemoryError instead)." >&2
  echo "This session's detected limit: $(exo_mem_limit_gib) GiB. Resubmit via SLURM:" >&2
  echo "    sbatch --mem=48G $EXO_ROOT/bin/exomiser.sbatch $SAMPLE $VCF $PPKT $OUTDIR" >&2
fi

echo
echo "exome  rc=$rc1  -> $OUTDIR/exome/$SAMPLE.genes.tsv"
echo "genome rc=$rc2  -> $OUTDIR/genome/$SAMPLE.genes.tsv"
(( rc1 == 0 && rc2 == 0 )) || exit 1
