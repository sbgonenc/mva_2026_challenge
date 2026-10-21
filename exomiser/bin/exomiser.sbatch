#!/usr/bin/env bash
# Two-pass Exomiser + Genomiser as a SLURM batch job.
#
#   sbatch --mem=48G bin/exomiser.sbatch <sample> <vcf.gz> <phenopacket> <outdir>
#
# #SBATCH values below are DEFAULTS and can be overridden at submit time, e.g.
#   sbatch --mem=64G --time=12:00:00 bin/exomiser.sbatch ...
#
# Why batch rather than an interactive/OOD session: a ~5M-variant genome pass
# needs 8-16 GiB of heap plus room for the ~54 GiB of memory-mapped MVStore
# files, and an interactive session's allocation is usually far smaller. A
# session that cannot supply it SIGKILLs the JVM (exit 137) with no Java error.
#
# No `set -e`: the conditional-flag idiom used downstream returns non-zero when
# a condition is false, which would abort the job under errexit.
#SBATCH --job-name=exomiser
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/exomiser_%j.out
#SBATCH --error=logs/exomiser_%j.err
#SBATCH --export=ALL
set -uo pipefail

EXO_ROOT=/fast/AG_Bunina/Berk/tools/exomiser

SAMPLE="${1:?usage: sbatch bin/exomiser.sbatch <sample> <vcf.gz> <phenopacket> <outdir>}"
VCF="${2:?vcf required}"
PPKT="${3:?phenopacket required}"
OUTDIR="${4:?outdir required}"

echo "host      : $(hostname)"
echo "job       : ${SLURM_JOB_ID:-<none>}"
echo "mem grant : ${SLURM_MEM_PER_NODE:-?} MB"
echo "cpus      : ${SLURM_CPUS_PER_TASK:-?}"
echo "sample    : $SAMPLE"
echo "vcf       : $VCF"

# Node-local scratch for JVM temp, so nothing lands on the shared filesystem.
export TMPDIR="${SLURM_TMPDIR:-/tmp}/exomiser.${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

exec "$EXO_ROOT/bin/run-exomiser.sh" "$SAMPLE" "$VCF" "$PPKT" "$OUTDIR"
