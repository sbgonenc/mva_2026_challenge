#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════
# prepare-patient-vcf.sh -- normalise a patient VCF for Exomiser 15.1.0
#
# Derived from bin/annotate_vcf.sh, but retargeted: Exomiser needs a
# *normalised* VCF, not an *annotated* one. It does its own transcript
# annotation (Jannovar) and reads frequency/pathogenicity/ClinVar from its own
# MVStore databases, so bcftools csq output and transferred INFO columns are
# ignored. What Exomiser cannot recover from is a representation mismatch:
# allele lookup is an EXACT match on (contig, pos, ref, alt) against
# gnomAD-derived, LEFT-ALIGNED keys. A mis-represented indel returns empty
# allele properties -- no frequency, no pathogenicity, no warning.
#
# Differences from annotate_vcf.sh, each deliberate:
#   - DROPS  bcftools norm -g/--gff-annot. Its documented behaviour is
#            "Follow HGVS 3'rule and RIGHT-align variants in transcripts".
#            That is the opposite of the left-aligned representation Exomiser
#            matches against, so it would silently break lookups for exactly
#            the indels that matter most: the ones inside genes.
#   - DROPS  bcftools csq       (Exomiser annotates with Jannovar)
#   - DROPS  bcftools annotate  (Exomiser reads its own ClinVar/frequency DBs)
#   - CHANGES --check-ref ws -> w, and FAILS on a high mismatch rate instead of
#            silently rewriting REF. Silent repair masks a wrong-assembly VCF,
#            which is the worst outcome in a diagnostic setting.
#   - ADDS   an explicit GRCh38 assembly gate (contig lengths + REF-match rate)
#   - ADDS   '*' spanning-deletion removal (matches nothing, wastes a rank slot)
#   - ADDS   per-run temp dir; the fixed /tmp/annotate_vcf/ collides between
#            concurrent jobs and between users on a shared cluster
#   - ADDS   variant counts at every stage, so silent loss is visible
#   - KEEPS  the getopts interface style, --threads plumbing, and -m - / --atomize
# ═══════════════════════════════════════════════════════════════════════════
set -euo pipefail

EXO_ROOT=/fast/AG_Bunina/Berk/tools/exomiser
# shellcheck source=/dev/null
[[ -r "$EXO_ROOT/bin/exomiser-env.sh" ]] && source "$EXO_ROOT/bin/exomiser-env.sh" || \
  export PATH="/fast/AG_Bunina/Berk/bin:$PATH"

in_vcf=""; reference=""; out_vcf=""; sample=""
n_cores=1; tmp_root="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
qual_filter=0; contig_style="auto"; max_ref_mismatch_pct=2

usage() {
  cat <<'USAGE'
Usage: prepare-patient-vcf.sh -i <in.vcf[.gz]> -r <GRCh38.fa> -o <out.vcf.gz> [options]

Required:
  -i FILE   Input patient VCF (plain or bgzipped)
  -r FILE   GRCh38 reference FASTA (must have a .fai; samtools faidx it if not)
  -o FILE   Output bgzipped VCF (.vcf.gz); a .tbi is written alongside

Options:
  -s NAME   Expected sample name. Verified against the VCF and reported; this
            must equal the phenopacket subject.id or Exomiser errors out.
  -c INT    Worker threads for bcftools            [1]
  -t DIR    Temp parent dir  [$SLURM_TMPDIR, else $TMPDIR, else /tmp]
  -q        Apply the Cooperstein et al. 2025 quality prefilter:
            drop hets outside 15-85% VAF and any genotype with GQ < 20.
            Skipped with a warning if AD/DP/GQ are absent.
  -C STYLE  Contig naming: auto | keep | chr | nochr   [auto]
            'auto' rewrites the VCF to match the REFERENCE's own style, which is
            what bcftools norm requires -- an Ensembl reference uses '13' while a
            GATK/UCSC VCF uses 'chr13', and norm then dies with
            'The sequence "chr13" was not found'. Renaming uses
            bcftools --rename-chrs, never sed on the header.
  -m PCT    Max tolerated REF-mismatch percent before aborting  [2]
  -h        This help
USAGE
  exit 1
}

while getopts "i:r:o:s:c:t:C:m:qh" opt; do
  case $opt in
    i) in_vcf="$OPTARG";;
    r) reference="$OPTARG";;
    o) out_vcf="$OPTARG";;
    s) sample="$OPTARG";;
    c) n_cores="$OPTARG";;
    t) tmp_root="$OPTARG";;
    C) contig_style="$OPTARG";;
    m) max_ref_mismatch_pct="$OPTARG";;
    q) qual_filter=1;;
    h) usage;;
    *) usage;;
  esac
done

[[ -n "$in_vcf" && -n "$reference" && -n "$out_vcf" ]] || { echo "ERROR: -i, -r and -o are required" >&2; usage; }
for f in "$in_vcf" "$reference"; do
  [[ -f "$f" ]] || { echo "ERROR: no such file: $f" >&2; exit 1; }
done
[[ -f "${reference}.fai" ]] || { echo "ERROR: missing ${reference}.fai -- run: samtools faidx '$reference'" >&2; exit 1; }
[[ "$out_vcf" == *.vcf.gz ]] || { echo "ERROR: -o must end in .vcf.gz (got '$out_vcf')" >&2; exit 1; }
case "$contig_style" in auto|keep|chr|nochr) ;; *) echo "ERROR: -C must be auto|keep|chr|nochr" >&2; exit 1;; esac
# Check the tools RUN, not just that they exist: a binary with unresolvable
# shared libraries satisfies `command -v` and then fails mid-run.
if declare -F exo_check_tools >/dev/null 2>&1; then
  exo_check_tools || exit 1
else
  for t in bcftools bgzip tabix; do
    command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH" >&2; exit 1; }
    "$t" --version >/dev/null 2>&1 || {
      echo "ERROR: $t at $(command -v "$t") will not run:" >&2
      "$t" --version 2>&1 | head -2 | sed 's/^/       /' >&2; exit 1; }
  done
fi

# Per-run temp dir. The original's fixed /tmp/annotate_vcf/ is a cross-user
# collision (same class of bug as Exomiser issue #541's shared /tmp/spring.log),
# and its unconditional `rm -fr` would delete a concurrent job's scratch.
TMP="$(mktemp -d "${tmp_root%/}/exo-prep.XXXXXXXX")"
cleanup() { rm -rf "$TMP"; }
trap cleanup EXIT
mkdir -p "$(dirname "$out_vcf")"

say() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }
nvar() { bcftools index -n "$1" 2>/dev/null || bcftools view -H "$1" | wc -l; }

say "input      : $in_vcf"
say "reference  : $reference"
say "output     : $out_vcf"
say "temp       : $TMP"
say "threads    : $n_cores"

# ── 0. Ensure a bgzipped, indexed input we can stream and count ───────────
if [[ "$in_vcf" == *.gz ]]; then
  cp -f "$in_vcf" "$TMP/in.vcf.gz"
else
  bgzip -c -@ "$n_cores" "$in_vcf" > "$TMP/in.vcf.gz"
fi
bcftools index -f --tbi --threads "$n_cores" "$TMP/in.vcf.gz" 2>/dev/null || \
  bcftools index -f --csi --threads "$n_cores" "$TMP/in.vcf.gz"

n_in=$(nvar "$TMP/in.vcf.gz")
say "variants in : $n_in"
[[ "$n_in" -gt 0 ]] || { echo "ERROR: input has no variant records" >&2; exit 1; }

# ── 1. Sample names ──────────────────────────────────────────────────────
mapfile -t samples < <(bcftools query -l "$TMP/in.vcf.gz" 2>/dev/null || bcftools view -h "$TMP/in.vcf.gz" | tail -1 | cut -f10-)
say "sample(s)   : ${samples[*]:-<none>}"
if [[ ${#samples[@]} -eq 0 ]]; then
  echo "ERROR: VCF has no sample columns; Exomiser needs genotypes" >&2; exit 1
fi
if [[ ${#samples[@]} -gt 1 ]]; then
  say "NOTE: multi-sample VCF. Exomiser requires a PED (--ped) and the PROBAND"
  say "      MUST be the first sample column, or gene-variant scores can come"
  say "      back as 0 for compound-het-compatible genes."
fi
if [[ -n "$sample" ]]; then
  printf '%s\n' "${samples[@]}" | grep -qxF "$sample" || {
    echo "ERROR: expected sample '$sample' not in VCF (${samples[*]})." >&2
    echo "       The phenopacket subject.id must match one of these exactly." >&2
    exit 1; }
  say "sample check: '$sample' present"
fi

# ── 2. CONTIG-NAME RECONCILIATION (must precede norm) ────────────────────
# bcftools norm resolves every contig against the reference by EXACT name. A
# 'chr13' VCF against an Ensembl '13' reference dies mid-normalisation with
# `The sequence "chr13" was not found`. Reconcile first, not afterwards.
say "--- contig-name reconciliation ---"
cut -f1 "${reference}.fai" | sort -u > "$TMP/ref.names"
bcftools view -H "$TMP/in.vcf.gz" | cut -f1 | uniq | sort -u > "$TMP/vcf.names"
ref_chr=0; vcf_chr=0
grep -qx 'chr[0-9XYM]*' "$TMP/ref.names" && ref_chr=1
grep -qx 'chr[0-9XYM]*' "$TMP/vcf.names" && vcf_chr=1
say "reference style: $( [ "$ref_chr" -eq 1 ] && echo 'chr-prefixed' || echo 'unprefixed' )"
say "VCF style      : $( [ "$vcf_chr" -eq 1 ] && echo 'chr-prefixed' || echo 'unprefixed' )"

# Resolve EACH VCF contig against the names that actually exist in the
# reference, rather than applying a blanket chr-prefix rule. A global rule
# cannot express the mitochondrial case: UCSC/GATK call it 'chrM', Ensembl
# calls it 'MT', so stripping 'chr' yields 'M' which exists in neither.
# Exomiser has a MITOCHONDRIAL inheritance mode, so chrM must be MAPPED, not
# dropped.
if [[ "$contig_style" != "keep" ]]; then
  : > "$TMP/rename.txt"
  while IFS= read -r c; do
    [[ -n "$c" ]] || continue
    base="${c#chr}"
    cands=()
    case "$contig_style" in
      chr)   cands=("chr$base" "$base" "$c") ;;
      nochr) cands=("$base" "chr$base" "$c") ;;
      auto)  if [[ "$ref_chr" -eq 1 ]]; then cands=("chr$base" "$base" "$c")
             else cands=("$base" "chr$base" "$c"); fi ;;
    esac
    # Mitochondrial synonyms, appended so an exact hit always wins first.
    case "$base" in
      M|MT) if [[ "$ref_chr" -eq 1 ]]; then cands+=("chrM" "chrMT" "MT" "M")
            else cands+=("MT" "M" "chrM" "chrMT"); fi ;;
    esac
    resolved=""
    for cand in "${cands[@]}"; do
      if grep -qxF -- "$cand" "$TMP/ref.names"; then resolved="$cand"; break; fi
    done
    # Unresolvable contigs are left alone; the orphan check below classifies them.
    if [[ -n "$resolved" && "$resolved" != "$c" ]]; then
      printf '%s\t%s\n' "$c" "$resolved" >> "$TMP/rename.txt"
    fi
  done < <(bcftools view -h "$TMP/in.vcf.gz" | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p')

  if [[ -s "$TMP/rename.txt" ]]; then
    say "renaming $(wc -l < "$TMP/rename.txt") contigs to match the reference"
    awk '{printf "        %s -> %s\n", $1, $2}' "$TMP/rename.txt" | head -6
    bcftools annotate --rename-chrs "$TMP/rename.txt" --threads "$n_cores" \
      -O z -o "$TMP/renamed.vcf.gz" "$TMP/in.vcf.gz"
    mv -f "$TMP/renamed.vcf.gz" "$TMP/in.vcf.gz"
    bcftools index -f --tbi --threads "$n_cores" "$TMP/in.vcf.gz" 2>/dev/null || \
      bcftools index -f --csi --threads "$n_cores" "$TMP/in.vcf.gz"
    bcftools view -H "$TMP/in.vcf.gz" | cut -f1 | uniq | sort -u > "$TMP/vcf.names"
  else
    say "contig names already match the reference; nothing to rename"
  fi
fi

# Every data-bearing contig must now exist in the reference, or norm dies.
# Classify the orphans, because the two cases need opposite handling:
#   - a PRIMARY chromosome missing  => naming or assembly problem, must abort
#   - only scaffolds/alts/decoys    => normal for a GRCh38 VCF against an
#     Ensembl primary assembly; drop them and carry on. Exomiser only analyses
#     primary chromosomes anyway, so aborting here would block real patients.
comm -23 "$TMP/vcf.names" "$TMP/ref.names" > "$TMP/orphan.names" || true
if [[ -s "$TMP/orphan.names" ]]; then
  grep -E '^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$' "$TMP/orphan.names" > "$TMP/orphan.primary" || true
  grep -vE '^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$' "$TMP/orphan.names" > "$TMP/orphan.other" || true

  if [[ -s "$TMP/orphan.primary" ]]; then
    echo "ERROR: PRIMARY chromosomes in the VCF are absent from the reference," >&2
    echo "       even after name resolution:" >&2
    sed 's/^/    /' "$TMP/orphan.primary" >&2
    echo "       Reference contigs available (primary only):" >&2
    grep -E '^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$' "$TMP/ref.names" \
      | paste -sd' ' | fold -s -w 70 | sed 's/^/           /' >&2
    echo "       Either the reference is a different assembly, or this contig has a" >&2
    echo "       naming convention not covered here. Force one with -C chr|nochr, or" >&2
    echo "       use a reference matching the VCF." >&2
    exit 1
  fi

  n_other=$(wc -l < "$TMP/orphan.other")
  say "NOTE: $n_other non-primary contig(s) absent from the reference; dropping them"
  head -5 "$TMP/orphan.other" | sed 's/^/        /'
  awk 'BEGIN{OFS="\t"} {print $1,1,$2}' "${reference}.fai" > "$TMP/keep.regions"
  bcftools view -R "$TMP/keep.regions" --threads "$n_cores" \
    -O z -o "$TMP/subset.vcf.gz" "$TMP/in.vcf.gz"
  mv -f "$TMP/subset.vcf.gz" "$TMP/in.vcf.gz"
  bcftools index -f --tbi --threads "$n_cores" "$TMP/in.vcf.gz" 2>/dev/null || \
    bcftools index -f --csi --threads "$n_cores" "$TMP/in.vcf.gz"
  n_sub=$(nvar "$TMP/in.vcf.gz")
  say "variants after dropping off-reference contigs : $n_sub (was $n_in)"
  bcftools view -H "$TMP/in.vcf.gz" | cut -f1 | uniq | sort -u > "$TMP/vcf.names"
fi
say "all $(wc -l < "$TMP/vcf.names") data-bearing contigs present in the reference"

# ── 2. ASSEMBLY GATE ─────────────────────────────────────────────────────
# The original script has no assembly verification at all, and its
# `--check-ref ws` would silently rewrite every mismatching REF -- turning a
# wrong-assembly VCF into a plausible-looking one. Two independent checks here,
# both fail-closed.
say "--- assembly verification (GRCh38) ---"

# 2a. Declared contig lengths vs the reference .fai, where the header has them.
awk '{print $1"\t"$2}' "${reference}.fai" | sed 's/^chr//' | sort -k1,1 > "$TMP/ref.len"
bcftools view -h "$TMP/in.vcf.gz" \
  | sed -n 's/^##contig=<ID=\([^,>]*\).*length=\([0-9]*\).*/\1\t\2/p' \
  | sed 's/^chr//' | sort -k1,1 > "$TMP/vcf.len" || true

if [[ -s "$TMP/vcf.len" ]]; then
  join -j1 "$TMP/vcf.len" "$TMP/ref.len" > "$TMP/len.join" || true
  n_cmp=$(wc -l < "$TMP/len.join"); n_bad=$(awk '$2!=$3' "$TMP/len.join" | wc -l)
  say "contig lengths compared: $n_cmp, mismatching: $n_bad"
  if [[ "$n_bad" -gt 0 ]]; then
    awk '$2!=$3 {printf "    %s: VCF=%s  REF=%s\n", $1,$2,$3}' "$TMP/len.join" | head -8 >&2
    echo "ERROR: contig lengths disagree with the reference. This VCF is not aligned" >&2
    echo "       to '$reference'. hg19 chr1=249250621 vs hg38 chr1=248956422 -- if you" >&2
    echo "       see that pair, the input is GRCh37 and must be lifted over, not patched." >&2
    exit 1
  fi
else
  say "WARNING: header declares no contig lengths; relying on the REF-match check below"
fi

# 2b. REF-allele agreement on a sample of records. A wrong-assembly VCF shows a
#     mismatch rate in the tens of percent; correct data is well under 1%.
bcftools view "$TMP/in.vcf.gz" | head -n 40000 > "$TMP/probe.vcf" || true
set +e
bcftools norm --check-ref w --fasta-ref "$reference" -O u "$TMP/probe.vcf" \
  > /dev/null 2> "$TMP/probe.err"
probe_rc=$?
set -e
# A non-zero exit, or a faidx "not found", means norm cannot process this VCF at
# all -- which is NOT the same as a REF-base mismatch and must not be counted as
# a clean probe. (Found by testing: a chr-prefixed VCF against an unprefixed
# reference reported 0% mismatch and then crashed inside the real norm call.)
if [[ "$probe_rc" -ne 0 ]] || grep -qiE 'was not found|faidx_fetch_seq failed|fai_build' "$TMP/probe.err"; then
  echo "ERROR: bcftools norm cannot process this VCF against '$reference' (rc=$probe_rc):" >&2
  head -6 "$TMP/probe.err" | sed 's/^/    /' >&2
  exit 1
fi
n_probe=$(grep -vc '^#' "$TMP/probe.vcf" || echo 0)
n_mis=$(grep -c 'REF_MISMATCH\|Reference allele mismatch' "$TMP/probe.err" || true)
n_mis=${n_mis:-0}
if [[ "$n_probe" -gt 0 ]]; then
  pct=$(( n_mis * 100 / n_probe ))
  say "REF mismatches on probe: $n_mis / $n_probe (${pct}%)  threshold ${max_ref_mismatch_pct}%"
  if [[ "$pct" -gt "$max_ref_mismatch_pct" ]]; then
    head -5 "$TMP/probe.err" >&2
    echo "ERROR: REF-allele mismatch rate ${pct}% exceeds ${max_ref_mismatch_pct}%." >&2
    echo "       Refusing to continue: this is the signature of a wrong assembly or a" >&2
    echo "       bad liftover. Raise with -m only if you know why the REFs differ." >&2
    exit 1
  fi
fi

# ── 3. Normalise: left-align, split multiallelics, atomize ───────────────
# NOTE: no -g/--gff-annot. See the header comment -- it right-aligns per the
# HGVS 3' rule, which breaks Exomiser's left-aligned exact-match lookup.
say "--- normalising (left-align, split multiallelics, atomize) ---"
bcftools norm \
  --fasta-ref "$reference" \
  --check-ref w \
  --multiallelics -any \
  --atomize \
  --rm-dup exact \
  --old-rec-tag EXO_ORIG \
  --threads "$n_cores" \
  -O z -o "$TMP/norm.vcf.gz" "$TMP/in.vcf.gz"
bcftools index -f --tbi --threads "$n_cores" "$TMP/norm.vcf.gz" 2>/dev/null || \
  bcftools index -f --csi --threads "$n_cores" "$TMP/norm.vcf.gz"
n_norm=$(nvar "$TMP/norm.vcf.gz")
say "variants after norm : $n_norm  (was $n_in; multiallelic splits raise this)"

# ── 4. Drop '*' spanning deletions ───────────────────────────────────────
say "--- removing '*' spanning-deletion alleles ---"
bcftools view -e 'ALT="*"' --threads "$n_cores" -O z -o "$TMP/nostar.vcf.gz" "$TMP/norm.vcf.gz"
bcftools index -f --tbi --threads "$n_cores" "$TMP/nostar.vcf.gz" 2>/dev/null || true
n_star=$(nvar "$TMP/nostar.vcf.gz")
say "variants after '*' removal : $n_star  (dropped $(( n_norm - n_star )))"
cur="$TMP/nostar.vcf.gz"

# ── 5. Optional quality prefilter ────────────────────────────────────────
if [[ "$qual_filter" -eq 1 ]]; then
  say "--- quality prefilter (VAF 15-85% for hets, GQ >= 20) ---"
  hdr=$(bcftools view -h "$cur")
  have_ad=0; have_gq=0; have_dp=0
  grep -q 'ID=AD,' <<<"$hdr" && have_ad=1
  grep -q 'ID=GQ,' <<<"$hdr" && have_gq=1
  grep -q 'ID=DP,' <<<"$hdr" && have_dp=1
  expr_parts=()
  if [[ "$have_ad" -eq 1 && "$have_dp" -eq 1 ]]; then
    expr_parts+=('(GT="het" && FMT/DP[0]>0 && (FMT/AD[0:1]/FMT/DP[0]<0.15 || FMT/AD[0:1]/FMT/DP[0]>0.85))')
  else
    say "WARNING: AD and/or DP absent -- skipping the VAF component"
  fi
  if [[ "$have_gq" -eq 1 ]]; then
    expr_parts+=('(FMT/GQ[0]<20)')
  else
    say "WARNING: GQ absent -- skipping the GQ component"
  fi
  if [[ "${#expr_parts[@]}" -gt 0 ]]; then
    ex=$(IFS='|'; printf '%s' "${expr_parts[*]}")
    ex=${ex//|/ || }
    bcftools view -e "$ex" --threads "$n_cores" -O z -o "$TMP/qual.vcf.gz" "$cur"
    bcftools index -f --tbi --threads "$n_cores" "$TMP/qual.vcf.gz" 2>/dev/null || true
    n_q=$(nvar "$TMP/qual.vcf.gz")
    say "variants after quality filter : $n_q  (dropped $(( n_star - n_q )))"
    cur="$TMP/qual.vcf.gz"
  else
    say "quality prefilter requested but no usable tags; skipped"
  fi
fi

# ── 6. Sort, write, index ────────────────────────────────────────────────
# bcftools sort, not `sort -V -k1,1 -k2,2`. The original's key spec sorts POS
# as a version string rather than numerically, and appending a separately
# sorted body to the header can yield an order the index does not agree with.
say "--- sorting and indexing ---"
bcftools sort -T "$TMP/sort" -O z -o "$out_vcf" "$cur"
bcftools index -f --tbi --threads "$n_cores" "$out_vcf" 2>/dev/null || \
  bcftools index -f --csi --threads "$n_cores" "$out_vcf"

n_out=$(nvar "$out_vcf")
say "--- done ---"
say "variants out : $n_out"
say "output       : $out_vcf"
say "index        : $(ls "${out_vcf}".{tbi,csi} 2>/dev/null | tr '\n' ' ')"
printf '\nSUMMARY  in=%s  normalised=%s  after_star=%s  out=%s  sample=%s\n' \
  "$n_in" "$n_norm" "$n_star" "$n_out" "${samples[0]}"
echo "Next: bin/run-exomiser.sh ${samples[0]} $out_vcf <phenopacket.yml> <outdir>"
