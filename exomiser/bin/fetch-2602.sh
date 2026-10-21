#!/usr/bin/env bash
# Fetch Exomiser data 2602 (hg38 + phenotype) and extract in place.
# Resumable: re-run to continue from completed chunks.
# 4-way parallelism, not 8 -- Zenodo/Globus dropped connections at 8.
set -uo pipefail
source /fast/AG_Bunina/Berk/tools/exomiser/bin/exomiser-env.sh || exit 1
U=https://g-879a9f.f5dc97.75bc.dn.glob.us/data
cd "$EXO_ROOT/data" || exit 1
for f in 2602_phenotype 2602_hg38; do
  echo "### $(date '+%F %T') fetching $f.zip"
  "$EXO_ROOT/bin/pget.sh" "$U/$f.zip" "$f.zip" 4 512 || { echo "### FAILED $f"; exit 1; }
  echo "### $(date '+%F %T') extracting $f.zip"
  unzip -q -o "$f.zip" -d "$EXO_ROOT/data" && rm -f "$f.zip"
  df -h "$EXO_ROOT/data" | tail -1
done
V=2602
for f in "${V}_hg38/${V}_hg38_variants.mv.db" "${V}_hg38/${V}_hg38_clinvar.mv.db" \
         "${V}_hg38/${V}_hg38_transcripts_ensembl.ser" "${V}_phenotype"; do
  if [ -e "$EXO_ROOT/data/$f" ]; then echo "OK      $f"; else echo "MISSING $f"; fi
done
echo "### $(date '+%F %T') 2602 READY - switch etc/application.properties 2512 -> 2602"
