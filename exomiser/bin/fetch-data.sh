#!/usr/bin/env bash
# Downloads Exomiser data 2602 (hg38 + phenotype) and ReMM v0.4 hg38.
# Resumable: safe to re-run after an interruption.
set -uo pipefail
source /fast/AG_Bunina/Berk/tools/exomiser/bin/exomiser-env.sh || exit 1

DATA_URL=https://g-879a9f.f5dc97.75bc.dn.glob.us/data
cd "$EXO_ROOT/data" || exit 1

for f in "${EXO_DATA_VER}_phenotype" "${EXO_DATA_VER}_hg38"; do
  echo "### $(date '+%F %T')  fetching $f.zip"
  wget -c -q --show-progress --progress=dot:giga "$DATA_URL/$f.zip" 2>&1 | tail -2
  echo "### $(date '+%F %T')  done $f.zip -> $(stat -c '%s' "$f.zip" 2>/dev/null) bytes"
done

echo "### $(date '+%F %T')  fetching ReMM v0.4 hg38"
cd "$EXO_ROOT/data/remm" || exit 1
wget -c -q --show-progress --progress=dot:giga -O ReMM.v0.4.hg38.tsv.gz \
  "https://zenodo.org/records/6576087/files/ReMM.v0.4.hg38.tsv.gz?download=1" 2>&1 | tail -2
wget -c -q -O ReMM.v0.4.hg38.tsv.gz.tbi \
  "https://zenodo.org/records/6576087/files/ReMM.v0.4.hg38.tsv.gz.tbi?download=1"
echo "### $(date '+%F %T')  ALL DOWNLOADS COMPLETE"
ls -l "$EXO_ROOT/data" "$EXO_ROOT/data/remm"
