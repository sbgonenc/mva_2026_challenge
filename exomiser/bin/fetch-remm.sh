#!/usr/bin/env bash
# ReMM v0.4 hg38 -- the one piece genuinely missing for Genomiser.
set -uo pipefail
source /fast/AG_Bunina/Berk/tools/exomiser/bin/exomiser-env.sh || exit 1
cd "$EXO_ROOT/data/remm" || exit 1
Z=https://zenodo.org/records/6576087/files
echo "### $(date '+%F %T') ReMM hg38 tabix index"
"$EXO_ROOT/bin/pget.sh" "$Z/ReMM.v0.4.hg38.tsv.gz.tbi?download=1" ReMM.v0.4.hg38.tsv.gz.tbi 2 64
echo "### $(date '+%F %T') ReMM hg38 scores (15.97 GB)"
"$EXO_ROOT/bin/pget.sh" "$Z/ReMM.v0.4.hg38.tsv.gz?download=1" ReMM.v0.4.hg38.tsv.gz 8 512
echo "### $(date '+%F %T') verifying MD5 (must be 2f8c8fcc75ea29f6798621d0446fe194)"
md5sum ReMM.v0.4.hg38.tsv.gz ReMM.v0.4.hg38.tsv.gz.tbi
echo "### $(date '+%F %T') REMM FETCH COMPLETE"
