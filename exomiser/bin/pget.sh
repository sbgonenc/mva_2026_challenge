#!/usr/bin/env bash
# pget.sh <url> <outfile> [nconn] [chunk_MB]
# Chunked parallel HTTPS downloader. Resumable: completed chunks are skipped
# on re-run. Exists because some of these hosts cap a single connection well
# below the link rate but serve HTTP range requests happily.
set -uo pipefail

URL="${1:?url required}"; OUT="${2:?outfile required}"
NCONN="${3:-8}"; CHUNK=$(( ${4:-512} * 1024 * 1024 ))
: "${SSL_CERT_FILE:=/etc/ssl/certs/ca-certificates.crt}"
export SSL_CERT_FILE CURL_CA_BUNDLE="$SSL_CERT_FILE"

SIZE=$(curl -sIL --max-time 60 "$URL" \
       | awk 'tolower($1)=="content-length:"{c=$2} END{printf "%d", c+0}')
if (( SIZE <= 0 )); then echo "ERROR: could not determine size of $URL" >&2; exit 1; fi

if [[ -s "$OUT" ]] && (( $(stat -c '%s' "$OUT") == SIZE )); then
  echo "already complete: $OUT ($SIZE bytes)"; exit 0
fi

PARTS="$OUT.parts"; mkdir -p "$PARTS"
N=$(( (SIZE + CHUNK - 1) / CHUNK ))
printf 'downloading -> %s\n  %d bytes, %d chunks, %d parallel\n' "$OUT" "$SIZE" "$N" "$NCONN"

: > "$PARTS/.todo"
for (( i=0; i<N; i++ )); do
  s=$(( i * CHUNK )); e=$(( s + CHUNK - 1 ))
  (( e > SIZE-1 )) && e=$(( SIZE-1 ))
  exp=$(( e - s + 1 ))
  f="$PARTS/p$(printf '%05d' "$i")"
  if [[ -s "$f" ]] && (( $(stat -c '%s' "$f") == exp )); then continue; fi
  echo "$i $s $e" >> "$PARTS/.todo"
done

todo=$(wc -l < "$PARTS/.todo")
echo "  $todo chunks to fetch ($(( N - todo )) already present)"

if (( todo > 0 )); then
  xargs -a "$PARTS/.todo" -P "$NCONN" -n 3 bash -c '
    i=$0; s=$1; e=$2
    f="'"$PARTS"'/p$(printf "%05d" "$i")"
    for attempt in 1 2 3 4 5; do
      curl -sS --max-time 1800 --retry 2 --retry-delay 3 \
           -r "${s}-${e}" -o "$f" "'"$URL"'" && \
        [[ $(stat -c "%s" "$f" 2>/dev/null || echo 0) -eq $((e-s+1)) ]] && exit 0
      sleep $(( attempt * 4 ))
    done
    echo "CHUNK FAILED: $i ($s-$e)" >&2; exit 1
  '
fi

# Verify every chunk before assembling: a short chunk silently corrupts the file.
bad=0
for (( i=0; i<N; i++ )); do
  s=$(( i * CHUNK )); e=$(( s + CHUNK - 1 )); (( e > SIZE-1 )) && e=$(( SIZE-1 ))
  f="$PARTS/p$(printf '%05d' "$i")"
  got=$(stat -c '%s' "$f" 2>/dev/null || echo 0)
  if (( got != e - s + 1 )); then echo "INCOMPLETE chunk $i: $got != $((e-s+1))" >&2; bad=1; fi
done
(( bad )) && { echo "ERROR: incomplete chunks; re-run to resume" >&2; exit 1; }

cat "$PARTS"/p[0-9]* > "$OUT" || exit 1
final=$(stat -c '%s' "$OUT")
if (( final != SIZE )); then
  echo "ERROR: assembled size $final != expected $SIZE" >&2; exit 1
fi
rm -rf "$PARTS"
echo "OK: $OUT ($final bytes)"
