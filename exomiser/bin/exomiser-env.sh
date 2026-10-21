#!/usr/bin/env bash
# Source this, do not execute it. Deliberately sets no shell options: it is
# sourced into interactive shells too, and `set -u` there is hostile.

export EXO_ROOT=/fast/AG_Bunina/Berk/tools/exomiser
export EXO_VERSION=15.1.0
export EXO_JAR="$EXO_ROOT/exomiser-cli-$EXO_VERSION/exomiser-cli-$EXO_VERSION.jar"
export EXO_CONF="$EXO_ROOT/etc/application.properties"
export EXO_EXAMPLES="$EXO_ROOT/exomiser-cli-$EXO_VERSION/examples"
export EXO_ENV=exomiser
export EXO_DATA_VER=2602

# ---------------------------------------------------------------------------
# This login environment ships SSL_CERT_FILE=/etc/ssl/certs/ca-bundle.crt,
# which does NOT exist on this host. Every curl/wget HTTPS call then dies with
# "curl: (77) error setting certificate file". Point it at the real bundle.
# ---------------------------------------------------------------------------
if [[ ! -r "${SSL_CERT_FILE:-/nonexistent}" ]]; then
  export SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt
fi
export CURL_CA_BUNDLE="$SSL_CERT_FILE"
export REQUESTS_CA_BUNDLE="$SSL_CERT_FILE"

# htslib tooling (bcftools/samtools/tabix/bgzip 1.24) is installed INSIDE the
# micromamba env from bioconda, deliberately. The group build in
# /fast/AG_Bunina/Berk/bin links against libdeflate.so.0 from
# /lib/x86_64-linux-gnu, which exists inside a Singularity container but NOT on
# the bare host -- there it dies with
#   bcftools: error while loading shared libraries: libdeflate.so.0
# So the group bin is appended as a fallback only, never preferred.
export PATH="$PATH:/fast/AG_Bunina/Berk/bin"

# ---------------------------------------------------------------------------
# Activate the micromamba env so `java` 21 is on PATH. Done here, not in the
# caller, so every run is self-contained. Use the binary, not the shell
# function: the function is only defined in interactive shells.
# ---------------------------------------------------------------------------
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-/fast/AG_Bunina/Berk/envs/mamba}"
_exo_mm="$(command -v micromamba 2>/dev/null)"
[[ -x "$_exo_mm" ]] || _exo_mm=/home/sgoenen/bin/micromamba
if [[ "${CONDA_DEFAULT_ENV:-}" != "$EXO_ENV" ]]; then
  if [[ -x "$_exo_mm" ]]; then
    eval "$("$_exo_mm" shell hook --shell bash)" 2>/dev/null
    micromamba activate "$EXO_ENV" 2>/dev/null
  fi
fi
# Belt and braces: if activation did not take, use the env prefix directly.
if ! command -v java >/dev/null 2>&1; then
  export PATH="$MAMBA_ROOT_PREFIX/envs/$EXO_ENV/bin:$PATH"
fi

# Make the env's own shared libraries resolvable, so anything linked against
# them works regardless of host vs container.
if [[ -d "$MAMBA_ROOT_PREFIX/envs/$EXO_ENV/lib" ]]; then
  export LD_LIBRARY_PATH="$MAMBA_ROOT_PREFIX/envs/$EXO_ENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi

# Fail loudly and early rather than letting Exomiser die on a class-file error.
_exo_java_major=$(java -version 2>&1 | sed -n '1s/[^"]*"\([0-9]*\).*/\1/p')
if [[ -z "$_exo_java_major" ]] || (( _exo_java_major < 21 )); then
  echo "ERROR: Exomiser $EXO_VERSION needs Java >= 21, found '${_exo_java_major:-none}'" >&2
  return 1 2>/dev/null || exit 1
fi

# Issue #541: Spring Boot logs to /tmp/spring.log by default. On a shared
# cluster the first user creates it 0644 and every other user's job dies with
# "openFile(/tmp/spring.log,true) call failed ... (Permission denied)".
# Give each run its own tmpdir and its own log.
exo_jvm_flags() {
  local tag="${1:-run}"
  local tmp="${TMPDIR:-/tmp}/exomiser.$USER.$$"
  mkdir -p "$tmp" "$EXO_ROOT/logs"
  printf '%s\n' \
    "-Djava.io.tmpdir=$tmp" \
    "-Dlogging.file.name=$EXO_ROOT/logs/exomiser.$tag.$$.log" \
    "-Dspring.config.location=$EXO_CONF"
}

# Verify the htslib tools actually EXECUTE, not merely exist on PATH. A binary
# whose shared libraries cannot be resolved passes `command -v` and then fails
# at the first real call -- which is exactly how the libdeflate.so.0 breakage
# surfaced, several steps into a patient run.
exo_check_tools() {
  local bad=0 t
  for t in bcftools samtools tabix bgzip; do
    if ! command -v "$t" >/dev/null 2>&1; then
      echo "ERROR: $t not on PATH" >&2; bad=1
    elif ! "$t" --version >/dev/null 2>&1; then
      echo "ERROR: $t is on PATH at $(command -v "$t") but will not run:" >&2
      "$t" --version 2>&1 | head -2 | sed 's/^/       /' >&2
      bad=1
    fi
  done
  return $bad
}

# ─────────────────────────────────────────────────────────────────────────────
# Effective memory limit for this process, in GiB.
# Walks the cgroup v2 tree and takes the tightest memory.max, because on an
# HPC the binding limit is the SLURM allocation, not the host's RAM. Falls
# back to MemAvailable when unconstrained.
# ─────────────────────────────────────────────────────────────────────────────
exo_mem_limit_gib() {
  local cg path p v best=""
  cg=$(awk -F: '{print $3; exit}' /proc/self/cgroup 2>/dev/null)
  path=/sys/fs/cgroup
  local IFS='/'
  for p in $cg; do
    [[ -z "$p" ]] && continue
    path="$path/$p"
    if [[ -r "$path/memory.max" ]]; then
      v=$(cat "$path/memory.max" 2>/dev/null)
      if [[ "$v" != "max" && -n "$v" ]]; then
        if [[ -z "$best" || "$v" -lt "$best" ]]; then best="$v"; fi
      fi
    fi
  done
  unset IFS
  # cgroup v1 fallback
  if [[ -z "$best" && -r /sys/fs/cgroup/memory/memory.limit_in_bytes ]]; then
    v=$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes 2>/dev/null)
    [[ -n "$v" && "$v" -lt 100000000000000 ]] && best="$v"
  fi
  if [[ -n "$best" ]]; then
    echo $(( best / 1024 / 1024 / 1024 ))
  else
    awk '/MemAvailable/{printf "%d", $2/1024/1024}' /proc/meminfo 2>/dev/null || echo 8
  fi
}

# Heap size for a pass, sized to the allocation instead of hardcoded.
#   exo_heap_gib exome|genome
# Reserves headroom for JVM non-heap plus the ~54 GiB of memory-mapped MVStore
# files Exomiser opens, which are charged to the same cgroup.
exo_heap_gib() {
  local mode="${1:-exome}" lim want reserve floor ceil
  lim=$(exo_mem_limit_gib)
  case "$mode" in
    genome) reserve=6; floor=6;  ceil=16 ;;
    *)      reserve=4; floor=3;  ceil=8  ;;
  esac
  want=$(( lim - reserve ))
  (( want > ceil )) && want=$ceil
  (( want < floor )) && want=$floor
  if (( lim < 12 && mode == "genome" )); then :; fi
  echo "$want"
}

# One-line advisory about whether this allocation can take a genome run.
exo_mem_report() {
  local lim; lim=$(exo_mem_limit_gib)
  echo "memory limit detected: ${lim} GiB  (exome heap $(exo_heap_gib exome)g, genome heap $(exo_heap_gib genome)g)"
  if (( lim < 16 )); then
    echo "WARNING: Exomiser documents 8 GiB minimum for a ~4.4M-variant genome," >&2
    echo "         and recommends 16 GiB total. With only ${lim} GiB the genome pass" >&2
    echo "         is likely to be killed (SIGKILL / exit 137). Submit it to SLURM" >&2
    echo "         with --mem=32G or more: see bin/exomiser.sbatch" >&2
  fi
}
