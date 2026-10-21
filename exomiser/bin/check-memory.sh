#!/usr/bin/env bash
# Report the memory limit that actually applies to THIS shell.
source /fast/AG_Bunina/Berk/tools/exomiser/bin/exomiser-env.sh >/dev/null 2>&1 || true
echo "host RAM        : $(free -g | awk '/^Mem:/{print $2" GiB total, "$7" GiB available"}')"
echo "SLURM job       : ${SLURM_JOB_ID:-<not in a SLURM job>}"
echo "SLURM mem grant : ${SLURM_MEM_PER_NODE:-<unset>} MB"
cg=$(awk -F: '{print $3; exit}' /proc/self/cgroup)
echo "cgroup          : $cg"
echo "effective limit : $(exo_mem_limit_gib) GiB"
d="/sys/fs/cgroup$cg"
[ -r "$d/memory.peak" ]   && echo "cgroup peak use : $(( $(cat $d/memory.peak)/1024/1024/1024 )) GiB"
[ -r "$d/memory.events" ] && { echo "OOM events:"; sed 's/^/  /' "$d/memory.events"; }
exo_mem_report
