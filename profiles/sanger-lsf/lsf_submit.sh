#!/usr/bin/env bash
# Submit one Snakemake job to LSF in the form the Sanger farm expects, and print its LSF job ID.
# Called by Snakemake (see config.yaml): lsf_submit.sh <queue> <group> <mem_mb> <runtime_min> <threads> <rule> <jobscript>
set -euo pipefail
queue=$1; group=$2; mem=$3; runtime=$4; threads=$5; rule=$6; jobscript=$7
mem=${mem%.*}   # whole MB
# A retry asking for more time than its queue allows moves up a queue (normal 12 h, long 48 h, basement 30 days)
runtime=${runtime%.*}
if [ "$queue" = normal ] && [ "$runtime" -gt 720 ]; then queue=long; fi
if [ "$queue" = long ] && [ "$runtime" -gt 2880 ]; then queue=basement; fi
# Make the farm's `module` command available to every shell in the job (needed for the modules in the config)
export BASH_ENV="$(cd "$(dirname "$0")" && pwd)/modules_init.sh"
# LSF logs go to <study_dir>/logs/lsf (set by the Snakefile)
logs=${PHYLOBUILD_LSF_LOGS:-logs/lsf}
mkdir -p "$logs"
out=$(bsub -q "$queue" -G "$group" -n "$threads" -W "$runtime" \
    -M "$mem" -R "select[mem>=${mem}] rusage[mem=${mem}] span[hosts=1]" \
    -J "smk_${rule}" -o "$logs/${rule}.%J.log" -e "$logs/${rule}.%J.log" \
    "$jobscript")
echo "$out" | grep -oE 'Job <[0-9]+>' | grep -oE '[0-9]+'
