#!/usr/bin/env bash
# Submit one Snakemake job to LSF in the form the Sanger farm expects, and print its LSF job ID.
# Called by Snakemake (see config.yaml): lsf_submit.sh <queue> <group> <mem_mb> <runtime_min> <threads> <rule> <jobscript>
set -euo pipefail
queue=$1; group=$2; mem=$3; runtime=$4; threads=$5; rule=$6; jobscript=$7
mem=${mem%.*}   # whole MB
# Make the farm's `module` command available to every shell in the job (needed for the modules in the config)
export BASH_ENV="$(cd "$(dirname "$0")" && pwd)/modules_init.sh"
mkdir -p logs/lsf
out=$(bsub -q "$queue" -G "$group" -n "$threads" -W "$runtime" \
    -M "$mem" -R "select[mem>=${mem}] rusage[mem=${mem}] span[hosts=1]" \
    -J "smk_${rule}" -o "logs/lsf/${rule}.%J.log" -e "logs/lsf/${rule}.%J.log" \
    "$jobscript")
echo "$out" | grep -oE 'Job <[0-9]+>' | grep -oE '[0-9]+'
