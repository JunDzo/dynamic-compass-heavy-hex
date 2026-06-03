#!/usr/bin/env bash
set -euo pipefail

out_dir=${1:?usage: scripts/pipeline/collect_stats.sh OUT_DIR [STATS_FILE_OR_DIR]}
cir_dir="$out_dir/circuits"
stats_target="${2:-$out_dir/stats.csv}"
if [[ "$stats_target" == *.csv ]]; then
  stats_file="$stats_target"
else
  mkdir -p "$stats_target"
  stats_file="$stats_target/stats.csv"
fi
procs="${SLURM_CPUS_PER_TASK:-12}"

mkdir -p "$(dirname "$stats_file")"

sinter collect \
    --circuits "$cir_dir"/*.stim \
    --save_resume_filepath "$stats_file" \
    --metadata_func auto \
    --decoders pymatching \
    --max_shots 10_000_000 \
    --max_errors 1000 \
    --processes "$procs"
