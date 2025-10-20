#!/bin/bash


out_dir=${1:?need out_dir}
circuits_dir="${out_dir}/circuits"
stats_file="${STATS_FILE:-${out_dir}/stats.csv}"   # <- allow override via env
procs="${SLURM_CPUS_PER_TASK:-12}"


python tesseractdecoder.py \
  --circuits "$circuits_dir"/*.stim \
  --save-resume "$stats_file" \
  --processes "$procs" \
  --max-shots 10_000_000 \
  --max-errors 1000

# python bpmatching.py \
#   --circuits "$circuits_dir"/*.stim \
#   --save-resume "$stats_file" \
#   --processes "$procs" \
#   --max_bp_iters 50 \
#   --max-shots 10_000_000 \
#   --max-errors 1000

# python mwpfdecoder.py \
#   --circuits "$circuits_dir"/*.stim \
#   --save-resume "$stats_file" \
#   --processes "$procs" \
#   --max-shots 10_000_000 \
#   --max-errors 1000