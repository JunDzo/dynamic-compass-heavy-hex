#!/bin/bash
set -euo pipefail
# cd /Volumes/jun/HH-Modified

out_dir="out-test"
cir_dir="$out_dir/circuits"
mkdir -p "$cir_dir"

# rm -rf /Volumes/jun/HH-Modified/out-test/circuits && echo "complete: removed old circuits output directory"


# PYTHONPATH=src parallel --ungroup tools/generate_circuit.py \
#     --out_dir "$cir_dir" \
#     --diameter {1} \
#     --rounds "auto" \
#     --noise_model si1000 \
#     --noise_strength {2} \
#     --b {3} \
#     --backend modified_no_hook.unflagged.reset \
#     ::: 5 \
#     ::: 0.001 \
#     ::: X \
#     # --table_file "data/analysis/tbls/breaker_table_d35.txt" \

PYTHONPATH=src parallel --ungroup tools/generate_circuit.py \
    --out_dir "$cir_dir" \
    --diameter {1} \
    --rounds "auto" \
    --noise_model uniform \
    --noise_strength {2} \
    --b {3} \
    --backend modified_no_hook.unflagged.no_reset \
    --table_file {4} \
    ::: 5 \
    ::: 0.001 \
    ::: Z X \
    ::: "data/analysis/tbls/breaker_table_d35.txt" \