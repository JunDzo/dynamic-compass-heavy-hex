#!/usr/bin/env bash
set -euo pipefail

out_dir=${1:?usage: scripts/pipeline/generate_uniform.sh OUT_DIR}

cir_dir="$out_dir/circuits"
mkdir -p "$cir_dir"

noise_strengths="${VAL:-${Z_VAL:-}}"
: "${noise_strengths:?set VAL or Z_VAL to one or more uniform noise strengths}"
: "${TABLE:?set TABLE to the breaker table file path}"
: "${ROUNDS:?set ROUNDS to auto or a space-separated rounds list}"
: "${BACKEND:?set BACKEND to the circuit backend, e.g. modified_no_hook.unflagged.no_reset}"

diameters="${DIAMETERS:-9 11 13 15 17}"
bases="${BASIS:-X}"
style="${STYLE:-heavy_hex}"
noise_model="${NOISE_MODEL:-uniform}"

PYTHONPATH=src parallel --ungroup tools/generate_circuit.py \
    --out_dir "$cir_dir" \
    --diameter {1} \
    --rounds {2} \
    --noise_model {3} \
    --noise_strength {4} \
    --b {5} \
    --backend {6} \
    --table_file {7} \
    --style {8} \
    ::: $diameters \
    ::: $ROUNDS \
    ::: $noise_model \
    ::: $noise_strengths \
    ::: $bases \
    ::: $BACKEND\
    ::: $TABLE \
    ::: $style

#  9 11 13 15 17
# PYTHONPATH=src parallel --ungroup tools/generate_circuit.py \
#     --out_dir "$cir_dir" \
#     --diameter {1} \
#     --rounds {2} \
#     --noise_model uniform \
#     --noise_strength {3} \
#     --b {4} \
#     --table_file {5} \
#     --backend {6} \
#     --style {7}} \
#     ::: $VAL \
#     ::: $ROUNDS \
#     ::: 1e-3 \
#     ::: Z X \
#     ::: $TABLE \
#     ::: $BACKEND\
#     ::: $STYLE\
