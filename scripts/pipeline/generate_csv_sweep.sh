#!/usr/bin/env bash
set -euo pipefail
out_dir=${1:?usage: scripts/pipeline/generate_csv_sweep.sh OUT_DIR}

# out_dir="out-test"
cir_dir="$out_dir/circuits"
mkdir -p "$cir_dir"

# Require Slurm-provided env values (array layer)
: "${Y:?set Y to the M/b2p probability}"
: "${X:?set X to the D2/b1p probability}"
: "${Z_VAL:?set Z_VAL to one or more idle sweep values}"
: "${BASIS:?set BASIS to X, Z, or both}"
: "${BACKEND:?set BACKEND to the circuit backend, e.g. modified_no_hook.unflagged.no_reset}"
: "${TABLE:?set TABLE to the breaker table file path}"

diameters="${DIAMETERS:-9 11 13 15 17}"
base_prob="${BASE_PROB:-2e-4}"
style="${STYLE:-heavy_hex}"

PYTHONPATH=src parallel --ungroup tools/generate_circuit.py \
    --out_dir "$cir_dir" \
    --diameter {1} \
    --rounds "auto" \
    --noise_model corr_mixing \
    --noise_strength {2} \
    --b {3} \
    --backend {4} \
    --table_file {5} \
    --base_prob {6} \
    --sweep_slot idle \
    --bias_slot_1 D2 \
    --bias_slot_1_prob "$X" \
    --bias_slot_2 M \
    --bias_slot_2_prob "$Y" \
    --style {7} \
    ::: $diameters \
    ::: $Z_VAL \
    ::: $BASIS \
    ::: $BACKEND \
    ::: $TABLE \
    ::: $base_prob \
    ::: $style

# PYTHONPATH=src parallel --ungroup tools/gen_circuits \
#     --out_dir "$cir_dir" \
#     --diameter {1} \
#     --rounds "auto" \
#     --noise_model corr_mixing \
#     --noise_strength {2} \
#     --style heavy_hex \
#     --b {3} \
#     --table_file "out/plot/t1-th/table.txt" \
#     --base_prob {4} \
#     --sweep_slot M \
#     --bias_slot_1 D1 \
#     --bias_slot_1_prob {5} \
#     --bias_slot_2 D2 \
#     --bias_slot_2_prob {6} \
#     ::: 13 \
#     ::: 1e-3 \
#     ::: X \
#     ::: 1.2e-2 \
#     ::: 1e-6 2e-6 \
#     ::: 1e-6 2e-6 \

# 1e-6 2e-6 3e-6 5e-6 7e-6 1e-5 2e-5 3e-5 5e-5 7e-5
# 1e-2 2e-2 3e-2 5e-2 7e-2 1e-1
