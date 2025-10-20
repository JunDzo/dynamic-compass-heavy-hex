#!/bin/bash

set -e
out_dir=$1

# out_dir="out"
cir_dir="$out_dir/circuits"
    
PYTHONPATH=src parallel --ungroup tools/gen_circuits \
    --out_dir "$cir_dir" \
    --diameter {1} \
    --rounds "auto" \
    --noise_model uniform \
    --noise_strength {2} \
    --style {3} \
    --b {4} \
    --table_file "/bucket/ElkoussU/jun/HH-Modified/plot/tb1-th/table.txt" \
    ::: {3..25} \
    ::: 1e-3 \
    ::: heavy_hex \
    ::: X Z \

# --table_file "$out_dir/table.txt" \

# PYTHONPATH=src parallel --ungroup tools/gen_circuits \
#      --out_dir "$cir_dir"  \
#      --diameter {1} \
#      --rounds "auto" \
#      --noise_model uniform \
#      --noise_strength {2} \
#      --style heavy_hex \
#      --b {3} \
#      --table_file "/bucket/ElkoussU/jun/HH-Modified/table.txt" \
#      ::: {3..13} \
#      ::: 1e-6 2e-6 3e-6 5e-6 7e-6 1e-5 2e-5 3e-5 5e-5 7e-5 1e-4 2e-4 3e-4 5e-4 7e-4 1e-3 2e-3 3e-3 5e-3 7e-3 1e-2 2e-2 3e-2 5e-2 7e-2 1e-1 \
#      ::: X Z \