#!/bin/bash
INX_DIR=data/analysis/swp_sug/sweep_suggestions_linear_axis_x.csv
INZ_DIR=data/analysis/swp_sug/sweep_suggestions_linear_axis_z.csv

OUTX_DIR=data/analysis/swp_sug/sweep_range_x.csv
OUTZ_DIR=data/analysis/swp_sug/sweep_range_z.csv

OUT_DIR=data/analysis/swp_sug/sweep_range.csv

python tools/expand_sweep_points.py \
  --in $INX_DIR \
  --out $OUTZ_DIR \
  --eta 6e-4 \
  --n 18 \
  --clip-low 0 \
  --diameters "7 9 11 13" \
  --base-prob 2e-4 

python tools/expand_sweep_points.py \
  --in $INX_DIR \
  --out $OUTX_DIR \
  --eta 6e-4 \
  --n 18 \
  --clip-low 0 \
  --diameters "7 9 11 13" \
  --base-prob 2e-4 

cat $OUTZ_DIR $OUTX_DIR | awk '!seen[$0]++' > $OUT_DIR

rm -rf $INX_DIR $INZ_DIR $OUTZ_DIR $OUTX_DIR
