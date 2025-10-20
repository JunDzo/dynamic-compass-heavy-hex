#!/bin/bash

CT_DIR=data/analysis/comm_th/common_thresholds_3.csv
OUTX_DIR=data/analysis/swp_sug/sweep_suggestions_linear_axis_x.csv
OUTZ_DIR=data/analysis/swp_sug/sweep_suggestions_linear_axis_z.csv

python tools/plot_threshold_surface_3d.py \
  --csv $CT_DIR \
  --basis Z \
  --fit linear \
  --emit-points 200 \
  --emit-points-mode axis \
  --plot-emitted \
  --emit-points-out $OUTZ_DIR \

