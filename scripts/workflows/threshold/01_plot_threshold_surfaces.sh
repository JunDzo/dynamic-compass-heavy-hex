#!/usr/bin/env bash
set -euo pipefail


# OUTX_DIR=data/analysis/swp_sug/sweep_suggestions_linear_axis_x.csv
# OUTZ_DIR=data/analysis/swp_sug/sweep_suggestions_linear_axis_z.csv

#------------------------------------- WITH NO-RESET -------------------------------------#
CT_DIR_NR=data/analysis/comm_th/com_th_nr.csv
OUT_PATH_X_NR=data/plots/surfaces/noreset/no_reset_surface_X_e14_a45.pdf
OUT_PATH_Z_NR=data/plots/surfaces/noreset/no_reset_surface_Z_e14_a45.pdf


python tools/threshold/plot_threshold_surface_3d.py \
  --csv "$CT_DIR_NR" \
  --basis X \
  --x-axis "b1p" \
  --y-axis "b2p" \
  --z-axis "p_th" \
  --fit-model "quadratic"\
  --fit-out "data/analysis/fit-surf/fit_nr_x.json" \
  --out "$OUT_PATH_X_NR" \
  #--mark-point 0.0041,0.054,0.012 \
  #--emit-points 200 \
  #--emit-points-mode axis \
  #--plot-emitted \
  #--emit-points-out $OUTX_DIR \



python tools/threshold/plot_threshold_surface_3d.py \
  --csv "$CT_DIR_NR" \
  --basis Z \
  --x-axis "b1p" \
  --y-axis "b2p" \
  --z-axis "p_th" \
  --fit-model "quadratic"\
  --ylim "0,0.118"\
  --fit-out "data/analysis/fit-surf/fit_nr_z.json" \
  --out "$OUT_PATH_Z_NR"
  #--mark-point 0.0041,0.054,0.012 \
  #--emit-points 200 \
  #--emit-points-mode axis \
  #--plot-emitted \
  #--emit-points-out $OUTZ_DIR \

#------------------------------------- WITH RESET -------------------------------------#
CT_DIR_R=data/analysis/comm_th/com_th_r.csv
OUT_PATH_X_R=data/plots/surfaces/reset/reset_surface_X_e14_a45.pdf
OUT_PATH_Z_R=data/plots/surfaces/reset/reset_surface_Z_e14_a45.pdf

python tools/threshold/plot_threshold_surface_3d.py \
  --csv "$CT_DIR_R" \
  --basis X \
  --x-axis "b1p" \
  --y-axis "b2p" \
  --z-axis "p_th" \
  --fit-model "quadratic"\
  --fit-out "data/analysis/fit-surf/fit_r_x.json" \
  --out "$OUT_PATH_X_R" \
  #--mark-point 0.0041,0.054,0.012 \
  #--emit-points 200 \
  #--emit-points-mode axis \
  #--plot-emitted \
  #--emit-points-out $OUTX_DIR \



python tools/threshold/plot_threshold_surface_3d.py \
  --csv "$CT_DIR_R" \
  --basis Z \
  --x-axis "b1p" \
  --y-axis "b2p" \
  --z-axis "p_th" \
  --fit-model "quadratic"\
  --ylim "0,0.118"\
  --fit-out "data/analysis/fit-surf/fit_r_z.json" \
  --out "$OUT_PATH_Z_R"
  #--mark-point 0.0041,0.054,0.012 \
  #--emit-points 200 \
  #--emit-points-mode axis \
  #--plot-emitted \
  #--emit-points-out $OUTZ_DIR \