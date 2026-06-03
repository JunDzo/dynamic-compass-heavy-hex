#!/usr/bin/env bash
set -euo pipefail

CT_DIR=data/analysis/comm_th/com_th_nr_2.csv

python tools/threshold/render_rotation.py --csv "$CT_DIR" --basis X --fit quadratic \
    --out threshold_rotate_X.mp4 --azim-start 140 --azim-end 220 --frames 200 --fps 25 --elev 25
