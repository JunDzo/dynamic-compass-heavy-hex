#!/usr/bin/env bash
set -euo pipefail

ST_DIR_NR=data/stats/stats_nr.csv
CT_DIR_NR=data/analysis/comm_th/com_th_nr.csv

python tools/threshold/find_approx_thresholds.py "$ST_DIR_NR" "$CT_DIR_NR" "idle" "linear"
python tools/threshold/count_thresholds.py --csv "$CT_DIR_NR"

################################################################################

ST_DIR_R=data/stats/stats_r.csv
CT_DIR_R=data/analysis/comm_th/com_th_r.csv

python tools/threshold/find_approx_thresholds.py "$ST_DIR_R" "$CT_DIR_R" "idle" "linear"
python tools/threshold/count_thresholds.py --csv "$CT_DIR_R"
