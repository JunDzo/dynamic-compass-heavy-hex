#!/bin/bash
ST_DIR=data/stats/stats_3.csv
CT_DIR=data/analysis/comm_th/common_thresholds_3.csv

python tools/find_common_thresholds.py $ST_DIR $CT_DIR
python tools/count_thresholds.py --csv $CT_DIR

