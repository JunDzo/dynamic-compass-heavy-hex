#!/usr/bin/env bash
set -euo pipefail

INX_DIR=data/analysis/swp_sug/sweep_suggestions_linear_axis_x.csv
INZ_DIR=data/analysis/swp_sug/sweep_suggestions_linear_axis_z.csv

OUTX_DIR=data/analysis/swp_sug/sweep_range_x.csv
OUTZ_DIR=data/analysis/swp_sug/sweep_range_z.csv

OUT_DIR=data/analysis/swp_sug/sweep_range.csv

GENERATED_OUTS=()

if [ -f "$INZ_DIR" ]; then
  python tools/threshold/expand_sweep_points.py \
    --in "$INZ_DIR" \
    --out "$OUTZ_DIR" \
    --eta 1e-3 \
    --n 21 \
    --clip-low 0 \
    --diameters "9 11 13 15 17" \
    --base-prob 2e-4
  GENERATED_OUTS+=("$OUTZ_DIR")
fi

if [ -f "$INX_DIR" ]; then
  python tools/threshold/expand_sweep_points.py \
    --in "$INX_DIR" \
    --out "$OUTX_DIR" \
    --eta 3e-4 \
    --n 9 \
    --clip-low 0 \
    --diameters "9 11 13 15 17" \
    --base-prob 2e-4
  GENERATED_OUTS+=("$OUTX_DIR")
fi

if [ ${#GENERATED_OUTS[@]} -eq 0 ]; then
  echo "No input found. Expected at least one of: $INX_DIR or $INZ_DIR" >&2
  exit 1
fi

cat "${GENERATED_OUTS[@]}" | awk '!seen[$0]++' > "$OUT_DIR"

rm -f "$INX_DIR" "$INZ_DIR" "$OUTZ_DIR" "$OUTX_DIR"
