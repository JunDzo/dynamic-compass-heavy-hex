#!/bin/bash
#SBATCH -J HeavyHexGenCsv
#SBATCH -t 1-00:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --array=0-69%30
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/slurm_common.sh"

hh_load_python_module

CODE_DIR="${CODE_DIR:-$(hh_project_root)}"
VENV_DIR="${VENV_DIR:-$CODE_DIR/.venv-slurm}"
hh_activate_venv_if_present "$VENV_DIR"
cd "$CODE_DIR"

CSV="${CSV:-data/analysis/swp_sug/sweep_range.csv}"
OUT_DIR="${OUT_DIR:-$CODE_DIR/out/slurm/no_resets}"
ROUNDS="${ROUNDS:-auto}"
BACKEND="${BACKEND:-modified_no_hook.unflagged.no_reset}"
# BACKEND="modified_no_hook.unflagged.reset"
TABLE="${TABLE:-$CODE_DIR/data/analysis/tbls/breaker_table_d35.txt}"
STYLE="${STYLE:-heavy_hex}"

mkdir -p "$CODE_DIR/logs/slurm/gen_csv" "$OUT_DIR"

NROWS=$(python - "$CSV" <<'PY'
import sys
import pandas as pd

df = pd.read_csv(sys.argv[1])
print(len(df))
PY
)

_get_col() {
python - "$CSV" "$1" "$2" <<'PY'
import math
import sys
import pandas as pd

csv, idx, col = sys.argv[1], int(sys.argv[2]), sys.argv[3]
df = pd.read_csv(csv)
if not (0 <= idx < len(df)):
    raise SystemExit(f"Row index {idx} out of range for {len(df)} rows")
val = df.iloc[idx].get(col, "")
if isinstance(val, float) and math.isnan(val):
    val = ""
print(val)
PY
}

TASKS=${SLURM_ARRAY_TASK_COUNT:-1}
IDX=${SLURM_ARRAY_TASK_ID:-0}
ROWS_PER_TASK=$(( (NROWS + TASKS - 1) / TASKS ))
START_IDX=$(( IDX * ROWS_PER_TASK ))
END_IDX=$(( START_IDX + ROWS_PER_TASK ))
if (( END_IDX > NROWS )); then END_IDX=$NROWS; fi
if (( START_IDX >= NROWS )); then
  echo "Task $IDX has no rows (START_IDX=$START_IDX >= NROWS=$NROWS). Exiting."
  exit 0
fi

echo "Task $IDX handling rows [$START_IDX, $END_IDX) of $NROWS; TASKS=$TASKS; ROWS_PER_TASK=$ROWS_PER_TASK"

for (( ROW_INDEX=START_IDX; ROW_INDEX<END_IDX; ROW_INDEX++ )); do
  X=$(_get_col "$ROW_INDEX" x)
  Y=$(_get_col "$ROW_INDEX" y)
  BASIS=$(_get_col "$ROW_INDEX" basis)
  Z_VAL=$(_get_col "$ROW_INDEX" z)
  export X Y Z_VAL BASIS BACKEND TABLE ROUNDS STYLE

  echo "[row $ROW_INDEX/$((NROWS-1))] X=$X Y=$Y BASIS=$BASIS"
  [ -n "$Z_VAL" ] && echo "[row $ROW_INDEX] Z list: $Z_VAL"

  bash scripts/pipeline/generate_csv_sweep.sh "$OUT_DIR"
  echo "Step 1: Circuits generated for row $ROW_INDEX"
done
