#!/bin/bash
#SBATCH -J HeavyHexGen
#SBATCH -t 1-00:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --array=0-9
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

OUT_DIR="${OUT_DIR:-$CODE_DIR/out/slurm/no_reset}"
TABLE="${TABLE:-$CODE_DIR/data/analysis/tbls/breaker_table_d35.txt}"
ROUNDS="${ROUNDS:-auto}"
BACKEND="${BACKEND:-modified_no_hook.unflagged.no_reset}"
# BACKEND="modified_no_hook.unflagged.reset"
STYLE="${STYLE:-heavy_hex}"
NOISE_MODEL="${NOISE_MODEL:-uniform}"

mkdir -p "$CODE_DIR/logs/slurm/gen" "$OUT_DIR"

VAL_LIST=(
  1.00e-4 1.20e-4 1.43e-4 1.71e-4 2.04e-4
  2.44e-4 2.91e-4 3.48e-4 4.16e-4 4.97e-4
)

TOTAL=${#VAL_LIST[@]}
N_ARRAY=${SLURM_ARRAY_TASK_COUNT:-1}
TASK=${SLURM_ARRAY_TASK_ID:-0}
CHUNK_SIZE=$(( (TOTAL + N_ARRAY - 1) / N_ARRAY ))
START=$(( TASK * CHUNK_SIZE ))
SUBARRAY=( "${VAL_LIST[@]:START:CHUNK_SIZE}" )

if (( ${#SUBARRAY[@]} == 0 )); then
  echo "Task $TASK has no VAL entries. Exiting."
  exit 0
fi

for VAL in "${SUBARRAY[@]}"; do
  export VAL TABLE ROUNDS BACKEND STYLE NOISE_MODEL
  echo "TASK=$TASK VAL=$VAL"
  bash scripts/pipeline/generate_uniform.sh "$OUT_DIR"
  echo "Step 1: Circuits generated for VAL=$VAL"
done
