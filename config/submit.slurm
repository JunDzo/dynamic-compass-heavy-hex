#!/bin/bash
#SBATCH --job-name=HeavyHexStats
#SBATCH --time=4-00:00:00
#SBATCH --mem=128G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --array=0-31
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/slurm_common.sh"

hh_load_python_module

CODE_DIR="${CODE_DIR:-$(hh_project_root)}"
VENV_DIR="${VENV_DIR:-$CODE_DIR/.venv-slurm}"
hh_activate_venv_if_present "$VENV_DIR"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

OUT_DIR="${OUT_DIR:-$CODE_DIR/out/slurm/no_resets}"
CIR_DIR="$OUT_DIR/circuits"
PARTS_DIR="$OUT_DIR/parts"
MANIFEST="$OUT_DIR/circuits_manifest.txt"

mkdir -p "$CODE_DIR/logs/slurm/stats" "$PARTS_DIR"

if [[ ! -d "$CIR_DIR" ]]; then
  echo "No circuit directory found at $CIR_DIR."
  echo "Run config/generate.slurm or config/generate_csv.slurm first, or set OUT_DIR to an existing run."
  exit 1
fi

if [[ ! -s "$MANIFEST" ]]; then
  find "$CIR_DIR" -maxdepth 1 -type f -name '*.stim' | sort > "$MANIFEST"
fi

TOTAL=$(wc -l < "$MANIFEST")
TASKS=${SLURM_ARRAY_TASK_COUNT:-1}
IDX=${SLURM_ARRAY_TASK_ID:-0}

PER_TASK=$(( (TOTAL + TASKS - 1) / TASKS ))
START=$(( IDX * PER_TASK + 1 ))
END=$(( START + PER_TASK - 1 ))
if (( END > TOTAL )); then END=$TOTAL; fi

if (( START > END )); then
  echo "Array task $IDX: nothing to do."
  exit 0
fi

mapfile -t CHUNK < <(sed -n "${START},${END}p" "$MANIFEST")

STATS_SHARD="$PARTS_DIR/stats_${IDX}.csv"
echo "Task $IDX: processing ${#CHUNK[@]} circuits -> $STATS_SHARD"

TMP_OUT="$OUT_DIR/_tmp_${IDX}"
mkdir -p "$TMP_OUT/circuits"
for f in "${CHUNK[@]}"; do
  ln -sf "$f" "$TMP_OUT/circuits/"
done

STATS_FILE="$STATS_SHARD" OUT_DIR="$TMP_OUT" \
  bash "$CODE_DIR/scripts/pipeline/collect_stats.sh" "$TMP_OUT" "$STATS_SHARD"

rm -rf "$TMP_OUT"
