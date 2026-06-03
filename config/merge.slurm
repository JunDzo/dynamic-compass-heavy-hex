#!/bin/bash
#SBATCH --job-name=HeavyHexMerge
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/slurm_common.sh"

CODE_DIR="${CODE_DIR:-$(hh_project_root)}"
OUT_DIR="${OUT_DIR:-$CODE_DIR/out/slurm/no_resets}"
PARTS_DIR="${PARTS_DIR:-$OUT_DIR/parts}"
FINAL="${FINAL:-$OUT_DIR/stats.csv}"

mkdir -p "$(dirname "$FINAL")"

shopt -s nullglob
parts=( "$PARTS_DIR"/stats_*.csv )
if (( ${#parts[@]} == 0 )); then
  echo "No shard CSV files found in $PARTS_DIR."
  echo "Run config/submit.slurm first, or set OUT_DIR/PARTS_DIR to the directory containing stats_*.csv files."
  exit 1
fi

first=1
: > "$FINAL"
for f in "${parts[@]}"; do
  if [[ $first -eq 1 ]]; then
    cat "$f" > "$FINAL"
    first=0
  else
    tail -n +2 "$f" >> "$FINAL"
  fi
done

echo "Merged ${#parts[@]} shard CSV files -> $FINAL"
