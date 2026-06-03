#!/bin/bash
#SBATCH --job-name=HeavyHexEnv
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/slurm_common.sh"

hh_load_python_module

CODE_DIR="${CODE_DIR:-$(hh_project_root)}"
VENV_DIR="${VENV_DIR:-$CODE_DIR/.venv-slurm}"
PYTHON_BIN="$(hh_python_bin)"

mkdir -p "$CODE_DIR/logs/slurm/env"
rm -rf "$VENV_DIR"
mkdir -p "$(dirname "$VENV_DIR")"
"$PYTHON_BIN" -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"

cd "$CODE_DIR"
python -m pip install --upgrade pip
python -m pip install -e .

echo "Done setting up the environment for Heavy Hex Modified."
