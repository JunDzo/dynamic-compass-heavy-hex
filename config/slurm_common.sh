#!/bin/bash

hh_project_root() {
  local script_dir
  script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  cd "$script_dir/.." && pwd
}

hh_load_python_module() {
  if ! command -v module >/dev/null 2>&1; then
    return 0
  fi

  module purge >/dev/null 2>&1 || true
  if [[ -n "${PYTHON_MODULE:-}" ]]; then
    module load "$PYTHON_MODULE"
  else
    module load python/3.11.4 >/dev/null 2>&1 || true
  fi
}

hh_python_bin() {
  if [[ -n "${PYTHON_BIN:-}" ]]; then
    printf '%s\n' "$PYTHON_BIN"
  elif command -v python3.11 >/dev/null 2>&1; then
    printf '%s\n' python3.11
  else
    printf '%s\n' python3
  fi
}

hh_activate_venv_if_present() {
  local venv_dir="$1"
  if [[ -f "$venv_dir/bin/activate" ]]; then
    source "$venv_dir/bin/activate"
  else
    echo "No virtual environment found at $venv_dir; using the current Python environment."
    echo "Run config/env.slurm first, or set VENV_DIR/PYTHON_BIN to use another environment."
  fi
}
