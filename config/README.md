# Slurm Configs

These scripts can run either with Slurm on a cluster or directly with `bash` on a local machine. They derive the repository root from the script location and avoid personal paths.

## Run Sequence

### 1. Set up the environment

Run this once before generation/decoding:

```bash
sbatch config/env.slurm
```

For local testing:

```bash
bash config/env.slurm
```

This creates `.venv-slurm` and installs the project with `pip install -e .`.

### 2A. Generate circuits for 2D plots

Use `generate.slurm`:

```bash
sbatch config/generate.slurm
```

The sweep values are controlled by `VAL_LIST` inside `generate.slurm`:

```bash
VAL_LIST=(
  1.00e-4 1.20e-4 1.43e-4
)
```

Edit this list to change the `noise_strength` values swept for 2D plots.

### 2B. Generate circuits for 3D plots

First generate the initial CSV sweep file:

```bash
python tools/threshold/generate_initial_sweep_range.py
```

This writes:

```text
data/analysis/swp_sug/sweep_range.csv
```

Then run:

```bash
sbatch config/generate_csv.slurm
```

`generate_csv.slurm` reads rows from the CSV and generates circuits for each row.

### 3. Decode generated circuits

Use `submit.slurm`:

```bash
sbatch config/submit.slurm
```

This uses a Slurm array, so each array task writes a separate CSV shard:

```text
OUT_DIR/parts/stats_<task_id>.csv
```

### 4. Merge shard CSV files

After all `submit.slurm` array jobs finish, merge the shard CSVs:

```bash
sbatch config/merge.slurm
```

This combines:

```text
OUT_DIR/parts/stats_*.csv
```

into:

```text
OUT_DIR/stats.csv
```

## Keep `OUT_DIR` Consistent

`OUT_DIR` is the most important path to keep consistent across generation, decoding, and merging.

For example:

```bash
OUT_DIR="$PWD/out/slurm/no_reset_2d" sbatch config/generate.slurm
OUT_DIR="$PWD/out/slurm/no_reset_2d" sbatch config/submit.slurm
OUT_DIR="$PWD/out/slurm/no_reset_2d" sbatch config/merge.slurm
```

For local testing:

```bash
OUT_DIR="$PWD/out/slurm/no_reset_2d" bash config/generate.slurm
OUT_DIR="$PWD/out/slurm/no_reset_2d" bash config/submit.slurm
OUT_DIR="$PWD/out/slurm/no_reset_2d" bash config/merge.slurm
```

## Cluster Settings To Modify

Each Slurm script has resource lines near the top. Change these to match your cluster and job size:

```bash
#SBATCH --job-name=HeavyHexStats
#SBATCH --time=4-00:00:00
#SBATCH --mem=128G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --array=0-31
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err
```

Common edits:

- `--job-name`: name shown in the queue.
- `--time`: wall-clock limit.
- `--mem`: memory request.
- `--cpus-per-task`: CPU count used by parallel generation or decoding.
- `--array`: number of Slurm array jobs. Increase for more parallelism.
- `--output` / `--error`: stdout/stderr log filename pattern.

These scripts do not set a fixed partition, account, or email. Add those if your cluster requires them:

```bash
sbatch --partition=compute --account=my_account config/submit.slurm
```

## Parameters In `generate.slurm`

Useful settings near the top of the script:

- `OUT_DIR`: output directory for generated circuits. Must match the `OUT_DIR` used by `submit.slurm` and `merge.slurm`.
- `TABLE`: breaker-table schedule file. For example, `data/analysis/tbls/breaker_table_d35.txt`. for detail, please check `data/analysis/tbls/README.md`.
- `ROUNDS`: number of rounds, usually `auto`.
- `BACKEND`: choose reset or no-reset circuits:
  - `modified_no_hook.unflagged.no_reset`
  - `modified_no_hook.unflagged.reset`
- `NOISE_MODEL`: noise model for this fixed-value sweep. Use `uniform` for the standard 2D workflow, or `si1000` if you want the SI1000 model supported by `tools/generate_circuit.py`.
- `STYLE`: style label written into circuit filename metadata. The default is `heavy_hex`, which is what the threshold tools expect.
- `VAL_LIST`: physical noise strengths to sweep for 2D plot data.
- `#SBATCH --array`: should cover the number of chunks you want for `VAL_LIST`.

## Parameters In `generate_csv.slurm`

Useful settings:

- `CSV`: input sweep CSV, defaulting to `data/analysis/swp_sug/sweep_range.csv`.
- `OUT_DIR`: output directory for generated circuits. Must match `submit.slurm` and `merge.slurm`.
- `TABLE`: breaker-table schedule file.
- `ROUNDS`: number of rounds, usually `auto`.
- `BACKEND`: choose reset or no-reset circuits:
  - `modified_no_hook.unflagged.no_reset`
  - `modified_no_hook.unflagged.reset`
- `STYLE`: style label written into circuit filename metadata. The default is `heavy_hex`, which is what the threshold tools expect. The 3D script uses the `corr_mixing` noise model through `scripts/pipeline/generate_csv_sweep.sh`. noise model detail can be found in `src/gen/_noise.py`
- `#SBATCH --array`: controls how the CSV rows are split across jobs.

CSV columns used by `generate_csv.slurm`:

- `x`: D2 / two-qubit depolarizing probability.
- `y`: measurement probability.
- `z`: one or more idle-noise values.
- `basis`: `X`, `Z`, or both depending on the row format.

## Parameters In `submit.slurm`

Useful settings:

- `OUT_DIR`: directory containing generated `circuits/`.
- `--cpus-per-task`: number of worker processes used by `sinter collect`.
- `--array`: number of shards used to split the circuit list.

The script writes per-task CSVs into `OUT_DIR/parts/`.

## Parameters In `merge.slurm`

Useful settings:

- `OUT_DIR`: the same output directory used by generation and submit.
- `PARTS_DIR`: directory containing `stats_*.csv`, defaulting to `OUT_DIR/parts`.
- `FINAL`: merged output CSV path, defaulting to `OUT_DIR/stats.csv`.

## Python Environment Notes

If your cluster uses environment modules, the scripts try to load `python/3.11.4` when `module` exists. Override this with:

```bash
PYTHON_MODULE=python/3.11 sbatch config/env.slurm
```

If your cluster does not use modules, the scripts use `python3.11` when available, otherwise `python3`. Override this with:

```bash
PYTHON_BIN=python3 bash config/env.slurm
```
