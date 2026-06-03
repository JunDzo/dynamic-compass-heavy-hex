# Scripts

This directory has two main workflows:

- `pipeline/`: data-production workflow. These scripts generate `.stim` circuits and collect decoder statistics.
- `workflows/threshold/`: threshold-analysis workflow. These scripts turn collected stats into threshold estimates, 3D surfaces, refined sweep points, and rotation renders.

The plotting helpers in `plot/` are used directly for 2D plots.

## Data-Production Pipeline

Use `pipeline/` to generate and collect data.

- `pipeline/generate_uniform.sh`: generates circuits for a fixed list of physical noise values. This is the generation path used for the paper's 2D plot data.
- `pipeline/generate_csv_sweep.sh`: generates circuits from CSV-provided sweep points. This is the generation path used for the paper's 3D threshold-surface data.
- `pipeline/collect_stats.sh`: runs `sinter collect` on generated `.stim` circuits and writes stats CSV files.

Older names are kept as compatibility wrappers:

- `pipeline/step1_generate_circuits.sh` -> `pipeline/generate_uniform.sh`
- `pipeline/step1_generate_circuits_csv.sh` -> `pipeline/generate_csv_sweep.sh`
- `pipeline/step2_collect_stats.sh` -> `pipeline/collect_stats.sh`

### Generate 2D Plot Data

```bash
export Z_VAL="0.001"
export TABLE="data/analysis/tbls/breaker_table_d35.txt"
export ROUNDS="auto"
export BACKEND="modified_no_hook.unflagged.no_reset"
bash scripts/pipeline/generate_uniform.sh data
```

Optional variables:

- `DIAMETERS`: defaults to `9 11 13 15 17`.
- `BASIS`: defaults to `X`.
- `VAL`: accepted as an alias for `Z_VAL`.

### Generate 3D Sweep Data

```bash
export X="0.001"
export Y="0.01"
export Z_VAL="0.001 0.002"
export BASIS="X Z"
export TABLE="data/analysis/tbls/breaker_table_d35.txt"
export BACKEND="modified_no_hook.unflagged.no_reset"
bash scripts/pipeline/generate_csv_sweep.sh data
```

Optional variables:

- `DIAMETERS`: defaults to `9 11 13 15 17`.
- `BASE_PROB`: defaults to `2e-4`.

### Collect Stats

```bash
bash scripts/pipeline/collect_stats.sh data data/stats/run.csv
```

If the second argument is omitted, stats are written to `OUT_DIR/stats.csv`. If the second argument is a directory, stats are written to `DIRECTORY/stats.csv`.

## Plotting And Threshold Workflow

For 2D plots, use the Python files in `plot/` directly:

- `plot/logical_error_rate.py`: logical-error-rate plotting helper.
- `plot/threshold_slice.py`: threshold-slice plotting helper.

For 3D threshold plots, run the threshold workflow from the repository root:

```bash
bash scripts/workflows/threshold/00_find_thresholds.sh
bash scripts/workflows/threshold/01_plot_threshold_surfaces.sh
```

The remaining threshold workflow steps are for iteration and presentation:

- `02_generate_sweep_ranges.sh`: estimates additional sweep points for an iterative, finer 3D plot.
- `03_render_threshold_rotation.sh`: renders the 3D threshold plot rotation as an MP4.

The threshold shell scripts call the implementation scripts in `tools/threshold/`; see `tools/threshold/README.md` for the role of each Python file.

## Slurm Configs

Project-local Slurm entrypoints are stored in `config/`.

- `config/generate.slurm`: fixed value-list generation for 2D plot data.
- `config/generate_csv.slurm`: CSV-driven generation for 3D threshold-surface data.
- `config/submit.slurm`: stat collection across Slurm array shards.
- `config/env.slurm`: project-local Slurm virtual environment setup.
