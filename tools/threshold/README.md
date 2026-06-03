# Threshold Tools

This folder contains the Python implementations used by the threshold workflow in `scripts/workflows/threshold/`, plus a few manual helper utilities.

## Workflow-Used Tools

These files are called directly by the shell workflow.

- `find_approx_thresholds.py`: reads a Sinter stats CSV and estimates approximate threshold crossings for each `(basis, b1p, b2p)` group. It writes a threshold summary CSV with stable columns such as `basis`, `b1p`, `b2p`, `p_th`, and `status`.
- `count_thresholds.py`: summarizes an approximate-threshold CSV by basis, counting attempted points, successful threshold estimates, and points where no approximate threshold was found.
- `plot_threshold_surface_3d.py`: reads the threshold summary CSV and plots 3D threshold surfaces. It can fit linear or quadratic surfaces, save fit coefficients, and emit fitted surface points for follow-up sweeps.
- `expand_sweep_points.py`: takes emitted/fitted surface points and expands each one into a finer local sweep range. This is used by workflow step `02_generate_sweep_ranges.sh` to suggest the next iteration of CSV sweep points.
- `render_rotation.py`: renders a rotating MP4 or GIF view of a 3D threshold surface, used by workflow step `03_render_threshold_rotation.sh`.

## Manual Helper Tools

These are not called by the current shell workflows, but can be useful during analysis.

- `generate_initial_sweep_range.py`: creates an initial coarse `data/analysis/swp_sug/sweep_range.csv` grid over `x`, `y`, `z`, and basis values.
- `check_th_range.py`: checks the maximum spread between `p_low` and `p_high` in `data/approx_thresholds.csv`.
- `plot_threshold.py`: helper for plotting a single 2D logical-error-rate threshold slice from Sinter stats at fixed `(b1p, b2p)`.
- `solve_surface.py`: CLI/library helper for solving fitted linear or quadratic threshold-surface models for one unknown variable.

## Typical 3D Workflow

Run these from the repository root after stats have been collected:

```bash
bash scripts/workflows/threshold/00_find_thresholds.sh
bash scripts/workflows/threshold/01_plot_threshold_surfaces.sh
```

For iterative refinement and presentation:

```bash
bash scripts/workflows/threshold/02_generate_sweep_ranges.sh
bash scripts/workflows/threshold/03_render_threshold_rotation.sh
```

Step `00` uses `find_approx_thresholds.py` and `count_thresholds.py`.
Step `01` uses `plot_threshold_surface_3d.py`.
Step `02` uses `expand_sweep_points.py`.
Step `03` uses `render_rotation.py`.
