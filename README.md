# Heavy-Hex Quantum Error Correction Threshold Analysis

This project implements tools for generating, simulating, and analyzing Heavy-Hex quantum error correction codes with various noise models to determine error correction thresholds.

## Project Structure

```
HH-Modified/
├── config/              # Slurm job configuration files
├── scripts/
│   ├── decoders/        # Decoder execution scripts (BeliefMatching, MWPF)
│   ├── pipeline/        # Circuit generation and data collection pipelines
│   └── visualization/   # Error rate plotting scripts
├── tools/               # Analysis and visualization tools
├── src/                 # Source code modules
│   ├── beliefmatching/  # Belief propagation decoder
│   ├── decoder/         # Decoder interfaces
│   ├── gen/             # Circuit generation utilities
│   ├── heavyhex/        # Heavy-Hex lattice implementation
│   └── utils/           # Utility functions
├── data/                # Output data directory
│   ├── circuits/        # Generated .stim circuit files
│   ├── stats/           # Decoder statistics CSV files
│   ├── plots/           # Output plots
│   │   ├── thresholds/  # Threshold analysis plots
│   │   ├── error_rates/ # Error rate plots
│   │   ├── surfaces/    # 3D threshold surfaces
│   │   └── breakers/    # Breaker selection visualizations
│   └── analysis/        # Analysis results (threshold CSVs)
│   │   ├── comm_th/     # from tools/find_common_thresholds.py
│   │   ├── swp_sug/     # from tools/plot_threshold_surface_3d.py
│   │   └── tbls/        # from tools/select_breakers_interactive.py
├── tests/               # Test files
└── out*/                # Legacy output directories
```

## Quick Start

### 1. Generate Circuits

Generate Heavy-Hex circuits with uniform noise:

```bash
python tools/generate_circuit.py \
    --diameter 5 \
    --rounds auto \
    --noise_model uniform \
    --noise_strength 0.001 \
    --b Z
```

Generate circuits with correlated mixing noise for 3D threshold studies:

```bash
python tools/generate_circuit.py \
    --diameter 7 \
    --rounds 28 \
    --noise_model corr_mixing \
    --noise_strength 0.001 \
    --base_prob 0.0002 \
    --sweep_slot idle \
    --bias_slot_1 D2 \
    --bias_slot_1_prob 0.0007 \
    --bias_slot_2 M \
    --bias_slot_2_prob 0.005 \
    --table_file data/analysis/breaker_table.txt
```

### 2. Run Decoders

Run MWPF decoder on generated circuits:

```bash
python scripts/decoders/run_mwpf_decoder.py \
    --circuits "data/circuits/*.stim" \
    --max-shots 10000000 \
    --max-errors 1000 \
    --processes 12 \
    --save-resume data/stats/stats.csv
```

Run BeliefMatching decoder:

```bash
python scripts/decoders/run_belief_matching_decoder.py \
    --circuits "data/circuits/*.stim" \
    --max-shots 10000000 \
    --max-errors 1000 \
    --max_bp_iters 20 \
    --processes 12 \
    --save-resume data/stats/stats.csv
```

### 3. Find Thresholds

Extract common thresholds from decoder statistics:

```bash
# Linear bracketing method
python tools/find_common_thresholds.py \
    data/stats/stats.csv \
    data/analysis/thresholds.csv \
    linear

# Quadratic collapse method
python tools/find_common_thresholds.py \
    data/stats/stats.csv \
    data/analysis/thresholds.csv \
    quadratic
```

### 4. Visualize Results

Plot 2D threshold slices:

```bash
python tools/plot_threshold_surface_2d.py
# (Edit settings in the file: CSV_PATH, BASIS, X_AXIS, Y_AXIS, etc.)
```

Plot 3D threshold surfaces:

```bash
python tools/plot_threshold_surface_3d.py \
    --csv data/analysis/thresholds.csv \
    --basis X \
    --fit linear \
    --out data/plots/surfaces/threshold_3d_X.png \
    --emit-points 100 \
    --emit-points-out data/analysis/surface_points_X.csv
```

Plot threshold slice for fixed parameters:

```bash
python tools/plot_threshold_slice.py
# (Edit settings in the file for specific (b1p, b2p) values)
```

## Tools Documentation

### Circuit Generation

**`tools/generate_circuit.py`** - Generate Heavy-Hex circuits with noise

- Supports uniform depolarizing and correlated mixing noise models
- Configurable code distance, rounds, and basis (X/Z)
- Optional breaker table for custom lattice configurations
- Outputs: `.stim` files in `data/circuits/`

### Threshold Analysis

**`tools/find_common_thresholds.py`** - Find common crossing thresholds

- Linear bracketing: Finds intervals where all distance pairs cross
- Quadratic collapse: Fits per-distance quadratics and minimizes variance
- Input: Sinter stats CSV from decoder runs
- Output: CSV with (b1p, b2p) → threshold mappings

**`tools/plot_threshold_slice.py`** - Plot logical error rate vs noise for fixed bias

- Visualizes threshold crossings for specific (b1p, b2p) slice
- Configurable log/linear scales
- Output: `threshold.pdf` in `data/plots/thresholds/`

**`tools/plot_threshold_surface_2d.py`** - Plot 2D threshold projections

- Visualize threshold as function of two parameters (fixing third)
- Supports error bars from threshold intervals
- Output: `threshold_slice_2d.png` in `data/plots/surfaces/`

**`tools/plot_threshold_surface_3d.py`** - Plot 3D threshold surfaces

- Full 3D visualization: (D2, p_m, p_idle)
- Linear and quadratic surface fitting
- Point emission for uniform sampling on fitted surface
- Output: PNG plots and optional CSV of surface points

### Breaker Configuration

**`tools/select_breakers_interactive.py`** - Interactive breaker selection GUI

- Click to select/deselect grid squares for breaker placement
- Supports multi-step configurations
- Output: Python literal file with complex coordinates

**`tools/visualize_breaker_table.py`** - Visualize breaker configurations

- Renders saved breaker tables as PDF plots
- Highlights selected breakers on grid
- Output: `visualized_step_*.pdf` in `data/plots/breakers/`

## Pipeline Scripts

### `scripts/pipeline/step1_generate_circuits.sh`

Generates a batch of circuits using GNU parallel. Called by Slurm scripts for HPC execution.

### `scripts/pipeline/step2_collect_stats.sh`

Runs Sinter decoder collection on all generated circuits. Uses PyMatching decoder by default.

### `scripts/pipeline/generate_error_rate_circuits.sh`

Specialized pipeline for uniform error rate studies with parameter sweeps.

## Configuration

### Slurm Job Files

**`config/heavy_hex.slurm`** - Cluster job configuration
- Multi-step pipeline: circuit generation → decoding → analysis
- Configured for HPC environments with module loading

**`config/heavy_hex_local.slurm`** - Local execution configuration
- CSV-based parameter sweep support
- Reads experiment parameters from `data/sweep_range.csv`
- Suitable for local development and testing

## Data Organization

All output data is organized under `data/`:

- **`data/circuits/`** - Generated Stim circuit files (`.stim`)
- **`data/stats/`** - Decoder statistics CSV files
- **`data/plots/thresholds/`** - Threshold slice plots
- **`data/plots/error_rates/`** - Logical error rate plots
- **`data/plots/surfaces/`** - 3D threshold surface visualizations
- **`data/plots/breakers/`** - Breaker configuration visualizations
- **`data/analysis/`** - Analysis results (threshold CSVs, surface points)

## Requirements

- Python 3.11+
- Stim (quantum circuit simulator)
- Sinter (Monte Carlo sampling framework)
- PyMatching / MWPF (decoders)
- NumPy, Pandas, Matplotlib, SciPy
- GNU parallel (for batch processing)

Install Python dependencies:
```bash
pip install -r requirements.txt  # (create this file if needed)
```

## Contributing

When adding new analysis tools:
1. Add comprehensive docstrings to all functions
2. Update output paths to use `data/` subdirectories
3. Add usage examples to this README
4. Document all CLI arguments with `--help` support

## License

[Specify license]

## Authors

Heavy-Hex Research Team

## References

[Add relevant publications and citations]
