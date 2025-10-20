# Tools Documentation

Analysis and visualization tools for Heavy-Hex quantum error correction threshold studies.

## Table of Contents

- [Circuit Generation](#circuit-generation)
- [Threshold Analysis](#threshold-analysis)
- [Visualization Tools](#visualization-tools)
- [Breaker Configuration](#breaker-configuration)
- [Output Conventions](#output-conventions)

---

## Circuit Generation

### `generate_circuit.py`

Generate Stim circuit files for Heavy-Hex codes with configurable noise models.

**Key Features**:
- Uniform depolarizing noise (single parameter)
- Correlated mixing noise (slot-specific error rates)
- Automatic rounds calculation (`--rounds auto` sets rounds = 4d)
- Custom breaker table support
- Metadata-rich filenames for tracking

**Basic Usage**:
```bash
python tools/generate_circuit.py \
    --diameter 5 \
    --rounds auto \
    --noise_model uniform \
    --noise_strength 0.001 \
    --b Z
```

**Advanced Usage** (3D threshold study):
```bash
python tools/generate_circuit.py \
    --diameter 11 \
    --rounds auto \
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

**Parameters**:
- `--diameter` (required): Code distance
- `--rounds`: Number of QEC rounds (int or "auto")
- `--noise_model`: "uniform" or "corr_mixing"
- `--noise_strength`: Primary noise parameter (swept variable)
- `--b, --basis`: Logical basis ("X" or "Z")
- `--table_file`: Path to breaker configuration
- `--out_dir`: Output directory (default: `data/circuits`)

**Noise Slots** (corr_mixing model):
- `D1`: Single-qubit gates
- `D2`: Two-qubit gates
- `M`: Measurements
- `idle`: Idle/storage errors

**Output**:
- Location: `data/circuits/`
- Format: `.stim` circuit files
- Naming: Encodes all parameters for metadata tracking

---

## Threshold Analysis

### `find_common_thresholds.py`

Extract common threshold crossings from decoder statistics.

**Purpose**: Find noise parameters where logical error rates for different code distances all intersect (indicating threshold behavior).

**Methods**:
1. **Linear Bracketing**: Scans adjacent p-intervals for pairwise order flips
2. **Quadratic Collapse**: Fits per-distance quadratics and minimizes variance

**Usage**:
```bash
# Linear method
python tools/find_common_thresholds.py \
    data/stats/stats.csv \
    data/analysis/thresholds.csv \
    linear

# Quadratic method
python tools/find_common_thresholds.py \
    data/stats/stats.csv \
    data/analysis/thresholds.csv \
    quadratic
```

**Arguments**:
1. Input stats CSV (from Sinter decoder runs)
2. Output thresholds CSV
3. Method: "linear" or "quadratic"

**Input Requirements**:
- Stats CSV must include columns: `json_metadata`, `shots`, `errors`
- Metadata must contain: `d`, `b`, `p` (or `p_idle`), `b1p`, `b2p`, `noise`, `basis`
- Requires ≥2 code distances and ≥2 p-values per group

**Output CSV Columns**:
- `basis`: X or Z
- `b1p`: Two-qubit error probability (D2 slot)
- `b2p`: Measurement error probability (M slot)
- `dists`: Comma-separated list of distances analyzed
- `p_low`, `p_high`: Threshold bracket (linear method)
- `p_th`: Point estimate (quadratic method stores in this field)
- `status`: "ok", "no_common_threshold", or "insufficient_data"
- `points_used`: Number of p-values in dataset

**Algorithm Notes**:

*Linear Method*:
- Finds intervals [p_low, p_high] where all distance pairs flip order
- Returns median of pairwise linear crossing estimates
- Fast and robust for well-sampled data

*Quadratic Method*:
- Fits f_d(p) = a_d + b_d p + c_d p² for each distance d
- Minimizes variance of {f_d(p_th)} to find collapse point
- Better for sparse data or smooth threshold regions

---

## Visualization Tools

### `plot_threshold_slice.py`

Plot logical error rate vs noise strength for a fixed (b1p, b2p) slice.

**Purpose**: Visualize how logical error rate varies with a single noise parameter for different code distances, revealing threshold crossings.

**Configuration** (edit within script):
```python
csv_path = "data/stats/stats.csv"
b1p = 0.0007  # Fix D2 error rate
b2p = 0.0027  # Fix measurement error rate
log_x = False  # Linear or log scale for x-axis
log_y = False  # Linear or log scale for y-axis
```

**Usage**:
```bash
# Edit settings in the file, then:
python tools/plot_threshold_slice.py
```

**Output**:
- File: `threshold.pdf`
- Recommended location: Move to `data/plots/thresholds/`

**Features**:
- Automatic grouping by code distance
- Configurable axes scales (log/linear)
- Binomial confidence intervals (from Sinter)
- Grid and legend for clarity

---

### `plot_threshold_surface_2d.py`

Plot 2D projections of threshold surface by fixing one parameter.

**Purpose**: Visualize threshold as a function of two noise parameters while holding a third fixed.

**Configuration** (edit within script):
```python
CSV_PATH = "data/analysis/thresholds.csv"
BASIS = "X"

# Choose two parameters for axes
X_AXIS = "b2p"    # e.g., measurement error
Y_AXIS = "p_idle"  # e.g., idle/storage error

# Fix the third parameter
FIX_AXIS = "b1p"   # e.g., two-qubit gate error
FIX_VALUE = 0.0004
REL_TOL = 0.05     # ±5% matching tolerance

# Scaling
LOG_X = False
LOG_Y = False
```

**Usage**:
```bash
# Edit settings, then:
python tools/plot_threshold_surface_2d.py
```

**Output**:
- File: `threshold_slice_2d.png`
- Recommended: `data/plots/surfaces/threshold_2d_<params>.png`

**Features**:
- Flexible axis assignment (any two of b1p, b2p, p_idle)
- Error bars from threshold intervals
- Tolerance-based slicing for fixed parameter
- Configurable log/linear scales

---

### `plot_threshold_surface_3d.py`

Create 3D visualizations of threshold surfaces with optional fitting.

**Purpose**: Full 3D representation of threshold surface: (D2, p_m, p_idle) with linear/quadratic fits and point sampling.

**Usage**:
```bash
python tools/plot_threshold_surface_3d.py \
    --csv data/analysis/thresholds.csv \
    --basis X \
    --fit linear \
    --out data/plots/surfaces/threshold_3d_X.png
```

**Advanced** (with surface sampling):
```bash
python tools/plot_threshold_surface_3d.py \
    --csv data/analysis/thresholds.csv \
    --basis X \
    --fit both \
    --emit-points 100 \
    --emit-points-mode hull \
    --emit-points-out data/analysis/surface_points_X.csv \
    --plot-emitted \
    --out data/plots/surfaces/threshold_3d_fitted_X.png
```

**Key Arguments**:
- `--csv`: Input threshold CSV
- `--basis`: X or Z
- `--fit`: "none", "linear", "quadratic", or "both"
- `--out`: Output PNG file
- `--emit-points N`: Sample N points on fitted surface
- `--emit-points-mode`: "hull" (data convex hull) or "axis" (intercept-based)
- `--emit-points-out`: CSV file for sampled (x,y,z,a,b,c) points
- `--plot-emitted`: Overlay sampled points on plot

**Output**:
- PNG plot in `data/plots/surfaces/`
- Optional CSV of surface points in `data/analysis/`
- Console output: fit coefficients, R², RMSE

**Fitting**:
- **Linear**: Z = a + bX + cY (plane fit)
- **Quadratic**: Z = a + bX + cY + dX² + eXY + fY² (paraboloid fit)

**Surface Sampling**:
Generates evenly-spaced (x,y,z) points on the fitted surface for:
- Uniform parameter exploration
- Experimental design guidance
- Interpolation validation

---

## Breaker Configuration

### `select_breakers_interactive.py`

Interactive GUI for selecting breaker configurations on Heavy-Hex grids.

**Purpose**: Manually select grid positions for "breaker" qubits that modify the Heavy-Hex lattice connectivity.

**Usage**:
```bash
python tools/select_breakers_interactive.py
```

**Interactive Workflow**:
1. Enter grid size `d` (code distance)
2. Enter number of steps (configurations) to create
3. For each step, a grid window appears:
   - Click squares to toggle selection (red = selected)
   - Only odd-parity squares (row+col odd) are selectable
   - Close window when done
4. Output saved to specified file (default: `table.txt`)

**Output**:
- File format: Python literal (list of lists of complex numbers)
- Example:
  ```python
  [
    [(1+0j), (2+1j), (3+2j)],  # Step 0
    [(1+1j), (3+1j), (5+3j)]   # Step 1
  ]
  ```
- Recommended location: `data/analysis/breaker_table_d<N>.txt`

**Use Case**:
Custom lattice configurations for studying:
- Connectivity requirements
- Decoder performance sensitivity
- Code variants

---

### `visualize_breaker_table.py`

Visualize saved breaker configurations as grid plots.

**Purpose**: Render breaker tables as PDF plots for verification and documentation.

**Usage**:
```bash
# Edit script to set input file and grid size
python tools/visualize_breaker_table.py
```

**Configuration** (edit within script):
```python
d = 7  # Grid size
visualize_table("data/analysis/breaker_table.txt")
```

**Output**:
- Files: `visualized_step_0.pdf`, `visualized_step_1.pdf`, ...
- Recommended location: Move to `data/plots/breakers/`

**Features**:
- Gray shading for invalid (even-parity) squares
- Red highlighting for selected breakers
- Grid lines and coordinates for reference
- One PDF per step in the table

---

## Output Conventions

### File Naming

**Circuit Files**:
```
r=<rounds>,d=<diameter>,p=<noise>,noise=<model>,c=<style>,b=<basis>,g=all.stim
```

**Plot Files**:
```
threshold_slice_<basis>_b1p<val>_b2p<val>.pdf
threshold_2d_<xparam>_vs_<yparam>_<basis>.png
threshold_3d_<basis>_<method>.png
error_rate_d<distance>_<noise_model>.pdf
```

**Data Files**:
```
stats_<experiment_id>.csv
thresholds_<analysis_type>.csv
surface_points_<basis>_<method>.csv
breaker_table_d<distance>.txt
```

### Directory Organization

All tools default to outputting in `data/` subdirectories:

```
data/
├── circuits/          # .stim files from generate_circuit.py
├── stats/             # decoder output CSVs
├── plots/
│   ├── thresholds/    # 1D threshold slices
│   ├── error_rates/   # logical error rate plots
│   ├── surfaces/      # 2D/3D threshold surfaces
│   └── breakers/      # breaker configuration visualizations
│── analysis/        # Analysis results (threshold CSVs)
│   ├── comm_th/     # from tools/find_common_thresholds.py
│   ├── swp_sug/     # from tools/plot_threshold_surface_3d.py
│   └── tbls/        # from tools/select_breakers_interactive.py
```

### CSV Metadata Standards

**Decoder Stats** (Sinter output):
- Required columns: `shots`, `errors`, `json_metadata`
- Metadata keys: `d`, `p`, `r`, `noise`, `c`, `b`, `b1p`, `b2p`, `sweep`, etc.

**Threshold Results**:
- Columns: `basis`, `b1p`, `b2p`, `dists`, `p_low`, `p_high`, `p_th`, `status`, `points_used`

**Surface Points**:
- Columns: `x` (b1p), `y` (b2p), `z` (p_idle), `basis`, `a`, `b`, `c` (fit coefficients)

---

## Tips and Best Practices

### Performance

1. **Generate circuits in parallel**: Use `scripts/pipeline/step1_generate_circuits.sh`
2. **Threshold finding**: Start with linear method (faster), refine with quadratic if needed
3. **Plot generation**: Use lower DPI for quick previews, high DPI for publications

### Workflow

1. Generate circuits → decode → find thresholds → visualize
2. Use consistent basis (X or Z) throughout a study
3. Keep threshold CSVs for reproducibility
4. Version control breaker tables with descriptive filenames

### Debugging

- Check circuit metadata: `stim circuit.stim | head -20`
- Validate stats CSV: Ensure all required metadata fields present
- Test threshold finding on small subset first
- Use `--help` on all tools for full argument lists

---

## See Also

- `../scripts/pipeline/README.md` - Pipeline execution documentation
- `../README.md` - Project overview and quick start
- `../config/` - Slurm job configuration

## Support

For issues or questions:
1. Check tool `--help` output
2. Review this documentation
3. Inspect example output files in `data/`
4. Contact: [research team contact]
