# Pipeline Scripts Documentation

This directory contains shell scripts for orchestrating the Heavy-Hex circuit generation and analysis pipeline.

## Overview

The pipeline consists of two main stages:

1. **Circuit Generation** (`step1_*.sh`) - Generate Stim circuit files with specified noise models
2. **Statistics Collection** (`step2_*.sh`) - Run decoders and collect error statistics

## Scripts

### `step1_generate_circuits.sh`

**Purpose**: Generate batches of Heavy-Hex circuits using GNU parallel

**Usage**:
```bash
bash scripts/pipeline/step1_generate_circuits.sh <output_directory>
```

**Environment Variables**:
- `PM` - Measurement error probability (p_m) for bias_slot_2
- `D2` - Two-qubit gate error probability for bias_slot_1

**Description**:
This script uses GNU parallel to generate multiple circuit files with different parameters simultaneously. It calls `tools/generate_circuit.py` for each parameter combination.

**Configuration**:
The script contains embedded parameter sweeps using parallel's syntax:
```bash
parallel --ungroup tools/generate_circuit.py \
    --diameter {1} \
    --rounds "auto" \
    --noise_strength {2} \
    --b {3} \
    ::: 7 9 11 13 \                    # code distances
    ::: 1e-5 2e-5 ... 2e-2 \           # noise strengths
    ::: X Z                             # bases
```

**Key Features**:
- Parallel execution for faster generation
- Automatic rounds calculation (4 × diameter)
- Support for correlated mixing noise model
- Reads breaker table from configured path
- Environment variable control for bias parameters

**Output**:
- Circuit files in `$OUT_DIR/circuits/`
- Filename format: `r=<rounds>,d=<diameter>,p=<noise>,...,b=<basis>.stim`

---

### `generate_error_rate_circuits.sh`

**Purpose**: Generate circuits for uniform error rate studies

**Usage**:
```bash
bash scripts/pipeline/generate_error_rate_circuits.sh <output_directory>
```

**Description**:
Specialized version of step1 for simple uniform noise studies. Generates circuits across a range of code distances with fixed noise parameters.

**Configuration**:
```bash
::: {3..25} \      # code distances from 3 to 25
::: 1e-3 \         # single noise strength
::: heavy_hex \    # code style
::: X Z            # both bases
```

**Use Case**:
- Initial error rate characterization
- Single-noise-parameter threshold finding
- Quick prototyping and testing

**Output**:
- Circuit files in `$OUT_DIR/circuits/`
- Simpler parameter space than step1_generate_circuits.sh

---

### `step2_collect_stats.sh`

**Purpose**: Run decoders on generated circuits and collect statistics

**Usage**:
```bash
bash scripts/pipeline/step2_collect_stats.sh <output_directory>
```

**Environment Variables**:
- `STATS_FILE` - Override default output CSV path (optional)
- `SLURM_CPUS_PER_TASK` - Number of parallel processes (defaults to 12)

**Description**:
Runs Sinter's collection framework with PyMatching decoder on all circuits in the output directory. Accumulates statistics across runs using the `--save_resume_filepath` feature.

**Key Parameters**:
```bash
--max_shots 10_000_000      # Maximum shots per circuit
--max_errors 1000           # Stop after 1000 errors observed
--processes 12              # Parallel decoding processes
--decoders pymatching       # Decoder to use
```

**Output**:
- Statistics CSV at `$OUT_DIR/stats.csv` (or `$STATS_FILE`)
- Resume-capable: can restart if interrupted
- Columns: circuit metadata, shots, errors, error rates, etc.

**Notes**:
- Uses Sinter's automatic metadata extraction (`--metadata_func auto`)
- Processes all `.stim` files in `$circuits_dir`
- Robust to interruption (can resume from checkpoint)

---

### `step2_collect_stats_tesseract.sh`

**Purpose**: Specialized statistics collection for Tesseract code studies

**Usage**:
```bash
bash scripts/pipeline/step2_collect_stats_tesseract.sh <output_directory>
```

**Description**:
Similar to `step2_collect_stats.sh` but may include Tesseract-specific decoder configurations or circuit filtering. Check the script for specific differences.

---

## Typical Workflow

### 1. Simple Uniform Noise Study

```bash
# Set output directory
export OUT_DIR="data/uniform_study"

# Generate circuits
bash scripts/pipeline/generate_error_rate_circuits.sh "$OUT_DIR"

# Run decoders
bash scripts/pipeline/step2_collect_stats.sh "$OUT_DIR"

# Analyze thresholds
python tools/find_common_thresholds.py \
    "$OUT_DIR/stats.csv" \
    data/analysis/thresholds.csv \
    linear
```

### 2. Correlated Mixing 3D Study

```bash
# Set output and environment parameters
export OUT_DIR="data/corr_study"
export PM=5e-4
export D2=7e-4

# Generate circuits
bash scripts/pipeline/step1_generate_circuits.sh "$OUT_DIR"

# Run decoders (use more shots for threshold precision)
bash scripts/pipeline/step2_collect_stats.sh "$OUT_DIR"

# Analyze thresholds
python tools/find_common_thresholds.py \
    "$OUT_DIR/stats.csv" \
    data/analysis/thresholds_3d.csv \
    quadratic

# Visualize 3D surface
python tools/plot_threshold_surface_3d.py \
    --csv data/analysis/thresholds_3d.csv \
    --basis X \
    --fit both \
    --out data/plots/surfaces/threshold_3d_X.png
```

### 3. HPC Execution with Slurm

```bash
# Submit cluster job
sbatch config/heavy_hex.slurm

# Or run locally with custom parameters
bash config/heavy_hex_local.slurm
```

## Integration with Slurm

The pipeline scripts are designed to be called from Slurm job files (`config/heavy_hex*.slurm`):

```bash
#!/bin/bash
#SBATCH --mem=32G
#SBATCH -c 128
# ... other Slurm directives ...

cd /path/to/HH-Modified
OUT_DIR=~/HH/out

# Step 1: Generate circuits
bash scripts/pipeline/step1_generate_circuits.sh "$OUT_DIR"

# Step 2: Collect stats
bash scripts/pipeline/step2_collect_stats.sh "$OUT_DIR"

# Optional: Post-processing
python tools/find_common_thresholds.py ...
```

## Performance Tips

1. **Parallel Generation**: Adjust parallel workers based on available cores
   ```bash
   export PARALLEL="-j 64"  # Use 64 parallel jobs
   ```

2. **Decoder Processes**: Match to available CPUs
   ```bash
   --processes $SLURM_CPUS_PER_TASK
   ```

3. **Checkpoint Resumption**: Always use `--save_resume_filepath` for long runs
   - Interrupted jobs can resume without starting over
   - Multiple runs accumulate statistics in the same file

4. **Noise Parameter Sweeps**: Adjust ranges in scripts for finer/coarser grids
   ```bash
   ::: 1e-5 2e-5 5e-5 1e-4 2e-4 5e-4 1e-3  # Logarithmic sweep
   ```

## Troubleshooting

**Problem**: `parallel: command not found`
- Install GNU parallel: `brew install parallel` (macOS) or `apt-get install parallel` (Linux)

**Problem**: Python modules not found
- Ensure PYTHONPATH includes `src/`: `export PYTHONPATH=src:$PYTHONPATH`

**Problem**: Out of memory during decoding
- Reduce `--processes` parameter
- Process fewer circuits at once by splitting into batches

**Problem**: Stats file corrupted
- Remove the checkpoint file and restart: `rm $OUT_DIR/stats.csv`

## Customization

To modify parameter sweeps, edit the parallel invocation in the scripts:

```bash
# Original
::: 7 9 11 13 \
```

To change noise models or add parameters, modify the `tools/generate_circuit.py` arguments:

```bash
--noise_model corr_mixing \
--base_prob {4} \
::: 1e-4 5e-4 1e-3 \  # Add sweep over base_prob
```

## See Also

- `config/README.md` - Slurm job configuration documentation
- `tools/generate_circuit.py` - Circuit generator CLI documentation
- Main project `README.md` - Overall workflow and examples
