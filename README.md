# Dynamic Compass Code Numerics

This repository contains the numerical code for the paper **"Low-valency scalable quantum error correction with a dynamic compass code"**.

The main way to use this code is through the Slurm scripts in `config/`. Those scripts set up the environment, generate circuits, run decoding, and merge the resulting CSV files.

## Attribution

This repository is a lightly modified version of code from:

**"Less Bacon, More Threshold"**  
arXiv:2305.12046  
Original code: <https://github.com/Strilanc/more-bacon-less-threshold>

All credit for the original implementation goes to the authors of the Less Bacon, More Threshold project. In particular, the code in `src/gen/` is from the original project.

## How To Run

Start with the workflow documented in:

```text
config/README.md
```

The usual sequence is:

```bash
sbatch config/env.slurm
sbatch config/generate.slurm        # for 2D plot data
sbatch config/generate_csv.slurm    # for 3D plot data
sbatch config/submit.slurm
sbatch config/merge.slurm
```

For 3D plot data, first generate the initial sweep CSV:

```bash
python tools/threshold/generate_initial_sweep_range.py
```

Then run `config/generate_csv.slurm`.

## Paper Figures

The figures and supporting stats used in the manuscript are in:

```text
figures/manuscript/
```

These manuscript outputs are generated from the Slurm workflow in `config/`.

## Important Folders

- `config/`: Slurm entrypoints and run instructions.
- `scripts/pipeline/`: circuit-generation and decoding helper scripts called by Slurm.
- `scripts/workflows/threshold/`: threshold plotting workflow for 3D surfaces and rotation videos.
- `tools/threshold/`: Python tools used by the threshold workflow.
- `data/analysis/tbls/`: schedule tables used in the paper.
- `figures/manuscript/`: manuscript figures and supporting stats.
