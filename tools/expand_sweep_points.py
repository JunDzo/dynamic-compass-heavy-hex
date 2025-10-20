
#!/usr/bin/env python3
"""

Example:
  python tools/expand_sweep_points.py \
    --in emitted_points.csv \
    --out expanded_jobs.csv \
    --eta 3e-5 \
    --n 9 \
    --clip-low 0 \
    --diameters "7 9 11 13" \
    --base-prob 2e-4 \
    --table-file /bucket/ElkoussU/jun/HH-Modified/out/plot/t1-th/table.txt

"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def linspace_inclusive(a: float, b: float, n: int) -> np.ndarray:
    """Like np.linspace but always includes both endpoints and returns float array."""
    if n <= 1:
        return np.array([float(a)], dtype=float)
    return np.linspace(float(a), float(b), int(n), dtype=float)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Expand z ranges for sweep from emitted points CSV",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # IO
    p.add_argument("--in", dest="in_path", required=True, help="Input CSV with x,y,z,basis")
    p.add_argument("--out", dest="out_path", required=True, help="Output CSV for GNU parallel")

    # Range spec
    p.add_argument("--eta", type=float, required=True, help="Half-width for z: create [z-eta, z+eta]")
    p.add_argument("--n", type=int, default=9, help="Number of points in the z range per row")
    p.add_argument("--clip-low", type=float, default=0.0, help="Minimum allowed z after expansion")

    # Column names (in case user renamed)
    p.add_argument("--x-col", default="x", help="Column name for x (D2)")
    p.add_argument("--y-col", default="y", help="Column name for y (PM)")
    p.add_argument("--z-col", default="z", help="Column name for z (idle)")
    p.add_argument("--basis-col", default="basis", help="Column name for basis")

    # Optional Cartesian expansion over diameters
    p.add_argument("--diameters", default=None, help="Space-separated list, e.g. '7 9 11 13'")

    # Convenience constants to carry through to parallel
    p.add_argument("--base-prob", type=float, default=None, help="Constant base_prob column to add")
    p.add_argument("--out-dir", default=None, help="Constant out_dir column to add")
    p.add_argument("--table-file", default=None, help="Constant table_file column to add")

    return p


def main() -> None:
    args = build_parser().parse_args()

    # Read and validate
    in_path = Path(args.in_path)
    if not in_path.exists():
        raise SystemExit(f"Input CSV not found: {in_path}")
    df = pd.read_csv(in_path)

    req = [args.x_col, args.y_col, args.z_col, args.basis_col]
    missing = [c for c in req if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns in input CSV: {missing}")

    x = df[args.x_col].to_numpy(float)
    y = df[args.y_col].to_numpy(float)
    z = df[args.z_col].to_numpy(float)
    basis = df[args.basis_col].astype(str).to_numpy()

    n = int(args.n)
    if n <= 0:
        raise SystemExit("--n must be >= 1")
    eta = float(args.eta)
    if eta < 0:
        raise SystemExit("--eta must be >= 0")
    zclip = float(args.clip_low)

    # Compute per-row ranges
    z_lo = np.maximum(z - eta, zclip)
    z_hi = z + eta

    # Expand rows
    rows = []
    for i in range(len(df)):
        z_grid = linspace_inclusive(z_lo[i], z_hi[i], n)
        # Format z as space-separated string without brackets/commas for GNU parallel loops
        z_str = " ".join(np.format_float_positional(v, trim='-') for v in z_grid)
        rows.append({
                "x": x[i],
                "y": y[i],
                "z": z_str,
                "basis": basis[i],
            })


    out_df = pd.DataFrame(rows, columns=["x", "y", "z", "basis"])  # stable order

    # Optional expansion over diameters
    if args.diameters:
        # Normalize to space-separated integers
        diam_str = " ".join(str(int(tok)) for tok in str(args.diameters).split())
        out_df["diameters"] = diam_str

    # Optional constant columns
    if args.base_prob is not None:
        out_df["base_prob"] = float(args.base_prob)
    if args.table_file is not None:
        out_df["table_file"] = str(args.table_file)

    # Write
    out_path = Path(args.out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)

    # Minimal summary
    print(f"Wrote {len(out_df)} rows to {out_path}")
    print("Columns:", ",".join(out_df.columns))


if __name__ == "__main__":
    main()