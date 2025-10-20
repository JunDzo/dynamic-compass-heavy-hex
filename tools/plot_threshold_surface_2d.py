

#!/usr/bin/env python3
"""
plot_threshold_surface_2d.py

A simple 2D plotter for threshold slices without CLI args.
You choose which two quantities to put on the axes and fix the third.

Data expectations (from find_common_thresholds.py output):
  columns: basis, b1p, b2p, p_low, p_high, (optional p_avg), status

Conventions:
  - b1p = D2 (two-qubit depolarizing)
  - b2p = p_m (measurement error)
  - p_idle is represented by the mean of [p_low, p_high]; we compute it as p_avg

Edit the SETTINGS block and run.
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =============================
# SETTINGS (edit these)
# =============================
CSV_PATH = "data/common_thresholds_3.csv"   # input table
BASIS = "X"                           # "X" or "Z" or None for all

# Choose axis names from {"b1p", "b2p", "p_idle"}
X_AXIS = "b2p"                        # e.g., "b1p" (D2) or "b2p" (p_m) or "p_idle"
Y_AXIS = "p_idle"                     # e.g., "p_idle" vs X

# Fix the third variable by value (only used if the third exists)
# Example: if X_AXIS=b1p and Y_AXIS=p_idle, then we fix b2p here.
FIX_AXIS = "b1p"                      # which variable to hold fixed
FIX_VALUE = 0.0003971485751819                    # target value for the fixed variable

# Match tolerance for the fixed variable (relative). 0.0 matches exactly.
REL_TOL = 0.05                         # 5% window by default

# Log scaling toggles (no CLI; just flip here)
LOG_X = False
LOG_Y = False

# Output
OUT_PNG = "data/plots/surfaces/threshold_slice_2d.png"
DPI = 300
# =============================

AXIS_LABELS = {
    "b1p": "D2 probability (b1p)",
    "b2p": "Measurement probability (p_m = b2p)",
    "p_idle": "Idle threshold (mean of [p_low, p_high])",
}


def _prepare_dataframe(csv_path: str, basis: str | None) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    # compute p_avg if missing
    if "p_avg" not in df.columns and {"p_low", "p_high"}.issubset(df.columns):
        df["p_avg"] = (df["p_low"] + df["p_high"]) / 2.0

    # keep only good rows
    req_cols = {"b1p", "b2p", "p_low", "p_high", "p_avg", "status"}
    missing = req_cols - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing columns: {sorted(missing)}")

    df = df[df["status"] == "ok"].copy()
    # basis filter if present
    if basis is not None and "basis" in df.columns:
        df = df[df["basis"] == basis].copy()

    # guard non-positive values
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["b1p", "b2p", "p_low", "p_high", "p_avg"])
    df = df[(df["b1p"] > 0) & (df["b2p"] > 0) & (df["p_low"] > 0) & (df["p_high"] > 0) & (df["p_avg"] > 0)]
    return df


def _apply_slice(df: pd.DataFrame, x_axis: str, y_axis: str, fix_axis: str | None, fix_value: float | None, rel_tol: float) -> pd.DataFrame:
    if x_axis == y_axis:
        raise ValueError("X_AXIS and Y_AXIS must be different.")

    # Map p_idle logical name to p_avg column
    def col_name(name: str) -> str:
        return "p_avg" if name == "p_idle" else name

    x_col = col_name(x_axis)
    y_col = col_name(y_axis)

    # If there is a third variable to fix, do it with a relative tolerance
    if fix_axis is not None:
        f_col = col_name(fix_axis)
        if f_col not in df.columns:
            raise ValueError(f"Unknown FIX_AXIS='{fix_axis}'. Use 'b1p', 'b2p', or 'p_idle'.")
        if fix_value is None:
            raise ValueError("FIX_VALUE must be set when FIX_AXIS is specified.")

        # Relative window around the target value
        lo = fix_value * (1.0 - rel_tol)
        hi = fix_value * (1.0 + rel_tol)
        df = df[(df[f_col] >= lo) & (df[f_col] <= hi)].copy()

    # Sort by x for pretty plotting
    df = df.sort_values(by=x_col)
    # Compute error bars for p_idle if used
    df["y"] = df[y_col]
    df["x"] = df[x_col]

    if y_axis == "p_idle":
        # vertical error bar from p_low to p_high
        df["yerr_low"] = df["y"] - df["p_low"]
        df["yerr_high"] = df["p_high"] - df["y"]
    else:
        df["yerr_low"] = 0.0
        df["yerr_high"] = 0.0

    if x_axis == "p_idle":
        # horizontal error bar if x is p_idle
        df["xerr_low"] = df["x"] - df["p_low"]
        df["xerr_high"] = df["p_high"] - df["x"]
    else:
        df["xerr_low"] = 0.0
        df["xerr_high"] = 0.0

    return df[["x", "y", "xerr_low", "xerr_high", "yerr_low", "yerr_high", "b1p", "b2p", "p_avg", "p_low", "p_high"]]


def main():
    df = _prepare_dataframe(CSV_PATH, BASIS)

    # Determine the fixed variable automatically if not set
    chosen = {X_AXIS, Y_AXIS}
    third = ({"b1p", "b2p", "p_idle"} - chosen)
    fix_axis = FIX_AXIS
    if fix_axis is None:
        if len(third) == 1:
            fix_axis = next(iter(third))
        else:
            raise ValueError("Ambiguous FIX_AXIS. Set FIX_AXIS explicitly.")

    data = _apply_slice(df, X_AXIS, Y_AXIS, fix_axis, FIX_VALUE, REL_TOL)

    if data.empty:
        raise SystemExit("No data after slicing. Adjust FIX_AXIS/FIX_VALUE/REL_TOL or axes.")

    fig, ax = plt.subplots(figsize=(8, 6))

    # Build asymmetrical error bars for x/y if p_idle used
    x = data["x"].to_numpy()
    y = data["y"].to_numpy()
    xerr = np.vstack([data["xerr_low"].to_numpy(), data["xerr_high"].to_numpy()])
    yerr = np.vstack([data["yerr_low"].to_numpy(), data["yerr_high"].to_numpy()])

    # Plot with error bars
    ax.errorbar(
        x, y,
        xerr=None if not np.any(xerr) else xerr,
        yerr=None if not np.any(yerr) else yerr,
        fmt="o", linewidth=1, capsize=3
    )

    # Axis labels
    ax.set_xlabel(AXIS_LABELS.get(X_AXIS, X_AXIS))
    ax.set_ylabel(AXIS_LABELS.get(Y_AXIS, Y_AXIS))

    # Log scaling if requested
    if LOG_X:
        ax.set_xscale("log")
    if LOG_Y:
        ax.set_yscale("log")

    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()

    # Ensure output directory exists
    import os
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)

    fig.savefig(OUT_PNG, dpi=DPI, bbox_inches='tight')
    print(f"✓ Saved 2D threshold slice to {OUT_PNG} with {len(data)} points.")
    plt.show()


if __name__ == "__main__":
    main()