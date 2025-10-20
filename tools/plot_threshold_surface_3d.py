#!/usr/bin/env python3
"""
plot_threshold_surface_3d.py

Read a CSV (e.g. produced by find_common_thresholds.py) and plot a 3D
threshold surface where:
  X = b1p (D2)
  Y = b2p (p_m)
  Z = p_idle (median of [p_low, p_high])
Vertical error bars show [p_low, p_high].

By default, data are plotted in log10 space (since mplot3d has no native log
axes). Tick labels are formatted back to scientific notation.

Usage:
  python plot_threshold_surface_3d.py \
      --csv common_thresholds.csv \
      --basis X \
      --out threshold_3d_X.png

Options:
  --csv   : input CSV (default: common_thresholds.csv)
  --basis : which basis to plot (X or Z). If omitted, plot all bases in one figure.
  --out   : output image file (PNG). If omitted, shows an interactive window.
"""

import argparse
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (needed to enable 3D)
import matplotlib.tri as mtri
from matplotlib.path import Path
from typing import Tuple, List, Optional

# --- Linear fit helpers ---
def _design_linear(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    X = np.asarray(X).ravel()
    Y = np.asarray(Y).ravel()
    return np.column_stack([
        np.ones_like(X),
        X, Y,
    ])

def _fit_linear_surface(X: np.ndarray, Y: np.ndarray, Z: np.ndarray) -> Tuple[np.ndarray, float, float]:
    """Return (coeffs, r2, rmse) for Z ≈ a + bX + cY."""
    M = _design_linear(X, Y)
    Z = np.asarray(Z).ravel()
    coef, residuals, rank, s = np.linalg.lstsq(M, Z, rcond=None)
    if Z.size > 0:
        y_hat = M @ coef
        ss_res = float(np.sum((Z - y_hat)**2))
        ss_tot = float(np.sum((Z - np.mean(Z))**2)) if Z.size > 1 else 0.0
        r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else 1.0
        rmse = float(np.sqrt(ss_res / max(1, Z.size)))
    else:
        r2, rmse = float("nan"), float("nan")
    return coef, r2, rmse

def _eval_linear_surface(coef: np.ndarray, X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    a, b, c = coef
    return a + b*X + c*Y
def _design_quadratic(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    X = np.asarray(X).ravel()
    Y = np.asarray(Y).ravel()
    return np.column_stack([
        np.ones_like(X),
        X, Y,
        X**2, X*Y, Y**2,
    ])

def _fit_quadratic_surface(X: np.ndarray, Y: np.ndarray, Z: np.ndarray) -> Tuple[np.ndarray, float, float]:
    """Return (coeffs, r2, rmse) for Z ≈ a + bX + cY + dX^2 + eXY + fY^2."""
    M = _design_quadratic(X, Y)
    Z = np.asarray(Z).ravel()
    # Least squares
    coef, residuals, rank, s = np.linalg.lstsq(M, Z, rcond=None)
    # Metrics
    if Z.size > 0:
        y_hat = M @ coef
        ss_res = float(np.sum((Z - y_hat)**2))
        ss_tot = float(np.sum((Z - np.mean(Z))**2)) if Z.size > 1 else 0.0
        r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else 1.0
        rmse = float(np.sqrt(ss_res / max(1, Z.size)))
    else:
        r2, rmse = float("nan"), float("nan")
    return coef, r2, rmse

def _eval_quadratic_surface(coef: np.ndarray, X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    a, b, c, d, e, f = coef
    return a + b*X + c*Y + d*X**2 + e*X*Y + f*Y**2


def parse_args():
    p = argparse.ArgumentParser(description="Plot 3D threshold points with error bars")
    p.add_argument("--csv", default="common_thresholds.csv", help="Input CSV")
    p.add_argument("--basis", choices=["X", "Z"], default="X", help="Filter to one basis")
    p.add_argument("--out", default=None, help="Output image filename (PNG)")
    p.add_argument(
        "--fit",
        choices=["none", "linear", "quadratic", "both"],
        default="linear",
        help="Fit a surface to the points in the current plotting space",
    )
    # Removed --log argument; always plot in linear scale.
    # Emit evenly distributed points on fitted surface
    p.add_argument("--emit-points", type=int, default=0, help="Number of evenly distributed (x,y,z) points to sample on the fitted surface (linear space)")
    p.add_argument("--emit-points-out", default=None, help="CSV file to write the sampled surface points")
    p.add_argument("--emit-points-mode", choices=["hull", "axis"], default="hull",
                   help="Sampling domain for --emit-points: 'hull' = convex hull of data; 'axis' = rectangle derived from fitted surface axis zero-crossings")
    p.add_argument("--emit-x-cap", type=float, default=None,
                   help="Optional positive cap for X when using --emit-points-mode axis (used if intercept is missing or infinite)")
    p.add_argument("--emit-y-cap", type=float, default=None,
                   help="Optional positive cap for Y when using --emit-points-mode axis (used if intercept is missing or infinite)")
    p.add_argument("--plot-emitted", action="store_true",
                   help="Overlay the emitted (x,y,z) points onto the 3D plot")
    return p.parse_args()

def _convex_hull(points: np.ndarray) -> np.ndarray:
    """Compute 2D convex hull (Andrew's monotone chain). Returns indices of hull in CCW order."""
    pts = np.asarray(points, float)
    if pts.shape[0] <= 3:
        return np.arange(pts.shape[0])
    pts_sorted = pts[np.lexsort((pts[:,1], pts[:,0]))]
    def cross(o, a, b):
        return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
    lower = []
    for p in pts_sorted:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(tuple(p))
    upper = []
    for p in reversed(pts_sorted):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(tuple(p))
    hull = np.array(lower[:-1] + upper[:-1])
    # map hull points back to original indices (handle duplicates by first match)
    idx = []
    used = set()
    for hp in hull:
        for i, pp in enumerate(pts):
            if i in used:
                continue
            if np.allclose(pp, hp):
                idx.append(i); used.add(i); break
    return np.array(idx, dtype=int)

def _mask_grid_to_hull(Xgrid: np.ndarray, Ygrid: np.ndarray, data_X: np.ndarray, data_Y: np.ndarray) -> np.ndarray:
    """Return boolean mask of points inside convex hull of (data_X, data_Y)."""
    P = np.column_stack([data_X.ravel(), data_Y.ravel()])
    hull_idx = _convex_hull(P)
    poly = P[hull_idx]
    path = Path(poly)
    pts = np.column_stack([Xgrid.ravel(), Ygrid.ravel()])
    inside = path.contains_points(pts)
    return inside.reshape(Xgrid.shape)

# --- Axis-intercept based range derivation (linear space) ---

def _axis_box_from_linear_fit(a: float, b: float, c: float,
                              x_cap: float=1.0,
                              y_cap: float=1.0,
                              ) -> Tuple[Tuple[float,float], Tuple[float,float]]:
    """Return [x_lo,x_hi], [y_lo,y_hi] rectangle where Z=a+bX+cY>0 along axes guidance.
    Strategy: use axis zero-crossings at Y=0 and X=0. If an intercept is missing or negative,
    fall back to user caps. If both are missing and caps absent, return NaNs.
    """
    x_lo, x_hi = float(0), float(1)
    y_lo, y_hi = float(0), float(1)

    # X-axis (Y=0): a + b X > 0
    if abs(b) < 1e-300: # cases for b = 0
        # No X dependence; if a<=0 then no positive region on axis
        if a > 0:
            x_lo, x_hi = 0.0 , 1.0
        else:
            x_lo, x_hi = float('nan'), float('nan')
    else:
        x0 = -a / b
        if b > 0:
            # X > x0
            x_lo = max(0.0, x0) 
            x_hi = x_cap
        else:
            # X < x0
            x_lo = 0.0
            x_hi = min(x0, x_cap)

    # Y-axis (X=0): a + c Y > 0
    if abs(c) < 1e-300:
        if a > 0:
            y_lo, y_hi = 0.0 , 1.0
        else:
            y_lo, y_hi = float('nan'), float('nan')
    else:
        y0 = -a / c
        if c > 0:
            y_lo = max(0.0, y0) 
            y_hi = y_cap
        else:
            y_lo = 0.0 
            y_hi = min(y0, y_cap)

    return (x_lo, x_hi), (y_lo, y_hi)


def _axis_box_from_quadratic_fit(q: np.ndarray,
                                 x_cap: Optional[float]=None,
                                 y_cap: Optional[float]=None,
                                 nonneg_only: bool=True) -> Tuple[Tuple[float,float], Tuple[float,float]]:
    """Rectangle guided by zeros of a+bX+cY+dX^2+eXY+fY^2 along axes.
    On X-axis: dX^2+bX+a=0; pick smallest positive root (if any) as upper bound.
    On Y-axis: fY^2+cY+a=0; pick smallest positive root (if any) as upper bound.
    If no positive root, use caps. Intersect with [0,∞) if nonneg_only.
    """
    a,b,c,d,e,f = [float(v) for v in q]
    # Roots
    xr = []
    if abs(d) < 1e-300:
        if abs(b) > 0:
            xr = [-a/b]
    else:
        D = b*b - 4*d*a
        if D >= 0:
            s = math.sqrt(D); xr = [(-b-s)/(2*d), (-b+s)/(2*d)]
    yr = []
    if abs(f) < 1e-300:
        if abs(c) > 0:
            yr = [-a/c]
    else:
        D2 = c*c - 4*f*a
        if D2 >= 0:
            s2 = math.sqrt(D2); yr = [(-c-s2)/(2*f), (-c+s2)/(2*f)]

    def _pick_pos(roots):
        pos = sorted([r for r in roots if r>0])
        return pos[0] if pos else None

    x_hi = _pick_pos(xr) or (x_cap if x_cap is not None else float('nan'))
    y_hi = _pick_pos(yr) or (y_cap if y_cap is not None else float('nan'))
    x_lo = 0.0 if nonneg_only else -math.inf
    y_lo = 0.0 if nonneg_only else -math.inf

    if not (isinstance(x_hi,float) and x_hi>0):
        x_lo, x_hi = float('nan'), float('nan')
    if not (isinstance(y_hi,float) and y_hi>0):
        y_lo, y_hi = float('nan'), float('nan')
    return (x_lo, x_hi), (y_lo, y_hi)

def main():
    args = parse_args()
    # Backwards compatibility: if --fitting is used and --fit left as "none", treat as quadratic
    df = pd.read_csv(args.csv)

    # Filter valid rows
    df = df.copy()
    df = df[df["status"] == "ok"]
    needed = {"b1p", "b2p", "p_low", "p_high"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing columns: {missing}")

    # basis filter    
    df = df[df["basis"] == args.basis]

    if df.empty:
        raise SystemExit("No rows to plot after filtering. Check basis/status.")

    # Compute z-mean, guard bad rows
    df = df.dropna(subset=["b1p", "b2p", "p_low", "p_high"])
    df = df[(df["p_low"] > 0) & (df["p_high"] > 0) & (df["b1p"] > 0) & (df["b2p"] > 0)]
    if df.empty:
        raise SystemExit("All rows invalid or non-positive. Nothing to plot.")

    x = df["b1p"].values  # D2
    y = df["b2p"].values  # p_m
    z_lo = df["p_low"].values
    z_hi = df["p_high"].values
    z_median = df["p_th"].values
    X, Y, Z, Zlo, Zhi = x, y, z_median, z_lo, z_hi

    X, Y, Z, Zlo, Zhi = x, y, z_median, z_lo, z_hi

    # Create figure
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Placeholders for optional emitted points
    emitted_pts_xy = None
    emitted_Zs = None
    sampled_df_emitted = None

    # Scatter the mean points
    ax.scatter(X, Y, Z, s=30, depthshade=True)

    # Draw vertical error bars (z only)
    for xi, yi, z0, z1 in zip(X, Y, Zlo, Zhi):
        # small horizontal tick at top/bottom is optional; we draw vertical segment only
        ax.plot([xi, xi], [yi, yi], [z0, z1], linewidth=1.5)

    # Optional: fit linear/quadratic surfaces in the current space
    if args.fit != "none":
        nx = ny = 40
        gx = np.linspace(X.min(), X.max(), nx)
        gy = np.linspace(Y.min(), Y.max(), ny)
        GX, GY = np.meshgrid(gx, gy)
        space = "linear space"

        if args.fit in ("linear", "both"):
            lcoef, lr2, lrmse = _fit_linear_surface(X, Y, Z)
            LZ = _eval_linear_surface(lcoef, GX, GY)
            mask = _mask_grid_to_hull(GX, GY, X, Y)
            Xg = GX[mask]; Yg = GY[mask]; Zg = LZ[mask]
            tri = mtri.Triangulation(Xg, Yg)
            ax.plot_trisurf(tri, Zg, alpha=0.25, antialiased=True, linewidth=0)
            print("Linear surface fit in", space)
            print("coefficients [a, b, c] =", lcoef)
            print(f"R^2 = {lr2:.4f}, RMSE = {lrmse:.4g}")

        if args.fit in ("quadratic", "both"):
            qcoef, qr2, qrmse = _fit_quadratic_surface(X, Y, Z)
            QZ = _eval_quadratic_surface(qcoef, GX, GY)
            mask = _mask_grid_to_hull(GX, GY, X, Y)
            Xg = GX[mask]; Yg = GY[mask]; Zg = QZ[mask]
            tri = mtri.Triangulation(Xg, Yg)
            ax.plot_trisurf(tri, Zg, alpha=0.4, antialiased=True, linewidth=0)
            print("Quadratic surface fit in", space)
            print("coefficients [a, b, c, d, e, f] =", qcoef)
            print(f"R^2 = {qr2:.4f}, RMSE = {qrmse:.4g}")

    # --- Precompute emitted surface points so we can plot them ---
    if getattr(args, 'emit_points', 0) and int(getattr(args, 'emit_points', 0)) > 0:
        n_pts = int(getattr(args, 'emit_points', 0))
        mode = getattr(args, 'emit_points_mode', 'hull')

        # Work in linear space for fitting and sampling
        Xlin, Ylin, Zlin = x.copy(), y.copy(), z_median.copy()

        # Linear fit used for reporting parameters (a,b,c) in output CSV
        lcoef_report, _, _ = _fit_linear_surface(Xlin, Ylin, Zlin)

        # Decide sampling rectangle
        if mode == 'axis':
            if args.fit in ("linear", "both"):
                lcoef_lin, _, _ = _fit_linear_surface(Xlin, Ylin, Zlin)
                (x_lo, x_hi), (y_lo, y_hi) = _axis_box_from_linear_fit(
                    float(lcoef_lin[0]), float(lcoef_lin[1]), float(lcoef_lin[2]),
                )
            elif args.fit == "quadratic":
                qcoef_lin, _, _ = _fit_quadratic_surface(Xlin, Ylin, Zlin)
                (x_lo, x_hi), (y_lo, y_hi) = _axis_box_from_quadratic_fit(
                    qcoef_lin,
                    x_cap=getattr(args, 'emit_x_cap', None),
                    y_cap=getattr(args, 'emit_y_cap', None),
                    nonneg_only=True,
                )
            else:
                x_lo=x_hi=y_lo=y_hi=float('nan')

            if not (isinstance(x_lo,float) and isinstance(x_hi,float) and x_lo<x_hi and x_lo>=0 and \
                    isinstance(y_lo,float) and isinstance(y_hi,float) and y_lo<y_hi and y_lo>=0):
                mode = 'hull'

        if mode == 'hull':
            xmin, xmax = float(Xlin.min()), float(Xlin.max())
            ymin, ymax = float(Ylin.min()), float(Ylin.max())

        # Build candidate XY grid
        if mode == 'axis':
            grid_side = max(2, int(np.sqrt(n_pts) + 0.5))
            gx = np.linspace(x_lo, x_hi, grid_side)
            gy = np.linspace(y_lo, y_hi, grid_side)
            GX, GY = np.meshgrid(gx, gy)
            cand = np.column_stack([GX.ravel(), GY.ravel()])
        else:
            grid_side = int(np.sqrt(n_pts) * 1.5)
            gx = np.linspace(xmin, xmax, grid_side)
            gy = np.linspace(ymin, ymax, grid_side)
            GX, GY = np.meshgrid(gx, gy)
            mask = _mask_grid_to_hull(GX, GY, Xlin, Ylin)
            cand = np.column_stack([GX.ravel(), GY.ravel()])[mask.ravel()]

        # Optional positivity filter at candidate stage for linear axis mode
        if mode == 'axis' and args.fit in ("linear", "both"):
            lcoef_lin_pf, _, _ = _fit_linear_surface(Xlin, Ylin, Zlin)
            a_pf, b_pf, c_pf = float(lcoef_lin_pf[0]), float(lcoef_lin_pf[1]), float(lcoef_lin_pf[2])
            # Keep only points satisfying a + b x + c y >= zmin
            zmin = 3e-5
            mask_pf = (a_pf + b_pf*cand[:,0] + c_pf*cand[:,1]) >= zmin
            cand = cand[mask_pf]

        if cand.shape[0] > 0:
            # Choose approximately evenly spaced subset
            if cand.shape[0] >= n_pts:
                idx = np.linspace(0, cand.shape[0]-1, n_pts).astype(int)
                emitted_pts_xy = cand[idx]
            else:
                emitted_pts_xy = cand

            # Evaluate Z on chosen model (in linear space)
            if args.fit in ("linear", "both"):
                coef, _, _ = _fit_linear_surface(Xlin, Ylin, Zlin)
                emitted_Zs = _eval_linear_surface(coef, emitted_pts_xy[:,0], emitted_pts_xy[:,1])
            elif args.fit == "quadratic":
                coef, _, _ = _fit_quadratic_surface(Xlin, Ylin, Zlin)
                emitted_Zs = _eval_quadratic_surface(coef, emitted_pts_xy[:,0], emitted_pts_xy[:,1])
            else:
                emitted_Zs = np.zeros(len(emitted_pts_xy))

            # Enforce Z >= emit_z_min
            zmin = 3e-5
            keep = np.isfinite(emitted_Zs) & (emitted_Zs >= zmin)
            if np.count_nonzero(keep) == 0:
                sampled_df_emitted = None
            else:
                emitted_pts_xy = emitted_pts_xy[keep]
                emitted_Zs = emitted_Zs[keep]
                sampled_df_emitted = pd.DataFrame({"x": emitted_pts_xy[:,0], "y": emitted_pts_xy[:,1], "z": emitted_Zs})
                # Add basis and linear fitting parameters (a,b,c)
                sampled_df_emitted["basis"] = args.basis
                a_rep, b_rep, c_rep = float(lcoef_report[0]), float(lcoef_report[1]), float(lcoef_report[2])
                sampled_df_emitted["a"] = a_rep
                sampled_df_emitted["b"] = b_rep
                sampled_df_emitted["c"] = c_rep

        # If requested, overlay on the current 3D plot
        if getattr(args, 'plot_emitted', False) and emitted_pts_xy is not None and emitted_Zs is not None:
            dispX, dispY, dispZ = emitted_pts_xy[:,0], emitted_pts_xy[:,1], emitted_Zs
            ax.scatter(dispX, dispY, dispZ, s=18, marker="x")

        # Emit sampled points to CSV or stdout preview (always emit if available)
        outpath = getattr(args, 'emit_points_out', None)
        if sampled_df_emitted is not None:
            if outpath:
                sampled_df_emitted.to_csv(outpath, index=False)
                print(f"Wrote {len(sampled_df_emitted)} evenly distributed surface points to {outpath}")
            else:
                print("Evenly distributed surface points (linear space):")
                print(sampled_df_emitted.head(min(20, len(sampled_df_emitted))).to_string(index=False))
        else:
            print("Warning: no emitted points to output.")

    # Axis labels
    ax.set_xlabel("D2 = b1p")
    ax.set_ylabel("p_m = b2p")
    ax.set_zlabel("p_idle (mean of [p_low, p_high])")

    # Linear mode: auto ticks; tighten limits with padding
    pad = 0.05
    ax.set_xlim(X.min() * (1 - pad), X.max() * (1 + pad))
    ax.set_ylim(Y.min() * (1 - pad), Y.max() * (1 + pad))
    ax.set_zlim(min(Zlo.min(), Z.min()), max(Zhi.max(), Z.max()))

    title_basis = ""
    if "basis" in df.columns:
        uniq_b = sorted(df["basis"].unique())
        if args.basis:
            title_basis = f" (basis={args.basis})"
        elif len(uniq_b) == 1:
            title_basis = f" (basis={uniq_b[0]})"
        else:
            title_basis = f" (bases={','.join(map(str, uniq_b))})"

    ax.set_title(f"Threshold surface: D2 vs p_m vs p_idle{title_basis}")

    plt.tight_layout()
    if args.out:
        fig.savefig(args.out, dpi=200)
        print(f"Saved 3D plot to {args.out}")
    else:
        plt.show()

if __name__ == "__main__":
    main()