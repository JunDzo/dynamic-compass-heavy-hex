#!/usr/bin/env python3
"""
plot_threshold_surface_3d.py

Read a CSV (e.g. produced by find_approx_thresholds.py) and plot a 3D
threshold surface where:
  X = b1p (D2)
  Y = b2p (p_m)
  Z = p_idle (p_th in CSV)

By default, data are plotted in log10 space (since mplot3d has no native log
axes). Tick labels are formatted back to scientific notation.

Usage:
  python plot_threshold_surface_3d.py \
      --csv approx_thresholds.csv \
      --basis X \
      --out threshold_3d_X.png

Options:
  --csv   : input CSV (default: approx_thresholds.csv)
  --basis : which basis to plot (X or Z). If omitted, plot all bases in one figure.
  --out   : output image file (PNG). If omitted, shows an interactive window.
"""

import argparse
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MultipleLocator
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (needed to enable 3D)
import matplotlib.tri as mtri
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.path import Path
from typing import Tuple, List, Optional
import json

def _design_quadratic(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    X = np.asarray(X).ravel()
    Y = np.asarray(Y).ravel()
    return np.column_stack([
        np.ones_like(X),
        X, Y,
        X**2, X*Y, Y**2,
    ])

def _fit_metrics(Z: np.ndarray, y_hat: np.ndarray) -> Tuple[float, float]:
    Z = np.asarray(Z).ravel()
    y_hat = np.asarray(y_hat).ravel()
    if Z.size == 0:
        return float("nan"), float("nan")
    ss_res = float(np.sum((Z - y_hat) ** 2))
    ss_tot = float(np.sum((Z - np.mean(Z)) ** 2)) if Z.size > 1 else 0.0
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    rmse = float(np.sqrt(ss_res / max(1, Z.size)))
    return r2, rmse


def _quadratic_axis_roots(coef: np.ndarray, eps: float = 1e-15) -> Tuple[List[float], List[float]]:
    """Return real roots on y=0 and x=0 for a+bX+cY+dX^2+eXY+fY^2."""
    a, b, c, d, _, f = [float(v) for v in coef]

    def _real_roots(q2: float, q1: float, q0: float) -> List[float]:
        if abs(q2) < eps:
            if abs(q1) < eps:
                return []
            return [-q0 / q1]
        disc = q1 * q1 - 4.0 * q2 * q0
        if disc < -eps:
            return []
        s = math.sqrt(max(0.0, disc))
        return [(-q1 - s) / (2.0 * q2), (-q1 + s) / (2.0 * q2)]

    return _real_roots(d, b, a), _real_roots(f, c, a)


def _positive_crossing_floor(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr) & (arr > 0.0)]
    if arr.size == 0:
        return 1e-12
    return max(float(np.max(arr)) * 1e-6, 1e-12)


def _quadratic_has_positive_axis_crossings(coef: np.ndarray, eps: float = 1e-12) -> bool:
    a, b, c, d, _, f = [float(v) for v in coef]
    disc_x = b * b - 4.0 * d * a
    disc_y = c * c - 4.0 * f * a
    if not (np.isfinite(a) and a > eps and disc_x >= -eps and disc_y >= -eps):
        return False
    x_roots, y_roots = _quadratic_axis_roots(coef, eps=eps)
    return any(r > eps for r in x_roots) and any(r > eps for r in y_roots)


def _positive_axis_guess(coef: np.ndarray, X: np.ndarray, Y: np.ndarray, Z: np.ndarray, x_floor: float, y_floor: float) -> np.ndarray:
    coef0 = np.asarray(coef, dtype=float).copy()
    eps = 1e-12
    if not np.isfinite(coef0[0]) or coef0[0] <= eps:
        positive_z = np.asarray(Z, dtype=float)
        positive_z = positive_z[np.isfinite(positive_z) & (positive_z > eps)]
        coef0[0] = float(np.min(positive_z)) if positive_z.size else eps

    x_roots, y_roots = _quadratic_axis_roots(coef0, eps=eps)
    x_pos = sorted(r for r in x_roots if np.isfinite(r) and r > x_floor)
    y_pos = sorted(r for r in y_roots if np.isfinite(r) and r > y_floor)
    x_cross = x_pos[0] if x_pos else float(np.max(X)) if np.asarray(X).size else 1.0
    y_cross = y_pos[0] if y_pos else float(np.max(Y)) if np.asarray(Y).size else 1.0
    x_cross = max(x_cross, x_floor)
    y_cross = max(y_cross, y_floor)

    # Make the starting point satisfy the positive root equalities.
    coef0[1] = -(coef0[0] + coef0[3] * x_cross * x_cross) / x_cross
    coef0[2] = -(coef0[0] + coef0[5] * y_cross * y_cross) / y_cross
    return np.r_[coef0, x_cross, y_cross]


def _fit_quadratic_surface_constrained(
    X: np.ndarray,
    Y: np.ndarray,
    Z: np.ndarray,
    p0: np.ndarray,
    maxiter: int = 200000,
) -> np.ndarray:
    try:
        from scipy.optimize import minimize
    except ImportError as exc:
        raise SystemExit("constrained quadratic fit requires scipy.optimize.minimize") from exc

    X = np.asarray(X, dtype=float).ravel()
    Y = np.asarray(Y, dtype=float).ravel()
    Z = np.asarray(Z, dtype=float).ravel()
    M = _design_quadratic(X, Y)
    x_floor = _positive_crossing_floor(X)
    y_floor = _positive_crossing_floor(Y)
    start = _positive_axis_guess(p0, X, Y, Z, x_floor=x_floor, y_floor=y_floor)
    eps = 1e-12

    def objective(v: np.ndarray) -> float:
        residuals = M @ v[:6] - Z
        return float(np.dot(residuals, residuals))

    def x_axis_eq(v: np.ndarray) -> float:
        a, b, _, d, _, _ = v[:6]
        x_cross = v[6]
        return float(a + b * x_cross + d * x_cross * x_cross)

    def y_axis_eq(v: np.ndarray) -> float:
        a, _, c, _, _, f = v[:6]
        y_cross = v[7]
        return float(a + c * y_cross + f * y_cross * y_cross)

    def disc_x(v: np.ndarray) -> float:
        a, b, _, d, _, _ = v[:6]
        return float(b * b - 4.0 * d * a)

    def disc_y(v: np.ndarray) -> float:
        a, _, c, _, _, f = v[:6]
        return float(c * c - 4.0 * f * a)

    constraints = [
        {"type": "eq", "fun": x_axis_eq},
        {"type": "eq", "fun": y_axis_eq},
        {"type": "ineq", "fun": disc_x},
        {"type": "ineq", "fun": disc_y},
    ]
    bounds = [
        (eps, None),
        (None, None),
        (None, None),
        (None, None),
        (None, None),
        (None, None),
        (x_floor, None),
        (y_floor, None),
    ]
    result = minimize(
        objective,
        start,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        options={"maxiter": int(maxiter), "ftol": 1e-15, "disp": False},
    )
    if not result.success:
        raise RuntimeError(f"constrained quadratic fit failed: {result.message}")
    coef = np.asarray(result.x[:6], dtype=float)
    if not _quadratic_has_positive_axis_crossings(coef):
        raise RuntimeError("constrained quadratic fit did not produce positive x/y/z axis crossings")
    return coef


def _fit_quadratic_surface(
    X: np.ndarray,
    Y: np.ndarray,
    Z: np.ndarray,
    method: str = "lstsq",
    p0: Optional[np.ndarray] = None,
    fit_bounds: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    maxfev: int = 200000,
) -> Tuple[np.ndarray, float, float]:
    """Return (coeffs, r2, rmse) for Z ≈ a + bX + cY + dX^2 + eXY + fY^2."""
    M = _design_quadratic(X, Y)
    Z = np.asarray(Z).ravel()
    coef, residuals, rank, s = np.linalg.lstsq(M, Z, rcond=None)
    if _quadratic_has_positive_axis_crossings(coef):
        fit_kind = "unconstrained"
    else:
        print("Unconstrained quadratic fit has no positive real crossing on all axes; retrying with constrained fit.")
        coef = _fit_quadratic_surface_constrained(X, Y, Z, p0=coef, maxiter=maxfev)
        fit_kind = "constrained"

    y_hat = _eval_quadratic_surface(coef, np.asarray(X), np.asarray(Y))
    r2, rmse = _fit_metrics(Z, y_hat)
    print(f"quadratic fit kind = {fit_kind}")
    return coef, r2, rmse


def _eval_quadratic_surface(coef: np.ndarray, X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    a, b, c, d, e, f = coef
    return a + b*X + c*Y + d*X**2 + e*X*Y + f*Y**2


def _fit_surface(
    model: str,
    X: np.ndarray,
    Y: np.ndarray,
    Z: np.ndarray,
    method: str = "lstsq",
    p0: Optional[np.ndarray] = None,
    fit_bounds: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    maxfev: int = 200000,
) -> Tuple[np.ndarray, float, float]:
    if model == "quadratic":
        return _fit_quadratic_surface(
            X, Y, Z, method=method, p0=p0, fit_bounds=fit_bounds, maxfev=maxfev
        )
    raise ValueError(f"Unsupported fit model: {model}. Only 'quadratic' is available.")


def _eval_surface(model: str, coef: np.ndarray, X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    if model == "quadratic":
        return _eval_quadratic_surface(coef, X, Y)
    raise ValueError(f"Unsupported fit model: {model}. Only 'quadratic' is available.")


def _parse_lim_pair(text: Optional[str], name: str) -> Optional[Tuple[float, float]]:
    if text is None:
        return None
    parts = [p.strip() for p in str(text).split(",")]
    if len(parts) != 2:
        raise SystemExit(f"--{name} must be 'min,max'")
    lo = float(parts[0])
    hi = float(parts[1])
    if not (np.isfinite(lo) and np.isfinite(hi) and hi > lo):
        raise SystemExit(f"--{name} must have finite values and max > min")
    return (lo, hi)


def _tick_scale_suffix(scale: float) -> str:
    if not np.isfinite(scale) or scale <= 0 or abs(scale - 1.0) < 1e-12:
        return ""
    exp = np.log10(scale)
    if abs(exp - round(exp)) < 1e-10:
        return rf" ($\times 10^{{{-int(round(exp))}}}$)"
    return rf" ($\times {scale:g}$)"


def _auto_tick_scale(lo: float, hi: float) -> float:
    vmax = max(abs(float(lo)), abs(float(hi)))
    if not np.isfinite(vmax) or vmax <= 0:
        return 1.0
    # Choose power-of-10 so the max shown value is roughly in [1, 10).
    exp = int(np.floor(np.log10(10.0 / vmax)))
    if exp < 2:
        exp = 2
    return float(10.0 ** exp)


def parse_args():
    p = argparse.ArgumentParser(description="Plot 3D threshold points")
    p.add_argument("--csv", default="approx_thresholds.csv", help="Input CSV")
    p.add_argument("--basis", choices=["X", "Z"], default="X", help="Filter to one basis")
    p.add_argument("--out", default=None, help="Output image filename (PNG)")
    p.add_argument("--x-axis", choices=["b1p", "b2p", "p_th"], default="b1p",
                   help="Quantity for X axis")
    p.add_argument("--y-axis", choices=["b1p", "b2p", "p_th"], default="b2p",
                   help="Quantity for Y axis")
    p.add_argument("--z-axis", choices=["b1p", "b2p", "p_th"], default="p_th",)
    p.add_argument("--x-label", default=None, help="Custom label for X axis (overrides default mapping)")
    p.add_argument("--y-label", default=None, help="Custom label for Y axis (overrides default mapping)")
    p.add_argument("--z-label", default=None, help="Custom label for Z axis (overrides default mapping)")
    p.add_argument("--xlim", "--x-lim", dest="xlim", default=None, help="X axis limits as 'min,max'")
    p.add_argument("--ylim", "--y-lim", dest="ylim", default=None, help="Y axis limits as 'min,max'")
    p.add_argument("--zlim", "--z-lim", dest="zlim", default=None, help="Z axis limits as 'min,max'")
    p.add_argument("--elev", type=float, default=14.0, help="3D view elevation angle in degrees")
    p.add_argument("--azim", type=float, default=45.0, help="3D view azimuth angle in degrees")
    p.add_argument("--tick-scale", type=float, default=0.0,
                   help="Global tick scale override; set <=0 to auto-scale each axis independently")
    # Removed --log argument; always plot in linear scale.
    p.add_argument("--fit-model", choices=["none", "quadratic"], default="quadratic",
                   help="Surface model used for fitting/plotting/emitting")
    p.add_argument("--fit-method", choices=["lstsq", "curve_fit"], default="lstsq",
                   help="Fitting backend to estimate coefficients")
    p.add_argument("--fit-init", default=None,
                   help="Optional initial coefficient guess for curve_fit, comma-separated")
    p.add_argument("--fit-maxfev", type=int, default=200000,
                   help="Maximum function evaluations for curve_fit")
    p.add_argument(
        "--fit-init-max-delta",
        type=float,
        default=None,
        help="Optional absolute coefficient bound around --fit-init for curve_fit: coef in [init-delta, init+delta]",
    )
    p.add_argument("--fit-out", default=None,
                   help="If provided, write fitted parameters and metrics to this JSON file. If omitted, a default name is derived.")
    
    # Emit evenly distributed points on fitted surface
    p.add_argument("--emit-points", type=int, default=0, help="Number of evenly distributed (x,y,z) points to sample on the fitted surface (linear space)")
    p.add_argument("--emit-points-out", default=None, help="CSV file to write the sampled surface points")
    p.add_argument("--emit-points-mode", choices=["axis"], default="axis",
                   help="Sampling domain for --emit-points: 'axis' = rectangle derived from fitted surface axis zero-crossings")
    p.add_argument("--emit-x-values", default=None,
                   help="Comma-separated X values to emit, e.g. '0.001,0.002,0.003'. Used with --emit-y-values.")
    p.add_argument("--emit-y-values", default=None,
                   help="Comma-separated Y values to emit, e.g. '0.01,0.02'. Used with --emit-x-values.")
    p.add_argument("--plot-emitted", action="store_true",
                   help="Overlay the emitted (x,y,z) points onto the 3D plot")
    return p.parse_args()


def _parse_float_list(values: str) -> List[float]:
    text = values.strip()
    if text.startswith("[") and text.endswith("]"):
        text = text[1:-1]
    if not text.strip():
        return []
    out = []
    for token in text.split(","):
        t = token.strip()
        if not t:
            continue
        out.append(float(t))
    return out


def _auto_major_ticks(
    lim: Tuple[float, float],
    scale: float,
    display_step: float,
) -> Tuple[np.ndarray, float]:
    lo, hi = [float(v) for v in lim]
    if not (np.isfinite(lo) and np.isfinite(hi) and hi > lo):
        return np.array([], dtype=float), 1.0
    scale = float(scale) if np.isfinite(scale) and scale > 0 else 1.0
    display_step = float(display_step)
    if not np.isfinite(display_step) or display_step <= 0:
        return np.array([], dtype=float), 1.0
    display_lo = lo * scale
    display_hi = hi * scale
    start = math.ceil((display_lo - 1e-12 * display_step) / display_step) * display_step
    stop = math.floor((display_hi + 1e-12 * display_step) / display_step) * display_step
    if stop < start:
        ticks_display = np.array([display_lo, display_hi], dtype=float)
    else:
        ticks_display = np.arange(start, stop + 0.5 * display_step, display_step, dtype=float)
    return ticks_display / scale, display_step / scale


def _apply_consistent_ticks(ax, xlim, ylim, zlim, sx, sy, sz) -> None:
    for axis_obj, set_ticks, lim, scale, display_step in (
        (ax.xaxis, ax.set_xticks, xlim, sx, 0.5),
        (ax.yaxis, ax.set_yticks, ylim, sy, 1.0),
        (ax.zaxis, ax.set_zticks, zlim, sz, 0.25),
    ):
        ticks, major_step = _auto_major_ticks(lim, scale, display_step)
        if ticks.size:
            set_ticks(ticks)
        if np.isfinite(major_step) and major_step > 0:
            axis_obj.set_minor_locator(MultipleLocator(major_step / 2.0))


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
    # Include boundary points with a tiny tolerance so edge-aligned slices
    # (e.g., b2p=0) are not visually clipped from the rendered surface.
    span = max(float(np.ptp(data_X)), float(np.ptp(data_Y)), 1.0)
    tol = 1e-9 * span
    inside = path.contains_points(pts, radius=tol)
    return inside.reshape(Xgrid.shape)

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


def _plot_domain_from_fit(
    model: str,
    coef: np.ndarray,
    X: np.ndarray,
    Y: np.ndarray,
    fallback_xlim: Tuple[float, float],
    fallback_ylim: Tuple[float, float],
) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """Choose a positive-quadrant plotting domain for the fitted surface."""
    x_domain = fallback_xlim
    y_domain = fallback_ylim
    if model == "quadratic":
        (x_lo, x_hi), (y_lo, y_hi) = _axis_box_from_quadratic_fit(
            coef,
            x_cap=float(fallback_xlim[1]),
            y_cap=float(fallback_ylim[1]),
            nonneg_only=True,
        )
        if np.isfinite(x_lo) and np.isfinite(x_hi) and x_hi > x_lo >= 0.0:
            x_domain = (float(x_lo), float(x_hi))
        if np.isfinite(y_lo) and np.isfinite(y_hi) and y_hi > y_lo >= 0.0:
            y_domain = (float(y_lo), float(y_hi))
    return x_domain, y_domain


def _limits_from_values(v: np.ndarray, pad_frac: float = 0.03) -> Tuple[float, float]:
    """Finite [min, max] with a small symmetric padding for display."""
    arr = np.asarray(v, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return (0.0, 1.0)
    lo = float(np.min(arr))
    hi = float(np.max(arr))
    if hi > lo:
        span = hi - lo
        pad = max(span * float(pad_frac), 1e-12)
        return (lo - pad, hi + pad)
    # Degenerate axis: create a tiny span around the value.
    center = lo
    delta = max(abs(center) * float(pad_frac), 1e-12)
    return (center - delta, center + delta)

def main():
    args = parse_args()
    df = pd.read_csv(args.csv)

    # Filter valid rows
    df = df.copy()
    df = df[df["status"] == "ok"]
    needed = {"b1p", "b2p", "p_th", "basis"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing columns: {missing}")

    # basis filter    
    df = df[df["basis"] == args.basis]

    if df.empty:
        raise SystemExit("No rows to plot after filtering. Check basis/status.")


    x = df[args.x_axis].values  # D2
    y = df[args.y_axis].values  # p_m
    z = df[args.z_axis].values
    X, Y, Z = x, y, z
    user_xlim = _parse_lim_pair(args.xlim, "xlim")
    user_ylim = _parse_lim_pair(args.ylim, "ylim")
    user_zlim = _parse_lim_pair(args.zlim, "zlim")
    xlim = user_xlim or (0.0, float(np.max(X)))
    ylim = user_ylim or (0.0, float(np.max(Y)))
    zlim = user_zlim or (0.0, float(np.max(Z)))

    # Create figure
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Placeholders for optional emitted points
    emitted_pts_xy = None
    emitted_Zs = None
    sampled_df_emitted = None

    # Storage for fit results for later export
    fit_result = None  # dict with keys: coef, r2, rmse

    # Scatter the mean points
    ax.scatter(X, Y, Z, s=30, depthshade=True)

    # For legend bookkeeping
    surf_fit = None
    fit_zero_contour = None
    fit_coef = None
    fit_init = None
    fit_bounds = None
    if args.fit_model != "none":
        # Surface fit in the current plotting space
        space = "linear space"

        if args.fit_init is not None:
            fit_init = np.asarray(_parse_float_list(args.fit_init), dtype=float)
        if args.fit_init_max_delta is not None:
            if fit_init is None:
                raise SystemExit("--fit-init-max-delta requires --fit-init")
            if args.fit_method != "curve_fit":
                raise SystemExit("--fit-init-max-delta is only used with --fit-method curve_fit")
            delta = float(args.fit_init_max_delta)
            if delta <= 0:
                raise SystemExit("--fit-init-max-delta must be > 0")
            fit_bounds = (fit_init - delta, fit_init + delta)
        fit_coef, fit_r2, fit_rmse = _fit_surface(
            args.fit_model, X, Y, Z,
            method=args.fit_method,
            p0=fit_init,
            fit_bounds=fit_bounds,
            maxfev=args.fit_maxfev,
        )
        fit_result = {"coef": [float(v) for v in fit_coef], "r2": float(fit_r2), "rmse": float(fit_rmse)}
        plot_xlim, plot_ylim = _plot_domain_from_fit(
            args.fit_model,
            fit_coef,
            X,
            Y,
            fallback_xlim=xlim,
            fallback_ylim=ylim,
        )
        nx = ny = 360
        gx = np.linspace(plot_xlim[0], plot_xlim[1], nx)
        gy = np.linspace(plot_ylim[0], plot_ylim[1], ny)
        GX, GY = np.meshgrid(gx, gy)
        FZ = _eval_surface(args.fit_model, fit_coef, GX, GY)
        if user_xlim is None:
            xlim = (0.0, float(max(plot_xlim[1], np.max(X)*1.05)))
        if user_ylim is None:
            ylim = (0.0, float(max(plot_ylim[1], np.max(Y)*1.05)))
        if user_zlim is None:
            fit_z = FZ[np.isfinite(FZ) & (FZ >= 0.0)]
            if fit_z.size:
                zlim = (0.0, float(max(np.max(fit_z), np.max(Z)*1.05)))
            else:
                zlim = (0.0, float(np.max(Z)))
        mask_lims = (
            (GX >= xlim[0]) & (GX <= xlim[1]) &
            (GY >= ylim[0]) & (GY <= ylim[1]) &
            (FZ >= 0.0) &
            np.isfinite(FZ)
        )
        mask = mask_lims
        FZ_plot = np.where(mask, FZ, np.nan)
        valid_cells = (
            np.isfinite(FZ_plot[:-1, :-1]) &
            np.isfinite(FZ_plot[1:, :-1]) &
            np.isfinite(FZ_plot[:-1, 1:]) &
            np.isfinite(FZ_plot[1:, 1:])
        )
        if np.any(valid_cells):
            surf_fit = ax.plot_surface(
                GX,
                GY,
                FZ_plot,
                alpha=0.4,
                antialiased=True,
                linewidth=0,
            )
        else:
            print("Warning: fitted surface has <3 in-range points after x/y/z clipping; skipping fit surface plot.")
        print(f"{args.fit_model.capitalize()} surface fit in", space)
        print(f"fit method = {args.fit_method}")
        print(f"coefficients ({len(fit_coef)} terms) =", fit_coef)
        print(f"R^2 = {fit_r2:.4f}, RMSE = {fit_rmse:.4g}")
    else:
        if args.fit_init is not None or args.fit_init_max_delta is not None:
            raise SystemExit("--fit-init/--fit-init-max-delta cannot be used with --fit-model none")

    # Optionally write fitted parameters/metrics to JSON
    summary = {
        "basis": str(args.basis),
        "axes": {"X": str(args.x_axis), "Y": str(args.y_axis), "Z": str(args.z_axis)},
        "data_summary": {
            "n_points": int(len(Z)),
            "x_min": float(np.min(X)), "x_max": float(np.max(X)),
            "y_min": float(np.min(Y)), "y_max": float(np.max(Y)),
            "z_min": float(np.min(Z)), "z_max": float(np.max(Z)),
        },
        "fit_model": str(args.fit_model),
        "fit_method": str(args.fit_method) if args.fit_model != "none" else None,
        "fit_init": None if fit_init is None else [float(v) for v in fit_init.tolist()],
        "fit_init_max_delta": None if args.fit_init_max_delta is None else float(args.fit_init_max_delta),
    }
    if args.fit_model != "none":
        summary[args.fit_model] = fit_result

    out_json = args.fit_out
    if not out_json or str(out_json).strip() == "":
        out_json = f"fit_{args.basis}_{args.x_axis}_{args.y_axis}_{args.z_axis}.json"
    try:
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)
        print(f"Wrote fitted parameters to {out_json}")
    except Exception as e:
        print(f"Warning: failed to write fit JSON to {out_json}: {e}")
    
    
    # --- Precompute emitted surface points so we can plot them ---
    has_manual_emit_lists = (args.emit_x_values is not None) or (args.emit_y_values is not None)
    if has_manual_emit_lists and (args.emit_x_values is None or args.emit_y_values is None):
        raise SystemExit("Both --emit-x-values and --emit-y-values must be provided together.")

    if has_manual_emit_lists or (getattr(args, 'emit_points', 0) and int(getattr(args, 'emit_points', 0)) > 0):
        if fit_coef is None:
            raise SystemExit("--emit-* options require a fitted model. Use --fit-model quadratic.")
        n_pts = int(getattr(args, 'emit_points', 0))
        mode = getattr(args, 'emit_points_mode', 'axis')

        # Work in linear space for fitting and sampling
        Xlin, Ylin, Zlin = x.copy(), y.copy(), z.copy()
        manual_emit = False
        if has_manual_emit_lists:
            x_vals = _parse_float_list(args.emit_x_values)
            y_vals = _parse_float_list(args.emit_y_values)
            if len(x_vals) == 0 or len(y_vals) == 0:
                raise SystemExit("--emit-x-values and --emit-y-values must each contain at least one value.")
            GX, GY = np.meshgrid(np.asarray(x_vals, dtype=float), np.asarray(y_vals, dtype=float))
            cand = np.column_stack([GX.ravel(), GY.ravel()])
            manual_emit = True
            print(f"Using manual emit grid from input lists: {len(x_vals)} x-values, {len(y_vals)} y-values ({len(cand)} points).")
        else:
            qcoef_lin = fit_coef
            (x_lo, x_hi), (y_lo, y_hi) = _axis_box_from_quadratic_fit(
                qcoef_lin,
                nonneg_only=True,
            )
            if not (isinstance(x_lo,float) and isinstance(x_hi,float) and x_lo<x_hi and x_lo>=0 and \
                    isinstance(y_lo,float) and isinstance(y_hi,float) and y_lo<y_hi and y_lo>=0):
                raise SystemExit("axis sampling could not find valid positive x/y zero-crossing bounds.")

            grid_side = max(2, int(np.sqrt(n_pts) + 0.5))
            gx = np.linspace(x_lo, x_hi, grid_side)
            gy = np.linspace(y_lo, y_hi, grid_side)
            GX, GY = np.meshgrid(gx, gy)
            cand = np.column_stack([GX.ravel(), GY.ravel()])

            # For axis-mode, keep only points where fitted quadratic surface stays above zmin.
            zmin = 3e-5
            qz_cand = _eval_quadratic_surface(fit_coef, cand[:, 0], cand[:, 1])
            cand = cand[qz_cand >= zmin]

        if cand.shape[0] > 0:
            # In manual mode, use all requested XY points.
            if manual_emit:
                emitted_pts_xy = cand
            else:
                # Choose approximately evenly spaced subset.
                if cand.shape[0] >= n_pts:
                    idx = np.linspace(0, cand.shape[0]-1, n_pts).astype(int)
                    emitted_pts_xy = cand[idx]
                else:
                    emitted_pts_xy = cand

            # Evaluate Z using the selected fitted model (in linear space)
            emitted_Zs = _eval_surface(args.fit_model, fit_coef, emitted_pts_xy[:,0], emitted_pts_xy[:,1])

            if manual_emit:
                # Keep all user-requested points except non-finite evaluations.
                keep = np.isfinite(emitted_Zs)
            else:
                # Enforce Z >= emit_z_min for auto-generated points.
                zmin = 3e-5
                keep = np.isfinite(emitted_Zs) & (emitted_Zs >= zmin)
            if np.count_nonzero(keep) == 0:
                sampled_df_emitted = None
            else:
                emitted_pts_xy = emitted_pts_xy[keep]
                emitted_Zs = emitted_Zs[keep]
                sampled_df_emitted = pd.DataFrame({"x": emitted_pts_xy[:,0], "y": emitted_pts_xy[:,1], "z": emitted_Zs})
                # Add basis tag to emitted points output
                sampled_df_emitted["basis"] = args.basis

        # If requested, overlay on the current 3D plot
        if getattr(args, 'plot_emitted', False) and emitted_pts_xy is not None and emitted_Zs is not None:
            dispX, dispY, dispZ = emitted_pts_xy[:,0], emitted_pts_xy[:,1], emitted_Zs
            ax.scatter(dispX, dispY, dispZ, s=18, marker="x")
            # When plotting emitted points, fit unspecified axes to emitted ranges.
            # Explicit --xlim/--ylim/--zlim always win.
            if user_xlim is None:
                xlim = _limits_from_values(dispX)
            if user_ylim is None:
                ylim = _limits_from_values(dispY)
            if user_zlim is None:
                zlim = _limits_from_values(dispZ)

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

    # Axis labels: map internal column names to noise-type labels by default
    noise_label_map = {
        "b1p": "2-qubit depolarizing ",
        "b2p": "measurement noise",
        "p_th": "idle noise",
    }

    xlab = args.x_label if args.x_label is not None else noise_label_map.get(args.x_axis, str(args.x_axis))
    ylab = args.y_label if args.y_label is not None else noise_label_map.get(args.y_axis, str(args.y_axis))
    zlab = args.z_label if args.z_label is not None else noise_label_map.get(args.z_axis, str(args.z_axis))

    scale_global = float(args.tick_scale)
    if np.isfinite(scale_global) and scale_global > 0:
        sx = sy = sz = scale_global
    else:
        sx = _auto_tick_scale(*xlim)
        sy = _auto_tick_scale(*ylim)
        sz = _auto_tick_scale(*zlim)

    ax.set_xlabel(f"{xlab}{_tick_scale_suffix(sx)}", fontsize=18, labelpad=12)
    ax.set_ylabel(f"{ylab}{_tick_scale_suffix(sy)}", fontsize=18, labelpad=12)
    ax.set_zlabel(f"{zlab}{_tick_scale_suffix(sz)}", fontsize=18, labelpad=12)
    ax.tick_params(axis='x', which='major', labelsize=14)
    ax.tick_params(axis='y', which='major', labelsize=14)
    ax.tick_params(axis='z', which='major', labelsize=14)
    ax.tick_params(axis='x', which='minor', labelsize=0)
    ax.tick_params(axis='y', which='minor', labelsize=0)
    ax.tick_params(axis='z', which='minor', labelsize=0)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda val, pos: f"{val * sx:g}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda val, pos: f"{val * sy:g}"))
    ax.zaxis.set_major_formatter(FuncFormatter(lambda val, pos: f"{val * sz:g}"))

    # Apply axis limits
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_zlim(*zlim)
    _apply_consistent_ticks(ax, xlim, ylim, zlim, sx, sy, sz)
    ax.view_init(elev=args.elev, azim=args.azim)

    title_basis = ""
    if "basis" in df.columns:
        uniq_b = sorted(df["basis"].unique())
        if args.basis:
            title_basis = f" (basis={args.basis})"
        elif len(uniq_b) == 1:
            title_basis = f" (basis={uniq_b[0]})"
        else:
            title_basis = f" (bases={','.join(map(str, uniq_b))})"

    # ax.set_title(f"Threshold surface: D2 vs p_m vs p_idle{title_basis}")

    # Legend: use proxy artists so it works reliably in 3D
    legend_handles = []
    legend_labels = []

    # # Data points
    # legend_handles.append(Line2D([0], [0], marker='o', linestyle='None', markersize=6, color='black'))
    # legend_labels.append('data points')

    # Fitted surfaces (only if created)
    # if 'surf_fit' in locals() and surf_fit is not None:
    #     legend_handles.append(Patch(alpha=0.4, color='black'))
    #     legend_labels.append(f'{args.fit_model} fit surface')

    # if fit_zero_contour is not None:
    #     legend_handles.append(Line2D([0], [0], color='crimson', linewidth=2.0))
    #     legend_labels.append('fitted z=0 contour')

    # Emitted points (only if requested and available)
    if getattr(args, 'plot_emitted', False) and emitted_pts_xy is not None and emitted_Zs is not None:
        legend_handles.append(Line2D([0], [0], marker='x', linestyle='None', markersize=6, color='black'))
        legend_labels.append('emitted surface points')

    if legend_handles:
        ax.legend(legend_handles, legend_labels, loc='best')

    plt.tight_layout()
    if args.out:
        fig.savefig(args.out, dpi=200)
        print(f"Saved 3D plot to {args.out}")
    else:
        plt.show()

if __name__ == "__main__":
    main()
