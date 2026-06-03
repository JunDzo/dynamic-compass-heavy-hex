#!/usr/bin/env python3
"""
render_rotation.py

Generate an auto-rotating video (MP4) or GIF of the same 3D threshold plot
as plot_threshold_surface_3d.py, by sweeping the azimuth angle.

Examples:
  # MP4 (requires ffmpeg)
  python render_rotation.py --csv approx_thresholds.csv --basis X --fit both \
    --out threshold_rotate.mp4 --fps 25 --frames 150 --elev 25

  # GIF (requires pillow; larger file)
  python render_rotation.py --csv approx_thresholds.csv --basis X --fit linear \
    --out threshold_rotate.gif --fps 20 --frames 120 --elev 25
"""

import argparse
import math
import shutil
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.path import Path
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.animation as animation


# -----------------------------
# Fit helpers (copied from your plot script)
# -----------------------------
def _design_linear(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    X = np.asarray(X).ravel()
    Y = np.asarray(Y).ravel()
    return np.column_stack([np.ones_like(X), X, Y])


def _fit_linear_surface(X: np.ndarray, Y: np.ndarray, Z: np.ndarray) -> Tuple[np.ndarray, float, float]:
    """Return (coeffs, r2, rmse) for Z ≈ a + bX + cY."""
    M = _design_linear(X, Y)
    Z = np.asarray(Z).ravel()
    coef, _, _, _ = np.linalg.lstsq(M, Z, rcond=None)
    if Z.size > 0:
        y_hat = M @ coef
        ss_res = float(np.sum((Z - y_hat) ** 2))
        ss_tot = float(np.sum((Z - np.mean(Z)) ** 2)) if Z.size > 1 else 0.0
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
        rmse = float(np.sqrt(ss_res / max(1, Z.size)))
    else:
        r2, rmse = float("nan"), float("nan")
    return coef, r2, rmse


def _eval_linear_surface(coef: np.ndarray, X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    a, b, c = coef
    return a + b * X + c * Y


def _design_quadratic(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    X = np.asarray(X).ravel()
    Y = np.asarray(Y).ravel()
    return np.column_stack([np.ones_like(X), X, Y, X**2, X * Y, Y**2])


def _fit_quadratic_surface(X: np.ndarray, Y: np.ndarray, Z: np.ndarray) -> Tuple[np.ndarray, float, float]:
    """Return (coeffs, r2, rmse) for Z ≈ a + bX + cY + dX^2 + eXY + fY^2."""
    M = _design_quadratic(X, Y)
    Z = np.asarray(Z).ravel()
    coef, _, _, _ = np.linalg.lstsq(M, Z, rcond=None)
    if Z.size > 0:
        y_hat = M @ coef
        ss_res = float(np.sum((Z - y_hat) ** 2))
        ss_tot = float(np.sum((Z - np.mean(Z)) ** 2)) if Z.size > 1 else 0.0
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
        rmse = float(np.sqrt(ss_res / max(1, Z.size)))
    else:
        r2, rmse = float("nan"), float("nan")
    return coef, r2, rmse


def _eval_quadratic_surface(coef: np.ndarray, X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    a, b, c, d, e, f = coef
    return a + b * X + c * Y + d * X**2 + e * X * Y + f * Y**2


# -----------------------------
# Geometry helpers (convex hull mask)
# -----------------------------
def _convex_hull(points: np.ndarray) -> np.ndarray:
    """Compute 2D convex hull indices (Andrew's monotone chain). Returns indices in CCW order."""
    pts = np.asarray(points, float)
    if pts.shape[0] <= 3:
        return np.arange(pts.shape[0])

    pts_sorted = pts[np.lexsort((pts[:, 1], pts[:, 0]))]

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

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

    idx = []
    used = set()
    for hp in hull:
        for i, pp in enumerate(pts):
            if i in used:
                continue
            if np.allclose(pp, hp):
                idx.append(i)
                used.add(i)
                break
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


# -----------------------------
# Plot building
# -----------------------------
def build_plot(
    csv_path: str,
    basis: str,
    x_axis: str,
    y_axis: str,
    z_axis: str,
    fit: str,
    error_bars: bool,
    figsize=(10, 8),
):
    df = pd.read_csv(csv_path)
    df = df.copy()
    df = df[df["status"] == "ok"]

    needed = {"b1p", "b2p", "p_th", "basis"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing columns: {missing}")

    df = df[df["basis"] == basis]
    if df.empty:
        raise SystemExit("No rows to plot after filtering. Check basis/status.")

    X = df[x_axis].values
    Y = df[y_axis].values
    Z = df[z_axis].values

    if error_bars:
        if ("p_low" not in df.columns) or ("p_high" not in df.columns):
            raise ValueError("--error-bars requested but CSV lacks p_low/p_high")
        Zlo = df["p_low"].values
        Zhi = df["p_high"].values
    else:
        Zlo = Zhi = None

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")

    # scatter
    ax.scatter(X, Y, Z, s=30, depthshade=True)

    # error bars
    if error_bars:
        for xi, yi, z0, z1 in zip(X, Y, Zlo, Zhi):
            ax.plot([xi, xi], [yi, yi], [z0, z1], linewidth=1.5)

    # fits
    surf_lin = None
    surf_quad = None
    if fit != "none":
        nx = ny = 40
        gx = np.linspace(X.min(), X.max(), nx)
        gy = np.linspace(Y.min(), Y.max(), ny)
        GX, GY = np.meshgrid(gx, gy)

        if fit in ("linear", "both"):
            lcoef, lr2, lrmse = _fit_linear_surface(X, Y, Z)
            LZ = _eval_linear_surface(lcoef, GX, GY)
            mask = _mask_grid_to_hull(GX, GY, X, Y)
            Xg = GX[mask]
            Yg = GY[mask]
            Zg = LZ[mask]
            tri = mtri.Triangulation(Xg, Yg)
            surf_lin = ax.plot_trisurf(tri, Zg, alpha=0.25, antialiased=True, linewidth=0)
            print("Linear fit: coef [a,b,c] =", lcoef, f"R^2={lr2:.4f}, RMSE={lrmse:.4g}")

        if fit in ("quadratic", "both"):
            qcoef, qr2, qrmse = _fit_quadratic_surface(X, Y, Z)
            QZ = _eval_quadratic_surface(qcoef, GX, GY)
            mask = _mask_grid_to_hull(GX, GY, X, Y)
            Xg = GX[mask]
            Yg = GY[mask]
            Zg = QZ[mask]
            tri = mtri.Triangulation(Xg, Yg)
            surf_quad = ax.plot_trisurf(tri, Zg, alpha=0.4, antialiased=True, linewidth=0)
            print("Quadratic fit: coef [a,b,c,d,e,f] =", qcoef, f"R^2={qr2:.4f}, RMSE={qrmse:.4g}")

    # labels (match your updated mapping)
    noise_label_map = {
        "b1p": "2 qubit depolarizing (p_2q)",
        "b2p": "measurement noise (p_m)",
        "p_th": "idle noise (p_idle)",
    }
    ax.set_xlabel(noise_label_map.get(x_axis, x_axis))
    ax.set_ylabel(noise_label_map.get(y_axis, y_axis))
    ax.set_zlabel(noise_label_map.get(z_axis, z_axis))

    # limits
    pad = 0.05
    x_min, x_max = float(np.min(X)), float(np.max(X))
    y_min, y_max = float(np.min(Y)), float(np.max(Y))
    z_min, z_max = float(np.min(Z)), float(np.max(Z))
    if error_bars:
        z_min = min(z_min, float(np.min(Zlo)))
        z_max = max(z_max, float(np.max(Zhi)))

    ax.set_xlim(x_min * (1 - pad), x_max * (1 + pad))
    ax.set_ylim(y_min * (1 - pad), y_max * (1 + pad))
    ax.set_zlim(z_min, z_max)

    ax.set_title(f"Threshold surface: p_2q vs p_m vs p_idle (basis={basis})")

    # legend (proxy artists for 3D)
    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="None", markersize=6, color="black"),
    ]
    legend_labels = ["data points"]
    if surf_lin is not None:
        legend_handles.append(Patch(alpha=0.25, color="black"))
        legend_labels.append("linear fit surface")
    if surf_quad is not None:
        legend_handles.append(Patch(alpha=0.4, color="black"))
        legend_labels.append("quadratic fit surface")
    ax.legend(legend_handles, legend_labels, loc="best")

    plt.tight_layout()
    return fig, ax


# -----------------------------
# Animation
# -----------------------------
def save_rotation(
    fig,
    ax,
    out_path: str,
    elev: float,
    azim_start: float,
    azim_end: float,
    n_frames: int,
    fps: int,
    dpi: int,
):
    ext = out_path.lower().split(".")[-1]
    if ext == "mp4":
        if shutil.which("ffmpeg") is None:
            raise SystemExit("ffmpeg not found. Install with: brew install ffmpeg")
        writer = animation.FFMpegWriter(fps=fps)
    elif ext == "gif":
        # PillowWriter is usually available if pillow is installed.
        writer = animation.PillowWriter(fps=fps)
    else:
        raise SystemExit("Unsupported output extension. Use .mp4 or .gif")

    def update(i):
        t = i / max(1, (n_frames - 1))
        azim = azim_start + (azim_end - azim_start) * t
        ax.view_init(elev=elev, azim=azim)
        return []

    ani = animation.FuncAnimation(fig, update, frames=n_frames, interval=1000 / fps, blit=False)
    ani.save(out_path, writer=writer, dpi=dpi)
    print(f"Saved rotation to {out_path}")


def parse_args():
    p = argparse.ArgumentParser(description="Render an auto-rotating 3D threshold plot to MP4/GIF.")
    p.add_argument("--csv", required=True, help="Input CSV (e.g., approx_thresholds.csv)")
    p.add_argument("--basis", choices=["X", "Z"], default="X", help="Basis filter")
    p.add_argument("--x-axis", choices=["b1p", "b2p", "p_th"], default="b1p")
    p.add_argument("--y-axis", choices=["b1p", "b2p", "p_th"], default="b2p")
    p.add_argument("--z-axis", choices=["b1p", "b2p", "p_th"], default="p_th")
    p.add_argument("--fit", choices=["none", "linear", "quadratic", "both"], default="linear")
    p.add_argument("--error-bars", action="store_true", help="Draw vertical error bars if p_low/p_high exist")
    p.add_argument("--out", required=True, help="Output file: .mp4 or .gif")

    p.add_argument("--elev", type=float, default=25.0, help="Elevation angle")
    p.add_argument("--azim-start", type=float, default=0.0, help="Starting azimuth angle")
    p.add_argument("--azim-end", type=float, default=360.0, help="Ending azimuth angle")

    p.add_argument("--frames", type=int, default=150, help="Number of frames")
    p.add_argument("--fps", type=int, default=25, help="Frames per second")
    p.add_argument("--dpi", type=int, default=200, help="Output DPI")

    return p.parse_args()


def main():
    args = parse_args()
    fig, ax = build_plot(
        csv_path=args.csv,
        basis=args.basis,
        x_axis=args.x_axis,
        y_axis=args.y_axis,
        z_axis=args.z_axis,
        fit=args.fit,
        error_bars=args.error_bars,
    )
    # Set initial view before rendering
    ax.view_init(elev=args.elev, azim=args.azim_start)

    save_rotation(
        fig=fig,
        ax=ax,
        out_path=args.out,
        elev=args.elev,
        azim_start=args.azim_start,
        azim_end=args.azim_end,
        n_frames=args.frames,
        fps=args.fps,
        dpi=args.dpi,
    )
    plt.close(fig)


if __name__ == "__main__":
    main()