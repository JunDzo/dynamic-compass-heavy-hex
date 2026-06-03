#!/usr/bin/env python3
"""
find_approx_thresholds.py
Discover an approximate threshold (crossing) across code distances for each (b1p, b2p) group
from a Sinter stats CSV. Two estimation styles are supported via CLI arg 3:

Output CSV columns (unchanged): basis, b1p, b2p, dists, p_th, status, points_used.

Usage:
    python find_approx_thresholds.py [stats_csv] [out_csv] [sweep_axis] [calc_mode]

Notes:
- The near-threshold region should contain ≥2 p-values and ≥2 distances.
- This script intentionally keeps output field names stable to ease downstream parsing.
"""

import sys
import math
import numpy as np
import pandas as pd
from sinter import read_stats_from_csv_files
from sinter._probability_util import fit_binomial, shot_error_rate_to_piece_error_rate
from statistics import median
from scipy.optimize import minimize_scalar

STATS_CSV = sys.argv[1] if len(sys.argv) > 1 else "out-test/stats.csv"
OUT_CSV   = sys.argv[2] if len(sys.argv) > 2 else "approx_thresholds.csv"
SWEEP_AXIS = sys.argv[3] if len(sys.argv) > 3 else "idle"  # 'idle' or 'M' or 'D2'
CALC_MODE = sys.argv[4] if len(sys.argv) > 4 else "linear"  # 'linear' or 'log'

if CALC_MODE not in {"linear", "log"}:
    raise SystemExit(f"Invalid calc_mode '{CALC_MODE}'. Use 'linear' or 'log'.")

# ------------ load & normalize ------------
stats = read_stats_from_csv_files(STATS_CSV)

rows = []
aggregated: dict[tuple[str, float, int, float, float], dict[str, int | float]] = {}

for s in stats:
    m = s.json_metadata
    if not (
        m.get('b') in {'X','Z'} and               
        m.get('noise') == 'corr_mixing' and
        m.get('sweep') == SWEEP_AXIS and
        m.get('c') == 'heavy_hex' and
        # m.get('d') in {17,23,29} and
        # m.get('r') == 52
        m.get('r') == m.get('d') * 4             
    ):
        continue

    p = m.get("p", None)
    if p is None:
        continue
    d = m.get("d", None)
    b1p = m.get("b1p", None)
    b2p = m.get("b2p", None)
    b = m.get("b", None)

    if d is None or b1p is None or b2p is None:
        continue
    # Sinter semantics: error rate is computed per *kept* shot.
    # kept_shots = shots - discards; best estimate (MLE) = errors / kept_shots
    raw_shots = int(s.shots)
    discards = int(getattr(s, 'discards', 0))
    kept_shots = raw_shots - discards
    if kept_shots <= 0:
        continue
    
    errors = int(s.errors)
    key = (b, float(p), int(d), float(b1p), float(b2p))

    if key not in aggregated:
            aggregated[key] = {
                'b': b,
                'p': float(p),
                'd': int(d),
                'b1p': float(b1p),
                'b2p': float(b2p),
                'shots': 0,
                'errors': 0,
            }
    aggregated[key]['shots'] += kept_shots
    aggregated[key]['errors'] += errors


for item in aggregated.values():
    kept_shots = int(item['shots'])
    errors = int(item['errors'])

    # Match Sinter's plot_error_rate: binomial fit on kept shots
    fit = fit_binomial(
        num_shots=kept_shots,
        num_hits=int(errors),
        max_likelihood_factor=1e3,  # same default you used when plotting
    )

    # Keep per-shot semantics (pieces=1, values=1)
    pieces = 1
    values = 1
    fail_best = shot_error_rate_to_piece_error_rate(
                        fit.best,
                        pieces=pieces,
                        values=values,
                    ) # pyright: ignore[reportCallIssue]
    if errors == 0:
        fail_best = float('nan')

    rows.append({
        "b": str(item["b"]),      # basis 'X' or 'Z'
        "b1p": float(item["b1p"]),  # 2-qubit depolarizing error rate
        "b2p": float(item["b2p"]),  # measurement noise error rate
        "p": float(item["p"]),      # idling error rate
        "d": int(item["d"]),        # code size
        "fail": float(fail_best),  # central estimate, matches Sinter's plotted value
    })

if not rows:
    print("No usable rows found. Check your CSV and json_metadata keys.")
    sys.exit(1)

df = pd.DataFrame(rows).sort_values(["b", "b1p", "b2p", "p", "d"])

def cross_detected(fail_lo: np.ndarray, fail_hi: np.ndarray) -> tuple[bool, tuple[int, int] | None]:
    """Return whether any pair of curves changes relative order between endpoints."""
    n = fail_lo.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            d_lo = fail_lo[i] - fail_lo[j]
            d_hi = fail_hi[i] - fail_hi[j]
            if d_lo * d_hi <= 0:
                return True, (i, j)
    return False, None


def segment_crossings(p0: float, p1: float, y0: list[float], y1: list[float], eps:float=1e-18):
    """Return one linear crossing x for two curves over [p0, p1], else None."""
    d0 = y0[0] - y0[1]
    d1 = y1[0] - y1[1]
    denom = d1 - d0
    t = -d0 / denom if abs(denom) >= eps else None
    if t is not None and 0.0 <= t <= 1.0:
        x = p0 + t * (p1 - p0)
    else:
        x = None
    return x

# Log-log interpolation crossing helper
def segment_crossings_loglog(p0: float, p1: float, y0: list[float], y1: list[float], eps: float = 1e-18):
    """Compute pairwise crossing points using linear interpolation in log-log space.

    This matches the visual interpolation more closely when both axes are shown
    on logarithmic scales. Only pairs with strictly positive x/y values are
    considered.
    """
    x= None
    if p0>0 and p1>0 and p1>0 and p0>0:
        log_p0 = np.log10(p0)
        log_p1 = np.log10(p1)
        d0 = np.log10(y0[0]) - np.log10(y0[1])
        d1 = np.log10(y1[0]) - np.log10(y1[1])
        denom = d1 - d0
        t = -d0 / denom if abs(denom) >= eps else None
        if t is not None and 0.0 <= t <= 1.0:
            log_x = log_p0 + t * (log_p1 - log_p0)
            x = 10 ** log_x
    if x is None:
        x = segment_crossings(p0, p1, y0, y1, eps)
    return x

def densest_linear_cluster(xs: list[float], span: float) -> list[float]:
    """Return the largest subset of sorted xs that fits inside a linear-width window."""
    if not xs:
        return []
    xs_sorted = sorted(xs)
    best_i = 0
    best_j = 0
    j = 0
    n = len(xs_sorted)
    for i in range(n):
        while j < n and xs_sorted[j] - xs_sorted[i] <= span:
            j += 1
        if (j - i) > (best_j - best_i):
            best_i, best_j = i, j
    return xs_sorted[best_i:best_j]


def densest_log_cluster(xs: list[float], log_span: float) -> list[float]:
    """Return the densest cluster in log10-space within a fixed window width."""
    if not xs:
        return []
    xs_sorted = sorted(xs)
    logs = np.log10(np.asarray(xs_sorted, dtype=float))
    best_i = 0
    best_j = 0
    j = 0
    n = len(xs_sorted)
    for i in range(n):
        while j < n and logs[j] - logs[i] <= log_span:
            j += 1
        if (j - i) > (best_j - best_i):
            best_i, best_j = i, j
    return xs_sorted[best_i:best_j]

def find_thresh_cluster(p_grid: np.ndarray, err_grid: np.ndarray,
                        min_hits: int = 3,
                        sat_cut: float = 0.3,
                        max_log_span: float = 0.1,
                        max_linear_span: float = 5e-4,
                        calc_mode: str = "linear",
                        eps: float = 1e-12) -> float | None:
    """Estimate an approximate threshold from pairwise linear crossings.

    The method scans adjacent p-grid segments, gathers crossing witnesses where
    curve-pair order flips, then keeps the densest cluster of witness x-values.
    A threshold is accepted only when enough witnesses exist and at least 70%
    of them lie inside a span of width `max_log_span`.

    Args:
        p_grid: Monotone 1D array of physical error rates.
        err_grid: 2D array shaped (num_p, num_distances) of failure rates.
        min_hits: Minimum witness crossings required before clustering.
        sat_cut: Skip a segment when its lower endpoint has min error > sat_cut.
        max_log_span: Logarithmic-width window used to extract the densest cluster.
        eps: Numeric tolerance for near-degenerate denominator / endpoint touch.

    Returns:
        Median crossing value from the accepted cluster, or None.
    """
    K, N = err_grid.shape
    if K < 2 or N < 2:
        return None
    
    witnesses: list[tuple[float, tuple[int, int]]] = []  # (x, (i,j))
    cross_d: list[tuple[int, int]] = []

    for k in range(K - 1):
        p0, p1 = p_grid[k], p_grid[k + 1]

        if min(err_grid[k, :]) > sat_cut:
            # print(f"Skipping segment [{p0:.3e}, {p1:.3e}] due to high failure rates")
            continue
        
        y0 = err_grid[k, :]
        y1 = err_grid[k + 1, :]

        if not cross_detected(y0, y1)[0]:
            continue

        for i in range(N):
            for j in range(i + 1, N):
                if (i,j) in cross_d:
                    continue
                y0 = [err_grid[k, i], err_grid[k, j]]
                y1 = [err_grid[k + 1, i], err_grid[k + 1, j]]
                d0 = y0[0] - y0[1]
                d1 = y1[0] - y1[1]
                if d0 * d1 < 0:
                    if calc_mode == "log":
                        x = segment_crossings_loglog(p0, p1, y0, y1, eps)
                    else:
                        x = segment_crossings(p0, p1, y0, y1, eps)
                    if x is None:
                        continue
                    # print(f"witness crossing at {x:.3e} between indices {i},{j}")
                    cross_d.append((i, j))
                    witnesses.append((x, (i, j)))
                    continue

                # endpoint touch at upper endpoint p1
                if abs(d1) < eps:
                    if k + 2 < K:
                            d2 = err_grid[k + 2, i] - err_grid[k + 2, j]
                            if d0 * d2 < 0:
                                # print(f"witness endpoint at {p1:.3e} between indices {i},{j}")
                                witnesses.append((p1, (i, j)))
        
    if not witnesses:
        return None

    xs = [w[0] for w in witnesses]
    if len(xs) < min_hits:
        return None

    xs_sorted = sorted(xs)
    if calc_mode == "log":
        cluster = densest_log_cluster(xs_sorted, log_span=max_log_span)
    else:
        cluster = densest_linear_cluster(xs_sorted, span=max_linear_span)

    if len(cluster)/len(xs_sorted) < 0.7:  # Require the best cluster to contain at least 70% of the witnesses
        return None

    threshold_value = float(np.median(cluster))
    if threshold_value < 0:
        return None

    return threshold_value


    
# ------------ main search per basis and (b1p, b2p) ------------
out_rows = []
for basis, df_basis in df.groupby("b", sort=False):
    for (b1p, b2p), grp in df_basis.groupby(["b1p", "b2p"], sort=False):

        grp = grp.sort_values("p", ascending=True)
        pivot = grp.pivot_table(index="p", columns="d", values="fail", aggfunc="first")
        pivot = pivot.dropna(axis=0, how="any") \
             .sort_index(axis=0) \
             .sort_index(axis=1)
                
        if pivot.shape[0] < 2 or pivot.shape[1] < 2:
            out_rows.append({"basis": basis, "b1p": b1p, "b2p": b2p, "status": "insufficient_data"})
            continue

        p_vals = pivot.index.values.astype(float)
        dists = pivot.columns.values
        err_grid = pivot.values.astype(float)  # shape K x N_d

        rng = find_thresh_cluster(
            p_vals,
            err_grid,
            min_hits=4,
            sat_cut=0.5,
            max_log_span=0.1,
            max_linear_span=0.75e-3,
            calc_mode=CALC_MODE,
            eps=1e-12,
        )
        
        if rng is None:
            status = "no_approx_threshold"
            p_th = None
        else:
            p_th = rng
            status = "ok"

        
        out_rows.append({
            "basis": basis,
            "b1p": b1p,
            "b2p": b2p,
            "dists": ",".join(map(str, dists)),
            "p_th": p_th,
            "status": status,
            "points_used": len(p_vals),
        })

out = pd.DataFrame(out_rows).sort_values(["basis", "b1p", "b2p"])
out.to_csv(OUT_CSV, index=False)
# print(f"Wrote {OUT_CSV} with {len(out)} rows.")
# print(out.head(12).to_string(index=False))
