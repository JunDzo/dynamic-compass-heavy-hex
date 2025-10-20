#!/usr/bin/env python3
"""
find_common_thresholds.py
Discover a common threshold (crossing) across code distances for each (b1p, b2p) group
from a Sinter stats CSV. Two estimation styles are supported via CLI arg 3:

  - "linear"    : scans adjacent p-intervals and checks pairwise order flips, returning a bracket
                  [p_low, p_high] and a median of pairwise linear crossing estimates.
  - "quadratic" : for each distance d, fits Pfail(p) ≈ a_d + b_d p + c_d p^2 near threshold,
                  then finds p_th minimizing the across-distance variance of {a_d + b_d p + c_d p^2}.
                  The CSV schema is unchanged; the point estimate is stored in p_th and p_low/p_high are None.

Output CSV columns (unchanged): basis, b1p, b2p, dists, p_low, p_high, p_th, status, points_used.

Usage:
    python find_common_thresholds.py [stats_csv] [out_csv] [linear|quadratic]

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
OUT_CSV   = sys.argv[2] if len(sys.argv) > 2 else "common_thresholds.csv"
FIT_STYLE = sys.argv[3] if len(sys.argv) > 3 else "linear"  # 'linear' or 'quadratic'

# ------------ load & normalize ------------
stats = read_stats_from_csv_files(STATS_CSV)

rows = []
for s in stats:
    m = s.json_metadata
    if not (
        m.get('b') in {'X','Z'} and               
        m.get('noise') == 'corr_mixing' and
        m.get('sweep') == 'idle' and
        m.get('c') == 'heavy_hex' and
        m.get('r') == m.get('d') * 4             
    ):
        continue
    # x-axis key: prefer 'p' else 'p_idle'
    p = m.get("p", None)
    if p is None:
        p = m.get("p_idle", None)
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
    kept_shots = max(1, raw_shots - discards)

    # Match Sinter's plot_error_rate: binomial fit on kept shots
    fit = fit_binomial(
        num_shots=kept_shots,
        num_hits=int(s.errors),
        max_likelihood_factor=1e3,  # same default you used when plotting
    )

    # Keep per-shot semantics (pieces=1, values=1)
    pieces = 1
    values = 1
    fail_low  = shot_error_rate_to_piece_error_rate(fit.low,  pieces=pieces, values=values)
    fail_best = shot_error_rate_to_piece_error_rate(fit.best, pieces=pieces, values=values)
    fail_high = shot_error_rate_to_piece_error_rate(fit.high, pieces=pieces, values=values)

    rows.append({
        "b": str(b),
        "b1p": float(b1p),
        "b2p": float(b2p),
        "p": float(p),
        "d": int(d),
        "fail": float(fail_best),  # central estimate, matches Sinter's plotted value
        "fail_low": float(fail_low),
        "fail_high": float(fail_high),
    })

if not rows:
    print("No usable rows found. Check your CSV and json_metadata keys.")
    sys.exit(1)

df = pd.DataFrame(rows).sort_values(["b", "b1p", "b2p", "p", "d"])


def _is_saturated_row(y_low: np.ndarray, y_high: np.ndarray, center: float = 0.3, tol: float = 0.05, min_frac: float = 0.75) -> bool:
    """Return True if a single p-row is saturated near 0.5.

    A row is considered saturated when at least `min_frac` of distances have
    failure rate within `tol` of `center` (default: within ±0.05 of 0.5 for ≥75% distances).
    """
    if y_low.size == 0:
        return False
    close = np.abs(y_low - y_high) <= tol
    uncer = np.abs(((y_low + y_high) / 2) - center) <= tol
    close = np.logical_and(close, uncer)
    return (np.count_nonzero(close) / y_low.size) >= min_frac

def cross_detected(fail_lo: np.ndarray, fail_hi: np.ndarray) -> tuple[bool, tuple[int, int] | None]:
    """
    Check for crossing between failure rate bounds.
    """
    n = fail_lo.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            d_lo = fail_lo[i] - fail_lo[j]
            d_hi = fail_hi[i] - fail_hi[j]
            if d_lo * d_hi <= 0:
                return True, (i, j)
    return False, None

def all_pairs_flip_sign_at_endpoints(fail_lo: np.ndarray, fail_hi: np.ndarray) -> bool:
    """
    Return True if, for every pair i<j, the difference (y_i - y_j) changes sign
    between the two endpoints (fail_lo vs fail_hi), or is exactly zero at one endpoint.
    This indicates a guaranteed crossing for that pair within the segment.
    """
    n = fail_lo.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            d_lo = fail_lo[i] - fail_lo[j]
            d_hi = fail_hi[i] - fail_hi[j]
            if d_lo == 0 and d_hi == 0:
                continue
            if d_lo * d_hi > 0:
                return False
    return True


def pairwise_linear_crossings_in_segment(p0: float, p1: float, fail0: np.ndarray, fail1: np.ndarray):
    """
    Compute linear interpolation crossing points p* where fail_i == fail_j within [p0, p1],
    for all pairs i<j, based on values at the segment endpoints.
    Returns a list of p* inside the segment.
    """
    xs = []
    n = fail0.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            d0 = fail0[i] - fail0[j]
            d1 = fail1[i] - fail1[j]
            denom = (d1 - d0)
            if abs(denom) < 1e-18:
                continue
            t = -d0 / denom  # fraction along [p0,p1]
            if 0.0 <= t <= 1.0:
                xs.append(p0 + t * (p1 - p0))
    return xs

def fit_quadratic_per_distance(p_values: np.ndarray, fail_values: np.ndarray):
    """
    Fit a quadratic Pfail(p) ≈ a + b p + c p^2 for a single code distance.

    Args:
        p_values: 1D array of physical error rates p (near threshold, sorted ascending).
        fail_values: 1D array of corresponding logical failure rates.

    Returns:
        (a, b, c) least-squares coefficients.
    """
    M = np.column_stack([np.ones_like(p_values), p_values, p_values**2])
    a, b, c = np.linalg.lstsq(M, fail_values, rcond=None)[0]
    return a, b, c

def estimate_threshold_from_collapse(per_distance_data: dict, p_bounds: tuple):
    """
    Estimate threshold p_th and intersection value B0 from per-distance quadratics.

    Args:
        per_distance_data: dict {distance: (p_values, fail_values)} near threshold.
        p_bounds: (p_min, p_max) bounded search interval for p_th.

    Returns:
        p_th_hat: float threshold estimate minimizing variance across distances.
        B0_hat:  float mean of Pfail at p_th across distances.
        abc:     dict {distance: (a, b, c)} per-distance quadratic coefficients.
    """
    abc = {}
    for dist, (p_values, fail_values) in per_distance_data.items():
        abc[dist] = fit_quadratic_per_distance(np.asarray(p_values), np.asarray(fail_values))

    dists_sorted = sorted(abc.keys())
    a = np.array([abc[d][0] for d in dists_sorted])
    b = np.array([abc[d][1] for d in dists_sorted])
    c = np.array([abc[d][2] for d in dists_sorted])

    def S(p):
        f = a + b*p + c*p*p
        return np.var(f) * len(f)  # proportional to SSE about mean

    res = minimize_scalar(S, bounds=p_bounds, method='bounded')
    p_th = res.x
    B0 = np.mean(a + b*p_th + c*p_th*p_th)
    return p_th, B0, abc

def estimate_threshold_from_quadratics(p_values: np.ndarray, fail_matrix: np.ndarray, p_bounds: tuple | None = None):
    """
    Estimate threshold p_th and intersection value B0 from a matrix of failure rates,
    by first fitting per-distance quadratics Pfail(p) ≈ a + b p + c p^2 and then
    minimizing the across-distance variance of the fitted polynomials evaluated at p.

    Args:
        p_values:   1D array of physical error rates p (monotone grid of length K).
        fail_matrix: 2D array of shape (K, N_dists) with Pfail per distance.
        p_bounds:   Optional (p_min, p_max) search interval. Defaults to hull of p_values.

    Returns:
        p_th: float threshold estimate.
        B0:  float intersection value at p_th.
        abc: dict {distance_index: (a, b, c)} per-distance quadratic coefficients.
    """
    K, N = fail_matrix.shape
    abc = {}
    a = np.empty(N)
    b = np.empty(N)
    c = np.empty(N)
    M = np.column_stack([np.ones_like(p_values), p_values, p_values**2])
    for idx in range(N):
        coeffs, _, _, _ = np.linalg.lstsq(M, fail_matrix[:, idx], rcond=None)
        a[idx], b[idx], c[idx] = coeffs
        abc[idx] = (a[idx], b[idx], c[idx])

    if p_bounds is None:
        p_bounds = (float(p_values.min()), float(p_values.max()))

    def S(p):
        f = a + b*p + c*p*p
        return np.var(f) * len(f)

    res = minimize_scalar(S, bounds=p_bounds, method='bounded')
    p_th = res.x
    B0 = np.mean(a + b*p_th + c*p_th*p_th)
    return p_th, B0, abc

def estimate_linear_crossing_interval(p_values: np.ndarray, fail_matrix: np.ndarray, width_rel_tol: int = 3):
    """
    Linear bracketing of the threshold region using pairwise order flips.

    Given a monotone grid of p-values (length K) and a matrix of shape (K, N_dists),
    the routine scans adjacent (and sometimes two-step) segments to find an interval
    [p_low, p_high] where every pair of distances flips order between endpoints.
    It returns a tight span of pairwise linear crossing estimates when available.

    Args:
        p_values: 1D monotone array of physical error rates p.
        fail_matrix: 2D array (K x N_dists) of logical failure rates.
        width_rel_tol: tolerance controlling the width of crossing intervals.

    Returns:
        (p_low, p_high, p_median) on success, or None if no common bracket is found.
    """
    K, N = fail_matrix.shape
    if K < 2 or N < 2:
        return None
    p_min, p_max = None, None
    for k in range(K - 1):
        p0, p1 = p_values[k], p_values[k + 1]
        y0 = fail_matrix[k, :]
        y1 = fail_matrix[k + 1, :]
        if _is_saturated_row(y0,y1):
            continue
        if not cross_detected(y0, y1)[0]:
            continue
        if all_pairs_flip_sign_at_endpoints(y0, y1):
            xs = pairwise_linear_crossings_in_segment(p0, p1, y0, y1)
            p_min, p_max = min(xs), max(xs)
            return(p_min, p_max)
        
        width_rel_tol = min(K-k-1, width_rel_tol)
        y2 = fail_matrix[k + width_rel_tol, :]
        yn1 = fail_matrix[k - 1, :] if k - 1 >= 0 else y0
        if all_pairs_flip_sign_at_endpoints(yn1, y2):
            xs = pairwise_linear_crossings_in_segment(p0, p1, y0, y1)
            if not xs:
                raise ValueError("Unexpected empty crossing list.")
            p_min, p_max = min(xs), max(xs)
            return(p_min, p_max)
        else:
            continue

    return (p_min, p_max)


    
# ------------ main search per basis and (b1p, b2p) ------------
out_rows = []
for basis, df_basis in df.groupby("b", sort=False):
    for (b1p, b2p), grp in df_basis.groupby(["b1p", "b2p"], sort=False):
        # pivot: rows=p, cols=d
        # 'fail' is the MLE per-shot estimate consistent with Sinter's plotting
        grp = grp.sort_values("p", ascending=True)
        pivot = grp.pivot_table(index="p", columns="d", values="fail", aggfunc="mean")
        pivot = pivot.dropna(axis=0, how="any") \
             .sort_index(axis=0) \
             .sort_index(axis=1)
                
        if pivot.shape[0] < 2 or pivot.shape[1] < 2:
            out_rows.append({"basis": basis, "b1p": b1p, "b2p": b2p, "status": "insufficient_data"})
            continue

        p_vals = pivot.index.values.astype(float)
        dists = pivot.columns.values
        y = pivot.values.astype(float)  # shape K x N_d

        rng = estimate_linear_crossing_interval(p_vals, y, width_rel_tol=4)
        if rng == (None, None) or rng is None:
            status = "no_common_threshold"
            p_low = p_high = None
        else:
            if FIT_STYLE == 'linear':
                p_low, p_high = rng
                status = "ok"
            if FIT_STYLE == 'quadratic':
                p_low, p_high = rng
                p_bounds = (float(p_vals.min()), float(p_vals.max()))
                p_th, B0, abc = estimate_threshold_from_quadratics(p_vals, y, p_bounds=p_bounds)
                status = "ok"
                
                p_median = p_th  # store the point estimate in this field for now

        
        out_rows.append({
            "basis": basis,
            "b1p": b1p,
            "b2p": b2p,
            "dists": ",".join(map(str, dists)),
            "p_low": p_low,
            "p_high": p_high,
            "p_th": p_low,
            "status": status,
            "points_used": len(p_vals),
        })

out = pd.DataFrame(out_rows).sort_values(["basis", "b1p", "b2p"])
out.to_csv(OUT_CSV, index=False)
print(f"Wrote {OUT_CSV} with {len(out)} rows.")
print(out.head(12).to_string(index=False))