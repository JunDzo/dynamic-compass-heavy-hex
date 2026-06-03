#!/usr/bin/env python
import matplotlib.pyplot as plt
from sinter import plot_error_rate ,read_stats_from_csv_files
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from statistics import median
from sinter._probability_util import fit_binomial, shot_error_rate_to_piece_error_rate
import json
from pathlib import Path

def load_err_rates(csv_path, meta_ok):
    """Load CSV and aggregate raw counts per (p, d) before computing logical error rates."""
    df = pd.read_csv(csv_path, sep=',', skipinitialspace=True)

    aggregated: dict[tuple[float, int], dict[str, int | float]] = {}

    for _, row in df.iterrows():
        m = json.loads(row['json_metadata'])
        if not meta_ok(m):
            continue

        p = m.get('p')
        d = m.get('d')
        if p is None or d is None:
            continue

        raw_shots = int(row['shots'])
        discards = int(row.get('discards', 0))
        kept_shots = raw_shots - discards
        if kept_shots <= 0:
            continue

        errors = int(row['errors'])
        key = (float(p), int(d))

        if key not in aggregated:
            aggregated[key] = {
                'p': float(p),
                'd': int(d),
                'shots': 0,
                'errors': 0,
            }

        aggregated[key]['shots'] += kept_shots
        aggregated[key]['errors'] += errors

    data = []
    for item in aggregated.values():
        kept_shots = int(item['shots'])
        errors = int(item['errors'])

        fit = fit_binomial(
            num_shots=kept_shots,
            num_hits=errors,
            max_likelihood_factor=1e3,
        )

        pieces = 1
        values = 1
        fail_best = shot_error_rate_to_piece_error_rate(
                        fit.best,
                        pieces=pieces,
                        values=values,
                    ) # pyright: ignore[reportCallIssue]
        if errors == 0:
            fail_best = float('nan')

        data.append({
            'p': item['p'],
            'd': item['d'],
            'error_rate': fail_best,
            'shots': kept_shots,
            'errors': errors,
        })

    return pd.DataFrame(data)

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

def denset_linear_cluster(xs: list[float], span: float) -> list[float]:
    """Return the densest cluster in linear space within a fixed window width."""
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

def segment_crossings(p0: float, p1: float, y0: list[float], y1: list[float], eps:float=1e-18):
    """Compute linear interpolation crossing points."""
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

# Example usage for debugging:
# xs_linear = segment_crossings(p0, p1, y0, y1)
# xs_loglog = segment_crossings_loglog(p0, p1, y0, y1)
# print("linear crossings:", xs_linear)
# print("log-log crossings:", xs_loglog)

def find_thresh_cluster(p_grid: np.ndarray, err_grid: np.ndarray,
                        min_hits: int = 3,
                        sat_cut: float = 0.3,
                        max_span: float = 0.1,
                        eps: float = 1e-12) -> dict | None:
    """
    Analyze failure rate matrix to identify a coherent threshold region.

    This routine scans adjacent p-segments and all pairs of curves (indexed by
    code distance) looking for valid crossings.  A crossing is either a strict
    sign change inside the segment or a confirmed endpoint touch that flips
    sign when extended one step further.  Each valid witness point is recorded
    along with the pair of curve indices that generated it.

    After collecting all witnesses we decide whether they form a tight cluster
    suitable to call a threshold.  The decision criteria are:

      * at least ``min_hits`` valid points exist
      * the log10 spread of the points is no larger than ``max_span``

    If the criteria are met the function returns a dictionary containing the
    median threshold, the raw crossing points list, count, log spread, witness
    interval, and a ``status`` string.  Otherwise ``None`` is returned.

    Args:
        p_grid (np.ndarray): monotonically increasing noise strengths of
            shape ``(K,)``.
        err_grid (np.ndarray): failure rates of shape ``(K, N)`` where the
            second axis indexes distinct code distances.
        min_hits (int): minimum number of valid witnesses required.
        max_span (float): allowable log10 span of witness points.
        eps (float): absolute tolerance for endpoint touching detection.

    Returns:
        Optional[dict]: see description above or ``None`` if no coherent
            threshold region was found.
    """
    K, N = err_grid.shape
    if K < 2 or N < 2:
        return None

    witnesses: list[tuple[float, tuple[int, int]]] = []  # (x, (i,j))
    cross_d: list[tuple[int, int]] = []

    for k in range(K - 1):
        p0 = p_grid[k]
        p1 = p_grid[k + 1]
        if min(err_grid[k, :]) > sat_cut:
            print(f"Skipping segment [{p0:.3e}, {p1:.3e}] due to high failure rates")
            continue
        for i in range(N):
            for j in range(i + 1, N):
                if (i,j) in cross_d:
                    continue
                y0 = [err_grid[k, i], err_grid[k, j]]
                y1 = [err_grid[k + 1, i], err_grid[k + 1, j]]
                d0 = y0[0] - y0[1]
                d1 = y1[0] - y1[1]

                # strict crossing inside interval
                if d0 * d1 < 0:
                    # x = segment_crossings(p0, p1, y0, y1)
                    x= segment_crossings_loglog(p0, p1, y0, y1)
                    if x is None:
                        continue
                    print(f"witness crossing at {x:.3e} between indices {i},{j}")
                    cross_d.append((i, j))
                    witnesses.append((x, (i, j)))
                    continue

                # endpoint touch at upper endpoint p1
                if abs(d1) < eps:
                    if k + 2 < K:
                            d2 = err_grid[k + 2, i] - err_grid[k + 2, j]
                            if d0 * d2 < 0:
                                print(f"witness endpoint at {p1:.3e} between indices {i},{j}")
                                witnesses.append((p1, (i, j)))
    if not witnesses:
        return None

    xs = [w[0] for w in witnesses]
    if len(xs) < min_hits:
        return None

    xs_sorted = sorted(xs)
    # Find the densest cluster in log10-space within the specified width
    cluster = densest_log_cluster(xs_sorted, log_span=max_span)
    # cluster = denset_linear_cluster(xs_sorted, span=max_span)

    if len(cluster)/len(xs_sorted) < 0.7:  # Require the best cluster to contain at least 70% of the witnesses
        return None


    spread_log10 = float(np.log10(max(cluster)) - np.log10(min(cluster))) if len(cluster) > 1 else 0.0

    threshold_value = float(np.median(cluster))
    witness_interval = (float(min(cluster)), float(max(cluster)))

    return {
        "threshold": threshold_value,
        "crossings": xs,
        "cluster": cluster,
        "n_crossings": len(xs),
        "cluster_size": len(cluster),
        "spread_log10": spread_log10,
        "witness_interval": witness_interval,
        "status": "threshold",
    }

def find_thresh_from_csv(csv_path, meta_ok):
    """Compute threshold using linear crossing method from CSV data."""
    df = load_err_rates(csv_path, meta_ok)
    if df.empty:
        return None
    
    # Pivot
    tbl = df.pivot_table(index='p', columns='d', values='error_rate', aggfunc='first')
    tbl = tbl.dropna(axis=0, how='any').sort_index(axis=0).sort_index(axis=1)
    
    if tbl.shape[0] < 2 or tbl.shape[1] < 2:
        return None
    
    print("Pivot table:")
    print(tbl)
    print("Aggregated counts by (p, d):")
    print(df[['p', 'd', 'shots', 'errors']].sort_values(['p', 'd']))
    
    p_grid = tbl.index.values.astype(float)
    err_grid = tbl.values.astype(float)
    
    result = find_thresh_cluster(p_grid, err_grid, min_hits=3, max_span=0.05, eps=1e-12)
    # result is either None or a dictionary with threshold info
    return result

def draw_thresh_seg(ax, thresh_info, y_min=None, y_max=None):
    """Draw a vertical threshold segment at x = threshold."""
    if not thresh_info or "threshold" not in thresh_info:
        return
    x = float(thresh_info["threshold"])
    if y_min is None or y_max is None:
        y_min, y_max = ax.get_ylim()
    ax.vlines(
        x=x,
        ymin=y_min,
        ymax=y_max,
        colors="#9F87FD",
        linestyles='--',
        linewidth=5.0,
        # label=f"p_th={x:.2e}",
    )


# ########################### PLOT THRESHOLD ########################### #




def group_func(stat):
    d = stat.json_metadata['d']
    return {
        'label': f'd={d}',
        'sort': d  # Tells sinter how to order the legend
    }

def make_filter_func(current_basis):
    def filter_func(stat):
        meta = stat.json_metadata
        return meta.get('b') == current_basis and meta.get('d') in [9, 11, 13, 15, 17]

    return filter_func

# and meta.get('b1p')==0.0 and meta.get('b2p')==0.0
# meta.get('d') in [7, 9, 11, 13, 15, 17, 31, 33, 35] and
# 7, 9, 11, 13, 15, 17,
# 31, 33, 35

def compute_threshold_info(current_basis):
    threshold_info = find_thresh_from_csv(
        csv_path,
        lambda m: m.get('b') == current_basis and m.get('d') in [9, 11, 13, 15, 17],
    )
    if threshold_info is not None:
        print(f"[{current_basis}] Threshold estimate: {threshold_info['threshold']:.2e}")
        print(f"[{current_basis}] Crossings ({threshold_info['n_crossings']}): {threshold_info['crossings']}")
        print(f"[{current_basis}] Spread log10: {threshold_info['spread_log10']:.3f}")
        print(f"[{current_basis}] Witness interval: {threshold_info['witness_interval']}")
    else:
        print(f"[{current_basis}] Could not estimate threshold.")
    return threshold_info


def make_x_ticks(threshold_info=None, show_xlabel=True):
    # ticks = [6e-4, 1e-3, 2e-3, 3e-3]
    ticks = [1e-3, 2e-3, 3e-3]

    if threshold_info and "threshold" in threshold_info:
        threshold = float(threshold_info["threshold"])
        if 6e-4 <= threshold <= 3e-3 and not any(np.isclose(threshold, tick, rtol=0, atol=5e-6) for tick in ticks):
            ticks.append(threshold)

    ticks = sorted(ticks)
    if show_xlabel:
        labels = [f"{tick * 1e3:.2f}".rstrip("0").rstrip(".") for tick in ticks]
    else:
        labels = []
        for tick in ticks:
            is_threshold = (
                threshold_info
                and "threshold" in threshold_info
                and np.isclose(float(threshold_info["threshold"]), tick, rtol=0, atol=5e-6)
            )
            labels.append(f"{tick * 1e3:.2f}".rstrip("0").rstrip(".") if is_threshold else "")
    return ticks, labels


def style_axis(ax, show_ylabel, show_xlabel, threshold_info=None):
    ax.set_box_aspect(1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_xlim(6e-4, 3e-3)
    ax.set_ylim(1e-3, 1)
    ax.set_xlim(1e-3, 2.1e-3)

    # ax.set_xlabel("Noise Strength (p)")
    # ax.xaxis.label.set_size(40)
    ax.tick_params(axis='both', which='major', labelsize=30, length=10, width=2)
    ax.tick_params(axis='both', which='minor', labelsize=30, length=7, width=1)

    if show_ylabel:
        ax.set_ylabel("Logical Error Rate (per shot)")
        ax.yaxis.label.set_size(30)
    else:
        ax.set_ylabel("")
        ax.tick_params(axis='y', which='both', left=True, labelleft=False)
        ax.spines['left'].set_visible(True)

    x_ticks, x_ticklabels = make_x_ticks(None, show_xlabel=show_xlabel)
    ax.xaxis.set_major_locator(ticker.FixedLocator(x_ticks))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter(x_ticklabels))
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.xaxis.set_minor_locator(
        ticker.LogLocator(base=10.0, subs=np.arange(2, 40) * 0.1, numticks=100)
        )   

    if show_xlabel:
        ax.tick_params(axis='x', which='both', bottom=True, labelbottom=True)
        ax.set_xlabel(r"Noise Strength ($\times 10^{-3}$)")
        ax.xaxis.label.set_size(30)
    else:
        ax.set_xlabel("")
        ax.tick_params(axis='x', which='major', bottom=True, labelbottom=True)
        ax.tick_params(axis='x', which='minor', bottom=True, labelbottom=False)
        ax.spines['bottom'].set_visible(True)


    ax.legend(loc='lower right', fontsize=30)
    ax.set_axisbelow(True)
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, subs=(1.0,), numticks=11))


def plot_basis(ax, current_basis, show_ylabel, show_xlabel):
    threshold_info = compute_threshold_info(current_basis)
    plot_error_rate(
        ax=ax,
        stats=stats,
        x_func=x_func,
        group_func=group_func,
        failure_units_per_shot_func=lambda stat: 1,
        failure_values_func=lambda stat: 1,
        filter_func=make_filter_func(current_basis),
        plot_args_func=lambda index, group_key, group_stats: {
            "color": colors[index % len(colors)],
            "markersize": 6,
        },
        highlight_max_likelihood_factor=1e3,
    )
    if threshold_info is not None:
        draw_thresh_seg(ax, threshold_info, y_min=1e-4, y_max=1)
    style_axis(ax, show_ylabel=show_ylabel, show_xlabel=show_xlabel, threshold_info=threshold_info)


def make_figure(current_basis,show_xlabel):
    if current_basis in ('XZ', 'ZX'):
        ordered_bases = list(current_basis)
        fig, axes = plt.subplots(1, 2, figsize=(20, 10), sharey=True)
        for index, single_basis in enumerate(ordered_bases):
            plot_basis(axes[index], single_basis, show_ylabel=index == 0, show_xlabel=show_xlabel)
        # fig.subplots_adjust(top=0.98)
        return fig

    fig, ax = plt.subplots(figsize=(11, 10))
    plot_basis(ax, current_basis, show_ylabel=True, show_xlabel=show_xlabel)
    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.12, top=0.96)
    return fig

plt.rcParams.update({'font.size': 18})
basis = 'Z'
colors = ["#83AAC4", "#B07AA1", "#FDCE87", "#EFABAA", "#A4A885"]
# colors = ["#155153", "#1C6C73", "#59999E", "#BBAC70", "#9EB3BA"]
# colors = ["#476883", "#A4A885", "#FDCE87", "#B27975", "#B07AA1"]

# 1) Load stats from CSV
csv_path = "figures/manuscript/no-reset-finer/stats.csv"
show_xlabel = True
stats = read_stats_from_csv_files(csv_path)

# 2) Define x_func, group_func, filter_func, etc.
x_func = lambda stat: stat.json_metadata['p']

fig = make_figure(basis, show_xlabel=show_xlabel)
output_path = Path(csv_path).resolve().parent / f"threshold_{basis}_S.pdf"
fig.savefig(output_path)
print(f"Wrote {output_path}")
plt.show()
