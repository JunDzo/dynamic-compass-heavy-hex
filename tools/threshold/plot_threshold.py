#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from sinter import plot_error_rate, read_stats_from_csv_files

plt.rcParams.update({'font.size': 18})

def plot_threshold(x_func,basis, b1p, b2p, csv_path="data/itaration_results/stats_3.csv", log_x=True, log_y=True):
    """
    Plot logical error rate vs x-axis variable for a fixed (b1p, b2p) slice.

    Parameters
    ----------
    x_func : callable
        Function mapping a sinter Stat -> float for the x-axis (e.g., lambda s: s.json_metadata['p_idle']).
    b1p : float
        Value for D2 (two-qubit depolarizing) used in the filter (meta['b1p']).
    b2p : float
        Value for M (measurement) used in the filter (meta['b2p']).
    csv_path : str, optional
        Path to the CSV file containing stats (default is "data/itaration_results/stats_4.csv").
    log_x : bool, optional
        Whether to use a logarithmic scale for the x-axis (default is True).
    log_y : bool, optional
        Whether to use a logarithmic scale for the y-axis (default is True).

    Returns
    -------
    (fig, ax) : matplotlib Figure and Axes
    """
    # 1) Load stats from CSV
    stats = read_stats_from_csv_files(csv_path)

    # 2) Groups (by code distance)
    def group_func(stat):
        d = stat.json_metadata['d']
        return {'label': f'd={d}', 'sort': d}

    # 3) Filter to a single slice in (D2, M), basis X, corr_mixing, and rounds=4d
    def filter_func(stat):
        m = stat.json_metadata
        return (
            m.get('b') == basis and
            m.get('noise') == 'corr_mixing' and
            m.get('r') == m.get('d') * 4 and
            # m.get('r') == 52 and
            m.get('b1') == 'D2' and
            m.get('b2') == 'M' and
            m.get('b2p') == b2p and
            m.get('b1p') == b1p
        )

    # 4) Make plot
    fig, ax = plt.subplots(figsize=(10, 10))

    plot_error_rate(
        ax=ax,
        stats=stats,
        x_func=x_func,
        group_func=group_func,
        failure_units_per_shot_func=lambda s: 1,
        # failure_units_per_shot_func= lambda s: s.json_metadata['r'],
        failure_values_func=lambda s: 1,
        filter_func=filter_func,
        plot_args_func=lambda idx, key, gstats: {},
        highlight_max_likelihood_factor=1e3,
    )

    # 5) Axes/log/limits
    if log_x:
        ax.set_xscale('log')
        ax.set_xlim(1e-5, 2e-2)
    else:
        ax.set_xscale('linear')
        ax.autoscale(axis='x', tight=True)
    if log_y:
        ax.set_yscale('log')
        ax.set_ylim(1e-5, 1)
    else:
        ax.set_yscale('linear')
        ax.autoscale(axis='y', tight=True)


    ax.set_xlabel("Noise Strength (p)")
    ax.set_ylabel("Logical Error Rate (per shot)")

    ax.legend(loc='best')
    ax.set_axisbelow(True)
    ax.grid(True, which='major', axis='both', color='gray', alpha=0.7, linestyle='-')
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, subs=(1.0,), numticks=11))
    ax.yaxis.set_minor_locator(ticker.LogLocator(subs=(2,3,4,5,6,7,8,9), numticks=30))
    ax.grid(True, which='minor', axis='both', color='gray', alpha=0.5, linewidth=0.5, linestyle=':')

    plt.tight_layout()
    return fig, ax

# Example usage:
if __name__ == "__main__":
    import os

    # Ensure output directory exists
    os.makedirs("data/plots/thresholds", exist_ok=True)

    # Generate plot with specified parameters
    fig, ax = plot_threshold(
        lambda s: s.json_metadata['p'],
        basis='Z',
        b1p=0.00025,
        b2p=0.0,
        csv_path="data/stats/stats_nr_3.csv",
        log_x=True,
        log_y=True
    )

    # Save to organized output location
    output_file = "data/plots/thresholds/threshold_slice.pdf"
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved plot to: {output_file}")

    plt.show()