#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from sinter import plot_error_rate, read_stats_from_csv_files

# ########################### PLOT ERROR RATE X ########################### #
# Update the global font size.
plt.rcParams.update({'font.size': 20})

# colors = ["#83AAC4", "#A4A885", "#FDCE87", "#EFABAA", "#B07AA1"]
colors = ["#83AAC4", "#EFABAA", "#B07AA1", "#A4A885", "#FDCE87"]


def custom_plot_args(index, group_key, group_stats):
    return {
        "color": colors[index % len(colors)],
        "markersize": 10,
    }


csv_path = 'figures/manuscript/various-d/stats.csv'
basis = 'XZ'

# Read the TaskStats list.
stats = read_stats_from_csv_files(csv_path)

# Print to verify
for stat in stats:
    print(
        f"Strong ID: {stat.strong_id}, Shots: {stat.shots}, Errors: {stat.errors}, "
        f"Metadata: {stat.json_metadata}"
    )


def x_func(stat):
    return stat.json_metadata['d']


def group_func(stat):
    meta = stat.json_metadata
    c = meta.get('c')
    return {
        'label': f'{c}',
        'sort': c,
    }


def make_filter_func(current_basis):
    def filter_func(stat):
        meta = stat.json_metadata
        return (
            meta.get('p') == 0.001
            and meta.get('noise') == 'uniform'
            and meta.get('b') == current_basis
            and meta.get('r') == meta.get('d') * 4
        )

    return filter_func


def style_axis(ax, show_ylabel):
    # ax.set_box_aspect(1)
    ax.set_xlim(0, 35)
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 1)
    ax.set_axisbelow(True)

    ax.set_xlabel("Grid Diameter (d)")
    ax.xaxis.label.set_size(40)
    ax.tick_params(axis='both', which='major', labelsize=40)
    ax.tick_params(axis='both', which='minor', labelsize=40)

    if show_ylabel:
        ax.set_ylabel("Logical Error Rate (per shot)")
        ax.yaxis.label.set_size(40)
    else:
        ax.set_ylabel("")
        ax.tick_params(axis='y', which='both', left=False, labelleft=False)
        ax.spines['left'].set_visible(True)

    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, subs=(1.0,), numticks=11))
    ax.yaxis.set_minor_locator(ticker.LogLocator(subs=(2, 3, 4, 5, 6, 7, 8, 9), numticks=30))
    ax.legend(loc='lower left', fontsize=35, frameon=True, ncol=1)


def plot_basis(ax, current_basis, show_ylabel):
    plot_error_rate(
        ax=ax,
        stats=stats,
        x_func=x_func,
        failure_units_per_shot_func=lambda stat: 1,
        failure_values_func=lambda stat: 1,
        group_func=group_func,
        filter_func=make_filter_func(current_basis),
        plot_args_func=custom_plot_args,
        highlight_max_likelihood_factor=1e3,
        line_fits=None,
        point_label_func=lambda stat: None,
    )
    style_axis(ax, show_ylabel=show_ylabel)


def make_figure(current_basis):
    if current_basis in ('XZ', 'ZX'):
        ordered_bases = list(current_basis)
        fig, axes = plt.subplots(1, 2, figsize=(20, 10), sharey=True)
        for index, single_basis in enumerate(ordered_bases):
            plot_basis(axes[index], single_basis, show_ylabel=index == 0)
        fig.subplots_adjust(top=0.98)
        return fig

    fig, ax = plt.subplots(figsize=(10, 10))
    plot_basis(ax, current_basis, show_ylabel=True)
    # fig.subplots_adjust(left=0.14, right=0.98, bottom=0.12, top=0.98)
    return fig


fig = make_figure(basis)
plt.show()

# Save the figure.
fig.savefig(f"figures/manuscript/various-d/error_rate_{basis}.pdf")
print(f"Wrote file://error_rate_{basis}.pdf")
