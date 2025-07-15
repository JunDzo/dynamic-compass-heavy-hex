#!/usr/bin/env python
import matplotlib.pyplot as plt
from sinter import plot_error_rate ,read_stats_from_csv_files
import matplotlib.ticker as ticker

# ########################### PLOT ERROR RATE X ########################### #
# Update the global font size.
# plt.rcParams.update({'font.size': 18})

# # Define a custom plot arguments function.
# def custom_plot_args(index, group_key, group_stats):
#     # Add additional plotting arguments here if needed.
#     return {}

# csv_path = 'out/stats.csv'

# # Read the TaskStats list.
# stats = read_stats_from_csv_files(csv_path)

# # Print to verify
# for stat in stats:
#     print(f"Strong ID: {stat.strong_id}, Shots: {stat.shots}, Errors: {stat.errors}, Metadata: {stat.json_metadata}")

# # Define helper functions using the pre-parsed metadata.
# x_func = lambda stat: stat.json_metadata['d']

# def group_func(stat):
#     meta = stat.json_metadata

#     if 'heavy_hex' not in meta.get('c'):
#         raise ValueError("Expected 'heavy_hex' in circuit metadata 'c'.")
#     else:
#         return "heavy_hex"

# failure_units_per_shot_func = lambda stat: stat.json_metadata['r']

# def filter_func(stat):
#     meta = stat.json_metadata
#     return (
#         meta.get('p') == 0.0005 and
#         meta.get('noise') == 'uniform' and
#         meta.get('b') == 'X' and
#         meta.get('r') == meta.get('d') * 4
#         # meta.get('r') == 3
#     )


# # Create a matplotlib figure and axes.
# fig, ax = plt.subplots(figsize=(10, 10))
# ax.set_ylim(1e-10, 1)
# ax.set_xlim(0, 35)

# # Call the sinter plotting function.
# plot_error_rate(
#     ax=ax,
#     stats=stats,
#     x_func=x_func,
#     failure_units_per_shot_func=failure_units_per_shot_func,
#     failure_values_func=lambda stat: 1,   # only "X" 
#     group_func=group_func,
#     filter_func=filter_func,
#     plot_args_func=custom_plot_args,
#     highlight_max_likelihood_factor=1e3,
#     line_fits=None,
#     point_label_func=lambda stat: None,
# )

# # Add plot titles and axis labels manually.
# # ax.set_title("basis=X, rounds=4d, noise=uniform, p=0.001", pad=15)
# ax.set_xlabel("Grid Diameter (d)")
# ax.set_ylabel("Logical Error Rate (per round)")
# # fig.suptitle("Logical Error Rate of Bacon Shor Code per Round vs Grid Diameter",x=0.5,y=0.96)  # Move subtitle slightly higher


# # Log scale
# ax.set_yscale('log')

# # Make sure the grid is behind the data
# ax.set_axisbelow(True)

# # Show major grid lines for both X and Y
# ax.grid(True, which='major', axis='both',
#          color='gray', alpha=0.7, linewidth=0.7, linestyle='-')

# ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, subs=(1.0,), numticks=11))
# ax.yaxis.set_minor_locator(ticker.LogLocator(subs=(2,3,4,5,6,7,8,9), numticks=30))

# # Now enable minor grid lines horizontally
# ax.grid(True, which='minor', axis='y',
#          color='gray', alpha=0.5, linewidth=0.5, linestyle=':')

# ax.legend(loc='upper right', fontsize=18, frameon=True, ncol=1)

# plt.tight_layout()
# plt.show()

# # Save the figure.
# fig.savefig("error_rate_x.pdf")
# print("Wrote file://error_rate_x.pdf")


# ########################### PLOT THRESHOLD ########################### #


plt.rcParams.update({'font.size': 18})

# 1) Load stats from CSV
csv_path = "out_heavyhex/stats.csv"
stats = read_stats_from_csv_files(csv_path)

# 2) Define x_func, group_func, filter_func, etc.
x_func = lambda stat: stat.json_metadata['p']

def group_func(stat):
    d = stat.json_metadata['d']
    return {
        'label': f'd={d}',
        'sort': d  # Tells sinter how to order the legend
    }

def filter_func(stat):
    """Only basis == 'X'. (Similarly, you can require decoder='pymatching', noise='uniform', etc. if needed.)"""
    meta = stat.json_metadata
    return meta.get('b') == 'X'

# 3) Create figure and axis
fig, ax = plt.subplots(figsize=(10, 10))

# 4) Plot with sinter
plot_error_rate(
    ax=ax,
    stats=stats,
    x_func=x_func,
    group_func=group_func,
    failure_units_per_shot_func=lambda stat: 1,  # or set to stat.json_metadata['r'] if desired
    failure_values_func=lambda stat: 1,
    filter_func=filter_func,
    plot_args_func=lambda index, group_key, group_stats: {},
    highlight_max_likelihood_factor=1e3,
)

# 5) Set log 
ax.set_xscale('log')          # Matches the CLI's "[log]Noise Strength"
ax.set_yscale('log')
ax.set_xlim(1e-6, 1e-1)       # Adjust as needed if your p range is in [1e-4, 1e-2], or remove to auto-scale
ax.set_ylim(1e-7, 1)          # Matches --ymin 1e-7

# 6) Titles, labels, and subtitle
# ax.set_title("basis=X, decoder=pymatching, noise=uniform", pad=15)
ax.set_xlabel("Noise Strength (p)")
ax.set_ylabel("Logical Error Rate (per shot)")
# fig.suptitle("Logical Error Rate vs Noise Strength", y=0.96)

# 7) Optional: legend
ax.legend(loc='best')

# 8) Optional: set grid, e.g.
ax.set_axisbelow(True)
ax.grid(True, which='major', axis='both', color='gray', alpha=0.7, linestyle='-')
ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, subs=(1.0,), numticks=11))
ax.yaxis.set_minor_locator(ticker.LogLocator(subs=(2,3,4,5,6,7,8,9), numticks=30))
ax.grid(True, which='minor', axis='both',
         color='gray', alpha=0.5, linewidth=0.5, linestyle=':')

# Show / save
plt.tight_layout()
fig.savefig("threshold.pdf")
plt.show()
