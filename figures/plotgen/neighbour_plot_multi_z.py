"""
Uses data files:

    + neighbour_analysis_gas_distance.npy
    + neighbour_analysis_dark_matter_distance.npy

This script makes a plot that is the fraction of mass from given components
as a function of radius, normalised by the virial radius.

Invoke this in script mode as
    
    python3 neighbour_plot_multi_z.py
"""

import sys

import matplotlib.pyplot as plt
import numpy as np

plt.style.use("mnras_flatiron")

# redshift: model: file path
sim_names = {
    0: {"NoJet": "s50nojetAHF", "Full": "s50j7kAHF"},
    2: {"NoJet": "s50nojetAHF_z2", "Full": "s50j7kAHF_z2"},
}
sims = {0: {"NoJet": {}, "Full": {}}, 2: {"NoJet": {}, "Full": {}}}

colours = {"Gas": "C0", "Stars": "C1", "Dark Matter": "C2"}
linestyles = {0: "solid", 2: "dashed"}

for redshift, simulations in sim_names.items():
    for name, path in simulations.items():
        sims[redshift][name]["Gas"] = np.load(
            f"{path}/neighbour_analysis_gas_distance.npy"
        )
        sims[redshift][name]["Stars"] = np.load(
            f"{path}/neighbour_analysis_star_distance.npy"
        )
        sims[redshift][name]["Dark Matter"] = np.load(
            f"{path}/neighbour_analysis_dark_matter_distance.npy"
        )


# Same bins as data gen script
bins_ratio = np.logspace(-4, 5, 256)
# These were given as kpc, but we'd like to plot mpc.
bins_distance = np.linspace(0, 15000, 256) / 1000
find_centers = lambda b: [0.5 * (x + y) for x, y in zip(b[:-1], b[1:])]
centers_ratio = find_centers(bins_ratio)
centers_distance = find_centers(bins_distance)

# Now we can move on to the simple distance histogram plot.

fig, axes = plt.subplots(
    2, 1, figsize=(3.321, 3.321 * 1.62), sharex=True, sharey=True, squeeze=True
)
axes = axes.flatten()


def trim_and_norm(data, bin_centers):
    """
    Trims out any bins are after the first zero bin, as well
    as normalising the whole thing (after trimming, of course).
    
    Assumes bins are equal widths.
    """

    # Stop when we get two bins next to each other that are zero.
    try:
        first_double_zero = [
            bin_number
            for bin_number, (x, y) in enumerate(zip(data[:-1], data[1:]))
            if x == 0 and y == 0
        ][0]
    except IndexError:
        # There is no double zero!
        first_double_zero = 2 * data.size

    new_data, new_centers, new_widths = [], [], []

    current_width = 0.0

    for bin_number in range(data.size):
        n_parts = data[bin_number]

        if bin_number == first_double_zero:
            # We've reached the end, friends.
            new_data.append(n_parts)
            new_centers.append(bin_centers[bin_number])

            break

        if n_parts == 0:
            # Skip me
            continue

        try:
            n_parts_next_door = data[bin_number + 1]
            last_bin = False
        except IndexError:
            n_parts_next_door = 0
            last_bin = True

        if n_parts_next_door == 0 and not last_bin:
            # Combine forces!
            total_n_parts = n_parts + n_parts_next_door
            center = 0.5 * (bin_centers[bin_number] + bin_centers[bin_number + 1])
        else:
            new_data.append(n_parts)
            new_centers.append(bin_centers[bin_number])

    width = bin_centers[1] - bin_centers[0]
    new_data = new_data / (sum(new_data) * width)

    return np.array(new_centers), np.array(new_data)


axes[0].semilogy()

# So we want to plot the simulations on different axes, showing the z=0 and z=2 results
# for both simultaneously.

legends = [{}, {}]
label_done = {"NoJet": False, "Full": False}

for redshift, simulation in sims.items():
    for (index, axis), legend, (model, data) in zip(
        enumerate(axes), legends, simulation.items()
    ):
        for ptype, distance_data in data.items():
            colour = colours[ptype]
            linestyle = linestyles[redshift]

            trimmed = trim_and_norm(distance_data, centers_distance)

            line = axis.plot(*trimmed, linestyle=linestyle, color=colour)[0]

            # Deal with legends... ==0 for explicity
            if index == 0:
                legend[f"$z={redshift}$"] = line
            elif index == 1:
                legend[ptype] = line

            # Add on simulation caption
            if not label_done[model]:
                axis.text(
                    1.0 - 0.025,
                    0.025,
                    model,
                    transform=axis.transAxes,
                    ha="right",
                    va="bottom",
                )
                label_done[model] = True


axes[0].set_xlim(0, 15)
axes[0].set_ylim(1e-7, 1e1)

ticks = [0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0]

axes[1].set_xticks(ticks)

ticks = [1e-7, 1e-5, 1e-3, 1e-1]
# Fix overlapping ticks
axes[1].set_yticks(ticks)

axes[0].set_ylabel("Probability density distribution", y=-0.0)
axes[1].set_xlabel("Spread metric $S$ [$h^{-1}$ Mpc]")

# Create the two legends.


axes[0].legend(
    legends[0].values(), legends[0].keys(), markerfirst=False, loc="upper right"
)
axes[1].legend(
    legends[1].values(), legends[1].keys(), markerfirst=False, loc="upper right"
)

fig.tight_layout()
fig.subplots_adjust(hspace=0)
fig.savefig("neighbour_analysis_hist_redshift.pdf")
