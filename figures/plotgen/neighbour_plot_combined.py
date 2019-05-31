"""
Uses data files:

    + neighbour_analysis_ratio_gas.npy
    + neighbour_analysis_ratio_star.npy
    + neighbour_analysis_gas_distance.npy
    + neighbour_analysis_star_distance.npy
    + neighbour_analysis_dark_matter_distance.npy
    + neighbour_analysis_gas_fb_agn.npy
    + neighbour_analysis_gas_fb_stellar.npy
    + neighbour_analysis_gas_fb_none.npy

This script makes a plot that is the fraction of mass from given components
as a function of radius, normalised by the virial radius.

You will need _three_ sets of data to do this. Pass the directories within
which these live to the script below in the following order:

    + Full model
    + NoJet
    + NonRadiative

Invoke this in script mode as
    
    python3 neighbour_plot.py <x> <y> <z>
"""

import sys

import matplotlib.pyplot as plt
import numpy as np

plt.style.use("mnras_flatiron")

sims = {"Full": {}, "NoJet": {}, "NonRadiative": {}}

filenames = sys.argv[1:]

for index, path in enumerate(filenames):
    name = list(sims.keys())[index]

    sims[name]["ratio_gas"] = np.load(f"{path}/neighbour_analysis_ratio_gas.npy")
    sims[name]["ratio_star"] = np.load(f"{path}/neighbour_analysis_ratio_star.npy")

    sims[name]["distance_gas"] = np.load(f"{path}/neighbour_analysis_gas_distance.npy")
    sims[name]["distance_star"] = np.load(
        f"{path}/neighbour_analysis_star_distance.npy"
    )
    sims[name]["distance_dark_matter"] = np.load(
        f"{path}/neighbour_analysis_dark_matter_distance.npy"
    )

    sims[name]["fb_agn"] = np.load(f"{path}/neighbour_analysis_gas_fb_agn.npy")
    sims[name]["fb_stellar"] = np.load(f"{path}/neighbour_analysis_gas_fb_stellar.npy")
    sims[name]["fb_none"] = np.load(f"{path}/neighbour_analysis_gas_fb_none.npy")


# Same bins as data gen script
bins_ratio = np.logspace(-4, 5, 256)
# These were given as kpc, but we'd like to plot mpc.
bins_distance = np.linspace(0, 15000, 256) / 1000
find_centers = lambda b: [0.5 * (x + y) for x, y in zip(b[:-1], b[1:])]
centers_ratio = find_centers(bins_ratio)
centers_distance = find_centers(bins_distance)

# Now we can move on to the simple distance histogram plot.

fig, axes = plt.subplots(
    1, 3, figsize=(6.972, 3.0), sharex=False, sharey=True, squeeze=True
)


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


fancy_names = {"gas": "Gas", "star": "Stars", "dark_matter": "Dark Matter"}

axes[0].semilogy()

fb_types = {"agn": "AGN", "stellar": "Stellar", "none": "None"}

fb_styles = {"agn": "dotted", "stellar": "dashed", "none": "solid"}

# Track the minimal value to set the bottom.
min_nonzero = 0

for index, path in enumerate(filenames):
    current_name = list(sims.keys())[index]
    ax = axes[index]

    for name, fancy_name in fb_types.items():
        raw_data = sims[current_name][f"fb_{name}"]
        style = fb_styles[name]

        trimmed = trim_and_norm(raw_data, centers_distance)
        min_nonzero = min(min(trimmed[1][trimmed[1] >= 1e-10]), 1e-10)
        print(trimmed)

        ax.plot(
            *trimmed,
            label=f"$x$ = Gas, $f$ = {fancy_name}",
            color="C0",
            linestyle=style,
        )

    c = 1
    for name, fancy_name in fancy_names.items():
        if name != "gas":
            raw_data = sims[current_name][f"distance_{name}"]

            trimmed = trim_and_norm(raw_data, centers_distance)
            min_nonzero = min(min(trimmed[1][trimmed[1] >= 1e-10]), 1e-10)

            ax.plot(*trimmed, label=f"$x$ = {fancy_name}", color=f"C{c}")
            c += 1

for ax in axes:
    ax.set_xlim(0, 15)
axes[0].set_ylim(1e-9, 1e0)
# axes[0].set_ylim(min_nonzero, axes[0].get_ylim()[1])
axes[0].semilogy()

ticks = [0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0]

for ax in axes[:2]:
    ax.set_xticks(ticks[:-1])
axes[2].set_xticks(ticks)

axes[0].set_ylabel("Normalised fraction of particles in bin")
axes[1].set_xlabel("$r_{x, f}$ [$h^{-1}$ Mpc]")

axes[1].legend()
fig.tight_layout()
fig.subplots_adjust(wspace=0)
fig.savefig("neighbour_analysis_feedback_histogram_combined.pdf")
