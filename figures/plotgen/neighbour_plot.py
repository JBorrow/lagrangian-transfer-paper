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

Invoke this in script mode as
    
    python3 neighbour_plot.py
"""

import matplotlib.pyplot as plt
import numpy as np

plt.style.use("mnras_flatiron")


ratio_gas = np.load("neighbour_analysis_ratio_gas.npy")
ratio_star = np.load("neighbour_analysis_ratio_star.npy")

distance_gas = np.load("neighbour_analysis_gas_distance.npy")
distance_star = np.load("neighbour_analysis_star_distance.npy")
distance_dark_matter = np.load("neighbour_analysis_dark_matter_distance.npy")

fb_agn = np.load("neighbour_analysis_gas_fb_agn.npy")
fb_stellar = np.load("neighbour_analysis_gas_fb_stellar.npy")
fb_none = np.load("neighbour_analysis_gas_fb_none.npy")


# Same bins as data gen script
bins_ratio = np.logspace(-4, 5, 256)
# These were given as kpc, but we'd like to plot mpc.
bins_distance = np.linspace(0, 15000, 256) / 1000.0
find_centers = lambda b: [0.5 * (x + y) for x, y in zip(b[:-1], b[1:])]
centers_ratio = find_centers(bins_ratio)
centers_distance = find_centers(bins_distance)

# Make the first histogram plot.

fig, ax = plt.subplots()

ax.loglog()

ax.fill_between(
    centers_ratio, np.ones_like(ratio_gas), ratio_gas, alpha=0.8, label="$x$ = Gas"
)
ax.fill_between(
    centers_ratio, np.ones_like(ratio_star), ratio_star, alpha=0.8, label="$x$ = Stars"
)

old_ylim = ax.get_ylim()
# Plot central position to guide the eye
ax.plot([1, 1], [1e-5, 1e8], color="white", ls="dashed")
ax.set_ylim(1, old_ylim[1])

ax.set_xlabel("$r_{x, f}/r_{DM, f}$")
ax.set_ylabel("Abundance of particles in box")

ax.legend()

fig.tight_layout()
fig.savefig("neighbour_analysis_ratio_histogram.pdf")

# Now we can move on to the simple distance histogram plot.

fig, ax = plt.subplots()


def trim_and_norm(data, bin_centers):
    """
    Trims out any bins are after the first zero bin, as well
    as normalising the whole thing (after trimming, of course).
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

    new_data, new_centers = [], []

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

    new_data = new_data / sum(new_data)

    return np.array(new_centers), np.array(new_data)


ax.semilogy()

fancy_names = {"gas": "Gas", "star": "Stars", "dark_matter": "Dark Matter"}

# Track the minimal value to set the bottom.
min_nonzero = 0

for name, fancy_name in fancy_names.items():
    raw_data = locals()[f"distance_{name}"]

    trimmed = trim_and_norm(raw_data, centers_distance)
    min_nonzero = min(min(trimmed[1][trimmed[1] != 0]), 0)

    ax.plot(*trimmed, label=f"$x$ = {fancy_name}")

ax.set_xlim(0, 15)
ax.set_ylim(min_nonzero, ax.get_ylim()[1])

ax.set_ylabel("Normalised fraction of particles in bin")
ax.set_xlabel("$r_{x, f}$ (Mpc/$h$)")

ax.legend()
fig.tight_layout()
fig.savefig("neighbour_analysis_simple_histogram.pdf")


### Now we can make the one that uses the binned data from feedback.

fig, ax = plt.subplots()

ax.semilogy()

fb_types = {"agn": "AGN", "stellar": "Stellar", "none": "No Feedback"}

fb_styles = {"agn": "dotted", "stellar": "dashed", "none": "solid"}

# Track the minimal value to set the bottom.
min_nonzero = 0

for name, fancy_name in fb_types.items():
    raw_data = locals()[f"fb_{name}"]
    style = fb_styles[name]

    trimmed = trim_and_norm(raw_data, centers_distance)
    min_nonzero = min(min(trimmed[1][trimmed[1] != 0]), 0)

    ax.plot(
        *trimmed, label=f"$x$ = Gas, $f$ = {fancy_name}", color="C0", linestyle=style
    )

for name, fancy_name in fancy_names.items():
    if name != gas:
        raw_data = locals()[f"distance_{name}"]

        trimmed = trim_and_norm(raw_data, centers_distance)
        min_nonzero = min(min(trimmed[1][trimmed[1] != 0]), 0)

        ax.plot(*trimmed, label=f"$x$ = {fancy_name}")

ax.set_xlim(0, 15)
ax.set_ylim(min_nonzero, ax.get_ylim()[1])
ax.semilogy()

ax.set_ylabel("Normalised fraction of particles in bin")
ax.set_xlabel("$r_{x, f}$ (Mpc/$h$)")

ax.legend()
fig.tight_layout()
fig.savefig("neighbour_analysis_feedback_histogram.pdf")
