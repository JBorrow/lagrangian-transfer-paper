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
lambda find_centers b: [0.5 * (x + y) for x, y in zip(b[:-1], b[1:])]
centers_ratio = find_centers(bins_ratio)
centers_distance = find_centers(bins_distance)

# Make the first histogram plot.

fig, ax = plt.subplots()

ax.loglog()

ax.fill_between(centers_ratio, np.ones_like(ratio_gas), ratio_gas, alpha=0.8, label="$x$ = gas")
ax.fill_between(centers_ratio, np.ones_like(ratio_star), ratio_star, alpha=0.8, label="$x$ = stars")

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

def trim_and_norm(data, min_particles=8):
    """
    Trims out any bins that have less than min_particles in them, as well
    as normalising the whole thing (after trimming, of course).
    """
    
    new_data = np.zeros_like(data)
    new_data[data >= min_particles] = data
    new_data /= new_data.sum()

    return new_data

ax.semilogy()

fancy_names = {
    "gas": "gas",
    "star": "stars",
    "dark_matter": "dark matter"
}

for name, fancy_name in fancy_names.items():
    raw_data = locals()[f"distance_{name}"]
    trimmed = trim_and_norm(raw_data)
    ax.plot(centers_distance, trimmed, label=f"$x$ = {fancy_name}")

ax.set_xlim(0, 15)

ax.set_ylabel("Normalised fraction of particles in bin")
ax.set_xlabel("$r_{x, f}$ (Mpc/$h$)")

ax.legend()
fig.tight_layout()
fig.savefig("neighbour_analysis_simple_histogram.pdf")


### Now we can make the one that uses the binned data from feedback.

fig, ax = plt.subplots()

ax.semilogy()

fb_types = {
    "agn": "AGN",
    "stellar": "Stellar",
    "none": "No Feedback"
}

fb_styles = {
    "agn": "dotted",
    "stellar": "dashed",
    "none": "solid"
}

for name, fancy_name in fb_types.items():
    raw_data = locals()[f"fb_{name}"]
    trimmed = trim_and_norm(raw_data)
    style = fb_styles[name]
    ax.plot(centers_distance, trimmed, label=f"$x$ = Gas, $f$ = {fancy_name}", color="C0", linestyle=style)

for name, fancy_name in fancy_names.items():
    if name != gas:
        raw_data = locals()[f"distance_{name}"]
        trimmed = trim_and_norm(raw_data)
        ax.plot(centers_distance, trimmed, label=f"$x$ = {fancy_name}")

ax.set_xlim(0, 15)
ax.semilogy()

ax.set_ylabel("Normalised fraction of particles in bin")
ax.set_xlabel("$r_{x, f}$ (Mpc/$h$)")

ax.legend()
fig.tight_layout()
fig.savefig("neighbour_analysis_feedback_histogram.pdf")
