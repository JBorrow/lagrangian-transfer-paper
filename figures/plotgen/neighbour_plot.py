"""
Uses data files:

    + neighbour_analysis_data_gas.npy
    + neighbour_analysis_data_dark_matter.npy
    + neighbour_analysis_data_star.npy
    + neighbour_analysis_radii_dm_g.npy
    + neighbour_analysis_radii_dm_s.npy
    + neighbour_analysis_radii_gas.npy
    + neighbour_analysis_radii_star.npy

This script makes a plot that is the fraction of mass from given components
as a function of radius, normalised by the virial radius.

Invoke this in script mode as
    
    python3 neighbour_plot.py
"""

import matplotlib.pyplot as plt
import numpy as np

plt.style.use("mnras_flatiron")

data_gas = np.load("neighbour_analysis_data_gas.npy")
data_dark_matter = np.load("neighbour_analysis_data_dark_matter.npy")
data_star = np.load("neighbour_analysis_data_star.npy")
radii_dm_g = np.load("neighbour_analysis_radii_dm_g.npy")
radii_dm_s = np.load("neighbour_analysis_radii_dm_s.npy")
radii_gas = np.load("neighbour_analysis_radii_gas.npy")
radii_star = np.load("neighbour_analysis_radii_star.npy")

# Make the first histogram plot.

fig, ax = plt.subplots()

ax.loglog()

ax.hist(
    radii_gas / radii_dm_g,
    bins=np.logspace(-4, 4, 256),
    range=(0, 5000),
    alpha=0.8,
    label="$x$ = Gas",
)
ax.hist(
    radii_star / radii_dm_s,
    bins=np.logspace(-4, 4, 256),
    range=(0, 5000),
    alpha=0.8,
    label="$x$ = Stars",
)

old_ylim = ax.get_ylim()
# Plot central position to guide the eye
ax.plot([1, 1], [1e-5, 1e8], color="white", ls="dashed")
ax.set_ylim(old_ylim)

ax.set_xlabel("$r_{x, f}/r_{DM, f}")
ax.set_ylabel("Abundance of particles in box")

fig.legend()

fig.tight_layout()
fig.savefig("neighbour_analysis_ratio_histogram.pdf")


# Now we can move on to the simple distance histogram plot.

fig, ax = plt.subplots()

ax.semilogy()

bin_max = max([data_gas.max(), data_star.max(), data_dm.max()])
bins = np.linspace(0, bin_max, 100)

ax.hist(
    data_gas[1] / 1000.0, bins=bins, density=True, histtype="step", label="$x$ = Gas"
)
ax.hist(
    data_star[1] / 1000.0, bins=bins, density=True, histtype="step", label="$x$ = Stars"
)
ax.hist(
    data_dark_matter[1] / 1000.0,
    bins=bins,
    density=True,
    histtype="step",
    label="$x$ = Dark Matter",
)

plt.ylabel("Normalised fraction of particles in bin")
plt.xlabel("$r_{x, f}$ (Mpc)")

fig.legend()
fig.tight_layout()
fig.savefig("neighbour_analysis_simple_histogram.pdf")
