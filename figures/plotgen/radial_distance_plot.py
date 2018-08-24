"""
Uses data files:

    + radial_distance_analysis_gas.npy
    + radial_distance_analysis_star.npy
    + radial_distance_analysis_gas_stderr.npy
    + radial_distance_analysis_star_stderr.npy
    + radial_distance_analysis_radial_bins.npy

This script makes a plot that is the fraction of mass from given components
as a function of radius, normalised by the virial radius.

Invoke this in script mode as
    
    python3 radial_distance_plotting.py
"""

import matplotlib.pyplot as plt
import numpy as np

plt.style.use("mnras_durham")

output = np.load("radial_distance_analysis_gas.npy")
output_star = np.load("radial_distance_analysis_star.npy")
output_std = np.load("radial_distance_analysis_gas_stderr.npy")
output_star_std = np.load("radial_distance_analysis_star_stderr.npy")
bin_centers = np.load("radial_distance_analysis_radial_bins.npy")


for data, errors, name, color in zip(
    output, output_std, ["Own LR", "Other LR", "Outside LR"], ["C0", "C1", "C2"]
):
    plt.fill_between(
        bin_centers, data - errors, data + errors, alpha=0.2, color=color, lw=0
    )
    plt.plot(bin_centers, data, label=name, color=color)

for data, errors, name, color in zip(
    output_star,
    output_star_std,
    ["Own LR", "Other LR", "Outside LR"],
    ["C0", "C1", "C2"],
):
    plt.fill_between(
        bin_centers, data - errors, data + errors, alpha=0.2, color=color, lw=0
    )
    plt.plot(bin_centers, data, color=color, ls="dashed")

plt.legend()
plt.xlabel(r"$r/r_{\rm vir}$")
plt.ylabel("Fraction of mass within bin")
plt.tight_layout()
plt.ylim(0, 0.9)
plt.xlim(0, 1)
plt.savefig("radial_distance_plot.pdf")
