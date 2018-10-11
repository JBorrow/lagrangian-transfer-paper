"""
Uses data files:

    + radial_distance_analysis_gas_fb.npy
    + radial_distance_analysis_radial_bins_fb.npy

This script makes a plot that is the fraction of mass from given components
as a function of radius, normalised by the virial radius.

Invoke this in script mode as
    
    python3 feedback_distance_plot.py
"""

import matplotlib.pyplot as plt
import numpy as np

plt.style.use("mnras_flatiron")

output = np.load("radial_distance_analysis_gas_fb.npy")
bins = np.load("radial_distance_analysis_radial_bins_fb.npy")

x = plt.rcParams["figure.figsize"][0]
fig, ax = plt.subplots(dpi=144, figsize=(x, x))

# Create the array to "add on" to
current_plot = np.zeros(len(bins), dtype=float)

order = np.array([1, 2, 0])

for data, name, color in zip(
    np.swapaxes(output, 0, 1)[order], np.array(["Own LR", "Other LR", "Outside LR"])[order], np.array(["C0", "C1", "C2"])[order]
):
    for x, style, fb, a in zip(data, ["", "/", "."], ["No FB", "Stellar", "AGN"], [0.6, 0.8, 1.0]):
        if style == "":
            label = name
        else:
            label=None
        
        ax.fill_between(
            bins, current_plot, current_plot+x, hatch=style,
            label=label, facecolor=color, edgecolor="white",
            alpha=a
        )
        
        current_plot += x


plt.legend(labelspacing=-2.25, loc=[0.7, 0.6])
plt.xlabel(r"$r/r_{\rm vir}$")
plt.ylabel("Fraction of mass within bin")
plt.tight_layout()
plt.ylim(0, 1)
plt.xlim(0, 1)
plt.savefig("feedback_distance_plot.pdf")
