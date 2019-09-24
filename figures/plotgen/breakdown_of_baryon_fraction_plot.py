"""
Uses data files:
    
    + baryon_fraction_vs_halo_mass_own.npy
    + baryon_fraction_vs_halo_mass_other.npy
    + baryon_fraction_vs_halo_mass_outside.npy
    + baryon_fraction_halo_masses.npy


This script makes figures that show the trends of lagrangian transfer as a
function of halo mass.

Invoke this in script mode as
    
    python3 breakdown_of_baryon_fraction_plot.py
"""

import matplotlib.pyplot as plt
import numpy as np

plt.style.use("mnras_flatiron")


def get_fancy_halo_mass(data):
    """
    Returns the halo mass 10^x'd and in physical units.
    """
    fancy_halo_mass = (1e10 / 0.7) * 10 ** (data)

    return fancy_halo_mass


own = np.load("baryon_fraction_vs_halo_mass_own.npy")
other = np.load("baryon_fraction_vs_halo_mass_other.npy")
outside = np.load("baryon_fraction_vs_halo_mass_outside.npy")
halo_mass = get_fancy_halo_mass(np.load("baryon_fraction_halo_masses.npy"))

switch = {"Outside LR": "outside", "Other LR": "other", "Own LR": "own"}

colors = ["C{}".format(x) for x in reversed(range(3))]

previous = [np.zeros_like(own)]

fig, ax = plt.subplots()

for (label, name), c in zip(switch.items(), colors):
    current_data = locals()["{}".format(name)]

    ax.fill_between(
        halo_mass,
        sum(previous),
        sum(previous) + current_data,
        color=c,
        label=label,
        alpha=0.6,
        linewidth=0,
    )

    previous.append(current_data)

# Now to draw the overall trend.

ax.plot(halo_mass, sum(previous), color="black", label="Overall")
ax.scatter(halo_mass, sum(previous), color="black", s=5)

ax.semilogx()
ax.set_ylabel("Baryon fraction $f_b / f_{b, c}$")
ax.set_xlabel("Halo mass [M$_\odot$]")

ax.set_xlim(halo_mass[0], halo_mass[-2])
ax.set_ylim(0, 1.0)

ax.legend()

fig.tight_layout()
fig.savefig("baryon_fraction_breakdown.pdf")
