"""
Uses data files:
    
    + component_fraction_vs_halo_mass_gas.npy
    + component_fraction_vs_halo_mass_stellar.npy
    + component_fraction_vs_halo_mass_both.npy

This script makes figures that show the trends of lagrangian transfer as a
function of halo mass.

Invoke this in script mode as
    
    python3 component_fraction_vs_halo_mass_plot.py
"""

import matplotlib.pyplot as plt
import numpy as np

plt.style.use("mnras_flatiron")

data_gas = np.load("component_fraction_vs_halo_mass_gas.npy").item()
data_stellar = np.load("component_fraction_vs_halo_mass_stellar.npy").item()
data_both = np.load("component_fraction_vs_halo_mass_both.npy").item()
try:
    data_dm = np.load("component_fraction_vs_halo_mass_dm.npy").item()
except:
    # No stuff :(
    pass

def get_fancy_halo_mass(data):
    """
    Returns the halo mass 10^x'd and in physical units.
    """
    fancy_halo_mass = (1e10 / 0.7) * 10 ** (data["halo_mass"])

    return fancy_halo_mass


switch = {
    "Own LR": "mass_fraction_from_lr",
    "Other LR": "mass_fraction_from_other_lr",
    "Outside LR": "mass_fraction_from_outside_lr",
}

colors = ["C{}".format(x) for x in range(3)]


# Make the three plots

for data_type in ["gas", "stellar", "both", "dm"]:
    try:
        this_data = locals()["data_{}".format(data_type)]
    except:
        # Skip this one
        continue
    halo_mass = get_fancy_halo_mass(this_data)

    for color, label in zip(colors, switch.keys()):
        name_of_item = switch[label]
        name_of_error = "{}_stddev".format(name_of_item)

        plt.fill_between(
            halo_mass,
            this_data[name_of_item] - this_data[name_of_error],
            this_data[name_of_item] + this_data[name_of_error],
            alpha=0.2,
            color=color,
            lw=0,
        )
        plt.plot(
            halo_mass, this_data[name_of_item], color=color, label=label
        )

    plt.semilogx()

    plt.xlim(halo_mass[0], halo_mass[-2])
    current_ylim = plt.ylim()[1]
    plt.ylim(0, min([current_ylim, 1]))

    plt.xlabel("Halo mass (M$_\odot$)")
    plt.ylabel("Fraction of mass in component")

    plt.legend()
    plt.tight_layout()

    plt.savefig("component_fraction_vs_halo_mass_{}.pdf".format(data_type))

    plt.close()
