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
try:
    data_inverse = np.load("inverse_component_fraction_vs_lr_mass.npy").item()
except:
    pass

def get_fancy_halo_mass(data):
    """
    Returns the halo mass 10^x'd and in physical units.
    """
    try:
        hm = data["halo_mass"]
    except KeyError:
        hm = data["lr_mass"]

    fancy_halo_mass = (1e10 / 0.7) * 10 ** (hm)

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
        plt.plot(halo_mass, this_data[name_of_item], color=color, label=label)

    plt.semilogx()

    plt.xlim(halo_mass[0], halo_mass[-1])
    current_ylim = plt.ylim()[1]
    plt.ylim(0, min([current_ylim, 1]))

    plt.xlabel("Halo mass (M$_\odot$)")
    plt.ylabel("Fraction of mass in component")

    plt.legend()
    plt.tight_layout()

    plt.savefig("component_fraction_vs_halo_mass_{}.pdf".format(data_type))

    plt.close()


# Now we make the plot of the three next to each other.

fig, ax = plt.subplots(1, 3, figsize=[6.973, 3], sharey=True, sharex=True)

for axis, data_type in zip(ax, ["both", "gas", "stellar"]):
    try:
        this_data = locals()["data_{}".format(data_type)]
    except:
        # Skip this one
        continue
    halo_mass = get_fancy_halo_mass(this_data)

    for color, label in zip(colors, switch.keys()):
        name_of_item = switch[label]
        name_of_error = "{}_stddev".format(name_of_item)

        axis.fill_between(
            halo_mass[:-1],
            (this_data[name_of_item] - this_data[name_of_error])[:-1],
            (this_data[name_of_item] + this_data[name_of_error])[:-1],
            alpha=0.2,
            color=color,
            lw=0,
        )
        axis.plot(halo_mass[:-1], this_data[name_of_item][:-1], color=color, label=label)

ax[1].semilogx()

current_ylim = max([axis.get_ylim()[1] for axis in ax])
ax[1].set_xlim(halo_mass[0], halo_mass[-2])
ax[1].set_ylim(0, min([current_ylim, 1]))

ax[1].set_xlabel("Halo mass [M$_\odot$]")
ax[0].set_ylabel("Fraction of mass in component")

ax[2].legend()

fig.tight_layout()
fig.subplots_adjust(wspace=0, hspace=0)

fig.savefig("component_fraction_mixed.pdf")

# Now for the inverse plot!

fig, ax = plt.subplots(1)

lr_mass = get_fancy_halo_mass(data_inverse)

halo_switch = {
    "Own Halo": "mass_fraction_to_halo",
    "Other Haloes": "mass_fraction_to_other_halo",
    "Outside Haloes": "mass_fraction_to_outside_halo" 
}

for color, label in zip(colors, halo_switch.keys()):
    name_of_item = halo_switch[label]
    name_of_error = "{}_stddev".format(name_of_item)

    ax.fill_between(
        lr_mass[:-1],
        (data_inverse[name_of_item] - data_inverse[name_of_error])[:-1],
        (data_inverse[name_of_item] + data_inverse[name_of_error])[:-1],
        alpha=0.2,
        color=color,
        lw=0,
    )
    ax.plot(lr_mass[:-1], data_inverse[name_of_item][:-1], color=color, label=label)

ax.semilogx()
ax.set_ylabel("Fraction of baryonic mass at $z=0$ from LR")
ax.set_xlabel("Mass of Lagrangian Region (LR) [M$_\odot$]")

current_ylim = ax.get_ylim()[1]
ax.set_xlim(lr_mass[0], lr_mass[-2])
ax.set_ylim(0, min([current_ylim, 1]))

ax.legend(loc=(0.4, 0.4))

fig.tight_layout()
fig.savefig("inverse_component_fraction.pdf")
