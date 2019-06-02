"""
Uses data files:
    
    + component_fraction_vs_halo_mass_gas.npy
    + component_fraction_vs_halo_mass_stellar.npy
    + component_fraction_vs_halo_mass_both.npy

from all of the various smoothings.

This script makes figures that show the trends of lagrangian transfer as a
function of halo mass.

Invoke this in script mode as
    
    python3 component_fraction_vs_halo_mass_plot.py
"""
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

plt.style.use("mnras_flatiron")

smoothings = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
alphas = np.linspace(1.0, 0.2, len(smoothings))

# To make the stuff at the top of the plot, we need a custom 'colourbar'
# in matplotlib parlance. We want to create a colour bar that shows
# the colours 0, 1.

from matplotlib.colors import LinearSegmentedColormap, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable

colors = [(83/255, 126/255, 186/255, a) for a in alphas]
cm = LinearSegmentedColormap.from_list("FlatironBlues", colors, N=len(smoothings))

def colorbar(mappable): 
    ax = mappable.axes 
    fig = ax.figure 
    divider = make_axes_locatable(ax) 
    cax = divider.append_axes("top", size="5%", pad=0) 
    return fig.colorbar(mappable, cax=cax, orientation="horizontal")

def update_colorbar(cbar):
    cbar.ax.set_xlabel(r"$R_{\rm LR} / R_{\rm vir}$")
    cbar.ax.xaxis.set_label_position("top") 
    cbar.ax.get_xaxis().set_ticks([])
 
    # Now we put the data ontop of the colorbar
    d_alpha = 1.0 / len(smoothings)
    alpha_0 = 0.0
    for n, smoothing in enumerate(smoothings):
        alpha = alpha_0 + (n + 0.5) * d_alpha
        cbar.ax.text(alpha, 0.4, smoothing, ha="center", va="center", alpha=0.8)
    
    return
        

data_gas = []
data_stellar = []
data_both = []
data_dm = []
data_inverse = []

for smoothing in smoothings:
    data_gas.append(np.load(f"{smoothing}/plots/component_fraction_vs_halo_mass_gas.npy").item())
    data_stellar.append(np.load(f"{smoothing}/plots/component_fraction_vs_halo_mass_stellar.npy").item())
    data_both.append(np.load(f"{smoothing}/plots/component_fraction_vs_halo_mass_both.npy").item())
    try:
        data_dm.append(np.load(f"{smoothing}/plots/component_fraction_vs_halo_mass_dm.npy").item())
    except:
        # No stuff :(
        pass
    try:
        data_inverse.append(np.load(f"{smoothing}/plots/inverse_component_fraction_vs_lr_mass.npy").item())
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

for data_type in ["gas", "stellar", "both"]:
    try:
        this_data = locals()["data_{}".format(data_type)]
    except:
        # Skip this one
        continue
    halo_mass = get_fancy_halo_mass(this_data[0])


    for color, label in zip(colors, switch.keys()):

        name_of_item = switch[label]

        for alpha, smooth_data in zip(alphas, this_data):
            if alpha == alphas[0]:
                this_label = label
            else:
                this_label = ""
            plt.plot(halo_mass, smooth_data[name_of_item], color=color, label=this_label, alpha=alpha)

    plt.semilogx()

    plt.xlim(halo_mass[0], halo_mass[-1])
    current_ylim = plt.ylim()[1]
    plt.ylim(0, min([current_ylim, 1]))


    plt.xlabel("Halo mass (M$_\odot$)")
    plt.ylabel("Fraction of mass in component")

    plt.legend()

    # Generate our colour map
    mappable = plt.scatter([1.0]*len(smoothings), [1.0]*len(smoothings), c=Normalize()(smoothings), cmap=cm)
    update_colorbar(colorbar(mappable))

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
    halo_mass = get_fancy_halo_mass(this_data[0])

    for color, label in zip(colors, switch.keys()):
        name_of_item = switch[label]

        for (i, alpha), smooth_data in zip(enumerate(alphas), this_data):
            if axis != ax[0]:
                if alpha == alphas[-1]:
                    this_label = label
                else:
                    this_label = ""
            else:
                if label == "Own LR":
                    this_label = smoothings[i]
                else:
                    this_label = ""


            axis.plot(halo_mass[:-1], smooth_data[name_of_item][:-1], color=color, label=this_label, alpha=alpha)

ax[1].semilogx()

current_ylim = max([axis.get_ylim()[1] for axis in ax])
ax[1].set_xlim(halo_mass[0], halo_mass[-2])
ax[1].set_ylim(0, min([current_ylim, 1]))

ax[1].set_xlabel("Halo mass [M$_\odot$]")
ax[0].set_ylabel("Fraction of mass in component")

ax[0].legend()
ax[2].legend()

fig.tight_layout()
fig.subplots_adjust(wspace=0, hspace=0)

fig.savefig("component_fraction_mixed.pdf")

# Now for the inverse plot!

fig, ax = plt.subplots(1)

lr_mass = get_fancy_halo_mass(data_inverse[0])

halo_switch = {
    "Own Halo": "mass_fraction_to_halo",
    "Other Halos": "mass_fraction_to_other_halo",
    "Outside Halos": "mass_fraction_to_outside_halo" 
}

for color, label in zip(colors, halo_switch.keys()):
    for color, label in zip(colors, switch.keys()):
        name_of_item = switch[label]

        for alpha, smooth_data in zip(alphas, this_data):
            if alpha == alphas[-1]:
                this_label = label
            else:
                this_label = ""
            ax.plot(halo_mass[:-1], smooth_data[name_of_item][:-1], color=color, label=this_label, alpha=alpha)

ax.semilogx()
ax.set_ylabel("Fraction of baryonic mass at $z=0$ from LR")
ax.set_xlabel("Mass of Lagrangian Region (LR) [M$_\odot$]")

current_ylim = ax.get_ylim()[1]
ax.set_xlim(lr_mass[0], lr_mass[-2])
ax.set_ylim(0, min([current_ylim, 1]))

ax.legend(loc=(0.4, 0.4))

fig.tight_layout()
fig.savefig("inverse_component_fraction.pdf")


