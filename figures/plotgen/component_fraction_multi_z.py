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

# redshift: model: file path
sim_names = {
    0: {"NoJet": "s50nojetAHF", "Full": "s50j7kAHF"},
    2: {"NoJet": "s50nojetAHF_z2", "Full": "s50j7kAHF_z2"},
}
sims = {0: {"NoJet": {}, "Full": {}}, 2: {"NoJet": {}, "Full": {}}}

colours = {"Gas": "C0", "Stars": "C1", "Dark Matter": "C2"}
linestyles = {0: "solid", 2: "dashed"}

for redshift, simulations in sim_names.items():
    for name, path in simulations.items():
        sims[redshift][name]["Gas"] = np.load(
            f"{path}/component_fraction_vs_halo_mass_gas.npy"
        ).item()


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

fig, axes = plt.subplots(
    2, 1, figsize=(3.321, 3.321 * 1.62), sharex=True, sharey=True, squeeze=True
)
axes = axes.flatten()

legends = [{}, {}]
label_done = {"NoJet": False, "Full": False}

for redshift, simulation in sims.items():
    for (index, axis), legend, (model, data) in zip(
        enumerate(axes), legends, simulation.items()
    ):
        for ptype, fraction_data in data.items():
            linestyle = linestyles[redshift]
            halo_mass = get_fancy_halo_mass(fraction_data)

            for color, label in zip(colors, switch.keys()):
                name_of_item = switch[label]
                name_of_error = "{}_stddev".format(name_of_item)
                
                if redshift == 2:
                    # Only create the error bars for z=2 results.
                    axis.fill_between(
                        halo_mass,
                        fraction_data[name_of_item] - fraction_data[name_of_error],
                        fraction_data[name_of_item] + fraction_data[name_of_error],
                        alpha=0.2,
                        color=color,
                        lw=0,
                    )


                line = axis.plot(
                    halo_mass,
                    fraction_data[name_of_item],
                    color=color,
                    linestyle=linestyle,
                )[0]

                # Deal with legends... ==0 for explicity
                if index == 0:
                    legend[f"$z={redshift}$"] = line
                elif index == 1:
                    if not label in legend.keys():
                        legend[label] = line

            # Add on simulation caption
            if not label_done[model]:
                axis.text(
                    0.025, 0.96, model, transform=axis.transAxes, ha="left", va="top"
                )
                label_done[model] = True


axes[0].semilogx()

axes[0].set_xlim(halo_mass[0], halo_mass[-1])
axes[0].set_ylim(0, 0.95)

axes[0].set_ylabel("Fraction of gas mass in halo", y=-0.0)
axes[1].set_xlabel("Halo mass [M$_\odot$]")

# Create the two legends.


axes[1].legend(
    legends[0].values(), legends[0].keys(), markerfirst=False, loc="center right"
)
axes[0].legend(
    legends[1].values(), legends[1].keys(), markerfirst=False, loc="center right"
)

fig.tight_layout()
fig.subplots_adjust(hspace=0)
fig.savefig("component_fraction_multi_z.pdf")
