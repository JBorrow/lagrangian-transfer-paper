"""
Generates data files:

    + baryon_fraction_vs_halo_mass_own.npy
    + baryon_fraction_vs_halo_mass_other.npy
    + baryon_fraction_vs_halo_mass_outside.npy
    + baryon_fraction_vs_halo_mass_own_stddev.npy
    + baryon_fraction_vs_halo_mass_other_stddev.npy
    + baryon_fraction_vs_halo_mass_outside_stddev.npy
    + baryon_fraction_halo_masses.npy

which are numpy arrays that can be analysed and properly described in the
accompanying plotting script.

Invoke this in script mode as
    
    python3 breakdown_of_baryon_fraction.py <location_of_lt_outputs.hdf5>
"""

import numpy as np
import ltcaesar as lt
from ltcaesar.analysis.plot import bin_x_by_y


def run_analysis(simulation: lt.objects.Simulation):
    """
    Runs the analysis and saves it; this is a function in case multiple
    analysis scripts need to be ran at once (to avoid double-loading data).
    """

    # Cosmic baryon fraction (This is a constant from Planck.)
    f_b = 0.15804394

    bins = np.linspace(1, 5, 30)

    # First, we need to mask out the halos with no gas or stellar content.

    use = ["gas", "stellar"]
    average_func = np.median

    bins = bins = np.linspace(1, 5, 30)

    masks = [getattr(simulation, "{}_mass_in_halo".format(x)) != 0 for x in use]
    mask = np.logical_and.reduce(masks)

    # Grab a bunch of masked arrays that are going to be used in the analysis
    masked_mass = (simulation.dark_matter_mass_in_halo)[mask]
    masked_lagrangian_mass = sum(
        [getattr(simulation, "{}_mass_in_halo_from_lagrangian".format(x))[mask] for x in use]
    )
    masked_other_lagrangian_mass = sum(
        [
            getattr(simulation, "{}_mass_in_halo_from_other_lagrangian".format(x))[mask]
            for x in use
        ]
    )
    masked_outside_lagrangian_mass = sum(
        [
            getattr(simulation, "{}_mass_in_halo_from_outside_lagrangian".format(x))[mask]
            for x in use
        ]
    )

    # Now reduce the data into fractions

    fraction_of_mass_from_lr = masked_lagrangian_mass / (masked_mass * f_b)
    fraction_of_mass_from_other_lr = masked_other_lagrangian_mass / (masked_mass * f_b)
    fraction_of_mass_from_outside_lr = masked_outside_lagrangian_mass / (masked_mass * f_b)

    masked_halo_masses = np.log10(simulation.dark_matter_mass_in_halo[mask])

    # Use local routine to fully reduce the data into a single line

    halo_mass, mass_fraction_from_own_lr, mass_fraction_from_own_lr_stddev = bin_x_by_y(
        masked_halo_masses, fraction_of_mass_from_lr, bins, average_func
    )

    _, mass_fraction_from_other_lr, mass_fraction_from_other_lr_stddev = bin_x_by_y(
        masked_halo_masses, fraction_of_mass_from_other_lr, bins, average_func
    )

    _, mass_fraction_from_outside_lr, mass_fraction_from_outside_lr_stddev = bin_x_by_y(
        masked_halo_masses, fraction_of_mass_from_outside_lr, bins, average_func
    )

    np.save("baryon_fraction_halo_masses.npy", halo_mass)

    for name in ["own", "outside", "other"]:
        filename = "baryon_fraction_vs_halo_mass_{}.npy".format(name)
        data = locals()["mass_fraction_from_{}_lr".format(name)]
        np.save(filename, data)

        filename = "baryon_fraction_vs_halo_mass_{}_stddev.npy".format(name)
        data = locals()["mass_fraction_from_{}_lr_stddev".format(name)]
        np.save(filename, data)

    return


if __name__ == "__main__":
    import sys

    filename = sys.argv[1]
    print("Reading data from {}".format(filename))

    sim = lt.read_data_from_file(filename)

    run_analysis(simulation=sim)
