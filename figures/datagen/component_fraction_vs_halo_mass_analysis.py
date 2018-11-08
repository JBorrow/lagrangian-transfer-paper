"""
Generates data files:

    + component_fraction_vs_halo_mass_stellar.npy
    + component_fraction_vs_halo_mass_gas.npy
    + component_fraction_vs_halo_mass_both.npy

which are numpy arrays that can be analysed and properly described in the
accompanying plotting script.

Invoke this in script mode as
    
    python3 component_fraction_vs_halo_mass_analysis.py <location_of_lt_outputs.hdf5>
"""

import numpy as np
import ltcaesar as lt


def run_analysis(simulation: lt.objects.Simulation):
    """
    Runs the analysis and saves it; this is a function in case multiple
    analysis scripts need to be ran at once (to avoid double-loading data).
    """

    data_stellar = lt.analysis.plot.mass_fraction_transfer_from_lr_data(
        simulation, bins=np.linspace(1, 5, 30), average_func=np.mean, use=["stellar"]
    )

    data_gas = lt.analysis.plot.mass_fraction_transfer_from_lr_data(
        simulation, bins=np.linspace(1, 5, 30), average_func=np.mean, use=["gas"]
    )

    data_both = lt.analysis.plot.mass_fraction_transfer_from_lr_data(
        simulation, bins=np.linspace(1, 5, 30), average_func=np.mean, use=["gas", "stellar"]
    )

    for name in ["stellar", "gas", "both"]:
        filename = "component_fraction_vs_halo_mass_{}.npy".format(name)
        data = locals()["data_{}".format(name)]

        np.save(filename, data)

    return


if __name__ == "__main__":
    import sys

    filename = sys.argv[1]
    print("Reading data from {}".format(filename))

    sim = lt.read_data_from_file(filename)

    run_analysis(simulation=sim)
