"""
Generates data files:

    + neighbour_analysis_data_gas.npy
    + neighbour_analysis_data_dark_matter.npy
    + neighbour_analysis_data_star.npy
    + neighbour_analysis_radii_dm_g.npy
    + neighbour_analysis_radii_dm_s.npy
    + neighbour_analysis_radii_gas.npy
    + neighbour_analysis_radii_star.npy

which are numpy arrays that can be analysed and properly described in the
accompanying plotting script.

Note that these arrays end up being quite large.

Invoke this in script mode as
    
    python3 neighbour_analysis.py <location_of_lt_outputs.hdf5>
"""

import numpy as np
import ltcaesar as lt


def run_analysis(simulation: lt.objects.Simulation):
    """
    Runs the analysis and saves it; this is a function in case multiple
    analysis scripts need to be ran at once (to avoid double-loading data).
    """

    star_data = lt.analysis.plot.find_distances_to_nearest_neighbours_data(
        simulation, "stars"
    )
    dark_matter_data = lt.analysis.plot.find_distances_to_nearest_neighbours_data(
        simulation, "dark_matter"
    )
    gas_data = lt.analysis.plot.find_distances_to_nearest_neighbours_data(
        simulation, "gas"
    )

    # Data for the histogram

    # First we need to figure out particles which have gas neighbours that
    # are also still alive in the final conditions
    unique_indicies_dm_gas = np.unique(dark_matter_data[3])[
        np.in1d(np.unique(dark_matter_data[3]), np.unique(gas_data[3]))
    ]

    unique_indicies_dm_star = np.unique(dark_matter_data[3])[
        np.in1d(np.unique(dark_matter_data[3]), np.unique(star_data[3]))
    ]

    ids_dm_g = np.take(sim.snapshot_end.dark_matter.ids, unique_indicies_dm_gas)
    ids_gas = np.take(sim.snapshot_end.dark_matter.ids, unique_indicies_dm_gas)
    ids_dm_s = np.take(sim.snapshot_end.dark_matter.ids, unique_indicies_dm_star)
    ids_star = np.take(sim.snapshot_end.dark_matter.ids, unique_indicies_dm_star)

    assert len(ids_gas) == len(ids_dm_g)
    assert len(ids_star) == len(ids_dm_s)

    indicies_dm_g = np.argsort(ids_dm_g)
    indicies_gas = np.argsort(ids_gas)
    indicies_dm_s = np.argsort(ids_dm_s)
    indicies_star = np.argsort(ids_star)

    radii_dm_g = dark_matter_data[1][indicies_dm_g]
    radii_gas = gas_data[1][indicies_gas]
    radii_dm_s = dark_matter_data[1][indicies_dm_s]
    radii_star = star_data[1][indicies_star]

    # Now we can dump these to arrays
    for x in ["dm_g", "gas", "dm_s", "star"]:
        full_name = "radii_{}".format(x)
        np.save("neighbour_analysis_{}.npy".format(full_name), locals()[full_name])

    # Move on to the more basic analysis
    for x in ["gas", "star", "dark_matter"]:
        full_name = "data_{}".format(x)
        np.save("neighbour_analysis_{}".format(full_name), locals()[full_name])

    return


if __name__ == "__main__":
    import sys

    filename = sys.argv[1]
    print("Reading data from {}".format(filename))

    sim = lt.read_data_from_file(filename)

    run_analysis(simulation=sim)
