"""
Generates data files:

+ distance_gas.npy
+ distance_dm.npy
+ distance_star.npy

which are numpy arrays that can be analysed and properly described in the
accompanying plotting script.

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

    for name in ["star", "dark_matter", "gas"]:
        data = locals()[f"{name}_data"]
        np.save(f"{name}_distance.npy", data)

    return


if __name__ == "__main__":
    import sys

    filename = sys.argv[1]
    print("Reading data from {}".format(filename))

    sim = lt.read_data_from_file(filename)

    run_analysis(simulation=sim)
