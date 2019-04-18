"""
Generates data files:

    + neighbour_analysis_gas_distance_ratio.npy
    + neighbour_analysis_star_distance_ratio.npy
    + neighbour_analysis_dark_matter_distance_ratio.npy

which are numpy arrays that can be analysed and properly described in the
accompanying plotting script.

Invoke this in script mode as
    
    python3 neighbour_analysis.py <location_of_lt_outputs.hdf5>
"""

import numpy as np
import ltcaesar as lt


def grab_feedback_numbers(simulation, ptype):
    """
    Returns a single array that tells you:

    + Which particles have been touched by AGN feedback (2 in array)
    + Which particle shave been touched by stellar feedback only (1)
    + Which particles have not been touched by any kind of feedback (0)
    """

    raw = simulation.snapshot_end.baryonic_matter.read_extra_array(
        "NWindLaunches", ptype
    )

    agn = raw >= 1000
    stellar = np.logical_and(raw < 1000, raw > 0)
    other = raw == 0

    output = np.empty(len(other), dtype="u1")
    output[other] = 0
    output[stellar] = 1
    output[agn] = 2

    return output


def get_rvir(simulation: lt.objects.Simulation, ptype, halo_catalogue):
    """
    Gets a numpy array of Rvir for that particle type based on their lagrangian
    regions.
    """

    if ptype == "dark_matter":
        lagrangian_regions = simulation.snapshot_end.dark_matter.halos
    elif ptype == "gas":
        lagrangian_regions = (
            simulation.snapshot_end.baryonic_matter.gas_lagrangian_regions
        )
    elif ptype == "stars":
        lagrangian_regions = (
            simulation.snapshot_end.baryonic_matter.star_lagrangian_regions
        )

    rvir = np.array([halo_catalogue[x].rvir for x in lagrangian_regions if x != -1])
    mask = lagrangian_regions != -1

    return rvir, mask


def run_analysis(simulation: lt.objects.Simulation, halo_catalogue):
    """
    Runs the analysis and saves it; this is a function in case multiple
    analysis scripts need to be ran at once (to avoid double-loading data).
    """

    rvir_star, mask_star = get_rvir(simulation, "stars", halo_catalogue)
    rvir_gas, mask_gas = get_rvir(simulation, "gas", halo_catalogue)
    rvir_dark_matter, mask_dark_matter = get_rvir(simulation, "dark_matter", halo_catalogue)

    star_data = lt.analysis.plot.find_distances_to_nearest_neighbours_data(
        simulation, "stars"
    )
    dark_matter_data = lt.analysis.plot.find_distances_to_nearest_neighbours_data(
        simulation, "dark_matter"
    )
    gas_data = lt.analysis.plot.find_distances_to_nearest_neighbours_data(
        simulation, "gas"
    )

    # Create the basic distance histograms
    for ptype in ["gas", "star", "dark_matter"]:
        data_name = f"{ptype}_data"
        this_data = locals()[data_name][1]
        this_rvir = locals()[f"rvir_{ptype}"]
        this_mask = locals()[f"mask_{ptype}"]
        distance = this_data[this_mask] / this_rvir
        np.save(f"neighbour_analysis_{ptype}_distance_ratio.npy", distance)

    return


if __name__ == "__main__":
    import sys
    import pickle

    filename = sys.argv[1]
    print("Reading data from {}".format(filename))

    with open("../halo_catalogue.pickle", "rb") as handle:
        halo_catalogue = pickle.load(handle)

    sim = lt.read_data_from_file(filename)

    run_analysis(simulation=sim, halo_catalogue=halo_catalogue.halos)
