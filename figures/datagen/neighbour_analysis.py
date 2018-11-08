"""
Generates data files:

    + neighbour_analysis_ratio_gas.npy
    + neighbour_analysis_ratio_star.npy
    + neighbour_analysis_gas_distance.npy
    + neighbour_analysis_star_distance.npy
    + neighbour_analysis_dark_matter_distance.npy
    + neighbour_analysis_gas_fb_agn.npy
    + neighbour_analysis_gas_fb_stellar.npy
    + neighbour_analysis_gas_fb_none.npy

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

    feedback_gas = grab_feedback_numbers(simulation, "gas")

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

    fb_gas = feedback_gas

    # Perform the binning, finally...
    bins_ratio = np.logspace(-4, 5, 256)
    bins_distance = np.linspace(0, 15000, 256)

    # Create the "ratio" histogram.
    ratio_gas, _ = np.histogram(radii_gas / radii_dm_g, bins=bins_ratio)
    ratio_star, _ = np.histogram(radii_star / radii_dm_s, bins=bins_ratio)

    # Dump ratio items
    for x in ["gas", "star"]:
        full_name = f"ratio_{x}"
        np.save(f"neighbour_analysis_{full_name}.npy", locals()[full_name])

    # Create the basic distance histograms
    for ptype in ["gas", "star", "dark_matter"]:
        data_name = f"{ptype}_data"
        this_data = locals()[data_name][1]
        distance, _ = np.histogram(this_data, bins=bins_distance)
        np.save(f"neighbour_analysis_{ptype}_distance.npy", distance)

    # Create the gas distances, binned by feedback method.

    gas = {
        "AGN": gas_data[1][fb_gas == 2],
        "Stellar": gas_data[1][fb_gas == 1],
        "None": gas_data[1][fb_gas == 0],
    }

    for name, data in gas.items():
        hist, _ = np.histogram(data, bins=bins_distance)
        np.save(f"neighbour_analysis_gas_fb_{name.lower()}.npy", hist)

    return


if __name__ == "__main__":
    import sys

    filename = sys.argv[1]
    print("Reading data from {}".format(filename))

    sim = lt.read_data_from_file(filename)

    run_analysis(simulation=sim)
