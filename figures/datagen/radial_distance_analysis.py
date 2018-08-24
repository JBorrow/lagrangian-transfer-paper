"""
Generates data files:

    + radial_distance_analysis_gas.npy
    + radial_distance_analysis_star.npy
    + radial_distance_analysis_gas_stderr.npy
    + radial_distance_analysis_star_stderr.npy
    + radial_distance_analysis_radial_bins.npy

which are numpy arrays that can be analysed and properly described in the
accompanying plotting script.

Invoke this in script mode as
    
    python3 radial_distance_analysis.py <location_of_lt_outputs.hdf5>
"""

import numpy as np
import ltcaesar as lt


def run_analysis(simulation: lt.objects.Simulation):
    """
    Runs the analysis and saves it; this is a function in case multiple
    analysis scripts need to be ran at once (to avoid double-loading data).
    """

    bin_edges = np.linspace(0, 1, 50)
    radial_bins = [[x, y] for x, y in zip(bin_edges[:-1], bin_edges[1:])]
    output, output_std = lt.analysis.radial.run_analysis_on_mass_bin(
        simulation, 1e12, 1e13, radial_bins
    )
    output_star, output_star_std = lt.analysis.radial.run_analysis_on_mass_bin(
        simulation, 1e12, 1e13, radial_bins, ptype="star"
    )

    bin_centers = [0.5 * (x + y) for x, y in zip(bin_edges[:-1], bin_edges[1:])]

    # Make the plot pretty
    bin_centers[0] = 0
    bin_centers[1] = 1

    np.save("radial_distance_analysis_gas.npy", output)
    np.save("radial_distance_analysis_star.npy", output_star)
    np.save("radial_distance_analysis_gas_stderr.npy", output_std)
    np.save("radial_distance_analysis_star_stderr.npy", output_star_std)
    np.save("radial_distance_analysis_radial_bins.npy", bin_centers)

    return


if __name__ == "__main__":
    import sys

    filename = sys.argv[1]
    print("Reading data from {}".format(filename))

    sim = lt.read_data_from_file(filename)

    run_analysis(simulation=sim)
