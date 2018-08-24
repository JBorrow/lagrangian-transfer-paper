"""
Runs all of the analysis scripts only loading the data once.

Invoke this in script mode as

    python3 all_analysis.py <location_of_lt_outputs.hdf5>

"""

import ltcaesar as lt

import radial_distance_analysis
import component_fraction_vs_halo_mass_analysis

import sys

filename = sys.argv[1]
print("Reading data from {}".format(filename))

sim = lt.read_data_from_file(filename)

print("Running radial_distance_analysis")
radial_distance_analysis.run_analysis(simulation=sim)

print("Running component_fraction_vs_halo_mass_analysis")
component_fraction_vs_halo_mass_analysis.run_analysis(simulation=sim)
