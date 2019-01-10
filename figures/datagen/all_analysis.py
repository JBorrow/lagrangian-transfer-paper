"""
Runs all of the analysis scripts only loading the data once.

Invoke this in script mode as

    python3 all_analysis.py <location_of_lt_outputs.hdf5>

"""

import ltcaesar as lt

import radial_distance_analysis
import component_fraction_vs_halo_mass_analysis
import neighbour_analysis
import split_by_feedback
import breakdown_of_baryon_fraction

import sys

filename = sys.argv[1]
print("Reading data from {}".format(filename))

sim = lt.read_data_from_file(filename)

print("Running baryon fraction analysis")
breakdown_of_baryon_fraction.run_analysis(simulation=sim)

print("Running component_fraction_vs_halo_mass_analysis")
component_fraction_vs_halo_mass_analysis.run_analysis(simulation=sim)

print("Running radial_distance_analysis")
radial_distance_analysis.run_analysis(simulation=sim)

print("Running neighbour_anaylsis")
neighbour_analysis.run_analysis(simulation=sim)

print("Running feedback analysis")
split_by_feedback.run_analysis(simulation=sim)
