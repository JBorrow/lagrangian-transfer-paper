import numpy as np
import ltcaesar as lt
from functools import lru_cache

@lru_cache(128)
def grab_feedback_numbers(halo, simulation, ptype):
    """
    Returns three boolean arrays that tell you:

    + Which particles have been touched by AGN feedback
    + Which particle shave been touched by stellar feedback only
    + Which particles have not been touched by any kind of feedback
    """
    
    raw = lt.analysis.radial.get_extra_baryonic_quantity(
        halo,
        simulation,
        "NWindLaunches",
        ptype
    )

    agn = raw >= 1000
    stellar = np.logical_and(raw < 1000, raw > 0)
    other = raw == 0

    return agn, stellar, other


def mass_fraction_binned_by_feedback(mask, halo, simulation, ptype):
    """
    Finds the fraction of mass from inside, outside, and from other
    lagrangian regions, in a given mass bin.
    """

    radii, masses, lagrangian_regions, halos = lt.analaysis.radial.get_relevant_baryonic_quantities(
        halo, simulation, ptype
    )

    agn, stellar, other = grab_feedback_numbers(halo, simulation, ptype)

    relevant_lr = lagrangian_regions[mask]
    relevant_halos = halos[mask]
    # Use fractional masses
    relevant_masses = masses[mask]
    # These need to be normalized on a bin-by-bin basis
    relevant_masses /= np.sum(relevant_masses)

    relevant_agn = agn[mask]
    relevant_stellar = stellar[mask]
    relevant_other = other[mask]

    # Create an indexing array based on those masks
    feedback_indicies = np.zeros_like(relevant_masses)
    feedback_indicies[relevant_agn] = 2
    feedback_indicies[relevant_stellar] = 1

    assert np.isclose(
        np.sum(relevant_masses), 1
    ), f"Sum: {np.sum(relevant_masses)}, length: {len(relevant_masses)}"

    from_own_lr = np.zeros(3, dtype=float)
    from_other_lr = np.zeros(3, dtype=float)
    from_outside_lr = np.zeros(3, dtype=float)

    for lr, halo, mass, fb in zip(relevant_lr, relevant_halos, relevant_masses, feedback_indicies):
        if lr == halo:
            from_own_lr[fb] += mass
        elif lr == -1:
            from_outside_lr[fb] += mass
        elif lr != -1:
            from_other_lr[fb] += mass
        else:
            raise Exception("Particle unassigned")

    sum_of_all = np.sum(from_own_lr + from_other_lr + from_outside_lr)

    assert np.isclose(sum_of_all, 1), (
        f"From Own: {from_own_lr}, From Other: {from_other_lr} ",
        f"From Outside: {from_outside_lr}, sum: {sum_of_all}.",
    )

    return from_own_lr, from_other_lr, from_outside_lr


def run_analysis(simulation: lt.objects.Simulation):
    """
    Runs the analysis and saves it; this is a function in case multiple
    analysis scripts need to be ran at once (to avoid double-loading data).
    """

    bin_edges = np.linspace(0, 1, 50)
    radial_bins = [[x, y] for x, y in zip(bin_edges[:-1], bin_edges[1:])]
    output, output_std = lt.analysis.radial.run_analysis_on_mass_bin(
        simulation, 1e12, 1e13, radial_bins, bin_func=mass_fraction_binned_by_feedback
    )

    bin_centers = [0.5 * (x + y) for x, y in zip(bin_edges[:-1], bin_edges[1:])]

    # Make the plot pretty
    bin_centers[0] = 0
    bin_centers[-1] = 1

    np.save("radial_distance_analysis_gas_fb.npy", output)
    np.save("radial_distance_analysis_gas_stderr_fb.npy", output_std)
    np.save("radial_distance_analysis_radial_bins_fb.npy", bin_centers)

    return


if __name__ == "__main__":
    import sys

    filename = sys.argv[1]
    print("Reading data from {}".format(filename))

    sim = lt.read_data_from_file(filename)

    run_analysis(simulation=sim)
