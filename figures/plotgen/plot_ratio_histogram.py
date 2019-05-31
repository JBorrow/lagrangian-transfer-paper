import matplotlib.pyplot as plt
import numpy as np

plt.style.use("mnras_flatiron")

data = dict(gas = np.load("neighbour_analysis_gas_distance_ratio.npy"),
            star = np.load("neighbour_analysis_star_distance_ratio.npy"),
            dm = np.load("neighbour_analysis_dark_matter_distance_ratio.npy"))

bins = np.logspace(-2.5, 2, 256)
bin_centers = [0.5 * (x + y) for x, y in zip(bins[1:], bins[:-1])]

hist = {k: np.histogram(v, bins=bins)[0] for k, v in data.items()}

names = dict(gas="Gas", star="Stars", dm="Dark Matter")
colours = dict(gas=0, star=1, dm=2)

for key in data.keys():
    plt.plot(bin_centers, hist[key]/hist[key].sum(), label=names[key], color=f"C{colours[key]}")

plt.legend()

plt.semilogx()

plt.xlabel("Ratio of $D_s$ to expected $R_{\\rm vir}$")
plt.ylabel("Fraction of particles in bin")

plt.xlim(bin_centers[0], bin_centers[-1])
plt.ylim(0, None)

plt.tight_layout()

plt.savefig("ratio_distance_histogram.pdf")

