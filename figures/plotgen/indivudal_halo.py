#%% [markdown]
# # Large Images
# We're going to test this notebook thing out over the internet.
# This script makes an image of the largest halo in the box and
# projects the fraction of gas from the lagrangian region


#%%
# Our imports
from ltcaesar import read_data_from_file
from sphviewer.tools import QuickView
import numpy as np
import matplotlib.pyplot as plt
from pickle import load
from matplotlib.colors import LogNorm, Normalize

#%%
# Setup our favourite stylesheet
plt.style.use("mnras_flatiron")

#%%
directory = "s50j7kAHF"
data = read_data_from_file(f"{directory}/lt/lt_outputs.hdf5")

#%% [markdown]
# We want to find the largest halo in the box and make a picture of it.
# Thankfully that should be reasonably easy to find because we have
# stored everything in our 'FakeCaesar' object.

#%%
with open(f"{directory}/lt/halo_catalogue.pickle", "rb") as handle:
    halo_catalogue = load(handle)


#%%
centers = [x.center for x in halo_catalogue.halos]
r_virs = [x.rvir for x in halo_catalogue.halos]

#%%
largest_rvir = r_virs[0]
largest_center = centers[0]

#%%
distance_to_plot = 2.5 * 1000  # kpc / h

#%%
gas_coordinates = data.snapshot_end.baryonic_matter.gas_coordinates

#%%
mask = np.logical_and.reduce(
    [
        gas_coordinates[:, 0] > largest_center[0] - distance_to_plot,
        gas_coordinates[:, 0] < largest_center[0] + distance_to_plot,
        gas_coordinates[:, 1] > largest_center[1] - distance_to_plot,
        gas_coordinates[:, 1] < largest_center[1] + distance_to_plot,
        gas_coordinates[:, 2] > largest_center[2] - distance_to_plot,
        gas_coordinates[:, 2] < largest_center[2] + distance_to_plot,
    ]
)

#%%
gas_coordinates_around_center = gas_coordinates[mask]

#%%
# First we'll set up some plotting commands
kwargs = dict(xsize=2048, ysize=2048, r="infinity", logscale=False, plot=False)

#%%
qv_nomass = QuickView(gas_coordinates_around_center, **kwargs)

#%%
fig, ax = plt.subplots()
ax.axis("off")
ax.imshow(
    qv_nomass.get_image(),
    norm=LogNorm(),
    extent=[-distance_to_plot, distance_to_plot, -distance_to_plot, distance_to_plot],
)

# Plot the virial radius of the cluster on top as a check
circle = plt.Circle((0, 0), largest_rvir, linestyle="dashed", color="white", fill=False)
ax.add_artist(circle)

#%% [markdown]
# Now that we've been able to get an image of our simulated cluster out, we
# should be able to plot the fraction of stuff from our LR and from outside.

#%%
baryons = data.snapshot_end.baryonic_matter
lagrangian_regions = baryons.gas_lagrangian_regions[mask]
# Will need to change this in the future if you want to plot a halo that's not 0
mask_in_lr = lagrangian_regions == 0

#%%
# Must use the hsml from previous image for consistency
qv_in_lr = QuickView(
    gas_coordinates_around_center[mask_in_lr], hsml=qv_nomass.get_hsml(), **kwargs
)

#%%
# Taking the ratio with the above will give the fraction of particles from own LR
# in each projected bin
fraction_from_in_lr = qv_in_lr.get_image() / qv_nomass.get_image()

#%%
fig, ax = plt.subplots()

pixels = np.linspace(-distance_to_plot, distance_to_plot, kwargs["xsize"])
x, y = np.meshgrid(pixels, pixels)

mappable = ax.pcolormesh(
    x,
    y,
    fraction_from_in_lr,
    norm=Normalize(vmin=0, vmax=1),
    cmap="BuPu",
    rasterized=True,
)

fig.colorbar(mappable, pad=0, label="Fraction of particles from own LR")

# Add a contour map of the actual density projection on top.
ax.contour(
    x,
    y,
    LogNorm()(qv_nomass.get_image()),
    colors="#333333",
    linewidths=0.35,
    antialiased=True,
)

# Plot the virial radius of the cluster on top as a check
circle = plt.Circle((0, 0), largest_rvir, linestyle="dashed", color="white", fill=False)
ax.add_artist(circle)

ax.tick_params(
    axis="both", bottom=False, left=False, labelbottom=False, labelleft=False
)

# Do not tight layout this, it will stretch the central circle to an oval.
# fig.tight_layout()

fig.savefig("fraction_from_own_lr_largest_halo.pdf")
