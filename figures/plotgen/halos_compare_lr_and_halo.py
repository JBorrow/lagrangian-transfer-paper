#%% [markdown]
# # Halo v.s. Lagrangian Regions
# Let's look at how large the difference is between where the gas comes
# from in the initial conditions to where the dark matter comes from.

#%%
# Our imports
from ltcaesar import read_data_from_file
from sphviewer.tools import QuickView
import h5py
import numpy as np
import matplotlib.pyplot as plt
from pickle import load
from matplotlib.colors import LogNorm, Normalize
from matplotlib.cm import Blues, Reds

#%%
# Setup our favourite stylesheet
plt.style.use("mnras_flatiron")

#%%
directory = "s50j7kAHF"
data = read_data_from_file(f"{directory}/lt/lt_outputs.hdf5")

#%%
baryons_ini = data.snapshot_ini.baryonic_matter
baryons_end = data.snapshot_end.baryonic_matter
dm_ini = data.snapshot_ini.dark_matter
kwargs = dict(xsize=2048, ysize=2048, r="infinity", logscale=False, plot=False)

#%% [markdown]
# We can grab the smoothing lengths out of the data.

#%%
with h5py.File(data.snapshot_ini.snapshot_filename, "r") as handle:
    mask = np.argsort(handle["PartType0/ParticleIDs"][...])

    hsml_ini_gas = handle["PartType0/SmoothingLength"][...][mask]


#%%
def lagrangian_region(lr):
    """
    Generate an image (returns a pixel grid) for the gas in the lagrangian region
    specified.
    """

    mask = data.lagrangian_regions == lr

    c = np.concatenate(
        [
            baryons_ini.gas_coordinates[mask],
            np.array([[0.0, 0.0, 0.0], [50000.0, 50000.0, 50000.0]]),
        ]
    )
    hsml = np.concatenate([hsml_ini_gas[mask], [0.0, 0.0]])

    qv = QuickView(c, hsml=hsml, **kwargs).get_image()

    return qv


#%%
def halo_gas(halo):
    """
    Generate an image (returns a pixel grid) for the gas in the halo
    specified.
    """

    baryonic_ids_in_halo_at_z0 = np.unique(
        np.concatenate(
            [
                baryons_end.gas_ids[baryons_end.gas_halos == halo],
                baryons_end.star_ids[baryons_end.star_halos == halo],
            ]
        )
    )

    _, mask, _ = np.intersect1d(
        baryons_ini.gas_ids,
        baryonic_ids_in_halo_at_z0,
        assume_unique=True,
        return_indices=True,
    )

    c = np.concatenate(
        [
            baryons_ini.gas_coordinates[mask],
            np.array([[0.0, 0.0, 0.0], [50000.0, 50000.0, 50000.0]]),
        ]
    )
    hsml = np.concatenate([hsml_ini_gas[mask], [0.0, 0.0]])

    qv = QuickView(c, hsml=hsml, **kwargs).get_image()

    return qv


#%%
def cmap_image(image, color=0):
    """
    Colormaps an image in the following way:
    
    + First, we LogNorm() the image
    + Pixels are all given 255 in the color axis if nonzero
    + Alpha is set as the LogNorm of hte image
    """

    normed = np.ma.filled(Normalize()(image), fill_value=0.0)
    colors = np.zeros_like(normed)

    colors[normed != 0.0] = 1.0

    out = np.zeros((kwargs["xsize"], kwargs["ysize"], 4), dtype=float)
    out[:, :, color] = normed
    out[:, :, 3] = normed

    print(out)
    print(out.shape)

    return out


#%%
def combine_images(top, bottom):
    """
    Combines two images by summing their colours and meaning their alphas.
    """

    return (top + bottom) * 0.5


#%%
def create_combined_image(halo_id):
    return combine_images(
        Blues(Normalize()(lagrangian_region(halo_id))),
        Reds(Normalize()(halo_gas(halo_id))),
    )


#%% [markdown]
# Fantastic, we can create the blended image! Now we just need to implement darker
# colour so that we can have multiple halos on at once. This just looks at the
# brightness of each pixel, and then chooses the 'brightest' one of the three.

#%%
def darker_colour(top, bottom):
    """
    Chooses the darker colour of the two for each pixel.
    """

    def lum(x):
        """ Brightness of the image """
        return x[:, :, :3].sum(axis=2) / 3.0

    # They're gonna have the same colour in a lot of places, let's make
    # sure that's included
    which_brightest = lum(top) <= lum(bottom)

    print(which_brightest.shape)

    out = bottom.copy()

    out[which_brightest] = top[which_brightest]

    return out


#%%
# Quickly import reduce here, idk why it's not part of the standard library
from functools import reduce
from tqdm import tqdm

#%%
group_ids_to_plot = [200, 212, 215, 219, 323, 324]  # [0, 431, 88, 299, 9]

#%%
images = map(create_combined_image, tqdm(group_ids_to_plot))

#%%
one_image = reduce(darker_colour, images)

#%%
fig, ax = plt.subplots(figsize=(8, 8))
ax.axis("off")
fig.subplots_adjust(0, 0, 1, 1)

ax.imshow(one_image, origin="lower")
fig.savefig("lagrangian_regions_combined_secondset.png", dpi=300)

#%% [markdown]
# We need to make the dark matter for the background!

#%%
qv = QuickView(data.snapshot_end.dark_matter.coordinates, **kwargs)

#%%
plt.imshow(qv.get_image, cmap="bone", norm=LogNorm())
plt.gca().axis("off")

#%% [markdown]
# On top of everyone's favourite dark matter image, we also want to show the
# halo virial radii.

#%%
with open(f"{directory}/lt/halo_catalogue.pickle", "rb") as handle:
    halo_catalogue = load(handle)

#%%
halo_masses = (
    sum([getattr(data, f"{x}_mass_in_halo") for x in ["gas", "stellar", "dark_matter"]])
    * 1e10
)  # Particle mass

#%%
def latex_float(f):
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str


#%%
# Make the plot of the halo centers
fig, ax = plt.subplots(figsize=(8, 8))
fig.subplots_adjust(0, 0, 1, 1)
ax.axis("off")

for group_id in group_ids_to_plot:
    halo = halo_catalogue.halos[group_id]

    circle = plt.Circle(
        halo.center[:-1], halo.rvir, color="black", fill=False, linestyle="dashed"
    )
    ax.add_artist(circle)

    # Place the tex tot top right of our halo with some offset.
    diff = halo.rvir * 0.5 ** (0.5)
    manual_offset = 200

    text_to_display = (
        f"Group {group_id}\n"
        f"$M_H = {latex_float(halo_masses[group_id])}$ $h^{{-1}}$M$_\\odot$"
    )

    ax.text(
        diff + manual_offset + halo.center[0],
        diff + manual_offset + halo.center[1],
        text_to_display,
        ha="left",
        va="bottom",
    )

ax.set_xlim(0, 50000)
ax.set_ylim(0, 50000)

fig.savefig("selected_halos_sizes_secondset.pdf")

#%%
