#%% [markdown]
# # Finding 'Interesting' Halos
# We want to find some interesting halos that have high fractions
# of transfer and outside material

#%%
from ltcaesar import read_data_from_file
import numpy as np
import matplotlib.pyplot as plt

#%%
# Setup our favourite stylesheet
plt.style.use("mnras_flatiron")

#%%
directory = "s50j7kAHF"
data = read_data_from_file(f"{directory}/lt/lt_outputs.hdf5")

#%% [markdown]
# We now need to calculate three things for each halo:
# + The fraction of baryonic mass from outside the halo
# + The fraction of baryonic mass from other halos
# + The fraction of baryonic mass from our own LR

#%%
def grab(x):
    return getattr(data, f"gas_{x}") + getattr(data, f"stellar_{x}")


baryonic_mass = grab("mass_in_halo")
from_outside = grab("mass_in_halo_from_outside_lagrangian")
from_other = grab("mass_in_halo_from_other_lagrangian")
from_own = grab("mass_in_halo_from_lagrangian")

#%%
# Now need to be careful if there are halos without baryons
mask = baryonic_mass != 0
ratio_outside = from_outside[mask] / baryonic_mass[mask]
ratio_other = from_other[mask] / baryonic_mass[mask]
ratio_own = from_own[mask] / baryonic_mass[mask]
group_ids = np.arange(len(from_own))[mask]


#%%
# Let's just make a test plot
for ratio, name in zip(
    [ratio_own, ratio_other, ratio_outside], ["From Own", "From Other", "From Outside"]
):
    plt.hist(ratio, label=name, histtype="step", density=True)

plt.xlabel("Ratio")
plt.ylabel("Fraction")
plt.ylim(0, 5)
plt.xlim(0, 1)
plt.legend()
plt.tight_layout()

#%% [markdown]
# We want to search for 'interesting' things within a halo mass
# range, let's for now say $10^{12} \rightarrow 10^{13}$

#%%
halo_mass_range = [1e12, 1e13]
unit_mass = 1e10
internal_mass_range = [x / unit_mass for x in halo_mass_range]

#%%
halo_masses = baryonic_mass + data.dark_matter_mass_in_halo
mass_mask = np.logical_and(
    halo_masses > internal_mass_range[0], halo_masses < internal_mass_range[1]
)
# Also need to include the nonzero baryonic mass condition
mass_mask = mass_mask[baryonic_mass != 0.0]

#%%
print(f"Found {mask_mass.sum():d} halos")

#%%
# Let's make a plot for these halos
for ratio, name in zip(
    [ratio_own, ratio_other, ratio_outside], ["From Own", "From Other", "From Outside"]
):
    plt.hist(ratio[mass_mask], label=name, histtype="step", density=True)

plt.title(f"For halo mass {halo_mass_range[0]:e} to {halo_mass_range[1]:e}")
plt.xlabel("Ratio")
plt.ylabel("Fraction")
plt.ylim(0, 5)
plt.xlim(0, 1)
plt.legend()
plt.tight_layout()

#%% [markdown]
# Let's choose some interesting halos, and some normal ones. We'll choose one
# with very high transfer from outside, one with very high transfer from others,
# and one with the median fraction from its own LR.

#%%
max_other = group_ids[mass_mask][ratio_other[mass_mask].argmax()]
max_outside = group_ids[mass_mask][ratio_outside[mass_mask].argmax()]
# Can't use .median() because we want a value that actually belongs to dataset
median_of_own = np.argsort(ratio_own[mass_mask])[np.sum(mass_mask) // 2]
median_own = group_ids[mass_mask][median_of_own]

#%%
print("These are halos:")
print(f"Max from other, id = {max_other}, mass = {halo_masses[max_other]}")
print(f"Max from outside, id = {max_outside}, mass = {halo_masses[max_outside]}")
print(f"Median from own, id = {median_own}, mass = {halo_masses[median_own]}")

#%% [markdown]
# Okay, great! We've found three really interesting halos. Let's take a look at
# how their values match up to each other.

#%%
ids = [max_other, max_outside, median_own]
names = ["Maximal Other", "Maximal Outside", "Median Own"]

fig, ax = plt.subplots(1, 3, sharey=True, figsize=(4, 2))

for id, name, a in zip(ids, names, ax):
    a.bar(
        [0.5, 1.5, 2.5],
        [ratio_own[id], ratio_other[id], ratio_outside[id]],
        width=1.0,
        color=["C0", "C1", "C2"],
    )
    a.set_xticks([0.5, 1.5, 2.5])
    a.set_xticklabels(["Own", "Other", "Outside"])
    a.set_title(name)
    a.set_xlim(0, 3)

fig.tight_layout()


#%%
print("IDs to use in the other code:")
print(ids)


#%%
