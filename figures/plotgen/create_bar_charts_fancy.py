#%% [markdown]
# We need to create bar charts for our fancy plot to show the fraction of stuff from
# each region. We'll do that here.

#%%
group_ids_to_plot = [0, 431, 88, 299, 9]

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
ratio_own = from_own / baryonic_mass
ratio_other = from_other / baryonic_mass
ratio_outside = from_outside / baryonic_mass

#%%
import os

os.mkdir("barcharts")


#%%
# Actually make the bar charts one by one
for id in group_ids_to_plot:
    fig, a = plt.subplots(figsize=(1, 1))

    a.bar(
        [0.5, 1.5, 2.5],
        [ratio_own[id], ratio_other[id], ratio_outside[id]],
        width=1.0,
        color=["C0", "C1", "C2"],
    )

    a.set_xticks([0.5, 1.5, 2.5])
    a.set_xticklabels(["Own", "Other", "Outside"])
    a.set_xlim(0, 3)
    a.set_ylim(0, 0.8)

    fig.tight_layout()

    fig.savefig("barcharts/group_{}.pdf".format(id))

#%%
