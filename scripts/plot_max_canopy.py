"""
    
    Plots from Energy.csv (produced by max_canopy.py)
    
    
"""
    
from context import name_dir, script_dir, json_dir
from icecream import ic
from file_funcs import check_folder, setup_script, cb_color_palete    

import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# get experiment names
# get grid resolution from our json dictionary
with open(str(json_dir) + "config.json") as f:
    config = json.load(f)
exps = config["exps"]["names"]
ttd = config["timing"]["ttd"]  # exttract the number of minutes between output

exdf = pd.DataFrame()

fig, axs = plt.subplots(1, 1,
                        gridspec_kw={"height_ratios":[1]} )

# folder to save figures (check its there, if not make)
fig_folder = check_folder(script_dir, "figures/Max_Vels")

# plot max winds of each experiment
    
dffile = "Canopy_Max_Velocities.csv"

df = pd.read_csv(script_dir + dffile)
print(df.head())

cols = df.columns.values.tolist()[1:-1]
ic(cols)

# extract times
Times = df["Times"].tolist() 
for ii in range(len(Times)):
    Times[ii] = Times[ii][11:]

# make a mins since start timing
mss = []
for ii in range(len(Times)):
    mss.append(ii * ttd)
    
# add to dataframe
df["MSS"] = mss

# dictionary of experiments and the desired plot colors
colors = cb_color_palete( exps )

# consistent font for plotting
hfont = {'fontname':'monospace'}

#plot
ax = df.plot(x="MSS", y=cols, kind="line",
        color=colors,
        figsize=(12, 12))

plt.legend(loc="upper left", fontsize=20, frameon=False)

plt.title("Max Canopy Height Wind Speed", fontsize=30, **hfont)

plt.ylim([0, 10])
plt.yticks(fontsize=15)
plt.ylabel("Wind Speed [m/s]", fontsize=20, **hfont)

plt.xlim([0, max(mss)])
ax.set_xticks(mss)
#ax.set_xticklabels(mss, **hfont)
plt.xlabel("Minutes since start [mins]", fontsize=20, **hfont)
# now set most labels as invisible
for label in ax.xaxis.get_ticklabels()[1::2]:
    label.set_visible(False)

plt.xticks(mss[0::4], fontsize=15)
plt.savefig(str(fig_folder) + "/Max_CanopyWind" )

print("Complete")


    