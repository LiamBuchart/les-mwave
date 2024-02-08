"""
    
    Plot the moisture flux divergence in a control volume over the crest of the slope and to the lee
    of the ridge
    
    lbuchart@eoas.ubc.ca
    January 15, 2024

"""

from context import name_dir, script_dir, json_dir
from icecream import ic
from file_funcs import check_folder, setup_script, cb_color_palete   

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

# get experiment names
# get grid resolution from our json dictionary
with open(str(json_dir) + "config.json") as f:
    config = json.load(f)
exps = config["exps"]["names"]
ttd = config["timing"]["ttd"]  # exttract the number of minutes between output
ic(exps)

moistfile = "leeslope_moisture_flux.csv"
compfile = "leeslope_moisture_component_"
components = ["x-moist", "y-moist", "z-moist"]

moist = pd.read_csv(script_dir + moistfile)

# extract times
Times = moist["Times"].tolist() 
for ii in range(len(Times)):
    Times[ii] = Times[ii][11:]

# make a mins since start timing
mss = []
for ii in range(len(Times)):
    mss.append(ii * ttd)
    
# add to dataframe
moist["MSS"] = mss

# folder to save figures (check its there, if not make)
fig_folder = check_folder(script_dir, "figures/moistureflux/")

# dictionary of experiments and the desired plot colors
colors = cb_color_palete( exps )

# just get column names of experiments
cols = moist.columns.values.tolist()[1:-2]
# drop the final Times column
moist = moist.drop(columns=["MSS", "Times"])

# normalize (min max operator)
ic(moist.head())
dfmin = min(moist.min())
dfmax = max(moist.max())
moist_full = pd.DataFrame(columns=exps)
for ex in exps:
    vals = moist[ex].to_numpy()
    diff = dfmax - dfmin
    moist_full[ex] = np.subtract(vals, dfmin) / (dfmax - dfmin)
ic(moist_full.head())

# add back to the to dataframe
moist_full["MSS"] = mss

# plot the values first individually then as subplots
moist_full.plot(x="MSS", y=cols, kind="line",
        color=colors)

plt.title(r'Control Volume Normalized Moisture Divergence')
plt.ylabel(r'Moisture Divergence [$gkm^{-1}ms^{-1}$]')

plt.xlim([0, 240])
plt.xticks(mss[0::3])

# save it
plt.savefig(str(fig_folder) + "AllMoistDivg")

# now plot the main components for each experiment
# dictionary of experiments and the desired plot colors
colors = cb_color_palete( components )
colors = list( colors.values() )

# create a list of latex names for each component
names = [r"$\frac{\partial q u}{\partial x}$", 
         r"$\frac{\partial q v}{\partial y}$", 
         r"$\frac{\partial q w}{\partial z}$"]

# grab dataframe which is single column
df_full = pd.read_csv(script_dir + "leeslope_moistflux_all_exps_comps.csv")
    
# normalize (min max operator)
dfmin = df_full["value"].min() 
dfmax = df_full["value"].max() 
ic(dfmin, dfmax)
df_full = df_full.assign( Norm = lambda x: ((x["value"] - dfmin) / ( dfmax - dfmin)) )
 
# set up the figure
fig, axs = plt.subplots(len(exps), 1, 
                         sharex=True,
                         gridspec_kw={"height_ratios":np.ones(len(exps))})

fig.suptitle("Normalized Moisture Divergence [ ]", fontsize=24)

# loop through all experiments and components to create plots
line = 0
for ex in exps: 
    df_ex = df_full.loc[ df_full["exp"] == ex]
    
    # loop through components   
    ii = 0 
    for comp in components:
        df_plot = df_ex.loc[ df_ex["component"] == comp]
        
        axs[line].plot(mss, df_plot["Norm"].tolist(), 
                       color=colors[ii],
                       label=names[ii])
        
        ii += 1
              
    axs[line].legend()
    axs[line].set_ylim([0, 1])
    # empty y axis label (instead put it on plot)
    axs[line].set_ylabel(" ")

    # place a text box bottom right corner
    axs[line].text(0.05, 0.45, ex, transform=axs[line].transAxes, fontsize=16,
                   verticalalignment='top') #, bbox=props)
    
    axs[line].set_xlim([0, 240])
    
    axs[line].get_legend().remove()
        
    line += 1

# label only final x-axis
axs[-1].set_xlabel("Minutes Since Start", fontsize=20)

# Shrink current axis by 20%
for ax in axs:
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.95, box.height])

handles, labels = axs[2].get_legend_handles_labels()
fig.legend(handles, labels, 
           loc="center right", bbox_to_anchor=(1, 0.5))
#plt.tight_layout()

# save    
plt.savefig(str(fig_folder) + "componentflux")

print("Complete")