"""
    
    Plots from Energy.csv (produced by tke_all)
    1. Mean energy density (with 25th and 75th percenile wiskers, shaded std)
    2. Total energy density
    3. Ep and Ek for each experiment
    
"""
    
from context import name_dir, script_dir, json_dir
from icecream import ic
from file_funcs import check_folder, setup_script    

import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# get experiment names
# get grid resolution from our json dictionary
with open(str(json_dir) + "config.json") as f:
    config = json.load(f)
exps = config["exps"]["names"]
ic(exps)

exdf = pd.DataFrame()

fig, axs = plt.subplots(int(len(exps)), 1, 
                        sharex=True,
                        gridspec_kw={"height_ratios":[1, 1, 1, 1, 1, 1]} )

# folder to save figures (check its there, if not make)
fig_folder = check_folder(script_dir, "figures/EnerDen")

# plot the energy density in indivdual sub plots
ii = 0
for ex in exps:
    
    dffile = "Energy_" + ex + ".csv"

    df = pd.read_csv(script_dir + dffile)
    print(df.head())
    cols = df.columns.values.tolist()
    cols = cols[1:]

    times = df[cols[-1]].to_numpy()
    Edmn = df[cols[0]].to_numpy()
    Estd = df[cols[2]].to_numpy()
    E25 = df[cols[3]].to_numpy()
    E75 = df[cols[5]].to_numpy()
    
    tt = range(0, len(df[cols[-1]]))
    
    # plot single experiment Energy density with (ensure a consistent axis)
    # have all plots on a common figure
    fig.suptitle("Mean Energy Density [J/kg]")
    
    axs[ii].plot(tt, df[cols[0]], 
                 'r-', label="Mean")
    axs[ii].plot(tt, df[cols[3]], 
                 'k-')
    axs[ii].plot(tt, df[cols[5]], 
                 'k-')
    
    axs[ii].fill_between(tt, df[cols[0]]-df[cols[2]], df[cols[0]]+df[cols[2]], 
                         alpha=0.4, facecolor='g', label="50% Interval")
    
    
    #axs[ii].set_ylim([0, 2000])
    
    axs[ii].legend()
    
    plt.savefig(str(fig_folder) + "/Mean_EnergyDensities")
    
    ii += 1


# plot just turbulent kinetic energy
fig, axs = plt.subplots(1, 1)

for ex in exps:
    
    dffile = "Energy_" + ex + ".csv"

    df = pd.read_csv(script_dir + dffile)
    cols = df.columns.values.tolist()
    cols = cols[1:]
    
    times = df[cols[-1]].to_numpy()
    Ek = df[cols[-3]].to_numpy()
    
    fig.suptitle("Mean Kinetic Energy [J/kg]")
    axs.plot(times, Ek, label=ex)
    
    axs.set_xticks(times)
    axs.set_xlabel( list(times[0::10]), rotation=45 )
    plt.setp(axs.get_xticklabels()[::2], visible=False)
    
    axs.legend(loc="lower center", ncol=len(exps), fancybox=True)
    
    plt.savefig(str(fig_folder) + "/Mean_KineticEnergy")


print("Complete")
    