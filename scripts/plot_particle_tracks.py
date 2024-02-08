"""
    
    Plot the output from the particle_tracking.py
    Input is two csv files x[z]_particle_trajectory_positions.csv
    
    February 2, 2024
    lbuchart@eoas.ubc.ca
    
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json

from file_funcs import cb_color_palete, check_folder
from context import script_dir, json_dir
from icecream import ic

xdf = pd.read_csv(script_dir + "x_particle_trajectory_positions.csv")
zdf = pd.read_csv(script_dir + "z_particle_trajectory_positions.csv")

edges = [0, 9500, 0, 44000]  # vertical and horizontal grid edges

# get experiment names from our json dictionary
with open(str(json_dir) + "config.json") as f:
    config = json.load(f)
exps = config["exps"]["names"]

# get height of release and x position from json file
xpos = int(config["trajectories"]["xpos"])
heights = config["trajectories"]["zpos"]

# dictionary of experiments and the desired plot colors
colors = cb_color_palete( exps )

# folder to save figures (check its there, if not make)
fig_folder = check_folder(script_dir, "figures/trajectory")
# release points
with open(str(json_dir) + "config.json") as f:
    config = json.load(f)
# get height of release and x position from json file
xpos = int(config["trajectories"]["xpos"])
heights = config["trajectories"]["zpos"]

for hh in heights:
    particles = 0
    
    # initialize a figure
    fig, ax = plt.subplots()
    
    plt.xlim([edges[2], edges[3]])
    #plt.ylim([edges[0], edges[1]])
    
    plt.title("Released at: " + str(hh) + "[m] AGL, " + str(xpos) + " east of grid edge", fontsize=20)
    ax.set_ylabel("Height [m] AGL", fontsize=20)
    ax.set_xlabel("Distance from Western Grid Edge [m]", fontsize=20)
    
    for exp in exps:
        ic(hh, exp)
        
        xx = xdf[exp + "_" + str(hh)]
        zz = zdf[exp + "_" + str(hh)]
        
        ic(len(xx), len(zz))
               
        plt.plot(xx, zz, color=colors[exp], label=exp)
        
        particles += 1
        
    plt.legend(fontsize=16)
    
    plt.savefig(str(fig_folder) + "/trajectories_height_" + str(hh))
    print("Next Release Point")
    
print("Complete")