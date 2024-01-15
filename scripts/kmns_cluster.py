"""
    
    Clustering on multiple datasets realted to the mwave experiments
    1. maximum canopy height winds vs max vertical velocity
    2. "all variables" including mean tke, max temp perturbation, max wind, max vert vel, momentum
    3. on variables in the hot dry wind index
    
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json

from context import script_dir, json_dir
from file_funcs import (cb_color_palete, check_folder)
from icecream import ic

from wrf import (getvar, xy, interp2dxy, interpline, 
                CoordPair, get_cartopy, to_np,
                extract_times)

from sklearn.cluster import KMeans

# first make a scatter plot of the max canopy velocities

can_file = "Canopy_Max_Velocities.csv"
mc_df = pd.read_csv(script_dir + can_file)
ic(mc_df.head())

vert_file = "Max_Vert_Velocity.csv"
wm_df = pd.read_csv(script_dir + vert_file)
ic(wm_df.head())

# remove first few rows as spin up (first 20 mins)
mc_df = mc_df.iloc[5:]
wm_df = wm_df.iloc[5:]

# get experiment names from our json dictionary
with open(str(json_dir) + "config.json") as f:
    config = json.load(f)
exps = config["exps"]["names"]

# folder to save figures (check its there, if not make)
fig_folder = check_folder(script_dir, "figures/Kmeans/")

# move all x and y values to their own lists
can = []
vert = []

# dictionary of experiments and the desired plot colors
colors = cb_color_palete( exps )
ic(colors)
ii = 0

plt.title("All Max Canopy and Vertical Velocities")
plt.ylabel("Max Vertical Vel [m/s]")
plt.xlabel("Max Canopy Height Wind [m/s]")
for ex in exps:
    # plot all values 10m vel (x-axis), max vert vel (y_axis)
    # as a single scatter plot
    x = mc_df[ex].to_list()
    y = wm_df[ex].to_list()
    
    
    plt.scatter(x, y, color=colors[ex], label=ex, s=100)
    
    ii += 1
    
    for jj in range(len(x)):
    
        can.append(x[jj])
        vert.append(y[jj])
    
plt.legend()
plt.savefig(fig_folder + "all_scatter")

# start our clustering check (elbow method)
data = list(zip(can, vert))
inertias = []

clust_nn = 11  # number of clusters
for ii in range(1, clust_nn):
    kmeans = KMeans(n_clusters=ii)
    kmeans.fit(data)
    inertias.append(kmeans.inertia_)

fig = plt.figure()    
plt.plot(range(1, clust_nn), inertias, marker='o')
plt.title("Elbow Method")
plt.xlabel("Number of Clusters")
plt.ylabel("Inertia")

plt.savefig(fig_folder + "inertias")

optim_clust = 3
kmeans = KMeans(n_clusters=optim_clust)
kmeans.fit(data)

fig = plt.figure()
plt.title("Max Velocities K-Means: " + str(optim_clust) + " Clusters")
plt.ylabel("Maximum Downward Vert. Vel. [m/s]")
plt.xlabel("Maximum Canopy Height Velocity [m/s]")
plt.scatter(can, vert, c=kmeans.labels_)

plt.savefig(fig_folder + "optimal_clustering")

print("Complete")