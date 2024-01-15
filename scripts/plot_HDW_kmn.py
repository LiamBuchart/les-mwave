"""
    
    Plot the Hot Dry Windy Index Timeseries
    Plots from HWD_" ".py (produced by hot_dry_windy.py)
    Do a kmeans clustering on HWD components
    
    
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

fig, axs = plt.subplots(1, 1,
                        gridspec_kw={"height_ratios":[1]} )

# folder to save figures (check its there, if not make)
fig_folder = check_folder(script_dir, "figures/HDW")

wind_file = "HDW_max_winds.csv"
vpd_file = "HDW_vapour_pressure_deficit.csv"

wdf = pd.read_csv(script_dir + wind_file)
vpd_df = pd.read_csv(script_dir + vpd_file)

# remove rows that will count as spin-up (first 20 mins)
wdf = wdf.iloc[5:]
vpd_df = vpd_df.iloc[5:]

cols = wdf.columns.values.tolist()[1:-1]
ic(cols)

# calculate the HDW place in single DataFrame
df = pd.DataFrame( )
for ex in exps:
    ls = wdf[ex].tolist()
    lvpd = vpd_df[ex].tolist()
    hdw = [a*b for a, b in zip(ls, lvpd)]
    df[ex] = hdw

# extract times
Times = wdf["Times"].tolist() 
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

#plot
ax = df.plot(x="MSS", y=cols, kind="line",
        color=colors,
        figsize=(12, 12))

plt.legend(loc="upper left", fontsize=20, frameon=False)

plt.title("Hot Dry Windy Index", fontsize=30)

#plt.ylim([0, 25])
plt.yticks(fontsize=15)
plt.ylabel("HDW [ ]", fontsize=20)

plt.xlim([0, max(mss)])
plt.xticks(fontsize=15)
ax.set_xticks(mss[5:])
#ax.set_xticklabels(mss, **hfont)
plt.xlabel("Minutes since start [mins]", fontsize=20)
# now set most labels as invisible
for label in ax.xaxis.get_ticklabels()[1::2]:
    label.set_visible(False)

plt.savefig(str(fig_folder) + "/HDW_Timeseries" )

## Do k-means clustering on the components of the HDW
from sklearn.cluster import KMeans
# move all x and y values to their own lists
mw = []
vpd = []

# dictionary of experiments and the desired plot colors
colors = cb_color_palete( exps )
ic(colors)
ii = 0

fig = plt.figure() 
plt.title("HWD_Components")
plt.ylabel("Max Vapor Pressure Deficit [hPa]")
plt.xlabel("Max Canopy Height Wind [m/s]")
for ex in exps:
    # plot all values 10m vel (x-axis), max vert vel (y_axis)
    # as a single scatter plot
    x = wdf[ex].to_list()
    y = vpd_df[ex].to_list()
     
    plt.scatter(x, y, color=colors[ex], label=ex, s=100)
    
    ii += 1
    
    for jj in range(len(x)):
    
        mw.append(x[jj])
        vpd.append(y[jj])
    
plt.legend()
plt.savefig(str(fig_folder) + "/all_scatter")

# start our clustering check (elbow method)
data = list(zip(mw, vpd))
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

plt.savefig(fig_folder + "/inertias")

optim_clust = 2
kmeans = KMeans(n_clusters=optim_clust)
kmeans.fit(data)

fig = plt.figure()
plt.title("Hot Dry Wind Indec K-Means: " + str(optim_clust) + " Clusters")
plt.ylabel("Max. Vapor Pressure Deficit [hPa]")
plt.xlabel("Max Winds [m/s]")
plt.scatter(mw, vpd, c=kmeans.labels_)

plt.savefig(fig_folder + "/optimal_clustering")


print("Complete")