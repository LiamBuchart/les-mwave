"""
    
    Plot vertical cross section of turbulent kinetic energy
    Include topography (add option to include wind barbs)
    
    Required: WRF output file directory (loop through many files)
    
    Output: plot for each time of the vertical cross section
    
    lbuchart@eoas.ubc.ca
    September 26, 2023
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from icecream import ic
from netCDF4 import Dataset
from context import name_dir, script_dir
from file_funcs import (setup_script, get_heights, weighted_mean_std_3d,
                        dist_from_ridge, fmt)

from wrf import (getvar, xy, interp2dxy, interpline, 
                CoordPair, get_cartopy, to_np,
                extract_times)
    
########## 

## USER INPUTS ##
# options ["real", "stable", "neutstab", "logT", "high_step", "low_step"]
exp = "high_step" + "/" # name of the experiment you want to plot
start = (0, 100)
end = (-1, 100)

## END USER INPUTS ##

path, save_path, relevant_files, wrfin = setup_script(exp)

# extract heights of all layers 
heights = getvar(wrfin[0], "height", units="m", meta=False)
ys = heights[:, 0, 0]  

all_tke = []
# loop through files list and make our plots
for ii in range(0, len(wrfin)):
    # import the file in a readable netcdf format
    print(relevant_files[ii])
    ncfile = wrfin[ii]
    
    # get the time in datetime format
    ct = extract_times(ncfile, timeidx=0)
    print(ct)
    
    # winds
    U = getvar(ncfile, "ua", units="m s-1",
               meta=True)
    V = getvar(ncfile, "va", units="m s-1",
               meta=True)
    W = getvar(ncfile, "wa", units="m s-1",
               meta=True)
    
    # means
    mnU = weighted_mean_std_3d(U, heights, all_stats=False)
    mnV = weighted_mean_std_3d(V, heights, all_stats=False)
    mnW = weighted_mean_std_3d(W, heights, all_stats=False)
    
    Up = np.squeeze(np.mean(U - mnU, axis=1))
    Vp = np.squeeze(np.mean(V - mnV, axis=1))
    Wp = np.squeeze(np.mean(W - mnW, axis=1))
    
    # tke
    tke =  0.5 * ((Up**2) + (Vp**2) + (Wp**2))
        
    # interpolation/cross section
    #tke_line = xy(tke, start_point=start, end_point=end)
    #tke_cross = to_np( interp2dxy(tke, tke_line) )
    
    # take mean along the y direction of tke
    tke_cross = tke # np.squeeze( np.mean(tke, axis=1) )
    
    # at 01:15 take a closer look at velocitites and perturbations (should be a rotor????)
    if ii == 50:
        ic(U.coords["bottom_top"].values)
        ic(tke_cross.coords["bottom_top"].values)
        
        f, a = plt.subplots(1, 1, sharex=True, gridspec_kw={"height_ratios":[1]})
        #a[0].imshow(np.mean(U, axis=1), aspect='auto')
        #a[1].imshow(np.mean(V, axis=1), aspect='auto')
        #a[2].imshow(np.mean(W, axis=1), aspect='auto')
        #a[3].imshow(Up, aspect='auto')
        #a[4].imshow(Vp, aspect='auto')
        a.imshow(tke_cross, aspect='auto')
        
        plt.savefig(save_path + "Vel_perturb_" + str(ct)[11:19])
        
        # just plot the xarray variable
        fig, ax = plt.subplots(figsize=(2, 2))
        tke_cross.plot()
        plt.savefig(save_path + "tke_wind" + str(ct)[11:19])
        
    # terrain
    ter = getvar(ncfile, "ter", meta=True)
    ter_cross = interpline(ter, 
                           start_point=CoordPair(start[0], start[1]),
                           end_point=CoordPair(end[0], end[1]))

    ridge_dist = dist_from_ridge(ter_cross)
    
    # get the cartopy projections
    proj = get_cartopy(W)
               
    # make the figure
    fig, axs = plt.subplots(2, 1,
                            sharex=True,
                            gridspec_kw={"height_ratios":[5, 1]})

    tke_contour = axs[0].imshow(tke_cross, cmap='cubehelix', vmin=0, vmax=400, aspect='auto', 
                                extent =[ridge_dist.min(), ridge_dist.max(), ys.max(), ys.min()],
                                interpolation="bicubic")

    ht_fill = axs[1].fill_between(ridge_dist, 0, to_np(ter_cross),
                                  facecolor="saddlebrown")
    
    # hide x labels on all but the bottom 
    for ax in axs:
        ax.label_outer()

    # make pretty titles and whatnot
    axs[0].set_xticks(np.arange(-6000, 40000, 6000))
    axs[0].xaxis.set_tick_params(top=True)
    axs[1].set_xticks(np.arange(-6000, 40000, 6000))

    axs[0].set_ylabel("Height AGL [m]", fontsize=10)
    axs[1].set_ylabel("Terrain Height [m]", fontsize=10)
    axs[1].set_xlabel("Distance from Ridge Top [m]", fontsize=10)
    #axs[0].yaxis.set_label_coords(-.1, .3)

    axs[0].set_yticks(np.arange(0, 11000, 1000))
    axs[0].set_ylim([0, 10000])
    
    axs[1].set_yticks(np.arange(0, 1250, 250))
    axs[1].set_ylim([0, 1200])

    # colorbar 
    #fig.tight_layout()  # call this before calling the colorbar and after calling 
    cbar = plt.colorbar(tke_contour, ax=axs)
    cbar.set_label("Turbulent Kinetic Energy [J/kg]", fontsize=10)
    
    plt.savefig(save_path + "2d_tke_" + str(ct)[11:19])
    plt.close()
    
    # concatenate over all time to get a mean picture
    tt = tke_cross.to_numpy()
    if ii == 5:
        all_tke = tt
    elif ii > 5: 
        all_tke = np.dstack((all_tke, tt))
        
        
# make a mean plot
fig, axs = plt.subplots(2, 1,
                        sharex=True,
                        gridspec_kw={"height_ratios":[5, 1]})

tke_contour = axs[0].imshow(np.mean(all_tke, axis=2), cmap='cubehelix', vmin=0, vmax=300, aspect='auto', 
                            extent =[ridge_dist.min(), ridge_dist.max(), ys.max(), ys.min()],
                            interpolation="bicubic")

ht_fill = axs[1].fill_between(ridge_dist, 0, to_np(ter_cross),
                              facecolor="saddlebrown")
    
# hide x labels on all but the bottom 
for ax in axs:
    ax.label_outer()

# make pretty titles and whatnot
axs[0].set_xticks(np.arange(-6000, 40000, 6000))
axs[0].xaxis.set_tick_params(top=True)
axs[1].set_xticks(np.arange(-6000, 40000, 6000))

axs[0].set_ylabel("Height AGL [m]", fontsize=10)
axs[1].set_ylabel("Terrain Height [m]", fontsize=10)
axs[1].set_xlabel("Distance from Ridge Top [m]", fontsize=10)

axs[0].set_yticks(np.arange(0, 11000, 1000))
axs[0].set_ylim([0, 10000])
    
axs[1].set_yticks(np.arange(0, 1250, 250))
axs[1].set_ylim([0, 1200])

# colorbar 
#fig.tight_layout()  # call this before calling the colorbar and after calling 
cbar = plt.colorbar(tke_contour, ax=axs)
cbar.set_label("Turbulent Kinetic Energy [J/kg]", fontsize=10)
    
plt.savefig(save_path + "2d_tke_exp_mean")
plt.close()
    
print("Complete")    
