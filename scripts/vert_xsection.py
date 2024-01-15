"""
    
    Plot vertical cross section of potential temperature and wind speed 
    through the center of the domain.
    Include topography (add option to include wind barbs)
    
    Required: WRF output file directory (loop through many files)
    
    Output: plot for each time of the vertical cross section
    
    lbuchart@eoas.ubc.ca
    August 12, 2022
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
  
from netCDF4 import Dataset
from context import name_dir, script_dir
from file_funcs import (setup_script, get_heights, 
                        dist_from_ridge, fmt)

from wrf import (getvar, xy, interp2dxy, interpline, 
                CoordPair, get_cartopy, to_np,
                extract_times)
    
########## 

## USER INPUTS ##
# options ["real", "stable", "neutstab", "logT", "high_step", "low_step"]
exp = "high_step/" # name of the experiment you want to plot
start = (0, 100)
end = (-1, 100)

## END USER INPUTS ##

path, save_path, relevant_files, wrfin = setup_script(exp)

# extract heights of all layers 
heights = get_heights(wrfin[0])

ys = heights[:, 0, 0]  
ys = ys[1:]

# loop through files list and make our plots
for ii in range(0, len(wrfin)):
    # import the file in a readable netcdf format
    print(relevant_files[ii])
    ncfile = wrfin[ii]
    
    # get the time in datetime format
    ct = extract_times(ncfile, timeidx=0)
    print(ct)
    
    # U winds
    U = getvar(ncfile, "ua", units="m s-1",
               meta=True)
    U_line = xy(U, start_point=start, end_point=end)
    U_cross = interp2dxy(U, U_line)
    print("Wind is in ", np.shape(U))
    
    # potential temperature
    theta = getvar(ncfile, "theta", units="k", 
                   meta=True)
    t_line = xy(theta, start_point=start, end_point=end)
    t_cross = interp2dxy(theta, t_line)
    print("temps are in")
    
    # terrain
    ter = getvar(ncfile, "ter", meta=True)
    ter_cross = interpline(ter, 
                           start_point=CoordPair(start[0], start[1]),
                           end_point=CoordPair(end[0], end[1]))
    print("terrain is in ")

    ridge_dist = dist_from_ridge(ter_cross)
    
    # get the cartopy projections
    proj = get_cartopy(U)


    # make the figure
    fig, axs = plt.subplots(2, 1, 
                            sharex=True,
                            gridspec_kw={"height_ratios":[5, 1]})

    temp_levels = np.arange(270, 370, 2)  # contour lines 
    wind_levels = np.arange(-30, 31, 1)

    xs = np.arange(0, U.shape[-1], 1)
    U_contour = axs[0].contourf(ridge_dist,
                                ys,
                                to_np(U_cross),
                                cmap="seismic",
                                extend="both",
                                levels=wind_levels)

    xs = np.arange(0, theta.shape[-1], 1)
    t_contour = axs[0].contour(ridge_dist, 
                               ys,
                               to_np(t_cross), 
                               colors="k",
                               levels=temp_levels)
    
    axs[0].clabel(t_contour, t_contour.levels[::2], inline=True, fmt=fmt, fontsize=10)

    ht_fill = axs[1].fill_between(ridge_dist, 0, to_np(ter_cross),
                                  facecolor="saddlebrown")
    
    # hide x labels on all but the bottom 
    for ax in axs:
        ax.label_outer()

    # make pretty titles and whatnot
    axs[0].set_xticks(np.arange(-6000, 30000, 6000))
    axs[1].set_xticks(np.arange(-6000, 30000, 6000))

    axs[0].set_ylabel("Height AGL [m]", fontsize=10)
    axs[1].set_ylabel("Terrain Height [m]", fontsize=10)
    #axs[0].yaxis.set_label_coords(-.1, .3)

    axs[0].set_yticks(np.arange(0, 11000, 1000))
    axs[0].set_ylim([0, 5000])
    
    axs[1].set_yticks(np.arange(0, 1250, 250))
    axs[1].set_ylim([0, 1200])

    # colorbar 
    fig.tight_layout()  # call this before calling the colorbar and after calling 
    cbar = plt.colorbar(U_contour, ax=axs)
    cbar.set_label("Wind Speed [m/s]", fontsize=10)
    
    plt.savefig(save_path + "potential_speed_" + str(ct)[11:19])
    plt.close()
    
    # concatenate over all time to get a mean picture
    UU = U_cross.to_numpy()
    t = t_cross.to_numpy()
    if ii == 5:
        all_spd = UU
        all_pot = t
    elif ii > 5: 
        all_spd = np.dstack((all_spd, UU))
        all_pot = np.dstack((all_pot, t))

        
# make the mean figure
fig, axs = plt.subplots(2, 1, 
                        sharex=True,
                        gridspec_kw={"height_ratios":[5, 1]})

temp_levels = np.arange(270, 370, 2)  # contour lines 
wind_levels = np.arange(-15, 26, 1)

xs = np.arange(0, U.shape[-1], 1)
U_contour = axs[0].contourf(ridge_dist,
                            ys,
                            np.mean(all_spd, axis=2),
                            cmap="seismic",
                            extend="both",
                            levels=wind_levels)

xs = np.arange(0, theta.shape[-1], 1)
t_contour = axs[0].contour(ridge_dist, 
                           ys,
                           np.mean(all_pot, axis=2), 
                           colors="k",
                           levels=temp_levels)
    
axs[0].clabel(t_contour, t_contour.levels[::2], inline=True, fmt=fmt, fontsize=10)

ht_fill = axs[1].fill_between(ridge_dist, 0, to_np(ter_cross),
                              facecolor="saddlebrown")
    
# hide x labels on all but the bottom 
for ax in axs:
    ax.label_outer()

# make pretty titles and whatnot
axs[0].set_xticks(np.arange(-6000, 42000, 6000))
axs[1].set_xticks(np.arange(-6000, 42000, 6000))

axs[0].set_ylabel("Height AGL [m]", fontsize=10)
axs[1].set_ylabel("Terrain Height [m]", fontsize=10)

axs[0].set_yticks(np.arange(0, 11000, 1000))
axs[0].set_ylim([0, 5000])
    
axs[1].set_yticks(np.arange(0, 1250, 250))
axs[1].set_ylim([0, 1200])

# colorbar 
fig.tight_layout()  # call this before calling the colorbar and after calling 
cbar = plt.colorbar(U_contour, ax=axs)
cbar.set_label("Wind Speed [m/s]", fontsize=10)
    
plt.savefig(save_path + "potential_speed_mean")
plt.close()
    
print("Complete")    
