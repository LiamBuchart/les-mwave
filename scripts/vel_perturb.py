"""
    
Plot the velocity perturbation onlong the 10m surface. Used to judge canopy height
wind speed pertubations - useful for tree heights in the southern half of BC and Alberta

lbuchart@eoas.ubc.ca
September 6, 2022    
    
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import metpy.calc as mpcalc
import json

from file_funcs import (setup_script, dist_from_ridge)

from wrf import (interpline, extract_times, get_cartopy,
                 getvar, to_np, latlon_coords, CoordPair,
                 interplevel, cartopy_xlim, cartopy_ylim)

##########

## USER INPUTS ##
# options ["real", "stable", "neutstab", "logT", "high_step", "low_step"]
exp = "high_step/"  # name of the experiment you are plotting
lev = 10  # height in m that you want wind velocity perturbation

start = (0, 100)
end = (-1, 100) 

## END USER INPUTS ##

path, save_path, relevant_files, wrfin = setup_script(exp=exp)

# loop through files to get velocity perturbation and plot
# plot alongside velocity??
for ii in range(0, len(wrfin)):
    # import the file in a readable netcdf format
    ncfile = wrfin[ii]
    
    # extract the time in datetime format
    ct = extract_times(ncfile, timeidx=0)
    print(ct)
    
    # get the 10m U-wind
    U = getvar(ncfile, "ua", units="m s-1",
               meta=True)
    
    height = getvar(ncfile, "height_agl", units="m",
                    meta=True)

    # get the terrain height
    ter = getvar(ncfile, "ter", meta=True)
    ter_cross = interpline(ter, 
                           start_point=CoordPair(start[0], start[1]),
                           end_point=CoordPair(end[0], end[1]))
    ridge_dist = dist_from_ridge(ter_cross)
    
    # interpolate to our height [m]
    spd = interplevel(U, height, lev,
                  meta=True)
    
    # calculate the velocity perturbation
    #print("Size: ", np.shape(spd))
    uus = np.mean(spd[:, 0])  # get an unperturbed upstream velocity values (taken from grid edge)
    u_perturb = (spd - uus) / uus
    
    # some values to make nice figures    
    # Get the latitude and longitude points
    lats, lons = latlon_coords(spd)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(spd)    


    # make the figure
    fig, ax = plt.subplots(constrained_layout=True)
    
    # levels to plot the velocity perturbation
    wind_levels = np.arange(-5, 5.1, 0.1)
    ter_levels = np.arange(0, 1000, 250)
    
    # Make the contours of terrain and wind perturbation
    ter_lines = plt.contour(ridge_dist, to_np(lats[:, 0]), to_np(ter), 
                            levels=ter_levels, colors="black")
    u_contour = plt.contourf(ridge_dist, to_np(lats[:, 0]), to_np(u_perturb),
                             levels=wind_levels, 
                             extend="both", cmap=get_cmap("seismic"))
    
    # make pretty titles and whatnot
    ax.set_xticks(np.arange(-6000, 42000, 6000))
    ax.set_xlabel("Distance from Ridge Top [m]")
    
    # colorbar 
    cbar = plt.colorbar(to_np(u_contour), ax=ax,
                        ticks=np.arange(wind_levels[0], wind_levels[-1], 1),
                        orientation="horizontal")
    cbar.set_label("Wind Perturbation []", fontsize=10)

    #ax.gridlines()
    
    plt.savefig(save_path + "perturb_u_" + str(ct)[11:19])
    plt.close()  
    
print("Complete")
