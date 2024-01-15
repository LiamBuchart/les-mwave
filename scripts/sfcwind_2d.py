"""
    
    Plot 2d surface winds at 10m AGL 
    Vectors to show direction colour contours for speed
    
    Required: WRF output file directory (loop through many files)
    
    Output: plot for each time of the vertical cross section
    
    lbuchart@eoas.ubc.ca
    February 17, 2023    
    
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import metpy.calc as mpcalc

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

# loop through files to get velocity 
for ii in range(0, len(wrfin)):
    # import the file in a readable netcdf format
    ncfile = wrfin[ii]
    
    # extract the time in datetime format
    ct = extract_times(ncfile, timeidx=0)
    print(ct)
    
    # get the 10m U-wind
    U = getvar(ncfile, "ua", units="m s-1",
               meta=True)
    V = getvar(ncfile, "va", units="m s-1",
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
    U10 = interplevel(U, height, lev,
                      meta=True)
    V10 = interplevel(V, height, lev,
                      meta=True)
    spd10 = np.sqrt(U10**2 + V10**2)
    
    # some values to make nice figures    
    # Get the latitude and longitude points
    lats, lons = latlon_coords(spd10)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(U10) 
    
    # make the figure
    fig, ax = plt.subplots(constrained_layout=True)
    
    # levels to plot the velocity perturbation
    wind_levels = np.arange(0, 10, 1)
    ter_levels = np.arange(0, 1000, 200) 
    
    # Make the contours of terrain and wind perturbation
    ter_lines = plt.contour(ridge_dist, to_np(lats[:, 0]), to_np(ter), 
                            levels=ter_levels, colors="black")
    spd_contour = plt.contourf(ridge_dist, to_np(lats[:, 0]), to_np(spd10),
                             levels=wind_levels, 
                             extend="max", cmap=get_cmap("Reds"))
    
    # classic looking wind barbs
    plt.quiver( to_np(lons[::10,::10]), to_np(lats[::10,::10]),
          to_np(U10[::10, ::10]), to_np(V10[::10, ::10]) )
    
    # make pretty titles and whatnot
    ax.set_xticks(np.arange(-6000, 40000, 6000))
    ax.set_xlabel("Distance from Ridge Top [m]")
    ax.set_xlim([-6000, 40000])
    
    # colorbar 
    cbar = plt.colorbar(to_np(spd_contour), ax=ax,
                        ticks=np.arange(wind_levels[0], wind_levels[-1], 1),
                        orientation="horizontal")
    cbar.set_label("Wind Speed [m/s]", fontsize=10)
    
    plt.savefig(save_path + "spd10_" + str(ct)[11:19])
    plt.close()  
    
    # concatenate over all time to get a mean picture
    SS = spd10.to_numpy()
    UU = U10.to_numpy()
    VV = V10.to_numpy()
    if ii == 5:
        all_spd = SS
        all_UU = UU
        all_VV = VV
    elif ii > 5: 
        all_spd = np.dstack((all_spd, SS))
        all_UU = np.dstack((all_UU, UU))
        all_VV = np.dstack((all_VV, VV))
        
# make the figure
fig, ax = plt.subplots(constrained_layout=True)
    
# levels to plot the velocity perturbation
wind_levels = np.arange(0, 8, 0.5)
ter_levels = np.arange(0, 1000, 200) 
U = np.mean(all_UU, axis=2)
V = np.mean(all_VV, axis=2)
spd = np.mean(all_spd, axis=2)
    
# Make the contours of terrain and wind perturbation
ter_lines = plt.contour(ridge_dist, to_np(lats[:, 0]), to_np(ter), 
                        levels=ter_levels, colors="black")
spd_contour = plt.contourf(ridge_dist, to_np(lats[:, 0]), to_np(spd),
                           levels=wind_levels, 
                           extend="max", cmap=get_cmap("Reds"))
    
# classic looking wind barbs
plt.quiver( to_np(lons[::10,::10]), to_np(lats[::10,::10]),
          to_np(U[::10, ::10]), to_np(V[::10, ::10]) )
    
# make pretty titles and whatnot
ax.set_xticks(np.arange(-6000, 40000, 6000))
ax.set_xlabel("Distance from Ridge Top [m]")
ax.set_xlim([-6000, 40000])
    
# colorbar 
cbar = plt.colorbar(to_np(spd_contour), ax=ax,
                    ticks=np.arange(wind_levels[0], wind_levels[-1], 1),
                    orientation="horizontal")
cbar.set_label("Wind Speed [m/s]", fontsize=10)
    
plt.savefig(save_path + "spd10_mean")
plt.close()  
    
print("Complete")