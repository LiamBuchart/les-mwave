"""

Plot a vertical cross section of the U velocity with the streamlines

lbuchart@eoas.ubc.ca
August 19, 2021

"""

import numpy as np
import matplotlib.pyplot as plt
import metpy.calc as mpcalc

from file_funcs import (setup_script, dist_from_ridge,
                        get_heights) 

from wrf import (xy, CoordPair, interpline, vinterp, interplevel,
                 extract_times, getvar, interp2dxy, to_np)

###########

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
    ncfile = wrfin[ii]
    
    # get the time in datetime format
    ct = extract_times(ncfile, timeidx=0)
    print(ct)
    
    # U winds
    U = getvar(ncfile, "ua", units="m s-1",
               meta=True)
    U_line = xy(U, start_point=start, end_point=end)
    
    # interpolate velocities to evenly spaced heights for streamlines
    new_heights = np.linspace(ys[0], ys[-1], 2*len(ys))
    
    # interpolate winds to these heights
    interp_levels = new_heights / 1000  # convert to kms
    # heights
    ht = getvar(ncfile, "z", units="m")
    
    # interpolate to evenly spaced heights for streamlines
    U_lev = interplevel(U, ht, new_heights, missing=-9999)
    print(np.shape(U_lev))

    # cross section
    U_cross = interp2dxy(U_lev, U_line)
    print("Wind is in")
    
    W = getvar(ncfile, "wa", units="m s-1",
               meta=True)
    W_line = xy(W, start_point=start, end_point=end)
    
    # interpolate vertical velocity in same manner
    W_lev = interplevel(W, ht, new_heights, missing=-9999)
    
    W_cross = interp2dxy(W_lev, W_line)
    print("Vert Vel is in")
    
    # terrain
    ter = getvar(ncfile, "ter", meta=True)
    ter_cross = interpline(ter, 
                           start_point=CoordPair(start[0], start[1]),
                           end_point=CoordPair(end[0], end[1]))
    print("terrain is in ")

    ridge_dist = dist_from_ridge(ter_cross)
    
    # make the figure
    fig, ax = plt.subplots(constrained_layout=True)
    # levels to plot winds on
    wind_levels = np.arange(-30, 31, 1)

    U_contour = ax.contourf(ridge_dist,
                                new_heights,
                                to_np(U_cross),
                                cmap="seismic",
                                extend="both",
                                levels=wind_levels)
    
    X, Y = np.meshgrid(ridge_dist, new_heights)

    streamline = ax.streamplot(ridge_dist, 
                                   new_heights, 
                                   U_cross, 
                                   W_cross,
                                   color="k",
                                   density=0.9)
    
    ht_fill = ax.fill_between(ridge_dist, 0, to_np(ter_cross),
                                  facecolor="saddlebrown")

    # make pretty titles and whatnot
    ax.set_xticks(np.arange(-6000, 46000, 6000))

    ax.set_ylabel("Height ASL [m]", fontsize=10)
    
    ax.set_yticks(np.arange(0, 11000, 500))
    ax.set_ylim([0, 5000])

    # colorbar 
    cbar = plt.colorbar(U_contour, ax=ax)
    cbar.set_label("Wind Speed [m/s]", fontsize=10)
    
    plt.savefig(save_path + "streamlines_" + str(ct)[11:19])
    plt.close()
    
    print("The max U wind is: ", np.max(U_cross), np.min(U_cross))
    
    # clear memory space by deleting large variables
    del U
    del U_cross
    del W
    del W_cross
