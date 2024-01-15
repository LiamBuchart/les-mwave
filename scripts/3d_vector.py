"""
    
    Plot 3d winds with quiver of 3 levels of winds
    10m, 50m, 250m
    Note: Attempting to add topography underlay
    
    Required: WRF output file directory (loop through many files)
    
    Output: plot for each output time step
    
    lbuchart@eoas.ubc.ca
    February 25, 2023     
    
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import metpy.calc as mpcalc
from context import json_dir

from file_funcs import (setup_script, dist_from_ridge, vert_dist_center)

from wrf import (interpline, extract_times, get_cartopy,
                 getvar, to_np, latlon_coords, CoordPair,
                 interplevel, cartopy_xlim, cartopy_ylim)

## USER INPUTS ##
exp = "high_step/"  # name of the experiment you are plotting
levs = [10, 50, 250] # height in m that you want wind velocity perturbation

start = (0, 100)
end = (-1, 100) 

## END USER INPUTS ##

path, save_path, relevant_files, wrfin = setup_script(exp=exp)

# get grid info from the json dictionary
with open(str(json_dir) + "config.json") as f:
    config = json.load(f)
dy = config["grid_dimensions"]["dy"]
dx = config["grid_dimensions"]["dx"]

# get the terrain heights
ter = getvar(wrfin[0], "ter", meta=True)
ter_cross = interpline(ter, 
                start_point=CoordPair(start[0], start[1]),
                end_point=CoordPair(end[0], end[1]))
ridge_dist = dist_from_ridge(ter_cross)

ter_y = interpline(ter, 
                   start_point=CoordPair(0, 0),
                   end_point=CoordPair(0, -1))
dy_dist = vert_dist_center(ter_y)

# create grid to plot on  - now configure lines to fit on this grid
xd, yd = np.meshgrid(np.linspace(ridge_dist[0], ridge_dist[-1], dx), 
                     np.linspace(dy_dist[0], dy_dist[-1], dy))

height = getvar(wrfin[0], "height_agl", units="m",
                meta=True)

# loop through files to get velocities
for ii in range(0, 5): #len(wrfin)):
    # import the file in a readable netcdf format
    ncfile = wrfin[ii]
    
    # extract the time in datetime format
    ct = extract_times(ncfile, timeidx=0)
    print(ct)
    
    # get desired winds
    U, V = getvar(ncfile, "ua", units="m s-1",
                  meta=True), getvar(ncfile, "va", units="m s-1",
                  meta=True)
    
    # interpolate to our height [m]
    # add loop here if you want more heights
    #Ulev = interplevel(U, height, levs[1],
    #                meta=True)
    print("The shape of U is: ", np.shape(U))
    Ulev = U[1, :, :]  # take second model level of wind
    #Vlev = interplevel(V, height, levs[1],
    #                meta=True)
    Vlev = V[1, :, :] # take second model level of wind
    
    # plot in 2d to get streamline properties
    fig_tmp, ax_tmp = plt.subplots()
    #print(np.shape(ridge_dist), np.shape(dy_dist), np.shape(Ulev), np.shape(Vlev))
    res = ax_tmp.streamplot(ridge_dist.T, dy_dist.T, Ulev, Vlev, color='k')
    plt.savefig(save_path + "streamlines_" + str(levs[1]) + "m_" + str(ct)[11:19])
    # extract the lines from the temporary figure
    lines = res.lines.get_paths()
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for line in lines:
        old_x = line.vertices.T[0]
        old_y = line.vertices.T[1]
        
        print(old_x, old_y)
        # apply 2d to 3d transformation
        new_z = np.exp(-(old_x ** 2 + old_y ** 2) / 40)
        new_x = 1.2 * old_x
        new_y = 0.8 * old_y
        
        ax.plot(new_x, new_y, new_z, 'k')
        ax.set_zlim([0, 0.2])
        exit # test
        
    plt.savefig(save_path + "3d_streamlines_" + str(levs[1]) + "m_" + str(ct)[11:19])
    plt.close()    