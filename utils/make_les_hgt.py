"""
Function which takes a terrain type input:
***** gaussian, sinusoidal, plateau, angular *****
and creates idealized terrain heights based on namelist input
information for your idealized WRF run,

Used to create topography for downslope wind storm experiments.
Note that all configurations with have a max height of 1000m.

lbuchart@eoas.ubc.ca
"""
import numpy as np
from context import json_dir, name_dir
import json

author="lbuchart"


terrain_type = "agnesi"  # what type of simple terrain (pick from list above in ***** ******
experiment = "neutstab"  # in the exp directory pick which experiment you will be making terrain for
outfile_name = "input_ht" 

# load grid dimensions which are saved in the context file
with open(str(json_dir) + "config.json") as f:
    config = json.load(f)

we_n = config["grid_dimensions"]["ndx"]
sn_n = config["grid_dimensions"]["ndy"]
mht = config["topo"]["mht"]  # topography height
dx = config["grid_dimensions"]["dx"]  # grid dimensions

# create a grid of zeros based on namelist size
hgt = np.empty((we_n, sn_n), np.float32)
grid_size = np.shape(hgt)
print(hgt.ndim, grid_size[0], grid_size[1])

# loop through the grid and add topography heights based on the chosen methods
# note that topography is added in the east-west direction
for ii in range(grid_size[1]):

    if terrain_type.__eq__("gaussian"):
        
        def gaussian(mht, x, mu, sig):
            # x = range over which gaussian hill is made
            # mu = offset of the guassian on domain
            # sig = steepness factors (>1 steeper)
            return  mht*np.exp( -np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
            
        mu = grid_size[0] / 4  # change gaussian width
        sig = 4
        x = np.arange(grid_size[0])
        
        hgt[:, ii] = gaussian(mht, x, mu, sig)  
    
    elif terrain_type.__eq__("agnesi"):
        
        def agnesi(mht, x, hw):
            # mht is the mountain/ridge maximum height
            # x = range over which the hill is made
            # hw = the half width of the mountain/ridge
            
            return mht * ( (hw**2) / ((hw**2) + (x**2)) )
        
        hw = 1000
        pos_peak = 6000   # position in m from the grid edge (x) that ridge peak is located
        x = (np.arange(grid_size[0]) * dx) - pos_peak
        
        hgt[:, ii] = agnesi(mht, x, hw)

# now add the required info and array to the outfile_name text file which we will create
dims = str(we_n) + " " + str(sn_n)
print(dims)

np.savetxt(outfile_name, hgt, header=dims, comments="", fmt="%i")
print("Saved the height array")