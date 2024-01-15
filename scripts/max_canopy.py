"""
    
    Plot the maximum canopy height instanteous winds.
    Ell experiments time series. 
    
    lbuchart@eoas.ubc.ca
    September 17, 2023
    
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import metpy.calc as mpcalc
import pandas as pd

from icecream import ic
from metpy.units import units
from context import name_dir, script_dir, json_dir
from file_funcs import (setup_script, get_heights, weighted_mean_std_3d, 
                        dist_from_ridge, fmt)

from wrf import (getvar, xy, interplevel, interpline, 
                CoordPair, get_cartopy, to_np,
                extract_times)

### USER INPUTS ###

lev = 10  # height to interpolate winds to

### END USER INPUTS ###

# get grid resolution from our json dictionary
with open(str(json_dir) + "config.json") as f:
    config = json.load(f)
exps = config["exps"]["names"]

# initialize data array to wind values 
df = pd.DataFrame( )
wm_df = pd.DataFrame( )

# loop through experiments to complete calculations
for ex in exps:  
    count = 0
    print("Experiment: ", ex)
    # get all required paths and files from the experiment
    path, save_path, relevant_files, wrfin = setup_script(ex)
    ic(path)
    
    # extract heights of all layers 
    height = getvar(wrfin[0], "height", units="m", meta=False)
    
    # loop through each of the output timesteps
    for ii in range(0, len(wrfin)):
        ncfile = wrfin[ii]
        
        ct = extract_times(ncfile, timeidx=0)
        #print(ct)
        
        # import all velocity directions
        U = getvar(ncfile, "ua", units="m s-1", meta=False)
        V = getvar(ncfile, "va", units="m s-1", meta=False)
        W = getvar(ncfile, "wa", units="m s-1", meta=False)
        
        # interpolate winds to desired lev height [m]
        Ulev = interplevel(U, height, lev, meta=True) 
        Vlev = interplevel(V, height, lev, meta=True) 
        Wlev = interplevel(W, height, lev, meta=True) 
        
        spd = np.array( np.max( np.sqrt( (np.square(Ulev)) + (np.square(Vlev)) + (np.square(Wlev))) ) )
        wmax = np.min(W)  # simply get the minimum vertical velocity (i.e. max downward)
        
        # place speeds in a vector for dataframe
        if count == 0:
            vel = np.array(spd)
            Times = np.array(ct)
            wm = np.array(wmax)
        else: 
            vel = np.append(vel, spd)
            Times = np.append(Times, ct)
            wm = np.append(wm, wmax)
            
        count += 1
        
    # append the speeds to the dataframe 
    ic(vel)
    df[str(ex)] = vel
    wm_df[str(ex)] = wm

df["Times"] = Times 
wm_df["Times"] = Times   
           
# same the dataframe 
df.to_csv("Canopy_Max_Velocities.csv", sep=",", header=True) 
wm_df.to_csv("Max_Vert_Velocity.csv", sep=",", header=True)  
    
print("Complete")       