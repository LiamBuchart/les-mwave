"""
    
    Calculate the hot dry windy index 
    Max wind * Vapour pressure deficit in lowest 50mb
    Save both components in dataframe
    
    November 28, 2023
    lbuchart@eoas.ubc.ca
    
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import metpy.calc as mpcalc
import pandas as pd

from metpy.units import units
from icecream import ic
from context import name_dir, script_dir, json_dir
from file_funcs import (setup_script, get_heights, weighted_mean_std_3d, 
                        dist_from_ridge)

from wrf import (getvar, xy, interplevel, interpline,
                extract_times)

# get experiment names from our json dictionary
with open(str(json_dir) + "config.json") as f:
    config = json.load(f)
exps = config["exps"]["names"]

# initialize data array to store HDW values 
mw_df = pd.DataFrame( )
vpd_df = pd.DataFrame( )

# loop through experiments to complete calculations
for ex in exps:  
    count = 0
    print("Experiment: ", ex)
    
    # get all required paths and files from the experiment
    path, save_path, relevant_files, wrfin = setup_script(ex)
    ic(path)
    
    # extract heights of all layers 
    height = getvar(wrfin[0], "height", units="m", meta=True)
    
    for ii in range(0, len(wrfin)):
        ncfile = wrfin[ii]
        
        ct = extract_times(ncfile, timeidx=0)
        
        # import all velocity directions
        U = getvar(ncfile, "ua", units="m s-1", meta=True) * units('m/s')
        V = getvar(ncfile, "va", units="m s-1", meta=True) * units('m/s')
        W = getvar(ncfile, "wa", units="m s-1", meta=True) * units('m/s')
        
        SPD = np.sqrt( (U**2) + (V**2) + (W**2) )
        
        # import pressure, mixing ratio and temperature
        P = getvar(ncfile, "pres", units="hPa", meta=True)
        T = getvar(ncfile, "temp", units="degC", meta=True)
        RH = getvar(ncfile, "rh", meta=True)
        
        # calculate mixing ratio, vapour pressure and saturation vapour pressure
        Q = mpcalc.mixing_ratio_from_relative_humidity(P * units('hPa'), T * units('degC'), RH) 
        
        VP = mpcalc.vapor_pressure(P * units.hPa, Q * units('g/kg'))  
        SVP = mpcalc.saturation_vapor_pressure(T * units.degC)  
        
        # vapour pressure deficit
        VPD = (SVP - VP) / 100  # convert to hPa
        #VPD.to('hPa')

        
        # now only get values in lowest 50mb
        PM = np.max(P)
        P_use = PM - 50
        
        SPD = np.where(P > P_use, SPD, -100)
        VPD = np.where(P > P_use, VPD, -100)
        
        # save the two values to
        if count == 0:
            mspd = np.array(np.max(SPD))
            Times = np.array(ct)
            mvpd = np.array(np.max(VPD))
        else: 
            mspd = np.append(mspd, np.max(SPD))
            Times = np.append(Times, ct)
            mvpd = np.append(mvpd, np.max(VPD))
            
        count += 1 
    
    # append the speeds to the dataframe 
    mw_df[str(ex)] = mspd
    vpd_df[str(ex)] = mvpd
    
mw_df["Times"] = Times 
vpd_df["Times"] = Times

# save dataframes
mw_df.to_csv("HDW_max_winds.csv", sep=",", header=True) 
vpd_df.to_csv("HDW_vapour_pressure_deficit.csv", sep=",", header=True)  
    
print("Complete")  