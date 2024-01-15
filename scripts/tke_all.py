"""

    Calculate the turbulent kinetic energy above the pbl
    for each experiment, plot time series (1/2(u'2+v'2+w'2))
    calculate potential energy (1/2(g/N)2*(bar(T'/T0)2))
    
    Clarification: Ek and Ep are calculated according to VanZandt (1985)
    according to their calculations on internal gravity waves
    
    lbuchart@eoas.ubc.ca
    May 4, 2023

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

from wrf import (getvar, xy, interp2dxy, interpline, 
                CoordPair, get_cartopy, to_np,
                extract_times)

#####

# get grid resolution from our json dictionary
with open(str(json_dir) + "config.json") as f:
    config = json.load(f)
exps = config["exps"]["names"]

# initialize data array to store energy values 
df = pd.DataFrame( )

# loop through experiments to complete calculations
for ex in exps:  
    count = 0
    print("Experiment: ", ex)
    # get all required paths and files from the experiment
    path, save_path, relevant_files, wrfin = setup_script(ex)
    
    # extract heights of all layers 
    heights = getvar(wrfin[0], "height", units="m", meta=False)
    hbf = getvar(wrfin[0], "height", units="m", meta=True)
    geopotential = getvar(wrfin[0], "geopotential", meta=True)
    terh = getvar(wrfin[0], "ter", units="m")
    
    # subtract terrain height from model grid heights [want AGL]
    for ll in range(0, np.shape(heights)[0]):
        heights[ll, :, :] = heights[ll, :, :] - terh
    
    # loop through each of the output timesteps
    for ii in range(0, len(wrfin)):
        
        ncfile = wrfin[ii]
        
        ct = extract_times(ncfile, timeidx=0)
        ic(ct)
        
        ## get pbl height - calculate as the maximum vertical gradient of the mean virtual potential temperature
        # import vapour mixing ratio
        mr = np.squeeze( ncfile['QVAPOR'][:, :, :] ) # grab values like this to convery to numpy array
        q2 = ncfile['Q2'][:, :, :]
        t2 = ncfile['T2'][:, :, :]

        # import the potential temperature
        tt = to_np( getvar(ncfile, "theta", units="K") )
        
        # calculate the virtual potential temperature as 
        # theta*(1+0.61[QVAPOUR]+[QVAPOUR_L]) -> unsaturated QVAPOUR_L=0
        tV = np.squeeze( tt*(1+0.61*(mr)) )
        tV2 = t2*(1+0.61*(q2))

        # take spatial mean, leaving vertical profile
        vtV = tV.mean(axis=(1, 2))    
        vh = heights.mean(axis=(1, 2))
        tV2m = np.mean(tV2)  
        
        # append 2m values to bottom of the virtual potential temperature vector
        vtV = np.append(tV2m, vtV)

        vh[0] = 0  # taking spatial mean alongside the height subtraction results in slightly negative first level
        
        # create a gradient vector of the virtual potential temperature
        gradtV = np.zeros(len(vtV)-1)
        for ii in range(0, len(vh)-1):
            dy = vh[ii+1] - vh[ii]
            dtV = vtV[ii+1] - vtV[ii]
            
            gradtV[ii] = dtV / dy
            
        #plt.figure(figsize=(12,12))
        #plt.plot(gradtV, vh, "k.")
        #plt.title("Vertical Gradient of Mean Virtual Potential Temperature")
        #plt.savefig(save_path + "grad_vpt_" + str(ct)[11:19] + ".png")
        
        # define max of the gradient as the top of the boundary layer
        # grab the eta level where this corresponds (to be used later)
        gradtV = gradtV[1:] # remove surface layer gradient
        gradtV = gradtV[:-7]  # remove all values near top of the troposphere (only for getting PHL height)
        ind = np.argmax(gradtV)
        
        # calculate kinetic anad potential internal wave energy 
        # above the boundary layer 

        # Ek [J/kg]
        U = getvar(ncfile, "ua", units="m s-1",
            meta=False)[ind+1:, :, :]
        Umn = weighted_mean_std_3d(U, heights[ind+1:, :, :], all_stats=False)
        V = getvar(ncfile, "va", units="m s-1",
            meta=False)[ind+1:, :, :]
        Vmn = weighted_mean_std_3d(V, heights[ind+1:, :, :], all_stats=False)
        W = getvar(ncfile, "wa", units="m s-1",
            meta=False)[ind+1:, :, :]
        Wmn = weighted_mean_std_3d(W, heights[ind+1:, :, :], all_stats=False)
        
        # anomaly values
        Up = U - Umn
        Vp = V - Vmn
        Wp = W - Wmn
        
        Ek = 0.5 * ( (Up)**2 + (Vp)**2 + (Wp)**2 )
        
        # Ep [J/kg]
        Th = getvar(ncfile, "theta", units="K", meta=True)
        temp = getvar(ncfile, "tk", meta=False)[ind+1:, :, :]  # temperature above boundary layer        
        T0 = weighted_mean_std_3d(temp, heights[ind+1:, :, :], all_stats=False)  # spatial mean temp above boundary layer  
        N = mpcalc.brunt_vaisala_frequency(hbf, Th)
        N = np.array(N)[ind+1:, :, :]  # conver to numpy array for easier use and take only above boundary layer
        
        Tp0 = ( (temp - T0) / T0 ) ** 2
        #Tp0 = weighted_mean_std_3d(Tp0, heights[ind+1:, :, :], all_stats=False)
        N0 =  weighted_mean_std_3d(N, heights[ind+1:, :, :], all_stats=False)
        ic(N0)
        
        ic(ind, 9.81, np.nanmean(N), N0, T0)
        Ep = 0.5 * ( (9.81/N0)**2 * np.ma.asarray(Tp0) )  # for consistent array typing used the .ma.asarray() method
        
        # calculate energ states
        
        Ed = Ep + Ek   # get the energy density
        
        # stats on energy density
        # calculated weighted variance and standard deviation (percentiles possible)
        Edmn, Edvar, Edstd, Edquant, quantiles = weighted_mean_std_3d(Ed, heights[ind+1:, :, :], all_stats=True)
        Epmn = weighted_mean_std_3d(Ep, heights[ind+1:, :, :], all_stats=False)
        Ekmn = weighted_mean_std_3d(Ek, heights[ind+1:, :, :], all_stats=False)
        
        #ic(Edmn, Ekmn, Epmn)
        
        # prepare vectors to save into dataframe for plotting
        if count == 0:
            Edv = np.array(Edmn)
            Ekv = np.array(Ekmn)
            Epv = np.array(Epmn)
            # sum of the energy density
            Edsumv = np.array(np.sum(Ed))
            Ed25 = np.array(Edquant[0])  # 25 percentile
            Ed50 = np.array(Edquant[1])  # 50 percentile
            Ed75 = np.array(Edquant[2])  # 75 percentile
            Edstdv = Edstd
            
            Times = np.array(ct)
        else: 
            Edv = np.append(Edv, Edmn)
            Ekv = np.append(Ekv, Ekmn)
            Epv = np.append(Epv, Epmn)
            
            Edsumv = np.append(Edsumv, np.sum(Ed))
            Ed25 = np.append(Ed25, Edquant[0])
            Ed50 = np.append(Ed50, Edquant[1])
            Ed75 = np.append(Ed75, Edquant[2])
            Edstdv = np.append(Edstdv, Edstd) 
            
            Times = np.append(Times, ct)
        
        count += 1
                
        ic(count)   
    
    # append arrays to dataframe df
    vals = { "Edmn": Edv, "Edsum": Edsumv, "Edstd": Edstdv,
            "Ed25": Ed25, "Ed50": Ed50, "Ed75": Ed75,
            "Ekmn": Ekv, "Epmn": Epv, "Times": Times}
    
    df = pd.DataFrame.from_dict(vals)
    
    print(script_dir)
    df.to_csv(script_dir + "/Energy_" + ex + ".csv", sep=",", header=True)
    
print("Complete")