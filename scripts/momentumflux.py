"""
    
    Calculate the vertical momentum flux through a volume placed
    to the lee of the slope near ridge top height (single term of momentum flux equation)    
    Normalized 
    
    lbuchart@eoas.ubc.ca
    October 8, 2023

"""

import json
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import metpy.calc as mpcalc
import pandas as pd
import math

from icecream import ic
from metpy.units import units
from context import name_dir, script_dir, json_dir
from file_funcs import (setup_script, mom_flux_names)

from wrf import (getvar, xy, interplevel, interpline, 
                CoordPair, get_cartopy, to_np,
                extract_times, vertcross)

### USER INPUTS ###

start = (0, 100)  # start of terrain cross section
end = (-1, 100)  # end of terrain cross section
top = 1500  # top of volume
bot = 1000   # bottom of volume
width = 20  # width of the control volume [# grid cells]
bl = (0, 0)  # bottom left of volume
br = (0, -1)  # bottom right of volume
tl = (-1, 0)  # top left of volume
tr = (-1, -1)  # top right of volume
vh = 10  # height to get sfc winds

### END USER INPUTS ###

# get grid resolution from our json dictionary
with open(str(json_dir) + "config.json") as f:
    config = json.load(f)
exps = config["exps"]["names"]
ndx = config["grid_dimensions"]["ndx"]
ndy = config["grid_dimensions"]["ndy"]
dx = config["grid_dimensions"]["dx"]
dy = config["grid_dimensions"]["dy"]

# grab the component names of the momentum flux 
components = ["x-mom", "y-mom", "z-mom", "x-press", "flux-divg"]
cc = ["xm", "ym", "zm", "px", "fd"]  # names of the python variables

# initialize dataframe to store flux values
df = pd.DataFrame( )
# dataframe to max velocity in the control volume
df_speed = pd.DataFrame( )
# a dataframe will all components in all experiments
df_full = pd.DataFrame( )

def airdensity(T, P):
    R = 287  # dry air constant [J K-1 kg-1]
    Rair = (P * R) / T
    
    return Rair

# define the mask of extract values from
path, save_path, relevant_files, wrfin = setup_script(exps[0])

ter = getvar(wrfin[0], "ter", meta=True)
ter_cross = interpline(ter, 
                       start_point=CoordPair(start[0], start[1]),
                       end_point=CoordPair(end[0], end[1]))

# find location ridge
max_ind = np.where(ter_cross == np.max(ter_cross))
minx = max_ind[0][0] 
maxx = max_ind[0][0] + width

# arbitariliy choose north south extent (just wide eough to encompass a good portion of the domain
# and avoid any edge effects)
miny = 35
maxy = ndy - 35

xx = [minx, maxx]
yy = [miny, maxy]

# loop through experiments to complete calculations
for ex in exps:  
    count = 0
    print("Experiment: ", ex)
    # initialize dataframe to store all momentum flux components
    df_comp = pd.DataFrame( )
    
    # get all required paths and files from the experiment
    path, save_path, relevant_files, wrfin = setup_script(ex)
    ic(path)
    
    # extract heights of all layers 
    height = getvar(wrfin[0], "height", units="m", meta=True)[:, xx[0]:xx[1], yy[0]:yy[1]]
    
    # loop through each of the output timesteps
    for ii in range(0, len(wrfin)):
        ncfile = wrfin[ii]
        
        ct = extract_times(ncfile, timeidx=0)
        
        # import all velocity directions
        U = getvar(ncfile, "ua", units="m s-1", meta=True)[:, xx[0]:xx[1], yy[0]:yy[1]]
        V = getvar(ncfile, "va", units="m s-1", meta=True)[:, xx[0]:xx[1], yy[0]:yy[1]]
        W = getvar(ncfile, "wa", units="m s-1", meta=True)[:, xx[0]:xx[1], yy[0]:yy[1]]
        
        # get air density
        P = getvar(ncfile, "pres", units="Pa", meta=True)[:, xx[0]:xx[1], yy[0]:yy[1]]  
        T = getvar(ncfile, "tk", meta=True)[:, xx[0]:xx[1], yy[0]:yy[1]]
        
        rho = airdensity(T, P)
        
        # interpolate rho, u, w to the top and bottom of the volume
        u_bot, u_top = interplevel(U, height, [bot, top], meta=False)  
        w_bot, w_top = interplevel(W, height, [bot, top], meta=False)         
        rho_bot, rho_top = interplevel(rho, height, [bot, top], meta=False)
        
        # interpolate rho, u, pres to the end of the volume (x-dir)
        # first define a sequence of heights to interpolate on in the vertical between top and bot of our volume
        h_levs = np.arange(bot, top, 10)
        
        rho_left = vertcross(rho, height, levels=h_levs, 
                             start_point=CoordPair(bl[0], bl[1]), 
                             end_point=CoordPair(tl[0], tl[1]))
        rho_right = vertcross(rho, height, levels=h_levs, 
                              start_point=CoordPair(br[0], br[1]), 
                              end_point=CoordPair(tr[0], tr[1]))
        
        u_left = vertcross(U, height, levels=h_levs, 
                           start_point=CoordPair(bl[0], bl[1]), 
                           end_point=CoordPair(tl[0], tl[1]))
        u_right = vertcross(U, height, levels=h_levs, 
                            start_point=CoordPair(br[0], br[1]), 
                            end_point=CoordPair(tr[0], tr[1]))
        
        P_left = vertcross(P, height, levels=h_levs, 
                           start_point=CoordPair(bl[0], bl[1]), 
                           end_point=CoordPair(tl[0], tl[1]))
        P_right = vertcross(P, height, levels=h_levs, 
                            start_point=CoordPair(br[0], br[1]), 
                            end_point=CoordPair(tr[0], tr[1]))
        
        # same interpolation for the north and south values (rho, u, v)
        rho_north = vertcross(rho, height, levels=h_levs, 
                              start_point=CoordPair(tl[0], tl[1]), 
                              end_point=CoordPair(tr[0], tr[1]))
        rho_south = vertcross(rho, height, levels=h_levs, 
                              start_point=CoordPair(bl[0], bl[1]), 
                              end_point=CoordPair(br[0], br[1]))
        
        u_north = vertcross(U, height, levels=h_levs, 
                            start_point=CoordPair(tl[0], tl[1]), 
                            end_point=CoordPair(tr[0], tr[1]))
        u_south = vertcross(U, height, levels=h_levs, 
                            start_point=CoordPair(bl[0], bl[1]), 
                            end_point=CoordPair(br[0], br[1]))
        
        v_north = vertcross(V, height, levels=h_levs, 
                            start_point=CoordPair(tl[0], tl[1]), 
                            end_point=CoordPair(tr[0], tr[1]))
        v_south = vertcross(V, height, levels=h_levs, 
                            start_point=CoordPair(bl[0], bl[1]), 
                            end_point=CoordPair(br[0], br[1]))
        
        # calculate all values of the momentum budget
        # excluding coriolis and the horz (x,y) flux divergence terms
        
        # x-momentum
        um_left = rho_left * u_left * u_left
        um_right = rho_right * u_right * u_right
        
        uxm = np.mean( (um_left - um_right) / (width * dx) )
        
        # v-momentum   
        vm_south = rho_south * u_south * v_south
        vm_north = rho_north * u_north * v_north
        
        vym = np.mean( (vm_north - vm_south) / ((maxy - miny) * dy) )
        
        # w-momentum
        wm_bot = rho_bot * u_bot * w_bot
        wm_top = rho_top * u_top * w_top
        
        wzm = np.mean( (wm_top - wm_bot) / (top - bot) )
        
        # pressure differential
        dpx = np.mean( (P_left - P_right) / (width * dx) )
        
        # vertical flux divergence (get the covariance at each height level)
        wp_top = w_top - np.mean(w_top)
        wp_bot = w_bot - np.mean(w_bot)
        
        up_top = u_top - np.mean(u_top)
        up_bot = u_bot - np.mean(u_bot)
        
        cov_top = np.mean(wp_top * up_top)
        cov_bot = np.mean(wp_bot * up_top)
        
        flux_div = (cov_top - cov_bot) / (top - bot)
        
        # flip the vertical coordinate (positive pointing towards the surface)
        # apply to to the w-momentum and flux divergence
        wzm = wzm * -1
        flux_div = flux_div * -1
        
        # full momentum
        M = uxm + vym + wzm + dpx + flux_div
        ic(np.max(M), np.min(M), np.mean(M))
        
        # now get the max sfc velocity beneath the control volume
        u10 = np.max( interplevel(U, height, vh, meta=False) )
        v10 = np.max( interplevel(V, height, vh, meta=False) )
        w10 = np.max( interplevel(W, height, vh, meta=False) )
        
        spd10 = np.sqrt(u10**2 + v10**2 + w10**2)
        ic(spd10)
        
        # place momentums in a vector for dataframe
        # components into vectors
        if count == 0:
            Mtm = np.array(M)
            spd = np.array(spd10)
            Times = np.array(ct)
            
            xm = np.array(uxm)
            ym = np.array(vym)
            zm = np.array(wzm)
            px = np.array(dpx)
            fd = np.array(flux_div)
        else: 
            Mtm = np.append(Mtm, M)
            spd = np.append(spd, spd10)
            Times = np.append(Times, ct)
            
            xm = np.append(xm, uxm)
            ym = np.append(ym, vym)
            zm = np.append(zm, wzm)
            px = np.append(px, dpx)
            fd = np.append(fd, flux_div)

        count += 1       
        
    # append each componenet of momentum budget to experiment dataframe
    for ii in range(len(cc)):  
        component = globals()[cc[ii]] 
        df_comp[str(components[ii])] = component
    
        data = {"value": component, "component": components[ii], "exp": ex, "Time": ct}  
        # fill our full dataframe with the above dictionary (easiest to normalize this way)
        df_full = pd.concat([df_full, pd.DataFrame.from_dict(data)])    
        
        df_comp["Times"] = Times
        
    # append the momentums and speeds to the dataframe 
    df[str(ex)] = Mtm
    df_speed[str(ex)] = spd
    
    # save the individual dataframe
    df_comp.to_csv("leeslope_momflux_component_" + str(ex) + ".csv", sep=",", header=True)
    df_full.to_csv("leeslope_all_exps_comps.csv", sep=",", header=True)

df["Times"] = Times
df_speed["Times"] = Times

# save the dateframe
df.to_csv("leeslope_momentum_flux.csv", sep=",", header=True)
df_speed.to_csv("leeslope_maxvel.csv", sep=",", header=True)

print("Complete")