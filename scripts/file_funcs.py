"""

Useful and reocurring functions for analyzing WRF output

"""

import os
import numpy as np
import json
import matplotlib.pyplot as plt

from icecream import ic
from netCDF4 import Dataset
from wrf import getvar
    
from context import name_dir, script_dir, json_dir

# # # # # # # # # #

def setup_script(exp, ):
    path = name_dir + exp + "/output/"

    save_path = script_dir + "figures/" + exp
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    # create list of file names (sorted) 
    all_files = sorted(os.listdir(path))

    # get wrfoutfiles
    relevant_files = get_wrfout(all_files)
    print(relevant_files)

    # import all datasets
    wrfin = [Dataset(path+x) for x in relevant_files]
    
    return path, save_path, relevant_files, wrfin

def get_wrfout(all_files):
    file_start = "wrfout"
    # function to just grab the wrfout files from a given directroy
    relevant_files = []
    for file_name in all_files: 
        if file_name.startswith(file_start):
            relevant_files.append(file_name)
    
    return relevant_files

def get_heights(file):
    # returns model grid heights for all levels in a 3d array (does not take into account terrain height)
    PH = getvar(file, "PH")
    PHB = getvar(file, "PHB")
    
    height = (PH + PHB) / 9.81
    
    return height

def dist_from_ridge(ter_cross):
    # use a cross section of terrain to get a 1d vector of distance from ridge high point

    # get grid resolution from our json dictionary
    with open(str(json_dir) + "config.json") as f:
        config = json.load(f)
    dx = config["grid_dimensions"]["dx"]

    # find location ridge
    max_ind = np.where(ter_cross == np.max(ter_cross))

    # create vector of distance
    ridge_dist = np.arange(1, len(ter_cross)+1, 1)
    for ii in range(0, len(ridge_dist)):
        ridge_dist[ii] = (ii - max_ind[0][0]) * dx

    return ridge_dist

def vert_dist_center(cross):
    # find distance from the center of the grid in the y direction
    
    # get grid info from the json dictionary
    with open(str(json_dir) + "config.json") as f:
        config = json.load(f)
    dy = config["grid_dimensions"]["dy"]
    
    # find mid point
    mid_ind = (len(cross) / 2) 
    
    # create vector distance
    dy_dist = np.arange(1, len(cross)+1, 1)
    for ii in range(0, len(dy_dist)):
        dy_dist[ii] = -(ii - mid_ind) * dy
        
    return dy_dist

def fmt(x):
    # This custom formatter removes trailing zeros, e.g. "1.0" becomes "1"
    # used for when wanting nice labels on a contour plot
    # from: https://matplotlib.org/stable/gallery/images_contours_and_fields/contour_label_demo.html
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s}" if plt.rcParams["text.usetex"] else f"{s}"

def extract_val(arr):
    new_var = ""
    for chars in range(0, len(arr)):
        vals = arr[chars]
        for item in vals:
            if item.isdigit():
                new_var = new_var + str(item)

    return new_var

def weighted_mean_std_3d(array, hh, all_stats,
                         array_sorted=False):
    # function which takes a weighted mean based off grid volume
    # get grid info from the json dictionary
    # only required when take a 3d mean of a wrf file
    # calculate std and desired quantiles (np.array [0-100%])
    with open(str(json_dir) + "config.json") as f:
        config = json.load(f)
    
    quantiles = [25, 50, 75]  # change for different values
    dx = config["grid_dimensions"]["dx"]
    dy = config["grid_dimensions"]["dy"]
    
    # get volume m^3 
    vol = hh * dx * dy
    volw = vol / np.sum(vol)  # relative volumetric weigth of each grid box
    
    # weighted mean
    wm = np.nansum(volw * array) / np.nansum(volw) 
    var = np.float( np.average((array-wm)**2, weights=volw) )
    
    if all_stats:
        af = np.ndarray.flatten(array)
        vf = np.ndarray.flatten(volw)
        if not array_sorted: 
            ind = np.argsort(af)
            d = af[ind]
            w = vf[ind]
        
        p = 1.*np.cumsum(w)/np.sum(w) * 100
        wq = np.interp(quantiles, p, d)
    
    # fix function to either just grab wm or all stats (all_stats=True or False)
    if all_stats: 
        print("Returning ", quantiles, " percentiles")
        return wm, var, np.sqrt(var), wq, quantiles
    else:
        return wm
    
def check_folder(path, name):
    # check if a folder exists
    # if does not exist create it
    # rreturns folder name (assuming you want to use it in your script!)
    
    new_dir = os.path.join(path, name)
    isExist = os.path.exists(new_dir)
    folder = path + name
    if not isExist:
        os.makedirs(new_dir)
        
    return folder

def cb_color_palete( key ):
    
    # a color blind friendly color cycle for all plots  
    # https://gist.github.com/thriveth/8560036
    CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#a65628', '#999999', '#984ea3',
                  '#e41a1c', '#dede00', '#f781bf']
    
    # use these colors combined with experiment names for definitve colors for each experiment
    # output as a dictionary       
    count = 0
    col_palete = {}    
    for ii in key:
        col_palete[ii] = CB_color_cycle[count]  
        count += 1     
        
    return col_palete

def mom_flux_names():
    
    components = ["x-mom", "y-mom", "z-mom", "x-press", "flux-divg"]
    
    return components 