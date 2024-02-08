"""
    Do some rudementary offline particle tracking based on the mean 2d cross section of u and w velocities.
    release a number of particles at different heights near the western edge of the domain. 
    heights every 500m

    lbuchart
    January 15, 2024
"""

# imports
import numpy as np
import json
import os
import pandas as pd

from scipy.spatial import KDTree, distance
from wrf import getvar, extract_times
from file_funcs import setup_script
from context import json_dir, name_dir, script_dir
from netCDF4 import Dataset
from icecream import ic

## Particle Tracking Class ##
class ParticleTracker:
    """
    First-Order offline particle tracking in 2D. Using WRF output.
    Processed using WRF-Python to create a experiment mean cross section
    
    Advection (callable): Function to advect the particles in space
    Nearest (callable): find nearest velocity point to the tracker
    setup_scripts (callable): setup files and path for WRF run (Compute Canada configured)
    WRF_proc (callable): function to process wrf output for use
    get_vel (callable): function to grab the velocities from WRF output
    x0z0 (tuple): initial position to release particle from
    u (float): horizontal velocity 
    w (float): vertical velocity
    time_step (float): time step we will march over 
    exp (string): name of experiment to process
    """
    
    def __init__(self, Advection, Nearest, setup_script, WRF_proc, get_vel, x0z0, u, w, time_step, exp):
        
        self.Advection = Advection  # first order linear advection
        self.Nearest = Nearest  # find nearest velocity points to calculate with
        self.setup_script = setup_script  # function to setup wrf processing
        self.WRF_proc = WRF_proc  # function to process and load WRF winds
        self.get_vel = get_vel  # function to grab the velocities from wrf output
        
        self.x0z0 = x0z0  # position of particle on grid
        self.u = u  # horz velocity
        self.w = w  # vertical velocity
        
        self.time_step = time_step  # advection timestep
        self.exp = exp  # name of the experiment to process
          
    def Calculate(self):
        # calculate the initial state
        #ic(self.x0z0)
        index = Nearest(self.x0z0, nodes)
        self.u, self.w = get_vel(index, U, W) 
        
    def Update(self):
        self.x0z0 = Advection(self.time_step, self.x0z0, self.u, self.w)   
        
def Advection(time_step, pos0, u, w):
        
    # linear advection in each direction
        
    x1 = pos0[1] + (u * time_step)
    z1 = pos0[0] + (w + time_step)
               
    return [z1, x1]
    
def Nearest(pos, nodes):
        
    # use a linear interpolation to find u  and w
    # velocity from wrf interpolated output
    nn = nodes[:, 0:2]  # positions
    iis = nodes[:, -2:]  # indexs
    
    # get the distance from and index of the nearest values in nodes from our position
    tree = KDTree(nn)  # create kd_tree to get nearest neightbours from
    dist, idx = tree.query(pos, k=1)  # get the closest position
    
    ids = iis[int(idx)]
    
    index = np.array(ids)
    #ic(dist, idx, index)
        
    return index
            
def setup_script(exp):
    path = name_dir + exp + "/output/"

    save_path = script_dir + "figures/" + exp
    
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    # create list of file names (sorted) 
    all_files = sorted(os.listdir(path))

    # get wrfoutfiles
    file_start = "wrfout"
    # function to just grab the wrfout files from a given directroy
    relevant_files = []
    for file_name in all_files: 
        if file_name.startswith(file_start):
            relevant_files.append(file_name)

    # import all datasets
    wrfin = [Dataset(path+x) for x in relevant_files]
    
    return path, wrfin
    
def WRF_proc(path, wrfin):
        
    # extract heights of all layers 
    heights = getvar(wrfin[0], "height_agl", units="m", meta=False)
    ys = heights[:, 0, 0] 

    # get grid resolution from our json dictionary
    with open(str(json_dir) + "config.json") as f:
        config = json.load(f)
    dx = config["grid_dimensions"]["dx"]
    ndx = config["grid_dimensions"]["ndx"]

    # create vector of distance
    horz_dist = np.arange(1, ndx+1, 1)
    for ii in range(0, len(horz_dist)):
        horz_dist[ii] = ii  * dx
    
    # make a (n, 2) array of every possible location on the grid by height and distance from edge
    nodes = np.array([])
    count = 0    
    for ii in range(len(ys)):
        for jj in range(len(horz_dist)):
            if count == 0:
                nodes = np.hstack(( nodes, np.array([ys[ii], horz_dist[jj], ii, jj]) ))
            else:
                nodes = np.vstack(( nodes, np.array([ys[ii], horz_dist[jj], ii, jj]) ))
        
            count += 1
    
    # loop through files list and make our plots
    all_U = []
    all_W = []
    
    for ii in range(20, len(wrfin)):
        # import the file in a readable netcdf format
        ncfile = wrfin[ii]
    
        # get the time in datetime format
        ct = extract_times(ncfile, timeidx=0)
    
        # winds
        W = getvar(ncfile, "wa", units="m s-1",
                   meta=True)
        U = getvar(ncfile, "ua", units="m s-1",
                   meta=True)

        # take mean over the north south distance
        W = np.mean(W, axis=1)
        U = np.mean(U, axis=1)

        # concatenate over all time to get a mean picture
        WW = W.to_numpy()
        UU = U.to_numpy()
        if ii == 5:
            all_W = WW
            all_U = UU
        elif ii > 5: 
            all_W = np.dstack((all_W, WW))
            all_U = np.dstack((all_U, UU))
            
    final_U = np.mean(all_U, axis=2)
    final_W = np.mean(all_W, axis=2)
    
    return final_U, final_W, nodes
    
def get_vel(index, U, W):
    # function to grab the actual velocity from the wrf output
    # based on the index values returned by the Nearest Function
        
    ii = int(index[0])  # the row
    jj = int(index[1])  # the column
    
    u_use = U[ii, jj]
    w_use = W[ii, jj]
        
    return u_use, w_use
    
def check_edges(pos, edges, is_in):
    # function to check that the position of the particles is inside
    # the edges of the grid
    maxz = edges[1]
    minz = edges[0]
    maxx = edges[3]
    minx = edges[2]
    
    x_pos = pos[1]
    z_pos = pos[0]
    
    if x_pos <= minx or x_pos >= maxx or z_pos <= minz or z_pos >= maxz:
        is_in = False
        
    return is_in
    
##
## simulation parameters ##
# time
num_steps = 500
time_step = 10  # time step [ss]
wrf_dt = 300  # the "time step" of wrf output 5mins [s]

# json file with hard coded variables
with open(str(json_dir) + "config.json") as f:
    config = json.load(f)
# get height of release and x position from json file
xpos = int(config["trajectories"]["xpos"])
heights = config["trajectories"]["zpos"]
ic(heights)

# "initial" velocities
u = 0
w = 0

# get experiment names from our json dictionary
exps = config["exps"]["names"]

edges = [0, 9500, 0, 44000]  # vertical and horizontal grid edges

# empty data frames to store particle trajectories
xdf = pd.DataFrame()
zdf = pd.DataFrame()

## END USER INPUTS ##

##
## initialize the particle tracking ##
# manualy make array with the position np.array(height[z], xpos)
pt = ParticleTracker(Advection, Nearest, setup_script, WRF_proc, get_vel, np.array([heights[0], xpos]), u, w, time_step, exps[0])

for exp in exps:
    print("Experiment to Process: ", exp)
    path, wrfin = setup_script(exp)
    U, W, nodes = WRF_proc(path, wrfin)
    for hh in heights:
        # empty nan vector for particle positions
        X = np.empty(num_steps+1)
        X[:] = np.nan
        Z = np.empty(num_steps+1)
        Z[:] = np.nan
        
        X[0] = xpos  # set the initial position in x-dir
        Z[0] = hh  # set the iniital position in z-dir
        
        print("Initial Particle Position: ", np.array([hh, xpos]))
        pt = ParticleTracker(Advection, Nearest, setup_script, WRF_proc, get_vel, np.array([hh, xpos]), u, w, time_step, exp)
        
        # boolean to show that the particle is the grid space
        is_in = True
    
        for ii in range(1, num_steps+1):
            pt.Calculate()
            pt.Update()
            pos = pt.x0z0
            is_in = check_edges(pos, edges, is_in)
                
            X[ii] = pos[1]
            Z[ii] = pos[0]
            
            if is_in != True:   
                break  # jump out of the loop if we find th point outside the grid
        
        ic(len(X), len(Z))
        
        ## Dictionary here to store trajectory for each height
        
        #xdf.at[hh, exp] = np.array(X)
        #zdf.at[hh, exp] = np.array(Z)
        
        xdf[exp + "_" + str(hh)] = X
        zdf[exp + "_" + str(hh)] = Z
        
        ic(zdf.head())
            
zdf.to_csv("z_particle_trajectory_positions.csv")
xdf.to_csv("x_particle_trajectory_positions.csv")                
            