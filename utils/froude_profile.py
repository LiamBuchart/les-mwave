"""
    
    Return the Froude number (F = Uh/N) for the wind profile 
    that is used in all simulations.
    Used to ensure F~=1 for idealized mountain wave simulations
    
    March 14, 2023
    lbuchart@eoas.ubc.ca    
    
"""

import pandas as pd
import metpy.calc as mpcalc
import numpy as np
import matplotlib.pyplot as plt
import json 
import xarray as xr

from metpy.units import units
from context import json_dir, script_dir, name_dir

##### USER INPUTS #####

exp = "low_step"
infile = "/input_sounding"

##### END USER INPUT #####

# load sounding header names and mountain height
with open(str(json_dir) + "config.json") as f:
    config = json.load(f)
        
headers = config["sounding"]["headers_exps"]
ht = config["topo"]["mht"]

# load the sounding
path = name_dir + exp 
df = pd.read_csv(path + infile, sep=' ')
df = df.dropna()  # get rid of first row which defined surface press, temp, and moisture

su = df[headers[3]].to_xarray() * units("m/s")  # grab the third column for U, see config.json for where U and V are placed in sounding

sh = df[headers[0]].to_xarray() * units("m")  # grab the heights, see config.json for where U and V are placed in sounding
st = df[headers[1]].to_xarray() * units("K")  # grab theta values

# calculate the brunt vaisale frequency
N = mpcalc.brunt_vaisala_frequency(sh, st)

print("Fr = ", (np.mean(su[10:20])/ht)/np.mean(N[10:20]))

print("Complete")