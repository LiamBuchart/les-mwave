"""
    
Function to "idealize" the Vernon Sounding data to better match 
historical downslope windstorm upstream soundings as in
https://doi.org/10.3390/atmos13050765  

Requires previously made input_sounding file (see gen_sounding)
then modified to reduce winds from 1000-5000 m AGL, only a small increase in
this layer before encountering th upper jet.
Output: input_sounding, save as gen_sounding, with modified winds

August 23, 2022
lbuchart@eoas.ubc.ca
"""

import numpy as np
import pandas as pd

from context import script_dir

##########

def make_ideal(heights, low, high, u):
    # idealize the sounding 
    ind4 = heights.index(5000)
    
    # make the wind speed constant from 1000 to 5000 m
    for ii in range(low, ind4):
        print(low, ind4)
        ll = u[low-1]
        hh = u[ind4]
        u[ii] = np.mean([ll, hh])
    
    # linearly increase from 5000-6000 m
    for ii in range(ind4, high+1):
        vals = len(np.arange(ind4, high+1))
        end = u[high+1]
        start = u[ind4]
        new = np.linspace(start, end, vals)
        u[ii] = new[ii-ind4]
    
    return u