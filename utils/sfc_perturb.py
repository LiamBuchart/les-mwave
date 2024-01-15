"""

Simple function to perturb the surface temperature in LES mode and create a turbulent boundary layer

Function comes from Chris Rodell

August 2, 2022
lbuchart@eoas.ubc.ca

"""

import json
import numpy as np
from context import json_dir, utils_dir

# load grid dimensions which are saved in the context file
with open(str(json_dir) + "config.json") as f:
    config = json.load(f)

ndx = config["grid_dimensions"]["ndx"]
ndy = config["grid_dimensions"]["ndy"]
sfc_T = config["sounding"]["T_sfc"]

# create the perturbed surface temperature to start convection
surface_T = (
    (np.random.rand(ndx, ndy) - 0.25) * 1.5) + sfc_T

# create the required headers for sfire
dim_header = ",".join(
    map(str, np.shape(surface_T)))

# save the output file in a format that sfire can deal with
np.savetxt(
    str(utils_dir) + "input_tsk", surface_T, header=dim_header, comments="", fmt="%1.2f"
)

print("save the surface temperature perturbation")