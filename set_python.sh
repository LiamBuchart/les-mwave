#!/bin/bash

# script that you run to load all the modules you need to run a python script
# also opens your virtual environment - assuming you called it ENV 
# and that it is in /bin and you activate it normally [active file present]
# though can specify path to your environment if youd like

# run with [source set_python.sh]

env_path=$PWD
env_path+="/ENV/bin/activate"
echo "$env_path"

module load python
module load scipy-stack
module load mpi4py

module load gdal
module load geos

source "$env_path"