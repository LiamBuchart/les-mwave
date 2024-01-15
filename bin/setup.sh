#!/bin/bash
#SBATCH -t 00:05:00
#SBATCH --mem-per-cpu=3000M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lbuchart@eoas.ubc.ca
#SBATCH --account=def-rstull

ml StdEnv/2020  intel/2020.1.217  openmpi/4.0.3
module load wrf/4.2.1

cd ../exps/high_step

mpirun -np 1 ./ideal.exe
