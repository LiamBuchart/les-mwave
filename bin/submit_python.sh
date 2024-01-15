#!/bin/bash
# ---------------------------------------------------------------------
# SLURM script for a multi-step job on a Compute Canada cluster. 
# ---------------------------------------------------------------------
#SBATCH --account=rrg-rstull
#SBATCH --time=01:00:00
#SBATCH --mem=10G  ## memory to ask for
#SBATCH --job-name="Wave Plot"
#SBATCH --mail-user=lbuchart@eoas.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --output=sim-%j.out
#SBATCH --error=sim-%j.err
# ---------------------------------------------------------------------
echo "Current working directory: $(pwd) "
echo "Starting run at:" "$(date)" 
# ---------------------------------------------------------------------
# Run

module load python scipy-stack mpi4py
source ~/projects/def-rstull/lbuchart/mwave/ENV/bin/activate

python ~/projects/def-rstull/lbuchart/mwave/scripts/gen_sounding.py
# ---------------------------------------------------------------------
echo "Job finished with exit code $? at:" "$(date)" 

