#!/bin/bash
#SBATCH --job-name=CROCO
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=40
#SBATCH --time=12:00:00              # Walltime (HH:MM:SS)
#SBATCH --output=CROCO.out           # Standard output file name
#SBATCH --error=CROCO.err            # Standard error file name
#SBATCH -A yourproject               # Project's quota used

cd ${SLURM_SUBMIT_DIR}

#===============================================================================

umask 022
set -u

# source your environment file
source myenv_mypath.sh

# launch the run script
./run_croco.bash >& run_croco.out

