#!/bin/bash
#SBATCH --job-name=Comp_wrf
#SBATCH --partition=slims
#SBATCH -n 8
#SBATCH --ntasks-per-node=4
#SBATCH --time=12:00:00              # temps d execution maximum demande (HH:MM:SS)
#SBATCH --output=comp_wrf.out
#SBATCH --error=comp_wrf.err

cd ${SLURM_SUBMIT_DIR}
source ../myenv_mypath.sh
./make_WRF_compil
