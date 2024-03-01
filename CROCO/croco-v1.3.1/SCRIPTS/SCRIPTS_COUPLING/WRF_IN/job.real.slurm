#!/bin/bash
#SBATCH --job-name=Comp_real
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=40
#SBATCH --time=12:00:00              # temps d execution maximum demande (HH:MM:SS)
#SBATCH --output=comp_wps.out
#SBATCH --error=comp_wps.err
#SBATCH -A yourproject

cd ${SLURM_SUBMIT_DIR}

#===============================================================================

umask 022
set -u

./run_real.bash configure.namelist.real 8

