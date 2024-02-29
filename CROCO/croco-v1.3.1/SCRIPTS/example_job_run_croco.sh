#!/bin/bash
#MSUB -r CROCO                # Request name
#MSUB -o CROCO.o              # Standard output file name
#MSUB -e CROCO.e              # Standard error file name
#MSUB -j oe                   # 
#MSUB -n 8                    # Total number of mpi tasks to use
#MSUB -T 7200                 # Walltime   
#MSUB -x
#MSUB -q skylake              # Standard for thin nodes, hybrid for GPU partition, xlarge for fat nodes 
#MSUB -A yourproject          # Project's quota used
#MSUB -m scratch,store,work
#===============================================================================

#===============================================================================

umask 022
set -u

# source your environment file
source myenv_mypath.sh

# launch the run script
./run_croco.bash >& run_croco.out
