#!/bin/bash
#MSUB -r Comp_wrf       # request name
#MSUB -o Comp_wrf.o
#MSUB -e Comp_wrf.e
#MSUB -j oe          #  VERIFIER que c'est OK sur Curie XXX
#MSUB -n 8        # Total number of mpi task to use
#MSUB -T 7200        
#MSUB -x
#MSUB -q skylake 
#MSUB -A yourproject
#MSUB -m scratch,store,work
#===============================================================================

umask 022
cd ${BRIDGE_MSUB_PWD}
source ../myenv_mypath.sh
./make_WRF_compil
