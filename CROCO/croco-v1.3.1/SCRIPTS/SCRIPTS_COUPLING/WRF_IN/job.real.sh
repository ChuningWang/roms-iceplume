#!/bin/bash
#MSUB -r Comp_real       # request name
#MSUB -o Comp_real.o
#MSUB -e Comp_real.e
#MSUB -j oe          #  VERIFIER que c'est OK sur Curie XXX
#MSUB -n 8        # Total number of mpi task to use
#MSUB -T 7200        
#MSUB -x
#MSUB -q skylake 
#MSUB -A yourproject
#MSUB -m scratch,store,work
#===============================================================================

#===============================================================================

umask 022
set -u

./run_real.bash configure.namelist.real 8
