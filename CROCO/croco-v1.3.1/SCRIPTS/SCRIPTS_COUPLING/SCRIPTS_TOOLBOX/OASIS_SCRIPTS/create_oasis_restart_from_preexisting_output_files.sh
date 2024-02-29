#!/bin/bash -e

## ----------------------------------------------------------------------------- #
## - Create restart file for oasis                                             - #
## - from pre-existing model file                                              - #
##                                                                             - #
## Mandatory inputs:                                                           - #
##  - the model file name (with full path)                                     - #
##  - the oasis restart file name (with full path)                             - #
##  - the model: wrf, croco, or ww3 cases are accepted                         - #
## Optional input:                                                             - #
##  - the grid levels: for croco: 0 for parent, 1 for child, etc               - #
##                     for wrf: examples: WRF_d01_EXT_d01 or WRF_d02_EXT_d01   - #
##                       domain 1 of WRF coupled with domain 1 of other model  - #
##                       domain 2 of WRF coupled with domain 1 of other model  - #
##                                                                             - #
## ----------------------------------------------------------------------------- #
#
# Further Information:   
# http://www.croco-ocean.org
#  
# This file is part of CROCOTOOLS
#
# CROCOTOOLS is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 2 of the License,
# or (at your option) any later version.
#
# CROCOTOOLS is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA
#
# Copyright (c) 2018 S. Jullien
# swen.jullien@ifremer.fr
## ----------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------ #

filein=$1
fileout=$2
model=$3
gridlevels=$4

# ------------------------------------------------------------------------------ #
echo '*******************************************************************'
echo 'START script create_oasis_restart_from_preexisting_output_files.sh'
echo '*******************************************************************'
echo ' '

# First check if inputs are ok
if [[ -z $filein ]] || [[ -z $fileout ]] || [[ -z $model ]] ; then
    echo 'ERROR: inputs are not correctly specified.'
    echo '       this script needs at least 3 inputs:'
    echo '       - the model file name (with full path)'
    echo '       - the oasis restart file name (with full path)'
    echo '       - the model: wrf croco or ww3 cases are accepted'
    echo ' Exit...'
    echo ' '
    exit 1
fi 
# ------------------------------------------------------------------------------ #

mydir=$(dirname "$fileout")
filetmp=$mydir/rst_tmp.nc

if [ $model == wrf ] ; then

    if [ -z $gridlevels ] ; then
      echo 'Default grid levels are assumed: WRF_d01_EXT_d01...'
      gridlevels='WRF_d01_EXT_d01'
    fi

    varlist=(${gridlevels}_TAUE \
            ${gridlevels}_TAUN \
            ${gridlevels}_TAUMOD \
            ${gridlevels}_WND_E_01 \
            ${gridlevels}_WND_N_01 \
            ${gridlevels}_SURF_NET_SOLAR \
            ${gridlevels}_SURF_NET_NON-SOLAR \
            ${gridlevels}_EVAP-PRECIP\
            ${gridlevels}_PSFC)

    dimtime=Time
    timerange=2

elif  [ $model == croco ] ; then

    if [ -z $gridlevels ] ; then
      echo 'Default grid level is assumed: 0 for parent...'
      gridlevels=''
    fi

    varlist=(CROCO_EOCE${gridlevels} \
            CROCO_NOCE${gridlevels} \
            CROCO_SST${gridlevels} \
            CROCO_SSH${gridlevels})

    dimtime=time
    timerange=1

elif  [ $model == ww3 ] ; then

#    varlist=(WW3_T0M1 \
#            WW3___HS \
#            WW3_CDIR \
#            WW3_SDIR \
#            WW3_TWOX \
#            WW3_TWOY \
#            WW3_TAWX \
#            WW3_TAWY \
#            WW3__CHA)

    varlist=(WW3_T0M1 \
            WW3___HS \
            WW3__DIR \
            WW3_TWOX \
            WW3_TWOY \
            WW3_TAWX \
            WW3_TAWY \
            WW3_ACHA)

    dimtime=time
    timerange=1

else
    echo 'ERROR: '$model' case is not implemented yet. Exit...'
    echo ' '
    exit 1
fi # model 

echo ' '
echo 'Varlist to proceed is '$varlist
echo '========================================================================='
lengthvar=${#varlist[@]}
for k in `seq 0 $(( ${lengthvar} - 1))` ; do

    var=${varlist[$k]}
    echo ' '
    echo '======================'
    echo 'Process '$var'...'
    echo '======================'

    # Extract or compute var
    echo '---> Extract or compute '$var
    ${SCRIPTDIR}/OASIS_SCRIPTS/from_${model}.sh $filein $filetmp $var $timerange $gridlevels

    if [ $model == wrf ] ; then
        # Put them on the stag grid
        echo '---> Put them on the stag grid' 
        ${SCRIPTDIR}/OASIS_SCRIPTS/to_wrf_stag_grid.sh $filetmp $filetmp
    fi

    # Remove time dimension
    echo '---> Remove time dimension...'
    ncwa -O -a $dimtime $filetmp $filetmp

    ncks -A -v $var $filetmp $fileout
    rm $filetmp

done # LOOP on varlist

# remove global attributes
ncatted -h -O -a ,global,d,, $fileout $fileout 

echo ' '    
