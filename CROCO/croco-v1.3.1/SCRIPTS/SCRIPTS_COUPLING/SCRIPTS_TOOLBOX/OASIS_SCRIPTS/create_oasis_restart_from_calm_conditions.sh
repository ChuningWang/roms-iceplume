#!/bin/bash -e

## ----------------------------------------------------------------------------- #
## - Create restart file for oasis                                             - #
## - with all variables set to 0                                               - #
##                                                                             - #
## Mandatory inputs:                                                           - #
##  - a file from this model containing the mask (with full path)              - #
##  - the oasis restart file name (with full path)                             - #
##  - the model: wrf, croco, or ww3 cases are accepted                         - #
##  - the list of variables that have to be generated in this restart file     - #
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
varlist=$4

# ------------------------------------------------------------------------------ #
echo '**********************************************************'
echo 'START script create_oasis_restart_from_calm_conditions.sh'
echo '**********************************************************'
echo ' '

# First check if inputs are ok
if [[ -z $filein ]] || [[ -z $fileout ]] || [[ -z $model ]] || [[ -z $varlist ]] ; then
    echo 'ERROR: inputs are not correctly specified.'
    echo '       this script needs 4 inputs:'
    echo '       - a file from this model containing the mask (with full path)'
    echo '       - the oasis restart file name (with full path)'
    echo '       - the model: wrf croco or ww3 cases are accepted'
    echo '       - the list of variables that have to be generated in this restart file'
    echo ' Exit...'
    echo ' '
    exit 1
fi 
# ------------------------------------------------------------------------------ #

mydir=$(dirname "$fileout")
filetmp=$mydir/rst_tmp.nc

echo 'Initialize restart file '$fileout' to 0 for variables: '$varlist
echo '========================================================================='

if [ $model == wrf ] ; then

    # Extract mask
    echo '---> Extract LANDMASK variable...'
    ncks -O -v LANDMASK -d Time,0 $filein ${filetmp}
    
    # Put it on the stag grid
    echo '---> Put it on the stag grid' 
    ${SCRIPTDIR}/OASIS_SCRIPTS/to_wrf_stag_grid.sh ${filetmp} ${filetmp}
 
    # remove time dimension
    echo '---> Remove time dimension...'
    ncwa -O -a Time ${filetmp} ${filetmp}
    
    # set the variable to 0 and rename it var0
    echo '---> Set the variable to 0 and rename it var0...'
    ncap2 -O -v -s "var0=double(LANDMASK*0)" ${filetmp} ${filetmp}

elif  [ $model == croco ] ; then

    # Extract dimensions
    xi_rho=`ncdump -h $filein | grep "xi_rho = " | cut -d '=' -f2 | cut -d ';' -f1`
    xi_rho=${xi_rho// /}
    eta_rho=`ncdump -h $filein | grep "eta_rho = " | cut -d '=' -f2 | cut -d ';' -f1`
    eta_rho=${eta_rho// /}

    # Extract mask
    echo '---> Extract mask_rho variable (only interior grid indices, i.e. in fortran convention : 2:end-1)...'
    ncks -O -F -d xi_rho,2,$((${xi_rho}-1)) -d eta_rho,2,$((${eta_rho}-1)) -v mask_rho $filein ${filetmp}

    # set the variable to 0 and rename it var0
    echo '---> Set the variable to 0 and rename it var0...'
    ncap2 -O -v -s "var0=double(mask_rho*0)" ${filetmp} ${filetmp}

elif  [ $model == ww3 ] ; then
    
    # Extract mask
    echo '---> Extract MAPSTA variable...'
    ncks -O -3 -v MAPSTA $filein ${filetmp}

    # set the variable to 0 and rename it var0
    echo '---> Set the variable to 0 and rename it var0...'
    ncap2 -O -v -s "var0=double(MAPSTA*0)" ${filetmp} ${filetmp}

elif  [ $model == toy ] ; then

    # Extract mask
    echo '---> Extract imask_t variable...'
    ncks -O -3 -v imask_t $filein ${filetmp}

    # set the variable to 0 and rename it var0
    echo '---> Set the variable to 0 and rename it var0...'
    ncap2 -O -v -s "var0=double(imask_t*0)" ${filetmp} ${filetmp}

else
    
    echo 'ERROR: '$model' case is not implemented yet. Exit...'
    echo ' '
    exit 1

fi # model 
    
# START LOOP on varlist #
#-----------------------#
for var in $varlist ; do

    echo ' '
    echo '================================'
    echo 'Create variable: '$var'...'
    echo '================================'
    ncks -A -v var0 ${filetmp} $fileout
    ncrename -v var0,$var $fileout

done
rm ${filetmp}

# remove global attributes
ncatted -h -O -a ,global,d,, $fileout $fileout

echo 'DONE for '$varlist' variables => initialized to 0'
echo ' '
