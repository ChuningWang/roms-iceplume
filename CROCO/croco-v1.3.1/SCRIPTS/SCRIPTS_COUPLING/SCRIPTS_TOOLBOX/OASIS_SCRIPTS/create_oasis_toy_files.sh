#!/bin/bash -e

## ----------------------------------------------------------------------------- #
## - Create data files for oasis toy                                           - #
## - from pre-existing model files                                             - #
##                                                                             - #
## Mandatory inputs:                                                           - #
##  - the model file name (with full path)                                     - #
##  - the data file name for the toy model                                     - #
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
timerange=$4
gridlevels=$5

# ------------------------------------------------------------------------------ #
echo '**************************************'
echo 'START script create_oasis_toy_files.sh'
echo '**************************************'
echo ' '

# First check if inputs are ok
if [[ -z $filein ]] || [[ -z $fileout ]] || [[ -z $timerange ]] || [[ -z $model ]] ; then
    echo 'ERROR: inputs are not correctly specified.'
    echo '       this script needs at least 4 inputs:'
    echo '       - the model file name (with full path)'
    echo '       - the data file name for the toy model'
    echo '       - the model: wrf croco or ww3 cases are accepted'
    echo '       - the time range to extract'
    echo ' Exit...'
    echo ' '
    exit 1
fi 
# ------------------------------------------------------------------------------ #

mydir=$(dirname "$fileout")
toytype=$(basename "$fileout")
toytype=`echo $toytype | cut -c5-7`  
gridtoy=$mydir/grid_${toytype}.nc
filetmp=$mydir/toy_tmp.nc

if [ $model == wrf ] ; then
    echo '==================================='
    echo 'Create toy model in behalf of wrf'
    echo '==================================='

    if [ -z $gridlevels ] ; then
      echo 'Default grid levels are assumed: WRF_d01_EXT_d01...'
      gridlevels='WRF_d01_EXT_d01'
    fi

    varlist=(${gridlevels}_TAUX \
            ${gridlevels}_TAUY \
            ${gridlevels}_TAUMOD \
            ${gridlevels}_WND_U_01 \
            ${gridlevels}_WND_V_01 \
            ${gridlevels}_SURF_NET_SOLAR \
            ${gridlevels}_SURF_NET_NON-SOLAR \
            ${gridlevels}_EVAP-PRECIP \
	    ${gridlevels}_PSFC)

    dimtime=Time

    modellon=XLONG
    modellat=XLAT
    modelmask=LANDMASK
    dimx=west_east
    dimy=south_north

    varlist_toy=(TOY_TAUX \
                 TOY_TAUY \
                 TOY_TAUM \
                 TOY_U_01 \
                 TOY_V_01 \
                 TOYSRFLX \
                 TOYSTFLX \
                 TOY__EMP \
		 TOY_PSFC)

elif  [ $model == croco ] ; then
    echo '==================================='
    echo 'Create toy model in behalf of croco'
    echo '==================================='

    if [ -z $gridlevels ] ; then
      echo 'Default grid level is assumed: 0 for parent...'
      gridlevels=''
    fi

    varlist=(CROCO_UOCE${gridlevels} \
            CROCO_VOCE${gridlevels} \
            CROCO_SST${gridlevels} \
            CROCO_SSH${gridlevels})

    dimtime=time

    modellon=lon_rho
    modellat=lat_rho
    modelmask=mask_rho
    dimx=xi_rho
    dimy=eta_rho

    varlist_toy=(TOY_UOCE \
                 TOY_VOCE \
                 TOY__SST \
                 TOY__SSH)

elif  [ $model == ww3 ] ; then
    echo '==================================='
    echo 'Create toy model in behalf of ww3'
    echo '==================================='

    varlist=(WW3_T0M1 \
            WW3__OHS \
            WW3__DIR \
            WW3_TWOX \
            WW3_TWOY \
            WW3_TAWX \
            WW3_TAWY \
            WW3_ACHA)

    dimtime=time

    modellon=longitude
    modellat=latitude
    modelmask=MAPSTA
    dimx=longitude
    dimy=latitude

    varlist_toy=(TOY_T0M1 \
                TOY___HS \
                TOY__DIR \
                TOY_TWOX \
                TOY_TWOY \
                TOY_TAWX \
                TOY_TAWY \
                TOY__CHA)
else
    echo 'ERROR: '$model' case is not implemented yet. Exit...'
    echo ' ' 
    exit 1
fi # model 

echo 'Varlist to proceed is '${varlist[*]}
echo '===================== '
lengthvar=${#varlist[@]}
for k in `seq 0 $(( ${lengthvar} - 1))` ; do

    var=${varlist[$k]}
    vartoy=${varlist_toy[$k]}
    echo ' '
    echo '====================='
    echo 'Process '$var'...'
    echo '====================='

    # Extract or compute var
    echo '---> Extract or compute '$var
    ${SCRIPTDIR}/OASIS_SCRIPTS/from_${model}.sh $filein $filetmp $var $timerange $gridlevels

    if [ $model == wrf ] ; then
        # Put them on the stag grid
        echo '---> Put them on the stag grid' 
        ${SCRIPTDIR}/OASIS_SCRIPTS/to_wrf_stag_grid.sh $filetmp $filetmp
        # rename dimensions
        echo '---> Rename dimensions...'
        # problem with some NCO versions, need to be in netcdf3
        ncks -O --3 $filetmp $filetmp
        ncrename -d $dimx,nlon -d $dimy,nlat -d $dimtime,time $filetmp
        ncks -O --4 $filetmp $filetmp
    fi

    ncks -A -v $var $filetmp $fileout
    rm $filetmp
    
    ncrename -v $var,$vartoy $fileout   
done # LOOP on varlist
# remove global attributes
ncatted -h -O -a ,global,d,, $fileout $fileout

echo ' '
echo '========================================================='
echo 'Create grid file '$gridtoy' for toy model'
echo '========================================================='

# Extract lon, lat, mask
echo '---> Extract lon, lat, mask variables...'
ncks -O -v $modellon,$modellat,$modelmask -d $dimtime,0 $filein $gridtoy

if [ $model == wrf ] ; then
    # Put them on the stag grid
    echo '---> Put them on the stag grid' 
    ${SCRIPTDIR}/OASIS_SCRIPTS/to_wrf_stag_grid.sh $gridtoy $gridtoy
fi

# remove time dimension
#echo '---> Remove time dimension...'
#ncwa -O -a $dimtime $gridtoy $gridtoy

# change mask from float to integer
#ncap2 -O -s "${modelmask}=int(${modelmask})" $gridtoy $gridtoy

# rename dimensions
echo '---> Rename dimensions...'
# problem with some NCO versions, need to be in netcdf3
ncks -O --3 $gridtoy $gridtoy
ncrename -d $dimx,nlon $gridtoy
ncrename -d $dimy,nlat $gridtoy
# rename variables
echo '---> Rename variables...'
ncrename -v ${modelmask},imask_t $gridtoy
if [ $model != ww3 ] ; then
 ncrename -v ${modellon},longitude -v ${modellat},latitude $gridtoy
fi
ncks -O --4 $gridtoy $gridtoy

# Modify mask values if necessary
if [ $model == wrf ] ; then
    ncap2 -O -s "imask_t=(imask_t-1)*(-1)" $gridtoy $gridtoy
elif [ $model == ww3 ] ; then
    ncap2 -O -s "where(imask_t != 1) imask_t=0" $gridtoy $gridtoy
elif [ $model == croco ] ; then
    Nnlon=`ncdump -h $gridtoy | grep "nlon = " | cut -d '=' -f2 | cut -d ';' -f1`
    Nnlon=${Nnlon// /}
    Nnlat=`ncdump -h $gridtoy | grep "nlat = " | cut -d '=' -f2 | cut -d ';' -f1`
    Nnlat=${Nnlat// /}
    ncks -O -F -d nlon,2,$((${Nnlon}-1)) -d nlat,2,$((${Nnlat}-1)) $gridtoy $gridtoy
fi
# remove global attributes
ncatted -h -O -a ,global,d,, $gridtoy $gridtoy 

echo '========================================================='
echo 'DONE creating grid file '$gridtoy' for toy model'
echo '========================================================='
echo ' '
