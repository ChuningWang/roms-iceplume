#!/bin/bash -e

## ----------------------------------------------------------------------------- #
## - Create grids.nc, masks.nc, files from WRF for oasis                       - #
## - because call to oasis_grid function not yet implemented in WRF            - #
##                                                                             - #
## Mandatory inputs:                                                           - #
##  - a file from WRF containing lon,lat,mask (with full path)                 - #
##  - the output destination directory                                         - #
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

gridfile=$1
mydir=$2

# ------------------------------------------------------------------------------ #

echo '*******************************************'
echo 'START script create_oasis_grids_for_wrf.sh'
echo '*******************************************'
echo ' '

# First check if inputs are ok
if [[ -z $gridfile ]] || [[ -z $mydir ]] ; then
    echo 'ERROR: inputs are not correctly specified.'
    echo '       this script needs 2 inputs:'
    echo '       - a file from WRF containing the mask, lon and lat (with full path)'
    echo '       - the output destination directory'
    echo ' Exit...'
    echo ' '
    exit 1
fi
# ------------------------------------------------------------------------------ #

mytmpgrd=$mydir/grd_tmp.nc
grdfile=$mydir/grids.wrf.nc
mskfile=$mydir/masks.wrf.nc

wrflon=XLONG
wrflat=XLAT
wrfmask=LANDMASK

# Extract lon,lat,mask
echo '---> Extract '${wrflon}', '${wrflat}', and '${wrfmask}' variables...'
ncks -O -v ${wrflon},${wrflat},${wrfmask} -d Time,0 $gridfile ${mytmpgrd}

# Put them on the stag grid
echo '---> Put them on the stag grid' 
${SCRIPTDIR}/OASIS_SCRIPTS/to_wrf_stag_grid.sh ${mytmpgrd} ${mytmpgrd}

# remove time dimension
echo '---> Remove time dimension...'
ncwa -O -a Time ${mytmpgrd} ${mytmpgrd}

# compute the last lon and lat
Nlon=`ncdump -h $gridfile | grep "west_east = " | cut -d '=' -f2 | cut -d ';' -f1`
Nlon=${Nlon// /}
Nlat=`ncdump -h $gridfile | grep "south_north = " | cut -d '=' -f2 | cut -d ';' -f1`
Nlat=${Nlat// /}
Nlonstag=$(($Nlon + 1))
Nlatstag=$(($Nlat + 1))
Nlonm1=$(($Nlon - 1))
Nlatm1=$(($Nlat - 1))
echo '---> compute the last lon...'
ncap2 -F -O -s "${wrflon}(:,$Nlonstag)=${wrflon}(:,$Nlon)+(${wrflon}(:,$Nlon)-${wrflon}(:,$Nlonm1))" ${mytmpgrd} ${mytmpgrd}
echo '---> compute the last lat...'
ncap2 -F -O -s "${wrflat}($Nlatstag,:)=${wrflat}($Nlat,:)+(${wrflat}($Nlat,:)-${wrflat}($Nlatm1,:))" ${mytmpgrd} ${mytmpgrd}

# change mask from float to integer
echo '---> Change mask from float to integer...'
ncap2 -O -s "${wrfmask}=int(${wrfmask})" ${mytmpgrd} ${mytmpgrd}
# problem with some NCO versions... Possible fix:
#ncks -O -v ${wrfmask} ${mytmpgrd} ${mytmpgrd}_mask
#ncap2 -O -v -s "tmpmask=int(${wrfmask})" ${mytmpgrd}_mask ${mytmpgrd}_mask
#ncrename -v tmpmask,${wrfmask} ${mytmpgrd}_mask
#ncks -O -x -v ${wrfmask} ${mytmpgrd} ${mytmpgrd}
#ncks -A -v ${wrfmask} ${mytmpgrd}_mask ${mytmpgrd}
#rm ${mytmpgrd}_mask 

# rename dimensions
echo '---> rename dimensions...'
# problem with some NCO versions, need to be in netcdf3
ncks -O --3 ${mytmpgrd} ${mytmpgrd}
ncrename -d west_east,x_atmt -d south_north,y_atmt ${mytmpgrd}
# rename variables
echo '---> Rename variables...'
ncrename -v ${wrfmask},atmt.msk -v ${wrflon},atmt.lon -v ${wrflat},atmt.lat ${mytmpgrd} 
# put the file back to netcdf4
ncks -O --4 ${mytmpgrd} ${mytmpgrd}

# create grid file
echo '---> Create grid file...'
echo '======================='
ncks -O -v atmt.lon,atmt.lat ${mytmpgrd} ${grdfile}
ncatted -h -O -a ,global,d,, ${grdfile} ${grdfile}
ncatted -h -O -a ,atmt.lon,d,, ${grdfile} ${grdfile}
ncatted -h -O -a ,atmt.lat,d,, ${grdfile} ${grdfile}

# create mask file
echo '---> Create mask file...'
echo '========================='
ncks -O -v atmt.msk ${mytmpgrd} ${mskfile}
ncatted -h -O -a ,global,d,, ${mskfile} ${mskfile}
ncatted -h -O -a ,atmt.msk,d,, ${mskfile} ${mskfile}

rm ${mytmpgrd}

echo 'DONE: grids.wrf.nc and masks.wrf.nc have been created in '$mydir
echo ' '
