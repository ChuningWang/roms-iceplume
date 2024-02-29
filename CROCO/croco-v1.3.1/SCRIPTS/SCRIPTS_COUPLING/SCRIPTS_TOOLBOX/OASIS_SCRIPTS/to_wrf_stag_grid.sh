#!/bin/bash -e

## ----------------------------------------------------------------------------- #
## - Function to put variables on WRF stag grid                                - #
##                                                                             - #
## Mandatory inputs:                                                           - #
##  - input file (from WRF with full path)                                     - #
##  - output file (with full path)                                             - #
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

# ------------------------------------------------------------------------------ #

echo '**************************'
echo 'ENTER to_wrf_stag_grid.sh'
echo '**************************'
echo ' '

# First check if inputs are ok
if [[ -z $filein ]] || [[ -z $fileout ]] ; then
    echo 'ERROR: inputs are not correctly specified.'
    echo '       this function needs 2 inputs:'
    echo '       - input file (from WRF with full path) '
    echo '       - output file (with full path)  '
    echo ' Exit...'
    echo ' '
    exit 1
fi 
# ------------------------------------------------------------------------------ #

mydir=$(dirname "$fileout")
mytmp=$mydir/to_wrf_stag_tmp.nc

    # Extract grid dimensions
    Nlon=`ncdump -h $filein | grep "west_east = " | cut -d '=' -f2 | cut -d ';' -f1`
    Nlon=${Nlon// /}
    Nlat=`ncdump -h $filein | grep "south_north = " | cut -d '=' -f2 | cut -d ';' -f1`
    Nlat=${Nlat// /}
    echo '---> Grid dimensions are:'
    echo 'west_east = '$Nlon
    echo 'south_north = '$Nlat
    echo ' ' 
    
    # need one more latitude: repeat the last T point
    echo '---> need one more latitude : repeat the last T point...'
    ncpdq -O -a south_north,Time $filein $mytmp
    ncks -F -O -d south_north,$Nlat $mytmp ${mytmp}2
    ncrcat -O $mytmp ${mytmp}2 $mytmp
    ncpdq -O -a Time,south_north $mytmp $mytmp
    # need one more longitude: repeat the last T point
    echo '---> need one more longitude: repeat the last T point...'
    ncpdq -O -a west_east,Time $mytmp $mytmp
    ncks -F -O -d west_east,$Nlon $mytmp ${mytmp}2
    ncrcat -O $mytmp ${mytmp}2 $mytmp
    ncpdq -O -a Time,west_east $mytmp $fileout

    rm $mytmp ${mytmp}2
 
