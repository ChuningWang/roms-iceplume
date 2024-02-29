#!/bin/bash -e

## ----------------------------------------------------------------------------- #
## - Function to extract or compute coupled variables from WW3                 - #
##                                                                             - #
## Mandatory inputs:                                                           - #
##  - the input file (from WW3 with full path)                                 - #
##  - the output file (with full path)                                         - #
##  - the coupled variable to compute                                          - #
##  - the time indice or range on which to extract the variable (F convention) - #
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
var=$3
timerange=$4

# ------------------------------------------------------------------------------ #

echo '******************'
echo 'ENTER from_ww3.sh'
echo '******************'
echo ' '

# First check if inputs are ok
if [[ -z $filein ]] || [[ -z $fileout ]] || [[ -z $var ]] || [[ -z $timerange ]] ; then
    echo 'ERROR: inputs are not correctly specified.'
    echo '       this function needs at least 4 inputs (and max 5):'
    echo '       - the input file (from WW3 with full path)'
    echo '       - the output file (with full path) '
    echo '       - the coupled variable to compute'
    echo '       - the time indice or range on which to extract the variable (F convention)'
    echo ' Exit...'
    echo ' '
    exit 1
fi
# ------------------------------------------------------------------------------ #

mydir=$(dirname "$fileout")
mytmp=$mydir/from_ww3_tmp.nc

    # Model variables
    if [ $var == WW3_ACHA ] ; then
      varin=MAPSTA
    elif [ $var == WW3__DIR ] ; then
      varin=dir
    elif [ $var == WW3__OHS ] ; then
      varin=hs
    elif [ $var == WW3_T0M1 ] ; then
      varin=t0m1
    elif [ $var == WW3_TWOX ] ; then
      varin=utwo
    elif [ $var == WW3_TWOY ] ; then
      varin=vtwo
    elif [ $var == WW3_TAWX ] ; then
      varin=utaw
    elif [ $var == WW3_TAWY ] ; then
      varin=vtaw
    else
      echo 'ERROR: '$var' variable not implemented yet'
      echo 'Exit...'
      echo ' '
      exit 1
    fi

    # Extract variable 
    echo '---> Extract '$varin'...'
    ncks -A -F -d time,$timerange -v $varin $filein $mytmp

    # compute or rename variable
    if [ $var == WW3_ACHA ] ; then
      echo '---> Create charnock coef = 0.0185...'
      ncap2 -A -v -s "WW3_ACHA=MAPSTA*0+0.0185" $mytmp $fileout

#    elif [ $var == WW3_SDIR ] ; then
#      echo '---> Compute sin of dir...'
#      ncap2 -A -v -s "WW3_SDIR=sin((270-dir)*3.1415926/180)" $mytmp $fileout
#
#    elif [ $var == WW3_CDIR ] ; then
#      echo '---> Compute cos of dir...'
#      ncap2 -A -v -s "WW3_CDIR=cos((270-dir)*3.1415926/180)" $mytmp $fileout
#
    else 
      echo '---> Rename '$varin' in '$var
      ncrename -v $varin,$var $mytmp
      ncks -A -v $var $mytmp $fileout
    fi

    rm $mytmp

