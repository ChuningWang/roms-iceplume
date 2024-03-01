#!/bin/bash -e

## ----------------------------------------------------------------------------- #
## - Function to extract or compute coupled variables from WRF                 - #
##                                                                             - #
## Mandatory inputs:                                                           - #
##  - the input file (from WRF with full path)                                 - #
##  - the output file (with full path)                                         - #
##  - the coupled variable to compute                                          - #
##  - the time indice or range on which to extract the variable (F convention) - #
## Optional input:                                                             - #
##  - the grid levels:                                                         - # 
##                       examples: WRF_d01_EXT_d01 or WRF_d02_EXT_d01          - #
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
var=$3
timerange=$4
gridlevels=$5

# ------------------------------------------------------------------------------ #

echo '******************'
echo 'ENTER from_wrf.sh'
echo '******************'
echo ' '

# First check if inputs are ok
if [[ -z $filein ]] || [[ -z $fileout ]] || [[ -z $var ]] || [[ -z $timerange ]] ; then
    echo 'ERROR: inputs are not correctly specified.'
    echo '       this function needs at least 4 inputs (and max 5):'
    echo '       - the input file (from WRF with full path)'
    echo '       - the output file (with full path) '
    echo '       - the coupled variable to compute'
    echo '       - the time indice or range on which to extract the variable (F convention)'
    echo ' Exit...'
    echo ' '
    exit 1
fi 
# ------------------------------------------------------------------------------ #

mydir=$(dirname "$fileout")
mytmp=$mydir/from_wrf_tmp.nc

    if [ -z $gridlevels ] ; then
      echo 'Default grid levels are assumed: WRF_d01_EXT_d01...'
      gridlevels='WRF_d01_EXT_d01'
    fi

    # Model variables 
    if [ $var == ${gridlevels}_TAUX ] ; then
      varin=UTAU
    elif [ $var == ${gridlevels}_TAUY ] ; then
      varin=VTAU
    elif [ $var == ${gridlevels}_TAUMOD ] ; then
      varin=TAUM
    elif [ $var == ${gridlevels}_TAUE ] ; then
      varin=UTAU,VTAU,COSALPHA,SINALPHA
    elif [ $var == ${gridlevels}_TAUN ] ; then
      varin=UTAU,VTAU,COSALPHA,SINALPHA
    elif [ $var == ${gridlevels}_WND_U_01 ] ; then
#      varin=U_01
      varin=U10
    elif [ $var == ${gridlevels}_WND_V_01 ] ; then
#      varin=V_01
      varin=V10
    elif [ $var == ${gridlevels}_WND_E_01 ] ; then
      varin=U10,V10,COSALPHA,SINALPHA
    elif [ $var == ${gridlevels}_WND_N_01 ] ; then
      varin=U10,V10,COSALPHA,SINALPHA
    elif [ $var == ${gridlevels}_SURF_NET_SOLAR ] ; then
      varin=SWDOWN
    elif [ $var == ${gridlevels}_SURF_NET_NON-SOLAR ] ; then
      varin=HFX,LH,GLW,SST
    elif [ $var == ${gridlevels}_EVAP-PRECIP ] ; then
      varin=QFX,RAINC,RAINNC,XTIME
    elif [ $var == ${gridlevels}_PSFC ] ; then
      varin=PSFC
    else
      echo 'ERROR: '$var' variable not implemented yet'
      echo 'Exit...'
      echo ' '
      exit 1
    fi

    echo ' '
    echo '==================='
    echo 'Extract '$varin 
    echo '==================='

    ncks -F -O -v $varin -d Time,$timerange $filein $mytmp

    # rename or compute variable
    if [ $var == ${gridlevels}_SURF_NET_SOLAR ] ; then
      echo '---> Compute variable: '$var'...'
      ncap2 -A -v -s "${gridlevels}_SURF_NET_SOLAR=SWDOWN*(1-0.08)" $mytmp $fileout

    elif [ $var == ${gridlevels}_SURF_NET_NON-SOLAR ] ; then
      echo '---> Compute variable: '$var'...'
      ncap2 -O -s "LW_out=5.67*10^(-8)*0.985*SST^4" $mytmp $mytmp
      ncap2 -A -v -s "${gridlevels}_SURF_NET_NON=GLW-LW_out-LH-HFX" $mytmp $mytmp
      ncrename -v ${gridlevels}_SURF_NET_NON,${gridlevels}_SURF_NET_NON-SOLAR $mytmp
      ncks -A -v ${gridlevels}_SURF_NET_NON-SOLAR $mytmp $fileout

    elif  [ $var == ${gridlevels}_EVAP-PRECIP ] ; then
      echo '---> Compute variable: '$var'...'
      ncap2 -O -v -s "${gridlevels}_EP=QFX-(RAINC+RAINNC)/(XTIME*60)" $mytmp $mytmp
      ncrename -v ${gridlevels}_EP,${gridlevels}_EVAP-PRECIP $mytmp
      ncks -A -v ${gridlevels}_EVAP-PRECIP $mytmp $fileout

    elif [ $var == ${gridlevels}_TAUE ] ; then
      echo '---> Compute variable: '$var'...'
      ncap2 -O -v -s "${gridlevels}_TE=UTAU*COSALPHA-VTAU*SINALPHA" $mytmp $mytmp
      ncrename -v ${gridlevels}_TE,${gridlevels}_TAUE $mytmp
      ncks -A -v ${gridlevels}_TAUE $mytmp $fileout

    elif [ $var == ${gridlevels}_TAUN ] ; then
      echo '---> Compute variable: '$var'...'
      ncap2 -O -v -s "${gridlevels}_TN=VTAU*COSALPHA+UTAU*SINALPHA" $mytmp $mytmp
      ncrename -v ${gridlevels}_TN,${gridlevels}_TAUN $mytmp
      ncks -A -v ${gridlevels}_TAUN $mytmp $fileout

    elif [ $var == ${gridlevels}_WND_E_01 ] ; then
      echo '---> Compute variable: '$var'...'
      ncap2 -O -v -s "${gridlevels}_WE=U10*COSALPHA-V10*SINALPHA" $mytmp $mytmp
      ncrename -v ${gridlevels}_WE,${gridlevels}_WND_E_01 $mytmp
      ncks -A -v ${gridlevels}_WND_E_01 $mytmp $fileout

    elif [ $var == ${gridlevels}_WND_N_01 ] ; then
      echo '---> Compute variable: '$var'...'
      ncap2 -O -v -s "${gridlevels}_WN=V10*COSALPHA+U10*SINALPHA" $mytmp $mytmp
      ncrename -v ${gridlevels}_WN,${gridlevels}_WND_N_01 $mytmp
      ncks -A -v ${gridlevels}_WND_N_01 $mytmp $fileout

    else
      echo '---> Rename variable: '$varin' to '$var
      ncrename -v $varin,$var $mytmp
      ncks -A -v $var $mytmp $fileout

    fi
    
    rm $mytmp    
