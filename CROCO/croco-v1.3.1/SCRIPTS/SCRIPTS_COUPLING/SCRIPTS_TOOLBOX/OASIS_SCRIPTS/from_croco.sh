#!/bin/bash -e

## ----------------------------------------------------------------------------- #
## - Function to extract or compute coupled variables from CROCO               - #
##                                                                             - #
## Mandatory inputs:                                                           - #
##  - the input file (from CROCO with full path)                               - #
##  - the output file (with full path)                                         - #
##  - the coupled variable to compute                                          - #
##  - the time indice or range on which to extract the variable (F convention) - #
## Optional input:                                                             - #
##  - the grid levels: 0 for parent, 1 for child, etc                          - # 
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
echo 'ENTER from_croco.sh'
echo '******************'
echo ' '

# First check if inputs are ok
if [[ -z $filein ]] || [[ -z $fileout ]] || [[ -z $var ]] || [[ -z $timerange ]] ; then
    echo 'ERROR: inputs are not correctly specified.'
    echo '       this function needs at least 4 inputs (and max 5):'
    echo '       - the input file (from CROCO with full path)'
    echo '       - the output file (with full path) '
    echo '       - the coupled variable to compute'
    echo '       - the time indice or range on which to extract the variable (F convention)'
    echo ' Exit...'
    echo ' '
    exit 1
fi
# ------------------------------------------------------------------------------ #

mydir=$(dirname "$fileout")
mytmp=$mydir/from_croco_tmp.nc

    if [ -z $gridlevels ] ; then
      echo 'Default grid level is assumed: 0 for parent...'
      gridlevels=''
    fi

    # Model variables
    if [ $var == CROCO_UOCE$gridlevels ] || [ $var == CROCO_EOCE$gridlevels ] ; then
      varin=u
    elif [ $var == CROCO_VOCE$gridlevels ] || [ $var == CROCO_NOCE$gridlevels ] ; then
      varin=v
    elif [ $var == CROCO_SSH$gridlevels ] ; then
      varin=zeta
    elif [ $var == CROCO_SST$gridlevels ] ; then
      varin=temp
    else
      echo 'ERROR: '$var' variable not implemented yet'
      echo 'Exit...'
      echo ' '
      exit 1
    fi

    # Extract number of vertical levels
    Ns_rho=`ncdump -h $filein | grep "s_rho = " | cut -d '=' -f2 | cut -d ';' -f1`
    Ns_rho=${Ns_rho// /}

    # Extract dimensions
    xi_rho=`ncdump -h $filein | grep "xi_rho = " | cut -d '=' -f2 | cut -d ';' -f1`
    xi_rho=${xi_rho// /}
    eta_rho=`ncdump -h $filein | grep "eta_rho = " | cut -d '=' -f2 | cut -d ';' -f1`
    eta_rho=${eta_rho// /}
    xi_u=`ncdump -h $filein | grep "xi_u = " | cut -d '=' -f2 | cut -d ';' -f1`
    xi_u=${xi_u// /}
    eta_u=`ncdump -h $filein | grep "eta_u = " | cut -d '=' -f2 | cut -d ';' -f1`
    eta_u=${eta_u// /}
    xi_v=`ncdump -h $filein | grep "xi_v = " | cut -d '=' -f2 | cut -d ';' -f1`
    xi_v=${xi_v// /}
    eta_v=`ncdump -h $filein | grep "eta_v = " | cut -d '=' -f2 | cut -d ';' -f1`
    eta_v=${eta_v// /}


    if [[ $var == "CROCO_EOCE${gridlevels}" ]] || [[ $var == "CROCO_NOCE${gridlevels}" ]]; then
        ll="u v angle"
        for vartmp in ${ll};do
            [[ ${vartmp} == "u" ]] && { dim1="xi_u,xi" ;dim2="eta_rho,eta" ;}
            [[ ${vartmp} == "v" ]] && { dim1="xi_rho,xi" ;dim2="eta_v,eta" ;}
            [[ ${vartmp} == "angle" ]] && { dim1="xi_rho,xi" ;dim2="eta_rho,eta" ;}

            ncks -O -F -d s_rho,$Ns_rho -d time,$timerange \
                   -d xi_rho,2,$((${xi_rho}-1)) -d eta_rho,2,$((${eta_rho}-1)) \
                   -d xi_u,1,$((${xi_u}-1))     -d eta_v,1,$((${eta_v}-1)) \
                   -v $vartmp $filein ${mytmp}.$vartmp
            [[ $vartmp == "angle" ]] && cp ${mytmp}.$vartmp angle.nc
            ncrename -O -d ${dim1} -d ${dim2} ${mytmp}.$vartmp ${mytmp}.$vartmp 
            ncap2 -O -v -s "${vartmp}=${vartmp}" ${mytmp}.$vartmp ${mytmp}.$vartmp
            [[ $vartmp == "angle" ]] && cp ${mytmp}.$vartmp angle.nc.2
            ncks -A ${mytmp}.${vartmp} ${mytmp}
            rm -f ${mytmp}.${vartmp}
        done
     
        if [[ $var == "CROCO_EOCE${gridlevels}" ]]; then
            ncap2 -O -s "u=u*cos(angle)-v*sin(angle)" ${mytmp} ${mytmp} 
            ncrename -O -d eta,eta_rho -d xi,xi_u ${mytmp} ${mytmp} 
        elif [[ $var == "CROCO_NOCE${gridlevels}" ]]; then
            ncap2 -O -s "v=u*sin(angle)+v*cos(angle)" ${mytmp} ${mytmp}
            ncrename -O -d xi,xi_rho -d eta,eta_v ${mytmp} ${mytmp}
        fi
    else
        ncks -O -F -d s_rho,$Ns_rho -d time,$timerange \
             -d xi_rho,2,$((${xi_rho}-1)) -d eta_rho,2,$((${eta_rho}-1)) \
             -d xi_u,1,$((${xi_u}-1))     -d eta_v,1,$((${eta_v}-1)) \
             -v $varin $filein $mytmp
    fi

    # remove vertical dimension
    if [ $varin != zeta ]; then
      echo '---> Remove vertical dimension...'
      ncwa -O -a s_rho $mytmp $mytmp 
    fi

    # rename variable
    echo '---> Rename '$varin' in '$var
    ncrename -v $varin,$var $mytmp
    ncks -A -v $var $mytmp $fileout

    rm $mytmp
