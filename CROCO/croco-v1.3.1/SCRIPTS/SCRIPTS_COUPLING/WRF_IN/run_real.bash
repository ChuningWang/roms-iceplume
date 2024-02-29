#!/bin/bash 

# --------------------------------------------------
#
#  Script to prepare WRF using the WRF real program
#
# --------------------------------------------------
#
# Run : ./run_real.bash configure.namelist.real NBPROCS
#
# Usage : 
#  - multiple nested domains (3 max currently)
# 
# Dependence : 
#  - need namelist.input.base
#  - read "configure.namelist.real" defining basic domain and 
#    run parameters
#  - Vtable : Vtable.AMIP or Vtable.GFS
#  - Vtable : Vtable.SSTROMS     
# 
# Source files: outputs from WPS 
#  - met_em.d...
#
# --------------------------------------------------
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
# adapted from J. Lefevre, IRD Noumea, NC
# and from S/ Masson, LOCEAN, Paris, France
# --------------------------------------------------

CONFIGURE_NAMELIST=$1
NBPROCS=$2

source $CONFIGURE_NAMELIST

#=========================================================================
#=======================  USER CHANGES  ==================================
#=========================================================================

#============= Define OPTIONS for running the script =================
#
# I/O parameters
#  his_interval_h  : history data interval time  [hours]
his_interval_h=3
#  his_frames      : nb of history records per file
his_frames=240
#  diag_int_m      : diagnoses interval output [min]
diag_int_m=180
#  diag_frames     : nb of diagnoses records per file
diag_frames=240
#  restart_flag    : .true. or .false.
rst="false"
#  rst_interval_h  : restart data interval       [hours]
rst_interval_h=24
#
# SWITCH fdda # 1 : grid nudging ; 2 : spectral nudging ; 0 : no nudging
switch_fdda=0
#
nudgedom="1"  # "1" is for the parent. For Parent+ nest "1 1", etc...
nudge_coef="0.0003" # Same than for nudgedom. Also need to have the same dimension
nudge_interval_m="360"
nudge_end_h="144"

#
#
#============= Set Environment ======================================= 
#
source ../myenv_mypath.sh
#
# Metgrid outputs path
export O_DATAROOT="${ATM_FILES_DIR}/WPS_DATA"
#
# real exe path
export REAL_EXE_DIR="${ATM}/exe_uncoupled"
#
# workdir path
export REAL_WORK_DIR="${CWORK}/outputs_real"
#
# real outputs path
export REAL_OUT_DIR="${ATM_FILES_DIR}"
#
# wrf in dir
export WRF_IN_DIR="${ATM_NAM_DIR}"
#
if [ -e $REAL_WORK_DIR ] ; then
 rm -f $REAL_WORK_DIR/*
else
 mkdir -p $REAL_WORK_DIR
fi
# MPI launch commands
# ------------------
if [ ${MACHINE} == "JEANZAY" ]; then
    export myMPI="srun -n $NBPROCS "
elif [ ${MACHINE} == "DATARMOR" ]; then
    export myMPI="$MPI_LAUNCH -np $NBPROCS "
elif [ ${MACHINE} == "IRENE" ]; then
    export myMPI="ccc_mprun -n $NBPROCS "
else
    echo "Define how to run the job in run_wps.bash"
    exit
fi
#
cp namelist.input.base.complete ${REAL_WORK_DIR}/.
# MPI parameters
nprocX=-1
nprocY=-1
niowrf=0
niogp=1
#

#=========================================================================
#=======================  END USER CHANGES  ==============================
#=========================================================================

echo "  "
echo "==============================="
echo "REAL_EXE_DIR=> "$REAL_EXE_DIR
echo "==============================="

echo "  "
echo "==============================="
echo "O_DATAROOT=> "$O_DATAROOT
echo "==============================="

echo "  "
echo "==============================="
echo "REAL_WORK_DIR=> "$REAL_WORK_DIR
echo "==============================="

echo "  "
echo "==============================="
echo "REAL_OUT_DIR=> "$REAL_OUT_DIR
echo "==============================="

echo " "
echo "************************************************************"
echo "*                                                          *"
echo "*                This script will run the                  *"
echo "*                    WRF real program                      *"
echo "*  (P. Marchesiello, J. Lefevre, F. Lemarie, S. Jullien)   *"
echo "*                                                          *"
echo "************************************************************"

#----------------------------------------------
# Rescale parameters from configure.namelist
#----------------------------------------------
export interval_s=`expr $obc_freq_h \* 3600`
export dx_d01=`echo "$dx*1000/1" | bc`
export dx_d02=`echo "$dx*1000/$refine_d02" | bc`
export dx_d03=`echo "$dx_d02/$refine_d03" | bc`
#export dx_d01=`expr $dx \* 1000`       # dx in meters
#export dx_d02=`expr $dx \* 1000 / $refine_d02`
#export dx_d03=`expr $dx_d02 \/ $refine_d03`

#------------------------------------------
# Create output data directory if needed
#-------------------------------------------

if ! [ -e $REAL_OUT_DIR ] ; then
  mkdir -p $REAL_OUT_DIR
fi

# O  O  O  O  O  O  O  O  O        START real.exe        O  O  O  O  O  O  O  O  O

  echo "         "
  echo "   BASH: "
  echo " =========================================================="
  echo "   REAL.EXE                                                "
  echo " =========================================================="
  echo "         "

  #************************
  # Check real.exe        *
  #************************
  echo "chdir $REAL_EXE_DIR"
  cd $REAL_EXE_DIR

  if [ -e real.exe ] ; then
    echo "   Executables are available "
    echo " "
    echo "cp real.exe $REAL_WORK_DIR"

#LEFEVRE J 2010-jul : 
# trouble avec real 3.2 et AMIP cdf2grib     
  cp real.exe $REAL_WORK_DIR
# cp /home/jlefevre/WRF/WRFV2.2/main/real.exe $REAL_EXE_DIR

  else
    echo "   Executables real.exe is missing ... "
    exit 1
  fi
  echo "chdir $REAL_WORK_DIR"
  cd $REAL_WORK_DIR
  echo `pwd`

  #**********************************************************
  # Create initial and boundary conditions for each domain  *
  #**********************************************************

  rm -f wrfinput* wrfbdy* wrflowinp* rsl.* met_em.*
  if [ $switch_fdda -ne 0 ]; then
   rm -f wrffdda*
   if  [ $switch_fdda -eq 1 ]; then
    echo "    Grid nudging activated     "
   else
    echo "    Spectral nudging activated "
   fi
  else
   echo "     No nudging activated "
  fi

  #----------------------------------------------
  # Rescale parameters from configure.namelist
  #----------------------------------------------
  export interval_s=`expr $obc_freq_h \* 3600`
  export rst_interval_m=`expr $rst_interval_h \* 60` # refresh interval in minutes

  for metfile in `ls $O_DATAROOT/met_em.*`
  do
    link_target=`basename $metfile`
    echo "Link $metfile to $link_target"           
    ln -s $metfile $link_target
  done

    echo " "
    echo "   First create namelist.input.prep     "
    echo " "
    if [ -e namelist.input ] ; then
      rm -f namelist.input
    fi

sed -e "s/<interval_s>/${interval_s}/g"         \
    -e "s/<sst_int_m>/${sst_interval_m}/g"      \
    -e "s/<max_domains>/${max_domains}/g" -e "s/<max_cpldom>/0/g"     \
    -e "s/<nbvertlev>/${nbvertlevel}/g"   -e  "s/<eta_lev>/${eta_levels}/g" -e "s/<ptop>/${ptop}/g"        \
    -e "s/<nbmetlev>/${nbmetlevel}/g"     -e "s/<nbmetsoil>/${nbmetsoil}/g"                                \
    -e "s/<xdim_d01>/${xdim_d01}/g"       -e "s/<xdim_d02>/${xdim_d02}/g"  -e "s/<xdim_d03>/${xdim_d03}/g" \
    -e "s/<ydim_d01>/${ydim_d01}/g"       -e "s/<ydim_d02>/${ydim_d02}/g"  -e "s/<ydim_d03>/${ydim_d03}/g" \
    -e "s/<dx_d01>/${dx_d01}/g"           -e "s/<dx_d02>/${dx_d02}/g"      -e "s/<dx_d03>/${dx_d03}/g"     \
    -e "s/<dy_d01>/${dx_d01}/g"           -e "s/<dy_d02>/${dx_d02}/g"      -e "s/<dy_d03>/${dx_d03}/g"     \
    -e "s/<i_str_d02>/${i_str_d02}/g"     -e "s/<i_str_d03>/${i_str_d03}/g"                                \
    -e "s/<j_str_d02>/${j_str_d02}/g"     -e "s/<j_str_d03>/${j_str_d03}/g"                                \
    -e "s/<coef_d02>/${refine_d02}/g"     -e "s/<coef_d03>/${refine_d03}/g"                                \
    -e "s/<isftcflx>/${isftcflx}/g" \
    namelist.input.base.complete > namelist.input.prep.${domain_name}

# Handle fdda for different domains
if [[ ${switch_fdda} != 0 ]]; then
    nbdom=$( echo "${nudgedom}" | wc -w)
    [[ ${nbdom} >  $( echo "${nudge_coef}" | wc -w) ]] && { echo "Missing values in nudge_coef for nest, we stop..."; exit ;}
    [[ ${nbdom} >  $( echo "${nudge_interval_m}" | wc -w) ]] && { echo "Missing values in nudge_interval_m for nest, we stop..."; exit ;}
    [[ ${nbdom} >  $( echo "${nudge_end_h}" | wc -w) ]] && { echo "Missing values in nudge_end_h for nest, we stop..."; exit ;}

    for lvl in `seq 1 ${nbdom}`; do
        dom=$(printf "d%02d" ${lvl})
        sed -e "s/<nudge_${dom}>/$( echo "${nudgedom}" | cut -d " " -f $(( ${lvl})) )/g"                 \
            -e "s/<nudge_coef_${dom}>/$( echo "${nudge_coef}" | cut -d " " -f $(( ${lvl})) )/g"       \
            -e "s/<nudge_end_h_${dom}>/$( echo "${nudge_end_h}" | cut -d " " -f $(( ${lvl})) )/g"      \
            -e "s/<nudge_int_m_${dom}>/$( echo "${nudge_interval_m}" | cut -d " " -f $(( ${lvl})) )/g" \
            namelist.input.prep.${domain_name} > namelist.tmp
            mv namelist.tmp namelist.input.prep.${domain_name}
   done

    sed -e "s/<nudge_d.*>/0/g" \
        -e "s/<nudge_coef_d.*>/0/g" \
        -e "s/<nudge_end_h_d.*>/0/g" \
        -e "s/<nudge_int_m_d.*>/0/g" \
        namelist.input.prep.${domain_name} > namelist.tmp
    mv namelist.tmp namelist.input.prep.${domain_name}
    chmod 755 namelist.input.prep.${domain_name}
else
    sed -e "s/<nudge_d.*>/0/g" \
        -e "s/<nudge_coef_d.*>/0/g" \
        -e "s/<nudge_end_h_d.*>/0/g" \
        -e "s/<nudge_int_m_d.*>/0/g" \
        namelist.input.prep.${domain_name} > namelist.tmp
    mv namelist.tmp namelist.input.prep.${domain_name}
    chmod 755 namelist.input.prep.${domain_name}
fi  
 
    echo " "
    echo "   Then create namelist.input adding run, time and output settings     "
    echo " "
    lmonth=( 1 3 5 7 8 10 12 )

# Remove heading 0 
while [ `echo ${start_m} | cut -b 1` -eq 0 ]; do
        start_m=`echo ${start_m} | cut -b 2-`
done
while [ `echo ${end_m} | cut -b 1` -eq 0 ]; do
        end_m=`echo ${end_m} | cut -b 2-`
done
#
for yy in `seq $start_y $end_y`; do
    [[ $yy == $start_y ]] && { mstart=$start_m ;} || { mstart=1 ;}
    [[ $yy == $end_y ]] && { mend=$end_m ;} || { mend=12 ;}

    for mm in `seq $mstart $mend`; do      
        [[ $yy == $start_y && $mm == $start_m ]] && { sday=$start_d ; shour=$start_h;} || { sday=01 ; shour=00;}

	if [[ $yy == $end_y && $mm == $end_m ]] ; then 
            eday=$end_d ; ehour=$end_h
            yy2=$yy
            emth=$mm
        elif [[ $mm == 12 ]]; then
            eday=01 ; ehour=00;
            yy2=$(( $yy + 1 ))
            emth=$(( $mm + 1 - 12 ))
        else
            eday=01 ; ehour=00;
            yy2=$yy
            emth=$(( $mm + 1 ))
        fi
	
        if [ $yy == $end_y ] && [ $mm == $end_m ] && [ $eday == $sday ] && [ $ehour == $shour ] ; then
                exit
        fi
 
        sed -e "s/<yr1>/$yy/g"   -e "s/<yr2>/$yy2/g"  \
            -e "s/<mo1>/$mm/g"   -e "s/<mo2>/$emth/g"  \
            -e "s/<dy1>/$sday/g"   -e "s/<dy2>/$eday/g"  \
            -e "s/<hr1>/$shour/g"   -e "s/<hr2>/$ehour/g"  \
            -e "s/<rst>/$rst/g"                    -e "s/<rst_int_h>/$rst_interval_h/g"   \
            -e "s/<his_int_h>/${his_interval_h}/g" -e "s/<his_nb_out>/${his_frames}/g"    \
            -e " s/<xtrm_int_m>/${diag_int_m}/g"    -e "s/<xtrm_nb_out>/${diag_frames}/g"  \
            -e "s/<nproc_x>/$nprocX/g"             -e "s/<nproc_y>/$nprocY/g"             \
            -e "s/<niotaskpg>/$niowrf/g"           -e "s/<niogp>/$niogp/g"                \
            -e "s/<dt>/${dt}/g"  \
            namelist.input.prep.${domain_name} > namelist.input

        chmod +x namelist.input
        cp namelist.input namelist.input.real.${domain_name}
        cp namelist.input.prep.${domain_name} $WRF_IN_DIR
        cp namelist.input.real.${domain_name} $WRF_IN_DIR 
 # RUN real.exe
        echo "   Run real.exe    "
        date
        ${myMPI}real.exe >& real.log
        date

        ls -l *

# Finalize with real outputs, rename and mov
        tdigit=$( printf "%02d" $mm)
        if [ -e wrfinput_d01 -a -e wrfbdy_d01 -a -e wrflowinp_d01 ] ; then
            echo "SUCCESS d01"
        # Create output data directory if needed
            if ! [ -e ${REAL_OUT_DIR} ] ; then
                mkdir ${REAL_OUT_DIR}
            fi 
            mv -f wrfinput_d01 ${REAL_OUT_DIR}/wrfinput_d01_${yy}_${tdigit}
            mv -f wrfbdy_d01 ${REAL_OUT_DIR}/wrfbdy_d01_${yy}_${tdigit}
            mv -f wrflowinp_d01 ${REAL_OUT_DIR}/wrflowinp_d01_${yy}_${tdigit}
            if  [ $switch_fdda -ne 0 ]; then
                mv -f wrffdda_d01 ${REAL_OUT_DIR}/wrffdda_d01_${yy}_${tdigit}
            fi  
        else
            echo "REAL ERROR d01"
        fi

        if [ $max_domains -ge 2 -a -e wrfinput_d02 -a -e wrflowinp_d02 ] ; then
            echo "SUCCESS d02"
            mv -f wrfinput_d02 ${REAL_OUT_DIR}/wrfinput_d02_${yy}_${tdigit}
            mv -f wrflowinp_d02 ${REAL_OUT_DIR}/wrflowinp_d02_${yy}_${tdigit}
        else
            echo "REAL ERROR or NON-EXISTENT d02"
        fi
        if [ $max_domains -eq 3 -a -e wrfinput_d03 -a -e wrflowinp_d03 ] ; then
            echo "SUCCESS d03"
            mv -f wrfinput_d03 ${REAL_OUT_DIR}/wrfinput_d03_${yy}_${tdigit}
            mv -f wrflowinp_d03 ${REAL_OUT_DIR}/wrflowinp_d03_${yy}_${tdigit}
        else
            echo "REAL ERROR or NON-EXISTENT d03"
        fi
    done
done
# O  O  O  O  O  O  O  O  O        END real.exe        O  O  O  O  O  O  O  O  O

