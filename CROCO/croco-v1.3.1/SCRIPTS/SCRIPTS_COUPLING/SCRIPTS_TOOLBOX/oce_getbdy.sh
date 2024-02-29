#!/bin/bash

module load $ncomod

if [[ ${bdy_ext} == *'clm'* ]]; then
    bryfile="croco_clm.nc"
    timevar="tclm_time"
else
    bryfile="croco_bry.nc"
    timevar="bry_time"
fi

# put 1-d stuff inside bdy file
cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
varlist="spherical,Vtransform,Vstretching,tstart,tend,theta_s,theta_b,Tcline,hc,sc_r,sc_w,Cs_r,Cs_w"
ncks -A -v "${varlist}"  ${OCE_FILES_DIR}/croco_${bdy_ext}_Y${cur_Y}M${cur_M}.nc ${bryfile}
#
options=("temp" "salt" "v2d" "v3d"  "zeta")
varlist="ssh tclm sclm uclm vclm temp salt v3d v2d zeta"
[[ ${bdy_ext} == *'bry'* ]] && varlist="bry ${varlist}"

# check if job remains in the same month or not
cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
while [[ `echo ${cur_M} | cut -b 1` -eq 0 ]]; do
    cur_M=`echo ${cur_M} | cut -b 2-`
done

mdy=$( valid_date ${MONTH_END_JOB} $(( ${DAY_END_JOB} + 1 )) ${YEAR_END_JOB} )
LOCAL_MTH_END=$( echo $mdy | cut -d " " -f 1 )


if [[ ${JOB_DUR_MTH} -eq 1 || ${LOCAL_MTH_END} -eq ${cur_M} ]]; then # Case 1 month or less
    echo "Job is one month long or less ---> Using netcdf of the current month"
    cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
    cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
    ln -sf ${OCE_FILES_DIR}/croco_${bdy_ext}_Y${cur_Y}M${cur_M}.nc ${bryfile}
else   
    if [[ ${JOB_DUR_MTH} -eq 0 && ${LOCAL_MTH_END} -ne ${cur_M} ]]; then
        echo "Job is less than a month BUT overlaps on next month ---> Concat netcdf of current and following month"
        nbloop=1
    else
        echo "Job is longer than one month ---> Concat netcdf of needed month"
        nbloop=$(( ${JOB_DUR_MTH}-1 ))     
    fi
    for i in `seq 0 $nbloop`; do
        mdy=$( valid_date $(( $MONTH_BEGIN_JOB + $i )) 1 $YEAR_BEGIN_JOB )
        cur_Y=$( printf "%04d\n"  $( echo $mdy | cut -d " " -f 3) )
        cur_M=$( printf "%02d\n"  $( echo $mdy | cut -d " " -f 1) )

        if [[ $i == 0 ]]; then
            tstart=1; tend=-2
        elif [[ $i == $nbloop ]]; then
            tstart=2; tend=""
        else
            tstart=2; tend=2
        fi

        for var in ${varlist}; do
            if [[ ${options[@]} =~ "$var" ]]; then
                extract=""
                for card in north south east west ;do
                    if [[ ${var} == "v2d" ]];then
                        extract="${extract}ubar_${card},vbar_${card},";
                    elif [[ ${var} == "v3d" ]];then
                        extract="${extract}u_${card},v_${card},";
                    else
                        extract="${extract}${var}_${card},"
                    fi
                done
                cnt=$( echo ${extract} | sed -e 's|\(.\)|\1\n|g' | grep ',' | wc -l )
                extract=$( echo "${extract}" | cut -d ',' -f -"${cnt}" )
            else
                extract="${var}_time"
            fi
            if [[ $i > 0 ]]; then
                ncks --mk_rec_dmn "${var}_time" -F -O -d "${var}_time",${tstart},${tend} -v "${extract}" ${OCE_FILES_DIR}/croco_${bdy_ext}_Y${cur_Y}M${cur_M}.nc tmp_${var}.nc
                ncks -O --mk_rec_dmn "${var}_time" ${bryfile} ${bryfile}
                ncrcat -A ${bryfile} tmp_${var}.nc ${bryfile}
                ncks -O --fix_rec_dmn "${var}_time" ${bryfile} ${bryfile}
            else
                ncks -F -O -d "${var}_time",${tstart},${tend} -v "${extract}" ${OCE_FILES_DIR}/croco_${bdy_ext}_Y${cur_Y}M${cur_M}.nc tmp_${var}.nc
                ncks -A -v "${extract}" tmp_${var}.nc ${bryfile}
            fi
            \rm tmp_${var}.nc
        done
    done
fi
module unload $ncomod
