#!/bin/bash

echo ' '
echo '-- WW3 inputs --------------'
echo 'copy and fill ww3 settings files *.inp'

cp ${WAV_NAM_DIR}/*.inp .

ms=$( printf "%02d"  ${MONTH_BEGIN_JOB} )
me=$( printf "%02d"  ${MONTH_END_JOB} )
ds=$( printf "%02d"  ${DAY_BEGIN_JOB} )
de=$( printf "%02d"  ${DAY_END_JOB} )

 ## - Fill ww3_grid.inp file -##
sed -e "s/<wavdt>/${DT_WAV}/g" \
    -e "s/<wavdtPRO>/${DT_WW_PRO}/g"  -e "s/<wavdtREF>/${DT_WW_REF}/g"  -e "s/<wavdtSRC>/${DT_WW_SRC}/g"  \
    -e "s/<wavnx>/${wavnx}/g"   -e "s/<wavny>/${wavny}/g"  \
    -e "s/<hmin>/${hmin}/g" \
    -e "s/<CEXPER>/${CEXPER}/g" \
    ${WAV_NAM_DIR}/ww3_grid.inp.base > ./ww3_grid.inp

 ## - Fill ww3_ounf.inp file -##
sed -e "s/<wav_int>/${wav_int}/g" \
    -e "s/<yr1>/${YEAR_BEGIN_JOB}/g"  -e "s/<mo1>/${ms}/g"  -e "s/<dy1>/${ds}/g"  -e "s/<hr1>/00/g" \
    ${WAV_NAM_DIR}/ww3_ounf.inp.base > ./ww3_ounf.inp

## - Fill ww3_shel.inp file -##
sed -e "s/<yr1>/${YEAR_BEGIN_JOB}/g"  -e "s/<mo1>/${ms}/g"  -e "s/<dy1>/${ds}/g"  -e "s/<hr1>/00/g"  \
    -e "s/<yr2>/${YEAR_END_JOB}/g"  -e "s/<mo2>/${me}/g"  -e "s/<dy2>/${de}/g"  -e "s/<hr2>/24/g" \
    -e "s/<wav_int>/${wav_int}/g"  -e "s/<wav_rst>/$(( ${TOTAL_JOB_DUR} * 24 * 3600))/g" \
    -e "s/<wavdt>/${DT_WAV}/g" \
    ${WAV_NAM_DIR}/ww3_shel.inp.base.${RUNtype} > ./ww3_shel.inp


echo 'link ww3 input files and copy associated settings files'
lengthforc=${#forcww3[@]}

cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
while [[ `echo ${cur_M} | cut -b 1` -eq 0 ]]; do
    cur_M=`echo ${cur_M} | cut -b 2-`
done
mdy=$( valid_date ${MONTH_END_JOB} $(( ${DAY_END_JOB} + 1 )) ${YEAR_END_JOB} )
LOCAL_MTH_END=$( echo $mdy | cut -d " " -f 1 )

for k in `seq 0 $(( ${lengthforc} - 1))` ; do
    module load ${ncomod}
    if [[ ${JOB_DUR_MTH} -eq 1 || ${LOCAL_MTH_END} -eq ${cur_M} ]]; then # Case 1 month or less
        echo "Job is one month long or less ---> Using netcdf of the current month"
        cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
        cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
        [[ ! -f ${WAV_FILES_DIR}/${forcin[$k]}_Y${cur_Y}M${cur_M}.nc ]] && { echo "Missing ${WAV_FILES_DIR}/${forcin[$k]}_Y${cur_Y}M${cur_M}.nc to build ${forcww3[$k]}.nc file."; exit ;}
        ln -sf ${WAV_FILES_DIR}/${forcin[$k]}_Y${cur_Y}M${cur_M}.nc ./${forcww3[$k]}.nc
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
            [[ ! -f ${WAV_FILES_DIR}/${forcin[$k]}_Y${cur_Y}M${cur_M}.nc ]] && { echo "Missing ${WAV_FILES_DIR}/${forcin[$k]}_Y${cur_Y}M${cur_M}.nc to build ${forcww3[$k]}.nc file."; exit ;}
            if [[ $i == 0 ]]; then
                ncrcat ${WAV_FILES_DIR}/${forcin[$k]}_Y${cur_Y}M${cur_M}.nc ./${forcww3[$k]}.nc
            else
                ncrcat -O ${forcww3[$k]}.nc ${WAV_FILES_DIR}/${forcin[$k]}_Y${cur_Y}M${cur_M}.nc ./${forcww3[$k]}.nc
            fi
        done
    fi
    module unload ${ncomod}
done

if [ ! -z $bouncin ]; then
    echo "link ww3 boundary files"
    echo "ln -sf ${WAV_FILES_DIR}/$bouncin* ./"
    ${io_getfile} ${WAV_FILES_DIR}/$bouncin* ./

    echo "cp -f ${WAV_NAM_DIR}/ww3_bounc.inp ./"
    cpfile ${WAV_NAM_DIR}/ww3_bounc.inp ./
fi

cp ${WAV_FILES_DIR}/*.inp ./.
