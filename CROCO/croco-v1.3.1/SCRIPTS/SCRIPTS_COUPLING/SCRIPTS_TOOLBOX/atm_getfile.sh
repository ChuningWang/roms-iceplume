#-------------------------------------------------------------------------------
#                                                          Configuration files
#-------------------------------------------------------------------------------
# link data files necessary for running wrf in a dedicated directory $wrf/data
if [ ! -d ${ATM_EXE_DIR}/../data ] ; then
 mkdir ${ATM_EXE_DIR}/../data
 ln -s ${ATM_EXE_DIR}/../run/* ${ATM_EXE_DIR}/../data/.
 # remove executables that could exist and namelist file
 rm -f ${ATM_EXE_DIR}/../data/*.exe
 rm -f ${ATM_EXE_DIR}/../data/namelist.input*
fi
echo 'link wrf data files'
echo "ln -sf ${ATM_EXE_DIR}/../data/* ."
ln -sf ${ATM_EXE_DIR}/../data/* .

#-------------------------------------------------------------------------------
#                                                          BDY
#-------------------------------------------------------------------------------
module load ${ncomod}


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
    [[ ! -f ${ATM_FILES_DIR}/wrfbdy_d01_${cur_Y}_${cur_M} ]] && { echo "Missing ${ATM_FILES_DIR}/wrfbdy_d01_${cur_Y}_${cur_M} to build atm bdy file."; exit ;}
    ${io_getfile} ${ATM_FILES_DIR}/wrfbdy_d01_${cur_Y}_${cur_M} wrfbdy_d01
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
         [[ ! -f ${ATM_FILES_DIR}/wrfbdy_d01_${cur_Y}_${cur_M} ]] && { echo "Missing ${ATM_FILES_DIR}/wrfbdy_d01_${cur_Y}_${cur_M} to build atm bdy file."; exit ;}
         if [[ $i == 0 ]]; then
             ncrcat ${ATM_FILES_DIR}/wrfbdy_d01_${cur_Y}_${cur_M} ./wrfbdy_d01
         else
             ncrcat -O wrfbdy_d01 ${ATM_FILES_DIR}/wrfbdy_d01_${cur_Y}_${cur_M} ./wrfbdy_d01
         fi
    done
fi
module unload ${ncomod}

    
#-------------------------------------------------------------------------------
#                                            Forcing fields (interannual case)
#-------------------------------------------------------------------------------
filelist='wrflowinp_d01'
if [ $NB_dom -ge 2 ] ; then
    filelist="$filelist wrflowinp_d02"
    if [ $NB_dom -eq 3 ] ; then
        filelist="$filelist wrflowinp_d03"
    fi
fi

module load ${ncomod}

for file in ${filelist}; do
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
        [[ ! -f ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ]] && { echo "Missing ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} to build atm bdy file."; exit ;}
        ${io_getfile} ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ${file}
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
            [[ ! -f ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ]] && { echo "Missing ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} to build atm bdy file."; exit ;}
            if [[ $i == 0 ]]; then
                ncrcat -O -F -d Time,1,-2 ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ${file}
            elif [[ $i == ${nbloop} ]]; then
                ncrcat  -O ${file} ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ${file}
            else
                ncrcat -O -F -d Time,1,-2 ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ${file}.tmp
                ncrcat -O ${file} ${file}.tmp ${file}
                rm -rf ${file}.tmp
            fi
        done
    fi
done
module unload ${ncomod}
#-------------------------------------------------------------------------------
#                                                        Spectral forcing (fdda)
#-------------------------------------------------------------------------------

if [[ ${switch_fdda} == 1 ]]; then
    filelist=''
    if [[ "$( echo ${nudgedom} | cut -d ' ' -f 1)" -ge 1 ]]; then
        filelist="$filelist wrffdda_d01" ;
        if [[ "$( echo ${nudgedom} | wc -w )" -ge 2 && "$( echo ${nudgedom} | cut -d ' ' -f 2)" -ge 1 ]]; then
            filelist="$filelist wrffdda_d02"
            if [[ "$( echo ${nudgedom} | wc -w )" -ge 3 && "$( echo ${nudgedom} | cut -d ' ' -f 3)" -ge 1 ]]; then
                filelist="$filelist wrffdda_d03"
            fi
       fi
    fi
    module load ${ncomod}

    for file in ${filelist}; do
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
            [[ ! -f ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ]] && { echo "Missing ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} to build atm bdy file."; exit ;}
            ${io_getfile} ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ${file}
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
                [[ ! -f ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ]] && { echo "Missing ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} to build atm fdda file."; exit ;}
                if [[ $i == 0 ]]; then
                    ncrcat ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ${file}
                else
                    ncrcat -O ${file} ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ${file}
                fi
            done
        fi
    done
    module unload ${ncomod}
fi
