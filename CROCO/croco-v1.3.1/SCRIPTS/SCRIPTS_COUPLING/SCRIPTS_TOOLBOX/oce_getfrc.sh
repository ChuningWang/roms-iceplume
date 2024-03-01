#!/bin/bash
#set -x

#
if [ ${interponline} -eq 1 ]; then 
    if [[ ${frc_ext} == *'AROME'* || ${frc_ext} == *'ARPEGE'* ]]; then
        vnames="${frc_ext}"
        ${io_getfile} ${OCE_FILES_ONLINEDIR}/${frc_ext} .
    else
        if [ ${frc_ext} == "ERA_ECMWF" ]; then
            vnames='T2M U10M V10M Q STRD SSR TP'
        else
            vnames='Temperature_height_above_ground Specific_humidity Precipitation_rate Downward_Short-Wave_Rad_Flux_surface Upward_Short-Wave_Rad_Flux_surface Downward_Long-Wave_Rad_Flux Upward_Long-Wave_Rad_Flux_surface U-component_of_wind V-component_of_wind'
        fi
#    
        printf "Creating link to data for the job duration\n"
#          
        echo "Checking if Previous month is needed"
        cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
        cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
        if [[ ${RESTART_FLAG} == "FALSE" ]]; then
            filefrom="${OCE_FILES_DIR}/croco_${ini_ext}_Y${cur_Y}M${cur_M}.nc"
        else
            filefrom="${RESTDIR_IN}/croco_rst_${DATE_END_JOBm1}.nc"
        fi
        # scrum_time of ini file
        tstartinsec=$( echo $( ncdump -v scrum_time ${filefrom} | grep 'scrum_time =' | cut -d '=' -f 2| cut -d ' ' -f 2 ))
        tstartinsec=`echo "scale=2; ${tstartinsec} + ${DT_OCE}*0.5" | bc ` # =0.5*dt like in croco 
        # Find first time value in forcing file
        fieldname=$( echo "$vnames" | awk '{print $1}' )
        ncdump -v time "${OCE_FILES_ONLINEDIR}/${fieldname}_Y${cur_Y}M${cur_M}.nc" | grep -n 'time =' > tmp$$
        ns=$( ncdump -v time ${OCE_FILES_ONLINEDIR}/${fieldname}_Y${cur_Y}M${cur_M}.nc | grep -c 'time =' )
        tstartfrc=`echo "scale=2; $( sed -n -e "${ns} p" tmp$$ | cut -d '=' -f 2 | cut -d ',' -f 1 ) * 86400" | bc`
        rm -rf tmp$$
        [[ $( echo "${tstartinsec}<=${tstartfrc}" | bc )>0 ]] && { echo "Previous month is needed!"; loopstrt=-1 ;} || { loopstrt=0 ;}      
#
        for i in `seq ${loopstrt} $(( ${JOB_DUR_MTH} ))`; do
            [ ${i} -eq -1 ] && printf "Adding link to the previous month (for temporal interpolation)\n"
            [ ${i} -eq ${JOB_DUR_MTH} ] && printf "Adding link to the following month (for temporal interpolation)\n"
 
            mdy=$( valid_date $(( $MONTH_BEGIN_JOB + $i )) $DAY_BEGIN_JOB $YEAR_BEGIN_JOB )
            cur_Y=$( printf "%04d\n"  $( echo $mdy | cut -d " " -f 3) )
            cur_M=$( printf "%02d\n"  $( echo $mdy | cut -d " " -f 1) )

            for varname in ${vnames} ; do
                [[ ! -f "${OCE_FILES_ONLINEDIR}/${varname}_Y${cur_Y}M${cur_M}.nc" ]] && { echo "File ${varname}_Y${cur_Y}M${cur_M}.nc is missing for online interpolation, we stop..." ; exit ;}
                ${io_getfile} ${OCE_FILES_ONLINEDIR}/${varname}_Y${cur_Y}M${cur_M}.nc ./
            done
        done

    # Check if next month is need when job duration is smaller than a month
        cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
        while [[ `echo ${cur_M} | cut -b 1` -eq 0 ]]; do
            cur_M=`echo ${cur_M} | cut -b 2-`
        done
        mdy=$( valid_date ${MONTH_END_JOB} $(( ${DAY_END_JOB} +1 )) ${YEAR_END_JOB} )
        LOCAL_MTH_END=$( echo $mdy | cut -d " " -f 1 )

        if [[ ${JOB_DUR_MTH} -eq 0 && ${LOCAL_MTH_END} -ne ${cur_M} ]]; then
            mdy=$( valid_date $(( ${MONTH_BEGIN_JOB} + 1 )) 1 ${YEAR_BEGIN_JOB} )
            cur_Y=$( printf "%04d\n"  $( echo $mdy | cut -d " " -f 3) )
            cur_M=$( printf "%02d\n"  $( echo $mdy | cut -d " " -f 1) )
            for varname in ${vnames} ; do
                ${io_getfile} ${OCE_FILES_ONLINEDIR}/${varname}_Y${cur_Y}M${cur_M}.nc ./
            done
        fi
    fi
#
else
    if [[ ${frc_ext} == *'frc'* ]]; then
        extend="frc"
        timevar="frc_time"
    else
        extend="blk"
        timevar="bulk_time"
    fi

    for nn in $( seq 0 ${AGRIFZ} ); do
        if [ ${nn} -gt 0 ]; then
            agrif_ext=".${nn}"
        else
            agrif_ext=""
        fi

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
            [[ ! -f ${OCE_FILES_DIR}/croco_${frc_ext}_Y${cur_Y}M${cur_M}.nc ]] && { echo "Missing ${OCE_FILES_DIR}/croco_${frc_ext}_Y${cur_Y}M${cur_M}.nc to build oce frc file."; exit ;}
	    ${io_getfile} ${OCE_FILES_DIR}/croco_${frc_ext}_Y${cur_Y}M${cur_M}.nc croco_${extend}.nc${agrif_ext}
        else
            if [[ ${JOB_DUR_MTH} -eq 0 && ${LOCAL_MTH_END} -ne ${cur_M} ]]; then
                echo "Job is less than a month BUT overlaps on next month ---> Concat netcdf of current and following month"
                nbloop=1
            else
                echo "Job is longer than one month ---> Concat netcdf of needed month"
                echo "Warning: Default overlap value used is 1 (before and after)"
               nbloop=$(( ${JOB_DUR_MTH}-1 ))
            fi
            for i in `seq 0 ${nbloop}`; do
                mdy=$( valid_date $(( $MONTH_BEGIN_JOB + $i )) 1 $YEAR_BEGIN_JOB )
                cur_Y=$( printf "%04d\n"  $( echo $mdy | cut -d " " -f 3) )
                cur_M=$( printf "%02d\n"  $( echo $mdy | cut -d " " -f 1) )

                if [[ $i == 0 ]]; then
                    tstart=1; tend=-2
                elif [[ $i == $nbloop ]]; then
                    tstart=2; tend=""
                else
                    tstart=2; tend=-2
                fi
                
                if [[ ${extend} == "blk" ]]; then
                    if [[ ${i} == 0 ]]; then
                        ncks -O -F --mk_rec_dmn bulk_time -d bulk_time,${tstart},${tend} ${OCE_FILES_DIR}/croco_${frc_ext}_Y${cur_Y}M${cur_M}.nc${agrif_ext} croco_${extend}.nc${agrif_ext}
                    else
                        ncrcat -F -d bulk_time,${tstart},${tend} croco_${extend}.nc${agrif_ext} ${OCE_FILES_DIR}/croco_${frc_ext}_Y${cur_Y}M${cur_M}.nc${agrif_ext} croco_${extend}.nc${agrif_ext}
                    fi
                elif  [[ ${extend} == "frc" ]]; then
                    varlist='sms shf swf srf sst sss wwv'
                    for var in ${varlist}; do
                        [[ ${var} == "sms" ]] && extract="sustr svstr"
                        [[ ${var} == "shf" ]] && extract="shflux"
                        [[ ${var} == "swf" ]] && extract="swflux"
                        [[ ${var} == "srf" ]] && extract="swrad"
                        [[ ${var} == "sst" ]] && extract="SST dQdSST"
                        [[ ${var} == "sss" ]] && extract="SSS"
                        [[ ${var} == "wwv" ]] && extract="Awave Dwave Pwave"
                        if [[ $i > 0 ]]; then
                            ncks --mk_rec_dmn "${var}_time" -F -O -d "${var}_time",${tstart},${tend} -v "${extract}" ${OCE_FILES_DIR}/croco_${frc_ext}_Y${cur_Y}M${cur_M}.nc tmp_${var}.nc
                            ncks -O --mk_rec_dmn "${var}_time" croco_${extend}.nc${agrif_ext} croco_${extend}.nc${agrif_ext}
                            ncrcat -A croco_${extend}.nc${agrif_ext} tmp_${var}.nc croco_${extend}.nc${agrif_ext}
                            ncks -O --fix_rec_dmn "${var}_time" croco_${extend}.nc${agrif_ext} croco_${extend}.nc${agrif_ext}
                        else
                            ncks -F -O -d "${var}_time",${tstart},${tend} -v "${extract}" ${OCE_FILES_DIR}/croco_${frc_ext}_Y${cur_Y}M${cur_M}.nc tmp_${var}.nc
                            ncks -A -v "${extract}" tmp_${var}.nc croco_${extend}.nc${agrif_ext}
                        fi
                        rm -rf tmp_${var}.nc
                    done
                else
                    echo "Did not understand your atmospheric forcing type, blk or frc are available ( need to be specified in input file name)"
                    exit
                fi
            done
        fi
    done
fi
