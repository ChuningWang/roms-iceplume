#!/bin/bash

if [ ${MACHINE} == "DATARMOR" ] ; then
    ${QSUB} -h -N ${ROOT_NAME_1} ${jobname}
    sub_ext=".pbs"
    launchcmd="-W depend=afterany:$( echo $( qselect -N ${ROOT_NAME_1}) | cut -c 1-7 )"
elif [ ${MACHINE} == "IRENE" ]; then
    ${QSUB} ${jobname}
    sub_ext=".sh"
    prejobname=${ROOT_NAME_1}
    echo "$( ccc_mstat -f -u $USER )" >tmp.text
    cnt=$(( $(wc -l tmp.text | awk '{ print $1 }') + 1 )) # check if there are already jobs running
    \rm -f tmp.text  
elif [ ${MACHINE} == "JEANZAY" ]; then
    ${QSUB} -H ${jobname}
    sub_ext=".sh"
    prejobname=${ROOT_NAME_1}
    echo "$( squeue --format="%.18i %.100j" -u $USER )" >tmp.text
    cnt=$(( $(wc -l tmp.text | awk '{ print $1 }') + 1 ))
    for linenb in `seq 1 $cnt`; do 
        line=$(sed -n "${linenb}p" tmp.text  )
	[[ $line == *"${prejobname}"* ]] && { prejobid=$( echo ${line} | cut -c 1-7 ); break; }       
    done
    firstjobid=${prejobid}
else
    printf "\n\n Chained job for your machine is not set up yet, we stop...\n\n" && exit 1

fi


cd ${SCRIPTDIR}/

newsdate=$( makedate $MONTH_BEGIN_JOBp1 $DAY_BEGIN_JOBp1 $YEAR_BEGIN_JOBp1 )
newedate=${newsdate}
months=$MONTH_BEGIN_JOBp1
days=$DAY_BEGIN_JOBp1
years=$YEAR_BEGIN_JOBp1
jobleft=$NBJOB
while [ ${newedate} -lt ${DATE_END_EXP} ] ; do
    #
    mdy=$( valid_date $(( $months + $JOB_DUR_MTH )) $(( $days + $JOB_DUR_DAY - 1 )) $years )
    monthe=$( echo $mdy | cut -d " " -f 1 )
    daye=$( echo $mdy | cut -d " " -f 2 )
    yeare=$( echo $mdy | cut -d " " -f 3 )
    newedate=$( makedate $monthe $daye $yeare )   
    jobleft=$(( $jobleft - 1))
    #
    future_date="${CEXPER}_${newsdate}_${newedate}${MODE_TEST}"   
    future_job="job_${future_date}${sub_ext}"
    cd ${JOBDIR_ROOT} 
    #
    sed -e "s/YEAR_BEGIN_JOB=${YEAR_BEGIN_JOB}/YEAR_BEGIN_JOB=${years}/" \
        -e "s/MONTH_BEGIN_JOB=${MONTH_BEGIN_JOB}/MONTH_BEGIN_JOB=${months}/" \
        -e "s/DAY_BEGIN_JOB=${DAY_BEGIN_JOB}/DAY_BEGIN_JOB=${days}/" \
        -e "s/export NBJOB=.*/export NBJOB=${jobleft}/" \
        -e "s/export RESTART_FLAG=.*/export RESTART_FLAG=\"TRUE\"/" \
        -e "s/export CHAINED_JOB=.*/export CHAINED_JOB=\"FALSE\"/" \
        ${jobname} > ${future_job}

    chmod 755 ${future_job}
    #
    if [ ${MACHINE} == "DATARMOR" ] ; then 
        newjobname="${CEXPER}_${newsdate}_${newedate}${MODE_TEST}"
        ${QSUB} -N ${newjobname} ${launchcmd} ${future_job} 
        launchcmd="-W depend=afterany:$( echo $( qselect -N ${newjobname}) | cut -c 1-7 )"
        cd ${SCRIPTDIR}/
#
    elif [ ${MACHINE} == "IRENE" ] ; then
        # 
        newjobname="${CEXPER}_${newsdate}_${newedate}${MODE_TEST}"
        sed -e "s/#MSUB -r .*/#MSUB -r ${newjobname}/" \
            -e "s/#MSUB -o .*/#MSUB -o ${newjobname}.jobid_%I.o/" \
            -e "s/#MSUB -e .*/#MSUB -e ${newjobname}.jobid_%I.e/" \
        ${future_job} > ${future_job}.tmp
        mv ${future_job}.tmp ${future_job}
        chmod 755 ${future_job}
        #
        echo "$( ccc_mstat -f -u $USER )" >tmp.text
        for linenb in `seq 3 $cnt`; do
            line=$(sed -n "${linenb}p" tmp.text  )
            [[ $line == *"${prejobname}"* ]] && { prejobid=$( echo ${line} | cut -c 1-7 ); break; }
        done
        \rm -f tmp.text
        #
        cnt=$(( ${cnt} + 1 ))
        ${QSUB} -a ${prejobid} ${future_job} 
        prejobname="${newjobname}"
        cd ${SCRIPTDIR}/
#
    elif [ ${MACHINE} == "JEANZAY" ] ; then
	newjobname="${CEXPER}_${newsdate}_${newedate}${MODE_TEST}"
	sed -e "s/#SBATCH --job-name=.*/#SBATCH --job-name=${newjobname}/" \
	    -e "s/#SBATCH --output=.*/#SBATCH --output=${newjobname}.out/" \
	    -e "s/#SBATCH --error=.*/#SBATCH --error=${newjobname}.out/" \
	${future_job} > ${future_job}.tmp
	mv ${future_job}.tmp ${future_job}
	chmod 755 ${future_job}
	echo "$( squeue --format="%.18i %.100j" -u $USER )" >tmp.text
        for linenb in `seq 1 $cnt`; do
            line=$(sed -n "${linenb}p" tmp.text  )
            [[ $line == *"${prejobname}"* ]] && { prejobid=$( echo ${line} | cut -c 1-7 ); break; }
        done
        \rm -f tmp.text
        #
        cnt=$(( ${cnt} + 1 ))
	launchcmd="--dependency=afterok:${prejobid}"
        ${QSUB} ${launchcmd} ${future_job}
	prejobname="${newjobname}"
        cd ${SCRIPTDIR}

    fi
#    
    mdy=$( valid_date $(( $monthe )) $(( $daye + 1 )) $yeare )
    months=$( echo $mdy | cut -d " " -f 1 )
    days=$( echo $mdy | cut -d " " -f 2 )
    years=$( echo $mdy | cut -d " " -f 3 )
    newsdate=$( makedate $months $days $years )
    
done
#
cd ${JOBDIR_ROOT}

[ ${MACHINE} == "DATARMOR" ] && { jobid=$( echo $( qselect -N ${ROOT_NAME_1}) | cut -c 1-7 ) ; qrls ${jobid} ; }
[ ${MACHINE} == "JEANZAY" ] && { scontrol release ${firstjobid} ; }
