

#===============================================================================
#  Print
#===============================================================================
#printf "\n Period of Experiment:\n"
#printf  "%12s%10d%22s%14d%17s\n" "BEGINING:" ${DATE_BEGIN_EXP} "(in ${CALTYPE} calendar)" ${JDAY_BEGIN_EXP} "(in Julian days)"
#printf  "%12s%10d%22s\n\n" "END:" ${DATE_END_EXP} "(in ${CALTYPE} calendar)"

#printf  "%16s%10d%20s%10d\n"     "YEAR_BEGIN_EXP:  " ${YEAR_BEGIN_EXP}   "YEAR_END_EXP:  " ${YEAR_END_EXP}
#printf  "%16s%10d%20s%10d\n"     "MONTH_BEGIN_EXP: " ${MONTH_BEGIN_EXP}  "MONTH_END_EXP: " ${MONTH_END_EXP}
#printf  "%16s%10d%20s%10d\n\n\n" "DAY_BEGIN_EXP:   " ${DAY_BEGIN_EXP}    "DAY_END_EXP:   " ${DAY_END_EXP}

printf "\n Period of Job:\n"
printf  "%12s%10d%22s%14d%17s\n" "BEGINING:" ${DATE_BEGIN_JOB} "(in ${CALTYPE} calendar)" ${JDAY_BEGIN_JOB} "(in Julian days)"
printf          "%44s%14d%17s\n" "Job duration" ${TOTAL_JOB_DUR} "(in Julian days)"
printf  "%12s%10d%22s%14d%17s\n\n\n" "END:" ${DATE_END_JOB} "(in ${CALTYPE} calendar)" ${JDAY_END_JOB} "(in Julian days)"

printf  "%16s%10d%20s%10d%20s%10d\n"     "YEAR_BEGIN_JOB:  " ${YEAR_BEGIN_JOB}   "YEAR_END_JOB:  " ${YEAR_END_JOB}    "YEAR_BEGIN_JOBp1:  " ${YEAR_BEGIN_JOBp1}
printf  "%16s%10d%20s%10d%20s%10d\n"     "MONTH_BEGIN_JOB: " ${MONTH_BEGIN_JOB}  "MONTH_END_JOB: " ${MONTH_END_JOB}   "MONTH_BEGIN_JOBp1: " ${MONTH_BEGIN_JOBp1}
printf  "%16s%10d%20s%10d%20s%10d\n\n\n" "DAY_BEGIN_JOB:   " ${DAY_BEGIN_JOB}    "DAY_END_JOB:   " ${DAY_END_JOB}    "DAY_BEGIN_JOBp1:   " ${DAY_BEGIN_JOBp1} 

printf "%-15s : %-80s\n" "CHOME" "${CHOME}"
printf "%-15s : %-80s\n" "CWORK" "${CWORK}"
printf "%-15s : %-80s\n" "SCRIPTDIR" "${SCRIPTDIR}"
printf "%-15s : %-80s\n\n" "JOBDIR" "${JOBDIR_ROOT}"

printf "%-15s : %-80s\n" "ATM_INPUTDIR" "${ATM_FILES_DIR}"
printf "%-15s : %-80s\n" "OCE_INPUTDIR" "${OCE_FILES_DIR}"
printf "%-15s : %-80s\n" "WAV_INPUTDIR" "${WAV_FILES_DIR}"
printf "%-15s : %-80s\n" "CPL_INPUTDIR" "${CPL_NAM_DIR}"
printf "%-15s : %-80s\n" "OUTPUTDIR_ROOT" "${OUTPUTDIR_ROOT}"
printf "%-15s : %-80s\n" "RESTDIR_ROOT" "${RESTDIR_ROOT}"

printf "%-15s : %-80s\n\n" "EXEDIR" "${EXEDIR_ROOT}"

printf "%-15s : %-80s\n\n" "LD_LIBRARY_PATH" "${LD_LIBRARY_PATH}"


