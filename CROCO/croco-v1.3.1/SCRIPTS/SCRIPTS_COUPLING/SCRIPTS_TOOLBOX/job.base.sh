
############################### from job.base.sh ###############################

        DATE1=`date "+%Y%m%d-%H:%M:%S"`
        OLD="old_${DATE1}"
        printf "\n date : ${DATE1}\n"


#===============================================================================
#  Step 1 commands:  get_file step  
#===============================================================================
if [ ${LOADL_STEP_NAME} == "get_file" ] || [ ${LOADL_STEP_NAME} == "XXX" ]; then
        
        #  Careful if modification in these vars, modif in the 3 steps of the job job.base.sh
        export JOBDIR="${JOBDIR_ROOT}/${ROOT_NAME_2}"
        export EXEDIR="${EXEDIR_ROOT}/${ROOT_NAME_2}"

# EXEDIR: Execution directory
        [ -d ${EXEDIR} ] && mv ${EXEDIR} ${EXEDIR}_${OLD}
        mkdir -p ${EXEDIR} && mkdir ${EXEDIR}/ls_l

# JOBDIR: job directory
        [ -d ${JOBDIR} ] && mv ${JOBDIR} ${JOBDIR}_${OLD}
        mkdir -p ${JOBDIR}

# some printings
. ${SCRIPTDIR}/common_printing.sh

cd ${EXEDIR}

# Copy job in EXEDIR
cp ${JOBDIR_ROOT}/${jobname} ./

#        printf "\n ************* SCRIPT files *****************\n\n"
#        for file in ${SCRIPTDIR}/*.sh ; do cpfile ${file} . ; done; chmod 755 *.sh

        printf "\n ************* EXECUTABLE files *****************\n\n"
            if [ ${USE_OCE}  -eq 1 ]; then
                if [ ${ONLINE_COMP} -eq 1 ]; then
                    . ${SCRIPTDIR}/oce_compile.sh 
                else
                    cpfile ${OCE_EXE_DIR}/croco.${RUNtype} crocox
                fi
                . ${SCRIPTDIR}/getversion.sh ${OCE}
            fi
#	    [ ${USE_OCE}  -eq 1 ] && cpfile ${OCE_EXE_DIR}/croco.${RUNtype} crocox
	    [ ${USE_ATM}  -eq 1 ] && { cpfile ${ATM_EXE_DIR}/wrf.exe wrfexe ; . ${SCRIPTDIR}/getversion.sh ${ATM} ; }
	    [ ${USE_WAV}  -eq 1 ] && { cp ${WAV_EXE_DIR}/ww3_* . && mv ww3_shel wwatch ; . ${SCRIPTDIR}/getversion.sh ${WAV} ;}
	    [ ${USE_XIOS} -ge 1 ] && cpfile ${XIOS_EXE_DIR}/xios_server.exe .
            # If toy is used #
            if [ ${USE_TOY}  -ge 1 ]; then
	        if [ ${nbtoy} -eq 1 ]; then
		    cpfile ${TOY_EXE_DIR}/toy_model toyexe
		else
		    for k in `seq 0 $(( ${nbtoy} - 1))`; do
			cpfile ${TOY_EXE_DIR}/toy_${toytype[$k]} toy${toytype[$k]}
		    done
		fi
	    fi
	    ##

        printf "\n ************* PARAMETER files *****************\n\n"
# if xios
	[[ ${ONLINE_XML} == "TRUE" ]] && { . ${SCRIPTDIR}/xios_process.sh ;}
        [ ${USE_XIOS} -ge 1 ] && {  cp ${XIOS_NAM_DIR}/*.xml . ; echo ""; }

# ocean/atmosphere input files (configuration, forcing, obc, levitus/restart...)

        DATE_END_JOBm1=$( makedate $MONTH_BEGIN_JOB $(( $DAY_BEGIN_JOB - 1 )) $YEAR_BEGIN_JOB )
        RESTDIR_IN=${RESTDIR_ROOT}/${DATE_END_JOBm1}

        [ ${USE_OCE} -eq 1 ] && printf "\n ************* get OCEAN CONFIGURATION, OBC, BLK... files *****************\n\n" 
        [ ${USE_OCE} -eq 1 ] && { . ${SCRIPTDIR}/oce_getfile.sh ; }

        [ ${USE_OCE} -eq 1 ] && printf "\n ************* get OCEAN RESTART files *****************\n\n" |tee ls_l/oce_getrst.txt
        [ ${USE_OCE} -eq 1 ] && printf "    see listing in ${EXEDIR}/ls_l/oce_getrst.txt \n"
        [ ${USE_OCE} -eq 1 ] && { . ${SCRIPTDIR}/oce_getrst.sh ; } >> ls_l/oce_getrst.txt

        [ ${USE_ATM} -eq 1 ] && printf "\n ************* get ATMOSPHERE CONFIGURATION LOWINPUT, BDY... files *****************\n\n" 
        [ ${USE_ATM} -eq 1 ] && { . ${SCRIPTDIR}/atm_getfile.sh ; }

	[ ${USE_ATM} -eq 1 ] && printf "\n ************* get ATMOSPHERE RESTART files *****************\n\n" |tee ls_l/atm_getrst.txt
	[ ${USE_ATM} -eq 1 ] && printf "    see listing in ${EXEDIR}/ls_l/atm_getrst.txt \n"
	[ ${USE_ATM} -eq 1 ] && { . ${SCRIPTDIR}/atm_getrst.sh ; } >> ls_l/atm_getrst.txt

        [ ${USE_WAV} -eq 1 ] && printf "\n ************* get WAVE CONFIGURATION files *****************\n\n"
        [ ${USE_WAV} -eq 1 ] && { . ${SCRIPTDIR}/wav_getfile.sh ; }

        [ ${USE_WAV} -eq 1 ] && printf "\n ************* get WAVE RESTART files *****************\n\n" |tee ls_l/wav_getrst.txt
        [ ${USE_WAV} -eq 1 ] && printf "    see listing in ${EXEDIR}/ls_l/wav_getrst.txt \n"
        [ ${USE_WAV} -eq 1 ] && { . ${SCRIPTDIR}/wav_getrst.sh ; } >> ls_l/wav_getrst.txt

	[ ${USE_CPL} -ge 1 ] && printf "\n ************* get OA3MCT RESTART files *****************\n\n"
	[ ${USE_CPL} -ge 1 ] && { . ${SCRIPTDIR}/cpl_getrst.sh ; }
        [ ${USE_CPL} -ge 1 ] && { . ${SCRIPTDIR}/cpl_getfile.sh ; }

        [ ${USE_TOY} -ge 1 ] && printf "\n ************* get TOY CONFIGURATION files *****************\n\n"
        [ ${USE_TOY} -ge 1 ] && { . ${SCRIPTDIR}/toy_getfile.sh ; }

# make the namelists from the namelist.base files
        printf "\n ************* make namelist files from namelist base files *****************\n\n"
	[ ${USE_OCE} -eq 1 ] && { . ${SCRIPTDIR}/oce_nam.sh ; }
	[ ${USE_ATM} -eq 1 ] && { . ${SCRIPTDIR}/atm_nam.sh ; }
        #[ ${USE_WAV} -eq 1 ] && echo "No namelist for WAV"
        [ ${USE_TOY} -ge 1 ] && { . ${SCRIPTDIR}/toy_nam.sh ; }
	[ ${USE_CPL} -ge 1 ] && { . ${SCRIPTDIR}/cpl_nam.sh ; }
        printf "\n date : `date "+%Y%m%d-%H:%M:%S"`\n"

fi # Step1
#===============================================================================
#  Step 2 commands:  run_model step (parallel)
#===============================================================================
if [ ${LOADL_STEP_NAME} == "run_model" ] || [ ${LOADL_STEP_NAME} == "XXX" ] ; then
        
        export JOBDIR="${JOBDIR_ROOT}/${ROOT_NAME_2}"
        export EXEDIR="${EXEDIR_ROOT}/${ROOT_NAME_2}"

cd ${EXEDIR} 

#       ls -l > ls_l/ls_pre_exe.txt
        printf "\n see ls -l in ${EXEDIR}/ls_l/ls_pre_exe.txt\n"
        printf "\n date : `date "+%Y%m%d-%H:%M:%S"`\n"
        printf "\n---------------  start   ---------------\n"
#      
        printf "\n ***************** make app.conf file for Multiple Program Multiple Data *****************\n\n"
        [ -f ${SCRIPTDIR}/MACHINE/${MACHINE}/launch_${MACHINE}.sh ] && { . ${SCRIPTDIR}/MACHINE/${MACHINE}/launch_${MACHINE}.sh ; } || { printf "\n Please create a file launch_${MACHINE}.sh \n" ; exit 0 ; }

	echo "launch run: $MPI_LAUNCH ${MPI_ext} app.conf " #${run_cmd} #$app.conf

##### RUN ######
time $MPI_LAUNCH ${MPI_ext} app.conf >& out_run.txt 
[ "$?" -gt "0" ] && { printf "ERROR during the RUN.\n Please check le log.files in ${EXEDIR} (out_run.txt, croco.log,...)"; exit ; }
################
#        fi
#	time ccc_mprun -f app.conf > out_run.txt 2>&1
#	time ccc_mprun -f app.conf 1>out_run.txt 2>err_run.txt

        printf "\n\n PWD: `pwd`\n"

        [ -f time.step ] && printf "\n  time.step: `cat time.step` \n"
        printf "\n---------------   end    ---------------\n"
        printf "\n date : `date "+%Y%m%d-%H:%M:%S"`\n"
        printf "\n see ls -l in ${EXEDIR}/ls_l/ls_post_exe.txt\n"
        ls -l > ls_l/ls_post_exe.txt
 
fi # Step2
 

#===============================================================================
#  Step 3 commands:  put_file step 
#===============================================================================
if [ ${LOADL_STEP_NAME} == "put_file" ] || [ "${LOADL_STEP_NAME}" == "XXX" ] ; then
     
        export JOBDIR="${JOBDIR_ROOT}/${ROOT_NAME_2}"
        export EXEDIR="${EXEDIR_ROOT}/${ROOT_NAME_2}"
        export OUTPUTDIR="${OUTPUTDIR_ROOT}/${ROOT_NAME_2}"
        export RESTDIR_OUT="${RESTDIR_ROOT}/${ROOT_NAME_3}"
        
# OUTPUTDIR: output directory
        ${MACHINE_STOCKAGE} ls ${OUTPUTDIR}  >  /dev/null  2>&1
        [ "$?" -eq "0" ] && ${MACHINE_STOCKAGE} mv ${OUTPUTDIR} ${OUTPUTDIR}_${OLD}
        ${MACHINE_STOCKAGE} mkdir -p ${OUTPUTDIR}

# RESTDIR_OUT: restart directory
        ${MACHINE_STOCKAGE} ls ${RESTDIR_OUT}  >  /dev/null  2>&1
        [ "$?" -eq "0" ] && ${MACHINE_STOCKAGE} mv ${RESTDIR_OUT} ${RESTDIR_OUT}_${OLD}
        ${MACHINE_STOCKAGE} mkdir -p ${RESTDIR_OUT}

cd ${EXEDIR} 
        printf "\n date : `date "+%Y%m%d-%H:%M:%S"`\n"
        [ ${USE_OCE} -eq 1 ] && printf "\n ************* put OCEAN OUTPUT/RESTART files *****************\n\n" |tee ls_l/oce_putfile.txt
        [ ${USE_OCE} -eq 1 ] && printf "    see listing in ${EXEDIR}/ls_l/oce_putfile.txt \n" 
        [ ${USE_OCE} -eq 1 ] && { . ${SCRIPTDIR}/oce_putfile.sh ; } >> ls_l/oce_putfile.txt

        [ ${USE_ATM} -eq 1 ] && printf "\n ************* put ATMOSPHERE OUTPUT/RESTART files *****************\n\n" |tee ls_l/atm_putfile.txt
        [ ${USE_ATM} -eq 1 ] && printf "    see listing in ${EXEDIR}/ls_l/atm_putfile.txt \n"
        [ ${USE_ATM} -eq 1 ] && { . ${SCRIPTDIR}/atm_putfile.sh ; } >> ls_l/atm_putfile.txt
      
        [ ${USE_WAV} -eq 1 ] && printf "\n ************* put WAVE OUTPUT/RESTART files *****************\n\n" |tee ls_l/wav_putfile.txt
        [ ${USE_WAV} -eq 1 ] && printf "    see listing in ${EXEDIR}/ls_l/wav_putfile.txt \n"
        [ ${USE_WAV} -eq 1 ] && { . ${SCRIPTDIR}/wav_putfile.sh ; } >> ls_l/wav_putfile.txt

        [ ${USE_CPL} -ge 1 ] && printf "\n ************* put OA3MCT OUTPUT/RESTART files *****************\n\n" 
        [ ${USE_CPL} -ge 1 ] && { . ${SCRIPTDIR}/cpl_putfile.sh ; }

        chmod -R ${permission} ${OUTPUTDIR} ${RESTDIR_OUT}
#-------------------------------------------------------------------------------
#  save output control ascii files in jobs directory
#-------------------------------------------------------------------------------

        printf "\n ************* save ascii job files in ${JOBDIR} *****************\n\n" 
# if ocean
        FILES_OCE="layout.dat ocean.output* namelist out_run.txt *time.step solver.stat* croco.log"        
        [ ${USE_OCE} -eq 1 ] && {  for file in ${FILES_OCE}; do cpfile2 ${file} ${JOBDIR}; done; echo ""; }
# if atmosphere
        FILES_ATM="namelist.input rsl.error.0000 rsl.out.0000"
        [ ${USE_ATM} -eq 1 ] && {  for file in ${FILES_ATM}; do cpfile2 ${file} ${JOBDIR}; done; echo ""; }
# if wave
        FILES_WAV="log.ww3 strt.out ounf.out prnc.*.out grid.out"
        [ ${USE_WAV} -eq 1 ] && {  for file in ${FILES_WAV}; do cpfile2 ${file} ${JOBDIR}; done; echo ""; }
# if coupler
        FILES_CPL="nout.000000 wrfexe.timers* crocox.timers* wwatch.timers*"
        [ ${USE_CPL} -ge 1 ] && {  for file in ${FILES_CPL}; do cpfile2 ${file} ${JOBDIR}; done; echo ""; }
# if agrif
        FILES_AGRIFZ="AGRIF_FixedGrids.in"
        [ ${AGRIFZ} -eq 1 ] && {  for file in ${FILES_AGRIFZ}; do cpfile2 ${file} ${JOBDIR}; done; echo ""; }
# job
        FILES_JOB="${jobname}"
        if [ ${COMPUTER} == "IRENE" ]; then
            FILES_JOB="${FILES_JOB} ${ROOT_NAME_1}*.o ${ROOT_NAME_1}*.e"
        elif [ ${COMPUTER} == "JEANZAY" ]; then
            FILES_JOB="${FILES_JOB} ${ROOT_NAME_1}.out"
        else
            FILES_JOB="${FILES_JOB} ${ROOT_NAME_1}.jobid_*.txt ${ROOT_NAME_1}.o*"
        fi
        
        cd ${JOBDIR_ROOT}; for file in ${FILES_JOB}; do mvfile2 ${file} ${JOBDIR}; done;  cd -; echo "";
#
#  mv output files in output
        printf "\n date : `date "+%Y%m%d-%H:%M:%S"`\n"

#-------------------------------------------------------------------------------
#  NEXT job!
#-------------------------------------------------------------------------------
        printf "\n *************  run the next job  *****************\n\n"

if [ "${MODE_TEST}" == "" ] ; then      #  en production
        sed -e "s/YEAR_BEGIN_JOB=.*/YEAR_BEGIN_JOB=${YEAR_BEGIN_JOBp1}/" \
            -e "s/MONTH_BEGIN_JOB=.*/MONTH_BEGIN_JOB=${MONTH_BEGIN_JOBp1}/" \
            -e "s/DAY_BEGIN_JOB=.*/DAY_BEGIN_JOB=${DAY_BEGIN_JOBp1}/" \
            -e "s/NBJOB=.*/NBJOB=$(( ${NBJOB} - 1 ))/" \
            -e "s/RESTART_FLAG=.*/RESTART_FLAG=\"TRUE\"/" \
            -e "s/CHAINED_JOB=.*/CHAINED_JOB=\"FALSE\"/" \
                ${SCRIPTDIR}/../myjob.sh > tata.sh
        mvfile2 ${SCRIPTDIR}/../myjob.sh ${JOBDIR}
        mvfile tata.sh ${SCRIPTDIR}/../myjob.sh
        chmod 755   ${SCRIPTDIR}/../myjob.sh

        cpfile ${SCRIPTDIR}/../mynamelist.sh ${JOBDIR}

	if [ ${DATE_END_JOB} -ne ${DATE_END_EXP} ]
	then
        cd ${SCRIPTDIR}/..
        chmod 755 submitjob.sh
        ./submitjob.sh
    else
	    printf "\n     run finished at date : ${DATE_END_EXP} \n\n\n\n"
	    exit 0
	fi

else # en test
        printf "\n\n\n\n  MODE_TEST=${MODE_TEST}  Test mode and non production => No job chaining.\n\n\n\n"
fi

fi # Step3

