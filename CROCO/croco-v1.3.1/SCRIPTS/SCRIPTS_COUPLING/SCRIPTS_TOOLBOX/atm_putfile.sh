#-------------------------------------------------------------------------------
#                                                                      Average
#-------------------------------------------------------------------------------


if [ ${USE_XIOS_ATM} -eq 1 ] ; then
    for file in ${ATM_XIOS_NAME}; do
	mv ${file}*.nc* ${OUTPUTDIR}/
    done
else
    module load $ncomod
    if [[ ${ATM_CASE} == "MOVING_NEST" ]]; then
        totnbdom=$(( ${NB_dom} + ${num_mv_nest} ))
    else
        totnbdom=${NB_dom}
    fi
    for dom in `seq 1 ${totnbdom}`; do
        ncrcat -O -F -d Time,1,-2 wrfout_d0${dom}_${YEAR_BEGIN_JOB}-* wrfout_d0${dom}_${DATE_BEGIN_JOB}_${DATE_END_JOB}.nc
        mv wrfout_d0${dom}_${DATE_BEGIN_JOB}_${DATE_END_JOB}.nc ${OUTPUTDIR}/.
        \rm wrfout_d0${dom}_${YEAR_BEGIN_JOB}-*
        mv wrfxtrm_d0${dom}_${YEAR_BEGIN_JOB}-* ${OUTPUTDIR}/wrfxtrm_d0${dom}_${DATE_BEGIN_JOB}_${DATE_END_JOB}.nc
   done
   module unload $ncomod
fi

#-------------------------------------------------------------------------------
#                                                                      Restart
#-------------------------------------------------------------------------------
year=$( printf "%04d"   ${YEAR_BEGIN_JOBp1} )
month=$( printf "%02d"   ${MONTH_BEGIN_JOBp1} )
day=$( printf "%02d"   ${DAY_BEGIN_JOBp1} )

date_rst="${year}-${month}-${day}"

for i in  wrfrst_d0?_${date_rst}_* 
   do
     ${io_putfile} $i ${RESTDIR_OUT}
done

