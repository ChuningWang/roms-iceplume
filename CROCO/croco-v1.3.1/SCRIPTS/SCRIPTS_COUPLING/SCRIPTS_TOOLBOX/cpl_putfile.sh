
#-------------------------------------------------------------------------------
#                                                                      Restart
#-------------------------------------------------------------------------------
if [ ${USE_ATM} -eq 1 ]; then
    for dom in $wrfcpldom; do
        if [ $dom == "d01" ]; then
            ${io_putfile} atm.nc ${RESTDIR_OUT}/atm_${CEXPER}_${DATE_END_JOB}.nc
        else
            ${io_putfile} atm${dom}.nc ${RESTDIR_OUT}/atm${dom}_${CEXPER}_${DATE_END_JOB}.nc
        fi
    done
fi
if [ ${USE_OCE} -eq 1 ]; then
   for nn in $( seq 0 ${AGRIFZ} ); do
        if [ ${nn} -gt 0 ];    then
            agrif_ext=".${nn}"
        else
            agrif_ext=""
        fi
        ${io_putfile} oce.nc${agrif_ext} ${RESTDIR_OUT}/oce_${CEXPER}_${DATE_END_JOB}.nc${agrif_ext}
   done
fi
[ ${USE_WAV} -eq 1 ] && ${io_putfile} wav.nc ${RESTDIR_OUT}/wav_${CEXPER}_${DATE_END_JOB}.nc

[[ ${USE_ATM} -eq 1 && ${USE_OCE} -eq 1 ]] && cp *atmt*_to_ocn* ${RESTDIR_OUT}/. && cp *ocn*_to_atmt* ${RESTDIR_OUT}/. 
[[ ${USE_ATM} -eq 1 && ${USE_WAV} -eq 1 ]] && cp *atmt*_to_ww3t* ${RESTDIR_OUT}/. && cp *ww3t_to_atmt* ${RESTDIR_OUT}/.
[[ ${USE_OCE} -eq 1 && ${USE_WAV} -eq 1 ]] && cp *ocn*_to_ww3t* ${RESTDIR_OUT}/. && cp *ww3t_to_ocn* ${RESTDIR_OUT}/. 

if [ ${USE_TOY} -ge 1 ] ; then
    for k in `seq 0 $(( ${nbtoy} - 1 ))`; do
        printf "move ${toytype[$k]}.nc"
        ${io_putfile} ${toytype[$k]}.nc ${RESTDIR_OUT}/${toytype[$k]}_${CEXPER}_${DATE_END_JOB}.nc
    done
    [ ${USE_OCE} -eq 1 ] && cp *toy*_to_ocn* ${RESTDIR_OUT}/. && cp *ocn*_to_toy* ${RESTDIR_OUT}/. 
    [ ${USE_WAV} -eq 1 ] && cp *toy*_to_ww3t* ${RESTDIR_OUT}/. && cp *ww3t_to_toy* ${RESTDIR_OUT}/.
    [ ${USE_ATM} -eq 1 ] && cp *toy*_to_atmt* ${RESTDIR_OUT}/. && cp *atmt_to_toy* ${RESTDIR_OUT}/.
    if [ ${nbtoy} -gt 1 ]; then
        cp *toy*_to_toy*
    fi
fi

cpfile2 grids.nc ${RESTDIR_OUT}/. && cpfile2 masks.nc ${RESTDIR_OUT}/. && cpfile2 areas.nc ${RESTDIR_OUT}/.




