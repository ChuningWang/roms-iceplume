#-------------------------------------------------------------------------------
#                                                                      Restart
#-------------------------------------------------------------------------------

if [[ ${RESTART_FLAG} == "FALSE" ]]; then
    cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
    cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
    for nn in $( seq 0 ${AGRIFZ} ); do
       if [ ${nn} -gt 0 ]; then
           agrif_ext=".${nn}"
       else
           agrif_ext=""
       fi
           ${io_getfile} ${OCE_FILES_DIR}/croco_${ini_ext}_Y${cur_Y}M${cur_M}.nc${agrif_ext} croco_ini.nc${agrif_ext}
    done
else
    for nn in $( seq 0 ${AGRIFZ} ); do
       if [ ${nn} -gt 0 ]; then
           agrif_ext=".${nn}"
       else
           agrif_ext=""
       fi
       ${io_getfile} ${RESTDIR_IN}/croco_rst_${DATE_END_JOBm1}.nc${agrif_ext} croco_ini.nc${agrif_ext}
    done
fi
#
if [[ ${USE_ATM} == 1 ]]; then
    if [[ -f "${OCE_FILES_DIR}/coupling_masks*.nc" ]]; then
        cp ${OCE_FILES_DIR}/coupling_masks*.nc .
    elif [[ $( echo "$wrfcpldom" | wc -w ) >1 ]]; then
        module load ${ncomod}
        for nn in `seq 0 $AGRIFZ`; do
            if [ ${nn} -gt 0 ]; then
                agrif_ext=".${nn}"
            else
                agrif_ext=""
            fi
            [ -f coupling_masks${nn}.nc ] && rm -rf coupling_masks${nn}.nc
            echo "Creating coupling mask for CROCO"
            ncks -v lon_rho,lat_rho croco_grd.nc${agrif_ext} coupling_masks${nn}.nc
            ncap2 -A -v -s "cplmask1=mask_rho*0+1" croco_grd.nc${agrif_ext} coupling_masks${nn}.nc
            maxdom=$( echo "$wrfcpldom" | wc -w )
            cnt=0
            Alevel=$(( $nn +1 ))
            while [ ${Alevel} -le ${maxdom} ]; do
                cnt=$(( $cnt +1 ))
                ncap2 -A -v -s "cplmask${cnt}=mask_rho*0+1" croco_grd.nc${agrif_ext} coupling_masks${nn}.nc
                if [[ ${RESTART_FLAG} == "FALSE" ]]; then
                    cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
                    cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
                    filefrom="${ATM_FILES_DIR}/wrfinput_d$( printf "%02d" ${Alevel})_${cur_Y}_${cur_M}"
                else
                    cur_Y=$( echo $DATE_END_JOBm1 | cut -c 1-4 )
                    cur_M=$( echo $DATE_END_JOBm1 | cut -c 5-6 )
                    filefrom="${RESTDIR_IN}/wrfrst_d0$( printf "%02d" ${Alevel})_${cur_Y}-${cur_M}"
                fi
                ncap2 -O -v -s 'latmin=XLAT.min();latmax=XLAT.max();lonmin=XLONG.min();lonmax=XLONG.max()' ${filefrom}* tmp.nc
                lonmin=$( ncdump -v lonmin tmp.nc  | grep "lonmin =" | cut -d ' ' -f 4)
                latmin=$( ncdump -v latmin tmp.nc  | grep "latmin =" | cut -d ' ' -f 4)
                lonmax=$( ncdump -v lonmax tmp.nc  | grep "lonmax =" | cut -d ' ' -f 4)
                latmax=$( ncdump -v latmax tmp.nc  | grep "latmax =" | cut -d ' ' -f 4)
                \rm tmp.nc
                for index in `seq 1 $cnt` ; do
                    if [[ $(( $index + $nn )) -lt $Alevel ]]; then
                        ncap2 -O -s "var_tmp=cplmask${index}; where( lat_rho > $latmin && lon_rho > $lonmin && lat_rho < $latmax && lon_rho < $lonmax ) var_tmp=0; cplmask${index}=var_tmp" coupling_masks${nn}.nc coupling_masks${nn}.nc
                    elif [[ $(( $index + $nn )) -eq $Alevel ]]; then
                        ncap2 -O -s "var_tmp=cplmask${index}; where( lat_rho < $latmin || lon_rho < $lonmin || lat_rho > $latmax || lon_rho > $lonmax ) var_tmp=0; cplmask${index}=var_tmp" coupling_masks${nn}.nc coupling_masks${nn}.nc
                    fi
                done
                ncks -O -v var_tmp -x coupling_masks${nn}.nc coupling_masks${nn}.nc
                Alevel=$(( ${Alevel} +1 ))
            done
        done
        module unload ${ncomod}
    fi
fi

