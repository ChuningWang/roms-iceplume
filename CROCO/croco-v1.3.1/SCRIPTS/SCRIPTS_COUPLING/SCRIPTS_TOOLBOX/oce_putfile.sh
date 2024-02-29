
#-------------------------------------------------------------------------------
#                                                                      Average
#-------------------------------------------------------------------------------

if [ ${USE_XIOS_OCE} -eq 1 ]; then 
    for file in ${OCE_XIOS_NAME}; do
        mv ${file}* ${OUTPUTDIR}/
    done
else
    module load $ncomod
    for nn in $( seq 0 ${AGRIFZ} ); do
        if [ ${nn} -gt 0 ];    then
            agrif_ext=".${nn}"
        else
            agrif_ext=""
        fi
        for ff in croco_his.nc${agrif_ext} croco_avg.nc${agrif_ext} ; do
            [ ${ff} == "croco_his.nc${agrif_ext}" ] && { name=croco_his ; ncrcat -O -F -d time,1,-2 croco_his.nc${agrif_ext} croco_his.nc${agrif_ext} ;} 
            [ ${ff} == "croco_avg.nc${agrif_ext}" ] && name=croco_avg
            mvfile2 ${ff} ${OUTPUTDIR}/${name}_${DATE_BEGIN_JOB}_${DATE_END_JOB}.nc${agrif_ext}
        done
    done
   module unload $ncomod
fi

#-------------------------------------------------------------------------------
#                                                                      Restart
#-------------------------------------------------------------------------------
for nn in $( seq 0 ${AGRIFZ} ); do
    if [ ${nn} -gt 0 ];    then
        agrif_ext=".${nn}"
    else
        agrif_ext=""
    fi
    for ff in croco_rst.nc ; do
#    fff=${CEXPER}_${DATE_END_JOB}_rst${ff/rst.????/}
        fff=croco_rst_${DATE_END_JOB}.nc${agrif_ext}
        ${io_putfile} $ff${agrif_ext} ${RESTDIR_OUT}/$fff
    done
done
