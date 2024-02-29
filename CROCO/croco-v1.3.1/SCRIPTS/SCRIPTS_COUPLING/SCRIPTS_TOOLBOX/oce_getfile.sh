#!/bin/bash
#-------------------------------------------------------------------------------
#                                                          Configuration files
#-------------------------------------------------------------------------------

# Get grid file
for nn in $( seq 0 ${AGRIFZ} ); do
    if [ ${nn} -gt 0 ];    then
        agrif_ext=".${nn}"
    else
        agrif_ext=""
    fi
    ${io_getfile} ${OCE_FILES_DIR}/croco_grd.nc${agrif_ext}                croco_grd.nc${agrif_ext}
done

# Get boundary foring
. ${SCRIPTDIR}/oce_getbdy.sh

# Get surface forcing if needed
[ ${surfrc_flag} == "TRUE" ] && . ${SCRIPTDIR}/oce_getfrc.sh

# Get tide forcing if any
for nn in $( seq 0 ${AGRIFZ} ); do
    if [ ${nn} -gt 0 ];    then
        agrif_ext=".${nn}"
    else
        agrif_ext=""
    fi
    [ ${tide_flag} == "TRUE" ] && ${io_getfile} ${OCE_FILES_DIR}/croco_frc.nc${agrif_ext} ./ 
    [ ${river_flag} == "TRUE" ] && ${io_getfile} ${OCE_FILES_DIR}/croco_runoff.nc${agrif_ext} ./
done

          
	       














