#-------------------------------------------------------------------------------
#                                                                   mozaic files
#-------------------------------------------------------------------------------
if [[ ${WEIGHT_FLAG} == 1 ]]; then
    weight_files="${weight_o2a} ${weight_a2o}"
    for file in ${weight_files}; do
        ${io_getfile}  ${CPL_FILES_DIR}/${file} .
    done
fi
