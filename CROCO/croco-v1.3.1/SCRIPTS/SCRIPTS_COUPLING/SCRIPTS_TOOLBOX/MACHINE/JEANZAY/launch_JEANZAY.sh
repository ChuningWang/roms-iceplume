
############ CREATE app.conf for JEANZAY ###########
mystartproc=0

if [ ${USE_ATM} -eq 1 ]; then
    myendproc=$(( ${NP_ATM} - 1 ))
    mod_Str=$mystartproc"-"$myendproc
    echo "$mod_Str ./wrfexe" >> app.conf
    if [ ${USE_XIOS_ATM} -eq 1 ]; then
        mystartproc=$(( ${myendproc} + 1 ))
        myendproc=$(( ${mystartproc} + ${NP_XIOS_ATM} - 1 ))
        mod_Str=$mystartproc"-"$myendproc
        echo "${mod_Str} ./xios_server.exe" >> app.conf
    fi

fi
#
if [ ${USE_OCE} -eq 1 ]; then
    if [[ ${MPI_NOLAND} == "TRUE" ]]; then
        OCE_PROCS=${MY_NODES}
    else
        OCE_PROCS=$(( ${NP_OCEX}*${NP_OCEY} ))
    fi 
    [ ${USE_ATM} -eq 1 ] && { mystartproc=$(( ${myendproc} + 1 )) ; myendproc=$(( ${mystartproc} + ${OCE_PROCS} - 1 )); } || { myendproc=$(( ${OCE_PROCS} - 1 )) ; }
    mod_Str=$mystartproc"-"$myendproc
    echo "$mod_Str ./crocox" >> app.conf
    if [ ${USE_XIOS_OCE} -eq 1 ]; then
        mystartproc=$(( ${myendproc} + 1 ))
        myendproc=$(( ${mystartproc} + ${NP_XIOS_OCE} - 1 ))
        mod_Str=$mystartproc"-"$myendproc
        echo "${mod_Str} ./xios_server.exe" >> app.conf
    fi
fi
#
if [ ${USE_WAV} -eq 1 ]; then
    if [ ${USE_ATM} -eq 1 ] || [ ${USE_OCE} -eq 1 ]; then
        mystartproc=$(( ${myendproc} + 1 ))
        myendproc=$(( ${mystartproc} + ${NP_WAV} - 1 ))
    else
        myendproc=$(( ${NP_WAV} - 1 ))
    fi
    mod_Str=$mystartproc"-"$myendproc
    echo "${mod_Str} ./wwatch" >> app.conf
fi

if [ ${USE_TOY} -ge 1 ]; then
    if [ ${USE_ATM} -eq 1 ] || [ ${USE_OCE} -eq 1 ] || [ ${USE_WAV} -eq 1 ]; then
        mystartproc=$(( ${myendproc} + 1 ))
        myendproc=$(( ${mystartproc} + ${NP_TOY} - 1 ))
    fi
    if [ ${nbtoy} -eq 1 ]; then
        mod_Str=$mystartproc"-"$myendproc
        echo "${mod_Str} ./toyexe" >> app.conf
     else
        for k in `seq 0 $(( ${nbtoy} - 1 ))`; do
            mod_Str=$mystartproc"-"$myendproc
            echo "${mod_Str} ./toy${toytype[$k]}" >> app.conf
            mystartproc=$(( ${myendproc} + 1 ))
            myendproc=$(( ${mystartproc} + ${NP_TOY} - 1 ))
        done
     fi
fi



