################################################################################
############################ END USER CHANGE ###################################
################################################################################

# Avoid error when a model is not used
[[ -z ${USE_ATM+x} ]] && export USE_ATM=0
[[ -z ${USE_OCE+x} ]] && export USE_OCE=0
[[ -z ${USE_WAV+x} ]] && export USE_WAV=0
[[ -z ${USE_XIOS_ATM+x} ]] && export USE_XIOS_ATM=0
[[ -z ${USE_XIOS_OCE+x} ]] && export USE_XIOS_OCE=0
[[ -z ${USE_TOYATM+x} ]] && export USE_TOYATM=0
[[ -z ${USE_TOYOCE+x} ]] && export USE_TOYOCE=0
[[ -z ${USE_TOYWAV+x} ]] && export USE_TOYWAV=0
[[ -z ${istoy+x} ]] && export istoy=""
[[ -z ${ONLINE_XML+x} ]] && export ONLINE_XML="FALSE"


if [ ${USE_ATM} == 0 ]; then
    export DT_ATM=1
fi

if [ ${USE_OCE} == 0 ]; then
    export DT_OCE=1
fi

if [ ${USE_WAV} == 0 ]; then
    export DT_WAV=1
fi

# KEY for XIOS and TOY #
export USE_XIOS=$(( ${USE_XIOS_ATM} + ${USE_XIOS_OCE} ))
export USE_TOY=$(( ${USE_TOYATM} + ${USE_TOYOCE} + ${USE_TOYWAV} ))
[ ${USE_TOY} -ge 1 ] && export USE_CPL=1 || export USE_CPL=$(( ${USE_ATM} * ${USE_OCE} + ${USE_ATM} * ${USE_WAV} + ${USE_OCE} * ${USE_WAV} ))

### TOY ###
[ -z ${toytype+x} ] && export toytype=()

export nbtoy=${#toytype[@]}
export model_to_toy=()
export toynamelist=()
export DT_TOY=()
if [[ ${USE_TOY} -ge 1 ]]; then
[[ ${USE_TOY} != $nbtoy ]] && { echo "There is an incoherence between the number of USE_TOY_* activated and the number of arguement in toytype. Make Sure they are coherent." ; exit;}
for k in `seq 0 $(( ${nbtoy} - 1))` ; do
    [ ${toytype[$k]} == "oce" ] && model_to_toy+=("croco")
    [ ${toytype[$k]} == "atm" ] && model_to_toy+=("wrf")
    [ ${toytype[$k]} == "wav" ] && model_to_toy+=("ww3")
    if [ ${nbtoy} -eq 1 ]; then
        toynamelist+=("TOYNAMELIST.nam.${toytype[$k]}.${RUNtype}")
    elif [ ${nbtoy} -eq 2 ]; then
        toynamelist+=("TWOTOYNAMELIST.nam.${toytype[$k]}.${RUNtype}")
    else
        printf "\n No namelist available for three namelist toy\n" ; exit 0
    fi
    if [ ${toytype[$k]} == "oce" ]; then
        targ=$( ncdump -v time ${toyfile[$k]} | grep "time = " | sed -n '2p' | cut -d ' ' -f 4-5 )
        tsp=$(( $( echo "${targ}" | cut -d ',' -f 2) - $( echo "${targ}" | cut -d ',' -f 1) ))
        DT_TOY+=("${tsp}")
    elif [ ${toytype[$k]} == "atm" ]; then
        targ=$( ncks -v XTIME ${toyfile[$k]} | grep "XTIME ="| cut -d '=' -f 2 | cut -d ',' -f 1-2)
        tsp=$(( $( echo "${targ}" | cut -d ',' -f 2)*60 - $( echo "${targ}" | cut -d ',' -f 1)*60 ))
        DT_TOY+=("${tsp}")
    elif [ ${toytype[$k]} == "wav" ]; then
        targ=$( ncdump -v time ${toyfile[$k]} | grep "time =" | sed -n '2p' | cut -d ' ' -f 4-5)
        arg1=$( echo "${targ}" | cut -d ',' -f 2 )
        arg2=$( echo "${targ}" | cut -d ',' -f 1 )
        tsp=$( echo "(${arg1}  - ${arg2})*86400" | bc | cut -d '.' -f 1)
        DT_TOY+=("${tsp}")
    fi
done
fi
######

