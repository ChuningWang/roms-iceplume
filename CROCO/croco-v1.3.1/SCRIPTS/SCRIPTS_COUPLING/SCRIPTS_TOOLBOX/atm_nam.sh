#!/bin/bash
#set -ue
#set -vx
#
. ${SCRIPTDIR}/caltools.sh
##
##======================================================================
##----------------------------------------------------------------------
##    II. modify namelist
##----------------------------------------------------------------------
##======================================================================
##
#

if [[ ${RESTART_FLAG} == "FALSE" ]]; then
  rst="false"
  file="wrfinput_d01"
else
  rst="true"
  file="wrfrst_d01_*"
fi

# look for starting hour in input file

module load ${ncomod}
if [[ ${RESTART_FLAG} == "FALSE" ]]; then
    YY=$( printf "%04d" ${YEAR_BEGIN_JOB} )
    MM=$( printf "%02d" ${MONTH_BEGIN_JOB} )
    DD=$( printf "%02d" ${DAY_BEGIN_JOB} )
    file="wrfinput_d01"
    fulldate=$( echo "$( ncdump -v Times ${file} )" | grep -m 1 "${YY}-${MM}-${DD}" | cut -d '=' -f 2)
    hh=$( echo ${fulldate} | cut -d '_' -f 2 | cut -d ':' -f 1 )
else
    YY=$( printf "%04d" ${YEAR_BEGIN_JOB} )
    MM=$( printf "%02d" ${MONTH_BEGIN_JOB} )
    DD=$( printf "%02d" ${DAY_BEGIN_JOB} )
    file="wrfrst_d01_*"
    fulldate=$( echo "$( ncdump -v Times ${file} )" | grep -m 2 "${YY}-${MM}-${DD}" | cut -d '{' -f 2 | cut -d '"' -f 2 )
    hh=$( echo ${fulldate} | cut -d '_' -f 2 | cut -d ':' -f 1 )
fi
module unload ${ncomod}
#

sed -e "s/<yr1>/${YEAR_BEGIN_JOB}/g"   -e "s/<yr2>/${YEAR_END_JOB}/g"  \
    -e "s/<mo1>/${MONTH_BEGIN_JOB}/g"   -e "s/<mo2>/${MONTH_END_JOB}/g"  \
    -e "s/<dy1>/${DAY_BEGIN_JOB}/g"   -e "s/<dy2>/${DAY_END_JOB}/g"  \
    -e "s/<hr1>/${hh}/g"   -e "s/<hr2>/24/g"  \
    -e "s/<rst>/$rst/g"              -e "s/<rst_int_h>/$(( ${TOTAL_JOB_DUR} * 24 ))/g"            \
    -e "s/<his_int_h>/${atm_his_h}/g"         -e "s/<his_nb_out>/${atm_his_frames}/g"    \
    -e "s/<xtrm_int_m>/${atm_diag_int_m}/g"   -e "s/<xtrm_nb_out>/${atm_diag_frames}/g"  \
    -e "s/<nproc_x>/${atm_nprocX}/g"            -e "s/<nproc_y>/${atm_nprocY}/g"             \
    -e "s/<niotaskpg>/${atm_niotaskpg}/g"       -e "s/<niogp>/${atm_niogp}/g"                \
    -e "s/time_step                           =.*/time_step                           =${DT_ATM} ,/g" \
    -e "s/<interval_s>/${interval_seconds}/g" \
    -e "s/<nbmetlev>/${nbmetlevel}/g"     -e "s/<nbmetsoil>/${nbmetsoil}/g"  \
    $ATM_NAM_DIR/${atmnamelist} > ./namelist.input

for domnb in `seq 1 $NB_dom` ; do
    dom="d$( printf "%02d" ${domnb})"
    if [[ ${RESTART_FLAG} == "FALSE" ]]; then
        file="wrfinput_${dom}"
    else
        file="wrfrst_${dom}*"
    fi
#
    searchf=("<xdim_${dom}>" "<ydim_${dom}>" "<nbvertlev>" "<dx_${dom}>" "<dy_${dom}>" "<i_str_${dom}>" "<j_str_${dom}>" "<coef_${dom}>" "<ptop>")
#
    dimx=$( ncdump -h  $file  | grep "west_east_stag =" | cut -d ' ' -f 3)
    dimy=$( ncdump -h  $file  | grep "south_north_stag =" | cut -d ' ' -f 3)
    dimz=$( ncdump -h  $file  | grep "bottom_top_stag =" | cut -d ' ' -f 3)
    dx=$( ncdump -h  $file | grep "DX =" | cut -d ' ' -f 3 | cut -d '.' -f 1)
    dy=$( ncdump -h  $file  | grep "DY =" | cut -d ' ' -f 3 | cut -d '.' -f 1)
    istrt=$( ncdump -h  $file  | grep "I_PARENT_START =" | cut -d ' ' -f 3)
    jstrt=$( ncdump -h  $file  | grep "J_PARENT_START =" | cut -d ' ' -f 3)
    coef=$( ncdump -h  $file  | grep "PARENT_GRID_RATIO =" | cut -d ' ' -f 3)
    ptop=$( ncdump -v P_TOP $file | grep "P_TOP =" | cut -d '=' -f 2 | cut -d ' ' -f 2  ) 

    sed -e "s/${searchf[0]}/${dimx}/g" \
        -e "s/${searchf[1]}/${dimy}/g" \
        -e "s/${searchf[2]}/${dimz}/g" \
        -e "s/${searchf[3]}/${dx}/g" \
        -e "s/${searchf[4]}/${dy}/g" \
        -e "s/${searchf[5]}/${istrt}/g" \
        -e "s/${searchf[6]}/${jstrt}/g" \
        -e "s/${searchf[7]}/${coef}/g" \
        -e "s/${searchf[8]}/${ptop}/g" \
        ./namelist.input > ./namelist.input.tmp
        mv namelist.input.tmp namelist.input
    chmod 755 namelist.input
done
#
if [[ ${switch_fdda} == 1 ]]; then
    nbdom=$( echo "${nudgedom}" | wc -w)
    [[ ${nbdom} >  $( echo "${nudge_coef}" | wc -w) ]] && { echo "Missing values in nudge_coef for nest, we stop..."; exit ;}
    [[ ${nbdom} >  $( echo "${nudge_interval_m}" | wc -w) ]] && { echo "Missing values in nudge_interval_m for nest, we stop..."; exit ;}
    [[ ${nbdom} >  $( echo "${nudge_end_h}" | wc -w) ]] && { echo "Missing values in nudge_end_h for nest, we stop..."; exit ;}

    for lvl in `seq 1 ${nbdom}`; do
        dom=$(printf "d%02d" ${lvl})
        sed -e "s/<nudge_${dom}>/$( echo "${nudgedom}" | cut -d " " -f $(( ${lvl})) )/g"                 \
            -e "s/<nudge_coef_${dom}>/$( echo "${nudge_coef}" | cut -d " " -f $(( ${lvl})) )/g"       \
            -e "s/<nudge_end_h_${dom}>/$( echo "${nudge_end_h}" | cut -d " " -f $(( ${lvl})) )/g"      \
            -e "s/<nudge_int_m_${dom}>/$( echo "${nudge_interval_m}" | cut -d " " -f $(( ${lvl})) )/g" \
            namelist.input > namelist.tmp
            mv namelist.tmp namelist.input
   done

    sed -e "s/<nudge_d.*>/0/g" \
        -e "s/<nudge_coef_d.*>/0/g" \
        -e "s/<nudge_end_h_d.*>/0/g" \
        -e "s/<nudge_int_m_d.*>/0/g" \
        namelist.input > namelist.tmp
    mv namelist.tmp namelist.input
    chmod 755 namelist.input
else
    sed -e "s/<nudge_d.*>/0/g" \
        -e "s/<nudge_coef_d.*>/0/g" \
        -e "s/<nudge_end_h_d.*>/0/g" \
        -e "s/<nudge_int_m_d.*>/0/g" \
        namelist.input> namelist.tmp
    mv namelist.tmp namelist.input
    chmod 755 namelist.input
fi
#
if [[ ${nestfeedback} == "FALSE" ]]; then
    sed -e "s/feedback                            = 1,/feedback                            = 0,/g" \
    ./namelist.input > ./namelist.input.tmp
    mv namelist.input.tmp namelist.input
fi
###### handle cpl dom ######
maxatmdom=$( echo "$wrfcpldom" | wc -w )
if [[ ${maxatmdom} > 0 ]] ; then
    max_cpldom=$( echo "${wrfcpldom}" | cut -d ' ' -f ${maxatmdom} | cut -d '0' -f 2)
    [[ ${max_cpldom} > 2 ]] && { echo "In the current state WRF can not couple more than 2 domains, we stop..."; exit; }
else
    max_cpldom=0
fi
#
sed -e "s/<max_cpldom>/${max_cpldom}/g" \
    ./namelist.input > ./namelist.input.tmp
mv namelist.input.tmp namelist.input

if [[ $maxatmdom == 1 && $AGRIFZ > 1 ]];then 
    numextmod=$(( $AGRIFZ + 1 ))
else
    [[ $maxatmdom > $(( $AGRIFZ +1 )) ]] && { numextmod=$maxatmdom ;} || { numextmod=$(( $AGRIFZ + 1 )) ;}
fi
sed -e "s/num_ext_model_couple_dom            = 1,/num_ext_model_couple_dom            = $numextmod,/g" \
./namelist.input > ./namelist.input.tmp
mv namelist.input.tmp namelist.input
####

sed -e "s/<isftcflx>/${isftcflx}/g" \
    ./namelist.input > ./namelist.input.tmp
mv namelist.input.tmp namelist.input

if [ $USE_WAV -eq 1 ] || [ $USE_TOYWAV -eq 1 ]; then
    [[ ${isftcflx} -ne 5 ]] && { echo "ERROR... Atm parameter (in namelist.input) isftcflx need to be 5 when coupling with WAV..." ; exit ; }
fi
#
if [[ ${ATM_CASE} == "MOVING_NEST" ]]; then
    [[ ${num_mv_nest} >  $( echo "${ref_coef}" | wc -w) ]] && { echo "Missing values in ref_coef  for nest, we stop..."; exit ;}
    [[ ${num_mv_nest} >  $( echo "${ew_size}" | wc -w) ]] && { echo "Missing values in ew_size for nest, we stop..."; exit ;}
    [[ ${num_mv_nest} >  $( echo "${ns_size}" | wc -w) ]] && { echo "Missing values in ns_size for nest, we stop..."; exit ;}
    [[ ${num_mv_nest} >  $( echo "${i_prt_strt}" | wc -w) ]] && { echo "Missing values in i_prt_str for nest, we stop..."; exit ;}
    [[ ${num_mv_nest} >  $( echo "${j_prt_strt}" | wc -w) ]] && { echo "Missing values in j_prt_str for nest, we stop..."; exit ;}

    for nest_nb in `seq 1 ${num_mv_nest}`; do
        dom=$( printf "d%02d" $(( ${nest_nb} +1 )) )
        rcoef=$( echo "${ref_coef}" | cut -d " " -f $(( ${nest_nb})) )
        xdim=$( echo "${ew_size}" | cut -d " " -f $(( ${nest_nb})) )    
        ydim=$( echo "${ns_size}" | cut -d " " -f $(( ${nest_nb})) )
        istrt=$( echo "${i_prt_strt}" | cut -d " " -f $(( ${nest_nb})) )
        jstrt=$( echo "${j_prt_strt}" | cut -d " " -f $(( ${nest_nb})) )
        echo "Computing nest dx/dy from its parent's grid"
        dx=`echo "${dx}/${rcoef}" | bc`
        dy=`echo "${dy}/${rcoef}" | bc`
    
        searchf=("<xdim_${dom}>" "<ydim_${dom}>" "<dx_${dom}>" "<dy_${dom}>" "<i_str_${dom}>" "<j_str_${dom}>" "<coef_${dom}>")
        sed -e "s/${searchf[0]}/${xdim}/g" \
            -e "s/${searchf[1]}/${ydim}/g" \
            -e "s/${searchf[2]}/${dx}/g" \
            -e "s/${searchf[3]}/${dy}/g" \
            -e "s/${searchf[4]}/${istrt}/g" \
            -e "s/${searchf[5]}/${jstrt}/g" \
            -e "s/${searchf[6]}/${rcoef}/g" \
            ./namelist.input > ./namelist.input.tmp
        mv namelist.input.tmp namelist.input
        chmod 755 namelist.input
    done
# add vortex_interval,max_vortex_speed,corral_dist,track_level,time_to_move
    linetoadd=("numtiles" "vortex_interval" "max_vortex_speed" "corral_dist" "track_level" "time_to_move")
    nbline=${#linetoadd[@]}

    for lineadd in `seq 0 $(( ${nbline} - 2))` ; do
        line1=${linetoadd[${lineadd}]}
        linep1=${linetoadd[$(( ${lineadd} +1 ))]}
        [[ ${linep1} == "vortex_interval" ]] && newline="vortex_interval                     = ${vortex_interval},${vortex_interval},${vortex_interval},"
        [[ ${linep1} == "max_vortex_speed" ]] && newline="max_vortex_speed                    = ${max_vortex_speed},${max_vortex_speed},${max_vortex_speed},"
        [[ ${linep1} == "corral_dist" ]] && newline="corral_dist                         = ${corral_dist},${corral_dist},${corral_dist},"
        [[ ${linep1} == "track_level" ]] && newline="track_level                         = ${track_level},"
        [[ ${linep1} == "time_to_move" ]] && newline="time_to_move                        = ${time_to_move},${time_to_move},${time_to_move},"

        sed -e "/${line1}.*/ a ${newline}" \
            ./namelist.input > ./namelist.input.tmp
        mv namelist.input.tmp namelist.input
        chmod 755 namelist.input
    done 

fi

if [[ ${ATM_CASE} == "MOVING_NEST" ]]; then
    sed -e "s|input_from_file                     = .true.,.true.,.true.,|input_from_file                     = .true.,.false.,.false.,|g" \
        -e "s|auxinput4_interval                  = <sst_int_m>,<sst_int_m>,<sst_int_m>,|auxinput4_interval                  = ${auxinput4_interval},|"\
        -e "s/<max_domains>/$(( ${NB_dom}+ ${num_mv_nest} ))/g" \
        ./namelist.input > ./namelist.input.tmp
     mv namelist.input.tmp namelist.input
     chmod 755 namelist.input
else
    sed -e "s|<sst_int_m>|${auxinput4_interval}|g"\
        -e "s/<max_domains>/${NB_dom}/g" \
        ./namelist.input > ./namelist.input.tmp
     mv namelist.input.tmp namelist.input
     chmod 755 namelist.input
fi



sed -e "s/<xdim_d.*>/1/g"  -e "s/<ydim_d.*>/1/g"\
    -e "s/<dx_d.*>/1/g"  -e "s/<dy_d.*>/1/g" \
    -e "s/<i_str_d.*>/1/g"  -e "s/<j_str_d.*>/1/g" \
    -e "s/<coef_d.*>/1/g"  -e "s/<coef_d.*>/1/g" \
./namelist.input > ./namelist.input.tmp
mv namelist.input.tmp namelist.input
chmod 755 namelist.input

