#-------------------------------------------------------------------------------
#                                                                      Restart
#-------------------------------------------------------------------------------
if [[ ${RESTART_FLAG} == "FALSE" ]]
then

 filelist='wrfinput_d01' 
 if [ $NB_dom -ge 2 ] ; then
  filelist="$filelist wrfinput_d02"
  if [ $NB_dom -eq 3 ] ; then
   filelist="$filelist wrfinput_d03"
  fi
 fi
 cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
 cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
 for file in $filelist
  do
   echo "ln -sf ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ./$file"
   ln -sf ${ATM_FILES_DIR}/${file}_${cur_Y}_${cur_M} ./$file
  done
  
  if [[ ${onlinecplmask} == "TRUE" ]]; then
      echo "EDITING CPLMASK ONLINE"  
      for dom in $wrfcpldom ; do
          if [[ ${dom} == "d01" ]]; then
              echo "set CPLMASK to 1 in coupled domain $dom"
              echo "ncap2 -O -s \"CPLMASK(:,0,:,:)=(LANDMASK+LAKEMASK-1)*(-1)\" ./wrfinput_$dom ./wrfinput_$dom"
              module load $ncomod
              ncap2 -O -s "CPLMASK(:,0,:,:)=(LANDMASK+LAKEMASK-1)*(-1)" ./wrfinput_$dom ./wrfinput_$dom
              if [[ $(echo ${wrfcpldom} | wc -w) == 1 && $AGRIFZ > 0 ]]; then
                  ncpdq -O -d num_ext_model_couple_dom_stag,0 -v CPLMASK -a num_ext_model_couple_dom_stag,Time wrfinput_$dom tmp.nc
                  for nn in `seq 1 $AGRIFZ`; do
                      [[ ${nn} -gt 0 ]] && { agrif_ext=".${nn}" ;} || { agrif_ext="" ;}
                      ncap2 -O -v -s 'latmin=lat_rho.min();latmax=lat_rho.max();lonmin=lon_rho.min();lonmax=lon_rho.max()' ${OCE_FILES_DIR}/croco_grd.nc${agrif_ext} tmp.nc.1
                      lonmin=$( ncdump -v lonmin tmp.nc.1  | grep "lonmin =" | cut -d ' ' -f 4)
                      latmin=$( ncdump -v latmin tmp.nc.1 | grep "latmin =" | cut -d ' ' -f 4)
                      lonmax=$( ncdump -v lonmax tmp.nc.1  | grep "lonmax =" | cut -d ' ' -f 4)
                      latmax=$( ncdump -v latmax tmp.nc.1  | grep "latmax =" | cut -d ' ' -f 4)
                      \rm tmp.nc.1

                      ncap2 -O -F -s "CPLMASK(1,:,:,:)=(LANDMASK+LAKEMASK-1)*(-1)" tmp.nc tmp.nc.1
                      ncap2 -O -F -s "var_tmp=CPLMASK(1,:,:,:); where( XLAT < $latmin || XLONG < $lonmin || XLAT > $latmax || XLONG > $lonmax ) var_tmp=0; CPLMASK(1,:,:,:)=var_tmp" tmp.nc.1 tmp.nc.1
                      ncks -O -v var_tmp -x tmp.nc.1 tmp.nc.1
                      [[ ! -f tmp.nc.2 ]] && { ncap2 -O -F -s "CPLMASK(1,:,:,:)=(LANDMASK+LAKEMASK-1)*(-1)" tmp.nc tmp.nc.2 ; }
                      ncrcat -O tmp.nc.2 tmp.nc.1 tmp.nc.2
                      num_ext_mod=$( ncdump -h tmp.nc.2 | grep "num_ext_model_couple_dom_stag = " | cut -d ' ' -f 6| cut -c 2)
                      for pp in `seq 1 $(( $num_ext_mod - 1 ))`; do
                          ncap2 -O -F -s "var_tmp=CPLMASK($pp,1,:,:); where( XLAT > $latmin && XLONG > $lonmin && XLAT < $latmax && XLONG < $lonmax ) var_tmp=0; CPLMASK($pp,1,:,:)=var_tmp" tmp.nc.2 tmp.nc.2
                          ncks -O -v var_tmp -x tmp.nc.2 tmp.nc.2
                      done
                  done
              ncpdq -O -a Time,num_ext_model_couple_dom_stag tmp.nc.2 tmp.nc.2
              ncks -A -v CPLMASK -x wrfinput_$dom tmp.nc.2
              mv tmp.nc.2 wrfinput_$dom
              rm -rf tmp.nc tmp.nc.1
              fi
          module unload $ncomod
          else
              module load $ncomod
              echo "set CPLMASK to 1 in coupled domain $dom"
              num_ext_mod=$( ncdump -h wrfinput_d01 | grep "num_ext_model_couple_dom_stag = " | cut -d ' ' -f 3)
              if [[ ${num_ext_mod} -lt $(echo ${wrfcpldom} | wc -w) ]];then
                  echo "Increase size of mum_ext_model by one for ${dom} (value in initial netcdf to small)"
                  ncpdq -O -v CPLMASK -a num_ext_model_couple_dom_stag,Time wrfinput_d01 tmp.nc
                  ncpdq -O -d num_ext_model_couple_dom_stag,0 -v CPLMASK -a num_ext_model_couple_dom_stag,Time wrfinput_d01 tmp2.nc
                  ncrcat -O tmp.nc tmp2.nc tmp2.nc
                  ncpdq -O -a Time,num_ext_model_couple_dom_stag tmp2.nc tmp2.nc
                  ncks -A wrfinput_d01 tmp2.nc
                  mv tmp2.nc wrfinput_d01
                 rm -rf tmp.nc
              fi
              num_ext_mod=$( ncdump -h wrfinput_d01 | grep "num_ext_model_couple_dom_stag = " | cut -d ' ' -f 3)
              echo "Find limits for domain $dom"
              ncap2 -O -v -s 'latmin=XLAT.min();latmax=XLAT.max();lonmin=XLONG.min();lonmax=XLONG.max()' wrfinput_${dom} tmp.nc
              lonmin=$( ncdump -v lonmin tmp.nc  | grep "lonmin =" | cut -d ' ' -f 4)
              latmin=$( ncdump -v latmin tmp.nc  | grep "latmin =" | cut -d ' ' -f 4)
              lonmax=$( ncdump -v lonmax tmp.nc  | grep "lonmax =" | cut -d ' ' -f 4)
              latmax=$( ncdump -v latmax tmp.nc  | grep "latmax =" | cut -d ' ' -f 4)
              rm -rf tmp.nc
              printf "Limits for domain ${dom} are:\n Lon min:$lonmin \n Lat min:$latmin \n Lon max:$lonmax \n Lat max:$latmax \n"
              ncap2 -F -O -s "var_tmp=CPLMASK(:,1,:,:); where( XLAT < $latmin || XLONG < $lonmin || XLAT > $latmax || XLONG > $lonmax ) var_tmp=0; CPLMASK(:,${num_ext_mod},:,:)=var_tmp" wrfinput_d01 wrfinput_d01.tmp
              if [[  ${dom} != "d01" ]]; then
                  ncap2 -F -O -s "var_tmp=CPLMASK(:,1,:,:); where( XLAT > $latmin && XLONG > $lonmin && XLAT < $latmax && XLONG < $lonmax ) var_tmp=0; CPLMASK(:,1,:,:)=var_tmp" wrfinput_d01.tmp wrfinput_d01.tmp
              fi
              ncks -O -v var_tmp -x wrfinput_d01.tmp wrfinput_d01
              rm -rf wrfinput_d01.tmp 
              module unload $ncomod
          fi
      done
  fi
else
    touch ls_l/getfile_atm_restarts.txt
    for file in `${MACHINE_STOCKAGE} ls ${RESTDIR_IN}/wrfrst_d0?_*`
# for i in ${RESTDIR_IN}/wrfrst_d01_*
    do
	${io_getfile} ${file} . >> ls_l/getfile_atm_restarts.txt
    done
fi
