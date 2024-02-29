#!/bin/bash

#set -x
##set -e
set -u

#echo '============================================================='
echo 'Create the link between test_repro sources and test dir'

source "$CVTKHOME/CONFIGURE_ANA"

rm -Rf $dir_test
mkdir -p $dir_test/Junk
[[ ! -d  $dir_web ]] && mkdir -p $dir_web

#
\cp -rf $CI_PROJECT_DIR/TEST_CASES/* $CVTKHOME/../common/TEST_CASES_CVTK/.

for file in $(ls $CVTKHOME/../common/TEST_CASES_CVTK/croco.in*)
do

  line=$(($(grep -n 'history:' $file  |  awk -F ':' '{print $1}') +1))
  [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
  [ ! -z $toto ] && sed -e "${line} s/$toto/6/" $file > tmp.txt && \mv tmp.txt $file

  line=$(($(grep -n 'restart:' $file  |  awk -F ':' '{print $1}') +1))
  [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
  [ ! -z $toto ] && sed -e "${line} s/$toto/6/" $file > tmp.txt && \mv tmp.txt $file

  line=$(($(grep -n 'time_stepping:' $file  |  awk -F ':' '{print $1}') +1))
  [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
  [ ! -z $line ] && mydt=$(sed -n ${line}p   $file   | awk '{print $2}')
  [ ! -z $toto ] && sed -e "${line} s/$toto/6/" $file > tmp.txt && \mv tmp.txt $file

  nb_sec=$( echo $mydt*6|bc )

  # initialisation
  start_date=''
  end_date=''
  lineend=''
  new_end_date=''
  #
  ffstart=$(grep -n 'start_date:' $file)
  if [ ! -z $ffstart ]; then
      linestr=$(($(grep -n 'start_date:' $file  |  awk -F ':' '{print $1}') +1))
      start_date=$(sed -n ${linestr}p   $file)
      nb_sec2=$(echo $nb_sec | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}') # arrondi a l'entier sup
      nb_sec3=$(( nb_sec2 ))
      nb_sec4=$(printf %02d $nb_sec3)
      new_end_date=$( echo  $start_date |cut -c1-17 )$nb_sec4
  fi
  #
  ffend=$(grep -n 'end_date:' $file)
  if [ ! -z $ffend ]; then
      lineend=$(($(grep -n 'end_date:' $file  |  awk -F ':' '{print $1}') +1))
      end_date=$(sed -n ${lineend}p   $file)
      sed -e "${lineend} s%${end_date}%${new_end_date}%g" $file > tmp.txt && \mv tmp.txt $file
  fi
  #
done
 
# CI common scripts and programms
ln -sf $CVTKHOME/../common/CONFIGURE_GLOBAL $dir_test/
ln -sf $CVTKHOME/../common/gitinfo.sh $dir_test/
ln -sf $CVTKHOME/../common/git_process.bash $dir_test/
ln -sf $CVTKHOME/../common/mk_CLEANALL.bash $dir_test/
ln -sf $CVTKHOME/../common/mk_CHECKALL.bash $dir_test/
ln -sf $CVTKHOME/../common/jobcomp_rvtk.bash $dir_test/
ln -sf $CVTKHOME/../common/print/* $dir_test/

# test repro specific + ana specific
ln -sf $dir_web $dir_test/
ln -sf $dir_home/Configure_Test_ana $dir_test/
ln -sf $dir_home/../CONFIGURE_ANA $dir_test/
ln -sf $dir_home/../mk_TestDIR.bash $dir_test/
ln -sf $dir_home/../mk_TESTALL.bash $dir_test/
ln -sf $dir_home/../gather_recap.bash $dir_test/

# cleaning
rm -Rf $dir_test/Configure_Test; 
cd $dir_test	
ln -sf Configure_Test_ana Configure_Test
cd -
#
echo 'Well done: Finish linking'

