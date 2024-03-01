#!/bin/bash

##set -x
##set -e
set -u

#echo '============================================================='
echo 'Create the link between test_repro sources and test dir'

source "$CVTKHOME/CONFIGURE_REG"

rm -Rf $dir_test
mkdir -p $dir_test/Junk
[[ ! -d  $dir_web ]] && mkdir -p $dir_web

#
\cp -rf $CI_PROJECT_DIR/OCEAN/croco.in* $CVTKHOME/../common/TEST_CASES_CVTK/VHR/
\cp -rf $CI_PROJECT_DIR/OCEAN/AGRIF_FixedGrids.in $CVTKHOME/../common/TEST_CASES_CVTK/VHR/

for file in $(ls $CVTKHOME/../common/TEST_CASES_CVTK/VHR/croco.in)
do 
  line=$(($(grep -n 'history:' $file  |  awk -F ':' '{print $1}') +1))
  [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
  [ ! -z $toto ] && sed -e "${line} s/$toto/3/" $file > tmp.txt && \mv tmp.txt $file

  line=$(($(grep -n 'restart:' $file  |  awk -F ':' '{print $1}') +1))
  [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
  [ ! -z $toto ] && sed -e "${line} s/$toto/3/" $file > tmp.txt && \mv tmp.txt $file
  
  line=$(($(grep -n 'time_stepping:' $file  |  awk -F ':' '{print $1}') +1))
  [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
  [ ! -z $toto ] && sed -e "${line} s/$toto/3/" $file > tmp.txt && \mv tmp.txt $file

  line=$(($(grep -n 'time_stepping:' $file  |  awk -F ':' '{print $1}') +1))
  [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
  [ ! -z $toto ] && sed -e "${line} s/$toto/900/" $file > tmp.txt && \mv tmp.txt $file
  
done

for file in $(ls $CVTKHOME/../common/TEST_CASES_CVTK/VHR/croco.in.1)
do 
  line=$(($(grep -n 'history:' $file  |  awk -F ':' '{print $1}') +1))
  [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
  [ ! -z $toto ] && sed -e "${line} s/$toto/9/" $file > tmp.txt && \mv tmp.txt $file

  line=$(($(grep -n 'restart:' $file  |  awk -F ':' '{print $1}') +1))
  [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $1}')
  [ ! -z $toto ] && sed -e "${line} s/$toto/9/" $file > tmp.txt && \mv tmp.txt $file

    line=$(($(grep -n 'time_stepping:' $file  |  awk -F ':' '{print $1}') +1))
  [ ! -z $line ] && toto=$(sed -n ${line}p   $file   | awk '{print $2}')
  [ ! -z $toto ] && sed -e "${line} s/$toto/300/" $file > tmp.txt && \mv tmp.txt $file
  
done

for file in $(ls $CVTKHOME/../common/TEST_CASES_CVTK/VHR/AGRIF_FixedGrids.in )
do 
    #sed '2c 79 137 37 117 3 3 3 3' $file > tmp.txt && \mv tmp.txt $file
    sed "2c $nest_position_reg" $file > tmp.txt && \mv tmp.txt $file
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
ln -sf $dir_home/Configure_Test_reg $dir_test/
ln -sf $dir_home/../CONFIGURE_REG $dir_test/
ln -sf $dir_home/../mk_TestDIR.bash $dir_test/
ln -sf $dir_home/../mk_TESTALL.bash $dir_test/
ln -sf $dir_home/../gather_recap.bash $dir_test/

# cleaning
rm -Rf $dir_test/Configure_Test; 
cd $dir_test	
ln -sf Configure_Test_reg Configure_Test
cd -
#
echo 'Well done: Finish linking'
