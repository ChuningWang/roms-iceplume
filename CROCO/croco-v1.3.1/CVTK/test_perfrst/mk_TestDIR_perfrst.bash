#!/bin/bash

##set -x
##set -e
###DO NOT USE SET -E !
##set -u

# $1 => CONFIGURE_ANA or CONFIGURE_VORT
# $2 => vort ou reg
# $3 => BASIN, CANYON, .... (NoAgrif, AGRIF1W, AGRIF2W)

configfile=$1 # $1 => CONFIGURE_ANA or CONFIGURE_VORT
scripttype=$2 # $2 => _vort ou _reg
testcase=$3   # $3 => BASIN, CANYON, .... (NoAgrif, AGRIF1W, AGRIF2W)

#usage : ./mk_TestDIR_perfrst.bash CONFIGURE_ANA _ana BASIN
#############################################

source $configfile
echo "   - SOURCE_CROCO="$SOURCE_CROCO

[ ! -d gitinfos ] && ./git_process.bash

numrev0=`sed -n '/revision/{n;p;}' gitinfos`
numrev=`echo $numrev0 | tr -d [:blank:]`
echo  "   - Testing CROCO Rev$numrev"
# Directory creation
#echo "======================================================"
[ ! -d $testcase ] && mkdir $testcase
echo '   - Create and setup dir. : ' $testcase

rm -Rf Configure_Test ; ln -sf  Configure_Test_${scripttype}_perfrst Configure_Test

cp Configure_Test/$testcase $testcase/configure_file 
cp gitinfos $testcase
ln -sf $dir_home/create_link_for_Test_${scripttype}_perfrst.sh $testcase
#===========================================================
cd $testcase ;  
export mytest=$testcase
./create_link_for_Test_${scripttype}_perfrst.sh

./test_croco_perfrst.sh
cd - >/dev/null

#echo "======================================================"
#echo "  "

