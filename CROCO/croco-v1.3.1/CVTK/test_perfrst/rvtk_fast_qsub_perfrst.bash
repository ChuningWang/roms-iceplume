#!/bin/bash
#
#======================================================================
# CROCO is a branch of ROMS developped at IRD and INRIA, in France
# The two other branches from UCLA (Shchepetkin et al) 
# and Rutgers University (Arango et al) are under MIT/X style license.
# CROCO specific routines (nesting) are under CeCILL-C license.
# 
# CROCO website : http://www.croco-ocean.org
#======================================================================
#
#---------------------------------------------------------------------
# Script to Run PERFRST DEBUG procedure managing
#--------------------------------------------------------------------
##set -x
##set -e
###DO NOT USE SET -E !

echo "=============================================="
echo "=> CONFIG "$mytest

echo "Remove *.exe* *.log* "
[ ! -z "$(ls *.exe* 2>/dev/null)" ] && /bin/rm *.exe*
[ ! -z "$(ls *.log* 2>/dev/null)" ] && /bin/rm *.log*
echo "Remove the CHECKFILE"
[ -f check_file ] && /bin/rm check_file


#=============================================================================================
#===
source CONFIGURE_GLOBAL_PERFRST
#===
echo " "
echo "=> MPIRUN COMMAND: "$MPIRUN
#
SOURCE_CVTK=${SOURCE_CROCO}/../CVTK/test_perfrst
echo 'Sources CVTK tests: '$SOURCE_CVTK
###\cp ${SOURCE_CROCO}/param.h param.h.exactrestart

#List of test cases with only one points in one direction
LIST_2DV_X='GRAV_ADJ IGW INNERSHELF INTERNAL SHOREFACE SWASH THACKER TANK ISOLITON KH_INST SANDBAR'
LIST_2DV_Y='OVERFLOW SHELFRONT'

source configure_file

Is2DV_X=0
Is2DV_Y=0
[ -n "$(echo $LIST_2DV_X |grep "${CONFIG_NAME}")" ] && Is2DV_X=1
[ -n "$(echo $LIST_2DV_Y |grep "${CONFIG_NAME}")" ] && Is2DV_Y=1

##############################################################################
#
# FILL THE CPPDEFS.H

\cp ${SOURCE_CROCO}/cppdefs.h cppdefs.h.exactrestart

# Title
echo TESTS OF $CONFIG_NAME
#
# 1- UNDEF ALL THE KEYS
#
echo 'undef '$LIST_KEY0
for EXAMPLE in $LIST_KEY0 ; do
    sed '/'${EXAMPLE}[[:graph:]]'/! s/'define\ \ \*$EXAMPLE'/'undef\ $EXAMPLE'/' < cppdefs.h.exactrestart > cppdefs.h.exactrestart.tmp
    \mv cppdefs.h.exactrestart.tmp cppdefs.h.exactrestart
done

# 2- DEFINE THE TYPE OF DEBUG TEST
#2.1 RVTK_DEBUG
sed '/'${KEY_DEBUG}[[:graph:]]'/! s/'undef\ \ \*$KEY_DEBUG'/'define\ $KEY_DEBUG'/' < cppdefs.h.exactrestart > cppdefs.h.exactrestart.tmp
\mv cppdefs.h.exactrestart.tmp cppdefs.h.exactrestart

#2.2
# RVTK_DEBUG_PERFRST key
# => trigger #define EXACT_RESTART
sed '/'RVTK_DEBUG_PERFRST[[:graph:]]'/! s/'undef\ \ \*RVTK_DEBUG_PERFRST'/'define\ RVTK_DEBUG_PERFRST'/' < cppdefs.h.exactrestart > cppdefs.h.exactrestart.tmp
mv cppdefs.h.exactrestart.tmp cppdefs.h.exactrestart

#
# MPI key
#
sed '/'MPI[[:graph:]]'/! s/'undef\ \ \*MPI'/'define\ MPI'/' < cppdefs.h.exactrestart > cppdefs.h.exactrestart.tmp
mv cppdefs.h.exactrestart.tmp cppdefs.h.exactrestart

#
# 3- DEFINE THE NAME OF THE CONFIG
sed '/'${EXAMPLE}[[:graph:]]'/! s/'undef\ \*BENGUELA_LR'/'define\ $CONFIG_NAME'/' < cppdefs.h.exactrestart > cppdefs.h.exactrestart.tmp
\mv cppdefs.h.exactrestart.tmp cppdefs.h.exactrestart

# 4- DEFINE THE VARIOUS CPPKEYS
#=4.1
for EXAMPLE in $LIST_KEY_PHYS
do
    sed '/'${EXAMPLE}[[:graph:]]'/! s/'undef\ \ \*$EXAMPLE'/'define\ $EXAMPLE'/' < cppdefs.h.exactrestart > cppdefs.h.exactrestart.tmp
    \mv cppdefs.h.exactrestart.tmp cppdefs.h.exactrestart
done
#==4.2
for EXAMPLE in $LIST_KEY_NEST
do
    echo $EXAMPLE
    sed '/'${EXAMPLE}[[:graph:]]'/! s/'undef\ \ \*$EXAMPLE'/'define\ $EXAMPLE'/' < cppdefs.h.exactrestart > cppdefs.h.exactrestart.tmp
    \mv cppdefs.h.exactrestart.tmp cppdefs.h.exactrestart
done
#
# write and read version
#
mv cppdefs.h.exactrestart cppdefs.h.exactrestart.write
cp cppdefs.h.exactrestart.write cppdefs.h.exactrestart.read
sed '/'XXXRVTK_DEBUG_READ[[:graph:]]'/! s/'define\ \ \*XXXRVTK_DEBUG_READ'/'define\ RVTK_DEBUG_READ'/' < cppdefs.h.exactrestart.read > cppdefs.h.exactrestart.read.tmp
mv cppdefs.h.exactrestart.read.tmp cppdefs.h.exactrestart.read

#echo ' '
echo '==============================='
echo 'START TESTING ...             '
#echo '==============================='

SUCCESS=0
SUCCESS_COMP=0
SUCCESS_EXEC_WRITE=0
SUCCESS_RESTART_READ=0
if [ ! -f ${TEST_NAME}_steps ]; then 
    echo 'Y' > ${TEST_NAME}_steps
    echo 'Y' >> ${TEST_NAME}_steps
    echo 'Y' >> ${TEST_NAME}_steps
fi

#  WRITE case compilation and run
[ -e  cppdefs.h.OK ] && \rm cppdefs.h.OK

##########################################
./comp_run_mpi_perfrst.bash
myreturn=$?
echo "myreturn is $myreturn"
echo $myreturn

SUCCESS=$(($SUCCESS+$myreturn))

if [ "$myreturn" -eq 1 ]; then
    SUCCESS_COMP=$(($SUCCESS_COMP+1))
    sed -e '1c N' ${TEST_NAME}_steps > tmp.txt 
    \mv tmp.txt ${TEST_NAME}_steps 
    sed -e '2c ?' ${TEST_NAME}_steps > tmp.txt 
    \mv tmp.txt ${TEST_NAME}_steps 
fi  
if [ "$myreturn" -eq 2 ]; then
    SUCCESS_EXEC=$(($SUCCESS_EXEC+1))
    sed -e '2c N' ${TEST_NAME}_steps > tmp.txt 
    \mv tmp.txt ${TEST_NAME}_steps 
fi  
##############################################################################

echo "&&&&&&&&&&&&&&&&&&&&&&&&&"
echo "Final SUCCESS is "$SUCCESS
echo " "
if [  "$SUCCESS" -ne 0 ]; then
    #sed not needed 
    sed -e '3c ?' ${TEST_NAME}_steps > tmp.txt ; \mv tmp.txt ${TEST_NAME}_steps
    #echo
    echo "Final SUCCESS -ne 0 => "
    echo "      SOMETHING WRONG HAPPENED WITH ${CONFIG_NAME}"
fi
#########################################################################################################


# #########################################################################################################
# # 6 - Extract results
# ##############################################################################
# #  runs
#echo "EXTRACTION $mytest"
#Fextract_results_perfrst $FLAG_MPI $FLAG_OPENMP
./extract_results_croco_perfrst.bash
