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
# Script to Run CVTK DEBUG procedure managing parallelization type 
# AND AGRIF nesting type (No nesting, Nesting 1-way, Nesting 2-ways) : 
# VORTEX and REGIONAL
#--------------------------------------------------------------------
##set -x
##set -e

echo "=============================================="
echo "=> CONFIG "$mytest

echo "Remove *.exe* *.log* "
[ ! -z "$(ls *.exe* 2>/dev/null)" ] && /bin/rm *.exe*
[ ! -z "$(ls *.log* 2>/dev/null)" ] && /bin/rm *.log*
echo "Remove the CHECKFILE"
[ -f check_file ] && /bin/rm check_file


#=============================================================================================
#===
source CONFIGURE_GLOBAL
#===
echo " "
echo "=> MPIRUN COMMAND: "$MPIRUN
#
SOURCE_CVTK=${SOURCE_CROCO}/../CVTK/test_repro
echo 'Sources CVTK tests: '$SOURCE_CVTK

source configure_file

#
# Create param.h and cppdefs.h for SERIAL, OPENMP and MPI
#
\cp ${SOURCE_CROCO}/cppdefs.h cppdefs.h.SERIAL
\cp ${SOURCE_CROCO}/param.h param.h.SERIAL

\cp ${SOURCE_CROCO}/cppdefs.h cppdefs.h.OPENMP
\cp ${SOURCE_CROCO}/param.h param.h.OPENMP

\cp ${SOURCE_CROCO}/cppdefs.h cppdefs.h.MPI
\cp ${SOURCE_CROCO}/param.h param.h.MPI


#List of test cases with only one points in one direction
LIST_2DV_X='GRAV_ADJ IGW INNERSHELF INTERNAL SHOREFACE SWASH THACKER TANK ISOLITON KH_INST SANDBAR'
LIST_2DV_Y='OVERFLOW SHELFRONT'

Is2DV_X=0
Is2DV_Y=0
[ -n "$(echo $LIST_2DV_X |grep "${CONFIG_NAME}")" ] && Is2DV_X=1
[ -n "$(echo $LIST_2DV_Y |grep "${CONFIG_NAME}")" ] && Is2DV_Y=1

##############################################################################
#
#   FILL THE CPPDEFS.H 
# Title
echo TESTS OF $CONFIG_NAME
#
# 1- UNDEF ALL THE KEYS
#
echo 'undef '$LIST_KEY0
for par in SERIAL OPENMP MPI ; do 
    echo 'PARA=' $par
    for EXAMPLE in $LIST_KEY0
    do
	sed '/'${EXAMPLE}[[:graph:]]'/! s/'define\ \ \*$EXAMPLE'/'undef\ $EXAMPLE'/' < cppdefs.h.$par > cppdefs.h.$par.tmp
	\mv cppdefs.h.$par.tmp cppdefs.h.$par
    done
    
    # 2- DEFINE THE TYPE OF DEBUG TEST 
    sed '/'${KEY_DEBUG}[[:graph:]]'/! s/'undef\ \ \*$KEY_DEBUG'/'define\ $KEY_DEBUG'/' < cppdefs.h.$par > cppdefs.h.$par.tmp
    \mv cppdefs.h.$par.tmp cppdefs.h.$par
    
    
    # 3- DEFINE THE NAME OF THE CONFIG
    sed '/'${EXAMPLE}[[:graph:]]'/! s/'undef\ \*BENGUELA_LR'/'define\ $CONFIG_NAME'/' < cppdefs.h.$par > cppdefs.h.$par.tmp
    \mv cppdefs.h.$par.tmp cppdefs.h.$par
     
    # 4- DEFINE THE VARIOUS CPPKEYS
    #=4.1
    for EXAMPLE in $LIST_KEY_PHYS
    do
	sed '/'${EXAMPLE}[[:graph:]]'/! s/'undef\ \ \*$EXAMPLE'/'define\ $EXAMPLE'/' < cppdefs.h.$par > cppdefs.h.$par.tmp
	\mv cppdefs.h.$par.tmp cppdefs.h.$par
    done
    #==4.2
    for EXAMPLE in $LIST_KEY_NEST
    do
	echo $EXAMPLE
	sed '/'${EXAMPLE}[[:graph:]]'/! s/'undef\ \ \*$EXAMPLE'/'define\ $EXAMPLE'/' < cppdefs.h.$par > cppdefs.h.$par.tmp
	\mv cppdefs.h.$par.tmp cppdefs.h.$par
    done
done

#=====================================================================================================
#=====================================================================================================
#####
echo 'FLAG OPENMP' ${FLAG_OPENMP}
echo 'FLAG MPI' ${FLAG_MPI}
####
#echo ' '
echo '==============================='
echo 'START TESTING ...             '
#echo '==============================='

SUCCESS=0
SUCCESS_COMP=0
SUCCESS_COMP_SERIAL=0
SUCCESS_COMP_OPENMP=0
SUCCESS_COMP_MPI=0
SUCCESS_EXEC=0
SUCCESS_EXEC_SERIAL=0
SUCCESS_EXEC_OPENMP=0
SUCCESS_EXEC_MPI=0

if [ ! -f ${TEST_NAME}_steps ]; then 
    echo 'Y' > ${TEST_NAME}_steps
    echo 'Y' >> ${TEST_NAME}_steps
    echo 'Y' >> ${TEST_NAME}_steps
    echo 'Y' >> ${TEST_NAME}_steps
fi


##############################################################################
# Serial runs
##############################################################################
echo "KEYS TESTED : "$LIST_KEY_PHYS
#echo ''
par1='SERIAL'
echo "SERIAL NPP=1 TEST $mytest"
[ -f check_file ] && /bin/rm check_file
sed 's/'NSUB_X=1,\ \ \*NSUB_E=NPP'/'NSUB_X=1,\ NSUB_E=1'/' < param.h.$par1 > param.h.$par1.tmp
sed 's/'NPP=1'/'NPP=1'/' < param.h.$par1 > param.h.$par1
\mv param.h.$par1.tmp param.h.$par1

[ -e  param.h.OK ] && \rm param.h.OK
[ -e  param.h.OK ] && \rm cppdefs.h.OK

\cp param.h.$par1 param.h.OK
\cp cppdefs.h.$par1 cppdefs.h.OK

Fqsub_serial
myreturn=$?
echo "myreturn is "$myreturn
#echo "SUCCESS="$SUCCESS
SUCCESS=$(($SUCCESS+$myreturn))
#echo "SUCCESS="$SUCCESS
if [ "$myreturn" -eq 1 ]; then
    SUCCESS_COMP_SERIAL=$(($SUCCESS_COMP_SERIAL+1))
    sed -e '1c N' ${TEST_NAME}_steps > tmp.txt 
    \mv tmp.txt ${TEST_NAME}_steps 
    sed -e '2c ?' ${TEST_NAME}_steps > tmp.txt 
    \mv tmp.txt ${TEST_NAME}_steps 
fi  
if [ "$myreturn" -eq 2 ]; then
    SUCCESS_EXEC_SERIAL=$(($SUCCESS_EXEC_SERIAL+1))
    sed -e '2c N' ${TEST_NAME}_steps > tmp.txt 
    \mv tmp.txt ${TEST_NAME}_steps 
fi  

# 4- 
##############################################################################
# Openmp runs
##############################################################################
if [ ${FLAG_OPENMP} -eq 1 ]; then 
    
    par1='OPENMP'
    if [ $Is2DV_Y == 1 ]; then
	echo " "
	echo "OPEN-MP 1X2 NPP=2 TEST $mytest"
	sed 's/'NSUB_X=1,\ \ \*NSUB_E=NPP'/'NSUB_X=1,\ NSUB_E=2'/' < param.h.$par1 > param.h.$par1.tmp
	sed 's/'NPP=4'/'NPP=2'/' < param.h.$par1.tmp > param.h.$par1
    elif [ $Is2DV_X == 1 ]; then
	echo " "
	echo "OPEN-MP 2x1 NPP=2 TEST $mytest"
	sed 's/'NSUB_X=1,\ \ \*NSUB_E=NPP'/'NSUB_X=2,\ NSUB_E=1'/' < param.h.$par1 > param.h.$par1.tmp
	sed 's/'NPP=4'/'NPP=2'/' < param.h.$par1.tmp > param.h.$par1
    else
	echo " "
	echo "OPEN-MP ${NBPROCS_X}X${NBPROCS_Y} TEST $mytest"
	sed 's/'NSUB_X=1,\ \ \*NSUB_E=NPP'/'NSUB_X=${NBPROCS_X},\ NSUB_E=${NBPROCS_Y}'/' < param.h.$par1 > param.h.$par1.tmp
	sed 's/'NPP=4'/'NPP=$(( $NBPROCS_X * $NBPROCS_Y ))'/' < param.h.$par1.tmp > param.h.$par1
    fi
    #
    sed '/'${par1}[[:graph:]]'/!s/'undef\ \ \*${par1}'/'define\ ${par1}'/' < cppdefs.h.$par1 > cppdefs.h.$par1.tmp
    \mv cppdefs.h.$par1.tmp cppdefs.h.$par1
    #
    [ -e  param.h.OK ] && \rm param.h.OK
    [ -e  cppdefs.h.OK ] && \rm cppdefs.h.OK
    \cp param.h.$par1 param.h.OK
    \cp cppdefs.h.$par1 cppdefs.h.OK
    #
    Fqsub_openmp
    myreturn=$?
    echo "myreturn is "$myreturn
#    echo "SUCCESS="$SUCCESS
    SUCCESS=$(($SUCCESS+$myreturn))
#    echo "SUCCESS="$SUCCESS
    if [ "$myreturn" -eq 1 ]; then
	SUCCESS_COMP_OPENMP=$(($SUCCESS_COMP_OPENMP+1))
	sed -e '1c N' ${TEST_NAME}_steps > tmp.txt
	\mv tmp.txt ${TEST_NAME}_steps 
	sed -e '2c ?' ${TEST_NAME}_steps > tmp.txt 
	\mv tmp.txt ${TEST_NAME}_steps 
    fi  
    if [ "$myreturn" -eq 2 ]; then
	SUCCESS_EXEC_OPENMP=$(($SUCCESS_EXEC_OPENMP+1))
	sed -e '2c N' ${TEST_NAME}_steps > tmp.txt 
	\mv tmp.txt ${TEST_NAME}_steps 
    fi  

    #else
    #  sed -e '2c ?' ${TEST_NAME}_steps > tmp.txt 
    #  \mv tmp.txt ${TEST_NAME}_steps     
fi 

###############################################################################
# 4- RVTK_DEBUG_REG_DEV
##############################################################################
# Mpi runs
##############################################################################
if [ ${FLAG_MPI} -eq 1 ]; then 
    
    par1='MPI'
    if [ $Is2DV_Y == 1 ]; then
	echo " "
	echo "MPI 1X2 TEST $mytest"
	sed 's/'NP_XI=1,\ \ \*NP_ETA=4'/'NP_XI=1,\ NP_ETA=2'/' < param.h.$par1 > param.h.$par.tmp
	
    elif [ $Is2DV_X == 1 ]; then
	echo " "
	echo "MPI 2X1 TEST $mytest"
	sed 's/'NP_XI=1,\ \ \*NP_ETA=4'/'NP_XI=2,\ NP_ETA=1'/' < param.h.$par1 > param.h.$par.tmp
    else
	echo " "
	echo "MPI ${NBPROCS_X}X${NBPROCS_Y} TEST $mytest"
	sed 's/'NP_XI=1,\ \ \*NP_ETA=4'/'NP_XI=${NBPROCS_X},\ NP_ETA=${NBPROCS_Y}'/' < param.h.$par1 > param.h.$par1.tmp
    fi
    \mv param.h.$par1.tmp param.h.$par1
    #
    sed '/'${par1}[[:graph:]]'/!s/'undef\ \ \*$par1'/'define\ $par1'/' < cppdefs.h.$par1 > cppdefs.h.$par1.tmp
    \mv cppdefs.h.$par1.tmp cppdefs.h.$par1
    #
    [ -e  param.h.OK ] && \rm param.h.OK
    [ -e  cppdefs.h.OK ] && \rm cppdefs.h.OK

    \cp param.h.$par1 param.h.OK
    \cp cppdefs.h.$par1 cppdefs.h.OK
    #
    Fqsub_mpi
    myreturn=$?
    echo "myreturn is "$myreturn
#    echo "SUCCESS="$SUCCESS
    SUCCESS=$(($SUCCESS+myreturn))
#    echo "SUCCESS="$SUCCESS
    if [ "$myreturn" -eq 1 ]; then
	SUCCESS_COMP_MPI=$(($SUCCESS_COMP_MPI+1))
	sed -e '1c N' ${TEST_NAME}_steps > tmp.txt 
	\mv tmp.txt ${TEST_NAME}_steps 
	sed -e '2c ?' ${TEST_NAME}_steps > tmp.txt 
	\mv tmp.txt ${TEST_NAME}_steps 
    fi  
    if [ "$myreturn" -eq 2 ]; then
	SUCCESS_EXEC_MPI=$(($SUCCESS_EXEC_MPI+1))
	sed -e '2c N' ${TEST_NAME}_steps > tmp.txt 
	\mv tmp.txt ${TEST_NAME}_steps 
    fi  
fi

##############################################################################
echo "&&&&&&&&&&&&&&&&&&&&&&&&&"
echo "Final SUCCESS is "$SUCCESS
echo " "
if [  "$SUCCESS" -ne 0 ]; then
    #sed not needed 
    sed -e '3c ?' ${TEST_NAME}_steps > tmp.txt ; \mv tmp.txt ${TEST_NAME}_steps
    sed -e '4c ?' ${TEST_NAME}_steps > tmp.txt ; \mv tmp.txt ${TEST_NAME}_steps
    #echo
    echo "Final SUCCESS -ne 0 => "
    echo "      SOMETHING WRONG HAPPENED WITH ${CONFIG_NAME}"
    #echo "EXITING ..."
    # echo
    #echo  | tee -a mylog.txt
    #echo -e "Final SUCESS is "$SUCCESS | tee -a mylog.txt
    #echo -e "Final SUCESS_COMP is "$SUCCESS_COMP | tee -a mylog.txt
    #echo -e "Final SUCESS_EXE is "$SUCCESS_EXE | tee -a mylog.txt
    #echo -e "${FMT_REDBLD}SOMETHING WRONG HAPPENED WITH ${CONFIG_NAME} ${FMT_ORD}" | tee -a mylog.txt
    #echo -e "${FMT_REDBLD}EXITING ...${FMT_ORD}"  | tee -a mylog.txt
    # if [ "$SUCCESS_COMP" -ne 0 ]; then
    # 	echo -e "${FMT_REDBLD}A COMPILATION ERROR WITH ${CONFIG_NAME} ${FMT_ORD}" | tee -a mylog.txt
    # fi
    # if [ "$SUCCESS_COMP" -eq 0 ] &&  [ "$SUCCESS_EXE" -ne 0 ]; then
    # 	echo -e "${FMT_REDBLD}COMPILATION IS OK WITH ${CONFIG_NAME} ${FMT_ORD}" | tee -a mylog.txt
    # 	echo -e "${FMT_REDBLD}AN EXECUTION ERROR WITH ${CONFIG_NAME} ${FMT_ORD}" | tee -a mylog.txt
    # fi
  
    ##echo  | tee -a mylog.txt
    ##exit  1
fi
#########################################################################################################


#########################################################################################################
# 6 - Extract results
##############################################################################
#  runs
echo "EXTRACTION $mytest"
Fextract_results $FLAG_MPI $FLAG_OPENMP
