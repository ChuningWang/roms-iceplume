#!/bin/bash
#===================================

set -u
###set -x
###DO NOT USE SET -E !

source CONFIGURE_GLOBAL_PERFRST
source configure_file

#cd $SUBMIT_DIR
#echo "$SUBMIT_DIR"

mkdir CROCO_OUTFILES
NBPROCS=4 #! see param.h in croco sources 1x4 for REGIONAL
#===================================
par1='write'
\cp cppdefs.h.exactrestart.${par1} cppdefs.h.OK
# 1.1 Compile
msg1="- Compilation failure for ${TEST_NAME} : ${par1}..."
msg2="${FMT_REDBLD}${msg1}${FMT_ORD}"
./jobcomp_rvtk.bash Compile_$par1 > jobcomp_${par1}_${TEST_NAME}.log  2>&1 || { echo -e "   $msg2" | tee -a mylog.txt ; echo -e $msg1 ; exit 1 ; }
/bin/mv croco croco_${par1}.exe

# 1.2 Run
msg1="- Execution failure for ${TEST_NAME} : ${par1}..."
msg2="${FMT_REDBLD}${msg1}${FMT_ORD}"

#msg3="- Parallel repro. failure for ${TEST_NAME} : ${par1}..."
#msg4="${FMT_REDBLD}${msg3}${FMT_ORD}"

$MPIRUN -np $NBPROCS ./croco_${par1}.exe croco.in.write > mpi_${NBPROCS}_${TEST_NAME}_$par1.log 2>&1  || { echo -e "   $msg2" | tee -a mylog.txt ; echo -e $msg1 ; exit 2 ; }
write_perfrst_exec_status=$?
echo "write_perfrst_exec_status is "$write_perfrst_exec_status
#  =0, OK or  clean stop before the end (bugbin or blow up)
# !=0, KO and bad stop before the end (input problem)


#=========================================================================================
#=========================================================================================

par1='read'
\cp cppdefs.h.exactrestart.${par1} cppdefs.h.OK
# 2.1 Compile
msg1="- Compilation failure for ${TEST_NAME} : ${par1}..."
msg2="${FMT_REDBLD}${msg1}${FMT_ORD}"
./jobcomp_rvtk.bash Compile_$par1 > jobcomp_${par1}_${TEST_NAME}.log  2>&1 || { echo -e "   $msg2" | tee -a mylog.txt ; echo -e $msg1 ; exit 3 ; }
/bin/mv croco croco_${par1}.exe

# 2.2 Run
msg1="- Execution failure for ${TEST_NAME} : ${par1}..."
msg2="${FMT_REDBLD}${msg1}${FMT_ORD}"

msg3="- Perfrst failure for ${TEST_NAME} : ${par1}..."
msg4="${FMT_REDBLD}${msg3}${FMT_ORD}"

$MPIRUN -np $NBPROCS ./croco_${par1}.exe croco.in.read > mpi_${NBPROCS}_${TEST_NAME}_$par1.log 2>&1  || { echo -e "   $msg2" | tee -a mylog.txt ; echo -e $msg1 ; exit 4 ; }
read_perfrst_exec_status=$?
echo "read_perfrst_exec_status is "$read_perfrst_exec_status
#  =0, OK or  clean stop before the end (bugbin or blow up)
# !=0, KO and bad stop before the end (input problem)

# 2.3 : bugbin detection in the read part
#
grep 'BUGBIN' mpi_${NBPROCS}_${TEST_NAME}_read.log > /dev/null 2>&1
bugbin_detec=$?
echo "bugbin detection flag is " $bugbin_detec
# bugbin_detection = 0  mean we have a problem; KO
# bugbin_detection !=0  mean the test is OK

if [ $read_perfrst_exec_status -eq 0 ] && [ $read_perfrst_exec_status -eq 0 ] && [ $bugbin_detec -eq 0 ] ; then
#if [ $bugbin_detec -eq 0 ] ; then
    # sortie car erreur de repro (bugbin detection)
    echo -e "   $msg4" | tee -a mylog.txt
    echo -e $msg3
    exit 3 
fi
#echo 'Output message status is' $?

#===
# # FAKE CHANGE #2
# # Additional check in case of clean stop before the end, mainly Parallel repro failure case
# #     Clean stop in case of blow up or bugbin => exit status =0
# #     Bad stop in case of bad input or else  => exit status !=0

# SUCCESS_TMP=1
# grep 'MAIN: DONE'  mpi_${NBPROCS}_${TEST_NAME}.log || SUCCESS_TMP=0
# if [  "$SUCCESS_TMP" -eq 0 ]; then
#   echo -e "   $msg4" | tee -a mylog.txt
#   echo -e $msg3 
#   exit 2 
# fi	
#===
