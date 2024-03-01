#!/bin/bash
#===================================
set -u
##set -x

source configure_file
NBPROCS=4 #! see param.h in croco sources 1x4 for REGIONAL



numrev0=`sed -n '/revision/{n;p;}' gitinfos`
numrev=`echo $numrev0 | tr -d [:blank:]`
today=`date +%Y%m%d`
mv Recap_${TEST_NAME}.git${numrev} Recap_${TEST_NAME}_${today}.git${numrev}

#============================================
filein=mpi_${NBPROCS}_${TEST_NAME}_read.log
fileout=checkbugbin_${TEST_NAME}.txt
rm -Rf $fileout
touch $fileout
GREP_CMD='grep -m 1'
execflag=`sed -n 2p ${TEST_NAME}_steps`
echo "execflag is "$execflag

if [ -s $filein ]; then   
    echo "===================" > $fileout
    echo 'CHECK (BUGBIN detection)' >> $fileout
    ${GREP_CMD} BUGBIN $filein >> $fileout
    res=`${GREP_CMD} BUGBIN $filein`
    echo 'res='$res >> $fileout
    if [ -n "$res" ]  ; then 
	echo 'check [Restartability failed]'  >> $fileout
	sed -e '4c N' ${TEST_NAME}_steps > tmp.txt 
    else
	if [[ $execflag == 'Y' ]]; then 
	    echo 'check [Restartability passed]'  >> $fileout
	    sed -e '4c Y' ${TEST_NAME}_steps > tmp.txt
	else
	    echo '... Restartability unknown'  >> $fileout
	    sed -e '4c ?__exec_failure' ${TEST_NAME}_steps > tmp.txt
	fi
    fi
    \mv tmp.txt ${TEST_NAME}_steps
else
    #sed -e '4c ?__comp_failure' ${TEST_NAME}_steps > tmp.txt
    sed -e '4c ?' ${TEST_NAME}_steps > tmp.txt
    \mv tmp.txt ${TEST_NAME}_steps
fi

# if [[ -d $filein_openmp && ! -z "$res_omp" ]] || [[ -d $filein_mpi  && ! -z "$res_mpi" ]] ; then
#     sed -e '3c N' ${TEST_NAME}_steps > tmp.txt 
#     \mv tmp.txt ${TEST_NAME}_steps
#     msg1="      => Repro failure for ${TEST_NAME} ..."
#     msg2="${FMT_REDBLD}${msg1}${FMT_ORD}"
#     echo -e "   $msg2" | tee -a mylog.txt
# else
#     if [ $FLAG_MPI -eq 1 -o  $FLAG_OPENMP -eq 1 ]; then
#  	sed -e '3c Y' ${TEST_NAME}_steps > tmp.txt 
#  	\mv tmp.txt ${TEST_NAME}_steps
#     else
#  	sed -e '3c ?' ${TEST_NAME}_steps > tmp.txt 
#  	\mv tmp.txt ${TEST_NAME}_steps
#     fi
# fi


##cat $fileout checkbugbin_${TEST_NAME}.txt
#============================================
file=Recap_${TEST_NAME}_${today}.git${numrev}
#echo '  ' >> $file
echo ' '
echo '>> oooooooooooooo <<' >> $file
echo 'REVISION GIT :' $numrev >> $file
echo 'DATE         :' $today >> $file
cat checkbugbin_${TEST_NAME}.txt >> $file 
echo '---------------------------' >> $file
echo ' ' >> $file

# Cleaning
rm -f  $fileout
