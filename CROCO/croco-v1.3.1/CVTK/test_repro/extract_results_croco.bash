#!/bin/bash
#===================================
set -u
##set -x

source configure_file



numrev0=`sed -n '/revision/{n;p;}' gitinfos`
numrev=`echo $numrev0 | tr -d [:blank:]`
today=`date +%Y%m%d`
mv Recap_${TEST_NAME}.git${numrev} Recap_${TEST_NAME}_${today}.git${numrev}

#============================================
filein_openmp=openmp_${NBPROCS}_${TEST_NAME}.log
filein_mpi=mpi_${NBPROCS}_${TEST_NAME}.log

fileout_openmp=openmp_checkbugbin_${TEST_NAME}.txt
fileout_mpi=mpi_checkbugbin_${TEST_NAME}.txt

rm -Rf $fileout_openmp $fileout_mpi

touch $fileout_openmp
touch $fileout_mpi

GREP_CMD='grep -m 1'

execflag=`sed -n 2p ${TEST_NAME}_steps`
echo "execflag is "$execflag

if [ $FLAG_OPENMP -eq 1 ]; then
    if [ -s $filein_openmp ]; then
	#=============================================
	# OPENMP
	echo "===================" > $fileout_openmp
	echo 'OPENMP CHECK (BUGBIN detection)' >> $fileout_openmp
	${GREP_CMD} BUGBIN $filein_openmp >> $fileout_openmp
	res_omp=`${GREP_CMD} BUGBIN $filein_openmp`
	echo 'res_omp='$res_omp >> $fileout_openmp
	if [ -n "$res_omp" ] ; then 
	    echo 'check [Parallel reproducibility failed]'  >> $fileout_openmp
	    sed -e '3c N' ${TEST_NAME}_steps > tmp.txt 
	else
	    if [[ $execflag == 'Y' ]]; then 
		echo 'check [Parallel reproducibility passed]'  >> $fileout_openmp
		sed -e '3c Y' ${TEST_NAME}_steps > tmp.txt 	
	    else
		echo '... Parallel reproducibility unknown'  >> $fileout_openmp
		sed -e '3c ?__exec_failure' ${TEST_NAME}_steps > tmp.txt
	    fi
	fi
	\mv tmp.txt ${TEST_NAME}_steps
   else
        #sed -e '3c ?__comp_failure' ${TEST_NAME}_steps > tmp.txt
        sed -e '3c ?' ${TEST_NAME}_steps > tmp.txt
        \mv tmp.txt ${TEST_NAME}_steps
   fi
else
    res_omp=""
    sed -e '3c ?_no_omp_test' ${TEST_NAME}_steps > tmp.txt 
    #sed -e '3c ?' ${TEST_NAME}_steps > tmp.txt 
    \mv tmp.txt ${TEST_NAME}_steps	
fi

#

if [ $FLAG_MPI -eq 1 ]; then
    if [ -s $filein_mpi ]; then
	#MPI
	echo "===================" > $fileout_mpi
	echo 'MPI CHECK (BUGBIN detection)' >> $fileout_mpi
	${GREP_CMD} BUGBIN $filein_mpi >> $fileout_mpi
	res_mpi=`${GREP_CMD} BUGBIN $filein_mpi`
	echo 'res_mpi='$res_mpi >> $fileout_mpi
	if [ -n "$res_mpi" ]  ; then 
	    echo 'check mpi [Parallel reproducibility failed]'  >> $fileout_mpi
	    sed -e '4c N' ${TEST_NAME}_steps > tmp.txt 
	else
	    if [[ $execflag == 'Y' ]]; then 
		echo 'check mpi [Parallel reproducibility passed]'  >> $fileout_mpi
		sed -e '4c Y' ${TEST_NAME}_steps > tmp.txt
	    else
		echo '... Parallel reproducibility unknown'  >> $fileout_mpi
		sed -e '4c ?__exec_failure' ${TEST_NAME}_steps > tmp.txt
	    fi
	fi
	\mv tmp.txt ${TEST_NAME}_steps
    else
        #sed -e '4c ?__comp_failure' ${TEST_NAME}_steps > tmp.txt
        sed -e '4c ?' ${TEST_NAME}_steps > tmp.txt
        \mv tmp.txt ${TEST_NAME}_steps
    fi
else
    res_mpi=""
    sed -e '4c ?_no_mpi_test' ${TEST_NAME}_steps > tmp.txt 
    #sed -e '4c ?' ${TEST_NAME}_steps > tmp.txt 
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

# MERGE
cat $fileout_mpi >>  $fileout_openmp
mv $fileout_openmp checkbugbin_${TEST_NAME}.txt

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
rm -f $fileout_openmp $fileout_mpi
