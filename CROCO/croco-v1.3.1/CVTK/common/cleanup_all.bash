#!/bin/bash
###===================================
###PBS -q sequentiel
###PBS -l walltime=02:00:00
###PBS -j oe
###PBS -M gildas.cambon@ird.fr -m abe
##cd $PBS_O_WORKDIR
##echo $PBS_O_LOGNAME
###===================================

#set -x

cd ${TESTROOTDIR}/REG
./mk_CLEANALL.bash
cd -

cd ${TESTROOTDIR}/VORT
./mk_CLEANALL.bash
cd -

cd ${TESTROOTDIR}/KTEST
./mk_CLEANALL.bash
cd -
