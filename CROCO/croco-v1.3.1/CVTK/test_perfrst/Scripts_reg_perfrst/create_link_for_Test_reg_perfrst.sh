#!/bin/bash
#echo '============================================================='

##set -x
##set -e
set -u

echo 'Create the link between TESTROOT . and '$PWD

source ../CONFIGURE_REG_PERFRST

#echo 'common scripts
ln -sf ../CONFIGURE_GLOBAL_PERFRST .
ln -sf ../CONFIGURE_REG_PERFRST .

[ -d TEST_CASES ] && rm TEST_CASES
ln -sf $PERFRSTHOME/../common/TEST_CASES_CVTK TEST_CASES
ln -sf $PERFRSTHOME/../common/TEST_CASES_CVTK/VHR/croco.in.read .
ln -sf $PERFRSTHOME/../common/TEST_CASES_CVTK/VHR/croco.in.write .
#ln -sf $dir_home/../TEST_CASES_CVTK/VHR/croco.in.1 .
ln -sf $PERFRSTHOME/../common/TEST_CASES_CVTK/VHR/AGRIF_FixedGrids.in .
ln -sf $PERFRSTHOME/../common/jobcomp_rvtk.bash .

# data for testrepro REG
ln -sf $dir_datafile/CROCO_FILES .
ln -sf $dir_datafile/DATA .

# test perfrst specific scripts
ln -sf $dir_home/../extract_results_croco_perfrst.bash .
ln -sf $dir_home/../comp_run_mpi_perfrst.bash .
ln -sf $dir_home/../test_croco_perfrst.sh .
ln -sf $dir_home/../rvtk_fast_qsub_perfrst.bash .

#echo 'Process namelist files'
cp -Rf jobcomp_rvtk.bash jobcomp_rvtk.bash.BACK
