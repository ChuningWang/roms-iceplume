#!/bin/bash
#echo '============================================================='

##set -x
##set -e
echo 'Create the link between TESTROOT . and '$PWD

source ../CONFIGURE_VORT

#echo 'common scripts
ln -sf ../CONFIGURE_GLOBAL .
ln -sf ../CONFIGURE_VORT .

[ -d TEST_CASES ] && rm TEST_CASES
ln -sf $CVTKHOME/../common/TEST_CASES_CVTK TEST_CASES
ln -sf $CVTKHOME/../common/TEST_CASES_CVTK/croco.in.Vortex croco.in
ln -sf $CVTKHOME/../common/TEST_CASES_CVTK/croco.in.Vortex.1 croco.in.1
ln -sf $CVTKHOME/../common/TEST_CASES_CVTK/AGRIF_FixedGrids.in.VORTEX AGRIF_FixedGrids.in
ln -sf $CVTKHOME/../common/TEST_CASES_CVTK/vortex_grd.nc* .
ln -sf $CVTKHOME/../common/TEST_CASES_CVTK/vortex_ini.nc* .
ln -sf $CVTKHOME/../common/TEST_CASES_CVTK/vortex_clm.nc .
ln -sf $CVTKHOME/../common/jobcomp_rvtk.bash .

# test repro specific scripts
ln -sf $dir_home/../extract_results_croco.bash .
ln -sf $dir_home/../comp_run_*.bash .
ln -sf $dir_home/../test_croco.sh .
ln -sf $dir_home/../rvtk_fast_qsub.bash .

#echo 'Process namelist files'
cp -Rf jobcomp_rvtk.bash jobcomp_rvtk.bash.BACK
