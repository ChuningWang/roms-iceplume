#!/bin/bash

#- dependancies
REQUIRE="matlab pdfcrop gs"
for i in $REQUIRE
do 
 has_it=$(which $i)
 if [ -z $has_it ]; then
 echo -e "\033[1;31m $i NOT available ... We quit \033[0m" && exit 1
 fi
done

b_n=$(basename ${0})
OPTIND=1

x_n='BASIN CANYON EQUATOR INNERSHELF INTERNAL IGW RIVER SEAMOUNT SHELFRONT SOLITON THACKER OVERFLOW UPWELLING VORTEX JET SHOREFACE SANDBAR RIP SWASH TANK GRAV_ADJ ISOLITON KH_INST TIDAL_FLAT DUNE'

#SINGLE_COLUMN PLUME MOVING_BATHY ACOUSTIC TS_HADV_TEST SED_TOY 

x_d=$(dirname $(dirname $PWD))
x_p="NO"
x_m="1"
x_r="Run"
x_s="mpirun"

#- Choice of the options ---
while getopts :hn:d:p:u:m:s:r: V
do
  case $V in
        (h) x_h=${OPTARG};
        echo "Usage      : "${b_n} \
            " [-h] [-n EXAMPLE] [-d ROOT_DIR] [-r RUNDIR] [-p PARALLEL] [-m MAX_PROC]";
        echo " -h               : help";       
        echo " -n EXAMPLE       : TEST name, as listed in cppdefs.h, default : all";
        echo " -d ROOTDIR       : Root of the git repository, default : same as CVTK";
        echo " -r RUNDIR        : Run repository for cppdefs.hand modified .F files, default : Run";
        echo " -p PARALLEL      : Type of parallelism (MPI or OPENMP), default : no";
        echo " -m MAX_PROC      : Max number of cpus available, default : 1";
        echo " -s MPI_RUN       : mpirun command, default : mpirun";
        echo "";
        exit 0;;
        (n)  x_n=${OPTARG};;
        (d)  x_d=${OPTARG};;
        (p)  x_p=${OPTARG};;
        (m)  x_m=${OPTARG};;
        (r)  x_r=${OPTARG};;
        (s)  x_s=${OPTARG};;
        (:)  echo ${b_n}" : -"${OPTARG}" option : missing value" 1>&2;
        exit 2;;
        (\?) echo ${b_n}" : -"${OPTARG}" option : not supported" 1>&2;
        exit 2;;
  esac
done
shift $(($OPTIND-1));

LIST_EXAMPLE=$x_n
ROOTDIR=$x_d
PARALLEL=$x_p
MAX_PROC=$x_m
MPIRUN=$x_s
RUNDIR=$x_r
    
[ ! -d ${ROOTDIR}/${RUNDIR} ] && ../../create_config.bash -f -n $RUNDIR
    
[ -f cppdefs.h ] && \rm cppdefs.h && \cp ${ROOTDIR}/OCEAN/cppdefs.h .
[ -f param.h ]    && \rm param.h   && cp ${ROOTDIR}/OCEAN/param.h .
[ -d TESTCASES ] && \rm -rf TESTCASES 
cp -r ${ROOTDIR}/TEST_CASES .
[ ! -d LOG ] && mkdir LOG

# Number of cases
NB_TEST=$(echo $LIST_EXAMPLE |wc -w )

i=0
for EXAMPLE in $LIST_EXAMPLE
  do
    ((i=$i+1))

    echo '----------------------'
    echo 
    echo "RUNNING TEST Number "${i} / ${NB_TEST}" : $EXAMPLE"
    [ -f LOG/${EXAMPLE}_run.log ] && \rm LOG/${EXAMPLE}_run.log
     ./run_real.sh -n $EXAMPLE -d $ROOTDIR -p $PARALLEL -m $MAX_PROC -r $RUNDIR -s $MPIRUN > LOG/${EXAMPLE}_run.log  2>&1  || { echo COMPILING or RUNNING TEST $EXAMPLE failed... EXITING... && exit ; }
 #   [ ! -f  ${EXAMPLE}_his.nc ] && { echo RUNNING TEST $EXAMPLE failed... EXITING... && exit ; }
    echo 

    echo '----------------------'
    echo 
    echo "PLOTTING TEST $EXAMPLE"
    [ -f LOG/${EXAMPLE}_plot.log ] && \rm LOG/${EXAMPLE}_plot.log
    ./plot_real.sh -n $EXAMPLE -d $ROOTDIR $i > LOG/${EXAMPLE}_plot.log  2>&1  || { echo PLOTTING TEST $EXAMPLE failed... EXITING... && exit ; } 
    echo 
  done

