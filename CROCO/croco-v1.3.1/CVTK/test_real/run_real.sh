#!/bin/bash
b_n=$(basename ${0})
OPTIND=1

x_n='BASIN CANYON EQUATOR INNERSHELF INTERNAL IGW RIVER SEAMOUNT SHELFRONT SOLITON THACKER OVERFLOW UPWELLING VORTEX JET SHOREFACE SANDBAR RIP SWASH TANK GRAV_ADJ ISOLITON KH_INST TIDAL_FLAT DUNE'
x_d=$(dirname $(dirname $PWD))
x_p="NO"
x_m="1"
x_r="Run"
x_s="mpirun"

while getopts :hn:d:p:u:m:s:r: V
do
  case $V in
        (h) x_h=${OPTARG};
        echo "Usage      : "${b_n} \
            " [-h] [-n EXAMPLE] [-d ROOT_DIR] [-r RUNDIR] [-p PARALLEL] [-m MAX_PROC]";
        echo " -h               : help";
        echo " -n EXAMPLE       : TEST name, as listed in cppdefs.h, default : all";
        echo " -d ROOTDIR       : Root of the git repository, default : same as CVTK";
        echo " -r RUNDIR        : Run repository for cppdefs.h and modified .F files, default : Run";
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
#MPIRUN="mpirun --mca btl self,sm"

[ ! -d ${ROOTDIR}/${RUNDIR} ] && ../../create_config.bash -f -n $RUNDIR
\cp ${ROOTDIR}/${RUNDIR}/jobcomp .
\cp ${ROOTDIR}/${RUNDIR}/*.h .
\cp ${ROOTDIR}/${RUNDIR}/*.F .
\cp ${ROOTDIR}/${RUNDIR}/*.F90 .
[ -f cppdefs.h ] && \rm cppdefs.h && \cp ${ROOTDIR}/${RUNDIR}/cppdefs.h .
[ -f param.h ]   && \rm param.h   && cp ${ROOTDIR}/${RUNDIR}/param.h .
[ -d TEST_CASES ] && \rm -rf TEST_CASES && \cp -r ${ROOTDIR}/${RUNDIR}/TEST_CASES .

#- undef everything
sed  '1,/^#if defined REGIONAL/ s/define /undef /g' cppdefs.h > tmp.txt && \mv tmp.txt cppdefs.h
#- undef MPI
sed  "s/^#\(.*\)define\(.*\)MPI\(.*\)/#undef MPI/ " cppdefs.h > tmp.txt && \mv tmp.txt cppdefs.h
#- undef OPENMP
sed "s/^#\(.*\)define\(.*\)OPENMP\(.*\)/#undef OPENMP/ " cppdefs.h > tmp.txt && \mv tmp.txt cppdefs.h
#- undef NBQ
#sed "s/^#\(.*\)define\(.*\)NBQ\(.*\)/#undef NBQ/ " cppdefs.h > tmp.txt && \mv tmp.txt cppdefs.h

#- suppress totally MPI or OPENMP (mandatory)
sed  -e "s/^#\(.*\)undef\(.*\)MPI$//g" -e "s/^#\(.*\)undef\(.*\)OPENMP$//g" cppdefs.h > tmp.txt && \mv tmp.txt cppdefs.h

#- define MPI if needed
if [ "$PARALLEL" == "MPI" ]; then
  sed   /"defined\(.*\)${EXAMPLE}"'/a\'$'\n''#define MPI'$'\n' cppdefs.h > tmp.txt && \mv tmp.txt cppdefs.h
  sed   /"defined\(.*\)${EXAMPLE}"'/a\'$'\n''#define NC4PAR'$'\n' cppdefs.h > tmp.txt && \mv tmp.txt cppdefs.h
fi

#- define OPENMP if need
if [ "$PARALLEL" == "OPENMP" ]; then
  sed   /"defined\(.*\)${EXAMPLE}"'/a\'$'\n''#define OPENMP'$'\n' cppdefs.h > tmp.txt && \mv tmp.txt cppdefs.h
fi

# proper path in jobcomp
sed  "s:.*SOURCE=.*:SOURCE=$ROOTDIR/OCEAN:g" jobcomp


# main loop
for EXAMPLE in $LIST_EXAMPLE
do
#- define config
  sed  "s/^#\(.*\)undef\(.*\)$EXAMPLE\(.*\)/#define $EXAMPLE/" cppdefs.h  > tmp.txt && \mv tmp.txt cppdefs.h

#- parallel of not
#-  find EXAMPLE dimensions
  if [ "$PARALLEL" == "MPI" -o "$PARALLEL" == "OPENMP" ]; then
 #   # caution here : respect the case for LLm0 and MMm0 in param.h, assumed to be on the same line etc
 #   old stuff
 #   NX=$(awk  "/(.*)$EXAMPLE/ {c=1; next} c-->0  " param.h  |sed -e 's/ *//g' -e 's/.*LLm0=\(.*\),\(.*\)MMm0.*/\1/' )
 #   NY=$(awk  "/(.*)$EXAMPLE/ {c=1; next} c-->0  " param.h  |sed -e 's/ *//g' -e 's/.*MMm0*=\(.*\),\(.*\)N.*/\1/' )
    NX=$(cpp -P -D$EXAMPLE -I${ROOTDIR}/OCEAN -include cppdefs.h param.h |grep LLm0 |sed -e '/^! */d' -e '/^ *!/d' -e 's/ *//g' -e 's/.*LLm0=\(.*\),\(.*\)MMm0.*/\1/' | sed  -n 2p )
    NY=$(cpp -P -D$EXAMPLE -I${ROOTDIR}/OCEAN -include cppdefs.h param.h |grep MMm0 |sed -e '/^! */d' -e '/^ *!/d' -e 's/ *//g' -e 's/.*MMm0=\(.*\),\(.*\)N.*/\1/' | sed  -n 2p   )

    PROC_X=1
    PROC_Y=1

#-  find maximum possible cpu in X
    while [ $( echo $NX/$PROC_X | bc ) -gt 16  -a $PROC_X -lt $MAX_PROC ];
    do
      ((PROC_X=$PROC_X+1))
    done
    echo $MAX_PROC $PROC_X
    ((MAX_PROC_Y=$MAX_PROC - $PROC_X ))
#-  find maximum possible cpu in Y
    while [ $( echo $NY/$PROC_Y | bc ) -gt 16  -a $MAX_PROC_Y -ge $PROC_X ];
    do
    ((PROC_Y=$PROC_Y+1))
    ((MAX_PROC_Y=$MAX_PROC - $PROC_X*$PROC_Y ))
    done
    ((NB_PROC=$PROC_Y * $PROC_X ))

    echo
    echo "Number of cpu (TOT/X/Y):"
    echo $NB_PROC $PROC_X $PROC_Y
    echo

    if [ "$PARALLEL" == "MPI" ] ; then
      NEW_LINE="parameter (NP_XI=$PROC_X,NP_ETA=$PROC_Y,NNODES=NP_XI*NP_ETA)"
      sed  -e /"^ *parameter\(.*\)NP_XI*="'/c\'$'\n'"\      ${NEW_LINE}" param.h > tmp.txt && \mv tmp.txt param.h
    else
      NEW_LINE="parameter (NSUB_X=$PROC_X,NSUB_E=$PROC_Y)"
      sed  -e /"^\ *parameter\(.*\)NSUB_X*="'/c\'$'\n'"\      ${NEW_LINE}" param.h > tmp.txt && \mv tmp.txt param.h
      NEW_LINE="parameter (NPP=$NB_PROC)"
      sed  -e /"^\ *parameter\(.*\)NPP*="'/c\'$'\n'"\      ${NEW_LINE}" param.h > tmp.txt && \mv tmp.txt param.h
    fi
  fi

#- input if needed
#  if [ "$EXAMPLE" == "VORTEX" -o "$EXAMPLE" == "JET" ]; then
  if [ "$EXAMPLE" == "VORTEX" ]; then
    example=$(echo $EXAMPLE |tr '[:upper:]' '[:lower:]')
    if [ ! -f  ${example}_grd.nc ]; then
      mymake="make_${example}"
      matlab -nodesktop  -nosplash -nodisplay -r "addpath ./TEST_CASES; ${mymake};exit"
    fi
  fi

#- compile and run
  ./jobcomp || exit 1
  if [ "${PARALLEL}" == "MPI" ]; then
    $MPIRUN -np $NB_PROC ./croco || exit 1
  elif [ "${PARALLEL}" == "OPENMP" ]; then
    export OMP_NUM_THREADS=$NB_PROC ; ./croco || exit 1
  else
    ./croco || exit 1
  fi

#- revert stuff  for safety
  sed  -e "s/#define $EXAMPLE /#undef  $EXAMPLE/" cppdefs.h   > tmp.txt && \mv tmp.txt cppdefs.h
  NEW_LINE="parameter (NP_XI=1,NP_ETA=1_Y,NNODES=NP_XI*NP_ETA)"
  sed  -e /"^ *parameter\(.*\)NP_XI*="'/c\'$'\n'"\      ${NEW_LINE}" param.h > tmp.txt && \mv tmp.txt param.h
  NEW_LINE="parameter (NSUB_X=1,NSUB_E=1)"
  sed  -e /"^\ *parameter\(.*\)NSUB_X*="'/c\'$'\n'"\      ${NEW_LINE}" param.h > tmp.txt && \mv tmp.txt param.h
  NEW_LINE="parameter (NPP=1)"
  sed  -e /"^\ *parameter\(.*\)NPP*="'/c\'$'\n'"\      ${NEW_LINE}" param.h >  tmp.txt && \mv tmp.txt param.h
done

# qq sed utiles si besoin :

# selectionne jusque non inclus
#sed -n '1,/^#if defined REGIONAL/ {/^#if defined REGIONAL/!p;} ' cppdefs.h

# selectionne jusque non inclus et undef
# sed -ne '1,/^#if defined REGIONAL/ s/define /undef /g ;  1,/^#if defined REGIONAL/ {/^#if defined REGIONAL/!p;} ; /\!/d' cppdefs.h

# extrait toutes les configs
#sed  -e  '1,/^#if defined REGIONAL/ s/define /undef /g  ; /^#if defined REGIONAL/,$d ; /^$/d; /^!/d; /^\//d; /^\*/d; /^\ /d' cppdefs.h | awk -F' ' '{print $2}'
