#!/bin/bash
#
########################################################
#  Define files and run parameters
########################################################
#
# Name used for the input files. For example croco_grd.nc
MODEL=croco

# Scratch directory where the model is run
SCRATCHDIR=`pwd`/SCRATCH

# Input directory where the croco_inter.in input file is located
INPUTDIR=`pwd`/CROCO_IN  # prod architecture
#INPUTDIR=`pwd`          # dev architecture

# AGRIF input file which defines the position of child grids
AGRIF_FILE=AGRIF_FixedGrids.in

# Directory where the croco input NetCDF files (croco_grd.nc, ...) are stored
MSSDIR=`pwd`/CROCO_FILES

# Directory where the croco output and restart NetCDF files (croco_his.nc, ...) are stored
MSSOUT=$SCRATCHDIR

# CROCO executable
CODFILE=croco

# number of processors for MPI run
NBPROCS=8

# command for running the mode : ./ for sequential job, mpirun -np NBPROCS for mpi run
# WARNING: for mpi run command, it is needed to add a space at the end!
RUNCMD='./'
#RUNCMD="mpirun -np $NBPROCS "
#RUNCMD="$MPI_LAUNCH "
#RUNCMD='srun '

#  Define environment variables for OPENMP
OMP_SCHEDULE=static
OMP_NUM_THREADS=2
OMP_DYNAMIC=false
OMP_NESTED=false
KMP_LIBRARY=throughput
KMP_STACKSIZE=2m
KMP_DUPLICATE_LIB_OK=TRUE

# Define which type of inputs are used
BULK_FILES=1
FORCING_FILES=0
CLIMATOLOGY_FILES=0
BOUNDARY_FILES=1
RUNOFF_FILES=0

# Atmospheric surface forcing dataset used for the bulk formula (NCEP)
ATMOS_BULK=ERA5
# Atmospheric surface forcing dataset used for the wind stress (NCEP, QSCAT)
ATMOS_FRC=QSCAT
# Oceanic boundary and initial dataset (SODA, ECCO,...)
OGCM=SODA
# Runoff dataset (Daie and Trenberth,...)
RUNOFF_DAT=DAI

# Model time step [seconds]
DT=3600
# Number of barotropic time steps within one baroclinic time step [number], NDTFAST in croco.in
NFAST=60

# Number total of grid levels (1: No child grid)
NLEVEL=1
# AGRIF nesting refinement coefficient
AGRIF_REF=3

# Start and End year 
NY_START=2005
NY_END=2005
# Start and End month
NM_START=1
NM_END=3
# Set month format at 1 or 2 digits (for input and output files): "%01d" = 1 digit/ "%02d" = 2 digit  
MTH_FORMAT="%02d"
#  Time Schedule  -  TIME_SCHED=0 --> yearly files
#                    TIME_SCHED=1 --> monthly files
TIME_SCHED=1

# Number of year that are considered to be part of the spin-up (i.e. 365 days per year)
NY_SPIN=0

# Output frequency [days]
#   average
ND_AVG=3
#   history (if = -1 set equal to NUMTIMES)
ND_HIS=-1
#   restart (if = -1 set equal to NUMTIMES)
ND_RST=-1

#  Restart file - RSTFLAG=0 --> No Restart
#		  RSTFLAG=1 --> Restart
RSTFLAG=0     
#  Exact restart - EXACT_RST=0 --> Exact restart OFF
#                - EXACT_RST=1 --> Exact restart ON
EXACT_RST=0

#unalias cp
#unalias mv
#limit coredumpsize unlimited
CP=/bin/cp
MV=/bin/mv
LN=/bin/ln
#
########################################################
#  END USER CHANGE
########################################################
#
if [[ $TIME_SCHED == 0 ]]; then
  NM_START=1979
  NM_END=1979
fi
#
# netcdf file prefixes
#
GRDFILE=${MODEL}_grd
FRCFILE=${MODEL}_frc
BLKFILE=${MODEL}_blk
INIFILE=${MODEL}_ini
CLMFILE=${MODEL}_clm
BRYFILE=${MODEL}_bry
RNFFILE=${MODEL}_runoff
#
if [ ! -e $MSSOUT ] ; then
 mkdir $MSSOUT
fi
#
if [[ $RSTFLAG != 0 ]]; then
  NY=$NY_START
  NM=$NM_START
  if [[ $TIME_SCHED == 0 ]]; then
    NY=$((NY - 1))
    TIME=Y${NY}
  else
    NM=$((NM - 1))
    if [[ $NM == 0 ]]; then
      NM=12
      NY=$((NY - 1))
    fi
    TIME=Y${NY}M$( printf ${MTH_FORMAT} ${NM})
  fi
  RSTFILE=${MODEL}_rst_${TIME}
fi
#
if [[ $TIME_SCHED == 0 ]]; then
  TIME=Y${NY_START}
else
  TIME=Y${NY_START}M$( printf ${MTH_FORMAT} ${NM_START})
fi
#
# Get the code
#
if [ ! -e $SCRATCHDIR ] ; then
 mkdir $SCRATCHDIR
fi
cd $SCRATCHDIR
echo "Getting $CODFILE from $INPUTDIR"
$CP -f $INPUTDIR/$CODFILE $SCRATCHDIR
chmod u+x $CODFILE
echo "Getting $AGRIF_FILE from $INPUTDIR"
$CP -f $INPUTDIR/$AGRIF_FILE $SCRATCHDIR
#
# Get the netcdf files
#
LEVEL=0
while [ $LEVEL != $NLEVEL ]; do
  if [[ ${LEVEL} == 0 ]]; then
    ENDF=
  else
    ENDF=.${LEVEL}
  fi
  echo "Getting ${GRDFILE}.nc${ENDF} from $MSSDIR"
  $LN -sf $MSSDIR/${GRDFILE}.nc${ENDF} $SCRATCHDIR
  echo "Getting ${MODEL}_inter.in${ENDF} from $INPUTDIR"
  $CP -f $INPUTDIR/${MODEL}_inter.in${ENDF} $SCRATCHDIR
  if [[ $RSTFLAG == 0 ]]; then
    echo "Getting ${INIFILE}_${OGCM}_${TIME}.nc${ENDF} from $MSSDIR"
    $CP -f $MSSDIR/${INIFILE}_${OGCM}_${TIME}.nc${ENDF} $SCRATCHDIR
    $CP -f ${INIFILE}_${OGCM}_${TIME}.nc${ENDF} ${INIFILE}.nc${ENDF}
  else
    echo "Getting ${RSTFILE}.nc${ENDF} from $MSSOUT"
    $CP -f $MSSOUT/${RSTFILE}.nc${ENDF} $SCRATCHDIR
    $CP -f ${RSTFILE}.nc${ENDF} ${INIFILE}.nc${ENDF}
  fi
  LEVEL=$((LEVEL + 1))
done
###########################################################
#  Compute
###########################################################
#
NY_END=$((NY_END + 1))
NM_END=$((NM_END + 1))
NY=$NY_START
while [ $NY != $NY_END ]; do
  if [[ $NY == $NY_START ]]; then
    NM=$NM_START
  else 
     NM=1
  fi
   MY_YEAR=$NY
   MY_YEAR=$((MY_YEAR + 1))
  if [[ $MY_YEAR == $NY_END ]]; then
     MONTH_END=$NM_END
  else 
     MONTH_END=13
  fi
  if [[ $TIME_SCHED == 0 ]]; then
     MONTH_END=2
  fi
  while [ $NM != $MONTH_END ]; do
    if [[ $TIME_SCHED == 0 ]]; then
      TIME=Y${NY}
      echo "Computing YEAR $NY"
    else
	TIME=Y${NY}M$( printf ${MTH_FORMAT} ${NM})
	echo "Computing YEAR $NY MONTH $( printf ${MTH_FORMAT} ${NM})"
    fi
#
# Get forcing and clim for this time
#
    LEVEL=0
    while [ $LEVEL != $NLEVEL ]; do
      if [[ ${LEVEL} == 0 ]]; then
        ENDF=
      else
        ENDF=.${LEVEL}
      fi
      if [[ ${FORCING_FILES} == 1 ]]; then
        echo "Getting ${FRCFILE}_${ATMOS_FRC}_${TIME}.nc${ENDF} from $MSSDIR"
        $LN -sf $MSSDIR/${FRCFILE}_${ATMOS_FRC}_${TIME}.nc${ENDF} ${FRCFILE}.nc${ENDF}
      fi
      if [[ ${BULK_FILES} == 1 ]]; then
        echo "Getting ${BLKFILE}_${ATMOS_BULK}_${TIME}.nc${ENDF} from $MSSDIR"
        $LN -sf $MSSDIR/${BLKFILE}_${ATMOS_BULK}_${TIME}.nc${ENDF} ${BLKFILE}.nc${ENDF}
      fi
     if [[ ${RNF_FILES} == 1 ]]; then
        echo "Getting ${RNFFILE}_${RUNOFF_DAT}_${TIME}.nc${ENDF} from $MSSDIR"
        $LN -sf $MSSDIR/${RNFFILE}_${RUNOFF_DAT}_${TIME}.nc${ENDF} ${RNFFILE}.nc${ENDF}
      fi
      
      LEVEL=$((LEVEL + 1))
    done
#
# No child climatology or boundary files
#
    if [[ ${CLIMATOLOGY_FILES} == 1 ]]; then
      echo "Getting ${CLMFILE}_${OGCM}_${TIME}.nc from $MSSDIR"
      $LN -sf $MSSDIR/${CLMFILE}_${OGCM}_${TIME}.nc ${CLMFILE}.nc
    fi
    if [[ ${BOUNDARY_FILES} == 1 ]]; then
      echo "Getting ${BRYFILE}_${OGCM}_${TIME}.nc from $MSSDIR"
      $LN -sf $MSSDIR/${BRYFILE}_${OGCM}_${TIME}.nc ${BRYFILE}.nc
    fi
#
# Set the number of time steps for each month 
# (30 or 31 days + 28 or 29 days for february)
#
    NUMTIMES=0
#
    if [[ ${NM} == 1 || ${NM} == 3 || ${NM} == 5 || ${NM} == 7 || ${NM} == 8 || ${NM} == 10 || ${NM} == 12 ]]; then
      NDAYS=31
    else
      NDAYS=30
      if [[ ${NM} == 2 ]]; then
        NDAYS=28
# February... check if it is a leap year

        B4=0
        B100=0
        B400=0

        B4=$((4 * ( NY / 4 )))
        B100=$((100 * ( NY / 100 )))
        B400=$((400 * ( NY / 400 )))

	
        if [[ $NY == $B4 && ((!($NY == $B100))||($NY == $B400)) ]]; then
#
          BSPIN=$(( NY - NY_START + 1 ))
          if [[ $BSPIN -gt $NY_SPIN ]]; then
	     echo Leap Year - $NY $B4 $B100 $B400
             NDAYS=29
          else
#.........   SPINUP!!!! In case of spinup I cant have leap years.
	     echo year $NY should be a Leap Year     
	     echo 'BUT : Spinup case: no leap year'
             NDAYS=28
          fi
#
        else
	  echo Not a Leap Year - $NY $B4 $B100 $B400
          NDAYS=28	  		  
        fi
      fi
    fi
    #
    # Put the number of time steps in the .in files
    #
    echo "YEAR = $NY MONTH = $NM DAYS = $NDAYS DT = $DT NTIMES = $NUMTIMES"
    NUMTIMES=$((NDAYS * 24 * 3600))
    NUMTIMES=$((NUMTIMES / DT))

    DT0=$DT
    LEVEL=0
    while [[ $LEVEL != $NLEVEL ]]; do
	if [[ ${LEVEL} == 0 ]]; then
            ENDF=
	else
            ENDF=.${LEVEL}
	    NUMTIMES=$((AGRIF_REF * NUMTIMES))
	    DT=$((DT / AGRIF_REF))
	fi
	NUMAVG=$((ND_AVG * 86400 / DT ))
	if [[ ${ND_HIS} -ne -1 ]]; then
	    NUMHIS=$((ND_HIS * 86400 / DT ))
	else
	    NUMHIS=$NUMTIMES
	fi
	if [[ ${ND_RST} -ne -1 ]]; then
	    NUMRST=$((ND_RST * 86400 / DT ))
	else
	    NUMRST=$NUMTIMES
	fi

	if [[ $EXACT_RST == 1 ]]; then
	    echo "Exact restart defined"
	    if [[ $NY == $NY_START && $NM == $NM_START ]]; then
		NUMRECINI=1
		echo "set NUMRECINI = $NUMRECINI"
	    else
		NUMRECINI=2
		echo "set NUMRECINI = $NUMRECINI"
	    fi
	else  # no exact restart
	    echo "No exact restart"
	    NUMRECINI=1
	    echo "set NUMRECINI = $NUMRECINI"
	fi
	
	echo " "
	echo "Writing in ${MODEL}_inter.in${ENDF}"
	echo "USING DT       = $DT"
	echo "USING NFAST    = $NFAST"
	echo "USING NUMTIMES = $NUMTIMES"
	echo "USING NUMAVG   = $NUMAVG"
	echo "USING NUMHIS   = $NUMHIS"
	echo "USING NUMRST   = $NUMRST"
	echo "USING NUMRECINI = $NUMRECINI"
	
	if [ ! -f ${MODEL}_inter.in${ENDF} ]; then
	    echo "=="
	    echo "=> ERROR : miss the ${MODEL}_inter.in${ENDF} file"
	  echo "=="
	  exit 1
	fi
	sed -e 's/NUMTIMES/'$NUMTIMES'/' -e 's/TIMESTEP/'$DT'/' -e 's/NFAST/'$NFAST'/' \
	    -e 's/\bNUMAVG\b/'$NUMAVG'/' -e 's/\bNUMHIS\b/'$NUMHIS'/' -e 's/\bNUMRST\b/'$NUMRST'/' \
	    -e 's/NUMRECINI/'$NUMRECINI'/' \
	    -e 's/NYONLINE/'$NY'/' -e 's/NMONLINE/'$NM'/' < ${MODEL}_inter.in${ENDF} > ${MODEL}_${TIME}_inter.in${ENDF}
	
	LEVEL=$((LEVEL + 1))
    done
    DT=$DT0
    #
    #  COMPUTE
    #
    echo " "
    echo "Computing for $TIME"
    date
    ${RUNCMD}$CODFILE  ${MODEL}_${TIME}_inter.in > ${MODEL}_${TIME}.out
    date
    #
    
    # Test if the month has finised properly
    echo "Test ${MODEL}_${TIME}.out"
    status=`tail -2 ${MODEL}_${TIME}.out | grep DONE | wc -l`
    if [[ $status == 1 ]]; then
      echo "All good"
      echo "XXXX${MYTEST}XXXX"
    else
      echo
      echo "Warning: month not finished properly"
      echo
      tail -20 ${MODEL}_${TIME}.out
      echo
      echo "Month ${TIME} did not work"
      echo
      exit 1
    fi
#
#  Archive
#
    LEVEL=0
    while [[ $LEVEL != $NLEVEL ]]; do
      if [[ ${LEVEL} == 0 ]]; then
        ENDF=
      else
        ENDF=.${LEVEL}
      fi
	  $CP -f ${MODEL}_rst.nc${ENDF} ${INIFILE}.nc${ENDF}
	  $MV -f ${MODEL}_his.nc${ENDF} ${MSSOUT}/${MODEL}_his_${TIME}.nc${ENDF}
	  $MV -f ${MODEL}_rst.nc${ENDF} ${MSSOUT}/${MODEL}_rst_${TIME}.nc${ENDF}
	  $MV -f ${MODEL}_avg.nc${ENDF} ${MSSOUT}/${MODEL}_avg_${TIME}.nc${ENDF}
      LEVEL=$((LEVEL + 1))
    done
    NM=$((NM + 1))
  done
  NY=$((NY + 1))
done
#
#############################################################












