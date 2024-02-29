############################# from mynamelist.sh ###############################

################################################################################
############################ USER CHANGES ######################################
################################################################################
#
export CEXPER=BENGUELA
export RUNtype=owa
#
export USE_ATM=1
export USE_TOYATM=0
export USE_XIOS_ATM=0
#
export USE_OCE=1
export USE_TOYOCE=0
export USE_XIOS_OCE=0
#
export USE_WAV=1
export USE_TOYWAV=0
#
[ ${USE_TOYATM}  -eq 1 ] && istoy=.toyatm  || istoy=''
[ ${USE_TOYWAV}  -eq 1 ] && istoy=.toywav  || istoy=''
[ ${USE_TOYOCE}  -eq 1 ] && istoy=.toyoce  || istoy=''
#
#-------------------------------------------------------------------------------
# Exe paths
# ------------------------------------------------------------------------------
export ATM_EXE_DIR=$ATM/exe_coupled
export OCE_EXE_DIR=$CHOME/croco_in
export WAV_EXE_DIR=$WAV/exe_${RUNtype}_${CEXPER}
export TOY_EXE_DIR=$CHOME/toy_in
export XIOS_EXE_DIR=$XIOS

#-------------------------------------------------------------------------------
# RUN_DIR
#-------------------------------------------------------------------------------
export EXEDIR_ROOT="$CWORK/rundir/${CEXPER}_execute"
export OUTPUTDIR_ROOT="$CWORK/rundir/${CEXPER}_outputs"
export RESTDIR_ROOT="$CWORK/rundir/${CEXPER}_restarts"

export  JOBDIR_ROOT=${CHOME}/jobs_${CEXPER}

#-------------------------------------------------------------------------------
# Model settings
# ------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# CPL
#-------------------------------------------------------------------------------
# namelist
export namcouplename=namcouple.base.${RUNtype}${istoy}

# coupling frequency
export CPL_FREQ=21600

#-------------------------------------------------------------------------------
# OCE
#-------------------------------------------------------------------------------
# namelist

# Online Compilation
export ONLINE_COMP=1

# Time steps
export TSP_OCE=1200
export TSP_OCEF=60

# Grid sizes
export ocenx=41 ; export oceny=42
export hmin=75; # minimum water depth in CROCO, delimiting coastline in WW3 

# domains
export AGRIFZ=0
export AGRIF_NAMES=""

# forcing files
export ini_ext='ini_SODA' # ini extension file (ini_SODA,...)
export bry_ext='bry_SODA' # bry extension file (bry_SODA,...)
export surfrc_flag="FALSE" # Flag if surface forcing is needed (FALSE if cpl)
export interponline=0 # switch (1=on, 0=off) for online surface interpolation
export frc_ext='blk_CFSR' # surface forcing extension(blk_CFSR, frc_CFSR,...). If interponline=1 just precise the type (ECMWF, CFSR,AROME,...)
export tide_flag="FALSE" # the forcing extension must be blk_??? otherwise tide forcing overwrites it 

# output settings
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                                          WARNING                                       ! 
# When XIOS is activated the following values (for the model) are not taken into account !
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
export oce_nhis=18     # history output interval (in number of timesteps) 
export oce_navg=18     # average output interval (in number of timesteps) 

#-------------------------------------------------------------------------------
# ATM
#-------------------------------------------------------------------------------
# namelist
export atmnamelist=namelist.input.prep.${CEXPER}.${RUNtype}

# Time steps
export TSP_ATM=100 #100   # 100 90 75 72 60 45

# Grid size
export atmnx=56 ; export atmny=50

# domains
export NB_dom=1 # Number of coupled domains
export wrfcpldom='d01'

# output settings
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                                          WARNING                                       ! 
# When XIOS is activated the following values (for the model) are not taken into account !
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
export atm_his_h=6                        # output interval (h)
export atm_his_frames=1000 # $((31*24))          # nb of outputs per file
export atm_diag_int_m=$((${atm_his_h}*60))  # diag output interval (m)
export atm_diag_frames=1000     # nb of diag outputs per file

#-------------------------------------------------------------------------------
# WAV
#-------------------------------------------------------------------------------
# namelist

# Time steps
export TSP_WAV=3600     # TMAX = 3*TCFL
export TSP_WW_PRO=1200  # TCFL --> ww3.grid to see the definition
export TSP_WW_REF=1800  # TMAX / 2
export TSP_WW_SRC=10

# Grid size
export wavnx=41 ; export wavny=42

# forcing files
export forcin=() # forcing file(s) list (leave empty if none)
export forcww3=() # name of ww3_prnc.inp extension/input file

# output settings
export flagout="TRUE" # Keep (TRUE) or not (FALSE) ww3 full output binary file (out_grd.ww3)
export wav_int=21600            # output interval (s)
# ww3 file to be used for creating restart file for oasis 
export wavfile=$CWORK/outputs_frc_ww3_CFSR/ww3.200501.nc # Usually done by running a frc mode on the area

#-------------------------------------------------------------------------------
# TOY
#-------------------------------------------------------------------------------
# type
export toytype=("wav" "atm") #oce,atm,wav

# namelist

# Time steps
### TOY ###
#!!!! WARNING !!!!
# When using toy, please put the output frequency of "toyfile" in the TSP of the model
# example: ww3.200501.nc output_frequency=3600 -----> TSP_WW3=3600
#!!!!!!!!!!!!!!!!!

# forcing files
export toyfile=("$CWORK/toy_files/ww3.frc.200501.nc" "$CWORK/toy_files/wrfout_d01_20050101_20050131.nc")
export timerange=('1,745' '2,125')

#-------------------------------------------------------------------------------
# XIOS
#-------------------------------------------------------------------------------
# namelist
export FILIN_XIOS="iodef.xml context_roms.xml context_wrf.xml field_def.xml domain_def.xml file_def_croco.xml file_def_wrf.xml" # files needed for xios. Need to be in INPUTDIRX (cf header.*)

# files
export ATM_XIOS_NAME="wrfout_inst" # All the names you have defined in your xml file
export OCE_XIOS_NAME="${CEXPER}_1d_aver ${CEXPER}_1h_inst_surf ${CEXPER}_6h_his"

################################################################################
############################ END USER CHANGE ###################################
################################################################################

# KEY for XIOS and TOY #
export USE_XIOS=$(( ${USE_XIOS_ATM} + ${USE_XIOS_OCE} ))
export USE_TOY=$(( ${USE_TOYATM} + ${USE_TOYOCE} + ${USE_TOYWAV} ))
[ ${USE_TOY} -ge 1 ] && export USE_CPL=1 || export USE_CPL=$(( ${USE_ATM} * ${USE_OCE} + ${USE_ATM} * ${USE_WAV} + ${USE_OCE} * ${USE_WAV} ))

### TOY ###
export nbtoy=${#toytype[@]}
export model_to_toy=()
export toynamelist=()
export TSP_TOY=()
for k in `seq 0 $(( ${nbtoy} - 1))` ; do
    [ ${toytype[$k]} == "oce" ] && model_to_toy+=("croco")
    [ ${toytype[$k]} == "atm" ] && model_to_toy+=("wrf")
    [ ${toytype[$k]} == "wav" ] && model_to_toy+=("ww3")
    if [ ${nbtoy} -eq 1 ]; then
        toynamelist+=("TOYNAMELIST.nam.${toytype[$k]}.${RUNtype}")
    elif [ ${nbtoy} -eq 2 ]; then
        toynamelist+=("TWOTOYNAMELIST.nam.${toytype[$k]}.${RUNtype}")
    else
        printf "\n No namelist available for three namelist toy\n" ; exit 0
    fi
    [ ${toytype[$k]} == "oce" ] && TSP_TOY+=("${TSP_OCE}")
    [ ${toytype[$k]} == "atm" ] && TSP_TOY+=("${TSP_ATM}")
    [ ${toytype[$k]} == "wav" ] && TSP_TOY+=("${TSP_WAV}")
done
######


