#-------------------------------------------------------------------------------
# WAV
#-------------------------------------------------------------------------------
# namelist

# Time steps
export DT_WAV=3600     # TMAX = 3*TCFL
export DT_WW_PRO=1200  # TCFL = 0.8 x dx/(g/fmin4pi) with fmin=0.0373 => 3-4 % of dx
export DT_WW_REF=1800  # TMAX / 2
export DT_WW_SRC=10    # TSRC = usually 10s  (could be between 5s and 60s)

# Grid size
export wavnx=41 ; export wavny=42

# forcing files
export forcin=() # forcing file(s) PREFIX list (leave empty if none), input filenames are supposed to be in the form: PREFIX_Y????M??.nc
export forcww3=() # name of ww3_prnc.inp extension/input file

#boundary files
export bouncin= # prefix for boundary files (leave empty is none)

# output settings
export flagout="TRUE" # Keep (TRUE) or not (FALSE) ww3 full output binary file (out_grd.ww3)
export wav_int=21600            # output interval (s)
# ww3 file to be used for creating restart file for oasis 
export wavfile=$CWORK/rundir/BENGUELA_LR_exp1_Wfrc_outputs/20050101_20050131/ww3_20050101_20050131.nc # Usually done by running a frc mode on the area

