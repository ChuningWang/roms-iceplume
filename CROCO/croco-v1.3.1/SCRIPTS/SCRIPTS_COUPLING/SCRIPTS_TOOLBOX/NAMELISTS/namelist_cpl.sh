#-------------------------------------------------------------------------------
# CPL
#-------------------------------------------------------------------------------
# namelist
export namcouplename=namcouple.base.${RUNtype} 
# Note: namelist example files are provided in OASIS_IN/
# if you want to use a pre-built weight file for grid interpolations, point to 
# e.g. namcouple.base.oa.smtho2a 

# coupling frequency
export CPL_FREQ=3600

# Start from existing condition
export CPL_restart="FALSE" # If TRUE: initialize CPL field with an history field
 
# Files to create oasis restart (leave empty if none). File need to be in *FILES_DIR. Default 1st time step
export oce_rst_file=""   
export atm_rst_file=""
export wav_rst_file=""

# Use pre-built weight file for grid interpolations
# The weight file should be placed in OASIS_FILES
# Note: for the moment only works with 1-atm and 1-oce domain
export WEIGHT_FLAG=0 
export weight_o2a=""
export weight_a2o=""

