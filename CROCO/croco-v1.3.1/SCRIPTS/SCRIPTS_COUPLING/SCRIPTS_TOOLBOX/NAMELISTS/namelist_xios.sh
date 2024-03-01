#-------------------------------------------------------------------------------
# XIOS
#-------------------------------------------------------------------------------

export ONLINE_XML="FALSE" # Process xml online

# name of xml files (defined in file_def_*). Here are default names in xml files.
export ATM_XIOS_NAME="wrfout wrf3d_1D wrf3d_1H" # All the names you have defined in your xml file
export OCE_OUTPUT_PREFIX="croco"
export OCE_XIOS_NAME="${OCE_OUTPUT_PREFIX}_3h_inst ${OCE_OUTPUT_PREFIX}_1h_avg_3d ${OCE_OUTPUT_PREFIX}_1h_inst_surf ${OCE_OUTPUT_PREFIX}_5d_aver" # Relative to what is in the file_def_croco.xml

