export USE_TOYATM=0
export USE_TOYOCE=0
export USE_TOYWAV=0
#
[ ${USE_TOYATM}  -eq 1 ] && istoy=".toyatm"  || istoy=""
[ ${USE_TOYWAV}  -eq 1 ] && istoy="${istoy}.toywav"  || istoy="${istoy}"
[ ${USE_TOYOCE}  -eq 1 ] && istoy="${istoy}.toyoce"  || istoy="${istoy}"
#
