! This is include file "diags_ek.h"
!  ==== == ======= ==== ==========
!

# if defined DIAGNOSTICS_EDDY && ! defined XIOS

# ifdef AVERAGES
      real timediags_eddy_avg
      real eddyzz_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE eddyzz_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyuu_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddyuu_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyvv_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddyvv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyuv_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddyuv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyub_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddyub_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyvb_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddyvb_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddywb_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddywb_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyuw_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddyuw_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyvw_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddyvw_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyubu_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE eddyubu_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyvbv_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE eddyvbv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyusu_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE eddyusu_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyvsv_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE eddyvsv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyugsu_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE eddyugsu_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyvgsv_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE eddyvgsv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif

# ifdef AVERAGES
      common /diag_timediags_eddy_avg/timediags_eddy_avg
      common /diag_eddyzz_avg/eddyzz_avg
     &       /diag_eddyuu_avg/eddyuu_avg
     &       /diag_eddyvv_avg/eddyvv_avg
     &       /diag_eddyuv_avg/eddyuv_avg
     &       /diag_eddyub_avg/eddyub_avg
     &       /diag_eddyvb_avg/eddyvb_avg
     &       /diag_eddywb_avg/eddywb_avg
     &       /diag_eddyuw_avg/eddyuw_avg
     &       /diag_eddyvw_avg/eddyvw_avg
     &       /diag_eddyubu_avg/eddyubu_avg
     &       /diag_eddyvbv_avg/eddyvbv_avg
     &       /diag_eddyusu_avg/eddyusu_avg
     &       /diag_eddyvsv_avg/eddyvsv_avg
     &       /diag_eddyugsu_avg/eddyugsu_avg
     &       /diag_eddyvgsv_avg/eddyvgsv_avg   
# endif
#endif /* DIAGNOSTICS_EDDY*/
