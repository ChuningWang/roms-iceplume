! This is include file "diags_ek.h"
!  ==== == ======= ==== ==========
!

#if defined OUTPUTS_SURFACE && defined AVERAGES && ! defined XIOS

      real timesurf_avg
      real surft_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE surft(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real surfs_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE surft(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real surfz_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE surft(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real surfu_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE surft(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real surfv_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE surft(BLOCK_PATTERN,*) BLOCK_CLAUSE

      common /surf_timesurf_avg/timesurf_avg
      common /surft_avg/surft_avg
     &       /surfs_avg/surfs_avg
     &       /surfz_avg/surfz_avg
     &       /surfu_avg/surfu_avg
     &       /surfv_avg/surfv_avg

# endif



 
