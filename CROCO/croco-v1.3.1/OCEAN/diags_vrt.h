! This is include file "diags_vrt.h"
!  ==== == ======= ==== ==========
!

#ifdef DIAGNOSTICS_VRT

      real vrtXadv(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtXadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtYadv(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtYadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtHdiff(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtHdiff(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtCor(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtCor(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtPrsgrd(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtPrsgrd(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtHmix(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtHmix(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtVmix(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtVmix(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtrate(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtrate(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtVmix2(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtVmix2(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtWind(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtWind(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtDrag(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtDrag(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real wrkWind(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE wrkWind(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real wrkDrag(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE wrkDrag(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined DIAGNOSTICS_BARO
      real vrtBaro(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtBaro(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
# if defined M3FAST
      real vrtfast(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtfast(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
# ifdef AVERAGES
      real timediags_vrt_avg
      real vrtXadv_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtXadv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtYadv_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtYadv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtHdiff_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtHdiff_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtCor_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtCor_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtPrsgrd_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtPrsgrd_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtHmix_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtHmix_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtVmix_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtVmix_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtrate_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtrate_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtVmix2_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtVmix2_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtWind_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtWind_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrtDrag_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtDrag_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined DIAGNOSTICS_BARO
      real vrtBaro_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtBaro_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
# if defined M3FAST
      real vrtfast_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtfast_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
# endif



      common /diag_vrtXadv/vrtXadv   
     &       /diag_vrtYadv/vrtYadv
     &       /diag_vrtHdiff/vrtHdiff
     &       /diag_vrtCor/vrtCor
     &       /diag_vrtPrsgrd/vrtPrsgrd
     &       /diag_vrtHmix/vrtHmix
     &       /diag_vrtVmix/vrtVmix
     &       /diag_vrtrate/vrtrate
     &       /diag_vrtVmix2/vrtVmix2
     &       /diag_vrtWind/vrtWind
     &       /diag_vrtDrag/vrtDrag
     &       /diag_wrkWind/wrkWind
     &       /diag_wrkDrag/wrkDrag
# if defined DIAGNOSTICS_BARO
     &       /diag_vrtBaro/vrtBaro
# endif
# if defined M3FAST
     &       /diag_vrtfast/vrtfast
# endif

# ifdef AVERAGES
      common /diag_timediags_vrt_avg/timediags_vrt_avg
      common /diag_vrtXadv_avg/vrtXadv_avg
     &       /diag_vrtYadv_avg/vrtYadv_avg
     &       /diag_vrtHdiff_avg/vrtHdiff_avg
     &       /diag_vrtCor_avg/vrtCor_avg
     &       /diag_vrtPrsgrd_avg/vrtPrsgrd_avg
     &       /diag_vrtHmix_avg/vrtHmix_avg
     &       /diag_vrtVmix_avg/vrtVmix_avg
     &       /diag_vrtrate_avg/vrtrate_avg
     &       /diag_vrtVmix2_avg/vrtVmix2_avg
     &       /diag_vrtWind_avg/vrtWind_avg
     &       /diag_vrtDrag_avg/vrtDrag_avg
# if defined DIAGNOSTICS_BARO
     &       /diag_vrtBaro_avg/vrtBaro_avg
# endif
# if defined M3FAST
     &       /diag_vrtfast_avg/vrtfast_avg
# endif
# endif      


#ifndef DIAGNOSTICS_UV
      real wrkXadv(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE wrkXadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real wrkYadv(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE wrkYadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real wrkHdiff(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE wrkHdiff(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real wrkCor(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE wrkCor(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real wrkPrsgrd(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE wrkPrsgrd(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real wrkHmix(GLOBAL_2D_ARRAY,2,2)
!CSDISTRIBUTE_RESHAPE wrkHmix(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real wrkVmix(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE wrkVmix(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real wrkrate(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE wrkrate(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real wrkVmix2(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE wrkVmix2(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined DIAGNOSTICS_BARO
      real wrkBaro(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE wrkBaro(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
# if defined M3FAST
      real wrkfast(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE wrkfast(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
      common /diag_wrkXadv/wrkXadv   
     &       /diag_wrkYadv/wrkYadv
     &       /diag_wrkHdiff/wrkHdiff
     &       /diag_wrkCor/wrkCor
     &       /diag_wrkPrsgrd/wrkPrsgrd
     &       /diag_wrkHmix/wrkHmix
     &       /diag_wrkVmix/wrkVmix
     &       /diag_wrkrate/wrkrate
     &       /diag_wrkVmix2/wrkVmix2
# if defined DIAGNOSTICS_BARO
     &       /diag_wrkBaro/wrkBaro
# endif
# if defined M3FAST
     &       /diag_wrkfast/wrkfast
# endif
#endif /*ifndef DIAGNOSTICS_UV */

#endif /* DIAGNOSTICS_VRT */



 
