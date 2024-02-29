! This is include file "diags_ek.h"
!  ==== == ======= ==== ==========
!

#ifdef DIAGNOSTICS_EK

      real ekHadv(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekHadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekHdiff(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekHdiff(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVadv(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekVadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekCor(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekCor(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekPrsgrd(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekPrsgrd(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekHmix(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekHmix(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVmix(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekVmix(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekrate(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekrate(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekvol(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekvol(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVmix2(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekVmix2(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekWind(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekWind(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekDrag(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekDrag(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined DIAGNOSTICS_BARO
      real ekBaro(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekBaro(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
# if defined M3FAST
      real ekfast(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekfast(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
# ifdef AVERAGES
      real timediags_ek_avg
      real ekHadv_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekHadv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekHdiff_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekHdiff_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVadv_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekVadv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekCor_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekCor_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekPrsgrd_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekPrsgrd_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekHmix_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekHmix_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVmix_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekVmix_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekrate_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekrate_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekvol_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekvol_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVmix2_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekVmix2_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekWind_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekWind_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekDrag_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekDrag_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
#  if defined DIAGNOSTICS_BARO
      real ekBaro_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekBaro_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
#  endif
#  if defined M3FAST
      real ekfast_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekfast_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
#  endif
# endif /* AVERAGES */

      common /diag_ekHadv/ekHadv   
     &       /diag_ekHdiff/ekHdiff
     &       /diag_ekVadv/ekVadv     
     &       /diag_ekCor/ekCor
     &       /diag_ekPrsgrd/ekPrsgrd
     &       /diag_ekHmix/ekHmix
     &       /diag_ekVmix/ekVmix
     &       /diag_ekrate/ekrate
     &       /diag_ekvol/ekvol
     &       /diag_ekVmix2/ekVmix2
     &       /diag_ekWind/ekWind
     &       /diag_ekDrag/ekDrag
# if defined DIAGNOSTICS_BARO
     &       /diag_ekBaro/ekBaro
# endif
# if defined M3FAST
     &       /diag_ekfast/ekfast
# endif

# ifdef AVERAGES
      common /diag_timediags_ek_avg/timediags_ek_avg
      common /diag_ekHadv_avg/ekHadv_avg
     &       /diag_ekHdiff_avg/ekHdiff_avg
     &       /diag_ekVadv_avg/ekVadv_avg     
     &       /diag_ekCor_avg/ekCor_avg
     &       /diag_ekPrsgrd_avg/ekPrsgrd_avg
     &       /diag_ekHmix_avg/ekHmix_avg
     &       /diag_ekVmix_avg/ekVmix_avg
     &       /diag_ekrate_avg/ekrate_avg
     &       /diag_ekvol_avg/ekvol_avg
     &       /diag_ekVmix2_avg/ekVmix2_avg
     &       /diag_ekWind_avg/ekWind_avg
     &       /diag_ekDrag_avg/ekDrag_avg
#  if defined DIAGNOSTICS_BARO
     &       /diag_ekBaro_avg/ekBaro_avg
#  endif
#  if defined M3FAST
     &       /diag_ekfast_avg/ekfast_avg
#  endif
# endif  /* AVERAGES */   

# ifdef DIAGNOSTICS_EK_MLD
      real ekHadv_mld(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekHadv_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekHdiff_mld(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekHdiff_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVadv_mld(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekVadv_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekCor_mld(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekCor_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekPrsgrd_mld(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekPrsgrd_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekHmix_mld(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekHmix_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVmix_mld(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekVmix_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekrate_mld(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekrate_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekvol_mld(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekvol_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVmix2_mld(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekVmix2_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
#  if defined DIAGNOSTICS_BARO
      real ekBaro_mld(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekBaro_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
#  endif

#  ifdef AVERAGES
      real ekHadv_mld_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekHadv_mld_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekHdiff_mld_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekHdiff_mld_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVadv_mld_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekVadv_mld_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekCor_mld_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekCor_mld_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekPrsgrd_mld_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekPrsgrd_mld_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekHmix_mld_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekHmix_mld_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVmix_mld_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekVmix_mld_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekrate_mld_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekrate_mld_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekvol_mld_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekvol_mld_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVmix2_mld_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekVmix2_mld_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
#   if defined DIAGNOSTICS_BARO
      real ekBaro_mld_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekBaro_mld_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
#   endif
#  endif /* AVERAGES */

      common /diag_ekHadv_mld/ekHadv_mld
     &       /diag_ekHdiff_mld/ekHdiff_mld
     &       /diag_ekVadv_mld/ekVadv_mld
     &       /diag_ekCor_mld/ekCor_mld
     &       /diag_ekPrsgrd_mld/ekPrsgrd_mld
     &       /diag_ekHmix_mld/ekHmix_mld
     &       /diag_ekVmix_mld/ekVmix_mld
     &       /diag_ekrate_mld/ekrate_mld
     &       /diag_ekvol_mld/ekvol_mld
     &       /diag_ekVmix2_mld/ekVmix2_mld
#  if defined DIAGNOSTICS_BARO
     &       /diag_ekBaro_mld/ekBaro_mld
#  endif
#  ifdef AVERAGES
      common /diag_ekHadv_mld_avg/ekHadv_mld_avg
     &       /diag_ekHdiff_mld_avg/ekHdiff_mld_avg
     &       /diag_ekVadv_mld_avg/ekVadv_mld_avg     
     &       /diag_ekCor_mld_avg/ekCor_mld_avg
     &       /diag_ekPrsgrd_mld_avg/ekPrsgrd_mld_avg
     &       /diag_ekHmix_mld_avg/ekHmix_mld_avg
     &       /diag_ekVmix_mld_avg/ekVmix_mld_avg
     &       /diag_ekrate_mld_avg/ekrate_mld_avg
     &       /diag_ekvol_mld_avg/ekvol_mld_avg
     &       /diag_ekVmix2_mld_avg/ekVmix2_mld_avg
#   if defined DIAGNOSTICS_BARO
     &       /diag_ekBaro_mld_avg/ekBaro_mld_avg
#   endif
#  endif
# endif  /* DIAGNOSTICS_EK_MLD */


      real ekwrkHadv(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkHadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkHdiff(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkHdiff(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkVadv(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkVadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkCor(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkCor(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkPrsgrd(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkPrsgrd(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkHmix(GLOBAL_2D_ARRAY,2,2)
!CSDISTRIBUTE_RESHAPE ekwrkHmix(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkVmix(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkVmix(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkrate(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkrate(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkvol(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkvol(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkVmix2(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkVmix2(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkwind(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkwind(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkdrag(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkdrag(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined DIAGNOSTICS_BARO
      real ekwrkBaro(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkBaro(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
# if defined M3FAST
      real ekwrkfast(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkfast(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
      common /diag_ekwrkHadv/ekwrkHadv   
     &       /diag_ekwrkHdiff/ekwrkHdiff
     &       /diag_ekwrkVadv/ekwrkVadv     
     &       /diag_ekwrkCor/ekwrkCor
     &       /diag_ekwrkPrsgrd/ekwrkPrsgrd
     &       /diag_ekwrkHmix/ekwrkHmix
     &       /diag_ekwrkVmix/ekwrkVmix
     &       /diag_ekwrkrate/ekwrkrate
     &       /diag_ekwrkvol/ekwrkvol
     &       /diag_ekwrkVmix2/ekwrkVmix2
     &       /diag_ekwrkwind/ekwrkwind
     &       /diag_ekwrkdrag/ekwrkdrag
# if defined DIAGNOSTICS_BARO
     &       /diag_ekwrkBaro/ekwrkBaro
# endif
# if defined M3FAST
     &       /diag_ekwrkfast/ekwrkfast
# endif
# ifdef DIAGNOSTICS_EK_MLD
      real ekwrkHadv_mld(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkHadv_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkHdiff_mld(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkHdiff_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkVadv_mld(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkVadv_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkCor_mld(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkCor_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkPrsgrd_mld(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkPrsgrd_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkHmix_mld(GLOBAL_2D_ARRAY,2,2)
!CSDISTRIBUTE_RESHAPE ekwrkHmix_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkVmix_mld(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkVmix_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkrate_mld(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkrate_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkvol_mld(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkvol_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkVmix2_mld(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkVmix2_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
#  if defined DIAGNOSTICS_BARO
      real ekwrkBaro_mld(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkBaro_mld(BLOCK_PATTERN,*) BLOCK_CLAUSE
#  endif

      common /diag_ekwrkHadv_mld/ekwrkHadv_mld
     &       /diag_ekwrkHdiff_mld/ekwrkHdiff_mld
     &       /diag_ekwrkVadv_mld/ekwrkVadv_mld
     &       /diag_ekwrkCor_mld/ekwrkCor_mld
     &       /diag_ekwrkPrsgrd_mld/ekwrkPrsgrd_mld
     &       /diag_ekwrkHmix_mld/ekwrkHmix_mld
     &       /diag_ekwrkVmix_mld/ekwrkVmix_mld
     &       /diag_ekwrkrate_mld/ekwrkrate_mld
     &       /diag_ekwrkvol_mld/ekwrkvol_mld
     &       /diag_ekwrkVmix2_mld/ekwrkVmix2_mld
#  if defined DIAGNOSTICS_BARO
     &       /diag_ekwrkBaro_mld/ekwrkBaro_mld
#  endif
# endif /* DIAGNOSTICS_EK_MLD */

# if defined DIAGNOSTICS_EK_FULL && ! defined DIAGNOSTICS_UV
      real MXadv(GLOBAL_2D_ARRAY,N,2)
      real MYadv(GLOBAL_2D_ARRAY,N,2)
      real MHdiff(GLOBAL_2D_ARRAY,N,2)
      real MVadv(GLOBAL_2D_ARRAY,N,2)
      real MCor(GLOBAL_2D_ARRAY,N,2)
      real MPrsgrd(GLOBAL_2D_ARRAY,N,2)
      real MHmix(GLOBAL_2D_ARRAY,N,2,2)
      real MVmix(GLOBAL_2D_ARRAY,N,2)
      real Mbody(GLOBAL_2D_ARRAY,N,2)
#  if defined DIAGNOSTICS_BARO
      real MBaro(GLOBAL_2D_ARRAY,N,2)
#  endif
#  if defined M3FAST
      real Mfast(GLOBAL_2D_ARRAY,N,2)
#  endif

      common /diag_MXadv/MXadv
     &       /diag_MYadv/MYadv
     &       /diag_MHdiff/MHdiff
     &       /diag_MVadv/MVadv
     &       /diag_MCor/MCor
     &       /diag_MPrsgrd/MPrsgrd
     &       /diag_MHmix/MHmix
     &       /diag_MVmix/MVmix
     &       /diag_Mbody/Mbody
#  if defined DIAGNOSTICS_BARO
     &       /diag_MBaro/MBaro
#  endif
#  if defined M3FAST
     &       /diag_Mfast/Mfast
#  endif
# endif /* DIAGNOSTICS_EK_FULL */

# endif /* DIAGNOSTICS_EK */



 
