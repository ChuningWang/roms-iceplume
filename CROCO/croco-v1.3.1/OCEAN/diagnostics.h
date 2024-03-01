! $Id: diagnostics.h 1458 2014-02-03 15:01:25Z gcambon $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
! This is include file "diagnostics.h": tracer equation terms
! for output purposes:
!
!
#ifdef DIAGNOSTICS_TS
      real TXadv(GLOBAL_2D_ARRAY,N,NT)
      real TYadv(GLOBAL_2D_ARRAY,N,NT)
      real TVadv(GLOBAL_2D_ARRAY,N,NT)
      real THmix(GLOBAL_2D_ARRAY,N,NT)
      real TVmix(GLOBAL_2D_ARRAY,N,NT)
#ifdef DIAGNOSTICS_TSVAR
      real TVmixt(GLOBAL_2D_ARRAY,N,NT)
#endif
      real TForc(GLOBAL_2D_ARRAY,N,NT)
      real Trate(GLOBAL_2D_ARRAY,N,NT)
!
# ifdef DIAGNOSTICS_TS_MLD
      real TXadv_mld(GLOBAL_2D_ARRAY,NT)
      real TYadv_mld(GLOBAL_2D_ARRAY,NT)
      real TVadv_mld(GLOBAL_2D_ARRAY,NT)
      real THmix_mld(GLOBAL_2D_ARRAY,NT)
      real TVmix_mld(GLOBAL_2D_ARRAY,NT)
      real TForc_mld(GLOBAL_2D_ARRAY,NT)
      real Trate_mld(GLOBAL_2D_ARRAY,NT)
      real Tentr_mld(GLOBAL_2D_ARRAY,NT)
      integer kbl_nstp(GLOBAL_2D_ARRAY)
# endif
# ifdef AVERAGES
      real timedia_avg
      real TXadv_avg(GLOBAL_2D_ARRAY,N,NT)
      real TYadv_avg(GLOBAL_2D_ARRAY,N,NT)
      real TVadv_avg(GLOBAL_2D_ARRAY,N,NT)
      real THmix_avg(GLOBAL_2D_ARRAY,N,NT)
      real TVmix_avg(GLOBAL_2D_ARRAY,N,NT)
#ifdef DIAGNOSTICS_TSVAR
      real TVmixt_avg(GLOBAL_2D_ARRAY,N,NT)
# endif
      real TForc_avg(GLOBAL_2D_ARRAY,N,NT)
      real Trate_avg(GLOBAL_2D_ARRAY,N,NT)
!
#  ifdef DIAGNOSTICS_TS_MLD
      real TXadv_mld_avg(GLOBAL_2D_ARRAY,NT)
      real TYadv_mld_avg(GLOBAL_2D_ARRAY,NT)
      real TVadv_mld_avg(GLOBAL_2D_ARRAY,NT)
      real THmix_mld_avg(GLOBAL_2D_ARRAY,NT)
      real TVmix_mld_avg(GLOBAL_2D_ARRAY,NT)
      real TForc_mld_avg(GLOBAL_2D_ARRAY,NT)
      real Trate_mld_avg(GLOBAL_2D_ARRAY,NT)
      real Tentr_mld_avg(GLOBAL_2D_ARRAY,NT)
#  endif
# endif	
      common /diag_TXadv/TXadv   
      common /diag_TYadv/TYadv
      common /diag_TVadv/TVadv  
      common /diag_THmix/THmix
      common /diag_TVmix/TVmix
#ifdef DIAGNOSTICS_TSVAR
      common /diag_TVmixt/TVmixt
# endif
      common /diag_TForc/TForc
      common /diag_Trate/Trate
!
# ifdef DIAGNOSTICS_TS_MLD
      common /diag_TXadv_mld/TXadv_mld
      common /diag_TYadv_mld/TYadv_mld
      common /diag_TVadv_mld/TVadv_mld
      common /diag_THmix_mld/THmix_mld
      common /diag_TVmix_mld/TVmix_mld
      common /diag_TForc_mld/TForc_mld
      common /diag_Trate_mld/Trate_mld
      common /diag_Tentr_mld/Tentr_mld	  
      common /diag_kbl_nstp/kbl_nstp
# endif
# ifdef AVERAGES
      common /diag_timedia_avg/timedia_avg
      common /diag_TXadv_avg/TXadv_avg  
      common /diag_TYadv_avg/TYadv_avg
      common /diag_TVadv_avg/TVadv_avg   
      common /diag_THmix_avg/THmix_avg
      common /diag_TVmix_avg/TVmix_avg
#ifdef DIAGNOSTICS_TSVAR
      common /diag_TVmixt_avg/TVmixt_avg
# endif
      common /diag_TForc_avg/TForc_avg
      common /diag_Trate_avg/Trate_avg
!
#  ifdef DIAGNOSTICS_TS_MLD
      common /diag_TXadv_mld_avg/TXadv_mld_avg
      common /diag_TYadv_mld_avg/TYadv_mld_avg
      common /diag_TVadv_mld_avg/TVadv_mld_avg
      common /diag_THmix_mld_avg/THmix_mld_avg
      common /diag_TVmix_mld_avg/TVmix_mld_avg
      common /diag_TForc_mld_avg/TForc_mld_avg
      common /diag_Trate_mld_avg/Trate_mld_avg
      common /diag_Tentr_mld_avg/Tentr_mld_avg
#  endif       	
# endif       	
#endif /* DIAGNOSTICS_TS */
!
#ifdef DIAGNOSTICS_UV
      real MXadv(GLOBAL_2D_ARRAY,N,2)
      real MYadv(GLOBAL_2D_ARRAY,N,2)
      real MVadv(GLOBAL_2D_ARRAY,N,2)
      real MCor(GLOBAL_2D_ARRAY,N,2)
      real MPrsgrd(GLOBAL_2D_ARRAY,N,2)
      real MHmix(GLOBAL_2D_ARRAY,N,2,2)
      real MHdiff(GLOBAL_2D_ARRAY,N,2)
      real MVmix(GLOBAL_2D_ARRAY,N,2)
      real MVmix2(GLOBAL_2D_ARRAY,N,2)
      real Mrate(GLOBAL_2D_ARRAY,N,2)
      real Mbody(GLOBAL_2D_ARRAY,N,2)
# if defined DIAGNOSTICS_BARO
      real MBaro(GLOBAL_2D_ARRAY,N,2)
# endif
# if defined M3FAST
      real Mfast(GLOBAL_2D_ARRAY,N,2)
# endif
# ifdef MRL_WCI
      real Mvf(GLOBAL_2D_ARRAY,N,2)
      real Mbrk(GLOBAL_2D_ARRAY,N,2)
      real MStCo(GLOBAL_2D_ARRAY,N,2)
      real MVvf(GLOBAL_2D_ARRAY,N,2)
      real MPrscrt(GLOBAL_2D_ARRAY,N,2)
      real Msbk(GLOBAL_2D_ARRAY,N,2) 
      real Mbwf(GLOBAL_2D_ARRAY,N,2)
      real Mfrc(GLOBAL_2D_ARRAY,N,2)
# endif      
# ifdef AVERAGES
      real timediaM_avg
      real MXadv_avg(GLOBAL_2D_ARRAY,N,2)
      real MYadv_avg(GLOBAL_2D_ARRAY,N,2)
      real MVadv_avg(GLOBAL_2D_ARRAY,N,2)
      real MCor_avg(GLOBAL_2D_ARRAY,N,2)
      real MPrsgrd_avg(GLOBAL_2D_ARRAY,N,2)
      real MHmix_avg(GLOBAL_2D_ARRAY,N,2)
      real MHdiff_avg(GLOBAL_2D_ARRAY,N,2)
      real MVmix_avg(GLOBAL_2D_ARRAY,N,2)
      real MVmix2_avg(GLOBAL_2D_ARRAY,N,2)
      real Mrate_avg(GLOBAL_2D_ARRAY,N,2)
# if defined DIAGNOSTICS_BARO
      real MBaro_avg(GLOBAL_2D_ARRAY,N,2)
# endif
# if defined M3FAST
      real Mfast_avg(GLOBAL_2D_ARRAY,N,2)
# endif
#  ifdef MRL_WCI
      real Mvf_avg(GLOBAL_2D_ARRAY,N,2)
      real Mbrk_avg(GLOBAL_2D_ARRAY,N,2)
      real MStCo_avg(GLOBAL_2D_ARRAY,N,2)
      real MVvf_avg(GLOBAL_2D_ARRAY,N,2)
      real MPrscrt_avg(GLOBAL_2D_ARRAY,N,2)
      real Msbk_avg(GLOBAL_2D_ARRAY,N,2) 
      real Mbwf_avg(GLOBAL_2D_ARRAY,N,2)
      real Mfrc_avg(GLOBAL_2D_ARRAY,N,2)
#  endif      
# endif	
      common /diag_MXadv/MXadv   
      common /diag_MYadv/MYadv
      common /diag_MHdiff/MHdiff
      common /diag_MVadv/MVadv  
      common /diag_MCor/MCor
      common /diag_MPrsgrd/MPrsgrd
      common /diag_MHmix/MHmix
      common /diag_MVmix/MVmix
      common /diag_MVmix2/MVmix2
      common /diag_Mrate/Mrate
      common /diag_Mbody/Mbody
# if defined DIAGNOSTICS_BARO
      common /diag_MBaro/MBaro
# endif
# if defined M3FAST
      common /diag_Mfast/Mfast
# endif
# ifdef MRL_WCI       
      common /diag_Mvf/Mvf
      common /diag_Mbrk/Mbrk
      common /diag_MStCo/MStCo
      common /diag_MVvf/MVvf
      common /diag_MPrscrt/MPrscrt
      common /diag_Msbk/Msbk
      common /diag_Mbwf/Mbwf
      common /diag_Mfrc/Mfrc
# endif      
# ifdef AVERAGES
      common /diag_timediaM_avg/timediaM_avg
      common /diag_MXadv_avg/MXadv_avg
      common /diag_MYadv_avg/MYadv_avg
      common /diag_MVadv_avg/MVadv_avg
      common /diag_MCor_avg/MCor_avg
      common /diag_MPrsgrd_avg/MPrsgrd_avg
      common /diag_MHmix_avg/MHmix_avg
      common /diag_MHdiff_avg/MHdiff_avg
      common /diag_MVmix_avg/MVmix_avg
      common /diag_MVmix2_avg/MVmix2_avg
      common /diag_Mrate_avg/Mrate_avg
# if defined DIAGNOSTICS_BARO
      common /diag_MBaro_avg/MBaro_avg
# endif
# if defined M3FAST
      common /diag_Mfast_avg/Mfast_avg
# endif
#  ifdef MRL_WCI       
      common /diag_Mvf_avg/Mvf_avg
      common /diag_Mbrk_avg/Mbrk_avg
      common /diag_MStCo_avg/MStCo_avg
      common /diag_MVvf_avg/MVvf_avg
      common /diag_MPrscrt_avg/MPrscrt_avg
      common /diag_Msbk_avg/Msbk_avg
      common /diag_Mbwf_avg/Mbwf_avg
      common /diag_Mfrc_avg/Mfrc_avg
#  endif      
# endif       	
#endif /* DIAGNOSTICS_UV */
#ifdef DIAGNOSTICS_BIO
# ifdef PISCES 
#  ifdef key_trc_diaadd
      real bioFlux(GLOBAL_2D_ARRAY,N,NumFluxTerms)
      real bioVSink(GLOBAL_2D_ARRAY,NumVSinkTerms)
#  endif
# else
      real bioFlux(GLOBAL_2D_ARRAY,N,NumFluxTerms)
      real bioVSink(GLOBAL_2D_ARRAY,0:N,NumVSinkTerms)
#  if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS 
      real GasExcFlux(GLOBAL_2D_ARRAY,NumGasExcTerms)
#  endif
# endif
# ifdef AVERAGES
#  ifdef PISCES 
#    ifdef key_trc_diaadd
      real bioFlux_avg(GLOBAL_2D_ARRAY,N,NumFluxTerms)
      real bioVSink_avg(GLOBAL_2D_ARRAY,NumVSinkTerms)
#    endif
#  else
      real bioFlux_avg(GLOBAL_2D_ARRAY,N,NumFluxTerms)
      real bioVSink_avg(GLOBAL_2D_ARRAY,0:N,NumVSinkTerms)
#  if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS 
      real GasExcFlux_avg(GLOBAL_2D_ARRAY,NumGasExcTerms)
#    endif
#  endif
      real timediabio_avg
# endif
# ifdef PISCES 
#    ifdef key_trc_diaadd
      common /diag_bioFlux/bioFlux
      common /diag_bioVSink/bioVSink
#    endif
# else
      common /diag_bioFlux/bioFlux
      common /diag_bioVSink/bioVSink
#  if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS 
      common /diag_GasFlux/GasExcFlux
#   endif
# endif
# ifdef AVERAGES
#  ifdef PISCES 
#   ifdef key_trc_diaadd
      common /diag_bioFlux_avg/bioFlux_avg
      common /diag_bioVSink_avg/bioVSink_avg
#   endif
#  else
      common /diag_bioFlux_avg/bioFlux_avg
      common /diag_bioVSink_avg/bioVSink_avg
#  if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      common /diag_GasFlux_avg/GasExcFlux_avg
#   endif 
#  endif
      common /diag_timediabio_avg/timediabio_avg
# endif
#endif /* DIAGNOSTICS_BIO */

