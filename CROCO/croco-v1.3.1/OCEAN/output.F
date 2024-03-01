! $Id: output.F 1586 2014-07-30 14:57:11Z marchesiello $
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

#include "cppdefs.h"
      subroutine output
!
!=======================================================================
!                                                                      !
! Output block: write restart/history files.
!                                                                      !
!=======================================================================
!
#ifdef AGRIF
      use Agrif_Util
#endif
      implicit none
      integer ilc
#include "param.h"
#include "scalars.h"
#include "ncscrum.h"
#ifdef FLOATS
# include "floats.h"
# include "ncscrum_floats.h"
#endif
#ifdef STATIONS
# include "sta.h"
# include "nc_sta.h"
#endif
!
                integer ntrds,trd,range,
     &        my_first,my_last, tile, my_iif

      integer omp_get_num_threads, omp_get_thread_num
      ilc=1+iic-ntstart   ! number of time step since restart

!
! XIOS Case
!
#ifdef XIOS
      ntrds=omp_get_num_threads()
      trd=omp_get_thread_num()
      range=(NSUB_X*NSUB_E+ntrds-1)/ntrds
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)
      
      do tile=my_first,my_last
          call send_xios_diags(tile)
      end do
#endif
!
!
! Conditions for writing output
!
      
      if ( iic.gt.ntstart .and. ( mod(ilc-1,nrst).eq.0
#ifdef EXACT_RESTART
     &                    .or.      mod(ilc,nrst).eq.0
#endif
     &     .or. (mod(ilc-1,nwrt).eq.0 .and. wrthis(indxTime))
#ifdef AVERAGES
     &    .or. (mod(ilc-ntsavg,navg).eq.0 .and. wrtavg(indxTime))
#endif
#if defined DIAGNOSTICS_TS
     &    .or. (mod(ilc-1,nwrtdia).eq.0 .and. wrtdia3D(NT+1))
# ifdef AVERAGES
     &    .or. (mod(ilc-ntsdia_avg,nwrtdia_avg).eq.0 .and. wrtdia3D_avg(NT+1))
# endif
# if defined DIAGNOSTICS_TS_MLD
     &    .or. (mod(ilc-1,nwrtdia).eq.0 .and. wrtdia2D(NT+1))
#  ifdef AVERAGES
     &    .or. (mod(ilc-ntsdia_avg,nwrtdia_avg).eq.0 .and. wrtdia2D_avg(NT+1))
#  endif
# endif
#endif
#if defined DIAGNOSTICS_UV
     &    .or. (mod(ilc-1,nwrtdiaM).eq.0 .and. wrtdiaM(3))
# ifdef AVERAGES
     &    .or. (mod(ilc-ntsdiaM_avg,nwrtdiaM_avg).eq.0 .and. wrtdiaM_avg(3))
# endif
#endif
#if defined DIAGNOSTICS_VRT
     &    .or. (mod(iic-ntstart,nwrtdiags_vrt).eq.0
     &    .and. wrtdiags_vrt(indxTime))
# ifdef AVERAGES
     &    .or. (mod(ilc-ntsdiags_vrt_avg,nwrtdiags_vrt_avg).eq.0
     &    .and. wrtdiags_vrt_avg(indxTime))
# endif
# endif
#if defined DIAGNOSTICS_EK
     &    .or. (mod(iic-ntstart,nwrtdiags_ek).eq.0
     &    .and. wrtdiags_ek(indxTime))
# ifdef AVERAGES
     &    .or. (mod(ilc-ntsdiags_ek_avg,nwrtdiags_ek_avg).eq.0
     &    .and. wrtdiags_ek_avg(indxTime))
# endif
# endif
#if defined DIAGNOSTICS_PV
     &    .or. (mod(iic-ntstart,nwrtdiags_pv).eq.0
     &    .and. wrtdiags_pv(indxTime))
# ifdef AVERAGES
     &    .or. (mod(ilc-ntsdiags_pv_avg,nwrtdiags_pv_avg).eq.0
     &    .and. wrtdiags_pv_avg(indxTime))
# endif
# endif
#if defined DIAGNOSTICS_EDDY && ! defined XIOS
# ifdef AVERAGES
     &    .or. (mod(ilc-ntsdiags_eddy_avg,nwrtdiags_eddy_avg).eq.0
     &    .and. wrtdiags_eddy_avg(indxTime))
# endif
# endif
#if defined OUTPUTS_SURFACE && ! defined XIOS
     &    .or. (mod(iic-ntstart,nwrtsurf).eq.0
     &    .and. wrtsurf(indxTime))
# ifdef AVERAGES
     &    .or. (mod(ilc-ntssurf_avg,nwrtsurf_avg).eq.0
     &    .and. wrtsurf_avg(indxTime))
# endif
# endif
#ifdef DIAGNOSTICS_BIO
     &    .or. (mod(ilc-1,nwrtdiabio).eq.0 .and. 
     &    (      wrtdiabioFlux(NumFluxTerms+1)
     &      .or. wrtdiabioVSink(NumVSinkTerms+1)
     &      .or. wrtdiabioGasExc(NumGasExcTerms+1) ) )

# ifdef AVERAGES
     &    .or. (mod(ilc-1,nwrtdiabio_avg).eq.0 .and. 
     &    (      wrtdiabioFlux_avg(NumFluxTerms+1)
     &      .or. wrtdiabioVSink_avg(NumVSinkTerms+1)
     &      .or. wrtdiabioGasExc_avg(NumGasExcTerms+1) ))
# endif
#endif /*DIAGNOSTICS_BIO*/
#ifdef FLOATS
# ifdef AGRIF
     &    .or. (Agrif_Root() .and. mod(ilc-1,nflt).eq.0
# else
     &    .or. (mod(ilc-1,nflt).eq.0
# endif
     &                             .and. nfloats.gt.0)
#endif
#ifdef STATIONS
# ifdef AGRIF
     &    .or. (Agrif_Root() .and. mod(ilc-1,nsta).eq.0
# else
     &    .or. (mod(ilc-1,nsta).eq.0
# endif
     &                             .and. nstas0.gt.0)
#endif
     &                                               )) then
!
! Get record number
!
        if (mod(ilc-1,nrst).eq.0
#ifdef EXACT_RESTART
     &                      .or. mod(ilc,nrst).eq.0
#endif
     &                      )     nrecrst=nrecrst+1
        if (mod(ilc-1,nwrt).eq.0) nrechis=nrechis+1
#ifdef AVERAGES
        if (mod(ilc-1,navg).eq.0) nrecavg=nrecavg+1
#endif
#if defined DIAGNOSTICS_TS
        if (mod(ilc-1,nwrtdia).eq.0) nrecdia=nrecdia+1
# ifdef AVERAGES
        if (mod(ilc-ntsdia_avg,nwrtdia_avg).eq.0) nrecdia_avg=nrecdia_avg+1
# endif
#endif
#if defined DIAGNOSTICS_UV
        if (mod(ilc-1,nwrtdiaM).eq.0) nrecdiaM=nrecdiaM+1
# ifdef AVERAGES
        if (mod(ilc-ntsdiaM_avg,nwrtdiaM_avg).eq.0) nrecdiaM_avg=nrecdiaM_avg+1
# endif
#endif
#if defined DIAGNOSTICS_VRT
        if (mod(iic-ntstart,nwrtdiags_vrt) .eq.0)
     &                                  nrecdiags_vrt=nrecdiags_vrt+1
# ifdef AVERAGES
        if (mod(ilc-ntsdiags_vrt_avg,nwrtdiags_vrt_avg) .eq.0)
     &                          nrecdiags_vrt_avg=nrecdiags_vrt_avg+1
# endif
#endif
#if defined DIAGNOSTICS_EK
        if (mod(iic-ntstart,nwrtdiags_ek) .eq.0)
     &                                    nrecdiags_ek=nrecdiags_ek+1
# ifdef AVERAGES
        if (mod(ilc-ntsdiags_ek_avg,nwrtdiags_ek_avg) .eq.0)
     &                            nrecdiags_ek_avg=nrecdiags_ek_avg+1
# endif
#endif
#if defined DIAGNOSTICS_PV
        if (mod(iic-ntstart,nwrtdiags_pv) .eq.0)
     &                                    nrecdiags_pv=nrecdiags_pv+1
# ifdef AVERAGES
        if (mod(ilc-ntsdiags_pv_avg,nwrtdiags_pv_avg) .eq.0)
     &                            nrecdiags_pv_avg=nrecdiags_pv_avg+1
# endif
#endif
#if defined DIAGNOSTICS_EDDY && ! defined XIOS
# ifdef AVERAGES
        if (mod(ilc-ntsdiags_eddy_avg,nwrtdiags_eddy_avg) .eq.0)
     &                     nrecdiags_eddy_avg=nrecdiags_eddy_avg+1
# endif
#endif
#if defined OUTPUTS_SURFACE && ! defined XIOS
        if (mod(iic-ntstart,nwrtsurf) .eq.0)
     &                                    nrecsurf=nrecsurf+1
# ifdef AVERAGES
        if (mod(ilc-ntssurf_avg,nwrtsurf_avg) .eq.0)
     &                            nrecsurf_avg=nrecsurf_avg+1
# endif
#endif
#ifdef DIAGNOSTICS_BIO
        if (mod(ilc-1,nwrtdiabio).eq.0) nrecdiabio=nrecdiabio+1
# ifdef AVERAGES
        if (mod(ilc-1,nwrtdiabio_avg).eq.0)
     &                          nrecdiabio_avg=nrecdiabio_avg+1
# endif
#endif /*DIAGNOSTICS_BIO*/
#ifdef FLOATS
# ifdef AGRIF
        if (Agrif_Root() .and. mod(ilc-1,nflt).eq.0)
# else
        if (mod(ilc-1,nflt).eq.0)
# endif
     &     nrecflt=nrecflt+1
#endif
#ifdef STATIONS
# ifdef AGRIF
        if (Agrif_Root() .and. mod(ilc-1,nsta).eq.0)
# else
        if (mod(ilc-1,nsta).eq.0)
# endif
     &     nrecsta=nrecsta+1
#endif
!
! Call output routines
!
!        write(*,*) "wrtdia3D ",wrtdia3D
!        write(*,*) "wrtdia2D ",wrtdia2D
!        write(*,*) "wrtdia_avg3D ",wrtdia3D_avg
!        write(*,*) "wrtdia_avg2D ",wrtdia2D_avg
!!
!        write(*,*) "wrtdiaM ",wrtdiaM
!        write(*,*) "wrtdiaM_avg ",wrtdiaM_avg
!
        if (mod(ilc-1,nrst).eq.0
#ifdef EXACT_RESTART
     &       .or. mod(ilc,nrst).eq.0
#endif
     &       ) call wrt_rst
#ifndef XIOS
        if (mod(ilc-1,nwrt).eq.0 .and. wrthis(indxTime)) call wrt_his
#  ifdef AVERAGES
        if (mod(ilc-ntsavg,navg).eq.0 .and. wrtavg(indxTime)) call wrt_avg
#  endif

#if defined DIAGNOSTICS_TS
        if (mod(ilc-1,nwrtdia).eq.0
     &      .and. ldefdia
     &      .and. (wrtdia3D(NT+1) .or.
     &      wrtdia2D(NT+1)))                              
     &                                                 call wrt_diags
# ifdef AVERAGES
        if (ldefdia_avg .and. mod(ilc-ntsdia_avg,nwrtdia_avg).eq.0
     &     .and. (wrtdia3D_avg(NT+1) .or. wrtdia2D_avg(NT+1))) 
     &                                             call wrt_diags_avg
# endif
#endif

#if defined DIAGNOSTICS_UV
        if (mod(ilc-1,nwrtdiaM).eq.0
     &      .and. ldefdiaM
     &      .and. wrtdiaM(3))
     &                                                call wrt_diagsM
# ifdef AVERAGES
        if (ldefdiaM_avg .and. mod(ilc-ntsdiaM_avg,nwrtdiaM_avg).eq.0
     &                 .and. wrtdiaM_avg(3))       call wrt_diagsM_avg
# endif
#endif

#if defined DIAGNOSTICS_VRT
        if (mod(iic-ntstart,nwrtdiags_vrt).eq.0
     &      .and. ldefdiags_vrt 
     &      .and. wrtdiags_vrt(indxTime)) 
     &                  call wrt_diags_vrt
# ifdef AVERAGES
        if (mod(ilc-ntsdiags_vrt_avg,nwrtdiags_vrt_avg) .eq.0
     &      .and. iic.gt.ntstart+1
     &      .and. ldefdiags_vrt_avg 
     &      .and. wrtdiags_vrt_avg(indxTime))
     &      call wrt_diags_vrt_avg
# endif
#endif

#if defined DIAGNOSTICS_EK
        if (mod(iic-ntstart,nwrtdiags_ek).eq.0
     &      .and. ldefdiags_ek 
     &      .and. wrtdiags_ek(indxTime)) 
     &                  call wrt_diags_ek
# ifdef AVERAGES
        if (mod(ilc-ntsdiags_ek_avg,nwrtdiags_ek_avg) .eq.0
     &      .and. iic.gt.ntstart+1
     &      .and. ldefdiags_ek_avg 
     &      .and. wrtdiags_ek_avg(indxTime))
     &      call wrt_diags_ek_avg
# endif
#endif

#if defined DIAGNOSTICS_PV
        if (mod(iic-ntstart,nwrtdiags_pv).eq.0
     &      .and. ldefdiags_pv 
     &      .and. wrtdiags_pv(indxTime)) 
     &                  call wrt_diags_pv
# ifdef AVERAGES
        if (mod(ilc-ntsdiags_pv_avg,nwrtdiags_pv_avg) .eq.0
     &      .and. iic.gt.ntstart+1
     &      .and. ldefdiags_pv_avg 
     &      .and. wrtdiags_pv_avg(indxTime))
     &      call wrt_diags_pv_avg
# endif
#endif

#if defined DIAGNOSTICS_EDDY && ! defined XIOS
# ifdef AVERAGES
        if (mod(ilc-ntsdiags_eddy_avg,nwrtdiags_eddy_avg) .eq.0
     &      .and. iic.gt.ntstart+1
     &      .and. ldefdiags_eddy_avg
     &      .and. wrtdiags_eddy_avg(indxTime))
     &      call wrt_diags_eddy_avg
# endif
#endif

#if defined OUTPUTS_SURFACE && ! defined XIOS
        if (mod(iic-ntstart,nwrtsurf).eq.0
     &      .and. ldefsurf 
     &      .and. wrtsurf(indxTime)) 
     &                  call wrt_surf
# ifdef AVERAGES
        if (mod(ilc-ntssurf_avg,nwrtsurf_avg) .eq.0
     &      .and. iic.gt.ntstart+1
     &      .and. ldefsurf_avg 
     &      .and. wrtsurf_avg(indxTime))
     &      call wrt_surf_avg
# endif
#endif

#if defined DIAGNOSTICS_BIO && ! defined XIOS
        if (mod(ilc-1,nwrtdiabio).eq.0 .and. 
     &       (      wrtdiabioFlux(NumFluxTerms+1)
     &       .or. wrtdiabioVSink(NumVSinkTerms+1)
     &       .or. wrtdiabioGasExc(NumGasExcTerms+1) ) ) then
          call wrt_bio_diags
        endif
# ifdef AVERAGES
        if (ldefdiabio_avg .and. mod(ilc-1,nwrtdiabio_avg).eq.0.and.   
     &       (      wrtdiabioFlux_avg(NumFluxTerms+1)
     &       .or. wrtdiabioVSink_avg(NumVSinkTerms+1)
     &       .or. wrtdiabioGasExc_avg(NumGasExcTerms+1) ) ) then
          call wrt_bio_diags_avg
        endif
# endif
#endif

#endif


#ifdef FLOATS
# ifdef AGRIF
        if (Agrif_Root()) then
# endif
          if (mod(ilc-1,nflt).eq.0
     &                .and.nfloats.gt.0) call wrt_floats
          diagfloats=.false.
# ifdef AGRIF
        endif
# endif
#endif
#ifdef STATIONS
# ifdef AGRIF
        if (Agrif_Root()) then
# endif
          if (mod(ilc-1,nsta).eq.0
     &                  .and.nstas0.gt.0) call wrt_sta
          diagsta=.false.
# ifdef AGRIF
        endif
# endif
#endif
!
      endif      ! conditions for writing output
!
#ifdef FLOATS
!
! Flag for diagnostic computation (writing at next time step)
!
# ifdef AGRIF
      if (Agrif_Root() .and. mod(ilc,nflt).eq.0) then
# else
      if (mod(ilc,nflt).eq.0) then
# endif
        diagfloats=.true.
      endif
#endif
#ifdef STATIONS
# ifdef AGRIF
      if (Agrif_Root() .and. mod(ilc,nsta).eq.0) then
# else
      if (mod(ilc,nsta).eq.0) then
# endif
        diagsta=.true.
      endif
#endif

      return
      end

