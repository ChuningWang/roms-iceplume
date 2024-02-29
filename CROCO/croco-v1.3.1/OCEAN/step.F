! $Id: step.F 1618 2014-12-18 14:39:51Z rblod $
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
!
!====================================================================
!                         subroutine step
!====================================================================
!
#include "cppdefs.h"
      subroutine step()
      implicit none
#include "param.h"
#include "scalars.h"
#include "zoom.h"
#include "grid.h"
#include "coupling.h"
#include "ocean3d.h"
#include "ocean2d.h"
#include "mpi_cpl.h"
#ifdef MUSTANG
# include "coupler_define_MUSTANG.h"
#endif

#ifdef AGRIF
      IF (agrif_fixed().NE.sortedint(nbtimes)) return
      nbtimes = nbtimes + 1
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef SOLVE3D

#if defined OA_COUPLING || defined OW_COUPLING
!--Get phase of OASIS3
       if ( (iif==-1).and.(oasis_time>=0).and.(nbstep3d<ntimes-ntstart+1) ) then
!         print *, 'oasis time before get = ', oasis_time
         call cpl_prism_get (oasis_time)
#ifdef AGRIF
       else if ((.not.agrif_root()).and.(iif==-1).and.
     &          (oasis_time>=0).and.(nbtimes<=Agrif_MaxLevel()+1).and.
     &          (nbprttime<nbmaxprttime)) then
         call cpl_prism_get (oasis_time)
#endif
       endif 
#endif /* OA_COUPLING */

       iif = iif+1
      IF (iif .EQ. 0) then
#endif
C$OMP PARALLEL
        call prestep3D_thread()
C$OMP END PARALLEL
#ifdef SOLVE3D
      endif
      IF ((iif.GT.0).AND.(iif.LT.(nfast+1))) THEN
#endif
C$OMP PARALLEL
        call step2d_thread()
C$OMP END PARALLEL
#ifdef SOLVE3D
      ENDIF
      IF (iif .EQ. (nfast+1)) then
#endif
C$OMP PARALLEL
        call step3D_uv_thread()
C$OMP END PARALLEL
#ifdef SOLVE3D
      endif
      IF (iif .EQ. (nfast + 2)) then
#endif
C$OMP PARALLEL
        call step3D_t_thread()
C$OMP END PARALLEL
#ifdef SOLVE3D
        iif = -1
        nbstep3d = nbstep3d + 1
      endif
#if defined OA_COUPLING || defined OW_COUPLING
!--Get phase of OASIS3
       if (iif==-1) then
         if (oasis_time>=0.and.(nbstep3d<=ntimes-ntstart+1)) then
           call cpl_prism_put (oasis_time)
#ifdef AGRIF
         else if (oasis_time>=0.and.(.not.agrif_root()).and.
     &           ((nbtimes>nbmaxtimes+1-Agrif_MaxLevel())).and.
     &           (nbprttime<nbmaxprttime) ) then
           call cpl_prism_put (oasis_time)
#endif
         endif
         oasis_time = oasis_time + dt
       endif
#endif /* OA_COUPLING */
#endif
!#ifdef ONLINE_ANALYSIS
!     if (mynode.eq.0) write(6,*) 'appel',iic,iif
!     call croco_oa(1)
!#endif
      return
      end
!
!====================================================================
!                subroutine prestep3D_thread
!====================================================================
!
      subroutine prestep3D_thread()
#ifdef AGRIF
      use Agrif_Util
# ifdef XIOS
# include "zoom.h"
# endif
#endif
#ifdef ONLINE_ANALYSIS
      use module_interface_oa, only : if_oa
#endif
      implicit none
#include "param.h"
#include "scalars.h"
#include "ncscrum.h"
#ifdef M3FAST
# include "nbq.h"
#endif
#ifdef FLOATS
# include "floats.h"
#endif
#include "ocean2d.h"
#include "ocean3d.h"
#if defined RVTK_DEBUG && ( defined BULK_FLUX || defined MRL_WCI )
# include "forces.h"
#endif
#if ( defined RVTK_DEBUG || defined  RVTK_DEBUG_PERFRST )
# include "mixing.h"
#endif

#ifdef STATIONS
# include "sta.h"
#endif
#if defined WKB_WWAVE
# include "wkb_wwave.h"
      integer iif_wave,winterp
#endif
      integer range
      integer ntrds,trd,
     &        my_first,my_last, tile, my_iif
      integer itrc

#ifdef FLOATS
      integer chunk_size_flt, Lstr,Lend, flt_str
      common /floats_step/ flt_str
#endif
#ifdef STATIONS
      integer chunk_size_sta, Mcpstr,Mcpend, sta_str
      common /sta_step/ sta_str
#endif
      integer omp_get_num_threads, omp_get_thread_num

      ntrds=omp_get_num_threads()
      trd=omp_get_thread_num()
C$OMP BARRIER
      range=(NSUB_X*NSUB_E+ntrds-1)/ntrds
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)

C$OMP MASTER
#ifdef AUTOTILING
      call begin_timing
#endif
      time=time_start+dt*float(iic-ntstart)
      tdays=time*sec2day
#ifdef USE_CALENDAR
      time_mars=time_mars+dt
      tdays=(time_mars-time_start)*sec2day
#endif

#ifdef SOLVE3D
      nstp=1+mod(iic-ntstart,2)
      nrhs=nstp
      nnew=3
#endif

#ifdef AGRIF
# ifdef XIOS
!      if( nbcoarse == 1 )then
        call iom_swap( "crocox" )
!      endif
# endif
#endif

#ifdef FLOATS
      nfp1=MOD(nfp1+1,NFT+1)  ! Shift time indices
      nf  =MOD(nf  +1,NFT+1)  ! for floats
      nfm1=MOD(nfm1+1,NFT+1)
      nfm2=MOD(nfm2+1,NFT+1)
      nfm3=MOD(nfm3+1,NFT+1)
      flt_str=0
#endif
#ifdef STATIONS
      sta_str=0
#endif

!
! Model input block: read forcing/climatology data from netCDF files.
!------ ----- ------ ---- ------------------- ---- ---- ------ ------
!
      if (synchro_flag) then
#ifdef AGRIF
        IF (Agrif_Root()) THEN
#endif
#if defined SOLVE3D && defined TCLIMATOLOGY && !defined ANA_TCLIMA
            call get_tclima
#endif
#if defined SOLVE3D && defined M3CLIMATOLOGY && !defined ANA_M3CLIMA \
 || (defined M2CLIMATOLOGY && !defined ANA_M2CLIMA)
            call get_uclima
#endif
#if defined ZCLIMATOLOGY && !defined ANA_SSH
            call get_ssh
#endif
#if defined FRC_BRY && !defined ANA_BRY
            call get_bry
# ifdef BIOLOGY
            call get_bry_bio
# endif
#endif
#if defined WKB_WWAVE && !defined ANA_BRY_WKB
            call get_bry_wkb
#endif
#ifdef AGRIF
        ENDIF
#endif
        call get_vbc
        synchro_flag=.false.
      endif  !<-- synchro_flag
C$OMP END MASTER
C$OMP BARRIER

#if defined RVTK_DEBUG_PERFRST
       call check_tab3d(Hz_bak(:,:,:),'Hz_bak beg prestep3d','r')
#endif

      if (may_day_flag.ne.0) return  !-->  EXIT

#ifdef SOLVE3D

      do tile=my_first,my_last
# ifdef AGRIF
        IF (Agrif_Root()) THEN
# endif
# ifdef TCLIMATOLOGY
          call ana_tclima (tile)   ! analytical values are given
#  ifndef ANA_TCLIMA
          call set_tclima (tile)   ! if data are missing from clim files
#  endif
# endif
# if defined M2CLIMATOLOGY || defined M3CLIMATOLOGY
#  if defined ANA_M2CLIMA || defined ANA_M3CLIMA
          call ana_uclima (tile)
#  else
          call set_uclima (tile)
#  endif
# endif
# ifdef ZCLIMATOLOGY
#  ifdef ANA_SSH
          call ana_ssh (tile)
#  else
          call set_ssh (tile)
#  endif
# endif
# ifdef FRC_BRY
#  ifdef ANA_BRY
          call ana_bry (tile)
#  else
          call set_bry (tile)
#  endif
#  ifdef BIOLOGY
          call ana_bry_bio (tile)  ! analytical values are given
#   ifndef ANA_BRY
          call set_bry_bio (tile)  ! if data are missing from bry files
#   endif
#  endif
# endif
# if defined WKB_WWAVE && !defined ANA_BRY_WKB
          call set_bry_wkb (tile)
# endif

# ifdef NBQCLIMATOLOGY
          call ana_nbq_clima (tile)
# elif defined NBQ_FRC_BRY
          call ana_nbq_bry (tile)
# endif

# ifdef AGRIF
        ENDIF
# endif
# ifdef ANA_WWAVE
        call ana_wwave (tile)
# endif
        call rho_eos (tile)
        call set_HUV (tile)
# if defined BIOLOGY
        call bio_diag (tile)
# else
        call diag(tile)
# endif
      enddo
C$OMP BARRIER
      do tile=my_first,my_last
        call set_vbc (tile)
# if defined  BBL && !defined M3FAST
#   ifdef AGRIF
        IF (Agrif_Fixed().GE.Agrif_lev_sedim) THEN
#   endif
        call bblm (tile)
#  ifdef AGRIF
        ENDIF
#  endif
# endif
# if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
#  ifdef AGRIF
        IF (Agrif_Root()) call clm_tides (tile)
#  else
        call clm_tides (tile)
#  endif
# endif
      enddo
C$OMP BARRIER

#  if defined RVTK_DEBUG && defined RVTK_DEBUG_ADVANCED && defined BULK_FLUX
C$OMP MASTER
      call check_tab2d(uwndg(:,:,1),'uwndg1 prestep3d_befbulk','u')
      call check_tab2d(uwndg(:,:,2),'uwndg2 prestep3d_befbulk','u')
      call check_tab2d(vwndg(:,:,1),'vwndg1 prestep3d_befbulk','v')
      call check_tab2d(vwndg(:,:,2),'vwndg2 prestep3d_befbulk','v')
      call check_tab2d(tair(:,:),'tair prestep3d_befbulk','r')
      call check_tab2d(rhum(:,:),'rhum prestep3d_befbulk','r')
      call check_tab2d(prate(:,:),'prate prestep3d_befbulk','r')
      call check_tab2d(radlw(:,:),'radlw prestep3d_befbulk','r')
      call check_tab2d(radsw(:,:),'radsw prestep3d_befbulk','r')
      call check_tab2d(wspd(:,:),'wspd prestep3d_befbulk','r')
      call check_tab2d(uwnd(:,:),'uwnd prestep3d_befbulk','u')
      call check_tab2d(vwnd(:,:),'vwnd prestep3d_befbulk','v')
C$OMP END MASTER
# endif

# ifdef BULK_FLUX
      do tile=my_first,my_last
        call bulk_flux (tile)   
      enddo
C$OMP BARRIER
#  if defined RVTK_DEBUG && defined RVTK_DEBUG_ADVANCED && defined BULK_FLUX
C$OMP MASTER
      call check_tab2d(uwndg(:,:,1),'uwndg1 prestep3d_aftbulk','u')
      call check_tab2d(uwndg(:,:,2),'uwndg2 prestep3d_aftbulk','u')
      call check_tab2d(vwndg(:,:,1),'vwndg1 prestep3d_aftbulk','v')
      call check_tab2d(vwndg(:,:,2),'vwndg2 prestep3d_aftbulk','v')

      call check_tab2d(tair(:,:),'tair prestep3d_aftbulk','r')
      call check_tab2d(rhum(:,:),'rhum prestep3d_aftbulk','r')
      call check_tab2d(prate(:,:),'prate prestep3d_aftbulk','r')
      call check_tab2d(radlw(:,:),'radlw prestep3d_aftbulk','r')
      call check_tab2d(radsw(:,:),'radsw prestep3d_aftbulk','r')
      call check_tab2d(wspd(:,:),'wspd prestep3d_aftbulk','r')
      call check_tab2d(uwnd(:,:),'uwnd prestep3d_aftbulk','u')
      call check_tab2d(vwnd(:,:),'vwnd prestep3d_aftbulk','v')
C$OMP END MASTER
C$OMP BARRIER
#  endif
# endif /* BULK_FLUX */

C$OMP BARRIER
# if defined RVTK_DEBUG_ADVANCED && defined BULK_FLUX
C$OMP MASTER      
      call check_tab2d(sustr(:,:),'sustr prestep3d_aftbulk','u')
      call check_tab2d(svstr(:,:),'svstr prestep3d_aftbulk','v')
C$OMP END MASTER
C$OMP BARRIER
# endif

      if (may_day_flag.ne.0) return  !-->  EXIT

# if defined WKB_WWAVE && !defined WKB_STEADY
      wkb_agrif_done=.FALSE.
#  ifdef MRL_CEW
      do tile=my_tile_range
        call wkb_cew_prep (tile)
      enddo
C$OMP BARRIER
C$OMP SINGLE
      wint=0
C$OMP END SINGLE
      do winterp=1,interp_max
C$OMP SINGLE
        wint=wint+1
        if (wint.gt.2) wint=1
C$OMP END SINGLE
        do tile=my_tile_range
          call wkb_uvfield (tile, winterp)
        enddo
C$OMP BARRIER
      enddo
#  endif
      do iif_wave=1,ndtfast         ! WKB ray equation w/ barotropic time stepping
C$OMP SINGLE
        wstp=wnew
        wnew=wstp+1
        if (wnew.ge.3) wnew=1
C$OMP END SINGLE

#  ifdef MRL_CEW
        do tile=my_tile_range
          call wkb_cew_finalize (tile)
        enddo
C$OMP BARRIER
#  endif
        do tile=my_tile_range
          call wkb_wwave (tile)
        enddo
C$OMP BARRIER
      enddo
#  if defined RVTK_DEBUG || defined RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
      call check_tab2d(wac(:,:,wnew),'wac in prestep ','r')
#   ifdef  WAVE_ROLLER
      call check_tab2d(war(:,:,wnew),'war in prestep ','r')
#   endif
C$OMP END MASTER
#  endif
# endif

# ifdef MRL_WCI
      do tile=my_tile_range
        call  mrl_wci (tile)
      enddo
C$OMP BARRIER
#  if defined RVTK_DEBUG || defined RVTK_DEBUG_ADVANCED
C$OMP MASTER
      call check_tab2d(ust2d(:,:),'ust2d in prestep ','u')
      call check_tab2d(vst2d(:,:),'vst2d in prestep ','v')
C$OMP END MASTER
#  endif
# endif

      do tile=my_first,my_last
# if defined ANA_VMIX
        call ana_vmix (tile)
# elif defined LMD_MIXING
        call lmd_vmix (tile)
# elif defined BVF_MIXING
        call bvf_mix (tile)
# endif
        call omega (tile)
# ifdef VIS_COEF_3D
        call hvisc_coef (tile)
# endif
      enddo
C$OMP BARRIER

      do tile=my_first,my_last
        call prsgrd (tile)
        call rhs3d (tile)
# ifdef NBQ
#  if defined AGRIF && defined NBQ_CHILD_ONLY
        if (.Not.Agrif_Root()) call rhs3d_w_nh (tile)
#  else
        call rhs3d_w_nh (tile)
#  endif
# endif
        call pre_step3d (tile)
# ifdef AGRIF
      if (.Not.Agrif_Root()) then
        call uv3dpremix (tile)
      endif
# endif
      enddo
C$OMP BARRIER

      do tile=my_first,my_last
# if defined UV_VIS2 || defined UV_VIS4
        call uv3dmix (tile)
# endif
# if defined SPONGE_VIS2 && !defined UV_VIS2
        call uv3dmix_spg (tile)
# endif
# ifdef SST_SKIN
        call sstskin (tile)
# endif

!---------------------------
# if defined DIAGNOSTICS_VRT
        call set_diags_vrt (tile)
# endif
# if defined DIAGNOSTICS_EK
        call set_diags_ek (tile)
# endif
# if defined DIAGNOSTICS_PV
        call set_diags_pv (tile)
# endif

# if defined AVERAGES && ! defined XIOS
        call set_avg (tile)
#  if defined DIAGNOSTICS_TS
        if (ldefdia_avg
     &      .and. (wrtdia3D_avg(NT+1) .or. wrtdia2D_avg(NT+1)))
     &      call set_diags_avg (tile)
#  endif
#  if defined DIAGNOSTICS_UV
        if (ldefdiaM_avg .and. wrtdiaM_avg(3))
     &      call set_diagsM_avg (tile)
#  endif
#  if defined DIAGNOSTICS_VRT
        if (ldefdiags_vrt_avg .and. wrtdiags_vrt_avg(indxTime))
     &      call set_diags_vrt_avg (tile)
#  endif
#  if defined DIAGNOSTICS_EK
        if (ldefdiags_ek_avg .and. wrtdiags_ek_avg(indxTime))
     &      call set_diags_ek_avg (tile)
#  endif
#  if defined DIAGNOSTICS_PV
        if (ldefdiags_pv_avg .and. wrtdiags_pv_avg(indxTime))
     &      call set_diags_pv_avg (tile)
#  endif
#  if defined DIAGNOSTICS_EDDY
        if (ldefdiags_eddy_avg .and. wrtdiags_eddy_avg(indxTime))
     &      call set_diags_eddy_avg (tile)
#  endif
#  ifdef DIAGNOSTICS_BIO
        call set_bio_diags_avg (tile)
#  endif
#  if defined OUTPUTS_SURFACE
        if (ldefsurf_avg .and. wrtsurf_avg(indxTime))
     &      call set_surf_avg (tile)
#  endif
# endif
      enddo

C$OMP BARRIER
C$OMP MASTER

#if defined RVTK_DEBUG || defined RVTK_DEBUG_PERFRST
      do itrc=1,NT
       call check_tab3d(t(:,:,:,nnew,itrc),'t prestep3d','r')
      enddo
      call check_tab3d(u(:,:,:,nnew),'u prestep3d','u')
      call check_tab3d(v(:,:,:,nnew),'v prestep3d','v')
#endif

      nrhs=3
      nnew=3-nstp

#ifdef ONLINE_ANALYSIS
! BLXD TODO moving OA calls just after the end of diags
      if ( if_oa.eqv..true. ) then
C$OMP END MASTER
C$OMP PARALLEL DO PRIVATE(tile)
         do tile=my_first,my_last
           ! BLXD call croco_oa(0)
           call online_spectral_diags(tile,1)
         enddo
         call output_oa
C$OMP MASTER
      endif
#endif

!
! Output block: write restart/history files.
!
      call output

C$OMP END MASTER
C$OMP BARRIER
#endif  /* SOLVE3D */

      if (may_day_flag.ne.0) return  !-->  EXIT


      return
      end
!
!====================================================================
!                   subroutine step2D_thread
!====================================================================
!
      subroutine step2D_thread()

#ifdef AGRIF
      use Agrif_Util
#endif
      implicit none
#include "param.h"
#include "scalars.h"
#include "ncscrum.h"
#include "ocean2d.h"
!#if ( defined RVTK_DEBUG || defined RVTK_DEBUG_PERFRST ) && defined MRL_WCI
!# include "forces.h"
!#endif
#ifdef FLOATS
# include "floats.h"
#endif
#ifdef STATIONS
# include "sta.h"
#endif
#if defined BBL && defined BSTRESS_FAST
# include "nbq.h"
#elif ( defined RVTK_DEBUG || defined  RVTK_DEBUG_PERFRST ) && defined M3FAST
# include "nbq.h"
#endif
#if defined RVTK_DEBUG || defined  RVTK_DEBUG_PERFRST
# include "forces.h"
#endif
      integer range
      integer ntrds,trd,
     &        my_first,my_last, tile, my_iif

#ifdef FLOATS
      integer chunk_size_flt, Lstr,Lend, flt_str
      common /floats_step/ flt_str
#endif
#ifdef STATIONS
      integer chunk_size_sta, Mcpstr,Mcpend, sta_str
      common /sta_step/ sta_str
#endif
      integer omp_get_num_threads, omp_get_thread_num

      ntrds=omp_get_num_threads()
      trd=omp_get_thread_num()
C$OMP BARRIER
      range=(NSUB_X*NSUB_E+ntrds-1)/ntrds
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)

!
! Solve the 2D primitive equations for the barotropic mode.
!-------------------------------------------------------------
! Note that no actual fast-time-step is performed during the
! auxiliary (nfast+1)th step. It is needed to finalize the fast-time
! averaging procedure, if any, and to compute the total depth of
! water column, as well as the new vertical coordinate transformation
! profiles, z_r, z_w, because the free surface elevation has been
! changed. All these operations are done during predictor time step.
! Thus, there is no need for the corrector step during the auxiliary
! time step.
!

#ifndef SOLVE3D
# if defined WKB_WWAVE && !defined WKB_STEADY
        wstp=wnew
        wnew=wstp+1
        if (wnew.ge.3) wnew=1
        do tile=my_tile_range
          call wkb_wwave (tile)
        enddo
C$OMP BARRIER
# endif
# ifdef MRL_WCI
        do tile=my_tile_range
          call mrl_wci (tile)
        enddo
C$OMP BARRIER
#  if defined RVTK_DEBUG || defined RVTK_DEBUG_ADVANCED
C$OMP MASTER
      call check_tab2d(ust2d(:,:),'ust2d in prestep ','u')
      call check_tab2d(vst2d(:,:),'vst2d in prestep ','v')
C$OMP END MASTER
#  endif
# endif
#endif

C$OMP MASTER
        kstp=knew                   ! This might look a bit silly,
        knew=kstp+1                 ! since both branches of this
        if (knew.gt.4) knew=1       ! if statement are identical.
C$OMP END MASTER
C$OMP BARRIER
!
! Model input block for 2D configurations only !!!
!
#ifndef SOLVE3D
# ifdef AGRIF
        IF (Agrif_Root()) THEN
# endif
          do tile=my_first,my_last
# if defined M2CLIMATOLOGY
#  if defined ANA_M2CLIMA
            call ana_uclima (tile)
#  else
            call set_uclima (tile)
#  endif
# endif
# ifdef ZCLIMATOLOGY
#  ifdef ANA_SSH
            call ana_ssh (tile)
#  else
            call set_ssh (tile)
#  endif
# endif
# if defined Z_FRC_BRY  || defined M2_FRC_BRY
#  ifdef ANA_BRY
            call ana_bry (tile)
#  else
            call set_bry (tile)
#  endif
# endif
# if defined WKB_WWAVE && !defined ANA_BRY_WKB
            call set_bry_wkb (tile)
# endif
          enddo
# ifdef AGRIF
        ENDIF
# endif
C$OMP BARRIER
        do tile=my_first,my_last
          call set_vbc (tile)
# ifdef BBL
          call bblm (tile)
# endif
# if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
#   ifdef AGRIF
          IF (Agrif_Root()) call clm_tides (tile)
#   else
          call clm_tides (tile)
#   endif
# endif
# ifdef AVERAGES
          call set_avg (tile)
# endif
# if defined UV_VIS2 && defined UV_VIS_SMAGO
        call hvisc_coef (tile)
# endif
        enddo
C$OMP BARRIER
C$OMP MASTER
        call output
C$OMP END MASTER
C$OMP BARRIER
#endif  /* !defined SOLVE3D */

#if defined BBL && defined BSTRESS_FAST
        if (mod(iif-1,inc_faststep).eq.0) then
          do tile=my_first,my_last
            call bblm (tile)
          enddo
        endif
C$OMP BARRIER
#endif
                  ! This might look a bit silly,
                  ! since both branches of this
                  ! if statement are identical.
                  ! Nevertheless, it makes sense,
                  ! since mpc will reverse one of
                  ! these loops to make zig-zag
                  ! tile-processing sequence.
        if (mod(knew,2).eq.0) then
          do tile=my_first,my_last
#ifdef M3FAST
# if defined AGRIF && defined NBQ_CHILD_ONLY
            if (Agrif_Root()) then
              call   step2d (tile)
            else
              call   step3d_fast (tile)
            endif
# else
            call   step3d_fast (tile)
# endif
#else
            call     step2d (tile)
#endif
          enddo
        else
          do tile=my_first,my_last
#ifdef M3FAST
# if defined AGRIF && defined NBQ_CHILD_ONLY
            if (Agrif_Root()) then
              call   step2d (tile)
            else
              call   step3d_fast (tile)
            endif
# else
            call   step3d_fast (tile)
# endif
#else
            call     step2d (tile)
#endif
          enddo
        endif

#if defined AGRIF && defined AGRIF_2WAY
        if (.not.agrif_root()) then
            call update2d()
        endif
#endif

#ifdef RVTK_DEBUG
C$OMP BARRIER
C$OMP MASTER
# ifdef M3FAST
       call check_tab2d(zeta(:,:,knew),'zeta step3d_fast','r')
       call check_tab2d(ubar(:,:,knew),'ubar step3d_fast','u')
       call check_tab2d(vbar(:,:,knew),'vbar step3d_fast','v')
#  ifdef RVTK_DEBUG_ADVANCED
        call check_tab2d(bustr,'bustr step3d_fast','u')
        call check_tab2d(bvstr,'bvstr step3d_fast','v')
        !!call check_tab3d(qdmu_nbq,'qdmu_nbq step3d_fast','u')
        !!call check_tab3d(qdmv_nbq,'qdmv_nbq step3d_fast','v')
#  endif
# else
       call check_tab2d(zeta(:,:,knew),'zeta step2d','r')
       call check_tab2d(ubar(:,:,knew),'ubar step2d','u')
       call check_tab2d(vbar(:,:,knew),'vbar step2d','v')
# endif
C$OMP END MASTER
#endif

      return
      end
!
!====================================================================
!                  subroutine step3D_uv_thread
!====================================================================
!
      subroutine step3D_uv_thread()
#ifdef AGRIF
      use Agrif_Util
#endif
      implicit none
#include "param.h"
#include "scalars.h"
#include "ncscrum.h"
#ifdef M3FAST
# include "nbq.h"
#endif
#ifdef FLOATS
# include "floats.h"
#endif
#ifdef STATIONS
# include "sta.h"
#endif
#include "zoom.h"
#include "ocean3d.h"
#include "ocean2d.h"
#include "coupling.h"
#ifdef AUTOTILING
# include "autotiling.h"
#endif
#if defined LMD_MIXING
# include "mixing.h"
#elif defined RVTK_DEBUG_PERFRST && defined TS_MIX_ISO_FILT
# include "mixing.h"
#endif
       integer range
      integer ntrds,trd,
     &        my_first,my_last, tile, my_iif
      integer i,j,k

#ifdef FLOATS
      integer chunk_size_flt, Lstr,Lend, flt_str
      common /floats_step/ flt_str
#endif
#ifdef STATIONS
      integer chunk_size_sta, Mcpstr,Mcpend, sta_str
      common /sta_step/ sta_str
#endif
      integer omp_get_num_threads, omp_get_thread_num
      integer itrc

      real t1,t2,t3

      ntrds=omp_get_num_threads()
      trd=omp_get_thread_num()
C$OMP BARRIER
      range=(NSUB_X*NSUB_E+ntrds-1)/ntrds
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)

C$OMP BARRIER

#ifdef SOLVE3D
# if defined NBQ && defined AGRIF && defined NBQ_CHILD_ONLY
      if (Agrif_Root()) then
        do tile=my_first,my_last
          call set_depth (tile)
        enddo
      endif
# elif !defined NBQ
      do tile=my_first,my_last
        call set_depth (tile)
      enddo
# endif
C$OMP BARRIER

      do tile=my_first,my_last
        call set_HUV2 (tile)
      enddo
C$OMP BARRIER

      do tile=my_first,my_last

        call omega (tile)

        call rho_eos (tile)

      enddo

C$OMP BARRIER

      do tile=my_first,my_last
        call prsgrd (tile)
        call rhs3d (tile)
# ifdef NBQ
#  if defined AGRIF && defined NBQ_CHILD_ONLY
        if (.Not.Agrif_Root()) call rhs3d_w_nh (tile)
#  else
        call rhs3d_w_nh (tile)
#  endif
# endif
        call step3d_uv1 (tile)
# ifdef NBQ
#  if defined AGRIF && defined NBQ_CHILD_ONLY
        if (.Not.Agrif_Root()) call step3d_w (tile)
#  else
        call step3d_w (tile)
#  endif
# endif
      enddo
C$OMP BARRIER

      do tile=my_first,my_last
        call step3d_uv2 (tile)
# ifdef DIF_COEF_3D
        call hdiff_coef (tile)
# endif
      enddo
C$OMP BARRIER

#endif /* SOLVE3D */
!
#if defined AGRIF && defined AGRIF_2WAY
!
! Update the outer domain after the last child step
! in case of 2-way nesting.
!
C$OMP BARRIER
C$OMP MASTER
      if (.Not.Agrif_Root()) then
        call Agrif_update_uv_np1
      endif
C$OMP END MASTER
#endif /*AGRIF && AGRIF_2WAY*/

      return
      end

!
!====================================================================
!                  subroutine step3D_t_thread
!====================================================================
!
      subroutine step3D_t_thread()
#ifdef AGRIF
      use Agrif_Util
#endif
#ifdef MUSTANG
      use plug_MUSTANG_CROCO, only: mustang_update_main
      use plug_MUSTANG_CROCO, only: mustang_deposition_main
# ifdef MORPHODYN
      use plug_MUSTANG_CROCO, only: mustang_morpho_main
      use comMUSTANG, only : l_morphocoupl, t_morpho
# endif
# if defined RVTK_DEBUG_AVANCED && defined SOLVE3D
      use comMUSTANG, only : c_sedtot
      use comMUSTANG, only : dzs 
# if defined key_MUSTANG_bedload && defined key_MUSTANG_V2
      use comMUSTANG , only : flx_bx, flx_by 
      use comsubstance, only : ibedload1, ibedload2
# endif
#endif
# include "coupler_define_MUSTANG.h"
#endif
      implicit none
#include "param.h"
#include "scalars.h"
#include "ncscrum.h"
#ifdef FLOATS
# include "floats.h"
#endif
#ifdef STATIONS
# include "sta.h"
#endif
#include "zoom.h"
#include "ocean3d.h"
#include "ocean2d.h"
#include "coupling.h"
#include "grid.h"
#ifdef AUTOTILING
# include "autotiling.h"
#endif
#if ( defined RVTK_DEBUG || defined RVTK_DEBUG_PERFRST )
# include "mixing.h"
#endif
      integer range
      integer ntrds,trd,
     &        my_first,my_last, tile, my_iif
      integer i,j,k

#ifdef FLOATS
      integer chunk_size_flt, Lstr,Lend, flt_str
      common /floats_step/ flt_str
#endif
#ifdef STATIONS
      integer chunk_size_sta, Mcpstr,Mcpend, sta_str
      common /sta_step/ sta_str
#endif
      integer omp_get_num_threads, omp_get_thread_num
      integer itrc, ilay

      real t1,t2,t3

      ntrds=omp_get_num_threads()
      trd=omp_get_thread_num()
C$OMP BARRIER
      range=(NSUB_X*NSUB_E+ntrds-1)/ntrds
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)
C$OMP BARRIER

#ifdef SOLVE3D

      do tile=my_first,my_last

        call omega (tile)

# ifdef GLS_MIXING
        call gls_mixing (tile)
# endif

# ifdef MUSTANG
        call mustang_update_main (tile)
# endif
        call step3d_t (tile)

# ifdef MUSTANG
        call mustang_deposition_main (tile)
# endif

# ifdef SEDIMENT
#  ifdef AGRIF
        IF (Agrif_Fixed().GE.Agrif_lev_sedim)
     &  call sediment (tile)
#  else
        call sediment (tile)
#  endif
# endif
# ifdef AGRIF
        if (.Not.Agrif_Root()) then
          call t3dpremix (tile)
        endif
# endif
      enddo
C$OMP BARRIER

# if defined RVTK_DEBUG_AVANCED && defined SOLVE3D
C$OMP MASTER
      do itrc = 1, NT
        call check_tab3d(Akt(:,:,1:N,itrc),'Akt step3d_t','r')
      enddo
#  if defined GLS_MIXING
        call check_tab3d(trb(:,:,1:N,nstp,igls),'gls step3d_t','r')
        call check_tab3d(trb(:,:,1:N,nstp,itke),'tke step3d_t','r')
#  endif
# if defined MUSTANG 
        do ilay = ksdmin, ksdmax
            call check_tab2d(c_sedtot(ilay, :, :),'c_sedtotx step3d_t','r')
            call check_tab2d(dzs(ilay, :, :),'dzs step3d_t','r')
        enddo
# if defined key_MUSTANG_bedload && defined key_MUSTANG_V2
        do itrc = ibedload1, ibedload2
            call check_tab2d(flx_bx(itrc, :, :),'flx_bx step3d_t','r')
            call check_tab2d(flx_by(itrc, :, :),'flx_by step3d_t','r')
        enddo
# endif
# endif
C$OMP END MASTER
# endif

# if defined TRACERS
      do tile=my_first,my_last
#  if defined TS_DIF2 || defined TS_DIF4
        call t3dmix (tile)
#  endif
#  if defined SPONGE_DIF2
        call t3dmix_spg (tile)
#  endif
      enddo
C$OMP BARRIER

# endif /* TRACERS */

# ifdef MORPHODYN
      do tile=my_first,my_last
#  ifdef SEDIMENT
#   ifdef ANA_MORPHODYN
        call ana_morphodyn (tile)
#   endif
#  elif defined MUSTANG
        if (l_morphocoupl .AND. CURRENT_TIME .GE. t_morpho ) then
            call mustang_morpho_main (tile)
        endif
#  endif
        call set_depth_morphodyn (tile)
      enddo
C$OMP BARRIER

# endif

#endif /* SOLVE3D */

#ifdef FLOATS
      chunk_size_flt=32
      do while (flt_str.lt.nfloats)
C$OMP CRITICAL
        Lstr=flt_str+1
        flt_str=Lstr+chunk_size_flt-1
C$OMP END CRITICAL
        Lend=min(Lstr+chunk_size_flt-1,nfloats)
        call step_floats (Lstr,Lend)
      enddo
c        call step_floats (1,nfloats) ! serial version for debugging
#endif /*FLOATS*/

#ifdef STATIONS
      chunk_size_sta=32
      do while (sta_str.lt.nstas)
C$OMP CRITICAL
        Mcpstr=sta_str+1
        sta_str=Mcpstr+chunk_size_sta-1
C$OMP END CRITICAL
        Mcpend=min(Mcpstr+chunk_size_sta-1,nstas)
        call step_sta (Mcpstr,Mcpend)
      enddo
#endif /*STATIONS*/
!
#if defined AGRIF && defined AGRIF_2WAY
!
! Update the outer domain after the last child step
! in case of 2-way nesting.
!
C$OMP BARRIER
C$OMP MASTER
      if ((.Not.Agrif_Root()).and.
     &    (nbcoarse == Agrif_Irhot())) then
        call Agrif_update_np1
      endif
C$OMP END MASTER
#endif /*AGRIF && AGRIF_2WAY*/

#if ( defined RVTK_DEBUG || defined RVTK_DEBUG_PERFRST ) && defined SOLVE3D
C$OMP BARRIER
C$OMP MASTER
      do itrc=1,NT
        call check_tab3d(t(:,:,:,nnew,itrc),'t step3d','r')
        !!!call check_tab3d(Akt(:,:,1:N,itrc),'Akt aft step3d_t','r')
      enddo
      call check_tab3d(u(:,:,:,nnew),'u step3d','u')
      call check_tab3d(v(:,:,:,nnew),'v step3d','v')
C$OMP END MASTER
#endif

C$OMP MASTER
#ifdef AUTOTILING
      call end_timing
      if (iic-lastiic == nbsampling) call auto_tiling
#endif

      iic=iic + 1

#ifdef AGRIF
      nbcoarse = 1 + mod(nbcoarse, Agrif_IRhot())
#endif
C$OMP END MASTER
C$OMP BARRIER

      return
      end
