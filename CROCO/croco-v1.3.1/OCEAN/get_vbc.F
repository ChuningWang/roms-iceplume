! $Id: get_vbc.F 1615 2014-12-17 13:27:07Z rblod $
!
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"

      subroutine get_vbc
!
!--------------------------------------------------------------------
!  This subroutine reads various forcing fields fom NetCDF files and
! save them as globally accessable arrays (declared in common blocks
! in file "forces.h"):
!
!   sustrg    kinematic surface momentum flux (wind stress) in
!                 the XI-direction [m^2/s^2].
!   svstrg    kinematic surface momentum flux (wind stress) in
!                 the ETA-direction [m^2/s^2].
!   srflxg    kinematic surface shortwave solar radiation flux
!                 [degC m/s].
!   stflxg    kinematic surface flux of tracer type variables
!                 [temperature: degC m/s; salinity: PSU m/s].
!-------------------------------------------------------------------
!
      implicit none
#include "param.h"
#include "sources.h"
!

!---------------------------------------------------------------
! River discharge, temperature and salinity
!---------------------------------------------------------------
!
#ifdef PSOURCE_NCFILE
        call get_psource
# ifdef PSOURCE_NCFILE_TS
        call get_psource_ts
# endif
#endif
!
#ifndef OW_COUPLING
!
!---------------------------------------------------------------
! Initial wave periode, amplitude, direction and dissipation
!---------------------------------------------------------------
!
# if (defined MUSTANG || defined BBL || defined MRL_WCI) && defined WAVE_OFFLINE
#  if defined AGRIF && defined BBL
      if (Agrif_lev_sedim.EQ.0) call get_wwave
#  else
      call get_wwave
#  endif
# endif
#endif /* !OW_COUPLING */
!
#ifndef OA_COUPLING
!
!---------------------------------------------------------------
! Kinematic surface momentum flux (wind stress) components
!  "sustrg" and "svstrg" [m^2/s^2].
!---------------------------------------------------------------
!
# ifndef ANA_SMFLUX
#  if (!defined BULK_FLUX && !defined OW_COUPLING) \
         || (defined OW_COUPLING && !defined WAVE_SMFLUX)
      call get_smflux
#  endif
# endif
!
# ifdef SOLVE3D
!
!---------------------------------------------------------------
! Kinematic surface temperature (heat) flux [degC m/s].
!---------------------------------------------------------------
!
#  ifdef BULK_FLUX
#   ifdef ONLINE
      ! Load data for ONLINE bulk forcing
      call get_bulk_online
#   else
      ! Load data for OFFLINE bulk forcing
      call get_bulk
#   endif
#  elif !defined ANA_STFLUX
#   if defined TEMPERATURE
      call get_stflux (itemp)
#   endif
#  endif
!
!---------------------------------------------------------------
! Flux correction to surface net heat flux.
!---------------------------------------------------------------
!
#  if defined QCORRECTION && !defined ANA_SST && defined TEMPERATURE
      call get_sst
#  endif
!
!---------------------------------------------------------------
! Kinematic surface freshwater flux (E-P) flux [PSU m/s].
!---------------------------------------------------------------
!
#  ifndef BULK_FLUX
#   if defined SALINITY && !defined ANA_SSFLUX
       call get_stflux (isalt)
#   endif
#  endif
!
!---------------------------------------------------------------
! Flux correction to surface salt flux.
!---------------------------------------------------------------
!
#  if defined SALINITY && defined SFLX_CORR && !defined ANA_SSS
      call get_sss
#  endif
!
!---------------------------------------------------------------
! Kinematic surface solar shortwave radiation flux [degC m/s].
!---------------------------------------------------------------
!
#  ifndef BULK_FLUX
#   if defined LMD_SKPP || defined LMD_BKPP || defined GLS_MIXING
#    if !defined ANA_SRFLUX && defined TEMPERATURE
       call get_srflux
#    endif
#   endif
#  endif /* !BULK_FLUX */
# endif /* SOLVE3D */
#endif /* !OA_COUPLING */
!
!---------------------------------------------------------------
! Kinematic surface solar shortwave radiation flux [degC m/s].
!---------------------------------------------------------------
!
#if ( defined SOLVE3D && defined TEMPERATURE && defined BHFLUX \
     && !defined ANA_BTFLUX )
      call get_btflux (itemp)
#endif
!
!---------------------------------------------------------------
! Kinematic bottom salt flux [PSU m/s].
! (Analytical bottom salt flux is usually set to zero.)
!---------------------------------------------------------------
!
#if ( defined SOLVE3D && defined SALINITY && defined BWFLUX \
     && !defined ANA_BTFLUX )
      call get_btflux (isalt)
#endif

      return
      end
!
!====================================================================
!     subroutine set_vbc
!====================================================================
!
      subroutine set_vbc (tile)
      implicit none
#include "param.h"
#include "private_scratch.h"
      integer tile, trd, omp_get_thread_num
#include "compute_tile_bounds.h"
!
      trd=omp_get_thread_num()
      call set_vbc_tile (Istr,Iend,Jstr,Jend,
     &         A2d(1,1,trd),A2d(1,2,trd),A2d(1,3,trd))
      return
      end

      subroutine set_vbc_tile (Istr,Iend,Jstr,Jend, wrk,wrk1,wrk2)
!
!----------------------------------------------------------------
!
!  This subroutine sets the vertical boundary conditions for
!  momentum and tracers.
!
!  On Output:  stored in common blocks (see include file forces.h)
!
!  sustr       kinematic surface momentum flux (wind stress) in
!                                      the XI-direction (m^2/s^2).
!  svstr       kinematic surface momentum flux (wind stress) in
!                                     the ETA-direction (m^2/s^2).
!  srflx       kinematic surface shortwave solar radiation flux
!                                                      (degC m/s).
!  stflx       kinematic surface flux of tracer type variables
!                      (temperature: degC m/s; salinity: PSU m/s).
!  btflx       Kinematic bottom flux of tracer type variables
!                      (temperature: degC m/s; salinity: PSU m/s).
!  bustr       kinematic bottom momentum flux in
!                                      the XI-direction (m^2/s^2).
!  bvstr       kinematic bottom momentum flux in
!                                     the ETA-direction (m^2/s^2).
!----------------------------------------------------------------
!
#ifdef AGRIF
      use Agrif_Util
#endif
      implicit none
#include "param.h"
#include "grid.h"
#include "ocean3d.h"
#include "forces.h"
#include "scalars.h"
#include "climat.h"
#include "sources.h"
#include "mixing.h"
#include "mpi_cpl.h"
#include "sediment.h"
#include "bbl.h"
      integer i,j,is, Istr,Iend,Jstr,Jend
      integer ISTRm,ISTRp,IENDp,JSTRm,JSTRp,JENDp
      real cff,cff1, Umag,Ucur,Vcur,Utur,znots
      real wrk(PRIVATE_2D_SCRATCH_ARRAY)
      real wrk1(PRIVATE_2D_SCRATCH_ARRAY)
      real wrk2(PRIVATE_2D_SCRATCH_ARRAY)
!
#include "compute_auxiliary_bounds.h"
!
#ifndef OA_COUPLING
!
!---------------------------------------------------------------
! Kinematic surface momentum flux (wind stress) components
! "sustr" and "svstr" [m^2/s^2].
!---------------------------------------------------------------
!
# ifdef ANA_SMFLUX
      call ana_smflux_tile (Istr,Iend,Jstr,Jend)
# else
#  if (!defined BULK_FLUX && !defined OW_COUPLING) \
         || (defined OW_COUPLING && !defined WAVE_SMFLUX)
      call set_smflux_tile (Istr,Iend,Jstr,Jend)
#  endif
# endif

# ifdef SOLVE3D
!
!---------------------------------------------------------------
! Kinematic surface temperature (heat) flux [degC m/s] and
!           surface solar shortwave radiation flux [degC m/s]
!           surface freshwater (E-P) flux [PSU m/s].
!
! --- This sets the input variables for Bulk formulation ---
!  --- The Bulk routine is actually called from step.F ---
!---------------------------------------------------------------
!
#  ifdef TEMPERATURE
#   ifdef BULK_FLUX
#    ifdef ONLINE
      call set_bulk_tile_online (Istr,Iend,Jstr,Jend)
#    else
      call set_bulk_tile (Istr,Iend,Jstr,Jend)
#    endif
#   else
!
!---------------------------------------------------------------
! Kinematic surface temperature (heat) flux [degC m/s].
!---------------------------------------------------------------
!
#    ifdef ANA_STFLUX
      call ana_stflux_tile (Istr,Iend,Jstr,Jend, itemp)
#    else
      call set_stflux_tile (Istr,Iend,Jstr,Jend, itemp)
#    endif
!
!---------------------------------------------------------------
! Kinematic surface solar shortwave radiation flux [degC m/s].
!---------------------------------------------------------------
!
#    if defined LMD_SKPP || defined LMD_BKPP || defined GLS_MIXING
#     ifdef ANA_SRFLUX
      call ana_srflux_tile (Istr,Iend,Jstr,Jend)
#     else
      call set_srflux_tile (Istr,Iend,Jstr,Jend)
#     endif
#    endif
#   endif /* BULK_FLUX */
#  endif /* TEMPERATURE */
!
!---------------------------------------------------------------
! Simplest sea ice correction to T,S fluxes
!---------------------------------------------------------------
!
#ifdef SEA_ICE_NOFLUX
      do j=JstrR,JendR
        do i=IstrR,IendR
              if( t(i,j,N,nrhs,itemp) .le. -1.8 ) then
                   stflx(i,j,itemp)=0.
#   if defined LMD_SKPP || defined LMD_BKPP || defined GLS_MIXING
                   srflx(i,j)=0.
#    endif
              endif
        enddo
      enddo
#endif
!
!---------------------------------------------------------------
! Flux correction to surface net heat flux.
!---------------------------------------------------------------
!
#  if defined QCORRECTION && defined TEMPERATURE
#   ifdef ANA_SST
      call ana_sst_tile (Istr,Iend,Jstr,Jend)
#   else
      call set_sst_tile (Istr,Iend,Jstr,Jend)
#   endif
#   ifndef BULK_FLUX
!
! If BULK_FLUX is defined, correction is done in bulk_flux
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          stflx(i,j,itemp)=stflx(i,j,itemp)+
     &                  dqdt(i,j)*(t(i,j,N,nrhs,itemp)-sst(i,j))
        enddo
      enddo
#   endif
#  endif /* QCORRECTION */
!
!---------------------------------------------------------------
! Kinematic surface freshwater flux (E-P) flux [PSU m/s].
!---------------------------------------------------------------
!
#  if defined SALINITY && !defined BULK_FLUX
#   ifdef ANA_SSFLUX
      call ana_stflux_tile (Istr,Iend,Jstr,Jend, isalt)
#   else
      call set_stflux_tile (Istr,Iend,Jstr,Jend, isalt)
#   endif
!
! Multiply flux by surface salinity.
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          stflx(i,j,isalt)=stflx(i,j,isalt)*t(i,j,N,nrhs,isalt)
        enddo
      enddo
#  endif
!
!--------------------------------------------------------------
! Flux correction to surface salt flux.
!--------------------------------------------------------------
!
#  if defined SALINITY && defined SFLX_CORR
#   ifdef ANA_SSS
      call ana_sss_tile (Istr,Iend,Jstr,Jend)
#   else
      call set_sss_tile (Istr,Iend,Jstr,Jend)
#   endif
#   ifndef BULK_FLUX
!
! If BULK_FLUX is defined, correction is done in bulk_flux
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          stflx(i,j,isalt)=stflx(i,j,isalt)
# ifdef SFLX_CORR_COEF
     &  -Hz(i,j,N)/(dSdt*day2sec)*(t(i,j,N,nrhs,isalt)-sss(i,j))
# else
     &                 +dqdt(i,j)*(t(i,j,N,nrhs,isalt)-sss(i,j))
# endif
        enddo
      enddo
#   endif
#  endif /* SFLX_CORR */
!
!---------------------------------------------------------------
! Diurnal modulation of surface solar shortwave radiation flux
!---------------------------------------------------------------
!
#  ifdef TEMPERATURE
#   if defined LMD_SKPP || defined LMD_BKPP || defined GLS_MIXING
#    ifdef ANA_DIURNAL_SW
      call ana_diurnal_sw_tile (Istr,Iend,Jstr,Jend)
#    endif
#   endif
#  endif
!
!---------------------------------------------------------------
! Ensure that drying water does not receive heat/salt fluxes
!  to avoid unrealistic T,S values
!---------------------------------------------------------------
!
#  ifdef WET_DRY
      do j=JstrR,JendR
        do i=IstrR,IendR
#   ifdef TEMPERATURE
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask_wet(i,j)
          srflx(i,j)=srflx(i,j)*rmask_wet(i,j)
#   endif
!          sustr(i,j)=sustr(i,j)*umask_wet(i,j)
!          svstr(i,j)=svstr(i,j)*vmask_wet(i,j)
#   ifdef SALINITY
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask_wet(i,j)
#   endif
        enddo
      enddo
#  endif
# endif /* SOLVE3D */
#endif /* !OA_COUPLING */

#if defined SOLVE3D && defined TRACERS
!
!---------------------------------------------------------------
!     Simplest sea ice correction to T,S fluxes
!     in case of OA_COUPLING
!---------------------------------------------------------------
!      
#  if defined OA_COUPLING && defined SEA_ICE_NOFLUX
      do j=JstrR,JendR
        do i=IstrR,IendR
          if ( t(i,j,N,nrhs,itemp) .le. -1.8 ) then
                   stflx(i,j,itemp)=0.
#   if defined LMD_SKPP || defined LMD_BKPP || defined GLS_MIXING
                   srflx(i,j)=0.
#   endif
          endif
        enddo
      enddo
#endif  /* OA_COUPLING && defined SEA_ICE_NOFLUX */
!
!---------------------------------------------------------------
!  Kinematic bottom temperature (heat) flux [degC m/s].
!  (Analytical bottom heat flux is usually set to zero.)
!---------------------------------------------------------------
!
# ifdef TEMPERATURE
#  ifdef ANA_BTFLUX
      call ana_btflux_tile (Istr,Iend,Jstr,Jend, itemp)
#  else
      call set_btflux_tile (Istr,Iend,Jstr,Jend, itemp)
#  endif
# endif
!
!---------------------------------------------------------------
! Kinematic bottom salt flux [PSU m/s].
! (Analytical bottom salt flux is usually set to zero.)
!---------------------------------------------------------------
!
# ifdef SALINITY
#  ifdef ANA_BSFLUX
      call ana_btflux_tile (Istr,Iend,Jstr,Jend, isalt)
#  else
      call set_btflux_tile (Istr,Iend,Jstr,Jend, isalt)
#  endif
!
!---------------------------------------------------------------
! Multiply flux by bottom salinity.
!---------------------------------------------------------------
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          btflx(i,j,isalt)=btflx(i,j,isalt)*t(i,j,1,nrhs,isalt)
        enddo
      enddo
# endif /* SALINITY */
#endif   /* SOLVE3D */
!
!---------------------------------------------------------------
! Kinematic bottom momentum flux [m^2/s^2]
!
! if BBL is defined, bottom flux (bustr,bvstr) is computed in
! subroutine bbl (called in step after call to set_vbc).
!---------------------------------------------------------------
!
#if defined BBL && defined AGRIF
      IF (Agrif_Fixed().LT.Agrif_lev_sedim) THEN
#endif

#if (!defined BBL || defined AGRIF)
# ifdef ANA_BMFLUX

      call ana_bmflux_tile (Istr,Iend,Jstr,Jend)

# elif defined SOLVE3D

#  ifdef BSTRESS_FAST
!
!   --> Compute bustr & bvstr in step3d_fast
!
#  else
!
!  Set bottom stress using logarithmic or linear
!  and/or quadratic formulation.
!
#   define CDRAG wrk
      if (maxval(Zob).ne.0.) then
        do j=JstrV-1,Jend     !! for currents
          do i=IstrU-1,Iend
            cff=vonKar/LOG((z_r(i,j,1)-z_w(i,j,0))/Zob(i,j))
            CDRAG(i,j)=MIN(Cdb_max,MAX(Cdb_min,cff*cff))
          enddo
        enddo
        do j=Jstr,Jend
          do i=IstrU,Iend
            cff=0.25*(v(i  ,j,1,nrhs)+v(i  ,j+1,1,nrhs)+
     &                v(i-1,j,1,nrhs)+v(i-1,j+1,1,nrhs))
            Umag=SQRT(cff*cff+u(i,j,1,nrhs)*u(i,j,1,nrhs))
            bustr(i,j)=0.5*(CDRAG(i-1,j)+CDRAG(i,j))*Umag*u(i,j,1,nrhs)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            cff=0.25*(u(i,j  ,1,nrhs)+u(i+1,j  ,1,nrhs)+
     &                u(i,j-1,1,nrhs)+u(i+1,j-1,1,nrhs))
            Umag=SQRT(cff*cff+v(i,j,1,nrhs)*v(i,j,1,nrhs))
            bvstr(i,j)=0.5*(CDRAG(i,j-1)+CDRAG(i,j))*Umag*v(i,j,1,nrhs)
          enddo
        enddo
#   undef CDRAG
      elseif (rdrg2.gt.0.) then
        do j=JstrV,Jend
          do i=Istr,Iend
            cff=0.25*(v(i  ,j,1,nrhs)+v(i  ,j+1,1,nrhs)+
     &                v(i-1,j,1,nrhs)+v(i-1,j+1,1,nrhs))
            Umag=SQRT(cff*cff+u(i,j,1,nrhs)*u(i,j,1,nrhs))
            bustr(i,j)=rdrg2*Umag*u(i,j,1,nrhs)
          enddo
        enddo
        do j=Jstr,Jend
          do i=IstrU,Iend
            cff=0.25*(u(i,j  ,1,nrhs)+u(i+1,j  ,1,nrhs)+
     &                u(i,j-1,1,nrhs)+u(i+1,j-1,1,nrhs))
            Umag=SQRT(cff*cff+v(i,j,1,nrhs)*v(i,j,1,nrhs))
            bvstr(i,j)=rdrg2*Umag*v(i,j,1,nrhs)
          enddo
        enddo
      else
        do j=Jstr,Jend
          do i=Istr,Iend
            bustr(i,j)=rdrg*u(i,j,1,nrhs)
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            bvstr(i,j)=rdrg*v(i,j,1,nrhs)
          enddo
        enddo
      endif

#   ifdef LIMIT_BSTRESS
!
!  From J. Warner's code:
!  Set limiting factor for bottom stress. The bottom stress is adjusted
!  to not change the direction of momentum. It only should slow down
!  to zero. The value of 0.75 is arbitrary limitation assigment.
!
      cff=0.75/dt
      do j=Jstr,Jend
        do i=IstrU,Iend
          cff1=cff*0.5*(Hz(i-1,j,1)+Hz(i,j,1))
          bustr(i,j)=SIGN(1.D0, bustr(i,j))*
     &               MIN(ABS(bustr(i,j)),
     &                   ABS(u(i,j,1,nrhs))*cff1)
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          cff1=cff*0.5*(Hz(i,j-1,1)+Hz(i,j,1))
          bvstr(i,j)=SIGN(1.D0, bvstr(i,j))*
     &               MIN(ABS(bvstr(i,j)),
     &                   ABS(v(i,j,1,nrhs))*cff1)
        enddo
      enddo
#   endif
#  endif /* BSTRESS_FAST */

#  ifdef SEDIMENT
!
!  Compute skin stress components at rho points
!  (for the sediment model)
!
        do j=Jstr,Jend
          do i=Istr,Iend
#   ifdef DUNE
            znots=Zobt
#   else
            znots=Ssize(i,j)/12.
#   endif
            cff=vonKar/LOG((z_r(i,j,1)-z_w(i,j,0))/znots)
            cff=MIN(Cdb_max,MAX(Cdb_min,cff*cff))
            Ucur=0.5*(u(i,j,1,nrhs)+u(i+1,j,1,nrhs))
            Vcur=0.5*(v(i,j,1,nrhs)+v(i,j+1,1,nrhs))
#   if defined GLS_MIXING && defined SANDBAR
            Utur=15.*SQRT(trb(i,j,1,nrhs,itke))
     &                           *exp(-0.01*(xr(i,j)-142.)**2)
#   else
            Utur=0.
#   endif
            Umag=SQRT(Ucur*Ucur+Vcur*Vcur+Utur*Utur)
            bustrw(i,j)=cff*Umag*Ucur
            bvstrw(i,j)=cff*Umag*Vcur
          enddo
        enddo
#  endif /* SEDIMENT */
!
! Set boundary conditions
!
#  ifndef EW_PERIODIC
      IF (EASTERN_EDGE) THEN
        DO j=Jstr,Jend
          bustr(Iend+1,j)=bustr(Iend,j)
        END DO
        DO j=JstrV,Jend
          bvstr(Iend+1,j)=bvstr(Iend,j)
        END DO
#   ifdef SEDIMENT
        DO j=Jstr,Jend
          bustrw(Iend+1,j)=bustrw(Iend,j)
          bvstrw(Iend+1,j)=bvstrw(Iend,j)
        END DO
#   endif
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=Jstr,Jend
          bustr(IstrU-1,j)=bustr(IstrU,j)
        END DO
        DO j=JstrV,Jend
          bvstr(Istr-1,j)=bvstr(Istr,j)
        END DO
#   ifdef SEDIMENT
        DO j=Jstr,Jend
          bustrw(Istr-1,j)=bustrw(Istr,j)
          bvstrw(Istr-1,j)=bvstrw(Istr,j)
        END DO
#   endif
      END IF
#  endif
#  ifndef NS_PERIODIC
      IF (NORTHERN_EDGE) THEN
        DO i=IstrU,Iend
          bustr(i,Jend+1) =bustr(i,Jend)
        END DO
        DO i=Istr,Iend
          bvstr(i,Jend+1) =bvstr(i,Jend)
        END DO
#   ifdef SEDIMENT
        DO i=Istr,Iend
          bustrw(i,Jend+1) =bustrw(i,Jend)
          bvstrw(i,Jend+1) =bvstrw(i,Jend)
        END DO
#   endif
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=IstrU,Iend
          bustr(i,Jstr-1)=bustr(i,Jstr)
        END DO
        DO i=Istr,Iend
          bvstr(i,JstrV-1)=bvstr(i,JstrV)
        END DO
#   ifdef SEDIMENT
        DO i=Istr,Iend
          bustrw(i,Jstr-1)=bustrw(i,Jstr)
          bvstrw(i,Jstr-1)=bvstrw(i,Jstr)
        END DO
#   endif
      END IF
#  endif
#  if !defined EW_PERIODIC && !defined NS_PERIODIC
      ISTRm=Istr-1
      ISTRp=Istr+1
      IENDp=Iend+1
      JSTRm=Jstr-1
      JSTRp=Jstr+1
      JENDp=Jend+1
      IF (SOUTHERN_EDGE.and.WESTERN_EDGE) THEN
        bustr(Istr,JSTRm)  =0.5*(bustr(ISTRp,JSTRm)+bustr(Istr,Jstr))
        bvstr(ISTRm,Jstr)  =0.5*(bvstr(Istr,Jstr)+bvstr(ISTRm,JSTRp))
#   ifdef SEDIMENT
        bustrw(ISTRm,JSTRm)=0.5*(bustrw(Istr,JSTRm)+bustrw(ISTRm,Jstr))
        bvstrw(ISTRm,JSTRm)=0.5*(bvstrw(Istr,JSTRm)+bvstrw(ISTRm,Jstr))
#   endif
      ENDIF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        bustr(IENDp,JSTRm) =0.5*(bustr(IENDp,Jstr)+bustr(Iend,JSTRm))
        bvstr(IENDp,Jstr)  =0.5*(bvstr(IENDp,JSTRp)+bvstr(Iend,Jstr))
#   ifdef SEDIMENT
        bustrw(IENDp,JSTRm)=0.5*(bustrw(IENDp,Jstr)+bustrw(Iend,JSTRm))
        bvstrw(IENDp,JSTRm)=0.5*(bvstrw(IENDp,Jstr)+bvstrw(Iend,JSTRm))
#   endif
      ENDIF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        bustr(Istr,JENDp)  =0.5*(bustr(Istr,Jend)+bustr(ISTRp,JENDp))
        bvstr(ISTRm,JENDp) =0.5*(bvstr(ISTRm,Jend)+bvstr(Istr,JENDp))
#   ifdef SEDIMENT
        bustrw(ISTRm,JENDp)=0.5*(bustrw(ISTRm,Jend)+bustrw(Istr,JENDp))
        bvstrw(ISTRm,JENDp)=0.5*(bvstrw(ISTRm,Jend)+bvstrw(Istr,JENDp))
#   endif
      ENDIF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        bustr(IENDp,JENDp) =0.5*(bustr(IENDp,Jend)+bustr(Iend,JENDp))
        bvstr(IENDp,JENDp) =0.5*(bvstr(IENDp,Jend)+bvstr(Iend,JENDp))
#   ifdef SEDIMENT
        bustrw(IENDp,JENDp)=0.5*(bustrw(IENDp,Jend)+bustrw(Iend,JENDp))
        bvstrw(IENDp,JENDp)=0.5*(bvstrw(IENDp,Jend)+bvstrw(Iend,JENDp))
#   endif
      ENDIF
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend, bustr)
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend, bvstr)
#   ifdef SEDIMENT
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend, bustrw)
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend, bvstrw)
#   endif
#  endif

# endif /* SOLVE3D */
#endif /* BBL */

#if defined BBL && defined AGRIF
      ENDIF
#endif

#ifdef SOLVE3D
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,sustr(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,svstr(START_2D_ARRAY))
#  ifdef TEMPERATURE
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                                      stflx(START_2D_ARRAY,itemp))
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,srflx(START_2D_ARRAY))
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                                      btflx(START_2D_ARRAY,itemp))
#  endif
#  ifdef SALINITY
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                                     stflx(START_2D_ARRAY,isalt))
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                                     btflx(START_2D_ARRAY,isalt))
#  endif
# endif
#endif
#ifdef PSOURCE_NCFILE
      call set_psource_tile(Istr,Iend,Jstr,Jend)
      do is=1,Nsrc
        Qbar(is)=qbardir(is)*Qbar(is)
      enddo
# ifdef PSOURCE_NCFILE_TS
      call set_psource_ts_tile(Istr,Iend,Jstr,Jend)
# endif
#endif
#if defined PSOURCE && defined ANA_PSOURCE
C$OMP BARRIER
C$OMP MASTER
      call ana_psource_tile (Istr,Iend,Jstr,Jend)
C$OMP END MASTER
#endif
#ifndef OW_COUPLING
# if (defined MUSTANG || defined BBL || defined MRL_WCI) && defined WAVE_OFFLINE
      call set_wwave_tile(Istr,Iend,Jstr,Jend)
# endif
#endif

#if defined RVTK_DEBUG_ADVANCED && defined SOLVE3D
C$OMP BARRIER
C$OMP MASTER
      call check_tab3d(u(:,:,:,nrhs),'#end u in set_vbc','u')
      call check_tab3d(v(:,:,:,nrhs),'#end v in set_vbc','v')

# ifdef TEMPERATURE
      call check_tab3d(t(:,:,:,nrhs,itemp),'#end temp in set_vbc','r')
# endif
# ifdef SALINITY
      call check_tab3d(t(:,:,:,nrhs,isalt),'#end salt in set_vbc','r')
# endif
C$OMP END MASTER
#endif
C$OMP BARRIER
#if defined RVTK_DEBUG_ADVANCED && defined BULK_FLUX
C$OMP MASTER
      call check_tab2d(uwndg(:,:,1),'uwndg1 set_vbc','u')
      call check_tab2d(uwndg(:,:,2),'uwndg2 set_vbc','u')
      call check_tab2d(vwndg(:,:,1),'vwndg1 set_vbc','v')
      call check_tab2d(vwndg(:,:,2),'vwndg2 set_vbc','v')
      call check_tab2d(tair(:,:),'tair set_vbc','r')
      call check_tab2d(rhum(:,:),'rhum set_vbc','r')
      call check_tab2d(prate(:,:),'prate set_vbc','r')
      call check_tab2d(radlw(:,:),'radlw set_vbc','r')
      call check_tab2d(radsw(:,:),'radsw set_vbc','r')
      call check_tab2d(wspd(:,:),'wspd set_vbc','r')
      call check_tab2d(uwnd(:,:),'uwnd set_vbc','u')
      call check_tab2d(vwnd(:,:),'vwnd set_vbc','v')
C$OMP END MASTER
#endif

      return
      end
