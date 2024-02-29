! $Id: online_set_bulk.F 1458 2014-02-03 15:01:25Z gcambon $
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
!  This is the "online_set_bulk.F" script
!------------------------------------------------------------------------------
!  This file contains the subfunctions enabling the online extraction of the 
!  forcing from a NCEP/CFSR dataset. A spatial and time interpolation are applied 
!  to the extracted data in order to adapt these to the considered simulation 
!  domain as well as the associated MPI/OPENMP discretisation (if defined MPI/
!  OPENMP).
!------------------------------------------------------------------------------
#include "cppdefs.h"

#if defined BULK_FLUX && defined ONLINE 

!*******************************************************************************
      subroutine set_bulk_tile_online (Istr,Iend,Jstr,Jend)
!------------------------------------------------------------------------------
!  This subfunction enables the time interpolation of the roms bulk forcings 
!  (tair, rhum, prate, radlw, radsw, uwnd, vwnd) on the simulation domain
!  for the ongoing tile.
!------------------------------------------------------------------------------
!  The steps in this interpolation is:
!  - Determination of forcings (tair, rhum, prate, radlw, radsw, uwnd, 
!  vwnd) at the first time step of the simulation from the known (tairg, rhumg, 
!  prateg, radlwg, radswg, uwndg, vwndg) fields.
!  - Time interpolation of the forcings (tair, rhum, prate, radlw, radsw,  
!  uwnd, vwnd) at cff from two known spatially interpolated fields (tairg,rhumg,
!  prateg, radlwg, radswg, uwndg, vwndg) at cff1 and cff2.
!
!  In this subroutine the aformentioned steps are completed for every bulk 
!  forcing iterating on the various blkvar_id.
      implicit none
# include "param.h"
# include "forces.h"
# include "scalars.h"
# include "grid.h"
# include "online.h"
# ifdef CFB_WIND_TRA
#  include "ocean3d.h"
#  include "params_bulk.h"
# endif
      integer Istr,Iend,Jstr,Jend, i,j, it1,it2, blkvar_id
      real cff,cff1,cff2, cff3,cff4
# ifdef SALINITY
      real cff5,cff6
# endif      
# ifdef CFB_WIND_TRA
      real cffu1,cffu2,cffv1,cffv2
      real uwnd_r, vwnd_r
# endif
!     
# include "compute_extended_bounds.h"
!
       
      !===== == === === ======= ========= 
      ! Loop on all the forcing datasets:
      !===== == === === ======= =========  
      blkvar_id=0
  10  blkvar_id = blkvar_id+1
!
! Do not process #7 :  upward longwave is obsolete
!
        if (blkvar_id.eq.7) goto 10
#  if defined ERA_ECMWF || defined AROME
!
! In the case of ERA the net short wave is used
! Net Short wave - 5 downward short wave not used
!
        if (blkvar_id.eq.4) goto 10
#  endif

      if(blkvar_id.gt.nblkvrs) return

      it1=3-itbulkO(blkvar_id)
      it2=itbulkO(blkvar_id)
      ! --------------------------
      ! Times of the interpolation
      ! --------------------------
      cff=time+0.5*dt
      cff1=bulk_timeO(it2,blkvar_id)-cff
      cff2=cff-bulk_timeO(it1,blkvar_id)


      !===== ============= ======= === ==== ====== 
      ! Time interpolation between two time steps:
      !===== ============= ======= === ==== ====== 
      ! The roms bulk forcings at cff are interpolated between two time steps 
      ! cff1 and cff2
      if (cff1.ge.0. .and. cff2.ge.0.) then
        if (ZEROTH_TILE .and. cff1.lt.dt) synchro_flag=.TRUE.
                                    !note cff order maters
        cff=srf_scale/(cff1+cff2)
        cff3=cff1*cff
        cff4=cff2*cff
# ifdef SALINITY
        cff=stf_scale(isalt)/(cff1+cff2)
        cff5=cff1*cff
        cff6=cff2*cff
# endif
        cff=1./(cff1+cff2)
        cff1=cff1*cff
        cff2=cff2*cff
        if (blkvar_id.eq.1) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              tair(i,j)=cff1*tairg(i,j,it1)+cff2*tairg(i,j,it2)
            enddo
          enddo
        elseif (blkvar_id.eq.2) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              rhum(i,j)=cff1*rhumg(i,j,it1)+cff2*rhumg(i,j,it2)
            enddo
          enddo
# ifdef SALINITY
        elseif (blkvar_id.eq.3) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              prate(i,j)=cff5*prateg(i,j,it1)+cff6*prateg(i,j,it2)
            enddo
          enddo
# endif
        elseif (blkvar_id.eq.5) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              radsw(i,j)=cff3*radswg(i,j,it1)+cff4*radswg(i,j,it2)
              srflx(i,j)=radsw(i,j)
            enddo
          enddo
        elseif (blkvar_id.eq.6) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              radlw(i,j)=cff3*radlwg(i,j,it1)+cff4*radlwg(i,j,it2)
            enddo
          enddo
        elseif (blkvar_id.eq.9) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              uwnd(i,j)=cff1*uwndg(i,j,it1)+cff2*uwndg(i,j,it2)
              vwnd(i,j)=cff1*vwndg(i,j,it1)+cff2*vwndg(i,j,it2)
              wspd(i,j)=cff1*wspdg(i,j,it1)+cff2*wspdg(i,j,it2)
            enddo
          enddo
# ifdef CFB_WIND_TRA
          cff      = 1.-swparam    ! current-wind coupling parameter: Ua => Ua-(1-sw)Uo
          
          do j=JstrR,min(JendR,Mm)
            do i=IstrR,min(IendR,Lm)
              cffu1=cff1*uwndg(i  ,j,it1)+cff2*uwndg(i  ,j,it2)
              cffu2=cff1*uwndg(i+1,j,it1)+cff2*uwndg(i+1,j,it2)
              uwnd_r   = 0.5*( cffu1 + cffu2 
     &                 - cff*( u(i+1,j,N,nrhs)+u(i,j,N,nrhs))
     &                   )

              cffv1=cff1*vwndg(i,j  ,it1)+cff2*vwndg(i,j  ,it2)
              cffv2=cff1*vwndg(i,j+1,it1)+cff2*vwndg(i,j+1,it2)
              vwnd_r   = 0.5*( cffv1 + cffv2
     &                 - cff*( v(i,j+1,N,nrhs)+v(i,j,N,nrhs))
     &                   )

              wspd_cfb(i,j) = SQRT( uwnd_r*uwnd_r+vwnd_r*vwnd_r )
            enddo
          enddo
#  ifndef EW_PERIODIC
            if (WESTERN_EDGE) then
              do j=Jstr,Jend
                wspd_cfb(Istr-1,j)=wspd_cfb(Istr,j)
              enddo
            endif
            if (EASTERN_EDGE) then
              do j=Jstr,Jend
                wspd_cfb(Iend+1,j)=wspd_cfb(Iend,j)
              enddo
            endif
#  endif
#  ifndef NS_PERIODIC
            if (SOUTHERN_EDGE) then
              do i=Istr,Iend
                wspd_cfb(i,Jstr-1)=wspd_cfb(i,Jstr)
              enddo
            endif
            if (NORTHERN_EDGE) then
              do i=Istr,Iend
                wspd_cfb(i,Jend+1)=wspd_cfb(i,Jend)
              enddo
            endif
#   ifndef EW_PERIODIC
            if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
              wspd_cfb(Istr-1,Jstr-1)=wspd_cfb(Istr,Jstr)
            endif
            if (WESTERN_EDGE .and. NORTHERN_EDGE) then
              wspd_cfb(Istr-1,Jend+1)=wspd_cfb(Istr,Jend)
            endif
            if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
              wspd_cfb(Iend+1,Jstr-1)=wspd_cfb(Iend,Jstr)
            endif
            if (EASTERN_EDGE .and. NORTHERN_EDGE) then
              wspd_cfb(Iend+1,Jend+1)=wspd_cfb(Iend,Jend)
            endif
#   endif
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
            call exchange_r2d_tile(Istr,Iend,Jstr,Jend,
     &                             wspd_cfb(START_2D_ARRAY))
#  endif
# endif
# ifdef READ_PATM
        elseif (blkvar_id.eq.10) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              patm2d(i,j)=cff1*patmg(i,j,it1)+cff2*patmg(i,j,it2)
            enddo
          enddo
# endif
        endif

      !====== == ==== ============== 
      ! Error in time interpolation:
      !====== == ==== ============== 
      else
        if (ZEROTH_TILE) then
          write(stdout,1) 'bulk_timeO',tdays,bulk_timeO(it2,blkvar_id)*sec2day
  1       format(/,' SET_BULK - current model time exceeds ending',
     &           1x,'value for variable: ',a,/,11x,'TDAYS = ',g12.4,
     &           2x,'TEND = ',g12.4)
          may_day_flag=2
        endif
      endif

      goto 10

      return
      end

#endif /* BULK_FLUX && ONLINE */
