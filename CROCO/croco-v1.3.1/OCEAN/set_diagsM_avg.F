! $Id: set_diagsM_avg.F 1458 2014-02-03 15:01:25Z gcambon $
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
#if (defined DIAGNOSTICS_UV && defined AVERAGES)

      subroutine set_diagsM_avg(tile)

!      USE param

      implicit none
      integer tile
# include "param.h"
# ifdef SOLVE3D
#  include "work.h"
#  include "ncscrum.h"
# endif
# include "compute_tile_bounds.h"
# ifdef TENDENCY
      call Wvlcty (tile, workr)
!!!      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
!!!     &                        workr(START_2D_ARRAY,1))
# endif
      call set_diagsM_avg_tile(Istr,Iend,Jstr,Jend)
      return
      end

      subroutine set_diagsM_avg_tile(Istr,Iend,Jstr,Jend)
!
! Compute time-averaged fields within a tile.
! ------- ------------- ------ ------ - -----
! Because of syncronization issues, the delayed mode averaging
! procedure is used. This procedure implies that all fields to be
! averaged are sampled during the next time step, rather than at
! the end of the time step when they were computed.
!
! Although this algorithm results in somewhat awkward controlling
! logic it has the advantage that all fields to be sampled
! correspond to exactly the same time, which is time step "n".
! Particularly, this is done this way because vertical velocity
! corresponding to the newly computed horizontal velocities
! becomes available only during the following time step.
! The same applies to the density field.
!
! The algorithm consists of three logical blocks: (1) initialization
! of the averages arrays: when mod(ilc-1,navg).eq.1 the target arrays
! are set to the first contribution; (2) accumulation of averaged
! data, when mod(ilc-1,navg).gt.1; and (3) adding the last
! contribution and scaling.
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j, iflux, k, ilc
      real    cff, cff1
      real    cffux, cffuy, cffvx, cffvy, cffM
      real    cffwx, cffwy, cffuz, cffvz

# include "param.h"
#  include "work.h"
# include "scalars.h"
# include "ncscrum.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "averages.h"
# include "diagnostics.h"
!
# include "compute_auxiliary_bounds.h"
!
      ilc=1+iic-ntstart   ! number of time step since restart
!
! calculation of averaged fluxes will only be performed if ntsdiaM_avg
! is a positive number
!
      if (ilc.gt.ntsdiaM_avg) then

        if (mod(ilc-ntsdiaM_avg,navg).eq.1) then
          cff =1.0
          cff1=0.0
          if (ZEROTH_TILE) then
            timediaM_avg=time
!            MPI_master_only write(*,*) 'started averaging M',
!     &                          iic,ntsdiaM_avg,nwrtdiaM_avg
          endif
        elseif (mod(ilc-ntsdiaM_avg,navg).gt.1) then
          cff =1.0
          cff1=1.0
          if (ZEROTH_TILE) then
           timediaM_avg=timediaM_avg+time
          endif
        elseif (mod(ilc-ntsdiaM_avg,navg).eq.0) then
          cff =1.0/float(nwrtdiaM_avg)
          cff1=1.0
          if (ZEROTH_TILE) then
            timediaM_avg=cff*(timediaM_avg+time)
!            MPI_master_only write(*,*) 'finished averaging M',
!     &                    iic,ntsdiaM_avg,nwrtdiaM_avg
          endif
        endif

        do iflux=1,2
          do k=1,N
            do j=JstrR,JendR
              do i=IstrR,IendR
# if defined TENDENCY
! compute u,v,w gradients on rho grid

            cffux = (u(i+1,j,k,nstp) -  u(i,j,k,nstp))
     &            * 0.5 * (pm_u(i+1,j) + pm_u(i,j))

            cffuy = 0.25 *(u(i,j+1,k,nstp) -  u(i,j,k,nstp))
     &            * 0.5 * (pn_u(i,j+1) + pn_u(i,j))
     &            + 0.25 * (u(i,j,k,nstp) -  u(i,j-1,k,nstp))
     &            * 0.5 * (pn_u(i,j) + pn_u(i,j-1))
     &           + 0.25 *(u(i+1,j+1,k,nstp) -  u(i+1,j,k,nstp))
     &            * 0.5 * (pn_u(i+1,j+1) + pn_u(i+1,j))
     &          + 0.25 * (u(i+1,j,k,nstp) -  u(i+1,j-1,k,nstp))
     &            * 0.5 * (pn_u(i+1,j) + pn_u(i+1,j-1))

            cffvx = 0.25 *(v(i+1,j,k,nstp) -  v(i,j,k,nstp))
     &            * 0.5 * (pm_v(i+1,j) + pm_v(i,j))
     &            + 0.25 * (v(i,j,k,nstp) -  v(i-1,j,k,nstp))
     &            * 0.5 * (pm_v(i,j) + pm_v(i-1,j))
     &            + 0.25 *(v(i+1,j+1,k,nstp) -  v(i,j+1,k,nstp))
     &            * 0.5 * (pm_v(i+1,j+1) + pm_v(i,j+1))
     &            + 0.25 * (v(i,j+1,k,nstp) -  v(i-1,j+1,k,nstp))
     &            * 0.5 * (pm_v(i,j+1) + pm_v(i-1,j+1))

            cffvy = (v(i,j+1,k,nstp) - v(i,j,k,nstp))
     &            * 0.5 * (pn_v(i,j+1) + pn_v(i,j))

            cffwx = 0.5 *(workr(i+1,j,k) -  workr(i,j,k))
     &            * 0.5 * (pm(i+1,j) + pm(i,j))
     &            + 0.5 * (workr(i,j,k) -  workr(i-1,j,k))
     &            * 0.5 * (pm(i,j) + pm(i-1,j))

            cffwy = 0.5 *(workr(i,j+1,k) -  workr(i,j,k))
     &            * 0.5 * (pn(i,j+1) + pn(i,j))
     &            + 0.5 * (workr(i,j,k) -  workr(i,j-1,k))
     &            * 0.5 * (pn(i,j) + pn(i,j-1))


               if (iflux.eq.1) then

!               var = -1*(ux * ux * ux + uy * uy * ux + ux * vx * uy + uy * vy * uy)
!                     -1*(vx * ux * vx + vy * uy * vx + vx * vx * vy + vy * vy * vy)


                cffM = -1*(cffux * cffux * cffux
     &               + cffuy * cffuy * cffux)

!                cffM = ( MXadv(i+1,j,k,iflux)/pm_u(i+1,j)
!     &               - MXadv(i,j,k,iflux)/pm_u(i,j) ) * cffux
!     &               + 0.5*(MXadv(i+1,j+1,k,iflux)/pm_u(i+1,j+1)
!     &               -  MXadv(i+1,j,k,iflux)  /pm_u(i+1,j)
!     &               +  MXadv(i,j+1,k,iflux) /pm_u(i,j+1)
!     &               -  MXadv(i,j,k,iflux)  /pm_u(i,j) ) * cffuy

                MXadv_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MXadv_avg(i,j,k,iflux) +
     &                            cffM )

                cffM = -1*(cffux * cffvx * cffuy
     &               + cffuy * cffvy * cffuy)

!                cffM = ( MYadv(i+1,j,k,iflux)/pm_u(i+1,j)
!     &               - MYadv(i,j,k,iflux)/pm_u(i,j) ) * cffux
!     &               + 0.5*(MYadv(i+1,j+1,k,iflux)/pm_u(i+1,j+1)
!     &               -  MYadv(i+1,j,k,iflux)  /pm_u(i+1,j)
!     &               +  MYadv(i,j+1,k,iflux) /pm_u(i,j+1)
!     &               -  MYadv(i,j,k,iflux)  /pm_u(i,j) ) * cffuy

                MYadv_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MYadv_avg(i,j,k,iflux) +
     &                            cffM )

!                cffM = ( MVadv(i+1,j,k,iflux)/pm_u(i+1,j)
!     &               - MVadv(i,j,k,iflux)/pm_u(i,j) ) * cffux
!     &               + 0.5*(MVadv(i+1,j+1,k,iflux)/pm_u(i+1,j+1)
!     &               -  MVadv(i+1,j,k,iflux)  /pm_u(i+1,j)
!     &               +  MVadv(i,j+1,k,iflux) /pm_u(i,j+1)
!     &               -  MVadv(i,j,k,iflux)  /pm_u(i,j) ) * cffuy


                if (k.eq.N) then

                cffuz = 0.5*(u(i,j,k,nstp) -  u(i,j,k-1,nstp)
     &                 +   u(i+1,j,k,nstp) -  u(i+1,j,k-1,nstp))
     &                  / (z_r(i,j,k) - z_r(i,j,k-1))

                elseif (k.eq.1) then

                cffuz = 0.5*(u(i,j,k+1,nstp) -  u(i,j,k,nstp)
     &                 +   u(i+1,j,k+1,nstp) -  u(i+1,j,k,nstp))
     &                  / (z_r(i,j,k+1) - z_r(i,j,k))
                else

                cffuz = 0.5*(u(i,j,k+1,nstp) -  u(i,j,k-1,nstp)
     &                 +   u(i+1,j,k+1,nstp) -  u(i+1,j,k-1,nstp))
     &                  / (z_r(i,j,k+1) - z_r(i,j,k-1))
                endif

                cffM = -1*(cffuz * cffwx * cffux
     &               + cffuz * cffwy * cffuy)

                MVadv_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MVadv_avg(i,j,k,iflux) +
     &                            cffM )

            cffM = ( MCor(i+1,j,k,iflux)-MCor(i,j,k,iflux) )
     &             * 0.5 * (pm_u(i+1,j) + pm_u(i,j)) * cffux
     & + (0.25 *(MCor(i,j+1,k,iflux)-MCor(i,j,k,iflux))
     &               * 0.5 * (pn_u(i,j+1) + pn_u(i,j))
     & + 0.25 *(MCor(i,j,k,iflux)-MCor(i,j-1,k,iflux))
     &               * 0.5 * (pn_u(i,j) + pn_u(i,j-1))
     & + 0.25 *(MCor(i+1,j+1,k,iflux)-MCor(i+1,j,k,iflux))
     &               * 0.5 * (pn_u(i+1,j+1) + pn_u(i+1,j))
     & + 0.25 *(MCor(i+1,j,k,iflux)-MCor(i+1,j-1,k,iflux))
     &               * 0.5 * (pn_u(i+1,j) + pn_u(i+1,j-1))) * cffuy

                MCor_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MCor_avg(i,j,k,iflux) +
     &                            cffM )

            cffM = ( MPrsgrd(i+1,j,k,iflux)-MPrsgrd(i,j,k,iflux) )
     &             * 0.5 * (pm_u(i+1,j) + pm_u(i,j)) * cffux
     & + (0.25 *(MPrsgrd(i,j+1,k,iflux)-MPrsgrd(i,j,k,iflux))
     &               * 0.5 * (pn_u(i,j+1) + pn_u(i,j))
     & + 0.25 *(MPrsgrd(i,j,k,iflux)-MPrsgrd(i,j-1,k,iflux))
     &               * 0.5 * (pn_u(i,j) + pn_u(i,j-1))
     & + 0.25 *(MPrsgrd(i+1,j+1,k,iflux)-MPrsgrd(i+1,j,k,iflux))
     &               * 0.5 * (pn_u(i+1,j+1) + pn_u(i+1,j))
     & + 0.25 *(MPrsgrd(i+1,j,k,iflux)-MPrsgrd(i+1,j-1,k,iflux))
     &               * 0.5 * (pn_u(i+1,j) + pn_u(i+1,j-1))) * cffuy

                MPrsgrd_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MPrsgrd_avg(i,j,k,iflux) +
     &                            cffM )

       cffM = ( MHmix(i+1,j,k,iflux,nstp)-MHmix(i,j,k,iflux,nstp) )
     &             * 0.5 * (pm_u(i+1,j) + pm_u(i,j)) * cffux
     & + (0.25 *(MHmix(i,j+1,k,iflux,nstp)-MHmix(i,j,k,iflux,nstp))
     &               * 0.5 * (pn_u(i,j+1) + pn_u(i,j))
     & + 0.25 *(MHmix(i,j,k,iflux,nstp)-MHmix(i,j-1,k,iflux,nstp))
     &               * 0.5 * (pn_u(i,j) + pn_u(i,j-1))
     & + 0.25 *(MHmix(i+1,j+1,k,iflux,nstp)-MHmix(i+1,j,k,iflux,nstp))
     &               * 0.5 * (pn_u(i+1,j+1) + pn_u(i+1,j))
     & + 0.25 *(MHmix(i+1,j,k,iflux,nstp)-MHmix(i+1,j-1,k,iflux,nstp))
     &               * 0.5 * (pn_u(i+1,j) + pn_u(i+1,j-1))) * cffuy

                MHmix_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MHmix_avg(i,j,k,iflux) +
     &                            cffM )

            cffM = ( MHdiff(i+1,j,k,iflux)-MHdiff(i,j,k,iflux) )
     &             * 0.5 * (pm_u(i+1,j) + pm_u(i,j)) * cffux
     & + (0.25 *(MHdiff(i,j+1,k,iflux)-MHdiff(i,j,k,iflux))
     &               * 0.5 * (pn_u(i,j+1) + pn_u(i,j))
     & + 0.25 *(MHdiff(i,j,k,iflux)-MHdiff(i,j-1,k,iflux))
     &               * 0.5 * (pn_u(i,j) + pn_u(i,j-1))
     & + 0.25 *(MHdiff(i+1,j+1,k,iflux)-MHdiff(i+1,j,k,iflux))
     &               * 0.5 * (pn_u(i+1,j+1) + pn_u(i+1,j))
     & + 0.25 *(MHdiff(i+1,j,k,iflux)-MHdiff(i+1,j-1,k,iflux))
     &               * 0.5 * (pn_u(i+1,j) + pn_u(i+1,j-1))) * cffuy

                MHdiff_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MHdiff_avg(i,j,k,iflux) +
     &                            cffM )

            cffM = ( MVmix(i+1,j,k,iflux)-MVmix(i,j,k,iflux) )
     &             * 0.5 * (pm_u(i+1,j) + pm_u(i,j)) * cffux
     & + (0.25 *(MVmix(i,j+1,k,iflux)-MVmix(i,j,k,iflux))
     &               * 0.5 * (pn_u(i,j+1) + pn_u(i,j))
     & + 0.25 *(MVmix(i,j,k,iflux)-MVmix(i,j-1,k,iflux))
     &               * 0.5 * (pn_u(i,j) + pn_u(i,j-1))
     & + 0.25 *(MVmix(i+1,j+1,k,iflux)-MVmix(i+1,j,k,iflux))
     &               * 0.5 * (pn_u(i+1,j+1) + pn_u(i+1,j))
     & + 0.25 *(MVmix(i+1,j,k,iflux)-MVmix(i+1,j-1,k,iflux))
     &               * 0.5 * (pn_u(i+1,j) + pn_u(i+1,j-1))) * cffuy

                MVmix_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MVmix_avg(i,j,k,iflux) +
     &                            cffM )

            cffM = ( MVmix2(i+1,j,k,iflux)-MVmix2(i,j,k,iflux) )
     &             * 0.5 * (pm_u(i+1,j) + pm_u(i,j)) * cffux
     & + (0.25 *(MVmix2(i,j+1,k,iflux)-MVmix2(i,j,k,iflux))
     &               * 0.5 * (pn_u(i,j+1) + pn_u(i,j))
     & + 0.25 *(MVmix2(i,j,k,iflux)-MVmix2(i,j-1,k,iflux))
     &               * 0.5 * (pn_u(i,j) + pn_u(i,j-1))
     & + 0.25 *(MVmix2(i+1,j+1,k,iflux)-MVmix2(i+1,j,k,iflux))
     &               * 0.5 * (pn_u(i+1,j+1) + pn_u(i+1,j))
     & + 0.25 *(MVmix2(i+1,j,k,iflux)-MVmix2(i+1,j-1,k,iflux))
     &               * 0.5 * (pn_u(i+1,j) + pn_u(i+1,j-1))) * cffuy

                MVmix2_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MVmix2_avg(i,j,k,iflux) +
     &                            cffM )

            cffM = ( Mrate(i+1,j,k,iflux)-Mrate(i,j,k,iflux) )
     &             * 0.5 * (pm_u(i+1,j) + pm_u(i,j)) * cffux
     & + (0.25 *(Mrate(i,j+1,k,iflux)-Mrate(i,j,k,iflux))
     &               * 0.5 * (pn_u(i,j+1) + pn_u(i,j))
     & + 0.25 *(Mrate(i,j,k,iflux)-Mrate(i,j-1,k,iflux))
     &               * 0.5 * (pn_u(i,j) + pn_u(i,j-1))
     & + 0.25 *(Mrate(i+1,j+1,k,iflux)-Mrate(i+1,j,k,iflux))
     &               * 0.5 * (pn_u(i+1,j+1) + pn_u(i+1,j))
     & + 0.25 *(Mrate(i+1,j,k,iflux)-Mrate(i+1,j-1,k,iflux))
     &               * 0.5 * (pn_u(i+1,j) + pn_u(i+1,j-1))) * cffuy

                Mrate_avg(i,j,k,iflux) = cff *
     &                     ( cff1*Mrate_avg(i,j,k,iflux) +
     &                            cffM )

# if defined DIAGNOSTICS_BARO
            cffM = ( MBaro(i+1,j,k,iflux)-MBaro(i,j,k,iflux) )
     &             * 0.5 * (pm_u(i+1,j) + pm_u(i,j)) * cffux
     & + (0.25 *(MBaro(i,j+1,k,iflux)-MBaro(i,j,k,iflux))
     &               * 0.5 * (pn_u(i,j+1) + pn_u(i,j))
     & + 0.25 *(MBaro(i,j,k,iflux)-MBaro(i,j-1,k,iflux))
     &               * 0.5 * (pn_u(i,j) + pn_u(i,j-1))
     & + 0.25 *(MBaro(i+1,j+1,k,iflux)-MBaro(i+1,j,k,iflux))
     &               * 0.5 * (pn_u(i+1,j+1) + pn_u(i+1,j))
     & + 0.25 *(MBaro(i+1,j,k,iflux)-MBaro(i+1,j-1,k,iflux))
     &               * 0.5 * (pn_u(i+1,j) + pn_u(i+1,j-1))) * cffuy

                MBaro_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MBaro_avg(i,j,k,iflux) +
     &                            cffM )
# endif

# if defined M3FAST
            cffM = ( Mfast(i+1,j,k,iflux)-Mfast(i,j,k,iflux) )
     &             * 0.5 * (pm_u(i+1,j) + pm_u(i,j)) * cffux
     & + (0.25 *(Mfast(i,j+1,k,iflux)-Mfast(i,j,k,iflux))
     &               * 0.5 * (pn_u(i,j+1) + pn_u(i,j))
     & + 0.25 *(Mfast(i,j,k,iflux)-Mfast(i,j-1,k,iflux))
     &               * 0.5 * (pn_u(i,j) + pn_u(i,j-1))
     & + 0.25 *(Mfast(i+1,j+1,k,iflux)-Mfast(i+1,j,k,iflux))
     &               * 0.5 * (pn_u(i+1,j+1) + pn_u(i+1,j))
     & + 0.25 *(Mfast(i+1,j,k,iflux)-Mfast(i+1,j-1,k,iflux))
     &               * 0.5 * (pn_u(i+1,j) + pn_u(i+1,j-1))) * cffuy

                Mfast_avg(i,j,k,iflux) = cff *
     &                     ( cff1*Mfast_avg(i,j,k,iflux) +
     &                            cffM )
# endif

               elseif (iflux.eq.2) then


                cffM = -1*(cffvx * cffux * cffvx
     &               + cffvy * cffuy * cffvx)

!                cffM = ( MXadv(i,j+1,k,iflux)/pm_v(i+1,j)
!     &               - MXadv(i,j,k,iflux)/pm_v(i,j) ) * cffvy
!     &               + 0.5*(MXadv(i+1,j+1,k,iflux)/pm_v(i+1,j+1)
!     &               -  MXadv(i,j+1,k,iflux)  /pm_v(i,j+1)
!     &               +  MXadv(i+1,j,k,iflux) /pm_v(i+1,j)
!     &               -  MXadv(i,j,k,iflux)  /pm_v(i,j) ) * cffuy

                MXadv_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MXadv_avg(i,j,k,iflux) +
     &                            cffM )

                cffM = -1*(cffvx * cffvx * cffvy
     &               + cffvy * cffvy * cffvy)

!                cffM = ( MYadv(i,j+1,k,iflux)/pm_v(i+1,j)
!     &               - MYadv(i,j,k,iflux)/pm_v(i,j) ) * cffvy
!     &               + 0.5*(MYadv(i+1,j+1,k,iflux)/pm_v(i+1,j+1)
!     &               -  MYadv(i,j+1,k,iflux)  /pm_v(i,j+1)
!     &               +  MYadv(i+1,j,k,iflux) /pm_v(i+1,j)
!     &               -  MYadv(i,j,k,iflux)  /pm_v(i,j) ) * cffuy

                MYadv_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MYadv_avg(i,j,k,iflux) +
     &                            cffM )

                if (k.eq.N) then

                cffvz = 0.5*(v(i,j,k,nstp) -  v(i,j,k-1,nstp)
     &                 +   v(i,j+1,k,nstp) -  v(i,j+1,k-1,nstp))
     &                  / (z_r(i,j,k) - z_r(i,j,k-1))

                elseif (k.eq.1) then

                cffvz = 0.5*(v(i,j,k+1,nstp) -  v(i,j,k,nstp)
     &                 +   v(i,j+1,k+1,nstp) -  v(i,j+1,k,nstp))
     &                  / (z_r(i,j,k+1) - z_r(i,j,k))
                else

                cffvz = 0.5*(v(i,j,k+1,nstp) -  v(i,j,k-1,nstp)
     &                 +   v(i,j+1,k+1,nstp) -  v(i,j+1,k-1,nstp))
     &                  / (z_r(i,j,k+1) - z_r(i,j,k-1))
                endif

                cffM = -1*(cffvz * cffwx * cffvx
     &               + cffvz * cffwy * cffvy)

!                cffM = ( MVadv(i,j+1,k,iflux)/pm_v(i+1,j)
!     &               - MVadv(i,j,k,iflux)/pm_v(i,j) ) * cffvy
!     &               + 0.5*(MVadv(i+1,j+1,k,iflux)/pm_v(i+1,j+1)
!     &               -  MVadv(i,j+1,k,iflux)  /pm_v(i,j+1)
!     &               +  MVadv(i+1,j,k,iflux) /pm_v(i+1,j)
!     &               -  MVadv(i,j,k,iflux)  /pm_v(i,j) ) * cffuy

                MVadv_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MVadv_avg(i,j,k,iflux) +
     &                            cffM )



        cffM = (MCor(i,j+1,k,iflux) - MCor(i,j,k,iflux))
     &            * 0.5 * (pn_v(i,j+1) + pn_v(i,j)) * cffvy
     &  +(0.25*(MCor(i+1,j,k,iflux)-MCor(i,j,k,iflux))
     &            * 0.5 * (pm_v(i+1,j) + pm_v(i,j))
     &  + 0.25*(MCor(i,j,k,iflux)-MCor(i-1,j,k,iflux))
     &            * 0.5 * (pm_v(i,j) + pm_v(i-1,j))
     &  + 0.25*(MCor(i+1,j+1,k,iflux)-MCor(i,j+1,k,iflux))
     &            * 0.5 * (pm_v(i+1,j+1) + pm_v(i,j+1))
     &  + 0.25*(MCor(i,j+1,k,iflux)-MCor(i-1,j+1,k,iflux))
     &            * 0.5 * (pm_v(i,j+1) + pm_v(i-1,j+1))) * cffuy

                MCor_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MCor_avg(i,j,k,iflux) +
     &                            cffM )

        cffM = (MPrsgrd(i,j+1,k,iflux) - MPrsgrd(i,j,k,iflux))
     &            * 0.5 * (pn_v(i,j+1) + pn_v(i,j)) * cffvy
     &  +(0.25*(MPrsgrd(i+1,j,k,iflux)-MPrsgrd(i,j,k,iflux))
     &            * 0.5 * (pm_v(i+1,j) + pm_v(i,j))
     &  + 0.25*(MPrsgrd(i,j,k,iflux)-MPrsgrd(i-1,j,k,iflux))
     &            * 0.5 * (pm_v(i,j) + pm_v(i-1,j))
     &  + 0.25*(MPrsgrd(i+1,j+1,k,iflux)-MPrsgrd(i,j+1,k,iflux))
     &            * 0.5 * (pm_v(i+1,j+1) + pm_v(i,j+1))
     &  + 0.25*(MPrsgrd(i,j+1,k,iflux)-MPrsgrd(i-1,j+1,k,iflux))
     &            * 0.5 * (pm_v(i,j+1) + pm_v(i-1,j+1))) * cffuy

                MPrsgrd_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MPrsgrd_avg(i,j,k,iflux) +
     &                            cffM )


        cffM = (MVmix(i,j+1,k,iflux) - MVmix(i,j,k,iflux))
     &            * 0.5 * (pn_v(i,j+1) + pn_v(i,j)) * cffvy
     &  +(0.25*(MVmix(i+1,j,k,iflux)-MVmix(i,j,k,iflux))
     &            * 0.5 * (pm_v(i+1,j) + pm_v(i,j))
     &  + 0.25*(MVmix(i,j,k,iflux)-MVmix(i-1,j,k,iflux))
     &            * 0.5 * (pm_v(i,j) + pm_v(i-1,j))
     &  + 0.25*(MVmix(i+1,j+1,k,iflux)-MVmix(i,j+1,k,iflux))
     &            * 0.5 * (pm_v(i+1,j+1) + pm_v(i,j+1))
     &  + 0.25*(MVmix(i,j+1,k,iflux)-MVmix(i-1,j+1,k,iflux))
     &            * 0.5 * (pm_v(i,j+1) + pm_v(i-1,j+1))) * cffuy

                MVmix_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MVmix_avg(i,j,k,iflux) +
     &                            cffM )


        cffM = (MHmix(i,j+1,k,iflux,nstp) - MHmix(i,j,k,iflux,nstp))
     &            * 0.5 * (pn_v(i,j+1) + pn_v(i,j)) * cffvy
     &  +(0.25*(MHmix(i+1,j,k,iflux,nstp)-MHmix(i,j,k,iflux,nstp))
     &            * 0.5 * (pm_v(i+1,j) + pm_v(i,j))
     &  + 0.25*(MHmix(i,j,k,iflux,nstp)-MHmix(i-1,j,k,iflux,nstp))
     &            * 0.5 * (pm_v(i,j) + pm_v(i-1,j))
     &  + 0.25*(MHmix(i+1,j+1,k,iflux,nstp)-MHmix(i,j+1,k,iflux,nstp))
     &            * 0.5 * (pm_v(i+1,j+1) + pm_v(i,j+1))
     &  + 0.25*(MHmix(i,j+1,k,iflux,nstp)-MHmix(i-1,j+1,k,iflux,nstp))
     &            * 0.5 * (pm_v(i,j+1) + pm_v(i-1,j+1))) * cffuy


                MHmix_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MHmix_avg(i,j,k,iflux) +
     &                            cffM )


        cffM = (MHdiff(i,j+1,k,iflux) - MHdiff(i,j,k,iflux))
     &            * 0.5 * (pn_v(i,j+1) + pn_v(i,j)) * cffvy
     &  +(0.25*(MHdiff(i+1,j,k,iflux)-MHdiff(i,j,k,iflux))
     &            * 0.5 * (pm_v(i+1,j) + pm_v(i,j))
     &  + 0.25*(MHdiff(i,j,k,iflux)-MHdiff(i-1,j,k,iflux))
     &            * 0.5 * (pm_v(i,j) + pm_v(i-1,j))
     &  + 0.25*(MHdiff(i+1,j+1,k,iflux)-MHdiff(i,j+1,k,iflux))
     &            * 0.5 * (pm_v(i+1,j+1) + pm_v(i,j+1))
     &  + 0.25*(MHdiff(i,j+1,k,iflux)-MHdiff(i-1,j+1,k,iflux))
     &            * 0.5 * (pm_v(i,j+1) + pm_v(i-1,j+1))) * cffuy


                MHdiff_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MHdiff_avg(i,j,k,iflux) +
     &                            cffM )


        cffM = (MVmix2(i,j+1,k,iflux) - MVmix2(i,j,k,iflux))
     &            * 0.5 * (pn_v(i,j+1) + pn_v(i,j)) * cffvy
     &  +(0.25*(MVmix2(i+1,j,k,iflux)-MVmix2(i,j,k,iflux))
     &            * 0.5 * (pm_v(i+1,j) + pm_v(i,j))
     &  + 0.25*(MVmix2(i,j,k,iflux)-MVmix2(i-1,j,k,iflux))
     &            * 0.5 * (pm_v(i,j) + pm_v(i-1,j))
     &  + 0.25*(MVmix2(i+1,j+1,k,iflux)-MVmix2(i,j+1,k,iflux))
     &            * 0.5 * (pm_v(i+1,j+1) + pm_v(i,j+1))
     &  + 0.25*(MVmix2(i,j+1,k,iflux)-MVmix2(i-1,j+1,k,iflux))
     &            * 0.5 * (pm_v(i,j+1) + pm_v(i-1,j+1))) * cffuy


                MVmix2_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MVmix2_avg(i,j,k,iflux) +
     &                            cffM )


        cffM = (Mrate(i,j+1,k,iflux) - Mrate(i,j,k,iflux))
     &            * 0.5 * (pn_v(i,j+1) + pn_v(i,j)) * cffvy
     &  +(0.25*(Mrate(i+1,j,k,iflux)-Mrate(i,j,k,iflux))
     &            * 0.5 * (pm_v(i+1,j) + pm_v(i,j))
     &  + 0.25*(Mrate(i,j,k,iflux)-Mrate(i-1,j,k,iflux))
     &            * 0.5 * (pm_v(i,j) + pm_v(i-1,j))
     &  + 0.25*(Mrate(i+1,j+1,k,iflux)-Mrate(i,j+1,k,iflux))
     &            * 0.5 * (pm_v(i+1,j+1) + pm_v(i,j+1))
     &  + 0.25*(Mrate(i,j+1,k,iflux)-Mrate(i-1,j+1,k,iflux))
     &            * 0.5 * (pm_v(i,j+1) + pm_v(i-1,j+1))) * cffuy


                Mrate_avg(i,j,k,iflux) = cff *
     &                     ( cff1*Mrate_avg(i,j,k,iflux) +
     &                            cffM )

# if defined DIAGNOSTICS_BARO
        cffM = (MBaro(i,j+1,k,iflux) - MBaro(i,j,k,iflux))
     &            * 0.5 * (pn_v(i,j+1) + pn_v(i,j)) * cffvy
     &  +(0.25*(MBaro(i+1,j,k,iflux)-MBaro(i,j,k,iflux))
     &            * 0.5 * (pm_v(i+1,j) + pm_v(i,j))
     &  + 0.25*(MBaro(i,j,k,iflux)-MBaro(i-1,j,k,iflux))
     &            * 0.5 * (pm_v(i,j) + pm_v(i-1,j))
     &  + 0.25*(MBaro(i+1,j+1,k,iflux)-MBaro(i,j+1,k,iflux))
     &            * 0.5 * (pm_v(i+1,j+1) + pm_v(i,j+1))
     &  + 0.25*(MBaro(i,j+1,k,iflux)-MBaro(i-1,j+1,k,iflux))
     &            * 0.5 * (pm_v(i,j+1) + pm_v(i-1,j+1))) * cffuy

                MBaro_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MBaro_avg(i,j,k,iflux) +
     &                            cffM )
# endif

# if defined M3FAST
        cffM = (Mfast(i,j+1,k,iflux) - Mfast(i,j,k,iflux))
     &            * 0.5 * (pn_v(i,j+1) + pn_v(i,j)) * cffvy
     &  +(0.25*(Mfast(i+1,j,k,iflux)-Mfast(i,j,k,iflux))
     &            * 0.5 * (pm_v(i+1,j) + pm_v(i,j))
     &  + 0.25*(Mfast(i,j,k,iflux)-Mfast(i-1,j,k,iflux))
     &            * 0.5 * (pm_v(i,j) + pm_v(i-1,j))
     &  + 0.25*(Mfast(i+1,j+1,k,iflux)-Mfast(i,j+1,k,iflux))
     &            * 0.5 * (pm_v(i+1,j+1) + pm_v(i,j+1))
     &  + 0.25*(Mfast(i,j+1,k,iflux)-Mfast(i-1,j+1,k,iflux))
     &            * 0.5 * (pm_v(i,j+1) + pm_v(i-1,j+1))) * cffuy

                Mfast_avg(i,j,k,iflux) = cff *
     &                     ( cff1*Mfast_avg(i,j,k,iflux) +
     &                            cffM )
# endif

               endif


#else
                MXadv_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MXadv_avg(i,j,k,iflux) +
     &                            MXadv(i,j,k,iflux) )
                MYadv_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MYadv_avg(i,j,k,iflux) +
     &                            MYadv(i,j,k,iflux) )
                MVadv_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MVadv_avg(i,j,k,iflux) +
     &                            MVadv(i,j,k,iflux) )
                MCor_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MCor_avg(i,j,k,iflux) +
     &                            MCor(i,j,k,iflux) )
                MPrsgrd_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MPrsgrd_avg(i,j,k,iflux) +
     &                            MPrsgrd(i,j,k,iflux) )
                MHmix_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MHmix_avg(i,j,k,iflux) +
     &                            MHmix(i,j,k,iflux,nstp) )
                MHdiff_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MHdiff_avg(i,j,k,iflux) +
     &                            MHdiff(i,j,k,iflux) )
                MVmix_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MVmix_avg(i,j,k,iflux) +
     &                            MVmix(i,j,k,iflux) )
                MVmix2_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MVmix2_avg(i,j,k,iflux) +
     &                            MVmix2(i,j,k,iflux) )
                Mrate_avg(i,j,k,iflux) = cff *
     &                     ( cff1*Mrate_avg(i,j,k,iflux) + 
     &                            Mrate(i,j,k,iflux) )
# if defined DIAGNOSTICS_BARO
                MBaro_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MBaro_avg(i,j,k,iflux) +
     &                            MBaro(i,j,k,iflux) )
# endif
# if defined M3FAST
                Mfast_avg(i,j,k,iflux) = cff *
     &                     ( cff1*Mfast_avg(i,j,k,iflux) +
     &                            Mfast(i,j,k,iflux) )
# endif
# endif
# if defined MRL_WCI
                Mvf_avg(i,j,k,iflux)= cff *
     &                     ( cff1*Mvf_avg(i,j,k,iflux) + 
     &                            Mvf(i,j,k,iflux) )
                Mbrk_avg(i,j,k,iflux)= cff *
     &                     ( cff1*Mbrk_avg(i,j,k,iflux) + 
     &                            Mbrk(i,j,k,iflux) )
                MStCo_avg(i,j,k,iflux)= cff *
     &                     ( cff1*MStCo_avg(i,j,k,iflux) + 
     &                            MStCo(i,j,k,iflux) )
                MVvf_avg(i,j,k,iflux)= cff *
     &                     ( cff1*MVvf_avg(i,j,k,iflux) + 
     &                            MVvf(i,j,k,iflux) )
                MPrscrt_avg(i,j,k,iflux)= cff *
     &                     ( cff1*MPrscrt_avg(i,j,k,iflux) + 
     &                            MPrscrt(i,j,k,iflux) )
                Msbk_avg(i,j,k,iflux)= cff *
     &                     ( cff1*Msbk_avg(i,j,k,iflux) + 
     &                            Msbk(i,j,k,iflux) )
                Mbwf_avg(i,j,k,iflux)= cff *
     &                     ( cff1*Mbwf_avg(i,j,k,iflux) + 
     &                            Mbwf(i,j,k,iflux) )
                Mfrc_avg(i,j,k,iflux)= cff *
     &                     ( cff1*Mfrc_avg(i,j,k,iflux) + 
     &                            Mfrc(i,j,k,iflux) )
# endif
              enddo
            enddo
          enddo
        enddo

      endif

      return
      end
#else /* DIAGNOSTICS_UV && AVERAGES */
      subroutine set_diagsM_avg_empty
      end
#endif /* DIAGNOSTICS_UV && AVERAGES */
