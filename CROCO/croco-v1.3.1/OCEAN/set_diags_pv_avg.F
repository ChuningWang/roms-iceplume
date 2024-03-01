!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
! 
! ROMS_AGRIF website : http://www.romsagrif.org
!======================================================================
!
#include "cppdefs.h"
#if (defined DIAGNOSTICS_PV && defined AVERAGES)

      subroutine set_diags_pv_avg(tile)

!      USE param

      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call set_diags_pv_avg_tile(Istr,Iend,Jstr,Jend)
      return
      end

      subroutine set_diags_pv_avg_tile(Istr,Iend,Jstr,Jend)
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
      integer Istr,Iend,Jstr,Jend, i,j, itrc, iflux, k, ilc
      real    cff, cff1
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "averages.h"
# include "diags_pv.h"
!
# include "compute_auxiliary_bounds.h"
!
      ilc=1+iic-ntstart   ! number of time step since restart
!
! calculation of averaged fluxes will only be performed if ntsdiags_pv_avg
! is a positive number
!
      if (ilc.gt.ntsdiags_pv_avg) then

        if (mod(ilc-ntsdiags_pv_avg,nwrtdiags_pv_avg).eq.1) then
          cff =1.0
          cff1=0.0
          if (ZEROTH_TILE) then
            timediags_pv_avg=time
!             MPI_master_only write(*,*) 'started averaging pv',iic,
!     &                             ntsdiags_pv_avg,nwrtdiags_pv_avg
          endif
        elseif (mod(ilc-ntsdiags_pv_avg,nwrtdiags_pv_avg).gt.1) then
          cff =1.0
          cff1=1.0
          if (ZEROTH_TILE) then
           timediags_pv_avg=timediags_pv_avg+time
          endif
        elseif (mod(ilc-ntsdiags_pv_avg,nwrtdiags_pv_avg).eq.0) then
          cff =1.0/float(nwrtdiags_pv_avg)
          cff1=1.0
          if (ZEROTH_TILE) then
            timediags_pv_avg=cff*(timediags_pv_avg+time)
!            MPI_master_only write(*,*) 'finished averaging pv',
!     &                 iic,ntsdiags_pv_avg,nwrtdiags_pv_avg
          endif
        endif

        do iflux=1,2
          do k=1,N
            do j=JstrR,JendR
              do i=IstrR,IendR 
                Mrhs_avg(i,j,k,iflux) = cff *
     &                     ( cff1*Mrhs_avg(i,j,k,iflux) +
     &                            Mrhs(i,j,k,iflux) )
              enddo
            enddo
          enddo
        enddo

        do itrc=1,NTA
          do k=1,N
            do j=JstrR,JendR
              do i=IstrR,IendR 
                Trhs_avg(i,j,k,itrc) = cff *
     &                     ( cff1*Trhs_avg(i,j,k,itrc) +
     &                            Trhs(i,j,k,itrc) )
              enddo
            enddo
          enddo
        enddo

# if defined DIAGNOSTICS_PV_FULL

          do k=0,N
            do j=JstrR,JendR
              do i=IstrR,IendR 
                pv_avg(i,j,k) = cff *
     &                     ( cff1*pv_avg(i,j,k) +
     &                            pv(i,j,k) )
                pvd_avg(i,j,k) = cff *
     &                     ( cff1*pvd_avg(i,j,k) +
     &                            pvd(i,j,k) )
              enddo
            enddo
          enddo
# endif

      endif

      return
      end
#else /* DIAGNOSTICS_PV && AVERAGES */
      subroutine set_diags_pv_avg_empty
      end
#endif /* DIAGNOSTICS_PV && AVERAGES */
