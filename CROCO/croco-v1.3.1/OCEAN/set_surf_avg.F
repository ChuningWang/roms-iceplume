#include "cppdefs.h"
#if defined OUTPUTS_SURFACE && defined AVERAGES && ! defined XIOS

 
      subroutine set_surf_avg (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call set_surf_avg_tile (istr,iend,jstr,jend)
      return
      end
 
      subroutine set_surf_avg_tile (istr,iend,jstr,jend)
!
! Compute time-averaged fields within a tile.
! ------- ------------- ------ ------ - -----
! Because of syncronization issues, the delayed mode averaging
! procedure is used. This procedure implies that all fields to be
! averaged are sampled during the next time step, rather than at
! the end of the time step when they were computed.
!
! Thought this algorithm results in somewhat ackwad controlling
! logic it has the advantage that that all fields to be sampled
! correspond to exactly the same time, which is time step "n".
! Particularly, this is done this way because vertical velocity
! corresponding to the newly computed horizontal velocities
! becomes available only during the following time step.
! The same applies to the density field.
!
! The algorithm consists of three logical blocks: (1) initialization
! of the averages arrays: when mod(iic-1,navg).eq.1 the target arrays
! are set to the first contribution; (2) accumulation of averaged
! data, when mod(iic-1,navg).gt.1; and (3) adding the last
! contribution and scaling.
!
      implicit none
      integer istr,iend,jstr,jend, i,j, ilc, iflux
# ifdef SOLVE3D
     &                       , itrc, k
# endif
      real cff
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"

# include "averages.h"
# include "surf.h"

!
#include "compute_auxiliary_bounds.h"
!
      ilc=1+iic-ntstart  ! number of time step since restart
!
! calculation of averaged fluxes will only be performed if ntssurf_avg
! is a positive number
!


      if (ilc.gt.ntssurf_avg) then

        if (mod(ilc-ntssurf_avg,nwrtsurf_avg).eq.1) then

          if (ZEROTH_TILE) then
            timesurf_avg=time
!            MPI_master_only write(*,*) 'started averaging surf',iic,
!     &                                      ntssurf_avg,nwrtsurf_avg
          endif


             do j=JstrR,JendR
                do i=IstrR,IendR 
                  surft_avg(i,j)=t(i,j,N,nstp,1)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  surfs_avg(i,j)=t(i,j,N,nstp,2)
           enddo
          enddo           
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  surfz_avg(i,j)=zeta(i,j,fast_indx_out)
           enddo         
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  surfu_avg(i,j)=u(i,j,N,nstp)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  surfv_avg(i,j)=v(i,j,N,nstp)
           enddo
          enddo


        elseif (mod(ilc-ntssurf_avg,nwrtsurf_avg).gt.1) then

          if (ZEROTH_TILE) timesurf_avg=timesurf_avg+time

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  surft_avg(i,j) = 
     &                       surft_avg(i,j) +
     &                       t(i,j,N,nstp,1)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  surfs_avg(i,j) = 
     &                       surfs_avg(i,j) +
     &                       t(i,j,N,nstp,2)
           enddo
          enddo           
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  surfz_avg(i,j) = 
     &                       surfz_avg(i,j) +
     &                       zeta(i,j,fast_indx_out)
           enddo       
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  surfu_avg(i,j) =
     &                       surfu_avg(i,j) +
     &                       u(i,j,N,nstp)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  surfv_avg(i,j) =
     &                       surfv_avg(i,j) +
     &                       v(i,j,N,nstp)
           enddo
          enddo
        elseif (mod(ilc-ntssurf_avg,nwrtsurf_avg).eq.0) then

          cff=1./float(nwrtsurf_avg)
          if (ZEROTH_TILE) then
            timesurf_avg=cff*(timesurf_avg+time)
!            MPI_master_only write(*,*) 'finish aver. surf',iic,
!     &                                 ntssurf_avg,nwrtsurf_avg
          endif


             do j=JstrR,JendR
                do i=IstrR,IendR 
                  surft_avg(i,j) = cff *
     &                       ( surft_avg(i,j) +
     &                         t(i,j,N,nstp,1) )
             enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  surfs_avg(i,j) = cff *
     &                       ( surfs_avg(i,j) +
     &                         t(i,j,N,nstp,2) )
           enddo
            enddo         
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  surfz_avg(i,j) = cff *
     &                       ( surfz_avg(i,j) +
     &                         zeta(i,j,fast_indx_out) )
           enddo           
          enddo
              do j=JstrR,JendR
                do i=IstrR,IendR
                  surfu_avg(i,j) = cff *
     &                       ( surfu_avg(i,j) +
     &                         u(i,j,N,nstp) )
           enddo
          enddo
              do j=JstrR,JendR
                do i=IstrR,IendR
                  surfv_avg(i,j) = cff *
     &                       ( surfv_avg(i,j) +
     &                         v(i,j,N,nstp) )
           enddo
          enddo

        endif
      endif       !<-- iic.gt.ntsavg





      return
      end
#else /* DIAGNOSTICS_EK && AVERAGES */
      subroutine set_surf_avg_empty
      end
#endif /* DIAGNOSTICS_EK && AVERAGES */
