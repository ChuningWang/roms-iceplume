#include "cppdefs.h"
#if (defined DIAGNOSTICS_VRT && defined AVERAGES)
 
      subroutine set_diags_vrt_avg (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call set_diags_vrt_avg_tile (istr,iend,jstr,jend)
      return
      end
 
      subroutine set_diags_vrt_avg_tile (istr,iend,jstr,jend)
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
#include "diags_vrt.h"

!
#include "compute_auxiliary_bounds.h"
!
      ilc=1+iic-ntstart  ! number of time step since restart
!
! calculation of averaged fluxes will only be performed if ntsdiags_vrt_avg
! is a positive number
!


      if (ilc.gt.ntsdiags_vrt_avg) then
        if (mod(ilc-ntsdiags_vrt_avg,nwrtdiags_vrt_avg).eq.1) then

          if (ZEROTH_TILE) then
            timediags_vrt_avg=time
      MPI_master_only write(*,*) 'started averaging vrt',iic,
     &                    ntsdiags_vrt_avg,nwrtdiags_vrt_avg
          endif


             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtXadv_avg(i,j)=vrtXadv(i,j)  
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtYadv_avg(i,j)=vrtYadv(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtHdiff_avg(i,j)=vrtHdiff(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtCor_avg(i,j)=vrtCor(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtPrsgrd_avg(i,j)=vrtPrsgrd(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtHmix_avg(i,j)=vrtHmix(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtVmix_avg(i,j)=vrtVmix(i,j)
            enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtVmix2_avg(i,j)=vrtVmix2(i,j)
            enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtWind_avg(i,j)=vrtWind(i,j)
            enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtDrag_avg(i,j)=vrtDrag(i,j)
            enddo
          enddo


             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtrate_avg(i,j)=vrtrate(i,j)
           enddo
          enddo

# if defined DIAGNOSTICS_BARO
             do j=JstrR,JendR
                do i=IstrR,IendR
                  vrtBaro_avg(i,j)=vrtBaro(i,j)
           enddo
          enddo
# endif

# if defined M3FAST
             do j=JstrR,JendR
                do i=IstrR,IendR
                  vrtfast_avg(i,j)=vrtfast(i,j)
           enddo
          enddo
# endif

        elseif (mod(ilc-ntsdiags_vrt_avg,nwrtdiags_vrt_avg).gt.1) then

          if (ZEROTH_TILE) timediags_vrt_avg=timediags_vrt_avg+time

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtXadv_avg(i,j) = 
     &                       vrtXadv_avg(i,j) +
     &                       vrtXadv(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtYadv_avg(i,j) = 
     &                       vrtYadv_avg(i,j) +
     &                       vrtYadv(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtHdiff_avg(i,j) =
     &                       vrtHdiff_avg(i,j) +
     &                       vrtHdiff(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtCor_avg(i,j) =
     &                       vrtCor_avg(i,j) +
     &                       vrtCor(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtPrsgrd_avg(i,j) =
     &                       vrtPrsgrd_avg(i,j) +
     &                       vrtPrsgrd(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtHmix_avg(i,j) =
     &                       vrtHmix_avg(i,j) +
     &                       vrtHmix(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtVmix_avg(i,j) = 
     &                       vrtVmix_avg(i,j) +
     &                       vrtVmix(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtVmix2_avg(i,j) = 
     &                       vrtVmix2_avg(i,j) +
     &                       vrtVmix2(i,j)
           enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtWind_avg(i,j) = 
     &                       vrtWind_avg(i,j) +
     &                       vrtWind(i,j)
           enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtDrag_avg(i,j) = 
     &                       vrtDrag_avg(i,j) +
     &                       vrtDrag(i,j)
           enddo
          enddo


             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtrate_avg(i,j) = 
     &                    vrtrate_avg(i,j) + 
     &                    vrtrate(i,j)  
           enddo
          enddo

# if defined DIAGNOSTICS_BARO
             do j=JstrR,JendR
                do i=IstrR,IendR
                  vrtBaro_avg(i,j) =
     &                       vrtBaro_avg(i,j) +
     &                       vrtBaro(i,j)
           enddo
          enddo
# endif

# if defined M3FAST
             do j=JstrR,JendR
                do i=IstrR,IendR
                  vrtfast_avg(i,j) =
     &                       vrtfast_avg(i,j) +
     &                       vrtfast(i,j)
           enddo
          enddo
# endif

        elseif (mod(ilc-ntsdiags_vrt_avg,nwrtdiags_vrt_avg).eq.0) then

          cff=1./float(nwrtdiags_vrt_avg)
          if (ZEROTH_TILE) then
                        timediags_vrt_avg=cff*(timediags_vrt_avg+time)
      MPI_master_only write(*,*) 'finish aver. vrt',iic,
     &                       ntsdiags_vrt_avg,nwrtdiags_vrt_avg
          endif


             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtXadv_avg(i,j) = cff *
     &                       ( vrtXadv_avg(i,j) +
     &                         vrtXadv(i,j) )
             enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtYadv_avg(i,j) = cff *
     &                       ( vrtYadv_avg(i,j) +
     &                         vrtYadv(i,j) )
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtHdiff_avg(i,j) = cff *
     &                       ( vrtHdiff_avg(i,j) +
     &                         vrtHdiff(i,j) )
           enddo
          enddo
              do j=JstrR,JendR
                do i=IstrR,IendR
                  vrtCor_avg(i,j) = cff *
     &                       ( vrtCor_avg(i,j) +
     &                         vrtCor(i,j) )
           enddo
          enddo
              do j=JstrR,JendR
                do i=IstrR,IendR
                  vrtPrsgrd_avg(i,j) = cff *
     &                       ( vrtPrsgrd_avg(i,j) +
     &                         vrtPrsgrd(i,j) )
           enddo
          enddo
              do j=JstrR,JendR
                do i=IstrR,IendR
                  vrtHmix_avg(i,j) = cff *
     &                       ( vrtHmix_avg(i,j) +
     &                         vrtHmix(i,j) )
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtVmix_avg(i,j) = cff *
     &                       ( vrtVmix_avg(i,j) +
     &                         vrtVmix(i,j) )
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtVmix2_avg(i,j) = cff *
     &                       ( vrtVmix2_avg(i,j) +
     &                         vrtVmix2(i,j) )
           enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtWind_avg(i,j) = cff *
     &                       ( vrtWind_avg(i,j) +
     &                         vrtWind(i,j) )
           enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtDrag_avg(i,j) = cff *
     &                       ( vrtDrag_avg(i,j) +
     &                         vrtDrag(i,j) )
           enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  vrtrate_avg(i,j) = cff *
     &                    ( vrtrate_avg(i,j) + 
     &                      vrtrate(i,j) )
           enddo
          enddo

# if defined DIAGNOSTICS_BARO
              do j=JstrR,JendR
                do i=IstrR,IendR
                  vrtBaro_avg(i,j) = cff *
     &                       ( vrtBaro_avg(i,j) +
     &                         vrtBaro(i,j) )
           enddo
          enddo
# endif

# if defined M3FAST
              do j=JstrR,JendR
                do i=IstrR,IendR
                  vrtfast_avg(i,j) = cff *
     &                       ( vrtfast_avg(i,j) +
     &                         vrtfast(i,j) )
           enddo
          enddo
# endif

        endif
      endif       !<-- ilc.gt.ntsavg

      return
      end
#else /* DIAGNOSTICS_VRT && AVERAGES */
      subroutine set_diags_vrt_avg_empty
      end
#endif /* DIAGNOSTICS_VRT && AVERAGES */
