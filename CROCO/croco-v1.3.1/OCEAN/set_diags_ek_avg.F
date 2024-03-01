#include "cppdefs.h"
#if (defined DIAGNOSTICS_EK && defined AVERAGES)
 
      subroutine set_diags_ek_avg (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call set_diags_ek_avg_tile (istr,iend,jstr,jend)
      return
      end
 
      subroutine set_diags_ek_avg_tile (istr,iend,jstr,jend)
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
#include "diags_ek.h"

!
#include "compute_auxiliary_bounds.h"
!
      ilc=1+iic-ntstart  ! number of time step since restart
!
! calculation of averaged fluxes will only be performed if ntsdiags_ek_avg
! is a positive number
!

      
      if (ilc.gt.ntsdiags_ek_avg) then
        if (mod(ilc-ntsdiags_ek_avg,nwrtdiags_ek_avg).eq.1) then

          if (ZEROTH_TILE) then
            timediags_ek_avg=time
!            MPI_master_only write(*,*) 'started averaging ek',iic,
!     &                            ntsdiags_ek_avg,nwrtdiags_ek_avg
          endif


             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekHadv_avg(i,j)=ekHadv(i,j)  
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekHdiff_avg(i,j)=ekHdiff(i,j)
           enddo
          enddo           
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVadv_avg(i,j)=ekVadv(i,j)
           enddo         
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekCor_avg(i,j)=ekCor(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekPrsgrd_avg(i,j)=ekPrsgrd(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekHmix_avg(i,j)=ekHmix(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVmix_avg(i,j)=ekVmix(i,j)
            enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVmix2_avg(i,j)=ekVmix2(i,j)
            enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekWind_avg(i,j)=ekWind(i,j)
            enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekDrag_avg(i,j)=ekDrag(i,j)
            enddo
          enddo
# if defined DIAGNOSTICS_BARO
             do j=JstrR,JendR
                do i=IstrR,IendR
                  ekBaro_avg(i,j)=ekBaro(i,j)
           enddo
          enddo
# endif

# if defined M3FAST
             do j=JstrR,JendR
                do i=IstrR,IendR
                  ekfast_avg(i,j)=ekfast(i,j)
           enddo
          enddo
# endif

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekrate_avg(i,j)=ekrate(i,j)
           enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekvol_avg(i,j)=ekvol(i,j)
           enddo
          enddo

        elseif (mod(ilc-ntsdiags_ek_avg,nwrtdiags_ek_avg).gt.1) then

          if (ZEROTH_TILE) timediags_ek_avg=timediags_ek_avg+time

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekHadv_avg(i,j) = 
     &                       ekHadv_avg(i,j) +
     &                       ekHadv(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekHdiff_avg(i,j) = 
     &                       ekHdiff_avg(i,j) +
     &                       ekHdiff(i,j)
           enddo
          enddo           
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVadv_avg(i,j) = 
     &                       ekVadv_avg(i,j) +
     &                       ekVadv(i,j)
           enddo       
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekCor_avg(i,j) =
     &                       ekCor_avg(i,j) +
     &                       ekCor(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekPrsgrd_avg(i,j) =
     &                       ekPrsgrd_avg(i,j) +
     &                       ekPrsgrd(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekHmix_avg(i,j) =
     &                       ekHmix_avg(i,j) +
     &                       ekHmix(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVmix_avg(i,j) = 
     &                       ekVmix_avg(i,j) +
     &                       ekVmix(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVmix2_avg(i,j) = 
     &                       ekVmix2_avg(i,j) +
     &                       ekVmix2(i,j)
           enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekWind_avg(i,j) = 
     &                       ekWind_avg(i,j) +
     &                       ekWind(i,j)
           enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekDrag_avg(i,j) = 
     &                       ekDrag_avg(i,j) +
     &                       ekDrag(i,j)
           enddo
          enddo

# if defined DIAGNOSTICS_BARO
             do j=JstrR,JendR
                do i=IstrR,IendR
                  ekBaro_avg(i,j) =
     &                       ekBaro_avg(i,j) +
     &                       ekBaro(i,j)
           enddo
          enddo
# endif

# if defined M3FAST
             do j=JstrR,JendR
                do i=IstrR,IendR
                  ekfast_avg(i,j) =
     &                       ekfast_avg(i,j) +
     &                       ekfast(i,j)
           enddo
          enddo
# endif

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekrate_avg(i,j) = 
     &                    ekrate_avg(i,j) + 
     &                    ekrate(i,j)  
           enddo
          enddo


             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekvol_avg(i,j) =
     &                    ekvol_avg(i,j) +
     &                    ekvol(i,j)
           enddo
          enddo

        elseif (mod(ilc-ntsdiags_ek_avg,nwrtdiags_ek_avg).eq.0) then

          cff=1./float(nwrtdiags_ek_avg)
          if (ZEROTH_TILE) then
            timediags_ek_avg=cff*(timediags_ek_avg+time)
!            MPI_master_only write(*,*) 'finish aver. ek',iic,
!     &                       ntsdiags_ek_avg,nwrtdiags_ek_avg
          endif


             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekHadv_avg(i,j) = cff *
     &                       ( ekHadv_avg(i,j) +
     &                         ekHadv(i,j) )
             enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekHdiff_avg(i,j) = cff *
     &                       ( ekHdiff_avg(i,j) +
     &                         ekHdiff(i,j) )
           enddo
            enddo         
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVadv_avg(i,j) = cff *
     &                       ( ekVadv_avg(i,j) +
     &                         ekVadv(i,j) )
           enddo           
          enddo
              do j=JstrR,JendR
                do i=IstrR,IendR
                  ekCor_avg(i,j) = cff *
     &                       ( ekCor_avg(i,j) +
     &                         ekCor(i,j) )
           enddo
          enddo
              do j=JstrR,JendR
                do i=IstrR,IendR
                  ekPrsgrd_avg(i,j) = cff *
     &                       ( ekPrsgrd_avg(i,j) +
     &                         ekPrsgrd(i,j) )
           enddo
          enddo
              do j=JstrR,JendR
                do i=IstrR,IendR
                  ekHmix_avg(i,j) = cff *
     &                       ( ekHmix_avg(i,j) +
     &                         ekHmix(i,j) )
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVmix_avg(i,j) = cff *
     &                       ( ekVmix_avg(i,j) +
     &                         ekVmix(i,j) )
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVmix2_avg(i,j) = cff *
     &                       ( ekVmix2_avg(i,j) +
     &                         ekVmix2(i,j) )
           enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekWind_avg(i,j) = cff *
     &                       ( ekWind_avg(i,j) +
     &                         ekWind(i,j) )
           enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekDrag_avg(i,j) = cff *
     &                       ( ekDrag_avg(i,j) +
     &                         ekDrag(i,j) )
           enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekrate_avg(i,j) = cff *
     &                    ( ekrate_avg(i,j) + 
     &                      ekrate(i,j) )
           enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekvol_avg(i,j) = cff *
     &                    ( ekvol_avg(i,j) +
     &                      ekvol(i,j) )
           enddo
          enddo
# if defined DIAGNOSTICS_BARO
              do j=JstrR,JendR
                do i=IstrR,IendR
                  ekBaro_avg(i,j) = cff *
     &                       ( ekBaro_avg(i,j) +
     &                         ekBaro(i,j) )
           enddo
          enddo
# endif

# if defined M3FAST
              do j=JstrR,JendR
                do i=IstrR,IendR
                  ekfast_avg(i,j) = cff *
     &                       ( ekfast_avg(i,j) +
     &                         ekfast(i,j) )
           enddo
          enddo
# endif

        endif
      endif       !<-- iic.gt.ntsavg





# ifdef DIAGNOSTICS_EK_MLD

      if (ilc.gt.ntsdiags_ek_avg) then
        if (mod(ilc-ntsdiags_ek_avg,nwrtdiags_ek_avg).eq.1) then

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekHadv_mld_avg(i,j)=ekHadv_mld(i,j)  
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekHdiff_mld_avg(i,j)=ekHdiff_mld(i,j)
           enddo
          enddo           
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVadv_mld_avg(i,j)=ekVadv_mld(i,j)
           enddo         
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekCor_mld_avg(i,j)=ekCor_mld(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekPrsgrd_mld_avg(i,j)=ekPrsgrd_mld(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekHmix_mld_avg(i,j)=ekHmix_mld(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVmix_mld_avg(i,j)=ekVmix_mld(i,j)
            enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVmix2_mld_avg(i,j)=ekVmix2_mld(i,j)
            enddo
          enddo


             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekrate_mld_avg(i,j)=ekrate_mld(i,j)
           enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekvol_mld_avg(i,j)=ekvol_mld(i,j)
           enddo
          enddo
# if defined DIAGNOSTICS_BARO
             do j=JstrR,JendR
                do i=IstrR,IendR
                  ekBaro_mld_avg(i,j)=ekBaro_mld(i,j)
           enddo
          enddo
# endif
        elseif (mod(ilc-ntsdiags_ek_avg,nwrtdiags_ek_avg).gt.1) then

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekHadv_mld_avg(i,j) = 
     &                       ekHadv_mld_avg(i,j) +
     &                       ekHadv_mld(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekHdiff_mld_avg(i,j) = 
     &                       ekHdiff_mld_avg(i,j) +
     &                       ekHdiff_mld(i,j)
           enddo
          enddo           
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVadv_mld_avg(i,j) = 
     &                       ekVadv_mld_avg(i,j) +
     &                       ekVadv_mld(i,j)
           enddo       
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekCor_mld_avg(i,j) =
     &                       ekCor_mld_avg(i,j) +
     &                       ekCor_mld(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekPrsgrd_mld_avg(i,j) =
     &                       ekPrsgrd_mld_avg(i,j) +
     &                       ekPrsgrd_mld(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekHmix_mld_avg(i,j) =
     &                       ekHmix_mld_avg(i,j) +
     &                       ekHmix_mld(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVmix_mld_avg(i,j) = 
     &                       ekVmix_mld_avg(i,j) +
     &                       ekVmix_mld(i,j)
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVmix2_mld_avg(i,j) = 
     &                       ekVmix2_mld_avg(i,j) +
     &                       ekVmix2_mld(i,j)
           enddo
          enddo


             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekrate_mld_avg(i,j) = 
     &                    ekrate_mld_avg(i,j) + 
     &                    ekrate_mld(i,j)  
           enddo
          enddo


             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekvol_mld_avg(i,j) =
     &                    ekvol_mld_avg(i,j) +
     &                    ekvol_mld(i,j)
           enddo
          enddo
# if defined DIAGNOSTICS_BARO
             do j=JstrR,JendR
                do i=IstrR,IendR
                  ekBaro_mld_avg(i,j) =
     &                       ekBaro_mld_avg(i,j) +
     &                       ekBaro_mld(i,j)
           enddo
          enddo
# endif
        elseif (mod(ilc-ntsdiags_ek_avg,nwrtdiags_ek_avg).eq.0) then

          cff=1./float(nwrtdiags_ek_avg)

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekHadv_mld_avg(i,j) = cff *
     &                       ( ekHadv_mld_avg(i,j) +
     &                         ekHadv_mld(i,j) )
             enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekHdiff_mld_avg(i,j) = cff *
     &                       ( ekHdiff_mld_avg(i,j) +
     &                         ekHdiff_mld(i,j) )
           enddo
            enddo         
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVadv_mld_avg(i,j) = cff *
     &                       ( ekVadv_mld_avg(i,j) +
     &                         ekVadv_mld(i,j) )
           enddo           
          enddo
              do j=JstrR,JendR
                do i=IstrR,IendR
                  ekCor_mld_avg(i,j) = cff *
     &                       ( ekCor_mld_avg(i,j) +
     &                         ekCor_mld(i,j) )
           enddo
          enddo
              do j=JstrR,JendR
                do i=IstrR,IendR
                  ekPrsgrd_mld_avg(i,j) = cff *
     &                       ( ekPrsgrd_mld_avg(i,j) +
     &                         ekPrsgrd_mld(i,j) )
           enddo
          enddo
              do j=JstrR,JendR
                do i=IstrR,IendR
                  ekHmix_mld_avg(i,j) = cff *
     &                       ( ekHmix_mld_avg(i,j) +
     &                         ekHmix_mld(i,j) )
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVmix_mld_avg(i,j) = cff *
     &                       ( ekVmix_mld_avg(i,j) +
     &                         ekVmix_mld(i,j) )
           enddo
          enddo
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekVmix2_mld_avg(i,j) = cff *
     &                       ( ekVmix2_mld_avg(i,j) +
     &                         ekVmix2_mld(i,j) )
           enddo
          enddo


             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekrate_mld_avg(i,j) = cff *
     &                    ( ekrate_mld_avg(i,j) + 
     &                      ekrate_mld(i,j) )
           enddo
          enddo

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  ekvol_mld_avg(i,j) = cff *
     &                    ( ekvol_mld_avg(i,j) +
     &                      ekvol_mld(i,j) )
           enddo
          enddo
# if defined DIAGNOSTICS_BARO
              do j=JstrR,JendR
                do i=IstrR,IendR
                  ekBaro_mld_avg(i,j) = cff *
     &                       ( ekBaro_mld_avg(i,j) +
     &                         ekBaro_mld(i,j) )
           enddo
          enddo
# endif
        endif
      endif       !<-- ilc.gt.ntsavg






# endif



      return
      end
#else /* DIAGNOSTICS_EK && AVERAGES */
      subroutine set_diags_ek_avg_empty
      end
#endif /* DIAGNOSTICS_EK && AVERAGES */
