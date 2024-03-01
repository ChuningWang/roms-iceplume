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
#if defined DIAGNOSTICS_EDDY && defined AVERAGES  && ! defined XIOS

      subroutine set_diags_eddy_avg(tile)

!      USE param

      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call set_diags_eddy_avg_tile(Istr,Iend,Jstr,Jend)
      return
      end

      subroutine set_diags_eddy_avg_tile(Istr,Iend,Jstr,Jend)
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
      integer trd, omp_get_thread_num
      real    GRho, cff, cff1, cffu, cffv, cffb, cffug, cffvg
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "forces.h"
# include "averages.h"
# include "diags_eddy.h"
# include "coupling.h"
!
# include "work.h"
# include "private_scratch.h"
!
# include "compute_auxiliary_bounds.h"
!
      ilc=1+iic-ntstart   ! number of time step since restart
!
! calculation of averaged fluxes will only be performed if ntsdiags_eddy_avg
! is a positive number
!
      if (ilc.gt.ntsdiags_eddy_avg) then


!###############################################

       ! Compute true vertical velocity (m/s).
        trd=omp_get_thread_num()
        call Wvlcty_tile (Istr,Iend,Jstr,Jend, workr, A2d(1,1,trd),
     &                                  A2d(1,1,trd), A2d(1,2,trd))

!###############################################

        if (mod(ilc-ntsdiags_eddy_avg,nwrtdiags_eddy_avg).eq.1) then
          cff =1.0
          cff1=0.0
          if (ZEROTH_TILE) then
            timediags_eddy_avg=time
!            MPI_master_only write(*,*) 'started averaging eddy',iic,
!     &                  ntsdiags_eddy_avg,nwrtdiags_eddy_avg
          endif
        elseif (mod(ilc-ntsdiags_eddy_avg,nwrtdiags_eddy_avg).gt.1) then
          cff =1.0
          cff1=1.0
          if (ZEROTH_TILE) then
           timediags_eddy_avg=timediags_eddy_avg+time
          endif
        elseif (mod(ilc-ntsdiags_eddy_avg,nwrtdiags_eddy_avg).eq.0) then
          cff =1.0/float(nwrtdiags_eddy_avg)
          cff1=1.0
          if (ZEROTH_TILE) then
            timediags_eddy_avg=cff*(timediags_eddy_avg+time)
!            MPI_master_only write(*,*) 'finished averaging eddy',
!     &                  iic,ntsdiags_eddy_avg,nwrtdiags_eddy_avg
          endif
        endif

!###############################################
          
        do j=JstrR,JendR
          do i=IstrR,IendR
            eddyzz_avg(i,j)= cff *
     &                       ( cff1*eddyzz_avg(i,j)
     &                                +Zt_avg1(i,j)**2 )
          enddo
        enddo

        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR

              GRho=-1*g/rho0

              cffu = 0.5*(u(i+1,j,k,nstp)+u(i,j,k,nstp))
              cffv = 0.5*(v(i,j+1,k,nstp)+v(i,j,k,nstp))

#  ifdef SPLIT_EOS
              cffb =rho1(i,j,k)+qp1(i,j,k)
     &                          *(z_w(i,j,N)-z_r(i,j,k))
#  else
              cffb =rho(i,j,k)
#  endif

              if (k.eq.1) then
                eddyubu_avg(i,j)= cff *
     &                           ( cff1*eddyubu_avg(i,j)
     &                             +cffu * bustr(i,j) )
                eddyvbv_avg(i,j)= cff *
     &                           ( cff1*eddyvbv_avg(i,j)
     &                             +cffv * bvstr(i,j) )
              endif

              if (k.eq.N) then
                  cffug = -g/f(i,j)*(
# ifdef MASKING
     &                    vmask(i,j+1)*
# endif
     &               (zeta(i,j+1,fast_indx_out)-zeta(i,j,fast_indx_out))
     &                   *0.5*(pm(i,j+1)+pm(i,j))
     &                    )
                  cffvg = g/f(i,j)*(
# ifdef MASKING
     &                    umask(i+1,j)*
# endif
     &               (zeta(i+1,j,fast_indx_out)-zeta(i,j,fast_indx_out))
     &                   *0.5*(pn(i+1,j)+pn(i,j))
     &                    )
                eddyusu_avg(i,j)= cff *
     &                           ( cff1*eddyusu_avg(i,j)
     &                             +cffu * sustr(i,j) )
                eddyvsv_avg(i,j)= cff *
     &                           ( cff1*eddyvsv_avg(i,j)
     &                             +cffv * svstr(i,j) )
                eddyugsu_avg(i,j)= cff *
     &                           ( cff1*eddyugsu_avg(i,j)
     &                             +cffug * sustr(i,j) )
                eddyvgsv_avg(i,j)= cff *
     &                           ( cff1*eddyvgsv_avg(i,j)
     &                             +cffvg * svstr(i,j) )
              endif

              eddyuu_avg(i,j,k)= cff *
     &                       ( cff1*eddyuu_avg(i,j,k)
     &                         +cffu**2 )
              eddyvv_avg(i,j,k)= cff *
     &                       ( cff1*eddyvv_avg(i,j,k)
     &                         +cffv**2 )
              eddyuv_avg(i,j,k)= cff *
     &                       ( cff1*eddyuv_avg(i,j,k)
     &                         +cffu*cffv )
              eddyub_avg(i,j,k)= cff *
     &                       ( cff1*eddyub_avg(i,j,k)
     &                     +GRho*cffu*cffb )
              eddyvb_avg(i,j,k)= cff *
     &                       ( cff1*eddyvb_avg(i,j,k)
     &                     +GRho*cffv*cffb )
              eddywb_avg(i,j,k)= cff *
     &                       ( cff1*eddywb_avg(i,j,k)
     &                     +GRho*workr(i,j,k)*cffb )
              eddyuw_avg(i,j,k)= cff *
     &                       ( cff1*eddyuw_avg(i,j,k)
     &                          +workr(i,j,k)*cffu )
              eddyvw_avg(i,j,k)= cff *
     &                       ( cff1*eddyvw_avg(i,j,k)
     &                          +workr(i,j,k)*cffv )

            enddo
          enddo
        enddo

      endif

      return
      end
#else /* DIAGNOSTICS_EDDY && AVERAGES */
      subroutine set_diags_eddy_avg_empty
      end
#endif /* DIAGNOSTICS_EDDY && AVERAGES */
