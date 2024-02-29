! $Id: diag.F 1458 2014-02-03 15:01:25Z gcambon $
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

      subroutine diag (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
#include "param.h"
#include "private_scratch.h"
#include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call diag_tile (Istr,Iend,Jstr,Jend, A2d(1,1,trd),A2d(1,2,trd))
      return
      end

      subroutine diag_tile (Istr,Iend,Jstr,Jend, ke2d,pe2d)
      implicit none     
#include "param.h" 
#include "grid.h"
#ifdef SOLVE3D
# include "ocean3d.h"
#else
# include "ocean2d.h"
#endif
#include "scalars.h"
#include "mpi_cpl.h"
#if defined NBQ && defined DIAG_CFL
# include "nbq.h"
#endif 
      integer Istr,Iend,Jstr,Jend,    i,j,k, NSUB,
     &                                trd, omp_get_thread_num
      real ke2d(PRIVATE_2D_SCRATCH_ARRAY),
     &     pe2d(PRIVATE_2D_SCRATCH_ARRAY)
      real*QUAD cff, my_avgke, my_avgpe, my_volume
      character echar*8
#ifdef MPI
      include 'mpif.h'
      integer size, step, status(MPI_STATUS_SIZE), ierr
      real*QUAD buff(9)
      common /xyz/ buff
#endif

#ifdef DIAG_CFL
      real my_Cu_Adv, my_Cu_W, ciV, cx, cw
# ifdef NBQ
      real my_Cu_Nbq_X, my_Cu_Nbq_Y, my_Cu_Nbq_Z
      real cx_nbq,cy_nbq,cz_nbq, alphax,alphay
      real my_Cu_Nbq
      real Cu_Nbq
# endif
      integer my_i_cmax, my_j_cmax, my_k_cmax
#endif

#ifndef NBQ_MASS
# define Hzr_half_nbq Hz
#else
# define Hzr_half_nbq Hzr
#endif
!
#ifdef DIAG_CFL
      my_Cu_Adv=0.  ;  my_Cu_W=0.
# ifdef NBQ
      my_Cu_Nbq_X=0.; my_Cu_Nbq_Y=0.; my_Cu_Nbq_Z=0.
      my_Cu_Nbq =0.
# endif
      if (mod(iic-1,ninfo).eq.0) then

        do j=jstr,jend
          do k=N,1,-1
            do i=istr,iend
# ifdef MASKING
              ciV=dt*rmask(i,j)*pm(i,j)*pn(i,j)/Hz(i,j,k)
# else
              ciV=dt*pm(i,j)*pn(i,j)/Hz(i,j,k)
# endif
              cw=ciV*( max( We(i,j,k  )
# ifdef VADV_ADAPT_IMP     
     &                     +Wi(i,j,k  )
# endif     
     &                          , 0.)
     &                 -min(We(i,j,k-1)
# ifdef VADV_ADAPT_IMP 
     &                     +Wi(i,j,k-1)
# endif        
     &                            , 0.) )

              cx=cw+ciV*( max(HUon(i+1,j,k), 0.)-min(HUon(i,j,k), 0.)
     &                   +max(HVom(i,j+1,k), 0.)-min(HVom(i,j,k), 0.))

              if (cx .gt. my_Cu_Adv) then
                my_Cu_Adv=cx ; my_Cu_W=cw
              endif

# ifdef NBQ
              cz_nbq=(soundspeed_nbq(i,j)*dtnbq/max(1.D-30,Hzr(i,j,k)))
#  ifdef MASKING
     &                *rmask(i,j) 
#  endif    
              my_Cu_Nbq_Z=max(my_Cu_Nbq_Z,cz_nbq)
            
              cx_nbq=soundspeed_nbq(i,j)*dtnbq*pm_u(i,j)
#  ifdef MASKING
     &                *umask(i,j) 
#  endif
              my_Cu_Nbq_X=max(my_Cu_Nbq_X,cx_nbq)
          
              cy_nbq=soundspeed_nbq(i,j)*dtnbq*pn_v(i,j)
#  ifdef MASKING
     &                *vmask(i,j) 
#  endif

              my_Cu_Nbq_Y=max(my_Cu_Nbq_Y,cy_nbq)
          
              alphax=(1/max(1.D-30,Hzr(i,j,k)))*
     &               (0.5*(z_r(i+1,j,k)-z_r(i-1,j,k))*pm(i,j))
#  ifdef MASKING
     &                *rmask(i,j) 
#  endif 

              alphay=(1/max(1.D-30,Hzr(i,j,k)))*
     &               (0.5*(z_r(i,j+1,k)-z_r(i,j-1,k))*pn(i,j))
#  ifdef MASKING
     &                *rmask(i,j) 
#  endif 
        
              my_Cu_Nbq=max(my_Cu_Nbq,
     &                (alphax*cz_nbq-cx_nbq)**2+(alphay*cz_nbq-cy_nbq)**2)
              if (Hzr_half_nbq(i,j,k).lt.1.D-10) my_Cu_Nbq_Z=0.D0
# endif /* NBQ */
            enddo
          enddo
        enddo            
      endif
# endif /* DIAG_CFL */

!
! Compute and report volume averaged kinetic, potential and total
!-----------------------------------------------------------------
! energy densities for either two- (shallow water) or three-
! dimensional versions of the model.
!
! At first, compute kinetic and potential energies, as well as total
! volume within the tile [subdomain of indices (Istr:Iend,Jstr:Jend)]
! by individual threads. In the case of three dimensions also perform
! verical summation at this stage.
!
      if (mod(iic-1,ninfo).eq.0) then
        do j=Jstr,Jend
#ifdef SOLVE3D
          do i=Istr,Iend
            ke2d(i,j)=0.
            pe2d(i,j)=0.5*g*z_w(i,j,N)*z_w(i,j,N)
                                     !!!zeta(i,j,kstp)*zeta(i,j,kstp)
          enddo
          cff=g/rho0
          do k=N,1,-1
            do i=Istr,Iend
             ke2d(i,j)=ke2d(i,j)+HZR(i,j,k)*0.25*(
     &                               u(i  ,j,k,nstp)*u(i,j,k,nstp)+
     &                               u(i+1,j,k,nstp)*u(i+1,j,k,nstp)+
     &                               v(i,j  ,k,nstp)*v(i,j  ,k,nstp)+
     &                               v(i,j+1,k,nstp)*v(i,j+1,k,nstp))
             pe2d(i,j)=pe2d(i,j)+cff*HZR(i,j,k)*rho(i,j,k)
     &                                       *(z_r(i,j,k)-z_w(i,j,0))
            enddo
          enddo
#else
          cff=0.5*g
          do i=Istr,Iend
            ke2d(i,j)=(zeta(i,j,krhs)+h(i,j))*0.25*(
     &                             ubar(i  ,j,krhs)*ubar(i  ,j,krhs)+
     &                             ubar(i+1,j,krhs)*ubar(i+1,j,krhs)+
     &                             vbar(i,j  ,krhs)*vbar(i,j  ,krhs)+
     &                             vbar(i,j+1,krhs)*vbar(i,j+1,krhs))
            pe2d(i,j)=cff*zeta(i,j,krhs)*zeta(i,j,krhs)
          enddo
#endif /* SOLVE3D */
        enddo
!
! After that integrate horizontally within the subdomain tile. Here,
! in order to reduce the round-off errors, the summation is performed
! in two stages, first the index j is collapsed, then in index i.
! In this order the partial sums consist on much fewer number of
! terms than if it would be the case of a straightforward two-
! dimensional summation. Thus adding numbers which are orders of
! magnitude apart is avoided. Also note that the partial sums are
! stored as quadro precision numbers for the same purpose.
!
        do i=Istr,Iend
          pe2d(i,Jend+1)=0.D0
          pe2d(i,Jstr-1)=0.D0
          ke2d(i,Jstr-1)=0.D0
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            cff=1./(pm(i,j)*pn(i,j))
#ifdef SOLVE3D
            pe2d(i,Jend+1)=pe2d(i,Jend+1)+cff*(z_w(i,j,N)-z_w(i,j,0))
#else
            pe2d(i,Jend+1)=pe2d(i,Jend+1)+cff*(zeta(i,j,krhs)+h(i,j))
#endif
            pe2d(i,Jstr-1)=pe2d(i,Jstr-1)+cff*pe2d(i,j)
            ke2d(i,Jstr-1)=ke2d(i,Jstr-1)+cff*ke2d(i,j)

          enddo
        enddo

        my_volume=0.
        my_avgpe=0.
        my_avgke=0.
        do i=Istr,Iend
          my_volume=my_volume+pe2d(i,Jend+1)
          my_avgpe =my_avgpe +pe2d(i,Jstr-1)
          my_avgke =my_avgke +ke2d(i,Jstr-1)
        enddo
        if (SINGLE_TILE_MODE) then
          NSUB=1
        else
          NSUB=NSUB_X*NSUB_E
        endif
!
! Perform global summation: whoever gets first to the critical region
! resets global sums before global summation starts; after the global
! summation is completed, thread, which is the last one to enter the
! critical region, finalizes the computation of diagnostics and
! prints them out. 
!
C$OMP CRITICAL (diag_cr_rgn)
          if (tile_count.eq.0) then
            volume=QuadZero               ! <-- Reset global sums for
            avgke= QuadZero               ! <-- multithreaded (shared
            avgpe= QuadZero               ! <-- memory) summation.
          endif 
          volume=volume+my_volume         ! Perform global 
          avgke =avgke +my_avgke          ! summation among 
          avgpe =avgpe +my_avgpe          ! the threads

#ifdef DIAG_CFL
          if (tile_count.eq.0) then
             Cu_Adv3d=my_Cu_Adv
             Cu_W=my_Cu_W
# ifdef NBQ
             Cu_Nbq_X=my_Cu_Nbq_X
             Cu_Nbq_Y=my_Cu_Nbq_Y
             Cu_Nbq_Z=my_Cu_Nbq_Z
             Cu_Nbq  = my_Cu_Nbq
# endif             
          else
            if (my_Cu_Adv.gt.Cu_Adv3d) then
               Cu_Adv3d=my_Cu_Adv
               Cu_W=my_Cu_W
            endif
# ifdef NBQ
            if(my_Cu_Nbq_X.gt.Cu_Nbq_X)  Cu_Nbq_X = my_Cu_Nbq_X
            if(my_Cu_Nbq_Y.gt.Cu_Nbq_Y)  Cu_Nbq_X = my_Cu_Nbq_Y
            if(my_Cu_Nbq_Z.gt.Cu_Nbq_Z)  Cu_Nbq_X = my_Cu_Nbq_Z    
            if(my_Cu_Nbq.gt.Cu_Nbq)      Cu_Nbq   = my_Cu_Nbq
# endif             
          endif
#endif /* DIAG_CFL */
          tile_count=tile_count+1         ! This counter identifies
          if (tile_count.eq.NSUB) then    ! the last thread, whoever
            tile_count=0                  ! it is, not always master.
#ifdef MPI
            if (NNODES.gt.1) then         ! Perform global summation 
              size=NNODES                 ! among MPI processes
   1           step=(size+1)/2 
                if (mynode.ge.step .and. mynode.lt.size) then
                  buff(1)=volume
                  buff(2)=avgke           ! This is MPI_Reduce
                  buff(3)=avgpe
# ifdef SOLVE3D
#  ifdef DIAG_CFL
                  buff(4)=Cu_Adv3d
                  buff(5)=Cu_W
#   ifdef NBQ                  
                  buff(6)=Cu_Nbq_X
                  buff(7)=Cu_Nbq_Y
                  buff(8)=Cu_Nbq_Z
                  buff(9)=Cu_Nbq
#   endif                  
#  endif
# endif
                  call MPI_Send (buff,  9, MPI_DOUBLE_PRECISION,
     &                 mynode-step, 17, MPI_COMM_WORLD,      ierr)
                elseif (mynode .lt. size-step) then
                  call MPI_Recv (buff,  9, MPI_DOUBLE_PRECISION,
     &                 mynode+step, 17, MPI_COMM_WORLD, status, ierr)
                  volume=volume+buff(1)
                  avgke=avgke+  buff(2)
                  avgpe=avgpe+  buff(3)
# ifdef SOLVE3D
#  ifdef DIAG_CFL
                  if (buff(4).gt.Cu_Adv3d) then
                     Cu_Adv3d=buff(4)
                     Cu_W=buff(5)
                  endif
#   ifdef NBQ  
                  if(buff(6).gt.Cu_Nbq_X)  Cu_Nbq_X = buff(6)             
                  if(buff(7).gt.Cu_Nbq_Y)  Cu_Nbq_Y = buff(7)   
                  if(buff(8).gt.Cu_Nbq_Z)  Cu_Nbq_Z = buff(8)
                  if(buff(9).gt.Cu_Nbq  )  Cu_Nbq   = buff(9)                  
#   endif                
#  endif
# endif
                endif
               size=step
              if (size.gt.1) goto 1
            endif
            if (mynode.eq.0) then
#endif /* MPI */

              avgke=avgke/volume          ! Compute and print global
              avgpe=avgpe/volume          ! diagnostics (last thread
              avgkp=avgke+avgpe           ! of master MPI process)

              if (first_time.eq.0) then
                first_time=1
                write(stdout,2) 'STEP','time[DAYS]','KINETIC_ENRG',
     &                  'POTEN_ENRG','TOTAL_ENRG','NET_VOLUME','trd'
   2            format(1x,A4,3x,A10,1x,A12,4x,A10,4x,A10,4x,A10,3x,A3)
              endif
              trd=omp_get_thread_num()
              write(stdout,3)iic-1,tdays,avgke,avgpe,avgkp,volume,trd
#ifdef CALENDAR 
              write(stdout,'(a)') tool_sectodat(time)
#endif
   3          format(I8, F12.5, 1PE16.9, 3(1PE14.7), I3)
#ifdef DIAG_CFL
              if (mod(iic-1,ninfo).eq.0) then
              write(stdout,4) 'CFL --  INT_3DADV  = ',Cu_adv3d,
     &                             '   EXT_GWAVES = ',Cu_max
              endif
   4          format(10x,A,F6.2,A,F6.2)
# ifdef NBQ
              if (mod(iic-1,ninfo).eq.0) then
              write(stdout,5) '        NBQ_HADV   = ',
     &                                 sqrt(Cu_Nbq_X**2+Cu_Nbq_Y**2),
     &                             '   NBQ_IMP    = ',sqrt(Cu_Nbq),
     &                             '   NBQ_VSOUND = ',Cu_Nbq_Z
 
              endif
   5          format(10x,A,F6.2,A,F6.2,A,F8.2)
# endif   
#endif   
!
! Raise may_day_flag to stop computations in the case of blowing up.
! [Criterion for blowing up here is the numerical overflow, so that
! avgkp is 'INF' or 'NAN' (any mix of lover and uppercase letters),
! therefore it is sufficient to check for the presence of letter 'N'.
!
              write(echar,'(1PE8.1)') avgkp
              do i=1,8
               if (echar(i:i).eq.'N' .or. echar(i:i).eq.'n'
     &                      .or. echar(i:i).eq.'*') may_day_flag=1 
              enddo
#ifdef MPI
            endif    ! <-- mynode.eq.0
#endif
/*
            if (ninfo.eq.1 .and. iic.gt.   8) ninfo=2
            if (ninfo.eq.2  .and. iic.gt. 16) ninfo=4
            if (ninfo.eq.4  .and. iic.gt. 32) ninfo=8
            if (ninfo.eq.8  .and. iic.gt. 64) ninfo=16
            if (ninfo.eq.16 .and. iic.gt.128) ninfo=32
            if (ninfo.eq.32 .and. iic.gt.256) ninfo=64
*/
          endif

#if defined MPI && defined DEBUG && !defined AGRIF
      call flush
#endif
  
C$OMP END CRITICAL (diag_cr_rgn)
      endif
      return
      end

