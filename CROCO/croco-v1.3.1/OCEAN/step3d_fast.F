! $Id: step3d_fast.F 1615 2014-12-17 13:27:07Z rblod $
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

#ifdef M3FAST

      subroutine step3d_fast (tile)
      implicit none
      integer tile, trd
# include "param.h"
# include "private_scratch.h"
C$    integer omp_get_thread_num
# include "compute_tile_bounds.h"
      trd=0
C$    trd=omp_get_thread_num()
!      
!======================================================================
!                       step3d_fast
!======================================================================
!
      call step3d_fast_tile ( Istr,Iend,Jstr,Jend
     &                   , A2d(1, 1,trd), A2d(1, 2,trd), A2d(1, 3,trd)
     &                   , A2d(1, 4,trd), A2d(1, 5,trd), A2d(1, 6,trd)
     &                   , A2d(1, 7,trd), A2d(1, 8,trd)
     &                   , A2d(1, 9,trd), A2d(1,10,trd)
     &                   , A2d(1,11,trd), A2d(1,12,trd)
     &                   , A2d(1,13,trd), A2d(1,14,trd)
     &                   , A2d(1,15,trd), A2d(1,16,trd), A2d(1,17,trd)
# ifdef NBQ
     &                   , A3d(1, 1,trd), A3d(1, 2,trd)
     &                   , A3d(1, 3,trd), A3d(1, 4,trd)
     &                   , A3d(1, 5,trd), A3d(1, 6,trd)
# endif
     &                   )

      return
      end
! 
!======================================================================
!                      step3d_fast_tile
!======================================================================
!
      subroutine step3d_fast_tile (Istr,Iend,Jstr,Jend
     &                          ,Dnew,rubar,rvbar
     &                          ,Drhs, UFx,UFe
     &                          ,VFx,VFe
     &                          ,urhs,vrhs 
     &                          ,DUon,DVom
     &                          ,ru_ext_nbq_sum, rv_ext_nbq_sum
     &                          ,ru_ext_nbq_old, rv_ext_nbq_old,work
# ifdef NBQ
     &                          ,Hzw_half_nbq_inv,  Hzr_half_nbq_inv 
     &                          ,Hzw_half_nbq_inv_u,Hzw_half_nbq_inv_v
     &                          ,Hzu_half_qdmu,     Hzv_half_qdmv
# endif
     &                          )
!
!
!***********************************************************************
!
!                    SOLVE FAST MODE 3D EQUATIONS
!
!***********************************************************************
!
!  This routines:
!  1- Computes non-NBQ RHS forcing terms of momentum equations. First
!     computes the barotropic (external) RHS forcing term (rubar,rvbar)
!     then adds it to the internal RHS forcing (computed in pre_step3d). 
!  2- Solves the 3D momentum conservation equations for fast-mode 
!     components (qdmu_nbq, qdmv_nbq, qdmw_nbq) by time integration of 
!     all forces:
!       Compressible pressure force + second viscosity + gravity
!       + NT Coriolis force + restoring force + non-NBQ RHS forces
!  3- Solves mass conservation equation, i.e., computes compressible 
!     density rho_nbq by time integration of momentum divergence
!
!  In this version, a first guest of zeta is derived from the surface 
!  vertical velocity (surface characteristic relation) instead of the 
!  depth-averaged conservation of mass. This satisfies dynamical coupling 
!  with the surface layer. After solving the 3D momentum equations, a 
!  final zeta field is diagnozed from mass conservation (then Hz is also
!  corrected for the internal time step).
!
!  W-momentum equation is solved with explicit or implicit methods:
!  - Explicit scheme: w-momentum is updated right after (and the same 
!                     way as) u- and v-momentum.
!  - Implicit scheme: horizontal component of divergence is first 
!                     precomputed (as required by fast-mode mass 
!                     conservation) before tridiagonal Gauss Elimination 
!                     is carried out for qdmw_nbq(m).
!
!  For all components, a Forward-Backward scheme is implemented:
!  - Explicit scheme: Forward:  zeta, qdmu_nbq, qdmw_nbq. 
!                     Backward:       rho_nbq.
!  - Implicit scheme: Forward:  zeta, qdmu_nbq. 
!                     Backward:       qdmw_nbq, rho_nbq.
!
!  In the NBQ_PERF option, the vertical grid is not evolving at fast
!  time step to gain computational time.
!
!***********************************************************************
!
      implicit none
# include "param.h"
      integer Istr,Iend,Jstr,Jend, i,j,k, kbak,kold,
     &        imin,imax,jmin,jmax,
     &        k1, k2, kp1
# ifdef MPI
     &       ,err
# endif
# ifdef PSOURCE
     &       ,is
# endif
      real    mybeta,myalpha,myepsilon,mygamma,
     &        VMAX,VMAXL, cff,cff0,cff1,cff2,cff3, 
     &        DUnew,DVnew, dum_s

# ifdef NBQ_THETAIMP
      real,parameter :: thetaimp_nbq = 0.5 
# else
      real,parameter :: thetaimp_nbq = 1.
# endif
!
! Parameters for NBQ pressure gradient scheme
!
      real,parameter :: gammau=1.0  ! 1.0  --> second-order
                                    ! 0.45 --> increased stability range
      real,parameter :: gammau_2=(1./3.)*(1.-gammau)

      real
     &                Dnew(PRIVATE_2D_SCRATCH_ARRAY),
     &               rubar(PRIVATE_2D_SCRATCH_ARRAY),
     &               rvbar(PRIVATE_2D_SCRATCH_ARRAY),
     &                Drhs(PRIVATE_2D_SCRATCH_ARRAY),
     &                 UFx(PRIVATE_2D_SCRATCH_ARRAY),
     &                 UFe(PRIVATE_2D_SCRATCH_ARRAY),
     &                 VFx(PRIVATE_2D_SCRATCH_ARRAY),
     &                 VFe(PRIVATE_2D_SCRATCH_ARRAY),
     &                urhs(PRIVATE_2D_SCRATCH_ARRAY),  
     &                vrhs(PRIVATE_2D_SCRATCH_ARRAY),
     &                DUon(PRIVATE_2D_SCRATCH_ARRAY),
     &                DVom(PRIVATE_2D_SCRATCH_ARRAY)
      real
     &       ru_ext_nbq_sum(PRIVATE_2D_SCRATCH_ARRAY),
     &       rv_ext_nbq_sum(PRIVATE_2D_SCRATCH_ARRAY),
     &       ru_ext_nbq_old(PRIVATE_2D_SCRATCH_ARRAY),
     &       rv_ext_nbq_old(PRIVATE_2D_SCRATCH_ARRAY),
     &                 work(PRIVATE_2D_SCRATCH_ARRAY)
# ifdef NBQ
      real    Hzu_half_qdmu(PRIVATE_2D_SCRATCH_ARRAY,0:N), 
     &        Hzv_half_qdmv(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &     Hzw_half_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &     Hzr_half_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,  N),
     &   Hzw_half_nbq_inv_u(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &   Hzw_half_nbq_inv_v(PRIVATE_2D_SCRATCH_ARRAY,0:N)
      real 
     &    dthetadiv_nbqdz_u(PRIVATE_2D_SCRATCH_ARRAY,2),
     &    dthetadiv_nbqdz_v(PRIVATE_2D_SCRATCH_ARRAY,2),
     &    dthetadiv_nbqdz_w(PRIVATE_2D_SCRATCH_ARRAY,2),
     &              dZdxq_u(PRIVATE_2D_SCRATCH_ARRAY,2),
     &              dZdyq_v(PRIVATE_2D_SCRATCH_ARRAY,2),
     &                   FX(PRIVATE_2D_SCRATCH_ARRAY),
     &                   FY(PRIVATE_2D_SCRATCH_ARRAY),
     &                   FC(PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &                   DC(PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &                   CF(PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &         qdmw_nbq_old(PRIVATE_2D_SCRATCH_ARRAY,0:N)
# endif
# ifdef UV_COR_NT
      real
     &               ntcoru(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &               ntcorv(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &               ntcorw(PRIVATE_2D_SCRATCH_ARRAY,0:N)
# endif

# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "coupling.h"
# include "forces.h"
# ifdef MRL_WCI
      real vstu,ustv,dudx,dvdx,dude,dvde
# endif
# include "mixing.h"
# include "climat.h"
# include "scalars.h"
# include "sources.h"
# ifdef AGRIF
#  include "zoom.h"
        integer irhox, irhoy, irhot
# endif
# if defined INTERNAL || defined BODYTIDE
      real,parameter :: U0 = 0.02
      real,parameter :: omega = 2.*pi/(12.4*3600.)
# endif
# ifdef ACOUSTIC
      real dist_d  
# endif
# ifdef WET_DRY
      real cff1_WD,cff2_WD
# endif
# include "nbq.h"
!
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
!
# ifdef MPI
#  include "mpi_cpl.h"
      include 'mpif.h'
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif   

# include "compute_auxiliary_bounds.h"

# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif
!
# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
!
!*********************************************************************
! Set several approximations to speed-up computation
!*********************************************************************
!
# ifdef NBQ
!----------------------------------------------------------------------
! d./ds terms update frequency
!----------------------------------------------------------------------
!
#  ifdef NBQ_GRID_SLOW
#   define NSTEP_DS mod(iic,1).eq.1 .and. iif.eq.1
#  else
#   define NSTEP_DS .true. 
#  endif
!
!----------------------------------------------------------------------
! Grid update frequency (Hz ...)
!----------------------------------------------------------------------
!
#  ifdef NBQ_GRID_SLOW
#   define NSTEP_GRID .or. iif.eq.nfast
#  else
#   define NSTEP_GRID .or. mod(iif,1).eq.0
#  endif

# endif /* NBQ */
!
!*********************************************************************
!
!                 EXTERNAL (2D) FAST MODE PROCESSING
!
!  Compute external (2D) rhs terms rubar,rvbar of 2D fast-mode equations
!  using Generalized Forward-Backward AB3-AM4 algorithm, and update 
!  internal and external forcing terms for 3D fast mode equations: 
!  ru_int, ru_ext~(rufrc+rubar).
!
!  1- AB3 forward step for D,ubar,vbar
!  2- Advance zeta(m+1)
!  3- AM4 backward step for rubar,rvbar 
!  4- update rufrc,rvfrc
!  5- update ru_int,rv_int
!  6- make some backups
!
! Reference for Generalized FB scheme:
! ------------------------------------
! Shchepetkin, A.F., and J.C. McWilliams, 2009: Computational kernel 
! algorithms for fine-scale, multiprocess, longtime oceanic simulations. 
! Pp. 119–182 in Handbook of Numerical Analysis: Computational Methods 
! for the Atmosphere and Oceans. R.M. Teman and J.J. Tribbia, eds, 
! Elsevier Science.
!
!*********************************************************************

!
!=====================================================================
! Store rho.h at first slow and fast time-step
!=====================================================================
!
# ifdef NBQ
      if (FIRST_TIME_STEP.and.FIRST_FAST_STEP) then
         do k=1,N
           do j=JstrV-2,Jend+1
             do i=IstrU-2,Iend+1
               rho_bak(i,j,k)=rho(i,j,k) 
             enddo
            enddo
          enddo
      endif

      ! Extrapolation in time:
      do j=JstrR,JendR
          do k=1,N
            do i=IstrR,IendR
              rho_grd(i,j,k)=(rho(i,j,k)*1.5-rho_bak(i,j,k)*0.5)/rho0
            ! rho_grd(i,j,k)=rho(i,j,k)/rho0
            enddo
          enddo
      enddo
# endif
 
# ifdef M3FAST_REINIT
!
!  Reinitialise fast mode at each first fast step
!
      if (FIRST_FAST_STEP) then
        do k=1,N
          do j=Jstr-1,Jend+1
            do i=IstrU-1,Iend+1
              qdmu_nbq(i,j,k)=0.5*u(i,j,k,nstp)*(Hz(i,j,k)+Hz(i-1,j,k))
            enddo
          enddo 
        enddo
        do k=1,N
          do j=JstrV-1,Jend+1
            do i=Istr-1,Iend+1
              qdmv_nbq(i,j,k)=0.5*v(i,j,k,nstp)*(Hz(i,j,k)+Hz(i,j-1,k))
            enddo
          enddo 
        enddo
#  ifdef NBQ
        do j=JstrV-1,Jend
          do k=1,N-1
            do i=IstrU-1,Iend
              qdmw_nbq(i,j,k)=0.5*wz(i,j,k,nstp)*(Hz(i,j,k)+Hz(i,j,k+1))
            enddo
          enddo
          do i=IstrU-1,Iend
            qdmw_nbq(i,j,0)=0.5*wz(i,j,0,nstp)*Hz(i,j,1)
          enddo
          do i=IstrU-1,Iend
            qdmw_nbq(i,j,N)=0.5*wz(i,j,N,nstp)*Hz(i,j,N)
          enddo
        enddo
#  endif
      endif    ! FIRST_FAST_STEP
# endif


!
!=====================================================================
! AB3 Forward Step: compute total depth of water column and vertically
! --- ------- ----  integrated mass fluxes which are needed to compute 
! rhs terms of the barotropic momentum equations (rubar,rvbar).
!=====================================================================
!
!----------------------------------------------------------------------
!  Set indices to extrapolate (D,ubar,vbar) at m+1/2 (AB3)
!----------------------------------------------------------------------
!
      mybeta=0.281105 ! parameter for AB3 extrapolation

       if (FIRST_FAST_STEP) then     
                                      ! Meaning of temporal indices
        kbak=kstp                     ! ------- -- -------- -------
        kold=kstp                     ! m-2     m-1      m      m+1
        cff1= 1.0                     ! kold    kbak     kstp   knew
        cff2= 0.0
        cff3= 0.0
       elseif (FIRST_FAST_STEP+1) then  
        kbak=kstp-1                   ! AB2 forward scheme
        if (kbak.lt.1) kbak=4
        kold=kbak
        cff1= 1.5
        cff2=-0.5
        cff3= 0.0
      else                             ! AB3 forward scheme
        kbak=kstp-1 
        if (kbak.lt.1) kbak=4
        kold=kbak-1
        if (kold.lt.1) kold=4
        cff1= 1.5+mybeta
        cff2=-2.0*mybeta-0.5
        cff3= mybeta
      endif
!
# ifdef RVTK_DEBUG
C$OMP BARRIER
C$OMP MASTER
       call check_tab2d(zeta(:,:,kstp),'zeta step3d_fast #0','r')
C$OMP END MASTER
# endif  
!
!----------------------------------------------------------------------
!  Extrapolate (D,ubar,vbar) at m+1/2
!----------------------------------------------------------------------
!
! Total depth/mass at m+1/2
!
      do j=JstrV-2,Jend+1
        do i=IstrU-2,Iend+1
# ifdef NBQ_MASS
          Drhs(i,j)=cff1*(zeta(i,j,kstp)+h(i,j))*rhobar_nbq(i,j,kstp)
     &             +cff2*(zeta(i,j,kbak)+h(i,j))*rhobar_nbq(i,j,kbak)
     &             +cff3*(zeta(i,j,kold)+h(i,j))*rhobar_nbq(i,j,kold)
# else
          Drhs(i,j)=cff1*(zeta(i,j,kstp)+h(i,j))
     &             +cff2*(zeta(i,j,kbak)+h(i,j))
     &             +cff3*(zeta(i,j,kold)+h(i,j))
# endif /* NBQ_MASS */
        enddo
      enddo
      
!
! Depth-average ubar velocity at m+1/2
!
      do j=Jstr-1,Jend+1
        do i=IstrU-1,Iend+1
          urhs(i,j)=cff1*ubar(i,j,kstp) 
     &             +cff2*ubar(i,j,kbak)
     &             +cff3*ubar(i,j,kold)
# if defined MRL_WCI && defined MASKING
          urhs(i,j)=urhs(i,j)*umask(i,j)+ust2d(i,j)*(umask(i,j)-1.0)
# endif
          DUon(i,j)=0.5*(Drhs(i,j)+Drhs(i-1,j))*on_u(i,j)*( urhs(i,j)
# ifdef MRL_WCI
     &                                                   + ust2d(i,j)
# endif
     &                                                              )
        enddo
      enddo
!
! Depth-average ubar velocity at m+1/2
!
      do j=JstrV-1,Jend+1
        do i=Istr-1,Iend+1
          vrhs(i,j)=cff1*vbar(i,j,kstp)
     &             +cff2*vbar(i,j,kbak)
     &             +cff3*vbar(i,j,kold)
# if defined MRL_WCI && defined MASKING
          vrhs(i,j)=vrhs(i,j)*vmask(i,j)+vst2d(i,j)*(vmask(i,j)-1.0)
# endif
          DVom(i,j)=0.5*(Drhs(i,j)+Drhs(i,j-1))*om_v(i,j)*( vrhs(i,j)
# ifdef MRL_WCI
     &                                                   + vst2d(i,j)
# endif
     &                                                              )
        enddo
      enddo

# ifdef OBC_VOLCONS
      call set_DUV_bc_tile (Istr,Iend,Jstr,Jend, Drhs, DUon,DVom)
# endif

# ifdef RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
      call check_tab2d(zeta(:,:,kstp),'zeta step3d_fast #1','r')
C$OMP END MASTER
# endif 
!
!-----------------------------------------------------------------------
! Load RHS values into additional AGRIF shared array for nesting
!-----------------------------------------------------------------------
! 
# ifdef AGRIF
      irhox = Agrif_Irhox()
      irhoy = Agrif_Irhoy()
      irhot = Agrif_Irhot()

      if (FIRST_FAST_STEP) then
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            Zt_avg3(i,j,0)=zeta(i,j,kstp)       
          enddo
        enddo 
        do j=JstrR,JendR
          do i=Istr,IendR
            DU_avg3(i,j,0)=DUon(i,j)
          enddo
        enddo 
        do j=Jstr,JendR
          do i=IstrR,IendR
            DV_avg3(i,j,0)=DVom(i,j)
          enddo
        enddo 
      endif

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        Zt_avg3(START_2D_ARRAY,0))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        DU_avg3(START_2D_ARRAY,0))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        DV_avg3(START_2D_ARRAY,0))
#  endif

#  ifdef RVTK_DEBUG_ADVANCED
      if (.not.agrif_Root()) then
C$OMP BARRIER
C$OMP MASTER
        call check_tab2d(Zt_avg3(:,:,0),'Zt_avg3 (0) step3d_fast','r')
        call check_tab2d(DU_avg3(:,:,0),'DU_avg3 (0) step3d_fast','u')
        call check_tab2d(DV_avg3(:,:,0),'DV_avg3 (0) step3d_fast','v')
C$OMP END MASTER  
      endif
#  endif  
# endif /* AGRIF */
!
!=====================================================================
!
!  Compute zeta(m+1) -- first guest -- and update grid
!  
!=====================================================================
!
!----------------------------------------------------------------------
!  Get derived z grid variables at first fast step
!----------------------------------------------------------------------
!
# ifdef NBQ
      if (FIRST_FAST_STEP) then
        call grid_nbq_tile(Istr,Iend,Jstr,Jend,
     &                     Hzw_half_nbq_inv,   Hzr_half_nbq_inv,
     &                     Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v,
     &                     Hzu_half_qdmu,      Hzv_half_qdmv)
      endif
!
!-------------------------------------------------------------------
!  zeta(m+1) is derived here from the surface vertical velocity 
!  (surface characteristic relation) instead of depth-averaged 
!  conservation of mass. This satisfies dynamical coupling with 
!  the surface layer.
!
!  The surface kinematic relation is:
!
!  zeta(m+1) = zeta(m) + dtfast * [ wsurf_nbq(m)
!                                  - usurf_nbq(m)*dzeta/dx(m+0.5)
!                                  - vsurf_nbq(m)*dzeta/dy(m+0.5) ]
!
!  with zeta(m+0.5)=(1.5+beta)*zeta(m)-(0.5+2*beta)*zeta(m-1)
!                                             +beta*zeta(m-2) --> AB3
!-------------------------------------------------------------------
!  Computes surface velocities
!-------------------------------------------------------------------
!
      if (IstrU.le.Iend) then
        do j=Jstr,Jend
          do i=IstrU-1,Iend+1     
            usurf_nbq(i,j)=(qdmu_nbq(i,j,N)  
     &                     *2./(Hz(i,j,N)+Hz(i-1,j,N))
#  ifdef MRL_WCI
     &                     +ust(i,j,N)
#   ifdef WET_DRY
     &                      *umask_wet(i,j)
#  endif
# endif
     &                     )
#  ifdef MASKING
     &                     *umask(i,j) 
#  endif
          enddo 
        enddo 
      endif

      if (JstrV.le.Jend) then
        do j=JstrV-1,Jend+1
          do i=Istr,Iend     
            vsurf_nbq(i,j)=(qdmv_nbq(i,j,N)              
     &                     *2./(Hz(i,j,N)+Hz(i,j-1,N))
#  ifdef MRL_WCI
     &                     +vst(i,j,N)
#   ifdef WET_DRY
     &                      *vmask_wet(i,j)
#   endif
#  endif
     &                      )
#  ifdef MASKING
     &                     *vmask(i,j) 
#  endif
          enddo
        enddo 
      endif

      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
#  ifdef NBQ
          wsurf_nbq(i,j)=(qdmw_nbq(i,j,N)         
#   ifdef NBQ_MASS     
     &                       /(1.+rho_grd(i,j,N))
#   endif
     &                       *Hzw_half_nbq_inv(i,j,N)
#  else
          wsurf_nbq(i,j)=(We(i,j,N)*pm(i,j)*pn(i,j)
#  endif
#  ifdef MRL_WCI
     &                     +wst(i,j,N)*rmask_wet(i,j)
#  endif
     &                    )
#  ifdef MASKING
     &                       *rmask(i,j) 
#  endif
        enddo
      enddo 

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI  
      if (IstrU.le.Iend) then
        call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                          usurf_nbq(START_2D_ARRAY))
      endif
      if (JstrV.le.Jend) then
        call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                          vsurf_nbq(START_2D_ARRAY))
      endif
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        wsurf_nbq(START_2D_ARRAY))
#  endif
!
!-------------------------------------------------------------------
! Apply point sources for river runoff simulations
!-------------------------------------------------------------------
!
#  ifdef PSOURCE
      do is=1,Nsrc 
#   ifdef MPI
        i=Isrc_mpi(is,mynode)
        j=Jsrc_mpi(is,mynode)
#   else
        i=Isrc(is)
        j=Jsrc(is)
#   endif
        if (IstrR.le.i .and. i.le.IendR .and.
     &      JstrR.le.j .and. j.le.JendR) then
          if (Dsrc(is).eq.0) then
            urhs(i,j)=2.*Qbar(is)/( on_u(i,j)
     &                             *(Drhs(i-1,j)+Drhs(i,j)) )
            DUon(i,j)=Qbar(is)
            usurf_nbq(i,j)=2.*Qsrc(is,N)
     &                     /( on_u(i,j)*(Hz(i,j,N)+Hz(i-1,j,N)) )
          else
            vrhs(i,j)=2.*Qbar(is)/( om_v(i,j)
     &                             *(Drhs(i,j-1)+Drhs(i,j)) )
            DVom(i,j)=Qbar(is)
            vsurf_nbq(i,j)=2.*Qsrc(is,N)
     &                     /( om_v(i,j)*(Hz(i,j,N)+Hz(i,j-1,N)) )
          endif
        endif
      enddo
#  endif
!
!-------------------------------------------------------------------
! Advance zeta at m+1
!-------------------------------------------------------------------
!
#  define zab3 UFx
      if (FIRST_FAST_STEP) then
        cff1 = 1.0
        cff2 = 0.0
        cff3 = 0.0
      elseif (FIRST_FAST_STEP+1) then  ! AB2
        cff1 = 1.5
        cff2 =-0.5
        cff3 = 0.0
      else                             ! AB3
        cff1 = 1.5+mybeta
        cff2 =-2.0*mybeta-0.5
        cff3 = mybeta
      endif

      do j=JstrV-2,Jend+1
        do i=IstrU-2,Iend+1
          zab3(i,j) =   cff1 * zeta(i,j,kstp)
     &                + cff2 * zeta(i,j,kbak)
     &                + cff3 * zeta(i,j,kold)
        enddo
      enddo

      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta(i,j,knew)=(zeta(i,j,kstp) + dtfast*( wsurf_nbq(i,j) 
     &                     -0.5*(usurf_nbq(i  ,j)
     &                                     *(zab3(i  ,j)
     &                                      -zab3(i-1,j))*pm_u(i,j)
     &                          +usurf_nbq(i+1,j)
     &                                     *(zab3(i+1,j)
     &                                      -zab3(i  ,j))*pm_u(i+1,j) 
     &                          )
     &                     -0.5*(vsurf_nbq(i  ,j)
     &                                     *(zab3(i,j  ) 
     &                                      -zab3(i,j-1))*pn_v(i,j)
     &                          +vsurf_nbq(i,j+1)
     &                                     *(zab3(i,j+1)
     &                                      -zab3(i,j  ))*pn_v(i,j+1) 
     &                          )
     &                                                ) )
#  ifdef MASKING
     &                                                  *rmask(i,j)
#  endif
        enddo
      enddo
#  undef zab3



# else /* NBQ */
!
!********************************************************************
!
!  Advance zeta at m+1 from continuity for the hydrostatic case *****
!
!********************************************************************
!
!  First, apply point sources
!
#  ifdef PSOURCE
      do is=1,Nsrc 
#   ifdef MPI
        i=Isrc_mpi(is,mynode)
        j=Jsrc_mpi(is,mynode)
#   else
        i=Isrc(is)
        j=Jsrc(is)
#   endif
        if (IstrR.le.i .and. i.le.IendR .and.
     &      JstrR.le.j .and. j.le.JendR) then
          if (Dsrc(is).eq.0) then
            urhs(i,j)=2.*Qbar(is)/( on_u(i,j)
     &                             *(Drhs(i-1,j)+Drhs(i,j)) )
            DUon(i,j)=Qbar(is)
          else
            vrhs(i,j)=2.*Qbar(is)/( om_v(i,j)
     &                             *(Drhs(i,j-1)+Drhs(i,j)) )
            DVom(i,j)=Qbar(is)
          endif
        endif
      enddo
#  endif
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta(i,j,knew)=(zeta(i,j,kstp) + dtfast*pm(i,j)*pn(i,j)
     &                                   *(DUon(i,j)-DUon(i+1,j  )
     &                                    +DVom(i,j)-DVom(i  ,j+1)))
#  ifdef MASKING
     &                                                   *rmask(i,j)
#  endif
        enddo
      enddo
      

# endif /* NBQ */
!
!-----------------------------------------------------------------------
!  Add nudging terms
!-----------------------------------------------------------------------
!
# ifdef ZNUDGING
#  ifdef ZONAL_NUDGING
      if (FIRST_TIME_STEP .or. mod(iic,10).eq.0) then
        if (FIRST_FAST_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   zeta(START_2D_ARRAY,knew),zetazon)
        endif
      endif
      if (FIRST_TIME_STEP) then
        if (FIRST_FAST_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   ssh(START_2D_ARRAY),sshzon)
        endif
      endif
#  endif  /* ZONAL_NUDGING */
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta(i,j,knew)=zeta(i,j,knew) + dtfast*Znudgcof(i,j)
#  ifdef ZONAL_NUDGING
     &                                 *(sshzon(j)-zetazon(j))
#  else
     &                              *(ssh(i,j)-zeta(i,j,knew))
#  endif /* ZONAL_NUDGING */
#  ifdef MASKING
     &                                             *rmask(i,j)
#  endif
        enddo
      enddo
# endif /* ZNUDGING */
!
!-----------------------------------------------------------------------
!  Compute wet/dry masks
!-----------------------------------------------------------------------
!
! First: modify new free-surface to ensure that depth 
!        is > Dcrit in masked cells.
!
# if defined WET_DRY && defined MASKING
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          cff=0.5+SIGN(0.5,Dcrit(i,j)-h(i,j))
          zeta(i,j,knew)=zeta(i,j,knew)+ 
     &                   cff*(Dcrit(i,j)-h(i,j))*(1.-rmask(i,j))
        enddo
      enddo 
# endif
!
! Then compute wet/dry masks
!
# ifdef WET_DRY
      call wetdry_tile (Istr,Iend,Jstr,Jend)
# endif
!
!-----------------------------------------------------------------------
!  Set boundary conditions for the free-surface
!-----------------------------------------------------------------------
!
      call zetabc_tile (Istr,Iend,Jstr,Jend)
!
!-----------------------------------------------------------------------
!  Perform exchanges
!-----------------------------------------------------------------------
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        zeta(START_2D_ARRAY,knew))
# endif
!
!-----------------------------------------------------------------------
!  Update vertical grid
!
!  As soon as zeta(m+1) is known, the grid can be updated at fast
!  step m+1. In PERF option, the grid is updated at lower frequency 
!  and only during re-evaluation of zeta(m+1).
!
!  Caution: Hz_bak must be assigned only once in set_depth. The
!  following code must thus be consistant with set_depth routine.
!-----------------------------------------------------------------------
!
# if defined NBQ && !defined NBQ_GRID_SLOW
      if ( (FIRST_TIME_STEP .and. FIRST_FAST_STEP) NSTEP_GRID ) then
!
! Update main grid parameters
!
        call set_depth_tile(Istr,Iend,Jstr,Jend)
!
! Update derived grid variables
!
        call grid_nbq_tile(Istr,Iend,Jstr,Jend,
     &                     Hzw_half_nbq_inv,   Hzr_half_nbq_inv,
     &                     Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v,
     &                     Hzu_half_qdmu,      Hzv_half_qdmv )
      endif
# endif /* NBQ_GRID_SLOW */

# ifndef NOT_NBQ_AM4
! 
!=====================================================================
! AM4 backward step: compute depth-averaged RHS of 2D fast-mode 
! --- -------- ----  momentum equations
!=====================================================================
!
!----------------------------------------------------------------------
! Surface pressure gradient
!----------------------------------------------------------------------
!
! Interpolate zeta fields half-step backward (AM4) for the subsequent 
! computation of barotropic pressure-gradient
!
      myalpha   = 0.25
      myepsilon = 0.00976186 - 0.13451357*myalpha
      mygamma   = 0.08344500 - 0.51358400*myalpha
!
      if (FIRST_FAST_STEP) then
        cff0=0.                  !---> Compute pressure-gradient
        cff1=1.                  !     terms using just zeta(:,:,kstp)
        cff2=0.
        cff3=0.
      elseif (FIRST_FAST_STEP+1) then
        cff0= 1.0833333333333    ! AM3 backward scheme
        cff1=-0.1666666666666    ! with coefficients chosen for
        cff2= 0.0833333333333    ! maximum stability, while maintaining
        cff3= 0.                 ! third-accuracy; alpha_max=1.73
      else
        cff0=0.5+2.*myepsilon+mygamma+2.*myalpha  ! AM4 backward scheme
        cff1=1.-cff0-mygamma-myepsilon            ! with implicit diffusion
        cff2=mygamma                              ! given by myalpha
        cff3=myepsilon
      endif
# endif
!
# define zwrk    UFx
# define rzeta   UFe
# define rzeta2  VFe
# define rzetaSA VFx

      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
# ifdef NOT_NBQ_AM4
          zwrk(i,j)=zeta(i,j,kstp)
# else
          zwrk(i,j)=(cff0*zeta(i,j,knew)+cff1*zeta(i,j,kstp)
     &              +cff2*zeta(i,j,kbak)+cff3*zeta(i,j,kold))
# endif
# ifdef NBQ_MASS
     &                                *rhobar_nbq(i,j,kstp) 
# endif
# ifdef VAR_RHO_2D
          rzeta(i,j)=(1.+rhoS(i,j))*zwrk(i,j)
          rzeta2(i,j)=rzeta(i,j)*zwrk(i,j)
          rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
# else
          rzeta(i,j)=zwrk(i,j)
          rzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
# endif
        enddo
      enddo
!
! Compute surface pressure gradient
!
      cff=0.5*g
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=cff*on_u(i,j)*( (h(i-1,j)+h(i,j))*(rzeta(i-1,j)
     &                        -rzeta(i,j)) +rzeta2(i-1,j)-rzeta2(i,j)
# ifdef VAR_RHO_2D
     &              +(h(i-1,j)-h(i,j))*( rzetaSA(i-1,j)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i-1,j)-rhoA(i,j))
     &                                      *(zwrk(i-1,j)-zwrk(i,j)))
# endif
# ifdef MRL_WCI
     &                  + ( h(i-1,j)+h(i,j)+rzeta(i-1,j)+rzeta(i,j) )
     &                                       *( sup(i,j)-sup(i-1,j) )
# endif
     &                                                              )
 
          rvbar(i,j)=cff*om_v(i,j)*( (h(i,j-1)+h(i,j))*(rzeta(i,j-1)
     &                        -rzeta(i,j)) +rzeta2(i,j-1)-rzeta2(i,j)
# ifdef VAR_RHO_2D
     &              +(h(i,j-1)-h(i,j))*( rzetaSA(i,j-1)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i,j-1)-rhoA(i,j))
     &                                      *(zwrk(i,j-1)-zwrk(i,j)))
# endif
# ifdef MRL_WCI
     &                  + ( h(i,j-1)+h(i,j)+rzeta(i,j-1)+rzeta(i,j) )
     &                                       *( sup(i,j)-sup(i,j-1) )
# endif
     &                                                              )
        enddo
      enddo            !--> discard  zwrk, rzeta, rzeta2, rzetaSA

# undef rzetaSA
# undef rzeta2
# undef rzeta
# undef zwrk
!
# ifdef UV_ADV
!----------------------------------------------------------------------
! Compute horizontal advection terms for momentum equations (2D only)
!-------- ---------- --------- ----- --- -------- --------- --- -----
!
! Centered second order advection scheme
!
! Numerical diffusion of momentum is implicitely added through 3D
! forcing of advection in rufrc and rvfrc (i.e., diffusion is
! at slow time scale)
!
! NOTE: mathematically necessary (minimal) index ranges for momentum-
! flux components are 
!
!      UFx(IstrU-1:Iend,Jstr:Jend)   VFx(Istr:Iend+1,JstrV:Jend)
!      UFe(IstrU:Iend,Jstr:Jend+1)   VFe(Istr,Iend,JstrV-1,Jend)
!
! however, for computational efficiency, these ranges are
! unified by suppressing U,V-suffices in order to allow fusion of the
! consecutive loops. This leads to slight increase of the redundant
! computations near western and southern boundaries in non-periodic
! directions.
!----------------------------------------------------------------------
!
      do j=Jstr,Jend
        do i=Istr-1,Iend
          UFx(i,j)=0.25*(DUon(i,j)+DUon(i+1,j))
     &                 *(urhs(i,j)+urhs(i+1,j))

          VFx(i+1,j)=0.25*(DUon(i+1,j)+DUon(i+1,j-1))
     &                   *(vrhs(i+1,j)+vrhs(i,j))
#  ifdef MASKING
     &                                 *pmask(i+1,j)
#  endif
        enddo
      enddo
      do j=Jstr-1,Jend
        do i=Istr,Iend
          VFe(i,j)=0.25*(DVom(i,j)+DVom(i,j+1))
     &                 *(vrhs(i,j)+vrhs(i,j+1))

          UFe(i,j+1)=0.25*(DVom(i,j+1)+DVom(i-1,j+1))
     &                   *(urhs(i,j+1)+urhs(i,j))
#  ifdef MASKING
     &                                 *pmask(i,j+1)
#  endif
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=rubar(i,j)-UFx(i,j)+UFx(i-1,j)
     &                         -UFe(i,j+1)+UFe(i,j)

          rvbar(i,j)=rvbar(i,j)-VFx(i+1,j)+VFx(i,j)
     &                         -VFe(i,j)+VFe(i,j-1)
        enddo
      enddo    !--> discard UFx,VFe,UFe,VFx, DUon,DVom
# endif /* UV_ADV */
!
!-----------------------------------------------------------------------
! Compute Coriolis (2D and 3D) term and advective curvilinear metric
! terms (2D only).
!-----------------------------------------------------------------------
!
# if defined UV_COR || (defined CURVGRID && defined UV_ADV)
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          cff=Drhs(i,j)*(
#  ifdef UV_COR
     &                   fomn(i,j)
#  endif
#  if (defined CURVGRID && defined UV_ADV)
     &          +0.5*( dndx(i,j)*(vrhs(i,j)+vrhs(i,j+1))
     &                -dmde(i,j)*(urhs(i,j)+urhs(i+1,j)))
#  endif
     &                   )
#  ifdef MRL_WCI
#   if defined CURVGRID && defined UV_ADV
          cff1 = Drhs(i,j)*(
     &    0.5*( dndx(i,j)*(vst2d(i,j)+vst2d(i,j+1))
     &         -dmde(i,j)*(ust2d(i,j)+ust2d(i+1,j)) ))
#   else
          cff1 = 0.0
#   endif
          UFx(i,j)=(cff+cff1)*(vrhs(i,j)+vrhs(i,j+1))
     &                 +cff*(vst2d(i,j)+vst2d(i,j+1))
          VFe(i,j)=(cff+cff1)*(urhs(i,j)+urhs(i+1,j))
     &                 +cff*(ust2d(i,j)+ust2d(i+1,j))
#  else
          UFx(i,j)=cff*(vrhs(i,j)+vrhs(i,j+1))
          VFe(i,j)=cff*(urhs(i,j)+urhs(i+1,j))
#  endif 
        enddo
      enddo
      do j=Jstr,Jend
        do i=IstrU,Iend
          rubar(i,j)=rubar(i,j)+0.25*(UFx(i,j)+UFx(i-1,j))
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          rvbar(i,j)=rvbar(i,j)-0.25*(VFe(i,j)+VFe(i,j-1))
        enddo 
      enddo
# endif /* UV_COR */
!
!-----------------------------------------------------------------------
! Linear and/or quadratic bottom stress.
!-----------------------------------------------------------------------
!
# ifndef BSTRESS_FAST
#  ifdef BBL
        do j=Jstr,Jend
          do i=IstrU,Iend
            rubar(i,j)=rubar(i,j) - bustr(i,j)
     &                             *om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rvbar(i,j)=rvbar(i,j) - bvstr(i,j)
     &                             *om_v(i,j)*on_v(i,j)
          enddo
        enddo
#  else
      if (rdrg2.gt.0.) then
        do j=Jstr,Jend
          do i=IstrU,Iend
            cff=0.25*( vbar(i  ,j,kstp)+vbar(i  ,j+1,kstp)
     &                +vbar(i-1,j,kstp)+vbar(i-1,j+1,kstp))
 
            rubar(i,j)=rubar(i,j) - ubar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(ubar(i,j,kstp)*ubar(i,j,kstp)+cff*cff)
     &                               )*om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            cff=0.25*( ubar(i,j  ,kstp)+ubar(i+1,j  ,kstp)
     &                +ubar(i,j-1,kstp)+ubar(i+1,j-1,kstp))
 
            rvbar(i,j)=rvbar(i,j) - vbar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(cff*cff+vbar(i,j,kstp)*vbar(i,j,kstp))
     &                               )*om_v(i,j)*on_v(i,j)
          enddo
        enddo
      else if (rdrg.gt.0.0) then
        do j=Jstr,Jend
          do i=IstrU,Iend
            rubar(i,j)=rubar(i,j) - rdrg*ubar(i,j,kstp)
     &                             *om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rvbar(i,j)=rvbar(i,j) - rdrg*vbar(i,j,kstp)
     &                             *om_v(i,j)*on_v(i,j)
          enddo
        enddo
      endif
#  endif
# endif /* BSTRESS_FAST */
!
!-----------------------------------------------------------------------
! Add 2D vortex-force terms combined with advection terms
!-----------------------------------------------------------------------
!
# ifdef MRL_WCI
      do j=Jstr,Jend
        do i=IstrU,Iend
          vstu = 0.25*( vst2d(i  ,j)+vst2d(i  ,j+1)
     &                 +vst2d(i-1,j)+vst2d(i-1,j+1) )
          dudx = 0.50*( urhs(i+1,j)-urhs(i-1,j  ) )
          dvdx = 0.50*( vrhs(i  ,j)-vrhs(i-1,j  )
     &                 +vrhs(i,j+1)-vrhs(i-1,j+1) )
          rubar(i,j) = rubar(i,j) + 0.5*on_u(i,j)*
     &                      ( Drhs(i-1,j)+Drhs(i,j) )
     &               *( ust2d(i,j)*dudx + vstu*dvdx )
        enddo
      enddo

      do j=JstrV,Jend
        do i=Istr,Iend
          ustv = 0.25*( ust2d(i,j  )+ust2d(i+1,j  )
     &                 +ust2d(i,j-1)+ust2d(i+1,j-1) )
          dude = 0.50*( urhs(i,j  )-urhs(i  ,j-1)
     &                 +urhs(i+1,j)-urhs(i+1,j-1) )
          dvde = 0.50*( vrhs(i,j+1)-vrhs(i  ,j-1) )
          rvbar(i,j) = rvbar(i,j) + 0.5*om_v(i,j)*
     &                      ( Drhs(i,j-1)+Drhs(i,j) )
     &               *( ustv*dude + vst2d(i,j)*dvde )
        enddo
      enddo
# endif
!
!=====================================================================
! Compute rufrc & rvfrc: internal mode forcing for barotropic fast mode
!
! During the first fast time step convert rufrc & fvfrc into forcing
! terms by subtracting the fast-time "rubar" and "rvbar" from them;
! These forcing terms are then extrapolated forward in time using
! optimized Adams-Bashforth weights, so that the resultant rufrc
! and rvfrc are centered effectively at time n+1/2. From now on,
! these newly computed forcing terms will remain constant during
! the fast time stepping and will be added to "rubar" and "rvbar"
! during all subsequent fast time steps.
!=====================================================================
!
      if (FIRST_FAST_STEP) then
!
        if (FIRST_TIME_STEP) then
          cff3=0.                        ! This version is designed
          cff2=0.                        ! for coupling during 3D
          cff1=1.                        ! predictor sub-step: here
        elseif (FIRST_TIME_STEP+1) then  ! forcing term "rufrc" is
          cff3= 0.                       ! computed as instantaneous
          cff2=-0.5                      ! value at 3D time step
          cff1= 1.5                      ! "nstp" first, and then
        else                             ! extrapolated half-step
          cff3= 0.281105                 ! forward using  AM3-like
          cff2=-0.5-2.*cff3              ! weights optimized for
          cff1= 1.5+cff3                 ! maximum stability (with
        endif                            ! special care for startup)
        do j=Jstr,Jend
         do i=IstrU,Iend
           cff=rufrc(i,j)-rubar(i,j)
           rufrc(i,j)=cff1*cff + cff2*rufrc_bak(i,j,3-nstp)
     &                         + cff3*rufrc_bak(i,j,nstp)
           rufrc_bak(i,j,nstp)=cff
         enddo
        enddo
        do j=JstrV,Jend
         do i=Istr,Iend
           cff=rvfrc(i,j)-rvbar(i,j)
           rvfrc(i,j)=cff1*cff + cff2*rvfrc_bak(i,j,3-nstp)
     &                         + cff3*rvfrc_bak(i,j,nstp)
           rvfrc_bak(i,j,nstp)=cff
         enddo
        enddo
      endif   !<-- FIRST_FAST_STEP
!
!=====================================================================
!   Update internal and external forcing terms for NBQ mode
!
!   Compute external forcing terms ru_ext_nbq and updated internal 
!   forcing terms ru_int_nbq for NBQ equations
!
!   ru_int_nbq     : RHS (3D) ( *mask & 2D correction)
!   ru_ext_nbq     : RHS (2D)
!   ru_ext_nbq_old : RHS (2D) at previous time-step
!   ru_ext_nbq_sum : time-integrated RHS (2D)
!=====================================================================
!
!-----------------------------------------------------------------------
!  First fast time step only
!-----------------------------------------------------------------------
!
      if (FIRST_FAST_STEP) then
!
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru_ext_nbq_sum(i,j)=0.
            ru_ext_nbq_old(i,j)=0.
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rv_ext_nbq_sum(i,j)=0.
            rv_ext_nbq_old(i,j)=0.
          enddo
        enddo 
# ifdef MASKING
        do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend
              ru_int_nbq(i,j,k)=ru_int_nbq(i,j,k)*umask(i,j)
            enddo
          enddo
        enddo
        do k=1,N
          do j=JstrV,Jend
            do i=Istr,Iend
              rv_int_nbq(i,j,k)=rv_int_nbq(i,j,k)*vmask(i,j)
            enddo
          enddo
        enddo
# endif
      endif ! FIRST_FAST_STEP
!
!-----------------------------------------------------------------------
!  All fast time steps
!-----------------------------------------------------------------------
!
# ifndef NBQ
      if (mod(iif-1,inc_faststep).eq.0) then
        nb_faststep = 0
        do j=Jstr,Jend
          do i=IstrU,Iend
            rubar_sum(i,j) = 0.
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rvbar_sum(i,j) = 0.
          enddo
        enddo
      endif
      nb_faststep = nb_faststep + 1
# endif

# define ru_ext_nbq UFx
      do j=Jstr,Jend
        do i=IstrU,Iend
# ifdef NBQ
          ru_ext_nbq(i,j)=(rufrc(i,j)+rubar(i,j))
     &                                   *pm_u(i,j)*pn_u(i,j)
     &                                  /(Drhs(i,j)+Drhs(i-1,j))
#  ifdef MASKING
     &                                                *umask(i,j)
#  endif
          ru_ext_nbq_old(i,j)=ru_ext_nbq(i,j)-ru_ext_nbq_old(i,j)
# else
          rubar_nbq(i,j)=(rufrc(i,j)+rubar(i,j))*pm_u(i,j)*pn_u(i,j)
#  ifdef MASKING
     &                                                *umask(i,j)
#  endif
          ru_ext_nbq(i,j)=rubar_nbq(i,j)/(Drhs(i,j)+Drhs(i-1,j))
          rubar_sum(i,j)=rubar_sum(i,j)+ru_ext_nbq(i,j)
# endif
          ru_ext_nbq_sum(i,j)=ru_ext_nbq_sum(i,j)+ru_ext_nbq(i,j)
        enddo
      enddo
# ifdef NBQ     
      do k=1,N
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru_int_nbq(i,j,k)=ru_int_nbq(i,j,k)+
     &                        ru_ext_nbq_old(i,j)
     &                        *(Hz(i-1,j,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=IstrU,Iend
          ru_ext_nbq_old(i,j)=ru_ext_nbq(i,j)
        enddo
      enddo
# endif
# undef ru_ext_nbq

# define rv_ext_nbq UFx
      do j=JstrV,Jend
        do i=Istr,Iend
# ifdef NBQ
          rv_ext_nbq(i,j)=(rvfrc(i,j)+rvbar(i,j))
     &                                   *pm_v(i,j)*pn_v(i,j)
     &                                   /(Drhs(i,j)+Drhs(i,j-1))
#  ifdef MASKING
     &                                                *vmask(i,j)
#  endif
          rv_ext_nbq_old(i,j)=rv_ext_nbq(i,j)-rv_ext_nbq_old(i,j)
# else
          rvbar_nbq(i,j)=(rvfrc(i,j)+rvbar(i,j))*pm_v(i,j)*pn_v(i,j)
#  ifdef MASKING
     &                                                *vmask(i,j)
#  endif
          rv_ext_nbq(i,j)=rvbar_nbq(i,j)/(Drhs(i,j)+Drhs(i,j-1))
          rvbar_sum(i,j)=rvbar_sum(i,j)+rv_ext_nbq(i,j)
# endif
          rv_ext_nbq_sum(i,j)=rv_ext_nbq_sum(i,j)+rv_ext_nbq(i,j)
        enddo
      enddo
# ifdef NBQ        
      do k=1,N
        do j=JstrV,Jend
          do i=Istr,Iend
            rv_int_nbq(i,j,k)=rv_int_nbq(i,j,k)+
     &                        rv_ext_nbq_old(i,j)
     &                        *(Hz(i,j-1,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          rv_ext_nbq_old(i,j)=rv_ext_nbq(i,j)
        enddo
      enddo
# endif
# undef rv_ext_nbq 

# ifdef NBQ
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,  
     &                        ru_int_nbq(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,  
     &                        rv_int_nbq(START_2D_ARRAY,1))
#  endif 
# endif

# if defined RVTK_DEBUG_ADVANCED && defined NBQ
C$OMP BARRIER
C$OMP MASTER
       call check_tab3d(ru_int_nbq,'ru_int_nbq (A)','u')
       call check_tab3d(rv_int_nbq,'rv_int_nbq (A)','v')
C$OMP END MASTER     
# endif  
!
!=====================================================================
!  Initializations and backups for 3D NBQ equations
!=====================================================================
!
!-----------------------------------------------------------------------
!  Initializations of ru_nbq_avg2 (as qdmu_nbq) at first fast step
!  It is used at last fast step for computation of new ru_nbq_avg2
!-----------------------------------------------------------------------
!
      if (FIRST_FAST_STEP) then
        do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend
              ru_nbq_avg2(i,j,k)=qdmu_nbq(i,j,k)
            enddo
          enddo 
        enddo
        do k=1,N
          do j=JstrV,Jend
            do i=Istr,Iend
              rv_nbq_avg2(i,j,k)=qdmv_nbq(i,j,k)
            enddo
          enddo 
        enddo
# ifdef NBQ
        do k=0,N
          do j=Jstr,Jend
            do i=Istr,Iend
              rw_nbq_avg2(i,j,k)=qdmw_nbq(i,j,k)
            enddo
          enddo 
        enddo
# endif
      endif    ! FIRST_FAST_STEP
!
!-----------------------------------------------------------------------
!  zw_nbq backups
!-----------------------------------------------------------------------
!
# ifdef NBQ
      if (FIRST_TIME_STEP .and. FIRST_FAST_STEP) then
        do j=JstrR,JendR
          do k=0,N
            do i=IstrR,IendR
              zw_nbq (i,j,k,:)=z_w(i,j,k)
            enddo
          enddo
        enddo
      endif
# endif
!
!-------------------------------------------------------------------
!  Acoustic wave emission
!-------------------------------------------------------------------
!
# ifdef ACOUSTIC 
      if (FIRST_TIME_STEP .and. FIRST_FAST_STEP) then
        period_exp = 0.025/2.
        for_a_exp  = 2.5
        amp_exp    = 1.e-2
        hmax_exp   = 128.
        dg_exp     = 2.
        time_nbq   = 0. 
      endif
# endif
!
!------------------------------------------------------------------
!  Fast bottom friction
!
!  Set bottom stress using logarithmic or linear 
!  and/or quadratic formulation. 
!-------------------------------------------------------------------
!
# if defined BSTRESS_FAST && !defined BBL
#  define UBOT UFx
#  define VBOT VFe

      if (mod(iif-1,inc_faststep).eq.0) then

        if (maxval(Zob).ne.0.) then
          do j=JstrV-1,Jend+1
            do i=IstrU-1,Iend+1
              UBOT(i,j)=2.*qdmu_nbq(i,j,1)/(Hz(i,j,1)+Hz(i-1,j,1))
              VBOT(i,j)=2.*qdmv_nbq(i,j,1)/(Hz(i,j,1)+Hz(i,j-1,1))
            enddo
          enddo
          do j=JstrV-1,Jend
            do i=IstrU-1,Iend
              !cff1=MAX(z_r(i,j,1)-z_w(i,j,0),Zob(i,j)+1.e-4)
              cff1=z_r(i,j,1)-z_w(i,j,0)
              cff=vonKar/LOG(cff1/Zob(i,j))
              work(i,j)=MIN(Cdb_max,MAX(Cdb_min,cff*cff))
            enddo
          enddo
          do j=Jstr,Jend
            do i=IstrU,Iend
              cff=0.25*(VBOT(i  ,j)+VBOT(i  ,j+1)+
     &                  VBOT(i-1,j)+VBOT(i-1,j+1))
              bustr(i,j)=0.5*(work(i-1,j)+work(i,j))*UBOT(i,j)*
     &                   SQRT(UBOT(i,j)*UBOT(i,j)+cff*cff)
            enddo
          enddo
          do j=JstrV,Jend
            do i=Istr,Iend
              cff=0.25*(UBOT(i,j  )+UBOT(i+1,j  )+
     &                  UBOT(i,j-1)+UBOT(i+1,j-1))
              bvstr(i,j)=0.5*(work(i,j-1)+work(i,j))*VBOT(i,j)*
     &                   SQRT(VBOT(i,j)*VBOT(i,j)+cff*cff)
            enddo
          enddo
          
        elseif (rdrg2.gt.0.) then
          do j=JstrV-1,Jend+1
            do i=IstrU-1,Iend+1
              UBOT(i,j)=2.*qdmu_nbq(i,j,1)/(Hz(i,j,1)+Hz(i-1,j,1))
              VBOT(i,j)=2.*qdmv_nbq(i,j,1)/(Hz(i,j,1)+Hz(i,j-1,1))
            enddo
          enddo
          do j=JstrV,Jend
            do i=Istr,Iend
              cff=0.25*(VBOT(i  ,j)+VBOT(i  ,j+1)+
     &                  VBOT(i-1,j)+VBOT(i-1,j+1))
              bustr(i,j)=rdrg2*UBOT(i,j)*
     &                   SQRT(UBOT(i,j)*UBOT(i,j)+cff*cff)
            enddo
          enddo
          do j=Jstr,Jend
            do i=IstrU,Iend
              cff=0.25*(UBOT(i,j  )+UBOT(i+1,j  )+
     &                  UBOT(i,j-1)+UBOT(i+1,j-1))
              bvstr(i,j)=rdrg2*VBOT(i,j)*
     &                    SQRT(VBOT(i,j)*VBOT(i,j)+cff*cff)
            enddo
          enddo
        else
          do j=Jstr,Jend
            do i=IstrU,Iend
              bustr(i,j)=rdrg*2.*qdmu_nbq(i,j,1)/(Hz(i,j,1)+Hz(i-1,j,1))
            enddo
          enddo
          do j=JstrV,Jend
            do i=Istr,Iend
              bvstr(i,j)=rdrg*2.*qdmv_nbq(i,j,1)/(Hz(i,j,1)+Hz(i,j-1,1))
            enddo
          enddo
        endif

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
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=Jstr,Jend
          bustr(IstrU-1,j)=bustr(IstrU,j)
        END DO
        DO j=JstrV,Jend
          bvstr(Istr-1,j)=bvstr(Istr,j)
        END DO
      END IF
#  endif
#  ifndef NS_PERIODIC
      IF (NORTHERN_EDGE) THEN
        DO i=IstrU,Iend
          bustr(i,Jend+1)=bustr(i,Jend)
        END DO
        DO i=Istr,Iend
          bvstr(i,Jend+1)=bvstr(i,Jend)
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=IstrU,Iend
          bustr(i,Jstr-1)=bustr(i,Jstr)
        END DO
        DO i=Istr,Iend
          bvstr(i,JstrV-1)=bvstr(i,JstrV)
        END DO
      END IF
#  endif
#  if !defined EW_PERIODIC && !defined NS_PERIODIC
      IF (SOUTHERN_EDGE.and.WESTERN_EDGE) THEN
        bustr(Istr,Jstr-1)=0.5*(bustr(Istr+1,Jstr-1)+bustr(Istr,Jstr))
        bvstr(Istr-1,Jstr)=0.5*(bvstr(Istr,Jstr)+bvstr(Istr-1,Jstr+1))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        bustr(Iend+1,Jstr-1)=0.5*(bustr(Iend+1,Jstr)+bustr(Iend,Jstr-1))
        bvstr(Iend+1,Jstr)=0.5*(bvstr(Iend+1,Jstr+1)+bvstr(Iend,Jstr))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        bustr(Istr,Jend+1)=0.5*(bustr(Istr,Jend)+bustr(Istr+1,Jend+1))
        bvstr(Istr-1,Jend+1)=0.5*(bvstr(Istr-1,Jend)+bvstr(Istr,Jend+1))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        bustr(Iend+1,Jend+1)=0.5*(bustr(Iend+1,Jend)+bustr(Iend,Jend+1))
        bvstr(Iend+1,Jend+1)=0.5*(bvstr(Iend+1,Jend)+bvstr(Iend,Jend+1))
      END IF
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,bustr(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,bvstr(START_2D_ARRAY))
#  endif


      endif
#  undef UBOT
#  undef VBOT
# endif
!
!------------------------------------------------------------------
!  Implicit part: system setup
!-------------------------------------------------------------------
!    
      do j=Jstr,Jend
        do i=Istr,Iend
# ifdef NBQ
          work(i,j)=pm(i,j)*pn(i,j)
# else
          if (LAST_FAST_STEP) work(i,j)=pm(i,j)*pn(i,j)
# endif
          DU_nbq(i,j)=0.
          DV_nbq(i,j)=0.
        enddo
      enddo
!	
!------------------------------------------------------------------
!  Store qdmw_nbq into working array
!------------------------------------------------------------------
!
# ifdef NBQ
      if (LAST_FAST_STEP) then
        do k=0,N 
          do j=Jstr,Jend             
            do i=Istr,Iend
               rw_nbq(i,j,k)=qdmw_nbq(i,j,k) 
            enddo
          enddo
        enddo
      endif
# endif
!	
!------------------------------------------------------------------
!  Recover Hz at first fast step (if final Hz correction needed)
!------------------------------------------------------------------
!
# ifdef NBQ_HZCORRECT
      if (iic.gt.ntstart .and. FIRST_FAST_STEP) then
        do k=1,N
          do j=JstrV-2,Jend+1
            do i=IstrU-2,Iend+1
              Hz(i,j,k) = Hz_correct(i,j,k)
            enddo
          enddo
        enddo
      endif
# endif
!
!----------------------------------------------------------------------
!  Store boundary values of nbq variables at previous 
!  time-step for use in radiation boundary conditions
!----------------------------------------------------------------------
!
# if defined OBC_NBQ && defined OBC_NBQORLANSKI
      call nbq_bry_store_tile (Istr,Iend,Jstr,Jend)
# endif
!
!
!***********************************************************************
!
!                  SOLVE FAST MODE 3D NBQ EQUATIONS
!
!***********************************************************************
!
!  Fast mode 3D momentum and mass-conservation equations can now be solved.
!
!  W-momentum equation is solved with explicit or implicit methods:
!  - Explicit scheme: w-momentum is updated right after (and the same 
!                     way as) u- and v-momentum.
!  - Implicit scheme: horizontal component of divergence is first 
!                     precomputed (as required by fast-mode mass 
!                     conservation) before tridiagonal Gauss Elimination 
!                     is carried out for qdmw_nbq(m).
!
!  Compressible pressure-force and second viscosity are calculated 
!  in thetadiv_nbq. Caution: this variable contains theta in the 
!  first part of the algorithm and momentum divergence in the remaining.
!
!  A Forward-backward scheme is implemented:
!  - Explicit scheme: Forward:  zeta, qdmu_nbq, qdmw_nbq. 
!                     Backward:       rho_nbq.
!  - Implicit scheme: Forward:  zeta, qdmu_nbq. 
!                     Backward:       qdmw_nbq, rho_nbq.
!
!***********************************************************************
!
# ifdef UV_COR_NT
!=====================================================================
!  Compute Non-traditional Coriolis pseudo-force
!
!    ntcoru = -e W cos(a) * Hz      [m2/s2]
!    ntcorv = +e W sin(a) * Hz
!    ntcorw = +e (U cos(a) - V sin(a) ) * Hz
!  
!  with e = 2 Omega cos(Phi)
!       a = angle between North and meridional grid axis
!
!=====================================================================
!
#  define Huw ntcoru
#  define Hvw ntcorv
      do j=Jstr,Jend
        do i=Istr,Iend+1
          Huw(i,j,0)=0.5*qdmu_nbq(i,j,1)
          Huw(i,j,N)=0.5*qdmu_nbq(i,j,N)
        enddo
      enddo
      do j=Jstr,Jend+1
        do i=Istr,Iend
          Hvw(i,j,0)=0.5*qdmv_nbq(i,j,1)
          Hvw(i,j,N)=0.5*qdmv_nbq(i,j,N)
        enddo
      enddo
      do k=1,N-1
        do j=Jstr,Jend
          do i=Istr,Iend+1
            Huw(i,j,k)=0.5*(qdmu_nbq(i,j,k)+qdmu_nbq(i,j,k+1))
          enddo
        enddo
        do j=Jstr,Jend+1
          do i=Istr,Iend
            Hvw(i,j,k)=0.5*(qdmv_nbq(i,j,k)+qdmv_nbq(i,j,k+1))
          enddo
        enddo
      enddo
      do k=0,N
        do j=Jstr,Jend
          do i=Istr,Iend
            ntcorw(i,j,k) = 0.5*e(i,j)*(
     &          cosa(i,j)*(Huw(i,j,k)+Huw(i+1,j,k))
     &        - sina(i,j)*(Hvw(i,j,k)+Hvw(i,j+1,k))
     &                                 )
          enddo
        enddo
      enddo
#  undef Huw
#  undef Hvw
!
      do k=1,N
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
            UFx(i,j)=-0.5*e(i,j)*cosa(i,j)
     &                            *( qdmw_nbq(i,j,k-1)
     &                              +qdmw_nbq(i,j,k  )) 
            VFe(i,j)=+0.5*e(i,j)*sina(i,j)
     &                            *( qdmw_nbq(i,j,k-1)
     &                              +qdmw_nbq(i,j,k)  )
          enddo
        enddo
        do j=Jstr,Jend
          do i=IstrU,Iend
            ntcoru(i,j,k)=0.5*(UFx(i,j)+UFx(i-1,j))
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            ntcorv(i,j,k)=0.5*(VFe(i,j)+VFe(i,j-1))
          enddo
        enddo
      enddo
# endif
!
!=====================================================================
!  Compute "Pressure - Viscosity" component (theta)
!=====================================================================
!
# ifdef NBQ
      do k=1,N
        do j=JstrV-2,Jend+1
          do i=IstrU-2,Iend+1
            thetadiv_nbq(i,j,k)=( -visc2_nbq*( thetadiv_nbq(i,j,k)
     &                                       +thetadiv3_nbq(i,j,k))
     &                         +soundspeed2_nbq(i,j)*rho_nbq(i,j,k) )
     &                                     *Hzr_half_nbq_inv(i,j,k)  
          enddo
        enddo
      enddo
!
!=====================================================================
!  Integrate fast-mode momentum equations
!
!  dqdm/dt=dtfast*(Compressible pressure force + second viscosity
!                  + gravity + non-NBQ forces + nudging)
!
!  If explicit scheme: all (x,y,z) mom components are done here
!  If implicit scheme: z-mom component is done after
!=====================================================================
!
!-----------------------------------------------------------------------
!  Pressure-Viscosity forces in XI- and ETA-Directions
!-----------------------------------------------------------------------
!
      k2 = 1
      do k=0,N
        k1=k2
        k2=3-k1

#  ifdef NBQ_GRID_SLOW
        if (NSTEP_DS) then
#  endif

          if (k.eq.0) then ! Bottom Boundary conditions

            do j=Jstr,Jend
              do i=Istr-1,Iend
                dthetadiv_nbqdz_u(i,j,k2)=0. 
              enddo
            enddo

            do j=Jstr-1,Jend
              do i=Istr,Iend
                dthetadiv_nbqdz_v(i,j,k2)=0.
              enddo
            enddo 

          else

            if (k.eq.N) then ! Top Boundary conditions

              do j=Jstr-1,Jend
                do i=Istr-1,Iend
#  ifdef NBQ_GRID_SLOW
                  dthetadiv_nbqdz(i,j,k,1)= - thetadiv_nbq(i,j,k)
#  else
                  dthetadiv_nbqdz(i,j)    = - thetadiv_nbq(i,j,k)
#  endif
                enddo
              enddo

            else

              do j=Jstr-1,Jend
                do i=Istr-1,Iend
#  ifdef NBQ_GRID_SLOW
                  dthetadiv_nbqdz(i,j,k,1)=thetadiv_nbq(i,j,k+1)
     &                                   - thetadiv_nbq(i,j,k)
#  else
                  dthetadiv_nbqdz(i,j)    =thetadiv_nbq(i,j,k+1) 
     &                                   - thetadiv_nbq(i  ,j,k)
#  endif
                enddo
              enddo

            endif
  
            do j=Jstr,Jend
              do i=Istr,Iend
#  ifdef NBQ_GRID_SLOW
                dthetadiv_nbqdz_u(i,j,k2)=Hzw_half_nbq_inv_u(i,j,k)*(
     &                                    dthetadiv_nbqdz(i,j,k,1)
     &                                   +dthetadiv_nbqdz(i-1,j,k,1))
#  else
                dthetadiv_nbqdz_u(i,j,k2)=Hzw_half_nbq_inv_u(i,j,k)*(
     &                                    dthetadiv_nbqdz(i,j)
     &                                   +dthetadiv_nbqdz(i-1,j))
#  endif          
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend
#  ifdef NBQ_GRID_SLOW
                 dthetadiv_nbqdz_v(i,j,k2)=Hzw_half_nbq_inv_v(i,j,k)*(
     &                                     dthetadiv_nbqdz(i,j,k,1)
     &                                    +dthetadiv_nbqdz(i,j-1,k,1))
#  else
                 dthetadiv_nbqdz_v(i,j,k2)=Hzw_half_nbq_inv_v(i,j,k)*(
     &                                     dthetadiv_nbqdz(i,j)
     &                                    +dthetadiv_nbqdz(i,j-1))
#  endif
              enddo
            enddo    
  
          endif    

#  ifdef NBQ_GRID_SLOW
        endif ! NSTEP_DS
#  endif

        if (k.gt.0) then
!
!-----------------------------------------------------------------------
!  Fast-mode U-momentum: qdmu_nbq
!-----------------------------------------------------------------------
!
          do j=Jstr,Jend
            do i=IstrU,Iend
              if (k.gt.1.and.k.lt.N) then 
#  ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,1)=(z_r(i  ,j,k)
     &                                     -z_r(i-1,j,k))  
     &              *(dthetadiv_nbqdz_u(i,j,k2)+
     &                dthetadiv_nbqdz_u(i,j,k1)) ! dZdx * (d(delta p)dz)_u
                endif 
                dum_s=dthetadiv_nbqdz(i,j,k,1)
#  else
                dum_s=(z_r(i,j,k)-z_r(i-1,j,k))                   
     &              *(dthetadiv_nbqdz_u(i,j,k2)+
     &                dthetadiv_nbqdz_u(i,j,k1)) ! dZdx * (d(delta p)dz)_u
#  endif
              elseif (k.gt.1) then
#  ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,1)=(z_r(i ,j,k)
     &                                     -z_r(i-1,j,k))       
     &                 *dthetadiv_nbqdz_u(i,j,k1) ! dZdx * (d(delta p)dz)_u
     &                 +(z_w(i,j,N)-z_w(i-1,j,N))
     &                 *dthetadiv_nbqdz_u(i,j,k2)
                endif
                dum_s=dthetadiv_nbqdz(i,j,k,1)
#  else
                dum_s=(z_r(i,j,k)-z_r(i-1,j,k))                      
     &                *dthetadiv_nbqdz_u(i,j,k1) ! dZdx * (d(delta p)dz)_u
     &                +(z_w(i,j,N)-z_w(i-1,j,N))   
     &                *dthetadiv_nbqdz_u(i,j,k2)
#  endif
              else
#  ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,1)=(z_r(i ,j,k)
     &                                     -z_r(i-1,j,k))       
     &             *2.*dthetadiv_nbqdz_u(i,j,k2) ! dZdx * (d(delta p)dz)_u
                endif
                dum_s=dthetadiv_nbqdz(i,j,k,1)
#  else
                dum_s=(z_r(i,j,k)-z_r(i-1,j,k))                      
     &              *2.*dthetadiv_nbqdz_u(i,j,k2) ! dZdx * (d(delta p)dz)_u
#  endif
              endif
#  ifdef MASKING
              if (umask(i-1,j)*umask(i+1,j) .eq. 0.) then
                dum_s=dum_s-(thetadiv_nbq(i  ,j,k)-
     &                       thetadiv_nbq(i-1,j,k))  ! - d(delta p)dx
              else
#  endif
                dum_s=dum_s
     &                -(gammau  *thetadiv_nbq(i  ,j,k)+
     &                  gammau_2*thetadiv_nbq(i+1,j,k)-
     &                  gammau  *thetadiv_nbq(i-1,j,k)-
     &                  gammau_2*thetadiv_nbq(i-2,j,k)) ! - d(delta p)dx
#  ifdef MASKING
              endif
#  endif
              dum_s=dum_s*Hzu_half_qdmu(i,j,k)
#  ifdef UV_COR_NT
              dum_s=dum_s+ntcoru(i,j,k)
#  endif
#  if defined INTERNAL || defined BODYTIDE
              dum_s=dum_s+0.5*omega*U0*cos(omega*time)
     &                       *(Hz(i-1,j,k)+Hz(i,j,k))
#  endif
#  ifdef BSTRESS_FAST
              if (k.eq.1) dum_s=dum_s-bustr(i,j)
#  endif
              qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k) + dtnbq * (
     &                        dum_s + ru_int_nbq(i,j,k) )  
#  ifdef MASKING
              qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k)*umask(i,j)
#  endif
              DU_nbq(i,j)=DU_nbq(i,j)+qdmu_nbq(i,j,k)
              ru_nbq(i,j,k)=dum_s/work(i,j) 

#  if defined NBQ_NUDGING && defined NBQCLIMATOLOGY
              qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k)*(1.-NBQnudgcof(i,j))
     &                        +u(i,j,k,nrhs)*Hzu_half_qdmu(i,j,k)
     &                                           *NBQnudgcof(i,j)
#  endif
            enddo 
          enddo
!
!-----------------------------------------------------------------------
!  Fast-mode V-momentum: qdmv_nbq
!-----------------------------------------------------------------------
!
          do j=JstrV,Jend
            do i=Istr,Iend
              if (k.gt.1.and.k.lt.N) then 
#  ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,2)=(z_r(i,j  ,k)
     &                                     -z_r(i,j-1,k)) 
     &                *(dthetadiv_nbqdz_v(i,j,k2)+
     &                  dthetadiv_nbqdz_v(i,j,k1)) ! dZdy * (d(delta p)dz)_v
                endif
                dum_s=dthetadiv_nbqdz(i,j,k,2)
#  else
                dum_s=(z_r(i,j,k)-z_r(i,j-1,k)) 
     &                *(dthetadiv_nbqdz_v(i,j,k2)+
     &                  dthetadiv_nbqdz_v(i,j,k1)) ! dZdy * (d(delta p)dz)_v
#  endif
              elseif (k.gt.1) then
#  ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,2)=(z_r(i,j  ,k)
     &                                     -z_r(i,j-1,k)) 
     &                *dthetadiv_nbqdz_v(i,j,k1) ! dZdy * (d(delta p)dz)_v
     &                +(z_w(i,j,N)-z_w(i,j-1,N))        
     &                *dthetadiv_nbqdz_v(i,j,k2)
                endif
                dum_s=dthetadiv_nbqdz(i,j,k,2)
#  else
                dum_s=(z_r(i,j,k)-z_r(i,j-1,k))            
     &                *dthetadiv_nbqdz_v(i,j,k1) ! dZdy * (d(delta p)dz)_v
     &                +(z_w(i,j,N)-z_w(i,j-1,N))      
     &                *dthetadiv_nbqdz_v(i,j,k2)
#  endif
              else
#  ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,2)=(z_r(i,j  ,k)
     &                                     -z_r(i,j-1,k)) 
     &                *2.*dthetadiv_nbqdz_v(i,j,k2) ! dZdy * (d(delta p)dz)_v
                endif
                dum_s=dthetadiv_nbqdz(i,j,k,2) 
#  else
                dum_s=(z_r(i,j,k)-z_r(i,j-1,k)) 
     &                *2.*dthetadiv_nbqdz_v(i,j,k2) ! dZdy * (d(delta p)dz)_v
#  endif
              endif
                
#  ifdef MASKING
              if (vmask(i,j-1)*vmask(i,j+1) .eq. 0.) then
              dum_s=dum_s
     &             -(thetadiv_nbq(i,j  ,k)-
     &               thetadiv_nbq(i,j-1,k)) ! - d(delta p)dy
              else
#  endif
              dum_s=dum_s
     &             -(gammau  *thetadiv_nbq(i,j  ,k)+
     &               gammau_2*thetadiv_nbq(i,j+1,k)-
     &               gammau  *thetadiv_nbq(i,j-1,k)-
     &               gammau_2*thetadiv_nbq(i,j-2,k)) ! - d(delta p)dy
#  ifdef MASKING
              endif
#  endif
              dum_s=dum_s*Hzv_half_qdmv(i,j,k)
#  ifdef UV_COR_NT
              dum_s=dum_s+ntcorv(i,j,k)
#  endif
#  if defined INTERNAL || defined BODYTIDE
              dum_s=dum_s+0.25*(f(i,j)+f(i,j-1))*
     &                       U0*sin(omega*time)*
     &                       (Hz(i,j,k)+Hz(i,j-1,k))
#  endif
#  ifdef BSTRESS_FAST
              if (k.eq.1) dum_s=dum_s-bvstr(i,j)
#  endif
              qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k) + dtnbq * (
     &                        dum_s + rv_int_nbq(i,j,k) )
#  ifdef MASKING
              qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k)*vmask(i,j)
#  endif
              DV_nbq(i,j)=DV_nbq(i,j)+qdmv_nbq(i,j,k)
              rv_nbq(i,j,k)=dum_s/work(i,j)  

#  if defined NBQ_NUDGING && defined NBQCLIMATOLOGY
              qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k)*(1.-NBQnudgcof(i,j))
     &                        +v(i,j,k,nrhs)*Hzv_half_qdmv(i,j,k)
     &                                           *NBQnudgcof(i,j)
#  endif
            enddo
          enddo
        endif  !<-- k>0
      enddo   !<-- k=0,N   

# else /* NBQ */

!
!********************************************************************
!
!  Advance fast 3D u,v for the hydrostatic case *****
!
!********************************************************************
!
#  ifdef BSTRESS_FAST
      do j=Jstr,Jend
        do i=IstrU,Iend
          rubar_nbq(i,j)=rubar_nbq(i,j)-bustr(i,j)
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend      
          rvbar_nbq(i,j)=rvbar_nbq(i,j)-bvstr(i,j)
        enddo
      enddo
#  endif
      if ((mod(iif-1,inc_faststep) .eq. inc_faststep-1) .OR.
     &    (LAST_FAST_STEP)) then
        do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend
              dum_s=0.
#  ifdef BSTRESS_FAST
              if (k.eq.1) dum_s=dum_s-bustr(i,j)
#  endif
              qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k)+dtfast * (
     &                      nb_faststep*(dum_s + ru_int_nbq(i,j,k))
     &                     +rubar_sum(i,j)*(Hz(i-1,j,k)+Hz(i,j,k)))    
#  ifdef MASKING
     &                                        *umask(i,j)
#  endif
              DU_nbq(i,j)=DU_nbq(i,j)+qdmu_nbq(i,j,k)
              if (LAST_FAST_STEP) ru_nbq(i,j,k)=dum_s/work(i,j) 
            enddo 
          enddo
          do j=JstrV,Jend
            do i=Istr,Iend
              dum_s=0.
#  ifdef BSTRESS_FAST
              if (k.eq.1) dum_s=dum_s-bvstr(i,j)
#  endif
              qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k)+dtfast * (
     &                       nb_faststep*(dum_s + rv_int_nbq(i,j,k))
     &                      +rvbar_sum(i,j)*(Hz(i,j-1,k)+Hz(i,j,k)))
#  ifdef MASKING
     &                                        *vmask(i,j)
#  endif
              DV_nbq(i,j)=DV_nbq(i,j)+qdmv_nbq(i,j,k)
              if (LAST_FAST_STEP) rv_nbq(i,j,k)=dum_s/work(i,j)
            enddo
          enddo
        enddo
      endif

# endif /* NBQ */

!
!-----------------------------------------------------------------------
! Apply point sources for river runoff simulations
!-----------------------------------------------------------------------
!
# ifdef PSOURCE
      do is=1,Nsrc 
#  ifdef MPI
        i=Isrc_mpi(is,mynode)
        j=Jsrc_mpi(is,mynode)
#  else
        i=Isrc(is)
        j=Jsrc(is)
#  endif
        if (IstrR.le.i .and. i.le.IendR .and.
     &      JstrR.le.j .and. j.le.JendR) then
          if (Dsrc(is).eq.0) then
            DU_nbq(i,j)=0.
#  ifdef WET_DRY
            if (Qbar(is).gt.0) then
              umask_wet(i,j)= + 1
            else
              umask_wet(i,j)= - 1
            endif
#  endif
            do k=1,N
              qdmu_nbq(i,j,k)=Qsrc(is,k)*pn_u(i,j)
              DU_nbq(i,j)=DU_nbq(i,j)+qdmu_nbq(i,j,k)
            enddo
          else
            DV_nbq(i,j)=0.
#  ifdef WET_DRY
            if (Qbar(is).gt.0) then
              vmask_wet(i,j)= + 1
            else
              vmask_wet(i,j)= - 1
            endif
#  endif
            do k=1,N
              qdmv_nbq(i,j,k)=Qsrc(is,k)*pm_v(i,j)
              DV_nbq(i,j)=DV_nbq(i,j)+qdmv_nbq(i,j,k)
            enddo
          endif
        endif
      enddo
# endif
!
!-----------------------------------------------------------------------
!  U & V momentum wet mask
!-----------------------------------------------------------------------
!
# ifdef WET_DRY
      do j=Jstr,Jend
        do i=IstrU,Iend
          cff1_WD=ABS(ABS(umask_wet(i,j))-1.)
          cff2_WD=0.5+SIGN(0.5,DU_nbq(i,j))*umask_wet(i,j)
          umask_wet(i,j)=0.5*umask_wet(i,j)*cff1_WD
     &                         +cff2_WD*(1.-cff1_WD)
          DU_nbq(i,j)=DU_nbq(i,j)*umask_wet(i,j)
#  ifdef MRL_WCI
          ust2d(i,j)=ust2d(i,j)*umask_wet(i,j)
#  endif
          do k=1,N
            qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k)*umask_wet(i,j)
          enddo
        enddo 
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          cff1_WD=ABS(ABS(vmask_wet(i,j))-1.)
          cff2_WD=0.5+SIGN(0.5,DV_nbq(i,j))*vmask_wet(i,j)
          vmask_wet(i,j)=0.5*vmask_wet(i,j)*cff1_WD
     &                         +cff2_WD*(1.-cff1_WD)
          DV_nbq(i,j)=DV_nbq(i,j)*vmask_wet(i,j)
#  ifdef MRL_WCI
          vst2d(i,j)=vst2d(i,j)*vmask_wet(i,j)
#  endif
          do k=1,N
            qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k)*vmask_wet(i,j)
          enddo
        enddo 
      enddo
# endif

# ifdef RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
      call check_tab3d(qdmu_nbq,'qdmu_nbqint','uint')
      call check_tab3d(qdmv_nbq,'qdmv_nbqint','vint')
C$OMP END MASTER
# endif
!
!-----------------------------------------------------------------------
!  U & V momentum open boundary conditions
!-----------------------------------------------------------------------
!
      do j=Jstr,Jend
        do i=IstrU,Iend
          ubar(i,j,knew)=urhs(i,j)
        enddo 
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          vbar(i,j,knew)=vrhs(i,j)
        enddo 
      enddo
# ifdef NBQ
      M2bc_nbq_flag=.true. ! apply boundary wet/dry conditions
                           ! and compute DU_nbq
      call u2dbc_tile   (Istr,Iend,Jstr,Jend, work)
      call v2dbc_tile   (Istr,Iend,Jstr,Jend, work)
# endif
      call unbq_bc_tile (Istr,Iend,Jstr,Jend, work)
      call vnbq_bc_tile (Istr,Iend,Jstr,Jend, work)
!
!-----------------------------------------------------------------------
! Exchange periodic boundaries and computational margins
!-----------------------------------------------------------------------
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI  
#  ifndef NBQ
      if ((mod(iif-1,inc_faststep) .eq. inc_faststep-1) .OR.
     &     (LAST_FAST_STEP)) then
#  endif
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        qdmu_nbq(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        qdmv_nbq(START_2D_ARRAY,1))
#  ifndef NBQ
      endif
#  endif
!      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
!     &                        DU_nbq(START_2D_ARRAY))
!      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
!     &                        DV_nbq(START_2D_ARRAY))
# endif
!$$$# ifdef RVTK_DEBUG_ADVANCED
!$$$C$OMP BARRIER
!$$$C$OMP MASTER
!$$$      call check_tab3d(qdmu_nbq,'qdmu_nbq','u')
!$$$      call check_tab3d(qdmv_nbq,'qdmv_nbq','v')
!$$$C$OMP END MASTER
!$$$# endif  


# ifdef NBQ
!
!-----------------------------------------------------------------------
!      Fast mode W-Momentum equation: qdmw_nbq
!                 ==> EXPLICIT scheme
!-----------------------------------------------------------------------
!
#  ifndef NBQ_IMP
      do j=Jstr,Jend
        do k=1,N-1
          do i=Istr,Iend   
            dum_s = thetadiv_nbq(i,j,k) - thetadiv_nbq(i,j,k+1)   
#   ifdef UV_COR_NT
            dum_s = dum_s + ntcorw(i,j,k)
#   endif
            qdmw_nbq(i,j,k)=qdmw_nbq(i,j,k)   
     &                      + dtnbq * ( dum_s + rw_int_nbq(i,j,k) )
#   ifdef NBQ_GRAV
     &               -0.25*(rho_nbq(i,j,k  )*Hzr_half_nbq_inv(i,j,k  )
     &                     +rho_nbq(i,j,k+1)*Hzr_half_nbq_inv(i,j,k+1))
     &                    *(Hzr(i,j,k)+Hzr(i,j,k+1))
     &                    *g*dtnbq
#   endif
#   ifdef MASKING
          qdmw_nbq(i,j,k)=qdmw_nbq(i,j,k) * rmask(i,j)
#   endif
          enddo             
        enddo

        k=N
        do i=Istr,Iend
          dum_s =   thetadiv_nbq(i,j,N)
#   ifdef UV_COR_NT
          dum_s = dum_s + ntcorw(i,j,k)
#   endif
          qdmw_nbq(i,j,N)=qdmw_nbq(i,j,N)   
     &                   + dtnbq * ( dum_s + rw_int_nbq(i,j,N) )
#   ifdef NBQ_GRAV
     &                   -rho_nbq(i,j,N)*0.5
     &                   *g*dtnbq
#   endif
#   ifdef MASKING
          qdmw_nbq(i,j,k)=qdmw_nbq(i,j,k) * rmask(i,j)
#   endif
        enddo
      enddo !<-- j loop
!
!-----------------------------------------------------------------------
! Vertical momentum open boundary conditions
!-----------------------------------------------------------------------
!
      call wnbq_bc_tile (Istr,Iend,Jstr,Jend, work)

#   ifdef UV_COR_NT
#    if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        qdmw_nbq(START_2D_ARRAY,0))
#    endif
#   endif

#  endif  /* ! NBQ_IMP */
!
!
!=====================================================================
!      Fast-mode conservation of mass
!  
!  ... and integration of W-momentum with IMPLICIT scheme
!
!  From now on, thetadiv_nbq array is used for divergence (not theta):
!
!=====================================================================
!
      if (IstrU.gt.Iend) then
        do j=Jstr,Jend
          do i=Istr,Iend+1
            FX(i,j)=0.
          enddo
        enddo
      endif
      if (JstrV.gt.Jend) then
        do j=Jstr,Jend+1
          do i=Istr,Iend 
            FY(i,j)=0.
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!  X-component dZdx*qdmu for horizontal divergence
!-----------------------------------------------------------------------
!
      k2 = 1
      do k=0,N
        k1=k2
        k2=3-k1

#  ifdef NBQ_GRID_SLOW
        if (NSTEP_DS) then
#  endif
          if (k.lt.N) then
            kp1 = k + 1
            do j=Jstr,Jend
              do i=Istr,Iend+1
                dZdxq_u(i,j,k2)=(z_r(i,j,kp1)-z_r(i-1,j,kp1)) 
     &                           *qdmu_nbq(i,j,kp1)  ! (dZdx * (rho u))_u
              enddo
            enddo
          endif

          if (k.eq.0) then  ! Bottom boundary conditions

#  ifdef NBQ_FREESLIP
            do j=Jstr,Jend
              do i=Istr,Iend+1 
#   ifdef NBQ_GRID_SLOW
                dZdxq_w(i,j,k ) = (z_w(i,j,0)-z_w(i-1,j,0))
     &                             *qdmu_nbq(i,j,1)  
     &                             /(Hzr(i,j,1)+Hzr(i-1,j,1))
#   else
                dZdxq_w(i,j,k2)= (z_w(i,j,0)-z_w(i-1,j,0))
     &                            *qdmu_nbq(i,j,1)  
     &                            /(Hzr(i,j,1)+Hzr(i-1,j,1))
#   endif
              enddo
            enddo

            do j=Jstr,Jend
              do i=Istr,Iend    
#   ifdef NBQ_GRID_SLOW
                qdmw_nbq(i,j,0)=0.5*(dZdxq_w(i  ,j,k)*pm_u(i  ,j)
     &                              +dZdxq_w(i+1,j,k)*pm_u(i+1,j) )
     &                                                 * Hzr(i,j,1)
#   else
                qdmw_nbq(i,j,0)=0.5*(dZdxq_w(i  ,j,k2)*pm_u(i  ,j)
     &                              +dZdxq_w(i+1,j,k2)*pm_u(i+1,j) )
     &                                                  * Hzr(i,j,1)
#   endif
#   ifdef MASKING
                qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)*rmask(i,j)
#   endif 
              enddo
            enddo 

#  else /* NBQ_FREESLIP */

#   ifndef NBQ_GRID_SLOW
            do j=Jstr,Jend
              do i=Istr,Iend+1
                dZdxq_w(i,j,k2)=0.  
                qdmw_nbq(i,j,0)=0.     
              enddo
            enddo 
#   endif
 
#  endif /* NBQ_FREESLIP */

          elseif (k .eq. N) then ! Top boundary conditions
           
            do j=Jstr,Jend
              do i=Istr,Iend+1
#  ifdef NBQ_GRID_SLOW
                dZdxq_w(i,j,k )= (z_w(i,j,N)-z_w(i-1,j,N))
     &                            *qdmu_nbq(i,j,N)
     &                            /(Hzr(i,j,N)+Hzr(i-1,j,N))
#  else
                dZdxq_w(i,j,k2)= (z_w(i,j,N)-z_w(i-1,j,N))
     &                            *qdmu_nbq(i,j,N)
     &                            /(Hzr(i,j,N)+Hzr(i-1,j,N))
#  endif
              enddo
            enddo  
 
          else ! k<>0 & k<>N   ! Inner domain   

            do j=Jstr,Jend
              do i=Istr,Iend+1
#  ifdef NBQ_GRID_SLOW
                dZdxq_w(i,j,k )=Hzw_half_nbq_inv_u(i,j,k)
     &                                 *(dZdxq_u(i,j,k1)+
     &                                   dZdxq_u(i,j,k2))
#  else
                dZdxq_w(i,j,k2)=Hzw_half_nbq_inv_u(i,j,k)
     &                                 *(dZdxq_u(i,j,k1)+
     &                                   dZdxq_u(i,j,k2))
#  endif
              enddo 
            enddo

          endif ! k<>0 , k<>N ,  Inner domain 

#  ifdef NBQ_GRID_SLOW

        else   ! NSTEP_DS: Update d./ds terms

          if (k.eq.0) then  ! Bottom boundary conditions

#   ifdef NBQ_FREESLIP
            do j=Jstr,Jend
              do i=Istr,Iend+1 
                dZdxq_w(i,j,k)=(z_w(i,j,0)-z_w(i-1,j,0))
     &                          *qdmu_nbq(i,j,1)  
     &                          /(Hzr(i,j,1)+Hzr(i-1,j,1))
              enddo
            enddo

            do j=Jstr,Jend
              do i=Istr,Iend    
                qdmw_nbq(i,j,0)=0.5*(dZdxq_w(i  ,j,k)*pm_u(i  ,j) 
     &                              +dZdxq_w(i+1,j,k)*pm_u(i+1,j) ) 
     &                                                  *Hzr(i,j,1)    
#    ifdef MASKING
                qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)*rmask(i,j)
#    endif 
              enddo
            enddo 

#   endif /* NBQ_FREESLIP */

          endif ! k.eq.0
        endif ! NSTEP_DS: Update d./ds terms

#  else  /* NBQ_GRID_SLOW */

        if (k.eq.0) then ! Bottom boundary conditions

#   ifdef NBQ_FREESLIP
          do j=Jstr,Jend
            do i=Istr,Iend+1 
              dZdxq_w(i,j,k2)=(z_w(i,j,0)-z_w(i-1,j,0))
     &                         *qdmu_nbq(i,j,1)  
     &                         /(Hzr(i,j,1)+Hzr(i-1,j,1))
            enddo
          enddo

          do j=Jstr,Jend
            do i=Istr,Iend    
              qdmw_nbq(i,j,0)=0.5*(dZdxq_w(i  ,j,k2)*pm_u(i  ,j) 
     &                            +dZdxq_w(i+1,j,k2)*pm_u(i+1,j) ) 
     &                                                 *Hzr(i,j,1)    
#    ifdef MASKING
              qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)*rmask(i,j)
#    endif 
            enddo
          enddo 

#   endif /* NBQ_FREESLIP */ 

        endif ! k.eq.0

#  endif  /* NBQ_GRID_SLOW */

        if (k.gt.0) then
          if (IstrU.le.Iend) then
            do j=Jstr,Jend
              do i=Istr,Iend+1
#  ifdef NBQ_GRID_SLOW
                FX(i,j)=-pm_u(i,j)*(dZdxq_w(i,j,k )-dZdxq_w(i,j,k-1))
#  else
                FX(i,j)=-pm_u(i,j)*(dZdxq_w(i,j,k2)-dZdxq_w(i,j,k1))
#  endif
#  ifdef MASKING
                FX(i,j)=FX(i,j)*umask(i,j)
#  endif                
              enddo
            enddo

            do j=Jstr,Jend
              do i=Istr,Iend
                thetadiv_nbq(i,j,k)=FX(i,j)+FX(i+1,j)
              enddo
            enddo

          else ! IstrU.gt.Iend

            do j=Jstr,Jend
              do i=Istr,Iend
                thetadiv_nbq(i,j,k)=0.
              enddo
            enddo

          endif 
        endif
      enddo ! k=0,N
!
!-----------------------------------------------------------------------
!  Y-component dZdy*qdmv for horizontal divergence
!-----------------------------------------------------------------------
!
      k2 = 1
      do k=0,N  !<-- k loop
        k1=k2
        k2=3-k1

#  ifdef NBQ_GRID_SLOW
        if (NSTEP_DS) then
#  endif
          if (k.lt.N) then
            kp1 = k + 1
            do j=Jstr,Jend+1
              do i=Istr,Iend
                dZdyq_v(i,j,k2)=(z_r(i,j,kp1)-z_r(i,j-1,kp1)) 
     &                           *qdmv_nbq(i,j,kp1)  ! (dZdy * (rho v))_v
              enddo
            enddo
          endif

          if (k.eq.0) then  ! Bottom boundary conditions

#  ifdef NBQ_FREESLIP
            do j=Jstr,Jend+1
              do i=Istr,Iend
#   ifdef NBQ_GRID_SLOW
                 dZdyq_w(i,j,k )= (z_w(i,j,0)-z_w(i,j-1,0))
     &                             *qdmv_nbq(i,j,1) 
     &                             /(Hzr(i,j,1)+Hzr(i,j-1,1))
#   else
                 dZdyq_w(i,j,k2)= (z_w(i,j,0)-z_w(i,j-1,0))
     &                             *qdmv_nbq(i,j,1) 
     &                             /(Hzr(i,j,1)+Hzr(i,j-1,1))
#   endif 
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend   
#   ifdef NBQ_GRID_SLOW
                qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)
     &                          +0.5*(dZdyq_w(i,j  ,k)*pn_v(i,j  )  
     &                               +dZdyq_w(i,j+1,k)*pn_v(i,j+1) )
     &                                                  * Hzr(i,j,1)
#   else
                qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)
     &                          +0.5*(dZdyq_w(i,j  ,k2)*pn_v(i,j  )  
     &                               +dZdyq_w(i,j+1,k2)*pn_v(i,j+1) )
     &                                                   * Hzr(i,j,1)
#   endif 
#   ifdef MASKING
                qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)*rmask(i,j)
#   endif 
              enddo
            enddo

#  else  /* NBQ_FREESLIP */

#   ifndef NBQ_GRID_SLOW   
            do j=Jstr,Jend +1
              do i=Istr,Iend  
                dZdyq_w(i,j,k2)=0.  
                qdmw_nbq(i,j,0)=0.
              enddo
            enddo 
#   endif
#  endif /* NBQ_FREESLIP */

          elseif (k .eq. N) then ! Top boundary conditions

            do j=Jstr,Jend+1
              do i=Istr,Iend
#  ifdef NBQ_GRID_SLOW
                dZdyq_w(i,j,k )= (z_w(i,j,N)-z_w(i,j-1,N))
     &                           *qdmv_nbq(i,j,N)
     &                           /(Hzr(i,j,N)+Hzr(i,j-1,N))
#  else
                dZdyq_w(i,j,k2)= (z_w(i,j,N)-z_w(i,j-1,N))
     &                            *qdmv_nbq(i,j,N)
     &                            /(Hzr(i,j,N)+Hzr(i,j-1,N))
#  endif 
              enddo
            enddo

          else

            do j=Jstr,Jend+1
              do i=Istr,Iend
#  ifdef NBQ_GRID_SLOW
                dZdyq_w(i,j,k )=Hzw_half_nbq_inv_v(i,j,k)
     &                           *(dZdyq_v(i,j,k1)+
     &                             dZdyq_v(i,j,k2)) ! (dZdy * (rho v))_uw/Hzw_v
#  else
                dZdyq_w(i,j,k2)=Hzw_half_nbq_inv_v(i,j,k)
     &                           *(dZdyq_v(i,j,k1)+
     &                             dZdyq_v(i,j,k2)) ! (dZdy * (rho v))_uw/Hzw_v
#  endif 
              enddo 
            enddo

          endif

#  ifdef NBQ_GRID_SLOW

        else  ! NSTEP_DS

          if (k.eq.0) then ! Bottom boundary conditions

#   ifdef NBQ_FREESLIP
            do j=Jstr,Jend+1
              do i=Istr,Iend
                dZdyq_w(i,j,k)= (z_w(i,j,0)-z_w(i,j-1,0))
     &                           *qdmv_nbq(i,j,1) 
     &                           /(Hzr(i,j,1)+Hzr(i,j-1,1))
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend   
                qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)
     &                         +0.5*(dZdyq_w(i,j  ,k)*pn_v(i,j  )  
     &                              +dZdyq_w(i,j+1,k)*pn_v(i,j+1) )
     &                                                 * Hzr(i,j,1) 
#    ifdef MASKING
                qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)*rmask(i,j)
#    endif 
              enddo
            enddo
#   endif
         endif
        endif ! NSTEP_DS

#  else /* NBQ_GRID_SLOW */

        if (k.eq.0) then  ! Bottom boundary conditions

#   ifdef NBQ_FREESLIP
          do j=Jstr,Jend+1
            do i=Istr,Iend
              dZdyq_w(i,j,k2)= (z_w(i,j,0)-z_w(i,j-1,0))
     &                          *qdmv_nbq(i,j,1) 
     &                          /(Hzr(i,j,1)+Hzr(i,j-1,1))
            enddo
          enddo
          do j=Jstr,Jend
            do i=Istr,Iend   
               qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)
     &                         +0.5*(dZdyq_w(i,j  ,k2)*pn_v(i,j  )  
     &                              +dZdyq_w(i,j+1,k2)*pn_v(i,j+1) )
     &                                                  * Hzr(i,j,1) 
#    ifdef MASKING
               qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)*rmask(i,j)
#    endif 
            enddo
          enddo
#   endif  /* NBQ_FREESLIP */
        endif

#  endif /* NBQ_GRID_SLOW */

        if (k.gt.0) then
          if (JstrV.le.Jend) then
            do j=Jstr,Jend+1
              do i=Istr,Iend 
#  ifdef NBQ_GRID_SLOW
                FY(i,j)=-pn_v(i,j)*(dZdyq_w(i,j,k)-dZdyq_w(i,j,k-1))
#  else
                FY(i,j)=-pn_v(i,j)*(dZdyq_w(i,j,k2)-dZdyq_w(i,j,k1))
#  endif 
#  ifdef MASKING
                FY(i,j)=FY(i,j)*vmask(i,j)
#  endif                 
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend
                thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                               +FY(i,j)+FY(i,j+1)
              enddo
            enddo
          endif
        endif

      enddo  !<-- k=0,N
!
!-----------------------------------------------------------------------
!  Compute total horizontal Divergence divH(qdmH(m+1))
!-----------------------------------------------------------------------
!
      do k=1,N !<-- k loop

        if (IstrU.le.Iend) then
          do j=Jstr,Jend
            do i=Istr,Iend+1
              FX(i,j)=on_u(i,j)*qdmu_nbq(i,j,k)
            enddo
          enddo
        endif
        if (JstrV.le.Jend) then
          do j=Jstr,Jend+1  
            do i=Istr,Iend
              FY(i,j)=om_v(i,j)*qdmv_nbq(i,j,k)
            enddo
          enddo
        endif
        if (IstrU.gt.Iend) then
          do j=Jstr,Jend
            do i=Istr,Iend  
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                           +pm(i,j)*pn(i,j)*(FY(i,j+1)-FY(i,j))  
#  ifdef MASKING
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)*rmask(i,j)
#  endif                              
            enddo
          enddo
        elseif (JstrV.gt.Jend) then
          do j=Jstr,Jend  
            do i=Istr,Iend   
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                            +pm(i,j)*pn(i,j)*(FX(i+1,j)-FX(i,j))  
#  ifdef MASKING
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)*rmask(i,j)
#  endif                              
            enddo
          enddo
        else
          do j=Jstr,Jend
            do i=Istr,Iend
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                            +pm(i,j)*pn(i,j)*(FX(i+1,j)-FX(i,j)+
     &                                              FY(i,j+1)-FY(i,j))  
#  ifdef MASKING
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)*rmask(i,j)
#  endif                              
            enddo
          enddo
        endif  

      enddo ! <-- k=1,N
!
!-----------------------------------------------------------------------
!  thetadiv2_nbq: complet time-corrective term  (dh/dt included) 
!  thetadiv3_nbq: reduced time-corrective term  (no dh/dt)
!-----------------------------------------------------------------------
!
#  ifdef NBQ_GRID_SLOW
       if (FIRST_FAST_STEP) then
#  endif

        do j=JstrR,JendR
          do k=0,N
            do i=IstrR,IendR
              zw_nbq(i,j,k,knew)=z_w(i,j,k)
            enddo
          enddo
        enddo
        do j=Jstr,Jend !<-- j loop

          do i=Istr,Iend
            FC(i,0)=0.    ! Bottom boundary condition
            CF(i,0)=0.
          enddo

          do k=1,N-1
            do i=Istr,Iend
              FC(i,k)=   
     &          -(zw_nbq(i,j,k,knew)-zw_nbq(i,j,k,kstp))/dtgrid_nbq
     &           *0.5*( (1.+rho_grd(i,j,k  ))
     &                  +rho_nbq(i,j,k  )*Hzr_half_nbq_inv(i,j,k  )  
     &                  +(1.+rho_grd(i,j,k+1))
     &                  +rho_nbq(i,j,k+1)*Hzr_half_nbq_inv(i,j,k+1))
              CF(i,k)=   
     &          -(zw_nbq(i,j,k,knew)-zw_nbq(i,j,k,kstp))/dtgrid_nbq 
     &          *0.5*(    rho_grd(i,j,k  ) 
     &                   +rho_nbq(i,j,k  )*Hzr_half_nbq_inv(i,j,k  ) 
     &                   +rho_grd(i,j,k+1)
     &                   +rho_nbq(i,j,k+1)*Hzr_half_nbq_inv(i,j,k+1) )

              thetadiv3_nbq(i,j,k)=CF(i,k)-CF(i,k-1) 
              thetadiv2_nbq(i,j,k)=FC(i,k)-FC(i,k-1) 
            enddo
          enddo
    
           do i=Istr,Iend
             FC(i,N)=   
     &        -(zw_nbq(i,j,N,knew)-zw_nbq(i,j,N,kstp))/dtgrid_nbq
     &         *( 1.+rho_grd(i,j,N)
     &          )

            CF(i,N)=  
     &        -(zw_nbq(i,j,N,knew)-zw_nbq(i,j,N,kstp))/dtgrid_nbq
     &         *( rho_grd(i,j,N)
     &          )
            thetadiv3_nbq(i,j,N)=CF(i,N)-CF(i,N-1) 
            thetadiv2_nbq(i,j,N)=FC(i,N)-FC(i,N-1) 
          enddo
        enddo !<-- j loop

#  ifdef NBQ_GRID_SLOW
       endif !<-- FIRST_FAST_STEP
#  endif
!
!-------------------------------------------------------------------
!  Solve implicit W-momentum equation
!-------------------------------------------------------------------
!
#  ifdef NBQ_IMP
      do k=0,N
        do j=Jstr,Jend
          do i=Istr,Iend
            qdmw_nbq_old(i,j,k)=qdmw_nbq(i,j,k)
          enddo
        enddo
      enddo
  
      do j=Jstr,Jend ! <-- j loop

        do k=1,N
          do i=Istr,Iend
            FC(i,k)= soundspeed2_nbq(i,j)*rho_nbq(i,j,k)
     &              -(thetaimp_nbq*soundspeed2_nbq(i,j)*dtnbq+visc2_nbq)
     &                       *(thetadiv_nbq(i,j,k)+thetadiv3_nbq(i,j,k))
    
#   ifdef NBQ_THETAIMP 
            FC(i,k)=FC(i,k)
     &       -thetaimp_nbq*(1.-thetaimp_nbq)*soundspeed2_nbq(i,j)*dtnbq*
     &                   ( Hzw_half_nbq_inv(i,j,k  )*qdmw_nbq(i,j,k)
     &                    -Hzw_half_nbq_inv(i,j,k-1)*qdmw_nbq(i,j,k-1) )
#   endif

            FC(i,k)=FC(i,k)*Hzr_half_nbq_inv(i,j,k) 
          enddo
        enddo
     
!.........Inner layers

        do k=1,N-1
          do i=Istr,Iend   
            dum_s = FC(i,k) - FC(i,k+1)
#   ifdef UV_COR_NT
            dum_s = dum_s + ntcorw(i,j,k)
#   endif
            qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k)   
     &                        + dtnbq * ( dum_s + rw_int_nbq(i,j,k) )
#   ifdef MASKING
            qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k) * rmask(i,j)
#   endif               
          enddo             
        enddo

!.........Surface BC

        k=N
        do i=Istr,Iend
          dum_s =   FC(i,k)
#   ifdef UV_COR_NT
          dum_s = dum_s + ntcorw(i,j,k)
#   endif
          qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k)   
     &                      + dtnbq * ( dum_s + rw_int_nbq(i,j,k) )
#   ifdef MASKING
          qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k) * rmask(i,j)
#   endif              
        enddo   

      enddo !<-- j loop


!---------------------
! Gaussian Elimination
!---------------------

!.......Compute coef.

      do j=Jstr,Jend  !<-- j loop

!..........Bottom BC
        k=1
        do i=Istr,Iend
          cff1=1./(dtnbq*(thetaimp_nbq**2*soundspeed2_nbq(i,j)
     &            *dtnbq+visc2_nbq)) 
          cff=(cff1+Hzw_half_nbq_inv(i,j,1)*(Hzr_half_nbq_inv(i,j,1)
     &             +Hzr_half_nbq_inv(i,j,2)))
          CF(i,1)=(-Hzw_half_nbq_inv(i,j,2)*Hzr_half_nbq_inv(i,j,2))
     &                                                          /cff
          DC(i,1)=qdmw_nbq(i,j,1)*cff1/cff   
     &           +qdmw_nbq(i,j,0)/cff*Hzw_half_nbq_inv(i,j,0)
     &                               *Hzr_half_nbq_inv(i,j,1)
        enddo

!..........Inner layers
        do k=2,N-1
          do i=Istr,Iend
            cff1=1./(dtnbq*(thetaimp_nbq**2*soundspeed2_nbq(i,j)
     &              *dtnbq+visc2_nbq)) 
            cff=(cff1+
     &               Hzw_half_nbq_inv(i,j,k)*(Hzr_half_nbq_inv(i,j,k)
     &              +Hzr_half_nbq_inv(i,j,k+1))   
     &              +Hzw_half_nbq_inv(i,j,k-1)*Hzr_half_nbq_inv(i,j,k)
     &              *CF(i,k-1))
            CF(i,k)=(-Hzw_half_nbq_inv(i,j,k+1)
     &               *Hzr_half_nbq_inv(i,j,k+1))/cff
            DC(i,k)=(qdmw_nbq(i,j,k)*cff1+Hzw_half_nbq_inv(i,j,k-1)
     &               *Hzr_half_nbq_inv(i,j,k)*DC(i,k-1)) /cff           
          enddo            
        enddo

!..........Surface BC
        k=N
        do i=Istr,Iend
          cff1=1./(dtnbq*(thetaimp_nbq**2*soundspeed2_nbq(i,j)
     &            *dtnbq+visc2_nbq)) 
          cff=(cff1+Hzw_half_nbq_inv(i,j,N)*Hzr_half_nbq_inv(i,j,N) 
     &             +Hzw_half_nbq_inv(i,j,N-1)
     &             *Hzr_half_nbq_inv(i,j,N)*CF(i,N-1))  
          CF(i,N)=0. 
          DC(i,k)=(qdmw_nbq(i,j,N)*cff1+Hzw_half_nbq_inv(i,j,N-1)
     &             *Hzr_half_nbq_inv(i,j,N)*DC(i,N-1))/cff
        enddo 

!..........Solves tri-diag system
        do i=Istr,Iend
          qdmw_nbq(i,j,N)=DC(i,k)   
#   ifdef NBQ_GRAV
     &                    -rho_nbq(i,j,N)*0.5*g*dtnbq
#   endif     
        enddo

        do k=N-1,1,-1
          do i=Istr,Iend
            qdmw_nbq(i,j,k)=DC(i,k)-CF(i,k)*qdmw_nbq(i,j,k+1)
#   ifdef NBQ_GRAV
     &            -0.25*(rho_nbq(i,j,k)*Hzr_half_nbq_inv(i,j,k)
     &                 +rho_nbq(i,j,k+1)*Hzr_half_nbq_inv(i,j,k+1))
     &                *(Hzr(i,j,k)+Hzr(i,j,k+1))
     &                 *g*dtnbq
#   endif
          enddo            
        enddo                        

      enddo   !<-- j loop 
!
!-------------------------------------------------------------------
!  W-momentum open boundary conditions
!-------------------------------------------------------------------
!
#   ifdef OBC_NBQ
      call wnbq_bc_tile (Istr,Iend,Jstr,Jend, work)
#   endif

#   ifdef UV_COR_NT
#    if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        qdmw_nbq(START_2D_ARRAY,0))
#    endif
#   endif

#  endif /* NBQ_IMP */
!
!-------------------------------------------------------------------
!  Mass conservation equation
!
!      div(m+1) = div(m+1) + divZ(qdmw_nbq(m+1))
!-------------------------------------------------------------------
!		
!.......Computes fluxes:  
!
      do j=Jstr,Jend  !<-- j loop
            
#  ifdef NBQ_FREESLIP
        do i=Istr,Iend
          FC(i,0)=Hzw_half_nbq_inv(i,j,0) *
     &            (thetaimp_nbq*qdmw_nbq(i,j,0)
     &           +(1.-thetaimp_nbq)*qdmw_nbq_old(i,j,0)) ! Bottom BC
        enddo
#  else
        do i=Istr,Iend
          FC(i,0)=0.   ! Bottom BC
        enddo
#  endif
        do k=1,N-1
          do i=Istr,Iend
            FC(i,k)=Hzw_half_nbq_inv(i,j,k) *
     &              (thetaimp_nbq*qdmw_nbq(i,j,k)
     &              +(1.-thetaimp_nbq)*qdmw_nbq_old(i,j,k))
            thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                         +FC(i,k)-FC(i,k-1)  
          enddo
        enddo
        do i=Istr,Iend
          FC(i,N)=Hzw_half_nbq_inv(i,j,N) *
     &            (thetaimp_nbq*qdmw_nbq(i,j,N)
     &           +(1.-thetaimp_nbq)*qdmw_nbq_old(i,j,N))
          thetadiv_nbq(i,j,N)=thetadiv_nbq(i,j,N)
     &                        +FC(i,N)-FC(i,N-1)    
        enddo

      enddo  !<-- j loop
!                 
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        thetadiv_nbq(START_2D_ARRAY,1))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        thetadiv2_nbq(START_2D_ARRAY,1))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        thetadiv3_nbq(START_2D_ARRAY,1))
#  endif
!
!=====================================================================
!
!      Solve rho_nbq from mass conservation 
!
!      using a Forward-Backward scheme:
!          rho_nbq(m+1) = rho_nbq(m) - dtfast*DIV(m+1)
!
!=====================================================================
!
      do k=1,N
        do j=JstrV-2,Jend+1
          do i=IstrU-2,Iend+1
            rho_nbq(i,j,k) = rho_nbq(i,j,k)  
     &                       - dtfast*(thetadiv_nbq(i,j,k)+
     &                                 thetadiv3_nbq(i,j,k))
          enddo
        enddo
      enddo

!     call rnbq_bc_tile(Istr,Iend,Jstr,Jend, work) !!! leave commented
!
!# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
!      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
!     &                       rho_nbq(START_2D_ARRAY,1))
!# endif
!
!-------------------------------------------------------------------
!  Acoustic wave emission
!-------------------------------------------------------------------
!
#  ifdef ACOUSTIC
        time_nbq = time_nbq + dtfast
        do k=1,N
          do j=JstrV-2,Jend+1
            do i=IstrU-2,Iend+1
              dist_d=sqrt((xr(i,j)-xl/2.)**2+(0.*(yr(i,j)-el/2.))**2
     &                             +(abs(z_w(i,j,k))-hmax_exp/2.)**2)

              rho_nbq(i,j,k)=rho_nbq(i,j,k)+amp_exp
     &                         *sin(2*acos(-1.)*time_nbq/period_exp)
     &                                 *exp(-dist_d**2/for_a_exp**2)
     &                                                   *Hzr(i,j,k)
            enddo
          enddo
        enddo
#  endif
!
!-------------------------------------------------------------------
!  rhobar_nbq: depth-mean density (/rho0)
!-------------------------------------------------------------------
!		
#  ifdef NBQ_MASS
!
! Compute rhobar_nbq(m+1) used in zeta diagnostic from 
! depth-integrated continuity equation
!
      do j=Jstr-1,Jend+1
        do i=Istr-1,Iend+1
          rhobar_nbq(i,j,knew)=0.
        enddo  
      enddo
      do k=1,N
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            rhobar_nbq(i,j,knew)= rhobar_nbq(i,j,knew)
     &                             +rho_nbq(i,j,k)
     &                             +rho_grd(i,j,k)*Hzr(i,j,k)
          enddo  
        enddo  
      enddo
      do j=Jstr-1,Jend+1
        do i=Istr-1,Iend+1
          rhobar_nbq(i,j,knew) = 1.+(rhobar_nbq(i,j,knew)) 
     &                            / (z_w(i,j,N)-z_w(i,j,0))
        enddo
      enddo
!      
! LAURENT: the next exchange should not be needed
!
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend, 
     &                           rhobar_nbq(START_2D_ARRAY,knew))
#   endif
#   ifdef RVTK_DEBUG_ADVANCED
!      call check_tab2d(rhobar_nbq(:,:,knew),'rhobar_nbq','r')
#   endif    
#  endif /* NBQ_MASS */

# endif /* NBQ */
!
!
!***********************************************************************
!
!
!                     FAST-MODE POST-PROCESSING
!
!
!***********************************************************************
!
!=====================================================================
!  Update total mass of water volume Dnew 
!=====================================================================
!
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          Dnew(i,j)=(zeta(i,j,knew)+h(i,j))
# ifdef NBQ_MASS
     &               *rhobar_nbq(i,j,knew)
# endif
        enddo
      enddo
            
# if defined RVTK_DEBUG_ADVANCED && defined NBQ
C$OMP BARRIER
C$OMP MASTER
      call check_tab3d(rho_nbq,'rho_nbq','r')
C$OMP END MASTER
# endif    

!=====================================================================
!  Update W NBQ forcing for internal mode (rw_nbq ~ rw_nbq_avg1)
!  Note: here rw_nbq contains qdmw_nbq(m)
!=====================================================================
!
# ifdef NBQ
      if (LAST_FAST_STEP) then
        do k=0,N 
          do j=Jstr,Jend              
            do i=Istr,Iend
              rw_nbq(i,j,k)=((qdmw_nbq(i,j,k)-rw_nbq(i,j,k))
     &                              /dtnbq-ndtnbq*rw_int_nbq(i,j,k))
     &                                            /(pm(i,j)*pn(i,j))
            enddo
          enddo
        enddo
      endif

!# if defined RVTK_DEBUG && defined NBQ
!      call check_tab3d(rw_nbq,'rw_nbq','rint')
!#  endif  
# endif
!
!=====================================================================
! Get filtered RHS terms
! and multiply by dx*dy to get units of rho*Hz*dx*dy*ru
!=====================================================================
!
      if (LAST_FAST_STEP) then
!
!-----------------------------------------------------------------------
! Store average fields AVG1 of rho and rhobar
!-----------------------------------------------------------------------
!     
# ifdef NBQ_MASS
        do j=Jstr,Jend
          do i=Istr,Iend
            rhobar_nbq_avg1(i,j)=rhobar_nbq(i,j,knew)
          enddo
        enddo 
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
              rho_nbq_avg1(i,j,k)=1.d0
     &                            + ( rho_nbq(i,j,k)/Hzr(i,j,k)
     &                               +rho_grd(i,j,k) )  
            enddo
          enddo 
        enddo
# endif /* NBQ_MASS */
!
!-----------------------------------------------------------------------
!  Compute average fields AVG2 of RHS NBQ forcing
!  Note: here ru_nbq_avg2, ru_nbq_2d_old ... are working arrays
!-----------------------------------------------------------------------
!
        do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend
!             ru_nbq_avg1(i,j,k)= ru_nbq(i,j,k)
# ifdef NBQ
              ru_int_nbq(i,j,k) = ru_int_nbq(i,j,k)
     &                      -ru_ext_nbq_old(i,j)*(Hz(i-1,j,k)+Hz(i,j,k))
# endif
              ru_nbq_avg2(i,j,k)=
     &                          ((qdmu_nbq(i,j,k)-ru_nbq_avg2(i,j,k))/dt
     &                   -ru_int_nbq(i,j,k)-(ru_ext_nbq_sum(i,j)/nfast)*
     &                      (Hz(i,j,k)+Hz(i-1,j,k)))*on_u(i,j)*om_u(i,j)
            enddo
          enddo 
        enddo          

        do k=1,N
          do j=JstrV,Jend
            do i=Istr,Iend             
!             rv_nbq_avg1(i,j,k)= rv_nbq(i,j,k) 
# ifdef NBQ
              rv_int_nbq(i,j,k) = rv_int_nbq(i,j,k)
     &                      -rv_ext_nbq_old(i,j)*(Hz(i,j-1,k)+Hz(i,j,k))
# endif
              rv_nbq_avg2(i,j,k)=
     &                          ((qdmv_nbq(i,j,k)-rv_nbq_avg2(i,j,k))/dt
     &                   -rv_int_nbq(i,j,k)-(rv_ext_nbq_sum(i,j)/nfast)*
     &                      (Hz(i,j,k)+Hz(i,j-1,k)))*on_v(i,j)*om_v(i,j)
            enddo
          enddo 
        enddo

# ifdef NBQ
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
!             rw_nbq_avg1(i,j,k)= rw_nbq(i,j,k)
              rw_nbq_avg2(i,j,k)=
     &                          ((qdmw_nbq(i,j,k)-rw_nbq_avg2(i,j,k))/dt
     &                           -rw_int_nbq(i,j,k))*on_r(i,j)*om_r(i,j)
            enddo
          enddo 
        enddo
# endif

      endif ! LAST_FAST_STEP
!
!-----------------------------------------------------------------------
! Dismiss coupling of NBQ, NBQ2EXT & NBQ2INT for testing
!-----------------------------------------------------------------------
!
# ifdef NBQ_NOCOUPLING 
      ru_nbq      =0.   ! 3D
      rv_nbq      =0.
      ru_nbq_avg2 =0.
      rv_nbq_avg2 =0.
#  ifdef NBQ
      rw_nbq      =0.
      rw_nbq_avg2 =0.
#  endif
#  ifdef NBQ_MASS
      rhobar_nbq  =1.
      rho_nbq     =0.
      rho_nbq_avg1=1.
#  endif
# endif
!
!-----------------------------------------------------------------------
! Exchange NBQ coupling
!-----------------------------------------------------------------------
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
#  ifdef NBQ_MASS
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        rhobar_nbq_avg1(START_2D_ARRAY))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        rho_nbq_avg1(START_2D_ARRAY,1))
#  endif
#  ifndef NBQ
      IF (LAST_FAST_STEP) then
#  endif
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,  
     &                        ru_nbq_avg2(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,  
     &                        rv_nbq_avg2(START_2D_ARRAY,1))
#  ifndef NBQ
      endif
#  endif
#  ifdef NBQ
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,  
     &                        rw_nbq_avg2(START_2D_ARRAY,0))
#  endif 
# endif
!
!=====================================================================
!  Depth-averaged velocity & forcing from fast mode
!=====================================================================
!
! Output: (ubar,vbar), (DU_avg1,DV_avg1)
!
# define Dstp DUon

      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          Dstp(i,j)=zeta(i,j,kstp)+h(i,j)
        enddo
      enddo
      
      cff=2*dtfast
      cff1=0.5*weight(1,iif)
      cff2=0.5*weight(2,iif)

      if (iif.eq.1) then
         DU_avg1=0.
         DV_avg1=0.
         DU_avg2=0.
         DV_avg2=0.
      endif

      do j=Jstr,Jend
        do i=IstrU,Iend
# ifdef NBQ
          DUnew=DU_nbq(i,j) *2.
# else
          DUnew=( (Dstp(i,j)+Dstp(i-1,j))*ubar(i,j,kstp)
     &                             + cff * rubar_nbq(i,j) )
#ifdef MASKING
     &                                     *umask(i,j)
#endif
# endif
# ifdef WET_DRY
     &                                     *umask_wet(i,j)
# endif
          ubar(i,j,knew)=DUnew/(Dnew(i,j)+Dnew(i-1,j))
          DU_avg1(i,j,nnew)=DU_avg1(i,j,nnew) 
     &                             +cff1*on_u(i,j)*( DUnew
# ifdef MRL_WCI
     &                +(Dnew(i,j)+Dnew(i-1,j))*ust2d(i,j)
#  ifdef WET_DRY
     &  *umask_wet(i,j)
#  endif
# endif
     &                                                   )
          DU_avg2(i,j)=DU_avg2(i,j)+cff2*on_u(i,j)*( DUnew 
# ifdef MRL_WCI
     &                +(Dnew(i,j)+Dnew(i-1,j))*ust2d(i,j)
#  ifdef WET_DRY
     & *umask_wet(i,j)
#  endif
# endif
     &                                                   )
        enddo
      enddo 
  
      do j=JstrV,Jend
        do i=Istr,Iend
# ifdef NBQ
          DVnew=DV_nbq(i,j) *2.
# else
          DVnew=( (Dstp(i,j)+Dstp(i,j-1))*vbar(i,j,kstp)
     &                             + cff * rvbar_nbq(i,j) )
#ifdef MASKING
     &                                     *vmask(i,j)
#endif
# endif
# ifdef WET_DRY
     &                                     *vmask_wet(i,j)
# endif
          vbar(i,j,knew)=DVnew/(Dnew(i,j)+Dnew(i,j-1))
          DV_avg1(i,j,nnew)=DV_avg1(i,j,nnew) 
     &                              +cff1*om_v(i,j)*(DVnew
# ifdef MRL_WCI
     &                +(Dnew(i,j)+Dnew(i,j-1))*vst2d(i,j)
#  ifdef WET_DRY
     & *vmask_wet(i,j)
#  endif
# endif
     &                                                   )
          DV_avg2(i,j)=DV_avg2(i,j)+cff2*om_v(i,j)*( DVnew
# ifdef MRL_WCI
     &                +(Dnew(i,j)+Dnew(i,j-1))*vst2d(i,j)
#  ifdef WET_DRY
     & * vmask_wet(i,j)
#  endif
# endif
     &                                                   )
        enddo
      enddo
!
!-----------------------------------------------------------------------
! Apply point sources for hydrostatic case
!-----------------------------------------------------------------------
!
# if defined PSOURCE && !defined NBQ
      do is=1,Nsrc 
#  ifdef MPI
        i=Isrc_mpi(is,mynode)
        j=Jsrc_mpi(is,mynode)
#  else
        i=Isrc(is)
        j=Jsrc(is)
#  endif
        if (IstrR.le.i .and. i.le.IendR .and.
     &      JstrR.le.j .and. j.le.JendR) then
          if (Dsrc(is).eq.0) then
            ubar(i,j,knew)=2.*Qbar(is)/( on_u(i,j)
     &                       *(Dnew(i-1,j)+Dnew(i,j)) )
            DU_avg1(i,j,nnew)=Qbar(is)
          else
            vbar(i,j,knew)=2.*Qbar(is)/( om_v(i,j)
     &                       *(Dnew(i,j-1)+Dnew(i,j)) )
            DV_avg1(i,j,nnew)=Qbar(is)
          endif
        endif
      enddo
# endif
!
!-----------------------------------------------------------------------
!  Set 2D Momemtum nudging
!-----------------------------------------------------------------------
!
# if defined M2NUDGING && defined M2CLIMATOLOGY

#  ifdef ZONAL_NUDGING
      if (FIRST_TIME_STEP .or. mod(iic,10).eq.0) then
        if (FIRST_FAST_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   ubar(START_2D_ARRAY,knew),ubzon)
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   vbar(START_2D_ARRAY,knew),vbzon)
        endif
      endif
      if (FIRST_TIME_STEP) then
        if (FIRST_FAST_STEP) then 
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   ubclm(START_2D_ARRAY),ubclmzon)
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   vbclm(START_2D_ARRAY),vbclmzon)
        endif
      endif
#  endif /* ZONAL_NUDGING */

      do j=Jstr,Jend
        do i=IstrU,Iend
#  ifdef ZONAL_NUDGING        
          DUnew = dtfast*M2nudgcof(i,j)*(ubclmzon(j)-ubzon(j))
#  else          
          DUnew = dtfast*M2nudgcof(i,j)*(ubclm(i,j)-ubar(i,j,knew))
#  endif
#  ifdef MASKING
     &                 * umask(i,j)
#  endif 
#  ifdef WET_DRY
     &                 * umask_wet(i,j)
#  endif
          ubar(i,j,knew)=ubar(i,j,knew) + DUnew
          DU_avg1(i,j,nnew)=DU_avg1(i,j,nnew) +cff1*DUnew*
     &                         (Dnew(i,j)+Dnew(i-1,j))*on_u(i,j)
          DU_avg2(i,j)     =DU_avg2(i,j)      +cff2*DUnew*
     &                         (Dnew(i,j)+Dnew(i-1,j))*on_u(i,j)
        enddo
      enddo
      
      do j=JstrV,Jend
        do i=Istr,Iend
#  if defined ZONAL_NUDGING  
          DVnew = dtfast*M2nudgcof(i,j)*(vbclmzon(j)-vbzon(j)) 
#  else
          DVnew = dtfast*M2nudgcof(i,j)*(vbclm(i,j)-vbar(i,j,knew))    
#  endif      
#  ifdef MASKING
     &                 * vmask(i,j)
#  endif 
#  ifdef WET_DRY
     &                 * vmask_wet(i,j)
#  endif
          vbar(i,j,knew)=vbar(i,j,knew) + DVnew
          DV_avg1(i,j,nnew)=DV_avg1(i,j,nnew) +cff1*DVnew*
     &                         (Dnew(i,j)+Dnew(i,j-1))*om_v(i,j)
          DV_avg2(i,j)     =DV_avg2(i,j)      +cff2*DVnew*
     &                         (Dnew(i,j)+Dnew(i,j-1))*om_v(i,j)
        enddo
      enddo
# endif /* M2NUDGING */
!
!----------------------------------------------------------------------- 
! Set boundary conditions and compute integral mass flux accross
! all open boundaries, if any.
!-----------------------------------------------------------------------
!
# ifdef NBQ
      M2bc_nbq_flag=.false.  ! skip wet/dry conditions
                             ! and DU_nbq computation
# endif
      call u2dbc_tile (Istr,Iend,Jstr,Jend, work) 
      call v2dbc_tile (Istr,Iend,Jstr,Jend, work)

# ifdef WET_DRY
#  ifndef EW_COM_PERIODIC
      if (WESTERN_EDGE) then
        DO j=Jstr,Jend
          ubar(Istr,j,knew)=ubar(Istr,j,knew)*umask_wet(Istr,j)
#   ifdef MRL_WCI
          ust2d(Istr,j)=ust2d(Istr,j)*umask_wet(Istr,j)
#   endif
        END DO
        DO j=JstrV,Jend
          vbar(Istr-1,j,knew)=vbar(Istr-1,j,knew)*vmask_wet(Istr-1,j)
#   ifdef MRL_WCI
          vst2d(Istr-1,j)=vst2d(Istr-1,j)*vmask_wet(Istr-1,j)
#   endif
        END DO
      END IF
      if (EASTERN_EDGE) then
        DO j=Jstr,Jend
          ubar(Iend+1,j,knew)=ubar(Iend+1,j,knew)*umask_wet(Iend+1,j)
#   ifdef MRL_WCI
          ust2d(Iend+1,j)=ust2d(Iend+1,j)*umask_wet(Iend+1,j)
#   endif
        END DO
        DO j=JstrV,Jend
          vbar(Iend+1,j,knew)=vbar(Iend+1,j,knew)*vmask_wet(Iend+1,j)
#   ifdef MRL_WCI
          vst2d(Iend+1,j)=vst2d(Iend+1,j)*vmask_wet(Iend+1,j)
#   endif
        END DO
      END IF
#  endif
#  ifndef NS_COM_PERIODIC
      if (SOUTHERN_EDGE) then
        DO i=IstrU,Iend
          ubar(i,Jstr-1,knew)=ubar(i,Jstr-1,knew)*umask_wet(i,Jstr-1)
#   ifdef MRL_WCI
          ust2d(i,Jstr-1)=ust2d(i,Jstr-1)*umask_wet(i,Jstr-1)
#   endif
        END DO
        DO i=IstrU,Iend
          vbar(i,Jstr,knew)=vbar(i,Jstr,knew)*vmask_wet(i,Jstr)
#   ifdef MRL_WCI
          vst2d(i,Jstr)=vst2d(i,Jstr)*vmask_wet(i,Jstr)
#   endif
        END DO
      END IF
      if (NORTHERN_EDGE) then
        DO i=Istr,Iend
          ubar(i,Jend+1,knew)=ubar(i,Jend+1,knew)*umask_wet(i,Jend+1)
#   ifdef MRL_WCI
          ust2d(i,Jend+1)=ust2d(i,Jend+1)*umask_wet(i,Jend+1)
#   endif
        END DO
        DO i=Istr,Iend
          vbar(i,Jend+1,knew)=vbar(i,Jend+1,knew)*vmask_wet(i,Jend+1)
#   ifdef MRL_WCI
          vst2d(i,Jend+1)=vst2d(i,Jend+1)*vmask_wet(i,Jend+1)
#   endif
        END DO
      END IF
#  endif
# endif
!
! zeta vill be recomputed via depth-integrated continuity equation
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        ubar(START_2D_ARRAY,knew))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        vbar(START_2D_ARRAY,knew))
# endif
          
# ifdef OBC_VOLCONS
      call obc_flux_tile (Istr,Iend,Jstr,Jend)
# endif

!      call unbq_bc_tile (Istr,Iend,Jstr,Jend, work)
!      call vnbq_bc_tile (Istr,Iend,Jstr,Jend, work)
!#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI  
!      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
!     &                        qdmu_nbq(START_2D_ARRAY,1))
!      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
!     &                        qdmv_nbq(START_2D_ARRAY,1))
!#  endif
!
!-----------------------------------------------------------------------
! Compute fast-time averaged barotropic mass fluxes along physical
! boundaries.
!-----------------------------------------------------------------------
!
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=Jstr-1,JendR
          Dnew(Istr-1,j)=(h(Istr-1,j)+zeta(Istr-1,j,knew))
#  ifdef NBQ_MASS
     &                          *rhobar_nbq(Istr-1,j,knew)
#  endif
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=Jstr-1,JendR
          Dnew(Iend+1,j)=(h(Iend+1,j)+zeta(Iend+1,j,knew))
#  ifdef NBQ_MASS
     &                          *rhobar_nbq(Iend+1,j,knew)
#  endif
        enddo
      endif
# endif /* !EW_PERIODIC */

# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=Istr-1,IendR
          Dnew(i,Jstr-1)=(h(i,Jstr-1)+zeta(i,Jstr-1,knew))
#  ifdef NBQ_MASS
     &                          *rhobar_nbq(i,Jstr-1,knew)
#  endif
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=Istr-1,IendR
          Dnew(i,Jend+1)=(h(i,Jend+1)+zeta(i,Jend+1,knew))
#  ifdef NBQ_MASS
     &                          *rhobar_nbq(i,Jend+1,knew)
#  endif
        enddo
      endif
# endif /* !NS_PERIODIC */

      cff1=0.5*weight(1,iif)
      cff2=0.5*weight(2,iif)
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=JstrR,JendR
          cff=(Dnew(IstrU-1,j)+Dnew(IstrU-2,j))*(ubar(IstrU-1,j,knew)
# ifdef MRL_WCI
     &                                              +ust2d(IstrU-1,j)
# endif
     &                                              )*on_u(IstrU-1,j)
          DU_avg1(IstrU-1,j,nnew)=DU_avg1(IstrU-1,j,nnew)+cff1*cff
          DU_avg2(IstrU-1,j)=DU_avg2(IstrU-1,j)+cff2*cff
        enddo
        do j=JstrV,Jend
          cff=(Dnew(Istr-1,j)+Dnew(Istr-1,j-1) )*(vbar(Istr-1,j,knew)
# ifdef MRL_WCI
     &                                               +vst2d(Istr-1,j)
# endif
     &                                               )*om_v(Istr-1,j)
          DV_avg1(Istr-1,j,nnew)=DV_avg1(Istr-1,j,nnew)+cff1*cff
          DV_avg2(Istr-1,j)=DV_avg2(Istr-1,j)+cff2*cff
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=JstrR,JendR
          cff=(Dnew(Iend+1,j)+Dnew(Iend,j))*(ubar(Iend+1,j,knew)
# ifdef MRL_WCI
     &                                          +ust2d(Iend+1,j)
# endif
     &                                          )*on_u(Iend+1,j)
          DU_avg1(Iend+1,j,nnew)=DU_avg1(Iend+1,j,nnew)+cff1*cff
          DU_avg2(Iend+1,j)=DU_avg2(Iend+1,j)+cff2*cff
        enddo
        do j=JstrV,Jend
          cff=(Dnew(Iend+1,j)+Dnew(Iend+1,j-1))*(vbar(Iend+1,j,knew)
# ifdef MRL_WCI
     &                                              +vst2d(Iend+1,j)
# endif
     &                                              )*om_v(Iend+1,j)
          DV_avg1(Iend+1,j,nnew)=DV_avg1(Iend+1,j,nnew)+cff1*cff
          DV_avg2(Iend+1,j)=DV_avg2(Iend+1,j)+cff2*cff
        enddo
      endif
# endif /* !EW_PERIODIC */
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=IstrU,Iend
          cff=(Dnew(i,Jstr-1)+Dnew(i-1,Jstr-1))*(ubar(i,Jstr-1,knew)
#  ifdef MRL_WCI
     &                                              +ust2d(i,Jstr-1)
#  endif
     &                                              )*on_u(i,Jstr-1)
          DU_avg1(i,Jstr-1,nnew)=DU_avg1(i,Jstr-1,nnew)+cff1*cff
          DU_avg2(i,Jstr-1)=DU_avg2(i,Jstr-1)+cff2*cff
        enddo
        do i=IstrR,IendR
          cff=(Dnew(i,JstrV-1)+Dnew(i,JstrV-2))*(vbar(i,JstrV-1,knew)
#  ifdef MRL_WCI
     &                                              +vst2d(i,JstrV-1)
#  endif
     &                                              )*om_v(i,JstrV-1)
          DV_avg1(i,JstrV-1,nnew)=DV_avg1(i,JstrV-1,nnew)+cff1*cff
          DV_avg2(i,JstrV-1)=DV_avg2(i,JstrV-1)+cff2*cff
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=IstrU,Iend
          cff=(Dnew(i,Jend+1)+Dnew(i-1,Jend+1))*(ubar(i,Jend+1,knew)
#  ifdef MRL_WCI
     &                                              +ust2d(i,Jend+1)
#  endif
     &                                              )*on_u(i,Jend+1)
          DU_avg1(i,Jend+1,nnew)=DU_avg1(i,Jend+1,nnew)+cff1*cff
          DU_avg2(i,Jend+1)=DU_avg2(i,Jend+1)+cff2*cff
        enddo
        do i=IstrR,IendR
          cff=(Dnew(i,Jend+1)+Dnew(i,Jend))*(vbar(i,Jend+1,knew)
#  ifdef MRL_WCI
     &                                          +vst2d(i,Jend+1)
#  endif
     &                                          )*om_v(i,Jend+1)
          DV_avg1(i,Jend+1,nnew)=DV_avg1(i,Jend+1,nnew)+cff1*cff
          DV_avg2(i,Jend+1)=DV_avg2(i,Jend+1)+cff2*cff
        enddo
      endif
# endif /* !NS_PERIODIC */
!
!=====================================================================
!
!  Adjust ZETA using depth-integrated continuity equation 
!  and update grid
!
!  Once rhobar_nbq and depth-averaged momentum is updated, surface anomalies
!  can be adjusted to satisfy the low-frequency mass conservation 
!  equation. However, this adjustment do not satisfy conservation at machine 
!  precision and a final correction is needed.
!
!  This recomputation of zeta(m+1) using div(ubar) is only done for PRECISE
!  option. In PERF option, zeta(m+1) is computed only once based on the 
!  surface characteristic relation. But the final numerical correction 
!  is applied in all cases.
!
!=====================================================================
!
# ifdef NBQ
!-----------------------------------------------------------------------
!  Update zeta(m+1)
!-----------------------------------------------------------------------
!
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            Dnew(i,j)=(zeta(i,j,knew)+h(i,j))
#  ifdef NBQ_MASS
     &                                          *rhobar_nbq(i,j,knew)
#  endif
          enddo
        enddo

        do j=Jstr,Jend
          do i=Istr,Iend
            zeta(i,j,knew)=( (h(i,j)+zeta(i,j,kstp))
#  ifdef NBQ_MASS
     &                                          *rhobar_nbq(i,j,kstp)
#  endif
     &                                + (dtfast*pm(i,j)*pn(i,j)*0.5*(
     &      (Dnew(i  ,j)+Dnew(i-1,j))*(ubar(i  ,j,knew)
#  ifdef MRL_WCI
     &                               +ust2d(i  ,j)
#  endif  
     &                                                 )*on_u(i  ,j)
     &     -(Dnew(i+1,j)+Dnew(i  ,j))*(ubar(i+1,j,knew)
#  ifdef MRL_WCI
     &                               +ust2d(i+1,j)
#  endif 
     &                                                 )*on_u(i+1,j)
     &     +(Dnew(i,j  )+Dnew(i,j-1))*(vbar(i,j  ,knew)
#  ifdef MRL_WCI
     &                               +vst2d(i,j  )
#  endif
     &                                                 )*om_v(i,j  )
     &     -(Dnew(i,j+1)+Dnew(i,j  ))*(vbar(i,j+1,knew)
#  ifdef MRL_WCI
     &                               +vst2d(i,j+1)
#  endif
     &                                                 )*om_v(i,j+1))) )
#  ifdef NBQ_MASS
     &                                             /rhobar_nbq(i,j,knew)
#  endif
     &                                                          - h(i,j)
          enddo
        enddo
        

!
!-----------------------------------------------------------------------
!  Set masking for zeta, including wet/dry conditions
!-----------------------------------------------------------------------
!
#  ifdef MASKING
      do j=Jstr,Jend
        do i=Istr,Iend
          zeta(i,j,knew)=zeta(i,j,knew)*rmask(i,j)
#   ifdef WET_DRY
!    modify new free-surface to ensure that depth 
!    is > Dcrit in masked cells.
          cff=0.5+SIGN(0.5,Dcrit(i,j)-h(i,j))
          zeta(i,j,knew)=zeta(i,j,knew)+ 
     &                   cff*(Dcrit(i,j)-h(i,j))*(1.-rmask(i,j))
#   endif
        enddo
      enddo 
#  endif /* MASKING */
!
!-----------------------------------------------------------------------
!  Set boundary conditions for the free-surface
!  --> ensure closed boundaries
!-----------------------------------------------------------------------
!
#  ifndef OBC_NBQ
      call zetabc_tile (Istr,Iend,Jstr,Jend)
#  endif

# endif /* NBQ */
!
!-----------------------------------------------------------------------
!  Update Zt_avg1 at last fast step
!-----------------------------------------------------------------------
!
      if (LAST_FAST_STEP) then
        do j=JstrR,JendR
          do i=IstrR,IendR
            Zt_avg1(i,j)=zeta(i,j,knew)
          enddo
        enddo
!#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
!        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
!     &                          Zt_avg1(START_2D_ARRAY))
!#  endif
      endif
!
!-----------------------------------------------------------------------
! Update grid parameters at m+1: Hz, z_r, z_w
! in prognostic or diagnostic way
!-----------------------------------------------------------------------
!
# ifdef NBQ
#  ifdef NBQ_GRID_SLOW
      if (LAST_FAST_STEP) then
#  endif

#  ifdef NBQ_HZ_PROGNOSTIC
!
!  Prognostic evaluation using momentum divergence
!
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
            do k=1,N
              Hz(i,j,k)=Hz_bak2(i,j,k) - dtfast*(thetadiv_nbq(i,j,k)
     &                                         +thetadiv2_nbq(i,j,k))
              Hzr(i,j,k)=(Hz(i,j,k)-rho_nbq(i,j,k))/(1.+rho(i,j,k)/rho0)
              z_w(i,j,k)=z_w(i,j,k-1)+Hzr(i,j,k)
              z_r(i,j,k)=0.5*(z_w(i,j,k)+z_w(i,j,k-1))
            enddo
          enddo
        enddo
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI      
        call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                          z_w(START_2D_ARRAY,0))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          z_r(START_2D_ARRAY,1))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          Hz(START_2D_ARRAY,1))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          Hzr(START_2D_ARRAY,1))
#   endif
#  else
!
!  Diagnostic evaluation from zeta(m+1)
!
        call set_depth_tile(Istr,Iend,Jstr,Jend)
#  endif
!
! Compute derived grid parameters if fast update
!
#  ifndef NBQ_GRID_SLOW
        call grid_nbq_tile(Istr,Iend,Jstr,Jend,
     &                     Hzw_half_nbq_inv,   Hzr_half_nbq_inv,
     &                     Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v,
     &                     Hzu_half_qdmu,      Hzv_half_qdmv)
#  endif

#  ifdef NBQ_GRID_SLOW
      endif !<-- LAST_FAST_STEP
#  endif
# endif /* NBQ */
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        zeta(START_2D_ARRAY,knew))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        DU_avg1(START_2D_ARRAY,nnew))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        DV_avg1(START_2D_ARRAY,nnew))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        DU_avg2(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        DV_avg2(START_2D_ARRAY))
#  if defined MRL_WCI && defined WET_DRY
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        ust2d(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        vst2d(START_2D_ARRAY))
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Apply conservation requirements for nesting
!-----------------------------------------------------------------------
!
# ifdef AGRIF
      if (.NOT.Agrif_Root()) THEN
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            Zt_avg3(i,j,iif)=zeta(i,j,knew)
          enddo
        enddo 
        do j=JstrR,JendR
          do i=Istr,IendR
            DU_avg3(i,j,iif) = 0.5*(h(i  ,j)+zeta(i  ,j,knew)+
     &                              h(i-1,j)+zeta(i-1,j,knew)) 
     &                               *on_u(i,j)*ubar(i,j,knew)
          enddo
        enddo 
        do j=Jstr,JendR
          do i=IstrR,IendR
            DV_avg3(i,j,iif) = 0.5*(h(i,j  )+zeta(i,j  ,knew)+
     &                              h(i,j-1)+zeta(i,j-1,knew)) 
     &                               *om_v(i,j)*vbar(i,j,knew)
          enddo
        enddo
 
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                          Zt_avg3(START_2D_ARRAY,iif))
        call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                          DU_avg3(START_2D_ARRAY,iif))
        call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                          DV_avg3(START_2D_ARRAY,iif))
#  endif

#  ifdef RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
        call check_tab2d(Zt_avg3(:,:,iif),'Zt_avg3 step3d_fast','r')
        call check_tab2d(DU_avg3(:,:,iif),'DU_avg3 step3d_fast','u')
        call check_tab2d(DV_avg3(:,:,iif),'DV_avg3 step3d_fast','v')
C$OMP END MASTER
#  endif   
      endif

#  ifdef AGRIF_CONSERV_VOL
      if (iif.eq.nfast) then
        if (agrif_root()) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              DU_avg1(i,j,5) = dt * DU_avg2(i,j)
              DV_avg1(i,j,5) = dt * DV_avg2(i,j)
            enddo
          enddo
        else
          do j=JstrR,JendR
            do i=IstrR,IendR
              DU_avg1(i,j,5) = dt * DU_avg2(i,j)
              DV_avg1(i,j,5) = dt * DV_avg2(i,j)
              DU_avg1(i,j,4) = DU_avg1(i,j,4) + DU_avg1(i,j,5)
              DV_avg1(i,j,4) = DV_avg1(i,j,4) + DV_avg1(i,j,5)
            enddo
          enddo       
        endif
      endif
#  endif
# endif /* AGRIF */
!
!-----------------------------------------------------------------------
!  Copy density for extrapolation: 
!-----------------------------------------------------------------------
!
# ifdef NBQ_MASS
      if (LAST_FAST_STEP) then
        do k=1,N
          do j=JstrV-2,Jend+1
            do i=IstrU-2,Iend+1
              rho_bak(i,j,k)=rho(i,j,k)
            enddo
           enddo
         enddo
      endif
# endif
!
!-----------------------------------------------------------------------
!  Correct Hz(m+1) for internal mode
!  by inverting internal continuity equation
!-----------------------------------------------------------------------
!
# ifdef NBQ_HZCORRECT
     
      if (LAST_FAST_STEP) then

        do k=1,N
          do j=JstrV-2,Jend+1
            do i=IstrU-2,Iend+1
              Hz_correct(i,j,k)=Hz(i,j,k)
            enddo
          enddo
        enddo

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
        call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                          DU_avg2(START_2D_ARRAY))
        call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                          DV_avg2(START_2D_ARRAY))
#  endif

        do j=Jstr,Jend
          do i=Istr,Iend
            dum_s=0.
            do k=1,N
              dum_s=dum_s+Hz(i,j,k)-Hz_bak(i,j,k)
            enddo
            do k=1,N
              Hz(i,j,k)=Hz(i,j,k)
     &                -(dum_s
     &             + (DU_avg2(i+1,j)-DU_avg2(i,j)
     &               +DV_avg2(i,j+1)-DV_avg2(i,j)
     &               )*pm(i,j)*pn(i,j)
     &               *dt)
     &               /(z_w(i,j,N)-z_w(i,j,0))
     &               *(z_w(i,j,k)-z_w(i,j,k-1))
#  ifdef MASKING
     &               *rmask(i,j)
#  endif

#  ifdef NBQ_HZCORR_DEBUG
              Hz_corr(i,j,k)=Hz(i,j,k)-Hz_correct(i,j,k)
#  endif

!# ifdef NBQ_MASS
!              Hzr(i,j,k)=(Hz(i,j,k)-rho_nbq(i,j,k)) ! adjust Zw,Zr?
!     &                        /(1.+rho(i,j,k)/rho0)
!# endif
!              z_w(i,j,k)=z_w(i,j,k-1)+Hzr(i,j,k)
!              z_r(i,j,k)=0.5*(z_w(i,j,k)+z_w(i,j,k-1))

            enddo
          enddo
        enddo

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          Hz(START_2D_ARRAY,1))
!# ifdef NBQ_MASS
!         call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
!     &                           Hzr(START_2D_ARRAY,1))
!# endif
!        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
!     &                          z_w(START_2D_ARRAY,0))
!        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
!     &                          z_r(START_2D_ARRAY,1))
#  endif 

      endif ! iif.eq.nfast

# endif  /* NBQ_HZCORRECT */
!
!-----------------------------------------------------------------------
!  TEST FOR CFL VIOLATION. IF SO, PRINT AND STOP
!-----------------------------------------------------------------------
!
      VMAXL=100.
      VMAX=0.
      do j=Jstr,Jend
        do i=Istr,Iend
          cff1=ubar(i,j,knew)
          cff2=vbar(i,j,knew)
          cff=max(abs(cff1),abs(cff2))
          IF (cff.GE.VMAX .or. cff1.ne.cff1 .or. cff2.ne.cff2) THEN
            IF (cff.GE.VMAX .and. cff1.eq.cff1 .and. cff2.eq.cff2) THEN 
              VMAX=cff
            ELSE
              VMAX=666.
            ENDIF
# ifdef MPI
            imax=i+iminmpi-1
            jmax=j+jminmpi-1
# else
            imax=i
            jmax=j
# endif
          ENDIF
        enddo
      enddo
      IF (VMAX.GT.VMAXL) THEN
        write(stdout,'(9(A/))')
     &     '                                         ',
     &     '                                         ',
     &     ' ======================================= ',
     &     ' =                                     = ',
     &     ' =  STEP3D_FAST:   ABNORMAL JOB END    = ',
     &     ' =                 BLOW UP             = ',
     &     ' =                                     = ',
     &     ' ======================================= ',
     &     '                                         '

# ifdef MPI
        write(stdout,'(A,I4)') '      mynode =',mynode
# endif

        if (VMAX.eq.666.) then
          write(stdout,'(A,F10.2)')
     &                                            '  VMAX (M/S) =   NaN'
        else
          write(stdout,'(A,F10.2)') 
     &                                            '  VMAX (M/S) =',VMAX
        endif
        write(stdout,'(A,2I6)')  
     &                                       '  IMAX JMAX   =',imax,jmax
        write(stdout,'(A,2I6/)')
     &                                       '  IINT IFAST  =',iic,iif

!     GC quik fix Comment this because it blocks the correct mpi_abort
!     croco do not hang out  
!!      call wrt_his            ! get output during blow-up

        may_day_flag=1
# ifdef MPI
        call mpi_abort (MPI_COMM_WORLD, err)
# else
        stop                            !-->  EXIT
# endif
      ENDIF

# undef zwrk 
# undef rzeta  
# undef rzeta2  
# undef rzetaSA 

      return
      end
#else
      subroutine step3d_fast_empty
      end
#endif /* M3FAST */
