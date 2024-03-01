! $Id: gls_mixing.F 1524 2014-04-14 17:00:06Z gcambon $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
#include "cppdefs.h"

#if defined SOLVE3D && defined GLS_MIXING

# define PSI_SMOOTH
# define TKE_SMOOTH

# define GLS_HADV_UP1
# define GLS_VADV_UP1

# ifdef GLS_MIXING_3D
#  define GLS_SHEAR_3D
#  define GLS_STRAIN_LIM
# endif
!
!=======================================================================

      SUBROUTINE gls_mixing (tile)

!=======================================================================
!
      IMPLICIT NONE
      INTEGER         :: tile, trd
      INTEGER         :: omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "ocean3d.h"
# include "compute_tile_bounds.h"

      trd=omp_get_thread_num()

      call gls_mixing_tile ( Istr, Iend, Jstr, Jend,
     &                    A3d(1, 1,trd), A3d(1, 2,trd), A3d(1, 3,trd),
     &                    A3d(1, 4,trd), A3d(1, 5,trd), A2d(1, 2,trd), 
     &                    A2d(1, 3,trd), A2d(1, 4,trd), A2d(1, 5,trd),
     &                    A2d(1, 6,trd), A2d(1, 7,trd), A2d(1, 8,trd))

      RETURN
      END
      SUBROUTINE gls_mixing_tile ( Istr, Iend, Jstr, Jend,
     &                          tke_new, gls_new, shear2, vort2, work,
     &                                diss, ustar_sfc_sq,ustar_bot_sq,
     &                                                 DC, FC, CF, RH)

!
!======================================================================
!                 ***  SUBROUTINE  gls_mixing  ***
!
!   Vertical mixing coefficients computed from the gls turbulent
!  closure parameterization (3 schemes: k-epsilon, k-omega, k-gen)
!
! *** References :
!
! Umlauf, L., and H. Burchard, A generic length-scale equation
! for geophysical turbulence models, J. Mar. Res., 61, 235-265, 2003
!
! Umlauf, L., and H. Burchard, Second-order turbulence closure models
!                        Cont. Shelf. Res., 25, 795-827, 2005
!
! Duhaut, T., Notes sur les schemas de turbulence a deux equations,
!                     Technical Report. IFREMER report , 2009
!
! *** History :
!   ! 2016-11 (F. Lemarié) : original code inspired by the Mars3D
!                            implementation by S. Petton, A.-C. Bennis,
!                             T. Duhaut, V. Garnier, F. Dumas
!   ! 2020 (Marchesiello, Benshila): added 3D shear + advection
!                        + wave forcing (Kumar et al COAWST paper 2012)
!   ! 2021 (Marchesiello): MPI bug fixes (advection, wave forcing
!                                                       (wkb_wwave.F));
!      + revision of parameters (background visc) including for SANDBAR
!      + fix OPENMP parallelization
!
!======================================================================
!
      IMPLICIT NONE
# include "param.h"
# include "mixing.h"
# include "scalars.h"

      ! Local integers
      INTEGER         ::   Istr, Iend, Jstr, Jend
      INTEGER         ::   i,       j,    k, tind, kref
      INTEGER         ::   imin, imax, jmin, jmax
      INTEGER         ::   ig,    ig1,  ig2
      ! Local arrays
      REAL            ::  tke_new     (PRIVATE_2D_SCRATCH_ARRAY,0:N  )
      REAL            ::  gls_new     (PRIVATE_2D_SCRATCH_ARRAY,0:N  )
      REAL            ::  work        (PRIVATE_2D_SCRATCH_ARRAY,0:N  )
      REAL            ::  shear2      (PRIVATE_2D_SCRATCH_ARRAY,0:N-1)
      REAL            ::  vort2       (PRIVATE_2D_SCRATCH_ARRAY,0:N-1)
      REAL            ::  diss        (PRIVATE_1D_SCRATCH_ARRAY,1:N-1)
      REAL            ::  ustar_sfc_sq(PRIVATE_2D_SCRATCH_ARRAY      )
      REAL            ::  ustar_bot_sq(PRIVATE_2D_SCRATCH_ARRAY      )
      REAL            ::  DC          (PRIVATE_1D_SCRATCH_ARRAY,0:N  )
      REAL            ::  FC          (PRIVATE_1D_SCRATCH_ARRAY,0:N  )
      REAL            ::  CF          (PRIVATE_1D_SCRATCH_ARRAY,1:N-1)
      REAL            ::  RH          (PRIVATE_1D_SCRATCH_ARRAY,1:N-1)
      ! Local scalars
      REAL            ::  cff , cff1 , cff2, cff3m, cff3p
      REAL            ::  invk, invG, Bprod, Sprod, epsilon
      REAL            ::  alpha_n, alpha_m, c_mu, c_mu_prim
      REAL            ::  alpha_n_min, alpha_m_max, cm0inv2, gls
      REAL            ::  flux_top, flux_bot, lgthsc, L_lim, du,dv,dw
      REAL            ::  trb_new, trb_sfc, trb_bot, z0_s, z0_b, gls_min
      REAL            ::  HUon_w, HVom_w, trb_min(2), Denom
      REAL            ::  su_r,sv_r
# ifdef MRL_WCI
      REAL            ::  cmu_fac, cb_wallE, gls_sigp_cb
# endif
      ! Parameter values
      REAL, PARAMETER ::  eps_min =  1.0E-12 ! min value for dissipation rate
      REAL, PARAMETER ::  tke_min =  1.0E-10 ! min value for TKE
# ifdef GLS_MIXING_3D
      REAL, PARAMETER ::  eps     =  1.0E-12 ! min stable stratification
      REAL, PARAMETER ::  nuws    =  0.1E-06 ! background diffusivity
      REAL, PARAMETER ::  nuwm    =  1.0E-06 ! background viscosity
# else
      REAL, PARAMETER ::  eps     =  1.0E-14 ! min stable stratification
      REAL, PARAMETER ::  nuws    =  0.1E-05 ! background diffusivity
      REAL, PARAMETER ::  nuwm    =  1.0E-05 ! background viscosity
# endif
      REAL, PARAMETER ::  galp    =  0.53    ! parameter for Galperin
                                             ! mixing length limitation
      REAL, PARAMETER ::  chk     =  1400./g ! charnock coefficient
# if defined SANDBAR || defined SHOREFACE
      REAL, PARAMETER ::  Zosmin  =  1.e-2   ! min surface roughness length
      REAL, PARAMETER ::  Zobmin  =  5.e-2   ! min bottom  roughness length
# else
      REAL, PARAMETER ::  Zosmin  =  1.e-2   ! min surface roughness length
      REAL, PARAMETER ::  Zobmin  =  1.e-4   ! min bottom  roughness length
# endif
# ifdef MRL_WCI
      REAL, PARAMETER ::  zalp    =  0.05    ! fraction of wave energy
# endif

      ! Choice of GLS model
      real            :: rp,    rm,    rn              !<-- n,m and p exponents
      real            :: beta1, beta2, beta3m, beta3p  !<-- beta terms for psi equation
      real            :: OneOverSig(2)   !<-- inverse of Schmidt number for tke and psi
#  if defined GLS_KOMEGA                            /* K-omega model */
      parameter( rp    = -1.0 , rm    = 0.5  , rn     = -1.0           )
      parameter( beta1 = 0.555, beta2 = 0.833, beta3m = -0.6, beta3p = 1.0)
      parameter( OneOverSig = (/ 0.5, 0.5 /) )       
#  elif defined GLS_KEPSILON                      /* K-epsilon model */
      parameter( rp    = 3.0 , rm    = 1.5 , rn     = -1.0             )
      parameter( beta1 = 1.44, beta2 = 1.92, beta3m = -0.4, beta3p = 1.0)
      parameter( OneOverSig = (/ 1.0, 0.7692 /) )  
#  else                                                 /* GEN model */
      parameter( rp    = 0.0, rm    = 1.0 , rn     = -0.67             )
      parameter( beta1 = 1.0, beta2 = 1.22, beta3m =  0.05, beta3p = 1.0)
      parameter( OneOverSig = (/ 1.25, 0.9345 /) ) 
#  endif
      REAL, PARAMETER :: e1 =  3.0 + 1.*rp / rn
      REAL, PARAMETER :: e2 =  1.5 + 1.*rm / rn
      REAL, PARAMETER :: e3 = -1.0 / rn
# if defined PSI_SMOOTH || defined TKE_SMOOTH
      ! 9-point isotropic laplacian filter
      REAL, PARAMETER :: smth_a = 1./12.
      REAL, PARAMETER :: smth_b = 3./16.
# endif
# ifdef GLS_STRAIN_LIM
      REAL, PARAMETER :: lars = 5.e-3*beta2/beta1  ! parameter for Larsen
                                                   ! dissipation limitation
# endif

# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
!
!--------------------------------------------------------------------
! Choice of Stability function (Default : Canuto A)
!--------------------------------------------------------------------
!
# include "gls_stab_func.h"
!
!--------------------------------------------------
# include "ocean3d.h"
# include "forces.h"
# include "grid.h"
# include "averages.h"
# ifdef BBL
#  include "bbl.h"
# endif
!
!--------------------------------------------------------------------
! Compute extended bounds (necessary because of
!                          the TKE/GLS smoothing)
!--------------------------------------------------------------------
!
# ifdef EW_PERIODIC
#  define I_EXT_RANGE Istr-1,Iend+1
      imin = Istr-1
      imax = Iend+1
# else
      if (WESTERN_EDGE) then
        imin=Istr
      else
        imin=Istr-1
      endif
      if (EASTERN_EDGE) then
        imax=Iend
      else
        imax=Iend+1
      endif
#  define I_EXT_RANGE imin,imax
# endif
# ifdef NS_PERIODIC
#  define J_EXT_RANGE Jstr-1,Jend+1
      jmin = Jstr-1
      jmax = Jend+1
# else
      if (SOUTHERN_EDGE) then
        jmin=Jstr
      else
        jmin=Jstr-1
      endif
      if (NORTHERN_EDGE) then
        jmax=Jend
      else
        jmax=Jend+1
      endif
#  define J_EXT_RANGE jmin,jmax
# endif
!
!--------------------------------------------------------------------
! Initialization of various constants
!--------------------------------------------------------------------
!
! Compute cmu0 and inverse of cmu0 squared
!
      cm0     =  ( (a2*a2 - 3.0*a3*a3 + 3.0*a1*nn)/(3.0*nn*nn) )**0.25
      cm0inv2 = 1./cm0**2
!
! Minimum value of alpha_n to ensure that alpha_m is positive
!
      alpha_n_min = 0.5*( - ( sf_d1 + sf_nb0 )
     &             + sqrt(  ( sf_d1 + sf_nb0 )**2
     &            - 4. * sf_d0 *( sf_d4 + sf_nb1 ) ) )
     &                            / ( sf_d4 + sf_nb1 )
!
! Compute gls_min consistently with eps_min/tke_min
!
      cff     = (cm0**3 )*(tke_min**1.5) / eps_min
      gls_min = (cm0**rp)*(tke_min**rm ) * ( cff**rn )

      trb_min(itke) = tke_min
      trb_min(igls) = gls_min

# ifdef MRL_WCI
      cmu_fac=(1.5/OneOverSig(1))**(1./3.)
     &             /(cm0**(4./3.))
      cb_wallE=1.
      gls_sigp_cb=vonKar**2/(cm0**2*beta2*cb_wallE)
     &                 *(rn**2-cff*rn/3.*(4.*rm+1.)
     &                    +cff**2*rm/9.*(2.+4.*rm))
# endif
!
!--------------------------------------------------------------------
! Horizontal advection (first-order upwind)
!--------------------------------------------------------------------
!
# ifdef GLS_HADV_UP1
#  define FX ustar_sfc_sq
#  define FE ustar_bot_sq
      ig=itke         ! <-- TKE
      DO k=1,N-1
        DO j=jmin,jmax
          DO i=imin,imax+1
            HUon_w = 0.5*(Huon(i,j,k)+Huon(i,j,k+1)) SWITCH umask(i,j)
            FX(i,j)=
     &            trb(i-1,j,k,nstp,ig)*max(HUon_w,0.)
     &           +trb(i  ,j,k,nstp,ig)*min(HUon_w,0.)
          ENDDO
        ENDDO
        DO j=jmin,jmax+1
          DO i=imin,imax
            HVom_w = 0.5*(Hvom(i,j,k)+Hvom(i,j,k+1)) SWITCH vmask(i,j)
            FE(i,j)=
     &            trb(i,j-1,k,nstp,ig)*max(HVom_w,0.)
     &           +trb(i,j  ,k,nstp,ig)*min(HVom_w,0.)
          ENDDO
        ENDDO
        DO j=jmin,jmax
          DO i=imin,imax
            cff    = 2. / ( Hz(i,j,k)+ Hz(i,j,k+1) )
            tke_new(i,j,k)= trb(i,j,k,nstp,ig)
     &                     -dt*cff*pm(i,j)*pn(i,j)
     &                       *(  FX(i+1,j)-FX(i,j)
     &                         + FE(i,j+1)-FE(i,j) )
          ENDDO
        ENDDO
      ENDDO
!
      ig=igls         ! <-- GLS
      DO k=1,N-1
        DO j=jmin,jmax
          DO i=imin,imax+1
            HUon_w = 0.5*(Huon(i,j,k)+Huon(i,j,k+1)) SWITCH umask(i,j)
            FX(i,j)=
     &            trb(i-1,j,k,nstp,ig)*max(HUon_w,0.)
     &           +trb(i  ,j,k,nstp,ig)*min(HUon_w,0.)
          ENDDO
        ENDDO
        DO j=jmin,jmax+1
          DO i=imin,imax
            HVom_w = 0.5*(Hvom(i,j,k)+Hvom(i,j,k+1)) SWITCH vmask(i,j)
            FE(i,j)=
     &            trb(i,j-1,k,nstp,ig)*max(HVom_w,0.)
     &           +trb(i,j  ,k,nstp,ig)*min(HVom_w,0.)
          ENDDO
        ENDDO
        DO j=jmin,jmax
          DO i=imin,imax
            cff    = 2. / ( Hz(i,j,k)+ Hz(i,j,k+1) )
            gls_new(i,j,k)= trb(i,j,k,nstp,ig)
     &                     -dt*cff*pm(i,j)*pn(i,j)
     &                       *(  FX(i+1,j)-FX(i,j)
     &                         + FE(i,j+1)-FE(i,j) )
          ENDDO
        ENDDO
      ENDDO
#  undef FX
#  undef FE
# else
      DO k=1,N-1
        DO j=jmin,jmax
          DO i=imin,imax
            tke_new(i,j,k)=trb(i,j,k,nstp,itke)
            gls_new(i,j,k)=trb(i,j,k,nstp,igls)
          ENDDO
        ENDDO
      ENDDO
# endif /* GLS_HADV_UP1 */
!
!--------------------------------------------------------------------
! Vertical advection (first-order upwind)
!--------------------------------------------------------------------
!
# ifdef GLS_VADV_UP1
      ig=itke         ! <-- TKE
      DO j=jmin,jmax
        DO k=1,N
          DO i=imin,imax
            cff=0.5*(We(i,j,k)+We(i,j,k-1))
            FC(i,k)=trb(i,j,k-1,nstp,ig)*max(cff,0.)
     &             +trb(i,j,k  ,nstp,ig)*min(cff,0.)
          ENDDO
        ENDDO
        DO k=1,N-1
          DO i=imin,imax
            cff=2./(Hz(i,j,k)+Hz(i,j,k+1))
            tke_new(i,j,k)= tke_new(i,j,k)
     &                     -dt*cff*pm(i,j)*pn(i,j)
     &                        *(FC(i,k+1)-FC(i,k))
          ENDDO
        ENDDO
      ENDDO
!
      ig=igls         ! <-- GLS
      DO j=jmin,jmax
        DO k=1,N
          DO i=imin,imax
            cff=0.5*(We(i,j,k)+We(i,j,k-1))
            FC(i,k)=trb(i,j,k-1,nstp,ig)*max(cff,0.)
     &             +trb(i,j,k  ,nstp,ig)*min(cff,0.)
          ENDDO
        ENDDO
        DO k=1,N-1
          DO i=imin,imax
            cff=2./(Hz(i,j,k)+Hz(i,j,k+1))
            gls_new(i,j,k)= gls_new(i,j,k)
     &                     -dt*cff*pm(i,j)*pn(i,j)
     &                        *(FC(i,k+1)-FC(i,k))
          ENDDO
        ENDDO
      ENDDO
# endif /* GLS_VADV_UP1 */
!
!--------------------------------------------------------------------
!  Compute the vertical shear (or complete deformation rate D)
!--------------------------------------------------------------------
!
!    D^2 = 2 DijDij      -- Dij = 0.5*(dui/dxj + duj/dxi) --
!        = (dudz + dwdx)^2 + (dvdz + dwdy)^2 + (dvdx + dudy)^2
!          + 2 dudx^2 + 2 dvdy^2 + 2 dwdz^2
!
!  Also compute enstrophy to identify nearly potential flow regions 
!  (Larsen & Fuhrman, 2018) : 
!    V^2 = 2 VijVij      -- Vij = 0.5*(dui/dxj - duj/dxi) --
!        = (dudz - dwdx)^2 + (dvdz - dwdy)^2 + (dvdx - dudy)^2
!--------------------------------------------------------------------
!
      tind  = nrhs
      DO k=1,N-1
         DO j=jmin,jmax
            DO i=imin,imax
               cff = 1. / ( Hz( i, j, k ) + Hz( i, j, k+1 ) )
               du  = cff*( u(i, j, k+1,tind)+u(i+1, j, k+1,tind)     ! dudz
     &                    -u(i, j, k  ,tind)-u(i+1, j, k  ,tind))
               dv  = cff*( v(i, j, k+1,tind)+v(i, j+1, k+1,tind)     ! dvdz
     &                    -v(i, j, k  ,tind)-v(i, j+1, k  ,tind))

# ifdef GLS_SHEAR_3D
#  ifdef NBQ
               dw  = 0.5*pm(i,j)
     &                  *( wz(i+1,j,k,tind)-wz(i-1,j,k,tind))        ! dwdx
               shear2(i,j,k) = (du+dw)*(du+dw)
#   ifdef GLS_STRAIN_LIM
               vort2(i,j,k)  = (du-dw)*(du-dw)
#   endif
               dw  = 0.5*pn(i,j)
     &                  *( wz(i,j+1,k,tind)-wz(i,j-1,k,tind))        ! dwdy
               shear2(i,j,k) = shear2(i,j,k) + (dv+dw)*(dv+dw)
#   ifdef GLS_STRAIN_LIM
               vort2(i,j,k)  = vort2(i,j,k)  + (dv-dw)*(dv-dw)
#   endif
#  else
               shear2(i,j,k) = du*du + dv*dv
#   ifdef GLS_STRAIN_LIM
               vort2(i,j,k)  = shear2(i,j,k)
#   endif
#  endif
               du  = 0.5*pm(i,j)
     &                  *( u(i+1, j, k+1,tind)-u(i, j, k+1,tind)     ! dudx
     &                    +u(i+1, j, k  ,tind)-u(i, j, k  ,tind))
               dv  = 0.5*pn(i,j)
     &                  *( v(i, j+1, k+1,tind)-v(i, j, k+1,tind)     ! dvdy
     &                    +v(i, j+1, k  ,tind)-v(i, j, k  ,tind))
#  ifdef NBQ
               dw  = cff*( wz(i,j,k+1,tind)-wz(i,j,k-1,tind))        ! dwdz
#  else
               dw = 0.
#  endif
               shear2(i,j,k) = shear2(i,j,k) +
     &                         2.*(du*du + dv*dv + dw*dw)
               du  = 0.125*pn(i,j)*
     &                        (u(i,j+1,k  ,tind)+u(i+1,j+1,k  ,tind) ! dudy
     &                        -u(i,j-1,k  ,tind)-u(i+1,j-1,k  ,tind)
     &                        +u(i,j+1,k+1,tind)+u(i+1,j+1,k+1,tind)
     &                        -u(i,j-1,k+1,tind)-u(i+1,j-1,k+1,tind))
               dv  = 0.125*pm(i,j)*
     &                        (v(i+1,j,k  ,tind)+v(i+1,j+1,k  ,tind) ! dvdx
     &                        -v(i-1,j,k  ,tind)-v(i-1,j+1,k  ,tind)
     &                        +v(i+1,j,k+1,tind)+v(i+1,j+1,k+1,tind)
     &                        -v(i-1,j,k+1,tind)-v(i-1,j+1,k+1,tind))
               shear2(i,j,k) = (du+dv)*(du+dv) + shear2(i,j,k)
#  ifdef GLS_STRAIN_LIM
               vort2(i,j,k)  = vort2(i,j,k)  + (dv-du)*(dv-du)
#  endif
# else
               shear2(i,j,k) = du*du + dv*dv
#  ifdef GLS_STRAIN_LIM
               vort2(i,j,k)  = shear2(i,j,k)
#  endif
# endif
           ENDDO
         ENDDO
      ENDDO
!
!--------------------------------------------------------------------
!  Compute ustar squared at the surface and at the bottom
!--------------------------------------------------------------------
!
      DO j=jmin,jmax
         DO i=imin,imax 
# ifdef OA_COUPLING
            ustar_sfc_sq(i,j)=smstr(i,j)
# else
            su_r=0.5*(sustr(i,j)+sustr(i+1,j))
            sv_r=0.5*(svstr(i,j)+svstr(i,j+1))
            ustar_sfc_sq( i, j ) = sqrt(su_r**2+sv_r**2)
# endif
            ustar_bot_sq( i, j ) = 
     &                        sqrt( (0.5*(bustr(i,j)+bustr(i+1,j)))**2
     &                             +(0.5*(bvstr(i,j)+bvstr(i,j+1)))**2 )
         ENDDO
      ENDDO

      !--------------------------------------------------
      DO j=jmin,jmax         !<-- j-outer loop
      !--------------------------------------------------

         ! Compute the dissipation rate
         DO i=imin,imax
            DO k=1,N-1
               cff       = (cm0**e1) * ( trb( i,j,k,nstp,itke )**e2 )
     &                               * ( trb( i,j,k,nstp,igls )**e3 )
               diss(i,k) = MAX( cff , eps_min )
            ENDDO
         ENDDO !<-- terminate i-loop

         !--------------------------------------------------
         DO ig = 1,ngls     ! ig = 2 for gls and = 1 for tke
         !--------------------------------------------------

            ! Off-diagonal terms for the tridiagonal problem
            cff=-0.5*dt
            DO k=2,N-1
               DO i=imin,imax
                  FC(i,k)= cff* OneOverSig(ig)*
     &                   ( Akv_old(i,j,k)+Akv_old(i,j,k-1) ) / Hz(i,j,k)
               ENDDO
            ENDDO

            DO i=imin,imax
               FC(i,1)=0.
               FC(i,N)=0.
            ENDDO

            ! Production/Dissipation terms and diagonal term
            DO k=1,N-1
               DO i=imin,imax
                  ig1   = (igls-ig)
                  ig2   = (ig-itke)
                  invk  =     1. / trb( i,j,k,nstp,itke )
                  gls   =          trb( i,j,k,nstp,igls )
                  ! invG = 1 for tke invg=1/psi for gls
                  invG  =  ig1+ig2*(1./gls)
                  cff1  =  ig1+ig2*beta1   * invk*gls
                  cff2  = (ig1+ig2*beta2 ) * invk
                  cff3m =  ig1+ig2*beta3m  * invk*gls
                  cff3p =  ig1+ig2*beta3p  * invk*gls
                  ! Shear and buoyancy production
                  Sprod =  cff1*Akv_old(i,j,k) * shear2(i,j,k)    
# if defined TEMPERATURE || defined SALINITY
                  Bprod = -Akt_old(i,j,k)*( cff3m*MAX(bvf(i,j,k),0.)
     &                                   +  cff3p*MIN(bvf(i,j,k),0.) )
# else
                  Bprod = 0.
# endif
                  ! Patankar trick to ensure non-negative solutions
                  cff   =       0.5*(Hz(i,j,k)+Hz(i,j,k+1))
                  IF( ig == itke ) THEN
                    trb_new=tke_new(i,j,k)
                  ELSE
                    trb_new=gls_new(i,j,k)
                  ENDIF
                  IF( (Bprod + Sprod) .gt. 0.) THEN
                     RH(i,k) = cff*( trb_new + dt*(Bprod+Sprod) )
                     DC(i,k) = cff*(1.+dt*cff2*diss(i,k))
     &                                                -FC(i,k)-FC(i,k+1)
                  ELSE
                     RH(i,k) = cff*( trb_new + dt*Sprod  )
                     DC(i,k) = cff*(1.+dt*(cff2*diss(i,k)
     &                              -invG*Bprod)) - FC(i,k) - FC(i,k+1)
                  ENDIF
               ENDDO
            ENDDO

            ! Boundary conditions
            IF( ig == itke ) THEN
               DO i=imin,imax
                  ! surface
# ifdef MRL_WCI
                  trb_sfc = MAX( tke_min, cmu_fac*
     &                                       (zalp*wepb(i,j))**(2./3.) )
                  flux_top=zalp*wepb(i,j)
# else
                  trb_sfc = MAX( tke_min, cm0inv2*ustar_sfc_sq(i,j) )
                  flux_top = 0.
# endif
                  ! bottom
                  trb_bot = MAX( tke_min, cm0inv2*ustar_bot_sq(i,j) )
                  flux_bot = 0.
                  ! finalize
                  RH(i,1  ) = RH(i,  1) + dt*flux_bot
                  RH(i,N-1) = RH(i,N-1) + dt*flux_top
                  tke_new(i,j,N) = trb_sfc
                  tke_new(i,j,0) = trb_bot
               ENDDO
            ELSE
               DO i=imin,imax
                  ! surface
# ifdef MRL_WCI
                  z0_s=MAX( Zosmin , 0.5*whrm(i,j) )
# else
                  z0_s=MAX( Zosmin , chk*ustar_sfc_sq(i,j) )   !<-- Charnock
# endif
                  cff = 0.5*(tke_new(i,j,N-1)+tke_new(i,j,N))
                  lgthsc = vonKar*(0.5*Hz(i,j,N)+z0_s)
                  trb_sfc = MAX(gls_min,(cm0**rp)*(lgthsc**rn)
     &                                              *(cff**rm))
                  flux_top = -rn*cm0**(rp+1.)
     &                             *vonKar*OneOverSig(igls)
     &                             *(cff**(rm+0.5))*(lgthsc**rn)
# ifdef MRL_WCI
                  flux_top=flux_top - OneOverSig(itke)*OneOverSig(igls)
     &                                  *   cm0**(rp+1.)
     &                                  *   rm
     &                                  *   cff**(rm-0.5)
     &                                  *   lgthsc** (rn+1.)
     &                                  *   zalp*   wepb(i,j)
# endif
                  ! bottom
#  ifdef BBL
                  z0_b = MAX( Zbapp(i,j) , Zobmin )
#  else
                  z0_b = MAX( Zob(i,j) , Zobmin )
#  endif
                  cff = 0.5*(tke_new(i,j,1)+tke_new(i,j,0))
                  lgthsc = vonKar*(0.5*Hz(i,j,1)+z0_b)
                  trb_bot = MAX(gls_min,(cm0**rp)*(lgthsc**rn)
     &                                   *(tke_new(i,j,0)**rm))
                  flux_bot =-rn*cm0**(rp+1.)
     &                          *vonKar*OneOverSig(igls)
     &                          *(cff**(rm+0.5))*(lgthsc**rn)
                  ! finalize
                  RH(i,  1) = RH(i,  1) + dt*flux_bot
                  RH(i,N-1) = RH(i,N-1) + dt*flux_top
                  gls_new(i,j,N) = trb_sfc
                  gls_new(i,j,0) = trb_bot
               ENDDO
            ENDIF

            ! tridiagonal resolution
            DO i=imin,imax
               cff       =  1./DC(i,N-1)
               CF(i,N-1) = cff*FC(i,N-1)
               RH(i,N-1) = cff*RH(i,N-1)
            ENDDO
            DO k=N-2,1,-1
               DO i=imin,imax
                  cff     =   1./(DC(i,k)-CF(i,k+1)*FC(i,k+1))
                  CF(i,k) = cff*FC(i,k)
                  RH(i,k) = cff*( RH(i,k)-FC(i,k+1)*RH(i,k+1))
               ENDDO
            ENDDO
            IF( ig == itke ) THEN
              DO i=imin,imax
                tke_new(i,j,1) = MAX( RH(i,1), trb_min(ig) )
              ENDDO
              DO k=2,N-1
                DO i=imin,imax
                  RH(i,k) = RH(i,k)-CF(i,k)*RH(i,k-1)
                  tke_new(i,j,k) = MAX( RH(i,k), trb_min(ig) )
                ENDDO
              ENDDO
            ELSE
              DO i=imin,imax
                gls_new(i,j,1) = MAX( RH(i,1), trb_min(ig) )
              ENDDO
              DO k=2,N-1
                DO i=imin,imax
                  RH(i,k) = RH(i,k)-CF(i,k)*RH(i,k-1)
                  gls_new( i,j,k) = MAX( RH(i,k), trb_min(ig) )
                ENDDO
              ENDDO
            ENDIF

         !--------------------------------------------------
         ENDDO     ! ig loop
         !--------------------------------------------------

      !--------------------------------------------------
      ENDDO     ! j loop
      !--------------------------------------------------

!
!--------------------------------------------------------------------
!    GLS Smoothing
!--------------------------------------------------------------------
!
# ifdef PSI_SMOOTH
      ig = igls
#  define TRB_NEW gls_new
#  define FX ustar_sfc_sq
#  define FE ustar_bot_sq
#  define FE1 work
      DO k=0,N
#        include "gls_smooth.h"
      ENDDO
#  undef FX
#  undef FE
#  undef FE1
#  undef TRB_NEW
# endif

# ifdef TKE_SMOOTH
        ig = itke
#  define TRB_NEW tke_new
#  define FX ustar_sfc_sq
#  define FE ustar_bot_sq
#  define FE1 work
      DO k=0,N
#        include "gls_smooth.h"
      ENDDO
#  undef FX
#  undef FE
#  undef FE1
#  undef TRB_NEW
# endif

# undef I_EXT_RANGE
# undef J_EXT_RANGE

!
!--------------------------------------------------------------------
!    Compute Akv & Lscale
!--------------------------------------------------------------------
!
      !--------------------------------------------------
      DO j=jstr,jend  !<-- j-outer loop
      !--------------------------------------------------
         DO k=1,N-1
            DO i=istr,iend
               !
               ! Galperin limitation : l <= l_lim
               L_lim = galp * sqrt( 2.* trb(i,j,k,nnew,itke)) /
     &                            ( sqrt(max(eps, bvf(i,j,k)))  )
               !
               ! Limitation on psi (use MAX because rn is negative)
               cff = (cm0**rp)*(L_lim**rn)*(trb(i,j,k,nnew,itke)**rm)
               trb( i,j,k,nnew,igls ) = MAX( trb( i,j,k,nnew,igls ),cff)
               !
               ! Dissipation rate epsilon
               epsilon = (cm0**e1) * ( trb( i,j,k,nnew,itke )**e2 )
     &                             * ( trb( i,j,k,nnew,igls )**e3 )
               !
               ! Limitation of epsilon
# ifdef GLS_STRAIN_LIM
               ! increase diss in potential flow (Larsen & Fuhrman, 2018)
               cff = lars * (shear2(i,j,k)/MAX(vort2(i,j,k),1.e-10))**2
               epsilon = MAX(epsilon, cff*epsilon)
# endif
               epsilon = MAX(epsilon, eps_min)
               !
               ! Compute alpha_n and alpha_m
               cff     = ( trb(i,j,k,nnew,itke)/epsilon )**2
               alpha_m     = cff*  shear2(i,j,k)
               alpha_n     = cff*     bvf(i,j,k)
               !
               ! Limitation of alpha_n and alpha_m
               alpha_n     = MIN(  MAX( 0.73*alpha_n_min , alpha_n ) , 1.0e10 )
               alpha_m_max = ( lim_am0 + lim_am1 * alpha_n
     &                                 + lim_am2 * alpha_n**2
     &                                 + lim_am3 * alpha_n**3) /
     &                         ( lim_am4 + lim_am5 * alpha_n
     &                                   + lim_am6 * alpha_n**2 )
               alpha_m = MIN(alpha_m , alpha_m_max)
               !
               ! Compute stability functions
               Denom = sf_d0  + sf_d1*alpha_n +  sf_d2*alpha_m
     &                        + sf_d3*alpha_n*alpha_m
     &                        + sf_d4*alpha_n**2 + sf_d5*alpha_m**2
               c_mu      = (sf_n0  + sf_n1*alpha_n  + sf_n2*alpha_m)
     &                                                        /Denom
               c_mu_prim = (sf_nb0 + sf_nb1*alpha_n + sf_nb2*alpha_m)
     &                                                        /Denom
               !
               ! Finalize the computation of Akv and Akt
               cff = trb( i,j,k,nnew,itke )**2 / epsilon
               Akv(i,j,k) = MAX( cff*c_mu , nuwm ) SWITCH rmask(i,j)
# ifdef TEMPERATURE
               Akt(i,j,k,itemp)= MAX( cff*c_mu_prim , nuws )
     &                                             SWITCH rmask(i,j)
# endif
# ifdef SALINITY
               Akt(i,j,k,isalt)= MAX( cff*c_mu_prim , nuws )
     &                                             SWITCH rmask(i,j)
# endif
               !
               ! Computate and store Lscale and epsilon
               Lscale( i, j , k ) =  cm0 * cm0 * cm0 * cff
     &                              / sqrt( trb( i,j,k,nnew,itke ) )
     &                                             SWITCH rmask(i,j)
               Eps_gls(i,j,k) = epsilon SWITCH rmask(i,j)
            ENDDO
         ENDDO

         DO i=istr,iend
           Akv(i,j,N) = MAX( 1.5*Akv(i,j,N-1)
     &                      -0.5*Akv(i,j,N-2), nuwm)
           Akv(i,j,0) = MAX( 1.5*Akv(i,j,  1)
     &                      -0.5*Akv(i,j,  2), nuwm)
# ifdef TEMPERATURE
           Akt(i,j,N,itemp) = MAX(  1.5*Akt(i,j,N-1,itemp)
     &                             -0.5*Akt(i,j,N-2,itemp), nuws )
           Akt(i,j,0,itemp) = MAX(  1.5*Akt(i,j,  1,itemp)
     &                             -0.5*Akt(i,j,  2,itemp), nuws )
# endif
# ifdef SALINITY
           Akt(i,j,N,isalt)=Akt(i,j,N,itemp)
           Akt(i,j,0,isalt)=Akt(i,j,0,itemp)
# endif
         ENDDO
      !--------------------------------------------------
      ENDDO  !<-- end j-outer loop
      !--------------------------------------------------

!
!------------------------------------------------
!  Apply boundary conditions
!------------------------------------------------
!
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstr,jend
          do k=0,N
            trb(istr-1,j,k,nnew,itke)=trb(istr,j,k,nnew,itke)
            trb(istr-1,j,k,nnew,igls)=trb(istr,j,k,nnew,igls)
            Akv(istr-1,j,k      )=Akv(istr,j,k      )
#  ifdef TEMPERATURE
            Akt(istr-1,j,k,itemp)=Akt(istr,j,k,itemp)
#  endif
#  ifdef SALINITY
            Akt(istr-1,j,k,isalt)=Akt(istr,j,k,isalt)
#  endif
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr,jend
          do k=0,N
            trb(iend+1,j,k,nnew,itke)=trb(iend,j,k,nnew,itke)
            trb(iend+1,j,k,nnew,igls)=trb(iend,j,k,nnew,igls)
            Akv(iend+1,j,k      )=Akv(iend,j,k      )
#  ifdef TEMPERATURE
            Akt(iend+1,j,k,itemp)=Akt(iend,j,k,itemp)
#  endif
#  ifdef SALINITY
            Akt(iend+1,j,k,isalt)=Akt(iend,j,k,isalt)
#  endif
          enddo
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          do k=0,N
            trb(i,jstr-1,k,nnew,itke)=trb(i,jstr,k,nnew,itke)
            trb(i,jstr-1,k,nnew,igls)=trb(i,jstr,k,nnew,igls)
            Akv(i,jstr-1,k      )=Akv(i,jstr,k      )
#  ifdef TEMPERATURE
            Akt(i,jstr-1,k,itemp)=Akt(i,jstr,k,itemp)
#  endif
#  ifdef SALINITY
            Akt(i,jstr-1,k,isalt)=Akt(i,jstr,k,isalt)
#  endif
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend
          do k=0,N
            trb(i,jend+1,k,nnew,itke)=trb(i,jend,k,nnew,itke)
            trb(i,jend+1,k,nnew,igls)=trb(i,jend,k,nnew,igls)
            Akv(i,jend+1,k)=Akv(i,jend,k)
#  ifdef TEMPERATURE
            Akt(i,jend+1,k,itemp)=Akt(i,jend,k,itemp)
#  endif
#  ifdef SALINITY
            Akt(i,jend+1,k,isalt)=Akt(i,jend,k,isalt)
#  endif
          enddo
        enddo
      endif
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=0,N
          trb(istr-1,jstr-1,k,nnew,itke)=trb(istr,jstr,k,nnew,itke)
          trb(istr-1,jstr-1,k,nnew,igls)=trb(istr,jstr,k,nnew,igls)
          Akv(istr-1,jstr-1,k      )=Akv(istr,jstr,k      )
#  ifdef TEMPERATURE
          Akt(istr-1,jstr-1,k,itemp)=Akt(istr,jstr,k,itemp)
#  endif
#  ifdef SALINITY
          Akt(istr-1,jstr-1,k,isalt)=Akt(istr,jstr,k,isalt)
#  endif
        enddo
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=0,N
          trb(istr-1,jend+1,k,nnew,itke)=trb(istr,jend,k,nnew,itke)
          trb(istr-1,jend+1,k,nnew,igls)=trb(istr,jend,k,nnew,igls)
          Akv(istr-1,jend+1,k      )=Akv(istr,jend,k      )
#  ifdef TEMPERATURE
          Akt(istr-1,jend+1,k,itemp)=Akt(istr,jend,k,itemp)
#  endif
#  ifdef SALINITY
          Akt(istr-1,jend+1,k,isalt)=Akt(istr,jend,k,isalt)
#  endif
        enddo
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=0,N
          trb(iend+1,jstr-1,k,nnew,itke)=trb(iend,jstr,k,nnew,itke)
          trb(iend+1,jstr-1,k,nnew,igls)=trb(iend,jstr,k,nnew,igls)
          Akv(iend+1,jstr-1,k      )=Akv(iend,jstr,k      )
#  ifdef TEMPERATURE
          Akt(iend+1,jstr-1,k,itemp)=Akt(iend,jstr,k,itemp)
#  endif
#  ifdef SALINITY
          Akt(iend+1,jstr-1,k,isalt)=Akt(iend,jstr,k,isalt)
#  endif
        enddo
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=0,N
          trb(iend+1,jend+1,k,nnew,itke)=trb(iend,jend,k,nnew,itke)
          trb(iend+1,jend+1,k,nnew,igls)=trb(iend,jend,k,nnew,igls)
          Akv(iend+1,jend+1,k)=Akv(iend,jend,k)
#  ifdef TEMPERATURE
          Akt(iend+1,jend+1,k,itemp)=Akt(iend,jend,k,itemp)
#  endif
#  ifdef SALINITY
          Akt(iend+1,jend+1,k,isalt)=Akt(iend,jend,k,isalt)
#  endif
        enddo
      endif
#  endif
# endif

!
!  Compute mixed layer depth
!
      do j=Jstr,Jend
        do i=Istr,Iend
          kbl(i,j)=1
          hbl(i,j)=z_w(i,j,N)-z_r(i,j,1)
# ifdef MASKING
          hbl(i,j)=hbl(i,j)*rmask(i,j)
# endif
        enddo
      enddo
      kref=max(1,N-3)
      do k=1,kref
        do j=Jstr,Jend
          do i=Istr,Iend
            cff=rho1(i,j,k)-rho1(i,j,kref)
            if (cff.gt.0.01) then
              kbl(i,j)=k
              hbl(i,j)=z_w(i,j,N)-z_r(i,j,k)
# ifdef MASKING
              hbl(i,j)=hbl(i,j)*rmask(i,j)
# endif
            endif
          enddo
        enddo
      enddo
!
! Exchange computational margines and/or periodic boundaries:
!--------- ------------- -------- ------ -------- -----------
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        trb(START_2D_ARRAY,0,nnew,itke))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        trb(START_2D_ARRAY,0,nnew,igls))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend, Akv)
#  ifdef TEMPERATURE
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        Akt(START_2D_ARRAY,0,itemp))
#  endif
#  ifdef SALINITY
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        Akt(START_2D_ARRAY,0,isalt))
#  endif
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        Lscale(START_2D_ARRAY,0))
!      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
!     &                        Eps_gls(START_2D_ARRAY,0))
# endif
!======================================================================
      return
      end
!======================================================================
#endif  /* SOLVE3D & GLS_MIXING */
