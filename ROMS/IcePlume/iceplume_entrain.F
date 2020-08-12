#include "cppdefs.h"
#ifdef ICEPLUME
!
! ==================================================================!
!                                                                   !
! These are the core functions of iceplume entrainment model.       !
!                                                                   !
! ==================================================================!
!
SUBROUTINE ICEPLUME_ENTRAIN(ng, I, iceDepthK,                           &
     &                      fIni, tIni, sIni,                           &
     &                      sgTyp, sgLen)
!
  USE mod_iceplume
  implicit none
!
! ==================================================================!
!                                                                   !
! Input variables:                                                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! ng          - grid identifier                                     !
! I           - point source index identifier                       !
! iceDepthK   - vertical grid index of grounding line depth         !
! fIni        - subglacial runoff volume flux [m^3 s^-1]            !
! tIni        - subglacial runoff temperature [degC]                !
! sIni        - subglacial runoff salinity [PSU]                    !
! sgTyp       - subglacial runoff type identifier                   !
! sgLen       - subglacial runoff discharge length [m]              !
!                                                                   !
! ==================================================================!
!                                                                   !
! Local variables:                                                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! RHO         - function to compute density                         !
! rhoA/tA/sA  - ambient density/temperature/salinity                !
! zIn/zOut    - in/out depth of the model steps                     !
!                                                                   !
! ==================================================================!
!
! In/out variables
!
  integer, intent(in) :: ng, I, iceDepthK
  real(r8), intent(in) :: fIni, tIni, sIni
  integer, intent(in) :: sgTyp
  real(r8), intent(in) :: sgLen
!
! Local variables declaration
!
  real(r8) :: RHO
  real(r8) :: rhoA, tA, sA
  real(r8) :: zIn, zOut
  integer :: K
!
! Local variables for ODEPACK
!
  integer :: IOPT, ISTATE, ITASK, ITOL, IWORK(20), LIW, LRW, MF, NEQ
  real(r8), parameter :: RTOL = 1.0E-5_r8
  real(r8), parameter :: ATOL = 1.0E-5_r8
  real(r8) :: RWORK(116), Y(10)
  external GENERAL_ENTRAIN_MODEL, JEX
!
! ==================================================================!
!                                                                   !
! Y is input/output vector for DLSODE                               !
!   Y(1)  - plume volume flux                                       !
!   Y(2)  - plume velocity                                          !
!   Y(3)  - plume temperature                                       !
!   Y(4)  - plume salinity                                          !
!   Y(5)  - integrated plume/glacier contact area                   !
!   Y(6)  - area integrated melt                                    !
!   Y(7)  - fake variable to pass in ng                             !
!   Y(8)  - fake variable to pass in I                              !
!   Y(9)  - fake variable to pass in sgTyp                          !
!   Y(10) - fake variable to pass in sgLen                          !
!                                                                   !
! ==================================================================!
!                                                                   !
! For ODEPACK solver. See ODEPACK documentation and source code in  !
! Cowton et al. 2015.                                               !
!                                                                   !
! ==================================================================!
!
  NEQ      = 6
  LRW      = 116
  LIW      = 116
!
  ITOL     = 1
  ITASK    = 1
  ISTATE   = 1
  IOPT     = 0
  MF       = 10
  IWORK(7) = 2  ! To limit # of repeat error messages are printed
!
! Initial conditions
!
  Y(1)  = fIni  ! initial plume volume flux
  Y(2)  = wIni  ! initial vertical velocity
  Y(3)  = tIni  ! initial temperature
  Y(4)  = sIni  ! initial salinity
  Y(5)  = 0.0   ! integrated contact area
  Y(6)  = 0.0   ! integrated melt rate
!
  Y(7)  = REAL(ng)     ! fake variable to pass in ng
  Y(8)  = REAL(I)      ! fake variable to pass in I
  Y(9)  = REAL(sgTyp)  ! fake variable to pass in sgTyp
  Y(10) = sgLen        ! fake variable to pass in sgLen
!
! Start at bottom of ice face
!
  zIn = PLUME(ng) % zW(I, iceDepthK)
  zOut = PLUME(ng) % zW(I, iceDepthK+1)
!
! Set initial conditions
!
  PLUME(ng) % f(I, iceDepthK)    = Y(1)
  PLUME(ng) % w(I, iceDepthK)    = Y(2)
  PLUME(ng) % t(I, iceDepthK)    = Y(3)
  PLUME(ng) % s(I, iceDepthK)    = Y(4)
  PLUME(ng) % a(I, iceDepthK)    = Y(5)
  PLUME(ng) % mInt(I, iceDepthK) = Y(6)
  PLUME(ng) % rho(I, iceDepthK)  = RHO(Y(3), Y(4), zIn)
!
  CALL PLUME_METRICS (sgTyp, Y(1)/Y(2), sgLen,                          &
     &                PLUME(ng) % lm(I, iceDepthK),                     &
     &                PLUME(ng) % lc(I, iceDepthK))
!
! Move up through water column from lowest layer
!
  DO K = iceDepthK+1, N(ng)
!
! ==================================================================!
!                                                                   !
! Use DLSODE to solve plume properties.                             !
!                                                                   !
! ==================================================================!
!
! Check to make sure plume has not reached neutral buoyancy in a lower
! layer
!
    IF (ISTATE .GT. -1) THEN
      CALL DLSODE (GENERAL_ENTRAIN_MODEL,                               &
     &             NEQ, Y, zIn, zOut,                                   &
     &             ITOL, RTOL, ATOL, ITASK,                             &
     &             ISTATE, IOPT, RWORK, LRW, IWORK,                     &
     &             LIW, JEX, MF)
!
      CALL PLUME_METRICS (sgTyp, Y(1)/Y(2), sgLen,                      &
     &                    PLUME(ng) % lm(I, K),                         &
     &                    PLUME(ng) % lc(I, K))
!
! ==================================================================!
!                                                                   !
! Test to see if neutral buoyancy has now been reached. If solver   !
! returns ISTATE = -1, then it has been unable to meet required     !
! tolerances at this level. This generally occurs because plume     !
! has reached neutral buoyancy and run out of momentum, and so is   !
! no longer rising. At this point, we therefore end the call to the !
! plume model. Our aim is to catch the plume at the point of        !
! neutral buoyancy. We therefore perform a manual comparrison of    !
! ambient and plume density. If plume density >= ambient density we !
! assign ISTATE = -1, again ending the call to the plume model.     !
!                                                                   !
! ==================================================================!
!
! Calculate plume density (rho = RHO(temp, salt, depth))
!
      PLUME(ng) % rho(I, K) = RHO(Y(3), Y(4), zIn)
# ifdef ICEPLUME_SHOOT_NO
!
! Calculate ambient density
!
      IF (K .EQ. 0) THEN
        tA = PLUME(ng) % tAm(I, 1)
        sA = PLUME(ng) % sAm(I, 1)
      ELSEIF (K .EQ. N(ng)) THEN
        tA = PLUME(ng) % tAm(I, N(ng))
        sA = PLUME(ng) % sAm(I, N(ng))
      ELSE
        tA = 0.5*(PLUME(ng) % tAm(I, K) + PLUME(ng) % tAm(I, K+1))
        sA = 0.5*(PLUME(ng) % sAm(I, K) + PLUME(ng) % sAm(I, K+1))
      ENDIF
      rhoA = RHO(tA, sA, zIn)
      IF ((PLUME(ng) % rho(I, K) .GT. rhoA) .OR. (K .EQ. N(ng))) THEN
# else
!
      IF ((Y(2) .LE. 0.0) .OR. (K .EQ. N(ng))) THEN
# endif
        ISTATE = -1
      ENDIF
!
! If the plume has not reached neutral buoyancy, then we assign a depth
! at which to calculate the next value and loop round to call the plume
! model again. Make sure we are not at the surface
!
      IF (ISTATE .GT. -1) THEN
        zIn = zOut
        zOut = PLUME(ng) % zW(I, K+1)
      ENDIF
!
    ELSE  ! (ISTATE .LE. -1)
!
! This section is entered once the plume has reached neutral buoyancy
! once plume has reached neutral buoyancy, no plume values
!
      Y(1) = 0.0
      Y(2) = 0.0
      Y(3) = 0.0
      Y(4) = 0.0
      Y(5) = 0.0
      Y(6) = 0.0
    ENDIF
!
! Save results.
!
    PLUME(ng) % f(I, K) = Y(1)
    PLUME(ng) % w(I, K) = Y(2)
    PLUME(ng) % t(I, K) = Y(3)
    PLUME(ng) % s(I, K) = Y(4)
    PLUME(ng) % a(I, K) = Y(5)
    PLUME(ng) % mInt(I, K) = Y(6)
!
  ENDDO
  RETURN
END SUBROUTINE ICEPLUME_ENTRAIN
!
! ==================================================================!
!                                                                   !
! Below is the generalized plume model.                             !
!                                                                   !
! ==================================================================!
!
SUBROUTINE GENERAL_ENTRAIN_MODEL (NEQ, T, Y, YDOT)
!
  USE mod_iceplume
  implicit none
!
  integer, intent(in) :: NEQ
  integer :: ng, I, sgTyp
  real(r8) :: sgLen
  real(r8) :: T, Y(10), YDOT(6)
  real(r8) :: RHO
  real(r8) :: tA, sA, tB, sB, rhoP, rhoA, gRed
  real(r8) :: mdot, a, b, c
  real(r8) :: Lm, Lc
!
! Check if velcoity is positive
!
  IF (Y(2) .LE. 0.0_r8) THEN
    YDOT(1) = 0.0_r8
    YDOT(2) = 0.0_r8
    YDOT(3) = 0.0_r8
    YDOT(4) = 0.0_r8
    YDOT(5) = 0.0_r8
    YDOT(6) = 0.0_r8
  ELSE
!
! Interpolate from imposed ambient profiles
!
    ng    = NINT(Y(7))
    I     = NINT(Y(8))
    sgTyp = NINT(Y(9))
    sgLen = Y(10)
!
    IF (T .LE. PLUME(ng) % zR(I, 1)) THEN
      tA = PLUME(ng) % tAm(I, 1)
      sA = PLUME(ng) % sAm(I, 1)
    ELSEIF (T .GE. PLUME(ng) % zR(I, N(ng))) THEN
      tA = PLUME(ng) % tAm(I, N(ng))
      sA = PLUME(ng) % sAm(I, N(ng))
    ELSE
      CALL LININT(N(ng), PLUME(ng) % zR(I, :),                          &
     &            PLUME(ng) % tAm(I, :), T, tA)
      CALL LININT(N(ng), PLUME(ng) % zR(I, :),                          &
     &            PLUME(ng) % sAm(I, :), T, sA)
    ENDIF
!
! Calculate reduced gravity
!
    rhoP = RHO(Y(3), Y(4), T)
    rhoA = RHO(tA, sA, T)
    gRed  = g*(rhoA-rhoP)/rhoRef
!
! Calculate plume metrics
!
    CALL PLUME_METRICS (sgTyp, Y(1)/Y(2), sgLen, Lm, Lc)
!
! Equations for sB, tB and mdot
!
    a = lambda1*(GamT*cW - GamS*cI)
    b = GamS*cI*(lambda1*Y(4) - lambda2 - lambda3*T +                   &
     &           tIce - (L/cI)) -                                       &
     &           GamT*cW*(Y(3) - lambda2 - lambda3*T + sIce)
    c = GamS*Y(4)*(cI*(lambda2 + lambda3*T - tIce) + L) +               &
     &  GamT*sIce*cW*(Y(3) - lambda2 - lambda3*T)
!
    sB   = (1.0/(2.0*a))*(-b - SQRT(b**2 - 4.0*a*c))
    tB   = lambda1*sB + lambda2 + lambda3*T
    mdot = GamS*SQRT(Cd)*Y(2)*(Y(4) - sB)/sB
!
! Plume volume flux
!
    YDOT(1) = alpha*Lc*Y(2) + Lm*mdot
!
! Plume vertical velocity
!
    YDOT(2) = (1./Y(1))*(-Y(2)*YDOT(1) + gRed*Y(1)/Y(2) -               &
     &                   Cd*Lm*Y(2)**2)
!
! Plume temperature
!
    YDOT(3) = (1./Y(1))*(-Y(3)*YDOT(1) + alpha*Lc*Ta*Y(2) +             &
     &                   Lm*mdot*tB -                                   &
     &                   SQRT(Cd)*GamT*Lm*Y(2)*(Y(3)-tB))
!
! Plume salinity
!
    YDOT(4) = (1./Y(1))*(-Y(4)*YDOT(1) + alpha*Lc*Sa*Y(2) +             &
     &                   Lm*mdot*sB -                                   &
     &                   SQRT(Cd)*GamS*Lm*Y(2)*(Y(4)-sB))
!
! Along-plume integrated contact area and melt rate
!
    YDOT(5) = Lm
    YDOT(6) = Lm*mdot
  ENDIF
!
END SUBROUTINE GENERAL_ENTRAIN_MODEL
!
! ==================================================================!
!                                                                   !
! Use this function to calculate plume metrics.                     !
!                                                                   !
! ==================================================================!
!
SUBROUTINE PLUME_METRICS (sgTyp, area, Ls, Lm, Lc)
!
! Calculate Lm (plume contact length with ice) and C (plume contact
! length with water), based on various of plume metrics.
!
! For each method, refer to the manual for details.
!
! If you wish to add new plume metrics, assign a new sgTyp number
! and code it to a new CASE.
!
  USE mod_iceplume
  implicit none
!
  integer, intent(in) :: sgTyp
  real(r8), intent(in) :: area, Ls
  real(r8), intent(out) :: Lm, Lc
  real(r8) :: gam, a, b, h
!
  SELECT CASE (sgTyp)
    CASE (2)
!
! Halfcone (Cowton)
!
      Lc = SQRT(2.*pi*area)
      Lm = 2./pi*Lc
!
    CASE (3)
!
! Finite line source
!
      Lc = SQRT(Ls**2 + 2.*pi*area)
      Lm = Ls + 2./pi*(Lc - Ls)
!
    CASE (4)
!
! Sheet (Jenkins)
!
      Lc = Ls
      Lm = Ls
!
    CASE (5)
!
! Fullcone
!
      Lc = 2.*SQRT(2.*pi*area)
      Lm = 0.
!
    CASE (6)
!
! Elipse (Aspect ratio 0.5)
!
      gam = 0.5_r8
      b = SQRT(2.*area/(gam*pi))
      a = gam*b
      h = ((a-b)/(a+b))**2
      Lc = 0.5*pi*(a+b)*(1.+3.*h/(10.+SQRT(4.-3.*h)))
      Lm = 2.*b
!
    CASE DEFAULT
!
! Halfcone (Cowton)
!
      Lc = SQRT(2.*pi*area)
      Lm = 2./pi*Lc
!
  END SELECT
!
END SUBROUTINE PLUME_METRICS
!
! ==================================================================!
!                                                                   !
! Use this function to calculate seawater density.                  !
!                                                                   !
! ==================================================================!
!
DOUBLE PRECISION FUNCTION RHO(T,S,z)
!
! Equation of state (UNESCO 1983)
!     T = temperature (deg C)
!     S = salinity (PSU)
!     z = depth (m)
!
  DOUBLE PRECISION T,S,z
  DOUBLE PRECISION rho_0, g, P
  DOUBLE PRECISION kw, Aw, Bw, k0
  DOUBLE PRECISION bulk_modulus
  DOUBLE PRECISION A, B, C, rho_w,rho_zero
!
  PARAMETER(rho_0=1027)
  PARAMETER(g=9.81)
!
  P= rho_0*g*abs(z)*1.0E-5
!
! RHO_1 (in situ)
  kw= 19652.21+ 148.4206*T- 2.327105*T**2+ &
   &    1.360477e-2*(T**3)-5.155288e-5*(T**4)
  Aw= 3.239908+ 1.43713e-3*T+ 1.16092e-4*T**2- &
   &    5.77905e-7*T**3
  Bw= 8.50935e-5- 6.12293e-6*T + 5.2787e-8*(T**2)
  k0= kw + (54.6746- 0.603459*T+ 1.09987e-2*(T**2) &
   &    -6.1670e-5*(T**3))*S +(7.944e-2 + 1.6483e-2* &
   &    T- 5.3009e-4*(T**2))*(S**1.5)
  A=  Aw+ (2.2838e-3- 1.0981e-5*T- 1.6078e-6*(T**2)) &
   &    *S+ 1.91075e-4*(S**1.5)
  B= Bw+ (-9.9348e-7+ 2.0816e-8*T+ 9.1697e-10*T**2)*S
  bulk_modulus= k0+ A*P+ B*P**2
!
  A= 8.24493e-1- 4.0899e-3*T+ 7.6438e-5*T**2- &
   &   8.2467e-7*T**3+5.3875e-9*T**4
  B= -5.72466e-3 + 1.0227e-4*T- 1.6546e-6*T**2
  C= 4.8314e-4
  rho_w= 999.842594 + 6.793952e-2*T- 9.095290e-3*T**2+ &
   &       1.001685e-4*T**3-1.120083e-6*T**4+ &
   &       6.536336e-9*T**5
  rho_zero= rho_w+ A*S + B*(S**1.5)+ C*(S**2)
!
  RHO= rho_zero/(1- (P/bulk_modulus))
!
END
!
! ==================================================================!
!                                                                   !
! Use this function to do 1D interpolation.                         !
!                                                                   !
! ==================================================================!
!
SUBROUTINE LININT(nx,xtab,ytab,x,y)
!
! Given a value of x return a value of y based on interpolation
! within a table of y values (ytab) corresponding to the x values
! contained in the array xtab.  The subroutine assumes that the
! values in xtab increase monotonically
!
! John Mahaffy 2/12/95
! Modified slightly TRC 2014
!
  integer nx
  double precision xtab(nx), ytab(nx), x, y
!
! local variables
!
  integer i, i1
  double precision  wx
!
  if (x.lt.(xtab(1)).or.x.GT.(xtab(nx))) then
    write(6, *) 'x = ', x, '  is out of table range'
    stop
  endif
  do 100 i=2,nx
       if (x.le.xtab(i)) go to 200
  100 continue
  200 i1=i-1
!
  wx=(x-xtab(i1))/(xtab(i1+1)-xtab(i1))
  y=(1-wx)*ytab(i1)+wx*ytab(i1+1)
!
END
!
! ==================================================================!
!                                                                   !
! Dummy routine for ODEPACK. Necessary for Jacobian matrix if       !
! stiff ODEs.                                                       !
!                                                                   !
! ==================================================================!
!
SUBROUTINE jex()
  RETURN
END
#endif
