#include "cppdefs.h"
#ifdef ICEPLUME
!
! ==================================================================!
!                                                                   !
! This subroutine call the iceplume core functions to calcualte     !
! plume status/meltr ates/etc.                                      !
!                                                                   !
! ==================================================================!
!
SUBROUTINE ICEPLUME_CALC(ng, I,                                         &
     &                   fIni, tIni, sIni,                              &
     &                   sgTyp, sgDep, sgLen)
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
! ng      - grid identifier                                         !
! I       - point source index identifier                           !
! fIni    - subglacial runoff volume flux [m^3 s^-1]                !
! tIni    - subglacial runoff temperature [degC]                    !
! sIni    - subglacial runoff salinity [PSU]                        !
! sgTyp   - subglacial runoff type identifier                       !
! sgDep   - subglacial runoff inject depth [m]                      !
! sgLen   - subglacial runoff discharge length [m]                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! Local variables:                                                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! iceDepthK   - vertical grid index of grounding line depth         !
! plumeDepthK - vertical grid index of plume depth                  !
! osDepthK    - vertical grid index of overshoot depth              !
! areaP       - surface area of plume in contact with ice in        !
!               that cell [m^2]                                     !
! areaB       - surface area of plume in contact with ice out       !
!               that cell [m^2]                                     !
! meanVel     - ice tangental velocity [m s^-1]                     !
! det         - total detrainment volume [m^3 s^-1]                 !
! mB          - background meltrate [m s^-1]                        !
!                                                                   !
! ==================================================================!
!
! In/out variables
!
  integer, intent(in) :: ng, I
  real(r8), intent(inout) :: fIni, tIni, sIni
  integer, intent(in) :: sgTyp
  real(r8), intent(in) :: sgDep, sgLen
!
! Local variables declaration
!
  integer :: iceDepthK, plumeDepthK, osDepthK
  real(r8) :: areaP, areaB, meanVel, det, mB, tB, sB
  real(r8) :: RHO
!
  integer :: K, itrc
  real(r8) :: cff, cff1, cff2
!
! If discharge input salinity is less or equal than zero, set it to
! a small number
!
  IF (sIni .LE. 0.0) THEN
    sIni = 0.001
  ENDIF
!
! ==================================================================!
!                                                                   !
! Find iceDepthK                                                    !
!                                                                   !
! ==================================================================!
!
  iceDepthK = -1
  IF (sgDep .GE. 0.0) THEN
    iceDepthK = 0
  ELSE
    DO K = 0, N(ng)
      IF (PLUME(ng) % zW(I, K) .GE. sgDep) THEN
        iceDepthK = K
        EXIT
      ENDIF
    ENDDO
  ENDIF
!
! ==================================================================!
!                                                                   !
! Use entrainment plume model to calculate T, S, & flux in plume    !
!                                                                   !
! ==================================================================!
!
! Initiate
!
  plumeDepthK = iceDepthK+1
  det = 0.0
!
  PLUME(ng) % f(I, 0)    = 0.0
  PLUME(ng) % w(I, 0)    = 0.0
  PLUME(ng) % t(I, 0)    = 0.0
  PLUME(ng) % s(I, 0)    = 0.0
  PLUME(ng) % a(I, 0)    = 0.0
  PLUME(ng) % mInt(I, 0) = 0.0
!
  DO k = 1, N(ng)
    PLUME(ng) % f(I, K)    = 0.0
    PLUME(ng) % w(I, K)    = 0.0
    PLUME(ng) % t(I, K)    = 0.0
    PLUME(ng) % s(I, K)    = 0.0
    PLUME(ng) % a(I, K)    = 0.0
    PLUME(ng) % mInt(I, K) = 0.0
!
    PLUME(ng) % ent(I, K) = 0.0
    PLUME(ng) % det(I, K) = 0.0
!
    PLUME(ng) % detI(I, K)    = 0
    PLUME(ng) % detFrac(I, K) = 0.0
  ENDDO
!
  IF (fIni .GT. 0.0) THEN
!
! Run the plume model
!
    CALL ICEPLUME_ENTRAIN(ng, I, iceDepthK,                             &
     &                    fIni, tIni, sIni,                             &
     &                    sgTyp, sgLen)
!
! Calculate volume flux differential to give entrainment / extrainment
! First clear ent / det
!
    DO K = iceDepthK+1, N(ng)
      PLUME(ng) % ent(I, K) =                                           &
     &    -(PLUME(ng) % f(I, K) - PLUME(ng) % f(I, K-1))
      PLUME(ng) % ent(I, K) = PLUME(ng) % ent(I, K) +                   &
     &    MAX(PLUME(ng) % mInt(I, K) -                                  &
     &        PLUME(ng) % mInt(I, K-1), 0.0)
      PLUME(ng) % ent(I, K) = MIN(PLUME(ng) % ent(I, K), 0.0)
    ENDDO
!
! Separate entrainment / detrainment
!
    DO K = iceDepthK+1, N(ng)
      IF (PLUME(ng) % f(I, K) .GT. 0.0) THEN
        osDepthK = K
        det = PLUME(ng) % f(I, K)
      ENDIF
    ENDDO
!
! Replace the last value of lc and lm
!
    PLUME(ng) % lm(I, osDepthK) = PLUME(ng) % lm(I, osDepthK-1)
    PLUME(ng) % lc(I, osDepthK) = PLUME(ng) % lc(I, osDepthK-1)
!
! Find detrainment depth index, Calculate total detrainment volume
! and thickness
!
# if defined ICEPLUME_OVERSHOOT_NO && defined ICEPLUME_DETRAIN_NO
    plumeDepthK = osDepthK
# else
    DO K = osDepthK, iceDepthK+1, -1
      cff = RHO(PLUME(ng) % t(I, osDepthK),                             &
     &          PLUME(ng) % s(I, osDepthK),                             &
     &          PLUME(ng) % zR(I, K))
      IF (K .EQ. 1) THEN
        cff1 = PLUME(ng) % rhoAm(I, 1)
        cff2 = 0.5*(PLUME(ng) % rhoAm(I, 1)+                            &
     &              PLUME(ng) % rhoAm(I, 2))
      ELSEIF (K .EQ. N(ng)) THEN
        cff1 = 0.5*(PLUME(ng) % rhoAm(I, N(ng)-1)+                      &
     &              PLUME(ng) % rhoAm(I, N(ng)))
        cff2 = PLUME(ng) % rhoAm(I, N(ng))
      ELSE
        cff1 = 0.5*(PLUME(ng) % rhoAm(I, K-1)+PLUME(ng) % rhoAm(I, K))
        cff2 = 0.5*(PLUME(ng) % rhoAm(I, K)+PLUME(ng) % rhoAm(I, K+1))
      ENDIF
      IF ((K .EQ. N(ng)) .AND. (cff .LT. cff2)) THEN
        plumeDepthK = N(ng)
        EXIT
      ENDIF
      IF ((cff .LT. cff1) .AND. (cff .GT. cff2)) THEN
        plumeDepthK = K
        EXIT
      ENDIF
    ENDDO
# endif
  ENDIF
!
! ==================================================================!
!                                                                   !
! Calculate melt rates.                                             !
!                                                                   !
! ==================================================================!
!
! Initiate
!
  DO K = 1, N(ng)
    PLUME(ng) % m(I, K)  = 0.0
    PLUME(ng) % mB(I, K) = 0.0
  ENDDO
!
  DO K = iceDepthK+1, N(ng)
!
! ==================================================================!
!                                                                   !
! If there is a plume in that cell, then need to calculate plume    !
! meltrate distinct to background melt rate. Plume melt rate is     !
! already encorporated in the plume model, and taken into account   !
! in the temperature and salinity of the plume outflow. It is       !
! useful though to have it available as a diagnostic.               !
!                                                                   !
! ==================================================================!
!
    areaP = PLUME(ng) % a(I, K) - PLUME(ng) % a(I, K-1)
    areaP = MAX(areaP, 0.0)
    areaB = PLUME(ng) % dz(I, K)*PLUME(ng) % dy(I) - areaP
    areaB = MAX(areaB, 0.0)
!
    IF ( areaP .GT. 0.0 ) THEN
      PLUME(ng) % m(I, K) = PLUME(ng) % mInt(I, K) -                    &
     &                      PLUME(ng) % mInt(I, K-1)
    ENDIF
# ifdef ICEPLUME_MELT
!
! ==================================================================!
!                                                                   !
! Calculate the background melt rate (not generated by plumes).     !
! This will then be used to update the temperature and salinity in  !
! the adjacent cells. Velocities are calculated at cell faces -     !
! find averages for cell centres. Does not include velocity         !
! perpendicular to ice - this differs depending on orientation of   !
! ice front.                                                        !
!                                                                   !
! ==================================================================!
!
    IF ( areaB .GT. 0.0 )THEN
!
      meanVel = SQRT((PLUME(ng) % vAm(I, K))**2 +                       &
     &               (PLUME(ng) % wAm(I, K))**2)
      CALL ICEPLUME_MELTRATE(PLUME(ng) % tAm(I, K),                     &
     &                       PLUME(ng) % sAm(I, K),                     &
     &                       meanVel,                                   &
     &                       PLUME(ng) % zR(I, K),                      &
     &                       mB, tB, sB)
!
! Calculate integrated background melt rate [m3 s-1]
!
      PLUME(ng) % mB(I, K) = mB*areaB
    ENDIF
# endif
  ENDDO
!
! ==================================================================!
!                                                                   !
! Calculate active and passive tracers                              !
!                                                                   !
! Update log                                                        !
! Add active tracers, T and S into this scheme. This is necessary   !
! when distributing detrainment in several layers, since writting   !
! T & S in PLUME directly is in general not a good idea.            !
! 2018/4/20 Chuning Wang                                            !
!                                                                   !
! ==================================================================!
!
! Clear local plume tracer variables
!
  DO itrc = 1, NT(ng)
    PLUME(ng) % trc(I, itrc)    = 0.0
    PLUME(ng) % trcCum(I, itrc) = 0.0
!
    PLUME(ng) % trcB(I, itrc)   = 0.0
  ENDDO
!
  IF (fIni .GT. 0.0) THEN
!
! Write the active tracer (T & S) concentration to PLUME (ng) % trc
!
    PLUME(ng) % trc(I, itemp) = PLUME(ng) % t(I, osDepthK)
    PLUME(ng) % trc(I, isalt) = PLUME(ng) % s(I, osDepthK)
!
# ifdef ICEPLUME_TRACER
!
! If use passive tracers, calculate passive tracer concentrations
! in the runoff
!
    DO itrc = 3, NT(ng)
      PLUME(ng) % trcCum(I, itrc) = PLUME(ng) % trcIni(I, itrc) * fIni
    ENDDO
!
! Add up total sum of each tracer in plume
!
    DO K = iceDepthK+1, osDepthK
      IF (PLUME(ng) % ent(I, K) .LT. 0.0) THEN
        DO itrc = 3, NT(ng)
          PLUME(ng) % trcCum(I, itrc) =                                 &
     &        PLUME(ng) % trcCum(I, itrc) +                             &
     &        (-PLUME(ng) % ent(I, K) *                                 &
     &          PLUME(ng) % trcAm(I, K, itrc))
        ENDDO
      ENDIF
    ENDDO
!
! Calculate concentration of tracer in outflow 
!
    DO itrc = 3, NT(ng)
      PLUME(ng) % trc(I, itrc) = PLUME(ng) % trcCum(I, itrc) / det
    ENDDO
#  ifdef ICEPLUME_MELT_TRACER
!
! If use melt water tracer, rewrite tracer concentration
! Calculate accumulated trcer amount from base to detrain depth
!
    PLUME(ng) % trcCum(I, NT(ng)-1) = 0.0
    PLUME(ng) % trcCum(I, NT(ng)  ) = 0.0
    DO K = iceDepthK+1, osDepthK
      PLUME(ng) % trcCum(I, NT(ng)-1) =                                 &
     &    PLUME(ng) % trcCum(I, NT(ng)-1) +                             &
     &    (-PLUME(ng) % ent(I, K)*PLUME(ng) % trcAm(I, K, NT(ng)-1)) +  &
     &    PLUME(ng) % m(I, K)*PLUME(ng) % trcIni(I, NT(ng)-1)
      PLUME(ng) % trcCum(I, NT(ng)  ) =                                 &
     &    PLUME(ng) % trcCum(I, NT(ng)  ) +                             &
     &    (-PLUME(ng) % ent(I, K)*PLUME(ng) % trcAm(I, K, NT(ng)  ))
    ENDDO
!
! Calculate tracer concentration
!
    PLUME(ng) % trc(I, NT(ng)-1) = PLUME(ng) % trcCum(I, NT(ng)-1)/det
    PLUME(ng) % trc(I, NT(ng)  ) = PLUME(ng) % trcCum(I, NT(ng)  )/det
#  endif
# endif
  ENDIF
# ifdef ICEPLUME_MELT
!
! Rewrite tracer concentration in background melt water 
!
  PLUME(ng) % trcB(I, itemp) = (cI*tIce - L)/cW
  PLUME(ng) % trcB(I, isalt) = sIce
#  ifdef ICEPLUME_MELT_TRACER
!
! If use melt water tracer, rewrite the last tracer type.
!
  PLUME(ng) % trcB(I, NT(ng)) = PLUME(ng) % trcIni(I, NT(ng))
#  endif
# endif
!
! ==================================================================!
!                                                                   !
! Use a detrainment model to distribute the detrainment plume in    !
! several layers. Details in iceplume_detrain.F.                    !
!                                                                   !
! *** This part is completely rewritten since Feb 20, 2019.         !
!                                                                   !
! ==================================================================!
!
  IF (fIni .GT. 0) THEN
!
! Write detrainment in outputs.
!
# ifdef ICEPLUME_DETRAIN_NO
    PLUME(ng) % det(I, plumeDepthK) = det
    PLUME(ng) % detI(I, plumeDepthK) = 1
    PLUME(ng) % detFrac(I, plumeDepthK) = 1.0
# else
    CALL ICEPLUME_DETRAIN(ng, I,                                        &
     &                    iceDepthK, plumeDepthK, osDepthK,             &
     &                    PLUME(ng) % lc(I, osDepthK),                  &
     &                    det)
# endif
  ENDIF  ! fini .GT. 0.0
!
! Calculate depth integrated volume transport [m3 s-1]
!
  PLUME(ng) % trs(I) = 0.0
  DO K = 1, N(ng)
    PLUME(ng) % trs(I) = PLUME(ng) % trs(I) +                           &
     & PLUME(ng) % det(I, K) + PLUME(ng) % ent(I, K) +                  &
     & PLUME(ng) % mB(I, K)
  ENDDO
END SUBROUTINE ICEPLUME_CALC
!
! ==================================================================!
!                                                                   !
! This function calculates background melt rate.                    !
!                                                                   !
! ==================================================================!
!
SUBROUTINE ICEPLUME_MELTRATE(temp, salt, vel, depth, mdot, sB, tB)
!
  USE mod_iceplume
  implicit none
!
  real(r8) :: temp, salt, vel
  real(r8) :: depth, absVelocity
  real(r8) :: a, b, c
  real(r8), intent(inout) :: mdot, sB, tB
!
! Routine cannot cope with zero velocity.
! Unlikely to occur anyway with currents, waves, convection etc.
! This is not very physical, but will do for now.
!
  IF ( ABS(vel) .LT. velBkg ) vel = velBkg
!
! Calculate melt rate from 3 equation formualtion (as for plume
! models)
! Equations for Sb, Tb and mdot
!
  a = lambda1*(GamT*cW - GamS*cI)
  b = GamS*cI*(lambda1*salt - lambda2 - lambda3*depth +                 &
     &         tIce - (L/cI)) -                                         &
     &         GamT*cW*(temp - lambda2 - lambda3*depth + sIce)
  c = GamS*salt*(cI*(lambda2 + lambda3*depth - tIce) + L) +             &
     & GamT*sIce*cW*(temp - lambda2 - lambda3*depth)
!
  sB   = (1.0/(2.0*a))*(-b - SQRT(b**2 - 4.0*a*c))
  tB   = lambda1*sB + lambda2 + lambda3*depth
  mdot = GamS*SQRT(CdBkg)*ABS(vel)*(salt - sB)/sB
!
END SUBROUTINE ICEPLUME_MELTRATE
!
! ==================================================================!
!                                                                   !
! These subroutines are copied from MITgcm.                         !
!                                                                   !
! ==================================================================!
!
SUBROUTINE SW_TEMP(S, T, P, PR, rv)
!
! *=============================================================*
! | S/R  SW_TEMP
! | o compute in-situ temperature from potential temperature
! *=============================================================*
!
! REFERENCES:
! Fofonoff, P. and Millard, R.C. Jr
! Unesco 1983. Algorithms for computation of fundamental properties of
! seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
! Eqn.(31) p.39
! 
! Bryden, H. 1973.
! New Polynomials for thermal expansion, adiabatic temperature
! gradient and potential temperature of sea water.
! DEEP-SEA RES., 1973, Vol20,401-408.
!
  USE mod_kinds
  IMPLICIT NONE
!
! === Global variables ===
! 
! !INPUT/OUTPUT PARAMETERS:
! === Routine arguments ===
! S      :: salinity
! T      :: potential temperature
! P      :: pressure
! PR     :: reference pressure
! rv :: return value (in-situ temeparture in degree C)
!
  real(r8) :: S, T, P, PR
  real(r8) :: rv
!
! !LOCAL VARIABLES:
! === local variables ===
!
  CALL SW_PTMP  (S, T, PR, P, rv)
!
END
!
SUBROUTINE SW_PTMP  (S, T, P, PR, rv)
!
! !DESCRIPTION: \bv
! *=============================================================*
! | S/R  SW_PTMP
! | o compute potential temperature as per UNESCO 1983 report.
! *=============================================================*
! \ev
! started:
!          Armin Koehl akoehl@ucsd.edu
!
! ==================================================================
! SUBROUTINE SW_PTMP
!
! ==================================================================
!
  USE mod_kinds
  IMPLICIT NONE
!
! === Global variables ===
!
! !INPUT/OUTPUT PARAMETERS:
! === Routine arguments ===
! S  :: salinity    [psu      (PSS-78) ]
! T  :: temperature [degree C (IPTS-68)]
! P  :: pressure    [db]
! PR :: Reference pressure  [db]
! rv :: return value (potential temeparture in degree C)
!
  real(r8) :: S, T, P, PR
  real(r8) :: rv
!
! !LOCAL VARIABLES
! === local variables ===
!
  real(r8) :: del_P ,del_th, th, q
  real(r8) :: onehalf, two, three
  parameter ( onehalf = 0.5d0, two = 2.d0, three = 3.d0 )
  real(r8) :: adtg_val
!
! theta1
!
  del_P   = PR - P
  call sw_adtg(S, T, P, adtg_val)
  del_th  = del_P*adtg_val
  th      = T + onehalf*del_th
  q       = del_th
!
! theta2
!
  call sw_adtg(S, th, P+onehalf*del_P, adtg_val)
  del_th  = del_P*adtg_val
  th      = th + (1 - 1/sqrt(two))*(del_th - q)
  q       = (two-sqrt(two))*del_th + (-two+three/sqrt(two))*q
!
! theta3
!
  call sw_adtg(S, th, P+onehalf*del_P, adtg_val)
  del_th  = del_P*adtg_val
  th      = th + (1 + 1/sqrt(two))*(del_th - q)
  q       = (two + sqrt(two))*del_th + (-two-three/sqrt(two))*q
!
! theta4
!
  call sw_adtg(S, th, P+del_P, adtg_val)
  del_th  = del_P*adtg_val
  rv      = th + (del_th - two*q)/(two*three)
!
END
!
SUBROUTINE SW_ADTG  (S,T,P, rv)
!
! !DESCRIPTION: \bv
! *=============================================================*
! | S/R  SW_ADTG
! | o compute adiabatic temperature gradient as per UNESCO 1983 routines.
! *=============================================================*
! \ev
!
! started:
!          Armin Koehl akoehl@ucsd.edu
!
! !USES:
!
  USE mod_kinds
  IMPLICIT NONE
!
! === Global variables ===
!
! !INPUT/OUTPUT PARAMETERS:
! === Routine arguments ===
!
  real(r8) :: S,T,P
  real(r8) :: rv
!
! !LOCAL VARIABLES:
! === local variables ===
!
  real(r8) :: a0,a1,a2,a3,b0,b1,c0,c1,c2,c3,d0,d1,e0,e1,e2
  real(r8) :: sref
!
  sref = 35.d0
  a0 =  3.5803d-5
  a1 = +8.5258d-6
  a2 = -6.836d-8
  a3 =  6.6228d-10
!
  b0 = +1.8932d-6
  b1 = -4.2393d-8
!
  c0 = +1.8741d-8
  c1 = -6.7795d-10
  c2 = +8.733d-12
  c3 = -5.4481d-14
!
  d0 = -1.1351d-10
  d1 =  2.7759d-12
!
  e0 = -4.6206d-13
  e1 = +1.8676d-14
  e2 = -2.1687d-16
!
  rv =  a0 + (a1 + (a2 + a3*T)*T)*T                                     &
     &  + (b0 + b1*T)*(S-sref)                                          &
     &  + ( (c0 + (c1 + (c2 + c3*T)*T)*T) + (d0 + d1*T)*(S-sref) )*P    &
     &  + (  e0 + (e1 + e2*T)*T )*P*P
!
END
#endif
