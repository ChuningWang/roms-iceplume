#include "cppdefs.h"
#ifdef ICEPLUME
!
! ==================================================================!
!                                                                   !
! Use this function to calculate the velocity of detrainment.       !
!                                                                   !
! ==================================================================!
!
SUBROUTINE ICEPLUME_DETRAIN(ng, I,                                      &
     &                      iceDepthK, plumeDepthK, osDepthK,           &
     &                      lc, det)
  USE mod_iceplume
  implicit none
!
! ==================================================================!
!                                                                   !
! Input/Output variables:                                           !
!                                                                   !
! ==================================================================!
!                                                                   !
! ng          - grid identifier                                     !
! I           - point source index identifier                       !
! iceDepthK   - glacier grounding line depth vertical grid index    !
! plumeDepthK - plume detrainment depth vertical grid index         !
! osDepthK    - plume overshoot depth vertical grid index           !
! lc          - plume/ocean contact length [m]                      !
! det         - detrainment volume flux [m^3 s^-1]                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! Local variables:                                                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! detDz       - detrainment layer total thickness [m]               !
! minDz       - minimum detrainment layer thickness [m]             !
! rhoP        - plume density [kg m^-3]                             !
! potE        - potential energy of plume/ambient water             !
!                                                                   !
! detVel      - detrainment velocity [m s^-1]                       !
! detWeight   - detrainment weight function                         !
! detSum      - detrainment volume flux before normalize [m^3 s^-1] !
!                                                                   !
! isSurface   - logical, if the plume discharge at surface          !
! isBottom    - logical, if the plume discharge at bottom           !
! udSwitch    - switch to determine search one layer up(1)/down(-1) !
!                                                                   !
! maxVel      - maximum velocity based on Richardson number         !
!               criterion [m s^-1]                                  !
! Fr          - Froude number                                       !
! rho1        - upper layer density [kg m^-3]                       !
! rho2        - lower layer density [kg m^-3]                       !
! h1          - upper layer thickness [m]                           !
! h2          - lower layer thickness [m]                           !
! gRed        - reduced gravity [m s^-2]                            !
# ifdef ICEPLUME_DET_NEUTRAL
!                                                                   !
! volR        - detrainment volume of runoff [m^3]                  !
! volM        - detrainment volume of meltwater                     !
! volF        - detrainment volume of total freshwater              !
! volE        - detrainment volume of entrainment                   !
! volP        - detrainment volume of plume                         !
!                                                                   !
! tF          - temperature of total freshwater [degC]              !
! tE          - temperature of entrainment                          !
! tP          - temperature of plume                                !
!                                                                   !
! sF          - salinity of total freshwater [PSU]                  !
! sE          - salinity of entrainment                             !
! sP          - salinity of plume                                   !
!                                                                   !
! rhoF        - density of total freshwater [kg m^-3]               !
! rhoE        - density of entrainment                              !
! rhoP        - density of plume                                    !
# endif
!                                                                   !
! ==================================================================!
!
! In/out variables
!
  integer, intent(in) :: ng, I
  integer, intent(in) :: iceDepthK, plumeDepthK, osDepthK
  real(r8), intent(in) :: lc, det
!
! Local variables declaration
!
  real(r8) :: detDz, minDz, detVel, detWeight, potE
  real(r8) :: rhoP
  integer  :: udSwitch
  logical  :: isSurface = .FALSE.
  logical  :: isBottom = .FALSE.
  real(r8) :: RHO
!
! For detrainment weight function
!
  real(r8) :: detSum
!
! For calculation of internal wave speed
!
  real(r8) :: maxVel, Fr, Ri
  real(r8) :: rho1, rho2, h1, h2, gRed
  real(r8) :: rhoH, h01, h02
# ifdef ICEPLUME_DET_NEUTRAL
!
! For variable detrainment tracer concentration
!
  real(r8) :: volR, volM, volF, volE, volP
  real(r8) :: tF, tE, tP, sF, sE, sP, rhoF, rhoE
# endif
!
! Other local variables
!
  integer  :: KI, K, K1, K2, itrc, counter
  real(r8) :: cff, cff1, cff2, cff3, cff4
# ifdef ICEPLUME_DET_NEUTRAL
!
! Initialize.
!
  DO K = 1, N(ng)
    PLUME(ng) % detF(I, K) = 0.0
    PLUME(ng) % detE(I, K) = 0.0
  ENDDO
  DO itrc = 1, NT(ng)
    DO K = 1, N(ng)
      PLUME(ng) % detTrc(I, K, itrc) = 0.0
    ENDDO
  ENDDO
# endif
!
! ==================================================================!
!                                                                   !
! PART I - Find minimum detrainment depth                           !
!                                                                   !
! ==================================================================!
!
! Initiate logical variables
!
  IF (plumeDepthK .EQ. N(ng)) THEN
    isSurface = .TRUE.
  ELSE
    isSurface = .FALSE.
  ENDIF
  IF (plumeDepthK .EQ. 1) THEN
    isBottom = .TRUE.
  ELSE
    isBottom = .FALSE.
  ENDIF
!
! Compute minimum detrainment thickness.
!
  rhoP = RHO(PLUME(ng) % t(I, osDepthK), PLUME(ng) % s(I, osDepthK),    &
     &       PLUME(ng) % zR(I, plumeDepthK))
!
! Use the Parameterization of Ching et al. (1993) to determine the bulk
! velocity and nt in several layers. First, calculate BV frequency.
!
  IF (isSurface) THEN
!
! Strong density jump. Always use momentum conservation.
!
    maxVel = 0.95*PLUME(ng) % w(I, N(ng))
  ELSE
!
! Use piecewise function.
!
    KI = plumeDepthK
    rhoH = 0.0
    h1 = 0.0
    DO K = KI+1, N(ng)
      h1 = h1+PLUME(ng) % dz(I, K)
      rhoH = rhoH + PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
    ENDDO
    rho1 = rhoH/h1
!
    rhoH = 0.0
    h2 = 0.0
    DO K = 1, KI 
      h2 = h2+PLUME(ng) % dz(I, K)
      rhoH = rhoH + PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
    ENDDO
    rho2 = rhoH/h2
!
    gRed = MAX(-g*(rho1-rho2)/rhoRef, gRedBkg)
    KI = MAXLOC(PLUME(ng) % w(I, :), 1)
    cff1 = PLUME(ng) % f(I, KI)/PLUME(ng) % w(I, KI)/lc
    Ri = gRed*cff1/(PLUME(ng) % w(I, KI)**2)
    IF (Ri .LT. 6.0) THEN
      maxVel = (0.7*(Ri**0.17))*PLUME(ng) % w(I, KI)
    ELSE
      maxVel = 0.95*PLUME(ng) % w(I, KI)
    ENDIF
  ENDIF
  maxVel = maxVel*(SQRT(pi/2)*detSigma)
  minDz = det/lc/maxVel
!
! Determine the vertical extent of detrainment plume.
!
  KI = plumeDepthK
  K1 = KI
  K2 = KI
  detDz = 0.0
!
! Add layers around plumeDepthK until it reaches critical thickness.
!
  counter = 0
  PLUME(ng) % detI(I, KI) = 1
  PLUME(ng) % detFrac(I, KI) = 1.0
  detDz = PLUME(ng) % dz(I, KI)
  potE = g*PLUME(ng) % dz(I, KI)*(PLUME(ng) % rhoAm(I, KI)-rhoP)
  IF (potE .GT. 0.0) THEN
    udSwitch = 1
  ELSE
    udSwitch = -1
  ENDIF
  DO WHILE ( (detDz .LT. minDz) .AND. (counter .LT. 100) )
    counter = counter + 1
    IF (udSwitch .EQ. 1) THEN
!
! Search one layer up
!
      DO K = 2, N(ng)
        IF ((PLUME(ng) % detI(I, K  ) .EQ. 0) .AND.                     &
     &      (PLUME(ng) % detI(I, K-1) .EQ. 1)) THEN
          KI = K
        ENDIF
      ENDDO
      K1 = KI
    ELSEIF (udSwitch .EQ. -1) THEN
!
! Search one layer down
!
      DO K = 1, N(ng)-1
        IF ((PLUME(ng) % detI(I, K+1) .EQ. 1) .AND.                     &
     &      (PLUME(ng) % detI(I, K  ) .EQ. 0)) THEN
          KI = K
        ENDIF
      ENDDO
      K2 = KI
    ENDIF
    IF (KI .EQ. N(ng)) isSurface = .TRUE.
    IF (KI .EQ. 1) isBottom = .TRUE.
!
! Update detrainment layer flag and detrainment depth.
!
    PLUME(ng) % detI(I, KI) = 1
    detDz = detDz + PLUME(ng) % dz(I, KI)
    PLUME(ng) % detFrac(I, KI) = 1.0
!
! Determine next layer to minimize potential energy anomaly.
!
    potE = potE +                                                       &
     &     g*PLUME(ng) % dz(I, KI)*(PLUME(ng) % rhoAm(I, KI)-rhoP)
    IF (potE .GT. 0.0) THEN
      udSwitch = 1
    ELSE
      udSwitch = -1
    ENDIF
!
! IF layer reached surface/bottom, force udSwitch to -1/1
!
    IF     ( (       isSurface ) .AND. ( .NOT. isBottom ) ) THEN
      udSwitch = -1
    ELSEIF ( ( .NOT. isSurface ) .AND. (       isBottom ) ) THEN
      udSwitch = 1
    ELSEIF ( (       isSurface ) .AND. (       isBottom ) ) THEN
      udSwitch = 0
    ENDIF
!
! Other exit loop criteria.
!
    IF ( (minDz .EQ. (minDz+1.0)) .OR. (udSwitch .EQ. 0) ) EXIT
  ENDDO
!
! Only take partial grid from top/bottom layer to satisfy potE=0
!
  IF ( (.NOT. isSurface) .AND. (counter .GT. 0) ) THEN
    IF ( potE .GT. 0.0 ) THEN
      KI = K2
    ELSE
      KI = K1
    ENDIF
    cff = 1 - potE / (g*(PLUME(ng) % rhoAm(I, KI)-rhoP)) /              &
     &        PLUME(ng) % dz(I, KI)
    IF ( (cff .GE. 0) .AND. (cff .LT. 1) ) THEN
      PLUME(ng) % detFrac(I, KI) = cff
      detDz = detDz - PLUME(ng) % dz(I, KI)*(1 - cff)
    ENDIF
  ENDIF
!
! ==================================================================!
!                                                                   !
! PART II - Distribute detrainmnet in several layers                !
!                                                                   !
! Update log                                                        !
! Use a Gause function to smooth the distribution.                  !
! 2018/06/08 Chuning Wang                                           !
!                                                                   !
! ==================================================================!
!
! Update detrainment volume flux
!
  detVel = det / (lc * detDz)
  DO K = 1,N(ng)
    IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
!
! First, calculate weight using Gaussian function.
!
      IF ( K .GT. plumeDepthK ) THEN
        cff = (PLUME(ng) % zR(I, K) -PLUME(ng) % zR(I, plumeDepthK)) /  &
     &    (PLUME(ng) % zR(I, K1)-PLUME(ng) % zR(I, plumeDepthK)+0.0001)
      ELSE
        cff = (PLUME(ng) % zR(I, K) -PLUME(ng) % zR(I, plumeDepthK)) /  &
     &    (PLUME(ng) % zR(I, K2)-PLUME(ng) % zR(I, plumeDepthK)+0.0001)
      ENDIF
      detWeight = EXP(-0.5 * (cff/detSigma)**2)
!
! Distribute the detrainment.
!
      PLUME(ng) % det(I, K) = detWeight * detVel * lc *                 &
     &  PLUME(ng) % detFrac(I, K) * PLUME(ng) % dz(I, K)
    ENDIF
  ENDDO
!
! Normalize
!
  detSum = SUM(PLUME(ng) % det(I, :))
  DO K = 1,N(ng)
    PLUME(ng) % det(I, K) = PLUME(ng) % det(I, K)*det/detSum
  ENDDO
# ifdef ICEPLUME_DET_NEUTRAL
!
! Calculate proportion of each water mass.
! There are three main water masses -
! 1. Runoff (R),
! 2. Melt (M),
! 3. Entrainment (E).
! R and M combined forms a new water mass,
! 4. Freshwater (F).
! M and F combined froms the final water mass,
! 5. Plume (P).
!
  volP = PLUME(ng) % f(I, osDepthK)
  volR = PLUME(ng) % f(I, iceDepthK)
  volM = PLUME(ng) % mInt(I, osDepthK)
  volF = volR + volM 
  volE = volP - volF
!
! Temperature. Note that tGade is effective temperature of Gade line.
!
  tP = PLUME(ng) % trc(I, itemp)
  tF = (tGade*volM + PLUME(ng) % t(I, iceDepthK)*volR)/volF
  tE = (tP*volP-tF*volF)/volE
!
! Salinity.
!
  sP = PLUME(ng) % trc(I, isalt)
  sF = (sIce*volM  + PLUME(ng) % s(I, iceDepthK)*volR)/volF
  sE = (sP*volP-sF*volF)/volE
!
! Density.
!
  rhoF = RHO(tF, sF, PLUME(ng) % zR(I, plumeDepthK))
  rhoE = (rhoP*volP-rhoF*volF)/volE
!
! Calculate proportion from mixing line.
!
  DO K = 1, N(ng)
    IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
      cff = (PLUME(ng) % zR(I, K) -                                     &
     &       PLUME(ng) % zR(I, plumeDepthK))*0.005
      cff1 = MIN(MAX((PLUME(ng) % rhoAm(I, K)-rhoE)/(rhoF-rhoE),        &
     &               0.01), 0.99)
      cff1 = ((PLUME(ng) % rhoAm(I, K)+cff)-rhoE)/(rhoF-rhoE)
      cff1 = MIN(MAX(cff1, 0.01), 0.99)
      cff2 = 1-cff1
      PLUME(ng) % detF(I, K) = PLUME(ng) % det(I, K)*cff1
      PLUME(ng) % detE(I, K) = PLUME(ng) % det(I, K)*cff2
    ENDIF
  ENDDO
!
! Make corrections to make sure each water mass sums to its previous
! volume.
!
  cff1 = volF/SUM(PLUME(ng) % detF(I, :))
  cff2 = volE/SUM(PLUME(ng) % detE(I, :))
  DO K = 1, N(ng)
    IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
      PLUME(ng) % detF(I, K) = PLUME(ng) % detF(I, K)*cff1
      PLUME(ng) % detE(I, K) = PLUME(ng) % detE(I, K)*cff2
    ENDIF
  ENDDO
!
! Rewrite detrainment volume.
!
  DO K = 1, N(ng)
    IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
      PLUME(ng) % det(I, K) = PLUME(ng) % detF(I, K) +                  &
     &                        PLUME(ng) % detE(I, K)
    ENDIF
  ENDDO
!
! ==================================================================!
!                                                                   !
! Update tracer concentration in detrainment model.                 !
!                                                                   !
! ==================================================================!
!
! Active tracers.
!
  DO K = 1, N(ng)
    IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
      PLUME(ng) % detTrc(I, K, itemp) =                                 &
     &    (tF*PLUME(ng) % detF(I, K)+tE*PLUME(ng) % detE(I, K)) /       &
     &    (PLUME(ng) % detF(I, K)+PLUME(ng) % detE(I, K))
      PLUME(ng) % detTrc(I, K, isalt) =                                 &
     &    (sF*PLUME(ng) % detF(I, K)+sE*PLUME(ng) % detE(I, K)) /       &
     &    (PLUME(ng) % detF(I, K)+PLUME(ng) % detE(I, K))
    ENDIF
  ENDDO
!
! Passive tracers.
!
  DO itrc = 3, NT(ng)
    cff1 = PLUME(ng) % trcIni(I, itrc)*volR/volF
#  ifdef ICEPLUME_MELT_TRACER
    IF ( itrc .EQ. NT(ng)-1 ) THEN
      cff1 = PLUME(ng) % trcIni(I, itrc)*volM/volF
    ELSEIF ( itrc .EQ. NT(ng) ) THEN
      cff1 = 0.0
    ENDIF
#  endif
    cff2 = (PLUME(ng) % trc(I, itrc)*volP - cff1*volF)/volE
    DO K = 1, N(ng)
      IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
        PLUME(ng) % detTrc(I, K, itrc) =                                &
     &    (cff1*PLUME(ng) % detF(I, K)+cff2*PLUME(ng) % detE(I, K)) /   &
     &    (PLUME(ng) % detF(I, K)+PLUME(ng) % detE(I, K))
      ENDIF
    ENDDO
  ENDDO
# endif
END SUBROUTINE ICEPLUME_DETRAIN
#endif
