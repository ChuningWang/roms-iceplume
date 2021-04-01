#include "cppdefs.h"
      MODULE mod_iceplume
#ifdef ICEPLUME
! =====================================================================!
!                                                                      !
! These are the module functions of iceplume model.                    !
!                                                                      !
! =====================================================================!
!
! This module stores all global variables.
!
        USE mod_kinds
        USE mod_param, ONLY: Ngrids, N, NT
        USE mod_scalars, ONLY: isalt, itemp, dt, ntstart, iic
        USE mod_sources, ONLY: Nsrc
!
        implicit none
!
! =====================================================================!
!                                                                      !
! Model parameters                                                     !
!                                                                      !
! alpha     - entrainment rate                                         !
! tIce      - ice temperature [degC]                                   !
! sIce      - ice salinity [PSU]                                       !
! rhoRef    - reference density [kg m^-3]                              !
! rhoAir    - air density [kg m^-3]                                    !
! g         - gravity acceleration [m s^-2]                            !
! cW        - heat capacity of water [J kg^-1 degC^-1]                 !
! cI        - heat capacity of ice [J kg^-1 degC^-1]                   !
! L         - latent heat of melting [J kg^-1]                         !
! tGade     - equavilent temperature of meltwater (Gade Line) [degC]   !
! lambda1   - freezing point slope [degC PSU^-1]                       !
! lambda2   - freezing point offset [degC]                             !
! lambda3   - freezing point depth slope [degC m^-1]                   !
!                                                                      !
! GamT      - thermal turbulent transfer coefficient                   !
! GamS      - salt turbulent transfer coefficient                      !
! Cd        - ice-plume drag coefficient                               !
!                                                                      !
! RiB       - critical Richardson number                               !
! gRedBkg   - background reduced gravity                               !
! CdBkg     - background ice-plume drag coefficient                    !
! velBkg    - background velocity [m s^-1]                             !
! wIni      - initial (discharge) velocity [m s^-1]                    !
!                                                                      !
! detSigma  - detrainment Gaussian profile standard deviation          !
# ifdef ICEPLUME_DET_AVERAGE
! avgExp    - logical switch to activate temporal exponential average  !
! avgAlpha  - average weight for exponential smoothing (0 to 1)        !
# endif
!                                                                      !
! =====================================================================!
!
        real(r8), parameter :: pi = 4.0d0*ATAN(1.0d0)    ! Pi
!
        real(r8), parameter :: alpha      = 0.1_r8
        real(r8), parameter :: tIce       = -10.0_r8
        real(r8), parameter :: sIce       = 0.0_r8
        real(r8), parameter :: rhoRef     = 1020.0_r8
        real(r8), parameter :: rhoAir     = 1.225_r8
        real(r8), parameter :: g          = 9.81_r8
        real(r8), parameter :: cW         = 3974.0_r8
        real(r8), parameter :: cI         = 2000.0_r8
        real(r8), parameter :: L          = 335000.0_r8
        real(r8), parameter :: tGade      = -(L - cI*tIce)/cW
        real(r8), parameter :: lambda1    = -0.0573_r8
        real(r8), parameter :: lambda2    = 0.0832_r8
        real(r8), parameter :: lambda3    = 0.000761_r8
!
        real(r8), parameter :: GamT       = 0.0220_r8
        real(r8), parameter :: GamS       = 0.000620_r8
        real(r8), parameter :: Cd         = 0.065_r8
! Original value from Cowton et al. 2015
!        real(r8), parameter :: Cd         = 0.0025_r8
!
        real(r8), parameter :: RiB        = 1.0_r8
        real(r8), parameter :: gRedBkg    = 0.01_r8
        real(r8), parameter :: CdBkg      = 0.065_r8
! Original value from Cowton et al. 2015
!        real(r8), parameter :: CdBkg      = 0.0025_r8
        real(r8), parameter :: velBkg     = 0.03_r8
        real(r8), parameter :: wIni       = 0.1_r8
!
        real(r8), parameter :: detSigma   = 0.5_r8
# ifdef ICEPLUME_DET_AVERAGE
        logical :: avgExp
        real(r8), parameter :: avgAlpha   = 1.0_r8
# endif
!
! =====================================================================!
!                                                                      !
! PLUME Type variables                                                 !
!                                                                      !
! Isrc, Jsrc - grid indices acquired from river file.                  !
! dx, dy     - grid lengths in cross/along glacier direction.          !
! dir        - direction of plume. +1 for positve direction and        !
!              -1 for negative direction. 0 for other situation.       !
! trs        - depth integrated total transport [m^3 s^-1]             !
# ifdef ICEPLUME_DET_AVERAGE
! avgIstr    - Location indices (Istr) when averaging density profiles.!
! avgIend    - Location indices (Iend) when averaging density profiles.!
! avgJstr    - Location indices (Jstr) when averaging density profiles.!
! avgJend    - Location indices (Jend) when averaging density profiles.!
! avgArea    - surface area of average span for each plume source.     !
# endif
!                                                                      !
! Profiles                                                             !
!                                                                      !
! zW         - depth [m]                                               !
! f          - plume vertical volume flux [m^3 s^-1]                   !
! w          - plume vertical velocity [m s^-1]                        !
! t          - plume temperature [degC]                                !
! s          - plume salinity [PSU]                                    !
! a          - plume area integrated [m^2]                             !
! mInt       - plume area integrated melt [m^3 s^-1]                   !
! rho        - plume density [kg m^-3]                                 !
!                                                                      !
! zR         - depth at Rho points [m]                                 !
! sAm        - ambient salinity [PSU]                                  !
! tAm        - ambient temperature [degC]                              !
! vAm        - horizontal velocity parallel to glacier [m s^-1]        !
! wAm        - vertical velocity parallel to glacier [m s^-1]          !
! tpAm       - ambient potential temperature [degC]                    !
! rhoAm      - ambient density [kg m^-3]                               !
# ifdef ICEPLUME_DET_AVERAGE
! rhoAm0     - ambient density of the previous step [kg m^-3]          !
# endif
!                                                                      !
! lm         - plume/glacier contact length [m]                        !
! lc         - plume/water contact length [m]                          !
!                                                                      !
! ent        - entrainment rate [m^3 s^-1]                             !
! det        - detrainment rate [m^3 s^-1]                             !
! detI       - detrainment flag                                        !
! detFrac    - fraction of detrainment in vertical direction           !
# ifdef ICEPLUME_DET_NEUTRAL
!                                                                      !
! detF       - detrainment rate of freshwater [m^3 s^-1]               !
! detE       - detrainment rate of entrainment [m^3 s^-1]              !
! detTrc     - detrainment tracer concentration                        !
# endif
!                                                                      !
! m          - plume melt rate [m^3 s^-1]                              !
! mB         - background melt rate [m^3 s^-1]                         !
!                                                                      !
! dz         - RHO layer thickness [m]                                 !
!                                                                      !
! Passive tracers                                                      !
!                                                                      !
! trc        - tracer concentration                                    !
! trcAm      - ambient tracer concentration                            !
! trcB       - background meltwater tracer concentration               !
! trcCum     - accumulative tracer concentration                       !
! trcIni     - initial tracer concentration in discharge               !
# ifdef ICEPLUME_MIX
!                                                                      !
! trsW       - depth integrated total transport through                !
!              con/divergence [m^3 s^-1]                               !
! entW       - entrainment rate through con/divergence [m^3 s^-1]      !
! detW       - detrainment rate through con/divergence [m^3 s^-1]      !
! mBW        - background melt through con/divergence [m^3 s^-1]       !
!                                                                      !
! For calculation of virtual tracer flux                               !
!                                                                      !
! dzNew      - new layer thickness [m]                                 !
! zWNew      - new omega surface depth [m]                             !
! zLNew      - new box vertical boundary [m]                           !
! dzWei      - new layer weight function                               !
! trcL       - new box tracer concentration                            !
! trcNew     - updated tracer concentration                            !
# endif
# ifdef ICEPLUME_WRT_AVERAGE
!                                                                      !
! If write plume states into average output files, define average      !
!   fields                                                             !
!                                                                      !
! sument     - total plume entrainment volume flux                     !
! summ       - total plume melt rate                                   !
! summB      - total ambient melt rate                                 !
! sumtrcEnt  - total plume entrainment tracer flux                     !
!                                                                      !
! avgent     - average plume entrainment volume flux                   !
! avgm       - average plume melt rate                                 !
! avgmB      - average ambient melt rate                               !
! avgtrcEnt  - average plume entrainment tracer flux                   !
!                                                                      !
# endif
!                                                                      !
! =====================================================================!
!
        TYPE T_PLUME
!
! Variables.
!
!
! Depth independent variables.
!
          integer,  pointer :: Isrc(:)
          integer,  pointer :: Jsrc(:)
          real(r8), pointer :: dx(:)
          real(r8), pointer :: dy(:)
          real(r8), pointer :: dir(:)
          real(r8), pointer :: trs(:)
# ifdef ICEPLUME_DET_AVERAGE
          integer,  pointer :: avgIstr(:)
          integer,  pointer :: avgIend(:)
          integer,  pointer :: avgJstr(:)
          integer,  pointer :: avgJend(:)
          real(r8), pointer :: avgArea(:, :)
# endif
!
! Plume state (omega surface).
!
          real(r8), pointer :: zW(:, :)
          real(r8), pointer :: f(:, :)
          real(r8), pointer :: w(:, :)
          real(r8), pointer :: t(:, :)
          real(r8), pointer :: s(:, :)
          real(r8), pointer :: a(:, :)
          real(r8), pointer :: mInt(:, :)
          real(r8), pointer :: rho(:, :)
!
! Grid, ambient state (rho surface).
!
          real(r8), pointer :: zR(:, :)
          real(r8), pointer :: tAm(:, :)
          real(r8), pointer :: sAm(:, :)
          real(r8), pointer :: vAm(:, :)
          real(r8), pointer :: wAm(:, :)
          real(r8), pointer :: tpAm(:, :)
          real(r8), pointer :: rhoAm(:, :)
# ifdef ICEPLUME_DET_AVERAGE
          real(r8), pointer :: rhoAm0(:, :)
# endif
!
! Plume shape parameters.
!
          real(r8), pointer :: lm(:, :)
          real(r8), pointer :: lc(:, :)
!
! Volume fluxes.
!
          real(r8), pointer :: ent(:, :)
          real(r8), pointer :: det(:, :)
          integer,  pointer :: detI(:, :)
          real(r8), pointer :: detFrac(:, :)
# ifdef ICEPLUME_DET_NEUTRAL
!
! For the neutral buoyancy detrainment model.
!
          real(r8), pointer :: detF(:, :)
          real(r8), pointer :: detE(:, :)
          real(r8), pointer :: detTrc(:, :, :)
# endif
!
! Melt rate, freshwater and heat fluxes.
!
          real(r8), pointer :: m(:, :)
          real(r8), pointer :: mB(:, :)
!
! Other profiles.
!
          real(r8), pointer :: dz(:, :)
!
! Passive tracer concentration.
!
          real(r8), pointer :: trcAm(:, :, :)
          real(r8), pointer :: trcB(:, :)
          real(r8), pointer :: trc(:, :)
          real(r8), pointer :: trcCum(:, :)
          real(r8), pointer :: trcIni(:, :)
# ifdef ICEPLUME_MIX
!
! Volume fluxes through convergence/divergence.
!
          real(r8), pointer :: trsW(:)
          real(r8), pointer :: entW(:, :)
          real(r8), pointer :: detW(:, :)
          real(r8), pointer :: mBW (:, :)
!
! For calculation of virtual tracer flux.
!
          real(r8), pointer :: dzNew(:, :)
          real(r8), pointer :: zWNew(:, :)
          real(r8), pointer :: zLNew(:, :)
          real(r8), pointer :: dzWei(:, :)
          real(r8), pointer :: trcL(:, :, :)
          real(r8), pointer :: trcNew(:, :, :)
# endif
# ifdef ICEPLUME_WRT_AVERAGE
!
! For adding up model states.
!
          real(r8), pointer :: sument(:)
          real(r8), pointer :: summ(:)
          real(r8), pointer :: summB(:)
          real(r8), pointer :: sumtrcEnt(:, :)
!
! For storage of model averages.
!
          real(r8), pointer :: avgent(:)
          real(r8), pointer :: avgm(:)
          real(r8), pointer :: avgmB(:)
          real(r8), pointer :: avgtrcEnt(:, :)
# endif
        END TYPE T_PLUME
!
        TYPE (T_PLUME), allocatable :: PLUME(:)
!
      CONTAINS
!
      SUBROUTINE allocate_iceplume(ng)
!
! =====================================================================!
!                                                                      !
! Allocate PLUME Type variables.                                       !
!                                                                      !
! =====================================================================!
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer :: ng
!
!  Local variable declarations.
!
      integer :: is, itrc, k
      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng .EQ. 1) allocate( PLUME(Ngrids) )
!
!  Allocate iceplume Sources/Sinks variables.
!
      allocate ( PLUME(ng) % Isrc (Nsrc(ng)) )
      allocate ( PLUME(ng) % Jsrc (Nsrc(ng)) )
      allocate ( PLUME(ng) % dx   (Nsrc(ng)) )
      allocate ( PLUME(ng) % dy   (Nsrc(ng)) )
      allocate ( PLUME(ng) % dir  (Nsrc(ng)) )
      allocate ( PLUME(ng) % trs  (Nsrc(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*6,r8)
# ifdef ICEPLUME_DET_AVERAGE
      allocate ( PLUME(ng) % avgIstr (Nsrc(ng)) )
      allocate ( PLUME(ng) % avgIend (Nsrc(ng)) )
      allocate ( PLUME(ng) % avgJstr (Nsrc(ng)) )
      allocate ( PLUME(ng) % avgJend (Nsrc(ng)) )
      allocate ( PLUME(ng) % avgArea (Nsrc(ng), N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*4+Nsrc(ng)*N(ng),r8)
# endif
!
      allocate ( PLUME(ng) % zW   (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % f    (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % w    (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % t    (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % s    (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % a    (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % mInt (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % rho  (Nsrc(ng), 0:N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*(N(ng)+1)*8,r8)
!
      allocate ( PLUME(ng) % zR    (Nsrc(ng), N(ng)  ) )
      allocate ( PLUME(ng) % tAm   (Nsrc(ng), N(ng)  ) )
      allocate ( PLUME(ng) % sAm   (Nsrc(ng), N(ng)  ) )
      allocate ( PLUME(ng) % vAm   (Nsrc(ng), N(ng)  ) )
      allocate ( PLUME(ng) % wAm   (Nsrc(ng), N(ng)  ) )
      allocate ( PLUME(ng) % tpAm  (Nsrc(ng), N(ng)  ) )
      allocate ( PLUME(ng) % rhoAm (Nsrc(ng), N(ng)+1) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng)*7+1,r8)
# ifdef ICEPLUME_DET_AVERAGE
      allocate ( PLUME(ng) % rhoAm0 (Nsrc(ng), N(ng)+1) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng)+1,r8)
# endif
!
      allocate ( PLUME(ng) % lm (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % lc (Nsrc(ng), 0:N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*(N(ng)+1)*2,r8)
!
      allocate ( PLUME(ng) % ent     (Nsrc(ng), N(ng)) )
      allocate ( PLUME(ng) % det     (Nsrc(ng), N(ng)) )
      allocate ( PLUME(ng) % detI    (Nsrc(ng), N(ng)) )
      allocate ( PLUME(ng) % detFrac (Nsrc(ng), N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng)*4,r8)
# ifdef ICEPLUME_DET_NEUTRAL
!
      allocate ( PLUME(ng) % detF    (Nsrc(ng), N(ng)) )
      allocate ( PLUME(ng) % detE    (Nsrc(ng), N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng)*2,r8)
!
      allocate ( PLUME(ng) % detTrc  (Nsrc(ng), N(ng), NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng)*NT(ng),r8)
# endif
!
      allocate ( PLUME(ng) % m  (Nsrc(ng), N(ng)) )
      allocate ( PLUME(ng) % mB (Nsrc(ng), N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng)*2,r8)
!
      allocate ( PLUME(ng) % dz (Nsrc(ng), N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng),r8)
!
      allocate ( PLUME(ng) % trcAm  (Nsrc(ng), N(ng), NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng)*NT(ng),r8)
!
      allocate ( PLUME(ng) % trcB   (Nsrc(ng), NT(ng)) )
      allocate ( PLUME(ng) % trc    (Nsrc(ng), NT(ng)) )
      allocate ( PLUME(ng) % trcCum (Nsrc(ng), NT(ng)) )
      allocate ( PLUME(ng) % trcIni (Nsrc(ng), NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*NT(ng)*4,r8)
# ifdef ICEPLUME_MIX
!
      allocate ( PLUME(ng) % trsW   (Nsrc(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng),r8)
!
      allocate ( PLUME(ng) % entW   (Nsrc(ng), N(ng)) )
      allocate ( PLUME(ng) % detW   (Nsrc(ng), N(ng)) )
      allocate ( PLUME(ng) % mBW    (Nsrc(ng), N(ng)) )
      allocate ( PLUME(ng) % dzNew  (Nsrc(ng), N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng)*4,r8)
!
      allocate ( PLUME(ng) % zWNew  (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % zLNew  (Nsrc(ng), 0:N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*(N(ng)+1)*2,r8)
!
      allocate ( PLUME(ng) % dzWei  (Nsrc(ng), N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng),r8)
!
      allocate ( PLUME(ng) % trcL   (Nsrc(ng), N(ng), NT(ng)) )
      allocate ( PLUME(ng) % trcNew (Nsrc(ng), N(ng), NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng)*NT(ng)*2,r8)
# endif
# ifdef ICEPLUME_WRT_AVERAGE
!
      allocate ( PLUME(ng) % sument    (Nsrc(ng)) )
      allocate ( PLUME(ng) % summ      (Nsrc(ng)) )
      allocate ( PLUME(ng) % summB     (Nsrc(ng)) )
      allocate ( PLUME(ng) % sumtrcEnt (Nsrc(ng), NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*3,r8)
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*NT(ng),r8)
!
      allocate ( PLUME(ng) % avgent    (Nsrc(ng)) )
      allocate ( PLUME(ng) % avgm      (Nsrc(ng)) )
      allocate ( PLUME(ng) % avgmB     (Nsrc(ng)) )
      allocate ( PLUME(ng) % avgtrcEnt (Nsrc(ng), NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*3,r8)
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*NT(ng),r8)
# endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DO is=1,Nsrc(ng)
        PLUME(ng) % Isrc (is) = 0
        PLUME(ng) % Jsrc (is) = 0
        PLUME(ng) % dx   (is) = IniVal
        PLUME(ng) % dy   (is) = IniVal
        PLUME(ng) % dir  (is) = IniVal
        PLUME(ng) % trs  (is) = IniVal
# ifdef ICEPLUME_DET_AVERAGE
        PLUME(ng) % avgIstr (is) = 0
        PLUME(ng) % avgIend (is) = 0
        PLUME(ng) % avgJstr (is) = 0
        PLUME(ng) % avgJend (is) = 0
        DO k=1,N(ng)
          PLUME(ng) % avgArea (is, k) = IniVal
        ENDDO
# endif
!
        DO k=0,N(ng)
          PLUME(ng) % zW   (is, k) = IniVal
          PLUME(ng) % f    (is, k) = IniVal
          PLUME(ng) % w    (is, k) = IniVal
          PLUME(ng) % t    (is, k) = IniVal
          PLUME(ng) % s    (is, k) = IniVal
          PLUME(ng) % a    (is, k) = IniVal
          PLUME(ng) % mInt (is, k) = IniVal
          PLUME(ng) % rho  (is, k) = IniVal
        ENDDO
!
        DO k=1,N(ng)
          PLUME(ng) % zR    (is, k) = IniVal
          PLUME(ng) % tAm   (is, k) = IniVal
          PLUME(ng) % sAm   (is, k) = IniVal
          PLUME(ng) % vAm   (is, k) = IniVal
          PLUME(ng) % wAm   (is, k) = IniVal
          PLUME(ng) % tpAm  (is, k) = IniVal
          PLUME(ng) % rhoAm (is, k) = IniVal
        ENDDO
        PLUME(ng) % rhoAm (is, N(ng)+1) = IniVal
# ifdef ICEPLUME_DET_AVERAGE
        DO k=1,N(ng)
          PLUME(ng) % rhoAm0 (is, k) = IniVal
        ENDDO
        PLUME(ng) % rhoAm0 (is, N(ng)+1) = IniVal
# endif
!
        DO k=0,N(ng)
          PLUME(ng) % lm (is, k) = IniVal
          PLUME(ng) % lc (is, k) = IniVal
        ENDDO
!
        DO k=1,N(ng)
          PLUME(ng) % ent     (is, k) = IniVal
          PLUME(ng) % det     (is, k) = IniVal
          PLUME(ng) % detI    (is, k) = 0
          PLUME(ng) % detFrac (is, k) = IniVal
        ENDDO
# ifdef ICEPLUME_DET_NEUTRAL
!
        DO k=1,N(ng)
          PLUME(ng) % detF (is, k) = IniVal
          PLUME(ng) % detE (is, k) = IniVal
!
          DO itrc=1,NT(ng)
            PLUME(ng) % detTrc (is, k, itrc) = IniVal
          ENDDO
        ENDDO
# endif
!
        DO k=1,N(ng)
          PLUME(ng) % m  (is, k) = IniVal
          PLUME(ng) % mB (is, k) = IniVal
!
          PLUME(ng) % dz (is, k) = IniVal
!
          DO itrc=1,NT(ng)
            PLUME(ng) % trcAm (is, k, itrc) = IniVal
          ENDDO
        ENDDO
!
        DO itrc=1,NT(ng)
          PLUME(ng) % trcB   (is, itrc) = IniVal
          PLUME(ng) % trc    (is, itrc) = IniVal
          PLUME(ng) % trcCum (is, itrc) = IniVal
          PLUME(ng) % trcIni (is, itrc) = IniVal
        ENDDO
# ifdef ICEPLUME_MIX
!
        PLUME(ng) % trsW (is) = IniVal
!
        DO k=1,N(ng)
          PLUME(ng) % entW  (is, k) = IniVal
          PLUME(ng) % detW  (is, k) = IniVal
          PLUME(ng) % mBW   (is, k) = IniVal
          PLUME(ng) % dzNew (is, k) = IniVal
        ENDDO
!
        DO k=0,N(ng)
          PLUME(ng) % zWNew (is, k) = IniVal
          PLUME(ng) % zLNew (is, k) = IniVal
        ENDDO
!
        DO k=1,N(ng)
          PLUME(ng) % dzWei (is, k) = IniVal
!
          DO itrc=1,NT(ng)
            PLUME(ng) % trcL   (is, k, itrc) = IniVal
            PLUME(ng) % trcNew (is, k, itrc) = IniVal
          ENDDO
        ENDDO
# endif
# ifdef ICEPLUME_WRT_AVERAGE
!
        PLUME(ng) % sument (is) = IniVal
        PLUME(ng) % summ   (is) = IniVal
        PLUME(ng) % summB  (is) = IniVal
        DO itrc=1,NT(ng)
          PLUME(ng) % sumtrcEnt (is, itrc) = IniVal
        ENDDO
!
        PLUME(ng) % avgent (is) = IniVal
        PLUME(ng) % avgm   (is) = IniVal
        PLUME(ng) % avgmB  (is) = IniVal
        DO itrc=1,NT(ng)
          PLUME(ng) % avgtrcEnt (is, itrc) = IniVal
        ENDDO
# endif
      ENDDO
!
      RETURN
      END SUBROUTINE allocate_iceplume
#endif
      END MODULE mod_iceplume
