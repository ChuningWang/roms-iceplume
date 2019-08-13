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
# ifdef ICEPLUME_DETRAIN_AVERAGE
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
  real(r8), parameter :: lambda1    = -0.0573_r8
  real(r8), parameter :: lambda2    = 0.0832_r8
  real(r8), parameter :: lambda3    = 0.000761_r8
!
  real(r8), parameter :: GamT       = 0.0220_r8
  real(r8), parameter :: GamS       = 0.000620_r8
  real(r8), parameter :: Cd         = 0.065_r8
! Original value from Cowton et al. 2015
!  real(r8), parameter :: Cd         = 0.0025_r8
!
  real(r8), parameter :: RiB        = 1.0_r8
  real(r8), parameter :: gRedBkg    = 0.03_r8
  real(r8), parameter :: CdBkg      = 0.0025_r8
  real(r8), parameter :: velBkg     = 0.03_r8
  real(r8), parameter :: wIni       = 1.0_r8
!
  real(r8), parameter :: detSigma   = 0.5_r8
# ifdef ICEPLUME_DETRAIN_AVERAGE
  real(r8), parameter :: avgAlpha   = 1.0_r8
# endif
!
! =====================================================================!
!                                                                      !
! PLUME Type variables                                                 !
!                                                                      !
! dir           - direction of plume. +1 for positve direction and     !
!                 -1 for negative direction. 0 for other situation.    !
! trs           - depth integrated total transport [m^3 s^-1]          !
!                                                                      !
! Profiles                                                             !
!                                                                      !
! zW            - depth [m]                                            !
! f             - plume vertical volume flux [m^3 s^-1]                !
! w             - plume vertical velocity [m s^-1]                     !
! t             - plume temperature [degC]                             !
! s             - plume salinity [PSU]                                 !
! a             - plume area integrated [m^2]                          !
! mInt          - plume area integrated melt [m^3 s^-1]                !
! rho           - plume density [kg m^-3]                              !
!                                                                      !
! zR            - depth at Rho points [m]                              !
! sAm           - ambient salinity [PSU]                               !
! tAm           - ambient temperature [degC]                           !
! vAm           - horizontal velocity parallel to glacier [m s^-1]     !
! wAm           - vertical velocity parallel to glacier [m s^-1]       !
! tpAm          - ambient potential temperature [degC]                 !
! rhoAm         - ambient density [kg m^-3]                            !
!                                                                      !
! lm            - plume/glacier contact length [m]                     !
! lc            - plume/water contact length [m]                       !
!                                                                      !
! ent           - entrainment rate [m^3 s^-1]                          !
! det           - detrainment rate [m^3 s^-1]                          !
! detI          - detrainment flag                                     !
! detFrac       - fraction of detrainment in vertical direction        !
# ifdef ICEPLUME_DETRAIN_NEUTRAL
!                                                                      !
! detF          - detrainment rate of freshwater [m^3 s^-1]            !
! detE          - detrainment rate of entrainment [m^3 s^-1]           !
! detTrc        - detrainment tracer concentration                     !
# endif
!                                                                      !
! m             - plume melt rate [m^3 s^-1]                           !
! mB            - background melt rate [m^3 s^-1]                      !
!                                                                      !
! dz            - RHO layer thickness [m]                              !
!                                                                      !
! Passive tracers                                                      !
!                                                                      !
! trc           - tracer concentration                                 !
! trcAm         - ambient tracer concentration                         !
! trcB          - background meltwater tracer concentration            !
! trcCum        - accumulative tracer concentration                    !
! trcIni        - initial tracer concentration in discharge            !
!                                                                      !
# ifdef ICEPLUME_VIRTUAL_MIX
!                                                                      !
! For calculation of virtual tracer flux                               !
! dzNew         - new layer thickness [m]                              !
! zWNew         - new omega surface depth [m]                          !
! zLNew         - new box vertical boundary [m]                        !
! dzWei         - new layer weight function                            !
! trcL          - new box tracer concentration                         !
! trcNew        - updated tracer concentration                         !
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
    real(r8), pointer :: dir(:)
    real(r8), pointer :: trs(:)
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
    integer(r8), pointer :: detI(:, :)
    real(r8), pointer :: detFrac(:, :)
# ifdef ICEPLUME_DETRAIN_NEUTRAL
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
# ifdef ICEPLUME_VIRTUAL_MIX
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
  END TYPE T_PLUME
!
  TYPE (T_PLUME), allocatable :: PLUME(:)
!
! =====================================================================!
!                                                                      !
! Allocate PLUME Type variables.                                       !
!                                                                      !
! =====================================================================!
!
  CONTAINS
    SUBROUTINE allocate_iceplume(ng)
      integer :: ng
      IF (ng .EQ. 1) allocate( PLUME(Ngrids) )
!
! Allocate profiles
!
      allocate( PLUME(ng) % dir(1:Nsrc(ng)) )
      allocate( PLUME(ng) % trs(1:Nsrc(ng)) )
!
      allocate( PLUME(ng) % zW   (Nsrc(ng), 0:N(ng)) )
      allocate( PLUME(ng) % f    (Nsrc(ng), 0:N(ng)) )
      allocate( PLUME(ng) % w    (Nsrc(ng), 0:N(ng)) )
      allocate( PLUME(ng) % t    (Nsrc(ng), 0:N(ng)) )
      allocate( PLUME(ng) % s    (Nsrc(ng), 0:N(ng)) )
      allocate( PLUME(ng) % a    (Nsrc(ng), 0:N(ng)) )
      allocate( PLUME(ng) % mInt (Nsrc(ng), 0:N(ng)) )
      allocate( PLUME(ng) % rho  (Nsrc(ng), 0:N(ng)) )
!
      allocate( PLUME(ng) % zR    (Nsrc(ng), N(ng)  ) )
      allocate( PLUME(ng) % tAm   (Nsrc(ng), N(ng)  ) )
      allocate( PLUME(ng) % sAm   (Nsrc(ng), N(ng)  ) )
      allocate( PLUME(ng) % vAm   (Nsrc(ng), N(ng)  ) )
      allocate( PLUME(ng) % wAm   (Nsrc(ng), N(ng)  ) )
      allocate( PLUME(ng) % tpAm  (Nsrc(ng), N(ng)  ) )
      allocate( PLUME(ng) % rhoAm (Nsrc(ng), N(ng)+1) )
!
      allocate( PLUME(ng) % lm (Nsrc(ng), 0:N(ng)) )
      allocate( PLUME(ng) % lc (Nsrc(ng), 0:N(ng)) )
!
      allocate( PLUME(ng) % ent  (Nsrc(ng), N(ng)) )
      allocate( PLUME(ng) % det  (Nsrc(ng), N(ng)) )
      allocate( PLUME(ng) % detI (Nsrc(ng), N(ng)) )
      allocate( PLUME(ng) % detFrac (Nsrc(ng), N(ng)) )
# ifdef ICEPLUME_DETRAIN_NEUTRAL
!
      allocate( PLUME(ng) % detF    (Nsrc(ng), N(ng)) )
      allocate( PLUME(ng) % detE    (Nsrc(ng), N(ng)) )
      allocate( PLUME(ng) % detTrc  (Nsrc(ng), N(ng), NT(ng)) )
# endif
!
      allocate( PLUME(ng) % m        (Nsrc(ng), N(ng)) )
      allocate( PLUME(ng) % mB       (Nsrc(ng), N(ng)) )
!
      allocate( PLUME(ng) % dz (Nsrc(ng), N(ng)) )
!
      allocate( PLUME(ng) % trcAm  (Nsrc(ng), N(ng), NT(ng)) )
      allocate( PLUME(ng) % trcB   (Nsrc(ng), NT(ng)) )
      allocate( PLUME(ng) % trc    (Nsrc(ng), NT(ng)) )
      allocate( PLUME(ng) % trcCum (Nsrc(ng), NT(ng)) )
      allocate( PLUME(ng) % trcIni (Nsrc(ng), NT(ng)) )
# ifdef ICEPLUME_VIRTUAL_MIX
!
      allocate( PLUME(ng) % dzNew  (Nsrc(ng), N(ng)) )
      allocate( PLUME(ng) % zWNew  (Nsrc(ng), 0:N(ng)) )
      allocate( PLUME(ng) % zLNew  (Nsrc(ng), 0:N(ng)) )
      allocate( PLUME(ng) % dzWei  (Nsrc(ng), N(ng)) )
      allocate( PLUME(ng) % trcL   (Nsrc(ng), N(ng), NT(ng)) )
      allocate( PLUME(ng) % trcNew (Nsrc(ng), N(ng), NT(ng)) )
# endif
    END SUBROUTINE allocate_iceplume
#endif
END
