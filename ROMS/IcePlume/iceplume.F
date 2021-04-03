#include "cppdefs.h"
#ifdef ICEPLUME
!
!=======================================================================
!                                                                      !
!  This script contains -                                              !
!    (1) A function to initiate and exchange variables.                !
!    (2) A wraper for the iceplume module.                             !
!    (3) A function to spread plume fluxes in a few adjoint grids.     !
!    (4) A function (coupler) to mix the fluxes into a grid instead of !
!      applying them as boundary conditions.                           !
!                                                                      !
!=======================================================================
!                                                                      !
!  PART I - Initiate the iceplume model.                               !
!    This function read in some position variables from SOURCES, and   !
!    pass them to ICEPLUME.                                            !
!                                                                      !
!=======================================================================
!
!***********************************************************************
      SUBROUTINE iceplume_init (ng)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: tile
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3, __LINE__, __FILE__)
# endif
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL iceplume_init_tile (ng, tile)
      ENDDO
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3, __LINE__, __FILE__)
# endif

      RETURN
      END SUBROUTINE iceplume_init
!
!***********************************************************************
      SUBROUTINE iceplume_init_tile (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars, ONLY : LuvSrc, LwSrc
      USE mod_grid, ONLY : GRID
      USE mod_sources
      USE mod_iceplume
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_collect
!
      integer :: Ispv
      real(r8) :: Fspv
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: is, Isrc, Jsrc, Dsrc, D2src
# ifdef ICEPLUME_DET_AVERAGE
      integer :: i, j
# endif
!
# include "set_bounds.h"
!
# ifdef ICEPLUME_DET_AVERAGE
!
!  Check if use temporal exponential averaging.
!
      IF ((0.0_r8 .LT. avgAlpha) .AND. (avgAlpha .LT. 1.0_r8)) THEN
        avgExp = .TRUE.
      ELSE
        avgExp = .FALSE.
      ENDIF
# endif
      DO is = 1, Nsrc(ng)
        Isrc = SOURCES(ng) % Isrc(is)
        Jsrc = SOURCES(ng) % Jsrc(is)
        IF ( ((IstrR .LE. Isrc) .AND. (Isrc .LE. IendR)) .AND.          &
     &       ((JstrR .LE. Jsrc) .AND. (Jsrc .LE. JendR)) ) THEN
          Dsrc = NINT(SOURCES(ng) % Dsrc(is))
          IF (LuvSrc(ng)) THEN
            IF (Dsrc .EQ. 0) THEN
              IF ((GRID(ng)%rmask(Isrc  , Jsrc) .LT. 0.5_r8) .AND.      &
     &            (GRID(ng)%rmask(Isrc-1, Jsrc) .GT. 0.5_r8)) THEN
                Isrc = Isrc-1
                D2src = -1
              ELSEIF ((GRID(ng)%rmask(Isrc  , Jsrc) .GT. 0.5_r8) .AND.  &
     &                (GRID(ng)%rmask(Isrc-1, Jsrc) .LT. 0.5_r8)) THEN
                D2src = 1
              ELSE
                D2src = 0
              ENDIF
            ELSEIF (Dsrc .EQ. 1) THEN
              IF ((GRID(ng)%rmask(Isrc, Jsrc  ) .LT. 0.5_r8) .AND.      &
     &            (GRID(ng)%rmask(Isrc, Jsrc-1) .GT. 0.5_r8)) THEN
                Jsrc = Jsrc-1
                D2src = -1
              ELSEIF ((GRID(ng)%rmask(Isrc, Jsrc  ) .GT. 0.5_r8) .AND.  &
     &                (GRID(ng)%rmask(Isrc, Jsrc-1) .LT. 0.5_r8)) THEN
                D2src = 1
              ELSE
                D2src = 0
              ENDIF
            ENDIF
          ELSEIF (LwSrc(ng)) THEN
            D2src = 0
          ENDIF
          PLUME(ng) % Isrc(is) = Isrc
          PLUME(ng) % Jsrc(is) = Jsrc
          PLUME(ng) % dir(is) = REAL(D2src)
!
!  Get grid dimensions.
!
          IF (Dsrc .EQ. 0) THEN
            PLUME(ng) % dx(is) = ABS(GRID(ng) % om_r(Isrc, Jsrc))
            PLUME(ng) % dy(is) = ABS(GRID(ng) % on_r(Isrc, Jsrc))
          ELSE
            PLUME(ng) % dx(is) = ABS(GRID(ng) % on_r(Isrc, Jsrc))
            PLUME(ng) % dy(is) = ABS(GRID(ng) % om_r(Isrc, Jsrc))
          ENDIF
!
!  Assign surface density (air density).
!
          PLUME(ng) % rhoAm(is, N(ng)+1) = rhoAir
        ENDIF
# ifdef ICEPLUME_DET_AVERAGE
!
!  Read in iceplume averaging span.
!
        PLUME(ng) % avgIstr(is) = NINT(SOURCES(ng) % SGIrange(is, 1))
        PLUME(ng) % avgIend(is) = NINT(SOURCES(ng) % SGIrange(is, 2))
        PLUME(ng) % avgJstr(is) = NINT(SOURCES(ng) % SGJrange(is, 1))
        PLUME(ng) % avgJend(is) = NINT(SOURCES(ng) % SGJrange(is, 2))
# endif
      ENDDO
# ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes.
!
      K = N(ng)+1
      CALL mp_collect (ng,iNLM,Nsrc(ng),Ispv,PLUME(ng)%Isrc(1:Nsrc(ng)))
      CALL mp_collect (ng,iNLM,Nsrc(ng),Ispv,PLUME(ng)%Jsrc(1:Nsrc(ng)))
      CALL mp_collect (ng,iNLM,Nsrc(ng),Fspv,PLUME(ng)%dir(1:Nsrc(ng)))
      CALL mp_collect (ng,iNLM,Nsrc(ng),Fspv,PLUME(ng)%dx(1:Nsrc(ng)))
      CALL mp_collect (ng,iNLM,Nsrc(ng),Fspv,PLUME(ng)%dy(1:Nsrc(ng)))
# endif
      RETURN
      END SUBROUTINE iceplume_init_tile
!
!=======================================================================
!                                                                      !
!  PART II - Wrapper for the iceplume module.                          !
!    This function extract data from the ocean model and run the plume !
!    model.                                                            !
!                                                                      !
!=======================================================================
!
!***********************************************************************
      SUBROUTINE iceplume (ng, tile)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 4, __LINE__, __FILE__)
# endif
      CALL iceplume_tile (ng, tile)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4, __LINE__, __FILE__)
# endif

      RETURN
      END SUBROUTINE iceplume
!
!***********************************************************************
      SUBROUTINE iceplume_tile (ng, tile)
!***********************************************************************
!
      USE mod_grid, ONLY : GRID
      USE mod_ocean, ONLY : OCEAN
      USE mod_stepping, ONLY : nstp
      USE mod_sources
      USE mod_scalars, ONLY : dt, iic
      USE mod_iceplume
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_assemble
!
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: is, K, itrc, Isrc, Jsrc, Dsrc
      real(r8), parameter :: prRef = 101.d3*1.d-4
      real(r8) :: pr, rhoAm
      real(r8) :: dx, dy
      real(r8) :: RHO
      real(r8), parameter :: spv = 0.0_r8
# ifdef ICEPLUME_DET_AVERAGE
      integer :: i, j
      real(r8) :: area
#  ifdef DISTRIBUTE
      integer :: Npts
#  endif
# endif
!
# include "set_bounds.h"
!
# if defined ICEPLUME_DET_AVERAGE && defined DISTRIBUTE
!
!  Initialize arrays to facilitate collective communications.
!
      PLUME(ng) % zR = spv
# endif
!
! ==================================================================!
!                                                                   !
! Read in and prepare profiles from OCEAN.                          !
!                                                                   !
! ==================================================================!
!
      DO is = 1, Nsrc(ng)
        Isrc = PLUME(ng) % Isrc(is)
        Jsrc = PLUME(ng) % Jsrc(is)
        IF ( ((IstrR .LE. Isrc) .AND. (Isrc .LE. IendR)) .AND.          &
     &       ((JstrR .LE. Jsrc) .AND. (Jsrc .LE. JendR)) ) THEN
          Dsrc = NINT(SOURCES(ng) % Dsrc(is))
!
          DO K = 0, N(ng)
            PLUME(ng) % zW(is, K) = GRID(ng) % z_w(Isrc, Jsrc, K)
          ENDDO
!
          DO K = 1, N(ng)
            PLUME(ng) % zR(is, K) = GRID(ng) % z_r(Isrc, Jsrc, K)
            PLUME(ng) % dz(is, K) =                                     &
     &          PLUME(ng) % zW(is, K) - PLUME(ng) % zW(is, K-1)
!
!  Temp, salt.
!
            PLUME(ng) % sAm(is, K) =                                    &
     &          OCEAN(ng) % t(Isrc, Jsrc, K, nstp(ng), isalt)
            PLUME(ng) % tpAm(is, K) =                                   &
     &          OCEAN(ng) % t(Isrc, Jsrc, K, nstp(ng), itemp)
# ifndef ICEPLUME_DET_AVERAGE
!
!  Density.
!
            PLUME(ng) % rhoAm(is, K) = OCEAN(ng) % rho(Isrc, Jsrc, K)
# endif
!
!  Convert potential temp to in-situ temp.
!
            pr = prRef+(ABS(PLUME(ng) % zR(is, K))*rhoRef*g)*1.d-4
            CALL SW_TEMP(PLUME(ng) % sAm(is, K),                        &
     &                   PLUME(ng) % tpAm(is, K),                       &
     &                   pr, prRef, PLUME(ng) % tAm(is, K))
!
!  Tracers.
!
            DO itrc = 1, NT(ng)
              PLUME(ng) % trcAm(is, K, itrc) =                          &
     &            OCEAN(ng) % t(Isrc, Jsrc, K, nstp(ng), itrc)
            ENDDO
!
!  U/V, W.
!
            IF (Dsrc .EQ. 0) THEN
              PLUME(ng) % vAm(is, K) = 0.5_r8 * (                       &
     &            OCEAN(ng) % v(Isrc  , Jsrc, K, nstp(ng)) +            &
     &            OCEAN(ng) % v(Isrc+1, Jsrc, K, nstp(ng)) )
            ELSEIF (Dsrc .EQ. 1) THEN
              PLUME(ng) % vAm(is, K) = 0.5_r8 * (                       &
     &            OCEAN(ng) % u(Isrc, Jsrc  , K, nstp(ng)) +            &
     &            OCEAN(ng) % u(Isrc, Jsrc+1, K, nstp(ng)) )
            ENDIF
            PLUME(ng) % wAm(is, K) = 0.5_r8 * (                         &
     &          OCEAN(ng) % wvel(Isrc, Jsrc, K-1) +                     &
     &          OCEAN(ng) % wvel(Isrc, Jsrc, K  ) )
          ENDDO
!
!  Discharge tracer concentration.
!
          DO itrc = 1, NT(ng)
# ifdef ICEPLUME_TRACER
            PLUME(ng) % trcIni(is, itrc) = SOURCES(ng) % SGtrc(is, itrc)
# else
            PLUME(ng) % trcIni(is, itrc) = 0.0_r8
# endif
          ENDDO
        ENDIF
      ENDDO
# ifdef ICEPLUME_DET_AVERAGE
!
! ==================================================================!
!                                                                   !
! Average density profiles in a box to surpress oscillation.        !
!                                                                   !
! ==================================================================!
!
#  ifdef DISTRIBUTE
      Npts = Nsrc(ng)*N(ng)
!
!  Exchange data between all parallel nodes.
!
      CALL mp_assemble (ng, iNLM, Npts, spv,                            &
     &                  PLUME(ng) % zR(1:Nsrc(ng), 1:N(ng)))
#  endif
!
      DO is = 1, Nsrc(ng)
        IF ( NINT(SOURCES(ng) % SGtyp(is)) .GT. 1 ) THEN
          DO K = 1, N(ng)
            PLUME(ng) % rhoAm0(is, K) = PLUME(ng) % rhoAm(is, K)
            PLUME(ng) % rhoAm(is, K) = spv
            PLUME(ng) % avgArea(is, K) = spv
          ENDDO
          DO j = PLUME(ng) % avgJstr(is), PLUME(ng) % avgJend(is)
            DO i = PLUME(ng) % avgIstr(is), PLUME(ng) % avgIend(is)
              IF ( ((IstrR .LE. i) .AND. (i .LE. IendR)) .AND.          &
      &            ((JstrR .LE. j) .AND. (j .LE. JendR)) .AND.          &
      &            (GRID(ng) % rmask(i, j) .GT. 0.5_r8) ) THEN
                area = ABS(GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j))
                DO K = 1, N(ng)
                  IF ( (PLUME(ng) % zR(is, K) .GT.                      &
      &                 GRID(ng) % z_r(i, j, 1)) ) THEN
!
!  Interpolate density to the right depth.
!
                    CALL LININT(N(ng),                                  &
     & GRID(ng) % z_r(i, j, :), OCEAN(ng) % rho(i, j, :),               &
     & PLUME(ng) % zR(is, K),   rhoAm)
!
!  Integrate density and averaging area.
!
                    PLUME(ng) % avgArea(is, K) =                        &
     &                PLUME(ng) % avgArea(is, K) + area
                    PLUME(ng) % rhoAm(is, K) =                          &
     &                PLUME(ng) % rhoAm(is, K) + rhoAm*area
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes.
!
      CALL mp_assemble (ng, iNLM, Npts, spv,                            &
     &                  PLUME(ng) % rhoAm(1:Nsrc(ng), 1:N(ng)))
      CALL mp_assemble (ng, iNLM, Npts, spv,                            &
     &                  PLUME(ng) % avgArea(1:Nsrc(ng), 1:N(ng)))
#  endif
!
      DO is = 1, Nsrc(ng)
        Isrc = PLUME(ng) % Isrc(is)
        Jsrc = PLUME(ng) % Jsrc(is)
        IF ( ((IstrR .LE. Isrc) .AND. (Isrc .LE. IendR)) .AND.          &
     &       ((JstrR .LE. Jsrc) .AND. (Jsrc .LE. JendR)) ) THEN
          DO K = 1, N(ng)
            PLUME(ng) % rhoAm(is, K) =                                  &
     &        PLUME(ng) % rhoAm(is, K) / PLUME(ng) % avgArea(is, K)
!
!  Convert from density anomaly to absolute density.
!
            PLUME(ng) % rhoAm(is, K) =                                  &
     &        PLUME(ng) % rhoAm(is, K)+1000.0_r8
!
!  Exponentially average over time.
!
            IF ( (avgExp) ) THEN
              IF ( iic(ng) .GT. ntstart(ng) ) THEN
                PLUME(ng) % rhoAm(is, K) =                              &
     &            PLUME(ng) % rhoAm0(is, K)*(1.0-avgAlpha) +            &
     &            PLUME(ng) % rhoAm(is, K) *avgAlpha
               ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
# endif
!
! ==================================================================!
!                                                                   !
! Call the main function here.                                      !
!                                                                   !
! ==================================================================!
!
      DO is = 1, Nsrc(ng)
        Isrc = PLUME(ng) % Isrc(is)
        Jsrc = PLUME(ng) % Jsrc(is)
        IF ( ((IstrR .LE. Isrc) .AND. (Isrc .LE. IendR)) .AND.          &
     &       ((JstrR .LE. Jsrc) .AND. (Jsrc .LE. JendR)) ) THEN
          IF (NINT(SOURCES(ng) % SGtyp(is)) .NE. 0) THEN
            CALL iceplume_calc(ng, is,                                  &
     &                         ABS(SOURCES(ng) % SGbar(is)),            &
     &                         SOURCES(ng) % SGtrc(is, itemp),          &
     &                         SOURCES(ng) % SGtrc(is, isalt),          &
     &                         NINT(SOURCES(ng) % SGtyp(is)),           &
     &                         SOURCES(ng) % SGdep(is),                 &
     &                         SOURCES(ng) % SGlen(is))
          ENDIF
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE iceplume_tile
# ifdef ICEPLUME_SPREAD
!
!=======================================================================
!                                                                      !
!  PART III - Spread the plume horizontally in a few grids.            !
!    This section is still being tested.                               !
!                                                                      !
!=======================================================================
!
!***********************************************************************
      SUBROUTINE iceplume_spread (ng, tile)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 4, __LINE__, __FILE__)
# endif
      CALL iceplume_spread_tile (ng, tile)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4, __LINE__, __FILE__)
# endif

      RETURN
      END SUBROUTINE iceplume_spread
!
!***********************************************************************
      SUBROUTINE iceplume_spread_tile (ng, tile)
!***********************************************************************
!
      USE mod_ocean, ONLY : OCEAN
      USE mod_stepping, ONLY : nstp
      USE mod_sources
      USE mod_iceplume
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: is, is0, is1, Isrc, Jsrc, K, itrc
      real(r8) :: fIni, trcIni, trcM
      real(r8) :: cff, cff1
!
# include "set_bounds.h"
!
      is  = 2
      Isrc = PLUME(ng) % Isrc(is)
      Jsrc = PLUME(ng) % Jsrc(is)
      IF ( ((IstrR .LE. Isrc) .AND. (Isrc .LE. IendR)) .AND.            &
     &     ((JstrR .LE. Jsrc) .AND. (Jsrc .LE. JendR)) ) THEN
        is0 = is-1
        is1 = is+1
        cff = PLUME(ng) % trs(is)/3.0_r8
        PLUME(ng) % trs(is0) = cff
        PLUME(ng) % trs(is1) = cff
        PLUME(ng) % trs(is ) = cff
        DO K = 1, N(ng)
          cff = PLUME(ng) % ent(is, K)/3.0_r8
          PLUME(ng) % ent(is0, K) = cff
          PLUME(ng) % ent(is1, K) = cff
          PLUME(ng) % ent(is , K) = cff
!
          cff = PLUME(ng) % det(is, K)/3.0_r8
          PLUME(ng) % det(is0, K) = cff
          PLUME(ng) % det(is1, K) = cff
          PLUME(ng) % det(is , K) = cff
!
          cff = PLUME(ng) % m(is, K)/3.0_r8
          PLUME(ng) % m(is0, K) = cff
          PLUME(ng) % m(is1, K) = cff
          PLUME(ng) % m(is , K) = cff
        ENDDO
!
!  Tracers.
!
        fIni = ABS(SOURCES(ng) % SGbar(is))/3.0_r8
        DO itrc = 1, NT(ng)
          trcIni = PLUME(ng) % trcIni(is, itrc)
          IF (itrc .EQ. isalt) THEN
            trcM = sIce
          ELSEIF (itrc .EQ. itemp) THEN
            trcM = tGade
          ELSE
            trcM = 0.0_r8
          ENDIF
# ifdef ICEPLUME_MELT_TRACER
          IF (itrc .EQ. NT(ng)-1) THEN
            trcIni = 0.0_r8
            trcM   = PLUME(ng) % trcIni(is, itrc)
          ELSEIF (itrc .EQ. NT(ng)) THEN
            trcIni = 0.0_r8
            trcM   = 0.0_r8
          ENDIF
# endif
          cff  = fIni
          cff1 = fIni*trcIni
          DO K = 1, N(ng)
            cff  = cff -PLUME(ng) % ent(is, K)+PLUME(ng) % m(is, K)
            cff1 = cff1-                                                &
     &        PLUME(ng) % ent(is, K)*PLUME(ng) % trcAm(is, K, itrc)+    &
     &        PLUME(ng) % m(is, K)*trcM
          ENDDO
          PLUME(ng) % trc(is, itrc) = cff1/cff
!
          cff  = fIni
          cff1 = fIni*trcIni
          DO K = 1, N(ng)
            cff  = cff -PLUME(ng) % ent(is0, K)+PLUME(ng) % m(is0, K)
            cff1 = cff1-                                                &
     &        PLUME(ng) % ent(is0, K)*PLUME(ng) % trcAm(is0, K, itrc)+  &
     &        PLUME(ng) % m(is0, K)*trcM
          ENDDO
          PLUME(ng) % trc(is1, itrc) = cff1/cff
!
          cff  = fIni
          cff1 = fIni*trcIni
          DO K = 1, N(ng)
            cff  = cff -PLUME(ng) % ent(is1, K)+PLUME(ng) % m(is1, K)
            cff1 = cff1-                                                &
     &        PLUME(ng) % ent(is1, K)*PLUME(ng) % trcAm(is1, K, itrc)+  &
     &        PLUME(ng) % m(is1, K)*trcM
          ENDDO
          PLUME(ng) % trc(is1, itrc) = cff1/cff
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE iceplume_spread_tile
# endif
# ifdef ICEPLUME_MIX
!
!=======================================================================
!                                                                      !
!  PART IV - Virtually mix the plume transport into the ocean grid.    !
!                                                                      !
!=======================================================================
!
!***********************************************************************
      SUBROUTINE iceplume_mix (ng, tile)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 4, __LINE__, __FILE__)
# endif
      CALL iceplume_mix_tile (ng, tile)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4, __LINE__, __FILE__)
# endif

      RETURN
      END SUBROUTINE iceplume_mix
!
!***********************************************************************
      SUBROUTINE iceplume_mix_tile (ng, tile)
!***********************************************************************
!
      USE mod_grid, ONLY : GRID
      USE mod_ocean, ONLY : OCEAN
      USE mod_stepping, ONLY : nstp, kstp
      USE mod_sources
      USE mod_scalars, ONLY : dt, Vtransform, hc, SCALARS
      USE mod_iceplume
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: is, Isrc, Jsrc, Dsrc, D2src
      integer :: K, K2, itrc
      real(r8) :: area
      real(r8) :: cff, cff1, cff2, cff3
!
!  From grid construction functions.
!
      real(r8) :: cff_r, cff1_r, cff2_r, cff_w, cff1_w, cff2_w
      real(r8) :: hinv, hwater, z_r0, z_w0
      real(r8) :: C2_r, C2_w, hh2, vert_n1, vert_a, vert_h0, vert_s0
!
# include "set_bounds.h"
!
      DO is = 1, Nsrc(ng)
        Isrc = PLUME(ng) % Isrc(is)
        Jsrc = PLUME(ng) % Jsrc(is)
        IF ( ((IstrR .LE. Isrc) .AND. (Isrc .LE. IendR)) .AND.          &
     &       ((JstrR .LE. Jsrc) .AND. (Jsrc .LE. JendR)) ) THEN
#  ifdef ICEPLUME_SLOPE
!
! ==================================================================!
!                                                                   !
! Separate the detrainment in horizontal and vertical direction.    !
!                                                                   !
! ==================================================================!
!
          cff = COS(pi*SOURCES(ng) % SGang(is)/180.0_r8)
#  else
          cff = 0.0_r8
#  endif
          cff1 = 1.0_r8 - cff
          PLUME(ng) % trsW(is) = PLUME(ng) % trs(is)*cff1
          PLUME(ng) % trs(is)  = PLUME(ng) % trs(is)*cff
          DO K = 1, N(ng)
            PLUME(ng) % entW(is, K) = PLUME(ng) % ent(is, K)*cff1
            PLUME(ng) % detW(is, K) = PLUME(ng) % det(is, K)*cff1
            PLUME(ng) % mBW(is, K)  = PLUME(ng) % mB(is, K)*cff1
            PLUME(ng) % ent(is, K)  = PLUME(ng) % ent(is, K)*cff
            PLUME(ng) % det(is, K)  = PLUME(ng) % det(is, K)*cff
            PLUME(ng) % mB(is, K)   = PLUME(ng) % mB(is, K)*cff
          ENDDO
!
! ==================================================================!
!                                                                   !
! Calculate virtual vertical mixing if grid is not adjoint to land. !
!                                                                   !
! ==================================================================!
!
! Calculate new box vertical boundary.
!
          area = PLUME(ng) % dx(is)*PLUME(ng) % dy(is)
          PLUME(ng) % zLNew(is, 0) = PLUME(ng) % zW(is, 0)
          DO K = 1, N(ng)
            PLUME(ng) % dzNew(is, K) = PLUME(ng) % dz(is, K) +          &
     &                 (PLUME(ng) % detW(is, K) +                       &
     &                  PLUME(ng) % entW(is, K) +                       &
     &                  PLUME(ng) % mBW(is, K))*dt(ng)/area
            PLUME(ng) % zLNew(is, K) =                                  &
     &          PLUME(ng) % zLNew(is, K-1) +                            &
     &          PLUME(ng) % dzNew(is, K)
          ENDDO
!
!  Calculate tracer concentration.
!
          DO K = 1, N(ng)
            cff  = area*PLUME(ng) % dz(is, K)
            cff1 = (PLUME(ng) % detW(is, K) +                           &
     &              PLUME(ng) % entW(is, K) +                           &
     &              PLUME(ng) % mBW(is, K))*dt(ng)
            DO itrc = 1, NT(ng)
              cff2 = (PLUME(ng) % detW(is, K)*                          &
#  ifdef ICEPLUME_DET_NEUTRAL
     &                PLUME(ng) % detTrc(is, K, itrc)+                  &
#  else
     &                PLUME(ng) % trc(is, itrc)+                        &
#  endif
     &                PLUME(ng) % entW(is, K)*                          &
     &                PLUME(ng) % trcAm(is, K, itrc)+                   &
     &                PLUME(ng) % mBW(is, K)*                           &
     &                PLUME(ng) % trcB(is, K, itrc))*dt(ng)
              PLUME(ng) % trcL(is, K, itrc) =                           &
     &            (PLUME(ng) % trcAm(is, K, itrc)*cff + cff2) /         &
     &            (cff + cff1)
            ENDDO
          ENDDO
!
!  Calculate new W-point location. This is copied from set_depth.
!
          cff = OCEAN(ng) % zeta(Isrc, Jsrc, kstp(ng))+                 &
     &          PLUME(ng)%trsW(is)*dt(ng)/area
!
          IF (Vtransform(ng).eq.3) THEN
            vert_s0 = 90./120.
            vert_n1 = 2./3. !3./5.
            vert_h0 = 60./(vert_s0)**(vert_n1)
            vert_a = 1./(vert_s0-1)
          END IF
!
          PLUME(ng) % zWNew(is, 0) = PLUME(ng) % zW(is, 0)
          hwater = GRID(ng) % h(Isrc, Jsrc)
          IF (Vtransform(ng).eq.1) THEN
            hinv = 1.0_r8/hwater
            DO k=1,N(ng)
              cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
              cff1_w=SCALARS(ng)%Cs_w(k)
              z_w0=cff_w+cff1_w*hwater
              PLUME(ng) % zWNew(is, k)=z_w0+cff*(1.0_r8+z_w0*hinv)
            ENDDO
          ELSE IF (Vtransform(ng).eq.3) THEN
            hinv=1.0_r8/hwater
            hh2=(min(vert_h0,hwater))*hinv
            DO k=1,N(ng)
              IF (SCALARS(ng)%sc_w(k).gt.-vert_s0) THEN 
                C2_w=-hh2*(-SCALARS(ng)%sc_w(k))**(vert_n1) 
                cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-C2_w)
                cff1_w=C2_w
                z_w0=cff_w+cff1_w*hwater
                PLUME(ng) % zWNew(is, k)=z_w0+cff*(1.0_r8+z_w0*hinv)
              ELSE
                C2_w=-hh2*(-SCALARS(ng)%sc_w(k))**(vert_n1)             &
     &            -(1-hh2)*(vert_a*(SCALARS(ng)%sc_w(k)+vert_s0))**2
                cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-C2_w)
                cff1_w=C2_w
                PLUME(ng) % zWNew(is, k)=z_w0+cff*(1.0_r8+z_w0*hinv)
              ENDIF
            ENDDO
          ELSE IF (Vtransform(ng).eq.2) THEN
            hinv=1.0_r8/(hc(ng)+hwater)
            DO k=1,N(ng)
              cff_w=hc(ng)*SCALARS(ng)%sc_w(k)
              cff1_w=SCALARS(ng)%Cs_w(k)
              cff2_w=(cff_w+cff1_w*hwater)*hinv
              PLUME(ng) % zWNew(is, k)=cff+(cff+hwater)*cff2_w
            ENDDO
          ENDIF
!
!  Calculate weight function.
!
          DO K = 1, N(ng)
            cff = PLUME(ng) % zWNew(is, K) -                            &
     &            PLUME(ng) % zWNew(is, K-1)
            DO K2 = 1, N(ng)
              cff1 = PLUME(ng) % zWNew(is, K) -                         &
     &               PLUME(ng) % zLNew(is, K2)
              cff1 = MIN(MAX(cff1, 0.0_r8), cff)
              cff2 = PLUME(ng) % zWNew(is, K-1) -                       &
     &               PLUME(ng) % zLNew(is, K2-1)
              cff2 = MIN(-MIN(cff2, 0.0_r8), cff)
              PLUME(ng) % dzWei(is, K2) = (cff-cff1-cff2)/cff
            ENDDO
            DO itrc = 1, NT(ng)
              PLUME(ng) % trcNew(is, K, itrc) = 0.0_r8
              DO K2 = 1, N(ng)
                PLUME(ng) % trcNew(is, K, itrc) =                       &
     &              PLUME(ng) % trcNew(is, K, itrc) +                   &
     &              PLUME(ng) % dzWei(is, K2)*                          &
     &              PLUME(ng) % trcL(is, K2, itrc)
              ENDDO
            ENDDO
          ENDDO
!
!  Substitute new zeta values (2-D field) back to OCEAN.
!
          OCEAN(ng) % zeta(Isrc, Jsrc, kstp(ng)) =                      &
     &      OCEAN(ng) % zeta(Isrc, Jsrc, kstp(ng))+                     &
     &      PLUME(ng)%trsW(is)*dt(ng)/area
!
!  Substitute new tracer values (3-D field) back to OCEAN.
!
          DO K = 1, N(ng)
            DO itrc = 1, NT(ng)
              OCEAN(ng) % t(Isrc, Jsrc, K, nstp(ng), itrc) =            &
     &            PLUME(ng) % trcNew(is, K, itrc)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE iceplume_mix_tile
# endif
#endif
