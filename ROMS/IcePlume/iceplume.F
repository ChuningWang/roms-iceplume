#include "cppdefs.h"
      MODULE iceplume_mod
#ifdef ICEPLUME
!
      implicit none
!
      PRIVATE
      PUBLIC  :: iceplume_init, iceplume
# ifdef ICEPLUME_MIX
      PUBLIC  :: iceplume_mix
# endif
# ifdef ICEPLUME_SPREAD
      PUBLIC  :: iceplume_spread
# endif
!
      CONTAINS
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
# ifdef MASKING
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
# else
          D2src = 0
# endif
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
!=======================================================================
!                                                                      !
! Read in and prepare profiles from OCEAN.                             !
!                                                                      !
!=======================================================================
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
            CALL sw_temp(PLUME(ng) % sAm(is, K),                        &
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
!=======================================================================
!                                                                      !
! Average density profiles in a box to surpress oscillation.           !
!                                                                      !
!=======================================================================
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
#  ifdef MASKING
              IF ( ((IstrR .LE. i) .AND. (i .LE. IendR)) .AND.          &
      &            ((JstrR .LE. j) .AND. (j .LE. JendR)) .AND.          &
      &            (GRID(ng) % rmask(i, j) .GT. 0.5_r8) ) THEN
#  else
              IF ( ((IstrR .LE. i) .AND. (i .LE. IendR)) .AND.          &
      &            ((JstrR .LE. j) .AND. (j .LE. JendR)) ) THEN
#  endif
                area = ABS(GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j))
                DO K = 1, N(ng)
                  IF ( (PLUME(ng) % zR(is, K) .GT.                      &
      &                 GRID(ng) % z_r(i, j, 1)) ) THEN
!
!  Interpolate density to the right depth.
!
                    CALL linint(N(ng),                                  &
     &  GRID(ng) % z_r(i, j, :), OCEAN(ng) % rho(i, j, :),              &
     &  PLUME(ng) % zR(is, K),   rhoAm)
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
!=======================================================================
!                                                                      !
! Call the main function here.                                         !
!                                                                      !
!=======================================================================
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
!=======================================================================
!                                                                      !
! Separate the detrainment in horizontal and vertical direction.       !
!                                                                      !
!=======================================================================
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
!=======================================================================
!                                                                      !
! Calculate virtual vertical mixing if grid is not adjoint to land.    !
!                                                                      !
!=======================================================================
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
            cff  = cff  - PLUME(ng) % ent(is, K)+PLUME(ng) % m(is, K)
            cff1 = cff1 -                                               &
     &        PLUME(ng) % ent(is, K)*PLUME(ng) % trcAm(is, K, itrc) +   &
     &        PLUME(ng) % m(is, K)*trcM
          ENDDO
          PLUME(ng) % trc(is, itrc) = cff1/cff
!
          cff  = fIni
          cff1 = fIni*trcIni
          DO K = 1, N(ng)
            cff  = cff  - PLUME(ng) % ent(is0, K)+PLUME(ng) % m(is0, K)
            cff1 = cff1 -                                               &
     &        PLUME(ng) % ent(is0, K)*PLUME(ng) % trcAm(is0, K, itrc) + &
     &        PLUME(ng) % m(is0, K)*trcM
          ENDDO
          PLUME(ng) % trc(is1, itrc) = cff1/cff
!
          cff  = fIni
          cff1 = fIni*trcIni
          DO K = 1, N(ng)
            cff  = cff  - PLUME(ng) % ent(is1, K)+PLUME(ng) % m(is1, K)
            cff1 = cff1 -                                               &
     &        PLUME(ng) % ent(is1, K)*PLUME(ng) % trcAm(is1, K, itrc) + &
     &        PLUME(ng) % m(is1, K)*trcM
          ENDDO
          PLUME(ng) % trc(is1, itrc) = cff1/cff
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE iceplume_spread_tile
# endif
!
!=======================================================================
!                                                                      !
! This subroutine call the iceplume core functions to calcualte plume  !
!   status/meltrates/etc.                                              ! 
!                                                                      !
!=======================================================================
!
      SUBROUTINE iceplume_calc(ng, I,                                   &
     &                         fIni, tIni, sIni0,                       &
     &                         sgTyp, sgDep, sgLen)
!
      USE mod_iceplume
!
!=======================================================================
!                                                                      !
! Input variables:                                                     !
!                                                                      !
!=======================================================================
!                                                                      !
! ng      - grid identifier                                            !
! I       - point source index identifier                              !
! fIni    - subglacial runoff volume flux [m^3 s^-1]                   !
! tIni    - subglacial runoff temperature [degC]                       !
! sIni    - subglacial runoff salinity [PSU]                           !
! sgTyp   - subglacial runoff type identifier                          !
! sgDep   - subglacial runoff inject depth [m]                         !
! sgLen   - subglacial runoff discharge length [m]                     !
!                                                                      !
!=======================================================================
!                                                                      !
! Local variables:                                                     !
!                                                                      !
!=======================================================================
!                                                                      !
! iceDepthK   - vertical grid index of grounding line depth            !
! plumeDepthK - vertical grid index of plume depth                     !
! osDepthK    - vertical grid index of overshoot depth                 !
! areaP       - surface area of plume in contact with ice in           !
!               that cell [m^2]                                        !
! areaB       - surface area of plume in contact with ice out          !
!               that cell [m^2]                                        !
! meanVel     - ice tangental velocity [m s^-1]                        !
! det         - total detrainment volume [m^3 s^-1]                    !
! mB          - background meltrate [m s^-1]                           !
!                                                                      !
!=======================================================================
!
! In/out variables
!
      integer, intent(in)  :: ng, I
      real(r8), intent(in) :: fIni, tIni, sIni0
      integer, intent(in)  :: sgTyp
      real(r8), intent(in) :: sgDep, sgLen
!
! Local variables declaration
!
      integer :: iceDepthK, plumeDepthK, osDepthK
      real(r8) :: sIni
      real(r8) :: areaP, areaB, meanVel, det, mB, sB
!
      integer :: K, itrc
      real(r8) :: cff, cff1, cff2
!
! If discharge input salinity is less or equal than zero, set it to
! a small number
!
      IF (sIni0 .LE. 0.0) THEN
        sIni = 0.001
      ELSE
        sIni = sIni0
      ENDIF
!
!=======================================================================
!                                                                      !
! Find iceDepthK                                                       !
!                                                                      !
!=======================================================================
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
!=======================================================================
!                                                                      !
! Use entrainment plume model to calculate T, S, & flux in plume       !
!                                                                      !
!=======================================================================
!
! Initialize and prepare profiles
!
      plumeDepthK = iceDepthK+1
      det = 0.0
!
      DO k = 0, N(ng)
        PLUME(ng) % f(I, K)    = 0.0
        PLUME(ng) % w(I, K)    = 0.0
        PLUME(ng) % t(I, K)    = 0.0
        PLUME(ng) % s(I, K)    = 0.0
        PLUME(ng) % a(I, K)    = 0.0
        PLUME(ng) % mInt(I, K) = 0.0
        PLUME(ng) % rho(I, K)  = 0.0
        PLUME(ng) % lm(I, K)   = 0.0
        PLUME(ng) % lc(I, K)   = 0.0
      ENDDO
!
      DO k = 1, N(ng)
        PLUME(ng) % ent(I, K) = 0.0
        PLUME(ng) % det(I, K) = 0.0
!
        PLUME(ng) % detI(I, K)    = 0
        PLUME(ng) % detFrac(I, K) = 0.0
      ENDDO
!
      IF ( (fIni .GT. 0.0) .AND. (sgTyp .NE. 1) ) THEN
!
! Run the plume model
!
        CALL iceplume_entrain(ng, I, iceDepthK,                         &
     &                        fIni, tIni, sIni,                         &
     &                        sgTyp, sgLen)
!
! Calculate volume flux differential to give entrainment/extrainment
! First clear ent/det
!
        DO K = iceDepthK+1, N(ng)
          PLUME(ng) % ent(I, K) =                                       &
     &        -(PLUME(ng) % f(I, K) - PLUME(ng) % f(I, K-1))
          PLUME(ng) % ent(I, K) = PLUME(ng) % ent(I, K) +               &
     &        MAX(PLUME(ng) % mInt(I, K) -                              &
     &            PLUME(ng) % mInt(I, K-1), 0.0)
          PLUME(ng) % ent(I, K) = MIN(PLUME(ng) % ent(I, K), 0.0)
        ENDDO
!
! Separate entrainment/detrainment
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
# if defined ICEPLUME_SHOOT_NO && defined ICEPLUME_DET_NO
        plumeDepthK = osDepthK
# else
        DO K = osDepthK, iceDepthK+1, -1
          cff = rho(PLUME(ng) % t(I, osDepthK),                         &
     &              PLUME(ng) % s(I, osDepthK),                         &
     &              PLUME(ng) % zR(I, K))
          IF (K .EQ. 1) THEN
            cff1 = PLUME(ng) % rhoAm(I, 1)
            cff2 = 0.5*(PLUME(ng) % rhoAm(I, 1)+                        &
     &                  PLUME(ng) % rhoAm(I, 2))
          ELSEIF (K .EQ. N(ng)) THEN
            cff1 = 0.5*(PLUME(ng) % rhoAm(I, N(ng)-1)+                  &
     &                  PLUME(ng) % rhoAm(I, N(ng)))
            cff2 = PLUME(ng) % rhoAm(I, N(ng))
          ELSE
            cff1 = 0.5*(PLUME(ng)%rhoAm(I, K-1)+PLUME(ng)%rhoAm(I, K))
            cff2 = 0.5*(PLUME(ng)%rhoAm(I, K)+PLUME(ng)%rhoAm(I, K+1))
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
!
!=======================================================================
!                                                                      !
! Use a detrainment model to distribute the detrainment plume in       !
! several layers. Details in iceplume_detrain.F.                       !
!                                                                      !
! *** This part is completely rewritten since Feb 20, 2019.            !
!                                                                      !
!=======================================================================
!
# ifdef ICEPLUME_DET_NO
        PLUME(ng) % det(I, plumeDepthK) = det
        PLUME(ng) % detI(I, plumeDepthK) = 1
        PLUME(ng) % detFrac(I, plumeDepthK) = 1.0
# else
        CALL iceplume_detrain(ng, I,                                    &
     &                        iceDepthK, plumeDepthK, osDepthK,         &
     &                        PLUME(ng) % lc(I, osDepthK),              &
     &                        det)
# endif
      ENDIF
# ifdef ICEPLUME_MELT
!
!=======================================================================
!                                                                      !
! Calculate melt rates.                                                !
!                                                                      !
!=======================================================================
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
!=======================================================================
!                                                                      !
! If there is a plume in that cell, then need to calculate plume       !
! meltrate distinct to background melt rate. Plume melt rate is        !
! already encorporated in the plume model, and taken into account      !
! in the temperature and salinity of the plume outflow. It is  useful  !
! though to have it available as a diagnostic.                         !
!                                                                      !
!=======================================================================
!
        areaP = PLUME(ng) % a(I, K) - PLUME(ng) % a(I, K-1)
        areaP = MAX(areaP, 0.0)
        areaB = PLUME(ng) % dz(I, K)*PLUME(ng) % dy(I) - areaP
        areaB = MAX(areaB, 0.0)
!
        IF ( areaP .GT. 0.0 ) THEN
          PLUME(ng) % m(I, K) = PLUME(ng) % mInt(I, K) -                &
     &                          PLUME(ng) % mInt(I, K-1)
        ENDIF
!
!=======================================================================
!                                                                      !
! Calculate the background melt rate (not generated by plumes).        !
! This will then be used to update the temperature and salinity in     !
! the adjacent cells. Velocities are calculated at cell faces -        !
! find averages for cell centres. Does not include velocity            !
! perpendicular to ice-this differs depending on orientation of        !
! ice front.                                                           !
!                                                                      !
!=======================================================================
!
        IF ( areaB .GT. 0.0 )THEN
!
          meanVel = SQRT((PLUME(ng) % vAm(I, K))**2 +                   &
     &                   (PLUME(ng) % wAm(I, K))**2)
          CALL iceplume_meltrate(PLUME(ng) % tAm(I, K),                 &
     &                           PLUME(ng) % sAm(I, K),                 &
     &                           meanVel,                               &
     &                           PLUME(ng) % zR(I, K),                  &
     &                           mB, PLUME(ng) % tB(I, K), sB)
!
! Calcate integrated background melt rate [m3 s-1]
!
          PLUME(ng) % mB(I, K) = mB*areaB
        ENDIF
      ENDDO
# endif
!
!=======================================================================
!                                                                      !
! Calculate active and passive tracers                                 !
!                                                                      !
! Update log                                                           !
! Add active tracers, T and S into this scheme. This is necessary      !
! when distributing detrainment in several layers, since writting      !
! T & S in PLUME directly is in general not a good idea.               !
! 2018/4/20 Chuning Wang                                               !
!                                                                      !
!=======================================================================
!
! Clear local plume tracer variables
!
      DO itrc = 1, NT(ng)
        PLUME(ng) % trc(I, itrc)    = 0.0
        PLUME(ng) % trcCum(I, itrc) = 0.0
!
        DO K = iceDepthK+1, N(ng)
          PLUME(ng) % trcB(I, K, itrc) = 0.0
        ENDDO
      ENDDO
!
      IF ( (fIni .GT. 0.0) .AND. (sgTyp .NE. 1) ) THEN
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
          PLUME(ng)%trcCum(I, itrc) = PLUME(ng)%trcIni(I, itrc)*fIni
        ENDDO
!
! Add up total sum of each tracer in plume
!
        DO K = iceDepthK+1, osDepthK
          IF (PLUME(ng) % ent(I, K) .LT. 0.0) THEN
            DO itrc = 3, NT(ng)
              PLUME(ng) % trcCum(I, itrc) =                             &
     &            PLUME(ng) % trcCum(I, itrc) +                         &
     &            (-PLUME(ng) % ent(I, K) *                             &
     &              PLUME(ng) % trcAm(I, K, itrc))
            ENDDO
          ENDIF
        ENDDO
!
! Calculate concentration of tracer in outflow 
!
        DO itrc = 3, NT(ng)
          PLUME(ng)%trc(I, itrc) = PLUME(ng)%trcCum(I, itrc)/det
        ENDDO
#  ifdef ICEPLUME_MELT_TRACER
!
! If use melt water tracer, rewrite tracer concentration
! Calculate accumulated trcer amount from base to detrain depth
!
        PLUME(ng) % trcCum(I, NT(ng)-1) = 0.0
        PLUME(ng) % trcCum(I, NT(ng)  ) = 0.0
        DO K = iceDepthK+1, osDepthK
          PLUME(ng) % trcCum(I, NT(ng)-1) =                             &
     & PLUME(ng) % trcCum(I, NT(ng)-1) +                                &
     & (-PLUME(ng) % ent(I, K)*PLUME(ng) % trcAm(I, K, NT(ng)-1)) +     &
     & PLUME(ng) % m(I, K)*PLUME(ng) % trcIni(I, NT(ng)-1)
           PLUME(ng) % trcCum(I, NT(ng)  ) =                            &
     & PLUME(ng) % trcCum(I, NT(ng)  ) +                                &
     & (-PLUME(ng) % ent(I, K)*PLUME(ng) % trcAm(I, K, NT(ng)  ))
        ENDDO
!
! Calculate tracer concentration
!
        PLUME(ng)%trc(I,NT(ng)-1) = PLUME(ng)%trcCum(I,NT(ng)-1)/det
        PLUME(ng)%trc(I,NT(ng)  ) = PLUME(ng)%trcCum(I,NT(ng)  )/det
#  endif
# endif
      ENDIF  ! (fIni .GT. 0.0) .AND. (sgTyp .NE. 1)
# ifdef ICEPLUME_MELT
!
! Rewrite tracer concentration in background melt water 
!
      DO K = iceDepthK+1, N(ng)
        PLUME(ng) % trcB(I, K, itemp) =                                 &
     &     (cI*(tIce-PLUME(ng) % tB(I, K)) - L)/cW
        PLUME(ng) % trcB(I, K, isalt) = sIce
#  ifdef ICEPLUME_MELT_TRACER
!
! If use melt water tracer, rewrite the last tracer type.
!
        PLUME(ng) % trcB(I, K, NT(ng)) = PLUME(ng) % trcIni(I, NT(ng))
#  endif
      ENDDO
# endif
!
! Calculate depth integrated volume transport [m3 s-1]
!
      PLUME(ng) % trs(I) = 0.0
      DO K = 1, N(ng)
        PLUME(ng) % trs(I) = PLUME(ng) % trs(I) +                       &
     &     PLUME(ng) % det(I, K) + PLUME(ng) % ent(I, K) +              &
     &     PLUME(ng) % mB(I, K)
      ENDDO
END SUBROUTINE iceplume_calc
!
!=======================================================================
!                                                                      !
! This function calculates background melt rate.                       !
!                                                                      !
!=======================================================================
!
      SUBROUTINE iceplume_meltrate(temp, salt, vel, depth, mdot, sB, tB)
!
      USE mod_iceplume
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
      IF ( ABS(vel) .LT. velBkg )  vel = velBkg
!
! Calculate melt rate from 3 equation formualtion (as for plume
! models)
! Equations for Sb, Tb and mdot
!
      a = lambda1*(GamT*cW-GamS*cI)
      b = GamS*cI*(lambda1*salt-lambda2-lambda3*depth +                 &
     &             tIce - (L/cI)) -                                     &
     &             GamT*cW*(temp-lambda2-lambda3*depth+sIce)
      c = GamS*salt*(cI*(lambda2+lambda3*depth-tIce) + L) +             &
     &     GamT*sIce*cW*(temp-lambda2-lambda3*depth)
!
      sB   = (1.0/(2.0*a))*(-b-SQRT(b**2 - 4.0*a*c))
      tB   = lambda1*sB+lambda2+lambda3*depth
      mdot = GamS*SQRT(CdBkg)*ABS(vel)*(salt-sB)/sB
!
      END SUBROUTINE iceplume_meltrate
!
!=======================================================================
!                                                                      !
! These are the core functions of iceplume entrainment model.          !
!                                                                      !
!=======================================================================
!
      SUBROUTINE iceplume_entrain(ng, I, iceDepthK,                     &
     &                            fIni, tIni, sIni,                     &
     &                            sgTyp, sgLen)
!
      USE mod_iceplume
!
!=======================================================================
!                                                                      !
! Input variables:                                                     !
!                                                                      !
!=======================================================================
!                                                                      !
! ng          - grid identifier                                        !
! I           - point source index identifier                          !
! iceDepthK   - vertical grid index of grounding line depth            !
! fIni        - subglacial runoff volume flux [m^3 s^-1]               !
! tIni        - subglacial runoff temperature [degC]                   !
! sIni        - subglacial runoff salinity [PSU]                       !
! sgTyp       - subglacial runoff type identifier                      !
! sgLen       - subglacial runoff discharge length [m]                 !
!                                                                      !
!=======================================================================
!                                                                      !
! Local variables:                                                     !
!                                                                      !
!=======================================================================
!                                                                      !
! rhoA/tA/sA  - ambient density/temperature/salinity                   !
! zIn/zOut    - in/out depth of the model steps                        !
!                                                                      !
!=======================================================================
!
! In/out variables
!
      integer, intent(in)  :: ng, I, iceDepthK
      real(r8), intent(in) :: fIni, tIni, sIni
      integer, intent(in)  :: sgTyp
      real(r8), intent(in) :: sgLen
!
! Local variables declaration
!
      real(r8) :: rhoA, tA, sA
      real(r8) :: zIn, zOut
      integer  :: K
!
! Local variables for ODEPACK
!
      integer :: IOPT, ISTATE, ITASK, ITOL, IWORK(20), LIW, LRW, MF, NEQ
      real(r8), parameter :: RTOL = 1.0E-5_r8
      real(r8), parameter :: ATOL = 1.0E-5_r8
      real(r8) :: RWORK(116), Y(10)
!
!=======================================================================
!                                                                      !
! Y is input/output vector for DLSODE                                  !
!   Y(1)  - plume volume flux                                          !
!   Y(2)  - plume velocity                                             !
!   Y(3)  - plume temperature                                          !
!   Y(4)  - plume salinity                                             !
!   Y(5)  - integrated plume/glacier contact area                      !
!   Y(6)  - area integrated melt                                       !
!   Y(7)  - fake variable to pass in ng                                !
!   Y(8)  - fake variable to pass in I                                 !
!   Y(9)  - fake variable to pass in sgTyp                             !
!   Y(10) - fake variable to pass in sgLen                             !
!                                                                      !
!=======================================================================
!                                                                      !
! For ODEPACK solver. See ODEPACK documentation and source code in     !
! Cowton et al. 2015.                                                  !
!                                                                      !
!=======================================================================
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
      PLUME(ng) % rho(I, iceDepthK)  = rho(Y(3), Y(4), zIn)
!
      CALL plume_metrics (sgTyp, Y(1)/Y(2), sgLen,                      &
     &                    PLUME(ng) % lm(I, iceDepthK),                 &
     &                    PLUME(ng) % lc(I, iceDepthK))
!
! Move up through water column from lowest layer
!
      DO K = iceDepthK+1, N(ng)
!
!=======================================================================
!                                                                      !
! Use DLSODE to solve plume properties.                                !
!                                                                      !
!=======================================================================
!
! Check to make sure plume has not reached neutral buoyancy in a lower
! layer
!
        IF (ISTATE .GT. -1) THEN
          CALL DLSODE (GENERAL_ENTRAIN_MODEL,                           &
     &                 NEQ, Y, zIn, zOut,                               &
     &                 ITOL, RTOL, ATOL, ITASK,                         &
     &                 ISTATE, IOPT, RWORK, LRW, IWORK,                 &
     &                 LIW, JEX, MF)
!
          CALL plume_metrics (sgTyp, Y(1)/Y(2), sgLen,                  &
     &                        PLUME(ng) % lm(I, K),                     &
     &                        PLUME(ng) % lc(I, K))
!
!=======================================================================
!                                                                      !
! Test to see if neutral buoyancy has now been reached. If solver      !
! returns ISTATE = -1, then it has been unable to meet required        !
! tolerances at this level. This generally occurs because plume        !
! has reached neutral buoyancy and run out of momentum, and so is      !
! no longer rising. At this point, we therefore end the call to the    !
! plume model. Our aim is to catch the plume at the point of           !
! neutral buoyancy. We therefore perform a manual comparrison of       !
! ambient and plume density. If plume density >= ambient density we    !
! assign ISTATE = -1, again ending the call to the plume model.        !
!                                                                      !
!=======================================================================
!
! Calculate plume density
!
          PLUME(ng) % rho(I, K) = rho(Y(3), Y(4), zIn)
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
          rhoA = rho(tA, sA, zIn)
          IF ((PLUME(ng)%rho(I, K) .GT. rhoA) .OR. (K .EQ. N(ng))) THEN
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
      END SUBROUTINE iceplume_entrain
!
!=======================================================================
!                                                                      !
! Below is the generalized plume model.                                !
!                                                                      !
!=======================================================================
!
      SUBROUTINE general_entrain_model (NEQ, T, Y, YDOT)
!
      USE mod_iceplume
!
      integer, intent(in) :: NEQ
      integer :: ng, I, sgTyp
      real(r8) :: sgLen
      real(r8) :: T, Y(10), YDOT(6)
!
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
        CALL linint(N(ng), PLUME(ng) % zR(I, :),                        &
     &              PLUME(ng) % tAm(I, :), T, tA)
        CALL linint(N(ng), PLUME(ng) % zR(I, :),                        &
     &              PLUME(ng) % sAm(I, :), T, sA)
!
! Calculate reduced gravity
!
        rhoP = rho(Y(3), Y(4), T)
        rhoA = rho(tA, sA, T)
        gRed  = g*(rhoA-rhoP)/rhoRef
!
! Calculate plume metrics
!
        CALL plume_metrics (sgTyp, Y(1)/Y(2), sgLen, Lm, Lc)
!
! Equations for sB, tB and mdot
!
        a = lambda1*(GamT*cW - GamS*cI)
        b = GamS*cI*(lambda1*Y(4) - lambda2 - lambda3*T +               &
     &               tIce - (L/cI)) -                                   &
     &               GamT*cW*(Y(3) - lambda2 - lambda3*T + sIce)
        c = GamS*Y(4)*(cI*(lambda2 + lambda3*T - tIce) + L) +           &
     &      GamT*sIce*cW*(Y(3) - lambda2 - lambda3*T)
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
        YDOT(2) = (1./Y(1))*(-Y(2)*YDOT(1) + gRed*Y(1)/Y(2) -           &
     &                       Cd*Lm*Y(2)**2)
!
! Plume temperature
!
        YDOT(3) = (1./Y(1))*(-Y(3)*YDOT(1) + alpha*Lc*Ta*Y(2) +         &
     &                       Lm*mdot*tB -                               &
     &                       SQRT(Cd)*GamT*Lm*Y(2)*(Y(3)-tB))
!
! Plume salinity
!
        YDOT(4) = (1./Y(1))*(-Y(4)*YDOT(1) + alpha*Lc*Sa*Y(2) +         &
     &                       Lm*mdot*sB -                               &
     &                       SQRT(Cd)*GamS*Lm*Y(2)*(Y(4)-sB))
!
! Along-plume integrated contact area and melt rate
!
        YDOT(5) = Lm
        YDOT(6) = Lm*mdot
      ENDIF
!
      END SUBROUTINE general_entrain_model
!
!=======================================================================
!                                                                      !
! Use this function to calculate plume metrics.                        !
!                                                                      !
!=======================================================================
!
      SUBROUTINE plume_metrics (sgTyp, area, Ls, Lm, Lc)
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
!
      integer, intent(in)   :: sgTyp
      real(r8), intent(in)  :: area, Ls
      real(r8), intent(out) :: Lm, Lc
!
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
      END SUBROUTINE plume_metrics
!
!=======================================================================
!                                                                      !
! Dummy routine for ODEPACK. Necessary for Jacobian matrix if          !
! stiff ODEs.                                                          !
!                                                                      !
!=======================================================================
!
      SUBROUTINE jex()
        RETURN
      END SUBROUTINE jex
!
!=======================================================================
!                                                                      !
! Use this function to calculate the velocity of detrainment.          !
!                                                                      !
!=======================================================================
!
      SUBROUTINE iceplume_detrain(ng, I,                                &
     &                            iceDepthK, plumeDepthK, osDepthK,     &
     &                            lc, det)
!
      USE mod_iceplume
!
!=======================================================================
!                                                                      !
! Input/Output variables:                                              !
!                                                                      !
!=======================================================================
!                                                                      !
! ng          - grid identifier                                        !
! I           - point source index identifier                          !
! iceDepthK   - glacier grounding line depth vertical grid index       !
! plumeDepthK - plume detrainment depth vertical grid index            !
! osDepthK    - plume overshoot depth vertical grid index              !
! lc          - plume/ocean contact length [m]                         !
! det         - detrainment volume flux [m^3 s^-1]                     !
!                                                                      !
!=======================================================================
!                                                                      !
! Local variables:                                                     !
!                                                                      !
!=======================================================================
!                                                                      !
! detDz       - detrainment layer total thickness [m]                  !
! minDz       - minimum detrainment layer thickness [m]                !
! rhoP        - plume density [kg m^-3]                                !
! potE        - potential energy of plume/ambient water                !
!                                                                      !
! detVel      - detrainment velocity [m s^-1]                          !
! detWeight   - detrainment weight function                            !
! detSum      - detrainment volume flux before normalize [m^3 s^-1]    !
!                                                                      !
! isSurface   - logical, if the plume discharge at surface             !
! isBottom    - logical, if the plume discharge at bottom              !
! udSwitch    - switch to determine search one layer up(1)/down(-1)    !
!                                                                      !
! maxVel      - maximum velocity based on Richardson number            !
!               criterion [m s^-1]                                     !
! Fr          - Froude number                                          !
! rho1        - upper layer density [kg m^-3]                          !
! rho2        - lower layer density [kg m^-3]                          !
! h1          - upper layer thickness [m]                              !
! h2          - lower layer thickness [m]                              !
! gRed        - reduced gravity [m s^-2]                               !
# ifdef ICEPLUME_DET_NEUTRAL
!                                                                      !
! volR        - detrainment volume of runoff [m^3]                     !
! volM        - detrainment volume of meltwater                        !
! volF        - detrainment volume of total freshwater                 !
! volE        - detrainment volume of entrainment                      !
! volP        - detrainment volume of plume                            !
!                                                                      !
! tF          - temperature of total freshwater [degC]                 !
! tE          - temperature of entrainment                             !
! tP          - temperature of plume                                   !
!                                                                      !
! sF          - salinity of total freshwater [PSU]                     !
! sE          - salinity of entrainment                                !
! sP          - salinity of plume                                      !
!                                                                      !
! rhoF        - density of total freshwater [kg m^-3]                  !
! rhoE        - density of entrainment                                 !
! rhoP        - density of plume                                       !
# endif
!                                                                      !
!=======================================================================
!
! In/out variables
!
      integer, intent(in)  :: ng, I
      integer, intent(in)  :: iceDepthK, plumeDepthK, osDepthK
      real(r8), intent(in) :: lc, det
!
! Local variables declaration
!
      real(r8) :: detDz, minDz, detVel, detWeight, potE
      real(r8) :: rhoP
      integer  :: udSwitch
      logical  :: isSurface = .FALSE.
      logical  :: isBottom = .FALSE.
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
# ifdef ICEPLUME_DET_SMOOTH
      real(r8) :: potEi, potEp, potEn, potERatio
# endif
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
!=======================================================================
!                                                                      !
! PART I - Find minimum detrainment depth                              !
!                                                                      !
!=======================================================================
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
      rhoP = rho(PLUME(ng) % t(I, osDepthK),                            &
     &           PLUME(ng) % s(I, osDepthK),                            &
     &           PLUME(ng) % zR(I, plumeDepthK))
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
# ifdef ICEPLUME_DET_SMOOTH
!
! Add a correction term to generate a smooth transition from subsurface
! plume to surface plume. The idea is to nudge the Richardson number
! towards 6.0 (thus the outflow behaves more like a surface plume) when
! the vertical momentum at surface is non-zero.
!
! This is an engineering fix in order to address Becca's concerns when
! running iceplume with a realistic topography.
!
! Chuning Wang 2021-05-17
!
        cff = 0.0
        DO K = 0,N(ng)
          IF (PLUME(ng) % w(I, K) .gt. cff) cff = PLUME(ng) % w(I, K)
        ENDDO
        cff1 = PLUME(ng) % w(I, N(ng)) / cff
        Ri = Ri + cff1*(6.0 - Ri)
# endif
!
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
            IF ((PLUME(ng) % detI(I, K  ) .EQ. 0) .AND.                 &
     &          (PLUME(ng) % detI(I, K-1) .EQ. 1)) THEN
              KI = K
            ENDIF
          ENDDO
          K1 = KI
        ELSEIF (udSwitch .EQ. -1) THEN
!
! Search one layer down
!
          DO K = 1, N(ng)-1
            IF ((PLUME(ng) % detI(I, K+1) .EQ. 1) .AND.                 &
     &          (PLUME(ng) % detI(I, K  ) .EQ. 0)) THEN
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
        potE = potE +                                                   &
     &         g*PLUME(ng) % dz(I, KI)*(PLUME(ng) % rhoAm(I, KI)-rhoP)
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
        cff = 1 - potE / (g*(PLUME(ng) % rhoAm(I, KI)-rhoP)) /          &
     &            PLUME(ng) % dz(I, KI)
        IF ( (cff .GE. 0) .AND. (cff .LT. 1) ) THEN
          PLUME(ng) % detFrac(I, KI) = cff
          detDz = detDz - PLUME(ng) % dz(I, KI)*(1 - cff)
        ENDIF
      ENDIF
# ifdef ICEPLUME_DET_SMOOTH
!
! Make correction for plume shape towards surface.
!
! Chuning Wang 2021-05-17
!
      potEp = 0.0
      potEn = 0.0
      DO K = K2,K1
        potEi = g*PLUME(ng) % detFrac(I, KI)*PLUME(ng) % dz(I, K)*      &
     &          (PLUME(ng) % rhoAm(I, K)-rhoP)
        IF (potEi .gt. 0.0) THEN
          potEp = potEp + potEi
        ELSE
          potEn = potEn - potEi
        ENDIF
      ENDDO
      IF ((potEp .gt. 0.0) .and. (potEn .gt. 0.0)) THEN
        potERatio = potEp/potEn
      ELSE
        potERatio = 1.0
      ENDIF
# endif
!
!=======================================================================
!                                                                      !
! PART II - Distribute detrainmnet in several layers                   !
!                                                                      !
! Update log                                                           !
! Use a Gause function to smooth the distribution.                     !
! 2018/06/08 Chuning Wang                                              !
!                                                                      !
!=======================================================================
!
! Update detrainment volume flux
!
      detVel = det / (lc * detDz)
      DO K = 1,N(ng)
        IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
!
! First, calculate weight using Gaussian function.
!
          h1 = PLUME(ng)%zR(I, K1)-PLUME(ng)%zR(I, plumeDepthK)+0.0001
          h2 = PLUME(ng)%zR(I, K2)-PLUME(ng)%zR(I, plumeDepthK)+0.0001
# ifdef ICEPLUME_DET_SMOOTH
          h1 = h1*potERatio
# endif
          IF ( K .GT. plumeDepthK ) THEN
            cff = (PLUME(ng)%zR(I, K)-PLUME(ng)%zR(I, plumeDepthK))/h1
          ELSE
            cff = (PLUME(ng)%zR(I, K)-PLUME(ng)%zR(I, plumeDepthK))/h2
          ENDIF
          detWeight = EXP(-0.5 * (cff/detSigma)**2)
!
! Distribute the detrainment.
!
          PLUME(ng) % det(I, K) = detWeight * detVel * lc *             &
     &      PLUME(ng) % detFrac(I, K) * PLUME(ng) % dz(I, K)
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
      rhoF = rho(tF, sF, PLUME(ng) % zR(I, plumeDepthK))
      rhoE = (rhoP*volP-rhoF*volF)/volE
!
! Calculate proportion from mixing line.
!
      DO K = 1, N(ng)
        IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
          cff = (PLUME(ng) % zR(I, K) -                                 &
     &           PLUME(ng) % zR(I, plumeDepthK))*0.005
          cff1 = MIN(MAX((PLUME(ng) % rhoAm(I, K)-rhoE)/(rhoF-rhoE),    &
     &                   0.01), 0.99)
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
          PLUME(ng) % det(I, K) = PLUME(ng) % detF(I, K) +              &
     &                            PLUME(ng) % detE(I, K)
        ENDIF
      ENDDO
!
!=======================================================================
!                                                                      !
! Update tracer concentration in detrainment model.                    !
!                                                                      !
!=======================================================================
!
! Active tracers.
!
      DO K = 1, N(ng)
        IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
          PLUME(ng) % detTrc(I, K, itemp) =                             &
     &        (tF*PLUME(ng) % detF(I, K)+tE*PLUME(ng) % detE(I, K)) /   &
     &        (PLUME(ng) % detF(I, K)+PLUME(ng) % detE(I, K))
          PLUME(ng) % detTrc(I, K, isalt) =                             &
     &        (sF*PLUME(ng) % detF(I, K)+sE*PLUME(ng) % detE(I, K)) /   &
     &        (PLUME(ng) % detF(I, K)+PLUME(ng) % detE(I, K))
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
            PLUME(ng) % detTrc(I, K, itrc) =                            &
     &        (cff1*PLUME(ng)%detF(I, K)+cff2*PLUME(ng)%detE(I, K)) /   &
     &        (PLUME(ng)%detF(I, K)+PLUME(ng)%detE(I, K))
          ENDIF
        ENDDO
      ENDDO
# endif
      END SUBROUTINE iceplume_detrain
!
!=======================================================================
!                                                                      !
! Use this function to calculate seawater density.                     !
!                                                                      !
!=======================================================================
!
      double precision FUNCTION rho(T,S,z)
!
! Equation of state (UNESCO 1983)
!     T = temperature (deg C)
!     S = salinity (PSU)
!     z = depth (m)
!
      double precision :: T,S,z
      double precision :: P
      double precision :: kw, Aw, Bw, k0
      double precision :: bulk_modulus
      double precision :: A, B, C, rho_w,rho_zero
!
      double precision, parameter :: rho_0 = 1027.0
      double precision, parameter :: g = 9.81
!

      P = rho_0*g*abs(z)*1.0e-5
!
! calculate in-situ rho
      kw = 19652.21 + 148.4206*T - 2.327105*T**2 + 1.360477e-2*(T**3) - &
     &     5.155288e-5*(T**4)
      Aw = 3.239908 + 1.43713e-3*T + 1.16092e-4*T**2 - 5.77905e-7*T**3
      Bw = 8.50935e-5 - 6.12293e-6*T + 5.2787e-8*(T**2)
      k0 = kw + (54.6746 - 0.603459*T + 1.09987e-2*(T**2) -             &
     &           6.1670e-5*(T**3))*S +                                  &
     &     (7.944e-2 + 1.6483e-2*T - 5.3009e-4*(T**2))*(S**1.5)
      A = Aw + (2.2838e-3 - 1.0981e-5*T - 1.6078e-6*(T**2))*S +         &
     &    1.91075e-4*(S**1.5)
      B = Bw + (-9.9348e-7 + 2.0816e-8*T + 9.1697e-10*T**2)*S
      bulk_modulus = k0 + A*P + B*P**2
!
      A = 8.24493e-1 - 4.0899e-3*T + 7.6438e-5*T**2 - 8.2467e-7*T**3 +  &
     &    5.3875e-9*T**4
      B = -5.72466e-3 + 1.0227e-4*T - 1.6546e-6*T**2
      C = 4.8314e-4
      rho_w = 999.842594 + 6.793952e-2*T - 9.095290e-3*T**2+            &
     &        1.001685e-4*T**3 - 1.120083e-6*T**4 + 6.536336e-9*T**5
      rho_zero = rho_w + A*S + B*(S**1.5) + C*(S**2)
!
      rho = rho_zero/(1.0 - (P/bulk_modulus))
!
      END FUNCTION rho
!
!=======================================================================
!                                                                      !
! These subroutines are modified from MITgcm.                          !
!                                                                      !
!=======================================================================
!
      SUBROUTINE sw_temp(S, T, P, PR, rv)
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
! DEEP-SEA RES., 1973, Vol20, 401-408.
!
      USE mod_kinds
!
! === Global variables ===
! 
! !INPUT/OUTPUT PARAMETERS:
! === Routine arguments ===
! S      :: salinity
! T      :: potential temperature
! P      :: pressure
! PR     :: reference pressure
! rv:: return value (in-situ temeparture in degree C)
!
      real(r8) :: S, T, P, PR
      real(r8) :: rv
!
! !LOCAL VARIABLES:
! === local variables ===
!
      CALL sw_ptmp  (S, T, PR, P, rv)
!
      END SUBROUTINE sw_temp
!
      SUBROUTINE sw_ptmp  (S, T, P, PR, rv)
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
!
! === Global variables ===
!
! !INPUT/OUTPUT PARAMETERS:
! === Routine arguments ===
! S  :: salinity    [psu      (PSS-78) ]
! T  :: temperature [degree C (IPTS-68)]
! P  :: pressure    [db]
! PR:: Reference pressure  [db]
! rv:: return value (potential temeparture in degree C)
!
      real(r8) :: S, T, P, PR
      real(r8) :: rv
!
! !LOCAL VARIABLES
! === local variables ===
!
      real(r8) :: del_P, del_th, th, q
      real(r8) :: onehalf, two, three
      parameter ( onehalf = 0.5, two = 2.0, three = 3.0 )
      real(r8) :: adtg_val
!
! theta1
!
      del_P   = PR-P
      CALL sw_adtg(S, T, P, adtg_val)
      del_th  = del_P*adtg_val
      th      = T+onehalf*del_th
      q       = del_th
!
! theta2
!
      CALL sw_adtg(S, th, P+onehalf*del_P, adtg_val)
      del_th  = del_P*adtg_val
      th      = th + (1-1/SQRT(two))*(del_th-q)
      q       = (two-sqrt(two))*del_th + (-two+three/SQRT(two))*q
!
! theta3
!
      CALL sw_adtg(S, th, P+onehalf*del_P, adtg_val)
      del_th  = del_P*adtg_val
      th      = th + (1+1/SQRT(two))*(del_th-q)
      q       = (two+SQRT(two))*del_th + (-two-three/SQRT(two))*q
!
! theta4
!
      CALL sw_adtg(S, th, P+del_P, adtg_val)
      del_th  = del_P*adtg_val
      rv      = th + (del_th-two*q)/(two*three)
!
      END SUBROUTINE sw_ptmp
!
      SUBROUTINE sw_adtg  (S, T, P, rv)
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
!
! === Global variables ===
!
! !INPUT/OUTPUT PARAMETERS:
! === Routine arguments ===
!
      real(r8) :: S, T, P
      real(r8) :: rv
!
! !LOCAL VARIABLES:
! === local variables ===
!
      real(r8) :: a0, a1, a2, a3, b0, b1, c0, c1, c2, c3
      real(r8) :: d0, d1, e0, e1, e2
      real(r8) :: sref
!
      sref = 35.0_r8
      a0 =  3.5803e-05_r8
      a1 =  8.5258e-06_r8
      a2 =  -6.836e-08_r8
      a3 =  6.6228e-10_r8
!
      b0 =  1.8932e-06_r8
      b1 = -4.2393e-08_r8
!
      c0 =  1.8741e-08_r8
      c1 = -6.7795e-10_r8
      c2 =   8.733e-12_r8
      c3 = -5.4481e-14_r8
!
      d0 = -1.1351e-10_r8
      d1 =  2.7759e-12_r8
!
      e0 = -4.6206e-13_r8
      e1 =  1.8676e-14_r8
      e2 = -2.1687e-16_r8
!
      rv =  a0 + (a1 + (a2+a3*T)*T)*T +                                 &
     &      (b0+b1*T)*(S-sref) +                                        &
     &      ( (c0 + (c1 + (c2+c3*T)*T)*T) + (d0+d1*T)*(S-sref) )*P +    &
     &      (  e0 + (e1+e2*T)*T )*P*P
!
      END SUBROUTINE sw_adtg
!
!=======================================================================
!                                                                      !
! Use this function to do 1D interpolation.                            !
!                                                                      !
!=======================================================================
!
      SUBROUTINE linint(nx, xtab, ytab, x, y)
!
! Given a value of x return a value of y based on interpolation
! within a table of y values (ytab) corresponding to the x values
! contained in the array xtab. The subroutine assumes that the
! values in xtab increase monotonically
!
! John Mahaffy 2/12/95
! Modified by CW 2020
!
      integer, intent(in)           :: nx
      double precision, intent(in)  :: xtab(nx), ytab(nx)
      double precision, intent(in)  :: x
      double precision, intent(out) :: y
!
! local variables
!
      integer :: i
      double precision :: wx
!
      IF (x .LE. xtab(1)) THEN
        y = ytab(1)
      ELSEIF (x .GE. xtab(nx)) THEN
        y = ytab(nx)
      ELSE
        DO i = 2, nx
          IF (x .LE. xtab(i)) THEN
            wx = (x-xtab(i-1))/(xtab(i)-xtab(i-1))
            y  = (1-wx)*ytab(i-1)+wx*ytab(i)
            EXIT
          ENDIF
        ENDDO
      ENDIF
!
      END SUBROUTINE linint
#endif
      END MODULE iceplume_mod
