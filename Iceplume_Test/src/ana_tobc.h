      SUBROUTINE ana_tobc (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets tracer-type variables open boundary conditions    !
!  using analytical expressions.                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_tobc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp(ng),                                     &
     &                    GRID(ng) % z_r,                               &
     &                    OCEAN(ng) % t)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(34)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_tobc
!
!***********************************************************************
      SUBROUTINE ana_tobc_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp,                                   &
     &                          z_r, t)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_boundary
      USE mod_ncparam
      USE mod_ocean
#ifdef SEDIMENT
      USE mod_sediment
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#else
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: i, ised, itrc, j, k
      real(r8) :: cff

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Tracers open boundary conditions.
!-----------------------------------------------------------------------
!
#ifdef ICEPLUME_TEST
      IF (ANY(LBC(ieast,isTvar(:),ng)%acquire).and.                     &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO j=JstrT,JendT
            BOUNDARY(ng)%t_east(j,k,itemp)=T0(ng)
# ifdef SALINITY
            BOUNDARY(ng)%t_east(j,k,isalt)=S0(ng)-(S0(ng)-USER(2))*     &
     &        0.5_r8*(TANH(0.2_r8*pi*(z_r(Iend+1,j,k)+USER(1)))+1.0_r8)
# endif
# ifdef SEDIMENT
            DO ised=1,NST
              BOUNDARY(ng)%t_east(j,k,idsed(ised))=0.0_r8
            END DO
# endif
          END DO
        END DO
      END IF

      IF (ANY(LBC(iwest,isTvar(:),ng)%acquire).and.                     &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO k=1,N(ng)
          DO j=JstrT,JendT
            BOUNDARY(ng)%t_west(j,k,itemp)=T0(ng)
# ifdef SALINITY
            BOUNDARY(ng)%t_west(j,k,isalt)=S0(ng)-(S0(ng)-USER(2))*     &
     &        0.5_r8*(TANH(0.2_r8*pi*(z_r(0,j,k)+USER(1)))+1.0_r8)
# endif
# ifdef SEDIMENT
            DO ised=1,NST
              BOUNDARY(ng)%t_west(j,k,idsed(ised))=0.0_r8
            END DO
# endif
          END DO
        END DO
      END IF

      IF (ANY(LBC(isouth,isTvar(:),ng)%acquire).and.                    &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO i=IstrT,IendT
            BOUNDARY(ng)%t_south(i,k,itemp)=T0(ng)
# ifdef SALINITY
            BOUNDARY(ng)%t_south(i,k,isalt)=S0(ng)-(S0(ng)-USER(2))*    &
     &        0.5_r8*(TANH(0.2_r8*pi*(z_r(i,0,k)+USER(1)))+1.0_r8)
# endif
# ifdef SEDIMENT
            DO ised=1,NST
              BOUNDARY(ng)%t_south(i,k,idsed(ised))=0.0_r8
            END DO
# endif
          END DO
        END DO
      END IF

      IF (ANY(LBC(inorth,isTvar(:),ng)%acquire).and.                    &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO i=IstrT,IendT
            BOUNDARY(ng)%t_north(i,k,itemp)=T0(ng)
# ifdef SALINITY
            BOUNDARY(ng)%t_north(i,k,isalt)=S0(ng)-(S0(ng)-USER(2))*    &
     &        0.5_r8*(TANH(0.2_r8*pi*(z_r(i,Jend+1,k)+USER(1)))+1.0_r8)
# endif
# ifdef SEDIMENT
            DO ised=1,NST
              BOUNDARY(ng)%t_north(i,k,idsed(ised))=0.0_r8
            END DO
# endif
          END DO
        END DO
      END IF

#else

      IF (ANY(LBC(ieast,isTvar(:),ng)%acquire).and.                     &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrT,JendT
              BOUNDARY(ng)%t_east(j,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END IF

      IF (ANY(LBC(iwest,isTvar(:),ng)%acquire).and.                     &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrT,JendT
              BOUNDARY(ng)%t_west(j,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END IF

      IF (ANY(LBC(isouth,isTvar(:),ng)%acquire).and.                    &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=IstrT,IendT
              BOUNDARY(ng)%t_south(i,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END IF

      IF (ANY(LBC(inorth,isTvar(:),ng)%acquire).and.                    &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=IstrT,IendT
              BOUNDARY(ng)%t_north(i,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END IF

#endif

      RETURN
      END SUBROUTINE ana_tobc_tile
