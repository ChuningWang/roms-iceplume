#include "cppdefs.h"
      MODULE uv3dmix2_geo_mod
#if defined SOLVE3D && defined UV_VIS2 && defined MIX_GEO_UV
# ifdef EW_PERIODIC
#  define IV_RANGE Istr-1,Iend+1
#  define IU_RANGE Istr-1,Iend+1
# else
#  define IV_RANGE MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
#  define IU_RANGE MAX(2,IstrU-1),MIN(Iend+1,Lm(ng))
# endif
# ifdef NS_PERIODIC
#  define JU_RANGE Jstr-1,Jend+1
#  define JV_RANGE Jstr-1,Jend+1
# else
#  define JU_RANGE MAX(1,Jstr-1),MIN(Jend+1,Mm(ng))
#  define JV_RANGE MAX(2,JstrV-1),MIN(Jend+1,Mm(ng))
# endif
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes harmonic mixing of momentum, rotated along    !
!  geopotentials,  from the  horizontal  divergence  of the  stress    !
!  tensor.  A transverse  isotropy is assumed so the stress tensor     !
!  is split into vertical and horizontal subtensors.                   !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!      Wajsowicz, R.C, 1993: A consistent formulation of the           !
!         anisotropic stress tensor for use in models of the           !
!         large-scale ocean circulation, JCP, 105, 333-338.            !
!                                                                      !
!      Sadourny, R. and K. Maynard, 1997: Formulations of              !
!         lateral diffusion in geophysical fluid dynamics              !
!         models, In Numerical Methods of Atmospheric and              !
!         Oceanic Modelling. Lin, Laprise, and Ritchie,                !
!         Eds., NRC Research Press, 547-556.                           !
!                                                                      !
!      Griffies, S.M. and R.W. Hallberg, 2000: Biharmonic              !
!         friction with a Smagorinsky-like viscosity for               !
!         use in large-scale eddy-permitting ocean models,             !
!         Monthly Weather Rev., 128, 8, 2935-2946.                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: uv3dmix2_geo

      CONTAINS
!
!***********************************************************************
      SUBROUTINE uv3dmix2_geo (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
# ifdef DIAGNOSTICS_UV
      USE mod_diags
# endif
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, 31)
# endif
      CALL uv3dmix2_geo_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        nrhs(ng), nnew(ng),                       &
# ifdef MASKING
     &                        GRID(ng) % pmask,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % vmask,                         &
# endif
     &                        GRID(ng) % Hz,                            &
     &                        GRID(ng) % om_p,                          &
     &                        GRID(ng) % om_r,                          &
     &                        GRID(ng) % om_u,                          &
     &                        GRID(ng) % om_v,                          &
     &                        GRID(ng) % on_p,                          &
     &                        GRID(ng) % on_r,                          &
     &                        GRID(ng) % on_u,                          &
     &                        GRID(ng) % on_v,                          &
     &                        GRID(ng) % pm,                            &
     &                        GRID(ng) % pn,                            &
     &                        GRID(ng) % z_r,                           &
     &                        MIXING(ng) % visc2_p,                     &
     &                        MIXING(ng) % visc2_r,                     &
# ifdef DIAGNOSTICS_UV
     &                        DIAGS(ng) % DiaRUfrc,                     &
     &                        DIAGS(ng) % DiaRVfrc,                     &
     &                        DIAGS(ng) % DiaU3wrk,                     &
     &                        DIAGS(ng) % DiaV3wrk,                     &
# endif
     &                        COUPLING(ng) % rufrc,                     &
     &                        COUPLING(ng) % rvfrc,                     &
     &                        OCEAN(ng) % u,                            &
     &                        OCEAN(ng) % v)
# ifdef PROFILE
      CALL wclock_off (ng, 31)
# endif
      RETURN
      END SUBROUTINE uv3dmix2_geo
!
!***********************************************************************
      SUBROUTINE uv3dmix2_geo_tile (ng, Istr, Iend, Jstr, Jend,         &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              nrhs, nnew,                         &
# ifdef MASKING
     &                              pmask, rmask, umask, vmask,         &
# endif
     &                              Hz,                                 &
     &                              om_p, om_r, om_u, om_v,             &
     &                              on_p, on_r, on_u, on_v,             &
     &                              pm, pn, z_r,                        &
     &                              visc2_p, visc2_r,                   &
# ifdef DIAGNOSTICS_UV
     &                              DiaRUfrc, DiaRVfrc,                 &
     &                              DiaU3wrk, DiaV3wrk,                 &
# endif
     &                              rufrc, rvfrc, u, v)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs, nnew

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_p(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_p(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: visc2_p(LBi:,LBj:)
      real(r8), intent(in) :: visc2_r(LBi:,LBj:)

#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRUfrc(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: DiaRVfrc(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: DiaU3wrk(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: DiaV3wrk(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: rufrc(LBi:,LBj:)
      real(r8), intent(inout) :: rvfrc(LBi:,LBj:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: visc2_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: visc2_r(LBi:UBi,LBj:UBj)

#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRUfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
      real(r8), intent(inout) :: DiaRVfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
      real(r8), intent(inout) :: DiaU3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
      real(r8), intent(inout) :: DiaV3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
#  endif
      real(r8), intent(inout) :: rufrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rvfrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k, k1, k2

      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff6, cff7, cff8
      real(r8) :: dmUdz, dnUdz, dmVdz, dnVdz, fac

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: UFe
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: VFe
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: UFx
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: VFx

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: UFs
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: VFs
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: dmUde
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: dmVde
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: dnUdx
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: dnVdx
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: dUdz
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: dVdz
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: dZde_p
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: dZde_r
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: dZdx_p
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: dZdx_r

# include "set_bounds.h"
!
!--------------------------------------------------------------------
!  Compute horizontal harmonic viscosity along geopotential surfaces.
!--------------------------------------------------------------------
!
!  Compute horizontal and vertical gradients.  Notice the recursive
!  blocking sequence.  The vertical placement of the gradients is:
!
!    dZdx_r, dZde_r, dnUdx, dmVde(:,:,k1) k      rho-points
!    dZdx_r, dZde_r, dnUdx, dmVde(:,:,k2) k+1    rho-points
!    dZdx_p, dZde_p, dnVdx, dmUde(:,:,k1) k      psi-points
!    dZdx_p, dZde_p, dnVdx, dmUde(:,:,k2) k+1    psi-points
!                       UFs, dUdz(:,:,k1) k-1/2  WU-points
!                       UFs, dUdz(:,:,k2) k+1/2  WU-points
!                       VFs, dVdz(:,:,k1) k-1/2  WV-points
!                       VFs, dVdz(:,:,k2) k+1/2  WV-points
!
      k2=1
      K_LOOP : DO k=0,N(ng)
        k1=k2
        k2=3-k1
        IF (k.lt.N(ng)) THEN
!
!  Compute slopes (nondimensional) at RHO- and PSI-points.
!
          DO j=-1+JU_RANGE+1
            DO i=IV_RANGE+1
              UFx(i,j)=0.5_r8*(pm(i-1,j)+pm(i,j))*                      &
     &                 (z_r(i,j,k+1)-z_r(i-1,j,k+1))
# ifdef MASKING
              UFx(i,j)=UFx(i,j)*umask(i,j)
# endif
            END DO
            DO i=IV_RANGE
              dZdx_r(i,j,k2)=0.5_r8*(UFx(i,j)+UFx(i+1,j))
            END DO
          END DO
          DO j=JU_RANGE+1
            DO i=-1+IV_RANGE+1
              VFe(i,j)=0.5_r8*(pn(i,j-1)+pn(i,j))*                      &
     &                 (z_r(i,j,k+1)-z_r(i,j-1,k+1))
# ifdef MASKING
              VFe(i,j)=VFe(i,j)*vmask(i,j)
# endif
            END DO
            DO i=IV_RANGE+1
              dZde_p(i,j,k2)=0.5_r8*(VFe(i-1,j)+VFe(i,j))
              dZdx_p(i,j,k2)=0.5_r8*(UFx(i,j-1)+UFx(i,j))
            END DO
          END DO
          DO j=JU_RANGE
            DO i=-1+IV_RANGE+1
              dZde_r(i,j,k2)=0.5_r8*(VFe(i,j)+VFe(i,j+1))
            END DO
          END DO
!
!  Compute momentum horizontal (1/m/s) and vertical (1/s) gradients.
!
          DO j=-1+JV_RANGE
            DO i=-1+IU_RANGE
              dnUdx(i,j,k2)=0.5_r8*pm(i,j)*                             &
     &                      ((pn(i  ,j)+pn(i+1,j))*u(i+1,j,k+1,nrhs)-   &
     &                       (pn(i-1,j)+pn(i  ,j))*u(i  ,j,k+1,nrhs))
# ifdef MASKING
              dnUdx(i,j,k2)=dnUdx(i,j,k2)*rmask(i,j)
# endif
              dmVde(i,j,k2)=0.5_r8*pn(i,j)*                             &
     &                      ((pm(i,j  )+pm(i,j+1))*v(i,j+1,k+1,nrhs)-   &
     &                       (pm(i,j-1)+pm(i,j  ))*v(i,j  ,k+1,nrhs))
# ifdef MASKING
              dmVde(i,j,k2)=dmVde(i,j,k2)*rmask(i,j)
# endif
            END DO
          END DO
          DO j=JU_RANGE+1
            DO i=IV_RANGE+1
              dmUde(i,j,k2)=0.125_r8*(pn(i-1,j  )+pn(i,j  )+            &
     &                                pn(i-1,j-1)+pn(i,j-1))*           &
     &                      ((pm(i-1,j  )+pm(i,j  ))*u(i,j  ,k+1,nrhs)- &
     &                       (pm(i-1,j-1)+pm(i,j-1))*u(i,j-1,k+1,nrhs))
# ifdef MASKING
              dmUde(i,j,k2)=dmUde(i,j,k2)*pmask(i,j)
# endif
              dnVdx(i,j,k2)=0.125_r8*(pm(i-1,j  )+pm(i,j  )+            &
     &                                pm(i-1,j-1)+pm(i,j-1))*           &
     &                      ((pn(i  ,j-1)+pn(i  ,j))*v(i  ,j,k+1,nrhs)- &
     &                       (pn(i-1,j-1)+pn(i-1,j))*v(i-1,j,k+1,nrhs))
# ifdef MASKING
              dnVdx(i,j,k2)=dnVdx(i,j,k2)*pmask(i,j)
# endif
            END DO
          END DO
        END IF
        IF ((k.eq.0).or.(k.eq.N(ng))) THEN
          DO j=JU_RANGE
            DO i=-1+IU_RANGE+1
              dUdz(i,j,k2)=0.0_r8
              UFs(i,j,k2)=0.0_r8
            END DO
          END DO
          DO j=-1+JV_RANGE+1
            DO i=IV_RANGE
              dVdz(i,j,k2)=0.0_r8
              VFs(i,j,k2)=0.0_r8
            END DO
          END DO
        ELSE
          DO j=JU_RANGE
            DO i=-1+IU_RANGE+1
              dUdz(i,j,k2)=(u(i,j,k+1,nrhs)-u(i,j,k,nrhs))/             &
     &                     (0.5_r8*(z_r(i-1,j,k+1)-z_r(i-1,j,k)+        &
     &                              z_r(i  ,j,k+1)-z_r(i  ,j,k)))
            END DO
          END DO
          DO j=-1+JV_RANGE+1
            DO i=IV_RANGE
              dVdz(i,j,k2)=(v(i,j,k+1,nrhs)-v(i,j,k,nrhs))/             &
     &                     (0.5_r8*(z_r(i,j-1,k+1)-z_r(i,j-1,k)+        &
     &                              z_r(i,j  ,k+1)-z_r(i,j  ,k)))
            END DO
          END DO
        END IF
!
!  Compute components of the rotated viscous flux (m5/s2) along
!  geopotential surfaces in the XI- and ETA-directions.
!
        IF (k.gt.0) THEN
          DO j=-1+JV_RANGE
            DO i=-1+IU_RANGE
              cff=visc2_r(i,j)*Hz(i,j,k)*                               &
     &            (on_r(i,j)*(dnUdx(i,j,k1)-0.5_r8*pn(i,j)*             &
     &                        (MIN(dZdx_r(i,j,k1),0.0_r8)*              &
     &                             (dUdz(i,j,k1)+dUdz(i+1,j,k2))+       &
     &                         MAX(dZdx_r(i,j,k1),0.0_r8)*              &
     &                             (dUdz(i,j,k2)+dUdz(i+1,j,k1))))-     &
     &             om_r(i,j)*(dmVde(i,j,k1)-0.5_r8*pm(i,j)*             &
     &                        (MIN(dZde_r(i,j,k1),0.0_r8)*              &
     &                             (dVdz(i,j,k1)+dVdz(i,j+1,k2))+       &
     &                         MAX(dZde_r(i,j,k1),0.0_r8)*              &
     &                             (dVdz(i,j,k2)+dVdz(i,j+1,k1)))))
# ifdef MASKING
              cff=cff*rmask(i,j)
# endif
              UFx(i,j)=on_r(i,j)*on_r(i,j)*cff
              VFe(i,j)=om_r(i,j)*om_r(i,j)*cff
            END DO
          END DO
          DO j=JU_RANGE+1
            DO i=IV_RANGE+1
              cff=visc2_p(i,j)*0.25_r8*(Hz(i-1,j  ,k)+Hz(i,j  ,k)+      &
     &                                  Hz(i-1,j-1,k)+Hz(i,j-1,k))*     &
     &            (on_p(i,j)*(dnVdx(i,j,k1)-                            &
     &                        0.125_r8*(pn(i-1,j-1)+pn(i-1,j)+          &
     &                                  pn(i  ,j-1)+pn(i  ,j))*         &
     &                        (MIN(dZdx_p(i,j,k1),0.0_r8)*              &
     &                             (dVdz(i-1,j,k1)+dVdz(i,j,k2))+       &
     &                         MAX(dZdx_p(i,j,k1),0.0_r8)*              &
     &                             (dVdz(i-1,j,k2)+dVdz(i,j,k1))))+     &
     &             om_p(i,j)*(dmUde(i,j,k1)-                            &
     &                        0.125_r8*(pm(i-1,j-1)+pm(i-1,j)+          &
     &                                  pm(i  ,j-1)+pm(i  ,j))*         &
     &                        (MIN(dZde_p(i,j,k1),0.0_r8)*              &
     &                             (dUdz(i,j-1,k1)+dUdz(i,j,k2))+       &
     &                         MAX(dZde_p(i,j,k1),0.0_r8)*              &
     &                             (dUdz(i,j-1,k2)+dUdz(i,j,k1)))))
# ifdef MASKING
              cff=cff*pmask(i,j)
# endif
              UFe(i,j)=om_p(i,j)*om_p(i,j)*cff
              VFx(i,j)=on_p(i,j)*on_p(i,j)*cff
            END DO
          END DO
!
!  Compute vertical flux (m2/s2) due to sloping terrain-following
!  surfaces.
!
          IF (k.lt.N(ng)) THEN
            DO j=Jstr,Jend
              DO i=IstrU,Iend
                cff=0.5_r8*(pn(i-1,j)+pn(i,j))
                dnUdz=cff*dUdz(i,j,k2)
                dnVdz=cff*0.25_r8*(dVdz(i-1,j+1,k2)+dVdz(i,j+1,k2)+     &
     &                             dVdz(i-1,j  ,k2)+dVdz(i,j  ,k2))
                cff=0.5_r8*(pm(i-1,j)+pm(i,j))
                dmUdz=cff*dUdz(i,j,k2)
                dmVdz=cff*0.25_r8*(dVdz(i-1,j+1,k2)+dVdz(i,j+1,k2)+     &
     &                             dVdz(i-1,j  ,k2)+dVdz(i,j  ,k2))
                cff1=MIN(dZdx_r(i-1,j,k1),0.0_r8)
                cff2=MIN(dZdx_r(i  ,j,k2),0.0_r8)
                cff3=MAX(dZdx_r(i-1,j,k2),0.0_r8)
                cff4=MAX(dZdx_r(i  ,j,k1),0.0_r8)
                UFs(i,j,k2)=on_u(i,j)*                                  &
     &                      0.25_r8*(visc2_r(i-1,j)+visc2_r(i,j))*      &
     &                      (cff1*(cff1*dnUdz-dnUdx(i-1,j,k1))+         &
     &                       cff2*(cff2*dnUdz-dnUdx(i  ,j,k2))+         &
     &                       cff3*(cff3*dnUdz-dnUdx(i-1,j,k2))+         &
     &                       cff4*(cff4*dnUdz-dnUdx(i  ,j,k1)))
                cff1=MIN(dZde_p(i,j  ,k1),0.0_r8)
                cff2=MIN(dZde_p(i,j+1,k2),0.0_r8)
                cff3=MAX(dZde_p(i,j  ,k2),0.0_r8)
                cff4=MAX(dZde_p(i,j+1,k1),0.0_r8)
                UFs(i,j,k2)=UFs(i,j,k2)+om_u(i,j)*                      &
     &                      0.25_r8*(visc2_r(i-1,j)+visc2_r(i,j))*      &
     &                      (cff1*(cff1*dmUdz-dmUde(i,j  ,k1))+         &
     &                       cff2*(cff2*dmUdz-dmUde(i,j+1,k2))+         &
     &                       cff3*(cff3*dmUdz-dmUde(i,j  ,k2))+         &
     &                       cff4*(cff4*dmUdz-dmUde(i,j+1,k1)))
                cff1=MIN(dZde_p(i,j  ,k1),0.0_r8)
                cff2=MIN(dZde_p(i,j+1,k2),0.0_r8)
                cff3=MAX(dZde_p(i,j  ,k2),0.0_r8)
                cff4=MAX(dZde_p(i,j+1,k1),0.0_r8)
                cff5=MIN(dZdx_p(i,j  ,k1),0.0_r8)
                cff6=MIN(dZdx_p(i,j+1,k2),0.0_r8)
                cff7=MAX(dZdx_p(i,j  ,k2),0.0_r8)
                cff8=MAX(dZdx_p(i,j+1,k1),0.0_r8)
                UFs(i,j,k2)=UFs(i,j,k2)+on_u(i,j)*                      &
     &                      0.25_r8*(visc2_r(i-1,j)+visc2_r(i,j))*      &
     &                      (cff1*(cff5*dnVdz-dnVdx(i,j  ,k1))+         &
     &                       cff2*(cff6*dnVdz-dnVdx(i,j+1,k2))+         &
     &                       cff3*(cff7*dnVdz-dnVdx(i,j  ,k2))+         &
     &                       cff4*(cff8*dnVdz-dnVdx(i,j+1,k1)))
                cff1=MIN(dZdx_r(i-1,j,k1),0.0_r8)
                cff2=MIN(dZdx_r(i  ,j,k2),0.0_r8)
                cff3=MAX(dZdx_r(i-1,j,k2),0.0_r8)
                cff4=MAX(dZdx_r(i  ,j,k1),0.0_r8)
                cff5=MIN(dZde_r(i-1,j,k1),0.0_r8)
                cff6=MIN(dZde_r(i  ,j,k2),0.0_r8)
                cff7=MAX(dZde_r(i-1,j,k2),0.0_r8)
                cff8=MAX(dZde_r(i  ,j,k1),0.0_r8)
                UFs(i,j,k2)=UFs(i,j,k2)-om_u(i,j)*                      &
     &                      0.25_r8*(visc2_r(i-1,j)+visc2_r(i,j))*      &
     &                      (cff1*(cff5*dmVdz-dmVde(i-1,j,k1))+         &
     &                       cff2*(cff6*dmVdz-dmVde(i  ,j,k2))+         &
     &                       cff3*(cff7*dmVdz-dmVde(i-1,j,k2))+         &
     &                       cff4*(cff8*dmVdz-dmVde(i  ,j,k1)))
              END DO
            END DO
!
            DO j=JstrV,Jend
              DO i=Istr,Iend
                cff=0.5_r8*(pn(i,j-1)+pn(i,j))
                dnUdz=cff*0.25_r8*(dUdz(i,j  ,k2)+dUdz(i+1,j  ,k2)+     &
     &                             dUdz(i,j-1,k2)+dUdz(i+1,j-1,k2))
                dnVdz=cff*dVdz(i,j,k2)
                cff=0.5_r8*(pm(i,j-1)+pm(i,j))
                dmUdz=cff*0.25_r8*(dUdz(i,j  ,k2)+dUdz(i+1,j  ,k2)+     &
     &                             dUdz(i,j-1,k2)+dUdz(i+1,j-1,k2))
                dmVdz=cff*dVdz(i,j,k2)
                cff1=MIN(dZdx_p(i  ,j,k1),0.0_r8)
                cff2=MIN(dZdx_p(i+1,j,k2),0.0_r8)
                cff3=MAX(dZdx_p(i  ,j,k2),0.0_r8)
                cff4=MAX(dZdx_p(i+1,j,k1),0.0_r8)
                VFs(i,j,k2)=on_v(i,j)*                                  &
     &                      0.25_r8*(visc2_r(i,j-1)+visc2_r(i,j))*      &
     &                      (cff1*(cff1*dnVdz-dnVdx(i  ,j,k1))+         &
     &                       cff2*(cff2*dnVdz-dnVdx(i+1,j,k2))+         &
     &                       cff3*(cff3*dnVdz-dnVdx(i  ,j,k2))+         &
     &                       cff4*(cff4*dnVdz-dnVdx(i+1,j,k1)))
                cff1=MIN(dZde_r(i,j-1,k1),0.0_r8)
                cff2=MIN(dZde_r(i,j  ,k2),0.0_r8)
                cff3=MAX(dZde_r(i,j-1,k2),0.0_r8)
                cff4=MAX(dZde_r(i,j  ,k1),0.0_r8)
                VFs(i,j,k2)=VFs(i,j,k2)+om_v(i,j)*                      &
     &                      0.25_r8*(visc2_r(i,j-1)+visc2_r(i,j))*      &
     &                      (cff1*(cff1*dmVdz-dmVde(i,j-1,k1))+         &
     &                       cff2*(cff2*dmVdz-dmVde(i,j  ,k2))+         &
     &                       cff3*(cff3*dmVdz-dmVde(i,j-1,k2))+         &
     &                       cff4*(cff4*dmVdz-dmVde(i,j  ,k1)))
                cff1=MIN(dZde_r(i,j-1,k1),0.0_r8)
                cff2=MIN(dZde_r(i,j  ,k2),0.0_r8)
                cff3=MAX(dZde_r(i,j-1,k2),0.0_r8)
                cff4=MAX(dZde_r(i,j  ,k1),0.0_r8)
                cff5=MIN(dZdx_r(i,j-1,k1),0.0_r8)
                cff6=MIN(dZdx_r(i,j  ,k2),0.0_r8)
                cff7=MAX(dZdx_r(i,j-1,k2),0.0_r8)
                cff8=MAX(dZdx_r(i,j  ,k1),0.0_r8)
                VFs(i,j,k2)=VFs(i,j,k2)-on_v(i,j)*                      &
     &                      0.25_r8*(visc2_r(i,j-1)+visc2_r(i,j))*      &
     &                      (cff1*(cff5*dnUdz-dnUdx(i,j-1,k1))+         &
     &                       cff2*(cff6*dnUdz-dnUdx(i,j  ,k2))+         &
     &                       cff3*(cff7*dnUdz-dnUdx(i,j-1,k2))+         &
     &                       cff4*(cff8*dnUdz-dnUdx(i,j  ,k1)))
                cff1=MIN(dZdx_p(i  ,j,k1),0.0_r8)
                cff2=MIN(dZdx_p(i+1,j,k2),0.0_r8)
                cff3=MAX(dZdx_p(i  ,j,k2),0.0_r8)
                cff4=MAX(dZdx_p(i+1,j,k1),0.0_r8)
                cff5=MIN(dZde_p(i  ,j,k1),0.0_r8)
                cff6=MIN(dZde_p(i+1,j,k2),0.0_r8)
                cff7=MAX(dZde_p(i  ,j,k2),0.0_r8)
                cff8=MAX(dZde_p(i+1,j,k1),0.0_r8)
                VFs(i,j,k2)=VFs(i,j,k2)+om_v(i,j)*                      &
     &                      0.25_r8*(visc2_r(i,j-1)+visc2_r(i,j))*      &
     &                      (cff1*(cff5*dmUdz-dmUde(i  ,j,k1))+         &
     &                       cff2*(cff6*dmUdz-dmUde(i+1,j,k2))+         &
     &                       cff3*(cff7*dmUdz-dmUde(i  ,j,k2))+         &
     &                       cff4*(cff8*dmUdz-dmUde(i+1,j,k1)))
              END DO
            END DO
          END IF
!
! Time-step harmonic, geopotential viscosity term.  Notice that
! momentum at this stage is HzU and HzV and has m2/s units.  Add
! contribution for barotropic forcing terms.
!
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              cff1=0.5_r8*((pn(i-1,j)+pn(i,j))*(UFx(i,j  )-UFx(i-1,j))+ &
     &                     (pm(i-1,j)+pm(i,j))*(UFe(i,j+1)-UFe(i  ,j)))
              cff2=0.25_r8*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
              cff3=UFs(i,j,k2)-UFs(i,j,k1)
              fac=dt(ng)*(cff1*cff2+cff3)
              rufrc(i,j)=rufrc(i,j)+cff1+cff3
              u(i,j,k,nnew)=u(i,j,k,nnew)+fac
# ifdef DIAGNOSTICS_UV
              DiaRUfrc(i,j,3,M2hvis)=DiaRUfrc(i,j,3,M2hvis)+cff1+cff3
              DiaU3wrk(i,j,k,M3hvis)=fac
# endif
            END DO
          END DO
          DO j=JstrV,Jend
            DO i=Istr,Iend
              cff1=0.5_r8*((pn(i,j-1)+pn(i,j))*(VFx(i+1,j)-VFx(i,j  ))- &
     &                     (pm(i,j-1)+pm(i,j))*(VFe(i  ,j)-VFe(i,j-1)))
              cff2=0.25_8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
              cff3=VFs(i,j,k2)-VFs(i,j,k1)
              fac=dt(ng)*(cff1*cff2+cff3)
              rvfrc(i,j)=rvfrc(i,j)+cff1+cff3
              v(i,j,k,nnew)=v(i,j,k,nnew)+fac
# ifdef DIAGNOSTICS_UV
              DiaRVfrc(i,j,3,M2hvis)=DiaRVfrc(i,j,3,M2hvis)+cff1+cff3
              DiaV3wrk(i,j,k,M3hvis)=fac
# endif
            END DO
          END DO
        END IF
      END DO K_LOOP
# undef IU_RANGE
# undef IV_RANGE
# undef JU_RANGE
# undef JV_RANGE
      RETURN
      END SUBROUTINE uv3dmix2_geo_tile
#endif
      END MODULE uv3dmix2_geo_mod
