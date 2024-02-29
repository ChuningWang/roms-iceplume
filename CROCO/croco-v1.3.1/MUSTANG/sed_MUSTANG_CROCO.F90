!------------------------------------------------------------------------------
 MODULE sed_MUSTANG_CROCO
!------------------------------------------------------------------------------

#include "cppdefs.h"
#if defined MUSTANG 

!&E=========================================================================
!&E                   ***  MODULE  sed_MUSTANG_CROCO  ***
!&E
!&E ** Purpose : concerns subroutines related to sediment dynamics link to 
!&E              hydrodynamic model to be used in CROCO system
!&E 
!&E ** Description :
!&E     subroutine sed_MUSTANG_settlveloc ! settling velocity in the water 
!&E                                         column
!&E     subroutine sed_gradvit            ! calcul gradient de vitesse, u*
!&E     subroutine sed_skinstress         ! computes the skin stress
!&E     subroutine sed_bottom_slope
!&E     subroutine sedinit_fromfile  ! reads filrepsed where all results of 
!&E                                    sediment dyn. are stored (depends on 
!&E                                    hydro model)
!&E     subroutine sed_exchange_w2s ! MPI treatment of slip deposit fluxes
!&E     subroutine sed_exchange_s2w ! MPI treatment of lateral erosion 
!&E     subroutine sed_exchange_flxbedload ! MPI treatment of bedload fluxes
!&E     subroutine sed_exchange_maskbedload ! MPI exchange of mask for 
!&E                                           bedload
!&E     subroutine sed_exchange_corflu ! MPI treatment of corflu fluxes
!&E     subroutine sed_obc_corflu ! corflu fluxes at boundaries
!&E     subroutine sed_meshedges_corflu ! corflu fluxes interpolation at 
!&E                                       mesh edges
!&E
!&E==========================================================================

#include "coupler_define_MUSTANG.h"

    !! * Modules used
    USE comMUSTANG
    USE comsubstance
    USE module_MUSTANG
    USE module_substance
# if defined key_MUSTANG_flocmod
    USE flocmod, ONLY : f_ws
#endif
    IMPLICIT NONE

    !! * Accessibility 
    PUBLIC sedinit_fromfile
    PUBLIC sed_skinstress
    PUBLIC sed_gradvit
    PUBLIC sed_MUSTANG_settlveloc
#ifdef key_MUSTANG_bedload
    PUBLIC sed_bottom_slope
#if defined MPI 
    PUBLIC sed_exchange_flxbedload
    PUBLIC sed_exchange_maskbedload
#endif
#endif
#if defined MPI  && defined key_MUSTANG_slipdeposit
    PUBLIC sed_exchange_w2s
#endif
#if defined MPI  && defined key_MUSTANG_lateralerosion
    PUBLIC sed_exchange_s2w
#endif
#if defined MUSTANG_CORFLUX
    PUBLIC sed_obc_corflu
    PUBLIC sed_meshedges_corflu
#if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
    PUBLIC sed_exchange_corflu
#endif
#endif


PRIVATE

CONTAINS
 
!!=============================================================================
SUBROUTINE sed_MUSTANG_settlveloc(ifirst, ilast, jfirst, jlast,   &
                                    WATER_CONCENTRATION) 
!&E----------------------------------------------------------------------------
!&E                 ***  ROUTINE sed_MUSTANG_settlveloc  ***
!&E
!&E ** Purpose : settling velocity computation
!&E
!&E ** Description : use arguments and common variable 
!&E  arguments IN : 
!&E         WATER_CONCENTRATION = t : WATER_CONCENTRATION 
!&E  arguments OUT:
!&E         ws_part : settling velocities for CROCO
!&E         ws3_bottom_MUSTANG: settling velocities in  bottom cell
!&E
!&E  need to be know via coupler_define_MUSTANG.h:
!&E         GRAVITY
!&E         kmax=NB_LAYER_WAT
!&E          
!&E  need to be know by code treated substance 
!&E  (if not ==> coupler_MUSTANG.F90)
!&E         imud1, nvpc, nvp, nv_adv, isand1, isand2
!&E         f_ws(iv) (if key_MUSTANG_flocmod)
!&E         ws_free_opt, ws_free_para, ws_free_min, ws_free_max,
!&E         ws_hind_opt, ws_hind_para   
!&E     
!&E  use module MUSTANG variables  :
!&E         ros(iv)
!&E         ws_sand(iv)
!&E        
!&E ** Called by :  MUSTANG_update
!&E
!&E----------------------------------------------------------------------------

!! * Arguments
INTEGER, INTENT(IN) :: ifirst, ilast, jfirst, jlast
REAL(KIND=rsh), DIMENSION(ARRAY_WATER_CONC), INTENT(IN) :: WATER_CONCENTRATION  
!! CROCO : WATER_CONCENTRATION  is directly t 
   
!! * Local declarations
INTEGER                    :: iv, k, ivpc, i, j
REAL(KIND=rsh)             :: cmes, phi, phiv, De, WSfree, Hind
REAL(KIND=rsh), PARAMETER  :: nuw = 0.00000102_rsh

!! * Executable part
DO j = jfirst, jlast
DO i = ifirst, ilast
           
    IF (htot(i, j) > h0fond) THEN
        DO k = 1, NB_LAYER_WAT
            cmes = 0.0_rsh
            DO ivpc = imud1, nvpc
                cmes = cmes + WATER_CONCENTRATION(i, j, k, nstp, itemp + ntrc_salt + ivpc)
            ENDDO
            cmes = MAX(0.0_rsh, cmes)

            ! first calculating sand settling velocity
            DO iv = isand1, isand2
                ws_part(i, j, k, itemp + ntrc_salt + iv) = ws_sand(iv)
            ENDDO
            
            ! next mud settling velocity 
#ifdef key_MUSTANG_flocmod   
            ws_part(i, j, k, itemp + ntrc_salt + imud1 : itemp + ntrc_salt + imud2 ) = f_ws(1:nv_mud)  
#else
            DO iv = imud1, nvp
                ! Free settling velocity - flocculation
                IF(ws_free_opt(iv) == 0) THEN ! constant settling velocity
                    WSfree = ws_free_min(iv)
                ELSEIF (ws_free_opt(iv) == 1) THEN ! Van Leussen 1994
                    WSfree = ws_free_para(1, iv) * cmes**ws_free_para(2, iv) &
                        * (1._rsh + ws_free_para(3, iv) * gradvit(k, i, j)) / &
                        (1._rsh + ws_free_para(4, iv) * gradvit(k, i, j)**2._rsh)
                ELSEIF (ws_free_opt(iv) == 2) THEN ! Winterwerp 1999
                    De = ws_free_para(1, iv) &
                        + ws_free_para(2, iv) &
                        * cmes / (ws_free_para(3, iv) * sqrt(gradvit(k, i, j)))
                    IF (De .GT. sqrt(nuw / gradvit(k, i, j))) THEN 
                        De = sqrt(nuw / gradvit(k, i, j)) 
                        ! in case of large C/low G limit floc size to kolmogorov microscale
                    ENDIF
                    WSfree = (ros(iv) - RHOREF) * GRAVITY / (18._rsh * RHOREF * nuw)  &
                        * ws_free_para(1, iv)**(3._rsh - ws_free_para(4, iv))  &
                        * De**(ws_free_para(4, iv) - 1._rsh)
                ELSEIF (ws_free_opt(iv) == 3) THEN ! Wolanski et al., 1989
                    WSfree = ws_free_para(1, iv) * cmes**ws_free_para(2, iv)
                ENDIF

                ! Hindered settling
                ! if ws_hind_opt.EQ.0 : no hindered settling... Hind = 1
                Hind = 1._rsh
                IF (ws_hind_opt(iv) == 1) THEN ! Scott, 1984
                    phi = MIN(1.0_rsh, cmes / ws_hind_para(1, iv))
                    Hind = (1._rsh - phi)**ws_hind_para(2, iv)
                ELSEIF (ws_hind_opt(iv) == 2) THEN ! Winterwerp, 2002 
                    ! WARNING : ros(iv) must be the same for all MUDS variables
                    phi = cmes / ros(iv)
                    IF (ws_free_opt(iv) == 2) THEN
                        phiv = phi * (De / ws_free_para(1, iv))**(3._rsh - ws_free_para(4, iv))
                    ELSE
                        phiv = cmes / ws_hind_para(1, iv)
                    ENDIF
                    Hind = (1._rsh - phiv)**ws_hind_para(2, iv) * &
                        (1._rsh - phi) / (1._rsh + 2.5_rsh * phiv)
                ELSEIF (ws_hind_opt(iv) == 3) THEN ! wolanski et al., 1989
                    IF (ws_free_opt(iv) == 3) THEN
                        Hind = 1._rsh / (cmes**2._rsh + ws_hind_para(1, iv)**2._rsh)**ws_hind_para(2, iv)
                    ENDIF  ! ws_hind_opt(iv) == 3 only if ws_free_opt(iv) == 3
                ENDIF

                ! limiting ws with min/max values...
                ! necessary if ndt_part not updated during the simulation, 
                ! ndt_part calculated in t3dmix_tridiagonal_settling from max ws
                ws_part(i, j, k, itemp + ntrc_salt + iv) = max(ws_free_min(iv), &
                    min(ws_free_max(iv), WSfree * Hind))
            ENDDO

#endif  /* key_MUSTANG_flocmod */

            DO iv = nvpc+1, nvp
                IF(irkm_var_assoc(iv) < imud1 .AND. irkm_var_assoc(iv) > 0) THEN    
                    ! sorbed substances on sands
                    ws_part(i, j, k, itemp + ntrc_salt + iv) = &
                        ws_part(i, j, k, itemp + ntrc_salt + irkm_var_assoc(iv))
                ENDIF
            ENDDO
            
            DO iv = nvp+1, nv_adv
                ws_part(i, j, k, itemp + ntrc_salt + iv) = 0.0_rsh
            ENDDO

        ENDDO ! do loop on k
        ws3_bottom_MUSTANG(1:nvp, i, j) = ws_part(i, j, 1, itemp+ntrc_salt+1:itemp+ntrc_salt+nvp)

    ELSE  ! htot(i,j) <= h0fond
        ws_part(i, j, :, :) = 0.0_rsh
        ws3_bottom_MUSTANG(:, i, j) = 0.0_rsh
    ENDIF

ENDDO
ENDDO

END SUBROUTINE sed_MUSTANG_settlveloc     

!!=============================================================================
SUBROUTINE sed_gradvit(ifirst, ilast, jfirst, jlast)
!&E--------------------------------------------------------------------------
!&E                 ***  ROUTINE sed_gradvit  ***
!&E
!&E ** Purpose : calculation of the turbulence energy G  
!&E
!&E ** Description : G= sqrt(turbulence dissipation/viscosity)
!&E                 to be programmed using hydrodynamic knowledge
!&E           using htot, RHOREF, sig, epn, nz ..
!&E
!&E     output : gradvit (in comMUSTANG)
!&E
!&E ** Called by :  MUSTANG_update
!&E
!&E--------------------------------------------------------------------------
!! * Modules used
#if defined key_MUSTANG_flocmod && defined SED_TOY_FLOC_0D
    USE flocmod, ONLY : flocmod_comp_g
#endif
#  include "mixing.h"
#  include "ocean3d.h"

!! * Arguments 
INTEGER, INTENT(IN)       :: ifirst, ilast, jfirst, jlast

!! * Local declarations
INTEGER        :: i, j, k
REAL(KIND=rsh) :: dist_surf_on_bottom, nuw, diss

nuw = 1.0e-6

DO j = jfirst, jlast
DO i = ifirst, ilast
    IF(htot(i, j) .GT. h0fond)  THEN
    DO k = 1, N
#if defined key_MUSTANG_flocmod && defined SED_TOY_FLOC_0D
        call flocmod_comp_g(gradvit(k, i, j), time-time_start)
#else
#if defined GLS_MIXING
        !
        ! Dissipation from turbulence clossure
        !
        if (k.eq.1) then
            diss = Eps_gls(i,j,k)
        elseif (k.eq.N) then
            diss = Eps_gls(i,j,k-1)
        else
            diss = 0.5*(Eps_gls(i,j,k-1)+Eps_gls(i,j,k))
        endif
        gradvit(k, i, j) = sqrt(diss/nuw)
#else
    ! gradvit : G=sqrt( turbulence dissipation rate/ vertical viscosity coefficient)
    ! if  turbulence dissipation rate has not been already evaluated: 
    ! use empirical formula from   Nezu and Nakawaga (1993)
    ! turbulence dissipation_rate = ustarbot**3 /Karman/Htot * (distance from surface/distance from bottom)
        dist_surf_on_bottom = ((z_w(i, j, N) - z_r(i, j, k)) / (z_r(i, j, k) - z_w(i, j, 0)))
        gradvit(k, i, j) = sqrt(ustarbot(i, j)**3._rsh / 0.4_rsh / htot(i, j) / &
                        (nuw + epsilon_MUSTANG) * dist_surf_on_bottom) 
#endif
#endif
    END DO
    ENDIF
ENDDO
ENDDO

END SUBROUTINE sed_gradvit

!!==============================================================================
  SUBROUTINE sed_skinstress(ifirst, ilast, jfirst, jlast)
  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE sed_skinstress  ***
  !&E
  !&E ** Purpose : computes  bottom shear stress
  !&E
  !&E ** Description :  
  !&E 
  !&E Compute bottom skin-friction stress due to combined maximum wave and current 
  !&E interaction
  !&E 
  !&E Available options to compute tauskin (combine current + wave (if WAVE_OFFLINE 
  !&E is defined)) and tauskin_x/tauskin_y (components / rho) :
  !&E - default : Soulsby formulation (with z0sed)
  !&E - BBL : d50 is constant (160microns) in the bustrw/bvstrw computation see 
  !&E bbl.F (and optionnal key_tauskin* do not apply)
  !&E 
  !&E Available options to compute tauskin_c
  !&E - default : compute at u,v points and then compute at rho point 
  !&E     with ubar if key_tauskin_c_ubar is defined, with bottom u(1) if not, 
  !&E     this case use 12 u,v points to compute tauc at rho point
  !&E     if key_tauskin_c_upwind is defined, x/y component used are tauskin 
  !&E     computed at uv point upwind from current
  !&E - key_tauskin_c_center : compute at rho point from immediate u,v value 
  !&E     with ubar if key_tauskin_c_ubar is definde, with bottom u(1) if not
  !&E
  !&E ** Called by :  MUSTANG_update
  !&E
  !&E--------------------------------------------------------------------------
   !! * Modules used
# ifdef WAVE_OFFLINE
   USE module_substance, ONLY : Uwave, Dwave, Pwave
# endif

#  ifdef BBL
#  include "bbl.h"
#  endif

#  include "grid.h"
#  include "ocean3d.h"
#  include "ocean2d.h"


  ! Arguments
  INTEGER, INTENT(IN)  :: ifirst, ilast, jfirst, jlast

  ! Local declaration
  INTEGER         :: i, j, k
  REAL(KIND=rsh)  :: speed ! current speed
  REAL(KIND=rsh)  :: speedu, speedv ! current speed*u, current speed*v, (m2.s-2)
  REAL(KIND=rsh)  :: urho, vrho ! current speed component (m.s-1)
  REAL(KIND=rsh),DIMENSION(ifirst-1: ilast+1, jfirst-1: jlast+1) :: Zr ! cell height (m)
  REAL(KIND=rsh),DIMENSION(ifirst-1: ilast+1, jfirst-1: jlast+1) :: ustar2_u ! (m2.s-2), compute at u point
  REAL(KIND=rsh),DIMENSION(ifirst-1: ilast+1, jfirst-1: jlast+1) :: ustar2_v ! (m2.s-2), compute at v point
# ifdef WAVE_OFFLINE
  REAL(KIND=rsh)  :: fws2ij, speedbar, alpha, beta, cosamb, sinamb, tauskin_cw
# endif
# ifdef key_tauskin_c_upwind   
  REAL(KIND=rsh)  :: cff1, cff2, cff3, cff4  ! temporary coefficients
# endif  

  ! Executable part


#  ifdef BBL /*warning, d50 is constant (160microns) in the bustrw/bvstrw computation see bbl.F */
  do j = jfirst, jlast
    do i = ifirst, ilast
      tauskin(i, j) = sqrt( bustrw(i, j)**2 + bvstrw(i, j)**2) * RHOREF
#  ifdef WET_DRY AND MASKING
      tauskin(i, j) = tauskin(i, j) * rmask_wet(i, j)
#  endif
#  ifdef key_MUSTANG_bedload
      urho = 0.5 * (u(i, j, 1, nnew) + u(i+1, j, 1, nnew))
      vrho = 0.5 * (v(i, j, 1, nnew) + v(i, j+1, 1, nnew))
      speed = SQRT(urho**2 + vrho**2)
      tauskin_c(i, j) = tauskin(i, j) ! used in eval_bedload
      tauskin_x(i, j) = urho / (speed + epsilon_MUSTANG) * tauskin_c(i, j)
      tauskin_y(i, j) = vrho / (speed + epsilon_MUSTANG) * tauskin_c(i, j)
      do j = jfirst, jlast
        do i = ifirst, ilast+1
          raphbx(i, j) = ABS(u(i, j, 1, nnew)) / (ABS(u(i, j, 1, nnew)) + epsilon_MUSTANG)
        enddo
      enddo  
      do j = jfirst , jlast+1
        do i = ifirst, ilast
          raphby(i, j) = ABS(v(i, j, 1, nnew)) / (ABS(v(i, j, 1, nnew)) + epsilon_MUSTANG)
        enddo
      enddo 
#  endif /* key_MUSTANG_bedload */

#  else /* else on #ifdef BBL */

      do j = jfirst-1, jlast+1
        do i = ifirst-1, ilast+1
          Zr(i, j) =  max(z_r(i, j, 1) - z_w(i, j, 0), z0sed(i, j) + 1.E-4)
        enddo
      enddo

# ifdef key_tauskin_c_center
  do j = jfirst, jlast
    do i = ifirst, ilast
#ifdef key_tauskin_c_ubar
#ifdef MPI
      if ((float(j + jj * Mm) .GE. 0) .and. (float(i + ii * Lm) .GE. 0) )then
#else
      if ((float(j       ) .GE. 0) .and. (float(i       ) .GE. 0)) then
#endif
        urho = 0.5 * (ubar(i, j, nnew) + ubar(i+1, j, nnew))
        vrho = 0.5 * (vbar(i, j, nnew) + vbar(i, j+1, nnew))  
      else
        urho = 0.
        vrho = 0.
      endif  
      speed = SQRT(urho**2 + vrho**2)
      tauskin_c(i, j) = 0.16_rsh * (LOG( ( z_w(i, j, N) - z_w(i, j, 0)) /   &
                    (z0sed(i, j) * 2.718)))**(-2) * speed**2                &
                    * (rho(i, j, 1) + rho0)
# else
      urho = 0.5 * (u(i, j, 1, nnew) + u(i+1, j, 1, nnew))
      vrho = 0.5 * (v(i, j, 1, nnew) + v(i, j+1, 1, nnew))
      speed = SQRT(urho**2 + vrho**2)
      tauskin_c(i, j) = 0.16_rsh * (LOG( ( Zr(i, j) ) /   &
                    (z0sed(i, j))))**(-2) * speed**2                &
                    * (rho(i, j, 1) + rho0)
# endif
#  ifdef key_MUSTANG_bedload     
      tauskin_x(i, j) = urho / (speed + epsilon_MUSTANG) * tauskin_c(i, j)
      tauskin_y(i, j) = vrho / (speed + epsilon_MUSTANG) * tauskin_c(i, j)
      do j = jfirst, jlast
        do i = ifirst, ilast+1
          raphbx(i, j) = ABS(u(i, j, 1, nnew)) / (ABS(u(i, j, 1, nnew)) + epsilon_MUSTANG)
        enddo
      enddo  
      do j = jfirst , jlast+1
        do i = ifirst, ilast
          raphby(i, j) = ABS(v(i, j, 1, nnew)) / (ABS(v(i, j, 1, nnew)) + epsilon_MUSTANG)
        enddo
      enddo 
# endif

# else /* else on #ifdef key_tauskin_c_center */

  do j = jfirst, jlast
    do i = ifirst, ilast+1
      raphbx(i, j) = ABS(u(i, j, 1, nnew)) / (ABS(u(i, j, 1, nnew)) + epsilon_MUSTANG)
#ifdef key_tauskin_c_ubar
      speedu = SQRT(0.0625_rsh * (vbar(  i, j, nnew) + vbar(   i, j+1, nnew) +   &
                                  vbar(i-1, j, nnew) + vbar( i-1, j+1, nnew))**2 &
                                + ubar(  i, j, nnew)**2) * ubar(i, j, nnew)
#ifdef MPI
      if (float(i + ii * Lm) .GE. 0) then
#else
      if (float(i          ) .GE. 0) then
#endif
        ustar2_u(i, j) = 0.16_rsh * (LOG( ( 0.5* (        &
              (z_w(i-1, j, N) - z_w(i-1, j, 0)) + (z_w(i, j, N) - z_w(i, j, 0))))  /   &
              ((0.5 * (z0sed(i-1, j) + z0sed(i, j))) 
              * 2.718)))**(-2) * speedu 
      else
        ustar2_u(i, j) = 0.
      endif
#else
      speedu = SQRT(0.0625_rsh*(v(  i, j, 1, nnew) + v(  i, j+1, 1, nnew) +   &
                                v(i-1, j, 1, nnew) + v(i-1, j+1, 1, nnew))**2 &
                              + u(  i, j, 1, nnew)**2) * u(i, j, 1, nnew)
      ustar2_u(i, j) = 0.16_rsh * (LOG(0.5 * (Zr(i-1, j) + Zr(i, j)) /   &
                      (0.5 * (z0sed(i-1, j) + z0sed(i, j)))              &
                        ))**(-2) * speedu 
#endif
    enddo
  enddo  
  DO j = jfirst, jlast+1
    DO i = ifirst, ilast
      raphby(i, j) = ABS(v(i, j, 1, nnew)) / (ABS(v(i, j, 1, nnew)) + epsilon_MUSTANG)
#ifdef key_tauskin_c_ubar
      speedv = SQRT(0.0625_rsh*(ubar(i,   j, nnew) + ubar(i+1,   j, nnew) +   &
                                ubar(i, j-1, nnew) + ubar(i+1, j-1, nnew))**2 &
                              + vbar(i, j, nnew)**2) * vbar(i, j, nnew) 
#ifdef MPI
      if (float(j + jj * Mm) .GE. 0) then
#else
      if (float(j          ) .GE. 0) then
#endif
        ustar2_v(i, j) = 0.16_rsh * (LOG( (0.5* (        &
              (z_w(i, j-1, N) - z_w(i, j-1, 0)) + (z_w(i, j, N) - z_w(i, j, 0)))) /   &
              ((0.5 * (z0sed(i, j-1) + z0sed(i, j)))                          &
              * 2.718)))**(-2) * speedv
      else
        ustar2_v(i, j) = 0.
      endif
#else
      speedv = SQRT(0.0625_rsh * (u(i,   j, 1, nnew) + u(i+1,   j, 1, nnew) +   &
                                  u(i, j-1, 1, nnew) + u(i+1, j-1, 1, nnew))**2 &
                                + v(i,   j, 1, nnew)**2) * v(i, j, 1, nnew)
      ustar2_v(i,j) = 0.16_rsh * (LOG(0.5 * (Zr(i, j-1) + Zr(i, j)) /   &
                      (0.5 * (z0sed(i, j-1) + z0sed(i, j)))              &
                      ))**(-2) * speedv        
#endif
    enddo
  enddo
  
  do j = jfirst, jlast
    do i = ifirst, ilast
# ifdef key_tauskin_c_upwind
      cff1 = 0.5 * (1.0 + SIGN(1.0, ustar2_u(i+1, j)))
      cff2 = 0.5 * (1.0 - SIGN(1.0, ustar2_u(i+1, j)))
      cff3 = 0.5 * (1.0 + SIGN(1.0, ustar2_u(i  , j)))
      cff4 = 0.5 * (1.0 - SIGN(1.0, ustar2_u(i  , j)))
      tauskin_x(i, j) = cff3 * (cff1 * ustar2_u(i, j) +                     &
                        cff2 * 0.5 * (ustar2_u(i, j) + ustar2_u(i+1, j))) + &
                        cff4 * (cff2 * ustar2_u(i+1, j) +                   &
                        cff1 * 0.5 * (ustar2_u(i, j) + ustar2_u(i+1, j)))

      cff1 = 0.5 * (1.0 + SIGN(1.0, ustar2_v(i, j+1)))
      cff2 = 0.5 * (1.0 - SIGN(1.0, ustar2_v(i, j+1)))
      cff3 = 0.5 * (1.0 + SIGN(1.0, ustar2_v(i, j)))
      cff4 = 0.5 * (1.0 - SIGN(1.0, ustar2_v(i, j)))
      tauskin_y(i, j) = cff3 * (cff1 * ustar2_v(i, j) +                     &
                        cff2 * 0.5 * (ustar2_v(i, j) + ustar2_v(i, j+1))) + &
                        cff4 * (cff2 * ustar2_v(i, j+1) +                   &
                        cff1 * 0.5 * (ustar2_v(i, j) + ustar2_v(i, j+1)))
# else
      tauskin_x(i, j) = (ustar2_u(i+1, j) * raphbx(i+1, j)+ &
                          ustar2_u(  i, j) * raphbx(  i, j)) / &
                          ( raphbx(i+1, j) + raphbx(  i, j) + epsilon_MUSTANG)

      tauskin_y(i, j) = (ustar2_v(i,   j) * raphby(i,   j) + &
                          ustar2_v(i, j+1) * raphby(i, j+1)) / &
                          ( raphby(i,   j) + raphby(i, j+1) + epsilon_MUSTANG)
# endif
      tauskin_c(i, j) = SQRT(tauskin_x(i, j)**2 + tauskin_y(i, j)**2) * (rho(i, j, 1) + rho0)
      tauskin_x(i, j) = tauskin_x(i, j) * (rho(i, j, 1) + rho0)
      tauskin_y(i, j) = tauskin_y(i, j) * (rho(i, j, 1) + rho0)

# endif /* end of #ifdef key_tauskin_c_center */

! combined wave and current interaction
! Note : tauskin direction is current direction (needed if bedload)
# ifdef WAVE_OFFLINE  
      fws2ij = fws2
      if (l_fricwave .AND. Uwave(i,j)*Pwave(i,j) > 0.0_rsh .AND.  &
                Pwave(i,j) > 0.001_rsh .AND. Uwave(i,j) > 0.001_rsh) then
        fws2ij = 0.5_rsh * 1.39_rsh * (Uwave(i,j) * Pwave(i,j) / &
                REAL(2.0_rlg * pi * z0sed(i, j),rsh))**(-0.52_rsh)
      endif

      ! calculation of shear stress due to waves tauskin_w 
      ! --------------------------
      tauskin_w(i, j) = ((rho(i, j, 1) + rho0) * fws2ij * Uwave(i, j)**2)

      speedbar = SQRT(ubar(i, j, nnew)**2 + vbar(i, j,nnew)**2)
      if (tauskin_c(i, j) > 0.0_rsh .AND. tauskin_w(i, j) > 0.0_rsh .AND. speedbar > 0.0_rsh) then

      ! calculation of shear stress with the formula of Soulsby (1995)
      ! ======================================================
      ! calculation of  tauskin_c (current) influenced by the waves
      ! ----------------------------------------------------
        tauskin_cw = tauskin_c(i, j) * (1 + (1.2 * (tauskin_w(i, j) / &
                     (tauskin_w(i, j) + tauskin_c(i, j)))**3.2))

      ! calculating the difference in angle between the direction of the waves and the current
      ! ---------------------------------------------------------------------------
      ! calculating the direction of the current relative to the north
        alpha = ACOS(vbar(i, j, nnew) / speedbar)   ! in radians
      ! calculation of wave orientation relative to north
        beta = Dwave(i, j)   ! beta and Dwave in radians
      ! calculation of cos(alpha-beta) and sin(alpha-beta)
        cosamb = ABS( COS(alpha - beta))
        sinamb = ABS( SIN(alpha - beta))

      ! calculation of tauskin (waves + current)
      ! -----------------------------------
        tauskin(i, j) = SQRT( (tauskin_cw + tauskin_w(i, j) * cosamb)**2 + &
                                           (tauskin_w(i, j) * sinamb)**2 )
    else
      tauskin(i, j) = tauskin_w(i, j) + tauskin_c(i, j)
    endif

# else /* WAVE_OFFLINE */
    tauskin(i, j) = tauskin_c(i, j)
# endif

# if defined WET_DRY && defined MASKING 
    tauskin(i, j) = tauskin(i, j) * rmask_wet(i, j)
# endif

# endif  /* end of #ifdef BBL */

      ! u* = ustarbot    
      if (htot(i, j) .GT. h0fond) then
        if (tauskin(i, j) < 0.) then
          ustarbot(i, j) = 0.0_rsh
        else
          ustarbot(i, j) = (tauskin(i, j) / RHOREF)**0.5_rsh
        endif
      endif
    enddo
  enddo

  END SUBROUTINE sed_skinstress


!!==============================================================================
#ifdef key_MUSTANG_bedload
  SUBROUTINE sed_bottom_slope(ifirst, ilast, jfirst, jlast, bathy)
   !&E--------------------------------------------------------------------------                         
   !&E                 ***  ROUTINE sed_bottom_slope  ***
   !&E
   !&E ** Purpose : evaluation of bottom slope
   !&E
   !&E ** Description : depend on host model grid, compute slope_dhdx and 
   !&E    slope_dhdy from bathy of neigbour cells if they are not masked
   !&E
   !&E ** Called by :  MUSTANG_update, MUSTANG_morpho 
   !&E
   !&E--------------------------------------------------------------------------
   !! * Arguments
   INTEGER, INTENT(IN)  :: ifirst, ilast, jfirst, jlast
   REAL(KIND=rsh),DIMENSION(ARRAY_BATHY_H0),INTENT(IN)  :: bathy  ! bathymetry (m)

   !! * Local declarations
   INTEGER :: i, j

      DO j = jfirst, jlast
        DO i = ifirst, ilast
          IF (bathy(i+1, j).LE. -valmanq .OR. bathy(i-1, j).LE. -valmanq) then
             slope_dhdx(i, j) = 0.0_rsh
          ELSE
             slope_dhdx(i, j) = -1.0_rsh*(-bathy(i+1, j)+bathy(i-1, j)) / (2.0_rsh * CELL_DX(i, j))
          ENDIF
          IF (bathy(i, j+1).LE. -valmanq .OR. bathy(i, j-1).LE. -valmanq) then
             slope_dhdy(i, j) = 0.0_rsh
          ELSE
             slope_dhdy(i, j) = -1.0_rsh*(-bathy(i, j+1)+bathy(i, j-1)) / (2.0_rsh * CELL_DY(i, j))
          ENDIF
        ENDDO
      ENDDO

  END SUBROUTINE sed_bottom_slope
!!=============================================================================
#if defined MPI
  SUBROUTINE sed_exchange_flxbedload(ifirst, ilast, jfirst, jlast)
    !&E--------------------------------------------------------------------------                         
    !&E                 ***  ROUTINE sed_exchange_flxbedload  ***
    !&E
    !&E ** Purpose : treatment of bedload fluxes
    !&E
    !&E ** Description : MPI exchange between processors
    !&E
    !&E ** Called by :  MUSTANG_update
    !&E
    !&E--------------------------------------------------------------------------
    ! * Arguments
    INTEGER, INTENT(IN)  :: ifirst, ilast, jfirst, jlast

    !! * Local declarations
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY) :: workexch
    INTEGER :: iv

    do iv=ibedload1,ibedload2
        workexch(:,:) = flx_bx(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
             &          workexch(START_2D_ARRAY))
        flx_bx(iv,:,:) = workexch(:,:)
   
        workexch(:,:) = flx_by(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
             &          workexch(START_2D_ARRAY))
        flx_by(iv,:,:) = workexch(:,:)
    enddo

  END SUBROUTINE sed_exchange_flxbedload
#endif /* MPI */
#endif /* key_MUSTANG_bedload */

!!=============================================================================
  SUBROUTINE sedinit_fromfile(BATHY_H0)
 
   !&E-------------------------------------------------------------------------
   !&E                 ***  ROUTINE sedinit_fromfile  ***
   !&E
   !&E ** Purpose : manages the fields to be re-initialized in case of the 
   !&E              continuation of a previous run
   !&E
   !&E ** Description : open and read a netcdf file, written during a previous run 
   !&E
   !&E ** Called by :  MUSTANG_init 
   !&E
   !&E-------------------------------------------------------------------------
   !! * Modules used
    implicit none

    !! * Arguments
    REAL(KIND=rsh),DIMENSION(ARRAY_BATHY_H0),INTENT(IN) :: BATHY_H0                         

# include "netcdf.inc"
    real time_scale
    integer iv, k, itrc, indWrk
    integer ncid, indx, varid,  ierr, lstr, lvar, latt, lenstr,       &
    start(2), count(2), ibuff(6), nf_fread, checkdims
    character units*180, nomcv*30
    character inised_name*180
    real tmp(GLOBAL_2D_ARRAY)
    real tmp3d(GLOBAL_2D_ARRAY, ksdmin:ksdmax)

    tmp(PROC_IN_ARRAY) = 0
    tmp3d(PROC_IN_ARRAY,ksdmin:ksdmax) = 0
    z0sed(PROC_IN_ARRAY) = z0seduni
    dzsmax(PROC_IN_ARRAY) = dzsmaxuni
    ksmi(PROC_IN_ARRAY) = ksmiuni
    ksma(PROC_IN_ARRAY) = 0
    hsed(PROC_IN_ARRAY) = -valmanq
    dzs(ksdmin:ksdmax,PROC_IN_ARRAY) = -valmanq
    cv_sed(-1:nv_tot,ksdmin:ksdmax,PROC_IN_ARRAY) = -valmanq
    c_sedtot(ksdmin:ksdmax,PROC_IN_ARRAY) = -valmanq
!
! Open initial conditions netCDF file for reading. Check that all
! spatial dimensions in that file are consistent with the model
! arrays, determine how many time records are available in the file
! and set record from which the dada will be read.
!
! The record is set as follows: (1) if only one time record is
! available in the file, then that record is used REGARDLESS of
! value of nrrec supplied from the parameter file; (2) if the
! file has multiple records and nrrec is positive, then nrrec is
! used, provided that nrrec is within the available records; 
! (3) if the file has multiple records and nrrec<0, then THE LAST 
! available record is used.

!      if (may_day_flag.ne.0) return      !-->  EXIT
      inised_name = filrepsed
      lstr = lenstr(inised_name)
      ierr = nf_open(inised_name(1:lstr), nf_nowrite, ncid)
      if (ierr.eq.nf_noerr) then
        ierr = checkdims(ncid, inised_name, lstr, indx)

        if (ierr.ne.nf_noerr) then
         goto 99
        elseif (indx.eq.0) then
          indx = 1
        elseif (indx.gt.0 .and. nrrec.gt.0 .and. nrrec.le.indx) then
          indx = nrrec
        elseif (indx.gt.0 .and. nrrec.gt.indx) then
          write(stdout,'(/1x,A,I4,A/16x,A,I4,A/16x,3A/)')   &
                 'SEDINIT_FROMFILE ERROR: requested restart time record', &
                  nrrec, ' exceeds',  'number of available records', &
                  indx,'in netCDF file', '''',ininame(1:lstr),'''.'
          goto 99                                        !--> ERROR
        endif
      else
        write(stdout,'(/1x,2A/15x,3A)') 'SEDINIT_FROMFILE ERROR: Cannot ', &
                    'open netCDF file', '''', ininame(1:lstr) ,'''.'
        goto 99                                           !--> ERROR
      endif


! ksmi, ksma
      ierr=nf_inq_varid (ncid,'ksmi', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (tmp, ncid, varid, indx, 0)
        ksmi(PROC_IN_ARRAY)=INT(tmp(PROC_IN_ARRAY))
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) 'ksmi', indx, inised_name(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) 'ksmi', inised_name(1:lstr)
        goto 99                                           !--> ERROR
      endif
      
       ierr=nf_inq_varid (ncid,'ksma', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (tmp, ncid, varid, indx, 0)
        ksma(PROC_IN_ARRAY)=INT(tmp(PROC_IN_ARRAY))

        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) 'ksma', indx, inised_name(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) 'ksma', inised_name(1:lstr)
        goto 99                                           !--> ERROR
      endif

     WHERE (BATHY_H0(PROC_IN_ARRAY) == -valmanq) ksmi(PROC_IN_ARRAY)=1
     WHERE (BATHY_H0(PROC_IN_ARRAY) == -valmanq) ksma(PROC_IN_ARRAY)=0
!  DZS
      ierr=nf_inq_varid (ncid,'DZS', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (tmp3d, ncid, varid, indx, 12)
         do k=ksdmin,ksdmax
            dzs(k,:,:)=tmp3d(:,:,k)
         enddo
        
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) 'DZS', indx, inised_name(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) 'DZS', inised_name(1:lstr)
        goto 99                                           !--> ERROR
      endif


! CVSED
 
      c_sedtot(:,:,:)=0.0_rsh
      do iv=-1,nv_tot 
       
       indWrk=indxT+ntrc_salt+ntrc_substot+iv+4+2
       nomcv=vname(1,indWrk)

       ierr=nf_inq_varid (ncid,nomcv, varid)
       if (ierr .eq. nf_noerr) then
        ierr=nf_fread (tmp3d, ncid, varid, indx, 12)

        do k=ksdmin,ksdmax
            cv_sed(iv,k,:,:)=tmp3d(:,:,k)
            c_sedtot(k,:,:)=c_sedtot(k,:,:)+cv_sed(iv,k,:,:)*typart(iv)  
        enddo

        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) nomcv, indx, inised_name(1:lstr)
          goto 99                                         !--> ERROR
        endif
       
       else
        if (iv > nvpc .and. iv < nvp+1 ) then
         ! not constitutive particulate  variables (Initial in M/Msed converted to M/m3 sed)
         IF (irkm_var_assoc(iv) >0) THEN
           do k=ksdmin,ksdmax
             cv_sed(iv,k,:,:)=cini_sed(iv)*cv_sed(irkm_var_assoc(iv),k,:,:)
           enddo
         ELSE
           do k=ksdmin,ksdmax
            cv_sed(iv,k,:,:)=cini_sed(iv)*c_sedtot(k,:,:)
           enddo
         END IF
         MPI_master_only  write(stdout,3) nomcv, inised_name(1:lstr)
        else if (iv > nvp) then
          ! dissolved variables (M/m3 EI)
          do k=ksdmin,ksdmax
            cv_sed(iv,k,:,:)=cini_sed(iv)             
          enddo
         MPI_master_only  write(stdout,3) nomcv,vname(1,indxT+ntrc_salt+iv), &
                                          inised_name(1:lstr)
        else
         MPI_master_only  write(stdout,1) nomcv, inised_name(1:lstr)
         goto 99                          !--> ERROR
        endif
       endif
      enddo
     
     do k=ksdmin,ksdmax
     WHERE (BATHY_H0(PROC_IN_ARRAY) == -valmanq) c_sedtot(k,PROC_IN_ARRAY)=-valmanq
     enddo

! Close input NetCDF file.
!
      ierr=nf_close(ncid)

  1   format(/1x,'SEDINIT_FROMFILE - unable to find variable:',    1x,A,  &
                                 /15x,'in input NetCDF file:',1x,A/)
  2   format(/1x,'SEDINIT_FROMFILE - error while reading variable:',1x, A, &
         2x,'at time record =',i4/15x,'in input NetCDF file:',1x,A/)
  3   format(/1x,'SEDINIT_FROMFILE - unable to find variable:',    1x,A,/10x,A, &
     &                            /15x,'in input NetCDF file:',1x,A,  &
         1x,'-> analytical value (cini_sed)'/)
      return
  99  may_day_flag=2
      return
      
  END SUBROUTINE sedinit_fromfile
!!=============================================================================

#if defined MPI && defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
    SUBROUTINE sed_exchange_maskbedload(ifirst, ilast, jfirst, jlast)
    !&E-------------------------------------------------------------------------
    !&E                 ***  ROUTINE sed_exchange_maskbedload ***
    !&E
    !&E ** Purpose : exchange MPI mask for bedload
    !&E
    !&E ** Description : exchange MPI mask for bedload
    !&E
    !&E ** Called by : MUSTANG_update
    !&E-------------------------------------------------------------------------

    !! * Arguments
    INTEGER,INTENT(IN) :: ifirst, ilast, jfirst, jlast

    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY) :: workexch
    workexch(:,:) = sedimask_h0plusxe(:,:)
    call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
              &           workexch(START_2D_ARRAY))
    sedimask_h0plusxe(:,:) = workexch(:,:)

    END SUBROUTINE sed_exchange_maskbedload
#endif /* defined MPI && defined key_MUSTANG_V2 && defined key_MUSTANG_bedload */

!!=============================================================================

#if defined MPI && defined key_MUSTANG_slipdeposit
    SUBROUTINE sed_exchange_w2s(ifirst, ilast, jfirst, jlast)
    !&E-------------------------------------------------------------------------
    !&E                 ***  ROUTINE sed_exchange_w2s ***
    !&E
    !&E ** Purpose : MPI exchange of slip deposit flux between processors
    !&E
    !&E ** Description : MPI exchange between processors
    !&E      used only if slopefac .NE. 0 (slip deposit if steep slope)
    !&E
    !&E ** Called by : MUSTANG_update
    !&E-------------------------------------------------------------------------

    !! * Arguments
    INTEGER,INTENT(IN) :: ifirst, ilast, jfirst, jlast

    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY) :: workexch
    INTEGER :: iv
      
    do iv = isand2+1, nvp
        workexch(:,:) = flx_w2s_corim1(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
                &          workexch(START_2D_ARRAY))
        flx_w2s_corim1(iv,:,:) = workexch(:,:)

        workexch(:,:) = flx_w2s_corip1(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
                &          workexch(START_2D_ARRAY))
        flx_w2s_corip1(iv,:,:) = workexch(:,:)

        workexch(:,:) = flx_w2s_corjm1(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
                &          workexch(START_2D_ARRAY))
        flx_w2s_corjm1(iv,:,:) = workexch(:,:)

        workexch(:,:) = flx_w2s_corjp1(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
                &          workexch(START_2D_ARRAY))
        flx_w2s_corjp1(iv,:,:) = workexch(:,:)

        workexch(:,:) = flx_w2s_corin(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
                &          workexch(START_2D_ARRAY))
        flx_w2s_corin(iv,:,:) = workexch(:,:)
    enddo
  
    END SUBROUTINE sed_exchange_w2s
#endif /* defined MPI && defined key_MUSTANG_slipdeposit */

!!=============================================================================

#if defined MPI && defined key_MUSTANG_lateralerosion
    SUBROUTINE sed_exchange_s2w(ifirst, ilast, jfirst, jlast)
    !&E-------------------------------------------------------------------------
    !&E                 ***  ROUTINE sed_exchange_s2w ***
    !&E
    !&E ** Purpose : MPI exchange of lateral erosion flux between processors
    !&E
    !&E ** Description : MPI exchange between processors
    !&E      used only if scoef_erolat .NE. 0 
    !&E
    !&E ** Called by : MUSTANG_update
    !&E-------------------------------------------------------------------------

    !! * Arguments
    INTEGER,INTENT(IN) :: ifirst, ilast, jfirst, jlast

    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY) :: workexch
    INTEGER :: iv
      
    do iv=-1,nv_adv
        workexch(:,:) = flx_s2w_corip1(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
                &          workexch(START_2D_ARRAY))
        flx_s2w_corip1(iv,:,:) = workexch(:,:)

        workexch(:,:) = flx_s2w_corim1(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
                &          workexch(START_2D_ARRAY))
        flx_s2w_corim1(iv,:,:) = workexch(:,:)

        workexch(:,:) = flx_s2w_corjp1(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
                &          workexch(START_2D_ARRAY))
        flx_s2w_corjp1(iv,:,:) = workexch(:,:)

        workexch(:,:) = flx_s2w_corjm1(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
                &          workexch(START_2D_ARRAY))
        flx_s2w_corjm1(iv,:,:) = workexch(:,:)
    enddo
  
    END SUBROUTINE sed_exchange_s2w
#endif /* defined MPI && defined key_MUSTANG_lateralerosion */

!!=============================================================================



#if defined MUSTANG_CORFLUX
#if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
  SUBROUTINE sed_exchange_corflu(ifirst, ilast, jfirst, jlast, type)
   !&E-------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_exchange_corflu ***
   !&E
   !&E ** Purpose : treatment of horizontal flow corrections for the transport 
   !&E              of sand in suspension
   !&E
   !&E ** Description : periodic borders and MPI exchange between processors
   !&E
   !&E ** Called by : MUSTANG_update
   !&E-------------------------------------------------------------------------

    !! * Arguments
    INTEGER,INTENT(IN) :: ifirst, ilast, jfirst, jlast, type

    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY) :: workexch
    INTEGER :: iv
    
    if (type .eq. 0) then  ! corflux and corfluy are still at rho point
        do iv = isand1, isand2
            workexch(:, :) = corflux(iv, :, :)
            call exchange_r2d_tile (ifirst, ilast, jfirst, jlast,  &
                &          workexch(START_2D_ARRAY))
            corflux(iv, :, :) = workexch(:, :)
    
            workexch(:, :) = corfluy(iv, :, :)
            call exchange_r2d_tile (ifirst, ilast, jfirst, jlast,  &
                &          workexch(START_2D_ARRAY))
            corfluy(iv, :, :) = workexch(:, :)
        enddo
    else  ! corflux and corfluy are at u,v point
        do iv = isand1, isand2
            workexch(:, :) = corflux(iv, :, :)
            call exchange_u2d_tile (ifirst, ilast, jfirst, jlast,  &
                &          workexch(START_2D_ARRAY))
            corflux(iv, :, :) = workexch(:, :)
    
            workexch(:, :) = corfluy(iv, :, :)
            call exchange_v2d_tile (ifirst, ilast, jfirst, jlast,  &
                &          workexch(START_2D_ARRAY))
            corfluy(iv, :, :) = workexch(:, :)
        enddo
    endif

   END SUBROUTINE sed_exchange_corflu
#endif /* defined EW_PERIODIC || defined NS_PERIODIC || defined MPI */
!!=============================================================================
   
   SUBROUTINE sed_obc_corflu(ifirst, ilast, jfirst, jlast)
 
    !&E------------------------------------------------------------------------
    !&E                 ***  ROUTINE sed_obc_corflu ***
    !&E
    !&E ** Purpose : treatment of horizontal flow corrections for the transport 
    !&E              of sand in suspension
    !&E
    !&E ** Description : extrapolation at borders
    !&E
    !&E ** Called by : MUSTANG_update
    !&E--------------------------------------------------------------------------
 
    !! * Arguments
    INTEGER,INTENT(IN) :: ifirst, ilast, jfirst, jlast
 
    !! * Local declarations
    INTEGER :: i, j, ivp

    do ivp = isand1, isand2

    !! * Executable part
#if defined MPI 
    if (float(ifirst + ii * Lm) .EQ. IMIN_GRID) then
#else
    if (float(ifirst) .EQ. IMIN_GRID) then
#endif
        corflux(ivp, ifirst, :) = corflux(ivp, ifirst+1, :)
        corfluy(ivp, ifirst, :) = corfluy(ivp, ifirst+1, :)
        corflux(ivp, ifirst-1, :) = corflux(ivp, ifirst+1, :)
        corfluy(ivp, ifirst-1, :) = corfluy(ivp, ifirst+1, :)
    endif
#if defined MPI 
    if (float(ilast + ii * Lm) .EQ. IMAX_GRID) then
#else
    if (float(ilast) .EQ. IMAX_GRID) then
#endif
        corflux(ivp, ilast, :) = corflux(ivp, ilast-1, :)
        corfluy(ivp, ilast, :) = corfluy(ivp, ilast-1, :)
        corflux(ivp, ilast+1, :) = corflux(ivp, ilast-1, :)
        corfluy(ivp, ilast+1, :) = corfluy(ivp, ilast-1, :)
    endif

#if defined MPI 
    if (float(jfirst + jj * Mm) .EQ. JMIN_GRID) then
#else
    if (float(jfirst) .EQ. JMIN_GRID) then
#endif
        corflux(ivp, :, jfirst) = corflux(ivp, :, jfirst+1)
        corfluy(ivp, :, jfirst) = corfluy(ivp, :, jfirst+1)
        corflux(ivp, :, jfirst-1) = corflux(ivp, :, jfirst+1)
        corfluy(ivp, :, jfirst-1) = corfluy(ivp, :, jfirst+1)
    endif
#if defined MPI 
    if (float(jlast + jj * Mm) .EQ. JMAX_GRID) then
#else
    if (float(jlast) .EQ. JMAX_GRID) then
#endif
        corflux(ivp, :, jlast) = corflux(ivp, :, jlast-1)
        corfluy(ivp, :, jlast) = corfluy(ivp, :, jlast-1)
        corflux(ivp, :, jlast+1) = corflux(ivp, :, jlast-1)
        corfluy(ivp, :, jlast+1) = corfluy(ivp, :, jlast-1)
    endif

! corners
#if defined MPI 
    if ((float(ifirst + ii * Lm) .EQ. IMIN_GRID) .and.  &
        (float(jfirst + jj * Mm) .EQ. JMIN_GRID)) then
#else
    if ((float(ifirst) .EQ. IMIN_GRID) .and.  &
        (float(jfirst) .EQ. JMIN_GRID)) then
#endif
        corflux(ivp, ifirst, jfirst) = corflux(ivp, ifirst+1, jfirst+1)
        corfluy(ivp, ifirst, jfirst) = corfluy(ivp, ifirst+1, jfirst+1)
    endif
#if defined MPI 
    if ((float(ifirst + ii * Lm) .EQ. IMIN_GRID) .and.  &
        (float(jlast + jj * Mm) .EQ. JMAX_GRID)) then
#else
    if ((float(ifirst) .EQ. IMIN_GRID) .and.  &
        (float(jlast) .EQ. JMAX_GRID)) then
#endif
        corflux(ivp, ifirst, jlast) = corflux(ivp, ifirst+1, jlast-1)
        corfluy(ivp, ifirst, jlast) = corfluy(ivp, ifirst+1, jlast-1)
    endif
#if defined MPI 
    if ((float(ilast + ii * Lm) .EQ. IMAX_GRID) .and.  &
        (float(jlast + jj * Mm) .EQ. JMAX_GRID)) then
#else
    if ((float(ilast) .EQ. IMAX_GRID) .and.  &
        (float(jlast) .EQ. JMAX_GRID)) then
#endif
        corflux(ivp, ilast, jlast) = corflux(ivp, ilast-1, jlast-1)
        corfluy(ivp, ilast, jlast) = corfluy(ivp, ilast-1, jlast-1)
    endif
#if defined MPI 
    if ((float(ilast + ii * Lm) .EQ. IMAX_GRID) .and.  &
        (float(jfirst + jj * Mm) .EQ. JMIN_GRID)) then
#else
    if ((float(ilast) .EQ. IMAX_GRID) .and.  &
        (float(jfirst) .EQ. JMIN_GRID)) then
#endif
        corflux(ivp, ilast, jfirst) = corflux(ivp, ilast-1, jfirst+1)
        corfluy(ivp, ilast, jfirst) = corfluy(ivp, ilast-1, jfirst+1)
    endif

    enddo ! ivp

    END SUBROUTINE sed_obc_corflu
!!=============================================================================
   
    SUBROUTINE sed_meshedges_corflu(ifirst, ilast, jfirst, jlast)
 
        !&E------------------------------------------------------------------------
        !&E       ***  ROUTINE sed_meshedges_corflu ***
        !&E
        !&E ** Purpose :  interpolate corflux and corfluy on mesh edges (in u & v) 
        !&E
        !&E ** Description :  make interpolation from rho point to u,v points
        !&E
        !&E ** Called by : MUSTANG_update
        !&E--------------------------------------------------------------------------
        !! * Arguments
        INTEGER,INTENT(IN) :: ifirst, ilast, jfirst, jlast

        !! * Local declarations
        INTEGER :: i,j,iv, ivp
        REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY) :: tmpx
        REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY) :: tmpy

        !! * Executable part
        DO iv = isand1, isand2           
            tmpy(:, :) = corfluy(iv, :, :)
            tmpx(:, :) = corflux(iv, :, :)
            DO j = jfirst, jlast+1
                DO i = ifirst, ilast+1
                    corflux(iv, i, j) = 0.5_rsh * (tmpx(i-1, j) + tmpx(i, j))
                    corfluy(iv, i, j) = 0.5_rsh * (tmpy(i, j-1) + tmpy(i, j))
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE sed_meshedges_corflu
!!=============================================================================
#endif /* MUSTANG_CORFLUX */

 
! **TODO** code for CROCO if needed
! SUBROUTINE bathy_actu_fromfile(h0)
! subroutine sed_exchange_hxe_MARS(iwhat,xh0,xssh) 
! SUBROUTINE sed_exchange_maskbedload_MARS
! SUBROUTINE sed_outsaverestart(h0)

#endif

END MODULE sed_MUSTANG_CROCO
