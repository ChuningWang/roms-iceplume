#include "cppdefs.h"

MODULE sedchem

   !!======================================================================
   !!                        ***  Module sedchem  ***
   !! sediment :   Variable for chemistry of the CO2 cycle
   !!======================================================================
# if defined key_pisces
   !!   modules used
   USE sed     ! sediment global variable
#ifdef AGRIF
   USE par_sed
#endif
   USE sedarr
   USE sms_pisces, ONLY : neos, rtrn

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC sed_chem
   PUBLIC ahini_for_at_sed     !
   PUBLIC solve_at_general_sed !

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"

   ! Maximum number of iterations for each method
   INTEGER, PARAMETER :: jp_maxniter_atgen    = 20
   REAL(wp), PARAMETER :: pp_rdel_ah_target = 1.E-4

   !! * Module variables
   REAL(wp) :: &
      calcon = 1.03E-2        ! mean calcite concentration [Ca2+] in sea water [mole/kg solution] 

   REAL(wp) ::   rgas   = 83.14472      ! universal gas constants

   ! coeff. for density of sea water (Millero & Poisson 1981) 
   REAL(wp), DIMENSION(5)  :: Adsw                       
   DATA Adsw/8.24493E-1, -4.0899E-3, 7.6438E-5 , -8.246E-7, 5.3875E-9 /

   REAL(wp), DIMENSION(3)  :: Bdsw 
   DATA Bdsw / -5.72466E-3, 1.0227E-4, -1.6546E-6 /

   REAL(wp)  :: Cdsw = 4.8314E-4

   REAL(wp), DIMENSION(6)  :: Ddsw                    
   DATA Ddsw / 999.842594 , 6.793952E-2 , -9.095290E-3, 1.001685E-4, -1.120083E-6, 6.536332E-9/

  REAL(wp) :: devk10  = -25.5
   REAL(wp) :: devk11  = -15.82
   REAL(wp) :: devk12  = -29.48
   REAL(wp) :: devk13  = -20.02
   REAL(wp) :: devk14  = -18.03
   REAL(wp) :: devk15  = -9.78
   REAL(wp) :: devk16  = -48.76
   REAL(wp) :: devk17  = -14.51
   REAL(wp) :: devk18  = -23.12
   REAL(wp) :: devk19  = -26.57
   REAL(wp) :: devk110  = -29.48
   !
   REAL(wp) :: devk20  = 0.1271
   REAL(wp) :: devk21  = -0.0219
   REAL(wp) :: devk22  = 0.1622
   REAL(wp) :: devk23  = 0.1119
   REAL(wp) :: devk24  = 0.0466
   REAL(wp) :: devk25  = -0.0090
   REAL(wp) :: devk26  = 0.5304
   REAL(wp) :: devk27  = 0.1211
   REAL(wp) :: devk28  = 0.1758
   REAL(wp) :: devk29  = 0.2020
   REAL(wp) :: devk210  = 0.1622
   !
   REAL(wp) :: devk30  = 0.
   REAL(wp) :: devk31  = 0.
   REAL(wp) :: devk32  = 2.608E-3
   REAL(wp) :: devk33  = -1.409e-3
   REAL(wp) :: devk34  = 0.316e-3
   REAL(wp) :: devk35  = -0.942e-3
   REAL(wp) :: devk36  = 0.
   REAL(wp) :: devk37  = -0.321e-3
   REAL(wp) :: devk38  = -2.647e-3
   REAL(wp) :: devk39  = -3.042e-3
   REAL(wp) :: devk310  = -2.6080e-3
   !
   REAL(wp) :: devk40  = -3.08E-3
   REAL(wp) :: devk41  = 1.13E-3
   REAL(wp) :: devk42  = -2.84E-3
   REAL(wp) :: devk43  = -5.13E-3
   REAL(wp) :: devk44  = -4.53e-3
   REAL(wp) :: devk45  = -3.91e-3
   REAL(wp) :: devk46  = -11.76e-3
   REAL(wp) :: devk47  = -2.67e-3
   REAL(wp) :: devk48  = -5.15e-3
   REAL(wp) :: devk49  = -4.08e-3
   REAL(wp) :: devk410  = -2.84e-3
   !
   REAL(wp) :: devk50  = 0.0877E-3
   REAL(wp) :: devk51  = -0.1475E-3
   REAL(wp) :: devk52  = 0.
   REAL(wp) :: devk53  = 0.0794E-3
   REAL(wp) :: devk54  = 0.09e-3
   REAL(wp) :: devk55  = 0.054e-3
   REAL(wp) :: devk56  = 0.3692E-3
   REAL(wp) :: devk57  = 0.0427e-3
   REAL(wp) :: devk58  = 0.09e-3
   REAL(wp) :: devk59  = 0.0714e-3
   REAL(wp) :: devk510  = 0.0

   !! $Id: sedchem.F90 10356 2018-11-22 11:39:19Z cetlod $
CONTAINS

   SUBROUTINE sed_chem( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_chem  ***
      !!
      !! ** Purpose :   set chemical constants
      !!
      !!   History :
      !!        !  04-10  (N. Emprin, M. Gehlen )  Original code
      !!        !  06-04  (C. Ethe)  Re-organization
      !!----------------------------------------------------------------------
      !!* Arguments
      INTEGER, INTENT(in) :: kt                     ! time step

      INTEGER  :: ji, jj
      REAL(wp) :: ztc, ztc2
      REAL(wp) :: zsal, zsal15  
      REAL(wp) :: zdens0, zaw, zbw, zcw    
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY,15) ::   zchem_data
      !!----------------------------------------------------------------------


      IF (lwp) WRITE(numsed,*) ' Getting Chemical constants from tracer model at time kt = ', kt
      IF (lwp) WRITE(numsed,*) ' '

      ! reading variables
      zchem_data(:,:,:) = rtrn

      IF (ln_sediment_offline) THEN
         CALL sed_chem_cst
      ELSE
         DO jj = JRANGE
            DO ji = IRANGE
               IF ( tmask(ji,jj,ikt) == 1 ) THEN
                  zchem_data(ji,jj,1) = ak13  (ji,jj,ikt)
                  zchem_data(ji,jj,2) = ak23  (ji,jj,ikt)
                  zchem_data(ji,jj,3) = akb3  (ji,jj,ikt)
                  zchem_data(ji,jj,4) = akw3  (ji,jj,ikt)
                  zchem_data(ji,jj,5) = aksp  (ji,jj,ikt)
                  zchem_data(ji,jj,6) = borat (ji,jj,ikt)
                  zchem_data(ji,jj,7) = ak1p3 (ji,jj,ikt)
                  zchem_data(ji,jj,8) = ak2p3 (ji,jj,ikt)
                  zchem_data(ji,jj,9) = ak3p3 (ji,jj,ikt)
                  zchem_data(ji,jj,10)= aksi3 (ji,jj,ikt)
                  zchem_data(ji,jj,11)= sio3eq(ji,jj,ikt)
                  zchem_data(ji,jj,12)= aks3  (ji,jj,ikt)
                  zchem_data(ji,jj,13)= akf3  (ji,jj,ikt)
                  zchem_data(ji,jj,14)= sulfat(ji,jj,ikt)
                  zchem_data(ji,jj,15)= fluorid(ji,jj,ikt)
               ENDIF
            ENDDO
         ENDDO

         CALL pack_arr ( jpoce, ak1s  (1:jpoce), zchem_data(PRIV_2D_BIOARRAY,1) , iarroce(1:jpoce) )
         CALL pack_arr ( jpoce, ak2s  (1:jpoce), zchem_data(PRIV_2D_BIOARRAY,2) , iarroce(1:jpoce) )
         CALL pack_arr ( jpoce, akbs  (1:jpoce), zchem_data(PRIV_2D_BIOARRAY,3) , iarroce(1:jpoce) )
         CALL pack_arr ( jpoce, akws  (1:jpoce), zchem_data(PRIV_2D_BIOARRAY,4) , iarroce(1:jpoce) )
         CALL pack_arr ( jpoce, aksps (1:jpoce), zchem_data(PRIV_2D_BIOARRAY,5) , iarroce(1:jpoce) )
         CALL pack_arr ( jpoce, borats(1:jpoce), zchem_data(PRIV_2D_BIOARRAY,6) , iarroce(1:jpoce) )
         CALL pack_arr ( jpoce, ak1ps (1:jpoce), zchem_data(PRIV_2D_BIOARRAY,7) , iarroce(1:jpoce) )
         CALL pack_arr ( jpoce, ak2ps (1:jpoce), zchem_data(PRIV_2D_BIOARRAY,8) , iarroce(1:jpoce) )
         CALL pack_arr ( jpoce, ak3ps (1:jpoce), zchem_data(PRIV_2D_BIOARRAY,9) , iarroce(1:jpoce) )
         CALL pack_arr ( jpoce, aksis (1:jpoce), zchem_data(PRIV_2D_BIOARRAY,10), iarroce(1:jpoce) )
         CALL pack_arr ( jpoce, sieqs (1:jpoce), zchem_data(PRIV_2D_BIOARRAY,11), iarroce(1:jpoce) )
         CALL pack_arr ( jpoce, aks3s (1:jpoce), zchem_data(PRIV_2D_BIOARRAY,12), iarroce(1:jpoce) )
         CALL pack_arr ( jpoce, akf3s (1:jpoce), zchem_data(PRIV_2D_BIOARRAY,13), iarroce(1:jpoce) )
         CALL pack_arr ( jpoce, sulfats(1:jpoce), zchem_data(PRIV_2D_BIOARRAY,14), iarroce(1:jpoce) )
         CALL pack_arr ( jpoce, fluorids(1:jpoce), zchem_data(PRIV_2D_BIOARRAY,15), iarroce(1:jpoce) )
      ENDIF

      DO ji = 1, jpoce
         ztc     = temp(ji)
         ztc2    = ztc * ztc
         ! zqtt    = ztkel * 0.01
         zsal    = salt(ji)
         zsal15  = SQRT( zsal ) * zsal

         ! Density of Sea Water - F(temp,sal) [kg/m3]
         zdens0 =  Ddsw(1) + Ddsw(2) * ztc + Ddsw(3) * ztc2 &
                  + Ddsw(4) * ztc * ztc2 + Ddsw(5) * ztc2 * ztc2 &
                  + Ddsw(6) * ztc * ztc2 * ztc2
         zaw =  Adsw(1) + Adsw(2) * ztc + Adsw(3)* ztc2 + Adsw(4) * ztc * ztc2 &
              + Adsw(5) * ztc2 * ztc2
         zbw =  Bdsw(1) + Bdsw(2) * ztc + Bdsw(3) * ztc2
         zcw =  Cdsw
         densSW(ji) = zdens0 + zaw * zsal + zbw * zsal15 + zcw * zsal * zsal
         densSW(ji) = densSW(ji) * 1E-3   ! to get dens in [kg/l]

         ak12s  (ji) = ak1s (ji) * ak2s (ji)
         ak12ps (ji) = ak1ps(ji) * ak2ps(ji)
         ak123ps(ji) = ak1ps(ji) * ak2ps(ji) * ak3ps(ji)

         calcon2(ji) = 0.01028 * ( salt(ji) / 35. ) * densSW(ji)
      ENDDO
       
   END SUBROUTINE sed_chem

   SUBROUTINE ahini_for_at_sed(p_hini)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE ahini_for_at  ***
      !!
      !! Subroutine returns the root for the 2nd order approximation of the
      !! DIC -- B_T -- A_CB equation for [H+] (reformulated as a cubic 
      !! polynomial) around the local minimum, if it exists.
      !! Returns * 1E-03_wp if p_alkcb <= 0
      !!         * 1E-10_wp if p_alkcb >= 2*p_dictot + p_bortot
      !!         * 1E-07_wp if 0 < p_alkcb < 2*p_dictot + p_bortot
      !!                    and the 2nd order approximation does not have 
      !!                    a solution
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(jpoce,jpksed), INTENT(OUT)  ::  p_hini
      INTEGER  ::   ji, jk
      REAL(wp)  ::  zca1, zba1
      REAL(wp)  ::  zd, zsqrtd, zhmin
      REAL(wp)  ::  za2, za1, za0
      REAL(wp)  ::  p_dictot, p_bortot, p_alkcb

      !
      DO jk = 1, jpksed
         DO ji = 1, jpoce
            p_alkcb  = pwcp(ji,jk,jwalk) / densSW(ji)
            p_dictot = pwcp(ji,jk,jwdic) / densSW(ji) 
            p_bortot = borats(ji) / densSW(ji)
            IF (p_alkcb <= 0.) THEN
                p_hini(ji,jk) = 1.e-3
            ELSEIF (p_alkcb >= (2.*p_dictot + p_bortot)) THEN
                p_hini(ji,jk) = 1.e-10
            ELSE
                zca1 = p_dictot/( p_alkcb + rtrn )
                zba1 = p_bortot/ (p_alkcb + rtrn )
           ! Coefficients of the cubic polynomial
                za2 = aKbs(ji)*(1. - zba1) + ak1s(ji)*(1.-zca1)
                za1 = ak1s(ji)*akbs(ji)*(1. - zba1 - zca1)    &
                &     + ak1s(ji)*ak2s(ji)*(1. - (zca1+zca1))
                za0 = ak1s(ji)*ak2s(ji)*akbs(ji)*(1. - zba1 - (zca1+zca1))
                                        ! Taylor expansion around the minimum
                zd = za2*za2 - 3.*za1   ! Discriminant of the quadratic equation
                                        ! for the minimum close to the root

                IF(zd > 0.) THEN        ! If the discriminant is positive
                  zsqrtd = SQRT(zd)
                  IF(za2 < 0) THEN
                    zhmin = (-za2 + zsqrtd)/3.
                  ELSE
                    zhmin = -za1/(za2 + zsqrtd)
                  ENDIF
                  p_hini(ji,jk) = zhmin + SQRT(-(za0 + zhmin*(za1 + zhmin*(za2 + zhmin)))/zsqrtd)
                ELSE
                  p_hini(ji,jk) = 1.e-7
                ENDIF
             !
            ENDIF
         END DO
      END DO
      !
   END SUBROUTINE ahini_for_at_sed

   !===============================================================================
   SUBROUTINE anw_infsup_sed( p_alknw_inf, p_alknw_sup )

   ! Subroutine returns the lower and upper bounds of "non-water-selfionization"
   ! contributions to total alkalinity (the infimum and the supremum), i.e
   ! inf(TA - [OH-] + [H+]) and sup(TA - [OH-] + [H+])

   ! Argument variables
   INTEGER :: jk
   REAL(wp), DIMENSION(jpoce,jpksed), INTENT(OUT) :: p_alknw_inf
   REAL(wp), DIMENSION(jpoce,jpksed), INTENT(OUT) :: p_alknw_sup

   DO jk = 1, jpksed
      p_alknw_inf(:,jk) =  -pwcp(:,jk,jwpo4) / densSW(:)
      p_alknw_sup(:,jk) =   (2. * pwcp(:,jk,jwdic) + 2. * pwcp(:,jk,jwpo4) + pwcp(:,jk,jwsil)     &
   &                        + borats(:) ) / densSW(:)
   END DO

   END SUBROUTINE anw_infsup_sed


   SUBROUTINE solve_at_general_sed( p_hini, zhi )

   ! Universal pH solver that converges from any given initial value,
   ! determines upper an lower bounds for the solution if required

   ! Argument variables
   !--------------------
   REAL(wp), DIMENSION(jpoce,jpksed), INTENT(IN)   :: p_hini
   REAL(wp), DIMENSION(jpoce,jpksed), INTENT(OUT)  :: zhi

   ! Local variables
   !-----------------
   INTEGER   ::  ji, jk, jn
   REAL(wp)  ::  zh_ini, zh, zh_prev, zh_lnfactor
   REAL(wp)  ::  zdelta, zh_delta
   REAL(wp)  ::  zeqn, zdeqndh, zalka
   REAL(wp)  ::  aphscale
   REAL(wp)  ::  znumer_dic, zdnumer_dic, zdenom_dic, zalk_dic, zdalk_dic
   REAL(wp)  ::  znumer_bor, zdnumer_bor, zdenom_bor, zalk_bor, zdalk_bor
   REAL(wp)  ::  znumer_po4, zdnumer_po4, zdenom_po4, zalk_po4, zdalk_po4
   REAL(wp)  ::  znumer_sil, zdnumer_sil, zdenom_sil, zalk_sil, zdalk_sil
   REAL(wp)  ::  znumer_so4, zdnumer_so4, zdenom_so4, zalk_so4, zdalk_so4
   REAL(wp)  ::  znumer_flu, zdnumer_flu, zdenom_flu, zalk_flu, zdalk_flu
   REAL(wp)  ::  zalk_wat, zdalk_wat
   REAL(wp)  ::  zfact, p_alktot, zdic, zbot, zpt, zst, zft, zsit
   LOGICAL   ::  l_exitnow
   REAL(wp), PARAMETER :: pz_exp_threshold = 1.0
   REAL(wp), DIMENSION(jpoce,jpksed) :: zalknw_inf, zalknw_sup, rmask, zh_min, zh_max, zeqn_absmin

      !  Allocate temporary workspace
   CALL anw_infsup_sed( zalknw_inf, zalknw_sup )

   rmask(:,:) = 1.0
   zhi(:,:)   = 0.

   ! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
   DO jk = 1, jpksed
      DO ji = 1, jpoce
         IF (rmask(ji,jk) == 1.) THEN
            p_alktot = pwcp(ji,jk,jwalk) / densSW(ji)
            aphscale = 1. + sulfats(ji)/aks3s(ji)
            zh_ini = p_hini(ji,jk)

            zdelta = (p_alktot-zalknw_inf(ji,jk))**2 + 4.*akws(ji) / aphscale

            IF(p_alktot >= zalknw_inf(ji,jk)) THEN
               zh_min(ji,jk) = 2.*akws(ji) /( p_alktot-zalknw_inf(ji,jk) + SQRT(zdelta) )
            ELSE
               zh_min(ji,jk) = aphscale * (-(p_alktot-zalknw_inf(ji,jk)) + SQRT(zdelta) ) / 2.
            ENDIF

            zdelta = (p_alktot-zalknw_sup(ji,jk))**2 + 4.*akws(ji) / aphscale

            IF(p_alktot <= zalknw_sup(ji,jk)) THEN
               zh_max(ji,jk) = aphscale * (-(p_alktot-zalknw_sup(ji,jk)) + SQRT(zdelta) ) / 2.
            ELSE
               zh_max(ji,jk) = 2.*akws(ji) /( p_alktot-zalknw_sup(ji,jk) + SQRT(zdelta) )
            ENDIF

            zhi(ji,jk) = MAX(MIN(zh_max(ji,jk), zh_ini), zh_min(ji,jk))
         ENDIF
      END DO
   END DO

   zeqn_absmin(:,:) = HUGE(1.)

   DO jn = 1, jp_maxniter_atgen
   DO jk = 1, jpksed
      DO ji = 1, jpoce
         IF (rmask(ji,jk) == 1.) THEN

            p_alktot = pwcp(ji,jk,jwalk) / densSW(ji)
            zdic  = pwcp(ji,jk,jwdic) / densSW(ji)
            zbot  = borats(ji) / densSW(ji)
            zpt =  pwcp(ji,jk,jwpo4) / densSW(ji)
            zsit = pwcp(ji,jk,jwsil) / densSW(ji)
            zst = sulfats(ji)
            zft = fluorids(ji)
            aphscale = 1. + sulfats(ji)/aks3s(ji)
            zh = zhi(ji,jk)
            zh_prev = zh

            ! H2CO3 - HCO3 - CO3 : n=2, m=0
            znumer_dic = 2.*ak1s(ji)*ak2s(ji) + zh*ak1s(ji)
            zdenom_dic = ak1s(ji)*ak2s(ji) + zh*(ak1s(ji) + zh)
            zalk_dic   = zdic * (znumer_dic/zdenom_dic)
            zdnumer_dic = ak1s(ji)*ak1s(ji)*ak2s(ji) + zh     &
                          *(4.*ak1s(ji)*ak2s(ji) + zh*ak1s(ji))
            zdalk_dic   = -zdic*(zdnumer_dic/zdenom_dic**2)


            ! B(OH)3 - B(OH)4 : n=1, m=0
            znumer_bor = akbs(ji)
            zdenom_bor = akbs(ji) + zh
            zalk_bor   = zbot * (znumer_bor/zdenom_bor)
            zdnumer_bor = akbs(ji)
            zdalk_bor   = -zbot*(zdnumer_bor/zdenom_bor**2)


            ! H3PO4 - H2PO4 - HPO4 - PO4 : n=3, m=1
            znumer_po4 = 3.*ak1ps(ji)*ak2ps(ji)*ak3ps(ji)  &
            &            + zh*(2.*ak1ps(ji)*ak2ps(ji) + zh* ak1ps(ji))
            zdenom_po4 = ak1ps(ji)*ak2ps(ji)*ak3ps(ji)     &
            &            + zh*( ak1ps(ji)*ak2ps(ji) + zh*(ak1ps(ji) + zh))
            zalk_po4   = zpt * (znumer_po4/zdenom_po4 - 1.) ! Zero level of H3PO4 = 1
            zdnumer_po4 = ak1ps(ji)*ak2ps(ji)*ak1ps(ji)*ak2ps(ji)*ak3ps(ji)  &
            &             + zh*(4.*ak1ps(ji)*ak1ps(ji)*ak2ps(ji)*ak3ps(ji)         &
            &             + zh*(9.*ak1ps(ji)*ak2ps(ji)*ak3ps(ji)                         &
            &             + ak1ps(ji)*ak1ps(ji)*ak2ps(ji)                                &
            &             + zh*(4.*ak1ps(ji)*ak2ps(ji) + zh * ak1ps(ji) ) ) )
            zdalk_po4   = -zpt * (zdnumer_po4/zdenom_po4**2)

            ! H4SiO4 - H3SiO4 : n=1, m=0
            znumer_sil = aksis(ji)
            zdenom_sil = aksis(ji) + zh
            zalk_sil   = zsit * (znumer_sil/zdenom_sil)
            zdnumer_sil = aksis(ji)
            zdalk_sil   = -zsit * (zdnumer_sil/zdenom_sil**2)

            ! HSO4 - SO4 : n=1, m=1
            aphscale = 1.0 + zst/aks3s(ji)
            znumer_so4 = aks3s(ji) * aphscale
            zdenom_so4 = aks3s(ji) * aphscale + zh
            zalk_so4   = zst * (znumer_so4/zdenom_so4 - 1.)
            zdnumer_so4 = aks3s(ji) * aphscale
            zdalk_so4   = -zst * (zdnumer_so4/zdenom_so4**2)

            ! HF - F : n=1, m=1
            znumer_flu =  akf3s(ji)
            zdenom_flu =  akf3s(ji) + zh
            zalk_flu   =  zft * (znumer_flu/zdenom_flu - 1.)
            zdnumer_flu = akf3s(ji)
            zdalk_flu   = -zft * (zdnumer_flu/zdenom_flu**2)

            ! H2O - OH
            zalk_wat   = akws(ji)/zh - zh/aphscale
            zdalk_wat  = -akws(ji)/zh**2 - 1./aphscale

            ! CALCULATE [ALK]([CO3--], [HCO3-])
            zeqn = zalk_dic + zalk_bor + zalk_po4 + zalk_sil   &
            &      + zalk_so4 + zalk_flu                       &
            &      + zalk_wat - p_alktot

            zalka = p_alktot - (zalk_bor + zalk_po4 + zalk_sil   &
            &       + zalk_so4 + zalk_flu + zalk_wat)

            zdeqndh = zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil &
            &         + zdalk_so4 + zdalk_flu + zdalk_wat

            ! Adapt bracketing interval
            IF(zeqn > 0.) THEN
               zh_min(ji,jk) = zh_prev
            ELSEIF(zeqn < 0.) THEN
               zh_max(ji,jk) = zh_prev
            ENDIF

            IF(ABS(zeqn) >= 0.5*zeqn_absmin(ji,jk)) THEN
            ! if the function evaluation at the current point is
            ! not decreasing faster than with a bisection step (at least linearly)
            ! in absolute value take one bisection step on [ph_min, ph_max]
            ! ph_new = (ph_min + ph_max)/2d0
            !
            ! In terms of [H]_new:
            ! [H]_new = 10**(-ph_new)
            !         = 10**(-(ph_min + ph_max)/2d0)
            !         = SQRT(10**(-(ph_min + phmax)))
            !         = SQRT(zh_max * zh_min)
               zh = SQRT(zh_max(ji,jk) * zh_min(ji,jk))
               zh_lnfactor = (zh - zh_prev)/zh_prev ! Required to test convergence below
            ELSE
            ! dzeqn/dpH = dzeqn/d[H] * d[H]/dpH
            !           = -zdeqndh * LOG(10) * [H]
            ! \Delta pH = -zeqn/(zdeqndh*d[H]/dpH) = zeqn/(zdeqndh*[H]*LOG(10))
            !
            ! pH_new = pH_old + \deltapH
            !
            ! [H]_new = 10**(-pH_new)
            !         = 10**(-pH_old - \Delta pH)
            !         = [H]_old * 10**(-zeqn/(zdeqndh*[H]_old*LOG(10)))
            !         = [H]_old * EXP(-LOG(10)*zeqn/(zdeqndh*[H]_old*LOG(10)))
            !         = [H]_old * EXP(-zeqn/(zdeqndh*[H]_old))

               zh_lnfactor = -zeqn/(zdeqndh*zh_prev)

               IF(ABS(zh_lnfactor) > pz_exp_threshold) THEN
                  zh          = zh_prev*EXP(zh_lnfactor)
               ELSE
                  zh_delta    = zh_lnfactor*zh_prev
                  zh          = zh_prev + zh_delta
               ENDIF

               IF( zh < zh_min(ji,jk) ) THEN
               ! if [H]_new < [H]_min
               ! i.e., if ph_new > ph_max then
               ! take one bisection step on [ph_prev, ph_max]
               ! ph_new = (ph_prev + ph_max)/2d0
               ! In terms of [H]_new:
               ! [H]_new = 10**(-ph_new)
               !         = 10**(-(ph_prev + ph_max)/2d0)
               !         = SQRT(10**(-(ph_prev + phmax)))
               !         = SQRT([H]_old*10**(-ph_max))
               !         = SQRT([H]_old * zh_min)
                  zh                = SQRT(zh_prev * zh_min(ji,jk))
                  zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below
               ENDIF

               IF( zh > zh_max(ji,jk) ) THEN
               ! if [H]_new > [H]_max
               ! i.e., if ph_new < ph_min, then
               ! take one bisection step on [ph_min, ph_prev]
               ! ph_new = (ph_prev + ph_min)/2d0
               ! In terms of [H]_new:
               ! [H]_new = 10**(-ph_new)
               !         = 10**(-(ph_prev + ph_min)/2d0)
               !         = SQRT(10**(-(ph_prev + ph_min)))
               !         = SQRT([H]_old*10**(-ph_min))
               !         = SQRT([H]_old * zhmax)
                  zh                = SQRT(zh_prev * zh_max(ji,jk))
                  zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below
               ENDIF
            ENDIF

            zeqn_absmin(ji,jk) = MIN( ABS(zeqn), zeqn_absmin(ji,jk))

            ! Stop iterations once |\delta{[H]}/[H]| < rdel
            ! <=> |(zh - zh_prev)/zh_prev| = |EXP(-zeqn/(zdeqndh*zh_prev)) -1| < rdel
            ! |EXP(-zeqn/(zdeqndh*zh_prev)) -1| ~ |zeqn/(zdeqndh*zh_prev)|
            ! Alternatively:
            ! |\Delta pH| = |zeqn/(zdeqndh*zh_prev*LOG(10))|
            !             ~ 1/LOG(10) * |\Delta [H]|/[H]
            !             < 1/LOG(10) * rdel

            ! Hence |zeqn/(zdeqndh*zh)| < rdel

            ! rdel <-- pp_rdel_ah_target
            l_exitnow = (ABS(zh_lnfactor) < pp_rdel_ah_target)

            IF(l_exitnow) THEN
               rmask(ji,jk) = 0.
            ENDIF

            zhi(ji,jk) =  zh

            IF(jn >= jp_maxniter_atgen) THEN
               zhi(ji,jk) = -1.
            ENDIF

         ENDIF
      END DO
   END DO
   END DO
   !

   END SUBROUTINE solve_at_general_sed

   SUBROUTINE sed_chem_cst
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sed_chem_cst  ***
      !!
      !! ** Purpose :   Sea water chemistry computed following MOCSY protocol
      !!                Computation is done at the bottom of the ocean only
      !!
      !! ** Method  : - ...
      !!---------------------------------------------------------------------
      INTEGER  ::   ji
      REAL(wp), DIMENSION(jpoce) :: saltprac, temps
      REAL(wp) ::   ztkel, ztkel1, zt , zsal  , zsal2 , zbuf1 , zbuf2
      REAL(wp) ::   ztgg , ztgg2, ztgg3 , ztgg4 , ztgg5
      REAL(wp) ::   zpres, ztc  , zcl   , zcpexp, zoxy  , zcpexp2
      REAL(wp) ::   zsqrt, ztr  , zlogt , zcek1, zc1, zplat
      REAL(wp) ::   zis  , zis2 , zsal15, zisqrt, za1, za2
      REAL(wp) ::   zckb , zck1 , zck2  , zckw  , zak1 , zak2  , zakb , zaksp0, zakw
      REAL(wp) ::   zck1p, zck2p, zck3p, zcksi, zak1p, zak2p, zak3p, zaksi
      REAL(wp) ::   zst  , zft  , zcks  , zckf  , zaksp1
      REAL(wp) ::   total2free, free2SWS, total2SWS, SWS2total
      !!---------------------------------------------------------------------
      !
      ! Computation of chemical constants require practical salinity
      ! Thus, when TEOS08 is used, absolute salinity is converted to
      ! practical salinity
      ! -------------------------------------------------------------
      IF (neos == -1) THEN
         saltprac(:) = salt(:) * 35.0 / 35.16504
      ELSE
         saltprac(:) = temp(:)
      ENDIF

      !
      ! Computations of chemical constants require in situ temperature
      ! Here a quite simple formulation is used to convert
      ! potential temperature to in situ temperature. The errors is less than
      ! 0.04°C relative to an exact computation
      ! ---------------------------------------------------------------------
         DO ji = 1, jpoce
            zpres = zkbot(ji) / 1000.
            za1 = 0.04 * ( 1.0 + 0.185 * temp(ji) + 0.035 * (saltprac(ji) - 35.0) )
            za2 = 0.0075 * ( 1.0 - temp(ji) / 30.0 )
            temps(ji) = temp(ji) - za1 * zpres + za2 * zpres**2
         END DO

      ! CHEMICAL CONSTANTS - DEEP OCEAN
      ! -------------------------------
      DO ji = 1, jpoce
         ! SET PRESSION ACCORDING TO SAUNDER (1980)
         zc1 = 5.92E-3 
         zpres = ((1-zc1)-SQRT(((1-zc1)**2)-(8.84E-6*zkbot(ji)))) / 4.42E-6
         zpres = zpres / 10.0

         ! SET ABSOLUTE TEMPERATURE
         ztkel   = temps(ji) + 273.15
         zsal    = saltprac(ji)
         zsqrt  = SQRT( zsal )
         zsal15  = zsqrt * zsal
         zlogt  = LOG( ztkel )
         ztr    = 1. / ztkel
         zis    = 19.924 * zsal / ( 1000.- 1.005 * zsal )
         zis2   = zis * zis
         zisqrt = SQRT( zis )
         ztc    = temps(ji)

         ! CHLORINITY (WOOSTER ET AL., 1969)
         zcl     = zsal / 1.80655

         ! TOTAL SULFATE CONCENTR. [MOLES/kg soln]
         zst     = 0.14 * zcl /96.062

         ! TOTAL FLUORIDE CONCENTR. [MOLES/kg soln]
         zft     = 0.000067 * zcl /18.9984

         ! DISSOCIATION CONSTANT FOR SULFATES on free H scale (Dickson 1990)
         zcks    = EXP(-4276.1 * ztr + 141.328 - 23.093 * zlogt         &
         &         + (-13856. * ztr + 324.57 - 47.986 * zlogt) * zisqrt &
         &         + (35474. * ztr - 771.54 + 114.723 * zlogt) * zis    &
         &         - 2698. * ztr * zis**1.5 + 1776.* ztr * zis2         &
         &         + LOG(1.0 - 0.001005 * zsal))

         ! DISSOCIATION CONSTANT FOR FLUORIDES on free H scale (Dickson and Riley 79)
         zckf    = EXP( 1590.2*ztr - 12.641 + 1.525*zisqrt   &
         &         + LOG(1.0d0 - 0.001005d0*zsal)            &
         &         + LOG(1.0d0 + zst/zcks))

         ! DISSOCIATION CONSTANT FOR CARBONATE AND BORATE
         zckb=  (-8966.90 - 2890.53*zsqrt - 77.942*zsal        &
         &      + 1.728*zsal15 - 0.0996*zsal*zsal)*ztr         &
         &      + (148.0248 + 137.1942*zsqrt + 1.62142*zsal)   &
         &      + (-24.4344 - 25.085*zsqrt - 0.2474*zsal)      &
         &      * zlogt + 0.053105*zsqrt*ztkel

         ! DISSOCIATION COEFFICIENT FOR CARBONATE ACCORDING TO
         ! MEHRBACH (1973) REFIT BY MILLERO (1995), seawater scale
         zck1    = -1.0*(3633.86*ztr - 61.2172 + 9.6777*zlogt  &
                   - 0.011555*zsal + 0.0001152*zsal*zsal)
         zck2    = -1.0*(471.78*ztr + 25.9290 - 3.16967*zlogt      &
                   - 0.01781*zsal + 0.0001122*zsal*zsal)

         ! PKW (H2O) (MILLERO, 1995) from composite data
         zckw    = -13847.26 * ztr + 148.9652 - 23.6521 * zlogt + ( 118.67 * ztr    &
                   - 5.977 + 1.0495 * zlogt ) * zsqrt - 0.01615 * zsal

         ! CONSTANTS FOR PHOSPHATE (MILLERO, 1995)
         zck1p    = -4576.752*ztr + 115.540 - 18.453*zlogt   &
         &          + (-106.736*ztr + 0.69171) * zsqrt       &
         &          + (-0.65643*ztr - 0.01844) * zsal

         zck2p    = -8814.715*ztr + 172.1033 - 27.927*zlogt  &
         &          + (-160.340*ztr + 1.3566)*zsqrt          &
         &          + (0.37335*ztr - 0.05778)*zsal

         zck3p    = -3070.75*ztr - 18.126                    &
         &          + (17.27039*ztr + 2.81197) * zsqrt       &
         &          + (-44.99486*ztr - 0.09984) * zsal

         ! CONSTANT FOR SILICATE, MILLERO (1995)
         zcksi    = -8904.2*ztr  + 117.400 - 19.334*zlogt   &
         &          + (-458.79*ztr + 3.5913) * zisqrt       &
         &          + (188.74*ztr - 1.5998) * zis           &
         &          + (-12.1652*ztr + 0.07871) * zis2       &
         &          + LOG(1.0 - 0.001005*zsal)

         ! APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE IN SEAWATER
         !       (S=27-43, T=2-25 DEG C) at pres =0 (atmos. pressure) (MUCCI 1983)
         zaksp0  = -171.9065 -0.077993*ztkel + 2839.319*ztr + 71.595*LOG10( ztkel )   &
         &         + (-0.77712 + 0.00284263*ztkel + 178.34*ztr) * zsqrt  &
         &         - 0.07711*zsal + 0.0041249*zsal15

         ! CONVERT FROM DIFFERENT PH SCALES
         total2free  = 1.0/(1.0 + zst/zcks)
         free2SWS    = 1. + zst/zcks + zft/(zckf*total2free)
         total2SWS   = total2free * free2SWS
         SWS2total   = 1.0 / total2SWS


         ! K1, K2 OF CARBONIC ACID, KB OF BORIC ACID, KW (H2O) (LIT.?)
         zak1    = 10**(zck1) * total2SWS
         zak2    = 10**(zck2) * total2SWS
         zakb    = EXP( zckb ) * total2SWS
         zakw    = EXP( zckw )
         zaksp1  = 10**(zaksp0)
         zak1p   = exp( zck1p )
         zak2p   = exp( zck2p )
         zak3p   = exp( zck3p )
         zaksi   = exp( zcksi )
         zckf    = zckf * total2SWS

         ! FORMULA FOR CPEXP AFTER EDMOND & GIESKES (1970)
         !        (REFERENCE TO CULBERSON & PYTKOQICZ (1968) AS MADE
         !        IN BROECKER ET AL. (1982) IS INCORRECT; HERE RGAS IS
         !        TAKEN TENFOLD TO CORRECT FOR THE NOTATION OF pres  IN
         !        DBAR INSTEAD OF BAR AND THE EXPRESSION FOR CPEXP IS
         !        MULTIPLIED BY LN(10.) TO ALLOW USE OF EXP-FUNCTION
         !        WITH BASIS E IN THE FORMULA FOR AKSPP (CF. EDMOND
         !        & GIESKES (1970), P. 1285-1286 (THE SMALL
         !        FORMULA ON P. 1286 IS RIGHT AND CONSISTENT WITH THE
         !        SIGN IN PARTIAL MOLAR VOLUME CHANGE AS SHOWN ON P. 1285))
         zcpexp  = zpres / (rgas*ztkel)
         zcpexp2 = zpres * zcpexp

         ! KB OF BORIC ACID, K1,K2 OF CARBONIC ACID PRESSURE
         !        CORRECTION AFTER CULBERSON AND PYTKOWICZ (1968)
         !        (CF. BROECKER ET AL., 1982)

         zbuf1  = -     ( devk10 + devk20 * ztc + devk30 * ztc * ztc )
         zbuf2  = 0.5 * ( devk40 + devk50 * ztc )
         ak1s(ji) = zak1 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk11 + devk21 * ztc + devk31 * ztc * ztc )
         zbuf2  = 0.5 * ( devk41 + devk51 * ztc )
         ak2s(ji) = zak2 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk12 + devk22 * ztc + devk32 * ztc * ztc )
         zbuf2  = 0.5 * ( devk42 + devk52 * ztc )
         akbs(ji) = zakb * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk13 + devk23 * ztc + devk33 * ztc * ztc )
         zbuf2  = 0.5 * ( devk43 + devk53 * ztc )
         akws(ji) = zakw * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk14 + devk24 * ztc + devk34 * ztc * ztc )
         zbuf2  = 0.5 * ( devk44 + devk54 * ztc )
         aks3s(ji) = zcks * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk15 + devk25 * ztc + devk35 * ztc * ztc )
         zbuf2  = 0.5 * ( devk45 + devk55 * ztc )
         akf3s(ji) = zckf * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk17 + devk27 * ztc + devk37 * ztc * ztc )
         zbuf2  = 0.5 * ( devk47 + devk57 * ztc )
         ak1ps(ji) = zak1p * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk18 + devk28 * ztc + devk38 * ztc * ztc )
         zbuf2  = 0.5 * ( devk48 + devk58 * ztc )
         ak2ps(ji) = zak2p * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         zbuf1  =     - ( devk110 + devk210 * ztc + devk310 * ztc * ztc )
         zbuf2  = 0.5 * ( devk410 + devk510 * ztc )
         aksis(ji) = zaksi * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         ! Convert to total scale
         ak1s(ji)  = ak1s(ji)  * SWS2total
         ak2s(ji)  = ak2s(ji)  * SWS2total
         akbs(ji)  = akbs(ji)  * SWS2total
         akws(ji)  = akws(ji)  * SWS2total
         ak1ps(ji) = ak1ps(ji) * SWS2total
         ak2ps(ji) = ak2ps(ji) * SWS2total
         ak3ps(ji) = ak3ps(ji) * SWS2total
         aksis(ji) = aksis(ji) * SWS2total
         akf3s(ji) = akf3s(ji) / total2free

         ! APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE
         !        AS FUNCTION OF PRESSURE FOLLOWING MILLERO
         !        (P. 1285) AND BERNER (1976)
         zbuf1  =     - ( devk16 + devk26 * ztc + devk36 * ztc * ztc )
         zbuf2  = 0.5 * ( devk46 + devk56 * ztc )
         aksps(ji) = zaksp1 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

         ! TOTAL F, S, and BORATE CONCENTR. [MOLES/L]
         borats(ji)   = 0.0002414 * zcl / 10.811
         sulfats(ji)  = zst
         fluorids(ji) = zft

         ! Iron and SIO3 saturation concentration from ...
         sieqs(ji) = EXP(  LOG( 10.) * ( 6.44 - 968. / ztkel )  ) * 1.e-6
      END DO
      !
   END SUBROUTINE sed_chem_cst


#endif

END MODULE sedchem
