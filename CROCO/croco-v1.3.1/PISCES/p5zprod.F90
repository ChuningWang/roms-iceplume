#include "cppdefs.h"

MODULE p5zprod
   !!======================================================================
   !!                         ***  MODULE p5zprod  ***
   !! TOP :  Growth Rate of the two phytoplanktons groups 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-05  (O. Aumont, C. Ethe) New parameterization of light limitation
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!   p5z_prod       :   Compute the growth Rate of the two phytoplanktons groups
   !!   p5z_prod_init  :   Initialization of the parameters for growth
   !!   p5z_prod_alloc :   Allocate variables for growth
   !!----------------------------------------------------------------------
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zlim
   USE p5zlim          !  Co-limitations of differents nutrients

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p5z_prod         ! called in p5zbio.F90
   PUBLIC   p5z_prod_init    ! called in trcsms_pisces.F90
   PUBLIC   p5z_prod_alloc

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"

   !! * Shared module variables
   REAL(wp), PUBLIC ::  pislopen        !:
   REAL(wp), PUBLIC ::  pislopep        !:
   REAL(wp), PUBLIC ::  pisloped        !:
   REAL(wp), PUBLIC ::  xadap           !:
   REAL(wp), PUBLIC ::  excretn         !:
   REAL(wp), PUBLIC ::  excretp         !:
   REAL(wp), PUBLIC ::  excretd         !:
   REAL(wp), PUBLIC ::  bresp           !:
   REAL(wp), PUBLIC ::  thetanpm        !:
   REAL(wp), PUBLIC ::  thetannm        !:
   REAL(wp), PUBLIC ::  thetandm        !:
   REAL(wp), PUBLIC ::  chlcmin         !:
   REAL(wp), PUBLIC ::  grosip          !:

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   zdaylen
   
   REAL(wp) :: r1_rday                !: 1 / rday
   REAL(wp) :: texcretn               !: 1 - excret 
   REAL(wp) :: texcretp               !: 1 - excretp 
   REAL(wp) :: texcretd               !: 1 - excret2        

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p5zprod.F90 10872 2019-04-15 12:32:09Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p5z_prod( kt , knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_prod  ***
      !!
      !! ** Purpose :   Compute the phytoplankton production depending on
      !!              light, temperature and nutrient availability
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: kt, knt
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zsilfac, znanotot, zpicotot, zdiattot, zconctemp, zconctemp2
      REAL(wp) ::   zration, zratiop, zratiof, zmax, zmax2, zsilim, ztn, zadap
      REAL(wp) ::   zpronmax, zpropmax, zprofmax, zrat
      REAL(wp) ::   zlim, zsilfac2, zsiborn, zprod, zprontot, zproptot, zprodtot
      REAL(wp) ::   zprnutmax, zdocprod, zprochln, zprochld, zprochlp
      REAL(wp) ::   zpislopen, zpislopep, zpisloped, thetannm_n, thetandm_n, thetanpm_n
      REAL(wp) ::   zrum, zcodel, zargu, zval, zfeup
      REAL(wp) ::   zfact, zrfact2, zmsk
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY) :: zmixnano, zmixpico, zmixdiat, zstrn
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zpislopeadn, zpislopeadp, zpislopeadd
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zprnut, zprmaxp, zprmaxn, zprmaxd
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zprbio, zprpic, zprdia, zysopt
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zprchln, zprchlp, zprchld
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zprorcan, zprorcap, zprorcad 
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zprofed, zprofep, zprofen
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zpronewn, zpronewp, zpronewd
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zproregn, zproregp, zproregd
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zpropo4n, zpropo4p, zpropo4d
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zprodopn, zprodopp, zprodopd
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zrespn, zrespp, zrespd
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zcroissn, zcroissp, zcroissd
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zmxl_fac, zmxl_chl
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zpligprod1, zpligprod2
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zw2d
      !!---------------------------------------------------------------------
      !
      zprorcan(:,:,:) = 0. ; zprorcap(:,:,:) = 0. ; zprorcad(:,:,:) = 0.
      zprofed (:,:,:) = 0. ; zprofep (:,:,:) = 0. ; zprofen (:,:,:) = 0.
      zpronewn(:,:,:) = 0. ; zpronewp(:,:,:) = 0. ; zpronewd(:,:,:) = 0.
      zproregn(:,:,:) = 0. ; zproregp(:,:,:) = 0. ; zproregd(:,:,:) = 0. 
      zpropo4n(:,:,:) = 0. ; zpropo4p(:,:,:) = 0. ; zpropo4d(:,:,:) = 0.
      zprdia  (:,:,:) = 0. ; zprpic  (:,:,:) = 0. ; zprbio  (:,:,:) = 0.
      zprodopn(:,:,:) = 0. ; zprodopp(:,:,:) = 0. ; zprodopd(:,:,:) = 0.
      zysopt  (:,:,:) = 0. ; zmxl_fac(:,:,:) = 0. ; zmxl_chl(:,:,:) = 0.
      zrespn  (:,:,:) = 0. ; zrespp  (:,:,:) = 0. ; zrespd  (:,:,:) = 0. 

      ! Computation of the optimal production
      zprnut (:,:,:) = 0.65 * r1_rday * tgfunc(:,:,:)
      zprmaxn(:,:,:) = ( 0.65 * (1. + zpsino3 * qnpmax ) ) * r1_rday * tgfunc(:,:,:)
      zprmaxp(:,:,:) = 0.5 / 0.65 * zprmaxn(:,:,:) 
      zprmaxd(:,:,:) = zprmaxn(:,:,:) 

      ! compute the day length depending on latitude and the day
      zrum = REAL( nday_year - 80, wp ) / REAL( nyear_len, wp )
      zcodel = ASIN(  SIN( zrum * rpi * 2. ) * SIN( rad * 23.5 )  )

      ! day length in hours
      zstrn(:,:) = 0.
      DO jj = JRANGE
         DO ji = IRANGE
            zargu = TAN( zcodel ) * TAN( gphit(ji,jj) * rad )
            zargu = MAX( -1., MIN(  1., zargu ) )
            zstrn(ji,jj) = MAX( 0.0, 24. - 2. * ACOS( zargu ) / rad / 15. )
         END DO
      END DO

         ! Impact of the day duration on phytoplankton growth
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  zval = MAX( 1., zstrn(ji,jj) )
                  IF( gdept_n(ji,jj,K) <= hmld(ji,jj) ) THEN
                     zval = zval * MIN(1., heup_01(ji,jj) / ( hmld(ji,jj) + rtrn ))
                  ENDIF
                  zmxl_chl(ji,jj,jk) = zval / 24.
                  zmxl_fac(ji,jj,jk) = 1.5 * zval / ( 12. + zval )
               ENDIF
            END DO
         END DO
      END DO

      zprbio(:,:,:) = zprmaxn(:,:,:) * zmxl_fac(:,:,:)
      zprdia(:,:,:) = zprmaxd(:,:,:) * zmxl_fac(:,:,:)
      zprpic(:,:,:) = zprmaxp(:,:,:) * zmxl_fac(:,:,:)


      ! Maximum light intensity
      zdaylen(:,:) = MAX(1., zstrn(:,:)) / 24.
      WHERE( zstrn(:,:) < 1.e0 ) zstrn(:,:) = 24.

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  ! Computation of the P-I slope for nanos and diatoms
                  ztn         = MAX( 0., tsn(ji,jj,K,jp_tem) - 15. )
                  zadap       = xadap * ztn / ( 2.+ ztn )
                  !
                  zpislopeadn(ji,jj,jk) = pislopen * trb(ji,jj,K,jpnch)    &
                  &                       /( trb(ji,jj,K,jpphy) * 12. + rtrn)
                  zpislopeadp(ji,jj,jk) = pislopep * ( 1. + zadap * EXP( -0.25 * epico(ji,jj,jk) ) )   &
                  &                       * trb(ji,jj,K,jppch)   &
                  &                       /( trb(ji,jj,K,jppic) * 12. + rtrn)
                  zpislopeadd(ji,jj,jk) = pisloped * trb(ji,jj,K,jpdch)    &
                     &                    /( trb(ji,jj,K,jpdia) * 12. + rtrn)
                  !
                  zpislopen = zpislopeadn(ji,jj,jk) / ( zprbio(ji,jj,jk) * rday * xlimphy(ji,jj,jk) + rtrn )
                  zpislopep = zpislopeadp(ji,jj,jk) / ( zprpic(ji,jj,jk) * rday * xlimpic(ji,jj,jk) + rtrn )
                  zpisloped = zpislopeadd(ji,jj,jk) / ( zprdia(ji,jj,jk) * rday * xlimdia(ji,jj,jk) + rtrn )

                  ! Computation of production function for Carbon
                  !  ---------------------------------------------
                  zprbio(ji,jj,jk) = zprbio(ji,jj,jk) * ( 1.- EXP( -zpislopen * enano(ji,jj,jk) )  )
                  zprpic(ji,jj,jk) = zprpic(ji,jj,jk) * ( 1.- EXP( -zpislopep * epico(ji,jj,jk) )  )
                  zprdia(ji,jj,jk) = zprdia(ji,jj,jk) * ( 1.- EXP( -zpisloped * ediat(ji,jj,jk) )  )

                  ! Computation of production function for Chlorophyll
                  !  -------------------------------------------------
                  zpislopen = zpislopen * zmxl_fac(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zpisloped = zpisloped * zmxl_fac(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zpislopep = zpislopep * zmxl_fac(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprchln(ji,jj,jk) = zprmaxn(ji,jj,jk) * ( 1.- EXP( -zpislopen * enanom(ji,jj,jk) )  )
                  zprchlp(ji,jj,jk) = zprmaxp(ji,jj,jk) * ( 1.- EXP( -zpislopep * epicom(ji,jj,jk) )  )
                  zprchld(ji,jj,jk) = zprmaxd(ji,jj,jk) * ( 1.- EXP( -zpisloped * ediatm(ji,jj,jk) )  )
               ENDIF
            END DO
         END DO
      END DO

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE

                IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !    Si/C of diatoms
                  !    ------------------------
                  !    Si/C increases with iron stress and silicate availability
                  !    Si/C is arbitrariliy increased for very high Si concentrations
                  !    to mimic the very high ratios observed in the Southern Ocean (silpot2)
                  zlim  = trb(ji,jj,K,jpsil) / ( trb(ji,jj,K,jpsil) + xksi1 )
                  zsilim = MIN( zprdia(ji,jj,jk) / ( zprmaxd(ji,jj,jk) + rtrn ), xlimsi(ji,jj,jk) )
                  zsilfac = 3.4 * EXP( -4.23 * zsilim ) * MAX( 0.e0, MIN( 1., 2.2 * ( zlim - 0.5 ) )  ) + 1.e0
                  zsiborn = trb(ji,jj,K,jpsil) * trb(ji,jj,K,jpsil)   &
                  &        * trb(ji,jj,K,jpsil)
                  IF (gphit(ji,jj) < -30 ) THEN
                    zsilfac2 = 1. + 2. * zsiborn / ( zsiborn + xksi2**3 )
                  ELSE
                    zsilfac2 = 1. +      zsiborn / ( zsiborn + xksi2**3 )
                  ENDIF
                  zysopt(ji,jj,jk) = grosip * zlim * zsilfac * zsilfac2
              ENDIF
            END DO
         END DO
      END DO

      !  Sea-ice effect on production                                                                               
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zprbio(ji,jj,jk)  = zprbio(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
               zprpic(ji,jj,jk)  = zprpic(ji,jj,jk) * ( 1. - fr_i(ji,jj) ) 
               zprdia(ji,jj,jk)  = zprdia(ji,jj,jk) * ( 1. - fr_i(ji,jj) ) 
               zprnut(ji,jj,jk)  = zprnut(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
            END DO
         END DO
      END DO

      ! Computation of the various production terms of nanophytoplankton 
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production terms for nanophyto.
                  zprorcan(ji,jj,jk) = zprbio(ji,jj,jk)  * xlimphy(ji,jj,jk) * trb(ji,jj,K,jpphy) * rfact2
                  !
                  zration = trb(ji,jj,K,jpnph) / ( trb(ji,jj,K,jpphy) + rtrn )
                  zratiop = trb(ji,jj,K,jppph) / ( trb(ji,jj,K,jpphy) + rtrn )
                  zratiof = trb(ji,jj,K,jpnfe) / ( trb(ji,jj,K,jpphy) + rtrn )
                  zprnutmax = zprnut(ji,jj,jk) * fvnuptk(ji,jj,jk) / rno3 * trb(ji,jj,K,jpphy) * rfact2
                  ! Uptake of nitrogen
                  zrat = MIN( 1., zration / (xqnnmax(ji,jj,jk) + rtrn) ) 
                  zmax = MAX(0., MIN(1., (1. - zrat)/ (1.05 - zrat) * 1.05))
                  zpronmax = zprnutmax * zmax * MAX(0., MIN(1., ( zratiop - xqpnmin(ji,jj,jk) )   &
                  &          / ( xqpnmax(ji,jj,jk) - xqpnmin(ji,jj,jk) + rtrn ), xlimnfe(ji,jj,jk) ) )
                  zpronewn(ji,jj,jk) = zpronmax * zdaylen(ji,jj) * xnanono3(ji,jj,jk)
                  zproregn(ji,jj,jk) = zpronmax * xnanonh4(ji,jj,jk)
                  ! Uptake of phosphorus
                  zrat = MIN( 1., zratiop / (xqpnmax(ji,jj,jk) + rtrn) )
                  zmax = MAX(0., MIN(1., (1. - zrat)/ (1.05 - zrat) * 1.05))
                  zpropmax = zprnutmax * zmax * xlimnfe(ji,jj,jk)
                  zpropo4n(ji,jj,jk) = zpropmax * xnanopo4(ji,jj,jk)
                  zprodopn(ji,jj,jk) = zpropmax * xnanodop(ji,jj,jk)
                  ! Uptake of iron
                  zrat = MIN( 1., zratiof / qfnmax )
                  zmax = MAX(0., MIN(1., (1. - zrat)/ (1.05 - zrat) * 1.05))
                  zprofmax = zprnutmax * qfnmax * zmax
                  zprofen(ji,jj,jk) = zprofmax * xnanofer(ji,jj,jk) * ( 3. - 2.4 * xlimnfe(ji,jj,jk)    &
                  &          / ( xlimnfe(ji,jj,jk) + 0.2 ) ) * (1. + 0.8 * xnanono3(ji,jj,jk) / ( rtrn  &
                  &          + xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk) ) * (1. - xnanofer(ji,jj,jk) ) )
               ENDIF
            END DO
         END DO
      END DO

      ! Computation of the various production terms of picophytoplankton 
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production terms for picophyto.
                  zprorcap(ji,jj,jk) = zprpic(ji,jj,jk)  * xlimpic(ji,jj,jk) * trb(ji,jj,K,jppic) * rfact2
                  !
                  zration = trb(ji,jj,K,jpnpi) / ( trb(ji,jj,K,jppic) + rtrn )
                  zratiop = trb(ji,jj,K,jpppi) / ( trb(ji,jj,K,jppic) + rtrn )
                  zratiof = trb(ji,jj,K,jppfe) / ( trb(ji,jj,K,jppic) + rtrn )
                  zprnutmax = zprnut(ji,jj,jk) * fvpuptk(ji,jj,jk) / rno3 * trb(ji,jj,K,jppic) * rfact2
                  ! Uptake of nitrogen
                  zrat = MIN( 1., zration / (xqnpmax(ji,jj,jk) + rtrn) )
                  zmax = MAX(0., MIN(1., (1. - zrat)/ (1.05 - zrat) * 1.05))
                  zpronmax = zprnutmax * zmax * MAX(0., MIN(1., ( zratiop - xqppmin(ji,jj,jk) )   &
                  &          / ( xqppmax(ji,jj,jk) - xqppmin(ji,jj,jk) + rtrn ), xlimpfe(ji,jj,jk) ) )
                  zpronewp(ji,jj,jk) = zpronmax * zdaylen(ji,jj) * xpicono3(ji,jj,jk) 
                  zproregp(ji,jj,jk) = zpronmax * xpiconh4(ji,jj,jk)
                  ! Uptake of phosphorus
                  zrat = MIN( 1., zratiop / (xqppmax(ji,jj,jk) + rtrn) )
                  zmax = MAX(0., MIN(1., (1. - zrat)/ (1.05 - zrat) * 1.05))
                  zpropmax = zprnutmax * zmax * xlimpfe(ji,jj,jk)
                  zpropo4p(ji,jj,jk) = zpropmax * xpicopo4(ji,jj,jk)
                  zprodopp(ji,jj,jk) = zpropmax * xpicodop(ji,jj,jk)
                  ! Uptake of iron
                  zrat = MIN( 1., zratiof / qfpmax )
                  zmax = MAX(0., MIN(1., (1. - zrat)/ (1.05 - zrat) * 1.05))
                  zprofmax = zprnutmax * qfpmax * zmax
                  zprofep(ji,jj,jk) = zprofmax * xpicofer(ji,jj,jk) * ( 3. - 2.4 * xlimpfe(ji,jj,jk)   &
                  &          / ( xlimpfe(ji,jj,jk) + 0.2 ) ) * (1. + 0.8 * xpicono3(ji,jj,jk) / ( rtrn   &
                  &          + xpicono3(ji,jj,jk) + xpiconh4(ji,jj,jk) ) * (1. - xpicofer(ji,jj,jk) ) )
               ENDIF
            END DO
         END DO
      END DO

      ! Computation of the various production terms of diatoms
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production terms for diatomees
                  zprorcad(ji,jj,jk) = zprdia(ji,jj,jk) * xlimdia(ji,jj,jk) * trb(ji,jj,K,jpdia) * rfact2
                  ! Computation of the respiration term according to pahlow 
                  ! & oschlies (2013)
                  !
                  zration = trb(ji,jj,K,jpndi) / ( trb(ji,jj,K,jpdia) + rtrn )
                  zratiop = trb(ji,jj,K,jppdi) / ( trb(ji,jj,K,jpdia) + rtrn )
                  zratiof = trb(ji,jj,K,jpdfe) / ( trb(ji,jj,K,jpdia) + rtrn )
                  zprnutmax = zprnut(ji,jj,jk) * fvduptk(ji,jj,jk) / rno3 * trb(ji,jj,K,jpdia) * rfact2
                  ! Uptake of nitrogen
                  zrat = MIN( 1., zration / (xqndmax(ji,jj,jk) + rtrn) )
                  zmax = MAX(0., MIN(1., (1. - zrat)/ (1.05 - zrat) * 1.05)) 
                  zpronmax = zprnutmax * zmax * MAX(0., MIN(1., ( zratiop - xqpdmin(ji,jj,jk) )   &
                  &          / ( xqpdmax(ji,jj,jk) - xqpdmin(ji,jj,jk) + rtrn ), xlimdfe(ji,jj,jk) ) )
                  zpronewd(ji,jj,jk) = zpronmax * zdaylen(ji,jj) * xdiatno3(ji,jj,jk)
                  zproregd(ji,jj,jk) = zpronmax * xdiatnh4(ji,jj,jk)
                  ! Uptake of phosphorus
                  zrat = MIN( 1., zratiop / (xqpdmax(ji,jj,jk) + rtrn) )
                  zmax = MAX(0., MIN(1., (1. - zrat)/ (1.05 - zrat) * 1.05)) 
                  zpropmax = zprnutmax * zmax * xlimdfe(ji,jj,jk)
                  zpropo4d(ji,jj,jk) = zpropmax * xdiatpo4(ji,jj,jk)
                  zprodopd(ji,jj,jk) = zpropmax * xdiatdop(ji,jj,jk)
                  ! Uptake of iron
                  zrat = MIN( 1., zratiof / qfdmax )
                  zmax = MAX(0., MIN(1., (1. - zrat)/ (1.05 - zrat) * 1.05))
                  zprofmax = zprnutmax * qfdmax * zmax
                  zprofed(ji,jj,jk) = zprofmax * xdiatfer(ji,jj,jk) * ( 3. - 2.4 * xlimdfe(ji,jj,jk)     &
                  &          / ( xlimdfe(ji,jj,jk) + 0.2 ) ) * (1. + 0.8 * xdiatno3(ji,jj,jk) / ( rtrn   &
                  &          + xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk) ) * (1. - xdiatfer(ji,jj,jk) ) )
               ENDIF
            END DO
         END DO
      END DO

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                     !  production terms for nanophyto. ( chlorophyll )
                  znanotot = enanom(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprod = rday * (zpronewn(ji,jj,jk) + zproregn(ji,jj,jk)) * zprchln(ji,jj,jk) * xlimphy(ji,jj,jk)
                  thetannm_n   = MIN ( thetannm, ( thetannm / (1. - 1.14 / 43.4 *tsn(ji,jj,K,jp_tem)))   &
                  &               * (1. - 1.14 / 43.4 * 20.))
                  zprochln = thetannm_n * zprod / ( zpislopeadn(ji,jj,jk) * znanotot + rtrn )
                  zprochln = MAX(zprochln, chlcmin * 12. * zprorcan (ji,jj,jk) )
                     !  production terms for picophyto. ( chlorophyll )
                  zpicotot = epicom(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprod = rday * (zpronewp(ji,jj,jk) + zproregp(ji,jj,jk)) * zprchlp(ji,jj,jk) * xlimpic(ji,jj,jk)
                  thetanpm_n   = MIN ( thetanpm, ( thetanpm / (1. - 1.14 / 43.4 *tsn(ji,jj,K,jp_tem)))   &
                  &               * (1. - 1.14 / 43.4 * 20.))
                  zprochlp = thetanpm_n * zprod / ( zpislopeadp(ji,jj,jk) * zpicotot + rtrn )
                  zprochlp = MAX(zprochlp, chlcmin * 12. * zprorcap(ji,jj,jk) )
                  !  production terms for diatomees ( chlorophyll )
                  zdiattot = ediatm(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprod = rday * (zpronewd(ji,jj,jk) + zproregd(ji,jj,jk)) * zprchld(ji,jj,jk) * xlimdia(ji,jj,jk)
                  thetandm_n   = MIN ( thetandm, ( thetandm / (1. - 1.14 / 43.4 *tsn(ji,jj,K,jp_tem)))   &
                  &               * (1. - 1.14 / 43.4 * 20.))
                  zprochld = thetandm_n * zprod / ( zpislopeadd(ji,jj,jk) * zdiattot + rtrn )
                  zprochld = MAX(zprochld, chlcmin * 12. * zprorcad(ji,jj,jk) )
                  !   Update the arrays TRA which contain the Chla sources and sinks
                  tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) + zprochln * texcretn
                  tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) + zprochld * texcretd
                  tra(ji,jj,jk,jppch) = tra(ji,jj,jk,jppch) + zprochlp * texcretp
               ENDIF
            END DO
         END DO
      END DO

      !   Update the arrays TRA which contain the biological sources and sinks
      DO jk = KRANGE
         DO jj = JRANGE
           DO ji = IRANGE
              zprontot = zpronewn(ji,jj,jk) + zproregn(ji,jj,jk)
              zproptot = zpronewp(ji,jj,jk) + zproregp(ji,jj,jk)
              zprodtot = zpronewd(ji,jj,jk) + zproregd(ji,jj,jk)
              zdocprod = excretd * zprorcad(ji,jj,jk) + excretn * zprorcan(ji,jj,jk)  &
              &          + excretp * zprorcap(ji,jj,jk)
              tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) - zpropo4n(ji,jj,jk) - zpropo4d(ji,jj,jk)  &
              &                     - zpropo4p(ji,jj,jk)
              tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - zpronewn(ji,jj,jk) - zpronewd(ji,jj,jk)  &
              &                     - zpronewp(ji,jj,jk)
              tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zproregn(ji,jj,jk) - zproregd(ji,jj,jk)  &
              &                     - zproregp(ji,jj,jk)
              tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zprorcan(ji,jj,jk) * texcretn    &
                 &                  - zpsino3 * zpronewn(ji,jj,jk) - zpsinh4 * zproregn(ji,jj,jk)   &
                 &                  - zrespn(ji,jj,jk) 
              zcroissn(ji,jj,jk) = tra(ji,jj,jk,jpphy) / rfact2/ (trb(ji,jj,K,jpphy) + rtrn)
              tra(ji,jj,jk,jpnph) = tra(ji,jj,jk,jpnph) + zprontot * texcretn
              tra(ji,jj,jk,jppph) = tra(ji,jj,jk,jppph) + zpropo4n(ji,jj,jk) * texcretn   &
              &                     + zprodopn(ji,jj,jk) * texcretn
              tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) + zprofen(ji,jj,jk) * texcretn
              tra(ji,jj,jk,jppic) = tra(ji,jj,jk,jppic) + zprorcap(ji,jj,jk) * texcretp     &
                 &                  - zpsino3 * zpronewp(ji,jj,jk) - zpsinh4 * zproregp(ji,jj,jk)   &
                 &                  - zrespp(ji,jj,jk) 
              zcroissp(ji,jj,jk) = tra(ji,jj,jk,jppic) / rfact2/ (trb(ji,jj,K,jppic) + rtrn)
              tra(ji,jj,jk,jpnpi) = tra(ji,jj,jk,jpnpi) + zproptot * texcretp
              tra(ji,jj,jk,jpppi) = tra(ji,jj,jk,jpppi) + zpropo4p(ji,jj,jk) * texcretp   &
              &                     + zprodopp(ji,jj,jk) * texcretp
              tra(ji,jj,jk,jppfe) = tra(ji,jj,jk,jppfe) + zprofep(ji,jj,jk) * texcretp
              tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) + zprorcad(ji,jj,jk) * texcretd   &
                 &                  - zpsino3 * zpronewd(ji,jj,jk) - zpsinh4 * zproregd(ji,jj,jk)   &
                 &                  - zrespd(ji,jj,jk) 
              zcroissd(ji,jj,jk) = tra(ji,jj,jk,jpdia) / rfact2 / (trb(ji,jj,K,jpdia) + rtrn)
              tra(ji,jj,jk,jpndi) = tra(ji,jj,jk,jpndi) + zprodtot * texcretd
              tra(ji,jj,jk,jppdi) = tra(ji,jj,jk,jppdi) + zpropo4d(ji,jj,jk) * texcretd   &
              &                     + zprodopd(ji,jj,jk) * texcretd
              tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) + zprofed(ji,jj,jk) * texcretd
              tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zprorcad(ji,jj,jk) * zysopt(ji,jj,jk) * texcretd
              tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + excretd * zprorcad(ji,jj,jk) + excretn * zprorcan(ji,jj,jk)  &
              &                     + excretp * zprorcap(ji,jj,jk)
              tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + excretd * zprodtot + excretn * zprontot   &
              &                     + excretp * zproptot
              tra(ji,jj,jk,jpdop) = tra(ji,jj,jk,jpdop) + excretd * zpropo4d(ji,jj,jk) + excretn * zpropo4n(ji,jj,jk)   &
              &    - texcretn * zprodopn(ji,jj,jk) - texcretd * zprodopd(ji,jj,jk) + excretp * zpropo4p(ji,jj,jk)     &
              &    - texcretp * zprodopp(ji,jj,jk)
              tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + o2ut * ( zproregn(ji,jj,jk) + zproregd(ji,jj,jk)   &
                 &                + zproregp(ji,jj,jk) ) + ( o2ut + o2nit ) * ( zpronewn(ji,jj,jk)           &
                 &                + zpronewd(ji,jj,jk) + zpronewp(ji,jj,jk) )   &
                 &                - o2ut * ( zrespn(ji,jj,jk) + zrespp(ji,jj,jk) + zrespd(ji,jj,jk) )
              zfeup = texcretn * zprofen(ji,jj,jk) + texcretd * zprofed(ji,jj,jk) + texcretp * zprofep(ji,jj,jk)
              tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zfeup
              tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) - texcretd * zprorcad(ji,jj,jk) * zysopt(ji,jj,jk)
              tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprorcan(ji,jj,jk) - zprorcad(ji,jj,jk) - zprorcap(ji,jj,jk)  &
              &                     + zpsino3 * zpronewn(ji,jj,jk) + zpsinh4 * zproregn(ji,jj,jk)   &
              &                     + zpsino3 * zpronewp(ji,jj,jk) + zpsinh4 * zproregp(ji,jj,jk)   &
              &                     + zpsino3 * zpronewd(ji,jj,jk) + zpsinh4 * zproregd(ji,jj,jk)  &
              &                     + zrespn(ji,jj,jk) + zrespd(ji,jj,jk) + zrespp(ji,jj,jk) 
              tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * ( zpronewn(ji,jj,jk) + zpronewd(ji,jj,jk)  &
              &                     + zpronewp(ji,jj,jk) ) - rno3 * ( zproregn(ji,jj,jk) + zproregd(ji,jj,jk)     &
              &                     + zproregp(ji,jj,jk) ) 
          END DO
        END DO
     END DO
     !
     IF( ln_ligand ) THEN
         zpligprod1(:,:,:) = 0.    ;    zpligprod2(:,:,:) = 0. 
         DO jk = KRANGE
            DO jj = JRANGE
              DO ji = IRANGE
                 zdocprod = excretd * zprorcad(ji,jj,jk) + excretn * zprorcan(ji,jj,jk) + excretp * zprorcap(ji,jj,jk)
                 zfeup    = texcretn * zprofen(ji,jj,jk) + texcretd * zprofed(ji,jj,jk) + texcretp * zprofep(ji,jj,jk)
                 tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) + zdocprod * ldocp  &
                 &       - zfeup * plig(ji,jj,jk) / ( rtrn + plig(ji,jj,jk) + 2.E3 * (1.0 - plig(ji,jj,jk) ) )
                 zpligprod1(ji,jj,jk) = zdocprod * ldocp
                 zpligprod2(ji,jj,jk) = zfeup * plig(ji,jj,jk) * lthet
              END DO
           END DO
        END DO
     ENDIF


     ! Total primary production per year
#if defined key_iomput
    IF( lk_iomput ) THEN
       IF( knt == nrdttrc ) THEN
          ALLOCATE( zw2d(PRIV_2D_BIOARRAY), zw3d(PRIV_3D_BIOARRAY) )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "PPPHYN" ) .OR. iom_use( "PPPHYD" ) .OR. iom_use( "PPPHYP" ) )  THEN
              zw3d(:,:,:) = zprorcan(:,:,:) * zfact * tmask(:,:,:)  ! primary production by nanophyto
              CALL iom_put( "PPPHYN"  , zw3d )
              !
              zw3d(:,:,:) = zprorcap(:,:,:) * zfact * tmask(:,:,:)  ! primary production by picophyto
              CALL iom_put( "PPPHYP"  , zw3d )
              !
              zw3d(:,:,:) = zprorcad(:,:,:) * zfact * tmask(:,:,:)  ! primary production by diatomes
              CALL iom_put( "PPPHYD"  , zw3d )
          ENDIF
          IF( iom_use( "PPNEWN" ) .OR. iom_use( "PPNEWD" ) .OR. iom_use( "PPNEWP" ) )  THEN
              zw3d(:,:,:) = zpronewn(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by nanophyto
              CALL iom_put( "PPNEWN"  , zw3d )
              !
              zw3d(:,:,:) = zpronewp(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by picophyto
              CALL iom_put( "PPNEWP"  , zw3d )
              !
              zw3d(:,:,:) = zpronewd(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by diatomes
              CALL iom_put( "PPNEWD"  , zw3d )
          ENDIF
          IF( iom_use( "PBSi" ) )  THEN
              zw3d(:,:,:) = zprorcad(:,:,:) * zfact * tmask(:,:,:) * zysopt(:,:,:) ! biogenic silica production
              CALL iom_put( "PBSi"  , zw3d )
          ENDIF
          IF( iom_use( "PFeN" ) .OR. iom_use( "PFeD" ) .OR. iom_use( "PFeP" ) )  THEN
              zw3d(:,:,:) = zprofen(:,:,:) * zfact * tmask(:,:,:)  ! biogenic iron production by nanophyto
              CALL iom_put( "PFeN"  , zw3d )
              !
              zw3d(:,:,:) = zprofep(:,:,:) * zfact * tmask(:,:,:)  ! biogenic iron production by picophyto
              CALL iom_put( "PFeP"  , zw3d )
              !
              zw3d(:,:,:) = zprofed(:,:,:) * zfact * tmask(:,:,:)  ! biogenic iron production by  diatomes
              CALL iom_put( "PFeD"  , zw3d )
          ENDIF
          IF( iom_use( "LPRODP" ) )  THEN
              zw3d(:,:,:) = zpligprod1(:,:,:) * 1e9 * zfact * tmask(:,:,:)
              CALL iom_put( "LPRODP"  , zw3d )
          ENDIF
          IF( iom_use( "LDETP" ) )  THEN
              zw3d(:,:,:) = zpligprod2(:,:,:) * 1e9 * zfact * tmask(:,:,:)
              CALL iom_put( "LDETP"  , zw3d )
          ENDIF
          IF( iom_use( "Mumax" ) )  THEN
              zw3d(:,:,:) = zprmaxn(:,:,:) * tmask(:,:,:)   ! Maximum growth rate
              CALL iom_put( "Mumax"  , zw3d )
          ENDIF
          IF( iom_use( "MuN" ) .OR. iom_use( "MuD" ) .OR. iom_use( "MuP" ) )  THEN
              zw3d(:,:,:) = zprbio(:,:,:) * xlimphy(:,:,:) * tmask(:,:,:)  ! Realized growth rate for nanophyto
              CALL iom_put( "MuN"  , zw3d )
              !
              zw3d(:,:,:) = zprpic(:,:,:) * xlimpic(:,:,:) * tmask(:,:,:)  ! Realized growth rate for picophyto
              CALL iom_put( "MuP"  , zw3d )
              !
              zw3d(:,:,:) =  zprdia(:,:,:) * xlimdia(:,:,:) * tmask(:,:,:)  ! Realized growth rate for diatoms
              CALL iom_put( "MuD"  , zw3d )
          ENDIF
          IF( iom_use( "LNlight" ) .OR. iom_use( "LDlight" ) .OR. iom_use( "LPlight" ) )  THEN
              zw3d(:,:,:) = zprbio (:,:,:) / (zprmaxn(:,:,:) + rtrn) * tmask(:,:,:) ! light limitation term
              CALL iom_put( "LNlight"  , zw3d )
              !
              zw3d(:,:,:) = zprpic (:,:,:) / (zprmaxp(:,:,:) + rtrn) * tmask(:,:,:) ! light limitation term
              CALL iom_put( "LPlight"  , zw3d )
              !
              zw3d(:,:,:) =  zprdia (:,:,:) / (zprmaxd(:,:,:) + rtrn) * tmask(:,:,:)  ! light limitation term
              CALL iom_put( "LDlight"  , zw3d )
          ENDIF
          IF( iom_use( "MunetN" ) .OR. iom_use( "MunetD" ) .OR. iom_use( "MunetP" ) )  THEN
              zw3d(:,:,:) = zcroissn(:,:,:) * tmask(:,:,:) ! ! Realized growth rate for nanophyto
              CALL iom_put( "MunetN"  , zw3d )
              !
              zw3d(:,:,:) = zcroissp(:,:,:) * tmask(:,:,:) ! ! Realized growth rate for picophyto
              CALL iom_put( "MunetP"  , zw3d )
              !
              zw3d(:,:,:) = zcroissd(:,:,:) * tmask(:,:,:) ! ! Realized growth rate for diatomes
              CALL iom_put( "MunetD"  , zw3d )
              !
          ENDIF
          !
          DEALLOCATE( zw2d, zw3d )
       ENDIF
     ENDIF
#endif
     !
#if defined key_trc_diaadd
      !   Supplementary diagnostics
     zfact = 1.e3 * rfact2r
     DO jk = KRANGE
        DO jj = JRANGE
          DO ji = IRANGE
             zmsk = zfact * tmask(ji,jj,K)
             trc3d(ji,jj,K,jp_pphy  )  = zprorcan(ji,jj,jk) * zmsk  ! primary production by nanophyto
             trc3d(ji,jj,K,jp_pphy2 )  = zprorcad(ji,jj,jk) * zmsk  ! primary production by diatom
             trc3d(ji,jj,K,jp_pnew  )  = zpronewn(ji,jj,jk) * zmsk ! new primary production by nanophyto
             trc3d(ji,jj,K,jp_pnew2 )  = zpronewd(ji,jj,jk) * zmsk ! new primary production by diatom
             trc3d(ji,jj,K,jp_pbsi  )  = zprorcad(ji,jj,jk) * zysopt(ji,jj,jk) * zmsk ! biogenic silica production
             trc3d(ji,jj,K,jp_pfed  )  = zprofed (ji,jj,jk) * zmsk  ! biogenic iron production by diatom
             trc3d(ji,jj,K,jp_pfen  )  = zprofen (ji,jj,jk) * zmsk!  biogenic iron production by nanophyto
             trc3d(ji,jj,K,jp_pnewo2)  = ( o2ut + o2nit ) &  ! Oxygen production by the New Produc.
                &                      * ( zpronewn(ji,jj,jk) + zpronewd(ji,jj,jk) + zpronewp(ji,jj,jk) ) * zmsk
             trc3d(ji,jj,K,jp_prego2)  = o2ut * &        ! Oxygen production by the Regen Produc.
                &                 (   zprorcan(ji,jj,jk) - zpronewn(ji,jj,jk)  &
                &                  + zprorcad(ji,jj,jk) - zpronewd(ji,jj,jk)   &
                &                  + zprorcap(ji,jj,jk) - zpronewp(ji,jj,jk) ) * zmsk
         END DO
        END DO
      END DO
#endif



      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('prod')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( charout, ltra='tra')
      ENDIF
      !
   END SUBROUTINE p5z_prod


   SUBROUTINE p5z_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p5z_prod_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton production parameters
      !!
      !! ** Method  :   Read the nampisprod namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampisprod
      !!----------------------------------------------------------------------
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!
      NAMELIST/namp5zprod/ pislopen, pislopep, pisloped, excretn, excretp, excretd,     &
         &                 thetannm, thetanpm, thetandm, chlcmin, grosip, bresp, xadap
      !!----------------------------------------------------------------------

      REWIND( numnatp_ref )              ! Namelist nampisprod in reference namelist : Pisces phytoplankton production
      READ  ( numnatp_ref, namp5zprod, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namp5zprod in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisprod in configuration namelist : Pisces phytoplankton production
      READ  ( numnatp_cfg, namp5zprod, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'namp5zprod in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, namp5zprod )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton growth, namp5zprod'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    mean Si/C ratio                           grosip       =', grosip
         WRITE(numout,*) '    P-I slope                                 pislopen     =', pislopen
         WRITE(numout,*) '    P-I slope  for diatoms                    pisloped     =', pisloped
         WRITE(numout,*) '    P-I slope  for picophytoplankton          pislopep     =', pislopep
         WRITE(numout,*) '    Acclimation factor to low light           xadap        =', xadap
         WRITE(numout,*) '    excretion ratio of nanophytoplankton      excretn      =', excretn
         WRITE(numout,*) '    excretion ratio of picophytoplankton      excretp      =', excretp
         WRITE(numout,*) '    excretion ratio of diatoms                excretd      =', excretd
         WRITE(numout,*) '    basal respiration in phytoplankton        bresp        =', bresp
         WRITE(numout,*) '    Maximum Chl/C in phytoplankton            chlcmin      =', chlcmin
         WRITE(numout,*) '    Minimum Chl/N in nanophytoplankton        thetannm     =', thetannm
         WRITE(numout,*) '    Minimum Chl/N in picophytoplankton        thetanpm     =', thetanpm
         WRITE(numout,*) '    Minimum Chl/N in diatoms                  thetandm     =', thetandm
      ENDIF
      !
      r1_rday   = 1. / rday 
      texcretn  = 1. - excretn
      texcretp  = 1. - excretp
      texcretd  = 1. - excretd
      !
   END SUBROUTINE p5z_prod_init


   INTEGER FUNCTION p5z_prod_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_prod_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( zdaylen(PRIV_2D_BIOARRAY), STAT = p5z_prod_alloc )
      !
      IF( p5z_prod_alloc /= 0 ) CALL ctl_warn( 'p5z_prod_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p5z_prod_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p5z_prod                    ! Empty routine
   END SUBROUTINE p5z_prod
#endif

   !!======================================================================
END MODULE p5zprod
