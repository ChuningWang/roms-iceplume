#include "cppdefs.h"

MODULE p4zprod
   !!======================================================================
   !!                         ***  MODULE p4zprod  ***
   !! TOP :  Growth Rate of the two phytoplanktons groups 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-05  (O. Aumont, C. Ethe) New parameterization of light limitation
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!   p4z_prod       : Compute the growth Rate of the two phytoplanktons groups
   !!   p4z_prod_init  : Initialization of the parameters for growth
   !!   p4z_prod_alloc : Allocate variables for growth
   !!----------------------------------------------------------------------
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p4zlim          ! Co-limitations of differents nutrients
!  USE prtctl_trc      ! print control for debugging
!  USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_prod         ! called in p4zbio.F90
   PUBLIC   p4z_prod_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_prod_alloc

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"

   REAL(wp), PUBLIC ::   pislopen     !:
   REAL(wp), PUBLIC ::   pisloped     !:
   REAL(wp), PUBLIC ::   xadap        !:
   REAL(wp), PUBLIC ::   excretn      !:
   REAL(wp), PUBLIC ::   excretd      !:
   REAL(wp), PUBLIC ::   bresp        !:
   REAL(wp), PUBLIC ::   chlcnm       !:
   REAL(wp), PUBLIC ::   chlcdm       !:
   REAL(wp), PUBLIC ::   chlcmin      !:
   REAL(wp), PUBLIC ::   fecnm        !:
   REAL(wp), PUBLIC ::   fecdm        !:
   REAL(wp), PUBLIC ::   grosip       !:

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotan   !: proxy of N quota in Nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotad   !: proxy of N quota in diatomee
   
   REAL(wp) ::   r1_rday    ! 1 / rday
   REAL(wp) ::   texcretn   ! 1 - excretn 
   REAL(wp) ::   texcretd   ! 1 - excretd        

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zprod.F90 11117 2019-06-17 08:50:02Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_prod( kt , knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod  ***
      !!
      !! ** Purpose :   Compute the phytoplankton production depending on
      !!              light, temperature and nutrient availability
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   !
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zsilfac, znanotot, zdiattot, zconctemp, zconctemp2
      REAL(wp) ::   zratio, zmax, zsilim, ztn, zadap, zlim, zsilfac2, zsiborn
      REAL(wp) ::   zprod, zproreg, zproreg2, zprochln, zprochld
      REAL(wp) ::   zmaxday, zdocprod, zpislopen, zpisloped
      REAL(wp) ::   zmxltst, zmxlday
      REAL(wp) ::   zrum, zcodel, zargu, zval, zfeup, chlcnm_n, chlcdm_n
      REAL(wp) ::   zfact, zmsk
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY) :: zw2d
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zw3d
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY) :: zstrn, zmixnano, zmixdiat
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zprmaxn,zprmaxd
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zpislopeadn, zpislopeadd, zysopt 
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zprdia, zprbio, zprdch, zprnch   
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zprorcan, zprorcad, zprofed, zprofen
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zpronewn, zpronewd
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zmxl_fac, zmxl_chl
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zpligprod1, zpligprod2
      !!---------------------------------------------------------------------
      !
      !  Allocate temporary workspace
      !
      zprorcan(:,:,:) = 0. ; zprorcad(:,:,:) = 0. ; zprofed (:,:,:) = 0.
      zprofen (:,:,:) = 0. ; zysopt  (:,:,:) = 0.
      zpronewn(:,:,:) = 0. ; zpronewd(:,:,:) = 0. ; zprdia  (:,:,:) = 0.
      zprbio  (:,:,:) = 0. ; zprdch  (:,:,:) = 0. ; zprnch  (:,:,:) = 0. 
      zmxl_fac(:,:,:) = 0. ; zmxl_chl(:,:,:) = 0. 

      ! Computation of the optimal production
      zprmaxn(:,:,:) = 0.8 * r1_rday * tgfunc(:,:,:)
      zprmaxd(:,:,:) = zprmaxn(:,:,:)

      ! compute the day length depending on latitude and the day
      zrum = FLOAT( nday_year - 80 ) / nyear_len
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

      ! Impact of the day duration and light intermittency on phytoplankton growth
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

      ! Maximum light intensity
      WHERE( zstrn(:,:) < 1.e0 ) zstrn(:,:) = 24.

      ! Computation of the P-I slope for nanos and diatoms
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  ztn         = MAX( 0., tsn(ji,jj,K,jp_tem) - 15. )
                  zadap       = xadap * ztn / ( 2.+ ztn )
                  zconctemp   = MAX( 0.e0 , trb(ji,jj,K,jpdia) - xsizedia )
                  zconctemp2  = trb(ji,jj,K,jpdia) - zconctemp
                  !
                  zpislopeadn(ji,jj,jk) = pislopen * ( 1.+ zadap  * EXP( -0.25 * enano(ji,jj,jk) ) )  &
                  &    * trb(ji,jj,K,jpnch) /( trb(ji,jj,K,jpphy) * 12. + rtrn)
                  !
                  zpislopeadd(ji,jj,jk) = (pislopen * zconctemp2 + pisloped * zconctemp)  &
                  &    / ( trb(ji,jj,K,jpdia) + rtrn )   &
                  &    * trb(ji,jj,K,jpdch) /( trb(ji,jj,K,jpdia) * 12. + rtrn)
               ENDIF
            END DO
         END DO
      END DO

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                   ! Computation of production function for Carbon
                   !  ---------------------------------------------
                   zpislopen = zpislopeadn(ji,jj,jk) / ( ( r1_rday + bresp * r1_rday ) &
                   &            * zmxl_fac(ji,jj,jk) * rday + rtrn)
                   zpisloped = zpislopeadd(ji,jj,jk) / ( ( r1_rday + bresp * r1_rday ) &
                   &            * zmxl_fac(ji,jj,jk) * rday + rtrn)
                   zprbio(ji,jj,jk) = zprbio(ji,jj,jk) * ( 1.- EXP( -zpislopen * enano(ji,jj,jk) )  )
                   zprdia(ji,jj,jk) = zprdia(ji,jj,jk) * ( 1.- EXP( -zpisloped * ediat(ji,jj,jk) )  )
                   !  Computation of production function for Chlorophyll
                   !--------------------------------------------------
                   zpislopen = zpislopeadn(ji,jj,jk) / ( zprmaxn(ji,jj,jk) * zmxl_chl(ji,jj,jk) * rday + rtrn )
                   zpisloped = zpislopeadd(ji,jj,jk) / ( zprmaxd(ji,jj,jk) * zmxl_chl(ji,jj,jk) * rday + rtrn )
                   zprnch(ji,jj,jk) = zprmaxn(ji,jj,jk) * ( 1.- EXP( -zpislopen * enanom(ji,jj,jk) ) )
                   zprdch(ji,jj,jk) = zprmaxd(ji,jj,jk) * ( 1.- EXP( -zpisloped * ediatm(ji,jj,jk) ) )
               ENDIF
            END DO
         END DO
      END DO

      !  Computation of a proxy of the N/C ratio
      !  ---------------------------------------
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
                zval = MIN( xnanopo4(ji,jj,jk), ( xnanonh4(ji,jj,jk) + xnanono3(ji,jj,jk) ) )   &
                &      * zprmaxn(ji,jj,jk) / ( zprbio(ji,jj,jk) + rtrn )
                quotan(ji,jj,jk) = MIN( 1., 0.2 + 0.8 * zval )
                zval = MIN( xdiatpo4(ji,jj,jk), ( xdiatnh4(ji,jj,jk) + xdiatno3(ji,jj,jk) ) )   &
                &      * zprmaxd(ji,jj,jk) / ( zprdia(ji,jj,jk) + rtrn )
                quotad(ji,jj,jk) = MIN( 1., 0.2 + 0.8 * zval )
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
                  zsilfac = 4.4 * EXP( -4.23 * zsilim ) * MAX( 0.e0, MIN( 1., 2.2 * ( zlim - 0.5 ) )  ) + 1.e0
                  zsiborn = trb(ji,jj,K,jpsil)**3
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

      !  Mixed-layer effect on production 
      !  Sea-ice effect on production

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zprbio(ji,jj,jk) = zprbio(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
               zprdia(ji,jj,jk) = zprdia(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
            END DO
         END DO
      END DO

      ! Computation of the various production terms 
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production terms for nanophyto. (C)
                  zprorcan(ji,jj,jk) = zprbio(ji,jj,jk)  * xlimphy(ji,jj,jk) * trb(ji,jj,K,jpphy) * rfact2
                  zpronewn(ji,jj,jk)  = zprorcan(ji,jj,jk)* xnanono3(ji,jj,jk) / ( xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk) + rtrn )
                  !
                  zratio = trb(ji,jj,K,jpnfe)   &
                  &       /(trb(ji,jj,K,jpphy) * fecnm + rtrn)
                  zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) ) 
                  zprofen(ji,jj,jk) = fecnm * zprmaxn(ji,jj,jk) * ( 1.0 - fr_i(ji,jj) )  &
                  &             * ( 4. - 4.5 * xlimnfe(ji,jj,jk) / ( xlimnfe(ji,jj,jk) + 0.5 ) )    &
                  &             * biron(ji,jj,jk) / ( biron(ji,jj,jk) + concnfe(ji,jj,jk) )  &
                  &             * zmax * trb(ji,jj,K,jpphy) * rfact2
                  !  production terms for diatoms (C)
                  zprorcad(ji,jj,jk) = zprdia(ji,jj,jk) * xlimdia(ji,jj,jk) * trb(ji,jj,K,jpdia) * rfact2
                  zpronewd(ji,jj,jk) = zprorcad(ji,jj,jk) * xdiatno3(ji,jj,jk) / ( xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk) + rtrn )
                  !
                  zratio = trb(ji,jj,K,jpdfe) &
                  &      / ( trb(ji,jj,K,jpdia) * fecdm + rtrn )
                  zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) ) 
                  zprofed(ji,jj,jk) = fecdm * zprmaxd(ji,jj,jk) * ( 1.0 - fr_i(ji,jj) )  &
                  &             * ( 4. - 4.5 * xlimdfe(ji,jj,jk) / ( xlimdfe(ji,jj,jk) + 0.5 ) )    &
                  &             * biron(ji,jj,jk) / ( biron(ji,jj,jk) + concdfe(ji,jj,jk) )  &
                  &             * zmax * trb(ji,jj,K,jpdia) * rfact2
               ENDIF
            END DO
         END DO
      END DO

      ! Computation of the chlorophyll production terms
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production terms for nanophyto. ( chlorophyll )
                  znanotot = enanom(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprod    = rday * zprorcan(ji,jj,jk) * zprnch(ji,jj,jk) * xlimphy(ji,jj,jk)
                  zprochln = chlcmin * 12. * zprorcan (ji,jj,jk)
                  chlcnm_n   = MIN ( chlcnm, ( chlcnm / (1. - 1.14 / 43.4 *tsn(ji,jj,K,jp_tem))) * (1. - 1.14 / 43.4 * 20.))
                  zprochln = zprochln + (chlcnm_n-chlcmin) * 12. * zprod / &
                                        & (  zpislopeadn(ji,jj,jk) * znanotot +rtrn)
                  !  production terms for diatoms ( chlorophyll )
                  zdiattot = ediatm(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprod    = rday * zprorcad(ji,jj,jk) * zprdch(ji,jj,jk) * xlimdia(ji,jj,jk)
                  zprochld = chlcmin * 12. * zprorcad(ji,jj,jk)
                  chlcdm_n   = MIN ( chlcdm, ( chlcdm / (1. - 1.14 / 43.4 * tsn(ji,jj,K,jp_tem))) * (1. - 1.14 / 43.4 * 20.))
                  zprochld = zprochld + (chlcdm_n-chlcmin) * 12. * zprod / &
                                        & ( zpislopeadd(ji,jj,jk) * zdiattot +rtrn )
                  !   Update the arrays TRA which contain the Chla sources and sinks
                  tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) + zprochln * texcretn
                  tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) + zprochld * texcretd
               ENDIF
            END DO
         END DO
      END DO

      !   Update the arrays TRA which contain the biological sources and sinks
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
              IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                 zproreg  = zprorcan(ji,jj,jk) - zpronewn(ji,jj,jk)
                 zproreg2 = zprorcad(ji,jj,jk) - zpronewd(ji,jj,jk)
                 zdocprod = excretd * zprorcad(ji,jj,jk) + excretn * zprorcan(ji,jj,jk)
                 tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) - zprorcan(ji,jj,jk) - zprorcad(ji,jj,jk)
                 tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - zpronewn(ji,jj,jk) - zpronewd(ji,jj,jk)
                 tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zproreg - zproreg2
                 tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zprorcan(ji,jj,jk) * texcretn
                 tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) + zprofen(ji,jj,jk) * texcretn
                 tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) + zprorcad(ji,jj,jk) * texcretd
                 tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) + zprofed(ji,jj,jk) * texcretd
                 tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zprorcad(ji,jj,jk) * zysopt(ji,jj,jk) * texcretd
                 tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zdocprod
                 tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + o2ut * ( zproreg + zproreg2) &
                 &                   + ( o2ut + o2nit ) * ( zpronewn(ji,jj,jk) + zpronewd(ji,jj,jk) )
                 !
                 zfeup = texcretn * zprofen(ji,jj,jk) + texcretd * zprofed(ji,jj,jk)
                 tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zfeup
                 tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) - texcretd * zprorcad(ji,jj,jk) * zysopt(ji,jj,jk)
                 tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprorcan(ji,jj,jk) - zprorcad(ji,jj,jk)
                 tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * ( zpronewn(ji,jj,jk) + zpronewd(ji,jj,jk) ) &
                 &                                         - rno3 * ( zproreg + zproreg2 )
              ENDIF
           END DO
        END DO
     END DO
     !
     IF( ln_ligand ) THEN
         zpligprod1(:,:,:) = 0.0    ;    zpligprod2(:,:,:) = 0.0
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                 IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                    zdocprod = excretd * zprorcad(ji,jj,jk) + excretn * zprorcan(ji,jj,jk)
                    zfeup    = texcretn * zprofen(ji,jj,jk) + texcretd * zprofed(ji,jj,jk)
                    tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) + zdocprod * ldocp   &
                    &       - zfeup * plig(ji,jj,jk) / ( rtrn + plig(ji,jj,jk) + 2.E3 * (1.0 - plig(ji,jj,jk) ) )
                    zpligprod1(ji,jj,jk) = zdocprod * ldocp
                    zpligprod2(ji,jj,jk) = zfeup * plig(ji,jj,jk) * lthet
                 ENDIF
              END DO
           END DO
        END DO
     ENDIF

#if defined key_iomput
     IF( lk_iomput ) THEN
       IF( knt == nrdttrc ) THEN
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "PPPHYN" ) .OR. iom_use( "PPPHYD" ) )  THEN
              DO jk = KRANGE
                 zw3d(:,:,jk) = zprorcan(:,:,jk) * zfact * tmask(:,:,jk)  ! primary production by nanophyto
              END DO
              CALL iom_put( "PPPHYN"  , zw3d )
              !
              DO jk = KRANGE
                 zw3d(:,:,jk) = zprorcad(:,:,jk) * zfact * tmask(:,:,jk)  ! primary production by diatomes
              END DO
              CALL iom_put( "PPPHYD"  , zw3d )
          ENDIF
          IF( iom_use( "PPNEWN" ) .OR. iom_use( "PPNEWD" ) )  THEN
              DO jk = KRANGE
                 zw3d(:,:,jk) = zpronewn(:,:,jk) * zfact * tmask(:,:,jk)  ! new primary production by nanophyto
              END DO
              CALL iom_put( "PPNEWN"  , zw3d )
              !
              DO jk = KRANGE
                 zw3d(:,:,jk) = zpronewd(:,:,jk) * zfact * tmask(:,:,jk)  ! new primary production by diatomes
              END DO
              CALL iom_put( "PPNEWD"  , zw3d )
          ENDIF
          IF( iom_use( "PBSi" ) )  THEN
              DO jk = KRANGE
                 zw3d(:,:,jk) = zprorcad(:,:,jk) * zfact * tmask(:,:,jk) * zysopt(:,:,jk) ! biogenic silica production
              END DO
              CALL iom_put( "PBSi"  , zw3d )
          ENDIF
          IF( iom_use( "PFeN" ) .OR. iom_use( "PFeD" ) )  THEN
              DO jk = KRANGE
                 zw3d(:,:,jk) = zprofen(:,:,jk) * zfact * tmask(:,:,jk)  ! biogenic iron production by nanophyto
                 END DO
              CALL iom_put( "PFeN"  , zw3d )
              !
              DO jk = KRANGE
                 zw3d(:,:,jk) = zprofed(:,:,jk) * zfact * tmask(:,:,jk)  ! biogenic iron production by  diatomes
              END DO
              CALL iom_put( "PFeD"  , zw3d )
          ENDIF
          IF( iom_use( "LPRODP" ) )  THEN
              DO jk = KRANGE
                 zw3d(:,:,jk) = zpligprod1(:,:,jk) * 1e9 * zfact * tmask(:,:,jk)
              END DO
              CALL iom_put( "LPRODP"  , zw3d )
          ENDIF
          IF( iom_use( "LDETP" ) )  THEN
              DO jk = KRANGE
                 zw3d(:,:,jk) = zpligprod2(:,:,jk) * 1e9 * zfact * tmask(:,:,jk)
              END DO
              CALL iom_put( "LDETP"  , zw3d )
          ENDIF
          IF( iom_use( "Mumax" ) )  THEN
              DO jk = KRANGE
                 zw3d(:,:,jk) = zprmaxn(:,:,jk) * tmask(:,:,jk)   ! Maximum growth rate
              END DO
              CALL iom_put( "Mumax"  , zw3d )
          ENDIF
          IF( iom_use( "MuN" ) .OR. iom_use( "MuD" ) )  THEN
              DO jk = KRANGE
                 zw3d(:,:,jk) = zprbio(:,:,jk) * xlimphy(:,:,jk) * tmask(:,:,jk)  ! Realized growth rate for nanophyto
              END DO
              CALL iom_put( "MuN"  , zw3d )
              !
              DO jk = KRANGE
                 zw3d(:,:,jk) =  zprdia(:,:,jk) * xlimdia(:,:,jk) * tmask(:,:,jk)  ! Realized growth rate for diatoms
              END DO
              CALL iom_put( "MuD"  , zw3d )
          ENDIF
          IF( iom_use( "LNlight" ) .OR. iom_use( "LDlight" ) )  THEN
              DO jk = KRANGE
                 zw3d(:,:,jk) = zprbio (:,:,jk) / (zprmaxn(:,:,jk) + rtrn) * tmask(:,:,jk) ! light limitation term
              END DO
              CALL iom_put( "LNlight"  , zw3d )
              !
              DO jk = KRANGE
                 zw3d(:,:,jk) = zprdia (:,:,jk) / (zprmaxd(:,:,jk) + rtrn) * tmask(:,:,jk)  ! light limitation term
              END DO
              CALL iom_put( "LDlight"  , zw3d )
          ENDIF
          IF( iom_use( "TPP" ) )  THEN
              DO jk = KRANGE
                 zw3d(:,:,jk) = ( zprorcan(:,:,jk) + zprorcad(:,:,jk) ) * zfact * tmask(:,:,jk)  ! total primary production
              END DO
              CALL iom_put( "TPP"  , zw3d )
          ENDIF
          IF( iom_use( "TPNEW" ) )  THEN
              DO jk = KRANGE
                 zw3d(:,:,jk) = ( zpronewn(:,:,jk) + zpronewd(:,:,jk) ) * zfact * tmask(:,:,jk)  ! total new production
              END DO
              CALL iom_put( "TPNEW"  , zw3d )
          ENDIF
          IF( iom_use( "TPBFE" ) )  THEN
              DO jk = KRANGE
                 zw3d(:,:,jk) = ( zprofen(:,:,jk) + zprofed(:,:,jk) ) * zfact * tmask(:,:,jk)  ! total biogenic iron production
              END DO
              CALL iom_put( "TPBFE"  , zw3d )
          ENDIF
          IF( iom_use( "INTPPPHYN" ) .OR. iom_use( "INTPPPHYD" ) ) THEN  
             zw2d(:,:) = 0.
             DO jk = KRANGE
               zw2d(:,:) = zw2d(:,:) + zprorcan(:,:,jk) * e3t_n(:,:,K) * zfact * tmask(:,:,jk)  ! vert. integrated  primary produc. by nano
             ENDDO
             CALL iom_put( "INTPPPHYN" , zw2d )
             !
             zw2d(:,:) = 0.
             DO jk = KRANGE
                zw2d(:,:) = zw2d(:,:) + zprorcad(:,:,jk) * e3t_n(:,:,K) * zfact * tmask(:,:,jk) ! vert. integrated  primary produc. by diatom
             ENDDO
             CALL iom_put( "INTPPPHYD" , zw2d )
          ENDIF
          IF( iom_use( "INTPP" ) ) THEN   
             zw2d(:,:) = 0.
             DO jk = KRANGE
                zw2d(:,:) = zw2d(:,:) + ( zprorcan(:,:,jk) + zprorcad(:,:,jk) ) * e3t_n(:,:,K) * zfact * tmask(:,:,jk) ! vert. integrated pp
             ENDDO
             CALL iom_put( "INTPP" , zw2d )
          ENDIF
          IF( iom_use( "INTPNEW" ) ) THEN    
             zw2d(:,:) = 0.
             DO jk = KRANGE
                zw2d(:,:) = zw2d(:,:) + ( zpronewn(:,:,jk) + zpronewd(:,:,jk) ) * e3t_n(:,:,K) * zfact * tmask(:,:,jk)  ! vert. integrated new prod
             ENDDO
             CALL iom_put( "INTPNEW" , zw2d )
          ENDIF
          IF( iom_use( "INTPBFE" ) ) THEN           !   total biogenic iron production  ( vertically integrated )
             zw2d(:,:) = 0.
             DO jk = KRANGE
                zw2d(:,:) = zw2d(:,:) + ( zprofen(:,:,jk) + zprofed(:,:,jk) ) * e3t_n(:,:,K) * zfact * tmask(:,:,jk) ! vert integr. bfe prod
             ENDDO
            CALL iom_put( "INTPBFE" , zw2d )
          ENDIF
          IF( iom_use( "INTPBSI" ) ) THEN           !   total biogenic silica production  ( vertically integrated )
             zw2d(:,:) = 0.
             DO jk = KRANGE
                zw2d(:,:) = zw2d(:,:) + zprorcad(:,:,jk) * zysopt(:,:,jk) * e3t_n(:,:,K) * zfact * tmask(:,:,jk)  ! vert integr. bsi prod
             ENDDO
             CALL iom_put( "INTPBSI" , zw2d )
          ENDIF
          !
       ENDIF
     ENDIF
#endif

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
                &                      * ( zpronewn(ji,jj,jk) + zpronewd(ji,jj,jk) ) * zmsk
             trc3d(ji,jj,K,jp_prego2)  = o2ut * &        ! Oxygen production by the Regen Produc.
                &                 (   zprorcan(ji,jj,jk) - zpronewn(ji,jj,jk)  &
                &                  + zprorcad(ji,jj,jk) - zpronewd(ji,jj,jk) ) * zmsk
         END DO
        END DO
      END DO
#endif

     IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('prod')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( charout, ltra='tra')
!         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
     ENDIF
      !
   END SUBROUTINE p4z_prod


   SUBROUTINE p4z_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_prod_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton production parameters
      !!
      !! ** Method  :   Read the nampisprod namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampisprod
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp4zprod/ pislopen, pisloped, xadap, bresp, excretn, excretd,  &
         &                 chlcnm, chlcdm, chlcmin, fecnm, fecdm, grosip
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_prod_init : phytoplankton growth'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampisprod in reference namelist : Pisces phytoplankton production
      READ  ( numnatp_ref, namp4zprod, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp4zprod in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampisprod in configuration namelist : Pisces phytoplankton production
      READ  ( numnatp_cfg, namp4zprod, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp4zprod in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, namp4zprod )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp4zprod'
         WRITE(numout,*) '      mean Si/C ratio                           grosip       =', grosip
         WRITE(numout,*) '      P-I slope                                 pislopen     =', pislopen
         WRITE(numout,*) '      Acclimation factor to low light           xadap        =', xadap
         WRITE(numout,*) '      excretion ratio of nanophytoplankton      excretn      =', excretn
         WRITE(numout,*) '      excretion ratio of diatoms                excretd      =', excretd
         WRITE(numout,*) '      basal respiration in phytoplankton        bresp        =', bresp
         WRITE(numout,*) '      Maximum Chl/C in phytoplankton            chlcmin      =', chlcmin
         WRITE(numout,*) '      P-I slope  for diatoms                    pisloped     =', pisloped
         WRITE(numout,*) '      Minimum Chl/C in nanophytoplankton        chlcnm       =', chlcnm
         WRITE(numout,*) '      Minimum Chl/C in diatoms                  chlcdm       =', chlcdm
         WRITE(numout,*) '      Maximum Fe/C in nanophytoplankton         fecnm        =', fecnm
         WRITE(numout,*) '      Minimum Fe/C in diatoms                   fecdm        =', fecdm
      ENDIF
      !
      r1_rday   = 1. / rday 
      texcretn  = 1. - excretn
      texcretd  = 1. - excretd
      !
   END SUBROUTINE p4z_prod_init


   INTEGER FUNCTION p4z_prod_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( quotan(PRIV_3D_BIOARRAY), quotad(PRIV_3D_BIOARRAY), STAT = p4z_prod_alloc )
      !
      IF( p4z_prod_alloc /= 0 ) CALL ctl_warn( 'p4z_prod_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p4z_prod_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_prod                    ! Empty routine
   END SUBROUTINE p4z_prod
#endif


   !!======================================================================
END MODULE p4zprod
