#include "cppdefs.h"

MODULE p4zopt
   !!======================================================================
   !!                         ***  MODULE p4zopt  ***
   !! TOP - PISCES : Compute the light availability in the water column
   !!======================================================================
   !! History :  1.0  !  2004     (O. Aumont) Original code
   !!            2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!            3.2  !  2009-04  (C. Ethe, G. Madec)  optimisation
   !!            3.4  !  2011-06  (O. Aumont, C. Ethe) Improve light availability of nano & diat
   !!----------------------------------------------------------------------
#if defined  key_pisces
   !!   p4z_opt       : light availability in the water column
   !!----------------------------------------------------------------------
   USE sms_pisces     ! Source Minus Sink of PISCES

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_opt        ! called in p4zbio.F90 module
   PUBLIC   p4z_opt_init   ! called in trcsms_pisces.F90 module
   PUBLIC   p4z_opt_alloc

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"


   !! * Shared module variables

   REAL(wp) ::   parlux      ! Fraction of shortwave as PAR
   REAL(wp) ::   xparsw      ! parlux/3
   REAL(wp) ::   xsi0r       ! 1. /rn_si0

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ekb, ekg, ekr  ! wavelength (Red-Green-Blue)

   INTEGER  ::   nksrp   ! levels below which the light cannot penetrate ( depth larger than 391 m)

   REAL(wp), DIMENSION(3,61) ::   xkrgb   ! tabulated attenuation coefficients for RGB absorption
   
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zopt.F90 10522 2019-01-16 08:35:15Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_opt( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_opt  ***
      !!
      !! ** Purpose :   Compute the light availability in the water column
      !!              depending on the depth and the chlorophyll concentration
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      INTEGER  ::   irgb
      REAL(wp) ::   zchl
      REAL(wp) ::   zc0 , zc1 , zc2, zc3, z1_dep
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zetmp5
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY) :: zdepmoy, zetmp1, zetmp2, zetmp3, zetmp4
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY) :: zqsr100, zqsr_corr
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zpar, ze0, ze1, ze2, ze3, zchl3d
      !!---------------------------------------------------------------------
      !
      IF( ln_p5z    )   ALLOCATE( zetmp5(PRIV_2D_BIOARRAY) )

      !     Initialisation of variables used to compute PAR
      !     -----------------------------------------------
      ze1(:,:,:) = 0.
      ze2(:,:,:) = 0.
      ze3(:,:,:) = 0.
      !
      !                                        !* attenuation coef. function of Chlorophyll and wavelength (Red-Green-Blue)
      !                                        !  --------------------------------------------------------
      DO jk = KRANGE   
         DO jj = JRANGE
            DO ji = IRANGE
                     zchl3d(ji,jj,jk) = trb(ji,jj,K,jpnch) + trb(ji,jj,K,jpdch)
                     IF( ln_p5z )   zchl3d(ji,jj,jk) = zchl3d(ji,jj,jk)    + trb(ji,jj,K,jppch)
            END DO
         END DO
      END DO

      !
      DO jk = KRANGE   
         DO jj = JRANGE
            DO ji = IRANGE
               zchl = ( zchl3d(ji,jj,jk) + rtrn ) * 1.e6
               zchl = MIN(  10. , MAX( 0.05, zchl )  )
               irgb = NINT( 41 + 20.* LOG10( zchl ) + rtrn )
               !                                                         
               ekb(ji,jj,jk) = xkrgb(1,irgb) * e3t_n(ji,jj,K)
               ekg(ji,jj,jk) = xkrgb(2,irgb) * e3t_n(ji,jj,K)
               ekr(ji,jj,jk) = xkrgb(3,irgb) * e3t_n(ji,jj,K)
            END DO
         END DO
      END DO
      !                                        !* Photosynthetically Available Radiation (PAR)
      !                                        !  --------------------------------------
      IF( l_trcdm2dc ) THEN                     !  diurnal cycle
         !
         DO jj = JRANGE
            DO ji = IRANGE
               zqsr_corr(ji,jj) = qsr_mean(ji,jj) / ( 1.-fr_i(ji,jj) + rtrn )
            END DO
         END DO
         !
         CALL p4z_opt_par( kt, zqsr_corr, ze1, ze2, ze3, pqsr100 = zqsr100 ) 
         !
         DO jk = KRANGE
            etot_ndcy(:,:,jk) =        ze1(:,:,jk) +        ze2(:,:,jk) +       ze3(:,:,jk)
            enano    (:,:,jk) =  1.85 * ze1(:,:,jk) + 0.69 * ze2(:,:,jk) + 0.46 * ze3(:,:,jk)
            ediat    (:,:,jk) =  1.62 * ze1(:,:,jk) + 0.74 * ze2(:,:,jk) + 0.63 * ze3(:,:,jk)
         END DO
         IF( ln_p5z ) THEN
            DO jk = KRANGE
              epico  (:,:,jk) =  1.94 * ze1(:,:,jk) + 0.66 * ze2(:,:,jk) + 0.4 * ze3(:,:,jk)
            END DO
         ENDIF
         !
         DO jj = JRANGE
            DO ji = IRANGE
               zqsr_corr(ji,jj) = qsr(ji,jj) / ( 1.-fr_i(ji,jj) + rtrn )
            END DO
         END DO
         !
         CALL p4z_opt_par( kt, zqsr_corr, ze1, ze2, ze3 ) 
         !
         DO jk = KRANGE
            etot(:,:,jk) =  ze1(:,:,jk) + ze2(:,:,jk) + ze3(:,:,jk)
         END DO
         !
      ELSE
         !
         DO jj = JRANGE
            DO ji = IRANGE
               zqsr_corr(ji,jj) = qsr(ji,jj) / ( 1.-fr_i(ji,jj) + rtrn )
            END DO
         END DO
         !
         CALL p4z_opt_par( kt, zqsr_corr, ze1, ze2, ze3, pqsr100 = zqsr100  ) 
         !
         DO jk = KRANGE
            etot (:,:,jk) =        ze1(:,:,jk) +        ze2(:,:,jk) +       ze3(:,:,jk)
            enano(:,:,jk) =  1.85 * ze1(:,:,jk) + 0.69 * ze2(:,:,jk) + 0.46 * ze3(:,:,jk)
            ediat(:,:,jk) =  1.62 * ze1(:,:,jk) + 0.74 * ze2(:,:,jk) + 0.63 * ze3(:,:,jk)
         END DO
         IF( ln_p5z ) THEN
            DO jk = KRANGE
              epico(:,:,jk) =  1.94 * ze1(:,:,jk) + 0.66 * ze2(:,:,jk) + 0.4 * ze3(:,:,jk)
            END DO
         ENDIF
         etot_ndcy(:,:,:) =  etot(:,:,:) 
      ENDIF


      IF( ln_qsr_bio ) THEN                    !* heat flux accros w-level (used in the dynamics)
         !                                     !  ------------------------
         DO jj = JRANGE
            DO ji = IRANGE
               zqsr_corr(ji,jj) = qsr(ji,jj)
            END DO
         END DO
         CALL p4z_opt_par( kt, zqsr_corr, ze1, ze2, ze3, pe0=ze0 )
         !
         etot3(:,:,1) =  qsr(:,:) * tmask(:,:,1)
         DO jk = KRANGE
            etot3(:,:,jk) =  ( ze0(:,:,jk) + ze1(:,:,jk) + ze2(:,:,jk) + ze3(:,:,jk) ) * tmask(:,:,jk)
         END DO
         !                                     !  ------------------------
      ENDIF
      !                                        !* Euphotic depth and level
      neln   (:,:) = 1                            !  ------------------------
      DO jj = JRANGE
         DO ji = IRANGE
            heup   (ji,jj) = gdepw_n(ji,jj,2)
            heup_01(ji,jj) = gdepw_n(ji,jj,2)
         END DO
      END DO

      DO jk = KRANGEL
         DO jj = JRANGE
           DO ji = IRANGE
              IF( etot_ndcy(ji,jj,jk) * tmask(ji,jj,jk) >=  zqsr100(ji,jj) )  THEN
                 neln(ji,jj) = jk+1                    ! Euphotic level : 1rst T-level strictly below Euphotic layer
                 !                                     ! nb: ensure the compatibility with nmld_trc definition in trd_mld_trc_zint
                 heup(ji,jj) = gdepw_n(ji,jj,jk+1)     ! Euphotic layer depth
              ENDIF
              IF( etot_ndcy(ji,jj,jk) * tmask(ji,jj,jk) >= 0.50 )  THEN
                 heup_01(ji,jj) = gdepw_n(ji,jj,jk+1)  ! Euphotic layer depth (light level definition)
              ENDIF
           END DO
        END DO
      END DO
      !
      heup   (:,:) = MIN( 300., heup   (:,:) )
      heup_01(:,:) = MIN( 300., heup_01(:,:) )
      !                                        !* mean light over the mixed layer
      zdepmoy(:,:)   = 0.e0                    !  -------------------------------
      zetmp1 (:,:)   = 0.e0
      zetmp2 (:,:)   = 0.e0

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( gdepw_n(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
                  zetmp1 (ji,jj) = zetmp1 (ji,jj) + etot     (ji,jj,jk) * e3t_n(ji,jj,K) ! remineralisation
                  zetmp2 (ji,jj) = zetmp2 (ji,jj) + etot_ndcy(ji,jj,jk) * e3t_n(ji,jj,K) ! production
                  zdepmoy(ji,jj) = zdepmoy(ji,jj) +                       e3t_n(ji,jj,K)
               ENDIF
            END DO
         END DO
      END DO
      !
      emoy(:,:,:) = etot(:,:,:)       ! remineralisation
      zpar(:,:,:) = etot_ndcy(:,:,:)  ! diagnostic : PAR with no diurnal cycle 
      !
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( gdepw_n(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
                  z1_dep = 1. / ( zdepmoy(ji,jj) + rtrn )
                  emoy (ji,jj,jk) = zetmp1(ji,jj) * z1_dep
                  zpar (ji,jj,jk) = zetmp2(ji,jj) * z1_dep
               ENDIF
            END DO
         END DO
      END DO
      !
      zdepmoy(:,:)   = 0.e0
      zetmp3 (:,:)   = 0.e0
      zetmp4 (:,:)   = 0.e0
      !
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( gdepw_n(ji,jj,jk+1) <= MIN(hmld(ji,jj), heup_01(ji,jj)) ) THEN
                  zetmp3 (ji,jj) = zetmp3 (ji,jj) + enano    (ji,jj,jk) * e3t_n(ji,jj,K) ! production
                  zetmp4 (ji,jj) = zetmp4 (ji,jj) + ediat    (ji,jj,jk) * e3t_n(ji,jj,K) ! production
                  zdepmoy(ji,jj) = zdepmoy(ji,jj) +                       e3t_n(ji,jj,K)
               ENDIF
            END DO
         END DO
      END DO
      enanom(:,:,:) = enano(:,:,:)
      ediatm(:,:,:) = ediat(:,:,:)
      !
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( gdepw_n(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
                  z1_dep = 1. / ( zdepmoy(ji,jj) + rtrn )
                  enanom(ji,jj,jk) = zetmp3(ji,jj) * z1_dep
                  ediatm(ji,jj,jk) = zetmp4(ji,jj) * z1_dep
               ENDIF
            END DO
         END DO
      END DO
      !
      IF( ln_p5z ) THEN
         zetmp5 (:,:) = 0.e0
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  IF( gdepw_n(ji,jj,jk+1) <= MIN(hmld(ji,jj), heup_01(ji,jj)) ) THEN
                     zetmp5(ji,jj)  = zetmp5 (ji,jj) + epico(ji,jj,jk) * e3t_n(ji,jj,K) ! production
                  ENDIF
               END DO
            END DO
         END DO
         !
         epicom(:,:,:) = epico(:,:,:)
         !
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  IF( gdepw_n(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
                     z1_dep = 1. / ( zdepmoy(ji,jj) + rtrn )
                     epicom(ji,jj,jk) = zetmp5(ji,jj) * z1_dep
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF

#if defined key_iomput
      IF( lk_iomput ) THEN
        IF( knt == nrdttrc ) THEN
           IF( iom_use( "Heup"  ) ) CALL iom_put( "Heup" , heup(:,:  ) * tmask(:,:,1) )  ! euphotic layer deptht
           IF( iom_use( "PARDM" ) ) CALL iom_put( "PARDM", zpar(:,:,:) * tmask(:,:,:) )  ! Photosynthetically Available Radiation
           IF( iom_use( "PAR"   ) ) CALL iom_put( "PAR"  , emoy(:,:,:) * tmask(:,:,:) )  ! Photosynthetically Available Radiation
        ENDIF
      ENDIF
#endif
      !
#if defined key_trc_diaadd
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               trc3d(ji,jj,K,jp_etot) = etot(ji,jj,jk) * tmask(ji,jj,jk)   ! PAR
            END DO
         END DO
      END DO
      !
      DO jj = JRANGE
         DO ji = IRANGE
            trc2d(ji,jj,jp_heup) = heup(ji,jj) * tmask(ji,jj,KSURF)   ! euphotic layer
         END DO
      END DO
#endif
      !
      IF( ln_p5z    )   DEALLOCATE( zetmp5 )
      !
   END SUBROUTINE p4z_opt


   SUBROUTINE p4z_opt_par( kt, pqsr, pe1, pe2, pe3, pe0, pqsr100 ) 
      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_opt_par  ***
      !!
      !! ** purpose :   compute PAR of each wavelength (Red-Green-Blue)
      !!                for a given shortwave radiation
      !!
      !!----------------------------------------------------------------------
#ifdef AGRIF
      USE ocean2pisces
#endif
      INTEGER                         , INTENT(in)              ::   kt                ! ocean time-step
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY)    , INTENT(in   )           ::   pqsr              ! shortwave
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY), INTENT(inout)           ::   pe1 , pe2 , pe3   ! PAR ( R-G-B)
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY), INTENT(inout), OPTIONAL ::   pe0               !
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY)    , INTENT(  out), OPTIONAL ::   pqsr100           !
      !
      INTEGER    ::   ji, jj, jk     ! dummy loop indices
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY) ::  zqsr   ! shortwave
      !!----------------------------------------------------------------------

      !  Real shortwave
      zqsr(:,:) = xparsw         * pqsr(:,:)
      
      !  Light at the euphotic depth 
      IF( PRESENT( pqsr100 ) )   pqsr100(:,:) = 0.01 * 3. * zqsr(:,:)

      IF( PRESENT( pe0 ) ) THEN     !  W-level
         !
         pe0(:,:,1) = pqsr(:,:) - 3. * zqsr(:,:)    !   ( 1 - 3 * alpha ) * q
         pe1(:,:,1) = zqsr(:,:)         
         pe2(:,:,1) = zqsr(:,:)
         pe3(:,:,1) = zqsr(:,:)
         !
         DO jk = KRANGEL
            DO jj = JRANGE
               DO ji = IRANGE
                  pe0(ji,jj,jk) = pe0(ji,jj,jk-1) * EXP( -e3t_n(ji,jj,KUP) * xsi0r )
                  pe1(ji,jj,jk) = pe1(ji,jj,jk-1) * EXP( -ekb  (ji,jj,jk-1 )        )
                  pe2(ji,jj,jk) = pe2(ji,jj,jk-1) * EXP( -ekg  (ji,jj,jk-1 )        )
                  pe3(ji,jj,jk) = pe3(ji,jj,jk-1) * EXP( -ekr  (ji,jj,jk-1 )        )
               END DO
              !
            END DO
            !
         END DO
        !
      ELSE   ! T- level
        !
        pe1(:,:,1) = zqsr(:,:) * EXP( -0.5 * ekb(:,:,1) )
        pe2(:,:,1) = zqsr(:,:) * EXP( -0.5 * ekg(:,:,1) )
        pe3(:,:,1) = zqsr(:,:) * EXP( -0.5 * ekr(:,:,1) )
        !
        DO jk = KRANGEL
           DO jj = JRANGE
              DO ji = IRANGE
                 pe1(ji,jj,jk) = pe1(ji,jj,jk-1) * EXP( -0.5 * ( ekb(ji,jj,jk-1) + ekb(ji,jj,jk) ) )
                 pe2(ji,jj,jk) = pe2(ji,jj,jk-1) * EXP( -0.5 * ( ekg(ji,jj,jk-1) + ekg(ji,jj,jk) ) )
                 pe3(ji,jj,jk) = pe3(ji,jj,jk-1) * EXP( -0.5 * ( ekr(ji,jj,jk-1) + ekr(ji,jj,jk) ) )
              END DO
           END DO
        END DO    
        !
      ENDIF
      ! 
   END SUBROUTINE p4z_opt_par

   SUBROUTINE p4z_opt_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_opt_init  ***
      !!
      !! ** Purpose :   Initialization of tabulated attenuation coef
      !!                and of the percentage of PAR in Shortwave
      !!
      !! ** Input   :   external ascii and netcdf files
      !!----------------------------------------------------------------------
      INTEGER :: numpar, ierr, ios   ! Local integer 
      !
      NAMELIST/nampisopt/ parlux
      !!----------------------------------------------------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_opt_init : '
         WRITE(numout,*) '~~~~~~~~~~~~ '
      ENDIF
      REWIND( numnatp_ref )              ! Namelist nampisopt in reference namelist : Pisces attenuation coef. and PAR
      READ  ( numnatp_ref, nampisopt, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampisopt in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampisopt in configuration namelist : Pisces attenuation coef. and PAR
      READ  ( numnatp_cfg, nampisopt, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nampisopt in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisopt )

      IF(lwp) THEN
         WRITE(numout,*) '   Namelist : nampisopt '
         WRITE(numout,*) '      Default value for the PAR fraction   parlux         = ', parlux
      ENDIF
      !
      xparsw = parlux / 3.0
      xsi0r  = 1.e0 / rn_si0
      !
!      CALL trc_oce_rgb( xkrgb )                  ! tabulated attenuation coefficients
      CALL trc_oce_rgb_read( xkrgb )     ! tabulated attenuation coefficients
      !
                         ekr      (:,:,:) = 0.
                         ekb      (:,:,:) = 0.
                         ekg      (:,:,:) = 0.
                         etot     (:,:,:) = 0.
                         etot_ndcy(:,:,:) = 0.
                         enano    (:,:,:) = 0.
                         ediat    (:,:,:) = 0.
      IF( ln_p5z     )   epico    (:,:,:) = 0.
      IF( ln_qsr_bio )   etot3    (:,:,:) = 0.
      ! 
   END SUBROUTINE p4z_opt_init


   INTEGER FUNCTION p4z_opt_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_opt_alloc  ***
      !!----------------------------------------------------------------------
      !
      ALLOCATE( ekb(PRIV_3D_BIOARRAY), ekr(PRIV_3D_BIOARRAY),  &
                ekg(PRIV_3D_BIOARRAY), STAT= p4z_opt_alloc  ) 
      !
      IF( p4z_opt_alloc /= 0 ) CALL ctl_warn( 'p4z_opt_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p4z_opt_alloc

   SUBROUTINE trc_oce_rgb( prgb )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_opt_init  ***
      !!
      !! ** Purpose :   Initialization of of the optical scheme
      !!
      !! ** Method  :   Set a look up table for the optical coefficients
      !!                i.e. the attenuation coefficient for R-G-B light
      !!                tabulated in Chlorophyll class (from JM Andre)
      !!
      !! ** Action  :   prgb(3,61) tabulated R-G-B attenuation coef.
      !!
      !! Reference  : Lengaigne et al. 2007, Clim. Dyn., V28, 5, 503-516.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(3,61), INTENT(out) ::   prgb   ! tabulated attenuation coefficient
      !!
      INTEGER  ::   jc     ! dummy loop indice
      INTEGER  ::   irgb   ! temporary integer
      REAL(wp) ::   zchl   ! temporary scalar
      REAL(wp), DIMENSION(4,61) ::   zrgb   ! tabulated attenuation coefficient (formerly read in 'kRGB61.txt')
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_oce_rgb : Initialisation of the optical look-up table'
         WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF
      !
      !  Chlorophyll        !     Blue attenuation     !     Green attenuation    !     Red attenuation      !
      zrgb(1, 1) =  0.010   ;   zrgb(2, 1) = 0.01618   ;   zrgb(3, 1) = 0.07464   ;   zrgb(4, 1) = 0.37807
      zrgb(1, 2) =  0.011   ;   zrgb(2, 2) = 0.01654   ;   zrgb(3, 2) = 0.07480   ;   zrgb(4, 2) = 0.37823
      zrgb(1, 3) =  0.013   ;   zrgb(2, 3) = 0.01693   ;   zrgb(3, 3) = 0.07499   ;   zrgb(4, 3) = 0.37840
      zrgb(1, 4) =  0.014   ;   zrgb(2, 4) = 0.01736   ;   zrgb(3, 4) = 0.07518   ;   zrgb(4, 4) = 0.37859
      zrgb(1, 5) =  0.016   ;   zrgb(2, 5) = 0.01782   ;   zrgb(3, 5) = 0.07539   ;   zrgb(4, 5) = 0.37879
      zrgb(1, 6) =  0.018   ;   zrgb(2, 6) = 0.01831   ;   zrgb(3, 6) = 0.07562   ;   zrgb(4, 6) = 0.37900
      zrgb(1, 7) =  0.020   ;   zrgb(2, 7) = 0.01885   ;   zrgb(3, 7) = 0.07586   ;   zrgb(4, 7) = 0.37923
      zrgb(1, 8) =  0.022   ;   zrgb(2, 8) = 0.01943   ;   zrgb(3, 8) = 0.07613   ;   zrgb(4, 8) = 0.37948
      zrgb(1, 9) =  0.025   ;   zrgb(2, 9) = 0.02005   ;   zrgb(3, 9) = 0.07641   ;   zrgb(4, 9) = 0.37976
      zrgb(1,10) =  0.028   ;   zrgb(2,10) = 0.02073   ;   zrgb(3,10) = 0.07672   ;   zrgb(4,10) = 0.38005
      zrgb(1,11) =  0.032   ;   zrgb(2,11) = 0.02146   ;   zrgb(3,11) = 0.07705   ;   zrgb(4,11) = 0.38036
      zrgb(1,12) =  0.035   ;   zrgb(2,12) = 0.02224   ;   zrgb(3,12) = 0.07741   ;   zrgb(4,12) = 0.38070
      zrgb(1,13) =  0.040   ;   zrgb(2,13) = 0.02310   ;   zrgb(3,13) = 0.07780   ;   zrgb(4,13) = 0.38107
      zrgb(1,14) =  0.045   ;   zrgb(2,14) = 0.02402   ;   zrgb(3,14) = 0.07821   ;   zrgb(4,14) = 0.38146
      zrgb(1,15) =  0.050   ;   zrgb(2,15) = 0.02501   ;   zrgb(3,15) = 0.07866   ;   zrgb(4,15) = 0.38189
      zrgb(1,16) =  0.056   ;   zrgb(2,16) = 0.02608   ;   zrgb(3,16) = 0.07914   ;   zrgb(4,16) = 0.38235
      zrgb(1,17) =  0.063   ;   zrgb(2,17) = 0.02724   ;   zrgb(3,17) = 0.07967   ;   zrgb(4,17) = 0.38285
      zrgb(1,18) =  0.071   ;   zrgb(2,18) = 0.02849   ;   zrgb(3,18) = 0.08023   ;   zrgb(4,18) = 0.38338
      zrgb(1,19) =  0.079   ;   zrgb(2,19) = 0.02984   ;   zrgb(3,19) = 0.08083   ;   zrgb(4,19) = 0.38396
      zrgb(1,20) =  0.089   ;   zrgb(2,20) = 0.03131   ;   zrgb(3,20) = 0.08149   ;   zrgb(4,20) = 0.38458
      zrgb(1,21) =  0.100   ;   zrgb(2,21) = 0.03288   ;   zrgb(3,21) = 0.08219   ;   zrgb(4,21) = 0.38526
      zrgb(1,22) =  0.112   ;   zrgb(2,22) = 0.03459   ;   zrgb(3,22) = 0.08295   ;   zrgb(4,22) = 0.38598
      zrgb(1,23) =  0.126   ;   zrgb(2,23) = 0.03643   ;   zrgb(3,23) = 0.08377   ;   zrgb(4,23) = 0.38676
      zrgb(1,24) =  0.141   ;   zrgb(2,24) = 0.03842   ;   zrgb(3,24) = 0.08466   ;   zrgb(4,24) = 0.38761
      zrgb(1,25) =  0.158   ;   zrgb(2,25) = 0.04057   ;   zrgb(3,25) = 0.08561   ;   zrgb(4,25) = 0.38852
      zrgb(1,26) =  0.178   ;   zrgb(2,26) = 0.04289   ;   zrgb(3,26) = 0.08664   ;   zrgb(4,26) = 0.38950
      zrgb(1,27) =  0.200   ;   zrgb(2,27) = 0.04540   ;   zrgb(3,27) = 0.08775   ;   zrgb(4,27) = 0.39056
      zrgb(1,28) =  0.224   ;   zrgb(2,28) = 0.04811   ;   zrgb(3,28) = 0.08894   ;   zrgb(4,28) = 0.39171
      zrgb(1,29) =  0.251   ;   zrgb(2,29) = 0.05103   ;   zrgb(3,29) = 0.09023   ;   zrgb(4,29) = 0.39294
      zrgb(1,30) =  0.282   ;   zrgb(2,30) = 0.05420   ;   zrgb(3,30) = 0.09162   ;   zrgb(4,30) = 0.39428
      zrgb(1,31) =  0.316   ;   zrgb(2,31) = 0.05761   ;   zrgb(3,31) = 0.09312   ;   zrgb(4,31) = 0.39572
      zrgb(1,32) =  0.355   ;   zrgb(2,32) = 0.06130   ;   zrgb(3,32) = 0.09474   ;   zrgb(4,32) = 0.39727
      zrgb(1,33) =  0.398   ;   zrgb(2,33) = 0.06529   ;   zrgb(3,33) = 0.09649   ;   zrgb(4,33) = 0.39894
      zrgb(1,34) =  0.447   ;   zrgb(2,34) = 0.06959   ;   zrgb(3,34) = 0.09837   ;   zrgb(4,34) = 0.40075
      zrgb(1,35) =  0.501   ;   zrgb(2,35) = 0.07424   ;   zrgb(3,35) = 0.10040   ;   zrgb(4,35) = 0.40270
      zrgb(1,36) =  0.562   ;   zrgb(2,36) = 0.07927   ;   zrgb(3,36) = 0.10259   ;   zrgb(4,36) = 0.40480
      zrgb(1,37) =  0.631   ;   zrgb(2,37) = 0.08470   ;   zrgb(3,37) = 0.10495   ;   zrgb(4,37) = 0.40707
      zrgb(1,38) =  0.708   ;   zrgb(2,38) = 0.09056   ;   zrgb(3,38) = 0.10749   ;   zrgb(4,38) = 0.40952
      zrgb(1,39) =  0.794   ;   zrgb(2,39) = 0.09690   ;   zrgb(3,39) = 0.11024   ;   zrgb(4,39) = 0.41216
      zrgb(1,40) =  0.891   ;   zrgb(2,40) = 0.10374   ;   zrgb(3,40) = 0.11320   ;   zrgb(4,40) = 0.41502
      zrgb(1,41) =  1.000   ;   zrgb(2,41) = 0.11114   ;   zrgb(3,41) = 0.11639   ;   zrgb(4,41) = 0.41809
      zrgb(1,42) =  1.122   ;   zrgb(2,42) = 0.11912   ;   zrgb(3,42) = 0.11984   ;   zrgb(4,42) = 0.42142
      zrgb(1,43) =  1.259   ;   zrgb(2,43) = 0.12775   ;   zrgb(3,43) = 0.12356   ;   zrgb(4,43) = 0.42500
      zrgb(1,44) =  1.413   ;   zrgb(2,44) = 0.13707   ;   zrgb(3,44) = 0.12757   ;   zrgb(4,44) = 0.42887
      zrgb(1,45) =  1.585   ;   zrgb(2,45) = 0.14715   ;   zrgb(3,45) = 0.13189   ;   zrgb(4,45) = 0.43304
      zrgb(1,46) =  1.778   ;   zrgb(2,46) = 0.15803   ;   zrgb(3,46) = 0.13655   ;   zrgb(4,46) = 0.43754
      zrgb(1,47) =  1.995   ;   zrgb(2,47) = 0.16978   ;   zrgb(3,47) = 0.14158   ;   zrgb(4,47) = 0.44240
      zrgb(1,48) =  2.239   ;   zrgb(2,48) = 0.18248   ;   zrgb(3,48) = 0.14701   ;   zrgb(4,48) = 0.44765
      zrgb(1,49) =  2.512   ;   zrgb(2,49) = 0.19620   ;   zrgb(3,49) = 0.15286   ;   zrgb(4,49) = 0.45331
      zrgb(1,50) =  2.818   ;   zrgb(2,50) = 0.21102   ;   zrgb(3,50) = 0.15918   ;   zrgb(4,50) = 0.45942
      zrgb(1,51) =  3.162   ;   zrgb(2,51) = 0.22703   ;   zrgb(3,51) = 0.16599   ;   zrgb(4,51) = 0.46601
      zrgb(1,52) =  3.548   ;   zrgb(2,52) = 0.24433   ;   zrgb(3,52) = 0.17334   ;   zrgb(4,52) = 0.47313
      zrgb(1,53) =  3.981   ;   zrgb(2,53) = 0.26301   ;   zrgb(3,53) = 0.18126   ;   zrgb(4,53) = 0.48080
      zrgb(1,54) =  4.467   ;   zrgb(2,54) = 0.28320   ;   zrgb(3,54) = 0.18981   ;   zrgb(4,54) = 0.48909
      zrgb(1,55) =  5.012   ;   zrgb(2,55) = 0.30502   ;   zrgb(3,55) = 0.19903   ;   zrgb(4,55) = 0.49803
      zrgb(1,56) =  5.623   ;   zrgb(2,56) = 0.32858   ;   zrgb(3,56) = 0.20898   ;   zrgb(4,56) = 0.50768
      zrgb(1,57) =  6.310   ;   zrgb(2,57) = 0.35404   ;   zrgb(3,57) = 0.21971   ;   zrgb(4,57) = 0.51810
      zrgb(1,58) =  7.079   ;   zrgb(2,58) = 0.38154   ;   zrgb(3,58) = 0.23129   ;   zrgb(4,58) = 0.52934
      zrgb(1,59) =  7.943   ;   zrgb(2,59) = 0.41125   ;   zrgb(3,59) = 0.24378   ;   zrgb(4,59) = 0.54147
      zrgb(1,60) =  8.912   ;   zrgb(2,60) = 0.44336   ;   zrgb(3,60) = 0.25725   ;   zrgb(4,60) = 0.55457
      zrgb(1,61) = 10.000   ;   zrgb(2,61) = 0.47804   ;   zrgb(3,61) = 0.27178   ;   zrgb(4,61) = 0.56870
      !
      prgb(:,:) = zrgb(2:4,:)
      !
      DO jc = 1, 61                         ! check
         zchl = zrgb(1,jc)
         irgb = NINT( 41 + 20.* LOG10( zchl ) + 1.e-15 )
         IF(lwp) WRITE(numout,*) '    jc =', jc, '  Chl = ', zchl, '  irgb = ', irgb
         IF( irgb /= jc ) THEN
            IF(lwp) WRITE(numout,*) '    jc =', jc, '  Chl = ', zchl, '  Chl class = ', irgb
            CALL ctl_stop( 'trc_oce_rgb : inconsistency in Chl tabulated attenuation coeff.' )
         ENDIF
      END DO
      !
   END SUBROUTINE trc_oce_rgb


   SUBROUTINE trc_oce_rgb_read( prgb )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_opt_init  ***
      !!
      !! ** Purpose :   Initialization of of the optical scheme
      !!
      !! ** Method  :   read the look up table for the optical coefficients
      !!
      !! ** input   :   xkrgb(61) precomputed array corresponding to the
      !!                          attenuation coefficient (from JM Andre)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(3,61), INTENT(out) ::   prgb   ! tabulated attenuation coefficient
      !!
      INTEGER  ::   jchl, jband   ! dummy loop indices
      INTEGER  ::   numlight
      REAL(wp) ::   zchl
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' trc_oce_rgb_read : optical look-up table read in kRGB61.txt file'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~'
      ENDIF
      !
      CALL ctl_opn( numlight, 'kRGB61.txt', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )
      DO jchl = 1, 61
         READ(numlight,*) zchl, ( prgb(jband,jchl), jband=1,3 )
      END DO
      CLOSE( numlight )
      !
      !
   END SUBROUTINE trc_oce_rgb_read

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                   No PISCES bio-model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE p4z_opt                   ! Empty routine
   END SUBROUTINE p4z_opt
#endif 

   !!======================================================================
END MODULE p4zopt
