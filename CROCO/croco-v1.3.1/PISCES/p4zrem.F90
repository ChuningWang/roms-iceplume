#include "cppdefs.h"

MODULE p4zrem
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP :   PISCES Compute remineralization/dissolution of organic compounds
   !!=========================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!   p4z_rem       :  Compute remineralization/dissolution of organic compounds
   !!   p4z_rem_init  :  Initialisation of parameters for remineralisation
   !!   p4z_rem_alloc :  Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zche          !  chemical model
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p4zlim
!   USE iom             !  I/O manager


   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_rem         ! called in p4zbio.F90
   PUBLIC   p4z_rem_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_rem_alloc

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"

   REAL(wp), PUBLIC ::   xremikc    !: remineralisation rate of DOC 
   REAL(wp), PUBLIC ::   xremikn    !: remineralisation rate of DON 
   REAL(wp), PUBLIC ::   xremikp    !: remineralisation rate of DOP 
   REAL(wp), PUBLIC ::   xremik     !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::   nitrif     !: NH4 nitrification rate 
   REAL(wp), PUBLIC ::   xsirem     !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::   xsiremlab  !: fast remineralisation rate of POC 
   REAL(wp), PUBLIC ::   xsilab     !: fraction of labile biogenic silica 
   REAL(wp), PUBLIC ::   feratb     !: Fe/C quota in bacteria
   REAL(wp), PUBLIC ::   xkferb     !: Half-saturation constant for bacteria Fe/C

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr   !: denitrification array

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zrem.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_rem( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic compounds
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zremik, zremikc, zremikn, zremikp, zsiremin, zfact 
      REAL(wp) ::   zsatur, zsatur2, znusil, znusil2, zdep, zdepmin, zfactdep
      REAL(wp) ::   zbactfer, zolimit, zrfact2, zmsk
      REAL(wp) ::   zammonic, zoxyremc, zoxyremn, zoxyremp
      REAL(wp) ::   zosil, ztem, zdenitnh4, zolimic, zolimin, zolimip, zdenitrn, zdenitrp
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY) :: ztempbac
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zdepbac, zolimi, zonitr, zdepprod,   &
      &                                        zfacsi, zfacsib, zdepeff, zfebact
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      !!---------------------------------------------------------------------
      !
      ! Initialisation of arrys
      zdepprod(:,:,:) = 1.
      zdepeff (:,:,:) = 0.3
      ztempbac(:,:)   = 0.
      zfacsib(:,:,:)  = xsilab / ( 1.0 - xsilab )
      zfebact(:,:,:)  = 0.
      zfacsi(:,:,:)   = xsilab

      ! Computation of the mean phytoplankton concentration as
      ! a crude estimate of the bacterial biomass
      ! this parameterization has been deduced from a model version
      ! that was modeling explicitely bacteria
      ! -------------------------------------------------------
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zdep = MAX( hmld(ji,jj), heup(ji,jj) )
               IF( gdept_n(ji,jj,K) < zdep ) THEN
                  zdepbac(ji,jj,jk) = MIN( 0.7 * ( trb(ji,jj,K,jpzoo)   &
                  &                      + 2.* trb(ji,jj,K,jpmes) ), 4.e-6 )
                  ztempbac(ji,jj)   = zdepbac(ji,jj,jk)
               ELSE
                  zdepmin = MIN( 1., zdep / gdept_n(ji,jj,K) )
                  zdepbac (ji,jj,jk) = zdepmin**0.683 * ztempbac(ji,jj)
                  zdepprod(ji,jj,jk) = zdepmin**0.273
                  zdepeff (ji,jj,jk) = zdepeff(ji,jj,jk) * zdepmin**0.3
               ENDIF
            END DO
         END DO
      END DO

      IF( ln_p4z ) THEN
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  ! DOC ammonification. Depends on depth, phytoplankton biomass
                  ! and a limitation term which is supposed to be a parameterization of the bacterial activity. 
                  zremik = xremik * xstep / 1.e-6 * xlimbac(ji,jj,jk) * zdepbac(ji,jj,jk) 
                  zremik = MAX( zremik, 2.74e-4 * xstep )
                  ! Ammonification in oxic waters with oxygen consumption
                  ! -----------------------------------------------------
                  zolimit = zremik * ( 1.- nitrfac(ji,jj,jk) )   &
                  &                * trb(ji,jj,K,jpdoc) 
                  zolimi(ji,jj,jk) = MIN( ( trb(ji,jj,K,jpoxy) - rtrn ) / o2ut, zolimit ) 
                  ! Ammonification in suboxic waters with denitrification
                  ! -------------------------------------------------------
                  zammonic = zremik * nitrfac(ji,jj,jk) * trb(ji,jj,K,jpdoc)
                  denitr(ji,jj,jk)  = zammonic * ( 1. - nitrfac2(ji,jj,jk) )
                  denitr(ji,jj,jk)  = MIN( ( trb(ji,jj,K,jpno3) - rtrn ) / rdenit, denitr(ji,jj,jk) )
                  zoxyremc          = zammonic - denitr(ji,jj,jk)
                  !
                  zolimi (ji,jj,jk) = MAX( 0.e0, zolimi (ji,jj,jk) )
                  denitr (ji,jj,jk) = MAX( 0.e0, denitr (ji,jj,jk) )
                  zoxyremc          = MAX( 0.e0, zoxyremc )

                  !
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zolimi (ji,jj,jk) + denitr(ji,jj,jk) + zoxyremc
                  tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zolimi (ji,jj,jk) + denitr(ji,jj,jk) + zoxyremc
                  tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - denitr (ji,jj,jk) * rdenit
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zolimi (ji,jj,jk) - denitr(ji,jj,jk) - zoxyremc
                  tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - zolimi (ji,jj,jk) * o2ut
                  tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zolimi (ji,jj,jk) + denitr(ji,jj,jk) + zoxyremc
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * ( zolimi(ji,jj,jk) + zoxyremc    &
                  &                     + ( rdenit + 1.) * denitr(ji,jj,jk) )
               END DO
            END DO
         END DO
      ELSE
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  ! DOC ammonification. Depends on depth, phytoplankton biomass
                  ! and a limitation term which is supposed to be a parameterization of the bacterial activity. 
                  ! -----------------------------------------------------------------
                  zremik = xstep / 1.e-6 * MAX(0.01, xlimbac(ji,jj,jk)) * zdepbac(ji,jj,jk) 
                  zremik = MAX( zremik, 2.74e-4 * xstep / xremikc )

                  zremikc = xremikc * zremik
                  zremikn = xremikn / xremikc
                  zremikp = xremikp / xremikc

                  ! Ammonification in oxic waters with oxygen consumption
                  ! -----------------------------------------------------
                  zolimit = zremikc * ( 1.- nitrfac(ji,jj,jk) ) * trb(ji,jj,K,jpdoc) 
                  zolimic = MAX( 0.e0, MIN( ( trb(ji,jj,K,jpoxy) - rtrn ) / o2ut, zolimit ) ) 
                  zolimi(ji,jj,jk) = zolimic
                  zolimin = zremikn * zolimic * trb(ji,jj,K,jpdon)   &
                  &       / ( trb(ji,jj,K,jpdoc) + rtrn )
                  zolimip = zremikp * zolimic * trb(ji,jj,K,jpdop)   &
                  &       / ( trb(ji,jj,K,jpdoc) + rtrn ) 

                  ! Ammonification in suboxic waters with denitrification
                  ! -------------------------------------------------------
                  zammonic = zremikc * nitrfac(ji,jj,jk) * trb(ji,jj,K,jpdoc)
                  denitr(ji,jj,jk)  = zammonic * ( 1. - nitrfac2(ji,jj,jk) )
                  denitr(ji,jj,jk)  = MAX(0., MIN(  ( trb(ji,jj,K,jpno3) - rtrn )   &
                  &                 / rdenit, denitr(ji,jj,jk) ) )
                  zoxyremc          = MAX(0., zammonic - denitr(ji,jj,jk))
                  zdenitrn  = zremikn * denitr(ji,jj,jk) * trb(ji,jj,K,jpdon)   &
                  &         / ( trb(ji,jj,K,jpdoc) + rtrn )
                  zdenitrp  = zremikp * denitr(ji,jj,jk) * trb(ji,jj,K,jpdop)   &
                  &         / ( trb(ji,jj,K,jpdoc) + rtrn )
                  zoxyremn  = zremikn * zoxyremc * trb(ji,jj,K,jpdon)   &
                  &         / ( trb(ji,jj,K,jpdoc) + rtrn )
                  zoxyremp  = zremikp * zoxyremc * trb(ji,jj,K,jpdop)   &
                  &         / ( trb(ji,jj,K,jpdoc) + rtrn )

                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zolimip + zdenitrp + zoxyremp
                  tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zolimin + zdenitrn + zoxyremn
                  tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - denitr(ji,jj,jk) * rdenit
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zolimic - denitr(ji,jj,jk) - zoxyremc
                  tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) - zolimin - zdenitrn - zoxyremn
                  tra(ji,jj,jk,jpdop) = tra(ji,jj,jk,jpdop) - zolimip - zdenitrp - zoxyremp
                  tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - zolimic * o2ut
                  tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zolimic + denitr(ji,jj,jk) + zoxyremc
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * ( zolimin + zoxyremn + ( rdenit + 1.) * zdenitrn )
               END DO
            END DO
         END DO
         !
      ENDIF

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               ! NH4 nitrification to NO3. Ceased for oxygen concentrations
               ! below 2 umol/L. Inhibited at strong light 
               ! ----------------------------------------------------------
               zonitr(ji,jj,jk)  = nitrif * xstep * trb(ji,jj,K,jpnh4) * ( 1.- nitrfac(ji,jj,jk) )  &
               &         / ( 1.+ emoy(ji,jj,jk) ) * ( 1. + fr_i(ji,jj) * emoy(ji,jj,jk) ) 
               zdenitnh4 = nitrif * xstep * trb(ji,jj,K,jpnh4) * nitrfac(ji,jj,jk)
               zdenitnh4 = MIN(  ( trb(ji,jj,K,jpno3) - rtrn ) / rdenita, zdenitnh4 ) 
               ! Update of the tracers trends
               ! ----------------------------
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zonitr(ji,jj,jk) - zdenitnh4
               tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + zonitr(ji,jj,jk) - rdenita * zdenitnh4
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2nit * zonitr(ji,jj,jk)
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2 * rno3 * zonitr(ji,jj,jk)+ rno3 * ( rdenita - 1. ) * zdenitnh4
            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem1')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( charout, ltra='tra')
       ENDIF

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE

               ! Bacterial uptake of iron. No iron is available in DOC. So
               ! Bacteries are obliged to take up iron from the water. Some
               ! studies (especially at Papa) have shown this uptake to be significant
               ! ----------------------------------------------------------
               zbactfer = feratb *  rfact2 * 0.6 / rday * tgfunc(ji,jj,jk) * xlimbacl(ji,jj,jk)   &
                  &              * trb(ji,jj,K,jpfer)                                             &
                  &              / ( xkferb + trb(ji,jj,K,jpfer) )                                &
                  &              * zdepprod(ji,jj,jk) * zdepeff(ji,jj,jk) * zdepbac(ji,jj,jk)
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zbactfer*0.33
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zbactfer*0.25
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zbactfer*0.08
               zfebact(ji,jj,jk)   = zbactfer * 0.33
               blim(ji,jj,jk)      = xlimbacl(ji,jj,jk)  * zdepbac(ji,jj,jk) / 1.e-6 * zdepprod(ji,jj,jk)
            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem2')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( charout, ltra='tra')
       ENDIF

      ! Initialization of the array which contains the labile fraction
      ! of bSi. Set to a constant in the upper ocean
      ! ---------------------------------------------------------------

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zdep     = MAX( hmld(ji,jj), heup_01(ji,jj) )
               zsatur   = MAX( rtrn, ( sio3eq(ji,jj,jk) - trb(ji,jj,K,jpsil) )   &
               &        / ( sio3eq(ji,jj,jk) + rtrn ) )
               zsatur2  = ( 1. + tsn(ji,jj,K,jp_tem) / 400.)**37
               znusil   = 0.225  * ( 1. + tsn(ji,jj,K,jp_tem) / 15.) * zsatur + 0.775 * zsatur2 * zsatur**9.25
               ! Remineralization rate of BSi depedant on T and saturation
               ! ---------------------------------------------------------
               IF ( gdept_n(ji,jj,K) > zdep ) THEN
                  zfacsib(ji,jj,jk) = zfacsib(ji,jj,jk-1) * EXP( -0.5 * ( xsiremlab - xsirem )  &
                  &                   * znusil * e3t_n(ji,jj,K) / wsbio4(ji,jj,jk) )
                  zfacsi(ji,jj,jk)  = zfacsib(ji,jj,jk) / ( 1.0 + zfacsib(ji,jj,jk) )
                  zfacsib(ji,jj,jk) = zfacsib(ji,jj,jk) * EXP( -0.5 * ( xsiremlab - xsirem )    &
                  &                   * znusil * e3t_n(ji,jj,K) / wsbio4(ji,jj,jk) )
               ENDIF
               zsiremin = ( xsiremlab * zfacsi(ji,jj,jk) + xsirem * ( 1. - zfacsi(ji,jj,jk) ) ) * xstep * znusil
               zosil    = zsiremin * trb(ji,jj,K,jpgsi)
               !
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) - zosil
               tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) + zosil
            END DO
         END DO
      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem3')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( charout, ltra='tra')
       ENDIF

#if defined key_iomput
     IF( lk_iomput ) THEN
      IF( knt == nrdttrc ) THEN
          zrfact2 = 1.e3 * rfact2r
          ALLOCATE( zw3d(PRIV_3D_BIOARRAY) )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "REMIN" ) )  THEN
              zw3d(:,:,:) = zolimi(:,:,:) * tmask(:,:,:) * zfact !  Remineralisation rate
              CALL iom_put( "REMIN"  , zw3d )
          ENDIF
          IF( iom_use( "DENIT" ) )  THEN
              zw3d(:,:,:) = denitr(:,:,:) * rdenit * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "DENIT"  , zw3d )
          ENDIF
          IF( iom_use( "BACT" ) )  THEN
               zw3d(:,:,:) = zdepbac(:,:,:) * 1.E6 * tmask(:,:,:)  ! Bacterial biomass
               CALL iom_put( "BACT", zw3d )
          ENDIF
          IF( iom_use( "FEBACT" ) )  THEN
               zw3d(:,:,:) = zfebact(:,:,:) * 1E9 * tmask(:,:,:) * zrfact2   ! Bacterial iron consumption
               CALL iom_put( "FEBACT" , zw3d )
          ENDIF
          !
          DEALLOCATE( zw3d )
       ENDIF
      ENDIF
#endif
      !
# if defined key_trc_diaadd
     zrfact2 = 1.e3 * rfact2r
     DO jk = KRANGE
        DO jj = JRANGE
          DO ji = IRANGE
             zmsk = zrfact2 * tmask(ji,jj,jk)
             trc3d(ji,jj,K,jp_remino2) = (-o2ut)  * zolimi(ji,jj,jk) * zmsk ! O2 consumption by nitrification
             trc3d(ji,jj,K,jp_nitrifo2)  = (-o2nit) * zonitr(ji,jj,jk) * zmsk ! O2 consumption by remin
         END DO
        END DO
      END DO
#endif

   END SUBROUTINE p4z_rem


   SUBROUTINE p4z_rem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_rem_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampisrem namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisrem
      !!
      !!----------------------------------------------------------------------
      NAMELIST/nampisrem/ xremik, nitrif, xsirem, xsiremlab, xsilab, feratb, xkferb, & 
         &                xremikc, xremikn, xremikp
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_rem_init : Initialization of remineralization parameters'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampisrem in reference namelist : Pisces remineralization
      READ  ( numnatp_ref, nampisrem, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampisrem in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampisrem in configuration namelist : Pisces remineralization
      READ  ( numnatp_cfg, nampisrem, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nampisrem in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, nampisrem )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist parameters for remineralization, nampisrem'
         IF( ln_p4z ) THEN
            WRITE(numout,*) '      remineralization rate of DOC              xremik    =', xremik
         ELSE
            WRITE(numout,*) '      remineralization rate of DOC              xremikc   =', xremikc
            WRITE(numout,*) '      remineralization rate of DON              xremikn   =', xremikn
            WRITE(numout,*) '      remineralization rate of DOP              xremikp   =', xremikp
         ENDIF
         WRITE(numout,*) '      remineralization rate of Si               xsirem    =', xsirem
         WRITE(numout,*) '      fast remineralization rate of Si          xsiremlab =', xsiremlab
         WRITE(numout,*) '      fraction of labile biogenic silica        xsilab    =', xsilab
         WRITE(numout,*) '      NH4 nitrification rate                    nitrif    =', nitrif
         WRITE(numout,*) '      Bacterial Fe/C ratio                      feratb    =', feratb
         WRITE(numout,*) '      Half-saturation constant for bact. Fe/C   xkferb    =', xkferb
      ENDIF
      !
      denitr(:,:,:) = 0.
      !
   END SUBROUTINE p4z_rem_init


   INTEGER FUNCTION p4z_rem_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( denitr(PRIV_3D_BIOARRAY), STAT=p4z_rem_alloc )
      !
      IF( p4z_rem_alloc /= 0 )   CALL ctl_warn( 'p4z_rem_alloc: failed to allocate arrays' )
      !
   END FUNCTION p4z_rem_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_rem                    ! Empty routine
   END SUBROUTINE p4z_rem
#endif 

   !!======================================================================
END MODULE p4zrem
