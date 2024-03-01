#include "cppdefs.h"

MODULE p4zflx
   !!======================================================================
   !!                         ***  MODULE p4zflx  ***
   !! TOP :   PISCES CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
   !!======================================================================
   !! History :   -   !  1988-07  (E. MAIER-REIMER) Original code
   !!             -   !  1998     (O. Aumont) additions
   !!             -   !  1999     (C. Le Quere) modifications
   !!            1.0  !  2004     (O. Aumont) modifications
   !!            2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                 !  2011-02  (J. Simeon, J. Orr) Include total atm P correction 
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!   p4z_flx       :   CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
   !!   p4z_flx_init  :   Read the namelist
   !!----------------------------------------------------------------------
   USE sms_pisces     !  PISCES Source Minus Sink variables
   USE p4zche         !  Chemical model
!   USE iom            !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_flx  
   PUBLIC   p4z_flx_init  
   PUBLIC   p4z_flx_alloc  

   !!* Substitution
#include "ocean2pisces.h90"
#include "top_substitute.h90"

   !                                 !!** Namelist  nampisext  **
   REAL(wp)          ::   atcco2      !: pre-industrial atmospheric [co2] (ppm) 	
   !                                  !!* nampisatm namelist (Atmospheric PRessure) *
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  patm, satmco2   !: atmospheric pco2 

   REAL(wp) ::   xconv  = 0.01 / 3600.   !: coefficients for conversion 

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zflx.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_flx ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_flx  ***
      !!
      !! ** Purpose :   CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
      !!
      !! ** Method  : 
      !!              - Include total atm P correction via Esbensen & Kushnir (1981) 
      !!              - Remove Wanninkhof chemical enhancement;
      !!              - Add option for time-interpolation of atcco2.txt  
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   !
      !
      INTEGER  ::   ji, jj, jm, iind, iindm1
      REAL(wp) ::   ztc, ztc2, ztc3, ztc4, zws, zkgwan
      REAL(wp) ::   zfld, zflu, zfld16, zflu16, zfact
      REAL(wp) ::   zvapsw, zsal, zfco2, zxc2, xCO2approx, ztkel, zfugcoeff
      REAL(wp) ::   zph, zdic, zsch_o2, zsch_co2
      REAL(wp) ::   zyr_dec, zdco2dt
      CHARACTER (len=25) ::   charout
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY) ::   zkgco2, zkgo2, zh2co3, &
      &                                          zoflx,  zpco2atm, oce_co2, zdpco2
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zw2d
      !!---------------------------------------------------------------------
      !
      ! SURFACE CHEMISTRY (PCO2 AND [H+] IN
      !     SURFACE LAYER); THE RESULT OF THIS CALCULATION
      !     IS USED TO COMPUTE AIR-SEA FLUX OF CO2

      DO jj = JRANGE
         DO ji = IRANGE
            ! DUMMY VARIABLES FOR DIC, H+, AND BORATE
            zfact = rhop(ji,jj,KSURF) / 1000. + rtrn
            zdic  = trb(ji,jj,KSURF,jpdic)
            zph   = MAX( hi(ji,jj,1), 1.e-10 ) / zfact
            ! CALCULATE [H2CO3]
            zh2co3(ji,jj) = zdic/(1. + ak13(ji,jj,1)/zph + ak13(ji,jj,1)*ak23(ji,jj,1)/zph**2)
         END DO
      END DO

      ! --------------
      ! COMPUTE FLUXES
      ! --------------

      ! FIRST COMPUTE GAS EXCHANGE COEFFICIENTS
      ! -------------------------------------------

      DO jj = JRANGE
         DO ji = IRANGE
            ztc  = MIN( 35., tsn(ji,jj,KSURF,jp_tem) )
            ztc2 = ztc * ztc
            ztc3 = ztc * ztc2 
            ztc4 = ztc2 * ztc2 
            ! Compute the schmidt Number both O2 and CO2
            zsch_co2 = 2116.8 - 136.25 * ztc + 4.7353 * ztc2 - 0.092307 * ztc3 + 0.0007555 * ztc4
            zsch_o2  = 1920.4 - 135.6  * ztc + 5.2122 * ztc2 - 0.109390 * ztc3 + 0.0009377 * ztc4
            !  wind speed 
            zws  = wndm(ji,jj)      &
            &    * wndm(ji,jj)
            ! Compute the piston velocity for O2 and CO2
            zkgwan = 0.251 * zws
            zkgwan = zkgwan * xconv * ( 1.- fr_i(ji,jj) ) * tmask(ji,jj,1)
            ! compute gas exchange for CO2 and O2
            zkgco2(ji,jj) = zkgwan * SQRT( 660./ zsch_co2 )
            zkgo2 (ji,jj) = zkgwan * SQRT( 660./ zsch_o2 )
         END DO
      END DO


      DO jj = JRANGE
         DO ji = IRANGE
            ztkel = tempis(ji,jj,1) + 273.15
            zsal  = salinprac(ji,jj,1) + ( 1.- tmask(ji,jj,1) ) * 35.
            zvapsw    = EXP(24.4543 - 67.4509*(100.0/ztkel) - 4.8489*LOG(ztkel/100) - 0.000544*zsal)
            zpco2atm(ji,jj) = satmco2(ji,jj) * ( patm(ji,jj) - zvapsw )
            zxc2      = ( 1.0 - zpco2atm(ji,jj) * 1E-6 )**2
            zfugcoeff = EXP( patm(ji,jj) * (chemc(ji,jj,2) + 2.0 * zxc2 * chemc(ji,jj,3) )   &
            &           / ( 82.05736 * ztkel ))
            zfco2 = zpco2atm(ji,jj) * zfugcoeff

            ! Compute CO2 flux for the sea and air
            zfld = zfco2 * chemc(ji,jj,1) * zkgco2(ji,jj)  ! (mol/L) * (m/s)
            zflu = zh2co3(ji,jj) * zkgco2(ji,jj)                                   ! (mol/L) (m/s) ?
            oce_co2(ji,jj) = ( zfld - zflu ) * rfact2 * tmask(ji,jj,1) * 1000.

            ! compute the trend
            tra(ji,jj,1,jpdic) = tra(ji,jj,1,jpdic) + ( zfld - zflu ) * rfact2 / e3t_n(ji,jj,KSURF) * tmask(ji,jj,1)

            ! Compute O2 flux 
            zfld16 = patm(ji,jj) * chemo2(ji,jj,1) * zkgo2(ji,jj)          ! (mol/L) * (m/s)
            zflu16 = trb(ji,jj,KSURF,jpoxy) * zkgo2(ji,jj)
            zoflx(ji,jj) = ( zfld16 - zflu16 ) * tmask(ji,jj,1)
            tra(ji,jj,1,jpoxy) = tra(ji,jj,1,jpoxy) + zoflx(ji,jj) * rfact2 / e3t_n(ji,jj,KSURF)
         END DO
      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('flx ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( charout, ltra='tra')
!         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
 
#if defined key_iomput
      IF( lk_iomput .AND. knt == nrdttrc ) THEN
         ALLOCATE( zw2d(PRIV_2D_BIOARRAY) )  
         IF( iom_use( "Cflx"  ) )  THEN
            zw2d(:,:) = oce_co2(:,:) * rfact2r
            CALL iom_put( "Cflx"     , zw2d )
         ENDIF

         IF( iom_use( "Oflx"  ) )  THEN
            zw2d(:,:) =  zoflx(:,:) * 1000 * tmask(:,:,1)
            CALL iom_put( "Oflx" , zw2d )
         ENDIF
         IF( iom_use( "Kg"    ) )  THEN
            zw2d(:,:) =  zkgco2(:,:) * tmask(:,:,1)
            CALL iom_put( "Kg"   , zw2d )
         ENDIF
         IF( iom_use( "Dpco2" ) ) THEN
           zw2d(:,:) = ( zpco2atm(:,:) - zh2co3(:,:) / ( chemc(:,:,1) + rtrn ) ) * tmask(:,:,1)
           CALL iom_put( "Dpco2" ,  zw2d )
         ENDIF
         IF( iom_use( "Dpo2" ) )  THEN
           zw2d(:,:) = ( atcox * patm(:,:) - atcox * trb(:,:,KSURF,jpoxy) / ( chemo2(:,:,1) + rtrn ) ) * tmask(:,:,1)
           CALL iom_put( "Dpo2"  , zw2d )
         ENDIF
         !
         DEALLOCATE( zw2d )
      ENDIF
#endif
      !
#if defined key_trc_diaadd
      DO jj = JRANGE
         DO ji = IRANGE
            ! Save diagnostics
            zdpco2(ji,jj) = (zpco2atm(ji,jj) - zh2co3(ji,jj) &
             &          / ( chemc(ji,jj,1) + rtrn ) ) * tmask(ji,jj,1)
            !
            trc2d(ji,jj,jp_flxco2) =  oce_co2(ji,jj)          !  carbon flux
            trc2d(ji,jj,jp_flxo2 ) =  zoflx(ji,jj)  * 1000.   !  O2 flux
            trc2d(ji,jj,jp_kgco2 ) =  zkgco2(ji,jj)           !  gas exchange for CO2
            trc2d(ji,jj,jp_dpco2 ) =  zdpco2(ji,jj)           ! delta pco2
         END DO
      END DO
#endif

   END SUBROUTINE p4z_flx


   SUBROUTINE p4z_flx_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_flx_init  ***
      !!
      !! ** Purpose :   Initialization of atmospheric conditions
      !!
      !! ** Method  :   Read the nampisext namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampisext
      !!----------------------------------------------------------------------
      INTEGER ::   jm, ios   ! Local integer 
      !!
      NAMELIST/nampisext/atcco2
      !!----------------------------------------------------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' p4z_flx_init : atmospheric conditions for air-sea flux calculation'
         WRITE(numout,*) ' ~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampisext in reference namelist : Pisces atm. conditions
      READ  ( numnatp_ref, nampisext, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampisext in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampisext in configuration namelist : Pisces atm. conditions
      READ  ( numnatp_cfg, nampisext, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nampisext in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisext )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : nampisext --- parameters for air-sea exchange'
         WRITE(numout,*) '         Constant Atmospheric pCO2 value               atcco2    =', atcco2
      ENDIF
      satmco2(:,:)  = atcco2      ! Initialisation of atmospheric pco2
      patm(:,:) = 1.0
      !
   END SUBROUTINE p4z_flx_init

   INTEGER FUNCTION p4z_flx_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_flx_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( satmco2(PRIV_2D_BIOARRAY), patm(PRIV_2D_BIOARRAY), STAT=p4z_flx_alloc )
      !
      IF( p4z_flx_alloc /= 0 )   CALL ctl_warn( 'p4z_flx_alloc : failed to allocate arrays' )
      !
   END FUNCTION p4z_flx_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_flx                   ! Empty routine
   END SUBROUTINE p4z_flx
#endif 

   !!======================================================================
END MODULE p4zflx
