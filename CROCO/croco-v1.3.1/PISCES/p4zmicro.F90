#include "cppdefs.h"

MODULE p4zmicro
   !!======================================================================
   !!                         ***  MODULE p4zmicro  ***
   !! TOP :   PISCES Compute the sources/sinks for microzooplankton
   !!======================================================================
#if defined key_pisces
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
   !!   p4z_micro      : Compute the sources/sinks for microzooplankton
   !!   p4z_micro_init : Initialize and read the appropriate namelist
   !!----------------------------------------------------------------------
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p4zlim          ! Co-limitations
   USE p4zprod         ! production
!   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_micro         ! called in p4zbio.F90
   PUBLIC   p4z_micro_init    ! called in trcsms_pisces.F90

   !!* Substitution
#  include "top_substitute.h90"
#  include "ocean2pisces.h90"

   REAL(wp), PUBLIC ::   part        !: part of calcite not dissolved in microzoo guts
   REAL(wp), PUBLIC ::   xprefc      !: microzoo preference for POC 
   REAL(wp), PUBLIC ::   xprefn      !: microzoo preference for nanophyto
   REAL(wp), PUBLIC ::   xprefd      !: microzoo preference for diatoms
   REAL(wp), PUBLIC ::   xthreshdia  !: diatoms feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::   xthreshphy  !: nanophyto threshold for microzooplankton 
   REAL(wp), PUBLIC ::   xthreshpoc  !: poc threshold for microzooplankton 
   REAL(wp), PUBLIC ::   xthresh     !: feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::   resrat      !: exsudation rate of microzooplankton
   REAL(wp), PUBLIC ::   mzrat       !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::   grazrat     !: maximal microzoo grazing rate
   REAL(wp), PUBLIC ::   xkgraz      !: Half-saturation constant of assimilation
   REAL(wp), PUBLIC ::   unass       !: Non-assimilated part of food
   REAL(wp), PUBLIC ::   sigma1      !: Fraction of microzoo excretion as DOM 
   REAL(wp), PUBLIC ::   epsher      !: growth efficiency for grazing 1 
   REAL(wp), PUBLIC ::   epshermin   !: minimum growth efficiency for grazing 1

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zmicro.F90 10374 2018-12-06 09:49:35Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_micro( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_micro  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for microzooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time step
      INTEGER, INTENT(in) ::   knt   ! ??? 
      !
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompadi, zcompaz , zcompaph, zcompapoc
      REAL(wp) :: zgraze  , zdenom, zdenom2
      REAL(wp) :: zfact   , zfood, zfoodlim, zbeta
      REAL(wp) :: zepsherf, zepshert, zepsherv, zgrarsig, zgraztotc, zgraztotn, zgraztotf
      REAL(wp) :: zgrarem, zgrafer, zgrapoc, zprcaca, zmortz
      REAL(wp) :: zrespz, ztortz, zgrasrat, zgrasratn
      REAL(wp) :: zgrazp, zgrazm, zgrazsd
      REAL(wp) :: zgrazmf, zgrazsf, zgrazpf
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zgrazing, zfezoo
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zw3d, zzligprod
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF (ln_ligand) THEN
         ALLOCATE( zzligprod(PRIV_3D_BIOARRAY) )
         zzligprod(:,:,:) = 0.
      ENDIF
      !
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zcompaz = MAX( ( trb(ji,jj,K,jpzoo) - 1.e-9 ), 0.e0 )
               zfact   = xstep * tgfunc2(ji,jj,jk) * zcompaz

               !  Respiration rates of both zooplankton
               !  -------------------------------------
               zrespz = resrat * zfact * trb(ji,jj,K,jpzoo)   &
                  &   / ( xkmort + trb(ji,jj,K,jpzoo) )       &
                  &   + resrat * zfact * 3. * nitrfac(ji,jj,jk)

               !  Zooplankton mortality. A square function has been selected with
               !  no real reason except that it seems to be more stable and may mimic predation.
               !  ---------------------------------------------------------------
               ztortz = mzrat * 1.e6 * zfact * trb(ji,jj,K,jpzoo) * (1. - nitrfac(ji,jj,jk))

               zcompadi  = MIN( MAX( ( trb(ji,jj,K,jpdia) - xthreshdia ), 0.e0 ), xsizedia )
               zcompaph  = MAX( ( trb(ji,jj,K,jpphy) - xthreshphy ), 0.e0 )
               zcompapoc = MAX( ( trb(ji,jj,K,jppoc) - xthreshpoc ), 0.e0 )
               
               !     Microzooplankton grazing
               !     ------------------------
               zfood     = xprefn * zcompaph + xprefc * zcompapoc + xprefd * zcompadi
               zfoodlim  = MAX( 0. , zfood - min(xthresh,0.5*zfood) )
               zdenom    = zfoodlim / ( xkgraz + zfoodlim )
               zdenom2   = zdenom / ( zfood + rtrn )
               zgraze    = grazrat * xstep * tgfunc2(ji,jj,jk)    &
               &           * trb(ji,jj,K,jpzoo) * (1. - nitrfac(ji,jj,jk))

               zgrazp    = zgraze  * xprefn * zcompaph  * zdenom2 
               zgrazm    = zgraze  * xprefc * zcompapoc * zdenom2 
               zgrazsd   = zgraze  * xprefd * zcompadi  * zdenom2 

               zgrazpf   = zgrazp  * trb(ji,jj,K,jpnfe)   &
               &          / (trb(ji,jj,K,jpphy) + rtrn)
               zgrazmf   = zgrazm  * trb(ji,jj,K,jpsfe)   &
               &          / (trb(ji,jj,K,jppoc) + rtrn)
               zgrazsf   = zgrazsd * trb(ji,jj,K,jpdfe)   &
               &          / (trb(ji,jj,K,jpdia) + rtrn)
               !
               zgraztotc = zgrazp  + zgrazm  + zgrazsd 
               zgraztotf = zgrazpf + zgrazsf + zgrazmf 
               zgraztotn = zgrazp * quotan(ji,jj,jk) + zgrazm + zgrazsd * quotad(ji,jj,jk)

               ! Grazing by microzooplankton
               zgrazing(ji,jj,jk) = zgraztotc

               !    Various remineralization and excretion terms
               !    --------------------------------------------
               zgrasrat  = ( zgraztotf + rtrn ) / ( zgraztotc + rtrn )
               zgrasratn = ( zgraztotn + rtrn ) / ( zgraztotc + rtrn )
               zepshert  =  MIN( 1., zgrasratn, zgrasrat / ferat3)
               zbeta     = MAX(0., (epsher - epshermin) )
               zepsherf  = epshermin + zbeta / ( 1.0 + 0.04E6 * 12. * zfood * zbeta )
               zepsherv  = zepsherf * zepshert 

               zgrafer   = zgraztotc * MAX( 0. , ( 1. - unass ) * zgrasrat - ferat3 * zepsherv ) 
               zgrarem   = zgraztotc * ( 1. - zepsherv - unass )
               zgrapoc   = zgraztotc * unass

               !  Update of the TRA arrays
               !  ------------------------
               zgrarsig  = zgrarem * sigma1
               tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zgrarsig
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zgrarsig
               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zgrarem - zgrarsig
               !
               IF( ln_ligand ) THEN
                  tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) + (zgrarem - zgrarsig) * ldocz
                  zzligprod(ji,jj,jk) = (zgrarem - zgrarsig) * ldocz
               ENDIF
               !
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2ut * zgrarsig
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zgrafer
               zfezoo(ji,jj,jk)    = zgrafer
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zgrapoc
               prodpoc(ji,jj,jk)   = prodpoc(ji,jj,jk) + zgrapoc
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zgraztotf * unass
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrarsig
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zgrarsig
               !   Update the arrays TRA which contain the biological sources and sinks
               !   --------------------------------------------------------------------
               zmortz = ztortz + zrespz
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - zmortz + zepsherv * zgraztotc 
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zgrazp
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zgrazsd
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch)   &
               &       - zgrazp  * trb(ji,jj,K,jpnch)/(trb(ji,jj,K,jpphy)+rtrn)
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch)   &
               &       - zgrazsd * trb(ji,jj,K,jpdch)/(trb(ji,jj,K,jpdia)+rtrn)
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi)   &
               &       - zgrazsd * trb(ji,jj,K,jpdsi)/(trb(ji,jj,K,jpdia)+rtrn)
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) + zgrazsd   &
               &       * trb(ji,jj,K,jpdsi)/(trb(ji,jj,K,jpdia)+rtrn)
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zgrazpf
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zgrazsf
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortz - zgrazm
               prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zmortz
               conspoc(ji,jj,jk) = conspoc(ji,jj,jk) - zgrazm
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + ferat3 * zmortz - zgrazmf
               !
               ! calcite production
               zprcaca = xfracal(ji,jj,jk) * zgrazp
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
               !
               zprcaca = part * zprcaca
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2. * zprcaca
               tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) + zprcaca
            END DO
         END DO
      END DO
      !
#if defined key_iomput
      IF( lk_iomput ) THEN
         IF( knt == nrdttrc ) THEN
           ALLOCATE( zw3d(PRIV_3D_BIOARRAY) )
           IF( iom_use( "GRAZ1" ) ) THEN
              zw3d(:,:,:) = zgrazing(:,:,:) * 1.e+3 * rfact2r * tmask(:,:,:)  !  Total grazing of phyto by zooplankton
              CALL iom_put( "GRAZ1", zw3d )
           ENDIF
           IF( iom_use( "FEZOO" ) ) THEN
              zw3d(:,:,:) = zfezoo(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:)   !
              CALL iom_put( "FEZOO", zw3d )
           ENDIF
           IF( iom_use( "LPRODZ" ) .AND. ln_ligand )  THEN
              zw3d(:,:,:) = zzligprod(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:)
              CALL iom_put( "LPRODZ"  , zw3d )
           ENDIF
           DEALLOCATE( zw3d )
         ENDIF
      ENDIF
#endif
      !
#if defined key_trc_diaadd
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               trc3d(ji,jj,K,jp_grapoc) = zgrazing(ji,jj,jk) * 1.e+3 * rfact2r * tmask(ji,jj,jk) !  grazing of phyto by microzoo
            END DO
         END DO
      END DO

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               trc3d(ji,jj,K,jp_mico2) = zgrazing(ji,jj,jk) * ( 1.- epsher - unass ) &
                  &                      * (-o2ut) * sigma1 * 1.e+3 * rfact2r * tmask(ji,jj,jk)   ! o2 consumption by Microzoo
            END DO
         END DO
      END DO
#endif

      IF (ln_ligand)  DEALLOCATE( zzligprod )
      !
      IF(ln_ctl) THEN      ! print mean trends (used for debugging)
         WRITE(charout, FMT="('micro')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( charout, ltra='tra')
!        CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
   END SUBROUTINE p4z_micro


   SUBROUTINE p4z_micro_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_micro_init  ***
      !!
      !! ** Purpose :   Initialization of microzooplankton parameters
      !!
      !! ** Method  :   Read the nampiszoo namelist and check the parameters
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampiszoo
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp4zzoo/ part, grazrat, resrat, mzrat, xprefn, xprefc, &
         &                xprefd,  xthreshdia,  xthreshphy,  xthreshpoc, &
         &                xthresh, xkgraz, epsher, epshermin, sigma1, unass
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*) 
         WRITE(numout,*) 'p4z_micro_init : Initialization of microzooplankton parameters'
         WRITE(numout,*) '~~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampiszoo in reference namelist : Pisces microzooplankton
      READ  ( numnatp_ref, namp4zzoo, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp4zzoo in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampiszoo in configuration namelist : Pisces microzooplankton
      READ  ( numnatp_cfg, namp4zzoo, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp4zzoo in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, namp4zzoo )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp4zzoo'
         WRITE(numout,*) '      part of calcite not dissolved in microzoo guts  part        =', part
         WRITE(numout,*) '      microzoo preference for POC                     xprefc      =', xprefc
         WRITE(numout,*) '      microzoo preference for nano                    xprefn      =', xprefn
         WRITE(numout,*) '      microzoo preference for diatoms                 xprefd      =', xprefd
         WRITE(numout,*) '      diatoms feeding threshold  for microzoo         xthreshdia  =', xthreshdia
         WRITE(numout,*) '      nanophyto feeding threshold for microzoo        xthreshphy  =', xthreshphy
         WRITE(numout,*) '      poc feeding threshold for microzoo              xthreshpoc  =', xthreshpoc
         WRITE(numout,*) '      feeding threshold for microzooplankton          xthresh     =', xthresh
         WRITE(numout,*) '      exsudation rate of microzooplankton             resrat      =', resrat
         WRITE(numout,*) '      microzooplankton mortality rate                 mzrat       =', mzrat
         WRITE(numout,*) '      maximal microzoo grazing rate                   grazrat     =', grazrat
         WRITE(numout,*) '      non assimilated fraction of P by microzoo       unass       =', unass
         WRITE(numout,*) '      Efficicency of microzoo growth                  epsher      =', epsher
         WRITE(numout,*) '      Minimum efficicency of microzoo growth          epshermin   =', epshermin
         WRITE(numout,*) '      Fraction of microzoo excretion as DOM           sigma1      =', sigma1
         WRITE(numout,*) '      half sturation constant for grazing 1           xkgraz      =', xkgraz
      ENDIF
      !
   END SUBROUTINE p4z_micro_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_micro                    ! Empty routine
   END SUBROUTINE p4z_micro
#endif 

   !!======================================================================
END MODULE p4zmicro
