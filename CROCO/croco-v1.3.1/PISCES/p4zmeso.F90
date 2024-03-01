#include "cppdefs.h"

MODULE p4zmeso
   !!======================================================================
   !!                         ***  MODULE p4zmeso  ***
   !! TOP :   PISCES Compute the sources/sinks for mesozooplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!   p4z_meso       : Compute the sources/sinks for mesozooplankton
   !!   p4z_meso_init  : Initialization of the parameters for mesozooplankton
   !!----------------------------------------------------------------------
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p4zprod         ! production
!  USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_meso              ! called in p4zbio.F90
   PUBLIC   p4z_meso_init         ! called in trcsms_pisces.F90

!!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"

   REAL(wp), PUBLIC ::  part2        !: part of calcite not dissolved in mesozoo guts
   REAL(wp), PUBLIC ::  xpref2d      !: mesozoo preference for diatoms
   REAL(wp), PUBLIC ::  xpref2n      !: mesozoo preference for nanophyto
   REAL(wp), PUBLIC ::  xpref2z      !: mesozoo preference for microzooplankton
   REAL(wp), PUBLIC ::  xpref2c      !: mesozoo preference for POC 
   REAL(wp), PUBLIC ::  xthresh2zoo  !: zoo feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2dia  !: diatoms feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2phy  !: nanophyto feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2poc  !: poc feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2     !: feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  resrat2      !: exsudation rate of mesozooplankton
   REAL(wp), PUBLIC ::  mzrat2       !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::  grazrat2     !: maximal mesozoo grazing rate
   REAL(wp), PUBLIC ::  xkgraz2      !: non assimilated fraction of P by mesozoo 
   REAL(wp), PUBLIC ::  unass2       !: Efficicency of mesozoo growth 
   REAL(wp), PUBLIC ::  sigma2       !: Fraction of mesozoo excretion as DOM 
   REAL(wp), PUBLIC ::  epsher2      !: growth efficiency
   REAL(wp), PUBLIC ::  epsher2min   !: minimum growth efficiency at high food for grazing 2
   REAL(wp), PUBLIC ::  grazflux     !: mesozoo flux feeding rate

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zmeso.F90 10367 2018-12-03 11:35:19Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_meso( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_meso  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for mesozooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step and ???
      !
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompadi, zcompaph, zcompapoc, zcompaz, zcompam
      REAL(wp) :: zgraze2 , zdenom, zdenom2
      REAL(wp) :: zfact   , zfood, zfoodlim, zproport, zbeta
      REAL(wp) :: zmortzgoc, zfrac, zfracfe, zratio, zratio2, zfracal, zgrazcal
      REAL(wp) :: zepsherf, zepshert, zepsherv, zgrarsig, zgraztotc, zgraztotn, zgraztotf
      REAL(wp) :: zgrarem2, zgrafer2, zgrapoc2, zprcaca, zmortz, zgrasrat, zgrasratn
      REAL(wp) :: zrespz, ztortz, zgrazd, zgrazz, zgrazpof
      REAL(wp) :: zgrazn, zgrazpoc, zgraznf, zgrazf
      REAL(wp) :: zgrazfffp, zgrazfffg, zgrazffep, zgrazffeg
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zgrazing, zfezoo2
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zw3d, zz2ligprod
      !!---------------------------------------------------------------------
      !
      zgrazing(:,:,:) = 0.
      zfezoo2 (:,:,:) = 0.
      !
      IF (ln_ligand) THEN
         ALLOCATE( zz2ligprod(PRIV_3D_BIOARRAY) )
         zz2ligprod(:,:,:) = 0.
      ENDIF
      !
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zcompam   = MAX( ( trb(ji,jj,K,jpmes) - 1.e-9 ), 0.e0 )
               zfact     = xstep * tgfunc2(ji,jj,jk) * zcompam

               !  Respiration rates of both zooplankton
               !  -------------------------------------
               zrespz    = resrat2 * zfact * ( trb(ji,jj,K,jpmes) / ( xkmort  &
               &           + trb(ji,jj,K,jpmes) )  &
               &           + 3. * nitrfac(ji,jj,jk) )

               !  Zooplankton mortality. A square function has been selected with
               !  no real reason except that it seems to be more stable and may mimic predation
               !  ---------------------------------------------------------------
               ztortz    = mzrat2 * 1.e6 * zfact * trb(ji,jj,K,jpmes)  * (1. - nitrfac(ji,jj,jk) )
               !
               zcompadi  = MAX( ( trb(ji,jj,K,jpdia) - xthresh2dia ), 0.e0 )
               zcompaz   = MAX( ( trb(ji,jj,K,jpzoo) - xthresh2zoo ), 0.e0 )
               zcompapoc = MAX( ( trb(ji,jj,K,jppoc) - xthresh2poc ), 0.e0 )
               ! Size effect of nanophytoplankton on grazing : the smaller it is, the less prone
               ! it is to predation by mesozooplankton
               ! -------------------------------------------------------------------------------
               zcompaph  = MAX( ( trb(ji,jj,K,jpphy) - xthresh2phy ), 0.e0 ) &
                  &      * MIN(1., MAX( 0., ( quotan(ji,jj,jk) - 0.2) / 0.3 ) )

               !   Mesozooplankton grazing
               !   ------------------------
               zfood     = xpref2d * zcompadi + xpref2z * zcompaz + xpref2n * zcompaph + xpref2c * zcompapoc 
               zfoodlim  = MAX( 0., zfood - MIN( 0.5 * zfood, xthresh2 ) )
               zdenom    = zfoodlim / ( xkgraz2 + zfoodlim )
               zdenom2   = zdenom / ( zfood + rtrn )
               zgraze2   = grazrat2 * xstep * tgfunc2(ji,jj,jk) * trb(ji,jj,K,jpmes) &
               &           * (1. - nitrfac(ji,jj,jk)) 

               zgrazd    = zgraze2  * xpref2d  * zcompadi  * zdenom2 
               zgrazz    = zgraze2  * xpref2z  * zcompaz   * zdenom2 
               zgrazn    = zgraze2  * xpref2n  * zcompaph  * zdenom2 
               zgrazpoc  = zgraze2  * xpref2c  * zcompapoc * zdenom2 

               zgraznf   = zgrazn   * trb(ji,jj,K,jpnfe)   &
               &          / ( trb(ji,jj,K,jpphy) + rtrn)
               zgrazf    = zgrazd   * trb(ji,jj,K,jpdfe)   &
               &          / ( trb(ji,jj,K,jpdia) + rtrn)
               zgrazpof  = zgrazpoc * trb(ji,jj,K,jpsfe)   &
               &          / ( trb(ji,jj,K,jppoc) + rtrn)

               !  Mesozooplankton flux feeding on GOC
               !  ----------------------------------
               zgrazffeg = grazflux  * xstep * wsbio4(ji,jj,jk)      &
               &           * tgfunc2(ji,jj,jk) * trb(ji,jj,K,jpgoc)  &
               &           * trb(ji,jj,K,jpmes)                      &
               &           * (1. - nitrfac(ji,jj,jk))
               zgrazfffg = zgrazffeg * trb(ji,jj,K,jpbfe)   &
               &           / (trb(ji,jj,K,jpgoc) + rtrn)
               zgrazffep = grazflux  * xstep *  wsbio3(ji,jj,jk)      &
               &           * tgfunc2(ji,jj,jk) * trb(ji,jj,K,jppoc)   &
               &           * trb(ji,jj,K,jpmes)                       &
               &           * (1. - nitrfac(ji,jj,jk))
               zgrazfffp = zgrazffep * trb(ji,jj,K,jpsfe)   &
               &           / (trb(ji,jj,K,jppoc) + rtrn)
               !
               zgraztotc = zgrazd + zgrazz + zgrazn + zgrazpoc + zgrazffep + zgrazffeg
               ! Compute the proportion of filter feeders
               zproport  = (zgrazffep + zgrazffeg)/(rtrn + zgraztotc)
               ! Compute fractionation of aggregates. It is assumed that 
               ! diatoms based aggregates are more prone to fractionation
               ! since they are more porous (marine snow instead of fecal pellets)
               zratio    = trb(ji,jj,K,jpgsi) / ( trb(ji,jj,K,jpgoc) + rtrn )
               zratio2   = zratio * zratio
               zfrac     = zproport * grazflux  * xstep * wsbio4(ji,jj,jk)      &
               &          * trb(ji,jj,K,jpgoc) * trb(ji,jj,K,jpmes)          &
               &          * ( 0.2 + 3.8 * zratio2 / ( 1.**2 + zratio2 ) )
               zfracfe   = zfrac * trb(ji,jj,K,jpbfe) / (trb(ji,jj,K,jpgoc) + rtrn)

               zgrazffep = zproport * zgrazffep
               zgrazffeg = zproport * zgrazffeg
               zgrazfffp = zproport * zgrazfffp
               zgrazfffg = zproport * zgrazfffg
               zgraztotc = zgrazd + zgrazz + zgrazn + zgrazpoc + zgrazffep + zgrazffeg
               zgraztotn = zgrazd * quotad(ji,jj,jk) + zgrazz + zgrazn * quotan(ji,jj,jk)   &
               &   + zgrazpoc + zgrazffep + zgrazffeg
               zgraztotf = zgrazf + zgraznf + zgrazz * ferat3 + zgrazpof + zgrazfffp + zgrazfffg

               ! Total grazing ( grazing by microzoo is already computed in p4zmicro )
               zgrazing(ji,jj,jk) = zgraztotc

               !    Mesozooplankton efficiency
               !    --------------------------
               zgrasrat  =  ( zgraztotf + rtrn )/ ( zgraztotc + rtrn )
               zgrasratn =  ( zgraztotn + rtrn )/ ( zgraztotc + rtrn )
               zepshert  = MIN( 1., zgrasratn, zgrasrat / ferat3)
               zbeta     = MAX(0., (epsher2 - epsher2min) )
               zepsherf  = epsher2min + zbeta / ( 1.0 + 0.04E6 * 12. * zfood * zbeta ) 
               zepsherv  = zepsherf * zepshert 

               zgrarem2  = zgraztotc * ( 1. - zepsherv - unass2 ) &
               &         + ( 1. - epsher2 - unass2 ) / ( 1. - epsher2 ) * ztortz
               zgrafer2  = zgraztotc * MAX( 0. , ( 1. - unass2 ) * zgrasrat - ferat3 * zepsherv )    &
               &         + ferat3 * ( ( 1. - epsher2 - unass2 ) /( 1. - epsher2 ) * ztortz )
               zgrapoc2  = zgraztotc * unass2

               !   Update the arrays TRA which contain the biological sources and sinks
               zgrarsig  = zgrarem2 * sigma2
               tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zgrarsig
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zgrarsig
               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zgrarem2 - zgrarsig
               !
               IF( ln_ligand ) THEN 
                  tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) + (zgrarem2 - zgrarsig) * ldocz
                  zz2ligprod(ji,jj,jk) = (zgrarem2 - zgrarsig) * ldocz
               ENDIF
               !
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2ut * zgrarsig
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zgrafer2
               zfezoo2(ji,jj,jk)   = zgrafer2
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrarsig
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zgrarsig              

               zmortz = ztortz + zrespz
               zmortzgoc = unass2 / ( 1. - epsher2 ) * ztortz + zrespz
               tra(ji,jj,jk,jpmes) = tra(ji,jj,jk,jpmes) - zmortz + zepsherv * zgraztotc 
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zgrazd
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - zgrazz
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zgrazn
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zgrazn * trb(ji,jj,K,jpnch)   &
               &                   / ( trb(ji,jj,K,jpphy) + rtrn )
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - zgrazd * trb(ji,jj,K,jpdch)   &
               &                   / ( trb(ji,jj,K,jpdia) + rtrn )
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) - zgrazd * trb(ji,jj,K,jpdsi)   &
               &                   / ( trb(ji,jj,K,jpdia) + rtrn )
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) + zgrazd * trb(ji,jj,K,jpdsi)   &
               &                   / ( trb(ji,jj,K,jpdia) + rtrn )
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zgraznf
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zgrazf

               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zgrazpoc - zgrazffep + zfrac
               prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zfrac
               conspoc(ji,jj,jk) = conspoc(ji,jj,jk) - zgrazpoc - zgrazffep
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zmortzgoc - zgrazffeg + zgrapoc2 - zfrac
               prodgoc(ji,jj,jk) = prodgoc(ji,jj,jk) + zmortzgoc + zgrapoc2
               consgoc(ji,jj,jk) = consgoc(ji,jj,jk) - zgrazffeg - zfrac
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) - zgrazpof - zgrazfffp + zfracfe
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + ferat3 * zmortzgoc - zgrazfffg     &
                 &                + zgraztotf * unass2 - zfracfe
               zfracal = trb(ji,jj,K,jpcal) / (trb(ji,jj,K,jppoc)   &
               &       + trb(ji,jj,K,jpgoc) + rtrn )
               zgrazcal = (zgrazffeg + zgrazpoc) * (1. - part2) * zfracal
               ! calcite production
               zprcaca = xfracal(ji,jj,jk) * zgrazn
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
               !
               zprcaca = part2 * zprcaca
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrazcal - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2. * ( zgrazcal + zprcaca )
               tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) - zgrazcal + zprcaca
            END DO
         END DO
      END DO
      !
#if defined key_iomput
      IF( lk_iomput .AND. knt == nrdttrc ) THEN
         ALLOCATE( zw3d(PRIV_3D_BIOARRAY) )
         IF( iom_use( "GRAZ2" ) ) THEN
            zw3d(:,:,:) = zgrazing(:,:,:) * 1.e+3 * rfact2r * tmask(:,:,:)  !   Total grazing of phyto by zooplankton
            CALL iom_put( "GRAZ2", zw3d )
         ENDIF
         IF( iom_use( "PCAL" ) ) THEN
            zw3d(:,:,:) = prodcal(:,:,:) * 1.e+3 * rfact2r * tmask(:,:,:)   !  Calcite production
            CALL iom_put( "PCAL", zw3d )  
         ENDIF
         IF( iom_use( "FEZOO2" ) ) THEN
            zw3d(:,:,:) = zfezoo2(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:)   !
            CALL iom_put( "FEZOO2", zw3d )
         ENDIF
         IF( iom_use( "LPRODZ2" ) .AND. ln_ligand )  THEN
            zw3d(:,:,:) = zz2ligprod(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:)
            CALL iom_put( "LPRODZ2"  , zw3d )
         ENDIF
         DEALLOCATE( zw3d )
      ENDIF
#endif
      !
#if defined key_trc_diaadd
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               trc3d(ji,jj,K,jp_grapoc2) = zgrazing(ji,jj,jk) * 1.e+3 * rfact2r * tmask(ji,jj,jk) !  grazing of phyto by mesozoo 
            END DO
         END DO
      END DO

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               trc3d(ji,jj,K,jp_meso2) = zgrazing(ji,jj,jk) * ( 1. - epsher2 - unass2 ) &
                     &                   * (-o2ut) * sigma2 * 1.e+3 * rfact2r * tmask(ji,jj,jk) ! o2 consumption by Mesozoo
            END DO
         END DO
      END DO
#endif

      IF (ln_ligand)  DEALLOCATE( zz2ligprod )
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('meso')")
        CALL prt_ctl_trc_info(charout)
        CALL prt_ctl_trc( charout, ltra='tra')
!        CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
   END SUBROUTINE p4z_meso


   SUBROUTINE p4z_meso_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_meso_init  ***
      !!
      !! ** Purpose :   Initialization of mesozooplankton parameters
      !!
      !! ** Method  :   Read the nampismes namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampismes
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp4zmes/ part2, grazrat2, resrat2, mzrat2, xpref2n, xpref2d, xpref2z,   &
         &                xpref2c, xthresh2dia, xthresh2phy, xthresh2zoo, xthresh2poc, &
         &                xthresh2, xkgraz2, epsher2, epsher2min, sigma2, unass2, grazflux
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*) 
         WRITE(numout,*) 'p4z_meso_init : Initialization of mesozooplankton parameters'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampismes in reference namelist : Pisces mesozooplankton
      READ  ( numnatp_ref, namp4zmes, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp4zmes in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampismes in configuration namelist : Pisces mesozooplankton
      READ  ( numnatp_cfg, namp4zmes, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp4zmes in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, namp4zmes )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp4zmes'
         WRITE(numout,*) '      part of calcite not dissolved in mesozoo guts  part2        =', part2
         WRITE(numout,*) '      mesozoo preference for phyto                   xpref2n      =', xpref2n
         WRITE(numout,*) '      mesozoo preference for diatoms                 xpref2d      =', xpref2d
         WRITE(numout,*) '      mesozoo preference for zoo                     xpref2z      =', xpref2z
         WRITE(numout,*) '      mesozoo preference for poc                     xpref2c      =', xpref2c
         WRITE(numout,*) '      microzoo feeding threshold  for mesozoo        xthresh2zoo  =', xthresh2zoo
         WRITE(numout,*) '      diatoms feeding threshold  for mesozoo         xthresh2dia  =', xthresh2dia
         WRITE(numout,*) '      nanophyto feeding threshold for mesozoo        xthresh2phy  =', xthresh2phy
         WRITE(numout,*) '      poc feeding threshold for mesozoo              xthresh2poc  =', xthresh2poc
         WRITE(numout,*) '      feeding threshold for mesozooplankton          xthresh2     =', xthresh2
         WRITE(numout,*) '      exsudation rate of mesozooplankton             resrat2      =', resrat2
         WRITE(numout,*) '      mesozooplankton mortality rate                 mzrat2       =', mzrat2
         WRITE(numout,*) '      maximal mesozoo grazing rate                   grazrat2     =', grazrat2
         WRITE(numout,*) '      mesozoo flux feeding rate                      grazflux     =', grazflux
         WRITE(numout,*) '      non assimilated fraction of P by mesozoo       unass2       =', unass2
         WRITE(numout,*) '      Efficiency of Mesozoo growth                   epsher2      =', epsher2
         WRITE(numout,*) '      Minimum Efficiency of Mesozoo growth           epsher2min  =', epsher2min
         WRITE(numout,*) '      Fraction of mesozoo excretion as DOM           sigma2       =', sigma2
         WRITE(numout,*) '      half sturation constant for grazing 2          xkgraz2      =', xkgraz2
      ENDIF
      !
   END SUBROUTINE p4z_meso_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_meso                    ! Empty routine
   END SUBROUTINE p4z_meso
#endif 

   !!======================================================================
END MODULE p4zmeso
