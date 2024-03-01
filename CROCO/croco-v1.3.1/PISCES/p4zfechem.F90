#include "cppdefs.h"

MODULE p4zfechem
   !!======================================================================
   !!                         ***  MODULE p4zfechem  ***
   !! TOP :   PISCES Compute iron chemistry and scavenging
   !!======================================================================
   !! History :   3.5  !  2012-07 (O. Aumont, A. Tagliabue, C. Ethe) Original code
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   p4z_fechem       : Compute remineralization/scavenging of iron
   !!   p4z_fechem_init  : Initialisation of parameters for remineralisation
   !!   p4z_fechem_alloc : Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p4zche          ! chemical model
   USE p4zsbc          ! Boundary conditions from sediments
!  USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_fechem        ! called in p4zbio.F90
   PUBLIC   p4z_fechem_init   ! called in trcsms_pisces.F90

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"

   LOGICAL          ::   ln_ligvar    !: boolean for variable ligand concentration following Tagliabue and voelker
   REAL(wp), PUBLIC ::   xlam1        !: scavenging rate of Iron 
   REAL(wp), PUBLIC ::   xlamdust     !: scavenging rate of Iron by dust 
   REAL(wp), PUBLIC ::   ligand       !: ligand concentration in the ocean 
   REAL(wp), PUBLIC ::   kfep         !: rate constant for nanoparticle formation

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zfechem.F90 10416 2018-12-19 11:45:43Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_fechem( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_fechem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of iron
      !!
      !! ** Method  :   A simple chemistry model of iron from Aumont and Bopp (2006)
      !!                based on one ligand and one inorganic form
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step
      !
      INTEGER  ::   ji, jj, jk, jic, jn
      REAL(wp) ::   zdep, zlam1a, zlam1b, zlamfac
      REAL(wp) ::   zkeq, zfeequi, zfesatur, zfecoll, fe3sol
      REAL(wp) ::   zdenom1, zscave, zaggdfea, zaggdfeb, zcoag
      REAL(wp) ::   ztrc, zdust
      REAL(wp) ::   zdenom2
      REAL(wp) ::   zzFeL1, zzFeL2, zzFe2, zzFeP, zzFe3, zzstrn2
      REAL(wp) ::   zrum, zcodel, zargu, zlight
      REAL(wp) ::   zkox, zkph1, zkph2, zph, zionic, ztligand
      REAL(wp) ::   za, zb, zc, zkappa1, zkappa2, za0, za1, za2
      REAL(wp) ::   zxs, zfunc, zp, zq, zd, zr, zphi, zfff, zp3, zq2
      REAL(wp) ::   ztfe, zoxy, zhplus, zxlam
      REAL(wp) ::   zaggliga, zaggligb
      REAL(wp) ::   dissol, zligco
      REAL(wp) :: zrfact2
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) ::   zTL1, zFe3, ztotlig, precip, zFeL1
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) ::   zcoll3d, zscav3d, zlcoll3d
      !!---------------------------------------------------------------------
      !
      zFe3 (:,:,:) = 0.
      zFeL1(:,:,:) = 0.
      zTL1 (:,:,:) = 0.

      ! Total ligand concentration : Ligands can be chosen to be constant or variable
      ! Parameterization from Tagliabue and Voelker (2011)
      ! -------------------------------------------------
      IF( ln_ligvar ) THEN
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  ztotlig(ji,jj,jk) =  0.09 * trb(ji,jj,K,jpdoc) * 1E6 + ligand * 1E9
                  ztotlig(ji,jj,jk) =  MIN( ztotlig(ji,jj,jk), 10. )
               END DO
            END DO
         END DO
      ELSE
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  IF( ln_ligand ) THEN  ;   ztotlig(ji,jj,jk) = trb(ji,jj,K,jplgw) * 1E9
                  ELSE                  ;   ztotlig(ji,jj,jk) = ligand * 1E9
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF

      ! ------------------------------------------------------------
      !  from Aumont and Bopp (2006)
      ! This model is based on one ligand and Fe' 
      ! Chemistry is supposed to be fast enough to be at equilibrium
      ! ------------------------------------------------------------
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zTL1(ji,jj,jk)  = ztotlig(ji,jj,jk)
               zkeq            = fekeq(ji,jj,jk)
               zfesatur        = zTL1(ji,jj,jk) * 1E-9
               ztfe            = trb(ji,jj,K,jpfer) 
               ! Fe' is the root of a 2nd order polynom
               zFe3 (ji,jj,jk) = ( -( 1. + zfesatur * zkeq - zkeq * ztfe )               &
                  &              + SQRT( ( 1. + zfesatur * zkeq - zkeq * ztfe )**2       &
                  &              + 4. * ztfe * zkeq) ) / ( 2. * zkeq )
               zFe3 (ji,jj,jk) = zFe3(ji,jj,jk) * 1E9
               zFeL1(ji,jj,jk) = MAX( 0., trb(ji,jj,K,jpfer) * 1E9 - zFe3(ji,jj,jk) )
           END DO
         END DO
      END DO
         !

      zdust = 0.         ! if no dust available
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               ! Scavenging rate of iron. This scavenging rate depends on the load of particles of sea water. 
               ! This parameterization assumes a simple second order kinetics (k[Particles][Fe]).
               ! Scavenging onto dust is also included as evidenced from the DUNE experiments.
               ! --------------------------------------------------------------------------------------
               zhplus  = max( rtrn, hi(ji,jj,jk) )
               fe3sol  = fesol(ji,jj,jk,1) * ( zhplus**3 + fesol(ji,jj,jk,2) * zhplus**2  &
               &         + fesol(ji,jj,jk,3) * zhplus + fesol(ji,jj,jk,4)     &
               &         + fesol(ji,jj,jk,5) / zhplus )
               !
               zfeequi = zFe3(ji,jj,jk) * 1E-9
               zhplus  = max( rtrn, hi(ji,jj,jk) )
               fe3sol  = fesol(ji,jj,jk,1) * ( zhplus**3 + fesol(ji,jj,jk,2) * zhplus**2  &
                  &         + fesol(ji,jj,jk,3) * zhplus + fesol(ji,jj,jk,4)     &
                  &         + fesol(ji,jj,jk,5) / zhplus )
               zfecoll = 0.5 * zFeL1(ji,jj,jk) * 1E-9
               ! precipitation of Fe3+, creation of nanoparticles
               precip(ji,jj,jk) = MAX( 0., ( zFe3(ji,jj,jk) * 1E-9 - fe3sol ) ) * kfep * xstep
               !
               ztrc   = ( trb(ji,jj,K,jppoc) + trb(ji,jj,K,jpgoc)   &
                  &   + trb(ji,jj,K,jpcal) + trb(ji,jj,K,jpgsi) ) * 1.e6 
               IF( ln_dust )  zdust  = dust(ji,jj) / ( wdust / rday ) * tmask(ji,jj,jk) &
               &  * EXP( -gdept_n(ji,jj,K) / 540. )
               IF (ln_ligand) THEN
                  zxlam  = xlam1 * MAX( 1.E-3, EXP(-2 * etot(ji,jj,jk) / 10. )   &
                     &   * (1. - EXP(-2 * trb(ji,jj,K,jpoxy) / 100.E-6 ) ))
               ELSE
                  zxlam  = xlam1 * 1.0
               ENDIF
               zlam1b = 3.e-5 + xlamdust * zdust + zxlam * ztrc
               zscave = zfeequi * zlam1b * xstep

               ! Compute the different ratios for scavenging of iron
               ! to later allocate scavenged iron to the different organic pools
               ! ---------------------------------------------------------
               zdenom1 = zxlam * trb(ji,jj,K,jppoc) / zlam1b
               zdenom2 = zxlam * trb(ji,jj,K,jpgoc) / zlam1b

               !  Increased scavenging for very high iron concentrations found near the coasts 
               !  due to increased lithogenic particles and let say it is unknown processes (precipitation, ...)
               !  -----------------------------------------------------------
               zlamfac = MAX( 0.e0, ( gphit(ji,jj) + 55.) / 30. )
               zlamfac = MIN( 1.  , zlamfac )
               zdep    = MIN( 1., 1000. / gdept_n(ji,jj,K) )
               zcoag   = 1E-4 * ( 1. - zlamfac ) * zdep * xstep * trb(ji,jj,K,jpfer)

               !  Compute the coagulation of colloidal iron. This parameterization 
               !  could be thought as an equivalent of colloidal pumping.
               !  It requires certainly some more work as it is very poorly constrained.
               !  ----------------------------------------------------------------
               zlam1a   = ( 0.369  * 0.3 * trb(ji,jj,K,jpdoc)                 &
                   &      + 102.4  * trb(ji,jj,K,jppoc) ) * xdiss(ji,jj,jk)   &
                   &      + ( 114.   * 0.3 * trb(ji,jj,K,jpdoc) )
               zaggdfea = zlam1a * xstep * zfecoll
               !
               zlam1b   = 3.53E3 * trb(ji,jj,K,jpgoc) * xdiss(ji,jj,jk)
               zaggdfeb = zlam1b * xstep * zfecoll
               !
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zscave - zaggdfea - zaggdfeb &
               &                     - zcoag - precip(ji,jj,jk)
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zscave * zdenom1 + zaggdfea
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zscave * zdenom2 + zaggdfeb
               zscav3d(ji,jj,jk)   = zscave
               zcoll3d(ji,jj,jk)   = zaggdfea + zaggdfeb
               !
            END DO
         END DO
      END DO
      !
      !  Define the bioavailable fraction of iron
      !  ----------------------------------------
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               biron(ji,jj,jk) = trb(ji,jj,K,jpfer) 
            END DO
         END DO
      END DO
      !
      IF( ln_ligand ) THEN
         !
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  zlam1a   = ( 0.369  * 0.3 * trb(ji,jj,K,jpdoc)                &
                      &    + 102.4  * trb(ji,jj,K,jppoc) ) * xdiss(ji,jj,jk)    &
                      &    + ( 114.   * 0.3 * trb(ji,jj,K,jpdoc) )
                  !
                  zlam1b   = 3.53E3 *   trb(ji,jj,K,jpgoc) * xdiss(ji,jj,jk)
                  zligco   = 0.5 * trb(ji,jj,K,jplgw)
                  zaggliga = zlam1a * xstep * zligco
                  zaggligb = zlam1b * xstep * zligco
                  tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) - zaggliga - zaggligb
                  zlcoll3d(ji,jj,jk)  = zaggliga + zaggligb
               END DO
            END DO
         END DO
         !
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  plig(ji,jj,jk) =  MAX( 0., ( ( zFeL1(ji,jj,jk) * 1E-9 ) / ( trb(ji,jj,K,jpfer) +rtrn ) ) )
               END DO
            END DO
         END DO
         !
      ENDIF
      !  Output of some diagnostics variables
      !     ---------------------------------
#if defined key_iomput
      IF( lk_iomput ) THEN
         IF( knt == nrdttrc ) THEN
            zrfact2 = 1.e3 * rfact2r  ! conversion from mol/L/timestep into mol/m3/s
            IF( iom_use("Fe3")    )  CALL iom_put("Fe3"    , zFe3   (:,:,:)       * tmask(:,:,:) )   ! Fe3+
            IF( iom_use("FeL1")   )  CALL iom_put("FeL1"   , zFeL1  (:,:,:)       * tmask(:,:,:) )   ! FeL1
            IF( iom_use("TL1")    )  CALL iom_put("TL1"    , zTL1   (:,:,:)       * tmask(:,:,:) )   ! TL1
            IF( iom_use("Totlig") )  CALL iom_put("Totlig" , ztotlig(:,:,:)       * tmask(:,:,:) )   ! TL
            IF( iom_use("Biron")  )  CALL iom_put("Biron"  , biron  (:,:,:)  * 1e9 * tmask(:,:,:) )   ! biron
            IF( iom_use("FESCAV") )  CALL iom_put("FESCAV" , zscav3d(:,:,:)  * 1e9 * tmask(:,:,:) * zrfact2 )
            IF( iom_use("FECOLL") )  CALL iom_put("FECOLL" , zcoll3d(:,:,:)  * 1e9 * tmask(:,:,:) * zrfact2 )
            IF( iom_use("LGWCOLL"))  CALL iom_put("LGWCOLL", zlcoll3d(:,:,:) * 1e9 * tmask(:,:,:) * zrfact2 )
         ENDIF
      ENDIF
#endif

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('fechem')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( charout, ltra='tra')
      ENDIF
      !
   END SUBROUTINE p4z_fechem


   SUBROUTINE p4z_fechem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_fechem_init  ***
      !!
      !! ** Purpose :   Initialization of iron chemistry parameters
      !!
      !! ** Method  :   Read the nampisfer namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisfer
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer 
      !!
      NAMELIST/nampisfer/ ln_ligvar, xlam1, xlamdust, ligand, kfep 
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_rem_init : Initialization of iron chemistry parameters'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )            ! Namelist nampisfer in reference namelist : Pisces iron chemistry
      READ  ( numnatp_ref, nampisfer, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampisfer in reference namelist', lwp )
      REWIND( numnatp_cfg )            ! Namelist nampisfer in configuration namelist : Pisces iron chemistry
      READ  ( numnatp_cfg, nampisfer, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nampisfer in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, nampisfer )

      IF(lwp) THEN                     ! control print
         WRITE(numout,*) '   Namelist : nampisfer'
         WRITE(numout,*) '      variable concentration of ligand          ln_ligvar    =', ln_ligvar
         WRITE(numout,*) '      scavenging rate of Iron                   xlam1        =', xlam1
         WRITE(numout,*) '      scavenging rate of Iron by dust           xlamdust     =', xlamdust
         WRITE(numout,*) '      ligand concentration in the ocean         ligand       =', ligand
         WRITE(numout,*) '      rate constant for nanoparticle formation  kfep         =', kfep
      ENDIF
      ! 
   END SUBROUTINE p4z_fechem_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_fechem                    ! Empty routine
   END SUBROUTINE p4z_fechem
#endif 

   
   !!======================================================================
END MODULE p4zfechem
