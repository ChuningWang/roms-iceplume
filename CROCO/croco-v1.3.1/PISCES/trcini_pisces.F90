#include "cppdefs.h"

MODULE trcini_pisces
   !!======================================================================
   !!                         ***  MODULE trcini_pisces  ***
   !! TOP :   initialisation of the PISCES biochemical model
   !!======================================================================
   !! History :    -   !  1988-07  (E. Maier-Reiner) Original code
   !!              -   !  1999-10  (O. Aumont, C. Le Quere)
   !!              -   !  2002     (O. Aumont)  PISCES
   !!             1.0  !  2005-03  (O. Aumont, A. El Moussaoui) F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec) from trcini.pisces.h90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !! trc_ini_pisces   : PISCES biochemical model initialisation
   !!----------------------------------------------------------------------
   USE sms_pisces      ! Source Minus Sink variables
   USE trcsms_pisces   !
   USE p4zsms          ! Main P4Z routine
   USE p4zche          !  Chemical model
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zopt          !  optical model
   USE p4zrem          !  Remineralisation of organic matter
   USE p4zflx          !  Gas exchange
   USE p4zlim          !  Co-limitations of differents nutrients
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p4zmicro        !  Sources and sinks of microzooplankton
   USE p4zmeso         !  Sources and sinks of mesozooplankton
   USE p4zmort         !  Mortality terms for phytoplankton
   USE p4zlys          !  Calcite saturation
   USE p4zsed          !  Sedimentation & burial
   USE p4zpoc          !  Remineralization of organic particles
   USE p4zligand       !  Remineralization of organic ligands
   USE p4zsbc          !  External supply of nutrients
   USE p4zfechem       !  Iron chemistry
   USE p5zlim          !  Co-limitations of differents nutrients
   USE p5zprod         !  Growth rate of the 2 phyto groups
   USE p5zmicro        !  Sources and sinks of microzooplankton
   USE p5zmeso         !  Sources and sinks of mesozooplankton
   USE p5zmort         !  Mortality terms for phytoplankton
   USE sedini          !  SEDIMENTS initialization routine
   USE sed


   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_pisces   ! called by trcini.F90 module
   PUBLIC   trc_nam_pisces   ! called by trcini.F90 module

   !! * Module variables
   REAL(wp), SAVE ::   sco2   =  2.312e-3
   REAL(wp), SAVE ::   alka0  =  2.426e-3
   REAL(wp), SAVE ::   oxyg0  =  177.6e-6
   REAL(wp), SAVE ::   po4    =  2.165e-6
   REAL(wp), SAVE ::   bioma0 =  1.000e-8
   REAL(wp), SAVE ::   silic1 =  91.51e-6
   REAL(wp), SAVE ::   no3    =  30.9e-6 * 7.625


CONTAINS

    SUBROUTINE trc_ini_pisces

      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_ini_pisces ***
      !!
      !! ** Purpose :   Initialisation of the PISCES biochemical model
      !!----------------------------------------------------------------------
      INTEGER  ::  ji, jj, jk, jn, ierr
      REAL(wp) ::  zcaralk, zbicarb, zco3
      REAL(wp) ::  ztmas, ztmas1


#if defined key_pisces_quota
      ln_p5z = .true.
      ln_p4z = .false.
#else
      ln_p4z = .true.
      ln_p5z = .false.
#endif
#if defined key_ligand
      ln_ligand = .true.
#else
      ln_ligand = .false.
#endif
#if defined key_sediment
      ln_sediment = .true.
#else
      ln_sediment = .false.
      ln_sed_2way = .false.
#endif

      IF(lwp) THEN
         WRITE(numout,*)
         IF( ln_p4z ) THEN
            WRITE(numout,*) 'p4z_ini :   PISCES biochemical model initialisation'
            WRITE(numout,*) '~~~~~~~'
         ELSE
            WRITE(numout,*) 'p5z_ini :   PISCES biochemical model initialisation'
            WRITE(numout,*) '~~~~~~~     With variable stoichiometry'
         ENDIF
      ENDIF

      neos = 0  !  standard salinity
      ierr =         sms_pisces_alloc()          
      ierr = ierr +  trc_alloc()          ! Start of PISCES-related alloc routines...
      ierr = ierr +  p4z_che_alloc()
      ierr = ierr +  p4z_sink_alloc()
      ierr = ierr +  p4z_opt_alloc()
      ierr = ierr +  p4z_flx_alloc()
      ierr = ierr +  p4z_sed_alloc()
      ierr = ierr +  p4z_lim_alloc()
      IF ( ln_p4z ) THEN
         ierr = ierr +  p4z_prod_alloc()
      ELSE
         ierr = ierr +  p5z_lim_alloc()
         ierr = ierr +  p5z_prod_alloc()
      ENDIF
      ierr = ierr +  p4z_rem_alloc()
      !
      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP in trc_ini_pisces : unable to allocate PISCES arrays' )

      IF( ln_ctl )  CALL prt_ctl_trc_ini

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE          ! masked grid volume
               tmask(ji,jj,jk) = tmask_i(ji,jj)
            END DO
         END DO
      END DO
      !
      CALL p4z_sms_init   ! Main routine
      !
      !                                            ! Time-step
      rfact   = rdt                                ! ---------
      rfactr  = 1. / rfact
      rfact2  = rfact / FLOAT( nrdttrc )
      rfact2r = 1. / rfact2
      xstep  = rfact2 / rday      ! Timestep duration for biology


      IF(lwp) WRITE(numout,*) 
      IF(lwp) WRITE(numout,*) '    Tracer  time step    rfact  = ', rfact, ' rdt = ', rdt
      IF(lwp) write(numout,*) '    Biology time step    rfact2 = ', rfact2
      IF(lwp) WRITE(numout,*)
 
      ! Set biological ratios
      ! ---------------------
      rno3   =   16.   / 122.
      po4r   =   1.e0  / 122.
      o2nit  =  32.    / 122.
      o2ut   = 133.    / 122.
      rdenit  =  ( ( o2ut + o2nit ) * 0.80 - rno3 - rno3 * 0.60 ) / rno3
      rdenita =   3. /  5.
      IF( ln_p5z ) THEN
         no3rat3 = no3rat3 / rno3
         po4rat3 = po4rat3 / po4r
      ENDIF

      CALL p4z_che        ! initialize the chemical constants

      ! Initialization of tracer concentration in case of  no restart 
      !--------------------------------------------------------------
      ln_rsttr = ( nrrec /= 0 ) 
#ifdef NEMO
      IF( .NOT. ln_rsttr ) THEN
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  trn(ji,jj,K,jpdic) = sco2
                  trn(ji,jj,K,jpdoc) = bioma0
                  trn(ji,jj,K,jptal) = alka0
                  trn(ji,jj,K,jpoxy) = oxyg0
                  trn(ji,jj,K,jpcal) = bioma0
                  trn(ji,jj,K,jppo4) = po4 / po4r
                  trn(ji,jj,K,jppoc) = bioma0
                  trn(ji,jj,K,jpgoc) = bioma0
                  trn(ji,jj,K,jpbfe) = bioma0 * 5.e-6
                  trn(ji,jj,K,jpsil) = silic1
                  trn(ji,jj,K,jpgsi) = bioma0 * 0.15
                  trn(ji,jj,K,jpdsi) = bioma0 * 0.15
                  trn(ji,jj,K,jpphy) = bioma0
                  trn(ji,jj,K,jpdia) = bioma0
                  trn(ji,jj,K,jpzoo) = bioma0
                  trn(ji,jj,K,jpmes) = bioma0
                  trn(ji,jj,K,jpfer) = 0.6E-9
                  trn(ji,jj,K,jpsfe) = bioma0 * 5.e-6
                  trn(ji,jj,K,jpdfe) = bioma0 * 2.e-5
                  trn(ji,jj,K,jpnfe) = bioma0 * 2.e-5
                  trn(ji,jj,K,jpnch) = bioma0 * 12. / 70.
                  trn(ji,jj,K,jpdch) = bioma0 * 12. / 70.
                  trn(ji,jj,K,jpno3) = no3
                  trn(ji,jj,K,jpnh4) = bioma0
                  IF( ln_ligand) THEN
                     trn(ji,jj,K,jplgw) = 0.6E-9
                  ENDIF
                  IF( ln_p5z ) THEN
                     trn(ji,jj,K,jpdon) = bioma0
                     trn(ji,jj,K,jpdop) = bioma0
                     trn(ji,jj,K,jppon) = bioma0
                     trn(ji,jj,K,jppop) = bioma0
                     trn(ji,jj,K,jpgon) = bioma0
                     trn(ji,jj,K,jpgop) = bioma0
                     trn(ji,jj,K,jpnph) = bioma0
                     trn(ji,jj,K,jppph) = bioma0
                     trn(ji,jj,K,jppic) = bioma0
                     trn(ji,jj,K,jpnpi) = bioma0
                     trn(ji,jj,K,jpppi) = bioma0
                     trn(ji,jj,K,jpndi) = bioma0
                     trn(ji,jj,K,jppdi) = bioma0
                     trn(ji,jj,K,jppfe) = bioma0 * 5.e-6
                     trn(ji,jj,K,jppch) = bioma0 * 12. / 55.
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF
#endif
      !
      
      ! initialize the half saturation constant for silicate
      ! ----------------------------------------------------
      DO jj = JRANGE
         DO ji = IRANGE
            xksi(ji,jj)    = 2.e-6
            xksimax(ji,jj) = xksi(ji,jj)
            fr_i(ji,jj) = 0.0
         ENDDO
      ENDDO

      IF (ln_p5z) THEN
         sized(:,:,:) = 1.0
         sizen(:,:,:) = 1.0
         sized(:,:,:) = 1.0
      ENDIF

      ! Initialization of chemical variables of the carbon cycle
      ! --------------------------------------------------------
      CALL ahini_for_at(hi)   !  set PH at kt=nit000
      !  
      IF(lwp) THEN               ! control print
         WRITE(numout,*)
         WRITE(numout,*)
         WRITE(numout,*) '          *** Total number of passive tracer jptra = ', jptra
      ENDIF

      CALL tracer_stat( nit000 )

      CALL p4z_opt_init       !  Optic: PAR in the water column
      CALL p4z_lim_init       !  Nutrient limitation
      IF ( ln_p4z ) THEN
         CALL p4z_lim_init
         CALL p4z_prod_init      !  Production initialization
      ELSE
         CALL p5z_lim_init       !  co-limitations by the various nutrients
         CALL p5z_prod_init      !  phytoplankton growth rate over the global ocean.
      ENDIF
      CALL p4z_sbc_init       !  External sources of nutrients
      CALL p4z_fechem_init       !  Iron chemistry
      CALL p4z_rem_init       !  remineralisation
      CALL p4z_poc_init          !  remineralisation of organic particles
      IF( ln_ligand ) &
         & CALL p4z_ligand_init  !  remineralisation of organic ligands

      IF ( ln_p4z ) THEN
         CALL p4z_mort_init      !  Phytoplankton mortality
         CALL p4z_micro_init     !  Microzooplankton
         CALL p4z_meso_init      !  Mesozooplankton
      ELSE
         CALL p5z_mort_init      !  phytoplankton mortality
         CALL p5z_micro_init     !  microzooplankton
         CALL p5z_meso_init      !  mesozooplankton
      ENDIF
      CALL p4z_lys_init       !  Carbonate chemistry
      CALL p4z_flx_init       !  gas exchange
      !
      ! Initialization of the sediment model
      IF( ln_sediment)   CALL sed_init

   END SUBROUTINE trc_ini_pisces


   SUBROUTINE trc_nam_pisces
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_sms_pisces_init  ***
      !!
      !! ** Purpose :   Initialization of PH variable
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      INTEGER  :: jn, ierr
#ifdef key_ligand
      TYPE(PTRACER), DIMENSION(jptra) :: tracer
#else
      TYPE(PTRACER), DIMENSION(jptra+1) :: tracer
#endif
      CHARACTER(LEN=20)::   clname

      NAMELIST/nampistrc/ tracer

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_pisces : read PISCES namelists'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'


      !                               ! Open the namelist file
      !                               ! ----------------------
      clname = 'namelist_pisces'
      CALL ctl_opn( numnatp_ref, TRIM( clname )//'_ref', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )
      CALL ctl_opn( numnatp_cfg, TRIM( clname )//'_cfg', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )

      IF(lwm) CALL ctl_opn( numonp     , 'output.namelist.pis' , &
      &   'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      
      ALLOCATE( ctrcnm(jptra), ctrcnl(jptra), ctrcnu(jptra), STAT = ierr  )  
      IF( ierr /= 0 )   CALL ctl_warn('trc_alloc: failed to allocate arrays')

      IF(lwp) WRITE(numout,*) 'number of tracer : ', jptra
      DO jn = 1, jptra
         WRITE( ctrcnm(jn),'("TR_",I2)'           ) jn
         WRITE( ctrcnl(jn),'("TRACER NUMBER ",I2)') jn
         ctrcnu(jn) = 'mmole/m3'
      END DO

      REWIND( numnatp_ref )
      READ  ( numnatp_ref, nampistrc, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampistrc in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampisrem in configuration namelist : Pisces remineralization
      READ  ( numnatp_cfg, nampistrc, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nampistrc in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, nampistrc )

      DO jn = 1, jptra
         ctrcnm(jn) = tracer(jn)%clsname
         ctrcnl(jn) = tracer(jn)%cllname
         ctrcnu(jn) = tracer(jn)%clunit
      END DO


      IF(lwp) THEN                   ! control print
         DO jn = 1, jptra
            WRITE(numout,*) '   tracer nb             : ', jn 
            WRITE(numout,*) '   short name            : ', TRIM(ctrcnm(jn))
            WRITE(numout,*) '   long name             : ', TRIM(ctrcnl(jn))
            WRITE(numout,*) '   unit                  : ', TRIM(ctrcnu(jn))
            WRITE(numout,*) ' '
         END DO
      ENDIF

   END SUBROUTINE trc_nam_pisces

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                            No PISCES biochemical model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_pisces             ! Empty routine
   END SUBROUTINE trc_ini_pisces
#endif

   !!======================================================================
END MODULE trcini_pisces
