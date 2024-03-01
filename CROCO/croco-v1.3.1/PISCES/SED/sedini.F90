#include "cppdefs.h"

MODULE sedini
   !!======================================================================
   !!              ***  MODULE  sedini  ***
   !! Sediment : define sediment variables
   !!=====================================================================
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   sed_init    : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   !! * Modules used
   USE sed     ! sediment global variable
   USE sed_oce
   USE sedarr
   USE sedadv
   USE sms_pisces, ONLY : rfact

   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC sed_init          ! routine called by opa.F90

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"


   !! Module variables
   REAL(wp)    ::  &
      sedzmin = 0.3    ,  &  !: Minimum vertical spacing
      sedhmax = 10.0   ,  &  !: Maximum depth of the sediment
      sedkth  = 5.0    ,  &  !: Default parameters
      sedacr  = 3.0          !: Default parameters
      
   REAL(wp)    ::  &
      porsurf =  0.95  ,  &  !: Porosity at the surface
      porinf  =  0.75  ,  &  !: Porosity at infinite depth
      rhox    =  2.0         !: Vertical length scale of porosity variation 

   REAL(wp)    ::  &
      rcopal  =   40.  ,  &  !: reactivity for si    [l.mol-1.an-1]
      dcoef   =  8.e-6       !: diffusion coefficient (*por)   [cm**2/s]

   REAL(wp), PUBLIC    ::  &
      redO2    =  172.  ,  &  !: Redfield coef for Oxygen
      redNo3   =   16.  ,  &  !: Redfield coef for Nitrate
      redPo4   =    1.  ,  &  !: Redfield coef for Phosphate
      redC     =  122.  ,  &  !: Redfield coef for Carbon
      redfep   =  0.175 ,  &  !: Ratio for iron bound phosphorus
      rcorgl   =   50.  ,  &  !: reactivity for POC/O2 [l.mol-1.an-1]    
      rcorgs   =   0.5  ,  &  !: reactivity of the semi-labile component
      rcorgr   =   1E-4 ,  &  !: reactivity of the refractory component
      rcnh4    =   10E6 ,  &  !: reactivity for O2/NH4 [l.mol-1.an-1]  
      rch2s    =   1.E5 ,  &  !: reactivity for O2/ODU [l.mol-1.an-1] 
      rcfe2    =   5.E8 ,  &  !: reactivity for O2/Fe2+ [l.mol-1.an-1]
      rcfeh2s  =   1.E4 ,  &  !: Reactivity for FEOH/H2S [l.mol-1.an-1]
      rcfes    =   1.E5 ,  &  !: Reactivity for FE2+/H2S [l.mol-1.an-1]
      rcfeso   =   3.E5 ,  &  !: Reactivity for FES/O2 [l.mol-1.an-1]
      xksedo2  =   5E-6 ,  &  !: half-sturation constant for oxic remin.
      xksedno3 =   5E-6 ,  &  !: half-saturation constant for denitrification
      xksedfeo =   0.6  ,  & !: half-saturation constant for iron remin
      xksedso4 =   2E-3       !: half-saturation constant for SO4 remin

   REAL(wp)    ::  &
      rccal   = 1000.,      & !: reactivity for calcite         [l.mol-1.an-1]
      rcligc  = 1.E-4         !: L/C ratio in POC

   REAL(wp), PUBLIC    ::  dbiot   = 15. , &  !: coefficient for bioturbation    [cm**2.(n-1)]
      dbtbzsc =  10.0  ,    &  !: Vertical scale of variation. If no variation, mixed layer in the sed [cm]
      xirrzsc = 2.0            !: Vertical scale of irrigation variation.
   REAL(wp)    ::  &
      ryear = 365. * 24. * 3600. !:  1 year converted in second

   REAL(wp), DIMENSION(jpwat), PUBLIC  :: diff1s
   DATA diff1s/4.59E-6, 1.104E-5, 4.81E-6 , 9.78E-6, 3.58E-6, 4.01E-6, 9.8E-6, 9.73E-6, 5.0E-6, 3.31E-6 /

   REAL(wp), DIMENSION(jpwat), PUBLIC  :: diff2s
   DATA diff2s/1.74E-7, 4.47E-7, 2.51E-7, 3.89E-7, 1.77E-7, 2.5E-7, 3.89E-7, 3.06E-7, 2.5E-7, 1.5E-7 /


   !! $Id: sedini.F90 10362 2018-11-30 15:38:17Z aumont $
CONTAINS


   SUBROUTINE sed_init
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_init  ***
      !!
      !! ** Purpose :  Initialization of sediment module
      !!               - Reading namelist
      !!               - Read the deepest water layer thickness
      !!                 ( using as mask ) in Netcdf file
      !!               - Convert unity if necessary
      !!               - sets initial sediment composition
      !!                 ( only clay or reading restart file )
      !!               - sets sediment grid, porosity and others constants
      !!
      !!   History :
      !!        !  04-10  (N. Emprin, M. Gehlen )  Original code
      !!        !  06-07  (C. Ethe)  Re-organization
      !!----------------------------------------------------------------------
      INTEGER :: ji, jj, ierr
      !!----------------------------------------------------------------------


      ! Reading namelist.sed variables
      !---------------------------------------

      CALL ctl_opn( numsed, 'sediment.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )

      IF (lwp) THEN
         WRITE(numsed,*)
         WRITE(numsed,*) '                 PISCES framework'
         WRITE(numsed,*) '                 SEDIMENT model'
         WRITE(numsed,*) '                version 3.0  (2018) '
         WRITE(numsed,*)
         WRITE(numsed,*)
      ENDIF

      IF(lwp) WRITE(numsed,*) ' sed_init : Initialization of sediment module  '
      IF(lwp) WRITE(numsed,*) ' '

      ! Read sediment Namelist
      !-------------------------
      CALL sed_init_nam

      ! Allocate SEDIMENT arrays
      ierr =        sed_alloc()
      ierr = ierr + sed_oce_alloc()
      ierr = ierr + sed_adv_alloc() 
      IF( ierr /= 0 )   CALL ctl_stop( 'sed_ini: unable to allocate sediment model arrays' )

      ! Determination of sediments number of points and allocate global variables
      epkbot(:,:) = 0.
      DO jj = JRANGE
         DO ji = IRANGE
            IF( tmask(ji,jj,ikt) == 1 ) epkbot(ji,jj) = e3t_n(ji,jj,KSED)
            gdepbot(ji,jj) = gdepw_n(ji,jj,ikt)
         ENDDO
      ENDDO

      ! computation of total number of ocean points
      !--------------------------------------------
      sedmask = 0.
      IF ( COUNT( epkbot(:,:) > 0. ) == 0 ) THEN 
          sedmask = 0.
      ELSE
          sedmask = 1.
      ENDIF 
      jpoce  = MAX( COUNT( epkbot(:,:) > 0. ) , 1 )

      ! Allocate memory size of global variables
      ALLOCATE( pwcp (jpoce,jpksed,jpwat) )  ;  ALLOCATE( pwcp0 (jpoce,jpksed,jpwat) ) ;  ALLOCATE( pwcp_dta  (jpoce,jpwat) )
      ALLOCATE( solcp(jpoce,jpksed,jpsol) )  ;  ALLOCATE( solcp0(jpoce,jpksed,jpsol) ) ;  ALLOCATE( rainrm_dta(jpoce,jpsol) )
      ALLOCATE( pwcp_avg(jpoce,jpksed,jpwat) )  ; ALLOCATE( solcp_avg(jpoce,jpksed,jpsol) )
      ALLOCATE( rainrm(jpoce,jpsol) )        ;  ALLOCATE( rainrg(jpoce,jpsol) )        ;  ALLOCATE( raintg(jpoce) ) 
      ALLOCATE( dzdep(jpoce) )               ;  ALLOCATE( iarroce(jpoce) )             ;  ALLOCATE( dzkbot(jpoce) )
      ALLOCATE( zkbot(jpoce) )               ;  ALLOCATE( db(jpoce,jpksed) )
      ALLOCATE( temp(jpoce) )                ;  ALLOCATE( salt(jpoce) )  
      ALLOCATE( diff(jpoce,jpksed,jpwat ) )  ;  ALLOCATE( irrig(jpoce, jpksed) )
      ALLOCATE( wacc(jpoce) )                ;  ALLOCATE( fecratio(jpoce) )
      ALLOCATE( press(jpoce) )               ;  ALLOCATE( densSW(jpoce) ) 
      ALLOCATE( hipor(jpoce,jpksed) )        ;  ALLOCATE( co3por(jpoce,jpksed) )
      ALLOCATE( dz3d(jpoce,jpksed) )         ;  ALLOCATE( volw3d(jpoce,jpksed) )       ;  ALLOCATE( vols3d(jpoce,jpksed) )
      ALLOCATE( sedligand(jpoce, jpksed) )
      ALLOCATE( trcsed_avg(jpoce, jpksed, jpdia3dsed) )  ; ALLOCATE( flxsed_avg(jpoce, jpdia2dsed) )

      ! Initialization of global variables
      pwcp  (:,:,:) = 0.   ;  pwcp0 (:,:,:) = 0.  ; pwcp_dta  (:,:) = 0.  
      solcp (:,:,:) = 0.   ;  solcp0(:,:,:) = 0.  ; rainrm_dta(:,:) = 0.
      pwcp_avg(:,:,:) = 0. ;  solcp_avg(:,:,:) = 0.
      rainrm(:,:  ) = 0.   ;  rainrg(:,:  ) = 0.  ; raintg    (:  ) = 0. 
      dzdep (:    ) = 0.   ;  iarroce(:   ) = 0   ; dzkbot    (:  ) = 0.
      temp  (:    ) = 0.   ;  salt   (:   ) = 0.  ; zkbot     (:  ) = 0.
      press (:    ) = 0.   ;  densSW (:   ) = 0.  ; db        (:,:) = 0. 
      hipor (:,:  ) = 0.   ;  co3por (:,: ) = 0.  ; irrig     (:,:) = 0. 
      dz3d  (:,:  ) = 0.   ;  volw3d (:,: ) = 0.  ; vols3d    (:,:) = 0. 
      fecratio(:)   = 1E-5 
      sedligand(:,:) = 0.6E-9
      trcsed_avg(:,:,:) = 0.0  ;  flxsed_avg(:,:) = 0.0

      ! Chemical variables      
      ALLOCATE( akbs  (jpoce) )  ;  ALLOCATE( ak1s   (jpoce) )  ;  ALLOCATE( ak2s  (jpoce) ) ;  ALLOCATE( akws  (jpoce) )     
      ALLOCATE( ak1ps (jpoce) )  ;  ALLOCATE( ak2ps  (jpoce) )  ;  ALLOCATE( ak3ps (jpoce) ) ;  ALLOCATE( aksis (jpoce) )    
      ALLOCATE( aksps (jpoce) )  ;  ALLOCATE( ak12s  (jpoce) )  ;  ALLOCATE( ak12ps(jpoce) ) ;  ALLOCATE( ak123ps(jpoce) )    
      ALLOCATE( borats(jpoce) )  ;  ALLOCATE( calcon2(jpoce) )  ;  ALLOCATE( sieqs (jpoce) ) 
      ALLOCATE( aks3s(jpoce) )   ;  ALLOCATE( akf3s(jpoce) )    ;  ALLOCATE( sulfats(jpoce) )
      ALLOCATE( fluorids(jpoce) )

      akbs  (:) = 0. ;   ak1s   (:) = 0. ;  ak2s  (:) = 0. ;   akws   (:) = 0.
      ak1ps (:) = 0. ;   ak2ps  (:) = 0. ;  ak3ps (:) = 0. ;   aksis  (:) = 0.
      aksps (:) = 0. ;   ak12s  (:) = 0. ;  ak12ps(:) = 0. ;   ak123ps(:) = 0.
      borats(:) = 0. ;   calcon2(:) = 0. ;  sieqs (:) = 0.
      aks3s(:)  = 0. ;   akf3s(:)   = 0. ;  sulfats(:) = 0. ;  fluorids(:) = 0.

      ! Mass balance calculation  
      ALLOCATE( fromsed(jpoce, jpsol) ) ; ALLOCATE( tosed(jpoce, jpsol) ) ;  ALLOCATE( rloss(jpoce, jpsol) )
      ALLOCATE( tokbot (jpoce, jpwat) ) 

      fromsed(:,:) = 0.    ;   tosed(:,:) = 0. ;  rloss(:,:) = 0.  ;   tokbot(:,:) = 0. 

      !  Restart flag
      ncidrstsed = -1
      ncidwrised = -1
      nrecsedpis_avg = 0
      nrecsedpis     = 0

      ! Initialization of sediment geometry
      !------------------------------------
      CALL sed_init_geom

   END SUBROUTINE sed_init

   SUBROUTINE sed_init_geom
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_init_geom  ***
      !!
      !! ** Purpose :  Initialization of sediment geometry
      !!               - Read the deepest water layer thickness
      !!                 ( using as mask ) in Netcdf file
      !!               - sets sediment grid, porosity and molecular weight
      !!                 and others constants
      !!
      !!   History :
      !!        !  06-07  (C. Ethe)  Original
      !!----------------------------------------------------------------------
      !! * Modules used
      !! * local declarations
      INTEGER  :: ji, jj, jk, jn
      REAL(wp) :: za0, za1, zt, zw, zsum, zsur, zprof, zprofw
      REAL(wp) :: ztmp1, ztmp2
      !---------------------------------------------------------- 

      IF(lwp) WRITE(numsed,*) ' sed_init_geom : Initialization of sediment geometry '
      IF(lwp) WRITE(numsed,*) ' '

      ! Computation of 1D array of sediments points
      indoce = 0
      DO jj = JRANGE
         DO ji = IRANGE
            IF (  epkbot(ji,jj) > 0. ) THEN
               indoce          = indoce + 1
               iarroce(indoce) = (jj - 1) * jpi + ji
            ENDIF
         END DO
      END DO

      IF ( indoce .EQ. 0 ) THEN
         indoce = 1
         iarroce(indoce) = 1
      ENDIF

      IF( indoce .NE. jpoce ) THEN
         CALL ctl_stop( 'sed_ini: number of ocean points indoce doesn''t match  number of point' )
      ELSE
         IF (lwp) WRITE(numsed,*) ' '
         IF (lwp) WRITE(numsed,*) ' total number of ocean points jpoce =  ',jpoce
         IF (lwp) WRITE(numsed,*) ' '
      ENDIF

      ! initialization of dzkbot in [cm]
      !------------------------------------------------    
      CALL pack_arr ( jpoce, dzkbot(1:jpoce), epkbot(PRIV_2D_BIOARRAY), iarroce(1:jpoce) )
      dzkbot(1:jpoce) = dzkbot(1:jpoce) * 1.e+2 
      CALL pack_arr ( jpoce, zkbot(1:jpoce), gdepbot(PRIV_2D_BIOARRAY), iarroce(1:jpoce) )

      ! Geometry and  constants 
      ! sediment layer thickness [cm]
      ! (1st layer= diffusive layer = pur water) 
      !------------------------------------------
      za1  = (  sedzmin - sedhmax / FLOAT(jpksed-1)  )                                                      &
         & / ( TANH((1-sedkth)/sedacr) - sedacr/FLOAT(jpksed-1) * (  LOG( COSH( (jpksed - sedkth) / sedacr) )      &
         &                                                   - LOG( COSH( ( 1  - sedkth) / sedacr) )  )  )
      za0  = sedzmin - za1 * TANH( (1-sedkth) / sedacr )
      zsur = - za0 - za1 * sedacr * LOG( COSH( (1-sedkth) / sedacr )  )

      profsedw(1) = 0.0
      profsed(1) = -dz(1) / 2.
      DO jk = 2, jpksed
         zw = REAL( jk , wp )
         zt = REAL( jk , wp ) - 0.5
         profsed(jk)  = ( zsur + za0 * zt + za1 * sedacr * LOG ( COSH( (zt-sedkth) / sedacr ) )  ) 
         profsedw(jk) = ( zsur + za0 * zw + za1 * sedacr * LOG ( COSH( (zw-sedkth) / sedacr ) )  )
      END DO

      dz(1) = 0.1
      DO jk = 2, jpksed
         dz(jk) = profsedw(jk) - profsedw(jk-1)
      END DO

      DO jk = 1, jpksed
         DO ji = 1, jpoce
            dz3d(ji,jk) = dz(jk)
         END DO
      ENDDO

      !  Porosity profile [0]
      !---------------------
      por(1) = 1.0
      DO jk = 2, jpksed
         por(jk) = porinf + ( porsurf-porinf) * exp(-rhox * (profsed(jk) ) )
      END DO
 
      ! inverse of  Porosity profile
      !-----------------------------
      por1(:) = 1. - por(:)

      ! Volumes of pore water and solid fractions (vector and array)
      !     WARNING : volw(1) and vols(1) are sublayer volums
      volw(:) = dz(:) * por(:)
      vols(:) = dz(:) * por1(:)

      ! temporary new value for dz3d(:,1) 
      dz3d(1:jpoce,1) = dzkbot(1:jpoce)

      ! WARNING : volw3d(:,1) and vols3d(:,1) are deepest water column volums
      DO jk = 1, jpksed
         volw3d(1:jpoce,jk) = dz3d(1:jpoce,jk) * por (jk)
         vols3d(1:jpoce,jk) = dz3d(1:jpoce,jk) * por1(jk)
      ENDDO

      ! Back to the old sublayer vlaue for dz3d(:,1)
      dz3d(1:jpoce,1) = dz(1)

      !---------------------------------------------
      ! Molecular weight [g/mol] for solid species
      !---------------------------------------------

      ! opal=sio2*0.4(h20)=28+2*16+0.4*(2+16)
      !---------------------------------------
      mol_wgt(jsopal) = 28. + 2. * 16. + 0.4 * ( 2. + 16. )  

      !  clay
      !  some kind of Illit (according to Pape)
      !  K0.58(Al 1.38 Fe(III)0.37Fe(II)0.04Mg0.34)[(OH)2|(Si3.41Al0.59)O10]
      !--------------------------------------------------------------------
      mol_wgt(jsclay) = 0.58 * 39. + 1.38 * 27. + ( 0.37 + 0.04 ) * 56.+ &
         &              0.34 * 24. + 2. * ( 16. + 1. ) + 3.41 * 38. +    &
         &              0.59 * 27. + 10. * 16.

      mol_wgt(jsfeo)  = 55.0 + 3.0 * ( 16.0 + 1.0)

      mol_wgt(jsfes)  = 55.0 + 32.0

      ! for chemistry Poc : C(122)H(244)O(86)N(16)P(1)
      ! But den sity of Poc is an Hydrated material (= POC + 30H2O)
      ! So C(122)H(355)O(120)N(16)P(1)
      !------------------------------------------------------------
      mol_wgt(jspoc) = ( 122. * 12. + 355. + 120. * 16.+ &
         &                16. * 14. + 31. ) / 122.
      mol_wgt(jspos) = mol_wgt(jspoc)
      mol_wgt(jspor) = mol_wgt(jspoc)

      ! CaCO3
      !---------
      mol_wgt(jscal) = 40. + 12. + 3. * 16.

      ! Density of solid material in sediment [g/cm**3]
      !------------------------------------------------
      denssol = 2.6

      ! Initialization of diffusion coefficient as function of porosity [cm**2/s]
      !--------------------------------------------------------------------
!      DO jn = 1, jpsol
!         DO jk = 1, jpksed
!            DO ji = 1, jpoce
!               diff(ji,jk,jn) = dcoef / ( 1.0 - 2.0 * log(por(jk)) )
!            END DO
!         END DO
!      END DO

      ! Accumulation rate from Burwicz et al. (2011). This is used to
      ! compute the flux of clays and minerals
      ! --------------------------------------------------------------
      DO ji = 1, jpoce
          ztmp1 = 0.117 / ( 1.0 + ( zkbot(ji) / 200.)**3 )
          ztmp2 = 0.006 / ( 1.0 + ( zkbot(ji) / 4000.)**10 )
          wacc(ji) = ztmp1 + ztmp2
      END DO


      ! Initialization of time step as function of porosity [cm**2/s]
      !------------------------------------------------------------------
   END SUBROUTINE sed_init_geom

   SUBROUTINE sed_init_nam
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_init_nam  ***
      !!
      !! ** Purpose :  Initialization of sediment geometry
      !!               - Reading namelist and defines constants variables
      !!
      !!   History :
      !!        !  06-07  (C. Ethe)  Original
      !!----------------------------------------------------------------------

      INTEGER ::   numnamsed_ref = -1           !! Logical units for namelist sediment
      INTEGER ::   numnamsed_cfg = -1           !! Logical units for namelist sediment
      INTEGER :: ios                 ! Local integer output status for namelist read
      CHARACTER(LEN=20)   ::   clname

      TYPE PSED
         CHARACTER(len = 20)  :: snamesed   !: short name
         CHARACTER(len = 80 ) :: lnamesed   !: long name
         CHARACTER(len = 20 ) :: unitsed    !: unit
      END TYPE PSED

      TYPE(PSED) , DIMENSION(jpsol     ) :: sedsol
      TYPE(PSED) , DIMENSION(jpwat     ) :: sedwat
      TYPE(PSED) , DIMENSION(jpdia3dsed) :: seddiag3d
      TYPE(PSED) , DIMENSION(jpdia2dsed) :: seddiag2d

      NAMELIST/nam_run/nrseddt,ln_sed_2way
      NAMELIST/nam_geom/jpksed, sedzmin, sedhmax, sedkth, sedacr, porsurf, porinf, rhox
      NAMELIST/nam_trased/sedsol, sedwat
      NAMELIST/nam_diased/seddiag3d, seddiag2d
      NAMELIST/nam_inorg/rcopal, dcoef, rccal, ratligc, rcligc
      NAMELIST/nam_poc/redO2, redNo3, redPo4, redC, redfep, rcorgl, rcorgs,  &
         &             rcorgr, rcnh4, rch2s, rcfe2, rcfeh2s, rcfes, rcfeso,  &
         &             xksedo2, xksedno3, xksedfeo, xksedso4
      NAMELIST/nam_btb/dbiot, ln_btbz, dbtbzsc, adsnh4, ln_irrig, xirrzsc
      NAMELIST/nam_rst/ln_rst_sed, cn_sedrst_indir, cn_sedrst_outdir, cn_sedrst_in, cn_sedrst_out
      NAMELIST/nam_output/ldefsedpis_avg,cn_sedwri_out, nrpfsedpis_avg, nwrtsedpis_avg, ntssedpis_avg

      INTEGER :: ji, jn, jn1
      !-------------------------------------------------------

      IF(lwp) WRITE(numsed,*) ' sed_init_nam : Read namelists '
      IF(lwp) WRITE(numsed,*) ' '

      ! ryear = 1 year converted in second
      !------------------------------------
      IF (lwp) THEN
         WRITE(numsed,*) ' '
         WRITE(numsed,*) 'number of seconds in one year : ryear = ', ryear
         WRITE(numsed,*) ' '     
      ENDIF

      ! Reading namelist.sed variables
      !---------------------------------
      clname = 'namelist_sediment'
      IF(lwp) WRITE(numsed,*) ' sed_init_nam : read SEDIMENT namelist'
      IF(lwp) WRITE(numsed,*) ' ~~~~~~~~~~~~~~'
      CALL ctl_opn( numnamsed_ref, TRIM( clname )//'_ref', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      CALL ctl_opn( numnamsed_cfg, TRIM( clname )//'_cfg', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )

      nitsed000 = nittrc000
      nitsedend = nitend
      ! Namelist nam_run
      REWIND( numnamsed_ref )              ! Namelist nam_run in reference namelist : Pisces variables
      READ  ( numnamsed_ref, nam_run, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_run in reference namelist', lwp )

      REWIND( numnamsed_cfg )              ! Namelist nam_run in reference namelist : Pisces variables
      READ  ( numnamsed_cfg, nam_run, IOSTAT = ios, ERR = 902)
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_run in configuration namelist', lwp )

      IF (lwp) THEN
         WRITE(numsed,*) ' namelist nam_run'
         WRITE(numsed,*) ' Nb of iterations for fast species    nrseddt = ', nrseddt
         WRITE(numsed,*) ' 2-way coupling between PISCES and Sed ln_sed_2way = ', ln_sed_2way
      ENDIF

      IF ( ln_p5z .AND. ln_sed_2way ) CALL ctl_stop( '2 ways coupling with sediment cannot be activated with PISCES-QUOTA' )

      REWIND( numnamsed_ref )              ! Namelist nam_geom in reference namelist : Pisces variables
      READ  ( numnamsed_ref, nam_geom, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_geom in reference namelist', lwp )

      REWIND( numnamsed_cfg )              ! Namelist nam_geom in reference namelist : Pisces variables
      READ  ( numnamsed_cfg, nam_geom, IOSTAT = ios, ERR = 904)
904   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_geom in configuration namelist', lwp )

      IF (lwp) THEN 
         WRITE(numsed,*) ' namelist nam_geom'
         WRITE(numsed,*) ' Number of vertical layers            jpksed  = ', jpksed
         WRITE(numsed,*) ' Minimum vertical spacing             sedzmin = ', sedzmin
         WRITE(numsed,*) ' Maximum depth of the sediment        sedhmax = ', sedhmax
         WRITE(numsed,*) ' Default parameter                    sedkth  = ', sedkth
         WRITE(numsed,*) ' Default parameter                    sedacr  = ', sedacr
         WRITE(numsed,*) ' Sediment porosity at the surface     porsurf = ', porsurf
         WRITE(numsed,*) ' Sediment porosity at infinite depth  porinf  = ', porinf
         WRITE(numsed,*) ' Length scale of porosity variation   rhox    = ', rhox
      ENDIF

      jpksedm1  = jpksed - 1
      dtsed = rfact

      REWIND( numnamsed_ref )              ! Namelist nam_trased in reference namelist : Pisces variables
      READ  ( numnamsed_ref, nam_trased, IOSTAT = ios, ERR = 905)
905   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_trased in reference namelist', lwp )

      REWIND( numnamsed_cfg )              ! Namelist nam_trased in reference namelist : Pisces variables
      READ  ( numnamsed_cfg, nam_trased, IOSTAT = ios, ERR = 906)
906   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_trased in configuration namelist', lwp )

      DO jn = 1, jpsol
         sedtrcd(jn) = sedsol(jn)%snamesed
         sedtrcl(jn) = sedsol(jn)%lnamesed
         sedtrcu(jn) = sedsol(jn)%unitsed
      END DO

      DO jn = 1, jpwat
         jn1 = jn + jpsol
         sedtrcd(jn1) = sedwat(jn)%snamesed
         sedtrcl(jn1) = sedwat(jn)%lnamesed
         sedtrcu(jn1) = sedwat(jn)%unitsed
      END DO

      IF (lwp) THEN
         WRITE(numsed,*) ' namelist nam_trased'
         WRITE(numsed,*) ' '
         DO jn = 1, jptrased
            WRITE(numsed,*) 'name of 3d output sediment field number :',jn,' : ',TRIM(sedtrcd(jn))
            WRITE(numsed,*) 'long name ', TRIM(sedtrcl(jn))
            WRITE(numsed,*) ' in unit = ', TRIM(sedtrcu(jn))
            WRITE(numsed,*) ' '
         END DO
         WRITE(numsed,*) ' '
      ENDIF

      REWIND( numnamsed_ref )              ! Namelist nam_diased in reference namelist : Pisces variables
      READ  ( numnamsed_ref, nam_diased, IOSTAT = ios, ERR = 907)
907   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diased in reference namelist', lwp )

      REWIND( numnamsed_cfg )              ! Namelist nam_diased in reference namelist : Pisces variables
      READ  ( numnamsed_cfg, nam_diased, IOSTAT = ios, ERR = 908)
908   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diased in configuration namelist', lwp )
      
      DO jn = 1, jpdia3dsed
         seddia3d(jn) = seddiag3d(jn)%snamesed
         seddia3l(jn) = seddiag3d(jn)%lnamesed
         seddia3u(jn) = seddiag3d(jn)%unitsed
      END DO

      DO jn = 1, jpdia2dsed
         seddia2d(jn) = seddiag2d(jn)%snamesed
         seddia2l(jn) = seddiag2d(jn)%lnamesed
         seddia2u(jn) = seddiag2d(jn)%unitsed
      END DO

      IF (lwp) THEN
         WRITE(numsed,*) ' namelist nam_diased'
         WRITE(numsed,*) ' '
         DO jn = 1, jpdia3dsed
            WRITE(numsed,*) 'name of 3D output diag number :',jn, ' : ', TRIM(seddia3d(jn))
            WRITE(numsed,*) 'long name ', TRIM(seddia3l(jn))
            WRITE(numsed,*) ' in unit = ',TRIM(seddia3u(jn))
            WRITE(numsed,*) ' '
         END DO

         DO jn = 1, jpdia2dsed
            WRITE(numsed,*) 'name of 2D output diag number :',jn, ' : ', TRIM(seddia2d(jn))
            WRITE(numsed,*) 'long name ', TRIM(seddia2l(jn))
            WRITE(numsed,*) ' in unit = ',TRIM(seddia2u(jn))
            WRITE(numsed,*) ' '
         END DO

         WRITE(numsed,*) ' '
      ENDIF

      !
      ! Additional biogenetic arrays
      !
        sname(1,1) = TRIM(seddia3d(1))
        sname(1,2) = TRIM(seddia3l(1))
        sname(1,3) = TRIM(seddia3u(1))
        sname(2,1) = TRIM(seddia3d(2))
        sname(2,2) = TRIM(seddia3l(2))
        sname(2,3) = TRIM(seddia3u(2))
        sname(3,1) = "dbioturb"
        sname(3,2) = "Bioturbation coefficient"
        sname(3,3) = "cm2/yr"
        sname(4,1) =  "irrig"
        sname(4,2) = "Irrigation coefficient"
        sname(4,3) = "-"
        sname(5,1) = "sedligand"
        sname(5,2) = "Fe ligands concentration"
        sname(5,3) = "mol/l"


      ! Inorganic chemistry parameters
      !----------------------------------
      REWIND( numnamsed_ref )              ! Namelist nam_inorg in reference namelist : Pisces variables
      READ  ( numnamsed_ref, nam_inorg, IOSTAT = ios, ERR = 909)
909   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_inorg in reference namelist', lwp )

      REWIND( numnamsed_cfg )              ! Namelist nam_inorg in reference namelist : Pisces variables
      READ  ( numnamsed_cfg, nam_inorg, IOSTAT = ios, ERR = 910)
910   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_inorg in configuration namelist', lwp )

      IF (lwp) THEN
         WRITE(numsed,*) ' namelist nam_inorg'
         WRITE(numsed,*) ' reactivity for Si      rcopal  = ', rcopal
         WRITE(numsed,*) ' diff. coef for por.    dcoef   = ', dcoef
         WRITE(numsed,*) ' reactivity for calcite rccal   = ', rccal
         WRITE(numsed,*) ' L/C ratio in POC       ratligc = ', ratligc
         WRITE(numsed,*) ' reactivity for ligands rcligc  = ', rcligc
         WRITE(numsed,*) ' '
      ENDIF

      ! Unity conversion to get saturation conc. psat in [mol.l-1]
      ! and reactivity rc in  [l.mol-1.s-1]
      !----------------------------------------------------------
      reac_sil   = rcopal / ryear     
      reac_ligc  = rcligc / ryear

      ! Additional parameter linked to POC/O2/No3/Po4
      !----------------------------------------------
      REWIND( numnamsed_ref )              ! Namelist nam_poc in reference namelist : Pisces variables
      READ  ( numnamsed_ref, nam_poc, IOSTAT = ios, ERR = 911)
911   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_poc in reference namelist', lwp )

      REWIND( numnamsed_cfg )              ! Namelist nam_poc in reference namelist : Pisces variables
      READ  ( numnamsed_cfg, nam_poc, IOSTAT = ios, ERR = 912)
912   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_poc in configuration namelist', lwp )

      IF (lwp) THEN
         WRITE(numsed,*) ' namelist nam_poc'
         WRITE(numsed,*) ' Redfield coef for oxy            redO2    = ', redO2
         WRITE(numsed,*) ' Redfield coef for no3            redNo3   = ', redNo3
         WRITE(numsed,*) ' Redfield coef for po4            redPo4   = ', redPo4
         WRITE(numsed,*) ' Redfield coef for carbon         redC     = ', redC
         WRITE(numsed,*) ' Ration for iron bound P          redfep   = ', redfep
         WRITE(numsed,*) ' reactivity for labile POC        rcorgl   = ', rcorgl
         WRITE(numsed,*) ' reactivity for semi-refract. POC rcorgs   = ', rcorgs
         WRITE(numsed,*) ' reactivity for refractory POC    rcorgr   = ', rcorgr
         WRITE(numsed,*) ' reactivity for NH4               rcnh4    = ', rcnh4
         WRITE(numsed,*) ' reactivity for H2S               rch2s    = ', rch2s
         WRITE(numsed,*) ' reactivity for Fe2+              rcfe2    = ', rcfe2
         WRITE(numsed,*) ' reactivity for FeOH/H2S          rcfeh2s  = ', rcfeh2s
         WRITE(numsed,*) ' reactivity for Fe2+/H2S          rcfes    = ', rcfes
         WRITE(numsed,*) ' reactivity for FeS/O2            rcfeso   = ', rcfeso
         WRITE(numsed,*) ' Half-sat. cste for oxic remin    xksedo2  = ', xksedo2
         WRITE(numsed,*) ' Half-sat. cste for denit.        xksedno3 = ', xksedno3
         WRITE(numsed,*) ' Half-sat. cste for iron remin    xksedfeo = ', xksedfeo
         WRITE(numsed,*) ' Half-sat. cste for SO4 remin     xksedso4 = ', xksedso4
         WRITE(numsed,*) ' '
      ENDIF


      so2ut  = redO2    / redC
      srno3  = redNo3   / redC
      spo4r  = redPo4   / redC
      srDnit = ( (redO2 + 32. ) * 0.8 - redNo3 - redNo3 * 0.6 ) / redC
      ! reactivity rc in  [l.mol-1.s-1]
      reac_pocl  = rcorgl / ryear
      reac_pocs  = rcorgs / ryear
      reac_pocr  = rcorgr / ryear
      reac_nh4   = rcnh4  / ryear
      reac_h2s   = rch2s  / ryear
      reac_fe2   = rcfe2  / ryear
      reac_feh2s = rcfeh2s/ ryear
      reac_fes   = rcfes  / ryear
      reac_feso  = rcfeso / ryear

      ! reactivity rc in  [l.mol-1.s-1]      
      reac_cal = rccal / ryear

      ! Bioturbation parameter
      !------------------------
      REWIND( numnamsed_ref )              ! Namelist nam_btb in reference namelist : Pisces variables
      READ  ( numnamsed_ref, nam_btb, IOSTAT = ios, ERR = 913)
913   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_btb in reference namelist', lwp )

      REWIND( numnamsed_cfg )              ! Namelist nam_btb in reference namelist : Pisces variables
      READ  ( numnamsed_cfg, nam_btb, IOSTAT = ios, ERR = 914)
914   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_btb in configuration namelist', lwp )

      IF (lwp) THEN
         WRITE(numsed,*) ' namelist nam_btb ' 
         WRITE(numsed,*) ' coefficient for bioturbation      dbiot    = ', dbiot
         WRITE(numsed,*) ' Depth varying bioturbation        ln_btbz  = ', ln_btbz
         WRITE(numsed,*) ' coefficient for btb attenuation   dbtbzsc  = ', dbtbzsc
         WRITE(numsed,*) ' Adsorption coefficient of NH4     adsnh4   = ', adsnh4
         WRITE(numsed,*) ' Bioirrigation in sediment         ln_irrig = ', ln_irrig
         WRITE(numsed,*) ' coefficient for irrig attenuation xirrzsc  = ', xirrzsc
         WRITE(numsed,*) ' '
      ENDIF

      ! Initial value (t=0) for sediment pore water and solid components
      !----------------------------------------------------------------
      REWIND( numnamsed_ref )              ! Namelist nam_rst in reference namelist : Pisces variables
      READ  ( numnamsed_ref, nam_rst, IOSTAT = ios, ERR = 915)
915   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_rst in reference namelist', lwp )

      REWIND( numnamsed_cfg )              ! Namelist nam_rst in reference namelist : Pisces variables
      READ  ( numnamsed_cfg, nam_rst, IOSTAT = ios, ERR = 916)
916   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_rst in configuration namelist', lwp )

      IF (lwp) THEN
         WRITE(numsed,*) ' namelist  nam_rst ' 
         WRITE(numsed,*) '  boolean term for restart (T or F) ln_rst_sed = ', ln_rst_sed 
         WRITE(numsed,*) ' '
      ENDIF

      ! Namelist corresponding to output
      !---------------------------------
      REWIND( numnamsed_ref )              ! Namelist nam_output in reference namelist : Pisces variables
      READ  ( numnamsed_ref, nam_output, IOSTAT = ios, ERR = 917)
917   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_output in reference namelist', lwp )

      REWIND( numnamsed_cfg )              ! Namelist nam_output in reference namelist : Pisces variables
      READ  ( numnamsed_cfg, nam_output, IOSTAT = ios, ERR = 918)
918   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_output in configuration namelist', lwp )

      IF (lwp) THEN
         WRITE(numsed,*) ' namelist  nam_output '
         WRITE(numsed,*) '  Write averaged output variables     ldefsedpis_avg = ', ldefsedpis_avg
         WRITE(numsed,*) '  Name of the ooutput file            cn_sedwri_out = ', cn_sedwri_out
         WRITE(numsed,*) '  Frequency of the averaged outputs   nfrecsedpis_avg = ', nrpfsedpis_avg
         WRITE(numsed,*) '  Frequency of the averaged outputs   nwrtsedpis_avg = ', nwrtsedpis_avg
         WRITE(numsed,*) '  Frequency of the averaged outputs   ntssedpis_avg = ', ntssedpis_avg
         WRITE(numsed,*) ' '
      ENDIF


      CLOSE( numnamsed_cfg )
      CLOSE( numnamsed_ref )

   END SUBROUTINE sed_init_nam

#endif

END MODULE sedini
