#include "cppdefs.h"

MODULE p4zsbc
   !!======================================================================
   !!                         ***  MODULE p4sbc  ***
   !! TOP :   PISCES surface boundary conditions of external inputs of nutrients
   !!======================================================================
   !! History :   3.5  !  2012-07 (O. Aumont, C. Ethe) Original code
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   p4z_sbc        :  Read and interpolate time-varying nutrients fluxes
   !!   p4z_sbc_init   :  Initialization of p4z_sbc
   !!----------------------------------------------------------------------
   USE sms_pisces      !  PISCES Source Minus Sink variables
!   USE iom             !  I/O manager
!   USE fldread         !  time interpolation
#ifdef AGRIF
      USE param, ONLY : Lmmpi,Mmmpi
!      USE trcini_pisces
#endif

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sbc
   PUBLIC   p4z_sbc_init   

   !! * Module variables
   LOGICAL , PUBLIC ::   ln_dust      !: boolean for dust input from the atmosphere
   LOGICAL , PUBLIC ::   ln_river     !: boolean for river input of nutrients
   LOGICAL , PUBLIC ::   ln_ndepo     !: boolean for atmospheric deposition of N
   LOGICAL , PUBLIC ::   ln_ironsed   !: boolean for Fe input from sediments
   REAL(wp), PUBLIC ::   sedfeinput   !: Coastal release of Iron
   REAL(wp), PUBLIC ::   dustsolub    !: Solubility of the dust
   REAL(wp), PUBLIC ::   mfrac        !: Mineral Content of the dust
   REAL(wp), PUBLIC ::   wdust        !: Sinking speed of the dust 
   REAL(wp), PUBLIC ::   nitrfix      !: Nitrogen fixation rate   
   REAL(wp), PUBLIC ::   diazolight   !: Nitrogen fixation sensitivty to light 
   REAL(wp), PUBLIC ::   concfediaz   !: Fe half-saturation Cste for diazotrophs 
   LOGICAL , PUBLIC ::   ll_sbc

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   dust             !: dust fields
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   no3dep, nh4dep   !: N depo fields
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   rivdic, rivalk   !: river input fields
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   rivdin, rivdip   !: river input fields
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   rivdon, rivdop   !: river input fields
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   rivdoc           !: river input fields
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   rivdsi           !: river input fields
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ironsed          !: Coastal supply of iron
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   dustmo, no3depmo, nh4depmo !: 2 consecutive set of dust fields 

   REAL(wp), PUBLIC :: sedsilfrac, sedcalfrac
   REAL(wp), PUBLIC :: year2daydta

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zsbc.F90 10868 2019-04-15 10:32:56Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_sbc( kt )
      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_sbc  ***
      !!
      !! ** purpose :   read and interpolate the external sources of nutrients
      !!
      !! ** method  :   read the files and interpolate the appropriate variables
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !
      INTEGER :: ji, jj, jk
      INTEGER, PARAMETER :: jpmois = 12
      INTEGER :: irec1, irec2, i15
      INTEGER :: nyear, nday, nmonth
      REAL    :: zpdtan, zpdtmo, zdemi, zt
      REAL    :: zxy, zjulian, zsec
      !!---------------------------------------------------------------------
      !
      ! Compute dust at nit000 or only if there is more than 1 time record in dust file
      IF( kt == nit000 .AND. lwp ) THEN
        WRITE(numout,*) ' '
        WRITE(numout,*) ' Number of days per year in file year2daydta = ', year2daydta
        WRITE(numout,*) ' '
      ENDIF
      !
      zpdtan = ( year2daydta * day2sec ) / rdt
      zpdtmo = zpdtan / float( jpmois )
      zdemi  = zpdtmo / 2.
      zt     = ( float( kt ) + zdemi) / zpdtmo


      !  recherche de l'indice des enregistrements
      !  du modele dynamique encadrant le pas de temps kt.
      !  --------------------------------------------------
      irec1 = int( zt )
      irec2 = irec1 + 1
      irec1 = MOD( irec1, jpmois )
      IF ( irec1 == 0 ) irec1 = jpmois
      irec2 = MOD( irec2, jpmois )
      IF ( irec2 == 0 ) irec2 = jpmois

      zxy = zt - float(int ( zt ) )
      !
      IF( ln_dust ) THEN
         DO jj = JRANGE
            DO ji = IRANGE
               dust(ji,jj) = ( 1. - zxy ) * dustmo(ji,jj,irec1) + zxy   * dustmo(ji,jj,irec2)
            END DO
          END DO
      ENDIF
      !
      IF( ln_ndepo ) THEN
         DO jj = JRANGE
            DO ji = IRANGE
               no3dep(ji,jj) = ( 1. - zxy ) * no3depmo(ji,jj,irec1) + zxy   * no3depmo(ji,jj,irec2)
               !
               nh4dep(ji,jj) = ( 1. - zxy ) * nh4depmo(ji,jj,irec1) + zxy  * nh4depmo(ji,jj,irec2)
            END DO
          END DO
      ENDIF

   END SUBROUTINE p4z_sbc


   SUBROUTINE p4z_sbc_init
      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_sbc_init  ***
      !!
      !! ** purpose :   initialization of the external sources of nutrients
      !!
      !! ** method  :   read the files and compute the budget
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
# include "netcdf.inc"
      INTEGER  :: ji, jj, jk, irec
      INTEGER :: ncid, varid, dimid, ierr, lstr, lenstr, nf_fread, nrec_dust
      INTEGER :: vartype, nvatts, latt, nvdims
      INTEGER :: vdims(5)
      INTEGER  :: ios                 ! Local integer output status for namelist read
      CHARACTER(len=16) :: varname, dimname, attname
      REAL(wp) :: zexpide, zdenitide, zmaskt

      REAL     ::  cycle_length
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  dustmp
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  no3deptmp, nh4deptmp
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  zcmask
#ifdef MPI
#define LOCALLM Lmmpi
#define LOCALMM Mmmpi
#else
#define LOCALLM Lm
#define LOCALMM Mm
#endif
      !
      NAMELIST/nampissbc/ ln_dust, ln_river, ln_ndepo, ln_ironsed, sedfeinput,   &
      &                   dustsolub, wdust, mfrac, nitrfix, diazolight, concfediaz 
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_sbc_init : initialization of the external sources of nutrients '
         WRITE(numout,*) '~~~~~~~~~~~~ '
      ENDIF
      !                            !* set file information
      REWIND( numnatp_ref )              ! Namelist nampissbc in reference namelist : Pisces external sources of nutrients
      READ  ( numnatp_ref, nampissbc, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampissbc in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampissbc in configuration namelist : Pisces external sources of nutrients
      READ  ( numnatp_cfg, nampissbc, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nampissbc in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampissbc )

      IF(lwp) THEN
         WRITE(numout,*) '   Namelist : nampissbc '
         WRITE(numout,*) '      dust input from the atmosphere           ln_dust     = ', ln_dust
         WRITE(numout,*) '      river input of nutrients                 ln_river    = ', ln_river
         WRITE(numout,*) '      atmospheric deposition of n              ln_ndepo    = ', ln_ndepo
         WRITE(numout,*) '      Fe input from sediments                  ln_ironsed  = ', ln_ironsed
         WRITE(numout,*) '      coastal release of iron                  sedfeinput  = ', sedfeinput
         WRITE(numout,*) '      solubility of the dust                   dustsolub   = ', dustsolub
         WRITE(numout,*) '      Mineral Fe content of the dust           mfrac       = ', mfrac
         WRITE(numout,*) '      sinking speed of the dust                wdust       = ', wdust
         WRITE(numout,*) '      nitrogen fixation rate                   nitrfix     = ', nitrfix
         WRITE(numout,*) '      nitrogen fixation sensitivty to light    diazolight  = ', diazolight
         WRITE(numout,*) '      Fe half-saturation cste for diazotrophs  concfediaz  = ', concfediaz
      END IF

      IF( ln_dust .OR. ln_river .OR. ln_ndepo ) THEN   ;   ll_sbc = .TRUE.
      ELSE                                             ;   ll_sbc = .FALSE.
      ENDIF

      ! coastal and island masks
      ! ------------------------
      IF( ln_ironsed ) THEN     
         ALLOCATE( zcmask(PRIV_3D_BIOARRAY), ironsed(PRIV_3D_BIOARRAY) )

         zcmask(:,:,:) = 0.0
         DO jj = JRANGE
            DO ji = IRANGE
               zcmask(ji,jj,jpk) = 1
            ENDDO
         ENDDO

         DO jk = 1, N
            DO jj = JRANGE
               DO ji = IRANGE
                  IF( tmask_i(ji,jj) /= 0. ) THEN
                     zmaskt = tmask_i(ji+1,jj) * tmask_i(ji-1,jj  ) * tmask_i(ji,jj+1)    &
                        &                      * tmask_i(ji  ,jj-1) * tmask_i(ji,jj  )
                     IF( zmaskt == 0. )   zcmask(ji,jj,jk ) = 0.1
                  ENDIF
               END DO
            END DO
         END DO
         !
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  zexpide   = MIN( 8.,( gdept_n(ji,jj,K) / 500. )**(-1.5) )
                  zdenitide = -0.9543 + 0.7662 * LOG( zexpide ) - 0.235 * LOG( zexpide )**2
                  zcmask(ji,jj,jk) = zcmask(ji,jj,jk) * MIN( 1., EXP( zdenitide ) / 0.5 )
               END DO
            END DO
         END DO
         ! Coastal supply of iron
         ! -------------------------
         ironsed(:,:,jpk) = 0.
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  ironsed(ji,jj,jk) = sedfeinput * zcmask(ji,jj,jk) / ( e3t_n(ji,jj,K) * rday )
               END DO
            END DO
         END DO
         DEALLOCATE( zcmask)
      ENDIF
      !
      !
!    READ DUST INPUT FROM ATMOSPHERE
!    -------------------------------------
!
      IF( ln_dust .OR. ln_ndepo ) THEN
         lstr = lenstr(bioname)
         ierr = nf_open (bioname(1:lstr), nf_nowrite, ncid)
         IF (ierr .NE. nf_noerr .AND. lwp ) THEN
            WRITE(numout,4) bioname
         ENDIF
         ierr = nf_inq_varid(ncid,"dust_time",varid)
! bug if compilation with gfortran
!         ierr =nf_inq_var (ncid, varid, varname, vartype, nvdims,  vdims,  nvatts) 
         ierr =nf_inq_varnatts (ncid, varid, nvatts) 
         year2daydta = year2day
         DO ji = 1, nvatts
            ierr = nf_inq_attname (ncid, varid, ji, attname)
            IF (ierr == nf_noerr) THEN
               latt = lenstr(attname)
               IF (attname(1:latt) == 'cycle_length') THEN
                  ierr = nf_get_att_FTYPE (ncid, varid, attname(1:latt), cycle_length)
                  IF (ierr == nf_noerr) THEN
                     year2daydta = cycle_length
                  ELSE
                     IF (lwp) write(numout,'(/1x,4A/)') 'SET_CYCLE ERROR while ', &
                     &        'reading attribute ''', attname(1:latt), '''.'
                  ENDIF
               ENDIF
            ENDIF
         END DO
      ENDIF

      IF ( ln_dust ) THEN
         lstr = lenstr(bioname)
         ierr = nf_open (bioname(1:lstr), nf_nowrite, ncid)
         IF ( ierr .NE. nf_noerr .AND. lwp ) THEN
            WRITE(numout,4) bioname
         ENDIF
         ierr = nf_inq_varid (ncid,"dust",varid)
         IF (ierr .NE. nf_noerr .AND. lwp ) THEN
            WRITE(numout,5) "dust", bioname
         ENDIF
         ierr = nf_inq_dimid(ncid,"dust_time",dimid)
         ierr = nf_inq_dimlen(ncid,dimid,nrec_dust)
         ALLOCATE( dustmp(GLOBAL_2D_ARRAY,nrec_dust), dustmo(GLOBAL_2D_ARRAY,12) )
         ALLOCATE( dust(PRIV_2D_BIOARRAY) )
         DO irec = 1, nrec_dust
            ierr = nf_fread(dustmp(START_2D_ARRAY,irec), ncid, varid, irec, r2dvar)
            IF (ierr .NE. nf_noerr .AND. lwp ) THEN
               WRITE(numout,6) "dust", irec
            ENDIF
         END DO
         ierr = nf_close(ncid)
         IF (lwp) WRITE(numout,*)
         IF (lwp) WRITE(numout,'(6x,A,1x,I4)') &
#ifdef MPI
         &                   'TRCINI_PISCES -- Read dust deposition ', mynode
#else
         &                   'TRCINI_PISCES -- Read dust deposition '
#endif

         DO irec = 1, nrec_dust
            DO jj = 1, LOCALMM
               DO ji = 1, LOCALLM
                  dustmo(ji,jj,irec) = dustmp(ji,jj,irec) / rmtss
               ENDDO
            ENDDO
         ENDDO
         !
         DEALLOCATE(dustmp)

      ENDIF
!
!    READ N DEPOSITION FROM ATMOSPHERE (use dust_time for time)
!    -------------------------------------
!
      IF (ln_ndepo) THEN
         lstr = lenstr(bioname)
         ierr = nf_open (bioname(1:lstr), nf_nowrite, ncid)
         IF (ierr .NE. nf_noerr .AND. lwp) THEN
            WRITE(numout,4) bioname
         ENDIF
         ierr = nf_inq_varid (ncid,"ndepo",varid)
         IF (ierr .NE. nf_noerr .AND. lwp ) THEN
            WRITE(numout,5) "ndepo", bioname
         ENDIF
         ierr = nf_inq_dimid(ncid,"dust_time",dimid)
         ierr = nf_inq_dimlen(ncid,dimid,nrec_dust)
         ALLOCATE( no3deptmp(GLOBAL_2D_ARRAY,nrec_dust), no3depmo(GLOBAL_2D_ARRAY,12) )
         ALLOCATE( no3dep(PRIV_2D_BIOARRAY) )
         DO irec = 1, nrec_dust
            ierr = nf_fread(no3deptmp(START_2D_ARRAY,irec), ncid, varid, irec, r2dvar)
            IF (ierr .NE. nf_noerr .AND. lwp ) THEN
               WRITE(numout,6) "ndepo", irec
            ENDIF
         END DO
         !
         ierr = nf_close(ncid)
         IF (lwp) WRITE(numout,*)
         IF (lwp) WRITE(numout,'(6x,A,1x,I4)') &
#ifdef MPI
         &                   'TRCINI_PISCES -- Read Nitrate deposition ', mynode
#else
         &                   'TRCINI_PISCES -- Read Nitrate deposition '
#endif

         DO irec = 1, nrec_dust
            DO jj = 1, LOCALMM
               DO ji = 1, LOCALLM
                  no3depmo(ji,jj,irec) = no3deptmp(ji,jj,irec)
               END DO
            END DO
         END DO
         !
         DEALLOCATE( no3deptmp )
         !
         lstr = lenstr(bioname)
         ierr = nf_open (bioname(1:lstr), nf_nowrite, ncid)
         IF (ierr .NE. nf_noerr .AND. lwp) THEN
            WRITE(numout,4) bioname
         ENDIF
         ierr = nf_inq_varid (ncid,"nhxdepo",varid)
         IF (ierr .NE. nf_noerr .AND. lwp ) THEN
            WRITE(numout,5) "nhxdepo", bioname
         ENDIF
         ierr = nf_inq_dimid(ncid,"dust_time",dimid)
         ierr = nf_inq_dimlen(ncid,dimid,nrec_dust)
         ALLOCATE( nh4deptmp(GLOBAL_2D_ARRAY,nrec_dust), nh4depmo(GLOBAL_2D_ARRAY,12) )
         ALLOCATE( nh4dep(PRIV_2D_BIOARRAY) )

         DO irec = 1, nrec_dust
            ierr = nf_fread(nh4deptmp(START_2D_ARRAY,irec), ncid, varid, irec, r2dvar)
            IF (ierr .NE. nf_noerr .AND. lwp ) THEN
               WRITE(numout,6) "nhxdepo", irec
            ENDIF
         END DO
         !
         ierr = nf_close(ncid)
         IF (lwp) WRITE(numout,*)
         IF (lwp) WRITE(numout,'(6x,A,1x,I4)') &
#ifdef MPI
         &                   'TRCINI_PISCES -- Read Ammoniun deposition ', mynode
#else
         &                   'TRCINI_PISCES -- Read Ammoniun deposition '
#endif

         DO irec = 1, nrec_dust
            DO jj= 1, LOCALMM
               DO ji =1, LOCALLM
                  nh4depmo(ji,jj,irec) = nh4deptmp(ji,jj,irec)
               END DO
            END DO
         END DO
         !
         DEALLOCATE( nh4deptmp )

      ENDIF

  4   FORMAT(/,' TRCINI_PISCES - unable to open forcing netCDF ',1x,A)
  5   FORMAT(/,' TRCINI_PISCES - unable to find forcing variable: ',A, &
     &                               /,14x,'in forcing netCDF  ',A)
  6   FORMAT(/,' TRCINI_PISCES - error while reading variable: ',A,2x, &
     &                                           ' at TIME index = ',i4)

 
      ! 
!      IF( ll_sbc ) CALL p4z_sbc( nit000 ) 
      !
      sedsilfrac = 0.03     ! percentage of silica loss in the sediments
      sedcalfrac = 0.6      ! percentage of calcite loss in the sediments
      !
   END SUBROUTINE p4z_sbc_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sbc                        ! Empty routine
   END SUBROUTINE p4z_sbc
#endif 


   !!======================================================================
END MODULE p4zsbc
