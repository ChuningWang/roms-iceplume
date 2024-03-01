#include "cppdefs.h"

MODULE sms_pisces   
   !!----------------------------------------------------------------------
   !!                     ***  sms_pisces.F90  ***  
   !! TOP :   PISCES Source Minus Sink variables
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2000-02 (O. Aumont) original code
   !!             3.2  !  2009-04 (C. Ethe & NEMO team) style
   !!----------------------------------------------------------------------

#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                         PISCES model
   !!----------------------------------------------------------------------
   USE par_pisces
   USE ocean2pisces
   USE trc

   IMPLICIT NONE
   PUBLIC

#include "ocean2pisces.h90"

   !!*  Time variables
   INTEGER  ::   numnatp_ref = -1           !! Logical units for namelist pisces
   INTEGER  ::   numnatp_cfg = -1           !! Logical units for namelist pisces
   INTEGER  ::   numonp      = -1           !! Logical unit for namelist pisces output

   !!* Model used
   LOGICAL  ::  ln_p2z            !: Flag to use LOBSTER model
   LOGICAL  ::  ln_p4z            !: Flag to use PISCES  model
   LOGICAL  ::  ln_p5z            !: Flag to use PISCES  quota model
   LOGICAL  ::  ln_ligand         !: Flag to enable organic ligands
   LOGICAL  ::  ln_sediment       !: Flag to enable sediment module

   !!*  Time variables
   INTEGER  ::   nrdttrc           !: ???
   INTEGER  ::   niter1max, niter2max           !: ???
   INTEGER  ::   neos
   REAL(wp) ::   rfact , rfactr    !: ???
   REAL(wp) ::   rfact2, rfact2r   !: ???
   REAL(wp) ::   xstep             !: Time step duration for biology

   !!*  Biological parameters 
   REAL(wp) ::   rno3              !: ???
   REAL(wp) ::   o2ut              !: ???
   REAL(wp) ::   po4r              !: ???
   REAL(wp) ::   rdenit            !: ???
   REAL(wp) ::   rdenita           !: ???
   REAL(wp) ::   o2nit             !: ???
   REAL(wp) ::   wsbio, wsbio2     !: ???
   REAL(wp) ::   wsbio2max         !: ???
   REAL(wp) ::   wsbio2scale       !: ???
   LOGICAL  ::   ln_sink_new       !: 
   REAL(wp) ::   xkmort            !: ???
   REAL(wp) ::   ferat3            !: ???
   REAL(wp) ::   ldocp             !: ???
   REAL(wp) ::   ldocz             !: ???
   REAL(wp) ::   lthet             !: ???
   REAL(wp) ::   no3rat3           !: ???
   REAL(wp) ::   po4rat3           !: ???

   !!* Damping 
   LOGICAL  ::   ln_pisdmp         !: relaxation or not of nutrients to a mean value
                                   !: when initialize from a restart file 
   LOGICAL  ::   ln_pisclo         !: Restoring or not of nutrients to initial value
                                   !: on close seas

   REAL(wp), DIMENSION(:), ALLOCATABLE ::   tra_ctl         !: previous trend values

   !!*  Biological fluxes for light
   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  neln  !: number of T-levels+ 1 in the euphotic layer
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  heup  !: euphotic layer depth
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  fr_i  !: euphotic layer depth
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  tmask
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  etot3
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  etot, etot_ndcy      !: PAR over 24h in case of diurnal cycle
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  enano, ediat   !: PAR for phyto, nano and diat
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  enanom, ediatm !: PAR for phyto, nano and diat 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  epico          !: PAR for pico
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  epicom         !: PAR for pico
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  emoy           !: averaged PAR in the mixed layer
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  heup_01 !: Absolute euphotic layer depth
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  xksi  !:  LOBSTER : zooplakton closure

   !!*  Biological fluxes for primary production
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:)  ::   xksimax    !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   biron      !: bioavailable fraction of iron
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   plig       !: proportion of iron organically complexed

   !!*  Sinking speed
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio3   !: POC sinking speed
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio4   !: GOC sinking speed

   !!*  SMS for the organic matter
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xfracal    !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   nitrfac    !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   nitrfac2   !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   orem       !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xdiss      !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prodcal    !: Calcite production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prodpoc    !: Calcite production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   conspoc    !: Calcite production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prodgoc    !: Calcite production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   consgoc    !: Calcite production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   blim       !: bacterial production factor
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sizen      !: size of diatoms
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sizep      !: size of diatoms
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sized      !: size of diatoms

   !!* Variable for chemistry of the CO2 cycle
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak13    !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak23    !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aksp    !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hi      !: 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   excess   !: 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aphscale   !:

   !!* Temperature dependancy of SMS terms
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfunc    !: Temp.  dependancy of various biological rates
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfunc2   !: Temp.  dependancy of mesozooplankton rates

   !! * Shared module variables

  REAL(wp), PARAMETER     ::  rtrn = 1.e-20
  INTEGER :: jip1 = 120
  INTEGER :: jjp1 = 60
  INTEGER :: jip2 = 38
  INTEGER :: jjp2 = 20
  INTEGER :: jkp = KSURF

   !!* Substitution
#  include "ocean2pisces.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.2 , LOCEAN-IPSL (2009) 
   !! $Id: sms_pisces.F90 1830 2010-04-12 13:03:51Z cetlod $ 
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!======================================================================   

CONTAINS

  INTEGER FUNCTION sms_pisces_alloc()
      !!----------------------------------------------------------------------
      !!        *** ROUTINE sms_pisces_alloc ***
      !!----------------------------------------------------------------------
      INTEGER ::   ierr(11)        ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !*  Biological fluxes for light
      ALLOCATE( etot   (PRIV_3D_BIOARRAY), neln(PRIV_2D_BIOARRAY), heup(PRIV_2D_BIOARRAY),    &
        &       heup_01(PRIV_2D_BIOARRAY), xksi(PRIV_2D_BIOARRAY), STAT=ierr(1) )
      !
      ALLOCATE( fr_i  (PRIV_2D_BIOARRAY),                                &
                tmask (PRIV_3D_BIOARRAY), etot3 (PRIV_3D_BIOARRAY) ,     &
                enano (PRIV_3D_BIOARRAY), ediat (PRIV_3D_BIOARRAY) ,     &
                enanom(PRIV_3D_BIOARRAY), ediatm(PRIV_3D_BIOARRAY) ,     &
                etot_ndcy(PRIV_3D_BIOARRAY), emoy(PRIV_3D_BIOARRAY),  STAT=ierr(2) )
      !
      !*  Biological fluxes for primary production
      ALLOCATE( xksimax(PRIV_2D_BIOARRAY), biron(PRIV_3D_BIOARRAY),  STAT=ierr(3) )
         !
      !*  SMS for the organic matter
      ALLOCATE( xfracal (PRIV_3D_BIOARRAY), orem    (PRIV_3D_BIOARRAY),   &
         &      nitrfac(PRIV_3D_BIOARRAY) , nitrfac2(PRIV_3D_BIOARRAY),   & 
         &      prodcal(PRIV_3D_BIOARRAY) , xdiss   (PRIV_3D_BIOARRAY),   &
         &      prodpoc(PRIV_3D_BIOARRAY) , conspoc(PRIV_3D_BIOARRAY) ,   &
         &      prodgoc(PRIV_3D_BIOARRAY) , consgoc(PRIV_3D_BIOARRAY) ,   &
         &      blim   (PRIV_3D_BIOARRAY) ,     STAT=ierr(4) )
         !
      !* Variable for chemistry of the CO2 cycle
      ALLOCATE( ak13(PRIV_3D_BIOARRAY)    ,                                  &
         &      ak23(PRIV_3D_BIOARRAY)    , aksp  (PRIV_3D_BIOARRAY) ,       &
         &      aphscale(PRIV_3D_BIOARRAY), excess(PRIV_3D_BIOARRAY) ,       &
         &      hi  (PRIV_3D_BIOARRAY)    ,    STAT=ierr(5) )
         !
      !* Temperature dependancy of SMS terms
      ALLOCATE( tgfunc(PRIV_3D_BIOARRAY)  , tgfunc2(PRIV_3D_BIOARRAY),   STAT=ierr(6) )
      !
      !* Sinking speed
      ALLOCATE( wsbio3 (PRIV_3D_BIOARRAY) , wsbio4 (PRIV_3D_BIOARRAY),     &
      &                             STAT=ierr(7) )
      !
      IF( ln_ligand ) THEN
         ALLOCATE( plig(PRIV_3D_BIOARRAY) ,                          STAT=ierr(9) )
      ENDIF
      IF( ln_p5z ) THEN
         !
         ALLOCATE( epico(PRIV_3D_BIOARRAY), epicom(PRIV_3D_BIOARRAY) ,   STAT=ierr(10) )

         !*  Size of phytoplankton cells
         ALLOCATE( sizen(PRIV_3D_BIOARRAY), sizep(PRIV_3D_BIOARRAY),         &
           &       sized(PRIV_3D_BIOARRAY), STAT=ierr(11) )
      ENDIF

      !
      sms_pisces_alloc = MAXVAL( ierr )
      !
      IF( sms_pisces_alloc /= 0 )   CALL ctl_warn('sms_pisces_alloc: failed to allocate arrays')
      !
   END FUNCTION sms_pisces_alloc


   SUBROUTINE tracer_stat( kt )
      !!----------------------------------------------------------------------
      !!                    ***  trc_rst_stat  ***
      !!
      !! ** purpose  :   Compute tracers statistics
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)  :: kt
      INTEGER  :: ji, jj, jk, jn
      REAL(wp) :: ztra, zmin, zmax, zmean, areatot, zcoef
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY,jptra)  :: ptra
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY)        :: zmask, zvol
      !!----------------------------------------------------------------------

      IF( lwp ) THEN
         WRITE(numout,*) 
         WRITE(numout,*) ' TRACER STAT at time-step kt = ', kt
         WRITE(numout,*) 
      ENDIF
      !
! to have coherent units when calling tracer_stat
      IF( kt .eq. nit000 ) THEN
        zcoef = 1.e-6
      ELSE
        zcoef = 1.
      ENDIF

      DO jn = 1, jptra
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  ptra(ji,jj,jk,jn) = trb(ji,jj,K,jn) * zcoef
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      areatot = 0.                                                           ! total volume of the ocean 
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE          ! masked grid volume
               zvol(ji,jj,jk)  = cvol(ji,jj,K)
               zmask(ji,jj,jk) = tmask(ji,jj,jk) * tmask_i(ji,jj) 
               areatot         = areatot + zvol(ji,jj,jk)
            ENDDO
        ENDDO
     ENDDO
     IF( lk_mpp )   CALL mpp_sum( areatot )     ! sum over the global domain  

     DO jn = 1, jptra
         ztra = 0.
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE          ! masked grid volume
                  ztra  = ztra + ptra(ji,jj,jk,jn) * zvol(ji,jj,jk) 
               ENDDO
            ENDDO
         ENDDO
         zmin  = MINVAL( ptra(:,:,:,jn), mask= ( zmask(:,:,:) /= 0. ) ) 
         zmax  = MAXVAL( ptra(:,:,:,jn), mask= ( zmask(:,:,:) /= 0. ) ) 
         IF( lk_mpp ) THEN
            CALL mpp_sum( ztra )      ! min over the global domain
            CALL mpp_min( zmin )      ! min over the global domain
            CALL mpp_max( zmax )      ! max over the global domain
         END IF
         zmean  = ztra / areatot
         IF(lwp) WRITE(numout,9000) jn, TRIM( ctrcnm(jn) ), zmean, zmin, zmax
      END DO
      WRITE(numout,*) 
9000  FORMAT(' tracer nb :',i2,'    name :',a10,'    mean :',e18.10,'    min :',e18.10, '    max :',e18.10 )
      !
   END SUBROUTINE tracer_stat

   SUBROUTINE prt_ctl_trc( charout, ltra )

      CHARACTER (len=*),  INTENT(in)           :: charout   ! information about the tab3d array
      CHARACTER (len=*),  INTENT(in), OPTIONAL :: ltra    ! information about the tab3d array
      INTEGER :: ji, jj, jk, jn
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY,jptra)        :: ztab
      REAL(wp)  :: zsum,  zvctl  


      IF( PRESENT( ltra ) ) THEN
         DO jn = 1, jptra
            DO jk = KRANGE
               DO jj = JRANGE
                  DO ji = IRANGE
                     ztab(ji,jj,jk,jn) = tra(ji,jj,jk,jn) * tmask(ji,jj,jk)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ELSE 
         DO jn = 1, jptra
            DO jk = KRANGE
               DO jj = JRANGE
                  DO ji = IRANGE
                     ztab(ji,jj,jk,jn) = trn(ji,jj,K,jn) * tmask(ji,jj,jk)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      WRITE(numout,*) charout

      IF( PRESENT( ltra ) ) THEN
          DO jn = 1, jptra
             zvctl  = tra_ctl(jn)
             zsum   = SUM( ztab(:,:,:,jn) )
             IF( lk_mpp ) CALL mpp_sum( zsum )      ! min over the global domain
             IF( lwp ) WRITE(numout,FMT="(3x,a10,' : ',D23.16)") TRIM(ctrcnm(jn)), zsum-zvctl
             tra_ctl(jn) = zsum
          END DO
       ELSE
          DO jn = 1, jptra
             zvctl  = tra_ctl(jn)
             zsum   = SUM( ztab(:,:,:,jn) )
             IF( lk_mpp ) CALL mpp_sum( zsum )      ! min over the global domain
             IF( lwp ) WRITE(numout,FMT="(3x,a10,' : ',D23.16)") TRIM(ctrcnm(jn)), zsum
          END DO
      ENDIF

   END SUBROUTINE prt_ctl_trc      

   SUBROUTINE prt_ctl_trc_ini

      ALLOCATE( tra_ctl(jptra) )
      tra_ctl(:) = 0.e0           ! Initialization to zero

   END SUBROUTINE prt_ctl_trc_ini

  SUBROUTINE prt_ctl_trc_info( clinfo )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE prt_ctl_trc_info  ***
      !!
      !! ** Purpose : - print information without any computation
      !!----------------------------------------------------------------------
      CHARACTER (len=*), INTENT(in) ::   clinfo      ! information to print
      !! 
      !
   END SUBROUTINE prt_ctl_trc_info



#else
   !!----------------------------------------------------------------------   
   !!  Empty module :                                     NO PISCES model
   !!----------------------------------------------------------------------
#endif
   
END MODULE sms_pisces    
