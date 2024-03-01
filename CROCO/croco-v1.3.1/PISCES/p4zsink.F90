#include "cppdefs.h"

MODULE p4zsink
   !!======================================================================
   !!                         ***  MODULE p4zsink  ***
   !! TOP :  PISCES  vertical flux of particulate matter due to gravitational sinking
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!             3.5  !  2012-07  (O. Aumont) Introduce potential time-splitting
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!   p4z_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   p4z_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE sms_pisces      !  PISCES Source Minus Sink variables
!   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sink         ! called in p4zbio.F90
   PUBLIC   p4z_sink_alloc

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"


   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinking, sinking2  !: POC sinking fluxes 
   !                                                          !  (different meanings depending on the parameterization)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkingn, sinking2n  !: POC sinking fluxes 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkingp, sinking2p  !: POC sinking fluxes 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkcal, sinksil   !: CaCO3 and BSi sinking fluxes
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkfer            !: Small BFe sinking fluxes
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkfer2           !: Big iron sinking fluxes

   INTEGER  :: ik100

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zsink.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   !!----------------------------------------------------------------------
   !!   'standard sinking parameterisation'                  ???
   !!----------------------------------------------------------------------

   SUBROUTINE p4z_sink ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      INTEGER  ::   ji, jj, jk, jit
      INTEGER  ::   iiter1, iiter2
      INTEGER  ::   ik1
      CHARACTER (len=25) :: charout
      REAL(wp) :: zmsk, zmax, zwsmax, zfact
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zw2d
      !!---------------------------------------------------------------------

      ! Initialization of some global variables
      ! ---------------------------------------
      prodpoc(:,:,:) = 0.
      conspoc(:,:,:) = 0.
      prodgoc(:,:,:) = 0.
      consgoc(:,:,:) = 0.

      !
      !    Sinking speeds of detritus is increased with depth as shown
      !    by data and from the coagulation theory
      !    -----------------------------------------------------------
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zmax  = MAX( heup_01(ji,jj), hmld(ji,jj) )
               zfact = MAX( 0., gdepw_n(ji,jj,jk+1) - zmax ) / wsbio2scale
               wsbio4(ji,jj,jk) = wsbio2 + MAX(0., ( wsbio2max - wsbio2 )) * zfact
            END DO
         END DO
      END DO

      ! limit the values of the sinking speeds to avoid numerical instabilities  
      wsbio3(:,:,:) = wsbio

      !
      !  Initializa to zero all the sinking arrays 
      !   -----------------------------------------
      sinking (:,:,:) = 0.e0
      sinking2(:,:,:) = 0.e0
      sinkcal (:,:,:) = 0.e0
      sinkfer (:,:,:) = 0.e0
      sinksil (:,:,:) = 0.e0
      sinkfer2(:,:,:) = 0.e0

      IF( .NOT. ln_sink_new ) THEN

!      LIMIT THE VALUES OF THE SINKING SPEEDS 
!      TO AVOID NUMERICAL INSTABILITIES

         !
         ! OA This is (I hope) a temporary solution for the problem that may 
         ! OA arise in specific situation where the CFL criterion is broken 
         ! OA for vertical sedimentation of particles. To avoid this, a time
         ! OA splitting algorithm has been coded. A specific maximum
         ! OA iteration number is provided and may be specified in the namelist 
         ! OA This is to avoid very large iteration number when explicit free
         ! OA surface is used (for instance). When niter?max is set to 1, 
         ! OA this computation is skipped. The crude old threshold method is 
         ! OA then applied. This also happens when niter exceeds nitermax.

         IF( MAX( niter1max, niter2max ) == 1 ) THEN
           iiter1 = 1
           iiter2 = 1
         ELSE
            iiter1 = 1
            iiter2 = 1
            DO jk = KRANGE
               DO jj = JRANGE
                  DO ji = IRANGE
                     IF( tmask(ji,jj,jk) == 1) THEN
                        zwsmax =  0.5 * e3t_n(ji,jj,K) / xstep
                        iiter1 =  MAX( iiter1, INT( wsbio3(ji,jj,jk) / zwsmax ) )
                        iiter2 =  MAX( iiter2, INT( wsbio4(ji,jj,jk) / zwsmax ) )
                     ENDIF
                  END DO
               END DO
            END DO
            IF( lk_mpp ) THEN
               CALL mpp_max( iiter1 )
               CALL mpp_max( iiter2 )
            ENDIF
            iiter1 = MIN( iiter1, niter1max )
            iiter2 = MIN( iiter2, niter2max )
         ENDIF

         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  IF( tmask(ji,jj,jk) == 1 ) THEN
                     zwsmax = 0.5 * e3t_n(ji,jj,K) / xstep
                     wsbio3(ji,jj,jk) = MIN( wsbio3(ji,jj,jk), zwsmax * FLOAT( iiter1 ) )
                     wsbio4(ji,jj,jk) = MIN( wsbio4(ji,jj,jk), zwsmax * FLOAT( iiter2 ) )
                  ENDIF
               END DO
            END DO
         END DO
!        Compute the sedimentation term using p4zsink2 for all
!        the sinking particles
!        -----------------------------------------------------
         DO jit = 1, iiter1
            CALL p4z_sink2_std( wsbio3, sinking , jppoc, iiter1 )
            CALL p4z_sink2_std( wsbio3, sinkfer , jpsfe, iiter1 )
         END DO

         DO jit = 1, iiter2
            CALL p4z_sink2_std( wsbio4, sinking2, jpgoc, iiter2 )
            CALL p4z_sink2_std( wsbio4, sinkfer2, jpbfe, iiter2 )
            CALL p4z_sink2_std( wsbio4, sinksil , jpgsi, iiter2 )
            CALL p4z_sink2_std( wsbio4, sinkcal , jpcal, iiter2 )
         END DO

         IF( ln_p5z ) THEN
            sinkingn (:,:,:) = 0.e0
            sinking2n(:,:,:) = 0.e0
            sinkingp (:,:,:) = 0.e0
            sinking2p(:,:,:) = 0.e0

            !   Compute the sedimentation term using p4zsink2 for all the sinking particles
            !   -----------------------------------------------------
            DO jit = 1, iiter1
               CALL p4z_sink2_std( wsbio3, sinkingn , jppon, iiter1 )
               CALL p4z_sink2_std( wsbio3, sinkingp , jppop, iiter1 )
            END DO


            DO jit = 1, iiter2
               CALL p4z_sink2_std( wsbio4, sinking2n, jpgon, iiter2 )
               CALL p4z_sink2_std( wsbio4, sinking2p, jpgop, iiter2 )
            END DO

         ENDIF

      ELSE

         CALL p4z_sink2_new( wsbio3, sinking , jppoc )
         CALL p4z_sink2_new( wsbio3, sinkfer , jpsfe )
         !
         CALL p4z_sink2_new( wsbio4, sinking2, jpgoc )
         CALL p4z_sink2_new( wsbio4, sinkfer2, jpbfe )
         CALL p4z_sink2_new( wsbio4, sinksil , jpdsi )
         CALL p4z_sink2_new( wsbio4 , sinkcal , jpcal )

         IF( ln_p5z ) THEN
            sinkingn (:,:,:) = 0.e0
            sinking2n(:,:,:) = 0.e0
            sinkingp (:,:,:) = 0.e0
            sinking2p(:,:,:) = 0.e0

            !   Compute the sedimentation term using p4zsink2 for all the sinking particles
            !   -----------------------------------------------------
            CALL p4z_sink2_new( wsbio3, sinkingn , jppon )
            CALL p4z_sink2_new( wsbio3, sinkingp , jppop )
            CALL p4z_sink2_new( wsbio4, sinking2n, jpgon )
            CALL p4z_sink2_new( wsbio4, sinking2p, jpgop )
         ENDIF

     ENDIF

     ik100 =10

#if defined key_iomput
     IF( lk_iomput ) THEN
       IF( knt == nrdttrc ) THEN
          ALLOCATE( zw2d(PRIV_2D_BIOARRAY), zw3d(PRIV_3D_BIOARRAY) )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "EPC100" ) )  THEN
              zw2d(:,:) = ( sinking(:,:,ik100) + sinking2(:,:,ik100) ) * zfact * tmask(:,:,1) ! Export of carbon at 100m
              CALL iom_put( "EPC100"  , zw2d )
          ENDIF
          IF( iom_use( "EPFE100" ) )  THEN
              zw2d(:,:) = ( sinkfer(:,:,ik100) + sinkfer2(:,:,ik100) ) * zfact * tmask(:,:,1) ! Export of iron at 100m
              CALL iom_put( "EPFE100"  , zw2d )
          ENDIF
          IF( iom_use( "EPCAL100" ) )  THEN
              zw2d(:,:) = sinkcal(:,:,ik100) * zfact * tmask(:,:,1) ! Export of calcite at 100m
              CALL iom_put( "EPCAL100"  , zw2d )
          ENDIF
          IF( iom_use( "EPSI100" ) )  THEN
              zw2d(:,:) =  sinksil(:,:,ik100) * zfact * tmask(:,:,1) ! Export of bigenic silica at 100m
              CALL iom_put( "EPSI100"  , zw2d )
          ENDIF
          IF( iom_use( "EXPC" ) )  THEN
              zw3d(:,:,:) = ( sinking(:,:,:) + sinking2(:,:,:) ) * zfact * tmask(:,:,:) ! Export of carbon in the water column
              CALL iom_put( "EXPC"  , zw3d )
          ENDIF
          IF( iom_use( "EXPFE" ) )  THEN
              zw3d(:,:,:) = ( sinkfer(:,:,:) + sinkfer2(:,:,:) ) * zfact * tmask(:,:,:) ! Export of iron 
              CALL iom_put( "EXPFE"  , zw3d )
          ENDIF
          IF( iom_use( "EXPCAL" ) )  THEN
              zw3d(:,:,:) = sinkcal(:,:,:) * zfact * tmask(:,:,:) ! Export of calcite 
              CALL iom_put( "EXPCAL"  , zw3d )
          ENDIF
          IF( iom_use( "EXPSI" ) )  THEN
              zw3d(:,:,:) = sinksil(:,:,:) * zfact * tmask(:,:,:) ! Export of bigenic silica
              CALL iom_put( "EXPSI"  , zw3d )
          ENDIF
          ! 
          DEALLOCATE( zw2d, zw3d )
        ENDIF
      ENDIF
#endif
      !
#if defined key_trc_diaadd
      zfact = 1.e3 * rfact2r
      ik1  = ik100 + 1
      DO jj = JRANGE
         DO ji = IRANGE
            zmsk = zfact * tmask(ji,jj,1)
            trc2d(ji,jj,jp_sinkco2) = ( sinking(ji,jj,ik1) + sinking2(ji,jj,ik1) ) * zmsk ! export of carbon at 100m
            trc2d(ji,jj,jp_sinkfer) = ( sinkfer(ji,jj,ik1) + sinkfer2(ji,jj,ik1) ) * zmsk ! export of biogenic iron
            trc2d(ji,jj,jp_sinkcal) =   sinkcal(ji,jj,ik1)  * zmsk   ! export of calcite
            trc2d(ji,jj,jp_sinksil) =   sinksil(ji,jj,ik1)  * zmsk   ! export of biogenic silica
         END DO
      END DO
#endif

      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( charout, ltra='tra')
!         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
   END SUBROUTINE p4z_sink

   SUBROUTINE p4z_sink2_std( pwsink, psinkflx, jp_tra, kiter )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink2  ***
      !!
      !! ** Purpose :   Compute the sedimentation terms for the various sinking
      !!     particles. The scheme used to compute the trends is based
      !!     on MUSCL.
      !!
      !! ** Method  : - this ROUTINE compute not exactly the advection but the
      !!      transport term, i.e.  div(u*tra).
      !!---------------------------------------------------------------------
#ifdef AGRIF
      USE ocean2pisces
#endif
      INTEGER , INTENT(in   )                         ::   jp_tra    ! tracer index index
      INTEGER , INTENT(in   )                         ::   kiter     ! number of iterations for time-splitting
      REAL(wp), INTENT(in   ), DIMENSION(PRIV_2D_BIOARRAY, jpk) ::   pwsink    ! sinking speed
      REAL(wp), INTENT(inout), DIMENSION(PRIV_2D_BIOARRAY, jpk+1) ::   psinkflx  ! sinking fluxe
      !!
      INTEGER  ::   ji, jj, jk, jn, jnt
      REAL(wp) ::   zigma,zew,zign, zflx, zstep
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY, jpk+1) ::  ztraz, zakz
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY, jpk+1) ::  zwsink2
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) ::  ztrn, ztmp, zdept, zmask
      !!---------------------------------------------------------------------

      zstep = rfact2 / FLOAT( kiter ) / 2.

      DO jk = 1, jpk
         DO jj = JRANGE
            DO ji = IRANGE
               zmask(ji,jj,jk) = tmask(ji,jj,jk)
               zdept(ji,jj,jk) = e3t_n(ji,jj,K)
               ztrn (ji,jj,jk) = trb(ji,jj,K,jp_tra)
           END DO
         END DO
      ENDDO

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               ztraz(ji,jj,jk) = 0.e0
               zakz (ji,jj,jk) = 0.e0
               ztmp(ji,jj,jk)  = ztrn(ji,jj,jk)
           END DO
         END DO
      END DO

      DO jk = 1, jpk-1
         DO jj = JRANGE
            DO ji = IRANGE
               zwsink2(ji,jj,jk+1) = -pwsink(ji,jj,jk) / rday * zmask(ji,jj,jk+1)
           END DO
         END DO
      END DO
      !
      DO jj = JRANGE
         DO ji = IRANGE
            zwsink2(ji,jj,1)   = 0.
            zwsink2(ji,jj,jpk+1) = 0.
         END DO
      END DO
      !
      ! Vertical advective flux
      DO jnt = 1, 2
         !  first guess of the slopes interior values
         DO jk = 2, jpk
            DO jj = JRANGE
               DO ji = IRANGE
                  ztraz(ji,jj,jk) = ( ztrn(ji,jj,jk-1) - ztrn(ji,jj,jk) ) * zmask(ji,jj,jk)
              END DO
            END DO
         END DO
         DO jj = JRANGE
            DO ji = IRANGE
               ztraz(ji,jj,1  ) = 0.0
               ztraz(ji,jj,jpk+1) = 0.0
            END DO
         END DO


         ! slopes
!         DO jk = KRANGEL
         DO jk = 2, jpk
            DO jj = JRANGE
               DO ji = IRANGE
                  zign = 0.25 + SIGN( 0.25, ztraz(ji,jj,jk) * ztraz(ji,jj,jk+1) )
                  zakz(ji,jj,jk) = ( ztraz(ji,jj,jk) + ztraz(ji,jj,jk+1) ) * zign
               END DO
            END DO
         END DO

         ! Slopes limitation
         DO jk = 2, jpk
            DO jj = JRANGE
               DO ji = IRANGE
                  zakz(ji,jj,jk) = SIGN( 1., zakz(ji,jj,jk) ) *        &
                     &             MIN( ABS( zakz(ji,jj,jk) ), 2. * ABS(ztraz(ji,jj,jk+1)), 2. * ABS(ztraz(ji,jj,jk) ) )
               END DO
            END DO
         END DO


         ! vertical advective flux
         DO jk = 1, jpk
            DO jj = JRANGE
               DO ji = IRANGE
                  zigma = zwsink2(ji,jj,jk+1) * zstep / e3w_n(ji,jj,jk+1)
                  zew   = zwsink2(ji,jj,jk+1)
                  psinkflx(ji,jj,jk+1) = -zew * ( ztrn(ji,jj,jk) - 0.5 * ( 1 + zigma ) * zakz(ji,jj,jk) ) * zstep
               END DO
            END DO
         END DO
         !
         ! Boundary conditions
         DO jj = JRANGE
            DO ji = IRANGE
               psinkflx(ji,jj,1  ) = 0.e0
               psinkflx(ji,jj,jpk+1) = 0.e0
            END DO
         END DO


         DO jk = 1, jpk
            DO jj = JRANGE
               DO ji = IRANGE
                  zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / zdept(ji,jj,jk)
                  ztrn(ji,jj,jk) = ztrn(ji,jj,jk) + zflx
               END DO
            END DO
         END DO


      ENDDO

      DO jk = 1, jpk
         DO jj = JRANGE
            DO ji = IRANGE
               zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / zdept(ji,jj,jk)
               ztmp(ji,jj,jk) = ztmp(ji,jj,jk) + 2. * zflx
            END DO
         END DO
      END DO

      DO jk = 1, jpk
         DO jj = JRANGE
            DO ji = IRANGE
               trn(ji,jj,K,jp_tra) = ztmp(ji,jj,jk)
               psinkflx(ji,jj,jk)   = 2. * psinkflx(ji,jj,jk)
            END DO
         END DO
      END DO

      !
   END SUBROUTINE p4z_sink2_std

   SUBROUTINE p4z_sink2_new( pwsink, psinkflx, jp_tra )
     !======================================================================
     !                                                                      !
     !  This routine computes the vertical settling (sinking) of suspended  !
     !  sediment via a semi-Lagrangian advective flux algorithm. It uses a  !
     !  parabolic,  vertical reconstructuion of the suspended particle in  !
     !  the water column with PPT/WENO constraints to avoid oscillations.   !
     !                                                                      !
     !  References:                                                         !
     !                                                                      !
     !  Colella, P. and P. Woodward, 1984: The piecewise parabolic method   !
     !    (PPM) for gas-dynamical simulations, J. Comp. Phys., 54, 174-201. !
     !                                                                      !
     !  Liu, X.D., S. Osher, and T. Chan, 1994: Weighted essentially        !
     !    nonoscillatory shemes, J. Comp. Phys., 115, 200-212.              !
     !                                                                      !
     !  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
     !    Arango, 2008:  Development of a three-dimensional,  regional,     !
     !    coupled wave, current, and sediment-transport model, Computers    !
     !    & Geosciences, 34, 1284-1306.                                     !
     !                                                                      !
     !=======================================================================
#ifdef AGRIF
      USE ocean2pisces
#endif
      INTEGER , INTENT(in   )                         ::   jp_tra    ! tracer index index
      REAL(wp), INTENT(in   ), DIMENSION(PRIV_2D_BIOARRAY, jpk) ::   pwsink    ! sinking speed
      REAL(wp), INTENT(inout), DIMENSION(PRIV_2D_BIOARRAY, jpk+1) ::   psinkflx  ! sinking fluxe
!
!  Local variable declarations.
!
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) ::  ztrn,  zdept, zmask

      INTEGER :: indx, ised, ks,ji,jj,jk
      REAL(wp) :: cff, cu, cffL, cffR, dltL, dltR
      INTEGER, DIMENSION(jpi,jpk) :: ksource
      REAL(wp), DIMENSION(jpi,jpk+1) :: FC
      REAL(wp), dimension(jpi,jpk) :: Hz_inv
      REAL(wp), dimension(jpi,jpk) :: Hz_inv2
      REAL(wp), dimension(jpi,jpk) :: Hz_inv3
      REAL(wp), dimension(jpi,jpk) :: qc
      REAL(wp), dimension(jpi,jpk) :: qR
      REAL(wp), dimension(jpi,jpk) :: qL
      REAL(wp), dimension(jpi,jpk) :: WR
      REAL(wp), dimension(jpi,jpk) :: WL


     ztrn(:,:,:) = 0.
     zdept(:,:,:) = 0.
     zmask(:,:,:) = 0.

     ksource(:,:) = 0.
     FC(:,:) = 0.
     Hz_inv(:,:) = 0.
     Hz_inv2(:,:) = 0.
     Hz_inv3(:,:) = 0.
     qc(:,:) = 0.
     qR(:,:) = 0.
     qL(:,:) = 0.
     WR(:,:) = 0.
     WL(:,:) = 0.


      DO jk = 1, jpk
         DO jj = JRANGE
            DO ji = IRANGE
               zmask(ji,jj,jk) = tmask(ji,jj,jk)
               zdept(ji,jj,jk) = e3t_n(ji,jj,K)
               ztrn (ji,jj,jk) = trn(ji,jj,K,jp_tra)
           END DO
         END DO
      ENDDO

!
!-----------------------------------------------------------------------
!  Add sediment vertical sinking (settling) term.
!-----------------------------------------------------------------------
!
!  Compute inverse thicknesses to avoid repeated divisions.
!
      DO jj = JRANGE

         DO jk= 1, jpk
            DO ji= IRANGE
               Hz_inv(ji,jk)=1./zdept(ji,jj,jk)
            ENDDO
          ENDDO
          DO jk = 2, jpk
             DO ji = IRANGE
                Hz_inv2(ji,jk)=1./(zdept(ji,jj,jk)+zdept(ji,jj,jk-1))
             ENDDO
           ENDDO
           DO jk = 2, jpk-1
              DO ji = IRANGE
                 Hz_inv3(ji,jk)=1./(zdept(ji,jj,jk+1)+zdept(ji,jj,jk)+zdept(ji,jj,jk-1))
              ENDDO
           ENDDO
!
!  Copy concentration of suspended sediment into scratch array "qc"
!  (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for sediment
!  concentration.
!
           DO ised=1,1
              DO jk = 1, jpk
                 DO ji = IRANGE
                    qc(ji,jk)=ztrn(ji,jj,jk)
                 ENDDO
              ENDDO
!
!
!-----------------------------------------------------------------------
!  Vertical sinking of suspended sediment.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of suspended sediment "qc" in terms
!  of a set of parabolic segments within each grid box. Then, compute
!  semi-Lagrangian flux due to sinking.
!
          DO jk = 2, jpk
            DO ji = IRANGE
              FC(ji,jk)=(qc(ji,jk-1)-qc(ji,jk))*Hz_inv2(ji,jk)
            END DO
          END DO
          DO jk = 2, jpk-1
            DO ji = IRANGE
              dltR=zdept(ji,jj,jk)*FC(ji,jk)
              dltL=zdept(ji,jj,jk)*FC(ji,jk+1)
              cff=zdept(ji,jj,jk+1)+2.*zdept(ji,jj,jk)+zdept(ji,jj,jk-1)
              cffR=cff*FC(ji,jk)
              cffL=cff*FC(ji,jk+1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
              IF ((dltR*dltL).LE.0.) THEN
                dltR=0.
                dltL=0.
              ELSE IF (ABS(dltR).GE.(cffL)) THEN
                dltR=cffL
              ELSE IF (ABS(dltL).GT.ABS(cffR)) THEN
                dltL=cffR
              ENDIF
!
!  Compute right and left side values (qR,qL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because qL(k+1)-qR(k) may still have different sign than
!        qc(k+1)-qc(k).  This possibility is excluded, after qL and qR
!        are reconciled using WENO procedure.
!
              cff=(dltR-dltL)*Hz_inv3(ji,jk)
              dltR=dltR-cff*zdept(ji,jj,jk-1)
              dltL=dltL+cff*zdept(ji,jj,jk+1)
              qR(ji,jk)=qc(ji,jk)+dltR
              qL(ji,jk)=qc(ji,jk)-dltL
              WR(ji,jk)=(2.*dltR-dltL)**2
              WL(ji,jk)=(dltR-2.*dltL)**2
            END DO
          END DO
          cff=1.e-14
          DO jk = 2, jpk-2
            DO ji = IRANGE
              dltL=MAX(cff,WL(ji,jk  ))
              dltR=MAX(cff,WR(ji,jk-1))
              qR(ji,jk)=(dltR*qR(ji,jk)+dltL*qL(ji,jk-1))/(dltR+dltL)
              qL(ji,jk-1)=qR(ji,jk)
            ENDDO
          ENDDO
          DO ji= IRANGE
            FC(ji,1)=0.              ! no-flux boundary condition
!
            qR(ji,1)=qc(ji,1)         ! default strictly monotonic
            qL(ji,1)=qc(ji,1)         ! conditions
            qR(ji,2)=qc(ji,1)
!
            qL(ji,jpk-1)=qc(ji,jpk)                 ! bottom grid boxes are
            qR(ji,jpk)  =qc(ji,jpk)                 ! re-assumed to be
            qL(ji,jpk)  =qc(ji,jpk)                 ! piecewise constant.
          ENDDO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
          DO jk=1, jpk
            DO ji= IRANGE
              dltR=qR(ji,jk)-qc(ji,jk)
              dltL=qc(ji,jk)-qL(ji,jk)
              cffR=2.*dltR
              cffL=2.*dltL
              IF ((dltR*dltL).LT.0.) THEN
                dltR=0.
                dltL=0.
              ELSE IF (ABS(dltR).GT.ABS(cffL)) THEN
                dltR=cffL
              ELSE IF (ABS(dltL).GT.ABS(cffR)) THEN
                dltL=cffR
              ENDIF
              qR(ji,jk)=qc(ji,jk)+dltR
              qL(ji,jk)=qc(ji,jk)-dltL
            ENDDO
          ENDDO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
          DO jk=1,jpk-1
            DO ji=IRANGE
              cff=rfact2 *ABS(pwsink(ji,jj,jk))/rday*zmask(ji,jj,jk)
              FC(ji,jk+1)=0.
              WL(ji,jk)=-gdepw_n(ji,jj,jk+1)+cff !!! ATTENTION
              WR(ji,jk)=zdept(ji,jj,jk)*qc(ji,jk)
              ksource(ji,jk)=jk
            ENDDO
          ENDDO

          DO jk=1, jpk
            DO ks=2,jk
              DO ji=IRANGE
                IF (WL(ji,jk).GT.-gdepw_n(ji,jj,ks)) THEN
                  ksource(ji,jk)=ks-1
                  FC(ji,jk+1)=FC(ji,jk+1)+WR(ji,ks)
                ENDIF
              ENDDO
            END DO
          END DO
!
!  Finalize computation of flux: add fractional part.
!
          DO jk=1,jpk-1
            DO ji=IRANGE
              ks=ksource(ji,jk)
              cu=MIN(1.,(WL(ji,jk)+gdepw_n(ji,jj,ks+1))*Hz_inv(ji,ks))
              FC(ji,jk+1)=FC(ji,jk+1)+                                      &
     &                  zdept(ji,jj,ks)*cu*                                  &
     &                  (qL(ji,ks)+                                      &
     &                   cu*(0.5*(qR(ji,ks)-qL(ji,ks))-                &
     &                       (1.5-cu)*                               &
     &                       (qR(ji,ks)+qL(ji,ks)-2.*qc(ji,ks))))
            END DO
          END DO
          DO ji=IRANGE
            DO jk=1,jpk-1
               trn(ji,jj,K,jp_tra)=qc(ji,jk)+(FC(ji,jk)-FC(ji,jk+1))*Hz_inv(ji,jk)
               psinkflx(ji,jj,jk)=FC(ji,jk)
!               psinkflx(ji,jj,jk)=(FC(ji,jk)-FC(ji,jk+1))*Hz_inv(ji,jk)
            ENDDO
          ENDDO
        END DO
      END DO

    !
   END SUBROUTINE p4z_sink2_new

   INTEGER FUNCTION p4z_sink_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(2)
      !!----------------------------------------------------------------------
      !
      ierr(:) = 0
      !
      ALLOCATE( sinking (PRIV_2D_BIOARRAY,1:jpk+1), sinking2(PRIV_2D_BIOARRAY,1:jpk+1),     &                
         &      sinkcal (PRIV_2D_BIOARRAY,1:jpk+1), sinksil (PRIV_2D_BIOARRAY,1:jpk+1),     &                
         &      sinkfer2(PRIV_2D_BIOARRAY,1:jpk+1), sinkfer (PRIV_2D_BIOARRAY,1:jpk+1), STAT=ierr(1) )
         !
      IF( ln_p5z    ) ALLOCATE( sinkingn(PRIV_2D_BIOARRAY,1:jpk+1), sinking2n(PRIV_2D_BIOARRAY,1:jpk+1),     &
         &                      sinkingp(PRIV_2D_BIOARRAY,1:jpk+1), sinking2p(PRIV_2D_BIOARRAY,1:jpk+1)   , STAT=ierr(2) )
      !
      p4z_sink_alloc = MAXVAL( ierr )
      IF( p4z_sink_alloc /= 0 ) CALL ctl_warn( 'p4z_sink_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p4z_sink_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sink                    ! Empty routine
   END SUBROUTINE p4z_sink
#endif 
   
   !!======================================================================
END MODULE p4zsink
