#include "cppdefs.h"

MODULE sedstp
   !!======================================================================
   !!                       ***  MODULE sedstp   ***
   !!   Sediment model : Sediment model time-stepping
   !!======================================================================
#if defined key_pisces
   USE sed      ! sediment global variables
   USE seddta   ! data read
   USE sedchem  ! chemical constant
   USE sedco3   ! carbonate in sediment pore water
   USE sedorg   ! Organic reactions and diffusion
   USE sedinorg ! Inorganic dissolution
   USE sedbtb   ! bioturbation
   USE sedadv   ! vertical advection
   USE sedmbc   ! mass balance calculation
   USE sedsfc   ! sediment surface data
   USE sedrst   ! restart
   USE sedwri   ! outputs
   USE setavg_sed
   USE sms_pisces, ONLY : rfact

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC sed_stp  ! called by step.F90

   !! $Id: sedstp.F90 10222 2018-10-25 09:42:23Z aumont $
CONTAINS

   SUBROUTINE sed_stp ( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sed_stp  ***
      !!
      !! ** Purpose :   Sediment time stepping
      !!                Simulation of pore water chemistry
      !!
      !! ** Action  :
      !!
      !!
      !!   History :
      !!        !  98-08 (E. Maier-Reimer, Christoph Heinze )  Original code
      !!        !  04-10 (N. Emprin, M. Gehlen ) coupled with PISCES
      !!        !  06-04 (C. Ethe)  Re-organization
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! number of iteration
      INTEGER :: ji,jk,js,jn,jw
      INTEGER :: ilc
      !!----------------------------------------------------------------------
        !
      dtsed  = rfact
!      dtsed2 = dtsed
      IF (kt /= nitsed000) THEN
         CALL sed_dta( kt )       ! Load  Data for bot. wat. Chem and fluxes
      ENDIF

      IF (sedmask == 1. ) THEN
         IF( kt /= nitsed000 )  THEN
           CALL sed_chem( kt )      ! update of chemical constant to account for salinity, temperature changes
         ENDIF

         CALL sed_btb( kt )         ! 1st pass of bioturbation at t+1/2
         CALL sed_org( kt )         ! Organic related reactions and diffusion
         CALL sed_inorg( kt )       ! Dissolution reaction
         CALL sed_btb( kt )         ! 2nd pass of bioturbation at t+1
         tokbot(:,:) = 0.0
         DO jw = 1, jpwat
            DO ji = 1, jpoce
               tokbot(ji,jw) = pwcp(ji,1,jw) * 1.e-3 * dzkbot(ji)
            END DO
         ENDDO
         CALL sed_adv( kt )         ! advection
         CALL sed_co3( kt )         ! pH actualization for saving
         ! This routine is commented out since it does not work at all
!         CALL sed_mbc( kt )         ! cumulation for mass balance calculation

         IF (ln_sed_2way) CALL sed_sfc( kt )         ! Give back new bottom wat chem to tracer model
      ENDIF
#if ! defined XIOS  && defined AVERAGES
      CALL set_avg_sed
      ilc = 1+iic-ntstart   ! number of time step since restart
      IF ( iic > ntstart .AND. mod(ilc-1,nwrtsedpis_avg) == 0 .AND. wrtavg(indxTime) ) THEN
         nrecsedpis_avg=nrecsedpis_avg+1
         CALL sed_wri
      ENDIF
#else
      CALL sed_wri
#endif
      IF( kt == nitsed000 ) THEN
          CALL iom_close( numrsr )       ! close input tracer restart file
      ENDIF
      IF (iic > ntstart .AND. MOD(ilc-1,nrst) == 0) CALL sed_rst_wri

      IF( kt == nitsedend )  CLOSE( numsed )

   END SUBROUTINE sed_stp

#endif

END MODULE sedstp
