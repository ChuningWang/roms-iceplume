#include "cppdefs.h"

MODULE sedadv
   !!======================================================================
   !!              ***  MODULE  sedadv  ***
   !!    Sediment : vertical advection and burial
   !!=====================================================================
#if defined key_pisces
   !! * Modules used
   !!----------------------------------------------------------------------
   !!   sed_adv :
   !!----------------------------------------------------------------------
   USE sed     ! sediment global variable

   IMPLICIT NONE
   PRIVATE

   PUBLIC sed_adv
   PUBLIC sed_adv_alloc

   !!* Substitution
#  include "ocean2pisces.h90"

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) :: dvolsp, dvolsm
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) :: c2por, ckpor

   REAL(wp) :: cpor
   REAL(wp) :: por1clay 
   REAL(wp) :: eps = 1.e-13

   !! $Id: sedadv.F90 10425 2018-12-19 21:54:16Z smasson $
CONTAINS

   SUBROUTINE sed_adv( kt )
      !!-------------------------------------------------------------------------
      !!                  ***  ROUTINE sed_adv  ***
      !!
      !! ** Purpose : vertical solid sediment advection and burial 
      !!
      !! ** Method  : At each grid point the 1-dimensional solid sediment column
      !!              is shifted according the rain added to the top layer and
      !!              the gaps produced through redissolution so that in the end
      !!              the original sediment mixed layer geometry is reestablished.
      !!
      !!
      !!   History :
      !!        !  98-08 (E. Maier-Reimer, Christoph Heinze )  Original code
      !!        !  04-10 (N. Emprin, M. Gehlen ) F90
      !!        !  06-04 (C. Ethe)  Re-organization
      !!-------------------------------------------------------------------------
      !!* Arguments
      INTEGER, INTENT(in) ::  &
         kt                     ! time step
      ! * local variables
      INTEGER :: ji, jk, js 
      INTEGER :: jn, ntimes, nztime, ikwneg
      
      REAL(wp), DIMENSION(jpksed,jpsol) :: zsolcpno
      REAL(wp), DIMENSION(jpksed)       :: zfilled, zfull, zfromup, zempty
      REAL(wp), DIMENSION(jpoce,jpksed) :: zgap, zwb
      REAL(wp), DIMENSION(jpoce,jpsol) :: zrainrf
      REAL(wp), DIMENSION(:  ), ALLOCATABLE :: zraipush

      REAL(wp) :: zkwnup, zkwnlo, zfrac,  zfromce, zrest, sumtot, zsumtot1

      !------------------------------------------------------------------------

!
      IF( kt == nitsed000 ) THEN
         IF (lwp) THEN
            WRITE(numsed,*) ' '
            WRITE(numsed,*) ' sed_adv : vertical sediment advection  '
            WRITE(numsed,*) ' '
         ENDIF
         por1clay = denssol * por1(jpksed) * dz(jpksed)
         cpor     = por1(jpksed) / por1(2)
         DO jk = 2, jpksed
            c2por(jk) = por1(2)      / por1(jk)
            ckpor(jk) = por1(jpksed) / por1(jk)
         ENDDO
         DO jk = jpksedm1, 2, -1
            dvolsp(jk) = vols(jk+1) / vols(jk)
         ENDDO
         DO jk = 3, jpksed
           dvolsm(jk) = vols(jk-1) / vols(jk)
        ENDDO
      ENDIF

      ! Initialization of data for mass balance calculation
      !---------------------------------------------------
      fromsed(:,:) = 0.
      tosed  (:,:) = 0. 
      rloss  (:,:) = 0.
      ikwneg = 1
      nztime = jpksed

      ALLOCATE( zraipush(nztime) )

      ! Initiate gap 
      !--------------
      zgap(:,:) = 0.
      DO js = 1, jpsol
         DO jk = 1, jpksed
            DO ji = 1, jpoce
               zgap(ji,jk) = zgap(ji,jk) + solcp(ji,jk,js)
            END DO
         ENDDO
      ENDDO

      zgap(1:jpoce,1:jpksed) = 1. - zgap(1:jpoce,1:jpksed)   

      ! Initiate burial rates
      !-----------------------
      zwb(:,:) = 0.
      DO jk = 2, jpksed
         zfrac =  dtsed / ( denssol * por1(jk) )     
         DO ji = 1, jpoce
            zwb(ji,jk) = zfrac * raintg(ji)
         END DO
      ENDDO


      DO ji = 1, jpoce
         zwb(ji,2) = zwb(ji,2) - zgap(ji,2) * dz(2)
      ENDDO

      DO jk = 3, jpksed
         zfrac = por1(jk-1) / por1(jk)
         DO ji = 1, jpoce
            zwb(ji,jk) = zwb(ji,jk-1) * zfrac - zgap(ji,jk) * dz(jk)
         END DO
      ENDDO

      zrainrf(:,:) = 0.
      DO ji = 1, jpoce
         IF( raintg(ji) /= 0. )  &
            &   zrainrf(ji,:) = rainrg(ji,:) / raintg(ji)
      ENDDO


      ! Computation of full and empty solid fraction in each layer
      ! for all 'burial' case
      !----------------------------------------------------------


      DO ji = 1, jpoce

         ! computation of total weight fraction in sediment
         !-------------------------------------------------
         zfilled(:) = 0.
         DO js = 1, jpsol
            DO jk = 2, jpksed
               zfilled(jk) = zfilled(jk) + solcp(ji,jk,js)
            ENDDO
         ENDDO
         
         DO js = 1, jpsol
            DO jk = 2, jpksed
               zsolcpno(jk,js) = solcp(ji,jk,js) / zfilled(jk)
            ENDDO
         ENDDO

         ! burial  3 cases: 
         ! zwb > 0 ==> rain > total rection loss 
         ! zwb = 0 ==> rain = 0
         ! zwb < 0 ==> rain > 0 and rain < total reaction loss
         !----------------------------------------------------------------

         IF( zwb(ji,jpksed) > 0. ) THEN

            zfull (jpksed) = zfilled(jpksed)
            zempty(jpksed) = 1. - zfull(jpksed)
            DO jk = jpksedm1, 2, -1
               zfull (jk) = zfilled(jk)
               zfull (jk) = zfull(jk) - zempty(jk+1) * dvolsp(jk)
               zempty(jk) = 1. - zfull(jk)
            ENDDO

            ! Computation of solid sediment species
            !--------------------------------------
            ! push entire sediment column downward to account rest of rain
            DO js = 1, jpsol
               DO jk = jpksed, 3, -1
                  solcp(ji,jk,js) = zfull(jk) * zsolcpno(jk,js) + zempty(jk) * zsolcpno(jk-1,js)
               ENDDO

               solcp(ji,2,js) = zfull(2) * zsolcpno(2,js) + zempty(2) * zrainrf(ji,js)

               DO jk = 2, jpksed
                  zsolcpno(jk,js) = solcp(ji,jk,js)
               END DO
            ENDDO

            zrest = zwb(ji,jpksed) * cpor
            ! what is remaining is less than dz(2)
            IF( zrest <= dz(2) ) THEN           

               zfromup(2) = zrest / dz(2)
               DO jk = 3, jpksed
                  zfromup(jk) = zwb(ji,jpksed) * ckpor(jk) / dz(jk)
               ENDDO
               DO js = 1, jpsol
                  zfromce = 1. - zfromup(2)
                  solcp(ji,2,js) = zfromce * zsolcpno(2,js) + zfromup(2) * zrainrf(ji,js)
                  DO jk = 3, jpksed
                     zfromce = 1. - zfromup(jk)
                     solcp(ji,jk,js) = zfromce * zsolcpno(jk,js) + zfromup(jk) * zsolcpno(jk-1,js)
                  ENDDO
                  fromsed(ji,js) = 0.
                  ! quantities to push in deeper sediment
                  tosed  (ji,js) = zsolcpno(jpksed,js) &
                     &           * zwb(ji,jpksed) * denssol * por1(jpksed)
               ENDDO

            ELSE ! what is remaining is great than dz(2)

               ntimes = INT( zrest / dz(2) ) + 1
               IF( ntimes > nztime ) THEN 
                  CALL ctl_stop( 'sed_adv : rest too large ' )
               ENDIF
               zraipush(1) = dz(2)
               zrest = zrest - zraipush(1)
               DO jn = 2, ntimes
                  IF( zrest >= dz(2) ) THEN
                     zraipush(jn) = dz(2)
                     zrest = zrest - zraipush(jn)
                  ELSE
                     zraipush(jn) = zrest
                     zrest = 0.
                  ENDIF
               ENDDO

               DO jn = 1, ntimes
                  DO js = 1, jpsol
                     DO jk = 2, jpksed
                        zsolcpno(jk,js) = solcp(ji,jk,js)
                     END DO
                  ENDDO
                  
                  zfromup(2) = zraipush(jn) / dz(2)
                  DO jk = 3, jpksed
                     zfromup(jk) = ( zraipush(jn) / dz(jk) ) * c2por(jk)
                  ENDDO

                  DO js = 1, jpsol
                     zfromce = 1. - zfromup(2)
                     solcp(ji,2,js) = zfromce * zsolcpno(2,js) + zfromup(2) * zrainrf(ji,js)
                     DO jk = 3, jpksed
                        zfromce = 1. - zfromup(jk)
                        solcp(ji,jk,js) = zfromce * zsolcpno(jk,js) + zfromup(jk) * zsolcpno(jk-1,js)
                     ENDDO
                     fromsed(ji,js) = 0.
                     tosed  (ji,js) = tosed(ji,js) + zsolcpno(jpksed,js) * zraipush(jn) &
                        &             * denssol * por1(2) 
                  ENDDO
               ENDDO
 
            ENDIF

         ELSE IF( raintg(ji) < eps ) THEN ! rain  = 0
!! Nadia    rloss(:,:) = rainrm(:,:)   bug ??????           

            rloss(ji,1:jpsol) = rainrm(ji,1:jpsol)

            zfull (2) = zfilled(2)
            zempty(2) = 1. - zfull(2)
            DO jk = 3, jpksed
               zfull (jk) = zfilled(jk)
               zfull (jk) = zfull (jk) - zempty(jk-1) * dvolsm(jk)
               zempty(jk) = 1. - zfull(jk)
            ENDDO

            ! fill boxes with weight fraction from underlying box
            DO js = 1, jpsol
               DO jk = 2, jpksedm1
                  solcp(ji,jk,js) = zfull(jk) * zsolcpno(jk,js) + zempty(jk) * zsolcpno(jk+1,js)
               END DO
               solcp(ji,jpksed,js) = zsolcpno(jpksed,js) * zfull(jpksed)
               tosed  (ji,js) = 0.
               fromsed(ji,js) = 0.
            ENDDO
            ! for the last layer, one make go up clay 
            solcp(ji,jpksed,jsclay) = solcp(ji,jpksed,jsclay) + zempty(jpksed) * 1.
            fromsed(ji,jsclay) = zempty(jpksed) * 1. * por1clay
         ELSE  ! rain > 0 and rain < total reaction loss


            DO jk = 2, jpksed
               zfull (jk) = zfilled(jk)
               zempty(jk) = 1. - zfull(jk)
            ENDDO

            ! Determination of indice of layer - ikwneg - where advection is reversed
            !------------------------------------------------------------------------
 iflag:     DO jk = 2, jpksed
               IF( zwb(ji,jk) < 0.  ) THEN
                  ikwneg = jk
                  EXIT iflag
               ENDIF
            ENDDO iflag

            ! computation of zfull and zempty 
            ! 3 cases : a/ ikwneg=2, b/ikwneg=3...jpksedm1, c/ikwneg=jpksed    
            !-------------------------------------------------------------      
            IF( ikwneg == 2 ) THEN ! advection is reversed in the first sediment layer

               zkwnup = rdtsed(ikwneg) * raintg(ji) / dz(ikwneg)
               zkwnlo = ABS( zwb(ji,ikwneg) ) / dz(ikwneg)
               zfull (ikwneg+1) = zfilled(ikwneg+1) - zkwnlo * dvolsm(ikwneg+1)
               zempty(ikwneg+1) = 1. - zfull(ikwneg+1)
               DO jk = ikwneg+2, jpksed
                  zfull (jk) = zfilled(jk) - zempty(jk-1) * dvolsm(jk)
                  zempty(jk) = 1. - zfull(jk)
               ENDDO
               DO js = 1, jpsol
                  solcp(ji,2,js) = zfull(2) * zsolcpno(2,js)+ zkwnlo * zsolcpno(3,js) &
                     &                                      + zkwnup * zrainrf(ji,js)
                  DO jk = 3, jpksedm1
                     solcp(ji,jk,js) = zfull(jk) * zsolcpno(jk,js) + zempty(jk) * zsolcpno(jk+1,js)
                  ENDDO
                  solcp(ji,jpksed,js) = zfull(jpksed) * zsolcpno(jpksed,js)
                  tosed(ji,js)   = 0.
                  fromsed(ji,js) = 0.
               ENDDO
               solcp(ji,jpksed,jsclay) =  solcp(ji,jpksed,jsclay) + zempty(jpksed) * 1.
               !! C. Heinze  fromsed(ji,jsclay) = zempty(jpksed) * 1. * denssol * por1(jpksed) / mol_wgt(jsclay)
               fromsed(ji,jsclay) = zempty(jpksed) * 1. * por1clay
               
            ELSE IF( ikwneg == jpksed ) THEN

               zkwnup = ABS( zwb(ji,ikwneg-1) ) * dvolsm(ikwneg) / dz(ikwneg)
               zkwnlo = ABS( zwb(ji,ikwneg) ) / dz(ikwneg)
               zfull (ikwneg-1) = zfilled(ikwneg-1) - zkwnup * dvolsp(ikwneg-1)
               zempty(ikwneg-1) = 1. - zfull(ikwneg-1)
               DO jk = ikwneg-2, 2, -1
                  zfull (jk) = zfilled(jk) - zempty(jk+1) * dvolsp(jk)
                  zempty(jk) = 1. - zfull(jk) 
               ENDDO
               DO  js = 1, jpsol
                  solcp(ji,2,js) = zfull(2) * zsolcpno(2,js) + zempty(2) * zrainrf(ji,js)
               ENDDO
               DO  js = 1, jpsol
                  DO jk = jpksedm1, 3, -1
                     solcp(ji,jk,js) = zfull(jk) * zsolcpno(jk,js) + zempty(jk) * zsolcpno(jk-1,js)
                  ENDDO
                  solcp(ji,jpksed,js) = zfull(jpksed) * zsolcpno(jpksed,js) &
                     &                       + zkwnup * zsolcpno(jpksedm1,js)
                  tosed(ji,js)   = 0.
                  fromsed(ji,js) = 0.
               ENDDO
               solcp(ji,jpksed,jsclay) = solcp(ji,jpksed,jsclay) + zkwnlo * 1.
               ! Heinze  fromsed(ji,jsclay) = zkwnlo * 1. * denssol * por1(jpksed) / mol_wgt(jsclay)
               fromsed(ji,jsclay) = zkwnlo * 1.* por1clay
            ELSE   ! 2 < ikwneg(ji) <= jpksedm1

               zkwnup = ABS( zwb(ji,ikwneg-1) ) * por1(ikwneg-1) / ( dz(ikwneg) * por1(ikwneg) )
               zkwnlo = ABS( zwb(ji,ikwneg) ) / dz(ikwneg)

               IF( ikwneg > 3 ) THEN

                  zfull (ikwneg-1) = zfilled(ikwneg-1) - zkwnup * dvolsp(ikwneg-1)
                  zempty(ikwneg-1) = 1. - zfull(ikwneg-1) 
                  DO jk = ikwneg-2, 2, -1
                     zfull (jk) = zfilled(jk) - zempty(jk+1) * dvolsp(jk)
                     zempty(jk)    = 1. - zfull(jk) 
                  ENDDO
                  DO js = 1, jpsol
                     solcp(ji,2,js) = zfull(2) * zsolcpno(2,js) + zempty(2) * zrainrf(ji,js)
                  ENDDO
                  DO js = 1, jpsol
                     DO jk = ikwneg-1, 3, -1
                        solcp(ji,jk,js) = zfull(jk) * zsolcpno(jk,js) + zempty(jk) * zsolcpno(jk-1,js)
                     ENDDO
                  ENDDO
               ELSE ! ikw = 3


                  zfull (2) = zfilled(2) - zkwnup * dvolsm(3)
                  zempty(2) = 1. - zfull(2)
                  DO js = 1, jpsol
                     solcp(ji,2,js) = zfull(2) * zsolcpno(2,js) + zempty(2) * zrainrf(ji,js)
                  ENDDO
               ENDIF

               IF( ikwneg < jpksedm1) THEN

                  zfull (ikwneg+1) = zfilled(ikwneg+1) - zkwnlo * dvolsm(ikwneg+1)
                  zempty(ikwneg+1) = 1. - zfull(ikwneg+1) 
                  DO jk = ikwneg+2, jpksed
                     zfull (jk) = zfilled(jk) - zempty(jk-1) * dvolsm(jk)
                     zempty(jk) = 1. - zfull(jk)
                  ENDDO
                  DO js = 1, jpsol
                     DO jk = ikwneg+1, jpksedm1
                        solcp(ji,jk,js) = zfull(jk) * zsolcpno(jk,js) + zempty(jk) * zsolcpno(jk+1,js) 
                     ENDDO
                     solcp(ji,jpksed,js) = zfull(jpksed) * zsolcpno(jpksed,js)
                  ENDDO
                  solcp(ji,jpksed,jsclay) = solcp(ji,jpksed,jsclay) + zempty(jpksed) * 1.
               ELSE

                  zfull (jpksed) = zfilled(jpksed) - zkwnlo * dvolsm(jpksed)
                  zempty(jpksed) = 1. - zfull(jpksed)                
                  DO js = 1, jpsol
                     solcp(ji,jpksed,js) = zfull(jpksed) * zsolcpno(jpksed,js)
                  ENDDO
                  solcp(ji,jpksed,jsclay) = solcp(ji,jpksed,jsclay) + zempty(jpksed) * 1.
               ENDIF ! jpksedm1
               
               ! ikwneg = jpksedm1  ; ikwneg+1 = jpksed ; ikwneg-1 = jpksed - 2
               DO js = 1, jpsol
                  solcp(ji,ikwneg,js) =  zfull(ikwneg) * zsolcpno(ikwneg  ,js) &
                     &                +  zkwnup        * zsolcpno(ikwneg-1,js) &
                     &                +  zkwnlo        * zsolcpno(ikwneg+1,js)   
                  tosed  (ji,js)   = 0.
                  fromsed(ji,js)   = 0.
               ENDDO
               ! Heinze  fromsed(ji,jsclay) = zempty * 1. * denssol * por1(jpksed) / mol_wgt(jsclay)
               fromsed(ji,jsclay) = zempty(jpksed) * 1. * por1clay

            ENDIF ! ikwneg(ji) = 2
         ENDIF  ! zwb > 0
      ENDDO  ! ji = 1, jpoce

      rainrm(:,:) = 0.
      rainrg(:,:) = 0.
      raintg(:)   = 0.

      DEALLOCATE( zraipush )

   END SUBROUTINE sed_adv


   INTEGER FUNCTION sed_adv_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( dvolsp(jpksed), dvolsm(jpksed), c2por(jpksed),         &
      &         ckpor(jpksed) ,           STAT = sed_adv_alloc )
      !
      IF( sed_adv_alloc /= 0 ) CALL ctl_stop( 'sed_adv_alloc : failed to allocate arrays.' )
      !
   END FUNCTION sed_adv_alloc

#endif 
END MODULE sedadv
