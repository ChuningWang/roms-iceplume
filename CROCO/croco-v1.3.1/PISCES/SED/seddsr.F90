#include "cppdefs.h"

MODULE seddsr
   !!======================================================================
   !!              ***  MODULE  seddsr  ***
   !!    Sediment : dissolution and reaction in pore water related 
   !!    related to organic matter
   !!=====================================================================
#if defined key_pisces
   !! * Modules used
   USE sms_pisces, ONLY : rtrn
   USE sed     ! sediment global variable
   USE sed_oce
   USE sedini

   IMPLICIT NONE
   PRIVATE

   PUBLIC sed_dsr

   !!* Substitution
#  include "ocean2pisces.h90"

   !! * Module variables

   REAL(wp) :: zadsnh4
   REAL(wp), DIMENSION(jpsol), PUBLIC      :: dens_mol_wgt  ! molecular density 
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zvolc    ! temp. variables


   !! $Id: seddsr.F90 10362 2018-11-30 15:38:17Z aumont $
CONTAINS
   
   SUBROUTINE sed_dsr( kt, knt ) 
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_dsr  ***
      !! 
      !!  ** Purpose :  computes pore water dissolution and reaction
      !!
      !!  ** Methode :  Computation of the redox reactions in sediment.
      !!                The main redox reactions are solved in sed_dsr whereas
      !!                the secondary reactions are solved in sed_dsr_redoxb.
      !!                A strand spliting approach is being used here (see 
      !!                sed_dsr_redoxb for more information). 
      !!
      !!   History :
      !!        !  98-08 (E. Maier-Reimer, Christoph Heinze )  Original code
      !!        !  04-10 (N. Emprin, M. Gehlen ) f90
      !!        !  06-04 (C. Ethe)  Re-organization
      !!        !  19-08 (O. Aumont) Debugging and improvement of the model.
      !!                             The original method is replaced by a 
      !!                              Strand splitting method which deals 
      !!                              well with stiff reactions.
      !!----------------------------------------------------------------------
      !! Arguments
      INTEGER, INTENT(in) ::   kt, knt       ! number of iteration
      ! --- local variables
      INTEGER :: ji, jk, js, jw, jn   ! dummy looop indices

      REAL(wp), DIMENSION(jpoce,jpksed) :: zrearat1, zrearat2, zrearat3    ! reaction rate in pore water
      REAL(wp), DIMENSION(jpoce,jpksed) :: zundsat    ! undersaturation ; indice jpwatp1 is for calcite   
      REAL(wp), DIMENSION(jpoce,jpksed) :: zkpoc, zkpos, zkpor, zlimo2, zlimno3, zlimso4, zlimfeo    ! undersaturation ; indice jpwatp1 is for calcite   
      REAL(wp), DIMENSION(jpoce)        :: zsumtot
      REAL(wp)  ::  zsolid1, zsolid2, zsolid3, zvolw, zreasat
      REAL(wp)  ::  zsatur, zsatur2, znusil, zkpoca, zkpocb, zkpocc
      REAL(wp)  ::  zratio, zgamma, zbeta, zlimtmp, zundsat2
      !!
      !!----------------------------------------------------------------------
!
      IF( kt == nitsed000 .AND. knt == 1 ) THEN
         IF (lwp) THEN
            WRITE(numsed,*) ' sed_dsr : Dissolution reaction '
            WRITE(numsed,*) ' '
         ENDIF
      ENDIF

     ! Initializations
     !----------------------
      
      zrearat1(:,:)   = 0.    ;   zundsat(:,:) = 0. ; zkpoc(:,:) = 0.
      zlimo2 (:,:)    = 0.    ;   zlimno3(:,:) = 0. ; zrearat2(:,:) = 0.
      zlimso4(:,:)    = 0.    ;   zkpor(:,:)   = 0. ; zrearat3(:,:) = 0.
      zkpos  (:,:)    = 0.
      zsumtot(:)      = rtrn
  
      ALLOCATE( zvolc(jpoce, jpksed, jpsol) )
      zvolc(:,:,:)    = 0.
      zadsnh4 = 1.0 / ( 1.0 + adsnh4 )

      ! Inhibition terms for the different redox equations
      ! --------------------------------------------------
      DO jk = 1, jpksed
         DO ji = 1, jpoce
            zkpoc(ji,jk) = reac_pocl 
            zkpos(ji,jk) = reac_pocs
            zkpor(ji,jk) = reac_pocr
         END DO
      END DO

      ! Conversion of volume units
      !----------------------------
      DO js = 1, jpsol
         DO jk = 1, jpksed
            DO ji = 1, jpoce
               zvolc(ji,jk,js) = ( vols3d(ji,jk) * dens_mol_wgt(js) ) /  &
                  &              ( volw3d(ji,jk) * 1.e-3 )
            ENDDO
         ENDDO
      ENDDO

      !----------------------------------------------------------
      ! 5.  Beginning of solid reaction
      !---------------------------------------------------------
      
      ! Definition of reaction rates [rearat]=sans dim 
      ! For jk=1 no reaction (pure water without solid) for each solid compo
      zrearat1(:,:) = 0.
      zrearat2(:,:) = 0.
      zrearat3(:,:) = 0.

      zundsat(:,:) = pwcp(:,:,jwoxy)

      DO jk = 2, jpksed
         DO ji = 1, jpoce
            zlimo2(ji,jk) = 1.0 / ( zundsat(ji,jk) + xksedo2 )
            zsolid1 = zvolc(ji,jk,jspoc)  * solcp(ji,jk,jspoc)
            zsolid2 = zvolc(ji,jk,jspos)  * solcp(ji,jk,jspos)
            zsolid3 = zvolc(ji,jk,jspor)  * solcp(ji,jk,jspor)
            zkpoca  = zkpoc(ji,jk) * zlimo2(ji,jk)
            zkpocb  = zkpos(ji,jk) * zlimo2(ji,jk)
            zkpocc  = zkpor(ji,jk) * zlimo2(ji,jk)
            zrearat1(ji,jk)  = ( zkpoc(ji,jk) * dtsed2 * zsolid1 ) / &
            &                 ( 1. + zkpoca * zundsat(ji,jk ) * dtsed2 )
            zrearat2(ji,jk)  = ( zkpos(ji,jk) * dtsed2 * zsolid2 ) / &
            &                 ( 1. + zkpocb * zundsat(ji,jk ) * dtsed2 )
            zrearat3(ji,jk)  = ( zkpor(ji,jk) * dtsed2 * zsolid3 ) / &
            &                 ( 1. + zkpocc * zundsat(ji,jk ) * dtsed2 )
         ENDDO
      ENDDO

      ! left hand side of coefficient matrix
!      DO jn = 1, 5
      DO jk = 2, jpksed
         DO ji = 1, jpoce
jflag1:     DO jn = 1, 10
               zsolid1 = zvolc(ji,jk,jspoc)  * solcp(ji,jk,jspoc)
               zsolid2 = zvolc(ji,jk,jspos)  * solcp(ji,jk,jspos)
               zsolid3 = zvolc(ji,jk,jspor)  * solcp(ji,jk,jspor)
               zbeta   = xksedo2 - pwcp(ji,jk,jwoxy) + so2ut * ( zrearat1(ji,jk)    &
               &         + zrearat2(ji,jk) + zrearat3(ji,jk) )
               zgamma = - xksedo2 * pwcp(ji,jk,jwoxy)
               zundsat2 = zundsat(ji,jk)
               zundsat(ji,jk) = ( - zbeta + SQRT( zbeta**2 - 4.0 * zgamma ) ) / 2.0
               zlimo2(ji,jk) = 1.0 / ( zundsat(ji,jk) + xksedo2 )
               zkpoca  = zkpoc(ji,jk) * zlimo2(ji,jk)
               zkpocb  = zkpos(ji,jk) * zlimo2(ji,jk)
               zkpocc  = zkpor(ji,jk) * zlimo2(ji,jk)
               zrearat1(ji,jk)  = ( zkpoc(ji,jk) * dtsed2 * zsolid1 ) / &
               &                 ( 1. + zkpoca * zundsat(ji,jk ) * dtsed2 )
               zrearat2(ji,jk)  = ( zkpos(ji,jk) * dtsed2 * zsolid2 ) / &
               &                 ( 1. + zkpocb * zundsat(ji,jk ) * dtsed2 )
               zrearat3(ji,jk)  = ( zkpor(ji,jk) * dtsed2 * zsolid3 ) / &
               &                 ( 1. + zkpocc * zundsat(ji,jk ) * dtsed2 )
               IF ( ABS( (zundsat(ji,jk)-zundsat2)/(zundsat2+rtrn)) < 1E-8 ) THEN
                  EXIT jflag1
               ENDIF
            END DO jflag1
         END DO
      END DO

      ! New solid concentration values (jk=2 to jksed) for each couple 
      DO jk = 2, jpksed
         DO ji = 1, jpoce
            zreasat = zrearat1(ji,jk) * zlimo2(ji,jk) * zundsat(ji,jk) / zvolc(ji,jk,jspoc)
            solcp(ji,jk,jspoc) = solcp(ji,jk,jspoc) - zreasat
            zreasat = zrearat2(ji,jk) * zlimo2(ji,jk) * zundsat(ji,jk) / zvolc(ji,jk,jspos)
            solcp(ji,jk,jspos) = solcp(ji,jk,jspos) - zreasat
            zreasat = zrearat3(ji,jk) * zlimo2(ji,jk) * zundsat(ji,jk) / zvolc(ji,jk,jspor)
            solcp(ji,jk,jspor) = solcp(ji,jk,jspor) - zreasat
         ENDDO
      ENDDO

      ! New pore water concentrations    
      DO jk = 2, jpksed
         DO ji = 1, jpoce
            ! Acid Silicic 
            pwcp(ji,jk,jwoxy)  = zundsat(ji,jk)
            zreasat = ( zrearat1(ji,jk) + zrearat2(ji,jk) + zrearat3(ji,jk) ) * zlimo2(ji,jk) * zundsat(ji,jk)    ! oxygen         
            ! For DIC
            pwcp(ji,jk,jwdic)  = pwcp(ji,jk,jwdic) + zreasat
            zsumtot(ji) = zsumtot(ji) + zreasat / dtsed2 * volw3d(ji,jk) * 1.e-3 * 86400. * 365. * 1E3
            ! For Phosphate (in mol/l)
            pwcp(ji,jk,jwpo4)  = pwcp(ji,jk,jwpo4) + zreasat * spo4r
            ! For iron (in mol/l)
            pwcp(ji,jk,jwfe2)  = pwcp(ji,jk,jwfe2) + fecratio(ji) * zreasat
            ! For alkalinity
            pwcp(ji,jk,jwalk)  = pwcp(ji,jk,jwalk) + zreasat * ( srno3 * zadsnh4 - 2.* spo4r )
            ! Ammonium
            pwcp(ji,jk,jwnh4)  = pwcp(ji,jk,jwnh4) + zreasat * srno3 * zadsnh4
            ! Ligands
            sedligand(ji,jk)   = sedligand(ji,jk) + ratligc * zreasat - reac_ligc * sedligand(ji,jk)
         ENDDO
      ENDDO

      !--------------------------------------------------------------------
      ! Begining POC denitrification and NO3- diffusion
      ! (indice n°5 for couple POC/NO3- ie solcp(:,:,jspoc)/pwcp(:,:,jwno3))
      !--------------------------------------------------------------------

      zrearat1(:,:) = 0.
      zrearat2(:,:) = 0.
      zrearat3(:,:) = 0.

      zundsat(:,:) = pwcp(:,:,jwno3)

      DO jk = 2, jpksed
         DO ji = 1, jpoce
            zlimno3(ji,jk) = ( 1.0 - pwcp(ji,jk,jwoxy) * zlimo2(ji,jk) ) / ( zundsat(ji,jk) + xksedno3 )
            zsolid1 = zvolc(ji,jk,jspoc) * solcp(ji,jk,jspoc)
            zsolid2 = zvolc(ji,jk,jspos) * solcp(ji,jk,jspos)
            zsolid3 = zvolc(ji,jk,jspor) * solcp(ji,jk,jspor)
            zkpoca = zkpoc(ji,jk) * zlimno3(ji,jk)
            zkpocb = zkpos(ji,jk) * zlimno3(ji,jk)
            zkpocc = zkpor(ji,jk) * zlimno3(ji,jk)
            zrearat1(ji,jk)  = ( zkpoc(ji,jk) * dtsed2 * zsolid1 ) / &
            &                 ( 1. + zkpoca * zundsat(ji,jk ) * dtsed2 )
            zrearat2(ji,jk)  = ( zkpos(ji,jk) * dtsed2 * zsolid2 ) / &
            &                 ( 1. + zkpocb * zundsat(ji,jk ) * dtsed2 )
            zrearat3(ji,jk)  = ( zkpor(ji,jk) * dtsed2 * zsolid3 ) / &
            &                 ( 1. + zkpocc * zundsat(ji,jk ) * dtsed2 )
        END DO
      END DO

!      DO jn = 1, 5
      DO jk = 2, jpksed
         DO ji = 1, jpoce
jflag2:    DO jn = 1, 10
               zlimtmp = ( 1.0 - pwcp(ji,jk,jwoxy) * zlimo2(ji,jk) )
               zsolid1 = zvolc(ji,jk,jspoc) * solcp(ji,jk,jspoc)
               zsolid2 = zvolc(ji,jk,jspos) * solcp(ji,jk,jspos)
               zsolid3 = zvolc(ji,jk,jspor) * solcp(ji,jk,jspor)
               zbeta   = xksedno3 - pwcp(ji,jk,jwno3) + srDnit * ( zrearat1(ji,jk)    &
               &         + zrearat2(ji,jk) + zrearat3(ji,jk) ) * zlimtmp
               zgamma = - xksedno3 * pwcp(ji,jk,jwno3)
               zundsat2 = zundsat(ji,jk)
               zundsat(ji,jk) = ( - zbeta + SQRT( zbeta**2 - 4.0 * zgamma ) ) / 2.0
               zlimno3(ji,jk) = ( 1.0 - pwcp(ji,jk,jwoxy) * zlimo2(ji,jk) ) / ( zundsat(ji,jk) + xksedno3 )
               zkpoca  = zkpoc(ji,jk) * zlimno3(ji,jk)
               zkpocb  = zkpos(ji,jk) * zlimno3(ji,jk)
               zkpocc  = zkpor(ji,jk) * zlimno3(ji,jk)
               zrearat1(ji,jk)  = ( zkpoc(ji,jk) * dtsed2 * zsolid1 ) / &
               &                 ( 1. + zkpoca * zundsat(ji,jk ) * dtsed2 )
               zrearat2(ji,jk)  = ( zkpos(ji,jk) * dtsed2 * zsolid2 ) / &
               &                 ( 1. + zkpocb * zundsat(ji,jk ) * dtsed2 )
               zrearat3(ji,jk)  = ( zkpor(ji,jk) * dtsed2 * zsolid3 ) / &
               &                 ( 1. + zkpocc * zundsat(ji,jk ) * dtsed2 )
               IF ( ABS( (zundsat(ji,jk)-zundsat2)/(zundsat2+rtrn)) < 1E-8 ) THEN
                  EXIT jflag2
               ENDIF
            END DO jflag2
         END DO
      END DO


      ! New solid concentration values (jk=2 to jksed) for each couple 
      DO jk = 2, jpksed
         DO ji = 1, jpoce
            zreasat = zrearat1(ji,jk) * zlimno3(ji,jk) * zundsat(ji,jk) / zvolc(ji,jk,jspoc)
            solcp(ji,jk,jspoc) = solcp(ji,jk,jspoc) - zreasat
            zreasat = zrearat2(ji,jk) * zlimno3(ji,jk) * zundsat(ji,jk) / zvolc(ji,jk,jspos)
            solcp(ji,jk,jspos) = solcp(ji,jk,jspos) - zreasat
            zreasat = zrearat3(ji,jk) * zlimno3(ji,jk) * zundsat(ji,jk) / zvolc(ji,jk,jspor)
            solcp(ji,jk,jspor) = solcp(ji,jk,jspor) - zreasat
         ENDDO
      ENDDO

      ! New dissolved concentrations
      DO jk = 2, jpksed
         DO ji = 1, jpoce
            ! For nitrates
            pwcp(ji,jk,jwno3)  =  zundsat(ji,jk)
            zreasat = ( zrearat1(ji,jk) + zrearat2(ji,jk) + zrearat3(ji,jk) ) * zlimno3(ji,jk) * zundsat(ji,jk)
            ! For DIC
            pwcp(ji,jk,jwdic)  = pwcp(ji,jk,jwdic) + zreasat
            zsumtot(ji) = zsumtot(ji) + zreasat / dtsed2 * volw3d(ji,jk) * 1.e-3 * 86400. * 365. * 1E3
            ! For Phosphate (in mol/l)
            pwcp(ji,jk,jwpo4)  = pwcp(ji,jk,jwpo4) + zreasat * spo4r            
            ! Ligands
            sedligand(ji,jk)   = sedligand(ji,jk) + ratligc * zreasat
            ! For iron (in mol/l)
            pwcp(ji,jk,jwfe2)  = pwcp(ji,jk,jwfe2) + fecratio(ji) * zreasat
            ! For alkalinity
            pwcp(ji,jk,jwalk)  = pwcp(ji,jk,jwalk) + zreasat * ( srDnit + srno3 * zadsnh4 - 2.* spo4r )           
            ! Ammonium
            pwcp(ji,jk,jwnh4)  = pwcp(ji,jk,jwnh4) + zreasat * srno3 * zadsnh4
         ENDDO
      ENDDO

      !--------------------------------------------------------------------
      ! Begining POC iron reduction
      ! (indice nï¿½5 for couple POFe(OH)3 ie solcp(:,:,jspoc)/pwcp(:,:,jsfeo))
      !--------------------------------------------------------------------

      zrearat1(:,:) = 0.
      zrearat2(:,:) = 0.
      zrearat3(:,:) = 0.

      zundsat(:,:) = solcp(:,:,jsfeo)

      DO jk = 2, jpksed
         DO ji = 1, jpoce
            zlimfeo(ji,jk) = ( 1.0 - pwcp(ji,jk,jwoxy) * zlimo2(ji,jk) ) * ( 1.0 - pwcp(ji,jk,jwno3)    &
            &                / ( pwcp(ji,jk,jwno3) + xksedno3 ) ) / ( zundsat(ji,jk) + xksedfeo )
            zsolid1 = zvolc(ji,jk,jspoc) * solcp(ji,jk,jspoc)
            zsolid2 = zvolc(ji,jk,jspos) * solcp(ji,jk,jspos)
            zsolid3 = zvolc(ji,jk,jspor) * solcp(ji,jk,jspor)
            zkpoca = zkpoc(ji,jk) * zlimfeo(ji,jk)
            zkpocb = zkpos(ji,jk) * zlimfeo(ji,jk)
            zkpocc = zkpor(ji,jk) * zlimfeo(ji,jk)
            zrearat1(ji,jk) = ( zkpoc(ji,jk) * dtsed2 * zsolid1 ) / &
            &                    ( 1. + zkpoca * zundsat(ji,jk) * dtsed2 )
            zrearat2(ji,jk) = ( zkpos(ji,jk) * dtsed2 * zsolid2 ) / &
            &                    ( 1. + zkpocb * zundsat(ji,jk) * dtsed2 )
            zrearat3(ji,jk) = ( zkpor(ji,jk) * dtsed2 * zsolid3 ) / &
            &                    ( 1. + zkpocc * zundsat(ji,jk) * dtsed2 )
         END DO
      END DO

!      DO jn = 1, 5
      DO jk = 2, jpksed
         DO ji = 1, jpoce
jflag3:     DO jn = 1, 10
               zlimtmp = ( 1.0 - pwcp(ji,jk,jwoxy) * zlimo2(ji,jk) ) * ( 1.0 - pwcp(ji,jk,jwno3)    &
               &                / ( pwcp(ji,jk,jwno3) + xksedno3 ) )
               zsolid1 = zvolc(ji,jk,jspoc) * solcp(ji,jk,jspoc)
               zsolid2 = zvolc(ji,jk,jspos) * solcp(ji,jk,jspos)
               zsolid3 = zvolc(ji,jk,jspor) * solcp(ji,jk,jspor)
               zreasat = ( zrearat1(ji,jk) + zrearat2(ji,jk) + zrearat3(ji,jk) ) / zvolc(ji,jk,jsfeo)
               zbeta   = xksedfeo - solcp(ji,jk,jsfeo) + 4.0 * zreasat * zlimtmp
               zgamma  = -xksedfeo * solcp(ji,jk,jsfeo)
               zundsat2 = zundsat(ji,jk)
               zundsat(ji,jk) = ( - zbeta + SQRT( zbeta**2 - 4.0 * zgamma ) ) / 2.0
               zlimfeo(ji,jk) = ( 1.0 - pwcp(ji,jk,jwoxy) * zlimo2(ji,jk) ) * ( 1.0 - pwcp(ji,jk,jwno3)    &
               &                / ( pwcp(ji,jk,jwno3) + xksedno3 ) ) / ( zundsat(ji,jk) + xksedfeo )
               zkpoca  = zkpoc(ji,jk) * zlimfeo(ji,jk)
               zkpocb  = zkpos(ji,jk) * zlimfeo(ji,jk)
               zkpocc  = zkpor(ji,jk) * zlimfeo(ji,jk)
               zrearat1(ji,jk) = ( zkpoc(ji,jk) * dtsed2 * zsolid1 ) / &
               &                    ( 1. + zkpoca * zundsat(ji,jk) * dtsed2 )
               zrearat2(ji,jk) = ( zkpos(ji,jk) * dtsed2 * zsolid2 ) / &
               &                    ( 1. + zkpocb * zundsat(ji,jk) * dtsed2 )
               zrearat3(ji,jk) = ( zkpor(ji,jk) * dtsed2 * zsolid3 ) / &
               &                    ( 1. + zkpocc * zundsat(ji,jk) * dtsed2 )
               IF ( ABS( (zundsat(ji,jk)-zundsat2)/( MAX(0.,zundsat2)+rtrn)) < 1E-8 ) THEN
                  EXIT jflag3
               ENDIF
            END DO jflag3
         END DO
      END DO


         ! New solid concentration values (jk=2 to jksed) for each couple 
      DO jk = 2, jpksed
         DO ji = 1, jpoce
            zreasat = zrearat1(ji,jk) * zlimfeo(ji,jk) * zundsat(ji,jk) / zvolc(ji,jk,jspoc)
            solcp(ji,jk,jspoc) = solcp(ji,jk,jspoc) - zreasat
            zreasat = zrearat2(ji,jk) * zlimfeo(ji,jk) * zundsat(ji,jk) / zvolc(ji,jk,jspos)
            solcp(ji,jk,jspos) = solcp(ji,jk,jspos) - zreasat
            zreasat = zrearat3(ji,jk) * zlimfeo(ji,jk) * zundsat(ji,jk) / zvolc(ji,jk,jspor)
            solcp(ji,jk,jspor) = solcp(ji,jk,jspor) - zreasat
         END DO
      END DO

      ! New dissolved concentrations
      DO jk = 2, jpksed
         DO ji = 1, jpoce
            zreasat = ( zrearat1(ji,jk) + zrearat2(ji,jk) + zrearat3(ji,jk) ) * zlimfeo(ji,jk) * zundsat(ji,jk)
            ! For FEOH
            solcp(ji,jk,jsfeo) = zundsat(ji,jk)
            ! For DIC
            pwcp(ji,jk,jwdic)  = pwcp(ji,jk,jwdic) + zreasat
            zsumtot(ji) = zsumtot(ji) + zreasat / dtsed2 * volw3d(ji,jk) * 1.e-3 * 86400. * 365. * 1E3
            ! For Phosphate (in mol/l)
            pwcp(ji,jk,jwpo4)  = pwcp(ji,jk,jwpo4) + zreasat * ( spo4r + 4.0 * redfep )
            ! Ligands
            sedligand(ji,jk)   = sedligand(ji,jk) + ratligc * zreasat
            ! For iron (in mol/l)
            pwcp(ji,jk,jwfe2)  = pwcp(ji,jk,jwfe2) + fecratio(ji) * zreasat
            ! For alkalinity
            pwcp(ji,jk,jwalk)  = pwcp(ji,jk,jwalk) + zreasat * ( srno3 * zadsnh4 - 2.* spo4r ) + 8.0 * zreasat
            ! Ammonium
            pwcp(ji,jk,jwnh4)  = pwcp(ji,jk,jwnh4) + zreasat * srno3 * zadsnh4
            pwcp(ji,jk,jwfe2)  = pwcp(ji,jk,jwfe2) + zreasat * 4.0
         ENDDO
      ENDDO

      !--------------------------------------------------------------------
      ! Begining POC denitrification and NO3- diffusion
      ! (indice nï¿½5 for couple POC/NO3- ie solcp(:,:,jspoc)/pwcp(:,:,jwno3))
      !--------------------------------------------------------------------

      zrearat1(:,:) = 0.
      zrearat2(:,:) = 0.
      zrearat3(:,:) = 0.

      zundsat(:,:) = pwcp(:,:,jwso4)

      DO jk = 2, jpksed
         DO ji = 1, jpoce
            zlimso4(ji,jk) = ( 1.0 - pwcp(ji,jk,jwoxy) * zlimo2(ji,jk) ) * ( 1.0 - pwcp(ji,jk,jwno3)    &
            &                / ( pwcp(ji,jk,jwno3) + xksedno3 ) ) * ( 1. - solcp(ji,jk,jsfeo)  &
            &                / ( solcp(ji,jk,jsfeo) + xksedfeo ) ) / ( zundsat(ji,jk) + xksedso4 )
            zsolid1 = zvolc(ji,jk,jspoc) * solcp(ji,jk,jspoc)
            zsolid2 = zvolc(ji,jk,jspos) * solcp(ji,jk,jspos)
            zsolid3 = zvolc(ji,jk,jspor) * solcp(ji,jk,jspor)
            zkpoca = zkpoc(ji,jk) * zlimso4(ji,jk)
            zkpocb = zkpos(ji,jk) * zlimso4(ji,jk)
            zkpocc = zkpor(ji,jk) * zlimso4(ji,jk)
            zrearat1(ji,jk)  = ( zkpoc(ji,jk) * dtsed2 * zsolid1 ) / &
            &                 ( 1. + zkpoca * zundsat(ji,jk ) * dtsed2 )
            zrearat2(ji,jk)  = ( zkpos(ji,jk) * dtsed2 * zsolid2 ) / &
            &                 ( 1. + zkpocb * zundsat(ji,jk ) * dtsed2 )
            zrearat3(ji,jk)  = ( zkpor(ji,jk) * dtsed2 * zsolid3 ) / &
            &                 ( 1. + zkpocc * zundsat(ji,jk ) * dtsed2 )
        END DO
      END DO
!
!      DO jn = 1, 5 
      DO jk = 2, jpksed
         DO ji = 1, jpoce
jflag4:     DO jn = 1, 10
               zlimtmp = ( 1.0 - pwcp(ji,jk,jwoxy) * zlimo2(ji,jk) ) * ( 1.0 - pwcp(ji,jk,jwno3)    &
               &         / ( pwcp(ji,jk,jwno3) + xksedno3 ) ) * ( 1. - solcp(ji,jk,jsfeo)  &
               &         / ( solcp(ji,jk,jsfeo) + xksedfeo ) ) 
               zsolid1 = zvolc(ji,jk,jspoc) * solcp(ji,jk,jspoc)
               zsolid2 = zvolc(ji,jk,jspos) * solcp(ji,jk,jspos)
               zsolid3 = zvolc(ji,jk,jspor) * solcp(ji,jk,jspor)
               zreasat = ( zrearat1(ji,jk) + zrearat2(ji,jk) + zrearat3(ji,jk) ) 
               zbeta   = xksedso4 - pwcp(ji,jk,jwso4) + 0.5 * zreasat * zlimtmp
               zgamma = - xksedso4 * pwcp(ji,jk,jwso4)
               zundsat2 = zundsat(ji,jk)
               zundsat(ji,jk) = ( - zbeta + SQRT( zbeta**2 - 4.0 * zgamma ) ) / 2.0
               zlimso4(ji,jk) = ( 1.0 - pwcp(ji,jk,jwoxy) * zlimo2(ji,jk) ) * ( 1.0 - pwcp(ji,jk,jwno3)    &
               &                / ( pwcp(ji,jk,jwno3) + xksedno3 ) ) * ( 1. - solcp(ji,jk,jsfeo)  &
               &                / ( solcp(ji,jk,jsfeo) + xksedfeo ) ) / ( zundsat(ji,jk) + xksedso4 )
               zkpoca  = zkpoc(ji,jk) * zlimso4(ji,jk)
               zkpocb  = zkpos(ji,jk) * zlimso4(ji,jk)
               zkpocc  = zkpor(ji,jk) * zlimso4(ji,jk)
               zrearat1(ji,jk)  = ( zkpoc(ji,jk) * dtsed2 * zsolid1 ) / &
               &                 ( 1. + zkpoca * zundsat(ji,jk ) * dtsed2 )
               zrearat2(ji,jk)  = ( zkpos(ji,jk) * dtsed2 * zsolid2 ) / &
               &                 ( 1. + zkpocb * zundsat(ji,jk ) * dtsed2 )
               zrearat3(ji,jk)  = ( zkpor(ji,jk) * dtsed2 * zsolid3 ) / &
               &                 ( 1. + zkpocc * zundsat(ji,jk ) * dtsed2 )
               IF ( ABS( (zundsat(ji,jk)-zundsat2)/(zundsat2+rtrn)) < 1E-8 ) THEN
                  EXIT jflag4
               ENDIF
            END DO jflag4
         END DO
      END DO

     ! New solid concentration values (jk=2 to jksed) for each couple 
      DO jk = 2, jpksed
         DO ji = 1, jpoce
            zreasat = zrearat1(ji,jk) * zlimso4(ji,jk) * zundsat(ji,jk) / zvolc(ji,jk,jspoc)
            solcp(ji,jk,jspoc) = solcp(ji,jk,jspoc) - zreasat
            zreasat = zrearat2(ji,jk) * zlimso4(ji,jk) * zundsat(ji,jk) / zvolc(ji,jk,jspos)
            solcp(ji,jk,jspos) = solcp(ji,jk,jspos) - zreasat
            zreasat = zrearat3(ji,jk) * zlimso4(ji,jk) * zundsat(ji,jk) / zvolc(ji,jk,jspor)
            solcp(ji,jk,jspor) = solcp(ji,jk,jspor) - zreasat
         ENDDO
      ENDDO
!
      ! New dissolved concentrations
      DO jk = 2, jpksed
         DO ji = 1, jpoce
            ! For sulfur
            pwcp(ji,jk,jwh2s)  = pwcp(ji,jk,jwh2s) - ( zundsat(ji,jk) - pwcp(ji,jk,jwso4) ) 
            pwcp(ji,jk,jwso4)  =  zundsat(ji,jk)
            zreasat = ( zrearat1(ji,jk) + zrearat2(ji,jk) + zrearat3(ji,jk) ) * zlimso4(ji,jk) * zundsat(ji,jk)
            ! For DIC
            pwcp(ji,jk,jwdic)  = pwcp(ji,jk,jwdic) + zreasat
            zsumtot(ji) = zsumtot(ji) + zreasat / dtsed2 * volw3d(ji,jk) * 1.e-3 * 86400. * 365. * 1E3
            ! For Phosphate (in mol/l)
            pwcp(ji,jk,jwpo4)  = pwcp(ji,jk,jwpo4) + zreasat * spo4r
            ! Ligands
            sedligand(ji,jk)   = sedligand(ji,jk) + ratligc * zreasat
            ! For iron (in mol/l)
            pwcp(ji,jk,jwfe2)  = pwcp(ji,jk,jwfe2) + fecratio(ji) * zreasat
            ! For alkalinity
            pwcp(ji,jk,jwalk)  = pwcp(ji,jk,jwalk) + zreasat * ( srno3 * zadsnh4 - 2.* spo4r ) + zreasat
            ! Ammonium
            pwcp(ji,jk,jwnh4)  = pwcp(ji,jk,jwnh4) + zreasat * srno3 * zadsnh4
         ENDDO
      ENDDO

      ! Oxydation of the reduced products. Here only ammonium and ODU is accounted for
      ! There are two options here: A simple time splitting scheme and a modified 
     ! Patankar scheme
      ! ------------------------------------------------------------------------------

      call sed_dsr_redoxb

      ! -------------------------------------------------------------- 
      !    4/ Computation of the bioturbation coefficient
      !       This parameterization is taken from Archer et al. (2002)
      ! --------------------------------------------------------------

      DO ji = 1, jpoce
         db(ji,:) = dbiot * zsumtot(ji) * pwcp(ji,1,jwoxy) / (pwcp(ji,1,jwoxy) + 20.E-6)
      END DO



      ! ------------------------------------------------------
      !    Vertical variations of the bioturbation coefficient
      ! ------------------------------------------------------
      IF (ln_btbz) THEN
         DO ji = 1, jpoce
            db(ji,:) = db(ji,:) * exp( -(profsedw(:) / dbtbzsc)**2 ) / (365.0 * 86400.0)
         END DO
      ELSE
         DO jk = 1, jpksed
            IF (profsedw(jk) > dbtbzsc) THEN
                db(:,jk) = 0.0
            ENDIF
         END DO
      ENDIF

      IF (ln_irrig) THEN
         DO jk = 1, jpksed
            DO ji = 1, jpoce
               irrig(ji,jk) = ( 7.63752 - 7.4465 * exp( -0.89603 * zsumtot(ji) ) ) * pwcp(ji,1,jwoxy)  &
               &             / (pwcp(ji,1,jwoxy) + 20.E-6) 
               irrig(ji,jk) = irrig(ji,jk) * exp( -(profsedw(jk) / xirrzsc) )
            END DO
         END DO
      ELSE
         irrig(:,:) = 0.0
      ENDIF

      DEALLOCATE( zvolc )
!      
   END SUBROUTINE sed_dsr

   SUBROUTINE sed_dsr_redoxb
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_dsr_redox  ***
      !! 
      !!  ** Purpose :  computes secondary redox reactions
      !!
      !!  ** Methode :  It uses Strand splitter algorithm proposed by 
      !!                Nguyen et al. (2009) and modified by Wang et al. (2018)
      !!                Basically, each equation is solved analytically when 
      !!                feasible, otherwise numerically at t+1/2. Then 
      !!                the last equation is solved at t+1. The other equations
      !!                are then solved at t+1 starting in the reverse order.
      !!                Ideally, it's better to start from the fastest reaction
      !!                to the slowest and then reverse the order to finish up
      !!                with the fastest one. But random order works well also.
      !!                The scheme is second order, positive and mass 
      !!                conserving. It works well for stiff systems.
      !! 
      !!   History :
      !!        !  18-08 (O. Aumont)  Original code
      !!----------------------------------------------------------------------
      !! Arguments
      ! --- local variables
      INTEGER   ::  ji, jk, jn   ! dummy looop indices

      REAL, DIMENSION(6)  :: zsedtrn, zsedtra
      REAL(wp)  ::  zalpha, zbeta, zgamma, zdelta, zepsi, zsedfer
      !!
      !!----------------------------------------------------------------------

      DO ji = 1, jpoce
         DO jk = 2, jpksed
            zsedtrn(1)  = pwcp(ji,jk,jwoxy)
            zsedtrn(2)  = MAX(0., pwcp(ji,jk,jwh2s) )
            zsedtrn(3)  = pwcp(ji,jk,jwnh4)
            zsedtrn(4)  = MAX(0., pwcp(ji,jk,jwfe2) - sedligand(ji,jk) )
            zsedfer     = MIN(0., pwcp(ji,jk,jwfe2) - sedligand(ji,jk) )
            zsedtrn(5)  = solcp(ji,jk,jsfeo) * zvolc(ji,jk,jsfeo)
            zsedtrn(6)  = solcp(ji,jk,jsfes) * zvolc(ji,jk,jsfes)
            zsedtra(:)  = zsedtrn(:) 

            ! First pass of the scheme. At the end, it is 1st order 
            ! -----------------------------------------------------
            ! Fe + O2
            zalpha = zsedtra(1) - 0.25 * zsedtra(4)
            zbeta  = zsedtra(4) + zsedtra(5) 
            zgamma = pwcp(ji,jk,jwalk) - 2.0 * zsedtra(4)
            zdelta = pwcp(ji,jk,jwpo4) - redfep * zsedtra(4)
            IF ( zalpha == 0. ) THEN
               zsedtra(4) = zsedtra(4) / ( 1.0 + zsedtra(4) * reac_fe2 * dtsed2 / 2.0 )
            ELSE
               zsedtra(4) = ( zsedtra(4) * zalpha ) / ( 0.25 * zsedtra(4) *   &
               &            ( exp( reac_fe2 * zalpha * dtsed2 / 2. ) - 1.0 )  &
               &            + zalpha * exp( reac_fe2 * zalpha * dtsed2 / 2. ) )
            ENDIF
            zsedtra(1) = zalpha + 0.25 * zsedtra(4) 
            zsedtra(5) = zbeta  - zsedtra(4)
            pwcp(ji,jk,jwalk) = zgamma + 2.0 * zsedtra(4)
            pwcp(ji,jk,jwpo4) = zdelta + redfep * zsedtra(4)
            ! H2S + O2
            zalpha = zsedtra(1) - 2.0 * zsedtra(2)
            zbeta  = pwcp(ji,jk,jwso4) + zsedtra(2)
            zgamma = pwcp(ji,jk,jwalk) - 2.0 * zsedtra(2)
            IF ( zalpha == 0. ) THEN
               zsedtra(2) = zsedtra(2) / ( 1.0 + zsedtra(2) * reac_h2s * dtsed2 / 2.0 )
            ELSE
               zsedtra(2) = ( zsedtra(2) * zalpha ) / ( 2.0 * zsedtra(2)   &
               &            * ( exp( reac_h2s * zalpha * dtsed2 / 2. ) - 1.0 )  &
               &            + zalpha * exp( reac_h2s * zalpha * dtsed2 / 2. ) )
            ENDIF
            zsedtra(1) = zalpha + 2.0 * zsedtra(2)
            pwcp(ji,jk,jwalk) = zgamma + 2.0 * zsedtra(2)
            pwcp(ji,jk,jwso4) = zbeta - zsedtra(2)
            ! NH4 + O2
            zalpha = zsedtra(1) - 2.0 * zsedtra(3)
            zgamma = pwcp(ji,jk,jwalk) - 2.0 * zsedtra(3)
            IF ( zalpha == 0. ) THEN
               zsedtra(3) = zsedtra(3) / ( 1.0 + zsedtra(3) * reac_nh4 * zadsnh4 * dtsed2 / 2.0 )
            ELSE
               zsedtra(3) = ( zsedtra(3) * zalpha ) / ( 2.0 * zsedtra(3)   &
               &            * ( exp( reac_nh4 * zadsnh4 * zalpha * dtsed2 / 2. ) - 1.0 )  &
               &            + zalpha * exp( reac_nh4 * zadsnh4 * zalpha * dtsed2 /2. ) )
            ENDIF
            zsedtra(1) = zalpha + 2.0 * zsedtra(3)
            pwcp(ji,jk,jwalk) = zgamma + 2.0 * zsedtra(3)
            ! FeS - O2
            zalpha = zsedtra(1) - 2.0 * zsedtra(6)
            zbeta  = zsedtra(4) + zsedtra(6)
            zgamma = pwcp(ji,jk,jwso4) + zsedtra(6)
            IF ( zalpha == 0. ) THEN
               zsedtra(6) = zsedtra(6) / ( 1.0 + zsedtra(6) * reac_feso * dtsed2 / 2.0 )
            ELSE
               zsedtra(6) = ( zsedtra(6) * zalpha ) / ( 2.0 * zsedtra(6)   &
               &            * ( exp( reac_feso * zalpha * dtsed2 / 2. ) - 1.0 )  &
               &            + zalpha * exp( reac_feso * zalpha * dtsed2 /2. ) )
            ENDIF
            zsedtra(1) = zalpha + 2.0 * zsedtra(6)
            zsedtra(4) = zbeta  - zsedtra(6)
            pwcp(ji,jk,jwso4) = zgamma - zsedtra(6)
!            ! Fe - H2S
            zalpha = zsedtra(2) - zsedtra(4)
            zbeta  = zsedtra(4) + zsedtra(6)
            zgamma = pwcp(ji,jk,jwalk) - 2.0 * zsedtra(4)
            IF ( zalpha == 0. ) THEN
               zsedtra(4) = zsedtra(4) / ( 1.0 + zsedtra(4) * reac_fes * dtsed2 / 2.0 )
            ELSE
               zsedtra(4) = ( zsedtra(4) * zalpha ) / ( zsedtra(4)  &
               &            * ( exp( reac_fes * zalpha * dtsed2 / 2. ) - 1.0 )  &
               &            + zalpha * exp( reac_fes * zalpha * dtsed2 /2. ) )
            ENDIF
            zsedtra(2) = zalpha + zsedtra(4)
            zsedtra(6) = zbeta  - zsedtra(4)
            pwcp(ji,jk,jwalk) = zgamma + 2.0 * zsedtra(4)
            ! FEOH + H2S
            zalpha = zsedtra(5) - 2.0 * zsedtra(2)
            zbeta  = zsedtra(5) + zsedtra(4)
            zgamma = pwcp(ji,jk,jwalk) - 2.0 * zsedtra(4)
            zdelta = pwcp(ji,jk,jwso4) + zsedtra(2)
            zepsi  = pwcp(ji,jk,jwpo4) + redfep * zsedtra(5)
            IF ( zalpha == 0. ) THEN
               zsedtra(2) = zsedtra(2) / ( 1.0 + zsedtra(2) * reac_feh2s * dtsed2 )
            ELSE
               zsedtra(2) = ( zsedtra(2) * zalpha ) / ( 2.0 * zsedtra(2)   &
               &            * ( exp( reac_feh2s * zalpha * dtsed2 ) - 1.0 )  &
               &            + zalpha * exp( reac_feh2s * zalpha * dtsed2 ) )
            ENDIF
            zsedtra(5) = zalpha + 2.0 * zsedtra(2)
            zsedtra(4) = zbeta  - zsedtra(5)
            pwcp(ji,jk,jwso4) = zdelta - zsedtra(2)
            pwcp(ji,jk,jwalk) = zgamma + 2.0 * zsedtra(4)
            pwcp(ji,jk,jwpo4) = zepsi - redfep * zsedtra(5)
            ! Fe - H2S
            zalpha = zsedtra(2) - zsedtra(4)
            zbeta  = zsedtra(4) + zsedtra(6)
            zgamma = pwcp(ji,jk,jwalk) - 2.0 * zsedtra(4)
            IF ( zalpha == 0. ) THEN
               zsedtra(4) = zsedtra(4) / ( 1.0 + zsedtra(4) * reac_fes * dtsed2 / 2.0 )
            ELSE
               zsedtra(4) = ( zsedtra(4) * zalpha ) / ( zsedtra(4) * ( exp( reac_fes * zalpha * dtsed2 / 2. ) - 1.0 )  &
               &            + zalpha * exp( reac_fes * zalpha * dtsed2 /2. ) )
            ENDIF
            zsedtra(2) = zalpha + zsedtra(4)
            zsedtra(6) = zbeta  - zsedtra(4)
            pwcp(ji,jk,jwalk) = zgamma + 2.0 * zsedtra(4)
            ! FeS - O2
            zalpha = zsedtra(1) - 2.0 * zsedtra(6)
            zbeta  = zsedtra(4) + zsedtra(6)
            zgamma = pwcp(ji,jk,jwso4) + zsedtra(6)
            IF (zalpha == 0.) THEN
               zsedtra(6) = zsedtra(6) / ( 1.0 + zsedtra(6) * reac_feso * dtsed2 / 2. )
            ELSE
               zsedtra(6) = ( zsedtra(6) * zalpha ) / ( 2.0 * zsedtra(6)         &
               &            * ( exp( reac_feso * zalpha * dtsed2 / 2. ) - 1.0 )  &
               &            + zalpha * exp( reac_feso * zalpha * dtsed2 /2. ) )
            ENDIF
            zsedtra(1) = zalpha + 2.0 * zsedtra(6)
            zsedtra(4) = zbeta  - zsedtra(6)
            pwcp(ji,jk,jwso4) = zgamma - zsedtra(6)
            ! NH4 + O2
            zalpha = zsedtra(1) - 2.0 * zsedtra(3)
            zgamma = pwcp(ji,jk,jwalk) - 2.0 * zsedtra(3)
            IF (zalpha == 0.) THEN
               zsedtra(3) = zsedtra(3) / ( 1.0 + zsedtra(3) * reac_nh4 * zadsnh4 * dtsed2 / 2.0) 
            ELSE
               zsedtra(3) = ( zsedtra(3) * zalpha ) / ( 2.0 * zsedtra(3)     &
               &            * ( exp( reac_nh4 * zadsnh4 * zalpha * dtsed2 / 2. ) - 1.0 )  &
               &            + zalpha * exp( reac_nh4 * zadsnh4 * zalpha * dtsed2 /2. ) )
            ENDIF
            zsedtra(1) = zalpha + 2.0 * zsedtra(3)
            pwcp(ji,jk,jwalk) = zgamma + 2.0 * zsedtra(3)
            ! H2S + O2
            zalpha = zsedtra(1) - 2.0 * zsedtra(2)
            zbeta  = pwcp(ji,jk,jwso4) + zsedtra(2)
            zgamma = pwcp(ji,jk,jwalk) - 2.0 * zsedtra(2)
            IF ( zalpha == 0. ) THEN
               zsedtra(2) = zsedtra(2) / ( 1.0 + zsedtra(2) * reac_h2s * dtsed2 / 2.0 )
            ELSE
               zsedtra(2) = ( zsedtra(2) * zalpha ) / ( 2.0 * zsedtra(2)   &
               &            * ( exp( reac_h2s * zalpha * dtsed2 / 2. ) - 1.0 )  &
               &            + zalpha * exp( reac_h2s * zalpha * dtsed2 / 2. ) )
            ENDIF
            zsedtra(1) = zalpha + 2.0 * zsedtra(2)
            pwcp(ji,jk,jwso4) = zbeta - zsedtra(2)
            pwcp(ji,jk,jwalk) = zgamma + 2.0 * zsedtra(2)
            ! Fe + O2
            zalpha = zsedtra(1) - 0.25 * zsedtra(4)
            zbeta  = zsedtra(4) + zsedtra(5)
            zgamma = pwcp(ji,jk,jwalk) - 2.0 * zsedtra(4)
            zdelta = pwcp(ji,jk,jwpo4) - redfep * zsedtra(4)
            IF ( zalpha == 0. ) THEN
               zsedtra(4) = zsedtra(4) / ( 1.0 + zsedtra(4) * reac_fe2 * dtsed2 / 2.0 )
            ELSE
               zsedtra(4) = ( zsedtra(4) * zalpha ) / ( 0.25 * zsedtra(4)  &
               &            * ( exp( reac_fe2 * zalpha * dtsed2 / 2. ) - 1.0 )  &
               &            + zalpha * exp( reac_fe2 * zalpha * dtsed2 / 2. ) )
            ENDIF
            zsedtra(1) = zalpha + 0.25 * zsedtra(4)
            zsedtra(5) = zbeta  - zsedtra(4)
            pwcp(ji,jk,jwpo4) = zdelta + redfep * zsedtra(4)
            pwcp(ji,jk,jwalk) = zgamma + 2.0 * zsedtra(4)
            pwcp(ji,jk,jwoxy)  = zsedtra(1)
            pwcp(ji,jk,jwh2s)  = zsedtra(2)
            pwcp(ji,jk,jwnh4)  = zsedtra(3)
            pwcp(ji,jk,jwfe2)  = zsedtra(4) + sedligand(ji,jk) + zsedfer
            pwcp(ji,jk,jwno3)  = pwcp(ji,jk,jwno3)  + ( zsedtrn(3) - pwcp(ji,jk,jwnh4) )
            solcp(ji,jk,jsfeo) = zsedtra(5) / zvolc(ji,jk,jsfeo)
            solcp(ji,jk,jsfes) = zsedtra(6) / zvolc(ji,jk,jsfes)
         END DO
      END DO

  END SUBROUTINE sed_dsr_redoxb

#endif

END MODULE seddsr
