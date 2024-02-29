#include "cppdefs.h"

MODULE sedinorg
   !!======================================================================
   !!              ***  MODULE  sedinorg  ***
   !!    Sediment : dissolution and reaction in pore water of 
   !!               inorganic species
   !!=====================================================================
#if defined key_pisces
   !! * Modules used
   USE sms_pisces, ONLY : rtrn
   USE sed     ! sediment global variable
   USE sed_oce
   USE sedmat  ! linear system of equations
   USE sedco3  ! carbonate ion and proton concentration 
   USE sedini
   USE seddsr

   IMPLICIT NONE
   PRIVATE

   PUBLIC sed_inorg

   !!* Substitution
#  include "ocean2pisces.h90"

   !! $Id: seddsr.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
CONTAINS
   
   SUBROUTINE sed_inorg( kt ) 
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_inorg  ***
      !! 
      !!  ** Purpose :  computes pore water dissolution and reaction
      !!
      !!  ** Methode :  implicit simultaneous computation of undersaturation
      !!               resulting from diffusive pore water transport and chemical
      !!               pore water reactions. Solid material is consumed according
      !!               to redissolution and remineralisation
      !!
      !!  ** Remarks :
      !!              - undersaturation : deviation from saturation concentration
      !!              - reaction rate   : sink of undersaturation from dissolution
      !!                                 of solid material 
      !!
      !!   History :
      !!        !  98-08 (E. Maier-Reimer, Christoph Heinze )  Original code
      !!        !  04-10 (N. Emprin, M. Gehlen ) f90
      !!        !  06-04 (C. Ethe)  Re-organization
      !!        !  19-08 (O. Aumont) Debugging and improvement of the model
      !!----------------------------------------------------------------------
      !! Arguments
      INTEGER, INTENT(in) ::   kt       ! number of iteration
      ! --- local variables
      INTEGER :: ji, jk, js, jw         ! dummy looop indices
      REAL(wp), DIMENSION(jpoce,jpksed) :: zrearat1, zrearat2    ! reaction rate in pore water
      REAL(wp), DIMENSION(jpoce,jpksed) :: zundsat    ! undersaturation ; indice jpwatp1 is for calcite   
      REAL(wp), DIMENSION(jpoce) :: zco3eq
      REAL(wp), DIMENSION(jpoce,jpksed,jpsol) :: zvolc    ! temp. variables
      REAL(wp), DIMENSION(jpoce) :: zsieq
      REAL(wp)  ::  zsolid1, zvolw, zreasat
      REAL(wp)  ::  zsatur, zsatur2, znusil, zsolcpcl, zsolcpsi
      !!
      !!----------------------------------------------------------------------

      IF( kt == nitsed000 ) THEN
         IF (lwp) THEN
            WRITE(numsed,*) ' sed_inorg : Dissolution reaction '
            WRITE(numsed,*) ' '
         ENDIF
!         ! 
      ENDIF

     ! Initializations
     !----------------------
      
      zrearat1(:,:) = 0.    ;   zundsat(:,:)  = 0. 
      zrearat2(:,:) = 0.    ;   zrearat2(:,:) = 0.
      zco3eq(:)     = rtrn
      zvolc(:,:,:)  = 0.

      ! -----------------------------------------------
      ! Computation of Si solubility
      ! Param of Ridgwell et al. 2002
      ! -----------------------------------------------

      DO ji = 1, jpoce
         zsolcpcl = 0.0
         zsolcpsi = 0.0
         DO jk = 1, jpksed
            zsolcpsi = zsolcpsi + solcp(ji,jk,jsopal) * dz(jk)
            zsolcpcl = zsolcpcl + solcp(ji,jk,jsclay) * dz(jk)
         END DO
         zsieq(ji) = sieqs(ji) * MAX(0.25, 1.0 - (0.045 * zsolcpcl / zsolcpsi )**0.58 )
         zsieq(ji) = MAX( rtrn, sieqs(ji) )
      END DO

      DO js = 1, jpsol
         DO jk = 1, jpksed
            DO ji = 1, jpoce    
               zvolc(ji,jk,js) = ( vols3d(ji,jk) * dens_mol_wgt(js) ) /  &
                  &              ( volw3d(ji,jk) * 1.e-3  )     
            ENDDO
         ENDDO
      ENDDO

      !----------------------------------------------------------
      ! 5.  Beginning of  Pore Water diffusion and solid reaction
      !---------------------------------------------------------
      
      !-----------------------------------------------------------------------------
      ! For jk=2,jpksed, and for couple 
      !  1 : jwsil/jsopal  ( SI/Opal )
      !  2 : jsclay/jsclay ( clay/clay ) 
      !  3 : jwoxy/jspoc   ( O2/POC )
      !  reaction rate is a function of solid=concentration in solid reactif in [mol/l] 
      !  and undersaturation in [mol/l].
      !  Solid weight fractions should be in ie [mol/l])
      !  second member and solution are in zundsat variable
      !-------------------------------------------------------------------------

      DO jk = 1, jpksed
         DO ji = 1, jpoce
            ! For Silicic Acid and clay
            zundsat(ji,jk) = zsieq(ji) - pwcp(ji,jk,jwsil)
         ENDDO
      ENDDO
      
      ! Definition of reaction rates [rearat]=sans dim 
      ! For jk=1 no reaction (pure water without solid) for each solid compo
      DO ji = 1, jpoce
         zrearat1(ji,:) = 0.
         zrearat2(ji,:) = 0.
      ENDDO

      ! left hand side of coefficient matrix
      DO jk = 2, jpksed
         DO ji = 1, jpoce
            zsolid1 = zvolc(ji,jk,jsopal) * solcp(ji,jk,jsopal)
            zsatur = MAX(0., zundsat(ji,jk) / zsieq(ji) )
            zsatur2 = (1.0 + temp(ji) / 400.0 )**37
            znusil = ( 0.225 * ( 1.0 + temp(ji) / 15.) + 0.775 * zsatur2 * zsatur**2.25 ) / zsieq(ji)
            zrearat1(ji,jk)  = ( reac_sil * znusil * dtsed * zsolid1 ) / &
               &                ( 1. + reac_sil * znusil * dtsed * zundsat(ji,jk) )
         ENDDO
      ENDDO

      CALL sed_mat( jwsil, jpoce, jpksed, zrearat1, zrearat2, zundsat, dtsed )

      ! New solid concentration values (jk=2 to jksed) for each couple 
      DO jk = 2, jpksed
         DO ji = 1, jpoce
            zreasat = zrearat1(ji,jk) * zundsat(ji,jk) / ( zvolc(ji,jk,jsopal) )
            solcp(ji,jk,jsopal) = solcp(ji,jk,jsopal) - zreasat
         ENDDO
      ENDDO

      ! New pore water concentrations    
      DO jk = 1, jpksed
         DO ji = 1, jpoce
            ! Acid Silicic 
            pwcp(ji,jk,jwsil)  = zsieq(ji) - zundsat(ji,jk)
         ENDDO
      ENDDO

      !---------------------------------------------------------------
      ! Performs CaCO3 particle deposition and redissolution (indice 9)
      !--------------------------------------------------------------

      ! computes co3por from the updated pwcp concentrations (note [co3por] = mol/l)

      CALL sed_co3( kt )

      ! *densSW(l)**2 converts aksps [mol2/kg sol2] into [mol2/l2] to get [undsat] in [mol/l]
      DO jk = 1, jpksed
         DO ji = 1, jpoce
            zco3eq(ji)     = aksps(ji) * densSW(ji) * densSW(ji) / ( calcon2(ji) + rtrn )
            zco3eq(ji)     = MAX( rtrn, zco3eq(ji) ) 
            zundsat(ji,jk) = MAX(0., zco3eq(ji) - co3por(ji,jk) )
         ENDDO
      ENDDO

      DO jk = 2, jpksed
         DO ji = 1, jpoce
            zsolid1 = zvolc(ji,jk,jscal) * solcp(ji,jk,jscal)
            zrearat1(ji,jk) = ( reac_cal * dtsed * zsolid1 / zco3eq(ji) ) / &
                  &               ( 1. + reac_cal * dtsed * zundsat(ji,jk) / zco3eq(ji) )
         END DO
      END DO

      ! solves tridiagonal system
      CALL sed_mat( jwdic, jpoce, jpksed, zrearat1, zrearat2, zundsat, dtsed )

      ! New solid concentration values (jk=2 to jksed) for cacO3
      DO jk = 2, jpksed
         DO ji = 1, jpoce
            zreasat = zrearat1(ji,jk) * zundsat(ji,jk) / zvolc(ji,jk,jscal)
            solcp(ji,jk,jscal) = solcp(ji,jk,jscal) - zreasat
         ENDDO
      ENDDO

      ! New dissolved concentrations
      DO jk = 1, jpksed
         DO ji = 1, jpoce
            zreasat = zrearat1(ji,jk) * zundsat(ji,jk)    
            ! For DIC
            pwcp(ji,jk,jwdic)  = pwcp(ji,jk,jwdic) + zreasat
            ! For alkalinity
            pwcp(ji,jk,jwalk)  = pwcp(ji,jk,jwalk) + 2.0 * zreasat 
         ENDDO
      ENDDO

      !-------------------------------------------------
      ! Beginning DIC, Alkalinity 
      !-------------------------------------------------
      
      DO jk = 1, jpksed
         DO ji = 1, jpoce      
            zundsat(ji,jk)   = pwcp(ji,jk,jwdic)
            zrearat1(ji,jk)  = 0.
         ENDDO
      ENDDO

      ! solves tridiagonal system
      CALL sed_mat( jwdic, jpoce, jpksed, zrearat1, zrearat2, zundsat, dtsed )

     ! New dissolved concentrations      
      DO jk = 1, jpksed
         DO ji = 1, jpoce                      
            pwcp(ji,jk,jwdic) = zundsat(ji,jk)
         ENDDO
      ENDDO            

      !-------------------------------------------------
      ! Beginning DIC, Alkalinity 
      !-------------------------------------------------

      DO jk = 1, jpksed
         DO ji = 1, jpoce
            zundsat(ji,jk) = pwcp(ji,jk,jwalk)
            zrearat1(ji,jk) = 0.
         ENDDO
      ENDDO
!
!      ! solves tridiagonal system
      CALL sed_mat( jwalk, jpoce, jpksed, zrearat1, zrearat2, zundsat, dtsed )
!
!      ! New dissolved concentrations      
      DO jk = 1, jpksed
         DO ji = 1, jpoce
            pwcp(ji,jk,jwalk) = zundsat(ji,jk)
         ENDDO
      ENDDO
      
      !----------------------------------
      !   Back to initial geometry
      !-----------------------------
      
      !---------------------------------------------------------------------
      !   1/ Compensation for ajustement of the bottom water concentrations
      !      (see note n° 1 about *por(2))
      !--------------------------------------------------------------------
      DO jw = 1, jpwat
         DO ji = 1, jpoce
            pwcp(ji,1,jw) = pwcp(ji,1,jw) + &
               &            pwcp(ji,2,jw) * dzdep(ji) * por(2) / dzkbot(ji)
         END DO
      ENDDO
      
      !-----------------------------------------------------------------------
      !    2/ Det of new rainrg taking account of the new weight fraction obtained 
      !      in dz3d(2) after diffusion/reaction (react/diffu are also in dzdep!)
      !      This new rain (rgntg rm) will be used in advection/burial routine
      !------------------------------------------------------------------------
      DO js = 1, jpsol
         DO ji = 1, jpoce
            rainrg(ji,js) = raintg(ji) * solcp(ji,2,js)
            rainrm(ji,js) = rainrg(ji,js) / mol_wgt(js)
         END DO
      ENDDO

      !  New raintg
      raintg(:) = 0.
      DO js = 1, jpsol
         DO ji = 1, jpoce
            raintg(ji) = raintg(ji) + rainrg(ji,js)
         END DO
      ENDDO
      
      !--------------------------------
      !    3/ back to initial geometry
      !--------------------------------
      DO ji = 1, jpoce
         dz3d  (ji,2) = dz(2)
         volw3d(ji,2) = dz3d(ji,2) * por(2)
         vols3d(ji,2) = dz3d(ji,2) * por1(2)
      ENDDO
!      
   END SUBROUTINE sed_inorg

#endif

END MODULE sedinorg
