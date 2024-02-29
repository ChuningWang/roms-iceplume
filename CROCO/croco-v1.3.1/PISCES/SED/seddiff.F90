#include "cppdefs.h"

MODULE seddiff
   !!======================================================================
   !!              ***  MODULE  seddsr  ***
   !!    Sediment : dissolution and reaction in pore water related 
   !!    related to organic matter
   !!=====================================================================
#if defined key_pisces
   !! * Modules used
   USE sed     ! sediment global variable
   USE sed_oce
   USE sedmat  ! linear system of equations
   USE sedini

   IMPLICIT NONE
   PRIVATE

   PUBLIC sed_diff

   !!* Substitution
#  include "ocean2pisces.h90"

   !! * Module variables

   !! $Id: seddsr.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
CONTAINS
   
   SUBROUTINE sed_diff( kt, knt ) 
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_diff  ***
      !! 
      !!  ** Purpose :  computes pore water diffusion
      !!
      !!  ** Methode :  implicit computation of undersaturation
      !!               resulting from diffusive pore water transport.
      !!
      !!  ** Remarks :
      !!              - undersaturation : deviation from saturation concentration
      !!   History :
      !!        !  98-08 (E. Maier-Reimer, Christoph Heinze )  Original code
      !!        !  04-10 (N. Emprin, M. Gehlen ) f90
      !!        !  06-04 (C. Ethe)  Re-organization
      !!        !  19-08 (O. Aumont) Debugging and improvement of the model
      !!----------------------------------------------------------------------
      !! Arguments
      INTEGER, INTENT(in) ::   kt, knt       ! number of iteration
      ! --- local variables
      INTEGER :: ji, jk, js   ! dummy looop indices

      REAL(wp), DIMENSION(jpoce,jpksed) :: zrearat1, zrearat2   ! reaction rate in pore water
      !!
      !!----------------------------------------------------------------------

      IF( kt == nitsed000 .AND. knt == 1 ) THEN
         IF (lwp) THEN
            WRITE(numsed,*) ' sed_diff : pore-water diffusion '
            WRITE(numsed,*) ' '
         ENDIF
      ENDIF

     ! Initializations
     !----------------------
      zrearat1(:,:)   = 0.
      zrearat2(:,:) = 0.

      !---------------------------
      ! Solves PO4 diffusion 
      !----------------------------

      ! solves tridiagonal system
      CALL sed_mat( jwpo4, jpoce, jpksed, zrearat1, zrearat2, pwcp(:,:,jwpo4), dtsed2 / 2.0 )

      !---------------------------
      ! Solves NH4 diffusion 
      !----------------------------

      ! solves tridiagonal system
      CALL sed_mat( jwnh4, jpoce, jpksed, zrearat1, zrearat2, pwcp(:,:,jwnh4), dtsed2 / 2.0 )

      !---------------------------
      ! Solves Fe2+ diffusion 
      !----------------------------

      ! solves tridiagonal system
      CALL sed_mat( jwfe2, jpoce, jpksed, zrearat1, zrearat2, pwcp(:,:,jwfe2), dtsed2 / 2.0 )

      !---------------------------
      ! Solves H2S diffusion 
      !----------------------------

      ! solves tridiagonal system
      CALL sed_mat( jwh2s, jpoce, jpksed, zrearat1, zrearat2, pwcp(:,:,jwh2s), dtsed2 / 2.0  )

      !---------------------------
      ! Solves SO4 diffusion 
      !----------------------------

      ! solves tridiagonal system
      CALL sed_mat( jwso4, jpoce, jpksed, zrearat1, zrearat2, pwcp(:,:,jwso4), dtsed2 / 2.0 )

      !---------------------------
      ! Solves O2 diffusion
      !----------------------------

      ! solves tridiagonal system
      CALL sed_mat( jwoxy, jpoce, jpksed, zrearat1, zrearat2, pwcp(:,:,jwoxy), dtsed2 / 2.0 )

      !---------------------------
      ! Solves NO3 diffusion
      !----------------------------

      ! solves tridiagonal system
      CALL sed_mat( jwno3, jpoce, jpksed, zrearat1, zrearat2, pwcp(:,:,jwno3), dtsed2 / 2.0 )

      CALL sed_mat( jwdic, jpoce, jpksed, zrearat1, zrearat2, sedligand(:,:), dtsed2 / 2.0 )

!      
   END SUBROUTINE sed_diff

#endif

END MODULE seddiff
