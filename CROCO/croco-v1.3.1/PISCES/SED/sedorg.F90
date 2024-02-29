#include "cppdefs.h"

MODULE sedorg
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
   USE seddiff
   USE seddsr

   IMPLICIT NONE
   PRIVATE

   PUBLIC sed_org

   !!* Substitution
#  include "ocean2pisces.h90"

   !! * Module variables

   REAL(wp) :: zadsnh4

   !! $Id: seddsr.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
CONTAINS
   
   SUBROUTINE sed_org( kt ) 
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_org  ***
      !! 
      !!  ** Purpose :  computes pore water diffusion and reaction
      !!
      !!  ** Methode :  Computation of the redox reactions in sediment.
      !!                The main redox reactions are solved in sed_dsr whereas
      !!                the secondary reactions are solved in sed_dsr_redoxb.
      !!                A strand spliting approach is being used here (see 
      !!                sed_dsr_redoxb for more information). 
      !!                Diffusive fluxes are computed in sed_diff
      !!
      !!   History :
      !!        !  98-08 (E. Maier-Reimer, Christoph Heinze )  Original code
      !!        !  04-10 (N. Emprin, M. Gehlen ) f90
      !!        !  06-04 (C. Ethe)  Re-organization
      !!        !  19-08 (O. Aumont) Debugging and improvement of the model.
      !!                             The original method is replaced by a 
      !!                             Strand splitting method which deals 
      !!                             well with stiff reactions.
      !!----------------------------------------------------------------------
      !! Arguments
      INTEGER, INTENT(in) ::   kt
      ! --- local variables
      INTEGER  :: ji, jk, js, jw, jnt   ! dummy looop indices
      REAL(wp) :: zadsnh4
      !!
      !!----------------------------------------------------------------------

!
      IF( kt == nitsed000 ) THEN
         IF (lwp) THEN
            WRITE(numsed,*) ' sed_org : Organic degradation related reactions and diffusion'
            WRITE(numsed,*) ' '
         ENDIF
!         ! 
         dens_mol_wgt(1:jpsol) = denssol / mol_wgt(1:jpsol)
         ! 
      ENDIF

      dtsed2 = dtsed / REAL( nrseddt, wp )

      ! 1. Change of geometry
      !    Increase of dz3d(2) thickness : dz3d(2) = dz3d(2)+dzdep
      !    Warning : no change for dz(2)
      !---------------------------------------------------------
      dz3d(1:jpoce,2) = dz3d(1:jpoce,2) + dzdep(1:jpoce)

      ! New values for volw3d(:,2) and vols3d(:,2)
      ! Warning : no change neither for volw(2) nor  vols(2)
      !------------------------------------------------------
      volw3d(1:jpoce,2) = dz3d(1:jpoce,2) * por(2)
      vols3d(1:jpoce,2) = dz3d(1:jpoce,2) * por1(2)

      ! 2. Change of previous solid fractions (due to volum changes) for k=2
      !---------------------------------------------------------------------

      DO js = 1, jpsol
         DO ji = 1, jpoce
            solcp(ji,2,js) = solcp(ji,2,js) * dz(2) / dz3d(ji,2)
         ENDDO
      END DO

      ! 3. New solid fractions (including solid rain fractions) for k=2
      !------------------------------------------------------------------
      DO js = 1, jpsol
         DO ji = 1, jpoce
            IF (raintg(ji) .ne. 0) THEN
               solcp(ji,2,js) = solcp(ji,2,js) + &
               &           ( rainrg(ji,js) / raintg(ji) ) * ( dzdep(ji) / dz3d(ji,2) )
               ! rainrm are temporary cancel
               rainrm(ji,js) = 0.
            ENDIF
         END DO
      ENDDO

      ! 4.  Adjustment of bottom water concen.(pwcp(1)):
      !     We impose that pwcp(2) is constant. Including dzdep in dz3d(:,2) we assume
      !     that dzdep has got a porosity of por(2). So pore water volum of jk=2 increase.
      !     To keep pwcp(2) cste we must compensate this "increase" by a slight adjusment
      !     of bottom water concentration.
      !     This adjustment is compensate at the end of routine
      !-------------------------------------------------------------
      DO jw = 1, jpwat
         DO ji = 1, jpoce
            pwcp(ji,1,jw) = pwcp(ji,1,jw) - &
               &            pwcp(ji,2,jw) * dzdep(ji) * por(2) / ( dzkbot(ji) + rtrn )
         END DO
      ENDDO

      zadsnh4 = 1.0 / ( 1.0 + adsnh4 )

      ! --------------------------------------------------
      ! Computation of the diffusivities
      ! --------------------------------------------------

      DO js = 1, jpwat
         DO jk = 1, jpksed
            DO ji = 1, jpoce
               diff(ji,jk,js) = ( diff1s(js) + diff2s(js) * temp(ji) ) / ( 1.0 - 2.0 * log( por(jk) ) )
            END DO
         END DO
      END DO

      ! Impact of bioirrigation and adsorption on diffusion
      ! ---------------------------------------------------

      diff(:,:,jwnh4) = diff(:,:,jwnh4) * ( 1.0 + irrig(:,:) ) * zadsnh4
      diff(:,:,jwsil) = diff(:,:,jwsil) * ( 1.0 + irrig(:,:) )
      diff(:,:,jwoxy) = diff(:,:,jwoxy) * ( 1.0 + irrig(:,:) )
      diff(:,:,jwdic) = diff(:,:,jwdic) * ( 1.0 + irrig(:,:) )
      diff(:,:,jwno3) = diff(:,:,jwno3) * ( 1.0 + irrig(:,:) )
      diff(:,:,jwpo4) = diff(:,:,jwpo4) * ( 1.0 + irrig(:,:) )
      diff(:,:,jwalk) = diff(:,:,jwalk) * ( 1.0 + irrig(:,:) )
      diff(:,:,jwh2s) = diff(:,:,jwh2s) * ( 1.0 + irrig(:,:) )
      diff(:,:,jwso4) = diff(:,:,jwso4) * ( 1.0 + irrig(:,:) )
      diff(:,:,jwfe2) = diff(:,:,jwfe2) * ( 1.0 + 0.2 * irrig(:,:) )

      DO jnt = 1, nrseddt
         CALL sed_diff( kt, jnt )        ! 1st pass in diffusion to get values at t+1/2
         CALL sed_dsr ( kt, jnt )        ! Dissolution reaction
         CALL sed_diff( kt, jnt )        ! 2nd pass in diffusion to get values at t+1
      END DO
!      
   END SUBROUTINE sed_org

#endif

END MODULE sedorg
