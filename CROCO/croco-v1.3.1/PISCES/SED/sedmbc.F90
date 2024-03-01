#include "cppdefs.h"

MODULE sedmbc
   !!======================================================================
   !!              ***  MODULE  sedmbc  ***
   !! Sediment : mass balance calculation
   !!=====================================================================
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   sed_mbc    : 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE sed     ! sediment global variable
   USE seddsr

   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC sed_mbc

   !!* Substitution
#  include "ocean2pisces.h90"

   !! * Module variables
   REAL(wp), DIMENSION(jpsol) :: rain_tot      ! total input rain
   REAL(wp), DIMENSION(jpsol) :: fromsed_tot   ! tota input from sediment
   REAL(wp), DIMENSION(jpsol) :: tosed_tot     ! total output from sediment
   REAL(wp), DIMENSION(jpsol) :: rloss_tot     ! total rain loss

   REAL(wp), DIMENSION(jpwat) :: diss_in_tot   ! total input in pore water
   REAL(wp), DIMENSION(jpwat) :: diss_out_tot  ! total output from pore water

   !! $Id: sedmbc.F90 10250 2018-10-29 13:19:44Z mathiot $
CONTAINS


   SUBROUTINE sed_mbc( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_mbc  ***
      !!
      !! ** Purpose :  computation of total tracer inventories for checking
      !!               mass conservation.
      !!
      !!
      !! ** Method   : tracer inventories of each reservoir are computed and added
      !!               subsequently.
      !!
      !!   History :
      !!        !  04-10  (N. Emprin, M. Gehlen )  Original code
      !!        !  06-07  (C. Ethe)  Re-organization
      !!----------------------------------------------------------------------

      !! Arguments
      INTEGER, INTENT(in) :: kt     ! time step

      !! local declarations
      INTEGER  :: ji,js, jw, jk
      REAL(wp) :: zinit, zfinal 
      REAL(wp) :: zinput, zoutput
      REAL(wp) :: zdsw, zvol
      REAL, DIMENSION(jpsol) :: zsolcp_inv_i, zsolcp_inv_f
      REAL, DIMENSION(jpwat) :: zpwcp_inv_i, zpwcp_inv_f
      REAL(wp) ::  zdelta_sil, zdelta_clay
      REAL(wp) ::  zdelta_co2, zdelta_fe
      REAL(wp) ::  zdelta_po4, zdelta_no3

      !!----------------------------------------------------------------------
      ! Initilization
      !---------------
!
      IF( kt == nitsed000 ) THEN

         DO js = 1, jpsol
            rain_tot   (js) = 0.
            fromsed_tot(js) = 0.
            tosed_tot  (js) = 0.
            rloss_tot  (js) = 0.
         ENDDO

         DO jw = 1, jpwat
            diss_in_tot (jw) = 0.
            diss_out_tot(jw) = 0.
         ENDDO

      ENDIF


      ! Calculation of the cumulativ input and output
      ! for mass balance check
      !----------------------------------------------

      ! cumulativ solid
      DO js = 1, jpsol
         DO ji = 1, jpoce
            ! input [mol]
            rain_tot   (js) = rain_tot   (js) + dtsed * rainrm_dta(ji,js)
            fromsed_tot(js) = fromsed_tot(js) + fromsed(ji,js) / mol_wgt(js)
            ! output [mol]
            tosed_tot  (js) = tosed_tot (js) + tosed(ji,js) / mol_wgt(js)
            rloss_tot  (js) = rloss_tot (js) + rloss(ji,js) / mol_wgt(js)
         ENDDO
      ENDDO

      ! cumulativ dissolved
      DO jw = 1, jpwat
         DO ji = 1, jpoce
            ! input [mol]
            diss_in_tot (jw) = diss_in_tot (jw) + pwcp_dta(ji,jw) * 1.e-3 * dzkbot(ji)
            ! output [mol]
            diss_out_tot(jw) = diss_out_tot(jw) + tokbot(ji,jw)
         ENDDO
      ENDDO

      ! Mass balance check
      !---------------------
      IF( kt == nitsedend ) THEN
         ! initial and final inventories for solid component (mole/dx.dy) in sediment
         zsolcp_inv_i(:) = 0.
         zsolcp_inv_f(:) = 0.
         zpwcp_inv_i (:) = 0.       
         zpwcp_inv_f (:) = 0.        
         DO js = 1, jpsol
            zdsw = denssol / mol_wgt(js)
            DO jk = 2, jpksed
               DO ji = 1, jpoce
                  zvol = vols3d(ji,jk) * zdsw
                  zsolcp_inv_i(js) = zsolcp_inv_i(js) + solcp0(ji,jk,js) * zvol 
                  zsolcp_inv_f(js) = zsolcp_inv_f(js) + solcp (ji,jk,js) * zvol
               ENDDO
            END DO
         ENDDO

         ! initial and final inventories for dissolved component (mole/dx.dy) in sediment
         DO jw = 1, jpwat
            DO jk = 2, jpksed
               DO ji = 1, jpoce 
                  zvol = volw3d(ji,jk) * 1.e-3
                  zpwcp_inv_i(jw) = zpwcp_inv_i(jw) + pwcp0(ji,jk,jw) * zvol
                  zpwcp_inv_f(jw) = zpwcp_inv_f(jw) + pwcp (ji,jk,jw) * zvol
               ENDDO
            END DO
         ENDDO

         ! mass balance for Silica/opal
         zinit      = zsolcp_inv_i(jsopal) + zpwcp_inv_i(jwsil)
         zfinal     = zsolcp_inv_f(jsopal) + zpwcp_inv_f(jwsil)
         zinput     = rain_tot    (jsopal) + diss_in_tot (jwsil)
         zoutput    = tosed_tot   (jsopal) + rloss_tot  (jsopal) + diss_out_tot(jwsil)
         zdelta_sil = ( zfinal + zoutput ) - ( zinit + zinput )


         ! mass balance for Clay
         zinit      = zsolcp_inv_i(jsclay) 
         zfinal     = zsolcp_inv_f(jsclay) 
         zinput     = rain_tot   (jsclay) + fromsed_tot(jsclay) 
         zoutput    = tosed_tot  (jsclay) + rloss_tot  (jsclay) 
         zdelta_clay= ( zfinal + zoutput ) - ( zinit + zinput )

         ! mass balance for carbon ( carbon in POC, CaCo3, DIC )
         zinit      = zsolcp_inv_i(jspoc) + zsolcp_inv_i(jspos) + zsolcp_inv_i(jspor) &
         &          + zsolcp_inv_i(jscal) + zpwcp_inv_i(jwdic)
         zfinal     = zsolcp_inv_f(jspoc) + zsolcp_inv_f(jspos) + zsolcp_inv_f(jspor) &
         &          + zsolcp_inv_f(jscal) + zpwcp_inv_f(jwdic)
         zinput     = rain_tot (jspoc) + rain_tot   (jspos)  +  rain_tot   (jspor) &
         &          + rain_tot (jscal) + diss_in_tot(jwdic)
         zoutput    = tosed_tot(jspoc) + tosed_tot(jspos) + tosed_tot(jspor) + tosed_tot(jscal) + diss_out_tot(jwdic) &
            &       + rloss_tot(jspoc) + rloss_tot(jspos) + rloss_tot(jspor) + rloss_tot(jscal) 
         zdelta_co2 = ( zfinal + zoutput ) - ( zinit + zinput )

         ! mass balance for Sulfur
         zinit      = zpwcp_inv_i(jwso4) + zpwcp_inv_i(jwh2s)   &
         &          + zsolcp_inv_i(jsfes) 
         zfinal     = zpwcp_inv_f(jwso4) + zpwcp_inv_f(jwh2s)   &
         &          + zsolcp_inv_f(jsfes)
         zinput     = diss_in_tot (jwso4) + diss_in_tot (jwh2s) &
         &          + rain_tot (jsfes)
         zoutput    = diss_out_tot(jwso4) + diss_out_tot(jwh2s) &
         &          + tosed_tot(jsfes)    + rloss_tot(jsfes)
         zdelta_no3 = ( zfinal + zoutput ) - ( zinit + zinput )

         ! mass balance for iron
         zinit      = zpwcp_inv_i(jwfe2)  + zsolcp_inv_i(jsfeo)   &
         &          + zsolcp_inv_i(jsfes)
         zfinal     = zpwcp_inv_f(jwfe2)  + zsolcp_inv_f(jsfeo)   &
         &          + zsolcp_inv_f(jsfes)
         zinput     = diss_in_tot (jwfe2) + rain_tot (jsfeo) &
         &          + rain_tot (jsfes)
         zoutput    = diss_out_tot(jwfe2) + tosed_tot(jsfeo) &
         &          + tosed_tot(jsfes)    + rloss_tot(jsfes) + rloss_tot(jsfeo)
         zdelta_fe  = ( zfinal + zoutput ) - ( zinit + zinput )


      END IF

      IF( kt == nitsedend) THEN 

         IF (lwp) THEN
         WRITE(numsed,*)
         WRITE(numsed,*)'==================    General mass balance   ==================  '
         WRITE(numsed,*)' '
         WRITE(numsed,*)' '
         WRITE(numsed,*)' Initial total solid Masses (mole/dx.dy)        '
         WRITE(numsed,*)' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numsed,*)'    Opal,      Clay,       POC,       POS,      POR,        CaCO3,     FeOH,     FeS'
         WRITE(numsed,'(8x,4(1PE10.3,2X))')zsolcp_inv_i(jsopal),zsolcp_inv_i(jsclay),zsolcp_inv_i(jspoc), &
            & zsolcp_inv_i(jspos),zsolcp_inv_i(jspor),zsolcp_inv_i(jscal),zsolcp_inv_i(jsfeo),zsolcp_inv_i(jsfes)
         WRITE(numsed,*)' '
         WRITE(numsed,*)' Initial total dissolved Masses (mole/dx.dy)    '
         WRITE(numsed,*)' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numsed,*)'    Si,         O2,         DIC,        Nit,         Phos,         Fe2+'
         WRITE(numsed,'(5x,5(1PE10.3,2X))') zpwcp_inv_i(jwsil), zpwcp_inv_i(jwoxy), &
            & zpwcp_inv_i(jwdic), zpwcp_inv_i(jwno3), zpwcp_inv_i(jwpo4), zpwcp_inv_i(jwfe2)
         WRITE(numsed,*)' '
         WRITE(numsed,*)'  Solid inputs :  Opale,      Clay,       POC,        CaCO3,        Fe'
         WRITE(numsed,'(A4,10X,5(1PE10.3,2X))')'Rain : ',rain_tot(jsopal),rain_tot(jsclay),rain_tot(jspoc)   &
            & + rain_tot(jspos) + rain_tot(jspor),&
            & rain_tot(jscal), rain_tot(jsfeo)
         WRITE(numsed,'(A12,6x,5(1PE10.3,2X))')' From Sed : ',fromsed_tot(jsopal), fromsed_tot(jsclay), &
            & fromsed_tot(jspoc)+fromsed_tot(jspos)+fromsed_tot(jspor), fromsed_tot(jscal),    &
            & fromsed_tot(jsfeo) + fromsed_tot(jsfes)
         WRITE(numsed,*)'Diss. inputs : Si,    O2,         DIC,         Nit,       Phos,      Fe'
         WRITE(numsed,'(A9,1x,6(1PE10.3,2X))')' From Pisc : ', diss_in_tot(jwsil), &
            & diss_in_tot(jwoxy), diss_in_tot(jwdic), diss_in_tot(jwno3), diss_in_tot(jwpo4), diss_in_tot(jwfe2)
         WRITE(numsed,*)' '
         WRITE(numsed,*)'Solid output : Opale,      Clay,       POC,        CaCO3,        Fe'
         WRITE(numsed,'(A6,8x,5(1PE10.3,2X))')'To sed', tosed_tot(jsopal),tosed_tot(jsclay),tosed_tot(jspoc) &
            & +tosed_tot(jspos)+tosed_tot(jspor),tosed_tot(jscal), tosed_tot(jsfeo)+tosed_tot(jsfes)
         WRITE(numsed,'(A5,9x,5(1PE10.3,2X))')'Perdu', rloss_tot(jsopal),rloss_tot(jsclay),rloss_tot(jspoc) &
            & +rloss_tot(jspos)+rloss_tot(jspor),rloss_tot(jscal),rloss_tot(jsfeo)+rloss_tot(jsfes)
         WRITE(numsed,*)'Diss. output : Si,     O2,        DIC,          Nit,       Phos,        Fe '  
         WRITE(numsed,'(A7,2x,6(1PE10.3,2X))')'To kbot', diss_out_tot(jwsil), &
            & diss_out_tot(jwoxy), diss_out_tot(jwdic), diss_out_tot(jwno3), diss_out_tot(jwpo4), diss_out_tot(jwfe2)
         WRITE(numsed,*)' '
         WRITE(numsed,*)'Final solid  Masses (mole/dx.dy) '
         WRITE(numsed,*)'    Opale,      Clay,       POC,        CaCO3,      Fe'
         WRITE(numsed,'(4x,5(1PE10.3,2X))')zsolcp_inv_f(jsopal),zsolcp_inv_f(jsclay),zsolcp_inv_f(jspoc)  &
            & +zsolcp_inv_f(jspos)+zsolcp_inv_f(jspor),zsolcp_inv_f(jscal),zsolcp_inv_f(jsfeo)+zsolcp_inv_f(jsfes)
         WRITE(numsed,*)' '
         WRITE(numsed,*)'Final dissolved  Masses (mole/dx.dy) (k=2-11)'
         WRITE(numsed,*)'    Si,        O2,         DIC,        Nit,        Phos,        Fe'
         WRITE(numsed,'(4x,6(1PE10.3,2X))') zpwcp_inv_f(jwsil), zpwcp_inv_f(jwoxy), &
            & zpwcp_inv_f(jwdic), zpwcp_inv_f(jwno3), zpwcp_inv_f(jwpo4), zpwcp_inv_f(jwfe2)
         WRITE(numsed,*)' '     
         WRITE(numsed,*)'Delta : Opale,      Clay,       C,         Fe,          S,'
         WRITE(numsed,'(7x,6(1PE11.3,1X))') zdelta_sil / ( zsolcp_inv_i(jsopal) + zpwcp_inv_i(jwsil) ) , &
         &            zdelta_clay / ( zsolcp_inv_i(jsclay) ) ,      & 
         &            zdelta_co2 / ( zsolcp_inv_i(jspoc) + zsolcp_inv_i(jspos) + zsolcp_inv_i(jspor) &
         &          + zsolcp_inv_i(jscal) + zpwcp_inv_i(jwdic) ),     &
         &            zdelta_fe / ( zpwcp_inv_i(jwfe2) + zsolcp_inv_i(jsfeo) + zsolcp_inv_i(jsfes) ) ,  &
         &            zdelta_no3 / ( zpwcp_inv_i(jwso4) + zpwcp_inv_i(jwh2s) + zsolcp_inv_i(jsfes) )
         WRITE(numsed,*)'=========================================================================='

      ENDIF
      ENDIF

   END SUBROUTINE sed_mbc

#endif

END MODULE sedmbc
