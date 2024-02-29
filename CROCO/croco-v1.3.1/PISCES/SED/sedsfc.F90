#include "cppdefs.h"

MODULE sedsfc
   !!======================================================================
   !!              ***  MODULE  sedsfc  ***
   !!    Sediment : Data at sediment surface
   !!=====================================================================
#if defined key_pisces
   !! * Modules used
   USE sed     ! sediment global variable
   USE sedini
   USE sedarr
   USE seddta

   PUBLIC sed_sfc

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"

   !! $Id: sedsfc.F90 10222 2018-10-25 09:42:23Z aumont $
CONTAINS

   SUBROUTINE sed_sfc( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sed_sfc ***
      !!
      !! ** Purpose :  Give data from sediment model to tracer model
      !!
      !!
      !!   History :
      !!        !  06-04 (C. Ethe)  Orginal code
      !!----------------------------------------------------------------------
      !!* Arguments
      INTEGER, INTENT(in) ::  kt              ! time step

      ! * local variables
      INTEGER :: ji, jj     ! dummy loop indices

      !------------------------------------------------------------------------
      ! reading variables

      CALL unpack_arr ( jpoce, trc_data(PRIV_2D_BIOARRAY,1), iarroce(1:jpoce), pwcp(1:jpoce,1,jwalk) )
      CALL unpack_arr ( jpoce, trc_data(PRIV_2D_BIOARRAY,2), iarroce(1:jpoce), pwcp(1:jpoce,1,jwdic) )
      CALL unpack_arr ( jpoce, trc_data(PRIV_2D_BIOARRAY,3), iarroce(1:jpoce), pwcp(1:jpoce,1,jwno3) )
      CALL unpack_arr ( jpoce, trc_data(PRIV_2D_BIOARRAY,4), iarroce(1:jpoce), pwcp(1:jpoce,1,jwpo4) )
      CALL unpack_arr ( jpoce, trc_data(PRIV_2D_BIOARRAY,5), iarroce(1:jpoce), pwcp(1:jpoce,1,jwoxy) )
      CALL unpack_arr ( jpoce, trc_data(PRIV_2D_BIOARRAY,6), iarroce(1:jpoce), pwcp(1:jpoce,1,jwsil) )
      CALL unpack_arr ( jpoce, trc_data(PRIV_2D_BIOARRAY,7), iarroce(1:jpoce), pwcp(1:jpoce,1,jwnh4) )
      CALL unpack_arr ( jpoce, trc_data(PRIV_2D_BIOARRAY,8), iarroce(1:jpoce), pwcp(1:jpoce,1,jwfe2) )


      DO jj = JRANGE
         DO ji = IRANGE
            IF ( tmask(ji,jj,ikt) == 1 ) THEN
               trb(ji,jj,KSED,jptal) = trc_data(ji,jj,1)
               trb(ji,jj,KSED,jpdic) = trc_data(ji,jj,2)
               trb(ji,jj,KSED,jpno3) = trc_data(ji,jj,3) * redC / redNo3
               trb(ji,jj,KSED,jppo4) = trc_data(ji,jj,4) * redC
               trb(ji,jj,KSED,jpoxy) = trc_data(ji,jj,5)
               trb(ji,jj,KSED,jpsil) = trc_data(ji,jj,6)
               trb(ji,jj,KSED,jpnh4) = trc_data(ji,jj,7) * redC / redNo3
               trb(ji,jj,KSED,jpfer) = trc_data(ji,jj,8)
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE sed_sfc

#endif

END MODULE sedsfc
