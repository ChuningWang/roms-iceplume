#include "cppdefs.h"

MODULE p4zbio
   !!======================================================================
   !!                         ***  MODULE p4zbio  ***
   !! TOP :   PISCES bio-model
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_bio        :   computes the interactions between the different
   !!                      compartments of PISCES
   !!----------------------------------------------------------------------
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zopt          !  optical model
   USE p4zlim          !  Co-limitations of differents nutrients
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p5zprod
   USE p4zmort         !  Mortality terms for phytoplankton
   USE p4zmicro        !  Sources and sinks of microzooplankton
   USE p4zmeso         !  Sources and sinks of mesozooplankton
   USE p5zlim          !  Co-limitations of differents nutrients
   USE p5zmort         !  Mortality terms for phytoplankton
   USE p5zmicro        !  Sources and sinks of microzooplankton
   USE p5zmeso         !  Sources and sinks of mesozooplankton
   USE p4zrem          !  Remineralisation of organic matter
   USE p4zpoc          !  Remineralization of organic particles
   USE p4zagg          !  Aggregation of particles
   USE p4zfechem
   USE p4zligand       !  Prognostic ligand model
  
   IMPLICIT NONE
   PRIVATE

   PUBLIC  p4z_bio    

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"

CONTAINS

   SUBROUTINE p4z_bio ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_bio  ***
      !!
      !! ** Purpose :   Ecosystem model in the whole ocean: computes the
      !!              different interactions between the different compartments
      !!              of PISCES
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      !
      INTEGER             :: ji, jj, jk, jn
      CHARACTER (len=25) :: charout
      LOGICAL :: ltra

      !!---------------------------------------------------------------------

      !     ASSIGN THE SHEAR RATE THAT IS USED FOR AGGREGATION
      !     OF PHYTOPLANKTON AND DETRITUS

      xdiss(:,:,:) = 1.
!!gm the use of nmld should be better here?
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               IF( gdepw_n(ji,jj,jk) > hmld(ji,jj) )   xdiss(ji,jj,jk) = 0.01
            END DO 
         END DO
      END DO

      CALL p4z_opt     ( kt, knt )     ! Optic: PAR in the water column
      CALL p4z_sink    ( kt, knt )     ! vertical flux of particulate organic matter
      CALL p4z_fechem  ( kt, knt )     ! Iron chemistry/scavenging
      IF ( ln_p4z ) THEN
         CALL p4z_lim     ( kt, knt )     ! co-limitations by the various nutrients
         CALL p4z_prod    ( kt, knt )     ! phytoplankton growth rate over the global ocean. 
         CALL p4z_mort    ( kt      )     ! phytoplankton mortality
         CALL p4z_micro   ( kt, knt )     ! microzooplankton
         CALL p4z_meso    ( kt, knt )     ! mesozooplankton
      ELSE
         CALL p5z_lim  ( kt, knt )     ! co-limitations by the various nutrients
         CALL p5z_prod ( kt, knt )     ! phytoplankton growth rate over the global ocean.
         !                             ! (for each element : C, Si, Fe, Chl )
         CALL p5z_mort ( kt      )     ! phytoplankton mortality
         !                             ! zooplankton sources/sinks routines
         CALL p5z_micro( kt, knt )           ! microzooplankton
         CALL p5z_meso ( kt, knt )           ! mesozooplankton
      ENDIF
      !                             ! (for each element : C, Si, Fe, Chl )
      CALL p4z_agg     ( kt, knt )     ! Aggregation of particles
      CALL p4z_rem     ( kt, knt )     ! remineralization terms of organic matter+scavenging of Fe
      CALL p4z_poc     ( kt, knt )     ! Remineralization of organic particles
      !
      IF( ln_ligand )  &
      & CALL p4z_ligand( kt, knt )

      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('bio ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( charout, ltra='tra')
!         CALL prt_ctl_trc(tab4d=trn, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
   END SUBROUTINE p4z_bio

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_bio                         ! Empty routine
   END SUBROUTINE p4z_bio
#endif 

   !!======================================================================
END MODULE  p4zbio
