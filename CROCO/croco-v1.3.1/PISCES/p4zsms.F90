#include "cppdefs.h"

MODULE p4zsms
        !!======================================================================
   !!                         ***  MODULE p4zsms  ***
   !! TOP :   PISCES Source Minus Sink manager
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   p4z_sms        : Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE sms_pisces
   USE p4zint          ! 
   USE p4zche          ! 
   USE p4zbio          ! 
   USE p4zsed          ! 
   USE p4zlys          ! 
   USE p4zflx          ! 
   USE p4zsbc          !
   USE trcwri_pisces          ! 
   USE sedmodel        ! Sediment model

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sms_init   ! called in p4zsms.F90
   PUBLIC   p4z_sms        ! called in p4zsms.F90

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xnegtr     ! Array used to indicate negative tracer values

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zsms.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_sms( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sms  ***
      !!
      !! ** Purpose :   Managment of the call to Biological sources and sinks
      !!              routines of PISCES bio-model
      !!
      !! ** Method  : - at each new day ...
      !!              - several calls of bio and sed ???
      !!              - ...
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !!
      INTEGER  :: ji, jj, jk, jn, jnt, jl
      REAL(wp) :: ztra
      CHARACTER (len=25) :: charout
      LOGICAL :: ltra
      !!---------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN
         !
         ALLOCATE( xnegtr(PRIV_3D_BIOARRAY) )
         !
         DO jn = jp_pcs0, jp_pcs1
           tra(:,:,:,jn) = 0.e0
         ENDDO
         !
      ENDIF

      IF( ll_sbc ) CALL p4z_sbc( kt )   ! external sources of nutrients

      CALL p4z_che          ! computation of chemical constants
      CALL p4z_int          ! computation of various rates for biogeochemistry
      !
      DO jnt = 1, nrdttrc          ! Potential time splitting if requested
         !
         CALL p4z_bio (kt, jnt )    ! Compute soft tissue production (POC)
         CALL p4z_lys( kt, jnt )    ! Compute CaCO3 saturation
         CALL p4z_sed (kt, jnt )    ! compute soft tissue remineralisation
         CALL p4z_flx( kt, jnt )    ! Compute surface fluxes
         !
         !                             ! test if tracers concentrations fall below 0.
         xnegtr(:,:,:) = 1.e0
         DO jn = jp_pcs0, jp_pcs1
            DO jk = KRANGE
               DO jj = JRANGE
                  DO ji = IRANGE
                     IF( ( trb(ji,jj,K,jn) + tra(ji,jj,jk,jn) ) < 0.e0 ) THEN 
                        ztra             = ABS(  ( trb(ji,jj,K,jn) - rtrn ) &
                                               / ( tra(ji,jj,jk,jn) + rtrn ) )
                        xnegtr(ji,jj,jk) = MIN( xnegtr(ji,jj,jk),  ztra )
                     ENDIF
                  END DO
               END DO
            END DO
         END DO
         !                                ! where at least 1 tracer concentration becomes negative
         !                                ! 
         DO jn = jp_pcs0, jp_pcs1
            DO jk = KRANGE
               DO jj = JRANGE
                  DO ji = IRANGE
                     trb(ji,jj,K,jn) = trb(ji,jj,K,jn) & 
                                     + xnegtr(ji,jj,jk) * tra(ji,jj,jk,jn)
                     tra(ji,jj,jk,jn) = 0.e0
                  END DO
               END DO
            END DO
         END DO
         !
      END DO

      IF( ln_sediment ) THEN
         !
         CALL sed_model( kt )     !  Main program of Sediment model
         !
      ENDIF

      IF( ln_ctl )  THEN
        CALL prt_ctl_trc( 'sms' ) 
        tra_ctl(:) = 0.
      ENDIF

      IF( kt - 1 == nitend )  CALL  tracer_stat( kt )

   END SUBROUTINE p4z_sms

   SUBROUTINE p4z_sms_init
      !!----------------------------------------------------------------------
      !!                     ***  p4z_sms_init  ***  
      !!
      !! ** Purpose :   read PISCES namelist
      !!
      !! ** input   :   file 'namelist.trc.s' containing the following
      !!             namelist: natext, natbio, natsms
      !!----------------------------------------------------------------------
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!
      NAMELIST/nampisbio/ nrdttrc, wsbio, xkmort, ferat3, wsbio2, wsbio2max, wsbio2scale,   &
                  &        ln_sink_new, niter1max, niter2max, ldocp, ldocz, lthet,  &
                  &        no3rat3, po4rat3
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_sms_init : PISCES initialization'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF

      REWIND( numnatp_ref )              ! Namelist nampisbio in reference namelist : Pisces variables
      READ  ( numnatp_ref, nampisbio, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampisbio in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampisbio in configuration namelist : Pisces variables
      READ  ( numnatp_cfg, nampisbio, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nampisbio in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, nampisbio )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : nampisbio'
         WRITE(numout,*) '      frequency for the biology                 nrdttrc     =', nrdttrc
         WRITE(numout,*) '      POC sinking speed                         wsbio       =', wsbio
         WRITE(numout,*) '      half saturation constant for mortality    xkmort      =', xkmort
         IF( ln_p5z ) THEN
            WRITE(numout,*) '      N/C in zooplankton                        no3rat3     =', no3rat3
            WRITE(numout,*) '      P/C in zooplankton                        po4rat3     =', po4rat3
         ENDIF
         WRITE(numout,*) '      Fe/C in zooplankton                       ferat3      =', ferat3
         WRITE(numout,*) '    Use of new sinking scheme (y/n)           ln_sink_new   =', ln_sink_new
         WRITE(numout,*) '      Big particles sinking speed               wsbio2      =', wsbio2
         WRITE(numout,*) '      Big particles maximum sinking speed       wsbio2max   =', wsbio2max
         WRITE(numout,*) '      Big particles sinking speed length scale  wsbio2scale =', wsbio2scale
         WRITE(numout,*) '    Maximum number of iterations for POC      niter1max   =', niter1max
         WRITE(numout,*) '    Maximum number of iterations for GOC      niter2max   =', niter2max
         IF( ln_ligand ) THEN
            IF( ln_p4z ) THEN
               WRITE(numout,*) '      Phyto ligand production per unit doc           ldocp  =', ldocp
               WRITE(numout,*) '      Zoo ligand production per unit doc             ldocz  =', ldocz
               WRITE(numout,*) '      Proportional loss of ligands due to Fe uptake  lthet  =', lthet
            ENDIF
         ENDIF
      ENDIF

   END SUBROUTINE p4z_sms_init


#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sms                   ! Empty routine
   END SUBROUTINE p4z_sms
#endif 

   !!======================================================================
END MODULE p4zsms
