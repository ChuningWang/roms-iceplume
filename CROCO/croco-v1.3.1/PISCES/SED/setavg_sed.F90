! $Id: set_avg.F 1458 2014-02-03 15:01:25Z gcambon $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"

MODULE setavg_sed

#if defined key_pisces

   !! * Modules used
   USE sed
   USE sedarr
   USE sms_pisces, ONLY : rtrn, rfact

   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC set_avg_sed         ! routine called by opa.F90

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "netcdf.inc"

CONTAINS

      SUBROUTINE set_avg_sed

         call set_avg_compute

      RETURN
      END SUBROUTINE set_avg_sed

      SUBROUTINE set_avg_compute
!
! Compute time-averaged fields within a tile.
! ------- ------------- ------ ------ - -----
! Because of syncronization issues, the delayed mode averaging
! procedure is used. This procedure implies that all fields to be
! averaged are sampled during the next time step, rather than at
! the end of the time step when they were computed.
!
! Although this algorithm results in somewhat awkward controlling
! logic it has the advantage that all fields to be sampled
! correspond to exactly the same time, which is time step "n".
! Particularly, this is done this way because vertical velocity
! corresponding to the newly computed horizontal velocities
! becomes available only during the following time step.
! The same applies to the density field.
!
! The algorithm consists of three logical blocks: (1) initialization
! of the averages arrays: when mod(ilc-1,nwrtsedpis_avg).eq.1 the target arrays
! are set to the first contribution; (2) accumulation of averaged
! data, when mod(ilc-1,nwrtsedpis_avg).gt.1; and (3) adding the last
! contribution and scaling.
!
      integer ji,jj, jk, jn, ilc, iout
      real cff, cff1, stf_cff
      REAL(wp) :: zflx, zrate, zsedph
      REAL(wp), DIMENSION(jpoce) :: zflxs
      parameter(stf_cff=86400/0.01)
      integer itrc,k
!
      ilc=1+iic-ntstart   ! number of time step since restart
!
      if (ilc.gt.ntssedpis_avg) then 

        if (mod(ilc-1,nwrtsedpis_avg).eq.1) then
          cff =1.0
          cff1=0.0
        elseif (mod(ilc-1,nwrtsedpis_avg).gt.1) then
          cff =1.0
          cff1=1.0
        elseif (mod(ilc-1,nwrtsedpis_avg).eq.0) then
          cff=1./float(nwrtsedpis_avg)
          cff1=1.0
          if (Istrp+Jstrp == 2) timesedpis_avg=timesedpis_avg+float(nwrtsedpis_avg)*dt
        endif

        do jn = 1, jpsol
           do jk = 1, jpksed
              DO ji = 1, jpoce
                 solcp_avg(ji,jk,jn)=cff*( cff1*solcp_avg(ji,jk,jn)  &
                 &                                +solcp(ji,jk,jn))
              END DO
           END DO
        END DO

        do jn = 1, jpwat
           do jk = 1, jpksed
              DO ji = 1, jpoce
                 pwcp_avg(ji,jk,jn)=cff*( cff1*pwcp_avg(ji,jk,jn)    &
                 &                                +pwcp(ji,jk,jn))
              END DO
           END DO
        END DO

        do jk = 1, jpksed
           do ji = 1, jpoce
              zsedph = -LOG10( hipor(ji,jk) / ( densSW(ji) + rtrn ) + rtrn )
              trcsed_avg(ji,jk,1) = cff*( cff1*trcsed_avg(ji,jk,1)    &
              &                                +zsedph)
              trcsed_avg(ji,jk,2) = cff*( cff1*trcsed_avg(ji,jk,2)    &
              &                                +co3por(ji,jk))
              trcsed_avg(ji,jk,3) = cff*( cff1*trcsed_avg(ji,jk,3)    &
              &                                +sedligand(ji,jk))
           END DO
        END DO 
        DO jn = 1, jpwat
           DO ji = 1, jpoce
              zflx = ( pwcp(ji,1,jn) - pwcp_dta(ji,jn) ) &
              &         * 1.e3 / 1.e2 * dzkbot(ji) / rfact 
              flxsed_avg(ji,jn) = cff*( cff1*flxsed_avg(ji,jn)    &
              &                                +zflx)
           END DO
        END DO

        zflxs(:) = 0.0
        DO jn = 1, jpsol
           zrate =  1.0 / ( denssol * por1(jpksed) ) / rfact
           DO ji = 1, jpoce
              zflxs(ji) = zflxs(ji) + ( tosed(ji,jn) - fromsed(ji,jn) ) * zrate
           ENDDO
        ENDDO

        DO ji = 1, jpoce
           flxsed_avg(ji,jpwatp1) = cff*( cff1*flxsed_avg(ji,jpwatp1)    &
           &                                +zflxs(ji))
        END DO

        DO ji = 1, jpoce
           flxsed_avg(ji,jpdia2dsed) = cff*( cff1*flxsed_avg(ji,jpdia2dsed)    &
           &                                +dzdep(ji)/dtsed)
        END DO

      endif
      return
      end SUBROUTINE set_avg_compute

#endif

END MODULE setavg_sed
