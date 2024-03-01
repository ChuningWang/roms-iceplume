! $Id: cpl_prism_put.F 1274 2013-09-16 09:24:57Z gcambon $
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
#if defined OA_COUPLING || defined OW_COUPLING

      SUBROUTINE cpl_prism_put (sndtime)
!     !---------------------------------------------------------------------
!     !              ***  ROUTINE cpl_prism_put  ***
!     !
!     ! ** Purpose : - 
!     !           
!     !
!     !----------------------------------------------------------------------
      USE mod_prism             ! OASIS3 prism module
# ifdef AGRIF
      use Agrif_Util
# endif
      IMPLICIT NONE

# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
# include "ocean2d.h"
# include "mpi_cpl.h"
# include "grid.h"
      
      INTEGER    :: sndtime, i, level, ja
      DOUBLE PRECISION, DIMENSION(Lmmpi,Mmmpi)  :: snd_field 
      CHARACTER*2 clevel
      
# ifdef AGRIF
      level=Agrif_Fixed()
      write(clevel,'(a,i1.1)') "_", level
# else
      level=0
      write(clevel,'(a)') "  "
# endif
!
!!$      if (level==0) then
!!$         MPI_master_only 
!!$     & write(*,*) 'CPL-CROCO : put for parent at time', sndtime
!!$      else		
!!$         MPI_master_only 
!!$     & write(*,*) 'CPL-CROCO : put for child:', level,
!!$     &  'at time', sndtime
!!$      endif
!      
      DO i = 1, nmaxfld
!-------------------------------------------------
!        SST       
         if (trim(ssnd_clname(i))=='CROCO_SST'//trim(clevel)) then
            snd_field = t(1:Lmmpi,1:Mmmpi,N,nnew,itemp) + 273.15    

            DO ja = 0, nmaxatm
               if ( ssnd_nid(ja,i) /= -1 ) then
                  call prism_put_proto(ssnd_nid(ja,i), sndtime, 
     &                  snd_field, comp_ierror)
                  if ( (comp_ierror == OASIS_Sent   ) .or.
     &                 (comp_ierror == OASIS_ToRest ) .or.
     &                 (comp_ierror == OASIS_SentOut) .or.
     &                 (comp_ierror == OASIS_ToRestOut) ) then
                     MPI_master_only write(*,*) 'CPL-CROCO : ',ja,
     &                trim(ssnd_clname(i)),' sent at time ',sndtime
                  endif
               endif  
            enddo
         
         endif
!-------------------------------------------------
!        SSH       
         if (trim(ssnd_clname(i))=='CROCO_SSH'//trim(clevel)) then
            snd_field = zeta(1:Lmmpi,1:Mmmpi,nnew)

            DO ja = 0, nmaxatm
               if ( ssnd_nid(ja,i) /= -1 ) then
                  call prism_put_proto(ssnd_nid(ja,i), sndtime,
     &                  snd_field, comp_ierror)
                  if ( (comp_ierror == OASIS_Sent   ) .or.
     &                 (comp_ierror == OASIS_ToRest ) .or.
     &                 (comp_ierror == OASIS_SentOut) .or.
     &                 (comp_ierror == OASIS_ToRestOut) ) then
                     MPI_master_only write(*,*) 'CPL-CROCO : ',ja,
     &                trim(ssnd_clname(i)),' sent at time ',sndtime
                  endif
               endif
            enddo

         endif
!-------------------------------------------------
!        U component of ocean surface current
!          at rho point
        if (trim(ssnd_clname(i))=='CROCO_UOCE'//trim(clevel)) then
           snd_field = 0.5 * (u(1:Lmmpi  ,1:Mmmpi,N,nnew)
     &                       +u(2:Lmmpi+1,1:Mmmpi,N,nnew))

            DO ja = 0, nmaxatm
               if ( ssnd_nid(ja,i) /= -1 ) then
                  call prism_put_proto(ssnd_nid(ja,i), sndtime,
     &                  snd_field, comp_ierror)
                  if ( (comp_ierror == OASIS_Sent   ) .or.
     &                 (comp_ierror == OASIS_ToRest ) .or.
     &                 (comp_ierror == OASIS_SentOut) .or.
     &                 (comp_ierror == OASIS_ToRestOut) ) then
                     MPI_master_only write(*,*) 'CPL-CROCO : ',ja,
     &                trim(ssnd_clname(i)),' sent at time ',sndtime
                  endif
               endif
            enddo

         endif
!-------------------------------------------------
!        V component of ocean surface current
!          at rho point
         if (trim(ssnd_clname(i))=='CROCO_VOCE'//trim(clevel)) then
           snd_field = 0.5 * (v(1:Lmmpi,1:Mmmpi  ,N,nnew)
     &                       +v(1:Lmmpi,2:Mmmpi+1,N,nnew))

            DO ja = 0, nmaxatm
               if ( ssnd_nid(ja,i) /= -1 ) then
                  call prism_put_proto(ssnd_nid(ja,i), sndtime,
     &                  snd_field, comp_ierror)
                  if ( (comp_ierror == OASIS_Sent   ) .or.
     &                 (comp_ierror == OASIS_ToRest ) .or.
     &                 (comp_ierror == OASIS_SentOut) .or.
     &                 (comp_ierror == OASIS_ToRestOut) ) then
                     MPI_master_only write(*,*) 'CPL-CROCO : ',ja,
     &                trim(ssnd_clname(i)),' sent at time ',sndtime
                  endif
               endif
            enddo
         
         endif
!-------------------------------------------------
!        Eastward component of ocean surface current
!         proceed rotation from xi/eta croco_grid coordinate
!         to East - North coordinate
         if (trim(ssnd_clname(i))=='CROCO_EOCE'//trim(clevel)) then
           snd_field = (0.5 * (u(1:Lmmpi  ,1:Mmmpi,N,nnew)
     &                       + u(2:Lmmpi+1,1:Mmmpi,N,nnew)) )
     &         * cos(angler(1:Lmmpi  ,1:Mmmpi))
     &               - (0.5 * (v(1:Lmmpi,1:Mmmpi  ,N,nnew)
     &                       + v(1:Lmmpi,2:Mmmpi+1,N,nnew)) )
     &         * sin(angler(1:Lmmpi  ,1:Mmmpi))

            DO ja = 0, nmaxatm
               if ( ssnd_nid(ja,i) /= -1 ) then
                  call prism_put_proto(ssnd_nid(ja,i), sndtime,
     &                  snd_field, comp_ierror)
                  if ( (comp_ierror == OASIS_Sent   ) .or.
     &                 (comp_ierror == OASIS_ToRest ) .or.
     &                 (comp_ierror == OASIS_SentOut) .or.
     &                 (comp_ierror == OASIS_ToRestOut) ) then
                     MPI_master_only write(*,*) 'CPL-CROCO : ',ja,
     &                trim(ssnd_clname(i)),' sent at time ',sndtime
                  endif
               endif
            enddo

         endif
!-------------------------------------------------
!        Northward component of ocean surface current
!         proceed rotation from xi/eta croco_grid coordinate
!         to East - North coordinate
         if (trim(ssnd_clname(i))=='CROCO_NOCE'//trim(clevel)) then
           snd_field =  (0.5 * (u(1:Lmmpi  ,1:Mmmpi,N,nnew)
     &                      +   u(2:Lmmpi+1,1:Mmmpi,N,nnew)) )
     &         *  sin(angler(1:Lmmpi  ,1:Mmmpi))
     &               + (0.5 * (v(1:Lmmpi,1:Mmmpi  ,N,nnew)
     &                       + v(1:Lmmpi,2:Mmmpi+1,N,nnew)) )
     &         * cos(angler(1:Lmmpi  ,1:Mmmpi))

            DO ja = 0, nmaxatm
               if ( ssnd_nid(ja,i) /= -1 ) then
                  call prism_put_proto(ssnd_nid(ja,i), sndtime,
     &                  snd_field, comp_ierror)
                  if ( (comp_ierror == OASIS_Sent   ) .or.
     &                 (comp_ierror == OASIS_ToRest ) .or.
     &                 (comp_ierror == OASIS_SentOut) .or.
     &                 (comp_ierror == OASIS_ToRestOut) ) then
                     MPI_master_only write(*,*) 'CPL-CROCO : ',ja,
     &                trim(ssnd_clname(i)),' sent at time ',sndtime
                  endif
               endif
            enddo

         endif
!
      END DO
      return
      end
#else 
!---------------------------------------------------------------------------------------
!     Dummy subroutine                                                           NO Coupling
!---------------------------------------------------------------------------------------
      subroutine cpl_prism_put_empty
      end
#endif 
