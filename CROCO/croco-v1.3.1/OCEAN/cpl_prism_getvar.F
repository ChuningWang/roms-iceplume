! $Id: cpl_prism_getvar.F 1574 2014-07-01 15:13:38Z gcambon $
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

      SUBROUTINE cpl_prism_getvar(prcv_id, msk, nmaxa, tab_t_tmp, 
     &                                                 rcvtime, iok)    
      
!     !---------------------------------------------------------------------
!     !              ***  ROUTINE cpl_prism_getvar  ***
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
      
      INTEGER,dimension(0:nmaxa):: prcv_id
      INTEGER :: nmaxa
      INTEGER :: rcvtime, comp_ierror, iok, ja
      real,dimension(Lmmpi,Mmmpi)  :: rcv_field 
      real,dimension(GLOBAL_2D_ARRAY) :: tab_t_tmp
      real,dimension(GLOBAL_2D_ARRAY,0:nmaxa) :: msk
!     
!     1.0. Get rcv_field from OASIS
      iok = 0
      DO ja = 0, nmaxa
        if ( prcv_id(ja) /= -1 ) then
          
          call prism_get_proto(prcv_id(ja), rcvtime,
     &         rcv_field(1:Lmmpi,1:Mmmpi), comp_ierror)
          
          if ( (comp_ierror == OASIS_Recvd      ) .or. 
     &         (comp_ierror == OASIS_FromRest   ) .or.  
     &         (comp_ierror == OASIS_RecvOut    ) .or.  
     &         (comp_ierror == OASIS_FromRestOut) ) then
            
!     1.a. Get the fields on the RHO grid
            rcv_field(1:Lmmpi,1:Mmmpi) = rcv_field(1:Lmmpi,1:Mmmpi)
     &           * msk(1:Lmmpi,1:Mmmpi,ja)
            if (iok == 0) then
              MPI_master_only write(*,*)
     &             'CPL-CROCO: getvar first',rcvtime,ja
              tab_t_tmp(1:Lmmpi,1:Mmmpi) = rcv_field(1:Lmmpi,1:Mmmpi)
              iok = comp_ierror
            else
              MPI_master_only write(*,*)
     &             'CPL-CROCO: getvar  next',rcvtime,ja
              tab_t_tmp(1:Lmmpi,1:Mmmpi) = tab_t_tmp(1:Lmmpi,1:Mmmpi)
     &             + rcv_field(1:Lmmpi,1:Mmmpi)
            endif
            
          endif
        endif
      enddo
        
      if ( iok /= 0 ) then
        
!     
!     1.b. Blend the borders of the fields
        if (WESTERN_EDGE) then
           tab_t_tmp(0,:)=2.*tab_t_tmp(1,:)-tab_t_tmp(2,:)
        endif
        if (EASTERN_EDGE) then
           tab_t_tmp(Lmmpi+1,:)=2.*tab_t_tmp(Lmmpi,:)
     &                            -tab_t_tmp(Lmmpi-1,:)
        endif
        if (SOUTHERN_EDGE) then
          tab_t_tmp(:,0)=2.*tab_t_tmp(:,1)-tab_t_tmp(:,2)
        endif
        if (NORTHERN_EDGE) then
          tab_t_tmp(:,Mmmpi+1)=2.*tab_t_tmp(:,Mmmpi)
     &                          - tab_t_tmp(:,Mmmpi-1)
        endif
        
!     !$
!     !$         WRITE(333,*) '****************'
!     !$         WRITE(333,*) 'prism_get_proto: ivarid '  , prcv_id
!     !$         WRITE(333,*) 'prism_get_proto:   kstep', rcvtime
!     !$         WRITE(333,*) 'prism_get_proto:   info ', comp_ierror
!     !$         WRITE(333,*) '     - Minimum value is ', MINVAL(rcv_field)
!     !$         WRITE(333,*) '     - Maximum value is ', MAXVAL(rcv_field)
!     !$         WRITE(333,*) '     -     Sum value is ', SUM(rcv_field)
!     !$         WRITE(333,*) '****************'

      endif
!     
      return
      end 
#else
!---------------------------------------------------------------------------------------
!     Dummy subroutine                                                           NO Coupling
!---------------------------------------------------------------------------------------
      subroutine cpl_prism_getvar_empty
      end
#endif


