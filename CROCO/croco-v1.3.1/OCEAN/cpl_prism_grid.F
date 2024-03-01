! $Id: cpl_prism_grid.F Joris Pianezze 9 aout 2016 $
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

      SUBROUTINE cpl_prism_grid(id_partition)

      USE mod_prism             ! OASIS3 prism module
# ifdef AGRIF
      use Agrif_Util
# endif
      IMPLICIT NONE
      include 'mpif.h'
      
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "mpi_cpl.h"
# include "set_global_definitions.h"
!
      integer, intent(in) :: id_partition
!
      integer error_flg
!
      real xbox, ybox
      real clonr(GLOBAL_2D_ARRAY, 4)
      real clonu(GLOBAL_2D_ARRAY, 4)
      real clonv(GLOBAL_2D_ARRAY, 4)
      real clatr(GLOBAL_2D_ARRAY, 4)
      real clatu(GLOBAL_2D_ARRAY, 4)
      real clatv(GLOBAL_2D_ARRAY, 4)
      real arear(GLOBAL_2D_ARRAY)
      real areau(GLOBAL_2D_ARRAY)
      real areav(GLOBAL_2D_ARRAY)
      integer rmask_inv(GLOBAL_2D_ARRAY)
      integer umask_inv(GLOBAL_2D_ARRAY)
      integer vmask_inv(GLOBAL_2D_ARRAY)
      integer j,i, level
      CHARACTER*2 clevel
!
!
!     !-------------------------------------------------------------------
!     !             ***  ROUTINE cpl_prism_grid  ***
!     !
!     ! ** Purpose :   Define grids/areas/corners/masks information for ocean
!     !                exchange between CROCO and COUPLER. (OASIS3-MCT software)
!     !
!     ! ** Method  :   OASIS3 MPI communication
!     !--------------------------------------------------------------------
# ifdef AGRIF
      level=Agrif_Fixed()
      write(clevel,'(a,i1.1)') "_", level
# else
      level=0
      write(clevel,'(a)') "  "
# endif
!------------------------------------------------------------------
!     Set the CROCO grid definition to PRISM system
!------------------------------------------------------------------
!
!    Need to use oasis3_mct_v3 !
!
      CALL prism_start_grids_writing(error_flg)

! Model grid longitudes and latitudes
!------------------------------------------------------------------
!     -- t-grid --
      CALL prism_write_grid('ocnt'//trim(clevel),LLm,MMm,
     &           lonr(1:Lmmpi,1:Mmmpi),latr(1:Lmmpi,1:Mmmpi),
     &           id_partition)
!     -- u-grid --
      CALL prism_write_grid('ocnu'//trim(clevel),LLm,MMm,
     &           lonu(1:Lmmpi,1:Mmmpi),latu(1:Lmmpi,1:Mmmpi),
     &           id_partition)
!     -- v-grid --
      CALL prism_write_grid('ocnv'//trim(clevel),LLm,MMm,
     &           lonv(1:Lmmpi,1:Mmmpi),latv(1:Lmmpi,1:Mmmpi),
     &           id_partition)

! Model grid cell corner longitudes and latitudes
!------------------------------------------------------------------
      DO j=1, Mmmpi
        DO i=1, Lmmpi
!
!       -- t-grid --
!
          clonr(i,j,1)=0.5*(lonu(i+1,j)+lonu(i+1,j+1))
          clonr(i,j,2)=0.5*(lonu(i,j)+lonu(i,j+1))
          clonr(i,j,3)=0.5*(lonu(i,j-1)+lonu(i,j))
          clonr(i,j,4)=0.5*(lonu(i+1,j-1)+lonu(i+1,j))
          clatr(i,j,1)=0.5*(latu(i+1,j)+latu(i+1,j+1))
          clatr(i,j,2)=0.5*(latu(i,j)+latu(i,j+1))
          clatr(i,j,3)=0.5*(latu(i,j-1)+latu(i,j))
          clatr(i,j,4)=0.5*(latu(i+1,j-1)+latu(i+1,j))
!
!       -- u-grid --
!
!          clonu(i,j,1)=0.5*(lonr(i,j)+lonr(i,j+1))
!          clonu(i,j,2)=0.5*(lonr(i-1,j)+lonr(i-1,j+1))
!          clonu(i,j,3)=0.5*(lonr(i-1,j-1)+lonr(i-1,j))
!          clonu(i,j,4)=0.5*(lonr(i,j-1)+lonr(i,j))
!          clatu(i,j,1)=0.5*(latr(i,j)+latr(i,j+1))
!          clatu(i,j,2)=0.5*(latr(i-1,j)+latr(i-1,j+1))
!          clatu(i,j,3)=0.5*(latr(i-1,j-1)+latr(i-1,j))
!          clatu(i,j,4)=0.5*(latr(i,j-1)+latr(i,j))
!
          clonu(i,j,1)=lonv(i,j+1)
          clonu(i,j,2)=lonv(i-1,j+1)
          clonu(i,j,3)=lonv(i-1,j)
          clonu(i,j,4)=lonv(i,j)
          clatu(i,j,1)=latv(i,j+1)
          clatu(i,j,2)=latv(i-1,j+1)
          clatu(i,j,3)=latv(i-1,j)
          clatu(i,j,4)=latv(i,j)
!
!       -- v-grid --
!
!          clonv(i,j,1)=0.5*(lonr(i,j)+lonr(i+1,j))
!          clonv(i,j,2)=0.5*(lonr(i-1,j)+lonr(i,j))
!          clonv(i,j,3)=0.5*(lonr(i-1,j-1)+lonr(i,j-1))
!          clonv(i,j,4)=0.5*(lonr(i,j-1)+lonr(i+1,j-1))
!          clatv(i,j,1)=0.5*(latr(i,j)+latr(i+1,j))
!          clatv(i,j,2)=0.5*(latr(i-1,j)+latr(i,j))
!          clatv(i,j,3)=0.5*(latr(i-1,j-1)+latr(i,j-1))
!          clatv(i,j,4)=0.5*(latr(i,j-1)+latr(i+1,j-1))
!
          clonv(i,j,1)=lonu(i+1,j)
          clonv(i,j,2)=lonu(i,j)
          clonv(i,j,3)=lonu(i,j-1)
          clonv(i,j,4)=lonu(i+1,j-1)
          clatv(i,j,1)=latu(i+1,j)
          clatv(i,j,2)=latu(i,j)
          clatv(i,j,3)=latu(i,j-1)
          clatv(i,j,4)=latu(i+1,j-1)
!
        ENDDO
      ENDDO
!
      CALL prism_write_corner('ocnt'//trim(clevel), LLm, MMm, 4,
     &             clonr(1:Lmmpi,1:Mmmpi,:), clatr(1:Lmmpi,1:Mmmpi,:),
     &             id_partition)
!
      CALL prism_write_corner('ocnu'//trim(clevel), LLm, MMm, 4,
     &             clonu(1:Lmmpi,1:Mmmpi,:), clatu(1:Lmmpi,1:Mmmpi,:),
     &             id_partition)
!
      CALL prism_write_corner('ocnv'//trim(clevel), LLm, MMm, 4,
     &             clonv(1:Lmmpi,1:Mmmpi,:), clatv(1:Lmmpi,1:Mmmpi,:),
     &             id_partition)
!
! Model grid cell area
!------------------------------------------------------------------
!
      DO j=1, Mmmpi
        DO i=1, Lmmpi
          arear(i,j)=om_r(i,j)*on_r(i,j)
          areau(i,j)=om_u(i,j)*on_u(i,j)
          areav(i,j)=om_v(i,j)*on_v(i,j)
        ENDDO
      ENDDO
!
      CALL prism_write_area('ocnt'//trim(clevel), LLm, MMm,
     &                         arear(1:Lmmpi,1:Mmmpi),
     &                         id_partition)

      CALL prism_write_area('ocnu'//trim(clevel), LLm, MMm,
     &                         areau(1:Lmmpi,1:Mmmpi),
     &                         id_partition)

      CALL prism_write_area('ocnv'//trim(clevel), LLm, MMm,
     &                         areav(1:Lmmpi,1:Mmmpi),
     &                         id_partition)

! Model grid mask : 0=sea / 1=land
!------------------------------------------------------------------
      rmask_inv=0
      umask_inv=0
      vmask_inv=0
      DO j=1, Mmmpi
        DO i=1, Lmmpi
          IF (rmask(i,j)==0) rmask_inv(i,j)=1
          IF (umask(i,j)==0) umask_inv(i,j)=1
          IF (vmask(i,j)==0) vmask_inv(i,j)=1
        ENDDO
      ENDDO

!     -- t-grid --
      CALL prism_write_mask('ocnt'//trim(clevel), LLm, MMm,
     &                     rmask_inv(1:Lmmpi,1:Mmmpi),
     &                     id_partition)

!     -- u-grid --
      CALL prism_write_mask('ocnu'//trim(clevel), LLm, MMm,
     &                     umask_inv(1:Lmmpi,1:Mmmpi),
     &                     id_partition)

!     -- v-grid --
      CALL prism_write_mask('ocnv'//trim(clevel), LLm, MMm,
     &                     vmask_inv(1:Lmmpi,1:Mmmpi),
     &                     id_partition)
      
! Terminate grid writing
!------------------------------------------------------------------
      CALL prism_terminate_grids_writing()
      CALL MPI_Barrier(ocean_grid_comm, error_flg)
!------------------------------------------------------------------
!     End of grid definition
!------------------------------------------------------------------
      return
      end
#else 
!---------------------------------------------------------------------------------------
!     Dummy subroutine                                                           NO Coupling
!---------------------------------------------------------------------------------------
      subroutine cpl_prism_grid_empty
      end
#endif 
