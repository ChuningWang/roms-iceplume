#include "cppdefs.h"
#ifdef ONLINE_ANALYSIS

      module module_parameter_oa

!======================================================================
!
!> @brief Croco interface for Online Analysis (OA=OnlineAnalysis)
!! - Allocation of arrays to store analyses
!
!! @details 
!!  More history/info in source module_interface_oa.F90
!
!> @authors  
!! - B. Lemieux-Dudon
!!  - Toward Croco Tile(-threads) compliant version (2021) 
!!  - Based on a preliminary Croco-OA interface version : spring 2020.
!! - More history (authors, comments) in source module_interface_oa.F90
!> @todo BLXD Scalogram arrays could be initialized elsewhere since they 
!!   do not depend on Croco spatial MPI subdomain dimensions 
!!   and/or tile dimensions.
!
!  REFERENCE:
!  
!======================================================================
        
      use scalars

      implicit none

      contains

      subroutine allocate_3D_glob_array_oa_cplx( arr_oa, nzupd3d_oa, kst )
      implicit none

      complex(8), dimension(:,:,:,:), allocatable, intent(out) :: arr_oa 
      integer, intent(in), optional :: kst
      integer, intent(in)           :: nzupd3d_oa
      integer                       :: kmin

        if ( present(kst) ) then
            kmin=kst
        else
            kmin=1
        endif

        if ( .not. allocated(arr_oa) ) then
            allocate( arr_oa( GLOBAL_2D_ARRAY, kmin:N, nzupd3d_oa ) )
            arr_oa(:,:,:,:) = ( 0.D0, 0.D0 )
        endif
        return

      end subroutine allocate_3D_glob_array_oa_cplx

      subroutine allocate_2D_glob_array_oa_cplx( arr_oa, nzupd2d_oa )
      implicit none
      complex(8), dimension(:,:,:), allocatable, intent(out) :: arr_oa
      integer, intent(in) :: nzupd2d_oa
       

        if ( .not. allocated(arr_oa) ) then
            allocate( arr_oa( GLOBAL_2D_ARRAY, nzupd2d_oa) )
            arr_oa(:,:,:) = ( 0.D0, 0.D0 )
        endif
        return

      end subroutine allocate_2D_glob_array_oa_cplx

!BLXD TODO scalogram arrays could be initialized elsewhere since they do not depend
!     on Croco spatial MPI subdomain dimensions and/or tile dimensions
     subroutine allocate_0D_sclg_array_oa_cplx( arr_oa, nper_sclg_oa, nzupd0d_oa )
     implicit none
     complex(8), dimension(:,:), allocatable, intent(out) :: arr_oa
     integer, intent(in) :: nzupd0d_oa, nper_sclg_oa
      
     
       if ( .not. allocated(arr_oa) ) then
           allocate( arr_oa( nper_sclg_oa, nzupd0d_oa) )
           arr_oa(:,:) = ( 0.D0, 0.D0 )
       endif
       return
     
     end subroutine allocate_0D_sclg_array_oa_cplx

     subroutine allocate_0D_sclg_array_oa_real( arr_oa, nper_sclg_oa, nzupd0d_oa )
     implicit none
     real(8), dimension(:,:), allocatable, intent(out) :: arr_oa
     integer, intent(in) :: nzupd0d_oa, nper_sclg_oa
      
       if ( .not. allocated(arr_oa) ) then
           allocate( arr_oa(  nper_sclg_oa, nzupd0d_oa) )
       endif
       return
     
     end subroutine allocate_0D_sclg_array_oa_real

     subroutine allocate_0D_sclg_array_oa_int( arr_oa, nzupd0d_oa )
     implicit none
     integer, dimension(:), allocatable, intent(out) :: arr_oa
     integer, intent(in) :: nzupd0d_oa
      
     
       if ( .not. allocated(arr_oa) ) then
           allocate( arr_oa( nzupd0d_oa) )
       endif
       return
     
     end subroutine allocate_0D_sclg_array_oa_int

end module module_parameter_oa

#else /* ONLINE_ANALYSIS */
      module module_parameter_oa_empty
      end module
#endif /* ONLINE_ANALYSIS */
