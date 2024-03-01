#include "cppdefs.h"
#ifdef ONLINE_ANALYSIS

      module module_grd_oa

!======================================================================
!
!> @brief Online Analysis (OA=OnlineAnalysis)
!! - Temporary use of the local derived type and pointers 
!!   association with the Croco grid arrays in order to  
!!   build the spatial structure of the OA state vector.
!
!! @details required when applying croco tile(-threads)
!!  with/without dual OpenMP-MPI parallelization of the horizontal domain.
!!  Some OA variable/parameters are tile(-threads) dependent,
!!  especially all variables depending on the ocean domain.
!!  More comments/info in source module_interface_oa.F90
!
!> @authors  
!! - B. Lemieux-Dudon
!!  - Croco Tile-thread compliant version (2021). 
!!  - derived data type for the grid array
!!  - fields of the local derived type have inherited from older 
!!    pieces of code of older OA versions
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version (2006)
!> @todo BLXD double precision type is obsolete.
!!  Module to define precision consistently with Croro grid arrays ?
!
!  REFERENCE:
!  
!======================================================================

      implicit none

      type grid_str_oa

          ! 2-Dimensional pointer to handle coordinate arrays
          double precision, dimension(:,:), pointer :: lon_t => null()
          double precision, dimension(:,:), pointer :: lon_u => null()
          double precision, dimension(:,:), pointer :: lon_v => null()
          double precision, dimension(:,:), pointer :: lon_f => null()
          double precision, dimension(:,:), pointer :: lat_t => null()
          double precision, dimension(:,:), pointer :: lat_u => null()
          double precision, dimension(:,:), pointer :: lat_v => null()
          double precision, dimension(:,:), pointer :: lat_f => null()
          double precision, dimension(:,:), pointer :: h_w   => null()
          double precision, dimension(:,:), pointer :: h_u   => null()
          double precision, dimension(:,:), pointer :: h_v   => null()
          double precision, dimension(:,:), pointer :: h_f   => null()

          !> 3-Dimensional pointer to handle mask
          integer, dimension(:,:,:), pointer :: mask_t => null()
          integer, dimension(:,:,:), pointer :: mask_u => null()
          integer, dimension(:,:,:), pointer :: mask_v => null()
          integer, dimension(:,:,:), pointer :: mask_f => null()

          integer :: imin, imax, jmin, jmax, kmin, kmax 

      end type

      type(grid_str_oa), allocatable, dimension(:), target  :: grd
      !BLXD_TILE_ISSUE rm pg from list of public module variables add grid_str_oa
      !type(grid_str_oa), pointer :: pg => null() ! croco tile-thread SHARED mod. var.

      private
      public :: associate_grid_oa_ptr, nullify_grid_oa_ptr, grd, grid_str_oa &
                ,allocate_tile_grid_tmp_str, deallocate_tile_grid_tmp_str    &
                ,get_3D_grid_mask, get_grid_mask, get_2D_subdomain_minmax

      CONTAINS 

      subroutine get_3D_grid_mask( i,j,k, grd_pt_code, pg &
     &                            ,h_g, lat_g, lon_g, msk_g )

        implicit none

        type(grid_str_oa), intent(in) :: pg
        integer, intent(in)           :: i,j,k 
        integer, intent(in)           :: grd_pt_code
        double precision, intent(out) :: h_g, lat_g, lon_g
        integer, intent(out)          :: msk_g 

        if (k/=pg%kmax) stop
        ! BLXD croco terrain following sigma coord.
        ! 2D rmask(i,j) can be used for w-grid point
        ! but avoid k=0
        if (grd_pt_code.eq.1) then
           h_g   = pg%h_w     (i,j)
           lat_g = pg%lat_t  (i,j)
           lon_g = pg%lon_t  (i,j)
           msk_g = pg%mask_t(i,j,pg%kmax) 

        elseif (grd_pt_code.eq.2) then
           h_g   = pg%h_u     (i,j)
           lat_g = pg%lat_u  (i,j)
           lon_g = pg%lon_u  (i,j)
           msk_g = pg%mask_u(i,j,pg%kmax) 
        elseif (grd_pt_code.eq.3) then
           h_g   = pg%h_v     (i,j)
           lat_g = pg%lat_v  (i,j)
           lon_g = pg%lon_v  (i,j)
           msk_g = pg%mask_v(i,j,pg%kmax) 
        else
           h_g   = pg%h_f     (i,j)
           lat_g = pg%lat_f  (i,j)
           lon_g = pg%lon_f  (i,j)
           msk_g = pg%mask_f(i,j,pg%kmax) 
        endif

        return
      end subroutine get_3D_grid_mask

      subroutine get_grid_mask( i,j,k, grd_pt_code, pg, msk_g )

        implicit none

        type(grid_str_oa), intent(in) :: pg
        integer, intent(in)           :: i,j,k 
        integer, intent(in)           :: grd_pt_code
        integer, intent(out)          :: msk_g 

        if (grd_pt_code.eq.1) then
           msk_g = pg%mask_t(i,j,pg%kmax) 
        elseif (grd_pt_code.eq.2) then
           msk_g = pg%mask_u(i,j,pg%kmax) 
        elseif (grd_pt_code.eq.3) then
           msk_g = pg%mask_v(i,j,pg%kmax) 
        else
           msk_g = pg%mask_f(i,j,pg%kmax) 
        endif

        return
      end subroutine get_grid_mask


      ! BLXD could be modified to meet up fortran 2008 standard
      subroutine allocate_tile_grid_tmp_str( tile_size )

        implicit none
        integer, intent(in) :: tile_size

        if ( .not. allocated(grd) ) allocate(grd(1:tile_size))

        return
      end subroutine allocate_tile_grid_tmp_str

      subroutine deallocate_tile_grid_tmp_str( tile_size )

        implicit none
        integer, intent(in), optional :: tile_size
        integer :: itil, itilsize

        if ( allocated(grd) ) then

            if ( present( tile_size ) ) then
                itilsize = size( grd(:) )
                if ( tile_size /= itilsize ) stop
            endif

            !do itil=1,tile_size
            !    call nullify_grid_oa_ptr(itil) 
            !enddo
            
            deallocate(grd)

        endif

        return
      end subroutine deallocate_tile_grid_tmp_str

                                                                                         
      subroutine associate_grid_oa_ptr( tile     &
                  ,imin, imax                    &
                  ,jmin, jmax                    &
                  ,kmin, kmax                    &
                  ,lon_t_oa                      &
                  ,lat_t_oa                      &
                  ,lon_u_oa                      &
                  ,lat_u_oa                      &
                  ,lon_v_oa                      &
                  ,lat_v_oa                      &
                  ,lon_f_oa                      &
                  ,lat_f_oa                      &
                  ,mask_t_oa                     &
                  ,mask_f_oa                     &
                  ,mask_u_oa                     &
                  ,mask_v_oa                     &
                  ,h_w_oa                        &
                  ,h_u_oa                        &
                  ,h_v_oa                        &
                  ,h_f_oa                        ) 

      implicit none

      !> Grid index range for the analysis of fields passed in argument to the OA module 
      integer, intent(in) :: tile &
           ,imin, imax            & 
           ,jmin, jmax            &
           ,kmin, kmax

      ! 2-Dimensional horizontal and vertical grid
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: lon_t_oa  !< Longitude at u-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: lon_u_oa  !< Longitude at u-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: lon_v_oa  !< Longitude at v-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: lon_f_oa  !< Longitude at f-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: lat_t_oa  !< Latitude at t-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: lat_u_oa  !< Latitude at u-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: lat_v_oa  !< Latitude at v-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: lat_f_oa  !< Latitude at f-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: h_w_oa    !< Ocean thickness at w-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: h_u_oa    !< Ocean thickness at u-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: h_v_oa    !< Ocean thickness at v-grid point
      double precision, dimension(imin:imax,jmin:jmax), intent(in), target :: h_f_oa    !< Ocean thickness at f-grid point           


      !> 3-Dimensional mask arrays passed from the calling code
      integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in), target :: mask_t_oa
      integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in), target :: mask_u_oa
      integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in), target :: mask_v_oa
      integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in), target :: mask_f_oa

      grd(tile)%lon_t(imin:imax,jmin:jmax) => lon_t_oa(:,:)
      grd(tile)%lon_u(imin:imax,jmin:jmax) => lon_u_oa(:,:)
      grd(tile)%lon_v(imin:imax,jmin:jmax) => lon_v_oa(:,:)
      grd(tile)%lon_f(imin:imax,jmin:jmax) => lon_f_oa(:,:)
      grd(tile)%lat_t(imin:imax,jmin:jmax) => lat_t_oa(:,:)
      grd(tile)%lat_u(imin:imax,jmin:jmax) => lat_u_oa(:,:)
      grd(tile)%lat_v(imin:imax,jmin:jmax) => lat_v_oa(:,:)
      grd(tile)%lat_f(imin:imax,jmin:jmax) => lat_f_oa(:,:)
      grd(tile)%h_w  (imin:imax,jmin:jmax) => h_w_oa  (:,:)
      grd(tile)%h_u  (imin:imax,jmin:jmax) => h_u_oa  (:,:)
      grd(tile)%h_v  (imin:imax,jmin:jmax) => h_v_oa  (:,:)
      grd(tile)%h_f  (imin:imax,jmin:jmax) => h_f_oa  (:,:)

      grd(tile)%mask_t(imin:imax,jmin:jmax,kmin:kmax) => mask_t_oa(:,:,:)
      grd(tile)%mask_u(imin:imax,jmin:jmax,kmin:kmax) => mask_u_oa(:,:,:)
      grd(tile)%mask_v(imin:imax,jmin:jmax,kmin:kmax) => mask_v_oa(:,:,:)
      grd(tile)%mask_f(imin:imax,jmin:jmax,kmin:kmax) => mask_f_oa(:,:,:)

      grd(tile)%imin = imin
      grd(tile)%imax = imax
      grd(tile)%jmin = jmin
      grd(tile)%jmax = jmax
      grd(tile)%kmin = kmin 
      grd(tile)%kmax = kmax

      end subroutine associate_grid_oa_ptr

      subroutine nullify_grid_oa_ptr( tile )

      implicit none 
      integer, intent(in) :: tile

      grd(tile)%lon_t => null()
      grd(tile)%lon_u => null()
      grd(tile)%lon_v => null()
      grd(tile)%lon_f => null()
      grd(tile)%lat_t => null()
      grd(tile)%lat_u => null()
      grd(tile)%lat_v => null()
      grd(tile)%lat_f => null()
      grd(tile)%h_w   => null()
      grd(tile)%h_u   => null()
      grd(tile)%h_v   => null()
      grd(tile)%h_f   => null()

      grd(tile)%mask_t => null()
      grd(tile)%mask_u => null()
      grd(tile)%mask_v => null()
      grd(tile)%mask_f => null()

      !print*,'OA out : associate '

      return
      end subroutine nullify_grid_oa_ptr

      subroutine get_2D_subdomain_minmax( imin, imax, jmin, jmax        &
     &                                   ,grd_pt_code, pg               &
     &                                   ,latmin, lonmin, latmax, lonmax )

        implicit none

        type(grid_str_oa), intent(in) :: pg
        integer, intent(in)           :: imin, imax, jmin, jmax
        integer, intent(in)           :: grd_pt_code
        double precision, intent(out) :: latmin, lonmin, latmax, lonmax 

        if (grd_pt_code.eq.1) then
           latmin = pg%lat_t  (imin,jmin)
           lonmin = pg%lon_t  (imin,jmin)
           latmax = pg%lat_t  (imax,jmax)
           lonmax = pg%lon_t  (imax,jmax)
        elseif (grd_pt_code.eq.2) then
           latmin = pg%lat_u  (imin,jmin)
           lonmin = pg%lon_u  (imin,jmin)
           latmax = pg%lat_u  (imax,jmax)
           lonmax = pg%lon_u  (imax,jmax)
        elseif (grd_pt_code.eq.3) then
           latmin = pg%lat_v  (imin,jmin)
           lonmin = pg%lon_v  (imin,jmin)
           latmax = pg%lat_v  (imax,jmax)
           lonmax = pg%lon_v  (imax,jmax)
        else
        !elseif (mod(tgv_oa(tv_oa(iv_g)),5).eq.3) then
           latmin = pg%lat_f  (imin,jmin)
           lonmin = pg%lon_f  (imin,jmin)
           latmax = pg%lat_f  (imax,jmax)
           lonmax = pg%lon_f  (imax,jmax)
        endif

        return
      end subroutine get_2D_subdomain_minmax

end module module_grd_oa

#else /* ONLINE_ANALYSIS */
      module module_grd_oa_empty
      end module
#endif /* ONLINE_ANALYSIS */

