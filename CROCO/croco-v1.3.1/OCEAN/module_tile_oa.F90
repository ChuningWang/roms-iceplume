#include "cppdefs.h"
#ifdef ONLINE_ANALYSIS

      module module_tile_oa

!======================================================================
!
!> @brief Online Analysis (OA=OnlineAnalysis)
!! - Derived data type to handle the vectorization of
!!   the spatial structure of the OA state vector
!!   compliant with Croco Tile(-threads) 
!
!! @details required when applying croco tile(-threads)
!!  with/without dual OpenMP-MPI parallelization of the horizontal domain.
!!  Some OA variable/parameters are tile(-threads) dependent,
!!  especially all variables depending on the ocean domain.
!!  More comments/info in source module_interface_oa.F90
!
!> @authors  
!! - B. Lemieux-Dudon
!!  - Toward a Croco Tile(-threads) compliant version (2021) 
!!  - More history (authors, comments) in source module_interface_oa.F90
!!  - Some fields of the local derived type are inherited from older 
!!    pieces of code of older OA version (F. Auclair and team, 2006).
!> @todo BLXD 
!!  - test with a Croco tile(-threads) working version 
!!  - quick patch including state vectors within the tile_space_str
!  REFERENCE:
!  
!======================================================================

      use module_oa_type

      implicit none

      integer :: ntiles  ! croco tile-thread THREADPRIVATE option/copying from master value bef. thread fork ?

      type tile_sclg_str
        integer, dimension(:), allocatable     :: v2sclg
        integer                   :: nsclg
        integer                   :: nper
      end type tile_sclg_str

      type tile_test_str
        integer, dimension(:,:,:), allocatable     :: vardp_test
        integer                            :: imin,imax,jmin,jmax
        integer                            :: kmin,kmax
      end type tile_test_str

      type tile_space_str

        integer                            :: imin,imax,jmin,jmax
        integer                            :: kmin,kmax  ! vertical dim added even though not tile dependent
        integer, dimension(:), allocatable :: begvs_oa   ! nzv_oa 
        integer, dimension(:), allocatable :: begvs3d_oa, kmin3d_oa ! nzvs_oa + 1
        integer, dimension(:), allocatable :: l2i_oa, l2j_oa ! nzvs_oa
        integer, dimension(:,:,:), allocatable :: ij2l_oa    ! imin:imax,jmin:jmax, nzv_oa

        ! BLDX removed from module_oa_variables since dependent on tile(-threads)
        !> Dimensions of the 2D and 3D spatial structures - dependent on tile-(threads) 
        integer :: nzvs3d_oa, nzvs_oa

        type(type_wf),dimension(:),allocatable::wf_oa

        type(type_level),dimension(:),allocatable::wlev_oa

        integer :: nzw_oa, nzlevel

      end type tile_space_str

      type(tile_sclg_str), allocatable, dimension(:), target :: scl

      type(tile_test_str), allocatable, dimension(:), target :: sts

      type(tile_space_str), allocatable, dimension(:), target :: st

      private

      public :: st, tile_space_str, ntiles, allocate_tile_space_str, deallocate_tile_space_str, &
                allocate_begvs_oa, allocate_begvs3d_oa, allocate_wf_oa, allocate_wlev_oa,       &
                deallocate_tile_varspace_oa,                                                    &
                scl, tile_sclg_str, allocate_tile_sclg_str, deallocate_tile_sclg_str,           &
                allocate_tile_sclg_oa, deallocate_tile_sclg_oa,                                 &
                sts, tile_test_str, allocate_tile_test_str, deallocate_tile_test_str,           &
                allocate_tile_test_oa, deallocate_tile_test_oa


      CONTAINS 

      subroutine allocate_tile_space_str( tile_num )

        implicit none

        integer, intent(in) :: tile_num

        if ( .not. allocated(st) ) then
            allocate(st(1:tile_num))
            ntiles = tile_num
        end if

        return
      end subroutine allocate_tile_space_str
                                                                                        
      subroutine allocate_tile_sclg_str( tile_num )
        implicit none
        integer, intent(in) :: tile_num

        !CR print*,'IN allocate_tile_sclg_str ', ntiles
        if ( .not. allocated(scl) ) then
            allocate(scl(1:tile_num))
            ntiles = tile_num
        end if
        !CR print*,'OUT allocate_tile_sclg_str ', ntiles
        return
      end subroutine allocate_tile_sclg_str

      subroutine allocate_tile_test_str( tile_num )
        implicit none
        integer, intent(in) :: tile_num

        if ( .not. allocated(sts) ) then
            allocate(sts(1:tile_num))
            ntiles = tile_num
        end if
        
        return
      end subroutine allocate_tile_test_str

      subroutine allocate_wf_oa( tile, nmsimult_oa )
        
        use module_oa_time, only : tallocated_oa
        implicit none
        
        integer, intent(in) :: tile
        integer, intent(in) :: nmsimult_oa 
        integer :: l_a

        if ( .not. allocated(st(tile)%wf_oa) ) then
            allocate ( st(tile)%wf_oa(nmsimult_oa) )
            do l_a=1,nmsimult_oa
               nullify(st(tile)%wf_oa(l_a)%coef)
            enddo
            ! see unset_convwind_oa in module_oa_interface
            !tallocated_oa(:)   = -1
            st(tile)%wf_oa(:)%t_indice  = -1
            st(tile)%wf_oa(:)%config    = -1
            st(tile)%wf_oa(:)%variable  = -1
            st(tile)%nzw_oa = nmsimult_oa
        endif

      end subroutine allocate_wf_oa

      subroutine allocate_wlev_oa( tile, nzlevel_oa )

        implicit none
        
        integer, intent(in) :: tile
        integer, intent(in) :: nzlevel_oa 
        integer :: l_a

        if ( .not. allocated(st(tile)%wlev_oa) ) then
            allocate ( st(tile)%wlev_oa(nzlevel_oa) )
            do l_a=1,nzlevel_oa
               nullify(st(tile)%wlev_oa(l_a)%rhp)
               nullify(st(tile)%wlev_oa(l_a)%z)
               nullify(st(tile)%wlev_oa(l_a)%k)
            enddo
            st(tile)%nzlevel = nzlevel_oa
        endif

      end subroutine allocate_wlev_oa


      subroutine allocate_begvs_oa( imin, imax, jmin, jmax, kmin, kmax, tile )

        use module_oa_variables, only : nzv_oa

        implicit none
        
        integer, intent(in) :: tile
        integer, intent(in) :: imin, imax, jmin, jmax, kmin, kmax

        st(tile)%imin = imin
        st(tile)%imax = imax
        st(tile)%jmin = jmin
        st(tile)%jmax = jmax
        st(tile)%kmin = kmin
        st(tile)%kmax = kmax

        if ( .not. allocated( st(tile)%begvs_oa ) ) then
            allocate( st(tile)%begvs_oa( nzv_oa + 1 ) )
        end if 

        if ( .not. allocated( st(tile)%ij2l_oa ) ) then
            allocate( st(tile)%ij2l_oa( imin:imax, jmin:jmax, nzv_oa ) )
        end if 

        return 
      end subroutine allocate_begvs_oa

      subroutine allocate_begvs3d_oa(tile)
        implicit none
        integer, intent(in) :: tile
        integer :: nzvs_oa

        nzvs_oa = st(tile)%nzvs_oa

        if ( .not. allocated( st(tile)%begvs3d_oa ) ) then
            allocate( st(tile)%begvs3d_oa( nzvs_oa + 1 ) )
        end if 

        ! BLXD check why kmin3d_oa as size nzvs_oa PLUS ONE
        if ( .not. allocated( st(tile)%kmin3d_oa ) ) then
            allocate( st(tile)%kmin3d_oa( nzvs_oa + 1 ) )
        end if 

        if ( .not. allocated( st(tile)%l2i_oa ) ) then
            allocate( st(tile)%l2i_oa( nzvs_oa ) )
        end if 

        if ( .not. allocated( st(tile)%l2j_oa ) ) then
            allocate( st(tile)%l2j_oa( nzvs_oa ) )
        end if 

        return
      end subroutine allocate_begvs3d_oa

      subroutine deallocate_tile_varspace_oa(tile)
        implicit none
        integer, intent(in) :: tile
        integer :: l_a
 
        if ( allocated( st(tile)%begvs_oa ) ) then
            deallocate( st(tile)%begvs_oa )
        end if 
        if ( allocated( st(tile)%begvs3d_oa ) ) then
            deallocate( st(tile)%begvs3d_oa )
        end if 
        if ( allocated( st(tile)%kmin3d_oa ) ) then
            deallocate( st(tile)%kmin3d_oa )
        end if 
        if ( allocated( st(tile)%ij2l_oa ) ) then
            deallocate( st(tile)%ij2l_oa )
        end if 
        if ( allocated( st(tile)%l2i_oa ) ) then
            deallocate( st(tile)%l2i_oa )
        end if 
        if ( allocated( st(tile)%l2j_oa ) ) then
            deallocate( st(tile)%l2j_oa )
        end if 

        if ( allocated(st(tile)%wf_oa) ) then
            do l_a=1,st(tile)%nzw_oa
               nullify(st(tile)%wf_oa(l_a)%coef)
            enddo
            deallocate ( st(tile)%wf_oa)
        endif

        if ( allocated(st(tile)%wlev_oa) ) then
            do l_a=1,st(tile)%nzlevel
               nullify(st(tile)%wlev_oa(l_a)%rhp)
               nullify(st(tile)%wlev_oa(l_a)%z)
               nullify(st(tile)%wlev_oa(l_a)%k)
            enddo
            deallocate ( st(tile)%wlev_oa)
        endif

        return
      end subroutine deallocate_tile_varspace_oa

      subroutine deallocate_tile_space_str( )
      implicit none
        ! BLXD removed since moved in rout. online_spectral_diags
        ! if ( allocated( st(tile)%begvs_oa ) ) then
        !    stop
        ! else
        if ( allocated( st ) ) then
            deallocate(st)
        end if
        return
      end subroutine deallocate_tile_space_str

      subroutine allocate_tile_sclg_oa( tile, nzv_oa, nper_sclg)
        implicit none
        integer, intent(in) :: tile, nzv_oa, nper_sclg

        ! BLXD this routine can eventually be called early in init_parameter_oa
        if ( .not. allocated( scl(tile)%v2sclg ) ) then
            allocate( scl(tile)%v2sclg( nzv_oa ) )
        end if 
        ! one scalogram per variable iv : -99 => not a iv config-var with scalogram
        scl(tile)%v2sclg( 1:nzv_oa ) = -99
        scl(tile)%nper    = nper_sclg
        scl(tile)%nsclg   = 0
        
        return
      end subroutine allocate_tile_sclg_oa

      subroutine set_tile_sclg_oa( tile, nper_sclg, nzupd0d )
        implicit none
        integer, intent(in) :: tile, nper_sclg, nzupd0d
        ! BLXD not useful anymore
        !if ( .not. allocated( scl(tile)%scal0d ) ) then
        !    allocate( scl(tile)%scal0d( nper_sclg, nzupd0d ) )
        !end if 
        scl(tile)%nper    = nper_sclg
        scl(tile)%nsclg   = nzupd0d
        return
      end subroutine set_tile_sclg_oa

      subroutine deallocate_tile_sclg_oa(tile)
        implicit none
        integer, intent(in) :: tile
 
        if ( allocated( scl(tile)%v2sclg ) ) then
            deallocate( scl(tile)%v2sclg )
        end if 

        !if ( allocated( scl(tile)%scal0d ) ) then
        !    deallocate( scl(tile)%scal0d )
        !end if 

        return
      end subroutine deallocate_tile_sclg_oa

      subroutine deallocate_tile_sclg_str
      implicit none
        if ( allocated( scl ) ) then
            deallocate(scl)
        end if
        return
      end subroutine deallocate_tile_sclg_str

      subroutine allocate_tile_test_oa( tile, imin, imax, jmin, jmax, kmin, kmax )
        implicit none
        integer, intent(in) :: tile, imin, imax, jmin, jmax, kmin, kmax

        if ( .not. allocated( sts(tile)%vardp_test ) ) then
            allocate( sts(tile)%vardp_test(imin:imax,jmin:jmax,kmin:kmax) )
            sts(tile)%imin = imin ; sts(tile)%imax = imax
            sts(tile)%jmin = jmin ; sts(tile)%jmax = jmax
            sts(tile)%kmin = kmin ; sts(tile)%kmax = kmax
        end if 
        return
      end subroutine allocate_tile_test_oa

      subroutine deallocate_tile_test_oa(tile)
        implicit none
        integer, intent(in) :: tile
 
        if ( allocated( sts(tile)%vardp_test ) ) then
            deallocate( sts(tile)%vardp_test )
        end if 
      end subroutine deallocate_tile_test_oa

      subroutine deallocate_tile_test_str
      implicit none
        if ( allocated( sts ) ) then
            deallocate(sts)
        end if
        return
      end subroutine deallocate_tile_test_str


end module module_tile_oa

#else /* ONLINE_ANALYSIS */
      module module_tile_oa_empty
      end module
#endif /* ONLINE_ANALYSIS */

