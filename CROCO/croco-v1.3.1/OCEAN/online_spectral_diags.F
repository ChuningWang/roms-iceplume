! $Id: online_spectral_diags.F 1458 2014-02-03 15:01:25Z blxd $
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
!======================================================================
!
!> @brief Croco interface for Online Analysis (OA=OnlineAnalysis)
!! - Performs online spectral and wavelet analyses.
!
!! @details written to support tile-threads and 
!!  possible dual OpenMP-MPI parallelization. 
!!  More history/info in source module_interface_oa.F90
!
!> @authors  
!! - B. Lemieux-Dudon
!!  - Twd a Croco Tile(-threads) compliant version (2021) 
!!    supporting tile-threads and possible dual OpenMP-MPI 
!!    parallelization of the horizontal domain
!!     => online_spectral_diags : 
!!        main routine called by Croco main.F and step.F 
!!     => module_tile_oa : tile-compliant OA state vector 
!!        spatial structure
!!     => output_oa      : OA updates global array looping on tiles 
!!        and finally send them to XIOS (out of the loop)
!!  - based on a preliminary Croco-OA interface version : spring 2020.
!! - More history (authors, comments) in source module_interface_oa.F90
!> @todo BLXD 
!!  - test with a Croco Tile(-threads) working version
!!  - isopycnal level tracking not yet tested with Croco-OA (2021)
!  REFERENCE:
!  
!======================================================================

#include "cppdefs.h"
! BLXD includes :
! - cppdefs_dev.h
! - set_global_definitions.h
!   with GLOBAL_2D_ARRAY
!        PRIVATE_SCRATCH_2D_ARRAY... 
#ifdef ONLINE_ANALYSIS

      subroutine online_spectral_diags (tile, icall)
      implicit none
      integer, intent(in) :: tile
      integer, intent(in) :: icall
C$    integer trd, omp_get_thread_num

# include "param.h"
# include "compute_tile_bounds.h"
! compute_tile_bounds.h : 
! 1) declaration
! 1) instructions for computing
!    Istr,Iend,Jstr,Jend
!    Using param.h, MPI subdom. (or not) and tile var. 
C$      trd=omp_get_thread_num()

      if ( icall==0 ) then

        call init_online_spectral_diags_tile (Istr,Iend
     &                                       ,Jstr,Jend
     &                                       ,tile) 
! No need of any openMP BARRIER since included at the end 
! of init_online_spectral_diags_tile

      else if ( icall==1 ) then

      !endif
      !if ( icall==0 .or. icall==1 ) then
        call run_online_spectral_diags_tile (Istr,Iend,Jstr,Jend,tile)
      end if
      return
      end subroutine online_spectral_diags
#else /* ONLINE_ANALYSIS */
      subroutine online_spectral_diags_empty
      return
      end subroutine online_spectral_diags_empty
#endif /* ONLINE_ANALYSIS */

!
!----------------------------------------------------------------------
!
!  Initialization of the online spectral analysis module
!
!----------------------------------------------------------------------
!
#ifdef ONLINE_ANALYSIS

      subroutine init_online_spectral_diags_tile (Istr,Iend,Jstr,Jend
     &                                           ,tile)
!
!
!#ifdef NBQ
!      use module_nh 
!      use module_nbq
!#endif
      use module_interface_oa, only : init_oa
      use module_grd_oa, only : deallocate_tile_grid_tmp_str
! BLXD TODO (eventually) TEST 
! use scalars => i) use param ii) include "scalars.h" 
! use param   => include "param.h"
! REQUIRES to comment
! - include "param.h"
! - include "scalars.h"
!      use scalars
      implicit none
      integer, intent(in) :: Istr, Iend, Jstr, Jend
      integer, intent(in) :: tile

! >> Include with variable declaration
# include "param.h"
# include "scalars.h"
! BLXD init_oa doesn't require (yet) 2D ocean fields 
! # include "ocean2d.h"
# include "ocean3d.h"
# include "grid.h"
#ifdef NBQ
# include "nbq.h"
#endif

! >> Put declaration bef. computations
! Local variables
      integer :: i,j,k    

      integer ::                        
     &  maskr_c (PRIVATE_2D_SCRATCH_ARRAY,1:N)  
     & ,masku_c (PRIVATE_2D_SCRATCH_ARRAY,1:N)  
     & ,maskv_c (PRIVATE_2D_SCRATCH_ARRAY,1:N)  
     & ,maskf_c (PRIVATE_2D_SCRATCH_ARRAY,1:N)  

      real ::                           
     &  hu_c(PRIVATE_2D_SCRATCH_ARRAY)           
     & ,hv_c(PRIVATE_2D_SCRATCH_ARRAY)          

#ifndef SPHERICAL
      real ::                          
     &   xu_c (PRIVATE_2D_SCRATCH_ARRAY)
     &  ,yu_c (PRIVATE_2D_SCRATCH_ARRAY)
     &  ,xv_c (PRIVATE_2D_SCRATCH_ARRAY)
     &  ,yv_c (PRIVATE_2D_SCRATCH_ARRAY)
#endif

      integer, dimension(1:2) :: ip_ran,jp_ran
     &                          ,ir_ran,jr_ran
     &                          ,iu_ran,ju_ran
     &                          ,iv_ran,jv_ran


! >> Include with computation instructions
# include "compute_auxiliary_bounds.h"
! From param.h, MPI subdom., tile 
! From compute_tile_bounds.h with Istr,Iend,Jstr,Jend
! Computes Physical Bdy indices IstrR, IendR, JstrR, JendR
!          Interior U or V points (start) IstrU, JstrV
!
! Does not require #include def_bounds.h
! LOCALLM, LOCALMM not needed at tile level

#  define I_EXT_RANGE Istr-1,Iend+1
#  define J_EXT_RANGE Jstr-1,Jend+1

#  define I_PRI_RANGE Istr-2,Iend+2
#  define J_PRI_RANGE Jstr-2,Jend+2

# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif

# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif

! BLXD OA diags must include Bdy
! - rho IstrR:IendR, JstrR:JendR
! - v   Istr :Iend , Jstr :JendR
! - u   Istr :IendR, Jstr :Jend  
! - p   Istr :IendR, Jstr :JendR

! BLXD In case of periodic Bdy two options regarding OA diags
! - calc inner domain points only then complete the bdy
! - calc inner domain points + bdy 


      ip_ran = (/I_PRI_RANGE/) 
      jp_ran = (/J_PRI_RANGE/)
      ir_ran = (/IstrR,IendR/) 
      jr_ran = (/JstrR,JendR/)
      iu_ran = (/IU_RANGE/)
      ju_ran = (/JU_RANGE/)
      iv_ran = (/IV_RANGE/)
      jv_ran = (/JV_RANGE/)


! BLXD rec. var3d_oa, var2d_oa calculated over :
! a) the only interior domain points 
! b) or also at physical Bdy points
! Option (b) chosen : requires to loop over indices
! rho : IstrR, IendR ...etc

! BLXD Warning Croco 2D masks (terrain-following coord) are of type real
!      OA masks are :
!      - 3D array of integer type 
!        to be used with geopotential vert. coord.
!      - an explicit int conversion is applied
!        below ifdef MASKING

      
      if ( size(hu_c, 1) .ne. ( ip_ran(2) - ip_ran(1)+1 )   .or.
     &     size(hu_c, 2) .ne. ( jp_ran(2) - jp_ran(1)+1 ) ) 
     & then
        MPI_master_only write(stdout,*) 
     &   'STOP : inconsistent private scratch size'
        stop
      endif
        
#ifndef MASKING
      maskr_c(PRIVATE_2D_SCRATCH_ARRAY,1:N) = 1
      masku_c(PRIVATE_2D_SCRATCH_ARRAY,1:N) = 1
      maskv_c(PRIVATE_2D_SCRATCH_ARRAY,1:N) = 1
      maskf_c(PRIVATE_2D_SCRATCH_ARRAY,1:N) = 1
#else

      maskr_c(PRIVATE_2D_SCRATCH_ARRAY,1:N) = 0
      masku_c(PRIVATE_2D_SCRATCH_ARRAY,1:N) = 0
      maskv_c(PRIVATE_2D_SCRATCH_ARRAY,1:N) = 0
      maskf_c(PRIVATE_2D_SCRATCH_ARRAY,1:N) = 0

! BLXD : change from only rho interior points to bdy points
!        compilation with -i4 -r8 options
!        if introducing kind precison INT(umask,kind=wp)
      do j=JU_RANGE
        do i=IU_RANGE
          masku_c(i,j,1:N) = INT(umask(i,j))
        enddo 
      enddo 
      do j=JV_RANGE
        do i=IV_RANGE
          maskv_c(i,j,1:N) = INT(vmask(i,j))
        enddo 
      enddo 

      do j=JV_RANGE
        do i=IU_RANGE
! BLXD ISSUE TO CHECK with the croco pmask at boundaries
! see grid.h setup_grid1.F
! pmask=(0=Land, 1=Sea, 1-gamma2 =boundary).
!      See grid.h pmask=(0=Land, 1=Sea, 1-gamma2 =boundary).
! using pmask instead
          maskf_c(i,j,1:N) = INT(pmask2(i,j))
        enddo 
      enddo 

      do j=JstrR,JendR
        do i=IstrR,IendR 
          maskr_c(i,j,1:N) = INT(rmask(i,j))
        enddo
      enddo

#endif

! BLXD add u-grid bdy points + bug with - sign 
      do j=JU_RANGE
        do i=IU_RANGE
          hu_c(i,j)=0.5*(h(i,j)+h(i-1,j))
        enddo 
      enddo 

! BLXD add v-grid bdy points + bug with - sign 
      do j=JV_RANGE
        do i=IV_RANGE
          hv_c(i,j)=0.5*(h(i,j)+h(i,j-1))
        enddo
      enddo

#ifndef SPHERICAL
! BLXD 2020 BUG change minus to plus sign
      do j=JU_RANGE
        do i=IU_RANGE
          xu_c(i,j)=0.5*(xr(i,j)+xr(i-1,j))
          yu_c(i,j)=0.5*(yr(i,j)+yr(i-1,j))
        enddo
      enddo
      do j=JV_RANGE
        do i=IV_RANGE
          xv_c(i,j)=0.5*(xr(i,j)+xr(i,j-1))
          yv_c(i,j)=0.5*(yr(i,j)+yr(i,j-1))
        enddo
      enddo
#endif

! BLXD GRID OPTION 2
! the grd structure could eventually be constructed here
! call allocate_tile_grid_tmp_str
! using its pointer fields
! grd(tile)%lon_t(PRIVATE_2D_SCRATCH_ARRAY) => lonr(PRIVATE_2D_SCRATCH_ARRAY)
! grd(tile)%mask_t(PRIVATE_2D_SCRATCH_ARRAY,1:N) => rmask(PRIVATE_2D_SCRATCH_ARRAY)
! however if the user namelist leads to "no analysis"
! this will be useless


      call init_oa( tile=tile+1
     & ,iic_oa=iic
     & ,imin=Istr-2, imax=Iend+2
     & ,jmin=Jstr-2, jmax=Jend+2
     & ,kmin=1,    kmax=N
#ifdef SPHERICAL
     & ,lon_t_oa=lonr(PRIVATE_2D_SCRATCH_ARRAY)
     & ,lat_t_oa=latr(PRIVATE_2D_SCRATCH_ARRAY)  
     & ,lon_u_oa=lonu(PRIVATE_2D_SCRATCH_ARRAY)  
     & ,lat_u_oa=latu(PRIVATE_2D_SCRATCH_ARRAY)  
     & ,lon_v_oa=lonv(PRIVATE_2D_SCRATCH_ARRAY)  
     & ,lat_v_oa=latv(PRIVATE_2D_SCRATCH_ARRAY) 
     & ,lon_f_oa=lonr(PRIVATE_2D_SCRATCH_ARRAY)    ! Wrong grid !
     & ,lat_f_oa=latr(PRIVATE_2D_SCRATCH_ARRAY)    ! Wrong grid !
#else
     & ,lon_t_oa=xr  (PRIVATE_2D_SCRATCH_ARRAY)
     & ,lat_t_oa=yr  (PRIVATE_2D_SCRATCH_ARRAY)    
     & ,lon_u_oa=xu_c(PRIVATE_2D_SCRATCH_ARRAY)  
     & ,lat_u_oa=yu_c(PRIVATE_2D_SCRATCH_ARRAY)  
     & ,lon_v_oa=xv_c(PRIVATE_2D_SCRATCH_ARRAY)  
     & ,lat_v_oa=yv_c(PRIVATE_2D_SCRATCH_ARRAY) 
     & ,lon_f_oa=xr  (PRIVATE_2D_SCRATCH_ARRAY)      ! Wrong grid !
     & ,lat_f_oa=yr  (PRIVATE_2D_SCRATCH_ARRAY)      ! Wrong grid !
#endif
     & ,mask_t_oa=maskr_c(PRIVATE_2D_SCRATCH_ARRAY,1:N)   
     & ,mask_f_oa=maskf_c(PRIVATE_2D_SCRATCH_ARRAY,1:N)  ! Wrong grid !
     & ,mask_u_oa=masku_c(PRIVATE_2D_SCRATCH_ARRAY,1:N)  
     & ,mask_v_oa=maskv_c(PRIVATE_2D_SCRATCH_ARRAY,1:N) 
     & ,h_w_oa=h         (PRIVATE_2D_SCRATCH_ARRAY)               
     & ,h_u_oa=hu_c      (PRIVATE_2D_SCRATCH_ARRAY)            
     & ,h_v_oa=hv_c      (PRIVATE_2D_SCRATCH_ARRAY)            
     & ,h_f_oa=h         (PRIVATE_2D_SCRATCH_ARRAY)      ! Wrong grid !
     & ,rhp_t  =rho   (PRIVATE_2D_SCRATCH_ARRAY,1:N)      
     & ,depth_t=z_r   (PRIVATE_2D_SCRATCH_ARRAY,1:N)     
     & )

C$OMP BARRIER
        call deallocate_tile_grid_tmp_str()

!CRBS print*, 'ooo OUT OF init_online_spectral_diags_tile ', mynode !CRBE

! ------- ------------- ------ ------ - -----
      return 
      end subroutine init_online_spectral_diags_tile

! BLXD GRID OPTION 2
! would require grd pointers twd null()
! grd structure deallocation

#else /* ONLINE_ANALYSIS */
      subroutine init_online_spectral_diags_empty
      end subroutine init_online_spectral_diags_empty
#endif /* ONLINE_ANALYSIS */
!
!----------------------------------------------------------------------
!
! Running online spectral analysis 
!
!----------------------------------------------------------------------
!
#ifdef ONLINE_ANALYSIS

      subroutine run_online_spectral_diags_tile (Istr,Iend,Jstr,Jend
     &                                          ,tile)
!
!
!#ifdef NBQ
!      use module_nh 
!      use module_nbq
!#endif
      use module_interface_oa, only : main_oa, scalogram_analysis
      use module_tile_oa, only : deallocate_tile_space_str
     &                          ,deallocate_tile_sclg_str
     &                          ,deallocate_tile_test_str
      use module_oa_variables, only : test_analysis

! BLXD WARNING TESTING 
! use scalars => i) use param ii) include "scalars.h" 
! use param   => include "param.h"
! REQUIRES to comment
! - include "param.h"
! - include "scalars.h"
!      use scalars
      implicit none
      integer, intent(in) :: Istr, Iend, Jstr, Jend
      integer, intent(in) :: tile
! Local variables
      integer :: i,j,k    

! >> Include with variable declaration
# include "param.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "grid.h"
#ifdef NBQ
# include "nbq.h"
#endif
! >> Include with computation instructions
# include "compute_auxiliary_bounds.h"
! From param.h, MPI subdom., tile 
! From compute_tile_bounds.h with Istr,Iend,Jstr,Jend
! Computes Physical Bdy indices IstrR, IendR, JstrR, JendR
!          Interior U or V points (start) IstrU, JstrV
!
! Does not require #include def_bounds.h
! LOCALLM, LOCALMM not needed at tile level

! BLXD removing "node" variables set at initialisation
!      which are common to each MPI subdomain
!      - mynode, if_print_mode... etc
!      removing simulation constant 
!      - dti, nt_max (cannot be adapative time step unless restarting)
!      Tile dimensions are already stored in the module_tile_oa 
!      see st%imin,st%imax,...
!      TODO check if removing imin,imax,... arguments 

!$OMP MASTER
#ifdef OA_TRACES
         print*,
     &   'eee ENTERING main_oa iic/node ',iic,mynode
#endif
!$OMP END MASTER


      call main_oa( tile=tile+1
     & ,ichoix=0
     & ,ivar_m=0
     & ,iic_oa=iic
     & ,imin=Istr-2, imax=Iend+2
     & ,jmin=Jstr-2, jmax=Jend+2
     & ,kmin=1,    kmax=N
     & ,rhp_t  =rho(PRIVATE_2D_SCRATCH_ARRAY,1:N)      
     & ,depth_t=z_r(PRIVATE_2D_SCRATCH_ARRAY,1:N)     
     & )

!$OMP MASTER
#ifdef OA_TRACES
         print*,
     &   'ooo OUT OF main_oa iic/node ',iic,mynode
#endif
!$OMP END MASTER

! At the very end icc = ntimes
! BLXD omp_barrier should be either here OR at the end of main_oa
C$OMP BARRIER
!        WARNING after calls to deallocate_tile_varspace_oa
!        in module_oa_interface with test on icc_oa as well
      if ( iic.eq.(ntimes+1) ) then
!$OMP MASTER
#ifdef OA_TRACES
         print*,'ooo END ONLINE_ANALYSIS icc/node ',iic,mynode
#endif
!$OMP END MASTER
         call deallocate_tile_space_str()
         if (scalogram_analysis) call deallocate_tile_sclg_str()
         if (test_analysis) call deallocate_tile_test_str()
         ! BLXD tile_isopycne_str not yet defined
      endif

! ------- ------------- ------ ------ - -----
      return 
      end subroutine run_online_spectral_diags_tile
#else /* ONLINE_ANALYSIS */
      subroutine run_online_spectral_diags_empty
      end subroutine run_online_spectral_diags_empty
#endif /* ONLINE_ANALYSIS */


