#include "cppdefs.h"
#ifdef ONLINE_ANALYSIS
!------------------------------------------------------------------------------
!                               NHOMS
!                Non Hydrostatic Ocean Modeling System      
!------------------------------------------------------------------------------
!
!> @note <a href="http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/index_snh_home.htm"> Main web documentation </a>
!!
!!     equipe d'oceanographie cotiere du laboratoire d'aerologie
!!     laboratoire d'aerologie
!!     cnrs - universite paul sabatier - observatoire midi pyrenees
!!     14 avenue edouard belin - 31400 toulouse - france
!
! DESCRIPTION: Module d'analyse spectrale et d'analyse en ondelettes en ligne pour Croco 
!
!> @brief OA Module : Online spectral and wavelet analysis module 
!
!> @details Procedures included in the module are:
!!
!! - init_parameter_oa       : reads the user namelist with global parameters (configuration, variable, domain, scales, frequency,..) 
!! - init_oa                 : initialization of the spatial, temporal and scale structures needed to stack/unstack analyses.
!! - rdnotebook_init_oa      : reading of the 1st namelist which defines the config-variables to analyse during the simulation.
!! - rdnotebook_oa           : reading of specific namelists each describing a perticular OA analysis to perform.
!! - init_configuration_oa   : example of initialization of "mixed" variables configuration (composite).
!! - var_grid_oa             : sets the grid point type for each defined OA variable (C grid).
!! - var_space_oa       : initializes the spatial structure of the state vector to analyse.
!! - var_per_oa         : initializes the frequency parameters of the state vector to analyse. 
!! - var_time_oa        : initializes the time parameters of the state vector to analyse.
!! - var_rec_oa         : calculates the "reconstruction factor" for Morlet wavelet, Fourier and Dirac (replaces inverse transformation).
!! - var_copy_oa        : used for "mixed" variables (composite), duplicates the variable information in another variable. 
!! - upd_init_oa        : prepares output arrays var2d_oa and var3d_oa and related parameters.
!! - main_oa            : main routine which performs the time convolution between the requested field and the OA atom (Fourier, wavelet, Dirac)
!!                        method without pointer, calls function var_oa.
!! - test_oa            : OA test functions application.
!! - psi_oa             : wavelet fonction
!! - psi_p2s_oa         : function which transforms the time period to the wavelet scale (s).
!! - var_oa             : function which returns the field value at one grid point i,j,k according to the analysis code
!! - var_upd_oa         : updates output arrays var2d, var3d with the analysis when available.
!! - box_oa             : heisenberg boxes (OA version prior to 2007)
!! - lev_init_oa        : initializes of the isopycne levels to analyse their positions.
!! - lev_upd_oa         : updates the wlev structured type with the curret position of the target density.
!! - update_level_oa    : searches the depth interval where the a target density is reached.
!! - subsave_oa         : calls var_upd_oa (outputs in file have been removed).
!! - allocate__*_oa     : dynamic variable allocation.
!! - deallocate_*_oa    : deallocation.
!! REMOVED module routines :
!! - struct_oa          : REMOVED (sauvegarde de la structure spatio-temporelle de la configuration).
!! - pressure_i_oa      : REMOVED (fonction pour le calcul de l'anomalie de pression).
!! - subsave_init_oa    : REMOVED (preparation des fichiers de sauvegarde)
!! - test_init_oa       : REMOVED (fonction test) 
!! ADDED module routines
!! - count_nmsimult_oa  : pre-calculates the state vector dimension as the number of simultaneously openened convolution windows.
!! - user_count_nmsimult_oa : namelist_oa nmsimult_oa_max /= -99 => user settings
!! - init_scalogram_oa  : initializes scalogram parameters to handle analyses outputs with possible distribution among MPI subdomains.  
!! - sclg_init_coords   : prepares the outputs of scalogram position and period of analysis handling MPI/tile(-threads) subdomain decomposition.
!! - init_temporal_box_oa : initialise Heisenberg boxes.
!! - calc_temporal_box_oa : calculate Heisenberg boxes (iterative, any atom).
!! - count_tracked_isopycne_oa
!! ADDED sources :
!! - output_oa          : OA updates global output array looping on tiles to finally send global array to XIOS
!! - var2d_oa, var3d_oa_out : called by output_oa
!! - module_parameter_oa: recycled module to allocate var3d_oa, var2d_oa to the croco GLOBAL_3D_ARRAY, 2D_ARRAY 
!! - online_spectral_diags : main routine called by Croco main.F and step.F 
!!                           supporting tile-threads and possible dual OpenMP-MPI parallelization of the horizontal domain
!! - module_tile_oa     : tile(-threads) compliant state vector structure
!! - module_grd_oa      : tile(-threads) compliant grid-derived parameters
!! REMOVED or fully modified related module :
!! - module_oa_upd, module_oa_stock, module_parameter_oa
!
! REVISION HISTORY:
!
!> @authors 
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version (2006)
!! - B. Lemieux-Dudon
!!  - Developments to perform online Scalogram analyses (06/2021)
!!   - from namelist options to calculation and XIOS netcdf outputs
!!   - scalograms distribution over MPI process/subdomains handled 
!!      either by XIOS or by MPI internal OA instructions (see if_mpi_oa)
!!   - Implementation of the Heisenberg theoritical and numerical uncertainties
!!      Onging. 
!!  - Toward a Croco-OA tile(-threads) compliant OA version (2021)
!!    supporting tile(-threads) and possible dual OpenMP-MPI parallelization 
!!    of the horizontal domain
!!    - Fields wf_oa%coef or wlev_oa%rhpt,... are domain dependant... 
!!      so tile(-threads) dependant.
!!    - wf_oa and wlev_oa now declared as field of the derived type st in module_tile_oa
!!      Croco tile(-thread) loop => array of derived type st(tile)% 
!!    - Croco grid passed as OA subroutine args are also tile-thread dependant
!!      as well as all the OA parameters related to the grid position of analyses
!!      Croco tile(-thread) loop => array of derived type grd(tile)% 
!!    See :
!!     => online_spectral_diags : main routine called by Croco main.F and step.F
!!     => module_tile_oa : tile(-threads) compliant state vector structure
!!     => module_grd_oa  : tile(-threads) compliant grid-derived parameters
!!     => openMP instructions
!!  - Handling netCDF outputs (of OA analyses) with specific routines which 
!!    are tile-thread compatible
!!     => output_oa : strictly updates array tiles to finally send the global array to XIOS
!!     => var2d_oa_out, var3d_oa_out : 2D using Croco working arrays work2d, work2d2
!!     => module_parameter_oa : allocating var3d_oa, var2d_oa to the croco GLOBAL_3D_ARRAY, 2D_ARRAY 
!!  - Memory optimization with pre-calculation of the state vector dimendion 
!!    to reduce the allocated memory (2021)
!!  - Croco-OA interface 1st version + XIOS OA outputs (2020)
!!    - Mix between the Stand-alone and Original OnlineA versions 
!!  - Memory optimization for the isopycnal level tracking (2015)
!!  - Headers, comments, cleaning, small changes from f77 to earlier standards  
!!  - Namelists (06/2014) + Stand-alone version + optimization (01/2015)
!
!> @date 2021
!> @todo BLXD 2020 Further Modify the Croco-OnlineAnalysis module interface?
!!   Interface is still a mix between the Stand-alone and Original OnlineA versions 
!!   => Part of the Croco fields are passed as arguments to the OA routines 
!!   using reduced dimension range arrays with memory duplication.
!> @todo BLXD 2020 general requirement to organize precision consistency
!!   beyond the use of the double precision compilation option.
!!    (ifort, gfortran jobcomp -r8 + mpc.F preproc not applied to F90)
!!    Module with wp declaration and use of KIND(wp) instruction
!!    to get rid of double precision type (obsolete), properly define 
!!    the cplx type precision,... etc.
!> @todo BLXD 2021 implement the Roms-Croco outputs based on wrt_his.F and 
!!    nf_write.F routines ?
!> @todo BLXD 2021 dev. to handle temporal convolution windows under calculation 
!!    in the case of simulation with restarts ?
!> @todo BLXD 2021 tests to perform with a working Croco tile(-thread) version
!> @todo BLXD 2021 stop on error with comments to review (printing subprocess)
!!       threads, standard output or io unit 
!> @todo BLXD 2021 tests scalogram with undef MPI cppkey 
!------------------------------------------------------------------------------

      module module_interface_oa

      implicit none

#ifdef MPI
 include 'mpif.h'
#endif


      !> MPI parameters for the Croco-OA interface
      integer, parameter :: root=0
      integer :: mynode, comm, comm_size
      logical :: if_print_node, if_mpi_oa

      !> Test if OA analyses requested 
      logical             :: if_oa

      !> Flag controling the type of OA outputs (forced to .true. august 2021 : XIOS version 2.5 only)
      logical             :: if_xios_output

      !> First and last simulation iteration indices
      integer :: kount0, nt_max                                      

      !> External time step
      double precision :: dti

      !> Shared Croco log io_unit
      integer :: io_unit
 
      !> directory_in_oa hardcoded to INPUT, directory out_oa set in the namelist_oa 
      character(len=200)  :: directory_in_oa, directory_out_oa
      integer             :: verbose_oa
      integer, parameter  :: io_nml=90, io_hist=91, io_hbg=92
      character(len=1), parameter :: txtslash='/'
      character(len=250)  :: file_hbg

      !> Namelist options including Heisenberg uncertainties  
      logical             :: if_sphe_deg2rad, if_extend_dom, if_chks
      logical             :: if_heisenberg_box
      real                :: unite_oa    !< namelist parameters for OA time unit 
      real                :: pi_oa       !< TODO? double precision, real(8)

      !> Test functions (BLXD TODO see if declared elsewhere, rm amp_test2_oa ?) 
      integer :: nper_test
      double precision, allocatable, dimension(:) :: period_test_oa, amp_test_oa, amp_test2_oa

      !> Maximum number of predefined configuration/var types
      integer, parameter :: maxtyp_oa = 100

      integer, parameter :: nmsgnv_oa = 99 ! #BLXD set the # of possible config. maxtyp_oa   

      !> Maximum number of analysis of the same type requested simultaneously
      integer, parameter :: max_idcfg_oa = 20                           
      
      !> Maximum number of mixed variables (composites)
      integer, parameter :: nmvar_oa = 5                              

      !> Flag for field dimension number (2D/3D) and type of Arakawa C-grid (brought from module_oa_space)
      integer,dimension(maxtyp_oa) ::                                  &  ! changed 100 to 200 jwf 20070619
            tgv3d_oa                                                   &  !
            ,tgv_oa
       
      character(len=5),dimension(maxtyp_oa) :: tgvnam_oa

      !> Number of variables associated to configuration (brought from module_oa_variable)
      integer:: nzvc_oa(maxtyp_oa)

      !> Last dimension of the allocatable array var?d_oa
      !! global dimension encompassing all the analysis shared by all tile/MPI process subdomains
      integer :: nzupd3d_oa, nzupd2d_oa, nzupd1d_oa, nzupd0d_oa
 
      !> Analysed 3D fields
      !! Dimensions are space dimensions i, j, k, and the index of the analysed variable has returned by tvar_oa
      complex(8), dimension(:,:,:,:), allocatable :: var3d_oa
      
      !> Analysed 2D fields
      !! Dimensions are space dimensions i, j, k, and the index of the analysed variable has returned by tvar_oa
      complex(8), dimension(:,:,:), allocatable   :: var2d_oa

      !> Returns the index of the analysed variable (required for var3d_oa and var2d_oa)
      !! dimension 1 : OA configuration-variable as requested in the namelist, tc_oa(ic), e.g., isopycne analysis has OA analysis index 20 
      !! dimension 2 : if several OA configuration-variable of the same type are requested, their rank, tvc_oa(ic) 
      !!               e.g. 1 for first OA analysis 20, 2 for the second OA analysis 20,... etc.
      !! dimension 3 : # of variables involved in the configuration (see composite config. with several variables)
      ! integer, dimension(1:maxtyp_oa,1:max_idcfg_oa,1:10)  :: tvar_oa
      integer, dimension(1:maxtyp_oa,1:max_idcfg_oa,1:nmvar_oa)  :: tvar_oa

      !> Scalogram analysis (BLXD TODO declared with test_analysis, isopycne_analysis?) 
      logical             :: scalogram_analysis

      real, parameter :: w0c=6.

#ifdef MPI
      !> Analysed 0D fields for scalogram global/local to tile-MPI process scal0d_oa
      complex(8), dimension(:,:), allocatable, target   :: scal0d_oa
      !> Array with Position/Scales of local to tile-MPI process scalogram
      real(8), dimension(:,:) , allocatable, target  :: per0d_oa
      integer, dimension(:)   , allocatable          :: iscal0d_oa
      integer, dimension(:)   , allocatable          :: jscal0d_oa
      integer, dimension(:)   , allocatable          :: kscal0d_oa
      integer, dimension(:)   , allocatable          :: index_s_oa
#endif
      integer, dimension(:) ,   allocatable  :: v2locsclg
      !integer, dimension(:) ,   allocatable  :: sclg2v !BLXD tmp keep

      complex(8), dimension(:,:), allocatable   :: scal0d_cr
      real(8), dimension(:,:)   , allocatable   :: per0d_cr
      integer, dimension(:)     , allocatable   :: iscal0d_cr
      integer, dimension(:)     , allocatable   :: jscal0d_cr
      integer, dimension(:)     , allocatable   :: kscal0d_cr
      integer, dimension(:)     , allocatable   :: index_s_cr


      !> Number of scalograms local to the current MPI process if any
      !! For now the total number of scalogram requested by the user is nzupd0d_oa
      !! BLXD TODO change the name eventually to nsclg_glo
      integer                                   :: nsclg_loc

      !> Firts global scalogram position on the current MPI process if any 
      !! (see XIOS outputs index_s indexing for scalogram distributed among MPI processes)
      integer                                   :: isclg_beg
  
      !> Namelist options for scalogram XIOS outputs (namelist_oa) 
      logical                                   :: if_record_sclg_ijpoints, if_xios_dom_grid

      ! BLXD bUG size moved from module variable to module_interface_oa
      !      integer, dimension(1:6000) :: des_oa  
      !
      ! Variable useful for NHOMS energy analysis if you wish to call OA at specific point in the code
      ! for specific variable that are updated at some irregular code place or time step (e.g., time splitting)
      ! depnds on iv_m so it should have nzv_oa size has cnb_oa 

      ! BLXD TODO integer, dimension(1:nzv_oa) :: io_traces

#ifdef OA_TEST_MPI_XIOS
      integer, dimension(:)     , allocatable   :: cntvt_cr
#endif

      !> Counting tiles for tile(-threads) and openMP ATOMIC tag
      integer :: tile_count_oa
      integer :: itgt_glob, jtgt_glob, ktgt

#ifdef MPI
      complex(8), dimension(:), pointer      :: buffx_s => null()
      real(8),    dimension(:), pointer      :: buffr_s => null()
      complex(8), dimension(:), allocatable  :: buff2x_s
      real(8),    dimension(:), allocatable  :: buff2r_s, buff2r_r
      !complex(8), dimension(:), pointer  :: buffx_r => null()
      integer,    dimension(:), pointer      :: buffind_s, buff
      integer,    dimension(1:3) :: buffij_s, buffij_r
      integer, dimension(1:MPI_STATUS_SIZE) :: mpi_status 
      integer, dimension(:), allocatable  :: sclg_request
      !integer, dimension(:), allocatable  :: sclgij_request, recv_request
      logical, parameter :: mpi_nonblock_send =.false.
      logical, parameter :: mpi_blocking_recv = .true.
      logical, parameter :: mpi_nobuff_for_per_recv = .true. !Success with .false. and blocking recv
                                                         
      logical, parameter :: mpi_use_test=.false.
      logical, parameter :: mpi_send_buff_ptr=.false.
#ifdef OA_TEST_MPI_XIOS
      ! BLXD test_mpi TODO REMOVE
      !logical, parameter :: if_test_mpi_blocking_send=.false., if_test_test_mpi_nonblocking_ssend=.false.
#endif
#endif

      contains

!----------------------------------------------------------------------
! PROCEDURE
!
!> @note Croco tile(-threads) : this routine must be called before any
!!       openMP parallel region fork and/or whitout loop over tiles
!!       so that the only MASTER thread initializes the implicitly 
!!       openMP SHARED module variables without any data race issue.
!
! DESCRIPTION: lecture des namelists pour initialiser les parametres
!  d'analyse (configutation, variable, frequence, echelles analysees).
!  connus de tous les process MPI et Croco tile(-threads)
!
!> @brief reads OA namelists and set global analysis parameters.
!
!> @details sets the type of OA analysis : atom (Fourier, wavelet,..)
!!  period/scale of analysis, ocean field target, composite or single
!!  variables, scalogram, reconstructed wavelet coeff,...
!!  characteristic of the time convolution windows. Set the targets :
!!  ocean field (composite or single variable), geographic area. 
!!  These parameters are global and known to all the MPI processes 
!!  and tile(-threads) if any. This routine is called in Croco main.F 
!!  before any openMP parallel region and outside any tile(-thread) loop.
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists + Stand-alone version + optimization (2015)
!!  - Headers, comments, cleaning, from f77 to earlier standards  
!!  - Croco-OA interface 1st working version (2020)
!!  - Toward Croco Tile(-threads) compliant interface and OA version (2021)
!!  - XIOS OA output (2021)
!!  - Developments to perform online Scalogram analyses (2021)
!!    - from namelist options to calculation and XIOS netcdf outputs
!!    - scalograms distribution over MPI process/subdomains handled 
!!      either by XIOS or by MPI internal OA instructions.
!!  - Implementation of the Heisenberg theoritical and numerical uncertainties. 
!!    Onging (2021). 
!> @date 2021
!> @todo BLXD
!! - Specific module needed to set consitent precisions with real(KIND=wp)...
!! - Outputs using Roms-Croco wrt_his/nf_write routines
!----------------------------------------------------------------------
!
      !subroutine finalize_parameter_oa
      !
      !use module_tile_oa, only : deallocate_tile_space_str
      !implicit none 
      !return
      !end subroutine finalize_parameter_oa

!     logical function init_parameter_oa(     &
      subroutine init_parameter_oa(     &
     &  io_unit_oa                      &
     & ,if_print_node_oa                &
     & ,mynode_oa                       &
     & ,comm_oa                         &
     & ,dti_oa                          & 
     & ,kount0_oa                       &
     & ,nt_max_oa                       &
     & ,ntiles)

      use module_oa_variables
      use module_oa_time
      use module_oa_space 
      use module_oa_periode
      !use module_oa_stock
      use module_tile_oa, only : allocate_tile_space_str, allocate_tile_sclg_str &
     &                          ,allocate_tile_test_str 
      use module_grd_oa, only  : allocate_tile_grid_tmp_str
      !use module_isopycne_oa, only  : allocate_tile_isopycne_str

      implicit none

      !> Time integration step (s)
      double precision, intent(in) :: dti_oa

!     !> Current model integration iteration
!     integer, intent(in) :: iic_oa 

      !> First and last simulation iteration indices
      integer, intent(in) :: kount0_oa, nt_max_oa                                       

      !> Number of Croco tiles
      integer, intent(in) :: ntiles                                       

      integer, intent(in) :: io_unit_oa, mynode_oa, comm_oa
      logical, intent(in) :: if_print_node_oa

      !integer, intent(inout) :: nper_oa
      integer :: nper_oa

      integer :: ic_s

#ifdef MPI
      integer :: ierr
#endif


      ! Parameter initializion
      ! TODO remove from module_oa_interface 
      pi_oa = acos(-1.)

      directory_in_oa  = "INPUT"

      ! Test configuration set to .false. by default
      test_analysis = .false.

      ! Isopycne configuration set to .false. by default
      isopycne_analysis = .false.
            
      ! Scalogram set to .false. by default
      scalogram_analysis = .false.

      ! BLXD pointers
      io_unit = io_unit_oa
      mynode  = mynode_oa
      comm    = comm_oa
      if_print_node = if_print_node_oa
      dti = dti_oa
      kount0 = kount0_oa 
      nt_max = nt_max_oa

      !open(unit=io_unit,file=trim(directory_out)//'/'//'trace.out')    
      if(if_print_node) write(io_unit,*) ' '
      if(if_print_node) write(io_unit,*) '..Begin Online Analysis OA LOG FILE'

!.....nombre de variables par configuration:
!     de 1 à 99 => configuration n'impliquant qu'une seule variable
!     sup a 100 => conf. a +s variables
!     BLXD TODO Specific init routine for nzvc_oa ?  
!     From 100 to parameter maxtyp_oa=200 => configuration including several variables

      do ic_s=1,nmsgnv_oa
         nzvc_oa(ic_s) = 1
      enddo
      nzvc_oa(100) = 2 

!..... BLXD removing des_oa initialization : 
!      WRONG type (integer) - WRONG array size 1:6000 instead of 1:nzv_oa
!..... No energy variable so far!
!      des_oa=0.

!.....lecture des nzc_oa configuration et calcul du nombre total de var. associe  nzv_oa
!     BLXD include composite configuration with several variables requires nzvc_oa
      call rdnotebook_init_oa


      if(if_print_node) write(io_unit,*) '...OA LOG FILE : # of requested analysis ',nzv_oa

!.....test si calcul oa ou pas
      if_no_analysis_requested : if(nzv_oa.eq.0) then

        if(if_print_node) write(io_unit,*) '...OA LOG FILE : WARNING zero analysis requested (nzv_oa=0)'

!.....flag returning false to croco, either function or module flag
!        init_parameter_oa = .false.
         if_oa         = .false.
         nper_oa       = 1  ! XIOS axis s_oa minimum size 1

        return

      else


        if(if_print_node) write(io_unit,*) '...OA LOG FILE : preliminary allocation '

!.....allocations dynamiques des parametres pour analyse ssi nzv_oa /= 0
        call allocate_part1_oa( )


!.....BLXD MOVING des_oa INITIALISATION HERE with correct type (integer) an array size (nzv_oa)
!     just allocated in allocate_part1
!       No energy variable so far!
        des_oa(1:nzv_oa) = 0 

!.....definition des ondelettes de morlet:

        ! BLXD appropriate place to init wavelet param, namelist can be a better place ?
        ! TODO use pi_oa or croco pi parameter
        fb_oa=2.                ! definition de l''ondelette de morlet complexe
        !fc_oa=6./2./3.14159274  ! par defaut fb_oa=2. et fc_oa=6./3.14159274/2.
        !fc_oa=6./2./pi_oa w0c redef. as module param
        !fc_oa=6.145/2./pi_oa w0c redef. as module param
        fc_oa=w0c/2./pi_oa

        if(if_print_node) write(io_unit,*) '...OA LOG FILE : Specific analysis config. soon read in namelists '

!.....lecture du notebook.oa: 
!      #BLXD removing the use of 2D grid parameter. Still contains the call to var_copy_oa
        call rdnotebook_oa( )

!.....Currently output can only be handled by XIOS - Croco-Roms netCDF output types to dev.
        if ( .not. if_xios_output ) then
            write(*,*) '...OA ERROR : outputs only handled by XIOS - Roms-Croco outputs to develope'
            stop
        endif

!.....Currently output can only be handled by XIOS version >= 2.5 - backward compat. to DEV
#ifndef XIOS2
        if ( if_xios_output ) then
            write(*,*) 'ERROR OA : XIOS version >= 2.5 required to output Online Analysis'
            stop
        endif   
#endif
     
!.....netCDF scalogram XIOS outputs applies a domain-based definition grid
!     (scalogram_period,scalogram_index) -> scalogram_value
!     XIOS2.5 axis-based definition grid to dev/test
      if (scalogram_analysis) then
        if ( if_xios_output .and. ( .not. if_xios_dom_grid ) ) then
            write(*,*) '...OA ERROR : XIOS grid axis not available => please change the OA namelist with if_xios_dom_grid set to T'
            stop
        endif
      endif

#ifndef MPI 
      if_mpi_oa = .false.
      if(if_print_node) write(io_unit,*) '...Non-MPI run => forcing to .false. flag if_mpi_oa '
#else
!.....BLXD set communicator size for further use     
      call MPI_COMM_SIZE(comm, comm_size, ierr)

      if (scalogram_analysis) then
        if ( if_mpi_oa ) then
            write(io_unit,*) '...OA WARNING : XIOS will not handle MPI scalogram distribution => check computing time' 
        else
            write(io_unit,*) '...OA WARNING : XIOS will handle MPI scalogram distribution' 
        endif
      endif
#endif

        call var_grid_oa

!.....BLXD here we know the # of tiles ntiles from module_tile_oa
        call allocate_tile_space_str(ntiles)
        call allocate_tile_grid_tmp_str(ntiles)

!.....BLXD scalogram analysis 
        if(scalogram_analysis) then
            if(if_print_node) write(io_unit,*) '...OA in : scalogram structure allocation ',mynode
            call allocate_tile_sclg_str(ntiles)
           if(if_print_node) write(io_unit,*) '...OA out : scalogram structure allocation ',mynode
        endif

!.....BLXD counting # of tracked isopycne (ispopycne_analysis true => nzlevel_oa>0)
        if ( isopycne_analysis ) then 
            call count_tracked_isopycne_oa()
        endif

!.....BLXD test analysis 
        !if(test_analysis) then
        if (test_analysis) then
            if(if_print_node) write(io_unit,*) '...OA in : test structure allocation ',mynode
            call allocate_tile_test_str(ntiles)
            if(if_print_node) write(io_unit,*) '...OA out : test structure allocation ',mynode
        endif

        call var_per_oa( .false., dti )
        
        call allocate_part3_oa   
        
        call var_per_oa( .true., dti )

        nper_oa       = max(nper_sclg_max,1)  ! BLXD TODO SEE if finally USED XIOS axis s_oa minimum size 1

        if(if_print_node) then
            if(scalogram_analysis) then
                if(if_print_node)then
                write(io_unit,*) '...OA LOG FILE : Scalogram Requested '
                write(io_unit,*) '...OA LOG FILE : Scalogram max # of period ', nper_sclg_max
                endif
            else
                if(if_print_node)then
                write(io_unit,*) '...OA LOG FILE : NO Scalogram  => max # of period set to 0      = ', nper_sclg_max
                write(io_unit,*) '...OA LOG FILE : NO Scalogram  => XIOS period axis dim set to 1 = ', nper_oa
                endif
            endif
        endif

        ! Counting tiles in either init_oa or main_oa
        tile_count_oa=0

        ! Pre-calculation of last dimension of var?d/scal0d_oa : nzupd?d_oa
        call upd_init_oa(.false.)
   
!.....if function returning true : init_parameter_oa = .true.
         
        if_oa         = .true.

      end if if_no_analysis_requested

      return 
      end subroutine init_parameter_oa
!      end function init_parameter_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note Croco tile(-threads) : this routine must be called before any
!!       openMP parallel region fork and/or whitout loop over tiles
!!       so that the only MASTER thread initializes the implicitly 
!!       openMP SHARED module variables without any data race issue.
!
! DESCRIPTION: Lecture de la namelist namelist_oa avec le nombre de configurations d'analyse
!  leur nom donnant acces aux namelists specifiques à chaque analyse.
!
!> @brief Reads a first namelist which defines the number and names of all the OA analysis to 
!! simultaneously perform during the simulation, sets global user options and 
!! enables the reading of specific analysis namelist.
!!
!> @details The namelist name is "namelist_oa". It provides the following parameter value:
!! - nzc_oa        \: number of simultaneous OA analysis to perform during the simulation,
!! - nzc_oa_names  \: a coma seperated list of the specific names of the requested analysis,
!! - nzc_oa_vartyp \: a coma seperated list of the specific OA configuration code corresponding to
!! - nzc_oa_scalogram \: a coma separeted list saying if configuration is a scalogram analysis
!!   the requested analysis.
!! The number and names of OA analysis to perform are applied to read the supplementary namelists 
!! that specifically describe each type of analysis. 
!! The names of the specific namelist are constructed as follows:
!! - "namelist_" \+ nzc_oa_names 
!
! REVISION HISTORY:
!
!> @authors
!! - B. Lemieux-Dudon
!!  - Namelists + Stand-alone OA version (2015)
!!  - Headers, comments, cleaning, from f77 to earlier standards  
!!  - Croco-OA 1st working interface version (2020)
!!  - XIOS outputs + Test functions (2020, 2021)
!!  - Scalogram analyses (analysing ocean field time series).
!!    - Two options included to handle the distribution of the scalograms among subdomains
!!      using either OA internal MPI instructions and/or XIOS facilities.
!!    - Heisenberg boxes.
!> @date 2021
!> @todo BLXD 
!! - clean/comment further namelist options 
!------------------------------------------------------------------------------

      subroutine rdnotebook_init_oa

      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      !use module_oa_stock

      implicit none

      integer             :: ic_r  !< Dummy configuration index, test function # of harmonics 
      character(len=250)  :: filin !< namelist filename

      namelist / oa_numbers / nzc_oa
      namelist / oa_names   / nzc_oa_names, nzc_oa_vartyp, nzc_oa_scalogram
      namelist / oa_test    / nper_test
      namelist / oa_test_signal / period_test_oa, amp_test_oa 
      namelist / oa_checks / if_chks, verbose_oa, directory_out_oa &
                            ,itgt_glob, jtgt_glob, ktgt
      namelist / oa_params / unite_oa, if_sphe_deg2rad, nmsimult_oa_max, if_extend_dom     &
                            ,if_mpi_oa, if_xios_output, if_xios_dom_grid             &
                            ,if_record_sclg_ijpoints                                 &
                            ,if_heisenberg_box

      ! Namelist default values 
      ! 1) oa_test_signal block for testing new dev. with specific cosine combination
      nper_test = 0

      ! 2) oa_checks BLXD see if moving directory_out_oa to oa_params block
      directory_out_oa = "OUTPUT" 
      verbose_oa       = 6
      if_chks          = .false.

      itgt_glob=-999 ; jtgt_glob=-999 ; ktgt=-999

      ! 2) oa_params DEFAULT settings    
      unite_oa         = 1
      if_sphe_deg2rad  = .false.
      if_mpi_oa        = .false.
      if_xios_dom_grid = .true.
      if_xios_output   = .true.
      if_record_sclg_ijpoints = .true.
      if_heisenberg_box = .false.
      nmsimult_oa_max   = -99

      !if_send_xios_2Dsclg=.false.
      !if_send_zero_size_arr=.true.
      !if_block_zero_size_arr=.false.

      ! BLXD TODO CLEAN see if removed
      if_extend_dom   = .false.

      filin= trim(directory_in_oa) // txtslash // 'namelist_oa'

      open(unit=io_nml, file = filin)
      read(unit=io_nml, nml  = oa_numbers)

      oa_analysis_requested : if (nzc_oa /= 0) then
    
      call allocate_namelist_oa 

      read(unit=io_nml, nml  = oa_names)
     
      nzv_oa = 0

      do ic_r = 1 , nzc_oa
       
         if(if_print_node) write(io_unit,*) 'OA analysis names and variable types' 
         if(if_print_node) write(io_unit,*) ic_r, nzc_oa, nzc_oa_names(ic_r), nzc_oa_vartyp(ic_r), nzc_oa_scalogram(ic_r)

!........A une config ic_r peut-etre associe plusieurs variables 
!        nzvc_oa est le nombre de variables impliquees dans la config


         nzv_oa = nzv_oa + nzvc_oa( nzc_oa_vartyp(ic_r) )

        ! If one of several var-config code set to 99 => Test function
        ! A single analytical function can be tested against different atoms
        if ( ( nzc_oa_vartyp(ic_r)==99 .or. nzc_oa_vartyp(ic_r)==98 ) .and. (nper_test==0) ) then
          read(unit=io_nml, nml  = oa_test)
        end if

      enddo


      oa_test_requested : if (nper_test > 0) then
          
          call allocate_namelist_test_oa(nper_test) 
          read(unit=io_nml, nml  = oa_test_signal)

      else if  (nper_test == -1 ) then      
          ! BLXD nper_test = -1 hardcoded test function 
          if(if_print_node) write(io_unit,*) 'OA hardcoded test function'
          nper_test = 1 
          call allocate_namelist_test_oa(nper_test) 
          amp_test_oa(1)    = 0.D0 
          period_test_oa(1) = 0.D0

      endif oa_test_requested 

      if (nper_test > 0) then

          do ic_r = 1, nper_test
            if(if_print_node) write(io_unit,*) 'OA test function (period,amplitude) ', period_test_oa(ic_r), amp_test_oa(ic_r)
          end do

      end if

      read(unit=io_nml, nml  = oa_checks)

      read(unit=io_nml, nml  = oa_params)

      if ( nmsimult_oa_max/=-99) then
          if(if_print_node) write(io_unit,*) 'OA user settings MAX. SIZE OF SIMULT. conv. window (ic,iv,lt) ', nmsimult_oa_max
          if (nmsimult_oa_max<1) then
            if(if_print_node) write(io_unit,*) 'ERROR : OA invalid MAX. SIZE OF SIMULT. conv. window'
            stop
          endif
      endif

      endif oa_analysis_requested 

      close(io_nml)

      return
      end subroutine rdnotebook_init_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: allocate/deallocate namelist parameters 
!
!> @brief allocation of names and type of requested OA analysis
!
!> @details variable allocation according to nzc_oa the number of type of 
!!  requested OA analysis :
!! - nzc_oa_names  \: names of the requested analysis,
!! - nzc_oa_vartyp \: configuration code of the requested analysis.
!! - nzc_oa_scalogram \: scalogram configuration analyses
!
! REVISION HISTORY:
!
!> @authors
!! - B. Lemieux-Dudon
!!  - Namelists + comments + Stand-alone OA version (2015)
!!  - Scalogram analyses (2021).
!> @date 2021
!> @todo
!------------------------------------------------------------------------------

      subroutine allocate_namelist_oa

      use module_oa_variables , only : nzc_oa_names, nzc_oa_vartyp, nzc_oa_scalogram
      use module_oa_periode , only : nzc_oa

      implicit none

          integer :: izc_oa
          allocate( nzc_oa_names(1:nzc_oa) )
          allocate( nzc_oa_vartyp(1:nzc_oa) )
          allocate( nzc_oa_scalogram(1:nzc_oa) )
          do izc_oa=1,nzc_oa
          nzc_oa_scalogram(izc_oa) = .false.
          enddo

      end subroutine allocate_namelist_oa

      subroutine deallocate_namelist_oa

      use module_oa_variables , only : nzc_oa_names, nzc_oa_vartyp, nzc_oa_scalogram

      implicit none

          deallocate( nzc_oa_names )
          deallocate( nzc_oa_vartyp )
          deallocate( nzc_oa_scalogram )

      end subroutine deallocate_namelist_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: allocation/desallocation des parametres pour les fonctions test. 
!
!> @brief allocate/deallocate test function parameters (spatially uniform cosine LC) 
!
!> @details linear combination of cosines, parameters partly user-defined 
!   namelist in namelist_oa :
!! - n              \: number of cosine (nper_test),
!! - period_test_oa \: array of cosine period in unit_oa (namelist_oa) 
!! - amp_test2_oa   \: array of cosine amplitude
!
! REVISION HISTORY:
!
!> @authors
!! - B. Lemieux-Dudon
!!  - Namelists + comments + stand-alone version (2015)
!!  - Test functions (2021).
!> @date 2021
!> @todo clean amp_test2_oa, add temporal section markers ?
!------------------------------------------------------------------------------

      subroutine allocate_namelist_test_oa(n)

      implicit none
      integer, intent(in) :: n

          if (.not. allocated(period_test_oa) ) allocate( period_test_oa(1:n) )
          if (.not. allocated(amp_test2_oa) )   allocate( amp_test2_oa(1:n) )
          if (.not. allocated(amp_test_oa) )    allocate( amp_test_oa(1:n) )

      end subroutine allocate_namelist_test_oa

      subroutine deallocate_namelist_test_oa

      implicit none

          if ( allocated(period_test_oa) ) deallocate( period_test_oa )
          if ( allocated(amp_test2_oa) )   deallocate( amp_test2_oa   )
          if ( allocated(amp_test_oa) )    deallocate( amp_test_oa    )

      end subroutine deallocate_namelist_test_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: allocation/desallocation des parametres de namelist
!              relatifs a une variable-configuration d'analyse 
!  a given variable-configuration analysis.
!
!> @brief allocate/deallocate namelist parameters related to
!  a given variable-configuration analysis (e.g., rho, M2)
!
!> @details parameters are array of size nzv_oa 
!!  
!
! REVISION HISTORY:
!
!> @authors
!! - B. Lemieux-Dudon (2014)
!!   Adapted from Symphonie/NHOMS initial version (2006)
!!   Bug tvc_oa
!> @date 2021
!> @todo
!------------------------------------------------------------------------------

      subroutine allocate_part1_oa( )

      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      !use module_oa_stock

      implicit none

!     #BLXD old NHOMS allocation
!     allocate (rhphat_oa_t(1:imax,1:jmax,1:kmax))
!     allocate (temhat_oa_t(1:imax,1:jmax,1:kmax))
!     allocate (salhat_oa_t(1:imax,1:jmax,1:kmax))       

!.....Allocation module_oa_variables:
      
      if ( .not. allocated(swt_d_oa) ) allocate(swt_d_oa(nzv_oa))          !< caracteristiques spatiales de la variable (voir notebook_oa)
      if ( .not. allocated(swt_t_oa) ) allocate(swt_t_oa(nzv_oa))          !< caracteristiques temporel  de la variable (voir notebook_oa)
      if ( .not. allocated(tv_oa)    ) allocate(tv_oa   (nzv_oa))          !< code variable associe 
      if ( .not. allocated(cnb_oa)   ) allocate(cnb_oa  (nzv_oa))          !< call number: pos. of call regarding time splitting
      if ( .not. allocated(updv_oa)  ) allocate(updv_oa (nzv_oa))          !< flag de remise a jour    
      if ( .not. allocated(save_oa)  ) allocate(save_oa (nzv_oa))          !< flag de sauvegarde
      if ( .not. allocated(tupd_oa)  ) allocate(tupd_oa (nzv_oa))          !< pour variables communes
    
      if ( .not. allocated(tpsi_oa)        ) allocate(tpsi_oa (nzv_oa))             !< type d atome utilise

      if ( .not. allocated(ltrec_fst_oa    )) allocate(ltrec_fst_oa (nzv_oa)   )    !< pour garder le 1er OA record de la simu
      if ( .not. allocated(if_first_rec_oa )) allocate(if_first_rec_oa (nzv_oa))    !<...
      if ( .not. allocated(ltrec_lst_oa    )) allocate(ltrec_lst_oa (nzv_oa)   )    !< pour garder le 1er OA record de la simu
      if ( .not. allocated(if_last_rec_oa  )) allocate(if_last_rec_oa (nzv_oa) )    !<...
    
!     Bug corr. des_oa allocated to wrong size +  tvc_oa must have the size of the # of configurations + tvar_oa known size
      if ( .not. allocated(des_oa     )) allocate( des_oa(nzv_oa)    ) 

!     Scalogram variables
      if ( .not. allocated(tv_sclg_oa )) allocate(tv_sclg_oa(nzv_oa) ) 

!     Structure temporelle du vecteur d'etat    
      if ( .not. allocated(begvt_oa )) allocate(begvt_oa(nzv_oa+1))

!     Structure spatiale du vecteur d'etat : begvs_oa
!     BLXD croco2021 tile-thread => spatial structure now handled in module_tile_oa.F90
!      if ( .not. allocated( ) allocate(    

!.....Configuration dependent variables
      if ( .not. allocated(tc_oa)   ) allocate(tc_oa (nzc_oa)      ) 
      if ( .not. allocated(tvc_oa)  ) allocate(tvc_oa(nzc_oa)      )    ! If several config. of the type is requested
      if ( .not. allocated(begc_oa) ) allocate(begc_oa (nzc_oa+1)  ) 

!.....Allocation module_oa_space:

      if ( .not. allocated(ptij_oa)) allocate(ptij_oa (2,nzv_oa)   )  ! point particulier demande par l utilisateur

      if ( .not. allocated(lat_oa )) allocate(lat_oa  (2,nzv_oa)   )  ! latitude  min, max de la structure 2d du vecteur d etat
      if ( .not. allocated(lon_oa )) allocate(lon_oa  (2,nzv_oa)   )  ! longitude  ...

      if ( .not. allocated(h_oa  )) allocate(h_oa    (2,nzv_oa)   )  ! profondeur ...  
      if ( .not. allocated(k_oa  )) allocate(k_oa    (2,nzv_oa)   )  ! niveaux verticaux min et max de la structure 3d du vecteur d etat

      if ( .not. allocated(dx_oa )) allocate(dx_oa   (nzv_oa)     ) ! resolution horizontale suivant x demandee par l utilisateur
      if ( .not. allocated(dy_oa )) allocate(dy_oa   (nzv_oa)     ) ! resolution horizontale suivant y demandee par l utilisateur
      if ( .not. allocated(dk_oa )) allocate(dk_oa   (nzv_oa)     ) ! resolution verticale   suivant z demandee par l utilisateur
    

!.....Allocation module_oa_periode et module_oa_time:

      if ( .not. allocated(nzpt_per_oa)   ) allocate(nzpt_per_oa(nzv_oa)      )   ! nombre de points de discretisation par periode
      if ( .not. allocated(kount_user_oa) ) allocate(kount_user_oa(3,nzv_oa)  )   ! description des periodes choisies par l utilisateur
      if ( .not. allocated(t0_oa)         ) allocate(t0_oa (nzv_oa)           )   ! date de la premiere sortie

      if ( .not. allocated(swt_wfpf_oa)   ) allocate(swt_wfpf_oa (nzv_oa)     )   ! calcul du coef wf ou du spectre pf (choix utilisateur)

      if ( .not. allocated(fl_rec_oa)     ) allocate(fl_rec_oa   (nzv_oa)     )   ! flag de reconstruction

      if ( .not. allocated(dori_oa)       ) allocate(dori_oa      (nzv_oa)    )   ! configuration frequentielle choisie par l utilisateur
      if ( .not. allocated(delta_t_oa)    ) allocate(delta_t_oa   (nzv_oa)    )   !< nombre de periodes pour calculer le produit de convolution


      if ( .not. allocated(per_oa)        ) allocate(per_oa (3,nzv_oa)        )   ! preriodes min, max, delta de chaque variable (ou configuration)

!     Pointeur de periode d'analyse
      if ( .not. allocated(begvp_oa)      ) allocate(begvp_oa (nzv_oa+1)      )   ! structure frequentielle du vecteur d etat


      return
      end subroutine allocate_part1_oa

      subroutine deallocate_part1_oa( )
      implicit none
      !deallocate()  
      return
      end subroutine deallocate_part1_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note Croco tile(-threads) : this routine must be called before any
!!       openMP parallel region fork and/or whitout loop over tiles
!!       so that the only MASTER thread initializes the implicitly 
!!       openMP SHARED module variables without any data race issue.
!
!
! DESCRIPTION: lecture des namelists specifiques à chaque configuration d'analyse
!              demandee par l'utilisateur dans la namelist principale namelist_oa 
!
!> @brief Read each specific namelist parameters associated to each requested OA analysis.
!!
!> @details Specific namelist associated to each configuration-variable-online-analysis, 
!! which contains the parameters describing the nature of the analysis 
!! (Fourier, Dirac, Wavelet atom,..), the horizontal and vertical domain to analyse, 
!! the temporal convolution window characteristics (temporal scales, frequency and sampling,...), 
!! ..., etc. The name of the specific namelists are constructed using the 
!! parameter nzc_oa_names as specified in the main namelist (namelist_oa): 
!! "namelist_" \+ nzc_oa_names \+ "_oa"
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists + Comments + Stand-alone version (2015)
!!  - Adapted to Croco parameters (grid, time stepping,...) (2020)
!!    - Add new Croco subroutine tool_datetosec to
!!      convert date YYYY, MM, DD, HH, MM, SS to the model time step index.
!> @date 2021
!> @todo BLXD 
!! - test new croco sub. tool_datetosec. with namelist swt_t_oa=4
!! - cleaning this old source would be nice...
!------------------------------------------------------------------------------

      subroutine rdnotebook_oa( )

      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
!     use module_oa_stock
!     use module_oa_upd
!     use module_nrj
#ifdef MPI
      use module_parameter_oa, only : iminmpi, jminmpi
#endif
      implicit none


      !> Local Namelist parameters
      integer :: TPSI_W,     &  !< Type of atom (0:Dirac, 1:Ondelette de Morlet, 2:Windowed Fourier, 3: Fourier)
                 SAVE_W,     &  !< Not applied (0,1:variable utilisateur,2:variable commune)
                 UPDV_W,     &  !< Variable remise à jour interactivement? Keep it to 2.
                 SWT_WFPF_W, &  !< Outputs: real part Wf (1), Energy Pf=|Wf|^2 (2) Absolute value |Wf| (4) Full complex analysis Wf
                 FL_REC_W,   &  !< Reconstruction flag 0: raw coefficient, 1: reconstruction
                 DORI_W,     &  !< Frequency analysis configuration: Discrete (1), Integration (2)
                 NZPT_PER_W, &  !< Number of points per period of analysis
                 CNB_W,      &  !< call number: position of the call in the baroclinic / barotropic time step. Not applied. Keep it to zero.
                 SWT_T_W,    &  !< 1: between two dates [T1:T2] , 2: [T1,end] , 3: single date T1, 4: all the simulation (#BLXD ONLY 4 TESTED)
                 SWT_D_W,    &  !< Domaine spatial (1: région [lat1,lat2],[lon1,lon2], 2: point (I,J), 3: Région [Hmin,Hmax] dans domaine (Lat,Lon)
                 DX_W,       &  ! Résolution en X (1 point sur...)
                 DY_W,       &  ! Résolution en Y (1 point sur...)
                 DK_W           ! Résolution en Z (1 point sur...)


      integer , dimension(1:2) :: K_W    !< Niveaux verticaux Min et Max K1,K2  (-99 pour colonne entiere)

      !> TODO test the routine to convert date YYYY,MM,DD,HH,MM,SS to the simulation time step index
      integer, dimension(1:6) :: DATE_DEB, & !< Annee, Mois, Jour, Heure, Min., Sec. --> ex : 2011,01,1,00,01,18
                                 DATE_FIN    !< Annee, Mois, Jour, Heure, Min., Sec. --> ex : 2011,01,1,00,15,15

      real  :: DT_W,  &         !< Temporal period in second at wich analysis are repeated (if SWT_T_W=1 ou 4)
               T0_W,  &         !< Time lag to launch the first analysis (SWT_T_W=4)
               DELTA_T_W        !< 1/2 Largeur de l'atome ou nombre de périodes par atome (1 for Dirac, 2 or more for Wavelet)


      real, dimension(1:2)  :: LAT_W,  & !< Latitude Min et Max
                               LON_W,  & !< Longitude Min et Max
                               HH_W,   & !< Profondeur Min et Max
                               PTIJ_W    !< Point particulier demande par l utilisateur

      real, dimension(1:3)  :: PER_W     !< Minimum, Maximum time period (s) (wavelets & fourier windows) or Length (Dirac brush), 0 for delta Dirac
                                         !! step determining # of periods to be examined (-99=optimum)

      !> Local 
      character(len=250)  :: filin !< namelist filename

      integer ::                                                      &
           iv_r                                                       &
          ,ip_r                                                       &
          ,ic_r                                                      

      real    :: dt_r

      double precision    :: time_from_croco_tref_in_sec
      integer             :: izv_oa

      namelist /oa_parameters/ TPSI_W, SAVE_W, UPDV_W, SWT_WFPF_W, FL_REC_W, PER_W,    &
                               DORI_W, DELTA_T_W, NZPT_PER_W, CNB_W, SWT_T_W, SWT_D_W, &
                               DATE_DEB, DATE_FIN, DT_W, T0_W, LAT_W, LON_W,           &
                               HH_W, K_W, DX_W, DY_W, DK_W, PTIJ_W
   
      oa_analysis_requested : if (nzc_oa /=0) then
    
!..... add unite_oa to namelist i) seconds unite_oa=1, ii) hours unite_oa=3600
      
      izv_oa = 0 ! Must be equal to nzv_oa at the end of the oa_config loop

      if(if_print_node) write(io_unit,*) '...Loop on Online Analysis configurations'
      oa_config_loop : do ic_r = 1 , nzc_oa
       
         if(if_print_node) write(io_unit,*) '...Configuration index ', ic_r

!..... namelist lue:
 
         filin = trim(directory_in_oa) // txtslash // 'namelist_' // trim( nzc_oa_names(ic_r) ) // '_oa'

         open(unit=io_nml, file=filin)
         read(unit=io_nml, nml  = oa_parameters)

!..... HYP : config. ic_r has a single variable ie, nzvc_oa( nzc_oa_vartyp(ic_r) ) = 1
         izv_oa = izv_oa + 1

!......lecture d'une configuration:
!
         tv_oa(izv_oa)      = nzc_oa_vartyp(ic_r)           ! type of variable (of config. izv_oa)
         tv_sclg_oa(izv_oa) = nzc_oa_scalogram(ic_r)        ! type of variable (of config. izv_oa)
         
         if(if_print_node) write(io_unit,*) '...Variable type/scalogram ', tv_oa(izv_oa),tv_sclg_oa(izv_oa)
  
!......mise a jour de la config.:

         tc_oa (ic_r)  = tv_oa (izv_oa)             ! type of variable(of config izv_oa)
         begc_oa(ic_r) = izv_oa                     ! config pointer : several variables per config

         tpsi_oa (izv_oa)    = TPSI_W          ! type of atom
         save_oa (izv_oa)    = SAVE_W      
         updv_oa (izv_oa)    = UPDV_W      
         swt_wfpf_oa(izv_oa) = SWT_WFPF_W
         fl_rec_oa(izv_oa)   = FL_REC_W
         per_oa (1,izv_oa)   = PER_W(1)       ! minimum period (h)
         per_oa (2,izv_oa)   = PER_W(2)       ! maximum period (h)
         per_oa (3,izv_oa)   = PER_W(3)       ! step determining # of periods to be examined : evenly spaced or optimal
         dori_oa(izv_oa)     = DORI_W
         delta_t_oa (izv_oa) = DELTA_T_W      ! # periods per atom
         nzpt_per_oa(izv_oa) = NZPT_PER_W     ! # calculation points per period
         cnb_oa(izv_oa)      = CNB_W          ! call number which may be used to define where a variable needs to be called.
!        pour les bilans energetiques cette variable est mise a jour automatiquement un peu plus loin.

         if ( per_oa(1,izv_oa).gt. per_oa(2,izv_oa) ) then
            write (*,*) " ERROR : per_oa (1,",izv_oa,") trop grand"
            if(if_print_node) write (io_unit,*) " ERROR : per_oa (1,",izv_oa,") trop grand"
            stop
         endif  

         if (per_oa (2,izv_oa).ne.per_oa (1,izv_oa).and.                        &
              per_oa (3,izv_oa).gt.(per_oa (2,izv_oa)-per_oa (1,izv_oa))        &
              .and.per_oa (3,izv_oa).ne.-99.) then                                                         
            write (*,*) " ERROR : per_oa (3,",izv_oa,") trop grand"
            if(if_print_node) write (io_unit,*) " ERROR : per_oa (3,",izv_oa,") trop grand"
            stop
         endif  
         
         do ip_r = 1 , 3
            per_oa  (ip_r,izv_oa) = per_oa(ip_r,izv_oa) * unite_oa 
         enddo

         swt_t_oa(izv_oa) = SWT_T_W
         swt_d_oa(izv_oa) = SWT_D_W

         if  ( tv_sclg_oa(izv_oa) .eqv. .true. ) then
             if ( dori_oa(izv_oa) /=1 ) then
                write (*,*) " ERROR OA : scalogram is possible with discrete periods only => DORI_W=1"
                if(if_print_node) write (io_unit,*) " ERROR OA : scalogram is possible with discrete periods only => DORI_W=1"
                stop
             else if ( swt_d_oa(izv_oa) /=2 ) then
                write (*,*) " ERROR OA : scalogram only possible for spatial config SWT_D_W = 2 (i,j) point"
                if(if_print_node) write (io_unit,*) " ERROR OA : scalogram only possible for spatial config SWT_D_W = 2 (i,j) point"
                stop
             else if ( (per_oa (2,izv_oa).eq.per_oa (1,izv_oa)) .or.                            &
                       ( per_oa (3,izv_oa).gt.(per_oa (2,izv_oa)-per_oa (1,izv_oa))  .and.      &
                         per_oa (3,izv_oa).ne.-99 ) ) then
                write (*,*) " ERROR OA : scalogram requires DP = PER_W(2)-PER_W(1)>0 and PER_W(3) = -99 or < DP"
                if(if_print_node) write (io_unit,*) " ERROR OA : scalogram requires DP = PER_W(2)-PER_W(1)>0 and PER_W(3) = -99 or < DP"
                stop
             endif
         endif

! #BLXD 2020 dev new croco sub. tool_datetosec
!            convert date YYYY, MM, DD, HH, MM, SS to the model time step index
!            TODO test with swt_t_oa=4 (even though not used ) then remove

         call tool_datetosec( DATE_DEB(1), DATE_DEB(2), DATE_DEB(3), &
          &                   DATE_DEB(4), DATE_DEB(5), DATE_DEB(6), &
          &                   time_from_croco_tref_in_sec)
         kount_user_oa(1,izv_oa) = int(time_from_croco_tref_in_sec/dti)
         
         call tool_datetosec( DATE_FIN(1), DATE_FIN(2), DATE_FIN(3), &
          &                   DATE_FIN(4), DATE_FIN(5), DATE_FIN(6), &
          &                   time_from_croco_tref_in_sec)
         kount_user_oa(2,izv_oa) = int(time_from_croco_tref_in_sec/dti)
         
!...>    (1) [T1,T2], (2) [T1,end], (3) single date T1             
         if  ( (swt_t_oa(izv_oa) == 1) .or. ( swt_t_oa(izv_oa) == 2)  .or. ( swt_t_oa(izv_oa) == 3) ) then 
         call tool_datetosec( DATE_DEB(1), DATE_DEB(2), DATE_DEB(3), &
          &                   DATE_DEB(4), DATE_DEB(5), DATE_DEB(6), &
          &                   time_from_croco_tref_in_sec)

             kount_user_oa(1,izv_oa) = int(time_from_croco_tref_in_sec/dti)
         end if
         if  (swt_t_oa(izv_oa) == 1) then 
         call tool_datetosec( DATE_FIN(1), DATE_FIN(2), DATE_FIN(3), &
          &                   DATE_FIN(4), DATE_FIN(5), DATE_FIN(6), &
          &                   time_from_croco_tref_in_sec)

             kount_user_oa(2,izv_oa) = int(time_from_croco_tref_in_sec/dti)
         end if
         ! #BLXD OnlineA 2020 change : nt_max-1 -> nt_max
         ! For wavelet and Fourier ( see var_time_oa )
         ! The nzvt_oa convolution windows are defined over the model time index 
         ! intervals [ kountv_oa(1,:), kountv_oa(2,:) ]
         ! Being given the width of the convolution window 2*dkount_tot_t
         ! and its time index center at k and being given the simulation time window :
         ! - the 1st possible convolution window can start at the 1st model time index
         ! which is the 'now' time index of the restart just before the 1st 
         ! model time step from now to after, ie kountv_oa(1, izvt_oa=1 ) = kount0
         ! - the last possible convolution window can have its last window time index
         ! at the last model time index, which is to say at 'after' state of the last
         ! - the last time index of the last possible convolution window si allowed to 
         !   match the 'after' index of the last time step of the model,
         !   kountv_oa( 2, izvt_oa=nzvt_oa ) = nt_max
         ! #BLXD TODO check with Dirac atom 
         !     a) if perv_oa = 0. => dkount_tot_t = 0 (OK) 
         !     b) else               dkount_tot_t =  INT(perv_oa / 2 / dti - 0.5) + 1

         !if  ( swt_t_oa(izv_oa) == 2) kount_user_oa(2,izv_oa) = nt_max-1
         if  ( swt_t_oa(izv_oa) == 2) kount_user_oa(2,izv_oa) = nt_max
         if  ( swt_t_oa(izv_oa) == 3) kount_user_oa(2,izv_oa) = kount_user_oa(1,izv_oa)
       
         if (swt_t_oa(izv_oa).eq.1.and.kount_user_oa(2,izv_oa).lt.kount_user_oa(1,izv_oa))   then
              if(if_print_node) write (io_unit,*) "probleme dates choisies t1 > t2" 
         end if
         !if (swt_t_oa(izv_oa).eq.1.and.kount_user_oa(1,izv_oa).gt.nt_max-1)                  then 
         if (swt_t_oa(izv_oa).eq.1.and.kount_user_oa(1,izv_oa).gt.nt_max)                  then 
              if(if_print_node) write (io_unit,*) "probleme dates choisies t1 > t final"
         end if
         if ( swt_t_oa(izv_oa) == 4) then
            kount_user_oa(1,izv_oa) = kount0
            ! BLD 2020 change kount_user_oa(2,izv_oa) = nt_max-1
            kount_user_oa(2,izv_oa) = nt_max
         endif

         if(if_print_node) write(io_unit,*) '=> SWT_T_OA option is ',swt_t_oa(izv_oa)
         if(if_print_node) write(io_unit,*) '   OA analysis over Time steps ',kount_user_oa(1,izv_oa),kount_user_oa(2,izv_oa)
         if(if_print_node) write(io_unit,*) '   Simu Time step interval including initial and last',kount0,nt_max
         if  ( (swt_t_oa(izv_oa) == 1) .or. ( swt_t_oa(izv_oa) == 2)  .or. ( swt_t_oa(izv_oa) == 3) ) then 
             write(*,*) 'ERROR OA : option not allowed (developments required), please set SWT_T_W to 4 in the OA namelist'  
             if(if_print_node) write(io_unit,*) 'ERRROR OA : option not allowed (to dev.), set SWT_T_W to 4 in the OA namelist'  
             stop
         end if                     
        
         dt_r = DT_W
         kount_user_oa(3,izv_oa) = max(1,int ( dt_r * unite_oa / dti ))

         if(if_print_node) write(io_unit,*) '   OA analysis time discretization ',kount_user_oa(3,izv_oa)
        
         t0_oa(izv_oa) = T0_W
         t0_oa(izv_oa) = t0_oa(izv_oa) * unite_oa
         kount_user_oa(1,izv_oa) = kount_user_oa(1,izv_oa) + int(t0_oa(izv_oa) / dti) 

         if(if_print_node) write(io_unit,*) '   OA analysis starting earlier with T0_W (sec) ',t0_oa(izv_oa)
         if(if_print_node) write(io_unit,*) '   OA analysis over Time steps ',kount_user_oa(1,izv_oa),kount_user_oa(2,izv_oa)

         lat_oa(1,izv_oa) = LAT_W(1)
         lat_oa(2,izv_oa) = LAT_W(2)
         lon_oa(1,izv_oa) = LON_W(1)
         lon_oa(2,izv_oa) = LON_W(2)
         h_oa(1,izv_oa) = HH_W(1)
         h_oa(2,izv_oa) = HH_W(2)

         k_oa(1,izv_oa) = K_W(1)
         k_oa(2,izv_oa) = K_W(2)


         dx_oa(izv_oa) = DX_W
         dy_oa(izv_oa) = DY_W
         dk_oa(izv_oa) = DK_W

#ifdef MPI
         ! BLXD MPI Bug
         ptij_oa(1,izv_oa)  =  PTIJ_W(1)  - iminmpi+1
         ptij_oa(2,izv_oa)  =  PTIJ_W(2)  - jminmpi+1
#else
         ptij_oa(1,izv_oa) = PTIJ_W(1)
         ptij_oa(2,izv_oa) = PTIJ_W(2)
#endif

         if ( ( swt_d_oa(izv_oa) == 2 ) .and. ( tv_sclg_oa(izv_oa) ) ) then
            if ( ( dk_oa(izv_oa) /= 1) .or. ( (k_oa(2,izv_oa)-k_oa(1,izv_oa))/=0 ) ) then       
                write (*,*) " ERROR OA : scalogram only possible with a single vertical grid point"
                if(if_print_node) then 
                write (io_unit,*) " ERROR OA : scalogram only possible with a single vertical grid point"
                write (io_unit,*) "            please set DK_W=1 and K_W=ktgt,ktgt where ktgt is your level target"
                endif
                stop
            endif
         endif

!  attention, test uniquement sur les points z et aps sur les points,x,y,rot
!  BLXD unit conv can be moved in init_oa
#ifdef SPHERICAL
         ! #BLXD 2020 change
         if ( if_sphe_deg2rad ) then
         if(if_print_node) write(io_unit,*) 'SPHERICAL (nml) : degrees converted to rads'

         lat_oa(1,izv_oa)=LAT_W(1)*pi_oa/180.D0
         lat_oa(2,izv_oa)=LAT_W(2)*pi_oa/180.D0     
         lon_oa(1,izv_oa)=LON_W(1)*pi_oa/180.D0
         lon_oa(2,izv_oa)=LON_W(2)*pi_oa/180.D0
         else 
         if(if_print_node) write(io_unit,*) 'SPHERICAL (nml) : degrees NOT converted to rads'
         endif
#else    
         
         if(if_print_node) write(io_unit,*) 'NOT SPHERICAL (nml) : coordinates xr, xp are in meters'
         
         lat_oa(1,izv_oa) = LAT_W(1)
         lat_oa(2,izv_oa) = LAT_W(2)
         lon_oa(1,izv_oa) = LON_W(1)
         lon_oa(2,izv_oa) = LON_W(2)
#endif

!.....NHOMS configuration treatment and nrj eliminated

!.....treatment of specific variable configuration

!.....TODO ? Move this in extend_domain. Not sure
 
         if ( tv_sclg_oa (izv_oa) ) then
            scalogram_analysis = .true.
            if(if_print_node) write (io_unit,*) "....Scalogram analysis requested"
         endif

         if ( tv_oa (izv_oa)==20 ) then
            isopycne_analysis = .true.
            if(if_print_node) write (io_unit,*) "....Isopycne analysis requested"
         else if ( (tv_oa (izv_oa)==99) .or. (tv_oa (izv_oa)==98) ) then
            test_analysis = .true.
            if(if_print_node) write (io_unit,*) "....Test analysis requested"
         endif

!.....treatment of "mixed" variable configuration
!.....HYP : config. ic_r has 1 variable ie, nzvc_oa( nzc_oa_vartyp(ic_r) ) = 1
         !if ( nzvc_oa( nzc_oa_vartyp(ic_r) ) >= 2 ) then
         if ( nzvc_oa( tv_oa(izv_oa) ) >= 2 ) then
            if(if_print_node) write(io_unit,*) 'Composite configuration : init_configuration_oa'
            call init_configuration_oa(izv_oa)
         endif

      enddo oa_config_loop

      if ( izv_oa /= nzv_oa ) then
        write (*,*) "ERROR in rdnotebook_oa : izv_oa should be equal to nzv_oa"
        if(if_print_node) write (io_unit,*) "ERROR in rdnotebook_oa : izv_oa should be equal to nzv_oa"
        stop
      endif

!.....desallocation tableau de lecture des namelists OA

      ! #BLXD DO NOT DEALLOCATE AND MAKE PUBLIC 
      ! nzc_oa_names (config name in namelist) => tv_nam_oa ?
      ! TODO useful to construct OA analysis variable field_def.xml ?
      ! write( vnam_oa, fmt='(a5,5,i3.3,a1,i3.3)') tgvnam_oa(tv_oa(iv)),'3d_r_',nzc_oa_names(ic),'_',tgvnam_oa(nzc_oa_vartyp(tv_oa(iv)))

      call deallocate_namelist_oa

!.....fin de la derniere configuration:

      begc_oa(nzc_oa + 1) = nzv_oa + 1

      close(io_nml)

!.....quelques controles:

      do ic_r=1,nzc_oa
         do iv_r=begc_oa(ic_r),begc_oa(ic_r+1)-1
            ! Fourier => analysis on the full simulation
            ! #BLXD why not choosing swt_t_oa(iv_r) = 1 or 2 ????
            if (tpsi_oa(iv_r).eq.3) then 
                if(if_print_node) write (io_unit,*) "WARNING : Fourier requested => changing SWT_T_OA to 4"
                swt_t_oa(iv_r) = 4
            endif
         enddo
      enddo

      endif oa_analysis_requested

      return

      end subroutine rdnotebook_oa


!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note WARNING user options hardcoded in this routine 
!
! DESCRIPTION: Example pour la construction d'analyse de variables composites. 
!
!> @brief Example of "mixed" variable configuration
!
!> @details Mixed variable configuration are configuration requiring several 
!! variables possibly available at different location of the calling code
!! (e.g., time-splitting). Note that :
!! - Configuration code from 1 to 99 corresponds to "simple" variable.
!! - Configuration code from 100 to XX corresponds to "mixed" variable configuration. 
!!   e.g., Kinetic energy calculated from zonal and meridional velocities.
!! See also the cnb_oa parameter set to zero for all variables by default. 
!! cnb_oa can be used to call main_oa at different code location 
!! In such case, change the cnb_oa parameter for particular variable, and adapt main_oa accordingly.
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Croco version (2020), Stand alone version (2015)
!> @date 2021
!> @todo BLXD var_oa_copy moved to module_oa_variable ? 
!------------------------------------------------------------------------------
      subroutine init_configuration_oa(izv_oa)

      use module_oa_variables, only : tv_oa

      implicit none

      integer, intent(inout) :: izv_oa ! counts for the final nzv_oa
      integer :: tmp_izv_oa

!.....configuration a multiples variables ie, nzvc_oa( nzc_oa_vartyp(ic_r) ) > 1

!-------------------------------------------
!     Example of mixed variable configuration (100)
!-------------------------------------------
!......> Composite variable with type index 100
!        Adding the variables, copying the config. and incrementing the total number of variable izv_oa  
         if (tv_oa(izv_oa).eq.100) then

            tmp_izv_oa = izv_oa

!------>1ere variable: velocity vel_u 
            tv_oa(izv_oa) = 1
            izv_oa = izv_oa + 1

!     copie toutes les infos pour la nouvelle variable:
            !call var_copy_oa ( izv_oa - 1 , izv_oa )
            call var_copy_oa ( tmp_izv_oa , izv_oa )

!------>2eme variable: velocity vel_v
            tv_oa(izv_oa) = 2
            izv_oa = izv_oa + 1

!     copie toutes les infos pour la nouvelle variable:
            !call var_copy_oa ( izv_oa - 2 , izv_oa )
            call var_copy_oa ( tmp_izv_oa , izv_oa )

!------>3eme variable: velocity vel_w
            tv_oa(izv_oa) = 3
            izv_oa = izv_oa + 1

!     copie toutes les infos pour la nouvelle variable:
            !call var_copy_oa ( izv_oa - 2 , izv_oa )
            call var_copy_oa ( tmp_izv_oa , izv_oa )

         else if (tv_oa(izv_oa)>100) then
            write (*,*) "ERROR OA : configuration above 100 is not yet hardcoded"
            if(if_print_node) write (io_unit,*) "ERROR OA : configuration above 100 is not yet hardcoded"
            stop
         endif

! Uncomment these lines and adapt them if you wish to call the OA treatment at specific
! point of the code (cf time splitting), for specific variables that are updated
! at some irregular time step or code place. 
!-------------------------------------------
!   Eventuellement point de sortie pour les variables
!      associees aux bilans energetiques
!-------------------------------------------
         ! if (tv_oa(izv_oa).eq.30) then
         !    outw_oa(5,2) =1
         !    outw_oa(5,5) =1
         !    cnb_oa(izv_oa)=-30
         !    des_oa(izv_oa)=-5005
         ! endif

      end subroutine init_configuration_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
! DESCRIPTION: 
!
!> @brief Handles "mixed" variables configuration. 
!
!> @details Configuration codes for mixed variables involves several variables,
!  and requires to copy out configuration informations to the variable parameters. 
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Scalogram (2021), Croco version (2020), stand alone version (2015)
!> @date 2021
!> @todo BLXD declare var_oa_copy into module_oa_variable ? 
!------------------------------------------------------------------------------

      subroutine var_copy_oa (                                       &
                  nzvold_c                                           &             
                 ,nzvnew_c   )

      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      !use module_oa_stock

      implicit none
      integer                                                         &
        nzvold_c                                                      &
       ,nzvnew_c

        tv_oa  (nzvnew_c)         = tv_oa  (nzvold_c)
        t0_oa  (nzvnew_c)         = t0_oa  (nzvold_c)
        tpsi_oa(nzvnew_c)         = tpsi_oa(nzvold_c)      
        save_oa(nzvnew_c)         = save_oa(nzvold_c)
        updv_oa(nzvnew_c)         = updv_oa(nzvold_c)
        swt_wfpf_oa(nzvnew_c)     = swt_wfpf_oa(nzvold_c)
        fl_rec_oa(nzvnew_c)       = fl_rec_oa(nzvold_c)
        per_oa (1,nzvnew_c)       = per_oa (1,nzvold_c)
        per_oa (2,nzvnew_c)       = per_oa (2,nzvold_c)
        per_oa (3,nzvnew_c)       = per_oa (3,nzvold_c) 
        delta_t_oa(nzvnew_c)      = delta_t_oa(nzvold_c)
        nzpt_per_oa(nzvnew_c)     = nzpt_per_oa(nzvold_c)
        cnb_oa(nzvnew_c)          = cnb_oa(nzvold_c)
        dori_oa(nzvnew_c)         = dori_oa(nzvold_c)
        ! BLXD croco2021 added
        des_oa(nzvnew_c)          = des_oa(nzvold_c)
        swt_t_oa(nzvnew_c)        = swt_t_oa(nzvold_c)
        swt_d_oa(nzvnew_c)        = swt_d_oa(nzvold_c)
        kount_user_oa(1,nzvnew_c) = kount_user_oa(1,nzvold_c)
        kount_user_oa(2,nzvnew_c) = kount_user_oa(2,nzvold_c)
        kount_user_oa(3,nzvnew_c) = kount_user_oa(3,nzvold_c)
        lat_oa(1,nzvnew_c)        = lat_oa(1,nzvold_c)
        lat_oa(2,nzvnew_c)        = lat_oa(2,nzvold_c)
        lon_oa(1,nzvnew_c)        = lon_oa(1,nzvold_c)
        lon_oa(2,nzvnew_c)        = lon_oa(2,nzvold_c)
        h_oa(1,nzvnew_c)          = h_oa(1,nzvold_c) 
        h_oa(2,nzvnew_c)          = h_oa(2,nzvold_c)
        k_oa(1,nzvnew_c)          = k_oa(1,nzvold_c)
        k_oa(2,nzvnew_c)          = k_oa(2,nzvold_c)
        dx_oa(nzvnew_c)           = dx_oa(nzvold_c)
        dy_oa(nzvnew_c)           = dy_oa(nzvold_c)
        dk_oa(nzvnew_c)           = dk_oa(nzvold_c)
        ptij_oa(1,nzvnew_c)       = ptij_oa(1,nzvold_c)
        ptij_oa(2,nzvnew_c)       = ptij_oa(2,nzvold_c)
        tv_sclg_oa  (nzvnew_c)    = tv_sclg_oa  (nzvold_c)
       return
       end subroutine var_copy_oa


!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note WARNING user options here hardcoded
!
! DESCRIPTION: a completer en fonction des variables a traiter  
!
!> @brief Specifies the number of dimensions of the variable to analyse, 
!! and at which grid point fields are defined.
!
!> @details This routine must be modified when adding a new variable:
!! - tgv3d_oa must be set to 3 or 2 for three and two-dimensional variable respectively.
!! - tgv_oa must be set to:
!!   - "1" for t-grid point,
!!   - "2" for u-grid point,
!!   - "3" for v-grid point,
!!   - "4" for f-grid point.
! 
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Croco version (2020)
!> @date 2021
!> @todo BLXD use croco grid definition ?
!------------------------------------------------------------------------------
      subroutine var_grid_oa
    
      !use module_oa_space , only : tgv3d_oa, tgv_oa, tgvnam_oa 

      implicit none

!.....vitesse vel_u:
      tgv_oa(1)   = 2 
      tgv3d_oa(1) = 3
      tgvnam_oa(1) = 'del_u'

!.....vitesse vel_v:
      tgv_oa(2)   = 3 
      tgv3d_oa(2) = 3
      tgvnam_oa(2) = 'del_v'

!     BLXD croco should avoid k=0  for w-grid point
!     2D rmask is appropriate for w-grid point
!.....vitesse vel_w:
      tgv_oa(3)   = 1 
      tgv3d_oa(3) = 3
      tgvnam_oa(3) = 'vel_w'

!.....vitesse vel_u:
      tgv_oa(4)   = 2 
      tgv3d_oa(4) = 3
      tgvnam_oa(4) = 'vel_u'

!.....vitesse vel_v:
      tgv_oa(5)   = 3 
      tgv3d_oa(5) = 3
      tgvnam_oa(5) = 'vel_v'

!.....vitesse u_bar:
      tgv_oa(6)   = 2 
      tgv3d_oa(6) = 2
      tgvnam_oa(6) = 'ubar_'

!.....vitesse v_bar:
      tgv_oa(7)   = 3 
      tgv3d_oa(7) = 2
      tgvnam_oa(7) = 'vbar_'

!.....variable temp:
      tgv_oa(8)   = 1 
      tgv3d_oa(8) = 3
      tgvnam_oa(8) = 'temp_'

!.....variable salt:
      tgv_oa(9)   = 1 
      tgv3d_oa(9) = 3
      tgvnam_oa(9) = 'salt_'

!.....variable density:
      tgv_oa(11)   = 1 
      tgv3d_oa(11) = 3
      tgvnam_oa(11) = 'rho__'

!.....variable rhp_t to track isopycne movement:
      tgv_oa(20)   = 1 
      tgv3d_oa(20) = 3
      tgvnam_oa(20) = 'isplv'

!.....variable scalogram ubar:
      tgv_oa(56)   = 1 
      tgv3d_oa(56) = 2
      tgvnam_oa(56) = 'scal_'

!.....variable scalogram density:
      tgv_oa(61)   = 1 
      tgv3d_oa(61) = 3
      tgvnam_oa(61) = 'scal_'

!.....test variable 2D vardp_test_oa :
      tgv_oa(97)   = 1 
      tgv3d_oa(97) = 2
      tgvnam_oa(97) = 'tst2_'

!.....test variable vardp_test_oa pour scalogram:
      tgv_oa(98)   = 1 
      tgv3d_oa(98) = 3
      tgvnam_oa(98) = 'scal_'

!.....test variable 3D vardp_test_oa:
      tgv_oa(99)   = 1 
      tgv3d_oa(99) = 3
      tgvnam_oa(99) = 'test_'

!.....horizontal kinetic energy at t-grid point:
      tgv_oa(100)   = 1 
      tgv3d_oa(100) = 3
      tgvnam_oa(100) = 'comp_'

      return
      end subroutine var_grid_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note This routine must be called under a Croco tile(-thread) loop 
!!       (see main.F).
!
! DESCRIPTION:  initilisation des structures spatiales, frequentielles et 
!!  temporelles du vecteur d'etat.
!
!> @brief construction of the spatial, temporal and scale structures needed
!!  to stack/unstack the analysis into/from the state vector. 
!
!> @details The OA state vector is a derived data type array whose dimension
!!   is the # of temporal convolution windows (time structure) 
!!   for all the user-requested configurations and variables and 
!!   whose field is the # of grid point to analyse (spatial structure). 
!!   OA initialization defines several arrays to vectorize the analysis.
!!   The spatial structure is Croco tile(-thread) dependent.
!!  
!! 1) Croco tile-thread version : OA module variables which are related to horizontal
!!    subdomain are now embeded in a derived type array having the dimension of the
!!    number of tiles. This routine is called within a tile loop in main.F.
!!    see, init_online_spectral_diags(Istr,Iend,Jstr,Jend,tile) 
!! 2) OA stand-alone/croco version WARNING : 
!!    Memory duplication when passing arguments to the init_oa subroutine :
!!    - the calling program should set the exact array dimension ranges 
!!      using lbound and ubound fonctions.   
!!    - for croco tile-thread version => use pointers (applied for grid arrays only)
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists + Stand-alone version + optimization (2015)
!!  - Headers, comments, cleaning, from f77 to earlier standards  
!!  - Croco-OA interface 1st version (2020)
!!  - XIOS outputs (2021)
!!  - Toward a Croco-OA tile(-threads) compliant interface and OA version (2021)
!!    supporting tile(-threads) and possible dual OpenMP-MPI parallelization 
!!    of the horizontal domain :
!!    - Variables with spatial dependance are tile(-thread) dependent and
!!      handled with an array of derived type which depends on tiles.
!!    - Croco tile(-thread) loop => array of derived type :
!!       - st(tile)% defined in module_tile_oa for perennial use over
!!         simulation
!!       - grd(tile)% defined in module_grd_oa for local temporary use 
!!         of the Croco grid parameters passed to OA as args.
!!    - openMP CRITICAL construct + last tile execution
!!  - Developments to perform online Scalogram analyses (2021)
!!    - from namelist options to calculation and XIOS netcdf outputs
!!    - scalograms distribution over MPI process/subdomains handled 
!!      either by XIOS or by MPI internal OA instructions (see if_mpi_oa).
!!    - Implementation of the Heisenberg theoritical and numerical uncertainties
!!      Ongoing. 
!! @todo BLXD
!! - The CROCO-OnlineAnalysis module interface must be modified : 
!!   Interface = mix between the Stand-alone and Original OnlineA versions 
!!   => Croco arrays are passed as arguments to the OnlineA routines with reduced 
!!   dimension range. It leads to memory duplication.
!! - specific module to set consitent precisions with real(KIND=wp)... etc 
!! - reintroduce pointers for isopycne tracking and/or function var_oa ?
!! - test/correct Croco Tile-thread compliant version with a 
!!    Croco tile-thread working vesion
!! - non XIOS outputs, handling temporal convolution windows and restarts. 
!------------------------------------------------------------------------------
!
      subroutine init_oa( tile          &
      ,iic_oa                           &
      ,imin, imax                       &
      ,jmin, jmax                       &
      ,kmin, kmax                       &
      ,lon_t_oa                         & !(-1:imax+2,-1:jmax+2)            
      ,lat_t_oa                         & !(-1:imax+2,-1:jmax+2)            
      ,lon_u_oa                         & !(0:imax+2,0:jmax+2)              
      ,lat_u_oa                         & !(0:imax+2,0:jmax+2)              
      ,lon_v_oa                         & !(0:imax+2,0:jmax+2)              
      ,lat_v_oa                         & !(0:imax+2,0:jmax+2)              
      ,lon_f_oa                         & !(0:imax+2,0:jmax+2)              
      ,lat_f_oa                         & !(0:imax+2,0:jmax+2)              
      ,mask_t_oa                        & !(-1:imax+2,-1:jmax+2,0:kmax+1)
      ,mask_f_oa                        & !(0:imax+1,0:jmax+1,0:kmax+1)
      ,mask_u_oa                        & !(0:imax+1,0:jmax+1,0:kmax+1)
      ,mask_v_oa                        & !(0:imax+1,0:jmax+1,0:kmax+1)
      ,h_w_oa                           & !(0:imax+1,0:jmax+1)              
      ,h_u_oa                           & !(0:imax+1,0:jmax+1)              
      ,h_v_oa                           & !(0:imax+1,0:jmax+1)              
      ,h_f_oa                           & !(0:imax+1,0:jmax+1)              
      ,rhp_t                            & 
      ,depth_t                          & !(-1:imax+2,-1:jmax+2,0:kmax+1) 
      ) 

      use module_grd_oa, only  : associate_grid_oa_ptr, nullify_grid_oa_ptr 
      use module_tile_oa, only : allocate_begvs_oa, allocate_begvs3d_oa       &
                                ,allocate_wf_oa, allocate_tile_sclg_oa        &
                                ,allocate_tile_test_oa ,ntiles
      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
!     use module_oa_stock

      implicit none

      !> Current model integration iteration
      integer, intent(in) :: iic_oa 

      !> Grid index range for the analysis of fields passed in argument to the OA module 
      !! Grid index are Roms-Croco tile(-thread) indices.
      !! If no tile(-thread) grid index are MPI subdomain index
      !! If no tile(-thread) and no MPI subdomain decomposition => full domain grid range
      integer, intent(in) ::       &
             imin, imax            & 
            ,jmin, jmax            &
            ,kmin, kmax

      !> tile parameter from croco
      integer, intent(in) :: tile

      !> 3-Dimensional density array at t-grid point
      double precision, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in), target :: rhp_t
                                                                          
      !> 3-Dimensional depth array at t-grid point
      double precision, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in), target :: depth_t

      ! 2-Dimensional array coordinate arrays passed from the calling code
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

      !> Lower and upper index bounds of the 3-Dimensional density array (see if pst%imin,...)
      integer, dimension(3)  ::  rhp_t_lbound, rhp_t_ubound

!.....History file:
!$OMP MASTER
#ifdef OA_TRACES
      if(if_print_node) write(io_unit,*) '...OA history file initialization'
#endif
      call history_oa(4,-1,-1,-1,-1, -1)
!$OMP END MASTER

!.....initialisation de la structure spatiale:
!
!     BLXD croco2021 tile-thread : 
!     - s(tile)%begvs_oa   : can be alloc. bef. the 1st var_space_oa call (nzv_oa size)
!     - s(tile)%begvs3d_oa : cannot be allocated. bef. size calculated in var_space_oa

!$OMP MASTER
#ifdef OA_TRACES
      if(if_print_node) write(io_unit,*) '...ALLOCATING allocate_begvs_oa size(nzv_oa+1)'
#endif
!$OMP END MASTER

      call allocate_begvs_oa(imin,imax,jmin,jmax,kmin,kmax,tile) 

!$OMP MASTER
#ifdef OA_TRACES
      if(if_print_node) write(io_unit,*) '...TEMPORARY Grid pointer association'
#endif
!$OMP END MASTER

      call associate_grid_oa_ptr( tile                    &
                        ,imin, imax                       &
                        ,jmin, jmax                       &
                        ,kmin, kmax                       &
                        ,lon_t_oa                         &
                        ,lat_t_oa                         &
                        ,lon_u_oa                         &
                        ,lat_u_oa                         &
                        ,lon_v_oa                         &
                        ,lat_v_oa                         &
                        ,lon_f_oa                         &
                        ,lat_f_oa                         &
                        ,mask_t_oa                        &
                        ,mask_f_oa                        &
                        ,mask_u_oa                        &
                        ,mask_v_oa                        &
                        ,h_w_oa                           &
                        ,h_u_oa                           &
                        ,h_v_oa                           &
                        ,h_f_oa                           )

      if (scalogram_analysis) then
        call allocate_tile_sclg_oa( tile, nzv_oa, nper_sclg_max )
      endif

!$OMP MASTER
#ifdef OA_TRACES
      if(if_print_node) write(io_unit,*) '...OA ENTERING var_space_oa .false. : calc. space dim'
#endif
!$OMP END MASTER

      call var_space_oa(.false., tile ) 

!$OMP MASTER
#ifdef OA_TRACES
      if(if_print_node) write(io_unit,*) '...OA ALLOCATING allocate_begvs3d_oa size(nzvs_oa+1)'
#endif
!$OMP END MASTER

      call allocate_begvs3d_oa(tile) 

!$OMP MASTER
#ifdef OA_TRACES
      if(if_print_node) write(io_unit,*) '...OA ENTERING var_space_oa .true. : set spatial variables'
#endif
!$OMP END MASTER

      call var_space_oa(.true., tile )

!.....BLXD horizontal grid information no longer needed : nullifying pointers 
       
      call nullify_grid_oa_ptr( tile ) 

! BLXD MPI subdomain without any work DO WE CONTINUE ?
! ic, iv, lt <-> lp,.. but nzvs_oa=0 and nzvs3d_oa=0 

!.....initialisation de la structure temporelle:
!$OMP CRITICAL (var_time_oa)

!      BLXD module variables (implicitly) SHARED for openMP
!      if tile-thread (CRITICAL REGION) or if sequential tiles
!      => all the tile will do the tile_count instruction one after the other

          tile_count_oa=tile_count_oa+1

!      The last tile will do the MPI receiving job

      last_tile : if ( tile_count_oa .eq. ntiles) then 

          call var_time_oa(.false., kount0, nt_max, dti )

          call allocate_part4_oa

          call var_time_oa(.true., kount0, nt_max, dti )

          tile_count_oa=0

      endif last_tile

!$OMP END CRITICAL (var_time_oa)

!.....initialisation boite d'heisenberg (module_oa_periodes):

      heisen : if (if_heisenberg_box) then

!---->Fichier de sortie:
!$OMP MASTER
        file_hbg = trim(directory_out_oa) // txtslash // 'heisenberg_oa.dat'

        if (mynode==0) open(unit=io_hbg,file=trim(file_hbg))   
!$OMP END MASTER

! implicit openMP barrier after critical since alloc/init_box_oa needed in var_rec_oa
!$OMP CRITICAL (heisenberg_init)

!      BLXD module variables (implicitly) SHARED for openMP
!      if tile-thread (CRITICAL REGION) or if sequential tiles
!      => all the tile will do the tile_count instruction one after the other

          tile_count_oa=tile_count_oa+1

!      The last tile will do the MPI receiving job

      last_tile0 : if ( tile_count_oa .eq. ntiles) then 

        call allocate_theo_box_oa
        
        call allocate_box_oa   
        
        call init_temporal_box_oa( nzvp_oa, psi_norm_l2, t0, sq_dt0, dt0 ) 

          tile_count_oa=0

      endif last_tile0

!$OMP END CRITICAL (heisenberg_init)

      endif heisen

!.....initialisation du facteur de reconstruction
!     BLXD heisenberg step 1

!$OMP CRITICAL (var_rec_nmsimult_oa)

!      BLXD module variables (implicitly) SHARED for openMP
!      if tile-thread (CRITICAL REGION) or if sequential tiles
!      => all the tile will do the tile_count instruction one after the other

          tile_count_oa=tile_count_oa+1

!      The last tile will do the MPI receiving job

      last_tile2 : if ( tile_count_oa .eq. ntiles) then 

      call var_rec_oa( dti )

      call unset_convwind_oa()

!.....BLXD set the max. simulation size of the module_oa_type array structure wf_oa
!     The nmsimult_oa size is a module variable (implicitely openMP SHARED)
!     Its value isn t tile(-threads) 

          if (nmsimult_oa_max==-99) then

            call count_nmsimult_oa( 0                          & ! ichoix
                                   ,1                          & ! ivar_m
                                   ,kount0                     &
                                   ,nt_max                     &
                                   ,dti                        &
                                   ,tile                       ) 
          else
      
            call user_count_nmsimult_oa( ) 

          endif
 
!         Pay attention that pst%nzw_oa is now set to nmsimult_oa in count_nmsimult_oa
!         user settings may underestimate the size required for wf_oa array
!         If simulation does not require restart (not yet DEV for OA)
!         To let the OA code calculate the required size set the nmsimult_oa_max to -99
!         in the oa_params namelist block (namelist_oa) 
         
          tile_count_oa=0

      endif last_tile2

! Implicit barrier with END CRITICAL region

!$OMP END CRITICAL (var_rec_nmsimult_oa)

!     BARRIER needed here since nmsimult_oa needed by all tile-threads 

      call allocate_wf_oa(tile,nmsimult_oa) 

!$OMP MASTER

      call unset_convwind_oa()

      if (if_heisenberg_box) then
        if (mynode==0) close(io_hbg)
      endif

!.....History file:        

      call history_oa(6,-1,-1,-1,-1, tile)

  
!.....ecriture du fichier de structure:

      !STDALONE call struct_oa

!.....History file:

      call history_oa(5,-1,-1,-1,-1, tile, iic_oa, nt_max )

!$OMP END MASTER

!.....initialisation eventuelle des variables de test
!     BLXD must be performed by all tile(-threads) in any order, concurrently or not 
      if (test_analysis) then
         call allocate_tile_test_oa( tile           &
                          ,imin=imin, imax=imax     &
                          ,jmin=jmin, jmax=jmax     &
                          ,kmin=kmin, kmax=kmax )
      endif

!.....initialisation eventuelle des "levels"
!     BLXD must be performed by all tile(-threads) in any order, concurrently or not 
      if (isopycne_analysis) then
          
        call allocate_lev_part1_oa(tile)

        rhp_t_lbound = (/imin,jmin,kmin/)
        rhp_t_ubound = (/imax,jmax,kmax/)

        call lev_init_oa( tile, rhp_t, rhp_t_lbound, rhp_t_ubound ) 

      endif

      ! OLD symphonie initialisations
      !......initialisation de variables symphonie:
      !    call var_upd_oa (-1,-1,-1)
      !......initialisation de l'energie
      !     call nrj_init_oa

!$OMP MASTER
      ! Output arrays var?d allocation using pre-calculated dimension nzupd?d_oa
      ! including the output scalogram option if_mpi_oa=T (mpi not handled by XIOS) 
      call upd_init_oa(.true.)
!$OMP END MASTER
! openMP barrier called few lines later if no CRITICAL region

      scalogram_init : if (scalogram_analysis) then 

!$OMP CRITICAL (init_tile_count_oa)
          ! CRITICAL => implicit BARRIER
          ! if tile-thread (CRITICAL REGION) or if sequential tiles
          ! => all tile will do this instruction one after the other
          tile_count_oa=tile_count_oa+1
          ! Here only the last tile will do the job
          last_tile3 : if ( tile_count_oa .eq. ntiles) then 
             call init_scalogram_oa()
           tile_count_oa=0
          end if last_tile3
!$OMP END CRITICAL (init_tile_count_oa)

      call sclg_init_coords(tile)

      else scalogram_init

!$OMP BARRIER

      endif scalogram_init

#ifdef OA_TEST_MPI_XIOS
      if ( if_chks) then
!$OMP SINGLE
        if ( allocated(cntvt_cr) ) deallocate(cntvt_cr)
!$OMP SINGLE NOWAIT
      endif
#ifdef MPI
      !if (if_chks) then
      !    ! BLXD TODO REMOVE these tests and flags
      !    if (if_test_mpi_blocking_send) call test_mpi_blocking_send() 
      !    if (if_test_test_mpi_nonblocking_ssend) call test_mpi_nonblocking_ssend(.true.)
      !endif
#endif
#endif

! BLXD IMPT barrier instruction : see init_online_spectral_diags_tile
!        call deallocate_tile_grid_tmp_str()
! DEALLOCATES THE GRID STRUCTURE, POINTERS HAVE ALREADY BEEN NULLIFIED
! C$OMP BARRIER

!$OMP MASTER
!$OMP END MASTER

      return
      end subroutine init_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note Croco openMP tile(-thread), this routine is called under a tile(-thread)
!!       loop. In the case of openMP threads, it must be called inside a 
!!       CRITICAL region (implicit BARRIER) by the last tile(-thread) only.
!
! DESCRIPTION:  dans le cas d'analyses en scalogram distribuees sur different
!               process MPI, initilisation des parametres permettant les sorties
!               XIOS et/ou les communications MPI realisees par le module OA. 
!
!> @brief in the cases of scalogram analysis distributed among different MPI
!!        process, parameters and variables initialization to ensure the XIOS
!!        distributed outputs and/or MPI communications within the OA module  
!!     
! REVISION HISTORY:
!
!> @authors 
!! B. Lemieux-Dudon
!> @date june-july 2021
!! @todo
!------------------------------------------------------------------------------
!
      subroutine init_scalogram_oa ()

      use module_oa_variables, only : nzv_oa, tupd_oa, tv_sclg_oa
      use module_oa_periode,  only  : nper_sclg_max, nzc_oa, begc_oa, tvc_oa, tc_oa
      use module_tile_oa, only      : scl, ntiles
      use module_oa_time, only      : nmsimult_oa

      implicit none

      integer                 :: iv, ic, itile
      integer                 :: isclg_loc, isclg_glo

#ifdef MPI
      integer, dimension(1:2) :: buff_s, buff_r
      logical, parameter :: snb=.false.
      integer :: tag, request, ierr
#endif

      if ( .not. allocated(v2locsclg) ) then
        allocate( v2locsclg(1:nzv_oa))
        v2locsclg(1:nzv_oa) = -99 
      end if
#ifdef MPI
#endif
#ifndef MPI
      if ( .not. allocated(index_s_cr) ) then
        if (nzupd0d_oa>0) then
        allocate( index_s_cr(0:nzupd0d_oa-1) )
        else
        allocate( index_s_cr(0) )
        endif
      end if
#else
      !if (if_mpi_oa) then
          if ( .not. allocated(index_s_cr) ) then
            if (nzupd0d_oa>0) then
            allocate( index_s_cr(0:nzupd0d_oa-1) )
            else
            allocate( index_s_cr(0) )
            endif
          end if
      !endif
#endif

      ! Counts the number of scalograms per MPI subprocess if any
       nsclg_loc = 0 ! nzupd0d_oa global number of scalogram
      ! Gives the correspondance from iv to the MPI 
      ! subdomain scalogram index v2locsclg(iv)
       isclg_loc = 0
      ! BLXD TODO The itile(-thread) loop order should be  controlled here
       do itile=1,ntiles
        nsclg_loc = nsclg_loc + scl(itile)%nsclg
        do iv = 1,nzv_oa
          ! it is a scalogram analysis
          if ( tv_sclg_oa(iv) ) then
            ! the iv scalogram analysis is included on 
            ! the itile region and therefore on the MPI process
            if ( scl(itile)%v2sclg(iv).ne.-99 ) then
               isclg_loc = isclg_loc + 1
               v2locsclg(iv) = isclg_loc
            endif
          endif 
        enddo
       enddo 
      if ( isclg_loc .ne. nsclg_loc ) then
        write(*,*)'ERROR OA : counts of MPI scalogram mismatch'
        if(if_print_node) write(io_unit,*)'ERROR OA : counts of MPI scalogram mismatch'
        stop
      else
      endif 
#ifdef MPI
      ! Only MPI node zero will know the global scal
      if ( .not. allocated( scal0d_oa) ) then
        allocate( scal0d_oa(nper_sclg_max,nsclg_loc) )
        scal0d_oa(:,:)=(0.D0,0.D0)
      endif 
      if ( .not. allocated( iscal0d_oa) ) then
        allocate( iscal0d_oa(nsclg_loc) )
      endif 
      if ( .not. allocated( jscal0d_oa) ) then
        allocate( jscal0d_oa(nsclg_loc) )
      endif 
      if ( .not. allocated( kscal0d_oa) ) then
        allocate( kscal0d_oa(nsclg_loc) )
      endif 
      if ( .not. allocated( per0d_oa) ) then
        allocate( per0d_oa(nper_sclg_max,nsclg_loc) )
        per0d_oa(:,:)=0.D0
      endif
      if ( .not. allocated( index_s_oa ) ) then
        if (nsclg_loc>0) then
        allocate( index_s_oa(0:nsclg_loc-1) )
        else
        allocate( index_s_oa(0) )
        endif
      endif 
      if ( .not. allocated( buff2r_s) ) then
        allocate( buff2r_s(nper_sclg_max) )
      endif
      if ( .not. allocated( buff2r_r) ) then
        allocate( buff2r_r(nper_sclg_max) )
      endif
#endif 
      isclg_loc  = 0
      isclg_glo  = 0
      do iv = 1,nzv_oa
         ! it is a scalogram analysis
         if ( tv_sclg_oa(iv) ) then
              isclg_glo = isclg_glo + 1
              !sclg2v(isclg_glo) = iv
#ifdef MPI
            if ( v2locsclg(iv).ne.-99 ) then
               ! Local mpi process scalogram index_s_oa only allocated ifMPI
               isclg_loc = isclg_loc + 1
               index_s_oa(isclg_loc-1) =  tupd_oa(iv)
            endif
             !eventually to test index_s_cr values and local values on mynode stored in another array ? 
             !if (isclg_glo <= nzupd0d_oa) then
             !   if (if_mpi_oa) index_s_cr(isclg_glo-1) =  tupd_oa(iv)
             !else
             !   if(if_print_node) write(io_unit,*)'ERROR OA : counts of tot. MPI scalogram mismatch for XIOS index_s_cr axis'
             !   stop
             !endif
#else
            
             if (isclg_glo <= nzupd0d_oa) then
                index_s_cr(isclg_glo-1) =  tupd_oa(iv)
             else
                write(*,*)'ERROR OA : counts of scalogram mismatch for XIOS index_s_cr axis'
                if(if_print_node) write(io_unit,*)'ERROR OA : counts of scalogram mismatch for XIOS index_s_cr axis'
                stop
            endif
#endif
         endif 
      enddo
      if ( isclg_loc .ne. nsclg_loc ) then
        write(*,*)'ERROR OA : counts of MPI scalogram mismatch for XIOS index_s_oa axis'
        if(if_print_node) write(io_unit,*)'ERROR OA : counts of MPI scalogram mismatch for XIOS index_s_oa axis'
        stop
      endif 


#ifdef MPI

!.....BLXD TODO TEST REMOVING communicator size setting since done earlier
!     in sub. init_parameter_oa 
      
      !call MPI_COMM_SIZE(comm, comm_size, ierr)
      recvfromrkminus1 : if ( (mynode>=1) .and. (mynode<=comm_size-1) ) then ! 0->1,...,6->7 rank recieving are [1,2,..,7]
        tag=mynode+mynode-1 ! rank 1 receiving : 1+0=1, rank 2 receiving : 2+1=3,..., rank 7 receiving : 7+6=13
        call MPI_Recv(buff_r,2, MPI_INTEGER, mynode-1, tag, comm, MPI_STATUS_IGNORE, ierr)
        if (ierr/=0) then
            print*,'ERROR OA : MPI_Recv buff_r mynode(<-mynode-1)/tag ',mynode,tag
            stop
        endif
        ! isclg_beg is the global position of the 1st scalogram hold by the MPI process mynode
        ! isclg_beg indexing regardless the real global scalogram index hold by arrays index_s_* 
        ! isclg_beg of mynode proc. = isclg_beg of proc (mynode-1) + nsclg_loc of proc. (mynode-1)
        if ( buff_r(2) > 0 ) then
            !isclg_beg = buff_r(1) + buff_r(2) - 1
            isclg_beg = buff_r(1) + buff_r(2)
        else 
            isclg_beg = buff_r(1)
        endif
        !if (if_mpi_oa) then
        tag=tag+1           ! rank 1 receiving : 1+0 + 1 = 2, rank 2 receiving : 1+2 + 1=4,...
        call MPI_Recv(index_s_cr,nzupd0d_oa, MPI_INTEGER, mynode-1, tag, comm, MPI_STATUS_IGNORE, ierr)
        if (ierr/=0) then
            print*,'ERROR OA : MPI_Recv index_s_cr mynode(<-mynode-1)/tag ',mynode,tag
            stop
        endif
        !endif
      endif recvfromrkminus1

      mpi_rank0 : if (mynode==0) then
        request=-99
        ! isclg_beg is the global position of the 1st scalogram hold by the MPI process mynode
        ! MPI process zero sends this index and the local MPI process # of scalograms
        tag=mynode+mynode+1 ! rank 0 sending : 0+1=1
        isclg_beg=0
        buff_s(1) = isclg_beg  
        buff_s(2) = nsclg_loc
        if (snb) then
            call MPI_ISSend(buff_s, 2, MPI_INTEGER, mynode+1, tag, comm, request, ierr)
            if (ierr/=0) then
                print*,'ERROR OA : MPI_ISsend buff_s mynode(->mynode+1)/tag ',mynode,tag
                stop
            endif
        else
            call MPI_Send(buff_s, 2, MPI_INTEGER, mynode+1, tag, comm, ierr)
            if (ierr/=0) then
                print*,'ERROR OA : MPI_Send buff_s mynode(->mynode+1)/tag ',mynode,tag
                stop
            endif
        endif
        !if (if_mpi_oa) then
        index_s_cr( isclg_beg : isclg_beg + nsclg_loc - 1) = index_s_oa( 0 : nsclg_loc - 1)
        !endif
        if (snb) then
            CALL MPI_Wait(request, MPI_STATUS_IGNORE, ierr)
            if (ierr/=0) then
                print*,'ERROR OA : MPI_Wait buff_s mynode(->mynode+1)/tag/request ',mynode,tag,request
                stop
            endif
        endif
        !if (if_mpi_oa) then
        tag=tag+1           ! rank 0 sending : 0+1 + 1 = 2
        call MPI_Send(index_s_cr,nzupd0d_oa, MPI_INTEGER, mynode+1, tag, comm, ierr)
        if (ierr/=0) then
            print*,'ERROR OA : MPI_Send index_s_cr mynode(->mynode+1)/tag ',mynode,tag
            stop
        endif
        !endif
      endif mpi_rank0 

      send2rankplus1 : if ( (mynode>=1) .and. (mynode<=comm_size-2) ) then ! 0->1,...,6->7 rank sending are [0,1,..,6]
        request=-98
        tag=mynode+mynode+1 ! rank 1 sending : 1+2=3, rank 2 sending : 2+3=5, ..., rank 6 sending : 6+7=13 
        buff_s(1) = isclg_beg  
        buff_s(2) = nsclg_loc
        if (snb) then
            call MPI_ISSend(buff_s,2, MPI_INTEGER, mynode+1, tag, comm, request,ierr)
            if (ierr/=0) then
                print*,'ERROR OA : MPI_Isend buff_s mynode(->mynode+1)/tag ',mynode,tag
                stop
            endif
        else
            call MPI_Send(buff_s, 2, MPI_INTEGER, mynode+1, tag, comm, ierr)
            if (ierr/=0) then
                print*,'ERROR OA : MPI_Send buff_s mynode(->mynode+1)/tag ',mynode,tag
                stop
            endif
        endif
        !if (if_mpi_oa) then
        index_s_cr( isclg_beg : isclg_beg + nsclg_loc - 1) = index_s_oa( 0 : nsclg_loc - 1)
        !endif
        if (snb) then
            CALL MPI_Wait(request, MPI_STATUS_IGNORE, ierr)
            if (ierr/=0) then
                print*,'ERROR OA : MPI_Wait buff_s mynode(->mynode+1)/tag/request ',mynode,tag,request
                stop
            endif
        endif
        !if (if_mpi_oa) then
        tag=tag+1           ! rank 1 sending : 1+2 + 1 = 4, rank 2 : 6,...
        call MPI_Send(index_s_cr,nzupd0d_oa, MPI_INTEGER, mynode+1, tag, comm, ierr)
        if (ierr/=0) then
            print*,'ERROR OA : MPI_Send index_s_cr mynode(->mynode+1)/tag ',mynode,tag
            stop
        endif
        !endif
      endif send2rankplus1


      !if (if_mpi_oa) then
      call MPI_Barrier(comm,ierr)
      if (ierr/=0) then
          print*,'ERROR OA : MPI_Barrier bef. index_s_cr Bcast mynode/ierr ',mynode,ierr
          stop
      endif
      ! If global index_s_cr needed
      call MPI_Bcast(index_s_cr,nzupd0d_oa, MPI_INTEGER, comm_size-1, comm, ierr)
      ! isclg_glo = 0
      ! cfg_loop2 : do ic = 1,nzc_oa
      ! var_loop2 : do iv = begc_oa(ic),begc_oa(ic)-1
      !   ! it is a scalogram analysis
      !   glo_scal2 : if ( tv_sclg_oa(iv) ) then
      !       !tupd_oa(iv)                                      = index_s_cr( isclg_glo )
      !       !tvar_oa(tc_oa(ic), tvc_oa(ic), iv-begc_oa(ic)+1) = index_s_cr( isclg_glo )
      !       sclg2v( index_s_cr(isclg_glo) ) = iv
      !       isclg_glo = isclg_glo + 1
      !   endif glo_scal2
      ! enddo var_loop2
      ! enddo cfg_loop2
      !endif

#endif

! BLXD TODO REMOVE USELESS STUFF >>>>
#ifdef OA_TEST_MPI_XIOS
      if (if_chks) then
          if ( .not. allocated( cntvt_cr) ) then
            allocate( cntvt_cr(nzv_oa) )
            cntvt_cr(1:nzv_oa)=-99
          endif 
      endif
#endif
! BLXD TODO REMOVE USELESS STUFF <<<<

! BLXD TODO REMOVE USELESS STUFF >>>>
#ifdef MPI
       if ( mpi_nonblock_send) then
       ! Allocating MPI request handle for non-blocking pt-pt communication
       ! size of the global # of scalograms so that all process know about the values 
       if ( .not. allocated(sclg_request) ) then
        allocate ( sclg_request(nzv_oa)   ) 
        sclg_request(min(1,nzv_oa):nzv_oa)=-99
       endif
       !if ( .not. allocated(sclgij_request) ) then
       ! allocate ( sclgij_request(nzv_oa)   ) 
       ! sclgij_request(min(1,nzv_oa):nzv_oa)=-99
       !endif
       ! Setting MPI request handle for non-blocking pt-pt communication
        do iv = 1,nzv_oa
          ! it is a scalogram analysis
          if ( tv_sclg_oa(iv) ) then
            ! the iv scalogram analysis is included on 
            ! the itile region and therefore on the MPI process
            sclg_request(iv)  = 1000+iv
            !sclgij_request(iv)= 2000+iv
            !recv_request(iv)  = 5000+iv
          endif 
        enddo
       endif
#endif
! BLXD TODO REMOVE USELESS STUFF <<<<


      return
      end subroutine init_scalogram_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note Croco openMP tile(-thread), this routine must be called under 
!! a tile(-thread) loop (see suboutine init_oa and Croco source main.F).
!
! DESCRIPTION:  initialise la structure spatiale du vecteur d'etat.
!
!> @brief Calculates the size of the horizontal and vertical domain 
!! over which analyses will be conducted i.e., nzvs_oa and nzvs3d_oa respectively
!! All the variable analyses are stored in a single vector, called state vector, 
!! which vectorizes the full set of domain grid-points positions to analyse. 
!! Variables enabling conversion from state vector index to domain indices 
!! are initialized.
!
!> @details Spatial domain options are set in the OA namelist with SWT_D_W parameter.  
!! - if set to 1: analysis is requested over the region [lat1,lat2],[lon1,lon2]. 
!! - if set to 2: analysis is requested at point (I,J).
!! - if set to 3: analysis is requested over the ocean colum [Hmin,Hmax] in the horizontal 
!!   domain the (Lat,Lon)
!! The domain to analyse is more over defined by the namelist parameters:
!! - dx_oa(iv_g), dy_oa(iv_g) : steps in i,j indices
!! - k_oa(1,iv_g),k_oa(2,iv_g),dk_oa(iv_g) : first, last k indices and k index step
!!
!! Variables:
!! - nzvs_oa,   begvs_oa    : 2d struture of the state vector
!! - nzvs3d_oa, begvs3d_oa  : 3d structure of the state vector
!! - ij2l_oa, l2i_oa, l2j_oa : conversion l <--> (i,j)
!!
!! The section of the state vector defined by indices begvs3d_oa(ls_l), begvs3d_oa(ls_l+1)-1
!! corresponds to the points to analyse over one given (i,j) ocean column for a given variable
!! with index iv_g. For each variable, the number of (i,j) points to analyse is given by index ls_l
!! which ranges from begvs_oa(iv_g) to begvs_oa(iv_g+1)-1. 
!!
!! For each variables: 
!! - the grid point is needed to get the mask value, i.e., tgv_oa.
!! - the 2d/3d attribute is required to construct the vertical part of the state vector, i.e., tgv3d_oa.
!!
!! The kmin3d_oa parameter stores the lowest analysed k index.
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - Benedicte Lemieux-Dudon
!!  - Namelists (06/2014) + comments + Stand-alone version + optimization (01/2015)
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!!  - Toward a Croco Tile-thread compliant interface and OA version (2021)
!!  - Developments to perform online Scalogram analyses (06/2021)
!!   - handles the scalograms distribution over MPI process/subdomains
!!      either with XIOS or using MPI internal OA instructions (see if_mpi_oa)
!!   - handles the scalogram distribution over Croco tile(-threads)
!> @date 2021
!> @todo BLXD
!! - Croco tile(-thread) to test with working Croco tile-thread version,
!! - continue the effort to synthesise redundant instructions with subroutine calls
!------------------------------------------------------------------------------
      subroutine var_space_oa( flag_s, tile )

      use module_grd_oa, only  : grd, grid_str_oa          &
     &                          ,get_3D_grid_mask          &
     &                          ,get_grid_mask             &
     &                          ,get_2D_subdomain_minmax
      use module_tile_oa, only : st, tile_space_str, scl, tile_sclg_str
      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      !use module_oa_stock

      implicit none

      logical, intent(in) :: flag_s
      integer, intent(in) :: tile
      
      ! BLXD_TILE_ISSUE 
      ! pst declaration moved from module_tile_oa to local subroutine
      ! => openMP PRIVATE variable with respect to Croco tile-thread loop
      type(tile_space_str), pointer :: pst => null()

      !BLXD_TILE_ISSUE rm pg from list of public module variables add grid_str_oa
      ! => openMP PRIVATE variable with respect to Croco tile-thread loop
      type(grid_str_oa), pointer :: pg => null()

      !BLXD_TILE_ISSUE rm pscl from list of public module variables add tile_sclg_str
      ! => openMP PRIVATE variable with respect to Croco tile-thread loop
      type(tile_sclg_str), pointer :: pscl => null()


      integer :: i, j, k, grd_pt_code                                 &                          
           ,iv_g                                                      &
           ,mykmin                                                    &
           ,msk_g                                                     &
           ,ls_g

      real    ::                                                      &
            lat_g                                                     &
           ,lon_g                                                     &
           ,h_g                                                         

      real    ::                                                      &
            latmin_g                                                  &
           ,lonmin_g                                                  &
           ,latmax_g                                                  &
           ,lonmax_g

      integer :: io_nodoa, num !, ii_glob, jj_glob
      integer          :: imin, imax, jmin, jmax, kmin, kmax
      integer, pointer :: nzvs_oa => null(), nzvs3d_oa => null()

!.... Grid and Space structure : temporary pointer association
      pg  => grd(tile)
      pst => st(tile) 
      if (scalogram_analysis) pscl => scl(tile)

!.... Initialisation of the size of the spatial structure
      pst%nzvs_oa   = 0
      pst%nzvs3d_oa = 0      

!.... Local pointers pointing toward 2D/3D the spatial structure size, resp.
!     just to handle the counts
      nzvs_oa => pst%nzvs_oa ; nzvs3d_oa => pst%nzvs3d_oa

!.... Checking grid size consistency

      if ( ( pg%imin /= pst%imin ) .or. &
           ( pg%imax /= pst%imax ) .or. &
           ( pg%jmin /= pst%jmin ) .or. &
           ( pg%jmax /= pst%jmax ) .or. &
           ( pg%kmin /= pst%kmin ) .or. &
           ( pg%kmax /= pst%kmax ) ) then
!$OMP ATOMIC
        write (*,*) 'ERROR OA : error with tile dim for grid and space structure tile/node', tile, mynode
        stop
      else
!.... Local variable (just to shorten loop instructions)
        imin=pg%imin ; imax=pg%imax ; kmin=pg%kmin
        jmin=pg%jmin ; jmax=pg%jmax ; kmax=pg%kmax
      endif

      spa_var_loop : do iv_g = 1 , nzv_oa

!.... Local/private openMP parameter => write log
         if (flag_s) then
             io_nodoa = set_io_nodoa(iv_g,mynode,8,7) ! [7/8](odd/even-iv)000+mynode
         else
             io_nodoa = set_io_nodoa(iv_g,mynode,6,5) ! [5/6](odd/even-iv)000+mynode
         end if

         if (flag_s.eq..false.) then


!.... Check user values : openMP tile-thread modifying openMP SHARED module var. implicit synchro END SINGLE
!$OMP SINGLE
            if (k_oa(1,iv_g).eq.-99) then
               k_oa(1,iv_g) = 1 
            endif
            if (k_oa(2,iv_g).eq.-99) then
               k_oa(2,iv_g) = kmax
            endif

            if ( k_oa(1,iv_g) .gt. k_oa(2,iv_g) ) then
               write (*,*) " k_oa(2,:) trop grand > kmax = N !"
               if(if_print_node) write (io_unit,*) " k_oa(2,:) trop grand > kmax = N !"
               stop
            endif
            if ( k_oa(2,iv_g).gt.kmax ) then
               write (*,*) " k_oa(2,:) trop grand > kmax = N !"
               if(if_print_node) write (io_unit,*) " k_oa(2,:) trop grand > kmax = N !"
               stop
            endif
            if ( k_oa(1,iv_g).lt.kmin ) then
               write (*,*) " k_oa(1,:) trop petit < kmin !"
               if(if_print_node) write (io_unit,*) " k_oa(1,:) trop petit < kmin !"
               stop
            endif

!.......... BLXD Croco T-grid indices 1:N, W-grid indices from 0:N and kmin=1,kmax=N 
            if ( k_oa(1,iv_g).eq.0 .and.       &
                (tgv_oa(tv_oa(iv_g)).eq.1      &
             .or.tgv_oa(tv_oa(iv_g)).eq.2      &
             .or.tgv_oa(tv_oa(iv_g)).eq.3      &
! Humm p grid point is the 2D vorticity grid point
! Should also start to 1 and not zero ( not w-grid point )
!            .or.tgv_oa(tv_oa(iv_g)).eq.3      &
               ) ) then
                k_oa(1,iv_g)=1
                if (k_oa(2,iv_g).eq.0) k_oa(2,iv_g)=1
            endif
!$OMP END SINGLE

            !if_domain_lonlat_1 : if (swt_d_oa(iv_g).eq.1.or.swt_d_oa(iv_g).eq.3) then
                ! BLDX TODO ADD here potential checks dependent of imin,imax,.. ?
                !if ( if_chks ) then
                !    call domain_checks()
                !end if
                !if ( if_extend_dom ) then
                !    grd_pt_code = mod(tgv_oa(tv_oa(iv_g)),5)
                !    call get_2D_subdomain_minmax( imin,imax,jmin, jmax &
                !                          ,grd_pt_code, pg               &
                !                          ,latmin_g, lonmin_g, latmax_g, lonmax_g )
                !end if
            !end if if_domain_lonlat_1

         else ! flag_s .true. 


            ! BLXD nzvs_oa => pst%nzvs_oa which is properly initialized to zero
            pst%begvs_oa( iv_g ) = nzvs_oa + 1

         end if

!.....configuration spatiale n°1 ou 3: [ lon,lat domain ] x [kmin,kmax,dk] or [hmin,hmax] x [kmin,kmax,dk]

         spatial_config1_3 : if (swt_d_oa(iv_g).eq.1.or.swt_d_oa(iv_g).eq.3) then


            jloop : do j = jmin, jmax, dy_oa(iv_g)
                iloop : do i = imin, imax, dx_oa(iv_g)

                  grd_pt_code = mod(tgv_oa(tv_oa(iv_g)),5)
                  call get_3D_grid_mask( i, j , kmax, grd_pt_code, pg, &
                                         h_g, lat_g, lon_g, msk_g )

#ifdef SPHERICAL
         if ( if_sphe_deg2rad ) then
                lat_g=lat_g*pi_oa/180.D0
                lon_g=lon_g*pi_oa/180.D0
         endif
#endif

                  if_spatial_dom_ok_loopij : if (                     &
                       msk_g.eq.1              .and.                  &
                       lat_g.ge.lat_oa(1,iv_g) .and.                  &
                       lat_g.le.lat_oa(2,iv_g) .and.                  &
                       lon_g.ge.lon_oa(1,iv_g) .and.                  &
                       lon_g.le.lon_oa(2,iv_g) .and.                  &
                       ( (h_g.ge.h_oa(1,iv_g)  .and.                  &
                       h_g.le.h_oa(2,iv_g) )                          &
                       .or.swt_d_oa(iv_g).ne.3) )                     &
                       then


!-------->structure 2d:

                     nzvs_oa          = nzvs_oa + 1
                     if (flag_s) then
                     pst%l2i_oa (nzvs_oa)  = i
                     pst%l2j_oa (nzvs_oa)  = j
                     pst%ij2l_oa(i,j,iv_g)= nzvs_oa
                     endif

!-------->structure 3d:
! #BLXD line_begin - line_end => SUBROUTINE  (step wth get_grid_mask)
                     mykmin=0
                     ! BLXD nzvs3d_oa and nzvs_oa are local pointers 
                     ! pointing twd pst%nzvs,nzvs3d_oa_oa 
                     ! - they have been properly initialized to zero
                     ! - they point twd tile (openMP thread PRIVATE) distinct memory address    
                     if (flag_s) pst%begvs3d_oa(nzvs_oa) = nzvs3d_oa + 1
                     if (tgv3d_oa(tv_oa(iv_g)).eq.2) then
                        nzvs3d_oa          = nzvs3d_oa + 1
                     else
                        do k = k_oa(1,iv_g),k_oa(2,iv_g),dk_oa(iv_g)
                          grd_pt_code = mod(tgv_oa(tv_oa(iv_g)),5)
                          call get_grid_mask( i, j , k, grd_pt_code, pg, msk_g )
                          if ( msk_g.eq.1 ) then
                             if ( mykmin.eq.0.and.flag_s ) then
                                pst%kmin3d_oa(nzvs_oa)=k            
                                mykmin=1
                             endif
                             nzvs3d_oa          = nzvs3d_oa + 1
                          endif
                        enddo
                     endif
! #BLXD line_end                  
                  endif if_spatial_dom_ok_loopij

               enddo iloop
            enddo jloop

         endif spatial_config1_3
         
!.....configuration spatiale n°2: [ lon,lat domain ] x [kmin,kmax,dk] or [hmin,hmax] x [kmin,kmax,dk]

         spatial_config2 : if (swt_d_oa(iv_g).eq.2) then
            i = ptij_oa(1,iv_g)
            j = ptij_oa(2,iv_g)


            if_spatial_dom_ok_pointij : if ( i>=imin .and. i<=imax .and. j>=jmin .and. j<=jmax ) then


            grd_pt_code = mod(tgv_oa(tv_oa(iv_g)),5)
            call get_grid_mask( i, j , kmax, grd_pt_code, pg, msk_g )

            if_mask_ok_pointij : if ( msk_g.eq.1 ) then
               nzvs_oa          = nzvs_oa + 1

               if (flag_s) then
                pst%l2i_oa (nzvs_oa)  = ptij_oa(1,iv_g)
                pst%l2j_oa (nzvs_oa)  = ptij_oa(2,iv_g)
                pst%ij2l_oa(i,j,iv_g)= nzvs_oa
                if ( tv_sclg_oa(iv_g) ) then
!.................. BLXD Scalogram Local tile(-thread)/MPI process counts 
!                        + convert iv_g to tile scalogram index
                    pscl%nsclg        =  pscl%nsclg + 1
                    pscl%v2sclg(iv_g) =  pscl%nsclg
                endif
               endif

!-------->structure 3d:
! #BLXD line_begin 2 TIMES THESE LINES => subroutine (step wth get_grid_mask)
               mykmin=0
               if (flag_s) pst%begvs3d_oa(nzvs_oa) = nzvs3d_oa + 1
               if (tgv3d_oa(tv_oa(iv_g)).eq.2) then
                  nzvs3d_oa          = nzvs3d_oa + 1
               else
                  do k = k_oa(1,iv_g),k_oa(2,iv_g),dk_oa(iv_g)
                    grd_pt_code = mod(tgv_oa(tv_oa(iv_g)),5)
                    call get_grid_mask( i, j , k, grd_pt_code, pg, msk_g )
                     if ( msk_g.eq.1 ) then
                        if ( mykmin.eq.0.and.flag_s ) then
                           pst%kmin3d_oa(nzvs_oa)=k            
                           mykmin=1
                        endif
                        nzvs3d_oa          = nzvs3d_oa + 1
                     endif
                  enddo
               endif
! #BLXD line_end

            endif if_mask_ok_pointij

            endif if_spatial_dom_ok_pointij

         endif spatial_config2

      enddo spa_var_loop

!.....last point:
      if (flag_s) then

       pst%begvs3d_oa( nzvs_oa + 1 ) = nzvs3d_oa + 1
       pst%begvs_oa  ( nzv_oa  + 1 ) = nzvs_oa + 1

   
      endif

!.....Local pointers nullified
      nzvs_oa   => null()
      nzvs3d_oa => null()

!.....Specific module_tile_oa and module_grd_oa pointer nullified
      pg  => null()
      pst => null() 
      if (scalogram_analysis) pscl => null()

      return
      end subroutine var_space_oa


!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note Croco tile(-threads) : this routine must be called before any
!!       openMP parallel region fork and whitout loop over tiles
!!       so that the only MASTER thread initializes the implicitly 
!!       openMP SHARED module variables without any data race issue.
!
! DESCRIPTION: 
!
!> @brief Sets the frequency parameters as requested in namelist 
!! which controls the "state vector" frequency structure.
!
!> @details No tile(-thread) loop => master thread updates shared global module parameters (no data race)
!! Variable :
!! - nzvp_oa is the total number of periods to analyse in the simulation (all variables included). 
!! - begvp_oa stores the number of periods per variables. It is refered as to the "frequency" structure
!!   of the state vector. 
!! Two options : (1) discrete period of analysis (including scalogram analysis), (2) integrated period of analysis. 
!!
!! Namelist parameters are:
!! PER_W(1), per_oa(1) : minimum period
!! PER_W(2), per_oa(2) : maximum period
!! PER_W(3), per_oa(3) : step determining the number of periods to be examined, dp in sec. (-99=optimum).
!! NZPT_PER_W, nzpt_per_oa : number of points in the convolution window.
!! 
!! Outputs:
!! - resv_oa : integer setting the temporal resolution (i.e., number of model time step) at which
!!   the convolution product will be calculated (e.g., one point every 4 model time steps)
!!   TODO check this definition in the case dori_oa equals 2.
!! - nzvp_oa and begvp_oa vector : parameters controling the frequential structure of the "state vector".
!! - perv_oa(1,:) : series of atom periods which will be examined for each cfg/var ( optimal steps if PER_W(3) = -99 ) 
!! - perv_oa(2,:) : reconstruction factor ( replaces inverse transformation ).
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - Benedicte Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - Croco-OA interface 1st version (2020)
!!  - Developments to perform online Scalogram analyses (2021)
!!    - from namelist options to calculation + XIOS netcdf outputs
!!    - scalograms distribution over MPI process/subdomains handled 
!!      either by XIOS or by MPI internal OA instructions.
!!  - Implementation of the Heisenberg theoritical and numerical uncertainties. 
!!    Onging (2021). 
!
!> @date 2021
!> @todo  BLXD
!! - remove goto syntax ? intent attribute ?
!! - modules with atoms specific paramters, functions, methods,..
!! - extend heisenberg boxes to other wavelet atoms.
!------------------------------------------------------------------------------
      subroutine var_per_oa( flag_p, dti )

      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      !use module_oa_stock

      implicit none

      double precision, intent(in) :: dti                                !< Integration time step

      integer                                                         &
         iv_p                                                         &  !< Variable index
        ,iv1_p

      real                                                            &
         ip_p

      real ::                                                         &   
         p_p                                                          &  !< Time Period, Scale
        ,t0c,dtc,w00c,dwc                                                 !< Wavelet resolution : Heisenberg parameters

      logical :: flag_p

      ! BLXD Morlet wavelet with sigma=1
      ! 1) Frequency center of the mother wavelet w0c = 2. * pi_oa * fc      
      !     => For the scaled wavelet w0c = 2. * pi_oa * fc / s
      !    Frequency spread of the mother wavelet dwc =  1./sqrt(2.) ~ 0.707106781
      !     => For the scaled wavelet dwc = 1. / sqrt(2.) / s     
      ! 2) Time Center of the mother wavelet t0c = 0. 
      !     => For the scaled wavelet consider the shift in time
      !    Time spread for the mother wavelet dtc = 1. / sqrt(2.)
      !     => For the scaled wavelet dtc = s / sqrt(2.)
      ! 3) Heisenberg box dtc x dwc = 1/2 = ( s * sigma / sqrt(2.) ) x ( 1. / ( s * sigma * sqrt(2.) ) )
      ! Note that the following period progression ratio is independent of s :
      ! (p2-p1)/(p2+p1) = dwc(scaled wavelet) / w0c(scaled wavelet) 
      !                 = dwc(mother wavelet) / w0c(mother wavelet)    
      !                 = 1. / ( 2. * pi_oa * fc_oa  * sqrt(2.) )
      !                 = 1. / ( 6.                  * sqrt(2.) )  
      ! 
      ! Declare a module specific to Morlet wavelet with declared cst and particular values
      ! Currently Heisenberg boxes only valid for the case of the Morlet wavelet
      ! Set gard-rail 

      t0c = 0.
      ! w0c redefined as module parameter 
      w00c = 6.0         
      dtc = 1./sqrt(2.) 
      dwc = dtc         
 
      nzvp_oa = 0

!.....precalcul des periodes blxd cas discret (1 ou +s per par var iv_p) ou integre :
      iv_loop1 : do iv_p = 1, nzv_oa

! cas resolution reelle
        if (per_oa(3,iv_p).eq.-99.*unite_oa) then
         
         if (flag_p) begvp_oa( iv_p ) = nzvp_oa + 1
         if (per_oa(1,iv_p).eq.0) per_oa(1,iv_p)=1. ! on commence pas a 0      
         p_p=per_oa(1,iv_p) 


         do while (p_p.le.per_oa(2,iv_p)) 
          nzvp_oa = nzvp_oa + 1

          if (flag_p) perv_oa (1,nzvp_oa) = p_p
          p_p=p_p*(1.+dwc/w0c)/(1.-dwc/w0c)
         enddo
         nzvp_oa = nzvp_oa + 1
         if (flag_p) perv_oa (1,nzvp_oa) =p_p

! cas resolution reguliere 
       else
         if (flag_p) begvp_oa( iv_p ) = nzvp_oa + 1
!---------------------------------------------------------
!.......periode entiere en secondes
!---------------------------------------------------------
        ! change for loop do ip_p = int(per_oa(1,iv_p)),int(per_oa(2,iv_p)),int(per_oa(3,iv_p))
        ! integer :: ip_p

        ip_p = per_oa(1,iv_p)

 100    continue
         nzvp_oa          = nzvp_oa + 1
         if (flag_p) perv_oa(1,nzvp_oa) = ip_p
         ip_p = ip_p + per_oa(3,iv_p)
        if (ip_p.le.per_oa(2,iv_p)) goto 100 

!---------------------------------------------------------
!       enddo
!---------------------------------------------------------
       endif

      enddo iv_loop1
      if (flag_p) begvp_oa(nzv_oa+1)     = nzvp_oa + 1

      if (.not.(flag_p)) return
      
!.....precalcul des resolutions:


      do iv_p = 1 , nzv_oa

!.......counting # of periods for scalogram variables
        nper_sclg(iv_p) = 0

        do iv1_p = begvp_oa(iv_p) , begvp_oa(iv_p+1)-1

!.........atome= dirac:
          if (tpsi_oa(iv_p).eq.0.and.flag_p) then
            resv_oa (iv1_p)    = 1 
          endif

!.........atome= ondelette de morlet ou windowed fourier ou fourier:
          wavelet : if (tpsi_oa(iv_p).ne.0.and.flag_p) then

           if ( tv_sclg_oa(iv_p) ) nper_sclg(iv_p) = nper_sclg(iv_p) + 1        

           single_discrete_per : if (dori_oa(iv_p).eq.1) then 
            ! Discrete
            ! BLXD 2020 consistent changes in var_time_oa
            ! resv_oa is an integer => automatic f2008 conv real to integer 
            !                          follows use of the intrinsic function INT
            if ( tv_sclg_oa(iv_p) ) then
                ! BLXD If scalogram better to calc. resolution using the smallest period
                !      PERV_W(1):PERV_W(2) same namelist => NZPT_PER_OA identical
                if ( MOD( perv_oa (1,begvp_oa(iv_p)), dti) /= 0 ) then   ! iv1_p replaced by begvp_oa(iv_p)
                    resv_oa (iv1_p)   = max ( (  int(                         &
                                     perv_oa (1,begvp_oa(iv_p))  /   dti      &
                                      ) + 1 )                                 &
                                   / nzpt_per_oa (iv_p)                       &
                                             , 1 )      
                else
                    resv_oa (iv1_p)   = max ( (  int(                         &
                                     perv_oa (1,begvp_oa(iv_p))  /   dti      &
                                      ) )                                     &
                                   / nzpt_per_oa (iv_p)                       &
                                             , 1 )      
                endif 

            else

                if ( MOD( perv_oa (1,iv1_p), dti) /= 0 ) then
                    resv_oa (iv1_p)   = max ( (  int(                         &
                                     perv_oa (1,iv1_p)  /   dti               &
                                      ) + 1 )                                 &
                                   / nzpt_per_oa (iv_p)                       &
                                             , 1 )      
                else
                    resv_oa (iv1_p)   = max ( (  int(                         &
                                     perv_oa (1,iv1_p)  /   dti               &
                                      ) )                                     &
                                   / nzpt_per_oa (iv_p)                       &
                                             , 1 )      
                endif 

            endif
   
           else single_discrete_per
            ! Integrated period
            ! BLXD 2020 consistent changes in var_time_oa
            ! resv_oa (iv1_p)   = max ( (  int(                         &
            !                  perv_oa (1,begvp_oa(iv_p))               & 
            !                /   dti                                    &
            !                   ) + 1 )                                 &
            !                           , 1)                            &
            !                / nzpt_per_oa (iv_p)
            ! #BLXD : BUG ?
            ! The above expression (from original code) is suspicious
            ! Why the MAX function should be here taken bef. division by nzpt_per_oa ?
            ! RESV_OA is supposed to be an integer with unit being a number of time steps
            !  
            ! BLXD PERV_W(1):PERV_W(2) same namelist => NZPT_PER_OA identical
            ! Integrated period : Again the best choice
            ! a) the smallest period of the integrated period range is used
            !   to calculate the resolution in terms of number of time steps 
            ! b) the largest period to calculate the convolution window
            if ( MOD( perv_oa (1,begvp_oa(iv_p)), dti) /= 0 ) then
                resv_oa (iv1_p)   = max ( (  int(                         &
                                 perv_oa (1,begvp_oa(iv_p))               & 
                               /   dti                                    &
                                  ) + 1 )                                 &
                               / nzpt_per_oa (iv_p)                       &
                                          , 1)                            

            else
                resv_oa (iv1_p)   = max ( (  int(                         &
                                 perv_oa (1,begvp_oa(iv_p))               &
                               /   dti                                    &
                                  ) )                                     &
                               / nzpt_per_oa (iv_p)                       &
                                         , 1 )      

            endif

           endif single_discrete_per
           endif wavelet
        enddo
       enddo

       nper_sclg_max = maxval( nper_sclg(:) )

      return

      end subroutine var_per_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note Croco tile(-threads) : this routine must be called before any
!!       openMP parallel region fork and/or whitout loop over tiles
!!       so that the only MASTER thread allocates the implicitly 
!!       openMP SHARED module variables without any data race issue.
!
! DESCRIPTION: allocation de variables globales definissant la structure frequentielle
!              du vecteur d'etat.
!
!> @brief allocates global variables defining the "state vector" frequency structure.
!
! REVISION HISTORY:
!
!> @authors B. Lemieux-Dudon
!!   Adapted from Symphonie/NHOMS initial version (2006)
!> @date 2021
!> @todo change subroutine name part3 -> per ? 
!------------------------------------------------------------------------------

      subroutine allocate_part3_oa

      ! BLXD change name allocate_nzvp_oa ? 
      use module_oa_variables
      use module_oa_time
      use module_oa_periode

      implicit none

      if ( .not. allocated(resv_oa) )                                 &
      allocate (                                                      &   
       resv_oa (nzvp_oa)                                              &  ! resolution temporelle pour le calcul de la convolution
       )

      if ( .not. allocated(perv_oa) )                                 &
      allocate (                                                      &  ! periodes associees a la structure vectorielle du vecteur d etat, 
       perv_oa (2,nzvp_oa)                                            &  ! facteurs de reconstruction associes a l'ondelette
       )

      !> Info sur la 1e/derniere analyse OA valide de la simu
      if ( .not. allocated(lper_fst_oa)     )  allocate ( lper_fst_oa (nzv_oa,nzvp_oa)       )
      if ( .not. allocated(if_first_per_oa) )  allocate ( if_first_per_oa (nzv_oa,nzvp_oa)   )
      if ( .not. allocated(lper_lst_oa)     )  allocate ( lper_lst_oa (nzv_oa,nzvp_oa)       )
      if ( .not. allocated(if_last_per_oa)  )  allocate ( if_last_per_oa (nzv_oa,nzvp_oa)    )

      !> Number of scalogram periods (i.e., scales)
      if ( .not. allocated(nper_sclg)       )  allocate( nper_sclg(nzv_oa) )

      return
      end subroutine allocate_part3_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note This routine is called under a Croco tile(-thread) loop and requires
!!  Preventing double allocation issues and openMP data race
!!  (all threads updating implicitly SHARED global module parameters),
!!  and Avoidingd multiple initialization of OA module parameters (tile-only simu),
!!  only the last tile(-threads) calls var_time_oa (with openMP critical construct)
!
! DESCRIPTION: 
!
!> @brief Sets the temporal convolution window parameters for all the requested analysis.
!
!> @details Variables :
!! - nzvt_oa is the total number of time windows requested for analysis (all variables included). 
!! - begvt_oa stores the number of time windows per variables. It is refered to as the "time" structure
!!   of the state vector.
!!
!! Corresponding namelist parameters :
!! - DELTA_T_W, delta_t_oa : set to 1 for Dirac analysis, set to 2 or more for wavelets or
!!   windowed Fourier. In the latter case, the convolution window is set 2 times the period 
!!   of interest (c.f., perv_oa(1,ip_t) ).
!! - SWT_T_W, swt_t_oa : defines the time domain where analysis are requested
!!   * 1 : from T1 to T2 #BLXD to test implemted subroutine tool_datetosec TODO
!!   * 2 : single date T1 at the end, NOT AVAILABLE YET TO TEST 
!!   * 3 : single date T1, NOT AVAILABLE YET TO TEST
!!   * 4 : the entire simulation.
!!   Check : the only namelist parameters are DATE_DEB, DATE_END (pattern is 2011, 01, 1, 00, 01, 18).
!!   kount_user_oa(1:2,nzv_oa) is calculated from DATE_DEB, DATE_END in terms of model time steps.
!! - DT_W : time period at which the analysis is repeated (s). Enables to deduce kount_user_oa(3,nzv_oa) 
!!   in terms of model time steps. 
!!   
!! - lt_p = begvt_oa(iv), begvt_oa(iv+1)-1
!!   is a "pointer" to all the requested analysis for a given variable-configuration
!!   including the requested number of simulation steps where to perform the analysis and all the requested period of analysis
!!   For a given lt_p, one can retreive the corresponding analysis time period per_t2p_oa(lt_p)
!!   To each lt_p corresponds :
!!   1) a convolution window [kountv_oa(1,lt_p), kountv_oa(2,lt_p)] centered at the simulation time step 
!!      [kountv_oa(2) + kountv_oa(1)]/2 and set according to :
!!   - the type of atom tpsi_oa(iv_p)
!!   - the analysis time period (delta_t_w) 
!!   2) a period of analysis recovered through ip_t = per_t2p_oa(lt_p)
!!      the requested period can be :
!!      - "discrete" and for each variable-cfg retrieved with ipt = begvp_oa(iv_t), begvp_oa(ivt+1)-1
!!         if dori_oa(iv_t) is set to 1 
!!
!! Outputs: 
!! - kountv_oa(1,nzvt_oa), kountv_oa(2,nzvt_oa) are delimiting each convolution window (for each variables).
!! - per_t2p_oa (nzvt_oa) : correspondance between "temporal" and "frequency" structure
!! - dkount_oa  (1:nzv_oa): duree d'echantillonnage,
!! - nztime_v_oa(1:nzv_oa): nombre reel de points d'echantillonnage par periode.    
!! - nzvt_oa, begvt_oa    : structure temporelle du vecteur d'etat.
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - Headers, comments, cleaning, from f77 to earlier standards  
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!!  - Toward Croco Tile-thread compliant OA version (2021)
!!  - Development of scalogram analysis (ocean field time series analysis)
!> @date 2021
!> @todo BLXD
!------------------------------------------------------------------------------
      subroutine var_time_oa( flag_t, kount0, nt_max, dti )
      
      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      !use module_oa_stock

      implicit none

      double precision, intent(in) :: dti                                !< Integration time step
      integer, intent(in)          :: kount0,                         &  !< First simulation iteration
                                      nt_max                             !< Last integration iteration

      logical, intent(in) :: flag_t

      integer                                                         &
            k                                                         &
           ,dkount_tot_t                                              &
           ,iv_t                                                      &
           ,ip_t                                                      &
           ,ip1_t                                                     &
           ,ip2_t

      integer, dimension(1:nzv_oa) :: comptt

      integer :: ncomptt

      nzvt_oa  = 0
      ncomptt  = 0


!---- >boucle sur toutes les variables:
      loop_on_variables : do iv_t = 1, nzv_oa

! #BLXD TODO simplify those variables
          if_first_rec_oa(iv_t)=0
          ltrec_fst_oa(iv_t) = -9999
          if_last_rec_oa(iv_t)=0
          ltrec_lst_oa(iv_t) = -9999

!*******************************************************************
!     preparation de la structure temporelle
!*******************************************************************
         
         if (flag_t) begvt_oa(iv_t)   = nzvt_oa + 1

         if (dori_oa(iv_t).eq.1) then 
            ! Discrete period (including scalogram discrete list of periods)
            ! lp_m pointer is given by per_t2p_oa(lt_m) and associated to each lt_m pointer
            ip1_t = begvp_oa( iv_t     ) 
            ip2_t = begvp_oa( iv_t + 1 ) - 1
         else
            ! Integrated period : 
            ! BLXD WARNING :
            ! - ip_t loop reduces to the last period pointer
            ! - per_t2p_oa(lt_m) is associated to the last period only
            ip1_t = begvp_oa( iv_t + 1 ) - 1 
            ip2_t = begvp_oa( iv_t + 1 ) - 1
         endif

!----->boucle sur toutes les periodes demandees pour cette variable:
         loop_on_requested_period : do ip_t = ip1_t , ip2_t

          ! only useful for discrete period
          if_first_per_oa(iv_t,ip_t) = 0
          lper_fst_oa    (iv_t,ip_t) = -9999
          if_last_per_oa (iv_t,ip_t) = 0
          lper_lst_oa    (iv_t,ip_t) = -9999

!----------------------------------------------------------------
!     precalcul de la largeur en kount de l'atome (from time periods perv_oa):
!     dkount_tot_t for 0: Dirac, 1: Wavelette, 2: Windowed Fourier
!----------------------------------------------------------------

!...........atome = dirac:
            if (tpsi_oa(iv_t).eq.0) then
               if (perv_oa(1,ip_t).eq.0.) then
                  dkount_tot_t    = 0
               else
                  dkount_tot_t    = int(perv_oa(1,ip_t)/2./dti-0.5)+1
               endif
            endif   

!...........atome = ondelette de morlet ou windowed fourier:     
            if_morlet_wind_fourier :if (tpsi_oa(iv_t).eq.1.or.tpsi_oa(iv_t).eq.2) then
               ! #BLXD 2020 WARNING : introduction of a special case when perv is diviqible by dti 
               ! The integer parameter dkount_tot_t defines the width of the convolution window 
               ! The convolution window centered at the k-time index is latter defined as
               ! the the time index intervals [k-dkount_tot_t,k+dkount_tot_t].
               ! Objective : setting dkount_tot_t in order to define a convolution window 
               !             at least as large as 2 times the time period requested by 
               !             the user ( ie, perv_oa )
               ! To meet this objective :
               ! dkount_tot_t is not taken as the integer part of the division of perv_oa by dti 
               ! (ie, p) but set to p+1 (ie, ceiling operation since perv_oa is positive):
               !       p * dti < perv_oa <= (p + 1) * dti
               ! However :
               ! Being given the user defined multiplicative factor delta_t_oa >= 1, 
               ! this implies "unexpected total number of wavelet analysis" 
               ! when the simulation lasts exactly p times the time period perv_oa
               ! Proposed modification : introduction of a special case when perv is divible by dti 
               ! TODO extend this tho var_per_oa (see resv_oa)
               if (dori_oa(iv_t).eq.1) then
                  ! If scalogram we use the same convolution window for all the periods
                  ! taking the limiting one which is the largest period conv. wind
                  if ( .not. tv_sclg_oa(iv_t)  ) then
                      if ( MOD( perv_oa (1,ip_t), dti) /= 0 ) then  
                      dkount_tot_t    = (int(                             & ! BLXD original
                           perv_oa (1,ip_t)  /   dti                      &
                           ) + 1)                                         &
                           * delta_t_oa(iv_t)
                      else 
                        dkount_tot_t    = int(                            & ! MOD = 0
                            perv_oa (1,ip_t) * delta_t_oa(iv_t) /   dti )                             
                      endif 
                  else  ! Salogram con wind from largest period associated to iv_t
                      if ( MOD( perv_oa (1,ip2_t), dti) /= 0 ) then
                      dkount_tot_t    = (int(                             &
                           perv_oa (1,ip2_t)  /   dti                     &
                           ) + 1)                                         &
                           * delta_t_oa(iv_t)
                      else 
                        dkount_tot_t    = int(                            &
                            perv_oa (1,ip2_t) * delta_t_oa(iv_t) /   dti )                             
                      endif 
                  endif
               else
                  if (ip2_t .ne. ip_t) then
                    write(*,*)'INTEGRATED PERIOD CASE in var_time_oa iv_t, ip_t',iv_t, ip_t
                    write(*,*)'ERROR OA : ip2_t msut be equal to ip_t ',ip2_t,ip_t
                    stop
                  endif
                  ! BLXD Useless same as dori_oa(iv_t).eq.1
                  ! For integrated period we use the same convolution window for all the periods
                  ! taking the limiting one which is the largest period conv. wind
                  ! since ip_t loop goes from ip2_t to ip2_t
                   if ( MOD( perv_oa (1,ip_t), dti) /= 0 ) then
                   dkount_tot_t    = (int(                             & ! BLXD original
                        perv_oa (1,ip_t)  /   dti                      &
                        ) + 1)                                         &
                        * delta_t_oa(iv_t)
                   else 
                     dkount_tot_t    = int(                            & ! MOD=0
                         perv_oa (1,ip_t) * delta_t_oa(iv_t) /   dti )                             
                   endif 
               endif

            endif if_morlet_wind_fourier

!...........atome = fourier :
            if (tpsi_oa(iv_t).eq.3) then
               dkount_tot_t = 99999
            endif   

!...........initialisation du compteur de la variable:
            comptt(iv_t)=0

!----------------------------------------------------------------
!     configuration temporelle numero 1: periode
!     configuration temporelle numero 4: simu
!     
!     atomes differents de fourier => Dirac, Morlet or windowed Fourier
!----------------------------------------------------------------
            not_fourier : if ((swt_t_oa(iv_t).eq.1.or.swt_t_oa(iv_t).eq.4).and.tpsi_oa(iv_t).ne.3 ) then

               simulation_window : do k = kount_user_oa(1,iv_t),kount_user_oa(2,iv_t),kount_user_oa(3,iv_t)
                  nzvt_oa             = nzvt_oa + 1

                  if (flag_t) then ! second pass
                   kountv_oa(1,nzvt_oa) = k - dkount_tot_t
                   kountv_oa(2,nzvt_oa) = k + dkount_tot_t
                   per_t2p_oa (nzvt_oa)  = ip_t

!----------------------------------------------------------------
                   ! #BLXD 2020 change 
                   ! The nzvt_oa convolution windows are defined over the model time index 
                   ! intervals [ kount_oa(1,:), kount_oa(2,:) ]
                   ! Being given the width of the convolution window 2*dkount_tot_t
                   ! and its time index center k and being given the simulation time window :
                   ! - the 1st possible convolution window can start at the 1st model time index
                   ! which is the 'now' time index of the restart just before the 1st 
                   ! model time step from now to after, ie kountv_oa(1, izvt_oa=1 ) = kount0
                   ! (see call init_oa followed by main_oa at k=kount0)
                   ! - the last time index of the last possible convolution window is allowed to 
                   !   match the 'after' index of the last time step of the model, i.e.,
                   !   kountv_oa( 2, izvt_oa=nzvt_oa ) = nt_max
                   ! => Excluded model time steps (flag -9999) should be 
                   !    therefore defined as follows (see var_rec_oa) 
                   ! ORIGINAL TESTS
                   !if (    kountv_oa(1,nzvt_oa).lt.kount0             &
                   !    .or.kountv_oa(1,nzvt_oa).gt.nt_max-1               &
                   !    .or.kountv_oa(2,nzvt_oa).lt.kount0             &
                   !    .or.kountv_oa(2,nzvt_oa).gt.nt_max-1 )             &
                   !    then
                   convwind_out_of_simulation_1 : if (    kountv_oa(1,nzvt_oa).lt.kount0  & 
                       .or.kountv_oa(1,nzvt_oa).gt.nt_max                                 &
                       .or.kountv_oa(2,nzvt_oa).lt.kount0                                 &
                       .or.kountv_oa(2,nzvt_oa).gt.nt_max )                               &
                       then


                     kountv_oa(1,nzvt_oa) = -9999
                     kountv_oa(2,nzvt_oa) = -9999

                     ! Cas des conf/var avec +s periodes discretes
                     ! les pointeurs temp. inlcus la localisation temporelle de l'analyse au cours de la simu
                     ! avec une correspondance biunivoque avec les pointeurs de periode
                     ! Les parametres ci-dessous permettent seulement de reperer le 1er rec pour la 1e periode
                     ! #BLXD testing tracking first/last conv. window
                     if ( ( if_first_rec_oa(iv_t)==1 ) .and. & 
                          ( if_last_rec_oa(iv_t)==0 ) ) then
                      ltrec_lst_oa(iv_t) = nzvt_oa-1   ! Periode discrete de 900s a 1800s le 1er lt qui ne passera plus est + grde per 1800s
                      if_last_rec_oa(iv_t) = 1         ! Periode integre iv_t plage de ip_t mais seul la plus grande est explore soit 1800s
                     endif

                     if ( ( if_first_per_oa(iv_t,ip_t)==1 ) .and. & 
                          ( if_last_per_oa (iv_t,ip_t)==0 ) ) then
                      lper_lst_oa   (iv_t,ip_t) = nzvt_oa-1 ! Periode integre iv_t plage de ip_t mais seul la plus grande est explore soit 1800s
                      if_last_per_oa(iv_t,ip_t) = 1         ! Periode discrete correspondance totale
                     endif

                   else convwind_out_of_simulation_1
!----------------------------------------------------------------

                     ! #BLXD testing tracking first/last conv. window
                     if ( if_first_rec_oa(iv_t)==0 ) then
                      ltrec_fst_oa(iv_t) = nzvt_oa   ! Periode discrete de 900s a 1800s => 1er lt de la plus petite periode  soit 900s
                      if_first_rec_oa(iv_t) = 1      ! Periode integre iv_t plage de ip_t mais seul la plus grande est explore soit 1800s
                     end if 

                     if ( if_first_per_oa(iv_t,ip_t)==0 ) then
                      lper_fst_oa   (iv_t,ip_t) = nzvt_oa ! Periode integre iv_t plage de ip_t mais seul la plus grande est explore soit 1800s
                      if_first_per_oa(iv_t,ip_t) = 1      ! Periode discrete correspondance totale
                     endif

                     ! #BLXD What is the purpose of updv_oa == 1 ?
                     ! if the above condition (if convwind_out_of_simulation_1 
                     ! is false non of the ORIGINAL tests
                     ! after updv_oa can be found true

                     if (updv_oa(iv_t).eq.1) then
                      if (kountv_oa(1,nzvt_oa).lt.kount0) then
                       kountv_oa(1,nzvt_oa) = kount0
                       kountv_oa(2,nzvt_oa) = kount0 + 2*dkount_tot_t
                      endif
                      if (kountv_oa(2,nzvt_oa).gt.nt_max-1) then
                       kountv_oa(1,nzvt_oa) = nt_max     - 2*dkount_tot_t
                       kountv_oa(2,nzvt_oa) = nt_max
                      endif
                     endif

                     comptt(iv_t)=comptt(iv_t)+1
                   endif convwind_out_of_simulation_1

                  endif ! second pass
               enddo simulation_window
            endif not_fourier


!----------------------------------------------------------------
!     configuration temporelle numero 4   : simu
!     
!     atome: fourier full simu time
!----------------------------------------------------------------
            ! BLXD I would prefer
            ! if ((swt_t_oa(iv_t).eq.1.or.swt_t_oa(iv_t).eq.4).and.tpsi_oa(iv_t).eq.3 ) then
            ! OK with init_parameter_oa : tpsi_oa == 3 <=> swt_t_oa(iv_t).eq.4)
            ! CHECK kountv_oa(1,)=ks-1, kountv_oa(2,)=ks-1
            fourier : if ( swt_t_oa(iv_t).eq.4 .and. tpsi_oa(iv_t).eq.3 ) then
               nzvt_oa             = nzvt_oa + 1

               if (flag_t) then
                kountv_oa(1,nzvt_oa) = nt_max-1                       &
                    - int ( int( real(nt_max-1-kount0+1) * dti        &
                           / perv_oa(1,ip_t) )                        &
                         * perv_oa(1,ip_t)/dti ) + 1             
                kountv_oa(2,nzvt_oa) = nt_max-1  
                per_t2p_oa (nzvt_oa) = ip_t
                ! One single Fourier analysis using all the simulation window for integration  
                comptt(iv_t)       = 1
               endif

            endif fourier


!----------------------------------------------------------------
!     configuration temporelle numero 2: une date a la fin de la simulation
!     0 : dirac, 1 : ondelette, 2 : windowed Fourier
!----------------------------------------------------------------
            if (swt_t_oa(iv_t).eq.2) then
               nzvt_oa             = nzvt_oa + 1
               if (flag_t) then
                kountv_oa(1,nzvt_oa) = nt_max - 1 - 2*dkount_tot_t
                kountv_oa(2,nzvt_oa) = nt_max - 1      
! BLXD TODO test with change : nt_max - 1  to nt_max
                !kountv_oa(1,nzvt_oa) = nt_max - 2*dkount_tot_t
                !kountv_oa(2,nzvt_oa) = nt_max      
                per_t2p_oa (nzvt_oa)  = ip_t
! BLXD TODO change nt_max - 1  to nt_max in the if instruction as well
                if (    kountv_oa(1,nzvt_oa).lt.kount0                 &
                    .or.kountv_oa(1,nzvt_oa).gt.nt_max-1               &
                    .or.kountv_oa(2,nzvt_oa).lt.kount0                 &
                    .or.kountv_oa(2,nzvt_oa).gt.nt_max-1 )             &
                    then
                  kountv_oa(1,nzvt_oa) = -9999
                  kountv_oa(2,nzvt_oa) = -9999

                else
! BLXD one variable a single date of conv. window at a single date

                   comptt(iv_t)       = 1
                endif
               endif
            endif

!----------------------------------------------------------------
!     configuration temporelle numero 3: date unique
!----------------------------------------------------------------
            if (swt_t_oa(iv_t).eq.3) then
               nzvt_oa             = nzvt_oa + 1

               if (flag_t) then
! BLXD TODO change nt_max - 1  to nt_max in the if instruction as well
               kountv_oa(1,nzvt_oa) = kount_user_oa(1,iv_t) - dkount_tot_t
               kountv_oa(2,nzvt_oa) = kount_user_oa(1,iv_t) + dkount_tot_t
               per_t2p_oa (nzvt_oa) = ip_t
               if (    kountv_oa(1,nzvt_oa).lt.kount0                  &
                    .or.kountv_oa(1,nzvt_oa).gt.nt_max - 1                 &
                    .or.kountv_oa(2,nzvt_oa).lt.kount0                 &
                    .or.kountv_oa(1,nzvt_oa).gt.nt_max - 1 )               &
                    then
                  kountv_oa(1,nzvt_oa) = -9999
                  kountv_oa(2,nzvt_oa) = -9999
               else
! BLXD one variable a single date of conv. window at a single date
                   comptt(iv_t)       = 1
               endif
               endif
            endif

!.....enddo ip_t

            if (flag_t) ncomptt   = ncomptt + comptt(iv_t)

         enddo loop_on_requested_period

         
!---- >dernier point:
         if (flag_t) begvt_oa( nzv_oa + 1 )   = nzvt_oa + 1

!.....enddo: boucle sur les variables...

      enddo loop_on_variables


      return
      end subroutine var_time_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note This routine is called under a Croco tile(-thread) loop and requires
!!  Preventing double allocation issues and openMP data race
!!  (all threads updating implicitly SHARED global module parameters),
!!  and Avoidingd multiple initialization of OA module parameters (tile-only simu),
!!  Only the last tile(-threads) calls allocate_part4_oa (with openMP critical construct)
!
! DESCRIPTION: allocation de variables globales definissant la structure temporelle
!              du vecteur d'etat.
!
!> @brief allocates global variables defining the "state vector" time structure.
!
! REVISION HISTORY:
!
!> @authors B. Lemieux-Dudon
!!   Adapted from Symphonie/NHOMS initial version (2006)
!> @date 2021
!> @todo change sub. name part4/time
!------------------------------------------------------------------------------

      subroutine allocate_part4_oa

      use module_oa_time
      use module_oa_periode
      !use module_oa_stock

      implicit none
      
      !> calcul des kounts de debut et de fin pour chaque variable (ou configuration)
      if ( .not. allocated(kountv_oa) )                               &
      allocate (                                                      &
       kountv_oa (2,nzvt_oa)                                          &
       ) 
     
      if ( .not. allocated(tallocated_oa) )                           &
      allocate (                                                      & 
       tallocated_oa(nzvt_oa)                                         &
       )

      if ( .not. allocated(per_t2p_oa) )                              &
      allocate (                                                      &
       per_t2p_oa (nzvt_oa)                                           &  ! transformation structure temporelle --> structure frequentielle
       )

      return
      end subroutine allocate_part4_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note This routine is called under a Croco tile(-thread) loop and requires
!!  Preventing double allocation issues and openMP data race
!!  (all threads updating implicitly SHARED global module parameters),
!!  and Avoidingd multiple initialization of OA module parameters (tile-only simu),
!!  only the last tile(-threads) calls allocate_box_oa (with openMP critical construct)
!
! DESCRIPTION: allocation de variables globale pour le calcul des incertitudes d'Heisenberg
!
!> @brief allocates global variables to calculate Heisenberg boxes.
!
! REVISION HISTORY:
!
!> @authors B. Lemieux-Dudon
!> @date 2021
!> @todo 
!------------------------------------------------------------------------------

      subroutine allocate_box_oa

      use module_oa_periode
      implicit none

      if ( .not. allocated(t0) )  allocate ( t0(nzvp_oa)) 
      if ( .not. allocated(dt0) ) allocate ( dt0(nzvp_oa)) 
      if ( .not. allocated(psi_norm_l2) ) allocate ( psi_norm_l2 (nzvp_oa) ) 

      if ( .not. allocated(w0) )  allocate ( w0(nzvp_oa)) 
      if ( .not. allocated(dw0) ) allocate ( dw0(nzvp_oa)) 

      if ( .not. allocated(ft_psi_norm_l2) ) allocate ( ft_psi_norm_l2 (nzvp_oa) ) 

      if ( .not. allocated(t0_bis) )  allocate ( t0_bis(nzvp_oa)) 
      if ( .not. allocated(dt0_bis) ) allocate ( dt0_bis(nzvp_oa)) 
      if ( .not. allocated(psi_norm_l2_bis) ) allocate ( psi_norm_l2_bis (nzvp_oa) ) 

      ! temporary used in count_nmsimult_oa NO BIS VALUE NEEDED 
      if ( .not. allocated(sq_dt0) ) allocate( sq_dt0(nzvp_oa) )
      if ( .not. allocated(sq_dw0) ) allocate( sq_dw0(nzvp_oa) )

      ! w0_apx, dw0_apx : approx heisenberg boxes with dt0
      if ( .not. allocated(w0_apx)  )      allocate ( w0_apx (nzvp_oa) ) 
      if ( .not. allocated(dw0_apx) )      allocate ( dw0_apx (nzvp_oa)) 

      ! #BLXD TODO bis quantities calculated in var_rec_oa
      !       TO BE INITIALIZED TO HUGE VALUES TO TRACK ERRORS
      !       OR TO SPECIFIC IMPOSSIBLE < 0 VALUES
      !       ESPECIALLY psi_norm_bis and dt0_bis 

      end subroutine allocate_box_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note This routine is called under a Croco tile(-thread) loop and requires
!!  to prevent double allocation issues and openMP data race
!!  (all threads updating implicitly SHARED global module parameters),
!!  and to avoid multiple initialization of OA module parameters (tile-only simu),
!!  only the last tile(-threads) calls allocate_theo_box_oa 
!!  (with openMP critical construct)
!
!
! DESCRIPTION: allocation/desallocation de variables globale pour le calcul theorique 
!              des incertitudes d'Heisenberg
!
!> @brief allocates/desallocates Heisenberg boxes global variables.
!
! REVISION HISTORY:
!
!> @authors B. Lemieux-Dudon
!> @date 2021
!> @todo
!------------------------------------------------------------------------------

      subroutine allocate_theo_box_oa

      use module_oa_periode
      implicit none

      if ( .not. allocated(t0_theo) ) allocate ( t0_theo (nzvp_oa)) 
      ! See functions to calculate on the fly the theoritical values
      !allocate ( dt0_theo (nzvp_oa)) 
      !allocate ( w0_theo (nzvp_oa), dw0_theo (nzvp_oa)) 

      end subroutine allocate_theo_box_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note This routine is called under a Croco tile(-thread) loop and requires
!!  to prevent double allocation issues and openMP data race
!!  (all threads updating implicitly SHARED global module parameters),
!!  and to avoid multiple initialization of OA module parameters (tile-only simu),
!!  only the last tile(-threads) calls calc_temporal_box_oa_oa 
!!  (using openMP critical construct)
!
!
! DESCRIPTION: estimation itérative des incertitudes d'Heisenberg numérique qqsoit l'atome
!
!> @brief evaluates numerical heisenberg uncertainties (any atom).
!> @details applies to any atoms by iterative calls to cumulate spectral power.
!
! REVISION HISTORY:
!
!> @authors B. Lemieux-Dudon
!> @date 2021
!> @todo
!------------------------------------------------------------------------------

      subroutine calc_temporal_box_oa( ndim, psi_norm, t0_psi, sq_dt0_psi, kpt, dti, norm_ndim, norm, ind ) 

      implicit none

      real*8, dimension(1:ndim), intent(inout) :: t0_psi, sq_dt0_psi
      real,   dimension(1:ndim), intent(inout) :: psi_norm
      real,   dimension(1:ndim), intent(in), optional :: norm_ndim
      real,                      intent(in), optional :: norm
      real*8,                    intent(in)    :: dti
      integer,                   intent(in)    :: kpt
      integer,                intent(in)    :: ndim
      integer,optional,       intent(in)    :: ind
      real    :: norm_r
      real*8  :: time
      integer :: ii,is,ie

      if ( present(ind) ) then
        if (ind<=ndim) then
            is = ind
            ie = ind
            if ( .not. present(norm) ) then
                write(*,*)'ERROR OA : calc_temporal_box_oa requires norm'
                if(if_print_node) write(io_unit,*)'ERROR OA : calc_temporal_box_oa requires norm'
                stop
            else 
                norm_r = norm
            endif
        else
            write(*,*)'ERROR OA : calc_temporal_box_oa max index is ', ndim
            if(if_print_node) write(io_unit,*)'ERROR OA : calc_temporal_box_oa max index is ', ndim
            stop
        endif
      else 
        is = 1
        ie = ndim
        if ( .not. present(norm_ndim) ) then
            write(*,*)'ERROR OA : calc_temporal_box_oa requires norm_ndim'
            if(if_print_node) write(io_unit,*)'ERROR OA : calc_temporal_box_oa requires norm_ndim'
            stop
        endif
      endif

    
      time = dble(kpt) * dti ! BLXD require precision module
      do ii=is,ie

        if ( present(norm_ndim) ) norm_r = norm_ndim(ii)
        !if (ii==2) then
        !    if(if_print_node) write(io_unit,*) time, norm_r
        !endif 

        !Psi L2-norm BLXD Note that there is no dti factor for integration
        psi_norm(ii)    = psi_norm(ii)     +  norm_r

        !Average t0 weighed by |psi_m|^2 BLXD Note that there is no dti factor for integration
        !consistently with other OA temporal integration 
        t0_psi(ii)      = t0_psi(ii)       +  norm_r * time 

        !Quadratic average for t0^2 weighed by |psi_m|^2 BLXD no dti factor for integration
        sq_dt0_psi(ii)  = sq_dt0_psi(ii)   +  norm_r * time**2
      enddo

      end subroutine calc_temporal_box_oa

      subroutine init_temporal_box_oa( ndim, psi_norm, t0_psi, sq_dt0_psi, std_dt0_psi ) 
                     
      implicit none

      real*8, dimension(1:ndim), intent(inout) :: t0_psi, sq_dt0_psi, std_dt0_psi
      real,   dimension(1:ndim), intent(inout) :: psi_norm
      integer,                   intent(in)    :: ndim

      psi_norm   (1:ndim) = 0.
      t0_psi     (1:ndim) = 0.d0
      sq_dt0_psi  (1:ndim) = 0.d0
      std_dt0_psi (1:ndim) = 0.d0

      end subroutine init_temporal_box_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note This routine is called under a Croco tile(-thread) loop and requires
!!  to prevent double allocation issues and openMP data race
!!  (all threads updating implicitly SHARED global module parameters),
!!  and to avoid multiple initialization of OA module parameters (tile-only simu),
!!  only the last tile(-threads) calls finalize_box_oa (with openMP critical construct)
!
!
! DESCRIPTION: finalisation du calcul itératif des incertitudes d'Heisenberg numérique.
!
!> @brief finalizes numerical heisenberg uncertainties with traces.
!
! REVISION HISTORY:
!
!> @authors B. Lemieux-Dudon
!> @date 2021
!> @todo
!------------------------------------------------------------------------------

      subroutine finalize_box_oa( ndim, psi_norm, psi_norm_tgt       &
                                      , t0_psi,   t0_psi_tgt         &
                                      , sq_dt0_psi, sq_dt0_psi_tgt   &
                                      , std_dt0_psi, std_dt0_psi_tgt &
                                      , myprec, myprec2, ind ) 

      implicit none

      real*8, dimension(1:ndim), intent(inout) :: t0_psi, sq_dt0_psi
      real,   dimension(1:ndim), intent(inout) :: psi_norm
      real*8, dimension(1:ndim), intent(inout) :: std_dt0_psi
      real*8,                    intent(in) :: t0_psi_tgt, sq_dt0_psi_tgt, std_dt0_psi_tgt
      real,                      intent(in) :: psi_norm_tgt
      real*8,                 intent(in)    :: myprec2
      real,                   intent(in)    :: myprec
      integer,                intent(in)    :: ndim
      integer,optional,       intent(in)    :: ind
      integer :: ii,is,ie
      logical :: ll_io_hbg
      !character(len=30), parameter :: fmt0="(1x,ES22.15E2)"
      !character(len=30), parameter :: fmt1="(1x,f22.6)"
      !character(len=*) :: myfmt1, myfmt2

      real   :: delta
      real*8 :: delta2

      if (mynode/=0) return

      if ( present(ind) ) then
        if (ind<=ndim) then
            is = ind
            ie = ind
        else
            if(if_print_node) write(io_hbg,*) 'ERROR OA : finalize_box_oa max index is ', ndim
            stop
        endif
      else 
        is = 1
        ie = ndim
      endif

      if (mynode==0) inquire(unit=io_hbg, opened=ll_io_hbg) 
      if ( .not. ll_io_hbg ) then
        print*,'ERROR OA : cannot write in heisenberg box file unit closed'
        stop
      endif
    
      do ii=is,ie,1

        if(if_print_node) then
           write(io_hbg,fmt="(a,i3)")   '    ** OA check Heisenberg box for evaluation # ',ii
           write(io_hbg,fmt="(a,2(f22.6,x))") '      => OA psi_oa L2-norm and target value',psi_norm(ii), psi_norm_tgt
        endif
        if ( abs( psi_norm(ii) - psi_norm_tgt ) > myprec )  then
            if(if_print_node) then
            write(io_hbg,fmt='(a,2(f22.6,x))') '         OA WARNING : psi_oa L2-norm does not meet the target value ',psi_norm(ii), psi_norm_tgt
            write(io_hbg,fmt='(a,1(f22.6))') '                      above threshold ',myprec
            endif
        endif

        if(if_print_node) write(io_hbg,fmt="(a,i3)")   '    ** OA Heisenberg box normalization using psi_oa L2-norm'

        t0_psi(ii)     = t0_psi(ii)/psi_norm(ii)
        sq_dt0_psi(ii) = sq_dt0_psi(ii)/psi_norm(ii)

        if(if_print_node) write(io_hbg,fmt="(a,2(1x,ES22.15E2))") '      => OA box temporal center and target value',t0_psi(ii), t0_psi_tgt
        if ( abs( t0_psi(ii) - t0_psi_tgt ) > myprec2 )  then
        if(if_print_node) then
           write(io_hbg,fmt='(a,2(1x,ES22.15E2))') '        OA WARNING : box temp. value does not meet the target value ',t0_psi(ii), t0_psi_tgt
           write(io_hbg,fmt='(a,1(f22.6))') '                     above threshold ',myprec2
        endif
        endif

        if(if_print_node) write(io_hbg,fmt="(a,2(1x,ES22.15E2))") '      => OA box temporal variance and target value',sq_dt0_psi(ii), sq_dt0_psi_tgt
        if ( abs( sq_dt0_psi(ii) - sq_dt0_psi_tgt ) > myprec2 )  then
            if(if_print_node) then
            write(io_hbg,fmt='(a,2(1x,ES22.15E2))') '        OA WARNING : box temp. value does not meet the target value ',sq_dt0_psi(ii), sq_dt0_psi_tgt
            write(io_hbg,fmt='(a,1(f22.6))') '                     above threshold ',myprec2
            endif
        endif

        std_dt0_psi(ii) = sqrt(sq_dt0_psi(ii))

        if(if_print_node) write(io_hbg,fmt="(a,2(1x,ES22.15E2))") '      => OA box temporal stdev and target value',std_dt0_psi(ii), std_dt0_psi_tgt
        if ( abs( std_dt0_psi(ii) - std_dt0_psi_tgt ) > myprec2 )  then
           if(if_print_node) then
           write(io_hbg,fmt='(a,2(1x,ES22.15E2))') '        OA WARNING : box temp. value does not meet the target value ',std_dt0_psi(ii), std_dt0_psi_tgt
           write(io_hbg,fmt='(a,1(f22.6))') '                     above threshold ',myprec2
           endif
        endif

      enddo

      ! BLXD only true in calls within var_rec_oa
      if ( (is==1) .and. (ie==2) ) then

        delta =  100. * 2. * ( psi_norm(2) - psi_norm(1) ) / ( psi_norm(2) + psi_norm(1) )

        if ( abs(delta) > myprec ) then
         if(if_print_node) then
         write(io_hbg,fmt='(a,1(f22.6))') '  OA WARNING : Psi L2-norm perc. ratio between the 2 eval. above threshold ',myprec
         write(io_hbg,fmt='(a,1(f22.6))') '            delta ',delta
         write(io_hbg,fmt='(a,1(f22.6))') '            norm1 ',psi_norm(1)
         write(io_hbg,fmt='(a,1(f22.6))') '            norm2 ',psi_norm(2)
         endif
         !stop
        end if

        delta2 =  100.d0 * 2.d0 * ( t0_psi(2) - t0_psi(1) ) / ( t0_psi(2) + t0_psi(1) )

        if ( abs(delta2) > myprec2 ) then
         if(if_print_node) then
         write(io_hbg,fmt='(a,1(f22.6))') '  OA WARNING : temp. box center perc. ratio between the 2 eval. above threshold ',myprec2
         write(io_hbg,fmt='(a,1(f22.6))') '            delta ',delta2
         write(io_hbg,fmt='(a,1(1x,ES22.15E2))') '            norm1 ',t0_psi(1)
         write(io_hbg,fmt='(a,1(1x,ES22.15E2))') '            norm2 ',t0_psi(2)
         endif
         !stop
        end if

        delta2 =  100.d0 * 2.d0 * ( std_dt0_psi(2) - std_dt0_psi(1) ) / ( std_dt0_psi(2) + std_dt0_psi(1) )

        if ( abs(delta2) > myprec2 ) then
         if(if_print_node) then
         write(io_hbg,fmt='(a,1(f22.6))') '  OA WARNING : temp. box stdev perc. ratio between the 2 eval. above threshold ',myprec2
         write(io_hbg,fmt='(a,1(f22.6))') '            delta ',delta2
         write(io_hbg,fmt='(a,1(1x,ES22.15E2))') '            norm1 ',std_dt0_psi(1)
         write(io_hbg,fmt='(a,1(1x,ES22.15E2))') '            norm2 ',std_dt0_psi(2)
         endif
         stop
        end if

      endif

      end subroutine finalize_box_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note This routine is called under a Croco tile(-thread) loop and requires
!!  to prevent double allocation issues and openMP data race
!!  (all threads updating implicitly SHARED global module parameters),
!!  and to avoid multiple initialization of OA module parameters (tile-only simu),
!!  only the last tile(-threads) calls unset_convwind_oa 
!!  (with openMP critical construct)
!
! DESCRIPTION: flag identifiant les fenetres de convolution temp. en cours de calcul 
!
!> @brief flag identifying convolution windows currently under calculation.
!
! REVISION HISTORY:
!
!> @authors B. Lemieux-Dudon
!> @date 2021
!> @todo
!------------------------------------------------------------------------------

      subroutine unset_convwind_oa

      use module_oa_time, only : tallocated_oa

      implicit none

      tallocated_oa(:)   = -1
      return
      end subroutine unset_convwind_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note croco tile(-thread) loop, all threads updating/allocating shared global module parameters 
!! => data race and double allocation issues.
!
! DESCRIPTION: allocation du vecteur d'etat w contenant les analyses.
!
!> @brief allocates global variables defining the "state vector" time structure.
!
! REVISION HISTORY:
!
!> @authors B. Lemieux-Dudon
!!   reducing allocation size to the max # of simultaneously openned convolution windows
!!   pre-calculated at the begining of the OA simulation (see, nmsimult_oa) 
!> @date 2021
!> @todo OBSOLOLETE REMOVE BLXD_TEST_TILE_ISSUE 
!------------------------------------------------------------------------------
      subroutine allocate_part5_oa(tile)
      !- st(tile)%wf_oa(l_a)%coef(dim_a) PRIVATE dim_a per tile (requires tile args in)
      !- rm pst from list of public module variables in module_tile_oa add tile_space_str
      use module_tile_oa, only : st, tile_space_str
      use module_oa_time
      implicit none
      integer, intent(in) :: tile

      type(tile_space_str), pointer :: pst => null()

      integer :: l_a

      pst => st(tile)
      if ( .not. allocated(pst%wf_oa) ) then
          allocate (                                                      &
           pst%wf_oa (nmsimult_oa)                                            &
           )
      
          do l_a=1,nmsimult_oa
             nullify(pst%wf_oa(l_a)%coef)
          enddo
      
          tallocated_oa(:)   = -1
          pst%wf_oa(:)%t_indice  = -1
          pst%wf_oa(:)%config    = -1
          pst%wf_oa(:)%variable  = -1
      endif
      pst => null()

      return
      end subroutine allocate_part5_oa

!----------------------------------------------------------------------
! PROCEDURE
!
!> @note This routine is called under a Croco tile(-thread) loop and requires
!!  Preventing double allocation issues and openMP data race
!!  (all threads updating implicitly SHARED global module parameters),
!!  and Avoidingd multiple initialization of OA module parameters (tile-only simu),
!!  only the last tile(-threads) calls user_count_nmsimult_oa 
!!  (using openMP critical construct)
!
! DESCRIPTION: 
! Optimization de la memoire pour allouer la structure temporelle
! du vecteur d'etat a sa juste taille qui correspond au nombre
! maximum de fenetres de convolution ouvertes a un meme pas de temps 
! de simulation.
! 
!> @brief memory optimization pre-calculating the size of 
!   temporal structure of the state vector (wf_oa dimension) 
!   using the exact same algo applied in main_oa.
!
!> @details nmsimult_oa based on the number of requested 
!!  config, variable, and frequency of analysis during the simulation
!   vectorize the analysis.
!!  The state vector is an derived data type array whose dimension
!!  have fields having dimension equal to the # of grid points 
!!  is the # of convolution windows (time structure), and which 
!!  to analyse (spatial structure). 
!!
! REVISION HISTORY:
!
!> @authors
!! - B. Lemieux-Dudon
!!  - Croco-OA interface 1st version (2020)
!!  - Croco Tile-thread compliant version (2021)
!!  - Memory optimization by reducing the state vector allocated size by
!!    pre-calculating the maximum number of temporal convolution windows
!!    simoultaneously opened. 
!!  - 2nd Heisenberg numerical uncertainty calculation.
!
!> @date 2021
!> @todo BLXD
!! - declare psi as real(KIND=wp) with specific module handl. prec. 
!! - revise as the algo for main_oa subroutine evolves, e.g.
!!   if convolution windows are simulatenously defined over
!!   the slow and fast modes.
!! - seperate Heisenberg boxes from nmsimult_oa calculation ?
!! - Heisenberg numerical boxes, still ongoing tests
!----------------------------------------------------------------------

      subroutine user_count_nmsimult_oa( ) 

      use module_tile_oa, only : st, ntiles
      use module_oa_time

      implicit none

      !> Tile count
      integer :: itile

      nmsimult_oa = nmsimult_oa_max

      do itile=1,ntiles
        st(itile)%nzw_oa = nmsimult_oa
      enddo

      return
      end subroutine user_count_nmsimult_oa

      subroutine count_nmsimult_oa( ichoix                 &
                          ,ivar_m                          &
                          ,kount0                          &
                          ,nt_max                          &
                          ,dti                             &
                          ,tile                            ) 

      use module_tile_oa, only : st, tile_space_str, ntiles
      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
!      use module_oa_stock

      implicit none

      !> Tile dependent calculation, ichoix set to zero in "regular" cases
      integer, intent(in) :: tile, ichoix

      !> To handle specific calls to main_oa with variables updated at different position in the calling code
      integer, intent(in), optional :: ivar_m

      !> Time integration step
      double precision, intent(in) :: dti

      !> First and last simulation iteration index
      integer, intent(in) :: nt_max, kount0                                        

      type(tile_space_str), pointer :: pst => null()

      !TODO declare flag here if only needed here
      logical :: ifl_test_composite

      integer ::                                                      &
         iv_m                                                         &
        ,ic_m                                                         &
        ,lt_m                                                         &
        ,ls_m                                                         &
        ,ls1_m

      integer ::                                                      &
         kpt_m                                                        &
        ,lp_m                                                         &
        ,lp1_m                                                        &
        ,lp2_m                                                        &
        ,lpreso_m                                                     &
        ,lpconw_m 

      !> For the model integration iteration
      integer :: iic_oa 

      !> Count of total size of the stock structure / process IO
      integer :: imsimult_oa, io_nodoa

      !> Tile count
      integer :: itile

      !BLXD not useful for now 
      ! For checking the count of spatial structure size coef
      ! might set the logical to false if to mudch time consuling at init time
      !logical, parameter :: if_count_space=.true. 
      !WARNING tile-(threads) dependent variables => loop over tile(-threads)
      !required (not only last tile)
      !integer, dimension(1:ntiles,1:nzc_oa,1:nzv_oa,1:nzvt_oa) ::     &
      !   intper_count_space                                           &
      !  ,sglper_count_space
      ! integer ::                                                     &
      ! i_m                                                            &
      !,j_m                                                            &
      !,k_m                                                            &

      integer, dimension(1:nzvt_oa) :: kount_tot
      integer                       :: kount_tot_low, kount_tot_up, kount_tot_hlf

! #BLXD : double precision ? TODO precision module required
      complex                                                         &
         psi_m

      real,   dimension(1:2) :: norm_r
      real                   :: norm
      real*8                 :: ut_m, std_dt0_theo

! #BLXD TODO precision module required 
      real,   parameter :: myprec  = 10.e-5
      real*8, parameter :: myprec2 = 10.e-10

! #BLXD TODO test .false. => if identical eval  use .false. 
      logical, parameter :: from_cplx_fun = .true.

#ifdef OA_TRACES
      integer, dimension(1:2), parameter :: io_trc=(100,200)
#endif

!.... The size of the module_oa_stock structure wf_oa is only dependent of
!     conf ic_m, var iv_m , and specific conv. window lt_m
!     It is independent of the tile which is related to the 2D horizontal domain
!     "tile-decomposition"
      imsimult_oa = 0 
      nmsimult_oa = 0 

!.... For each conf ic_m, var iv_m , and specific conv. window lt_m
!     One can check the # of analysed spatial points i,j,k 
!     Tile-(threads) dependent var.
!      sglper_count_space(tile, :,:,:) = 0
!      intper_count_space(tile, :,:,:) = 0

      oa_analysis_requested : if (nzv_oa.ne.0) then

!.... Pointer to tile structure
!     if (if_count_space) then
!      pst => st(tile)
!     endif

!.....flag variable composite 
      ifl_test_composite = .true.      

      loop_timestepping : do iic_oa = kount0, nt_max

      loop_ic : do ic_m = 1 , nzc_oa
      loop_iv : do iv_m = begc_oa(ic_m),begc_oa(ic_m+1)-1

!........ce test est utilise pour les appels cibles lors des calculs energetiques:
!        STDALONE to call main_oa at specific location for composite variable confguration
!        if (present(ivar_m)) then
!           ifl_test_composite = (tv_oa(iv_m).eq.ivar_m)
!        endif
!         if ( (ichoix.ge.0.or.ifl_test_composite) .and. ( tv_oa(iv_m).ne.20) ) then

!---------------------------------------------------------------------------
!.....periodes discretisees:
!---------------------------------------------------------------------------

     if_discrete_per : if (dori_oa(iv_m).eq.1) then

!....Morlet => Heisenberg box 
!....#BLXD TODO calculation for all requested period => possible optimization
!      if (tpsi_oa(iv_m)==1) then 
!        call init_temporal_box_oa( nzvp_oa, psi_norm_l2, t0, sq_dt0, dt0 ) 
!      end if

      time_of_analysis_loop : do lt_m  = begvt_oa(iv_m) , begvt_oa(iv_m+1) - 1

!.....for a given iv_m => set of lp_m pointers with bijective corr. to a set of lt_m pointers
         lp_m = per_t2p_oa(lt_m)

!.....test pour le calcul:
       tstp_match_convol_window : if ( iic_oa.ge.kountv_oa(1,lt_m) .and. iic_oa.le.kountv_oa(2,lt_m)                  &
            .and.               mod( int( (iic_oa-kountv_oa(1,lt_m)), kind=8 ), resv_oa(per_t2p_oa(lt_m)) ).eq.0      &
            .and.               ( ichoix.eq.cnb_oa(iv_m) .or. ichoix.lt.0 )                                           &
                             ) then

! #BLXD ( ichoix.eq.cnb_oa(iv_m) .or. ichoix.lt.0 ) = calling main_oa for particular variable


!........allocation if necessary (if variable isn't allocated yet, it is done.., at user-specified kount):


         if (tallocated_oa(lt_m).eq.-1) then
           imsimult_oa = imsimult_oa + 1
           nmsimult_oa = max(nmsimult_oa,imsimult_oa)
           !tallocated_oa(lt_m) = imsimult_oa
           tallocated_oa(lt_m) = lt_m 
           kount_tot(lt_m) = 0
         else 
           kount_tot(lt_m) = kount_tot(lt_m) + 1
         endif

!........time translation applied to the atom (s)
!        ut_m  = real((kountv_oa(1,lt_m)+kountv_oa(2,lt_m))/2.)*dti

!........precalcul de l'indice temporel prenant en compte la translation (u):
!
         kpt_m =  iic_oa - ( int((kountv_oa(1,lt_m)+kountv_oa(2,lt_m))/2) )

!........precalcul de psi:

         psi_m = psi_oa(                                                        &
                    tpsi_oa(iv_m)                                               &
                   ,psi_p2s_oa(                                                 &
                                     tpsi_oa(iv_m)                              &         
                                    ,perv_oa(1,per_t2p_oa(lt_m)),fb_oa,fc_oa )  &                                     
                   ,real(kpt_m)* dti                                            &
                   ,dti*resv_oa(per_t2p_oa(lt_m))                               &           
                   ,fb_oa,fc_oa  ) 

!........Heisenberg box
         if_hsbg_dp1 : if ( if_heisenberg_box ) then

!........BLXD TODO Heisenberg box only dev. for Morlet wavelet
            if_mrlt_psi_dp1 : if (tpsi_oa(iv_m)==1) then

!...........Psi_oa L2-norm 
            norm_r(1) =  real( conjg(psi_m) * psi_m )

            norm_r(2) = module_morlet_psi_oa (                                  &
                    tpsi_oa(iv_m)                                               &
                   ,psi_p2s_oa(                                                 &
                                     tpsi_oa(iv_m)                              &         
                                    ,perv_oa(1,per_t2p_oa(lt_m)),fb_oa,fc_oa )  &                                     
                   ,real(kpt_m)* dti                                            &
                   ,dti*resv_oa(per_t2p_oa(lt_m))                               &           
                   ,fb_oa, fc_oa, from_cplx_fun  )  

            if ( abs( norm_r(1) - norm_r(2) ) > myprec2 ) then
                if(if_print_node) then
                 write(io_unit,*) 'ERROR OA : count_nmsimult_oa psi norms are too different ', norm_r(1), norm_r(2)
                 write(io_unit,*) 'ERROR OA : SHOULD STOP module_morlet_psi_oa gives bad results !!!!'
                endif
                !stop
            else 
                norm = norm_r(1)
            end if

            ! BLXD args are full size array with dimension nzvp_oa, but processing only targets lp_m
            call calc_temporal_box_oa( ndim=nzvp_oa                 &
                                      ,psi_norm=psi_norm_l2         &
                                      ,t0_psi=t0                    &
                                      ,sq_dt0_psi=sq_dt0            &
                                      ,kpt=kpt_m                    &
                                      ,dti=dti                      &
                                      ,norm=norm                    &
                                      ,ind=lp_m ) 
            endif if_mrlt_psi_dp1
         endif if_hsbg_dp1


#ifdef OA_TRACES
         !if ( ( lt_m==ltrec_fst_oa(iv_m) ) ) then
         if ( ( lt_m==lper_fst_oa(iv_m,lp_m) ) ) then
         if (lp_m <= 99) then
            if ( (iv_m==1) .or. (iv_m==2) ) then
             io_nodoa = (io_trc(iv_m)+lp_m)*1000+mynode  !io_nodoa = 10000+mynode now 10X000+mynode
             write (io_nodoa,fmt='(i4,3(1x,ES22.15E2))') kpt_m, REAL(DBLE(psi_m)), REAL(DIMAG(psi_m)),  perv_oa(2,per_t2p_oa(lt_m))
            end if
           else
              print*,'WARNING OA : cannot save traces in count_nmsimult_oa'
           endif
         end if
#endif
!........boucle spatiale:

           !if (if_count_space) then
           !horizontal_space_loop : do ls_m = pst%begvs_oa(iv_m),pst%begvs_oa(iv_m+1)-1
           ! i_m = pst%l2i_oa(ls_m)
           ! j_m = pst%l2j_oa(ls_m)
           !
           ! vertical_space_loop : do ls1_m = pst%begvs3d_oa(ls_m),pst%begvs3d_oa(ls_m+1)-1
           ! k_m = pst%kmin3d_oa(ls_m) + (ls1_m - pst%begvs3d_oa(ls_m))* dk_oa(iv_m)
           ! sglper_count_space(tile, ic_m,iv_m,lt_m) = sglper_count_space(tile, ic_m,iv_m,lt_m) + 1
           !
           !
           ! enddo vertical_space_loop
           !enddo  horizontal_space_loop
           !endif

       else if ( (kountv_oa(2,lt_m).le.iic_oa.or.iic_oa.eq.nt_max)      &
     &        .and.(ichoix.eq.cnb_oa(iv_m).or.ichoix.lt.0)               &
     &        .and.(ichoix.eq.des_oa(iv_m)                               &
     &        .or.des_oa(iv_m).eq.0)                                     &
     &                        ) then
            if ( tallocated_oa(lt_m)==lt_m) then
               if ( imsimult_oa .lt. 1 ) then
                write(io_unit,*)'ERROR : isimult_oa cannot be samller than 1 ',mynode
               endif
!..............Updating current count of total opened conv window

               imsimult_oa = imsimult_oa - 1
               tallocated_oa(lt_m)=-1

!..............Heisenberg box
               if_hsbg_dp2 : if ( if_heisenberg_box ) then

!              BLXD TODO Heisenberg box only dev. for Morlet wavelet
               if_mrlt_psi_dp2 : if (tpsi_oa(iv_m)==1) then

               ! lp_m is a lost info since the model time step loop is the most outer loop
               lp_m           = per_t2p_oa(lt_m)
               lp1_m          = begvp_oa(iv_m)         ! scalog. smallest per used for reso  
               lp2_m          = begvp_oa(iv_m+1) -1    ! scalog. largest per  used for conv window

               ut_m           = dble((kountv_oa(1,lt_m)+kountv_oa(2,lt_m))/2.)*dti
               t0_theo(lp_m)  = ut_m   ! BLXD eliminate ut_m

               std_dt0_theo = morlet_psi_dt0_oa(                                                           & 
                          fb_p=fb_oa                                                                         &
                         ,fc_p=fc_oa                                                                         &
                         ,scale_p=morlet_psi_p2s_oa( perv_oa(1,lp_m), fb_oa, fc_oa )                       )
               !         ,scale_p=psi_p2s_oa( tpsi_oa(iv_p), perv_oa(1,lp_m), fb_oa, fc_oa )   ! use morlet_psi_p2s_oa here
               
               ! BLXD compare temporal integral values with theo after normalizing by psi L2-norm
               ! Note that psi_norm_l2 has the size of the # of periods 
               ! but only period corresponding to pointer lp_m is here updated
               call finalize_box_oa( nzvp_oa, psi_norm_l2,    1.          & !target value psi_norm = 1. theo
                                            , t0,      ut_m               & !target value t0_psi = 0.d0 theo (no time shift)
                                            , sq_dt0,  (std_dt0_theo**2)  & !target value sq_dt0_psi  "
                                            , dt0,     std_dt0_theo       & !target value std_dt0_psi "
                                            , myprec, myprec2, lp_m ) 

        
               ! BLXD TODO include Inverse Transform here to calculate sq_dw0 => dw0
               if ( tv_sclg_oa(iv_m) ) then
                lpreso_m = lp1_m ; lpconw_m = lp2_m
               else
                lpreso_m = lp_m ; lpconw_m = lp_m
               endif

               if ( MOD( perv_oa (1,lpconw_m), dti) /= 0 ) then ! for resv_oa lp_m would works, see var_per_oa
                kount_tot_hlf = int( delta_t_oa(iv_m)*( int( perv_oa (1,lpconw_m)/ dti ) + 1 ) / resv_oa(lpreso_m) )
               else
                kount_tot_hlf = int( delta_t_oa(iv_m)*  int( perv_oa (1,lpconw_m)/ dti )       / resv_oa(lpreso_m) )
               endif
               kount_tot_low = 2 * kount_tot_hlf
               kount_tot_up  = kount_tot_low + 1

#ifdef OA_TRACES                 
               if ( ( kount_tot(lt_m) < kount_tot_low ) .or. ( kount_tot(lt_m) > kount_tot_up ) ) then
                 if(if_print_node) then
                  write(io_unit,*)'ERROR OA : total discrete concolution points kount_tot is suspicious '
                  write(io_unit,*)'           from var_per,time_oa interval should be ',kount_tot_low,kount_tot_up 
                  write(io_unit,*)'           current estimate count is               ',kount_tot(lt_m)
                 endif
                 ! stop
               endif
#endif
               kount_tot_up  = int( (kountv_oa(2,lt_m)-kountv_oa(2,lt_m))/resv_oa(lpreso_m) )
#ifdef OA_TRACES                 
               if ( kount_tot_up  /=  kount_tot(lt_m) ) then
                 if(if_print_node) then
                  write(io_unit,*)'ERROR OA : total discrete concolution points kount_tot is suspicious '
                  write(io_unit,*)'           current estimate count is               ',kount_tot(lt_m)
                  write(io_unit,*)'           OA conv window divided by resv_oa       ',kount_tot_up
                 endif
               endif
#endif

               ! Inverse Transform => Integral 
               ! - with dwi_p = 1. / (temporal conv window)  
               ! - from -N/2 to N/2  with N * dwi = 1 / dti_eff =  1 / ( dti * resv_oa )  

               call moments_of_Morlet_Power_Spectrum(                          & 
                 tpsi_p  = tpsi_oa(iv_m)                                       &
                ,m0_p    = ft_psi_norm_l2(lp_m)                                &
                ,m1_p    = w0(lp_m)                                            &
                ,m2_p    = sq_dw0(lp_m)                                        &
                ,u_p     = ut_m                                                &
                ,scale_p = morlet_psi_p2s_oa( perv_oa(1,lp_m), fb_oa, fc_oa )  &
                ,dwi_p   = 1. / (2.*delta_t_oa(iv_m)*perv_oa (1,lpconw_m))     & 
                ,fb_p    = fb_oa                                               & 
                ,fc_p    = fc_oa                                               &
                ,nw      = kount_tot_low                                     )   ! even integer
               
                dw0(lp_m)    = sqrt( sq_dw0(lp_m) )

               ! Approximate quantity dw0_apx from Heisenberg equality, leave w0 center identical
   
               dw0_apx(lp_m) =  morlet_psi_dw0_oa(                                                         & 
                          fb_p=fb_oa                                                                       &
                         ,fc_p=fc_oa                                                                       &
                         ,scale_p=morlet_psi_p2s_oa( perv_oa(1,lp_m), fb_oa, fc_oa )                       &
                         ,from_dt0_eval=dt0(lp_m)                                                          )

               w0_apx(lp_m) = w0(lp_m)
      
               endif if_mrlt_psi_dp2
               endif if_hsbg_dp2 

            endif
       endif tstp_match_convol_window

      enddo time_of_analysis_loop
 
!.....periodes integrees:
!---------------------------------------------------------------------------
      elseif (dori_oa(iv_m).eq.2) then if_discrete_per
!---------------------------------------------------------------------------

!....Morlet => Heisenberg box
!      if (tpsi_oa(iv_m)==1) then 
!        call init_temporal_box_oa( nzvp_oa, psi_norm_l2, t0, sq_dt0, dt0 ) 
!      end if

      time_of_analysis_loop2 : do lt_m  = begvt_oa(iv_m) , begvt_oa(iv_m+1) - 1

      lp2_m = begvp_oa(iv_m+1)-1

      if ( ( per_t2p_oa(lt_m) - lp2_m ) /= 0 ) then
        write(*,*) &
     & 'ERROR OA : count_nmsimult_oa integrated period must have the time pointer connected to the last period pointer ',iv_m
        if(if_print_node)write(io_unit,*) &
     & 'ERROR OA : count_nmsimult_oa integrated period must have the time pointer connected to the last period pointer ',iv_m
        stop
      end if

!.....test pour le calcul:

       tstp_match_convol_window2 : if ( iic_oa.ge.kountv_oa(1,lt_m).and.                                     &
                               iic_oa.le.kountv_oa(2,lt_m).and.                                      &
                               mod( int((iic_oa-kountv_oa(1,lt_m)),kind=8),                          &
                                    resv_oa(per_t2p_oa(lt_m)) ).eq.0                                 &
                                     .and.                                                           &
                                     (ichoix.eq.cnb_oa(iv_m).or.ichoix.lt.0)                         &
                             ) then

!........allocation si necessaire:
          if (tallocated_oa(lt_m).eq.-1) then
            imsimult_oa = imsimult_oa + 1
            nmsimult_oa = max(nmsimult_oa,imsimult_oa)
            tallocated_oa(lt_m) = lt_m
          endif

!........atome time translation (s)
!         ut_m  = real((kountv_oa(1,lt_m)+kountv_oa(2,lt_m))/2.)*dti

!........precalcul de l'indice temporel prenant en compte la translation (u):
!        Convolution window set according to the largest period eg, pointer lp_m = bevp_oa(iv_m+1)-1

         kpt_m =  iic_oa - ( int((kountv_oa(1,lt_m)+kountv_oa(2,lt_m))/2) )

         period_of_analysis_loop2 : do lp_m = begvp_oa(iv_m) , begvp_oa(iv_m+1) - 1

!........precalcul de psi: 

         ! BLXD pointer lt_m is associated to (iv_m) and the last and only the last period pointer lp_m
         ! which is the pointer for the largest period of the integartion interval
         ! see per_t2p_oa within a single period pointer loop in var_time_oa
         ! resv_oa is identical whatever lp_m is (do not depend on lp_m), perv_oa depends on lp_m
                   
         psi_m = psi_oa(                                                        &
                tpsi_oa(iv_m)                                                   &
               ,psi_p2s_oa(                                                     &
                                 tpsi_oa(iv_m)                                  &
                                ,perv_oa(1,lp_m),fb_oa,fc_oa )                  &
               ,real(kpt_m)* dti                                                &
               ,dti*resv_oa(per_t2p_oa(lt_m))                                   &
               ,fb_oa,fc_oa   ) 

!........Heisenberg box
         if_hsbg_ip1 : if ( if_heisenberg_box ) then

!           BLXD TODO Heisenberg box only dev. for Morlet wavelet
            if_mrlt_psi_ip1 : if (tpsi_oa(iv_m)==1) then

!...........Psi_oa L2-norm 
            norm_r(1) =  real( conjg(psi_m) * psi_m )

            norm_r(2) = module_morlet_psi_oa (                                  &
                    tpsi_oa(iv_m)                                               &
                   ,psi_p2s_oa(                                                 &
                                     tpsi_oa(iv_m)                              &         
                                    ,perv_oa(1,lp_m),fb_oa,fc_oa )              &                                     
                   ,real(kpt_m)* dti                                            &
                   ,dti*resv_oa(per_t2p_oa(lt_m))                               &           
                   ,fb_oa, fc_oa, from_cplx_fun  )  

            if ( abs( norm_r(1) - norm_r(2) ) > myprec2 ) then
                if(if_print_node) then
                 write(io_unit,*) 'ERROR OA : count_nmsimult_oa psi norms are too different ', norm_r(1), norm_r(2)
                 write(io_unit,*) 'ERROR OA : SHOULD STOP module_morlet_psi_oa gives bad results !!!!'
                endif
                !stop
            else 
                norm = norm_r(1)
            end if

            ! BLXD args are full size array with dimension nzvp_oa, but processing only targets lp_m
            call calc_temporal_box_oa( ndim=nzvp_oa                 &
                                      ,psi_norm=psi_norm_l2         &
                                      ,t0_psi=t0                    &
                                      ,sq_dt0_psi=sq_dt0            &
                                      ,kpt=kpt_m                    &
                                      ,dti=dti                      &
                                      ,norm=norm                    &
                                      ,ind=lp_m ) 


         endif if_mrlt_psi_ip1
         endif if_hsbg_ip1

         enddo period_of_analysis_loop2

!........boucle spatiale:
           !if (if_count_space) then
           !horizontal_space_loop2 : do ls_m = pst%begvs_oa(iv_m),pst%begvs_oa(iv_m+1)-1
           ! i_m = pst%l2i_oa(ls_m)
           ! j_m = pst%l2j_oa(ls_m)
           !
           ! vertical_space_loop2 : do ls1_m = pst%begvs3d_oa(ls_m),pst%begvs3d_oa(ls_m+1)-1
           ! k_m = pst%kmin3d_oa(ls_m) + (ls1_m-pst%begvs3d_oa(ls_m))* dk_oa(iv_m)
           !
           ! intper_count_space(tile, ic_m,iv_m,lt_m) = intper_count_space(tile, ic_m,iv_m,lt_m) + 1
           !
           !
           ! enddo vertical_space_loop2
           !enddo horizontal_space_loop2
           !endif

       else if ( (kountv_oa(2,lt_m).le.iic_oa.or.iic_oa.eq.nt_max)                   &
             .and.(ichoix.eq.cnb_oa(iv_m).or.ichoix.lt.0)                            &
             .and.(ichoix.eq.des_oa(iv_m)                                            &
             .or.des_oa(iv_m).eq.0)                                                  &
                             ) then tstp_match_convol_window2
            
              if ( tallocated_oa(lt_m)==lt_m) then

               if ( imsimult_oa .lt. 1 ) then
                write(io_unit,*)'ERROR : isimult_oa cannot be samller than 1 ',mynode
               endif
!..............Updating current count of total opened conv window
               imsimult_oa = imsimult_oa - 1
               tallocated_oa(lt_m)=-1

!..............Morlet => Heisenberg box
               if_hsbg_ip2 : if ( if_heisenberg_box ) then

!              BLXD TODO Heisenberg box only dev. for Morlet wavelet
               if_mrlt_psi_ip2 : if (tpsi_oa(iv_m)==1) then

!..............Center of convolution window set according to the largest period eg, pointer lp_m = bevp_oa(iv_m+1)-1

               ut_m           = dble((kountv_oa(1,lt_m)+kountv_oa(2,lt_m))/2.)*dti

               lp1_m          = begvp_oa(iv_m)         ! scalog. smallest per used for reso  
               lp2_m          = begvp_oa(iv_m+1) -1    ! scalog. largest per  used for conv window

               ! BLXD TODO include Inverse Transform here to calculate sq_dw0 => dw0
               !if ( tv_sclg_oa(iv_m) ) then ! NOT POSSIBLE CURRENTLY AND BOTH CONDITION IDENTICAL !
               ! lpreso_m = lp2_m ; lpconw_m = lp1_m
               !else
               ! lpreso_m = lp2_m ; lpconw_m = lp1_m
               !endif

               period_of_analysis_loop3 : do lp_m = begvp_oa(iv_m) , begvp_oa(iv_m+1) - 1

                   t0_theo(lp_m)  = ut_m

                   std_dt0_theo = morlet_psi_dt0_oa(                                                           & 
                              fb_p=fb_oa                                                                         &
                             ,fc_p=fc_oa                                                                         &
                             ,scale_p=morlet_psi_p2s_oa( perv_oa(1,lp_m), fb_oa, fc_oa )                       )
                   !         ,scale_p=psi_p2s_oa( tpsi_oa(iv_p), perv_oa(1,lp_m), fb_oa, fc_oa )   ! use morlet_psi_p2s_oa here
                   
                   ! BLXD compare temporal integral values with theo after normalizing by psi L2-norm
                   ! Note that psi_norm_l2 has the size of the # of periods 
                   ! but only period corresponding to pointer lp_m is here updated
                   call finalize_box_oa( nzvp_oa, psi_norm_l2,    1.          & !target value psi_norm = 1. theo
                                                , t0,      ut_m               & !target value t0_psi = 0.d0 theo (no time shift)
                                                , sq_dt0,  (std_dt0_theo**2)  & !target value sq_dt0_psi  "
                                                , dt0,     std_dt0_theo       & !target value std_dt0_psi "
                                                , myprec, myprec2, lp_m ) 

                   ! BLXD TODO include Inverse Transform here to calculate sq_dw0 => dw0
                   ! Requires # of points used to calculate temporal convolution
                   ! In the case of integrated period :
                   ! - the smallest period of the interval is used to calc the reso resv_oa(lp_m) in # of time steps
                   ! - the largest period is used to define the convolution interval -> pointer lp2_m 
                   ! kount = ( conv wind  of largest period) / (reso of each particular period)
                   ! dwi_p = 1 / (conv wind  of largest period)
                   ! BLXD TODO it seems that main_oa with dori_oa = 2 => resv_oa always set to per_t2p_oa(lt_m)
                   if ( MOD( perv_oa (1,lp2_m), dti) /= 0 ) then
                    !kount_tot_hlf = int( delta_t_oa(iv_m)*( int( perv_oa (1,lp2_m)/ dti ) + 1 ) / resv_oa(lp_m) )
                    kount_tot_hlf = int( delta_t_oa(iv_m)*( int( perv_oa (1,lp2_m)/ dti ) + 1 ) / resv_oa(per_t2p_oa(lt_m)) )
                   else
                    !kount_tot_hlf = int( delta_t_oa(iv_m)*  int( perv_oa (1,lp2_m)/ dti )        / resv_oa(lp_m) )
                    kount_tot_hlf = int( delta_t_oa(iv_m)*  int( perv_oa (1,lp2_m)/ dti )        / resv_oa(per_t2p_oa(lt_m)) )
                   endif
                   kount_tot_low = 2 * kount_tot_hlf
                   kount_tot_up  = kount_tot_low + 1
                     
                   if ( ( kount_tot(lt_m) < kount_tot_low ) .or. ( kount_tot(lt_m) > kount_tot_up ) ) then
                     if(if_print_node) then
                      write(io_unit,*)'ERROR OA : total discrete concolution points kount_tot is suspicious '
                      write(io_unit,*)'           from var_per,time_oa interval should be ',kount_tot_low,kount_tot_up 
                      write(io_unit,*)'           current estimate count is               ',kount_tot(lt_m)
                     endif
                     ! stop
                   endif
                   kount_tot_up  = int( (kountv_oa(2,lt_m)-kountv_oa(2,lt_m))/resv_oa(lp_m) )
                   if ( kount_tot_up  /=  kount_tot(lt_m) ) then
                     if(if_print_node) then
                      write(io_unit,*)'ERROR OA : total discrete concolution points kount_tot is suspicious '
                      write(io_unit,*)'           current estimate count is               ',kount_tot(lt_m)
                      write(io_unit,*)'           OA conv window divided by resv_oa       ',kount_tot_up
                     endif
                   endif

                   ! Inverse Transform => Integral 
                   ! - with dwi_p = 1. / (temporal conv window)  
                   ! - from -N/2 to N/2  with N * dwi = 1 / dti_eff =  1 / ( dti * resv_oa )  

                   call moments_of_Morlet_Power_Spectrum(                          & 
                     tpsi_p  = tpsi_oa(iv_m)                                       &
                    ,m0_p    = ft_psi_norm_l2(lp_m)                                &
                    ,m1_p    = w0(lp_m)                                            &
                    ,m2_p    = sq_dw0(lp_m)                                        &
                    ,u_p     = ut_m                                                &
                    ,scale_p = morlet_psi_p2s_oa( perv_oa(1,lp_m), fb_oa, fc_oa )  &
                    ,dwi_p   = 1. / (2.*delta_t_oa(iv_m)*perv_oa (1,lp2_m))        & 
                    ,fb_p    = fb_oa                                               & 
                    ,fc_p    = fc_oa                                               &
                    ,nw      = kount_tot_low                                     )   ! even integer
                   
                    dw0(lp_m) = sqrt( sq_dw0(lp_m) )

                    ! Approximate quantity dw0_apx from Heisenberg equality, leave w0 center identical
        
                    dw0_apx(lp_m) =  morlet_psi_dw0_oa(                                                        & 
                              fb_p=fb_oa                                                                       &
                             ,fc_p=fc_oa                                                                       &
                             ,scale_p=morlet_psi_p2s_oa( perv_oa(1,lp_m), fb_oa, fc_oa )                       &
                             ,from_dt0_eval=dt0(lp_m)                                                          )

                    w0_apx(lp_m) = w0(lp_m)
               
               enddo period_of_analysis_loop3

               endif if_mrlt_psi_ip2
               endif if_hsbg_ip2 

             endif
       endif tstp_match_convol_window2

      enddo time_of_analysis_loop2
      endif if_discrete_per

!.....enddo associe a iv_m et ic_m:
      enddo loop_iv
      enddo loop_ic

      enddo loop_timestepping

      if (if_print_node) then
        write(io_unit,*)'... OA precalculation of wf_oa array size (ic,iv,lt) tile/nmsimult_oa = ',tile,nmsimult_oa
      endif

      do itile=1,ntiles
        st(itile)%nzw_oa = nmsimult_oa
      enddo

      !if (if_count_space) then
      !do ic_m = 1 , nzc_oa
      !do iv_m = begc_oa(ic_m),begc_oa(ic_m+1)-1
      !  if (dori_oa(iv_m).eq.1) then
      !      do lt_m  = begvt_oa(iv_m) , begvt_oa(iv_m+1) - 1
      !      if (if_print_node) write(io_unit,*)'tile,ic,tv_oa,SGLPER OUT = ',tile,ic_m, tv_oa(iv_m),sglper_count_space(tile, ic_m,iv_m,lt_m)
      !      end do
      !  elseif (dori_oa(iv_m).eq.2) then
      !      do lt_m  = begvt_oa(iv_m) , begvt_oa(iv_m+1) - 1
      !      if (if_print_node) write(io_unit,*)'tile,ic,tv_oa,INTPER OUT = ',tile,ic_m, tv_oa(iv_m),intper_count_space(tile, ic_m,iv_m,lt_m)
      !      end do
      !  endif
      !enddo
      !enddo
      !endif

!.....Pointer to tile struct. nullified

!     if (if_count_space) then
!      pst => null()
!     endif

      endif oa_analysis_requested

      return
      end subroutine count_nmsimult_oa

!------------------------------------------------------------------------------
! PROCEDURE
!------------------------------------------------------------------------------
!
!> @note This routine is called under a Croco tile(-thread) loop and requires
!!  to prevent double allocation issues and openMP data race
!!  (all threads updating implicitly SHARED global module parameters),
!!  and to avoid multiple initialization of OA module parameters (tile-only simu),
!!  only the last tile(-threads) calls var_time_oa (with openMP critical construct)
!
! DESCRIPTION: Initialise les coefficients de reconstruction pour chaque periode et
!  calcule les incertitudes théoriques et numériques d'Heisenberg
!
!> @brief calculates the Morlet wavelet reconstruction, estimates Heisenberg boxes. 

!> @details 
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014) + Stand-alone version + optimization (01/2015)
!!  - Headers, comments
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!!  - Calculation of the Heisenberg theoritical and numerical uncertainties, 
!> @date 2021
!> @todo BLXD
!! - seperate reconstruction coeff and Heisenberg boxes calculation ?
!! - Heisenberg numerical boxes, still ongoing tests
!------------------------------------------------------------------------------

      subroutine var_rec_oa( dti )

      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      !use module_oa_stock

      implicit none

      double precision, intent(in) :: dti    !< Integration time step

      integer :: kpt_m
      real    :: signal_r 
   
! BLXD : DBLE PRECISION ? 
      complex ::                                                      &
              signal_rec                                              &
             ,temp_r                                                  &
             ,psi_m

      integer, parameter     :: ndim=2
      real,   dimension(1:ndim) :: psi_norm, norm_r
      real*8, dimension(1:ndim) :: t0_psi, sq_dt0_psi, std_dt0_psi

      real*8                    :: std_dt0_theo 
      real                      :: scale_p ! #BLXD RM when prec test ok using morlet_psi_p2s_oa      

! #BLXD TODO precision module  required
      real,   parameter :: myprec  = 10.e-5
      real*8, parameter :: myprec2 = 10.e-10

! #BLXD TODO test .false. => if identical eval  use .false. 
      logical, parameter :: from_cplx_fun = .true.

      integer                                                         &
            ip_p                                                      &
           ,iv_p                                                      &
           ,lt_p                                                      &
           ,kount_p

      integer :: io_nodoa
#ifdef OA_TRACES
      integer, dimension(1:2), parameter :: io_trc=(1000,2000), io_trc_bis=(1100,2100)
      ! BLXD integer, dimension(1:nzv_oa) :: io_traces
#endif
!

      do iv_p = 1 , nzv_oa
!---------------------------------------------------------------------------
!.....initialisation (pour tous les atomes):
!---------------------------------------------------------------------------
         do ip_p = begvp_oa(iv_p) , begvp_oa(iv_p+1)-1
            perv_oa(2,ip_p) = 1.
         enddo

!.....pour certains atomes, reconstruction possible si demandee:
         reconstruct_requested : if ( (fl_rec_oa(iv_p).eq.1 ) .or. tv_sclg_oa(iv_p) ) then

!---------------------------------------------------------------------------
!.....dirac:
!---------------------------------------------------------------------------
!TODO replace some if-endif by if-else if-endif

            type_of_oa_atom : if (tpsi_oa(iv_p).eq.0) then

               do ip_p = begvp_oa(iv_p) , begvp_oa(iv_p+1)-1
                  if (perv_oa(1,ip_p).eq.0) then
                     perv_oa(2,ip_p) = 1.
                  else
                     perv_oa(2,ip_p)=real(2*(int(perv_oa(1,ip_p)      & 
                          /2./dti-0.5) )/1)                        &
                          +3.
                  endif
               enddo
!TODO replace some if-endif by if-else if-endif
!            endif


!---------------------------------------------------------------------------
!.....fourier:
!---------------------------------------------------------------------------
!TODO replace some if-endif by if-else if-endif
!            if (tpsi_oa(iv_p).eq.2.or.tpsi_oa(iv_p).eq.3) then 

            else if (tpsi_oa(iv_p).eq.2.or.tpsi_oa(iv_p).eq.3) then type_of_oa_atom 
               lt_p  = begvt_oa(iv_p)
               if (dori_oa(iv_p).eq.2) stop 'a programmer...'
               do ip_p = begvp_oa(iv_p) , begvp_oa(iv_p+1)-1
                  do while ((per_t2p_oa(lt_p).ne.ip_p.or.kountv_oa(1,lt_p).eq.-9999.or.kountv_oa(2,lt_p).eq.-9999).and.lt_p.lt.begvt_oa(iv_p+1) )
                     lt_p=lt_p+1
                  enddo
                  perv_oa(2,ip_p)=real((kountv_oa(2,lt_p)-kountv_oa(1,lt_p))/resv_oa(ip_p)+1)
               enddo
              
!TODO replace some if-endif by if-else if-endif
!            endif


!---------------------------------------------------------------------------
!.....ondelette de morlet:
!---------------------------------------------------------------------------
!TODO replace some if-endif by if-else if-endif
!            if (tpsi_oa(iv_p).eq.1) then

             else if (tpsi_oa(iv_p).eq.1) then type_of_oa_atom 

!.....periodes discretisees:
               single_discrete_period : if (dori_oa(iv_p).eq.1) then
                         
                 
                  single_discrete_per_loop : do ip_p = begvp_oa(iv_p) , begvp_oa(iv_p+1)-1

                     temp_r = 0.

                     if (if_heisenberg_box)  call init_temporal_box_oa( ndim, psi_norm, t0_psi, sq_dt0_psi, std_dt0_psi ) 


!     boucle temporelle
                     
                     ! BLXD begvt_oa(iv_p) points toward the set of user requested temporal convolution windows
                     ! Some conv. window will not be calculated since the convolution window is crop 
                     ! either due to restart or to the end of simulation (cf, COI cone of influence)
                     ! For single discrete period analysis the relationship between
                     ! temporal convolution window lt_p and period of analysis over this convolution window
                     ! is bijective/biunivoque : a given lt_p -> a single ip_p and reciprocally  
                     ! => one given reconstruction factors perv_oa(2,ip_p) for one given lt_p temporal convolution pointer
                      
                     lt_p = begvt_oa(iv_p)

                     do while (lt_p.lt.begvt_oa(iv_p+1).and.(per_t2p_oa(lt_p).ne.ip_p.or.kountv_oa(1,lt_p).eq.-9999.or.kountv_oa(2,lt_p).eq.-9999))
                     ! #BLXD issue with array per_t2p_oa bounds overflow with intel ifort
                     !      Even if lt_p >= begvt_oa(iv_p+1) is false the following conditions are tested (per_t2p_oa(lt_p)/=ip_p) etc
                     !      => run time failure 
                     ! Adding exit condition within the loop 
                     !do while ((per_t2p_oa(lt_p).ne.ip_p.or.kountv_oa(1,lt_p).eq.-9999.or.kountv_oa(2,lt_p).eq.-9999).and.lt_p.lt.begvt_oa(iv_p+1))
                        lt_p = lt_p + 1
                        if (lt_p.eq.begvt_oa(iv_p+1)) then
                            if(if_print_node) write(io_unit,*) 'ERROR : var_rec_oa disc. per. cannot lt_p find pointer ', iv_p, perv_oa(1,ip_p)
                            exit
                        end if
                     enddo
                 
                     if ( ( if_first_per_oa(iv_p,ip_p)==0 ) .and. ( lt_p.lt.begvt_oa(iv_p+1) ) ) then
                            write(*,*) 'ERROR OA : inconsistent 1st conv. wind.'
                            if(if_print_node) then
                            write(io_unit,*) 'ERROR : two methods to get the 1st conv. wind. for var iv_p and period ip_p' 
                            write(io_unit,*) '        but inconsitent result'
                            endif
                            stop
                     end if
                     ! #BLXD 1st temporal pointer (requested by the user) for which the convolution can be performed (if any)
                     ! found_lt_pointer_single_per : if (lt_p.lt.begvt_oa(iv_p+1)) then
                      found_lt_pointer_single_per : if (if_first_per_oa(iv_p,ip_p)==1) then

                        if ( lt_p /= lper_fst_oa(iv_p,ip_p) ) then
                            write(*,*) 'ERROR OA : inconsistent conv. wind. pointer'
                            if(if_print_node) then
                            write(io_unit,*) 'ERROR : found inconsistent lt_p convolution window pointer for lper_fst_oa ',lper_fst_oa(iv_p,ip_p)
                            write(io_unit,*) 'ERROR :  see iv_p, perv_oa = ', iv_p, perv_oa(1,ip_p)
                            endif
                            stop
                        endif
          
                        loop_on_convwind_single_per : do  kount_p =  kountv_oa(1,lt_p),kountv_oa(2,lt_p),resv_oa(per_t2p_oa(lt_p))
                           kpt_m = kount_p - ( int(( kountv_oa(1,lt_p)+kountv_oa(2,lt_p))/2) )

                           
!     precalcul de psi:
                           
                           psi_m = psi_oa(                            &
                                tpsi_oa(iv_p)                         &
                                ,psi_p2s_oa(                          &
                                tpsi_oa(iv_p)                         &
                                ,perv_oa(1,ip_p),fb_oa,fc_oa )        &
                                ,real(kpt_m)* dti                     &   !- tempschoisi*dti
                                ,dti*resv_oa(per_t2p_oa(lt_p))        &
                                ,fb_oa,fc_oa  )

!     cosinus signal test
                           signal_r=cos(2*pi_oa/perv_oa(1,ip_p) *real(kpt_m)* dti)  
!                                                              *real(kount_p)* dti)  

!     reconstruction factor
                           temp_r = temp_r + conjg(psi_m) * signal_r

                           if_hsbg_dp1 : if (if_heisenberg_box) then


!     Psi L2-norm converted to real type => BLXD prec module needed !! 
                           norm_r(1)        =  real( conjg(psi_m) * psi_m )


!     Other calculation of psi L2 norm and t0_psi from morlet_psi_oa real function :
                           
                           norm_r(2) = module_morlet_psi_oa (                 &
                                tpsi_oa(iv_p)                                 &
                                ,psi_p2s_oa(                                  &
                                tpsi_oa(iv_p)                                 &
                                ,perv_oa(1,ip_p),fb_oa,fc_oa )                &
                                ,real(kpt_m)*dti                              &
                                ,dti*resv_oa(per_t2p_oa(lt_p))                &
                                ,fb_oa, fc_oa, from_cplx_fun  )  

                        call calc_temporal_box_oa( ndim=ndim              &
                                                  ,norm_ndim=norm_r       &
                                                  ,psi_norm=psi_norm      &
                                                  ,t0_psi=t0_psi          &
                                                  ,sq_dt0_psi=sq_dt0_psi  &
                                                  ,kpt=kpt_m              &
                                                  ,dti=dti ) 
                           endif if_hsbg_dp1 

#ifdef OA_TRACES
                           !if ( ( lt_p==ltrec_fst_oa(iv_p) ) ) then
                           if ( ( lt_p==lper_fst_oa(iv_p,ip_p) ) ) then
                             if ( (iv_p==1) .or. (iv_p==2) ) then 
                              io_nodoa = (io_trc(iv_p)+ip_p)*1000+mynode  !io_nodoa = io_trc(1)+mynode !1000+mynode 100X000+mynode
                              write (io_nodoa,fmt='(i4,3(1x,ES22.15E2))') kpt_m, REAL(DBLE(psi_m)), REAL(DIMAG(psi_m)), signal_r
                              if (if_heisenberg_box) then 
                              io_nodoa = (io_trc_bis(iv_p)+ip_p)*1000+mynode  !io_nodoa = io_trc(1)+mynode !1000+mynode 100X000+mynode
                              write (io_nodoa,fmt='(i4,4(1x,ES22.15E2))') kpt_m, norm_r(1), norm_r(2), psi_norm(1), psi_norm(2)
                              endif
                             end if
                            end if
#endif

                        enddo loop_on_convwind_single_per

!     signal_rec=cos(2*pi_oa/perv_oa(1,ip_p)*real(kpt_m)* dti)
!    &  +(0.,1.)*sin(2*pi_oa/perv_oa(1,ip_p)*real(kpt_m)* dti) 
!     perv_oa(2,ip_p)=abs(perv_oa(2,ip_p)/signal_rec)

                        if ( .not. tv_sclg_oa(iv_p) ) perv_oa(2,ip_p)   = real(temp_r)

                        if_hsbg_dp2 : if (if_heisenberg_box) then

                        ! #BLXD morlet_psi_ps2_oa to prepare future mod with class/method
                        scale_p = morlet_psi_p2s_oa(                                   & 
                            per_p= perv_oa(1,ip_p)                                     &
                           ,fb_p = fb_oa                                               &
                           ,fc_p = fc_oa                                               )                                      
                        
                        if ( abs( scale_p - psi_p2s_oa( tpsi_oa(iv_p), perv_oa(1,ip_p), fb_oa, fc_oa) ) > myprec2 ) then
                           write(*,*) 'ERROR OA : var_rec_oa wavelet scale calculation mismatch'
                           if(if_print_node)write(io_unit,*) 'ERROR OA : var_rec_oa wavelet scale calculation mismatch'
                           stop
                        endif
                        
                        std_dt0_theo = morlet_psi_dt0_oa(                                                           & 
                                   fb_p=fb_oa                                                                         &
                                  ,fc_p=fc_oa                                                                         &
                                  ,scale_p=morlet_psi_p2s_oa( perv_oa(1,ip_p), fb_oa, fc_oa )                       )
                        !         ,scale_p=psi_p2s_oa( tpsi_oa(iv_p), perv_oa(1,ip_p), fb_oa, fc_oa )   ! use morlet_psi_p2s_oa here
                        
                        ! BLXD compare temporal integral values with theo after normalizing by psi L2-norm
                        call finalize_box_oa( ndim, psi_norm,    1.                 & !target value psi_norm = 1. theo
                                                  , t0_psi,      0.d0               & !target value t0_psi = 0.d0 theo (no time shift)
                                                  , sq_dt0_psi,  std_dt0_theo**2    & !target value sq_dt0_psi  "
                                                  , std_dt0_psi, std_dt0_theo       & !target value std_dt0_psi "
                                                  , myprec, myprec2 ) 
                        
                        ! #BLXD TODO Check traces to choose/correct eval 1/2 
                        t0_bis (ip_p)          = t0_psi(1)          ! [T]
                        dt0_bis(ip_p)          = std_dt0_psi(1)     ! [T]
                        psi_norm_l2_bis(ip_p)  = psi_norm(1)        ! Energy pdf

                        endif if_hsbg_dp2
#ifdef OA_TRACES
                        if(if_print_node) write (io_unit,*) '       => REC FACTOR ip_p ?  ',perv_oa(1,ip_p), perv_oa(2,ip_p)
                        !if ( ( lt_p==ltrec_fst_oa(iv_p) ) ) then
                        if ( ( lt_p==lper_fst_oa(iv_p,ip_p) ) ) then
                          if ( (iv_p==1) .or. (iv_p==2) ) then 
                           io_nodoa = (io_trc(iv_p)+ip_p)*1000+mynode  !io_nodoa = io_trc(1)+mynode !1000+mynode 100X000+mynode
                           write (io_nodoa,fmt='(i4,3(1x,ES22.15E2))') kpt_m, REAL(DBLE(psi_m)), REAL(DIMAG(psi_m)), perv_oa(2,ip_p)
                           if (if_heisenberg_box) then 
                           io_nodoa = (io_trc_bis(iv_p)+ip_p)*1000+mynode  !io_nodoa = io_trc(1)+mynode !1000+mynode 100X000+mynode
                           write (io_nodoa,fmt='(i4,4(1x,ES22.15E2))') kpt_m, t0_psi(1), t0_psi(2), psi_norm(1), psi_norm(2)
                           write (io_nodoa,fmt='(i4,4(1x,ES22.15E2))') kpt_m, std_dt0_psi(1), std_dt0_psi(2), psi_norm(1), psi_norm(2)
                           endif
                          end if
                        end if
#endif
                     else found_lt_pointer_single_per

                        ! t0 (ip_p)          = t0_psi(1)          ! [T]
                        ! dt0(ip_p)          = std_dt0_psi(1)     ! [T]
                        ! psi_norm_l2(ip_p)  = psi_norm(1)        ! Energy pdf
                        ! perv_oa(2,ip_p)    = 

                        if(if_print_node) then
                        write(io_unit,*)'WARNING var_rec_oa discrete period'
                        write(io_unit,*)'WARNING : the analysis cannot work, your simulation seems to be too short'
                        write(io_unit,*)'WARNING see iv_p, lt_p, perv_oa = ', iv_p, lt_p, perv_oa(1,ip_p)
                        write(io_unit,*)'WARNING'
                        endif

                     endif found_lt_pointer_single_per

                  enddo single_discrete_per_loop

!.....periodes integrees:
               elseif (dori_oa(iv_p).eq.2) then single_discrete_period


                  integrated_period_loop : do ip_p = begvp_oa(iv_p) , begvp_oa(iv_p+1)-1


                     temp_r = 0.

                     if (if_heisenberg_box) call init_temporal_box_oa( ndim, psi_norm, t0_psi, sq_dt0_psi, std_dt0_psi ) 

!     boucle temporelle
                     
                     lt_p = begvt_oa(iv_p)

                     ! BLXD for integrated period there are several period pointer for a given temporal pointer
                     !      Integrated : per_t2p_oa(lt_m) is  associated last period only
                     ! lt_p = begvt_oa(iv_p) should correspond to ip_p = begvp_oa(iv_p+1)-1 = per_t2p_oa(lt_p)

                     do while ( lt_p.lt.begvt_oa(iv_p+1)                 .and.  &
                              ( per_t2p_oa(lt_p).ne.(begvp_oa(iv_p+1)-1) .or.   &
                                kountv_oa(1,lt_p).eq.-9999               .or.   &
                                kountv_oa(2,lt_p).eq.-9999) )
                        lt_p = lt_p + 1
                        if (lt_p.eq.begvt_oa(iv_p+1)) then
                            if(if_print_node) write(io_unit,*) 'ERROR : var_rec_oa int. per. cannot lt_p find pointer ', iv_p, perv_oa(1,ip_p)
                            exit
                        end if
                     enddo

                     if ( ( if_first_per_oa(iv_p,ip_p)==0 ) .and. ( lt_p.lt.begvt_oa(iv_p+1) ) ) then
                            write(*,*) 'ERROR OA : inconsistent 1st conv. wind.'
                            if(if_print_node) then
                            write(io_unit,*) 'ERROR : two methods to get the 1st conv. wind. for var iv_p and period ip_p' 
                            write(io_unit,*) '        but inconsitent result'
                            endif
                            stop
                     end if

                     ! #BLXD 1st temporal pointer (requested by the user) for which the convolution can be performed (if any)
                     found_lt_pointer_integrated_per : if (if_first_per_oa(iv_p,ip_p)==1) then

                        if ( lt_p /= lper_fst_oa(iv_p,ip_p) ) then
                            write(*,*) 'ERROR OA : inconsistent conv. wind. pointer'
                            if(if_print_node) then
                            write(io_unit,*) 'ERROR : found inconsistent lt_p convolution window pointer for lper_fst_oa ',lper_fst_oa(iv_p,ip_p)
                            write(io_unit,*) 'ERROR :  see iv_p, perv_oa = ', iv_p, perv_oa(1,ip_p)
                            endif
                            stop
                        endif


                     loop_on_convwind_integrated_per : do  kount_p = kountv_oa(1,lt_p) , kountv_oa(2,lt_p),resv_oa(per_t2p_oa(lt_p))
                        
                        kpt_m = kount_p - ( int((kountv_oa(1,lt_p) +kountv_oa(2,lt_p))/2) )

!     precalcul de psi:
                        
                        psi_m = psi_oa(                          &
                             tpsi_oa(iv_p)                            &
                             ,psi_p2s_oa(                        &
                             tpsi_oa(iv_p)                            &
                             ,perv_oa(1,ip_p),fb_oa,fc_oa )           &
                             ,real(kpt_m)* dti                     &  !- tempschoisi*dti
                             ,dti*resv_oa(per_t2p_oa(lt_p))        &
                             ,fb_oa,fc_oa  )

!     cosinus signal test
                        signal_r=cos(2*pi_oa/perv_oa(1,ip_p)*real(kpt_m) * dti)  
                        
!     reconstruction factor
                        temp_r = temp_r + conjg(psi_m) * signal_r

                        if_hsbg_ip1 : if (if_heisenberg_box) then

!     L2-norm 
                        norm_r(1)         =  real( conjg(psi_m) * psi_m )

!     Other calculation of psi L2 norm and t0_psi from morlet_psi_oa real function :
                           
                        norm_r(2) = module_morlet_psi_oa (                 &
                             tpsi_oa(iv_p)                                 &
                             ,psi_p2s_oa(                                  &
                             tpsi_oa(iv_p)                                 &
                             ,perv_oa(1,ip_p),fb_oa,fc_oa )                &
                             ,real(kpt_m)* dti                             &   !- tempschoisi*dti
                             ,dti*resv_oa(per_t2p_oa(lt_p))                &
                             ,fb_oa, fc_oa, from_cplx_fun  )  

                        call calc_temporal_box_oa( ndim=ndim              &
                                                  ,norm_ndim=norm_r       &
                                                  ,psi_norm=psi_norm      &
                                                  ,t0_psi=t0_psi          &
                                                  ,sq_dt0_psi=sq_dt0_psi  &
                                                  ,kpt=kpt_m              &
                                                  ,dti=dti ) 
                        endif if_hsbg_ip1

!     fin boucle temporelle
                     enddo loop_on_convwind_integrated_per

!     signal_rec=cos(2*pi_oa/perv_oa(1,ip_p)*real(! tempschoisi)* dti)
!     &  +(0.,1.)*sin(2*pi_oa/perv_oa(1,ip_p)*real(! tempschoisi)* dti) 
!     perv_oa(2,ip_p)=real(perv_oa(2,ip_p)/signal_rec)

                     if ( .not. tv_sclg_oa(iv_p) ) perv_oa(2,ip_p)   = real(temp_r)

                     if_hsbg_ip2 : if (if_heisenberg_box) then

                     ! #BLXD morlet_psi_ps2_oa to prepare future mod with class/method
                     scale_p = morlet_psi_p2s_oa(                                   & 
                         per_p= perv_oa(1,ip_p)                                     &
                        ,fb_p = fb_oa                                               &
                        ,fc_p = fc_oa                                               )                                      

                     if ( abs( scale_p - psi_p2s_oa( tpsi_oa(iv_p), perv_oa(1,ip_p), fb_oa, fc_oa) ) > myprec2 ) then
                        write(*,*) 'ERROR OA : var_rec_oa wavelet scale calculation mismatch'
                        if(if_print_node)write(io_unit,*) 'ERROR OA : var_rec_oa wavelet scale calculation mismatch'
                        stop
                     endif

                     std_dt0_theo = morlet_psi_dt0_oa(                                                           & 
                                fb_p=fb_oa                                                                         &
                               ,fc_p=fc_oa                                                                         &
                               ,scale_p=morlet_psi_p2s_oa( perv_oa(1,ip_p), fb_oa, fc_oa )                       )
                     !         ,scale_p=psi_p2s_oa( tpsi_oa(iv_p), perv_oa(1,ip_p), fb_oa, fc_oa )   ! use morlet_psi_p2s_oa here

                     ! BLXD compare temporal integral values with theo after normalizing by psi L2-norm
                     call finalize_box_oa( ndim, psi_norm,    1.                 & !target value psi_norm = 1. theo
                                               , t0_psi,      0.d0               & !target value t0_psi = 0.d0 theo (no time shift)
                                               , sq_dt0_psi,  std_dt0_theo**2    & !target value sq_dt0_psi  "
                                               , std_dt0_psi, std_dt0_theo       & !target value std_dt0_psi "
                                               , myprec, myprec2 ) 

                     ! #BLXD TODO Check traces to choose/correct eval 1/2 
                     t0_bis (ip_p)          = t0_psi(1)          ! [T]
                     dt0_bis(ip_p)          = std_dt0_psi(1)     ! [T]
                     psi_norm_l2_bis(ip_p)  = psi_norm(1)        ! Energy pdf
                
                     endif if_hsbg_ip2

                     else found_lt_pointer_integrated_per

                        if(if_print_node) then
                        write(io_unit,*)'WARNING var_rec_oa integrated period'
                        write(io_unit,*)'WARNING : the analysis cannot work, your simulation seems to be too short'
                        write(io_unit,*)'WARNING see iv_p, lt_p, perv_oa = ', iv_p, lt_p, perv_oa(1,ip_p)
                        write(io_unit,*)'WARNING'
                        endif

                     endif found_lt_pointer_integrated_per

                  enddo integrated_period_loop

!.....endif associe au test sur dori_oa
               endif single_discrete_period

            endif type_of_oa_atom

!........endif associe au test: fl_rec_w=1
         endif reconstruct_requested

!......enddo associe a la variable iv_p:
      enddo
     
      return
      end subroutine var_rec_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note in case of openMP threads, execution by the MASTER thread only
!!       see, call ipd_init_oa
!
! DESCRIPTION: preparation des variables communes + scalograms.
!
!> @brief Allocate and initializes output simple variables, scalograms 
!!   and related parameters.
!
!> @details Dans la version "tile-thread compliant" les variables var3d et var2d 
!!   ont pour dimension GLOBAL_3D_ARRAY et GLOBAL_2D_ARRAY (voir module_parameter_oa)
!!   Ces variables sont traites par le module OA  par sous-sous domaines 
!!   horizontaux compatibles avec les boucles croco en "tile-threads"
!!   (double parallelisation possible pour les sous-domaines horizontaux)
!!   c-a-d que var2d/3d_oa sont mis a jour sur les dimensions 
!!   compatibles avec les PRIVATE_SCRATCH_ARRAY de set_global_definition.h.       
!!   Actuellement, on a imin,imax=Istr-2:Iend+2, jmin,jmax=Jstr-2,Jend+2
!!   Un ajustement est possible via ip_ran et I_PRI_RANGE, J_PRI_RANGE,..etc 
!!   (point interieur seulement avec ou sans les bdy points)
!!   voir online_spectral_diags.F
!!   Le traitement fait par var_space_oa, qui definit le pointeur 
!!   de domaine spatial begvs_oa  
!!   (structure vectoriel qui empile/desempile les analyses OA)
!!   est sensible a la cle cpp MASKING. 
!!   La conversion index de pointeur begsv_oa vers indices de domaine horizontal
!!   se fait via les variables ls2i_oa/l2sj_oa qui a actuellement pour
!!   domaine d'application {Istr-2:Iend+2, Jstr-2:Jend+2}.
!!   Au final, les points (i,j) mis a jour avec les analyses OA 
!!   et restitues dans var3d_oa/var2d_oa sont sensibles a la strategie adoptee 
!!   (bdy, points interieurs, points de frontiere ouverte,...)  
!!   Si var3d_oa, var2d_oa sont uniquement des variables diagnostics
!!   on peut tout de meme souhaiter avoir les analyses aus points frontieres 
!!   (forçage)
!!   Si var3d_oa, var2d_oa deviennent des variables actives 
!!   avec par exemple des derivees spatiales par exemple, il faudra envisager
!!   les echanges MPI necessaires sur les ghost points et faire attention
!!   a la compatibilite avec la cle cpp THREE_GHOST_POINTS et les schemas d'ordre eleve
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Stand-alone version + f77 syntax twd earlier standard (2015)
!!  - tvar_oa modified alloc dim., bug with tvc_oa alloc. dim nzc_oa <- nzv_oa
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!!  - Toward a Croco-OA tile(-threads) compliant interface and OA version (2021)
!!  - Scalogram analyses.
!> @date 2021
!> @todo BLXD
!! - Croco-OA module_parameter_oa not really needed for scalogram, modify allocation ?
!------------------------------------------------------------------------------

      subroutine upd_init_oa( flag_alloc )

      use module_oa_variables
      use module_parameter_oa, only : allocate_3D_glob_array_oa_cplx  &
                                     ,allocate_2D_glob_array_oa_cplx  &
                                     ,allocate_0D_sclg_array_oa_cplx  &
                                     ,allocate_0D_sclg_array_oa_real  &
                                     ,allocate_0D_sclg_array_oa_int
      use module_oa_periode


      implicit none

      logical, intent(in) :: flag_alloc

      integer                                                         &
       ic_u                                                           &
      ,iv_u                                                           &
      ,nz_u                                                           &
      ,nzc_u


     if_flag_alloc : if ( .not. flag_alloc ) then

!$OMP MASTER
      ! BLXD : preventing tile-threads data race when updating implicitly 
      !        openMP SHARED module variables
      ! This routine should be called by the MASTER thread only + synchronization

      nzupd3d_oa = 0
      nzupd2d_oa = 0
      nzupd1d_oa = 0
      nzupd0d_oa = 0

      tvar_oa(1:maxtyp_oa,1:max_idcfg_oa,1:nmvar_oa) = 0

      set_tvar_oa_loop_config : do ic_u = 1 , nzc_oa

       nzc_u = 1
!      nzc_u = 0
!100  continue
!....Search possible configuration of the same type ie, having the same config. code tc_oa
     do while (tvar_oa(tc_oa(ic_u),nzc_u,1).ne.0)
       nzc_u = nzc_u + 1
     enddo
!     if (tvar_oa(tc_oa(ic_u),nzc_u,1).ne.0) goto 100
!....tvc_oa : if several configuration of the same type are requested tvc_oa ordinates them.
!....BLXD BUG CORRECTION tvc_oa should be allocated to nzc_oa not to nzv_oa
      tvc_oa (ic_u) = nzc_u
      var_in_cfg_loop : do iv_u = begc_oa(ic_u),begc_oa(ic_u+1)-1
      if (updv_oa(iv_u).eq.2) then ! #BLXD "variable remise a jour interactivement" what if UPDV_W=1 ?
        if (tgv3d_oa(tv_oa(iv_u)).eq.3) then
          if ( .not. tv_sclg_oa(iv_u) )  then
          nzupd3d_oa    = nzupd3d_oa + 1
          tupd_oa(iv_u) = nzupd3d_oa 
          tvar_oa(tc_oa(ic_u),nzc_u,iv_u-begc_oa(ic_u)+1) = nzupd3d_oa
          else
          ! Currently scalogram only possible for a given point i,j,k
          ! Parameters for the global scal0d_cr(1:nper_sclg_max,1:nzupd0d) array
          nzupd0d_oa    = nzupd0d_oa + 1
          tupd_oa(iv_u) = nzupd0d_oa 
          tvar_oa(tc_oa(ic_u),nzc_u,iv_u-begc_oa(ic_u)+1) = nzupd0d_oa
          !if ( scl(tile)%v2sclg(iv_u) /= -99 ) then
          !   ! Retrieving global scal0d_cr(1:nper_sclg_max,1:nzupd0d) array second dimension
          !   scl(tile)%tupd(iv_u) = nzupd0d_oa
          !endif
          endif
        else if (tgv3d_oa(tv_oa(iv_u)).eq.2) then
          if ( .not. tv_sclg_oa(iv_u) )  then
          nzupd2d_oa    = nzupd2d_oa + 1
          tupd_oa(iv_u) = nzupd2d_oa 
          tvar_oa(tc_oa(ic_u),nzc_u,iv_u-begc_oa(ic_u)+1) = nzupd2d_oa
          else
          ! Currently scalogram only possible for a given point i,j,k
          nzupd0d_oa    = nzupd0d_oa + 1
          tupd_oa(iv_u) = nzupd0d_oa 
          tvar_oa(tc_oa(ic_u),nzc_u,iv_u-begc_oa(ic_u)+1) = nzupd0d_oa
          !if ( scl(tile)%v2sclg(iv_u) /= -99 ) then
          !   ! Retrieving global scal0d_cr(1:nper_sclg_max,1:nzupd0d) array second dimension
          !   scl(tile)%tupd(iv_u) = nzupd0d_oa
          !endif
          endif
        ! #BLXD TODO eventually see if 0D/1D analysis can be useful
        !       => var1d_oa, var0d_oa
        !else if (tgv3d_oa(tv_oa(iv_u)).eq.1) then
        !  nzupd1d_oa    = nzupd1d_oa + 1
        !  tupd_oa(iv_u) = nzupd1d_oa 
        !  tvar_oa(tc_oa(ic_u),nzc_u,iv_u-begc_oa(ic_u)+1) = nzupd1d_oa
        !else if (tgv3d_oa(tv_oa(iv_u)).eq.0) then
        !  nzupd0d_oa    = nzupd0d_oa + 1
        !  tupd_oa(iv_u) = nzupd0d_oa 
        !  tvar_oa(tc_oa(ic_u),nzc_u,iv_u-begc_oa(ic_u)+1) = nzupd0d_oa
        endif
      ! #BLXD 3rd dimension of tvar_oa is the number of variable per configuration
      !       in simple case it is always 1, but for composite variables... begc_oa(ic_u+1)-begc_oa(ic_u)-1
      ! ( tc_oa(ic_u), tvc_oa(ic_u), 1 )
      endif
      enddo var_in_cfg_loop

      enddo set_tvar_oa_loop_config

      !STDALONE !.....ajout d'une variable pour le calcul de l'energie:
      !STDALONE       if (tc_oa(ic_u).eq.112.and.updv_oa(begc_oa(ic_u)).eq.2) then
      !STDALONE           nzupd3d_oa    = nzupd3d_oa + 1
      !STDALONE           tvar_oa(tc_oa(ic_u),nzc_u,begc_oa(ic_u+1)-begc_oa(ic_u)+1) = nzupd3d_oa
      !STDALONE !         nzupd2d_oa    = nzupd2d_oa + 1
      !STDALONE !         tvar_oa(tc_oa(ic_u),nzc_u,begc_oa(ic_u+1)-begc_oa(ic_u)+2) = nzupd2d_oa
      !STDALONE       endif
      !STDALONE       if (tc_oa(ic_u).eq.110.and.updv_oa(begc_oa(ic_u)).eq.2) then
      !STDALONE           nzupd3d_oa    = nzupd3d_oa + 1
      !STDALONE           tvar_oa(tc_oa(ic_u),nzc_u,begc_oa(ic_u+1)-begc_oa(ic_u)+1) = nzupd3d_oa
      !STDALONE       endif
      !STDALONE       enddo  

!.....History file:

        call history_oa(11,-1,-1,-1,-1, -1)

!$OMP END MASTER

      else if_flag_alloc

!.....allocations dynamiques:


!$OMP MASTER
      ! BLXD : preventing tile or tile-threads to allocate twice the module arrays 
      !        tile-threads => allocation by MASTER thread only
      !        tiles        => if ( .not. allocated then allocate(...     
      ! This routine should be called by the MASTER thread only + synchronization

      call allocate_3D_glob_array_oa_cplx( var3d_oa, nzupd3d_oa )
      call allocate_2D_glob_array_oa_cplx( var2d_oa, nzupd2d_oa )

      if ( scalogram_analysis ) then
          ! BLXD tmp
          ! call allocate_0D_sclg_array_oa_int( sclg2v, nzupd0d_oa )

          ! BLXD TODO Scalogram
          ! Croco interface module_parameter_oa not really needed

#ifndef MPI
            call allocate_0D_sclg_array_oa_cplx( scal0d_cr, nper_sclg_max, nzupd0d_oa )
            call allocate_0D_sclg_array_oa_int( iscal0d_cr, nzupd0d_oa )
            call allocate_0D_sclg_array_oa_int( jscal0d_cr, nzupd0d_oa )
            call allocate_0D_sclg_array_oa_int( kscal0d_cr, nzupd0d_oa )
            call allocate_0D_sclg_array_oa_real( per0d_cr, nper_sclg_max, nzupd0d_oa )
#else
        if ( if_mpi_oa ) then
            call allocate_0D_sclg_array_oa_cplx( scal0d_cr, nper_sclg_max, nzupd0d_oa )
            call allocate_0D_sclg_array_oa_int( iscal0d_cr, nzupd0d_oa )
            call allocate_0D_sclg_array_oa_int( jscal0d_cr, nzupd0d_oa )
            call allocate_0D_sclg_array_oa_int( kscal0d_cr, nzupd0d_oa )
            call allocate_0D_sclg_array_oa_real( per0d_cr, nper_sclg_max, nzupd0d_oa )
        endif
#endif
      endif
!$OMP END MASTER

      endif if_flag_alloc


      return
      end subroutine upd_init_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note Croco openMP tile(-thread), this routine must be called under 
!! a tile(-thread) loop (see suboutine init_oa and Croco source main.F).
!
! DESCRIPTION: 
!
!> @brief
!
!> @details 
!
!
! REVISION HISTORY:
!
!> @authors 
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon : 
!!  - Memory optimization when tracking isopycnal levels (2014)
!!   - the structured type array wlev, sized to nzlevel_oa intead of nzv_oa 
!!     (the total number of OA analysis requested in the simulation).
!!  - Stand-alone version (2015) : isopycne values can be picked along the rhp_t
!!    field profile (according to namelist parameters).
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!!  - Toward Croco-OA tile-(threads) compliant interface and OA version (2021) 
!!    - wlev_oa now declared as field of the derived type st in module_oa_tile
!!      fields wlev_oa%rhpt,... are domain dependant... so tile(-threads) depen
!!      dant.
!!     
!> @date 2021
!> @todo
!! - reading isopycne targets in a parameter file ?
!------------------------------------------------------------------------------
      subroutine  lev_init_oa( tile, rhp_t, rhp_t_lbound, rhp_t_ubound )

      use module_tile_oa, only : st, tile_space_str
      use module_oa_variables
      use module_oa_space
      use module_oa_level

      implicit none

      !> Lower and upper index bounds of the 3-Dimensional density array
      integer, dimension(3), intent(in) ::  &
        rhp_t_lbound, rhp_t_ubound

      ! 3-Dimensional density array passed in arguments to the OA module
      double precision, dimension(rhp_t_lbound(1):rhp_t_ubound(1),rhp_t_lbound(2):rhp_t_ubound(2),rhp_t_lbound(3):rhp_t_ubound(3)) :: &
        rhp_t  !< Density array at t-grid point

      ! tile number
      integer, intent(in) :: tile 

      type(tile_space_str), pointer :: pst => null()

      integer                                                         &             
       i_l                                                            &
      ,j_l                                                            &
      ,k_l                                                            &
      ,ls_l                                                           &
      ,ls1_l                                                          &
      ,iv_l                                                           &
      ,dim_l                                                          &
      ,izlevel_oa

   !STDALONE   double precision                                                &
   !STDALONE    val_l(1:kmax+1)
   !STDALONE
   !STDALONE   integer ifl_l
   !STDALONE
   !STDALONE   character*30 file_l

!.....BLXD isopycne_analysis true (new flag) <=> ifl_l==1  (old falg)
!          isopycne_analysis true if at least one variable as type tv_oa(iv_l)==20 <=> nzlevel_oa > 0
      
      isopycne_requested : if ( nzlevel_oa>0 ) then

!.... Pointer to tile structure
          pst => st(tile)

   !BLXD Optimization the array of structured type wlev is now sized to nzlevel_oa intead of
   !     nzv_oa (the total number of OA analysis requested in the simulation)

          izlevel_oa = 0

    !STDALONE To read the target density in a file, adpat the following lines to your directories
    !STDALONE !.....lecture du fichier contenant les caracteristiques des isopycnes
    !STDALONE !     a suivre:
    !STDALONE       open(unit=10,file='DATA/isopycnes.in')
    !STDALONE       file_l ='DATA/'//dom_c//'isopycnes.dat'
    !STDALONE       open(unit=11,file=file_l)
    !STDALONE         iv_l = 1
    !STDALONE         ls_l  =begvs_oa(iv_l)
    !STDALONE         do ls1_l=begvs3d_oa(ls_l),begvs3d_oa(ls_l+1)-1
    !STDALONE           read (10,*) val_l(ls1_l-begvs3d_oa(ls_l)+1)
    !STDALONE           write(11,*) val_l(ls1_l-begvs3d_oa(ls_l)+1)
    !STDALONE         enddo
    !STDALONE       close(10)
    !STDALONE       close(11)

!.....construction de la structure de variables necessaire
!     au suivi des isopycnes:


      variable_loop : do iv_l = 1 , nzv_oa

        if (tv_oa(iv_l).eq.20) then

            izlevel_oa            = izlevel_oa + 1
!$OMP CRITICAL(isopycne_var_level)
            lev2v_oa (izlevel_oa) = iv_l
            v2lev_oa (iv_l)       = izlevel_oa
!$OMP END CRITICAL(isopycne_var_level)
            dim_l                 = pst%begvs3d_oa(pst%begvs_oa(iv_l+1)) &
                                  - pst%begvs3d_oa(pst%begvs_oa(iv_l))
            ! BLXD_TEST_TILE_ISSUE
            call allocate_wlev_oa_ptr( tile                           &
             ,izlevel_oa                                              &
             ,iv_l                                                    &
             ,dim_l )
  !           ,begvs3d_oa(begvs_oa(iv_l+1))-begvs3d_oa(begvs_oa(iv_l)) )
            do ls_l = pst%begvs_oa(iv_l),pst%begvs_oa(iv_l+1)-1
             i_l = pst%l2i_oa(ls_l)
             j_l = pst%l2j_oa(ls_l)
             do ls1_l=pst%begvs3d_oa(ls_l),pst%begvs3d_oa(ls_l+1)-1
              k_l = pst%kmin3d_oa(ls_l) + (ls1_l - pst%begvs3d_oa(ls_l))* dk_oa(iv_l)
              pst%wlev_oa(izlevel_oa)%rhp(ls1_l - pst%begvs3d_oa(pst%begvs_oa(iv_l)) + 1 ) &
                  = rhp_t(i_l,j_l,k_l)
    !STDALONE         = val_l(ls1_l - begvs3d_oa(ls_l)+1)
    ! BLXD TODO cehck precision conversion from module_oa_type where rhp and Z are real precision, here real(8)
             enddo
            enddo
         endif

       enddo variable_loop

!.....Specific module_tile_oa pointer nullified
      pst => null()

      endif isopycne_requested

      return
      end subroutine  lev_init_oa

!----------------------------------------------------------------------
! PROCEDURE
!
!> @note this routine must be called before outside of a tile(-thread) loop.
!
! DESCRIPTION: Optimization memoire en reduisant taille du vecteur d'etat
! pour le suivi d'isopycnes
! 
!> @brief memory optimization reducing the size of the state vector
!!   for isopycnal tracking 
!! (wlev_oa dimension now nzlevel_oa instead of nzv_oa, total variable #) 
!
!> @details
!!
! REVISION HISTORY:
!
!> @authors B. Lemieux-Dudon 
!> @date 2014
!------------------------------------------------------------------------------


      subroutine count_tracked_isopycne_oa

      use module_oa_variables
      use module_oa_level, only : nzlevel_oa

      implicit none

      integer :: iv_l

      nzlevel_oa = 0
      do iv_l = 1 , nzv_oa
         if (tv_oa(iv_l).eq.20) then
            nzlevel_oa = nzlevel_oa + 1
         endif
      enddo

      end subroutine count_tracked_isopycne_oa

!----------------------------------------------------------------------
! PROCEDURE
!
!> @note Croco tile(-thread) loop : double allocation, thread data race 
!
! DESCRIPTION: allocation du vecteur d'etat pour le suivi d'isopycnes
! 
!> @brief allocation of the state vector for isopycnal tracking 
!! BLXD2014 :wlev_oa dim. now nzlevel_oa instead of nzv_oa
!
!> @details
!!
! REVISION HISTORY:
!
!> @authors B. Lemieux-Dudon 
!!  - Toward Croco-OA tile-(threads) compliant interface and OA version (2021) 
!!    - Fields wlev_oa%rhpt,... are domain dependant... so tile(-threads) depen
!!      dant.
!!    - wlev_oa now declared as field of the derived type st in module_oa_tile
!!      Croco tile(-thread) loop => Introduce st(tile)% 
!> @date 2021
!> @todo
!------------------------------------------------------------------------------
      subroutine allocate_lev_part1_oa (tile)

      use module_tile_oa, only : st, tile_space_str, allocate_wlev_oa
      use module_oa_level, only : nzlevel_oa, lev2v_oa, v2lev_oa


      implicit none
      integer, intent(in) :: tile

      type(tile_space_str), pointer :: pst => null()

      pst => st(tile)
 
      !if (.not. allocated(pst%wlev_oa ) ) allocate(pst%wlev_oa (nzlevel_oa))
      call allocate_wlev_oa(tile,nzlevel_oa)
      if (.not. allocated(lev2v_oa) ) allocate(lev2v_oa(nzlevel_oa))
      if (.not. allocated(v2lev_oa) ) allocate(v2lev_oa(nzlevel_oa))

      pst => null()


      return
      end subroutine allocate_lev_part1_oa 

!----------------------------------------------------------------------
! PROCEDURE
!
!> @note to apply under a Croco tile(-thread) loop.
!
! DESCRIPTION: allocation de la partie spatiale du vecteur d'etat pour 
!              le suivi d'isopycnes
! 
!> @brief allocation of the spatial part of the isopycnal tracking state
!  vector
!
!> @details allocates pointers rhp,z,k of the st(tile)%wlev_oa structure
!!  - Fields wf_oa%coef or wlev_oa%rhpt,... are domain dependant... 
!!    so tile(-threads) dependant.
!!  - wlev_oa now declared as field of the derived type st in module_tile_oa
!!    Croco tile(-thread) loop => array of derived type st(tile)% 
!!
! REVISION HISTORY:
!
!> @authors B. Lemieux-Dudon 
!! - st(tile)%wlev_oa(l_a)%rhp(dim_a) PRIVATE dim_a per tile
!! - rm pst from list of public module variables add tile_space_str
!! - check double allocation fault, thread data race 
!> @date 2021
!> @todo 
!------------------------------------------------------------------------------


      subroutine allocate_wlev_oa_ptr (tile,lz_a,lv_a,dim_a)

      use module_tile_oa, only : st, tile_space_str

      implicit none
      integer, intent(in) :: tile
      integer, intent(in) :: dim_a, lz_a, lv_a

      type(tile_space_str), pointer :: pst => null()

      pst => st(tile)
      if (.not. associated(pst%wlev_oa(lz_a)%z   ) ) allocate( pst%wlev_oa(lz_a)%z  (dim_a) )
      if (.not. associated(pst%wlev_oa(lz_a)%rhp ) ) allocate( pst%wlev_oa(lz_a)%rhp(dim_a) )
      if (.not. associated(pst%wlev_oa(lz_a)%k   ) ) allocate( pst%wlev_oa(lz_a)%k  (dim_a) )
      pst%wlev_oa(lz_a)%k(:) = 0
      pst%wlev_oa(lz_a)%z(:) = 0.
      pst => null()

!.....History file:
!$OMP MASTER
      call history_oa(1,lz_a,lv_a,dim_a,-1, -1)
!$OMP END MASTER


      return
      end subroutine allocate_wlev_oa_ptr


!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note to apply under a Croco tile(-thread) loop.
!
! DESCRIPTION:  Initilisation des variables de periodes et de position (i,j,k)
!               associees aux scalograms demandes dans la namelist utilisateur 
!               ceci en gerant la distribution potentielle des scalograms 
!               sur differents process MPI (et tile(-threads)) pour
!               en particulier preparer les sorties XIOS (ou non XIOS).
!
!> @brief Initialization of the variables holding the scalogram position (i,j,k)
!!        and periods, in agreement with the user namelist to ensure
!!        the XIOS or non-XIOS outputs when scalogram are distributed (or not)
!!        among different MPI process (if any) and eventually distributed among 
!!        Croco tile(-threads) in the case of a double domain segmentaion(-parallelisation).
!!     
! REVISION HISTORY:
!
!> @authors 
!! B. Lemieux-Dudon
!!  - Developments to perform online Scalogram analyses (06/2021)
!!   - from namelist options to calculation and XIOS netcdf outputs
!!   - scalograms distribution over MPI process/subdomains handled 
!!      either by XIOS or by MPI internal OA instructions (see if_mpi_oa)
!!  - Toward a Croco-OA tile(-threads) compliant OA version (2021)
!!    supporting tile(-threads) and possible dual OpenMP-MPI parallelization 
!!    of the horizontal domain.
!!  - MPI + openMP instructions.
!> @date june-july 2021
!! @todo 
!------------------------------------------------------------------------------
!
      subroutine sclg_init_coords(  tile ) 

      use module_tile_oa, only : tile_space_str, st, ntiles
      use module_oa_variables
      use module_oa_periode
      use module_oa_space     !tgv3d_oa, dk_oa
#ifdef MPI
      use module_parameter_oa, only : iminmpi, jminmpi
#endif

      implicit none

      integer, intent(in) :: tile

      type(tile_space_str), pointer :: pst => null()

      integer :: ic_m, iv_m, lp_m, lps_m, ls_m, ls1_m, i_m, j_m, k_m, iu_glob, ju_glob
#ifdef MPI
      integer :: tag, tagij, ierr, is, myrank2, myrank, req
      logical :: flag
#endif


#ifdef MPI
!$OMP MASTER
      !BLXD TEST REMOVING communicator size setting since done earlier
      !in sub. init_parameter_oa 
      !call MPI_COMM_SIZE(comm, comm_size, ierr)
      if (if_chks) then
          call MPI_COMM_RANK(comm, myrank2, ierr)
          if ( mynode /= myrank2 ) then
            print*," ERROR OA : MPI comm RANK retuns myrank2 /= mynode ",comm_size,mynode, myrank2
            stop
          endif
      endif
!$OMP END MASTER
#endif

      oa_analysis_requested : if (nzv_oa.ne.0) then

!.....Pointer association to specific module_tile_oa
      pst => st(tile)

      cfg_loop : do ic_m = 1 , nzc_oa
        var_loop : do iv_m = begc_oa(ic_m),begc_oa(ic_m+1)-1

      if_discrete_per : if (dori_oa(iv_m).eq.1) then

!.....If scalogram initialize once period axis and i,j points for outputs
      if_scal_ini : if ( tv_sclg_oa(iv_m) ) then 
        std_var_ini : if (iv_m.ne.-1.and.updv_oa(iv_m).eq.2) then

#ifdef OA_TEST_MPI_XIOS
!$OMP ATOMIC
            if (if_chks) cntvt_cr( tupd_oa(iv_m) ) = begvp_oa(iv_m+1) - begvp_oa(iv_m)
#endif

!.....openMP PRIVATE variable myrank 
      myrank = -99

!.....Local pointer to module_tile_oa array structure
       
       !loop_on_periods : do lp_u=begvt_oa(iv_u),begvt_oa(iv_u+1)-1
       loop_on_periods : do lp_m=begvp_oa(iv_m),begvp_oa(iv_m+1)-1

         if (iv_m > 1) then
            lps_m = lp_m - begvp_oa(iv_m) + 1
         else
            lps_m = lp_m
         endif

#ifdef OA_TEST_MPI_XIOS
!$OMP ATOMIC
        if (if_chks) cntvt_cr( tupd_oa(iv_m) ) = cntvt_cr( tupd_oa(iv_m) ) - 1
#endif

!.....History file:
!$OMP MASTER
        call history_oa(15,lp_m,iv_m,-1,-1, -1)
!$OMP END MASTER
        
#ifdef MPI
        if ( v2locsclg(iv_m)/=-99 ) then
            ! This scalogram is local to the current MPI process mynode
            ! All threads of the MPI process which holds the scalogram associated to iv
            ! will pass here => the thread holding this scalogram must be the one updating the variable
            ! (unless setting a barrier forcing synchronization?)
!$OMP ATOMIC
            per0d_oa(lps_m,v2locsclg(iv_m)) = perv_oa(1,lp_m)
        end if
#else
        ! MASTER thread only + synchro ?
!$OMP ATOMIC
        per0d_cr(lps_m,tupd_oa(iv_m)) = perv_oa(1,lp_m)
#endif

       enddo loop_on_periods

#ifdef OA_TEST_MPI_XIOS
       if (if_chks) then
           if ( cntvt_cr( tupd_oa(iv_m) )/=0 ) then !BLXD check OK TODO remove ?
            print*,'WARNING tile-thread ISSUE ? period count down not reaching zero ',mynode, cntvt_cr( tupd_oa(iv_m) )
            stop
           endif
       endif
#endif

       var_number_of_dimensions : if (tgv3d_oa(tv_oa(iv_m)).eq.3) then      ! three-dimensional variables


!....... BLXD Note 
!        If variable iv is a scalogram, it corresponds to a unique scalogram index isclg.
!        Scalogram isclg is part of a given MPI subdomain process and part of a unique tile(-thread) 
!        For a given iv variable the ls_m loop is executed by a single MPI process and a single thread
!        openMP data race cannot occur below this loop, this is true for any of the updated values 
!        targetting global module variables (implicitly openMP SHARED) 
         do ls_m = pst%begvs_oa(iv_m),pst%begvs_oa(iv_m+1)-1
          i_m = pst%l2i_oa(ls_m)
          j_m = pst%l2j_oa(ls_m)
#ifdef MPI
          iu_glob = i_m + iminmpi-1 
          ju_glob = j_m + jminmpi-1 
#else
          iu_glob = i_m
          ju_glob = j_m
#endif

          do ls1_m = pst%begvs3d_oa(ls_m),pst%begvs3d_oa(ls_m+1)-1

           k_m = pst%kmin3d_oa(ls_m) + (ls1_m-pst%begvs3d_oa(ls_m))* dk_oa(iv_m)
#ifdef MPI
            iscal0d_oa(v2locsclg(iv_m)) = iu_glob
            jscal0d_oa(v2locsclg(iv_m)) = ju_glob
            kscal0d_oa(v2locsclg(iv_m)) = k_m
            ! myrank PRIVATE openMP var. changed if and only if I am the thread which holds the iv-var scalogram
            myrank    = mynode
#else
            iscal0d_cr(tupd_oa(iv_m)) = iu_glob 
            jscal0d_cr(tupd_oa(iv_m)) = ju_glob 
            kscal0d_cr(tupd_oa(iv_m)) = k_m 

#endif
          enddo
         enddo

       else var_number_of_dimensions                                         ! two-dimensional variables

         do ls_m = pst%begvs_oa(iv_m),pst%begvs_oa(iv_m+1)-1
          i_m = pst%l2i_oa(ls_m)
          j_m = pst%l2j_oa(ls_m)
#ifdef MPI
          iu_glob = i_m + iminmpi-1 
          ju_glob = j_m + jminmpi-1 
#else
          iu_glob = i_m
          ju_glob = j_m
#endif

          do ls1_m = pst%begvs3d_oa(ls_m),pst%begvs3d_oa(ls_m+1)-1
#ifdef MPI
            iscal0d_oa(v2locsclg(iv_m)) = iu_glob
            jscal0d_oa(v2locsclg(iv_m)) = ju_glob
            kscal0d_oa(v2locsclg(iv_m)) = -99 

            myrank    = mynode

#else
            iscal0d_cr(tupd_oa(iv_m)) = iu_glob 
            jscal0d_cr(tupd_oa(iv_m)) = ju_glob 
            kscal0d_cr(tupd_oa(iv_m)) = -99 

#endif
          enddo
         enddo

        endif var_number_of_dimensions

#ifdef MPI

        if_oa_handles_sclg_mpi : if ( if_mpi_oa ) then

!.....MPI_Send the MPI subdomain array per0d_oa for var iv_m -> one scal.

! BLXD
!     The subroutine variable myrank is a local PRIVATE thread variable :
!     myrank = -99                     => local tile(-thread) of the current MPI process has nothing to send
!     myrank = rank of the MPI process => local tile(-thread) of thr current MPI process has scalogram data to send
!     Note that all the tile-thread must test if mynode=myrank 
!     Do not restrict MPI communication to the last thread  here


!$OMP CRITICAL (send_ini_sclg_oa)
! BLXD TODO Check if openMP CRITICAL needed ? 
!           An openMP BARRIER at the end of condition mynode==myrank seems sufficiant 
! Either Croco sequential tiles or Croco tile-threads (with openMP CRITICAL REGION) 
! will process the following instructions one after one 
! CRITICAL region 
!   => there is an implicit BARRIER with synchronization of global module var. (implicitly SHARED)

           if_node_correct : if ( mynode == myrank ) then
            if_scalogram_not_on_root : if ( myrank /= root ) then

            ! MPI subprocess scalogram iv_m in buffer to send to root process
            ! with tag set to the global scalogram index 
            ! MPI tags for per and ij-points
            tag   = tupd_oa(iv_m)
            tagij = tupd_oa(iv_m)

            nptr :if (mpi_send_buff_ptr) then
                buffr_s(1:nper_sclg(iv_m)) => per0d_oa(1:nper_sclg(iv_m),v2locsclg(iv_m))
                call MPI_Send (buffr_s, nper_sclg(iv_m), MPI_DOUBLE_PRECISION, root, tag, comm, ierr)
            else nptr
                !buff2r_s(1:nper_sclg_max) = dble(per0d_oa(1:nper_sclg_max,v2locsclg(iv_m)))
                buff2r_s(1:nper_sclg_max) = per0d_oa(1:nper_sclg_max,v2locsclg(iv_m))
                call MPI_Send (buff2r_s, nper_sclg_max, MPI_DOUBLE_PRECISION, root, tag, comm, ierr)
            endif nptr
            if (ierr/=0) then
                if(if_print_node)write(io_unit,*)'ERROR OA : MPI_(I)Send scalogram period for node/ierr/var ',mynode, ierr, iv_m
                print*, 'ERROR OA : MPI_(I)Send scalogram period for node/ierr/var ',mynode, ierr, iv_m
                stop 
            else 
            endif
            if (if_record_sclg_ijpoints) then
                buffij_s(1)=iscal0d_oa(v2locsclg(iv_m))
                buffij_s(2)=jscal0d_oa(v2locsclg(iv_m))
                buffij_s(3)=jscal0d_oa(v2locsclg(iv_m))
                call MPI_Send (buffij_s,  3, MPI_INTEGER, root, tagij, comm, ierr)
                if (ierr/=0) then
                    if(if_print_node)write(io_unit,*)'ERROR OA : MPI_(I)Send scalogram period for node/ierr/var ',mynode, ierr, iv_m
                    print*, 'ERROR OA : MPI_(I)Send scalogram i,jscal0d_oa for node/ierr/var ',mynode, ierr, iv_m
                    stop 
                else 
                endif
            endif
            else if_scalogram_not_on_root
               ! Root rank process holds the saclogram then no MPI exchange !!!
                per0d_cr(1:nper_sclg(iv_m),tupd_oa(iv_m))= per0d_oa(1:nper_sclg(iv_m),v2locsclg(iv_m))
                iscal0d_cr(tupd_oa(iv_m))                = iscal0d_oa(v2locsclg(iv_m))
                jscal0d_cr(tupd_oa(iv_m))                = jscal0d_oa(v2locsclg(iv_m))
                kscal0d_cr(tupd_oa(iv_m))                = kscal0d_oa(v2locsclg(iv_m))
            endif if_scalogram_not_on_root

           endif if_node_correct

!$OMP END CRITICAL (send_ini_sclg_oa)


!..... If scalogram : iv_m -> one scal. hold by one MPI process (and eventually one tile/tile-thread)
!..... The MPI root process knows that iv is a scalogram and therefore waits for the scalogram sent by another process
!      => require an MPI Receive for arrays per0d_oa and i,j grid points 
!         except if the root process holds the scalogram itself

! BLXD TODO tile_ount_oa is necessary for Croco tile-only simulation 
!      (as opposed to Croco tile-threads with openMP)
!      as a result openMP CRITICAL region must replace MASTER region

!$OMP CRITICAL (recv_ini_sclg_oa)
!      if tile-thread (CRITICAL REGION) or if sequential tiles
!      => all the tile will do the tile_count instruction one after the other

          tile_count_oa=tile_count_oa+1

!      The last tile will do the MPI receiving job

          last_tile : if ( tile_count_oa .eq. ntiles) then 

          ! Receving the scalogram period and ijk-points sent by one MPI process node
          ! if tile-thread (CRITICAL REGION) or if sequential tiles
          ! => only the last tile performs the MPI comm
           if_node_root : if ( mynode == root ) then
            sending_rank_is_not_root : if ( myrank == -99 ) then


            ! MPI tags for per and ij-points
            tag   = tupd_oa(iv_m)
            tagij = tupd_oa(iv_m)

            ! BLXD TO TEST reduced buffer size !!!! 
            nptrr :if (mpi_send_buff_ptr) then
                call MPI_Recv( per0d_cr(1:nper_sclg(iv_m),tupd_oa(iv_m)), nper_sclg(iv_m), &
                        MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, comm, mpi_status, ierr )
            else nptrr
                call MPI_Recv( per0d_cr(1:nper_sclg_max,tupd_oa(iv_m)), nper_sclg_max,     &
                    MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, comm, mpi_status, ierr )
            endif nptrr
            if ( ierr /= 0 ) then
                write(io_unit,*)'ERROR OA : MPI_Recv scalogram per0d_cr ierr ',mynode, ierr, iv_m
                print*,'ERROR OA : MPI_Recv scalogram per0d_cr ierr ',mynode, ierr, iv_m
                stop
            endif


            if (if_record_sclg_ijpoints) then
                call MPI_Recv( buffij_r, 3, MPI_INTEGER, MPI_ANY_SOURCE, tag, comm, mpi_status, ierr )
                if ( ierr /= 0 ) then
                    write(io_unit,*)'ERROR OA : MPI_Recv ijscal0d_oa mynode/ierr/iv_m ',mynode, ierr, iv_m
                    print*,'ERROR OA : MPI_Recv ijscal0d_oa mynode/ierr/iv_m ',mynode, ierr, iv_m
                    stop
                endif
                
                ! POSSIBLEMENT TOUT bufferiser avec do ic do if if scal on MPI process
                ! using index_s_cr
                iscal0d_cr(tupd_oa(iv_m))=buffij_r(1)
                jscal0d_cr(tupd_oa(iv_m))=buffij_r(2)
                kscal0d_cr(tupd_oa(iv_m))=buffij_r(3)


            endif



            endif sending_rank_is_not_root

           endif if_node_root
           
           call MPI_Barrier(comm,ierr)
           if (ierr/=0) then
               print*,'ERROR OA : MPI_Barrier bef. per/ijkscal Bcast mynode/ierr ',mynode,ierr
               stop
           endif
           tag=tupd_oa(iv_m)
           call MPI_Bcast(per0d_cr(1:nper_sclg(iv_m),tupd_oa(iv_m)), nper_sclg(iv_m), &
                            MPI_DOUBLE_PRECISION, root, comm, ierr)
           if (if_record_sclg_ijpoints) then
               call MPI_Bcast(iscal0d_cr(tupd_oa(iv_m)), 1, MPI_INTEGER, root, comm, ierr)
               call MPI_Bcast(jscal0d_cr(tupd_oa(iv_m)), 1, MPI_INTEGER, root, comm, ierr)
               call MPI_Bcast(kscal0d_cr(tupd_oa(iv_m)), 1, MPI_INTEGER, root, comm, ierr)
           endif


! The last tile Sets back to zero the module tile count variable

           tile_count_oa=0

          endif last_tile

!$OMP END CRITICAL (recv_ini_sclg_oa)

        endif if_oa_handles_sclg_mpi
#endif /* MPI */

        !endif if_first_time_a_scalogram                        
        endif std_var_ini
      endif if_scal_ini

      endif if_discrete_per
      enddo var_loop
      enddo cfg_loop

      pst => null()

      endif oa_analysis_requested 

      return
      end subroutine sclg_init_coords

!----------------------------------------------------------------------
! PROCEDURE
!
!> @note this routine is called under a Croco tile(-thread) loop 
!!   (see croco sources main.F and step.F).
!
! DESCRIPTION: 
! subroutine appelee par le programme principal
! pour chaque variable ("configuration") main_oa teste
! si l'on se trouve dans une periode de "calcul" et si tel
! est le cas il gere le calcul du coef...
! outputs: wf_oa: variable contenant tous les coefs...
!
!> @brief performs the time convolution between numerical fields and 
!! the requested "atoms" (i.e., wavelets, Fourier, Dirac). 
!
!> @details For each OA variariable-configuration (fields, scale), at each time step 
!! cumulates the convolution product between the evolving numerical field 
!! and the requested "atom" (i.e., wavelets, Fourier, Dirac).
!! The resulting coefficient is stored in the wf_oa structure.
!!
!! In the case of Fourier or wavelets, two types of analysis are available:
!! - dori_oa set to 1 : no scale integration, analysis at a discrete given scale 
!! - dori_oa set to 2 : scale integration
!!
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014) + Stand-alone version + optimization (01/2015)
!!  - Headers, comments, cleaning, f77 to earlier standards (intent attr.,..)
!!  - Croco-OA interface 1st version (2020)
!!    Interface = mix between the Stand-alone and Original OnlineA versions 
!!  - Toward a Croco-OA tile(-threads) compliant OA version (2021)
!!    supporting tile(-threads) and possible dual OpenMP-MPI parallelization 
!!    of the horizontal domain.
!!  - MPI + openMP instructions, use of st(tile)% derived type for
!!    tile(-thread) dependant variables (state vector wf_oa%coef,...)
!!  - Developments to perform online Scalogram analyses (06/2021)
!!   - from namelist options to calculation and XIOS netcdf outputs
!!   - scalograms distribution over MPI process/subdomains handled 
!!     either by XIOS or by MPI internal OA instructions (see if_mpi_oa)
!
!> @date 2021
!> @todo BLXD Modify Croco-OnlineA module interface.
!  - Croco-OnlineA module interface, 1st version, Spring 2020
!!   Interface = mix between the Stand-alone and Original OnlineA versions 
!!   => Croco arrays are passed as arguments to the OnlineA routines with 
!!   reduced dimension range. Check memory duplication.
!! @todo BLXD 
!!  - specific module to set consitent precisions 
!!    with real(KIND=wp)... etc 
!!    eg, var_oa function defined as real => real(8)
!!    eg, complex =>
!!    double precision is an obsolete type TO REPLACE
!!  - introduce tiny instead of 1.e-30 ?
!!  - still required to test with a Croco tile(-thread) working version 
!----------------------------------------------------------------------

      subroutine main_oa(  tile                       &
                          ,ichoix                     &
                          ,ivar_m                     &
                          ,iic_oa                     &
                          ,imin, imax                 &
                          ,jmin, jmax                 &
                          ,kmin, kmax                 &
                          ,rhp_t                      & !(-1:imax+2,-1:jmax+2,0:kmax+1) 
                          ,depth_t                    ) 

      use module_tile_oa, only : tile_space_str, st, deallocate_tile_varspace_oa  &
                                ,deallocate_tile_sclg_oa, deallocate_tile_test_oa &
                                ,ntiles
      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
!      use module_oa_stock

      implicit none

      integer, intent(in) :: ichoix

      !> Current model integration iteration
      integer, intent(in) :: iic_oa 

      !> To handle specific calls to main_oa with variables updated at different position in the calling code
      integer, intent(in), optional :: ivar_m
     
      !> Grid index range for the analysis of fields passed in argument to the OA module 
      !! Grid index are Roms-Croco tile(-thread) indices.
      !! If no tile(-thread) grid index are MPI subdomain index
      !! If no tile(-thread) and no MPI subdomain decomposition => full domain grid range
      integer, intent(in) ::       &
             imin, imax            & 
            ,jmin, jmax            &
            ,kmin, kmax

      !> tile parameter from croco
      integer, intent(in) :: tile

      !> 3-Dimensional density array at t-grid point
      double precision, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in), target :: rhp_t
                                                                          
      !> 3-Dimensional depth array at t-grid point
      double precision, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in), target :: depth_t

      !> Lower and upper index bounds of the 3-Dimensional density array (see if pst%imin,...)
      integer, dimension(3)  ::  rhp_t_lbound, rhp_t_ubound

      type(tile_space_str), pointer :: pst => null()

      !logical, dimension(1:nzv_oa) :: isclg_send_recv
      !integer, dimension(1:nzv_oa) :: rank_sending  

      !complex :: psi_oa
      !real    :: psi_p2s_oa

      real    :: var_oa
      
      !TODO declare flag here if only needed here
      logical :: ifl_test_composite

      integer                                                         &
         iv_m                                                         &
        ,ic_m                                                         &
        ,lt_m                                                         &
        ,ls_m                                                         &
        ,ls1_m

      integer                                                         &
         i_m                                                          &
        ,j_m                                                          &
        ,k_m                                                          &
        ,kpt_m                                                        &
        ,lp_m                                                         &
        ,la_m
    
! BLXD TODO : what precision ? 
      complex                                                         &
         psi_m

! BLXD var_oa function real(4) ?
      real :: tmpvar

! BLXD introduce tiny instead of 1.e-30
! ....Comment : myrank  is a PRIVATE openMP var.
      integer :: io_nodoa, myrank


#ifdef MPI
      integer :: tag, tagij, ierr
#endif


!---------------------------------------------------------------------------
!.....periodes discretisees:
!---------------------------------------------------------------------------

      oa_analysis_requested : if (nzv_oa.ne.0) then

!.....Pointer association to specific module_tile_oa
      pst => st(tile)

!.....flag variable composite 
      ifl_test_composite = .true.      

!.....test variable
       if (test_analysis) then
        call test_oa( tile=tile, iic_oa=iic_oa )
         !,imin=imin, imax=imax     &
         !,jmin=jmin, jmax=jmax     &
         !,kmin=kmin, kmax=kmax )
       endif

      !STDALONE if (flag_nrj_oa.eq.1) call nrj_upd_oa(1)

      if ( isopycne_analysis ) then
        rhp_t_lbound = (/imin,jmin,kmin/)
        rhp_t_ubound = (/imax,jmax,kmax/)
        call lev_upd_oa(     tile                                          &
                             ,rhp_t                                        &
                             ,depth_t                                      & !(-1:imax+2,-1:jmax+2,0:kmax+1) 
                             ,rhp_t_lbound                                 &
                             ,rhp_t_ubound  )
      endif


      do ic_m = 1 , nzc_oa
      do iv_m = begc_oa(ic_m),begc_oa(ic_m+1)-1
         !STDALONE to call main_oa at specific location for composite varaible confguration
         if (present(ivar_m)) then
            ifl_test_composite = tv_oa(iv_m).eq.ivar_m
         endif

!........ce test est utilise pour les appels cibles lors des calculs energetiques:
         if_energy_calc : if (ichoix.ge.0.or.ifl_test_composite) then
!---------------------------------------------------------------------------
!.....periodes discretisees:
!---------------------------------------------------------------------------

     if_discrete_per : if (dori_oa(iv_m).eq.1) then

      time_of_analysis_loop : do lt_m  = begvt_oa(iv_m) , begvt_oa(iv_m+1) - 1


!.....test pour le calcul:
       tstp_match_convol_window : if ( iic_oa.ge.kountv_oa(1,lt_m) .and. iic_oa.le.kountv_oa(2,lt_m)                       &
        &  .and.               mod( int( (iic_oa-kountv_oa(1,lt_m)), kind=8 ), resv_oa(per_t2p_oa(lt_m)) ).eq.0            &
        &  .and.               ( ichoix.eq.cnb_oa(iv_m) .or. ichoix.lt.0 )                                                 &
                             ) then


! #BLXD ( ichoix.eq.cnb_oa(iv_m) .or. ichoix.lt.0 ) = calling main_oa for particular variable

!........allocation if necessary (if variable isn't allocated yet, it is done.., at user-specified kount):


          test_if_allocation_needed_dp : if (tallocated_oa(lt_m).eq.-1) then

            if ( .not. associated(pst) ) then
                print*, 'ERROR MAIN_OA : pst not associated alloc_win_oa',ic_m,iv_m,lt_m,mynode
                stop
            endif


            call allocate_win_oa ( tile                                                        &
                    ,lt_m, ic_m, iv_m                                                          &
                    ,pst%begvs3d_oa(pst%begvs_oa(iv_m+1)) - pst%begvs3d_oa(pst%begvs_oa(iv_m)) &
                    ,iic_oa, la_m )
            

           endif test_if_allocation_needed_dp          


!........precalcul de l'indice temporel prenant en compte la translation (u):

         kpt_m =  iic_oa - ( int((kountv_oa(1,lt_m)+kountv_oa(2,lt_m))/2) )

!........precalcul de psi:

          psi_m = psi_oa(                                                     &
                tpsi_oa(iv_m)                                                 &
               ,psi_p2s_oa(                                                   &
                                 tpsi_oa(iv_m)                                &         
                                ,perv_oa(1,per_t2p_oa(lt_m)),fb_oa,fc_oa )    &                                     
               ,real(kpt_m)* dti                                              &
               ,dti*resv_oa(per_t2p_oa(lt_m))                                 &           
               ,fb_oa,fc_oa  ) 
  
!........boucle spatiale:
!        BLXD No openMP data race in the below loop, to a given var. index is
!        associated a set of (i,j,[k]) analysis grid-points
!        Each tile-thread covers disjoint part of the set of analysed points

           horizontal_space_loop : do ls_m = pst%begvs_oa(iv_m),pst%begvs_oa(iv_m+1)-1
            i_m = pst%l2i_oa(ls_m)
            j_m = pst%l2j_oa(ls_m)

            vertical_space_loop : do ls1_m = pst%begvs3d_oa(ls_m),pst%begvs3d_oa(ls_m+1)-1
            k_m = pst%kmin3d_oa(ls_m) + (ls1_m - pst%begvs3d_oa(ls_m))* dk_oa(iv_m)
            tmpvar =  var_oa( tile                                               &
                                   ,tv_oa(iv_m)                                  &
                                   ,ichoix                                       &
                                   ,i_m                                          &
                                   ,j_m                                          &
                                   ,k_m                                          &
                                   ,iv_m                                         &
                                   ,ls1_m-pst%begvs3d_oa(pst%begvs_oa(iv_m))+1       ) 

             !BLXD_TILE_ISSUE one dim_a spatial dimension per tile(-thread)
             !wf_oa transformed into pst%wf_oa(tallocated_oa(lt_m))%coef ( dim_a ) 

             pst%wf_oa(tallocated_oa(lt_m))%coef (                               &
                      ls1_m-pst%begvs3d_oa(pst%begvs_oa(iv_m))+1 )   =           &
             pst%wf_oa(tallocated_oa(lt_m))%coef (                               &
                      ls1_m-pst%begvs3d_oa(pst%begvs_oa(iv_m))+1 )               &
                     + conjg(psi_m) * tmpvar                                     &
                     / max(perv_oa(2,per_t2p_oa(lt_m)),1.e-30)

            enddo vertical_space_loop
           enddo  horizontal_space_loop


       endif tstp_match_convol_window
      enddo time_of_analysis_loop
 
!.....periodes integrees:
!---------------------------------------------------------------------------
      elseif (dori_oa(iv_m).eq.2) then if_discrete_per
!---------------------------------------------------------------------------

      time_of_analysis_loop2 : do lt_m  = begvt_oa(iv_m) , begvt_oa(iv_m+1) - 1

!.....test pour le calcul:

       tstp_match_convol_window2 : if ( iic_oa.ge.kountv_oa(1,lt_m) .and. iic_oa.le.kountv_oa(2,lt_m)                      &
        &  .and.               mod( int( (iic_oa-kountv_oa(1,lt_m)), kind=8 ), resv_oa(per_t2p_oa(lt_m)) ).eq.0            &
        &  .and.               ( ichoix.eq.cnb_oa(iv_m) .or. ichoix.lt.0 )                                                 &
                             ) then


!........allocation si necessaire:

         if (tallocated_oa(lt_m).eq.-1) then
            call allocate_win_oa ( tile                                                        &
                    ,lt_m, ic_m, iv_m                                                          &
                    ,pst%begvs3d_oa(pst%begvs_oa(iv_m+1)) - pst%begvs3d_oa(pst%begvs_oa(iv_m)) &
                    ,iic_oa, la_m )
         endif

!........precalcul de l'indice temporel prenant en compte la translation (u):

         kpt_m =  iic_oa  - ( int((kountv_oa(1,lt_m)+kountv_oa(2,lt_m))/2) )
     
!........boucle sur toutes les periodes:

         period_of_analysis_loop2 : do lp_m = begvp_oa(iv_m) , begvp_oa(iv_m+1) - 1

!........precalcul de psi: 

         ! BLXD pointer lt_m is associated to (iv_m) and the last and only the last period pointer lp_m
         ! which is the pointer for the largest period of the integartion interval
         ! see per_t2p_oa within a single period pointer loop in var_time_oa
         ! resv_oa is identical whatever lp_m is (do not depend on lp_m), perv_oa depends on lp_m
          
         psi_m = psi_oa(                                                        &
                tpsi_oa(iv_m)                                                   &
               ,psi_p2s_oa(                                                     &
                                 tpsi_oa(iv_m)                                  &
                                ,perv_oa(1,lp_m),fb_oa,fc_oa )                  &
               ,real(kpt_m)* dti                                             &
               ,dti*resv_oa(per_t2p_oa(lt_m))                                &
               ,fb_oa,fc_oa   ) 

!........boucle spatiale:
           horizontal_space_loop2 : do ls_m = pst%begvs_oa(iv_m),pst%begvs_oa(iv_m+1)-1
            i_m = pst%l2i_oa(ls_m)
            j_m = pst%l2j_oa(ls_m)

            vertical_space_loop2 : do ls1_m = pst%begvs3d_oa(ls_m),pst%begvs3d_oa(ls_m+1)-1
            k_m = pst%kmin3d_oa(ls_m) + (ls1_m-pst%begvs3d_oa(ls_m))* dk_oa(iv_m)

            tmpvar = var_oa( tile                                               &
                                   ,tv_oa(iv_m)                                 &
                                   ,ichoix                                      &
                                   ,i_m                                         &
                                   ,j_m                                         &
                                   ,k_m                                         &
                                   ,iv_m                                        &
                                   ,ls1_m-pst%begvs3d_oa(pst%begvs_oa(iv_m))+1 )
             !BLXD_TILE_ISSUE
             !wf_oa(tallocated_oa(lt_m))%coef (                                  &
             !         ls1_m-pst%begvs3d_oa(pst%begvs_oa(iv_m))+1 ) =            &
             !wf_oa(tallocated_oa(lt_m))%coef (                                  &
             !         ls1_m-pst%begvs3d_oa(pst%begvs_oa(iv_m))+1 )              &
             !        + conjg(psi_m) * tmpvar                                    &
             !        / max(perv_oa(2,lp_m),1.e-30)

             pst%wf_oa(tallocated_oa(lt_m))%coef (                              &
                      ls1_m-pst%begvs3d_oa(pst%begvs_oa(iv_m))+1 ) =            &
             pst%wf_oa(tallocated_oa(lt_m))%coef (                              &
                      ls1_m-pst%begvs3d_oa(pst%begvs_oa(iv_m))+1 )              &
                     + conjg(psi_m) * tmpvar                                    &
                     / max(perv_oa(2,lp_m),1.e-30)

            enddo vertical_space_loop2
           enddo horizontal_space_loop2
         enddo period_of_analysis_loop2
       endif tstp_match_convol_window2

      enddo time_of_analysis_loop2
      endif if_discrete_per

!.....enddo associe a iv_m et ic_m:
      endif if_energy_calc
      enddo
      enddo

!.....eventuelle desalocation:

      do la_m = 1,nmsimult_oa
         ! BLXD_TILE_ISSUE 
         ! instructions done by all MPI processes and tile-(threads) since they all 
         ! know global module variables nmsimult_oa, kountv_oa, etc...
         ! and they all see the same temporal structure, and they all understand 
         ! that there is a convolution window ending


         if ( associated(pst%wf_oa(la_m)%coef) ) then

         if_ending_convwind : if ((kountv_oa(2,pst%wf_oa(la_m)%t_indice).le.iic_oa.or.iic_oa.eq.nt_max)      &
           .and.(ichoix.eq.cnb_oa(pst%wf_oa(la_m)%variable).or.ichoix.lt.0)          &
           .and.(ichoix.eq.des_oa(pst%wf_oa(la_m)%variable)                          &
                 .or.des_oa(pst%wf_oa(la_m)%variable).eq.0)                          &
                 ) then                

            lt_m = pst%wf_oa(la_m)%t_indice
            iv_m = pst%wf_oa(la_m)%variable
            ic_m = pst%wf_oa(la_m)%config



            call subsave_oa ( tile, la_m, myrank )                           
            call deallocate_win_oa(tile, lt_m)


#ifdef MPI

      if_oa_handles_sclg_mpi : if ( if_mpi_oa ) then

!.....If scalogram : MPI Send process subdomain array scal0d_oa for var iv_u -> one scal.
!     If range of period completed and if all tiles completed + If test in var_upd_oa ok
        std_var : if (iv_m.ne.-1.and.updv_oa(iv_m).eq.2) then
        if_sclg : if ( tv_sclg_oa(iv_m) )  then

            ! Check if this analysis la_m <-> ic_m, iv_m, lt_m is the last period of the scalogram
            if_sclg_complete : if ( per_t2p_oa(lt_m) == begvp_oa( iv_m + 1 ) - 1 )  then

!$OMP CRITICAL (recv_sclg_oa)
!      if tile-thread (CRITICAL REGION) or if sequential tiles
!      => all the tile will do the tile_count instruction one after the other

          tile_count_oa=tile_count_oa+1


          ! if tile-thread (CRITICAL REGION) or if sequential tiles
          ! => only the last tile performs the MPI comm
          last_tile : if ( tile_count_oa .eq. ntiles) then 

          if_node_correct2 : if ( mynode == root ) then


            sending_rank_is_not_root : if ( myrank == -99 ) then

             tag=tupd_oa(iv_m)

             sclg_br : if ( mpi_blocking_recv ) then


                 call MPI_Recv( scal0d_cr(1:nper_sclg(iv_m),tupd_oa(iv_m)), nper_sclg(iv_m), &
                                MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, tag, comm, mpi_status, ierr )

                 if ( ierr /= 0 ) then
                     write(io_unit,*)'ERROR OA : MPI_Recv scalogram per0d_cr ierr ',mynode, ierr, iv_m
                     stop
                 endif

               

             else sclg_br
                write(io_unit,*)'ERROR OA : After subsave_oa NON BLOCKING RECV TO DEV ',mynode, iv_m
                stop

             endif sclg_br

            else sending_rank_is_not_root

#ifdef OA_TRACES
                 if(if_print_node) write (io_unit,*) '*** MAIN_OA : after subsave_oa ROOT node with scalogram no MPI_Recv ',mynode,iv_m
#endif

            endif sending_rank_is_not_root


           endif if_node_correct2

           tag=tupd_oa(iv_m)
           call MPI_Barrier(comm,ierr)
           if (ierr/=0) then
               print*,'ERROR OA : MPI_Barrier bef. scal Bcast mynode ',mynode
               stop
           endif
           call MPI_Bcast(scal0d_cr(1:nper_sclg(iv_m),tupd_oa(iv_m)), nper_sclg(iv_m), &
                            MPI_DOUBLE_COMPLEX, root, comm, ierr)

! The last tile Sets back to zero the module tile count variable

           tile_count_oa=0

          endif last_tile

!$OMP END CRITICAL (recv_sclg_oa)

           endif if_sclg_complete 


        endif if_sclg 
        endif std_var

        endif if_oa_handles_sclg_mpi

#endif  /* MPI */

         endif if_ending_convwind
         endif
      enddo


!.....Specific module_tile_oa pointer nullyfied
      pst => null()

!.....If end of simulation Deallocate deallocate_tile_varspace_oa
!     consistently with the deallocate_tile_space_oa removing st
!     in online_spectral_diags after the $OMP BARRIER
      if (iic_oa.eq.(nt_max+1)) then
        call deallocate_tile_varspace_oa(tile)
        if (scalogram_analysis) call deallocate_tile_sclg_oa(tile)
        if (test_analysis) call deallocate_tile_test_oa(tile)
      endif

      endif oa_analysis_requested

      return
      
      end subroutine main_oa

 
!----------------------------------------------------------------------
! PROCEDURE
!
!> @note This routine must be called under a Croco tile(-thread)
!!       loop (see main_oa subroutine and Croco sources main.F/step.F).
!
! DESCRIPTION: sauvegarde des coefficients qui viennent d'etre calcules.
!
!> @brief 
!
!> @details
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - 2021 Scalograms + MPI + XIOS outputs
!!  - 2021 Twd Croco tile(-threads) compliant OA version
!!  - 2020 Croco interface 1st version. 
!!  - 2014-2015 optimization, f77 to earlier stds, Stand-alone version, Headers
!> @date 2021
!> @todo
!----------------------------------------------------------------------

      subroutine subsave_oa( tile, la_s, myrank )
 
      use module_oa_variables
      !use module_oa_stock
      use module_tile_oa, only : st, tile_space_str


      implicit none

      integer, intent(in)  :: tile, la_s    
      integer, intent(out) :: myrank
 
      type(tile_space_str), pointer :: pst => null()

      integer ::                                                      &
            iv_s                                                      &
           ,lt_s                                                      &
           ,ic_s                                                       

      pst => st(tile)
 
      ic_s = pst%wf_oa(la_s)%config
      iv_s = pst%wf_oa(la_s)%variable
      lt_s = pst%wf_oa(la_s)%t_indice

      if (updv_oa(iv_s).ne.0) then
       call var_upd_oa ( tile, ic_s, iv_s, lt_s, la_s, myrank )
      endif

      pst => null()

      return
      end subroutine subsave_oa

!------------------------------------------------------------------------------
! PROCEDURE
!------------------------------------------------------------------------------
!
!> @note This routine must called under a Croco tile(-thread)
!!       loop (see routine subsave_oa,main_oa and croco sources main.F/step.F).
!
! DESCRIPTION: Mise à jour des analyses lorsqu'elles sont disponibles, c-a-d
!  lorsque la fenetre de convolution temporelle se clot pour les analyses 
!  classiques et quand l'ensemble des periodes d'analyses des scalograms
!  ont ete traitees. 
!
!> @brief Updates output arrays var2d, var3d with the analysis when available.
!!  Updates scal0d_oa/_cr when a scalogram analyses is completed.
!
!> @details in the case of scalogram analysis, handles the possible 
!!  MPI/tile(-threads) subdomain decomposition (see if_mpi_oa).
!
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014) + Stand-alone version + optimization (01/2015)
!!  - Comments, headers, cleaning, f77 to earlier standards (intent attr.,..)
!!  - Croco-OA interface 1st version (2020)
!!    Interface = mix between the Stand-alone and Original OnlineA versions 
!!  - Toward a Croco-OA tile(-threads) compliant OA version (2021)
!!    supporting tile(-threads) and possible dual OpenMP-MPI parallelization 
!!    of the horizontal domain
!!  - MPI + openMP instructions, use of st(tile)% derived type for
!!    tile(-thread) dependant variables (state vector wf_oa%coef,...)
!!  - XIOS outputs (output_oa, var2d_oa_out, var3d_oa_out,...)
!!  - Developments to perform online Scalogram analyses (06/2021)
!!   - from namelist options to calculation and XIOS netcdf outputs
!!   - scalograms distribution over MPI process/subdomains handled 
!!     either by XIOS or by MPI internal OA instructions (see if_mpi_oa)
!
!> @date 2021
!> @todo BLXD
!!  - Croco 2021 non-cplx OA analyses have not been fully tested 
!!    i.e., options different from swt_wfpf_oa(iv_u)=4 
!------------------------------------------------------------------------------

      subroutine var_upd_oa( tile                                     &
      ,ic_u                                                           &
      ,iv_u                                                           &
      ,lt_u                                                           &
      ,la_u                                                           &
      ,myrank                                               )

      use module_tile_oa, only : st, tile_space_str, ntiles

      use module_oa_variables !updv_oa, kmin3d_oa, begvs3d_oa, begvsoa, l2i_oa, l2j_oa, tupd_oa
      use module_oa_time
      use module_oa_space     !tgv3d_oa, dk_oa
      use module_oa_periode   !swt_wfpf_oa
      !use module_oa_stock     !wf_oa, tallocated_oa
#ifdef MPI
      use module_parameter_oa, only : iminmpi, jminmpi
#endif
      implicit none

      integer, intent(in)  :: tile, ic_u, iv_u, lt_u, la_u 
      integer, intent(out) :: myrank

      !BLXD_TILE_ISSUE rm pst from list of public module variables add tile_space_str
      type(tile_space_str), pointer :: pst => null()

      integer    ::                                                   &
       i_u                                                            &
      ,j_u                                                            &
      ,k_u                                                            &
      ,ls_u                                                           &
      ,ls1_u                                                           
      
      integer, pointer :: pla_u => null(), ps1_u => null()
      integer :: io_nodoa, iu_glob, ju_glob, ierr, itag

 
      variable_commune : if (iv_u.ne.-1.and.updv_oa(iv_u).eq.2) then
!---------------------------------------------------------------------------
!.....mise a jour des variables communes:
!---------------------------------------------------------------------------

! BLXD local sub. variable => openMP tile-threads PRIVATE var.
      myrank = -99

!.....History file:

!$OMP MASTER

      call history_oa(14,lt_u,iv_u,-1,-1, -1)

!$OMP END MASTER

!.....Local pointer to module_tile_oa array structure
      ! BLXD TODO openMP SHARED issue; solution declare pst locally
      pst => st(tile) 

!.....Local pointer with iv_u, lt_u args in
      pla_u => tallocated_oa(lt_u)
      ps1_u => pst%begvs3d_oa(pst%begvs_oa(iv_u))

       var_number_of_dimensions : if (tgv3d_oa(tv_oa(iv_u)).eq.3) then  ! three-dimensional variables

!.......sortie de la partie reelle du coeff
        what_2_extract_3dvar : if (swt_wfpf_oa(iv_u).eq.1) then

!........BLXD For a given iv variable this loop is executed by a single MPI process and a single thread
!        openMP data race cannot occur below this loop 
!        in any of the value update targetting global module variables (implicitly openMP SHARED) 

         do ls_u = pst%begvs_oa(iv_u),pst%begvs_oa(iv_u+1)-1
          i_u = pst%l2i_oa(ls_u)
          j_u = pst%l2j_oa(ls_u)
          do ls1_u = pst%begvs3d_oa(ls_u),pst%begvs3d_oa(ls_u+1)-1
           k_u = pst%kmin3d_oa(ls_u) + (ls1_u-pst%begvs3d_oa(ls_u))* dk_oa(iv_u)
           !var3d_oa(i_u,j_u,k_u,tupd_oa(iv_u)) = real(wf_oa(tallocated_oa(lt_u))%coef (ls1_u-begvs3d_oa(begvs_oa(iv_u))+1 ) )
           if ( .not. tv_sclg_oa(iv_u) )  then
            var3d_oa(i_u,j_u,k_u,tupd_oa(iv_u)) = real( pst%wf_oa(pla_u)%coef( ls1_u-ps1_u+1 ) )
           else

#ifdef MPI

!...........BLXD iv is a scalogram -> unique icslg index -> hold by a unique process -> hold by a unique tile(-thread) 
!           if tile(-thread/process passes here => the tile(-thread) process is concerned by MPI sending the scalogram
!           tile/tile-thread will fill disjoint part of the sacl0d_oa array successively (no data race)
            scal0d_oa(pers_t2p_oa(iv_u,lt_u),v2locsclg(iv_u)) = real( pst%wf_oa(pla_u)%coef( ls1_u-ps1_u+1 ) )
            myrank    = mynode

#else
            scal0d_cr(pers_t2p_oa(iv_u,lt_u),tupd_oa(iv_u)) = real( pst%wf_oa(pla_u)%coef( ls1_u-ps1_u+1 ) )
#endif
           end if 
          enddo
         enddo
        ! BLXD exclusive options if => else if
        !endif

!.......sortie du coef au carre:
        else if (swt_wfpf_oa(iv_u).eq.2) then what_2_extract_3dvar

         do ls_u = pst%begvs_oa(iv_u),pst%begvs_oa(iv_u+1)-1
          i_u = pst%l2i_oa(ls_u)
          j_u = pst%l2j_oa(ls_u)
          do ls1_u = pst%begvs3d_oa(ls_u),pst%begvs3d_oa(ls_u+1)-1
           k_u = pst%kmin3d_oa(ls_u) + (ls1_u-pst%begvs3d_oa(ls_u))* dk_oa(iv_u)
           if ( .not. tv_sclg_oa(iv_u) )  then
            var3d_oa(i_u,j_u,k_u,tupd_oa(iv_u)) = abs( pst%wf_oa(pla_u)%coef( ls1_u-ps1_u+1 ) )**2
           else
#ifdef MPI
            scal0d_oa(pers_t2p_oa(iv_u,lt_u),v2locsclg(iv_u)) = abs( pst%wf_oa(pla_u)%coef( ls1_u-ps1_u+1 ) )**2
            myrank    = mynode
#else
            scal0d_cr(pers_t2p_oa(iv_u,lt_u),tupd_oa(iv_u)) = abs( pst%wf_oa(pla_u)%coef( ls1_u-ps1_u+1 ) )**2
#endif
           endif
          enddo
         enddo
        !endif

!.......sortie du valeur absolue du coef:
        else if (swt_wfpf_oa(iv_u).eq.3) then what_2_extract_3dvar

         do ls_u = pst%begvs_oa(iv_u),pst%begvs_oa(iv_u+1)-1
          i_u = pst%l2i_oa(ls_u)
          j_u = pst%l2j_oa(ls_u)
          do ls1_u = pst%begvs3d_oa(ls_u),pst%begvs3d_oa(ls_u+1)-1
           k_u = pst%kmin3d_oa(ls_u) + (ls1_u-pst%begvs3d_oa(ls_u))* dk_oa(iv_u)
           if (.not. tv_sclg_oa(iv_u) )  then
            var3d_oa(i_u,j_u,k_u,tupd_oa(iv_u)) = abs( pst%wf_oa(pla_u)%coef( ls1_u-ps1_u+1 ) )
           else
#ifdef MPI
            scal0d_oa(pers_t2p_oa(iv_u,lt_u),v2locsclg(iv_u)) = abs( pst%wf_oa(pla_u)%coef( ls1_u-ps1_u+1 ) )
            myrank    = mynode
#else
            scal0d_cr(pers_t2p_oa(iv_u,lt_u),tupd_oa(iv_u)) = abs( pst%wf_oa(pla_u)%coef( ls1_u-ps1_u+1 ) )
#endif
           endif
          enddo
         enddo
        !endif BLXD exclusive swt_wfpf_oa option for a single iv var

!.......sortie du coef complexe:
        else if (swt_wfpf_oa(iv_u).eq.4) then what_2_extract_3dvar

!........BLXD 
!        For a given iv variable the following loop is executed by a single MPI process 
!        and by a single tile(-thread)
!        For Croco tile-thread simu, openMP data race cannot occur in the below the loop
!        regarding any of the global module variables (implicitly openMP SHARED) updated
!        in the ls_u loop 
         do ls_u = pst%begvs_oa(iv_u),pst%begvs_oa(iv_u+1)-1
          i_u = pst%l2i_oa(ls_u)
          j_u = pst%l2j_oa(ls_u)


          do ls1_u = pst%begvs3d_oa(ls_u),pst%begvs3d_oa(ls_u+1)-1
           k_u = pst%kmin3d_oa(ls_u) + (ls1_u-pst%begvs3d_oa(ls_u))* dk_oa(iv_u)

           not_a_scalogram : if ( .not. tv_sclg_oa(iv_u) )  then
             var3d_oa(i_u,j_u,k_u,tupd_oa(iv_u)) = pst%wf_oa(pla_u)%coef( ls1_u-ps1_u+1 )
#if defined OA_TRACES
             if_verb3d : if (verbose_oa>=6) then
#ifdef MPI
                iu_glob = i_u + iminmpi-1 
                ju_glob = j_u + jminmpi-1 
#else
                iu_glob = i_u 
                ju_glob = j_u 
#endif
                 !if ( (lt_u==ltrec_fst_oa(iv_u)) ) then
# ifdef IGW
                if ( iu_glob ==itgt_glob .and. (k_u==ktgt) ) then
# else
#  ifdef MILES
                if ( iu_glob==itgt_glob .and. ju_glob==jtgt_glob .and. k_u==ktgt ) then
#  endif
# endif
                     io_nodoa = set_io_nodoa(iv_u,mynode,4,3)! [3/4](odd/even-iv)000+mynode
                     ! old script iv_u =1 => 30000+mynode iv_u=2 => 40000+mynode
                     write (io_nodoa,fmt='(i5,i5,i3,2(1x,ES22.15E2))')i_u,j_u,k_u &
                         ,REAL(DBLE( var3d_oa(i_u,j_u,k_u,tupd_oa(iv_u)) )),REAL(DIMAG( var3d_oa(i_u,j_u,k_u,tupd_oa(iv_u)) ))
                end if
             endif if_verb3d
#endif
          else not_a_scalogram

#ifdef MPI
!............BLXD iv is a scalogram -> unique icslg index -> hold by a unique process -> hold by a unique tile(-thread) 
!            if (tile-thread)/process passing here => the (tile-thread) process is concerned by MPI sending the scalogram
!            tile/tile-thread will fill disjoint part of the sacl0d_oa array successively (no data race)
             scal0d_oa(pers_t2p_oa(iv_u,lt_u),v2locsclg(iv_u)) = pst%wf_oa(pla_u)%coef( ls1_u-ps1_u+1 )
             myrank    = mynode
#else
             scal0d_cr(pers_t2p_oa(iv_u,lt_u),tupd_oa(iv_u)) = pst%wf_oa(pla_u)%coef( ls1_u-ps1_u+1 )
#endif
#if defined OA_TRACES
             if_verb0d : if (verbose_oa>=5) then
#ifdef MPI
                iu_glob = i_u + iminmpi-1 
                ju_glob = j_u + jminmpi-1 
#else
                iu_glob = i_u 
                ju_glob = j_u 
#endif
                if (iv_u==1) then
                 !if ( (lt_u==ltrec_fst_oa(iv_u)) ) then
                     io_nodoa = set_io_nodoa(pers_t2p_oa(iv_u,lt_u),mynode,4,3)! [3/4](odd/even-iv)000+mynode
                     ! old script iv_u =1 => 30000+mynode iv_u=2 => 40000+mynode
#ifdef MPI
            write (io_nodoa,fmt='(i5,i5,i3,2(1x,ES22.15E2))')i_u,j_u,k_u &
                ,REAL(DBLE( scal0d_oa(pers_t2p_oa(iv_u,lt_u),v2locsclg(iv_u)) )),REAL(DIMAG( scal0d_oa(pers_t2p_oa(iv_u,lt_u),v2locsclg(iv_u)) ))
#else
            write (io_nodoa,fmt='(i5,i5,i3,2(1x,ES22.15E2))')i_u,j_u,k_u &
                ,REAL(DBLE( scal0d_cr(pers_t2p_oa(iv_u,lt_u),tupd_oa(iv_u)) )),REAL(DIMAG( scal0d_cr(pers_t2p_oa(iv_u,lt_u),tupd_oa(iv_u)) ))
#endif
                endif
             endif if_verb0d
#endif

          endif not_a_scalogram
          enddo
         enddo

        endif what_2_extract_3dvar

       else var_number_of_dimensions                                                          ! two-dimensional variables

!.......sortie de la partie reelle du coef
        what_2_extract_2dvar : if (swt_wfpf_oa(iv_u).eq.1) then

         do ls_u = pst%begvs_oa(iv_u),pst%begvs_oa(iv_u+1)-1
          i_u = pst%l2i_oa(ls_u)
          j_u = pst%l2j_oa(ls_u)
          do ls1_u = pst%begvs3d_oa(ls_u),pst%begvs3d_oa(ls_u+1)-1
           if ( .not. tv_sclg_oa(iv_u) )  then
            var2d_oa(i_u,j_u,tupd_oa(iv_u)) = real(pst%wf_oa(pla_u)%coef(ls1_u-ps1_u+1 ) )
           else
#ifdef MPI
            scal0d_oa(pers_t2p_oa(iv_u,lt_u),v2locsclg(iv_u)) = real(pst%wf_oa(pla_u)%coef( ls1_u-ps1_u+1 ))
            myrank    = mynode
#else
            scal0d_cr(pers_t2p_oa(iv_u,lt_u),tupd_oa(iv_u)) = real( pst%wf_oa(pla_u)%coef( ls1_u-ps1_u+1 ) )

#endif
           endif
          enddo
         enddo
        ! BLXD exclusive options if => else if
        !endif

!.......sortie du coef au carre (energy):
        else if (swt_wfpf_oa(iv_u).eq.2) then what_2_extract_2dvar

         do ls_u = pst%begvs_oa(iv_u),pst%begvs_oa(iv_u+1)-1
          i_u = pst%l2i_oa(ls_u)
          j_u = pst%l2j_oa(ls_u)
          do ls1_u = pst%begvs3d_oa(ls_u),pst%begvs3d_oa(ls_u+1)-1
           if ( .not. tv_sclg_oa(iv_u) )  then
            var2d_oa(i_u,j_u,tupd_oa(iv_u)) = abs(pst%wf_oa(pla_u)%coef (ls1_u-ps1_u+1 ) )**2
           else
#ifdef MPI
            scal0d_oa(pers_t2p_oa(iv_u,lt_u),v2locsclg(iv_u)) = abs(pst%wf_oa(pla_u)%coef (ls1_u-ps1_u+1 ) )**2
            myrank    = mynode
#else
            scal0d_cr(pers_t2p_oa(iv_u,lt_u),tupd_oa(iv_u)) = abs(pst%wf_oa(pla_u)%coef (ls1_u-ps1_u+1 ) )**2

#endif
           endif
          enddo
         enddo
        !endif

!.......sortie de la valeur absolue du coef: 
        else if (swt_wfpf_oa(iv_u).eq.3) then what_2_extract_2dvar

         do ls_u = pst%begvs_oa(iv_u),pst%begvs_oa(iv_u+1)-1  
          i_u = pst%l2i_oa(ls_u)
          j_u = pst%l2j_oa(ls_u)
          do ls1_u = pst%begvs3d_oa(ls_u),pst%begvs3d_oa(ls_u+1)-1
           if ( .not. tv_sclg_oa(iv_u) )  then
            var2d_oa(i_u,j_u,tupd_oa(iv_u)) = abs(pst%wf_oa(pla_u)%coef(ls1_u-ps1_u+1 ) )
           else
#ifdef MPI
            scal0d_oa(pers_t2p_oa(iv_u,lt_u),v2locsclg(iv_u)) = abs(pst%wf_oa(pla_u)%coef(ls1_u-ps1_u+1 ) )
            myrank    = mynode
#else
            scal0d_cr(pers_t2p_oa(iv_u,lt_u),tupd_oa(iv_u)) = abs(pst%wf_oa(pla_u)%coef(ls1_u-ps1_u+1 ) )

#endif
           endif
          enddo
         enddo
        !endif

!.......sortie du coeff complexe:
        else if (swt_wfpf_oa(iv_u).eq.4) then what_2_extract_2dvar

         do ls_u = pst%begvs_oa(iv_u),pst%begvs_oa(iv_u+1)-1  
          i_u = pst%l2i_oa(ls_u)
          j_u = pst%l2j_oa(ls_u)
          do ls1_u = pst%begvs3d_oa(ls_u),pst%begvs3d_oa(ls_u+1)-1
           not_a_scalogram_2dvar : if ( .not. tv_sclg_oa(iv_u) )  then
            var2d_oa(i_u,j_u,tupd_oa(iv_u)) = (pst%wf_oa(pla_u)%coef (ls1_u-ps1_u+1 ) )

#if defined OA_TRACES
            if_verb2d : if (verbose_oa>=6) then
#ifdef MPI
                iu_glob = i_u + iminmpi-1 
                ju_glob = j_u + jminmpi-1 
#else
                iu_glob = i_u 
                ju_glob = j_u 
#endif
                !if ( (lt_u==ltrec_fst_oa(iv_u)) ) then
                if ( iu_glob ==itgt_glob .and. (ju_glob==jtgt_glob) ) then
                    io_nodoa = set_io_nodoa(iv_u,mynode,4,3)! [3/4](odd/even-iv)000+mynode
                    ! old script iv_u =1 => 30000+mynode iv_u=2 => 40000+mynode
                    write (io_nodoa,fmt='(i5,i5,2(1x,ES22.15E2))')i_u,j_u &
                        ,REAL(DBLE( var2d_oa(i_u,j_u,tupd_oa(iv_u)) )),REAL(DIMAG( var2d_oa(i_u,j_u,tupd_oa(iv_u)) ))
                end if
            endif if_verb2d
#endif
           else not_a_scalogram_2dvar
#ifdef MPI
            scal0d_oa(pers_t2p_oa(iv_u,lt_u),v2locsclg(iv_u)) = (pst%wf_oa(pla_u)%coef (ls1_u-ps1_u+1 ) )
            myrank    = mynode
#else
            scal0d_cr(pers_t2p_oa(iv_u,lt_u),tupd_oa(iv_u)) = (pst%wf_oa(pla_u)%coef (ls1_u-ps1_u+1 ) )

#endif
           endif not_a_scalogram_2dvar
          enddo
         enddo

        endif what_2_extract_2dvar

        endif var_number_of_dimensions

#ifdef MPI

        if_oa_handles_sclg_mpi : if ( if_mpi_oa ) then

!.....If scalogram : MPI Send MPI subdomain array scal0d_oa for var iv_u <-> one scal.
        if_sclg : if ( tv_sclg_oa(iv_u) )  then

!.....If range of period completed and if all tiles completed, if MPI

         if_nper_sclg_completed : if ( per_t2p_oa(lt_u) == (begvp_oa(iv_u+1) - 1) ) then

!$OMP CRITICAL (send_sclg_oa)
! BLXD
!     The subroutine variable myrank is a local PRIVATE thread variable :
!     myrank = -99                     => local tile(-thread) of the current MPI process has nothing to send
!     myrank = rank of the MPI process => local tile(-thread) of thr current MPI process has scalogram data to send
!     Note that all the tile-thread must test if mynode=myrank 
!     Do not restrict MPI communication to the last thread  here
!       i.e., last_tile : if ( tile_count_oa .eq. ntiles) then 
!
! Either Croco sequential tiles or Croco tile-threads (with openMP CRITICAL REGION) 
! will process the following instructions one after one 

! BLXD TODO Check if openMP CRITICAL needed ? 
!           An openMP BARRIER at the end of condition mynode==myrank seems sufficiant 

           if_node_correct : if ( mynode == myrank ) then

            if_scalogram_not_on_root : if ( myrank /= root ) then
                ! MPI subprocess scalogram iv_u in buffer to send to root process
                ! with tag set to the global scalogram index
                itag = tupd_oa(iv_u)
                buffx_s(1:nper_sclg(iv_u)) => scal0d_oa(1:nper_sclg(iv_u),v2locsclg(iv_u))
                bs_sclg : if ( mpi_nonblock_send) then
                    call MPI_ISend (buffx_s, nper_sclg(iv_u), MPI_DOUBLE_COMPLEX, root, itag, comm, sclg_request(iv_u), ierr)
                    if (ierr/=0) then
                        print*,'ERROR OA : NON BLOCKING MPI_ISend scalogram for mynode-var returns an error ',mynode,iv_u,ierr
                        stop
                    endif
                else bs_sclg
                    !call MPI_Send (buffx_s, nper_sclg_max, MPI_DOUBLE_COMPLEX, root, itag, comm, ierr)
                    call MPI_Send (buffx_s, nper_sclg(iv_u), MPI_DOUBLE_COMPLEX, root, itag, comm, ierr)
                    if (ierr/=0) then
                        print*,'ERROR OA : NON BLOCKING MPI_ISend scalogram for mynode-var returns an error ',mynode,iv_u,ierr
                        stop
                    endif
                endif bs_sclg
            else if_scalogram_not_on_root
               ! Root rank process holds the saclogram then no MPI exchange !!!
               ! var_upd_oa returns myrank /= -99
                scal0d_cr(1:nper_sclg(iv_u),tupd_oa(iv_u))=scal0d_oa(1:nper_sclg(iv_u),v2locsclg(iv_u))
            endif if_scalogram_not_on_root
           endif if_node_correct

!$OMP END CRITICAL (send_sclg_oa)

         endif if_nper_sclg_completed
        endif if_sclg 

        endif if_oa_handles_sclg_mpi
#endif /* MPI */

!.....Nullify pointer to module_oa_tile struct.
       
      pst => null()

!.....Nullify local pointer
      pla_u => null()
      ps1_u => null()

      endif variable_commune

      return
      end subroutine var_upd_oa


!----------------------------------------------------------------------
! PROCEDURE
!
!> @note This routine is called under a Croco tile(-thread)
!!       loop (see routine main_oa and croco sources main.F/step.F).
!
! DESCRIPTION: Allocation et desallocation des champs du type derivee
!  wf_oa (coef d'analyse en ondelette, spectrale,...etc) au gre des 
!  besoins d'ouverture/fermeture des fenetres de convolution temporelle.
!
!
!> @brief Allocates/deallocates array fields of derived type wf_oa
!
!> @details  wf_oa stores all the online analyses (wavelet,spectral coef.)
!!   before using them or dumping them for outputs in var2d/var3d/scal 
!!   variables.
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Toward a Croco-OA tile(-threads) compliant interface and OA version (2021)
!!  - Memory optimization with the pre-calculation of the size of wf_oa array,
!!    i.e., the max number of convolution windows simultaneously openened
!!    during the simulation (# of conv. wind. are namelist-dependent).
!!    This max number of conv. wind can be overwritten using the namelist parameter
!!    nmsimult_oa_max (default value is -99 => online pre-calculation). 
!> @date 2021
!> @todo BLXD
!!  - general requirement to organize precision consistency
!!    eventhough compiling with double prec option
!!    Module with wp declaration and use of KIND(wp) ?
!----------------------------------------------------------------------

      subroutine allocate_win_oa (tile, lti_a, lc_a, lv_a, dim_a, iic_oa, l_a)

      use module_tile_oa , only     : st, tile_space_str
      use module_oa_time , only : nmsimult_oa, tallocated_oa
      implicit none

      integer, intent(in) :: tile
      integer, intent(in) :: iic_oa

      integer, intent(in) :: &
        lc_a                 &  !< index related to a specific configuration of analysis 
       ,lv_a                 &  !< index related to a specific variable to analyse
       ,lti_a                   !< index related to a specific time of analysis (i.e., time convolution window)

      integer, intent(in)  :: dim_a !< Size of the spatial domain to analyse for variable lv_a 
      integer, intent(out) :: l_a

      type(tile_space_str), pointer :: pst => null()


      pst => st(tile)

      l_a=1 ! Shall we keep track of imsimult_oa ?
      do while ( (associated(pst%wf_oa(l_a)%coef)).and.l_a.lt.nmsimult_oa)
         l_a = l_a + 1
      enddo

      if (l_a.eq.nmsimult_oa ) then
!$OMP MASTER
         if(if_print_node) write (io_unit,*) 'Pre-calculated size for wf_oa now reached ',nmsimult_oa 
!$OMP END MASTER
      else if (l_a.gt.nmsimult_oa) then
         write (*,*) 'ERROR OA : above pre-calculated size for wf_oa ',nmsimult_oa
         if(if_print_node) write (io_unit,*) 'STOP : above pre-calculated size for wf_oa ',nmsimult_oa
         stop
      endif

      if (l_a.eq.nmsimult_oa.and.associated(pst%wf_oa(l_a)%coef) ) then
         write (*,*) 'ERROR OA : nmsimult_oa trop petit ! ', nmsimult_oa
         if(if_print_node) write (io_unit,*) 'ERROR OA : nmsimult_oa trop petit ! ', nmsimult_oa
         if(if_print_node) write (io_unit,*) 'si simu avec restart (to dev) fixez nmsimult_oa_max via la namelist_oa'
         stop
      endif

      ! BLXD TODO tile tile-thread loop => protect against MPI process double allocation
      ! Check if .not. associated instead allocated
      if (.not. associated(pst%wf_oa(l_a)%coef) ) allocate( pst%wf_oa(l_a)%coef(dim_a) )

!.....History file:
!$OMP MASTER
      call history_oa(2,lc_a,lv_a,l_a,lti_a,-1, iic_oa)
!$OMP END MASTER

      pst%wf_oa(l_a)%coef(:)               = (0.D0,0.D0) ! #BLDX double prec
      tallocated_oa(lti_a)                 = l_a
      pst%wf_oa(l_a)%t_indice              = lti_a 
      pst%wf_oa(l_a)%config                = lc_a 
      pst%wf_oa(l_a)%variable              = lv_a 

      pst => null()

      return
      end subroutine allocate_win_oa


      subroutine deallocate_win_oa ( tile, lt_a )    

      use module_tile_oa , only     : st, tile_space_str
      use module_oa_time , only : nmsimult_oa, tallocated_oa

      implicit none

      integer, intent(in) :: tile
      integer, intent(in) :: lt_a

      !BLXD_TILE_ISSUE rm pst from list of public module variables add tile_space_str
      type(tile_space_str), pointer :: pst => null()   
      
!.....History file:
!$OMP MASTER
      call history_oa(3,lt_a,-1,-1,-1, -1)
!$OMP END MASTER

      ! BLXD AUG2021 coef is an "allocated pointer"
      ! does not point twd a target
      ! deallocate (wf_oa(tallocated_oa(lt_a))%coef )


      pst => st(tile)

      deallocate (pst%wf_oa(tallocated_oa(lt_a))%coef )
      nullify( pst%wf_oa(tallocated_oa(lt_a))%coef )

      pst%wf_oa(tallocated_oa(lt_a))%t_indice   = -1 
      pst%wf_oa(tallocated_oa(lt_a))%config     = -1
      pst%wf_oa(tallocated_oa(lt_a))%variable   = -1 
      tallocated_oa(lt_a)                       = -1

      pst => null()

      return
      end subroutine deallocate_win_oa 

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
!
! DESCRIPTION: 
!
!> @brief returns the complex Morlet wavelet or windowed Fourier.
!
!> @details 
!
!
! REVISION HISTORY:
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - f77 to earlier fortran standard (intent attr., goto,...) 
!!  - Stand-alone (2015)
!> @date 2021
!> @todo BLXD
!   check precision requirements (cplx,..,)
!------------------------------------------------------------------------------
      complex function psi_oa(                                        & 
                 tpsi_p                                               & 
                ,scale_p                                              & 
                ,t_p                                                  & 
                ,dti_p                                                & 
                ,fb_p                                                 & 
                ,fc_p             )


      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      !use module_oa_stock

!....attention il faut encore tenir compte des normalisations...

      integer, intent(in) :: tpsi_p

      double precision, intent(in) :: dti_p, t_p

      real, intent(in)  ::                                            &
        scale_p                                                       & 
       ,fb_p                                                          & 
       ,fc_p
      
      real ::                                                         & 
        w0_p                                                          & 
       ,rec_p 

! attention exemple de psi incomplet pour windowed fourier...
! les coefs pour la reconstruction n'ont
! pas ete ajoutes.

!.....dirac:

      if (tpsi_p .eq. 0 ) then      
       psi_oa = 1.

!.....ondelette de morlet:
  
      else if (tpsi_p .eq. 1 ) then      
       psi_oa =                                                        & 
          exp( - ( t_p /scale_p ) ** 2 / fb_p )                        & 
        * exp( t_p / scale_p * fc_p*(2.*pi_oa) * (0.,1.) )             & 
!       / (2*pi_oa)**0.5                                               & 
        / (pi_oa)**0.25                                                & 
!       * sqrt( dti_p/ scale_p )    
        / sqrt( scale_p )    

!....."windowed fourier":
      else if (tpsi_p .eq. 2.or.tpsi_p.eq.3) then
       psi_oa =  exp(  2.* pi_oa* ( t_p - dti_p ) / scale_p * (0.,1.) )
      endif 

      end function psi_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note toward f2008 standard with object oriented structure (class, method)
!
!
! DESCRIPTION: 
!
!> @brief returns the complex Morlet wavelet value.
!
!> @details
!
!
! REVISION HISTORY:
!> @authors
!! - B. Lemieux-Dudon
!> @date 2021
!> @todo BLXD
!   check double precision, rm tpsi_p arg 
!------------------------------------------------------------------------------
      complex function morlet_psi_oa(                                  & 
                 tpsi_p                                                & 
                ,scale_p                                               & 
                ,t_p                                                   & 
                ,dti_p                                                 & 
                ,fb_p                                                  & 
                ,fc_p             )


      ! BLXD Tests on atom flag might be removed after tests
      integer, intent(in) :: tpsi_p

      ! useful for psi_oa_morlet FT
      double precision, intent(in) :: dti_p

      double precision, intent(in) :: t_p

      real, intent(in)  ::                                            &
        scale_p                                                       & 
       ,fb_p                                                          & 
       ,fc_p
!
!      real ::                                                         & 
!        w0_p                                                          & 
!       ,rec_p 


!.....ondelette de morlet seulement :
!     BLXD WARNING about the psi_oa fb_p parameter and the normalization
!     1) If fb_p = 2 => sigma =1. and everything is fine
!     2) If changing the psi_oa fb_p parameter 
!       => changing sigma values such as sigma^2 = fb_p / 2.
!     In case 2) the psi_oa formulas will have to be modified :   
!     - The C normalization factor being 
!       C = 1. / sqrt(sigma) / pi^(1./4.) = 1. / ( sigma^2 * pi )^(1./4.)
!       which writes in terms of fb_p
!       C = 1. / ( fb_p * pi / 2. )^(1./4.) 
!       => REPLACE / (pi_oa)**0.25 by / ( fb_p * pi_oa / 2. )**0.25 

  
!      if (tpsi_p .eq. 1 ) then      
       morlet_psi_oa =                                                 & 
          exp( - ( t_p /scale_p )** 2 / fb_p )                         & 
        * exp( t_p / scale_p * fc_p*(2.*pi_oa) * (0.,1.) )             & 
!       / (2*pi_oa)**0.5                                               & 
        / (pi_oa)**0.25                                                & 
!       * sqrt( dti_p/ scale_p )    
        / sqrt( scale_p )    

!      else if (tpsi_p .eq. 0 .or. tpsi_p .eq. 2 .or. tpsi_p.eq.3 ) then
!      if(if_print_node)write(io_unit,*) "ERROR : atom flag tpsi_p must correspond to Morlet wavelet" 
!      stop
!     else 
!      if(if_print_node)write(io_unit,*) "ERROR : unknown atom flag tpsi_p" 
!      stop
!     endif 

      end function morlet_psi_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note toward f2008 standard with object oriented structure (class, method)
!
!
! DESCRIPTION: 
!
!> @brief several Morlet wavelet results.
!
!> @details 
!
!
! REVISION HISTORY:
!> @authors
!! - B. Lemieux-Dudon
!> @date 2021
!> @todo BLXD function related to Morlet wavelet should be grouped in a Morelt
!        module with class and method 
!   check double precision 
!------------------------------------------------------------------------------


      real function morlet_psi_sigma2_oa(                          & 
                 fb_p                                              )
        implicit none
        real, intent(in) :: fb_p
        if (fb_p==2.) then
            morlet_psi_sigma2_oa = 1.
        else
            morlet_psi_sigma2_oa = fb_p/2.
        endif
      end function morlet_psi_sigma2_oa

      real function morlet_psi_normalization_oa(                   & 
                 sigma2                                            )
        implicit none
        real, intent(in) :: sigma2
        if (sigma2==1.) then
            morlet_psi_normalization_oa = 1. / ( pi_oa )**0.25
        else
            morlet_psi_normalization_oa = 1. / ( pi_oa * sigma2 )**0.25
        endif
      end function morlet_psi_normalization_oa

      real function morlet_psi_w0c_oa(                               & 
                 fc_p                                                &
                ,scale_p )   ! use morlet_psi_p2s_oa here
        implicit none
        real, intent(in) :: fc_p
        real, intent(in), optional :: scale_p

         if ( present(scale_p) ) then

            morlet_psi_w0c_oa  = fc_p*2.*pi_oa/scale_p   ! theo freq box oa  center
                                                        ! eta/a
         else 
                                                        ! mother wavelet theo freq center
            morlet_psi_w0c_oa  = fc_p*2.*pi_oa           ! eta

         endif

      end function morlet_psi_w0c_oa

      real function morlet_psi_dw0_oa(                              & 
                 fb_p                                               &
                ,fc_p                                               &
                ,scale_p                                            &   ! use morlet_psi_p2s_oa here
                ,from_dt0_eval                                      )
        implicit none
        real, intent(in) :: fb_p, fc_p
        real, intent(in), optional :: scale_p
        real*8, intent(in), optional :: from_dt0_eval
        real :: sigma

        if ( present(from_dt0_eval) ) then
            morlet_psi_dw0_oa = 1. / ( 2. * from_dt0_eval )
        else
         if ( present(scale_p) ) then
            if (fb_p==2.) then
                !sigma=1 not useful
                morlet_psi_dw0_oa=1./sqrt(2.)/scale_p
            else
                sigma = sqrt( morlet_psi_sigma2_oa(fb_p) )
                morlet_psi_dw0_oa=1./sqrt(2.)/scale_p/sigma            
            endif
         else
             if (fb_p==2.) then
                !sigma=1 not useful
                 morlet_psi_dw0_oa=1./sqrt(2.)
            else
                sigma = sqrt( morlet_psi_sigma2_oa(fb_p) )
                morlet_psi_dw0_oa=1./sqrt(2.)/sigma            
            endif
         endif           
        endif

      end function morlet_psi_dw0_oa

      real function morlet_psi_dt0_oa(                              & 
                 fb_p                                               &
                ,fc_p                                               &
                ,scale_p                                            &! use morlet_psi_p2s_oa here
                ,from_dw0_eval                                      )

        implicit none
        real, intent(in) :: fb_p, fc_p
        real, intent(in), optional :: scale_p
        real, intent(in), optional :: from_dw0_eval
        real :: sigma

        if ( present(from_dw0_eval) ) then
            morlet_psi_dt0_oa = 1. / ( 2. * from_dw0_eval )
        else
         if ( present(scale_p) ) then
            if (fb_p==2.) then
                !sigma=1 not useful
                morlet_psi_dt0_oa=scale_p/sqrt(2.)
            else
                sigma = sqrt( morlet_psi_sigma2_oa(fb_p) )
                morlet_psi_dt0_oa=scale_p*sigma/sqrt(2.)
            endif
         else
             if (fb_p==2.) then
                !sigma=1 not useful
                morlet_psi_dt0_oa=1./sqrt(2.)
            else
                sigma = sqrt( morlet_psi_sigma2_oa(fb_p) )
                morlet_psi_dt0_oa=sigma/sqrt(2.)            
            endif
         endif           
        endif
      end function morlet_psi_dt0_oa

      real function morlet_psi_p2s_oa(                              & 
                 per_p                                              &
                ,fb_p                                               &
                ,fc_p                                               )
        implicit none
        real, intent(in) :: per_p, fb_p, fc_p
        real :: sigma_p, wc_p

       if (fb_p==2.) then
            ! sigma_p = 1.
            ! identical to function psi_p2s_oa if wavelet tpsi_p == 1 
            morlet_psi_p2s_oa = per_p / ( 4. * pi_oa ) * (fc_p*2.*pi_oa + sqrt(2.+(fc_p*2.*pi_oa)**2))  
       else
            sigma_p  = sqrt( fb_p / 2. )
            wc_p   = fc_p*2.*pi_oa      ! eta
            morlet_psi_p2s_oa = per_p / ( 4. * pi_oa * sigma_p ) * ( sigma_p * wc_p + sqrt( 2.+ sigma_p**2 * wc_p**2 ) )  
       endif
      end function morlet_psi_p2s_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
!
! DESCRIPTION: 
!
!> @brief returns the complex Morlet wavelet value.
!
!> @details 
!
!
! REVISION HISTORY:
!> @authors
!! - B. Lemieux-Dudon
!> @date 2021
!> @todo BLXD
!   check double precision 
!------------------------------------------------------------------------------
      complex function FT_morlet_psi_oa(                               & 
                 tpsi_p                                                & 
                ,u_p                                                   &
                ,scale_p                                               & 
                ,w_p                                                   & 
                ,dwi_p                                                 & 
                ,fb_p                                                  & 
                ,fc_p             )


      ! BLXD Tests on atom flag might be removed after tests
      integer, intent(in) :: tpsi_p

      ! useful for psi_oa_morlet FT
      double precision, intent(in) :: dwi_p

      double precision, intent(in) :: w_p
      real*8          , intent(in) :: u_p

      real, intent(in)  ::                                            &
        scale_p                                                       & 
       ,fb_p                                                          & 
       ,fc_p

!      real ::                                                         & 
!        w0_p                                                          & 
!       ,rec_p 


!.....ondelette de morlet seulement :
!     BLXD WARNING about the psi_oa fb_p parameter and the normalization
!     1) If fb_p = 2 => sigma =1. and everything is fine
!     2) If changing the psi_oa fb_p parameter 
!       => changing sigma values such as sigma^2 = fb_p / 2.
!     In case 2) the psi_oa formulas will have to be modified :   
!     - The C normalization factor being 
!       C = 1. / sqrt(sigma) / pi^(1./4.) = 1. / ( sigma^2 * pi )^(1./4.)
!       which writes in terms of fb_p
!       C = 1. / ( fb_p * pi / 2. )^(1./4.) 
!       => REPLACE / (pi_oa)**0.25 by / ( fb_p * pi_oa / 2. )**0.25 
!     PAY ATTENTION TO Inverse Transforme :
!     - normalization 1. / (2*pi)
!     - summation over p= -N/2 to N/2 covering 1/dti with N*dti = 2.*DELTA_T_W*PER_W = 2*delta_oa_w*per_oa
!     - step dwi = 1. / (N*dti) = 1. / (2*delta_oa_w*per_oa) 
 
!      if (tpsi_p .eq. 1 ) then      
       FT_morlet_psi_oa =                                                    & 
          exp( - fb_p * ( scale_p * w_p - 2. * pi_oa * fc_p )**2 / 4.   )    & 
        * exp( u_p * w_p  * (0.,1.) )                                        & 
        * (2*pi_oa*fb_p)**0.25                                               & 
        * sqrt( scale_p )    



!      else if (tpsi_p .eq. 0 .or. tpsi_p .eq. 2 .or. tpsi_p.eq.3 ) then
!      if(if_print_node)write(io_unit,*) "ERROR : atom flag tpsi_p must correspond to Morlet wavelet" 
!      stop
!     else 
!      if(if_print_node)write(io_unit,*) "ERROR : unknown atom flag tpsi_p" 
!      stop
!     endif 

      end function FT_morlet_psi_oa


!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
!
! DESCRIPTION: 
!
!> @brief returns the module of the Morlet wavelet value.
!
!> @details two calculation ways to check
!
!
! REVISION HISTORY:
!> @authors
!! - B. Lemieux-Dudon
!> @date 2021
!> @todo BLXD
!   check double precision 
!------------------------------------------------------------------------------
      real function module_morlet_psi_oa(                              & 
                 tpsi_p                                                & 
                ,scale_p                                               & 
                ,t_p                                                   & 
                ,dti_p                                                 & 
                ,fb_p                                                  & 
                ,fc_p                                                  &
                ,from_cplx_fun                                         )


      ! BLXD Tests on atom flag might be removed after tests
      integer, intent(in) :: tpsi_p

      ! useful for psi_oa_morlet FT
      double precision, intent(in) :: dti_p

      double precision, intent(in) :: t_p

      logical, intent(in) :: from_cplx_fun

      real, intent(in)  ::                                            &
        scale_p                                                       & 
       ,fb_p                                                          & 
       ,fc_p

! BLXD TODO : what precision ? 
      complex :: psi_p

!
!      real ::                                                         & 
!        w0_p                                                          & 
!       ,rec_p 


!.....ondelette de morlet seulement :
!     BLXD WARNING about the psi_oa fb_p parameter and the normalization
!     1) If fb_p = 2 => sigma =1. and everything is fine
!     2) If changing the psi_oa fb_p parameter 
!       => changing sigma values such as sigma^2 = fb_p / 2.
!     In case 2) the psi_oa formulas will have to be modified :   
!     - The C normalization factor being 
!       C = 1. / sqrt(sigma) / pi^(1./4.) = 1. / ( sigma^2 * pi )^(1./4.)
!       which writes in terms of fb_p
!       C = 1. / ( fb_p * pi / 2. )^(1./4.) 
!       => REPLACE / (pi_oa)**0.25 by / ( fb_p * pi_oa / 2. )**0.25 

  
!      if (tpsi_p .eq. 1 ) then      
       if ( from_cplx_fun ) then

          psi_p = morlet_psi_oa(  &  !                             
                        tpsi_p    &  ! tpsi_oa(iv_m)                            
                       ,scale_p   &  !,psi_p2s_oa( tpsi_oa(iv_m),perv_oa(1,per_t2p_oa(lt_m)),fb_oa,fc_oa ) 
                       ,t_p       &  !,real(kpt_m)* dti                          
                       ,dti_p     &  !,dti*resv_oa(per_t2p_oa(lt_m))             
                       ,fb_p      &  !,fb_oa 
                       ,fc_p      )  !,fc_oa  ) 
                                                                             
          ! conv complex to real type prec module needed !!!
          module_morlet_psi_oa = real( psi_p * conjg(psi_p) )

        else

          ! BLXD : removing imaginary part ; abs value shouldn't be useful since factors are all positive
          module_morlet_psi_oa =                                 &
               exp( - 2.* ( t_p /scale_p ) ** 2 / fb_p )         & 
             / (pi_oa)**0.5                                      & 
             / scale_p    

        endif

!      else if (tpsi_p .eq. 0 .or. tpsi_p .eq. 2 .or. tpsi_p.eq.3 ) then
!      if(if_print_node)write(io_unit,*) "ERROR : atom flag tpsi_p must correspond to Morlet wavelet" 
!      stop
!     else 
!      if(if_print_node)write(io_unit,*) "ERROR : unknown atom flag tpsi_p" 
!      stop
!     endif 

      end function module_morlet_psi_oa

!------------------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
!
! DESCRIPTION: 
!
!> @brief returns the module complex Morlet wavelet value.
!
!> @details 
!
!
! REVISION HISTORY:
!> @authors
!! - B. Lemieux-Dudon
!> @date 2021
!> @todo BLXD
!   check double precision 
!------------------------------------------------------------------------------
      real function module_FT_morlet_psi_oa(                                  & 
                 tpsi_p                                                & 
                ,u_p                                                   &
                ,scale_p                                               & 
                ,w_p                                                   & 
                ,dwi_p                                                 & 
                ,fb_p                                                  & 
                ,fc_p                                                  & 
                ,from_cplx_fun                                         )



      ! BLXD Tests on atom flag might be removed after tests
      integer, intent(in) :: tpsi_p

      ! useful for psi_oa_morlet FT
      double precision, intent(in) :: dwi_p

      double precision, intent(in) :: w_p, u_p

      real, intent(in)  ::                                            &
        scale_p                                                       & 
       ,fb_p                                                          & 
       ,fc_p

      logical, intent(in) :: from_cplx_fun

      complex :: psi_p


!.....ondelette de morlet seulement :
!     BLXD WARNING about the psi_oa fb_p parameter and the normalization
!     1) If fb_p = 2 => sigma =1. and everything is fine
!     2) If changing the psi_oa fb_p parameter 
!       => changing sigma values such as sigma^2 = fb_p / 2.
!     In case 2) the psi_oa formulas will have to be modified :   
!     - The C normalization factor being 
!       C = 1. / sqrt(sigma) / pi^(1./4.) = 1. / ( sigma^2 * pi )^(1./4.)
!       which writes in terms of fb_p
!       C = 1. / ( fb_p * pi / 2. )^(1./4.) 
!       => REPLACE / (pi_oa)**0.25 by / ( fb_p * pi_oa / 2. )**0.25 
!     PAY ATTENTION TO Inverse Transforme :
!     - normalization 1. / (2*pi)
!     - summation over p= -N/2 to N/2 covering 1/dti with N*dti = 2.*DELTA_T_W*PER_W = 2*delta_oa_w*per_oa
!     - step dwi = 1. / (N*dti) = 1. / (2*delta_oa_w*per_oa) 
 
!      if (tpsi_p .eq. 1 ) then      

       if ( from_cplx_fun ) then

          psi_p = FT_morlet_psi_oa(  &  !                             
                        tpsi_p    &  ! tpsi_oa(iv_m)  
                       ,u_p       &                          
                       ,scale_p   &  !,psi_p2s_oa( tpsi_oa(iv_m),perv_oa(1,per_t2p_oa(lt_m)),fb_oa,fc_oa ) 
                       ,w_p       &  !,real(kpt_m)* dti                          
                       ,dwi_p     &  !,dti*resv_oa(per_t2p_oa(lt_m))             
                       ,fb_p      &  !,fb_oa 
                       ,fc_p      )  !,fc_oa  ) 

          ! conv cpx to real type prec module needed
          module_FT_morlet_psi_oa = real( psi_p * conjg(psi_p) )

        else
       
         module_FT_morlet_psi_oa =                                           & 
          exp( - fb_p * ( scale_p * w_p - 2. * pi_oa * fc_p ) ** 2 / 2. )    & 
        !* exp( u_p * w_p  * (0.,1.) )                                        & 
        * (2*pi_oa*fb_p)**0.5                                                & 
        * scale_p    


        endif

!      else if (tpsi_p .eq. 0 .or. tpsi_p .eq. 2 .or. tpsi_p.eq.3 ) then
!      if(if_print_node)write(io_unit,*) "ERROR : atom flag tpsi_p must correspond to Morlet wavelet" 
!      stop
!     else 
!      if(if_print_node)write(io_unit,*) "ERROR : unknown atom flag tpsi_p" 
!      stop
!     endif 

      end function module_FT_morlet_psi_oa


      subroutine moments_of_Morlet_Power_Spectrum(                     & 
                 tpsi_p                                                &
                ,m0_p                                                  &
                ,m1_p                                                  &
                ,m2_p                                                  &
                ,u_p                                                   &
                ,scale_p                                               & 
                ,dwi_p                                                 & 
                ,fb_p                                                  & 
                ,fc_p                                                  &
                ,nw                                                 )
    


      ! BLXD Tests on atom flag might be removed after tests

      integer, intent(in) :: tpsi_p
      integer, intent(in) :: nw
      real,    intent(in) ::                                          &
        scale_p                                                       & 
       ,fb_p                                                          & 
       ,fc_p

      double precision, intent(in) :: dwi_p, u_p

      real, intent(inout) :: m0_p, m1_p, m2_p

      double precision  :: sum0_p, sum1_p, sum2_p
      double precision  :: s1_p, s2_p
      real*8            :: w_p
      real              :: mps_p
      integer :: kw
      logical, parameter :: from_cplx_fun=.true.

!      real ::                                                         & 
!        w0_p                                                          & 
!       ,rec_p 


!.....ondelette de morlet seulement :
!      if (tpsi_p .eq. 1 ) then      

           sum0_p = 0.d0
           sum1_p = 0.d0
           sum2_p = 0.d0

           do kw=1,nw
           ! even nw : -nw/2 to nw/2 integral -W/2 -> W/2
               w_p = real(kw-int(nw/2)) * dwi_p
                
               s1_p = w_p - morlet_psi_w0c_oa( fc_p, scale_p )
               s2_p = s1_p**2

               mps_p = module_FT_morlet_psi_oa(        &  !                             
                                tpsi_p                 &  ! tpsi_oa(iv_m)  
                               ,u_p                    &                          
                               ,scale_p                &  !,psi_p2s_oa( tpsi_oa(iv_m),perv_oa(1,per_t2p_oa(lt_m)),fb_oa,fc_oa ) 
                               ,w_p                    &  !,real(kpt_m)* dti                          
                               ,dwi_p                  &  !,dti*resv_oa(per_t2p_oa(lt_m))             
                               ,fb_p                   &  !,fb_oa 
                               ,fc_p                   &
                               ,from_cplx_fun)  !,fc_oa  ) from_cplx_fun .true./.false 

              sum0_p = sum0_p + mps_p
              sum1_p = sum1_p + mps_p * s1_p
              sum2_p = sum2_p + mps_p * s2_p

           enddo
 
           ! Normalization
           sum1_p = sum1_p / sum0_p
           sum2_p = sum2_p / sum0_p

!      else if (tpsi_p .eq. 0 .or. tpsi_p .eq. 2 .or. tpsi_p.eq.3 ) then
!      if(if_print_node)write(io_unit,*) "ERROR : atom flag tpsi_p must correspond to Morlet wavelet" 
!      stop
!     else 
!      if(if_print_node)write(io_unit,*) "ERROR : unknown atom flag tpsi_p" 
!      stop
!     endif 

      end subroutine moments_of_Morlet_Power_Spectrum

      complex function conv_signal_wth_FT_morlet_psi_oa(               & 
                 tpsi_p                                                &
                ,u_p                                                   &
                ,scale_p                                               & 
                ,dwi_p                                                 & 
                ,fb_p                                                  & 
                ,fc_p                                                  &
                ,signal_p                                              &
                ,nw                                                    )
    


      ! BLXD Tests on atom flag might be removed after tests

      integer, intent(in) :: tpsi_p

      integer, intent(in) :: nw
      complex, dimension(1:nw), intent(in) :: signal_p

      double precision, intent(in) :: dwi_p, u_p
      real, intent(in)  ::                                            &
        scale_p                                                       & 
       ,fb_p                                                          & 
       ,fc_p

      complex :: sum_p, s_p
      real*8  :: w_p
      integer :: kw

!      real ::                                                         & 
!        w0_p                                                          & 
!       ,rec_p 


!.....ondelette de morlet seulement :
!      if (tpsi_p .eq. 1 ) then      

           sum_p = 0.d0
           do kw=1,nw
           ! even nw : -nw/2 to nw/2 integral -W/2 -> W/2
               w_p = real( kw-int(nw/2) ) * dwi_p
               s_p = signal_p(kw)
               !s_p = ( w_p - morlet_psi_w0c_oa( fc_p, scale_p ) )**2

               sum_p = sum_p + s_p * FT_morlet_psi_oa(  &  !                             
                                           tpsi_p       &  ! tpsi_oa(iv_m)  
                                          ,u_p          &                          
                                          ,scale_p      &  !,psi_p2s_oa( tpsi_oa(iv_m),perv_oa(1,per_t2p_oa(lt_m)),fb_oa,fc_oa ) 
                                          ,w_p          &  !,real(kpt_m)* dti                          
                                          ,dwi_p        &  !,dti*resv_oa(per_t2p_oa(lt_m))             
                                          ,fb_p         &  !,fb_oa 
                                          ,fc_p         )  !,fc_oa

           enddo 
           sum_p = sum_p * dwi_p / ( 2.d0 * pi_oa )
           conv_signal_wth_FT_morlet_psi_oa = sum_p

!      else if (tpsi_p .eq. 0 .or. tpsi_p .eq. 2 .or. tpsi_p.eq.3 ) then
!      if(if_print_node)write(io_unit,*) "ERROR : atom flag tpsi_p must correspond to Morlet wavelet" 
!      stop
!     else 
!      if(if_print_node)write(io_unit,*) "ERROR : unknown atom flag tpsi_p" 
!      stop
!     endif 

      end function conv_signal_wth_FT_morlet_psi_oa

!------------------------------------------------------------------------------
! PROCEDURE
!------------------------------------------------------------------------------
!
!> @note
!
! DESCRIPTION: 
!
!
! DESCRIPTION: 
!
!> @brief Transforms the requested time period of the analysis into the wavelet scale (s).
!
!> @details transformation periode --> echelle (s) pour l'ondelette choisie.
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - doxygen comments,  cleaning
!!  - intent in/out specification, changing if statements since tpsi_p has a single value
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!> @date 2021
!> @todo 
!
!------------------------------------------------------------------------------

      real function psi_p2s_oa(                                       &
                 tpsi_p                                               &
                ,per_p                                                &
                ,fb_p                                                 &
                ,fc_p                )

      use module_oa_variables
      implicit none

      integer, intent(in) :: tpsi_p

      real, intent(in)    ::                                          &
       per_p                                                          &
      ,fb_p                                                           &
      ,fc_p

!.....dirac:

      if (tpsi_p .eq. 0 ) then
       psi_p2s_oa = 0.

!.....ondelette de morlet:

      else if (tpsi_p .eq. 1 ) then
       psi_p2s_oa = per_p / ( 4. * pi_oa ) * (fc_p*2.*pi_oa + sqrt(2.+(fc_p*2.*pi_oa)**2))  
                                                 ! t/torrence central frequency

!....."windowed fourier":

      else if (tpsi_p .eq. 2.or.tpsi_p.eq.3 ) then
       psi_p2s_oa = per_p
      endif 

      return

      end function psi_p2s_oa

!------------------------------------------------------------------------------
! PROCEDURE
!------------------------------------------------------------------------------
!
!> @note Calculus must be adapted according to the wavelet shape (Morlet,...)
!
! DESCRIPTION: 
!
!> @brief
!> @details calcul de la taille des boites d'heisenberg et sauvegarde
!! attention reprendre calcul selon ondelette
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!> @date 2006
!> @todo BLXD
!! - exponential overflow underflow !
!------------------------------------------------------------------------------

      subroutine box_oa

      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      !use module_oa_stock

      implicit none

      integer ib

      double precision tb(20001),wb(10001),w0b,npsi2,dtb,dwb
      double precision w0bb,dtbb,dwbb,t0b 

      w0b=fc_oa*2*pi_oa
      dtb=0.1
      dwb=0.1


      tb(1)=-1000.
      do ib=1,20000
       tb(ib+1)=tb(ib)+dtb
      enddo

      wb(1)=0.
      do ib=1,10000
       wb(ib+1)=wb(ib)+dwb
      enddo

      npsi2=0.
      do ib=1,20001
       npsi2=npsi2+exp(-tb(ib)*tb(ib)/fb_oa*2)
      enddo
      npsi2=npsi2*dtb/(pi_oa**0.5) !torrence

!---------------------------------------------------------------
!      npsi2=npsi2*dtb/(2.*pi_oa)  !matlab
!---------------------------------------------------------------

      t0b=0.
!---------------------------------------------------------------
! parametres (ivane)
!
!      do ib=1,20001
!      t0b=t0b+tb(ib)*exp(-tb(ib)*tb(ib)/fb_oa*2.)
!      enddo
!      t0b=sqrt(t0b*dtb/npsi2/(2.*pi_oa)) !matlab
!!      t0b=sqrt(t0b*dtb/npsi2/(pi_oa**0.5)) !torrence
!---------------------------------------------------------------

      w0bb=0.
      do ib=1,10001
      w0bb=w0bb+wb(ib)*exp(-(wb(ib)-w0b)*(wb(ib)-w0b)*fb_oa/2.)*fb_oa/2. ! torrence
      enddo
      w0bb=w0bb*dwb/npsi2/(pi_oa**0.5) !torrence

!---------------------------------------------------------------
!      w0bb=w0bb*dwb/npsi2/(2.*pi_oa) !matlab
!---------------------------------------------------------------

      dtbb=0.
      do ib=1,20001
      dtbb=dtbb+(tb(ib)-t0b)*(tb(ib)-t0b)*exp(-tb(ib)*tb(ib)/fb_oa*2.)
      enddo
      dtbb=sqrt(dtbb*dtb/npsi2/(pi_oa**0.5)) !torrence

!---------------------------------------------------------------
!      dtbb=sqrt(dtbb*dtb/npsi2/(2.*pi_oa)) !matlab
!---------------------------------------------------------------

      dwbb=0.
      do ib=1,10001
      dwbb=dwbb+(wb(ib)-w0bb)*(wb(ib)-w0bb)*exp(-(wb(ib)-w0b)*(wb(ib)-w0b)*fb_oa/2.)*fb_oa/2.   ! torrence
      enddo

      dwbb=sqrt(dwbb*dwb/npsi2/(pi_oa**0.5)) !torrence

!---------------------------------------------------------------
!      dwbb=sqrt(dwbb*dwb/npsi2/(2.*pi_oa)) !matlab
!---------------------------------------------------------------
 
      return
      end subroutine box_oa      

!----------------------------------------------------------------------
! PROCEDURE
!
!> @note WARNING user options hardcoded in this routine 
!
! DESCRIPTION: fonctions tests pour le developpement 
!
!> @brief test variable is an homogeneous LC of cosine functions
!
!> @details  the test variable is an homogeneous 3D field 
!! varying as a LC of cosine functions with periods and amplitude
!! set in the OA namelist 
!!
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - 2020-2021 Croco version. 
!!  - 2014 namelists, doxygen comments, Stand-alone version, optimization
!> @date 2021
!> @todo BLXD CLEAN, REORGANIZE
!----------------------------------------------------------------------

      subroutine test_oa_back( tile                     & 
       ,iic_oa )
!      ,dti                        & 
!      ,imin, imax                 &
!      ,jmin, jmax                 &
!      ,kmin, kmax )

      !use module_oa_variables, only : vardp_test_oa
      use module_tile_oa, only : sts, tile_test_str

      implicit none

      ! BLXD TODO clean
      !!> Time integration step
      !!double precision, intent(in) :: dti

      !> Croco tile(-thread)
      integer, intent(in) :: tile 

      !> Current model integration iteration
      integer, intent(in) :: iic_oa 

      type(tile_test_str), pointer :: psts => null()

      double precision  :: time_t, var_t
      integer           :: i, j, k, iper, iper1, iper2

        psts => sts(tile)

        iper1=1
        iper2=nper_test

#ifdef IGW
       if ( scalogram_analysis ) then
           !iper1 = 1; iper2=nper_test
           if (iic_oa <= 52 ) then
            iper1=1 ; iper2=1
           else if ( (iic_oa > 52 ) .and. (iic_oa <= 96) ) then
            iper1=1
            iper2=nper_test
           else if (iic_oa > 96 ) then
            iper1 = 2; iper2=2
           !else
           endif
       endif
#endif

       time_t = dti * real(iic_oa)
     
       var_t=0.D0 
       do iper = iper1,iper2
         var_t = var_t + amp_test_oa(iper) * cos(2.d0*pi_oa/period_test_oa(iper)*time_t)
       end do

       ! BLXD removing temporary hardcoded option 
       !if (nper_test /= -1) then 
       do k=psts%kmin,psts%kmax
         do j=psts%jmin,psts%jmax
           do i=psts%imin,psts%imax !0,imax+1
            psts%vardp_test(i,j,k) = var_t
           enddo
         enddo
       enddo
       !end if

       psts => null() 

      return

      end subroutine test_oa_back

      subroutine test_oa( tile     & 
       ,iic_oa )
!       ichoix                     & 
!       iic_oa                     & 
!      ,dti                        & 
!      ,imin, imax                 &
!      ,jmin, jmax                 &
!      ,kmin, kmax )


      use module_tile_oa, only : sts, tile_test_str

      implicit none

      ! BLXD TODO clean
      !!> Time integration step
      !!double precision, intent(in) :: dti

      !> Croco tile(-thread)
      integer, intent(in) :: tile 

      !> Current model integration iteration
      integer, intent(in) :: iic_oa 

      type(tile_test_str), pointer :: psts => null()

      double precision  :: time_t, var_t
      integer           :: i, j, k, iper, iper1, iper2
      !logical           :: ll_ptra
#ifdef IGW
      real, parameter:: it_phy1=90.
      real              :: it_tr12, it_phy2, it_tr23, &
                           it_phy3, it_tr34, it_phy4, it_phy5, it_phy6
      real              :: delta_it, delta_amp

#endif

        psts => sts(tile)

        iper1=1
        iper2=nper_test

#ifdef IGW
      
       if_sclg : if ( scalogram_analysis ) then

           nper3 : if (nper_test==2) then

           !iper1 = 1; iper2=nper_test
           if (iic_oa <= 52 ) then
            iper1=1 ; iper2=1
           else if ( (iic_oa > 52 ) .and. (iic_oa <= 96) ) then
            iper1=1
            iper2=nper_test
           else if (iic_oa > 96 ) then
            iper1 = 2; iper2=2
           !else
           endif

           else if (nper_test==3) then nper3

           it_tr12=it_phy1+30.
           it_phy2=it_tr12+60.
           it_tr23=it_phy2+30.
           it_phy3=it_tr23+60.
           it_tr34=it_phy3+30.
           it_phy4=it_tr34+60.
           it_phy5=it_phy4+90.
           it_phy6=it_phy5+126.

           !iper1 = 1; iper2=nper_test
           if (iic_oa < it_phy1 ) then
            iper1=1 ; iper2=1
            amp_test2_oa(iper1) = amp_test_oa(iper1)
            do iper=iper2,nper_test
                amp_test2_oa(iper) = 0.
            enddo
           else if ( (iic_oa >= it_phy1 ) .and. (iic_oa <= it_tr12) ) then
            iper1=1; iper2=min(2,nper_test)
            delta_it  = (it_tr12 - it_phy1 )
            delta_amp = amp_test_oa(iper2)/delta_it
            amp_test2_oa(iper1) = amp_test_oa(iper1)
            amp_test2_oa(iper2) = delta_amp * (iic_oa - it_phy1)
            do iper=iper2+1,nper_test,1
                amp_test2_oa(iper) = 0.
            enddo
           else if ( (iic_oa > it_tr12 ) .and. (iic_oa <= it_phy2) ) then
            iper1=1; iper2=min(2,nper_test)
            amp_test2_oa(iper1) = amp_test_oa(iper1)
            amp_test2_oa(iper2) = amp_test_oa(iper2)
            do iper=iper2+1,nper_test,1
                amp_test2_oa(iper) = 0.
            enddo
           else if ( (iic_oa > it_phy2 ) .and. (iic_oa <= it_tr23) ) then
            iper1=1; iper2=min(2,nper_test)
            delta_it  = (it_tr23 - it_phy2 )
            delta_amp = - amp_test_oa(iper1)/delta_it
            amp_test2_oa(iper1) = amp_test_oa(iper1)+delta_amp*(iic_oa-it_phy2)
            amp_test2_oa(iper2) = amp_test_oa(iper2)
            do iper=iper2+1,nper_test,1
                amp_test2_oa(iper) = 0.
            enddo
           else if ( (iic_oa > it_tr23 ) .and. (iic_oa <= it_phy3) ) then
            iper1=min(2,nper_test); iper2=min(2,nper_test)
            amp_test2_oa(1) = 0. 
            amp_test2_oa(2) = amp_test_oa(2)
            amp_test2_oa(3) = 0. 
           else if ( (iic_oa > it_phy3 ) .and. (iic_oa <= it_tr34) ) then
            iper1=min(2,nper_test); iper2=min(3,nper_test)
            delta_it  = (it_tr34 - it_phy3 )
            delta_amp = amp_test_oa(3)/delta_it
            amp_test2_oa(1) = 0. 
            amp_test2_oa(2) = amp_test_oa(2)
            amp_test2_oa(3) = delta_amp * (iic_oa - it_phy3)
           else if ( (iic_oa > it_tr34 ) .and. (iic_oa <= it_phy4) ) then
            iper1=min(2,nper_test); iper2=min(3,nper_test)
            amp_test2_oa(1) = 0. 
            amp_test2_oa(iper1) = amp_test_oa(iper1)
            amp_test2_oa(iper2) = amp_test_oa(iper2)
           else if ( (iic_oa > it_phy4 ) .and. (iic_oa <= it_phy5) ) then
            iper1=min(2,nper_test); iper2=min(2,nper_test)
            amp_test2_oa(1) = 0. 
            amp_test2_oa(2) = amp_test_oa(2)
            amp_test2_oa(3) = 0. 
           else if (iic_oa > it_phy5 ) then
            iper1=1 ; iper2=1
            amp_test2_oa(iper1) = amp_test_oa(iper1)
            do iper=iper2+1,nper_test,1
                amp_test2_oa(iper) = 0.
            enddo
           endif
        endif nper3

        if (nper_test/=3) then
            do iper=iper1,iper2
                amp_test2_oa(iper) = amp_test_oa(iper)
            enddo
        endif

       endif if_sclg
#endif

       time_t = dti * real(iic_oa)
     
       var_t=0.D0 
       do iper = iper1,iper2
         var_t = var_t + amp_test2_oa(iper) * cos(2.d0*pi_oa/period_test_oa(iper)*time_t)
       end do

       ! BLXD removing temporary hardcoded option 
       !if (nper_test /= -1) then 
       do k=psts%kmin,psts%kmax
         do j=psts%jmin,psts%jmax
           do i=psts%imin,psts%imax !0,imax+1
            psts%vardp_test(i,j,k) = var_t
           enddo
         enddo
       enddo
       !end if

       psts => null()

      return

      end subroutine test_oa

!----------------------------------------------------------------------
! PROCEDURE
!
!> @note
!
!> @brief Updates the isopycne localization (closest level and depth).
!
!> @details lev_upd_oa calls update_level_oa which searches for current
!! depth position of the target density (configuration code 20 requested in the OA namelist)
!! and stores the closest level and depth in the wlev_oa structured type array.
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - 2020-2021 Croco version : croco inteface + scalogram + MPI + XIOS
!!     toward a tile(-thread) compliant version  
!!  - 2014 namelists, doxygen comments, Stand-alone version, optimization
!> @date 2021
!> @todo BLXD Croco tile-threads
!!       - wlev_oa -> st(tile)%wlev_oa (private tile dimension for field %rhp;%z,..)
!!       - pst must be a local subroutine variable pst instead of a global module var. 
!------------------------------------------------------------------------------

      subroutine  lev_upd_oa( tile                  & 
      ,rhp_t                                        &
      ,depth_t                                      & !(-1:imax+2,-1:jmax+2,0:kmax+1) 
      ,rhp_t_lbound                                 &
      ,rhp_t_ubound  )

      use module_tile_oa, only : tile_space_str, st
      use module_oa_space
      use module_oa_level

      implicit none

      !> Tile index
      integer, intent(in)  :: tile
 
      !> Lower and upper index bounds of the 3-Dimensional density array
      integer, dimension(3), intent(in) :: rhp_t_lbound, rhp_t_ubound

      !> 3D depth and density array at t-grid point
      double precision,                                                                                                  &
       dimension(rhp_t_lbound(1):rhp_t_ubound(1),rhp_t_lbound(2):rhp_t_ubound(2),rhp_t_lbound(3):rhp_t_ubound(3)), &
        intent(in) :: rhp_t, depth_t

      type(tile_space_str), pointer :: pst => null()

      ! Local indices 
      integer ::                                                      &
       i_l                                                            &
      ,j_l                                                            &
      ,k_l                                                            &
      ,ls_l                                                           &
      ,ls1_l                                                          &
      ,iv_l                                                           &
      ,il_l

      integer, pointer :: ps1_l => null()

      pst => st(tile)

      do il_l = 1 , nzlevel_oa
         iv_l = lev2v_oa(il_l)
           do ls_l = pst%begvs_oa(iv_l),pst%begvs_oa(iv_l+1)-1
            i_l = pst%l2i_oa(ls_l)
            j_l = pst%l2j_oa(ls_l)
            ps1_l=> pst%begvs3d_oa(pst%begvs_oa(iv_l))
            do ls1_l = pst%begvs3d_oa(ls_l),pst%begvs3d_oa(ls_l+1)-1
              k_l = pst%kmin3d_oa(ls_l) + (ls1_l - pst%begvs3d_oa(ls_l))* dk_oa(iv_l)
              call update_level_oa(                                   &
                i_l                                                   &
               ,j_l                                                   &
               ,pst%wlev_oa(il_l)%rhp(ls1_l-ps1_l+1)                  &
               ,pst%wlev_oa(il_l)%k  (ls1_l-ps1_l+1)                  &
               ,pst%wlev_oa(il_l)%z  (ls1_l-ps1_l+1)                  &
               ,pst%kmin, pst%kmax                                    &
               ,rhp_t                                                 &
               ,depth_t                                               &
               ,rhp_t_lbound                                          &
               ,rhp_t_ubound  )
            enddo
           enddo

      enddo
      ps1_l => null()
      pst => null()

      return
      end subroutine  lev_upd_oa

!------------------------------------------------------------------------------
! PROCEDURE
!------------------------------------------------------------------------------
!
!> @note
!
!
! DESCRIPTION: 
!
!> @brief
!
!> @details 
!
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - doxygen comments, cleaning
!!  - update_level_oa :
!!      - lower/upper k levels set to kmin,kmax instead of 1,kmax+1.
!!      - when the search fails, extrapolation treatment is removed (requires supplementary arguments).
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!> @date 2021
!> @todo update_level_oa
!! - extrapolation treatment?
!------------------------------------------------------------------------------

      subroutine update_level_oa (                  &
       i_c                                          &
      ,j_c                                          &
      ,rh_c                                         &
      ,k_c                                          &
      ,z_c                                          &
      ,kmin, kmax                                   &
      ,rhp_t                                        &
      ,depth_t                                      & !(-1:imax+2,-1:jmax+2,0:kmax+1) 
      ,rhp_t_lbound                                 &
      ,rhp_t_ubound  )


      implicit none

      !> Vertical grid index range for the analysis of density field passed in argument to the OA module 
      integer, intent(in) :: kmin, kmax

      !> Lower and upper index bounds of the 3-Dimensional density array
      integer, dimension(3), intent(in) :: rhp_t_lbound, rhp_t_ubound

      ! BLXD TODO optimize successive precision conversion from module_oa_type where rhp and Z are real precision, here real(8)
      !> 3-Dimensional density array at t-grid point
      double precision, dimension(rhp_t_lbound(1):rhp_t_ubound(1),rhp_t_lbound(2):rhp_t_ubound(2),rhp_t_lbound(3):rhp_t_ubound(3)) &
                      , intent(in) :: rhp_t, depth_t


      integer                                                         &
       i_c                                                            &
      ,j_c                                                            &
      ,k_c                                                            &
      ,iflag

      real                                                            &
       z_c                                                            &
      ,rh_c                                                           &
      ,z2_c

      integer                                                         &
       l_c                                                            &
      ,dk1_c                                                          &
      ,dk2_c                                                          &
      ,nk_c                                                           &
      ,k2_c


!.....diverses initialisations:
      nk_c = kmax
      !STDALONE k_c   = min(nk_c,max(1,k_c)) 
      k_c   = min(nk_c,max(kmin,k_c)) 
      dk1_c = 0
      dk2_c = 0
      l_c   = 0
      iflag = 0
      k2_c  = 0

!.....boucle sur l'indice l_c:
!                              l_c impair: on cherche au dessus,
!                              l_c pair  : on cherche en dessous.
!
!
!
! la sortie de la boucle est le point delicat. elle est actuellement
! basee sur 3 tests:
!
!
! test 1: la boucle est executee au moins une fois,
! test 2: 1ere partie: la recherche a ete concluante (une profondeur a ete calculee)
!         2eme partie: on sort des bornes au dessus et au dessous (iflag=2)... abandon
! test 3: la premiere recherche "en dessous" n'a rien donne et sort des bornes (<1), on cherche au dessus...
!

      do while (                                                      &
        ( dk1_c.eq.0 .and. dk2_c.eq.0 ).or.                           &
!STDALONE        (  (rh_c-rhp_t(i_c,j_c,min(nk_c,max(1,k_c+dk1_c))))*(rh_c-rhp_t(i_c,j_c,min(nk_c,max(1,k_c+dk2_c)))).gt.0   &
        (  (rh_c-rhp_t(i_c,j_c,min(nk_c,max(kmin,k_c+dk1_c))))*(rh_c-rhp_t(i_c,j_c,min(nk_c,max(kmin,k_c+dk2_c)))).gt.0   &
         .and.iflag.ne.2                                              &
         .and.k_c+dk1_c.ge.1                                          &
         .and.k_c+dk2_c.ge.1                                          &
         .and.k_c+dk1_c.le.nk_c                                       &
         .and.k_c+dk1_c.le.nk_c )                                     &
         .or.                                                         &
         k_c+dk1_c.lt.1                                               &
         )

!......incrementation de l_c l'indice de boucle:
       l_c = l_c + 1

!......calcul des increments verticaux rã©els:
       if (mod(l_c,2).ne.0) then
          dk1_c = -l_c/2-1
          dk2_c = -l_c/2
       else
          dk1_c = l_c/2-1
          dk2_c = l_c/2
       endif

!......bornes verticales:
!STDALONE      if (k_c+dk1_c.ge.1.and.k_c+dk2_c.le.nk_c) then
     if (k_c+dk1_c.ge.kmin.and.k_c+dk2_c.le.nk_c) then

!.........remise a zero de iflag (utilise pour les sorties de domaines):
          iflag = 0

!........calcul de la profondeur dans le cas ou la densite n'est pas constante.
     if (rhp_t(i_c,j_c,k_c+dk2_c)-rhp_t(i_c,j_c,k_c+dk1_c).ne.0.) then
          z2_c = depth_t(i_c,j_c,k_c+dk1_c) +                        &
          ( depth_t(i_c,j_c,k_c+dk2_c)-depth_t(i_c,j_c,k_c+dk1_c) )  &
          * ( rh_c                    -rhp_t(i_c,j_c,k_c+dk1_c) )    &  
          / ( rhp_t(i_c,j_c,k_c+dk2_c)-rhp_t(i_c,j_c,k_c+dk1_c) )     
     else
      z2_c=0.
         endif
!........calcul du niveau le plus proche (utilise le cas echeant pour demarrer
!        la recherche au prochain pas de temps).
         if ( abs(rh_c-rhp_t(i_c,j_c,k_c+dk1_c)).le. abs(rh_c-rhp_t(i_c,j_c,k_c+dk2_c))   ) then
            k2_c = k_c + dk1_c
         else 
                k2_c = k_c + dk2_c
         endif
       else
!.......dans le cas ou l'on sort du domaine (par le haut ou par le bas), 
!       on increment iflag. lorsque iflag=2, i.e. l'on sort successivement
!       par le haut et par le bas, on arrete la recherche... pas de solution! 
        iflag = iflag + 1
       endif
!.....sortie de la boucle principale:
      enddo 
!.....si l'on a quitte la boucle precedente normalement (i.e. sans etre sorti
!     du domaine...) on valide la solution z (et l'on conserve le k le plus proche
!     pour gagner du temps lors de la prochaine recherche). 
!STDALONE      if (k_c+dk1_c.ge.1.and.k_c+dk2_c.le.nk_c) then
      if (k_c+dk1_c.ge.kmin.and.k_c+dk2_c.le.nk_c) then
        z_c = z2_c
        k_c = k2_c
      else
!.......extrapolation eventuelle: attention il s'agit d'une solution tres discutable...
!       mais qui permet de fournir une solution tres approximative.
        !STDALONE
        !TODO extrapolation?

        !if (rh_c.le.rhp_t(i_c,j_c,nk_c).and.rhp_t(i_c,j_c,nk_c-1) -rhp_t(i_c,j_c,nk_c).ne.0.) then
        ! z_c = min((hssh_w(i_c,j_c,1)-h_w(i_c,j_c))*0. ,                            &
        !   depth_t(i_c,j_c,nk_c) +                                  &
        ! ( depth_t(i_c,j_c,nk_c-1)-depth_t(i_c,j_c,nk_c) )         &
        ! * ( rh_c                  -rhp_t(i_c,j_c,nk_c) )            & 
        ! / ( rhp_t(i_c,j_c,nk_c-1) -rhp_t(i_c,j_c,nk_c) )            &
        !          )
        ! k_c = nk_c
        !elseif (rhp_t(i_c,j_c,1) -rhp_t(i_c,j_c,2).ne.0.) then
        ! z_c = max(-h_w(i_c,j_c),                                    &
        !   depth_t(i_c,j_c,2) +                                     &
        ! ( depth_t(i_c,j_c,1)-depth_t(i_c,j_c,2) )                 &
        ! * ( rh_c                  -rhp_t(i_c,j_c,2) )               &
        ! / ( rhp_t(i_c,j_c,1) -rhp_t(i_c,j_c,2) )                    &
        !          )
        ! k_c = 1
        !else
          z_c = 0.
          k_c = 0
        !endif
      endif
            
      return
      end subroutine update_level_oa

!------------------------------------------------------------------------------
! PROCEDURE
!------------------------------------------------------------------------------
!
!> @note
!
! DESCRIPTION: 
!
!> @brief Handling OA module history output.
!
!> @details 
!
!
! REVISION HISTORY:
!!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Symphonie/NHOMS initial version 2006
!! - B. Lemieux-Dudon
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!!  - Croco tile(-threads) compliant version (2021)
!!  - Heisenberg theoritical uncertainties + Scalogram (june 2021)
!> @date 2021
!> @todo Croco tile(-threads)
!!       - pst must be a local subroutine variable pst instead of a global module var. 
!------------------------------------------------------------------------------

      subroutine history_oa ( ichoix, i1_h, i2_h, i3_h, i4_h, tile &
                              ,iic_oa, nt_max )

      !BLXD_TILE_ISSUE
      !use module_tile_oa, only : pst, st
      use module_tile_oa, only : st, tile_space_str
      use module_oa_variables
      use module_oa_time
      use module_oa_space
      use module_oa_periode
      !use module_oa_stock
#ifdef MPI
      use module_parameter_oa, only : iminmpi, jminmpi
#endif

      implicit none
 
      integer, intent(in) :: tile
      integer, intent(in) :: ichoix
      integer, intent(in) :: i1_h, i2_h, i3_h, i4_h

      !> Current model integration iteration
      integer, intent(in), optional :: iic_oa 

      !> Last simulation iteration index
      integer, intent(in), optional :: nt_max                                        

      type(tile_space_str), pointer :: pst => null()

      real*8  :: dt0_theo, dw0_theo, w0_theo
      real    :: dp, dp_r, dw_r_num, dw_r_num_bis, dw_r_theo
     
      integer ::                                                     &
          k_m                                                        &
         ,ir_o                                                       &
         ,l_a                                                        &
         ,lv_a                                                       &
         ,lc_a                                                       &
         ,lti_a                                                      &
         ,ic_o                                                       &
         ,ic_u                                                       &
         ,lt_u                                                       &
         ,lt_o                                                       &
         ,lp_o                                                       &
         ,iv_o                                                       &
         ,la_s                                                       &
         ,lt_a                                                       &
         ,dim_a                                                      &
         ,iv_u                                                       &
         ,iv_s

       integer            :: lp_o1
       character(len=250) :: file_hist

! BLXD TODO tile-thread ISSUE TO SOLVE
#ifdef MPI
      if (.not. if_print_node) return
#endif

!---->Fichier de sortie:

      file_hist = trim(directory_out_oa) // txtslash // 'history_oa.dat'

      if (ichoix.eq.1) then
!*******************************************************************************
! allocate_wlev_oa_ptr
!*******************************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
       !if (ifl_oa_out.eq.1) call sequential_begin()
       open(unit=io_hist,file=trim(file_hist),position='append')   
       write (io_hist,*) 'allocation: level n°',i1_h,',var n°',i2_h,',dim=',i3_h
       close(io_hist)
       !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif

      if (ichoix.eq.2) then
!*******************************************************************************
! allocate_win_oa  
!*******************************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
       !if (ifl_oa_out.eq.1) call sequential_begin()
       open(unit=io_hist,file=trim(file_hist),position='append')   
!       write (io_hist,*) 'allocation: fenetre n°',l_a,',iic_oa=',iic_oa       &
!                    ,',lc=',i1_h                                       &
!                    ,',lv=',i2_h                                       &
!                    ,',dim=',i3_h                                      &
!                    ,',lti=',i4_h
       write (io_hist,*) 'allocation: fenetre n°',i3_h,',iic_oa=',iic_oa       &
                    ,',lc=',i1_h                                       &
                    ,',lv=',i2_h                                       &
                    ,',dim=?'                                      &
                    ,',lti=',i4_h
       close(io_hist)
       !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif

      if (ichoix.eq.3) then
!*******************************************************************************
! deallocate_win_oa  
!*******************************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()
        open(unit=io_hist,file=trim(file_hist),position='append')   
        write (io_hist,*) 'desallocation: fenetre n°',tallocated_oa(i1_h)
        close(io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif


      if (ichoix.eq.4) then
!*******************************************************************************
! initial_oa
!*******************************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()
        open(unit=io_hist,file=trim(file_hist),position='append')   
        write(io_hist,*)
        write(io_hist,*) '*************************'
        write(io_hist,*) 'subroutine: initial_oa'
        write(io_hist,*) '*************************'
        close(io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif

      chx5 : if (ichoix.eq.5) then
!*******************************************************************************
! initial_oa   
!*******************************************************************************
!******************************************************************
!        caracteristiques d'une variable:      
!******************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()

         pst => st(tile)

         open(unit=io_hist,file=trim(file_hist),position='append')   
         loop_config : do ic_o=1,nzc_oa
            write(io_hist,*)
            write(io_hist,*)
            write(io_hist,*) '*************************************'
            write(io_hist,*) ' configuration de variables n°',ic_o
            write(io_hist,*) '*************************************'
            if (tc_oa(ic_o).ge.100) then
               write (io_hist,*) 'type de configuration:',tc_oa(ic_o)
               write (io_hist,*) 'WARNING ! composite configuration are hardcoded'
               write (io_hist,*) 'Composite configuration never tested with croco'
               write (io_hist,*) 'ERROR : check/change the configuration code, remove stop and recompile'
               stop 
            else
               write (io_hist,*) 'type de cfg-variable simple:',tc_oa(ic_o)
            endif

            loop_variable : do iv_o=begc_oa(ic_o),begc_oa(ic_o+1)-1
               write(io_hist,*)
               write(io_hist,*) '*************************************'
               write(io_hist,*) '        variable n°',iv_o
               write(io_hist,*) '*************************************'

               write(io_hist,*) '----------------------'
               write(io_hist,*) '---caracteristiques---'
               write(io_hist,*) '----------------------'
               write(io_hist,*) 'type            : ',tv_oa(iv_o) 
               if (swt_wfpf_oa(iv_o).eq.1) then
                  write(io_hist,*) 'calcul du coef. reel'
               endif 
               if (swt_wfpf_oa(iv_o).eq.2) then
                  write(io_hist,*) 'calcul du module au carre'
               endif 
               if (swt_wfpf_oa(iv_o).eq.3) then
                  write(io_hist,*) 'calcul du module'
               endif 
               if (swt_wfpf_oa(iv_o).eq.4) then
                  write(io_hist,*) 'calcul du coef. complexe'
               endif 
               if (save_oa(iv_o).eq.0) then
                  write(io_hist,*) 'variable non sauvegardee dans un fichier'
               endif 
               if (save_oa(iv_o).eq.1) then
                  write(io_hist,*) 'variable sauvegardee dans un fichier'
               endif 
               if (updv_oa(iv_o).eq.1) then
                  write(io_hist,*) 'variable utilisee pour une mise a jour',iv_s
               endif 

               write(io_hist,*) 
               write(io_hist,*) '----------------------------------'
               write(io_hist,*) '---caracteristiques de l''atome---'
               write(io_hist,*) '----------------------------------'

               write(io_hist,*) 'type de l''atome                : ',tpsi_oa(iv_o)
               
               if (tpsi_oa(iv_o).eq.0) then
                  write(io_hist,*) 'dirac'
               !endif
               else if (tpsi_oa(iv_o).eq.1) then
                  write(io_hist,*) 'ondelette morlet complexe fb_oa,fc_oa: ',fb_oa,fc_oa
               !endif
               else if (tpsi_oa(iv_o).eq.2) then
                  write(io_hist,*) 'windows fourier'
               !endif
               else if (tpsi_oa(iv_o).eq.3) then
                  write(io_hist,*) 'transformee de fourier classique'
               endif

               if (tpsi_oa(iv_o).ne.3) then
                  write(io_hist,*) 'largeur (h ou nbre de periodes): ',delta_t_oa(iv_o)
               endif
               write(io_hist,*) 'nombre de points par periode   : ',nzpt_per_oa(iv_o)

               write(io_hist,*) 
               write(io_hist,*) '------------------------------'
               write(io_hist,*) '---configuration temporelle---'
               write(io_hist,*) '------------------------------'
               if (unite_oa.eq.1.) then
                write(io_hist,*) 'unites                     : ','secondes'
               else
               if (unite_oa.eq.3600.) then
                 write(io_hist,*) 'unites                     : ','heures'
                else
                 write(io_hist,*) 'unites                     : ','autre'
                 endif
               endif
               write(io_hist,*) 'type echantillonnage       : ',swt_t_oa(iv_o)
               write(io_hist,*) 'premiere sortie (h/s)      : ',t0_oa(iv_o)
               write(io_hist,*) 'discretisation/integration : ',dori_oa(iv_o)

               if (swt_t_oa(iv_o).eq.1.or.swt_t_oa(iv_o).eq.4) then
                  write(io_hist,*) 'extraction: kount initial  : ',kount_user_oa(1,iv_o)
                  write(io_hist,*) '            kount final    : ',kount_user_oa(2,iv_o)
                  write(io_hist,*) '            delta(kount)   : ',kount_user_oa(3,iv_o)
               else if (swt_t_oa(iv_o).eq.3) then
                  write(io_hist,*) 'extraction: kount initial: ',kount_user_oa(1,iv_o)
               !else if (swt_t_oa(iv_o).eq.3) then #BLDX croco nt_max ?
                  write(io_hist,*) 'extraction: kount final  : ',nt_max-1
               endif 

               do lt_o  = begvt_oa(iv_o) , begvt_oa(iv_o+1) - 1
                  ! #BLXD Corrected Log. 
                  ! see call var_oa if (   kountv_oa(1,lt_o).ne.-9999.or.kountv_oa(2,lt_o).ne.-9999) then
                  if (   kountv_oa(1,lt_o).ne.-9999.and.kountv_oa(2,lt_o).ne.-9999) then
                     write(io_hist,*) 'BLXD localisation de l''atome n°',lt_o-begvt_oa(iv_o)+1,' [t1,t2]= ',kountv_oa(1,lt_o), kountv_oa(2,lt_o)
                  else
                     write(io_hist,*) 'pas d''analyse pour la localisation n°',lt_o-begvt_oa(iv_o)+1
                  endif
               enddo

               write(io_hist,*) 
               write(io_hist,*) '----------------------------'
               write(io_hist,*) '---configuration spatiale---'
               write(io_hist,*) '----------------------------'
               write(io_hist,*) 'type d''echantillonnage        : ',swt_d_oa(iv_o)  

               if (swt_d_oa(iv_o).eq.1.or.swt_d_oa(iv_o).eq.3) then
#ifdef SPHERICAL
                  write(io_hist,*) ' BLXD CHECK lon, lat min max'
         if ( if_sphe_deg2rad ) then
                  write(io_hist,*) 'degree    : lat min,lat max   : ',lat_oa(1,iv_o)*180/pi_oa,lat_oa(2,iv_o)*180/pi_oa
                  write(io_hist,*) '            lon min,lon max   : ',lon_oa(1,iv_o)*180/pi_oa,lon_oa(2,iv_o)*180/pi_oa
         else
                  write(io_hist,*) 'natif     : lat min,lat max   : ',lat_oa(1,iv_o),lat_oa(2,iv_o)
                  write(io_hist,*) '            lon min,lon max   : ',lon_oa(1,iv_o),lon_oa(2,iv_o)
         endif
#else
                  write(io_hist,*) 'coord.    : lat min,lat max   : ',lat_oa(1,iv_o),lat_oa(2,iv_o)
                  write(io_hist,*) '            lon min,lon max   : ',lon_oa(1,iv_o),lon_oa(2,iv_o)
#endif
               else if (swt_d_oa(iv_o).eq.2) then
                  write(io_hist,*) 'extraction: point loc. (i,j)  : ',ptij_oa(1,iv_o),ptij_oa(2,iv_o)
#ifdef MPI
                  write(io_hist,*) 'extraction: point glob (i,j)  : ',ptij_oa(1,iv_o)+ iminmpi-1,ptij_oa(2,iv_o)+ jminmpi-1
#endif

               endif

               if (swt_d_oa(iv_o).eq.3)                               &
                  write(io_hist,*) '            prof min,prof max : ',h_oa(1,iv_o),h_oa(2,iv_o)
               write(io_hist,*)    '            kmin,kmax         : ',k_oa(1,iv_o),k_oa(2,iv_o)
               if (swt_d_oa(iv_o).ne.2) then
                  write(io_hist,*) '            di,dj,dk          : ',dx_oa(iv_o),dy_oa(iv_o),dk_oa(iv_o)
               else
                  write(io_hist,*) '            dk                : ',dk_oa(iv_o)
               endif

               write(io_hist,*) 
               write(io_hist,*) '---------------------------------'
               write(io_hist,*) '---configuration frequentielle---'
               write(io_hist,*) '---------------------------------'
               if (dori_oa(iv_o).eq.1) then
                  write(io_hist,*) 'periodes (t) discretes'
               else
                  write(io_hist,*) 'periodes (t) integrees'
               endif
               write(io_hist,*)    'tmin,dt,tmax (h ou s)      : ',per_oa (1,iv_o)/unite_oa,per_oa (3,iv_o)/unite_oa,per_oa (2,iv_o)/unite_oa
               lp_o1 = begvp_oa(iv_o)
               per_loop : do lp_o = begvp_oa(iv_o) , begvp_oa(iv_o+1)-1
                  write (io_hist,*) 
                  write (io_hist,*) '=> Periode n°',lp_o-begvp_oa(iv_o)+1,'en h ou s     :',perv_oa(1,lp_o)/unite_oa
                  write (io_hist,*) '   coef. reconstruction  :',perv_oa(2,lp_o)
                  write (io_hist,*) '   resolution de l''atome :',resv_oa(lp_o)
                  write (io_hist,*) 
                  if_hsbg : if (if_heisenberg_box) then
                  if (tpsi_oa(iv_o).eq.1) then
                        w0_theo = morlet_psi_w0c_oa(                                                       & 
                          fc_p   =fc_oa                                                                         &
                         ,scale_p=morlet_psi_p2s_oa( perv_oa(1,lp_o), fb_oa, fc_oa )                       ) 

                        dw0_theo = morlet_psi_dw0_oa(                                                      & 
                          fb_p=fb_oa                                                                         &
                         ,fc_p=fc_oa                                                                         &
                         ,scale_p=morlet_psi_p2s_oa( perv_oa(1,lp_o), fb_oa, fc_oa )                       ) 

                        dt0_theo = morlet_psi_dt0_oa(                                                      & 
                          fb_p   =fb_oa                                                                         &
                         ,fc_p   =fc_oa                                                                         &
                         ,scale_p=morlet_psi_p2s_oa( perv_oa(1,lp_o), fb_oa, fc_oa )                       )

                      write (io_hist,*) '=> Incertitudes d''Heisenberg'
                      write (io_hist,*) '   ** Psi L2-norm            :',psi_norm_l2(lp_o), psi_norm_l2_bis(lp_o)
                      write (io_hist,*) '   ** Analysed Per./freq.    :',perv_oa(1,lp_o), (2.*pi_oa/perv_oa(1,lp_o))
                      write (io_hist,*) '   ** Temporal box estimate'
                      write (io_hist,*) '      time box center        :',t0(lp_o)         , t0_bis(lp_o)
                      write (io_hist,*) '      time box width         :',dt0(lp_o)        , dt0_bis(lp_o)
                      write (io_hist,*) '   ** Theoritical temporal box'
                      !write (io_hist,*) '      time box center        :',t0_theo(lp_o)
                      !write (io_hist,*) '      time box width         :',dt0_theo(lp_o)
                      write (io_hist,*) '      time box center        :',t0_theo(lp_o)
                      write (io_hist,*) '      time box width         :',dt0_theo
                      write (io_hist,*) '   ** Frequential box estimate'
                      write (io_hist,*) '      freq.box center        :',w0(lp_o)         , w0_apx(lp_o)
                      write (io_hist,*) '      freq.box width         :',dw0(lp_o)        , dw0_apx(lp_o)
                      write (io_hist,*) '   ** Theoritical frequential box'
                      !write (io_hist,*) '      freq.box center        :',w0_theo(lp_o) 
                      !write (io_hist,*) '      freq.box width         :',dw0_theo(lp_o)
                      write (io_hist,*) '      freq.box center        :',w0_theo
                      write (io_hist,*) '      freq.box width         :',dw0_theo
                      dp = perv_oa(1,lp_o) - perv_oa(1,lp_o1)
                      if ( dp > 0. ) then
                        dp_r = dp / ( perv_oa(1,lp_o) + perv_oa(1,lp_o1) )
                        dw_r_num      = dw0(lp_o)     / w0(lp_o)
                        dw_r_num_bis  = dw0_apx(lp_o) / w0_apx(lp_o)
                        dw_r_theo     = dw0_theo/ w0_theo
                        write (io_hist,*) '   ** User defined period progression of the scalogram  ',dp 
                        write (io_hist,*) '   ** User def. period progr. ratio vs numerical dw0/w0 ',dp_r, dw_r_num, dw_r_num_bis 
                        write (io_hist,*) '   ** User def. period progr. ratio vs numerical dw0/w0 ',dp_r, dw_r_theo
                      endif
                      lp_o1 = lp_o
                  endif
                  endif if_hsbg
               enddo per_loop

            enddo loop_variable
         enddo loop_config

         write(io_hist,*)
         write(io_hist,*)
         write(io_hist,*) '*************************************'
         write(io_hist,*) '          tailles en memoire'
         write(io_hist,*) '       des differentes structures'
         write(io_hist,*) '*************************************'
         write(io_hist,*)
! #BLXD
         write(io_hist,*) 'allocations simultanees (user/=-99)  : ',nmsimult_oa_max
         write(io_hist,*) 'allocations simultanees (set/calc.)  : ',nmsimult_oa
         write(io_hist,*) 'structure temporelle          : ',nzvt_oa
         !write(io_hist,*) 'structure spatiale (2d)       : ',nzvs_oa
         write(io_hist,*) 'BLXD CHECK tile structure spatiale (2d)       : ',pst%nzvs_oa
         !write(io_hist,*) 'structure spatiale (3d)       : ',nzvs3d_oa
         write(io_hist,*) 'BLXD CHECK tile structure spatiale (3d)       : ',pst%nzvs3d_oa
         write(io_hist,*) 'structure frequentielle       : ',nzvp_oa
         close(io_hist)
         !if (ifl_oa_out.eq.1) call sequential_end()
       !endif

         pst => null()

      endif chx5

      chx6 : if (ichoix.eq.6) then
!*******************************************************************************
! initial_oa   
!*******************************************************************************
!******************************************************************
!     alocations dynamiques:
!******************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()
      
        pst => st(tile)

        open(unit=io_hist,file=trim(file_hist),position='append')   
        write (io_hist,*)
        write(io_hist,*) '*************************'
        write (io_hist,*) 'subroutine: allocate::'
        write(io_hist,*) '*************************'
        write (io_hist,*)
! #BLXD
        write (io_hist,*) 'nmsimult_oa_max =',nmsimult_oa_max
        write (io_hist,*) 'nmsimult_oa =',nmsimult_oa

! #BLXD
        write (io_hist,*) 'nzv_oa      =',nzv_oa
        write (io_hist,*) 'BLXD nzvs_oa     =',pst%nzvs_oa
        write (io_hist,*) 'nzvt_oa     =',nzvt_oa
        write (io_hist,*) 'nzvp_oa     =',nzvp_oa
! #BLD typo 
!        write (io_hist,*) 'nzvc_oa     =',nzc_oa
         write (io_hist,*) 'nzvc_oa     =',nzvc_oa
        write (io_hist,*) 'BLXD nzvs3d_oa   =',pst%nzvs3d_oa
        write (io_hist,*)
        close(io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
        pst => null()

      endif chx6

!STDALONE struct_oa eliminated      if (ichoix.eq.7) then
!*******************************************************************************
! struct_oa   
!*******************************************************************************
!      !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
!       !if (ifl_oa_out.eq.1) call sequential_begin()
!       open(unit=io_hist,file=trim(file_hist),position='append')   
!       write(io_hist,*)
!       write(io_hist,*) '*************************'
!       write(io_hist,*) 'subroutine: struct_oa'
!       write(io_hist,*) '*************************'
!       close(io_hist)
!       !if (ifl_oa_out.eq.1) call sequential_end()
!      !endif
!     endif

!STDALONE struct_oa eliminated     if (ichoix.eq.8) then
!*******************************************************************************
! struct_oa   
!*******************************************************************************
!      !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
!       !if (ifl_oa_out.eq.1) call sequential_begin()
!       !open(unit=io_hist,file=trim(dir_history_exp)//trim(file_hist)//'_'//dom_c//'.out',position='append')    
!       open(unit=io_hist,file=trim(file_hist),position='append')   
!       write(3,*) '...fichier structure initialise.'
!       close(3)
!       !if (ifl_oa_out.eq.1) call sequential_end()
!      !endif
!     endif

!STDALONE subsave_init_oa eliminated     if (ichoix.eq.9) then
!*******************************************************************************
! subsave_init_oa   
!*******************************************************************************
!      !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
!       !if (ifl_oa_out.eq.1) call sequential_begin()
!       !open(unit=io_hist,file=trim(dir_history_exp)//trim(file_hist)//'_'//dom_c//'.out',position='append')    
!       open(unit=io_hist,file=trim(file_hist),position='append')   
!       write(io_hist,*)
!       write(io_hist,*) '*****************************'
!       write(io_hist,*) 'subroutine: save_init_oa'
!       write(io_hist,*) '*****************************'
!       close(io_hist)
!       !if (ifl_oa_out.eq.1) call sequential_end()
!      !endif
!     endif

      if (ichoix.eq.10) then
!*******************************************************************************
! subsave_oa   
!*******************************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()
        open(unit=io_hist,file=trim(file_hist),position='append')   
        write(io_hist,*) '*** sauvegarde de la fenetre n°',i1_h,'***'
        close(io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
        !endif
      endif

      if (ichoix.eq.11) then
!*******************************************************************************
! upd_init_oa   
!*******************************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()
        open(unit=io_hist,file=trim(file_hist),position='append')   
        write(io_hist,*)
        write(io_hist,*) '*************************************'
        write(io_hist,*) ' allocation dynamique des variables 2d/3d'
        write(io_hist,*) '  => nzupd0d_oa,nzupd2d_oa,nzupd3d_oa:',nzupd2d_oa,nzupd3d_oa
        write(io_hist,*) ' allocation dynamique des variables 0d-scalogram'
        write(io_hist,*) '  => nzupd0d_oa                      :',nzupd0d_oa
        write(io_hist,*) '*************************************'
        write(io_hist,*)
        close (io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif


      if (ichoix.eq.12) then
!*******************************************************************************
! var_upd_oa      
!*******************************************************************************
!STDALONE ichoix==12 NOT FOUND
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
       !if (ifl_oa_out.eq.1) call sequential_begin()
        open(unit=io_hist,file=trim(file_hist),position='append')   
        write(io_hist,*) 'initialisation des variables mises a jour'
        close(io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif

      if (ichoix.eq.13) then
!*******************************************************************************
! var_upd_oa      
!*******************************************************************************
!STDALONE ichoix==13 NOT FOUND
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()
        open(unit=io_hist,file=trim(file_hist),position='append')   
        write(io_hist,*) 'mise a jour des variables specifiques it=',i1_h,tallocated_oa(i1_h)
        close(io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif

      if (ichoix.eq.14) then
!*******************************************************************************
! var_upd_oa      
!*******************************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()
        open(unit=io_hist,file=trim(file_hist),position='append')   
        write(io_hist,*) 'mise a jour des variables specifiques it=',i1_h,tallocated_oa(i1_h)        
        write(io_hist,*) ' - tupd_oa    ',tupd_oa(i2_h)
        write(io_hist,*) ' - 2d-3d var  ',tgv3d_oa(tv_oa(i2_h))
        write(io_hist,*) ' - coeff type ',swt_wfpf_oa(i2_h)
        close(io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif

      if (ichoix.eq.15) then
!*******************************************************************************
! var_upd_oa      
!*******************************************************************************
       !if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
        !if (ifl_oa_out.eq.1) call sequential_begin()
        open(unit=io_hist,file=trim(file_hist),position='append')   
        write(io_hist,*) 'mise a jour des coord. des scalogram ip, perv =',i1_h,perv_oa(1,i1_h) 
        write(io_hist,*) ' - index global du scalogram tupd_oa    ',tupd_oa(i2_h)
        close(io_hist)
        !if (ifl_oa_out.eq.1) call sequential_end()
       !endif
      endif

!STDALONE energie_oaavelet removed      if (ichoix.eq.15) then
!*******************************************************************************
! var_upd_oa      
!*******************************************************************************
!      if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
!       if (ifl_oa_out.eq.1) call sequential_begin()
!       open(unit=io_hist,file=trim(dir_history_exp)//trim(file_hist)//'_'//dom_c//'.out',position='append')    
!       write (io_hist,*) 'appel de la routine energie_oaavelet: ','ic=',i1_h
!       close (io_hist)
!       if (ifl_oa_out.eq.1) call sequential_end()
!      endif
!      endif

!STDALONE energie_oaavelet removed      if (ichoix.eq.16) then
!*******************************************************************************
! var_upd_oa   
!*******************************************************************************
!      if (ifl_oa_out.eq.1.or.(ifl_oa_out.eq.2.and.par%rank.eq.0)) then
!       if (ifl_oa_out.eq.1) call sequential_begin()
!       open(unit=io_hist,file=trim(dir_history_exp)//trim(file_hist)//'_'//dom_c//'.out',position='append')    
!       write (io_hist,*) 'appel de la routine energie_oaavelet: ' ,'ic=',i1_h
!       close (io_hist)
!       if (ifl_oa_out.eq.1) call sequential_end()
!      endif
!      endif

       return

       end subroutine history_oa

       integer function pers_t2p_oa( iv, lt )
        use module_oa_periode, only : per_t2p_oa, begvp_oa
        implicit none
        integer, intent(in) :: iv, lt
        if ( iv == 1 ) then
            pers_t2p_oa = per_t2p_oa(lt)
        else if ( iv > 1 ) then
            pers_t2p_oa = per_t2p_oa(lt) - begvp_oa(iv) + 1
        else
           print*,'ERROR OA : issue with variable neg index iv ', iv
           stop
        endif
        return
       end function 

       function set_io_nodoa(iv, mynod, inteven, intodd)
        use module_oa_variables, only : nzv_oa
        implicit none
        integer, intent(in) :: mynod, inteven, intodd, iv
        integer :: ionodoa, set_io_nodoa
        if (nzv_oa >= 100) then
           ionodoa=1000
        else if (nzv_oa >= 10) then
           ionodoa=100
        else
           ionodoa=10
        endif 
        if ( MOD(iv,2)==0 ) then
         set_io_nodoa = (inteven*ionodoa+iv)*1000+mynod
        else
         set_io_nodoa = (intodd*ionodoa+iv)*1000+mynod
        end if
        return
       end function set_io_nodoa

       subroutine test_mpi_blocking_send()
       
       IMPLICIT NONE
       
       INTEGER :: ierror
       INTEGER :: size
       INTEGER, PARAMETER :: sender_rank = 0
       INTEGER, PARAMETER :: receiver_rank = 1
       INTEGER :: my_rank
       INTEGER :: buffer
       
       ! Get my rank and do the corresponding job
       ! CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)
       CALL MPI_Comm_rank(comm, my_rank, ierror)
       SELECT CASE (my_rank)
           CASE (sender_rank)
               buffer = 12345
               WRITE(*,'(A,I0,A,I0,A)') 'MPI process ', my_rank, ' sends value ', buffer, '.'
               CALL MPI_Send(buffer, 1, MPI_INTEGER, receiver_rank, 0, comm, ierror)
           CASE (receiver_rank)
               !CALL MPI_Recv(buffer, 1, MPI_INTEGER, sender_rank, 0, comm, MPI_STATUS_IGNORE, ierror)
               CALL MPI_Recv(buffer, 1, MPI_INTEGER, MPI_ANY_SOURCE, 0, comm, MPI_STATUS_IGNORE, ierror)
               WRITE(*,'(A,I0,A,I0,A)') 'MPI process ', my_rank, ' received value ', buffer, '.'
       END SELECT 
       
       end subroutine test_mpi_blocking_send
       
       subroutine test_mpi_nonblocking_ssend(do_things_outside_case)
       
       IMPLICIT NONE
       
       INTEGER :: ierror
       INTEGER :: size
       INTEGER, PARAMETER :: sender_rank = 0
       INTEGER, PARAMETER :: receiver_rank = 1
       INTEGER :: my_rank
       INTEGER :: buffer
       LOGICAL, INTENT(IN) :: do_things_outside_case
       INTEGER :: request, x
       
       ! Get my rank and do the corresponding job
       ! CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)
       CALL MPI_Comm_rank(comm, my_rank, ierror)
       
       x = 0
       request=-99
       SELECT CASE (my_rank)
         CASE (sender_rank)
             buffer = 6789
             WRITE(*,'(A,I0,A,I0,A)') 'MPI process ', my_rank, ' sends value with ISsend ', buffer, '.'
             !CALL MPI_Issend(buffer, 1, MPI_INTEGER, receiver_rank, 0, MPI_COMM_WORLD, request, ierror)
             CALL MPI_Issend(buffer, 1, MPI_INTEGER, receiver_rank, 0, comm, request, ierror)
             if ( .not. do_things_outside_case) then    
                 ! Do other things while the MPI_Isend completes
                 ! <...>
                 WRITE(*,'(A,I0,A,I0)') 'I am MPI process ', my_rank, ' and I do things before WAIT INSIDE CASE with req =',request
                 x = x + 1 
                 WRITE(*,'(A,I0,A,I0)') 'I am MPI process ', my_rank, ' and I do things before WAIT and x = ',x
                 ! Let's wait for the MPI_Isend to complete before progressing further.
                 CALL MPI_Wait(request, MPI_STATUS_IGNORE, ierror)
             endif
         CASE (receiver_rank)
             !CALL MPI_Recv(buffer, 1, MPI_INTEGER, sender_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
             CALL MPI_Recv(buffer, 1, MPI_INTEGER, sender_rank, 0, comm, MPI_STATUS_IGNORE, ierror)
             WRITE(*,'(A,I0,A,I0,A)') 'MPI process ', my_rank, ' received value with Recv ', buffer, '.'
       END SELECT 
       
       IF (my_rank==sender_rank) then
       if ( do_things_outside_case) then    
           ! Do other things while the MPI_Isend completes
           ! <...>
           WRITE(*,'(A,I0,A,I0)') 'I am MPI process ', my_rank, ' and I do things before WAIT OUTSIDE CASE with req =',request
           x = x + 1 
           WRITE(*,'(A,I0,A,I0)') 'I am MPI process ', my_rank, ' and I do things before WAIT and x = ',x
           ! Let's wait for the MPI_Isend to complete before progressing further.
           CALL MPI_Wait(request, MPI_STATUS_IGNORE, ierror)
       endif
       ENDIF
       
       end subroutine test_mpi_nonblocking_ssend

      end module module_interface_oa
 
#else /* ONLINE_ANALYSIS */
      module module_interface_oa_empty
      end module
#endif /* ONLINE_ANALYSIS */

