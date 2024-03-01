! $Id: main.F 1615 2014-12-17 13:27:07Z rblod $
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
!
      program main
!
!======================================================================
!
!                     OCEAN MODEL MAIN DRIVER
!
!    Advances forward the equations for all nested grids, if any.
!
!======================================================================
!
#if defined OA_COUPLING || defined OW_COUPLING
          USE mod_prism
#endif
#ifdef XIOS
          USE xios           ! XIOS module
#endif
#ifdef PISCES
          USE pisces_ini     ! PISCES modules
          USE trcini_pisces
#endif
#ifdef SUBSTANCE
          USE substance, ONLY : substance_read_alloc
          USE substance, ONLY : substance_surfcell
#endif
#ifdef MUSTANG
          USE plug_MUSTANG_CROCO, ONLY : mustang_init_main
#endif
#ifdef ONLINE_ANALYSIS
          USE module_interface_oa, only : init_parameter_oa
     &                                   ,if_oa
#endif
!     implicit none
      integer tile, subs, trd, ierr
#include "param.h"
#include "private_scratch.h"
#include "nbq.h"
#include "scalars.h"
#include "ncscrum.h"
#include "grid.h"
#include "ocean2d.h"
#include "mpi_cpl.h"
#ifdef FLOATS
# include "floats.h"
# include "ncscrum_floats.h"
#endif
#ifdef STATIONS
# include "sta.h"
# include "nc_sta.h"
#endif
#ifdef AGRIF
      Type(Agrif_pgrid),pointer  :: parcours
#endif
#ifdef MPI
      include 'mpif.h'
!     real*8 start_time2, start_time1, exe_time
#endif
      integer :: iifroot, iicroot
#ifdef WKB_WWAVE
      integer winterp
#endif
#ifdef AGRIF
      integer size_XI,size_ETA,se,sse, sz,ssz
      external :: step
# include "zoom.h"
# include "dynparam.h"
#endif
#ifdef WKB_WWAVE
# include "wkb_wwave.h"
#endif
#ifdef USE_CALENDAR
      character*19 :: tool_sectodat
!      real*8       :: tool_datosec
#endif
#ifdef AUTOTILING
! Initial values of NSUB_X, NSUB_E
      NSUB_X = 1
      NSUB_E = NPP
#endif
!
#include "dynderivparam.h"
!
!----------------------------------------------------------------------
!  Initialize communicators and subgrids decomposition:
!  MPI parallelization, XIOS server, OASIS coupling, AGRIF nesting,
!----------------------------------------------------------------------
!
#ifdef MPI
# if (!defined AGRIF && !defined OA_COUPLING && !defined OW_COUPLING)
      call MPI_Init (ierr)
# endif
!
!  XIOS, OASIS and AGRIF: split MPI communicator
!  (XIOS with OASIS is not done yet)
!
# if (defined XIOS && !defined AGRIF)
      call xios_initialize( "crocox",return_comm=MPI_COMM_WORLD )
# elif (defined OA_COUPLING && !defined AGRIF)
      call cpl_prism_init  ! If AGRIF --> call cpl_prism_init in zoom.F
# elif (defined OW_COUPLING && !defined AGRIF)
      call cpl_prism_init  ! In AGRIF case, cpl_prism_init is in zoom.F
# elif defined AGRIF
      call Agrif_MPI_Init(MPI_COMM_WORLD)
# endif
#endif /* MPI */
!
!  Initialize AGRIF nesting
!
#ifdef AGRIF
      call Agrif_Init_Grids()
      call declare_zoom_variables()
#endif
!
!  Initialize automatic tiling
!
#ifdef AUTOTILING
      call init_auto_tiling
#endif
!
!  Setup MPI domain decomposition
!
#ifdef MPI
!     start_time1=PMPI_Wtime()
      call MPI_Setup (ierr)
      if (ierr.ne.0) goto 100        !--> ERROR
#endif
!
!  Initialize debug procedure
!
#if defined RVTK_DEBUG || defined RVTK_DEBUG_ADVANCED || \
    defined RVTK_DEBUG_PERFRST

      call debug_ini
#endif
!
#define CR  !
!
!----------------------------------------------------------------------
!  Read in tunable model parameters in roms.in file
!----------------------------------------------------------------------
!
      call read_inp (ierr)
      if (ierr.ne.0) goto 100
!
!----------------------------------------------------------------------
!  Initialize global model parameters
!----------------------------------------------------------------------
!
!  Gobal scalar variables
!
      call init_scalars (ierr)
      if (ierr.ne.0) goto 100
!
#ifdef SOLVE3D
!
!  PISCES biogeochemeical model parameters
!
# if defined BIOLOGY && defined PISCES
      call trc_nam_pisces
# endif
!
!
!  Read sediment initial values and parameters from sediment.in file
!
# ifdef SEDIMENT
#  ifdef AGRIF
      if (Agrif_lev_sedim.EQ.0) call init_sediment
#  else
      call init_sediment
#  endif
# endif
#endif

#ifdef SUBSTANCE
!
!  Substance var need for MUSTANG and BIOLink
!
      call substance_read_alloc(may_day_flag,indxT,indxTsrc)
#endif

!
! Online spectral analysis module
!
# ifdef ONLINE_ANALYSIS
      CALL init_parameter_oa(
     &  io_unit_oa=stdout
#ifdef MPI
     & ,if_print_node_oa=(mynode==0)
#else
     & ,if_print_node_oa=.true.
#endif
#ifdef MPI
     & ,mynode_oa=mynode
     & ,comm_oa=MPI_COMM_WORLD
#else
     & ,mynode_oa=0
     & ,comm_oa=0
#endif
     & ,dti_oa=dt
     & ,kount0_oa=ntstart-1
     & ,nt_max_oa=ntimes
     & ,ntiles=NSUB_X*NSUB_E)
# endif
!
!
!----------------------------------------------------------------------
!  Create parallel threads; start timers for each thread;
!  initialize (FIRST-TOUCH) model global arrays (most of them
!  are just set to to zero).
!----------------------------------------------------------------------
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call start_timers()
        call init_arrays (tile)
      enddo
CR      write(*,*) '-11' MYID
!
!----------------------------------------------------------------------
!  Set horizontal grid, model bathymetry and Land/Sea mask
!----------------------------------------------------------------------
!
#ifdef ANA_GRID
!
!  Set grid analytically
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call ana_grid (tile)
      enddo
# if defined RVTK_DEBUG || defined RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
       call check_tab2d(h(:,:),'h initialisation #1','r')
C$OMP END MASTER
# endif
#else
!
!  Read grid from GRID NetCDF file
!
      call get_grid
      if (may_day_flag.ne.0) goto 99 !-->  EXIT
#endif
!
!  Compute various metric term combinations.
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call setup_grid1 (tile)
      enddo
CR      write(*,*) '-10' MYID
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call setup_grid2 (tile)
      enddo
CR      write(*,*) ' -9' MYID
!
!----------------------------------------------------------------------
!  Setup vertical grid variables setup vertical S-coordinates
!  and fast-time averaging for coupling of
!  split-explicit baroropic mode.
!----------------------------------------------------------------------
!
#ifdef SOLVE3D
!
!  Set vertical S-coordinate functions
!
      call set_scoord
!
!  Set fast-time averaging for coupling of split-explicit baroropic mode.
!

      call set_weights
CR      write(*,*) ' -8' MYID
!
!  Create three-dimensional S-coordinate system,
!  which may be needed by ana_initial
!  (here it is assumed that free surface zeta=0).
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_depth (tile)
      enddo
CR      write(*,*) ' -7' MYID
!
!  Make grid diagnostics
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call grid_stiffness (tile)
      enddo
#endif
!
!----------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables
!----------------------------------------------------------------------
!
!  Read from NetCDF file
!
#ifdef ANA_INITIAL
      if (nrrec.ne.0) then            ! read from NetCDF file
#endif
#ifdef EXACT_RESTART
        call get_initial (nrrec-1, 2) ! Set initial conditions
                                      ! in case of restart
C$OMP BARRIER
# ifdef SOLVE3D
        do tile=0,NSUB_X*NSUB_E-1
          call set_depth (tile)       !<-- needed to initialize Hz_bak
        enddo
C$OMP BARRIER
# endif
#endif
        call get_initial (nrrec, 1)   ! Set initial conditions
#ifdef ANA_INITIAL
      else  ! nrrec.eq.0
# if defined OA_COUPLING || defined OW_COUPLING
        call cpl_prism_define
        oasis_time = 0
        MPI_master_only write(*,*)'CPL-ROMS: OASIS_TIME',oasis_time
# endif
      endif
#endif
                                ! Set initial model clock: at this
      time=start_time           ! moment "start_time" (global scalar)
      tdays=time*sec2day        ! is set by get_initial or analytically
                                ! --> copy it into threadprivate "time"

#ifdef USE_CALENDAR
      ! compute start time and ntimes
      time=start_time +origin_date_in_sec
!# ifdef ANA_INITIAL
!      ! if no initial file set start_date from croco.in
!      if (nrrec.eq.0) time=tool_datosec(start_date)
!# endif
      time_end=tool_datosec(end_date)
      ntimes=int((time_end-time)/dt)
      MPI_master_only write(stdout,*)
     &     'Ntimes from date_start and date_end:',ntimes
      ntimes=ntimes+ntstart
#endif
!
!  Set initial conditions analytically for ideal cases
!  or for tracer variables not present in NetCDF file
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call ana_initial (tile)
      enddo
!
!----------------------------------------------------------------------
!  Initialize specific PISCES variables
!----------------------------------------------------------------------
!
#if defined BIOLOGY && defined PISCES
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call  pisces_ini_tile (tile)
      enddo
#endif
CR      write(*,*) ' -6' MYID
      if (may_day_flag.ne.0) goto 99 !-->  EXIT
!
!----------------------------------------------------------------------
!  Bottom sediment parameters for BBL or SEDIMENT model
!----------------------------------------------------------------------
!
#if (defined BBL && defined ANA_BSEDIM) || defined SEDIMENT
!
!  --- Set analytically ---
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1

# if defined BBL && defined ANA_BSEDIM
#  ifdef AGRIF
        if (Agrif_lev_sedim.EQ.0) call ana_bsedim (tile)
#  else
        call ana_bsedim (tile)
#  endif
# endif
!
# ifdef SEDIMENT
#  ifdef AGRIF
        if (Agrif_lev_sedim.EQ.0) call ana_sediment (tile)
#  else
        call ana_sediment (tile)
#  endif
# endif
      enddo
#endif
!
!  --- Read from NetCDF file (to do) ---
!
#if defined BBL && !defined ANA_BSEDIM && !defined SEDIMENT
# ifdef AGRIF
      if (Agrif_lev_sedim.EQ.0) call get_bsedim
#  else
      call get_bsedim
#  endif
#endif
!
#if defined SEDIMENT && !defined ANA_SEDIMENT
# ifdef AGRIF
      if (Agrif_lev_sedim.EQ.0) call get_sediment
# else
      call get_sediment
# endif
#endif


!----------------------------------------------------------------------
!  SUBSTANCE : computing cell surfaces need for MUSTANG and BIOLink
!----------------------------------------------------------------------
# ifdef SUBSTANCE
      call substance_surfcell
# endif
!
!----------------------------------------------------------------------
!  MUSTANG : initialization
!----------------------------------------------------------------------
!
#ifdef MUSTANG
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call MUSTANG_init_main (tile)
      enddo
CR      write(*,*)  ' -5' MYID0
#endif
!
!----------------------------------------------------------------------
!  Finalize grid setup
!----------------------------------------------------------------------
!
!  Finalize vertical grid now that zeta is knowned
!  zeta is also corrected here for Wetting/Drying
!  in both 2D and 3D cases
!
#if defined SOLVE3D || defined WET_DRY
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_depth (tile)
      enddo
CR      write(*,*)  ' -5' MYID
#endif
!
!----------------------------------------------------------------------
!  Initialize diagnostic fields: mass flux, rho, omega
!----------------------------------------------------------------------
!
#ifdef SOLVE3D
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_HUV (tile)
# ifdef RESET_RHO0
        call reset_rho0 (tile)
# endif
      enddo
CR      write(*,*)  ' -4' MYID
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call omega (tile)
        call rho_eos (tile)
      enddo
CR      write(*,*)  ' -3' MYID
#endif
!
!----------------------------------------------------------------------
!  Initialize surface wave variables
!----------------------------------------------------------------------
!
#ifdef MRL_WCI
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call mrl_wci (tile)
      enddo
CR      write(*,*)  ' -2' MYID
#endif
!
!----------------------------------------------------------------------
!  Set nudging coefficients
!  for sea surface height, momentum and tracers
!----------------------------------------------------------------------
!
# if defined TNUDGING  || defined ZNUDGING  \
  || defined M2NUDGING || defined M3NUDGING \
                       || defined SPONGE
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_nudgcof (tile)
      enddo
# endif
!
!----------------------------------------------------------------------
!  Read initial climatology fields from NetCDF file
!  (or any external oceanic forcing data)
!  in 3D interior or 2D boundary arrays
!----------------------------------------------------------------------
!
#if defined TCLIMATOLOGY && !defined ANA_TCLIMA
      call get_tclima
#endif
#if defined M2CLIMATOLOGY && !defined ANA_M2CLIMA ||\
   (defined M3CLIMATOLOGY && !defined ANA_M3CLIMA)
      call get_uclima
#endif
#if defined ZCLIMATOLOGY && !defined ANA_SSH
      call get_ssh
#endif
#if defined FRC_BRY && !defined ANA_BRY
      call get_bry
# ifdef BIOLOGY
      call get_bry_bio
# endif
#endif
#if !defined ANA_BRY_WKB && defined WKB_WWAVE
      call get_bry_wkb
# endif
!
!----------------------------------------------------------------------
!  Set analytical initial climatology fields
!  (or any external oceanic forcing data)
!  for sea surface height, momentum and tracers
!----------------------------------------------------------------------
!
#if (defined ZCLIMATOLOGY  && defined ANA_SSH)     || \
    (defined M2CLIMATOLOGY && defined ANA_M2CLIMA) || \
    (defined M3CLIMATOLOGY && defined ANA_M3CLIMA) || \
     defined TCLIMATOLOGY
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
# ifdef TCLIMATOLOGY
        call ana_tclima (tile)
# endif
# if defined M2CLIMATOLOGY && defined ANA_M2CLIMA ||\
    (defined M3CLIMATOLOGY && defined ANA_M3CLIMA)
        call ana_uclima (tile)
# endif
# if defined ZCLIMATOLOGY && defined ANA_SSH
        call ana_ssh (tile)
# endif
      enddo
#endif
!
CR      write(*,*) ' -2' MYID
!
!----------------------------------------------------------------------
!  Read surface forcing from NetCDF file
!----------------------------------------------------------------------
!
      call get_vbc
!
!----------------------------------------------------------------------
!  Read tidal harmonics from NetCDF file
!----------------------------------------------------------------------
!
#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
      call get_tides
#endif
!
CR      write(*,*) ' -1' MYID
!
!----------------------------------------------------------------------
! OA "Stand Alone" module : second initialization step (spatial domain)
!----------------------------------------------------------------------
!
#ifdef ONLINE_ANALYSIS
      if ( if_oa.eqv..true. ) then
C$OMP PARALLEL DO PRIVATE(tile)
         do tile=0,NSUB_X*NSUB_E-1
           call online_spectral_diags(tile,0)
         enddo
      endif
#endif
!
!----------------------------------------------------------------------
!  Initialize XIOS I/O server
!----------------------------------------------------------------------
!
#ifdef XIOS
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call init_xios(tile)
      enddo
#endif
!
      if (may_day_flag.ne.0) goto 99 !-->  EXIT
!
!----------------------------------------------------------------------
!  Initialization for Lagrangian floats
!
!  It is done here and not in init_scalars since it must be done only
!  once (whether child levels exist or not)
!----------------------------------------------------------------------
!
#ifdef FLOATS
      nrecflt=0
      ncidflt=-1
      flopspval=1.E15  ! spval is the nodata flag for float variables
      deltac2p=2.3     ! distance from the boundary at which a float
                       ! is transferred from child to parent
      deltap2c=2.5     ! same for transfer from parent to child

      call init_arrays_floats
      call init_floats
# ifdef SPHERICAL
      call interp_r2d_type_ini (lonr(START_2D_ARRAY), iflon)
      call interp_r2d_type_ini (latr(START_2D_ARRAY), iflat)
# else
      call interp_r2d_type_ini (  xr(START_2D_ARRAY), iflon)
      call interp_r2d_type_ini (  yr(START_2D_ARRAY), iflat)
# endif
# ifdef SOLVE3D
      call fill_ini ! fills in trackaux for ixgrd,iygrd,izgrd
                    ! and ifld (either izgrd or ifld is meaningful)
# endif
      if (ldefflt) call wrt_floats
#endif /* FLOATS */
!
!----------------------------------------------------------------------
!  Initialization for stations
!
!  It is done here and not in init_scalars since it must be done only
!  once (whether child levels exist or not)
!----------------------------------------------------------------------
!
#ifdef STATIONS
      nrecsta=0
      ncidsta=-1
      staspval=1.E15  ! nodata flag for float variables.
      stadeltap2c=2.5 ! distance from the boundary at which a
                      ! float is transfered from parent to child
      call init_arrays_sta
      call init_sta
# ifdef SPHERICAL
      call interp_r2d_sta_ini (lonr(START_2D_ARRAY), istalon)
      call interp_r2d_sta_ini (latr(START_2D_ARRAY), istalat)
# else
      call interp_r2d_sta_ini (  xr(START_2D_ARRAY), istalon)
      call interp_r2d_sta_ini (  yr(START_2D_ARRAY), istalat)
# endif
# ifdef SOLVE3D
      call fill_sta_ini ! fills in trackaux for ixgrd,iygrd,izgrd
                        ! and ifld (either izgrd or ifld is meaningful)
# endif
      if (ldefsta) call wrt_sta
#endif /* STATIONS */
!
!----------------------------------------------------------------------
!  WKB surface wave model:
!
!  initialization and spinup to equilibrium
!----------------------------------------------------------------------
!
#ifdef WKB_WWAVE
C$OMP BARRIER
C$OMP MASTER
        MPI_master_only write(stdout,'(/1x,A/)')
     &     'WKB: started steady wave computation.'
C$OMP END MASTER
      iic=0
      winfo=1
      iwave=1
      thwave=1.D+10
# ifndef ANA_BRY_WKB
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_bry_wkb (tile)   ! set boundary forcing
      enddo
# endif
# ifdef MRL_CEW
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call wkb_cew_prep (tile)  ! prepare coupling mode
      enddo
C$OMP BARRIER
      wint=0
      do winterp=1,interp_max
        wint=wint+1
        if (wint.gt.2) wint=1
C$OMP PARALLEL DO PRIVATE(tile)
        do tile=0,NSUB_X*NSUB_E-1
          call wkb_uvfield (tile, winterp)
        enddo
      enddo
C$OMP BARRIER
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call wkb_cew_finalize (tile)
      enddo
C$OMP BARRIER
# endif /* MRL_CEW */
!
!  Spinup: intergrate wave model to equilibrium
!
      do while (iwave.le.50000.and.thwave.ge.1.D-10)
        wstp=wnew
        wnew=wstp+1
        if (wnew.ge.3) wnew=1
C$OMP PARALLEL DO PRIVATE(tile)
        do tile=0,NSUB_X*NSUB_E-1
# ifdef WAVE_OFFLINE
          if (iwave.eq.1) call set_wwave(tile)
# endif
          call wkb_wwave (tile)
        enddo
        call wkb_diag (0)
        iwave=iwave+1
        thwave=max(av_wac,av_wkn)
      enddo
# if defined RVTK_DEBUG || defined RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
      call check_tab2d(wac(:,:,wnew),'wac initialisation #1','r')
C$OMP END MASTER
# endif
!
!  Re-initialize wave forcing terms
!
      first_time=0
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call mrl_wci (tile)
      enddo
C$OMP BARRIER
C$OMP MASTER
        MPI_master_only write(stdout,'(/1x,A/)')
     &     'WKB: completed steady wave computation.'
C$OMP END MASTER
#endif /* WKB_WWAVE */

!
!----------------------------------------------------------------------
!  Set initial non-Boussinesq (or fast 3D) parameters and variables
!----------------------------------------------------------------------
!
#ifdef M3FAST
# ifdef NBQ_MASS
!
! Re-evaluate Hz and Huon,Hvom using
! previously computed density rho
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_depth (tile)
      enddo

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call set_HUV (tile)
      enddo
# endif
!
! Set initial NBQ param. & var.
!
C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call initial_nbq(tile)
      enddo
#endif
!
!----------------------------------------------------------------------
! OA "Stand Alone" module
!----------------------------------------------------------------------
!
#ifdef ONLINE_ANALYSIS
      if ( if_oa.eqv..true. ) then
C$OMP PARALLEL DO PRIVATE(tile)
         do tile=0,NSUB_X*NSUB_E-1
           call online_spectral_diags(tile,1)
         enddo
        call output_oa
      endif
#endif
!
!----------------------------------------------------------------------
!  Write initial fields into history NetCDF files
!----------------------------------------------------------------------
!
#ifdef XIOS
      if (nrrec.eq.0) then
C$OMP PARALLEL DO PRIVATE(tile)
        do tile=0,NSUB_X*NSUB_E-1
          call send_xios_diags(tile)
        enddo
      endif
#else
      if (ldefhis .and. wrthis(indxTime)) call wrt_his
#endif

CR      write(*,*) '  0' MYID
      if (may_day_flag.ne.0) goto 99     !-->  EXIT
!
!
!**********************************************************************
!                                                                     *
!             *****   ********   *****   ******  ********             *
!            **   **  *  **  *  *   **   **   ** *  **  *             *
!            **          **    **   **   **   **    **                *
!             *****      **    **   **   **   *     **                *
!                 **     **    *******   *****      **                *
!            **   **     **    **   **   **  **     **                *
!             *****      **    **   **   **   **    **                *
!                                                                     *
!**********************************************************************
!
!
#undef CR
      MPI_master_only write(stdout,'(/1x,A27/)')
     &                'MAIN: started time-stepping.'
      next_kstp=kstp
      time_start=time
#ifdef USE_CALENDAR
      time_mars=time
      time_end=tool_datosec(end_date)
#endif

! XIOS (
!  Clean log output a bit
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
! XIOS )

#ifdef SOLVE3D
      iif = -1
      nbstep3d = 0
#endif
      iic = ntstart

#ifdef AGRIF
      iind = -1
      grids_at_level = -1
      sortedint = -1
      call computenbmaxtimes
#endif

#ifdef MPI_TIME
! Start Chrono
      start_time1=PMPI_Wtime()
#endif

      do iicroot=ntstart,ntimes+1
#ifdef MPI_TIME
        start_time1_1=PMPI_Wtime()
#endif

#ifdef USE_CALENDAR
        if (mod(iicroot-1,ninfo) .eq. 0) then
          MPI_master_only write(stdout,'(a)') tool_sectodat(time_mars)
        endif
        if (time .gt. time_end) goto 99
#endif

#ifdef SOLVE3D
# ifndef AGRIF
        do iifroot = 0,nfast+2
# else
        nbtimes = 0
        nbmaxprttime = ntimes+1
        nbprttime = iicroot
        do while (nbtimes.LE.nbmaxtimes)
# endif
#endif

#ifdef AGRIF
          call Agrif_Step(step)
#else
          call step()
#endif

#ifdef SOLVE3D
        enddo
#endif
#ifdef MPI_TIME
        start_time1_2=PMPI_Wtime()
        exe_time_1 = start_time1_2 - start_time1_1
          if (mynode.eq.0) write(97,*) exe_time_1
#endif
        if (may_day_flag.ne.0) goto 99     !-->  EXIT
      enddo                                !-->  end of time step

  99  continue                             ! SHUTDOWN:
#ifdef MPI_TIME
      start_time2=PMPI_Wtime()
      exe_time = start_time2 - start_time1
      MPI_master_only write(stdout,*) '*******************************'
      MPI_master_only write(stdout,*) 'time time-stepping : ',exe_time
      MPI_master_only write(stdout,*) '*******************************'
#endif
C$OMP PARALLEL DO PRIVATE(tile)            ! Stop timers and
      do tile=0,NSUB_X*NSUB_E-1            ! close netCDF files.
        call stop_timers()
      enddo
      call closecdf
#ifdef XIOS
      call iom_context_finalize( "crocox")   ! needed for XIOS+AGRIF
#endif
#ifdef AGRIF
!
!  Close the netcdf files also for the child grids
!
      parcours=>Agrif_Curgrid%child_list % first
      do while (associated(parcours))
        Call Agrif_Instance(parcours % gr)
        call closecdf
# ifdef XIOS
        call iom_context_finalize( "crocox") ! needed for XIOS+AGRIF
# endif
        parcours => parcours % next
      enddo
#endif

 100  continue

#ifdef MPI
      if (ierr.ne.0) call mpi_abort (MPI_COMM_WORLD, ierr)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! XIOS

!     start_time2=PMPI_Wtime()
!     exe_time = start_time2 - start_time1
!     if (mynode.eq.0) print *,'exe_time =',exe_time

# if defined XIOS
                                ! case XIOS + (OASIS / no OASIS)
                                !      > MPI finalize done by XIOS
      call xios_finalize()      !      > if OASIS, finalize is done by XIOS
                                !      > if AGRIF : done using iom_context_finalize
!!!#   if !defined AGRIF
!!!      call MPI_Finalize (ierr)  !       if No AGRIF : MPI_Finalize is needed
!!!#   endif

#  if  !defined OA_COUPLING && !defined OW_COUPLING  && !defined AGRIF
       call MPI_Finalize (ierr)      !  if No coupling + No AGRIF : MPI_Finalize is needed
#  endif

# elif defined OA_COUPLING || defined OW_COUPLING
                                          ! case no XIOS + OASIS
      call prism_terminate_proto(ierr)    !   > Finalize OASIS3 (without XIOS)
!!!      call MPI_Finalize (ierr)         !   > Finalize CROCO (without XIOS)
# else
                                          ! case no XIOS + no OASIS
      call MPI_Finalize (ierr)            !   > Finalize CROCO (without XIOS)
# endif
#endif

      stop
      end
