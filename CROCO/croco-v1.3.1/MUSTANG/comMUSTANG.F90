#include "cppdefs.h"

MODULE comMUSTANG

#ifdef MUSTANG

!!============================================================================
!! ***  MODULE  comMUSTANG  ***
!! Purpose : declare all common variables related to sediment dynamics
!!============================================================================

!! * Modules used
USE comsubstance ! for lchain, rsh, rlg, riosh

implicit none

! default
public

#include "coupler_define_MUSTANG.h"
  
    !! * Shared or public variables for MUSTANG 

    ! parameters
    REAL(kind=rlg), PARAMETER :: epsi30_MUSTANG = 1.e-30
    REAL(kind=rlg), PARAMETER :: epsilon_MUSTANG = 1.e-09 
    REAL(kind=rsh), PARAMETER :: valmanq = 999.0
    REAL(kind=riosh), PARAMETER :: rg_valmanq_io = 999.0

    ! namelists

    ! namsedim_init
    CHARACTER(len=19) :: date_start_dyninsed ! starting date for dynamic 
        ! processes in sediment; format '01/01/0000 00:00:00'
    CHARACTER(len=19) :: date_start_morpho ! starting date for morphodynamic; 
        ! format '01/01/0000 00:00:00'
    LOGICAL  :: l_repsed ! set to .true. if sedimentary variables are 
        ! initialized from a previous run
    LOGICAL  :: l_initsed_vardiss !set to .true. if initialization of dissolved
        ! variables, temperature and salinity in sediment (will be done with 
        ! concentrations in water at bottom (k=1))
    LOGICAL  :: l_unised !set to .true. for a uniform bottom initialization
    LOGICAL  :: l_init_hsed ! set to .true. if we want to adjust the sediment 
        ! thickness in order to be coherent with sediment 
        ! parameters (calculation of a new hseduni based on 
        ! cseduni, cvolmax values, and csed_ini of each sediment)
    CHARACTER(len=lchain) :: filrepsed ! file path from which the model is 
        ! initialized for the continuation of a previous run
    CHARACTER(len=lchain) :: fileinised ! file path for initialization (if 
        ! l_unised is False)
    REAL(KIND=rsh) :: cseduni ! initial sediment concentration (kg/m3)
    REAL(KIND=rsh) :: hseduni ! initial uniform sediment thickness (m) 
    REAL(KIND=rsh) :: csed_mud_ini ! real, mud concentration into initial 
        ! sediment (kg/m3) (if = 0. ==> csed_mud_ini = cfreshmud)
    INTEGER        :: ksmiuni ! lower grid cell index in the sediment
    INTEGER        :: ksmauni ! upper grid cell index in the sediment
    REAL(KIND=rsh) :: sini_sed ! initial interstitial water uniform salinity
    REAL(KIND=rsh) :: tini_sed ! initial interstitial water uniform temperature
    REAL(KIND=rsh) :: poro_mud_ini !only if key_MUSTANG_V2, initial porosity of
        ! mud fraction


    ! namsedim_layer
    LOGICAL  :: l_dzsmaxuni ! set to .true. dzsmax = dzsmaxuni , 
        ! if set to .false. then linearly computed in MUSTANG_sedinit
        ! from dzsmaxuni to dzsmaxuni/100 depending on water depth
    LOGICAL  :: l_dzsminuni !only if key_MUSTANG_V2, set to .false. if dzsmin 
        ! vary with sediment bed composition, else dzsmin =  dzsminuni
    REAL(KIND=rsh) :: dzsminuni !only if key_MUSTANG_V2, minimum sediment 
        ! layer thickness (m)
    REAL(KIND=rsh) :: dzsmin ! minimum sediment layer thickness (m)
    REAL(KIND=rsh) :: dzsmaxuni ! uniform maximum thickness for the superficial
        ! sediment layer (m), must be >0
    REAL(KIND=rsh) :: dzsmax_bottom ! maximum thickness of bottom layers 
        ! which result from the fusion when ksdmax is exceeded (m)
    REAL(KIND=rsh) :: k1HW97 !only if key_MUSTANG_V2, 
        ! ref value k1HW97 = 0.07, parameter to compute active layer 
        ! thickness (Harris and Wiberg, 1997)
    REAL(KIND=rsh) :: k2HW97 !only if key_MUSTANG_V2
        ! rref value k2HW97 = 6.0,parameter to compute active layer 
        ! thickness (Harris and Wiberg, 1997)
    REAL(KIND=rsh) :: fusion_para_activlayer !only if key_MUSTANG_V2
        ! criterion cohesiveness for fusion in active layer
        ! 0 : no fusion, 
        ! = 1 : frmudcr1, 
        ! > 1 : between frmudcr1 & frmudcr2
    INTEGER :: nlayer_surf_sed ! number of layers below the sediment surface 
        ! that can not be melted (max thickness = dzsmax)


    ! namsedim_bottomstress
    LOGICAL  :: l_z0seduni ! boolean, set to .false. for z0sed computation from
        ! sediment diameter (if true, z0seduni is used)
    REAL(KIND=rsh) :: z0seduni ! uniform bed roughness (m)
    REAL(KIND=rsh) :: z0sedmud ! mud (i.e.minimum) bed roughness (m) 
        ! (used only if l_unised is false)
    REAL(KIND=rsh) :: z0sedbedrock ! bed roughness for bedrock (no sediment) (m) 
        ! (used only if l_unised is false)
    LOGICAL :: l_fricwave ! boolean, set to .true. if using wave related friction 
        ! factor for bottom shear stress (from wave orbital velocity and period)
        ! if .false. then fricwav namelist value is used
    REAL(KIND=rsh) :: fricwav ! default value is 0.06, wave related friction 
        !factor (used for bottom shear stress computation)
    LOGICAL :: l_z0hydro_coupl_init ! boolean, set to .true. if evaluation of 
        ! z0 hydro depends on sediment composition at the beginning 
        ! of the simulation
    LOGICAL :: l_z0hydro_coupl ! boolean, set to .true. if evaluation of 
        ! z0 hydro depends on sediment composition along the run
    REAL(KIND=rsh) :: coef_z0_coupl ! parameter to compute z0hydro in the 
        ! first centimeter : z0hydro = coef_z0_coupl * sand diameter 
    REAL(KIND=rsh) :: z0_hydro_mud ! z0hydro if pure mud (m)
    REAL(KIND=rsh) :: z0_hydro_bed ! z0hydro if no sediment (m)


    ! namsedim_deposition
    REAL(KIND=rsh) :: cfreshmud ! fresh deposit concentration (kg/m3) 
        ! (must be around 100 if consolidation 
        ! or higher (300-500 if no consolidation)
    REAL(KIND=rsh) :: csedmin ! concentration of the upper layer under 
        ! which there is fusion with the underlying sediment cell (kg/m3)
    REAL(KIND=rsh) :: cmudcr ! critical relative concentration of the surface 
        ! layer above which no mixing is allowed with the underlying 
        ! sediment (kg/m3)
    REAL(KIND=rsh) :: aref_sand  ! parameter used in sandconcextrap,
        ! reference height above sediment, used for computing of 
        ! sand deposit. Parameter used for sand extrapolation on water 
        ! column and correct sand transport, value by default = 0.02 
        ! correspond to Van Rijn experiments 
        ! DO NOT CHANGED IF NOT EXPERT
    REAL(KIND=rsh) :: cvolmaxsort ! max volumic concentration of sorted sand
    REAL(KIND=rsh) :: cvolmaxmel ! maxvolumic concentration of mixed sediments
    REAL(KIND=rsh) :: slopefac !slope effect multiplicative on deposit 
        ! (only if key_MUSTANG_slipdeposit)


    ! namsedim_erosion
    REAL(KIND=rsh) :: activlayer ! active layer thickness (m)
    REAL(KIND=rsh) :: frmudcr2 ! critical mud fraction under which the 
        ! behaviour is intermediate between sand and mud and above which the 
        ! behavior is purely muddy
    REAL(KIND=rsh) :: coef_frmudcr1 ! to compute critical mud fract. frmudcr1 
        ! underwhich the behaviour is purely sandy 
        ! (frmudcr1=min(coef_frmudcr1*d50 sand,frmudcr2))
    REAL(KIND=rsh) :: x1toce_mud ! coef. for the formulation of the critical 
        ! erosion stress in mud behavior toce=x1toce*csed**x2toce
    REAL(KIND=rsh) :: x2toce_mud ! coef. for the formulation of the critical 
    ! erosion stress in mud behavior toce=x1toce*csed**x2toce
    REAL(KIND=rsh) :: E0_sand_para ! coefficient used to modulate erosion 
        ! flux for sand (=1 if no correction )
    REAL(KIND=rsh) :: n_eros_sand ! parameter for erosion flux for sand 
        ! (E0_sand*(tenfo/toce-1.)**n_eros_sand )
        ! WARNING : choose parameters compatible with E0_sand_option 
        ! (example : n_eros_sand=1.6 for E0_sand_option=1)
    REAL(KIND=rsh) :: E0_mud ! erosion flux for mud
    REAL(KIND=rsh) :: n_eros_mud ! E0_mud*(tenfo/toce-1.)**n_eros_mud
    INTEGER        :: ero_option ! choice of erosion formulation for mixing 
        ! sand-mud
        ! ero_option= 0 : pure mud behavior 
        ! ero_option= 1 : linear interpolation between sand and mud behavior, 
        !   depend on proportions of the mixture
        ! ero_option= 2 : formulation derived from that of J.Vareilles (2013)
        ! ero_option= 3 : formulations proposed by B. Mengual (2015) with 
        !   exponential coefficients depend on proportions of the mixture
    INTEGER        :: E0_sand_option ! integer, choice of formulation for 
        ! E0_sand evaluation :
        ! E0_sand_option = 0 E0_sand = E0_sand_Cst 
        ! E0_sand_option = 1 E0_sand evaluated with Van Rijn (1984) 
        ! E0_sand_option = 2 E0_sand evaluated with erodimetry 
        !    (min(0.27,1000*d50-0.01)*toce**n_eros_sand)
        ! E0_sand_option = 3 E0_sand evaluated with Wu and Lin (2014)
    REAL(KIND=rsh) :: xexp_ero !used only if ero_option=3 : adjustment on 
        ! exponential variation  (more brutal when xexp_ero high)
    REAL(KIND=rsh) :: E0_sand_Cst ! constant erosion flux for sand 
        ! (used if E0_sand_option= 0) 
    REAL(KIND=rsh) :: E0_mud_para_indep !only if key_MUSTANG_V2,
        ! parameter to correct E0_mud in case of erosion 
        ! class by class in non cohesive regime
    LOGICAL        :: l_peph_suspension !only if key_MUSTANG_V2,
        ! set to .true. if hindering / exposure processes in critical 
        ! shear stress estimate for suspension
    LOGICAL        :: l_eroindep_noncoh !only if key_MUSTANG_V2,
        ! set to .true. in order to activate independant erosion for 
        ! the different sediment classes sands and muds  
        ! set to .false. to have the mixture mud/sand eroded as in V1
    LOGICAL        :: l_eroindep_mud !only if key_MUSTANG_V2,
        ! set to .true. if mud erosion independant for sands erosion
        ! set to .false. if mud erosion proportionnal to total sand erosion
    LOGICAL        :: l_xexp_ero_cst !only if key_MUSTANG_V2, set to .true. 
        ! if xexp_ero estimated from empirical formulation, depending on 
        ! frmudcr1 
    INTEGER        :: tau_cri_option !only if key_MUSTANG_V2, 
        ! choice of critical stress formulation , 
        ! 0: Shields 1: Wu and Lin (2014)
    INTEGER        :: tau_cri_mud_option_eroindep !only if key_MUSTANG_V2
        ! choice of mud critical stress formulation 
        ! 0: x1toce_mud*cmudr**x2toce_mud
        ! 1: toce_meansan if somsan>eps (else->case0)
        ! 2: minval(toce_sand*cvsed/cvsed+eps) if >0 (else->case0)
        ! 3: min( case 0; toce(isand2) )


#ifdef key_MUSTANG_V2
    ! namsedim_poro 
    INTEGER :: poro_option ! choice of porosity formulation
        ! 1: Wu and Li (2017) (incompatible with consolidation))
        ! 2: mix ideal coarse/fine packing 
    REAL(KIND=rsh) :: Awooster ! parameter of the formulation of 
        ! Wooster et al. (2008) for estimating porosity associated to the 
        ! non-cohesive sediment see Cui et al. (1996) ref value = 0.42
    REAL(KIND=rsh) :: Bwooster ! parameter of the formulation of 
    ! Wooster et al. (2008) for estimating porosity associated to the 
    ! non-cohesive sediment see Cui et al. (1996) ref value = -0,458
    REAL(KIND=rsh) :: Bmax_wu ! maximum portion of the coarse sediment class 
        ! participating in filling , ref value = 0.65
    REAL(KIND=rsh) :: poro_min ! minimum porosity below which consolidation 
        ! is stopped
#endif


#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
    ! namsedim_bedload 
    LOGICAL :: l_peph_bedload ! set to .true. if hindering / exposure processes
        ! in critical shear stress estimate for bedload
    LOGICAL :: l_slope_effect_bedload ! set to .true. if accounting for slope 
        ! effects in bedload fluxes (Lesser formulation)
    LOGICAL :: l_fsusp ! limitation erosion fluxes of non-coh sediment in case 
        ! of simultaneous bedload transport, according to Wu & Lin formulations
        ! set to .true. if erosion flux is fitted to total transport 
        ! should be set to .false. if E0_sand_option=3 (Wu & Lin)
    REAL(KIND=rsh) :: alphabs ! coefficient for slope effects (default 
        ! coefficients Lesser et al. (2004), alphabs = 1.)
    REAL(KIND=rsh) :: alphabn ! coefficient for slope effects (default 
        ! coefficients Lesser et al. (2004), default alphabn is 1.5 but 
        ! can be higher, until 5-10 (Gerald Herling experience))
    REAL(KIND=rsh) :: hmin_bedload  ! no bedload in u/v directions if 
        ! h0+ssh <= hmin_bedload in neighbouring cells
#endif


!**TODO** put under cpp key #ifdef key_MUSTANG_lateralerosion
    ! namsedim_lateral_erosion 
    REAL(KIND=rsh) :: htncrit_eros ! critical water height so as to prevent 
        ! erosion under a given threshold (the threshold value is different for
        ! flooding or ebbing, cf. Hibma's PhD, 2004, page 78)
    REAL(KIND=rsh) :: coef_erolat ! slope effect multiplicative factor 
    REAL(KIND=rsh) :: coef_tauskin_lat ! parameter to evaluate the lateral 
        ! stress as a function of the average tangential velocity on the 
        ! vertical
    LOGICAL        :: l_erolat_wet_cell ! set to .true in order to take 
        ! into account wet cells lateral erosion
!**TODO** put under cpp key #key_MUSTANG_lateralerosion


    ! namsedim_consolidation 
    LOGICAL        :: l_consolid ! set to .true. if sediment consolidation is 
        ! accounted for
    REAL(KIND=rsh) :: dt_consolid ! time step for consolidation processes
    REAL(KIND=rlg) :: subdt_consol ! sub time step for consolidation and 
                                   ! particulate bioturbation  in sediment
    REAL(KIND=rsh) :: csegreg ! NOT CHANGE VALUE if not expert, default 250.0
    REAL(KIND=rsh) :: csandseg ! NOT CHANGE VALUE if not expert, default 1250.0
    REAL(KIND=rsh) :: xperm1 ! permeability=xperm1*d50*d50*voidratio**xperm2
    REAL(KIND=rsh) :: xperm2 ! permeability=xperm1*d50*d50*voidratio**xperm2
    REAL(KIND=rsh) :: xsigma1 ! parameter used in Merckelback & Kranenburg s 
        ! (2004) formulation NOT CHANGE VALUE if not expert, default 6.0e+05
    REAL(KIND=rsh) :: xsigma2 ! parameter used in Merckelback & Kranenburg s 
        ! (2004) formulation NOT CHANGE VALUE if not expert, default 6


    ! namsedim_diffusion     
    LOGICAL        :: l_diffused ! set to .true. if taking into account 
        ! dissolved diffusion in sediment and at the water/sediment interface
    REAL(KIND=rsh) :: dt_diffused ! time step for diffusion in sediment
    INTEGER        :: choice_flxdiss_diffsed ! choice for expression of 
        ! dissolved fluxes at sediment-water interface
        ! 1 : Fick law : gradient between Cv_wat at dz(1)/2
        ! 2 : Fick law : gradient between Cv_wat at distance epdifi
    REAL(KIND=rsh) :: xdifs1 ! diffusion coefficients within the sediment
    REAL(KIND=rsh) :: xdifsi1 ! diffusion coefficients at the water-sediment 
        ! interface
    REAL(KIND=rsh) :: epdifi ! diffusion thickness in the water at the 
        ! sediment-water interface
    REAL(KIND=rsh) :: fexcs ! factor of eccentricity of concentrations in 
        ! vertical fluxes evaluation (.5 a 1) 


    ! namsedim_bioturb  
    LOGICAL        :: l_bioturb ! set to .true. if taking into account 
        ! particulate bioturbation (diffusive mixing) in sediment
    LOGICAL        :: l_biodiffs ! set to .true. if taking into account 
        ! dissolved bioturbation diffusion in sediment
    REAL(KIND=rsh) :: dt_bioturb ! time step for bioturbation in sediment
    REAL(KIND=rsh) :: subdt_bioturb ! sub time step for bioturbation 
    REAL(KIND=rsh) :: xbioturbmax_part ! max particular bioturbation 
        ! coefficient by bioturbation Db (in surface)
    REAL(KIND=rsh) :: xbioturbk_part ! for part. bioturbation coefficient 
        ! between max Db at sediment surface and 0 at bottom
    REAL(KIND=rsh) :: dbiotu0_part ! max depth beneath the sediment 
        ! surface below which there is no bioturbation
    REAL(KIND=rsh) :: dbiotum_part ! sediment thickness where the 
        ! part-bioturbation coefficient Db is constant (max)
    REAL(KIND=rsh) :: xbioturbmax_diss ! max diffusion coeffient by 
        ! biodiffusion Db (in surface)
    REAL(KIND=rsh) :: xbioturbk_diss ! coef (slope) for biodiffusion 
        ! coefficient between max Db at sediment surface and 0 at bottom
    REAL(KIND=rsh) :: dbiotu0_diss ! max depth beneath the sediment 
        ! surface below which there is no bioturbation
    REAL(KIND=rsh) :: dbiotum_diss ! sediment thickness where the  
        ! diffsolved-bioturbation coefficient Db is constant (max)
    REAL(KIND=rsh) :: frmud_db_min ! mud fraction limit (min) below which 
        ! there is no Biodiffusion
    REAL(KIND=rsh) :: frmud_db_max ! mud fraction limit (max)above which 
        ! the biodiffusion coefficient Db is maximum (muddy sediment)


    ! namsedim_morpho   
    LOGICAL :: l_morphocoupl ! set to .true if coupling module morphodynamic  
    LOGICAL :: l_MF_dhsed ! set to .true. if morphodynamic applied with 
        ! sediment height variation amplification 
        ! (MF_dhsed = MF; then MF will be = 0)
        ! set to .false. if morphodynamic is applied with 
        ! erosion/deposit fluxes amplification (MF_dhsed not used)
    LOGICAL :: l_bathy_actu ! set to .true. if reading a new bathy issued a 
        ! previous run and saved in filrepsed (given in namelist namsedim_init)  
        !!! NOT IMPLEMENTED YET !!! **TODO**
    REAL(KIND=rsh) :: MF ! morphological factor : multiplication factor for 
        ! morphologicalevolutions, equivalent to a "time acceleration" 
        ! (morphological evolutions over a MF*T duration are assumed to be 
        ! equal to MF * the morphological evolutions over T). 
    REAL(KIND=rlg) :: dt_morpho ! time step for morphodynamic (s)


#if  ! defined key_noTSdiss_insed
    ! namtempsed  
    REAL(KIND=rsh) :: mu_tempsed1 ! parameters used to estimate thermic 
        ! diffusitiyfunction of mud fraction 
    REAL(KIND=rsh) :: mu_tempsed2 ! parameters used to estimate thermic
        ! diffusitiyfunction of mud fraction 
    REAL(KIND=rsh) :: mu_tempsed3 ! parameters used to estimate thermic
        ! diffusitiyfunction of mud fraction 
    REAL(KIND=rsh) :: epsedmin_tempsed ! sediment thickness limits for 
         ! estimation heat loss at bottom, if hsed < epsedmin_tempsed : 
        ! heat loss at sediment bottom = heat flux a sediment surface
    REAL(KIND=rsh) :: epsedmax_tempsed ! sediment thickness limits for
        ! estimation heat loss at bottom, if hsed > epsedmax_tempsed : 
        ! heat loss at sediment bottom = 0.
#endif


    ! namsedoutput  
    LOGICAL :: l_outsed_flx_WS_all ! set to .true. if output fluxes threw 
        ! interface Water/sediment (2 2D variables per 
        ! constitutive particulate variable)
    LOGICAL :: l_outsed_flx_WS_int ! set to .true. if output fluxes threw 
        ! interface Water/sediment (integration on all 
        ! constitutive particulate variables)
    LOGICAL :: l_outsed_saltemp ! set to .true. if output Salinity and 
        ! Temperature in sediment
    INTEGER :: nk_nivsed_out ! number of saved sediment layers 
        ! unused if choice_nivsed_out = 1                     
        ! <ksdmax if choice_nivsed_out = 2, 
        ! unused if choice_nivsed_out = 3
        !  <6 if choice_nivsed_out = 4, 
    INTEGER :: choice_nivsed_out ! choice of saving output  (1 to 4)
    REAL(KIND=rsh), DIMENSION(5) :: ep_nivsed_out ! 5 values of sediment 
        ! layer thickness (mm), beginning with surface layer 
        ! (used if choice_nivsed_out=4)
    REAL(KIND=rsh) :: epmax_nivsed_out ! maximum thickness (mm) for 
        ! output each layers of sediment (used if choice_nivsed_out=3). 
        ! Below the layer which bottom level exceed this thickness, 
        ! an addition layer is an integrative layer till bottom


#if defined key_MUSTANG_debug && defined key_MUSTANG_V2
    ! namsedim_debug 
    LOGICAL        :: l_debug_effdep ! set to .true. if print some 
        ! informations for debugging MUSTANG deposition
    LOGICAL        :: l_debug_erosion ! set to .true. if  print 
        ! informations for debugging  in erosion routines
    REAL(kind=rlg)   :: lon_debug, lat_debug ! define mesh location 
        ! where we print these informations
    INTEGER          :: i_MUSTANG_debug, j_MUSTANG_debug ! indexes 
        ! of the mesh where we print these informations 
        ! (only if lon_debug and lat_debug = 0.)
    CHARACTER(len=19):: date_start_debug ! starting date for write 
        ! debugging informations 
#endif



#ifdef key_MUSTANG_flocmod
    ! namflocmod  
    LOGICAL :: l_ASH ! set to .true. if aggregation by shear
    LOGICAL :: l_ADS ! set to .true. if aggregation by differential settling
    LOGICAL :: l_COLLFRAG ! set to .true. if fragmentation by collision
    INTEGER :: f_ero_iv ! fragment class (mud variable index corresponding to 
        ! the eroded particle size - typically 1)
    REAL(KIND=rsh) :: f_ater ! ternary fragmentation factor : proportion of 
        ! flocs fragmented as half the size of the initial binary fragments 
        ! (0.0 if full binary fragmentation, 0.5 if ternary fragmentation)
    REAL(KIND=rsh) :: f_dmin_frag ! minimum diameter for fragmentation 
        ! (default 10e-6 microns)
    REAL(KIND=rsh) :: f_ero_frac ! floc erosion (% of floc mass eroded) 
        ! (default 0.05)
    REAL(KIND=rsh) :: f_ero_nbfrag ! nb of fragments produced by erosion 
        ! (default 2.0)
    REAL(KIND=rsh) :: f_mneg_param ! negative mass after 
        ! flocculation/fragmentation allowed before redistribution 
        ! (default 0.001 g/l)
    REAL(KIND=rsh) :: f_collfragparam ! fraction of shear aggregation leading 
        ! to fragmentation by collision (default 0.0, must be less than 1.0)
    REAL(KIND=rsh) :: f_cfcst ! fraction of mass lost by flocs if fragmentation
        ! by collision .. (default : =3._rsh/16._rsh)
    REAL(KIND=rsh) :: f_fp ! relative depth of inter particle penetration  
        ! (default =0.1) (McAnally, 1999)
    REAL(KIND=rsh) :: f_fy ! floc yield strength  (default= 1.0e-10) 
        ! (Winterwerp, 2002)
    REAL(KIND=rsh) :: f_dp0 ! primary particle size (default 4.e-6 m)
    REAL(KIND=rsh) :: f_alpha ! flocculation efficiency parameter 
        ! (default 0.15)
    REAL(KIND=rsh) :: f_beta ! floc break up parameter (default 0.1)
    REAL(KIND=rsh) :: f_nb_frag ! nb of fragments of equal size by shear 
        ! fragmentation (default 2.0 as binary fragmentation)
    REAL(KIND=rsh) :: f_nf ! fractal dimension (default 2.0, usual range from 
        ! 1.6 to 2.8)
    REAL(KIND=rsh) :: f_clim ! min concentration below which flocculation 
        !processes are not calculated
#endif


! end namelist variables


    REAL(KIND=rsh) :: h0fond  ! RESIDUAL_THICKNESS_WAT see coupler_define (m)

    ! fwet =1 if not used  
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: fwet

    ! Initialization 
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: cini_sed
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: cv_sedini
    REAL(KIND=rsh) :: hsed_new

    ! Fluxes at the interface water-sediment
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_s2w
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_sum
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: EROS_FLUX_s2w
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: SETTL_FLUX_w2s
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: SETTL_FLUXSUM_w2s

    ! Sediment parameters
    REAL(KIND=rsh)            :: ros_sand_homogen 
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: typart
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: diamstar
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: ws_sand
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: rosmrowsros
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: stresscri0
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: tetacri0


    INTEGER :: nv_use
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ksmi
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ksma
    REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE   :: cv_sed
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: c_sedtot
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: poro
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: dzs       
    REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: dzsmax
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: gradvit       

    ! Sediment height
    REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: hsed
    REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: hsed0
    REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE       :: hsed_previous

    ! Bottom stress variables
    REAL(KIND=rsh) :: fws2  ! fricwav/2   
    REAL(KIND=rsh),DIMENSION(:,:),ALLOCATABLE     :: z0sed ! roughness (m)
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: tauskin ! max bottom stress due to the combinaison current/wave (N.m-2)
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: tauskin_c ! bottom stress due to current (N.m-2)
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: tauskin_w ! bottom stress due to wave (N.m-2)
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: tauskin_x ! bottom stress - component on x axis
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: tauskin_y ! bottom stress - component on y axis
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: ustarbot ! (m/s)
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: raphbx ! adim.
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: raphby ! adim.

    REAL(KIND=rlg), DIMENSION(:,:), ALLOCATABLE   :: phieau_s2w
    REAL(KIND=rlg), DIMENSION(:,:), ALLOCATABLE   :: phieau_s2w_consol
    REAL(KIND=rlg), DIMENSION(:,:), ALLOCATABLE   :: phieau_s2w_drycell

    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: htot
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: alt_cw1

    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: sal_bottom_MUSTANG
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: temp_bottom_MUSTANG
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: epn_bottom_MUSTANG
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: cw_bottom_MUSTANG
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: ws3_bottom_MUSTANG ! settling velocities in  bottom cell (m/s)
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE   :: roswat_bot

!**TODO** put under cppkey MUSTANG_CORFLUX
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: corflux
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE     :: corfluy

#ifdef key_sand2D
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: rouse2D ! Rouse2D number
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: sum_tmp ! SUM(dzcche*((htot-hzed)/hzed)**rouse) 
#endif


    ! Dynamic in sediment (consolidation/diffusion/bioturbation)
    REAL(KIND=rlg)   :: tstart_dyninsed ! time beginning dynamic in sediment
    REAL(KIND=rlg)   :: t_dyninsed      ! time of next dynamic in sediment step
    REAL(KIND=rlg)   :: dt_dyninsed     ! time step for dynamic in sediment (min of dt for each process)
    LOGICAL :: l_dyn_insed ! true if (l_consolid .OR. l_bioturb .OR. l_diffused .OR. l_biodiffs)
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: fludif
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: fluconsol
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: fluconsol_drycell
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flu_dyninsed
   
    ! Diffusion
    REAL(KIND=rsh) :: cexcs

    ! Morphodynamic
    REAL(KIND=rlg) :: tstart_morpho   ! time beginning morphodynamic
    REAL(KIND=rlg) :: t_morpho        ! time of next morphodynamic step
    REAL(KIND=rsh) :: MF_dhsed
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: morpho0

#ifdef key_MUSTANG_V2
    REAL(KIND=rsh) :: coeff_dzsmin
    LOGICAL,DIMENSION(:,:),ALLOCATABLE :: l_isitcohesive
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: psi_sed
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: poro_mud
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: crel_mud
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: sigmapsg
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: stateconsol
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: permeab
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: E0_sand
#ifdef  key_MUSTANG_bedload
        REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE  :: flx_bx
        REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE  :: flx_by
        REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE    :: slope_dhdx
        REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE    :: slope_dhdy
        REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE    :: sedimask_h0plusxe
#if defined MORPHODYN
            INTEGER :: it_morphoYes
#endif
#endif
#ifdef key_MUSTANG_debug
        REAL(KIND=rlg)   :: t_start_debug
#endif
#endif
   

    ! Sedim output
    INTEGER :: nv_out3Dk_specif
    INTEGER :: nv_out3Dnv_specif
    INTEGER :: nv_out2D_specif
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: ep_nivsed_outp1
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: nivsed_out
    REAL(KIND=riosh), DIMENSION(:,:,:,:), ALLOCATABLE  :: var3D_cvsed
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3D_dzs
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3D_TEMP
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE :: var3D_SAL
#if defined key_BLOOM_insed
    REAL(KIND=riosh), DIMENSION(:,:,:,:), ALLOCATABLE  :: var3D_diagsed
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE  :: var2D_diagsed
#endif
#ifdef key_MUSTANG_specif_outputs
    REAL(KIND=rsh), DIMENSION(:,:,:,:), ALLOCATABLE    :: varspecif3Dk_save
    REAL(KIND=rsh), DIMENSION(:,:,:,:), ALLOCATABLE    :: varspecif3Dnv_save
    REAL(KIND=riosh), DIMENSION(:,:,:,:), ALLOCATABLE  :: varspecif3Dnv_out
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE      :: varspecif2D_save
    REAL(KIND=riosh), DIMENSION(:,:,:), ALLOCATABLE    :: varspecif2D_out
    REAL(KIND=riosh), DIMENSION(:,:,:,:), ALLOCATABLE  :: var3D_specifout
#endif


!**TODO** put under cpp key #if defined key_MUSTANG_lateralerosion
!  used in erosion only but exchange and dimensions could depend on grid model 
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_s2w_corim1
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_s2w_corip1
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_s2w_corjm1
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: flx_s2w_corjp1
!**TODO** put under cpp key #if ! defined key_nofluxwat_IWS
        REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: phieau_s2w_corim1
        REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: phieau_s2w_corip1
        REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: phieau_s2w_corjm1
        REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: phieau_s2w_corjp1
!**TODO** put under cpp key #endif
!**TODO** put under cpp key #endif

! slipdeposit : **TODO** put under cpp key key_MUSTANG_slipdeposit
   !  used in accretion (settling) only bud exchange and dimensions could depend on grid model 
   REAL(KIND=rsh),DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_corin
   REAL(KIND=rsh),DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_corim1
   REAL(KIND=rsh),DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_corip1
   REAL(KIND=rsh),DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_corjm1
   REAL(KIND=rsh),DIMENSION(:,:,:), ALLOCATABLE :: flx_w2s_corjp1


#if ! defined key_noTSdiss_insed
    ! Temperature in sediment 
    REAL(KIND=rsh), DIMENSION(:,:), ALLOCATABLE :: phitemp_s
    INTEGER, DIMENSION(:), ALLOCATABLE  :: ivdiss
    INTEGER       , DIMENSION(:), ALLOCATABLE :: D0_funcT_opt
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: D0_m0
    REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: D0_m1
#endif

#if ! defined key_nofluxwat_IWS && ! defined key_noTSdiss_insed
    REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE :: WATER_FLUX_INPUTS ! not operationnal, stil to code **TODO**
#endif

#ifdef key_BLOOM_insed
    LOGICAL :: l_out_subs_diag_sed
#endif


    CONTAINS
 
#endif /* ifdef MUSTANG */

END MODULE comMUSTANG
