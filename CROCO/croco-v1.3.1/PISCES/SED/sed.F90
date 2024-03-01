#include "cppdefs.h"

MODULE sed
   !!======================================================================
   !!                        ***  sed  ***
   !! Sediment :   set sediment global variables
   !!======================================================================
   !! History :
   !!        !  06-12  (C. Ethe)  Orignal
   !!----------------------------------------------------------------------
#if defined key_pisces

   USE par_sed
   USE par_pisces
   USE ocean2pisces
   USE oce_sed


   IMPLICIT NONE
   PUBLIC

   PUBLIC sed_alloc

#include "ocean2pisces.h90"


   !! Namelist
   REAL(wp), PUBLIC               ::  reac_sil            !: reactivity of silicate in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_clay           !: reactivity of clay in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_ligc           !: reactivity of Ligands [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_pocl           !: reactivity of pocl in  [s-1]
   REAL(wp), PUBLIC               ::  reac_pocs           !: reactivity of pocs in  [s-1]
   REAL(wp), PUBLIC               ::  reac_pocr           !: reactivity of pocr in  [s-1]
   REAL(wp), PUBLIC               ::  reac_nh4            !: reactivity of NH4 in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_h2s            !: reactivity of ODU in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_fe2            !: reactivity of Fe2+ in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_feh2s          !: reactivity of Fe2+ in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_fes            !: reactivity of Fe with H2S in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_feso           !: reactivity of FeS with O2 in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_cal            !: reactivity of cal in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  adsnh4              !: adsorption coefficient of NH4
   REAL(wp), PUBLIC               ::  ratligc             !: C/L ratio in POC
   REAL(wp), PUBLIC               ::  so2ut 
   REAL(wp), PUBLIC               ::  srno3 
   REAL(wp), PUBLIC               ::  spo4r 
   REAL(wp), PUBLIC               ::  srDnit 
   REAL(wp), PUBLIC               ::  dtsed               !: sedimentation time step
   REAL(wp), PUBLIC               ::  dtsed2              !: sedimentation time step
   INTEGER , PUBLIC               ::  nitsed000
   INTEGER , PUBLIC               ::  nitsedend
   INTEGER, PUBLIC                ::  nrseddt
   REAL    , PUBLIC               ::  sedmask
   REAL(wp), PUBLIC               ::  denssol                !: density of solid material
   INTEGER , PUBLIC               ::  numrsr, numrsw   !: logical unit for sed restart (read and write)
   LOGICAL , PUBLIC               ::  lrst_sed       !: logical to control the trc restart write
   LOGICAL , PUBLIC               ::  ln_rst_sed  = .TRUE.     !: initialisation from a restart file or not
   LOGICAL , PUBLIC               ::  ln_btbz     = .FALSE.    !: Depth variation of the bioturbation coefficient
   LOGICAL , PUBLIC               ::  ln_irrig    = .FALSE.    !: iActivation of the bioirrigation
   LOGICAL , PUBLIC               ::  ln_sed_2way = .FALSE.    !: 2 way coupling with PISCES
   LOGICAL , PUBLIC               ::  ln_sediment_offline = .FALSE. !: Offline mode for sediment module
   INTEGER             , PUBLIC   ::  nn_rstsed      !: control of the time step ( 0 or 1 ) for pass. tr.
   INTEGER , PUBLIC               ::  nn_dtsed = 1   !: frequency of step on passive tracers
   CHARACTER(len = 80) , PUBLIC   ::  cn_sedrst_in   !: suffix of pass. tracer restart name (input)
   CHARACTER(len = 256), PUBLIC   ::  cn_sedrst_indir  !: restart input directory
   CHARACTER(len = 80) , PUBLIC   ::  cn_sedrst_out  !: suffix of pass. tracer restart name (output)
   CHARACTER(len = 256), PUBLIC   ::  cn_sedrst_outdir  !: restart output directory
   CHARACTER(len = 80) , PUBLIC   ::  cn_sedwri_out  !: suffix of pass. tracer restart name (output)
   INTEGER, DIMENSION(jptrased)   ::  rstsed, wrised
   INTEGER, DIMENSION(jpdia3dsed) ::  dia3wrised
   INTEGER, DIMENSION(jpdia2dsed) ::  dia2wrised
   INTEGER, DIMENSION(5)          ::  rstadd
   INTEGER                        ::  ncidrstsed, ncidwrised, nrecsedpis, nrecsedpis_avg
   INTEGER                        ::  nwrtsedpis_avg, ntssedpis_avg 
   INTEGER                        ::  nrpfsedpis_avg, sedTsteppis_avg,sedTimepis_avg, sedTime2pis_avg
   REAL(wp)                       ::  timesedpis_avg
   LOGICAL                        ::  ldefsedpis_avg

   !
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE ::  pwcp       !: pore water sediment data at given time-step
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE ::  pwcp0      !: pore water sediment data at initial time
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE ::  pwcp_avg   !: pore water sediment data at initial time
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE ::  solcp      !: solid sediment data at given time-step
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE ::  solcp0     !: solid sediment data at initial time
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE ::  solcp_avg  !: solid sediment data at initial time
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE ::  trc_dta
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE ::  diff
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE ::  trcsed_avg
   REAL(wp), PUBLIC, DIMENSION(:,:) , ALLOCATABLE  ::  flxsed_avg

   !! * Shared module variables
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  pwcp_dta   !: pore water data at given time-step
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  rainrm_dta !: rain data at at initial time
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  rainrm     !: rain data at given time-step
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  rainrg     !: rain of each solid component in [g/(cm**2.s)]
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  fromsed    !:
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  tosed      !:
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  rloss      !: 
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  tokbot        
   !
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  temp       !: temperature
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  salt       !: salinity
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  press      !: pressure
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  raintg     !: total massic flux rained in each cell (sum of sol. comp.)
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  fecratio   !: Fe/C ratio in falling particles to the sediments
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  dzdep      !: total thickness of solid material rained [cm] in each cell
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  zkbot      !: total thickness of solid material rained [cm] in each cell
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  wacc       !: total thickness of solid material rained [cm] in each cell
   !
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  hipor      !: [h+] in mol/kg*densSW 
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  co3por     !: [co3--]solid sediment at initial time
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  dz3d       !:  ???
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  volw3d     !:  ???
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  vols3d     !:  ???


   !! Chemistry
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  densSW 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  borats 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  calcon2
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  akbs  
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  ak1s 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  ak2s   
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  akws  
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  ak12s  
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  ak1ps 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  ak2ps  
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  ak3ps 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  ak12ps 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  ak123ps
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  aksis 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  aksps 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  sieqs
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  aks3s
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  akf3s
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  sulfats
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  fluorids

   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  mol_wgt    !: molecular weight of solid sediment data
 
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE ::  trc_data    !: tracer data to share with sediment model
   !! Geometry
   INTEGER , PUBLIC, SAVE                          ::  jpoce, indoce !: Ocean points ( number/indices )
   INTEGER , PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  iarroce       !: Computation of 1D array of sediments points
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  epkbot        !: ocean bottom layer thickness
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  gdepbot       !: Depth of the sediment
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  dzkbot        !: ocean bottom layer thickness in meters
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  dz            !: sediment layers thickness
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  por           !: porosity profile     
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  por1          !: 1-por 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  volw          !: volume of pore water cell fraction
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  vols          !: volume of solid cell fraction
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  db            !: bioturbation ceofficient
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  irrig        !: bioturbation ceofficient
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  rdtsed        !:  sediment model time-step
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE :: sedligand
   REAL(wp)  ::   dens               !: density of solid material
   !! Inputs / Outputs
   CHARACTER( len = 80 ), DIMENSION(jptrased  ) ::  sedtrcl
   CHARACTER( len = 20 ), DIMENSION(jptrased  ) ::  sedtrcd , sedtrcu
   CHARACTER( len = 80 ), DIMENSION(jpdia3dsed) ::  seddia3l 
   CHARACTER( len = 20 ), DIMENSION(jpdia3dsed) ::  seddia3d, seddia3u
   CHARACTER( len = 80 ), DIMENSION(jpdia2dsed) ::  seddia2l 
   CHARACTER( len = 20 ), DIMENSION(jpdia2dsed) ::  seddia2d, seddia2u
   CHARACTER( len = 40 ), DIMENSION(5,3       ) ::  sname
   !
   INTEGER, PUBLIC ::  numsed = 27    ! units

   !! $Id: sed.F90 10425 2018-12-19 21:54:16Z smasson $
CONTAINS

   INTEGER FUNCTION sed_alloc()
      !!-------------------------------------------------------------------
      !!                    *** ROUTINE sed_alloc ***
      !!-------------------------------------------------------------------
      !!-------------------------------------------------------------------
      !
      ALLOCATE( trc_data(PRIV_2D_BIOARRAY,jpdta)                          ,   &
         &      epkbot(PRIV_2D_BIOARRAY), gdepbot(PRIV_2D_BIOARRAY)       ,   &
         &      dz(jpksed)  , por(jpksed) , por1(jpksed)                  ,   &
         &      volw(jpksed), vols(jpksed), rdtsed(jpksed)                ,   &
         &      mol_wgt(jpsol),                                           STAT=sed_alloc )

      IF( sed_alloc /= 0 )   CALL ctl_warn( 'sed_alloc: failed to allocate arrays' )
      !
   END FUNCTION sed_alloc

#endif

END MODULE sed
