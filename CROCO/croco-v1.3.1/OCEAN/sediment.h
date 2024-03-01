! $Id: sediment.h 1458 2014-02-03 15:01:25Z gcambon $
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
#ifdef SEDIMENT
/*
** Include file "sediment.h".
************************************* Meinte Blaas, John C. Warner ***
** Copyright (c) 2002/2004 Rutgers, UCLA                            **
************************************************* Hernan G. Arango ***
**                                                                  **
**  Bthk(NLAY)     User-defined initial bed layer thickness (m)     **
**  Bpor(NLAY)     User-defined initial porosity of bed layer       **
**  Bfr(NLAY,NST)  User-defined initial vol.fraction of layer, class *
**                                                                  **
**  Hrip      User-defined initial ripple height from file          **
**  Lrip      User-defined initial ripple length from file          **
**                                                                  **
**  bed_thick Sediment bed layer thickness (m)                      **
**  bed_poros Sediment bed layer porosity (void/total bed vol.)     **
**  bed_frac  Volume fraction of size class in bed layer            **
**  bed_age   Mass of sediment per layer                            **
**  bed_mass  Mass of sediment per layer                            **
**  bot_thick Active layer thickness                                **
**                                                                  **
**  settling_flux   depostion flux (kg/m2)                          **
**  ero_flux        erosion flux   (kg/m2)                          **
**  bedldu          bedload flux in xsi direction (kg/m2)           **                                                 
**  bedldv          bedload flux in eta direction (kg/m2)           **
**                                                                  **
**  morph_fac  Morphodynamic factor (non dimensionnal)              **
**                                                                  **
** Parameters for sediment model:                                   **
**                                                                  **
**  Csed     Sediment concentration (kg/m3), used during analytical **
**             initialization.                                      **
**  Erate    Surface erosion rate (kg/m2/s).                        **
**  Sd       Sediment grain diameter per size class (m).            **
**  Srho     Sediment grain density (kg/m3).                        **
**  Wsed     Particle settling velocity (m/s).                      **
**  tau_ce   Kinematic critical shear for erosion (m2/s2).          **
**  tau_cd   Kinematic critical shear for deposition (m2/s2).       **
**                                                                  **
** Sediment tracers identification indices:                         **
**                                                                  **
**  idsed    Cohesive and noncohesive sediment indices.             **
**  idmud    Cohesive sediment indices.                             **
**  isand    Noncohesive sediment indices.                          **
**                                                                  **
**  Stitle   Name of sediment.in inputfile                          **
**********************************************************************
*/  
      real Csed(NST), Erate(NST), Sd(NST), Srho(NST),
     &     Wsed(NST), tau_ce(NST), tau_cd(NST)
      common /ssediment/
     &        Csed, Erate, Sd,
     &        Srho, Wsed, tau_ce, tau_cd

      real Bthk(NLAY), Bpor(NLAY)
      common /sediment_bedthk/ Bthk, Bpor
             
      real Bfr(NLAY,NST)
      common /sediment_bedfrc/ Bfr      
      
      real Hrip, Lrip
      common /sediment_bedrip/ Hrip, Lrip
      
      real bed_thick(GLOBAL_2D_ARRAY,NLAY),
     &     bed_poros(GLOBAL_2D_ARRAY,NLAY),
     &     bed_age  (GLOBAL_2D_ARRAY,NLAY),
     &     bed_bctr (GLOBAL_2D_ARRAY,NLAY),
     &     bed_mass (GLOBAL_2D_ARRAY,NLAY,2,NST),
     &     worksed_bed(GLOBAL_2D_ARRAY,NLAY)
      common /sediment_bed/ bed_thick, bed_poros,
     &                      bed_age  , bed_mass,
     &                      bed_bctr ,
     &                      worksed_bed

      real bed_frac(GLOBAL_2D_ARRAY,NLAY,NST),
     &     worksed_frac(GLOBAL_2D_ARRAY,NLAY)
      common /sediment_frac/ bed_frac, worksed_frac     

      real bot_thick(GLOBAL_2D_ARRAY)
      common /bot_thick/ bot_thick 

      real bot_dnet(GLOBAL_2D_ARRAY)
      common /bot_dnet/ bot_dnet 

      real bot_dprp(GLOBAL_2D_ARRAY)
      common /bot_dprp/ bot_dprp 

      real bedload_coeff, morph_fac
      common /bed_coeff/ bedload_coeff, morph_fac 

      real tau_ce_2d(GLOBAL_2D_ARRAY,NST)
      common /hidexp/ tau_ce_2d

      real minlayer_thick
      common /layer_thick/ minlayer_thick
  
# ifdef SUSPLOAD
      real settling_flux(GLOBAL_2D_ARRAY,NST)
      real ero_flux(GLOBAL_2D_ARRAY,NST)
      common /sed_settling/ settling_flux,ero_flux  
# endif
# ifdef BEDLOAD
      real bedldu(GLOBAL_2D_ARRAY,NST)
      real bedldv(GLOBAL_2D_ARRAY,NST)
      common /sed_bedload/ bedldu, bedldv
#  if defined BEDLOAD_WENO5 || defined BEDLOAD_UP5
      real bedload_FX(GLOBAL_2D_ARRAY)
      real bedload_FE(GLOBAL_2D_ARRAY)
      common /sed_bedlflx/ bedload_FX,bedload_FE
#  endif
# endif
# ifdef MORPHODYN
      real bed_thick_tot(GLOBAL_2D_ARRAY,2)
      common /sed_morph/ bed_thick_tot
# endif
  
# ifdef AVERAGES
      real bed_frac_avg(GLOBAL_2D_ARRAY,NLAY,NST)
      common /sediment_frac_avg/ bed_frac_avg
      real bed_thick_avg(GLOBAL_2D_ARRAY,NLAY)
      common /sediment_thick_avg/ bed_thick_avg
      real bed_poros_avg(GLOBAL_2D_ARRAY,NLAY)
      common /sediment_poros_avg/ bed_poros_avg
#  ifdef SUSPLOAD
      real settling_flux_avg(GLOBAL_2D_ARRAY,NST)
      real ero_flux_avg(GLOBAL_2D_ARRAY,NST)
      common /sed_settling_avg/ settling_flux_avg,ero_flux_avg
#  endif
#  ifdef BEDLOAD
      real bedldu_avg(GLOBAL_2D_ARRAY,NST)
      real bedldv_avg(GLOBAL_2D_ARRAY,NST)
      common /sed_bedload_avg/ bedldu_avg, bedldv_avg
#  endif
# endif

#if defined COHESIVE_BED || defined MIXED_BED
      real tcr_min, tcr_max, tcr_slp, tcr_off, tcr_tim
      common /sed_tcr/ tcr_min, tcr_max, tcr_slp, tcr_off, tcr_tim
#endif

#if defined MIXED_BED
      real transC, transN
      common /sed_trans/ transC, transN
#endif

# if defined SED_FLOCS
      real mud_frac_eq(NMUD)
      real t_dfloc
      common /sedfloc_bed/ mud_frac_eq, t_dfloc

      logical l_ASH
      logical l_ADS
      logical l_COLLFRAG
      logical l_testcase
      integer f_ero_iv
      real f_dp0,f_alpha,f_beta,f_nb_frag
      real f_dmax,f_ater,f_clim
      real f_ero_frac,f_ero_nbfrag
      real f_nf
      real f_frag
      real f_fter
      real f_collfragparam
      common /sedflocs1/ l_ASH,l_ADS,l_COLLFRAG,
     & l_testcase,f_dp0,f_alpha,f_beta,f_nb_frag,
     & f_dmax,f_ater,f_clim,f_ero_frac,f_ero_nbfrag,
     & f_nf, f_frag,f_fter,f_collfragparam,f_ero_iv

      real rhoref
      parameter  (rhoref = 1030.0_r8)

      real f_diam(NMUD)
      real f_vol(NMUD)
      real f_rho(NMUD)
      real f_cv(NMUD)
      real f_l3(NMUD)
      real f_mass(0:NMUD+1)
      real f_coll_prob_sh(NMUD,NMUD)
      real f_coll_prob_ds(NMUD,NMUD)
      real f_l1_sh(NMUD,NMUD)
      real f_l1_ds(NMUD,NMUD)
      real f_g3(NMUD,NMUD)
      real f_l4(NMUD,NMUD)
      real f_g1_sh(NMUD,NMUD,NMUD)
      real f_g1_ds(NMUD,NMUD,NMUD)
      real f_g4(NMUD,NMUD,NMUD)    
      common /sedflocs2/ f_diam,f_vol,f_rho,f_cv,f_l3,
     &   f_mass,f_coll_prob_sh,f_coll_prob_ds,
     &   f_l1_sh,f_l1_ds,f_g3,f_l4,f_g1_sh,f_g1_ds,f_g4 
# endif

      character*80 Stitle
      common /charseds/ Stitle

#endif /* SEDIMENT */

