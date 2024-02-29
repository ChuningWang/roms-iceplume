!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org/
!
! These Sediment Dynamics routines are based on the Community Sediment 
! Transport Model (CSTM) as implemented in ROMS by J. Warner (USGS)
!
! Adapted from John Warner's code to CROCO by R. Benshila, CNRS 2015
! -- modifs: 
!   P. Marchesiello, Dec 2015
!                    Jun 2020: Van der A (2013) bedload
!   P. Marchesiello & G. Morvan, Sep 2021: 5th order interpolation
!     of bedload fluxes; DUNE test case in common with MUSTANG model
!
!=====================================================================
!
#include "cppdefs.h"
#ifdef SEDIMENT

# ifdef BEDLOAD /* for bedload scheme */
#  define bnew nnew
# else
#  define bnew nstp      
# endif

# ifdef BEDLOAD
#  ifndef BEDLOAD_VANDERA
#   define BSTRESS_UPWIND
#  endif
#  ifdef ANA_DUNE
#   undef  BSTRESS_UPWIND
#   undef  BEDLOAD_SLUMP
#   undef  BEDLOAD_SHALLOW_LIMIT
#  endif
#  ifdef BEDLOAD_VANDERA /* wave half-cycle concept */
#   define BEDLOAD_VANDERA_ZEROCURR
#   ifndef BEDLOAD_VANDERA_ZEROCURR
#     define BEDLOAD_VANDERA_FD_MADSEN
#     undef  BEDLOAD_VANDERA_FD_SANTOSS
#   endif
#   define BEDLOAD_VANDERA_ASYM_LIMITS
#   undef  BEDLOAD_VANDERA_SURFACE_WAVE
#   undef  BEDLOAD_VANDERA_WAVE_AVG_STRESS
#   define BEDLOAD_SLUMP
#   define BEDLOAD_SHALLOW_LIMIT
#  endif
# endif
!
# undef  LINEAR_CONTINUATION /* for settling scheme */
# undef  NEUMANN
!
# define MORPH_FAC *morph_fac

      subroutine sediment (tile)
!
!=======================================================================
!                                                                      !
!  This routine it is the main driver for the sediment-transport       !
!  model. Currently, it includes calls to the following routines:      !
!                                                                      !
!  * Vertical settling of sediment in the water column.                !
!  * Erosive and depositional flux interactions of sediment            !
!    between water column and the bed.                                 !
!  * Transport of multiple grain sizes.                                !
!  * Bed layer stratigraphy.                                           !
!  * Bed morphology.                                                   !
!  * Bedload based on Meyer Peter Mueller.                             !
!  * Bedload based on Soulsby combined waves + currents                !
!    (p166 Soulsby 1997)                                               !
!  * Bedload slope term options: Nemeth et al, 2006, Coastal           !
!    Engineering, v 53, p 265-275; Lesser et al, 2004, Coastal         !
!    Engineering, v 51, p 883-915.                                     !
!                                                                      !
!  * Seawater/sediment vertical level distribution:                    !
!                                                                      !
!         W-level  RHO-level                                           !
!                                                                      !
!            N     _________                                           !
!                 |         |                                          !
!                 |    N    |                                          !
!          N-1    |_________|  S                                       !
!                 |         |  E                                       !
!                 |   N-1   |  A                                       !
!            2    |_________|  W                                       !
!                 |         |  A                                       !
!                 |    2    |  T                                       !
!            1    |_________|  E                                       !
!                 |         |  R                                       !
!                 |    1    |                                          !
!            0    |_________|_____ bathymetry                          !
!                 |/////////|                                          !
!                 |    1    |                                          !
!            1    |_________|  S                                       !
!                 |         |  E                                       !
!                 |    2    |  D                                       !
!            2    |_________|  I                                       !
!                 |         |  M                                       !
!                 |  NLAY-1 |  E                                       !
!        Nbed-1   |_________|  N                                       !
!                 |         |  T                                       !
!                 |  NLAY   |                                          !
!         Nbed    |_________|                                          !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================

      implicit none
# include "param.h"
# include "ocean3d.h"
# include "scalars.h"
# include "sediment.h"
# include "private_scratch.h"
      integer tile, trd, omp_get_thread_num, itrc
# include "compute_tile_bounds.h"
!
      trd=omp_get_thread_num()
# ifdef BEDLOAD
      call sed_bedload_tile (Istr,Iend,Jstr,Jend,
     &                       A2d(1,1,trd), A2d(1,2,trd), 
     &                       A2d(1,3,trd), A2d(1,4,trd),
     &                       A2d(1,5,trd),
#  ifdef BEDLOAD_VANDERA
     &                       A2d(1,6,trd), A2d(1,7,trd), 
     &                       A2d(1,8,trd)
#  else
     &                       A2d(1,6,trd), A2d(1,7,trd) 
#  endif
     &                      )
C$OMP BARRIER
# endif /* BEDLOAD */
# ifdef SUSPLOAD
#  ifdef SED_FLOCS
      call sed_flocmod_tile (Istr,Iend,Jstr,Jend)
#  endif
C$OMP BARRIER
      call sed_settling_tile (Istr,Iend,Jstr,Jend, 
     &                        A2d(1,1,trd), A2d(1,2,trd),
     &                        A2d(1,3,trd), A2d(1,4,trd),
     &                        A2d(1,5,trd), A2d(1,6,trd),
     &                        A2d(1,7,trd), A2d(1,8,trd),
     &                        A2d(1,9,trd), A2d(1,10,trd),
     &                        B2d(1,trd) )
C$OMP BARRIER
      call sed_fluxes_tile (Istr,Iend,Jstr,Jend,A2d(1,1,trd))
C$OMP BARRIER
# endif /* SUSPLOAD */

# if defined COHESIVE_BED || defined MIXED_BED
      call sed_bed_cohesive_tile (Istr,Iend,Jstr,Jend,
     &                   A2d(1,1,trd), A2d(1,2,trd) )
# else
      call sed_bed_tile (Istr,Iend,Jstr,Jend,
     &                   A2d(1,1,trd), A2d(1,2,trd) ,
     &                   A2d(1,3,trd), A2d(1,4,trd) )
# endif
C$OMP BARRIER
!
      call sed_surf_tile (Istr,Iend,Jstr,Jend)     
C$OMP BARRIER
# if defined TRACERS
      do itrc=1,NT
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
#  ifdef THREE_GHOST_POINTS_TS
        call  exchange_r3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &          t(START_2D_ARRAY,1,nnew,itrc))
#  else
        call  exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &          t(START_2D_ARRAY,1,nnew,itrc))
#  endif
# endif
      enddo
# endif

!
      return
      end

# ifdef BEDLOAD
!======================================================================
!                                                                      !
!  This routine computes sediment bedload transport using the Meyer-   !
!  Peter and Muller (1948) formulation for unidirectional flow,        !
!  or SANTOSS formulation (Van der A 2013), that accounts for waves by !
!  treating phase-lag and wave asymetry effects in half-cycles         !
!  (from Warner et al. implementation in COAWST)                       !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Van der A, D. A., Ribberink, J. S., van der Werf, J. J.,            !
!    O’Donoghue, T., Buijsrogge, R. H., and Kranenburg, W. M., 2013:   !
!    Practical sand transport formula for non-breaking waves and       !
!    currents. Coastal Engineering, 76, 26-42                          !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!  Meyer-Peter, E. and R. Muller, 1948: Formulas for bedload transport !
!    In: Report on the 2nd Meeting International Association Hydraulic !
!    Research, Stockholm, Sweden, pp 39-64.                            !
!                                                                      !
!  Wu W. & Lin Q., 2014: Nonuniform sediment transport under non-      !
!    breaking waves and currents. Coastal Engineering, 90, pp 1-11.    !
!                                                                      !
!======================================================================
!
      subroutine sed_bedload_tile (Istr, Iend, Jstr, Jend,
     &                             tau_w,FX,FE,FX_r,FE_r
#  ifdef BEDLOAD_VANDERA
     &                            ,tau_m,phic,phicw
#  else
     &                            ,angleu,anglev
#  endif      
     &                            )
!
      implicit none
      integer Istr,Iend,Jstr,Jend
      real  FX(PRIVATE_2D_SCRATCH_ARRAY)
     &     ,FE(PRIVATE_2D_SCRATCH_ARRAY)
     &     ,FX_r(PRIVATE_2D_SCRATCH_ARRAY)
     &     ,FE_r(PRIVATE_2D_SCRATCH_ARRAY)
     &     ,tau_w(PRIVATE_2D_SCRATCH_ARRAY)
#  ifdef BEDLOAD_VANDERA
     &     ,tau_m(PRIVATE_2D_SCRATCH_ARRAY)
     &     ,phic(PRIVATE_2D_SCRATCH_ARRAY)
     &     ,phicw(PRIVATE_2D_SCRATCH_ARRAY)
#  else
     &     ,angleu(PRIVATE_2D_SCRATCH_ARRAY)
     &     ,anglev(PRIVATE_2D_SCRATCH_ARRAY)
#  endif
!
      integer i, j, k, ised
      real eps, Dstp
      real cff, cff1, cff2, cff3, cff4, cff5,cff6
      real a_slopex, a_slopey, sed_angle
      real bedld, bedld_mass, dzdx, dzdy, dzdxdy
      real smgd, smgdr, osmgd, Umag
      real Ua, Ra, phi, Clim
      real kdmax,Kbh,Kbh2,Kdh,Fwave
      real rhs_bed, bed_change
      real twopi, otwopi, sqrt2
      integer imin,imax,jmin,jmax
      parameter (eps = 1.e-20)
      parameter (kdmax = 100.)
      real K1, K2, K3, K4, K5, K6
      parameter (K1=0.6666666666, K2=0.3555555555,
     &           K3=0.1608465608, K4=0.0632098765,
     &           K5=0.0217540484, K6=0.0065407983)
#  if defined BEDLOAD_WENO5
      real :: flux5_weno
#     define FLUX5 flux5_weno
#  elif defined BEDLOAD_UP5
      real :: flux5
#     define FLUX5 flux5
#  endif
#  ifdef BEDLOAD_VANDERA
      real wavecycle, alphac, alphaw
      real ksd_wbl, ustrc_wbl, thck_wbl, udelta_wbl
      real phic_sgwbl, fd_wbl
      real Hs, Td, depth, k_wn, c_w
      real d50, d50_mix, d90, rhos
      real urms, uhat, ahat, phi_curwave
      real theta, fd, ksw, eta, alpha, tau_wRe
      real ursell_no, RR_asymwave, beta_asymwave
      real r, Su, Au, Sk, Ak, T_cu, T_tu
      real ucrest_r, utrough_r, T_crest, T_trough
      real umax, umin, uhat_c, uhat_t, mag_uc, mag_ut
      real dsf_c, dsf_t, theta_c, theta_t
      real theta_cx, theta_cy, theta_tx, theta_ty
      real mag_theta_c, mag_theta_t, mag_bstrc
      real om_cc, om_tt, om_ct, om_tc
      real bedld_cx, bedld_cy, bedld_tx, bedld_ty
      real bedld_x, bedld_y
      real kh_calc    ! functions
      real mu_calc, ksw_calc, fw_calc, fwi_calc
      real dsf_calc, w_s_calc, w_sc_calc, theta_cr_calc
#  endif
#  ifdef BEDLOAD_SLUMP
      real sedslope_crit, slopefac, slopefac_local
#  endif
#  include "param.h"
#  include "scalars.h"
#  include "grid.h"
#  include "ocean3d.h"
#  include "ocean2d.h"
#  include "sediment.h"
#  include "forces.h"
#  ifdef WKB_WWAVE
#   include "wkb_wwave.h"
#  endif
#  include "bbl.h"
#  include "mixing.h"
#  include "compute_auxiliary_bounds.h"
!
!-----------------------------------------------------------------------
! Set boundary Indices
!-----------------------------------------------------------------------
!
# ifdef EW_PERIODIC
#  define I_EXT_RANGE Istr-1,Iend+1
# else
      if (WESTERN_EDGE) then
        imin=Istr
      else
        imin=Istr-1
      endif
      if (EASTERN_EDGE) then
        imax=Iend
      else
        imax=Iend+1
      endif
#  define I_EXT_RANGE imin,imax
# endif
# ifdef NS_PERIODIC
#  define J_EXT_RANGE Jstr-1,Jend+1
# else
      if (SOUTHERN_EDGE) then
        jmin=Jstr
      else
        jmin=Jstr-1
      endif
      if (NORTHERN_EDGE) then
        jmax=Jend
      else
        jmax=Jend+1
      endif
#  define J_EXT_RANGE jmin,jmax
# endif

# if defined BEDLOAD_WENO5 || defined BEDLOAD_UP5
#  define FX_r bedload_FX
#  define FE_r bedload_FE
# endif
!
!-----------------------------------------------------------------------
! Set parameters
!-----------------------------------------------------------------------
!
      twopi=2.*pi
      otwopi=1./twopi
      sqrt2=SQRT(2.)
!
#  ifdef BEDLOAD_VANDERA
      alphac=1.0   ! bedload scale fac for curr contrib
      alphaw=1.0   ! bedload scale fac for wave contrib
#  endif
!
!-----------------------------------------------------------------------
! Compute bottom stress components at rho points
!    tau_m: form stress for flow
!    tau_w: skin stress for sediments
!-----------------------------------------------------------------------
!
      do j=J_EXT_RANGE
        do i=I_EXT_RANGE
#  if defined BEDLOAD_VANDERA && !defined BEDLOAD_VANDERA_ZEROCURR
          tau_m(i,j)=0.5*SQRT( (bustr(i,j)+bustr(i+1,j))**2 
     &                        +(bvstr(i,j)+bvstr(i,j+1))**2 ) 
#  endif
          tau_w(i,j)=SQRT( bustrw(i,j)*bustrw(i,j)
     &                   + bvstrw(i,j)*bvstrw(i,j) )
#  if !defined BEDLOAD_VANDERA && defined BSTRESS_UPWIND
#   define tau_wX angleu
#   define tau_wE anglev
          cff1=0.5*(1.0+SIGN(1.0,bustr(i+1,j)))
          cff2=0.5*(1.0-SIGN(1.0,bustr(i+1,j)))
          cff3=0.5*(1.0+SIGN(1.0,bustr(i  ,j)))
          cff4=0.5*(1.0-SIGN(1.0,bustr(i  ,j)))
          tau_wX(i,j)=cff3*(cff1*bustr(i,j)+                          
     &                cff2*0.5*(bustr(i,j)+bustr(i+1,j)))+         
     &                cff4*(cff2*bustr(i+1,j)+                        
     &                cff1*0.5*(bustr(i,j)+bustr(i+1,j)))
          cff1=0.5*(1.0+SIGN(1.0,bvstr(i,j+1)))
          cff2=0.5*(1.0-SIGN(1.0,bvstr(i,j+1)))
          cff3=0.5*(1.0+SIGN(1.0,bvstr(i,j)))
          cff4=0.5*(1.0-SIGN(1.0,bvstr(i,j)))
          tau_wE(i,j)=cff3*(cff1*bvstr(i,j)+                          
     &                cff2*0.5*(bvstr(i,j)+bvstr(i,j+1)))+         
     &                cff4*(cff2*bvstr(i,j+1)+                        
     &                cff1*0.5*(bvstr(i,j)+bvstr(i,j+1)))
#  endif
        enddo
      enddo
!
!-----------------------------------------------------------------------
! Compute angles
!-----------------------------------------------------------------------
!      
! Compute some constant bed slope parameters.
!
      sed_angle=DTAN(33.*pi/180.)    
!
!  Compute angle between currents and waves (radians).
!
      do j=J_EXT_RANGE
       do i=I_EXT_RANGE
#  ifdef BEDLOAD_VANDERA
!
! Compute angle between currents and waves
! .. from current direction toward wave vector.
!
          cff1=0.5*(bustr(i,j)+bustr(i+1,j))
          cff2=0.5*(bvstr(i,j)+bvstr(i,j+1))
          if (cff1.eq.0.) then
            phic(i,j)=0.5*pi*SIGN(1.,cff2)
          else
            phic(i,j)=ATAN2(cff2,cff1)
          endif
          phicw(i,j)=Dwave(i,j)-phic(i,j)  ! Dwave assumed in CROCO 
                                           ! grid convention
!
! Convert phicw to angle between waves and current
!  .. from wave direction towards current vector
!
          phicw(i,j)=twopi-phicw(i,j)
!
#  elif !defined BSTRESS_UPWIND
          cff1=0.5*(bustr(i,j)+bustr(i+1,j))
          cff2=0.5*(bvstr(i,j)+bvstr(i,j+1))
          Umag=SQRT(cff1*cff1+cff2*cff2)+eps
          angleu(i,j)=cff1/Umag
          anglev(i,j)=cff2/Umag
#  endif /* BEDLOAD_VANDERA */
        enddo
      enddo
!
      do ised=1,NST   ! <============== sed loop
!
!-----------------------------------------------------------------------
! Compute magnitude of bedload at rho points
!
!  bedld has dimensions of kg m-1 s-1. 
!  FX_r and FE_r have dimensions of kg.
!-----------------------------------------------------------------------
!
#  ifdef BEDLOAD_MARIEU
        smgdr=Srho(ised)
#  else
        smgd=(Srho(ised)/rho0-1.)*g*Sd(ised)
        osmgd=1./smgd 
        smgdr=SQRT(smgd)*Sd(ised)*Srho(ised)
#  endif
!
#  ifdef BEDLOAD_VANDERA
        rhos=Srho(ised)          ! (kg/m3)
        d50=Sd(ised)             ! (m)
        d90=1.3*d50              ! (m)
#  endif
!
        do j=J_EXT_RANGE     ! <===== i,j loops
          do i=I_EXT_RANGE
!
#  ifdef BEDLOAD_VANDERA
!
! -----------------  Van der A formulation  ---------------------
!
! VA-2013 equation 1 is solved:
! The transport formula is based on the wave “half-cycle” concept: 
! The wave-averaged net sand transport rate in the wave boundary 
! layer is described as the difference between sand transported 
! during the positive “crest” and the negative “trough” half-cycles. 
! Phase-lag and wave asymetry effects are taken into account.
!
! Wave & sediment parameters
!
            depth=z_w(i,j,N)+h(i,j)             ! (m)
            Hs=2.*Awave(i,j)                    ! (m)
            Td=MAX(Pwave(i,j),1.)               ! (s)
            phi_curwave=phicw(i,j)              ! angle cw
            k_wn=kh_calc(Td,depth)/depth        ! Wave number 
            c_w=twopi/(k_wn*Td)                 ! Wave speed

            if(NST>1) then 
              d50_mix=Sdens(i,j)
            else
             d50_mix=d50
            endif
!
! Compute bed wave orbital velocity (m/s) & excursion amplitude
!
            urms=twopi/Td*Awave(i,j)/(SINH(k_wn*depth)+eps)
#   ifdef MASKING
     &                                          *rmask(i,j)
#   endif
#   ifdef WET_DRY
     &                                      *rmask_wet(i,j)
#   endif
            uhat=urms*sqrt2
            ahat=uhat*Td*otwopi
! 
! Magnitude of stress for computing current velocity 
! at the wave boundary layer
!
#   ifndef BEDLOAD_VANDERA_ZEROCURR
            mag_bstrc=tau_m(i,j)
#    ifdef MASKING
     &                             *rmask(i,j)
#    endif
#    ifdef WET_DRY
     &                         *rmask_wet(i,j)
#    endif
#   else
            mag_bstrc=0.
#   endif
!
! Compute near bottom current at wave boundary layer (wbl)
! Use wbl thickness based on Madsen to get near bottom current velocity. 
! Make sure that wbl is under total depth and greater than apparent 
! roughness.
!
#   ifdef BEDLOAD_VANDERA_ZEROCURR
            udelta_wbl=0.
            ksd_wbl=0.
            thck_wbl=0.
            ustrc_wbl=0.
#   else
#    ifdef BBL
            ksd_wbl=30.*Zbapp(i,j)
#    else
            ksd_wbl=30.*Zob(i,j)
#    endif
            thck_wbl=0.09*ksd_wbl*(ahat/ksd_wbl)**0.82
            !thck_wbl=0.1
            ustrc_wbl=sqrt(mag_bstrc)
            cff1=MAX(MIN(0.9*depth, thck_wbl), 1.1*ksd_wbl)
            cff2=LOG(30.*cff1/ksd_wbl)
            udelta_wbl=(ustrc_wbl/vonKar)*cff2
#   endif
!
! Ruessink et al. provides equations for calculating skewness parameters
! Uses Malarkey and Davies equations to get "r" and "phi"
! common to both crest and trough.
!
            call skewness_params(Hs, Td, depth, r, phi, ursell_no)
!        
! Abreu et al. use skewness params (r, phi) to get representative critical 
! orbital velocity for crest and trough cycles.
! 
            call abreu_points(r, phi, uhat, Td,
     &                        T_crest, T_trough,
     &                        T_cu, T_tu, umax, umin,
     &                        RR_asymwave, beta_asymwave)
!
! Representative crest half cycle water particle velocity
! as well as full cycle orbital velocity and excursion.
! VA2013 Equation 10, 11.
!
            uhat_c=umax
            uhat_t=umin
            ucrest_r=0.5*sqrt2*uhat_c
            utrough_r=0.5*sqrt2*uhat_t
! 
! Full wave cycle 
! 
            call full_wave_cycle_stress_factors(d50, d90, osmgd,
     &                                           Td, c_w, depth,
     &                                              phi_curwave,
     &                                  RR_asymwave, uhat, ahat,
     &                                    umax, umin, mag_bstrc,
     &                            T_crest, T_trough, T_cu, T_tu,
     &                                       ksd_wbl, ustrc_wbl,
     &                             thck_wbl, udelta_wbl, fd_wbl,
     &                                 alpha, eta, ksw, tau_wRe )
!
! Bed shear stress (Shields parameter) for Crest half cycle 
! alpha VA2013 Eqn. 19  
!
            call half_wave_cycle_stress_factors( T_cu, T_crest,
     &                                               ahat, ksw,
     &                                           fd_wbl, alpha,
     &                                          alphac, alphaw,
     &                                              d50, osmgd,
     &               ucrest_r, uhat_c, udelta_wbl, phi_curwave,
     &                                         tau_wRe,  dsf_c,
     &                         theta_cx, theta_cy, mag_theta_c )
!
! Compute sediment load entrained during each crest half cycle
!
            wavecycle=1.
            call sandload_vandera( wavecycle,
     &                             Hs, Td,  depth, RR_asymwave,
     &                                 d50, d50_mix, rhos, c_w,
     &                                              eta, dsf_c,
     &                                   T_crest, T_cu, uhat_c,
     &                               mag_theta_c, om_cc, om_ct )
!
! Bed shear stress (Shields parameter) for Trough half cycle 
! alpha VA2013 Eqn. 19  
!
            call half_wave_cycle_stress_factors( T_tu, T_trough, 
     &                                                ahat, ksw,
     &                                            fd_wbl, alpha,
     &                                           alphac, alphaw,
     &                                               d50, osmgd,
     &               utrough_r, uhat_t, udelta_wbl, phi_curwave,
     &                                           tau_wRe, dsf_t,
     &                          theta_tx, theta_ty, mag_theta_t )
!
! Compute sediment load entrained during each trough half cycle
!
            wavecycle=-1.
            call sandload_vandera( wavecycle,
     &                              Hs, Td,  depth, RR_asymwave,
     &                                  d50, d50_mix, rhos, c_w,
     &                                               eta, dsf_t,
     &                                   T_trough, T_tu, uhat_t,
     &                                mag_theta_t, om_tt, om_tc )
!
! Compute sediment load entrained during crest and trough half cycles
!
            cff1=MAX(0.5*T_crest/T_cu, 0.)
            cff2=sqrt(mag_theta_c)*(om_cc+cff1*om_tc)  
            cff3=(theta_cx/mag_theta_c)
            bedld_cx=cff2*cff3
            cff3=(theta_cy/mag_theta_c)
            bedld_cy=cff2*cff3
!
            cff1=MAX(0.5*T_trough/T_tu, 0.)
            cff2=sqrt(mag_theta_t)*(om_tt+cff1*om_ct)
            cff3=(theta_tx/mag_theta_t)
            bedld_tx=cff2*cff3
            cff3=(theta_ty/mag_theta_t)
            bedld_ty=cff2*cff3
!
! Total bedload (kg m-1 sec-1). Use velocity-load equation 1 of VA2013 
!
            bedld_x=smgdr*( bedld_cx*T_crest+
     &                      bedld_tx*T_trough )/Td
            bedld_y=smgdr*( bedld_cy*T_crest+
     &                      bedld_ty*T_trough )/Td
!
! Partition bedload into xi and eta directions, still at rho points
! (FX_r and FE_r have dimensions of kg).
! --> project on CROCO grid (bedld is wave-aligned in Van Der A)
!
            FX_r(i,j)=(bedld_x*COS(Dwave(i,j)) -
     &                 bedld_y*SIN(Dwave(i,j)))*on_r(i,j)*dt
            FE_r(i,j)=(bedld_x*SIN(Dwave(i,j)) +
     &                 bedld_y*COS(Dwave(i,j)))*om_r(i,j)*dt
!
! -----------------  END Van der A formulation  ------------------
!
#  else
!
!                    Use simpler formulations
!
! -------------  Meyer-Peter Muller or Wu & Lin  -----------------
!
#   ifdef BSTRESS_UPWIND
!
! Use partitions of stress from upwind direction, still at rho points.
!
#    ifdef BEDLOAD_MPM
            bedld=8.*(MAX((ABS(tau_wX(i,j))*osmgd-0.047),0.)**1.5)
#    elif defined BEDLOAD_WULIN
            cff=ABS(tau_wX(i,j))/tau_ce_2d(i,j,ised)-1.
            bedld=0.0053*(MAX(cff,0.)**2.2)
#    elif defined BEDLOAD_MARIEU
            bedld=(0.001*ubar(i,j,knew)**3)
#    endif
     &                                 *smgdr*SIGN(1.,tau_wX(i,j))
            FX_r(i,j)=bedld*on_r(i,j)*dt
#    ifdef BEDLOAD_MPM
            bedld=8.*(MAX((ABS(tau_wE(i,j))*osmgd-0.047),0.)**1.5)
#    elif defined BEDLOAD_WULIN
            cff=ABS(tau_wE(i,j))/tau_ce_2d(i,j,ised)-1
            bedld=0.0053*(MAX(cff,0.)**2.2)
#    elif defined BEDLOAD_MARIEU
            bedld=0.
#    endif
     &                                 *smgdr*SIGN(1.,tau_wE(i,j))
            FE_r(i,j)=bedld*om_r(i,j)*dt
#    undef tau_wX
#    undef tau_wE
#   else
!
! Use mean stress tau_w at rho points
!
#    ifdef BEDLOAD_MPM
            bedld=8.*(MAX((tau_w(i,j)*osmgd-0.047),0.)**1.5)*smgdr
#    elif defined BEDLOAD_WULIN
            cff=tau_w(i,j)/tau_ce_2d(i,j,ised)-1.
            bedld=0.0053*(MAX(cff,0.)**2.2)*smgdr
#    elif defined BEDLOAD_MARIEU
            bedld=(0.001*ubar(i,j,knew)**3)*smgdr
#    endif
!
! Partition bedld into xi and eta directions, still at rho points.
! (FX_r and FE_r have dimensions of kg).
!
            FX_r(i,j)=angleu(i,j)*bedld*on_r(i,j)*dt
            FE_r(i,j)=anglev(i,j)*bedld*om_r(i,j)*dt
#   endif /* BSTRESS_UPWIND */

#  endif /* BEDLOAD_VANDERA */
!
!-----------------------------------------------------------------------
! Correct for along-direction slope. 
! Limit slope to 0.9*sed angle.
!-----------------------------------------------------------------------
!
#  if defined SLOPE_NEMETH || defined SLOPE_LESSER
            cff1=0.5*(1.+SIGN(1.,FX_r(i,j)))
            cff2=0.5*(1.-SIGN(1.,FX_r(i,j)))
            cff3=0.5*(1.+SIGN(1.,FE_r(i,j)))
            cff4=0.5*(1.-SIGN(1.,FE_r(i,j)))
#  endif
#  ifdef SLOPE_NEMETH
            dzdx=(h(i+1,j)-h(i  ,j))/om_u(i+1,j)*cff1+       
     &           (h(i-1,j)-h(i  ,j))/om_u(i  ,j)*cff2
            dzdy=(h(i,j+1)-h(i,j  ))/on_v(i,j+1)*cff3+      
     &           (h(i,j-1)-h(i,j  ))/on_v(i  ,j)*cff4

            a_slopex=1.7*dzdx
            a_slopey=1.7*dzdy
!
! Add contribution of bed slope to bed load transport fluxes.
!
            FX_r(i,j)=FX_r(i,j)*(1.+a_slopex)
            FE_r(i,j)=FE_r(i,j)*(1.+a_slopey)
!
#  elif defined SLOPE_LESSER
            dzdx=MIN(((h(i+1,j)-h(i  ,j))/om_u(i+1,j)*cff1+      
     &                (h(i  ,j)-h(i-1,j))/om_u(i  ,j)*cff2),0.52)* 
     &                SIGN(1.,FX_r(i,j))
            dzdy=MIN(((h(i,j+1)-h(i,j  ))/on_v(i,j+1)*cff3+      
     &                (h(i,j  )-h(i,j-1))/on_v(i  ,j)*cff4),0.52)* 
     &                SIGN(1.,FE_r(i,j))
            cff=ATAN(dzdx)
            a_slopex=sed_angle/(COS(cff)*(sed_angle-dzdx))
            cff=ATAN(dzdy)
            a_slopey=sed_angle/(COS(cff)*(sed_angle-dzdy))
!
! Add contribution of bed slope to bed load transport fluxes.
!
            FX_r(i,j)=FX_r(i,j)*a_slopex
            FE_r(i,j)=FE_r(i,j)*a_slopey
!
#  endif /* SLOPE_NEMETH || SLOPE_LESSER */
!
!-----------------------------------------------------------------------
! Apply morphology factor.
!-----------------------------------------------------------------------
!
            FX_r(i,j)=FX_r(i,j) MORPH_FAC
            FE_r(i,j)=FE_r(i,j) MORPH_FAC
!
!-----------------------------------------------------------------------
! Apply bedload transport rate coefficient. 
! Also limit bedload to the fraction of each sediment class.
!-----------------------------------------------------------------------
!
            FX_r(i,j)=FX_r(i,j)*bedload_coeff*bed_frac(i,j,1,ised)
            FE_r(i,j)=FE_r(i,j)*bedload_coeff*bed_frac(i,j,1,ised)
!
!-----------------------------------------------------------------------
! Limit bed load to available bed mass.
!-----------------------------------------------------------------------
!
            bedld_mass=ABS(FX_r(i,j))+ABS(FE_r(i,j))+eps
            FX_r(i,j)=MIN(ABS(FX_r(i,j)),                      
     &                    bed_mass(i,j,1,nstp,ised)*           
     &                    om_r(i,j)*on_r(i,j)*ABS(FX_r(i,j))/    
     &                    bedld_mass)*                        
     &                    SIGN(1.,FX_r(i,j))
            FE_r(i,j)=MIN(ABS(FE_r(i,j)),                    
     &                    bed_mass(i,j,1,nstp,ised)*          
     &                    om_r(i,j)*on_r(i,j)*ABS(FE_r(i,j))/     
     &                    bedld_mass)*             
     &                    SIGN(1.,FE_r(i,j))
#  ifdef MASKING
            FX_r(i,j)=FX_r(i,j)*rmask(i,j)
            FE_r(i,j)=FE_r(i,j)*rmask(i,j)
#  endif
#  ifdef WET_DRY
            FX_r(i,j)=FX_r(i,j)*rmask_wet(i,j)
            FE_r(i,j)=FE_r(i,j)*rmask_wet(i,j)
#  endif
!
          enddo ! I_EXT_RANGE
        enddo  ! J_EXT_RANGE  <===== i,j loops
!
!-----------------------------------------------------------------------
! Apply boundary conditions
!-----------------------------------------------------------------------
!
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=J_EXT_RANGE
            FX_r(Istr-1,j)=FX_r(Istr,j)
            FE_r(Istr-1,j)=FE_r(Istr,j)
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=J_EXT_RANGE
            FX_r(Iend+1,j)=FX_r(Iend,j)
            FE_r(Iend+1,j)=FE_r(Iend,j)
          enddo
        endif
#  endif
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=I_EXT_RANGE
            FX_r(i,Jstr-1)=FX_r(i,Jstr)
            FE_r(i,Jstr-1)=FE_r(i,Jstr)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=I_EXT_RANGE
            FX_r(i,Jend+1)=FX_r(i,Jend)
            FE_r(i,Jend+1)=FE_r(i,Jend)
          enddo
        endif
#  endif

#  if !defined EW_PERIODIC && !defined NS_PERIODIC
        if ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) then
          FX_r(Istr-1,Jstr-1)=0.5*(FX_r(Istr  ,Jstr-1)+     
     &                             FX_r(Istr-1,Jstr  ))
          FE_r(Istr-1,Jstr-1)=0.5*(FE_r(Istr  ,Jstr-1)+       
     &                             FE_r(Istr-1,Jstr  ))
        endif
        if ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) then
          FX_r(Iend+1,Jstr-1)=0.5*(FX_r(Iend  ,Jstr-1)+      
     &                             FX_r(Iend+1,Jstr  ))
          FE_r(Iend+1,Jstr-1)=0.5*(FE_r(Iend  ,Jstr-1)+       
     &                             FE_r(Iend+1,Jstr  ))
        endif
        if ((NORTHERN_EDGE).and.(WESTERN_EDGE)) then
          FX_r(Istr-1,Jend+1)=0.5*(FX_r(Istr-1,Jend  )+       
     &                             FX_r(Istr  ,Jend+1))
          FE_r(Istr-1,Jend+1)=0.5*(FE_r(Istr-1,Jend  )+      
     &                             FE_r(Istr  ,Jend+1))
        endif
        if ((NORTHERN_EDGE).and.(EASTERN_EDGE)) then
          FX_r(Iend+1,Jend+1)=0.5*(FX_r(Iend+1,Jend  )+      
     &                             FX_r(Iend  ,Jend+1))
          FE_r(Iend+1,Jend+1)=0.5*(FE_r(Iend+1,Jend  )+    
     &                             FE_r(Iend  ,Jend+1))
        endif
#  endif

#  undef I_EXT_RANGE
#  undef J_EXT_RANGE
!
!-----------------------------------------------------------------------
! Compute face fluxes at u and v points before taking divergence.
!
!  --> Extended stencil in case of high-order interpolation.
!-----------------------------------------------------------------------
!
#  if defined BEDLOAD_WENO5 || defined BEDLOAD_UP5
C$OMP BARRIER
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
        if (IstrU.le.Iend) then
          call exchange_r2d_3pts_tile (Istr,Iend,Jstr,Jend,FX_r)
        endif
        if (JstrV.le.Jend) then
          call exchange_r2d_3pts_tile (Istr,Iend,Jstr,Jend,FE_r)
        endif
#   endif

#   ifdef EW_PERIODIC
        imin=1
        imax=LOCALLM+1
#   else
#    ifdef MPI
        if (WESTERN_EDGE) then
          imin=3
        else
          imin=1
        endif
        if (EASTERN_EDGE) then
          imax=Lmmpi-1
        else
          imax=Lmmpi+1
        endif
#    else
        imin=3
        imax=Lm-1 
#    endif
#   endif
#   ifdef NS_PERIODIC
        jmin=1
        jmax=LOCALMM+1
#   else
#    ifdef MPI
        if (SOUTHERN_EDGE) then
          jmin=3
        else
          jmin=1
        endif
        if (NORTHERN_EDGE) then
          jmax=Mmmpi-1
        else
          jmax=Mmmpi+1
        endif
#    else
        jmin=3
        jmax=Mm-1
#    endif
#   endif
#  endif /* BEDLOAD_WENO5 || BEDLOAD_UP5 */
!
        do j=Jstr,Jend
          do i=Istr,Iend+1
#  if defined BEDLOAD_WENO5 || defined BEDLOAD_UP5
            cff=SIGN(1.,FX_r(i,j))
            cff1=0.5*(1.+SIGN(1.,FX_r(i,j)))
            cff2=0.5*(1.-SIGN(1.,FX_r(i,j)))
            IF ( (( i.lt.imin ) .or. ( i.gt.imax )) 
     &                          .and.( imin.lt.imax ) ) THEN
               FX(i,j)=0.5*(1.+SIGN(1.,FX_r(i-1,j)))*
     &              (cff1*FX_r(i-1,j)+
     &               cff2*0.5*(FX_r(i-1,j)+FX_r(i,j)))+
     &                 0.5*(1.-SIGN(1.,FX_r(i-1,j)))*
     &              (cff2*FX_r(i  ,j)+
     &               cff1*0.5*(FX_r(i-1,j)+FX_r(i,j)))
            ELSE
               FX(i,j)=FLUX5(FX_r(i-3,j),FX_r(i-2,j),
     &                       FX_r(i-1,j),FX_r(i  ,j),
     &                       FX_r(i+1,j),FX_r(i+2,j),cff)
            ENDIF
#  elif defined BEDLOAD_UP1
            cff1=0.5*(1.+SIGN(1.,FX_r(i,j)))
            cff2=0.5*(1.-SIGN(1.,FX_r(i,j)))
            FX(i,j)=0.5*(1.+SIGN(1.,FX_r(i-1,j)))*     
     &              (cff1*FX_r(i-1,j)+                    
     &               cff2*0.5*(FX_r(i-1,j)+FX_r(i,j)))+   
     &              0.5*(1.-SIGN(1.,FX_r(i-1,j)))*    
     &              (cff2*FX_r(i  ,j)+                   
     &               cff1*0.5*(FX_r(i-1,j)+FX_r(i,j)))
#  else
            FX(i,j)=0.5*(FX_r(i-1,j)+FX_r(i,j))
#  endif
!
#  ifdef SLOPE_KIRWAN
            dzdx=(h(i,j)-h(i-1,j))/om_u(i,j)
            a_slopex=(MAX(0.,abs(dzdx)-0.01)
     &               *SIGN(1.,dzdx)*0.1)
     &               *om_r(i,j)*dt
            FX(i,j)=FX(i,j)+a_slopex MORPH_FAC
#  endif /* SLOPE_KIRWAN */
#  ifdef MASKING
            FX(i,j)=FX(i,j)*umask(i,j)
#  endif
#  ifdef WET_DRY
            FX(i,j)=FX(i,j)*umask_wet(i,j)
#  endif
          enddo
        enddo
        do j=Jstr,Jend+1
          do i=Istr,Iend
#  if defined BEDLOAD_WENO5 || defined BEDLOAD_UP5
            cff=SIGN(1.,FE_r(i,j))
            cff1=0.5*(1.+SIGN(1.,FE_r(i,j)))
            cff2=0.5*(1.-SIGN(1.,FE_r(i,j)))
            IF ( (( j.lt.jmin ) .or. ( j.gt.jmax ))
     &                          .and.( jmin.lt.jmax ) ) THEN
               FE(i,j)=0.5*(1.+SIGN(1.,FE_r(i,j-1)))*
     &              (cff1*FE_r(i,j-1)+
     &               cff2*0.5*(FE_r(i,j-1)+FE_r(i,j)))+
     &                 0.5*(1.-SIGN(1.,FE_r(i,j-1)))*
     &              (cff2*FE_r(i  ,j)+
     &               cff1*0.5*(FE_r(i,j-1)+FE_r(i,j)))
            ELSE
               FE(i,j)=FLUX5(
     &                      FE_r(i,j-3),FE_r(i,j-2),
     &                      FE_r(i,j-1),FE_r(i  ,j),
     &                      FE_r(i,j+1),FE_r(i,j+2),cff)
            ENDIF
#  elif defined BEDLOAD_UP1
            cff1=0.5*(1.+SIGN(1.,FE_r(i,j)))
            cff2=0.5*(1.-SIGN(1.,FE_r(i,j)))
            FE(i,j)=0.5*(1.+SIGN(1.,FE_r(i,j-1)))*     
     &              (cff1*FE_r(i,j-1)+                      
     &               cff2*0.5*(FE_r(i,j-1)+FE_r(i,j)))+       
     &              0.5*(1.-SIGN(1.,FE_r(i,j-1)))*      
     &              (cff2*FE_r(i  ,j)+                    
     &               cff1*0.5*(FE_r(i,j-1)+FE_r(i,j)))
#  else
            FE(i,j)=0.5*(FE_r(i,j-1)+FX_r(i,j))
#  endif
#  ifdef SLOPE_KIRWAN
            dzdy=(h(i,j)-h(i,j-1))/on_v(i,j)
            a_slopey=(MAX(0.,abs(dzdy)-0.01)
     &               *SIGN(1.,dzdy)*0.1)
     &               *on_r(i,j)*dt
            FE(i,j)=FE(i,j)+a_slopey MORPH_FAC
#  endif /* SLOPE_KIRWAN */
#  ifdef MASKING
            FE(i,j)=FE(i,j)*vmask(i,j)
#  endif
#  ifdef WET_DRY
            FE(i,j)=FE(i,j)*vmask_wet(i,j)
#  endif
          enddo
        enddo

#  if defined BEDLOAD_WENO5 || defined BEDLOAD_UP5
#   undef FX_r
#   undef FE_r
#  endif
!
#  ifdef BEDLOAD_SLUMP
!-----------------------------------------------------------------------
!  Sed slump computation to allow slumping (from COAWST).
!     sedslope is the critical slope to allow slumping.
!     slopefac is the scale factor for sediment movement.
!-----------------------------------------------------------------------
!
        sedslope_crit=0.1  ! Critical bed slope for slumping
        slopefac=0.1       ! Bedload bed slumping factor
!
!  U-direction slumping
!
        do j=Jstr,Jend
          do i=Istr,Iend+1
            dzdx=(h(i,j)-h(i-1,j))/om_u(i,j)
            dzdy=(h(i,j)-h(i,j-1))/on_v(i,j)
            dzdxdy=sqrt(dzdx**2.+dzdy**2.)
            cff=Srho(ised)*(1.-bed_poros(i,j,1))
            slopefac_local=cff*dt*slopefac MORPH_FAC
            cff=ABS(dzdxdy)-sedslope_crit
            cff1=(0.5+SIGN(0.5,cff))
            cff=0.5*cff*cff1/(pm(i,j)*pn(i,j))
            cff2=cff*slopefac_local*SIGN(1.,dzdx)
            FX(i,j)=FX(i,j)+cff2
          end do
        end do
!
!  V-direction slumping
!
        do j=Jstr,Jend+1
          do i=Istr,Iend
            dzdy=(h(i,j)-h(i,j-1))/on_v(i,j) ! positive is downhill
            dzdx=(h(i,j)-h(i-1,j))/om_u(i,j)
            dzdxdy=sqrt(dzdx**2.+dzdy**2.)
            cff=Srho(ised)*(1.-bed_poros(i,j,1))
            slopefac_local=cff*dt*slopefac MORPH_FAC
            cff=ABS(dzdxdy)-sedslope_crit
            cff1=(0.5+SIGN(0.5,cff))
            cff=0.5*cff*cff1/(pm(i,j)*pn(i,j))
            cff2=cff*slopefac_local*SIGN(1.,dzdy)
            FE(i,j)=FE(i,j)+cff2
          end do
        end do
#  endif /* BEDLOAD_SLUMP */
!
!-----------------------------------------------------------------------
! Limit fluxes to prevent bottom from breaking thru water surface.
!-----------------------------------------------------------------------
!
        do j=Jstr,Jend
          do i=Istr,Iend
!
! Total thickness available.
!
#  ifdef WET_DRY
            Dstp=MAX(z_w(i,j,N)-z_w(i,j,0)-1.5*Dcrit(i,j),0.)
#  else
            Dstp=MAX(z_w(i,j,N)-z_w(i,j,0)-0.01,0.)  !!!
#  endif
!
! Thickness change that wants to occur.
!
            rhs_bed=(FX(i+1,j)-FX(i,j)+                          
     &               FE(i,j+1)-FE(i,j))*pm(i,j)*pn(i,j)
            bed_change=rhs_bed/(Srho(ised)*(1.-bed_poros(i,j,1)))
!
! Limit that change to be less than available (if accretional)
!
            cff=MAX(bed_change-Dstp,0.)
            cff1=cff/ABS(bed_change+eps)
            cff1=cff1*0.5*(1.-SIGN(1.,bed_change))
            FX(i+1,j  )=FX(i+1,j  )*(1.-cff1)
            FX(i  ,j  )=FX(i  ,j  )*(1.-cff1)
            FE(i  ,j+1)=FE(i  ,j+1)*(1.-cff1)
            FE(i  ,j  )=FE(i  ,j  )*(1.-cff1)
          end do
        end do
!
!-----------------------------------------------------------------------
! Apply closed boundary conditions on (FX,FE) if any
!-----------------------------------------------------------------------
!
#  ifndef EW_PERIODIC
#   ifndef OBC_WEST
        if (WESTERN_EDGE) then
          do j=Jstr-1,Jend+1
            FX(Istr,j)=0.
          enddo
        endif
#   endif
#   ifndef OBC_EAST
        if (EASTERN_EDGE) then
          do j=Jstr-1,Jend+1
            FX(Iend+1,j)=0.
          enddo
        endif
#   endif
#  endif
#  ifndef NS_PERIODIC
#   ifndef OBC_SOUTH
        if (SOUTHERN_EDGE) then
          do i=Istr-1,Iend+1
            FE(i,Jstr)=0.
          enddo
        endif
#   endif
#   ifndef OBC_NORTH
        if (NORTHERN_EDGE) then
          do i=Istr-1,Iend+1
            FE(i,Jend+1)=0.
          enddo
        endif
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Determine flux divergence and evaluate change in bed properties.
!-----------------------------------------------------------------------
!
        do j=Jstr,Jend
          do i=Istr,Iend
            cff=(FX(i+1,j)-FX(i,j)+                           
     &           FE(i,j+1)-FE(i,j))*pm(i,j)*pn(i,j)
#  ifdef SPONGE_SED
     &                             *sed_sponge(i,j)
#  endif

#  ifdef BEDLOAD_SHALLOW_LIMIT
            cff=cff/max(1.,(0.5/(z_w(i,j,N)-z_w(i,j,0)
#   ifdef WET_DRY
     &                           -Dcrit(i,j)+eps
#   endif
     &            )**2))
     &       
#  endif

            bed_mass(i,j,1,nnew,ised)=MAX(bed_mass(i,j,1,nstp,ised)-   
     &                                    cff,0.)
#  if !defined SUSPLOAD
            if (NLAY.gt.1) then
              do k=2,NLAY
                bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k,nstp,ised)
              enddo
            endif
#  endif
            bed_thick(i,j,1)=MAX(bed_thick(i,j,1)-            
     &                           cff/(Srho(ised)*              
     &                           (1.-bed_poros(i,j,1))), 0.)
          enddo
        enddo
!
!-----------------------------------------------------------------------
!  Output bedload transports [kg/m2/s].
!-----------------------------------------------------------------------
!
        cff=1./dt
        do j=Jstr,Jend
          do i=Istr,Iend
            bedldu(i,j,ised)=(FX(i,j)-FX(i+1,j))*pm(i,j)*pn(i,j)*cff
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            bedldv(i,j,ised)=(FE(i,j)-FE(i,j+1))*pm(i,j)*pn(i,j)*cff
          enddo
        enddo
!
      enddo  ! ised =========
!
!-----------------------------------------------------------------------
!  Update mean surface properties.
!-----------------------------------------------------------------------
!
      do j=Jstr,Jend
        do i=Istr,Iend
          cff3=0.0
          do ised=1,NST
            cff3=cff3+bed_mass(i,j,1,nnew,ised)
          end do
          if (cff3.eq.0.0) then
            cff3=eps
          end if
          do ised=1,NST
            bed_frac(i,j,1,ised)=bed_mass(i,j,1,nnew,ised)/cff3
          end do
!
          cff1=1.
          cff2=1.
          cff3=1.
          cff4=1.
          do ised=1,NST
            cff1=cff1*tau_ce_2d(i,j,ised)**bed_frac(i,j,1,ised)
            cff2=cff2*Sd(ised)**bed_frac(i,j,1,ised)
            cff3=cff3*(wsed(ised)+eps)**bed_frac(i,j,1,ised)
            cff4=cff4*Srho(ised)**bed_frac(i,j,1,ised)
          enddo
          taucb(i,j)=cff1
          Ssize(i,j)=cff2
          w_set(i,j)=cff3
          Sdens(i,j)=MAX(cff4,1050.)
        enddo
      enddo
!
! Set boundary conditions
!
# ifndef EW_PERIODIC
      if (EASTERN_EDGE) then
        do j=Jstr,Jend
           taucb(Iend+1,j)=taucb(Iend,j)
           Ssize(Iend+1,j)=Ssize(Iend,j)
           w_set(Iend+1,j)=w_set(Iend,j)
           Sdens(Iend+1,j)=Sdens(Iend,j)
        enddo
      endif
      if (WESTERN_EDGE) then
        do j=Jstr,Jend
           taucb(Istr-1,j)=taucb(Istr,j)
           Ssize(Istr-1,j)=Ssize(Istr,j)
           w_set(Istr-1,j)=w_set(Istr,j)
           Sdens(Istr-1,j)=Sdens(Istr,j)
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (NORTHERN_EDGE) then
        do i=Istr,Iend
          taucb  (i,Jend+1)=taucb(i,Jend)
          Ssize  (i,Jend+1)=Ssize(i,Jend)
          w_set  (i,Jend+1)=w_set(i,Jend)
          Sdens  (i,Jend+1)=Sdens(i,Jend)
        enddo
      endif
      if (SOUTHERN_EDGE) then
        do i=Istr,Iend
          taucb  (i,Jstr-1)=taucb(i,Jstr)
          Ssize  (i,Jstr-1)=Ssize(i,Jstr)
          w_set  (i,Jstr-1)=w_set(i,Jstr)
          Sdens  (i,Jstr-1)=Sdens(i,Jstr)
        enddo
      endif
# endif
# if !defined EW_PERIODIC && !defined NS_PERIODIC
      if (SOUTHERN_EDGE.and.WESTERN_EDGE) then
        taucb(Istr-1,Jstr-1)=0.5*(taucb(Istr,Jstr-1)+taucb(Istr-1,Jstr))
        Ssize(Istr-1,Jstr-1)=0.5*(Ssize(Istr,Jstr-1)+Ssize(Istr-1,Jstr))
        w_set(Istr-1,Jstr-1)=0.5*(w_set(Istr,Jstr-1)+w_set(Istr-1,Jstr))
        Sdens(Istr-1,Jstr-1)=0.5*(Sdens(Istr,Jstr-1)+Sdens(Istr-1,Jstr))
      endif
      if (SOUTHERN_EDGE.and.EASTERN_EDGE) then
        taucb(Iend+1,Jstr-1)=0.5*(taucb(Iend,Jstr-1)+taucb(Iend+1,Jstr))
        Ssize(Iend+1,Jstr-1)=0.5*(Ssize(Iend,Jstr-1)+Ssize(Iend+1,Jstr))
        w_set(Iend+1,Jstr-1)=0.5*(w_set(Iend,Jstr-1)+w_set(Iend+1,Jstr))
        Sdens(Iend+1,Jstr-1)=0.5*(Sdens(Iend,Jstr-1)+Sdens(Iend+1,Jstr))
      endif
      if (NORTHERN_EDGE.and.WESTERN_EDGE) then
        taucb(Istr-1,Jend+1)=0.5*(taucb(Istr,Jend+1)+taucb(Istr-1,Jend))
        Ssize(Istr-1,Jend+1)=0.5*(Ssize(Istr,Jend+1)+Ssize(Istr-1,Jend))
        w_set(Istr-1,Jend+1)=0.5*(w_set(Istr,Jend+1)+w_set(Istr-1,Jend))
        Sdens(Istr-1,Jend+1)=0.5*(Sdens(Istr,Jend+1)+Sdens(Istr-1,Jend))
      endif
      if (NORTHERN_EDGE.and.EASTERN_EDGE) then
        taucb(Iend+1,Jend+1)=0.5*(taucb(Iend,Jend+1)+taucb(Iend+1,Jend))
        Ssize(Iend+1,Jend+1)=0.5*(Ssize(Iend,Jend+1)+Ssize(Iend+1,Jend))
        w_set(Iend+1,Jend+1)=0.5*(w_set(Iend,Jend+1)+w_set(Iend+1,Jend))
        Sdens(Iend+1,Jend+1)=0.5*(Sdens(Iend,Jend+1)+Sdens(Iend+1,Jend))
      endif
# endif

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 taucb(START_2D_ARRAY))
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 Ssize(START_2D_ARRAY))
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 w_set(START_2D_ARRAY))
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 Sdens(START_2D_ARRAY))
# endif
!
!-----------------------------------------------------------------------
! Set boundary conditions
!-----------------------------------------------------------------------
!
#  ifndef EW_PERIODIC
      if (EASTERN_EDGE) then
        do j=Jstr,Jend
           do k=1,NLAY
             do ised=1,NST
                bed_mass(Iend+1,j,k,nnew,ised)=
     &              bed_mass(Iend,j,k,nnew,ised)
                bed_frac(Iend+1,j,k,ised)=
     &              bed_frac(Iend,j,k,ised)
             enddo
           enddo
           bed_thick(Iend+1,j,1)=bed_thick(Iend,j,1)
         enddo
      endif
      if (WESTERN_EDGE) then
        do j=Jstr,Jend
          do k=1,NLAY
            do ised=1,NST
                bed_mass(Istr-1,j,k,nnew,ised)=
     &              bed_mass(Istr,j,k,nnew,ised)
                bed_frac(Istr-1,j,k,ised)=
     &              bed_frac(Istr,j,k,ised)
             enddo
           enddo
          bed_thick(Istr-1,j,1)=bed_thick(Istr,j,1)
        enddo
      endif
#  endif
#  ifndef NS_PERIODIC
      if (NORTHERN_EDGE) then
        do i=Istr,Iend
          do k=1,NLAY
            do ised=1,NST
               bed_mass(i,Jend+1,k,nnew,ised)=
     &              bed_mass(i,Jend,k,nnew,ised)
               bed_frac(i,Jend+1,k,ised)=
     &              bed_frac(i,Jend,k,ised)
            enddo
          enddo
          bed_thick(i,Jend+1,1)=bed_thick(i,Jend,1)
        enddo
      endif
      if (SOUTHERN_EDGE) then
        do i=Istr,Iend
          do k=1,NLAY
            do ised=1,NST
               bed_mass(i,Jstr-1,k,nnew,ised)=
     &              bed_mass(i,Jstr,k,nnew,ised)
               bed_frac(i,Jstr-1,k,ised)=
     &              bed_frac(i,Jstr,k,ised)
            enddo
          enddo
          bed_thick(i,Jstr-1,1)=bed_thick(i,Jstr,1)
        enddo
      endif
#  endif
#  if !defined EW_PERIODIC && !defined NS_PERIODIC
      if (SOUTHERN_EDGE.and.WESTERN_EDGE) then
        bed_mass(Istr-1,Jstr-1,:,nnew,:)=
     &  0.5*(bed_mass(Istr,Jstr-1,:,nnew,:)
     &      +bed_mass(Istr-1,Jstr,:,nnew,:))
        bed_frac(Istr-1,Jstr-1,:,:)=
     &  0.5*(bed_frac(Istr,Jstr-1,:,:)
     &      +bed_frac(Istr-1,Jstr,:,:))
        bed_thick(Istr-1,Jstr-1,1)=
     &  0.5*(bed_thick(Istr,Jstr-1,1)
     &      +bed_thick(Istr-1,Jstr,1))
      endif
      if (SOUTHERN_EDGE.and.EASTERN_EDGE) then
        bed_mass(Iend+1,Jstr-1,:,nnew,:)=
     &  0.5*(bed_mass(Iend,Jstr-1,:,nnew,:)
     &      +bed_mass(Iend+1,Jstr,:,nnew,:))
        bed_frac(Iend+1,Jstr-1,:,:)=
     &  0.5*(bed_frac(Iend,Jstr-1,:,:)
     &      +bed_frac(Iend+1,Jstr,:,:))
        bed_thick(Iend+1,Jstr-1,1)=
     &  0.5*(bed_thick(Iend,Jstr-1,1)
     &      +bed_thick(Iend+1,Jstr,1))
      endif
      if (NORTHERN_EDGE.and.WESTERN_EDGE) then
        bed_mass(Istr-1,Jend+1,:,nnew,:)=
     &  0.5*(bed_mass(Istr,Jend+1,:,nnew,:)
     &      +bed_mass(Istr-1,Jend,:,nnew,:))
        bed_frac(Istr-1,Jend+1,:,:)=
     &  0.5*(bed_frac(Istr,Jend+1,:,:)
     &      +bed_frac(Istr-1,Jend,:,:))
        bed_thick(Istr-1,Jend+1,1)=
     &  0.5*(bed_thick(Istr,Jend+1,1)
     &      +bed_thick(Istr-1,Jend,1))
      endif
      if (NORTHERN_EDGE.and.EASTERN_EDGE) then
        bed_mass(Iend+1,Jend+1,:,nnew,:)=
     &   0.5*(bed_mass(Iend,Jend+1,:,nnew,:)
     &       +bed_mass(Iend+1,Jend,:,nnew,:))
        bed_frac(Iend+1,Jend+1,:,:)=
     &   0.5*(bed_Frac(Iend,Jend+1,:,:)
     &       +bed_Frac(Iend+1,Jend,:,:))
        bed_thick(Iend+1,Jend+1,1)=
     &   0.5*(bed_thick(Iend,Jend+1,1)
     &       +bed_thick(Iend+1,Jend,1))
      endif
#  endif

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      do ised=1,NST
        do k=1,NLAY
          call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 bed_mass(START_2D_ARRAY,k,nnew,ised))
          call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 bed_frac(START_2D_ARRAY,k,ised))         
        enddo 
        call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                 bedldu(START_2D_ARRAY,ised))
        call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                 bedldv(START_2D_ARRAY,ised))
      enddo
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 bed_thick(START_2D_ARRAY,1))
#  endif

      return 
      end
# endif /* BEDLOAD */
      
# ifdef SUSPLOAD
!======================================================================
!                                                                      !
!  This routine computes the vertical settling (sinking) of suspended  !
!  sediment via a semi-Lagrangian advective flux algorithm. It uses a  !
!  parabolic,  vertical reconstructuion of the suspended  sediment in  !
!  the water column with PPT/WENO constraints to avoid oscillations.   !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Colella, P. and P. Woodward, 1984: The piecewise parabolic method   !
!    (PPM) for gas-dynamical simulations, J. Comp. Phys., 54, 174-201. !
!                                                                      !
!  Liu, X.D., S. Osher, and T. Chan, 1994: Weighted essentially        !
!    nonoscillatory shemes, J. Comp. Phys., 115, 200-212.              !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
!
      subroutine sed_settling_tile (Istr,Iend,Jstr,Jend,
     &                              Hz_inv,Hz_inv2,Hz_inv3,qc,
     &                              FC,qR,qL,WR,WL,tau_w,ksource) 
!
      implicit none
#  include "param.h"
#  include "ocean3d.h"
#  include "sediment.h"
#  include "bbl.h"
#  include "scalars.h"
#  include "grid.h"
#  include "mixing.h"
#  ifndef BBL
#   include "forces.h"
#  endif
#  ifdef DIAGNOSTICS_TS
#   include "diagnostics.h"
#  endif
!
      integer Istr,Iend,Jstr,Jend
      real Hz_inv(PRIVATE_1D_SCRATCH_ARRAY,N),  
     &    Hz_inv2(PRIVATE_1D_SCRATCH_ARRAY,N), 
     &    Hz_inv3(PRIVATE_1D_SCRATCH_ARRAY,N),  
     &         qc(PRIVATE_1D_SCRATCH_ARRAY,N),  
     &         FC(PRIVATE_1D_SCRATCH_ARRAY,0:N),  
     &         qR(PRIVATE_1D_SCRATCH_ARRAY,N), 
     &         qL(PRIVATE_1D_SCRATCH_ARRAY,N),  
     &         WR(PRIVATE_1D_SCRATCH_ARRAY,N), 
     &         WL(PRIVATE_1D_SCRATCH_ARRAY,N)
      real tau_w(PRIVATE_2D_SCRATCH_ARRAY)
      integer  ksource(PRIVATE_1D_SCRATCH_ARRAY,N),
     &         ised,indx,i,j,k,ks  
      real dltR, dltL,cff,cffR,cffL,cu
!
!-----------------------------------------------------------------------
!  Compute bottom skin-friction stress due to combined maximum wave
!  and current interaction
!-----------------------------------------------------------------------
!
      do j=Jstr,Jend
        do i=Istr,Iend
          tau_w(i,j)=sqrt(bustrw(i,j)**2
     &                   +bvstrw(i,j)**2)
#  ifdef WET_DRY
          tau_w(i,j)=tau_w(i,j)*rmask_wet(i,j)
#  endif
        enddo
      enddo
      do j=Jstr,Jend                        ! Auxiliary step: save   
        do k=1,N                            ! inverses of grid box 
          do i=Istr,Iend                    ! heights to avoid
            Hz_inv(i,k)=1./Hz(i,j,k)        ! repeated divisions
          enddo
        enddo
        do k=1,N-1
          do i=Istr,Iend
            Hz_inv2(i,k)=1./(Hz(i,j,k)+Hz(i,j,k+1))
          enddo
        enddo
        do k=2,N-1
          do i=Istr,Iend
            Hz_inv3(i,k)=1./(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          enddo
        enddo

!     
! Vertical sinking of suspended particles:   Copy concentration of
!========= ======= == ========= ==========   suspended sediment into
! scratch array "qc" (q-central, restrict it to be positive) which
! are hereafter interpreted as a set of grid-box averaged values for
! concentration. then reconstruct vertical profile of "qc" in terms
! of a set of parabolic segments within each grid box; and, finally,
! compute semi-Lagrangian flux due to sinking.
!
        do ised=1,NST
          indx=itrc_sed+ised-1
          do k=1,N
            do i=Istr,Iend
              qc(i,k)=t(i,j,k,nnew,indx)
            enddo
          enddo

          do k=N-1,1,-1
            do i=Istr,Iend
              FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
            enddo
          enddo
          do k=2,N-1
            do i=Istr,Iend
              dltR=Hz(i,j,k)*FC(i,k)
              dltL=Hz(i,j,k)*FC(i,k-1)
              cff=Hz(i,j,k-1)+2.*Hz(i,j,k)+Hz(i,j,k+1)
              cffR=cff*FC(i,k)
              cffL=cff*FC(i,k-1)            ! Apply PPM monotonicity
                                            ! constraint to prevent
              if (dltR*dltL .le. 0.) then   ! oscillation within the
                dltR=0.                     ! grid box
                dltL=0.
              elseif (abs(dltR) .gt. abs(cffL)) then
                dltR=cffL
              elseif (abs(dltL) .gt. abs(cffR)) then
                dltL=cffR
              endif                         ! Compute right and left
                                            ! side values qR,qL of 
              cff=(dltR-dltL)*Hz_inv3(i,k)  ! parabolic segments
              dltR=dltR-cff*Hz(i,j,k+1)     ! within grid box Hz(k)
              dltL=dltL+cff*Hz(i,j,k-1)     ! (WR,WL are measures of 
              qR(i,k)=qc(i,k)+dltR          ! quadratic variations).
              qL(i,k)=qc(i,k)-dltL
              WR(i,k)=( 2.*dltR-dltL )**2   ! NOTE: Although each 
              WL(i,k)=( dltR-2.*dltL )**2   ! parabolic segment is
            enddo                           ! monotone within its
          enddo          !--> discard FC    ! grid box, monotonicity 
                                            ! of the whole profile is
          cff=1.0E-14                       ! not guaranteed, because 
          do k=2,N-2                        ! qL(k+1)-qR(k) may still 
            do i=Istr,Iend                  ! have different sign 
              dltL=max(WL(i,k),   cff)      ! than qc(k+1)-qc(k)...
              dltR=max(WR(i,k+1), cff)
              qR(i,k)=(dltR*qR(i,k)+dltL*qL(i,k+1))/(dltR+dltL)
              qL(i,k+1)=qR(i,k)
            enddo                           !   ...this possibility
          enddo      !--> discard WR,WL     ! is excluded, after qL
                                            ! and qR are reconciled
          do i=Istr,Iend                    ! using WENO procedure.
            FC(i,N)=0.     !<-- no-flux BC
#  if defined LINEAR_CONTINUATION
            qL(i,N)=qR(i,N-1)
            qR(i,N)=2.*qc(i,N)-qL(i,N)
#  elif defined NEUMANN
            qL(i,N)=qR(i,N-1)
            qR(i,N)=1.5*qc(i,N)-0.5*qL(i,N)
#  else
            qR(i,N)=qc(i,N)                 ! Strictly monotone
            qL(i,N)=qc(i,N)                 ! version as the default:
            qR(i,N-1)=qc(i,N)               ! distributions at top...
#  endif
#  if defined LINEAR_CONTINUATION 
            qR(i,1)=qL(i,2)
            qL(i,1)=2.*qc(i,1)-qR(i,1)
#  elif defined NEUMANN
            qR(i,1)=qL(i,2)
            qL(i,1)=1.5*qc(i,1)-0.5*qR(i,1)
#  else  
            qL(i,2)=qc(i,1)                 !     ...and bottom grid
            qR(i,1)=qc(i,1)                 ! boxes re assumed to be
            qL(i,1)=qc(i,1)                 ! piecewise constant.
#  endif
          enddo

          do k=1,N                          ! Since the reconciled 
            do i=Istr,Iend                  ! interfacial values may
              dltR=qR(i,k)-qc(i,k)          ! cause non-monotonic
              dltL=qc(i,k)-qL(i,k)          ! behavior of parabolic
              cffR=2.*dltR                  ! segments inside grid
              cffL=2.*dltL                  ! box apply monotonicity
                                            ! constraint again.
              if (dltR*dltL .lt. 0.) then
                dltR=0.
                dltL=0.
              elseif (abs(dltR) .gt. abs(cffL)) then
                dltR=cffL
              elseif (abs(dltL) .gt. abs(cffR)) then
                dltL=cffR
              endif

              qR(i,k)=qc(i,k)+dltR
              qL(i,k)=qc(i,k)-dltL
            enddo
          enddo     !--> discard everything, except qR,qL
!
! After this moment reconstruction is considered complete. The next
! stage is to compute vertical advective fluxes FC. It is expected
! that sinking may occurs relatively fast, the algorithm is designed
! to be free of CFL criterion, which is achieved by allowing
! integration bounds for semi-Lagrangian advective flux to use as
! many grid boxes in upstream direction as necessary.
!
          cff=dt*abs(Wsed(ised))          ! In the two code segments
          do k=1,N                        ! WL is z-coordinate of the
            do i=Istr,Iend                ! departure point for grid
              FC(i,k-1)=0.                ! box interface z_w with
              WL(i,k)=z_w(i,j,k-1)+cff    ! the same indices;
#  ifdef WET_DRY
     &                *rmask_wet(i,j)
#  endif
              WR(i,k)=Hz(i,j,k)*qc(i,k)   ! FC is finite volume flux;
              ksource(i,k)=k              ! ksource(:,k) is index of
            enddo                         ! vertical grid box which
          enddo                           ! contains the departure
          do k=1,N                        ! point (restricted by N);
            do ks=k,N-1                   ! During the search: also
              do i=Istr,Iend
                if (WL(i,k) .gt. z_w(i,j,ks)) then
                  ksource(i,k)=ks+1
                  FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                endif
              enddo                       ! add in content of whole
            enddo                         ! grid boxes participating
          enddo       !--> discard WR     ! in FC.

          do k=1,N                        ! Finalize computation of
            do i=Istr,Iend                ! flux: add fractional part
              ks=ksource(i,k)
              cu=min(1.,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
              FC(i,k-1)=FC(i,k-1) + Hz(i,j,ks)*cu*( qL(i,ks)
     &                             +cu*( 0.5*(qR(i,ks)-qL(i,ks))
     &             -(1.5-cu)*(qR(i,ks)+qL(i,ks)-2.*qc(i,ks))  ))
            enddo
          enddo
          cff=1/tau_cd(ised) 
          do i=Istr,Iend
            !FC(i,0)=FC(i,0)*min(1.,REAL(int(cff*tau_w(i,j))))
            do k=1,N,+1
               t(i,j,k,nnew,indx)=qc(i,k)
     &                   +(FC(i,k)-FC(i,k-1))*Hz_inv(i,k) 
#  ifdef SPONGE_SED
     &                                       *sed_sponge(i,j)
#  endif
            enddo   
            settling_flux(i,j,ised)=FC(i,0)
#  ifdef SPONGE_SED
     &                                       *sed_sponge(i,j)
#  endif
          enddo   
!
!  Diagnostics for settling: use TForc array
!
#  ifdef DIAGNOSTICS_TS
          do k=1,N 
            do i=Istr,Iend
              TForc(i,j,k,indx)=(FC(i,k)-FC(i,k-1))/dt
#   ifdef MASKING
     &                                    * rmask(i,j)
#   endif
            enddo
          enddo
#  endif /* DIAGNOSTICS_TS */

        enddo   ! ised loop
      enddo     ! j loop

!
!  Update Trate diagnostics
!
#  ifdef DIAGNOSTICS_TS
      do ised=1,NST
        indx=itrc_sed+ised-1  
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
              Trate(i,j,k,indx)=(Hz(i,j,k)*t(i,j,k,nnew,indx)
     &                          -Hz_bak(i,j,k)*t(i,j,k,nstp,indx))
     &                                                         /dt
#   ifdef MASKING
     &                                                 *rmask(i,j)
#   endif
            enddo
          enddo
        enddo
      enddo
#  endif /* DIAGNOSTICS_TS */
!
!  Exchange for output purpose
!
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      do ised=1,NST
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                          settling_flux(START_2D_ARRAY,ised))
      enddo
#  endif

      return
      end
!
!=======================================================================
!                                                                      !
!  This computes sediment erosion flux and resuspended load.           !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
!
      subroutine sed_fluxes_tile (Istr,Iend,Jstr,Jend, tau_w) 
!
      implicit none
#  include "param.h"
#  include "grid.h"
#  include "bbl.h"
#  include "forces.h"
#  include "sediment.h"
#  include "ocean3d.h"
#  include "scalars.h"
#  include "mixing.h"
#  ifdef DIAGNOSTICS_TS
#   include "diagnostics.h"
#  endif
      integer Istr,Iend,Jstr,Jend
      integer i,j,ised,ised0
      integer indx
      integer imin,jmin,imax,jmax
      real cff,cff1,cff2,cff3,cff4,cff5,tau_water,Zr
      real tau_w(PRIVATE_2D_SCRATCH_ARRAY)
#  ifdef BED_ARMOR
      integer ised0
      real ph,pe
#  endif
!
!-----------------------------------------------------------------------
!  Compute bottom skin-friction stress due to combined maximum wave 
!  and current interaction
!-----------------------------------------------------------------------
!
      do j=Jstr,Jend
        do i=Istr,Iend
          tau_w(i,j)=sqrt(bustrw(i,j)**2
     &                   +bvstrw(i,j)**2)
#  ifdef WET_DRY
          tau_w(i,j)=tau_w(i,j)*rmask_wet(i,j)
#  endif
        enddo
      enddo
!
!-----------------------------------------------------------------------
!  Compute erosion flux and resuspended load
!-----------------------------------------------------------------------
!
      do j=jstr,jend
        do ised=1,NST
          do i=istr,iend
            indx=itrc_sed+ised-1     
! 
!  Calculate critical shear stress in Pa
!  => change for COHESIVE bed here
!

#  if defined COHESIVE_BED
            cff = rho0/bed_bctr(i,j,1)
#  elif defined MIXED_BED
            cff = MAX(bot_dprp(i,j)*bed_bctr(i,j,1)/rho0+          
     &            (1.-bot_dprp(i,j))*tau_ce(ised),           
     &            tau_ce(ised))
            cff=1./cff
#  else
            cff=1./tau_ce_2d(i,j,ised)
#  endif
!
!  Compute erosion, ero_flux (kg/m2).
!
            cff1=(1.0-bed_poros(i,j,1))*bed_frac(i,j,1,ised)
            cff2=dt*Erate(ised)*cff1
            cff3=Srho(ised)*cff1
            cff4=bed_mass(i,j,1,bnew,ised)
            cff5=settling_flux(i,j,ised)
#  if defined SED_TAU_CD_CONST || defined SED_TAU_CD_LIN
             if (tau_w(i,j).GT.tau_cd(ised) ) then
                cff5=0.0
             endif
#  endif
# if defined SED_TAU_CD_LIN
             if (tau_w(i,j).LE.tau_cd(ised).AND.                      
     &                       tau_cd(ised).GT.0.0) then
                cff5=(1.0-tau_w(i,j)/tau_cd(ised))*cff5
             endif
# endif            
            ero_flux(i,j,ised)=                                     
     &            MIN(MAX(0.,cff2*(cff*tau_w(i,j)-1.)),      
     &                MIN(cff3*bot_thick(i,j) MORPH_FAC ,cff4)+
     &                cff5)
            ero_flux(i,j,ised)=ero_flux(i,j,ised)
#  ifdef SPONGE_SED
     &                             *sed_sponge(i,j)
#  endif
!
!  Update global tracer variables (m Tunits for nnew indx, Tuints for 3)
!  for erosive flux.
!                        
            t(i,j,1,nnew,indx)=t(i,j,1,nnew,indx)+ero_flux(i,j,ised)
     &                                                    /Hz(i,j,1)
#  if defined SED_TAU_CD_CONST || defined SED_TAU_CD_LIN
     &                       + (settling_flux(i,j,ised)-cff5)/Hz(i,j,1)
            settling_flux(i,j,ised)=cff5
#  endif
          enddo  ! i loop

!
!  Update erosion/settling diagnostics (bottom level only)
!
#  ifdef DIAGNOSTICS_TS
          do i=Istr,Iend
            TForc(i,j,1,indx) = TForc(i,j,1,indx) +
     &                          ero_flux(i,j,ised)/dt
#   ifdef MASKING
     &                                   * rmask(i,j)
#   endif
          enddo
#  endif /* DIAGNOSTICS_TS */

        enddo    ! ised loop 
      enddo      ! j loop

!
!  Update Trate diagnostics (bottom level only)
!
#  ifdef DIAGNOSTICS_TS
      do ised=1,NST
        indx=itrc_sed+ised-1  
        do j=Jstr,Jend
          do i=Istr,Iend
            Trate(i,j,1,indx)=(Hz(i,j,1)*t(i,j,1,nnew,indx)
     &                        -Hz_bak(i,j,1)*t(i,j,1,nstp,indx))
     &                                                       /dt
#   ifdef MASKING
     &                                               *rmask(i,j)
#   endif
          enddo
        enddo
      enddo
#  endif /* DIAGNOSTICS_TS */
!
!  Exchange for output purpose
!
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      do ised=1,NST
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                          ero_flux(START_2D_ARRAY,ised))
      enddo
#   if defined SED_TAU_CD_CONST || defined SED_TAU_CD_LIN
      do ised=1,NST
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                          settling_flux(START_2D_ARRAY,ised))
      enddo  
#   endif 
#  endif

      return
      end
# endif /* SUSPLOAD */
!
!=======================================================================
!                                                                      !
!  This routine computes sediment bed layer stratigraphy.              !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
!
# if !defined COHESIVE_BED && !defined MIXED_BED

      subroutine sed_bed_tile (Istr,Iend,Jstr,Jend,
     &                         tau_w,dep_mass,Awrk1,Awrk2)
!
      implicit none
# include "param.h"
# include "bbl.h"
# include "forces.h"
# include "grid.h"
# include "ocean3d.h"
# include "scalars.h"
# include "sediment.h"
      integer Istr,Iend,Jstr,Jend
      integer i,j,k,ks,ised,Ksed, NLAYm1
      real tau_w(PRIVATE_2D_SCRATCH_ARRAY)
      real dep_mass(PRIVATE_1D_SCRATCH_ARRAY,NST)
      real Awrk1(PRIVATE_2D_SCRATCH_ARRAY)
      real Awrk2(PRIVATE_2D_SCRATCH_ARRAY)
      real eps,cff,cff1,cff2,cff3
      real thck_avail,thck_to_add,newlayer_thick
      integer imin,imax,jmin,jmax
      parameter ( eps=1.D-20, newlayer_thick=0.001 )
!
!-----------------------------------------------------------------------
!  Compute bottom skin-friction stress due to combined maximum wave 
!  and current interaction
!-----------------------------------------------------------------------
!
      do j=Jstr,Jend
        do i=Istr,Iend
          tau_w(i,j)=sqrt(bustrw(i,j)**2
     &                   +bvstrw(i,j)**2)
!
!-----------------------------------------------------------------------
!  Calculate active layer thickness
!     --> Harris and Wiberg 1997: za = k1 (tau-tc)rho + k2 Sd50
!-----------------------------------------------------------------------
!
          bot_thick(i,j)=MAX(newlayer_thick,            
     &                       0.007*(tau_w(i,j)-taucb(i,j))*rho0)+   
     &                       6.*Ssize(i,j)
        enddo
      enddo
!-----------------------------------------------------------------------
!  Update bed properties according to ero_flux and dep_flux.
!-----------------------------------------------------------------------
!
!
      NLAYm1=max(NLAY-1,1)  ! for NLAY=1 case
!
      do j=Jstr,Jend
!
        do ised=1,NST
          do i=Istr,Iend
            dep_mass(i,ised)=0.0

# ifdef SUSPLOAD
!  Suspended erosion/deposition. Apply morphology factor.
!
!  The deposition and resuspension of sediment on the bottom "bed"
!  is due to precipitation flux, already computed, and the
!  resuspension (erosion, hence called ero_flux). The resuspension is
!  applied to the bottom-most grid box value qc(:,1) so the total mass
!  is conserved. Restrict "ero_flux" so that "bed" cannot go negative
!  after both fluxes are applied.
!
            ero_flux(i,j,ised)=ero_flux(i,j,ised) MORPH_FAC
            settling_flux(i,j,ised)=settling_flux(i,j,ised) MORPH_FAC
# endif
!
! Net deposited sediment, including bedload
!
            cff=
# ifdef SUSPLOAD
     &          settling_flux(i,j,ised)-ero_flux(i,j,ised)
# endif
# ifdef BEDLOAD
     &          +bed_mass(i,j,1,nnew,ised)-bed_mass(i,j,1,nstp,ised)
# endif

            if ( cff .gt. 0. ) then
!
!  If first time step of deposit, then store deposit material in
!  temporary array, dep_mass.
!
!  Use initial resolution Bthk(1) to trigger new layer formation
!  (prevents stacking all layers near the surface - PM).
! 
              if ( time .gt. (bed_age(i,j,1)+1.1*dt) .and.
!     &             bed_thick(i,j,1) .gt. newlayer_thick ) then
     &             bed_thick(i,j,1) .gt. Bthk(1) ) then
                dep_mass(i,ised)=cff
              endif
              !bed_age(i,j,1)=time   !!!???
            endif
!
!  Update bed mass arrays.
!
            bed_mass(i,j,1,nnew,ised)=MAX(bed_mass(i,j,1,nstp,ised)+ 
     &                                    cff, 0.) 
            do k=2,NLAY
              bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k,nstp,ised)
            enddo
          enddo   ! i loop
        enddo     ! ised loop
!
!  If first time step of deposit, create new layer and combine bottom
!  two bed layers.
!
        do i=Istr,Iend
          cff=0.
!
!  Determine if deposition ocurred here.
!
          if (NLAY.gt.1) then
!
            do ised=1,NST
              cff=cff+dep_mass(i,ised)
            enddo
!
            if (cff.gt.0.) then
!     
!  Combine bottom layers.
!
              bed_poros(i,j,NLAY)=0.5*(bed_poros(i,j,NLAYm1)+        
     &                                 bed_poros(i,j,NLAY))
              bed_age  (i,j,NLAY)=0.5*(bed_age  (i,j,NLAYm1)+     
     &                                 bed_age  (i,j,NLAY))
              do ised=1,NST
                bed_mass(i,j,NLAY,nnew,ised)=                           
     &                             bed_mass(i,j,NLAYm1,nnew,ised)+      
     &                             bed_mass(i,j,NLAY  ,nnew,ised)
              enddo
!
!  Push layers down.
!
              do k=NLAY-1,2,-1
                bed_poros(i,j,k)=bed_poros(i,j,k-1)
                bed_age  (i,j,k)=bed_age  (i,j,k-1)
                do ised =1,NST
                  bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k-1,nnew,ised)
                enddo
              enddo
!
!  Set new top layer parameters.
!
              k=min(2,NLAY) ! for compilation if NLAY=1
              do ised=1,NST
                bed_mass(i,j,k,nnew,ised)=MAX(bed_mass(i,j,k,nnew,ised)-
     &                                        dep_mass(i,ised),0.)
                bed_mass(i,j,1,nnew,ised)=dep_mass(i,ised)
              enddo
              bed_age(i,j,1)=time  !!! update here
!
            endif ! new deposition cff>0
          endif  ! NLAY=1
!
!  Recalculate thickness and fractions for all layers.
!
          do k=1,NLAY
            cff3=0.
            do ised=1,NST
              cff3=cff3+bed_mass(i,j,k,nnew,ised)
            enddo
            if (cff3.eq.0.) then
              cff3=eps
            endif
            bed_thick(i,j,k)=0.
            do ised=1,NST
              bed_frac(i,j,k,ised)=bed_mass(i,j,k,nnew,ised)/cff3
              bed_thick(i,j,k)=MAX(bed_thick(i,j,k)+                    
     &                         bed_mass(i,j,k,nnew,ised)/               
     &                         (Srho(ised)*(1.-bed_poros(i,j,k))),0.)
            enddo  ! ised loop 
          enddo    ! k loop   
        enddo      ! i loop
      enddo        ! j loop
!
!  Ensure top bed layer thickness is greater or equal than active layer
!  thickness. If need to add sed to top layer, then entrain from lower
!  levels. Create new layers at bottom to maintain Nlay
!
      do j=Jstr,Jend
        do i=Istr,Iend
          if (bot_thick(i,j).gt.bed_thick(i,j,1)) then
            if (NLAY.eq.1) then
              bot_thick(i,j)=bed_thick(i,j,1)
            else
              thck_to_add=bot_thick(i,j)-bed_thick(i,j,1)
              thck_avail=0.
              Ksed=1       ! initialize
              do k=2,NLAY
                if (thck_avail.lt.thck_to_add) then
                  thck_avail=thck_avail+bed_thick(i,j,k)
                  Ksed=k
                endif
              enddo
!
!  Catch here if there was not enough bed material.
!
              if (thck_avail.lt.thck_to_add) then
                bot_thick(i,j)=bed_thick(i,j,1)+thck_avail
                thck_to_add=thck_avail
              endif
!
!  Update bed mass of top layer and fractional layer.
!
              cff2=MAX(thck_avail-thck_to_add,0.)/            
     &             MAX(bed_thick(i,j,Ksed),eps)
              do ised=1,NST
                cff1=0.
                do k=1,Ksed
                  cff1=cff1+bed_mass(i,j,k,nnew,ised)
                enddo
                cff3=cff2*bed_mass(i,j,Ksed,nnew,ised)
                bed_mass(i,j,1   ,nnew,ised)=cff1-cff3
                bed_mass(i,j,Ksed,nnew,ised)=cff3
              enddo
!
!  Update thickness of fractional layer ksource_sed
!
              bed_thick(i,j,Ksed)=MAX(thck_avail-thck_to_add,0.)
!
!  Upate bed fraction of top layer.
!
              cff3=0.
              do ised=1,NST
                cff3=cff3+bed_mass(i,j,1,nnew,ised)
              enddo
              if (cff3.eq.0.) then
                cff3=eps
              endif
              do ised=1,NST
                bed_frac(i,j,1,ised)=bed_mass(i,j,1,nnew,ised)/cff3
              enddo
!
!  Upate bed thickness of top layer.
!
              bed_thick(i,j,1)=bot_thick(i,j)
!
!  Pull all layers closer to the surface.
!
              ks=Ksed-2  ! ks is the number of new layers
              do k=Ksed,NLAY
                bed_thick(i,j,k-ks)=bed_thick(i,j,k)
                bed_poros(i,j,k-ks)=bed_poros(i,j,k)
                bed_age  (i,j,k-ks)=bed_age  (i,j,k)
                do ised=1,NST
                  bed_frac(i,j,k-ks,ised)=bed_frac(i,j,k,ised)
                  bed_mass(i,j,k-ks,nnew,ised)=bed_mass(i,j,k,nnew,ised)
                enddo
              enddo
!
!  Add new layers onto the bottom. Split what was in the bottom layer to
!  fill these new empty cells. 
!
              !ks=0
              cff=1./REAL(ks+1)
              do k=NLAY,NLAY-ks,-1
                bed_thick(i,j,k)=bed_thick(i,j,NLAY-ks)*cff
                bed_age  (i,j,k)=bed_age  (i,j,NLAY-ks)
                do ised=1,NST
                  bed_frac(i,j,k,ised)=bed_frac(i,j,NLAY-ks,ised)
                  bed_mass(i,j,k,nnew,ised)=                          
     &                             bed_mass(i,j,NLAY-ks,nnew,ised)*cff
                enddo   ! ised loop
              enddo     ! k loop
            endif       ! NLAY > 1
          endif         ! increase top bed layer
        enddo           ! i loop
      enddo             ! j loop
!
! Set boundary conditions
!
# ifndef EW_PERIODIC
      if (EASTERN_EDGE) then
        do j=Jstr,Jend
          do k=1,NLAY
            do ised=1,NST
            bed_mass(Iend+1,j,k,nnew,ised)=
     &                       bed_mass(Iend,j,k,nnew,ised)
            bed_frac(Iend+1,j,k,ised)=
     &                        bed_frac(Iend,j,k,ised)
             enddo
             bed_thick(Iend+1,j,k)=bed_thick(Iend,j,k)
             bed_poros(Iend+1,j,k)=bed_poros(Iend,j,k)
             bed_age(Iend+1,j,k)=bed_age(Iend,j,k)
           enddo
           bot_thick(Iend+1,j)=bot_thick(Iend,j)
         enddo
        endif
      if (WESTERN_EDGE) then
        do j=Jstr,Jend
          do k=1,NLAY
            do ised=1,NST
            bed_mass(Istr-1,j,k,nnew,ised)=
     &                       bed_mass(Istr,j,k,nnew,ised)
            bed_frac(Istr-1,j,k,ised)=
     &                        bed_frac(Istr,j,k,ised)
             enddo
             bed_thick(Istr-1,j,k)=bed_thick(Istr,j,k)
             bed_poros(Istr-1,j,k)=bed_poros(Istr,j,k)
             bed_age(Istr-1,j,k)=bed_age(Istr,j,k)
           enddo
           bot_thick(Istr-1,j)=bot_thick(Istr,j)
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (NORTHERN_EDGE) then
        do i=Istr,Iend
          do k=1,NLAY
            do ised=1,NST
              bed_mass (i,Jend+1,k,nnew,ised)=
     &                            bed_mass(i,Jend,k,nnew,ised)
              bed_frac (i,Jend+1,k,ised)=
     &                            bed_frac(i,Jend,k,ised)
            enddo
            bed_thick(i,Jend+1,k)=bed_thick(i,Jend,k)
            bed_poros(i,Jend+1,k)=bed_poros(i,Jend,k)
            bed_age  (i,Jend+1,k)=bed_age(i,Jend,k)
          enddo
          bot_thick  (i,Jend+1)=bot_thick(i,Jend)
        enddo
      endif
      if (SOUTHERN_EDGE) then
        do i=Istr,Iend
          do k=1,NLAY
            do ised=1,NST
              bed_mass (i,Jstr-1,k,nnew,ised)=
     &                            bed_mass(i,Jstr,k,nnew,ised)
              bed_frac (i,Jstr-1,k,ised)=
     &                            bed_frac(i,Jstr,k,ised)
            enddo
            bed_thick(i,Jstr-1,k)=bed_thick(i,Jstr,k)
            bed_poros(i,Jstr-1,k)=bed_poros(i,Jstr,k)
            bed_age  (i,Jstr-1,k)=bed_age(i,Jend,k)
          enddo
          bot_thick  (i,Jstr-1)=bot_thick(i,Jstr)
        enddo
      endif
# endif
# if !defined EW_PERIODIC && !defined NS_PERIODIC
      if (SOUTHERN_EDGE.and.WESTERN_EDGE) then
        bed_mass(Istr-1,Jstr-1,:,nnew,:)=
     &  0.5*(bed_mass(Istr,Jstr-1,:,nnew,:)
     &      +bed_mass(Istr-1,Jstr,:,nnew,:))
        bed_frac(Istr-1,Jstr-1,:,:)=
     &  0.5*(bed_frac(Istr,Jstr-1,:,:)
     &      +bed_frac(Istr-1,Jstr,:,:))
        bed_thick(Istr-1,Jstr-1,:)=
     &  0.5*(bed_thick(Istr,Jstr-1,:)
     &      +bed_thick(Istr-1,Jstr,:))
        bed_poros(Istr-1,Jstr-1,:)=
     &  0.5*(bed_poros(Istr,Jstr-1,:)
     &      +bed_poros(Istr-1,Jstr,:))
        bed_age(Istr-1,Jstr-1,:)=
     &  0.5*(bed_age(Istr,Jstr-1,:)
     &      +bed_age(Istr-1,Jstr,:))
      endif
      if (SOUTHERN_EDGE.and.EASTERN_EDGE) then
        bed_mass(Iend+1,Jstr-1,:,nnew,:)=
     &  0.5*(bed_mass(Iend,Jstr-1,:,nnew,:)
     &      +bed_mass(Iend+1,Jstr,:,nnew,:))
        bed_frac(Iend+1,Jstr-1,:,:)=
     &  0.5*(bed_frac(Iend,Jstr-1,:,:)
     &      +bed_frac(Iend+1,Jstr,:,:))
        bed_thick(Iend+1,Jstr-1,:)=
     &  0.5*(bed_thick(Iend,Jstr-1,:)
     &      +bed_thick(Iend+1,Jstr,:))
        bed_poros(Iend+1,Jstr-1,:)=
     &  0.5*(bed_poros(Iend,Jstr-1,:)
     &      +bed_poros(Iend+1,Jstr,:))
        bed_age(Iend+1,Jstr-1,:)=
     &  0.5*(bed_age(Iend,Jstr-1,:)
     &      +bed_age(Iend+1,Jstr,:))
      endif
      if (NORTHERN_EDGE.and.WESTERN_EDGE) then
        bed_mass(Istr-1,Jend+1,:,nnew,:)=
     &  0.5*(bed_mass(Istr,Jend+1,:,nnew,:)
     &      +bed_mass(Istr-1,Jend,:,nnew,:))
        bed_frac(Istr-1,Jend+1,:,:)=
     &  0.5*(bed_frac(Istr,Jend+1,:,:)
     &      +bed_frac(Istr-1,Jend,:,:))
        bed_thick(Istr-1,Jend+1,:)=
     &  0.5*(bed_thick(Istr,Jend+1,:)
     &      +bed_thick(Istr-1,Jend,:))
        bed_poros(Istr-1,Jend+1,:)=
     &  0.5*(bed_poros(Istr,Jend+1,:)
     &      +bed_poros(Istr-1,Jend,:))
        bed_age(Istr-1,Jend+1,:)=
     &  0.5*(bed_age(Istr,Jend+1,:)
     &      +bed_age(Istr-1,Jend,:))
      endif
      if (NORTHERN_EDGE.and.EASTERN_EDGE) then
        bed_mass(Iend+1,Jend+1,:,nnew,:)=
     &   0.5*(bed_mass(Iend,Jend+1,:,nnew,:)
     &       +bed_mass(Iend+1,Jend,:,nnew,:))
        bed_frac(Iend+1,Jend+1,:,:)=
     &   0.5*(bed_frac(Iend,Jend+1,:,:)
     &       +bed_frac(Iend+1,Jend,:,:))
        bed_thick(Iend+1,Jend+1,:)=
     &   0.5*(bed_thick(Iend,Jend+1,:)
     &       +bed_thick(Iend+1,Jend,:))
        bed_poros(Iend+1,Jend+1,:)=
     &   0.5*(bed_poros(Iend,Jend+1,:)
     &       +bed_poros(Iend+1,Jend,:))
        bed_age(Iend+1,Jend+1,:)=
     &   0.5*(bed_age(Iend,Jend+1,:)
     &       +bed_age(Iend+1,Jend,:))
      endif
# endif

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      do ised=1,NST
        do k=1,NLAY
          call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 bed_frac(START_2D_ARRAY,k,ised))
          call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 bed_mass(START_2D_ARRAY,k,nnew,ised))
        enddo 
      enddo
      do k=1,NLAY
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 bed_poros(START_2D_ARRAY,k))
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 bed_thick(START_2D_ARRAY,k))
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 bed_age(START_2D_ARRAY,k))
      enddo   
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 bot_thick(START_2D_ARRAY))
# endif
!
!-----------------------------------------------------------------------
!  Store total bed thickness
!-----------------------------------------------------------------------
!
# ifdef MORPHODYN
      do j=Jstr,Jend
        do i=Istr,Iend
          bed_thick_tot(i,j,nnew)=0.
          do k=1,NLAY
            bed_thick_tot(i,j,nnew)=bed_thick_tot(i,j,nnew)+
     &                              bed_thick(i,j,k)
          enddo
        enddo
      enddo
!
! Set boundary conditions
!
#  ifndef EW_PERIODIC
      if (EASTERN_EDGE) then
        do j=Jstr,Jend
          bed_thick_tot(Iend+1,j,nnew)=bed_thick_tot(Iend,j,nnew)
        enddo
      endif
      if (WESTERN_EDGE) then
        do j=Jstr,Jend
          bed_thick_tot(Istr-1,j,nnew)=bed_thick_tot(Istr,j,nnew)
        enddo
      endif
#  endif
#  ifndef NS_PERIODIC
      if (NORTHERN_EDGE) then
        do i=Istr,Iend
          bed_thick_tot(i,Jend+1,nnew)=bed_thick_tot(i,Jend,nnew)
        enddo
      endif
      if (SOUTHERN_EDGE) then
        do i=Istr,Iend
          bed_thick_tot(i,Jstr-1,nnew)=bed_thick_tot(i,Jstr,nnew)
        enddo
      endif
#  endif
#  if !defined EW_PERIODIC && !defined NS_PERIODIC
      if (SOUTHERN_EDGE.and.WESTERN_EDGE) then
        bed_thick_tot(Istr-1,Jstr-1,nnew)=
     &  0.5*(bed_thick_tot(Istr,Jstr-1,nnew)
     &      +bed_thick_tot(Istr-1,Jstr,nnew))
      endif
      if (SOUTHERN_EDGE.and.EASTERN_EDGE) then
        bed_thick_tot(Iend+1,Jstr-1,nnew)=
     &  0.5*(bed_thick_tot(Iend,Jstr-1,nnew)
     &      +bed_thick_tot(Iend+1,Jstr,nnew))
      endif
      if (NORTHERN_EDGE.and.WESTERN_EDGE) then
        bed_thick_tot(Istr-1,Jend+1,nnew)=
     &  0.5*(bed_thick_tot(Istr,Jend+1,nnew)
     &      +bed_thick_tot(Istr-1,Jend,nnew))
      endif
      if (NORTHERN_EDGE.and.EASTERN_EDGE) then
        bed_thick_tot(Iend+1,Jend+1,nnew)=
     &   0.5*(bed_thick_tot(Iend,Jend+1,nnew)
     &       +bed_thick_tot(Iend+1,Jend,nnew))
      endif
#  endif
!
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
       call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                         bed_thick_tot(START_2D_ARRAY,nnew))
#  endif
!
!-----------------------------------------------------------------------
!  Store bed evolution
!-----------------------------------------------------------------------
!
      do j=Jstr-1,Jend+1
        do i=Istr-1,Iend+1
          dh(i,j)=bed_thick_tot(i,j,nstp)-bed_thick_tot(i,j,nnew)
        enddo
      enddo
!
# endif /* MORPHODYN */

      return
      end

# else
**********************************************************************
      subroutine sed_bed_cohesive_tile (                       
     &                         Istr,Iend,Jstr,Jend,tau_w,dep_mass)
!***********************************************************************
!
      implicit none
# include "param.h"
# include "bbl.h"
# include "forces.h"
# include "grid.h"
# include "ocean3d.h"
# include "scalars.h"
# include "sediment.h"
      integer Istr,Iend,Jstr,Jend
      real tau_w(PRIVATE_2D_SCRATCH_ARRAY)
      real dep_mass(PRIVATE_1D_SCRATCH_ARRAY,NST)
      integer Ksed, i, ised, j, k, ks, nnn, indx
      real cff, cff1, cff2, cff3,cff4,eps
      real thck_avail, thck_to_add
      real nlysm(NST)   !!!!???(OpenMP)
      integer imin,imax,jmin,jmax
      real add,newlayer_thick

!     parameter (eps = 1.0E-14,newlayer_thick=0.0005) !BEDTOY
      parameter (eps = 1.0E-14,newlayer_thick=0.001)

# if defined MIXED_BED
      real pcoh 
# endif
# if defined COHESIVE_BED || defined MIXED_BED
      real alpha, bzactv, frt, tcb_top, tcb_bot
      real tcr(NLAY)  !!!!???(OpenMP)
# endif
# if defined COHESIVE_BED || defined MIXED_BED
      real bmz(NLAY)
# endif
# if defined SED_FLOCS && defined SED_DEFLOC
      real masseq(NMUD)
# endif

! KLUDGE alert
! minlayer_thick = 0.0005
      minlayer_thick = newlayer_thick
!
!-----------------------------------------------------------------------
! Compute sediment bed layer stratigraphy.
!-----------------------------------------------------------------------
!
! why MPM ???
# if defined BEDLOAD_MPM || defined SUSPLOAD
      do j=Jstr-1,Jend+1
        do i=Istr-1,Iend+1
          tau_w(i,j)=SQRT(bustrw(i,j)**2              
     &                   +bvstrw(i,j)**2)
#  ifdef WET_DRY
          tau_w(i,j)=tau_w(i,j)*rmask_wet(i,j)
#  endif
        enddo
      enddo
# endif
!
!-----------------------------------------------------------------------
!  Update bed properties according to ero_flux and dep_flux.
!-----------------------------------------------------------------------
!
# ifdef SUSPLOAD
      J_LOOP : do j=Jstr,Jend
!
!  The deposition and resuspension of sediment on the bottom "bed"
!  is due to precipitation flux FC(:,0), already computed, and the
!  resuspension (erosion, hence called ero_flux). The resuspension is
!  applied to the bottom-most grid box value qc(:,1) so the total mass
!  is conserved. Restrict "ero_flux" so that "bed" cannot go negative
!  after both fluxes are applied.
!
      do i=Istr,Iend
         SED_LOOP: do ised=1,NST  
            dep_mass(i,ised)=0.0
#  ifdef SED_MORPH
            ! Apply morphology factor.
             ero_flux(i,j,ised)=ero_flux(i,j,ised)   MORPH_FAC
             settling_flux(i,j,ised)=settling_flux(i,j,ised)  MORPH_FAC
#  endif

            !  Update bed mass arrays.
            bed_mass(i,j,1,nnew,ised)=MAX(bed_mass(i,j,1,bnew,ised)-
     &                                   (ero_flux(i,j,ised)-
     &                                    settling_flux(i,j,ised)), 0.)
            do k=2,NLAY
               bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k,nstp,ised)
            enddo
         enddo SED_LOOP
      
         cff3=0.0
         do ised=1,NST
            cff3=cff3+bed_mass(i,j,1,nnew,ised)
         enddo

         if (cff3.eq.0.0) then 
            cff3=eps 
         endif 
         bed_thick(i,j,1)=0.0
         do ised=1,NST
            bed_frac(i,j,1,ised)=bed_mass(i,j,1,nnew,ised)/cff3
            bed_thick(i,j,1)=MAX(bed_thick(i,j,1)+                      
     &        bed_mass(i,j,1,nnew,ised)/                                
     &        (Srho(ised)*                                         
     &        (1.0-bed_poros(i,j,1))),0.0)
         enddo
      enddo
      
      enddo J_LOOP 
# endif /* SUSPLOAD section */
!
!-----------------------------------------------------------------------
!  At this point, all deposition or erosion is complete, and
!  has been added/subtracted to top layer. Thickness has NOT been corrected.
!-----------------------------------------------------------------------
!
!!!!! taucb ????
      J_LOOP_CB : do j=Jstr,Jend
        I_LOOP_CB: do i=Istr,Iend
!       Calculate active layer thickness, bot_thick(i,j).
!       (trunk version allows this to be zero...this has minimum of 6*D50)
          bot_thick(i,j)=MAX(0.0,                                 
     &                          0.007*                               
     &                          (tau_w(i,j)-taucb(i,j))*rho0)+   
     &                          6.0*Ssize(i,j)
# if defined COHESIVE_BED 
          bot_thick(i,j)=MAX(0.0,                                 
     &                          0.007*                               
     &                          (tau_w(i,j)*rho0-bed_bctr(i,j,1)))+    
     &                          6.0*Ssize(i,j)
# endif
# if defined MIXED_BED
          cff1=             MAX(0.0,                                 
     &                          0.007*                               
     &                          (tau_w(i,j)*rho0-bed_bctr(i,j,1)))+    
     &                          6.0*Ssize(i,j)
          cff2=             MAX(0.0,                                 
     &                          0.007*                               
     &                          (tau_w(i,j)-taucb(i,j))*rho0)+   
     &                          6.0*Ssize(i,j)
          bot_thick(i,j)=MAX(cff1,cff2)
# endif
# ifdef MORPHODYN
           ! Apply morphology factor.
           bot_thick(i,j)=MAX( bot_thick(i,j)*morph_fac,      
     &          bot_thick(i,j))
# endif
# if defined COHESIVE_BED || defined MIXED_BED
!  Find first layer with tc > tb
!  Remember, the critical stresses apply to the TOP of each layer
          Ksed = 1
          bzactv = 0.0
          frt = 0.0 ! CRS
          cff1 = rho0*tau_w(i,j)
          tcb_top = bed_bctr(i,j,1)
          tcb_bot = bed_bctr(i,j,2)
#  if defined MIXED_BED
!         Calc. tau crit for bottom of next layer
!         Update cohesive property of seds in top layer
          cff3 = 0.0
          cff4 = 1.0
          do ised=imud1,imud2
             cff3=cff3+bed_frac(i,j,1,ised)
             cff4=cff4*tau_ce(ised)**bed_frac(i,j,1,ised)
          enddo
          do ised=isand1,isand2
             cff4=cff4*tau_ce(ised)**bed_frac(i,j,1,ised)
          enddo   
          pcoh=min(max((cff3-transN)/(transC-transN)        
     &      ,0.0),1.0)
          tcb_top = pcoh*bed_bctr(i,j,1)+(pcoh-1.0)*cff4*rho0   !BF
#  endif
          if(cff1 .GT. tcb_top)then
!            Calculate tcb_temp for next layer
             tcb_bot = bed_bctr(i,j,Ksed+1)
#  if defined MIXED_BED
!            Recalculate cohesive fraction and mean grain tau crit
!            Note that Ksed is used for grain props, and Ksed+1 for
!            bed tau crit at bottom of layer Ksed
             cff3 = 0.0
             cff4 = 1.0
             do ised=imud1,imud2
                cff3=cff3+bed_frac(i,j,Ksed,ised)
                cff4=cff4*tau_ce(ised)**bed_frac(i,j,Ksed,ised)
             enddo
             do ised=isand1,isand2
                cff4=cff4*tau_ce(ised)**bed_frac(i,j,Ksed,ised)
             enddo   
!            Calculate cohesive behavior and blended tau crit
             pcoh=min(max((cff3-transN)/(transC-transN),    
     &            0.0),1.0)
             tcb_bot =pcoh*bed_bctr(i,j,Ksed+1)+(pcoh-1.0)*cff4*rho0  !BF
#  endif
             do while ( (Ksed.LE.(NLAY-1)) .AND.                        
     &            (cff1 .GT. tcb_bot))
!ALA   Instead of adding entire 2nd layer, add just what is needed
!ALA!               Add entire layer
!ALA                bzactv = bzactv + bed_thick(i,j,Ksed)
!       Add thickness equivalent to difference between cff1 and tcb_bot
                if (cff1.GT.bed_bctr(i,j,Ksed+1)) then
                   bzactv = bzactv + bed_thick(i,j,Ksed)
                   tcb_top = tcb_bot
                   tcb_bot = bed_bctr(i,j,Ksed+1)
                else
                   bzactv = bzactv + MIN(bed_thick(i,j,Ksed),           
     &                      (MAX(0.0,0.007*(cff1-tcb_bot))))
                   tcb_top = cff1
                   tcb_bot = bed_bctr(i,j,Ksed+1)   
                endif
#  if defined MIXED_BED
!               Recalculate cohesive fraction and mean grain tau crit
                cff3 = 0.0
                cff4 = 1.0
                do ised=imud1,imud2
                   cff3=cff3+bed_frac(i,j,Ksed,ised)
                   cff4=cff4*tau_ce(ised)**bed_frac(i,j,Ksed,ised)
                enddo
                do ised=isand1,isand2
                   cff4=cff4*tau_ce(ised)**bed_frac(i,j,Ksed,ised)
                enddo
!               Calculate cohesive behavior and blended tau crit
                pcoh=min(max((cff3-transN)/(transC-transN), 
     &               0.0),1.0)
                tcb_bot = pcoh*bed_bctr(i,j,Ksed+1)+(pcoh-1.0)*cff4 
     &                         *rho0                                 !BF
#  endif
                Ksed = Ksed+1
             enddo
             frt = MAX(0.0,(cff1-tcb_top) ) /                        
     &            MAX( eps,                                             
     &            ( tcb_bot-tcb_top ))
          endif
          bzactv = bzactv+frt*bed_thick(i,j,Ksed)
          bzactv = MAX( bzactv, 6.0*Ssize(i,j) ) !CRS
          bot_thick(i,j)=min(bot_thick(i,j),bzactv) !?CRS
!          bot_thick(i,j)=max(bot_thick(i,j),bzactv) !?CRS
# endif /* defined COHESIVE_BED || defined MIXED_BED */
        enddo I_LOOP_CB
      enddo J_LOOP_CB

        J_LOOP2 : do j=Jstr,Jend
        do i=Istr,Iend

!
!          Calculate net deposition and erosion
           cff=0.0
           cff2=0.0
           do ised=1,NST
              cff=cff+settling_flux(i,j,ised)
              cff2=cff2+ero_flux(i,j,ised)
              dep_mass(i,ised)=0.0
              if ((ero_flux(i,j,ised)-settling_flux(i,j,ised)).lt.      
     &             0.0) then
                 dep_mass(i,ised)=settling_flux(i,j,ised)-              
     &                ero_flux(i,j,ised)
              endif
           enddo
           bot_dnet(i,j)=cff-cff2

           if ( cff-cff2.GT.0.0) then ! NET deposition
              !  Deposition. Determine if we need to create a new bed layer 

# if defined COHESIVE_BED || defined MIXED_BED
!
!  Calculate tau_crit of deposited bed
!
!  Calculate new mass in top layer
             bmz(1) = 0.0
             do ised=1,NST
                bmz(1) = bmz(1)+bed_mass(i,j,1,nnew,ised)
             enddo
             if (NLAY.GT.1) then
!  Average of cff1 and cff2, where
!    cff1 = linear extension of previous tcr slope to new surface
!    cff2 = minimum deposition
                cff1 = bed_bctr(i,j, 1) -                               
     &             bot_dnet(i,j) *                                  
     &             (bed_bctr(i,j,2)-bed_bctr(i,j,1)) /                
     &             (bmz(1)-bot_dnet(i,j))
                cff2 = MAX( rho0*tau_w(i,j) , tcr_min )
                bed_bctr(i,j,1) = MIN(bed_bctr(i,j,1),                
     &                                 MAX( 0.5*(cff1+cff2), cff2 ))
!                cff1 = bot_dnet(i,j)/MAX(bmz(1),eps)
!                cff2 = bed_bctr(i,j,1)+bed_bctr(i,j,2)-2*tcr_min
!                bed_bctr(i,j,1) =bed_bctr(i,j,1) - cff1*cff2 
             else
!  TODO: Not sure what mud tau_crit should be for single-layer bed.
!  Try weighted average of dep and old value
                cff1 = (bed_bctr(i,j,1)*(bmz(1)-bot_dnet(i,j))     
     &                  + cff2*bot_dnet(i,j)) / bmz(1)
                bed_bctr(i,j,1) = MAX( cff1, cff2 )
             endif
# endif
              ! (no test for age here)
              bed_age(i,j,1)=time
              if(bed_thick(i,j,1).gt.                                   
     &             MAX(bot_thick(i,j),newlayer_thick)) then
                 ! Top layer is too thick
                 if (NLAY.gt.2) then
                    if(bed_thick(i,j,2).lt.minlayer_thick) then
                    ! Layer 2 is smaller than minimum size
                    ! Instead of pushing down all layers, just combine top 2 layers
                       cff=0.0
                       cff1=0.0
                       cff2=0.0
                       do ised=1,NST
                          cff =cff +dep_mass(i,ised)
                          cff1=cff1+bed_mass(i,j,1,nnew,ised)
                          cff2=cff2+bed_mass(i,j,2,nnew,ised)
                       enddo
# if defined COHESIVE_BED || defined MIXED_BED
!
!   Assign new tau_crit at 2nd layer
!
                    bed_bctr(i,j,2)= bed_bctr(i,j,2) -                
     &       (cff1-cff)/(cff1+cff2)*(bed_bctr(i,j,3)-bed_bctr(i,j,1))
!                  taucb(i,j,2) = 0.5*( taucb(i,j,1)+         
!     &                                        taucb(i,j,2) )
# endif
!

                       !  Update bed mass
                       do ised=1,NST
                          bed_mass(i,j,2,nnew,ised)=                    
     &                         MAX(bed_mass(i,j,2,nnew,ised)+           
     &                         bed_mass(i,j,1,nnew,ised)-               
     &                         dep_mass(i,ised),0.0)
                          bed_mass(i,j,1,nnew,ised)=dep_mass(i,ised)
                       enddo
                       ! ALA - average time and porosity
                       ! ALA CHECK WITH CRS cff1 or cff1-cff for first layer
                       bed_age(i,j,2)=(bed_age(i,j,1)*cff1+         
     &                         bed_age(i,j,2)*cff2)/(cff1+cff2)
                       bed_age(i,j,1)=time
                       bed_poros(i,j,2)=(bed_poros(i,j,1)*cff1+         
     &                         bed_poros(i,j,2)*cff2)/(cff1+cff2)
                       ! ALA CHECK WITH CRS POROSITY OF 1ST LAYER
                       bed_poros(i,j,1)=bed_poros(i,j,1)
                    else
                    ! Layer 2 is > minlayer thick, need another layer
                    !  Combine bottom layers.
                       cff1=0.0
                       cff2=0.0
                       do ised=1,NST
                          cff1=cff1+bed_mass(i,j,NLAY-1,nnew,ised)
                          cff2=cff2+bed_mass(i,j,NLAY,nnew,ised)
                       enddo
                       bed_poros(i,j,NLAY)=                             
     &                      (bed_poros(i,j,NLAY-1)*cff1+                
     &                      bed_poros(i,j,NLAY)*cff2)/(cff1+cff2)
                       bed_age(i,j,NLAY)=                             
     &                      (bed_age(i,j,NLAY-1)*cff1+                
     &                      bed_age(i,j,NLAY)*cff2)/(cff1+cff2)
# if defined COHESIVE_BED || defined MIXED_BED
!
!   Assign tcrit at top of new bottom bed tcrit for NLAY-1
              bed_bctr(i,j,NLAY)= bed_bctr(i,j,NLAY-1)
!              taucb(i,j,NLAY)=taucb(i,j,NLAY-1) +                &
!     &            cff2*(taucb(i,j,NLAY)-taucb(i,j,NLAY-1))/cff1
# endif

                       do ised=1,NST
                          bed_mass(i,j,NLAY,nnew,ised)=                 
     &                         bed_mass(i,j,NLAY-1,nnew,ised)+          
     &                         bed_mass(i,j,NLAY  ,nnew,ised)
                       enddo
                       !
                       !  Push layers down.
                       do k=NLAY-1,2,-1
                          bed_poros(i,j,k)=bed_poros(i,j,k-1)
                          bed_age(i,j,k)=bed_age(i,j,k-1)
                          do ised =1,NST
                             bed_mass(i,j,k,nnew,ised)=                 
     &                            bed_mass(i,j,k-1,nnew,ised)
                          enddo
# if defined COHESIVE_BED || defined MIXED_BED
                bed_bctr(i,j,k)=bed_bctr(i,j,k-1)
# endif
                       enddo
                       !  Set new top parameters for top 2 layers
# if defined COHESIVE_BED || defined MIXED_BED
                       !  Tau_crit at top already determined
                       !  Set tau_crit at 2nd layer
                       cff=0.0
                       cff1=0.0
                       do ised=1,NST
                          cff =cff +dep_mass(i,ised)
                          cff1=cff1+bed_mass(i,j,1,nnew,ised)
                       enddo

                       cff2 = (bed_bctr(i,j,2)-bed_bctr(i,j,1))/cff1
                       bed_bctr(i,j,2) = bed_bctr(i,j,1)+cff*cff2
# endif

                       do ised=1,NST
                          bed_mass(i,j,2,nnew,ised)=                    
     &                         MAX(bed_mass(i,j,2,nnew,ised)-           
     &                         dep_mass(i,ised),0.0)
                          bed_mass(i,j,1,nnew,ised)=dep_mass(i,ised)
                       enddo
                    endif
                 elseIF (NLAY.eq.2) then 
                  ! NBED=2
                  !
# if defined COHESIVE_BED || defined MIXED_BED
                    !  Tau_crit at top already determined 
                    !  Set tau_crit at 2nd layer
                    cff=0.0
                    cff1=0.0
                    do ised=1,NST
                       cff =cff +dep_mass(i,ised)
                       cff1=cff1+bed_mass(i,j,1,nnew,ised)
                    enddo
                    cff2 = (bed_bctr(i,j,2)-bed_bctr(i,j,1))/cff1
                    bed_bctr(i,j,2) = bed_bctr(i,j,1)+cff*cff2
# endif
                    cff1=0.0
                    cff2=0.0
                    do ised=1,NST
                       cff1=cff1+bed_mass(i,j,1,nnew,ised)
                       cff2=cff2+bed_mass(i,j,2,nnew,ised)
                    enddo
                    do ised=1,NST
                       bed_mass(i,j,2,nnew,ised)=                       
     &                      MAX(bed_mass(i,j,2,nnew,ised)+              
     &                      bed_mass(i,j,1,nnew,ised)-                  
     &                      dep_mass(i,ised),0.0)
                       bed_mass(i,j,1,nnew,ised)=dep_mass(i,ised)
                    enddo
                    ! ALA - average time and porosity
                    bed_age(i,j,2)=(bed_age(i,j,1)*cff1+            
     &                      bed_age(i,j,2)*cff2)/(cff1+cff2)
                    bed_age(i,j,1)=time
                    bed_poros(i,j,2)=(bed_poros(i,j,1)*cff1+            
     &                      bed_poros(i,j,2)*cff2)/(cff1+cff2)
                    ! ALA CHECK WITH CRS POROSITY OF 1ST LAYER
                    bed_poros(i,j,1)=bed_poros(i,j,1)
                 else
                 ! NBED=1
                 endif
              else
                ! Net deposition has occurred, but no new bed layer was created 
             endif
           else
                ! Net erosion occurred
# if defined COHESIVE_BED || defined MIXED_BED
              bmz(1) = 0.0
              do ised=1,NST
                 bmz(1) = bmz(1)+bed_mass(i,j,1,nnew,ised)
              enddo
!  recalc tc for top of bed, based on linear
!  interpolation and mass removed / orig. mass in top layer
              bed_bctr(i,j,1)=bed_bctr(i,j,1)+                        
     &             MIN(1.0,-bot_dnet(i,j)/MAX(eps,bmz(1)))*      
     &             MAX(0.0,(bed_bctr(i,j,2)-bed_bctr(i,j,1)))     
# endif
                bed_age(i,j,1)=time
                if (NLAY.eq.2) then 
                 ! NBED=2
                    do ised=1,NST
                       bed_mass(i,j,2,nnew,ised)=                       
     &                      MAX(bed_mass(i,j,2,nnew,ised)+              
     &                      bed_mass(i,j,1,nnew,ised)-                  
     &                      dep_mass(i,ised),0.0)
                       bed_mass(i,j,1,nnew,ised)=dep_mass(i,ised)
                    enddo
                elseif (NLAY.eq.1) then
                ! ALF NO NEED TO do ANYTHING
                else
                endif
           endif

          ! Recalculate thickness and fractions for all layers.
           do k=1,NLAY
              cff3=0.0
              do ised=1,NST
                 cff3=cff3+bed_mass(i,j,k,nnew,ised)
              enddo
              if (cff3.eq.0.0) then 
                 cff3=eps 
              endif 
              bed_thick(i,j,k)=0.0
              do ised=1,NST
                 bed_frac(i,j,k,ised)=bed_mass(i,j,k,nnew,ised)/cff3
                 bed_thick(i,j,k)=MAX(bed_thick(i,j,k)+                 
     &                         bed_mass(i,j,k,nnew,ised)/               
     &                         (Srho(ised)*                          
     &                (1.0-bed_poros(i,j,k))),0.0)
              enddo
           enddo
        enddo
      enddo J_LOOP2

      J_LOOP3 : do j=Jstr,Jend
        I_LOOP3 : do i=Istr,Iend
          if (bot_thick(i,j).gt.bed_thick(i,j,1)) then
            if (NLAY.eq.1) then
              bot_thick(i,j)=bed_thick(i,j,1)
            else
              thck_to_add=bot_thick(i,j)-bed_thick(i,j,1)
              thck_avail=0.0
              Ksed=1                                        ! initialize
              do k=2,NLAY
                if (thck_avail.lt.thck_to_add) then
                  thck_avail=thck_avail+bed_thick(i,j,k)
                  Ksed=k
                endif
              enddo
!
!  Catch here if there was not enough bed material.
!
              if (thck_avail.lt.thck_to_add) then
                bot_thick(i,j)=bed_thick(i,j,1)+thck_avail
                thck_to_add=thck_avail
              endif
!
!  Update bed mass of top layer and fractional layer.
!
              cff2=MAX(thck_avail-thck_to_add,0.0)/                  
     &             MAX(bed_thick(i,j,Ksed),eps)         
              do ised=1,NST
                cff1=0.0
                do k=1,Ksed
                  cff1=cff1+bed_mass(i,j,k,nnew,ised)
                enddo
                cff3=cff2*bed_mass(i,j,Ksed,nnew,ised)
                bed_mass(i,j,1   ,nnew,ised)=cff1-cff3
                bed_mass(i,j,Ksed,nnew,ised)=cff3
              enddo
!
!  Update thickness of fractional layer ksource_sed.
!
              bed_thick(i,j,Ksed)=MAX(thck_avail-thck_to_add,0.0)

# if defined COHESIVE_BED || defined MIXED_BED
!
! Update tau_cr of fractional layer   !RB debug
!           IF (Ksed < NLAY ) then
            bed_bctr(i,j,Ksed) = bed_bctr(i,j,Ksed+1)-              
     &                 cff2*(bed_bctr(i,j,Ksed+1)-bed_bctr(i,j,Ksed))
!           endif
# endif

!
!  Update bed fraction of top layer.
!
              cff3=0.0
              do ised=1,NST
                cff3=cff3+bed_mass(i,j,1,nnew,ised)
              enddo
              if (cff3.eq.0.0) then
                 cff3=eps
              endif
              do ised=1,NST
                bed_frac(i,j,1,ised)=bed_mass(i,j,1,nnew,ised)/cff3
              enddo
!
!  Upate bed thickness of top layer.
!
              bed_thick(i,j,1)=bot_thick(i,j)
!
!  Pull all layers closer to the surface.
!
              do k=Ksed,NLAY
                ks=Ksed-2
                bed_thick(i,j,k-ks)=bed_thick(i,j,k)
                bed_poros(i,j,k-ks)=bed_poros(i,j,k)
                bed_age(i,j,k-ks)=bed_age(i,j,k)

# if defined COHESIVE_BED || defined MIXED_BED
                bed_bctr(i,j,k-ks)=bed_bctr(i,j,k)
# endif

                do ised=1,NST
                  bed_frac(i,j,k-ks,ised)=bed_frac(i,j,k,ised)
                  bed_mass(i,j,k-ks,nnew,ised)=bed_mass(i,j,k,nnew,ised)
                enddo
              enddo
!
!  Add new layers onto the bottom. Split what was in the bottom layer to
!  fill these new empty cells. ("ks" is the number of new layers).
!
              ks=Ksed-2
              cff=1.0/REAL(ks+1)
# if defined COHESIVE_BED || defined MIXED_BED
#  undef BF_TCR
#  if defined BF_TCR
          bmz(1) = 0.0
          do ised=1,NST
             bmz(1) = bmz(1)+bed_mass(i,j,1,nnew,ised)
          enddo
          do k=2,NLAY
             bmz(k) = bmz(k-1)
             do ised=1,NST
                bmz(k)=bmz(k)+bed_mass(i,j,k,nnew,ised)
             enddo
          enddo
          tcr(1) = tcr_min
          do k=2,NLAY
             tcr(k) = tcr_min
             if (bmz(k-1).GT.eps) then
!                tcr(k) = exp((log(bmz(k-1))-                            
!     &               bottom(i,j,idoff))/                                
!     &               bottom(i,j,idslp))
                tcr(k) = exp((log(bmz(k-1))-                            
     &                    tcr_off)/                                 
     &                    tcr_slp)
             endif
             tcr(k) = MIN( MAX( tcr(k), tcr_min), tcr_max )
          enddo
#  else
!alf              cff1 = cff*(tcr_max-taucb(i,j,NLAY-ks))
!alf          Use the value from the bottom layer as max.
              cff1 = cff*(bed_bctr(i,j,NLAY)-bed_bctr(i,j,NLAY-ks))
#  endif

! Interpolate bottom layer to tau_crit_max
! TODO: should be the reference profile and not tau_crit_max
                  do k=NLAY-1,NLAY-ks,-1
#  if defined BF_TCR
                     bed_bctr(i,j,k)=tcr(k)
#  else
                     bed_bctr(i,j,k)=bed_bctr(i,j,NLAY-ks)+           
     &                       REAL(k-NLAY+ks) * cff1
#  endif
                  enddo
# endif
              ! ALA CHECK WITH CRS about bed_frac 
              nnn=0
              do ised=1,NST
                 nlysm(ised)=newlayer_thick*REAL(ks+1)*          
     &                  (Srho(ised)*                                 
     &                  (1.0-bed_poros(i,j,NLAY-ks)))*               
     &                  bed_frac(i,j,NLAY-ks,ised)
                  IF (ks.gt.0) then
                        if (bed_mass(i,j,NLAY-ks,nnew,ised).gt.         
     &                       nlysm(ised)) then
                           nnn=nnn+1
                           nlysm(ised)=                                 
     &                           newlayer_thick*REAL(ks)*        
     &                          (Srho(ised)*                         
     &                          (1.0-bed_poros(i,j,NLAY-ks)))*       
     &                          bed_frac(i,j,NLAY-ks,ised)
                        endif
                    endif
              enddo
              if (nnn.eq.NST) then
                  bed_thick(i,j,NLAY)=bed_thick(i,j,NLAY-ks)-           
     &                                newlayer_thick*REAL(ks)
                  do ised=1,NST
                     bed_mass(i,j,NLAY,nnew,ised)=                      
     &                  bed_mass(i,j,NLAY-ks,nnew,ised)-nlysm(ised)
                  enddo
                  do k=NLAY-1,NLAY-ks,-1
                     bed_thick(i,j,k)=newlayer_thick
                     bed_age(i,j,k)=bed_age(i,j,NLAY-ks)
                     do ised=1,NST
                        bed_frac(i,j,k,ised)=bed_frac(i,j,NLAY-ks,ised)
                        bed_mass(i,j,k,nnew,ised)=                      
     &                           nlysm(ised)/REAL(ks)
                     enddo
                  enddo
              else
                 cff=1.0/REAL(ks+1)
                  do k=NLAY,NLAY-ks,-1
                     bed_thick(i,j,k)=bed_thick(i,j,NLAY-ks)*cff
                     bed_age(i,j,k)=bed_age(i,j,NLAY-ks)
                     do ised=1,NST
                        bed_frac(i,j,k,ised)=bed_frac(i,j,NLAY-ks,ised)
                        bed_mass(i,j,k,nnew,ised)=                      
     &                             bed_mass(i,j,NLAY-ks,nnew,ised)*cff
                     enddo
                  enddo
              endif
            endif  ! NLAY > 1
          endif  ! increase top bed layer
# if defined MIXED_BED
!
!         Update cohesive property of seds in top layer
          cff1 = 0.0
          do ised=imud1,imud2
             cff1=cff1+bed_frac(i,j,1,ised)
          enddo
          bot_dprp(i,j)=min(max((cff1-transN)/                  
     &     (transC-transN),0.0),1.0)
# endif
          enddo I_LOOP3 

# if defined SED_FLOCS && defined SED_DEFLOC
          I_LOOP_DFL: do i=Istr,Iend
!        Activate defloc in all layers but the top one
!            do k=2,NLAY 
!        Activate defloc in all layers
            do k=1,NLAY
              alpha = MIN((time-bed_age(i,j,k))/t_dfloc,      
     &                  1.0)
              cff3=0.0
              do ised=1,NST
                cff3=cff3+bed_mass(i,j,k,nnew,ised)
              enddo
              cff1 = 0.0
              do ised=imud1,imud2
                cff1 = cff1+bed_mass(i,j,k,nnew,ised)
              enddo
              do ised=imud1,imud2
                indx= imud1-NSAND 
                masseq(indx)=mud_frac_eq(indx)*cff1
                bed_mass(i,j,k,nnew,ised)=alpha*masseq(indx)+           
     &                  bed_mass(i,j,k,nnew,ised)*(1.0-alpha)
                bed_frac(i,j,k,ised)=bed_mass(i,j,k,nnew,ised)/cff3
              enddo
            enddo
          enddo I_LOOP_DFL
# endif

# if defined COHESIVE_BED || defined MIXED_BED
         I_LOOP_CB2: do i=Istr,Iend
!
!  Key to this algorithm: "mass available depth" (mad) [kg/m2]
!  present tau crit profile (tcp)
!  representative tau crit profile (tcr)          
!  Compute depth and mass depth of sediment column
!  Compute cumulative depth
!  Compute mass depth
          bmz(1) = 0.0
          do ised=1,NST
             bmz(1) = bmz(1)+bed_mass(i,j,1,nnew,ised)
          enddo
          do k=2,NLAY
             bmz(k) = bmz(k-1)
             do ised=1,NST
                bmz(k)=bmz(k)+bed_mass(i,j,k,nnew,ised)
             enddo
          enddo

!  Calculate representative critical shear stress profile
!  Note that the values are for the TOP of the layer...we
!  assume the bottom of the bottom layer has tcr = tcr_max
          tcr(1) = tcr_min
          do k=2,NLAY
             tcr(k) = tcr_min
             if (bmz(k-1).GT.eps) then
!                tcr(k) = exp((log(bmz(k-1))-                            &
!     &               bottom(i,j,idoff))/                                &
!     &               bottom(i,j,idslp))
                tcr(k) = exp((log(bmz(k-1))-                            
     &                    tcr_off)/                                 
     &                    tcr_slp)
             endif
             tcr(k) = MIN( MAX( tcr(k), tcr_min), tcr_max )
          enddo
!          ks=Ksed-2
!          do k=NLAY,NLAY-ks,-1
!             taucb(i,j,k)=tcr(k)
!          enddo
!  Relax tau crit profile taucb(i,j,k) 
!  towards representative profile tcr...100 x slower for "swelling"
!  than consolidation
!  TODO - make the factor an input parameter
          if( bed_bctr(i,j,1).LE.tcr(1)) then
!             alpha = MIN(dt/bottom(i,j,idtim),1.0)
             alpha = MIN(dt/tcr_tim,1.0)
          else
!             alpha = MIN(dt/(100.0*bottom(i,j,idtim)),1.0)
             alpha = MIN(dt/(100.0*tcr_tim),1.0)
          endif
          do k=1,NLAY
             bed_bctr(i,j,k)=bed_bctr(i,j,k)+                         
     &            alpha*(tcr(k)-bed_bctr(i,j,k))
             bed_bctr(i,j,k)=                                          
     &          MIN( MAX( bed_bctr(i,j,k), tcr_min), tcr_max )
          enddo
!ALA    Enforce last layer = tcr_max
!ALA          bed_bctr(i,j,NLAY)= tcr_max
        enddo I_LOOP_CB2
# endif /* cohesive bed */

      enddo J_LOOP3
!
!-----------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions to property arrays.
!-----------------------------------------------------------------------
!
! Set boundary conditions
!
# ifndef EW_PERIODIC
      if (EASTERN_EDGE) then
        do j=Jstr,Jend
          do k=1,NLAY
            do ised=1,NST
            bed_mass(Iend+1,j,k,nnew,ised)=
     &                       bed_mass(Iend,j,k,nnew,ised)
            bed_frac(Iend+1,j,k,ised)=
     &                        bed_frac(Iend,j,k,ised)
             enddo
             bed_thick(Iend+1,j,k)=bed_thick(Iend,j,k)
             bed_poros(Iend+1,j,k)=bed_poros(Iend,j,k)
             bed_age(Iend+1,j,k)=bed_age(Iend,j,k)
#  if defined MIXED_BED || defined COHESIVE_BED
             bed_bctr(Iend+1,j,k)=bed_bctr(Iend,j,k)
#  endif
           enddo
           bot_thick(Iend+1,j)=bot_thick(Iend,j)
         enddo
        endif
      if (WESTERN_EDGE) then
        do j=Jstr,Jend
          do k=1,NLAY
            do ised=1,NST
            bed_mass(Istr-1,j,k,nnew,ised)=
     &                       bed_mass(Istr,j,k,nnew,ised)
            bed_frac(Istr-1,j,k,ised)=
     &                        bed_frac(Istr,j,k,ised)
             enddo
             bed_thick(Istr-1,j,k)=bed_thick(Istr,j,k)
             bed_poros(Istr-1,j,k)=bed_poros(Istr,j,k)
             bed_age(Istr-1,j,k)=bed_age(Istr,j,k)
#  if defined MIXED_BED || defined COHESIVE_BED
             bed_bctr(Istr-1,j,k)=bed_bctr(Istr,j,k)
#  endif
           enddo
           bot_thick(Istr-1,j)=bot_thick(Istr,j)
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (NORTHERN_EDGE) then
        do i=Istr,Iend
          do k=1,NLAY
            do ised=1,NST
              bed_mass (i,Jend+1,k,nnew,ised)=
     &                            bed_mass(i,Jend,k,nnew,ised)
              bed_frac (i,Jend+1,k,ised)=
     &                            bed_frac(i,Jend,k,ised)
            enddo
            bed_thick(i,Jend+1,k)=bed_thick(i,Jend,k)
            bed_poros(i,Jend+1,k)=bed_poros(i,Jend,k)
            bed_age  (i,Jend+1,k)=bed_age(i,Jend,k)
#  if defined MIXED_BED || defined COHESIVE_BED
            bed_bctr (i,Jend+1,k)=bed_bctr(i,Jend,k)
#  endif
          enddo
          bot_thick  (i,Jend+1)=bot_thick(i,Jend)
        enddo
      endif
      if (SOUTHERN_EDGE) then
        do i=Istr,Iend
          do k=1,NLAY
            do ised=1,NST
              bed_mass (i,Jstr-1,k,nnew,ised)=
     &                            bed_mass(i,Jstr,k,nnew,ised)
              bed_frac (i,Jstr-1,k,ised)=
     &                            bed_frac(i,Jstr,k,ised)
            enddo
            bed_thick(i,Jstr-1,k)=bed_thick(i,Jstr,k)
            bed_poros(i,Jstr-1,k)=bed_poros(i,Jstr,k)
            bed_age  (i,Jstr-1,k)=bed_age(i,Jend,k)
#  if defined MIXED_BED || defined COHESIVE_BED
            bed_bctr (i,Jstr-1,k)=bed_bctr(i,Jend,k)
#  endif
          enddo
          bot_thick  (i,Jstr-1)=bot_thick(i,Jstr)
        enddo
      endif
# endif
# if !defined EW_PERIODIC && !defined NS_PERIODIC
      if (SOUTHERN_EDGE.and.WESTERN_EDGE) then
        bed_mass(Istr-1,Jstr-1,:,nnew,:)=
     &  0.5*(bed_mass(Istr,Jstr-1,:,nnew,:)
     &      +bed_mass(Istr-1,Jstr,:,nnew,:))
        bed_frac(Istr-1,Jstr-1,:,:)=
     &  0.5*(bed_frac(Istr,Jstr-1,:,:)
     &      +bed_frac(Istr-1,Jstr,:,:))
        bed_thick(Istr-1,Jstr-1,:)=
     &  0.5*(bed_thick(Istr,Jstr-1,:)
     &      +bed_thick(Istr-1,Jstr,:))
        bed_poros(Istr-1,Jstr-1,:)=
     &  0.5*(bed_poros(Istr,Jstr-1,:)
     &      +bed_poros(Istr-1,Jstr,:))
        bed_age(Istr-1,Jstr-1,:)=
     &  0.5*(bed_age(Istr,Jstr-1,:)
     &      +bed_age(Istr-1,Jstr,:))
#  if defined MIXED_BED || defined COHESIVE_BED
        bed_bctr(Istr-1,Jstr-1,:)=
     &  0.5*(bed_bctr(Istr,Jstr-1,:)
     &      +bed_bctr(Istr-1,Jstr,:))
# endif
      endif
      if (SOUTHERN_EDGE.and.EASTERN_EDGE) then
        bed_mass(Iend+1,Jstr-1,:,nnew,:)=
     &  0.5*(bed_mass(Iend,Jstr-1,:,nnew,:)
     &      +bed_mass(Iend+1,Jstr,:,nnew,:))
        bed_frac(Iend+1,Jstr-1,:,:)=
     &  0.5*(bed_frac(Iend,Jstr-1,:,:)
     &      +bed_frac(Iend+1,Jstr,:,:))
        bed_thick(Iend+1,Jstr-1,:)=
     &  0.5*(bed_thick(Iend,Jstr-1,:)
     &      +bed_thick(Iend+1,Jstr,:))
        bed_poros(Iend+1,Jstr-1,:)=
     &  0.5*(bed_poros(Iend,Jstr-1,:)
     &      +bed_poros(Iend+1,Jstr,:))
        bed_age(Iend+1,Jstr-1,:)=
     &  0.5*(bed_age(Iend,Jstr-1,:)
     &      +bed_age(Iend+1,Jstr,:))
#  if defined MIXED_BED || defined COHESIVE_BED
        bed_bctr(Iend+1,Jstr-1,:)=
     &  0.5*(bed_bctr(Iend,Jstr-1,:)
     &      +bed_bctr(Iend+1,Jstr,:))
#  endif
      endif
      if (NORTHERN_EDGE.and.WESTERN_EDGE) then
        bed_mass(Istr-1,Jend+1,:,nnew,:)=
     &  0.5*(bed_mass(Istr,Jend+1,:,nnew,:)
     &      +bed_mass(Istr-1,Jend,:,nnew,:))
        bed_frac(Istr-1,Jend+1,:,:)=
     &  0.5*(bed_frac(Istr,Jend+1,:,:)
     &      +bed_frac(Istr-1,Jend,:,:))
        bed_thick(Istr-1,Jend+1,:)=
     &  0.5*(bed_thick(Istr,Jend+1,:)
     &      +bed_thick(Istr-1,Jend,:))
        bed_poros(Istr-1,Jend+1,:)=
     &  0.5*(bed_poros(Istr,Jend+1,:)
     &      +bed_poros(Istr-1,Jend,:))
        bed_age(Istr-1,Jend+1,:)=
     &  0.5*(bed_age(Istr,Jend+1,:)
     &      +bed_age(Istr-1,Jend,:))
#  if defined MIXED_BED || defined COHESIVE_BED
        bed_bctr(Istr-1,Jend+1,:)=
     &  0.5*(bed_bctr(Istr,Jend+1,:)
     &      +bed_bctr(Istr-1,Jend,:))
#  endif
      endif
      if (NORTHERN_EDGE.and.EASTERN_EDGE) then
        bed_mass(Iend+1,Jend+1,:,nnew,:)=
     &   0.5*(bed_mass(Iend,Jend+1,:,nnew,:)
     &       +bed_mass(Iend+1,Jend,:,nnew,:))
        bed_frac(Iend+1,Jend+1,:,:)=
     &   0.5*(bed_frac(Iend,Jend+1,:,:)
     &       +bed_frac(Iend+1,Jend,:,:))
        bed_thick(Iend+1,Jend+1,:)=
     &   0.5*(bed_thick(Iend,Jend+1,:)
     &       +bed_thick(Iend+1,Jend,:))
        bed_poros(Iend+1,Jend+1,:)=
     &   0.5*(bed_poros(Iend,Jend+1,:)
     &       +bed_poros(Iend+1,Jend,:))
        bed_age(Iend+1,Jend+1,:)=
     &   0.5*(bed_age(Iend,Jend+1,:)
     &       +bed_age(Iend+1,Jend,:))
#  if defined MIXED_BED || defined COHESIVE_BED
        bed_bctr(Iend+1,Jend+1,:)=
     &   0.5*(bed_bctr(Iend,Jend+1,:)
     &       +bed_bctr(Iend+1,Jend,:))
#  endif
      endif
# endif

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      do ised=1,NST
        do k=1,NLAY
          call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 bed_frac(START_2D_ARRAY,k,ised))
          call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 bed_mass(START_2D_ARRAY,k,nnew,ised))
         enddo 
      enddo
      do k=1,NLAY
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 bed_poros(START_2D_ARRAY,k))
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 bed_thick(START_2D_ARRAY,k))
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 bed_age(START_2D_ARRAY,k))
#  if defined MIXED_BED || defined COHESIVE_BED
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 bed_bctr(START_2D_ARRAY,k))
#  endif
      enddo   
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 bot_thick(START_2D_ARRAY))
# endif
!
!-----------------------------------------------------------------------
!  Store total bed thickness
!-----------------------------------------------------------------------
!
# ifdef MORPHODYN
      do j=Jstr,Jend
        do i=Istr,Iend
          bed_thick_tot(i,j,nnew)=0.
          do k=1,NLAY
            bed_thick_tot(i,j,nnew)=bed_thick_tot(i,j,nnew)+
     &                              bed_thick(i,j,k)
          enddo
        enddo
      enddo
!
! Set boundary conditions
!
#  ifndef EW_PERIODIC
      if (EASTERN_EDGE) then
        do j=Jstr,Jend
          bed_thick_tot(Iend+1,j,nnew)=bed_thick_tot(Iend,j,nnew)
         enddo
      endif
      if (WESTERN_EDGE) then
        do j=Jstr,Jend
          bed_thick_tot(Istr-1,j,nnew)=bed_thick_tot(Istr,j,nnew)
        enddo
      endif
#  endif
#  ifndef NS_PERIODIC
      if (NORTHERN_EDGE) then
        do i=Istr,Iend
          bed_thick_tot(i,Jend+1,nnew)=bed_thick_tot(i,Jend,nnew)
        enddo
      endif
      if (SOUTHERN_EDGE) then
        do i=Istr,Iend
          bed_thick_tot(i,Jstr-1,nnew)=bed_thick_tot(i,Jstr,nnew)
        enddo
      endif
#  endif
#  if !defined EW_PERIODIC && !defined NS_PERIODIC
      if (SOUTHERN_EDGE.and.WESTERN_EDGE) then
        bed_thick_tot(Istr-1,Jstr-1,nnew)=
     &  0.5*(bed_thick_tot(Istr,Jstr-1,nnew)
     &      +bed_thick_tot(Istr-1,Jstr,nnew))
      endif
      if (SOUTHERN_EDGE.and.EASTERN_EDGE) then
        bed_thick_tot(Iend+1,Jstr-1,nnew)=
     &  0.5*(bed_thick_tot(Iend,Jstr-1,nnew)
     &      +bed_thick_tot(Iend+1,Jstr,nnew))
      endif
      if (NORTHERN_EDGE.and.WESTERN_EDGE) then
        bed_thick_tot(Istr-1,Jend+1,nnew)=
     &  0.5*(bed_thick_tot(Istr,Jend+1,nnew)
     &      +bed_thick_tot(Istr-1,Jend,nnew))
      endif
      if (NORTHERN_EDGE.and.EASTERN_EDGE) then
        bed_thick_tot(Iend+1,Jend+1,nnew)=
     &   0.5*(bed_thick_tot(Iend,Jend+1,nnew)
     &       +bed_thick_tot(Iend+1,Jend,nnew))
      endif
#  endif
!
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
       call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                         bed_thick_tot(START_2D_ARRAY,nnew))
#  endif
!
!-----------------------------------------------------------------------
!  Store bed evolution
!-----------------------------------------------------------------------
!
      do j=Jstr-1,Jend+1
        do i=Istr-1,Iend+1
          dh(i,j)=bed_thick_tot(i,j,nstp)-bed_thick_tot(i,j,nnew)
        enddo
      enddo
# endif
      return
      end

# endif
!=======================================================================
!                                                                      !
!  This routine computed sediment surface layer (sediment-water        !
!  interface) properties.                                              !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
!
      subroutine sed_surf_tile (Istr,Iend,Jstr,Jend)
!
      implicit none
# include "param.h"
# include "bbl.h"
# include "grid.h"
# include "scalars.h"
# include "sediment.h"
      integer Istr,Iend,Jstr,Jend
      integer i,j,ised
      real cff,cff1,cff2,cff3,cff4,eps
      parameter (eps=1.D-20)
#  if defined BED_HIDEXP || defined BED_ARMOR
      integer ised0
      real ph,pe, phik,phim,phis, invlog2
#  endif
!
# include "compute_auxiliary_bounds.h"

# define I_EXT_RANGE Istr-1,Iend+1
# define J_EXT_RANGE Jstr-1,Jend+1
!
!-----------------------------------------------------------------------
! Compute Hiding/Exposure and Armoring
!
! Hiding /Exposure: 
! Resuspension of smaller size classes may be reduced due to hiding, 
! and the coarser grains resuspend more easily due to exposure. 
!
! Armoring:
! In addition, bed erosion may become limited when selective entrainment 
! of the fine fraction causes "bed armoring": the remaining coarser 
! grains inhibit resuspension of the fine grains.
!-----------------------------------------------------------------------
!
#  ifdef BED_HIDEXP
#   ifdef BED_HIDEXP_WULIN
!
! Modulate critical stress tau_ce with hiding/exposure effects 
! using hidden and exposed probabilities of each sediment particle 
! in the bed material (Wu & Lin 2014)
!
      do j=J_EXT_RANGE
        do i=I_EXT_RANGE
          do ised=1,NST
            ph=0.
            pe=0.
            do ised0=1,NST
              cff =1./(Sd(ised)+Sd(ised0))
              ph=ph+bed_frac(i,j,1,ised0)*cff*Sd(ised0)  ! prob hiding
              pe=pe+bed_frac(i,j,1,ised0)*cff*Sd(ised )  ! prob exposure
            enddo
            tau_ce_2d(i,j,ised)=tau_ce(ised)*(pe/ph)**(-0.6)
          enddo
        enddo
      enddo
#   else
!
!  Hiding/exposure factor is function of sediment size:
!  (d(ised)/dm)^(-gamma), e.g., Wilcock & Crowe (2003)
!
      invlog2=1./log(2.)
      do j=J_EXT_RANGE
        do i=I_EXT_RANGE
          phim=0.
          do ised=1,NST
            phik=-log(Sd(ised))*invlog2
            phim=phim+phik*bed_frac(i,j,1,ised)
          enddo
          cff=2.**(phim)  ! 1/mean(Sd)
          do ised=1,NST
            cff1=1-0.67/(1+exp(1.5-Sd(ised)*cff)) ! or 1.
            tau_ce_2d(i,j,ised)=tau_ce(ised)*(Sd(ised)*cff)**(-cff1)
          enddo
        enddo
      enddo
#   endif
#  else
      do ised=1,NST
        do j=J_EXT_RANGE
          do i=I_EXT_RANGE
            tau_ce_2d(i,j,ised)=tau_ce(ised)
          enddo
        enddo
      enddo
#  endif
!
#  ifdef BED_ARMOR
!
! Bed erosion may become limited when selective entrainment of the fine 
! fraction causes "bed armoring": the remaining coarser grains inhibit 
! resuspension of the fine grains.
! Formulation in Blaas et al. (2007) from Garcia and Parker (1991) 
!
      invlog2=1./log(2.)
      do j=J_EXT_RANGE
        do i=I_EXT_RANGE
          phim=0.
          phis=0.
          do ised=1,NST
            phik=-log(Sd(ised))*invlog2
            phim=phim+phik*bed_frac(i,j,1,ised) ! mean
          enddo
          do ised=1,NST
            phik=-log(Sd(ised))*invlog2
            phis=phis+((phik-phim)**2)*bed_frac(i,j,1,ised) ! std
          enddo
          phis=SQRT(phis)
          cff=1./(1.-0.29*phis)**5. ! 1/lambda**5
          do ised=1,NST
            tau_ce_2d(i,j,ised)=tau_ce_2d(i,j,ised)*cff
          enddo
        enddo
      enddo
#  endif

#  undef DEBUG_ARMOR
#  ifdef DEBUG_ARMOR
      do ised=1,NST
        MPI_master_only write(*,*) '================================'
        MPI_master_only write(*,*) ' ised      :',ised
        MPI_master_only write(*,*) ' Sd        :',1.e3*Sd(ised)
        MPI_master_only write(*,*) ' bed_frac  :',bed_frac(25,1,1,ised)
        MPI_master_only write(*,*) ' tau_ce_2d :',rho0
     &                                            *tau_ce_2d(25,1,ised)
        MPI_master_only write(*,*) '================================'
      enddo
#  endif
!
!-----------------------------------------------------------------------
!  Update mean surface properties.
!  Sd50 must be positive definite, due to BBL routines.
!  Srho must be >1000, due to (s-1) in BBL routines
!-----------------------------------------------------------------------
!
      do j=Jstr-1,Jend+1
        do i=Istr-1,Iend+1
          cff1=1.
          cff2=1.
          cff3=1.
          cff4=1.
          do ised=1,NST
            cff1=cff1*tau_ce_2d(i,j,ised)**bed_frac(i,j,1,ised)
            cff2=cff2*Sd(ised)**bed_frac(i,j,1,ised)
            cff3=cff3*(wsed(ised)+eps)**bed_frac(i,j,1,ised)
            cff4=cff4*Srho(ised)**bed_frac(i,j,1,ised)
          enddo
          taucb(i,j)=cff1
          Ssize(i,j)=cff2
          w_set(i,j)=cff3
          Sdens(i,j)=MAX(cff4,1050.)
        enddo
      enddo
!
! Set boundary conditions
!
# ifndef EW_PERIODIC
      if (EASTERN_EDGE) then
        do j=Jstr,Jend
           taucb(Iend+1,j)=taucb(Iend,j)
           Ssize(Iend+1,j)=Ssize(Iend,j)
           w_set(Iend+1,j)=w_set(Iend,j)
           Sdens(Iend+1,j)=Sdens(Iend,j)
        enddo
      endif
      if (WESTERN_EDGE) then
        do j=Jstr,Jend
           taucb(Istr-1,j)=taucb(Istr,j)
           Ssize(Istr-1,j)=Ssize(Istr,j)
           w_set(Istr-1,j)=w_set(Istr,j)
           Sdens(Istr-1,j)=Sdens(Istr,j)
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (NORTHERN_EDGE) then
        do i=Istr,Iend
          taucb  (i,Jend+1)=taucb(i,Jend)
          Ssize  (i,Jend+1)=Ssize(i,Jend)
          w_set  (i,Jend+1)=w_set(i,Jend)
          Sdens  (i,Jend+1)=Sdens(i,Jend)
        enddo
      endif
      if (SOUTHERN_EDGE) then
        do i=Istr,Iend
          taucb  (i,Jstr-1)=taucb(i,Jstr)
          Ssize  (i,Jstr-1)=Ssize(i,Jstr)
          w_set  (i,Jstr-1)=w_set(i,Jstr)
          Sdens  (i,Jstr-1)=Sdens(i,Jstr)
        enddo
      endif
# endif
# if !defined EW_PERIODIC && !defined NS_PERIODIC
      if (SOUTHERN_EDGE.and.WESTERN_EDGE) then
        taucb(Istr-1,Jstr-1)=0.5*(taucb(Istr,Jstr-1)+taucb(Istr-1,Jstr))
        Ssize(Istr-1,Jstr-1)=0.5*(Ssize(Istr,Jstr-1)+Ssize(Istr-1,Jstr))
        w_set(Istr-1,Jstr-1)=0.5*(w_set(Istr,Jstr-1)+w_set(Istr-1,Jstr))
        Sdens(Istr-1,Jstr-1)=0.5*(Sdens(Istr,Jstr-1)+Sdens(Istr-1,Jstr))
      endif
      if (SOUTHERN_EDGE.and.EASTERN_EDGE) then
        taucb(Iend+1,Jstr-1)=0.5*(taucb(Iend,Jstr-1)+taucb(Iend+1,Jstr))
        Ssize(Iend+1,Jstr-1)=0.5*(Ssize(Iend,Jstr-1)+Ssize(Iend+1,Jstr))
        w_set(Iend+1,Jstr-1)=0.5*(w_set(Iend,Jstr-1)+w_set(Iend+1,Jstr))
        Sdens(Iend+1,Jstr-1)=0.5*(Sdens(Iend,Jstr-1)+Sdens(Iend+1,Jstr))
      endif
      if (NORTHERN_EDGE.and.WESTERN_EDGE) then
        taucb(Istr-1,Jend+1)=0.5*(taucb(Istr,Jend+1)+taucb(Istr-1,Jend))
        Ssize(Istr-1,Jend+1)=0.5*(Ssize(Istr,Jend+1)+Ssize(Istr-1,Jend))
        w_set(Istr-1,Jend+1)=0.5*(w_set(Istr,Jend+1)+w_set(Istr-1,Jend))
        Sdens(Istr-1,Jend+1)=0.5*(Sdens(Istr,Jend+1)+Sdens(Istr-1,Jend))
      endif
      if (NORTHERN_EDGE.and.EASTERN_EDGE) then
        taucb(Iend+1,Jend+1)=0.5*(taucb(Iend,Jend+1)+taucb(Iend+1,Jend))
        Ssize(Iend+1,Jend+1)=0.5*(Ssize(Iend,Jend+1)+Ssize(Iend+1,Jend))
        w_set(Iend+1,Jend+1)=0.5*(w_set(Iend,Jend+1)+w_set(Iend+1,Jend))
        Sdens(Iend+1,Jend+1)=0.5*(Sdens(Iend,Jend+1)+Sdens(Iend+1,Jend))
      endif
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 taucb(START_2D_ARRAY))
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 Ssize(START_2D_ARRAY))
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 w_set(START_2D_ARRAY))
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                 Sdens(START_2D_ARRAY))
# endif

#  undef I_EXT_RANGE
#  undef J_EXT_RANGE

      return
      end
!                 
!====================================================================
!
!  This routine reads in and processes sediment parameters and
!  some initial data from sediments.in file.
!
!=====================================================================
!
       subroutine init_sediment
# ifdef SED_FLOCS
        USE flocmod, ONLY : flocmod_alloc, flocmod_init
        USE flocmod, ONLY : f_ws
# endif
      implicit none
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "ncscrum.h"
# include "scalars.h"
# include "sediment.h"
!
      integer Nclass, i, icard, iunit, lstr, nl
      integer lenstr
      real rdum, cff1,cff2
      logical ldum
      parameter (iunit=50)
# ifdef SED_FLOCS
      real :: f_mneg_param, f_fp, f_fy, f_cfcst, f_dmin_frag
# endif
!
!---------------------------------------------------------------------
!  Read in initial float locations.
!---------------------------------------------------------------------
!
      lstr=lenstr(sedname)

      open(iunit,file=sedname(1:lstr),form='formatted',
     &  status='old', err=195)
!
!  Read input parameters according to their input card number.
!
      icard=0
      do while (icard.lt.99)
!
!  Read in sediment case title.
!
        if (icard.eq.1) then
          read(iunit,'(a)',err=60) Stitle
          lstr=lenstr(Stitle)
          write(stdout,10) Stitle(1:lstr)
  10      format(1x,'(',a,')',/)
!
!  Read in parameter and initial values per size class.
!
        elseif (icard.eq.2) then
          i=0
          Nclass=0
          do i=1,NST
           read(iunit,*,err=40) Sd(i), Csed(i), Srho(i), Wsed(i),
     &                          Erate(i), tau_ce(i), tau_cd(i), 
     &                          (Bfr(nl,i),nl=1,NLAY)

!
! Re-compute critical shear stress tau_ce (Wu & Lin, 2014)
!
# ifdef TAU_CRIT_SHIELDS
           cff1=10*Sd(i)*(g*(Srho(i)/rho0 -1.))**0.33        ! Sd_star (m)
           cff2=0.3/(1.+1.2*cff1)+0.055*(1.-exp(-0.02*cff1)) ! Critical shield
           tau_ce(i)=cff2*Sd(i)/1000*g*(Srho(i)-rho0)        ! tau_ce (N/m2)
# elif defined TAU_CRIT_WULIN
           cff2=0.03                                         ! Critical shield
           tau_ce(i)=cff2*Sd(i)/1000*g*(Srho(i)-rho0)        ! tau_ce (N/m2)
# endif
     
           MPI_master_only write(stdout,'(/A,2x,i2)')
     &     'Sediment parameters for grain-size class:', i
     
           MPI_master_only write(stdout,'(7(/f10.5,2x,A))')
     &     Sd(i),    'Sd      Median sediment grain diameter (mm).',
     &     Csed(i),  'Csed    Analyt. init. of sed. conc. (kg/m3).',
     &     Srho(i),  'Srho    Sediment grain density (kg/m3).',
     &     Wsed(i),  'Wsed    Particle settling velocity (mm/s).',
     &     Erate(i), 'Erate   Surface erosion rate (kg/(m2s)).',
     &     tau_ce(i),'tau_ce  Critical shear for erosion (N/m2).',
     &     tau_cd(i),'tau_cd  Critical shear for deposition (N/m2).'
          
           do nl=1,NLAY
            MPI_master_only write(stdout,'(f10.5,2x,A,2x,i2)')
     &      Bfr(nl,i), 'bed_frac   Volume fraction, layer:',nl
           enddo
 
           Sd(i)=Sd(i)/1000.         ! mm --> m
           Wsed(i)=Wsed(i)/1000.     ! mm/s --> m/s
           tau_ce(i)=tau_ce(i)/rho0  ! N/m2 --> m2/s2
           tau_cd(i)=tau_cd(i)/rho0  ! N/m2 --> m2/s2
           Nclass=Nclass+1
          enddo

        elseif (icard.eq.3) then
           MPI_master_only write(stdout,'(/A)')
     &     'Bed parameters for all bed levels'

          read(iunit,*,err=60) (Bthk(nl),nl=1,NLAY)
           do nl=1,NLAY
            MPI_master_only write(stdout,'(f10.5,2x,A,2x,i2,2x,A)')
     &      Bthk(nl), 'Bthk   Init. bed thickness, layer:',nl,'(m).'
           enddo

        elseif (icard.eq.4) then

          read(iunit,*,err=60) (Bpor(nl),nl=1,NLAY)
           do nl=1,NLAY
            MPI_master_only write(stdout,'(f10.5,2x,A,2x,i2,2x,A)')
     &      Bpor(nl), 'Bpor   Init. bed porosity, layer:',nl,'(m).'
           enddo

        elseif (icard.eq.5) then
          read(iunit,*,err=60) Hrip                     
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      Hrip, 'Hrip   Init. ripple roughness (m).'
        
        elseif (icard.eq.6) then
          read(iunit,*,err=60) Lrip                     
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      Lrip, 'Lrip   Init. ripple height (m).'

        elseif (icard.eq.7) then
          read(iunit,*,err=60) bedload_coeff                     
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      bedload_coeff, 'bedload_coeff   Bedload rate.'

        elseif (icard.eq.8) then
          read(iunit,*,err=60) morph_fac                     
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      morph_fac, 'morph_fac   Morphological factor.'

# ifdef MIXED_BED
        elseif (icard.eq.9) then
          read(iunit,*,err=60) transC                     
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      transC, 'transC   Cohesive transition.'

        elseif (icard.eq.10) then
          read(iunit,*,err=60) transN                     
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      transN, 'transN   Non cohesive transition.'
# else
        elseif (icard.eq.9) then
          read(iunit,*,err=60) rdum                     
        elseif (icard.eq.10) then
          read(iunit,*,err=60) rdum                     
# endif
# if defined COHESIVE_BED || defined MIXED_BED
        elseif (icard.eq.11) then
          read(iunit,*,err=60) tcr_min                     
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      tcr_min, 'tcr_min   Mud min taucr.'

        elseif (icard.eq.12) then
          read(iunit,*,err=60) tcr_max                     
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      tcr_max, 'tcr_max   Mud max taucr.'

        elseif (icard.eq.13) then
          read(iunit,*,err=60) tcr_slp                     
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      tcr_slp, 'tcr_slp   Mud taucr slope.'

        elseif (icard.eq.14) then
          read(iunit,*,err=60) tcr_off                     
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      tcr_off, 'tcr_off   Mud profile offset.'

        elseif (icard.eq.15) then
          read(iunit,*,err=60) tcr_tim                     
            MPI_master_only write(stdout,'(f10.1,2x,A)')
     &      tcr_tim, 'tcr_tim   Mud consolidation rate.'

# endif
# ifdef SED_FLOCS
        elseif (icard.eq.16) then
          read(iunit,*,err=60) l_ADS,l_ASH,l_COLLFRAG,l_testcase
            MPI_master_only write(stdout,'(L10,2x,A)')
     &      l_ADS, 'l_ADS   Mud differential settling.'
            MPI_master_only write(stdout,'(L10,2x,A)')
     &      l_ASH, 'l_ASH   Mud shear aggregation.'
            MPI_master_only write(stdout,'(L10,2x,A,1x,A)')
     &      l_COLLFRAG, 'l_COLLFRAG   Mud collision-induced',
     &                   'fragmentation.'
            MPI_master_only write(stdout,'(L10,2x,A,1x,A)')
     &      l_testcase, 'l_testcase   Mud boolean for values', 
     &                  'from Verney et al.'

        elseif (icard.eq.17) then
          read(iunit,*,err=60) f_dp0  ,f_nf   ,f_dmax,f_nb_frag,
     &                         f_alpha,f_beta ,f_ater,f_ero_frac,
     &                         f_ero_nbfrag,f_collfragparam,
     &                         f_clim,f_ero_iv                       
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      f_dp0, 'f_dp0   Mud Primary particle size(m).'
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      f_nf, 'f_nf   Mud floc fractal dimension.'
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      f_dmax, 'f_dmax   Mud maximum diameter.'
            MPI_master_only write(stdout,'(f10.5,2x,A,1x,A)')
     &      f_nb_frag, 'f_nb_frag   Mud Number of fragments by shear',
     &                             'erosion.'
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      f_alpha, 'f_alpha   Mud flocculation efficiency.'
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      f_beta, 'f_beta   Mud shear fragmentation rate.'
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      f_ater, 'f_ater   Mud ternary breakup.'
            MPI_master_only write(stdout,'(f10.5,2x,A,1x,A)')
     &      f_ero_frac, 'f_ero_frac   Mud fraction of shear',
     &             'fragmentation term transfered to shear erosion.'
            MPI_master_only write(stdout,'(f10.5,2x,A,1x,A)')
     &      f_ero_nbfrag, 'f_ero_nbfrag   Mud Number of fragments',
     &             'induced by shear erosion.'
            MPI_master_only write(stdout,'(f10.5,2x,A,1x,A)')
     &      f_collfragparam, 'f_collfragparam   Mud Fragmentation rate',
     &             'for collision-induced breakup.'
            MPI_master_only write(stdout,'(f10.5,2x,A,1x,A)')
     &      f_clim, 'f_clim   Mud Min concentration below which',
     &           'flocculation processes are not calc..'
            MPI_master_only write(stdout,'(I10,2x,A)')
     &      f_ero_iv, 'f_ero_iv   Mud Fragment size class.'

#  ifdef SED_DEFLOC
        elseif (icard.eq.18) then
          read(iunit,*,err=60) (mud_frac_eq(nl),nl=1,NMUD)
           do nl=1,NMUD
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      mud_frac_eq(nl),'mud_frac_eq   Equilibrium fractional class'
           enddo
        elseif (icard.eq.19) then
          read(iunit,*,err=60) t_dfloc                     
            MPI_master_only write(stdout,'(f10.5,2x,A)')
     &      t_dfloc, 't_dfloc   Floc decomposition rate.'
#  else
        elseif (icard.eq.18) then
          read(iunit,*,err=60) rdum
        elseif (icard.eq.19) then
          read(iunit,*,err=60) rdum
#  endif
# else
        elseif (icard.eq.16) then
          read(iunit,*,err=60) ldum                     
        elseif (icard.eq.17) then
          read(iunit,*,err=60) rdum                   
        elseif (icard.eq.18) then
          read(iunit,*,err=60) rdum
        elseif (icard.eq.19) then
          read(iunit,*,err=60) rdum                     
# endif
        endif ! icard 
!
!  Read until last input card ID.
!
        read(iunit,*,err=60) icard
      enddo
      goto 90
!
!  Error while reading input parameters.
!
  40  write(stdout,50) icard, i, Nclass, sedname
  50  format(/,' INIT_SEDIMENTS - error reading variables in card: ',
     &             i2, ', entry: ',i3,/,15x,
     &             'nclass:',i3, 'input script: ',a)
  
  60  write(stdout,80) icard, sedname
  80  format(/,' INIT_SEDIMENTS - error while reading input card: ',
     &       i2,15x,'in input script: ',a)

 195  write(stdout,205) sedname
 205  format(/,'sediment inputfile : ',A,/, 
     &   ' inputfile not found or bad sediment initialization')
      stop

  90  close(iunit)
      write(stdout,100) Nclass
 100  format(/,
     &       'Number of size classes in sediment computation:'i2)


       do i=1,NST
         tau_ce_2d(:,:,i)=tau_ce(i)
       enddo
# ifdef SED_FLOCS
        ! some parameters are not read in sediment.in
        ! f_mneg_param, f_fp, f_fy, f_cfcst, f_dmin_frag
        ! values from paraMUSTANG_default.txt
        ! && f_dmax not used but read in sediment.in
        f_mneg_param=0.001
        f_cfcst = 0.1875
        f_fp = 0.1
        f_fy = 1.0e-10 
        f_dmin_frag = 0.00001
        CALL flocmod_alloc(NMUD)
        CALL flocmod_init(l_ADS, l_ASH, l_COLLFRAG,             
     &      f_dp0, f_nf, f_nb_frag, f_alpha, f_beta, f_ater,    
     &      f_ero_frac, f_ero_nbfrag, f_ero_iv, f_mneg_param,   
     &      f_collfragparam, f_dmin_frag, f_cfcst, f_fp, f_fy,  
     &      f_clim, Sd(imud1:imud1+NMUD-1) , 
     &      Srho(imud1:imud1+NMUD-1),      
     &      rho0, l_testcase, stdout)
        Wsed(imud1:imud1+NMUD-1) = f_ws(1:NMUD)
# endif
      return
      end

# ifdef BEDLOAD_VANDERA
!======================================================================!
!======================================================================!
!                          VANDERA FUNCTIONS                           !
!                                                                      !
!   Subroutines and functions required for the sediment bedload        !
!   calculations using Van der A's formulations.                       !  
!                                                                      !
! Tarandeep S. Kalra; Chris Sherwood; John C. Warner                   !
! ... adapted by P. Marchesiello (2021) for CROCO
!======================================================================!
!======================================================================!
!
      subroutine sandload_vandera( wavecycle,
     &                                   Hs, Td,  depth, RR,
     &                              d50, d50_mix, rhos, c_w,
     &                                             eta, dsf,
     &                       T_i, T_iu, uhat_i, mag_theta_i,
     &                                         om_ii, om_iy )
!
      implicit none
#  include "param.h"
#  include "scalars.h"
      real wavecycle
      real Hs, Td, depth, RR
      real d50, d50_mix, rhos, c_w
      real eta, dsf
      real T_i, T_iu
      real uhat_i, mag_theta_i
      real om_ii, om_iy
!
! local variables
! 
      real m_fac,n_fac
      parameter(m_fac=11., n_fac=1.2)
      real alpha_fac
      parameter(alpha_fac=8.2)
      real xi
      parameter(xi=1.7) ! Based on Santoss_core.m
      real eps
      parameter(eps=1.0E-14)
      real eps_eff, P
      real om_i, theta_diff, theta_ieff, theta_cr
      real w_s, what, w_sc, delta
      real cff, cff1
      real w_s_calc, what_calc, theta_cr_calc
!
! Layer thickness
!
      if (eta.gt.0.) then
        delta=eta   ! For ripple regime
      else
        delta=dsf   ! For sheet flow regime
      endif
! 
! Find settling velocity based on Soulsby (1997). 
! VA2013 Use 0.8*d50 for settling velocity (text under equation 28).
!
      w_s=w_s_calc(0.8*d50, rhos)
!
! Correction by vertical advection (Stokes 2nd-order)
!    wavecycle=+1 => VA2013 Eq 29, for crest cycle
!    wavecycle=-1 => VA2013 Eq 30, for trough cycle
!
      what=what_calc(Hs, Td, depth, RR, delta)
      w_sc=MAX(w_s+wavecycle*what,eps)
!
! VA2013 Eq 27-28, Phase lag parameter
!
      cff=(1.-wavecycle*xi*uhat_i)/c_w
      cff1=1./MAX(2.*(T_i-T_iu)*w_sc, eps)
      P=alpha_fac*delta*cff*cff1 
!
! Hiding/exposure correction factor
!
!      eps_eff=(d50/d50_mix)**0.25 
      eps_eff=1. 
!
! CRS for multiple sed types
!
      theta_ieff=eps_eff*mag_theta_i
! 
! Find critical Shields parameters based on Soulsby (1997).
!
      theta_cr=theta_cr_calc(d50, rhos)
!
! VA2013 Eq 2: Sandload entrained during each half-cycle
!
      theta_diff=MAX((theta_ieff-theta_cr),0.)
      om_i=m_fac*(theta_diff)**n_fac
!
! VA2013 Eq 23-26: Sandload transported within the 
!        half-cycle (om_ii) or in the next (om_iy) 
! 
      if (P.lt.1) then
        om_ii=om_i
        om_iy=0.
      else
        om_ii=om_i/P
        om_iy=om_i-om_ii
      endif
!
      return
      end
!
!======================================================================
!
      subroutine full_wave_cycle_stress_factors( d50, d90, osmgd,
     &                                            Td, c_w, depth,
     &                                               phi_curwave,
     &                                            RR, uhat, ahat,
     &                                                umax, umin,
     &                                                 mag_bstrc,
     &                                      T_c, T_t, T_cu, T_tu,
     &                                                ksd, ustrc,
     &                                          delw, udelta, fd,
     &                                   alpha, eta, ksw, tau_wRe )
!
! Iterative solution to obtain current and wave related bed roughness,
! computed from converged time-averaged Shields parameter (VA2013 
! Apendix A). Also, compute full-cycle friction factor and wave Reynolds 
! stress
!
!  This subroutine returns the following:
!  eta                 : ripple height
!  udelta              : current velocity at the wave boundary layer
!  fd                  : current friction factor  
!  tau_wRe             : Wave averaged Reynolds stress
!  T_c, T_t, T_cu, T_tu: Updated time periods in half cycles 
!                        based on current velocity
! 
      implicit none
#  include "param.h"
#  include "scalars.h"
!
!  Imported variable declarations.
!  
! Input the crest or trough half cycle velocity
! d50 -- grain size in meters
! Different for crest and trough half cycles
!       
      real d50, d90, osmgd
      real Td, c_w, depth
      real phi_curwave
      real RR, uhat, ahat
      real umax, umin
      real mag_bstrc
      real T_c, T_t, T_cu, T_tu
      real udelta, delw, fd
      real alpha, eta, ksw, tau_wRe
!
!  Local variables
! 
      integer iter
      integer total_iters
      parameter (total_iters=15)
      real tol,von_k
      parameter(tol=0.001, von_k=0.41)
      real eps
      parameter(eps=1.0E-14)
      real crs_fac
      parameter(crs_fac=1.)
      real theta_avg_old, theta_avg, theta_hat_i
      real psi  ! maximum mobility number
      real rlen ! ripple length
      real ksd, ustrc
      real fw
      real alpha_w, fwd
      real ustarw
      real ksw_calc, ksd_calc, fw_calc
      real fd_calc_madsen, fd_calc_santoss, mu_calc
!
! Maximum mobility number at crest and trough
! For irregular waves, use Rayleigh distributed maximum value
! - VA2013, text under equation Appendix B.4
!
      psi=(1.27*uhat)**2*osmgd
!
! Use Appendix B eqn B.1 and B.2 to get ripple height and length
!
      call ripple_dim(psi, d50, eta, rlen)
!
      eta=eta*ahat
      rlen=rlen*ahat
!
! Initiliaze with theta_avg=0 and theta_hat_i=theta_avg
!
      theta_avg=0.
      theta_avg_old=0.
!
! Calculate friction factor fd
!
#  ifdef BEDLOAD_VANDERA_FD_MADSEN
      fd=fd_calc_madsen(udelta, mag_bstrc)
#  elif defined BEDLOAD_VANDERA_ZEROCURR
      fd=0.
#  endif
!
      do iter=1,total_iters !---- ! Iterative solution
!
! Wave related bed roughness - VA2013 A.5
!
        ksw=ksw_calc(d50, mu_calc(d50), theta_avg, eta, rlen)
!
! Full-cycle wave friction factor - VA2013 Appendix Eq A.4
!
        fw=fw_calc(ahat, ksw)
!
#  if defined BEDLOAD_VANDERA_FD_SANTOSS &&\
     !defined BEDLOAD_VANDERA_ZEROCURR
!
! Current-related bed roughness - VA2013 Appendix A.1
!
        ksd=ksd_calc(d50, d90, mu_calc(d50), theta_avg, eta, rlen)
!
! Friction factor from ksd, udelta, delta
!
        fd=fd_calc_santoss(udelta, delw, ksd)
!
#  endif
!
! Time-averaged absolute Shields stress - VA2013 Appendix Eq. A.3
!
#  ifdef BEDLOAD_VANDERA_ZEROCURR
        theta_avg=osmgd*0.25*fw*uhat**2
#  else
        theta_avg=osmgd*(0.50*fd*udelta**2.+
     &                   0.25*fw*uhat**2.)
#  endif
!
        if (ABS(theta_avg-theta_avg_old).lt.tol) then
          EXIT
        endif
        theta_avg_old=theta_avg
!
      enddo  ! end iterations ----
!
! Full-cycle current friction factor - VA2013 Eq 20
!
#  ifdef BEDLOAD_VANDERA_ZEROCURR
      fwd=fw
      alpha=0.
#  else
      alpha=udelta/(udelta+uhat)
      fwd=alpha*fd+(1.-alpha)*fw
#  endif
!
! Wave Reynolds stress - VA2013 Eq 22
! (contribution associated with progressive surface waves)
!
#  ifdef BEDLOAD_VANDERA_WAVE_AVG_STRESS
      alpha_w=0.424
      tau_wRe=MAX((rho0*fwd*alpha_w*uhat**3./(2.*c_w)),eps)
#  else 
      tau_wRe=0.
#  endif
!        
! Compute the change in time period based on converged udelta 
! (current velocity at wave boundary layer)
!
#  ifndef BEDLOAD_VANDERA_ZEROCURR
      call current_timeperiod(udelta, phi_curwave, umax, umin, RR,
     &                        T_c, T_t, Td)
#  endif
!
! Calculate the effect of surface waves 
!
#  ifdef BEDLOAD_VANDERA_SURFACE_WAVE
      call surface_wave_mod(Td, depth, uhat, T_c, T_cu, T_t, T_tu)
#  endif 
!
      return
      end
!
!======================================================================
!
      subroutine half_wave_cycle_stress_factors( T_iu, T_i, ahat, ksw,
     &                                                      fd, alpha,
     &                                                 alphac, alphaw,
     &                                                     d50, osmgd,
     &                              ui_r, uhat_i, udelta, phi_curwave,
     &                                                   tau_wRe, dsf,
     &                                theta_ix, theta_iy, mag_theta_i )
!
!  This subroutine returns the following:
!  dsf                 : sheetflow thickness
!  theta_ix, theta_iy  : Shields parameter in x and y dir. 
!  mag_theta_i         : Magnitude of Shields parameter for half cycle
! 
      implicit none
#  include "param.h"
#  include "scalars.h"
!  
! Input the crest or trough half cycle velocity
! d50 -- grain size in meters
! Different for crest and trough half cycles 
!       
      real T_iu, T_i, ahat, ksw
      real fd, alpha
      real alphac, alphaw
      real d50, osmgd
      real ui_r, uhat_i, udelta, phi_curwave
      real tau_wRe
      real dsf, theta_ix, theta_iy, mag_theta_i
      real fwi_calc, dsf_calc
!
!  Local variables
! 
      real eps
      parameter(eps = 1.0E-14)
      real fw_i, fwd_i
      real alpha_w, fwd, k, c_w
      real theta_hat_i
      real ui_rx, ui_ry, mag_ui
!        
! Wave friction factor for wave and crest half cycle - VA2013 Eq 21
! 
      fw_i=fwi_calc(T_iu, T_i, ahat, ksw)
!
! Wave current friction factor (Madsen and Grant) - VA2013 Eq 18
! Different for crest and trough 
!
#  ifdef BEDLOAD_VANDERA_ZEROCURR
      fwd_i=fw_i
#  else
      fwd_i=alpha*fd+(1.-alpha)*fw_i
#  endif
!
! Magnitude of Shields parameter - VA2013 Eq 17
! 
      theta_hat_i=0.5*fwd_i*uhat_i**2*osmgd
!
! Sheet flow thickness - VA2013 Appendix C C.1 
!
      dsf=dsf_calc(d50, theta_hat_i) ! this dsf is in m 
!
! Compute velocity magnitude based on representative velocities
!  - VA2013 Eq 12 
!
! Get the representative trough half cycle water particle velocity
! as well as full cycle orbital velocity and excursion
!
#  ifdef BEDLOAD_VANDERA_ZEROCURR
      ui_rx=ui_r*alphaw
      ui_ry=0.
#  else
      ui_rx=udelta*COS(phi_curwave)*alphac+ui_r*alphaw
      ui_ry=udelta*SIN(phi_curwave)*alphac
#  endif
!
! mag_ui is set to a min value to avoid non-zero division
!
#  ifdef BEDLOAD_VANDERA_ZEROCURR
      mag_ui=MAX( SQRT(ui_rx*ui_rx), eps )
#  else
      mag_ui=MAX( SQRT(ui_rx*ui_rx+ui_ry*ui_ry), eps )
#  endif
!
! Magnitude of Shields parameter - VA2013 Eq 17
! 
      mag_theta_i=MAX(0.5*fwd_i*osmgd*(mag_ui**2.), eps)
!
! Shields parameter in crest cycle
! rho0 is required for non-dimensionalizing 
!
      theta_ix=ABS(mag_theta_i)*ui_rx/(mag_ui)+tau_wRe*osmgd/rho0
#  ifdef BEDLOAD_VANDERA_ZEROCURR
      theta_iy=0.
#  else
      theta_iy=ABS(mag_theta_i)*ui_ry/(mag_ui)
#  endif
!
! mag_theta_i is set to a min value to avoid non-zero division
!
#  ifdef BEDLOAD_VANDERA_ZEROCURR
      mag_theta_i=MAX( sqrt(theta_ix*theta_ix+theta_iy*theta_iy),eps )
#  else
      mag_theta_i=MAX( sqrt(theta_ix*theta_ix),eps )
#  endif
!
      return
      end
!
!======================================================================
!
      subroutine current_timeperiod(unet, phi_curwave, umax, umin,
     &                              RR, T_c, T_t, Td)
!
! This subroutine returns the following:
! T_c, T_t  : Time period in crest and trough cycle
!
! Modify the crest and trough time periods based on current velocity
! This function was developed by Chris Sherwood and Tarandeep Kalra
!
! The basis for these formulations are formed from Appendix A.3 
! in SANTOSS report (Report nb: SANTOSS_UT_IR3 Date: January 2010)
!
      implicit none
#  include "param.h"
#  include "scalars.h"
      real unet, phi_curwave
      real umax, umin
      real RR, Td
      real T_c, T_t
!
!  Local variables
!
      real unet_xdir, udelta, delt
!
      unet_xdir=unet*cos(phi_curwave)
      if(RR.eq.0.5) then
        T_c=0.5*Td
        T_t=0.5*Td
        if(unet_xdir.ge.umax) then
          T_c=Td
          T_t=0.
        elseif(unet_xdir.le.umin) then
          T_c=0.
          T_t=Td
        elseif(unet_xdir.lt.0.0.and.unet_xdir.gt.umin) then
          delt=ASIN(-unet/umin)/pi
          T_t=T_t*(1.-2.*delt)
          T_c=Td-T_t
        elseif(unet_xdir.gt.0.0.and.unet_xdir.lt.umax) then
          delt=ASIN(unet_xdir/(-umax))/pi
          T_c=T_c*(1.-2.*delt)
          T_t=Td-T_c
        elseif(unet_xdir.eq.0.) then
          T_c=T_c
          T_t=T_t
        endif
      elseif(RR.gt.0.5) then
        T_c=T_c
        T_t=T_t
        if(unet_xdir.ge.umax) then
          T_c=Td
          T_t=0.
        elseif(unet_xdir<=umin) then
          T_c=0.
          T_t=Td
        elseif(unet_xdir.lt.0.0.and.unet_xdir.gt.umin) then
          delt=ASIN(-unet_xdir/umin)/pi
          T_t=T_t*(1.-2.*delt)
          T_c=Td-T_t
        elseif(unet_xdir.gt.0.0.and.unet_xdir.lt.umax) then
          delt=ASIN(unet_xdir/(-umax))/pi
          T_c=T_c*(1.-2.*delt)
          T_t=Td-T_c
        elseif(unet_xdir.eq.0.) then
          T_c=T_c
          T_t=T_t
        endif
      endif
      T_c=MAX(T_c,0.)
      T_t=MAX(T_t,0.)
!
      return
      end
!
!======================================================================
!
      subroutine surface_wave_mod(Td, depth, uhat, T_c, T_cu, T_t, T_tu)
! 
! This subroutine returns the following:
! T_c, T_cu, T_t, T_tu  : Change in time period in crest and 
!                         trough cycle due to particle displacement
!                         under surface waves. 
!
! Crest period extension for horizontal particle displacement.
! Tuning factor eps_eff = 0.55 from measurements GWK Schretlen 2010         
! Equations in Appendix B of SANTOSS Report (SANTOSS_UT_IR3, 2010) 
!
      implicit none
#  include "param.h"
#  include "scalars.h"
      real Td, depth, uhat
      real T_c, T_cu, T_t, T_tu
!
!  Local variables
!
      real eps
      parameter(eps = 1.0E-14)
      real k_wn, eps_eff, c
      real delta_Tc, delta_Tt
      real T_c_new, T_cu_new
      real T_t_new, T_tu_new
      real kh_calc
!
      k_wn=kh_calc(Td,depth)/depth
      c=2.*pi/(k_wn*Td)
!
      eps_eff=0.55
!
      delta_Tc=eps_eff*uhat/(c*pi-eps_eff*2.*uhat)
      T_c_new=T_c+delta_Tc
!
! Avoid non zero values for time periods 
!
      T_c_new=MAX( T_c_new, 0. )
      T_cu_new=MAX( T_cu*T_c_new/T_c, 0. )
!
      delta_Tt=eps_eff*uhat/(c*pi+eps_eff*2.*uhat)
      T_t_new=T_t-delta_Tt
      T_t_new=MAX( T_t_new, 0. )
      T_tu_new=MAX( T_tu*T_t_new/T_t, 0. )
!
      T_c=T_c_new
      T_cu=T_cu_new
      T_t=T_t_new
      T_tu=T_tu_new
!
      return
      end
!
!======================================================================
!
      subroutine ripple_dim(psi, d50, eta, rlen)
! 
! This subroutine returns the following:
! eta, rlen : Ripple dimensions: (height and length) 
!
! Calculate ripple dimensions of O'Donoghue et al. 2006
! based on VA2013 Appendix B
!        
      implicit none 
#  include "param.h"
#  include "scalars.h"
      real psi, d50
      real eta, rlen
      real d50_mm 
      real m_eta, m_lambda, n_eta, n_lambda 
      real eps
      parameter(eps=1.0E-14)
!     
      d50_mm=0.001*d50
      if(d50_mm.lt.0.22) then
        m_eta=0.55
        m_lambda=0.73
      elseif(d50_mm.ge.0.22.and.d50_mm.lt.0.30) then
        m_eta=0.55+(0.45*(d50_mm-0.22)/(0.30-0.22))
        m_lambda=0.73+(0.27*(d50_mm-0.22)/(0.30-0.22))
      else
        m_eta=1.
        m_lambda=1.
      endif
! 
! Smooth transition between ripple regime and bed sheet flow regime 
!
      if(psi.le.190.) then
        n_eta=1.
      elseif(psi.gt.190..and.psi.lt.240.) then
        n_eta=0.5*(1.+cos(pi*(psi-190.)/(50.)))
      elseif(psi.ge.240.) then
        n_eta=0.
      endif
      n_lambda=n_eta
!
      eta=MAX(0.,m_eta*n_eta*(0.275-0.022*psi**0.42))
      rlen=MAX(eps,m_lambda*n_lambda*(1.97-0.44*psi**0.21))
!
      return
      end
!
!======================================================================
!
      subroutine skewness_params( H_s, T, depth, r, phi, Ur )
!        
! Ruessink et al. (2012) provides equations for calculating skewness 
! parameters. Uses Malarkey and Davies equations to get "bb" and "r",
! given input of H_s, T and depth 
!
! r     - skewness/asymmetry parameter r=2b/(1+b^2)      [value]
! phi   - skewness/asymmetry parameter                   [value]
! Su    - umax/(umax-umin)                               [value]
! Au    - amax/(amax-amin)                               [value]
! alpha - tmax/pi                                        [value]
!
      implicit none
#  include "param.h"
#  include "scalars.h"
      real H_s, T, depth
      real Ur
      real r, phi
      real kh_calc
!
! Local variables 
! 
      real p1,p2,p3,p4,p5,p6
      parameter(p1=0.,p2=0.857,p3=-0.471,p4=0.297,p5=0.815,p6=0.672)
      real B, psi, bb
      real k_wn, cff
      real Su, Au
!
! Ruessink et al. (2012), Coastal Engineering 65:56-63.
!
! k is the local wave number computed with linear wave theory.
!
      k_wn=kh_calc(T,depth)/depth    
!
! Ursell number defined as 3/(32*pi^2) * Hs*L^2/h^3 
!               (Fourier coeff for NL terms)
!
      Ur=0.375*H_s/(k_wn*k_wn*depth**3.)
!
! Ruessink et al. (2012), Eq 9
!
      cff=EXP((p3-log10(Ur))/p4)
      B=p1+((p2-p1)/(1.+cff))
!
      psi=-90.*deg2rad*(1.-TANH(p5/Ur**p6))
!
#  ifdef BEDLOAD_VANDERA_ASYM_LIMITS
      B=MIN(B,0.8554) ! according to fig.2, max values w.r.t Ur=24
      psi=MAX(psi,-1.4233) 
#  endif     
! 
! Markaley and Davies, equation provides bb which is "b" in paper
! Check from where CRS found these equations
! 
      bb=sqrt(2.)*B/(sqrt(2.*B**2.+9.))
      r=2.*bb/(bb**2.+1.)
!
! Ruessink et al., 2012 under Equation 12.
!
      phi=-psi-0.5*pi
!
! Skewness Su and asymmetry Au 
!          
      Su=B*cos(psi)
      Au=B*sin(psi)
!
      return
      end
!
!======================================================================
!
      subroutine abreu_points( r, phi, Uw, T, T_c, T_t,
     &                            T_cu, T_tu, umax, umin, RR, beta )
! 
!  Calculate umax, umin, and phases of asymmetrical wave orbital velocity 
!
!  Use the asymmetry parameters from Ruessink et al, 2012
!  to get the umax, umin and phases of asymettrical wave 
!  orbital velocity to be used by Van Der A. 
!  T_c is duration of crest
!  T_cu Duration of accerating flow within crest half cycle
!
      implicit none
#  include "param.h"
#  include "scalars.h"
      real r, phi, Uw, T
      real T_c, T_t, T_cu, T_tu
      real umax, umin, RR, beta
!
! Local variables 
! 
      real b, c, ratio, tmt, tmc, tzd, tzu
      real omega, w, phi_new
      real P, F0, betar_0
!     real T_tu, T_cu, T_c, T_t
      real cff1, cff2, cff
      real Sk, Ak
!
      omega=2.*pi/T
!
      phi_new=-phi

! Malarkey and Davies (Under equation 16b) 
      P=SQRT(1.-r*r)
!
! Malarkey and Davies (Under equation 16b) 
!
      b=r/(1.+P)
!
! Appendix E of Malarkey and Davies 
!
      c=b*SIN(phi_new)
!
      cff1=4.*c*(b*b-c*c)+(1.-b*b)*(1.+b*b-2.*c*c)
      cff2=(1.+b*b)**2.-4.*c*c
      ratio=cff1/cff2
!
! These "if" conditionals prevent ASIN to be between [-1,1] and prevent NaNs
! Not a problem in the MATLAB code
!
      if(ratio.gt.1.)then
        ratio=1.
      endif
      if(ratio.lt.-1.)then
        ratio=-1.
      endif
      tmc=ASIN(ratio)
!
      cff1=4.*c*(b*b-c*c)-(1.-b*b)*(1.+b*b-2.*c*c)
      cff2=(1.+b*b)**2.-4.*c*c
      ratio=cff1/cff2
      if(ratio.gt.1.)then
        ratio=1.
      endif
      if(ratio.lt.-1.)then
        ratio=-1.
      endif
      tmt=ASIN(ratio)
!       
      if(tmc.lt.0.) then
        tmc=tmc+2.*pi
      endif
      if(tmt.lt.0.) then
        tmt=tmt+2.*pi
      endif
! 
! Non dimensional umax and umin, under E5 in Malarkey and Davies 
! 
      umax=1.+c
      umin=umax-2.
!
!     Dimensionalize
!
      umax=umax*Uw
      umin=umin*Uw
!
! phase of zero upcrossing and downcrossing (radians)
!
      tzu=ASIN(b*SIN(phi_new))
      tzd=2.*ACOS(c)+tzu
! 
! MD, equation 17
!
      RR=0.5*(1.+b*SIN(phi_new))
! 
! MD, under equation 18
! 
      if(r.le.0.5) then
        F0=1.0-0.27*(2.*r)**(2.1)
      else
        F0=0.59+0.14*(2.*r)**(-6.2)
      endif
!
! MD, Equation 15a,b 
!
      if(r.ge.0..and.r.lt.0.5)then
        betar_0=0.5*(1.+r)
      elseif(r.gt.0.5.and.r.lt.1.)then
        cff1=4.*r*(1.+r)
        cff2=cff1+1.
        betar_0=cff1/cff2
      endif
!
! MD, Equation 18
!
      cff=SIN((0.5*pi-ABS(phi_new))*F0)/SIN(0.5*pi*F0)
      beta=0.5+(betar_0-0.5)*cff
!
! MD, Table 1, get asymmetry parameterization
! using GSSO (10a,b)
!
      cff=SQRT(2.*(1.+b*b)**3.)
      Sk=3.*SIN(phi_new)/cff
      Ak=-3.*COS(phi_new)/cff
!
! These are the dimensional fractions of wave periods needed by Van der A eqn.
!
      w=1./omega
      T_c=(tzd-tzu)*w
      T_t=T-T_c
      T_cu=(tmc-tzu)*w
      T_tu=(tmt-tzd)*w
!
      return
      end

! 
!======================================================================
!
      FUNCTION kh_calc(Td,depth)
!
! Calculate wave number from Wave period and depth 
!
!#  define KH_HUNT
!
      implicit none
#  include "param.h"
#  include "scalars.h"
      real kh_calc
      real Td, depth
      real omega, Khd
      real cff, x, y, t
#  ifdef KH_HUNT
      real K1, K2, K3, K4, K5, K6
      parameter (K1=0.6666666666, K2=0.3555555555,
     &           K3=0.1608465608, K4=0.0632098765,
     &           K5=0.0217540484, K6=0.0065407983)
#  endif
!
      omega=2.*pi/Td
!
#  ifdef KH_HUNT
!
! Hunt (1979) "Direct solutions of wave dispersion equation"
!
      Khd=depth*omega*omega/g
      kh_calc= SQRT( Khd*Khd+
     &         Khd/(1.+Khd*(K1+Khd*(K2+Khd*(K3+Khd*(K4+
     &         Khd*(K5+K6*Khd)))))) )
#  else
!
! Soulsby (2006) "Simplified calculation of wave orbital velocities"
!
      if(depth.lt.0.) then  ! avoid negative wavenumber
        x=0.
      else 
        x=omega**2.*depth/g
      endif 
!
      if(x.lt.1.) then
        y=SQRT(x)
      else
        y=x
      endif
!      
      t=TANH(y)   ! Iteratively solving 3 times (Soulsby 1997)
      cff=(y*t-x)/(t+y*(1.-t*t))
      y=y-cff
!
      t=TANH(y)
      cff=(y*t-x)/(t+y*(1.-t*t))
      y=y-cff
!
      t=TANH(y)
      cff=(y*t-x)/(t+y*(1.-t*t))
      y=y-cff
!
      kh_calc=y
#  endif /* KH_HUNT */
!
      return
      end
! 
!======================================================================
!
      FUNCTION w_s_calc(d50, rhos)
!
! Critical Shields parameter from Soulsby (1997).
! Dynamics of Marine Sands 
!
      implicit none
#  include "param.h"
#  include "scalars.h"
      real w_s_calc
      real nu
      parameter(nu=1.36E-6)
      real d50, rhos
      real s, dstar
      real cff, cff1
!
      s=rhos/rho0
      dstar=(g*(s-1)/(nu*nu))**(1./3.)*d50
      cff=nu/d50
      cff1=10.36
      w_s_calc=cff*(sqrt(cff1*cff1+1.049*dstar**3.)-cff1)
!
      return
      end
! 
!======================================================================
!
      FUNCTION what_calc(Hs, Td, depth, RR, zws) 
! 
! Second order Stokes theory to get vertical velocity of water particle
! at a given elevation based on santoss_core.m
! 
      implicit none 
      real what_calc
      real eps_inv
      parameter (eps_inv=1.0E14)
      real pi
      parameter (pi=3.14159265358979323846)
      real Hs, Td, depth, RR, zws
      real cff, worb1, worb2, worb 
!      
      worb1=pi*Hs*zws/(Td*depth)
      worb2=worb1*2.*(RR+RR-1.)
!
!  Using the SANTOSS model formulation 
! 
      cff=0.125
      worb=cff*worb1*SQRT(64.-(-worb1+
     &          SQRT(worb1**2.+32.*
     &               worb2**2.))**2./(worb2**2.))+
     &               worb2*SIN(2.*ACOS(cff*(-worb1+
     &          SQRT(worb1**2.+32.*worb2**2.))/worb2))
!
! Prevent worb from going to Infinity when worb2=0.0 
! 
      worb=MIN(worb, eps_inv)
      what_calc=worb
!
      return  
      end
! 
!======================================================================
!
      FUNCTION mu_calc(d50)
!
! Bed roughness factor based on grain size - VA2013 Appendix A
! required for current and wave related bed roughness. 
!
      implicit none
      real mu_calc
      real d50, d50_mm
!
      d50_mm=d50*1000.
!
      if(d50_mm.le.0.15) then
        mu_calc=6.
      elseif(d50_mm.gt.0.15.and.d50_mm.lt.0.2) then
        mu_calc=6.-5.*((d50_mm-0.15)/(0.2-0.15))
      elseif(d50_mm.ge.0.2) then
        mu_calc=1.
      endif
!
      return
      end
! 
!======================================================================
!
      FUNCTION ksd_calc(d50, d90, mu, theta_avg, eta, rlen)
!
! Current-related bed roughness - VA2013 Appendix A.1
!
      implicit none
      real ksd_calc
      real d50, d90, mu, theta_avg, eta, rlen
      real ripple_fac
!
      rlen=MAX(rlen,d50)
      ripple_fac=0.4*eta**2./rlen
      ksd_calc=MAX( 3.*d90, d50*(mu+6.*(theta_avg-1.)) )
     &         + ripple_fac
!
      return
      end
! 
!======================================================================
!
      FUNCTION ksw_calc(d50, mu, theta_avg, eta, rlen)
!
! Wave related bed roughness - VA2013 Eq A.5
! 
      implicit none
      real ksw_calc
      real d50, mu, theta_avg, eta, rlen
      real ripple_fac, ksw
!
      rlen=MAX(rlen,d50)
      ripple_fac=0.4*eta**2./rlen
      ksw_calc=MAX( d50, d50*(mu+6.*(theta_avg-1.)) )
     &         + ripple_fac
!
      return
      end
! 
!======================================================================
!
      FUNCTION fw_calc(ahat, ksw)
!
! Full-cycle wave friction factor - VA2013 Eq A.4
!
      implicit none
      real fw_calc
      real ahat, ksw, ratio, fw
!
      ratio=ahat/ksw
      if(ratio.gt.1.587) then
        fw_calc=0.00251*EXP(5.21*(ratio)**(-0.19))
      else
        fw_calc=0.3
      endif
!
      return
      end
! 
!======================================================================
!
      FUNCTION fd_calc_santoss(udelta, delta, ksd)
!
! Current friction factor Assuming logarithmic velocity profile.
! - VA2013 Eq 20
!
      implicit none
      real fd_calc_santoss
      real udelta, delta, ksd
      real min_udelta
      parameter(min_udelta=1.0E-4)
      real von_k
      parameter(von_k=0.41)
 
      if(udelta.lt.min_udelta) then
        fd_calc_santoss=0.
      else
        fd_calc_santoss=2.*(von_k/LOG(30.*delta/ksd))**2.
      endif
     
      return
      end
! 
!======================================================================
!
      FUNCTION fd_calc_madsen(udelta, mag_bstrc)
!
! Current friction factor from current stresses. 
! 
      implicit none
#  include "param.h"
#  include "scalars.h"
!
      real fd_calc_madsen
      real eps
      parameter(eps=1.E-14)
      real min_udelta
      parameter(min_udelta=1.0E-4)
      real udelta, mag_bstrc
! 
      if(udelta.lt.min_udelta) then
        fd_calc_madsen=0.
      else
        fd_calc_madsen=MAX((mag_bstrc/(0.5*udelta*udelta)),eps)
!        fd_calc_madsen=MIN(fd_calc_madsen,2.) 
      endif
!     
      return
      end
! 
!======================================================================
!
      FUNCTION fwi_calc(T_iu, T_i, ahat, ksw)
!
! Wave friction factor for half cycles - VA2013 Eq 21
!       
      implicit none
#  include "param.h"
#  include "scalars.h"
      real fwi_calc
      real eps
      parameter(eps = 1.0E-14)
      real T_iu, T_i, ahat, ksw
      real c1, ratio, fwi
      real cff1, cff2, cff3
!
      fwi_calc=0.3
!      
      c1=2.6
      ratio=ahat/ksw
      if(ratio.gt.1.587) then
        cff1=MAX( (T_iu/T_i),0. )
        cff2=(2.*cff1)**c1
        cff3=cff2*ratio
!
! These "if" condition prevents arithematic overflow error
! when 0.0**-0.19
!
        if(cff3.le.0.) then 
          fwi_calc=0.
        else 
          fwi_calc=0.00251*EXP(5.21*(cff3)**(-0.19))
        end if
      end if 
!
      return
      end
! 
!======================================================================
!
      FUNCTION dsf_calc(d50, theta_i)
!
! Sheet flow thickness - VA2013 Appendix C.1.
!       
      implicit none
      real dsf_calc
      real d50, theta_i
      real d50_mm
      real cff
!
      d50_mm=d50*1000.
      if(d50_mm.le.0.15)then
        cff=25.*theta_i
      elseif(d50_mm.gt.0.15.and.d50_mm.lt.0.2)then
        cff=25.-(12.*(d50_mm-0.15)/0.05)
      elseif(d50_mm.ge.0.2)then
        cff=13.*theta_i
      endif
      dsf_calc=MAX(d50*cff,d50)
!
      return
      end
! 
!======================================================================
!
      FUNCTION theta_cr_calc(d50, rhos)
!
! Critical Shields parameter from Soulsby (1997).
!
      implicit none
#  include "param.h"
#  include "scalars.h"
      real theta_cr_calc
      real nu
      parameter(nu=1.36E-6)
      real d50, rhos
      real s, dstar
      real cff1, cff2
!
      s=rhos/rho0
      dstar=(g*(s-1)/(nu*nu))**(1./3.)*d50
      cff1=0.3/(1.+1.2*dstar)
      cff2=0.055*(1.-EXP(-0.02*dstar))
      theta_cr_calc=cff1+cff2
!
      return
      end
# endif /* BEDLOAD_VANDERA */
!======================================================================
!
# ifdef SED_FLOCS
      subroutine sed_flocmod_tile (Istr,Iend,Jstr,Jend)
        USE flocmod, ONLY : flocmod_main, flocmod_comp_g
      implicit none 
#  include "param.h"
#  include "scalars.h"
#  include "ocean3d.h"
#  include "mixing.h"
#  include "bbl.h"
#  include "sediment.h"
#  include "forces.h"

      integer Istr,Iend,Jstr,Jend
      integer i, j, k
!
!  Variable declarations for floc model

      real Gval, diss, nueau
      real epsilon8, ustr2, effecz

      if (l_testcase) then
        CALL flocmod_comp_g(Gval, time)
      else
        epsilon8=epsilon(1.0)
      endif

      do j=Jstr,Jend
        do i=Istr,Iend
          do k=1,N       
            if (.not. l_testcase) then
               
# if defined FLOC_TURB_DISS && defined GLS_MIXING
!
! Dissipation from turbulence clossure
!
              if (k.eq.1) then
                diss = Eps_gls(i,j,k)
              elseif (k.eq.N) then
                diss = Eps_gls(i,j,k-1)
              else 
                diss = 0.5*(Eps_gls(i,j,k-1)+Eps_gls(i,j,k))
              endif
#  ifdef FLOC_BBL_DISS
!
! Dissipation from wavecurrent bottom stress
! NOT READY FOR PRIME TIME: NEEDS VERTICAL DISTRIBUTION
! As first cut, use turbulence closure
!
              if (k.eq.1) then       
#   ifdef BBL
                ustr2 =sqrt((bustrw(i,j)**2                
     &                      +bvstrw(i,j)**2))
                effecz = (ustr2**0.5)*Pwave(i,j)*0.5/pi
                diss = MAX((ustr2**(1.5))/(vonKar*effecz),diss)
#   else
                ustr2 =sqrt((bustr(i,j)**2.0+bvstr(i,j)**2.0 ))
                diss = MAX((ustr2**(1.5))/(vonKar*                  
     &                       (z_r(i,j,1)-z_w(i,j,0))),diss)
#   endif
              endif
#  endif
# else 
!
! Constant tke dissipation
!
                diss = epsilon8         
# endif              
                nueau = 1.5E-6
                Gval=sqrt(diss/nueau)
            endif ! not l_testcase
            
            CALL flocmod_main( dt, 
     &            t(i,j,k,nnew,itrc_mud:itrc_mud+NMUD-1),  
     &            Gval )

          enddo
        enddo
      enddo  !J_LOOP

      return
      end subroutine sed_flocmod_tile
      
# endif /* SED_FLOCS */
!
!====================================================================
!
!  UP5 interpolation for bedload fluxes
!
!=====================================================================
!
# ifdef BEDLOAD_UP5
      function flux5(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua)
        implicit none 
        REAL    :: q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2
        REAL    :: ua
        REAL    :: flux5, flux6

        flux6 = ( 37.*(q_i+q_im1) - 8.*(q_ip1+q_im2)
     &       +(q_ip2+q_im3) )/60.0

        flux5= flux6-sign(1.,ua)*(
     &      (q_ip2-q_im3)-5.*(q_ip1-q_im2)+10.*(q_i-q_im1) )/60.0
        return
      end function flux5
# endif
#else       
      subroutine sediment_empty
      end
#endif



