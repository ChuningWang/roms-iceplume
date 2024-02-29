#include "cppdefs.h"

#if defined MUSTANG
/*
   !&E==========================================================================
   !&E                   ***  coupler_define_MUSTANG  ***
   !&E
   !&E
   !&E ** Purpose : definitions of dimensions, variables and parameters 
   !&E               for MUSTANG         
   !&E 
   !&E ** Description : must be completed by the user
   !&E          when coupling with a hydrodynamic model
   !&E
   !&E==========================================================================
*/
# define sed_MUSTANG_HOST sed_MUSTANG_CROCO

/* CROCO */
# define iscreenlog stdout
# define ierrorlog stdout
# define iwarnlog stdout
# define NAME_SUBS vname(1,indxT+ntrc_salt+isubs)

/* Directory where are namelists files */
# define REPFICNAMELIST 'MUSTANG_NAMELIST'

/* Spatial Grid limits definition  of loops inside the domain - except meshes at open boundaries */
# define IMIN_GRID 1
# define IMAX_GRID Lm
# define JMIN_GRID 1
# define JMAX_GRID Mm

/* dimensions table definition */
# define PROC_IN_ARRAY       GLOBAL_2D_ARRAY   
# define PROC_IN_ARRAY_m1p2  GLOBAL_2D_ARRAY
# define PROC_IN_ARRAY_m1p1  GLOBAL_2D_ARRAY
# define PROC_IN_ARRAY_m2p2  GLOBAL_2D_ARRAY
# define PROC_IN_ARRAY_0p1   GLOBAL_2D_ARRAY

/* dimensions of variables in hydro modele */
# define ARRAY_EROS_FLUX_s2w PROC_IN_ARRAY,1:NT
# define ARRAY_SETTL_FLUX_w2s PROC_IN_ARRAY,1:NT
# define ARRAY_SETTL_FLUXSUM_w2s PROC_IN_ARRAY,1:NT
# define ARRAY_WATER_FLUX_INPUTS PROC_IN_ARRAY,1:NB_LAYER_WAT
# define ARRAY_BATHY_H0 GLOBAL_2D_ARRAY
# define ARRAY_WATER_ELEVATION GLOBAL_2D_ARRAY,1:4
# define ARRAY_VELOCITY_U GLOBAL_2D_ARRAY,1:4
# define ARRAY_VELOCITY_V GLOBAL_2D_ARRAY,1:4
# define ARRAY_CELL_SURF GLOBAL_2D_ARRAY
# define ARRAY_CELL_DX GLOBAL_2D_ARRAY
# define ARRAY_CELL_DY GLOBAL_2D_ARRAY
/* dimensions of variables in hydro modele !*/
# define ARRAY_WATER_CONC GLOBAL_2D_ARRAY,N,3,NT
# define ARRAY_DHSED GLOBAL_2D_ARRAY 
# define ARRAY_FROFON GLOBAL_2D_ARRAY 
# define ARRAY_Z0HYDRO PROC_IN_ARRAY 

#if defined key_MUSTANG_V2 && defined key_MUSTANG_debug
# define ARRAY_LATLON GLOBAL_2D_ARRAY
#endif

#define  ARRAY_morpho GLOBAL_2D_ARRAY

/* general variable hydro , bathy, time ... defined in hydro model but using by MUSTANG
!*/
# define NUMBER_PI pi
# define NB_LAYER_WAT N
# define COORD_SIGMA sc_r
# define BATHY_H0 h
# ifdef WET_DRY
# define RESIDUAL_THICKNESS_WAT D_wetdry
# else
# define RESIDUAL_THICKNESS_WAT 0.
# endif
# define WATER_ELEVATION zeta
# define CELL_DX om_r
# define CELL_DY on_r
# define CELL_SURF surf_cell
# define BAROTROP_VELOCITY_U ubar
# define BAROTROP_VELOCITY_V vbar
# define TIME_STEP dt   /* in MARS :  time step declared in rlg and therefore also in MUSTANG */
# define TRANSPORT_TIME_STEP dt /* in MARS : solving equations every half time step (in rlg)*/
# define CURRENT_TIME time
# define RHOREF rho0
# define TEMPREF_LIN 10.0 
# define SALREF_LIN 35.0
# define GRAVITY g
# define BOTTOM_THICK_LAYER epn_bottom
# define Z0HYDRO zob
# define WATER_CONCENTRATION t  /* water concentration in hydro model (=cv_wat in MARS)*/
# define DHSED dh  
#if defined key_MUSTANG_V2 && defined key_MUSTANG_debug
# define LATITUDE latr
# define LONGITUDE lonr
#endif

/* surface elevation (i,j) and current could have different dimensions*/
# define SURF_ELEVATION_ij WATER_ELEVATION(i,j,3) 
# define CURRENTU_ij BAROTROP_VELOCITY_U(i+1,j,3) 
# define CURRENTV_ij BAROTROP_VELOCITY_V(i,j+1,3) 

# define CURRENTV_ip1jm1 BAROTROP_VELOCITY_V(i+1,j,3) 
# define CURRENTV_ip1j   BAROTROP_VELOCITY_V(i+1,j+1,3) 
# define CURRENTV_im1jm1 BAROTROP_VELOCITY_V(i-1,j,3)
# define CURRENTV_im1j   BAROTROP_VELOCITY_V(i-1,j+1,3)

# define CURRENTU_im1jp1 BAROTROP_VELOCITY_U(i,j+1,3) 
# define CURRENTU_ijp1   BAROTROP_VELOCITY_U(i+1,j+1,3) 
# define CURRENTU_im1jm1 BAROTROP_VELOCITY_U(i,j-1,3) 
# define CURRENTU_ijm1   BAROTROP_VELOCITY_U(i+1,j-1,3) 

/* name of fluxes exchange between MUSTANG and hydro model */
# define EROS_FLUX_s2w flx_s2w_CROCO /* Erosion flux from sediment to water */
# define SETTL_FLUX_w2s flx_w2s_CROCO  /* Tendance Flux de depot eau vers sediment (=flx_s2w in MARS) only for particulate */
# define SETTL_FLUXSUM_w2s flx_w2s_sum_CROCO /* effective deposit Flux (sum) from water to sediment in hydro model (=flx_w2s_sum in MARS)*/
/* to locate the number of variables simulated by MUSTANG in the host hydro model (used in coupleur to tranfer exchange arrays  */      
# define IV_HOSTMODEL itsubs1+iv-1
# define ITEMP_HOSTMODEL itemp
# define ISAL_HOSTMODEL itemp+1
# define WATER_FLUX_INPUT_BOTCELL phieau_CROCO(:,:,1) /* Flux d eau apporte dans la maille de fond in hydro model (=phieau in MARS) */

/* Lateral Erosion in neighboring cells (could depend on grid architecture) */
# define HTOT_NEAR_E htot(i+1,j)
# define HTOT_NEAR_W htot(i-1,j)
# define HTOT_NEAR_N htot(i,j+1)
# define HTOT_NEAR_S htot(i,j-1)
# define SURF_NEAR_E CELL_SURF(i+1,j)
# define SURF_NEAR_W CELL_SURF(i-1,j)
# define SURF_NEAR_N CELL_SURF(i,j+1)
# define SURF_NEAR_S CELL_SURF(i,j-1)
# define V_NEAR_E (CURRENTV_ip1jm1 + CURRENTV_ip1j)
# define V_NEAR_W (CURRENTV_im1jm1 + CURRENTV_im1j)
# define U_NEAR_N (CURRENTU_im1jp1 + CURRENTU_ijp1)
# define U_NEAR_S (CURRENTU_im1jm1 + CURRENTU_ijm1)

/* sliding proces of fluid mud slope of the bottom in the neighboring cells 
(could depend on grid architecture) */
# define SLOPE_W ((BATHY_H0(i-1,j)-BATHY_H0(i,j))/CELL_DX(i,j))
# define SLOPE_E ((BATHY_H0(i+1,j)-BATHY_H0(i,j))/CELL_DX(i,j))
# define SLOPE_N ((BATHY_H0(i,j+1)-BATHY_H0(i,j))/CELL_DY(i,j))
# define SLOPE_S ((BATHY_H0(i,j-1)-BATHY_H0(i,j))/CELL_DY(i,j))

#endif
