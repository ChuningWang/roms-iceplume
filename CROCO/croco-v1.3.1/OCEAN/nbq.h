! $Id:$
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
#ifdef M3FAST
!**********************************************************************
      logical M2bc_nbq_flag
      common /nbq_M2bc/ M2bc_nbq_flag

!********************************************************************** 
      integer iteration_nbq_max  
      common /nbq_var1/ iteration_nbq_max             
      integer iteration_nbq      
      common /nbq_var2/ iteration_nbq      
      integer ifl_nbq     
      common /nbq_var3/ ifl_nbq     
      integer slip_nbq   
      common /nbq_var4/ slip_nbq    

!**********************************************************************
      real soundspeed_nbq(GLOBAL_2D_ARRAY) 
      common /nbq_param1/ soundspeed_nbq                              
      real soundspeed2_nbq(GLOBAL_2D_ARRAY) 
      common /nbq_param2/ soundspeed2_nbq

      double precision time_nbq               
      common /nbq_param3/ time_nbq
      double precision csvisc1_nbq                
      common /nbq_param4/ csvisc1_nbq 
      double precision csvisc2_nbq            
      common /nbq_param5/ csvisc2_nbq 
      double precision cw_int_nbq       
      common /nbq_param6/ cw_int_nbq
      double precision ifl_imp_nbq
      common /nbq_param7/ ifl_imp_nbq

!**********************************************************************
      integer ndtnbq
      common /time_nbq1/ ndtnbq
      real dtnbq
      common /time_nbq2/ dtnbq 
      real csound_nbq
      common /nbq_csound/ csound_nbq
      real visc2_nbq
      common /nbq_visc2/ visc2_nbq

      real dtgrid_nbq
      common /nbq_dtgrid/ dtgrid_nbq

!**********************************************************************
      real qdmu_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_qdmu_nbq/ qdmu_nbq
      real qdmv_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_qdmv_nbq/ qdmv_nbq
# ifdef NBQ
      real qdmw_nbq(GLOBAL_2D_ARRAY,0:N)
      common /nbq_qdmw_nbq/ qdmw_nbq
# endif

!**********************************************************************       
# ifdef NBQ
      real thetadiv_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_thetadiv_nbq/ thetadiv_nbq
      real thetadiv2_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_thetadiv2_nbq/ thetadiv2_nbq
      real thetadiv3_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_thetadiv3_nbq/ thetadiv3_nbq
# endif

!**********************************************************************
      real ru_int_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_ruint/ ru_int_nbq
      real rv_int_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_rvint/ rv_int_nbq

      real ru_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_ru/ ru_nbq
      real rv_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_rv/ rv_nbq

      real ru_nbq_avg2(GLOBAL_2D_ARRAY,N)
      common /avg2_runbq/ ru_nbq_avg2
      real rv_nbq_avg2(GLOBAL_2D_ARRAY,N)
      common /avg2_rvnbq/ rv_nbq_avg2

# ifdef NBQ
      real rw_int_nbq(GLOBAL_2D_ARRAY,0:N)
      common /nbq_rwint/ rw_int_nbq   
      real rw_nbq(GLOBAL_2D_ARRAY,0:N)
      common /nbq_rw/ rw_nbq
      real rw_nbq_avg2(GLOBAL_2D_ARRAY,0:N)
      common /avg2_rwnbq/ rw_nbq_avg2
      real rho_nbq(GLOBAL_2D_ARRAY,N)
      common/nbq_rho_nbq/rho_nbq
# endif

!**********************************************************************
      integer inc_faststep
      common/nbq_inc_faststep/inc_faststep
      integer nb_faststep
      common/nbq_nb_faststep/nb_faststep

      real DU_nbq(GLOBAL_2D_ARRAY)
      common /nbq_DU_nbq/ DU_nbq
      real DV_nbq(GLOBAL_2D_ARRAY)
      common /nbq_DV_nbq/ DV_nbq

      real ru_int_nbq_2d (GLOBAL_2D_ARRAY)  
      common /nbq_ruint_2d/ru_int_nbq_2d
      real rv_int_nbq_2d (GLOBAL_2D_ARRAY)  
      common /nbq_rvint_2d/rv_int_nbq_2d

# ifdef NBQ
      real rho_grd(GLOBAL_2D_ARRAY,N)
      common/nbq_rho_grd/rho_grd
      real rho_bak(GLOBAL_2D_ARRAY,N)
      common/nbq_rho_bak/rho_bak
#  ifdef NBQ_MASS
      real rho_nbq_avg1(GLOBAL_2D_ARRAY,0:N)
      common /avg1_rhonbq/ rho_nbq_avg1
      real rhobar_nbq(GLOBAL_2D_ARRAY,4)
      common /nbq_rhobar/ rhobar_nbq
      real rhobar_nbq_avg1(GLOBAL_2D_ARRAY)
      common /nbq_rhobar_AVG1/ rhobar_nbq_avg1
#  endif
# else
      real rubar_nbq(GLOBAL_2D_ARRAY)
      common /nbq_rubar/ rubar_nbq
      real rvbar_nbq(GLOBAL_2D_ARRAY)
      common /nbq_rvbar/ rvbar_nbq

      real rubar_sum(GLOBAL_2D_ARRAY)
      common /nbq_rubar_sum/ rubar_sum
      real rvbar_sum(GLOBAL_2D_ARRAY)
      common /nbq_rvbar_sum/ rvbar_sum
# endif

!**********************************************************************
      real Hzw_half_nbq(GLOBAL_2D_ARRAY,0:N)
      common /grid_Hzw_half_nbq/ Hzw_half_nbq

# ifdef NBQ
      real zw_nbq(GLOBAL_2D_ARRAY,0:N,4)
      common /nbq_zw/ zw_nbq
# endif

# ifdef NBQ_HZCORRECT
       real Hz_correct(GLOBAL_2D_ARRAY,N)
       common /grid_Hz_correct/ Hz_correct
#  ifdef NBQ_HZCORR_DEBUG
      real  Hz_corr(GLOBAL_2D_ARRAY,N) 
      common/corr_Hz/Hz_corr
#  endif
# endif

# ifdef NBQ_HZ_PROGNOSTIC
      real Hz_bak2(GLOBAL_2D_ARRAY,1:N)
      common /nbq_H_bak2/ Hz_bak2
# endif

!**********************************************************************
# ifdef NBQ
#  ifdef NBQ_GRID_SLOW
      real dthetadiv_nbqdz(GLOBAL_2D_ARRAY,0:N,2)
      common /nbq_nods3/ dthetadiv_nbqdz
      real dZdxq_w(GLOBAL_2D_ARRAY,0:N+1)
      common /nbq_nods5/ dZdxq_w
      real dZdyq_w(GLOBAL_2D_ARRAY,0:N+1)
      common /nbq_nods7/ dZdyq_w
#  else
      real dthetadiv_nbqdz(GLOBAL_2D_ARRAY)
      common /nbq_nods3/ dthetadiv_nbqdz
      real dZdxq_w(GLOBAL_2D_ARRAY,0:N+1)
      common /nbq_nods5/ dZdxq_w
      real dZdyq_w(GLOBAL_2D_ARRAY,0:N+1)
      common /nbq_nods7/ dZdyq_w
#  endif /* NBQ_GRID_SLOW */
# endif

!**********************************************************************
# ifdef NBQ
      real wsurf_nbq(GLOBAL_2D_ARRAY)
      common /nbq_wsurf/ wsurf_nbq
      real usurf_nbq(GLOBAL_2D_ARRAY)
      common /nbq_usurf/ usurf_nbq
      real vsurf_nbq(GLOBAL_2D_ARRAY)
      common /nbq_vsurf/ vsurf_nbq
# endif

!**********************************************************************
# if defined OBC_NBQ && defined OBC_NBQORLANSKI
#  ifdef OBC_COM_WEST
      real qdmu_nbq_west(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_unbq_west/ qdmu_nbq_west
      real qdmv_nbq_west(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_vnbq_west/ qdmv_nbq_west
#   ifdef NBQ
      real qdmw_nbq_west(GLOBAL_1D_ARRAYETA,0:N,2)
      common /bry_wnbq_west/ qdmw_nbq_west
      real  rho_nbq_west(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_rnbq_west/ rho_nbq_west
#   endif
#  endif
#  ifdef OBC_COM_EAST
      real qdmu_nbq_east(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_unbq_east/ qdmu_nbq_east
      real qdmv_nbq_east(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_vnbq_east/ qdmv_nbq_east
#   ifdef NBQ
      real qdmw_nbq_east(GLOBAL_1D_ARRAYETA,0:N,2)
      common /bry_wnbq_east/ qdmw_nbq_east
      real  rho_nbq_east(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_rnbq_east/ rho_nbq_east
#   endif
#  endif
#  ifdef OBC_COM_SOUTH
      real qdmu_nbq_south(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_unbq_south/ qdmu_nbq_south
      real qdmv_nbq_south(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_vnbq_south/ qdmv_nbq_south
#   ifdef NBQ
      real qdmw_nbq_south(GLOBAL_1D_ARRAYXI,0:N,2)
      common /bry_wnbq_south/ qdmw_nbq_south
      real  rho_nbq_south(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_rnbq_south/ rho_nbq_south
#   endif
#  endif
#  ifdef OBC_COM_NORTH
      real qdmu_nbq_north(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_unbq_north/ qdmu_nbq_north
      real qdmv_nbq_north(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_vnbq_north/ qdmv_nbq_north
#   ifdef NBQ
      real qdmw_nbq_north(GLOBAL_1D_ARRAYXI,0:N,2)
      common /bry_wnbq_north/ qdmw_nbq_north
      real  rho_nbq_north(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_rnbq_north/ rho_nbq_north
#   endif
#  endif
# endif     

!**********************************************************************
# ifdef NBQ_NUDGING
      real NBQnudgcof(GLOBAL_2D_ARRAY)
      common /nbq_nudg/ NBQnudgcof
# endif

!**********************************************************************
# ifdef ACOUSTIC
      real  period_exp  
      common/ACOUS1/period_exp
      real  for_a_exp   
      common/ACOUS2/for_a_exp
      real  dg_exp     
      common/ACOUS3/dg_exp 
      real  hmax_exp    
      common/ACOUS4/hmax_exp
      real  amp_exp
      common/ACOUS4/amp_exp
# endif

#endif /* M3FAST */

  
