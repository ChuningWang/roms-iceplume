! $Id: bbl.h 1458 2014-02-03 15:01:25Z gcambon $
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
#if defined BBL || defined SEDIMENT

/*
** Include file "bbl.h"
**********************************************************************
** Copyright (c) 2003 Rutgers/UCLA                                  **
************************************************* Hernan G. Arango ***
****************************************** Christopher R. Sherwood ***
**************************************************** Meinte Blaas  ***
**                                                                  **
** Abed         wind-induced, bed wave excursion amplitude (m).     **
** Hripple      Bed ripple height (m).                              **
** Lripple      Bed ripple length (m).                              **
** w_set        Input settling velo (m/s) sediment only bbl(rho pts)**
** Sdens        Input sediment grain density (kg/m3) "   "   "   "  **
** Ssize        Input sediment grain diameter (m)    "   "   "   "  **
** taucb        Input threshold stress bedload(N/m^2)"   "   "   "  **
** Ubed         Wind-induced, bed wave orbital U-velocity (m/s).    **
** Vbed         Wind-induced, bed wave orbital V-velocity (m/s).    **
** Zbnot        Physical hydraulic bottom roughness  (m)            **
** Zbapp        Total apparent hydraulic bottom roughness (m).      **
** bustrw       Kinematic bottom skin stress (m2/s2) in the         **
**                XI-direction at horizontal Rho-points.            **
** bvstrw       Kinematic bottom skin stress (m2/s2) in the         **
**                ETA-direction at horizontal Rho-points.           **
**********************************************************************
*/

      real Abed(GLOBAL_2D_ARRAY)
      common /bbl_Abed/ Abed      

      real Hripple(GLOBAL_2D_ARRAY)
      common /bbl_Hripple/ Hripple

      real Lripple(GLOBAL_2D_ARRAY)
      common /bbl_Lripple/ Lripple

      real w_set(GLOBAL_2D_ARRAY)
      common /bbl_wset/ w_set

      real Sdens(GLOBAL_2D_ARRAY)
      common /bbl_Sdens/ Sdens

      real Ssize(GLOBAL_2D_ARRAY)
      common /bbl_Ssize/ Ssize

      real taucb(GLOBAL_2D_ARRAY)
      common /bbl_taucb/ taucb
      
      real Zbnot(GLOBAL_2D_ARRAY)
      common /bbl_Zbnot/ Zbnot

      real Zbapp(GLOBAL_2D_ARRAY)
      common /bbl_Zbapp/ Zbapp

      real bustrw(GLOBAL_2D_ARRAY)
      common /bbl_bustrw/ bustrw

      real bvstrw(GLOBAL_2D_ARRAY)
      common /bbl_bvstrw/ bvstrw
#endif
