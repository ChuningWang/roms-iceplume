
#include "cppdefs.h"

#if defined SOLVE3D && defined ONLINE_ANALYSIS

!------------------------------------------------------------------------------
! PROCEDURE
!------------------------------------------------------------------------------
!
!> @note
!
! DESCRIPTION: 
!
!> @brief Fonction donnant la valeur de la variable choisie
!!     en un point donne.
!
!> @details Attention lors de l'ajout d'une variable penser completer
!!     var_grid_oa pour specifier le type de point sur lequel
!!     se trouve la nouvelle variable.
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Initial version
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - newvar_oa : var_oa pointer version in case of long list
!!    of if-statement to determine which var-config the function
!!    var_oa should return (cf, repated call to construct the correlation
!!    product between psi_oa and the requested model fields which 
!!    are applied inside the time, period and space loops in main_oa).
!! - Francis Auclair, B. Lemieux-Dudon, C. Nguyen
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!> @date 2015 January
!> @todo BLXD Croco-OA interface 2020
!! WARNING only var-config type 11 and 99 have been tested
!------------------------------------------------------------------------------

      real function var_oa(   tile                                    &
           ,ivar_v                                                    &
           ,cnb_v                                                     &
           ,i_v                                                       &
           ,j_v                                                       &
           ,k_v                                                       &
           ,lv_v                                                      &
           ,ls1_v                                              ) 

!BLXD_TILE_ISSUE
!      use module_oa_variables, only : vardp_test_oa
      use module_tile_oa, only : st, tile_space_str, sts, tile_test_str 
! BLXD 2020 removing useless modules
!     use module_oa_time
!     use module_oa_space
!     use module_oa_periode
!     use module_oa_stock
      use module_oa_level, only : v2lev_oa
! BLXD 2020 module_oa_upd does not exist anymore
!     use module_oa_upd
!     use module_nrj
      use scalars

      implicit none
# include "ocean2d.h"
# include "ocean3d.h"
# include "grid.h"
#ifdef NBQ
# include "nbq.h"
#endif

      integer, intent(in) ::                                          &
            cnb_v                                                     &
           ,tile                                                      &
           ,ivar_v                                                    &
           ,i_v                                                       &
           ,j_v                                                       &
           ,k_v                                                       &
           ,lv_v                                                      &
           ,ls1_v

      ! BLXD_TILE_ISSUE 
      ! pst declaration moved from module_tile_oa to local subroutien
      ! => openMP PRIVATE variable with respect to Croco tile-thread loop
      type(tile_space_str), pointer :: pst => null()

      type(tile_test_str), pointer  :: psts => null()

!*********************************************************************
! Calling code "simple" variables
!*********************************************************************

      if (ivar_v.eq.1) then
         var_oa = u(i_v,j_v,k_v,nstp)-ubar(i_v,j_v,fast_indx_out)
      elseif (ivar_v.eq.2) then
         var_oa = v(i_v,j_v,k_v,nstp)-vbar(i_v,j_v,fast_indx_out)
#  ifdef NBQ
      elseif (ivar_v.eq.3) then
         var_oa = wz(i_v,j_v,k_v,nstp) 
#  endif
      elseif (ivar_v.eq.4) then
         var_oa = u(i_v,j_v,k_v,nstp)
      elseif (ivar_v.eq.5) then
         var_oa = v(i_v,j_v,k_v,nstp)
      elseif (ivar_v.eq.6) then
         var_oa = ubar(i_v,j_v,fast_indx_out)
      elseif (ivar_v.eq.7) then
         var_oa = vbar(i_v,j_v,fast_indx_out)
      elseif (ivar_v.eq.8) then
         var_oa = t(i_v,j_v,k_v,nstp,itemp)
# ifdef SALINITY
      elseif (ivar_v.eq.9) then
         var_oa = t(i_v,j_v,k_v,nstp,isalt)
# endif
      elseif (ivar_v.eq.11) then
         var_oa = rho(i_v,j_v,k_v)
      elseif (ivar_v.eq.20) then
         pst => st(tile)
         var_oa = pst%wlev_oa(v2lev_oa(lv_v))%z(ls1_v)
         pst => null()
      elseif (ivar_v.eq.56) then
         var_oa = ubar(i_v,j_v,fast_indx_out)
      elseif (ivar_v.eq.61) then
         var_oa = rho(i_v,j_v,k_v)
!*********************************************************************
! OA test variable
!*********************************************************************
      elseif (ivar_v.eq.97) then
         psts => sts(tile)
         ! var_oa = 0.5*vardp_test_oa(i_v,j_v,1) test var2d code 97 + var3d code 99
         var_oa = psts%vardp_test(i_v,j_v,1)
         psts => null()
      elseif (ivar_v.eq.98) then
         psts => sts(tile)
         ! var_oa = 0.5*vardp_test_oa(i_v,j_v,1) test var2d code 98 + var3d code 99
         var_oa = psts%vardp_test(i_v,j_v,k_v)
         psts => null()
      elseif (ivar_v.eq.99) then
         psts => sts(tile)
         var_oa = psts%vardp_test(i_v,j_v,k_v)
         psts => null()
!*********************************************************************
! Calling code "mixed" variables (composite)
!*********************************************************************
      elseif (ivar_v.eq.100) then
!........sum of the squared vel_u and vel_v at t-grid : 
!        "kinetic energy without rho factor" just to provide and example to handle "mixed" variables
!        var_oa = ( ( u(i_v,j_v,k_v) + u(i_v+1,j_v,k_v) )**2 + ( v(i_v,j_v,k_v) + v(i_v,j_v+1,k_v) )**2 ) / 8.0
! BLXD add a comment to notify termination
         stop
      endif


      return
      end function var_oa

# else
      real function var_oa(                                           &
            ivar_v                                                    &
           ,cnb_v                                                     &
           ,i_v                                                       &
           ,j_v                                                       &
           ,k_v                                                       &
           ,lv_v                                                      &
           ,ls1_v                                              ) 

      integer, intent(in) ::                                          &
            cnb_v                                                     &
           ,ivar_v                                                    &
           ,i_v                                                       &
           ,j_v                                                       &
           ,k_v                                                       &
           ,lv_v                                                      &
           ,ls1_v
      var_oa = var_oa
      return
      end function var_oa



#endif



