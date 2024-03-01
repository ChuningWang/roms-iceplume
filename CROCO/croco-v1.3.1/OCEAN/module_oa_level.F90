#include "cppdefs.h"
#ifdef ONLINE_ANALYSIS
!------------------------------------------------------------------------------
!                               NHOMS
!                Non Hydrostatic Ocean Modeling System      
!------------------------------------------------------------------------------
!
!> @note <a href="http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/index_snh_home.htm"> Main web documentation </a>
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
!!  - Initial version
!! - B. Lemieux-Dudon : 
!!  - modification in the tracking of isopycnal levels_:
!!      - ifl_l flag eliminated, replaced by counting analysis of type 20 in nzlevel_oa, 
!!        and testing if nzlevel_oa>0.
!!      - enables to diminish the size of the structured type array wlev, now sized 
!!        to nzlevel_oa intead of nzv_oa (the total number of OA analysis requested 
!!        in the simulation).
!! - Francis Auclair, B. Lemieux-Dudon, C. Nguyen
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!> @date 2015 January
!> @todo
!
!------------------------------------------------------------------------------

!************************************************************************
!.....module module_oa_level
!************************************************************************

      module module_oa_level

      use module_oa_type

      integer::nzlevel_oa                      

      !BLXD_TILE_ISSUE
      !type(type_level),dimension(:),allocatable::wlev_oa

      integer,dimension(:),allocatable::lev2v_oa
      integer,dimension(:),allocatable::v2lev_oa

      end module module_oa_level
#else
      module module_oa_level_empty
      end module
#endif
