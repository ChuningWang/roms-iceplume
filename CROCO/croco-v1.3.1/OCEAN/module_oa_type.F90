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
!! - Francis Auclair, B. Lemieux-Dudon, C. Nguyen
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!
!> @date 2015
!> @todo BLXD
!  organize consistent calc. with double precision complex (ifort, gfortran, -r8, croco preproc.)
!------------------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     equipe d'oceanographie cotiere du laboratoire d'aerologie
!     laboratoire d'aerologie
!     cnrs - universite paul sabatier - observatoire midi pyrenees
!     14 avenue edouard belin - 31400 toulouse - france
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!************************************************************************
!.....module module_oa_type: declaration de type
!************************************************************************

      module module_oa_type

!.....declaration des types:

!----> type ondelettes:

      type type_wf

      integer config,variable,t_indice
      complex(8),dimension(:),pointer::coef 

      end type type_wf

!----> type level:

      type type_level

      real   ,dimension(:),pointer::z
      real   ,dimension(:),pointer::rhp
      integer,dimension(:),pointer::k

      end type type_level

      end module module_oa_type

#else
      module module_oa_type_empty
      end module
#endif
