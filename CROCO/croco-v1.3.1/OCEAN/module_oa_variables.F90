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
!! - B. Lemieux-Dudon
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!!  - BLXD nzvc_oa variables moved to the module_oa_interface with
!!    parametrized dimension maxtyp_oa.
!!  - Max # of convolution windows simultaneously openned pre-calculated 
!!    nmsimult_oa
!> @date 2015
!> @todo
!
!------------------------------------------------------------------------------

!************************************************************************
!.....module module_oa_variables
!************************************************************************

      module module_oa_variables


!.....variables:

      integer::                                               &
           nzv_oa                                                ! nombre de variables (assoc. aux configurations)
           !,nzvs_oa                                           & ! nombre total de points (taille de la structure spatiale 2d du vecteur d etat)
           !,nzvs3d_oa                                           ! nombre total de points (taille de la structure spatiale 3d du vecteur d etat)

      integer,dimension(:),allocatable::                      & ! (nmv_oa)
           swt_d_oa                                           & ! caracteristiques spatiales de la variable (voir notebook_wavelet)
           ,swt_t_oa                                          & ! caracteristiques temporelles de la variable (voir notebook_wavelet)          
           ,tv_oa                                             & ! variable symphonie associee         
           ,cnb_oa                                            & ! call number: position of the call in the baroclinic / barotropic time step.
           ,updv_oa
!          ,tvc_oa                                            & ! configuration de la variable associee

!.....flag pour traces et deboggage
      integer,dimension(:),allocatable::                      & ! (nmv_oa)
            ltrec_fst_oa                                      & !  records 1st temporal OA in simulation 
           ,ltrec_lst_oa                                      & ! ...
           ,if_first_rec_oa                                   &
           ,if_last_rec_oa                                     

!.....flag pour variable de type scalogram
      logical,dimension(:),allocatable:: tv_sclg_oa

!.....Variable for NHOMS energy analysis if you wish to call OA at specific point in the code
!     for specific variable that are updated at some irregular code place or time step (e.g., time splitting)
!     depends on iv_m so it should have nzv_oa size has cnb_oa
!     BLXD bug correction declared with wrong size, ie
!     integer, dimension(1:6000) :: des_oa  

      integer,dimension(:),allocatable::                      &
          des_oa                                                ! (nmv_oa)

      integer,dimension(:),allocatable::                      &  ! (nmv_oa)
           tpsi_oa                                               ! type d atome utilise


      integer,dimension(:),allocatable::                      & ! (nmv_oa+1)
          begvt_oa                                              ! structure temporelle du vecteur d etat


      integer                                                         &
           flag_nrj_oa                                                  ! flag pour le calcul de l'energie
 
      integer,dimension(:),allocatable::                              & !(nmv_oa)
           save_oa                                                    & ! sauvegarde finale de la variable dans un fichier
          ,tupd_oa                                                      ! pour variables communes

      logical :: test_analysis
      logical :: isopycne_analysis

      !> Variables required to read OA namelists
      character(len=50), dimension(:), allocatable :: nzc_oa_names
      integer, dimension(:), allocatable           :: nzc_oa_vartyp
      logical, dimension(:), allocatable           :: nzc_oa_scalogram

      end module module_oa_variables

#else
      module module_oa_variables_empty
      end module
#endif
