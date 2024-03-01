#include "cppdefs.h"
#ifdef ONLINE_ANALYSIS 
!------------------------------------------------------------------------------
! Croco interface for Online Analysis (OnlineAnalysis = OA)
!------------------------------------------------------------------------------
!> @brief Croco interface for Online Analysis (OA=OnlineAnalysis)
!! - Performs online spectral and wavelet analysis.
!! - Preliminary Croco version : spring 2020.
!> @authors 
!! - Initial version : Francis Auclair , Jochem Floor and Ivane Pairaud:
!! - Stand-alone version + Namelists + optimization : B. Lemieux-Dudon
!! - Preliminary Croco version : F. Auclair, B. Lemieux-D., C. Nguyen
!> @todo BLXD Optimize the Ocean model / OnlineAnalysis interface for Croco.
!> @todo BLXD hdle intel/gnu dble prec. real, complex (-r8, croco mpc.F preproc)
!------------------------------------------------------------------------------


      subroutine croco_oa (icall)

!#ifdef NBQ
!      use module_nh 
!      use module_nbq
!#endif
      use module_interface_oa, only : init_oa, main_oa
      use scalars
      implicit none
! BLXD
! # include "param.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "grid.h"
# include "nbq.h"

#include "def_bounds.h"

      integer :: icall
      integer :: i,j,k    

      character (len=200) :: dum1_c,dum2_c

      integer ::                       &
        maskr_c (GLOBAL_2D_ARRAY,1:N)  &
       ,masku_c (GLOBAL_2D_ARRAY,1:N)  &
       ,maskv_c (GLOBAL_2D_ARRAY,1:N)  &
       ,maskf_c (GLOBAL_2D_ARRAY,1:N)   

      real ::                          &
        hu_c(GLOBAL_2D_ARRAY)          &
       ,hv_c(GLOBAL_2D_ARRAY)          

#ifndef SPHERICAL
      real ::                          &
        xu_c (GLOBAL_2D_ARRAY)         &
       ,yu_c (GLOBAL_2D_ARRAY)         &
       ,xv_c (GLOBAL_2D_ARRAY)         &
       ,yv_c (GLOBAL_2D_ARRAY)     
#endif

      !real :: u_test (GLOBAL_2D_ARRAY,1:N) 

      integer  ::                                                     &
         istr_oa                                                      &
        ,jstr_oa                                                      &
        ,iend_oa                                                      &
        ,jend_oa                                                       

      integer, parameter  :: verbose_oa=5

       istr_oa  = 1
       iend_oa  = Lm 
       jstr_oa  = 1
       jend_oa  = Mm 

      if (icall==0) then
!**********************************************************************
! Initializations
!**********************************************************************

       dum1_c="INPUT"
       dum2_c="OUTPUT"

#ifndef MASKING
       maskr_c = 1
       masku_c = 1
       maskv_c = 1
       maskf_c = 1
#else
       do i=istr_oa,iend_oa
       do j=jstr_oa,jend_oa
          maskr_c(i,j,1:N)  = INT(rmask(i,j))
          maskf_c(i,j,1:N)  = INT(pmask(i,j))
       enddo
       enddo

       do i=istr_oa,iend_oa
       do j=jstr_oa,jend_oa
          masku_c(i,j,1:N)=INT(umask(i,j))
       enddo
       enddo

       do i=istr_oa,iend_oa
       do j=jstr_oa,jend_oa
          maskv_c(i,j,1:N)=INT(vmask(i,j))
       enddo
       enddo
#endif

! #BLXD - -> +
       do i=istr_oa,iend_oa
       do j=jstr_oa,jend_oa
          hu_c(i,j)=0.5*(h(i,j)+h(i-1,j))
       enddo
       enddo

       do i=istr_oa,iend_oa
       do j=jstr_oa,jend_oa
          hv_c(i,j)=0.5*(h(i,j)+h(i,j-1))
       enddo
       enddo

#ifndef SPHERICAL
       ! #BLD 2020 changing minus to plus sign (check)
       do i=istr_oa,iend_oa
       do j=jstr_oa,jend_oa
          xu_c(i,j)=0.5*(xr(i,j)+xr(i-1,j))
          yu_c(i,j)=0.5*(yr(i,j)+yr(i-1,j))
       enddo
       enddo

       do i=istr_oa,iend_oa
       do j=jstr_oa,jend_oa
          xv_c(i,j)=0.5*(xr(i,j)+xr(i,j-1))
          yv_c(i,j)=0.5*(yr(i,j)+yr(i,j-1))
       enddo
       enddo
#endif

    !    pointer_approach=.false.  & 

      ! BLXD 2020 kount0 set to nstart-1
      ! check to clean
      ! write(*,*) ' IMPT nstart , icc = ', ntstart, iic

       call init_oa(                                   & 

       directory_in_oa =dum1_c                         &
      ,directory_out_oa=dum2_c                         &
      ,io_unit_oa=stdout                               &
#ifdef MPI
      ,if_print_node_oa=(mynode==0)                            &
#else
      ,if_print_node_oa=.true.                            &
#endif                                             
#ifdef MPI
      ,mynode_oa=mynode                                &
#else
      ,mynode_oa=0                                      &
#endif
      ,iic_oa=iic                                      &
      ,kount0=ntstart-1                                  &
      ,nt_max=ntimes                                   &   
      ,dti=dt                                          &

      ,imin=istr_oa, imax=iend_oa                      &
      ,jmin=jstr_oa, jmax=jend_oa                      &
      ,kmin=1,        kmax=N                           & 

#ifdef SPHERICAL
      ,lon_t=lonr(istr_oa:iend_oa,jstr_oa:jend_oa) & 
#else
      ,lon_t=xr(istr_oa:iend_oa,jstr_oa:jend_oa)   & 
#endif
      ,lon_t_lbound=(/istr_oa,jstr_oa/)              &
      ,lon_t_ubound=(/iend_oa,jend_oa/)              &
#ifdef SPHERICAL
      ,lat_t=latr(istr_oa:iend_oa,jstr_oa:jend_oa) & 
#else
      ,lat_t=yr(istr_oa:iend_oa,jstr_oa:jend_oa)   & 
#endif
      ,lat_t_lbound=(/istr_oa,jstr_oa/)              &
      ,lat_t_ubound=(/iend_oa,jend_oa/)              &
#ifdef SPHERICAL
      ,lon_u=lonu(istr_oa:iend_oa,jstr_oa:jend_oa) & 
#else
      ,lon_u=xu_c(istr_oa:iend_oa,jstr_oa:jend_oa) & 
#endif
      ,lon_u_lbound=(/istr_oa,jstr_oa/)              &
      ,lon_u_ubound=(/iend_oa,jend_oa/)              &
#ifdef SPHERICAL
      ,lat_u=latu(istr_oa:iend_oa,jstr_oa:jend_oa) & 
#else
      ,lat_u=yu_c(istr_oa:iend_oa,jstr_oa:jend_oa) & 
#endif
      ,lat_u_lbound=(/istr_oa,jstr_oa/)              &
      ,lat_u_ubound=(/iend_oa,jend_oa/)              &
#ifdef SPHERICAL
      ,lon_v=lonv(istr_oa:iend_oa,jstr_oa:jend_oa) & 
#else
      ,lon_v=xv_c(istr_oa:iend_oa,jstr_oa:jend_oa) & 
#endif
      ,lon_v_lbound=(/istr_oa,jstr_oa/)              &
      ,lon_v_ubound=(/iend_oa,jend_oa/)              &
#ifdef SPHERICAL
      ,lat_v=latv(istr_oa:iend_oa,jstr_oa:jend_oa) &
#else
      ,lat_v=yv_c(istr_oa:iend_oa,jstr_oa:jend_oa) &
#endif
      ,lat_v_lbound=(/istr_oa,jstr_oa/)              &
      ,lat_v_ubound=(/iend_oa,jend_oa/)              &
#ifdef SPHERICAL
      ,lon_f=lonr(istr_oa:iend_oa,jstr_oa:jend_oa) &   ! Wrong grid !
#else
      ,lon_f=xr(istr_oa:iend_oa,jstr_oa:jend_oa)   &   ! Wrong grid !
#endif
      ,lon_f_lbound=(/istr_oa,jstr_oa/)              &
      ,lon_f_ubound=(/iend_oa,jend_oa/)              &
#ifdef SPHERICAL
      ,lat_f=latr(istr_oa:iend_oa,jstr_oa:jend_oa) &   ! Wrong grid !
#else
      ,lat_f=yr(istr_oa:iend_oa,jstr_oa:jend_oa)   &   ! Wrong grid !
#endif
      ,lat_f_lbound=(/istr_oa,jstr_oa/)              &
      ,lat_f_ubound=(/iend_oa,jend_oa/)              &

      ,mask_t=maskr_c(istr_oa:iend_oa,jstr_oa:jend_oa,1:N) &  
      ,mask_t_lbound=(/istr_oa,jstr_oa,1/)                   &
      ,mask_t_ubound=(/iend_oa,jend_oa,N/)                   &
      ,mask_f=maskf_c(istr_oa:iend_oa,jstr_oa:jend_oa,1:N) & ! Wrong grid !
      ,mask_f_lbound=(/istr_oa,jstr_oa,1/)                   &
      ,mask_f_ubound=(/iend_oa,jend_oa,N/)                   &
      ,mask_u=masku_c(istr_oa:iend_oa,jstr_oa:jend_oa,1:N) & 
      ,mask_u_lbound=(/istr_oa,jstr_oa,1/)                   &
      ,mask_u_ubound=(/iend_oa,jend_oa,N/)                   &
      ,mask_v=maskv_c(istr_oa:iend_oa,jstr_oa:jend_oa,1:N) &
      ,mask_v_lbound=(/istr_oa,jstr_oa,1/)                   &
      ,mask_v_ubound=(/iend_oa,jend_oa,N/)                   &

      ,h_w=h(istr_oa:iend_oa,jstr_oa:jend_oa)              & 
      ,h_w_lbound=(/istr_oa,jstr_oa/)                        &
      ,h_w_ubound=(/iend_oa,jend_oa/)                        &
      ,h_u=hu_c(istr_oa:iend_oa,jstr_oa:jend_oa)           & 
      ,h_u_lbound=(/istr_oa,jstr_oa/)                        &
      ,h_u_ubound=(/iend_oa,jend_oa/)                        &
      ,h_v=hv_c(istr_oa:jstr_oa,jstr_oa:jend_oa)           & 
      ,h_v_lbound=(/istr_oa,jstr_oa/)                        &
      ,h_v_ubound=(/iend_oa,jend_oa/)                        &
      ,h_f=h(istr_oa:iend_oa,jstr_oa:jend_oa)              & ! Wrong grid !
      ,h_f_lbound=(/istr_oa,jstr_oa/)                        &
      ,h_f_ubound=(/iend_oa,jend_oa/)                        &

      ,rhp_t=rho(istr_oa:iend_oa,jstr_oa:jend_oa,1:N)      &
      ,rhp_t_lbound=(/istr_oa,jstr_oa,1/)                    &
      ,rhp_t_ubound=(/iend_oa,jend_oa,N/)                    &

      ,depth_t=z_r(istr_oa:iend_oa,jstr_oa:jend_oa,1:N)    & 
      ,depth_t_lbound=(/istr_oa,jstr_oa,1/)                  &
      ,depth_t_ubound=(/iend_oa,jend_oa,N/)              )

!#endif      

!     pointer_approach=.true.  

      ! BLXD check to clean
      ! write(*,*) ' IMPT nstart , icc = ', ntstart, iic

      call main_oa(                                                  &
                            ichoix=0                                 &
                           ,ivar_m=1                                 &
      ,io_unit_oa=stdout                                             &
#ifdef MPI                                                           
      ,if_print_node_oa=(mynode==0)                                  &
#else                                                                
      ,if_print_node_oa=.true.                                       &
#endif                                                               
#ifdef MPI                                                           
      ,mynode_oa=mynode                                              &
#else                                                                
      ,mynode_oa=0                                                   &
#endif
                           ,iic_oa=iic                               &
                           ,dti=dt                                   &
                           ,nt_max=ntimes                            &
                           ,imin=istr_oa, imax=iend_oa               &
                           ,jmin=jstr_oa, jmax=jend_oa               &
                           ,kmin=1,        kmax=N                    &

                           ,rhp_t=rho(istr_oa:iend_oa,jstr_oa:jend_oa,1:N)      &
                           ,rhp_t_lbound=(/istr_oa,jstr_oa,1/)                    &
                           ,rhp_t_ubound=(/iend_oa,jend_oa,N/)                    &

                           ,depth_t=z_r(istr_oa:iend_oa,jstr_oa:jend_oa,1:N)    & 
                           ,depth_t_lbound=(/istr_oa,jstr_oa,1/)                  &
                           ,depth_t_ubound=(/iend_oa,jend_oa,N/)      )


      elseif (icall==1) then
!**********************************************************************
! At every (internal-mode) time step
!**********************************************************************

! BLXD   Spring 2020
!        Currently, the croco_oa code (OA=OnlineAnalysis) has been set as a mix of :
!        - the Stand-alone OA code (argument-based inteface between Ocean Model and OA),
!        - the original OA code (OA imbedded/intricated winthin the Ocean model).
!        In the original and the stand-alone OnlineAnalysis code, 
!        the model fields to which applying the analysis are 
!        identified by the mean of OA configuration codes (see namelist_oa, nzc_oa_vartyp).
!        The correlation product between psi_oa and the model field is performed
!        thanks to repetitive calls to the external var_oa real scalar function.
!        Based on a series of if-statement which tests the OA configuration code, 
!        the var_oa function takes the value of a specific model field at a given location i,j,k.
!        The var_oa function is called by main_oa to cumulate the correlation product 
!        at each model (or computation) time step. Compared to the original OA code,
!        the stand-alone OnlineAnalysis code can also apply the var_oa external function.
!        In this latter case, there is no need to pass the user-requested model field 
!        (to be analysed) to any of the two main subroutines init_oa and main_oa
!        (since var_oa function plays the role of the interface).
!        However, for the purpose of tracking the isopycnal movements, 
!        the rhp_t argument was introduced in the stand-alone OnlineAnalysis version.
!        It must be stressed out that the rhp_t argument rhp IS NOT applied 
!        for performing any spectral analysis on the 3D model density field. 
!        BUG : as such, args rhp_t cannot be used to pass a user-defined test function
!        since it is ONLY applied to track isopycnal levels. To use a test function 
!        please choose the test variable option (eg, nzc_oa_vartyp=99) and use
!        the test_oa function and the vardp_test_oa variable (see, module_oa_inteface) 
!
!      u_test=23*cos(2.*3.14159/500.*dt*real(iic)) !!-cos(2*3.14159/50.*dt*real(iic))
!     do i=istr_oa,iend_oa
!     do j=jstr_oa,jend_oa
!     do k=1,N
!       !u_test(i,j,k)=u(i,j,k,nstp)-ubar(i,j,fast_indx_out)
!        u_test(i,j,k)=real(i)
!     enddo
!     enddo
!     enddo

      call main_oa(                                                    &
                            ichoix=0                                   &
                           ,ivar_m=1                                   &
      ,io_unit_oa=stdout                                               &
#ifdef MPI                                                             
      ,if_print_node_oa=(mynode==0)                                    &
#else                                                                  
      ,if_print_node_oa=.true.                                         &
#endif                                                                 
#ifdef MPI                                                             
      ,mynode_oa=mynode                                                &
#else                                                                  
      ,mynode_oa=0                                                     &
#endif
                           ,iic_oa=iic                                 &
                           ,dti=dt                                     &
                           ,nt_max=ntimes                              &
                           ,imin=istr_oa, imax=iend_oa                 &
                           ,jmin=jstr_oa, jmax=jend_oa                 &
                           ,kmin=1,        kmax=N                      &

                           ,rhp_t=rho(istr_oa:iend_oa,jstr_oa:jend_oa,1:N)   &
                           ,rhp_t_lbound=(/istr_oa,jstr_oa,1/)               &
                           ,rhp_t_ubound=(/iend_oa,jend_oa,N/)               &

                           ,depth_t=z_r(istr_oa:iend_oa,jstr_oa:jend_oa,1:N) & 
                           ,depth_t_lbound=(/istr_oa,jstr_oa,1/)             &
                           ,depth_t_ubound=(/iend_oa,jend_oa,N/)      )


      endif

      end subroutine croco_oa
#else
      subroutine croco_oa_empty
      end subroutine croco_oa_empty
#endif
