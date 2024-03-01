! $Id: var3d_oa_out.F
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
!======================================================================
!
!> @brief Croco interface for Online Analysis (OA=OnlineAnalysis)
!! - Prepares spectral and wavelet analysis outputs.
!!
!! @details called in output_oa.F 
!!  written to support tile-threads and 
!!  possible dual OpenMP-MPI parallelization. 
!!  See Online_spectral_diags called by main.F and step.F
!!  More history/info in source module_interface_oa.F90
!!
!> @authors  
!! - B. Lemieux-Dudon
!!  - Croco Tile-thread compliant version (2021) 
!!    supporting tile-threads and possible dual OpenMP-MPI 
!!    parallelization of the horizontal domain
!!     => output_oa : OA updates global array looping on tiles 
!!        and finally send them to XIOS (out of the loop)
!!     => var2d_oa_out, var3d_oa_out : 
!!        2D using Croco working arrays work2d, work2d2
!!     => module_parameter_oa : allocating var3d_oa, var2d_oa 
!!        to the croco GLOBAL_3D_ARRAY, 2D_ARRAY 
!!  - based on a preliminary Croco-OA interface version : spring 2020.
!! - More history (authors, comments) in source module_interface_oa.F90
!> @todo BLXD TODO test var3d_oa_out with work.h as applied in
!!  var2d_oa_out
!
!  REFERENCE:
!  
!======================================================================

#include "cppdefs.h"
! BLXD includes :
! - cppdefs_dev.h
! - set_global_definitions.h
!   with GLOBAL_2D_ARRAY
!        PRIVATE_SCRATCH_2D_ARRAY... 
#if defined ONLINE_ANALYSIS && defined SOLVE3D

      subroutine var3d_oa_out (tile, ic, iv, var3d_re, var3d_im
     &                             , vnam_oa_r, vnam_oa_i)

      implicit none
# include "param.h"

! "Global" 3D work array (full domain or MPI subdomain)
!  Argt Out calculated over the tile horizontal dimension
      real, intent(inout) :: var3d_re(GLOBAL_2D_ARRAY,1:N)
      real, intent(inout) :: var3d_im(GLOBAL_2D_ARRAY,1:N)
      integer, intent(in) :: iv, ic
! Cannot work, require that the loop over ic, iv is moved
! outside var3d, var3d_oa_out
      character(len=22), intent(inout) :: vnam_oa_r, vnam_oa_i 

      integer tile
C$    integer trd, omp_get_thread_num
! BLXD private scratch not needed since no A2d, A3d array 
! passed in argument of the _tile subroutine
! # include "private_scratch.h"
! compute_tile_nounds => Istr,Iend,Jstr,Jend assoc. to tile
# include "compute_tile_bounds.h"

C$    trd=omp_get_thread_num()
      call var3d_oa_out_tile (Istr,Iend,Jstr,Jend
     &                        ,var3d_re, var3d_im
     &                        ,vnam_oa_r, vnam_oa_i
     &                        ,ic, iv )


      end subroutine var3d_oa_out

#else /* ONLINE_ANALYSIS && SOLVE3D */

      subroutine var3d_oa_out_empty
      return
      end subroutine var3d_oa_out_empty

#endif /* ONLINE_ANALYSIS && SOLVE3D */

!
!----------------------------------------------------------------------
!
!  Processing the outputs of the online spectral analysis module
!  before outputs
!
!----------------------------------------------------------------------
!
#ifdef ONLINE_ANALYSIS

      subroutine var3d_oa_out_tile ( Istr,Iend,Jstr,Jend
     &                              ,var3d_re,  var3d_im
     &                              ,vnam_oa_r, vnam_oa_i
     &                              ,ic, iv )
!
      use module_interface_oa, only : var3d_oa
     & ,itgt_glob, jtgt_glob, ktgt
     & ,tvar_oa ,tgv3d_oa, tgv_oa, tgvnam_oa, nzvc_oa
     & ,verbose_oa
      use module_oa_variables,  only : tv_oa     
     & ,tupd_oa  ! BLXD added for test
      use module_oa_periode,    only : begc_oa, tvc_oa, tc_oa, nzc_oa
     & ,swt_wfpf_oa
      use scalars
      implicit none
! >> Include with variable declaration
! # include "param.h"
      integer, intent(in) :: Istr, Iend, Jstr, Jend
      ! BLXD each tile will complete the GLOBAL_2D_ARRAY var2d_re, var3d_im
      ! In case of OpenMP threads, it should require and $OMP BARRIER
      ! at the end of the tile loop, shouldn't it ?  
      real, intent(inout) :: var3d_re(GLOBAL_2D_ARRAY,1:N)
     &                      ,var3d_im(GLOBAL_2D_ARRAY,1:N) 
      integer, intent(in) :: iv, ic
      character(len=22), intent(inout)  :: vnam_oa_r, vnam_oa_i 

! Local variables
      integer :: i, j, k, ivc
      integer, parameter :: verbl=5
      !integer imin,imax,jmin,jmax
! BLXD not needed :
# include "grid.h"
!# include "ocean3d.h"
!# include "scalars.h"
! >> Include with computation instructions
# include "compute_auxiliary_bounds.h"
! From param.h, MPI subdom., tile 
! From compute_tile_bounds.h with Istr,Iend,Jstr,Jend
! Computes Physical Bdy indices IstrR, IendR, JstrR, JendR
!          Interior U or V points (start) IstrU, JstrV

! TODO BLXD instruction and now computations ?

#  define I_EXT_RANGE Istr-1,Iend+1
#  define J_EXT_RANGE Jstr-1,Jend+1

# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif

# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif

# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif

      
! 3Drho-point
        ivc = nzvc_oa( tc_oa(ic) )
        var3d_oa_grid_point : if ( tgv_oa( tv_oa(iv) ) == 1 ) then     

        if_oa_re_im_rho_3D : if ( swt_wfpf_oa(iv)==4 ) then
      ! #BLXD 2020 REAL/DIMAG for real/imag part of dble prec var3d_oa 
      ! Warning mpc.F parser REAL -> REAL real_type
      ! Warning : different gnu, intel standards 
            do k=1,N
!                do j=Jstr,Jend
!                    do i=Istr,Iend
                do j=JstrR,JendR
                    do i=IstrR,IendR
                        var3d_re( i, j, k) = 
     & REAL(
     & DBLE( var3d_oa( i,j,k, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ))))
     &              SWITCH rmask(i,j)
                       var3d_im( i, j, k) = 
     & REAL(
     & DIMAG( var3d_oa( i,j,k, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ))))
     &              SWITCH rmask(i,j)

#ifdef OA_TRACES && ( MILES || IGW)
                       if(verbose_oa>verbl) then
                       call traces_oa_ik(i,j,k,itgt_glob,ktgt) 
                       call trace2_oa_ik
     &                    (i,j,k,itgt_glob,ktgt,var3d_re(i,j,k)
     &                                         ,var3d_im(i,j,k) ) 
                       endif
#endif
                    enddo
                enddo
            enddo
            vnam_oa_r = get_vnam( vnamfmt='(a5,a5,i3.3,a1,i3.3,a5)'
     &              ,vnamdp='3d_r_',vnams='_real',verbose=verbose_oa)
            vnam_oa_i = get_vnam( vnamfmt='(a5,a5,i3.3,a1,i3.3,a5)'
     &              ,vnamdp='3d_r_',vnams='_imag',verbose=verbose_oa)
      
! if_oa_re_im_rho_3D
        else 

            do k=1,N
!                do j=Jstr,Jend
!                    do i=Istr,Iend
                do j=JstrR,JendR
                    do i=IstrR,IendR
                        var3d_re( i, j, k) = 
     & var3d_oa( i,j,k, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ))
     &              SWITCH rmask(i,j)

                    enddo
                enddo
            enddo

            vnam_oa_r = get_vnam( vnamfmt='(a5,a5,i3.3,a1,i3.3)'
     &              ,vnamdp='3d_r_',vnams='',verbose=verbose_oa)

        endif if_oa_re_im_rho_3D
     
! 3Du-point
        else if ( tgv_oa( tv_oa(iv) ) == 2 ) then !u-point
       
        if_oa_re_im_u_3D : if ( swt_wfpf_oa(iv)==4 ) then
! #BLXD REAL/DIMAG for real/imag part of dble prec var3d_oa 
! Warning mpc.F parser REAL -> REAL real_type
! Warning : different gnu, intel standards 
      
            do k=1,N
!                do j=Jstr,Jend
!                    do i=Istr,Iend
                do j=JU_RANGE
                    do i=IU_RANGE

                        var3d_re( i, j, k) = 
     & REAL(
     & REAL( var3d_oa( i,j,k, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ))))
     &              SWITCH umask(i,j)
                       var3d_im( i, j, k) = 
     & REAL(
     & DIMAG( var3d_oa( i,j,k, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ))))
     &              SWITCH umask(i,j)

                    enddo
                enddo
            enddo
            vnam_oa_r = get_vnam( vnamfmt='(a5,a5,i3.3,a1,i3.3,a5)'
     &              ,vnamdp='3d_u_',vnams='_real',verbose=verbose_oa)
            vnam_oa_i = get_vnam( vnamfmt='(a5,a5,i3.3,a1,i3.3,a5)'
     &              ,vnamdp='3d_u_',vnams='_imag',verbose=verbose_oa)

! if_oa_re_im_u_3D 
        else

            do k=1,N
!                do j=Jstr,Jend
!                    do i=Istr,Iend
                do j=JU_RANGE
                    do i=IU_RANGE
                        var3d_re( i, j, k) = 
     &  var3d_oa( i,j,k, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ))
     &              SWITCH umask(i,j)
                    enddo
                enddo
            enddo
            vnam_oa_r = get_vnam( vnamfmt='(a5,a5,i3.3,a1,i3.3)'
     &              ,vnamdp='3d_u_',vnams='',verbose=verbose_oa)

        endif if_oa_re_im_u_3D

! 3Dv-point
        else if ( tgv_oa( tv_oa(iv) ) == 3 ) then
       
        if_oa_re_im_v_3D : if ( swt_wfpf_oa(iv)==4 ) then
! #BLXD REAL/DIMAG for real/imag part of dble prec var3d_oa 
! Warning mpc.F parser REAL -> REAL real_type
! Warning : different gnu, intel standards 
            do k=1,N
!                do j=Jstr,Jend
!                    do i=Istr,Iend
                do j=JV_RANGE
                    do i=IV_RANGE
                         
                        var3d_re( i, j, k) =
     & REAL(
     & DBLE( var3d_oa( i,j,k, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ))))
     &              SWITCH vmask(i,j)
                       var3d_im( i, j, k) = 
     & REAL(
     & DIMAG( var3d_oa( i,j,k, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ))))
     &              SWITCH vmask(i,j)

                    enddo
                enddo
            enddo
            vnam_oa_r = get_vnam( vnamfmt='(a5,a5,i3.3,a1,i3.3,a5)'
     &              ,vnamdp='3d_v_',vnams='_real',verbose=verbose_oa)
            vnam_oa_i = get_vnam( vnamfmt='(a5,a5,i3.3,a1,i3.3,a5)'
     &              ,vnamdp='3d_v_',vnams='_imag',verbose=verbose_oa)

! if_oa_re_im_v_3D
        else

            do k=1,N
!                do j=Jstr,Jend
!                    do i=Istr,Iend
                do j=JV_RANGE
                    do i=IV_RANGE
                        var3d_re( i, j, k) = 
     & var3d_oa( i, j, k, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ) )
     &              SWITCH vmask(i,j)
                    enddo
                enddo
            enddo
            vnam_oa_r = get_vnam( vnamfmt='(a5,a5,i3.3,a1,i3.3)'
     &              ,vnamdp='3d_v_',vnams='',verbose=verbose_oa)

        endif if_oa_re_im_v_3D

! 3Df-point
        else if ( tgv_oa( tv_oa(iv) ) == 4 ) then

        if_oa_re_im_f_3D : if ( swt_wfpf_oa(iv)==4 ) then
! #BLXD REAL/DIMAG for real/imag part of dble prec var3d_oa 
! Warning mpc.F parser REAL -> REAL real_type
! Warning : different gnu, intel standards 
            do k=1,N
!                do j=Jstr,Jend
!                    do i=Istr,Iend
                do j=JV_RANGE
                    do i=IU_RANGE
                        var3d_re( i, j, k) =
! BLD WARNING replacing pmask by pmask2
     & REAL(
     & DBLE( var3d_oa( i, j, k, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ) ) ) )
     &              SWITCH pmask2(i,j)
                       var3d_im( i, j, k) = 
     & REAL(
     & DIMAG( var3d_oa( i, j, k, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ) ) ) )
     &              SWITCH pmask2(i,j)

                    enddo
                enddo
            enddo
            vnam_oa_r = get_vnam( vnamfmt='(a5,a5,i3.3,a1,i3.3,a5)'
     &              ,vnamdp='3d_p_',vnams='_real',verbose=verbose_oa)
            vnam_oa_i = get_vnam( vnamfmt='(a5,a5,i3.3,a1,i3.3,a5)'
     &              ,vnamdp='3d_p_',vnams='_imag',verbose=verbose_oa)

! if_oa_re_im_f_3D
        else

            do k=1,N
!                do j=Jstr,Jend
!                    do i=Istr,Iend
                do j=JV_RANGE
                    do i=IU_RANGE
                        var3d_re( i, j, k) = 
     & var3d_oa( i, j, k, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ) )
     &              SWITCH pmask2(i,j)
 
                    enddo
                enddo
            enddo

            vnam_oa_r = get_vnam( vnamfmt='(a5,a5,i3.3,a1,i3.3)'
     &              ,vnamdp='3d_p_',vnams='',verbose=verbose_oa)

        endif if_oa_re_im_f_3D

        endif var3d_oa_grid_point

      contains

      subroutine traces_oa_ik(i,j,k,iiglob_tgt,k_tgt)        
        implicit none
        integer, intent(in)    :: i, j , k
        integer, intent(in)    :: iiglob_tgt, k_tgt
        integer                :: io_node, ii_glob

        ii_glob = i + iminmpi-1 
     
        if ( ( iiglob_tgt == -999 ) .or. ( k_tgt == -999) ) return
        if ( (ii_glob == iiglob_tgt) .and. (k==k_tgt) ) then
            if (iv==1) then
            io_node = 3000+mynode
            write (io_node,*) 'var3d_oa iv1 =',j,var3d_oa(i,j,k,tupd_oa(iv))
            else if (iv==2) then
            io_node = 4000+mynode
            write (io_node,*) 'var3d_oa iv2 =',j,var3d_oa(i,j,k,tupd_oa(iv))
            end if
        end if

        return
      end subroutine traces_oa_ik

      subroutine trace2_oa_ik(i,j,k,iiglob_tgt,k_tgt,tab3d_re,tab3d_im)        
        implicit none
        real, intent(in) :: tab3d_re
        real, intent(in) :: tab3d_im
        integer, intent(in)    :: i, j , k
        integer, intent(in)    :: iiglob_tgt, k_tgt
        integer                :: io_node, ii_glob

        ii_glob = i + iminmpi-1 
      
        if ( ( iiglob_tgt == -999 ) .or. ( k_tgt == -999) ) return
        if ( (ii_glob == iiglob_tgt) .and. (k==k_tgt) ) then
            if (iv==1) then
            io_node = 5000+mynode
            write (io_node,*) 'Re/Im var3d_oa iv1 =',j,tab3d_re
     &                                                ,tab3d_im
            else if (iv==2) then
            io_node = 6000+mynode
            write (io_node,*) 'Re/Im var3d_oa iv2 =',j,tab3d_re
     &                                                ,tab3d_im
            end if
        endif
        return
      end subroutine trace2_oa_ik

!     function get_vnam( vnamfmt='(a5,a5,i3.3,a1,i3.3,a5)'
!    &                  ,vnamdp='3d_r_',vnams='_real', verbose=6)
      function get_vnam( vnamfmt
     &                  ,vnamdp,vnams,verbose)

        implicit none
        character(len=*), intent(in) :: vnamfmt, vnamdp, vnams
        integer, intent(in)           :: verbose
        !character(len=17)  intent(in) :: vnam_oa 
        !character(len=22)  intent(in) :: vnam_oa_r, vnam_oa_i 
        character(len=22)  :: get_vnam, tmpvnam !, tmpfmt 
        

!            write( vnam_oa_r, fmt='(a5,a5,i3.3,a1,i3.3,a5)') 
!& tgvnam_oa(tv_oa(iv)),'3d_r_',ic,'_',iv,'_real'
!            write( vnam_oa_i, fmt='(a5,a5,i3.3,a1,i3.3,a5)') 
!& tgvnam_oa(tv_oa(iv)),'3d_r_',ic,'_',iv,'_imag'

!            write( vnam_oa, fmt='(a5,a5,i3.3,a1,i3.3)') 
!& tgvnam_oa(tv_oa(iv)),'3d_r_',ic,'_',iv

        !tmpfmt=trim(vnamfmt)
        write( tmpvnam, fmt=trim(vnamfmt)) 
     &         tgvnam_oa(tv_oa(iv)),trim(vnamdp),ic,'_',iv,trim(vnams)
       
        if (verbose >=16 ) then
!        MPI_master_only write(stdout,fmt='(2(22a))') 
!    &   vnam_oa_r,vnam_oa_i
!        MPI_master_only write(stdout,fmt='(17a)') 
!    &   vnam_oa

         MPI_master_only write(stdout,fmt='(a22)') 
     &   tmpvnam
         MPI_master_only write(stdout,fmt='(i3.3)') 
     &   tvar_oa(tc_oa(ic),tvc_oa(ic),ivc)
         MPI_master_only write(stdout,fmt='(3(i3.3,1x))') 
     &   tc_oa(ic), tvc_oa(ic), ivc
        end if

        get_vnam = tmpvnam 

        return
      end function get_vnam


      
      end subroutine var3d_oa_out_tile

#else /* ONLINE_ANALYSIS */

      subroutine var3d_oa_out_tile_empty
      return
      end subroutine var3d_oa_out_tile_empty

#endif /* ONLINE_ANALYSIS */
