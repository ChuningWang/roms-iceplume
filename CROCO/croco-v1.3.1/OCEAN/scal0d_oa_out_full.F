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
!! @details called in output_oa.F, scalogram only
!!
!> @authors  
!! - B. Lemieux-Dudon
!!  - based on a preliminary Croco-OA interface version : spring 2020.
!! - More history (authors, comments) in source module_interface_oa.F90
!> @todo 
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

      subroutine scal0d_oa_out_full ( 
     &                               scal0d_re, scal0d_im
     &                             , vnam_oa_r, vnam_oa_i
     &                             , per0d
     &                             , pnam_oa
     &                             , inam_oa, jnam_oa
     &                             , knam_oa
     &                             , if_imag_part )

     !& ,scal0d_oa, scal0d_cr,verbose_oa, io_hist, directory_out_oa
     !& ,scalogram_analysis

      use module_interface_oa, only : scal0d_cr, per0d_cr
     & ,tvar_oa , tgv_oa, tgvnam_oa, nzvc_oa, nzupd0d_oa
     & ,verbose_oa
      use module_oa_variables,  only : tv_oa     
     & ,tupd_oa, tv_sclg_oa
      use module_oa_periode,    only : begc_oa, tvc_oa, tc_oa
     & ,nzc_oa, swt_wfpf_oa, nper_sclg_max, nper_sclg
      !use module_oa_space,      only : ptij_oa
      use scalars
      implicit none
! >> Include with variable declaration
! # include "param.h"
      real(8), intent(inout) :: scal0d_re(1:nper_sclg_max,nzupd0d_oa)
     &                         ,scal0d_im(1:nper_sclg_max,nzupd0d_oa)
     &                         ,per0d(1:nper_sclg_max,nzupd0d_oa)
      logical, intent(inout) :: if_imag_part
      character(len=9),  intent(inout)  :: vnam_oa_r, vnam_oa_i 
      character(len=5), intent(inout)   :: pnam_oa 
      character(len=7), intent(inout)   :: inam_oa, jnam_oa
     &                                    ,knam_oa

! Local variables
      integer :: iv, ic
      integer :: ip, ivc, io_node
      integer, parameter :: verbl=6

      if_imag_part = .false.
      scal0d_re(1:nper_sclg_max,1:nzupd0d_oa) = 0. 
      scal0d_im(1:nper_sclg_max,1:nzupd0d_oa) = 0.
      per0d    (1:nper_sclg_max,1:nzupd0d_oa) = 0.


      oa_config_loop : do ic = 1, nzc_oa
      
        ivc = 0
        do iv = begc_oa(ic), begc_oa(ic+1) -1
            ivc = ivc + 1
        end do
        oa_stop1 : if ( ivc /= nzvc_oa( tc_oa(ic) ) ) then
           MPI_master_only write(stdout,fmt='(a)') '...OA-ERR-iv-ic'
           stop
           !MPI_master_only write(stdout,*) 'Expected total # of var nzvc_oa( tc_oa(ic) ) ',nzvc_oa( tc_oa(ic) )
        end if oa_stop1

        !could enable to retrieve a name, requires to use tupd_oa
        ! # of variable per configuration ic
        ! do iv = 1, nzvc_oa( tc_oa(ic) )
        ivc = 0
        the_oa_var_loop : do iv = begc_oa(ic), begc_oa(ic+1) -1
           ivc = ivc + 1
        oa_stop2 : if ( 
     &     tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ) /= tupd_oa(iv) ) 
     &     then
           MPI_master_only write(stdout,fmt='(a)') '...OA-ERR tupd_oa'
           stop
        end if oa_stop2

        !ivc = nzvc_oa( tc_oa(ic) )

#ifdef OA_TRACES && ( MILES || IGW)
        call traces_oa1(verbose_oa)
#endif

          is_a_scalogram : if ( tv_sclg_oa(iv) ) then 

                  ! iv <-> ONE SCALOGRAM 
                  ! - complex 2D array scal0d_oa( 1: nper_sclg_max, tupd_oa(iv))
                  ! - dumped in 1D array w_re(1: nper_sclg_max)


        if_oa_re_im : if ( swt_wfpf_oa(iv)==4 ) then

      ! #BLXD 2020 REAL/DIMAG for real/imag part of dble prec scal0d_oa 
      ! Warning mpc.F parser REAL -> REAL real_type
      ! Warning : different gnu, intel standards 

        ! If one single variable is imag all MUST be
        if_imag_part = .true.

        loop_per_re_im : do ip=1,nper_sclg(iv)
#ifdef OA_TRACES_XIOS             
            MPI_master_only write(stdout,fmt='(a,2(i3))') 
     &      '...OA-Scalogram Part/ImPart for period ',ip,iv
            write(*,fmt='(a,4(i3))') 
     &      '>><<OA-Scalogram node/ic/iv/ip '
     &      ,mynode,ic,iv,ip
#endif
                       scal0d_re(ip,tupd_oa(iv)) = 
     & REAL(
     & DBLE( scal0d_cr( ip, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ))))

                       scal0d_im(ip,tupd_oa(iv)) = 
     & REAL(
     & DIMAG( scal0d_cr( ip, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ))))


        enddo loop_per_re_im

            ! scal_ + real or imag
            vnam_oa_r = get_vnam_full( vnamfmt='(a5,a4)'
     &              ,vnams='real',verbose=verbose_oa)
            vnam_oa_i = get_vnam_full( vnamfmt='(a5,a4)'
     &              ,vnams='imag',verbose=verbose_oa)

        else if_oa_re_im

        loop_per_re : do ip=1,nper_sclg(iv)
                       scal0d_re(ip,tupd_oa(iv)) = 
     &  scal0d_cr(ip, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ))
        enddo loop_per_re

            ! scal_ + real or imag
            vnam_oa_r = get_vnam_full( vnamfmt='(a5,a4)'
     &              ,vnams='',verbose=verbose_oa)

        endif if_oa_re_im

        loop_for_per_axis : do ip=1,nper_sclg(iv)

                       per0d(ip,tupd_oa(iv)) = 
     &        per0d_cr( ip, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ))

#ifdef OA_TRACES_XIOS             
            MPI_master_only write(stdout,fmt='(a,i3,f16.2)') 
     &      '...OA-Scalogram Period axis ',ip,per0d(ip,tupd_oa(iv))
#endif

#ifdef OA_TRACES && ( MILES || IGW)
                       if(verbose_oa>=verbl) then
                       call traces2_oa_scal0d(ip) 
                       !call traces_oa_scal0d(ip) 
                       endif
#endif

        enddo loop_for_per_axis


       pnam_oa = get_pnam_full( pnamfmt='(a5)'
     &              ,pnamdp='perv_',verbose=verbose_oa )

       inam_oa = get_vnam_full( vnamfmt='(a5,a1)'
     &              ,vnams='i',verbose=verbose_oa)
       jnam_oa = get_vnam_full( vnamfmt='(a5,a1)'
     &              ,vnams='j',verbose=verbose_oa)
       knam_oa = get_vnam_full( vnamfmt='(a5,a1)'
     &              ,vnams='k',verbose=verbose_oa)

          endif is_a_scalogram

        enddo the_oa_var_loop

      enddo oa_config_loop

      contains

      subroutine traces_oa1(verbose)        

        implicit none
        integer, intent(in)    :: verbose
        if ( verbose >=6 ) then
            MPI_master_only write(stdout,fmt='(a13,(3(i3.3,1x)))')     
     &     'Cfg-#var-Var ',ic,ivc,iv
        end if
        return
      end subroutine traces_oa1

      subroutine traces_oa_scal0d(ip)        
        implicit none
        integer, intent(in) :: ip

            if (ip<=9) then

                io_node = (10+ip)*1000+mynode
                write (io_node,fmt='(a,i3.3,2(1x,ES22.15E2))') 
     &          'ip, scal0d_cr = '
     &         ,ip
     &         ,REAL(DBLE ( scal0d_cr(ip,tupd_oa(iv)) ))
     &         ,REAL(DIMAG( scal0d_cr(ip,tupd_oa(iv)) ))
               
               io_node = (20+ip)*1000+mynode
               write (io_node,fmt='(a,i3.3,2(1x,ES22.15E2))') 
     &         'ip, scal0d_re/im = '
     &         ,ip
     &         ,scal0d_re(ip,tupd_oa(iv))
     &         ,scal0d_im(ip,tupd_oa(iv))

            endif
        return
      end subroutine traces_oa_scal0d

      subroutine traces2_oa_scal0d(ip)        
        implicit none
        integer, intent(in) :: ip

            if (ip<=9) then

                io_node = (10+ip)*1000+mynode
                write (io_node,fmt='(a,f16.2,2(1x,ES22.15E2))') 
     &          'ip, scal0d_cr = '
     &         ,per0d_cr(ip,tupd_oa(iv))
     &         ,REAL(DBLE ( scal0d_cr(ip,tupd_oa(iv)) ))
     &         ,REAL(DIMAG( scal0d_cr(ip,tupd_oa(iv)) ))
               
               io_node = (20+ip)*1000+mynode
               write (io_node,fmt='(a,f16.2,2(1x,ES22.15E2))') 
     &         'ip, scal0d_re/im = '
     &         ,per0d    (ip,tupd_oa(iv))
     &         ,scal0d_re(ip,tupd_oa(iv))
     &         ,scal0d_im(ip,tupd_oa(iv))

            endif
        return
      end subroutine traces2_oa_scal0d

      function get_vnam_full( vnamfmt, vnams,verbose)

        implicit none
        character(len=*), intent(in) :: vnamfmt, vnams
        integer, intent(in)           :: verbose
        character(len=9)  :: get_vnam_full, tmpvnam !, tmpfmt 
        
        write( tmpvnam, fmt=trim(vnamfmt)) 
     &         tgvnam_oa(tv_oa(iv)),trim(vnams)

#ifdef OA_TRACES       
        if (verbose >=16 ) then
         MPI_master_only write(stdout,fmt='(a9)') 
     &   tmpvnam
         MPI_master_only write(stdout,fmt='(i3.3)') 
     &   tvar_oa(tc_oa(ic),tvc_oa(ic),ivc)
         MPI_master_only write(stdout,fmt='(3(i3.3,1x))') 
     &   tc_oa(ic), tvc_oa(ic), ivc
        end if
#endif

        get_vnam_full = tmpvnam 

        return
      end function get_vnam_full

      function get_pnam_full( pnamfmt
     &                  ,pnamdp,verbose)

        implicit none
        character(len=*), intent(in) :: pnamfmt, pnamdp
        integer, intent(in)           :: verbose
        character(len=5)  :: get_pnam_full, tmppnam !, tmpfmt 

        write( tmppnam, fmt=trim(pnamfmt)) 
     &         trim(pnamdp)

        get_pnam_full = tmppnam 

      end function get_pnam_full

      function get_pnam( pnamfmt
     &                  ,pnamdp,verbose)

        implicit none
        character(len=*), intent(in) :: pnamfmt, pnamdp
        integer, intent(in)           :: verbose
        character(len=22)  :: get_pnam, tmppnam !, tmpfmt 

        write( tmppnam, fmt=trim(pnamfmt)) 
     &         trim(pnamdp),ic,'_',iv

        get_pnam = tmppnam 

      end function get_pnam

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
       
#ifdef OA_TRACES
        if (verbose >=16 ) then
         MPI_master_only write(stdout,fmt='(a22)') 
     &   tmpvnam
         MPI_master_only write(stdout,fmt='(i3.3)') 
     &   tvar_oa(tc_oa(ic),tvc_oa(ic),ivc)
         MPI_master_only write(stdout,fmt='(3(i3.3,1x))') 
     &   tc_oa(ic), tvc_oa(ic), ivc
        end if
#endif

        get_vnam = tmpvnam 

        return
      end function get_vnam
      
      end subroutine scal0d_oa_out_full 

#else /* ONLINE_ANALYSIS */

      subroutine scal0d_oa_out_full_empty
      return
      end subroutine scal0d_oa_out_full_empty

#endif /* ONLINE_ANALYSIS */
