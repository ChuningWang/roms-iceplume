! $Id: output_oa.F
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
!! - Outputs the online spectral and wavelet analyses with XIOS
!!
!! @details written to support tile-threads and 
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
!
#include "cppdefs.h"
! BLXD includes :
! - cppdefs_dev.h
! - set_global_definitions.h
!   with GLOBAL_2D_ARRAY
!        PRIVATE_SCRATCH_2D_ARRAY... 
#ifdef ONLINE_ANALYSIS

      subroutine output_oa
!
!=======================================================================
!                                                                      !
! Output for OA Online analysis block
!                                                                      !
!=======================================================================
!
      use xios_module
      use module_interface_oa, only : tvar_oa
     & ,iscal0d_cr, jscal0d_cr, kscal0d_cr, per0d_cr
#ifdef MPI
     & ,iscal0d_oa, jscal0d_oa, kscal0d_oa, per0d_oa
#endif
     & ,tgv3d_oa, tgv_oa, tgvnam_oa, nzvc_oa, nzupd0d_oa
     & ,verbose_oa, io_hist, directory_out_oa
     & ,scalogram_analysis
     & ,nsclg_loc, if_mpi_oa
     & ,if_record_sclg_ijpoints
      use module_oa_variables,  only : tv_oa     
     & ,tupd_oa, tv_sclg_oa
      use module_oa_periode,    only : begc_oa, tvc_oa, tc_oa, nzc_oa
     & ,swt_wfpf_oa, nper_sclg_max

      implicit none
!
#include "param.h"
#include "scalars.h"
#include "work.h"
!
                integer :: ntrds,trd,range,
     &        my_first,my_last, tile, ierr

      integer :: omp_get_num_threads, omp_get_thread_num

      character(len=22) :: vnamr, vnami, pnam 

#ifdef OA_WORK
      real(8) :: workre(GLOBAL_2D_ARRAY,N)
#endif
      real(8) :: workim(GLOBAL_2D_ARRAY,N)

      character(len=9) :: vvnamr, vvnami
      character(len=5) :: ppnam 
      character(len=7) :: inam, jnam, knam 

      real(8) :: ww_im(nper_sclg_max, nzupd0d_oa)
      real(8) :: ww_re(nper_sclg_max, nzupd0d_oa)
      real(8) :: pper (nper_sclg_max, nzupd0d_oa)

      real(8) :: wl_im(nper_sclg_max, nsclg_loc)
      real(8) :: wl_re(nper_sclg_max, nsclg_loc)
      real(8) :: perl (nper_sclg_max, nsclg_loc)

      logical :: if_imag_part
      integer :: ic, iv, ivc, isclg, ilc

      ilc=1+iic-ntstart   ! number of time step since rst
      ! Must be > 1 to get 'once' variables 
      call xios_update_calendar(max(ilc-1,1))

! Univoque relationship
! ONE isclg SCALOGRAM -> ONE iv variable
! But 
! ONE iv variable -> not sure this is a saclogram
! - if MPI handled by XIOS
!   complex 2D array scal0d_oa( 1: nper_sclg_max, v2locsclg(iv))
! - if no MPI OR MPI handled by OA module 
! - complex 2D array scal0d_cr( 1: nper_sclg_max, tupd_oa(iv))

      oa_config_loop : do ic = 1, nzc_oa
      
        ivc = 0
        do iv = begc_oa(ic), begc_oa(ic+1) -1
            ivc = ivc + 1
        end do
        oa_stop1 : if ( ivc /= nzvc_oa( tc_oa(ic) ) ) then
           MPI_master_only write(stdout,fmt='(a)') '...OA-ERR-iv-ic'
           stop
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
      
#ifdef OA_TRACES && ( MILES || IGW )
        call traces_oa1(verbose_oa)
#endif

! 3Dvar
        if_is_3dvar : if ( tgv3d_oa( tv_oa(iv) ) == 3 ) then

!
! XIOS Case only
!
#ifdef XIOS
          
          v3d_not_scalogram : if ( .not. tv_sclg_oa(iv) ) then 

              ntrds=omp_get_num_threads()
              trd=omp_get_thread_num()
              range=(NSUB_X*NSUB_E+ntrds-1)/ntrds
              my_first=trd*range
              my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)

              do tile=my_first,my_last
#ifdef OA_WORK
                call var3d_oa_out(tile,ic,iv,workre,workim,vnamr,vnami)
#else
                ! BLXD TODO test with common space work3d_r and workr
#ifdef OA_TRACES              
              MPI_master_only write(stdout,fmt='(a)') 
     &          '...OA var3d outputs using croco workr'
#endif              
                call var3d_oa_out(tile,ic,iv,workr,workim,vnamr,vnami)
#endif
              end do
C$OMP BARRIER
#ifdef OA_WORK
              call xios_send_field( trim(vnamr), real(workre)  )
#else
              call xios_send_field( trim(vnamr), real(workr)   )
#endif
              if ( swt_wfpf_oa(iv)==4 )
     &          call xios_send_field( trim(vnami), real(workim) )
          endif v3d_not_scalogram

#else
           MPI_master_only write(stdout,fmt='(a)') 
     &      '...Non XIOS outputs doesn''t work yet with OA'
           stop
#endif

! 2Dvar
        else if ( tgv3d_oa( tv_oa(iv) ) == 2 ) then if_is_3dvar

#ifdef XIOS
        v2d_not_scalogram : if ( .not. tv_sclg_oa(iv) ) then 

#ifdef OA_TRACES         
              MPI_master_only write(stdout,fmt='(a)') 
     &          '...OA var2d outputs using croco work2d/work2d2'
#endif              
              ntrds=omp_get_num_threads()
              trd=omp_get_thread_num()
              range=(NSUB_X*NSUB_E+ntrds-1)/ntrds
              my_first=trd*range
              my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)
              
              do tile=my_first,my_last
                  call var2d_oa_out(tile,ic,iv,work2d,work2d2,vnamr,vnami)
              end do

C$OMP BARRIER
              call xios_send_field( trim(vnamr), real(work2d)  )
              if ( swt_wfpf_oa(iv)==4 )
     &          call xios_send_field( trim(vnami), real(work2d2) )

          endif v2d_not_scalogram

#else
           MPI_master_only write(stdout,fmt='(a)') 
     &      '...Non XIOS outputs doesn''t work with OA'
           stop
#endif

        else if_is_3dvar
           MPI_master_only write(stdout,fmt='(a)') 
     &      '...ERROR output_oa : unrecogn. var. dimension'
           stop
        endif if_is_3dvar
!
        end do the_oa_var_loop
      end do oa_config_loop
!   !
        if_output_scalogram : if ( scalogram_analysis ) then

         if_nonloc_sclg : if ( if_mpi_oa ) then

              call scal0d_oa_out_full( ww_re,ww_im
     &                                ,vvnamr,vvnami
     &                                ,pper,ppnam
     &                                ,inam,jnam,knam
     &                                ,if_imag_part)

#ifdef OA_TRACES_XIOS             
                  write(*,fmt='(a,1(i5))') 
     &            '...OA-Scalogram : send re to XIOS subdom',mynode
#endif
              call xios_send_field( trim(vvnamr), 
     &                  real(ww_re(1:nper_sclg_max, 1:nzupd0d_oa)))
              if ( if_imag_part ) then
#ifdef OA_TRACES_XIOS             
                  write(*,fmt='(a,1(i5))') 
     &            '...OA-Scalogram : send im to XIOS subdom',mynode
#endif
              call xios_send_field( trim(vvnami), 
     &                  real(ww_im(1:nper_sclg_max, 1:nzupd0d_oa)))
              endif
#ifdef OA_TRACES_XIOS             
                  write(*,fmt='(a,1(i5))') 
     &            '...OA-Scalogram : send period to XIOS subdom',mynode
#endif
              call xios_send_field( trim(ppnam),
     &                  real(pper(1:nper_sclg_max, 1:nzupd0d_oa)))
!    &                    per0d_cr(1:nper_sclg_max, 1:nzupd0d_oa))

              if ( if_record_sclg_ijpoints ) then

#ifdef OA_TRACES_XIOS             
                  write(*,fmt='(a,1(i5))')
     &        '...OA-Scalogram : send i-pt to XIOS subdom',mynode
#endif
              call xios_send_field( trim(inam),
     &              real(iscal0d_cr(1:nzupd0d_oa)))
#ifdef OA_TRACES_XIOS             
                  write(*,fmt='(a,1(i5))')
     &        '...OA-Scalogram : send j-pt to XIOS subdom',mynode
#endif
              call xios_send_field( trim(jnam),
     &              real(jscal0d_cr(1:nzupd0d_oa)))
#ifdef OA_TRACES_XIOS             
                  write(*,fmt='(a,1(i5))')
     &        '...OA-Scalogram : send k-pt to XIOS subdom',mynode
#endif
              call xios_send_field( trim(knam),
     &              real(kscal0d_cr(1:nzupd0d_oa)))
#ifdef OA_TRACES_XIOS             
              write(*,fmt='(a,1(i5))')
     &        '...OA-Scalogram : sending END XIOS subdom',mynode
#endif
              endif

         else if_nonloc_sclg
#ifdef MPI              
              ! if nsclg_loc=0 isclg=0 => array range 0:0
              isclg=min(1,nsclg_loc)
              !if ( nsclg_loc > 0 ) then
              call scal0d_oa_out_loc(  wl_re,wl_im
     &                                ,vvnamr,vvnami
     &                                ,perl,ppnam
     &                                ,inam,jnam
     &                                ,knam
     &                                ,if_imag_part)
              !scal0d_oa_out_loc returns 
              ! - wl_*,perl arrays with zeros if nsclg_loc==0
              ! - if_imag_part true if ImPart updating necessary 
              !else
              !  wl_re(1:nper_sclg_max,0) = 0. 
              !  wl_im(1:nper_sclg_max,0) = 0.
              !  perl (1:nper_sclg_max,0) = 0.
              !endif
              loc_scal : if ( nsclg_loc>0) then

#ifdef OA_TRACES_XIOS             
                  write(*,fmt='(a,1(i5))') 
     &            '..OA-Scalogram-loc : send re to XIOS subdom',mynode
#endif
              call xios_send_field( trim(vvnamr), 
     &                  real(wl_re(1:nper_sclg_max, isclg:nsclg_loc)))

     !CRBS         write(stdout,fmt='(2i,1(ES22.15E2))')mynode,                !CRBE
     !CRBS &                                             isclg,wl_re(1,isclg)  !CRBE
              if ( if_imag_part ) then
#ifdef OA_TRACES_XIOS             
                  write(*,fmt='(a,1(i5))') 
     &            '..OA-Scalogram-loc : send im to XIOS subdom',mynode
#endif
              call xios_send_field( trim(vvnami), 
     &                  real(wl_im(1:nper_sclg_max, isclg:nsclg_loc)))
              endif
#ifdef OA_TRACES_XIOS             
                  write(*,fmt='(a,1(i5))') 
     &            '..OA-Scalogram : send period to XIOS subdom',mynode
#endif
              call xios_send_field( trim(ppnam),
     &                  real(perl(1:nper_sclg_max, isclg:nsclg_loc)))
     !&                    per0d_oa(1:nper_sclg_max, isclg:nsclg_loc))

              if (if_record_sclg_ijpoints) then
#ifdef OA_TRACES_XIOS             
                  write(*,fmt='(a,1(i5))') 
     &            '..OA-Scalogram : send ij-pts to XIOS subdom',mynode
#endif
              call xios_send_field( trim(inam),
     &                     real(iscal0d_oa(isclg:nsclg_loc)))
              call xios_send_field( trim(jnam),
     &                     real(jscal0d_oa(isclg:nsclg_loc)))
              call xios_send_field( trim(knam),
     &                     real(kscal0d_oa(isclg:nsclg_loc)))
              endif
#ifdef OA_TRACES_XIOS             
                  write(*,fmt='(a,1(i5))') 
     &            'OA-Scalogram-loc : sending END XIOS subdom',mynode
#endif
              else loc_scal

#ifdef OA_TRACES_XIOS             
                  write(*,fmt='(a,1(i5))') 
     &            '..OA-Scalogram-loc :send re40 to XIOS subdo',mynode
#endif
              call xios_send_field( trim(vvnamr), 
     &                  RESHAPE(
     &                  real(wl_re(1:nper_sclg_max, isclg:nsclg_loc))
     &                  ,(/nper_sclg_max,0/) ) )
 !CRBS #ifdef OA_TRACES_XIOS                                                !CRBE 
 !CRBS             write(*,fmt='(2i,1(ES22.15E2))')mynode,             !CRBE
 !CRBS    &                                             isclg,wl_re(1,isclg)!CRBE
 !CRBS #endif                                                               !CRBE
              if ( if_imag_part ) then
#ifdef OA_TRACES_XIOS             
                  write(*,fmt='(a,1(i5))') 
     &            '..OA-Scalogram-loc :send im40 to XIOS subdo',mynode
#endif
              call xios_send_field( trim(vvnami), 
     &                  RESHAPE(
     &                  real(wl_im(1:nper_sclg_max, isclg:nsclg_loc))
     &                  ,(/nper_sclg_max,0/) ) )

              endif
#ifdef OA_TRACES_XIOS             
                  write(*,fmt='(a,1(i5))') 
     &            '..OA-Scalogram :send period40 to XIOS subdo',mynode
#endif
              call xios_send_field( trim(ppnam),
     &                  RESHAPE(
     &                  real(perl(1:nper_sclg_max, isclg:nsclg_loc))
     &                  ,(/nper_sclg_max,0/) ) )

              if (if_record_sclg_ijpoints) then
#ifdef OA_TRACES_XIOS             
                  write(*,fmt='(a,1(i5))') 
     &            '..OA-Scalogram :send ij-pts40 to XIOS subdo',mynode
#endif
              call xios_send_field( trim(inam),
     &                     REAL( RESHAPE(
     &                     iscal0d_oa(isclg:nsclg_loc)
     &                    ,(/0/) ) ))
              call xios_send_field( trim(jnam),
     &                     REAL( RESHAPE(
     &                     jscal0d_oa(isclg:nsclg_loc)
     &                    ,(/0/) ) ))
              call xios_send_field( trim(knam),
     &                     REAL( RESHAPE(
     &                     kscal0d_oa(isclg:nsclg_loc)
     &                    ,(/0/) ) ))
              endif
#ifdef OA_TRACES_XIOS             
                  write(*,fmt='(a,1(i5))') 
     &            'OA-Scalogram-loc : sending END40 XIOS subdo',mynode
#endif
              endif loc_scal
#else /* MPI */
                  write(stdout,fmt='(a)') 
     &            '..OA-Sclg : XIOS disctribut. axis impossible'
#endif /* MPI */
         endif if_nonloc_sclg

        endif if_output_scalogram

      contains

!      subroutine traces_oa1(verbose,iic)        
      subroutine traces_oa1(verbose)

        implicit none
        integer, intent(in)    :: verbose
        if ( verbose >=15 ) then
            MPI_master_only write(stdout,fmt='(a13,(3(i3.3,1x)))')     
     &     'Cfg-#var-Var ',ic,ivc,iv
        end if
        return
      end subroutine traces_oa1



      end subroutine output_oa

#else /* ONLINE_ANALYSIS */

      subroutine output_oa_empty
      return
      end subroutine output_oa_empty

#endif /* ONLINE_ANALYSIS */

