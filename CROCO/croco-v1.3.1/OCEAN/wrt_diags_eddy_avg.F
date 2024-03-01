!
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
! 
! ROMS_AGRIF website : http://www.romsagrif.org
!======================================================================
!
#include "cppdefs.h"
#if defined DIAGNOSTICS_EDDY && defined AVERAGES && ! defined XIOS
!
!---------------------------------------------------------------
! Write diagnostics fields at requested levels into diagnostics 
! netCDF file.

!---------------------------------------------------------------
!
      subroutine wrt_diags_eddy_avg

      implicit none
      integer ierr,  record, lstr, lvar, lenstr
     &   , start(2), count(2), ibuff(4), nf_fwrite, type
#if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE), blank
#endif
#include "param.h"
#include "scalars.h"
#include "ncscrum.h"
#include "forces.h"
#include "grid.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "mixing.h"
#include "diags_eddy.h"
#include "mpi_cpl.h"
#ifdef SEDIMENT
# include "sediment.h"
#endif
#ifdef BBL
# include "bbl.h"
#endif
#ifdef SOLVE3D
      integer tile,itrc,i,j,k,ivar
# ifdef SEDIMENT
     &      , indxWrk
# endif
# include "work.h"
#endif
#include "netcdf.inc"
#if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1, 
     &                 1, MPI_COMM_WORLD, status, ierr) 
      endif
#endif
!
#undef DEBUG
!
!
! Create/open diagnostic file; write grid arrays, if so needed.
!
      call def_diags_eddy_avg (nciddiags_eddy_avg, nrecdiags_eddy_avg, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(diags_eddyname_avg)
!                                          !!! WARNING: Once time
! Set record within the file.              !!! stepping has been 
!     
                                           !!! started, it is assumed
      if (iic.eq.0) nrecdiags_eddy_avg=nrecdiags_eddy_avg+1  !!! that global history
      if (nrpfdiags_eddy_avg.eq.0) then            !!! record index "nrecdiags_eddy"
        record=nrecdiags_eddy_avg                  !!! is advanced by main.
      else
        record=1+mod(nrecdiags_eddy_avg-1, nrpfdiags_eddy_avg)
      endif
!
!---------------------------------------------------------------
! Write out evolving model variables:
!---------------------------------------------------------------
!
! Time step number and record numbers.
!
      type=filetype_diags_eddy_avg
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrecdiags_eddy_avg
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiags_eddy_avg, diags_eddyTstep_avg,
     &                       start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step_avg', record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_FTYPE (nciddiags_eddy_avg, diags_eddyTime_avg,
     &                       record, timediags_eddy_avg)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time2
!
      ierr=nf_put_var1_FTYPE (nciddiags_eddy_avg, diags_eddyTime2_avg,
     &                       record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime2))
        write(stdout,1) vname(1,indxTime2)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif

!---------------------------------------------------------------
! PV diagnostic variables.
!---------------------------------------------------------------
!
!

         itrc=1

       if (wrtdiags_eddy_avg(1)) then
!---------------------------------------------------------------
! indxeddyubu
!
        call fillvalue2d(eddyubu_avg,nciddiags_eddy_avg,diags_eddyubu_avg(1),
     &                      indxeddyubu,
     &                      record,u2dvar,type)

!---------------------------------------------------------------
! indxeddyvbv
!
        call fillvalue2d(eddyvbv_avg,nciddiags_eddy_avg,diags_eddyvbv_avg(1),
     &                      indxeddyvbv,
     &                      record,v2dvar,type)
!---------------------------------------------------------------
! indxeddyusu
!
        call fillvalue2d(eddyusu_avg,nciddiags_eddy_avg,diags_eddyusu_avg(1),
     &                      indxeddyusu,
     &                      record,u2dvar,type)

!---------------------------------------------------------------
! indxeddyvsv
!
        call fillvalue2d(eddyvsv_avg,nciddiags_eddy_avg,diags_eddyvsv_avg(1),
     &                      indxeddyvsv,
     &                      record,v2dvar,type)
!---------------------------------------------------------------
! indxeddyugsu
!
        call fillvalue2d(eddyugsu_avg,nciddiags_eddy_avg,diags_eddyugsu_avg(1),
     &                      indxeddyugsu,
     &                      record,u2dvar,type)

!---------------------------------------------------------------
! indxeddyvgsv
!
        call fillvalue2d(eddyvgsv_avg,nciddiags_eddy_avg,diags_eddyvgsv_avg(1),
     &                      indxeddyvgsv,
     &                      record,v2dvar,type)
!---------------------------------------------------------------
! indxeddyzz
!
        call fillvalue2d(eddyzz_avg,nciddiags_eddy_avg,diags_eddyzz_avg(1),
     &                      indxeddyzz,
     &                      record,r2dvar,type)

!---------------------------------------------------------------
! indxeddyuu
!
        call fillvalue3d(eddyuu_avg,nciddiags_eddy_avg,diags_eddyuu_avg(1),
     &                      indxeddyuu,
     &                      record,r3dvar,type)

!---------------------------------------------------------------
! indxeddyvv
!
        call fillvalue3d(eddyvv_avg,nciddiags_eddy_avg,diags_eddyvv_avg(1),
     &                      indxeddyvv,
     &                      record,r3dvar,type)
!---------------------------------------------------------------
! indxeddyuv
!
        call fillvalue3d(eddyuv_avg,nciddiags_eddy_avg,diags_eddyuv_avg(1),
     &                      indxeddyuv,
     &                      record,r3dvar,type)
!---------------------------------------------------------------
! indxeddyuw
!
        call fillvalue3d(eddyuw_avg,nciddiags_eddy_avg,diags_eddyuw_avg(1),
     &                      indxeddyuw,
     &                      record,r3dvar,type)
!---------------------------------------------------------------
! indxeddyvw
!
        call fillvalue3d(eddyvw_avg,nciddiags_eddy_avg,diags_eddyvw_avg(1),
     &                      indxeddyvw,
     &                      record,r3dvar,type)
!---------------------------------------------------------------
! indxeddywb
!
        call fillvalue3d(eddywb_avg,nciddiags_eddy_avg,diags_eddywb_avg(1),
     &                      indxeddywb,
     &                      record,r3dvar,type)
!---------------------------------------------------------------
! indxeddyub
!
        call fillvalue3d(eddyub_avg,nciddiags_eddy_avg,diags_eddyub_avg(1),
     &                      indxeddyub,
     &                      record,r3dvar,type)
!---------------------------------------------------------------
! indxeddyvb
!
        call fillvalue3d(eddyvb_avg,nciddiags_eddy_avg,diags_eddyvb_avg(1),
     &                      indxeddyvb,
     &                      record,r3dvar,type)
!---------------------------------------------------------------
       endif

!


  1   format(/1x,'WRT_DIAGS_EDDY_AVG ERROR while writing variable ''', A,
     &          ''' into diag file.', /11x, 'Time record:',
     &                   I6,3x,'netCDF error code',i4,3x,a,i4)
      goto 100 
  99  may_day_flag=3
 100  continue

!
!---------------------------------------------------------------
! Synchronize netCDF file to disk to allow other processes
! to access data immediately after it is written.
!---------------------------------------------------------------
!
#if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      ierr=nf_close (nciddiags_eddy_avg)
      if (nrpfdiags_eddy_avg.gt.0 .and. record.ge.nrpfdiags_eddy_avg)
     &            nciddiags_eddy_avg=-1
#else
      if (nrpfdiags_eddy_avg.gt.0 .and. record.ge.nrpfdiags_eddy_avg) then
        ierr=nf_close (nciddiags_eddy_avg)
        nciddiags_eddy_avg=-1
      else
        ierr=nf_sync(nciddiags_eddy_avg)
      endif
#endif


      if (ierr .eq. nf_noerr) then
        MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')
     &            'WRT_DIAGS_EDDY_AVG -- wrote',
     &            ' diag fields into time record =', record, '/',
     &             nrecdiags_eddy_avg  MYID
      else
        MPI_master_only write(stdout,'(/1x,2A/)') 
     &             'WRT_DIAGS_EDDY_AVG ERROR: Cannot ',
     &             'synchronize/close diag netCDF file.'
        may_day_flag=3
      endif

#if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode .lt. NNODES-1) then
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &                        1, MPI_COMM_WORLD,  ierr)
      endif
#endif
      return
      end

#else
      subroutine wrt_diags_eddy_avg_empty
      end
#endif /* (DIAGNOSTICS_EDDY) */
