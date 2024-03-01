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
#if defined DIAGNOSTICS_PV
!
!---------------------------------------------------------------
! Write diagnostics fields at requested levels into diagnostics 
! netCDF file.

!---------------------------------------------------------------
!
      subroutine wrt_diags_pv

      implicit none
      integer ierr,  record, lstr, lvar, lenstr
     &   , start(2), count(2), ibuff(4), nf_fwrite, type
#if defined MPI & !defined PARALLEL_FILES
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
#include "diags_pv.h"
#include "diagnostics.h"
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
      call def_diags_pv (nciddiags_pv, nrecdiags_pv, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(diags_pvname)
!                                          !!! WARNING: Once time
! Set record within the file.              !!! stepping has been 
!     
                                           !!! started, it is assumed
      if (iic.eq.0) nrecdiags_pv=nrecdiags_pv+1  !!! that global history
      if (nrpfdiags_pv.eq.0) then            !!! record index "nrecdiags_pv"
        record=nrecdiags_pv                  !!! is advanced by main.
      else
        record=1+mod(nrecdiags_pv-1, nrpfdiags_pv)
      endif
!
!---------------------------------------------------------------
! Write out evolving model variables:
!---------------------------------------------------------------
!
! Time step number and record numbers.
!
      type=filetype_diags_pv
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrecdiags_pv
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiags_pv, diags_pvTstep,
     &                                            start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step', record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_FTYPE (nciddiags_pv, diags_pvTime, record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time2
!
      ierr=nf_put_var1_FTYPE (nciddiags_pv, diags_pvTime2, record, time)
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
      if (wrtdiags_pv(1)) then


# if defined DIAGNOSTICS_PV_FULL

         itrc=1


         call fillvalue3d(pv,nciddiags_pv,diags_pvpv(itrc),
     &                      indxpvpv,
     &                      record,u3dvar,type)


         call fillvalue3d(pvd,nciddiags_pv,diags_pvpvd(itrc),
     &                      indxpvpvd,
     &                      record,v3dvar,type)

# endif
!---------------------------------------------------------------

        do itrc=1,2

          if (itrc.eq.1) then
            ivar=u3dvar
          else
            ivar=v3dvar
          endif
!
! indxMrhs
!


           workr=Mrhs(:,:,:,itrc)
           call fillvalue3d(workr,nciddiags_pv,diags_pvMrhs(itrc),
     &                      indxpvMrhs+itrc-1,
     &                      record,ivar,type)  

        enddo
!---------------------------------------------------------------

      endif

!---------------------------------------------------------------

        do itrc=1,NTA
          if (wrtdiags_pv(itrc)) then

!
! indxTrhs
!

           workr=Trhs(:,:,:,itrc)
!           workr=TXadv(:,:,:,itrc)
           call fillvalue3d(workr,nciddiags_pv,diags_pvTrhs(itrc),
     &                      indxpvTrhs+itrc-1,
     &                      record,r3dvar,type)


          endif
        enddo
!---------------------------------------------------------------


  1   format(/1x,'WRT_DIAGS_PV ERROR while writing variable ''', A,
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
      ierr=nf_close (nciddiags_pv)
      if (nrpfdiags_pv.gt.0 .and. record.ge.nrpfdiags_pv) nciddiags_pv=-1
#else
      if (nrpfdiags_pv.gt.0 .and. record.ge.nrpfdiags_pv) then
        ierr=nf_close (nciddiags_pv)
        nciddiags_pv=-1
      else
        ierr=nf_sync(nciddiags_pv)
      endif
#endif


      if (ierr .eq. nf_noerr) then
        MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')
     &            'WRT_DIAGS_PV -- wrote',
     &            ' diag fields into time record =', record, '/',
     &             nrecdiags_pv  MYID
      else
        MPI_master_only write(stdout,'(/1x,2A/)') 
     &             'WRT_DIAGS_PV ERROR: Cannot ',
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
      subroutine wrt_diags_pv_empty
      end
#endif /* (DIAGNOSTICS_PV) */
