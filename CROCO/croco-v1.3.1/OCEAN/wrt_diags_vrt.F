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
#if defined DIAGNOSTICS_VRT 
!
!---------------------------------------------------------------
! Write diagnostics fields at requested levels into diagnostics 
! netCDF file.

!---------------------------------------------------------------
!
      subroutine wrt_diags_vrt

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
#include "diags_vrt.h"
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
      call def_diags_vrt (nciddiags_vrt, nrecdiags_vrt, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(diags_vrtname)
!                                          !!! WARNING: Once time
! Set record within the file.              !!! stepping has been 
!     
                                           !!! started, it is assumed
      if (iic.eq.0) nrecdiags_vrt=nrecdiags_vrt+1  !!! that global history
      if (nrpfdiags_vrt.eq.0) then            !!! record index "nrecdiags_vrt"
        record=nrecdiags_vrt                  !!! is advanced by main.
      else
        record=1+mod(nrecdiags_vrt-1, nrpfdiags_vrt)
      endif
!
!---------------------------------------------------------------
! Write out evolving model variables:
!---------------------------------------------------------------
!
! Time step number and record numbers.
!
      type=filetype_diags_vrt
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrecdiags_vrt
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiags_vrt, diags_vrtTstep, start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step', record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_FTYPE (nciddiags_vrt, diags_vrtTime, record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time2
!
      ierr=nf_put_var1_FTYPE (nciddiags_vrt, diags_vrtTime2, record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime2))
        write(stdout,1) vname(1,indxTime2)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif

!---------------------------------------------------------------
! Momentum diagnostic variables.
!---------------------------------------------------------------
!
!
        itrc = 1
          if (wrtdiags_vrt(itrc)) then

            ivar=p2dvar

!
! indxvrtXadv
!
           work2d=vrtXadv(:,:)
           call fillvalue2d(work2d,nciddiags_vrt,diags_vrtXadv(itrc),
     &                      indxvrtXadv+itrc-1,
     &                      record,ivar,type)  
!
! indxvrtYadv
!
           work2d=vrtYadv(:,:)
           call fillvalue2d(work2d,nciddiags_vrt,diags_vrtYadv(itrc),
     &                      indxvrtYadv+itrc-1,
     &                      record,ivar,type)
!
! indxvrtHdiff
!
           work2d=vrtHdiff(:,:)
           call fillvalue2d(work2d,nciddiags_vrt,diags_vrtHdiff(itrc),
     &                      indxvrtHdiff+itrc-1,
     &                      record,ivar,type)
!
! indxvrtCor
!
           work2d=vrtCor(:,:)
           call fillvalue2d(work2d,nciddiags_vrt,diags_vrtCor(itrc),
     &                      indxvrtCor+itrc-1,
     &                      record,ivar,type)  
!
! indxvrtPrsgrd
!
           work2d=vrtPrsgrd(:,:)
           call fillvalue2d(work2d,nciddiags_vrt,diags_vrtPrsgrd(itrc),
     &                      indxvrtPrsgrd+itrc-1,
     &                      record,ivar,type)  
!
! indxvrtHmix
!
           work2d=vrtHmix(:,:)
           call fillvalue2d(work2d,nciddiags_vrt,diags_vrtHmix(itrc),
     &                      indxvrtHmix+itrc-1,
     &                      record,ivar,type)  
!
! indxvrtVmix
!
           work2d=vrtVmix(:,:)
           call fillvalue2d(work2d,nciddiags_vrt,diags_vrtVmix(itrc),
     &                      indxvrtVmix+itrc-1,
     &                      record,ivar,type)  

!
! indxvrtrate
!
           work2d=vrtrate(:,:)
           call fillvalue2d(work2d,nciddiags_vrt,diags_vrtrate(itrc),
     &                      indxvrtrate+itrc-1,
     &                      record,ivar,type)  
!

!
! indxvrtVmix2
!
           work2d=vrtVmix2(:,:)
           call fillvalue2d(work2d,nciddiags_vrt,diags_vrtVmix2(itrc),
     &                      indxvrtVmix2+itrc-1,
     &                      record,ivar,type)  

!
! indxvrtWind
!
           work2d=vrtWind(:,:)
           call fillvalue2d(work2d,nciddiags_vrt,diags_vrtWind(itrc),
     &                      indxvrtWind+itrc-1,
     &                      record,ivar,type)


!
! indxvrtDrag
!
           work2d=vrtDrag(:,:)
           call fillvalue2d(work2d,nciddiags_vrt,diags_vrtDrag(itrc),
     &                      indxvrtDrag+itrc-1,
     &                      record,ivar,type)
# if defined DIAGNOSTICS_BARO
!
! indxvrtBaro
!
           work2d=vrtBaro(:,:)
           call fillvalue2d(work2d,nciddiags_vrt,diags_vrtBaro(itrc),
     &                      indxvrtBaro+itrc-1,
     &                      record,ivar,type)
# endif

# if defined M3FAST
!
! indxvrtfast
!
           work2d=vrtfast(:,:)
           call fillvalue2d(work2d,nciddiags_vrt,diags_vrtfast(itrc),
     &                      indxvrtfast+itrc-1,
     &                      record,ivar,type)
# endif

        endif
!
#ifdef DEBUG
      i=5
      j=5
      write(*,*) 'Write diag vrt his at : i='
     &                      ,i,' j=',j
      write(*,*) 'vrtXadv(i,j)   ',vrtXadv(i,j)
      write(*,*) 'vrtYadv(i,j)   ',vrtYadv(i,j)
      write(*,*) 'vrtCor(i,j)    ',vrtCor(i,j)
      write(*,*) 'vrtPrsgrd(i,j) ',vrtPrsgrd(i,j)
      write(*,*) 'vrtHmix(i,j)   ',vrtHmix(i,j)
      write(*,*) 'vrtVmix(i,j)   ',vrtVmix(i,j)
# if defined DIAGNOSTICS_BARO
      write(*,*) 'vrtPrsgrd(i,j) ',vrtPrsgrd(i,j)
# endif
      write(*,*) 'vrtrate(i,j)   ',vrtrate(i,j)

      write(*,*) 'SumUVhis(i,j) = ',vrtXadv(i,j)
     &                                  +  vrtYadv(i,j)
     &                                  +  vrtCor(i,j)
     &                                  +  vrtPrsgrd(i,j)
     &                                  +  vrtHmix(i,j)
     &                                  +  vrtVmix(i,j)
     &                                  -  vrtrate(i,j)
      write(*,*) '-----------------------------------------'
#endif /* DEBUG */

  1   format(/1x,'WRT_DIAGS_VRT ERROR while writing variable ''', A,
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
      ierr=nf_close (nciddiags_vrt)
      if (nrpfdiags_vrt.gt.0 .and. record.ge.nrpfdiags_vrt) nciddiags_vrt=-1
#else
      if (nrpfdiags_vrt.gt.0 .and. record.ge.nrpfdiags_vrt) then
        ierr=nf_close (nciddiags_vrt)
        nciddiags_vrt=-1
      else
        ierr=nf_sync(nciddiags_vrt)
      endif
#endif


      if (ierr .eq. nf_noerr) then
        MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')
     &            'WRT_DIAGS_VRT -- wrote',
     &            ' diag fields into time record =', record, '/',
     &             nrecdiags_vrt  MYID
      else
        MPI_master_only write(stdout,'(/1x,2A/)') 
     &             'WRT_DIAGS_VRT ERROR: Cannot ',
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
      subroutine wrt_diags_vrt_empty
      end
#endif /* (DIAGNOSTICS_VRT) */
