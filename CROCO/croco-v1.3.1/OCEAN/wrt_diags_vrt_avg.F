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
#if defined DIAGNOSTICS_VRT && defined AVERAGES
!
!---------------------------------------------------------------
! Write diagnostics fields at requested levels into diagnostics 
! netCDF file.

!---------------------------------------------------------------
!
      subroutine wrt_diags_vrt_avg

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
      call def_diags_vrt_avg (nciddiags_vrt_avg, nrecdiags_vrt_avg,
     &                       ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(diags_vrtname_avg)
!                                          !!! WARNING: Once time
! Set record within the file.              !!! stepping has been 
!     
                                           !!! started, it is assumed
      if (iic.eq.0) nrecdiags_vrt_avg=nrecdiags_vrt_avg+1  !!! that global history
      if (nrpfdiags_vrt_avg.eq.0) then            !!! record index "nrecdiags_vrt"
        record=nrecdiags_vrt_avg                  !!! is advanced by main.
      else
        record=1+mod(nrecdiags_vrt_avg-1, nrpfdiags_vrt_avg)
      endif
!
!---------------------------------------------------------------
! Write out evolving model variables:
!---------------------------------------------------------------
!
! Time step number and record numbers.
!
      type=filetype_diags_vrt_avg
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrecdiags_vrt_avg
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiags_vrt_avg, diags_vrtTstep_avg,
     &               start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step_avg', record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_FTYPE (nciddiags_vrt_avg, diags_vrtTime_avg,
     &                                    record, timediags_vrt_avg)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time2
!
      ierr=nf_put_var1_FTYPE (nciddiags_vrt_avg, diags_vrtTime2_avg,
     &                                         record, time)
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
          if (wrtdiags_vrt_avg(itrc)) then

            ivar=p2dvar

!
! indxvrtXadv
!
           work2d=vrtXadv_avg(:,:)
           call fillvalue2d(work2d,nciddiags_vrt_avg,
     &                      diags_vrtXadv_avg(itrc),
     &                      indxvrtXadv+itrc-1,
     &                      record,ivar,type)  
!
! indxvrtYadv
!
           work2d=vrtYadv_avg(:,:)
           call fillvalue2d(work2d,nciddiags_vrt_avg,
     &                      diags_vrtYadv_avg(itrc),
     &                      indxvrtYadv+itrc-1,
     &                      record,ivar,type)
!
! indxvrtHdiff
!
           work2d=vrtHdiff_avg(:,:)
           call fillvalue2d(work2d,nciddiags_vrt_avg,
     &                      diags_vrtHdiff_avg(itrc),
     &                      indxvrtHdiff+itrc-1,
     &                      record,ivar,type)
!
! indxvrtCor
!
           work2d=vrtCor_avg(:,:)
           call fillvalue2d(work2d,nciddiags_vrt_avg,
     &                      diags_vrtCor_avg(itrc),
     &                      indxvrtCor+itrc-1,
     &                      record,ivar,type)  
!
! indxvrtPrsgrd
!
           work2d=vrtPrsgrd_avg(:,:)
           call fillvalue2d(work2d,nciddiags_vrt_avg,
     &                      diags_vrtPrsgrd_avg(itrc),
     &                      indxvrtPrsgrd+itrc-1,
     &                      record,ivar,type)  
!
! indxvrtHmix
!
           work2d=vrtHmix_avg(:,:)
           call fillvalue2d(work2d,nciddiags_vrt_avg,
     &                      diags_vrtHmix_avg(itrc),
     &                      indxvrtHmix+itrc-1,
     &                      record,ivar,type)  
!
! indxvrtVmix
!
           work2d=vrtVmix_avg(:,:)
           call fillvalue2d(work2d,nciddiags_vrt_avg,
     &                      diags_vrtVmix_avg(itrc),
     &                      indxvrtVmix+itrc-1,
     &                      record,ivar,type)  

!
! indxvrtrate
!
           work2d=vrtrate_avg(:,:)
           call fillvalue2d(work2d,nciddiags_vrt_avg,
     &                      diags_vrtrate_avg(itrc),
     &                      indxvrtrate+itrc-1,
     &                      record,ivar,type)  
!

!
! indxvrtVmix2
!
           work2d=vrtVmix2_avg(:,:)
           call fillvalue2d(work2d,nciddiags_vrt_avg,
     &                      diags_vrtVmix2_avg(itrc),
     &                      indxvrtVmix2+itrc-1,
     &                      record,ivar,type)  

!
! indxvrtWind
!
           work2d=vrtWind_avg(:,:)
           call fillvalue2d(work2d,nciddiags_vrt_avg,
     &                      diags_vrtWind_avg(itrc),
     &                      indxvrtWind+itrc-1,
     &                      record,ivar,type)


!
! indxvrtDrag
!
           work2d=vrtDrag_avg(:,:)
           call fillvalue2d(work2d,nciddiags_vrt_avg,
     &                      diags_vrtDrag_avg(itrc),
     &                      indxvrtDrag+itrc-1,
     &                      record,ivar,type)
# if defined DIAGNOSTICS_BARO
!
! indxvrtBaro
!
           work2d=vrtBaro_avg(:,:)
           call fillvalue2d(work2d,nciddiags_vrt_avg,
     &                      diags_vrtBaro_avg(itrc),
     &                      indxvrtBaro+itrc-1,
     &                      record,ivar,type)
# endif

# if defined M3FAST
!
! indxvrtfast
!
           work2d=vrtfast_avg(:,:)
           call fillvalue2d(work2d,nciddiags_vrt_avg,
     &                      diags_vrtfast_avg(itrc),
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
      write(*,*) 'vrtXadv_avg(i,j)   ',vrtXadv_avg(i,j)
      write(*,*) 'vrtYadv_avg(i,j)   ',vrtYadv_avg(i,j)
      write(*,*) 'vrtCor_avg(i,j)    ',vrtCor_avg(i,j)
      write(*,*) 'vrtPrsgrd_avg(i,j) ',vrtPrsgrd_avg(i,j)
      write(*,*) 'vrtHmix_avg(i,j)   ',vrtHmix_avg(i,j)
      write(*,*) 'vrtVmix_avg(i,j)   ',vrtVmix_avg(i,j)
# if defined DIAGNOSTICS_BARO
      write(*,*) 'vrtBaro_avg(i,j) ',vrtBaro_avg(i,j)
# endif
# if defined M3FAST
      write(*,*) 'vrtfast_avg(i,j) ',vrtfast_avg(i,j)
# endif
      write(*,*) 'vrtrate_avg(i,j)   ',vrtrate_avg(i,j)

      write(*,*) 'SumUVhis_avg(i,j) = ',vrtXadv_avg(i,j)
     &                                  +  vrtYadv_avg(i,j)
     &                                  +  vrtCor_avg(i,j)
     &                                  +  vrtPrsgrd_avg(i,j)
     &                                  +  vrtHmix_avg(i,j)
     &                                  +  vrtVmix_avg(i,j)
# if defined M3FAST
     &                                  +  vrtfast_avg(i,j)
# endif
     &                                  -  vrtrate_avg(i,j)
      write(*,*) '-----------------------------------------'
#endif /* DEBUG */

  1   format(/1x,'WRT_DIAGS_VRT_AVG ERROR while writing variable ''',A,
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
      ierr=nf_close (nciddiags_vrt_avg)
      if (nrpfdiags_vrt_avg.gt.0 .and. record.ge.nrpfdiags_vrt_avg)
     &            nciddiags_vrt_avg=-1
#else
      if (nrpfdiags_vrt_avg.gt.0 .and.
     &                       record.ge.nrpfdiags_vrt_avg) then
        ierr=nf_close (nciddiags_vrt_avg)
        nciddiags_vrt_avg=-1
      else
        ierr=nf_sync(nciddiags_vrt_avg)
      endif
#endif


      if (ierr .eq. nf_noerr) then
        MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')
     &            'WRT_DIAGS_VRT_AVG -- wrote',
     &            ' diag fields into time record =', record, '/',
     &             nrecdiags_vrt_avg  MYID
      else
        MPI_master_only write(stdout,'(/1x,2A/)') 
     &             'WRT_DIAGS_VRT_AVG ERROR: Cannot ',
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
      subroutine wrt_diags_vrt_avg_empty
      end
#endif /* (DIAGNOSTICS_VRT) */
