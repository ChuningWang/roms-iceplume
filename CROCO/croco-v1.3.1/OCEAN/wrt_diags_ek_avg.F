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
#if defined DIAGNOSTICS_EK && defined AVERAGES
!
!---------------------------------------------------------------
! Write diagnostics fields at requested levels into diagnostics 
! netCDF file.

!---------------------------------------------------------------
!
      subroutine wrt_diags_ek_avg

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
#include "diags_ek.h"
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
      call def_diags_ek_avg (nciddiags_ek_avg, nrecdiags_ek_avg, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(diags_ekname_avg)
!                                          !!! WARNING: Once time
! Set record within the file.              !!! stepping has been 
!     
                                           !!! started, it is assumed
      if (iic.eq.0) nrecdiags_ek_avg=nrecdiags_ek_avg+1  !!! that global history
      if (nrpfdiags_ek_avg.eq.0) then            !!! record index "nrecdiags_ek"
        record=nrecdiags_ek_avg                  !!! is advanced by main.
      else
        record=1+mod(nrecdiags_ek_avg-1, nrpfdiags_ek_avg)
      endif
!
!---------------------------------------------------------------
! Write out evolving model variables:
!---------------------------------------------------------------
!
! Time step number and record numbers.
!
      type=filetype_diags_ek_avg
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrecdiags_ek_avg
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiags_ek_avg, diags_ekTstep_avg,
     &                       start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step_avg', record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_FTYPE (nciddiags_ek_avg, diags_ekTime_avg,
     &                       record, timediags_ek_avg)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time2
!
      ierr=nf_put_var1_FTYPE (nciddiags_ek_avg, diags_ekTime2_avg,
     &                       record, time)
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
          if (wrtdiags_ek_avg(itrc)) then

            ivar=r2dvar

!
! indxekHadv
!
           work2d=ekHadv_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekHadv_avg(itrc),
     &                      indxekHadv+itrc-1,
     &                      record,ivar,type)  
!
! indxekHdiff
!
           work2d=ekHdiff_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekHdiff_avg(itrc),
     &                      indxekHdiff+itrc-1,
     &                      record,ivar,type)  
!
! indxekVadv
!
           work2d=ekVadv_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekVadv_avg(itrc),
     &                      indxekVadv+itrc-1,
     &                      record,ivar,type)  
!
! indxekCor
!
           work2d=ekCor_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekCor_avg(itrc),
     &                      indxekCor+itrc-1,
     &                      record,ivar,type)  
!
! indxekPrsgrd
!
           work2d=ekPrsgrd_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekPrsgrd_avg(itrc),
     &                      indxekPrsgrd+itrc-1,
     &                      record,ivar,type)  
!
! indxekHmix
!
           work2d=ekHmix_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekHmix_avg(itrc),
     &                      indxekHmix+itrc-1,
     &                      record,ivar,type)  
!
! indxekVmix
!
           work2d=ekVmix_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekVmix_avg(itrc),
     &                      indxekVmix+itrc-1,
     &                      record,ivar,type)  

!
! indxekrate
!
           work2d=ekrate_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekrate_avg(itrc),
     &                      indxekrate+itrc-1,
     &                      record,ivar,type)  
!

!
! indxekvol
!
           work2d=ekvol_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekvol_avg(itrc),
     &                      indxekvol+itrc-1,
     &                      record,ivar,type)  
!

!
! indxekVmix2
!
           work2d=ekVmix2_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekVmix2_avg(itrc),
     &                      indxekVmix2+itrc-1,
     &                      record,ivar,type)  

!
! indxekWind
!
           work2d=ekWind_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekWind_avg(itrc),
     &                      indxekWind+itrc-1,
     &                      record,ivar,type)


!
! indxekDrag
!
           work2d=ekDrag_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekDrag_avg(itrc),
     &                      indxekDrag+itrc-1,
     &                      record,ivar,type)

# if defined DIAGNOSTICS_BARO
!
! indxekBaro
!
           work2d=ekBaro_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekBaro_avg(itrc),
     &                      indxekBaro+itrc-1,
     &                      record,ivar,type)
# endif

# if defined M3FAST
!
! indxekfast
!
           work2d=ekfast_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekfast_avg(itrc),
     &                      indxekfast+itrc-1,
     &                      record,ivar,type)
# endif


# if defined DIAGNOSTICS_EK_MLD



!
! indxekHadv
!
           work2d=ekHadv_mld_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekHadv_mld_avg(itrc),
     &                      indxekHadv_mld+itrc-1,
     &                      record,ivar,type)  
!
! indxekHdiff
!
           work2d=ekHdiff_mld_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekHdiff_mld_avg(itrc),
     &                      indxekHdiff_mld+itrc-1,
     &                      record,ivar,type)  
!
! indxekVadv
!
           work2d=ekVadv_mld_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekVadv_mld_avg(itrc),
     &                      indxekVadv_mld+itrc-1,
     &                      record,ivar,type)  
!
! indxekCor
!
           work2d=ekCor_mld_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekCor_mld_avg(itrc),
     &                      indxekCor_mld+itrc-1,
     &                      record,ivar,type)  
!
! indxekPrsgrd
!
           work2d=ekPrsgrd_mld_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekPrsgrd_mld_avg(itrc),
     &                      indxekPrsgrd_mld+itrc-1,
     &                      record,ivar,type)  
!
! indxekHmix
!
           work2d=ekHmix_mld_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekHmix_mld_avg(itrc),
     &                      indxekHmix_mld+itrc-1,
     &                      record,ivar,type)  
!
! indxekVmix
!
           work2d=ekVmix_mld_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekVmix_mld_avg(itrc),
     &                      indxekVmix_mld+itrc-1,
     &                      record,ivar,type)  

!
! indxekrate
!
           work2d=ekrate_mld_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekrate_mld_avg(itrc),
     &                      indxekrate_mld+itrc-1,
     &                      record,ivar,type)  
!

!
! indxekvol
!
           work2d=ekvol_mld_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekvol_mld_avg(itrc),
     &                      indxekvol_mld+itrc-1,
     &                      record,ivar,type)  
!

!
! indxekVmix2
!
           work2d=ekVmix2_mld_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekVmix2_mld_avg(itrc),
     &                      indxekVmix2_mld+itrc-1,
     &                      record,ivar,type)
# if defined DIAGNOSTICS_BARO
!
! indxekBaro
!
           work2d=ekBaro_mld_avg(:,:)
           call fillvalue2d(work2d,nciddiags_ek_avg,
     &                      diags_ekBaro_mld_avg(itrc),
     &                      indxekBaro_mld+itrc-1,
     &                      record,ivar,type)
# endif

# endif

        endif
!
#ifdef DEBUG
      i=5
      j=5
      write(*,*) 'Write diag ek his at : i='
     &                      ,i,' j=',j
      write(*,*) 'ekHadv_avg(i,j)   ',ekHadv_avg(i,j)
      write(*,*) 'ekHdiff_avg(i,j)   ',ekHdiff_avg(i,j)
      write(*,*) 'ekVadv_avg(i,j)   ',ekVadv_avg(i,j)
      write(*,*) 'ekCor_avg(i,j)    ',ekCor_avg(i,j)
      write(*,*) 'ekPrsgrd_avg(i,j) ',ekPrsgrd_avg(i,j)
      write(*,*) 'ekHmix_avg(i,j)   ',ekHmix_avg(i,j)
      write(*,*) 'ekVmix_avg(i,j)   ',ekVmix_avg(i,j)
      write(*,*) 'ekvol_avg(i,j)   ',ekvol_avg(i,j)
# if defined DIAGNOSTICS_BARO
      write(*,*) 'ekBaro_avg(i,j) ',ekBaro_avg(i,j)
# endif
# if defined M3FAST
      write(*,*) 'ekfast_avg(i,j) ',ekfast_avg(i,j)
# endif
      write(*,*) 'ekrate_avg(i,j)   ',ekrate_avg(i,j)

      write(*,*) 'SumUVhis_avg(i,j) = ',ekHadv_avg(i,j)
     &                                  +  ekVadv_avg(i,j)
     &                                  +  ekCor_avg(i,j)
     &                                  +  ekPrsgrd_avg(i,j)
     &                                  +  ekHmix_avg(i,j)
     &                                  +  ekVmix_avg(i,j)
     &                                  +  ekvol_avg(i,j)
     &                                  -  ekrate_avg(i,j)
      write(*,*) '-----------------------------------------'
#endif /* DEBUG */

  1   format(/1x,'WRT_DIAGS_EK_AVG ERROR while writing variable ''', A,
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
      ierr=nf_close (nciddiags_ek_avg)
      if (nrpfdiags_ek_avg.gt.0 .and. record.ge.nrpfdiags_ek_avg)
     &            nciddiags_ek_avg=-1
#else
      if (nrpfdiags_ek_avg.gt.0 .and. record.ge.nrpfdiags_ek_avg) then
        ierr=nf_close (nciddiags_ek_avg)
        nciddiags_ek_avg=-1
      else
        ierr=nf_sync(nciddiags_ek_avg)
      endif
#endif


      if (ierr .eq. nf_noerr) then
        MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')
     &            'WRT_DIAGS_EK_AVG -- wrote',
     &            ' diag fields into time record =', record, '/',
     &             nrecdiags_ek_avg  MYID
      else
        MPI_master_only write(stdout,'(/1x,2A/)') 
     &             'WRT_DIAGS_EK_AVG ERROR: Cannot ',
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
      subroutine wrt_diags_ek_avg_empty
      end
#endif /* (DIAGNOSTICS_EK) */
