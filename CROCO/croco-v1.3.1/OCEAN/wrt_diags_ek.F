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
#ifdef DIAGNOSTICS_EK
!
!---------------------------------------------------------------
! Write diagnostics fields at requested levels into diagnostics 
! netCDF file.

!---------------------------------------------------------------
!
      subroutine wrt_diags_ek

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
      call def_diags_ek (nciddiags_ek, nrecdiags_ek, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(diags_ekname)
!                                          !!! WARNING: Once time
! Set record within the file.              !!! stepping has been 
!     
                                           !!! started, it is assumed
      if (iic.eq.0) nrecdiags_ek=nrecdiags_ek+1  !!! that global history
      if (nrpfdiags_ek.eq.0) then            !!! record index "nrecdiags_ek"
        record=nrecdiags_ek                  !!! is advanced by main.
      else
        record=1+mod(nrecdiags_ek-1, nrpfdiags_ek)
      endif
!
!---------------------------------------------------------------
! Write out evolving model variables:
!---------------------------------------------------------------
!
! Time step number and record numbers.
!
      type=filetype_diags_ek
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrecdiags_ek
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiags_ek, diags_ekTstep,
     &                                            start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step', record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_FTYPE (nciddiags_ek, diags_ekTime, record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time2
!
      ierr=nf_put_var1_FTYPE (nciddiags_ek, diags_ekTime2, record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime2))
        write(stdout,1) vname(1,indxTime2)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif

!---------------------------------------------------------------
! Energy diagnostic variables.
!---------------------------------------------------------------
!
!
        itrc = 1
          if (wrtdiags_ek(itrc)) then

            ivar=r2dvar

!
! indxekHadv
!
           work2d=ekHadv(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekHadv(itrc),
     &                      indxekHadv+itrc-1,
     &                      record,ivar,type)  
!
! indxekHdiff
!
           work2d=ekHdiff(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekHdiff(itrc),
     &                      indxekHdiff+itrc-1,
     &                      record,ivar,type)  
!
! indxekVadv
!
           work2d=ekVadv(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekVadv(itrc),
     &                      indxekVadv+itrc-1,
     &                      record,ivar,type)  
!
! indxekCor
!
           work2d=ekCor(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekCor(itrc),
     &                      indxekCor+itrc-1,
     &                      record,ivar,type)  
!
! indxekPrsgrd
!
           work2d=ekPrsgrd(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekPrsgrd(itrc),
     &                      indxekPrsgrd+itrc-1,
     &                      record,ivar,type)  
!
! indxekHmix
!
           work2d=ekHmix(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekHmix(itrc),
     &                      indxekHmix+itrc-1,
     &                      record,ivar,type)  
!
! indxekVmix
!
           work2d=ekVmix(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekVmix(itrc),
     &                      indxekVmix+itrc-1,
     &                      record,ivar,type)  

!
! indxekrate
!
           work2d=ekrate(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekrate(itrc),
     &                      indxekrate+itrc-1,
     &                      record,ivar,type)  
!

!
! indxekvol
!
           work2d=ekvol(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekvol(itrc),
     &                      indxekvol+itrc-1,
     &                      record,ivar,type)  
!

!
! indxekVmix2
!
           work2d=ekVmix2(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekVmix2(itrc),
     &                      indxekVmix2+itrc-1,
     &                      record,ivar,type)  

!
! indxekWind
!
           work2d=ekWind(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekWind(itrc),
     &                      indxekWind+itrc-1,
     &                      record,ivar,type)


!
! indxekDrag
!
           work2d=ekDrag(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekDrag(itrc),
     &                      indxekDrag+itrc-1,
     &                      record,ivar,type)
# if defined DIAGNOSTICS_BARO
!
! indxekBaro
!
           work2d=ekBaro(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekBaro(itrc),
     &                      indxekBaro+itrc-1,
     &                      record,ivar,type)
!
# endif

# if defined M3FAST
!
! indxekfast
!
           work2d=ekfast(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekfast(itrc),
     &                      indxekfast+itrc-1,
     &                      record,ivar,type)
!
# endif

#if defined DIAGNOSTICS_EK_MLD

!---------------------------------------------------------------
! Energy diagnostic variables.
!---------------------------------------------------------------
!



!
! indxekHadv
!
           work2d=ekHadv_mld(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekHadv_mld(itrc),
     &                      indxekHadv_mld+itrc-1,
     &                      record,ivar,type)  
!
! indxekHdiff
!
           work2d=ekHdiff_mld(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekHdiff_mld(itrc),
     &                      indxekHdiff_mld+itrc-1,
     &                      record,ivar,type)  
!
! indxekVadv
!
           work2d=ekVadv_mld(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekVadv_mld(itrc),
     &                      indxekVadv_mld+itrc-1,
     &                      record,ivar,type)  
!
! indxekCor
!
           work2d=ekCor_mld(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekCor_mld(itrc),
     &                      indxekCor_mld+itrc-1,
     &                      record,ivar,type)  
!
! indxekPrsgrd
!
           work2d=ekPrsgrd_mld(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekPrsgrd_mld(itrc),
     &                      indxekPrsgrd_mld+itrc-1,
     &                      record,ivar,type)  
!
! indxekHmix
!
           work2d=ekHmix_mld(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekHmix_mld(itrc),
     &                      indxekHmix_mld+itrc-1,
     &                      record,ivar,type)  
!
! indxekVmix
!
           work2d=ekVmix_mld(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekVmix_mld(itrc),
     &                      indxekVmix_mld+itrc-1,
     &                      record,ivar,type)  

!
! indxekrate
!
           work2d=ekrate_mld(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekrate_mld(itrc),
     &                      indxekrate_mld+itrc-1,
     &                      record,ivar,type)  
!

!
! indxekvol
!
           work2d=ekvol_mld(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekvol_mld(itrc),
     &                      indxekvol_mld+itrc-1,
     &                      record,ivar,type)  
!

!
! indxekVmix2
!
           work2d=ekVmix2_mld(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekVmix2_mld(itrc),
     &                      indxekVmix2_mld+itrc-1,
     &                      record,ivar,type)  
# if defined DIAGNOSTICS_BARO
!
! indxekBaro
!
           work2d=ekBaro_mld(:,:)
           call fillvalue2d(work2d,nciddiags_ek,diags_ekBaro_mld(itrc),
     &                      indxekBaro_mld+itrc-1,
     &                      record,ivar,type)
!
# endif

#endif

        endif


#ifdef DEBUG
      i=5
      j=5
      write(*,*) 'Write diag ek his at : i='
     &                      ,i,' j=',j
      write(*,*) 'ekHadv(i,j)   ',ekHadv(i,j)
      write(*,*) 'ekHdiff(i,j)   ',ekHdiff(i,j)
      write(*,*) 'ekVadv(i,j)   ',ekVadv(i,j)
      write(*,*) 'ekCor(i,j)    ',ekCor(i,j)
      write(*,*) 'ekPrsgrd(i,j) ',ekPrsgrd(i,j)
      write(*,*) 'ekHmix(i,j)   ',ekHmix(i,j)
      write(*,*) 'ekVmix(i,j)   ',ekVmix(i,j)
      write(*,*) 'ekvol(i,j)   ',ekvol(i,j)

      write(*,*) 'ekrate(i,j)   ',ekrate(i,j)
# if defined DIAGNOSTICS_BARO
      write(*,*) 'ekBaro(i,j) ',ekBaro(i,j)
# endif
# if defined M3FAST
      write(*,*) 'ekfast(i,j) ',ekfast(i,j)
# endif
      write(*,*) 'SumUVhis(i,j) = ',ekHadv(i,j)
     &                                  +  ekVadv(i,j)
     &                                  +  ekCor(i,j)
     &                                  +  ekPrsgrd(i,j)
     &                                  +  ekHmix(i,j)
     &                                  +  ekVmix(i,j)
     &                                  +  ekvol(i,j)
# if defined M3FAST
     &                                  +  ekfast(i,j)
# endif
     &                                  -  ekrate(i,j)
      write(*,*) '-----------------------------------------'
#endif /* DEBUG */

  1   format(/1x,'WRT_DIAGS_EK ERROR while writing variable ''', A,
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
      ierr=nf_close (nciddiags_ek)
      if (nrpfdiags_ek.gt.0 .and. record.ge.nrpfdiags_ek) nciddiags_ek=-1
#else
      if (nrpfdiags_ek.gt.0 .and. record.ge.nrpfdiags_ek) then
        ierr=nf_close (nciddiags_ek)
        nciddiags_ek=-1
      else
        ierr=nf_sync(nciddiags_ek)
      endif
#endif


      if (ierr .eq. nf_noerr) then
        MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')
     &            'WRT_DIAGS_EK -- wrote',
     &            ' diag fields into time record =', record, '/',
     &             nrecdiags_ek  MYID
      else
        MPI_master_only write(stdout,'(/1x,2A/)') 
     &             'WRT_DIAGS_EK ERROR: Cannot ',
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
      subroutine wrt_diags_ek_empty
      end
#endif /* (DIAGNOSTICS_EK) */
