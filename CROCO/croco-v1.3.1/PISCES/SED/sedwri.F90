! $Id: wrt_rst.F 1571 2014-07-01 12:38:05Z gcambon $
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
#include "cppdefs.h"

MODULE sedwri

#if defined key_pisces

   !! * Modules used
   USE sed
   USE sedarr
   USE sms_pisces, ONLY : rtrn, rfact
#ifdef AGRIF
      USE param, ONLY : Lmmpi,Mmmpi
#endif

   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC sed_wri         ! routine called by opa.F90

   !!* Substitution
#  include "ocean2pisces.h90"

#if ! defined XIOS

# define rec_per_file nrpfsedpis_avg
# define vidTime sedTimepis_avg
# define vidTime2 sedTime2pis_avg
# define vidTstep sedTsteppis_avg


CONTAINS         ! Write model prognostic

      SUBROUTINE sed_wri      ! variables into restart
                                  ! netCDF file.
# include "netcdf.inc"

      INTEGER :: ierr, record, lstr, lvar, lenstr   &
      &  , start(2), count(2), ibuff(4), nf_fwrite, itrc, itype
      INTEGER :: ji, jj, jk, jn
      REAL(wp) :: zrate
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: trcsedtmp, trcsedi, flxsedi3d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: flxsedtmp, flxsedi2d
      REAL(wp), DIMENSION(jpoce,jpksed)   :: zdta
      REAL(wp), DIMENSION(jpoce,jpwatp1)  :: zflx

#ifdef MPI
#define LOCALLM Lmmpi
#define LOCALMM Mmmpi
#else
#define LOCALLM Lm
#define LOCALMM Mm
#endif


#if defined MPI & !defined PARALLEL_FILES
      INCLUDE 'mpif.h'
      INTEGER status(MPI_STATUS_SIZE), blank
#endif

#if defined MPI & !defined PARALLEL_FILES
      IF (mynode > 0) THEN
         call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,      &
         &                 1, MPI_COMM_WORLD, status, ierr)
      ENDIF
#endif
!
! Create/open restart file; write grid arrays, if so needed.
!
      CALL def_wri_sed (ncidwrised, nrecsedpis_avg, ierr)
      IF (ierr .NE. nf_noerr) GOTO 99
!                                            !!! WARNING: Here it is
! Set record within the file.                !!! assumed that global
!                                            !!! restart record index 
      nrecsedpis_avg = max(nrecsedpis_avg,1)                 !!! nrecrst is already
      IF (nrpfsedpis_avg == 0) THEN                 !!! advanced by main.
         record = nrecsedpis_avg
      ELSE
         record = 1+mod(nrecsedpis_avg-1, nrpfsedpis_avg)
      ENDIF

!
! Write out evolving model variables:
! ----- --- -------- ----- ----------
!
! Time step number and record indices. 
!
      itype=filetype_avg

      ibuff(1) = iic
      ibuff(2) = nrecrst
      ibuff(3) = nrechis
      ibuff(4) = nrecsedpis_avg
      start(1) = 1
      start(2) = record
      count(1) = 4
      count(2) = 1
      ierr = nf_put_vara_int (ncidwrised, sedTsteppis_avg, start, count, ibuff)
      IF (ierr .NE. nf_noerr) THEN
         WRITE(stdout,1) 'time_step', record, ierr      
         GOTO 99                                           !--> ERROR
      ENDIF
!!
! Averaged Time.
!
      ierr = nf_put_var1_FTYPE (ncidwrised, sedTimepis_avg, record, timesedpis_avg)
      IF (ierr .NE. nf_noerr) THEN
         lvar = lenstr(vname(1,indxTime))
         WRITE(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
         GOTO 99                                           !--> ERROR
      ENDIF
!     
!     Averaged time2
!     
      ierr=nf_put_var1_FTYPE (ncidwrised, sedTime2pis_avg, record, timesedpis_avg)
      if (ierr .ne. nf_noerr) then
         lvar=lenstr(vname(1,indxTime2))
         write(stdout,1) vname(1,indxTime2)(1:lvar), record, ierr
         goto 99                !--> ERROR
      endif
!
! Tracer variables.
!
!
      ALLOCATE(trcsedtmp(GLOBAL_2D_ARRAY,jpksed,jptrased), trcsedi(PRIV_2D_BIOARRAY,jpksed,jptrased) )
      trcsedi(:,:,:,:)   = 0.0
      trcsedtmp(:,:,:,:) = 0.0

      ! Back to 2D geometry
      DO jn = 1, jpsol
         CALL unpack_arr( jpoce, trcsedi(PRIV_2D_BIOARRAY,1:jpksed,jn) , iarroce(1:jpoce), &
         &                       solcp_avg(1:jpoce,1:jpksed,jn ) )
      END DO

      DO jn = 1, jpwat
         CALL unpack_arr( jpoce, trcsedi(PRIV_2D_BIOARRAY,1:jpksed,jpsol+jn) , iarroce(1:jpoce), &
         &                       pwcp_avg(1:jpoce,1:jpksed,jn  )  )
      END DO
!
      DO itrc = 1, jptrased
         DO jk = 1, jpksed
            DO jj = 1, LOCALMM
               DO ji = 1, LOCALLM
                  trcsedtmp(ji,jj,jk,itrc) = trcsedi(ji,jj,jk,itrc)
               END DO
            END DO
         END DO
      END DO

      DO itrc = 1, jptrased
         ierr = nf_fwrite(trcsedtmp(START_2D_ARRAY,1,itrc), ncidwrised,   &
         &                             wrised(itrc), record, r3dsed)
        IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) TRIM(sedtrcd(itrc)), record, ierr
          GOTO 99                                         !--> ERROR
        ENDIF
        END DO

      DEALLOCATE(trcsedtmp, trcsedi )

!
!     3D diagnostics variables
!
      ALLOCATE(trcsedtmp(GLOBAL_2D_ARRAY,jpksed,jpdia3dsed), trcsedi(PRIV_2D_BIOARRAY,jpksed,jpdia3dsed) )
! pH
      zdta(:,:) = 0.
      DO jk = 1, jpksed
         DO ji = 1, jpoce
            zdta(ji,jk) = -LOG10( hipor(ji,jk) / ( densSW(ji) + rtrn ) + rtrn )
         ENDDO
      ENDDO

      CALL unpack_arr( jpoce, trcsedi(PRIV_2D_BIOARRAY,1:jpksed,1)  , iarroce(1:jpoce), &
         &                   zdta(1:jpoce,1:jpksed)  )

      CALL unpack_arr( jpoce, trcsedi(PRIV_2D_BIOARRAY,1:jpksed,2)  , iarroce(1:jpoce), &
         &                   co3por(1:jpoce,1:jpksed)  )

      CALL unpack_arr( jpoce, trcsedi(PRIV_2D_BIOARRAY,1:jpksed,3)  , iarroce(1:jpoce), &
         &                   sedligand(1:jpoce,1:jpksed)  )

      DO itrc = 1, jpdia3dsed
         DO jk = 1, jpksed
            DO jj = 1, LOCALMM
               DO ji = 1, LOCALLM
                  trcsedtmp(ji,jj,jk,itrc) = trcsedi(ji,jj,jk,itrc)
               END DO
            END DO
         END DO
      END DO

      DO itrc = 1, jpdia3dsed
         ierr = nf_fwrite(trcsedtmp(START_2D_ARRAY,1,itrc), ncidwrised,   &
         &                             dia3wrised(itrc), record, r3dsed)
         IF (ierr .NE. nf_noerr) THEN
            WRITE(stdout,1) TRIM(seddia3d(itrc)), record, ierr
            GOTO 99                                         !--> ERROR
         ENDIF
      END DO

      DEALLOCATE(trcsedtmp, trcsedi)

!
!     3D diagnostics variables
!
      ALLOCATE(flxsedtmp(GLOBAL_2D_ARRAY,jpdia2dsed), flxsedi2d(PRIV_2D_BIOARRAY,jpdia2dsed) )

      zflx(:,:) = 0.
      flxsedtmp(:,:,:) = 0.0
      flxsedi2d(:,:,:)  = 0.0
      ! Calculation of fluxes mol/cm2/s
      DO jn = 1, jpwat
         DO ji = 1, jpoce
            zflx(ji,jn) = ( pwcp(ji,1,jn) - pwcp_dta(ji,jn) ) &
               &         * 1.e3 / 1.e2 * dzkbot(ji) / rfact
         ENDDO
      ENDDO

      ! Calculation of accumulation rate per dt
      DO jn = 1, jpsol
         zrate =  1.0 / ( denssol * por1(jpksed) ) / rfact
         DO ji = 1, jpoce
            zflx(ji,jpwatp1) = zflx(ji,jpwatp1) + ( tosed(ji,jn) - fromsed(ji,jn) ) * zrate
         ENDDO
      ENDDO

      DO itrc = 1, jpdia2dsed - 1
         CALL unpack_arr( jpoce, flxsedi2d(PRIV_2D_BIOARRAY,itrc), iarroce(1:jpoce), zflx(1:jpoce,itrc)  )
      END DO
      zflx(:,1) = dzdep(:) / dtsed
      CALL unpack_arr( jpoce, flxsedi2d(PRIV_2D_BIOARRAY,jpdia2dsed), iarroce(1:jpoce), zflx(1:jpoce,1) )

      DO itrc = 1, jpdia2dsed
         DO jj = 1, LOCALMM
            DO ji = 1, LOCALLM
               flxsedtmp(ji,jj,itrc) = flxsedi2d(ji,jj,itrc)
            END DO
         END DO
      END DO

      DO itrc = 1, jpdia2dsed
         ierr = nf_fwrite(flxsedtmp(START_2D_ARRAY,itrc), ncidwrised,   &
         &                             dia2wrised(itrc), record, r2dsed)
         IF (ierr .NE. nf_noerr) THEN
            WRITE(stdout,1) TRIM(seddia2d(itrc)), record, ierr
            GOTO 99                                         !--> ERROR
         ENDIF
      END DO

      DEALLOCATE(flxsedtmp, flxsedi2d)

  1   FORMAT(/1x, 'SED_WRI ERROR while writing variable ''', A,   &
       &           ''' into output file.', /11x, 'Time record:', &
       &               i6, 3x, 'netCDF error code', i4, 3x, A,i4) 
      GOTO 100 
  99  may_day_flag=3
 100  CONTINUE

!
! Synchronize restart netCDF file to disk to allow other
! processes to access data immediately after it is written.
!
#if defined MPI & !defined PARALLEL_FILES
      ierr = nf_close (ncidwrised)
      IF (nrpfsedpis_avg > 0 .AND. record >= nrpfsedpis_avg) ncidwrised = -1
#else
      IF (nrpfsedpis_avg > 0 .AND. record >= nrpfsedpis_avg) THEN
        ierr = nf_close (ncidwrised)
        ncidwrised = -1
      ELSE
        ierr = nf_sync(ncidwrised)
      ENDIF
#endif
      IF (ierr == nf_noerr) THEN
         MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')    & 
         &            'SED_WRI -- wrote ',                          &
         &            'output fields into time record =', record, '/',  &
         &             nrecsedpis_avg
      ELSE
         MPI_master_only  write(stdout,'(/1x,2A/)')     & 
         &             'SED_WRI ERROR: Cannot ',        &
         &             'synchronize/close output netCDF file.'
         may_day_flag = 3
      ENDIF

#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      IF (mynode < NNODES-1) THEN
         CALL MPI_Send (blank, 1, MPI_INTEGER, mynode+1, 1, MPI_COMM_WORLD,  ierr)
      ENDIF
#endif
      RETURN
      END SUBROUTINE sed_wri

      SUBROUTINE def_wri_sed( ncid, total_rec, ierr)  ! restart netCDF

# include "netcdf.inc"

      logical :: create_new_file
      integer :: ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim    &
      &      , r2dgrd(3),  auxil(2),  checkdims                           &
      &      , r3dgrd(4),  u3dgrd(4), v3dgrd(4),  w3dgrd(4), itrc
      CHARACTER(len=20) :: cltra
      character(len=60) :: text

!
! Put time record index into file name. In  the case when model
! output is to be arranged into sequence of named files, the naming
! convention is as follows: 'rst_root.INDEX.[MPI_node.]nc', where
! INDEX is an integer number such that (i) it is divisible by the
! specified number of records per file; and (ii)
!
!      INDEX + record_within_the_file = total_record
!
! where, 1 =< record_within_the_file =< records_per_file, so that
! total_record changes continuously throughout the sequence of files.
!
      if (may_day_flag.ne.0) return      !-->  EXIT

      ierr = 0
      lstr = lenstr(TRIM(cn_sedwri_out))
      IF (rec_per_file > 0) THEN
         lvar=total_rec - (1+mod(total_rec-1, rec_per_file))
         call insert_time_index (TRIM(cn_sedwri_out), lstr, lvar, ierr)
         IF (ierr .NE. 0) GOTO 99
      ENDIF
!
! Decide whether to create a new file, or open existing one.
! Overall the whole code below is organized into 3-way switch,
!
! 10  if (create_new_file) then
!        .... create new file, save netCDF ids for all variables;
!     elseif (ncid.eq.-1) then
!        .... try to open existing file and check its dimensions
!       if (cannot be opened or rejected) then
!         create_new_file=.true.
!         goto 10
!       endif   and prepare
!        .... prepare the file for adding new data,
!        .... find and save netCDF ids for all variables
!     else
!        .... just open, no checking, all ids are assumed to be
!        .... already known (MPI single file output only).
!     endif
!
! which is designed to implement flexible opening policy:
! if ldefhis=.true., it forces creation of a new file [if the
! file already exists, it will be overwritten]; on the other hand,
! ldefhis=.false., it is assumed that the file already exists and
! an attempt to open it is made; if the attempt is successful, the
! file is prepared for appending hew data; if it fails, a new file
! is created.
!
      create_new_file = ldefsedpis_avg
      IF (ncid .NE. -1) create_new_file = .false.
#if defined MPI & !defined PARALLEL_FILES
      IF (mynode > 0) create_new_file = .false.
#endif
!
! Create new restart file:    Put global attributes
!======= === ======= =====    and define all variables.
!
  10  if (create_new_file) then

        ierr  = nf_create(TRIM(cn_sedwri_out),NF_CLOBBER, ncid)
        IF (ierr .NE. nf_noerr) THEN
           WRITE(stdout,'(/3(1x,A)/)') 'ERROR in DEF_WRI_SED: Cannot',    &
           &             'create output NetCDF file:', TRIM(cn_sedwri_out)
           may_day_flag = 3
           RETURN
        ENDIF
!
! Put global attributes.
! --- ------ -----------
!
        CALL put_global_atts (ncid, ierr)
        IF (ierr .NE. nf_noerr) THEN
           WRITE(stdout,*) TRIM(cn_sedwri_out)
           may_day_flag = 3
           RETURN                         !-->  EXIT
        ENDIF
!
! Define dimensions of staggered fields.
! ------ ---------- -- --------- -------
!
        ierr = nf_def_dim (ncid, 'xi_rho',   xi_rho,   r2dgrd(1))
        ierr = nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr = nf_def_dim (ncid, 'profsed',   jpksed,  r3dgrd(3))
        ierr = nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr = nf_def_dim (ncid, 'auxil',    4,        auxil(1))
        auxil(2)  = timedim

        r2dgrd(3) = timedim           ! Free surface
        r3dgrd(1) = r2dgrd(1)         !
        r3dgrd(2) = r2dgrd(2)         ! 3D RHO-type
        r3dgrd(4) = timedim           !
!
! Define evolving model variables:
! ------ -------- ----- ----------
!
!
! Time step number and time record numbers:
!
        ierr = nf_def_var (ncid, 'time_step', nf_int, 2, auxil,     &
        &       vidTstep)
        ierr = nf_put_att_text (ncid, sedTsteppis_avg, 'long_name', 48,    &
        &       'time step and record numbers from initialization')
!
! Time.
!
        lvar = lenstr(vname(1,indxTime))
        ierr = nf_def_var (ncid, vname(1,indxTime)(1:lvar),           &
        &                              NF_FOUT, 1, timedim, vidTime)
        text='avg'//vname(2,indxTime)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidTime, 'long_name',             &
        &       lvar, text(1:lvar))
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', lvar,       &
        &                                vname(2,indxTime)(1:lvar))

        lvar = lenstr(vname(3,indxTime))
        ierr = nf_put_att_text (ncid, vidTime, 'units',     lvar,     &
        &                                  vname(3,indxTime)(1:lvar))
        lvar = lenstr (vname(4,indxTime))
        ierr = nf_put_att_text(ncid, vidTime, 'field',     lvar,      &
        &                                  vname(4,indxTime)(1:lvar))

        CALL nf_add_attribute(ncid, vidTime, indxTime, 5, NF_FOUT, ierr)

! Time2.
!
        lvar = lenstr(vname(1,indxTime2))
        ierr = nf_def_var (ncid, vname(1,indxTime2)(1:lvar),            &
        &                              NF_FOUT, 1, timedim, vidTime2)
        text='avg'//vname(2,indxTime2)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidTime2, 'long_name',              &
        &       lvar, text(1:lvar))
        lvar=lenstr(vname(2,indxTime2))
        ierr=nf_put_att_text (ncid, vidTime2, 'long_name', lvar,        &
        &                                vname(2,indxTime2)(1:lvar))

        lvar = lenstr(vname(3,indxTime2))
        ierr = nf_put_att_text (ncid, vidTime2, 'units',     lvar,     &
        &                                  vname(3,indxTime2)(1:lvar))
        lvar = lenstr (vname(4,indxTime2))
        ierr = nf_put_att_text(ncid, vidTime2, 'field',     lvar,      &
        &                                  vname(4,indxTime2)(1:lvar))

        CALL nf_add_attribute(ncid, vidTime2, indxTime2, 5, NF_FOUT, ierr)
!
! Tracer variables.
!
        DO itrc = 1, jptrased
           cltra = TRIM(sedtrcd(itrc))
           ierr  = nf_def_var (ncid, cltra, NF_FOUT, 4, r3dgrd, wrised(itrc))
           text='averaged '//TRIM(sedtrcl(itrc))
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, wrised(itrc), 'long_name', lvar, text(1:lvar))
           lvar = lenstr(TRIM(sedtrcu(itrc)))
           ierr = nf_put_att_text (ncid, wrised(itrc), 'units', lvar, TRIM(sedtrcu(itrc)))
           CALL nf_add_attribute(ncid, wrised(itrc), itrc, 5, NF_FOUT,ierr)
        END DO
!
! 3D diagnostics variable
!
        DO itrc = 1, jpdia3dsed
           cltra = TRIM(seddia3d(itrc))
           ierr  = nf_def_var (ncid, cltra, NF_FOUT, 4, r3dgrd, dia3wrised(itrc))
           text='averaged '//TRIM(seddia3l(itrc))
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, dia3wrised(itrc), 'long_name', lvar, text(1:lvar))
           lvar = lenstr(TRIM(seddia3u(itrc)))
           ierr = nf_put_att_text (ncid, dia3wrised(itrc), 'units', lvar, TRIM(seddia3u(itrc)))
           CALL nf_add_attribute(ncid, dia3wrised(itrc), itrc, 5, NF_FOUT,ierr)
        END DO
!
! 2D diagnostics variable
!
        DO itrc = 1, jpdia2dsed
           cltra = TRIM(seddia2d(itrc))
           ierr  = nf_def_var (ncid, cltra, NF_FOUT, 3, r2dgrd, dia2wrised(itrc))
           text='averaged '//TRIM(seddia2l(itrc))
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, dia2wrised(itrc), 'long_name', lvar, text(1:lvar))
           lvar = lenstr(TRIM(seddia2u(itrc)))
           ierr = nf_put_att_text (ncid, dia2wrised(itrc), 'units', lvar, TRIM(seddia2u(itrc)))
           CALL nf_add_attribute(ncid, dia2wrised(itrc), itrc, 5, NF_FOUT,ierr)
        END DO
!
! Leave definition mode.                  Also initialize record
! ----- ---------- -----                  dimension size to zero.
!
        ierr = nf_enddef(ncid)
        WRITE(*,'(6x,4A,1x,A,i4)') 'DEF_WRI_SED - Created new ',        &
        &              'netCDF file ''', TRIM(cn_sedwri_out), '''.'
!
! Open an existing file and prepare for appending data.
! ==== == ======== ==== === ======= === ========= =====
! Check consistency of the dimensions of fields from the
! file with model dimensions. Determine the current size
! of unlimited dimension and set initial record [in the
! case of MPI serialized output, at this moment the last
! time record is assumed to be **partially** written by
! MPI processes with lower rank. Thus the next write is
! expected to be into the same record rather than next
! one (except MPI-master, who initializes the record).
!
! In the case when file is rejected (whether it cannot
! be opened, or something is wrong with its dimensions, 
! create new file. 
!
      ELSEIF (ncid == -1) THEN
        ierr = nf_open (TRIM(cn_sedwri_out), nf_write, ncid)
        IF (ierr == nf_noerr) THEN
           ierr = checkdims (ncid, TRIM(cn_sedwri_out), lstr, rec)
           IF (ierr == nf_noerr) THEN
              IF (rec_per_file == 0) THEN
                 ierr = rec+1 - total_rec
              ELSE
                 ierr = rec+1 - (1+mod(total_rec-1, rec_per_file))
              ENDIF
              IF (ierr > 0) THEN
                 MPI_master_only write( stdout,                              &
        &                 '(/1x,A,I5,1x,A/8x,3A,I5/8x,A,I5,1x,A/)'         &
        &           ) 'DEF_WRI_SED WARNING: Actual number of records', rec,    &
        &             'in netCDF file',  '''',  TRIM(cn_sedwri_out),       &
        &             ''' exceeds the record number from output data',    &
        &             rec+1-ierr,'/', total_rec,', restart is assumed.'
                 rec = rec-ierr
              ELSEIF (rec_per_file == 0) THEN
                 total_rec = rec+1           ! <-- set to the next record
#if defined MPI & !defined PARALLEL_FILES
                 IF (mynode > 0) total_rec = total_rec-1
#endif
              ENDIF
              ierr = nf_noerr
           ENDIF
        ENDIF

        IF (ierr .NE. nf_noerr) THEN
#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
           IF (mynode == 0) THEN
              create_new_file = .true.
              GOTO 10
           ELSE
              WRITE(stdout,'(/1x,4A, 1x,A,I4/)')     'DEF_WRI_SED ERROR: ',    &
              &     'Cannot open output netCDF file ''',TRIM(cn_sedwri_out),'''.'
              GOTO 99                                     !--> ERROR 
           ENDIF
#else
           create_new_file=.true.
           GOTO 10
#endif
        ENDIF
!
! Find netCDF IDs of evolving model variables:
! ---- ------ --- -- -------- ----- ----------
!
! Time step indices:
!
        ierr = nf_inq_varid (ncid, 'time_step', vidTstep)
        IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) 'time_step', TRIM(cn_sedwri_out)
          GOTO 99                                         !--> ERROR
        ENDIF
!
! Time.
!
        lvar = lenstr(vname(1,indxTime))
        ierr = nf_inq_varid (ncid, vname(1,indxTime)(1:lvar), vidTime)
        IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) vname(1,indxTime)(1:lvar), TRIM(cn_sedwri_out)
          GOTO 99                                         !--> ERROR
        ENDIF
!
! Time2.
!
        lvar = lenstr(vname(1,indxTime2))
        ierr = nf_inq_varid (ncid, vname(1,indxTime2)(1:lvar), vidTime2)
        IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) vname(1,indxTime2)(1:lvar), TRIM(cn_sedwri_out)
          GOTO 99                                         !--> ERROR
        ENDIF
!
! Tracer variables.
!
       DO itrc = 1, jptrased
          ierr = nf_inq_varid (ncid, TRIM(sedtrcd(itrc)), wrised(itrc))
          IF (ierr .NE. nf_noerr) THEN
             WRITE(stdout,1) TRIM(sedtrcd(itrc)), TRIM(cn_sedwri_out)
             GOTO 99                                       !--> ERROR
          ENDIF
       END DO
!
!  3D diagnostic variables
!
       DO itrc = 1, jpdia3dsed
          ierr = nf_inq_varid (ncid, TRIM(seddia3d(itrc)), dia3wrised(itrc))
          IF (ierr .NE. nf_noerr) THEN
             WRITE(stdout,1) TRIM(seddia3d(itrc)), TRIM(cn_sedwri_out)
             GOTO 99                                       !--> ERROR
          ENDIF
       END DO
!
!  2D diagnostic variables
!
       DO itrc = 1, jpdia2dsed
          ierr = nf_inq_varid (ncid, TRIM(seddia2d(itrc)), dia2wrised(itrc))
          IF (ierr .NE. nf_noerr) THEN
             WRITE(stdout,1) TRIM(seddia2d(itrc)), TRIM(cn_sedwri_out)
             GOTO 99                                       !--> ERROR
          ENDIF
       END DO
!
        MPI_master_only WRITE(*,'(6x,2A,i4,1x,A,i4)')              &
        &             'DEF_WRI_SED -- Opened ',                        &
        &             'existing output file,  record =', rec

#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      ELSE
         ierr = nf_open (TRIM(cn_sedwri_out), nf_write, ncid)
         IF (ierr .NE. nf_noerr) THEN
            MPI_master_only WRITE(stdout,'(/1x,4A, 1x,A,I4/)')        &
            &          'DEF_WRI_SED ERROR: Cannot',                       &
            &          'open output netCDF file ''', cn_sedwri_out(1:lstr), '''.'
            GOTO 99                                         !--> ERROR
         ENDIF
#endif
      ENDIF              !<-- create_new_file
      ierr=nf_set_fill (ncid, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
        write(*,'(6x,2A,i4,1x,A,i4)') 'DEF_WRI_SED ERROR: Cannot ',     &
        &    'switch to ''nf_nofill'' more; netCDF error code =', ierr
      endif

   1  FORMAT(/1x,'DEF_WRI_SED ERROR: Cannot find variable ''',        &
      &               A, ''' in netCDF file ''', A, '''.'/)

  99  RETURN                                              !--> ERROR
      END SUBROUTINE def_wri_sed

#else

CONTAINS         ! Write model prognostic

   SUBROUTINE sed_wri
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_wri  ***
      !!
      !! ** Purpose :  output of sediment passive tracer
      !!
      !!   History :
      !!        !  06-07  (C. Ethe)  original
      !!----------------------------------------------------------------------

      INTEGER  :: ji, jj, jk, js, jw, jn
      INTEGER  :: it
      CHARACTER(len = 20)  ::  cltra 
      REAL(wp)  :: zrate
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zdta, zflx
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: ztrcsedi, flxsedi3d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: trcsedi
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: flxsedi2d

      !!-------------------------------------------------------------------


      ! Initialisation
      ! ----------------- 
      ALLOCATE( zdta(jpoce,jpksed) )    ;   ALLOCATE( zflx(jpoce,jpwatp1) )
      ALLOCATE( ztrcsedi(PRIV_2D_BIOARRAY,jpksed,jptrased) )
      ALLOCATE( flxsedi2d(PRIV_2D_BIOARRAY,jpdia2dsed) )
      ALLOCATE( flxsedi3d(PRIV_2D_BIOARRAY,jpksed,jpdia3dsed) )
      ALLOCATE( trcsedi(GLOBAL_2D_ARRAY,jpksed,jptrased) )

      ! Initialize variables
      ! --------------------

      trcsedi(:,:,:,:)   = 0.0
      ztrcsedi(:,:,:,:)  = 0.0
      flxsedi3d(:,:,:,:) = 0.0
      flxsedi2d(:,:,:)   = 0.0

      ! 2.  Back to 2D geometry
      ! -----------------------------------------------------------------
      DO jn = 1, jpsol
         CALL unpack_arr( jpoce, ztrcsedi(PRIV_2D_BIOARRAY,1:jpksed,jn) , iarroce(1:jpoce), &
         &                       solcp(1:jpoce,1:jpksed,jn ) )
      END DO
      
      DO jn = 1, jpwat
         CALL unpack_arr( jpoce, ztrcsedi(PRIV_2D_BIOARRAY,1:jpksed,jpsol + jn) , iarroce(1:jpoce), &
         &                       pwcp(1:jpoce,1:jpksed,jn  )  )
      END DO      

      ! porosity
      zdta(:,:) = 0.
      DO jk = 1, jpksed
         DO ji = 1, jpoce
            zdta(ji,jk) = -LOG10( hipor(ji,jk) / ( densSW(ji) + rtrn ) + rtrn )
         ENDDO
      ENDDO

      CALL unpack_arr( jpoce, flxsedi3d(PRIV_2D_BIOARRAY,1:jpksed,1)  , iarroce(1:jpoce), &
         &                   zdta(1:jpoce,1:jpksed)  )
      
      CALL unpack_arr( jpoce, flxsedi3d(PRIV_2D_BIOARRAY,1:jpksed,2)  , iarroce(1:jpoce), &
         &                   co3por(1:jpoce,1:jpksed)  )
      
!      flxsedi3d = 0.
      zflx(:,:) = 0.    
      ! Calculation of fluxes mol/cm2/s
      DO jw = 1, jpwat
         DO ji = 1, jpoce
            zflx(ji,jw) = ( pwcp(ji,1,jw) - pwcp_dta(ji,jw) ) &
               &         * 1.e3 / 1.e2 * dzkbot(ji) / rfact
         ENDDO
      ENDDO

      ! Calculation of accumulation rate per dt
      DO js = 1, jpsol
         zrate =  1.0 / ( denssol * por1(jpksed) ) / rfact
         DO ji = 1, jpoce
            zflx(ji,jpwatp1) = zflx(ji,jpwatp1) + ( tosed(ji,js) - fromsed(ji,js) ) * zrate
         ENDDO
      ENDDO

      DO jn = 1, jpdia2dsed - 1 
         CALL unpack_arr( jpoce, flxsedi2d(PRIV_2D_BIOARRAY,jn), iarroce(1:jpoce), zflx(1:jpoce,jn)  )
      END DO
      zflx(:,1) = dzdep(:) / dtsed
      CALL unpack_arr( jpoce, flxsedi2d(PRIV_2D_BIOARRAY,jpdia2dsed), iarroce(1:jpoce), zflx(1:jpoce,1) )

       ! Start writing data
       ! ---------------------
       DO jn = 1, jptrased
          DO jk = 1, jpksed
             DO jj = JRANGE
                DO ji = IRANGE
                   trcsedi(ji,jj,jk,jn) = ztrcsedi(ji,jj,jk,jn)
                END DO
             END DO
          END DO
       END DO

       DO jn = 1, jptrased
          cltra = sedtrcd(jn) ! short title for 3D diagnostic
          CALL iom_put( cltra, trcsedi(:,:,:,jn) )
       END DO

       DO jn = 1, jpdia3dsed
          DO jk = 1, jpksed
             DO jj = JRANGE
                DO ji = IRANGE
                   trcsedi(ji,jj,jk,jn) = flxsedi3d(ji,jj,jk,jn)
                END DO
             END DO
          END DO
       END DO

       DO jn = 1, jpdia3dsed
          cltra = seddia3d(jn) ! short title for 3D diagnostic
          CALL iom_put( cltra, trcsedi(:,:,:,jn) )
       END DO

       DO jn = 1, jpdia2dsed
          DO jj = JRANGE
             DO ji = IRANGE
                trcsedi(ji,jj,1,jn) = flxsedi2d(ji,jj,jn)
             END DO
          END DO
       END DO
 
       DO jn = 1, jpdia2dsed
          cltra = seddia2d(jn) ! short title for 2D diagnostic
          CALL iom_put( cltra, trcsedi(:,:,1,jn) )
       END DO

      DEALLOCATE( zdta, zflx, flxsedi2d, flxsedi3d, trcsedi, ztrcsedi ) 

   END SUBROUTINE sed_wri

#endif
#endif


END MODULE sedwri
