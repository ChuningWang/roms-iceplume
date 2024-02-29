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

MODULE sedrst

#if defined key_pisces

   !! * Modules used
   USE sed
   USE sedarr
   USE sms_pisces, ONLY : rtrn
#ifdef AGRIF
      USE param, ONLY : Lmmpi,Mmmpi
#endif
   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC sed_rst_wri         ! routine called by opa.F90
   PUBLIC sed_rst_read

   !!* Substitution
#  include "ocean2pisces.h90"

CONTAINS         ! Write model prognostic

      SUBROUTINE sed_rst_wri      ! variables into restart
                                  ! netCDF file.
# include "netcdf.inc"

      INTEGER :: ierr, record, lstr, lvar, lenstr   &
      &  , start(2), count(2), ibuff(4), nf_fwrite, itrc  
      INTEGER :: ji, jj, jk, jn
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: trcsedtmp, trcsedi
      REAL(wp), DIMENSION(jpoce,jpksed)   :: zdta


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

#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      IF (mynode > 0) THEN
         call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,       &
         &                 1, MPI_COMM_WORLD, status, ierr)
      ENDIF
#endif
!
! Create/open restart file; write grid arrays, if so needed.
!
      CALL def_rst_sed (ncidrstsed, nrecrst, ierr)
      IF (ierr .NE. nf_noerr) GOTO 99
      lstr = lenstr(rstname)
!                                            !!! WARNING: Here it is
! Set record within the file.                !!! assumed that global
!                                            !!! restart record index 
      nrecrst = max(nrecrst,1)                 !!! nrecrst is already
      IF (nrpfrst == 0) THEN                 !!! advanced by main.
         record = nrecrst
      ELSE
         record = 1+mod(nrecrst-1, abs(nrpfrst))
      ENDIF

!
! Write out evolving model variables:
! ----- --- -------- ----- ----------
!
! Time step number and record indices. 
!
      ibuff(1) = iic
      ibuff(2) = nrecrst
      ibuff(3) = nrechis
#ifdef AVERAGES
      ibuff(4) = nrecavg
#else
      ibuff(4) = 0
#endif
      start(1) = 1
      start(2) = record
      count(1) = 4
      count(2) = 1
      ierr = nf_put_vara_int (ncidrstsed, rstTstep, start, count, ibuff)
      IF (ierr .NE. nf_noerr) THEN
         WRITE(stdout,1) 'time_step', record, ierr      
         GOTO 99                                           !--> ERROR
      ENDIF
!
! Time.
!
      ierr = nf_put_var1_FTYPE (ncidrstsed, rstTime, record, time)
      IF (ierr .NE. nf_noerr) THEN
         lvar = lenstr(vname(1,indxTime))
         WRITE(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
         GOTO 99                                           !--> ERROR
      ENDIF
!
! Tracer variables.
!
!
      ALLOCATE(trcsedtmp(GLOBAL_2D_ARRAY,jpksed,jptrased), trcsedi(PRIV_2D_BIOARRAY,jpksed,jptrased) )
      trcsedi(:,:,:,:)   = 0.0
      trcsedtmp(:,:,:,:) = 0.0
!
!      ! Back to 2D geometry
      DO jn = 1, jpsol
         CALL unpack_arr( jpoce, trcsedi(PRIV_2D_BIOARRAY,1:jpksed,jn) , iarroce(1:jpoce), &
         &                       solcp(1:jpoce,1:jpksed,jn ) )
      END DO
!
      DO jn = 1, jpwat
         CALL unpack_arr( jpoce, trcsedi(PRIV_2D_BIOARRAY,1:jpksed,jpsol+jn) , iarroce(1:jpoce), &
         &                       pwcp(1:jpoce,1:jpksed,jn  )  )
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
         ierr = nf_fwrite(trcsedtmp(START_2D_ARRAY,1,itrc), ncidrstsed,   &
         &                             rstsed(itrc), record, r3dsed)
        IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) TRIM(sedtrcd(itrc)), record, ierr
          GOTO 99                                         !--> ERROR
        ENDIF
        END DO

      DEALLOCATE(trcsedtmp, trcsedi )
!
! Additional variables.
!

      ALLOCATE(trcsedtmp(GLOBAL_2D_ARRAY,jpksed,5), trcsedi(PRIV_2D_BIOARRAY,jpksed,5) )

      trcsedi(:,:,:,:)   = 0.0
      trcsedtmp(:,:,:,:) = 0.0

      ! pH
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
      &                   db(1:jpoce,1:jpksed)  )

      CALL unpack_arr( jpoce, trcsedi(PRIV_2D_BIOARRAY,1:jpksed,4)  , iarroce(1:jpoce), &
      &                   irrig(1:jpoce,1:jpksed)  )

      CALL unpack_arr( jpoce, trcsedi(PRIV_2D_BIOARRAY,1:jpksed,5)  , iarroce(1:jpoce), &
      &                   sedligand(1:jpoce,1:jpksed)  )

      DO itrc = 1, 5
         DO jk = 1, jpksed
            DO jj = 1, LOCALMM
               DO ji = 1, LOCALLM
                  trcsedtmp(ji,jj,jk,itrc) = trcsedi(ji,jj,jk,itrc)
               END DO
            END DO
         END DO
      END DO

      DO itrc = 1, 5
         ierr=nf_fwrite(trcsedtmp(START_2D_ARRAY,1,itrc), ncidrstsed,   &
         &                             rstadd(itrc), record, r3dsed)
         IF (ierr .NE. nf_noerr) THEN
            WRITE(stdout,1) TRIM(sname(itrc,1)), record, ierr
            GOTO 99                                         !--> ERROR
         ENDIF
      END DO

      DEALLOCATE(trcsedtmp, trcsedi )


  1   FORMAT(/1x, 'SED_RST_WRI ERROR while writing variable ''', A,   &
       &           ''' into restart file.', /11x, 'Time record:', &
       &               i6, 3x, 'netCDF error code', i4, 3x, A,i4) 
      GOTO 100 
  99  may_day_flag=3
 100  CONTINUE

!
! Synchronize restart netCDF file to disk to allow other
! processes to access data immediately after it is written.
!
#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      ierr = nf_close (ncidrstsed)
      IF (nrpfrst > 0 .AND. record >= nrpfrst) ncidrstsed = -1
#else
      IF (nrpfrst > 0 .AND. record >= nrpfrst) THEN
        ierr = nf_close (ncidrstsed)
        ncidrstsed = -1
      ELSE
        ierr = nf_sync(ncidrstsed)
      ENDIF
#endif
      IF (ierr == nf_noerr) THEN
         MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')    & 
         &            'WRT_RST_SED -- wrote ',                          &
         &            'restart fields into time record =', record, '/',  &
         &             nrecrst  
      ELSE
         MPI_master_only  write(stdout,'(/1x,2A/)')     & 
         &             'WRT_RST_SED ERROR: Cannot ',        &
         &             'synchronize/close restart netCDF file.'
         may_day_flag = 3
      ENDIF

#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      IF (mynode < NNODES-1) THEN
         CALL MPI_Send (blank, 1, MPI_INTEGER, mynode+1, 1, MPI_COMM_WORLD,  ierr)
      ENDIF
#endif
      RETURN
      END SUBROUTINE sed_rst_wri 

      SUBROUTINE def_rst_sed( ncid, total_rec, ierr)  ! restart netCDF

# include "netcdf.inc"

      logical :: create_new_file
      integer :: ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim    &
      &      , r2dgrd(3),  auxil(2),  checkdims                           &
#ifdef NC4PAR
      &      , csize,cmode           &
#endif
      &      , r3dgrd(4),  u3dgrd(4), v3dgrd(4),  w3dgrd(4), itrc
      CHARACTER(len=20) :: cltra

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
      ierr = 0
      lstr = lenstr(cn_sedrst_out)
      IF (nrpfrst > 0) THEN
         lvar=total_rec - (1+mod(total_rec-1, nrpfrst))
         call insert_time_index (cn_sedrst_out, lstr, lvar, ierr)
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
      create_new_file = ldefhis
      IF (ncid .NE. -1) create_new_file = .false.
#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      IF (mynode > 0) create_new_file = .false.
#endif
!
! Create new restart file:    Put global attributes
!======= === ======= =====    and define all variables.
!
  10  if (create_new_file) then

#ifndef NC4PAR
        ierr  = nf_create(TRIM(cn_sedrst_out),NF_CLOBBER, ncid)
#else
        cmode = ior(nf_netcdf4,nf_classic_model)
        cmode = ior(cmode, nf_mpiio)
        csize = xi_rho*eta_rho/NNODES
        WRITE(stdout,*)'CREATE RST NC4 PARALLEL FILE'
        ierr  = nf_create_par(cn_sedrst_out(1:lstr),cmode, &
        &       MPI_COMM_WORLD,MPI_INFO_NULL,ncid)
#endif

        IF (ierr .NE. nf_noerr) THEN
           WRITE(stdout,'(/3(1x,A)/)') 'ERROR in DEF_RST_SED: Cannot',    &
           &             'create restart NetCDF file:', TRIM(cn_sedrst_out)
           GOTO 99                                         !--> ERROR
        ENDIF
        IF (nrpfrst == 0) total_rec = 0
!
! Put global attributes.
! --- ------ -----------
!
        CALL put_global_atts (ncid, ierr)
!
! Define dimensions of staggered fields.
! ------ ---------- -- --------- -------
!
        ierr = nf_def_dim (ncid, 'xi_rho',   xi_rho,  r2dgrd(1))
        ierr = nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr = nf_def_dim (ncid, 'profsed',    jpksed,   r3dgrd(3))
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
        &       rstTstep)
#ifdef NC4PAR
        ierr = nf_var_par_access(ncid,rstTstep,nf_collective)
#endif
        ierr = nf_put_att_text (ncid, rstTstep, 'long_name', 48,    &
        &       'time step and record numbers from initialization')
!
! Time.
!
        lvar = lenstr(vname(1,indxTime))
        ierr = nf_def_var (ncid, vname(1,indxTime)(1:lvar),           &
        &                              NF_DOUBLE, 1, timedim, rstTime)
#ifdef NC4PAR
        ierr = nf_var_par_access(ncid,rstTime,nf_collective)
#endif
        lvar = lenstr(vname(2,indxTime))
        ierr = nf_put_att_text (ncid, rstTime, 'long_name', lvar,     &
        &                                  vname(2,indxTime)(1:lvar))
        lvar = lenstr(vname(3,indxTime))
        ierr = nf_put_att_text (ncid, rstTime, 'units',     lvar,     &
        &                                  vname(3,indxTime)(1:lvar))
        lvar = lenstr (vname(4,indxTime))
        ierr = nf_put_att_text(ncid, rstTime, 'field',     lvar,      &
        &                                  vname(4,indxTime)(1:lvar))

        CALL nf_add_attribute(ncid, rstTime, indxTime, 5, NF_DOUBLE, ierr)
!
! Time2.
!
        lvar = lenstr(vname(1,indxTime2))
        ierr = nf_def_var (ncid, vname(1,indxTime2)(1:lvar),            &
        &                              NF_DOUBLE, 1, timedim, rstTime2)
#ifdef NC4PAR
        ierr = nf_var_par_access(ncid,rstTime2,nf_collective)
#endif
        lvar = lenstr(vname(2,indxTime2))
        ierr = nf_put_att_text (ncid, rstTime2, 'long_name', lvar,     &
        &                                  vname(2,indxTime2)(1:lvar))
        lvar = lenstr(vname(3,indxTime2))
        ierr = nf_put_att_text (ncid, rstTime2, 'units',     lvar,     &
        &                                  vname(3,indxTime2)(1:lvar))
        lvar = lenstr (vname(4,indxTime2))
        ierr = nf_put_att_text(ncid, rstTime2, 'field',     lvar,      &
        &                                  vname(4,indxTime2)(1:lvar))

        CALL nf_add_attribute(ncid, rstTime2, indxTime2, 5, NF_DOUBLE, ierr)
!
! Tracer variables.
!
        DO itrc = 1, jptrased
           cltra = TRIM(sedtrcd(itrc))
           ierr  = nf_def_var (ncid, cltra, NF_DOUBLE, 4, r3dgrd, rstsed(itrc))
#ifdef NC4PAR
           ierr = nf_var_par_access(ncid,rstsed(itrc),nf_collective)
#endif
           lvar = lenstr(TRIM(sedtrcl(itrc)))
           ierr = nf_put_att_text (ncid, rstsed(itrc), 'long_name',    &
           &                     lvar, TRIM(sedtrcl(itrc)))
           lvar = lenstr(TRIM(sedtrcu(itrc)))
           ierr = nf_put_att_text (ncid, rstsed(itrc), 'units', lvar, TRIM(sedtrcu(itrc)))
           CALL nf_add_attribute(ncid, rstsed(itrc), itrc, 5, NF_DOUBLE,ierr)
        END DO

        DO itrc = 1, 5
           cltra = sname(itrc,1)
           ierr = nf_def_var (ncid, cltra, NF_DOUBLE, 4, r3dgrd, rstadd(itrc))
#ifdef NC4PAR
           ierr = nf_var_par_access(ncid,rstadd(itrc),nf_collective)
#endif
           lvar = lenstr(sname(itrc,2))
           ierr = nf_put_att_text (ncid, rstadd(itrc), 'long_name',    &
           &                     lvar, sname(itrc,2) )
           lvar = lenstr(sname(itrc,3))
           ierr = nf_put_att_text (ncid, rstadd(itrc), 'units', lvar, sname(itrc,3))
           CALL nf_add_attribute(ncid, rstadd(itrc), itrc, 5, NF_DOUBLE,ierr)
        END DO

!
! Leave definition mode.                  Also initialize record
! ----- ---------- -----                  dimension size to zero.
!
        ierr = nf_enddef(ncid)
        WRITE(*,'(6x,4A,1x,A,i4)') 'DEF_RST_SED - Created new ',        &
        &              'netCDF file ''', TRIM(cn_sedrst_out), '''.'
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
#ifndef NC4PAR
        ierr = nf_open (TRIM(cn_sedrst_out), nf_write, ncid)
#else
        ierr = nf_open_par (TRIM(cn_sedrst_out), IOR(nf_write, nf_mpiio),   &
        &     MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
#endif
        IF (ierr == nf_noerr) THEN
           ierr = checkdims (ncid, cn_sedrst_out, lstr, rec)
           IF (ierr == nf_noerr) THEN
              IF (nrpfrst == 0) THEN
                 ierr = rec+1 - nrecrst
              ELSE
                 ierr = rec+1 - (1+mod(nrecrst-1, abs(nrpfrst)))
              ENDIF
              IF (ierr > 0) THEN
                 MPI_master_only write( stdout,                              &
        &                 '(/1x,A,I5,1x,A/8x,3A,I5/8x,A,I5,1x,A/)'         &
        &           ) 'DEF_RST_SED WARNING: Actual number of records', rec,    &
        &             'in netCDF file',  '''',  TRIM(cn_sedrst_out),       &
        &             ''' exceeds the record number from restart data',    &
        &             rec+1-ierr,'/', total_rec,', restart is assumed.'
                 rec = rec-ierr
              ELSEIF (nrpfrst == 0) THEN
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
              WRITE(stdout,'(/1x,4A, 1x,A,I4/)')     'DEF_RST_SED ERROR: ',    &
              &     'Cannot open restart netCDF file ''',TRIM(cn_sedrst_out),'''.'
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
        ierr = nf_inq_varid (ncid, 'time_step', rstTstep)
        IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) 'time_step', TRIM(cn_sedrst_out)
          GOTO 99                                         !--> ERROR
        ENDIF
!
! Time.
!
        lvar = lenstr(vname(1,indxTime))
        ierr = nf_inq_varid (ncid, vname(1,indxTime)(1:lvar), rstTime)
        IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) vname(1,indxTime)(1:lvar), TRIM(cn_sedrst_out)
          GOTO 99                                         !--> ERROR
        ENDIF
!
! Time2.
!
        lvar = lenstr(vname(1,indxTime2))
        ierr = nf_inq_varid (ncid, vname(1,indxTime2)(1:lvar), rstTime2)
        IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) vname(1,indxTime2)(1:lvar), TRIM(cn_sedrst_out)
          GOTO 99                                         !--> ERROR
        ENDIF
!
! Tracer variables.
!
       DO itrc = 1, jptrased
          ierr = nf_inq_varid (ncid, TRIM(sedtrcd(itrc)), rstsed(itrc))
          IF (ierr .NE. nf_noerr) THEN
             WRITE(stdout,1) TRIM(sedtrcd(itrc)), cn_sedrst_out(1:lstr)
             GOTO 99                                       !--> ERROR
          ENDIF
       END DO
!
! Additional variables.
!
       DO itrc = 1, 5
          ierr = nf_inq_varid (ncid, TRIM(sname(itrc,1)),rstadd(itrc))
          IF (ierr .NE. nf_noerr) THEN
            WRITE(stdout,1) TRIM(sname(itrc,1)),cn_sedrst_out(1:lstr)
            GOTO 99                                       !--> ERROR
          ENDIF
       ENDDO
!
        MPI_master_only WRITE(*,'(6x,2A,i4,1x,A,i4)')              &
        &             'DEF_RST_SED -- Opened ',                        &
        &             'existing restart file,  record =', rec

#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      ELSE
         ierr = nf_open (cn_sedrst_out(1:lstr), nf_write, ncid)
         IF (ierr .NE. nf_noerr) THEN
            MPI_master_only WRITE(stdout,'(/1x,4A, 1x,A,I4/)')        &
            &          'DEF_RST_SED ERROR: Cannot',                       &
            &          'open restart netCDF file ''', cn_sedrst_out(1:lstr), '''.'
            GOTO 99                                         !--> ERROR
         ENDIF
#endif
      ENDIF              !<-- create_new_file
   1  FORMAT(/1x,'DEF_RST_SED ERROR: Cannot find variable ''',        &
      &               A, ''' in netCDF file ''', A, '''.'/)
  99  RETURN                                              !--> ERROR
      END SUBROUTINE def_rst_sed

                              ! Read initial conditions for the
      SUBROUTINE sed_rst_read ! primitive variables from NetCDF
                              ! initialization file.

!======================================================
!
!======================================================

# include "netcdf.inc"

      real(wp) :: time_scale
      integer  :: itrc
      integer  :: ji, jj, jk, jn
      integer  :: ncid, indx, varid,  ierr, lstr, lvar, latt, lenstr,    &
      &        start(2), count(2), ibuff(6), nf_fread, checkdims
      character :: units*180
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: trcsedtmp
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: zdta
      REAL(wp), DIMENSION(jpoce,jpksed)         :: zhipor


#ifdef MPI
#define LOCALLM Lmmpi
#define LOCALMM Mmmpi
#else
#define LOCALLM Lm
#define LOCALMM Mm
#endif

!
! Open initial conditions netCDF file for reading. Check that all
! spatial dimensions in that file are consistent with the model
! arrays, determine how many time records are available in the file
! and set record from which the dada will be read.
!
! The record is set as follows: (1) if only one time record is
! available in the file, then that record is used REGARDLESS of
! value of nrrec supplied from the parameter file; (2) if the
! file has multiple records and nrrec is positive, then nrrec is
! available record is used.
!
      IF (may_day_flag .NE. 0) RETURN      !-->  EXIT
      lstr = lenstr(cn_sedrst_in)
      ierr = nf_open(TRIM(cn_sedrst_in), nf_nowrite, ncid)
      IF (ierr == nf_noerr) THEN
         IF (ierr .NE. nf_noerr) THEN
            GOTO 99
         ELSEIF (indx == 0) then
            indx = 1
         ELSEIF (indx > 0 .AND. nrrec > 0 .AND. nrrec <= indx) THEN
            indx = nrrec
         ELSEIF (indx > 0 .AND. nrrec > indx) THEN
            WRITE(stdout,'(/1x,A,I4,A/16x,A,I4,A/16x,3A/)')                   &
            &            'SED_RST_READ ERROR: requested restart time record',  &
            &             nrrec, ' exceeds',  'number of available records',  &
            &             indx,'in netCDF file', '''',TRIM(cn_sedrst_in),'''.'
            GOTO 99                                        !--> ERROR
         ENDIF
      ELSE
         WRITE(stdout,'(/1x,2A/15x,3A)') 'SED_RST_READ ERROR: Cannot ',      &
         &               'open netCDF file', '''', TRIM(cn_sedrst_in) ,'''.'
         GOTO 99                                           !--> ERROR
      ENDIF
!
! Read in evolving model variables:
! ---- -- -------- ----- ----------
!
! Time: find netCDF id, read value, read attribute 'units'
! and set starting time index and time clock in days.
!
      lvar = lenstr(vname(1,indxTime))
      ierr = nf_inq_varid (ncid, vname(1,indxTime)(1:lvar), varid)
      IF (ierr == nf_noerr) THEN
        ierr = nf_get_var1_FTYPE (ncid, varid, indx, time)
        IF (ierr == nf_noerr) THEN
          ierr = nf_get_att_text(ncid, varid, 'units', units)
          IF (ierr == nf_noerr) THEN
            latt = lenstr(units)
            IF (units(1:6) == 'second') THEN
               time_scale = 1.
            ELSEIF (units(1:3) == 'day') THEN
              time_scale = day2sec
            ELSE
              WRITE (stdout,'(/1x,4A/8x,3A/)') 'SED_RST_READ ',      &
       &              'ERROR: unknown units of for variable ''',      &
       &               vname(1,indxTime)(1:lvar), '''',               &
       &              'in netCDF file ''', TRIM(cn_sedrst_in),'''.'
              GOTO 99                                    !--> ERROR
            ENDIF
          ELSE
            WRITE (stdout,'(/1x,2A/8x,5A/)') 'SED_RST_READ ERROR: ',   &
       &             'cannot read attribute ''units'' for variable',  &
       &             '''', vname(1,indxTime)(1:lvar),                 &
       &             ''' in netCDF file ''',  TRIM(cn_sedrst_in), '''.'
            GOTO 99                                       !--> ERROR
          ENDIF
        ELSE
          MPI_master_only write(stdout,2) vname(1,indxTime)(1:lvar)  &
          &                                , indx, TRIM(cn_sedrst_in)
          GOTO 99                                         !--> ERROR
        ENDIF
      ELSE
        MPI_master_only write(stdout,1) vname(1,indxTime)(1:lvar), TRIM(cn_sedrst_in)
        GOTO 99                                           !--> ERROR
      ENDIF

      time = time*time_scale
      tdays = time*sec2day

      ierr = nf_inq_varid (ncid, 'time_step', varid)
      IF (ierr == nf_noerr) THEN
         start(1) = 1
         start(2) = indx
         count(1) = 4
         count(2) = 1
         ierr = nf_get_vara_int (ncid, varid, start, count, ibuff)
         IF (ierr == nf_noerr) THEN
            ntstart = ibuff(1)
            nrecrst = ibuff(2)
            nrechis = ibuff(3)

            MPI_master_only WRITE(stdout,                            &
            &     '(6x,A,G12.4,A,I2,A,I6,A,I3,A,I3,A)')              &
            &     'SED_RST_READ: Restarted from day =', tdays, ' rec =',   &
            &      indx, '(', ntstart, ',', nrecrst, ',', nrechis, ').'

         ELSE
            MPI_master_only write(stdout,'(/1x,2A/)')                     &
            &                            'SED_RST_READ ERROR: Cannot ',    &
            &                            'read time and record indices.'
            GOTO 99                                         !--> ERROR
         ENDIF
      ELSE
         ntstart = 1
         nrecrst = 0
         nrechis = 0
         MPI_master_only WRITE(stdout,'(6x,2A,G12.4,1x,A,I4)')      &
         &          'SED_RST_READ -- ',                              &
         &          'Processing data for time =', tdays, 'record =', indx
      ENDIF
      IF (ntstart < 1) ntstart = 1
      ntimes = ntstart+ntimes-1
!
! Tracer variables.
!
      ALLOCATE(trcsedtmp(GLOBAL_2D_ARRAY,jpksed,jptrased), zdta(PRIV_2D_BIOARRAY,jpksed,jptrased) )

      zdta(:,:,:,:) = 0.0
      trcsedtmp(:,:,:,:) = 0.0


      DO itrc = 1, jptrased
        ierr = nf_inq_varid (ncid, TRIM(sedtrcd(itrc)), varid)
        IF (ierr == nf_noerr) THEN
          ierr = nf_fread (trcsedtmp(START_2D_ARRAY,1,itrc), ncid,  varid, indx, r3dsed)
          IF (ierr .NE. nf_noerr) THEN
            MPI_master_only WRITE(stdout,2) TRIM(sedtrcd(itrc)), indx, TRIM(cn_sedrst_in)
            GOTO 99                                       !--> ERROR
          ENDIF
        ELSE
           MPI_master_only WRITE(stdout,3) TRIM(sedtrcd(itrc)), TRIM(cn_sedrst_in)
        ENDIF
      END DO

      DO itrc = 1, jptrased
         DO jk = 1, jpksed
            DO jj = 1, LOCALMM
               DO ji = 1, LOCALLM
                  zdta(ji,jj,jk,itrc) = trcsedtmp(ji,jj,jk,itrc)
               END DO
            END DO
         END DO
      END DO

      DO jn = 1, jpsol
         CALL pack_arr( jpoce, solcp(1:jpoce,1:jpksed,jn), &
         &              zdta(1:jpi,1:jpj,1:jpksed,jn), iarroce(1:jpoce) )
      END DO

      DO jn = 1, jpwat
         CALL pack_arr( jpoce, pwcp(1:jpoce,1:jpksed,jn), &
         &              zdta(1:jpi,1:jpj,1:jpksed,jpsol+jn), iarroce(1:jpoce) )
      END DO
      DEALLOCATE( zdta, trcsedtmp )

      ! Initialization of sediment composant only ie jk=2 to jk=jpksed
      ! ( nothing in jk=1)
      solcp(1:jpoce,1,:) = 0.
      pwcp (1:jpoce,1,:) = 0.
!
! Additional variables.
!
      ALLOCATE(trcsedtmp(GLOBAL_2D_ARRAY,jpksed,5), zdta(PRIV_2D_BIOARRAY,jpksed,5) )

      zdta(:,:,:,:) = 0.0
      trcsedtmp(:,:,:,:) = 0.0

      DO itrc = 1, 5
        ierr = nf_inq_varid (ncid, TRIM(sname(itrc,1)), varid)
        IF (ierr == nf_noerr) THEN
          ierr = nf_fread (trcsedtmp(START_2D_ARRAY,1,itrc), ncid,  varid,  &
          &                                               indx, r3dsed)
          IF (ierr .NE. nf_noerr) THEN
            MPI_master_only WRITE(stdout,2) TRIM(sname(itrc,1)), indx, TRIM(cn_sedrst_in)
            GOTO 99                                       !--> ERROR
          ENDIF
        ELSE
           MPI_master_only WRITE(stdout,3) TRIM(sname(itrc,1)), TRIM(cn_sedrst_in)
        ENDIF
      ENDDO

      DO itrc = 1, 5
         DO jk = 1, jpksed
            DO jj = 1, LOCALMM
               DO ji = 1, LOCALLM
                  zdta(ji,jj,jk,itrc) = trcsedtmp(ji,jj,jk,itrc)
               END DO
            END DO
         END DO
      END DO

      zhipor(:,:) = 0.
      CALL pack_arr( jpoce, zhipor(1:jpoce,1:jpksed), &
      &             zdta(PRIV_2D_BIOARRAY,1:jpksed,1), iarroce(1:jpoce) )

      ! Initialization of [h+] in mol/kg
      DO jk = 1, jpksed
         DO ji = 1, jpoce
            hipor (ji,jk) = 10.**( -1. * zhipor(ji,jk) )
         ENDDO
      ENDDO

      CALL pack_arr( jpoce, co3por(1:jpoce,1:jpksed), &
      &             zdta(PRIV_2D_BIOARRAY,1:jpksed,2), iarroce(1:jpoce) )

      CALL pack_arr( jpoce, db(1:jpoce,1:jpksed), &
      &             zdta(PRIV_2D_BIOARRAY,1:jpksed,3), iarroce(1:jpoce) )

      CALL pack_arr( jpoce, irrig(1:jpoce,1:jpksed), &
      &             zdta(PRIV_2D_BIOARRAY,1:jpksed,4), iarroce(1:jpoce) )

      CALL pack_arr( jpoce, sedligand(1:jpoce,1:jpksed), &
      &             zdta(PRIV_2D_BIOARRAY,1:jpksed,5), iarroce(1:jpoce) )

      DEALLOCATE( zdta, trcsedtmp )

!======================================================
! END MODIF_JG_2
!======================================================

!
!  Close input NetCDF file.
!
      ierr = nf_close(ncid)

  1   FORMAT(/1x,'SED_RST_READ - unable to find variable:',    1x,A,    &
      &                            /15x,'in input NetCDF file:',1x,A/)
  2   FORMAT(/1x,'SED_RST_READ - error while reading variable:',1x, A,  &
      &    2x,'at time record =',i4/15x,'in input NetCDF file:',1x,A/)
  3   FORMAT(/1x,'SED_RST_READ - unable to find variable:',    1x,A,    &
      &                            /15x,'in input NetCDF file:',1x,A,  &
      &    1x,'-> analytical value'/)
      RETURN
  99  may_day_flag = 2
      RETURN
      END SUBROUTINE sed_rst_read

#endif

END MODULE sedrst
