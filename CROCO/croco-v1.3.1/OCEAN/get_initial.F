! $Id: get_initial.F 1458 2014-02-03 15:01:25Z gcambon $
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
                                                ! Read initial conditions for the
      subroutine get_initial (req_rec, tindx)   ! primitive variables from NetCDF
                                                ! initialization file.
      implicit none
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "forces.h"
#ifdef SEDIMENT
# include "sediment.h"
#endif
#ifdef BBL
# include "bbl.h"
#endif
#if defined OA_COUPLING || defined OW_COUPLING
# include "mpi_cpl.h"
#endif
#if defined NBQ || defined MORPHODYN
# include "grid.h"
#endif
# include "mixing.h"

# if defined EXACT_RESTART
# include "coupling.h"
# if defined M3FAST
# include "nbq.h"
# endif
# endif
!======================================================
!
!======================================================

      real time_scale
      integer itrc
#ifdef EXACT_RESTART
      real time_bak
#endif
#if defined SEDIMENT || defined BBL
     &      ,  indxWrk
#endif
      integer ncid, req_rec, tindx, record, max_rec,
     &        indx, varid,  ierr, lstr, lvar, latt, lenstr,
     &        start(2), count(2), ibuff(6),   nf_fread, checkdims
      character units*180
      CHARACTER(len=19)          :: tool_sectodat
#include "netcdf.inc"
#define time illegal
#define tdays illegal
#define nrrec illegal

!     
! Open initial conditions netCDF file for reading. Check that all
! spatial dimensions in that file are consistent with the model
! arrays, determine how many time records are available in the file
! and determine number of record from which the data will be read.
! The record is set as follows: (1) if only one time record is
! available in the file, then that record is used REGARDLESS of value
! of "nrrec" supplied from the parameter file; (2) if the file has
! multiple records and "nrrec" is positive, then "nrrec" is used,
! provided that "nrrec" is within the available records; (3) if the
! file has multiple records and nrrec<0, then THE LAST available
! record is used.
!
      if (may_day_flag.ne.0) return      !-->  EXIT
      lstr=lenstr(ininame)
      ierr=nf_open(ininame(1:lstr), nf_nowrite, ncid)
      if (ierr .eq. nf_noerr) then
        ierr=checkdims (ncid, ininame, lstr, max_rec)
        if (ierr. eq. nf_noerr) then
          if (max_rec.gt.0) then
            if (req_rec.gt.0) then
              if (req_rec.le.max_rec) then
                record=req_rec
              else
       MPI_master_only write(stdout,'(/1x,2A,I4,1x,A/12x,A,I4,1x,3A/)')
     &           '### ERROR: get_initial :: requested restart time ',
     &           'record',  req_rec, 'exceeds number',  'of records',
     &                       max_rec,  'available in netCDF file ''',
     &                                        ininame(1:lstr), '''.'

                goto 99         !--> ERROR
              endif
            else
              record=max_rec
            endif
          else
            record=1
          endif
        else
          goto 99
        endif
      else
        MPI_master_only write(stdout,'(/1x,4A/12x,A/)')
     &               '### ERROR: get_initial :: ',
     &               'Cannot open netCDF file ''',   ininame(1:lstr),
     &                                  '''.',     nf_strerror(ierr)
        goto 99                                           !--> ERROR
      endif
!
! Read in evolving model variables:
!----- -- -------- ----- ----------
! Time: find netCDF id, read value, read attribute 'units' and set
! starting time index and time clock in days.
!
! Note that if EXACT_RESTART CPP-switch is defined, make a "soft"
! attempt to do exact restart, where "soft" means that exact restart
! is done only when file of initial conditions contains sufficient
! data, i.e. two consecutive time records of evolving fields one time
! step apart from each other.  Thus, in order to accept the file for
! exact restart, it must pass two consecutive checks: (i) ocean_time
! values in two consecutive values must differ by "dt" of current
! run, and (ii) "ntstart" from two consecutive records of netCDF
! structure "time_step" must differ by one.  If either check
! is fails, forward step is used as the initial time step.  Since  
! "get_initial" is expected to be called twice consecutively as
!
!           call get_initial (req_rec=rec-1, tindx=2)
!           call get_initial (req_rec=rec,   tindx=1)
!
! where "rec" is record number in netCDF file which contains fields 
! corresponding to time step "n" while "rec-1" corresponds to "n-1"
! (hence, making it possible to start time stepping with regular LF
! predictor step rather than forward), both checks are performed
! during the first call, tindx=2. If either check fails, the exact
! restart is cancelled and no reading of 2D and 3D fields will be
! performed for tindx=2.
!
! The possibility of exact restart is communicated with the rest
! of the code via integer variable "forw_start"  which is set
! exclussively by this routine and is used as part of CPP-macro
!
!             FIRST_TIME_STEP iic.eq.forw_start
!
! where the possibilities are as follows:
!
!       forw_start=1 means that "exact" restart, i.e., forward step
!                    is to be performed only during absolutely first
!                    time step, iic=1, and
!
!       forw_start=ntstart, means that restarted time stepping
!                           should also begin with forward step
!                           (this is approximate restart). 
!
! This mechanism of exact restart is designed to handle essentially
! three situations: (1) initial run using a 3rd-party file which
! contain initial time, but does not contain "time_step". In this
! case ntstart is set to 1, and forward step is assumed at the first
! step; (2) restart from restart file generated by this code, but
! with deactivated CPP-switch EXACT_RESTART.  This file contains both
! both time variable and "time_step" structure, but only one
! consecutive record. This situation is identified automatically and
! approximate restart is assumed. This is compatibility mode. This
! also incldes restart from a history file generated by this code.  
! (3) restart from a file created by this code with activated
! EXACT_RESTART.
!
      lvar=lenstr(vname(1,indxTime))
      ierr=nf_inq_varid (ncid, vname(1,indxTime)(1:lvar), varid)
      
!! GC
!! The following is done for backward compatibility: normally time
!! variable is named "ocean_time", but legacy startup files may name
!! it either "roms_time" or "scrum_time".
!!
!      if (ierr .ne. nf_noerr) then
!        ierr=nf_inq_varid (ncid, 'roms_time', varid)
!      endif
!      if (ierr .ne. nf_noerr) then
!        ierr=nf_inq_varid (ncid, 'scrum_time', varid)
!      endif
! GC

      if (ierr .eq. nf_noerr) then
        ierr=nf_get_var1_FTYPE (ncid, varid, record, start_time)
        if (ierr .eq. nf_noerr) then
          ierr=nf_get_att_text(ncid, varid, 'units', units)
          if (ierr .eq. nf_noerr) then
            latt=lenstr(units)
            if (units(1:6).eq.'second') then
              time_scale=1.
            elseif (units(1:3).eq.'day') then
              time_scale=day2sec
            else
              MPI_master_only write(stdout,'(/1x,4A/8x,3A/)')
     &              'GET_INITIAL ',
     &              'ERROR: unknown units of for variable ''', 
     &                               vname(1,indxTime)(1:lvar), '''',
     &                   'in netCDF file ''', ininame(1:lstr), '''.'
              goto 99                                     !--> ERROR
            endif
            start_time=start_time*time_scale

# if defined XIOS && !defined USE_CALENDAR
      start_date_run = tool_sectodat(start_time+xios_origin_date_in_sec)
      READ(start_date_run(1:4),fmt='(i4)') start_year
      READ(start_date_run(6:7),fmt='(i2)') start_month
      READ(start_date_run(9:10),fmt='(i2)') start_day
      READ(start_date_run(12:13),fmt='(i2)') start_hour
      READ(start_date_run(15:16),fmt='(i2)') start_minute
      READ(start_date_run(18:19),fmt='(i2)') start_second
      !
      READ(xios_origin_date(1:4),fmt='(i4)') origin_year
      READ(xios_origin_date(6:7),fmt='(i2)') origin_month
      READ(xios_origin_date(9:10),fmt='(i2)') origin_day
      READ(xios_origin_date(12:13),fmt='(i2)') origin_hour
      READ(xios_origin_date(15:16),fmt='(i2)') origin_minute
      READ(xios_origin_date(18:19),fmt='(i2)') origin_second
#endif

#if defined USE_CALENDAR
            call tool_origindate(ncid,varid,origin_date_in_sec)
      MPI_master_only write(stdout,'(/2x,A,1x,A/)')
     &    'Origin date of time axis (YYYY-MM-DD hh:mm:ss) :',
     &    tool_sectodat(origin_date_in_sec)
      origin_date=tool_sectodat(origin_date_in_sec)
      READ(origin_date(1:4),fmt='(i4)') origin_year
      READ(origin_date(6:7),fmt='(i2)') origin_month
      READ(origin_date(9:10),fmt='(i2)') origin_day
      READ(origin_date(12:13),fmt='(i2)') origin_hour
      READ(origin_date(15:16),fmt='(i2)') origin_minute
      READ(origin_date(18:19),fmt='(i2)') origin_second
      ! Modify units when using USE_CALENDAR
      vname(3,indxTime)='seconds since '
     &                     //tool_sectodat(origin_date_in_sec)
              ! Modify units when using USE_CALENDAR
      vname(3,indxTime2)='seconds since '
     &                     //tool_sectodat(origin_date_in_sec)
#  if defined XIOS
      MPI_master_only write(stdout,'(/3x,A,/)')
     &    '!!! WARNING: Using XIOS with USE_CALENDAR.',
     &    ' USING USE_CALENDAR origin date for XIOS ',
     &    ' (instead of xios_origin_date)'

      READ(origin_date(1:4),fmt='(i4)') start_year
      READ(origin_date(6:7),fmt='(i2)') start_month
      READ(origin_date(9:10),fmt='(i2)') start_day
      READ(origin_date(12:13),fmt='(i2)') start_hour
      READ(origin_date(15:16),fmt='(i2)') start_minute
      READ(origin_date(18:19),fmt='(i2)') start_second
#endif
#endif

#ifdef EXACT_RESTART
            if (tindx.eq.2) then
              forw_start=0
              if (record.lt.max_rec) then
                time_bak=start_time
                ierr=nf_get_var1_FTYPE (ncid, varid, record+1,
     &                                             start_time)
                if (ierr .eq. nf_noerr) then
                  start_time=start_time*time_scale
!
! Here expression "abs(start_time-time_bak-dt).lt.0.001*dt" below is
! a roundoff-error tolerant version of "start_time.eq.time_bak+dt".
!
                  if (abs(start_time-time_bak-dt) .lt. 0.001*dt) then
                    forw_start=1
                  else
              MPI_master_only write(stdout,'(/1x,2A,2I4/10x,4A/10x,A/)')
     &                'WARNING: Exact restart is requested, but ',
     &                'is not possible: records', record,record+1, 
     &                'in ''',   ininame(1:lstr),   ''' are not ',
     &                'consecutive time steps ==> proceeding ',
     &                          'with forward initial time step.' 
                  endif
                else
                  MPI_master_only write(stdout,2)
     &                       vname(1,indxTime)(1:lvar), record,
     &                       ininame(1:lstr), nf_strerror(ierr)
                  goto 99

                endif
              else
                MPI_master_only write(stdout,'(/1x,2A/10x,4A,A/)')
     &                  'WARNING: Exact restart ',
     &                  'is requested, but is not possible: initial',
     &                  'file ''',  ininame(1:lstr),  ''' does not ',
     &                                 'contain sufficient records.'
              endif
              if (forw_start.ne.1) return
              forw_start=0
            endif
#endif
          else
            MPI_master_only write (stdout,'(/1x,2A/8x,5A/)')
     &             'GET_INITIAL ERROR: ',
     &             'cannot read attribute ''units'' for variable',
     &             '''', vname(1,indxTime)(1:lvar),
     &             ''' in netCDF file ''',  ininame(1:lstr), '''.'
            goto 99                                       !--> ERROR
          endif
        else
          MPI_master_only write(stdout,2) vname(1,indxTime)(1:lvar)
     &                                , record, ininame(1:lstr)
     &                                , nf_strerror(ierr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only write(stdout,1) vname(1,indxTime)(1:lvar)
     &                                , ininame(1:lstr)
     &                                , nf_strerror(ierr)
        goto 99                                           !--> ERROR
      endif

#ifdef XIOS
      start_date_run = tool_sectodat(start_time+xios_origin_date_in_sec)
      READ(start_date_run(1:4),fmt='(i4)') start_year
      READ(start_date_run(6:7),fmt='(i2)') start_month
      READ(start_date_run(9:10),fmt='(i2)') start_day
      READ(start_date_run(12:13),fmt='(i2)') start_hour
      READ(start_date_run(15:16),fmt='(i2)') start_minute
      READ(start_date_run(18:19),fmt='(i2)') start_second
#endif

      ierr=nf_inq_varid (ncid, 'time_step', varid)
      if (ierr .eq. nf_noerr) then
        start(1)=1
        start(2)=record
        count(1)=4
        count(2)=1
        ierr=nf_get_vara_int (ncid, varid, start, count, ibuff)
        if (ierr .eq. nf_noerr) then
          ntstart=ibuff(1)
          nrecrst=ibuff(2)
#if defined NO_RESTART_HIS
          nrechis=0
#else
          nrechis=ibuff(3)
#endif
#ifdef AVERAGES
          nrecavg=ibuff(4)
#endif

#if defined RESTART_DIAGS && ! defined XIOS

#ifdef DIAGNOSTICS_VRT
          nrecdiags_vrt=ibuff(3)
# ifdef AVERAGES
          nrecdiags_vrt_avg=ibuff(4)
# endif
#endif
#ifdef DIAGNOSTICS_EK
          nrecdiags_ek=ibuff(3)
# ifdef AVERAGES
          nrecdiags_ek_avg=ibuff(4)
# endif
#endif
#ifdef DIAGNOSTICS_PV
          nrecdiags_pv=ibuff(3)
# ifdef AVERAGES
          nrecdiags_pv_avg=ibuff(4)
# endif
#endif

#if defined DIAGNOSTICS_EDDY && ! defined XIOS
          nrecdiags_eddy=ibuff(3)
# ifdef AVERAGES
          nrecdiags_eddy_avg=ibuff(4)
# endif
#endif
#if defined DIAGNOSTICS_TS
          nrecdia=ibuff(3)
# ifdef AVERAGES
          nrecdia_avg=ibuff(4)
# endif
#endif
#ifdef DIAGNOSTICS_UV
          nrecdiaM=ibuff(3)
# ifdef AVERAGES
          nrecdiaM_avg=ibuff(4)
# endif
#endif
          
#endif
 

#ifdef EXACT_RESTART
          if (tindx.eq.2 .and. record.lt.max_rec) then
            start(2)=record+1
            ierr=nf_get_vara_int (ncid, varid, start, count, ibuff)
            if (ierr .eq. nf_noerr) then
              if (ibuff(1) .eq. ntstart+1) then
                forw_start=1
              else
          MPI_master_only write(stdout,'(/1x,3A,2I4/10x,4A/10x,A)')
     &           'WARNING: ',
     &           'Exact restart is requested, but is not possible: ',
     &           'records', record,record+1, 'in ''',ininame(1:lstr),
     &           ''' are not consecutive time steps ==> proceeding ',
     &                         'with', 'forward initial time step.'
              endif
            else
             write(*,2) 'time_step', record, ininame(1:lstr),
     &                                     nf_strerror(ierr)
              goto 99
            endif
          endif 
          if (tindx.eq.2 .and. forw_start.ne.1) return
          if (tindx.eq.1) then

            if (forw_start.eq.1) then
              MPI_master_only write(stdout,6)
     &             start_time*sec2day, record,
     &             ntstart, nrecrst, nrechis
            else
#endif
              MPI_master_only  write(stdout,7)
     &             start_time*sec2day, record,
     &             ntstart, nrecrst, nrechis
#ifdef EXACT_RESTART
            endif
          endif
     
#endif
        else
          MPI_master_only write(stdout,'(/1x,2A/)') 
     &         'GET_INITIAL ERROR: Cannot ',
     &         'read time and record indices.'
          goto 99                                         !--> ERROR
        endif
# if defined OUTPUTS_SURFACE && ! defined XIOS
        start(1)=1
        start(2)=record
        count(1)=6
        count(2)=1
        ierr=nf_get_vara_int (ncid, varid, start, count, ibuff)
        if (ierr .eq. nf_noerr) then
#if defined RESTART_DIAGS
          nrecsurf=ibuff(5)
#ifdef AVERAGES
          nrecsurf_avg=ibuff(6)
#endif
#else
          nrecsurf=0
#ifdef AVERAGES
          nrecsurf_avg=0
#endif
#endif
        endif
#endif
      else
        ntstart=1
        nrecrst=0
        nrechis=0
#ifdef AVERAGES
        nrecavg=0
#endif
        MPI_master_only write(stdout,'(/1x,2A,G12.4,1x,A,I4)') 
     &          'GET_INITIAL -- ',
     &       'Processing data for time =', start_time*sec2day,
     &                                      'record =', record
      endif
#ifdef EXACT_RESTART
      if (tindx.eq.1 .and. forw_start.eq.0) forw_start=ntstart
#endif


#ifdef AVERAGES

!      if (tindx.eq.1) then
!        if (ntstart.gt.1) then
!          MPI_master_only write(stdout,'(/4x,A,I5,A,I5/)')
!     &         'Adjusting averaging phase ntsavg =',  ntsavg,
!     &                          ' ---> ', ntsavg + ntstart-1
!          ntsavg=ntsavg + ntstart-1
!        endif
!      endif

#ifdef DIAGNOSTICS_UV
          ntsdiaM_avg = ntsavg
#endif
#ifdef DIAGNOSTICS_VRT
          ntsdiags_vrt_avg = ntsavg
#endif
#ifdef DIAGNOSTICS_EK
          ntsdiags_ek_avg = ntsavg
#endif
#ifdef DIAGNOSTICS_TS
          ntsdia_avg = ntsavg
#endif
#ifdef DIAGNOSTICS_PV
          ntsdiags_pv_avg = ntsavg
#endif
#if defined DIAGNOSTICS_EDDY && ! defined XIOS
          ntsdiags_eddy_avg = ntsavg
#endif
#endif
!     
          if (tindx.eq.1) then
            if (ntstart.lt.1) ntstart=1
            ntimes=ntstart+ntimes-1
          endif

!          write(*,*)'==============='
!          write(*,*)'FORW_START=',forw_start
!          write(*,*)'==============='
!          write(*,*)'==============='
!          write(*,*)'RECORD=',record
!          write(*,*)'==============='

!
!
! Free-surface.
!
      lvar=lenstr(vname(1,indxZ))
      ierr=nf_inq_varid (ncid, vname(1,indxZ)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (zeta(START_2D_ARRAY,1), ncid, varid,
     &                                         record, r2dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxZ)(1:lvar)
     &                                ,record,ininame(1:lstr)

          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) vname(1,indxZ)(1:lvar)
     &                                ,ininame(1:lstr)

        goto 99                                           !--> ERROR
      endif
!
! 2D momentum component in the XI-direction.
!
      lvar=lenstr(vname(1,indxUb))
      ierr=nf_inq_varid (ncid, vname(1,indxUb)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (ubar(START_2D_ARRAY,1), ncid, varid,
     &                                         record, u2dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxUb)(1:lvar)
     &                                  , record, ininame(1:lstr)
          goto 99                                        !--> ERROR
        endif
      else
        MPI_master_only write(stdout,1) vname(1,indxUb)(1:lvar)
     &                                ,ininame(1:lstr)
        goto 99                                           !--> ERROR
      endif
!
! 2D momentum component in the ETA-direction.
!
      lvar=lenstr(vname(1,indxVb))
      ierr=nf_inq_varid (ncid, vname(1,indxVb)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (vbar(START_2D_ARRAY,1), ncid, varid,
     &                                         record, v2dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxVb)(1:lvar)
     &                                , record, ininame(1:lstr)
          goto 99                                          !--> ERROR
        endif
      else
        MPI_master_only write(stdout,1) vname(1,indxVb)(1:lvar)
     &                                ,ininame(1:lstr)
        goto 99                                           !--> ERROR
      endif

#ifdef MORPHODYN
!
! Time evolving bathymetry (only needed for restart)
!
      lvar=lenstr(vname(1,indxHm))
      ierr=nf_inq_varid (ncid, vname(1,indxHm)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (h, ncid, varid, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxHm)(1:lvar)
     &                                   , record,ininame(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,3) vname(1,indxHm)(1:lvar)
     &                                          ,ininame(1:lstr)
      endif
#endif

#ifdef BBL
!
! Ripple height
!
      got_inibed(1)=.true.
      lvar=lenstr(vname(1,indxHrip))
      ierr=nf_inq_varid (ncid, vname(1,indxHrip)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (Hripple, ncid, varid, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxHrip)(1:lvar)
     &                                    , record, ininame(1:lstr)
          goto 99                                        !--> ERROR
        endif
      else
        MPI_master_only write(stdout,3) vname(1,indxHrip)(1:lvar)
     &                                , ininame(1:lstr)
        got_inibed(1)=.false.
      endif
!
! Ripple length
!
      got_inibed(2)=.true.
      lvar=lenstr(vname(1,indxLrip))
      ierr=nf_inq_varid (ncid, vname(1,indxLrip)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (Lripple, ncid, varid, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxLrip)(1:lvar)
     &                                    , record, ininame(1:lstr)
          goto 99                                        !--> ERROR
        endif
      else
        MPI_master_only write(stdout,3) vname(1,indxLrip)(1:lvar)
     &                                , ininame(1:lstr)
        got_inibed(2)=.false.
      endif
#endif /* BBL */

#ifdef SOLVE3D
!
! 3D momentum component in the XI-direction.
!
      lvar=lenstr(vname(1,indxU))
      ierr=nf_inq_varid (ncid, vname(1,indxU)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (u(START_2D_ARRAY,1,tindx), ncid, varid,
     &                                         record, u3dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxU)(1:lvar)
     &                                 , record, ininame(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only write(stdout,1) vname(1,indxU)(1:lvar)
     &                                ,ininame(1:lstr)
        goto 99                                           !--> ERROR
      endif
!
! 3D momentum component in the ETA-direction.
!
      lvar=lenstr(vname(1,indxV))
      ierr=nf_inq_varid (ncid, vname(1,indxV)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (v(START_2D_ARRAY,1,tindx), ncid, varid,
     &                                         record, v3dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxV)(1:lvar)
     &                                 , record, ininame(1:lstr)
          goto 99                                        !--> ERROR
        endif
      else
        MPI_master_only write(stdout,1) vname(1,indxV)(1:lvar)
     &                                ,ininame(1:lstr)
        goto 99                                           !--> ERROR
      endif
!
! Tracer variables.
!
#ifdef TRACERS 
      do itrc=1,NT
        got_tini(itrc)=.true.
        lvar=lenstr(vname(1,indxV+itrc))
        ierr=nf_inq_varid (ncid, vname(1,indxV+itrc)(1:lvar),varid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_fread (t(START_2D_ARRAY,1,tindx,itrc), ncid, varid,
     &                                                record, r3dvar)
          if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) vname(1,indxV+itrc)(1:lvar)
     &                                        , record, ininame(1:lstr)
            goto 99                                       !--> ERROR
          endif
        else
          MPI_master_only write(stdout,3) vname(1,indxV+itrc)(1:lvar),
     &                                    ininame(1:lstr)
          got_tini(itrc)=.false.
        endif
      enddo
#endif /* TRACERS */

# ifdef SEDIMENT
!
! Bed thickness
!
      got_inised(1)=.true.
      lvar=lenstr(vname(1,indxBTHK))
      ierr=nf_inq_varid (ncid, vname(1,indxBTHK)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (bed_thick, ncid, varid, record, b3dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxBTHK)(1:lvar)
     &                                    , record, ininame(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only write(stdout,3) vname(1,indxBTHK)(1:lvar)
     &     ,ininame(1:lstr)
        got_inised(1)=.false.
      endif
!
! Bed porosity
!
      got_inised(2)=.true.
      lvar=lenstr(vname(1,indxBPOR))
      ierr=nf_inq_varid (ncid, vname(1,indxBPOR)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (bed_poros, ncid, varid, record, b3dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxBPOR)(1:lvar)
     &                                    , record, ininame(1:lstr)
          goto 99                                        !--> ERROR
        endif
      else
        MPI_master_only write(stdout,3) vname(1,indxBPOR)(1:lvar)
     &       , ininame(1:lstr)
     &       , nf_strerror(ierr)
        got_inised(2)=.false.
      endif
!
! Bed sediment fractions
!
      got_inised(3)=.true.
      do itrc=1,NST
         indxWrk=indxBFRA(1)+itrc-1
         lvar=lenstr(vname(1,indxWrk))
      ierr=nf_inq_varid (ncid, vname(1,indxWrk)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (bed_frac(START_2D_ARRAY,1,itrc), ncid, varid,
     &                                               record, b3dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxWrk)(1:lvar)
     &                                   , record, ininame(1:lstr)
          goto 99                                        !--> ERROR
        endif
      else
        MPI_master_only write(stdout,3) vname(1,indxWrk)(1:lvar)
     &                                , ininame(1:lstr)
        got_inised(3)=.false.
      endif
      enddo
# endif /* SEDIMENT */

# ifdef LMD_MIXING
!
!  hbl
!
#  ifdef LMD_SKPP
#   ifdef LMD_SKPP2005
      lvar=lenstr(vname(1,indxHbl))
      ierr=nf_inq_varid (ncid, vname(1,indxHbl)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (hbls(START_2D_ARRAY,1), ncid, varid,
!#   else
!        ierr=nf_fread (hbl(START_2D_ARRAY), ncid, varid,
!#   endif
     &                                         1, r2dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxHbl)(1:lvar)
     &                                   , record, ininame(1:lstr)
          goto 99                                        !--> ERROR
        endif
      else
        MPI_master_only write(stdout,4) vname(1,indxHbl)(1:lvar)
     &                                , ininame(1:lstr)
      endif
#   endif /* LMD_SKPP2005 */
#  endif/* LMD_SKPP */

#  ifdef LMD_BKPP
#   ifdef LMD_BKPP2005
!
!  hbbl
!
      lvar=lenstr(vname(1,indxHbbl))
      ierr=nf_inq_varid (ncid, vname(1,indxHbbl)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (hbbl(START_2D_ARRAY), ncid, varid,
     &                                         1, r2dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxHbbl)(1:lvar)
     &                                    , record, ininame(1:lstr)
          goto 99                                        !--> ERROR
        endif
      else
        MPI_master_only write(stdout,4) vname(1,indxHbbl)(1:lvar)
     &                                , ininame(1:lstr)
      endif
      
#   endif /* LMD_BKPP2005 */
#  endif /* LMD_BKPP */
# endif /* LMD_MIXING */
      
# ifdef GLS_MIXING
!
!  turbulent kinetic energy.
!
        lvar=lenstr(vname(1,indxTke))
        ierr=nf_inq_varid (ncid, vname(1,indxTke)(1:lvar), varid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_fread (trb(START_2D_ARRAY,0,tindx,itke), ncid,
     &                                    varid, record, w3dvar)
          if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) vname(1,indxTke)(1:lvar)
     &                                     , record, ininame(1:lstr)
            goto 99                                       !--> ERROR
          endif
        else
          MPI_master_only write(stdout,4) vname(1,indxTke)(1:lvar)
     &                                  , ininame(1:lstr)          
      endif
      
!
!  generic length scale.
!
        lvar=lenstr(vname(1,indxGls))
        ierr=nf_inq_varid (ncid, vname(1,indxGls)(1:lvar), varid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_fread (trb(START_2D_ARRAY,0,tindx,igls), ncid,
     &                                    varid, record, w3dvar)
          if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) vname(1,indxGls)(1:lvar)
     &                                     , record, ininame(1:lstr)
            goto 99                                       !--> ERROR
          endif
        else 
          MPI_master_only write(stdout,4) vname(1,indxGls)(1:lvar)
     &                                  , ininame(1:lstr)          
      endif
      
!
!  vertical viscosity coefficient
!
        lvar=lenstr(vname(1,indxAkv))
        ierr=nf_inq_varid (ncid, vname(1,indxAkv)(1:lvar), varid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_fread (Akv(START_2D_ARRAY,0), ncid,
     &                                    varid, record, w3dvar)
          if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) vname(1,indxAkv)(1:lvar)
     &                                     , record, ininame(1:lstr)
            goto 99                                       !--> ERROR
          endif
        else
          MPI_master_only write(stdout,4) vname(1,indxAkv)(1:lvar)
     &                                  , ininame(1:lstr)
        endif
        
!
!  vertical diffusion coefficient for potential temperature.
!
#  if defined TEMPERATURE
        lvar=lenstr(vname(1,indxAkt))
        ierr=nf_inq_varid (ncid, vname(1,indxAkt)(1:lvar), varid)
       if (ierr .eq. nf_noerr) then
          ierr=nf_fread (Akt(START_2D_ARRAY,0,itemp), ncid,
     &                                    varid, record, w3dvar)
          if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) vname(1,indxAkt)(1:lvar)
     &                                     , record, ininame(1:lstr)
            goto 99                                       !--> ERROR
          endif
        else
          MPI_master_only write(stdout,4) vname(1,indxAkt)(1:lvar)
     &                                  , ininame(1:lstr)
        endif        
#  endif /* TEMPERATURE */

#  ifdef SALINITY
!
!  vertical diffusion coefficient for salinity.
!
        lvar=lenstr(vname(1,indxAks))
        ierr=nf_inq_varid (ncid, vname(1,indxAks)(1:lvar), varid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_fread (Akt(START_2D_ARRAY,0,isalt), ncid,
     &                                    varid, record, w3dvar)
          if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) vname(1,indxAks)(1:lvar)
     &                                     , record, ininame(1:lstr)
            goto 99                                       !--> ERROR
          endif
        else
          MPI_master_only write(stdout,4) vname(1,indxAks)(1:lvar)
     &                                  , ininame(1:lstr)
        endif
        
#  endif /* SALINITY */
# endif  /* GLS_MIXING */

# if defined LMD_MIXING
#  if defined M3FAST
!
!  bustr
!
      lvar=lenstr(vname(1,indxBustr))
      ierr=nf_inq_varid (ncid, vname(1,indxBustr)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
      ierr=nf_fread (bustr(START_2D_ARRAY), ncid,
     &                                     varid, record, u2dvar)
        if (ierr .ne. nf_noerr) then
         MPI_master_only write(stdout,2) vname(1,indxBustr)(1:lvar)
     &                                    , record, ininame(1:lstr)
          goto 99                                        !--> ERROR
        endif
      else
        MPI_master_only write(stdout,4) vname(1,indxBustr)(1:lvar)
     &                                , ininame(1:lstr)
      endif
!  
!  bvstr
!
       lvar=lenstr(vname(1,indxBvstr))
       ierr=nf_inq_varid (ncid, vname(1,indxBvstr)(1:lvar), varid)
       if (ierr .eq. nf_noerr) then
       ierr=nf_fread (bvstr(START_2D_ARRAY), ncid,
     &                                     varid, record, v2dvar)
         if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxBvstr)(1:lvar)
     &                                     , record, ininame(1:lstr)
           goto 99                                        !--> ERROR
         endif
       else
         MPI_master_only write(stdout,4)  vname(1,indxBvstr)(1:lvar)
     &                                  , ininame(1:lstr)
        endif
#  endif        
# endif
                
# if defined EXACT_RESTART
!
!  rufrc
!
      lvar=lenstr(vname(1,indxrufrc))
      ierr=nf_inq_varid (ncid, vname(1,indxrufrc)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
      ierr=nf_fread (rufrc_bak(START_2D_ARRAY,3-tindx), ncid,
     &                                     varid, record, u2dvar)
        if (ierr .ne. nf_noerr) then
         MPI_master_only write(stdout,2) vname(1,indxrufrc)(1:lvar)
     &                                    , record, ininame(1:lstr)
          goto 99                                        !--> ERROR
        endif
      else
        MPI_master_only write(stdout,5) vname(1,indxrufrc)(1:lvar)
     &                                , ininame(1:lstr)
      endif
!  
!  rvfrc
!
       lvar=lenstr(vname(1,indxrvfrc))
       ierr=nf_inq_varid (ncid, vname(1,indxrvfrc)(1:lvar), varid)
       if (ierr .eq. nf_noerr) then
       ierr=nf_fread (rvfrc_bak(START_2D_ARRAY,3-tindx), ncid,
     &                                     varid, record, v2dvar)
         if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxrvfrc)(1:lvar)
     &                                     , record, ininame(1:lstr)
           goto 99                                        !--> ERROR
         endif
       else
         MPI_master_only write(stdout,5)  vname(1,indxrvfrc)(1:lvar)
     &                                  , ininame(1:lstr)
        endif
       
#  ifdef M3FAST
!
! ru_nbq.
!
      lvar=lenstr(vname(1,indxru_nbq))
      ierr=nf_inq_varid (ncid,
     &                  vname(1,indxru_nbq)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (ru_nbq(START_2D_ARRAY,1), ncid, varid,
     &                                         record, u3dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2)
     &                           vname(1,indxru_nbq)(1:lvar)
     &                         , record, ininame(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only write(stdout,5) vname(1,indxru_nbq)(1:lvar)
     &                                 ,ininame(1:lstr)
      endif
!
! rv_nbq.
!
      lvar=lenstr(vname(1,indxrv_nbq))
      ierr=nf_inq_varid (ncid,
     &                  vname(1,indxrv_nbq)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (rv_nbq(START_2D_ARRAY,1), ncid, varid,
     &                                         record, v3dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2)
     &                           vname(1,indxrv_nbq)(1:lvar)
     &                             , record, ininame(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only write(stdout,5) vname(1,indxrv_nbq)(1:lvar)
     &                                 ,ininame(1:lstr)
      endif
!
! ru_nbq_avg2.
!
      lvar=lenstr(vname(1,indxru_nbq_avg2))
      ierr=nf_inq_varid (ncid,
     &                  vname(1,indxru_nbq_avg2)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (ru_nbq_avg2(START_2D_ARRAY,1), ncid, varid,
     &                                         record, u3dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2)
     &                           vname(1,indxru_nbq_avg2)(1:lvar)
     &                                  , record, ininame(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only write(stdout,5) vname(1,indxru_nbq_avg2)(1:lvar)
     &                                 ,ininame(1:lstr)
      endif
!
! rv_nbq_avg2.
!
      lvar=lenstr(vname(1,indxrv_nbq_avg2))
      ierr=nf_inq_varid (ncid,
     &                  vname(1,indxrv_nbq_avg2)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (rv_nbq_avg2(START_2D_ARRAY,1), ncid, varid,
     &                                         record, v3dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2)
     &                           vname(1,indxrv_nbq_avg2)(1:lvar)
     &                         , record, ininame(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only write(stdout,5) vname(1,indxrv_nbq_avg2)(1:lvar)
     &                                 ,ininame(1:lstr)
      endif
!
! qdmu_nbq.
!
      lvar=lenstr(vname(1,indxqdmu_nbq))
      ierr=nf_inq_varid (ncid,
     &                  vname(1,indxqdmu_nbq)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (qdmu_nbq(START_2D_ARRAY,1), ncid, varid,
     &                                         record, u3dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2)
     &                           vname(1,indxqdmu_nbq)(1:lvar)
     &                         , record, ininame(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only write(stdout,5) vname(1,indxqdmu_nbq)(1:lvar)
     &                                 ,ininame(1:lstr)
      endif
!
! qdmv_nbq.
!
      lvar=lenstr(vname(1,indxqdmv_nbq))
      ierr=nf_inq_varid (ncid,
     &                  vname(1,indxqdmv_nbq)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (qdmv_nbq(START_2D_ARRAY,1), ncid, varid,
     &                                         record, v3dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2)
     &                           vname(1,indxqdmv_nbq)(1:lvar)
     &                         , record, ininame(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only write(stdout,5) vname(1,indxqdmv_nbq)(1:lvar)
     &                                 ,ininame(1:lstr)
      endif
      
#  endif  /* M3FAST */

#  ifdef TS_MIX_ISO_FILT
!
! density gradients for use in t3dmix
!
      lvar=lenstr(vname(1,indxdRdx))
      ierr=nf_inq_varid (ncid,
     &                  vname(1,indxdRdx)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (dRdx(START_2D_ARRAY,1), ncid, varid,
     &                                         record, u3dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2)
     &                           vname(1,indxdRdx)(1:lvar)
     &                           , record, ininame(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only write(stdout,5) vname(1,indxdRdx)(1:lvar)
     &                                 ,ininame(1:lstr)
      endif
!
!
      lvar=lenstr(vname(1,indxdRde))
      ierr=nf_inq_varid (ncid,
     &                  vname(1,indxdRde)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (dRde(START_2D_ARRAY,1), ncid, varid,
     &                                         record, v3dvar)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxdRde)(1:lvar)
     &                                    , record, ininame(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only write(stdout,5) vname(1,indxdRde)(1:lvar)
     &                                 ,ininame(1:lstr)
      endif
#  endif /* TS_MIX_ISO_FILT */

# endif  /* EXACT_RESTART */

!
!  Close input NetCDF file.
!
      ierr=nf_close(ncid)
#endif /* SOLVE3D */

#if defined OA_COUPLING || defined OW_COUPLING
#  if defined EXACT_RESTART
!     Do cpl_prism_define just on time. 
!     In case of EXACT_RESTART add a condition on tindx value
      if (tindx.eq.1) then
#  endif
        call cpl_prism_define
        oasis_time = 0
        MPI_master_only write(*,*)'CPL-ROMS: OASIS_TIME',oasis_time
#  if defined EXACT_RESTART   
      endif
#  endif
#endif

  1   format(/1x,'GET_INITIAL - unable to find variable:',    1x,A,
     &                            /15x,'in input NetCDF file:',1x,A/)
  2   format(/1x,'GET_INITIAL - error while reading variable:',1x, A,
     &    2x,'at time record =',i4/15x,'in input NetCDF file:',1x,A/)
  3   format(/1x,'GET_INITIAL - unable to find variable:',    1x,A,
     &                            /15x,'in input NetCDF file:',1x,A,
     &    1x,'-> analytical value'/)
  4   format(/1x,'GET_INITIAL -  unable to find variable:',1x,A,
     &    /15x' in input NetCDF file:',1x,A,
     &     1x,' ==> Initialized to zero state.',
     &    /18x' >> CAUTION in case of #define EXACT_RESTART << ',
     &    /21x ' If it is the case ', 
     &    /21x'    - OK if it is a ''cold start''',
     &        ' i.e coming from a 3rd-party initial file',
     &    /21x'    - otherwise if it is a ''hot start'' ',
     &         'i.e from a restart file produced by this code: ',
     &    /21x'      => problem: run is not restartable',
     &    /21x'      => check your initial file')            
  5   format(/1x,'GET_INITIAL -  unable to find variable:',1x,A,
     &    /15x' in input NetCDF file:',1x,A,
     &     1x,' ==> Initialized to zero state.',
     &    /18x' >> EXACT_RESTART is defined <<',
     &    /21x'    - OK if it is a ''cold start'' ',
     &         'i.e from a 3rd-party initial file',
     &    /21x'    - otherwise, if it a ''hot start'' ',
     &         'i.e from a restart file produced by this code ',
     &    /21x'      => problem: run is not restartable ',
     &    /21x'      => check your initial file')      
 6    format(/6x,'GET_INITIAL: Exact restart from day =',F12.4,
     &     '  rec =',I4,'(',I8,',',I4,','I4,').')
 7    format(/6x,'GET_INITIAL: Restarted from day =',F12.4,
     &     ' rec=',I4,'(',I8,',',I4,','I4').')
      
      return
  99  may_day_flag=2
      return
      end
