!     $Id: get_btflux.F gcambon $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers Universit (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"
#if defined SOLVE3D && ( defined BHFLUX || defined BWFLUX )
      subroutine get_btflux (itrc)
!
!  Read in point or grided bottom hydrothermal heat flux at the appropriate
!  time from forcing NetCDF file.
!
!
# define BTFLUX_DATA
      implicit none
# include "param.h"
# include "forces.h"
# include "scalars.h"
# include "ncscrum.h"
# include "netcdf.inc"
      real cff
      integer itrc, i,ierr, lstr,lvar,lenstr, nf_fread, advance_cycle
!
! Initialization: Inquire about the contents of forcing NetCDF file:
!================  variables and dimensions. Check for consistency.
!
      if (may_day_flag.ne.0) return      !-->  EXIT
      if (itbtf(itrc).eq.0 .or. iic.eq.0) then
        lstr=lenstr(btfname)
c**        call opencdf (btfname,N)
        if (may_day_flag.ne.0) return    !-->  EXIT
!
! If not opened yet, open forcing NetCDF file for reading.
! Find and save IDs for relevant variables, determine whether
! surface flux for tracer itrc is a field or scalar value.
!
        if (ncidbtf.eq.-1) then
          ierr=nf_open (btfname(1:lstr), nf_nowrite, ncidbtf)
          if (ierr .ne. nf_noerr) goto 4
        endif

        if (itrc.eq.itemp) then
          ierr=nf_inq_varid (ncidbtf, 'bhf_time', btf_tid(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,3) 'bhf_time', btfname(1:lstr)
            goto 99                                       !--> ERROR
          endif
# ifdef SALINITY
        elseif (itrc.eq.isalt) then
          ierr=nf_inq_varid (ncidbtf, 'bwf_time', btf_tid(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,3) 'bwf_time', btfname(1:lstr)
            goto 99                                       !--> ERROR
          endif
# endif
        endif
        lvar=lenstr(vname(1,indxBhflx+itrc-1))
        ierr=nf_inq_varid (ncidbtf, vname(1,indxBhflx+itrc-1),
     &                                          btf_id(itrc))
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidbtf, btf_id(itrc), i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lbtfgrd(itrc)=1
            else
              lbtfgrd(itrc)=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxBhflx+itrc-1), btfname(1:lstr)
          goto 99                                         !--> ERROR
        endif

!
! Determine whether there is cycling to reuse the input data and
! find cycling period "btf_cycle", set initial cycling index
! "btf_ncycle" and record index "btf_rec".
! Set initial value for time index "itbtf" and both time record
! bounds to large negative artificial values, so that it will
! trigger the logic in reading part below.
!
        call set_cycle (ncidbtf,  btf_tid(itrc),    ntbtf(itrc),
     &           btf_cycle(itrc), btf_ncycle(itrc), btf_rec(itrc))
        if (may_day_flag.ne.0) return    !-->  EXIT

        itbtf(itrc)=2
        btf_time(1,itrc)=-1.E+20
        btf_time(2,itrc)=-1.E+20
!
! Set scale factors to convert input fluxes to model units:
!
!     Heat  flux - convert from [Watts/m^2] to [Celsius m/s]
!     Fresh Water flux - convert from [cm/day] to [PSU m/s]
!
        btf_scale(itemp)=1./(rho0*Cp)
# ifdef SALINITY
        btf_scale(isalt)=0.01/86400.
# endif
      endif
!
! Reading data from the forcing file:  Get out, if model time is
!======== ==== ==== === ======= =====  already within the interval
! set by the past and future data times. Otherwise flip the time
! index, increment record and cyclin indices and read a new portion
! of data. Repeat it until model time is between the two times from
! data.
!
  1   i=3-itbtf(itrc)
       cff=time+0.5*dt
        if (btf_time(i,itrc).le.cff .and.
     &      cff.lt.btf_time(itbtf(itrc),itrc)) return
        ierr=advance_cycle (btf_cycle(itrc),  ntbtf(itrc),
     &                      btf_ncycle(itrc), btf_rec(itrc))
        if (ierr .ne. 0) goto 7                           !--> ERROR

        ierr=nf_get_var1_FTYPE(ncidbtf, btf_tid(itrc),
     &                             btf_rec(itrc), cff)
# ifdef USE_CALENDAR
        call tool_origindate(ncidbtf,btf_tid,
     &                       btf_origin_date_in_sec)
        cff=cff+btf_origin_date_in_sec*sec2day
# endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'btf_time', btf_rec(i)
          goto 99                                         !--> ERROR
        endif
        btf_time(i,itrc)=cff*day2sec+btf_cycle(itrc)*btf_ncycle(itrc)
        if (btf_time(itbtf(itrc),itrc).eq.-1.E+20)
     &      btf_time(itbtf(itrc),itrc)=btf_time(i,itrc)

        if (lbtfgrd(itrc).eq.1) then
          ierr=nf_fread(btflxg(START_2D_ARRAY,i,itrc), ncidbtf,
     &                     btf_id(itrc), btf_rec(itrc), r2dvar)
        else
          ierr=nf_get_var1_FTYPE(ncidbtf, btf_id(itrc),
     &                     btf_rec(itrc), btflxp(i,itrc))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'btflux', btf_rec 
          goto 99                                         !--> ERROR
        endif
       itbtf(itrc)=i
       MPI_master_only write(stdout,'(6x,A,1x,I2,1x,A,1x,g12.4,1x,I4)')
     &            'GET_BTFLUX -- Read surface flux of tracer', itrc,
     &                                             'for time =', cff
# ifdef USE_CALENDAR
     &  -btf_origin_date_in_sec*sec2day
# endif
# ifdef MPI
     &                                                      , mynode
#endif

      if (ntbtf(itrc).gt.1) goto 1
      if (ntbtf(itrc).eq.1) return

!
! Sort out error messages: The following portion of the code is
!===== === ===== ========= not accessed unless something goes wrong.
!
  3   format(/,' GET_BTFLUX - ERROR: unable to find forcing ',
     &              'variable: ',a,/,14x,'in forcing file: ',a)
  4   write(stdout,5) btfname(1:lstr)
  5   format(/,' GET_BTFLUX - ERROR: unable to open file: ', a)
      goto 99
  6   format(/,' GET_BTFLUX - ERROR while reading variable: ',
     &                                a,' at TIME index = ',i4)
  7   write(stdout,8) btf_rec(itrc),ntbtf(itrc),btfname(1:lstr),
     &                tdays, btf_time(itbtf(itrc),itrc)*sec2day
# ifdef USE_CALENDAR
     &  -btf_origin_date_in_sec*sec2day
# endif
  8   format(/,' GET_BTFLUX - ERROR: requested time record ',I4,
     &       1x,'exeeds the last available', /,14x,'record ',I4,
     &       1x,'in forcing file: ', a,/,14x,'TDAYS = ',
     &             g12.4,2x,'last available BTF_TIME = ',g12.4)
  99  may_day_flag=2
      return
      end

      subroutine set_btflux_tile (Istr,Iend,Jstr,Jend, itrc)
!
!--------------------------------------------------------------------
! Set-up surface tracer flux for current tile.
!--------------------------------------------------------------------
!
# define BTFLUX_DATA
      implicit none
      integer Istr,Iend,Jstr,Jend, itrc, i,j, it1,it2
      real cff, cff1, cff2
# include "param.h"
# include "forces.h"
# include "scalars.h"
!
# include "compute_extended_bounds.h"
!
!      write(*,*)'Enter set_btflux'
      it1=3-itbtf(itrc)
      it2=itbtf(itrc)
      cff=time+0.5*dt
      cff1=btf_time(it2,itrc)-cff
      cff2=cff-btf_time(it1,itrc)
!
! Load time invariant surface tracer flux.
! Time interpolation is not performed in this case.
!
      if (btf_cycle(itrc).lt.0.) then
        if (FIRST_RST_TIME_STEP) then
          if (lbtfgrd(itrc).eq.1) then
            cff=btf_scale(itrc)
            do j=JstrR,JendR
              do i=IstrR,IendR
                btflx(i,j,itrc)=cff*btflxg(i,j,itbtf(itrc),itrc)
!                if ( i==10 .and. j==10) then
!                  write(*,*)'I=',i,' J=',j
!                  write(*,*)'btflx(10,10,1)=',btflx(10,10,1)
!                endif
              enddo
            enddo
          else
            cff=btf_scale(itrc)*btflxp(itbtf(itrc),itrc)
            do j=JstrR,JendR
              do i=IstrR,IendR
                btflx(i,j,itrc)=cff
              enddo
            enddo
          endif
        endif
!
! Time-interpolate surface tracer flux from grided or point data.
! Check that for the next time step [when time=time+dt] time+dt is
! still between srf_tintrp(it1) and srf_tintrp(it2); and if not,
! set synchro_flag top signal that the new forcing data should be
! read from the netCDF input file (master thread only).
!
      elseif (cff1.ge.0. .and. cff2.ge.0.) then
        if (ZEROTH_TILE .and. cff1.lt.dt) synchro_flag=.TRUE.
        cff=btf_scale(itrc)/(cff1+cff2)
        cff1=cff1*cff
        cff2=cff2*cff
        if (lbtfgrd(itrc).eq.1) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              btflx(i,j,itrc)=cff1*btflxg(i,j,it1,itrc)
     &                       +cff2*btflxg(i,j,it2,itrc)
            enddo
          enddo
        else
          cff=cff1*btflxp(it1,itrc)+cff2*btflxp(it2,itrc)
          do j=JstrR,JendR
            do i=IstrR,IendR
              btflx(i,j,itrc)=cff
            enddo
          enddo
        endif
        !if ( i==10 .and. j==10) then
        !  write(*,*)'I=',i,' J=',j
        !  write(*,*)'btflx(10,10,1)=',btflx(10,10,1)
        !endif
!     
!     Unable to set-up surface tracer flux:
!     Complain about the error and signal to quit (ONE THRED ONLY).
!
      else
        if (ZEROTH_TILE) then
          write(stdout,1) 'btf_time',tdays,btf_time(it2,itrc)*sec2day
# ifdef USE_CALENDAR
     &  -btf_origin_date_in_sec*sec2day
# endif
  1       format(/,' SET_BTFLUX_TILE - current model time exceeds',
     &          1x,'ending value for variable: ',A8,/,14x,'TDAYS = ',
     &                                    g12.4, 2x, 'TBTF = ',g12.4)
          may_day_flag=2
        endif
      endif
      return
      end
#else
      subroutine get_btflux_empty
      return
      end
#endif /* defined SOLVE3D && (!defined ANA_BTFLUX || !defined ANA_BSFLUX) */
