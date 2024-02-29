! $Id: get_wwave.F 1458 2014-02-03 15:01:25Z gcambon $
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
#if defined  WAVE_OFFLINE && (defined MUSTANG || defined BBL || defined MRL_WCI)

                              ! Read point or gridded wind wave
      subroutine get_wwave    ! height, direction, period at approp.
                              ! time from forcing NetCDF file.
      implicit none
# include "param.h" 
# include "scalars.h"
# include "ncscrum.h"
# include "forces.h"
      real cff
      integer i,ierr, lstr,lvar, lenstr, nf_fread, advance_cycle
# include "netcdf.inc"
!
! Initialization: Inquire about the contents of forcing NetCDF file:
!================  variables and dimensions. Check for consistency.
!


# ifndef MUSTANG
      ncidwave=ncidfrc
      wave_file=frcname
# endif

      if (may_day_flag.ne.0) return      !-->  EXIT
      if (itww.eq.0 .or. iic.eq.0) then
        lstr=lenstr(wave_file)
        if (may_day_flag.ne.0) return    !-->  EXIT
!
! If not opened yet, open forcing NetCDF file for reading.
! Find and save IDs for relevant variables, determine whether
! wind wave data are field or scalar values. 
!

        if (ncidwave.eq.-1) then
          ierr=nf_open(wave_file(1:lstr), nf_nowrite, ncidwave)
          if (ierr .ne. nf_noerr) goto 4                !--> ERROR 
        endif

        ierr=nf_inq_varid (ncidwave, 'wwv_time', ww_tid)
        if (ierr .ne. nf_noerr) then
          write(stdout,3) 'wwv_time', wave_file(1:lstr)
          goto 99                                         !--> ERROR
        endif

# ifndef MUSTANG
        ncidfrc=ncidwave
# endif

! Check if wave amplitude present, and if so, if field or scalar series
      
        lvar=lenstr(vname(1,indxWWA))
        ierr=nf_inq_varid (ncidwave,vname(1,indxWWA)(1:lvar),wwa_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidwave, wwa_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              wwagrd=1.
            else
              wwagrd=0.
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxWWA)(1:lvar),wave_file(1:lstr)
          goto 99                                         !--> ERROR
        endif

! Check if wave direction present, and if so, if field or scalar series

        lvar=lenstr(vname(1,indxWWD))
        ierr=nf_inq_varid (ncidwave,vname(1,indxWWD)(1:lvar),wwd_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidwave, wwd_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              wwdgrd=1.
            else
              wwdgrd=0.
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxWWD)(1:lvar),wave_file(1:lstr)
          goto 99                                         !--> ERROR
        endif

     
! Check if wave period present, and if so, if field or scalar series
      
        lvar=lenstr(vname(1,indxWWP))
        ierr=nf_inq_varid (ncidwave,vname(1,indxWWP)(1:lvar),wwp_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidwave, wwp_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              wwpgrd=1.
            else
              wwpgrd=0.
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxWWP)(1:lvar),wave_file(1:lstr)
          goto 99                                         !--> ERROR
        endif

# if defined MUSTANG
! Check if wave orbital velocity present, and if so, if field or scalar series

        lvar=lenstr(vname(1,indxWWU))
        ierr=nf_inq_varid (ncidwave,vname(1,indxWWU)(1:lvar),wwu_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidwave, wwu_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              wwugrd=1.
            else
              wwugrd=0.
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxWWU)(1:lvar),wave_file(1:lstr)
          goto 99                                         !--> ERROR
        endif

# endif

! Check if breaking dissipation present, and if so, if field or scalar series

# if defined MRL_WCI && defined WAVE_OFFLINE_BREAKING   
        lvar=lenstr(vname(1,indxWEB))
        ierr=nf_inq_varid (ncidwave,vname(1,indxWEB)(1:lvar),wweb_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidwave, wweb_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              wwebgrd=1.
            else
              wwebgrd=0.
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxWEB)(1:lvar),wave_file(1:lstr)
          goto 99                                         !--> ERROR
        endif
# endif

! Check if frictional dissipation present, and if so, if field or scalar series

# if defined MRL_WCI && defined WAVE_OFFLINE_FRICTION  
        lvar=lenstr(vname(1,indxWED))
        ierr=nf_inq_varid (ncidwave,vname(1,indxWED)(1:lvar),wwed_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidwave, wwed_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              wwedgrd=1.
            else
              wwedgrd=0.
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxWED)(1:lvar),wave_file(1:lstr)
          goto 99                                         !--> ERROR
        endif
# endif

! Check if roller dissipation present, and if so, if field or scalar series

# if defined MRL_WCI && defined WAVE_OFFLINE_ROLLER
        lvar=lenstr(vname(1,indxWER))
        ierr=nf_inq_varid (ncidwave,vname(1,indxWER)(1:lvar),wwer_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidwave, wwer_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              wwergrd=1.
            else
              wwergrd=0.
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxWER)(1:lvar),wave_file(1:lstr)
          goto 99                                         !--> ERROR
        endif
# endif

!
! Determine whether cycling is defined to reuse the input data
! and find cycling period "ww_cycle", set initial cycling
! index "ww_ncycle" and record index "ww_rec".
! Set initial value for time index "itww" and both time record
! bounds to large negative artificial values, so that it will
! trigger the logic in the reading part below.
! Also set scale factor to convert input data to model units:
! convert wave direction from degrees to radians;
! no conversion for amplitude and period.
!
        call set_cycle (ncidwave,   ww_tid,    ntww,
     &                  ww_cycle, ww_ncycle, ww_rec)
        if (may_day_flag.ne.0) return    !-->  EXIT
      
        itww=2
        wwv_time(1)=-1.E+20
        wwv_time(2)=-1.E+20


        wwa_scale=1.
        wwd_scale=deg2rad
        wwp_scale=1.
        wweb_scale=1.
        wwed_scale=1.
        wwer_scale=1.
#  ifdef MUSTANG
        wwu_scale=1.
#  endif
     
      endif                   !<-- itww.eq.0 .or. iic.eq.0

!
! Reading data from the forcing file:  Get out, if model time is
!======== ==== ==== === ======= =====  already within the interval
! set by the past and future data times. Otherwise flip the time
! index, increment record and cyclin indices and read a new portion
! of data. Repeat it until model time is between the two times from
! data.
!
  1   i=3-itww
      cff=time+0.5*dt
      if (wwv_time(i).le.cff .and. cff.lt.wwv_time(itww)) return
      ierr=advance_cycle (ww_cycle, ntww, ww_ncycle, ww_rec)
      if (ierr.ne.0) goto 7                           !--> ERROR

      ierr=nf_get_var1_FTYPE(ncidwave, ww_tid, ww_rec, cff)
      if (ierr .ne. nf_noerr) then
        write(stdout,6) 'wwv_time', ww_rec 
        goto 99                                         !--> ERROR
      endif
# ifdef USE_CALENDAR
        call tool_origindate(ncidwave,ww_tid,
     &                       ww_origin_date_in_sec)
        cff=cff+ww_origin_date_in_sec*sec2day
# endif
      wwv_time(i)=cff*day2sec+ww_cycle*ww_ncycle
      if (wwv_time(itww).eq.-1.E+20) wwv_time(itww)=wwv_time(i)

! Read wave amplitude

      if (wwagrd.eq.1.) then
        ierr=nf_fread (wwag(START_2D_ARRAY,i), ncidwave, wwa_id,
     &                                           ww_rec, r2dvar)
      else
        ierr=nf_get_var1_FTYPE (ncidwave,wwa_id,ww_rec,wwap(i))
      endif

      if (ierr .ne. nf_noerr) then
        write(stdout,6) 'wind-wave amplitude', ww_rec
        goto 99                                         !--> ERROR
      endif 
      
! Read wave direction

      if (wwdgrd.eq.1.) then
        ierr=nf_fread (wwdg(START_2D_ARRAY,i), ncidwave, wwd_id,
     &                                           ww_rec, r2dvar)
      else
        ierr=nf_get_var1_FTYPE (ncidwave,wwd_id,ww_rec,wwdp(i))
      endif

      if (ierr .ne. nf_noerr) then
        write(stdout,6) 'wind-wave direction', ww_rec
        goto 99                                         !--> ERROR
      endif 
      
! Read wave period

      if (wwpgrd.eq.1.) then
        ierr=nf_fread (wwpg(START_2D_ARRAY,i), ncidwave, wwp_id,
     &                                           ww_rec, r2dvar)
      else
        ierr=nf_get_var1_FTYPE (ncidwave,wwp_id,ww_rec,wwpp(i))
      endif

      if (ierr .ne. nf_noerr) then
        write(stdout,6) 'wind-wave period', ww_rec
        goto 99                                         !--> ERROR
      endif 
      
! Read orbital velocity

# ifdef MUSTANG
      if (wwugrd.eq.1.) then
        ierr=nf_fread (wwug(START_2D_ARRAY,i), ncidwave, wwu_id,
     &                                           ww_rec, r2dvar)
      else
        ierr=nf_get_var1_FTYPE (ncidwave,wwu_id,ww_rec,wwup(i))
      endif

      if (ierr .ne. nf_noerr) then
        write(stdout,6) 'wave orbital velocity', ww_rec
        goto 99                                         !--> ERROR
      endif
# endif

! Read breaking dissipation

# if defined MRL_WCI && defined WAVE_OFFLINE_BREAKING
      if (wwebgrd.eq.1.) then
        ierr=nf_fread (wweb(START_2D_ARRAY,i), ncidwave, wweb_id,
     &                                           ww_rec, r2dvar)
      else
        ierr=nf_get_var1_FTYPE (ncidwave,wweb_id,ww_rec,wwebp(i))
      endif

      if (ierr .ne. nf_noerr) then
        write(stdout,6) 'wave breaking dissipation', ww_rec
        goto 99                                         !--> ERROR
      endif       
# endif
    
! Read frictional dissipation

# if defined MRL_WCI && defined WAVE_OFFLINE_FRICTION
      if (wwedgrd.eq.1.) then
        ierr=nf_fread (wwedd(START_2D_ARRAY,i), ncidwave, wwed_id,
     &                                           ww_rec, r2dvar)
      else
        ierr=nf_get_var1_FTYPE (ncidwave,wwed_id,ww_rec,wwedp(i))
      endif

      if (ierr .ne. nf_noerr) then
        write(stdout,6) 'wave frictional dissipation', ww_rec
        goto 99                                         !--> ERROR
      endif       
# endif
   
! Read frictional dissipation

# if defined MRL_WCI && defined WAVE_OFFLINE_ROLLER 
      if (wwergrd.eq.1.) then
        ierr=nf_fread (wwer(START_2D_ARRAY,i), ncidwave, wwer_id,
     &                                           ww_rec, r2dvar)
      else
        ierr=nf_get_var1_FTYPE (ncidwave,wwer_id,ww_rec,wwerp(i))
      endif

      if (ierr .ne. nf_noerr) then
        write(stdout,6) 'wave roller dissipation', ww_rec
        goto 99                                         !--> ERROR
      endif       
# endif

  
      itww=i
      MPI_master_only write(stdout,'(6x,A,1x,A,1x,F10.4,1x,I4)')
     &  'GET_WWAVE --','Read wind wave input data for time =', cff
# ifdef USE_CALENDAR
     &  -ww_origin_date_in_sec*sec2day
# endif
# ifdef MPI
     &                                                      , mynode
# endif
      if (ntww.gt.1) goto 1
      if (ntww.eq.1) return 
!
!
! Sort out error messages: The following portion of the code is
!===== === ===== ========= not accessed unless something goes wrong.
!
  3   format(/,' GET_WWAVE - unable to find forcing variable: ',a,
     &                          /,15x,'in forcing NetCDF file: ',a)
  4   write(stdout,5) wave_file(1:lstr)
  5   format(/,' GET_WWAVE - unable to open forcing NetCDF file:',
     &                                                        1x,a)
      goto 99

  6   format(/,' GET_WWAVE - error while reading variable: ',a,2x,
     &                                      ' at TIME index = ',i4)
  7   write(stdout,8) ww_rec, ntww, wave_file(1:lstr), tdays,
     &                                wwv_time(itww)*sec2day
# ifdef USE_CALENDAR
     &  -ww_origin_date_in_sec*sec2day
# endif
  8   format(/,' GET_WWAVE - ERROR: requested time record ',I4,
     &       1x,'exceeds the last available', /,14x,'record ',I4,
     &       1x,'in forcing NetCDF file: ',a, /,14x,'TDAYS = ',
     &       g12.4,2x,'last available wwv_time = ',g12.4)

  99  may_day_flag=2
    
      return
      end
!
!=====================================================================
!
!                             SET_WWAVE
!
!=====================================================================
!
      subroutine set_wwave (tile)
      implicit none
# include "param.h"
      integer tile 
# include "compute_tile_bounds.h"
!
      call set_wwave_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine set_wwave_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
!  Set-up wave amplitude, direction, period for current tile.
!--------------------------------------------------------------------
!
      implicit none

      integer Istr,Iend,Jstr,Jend, i,j,k, it1,it2
      real    fac,fac1,fac2,fac11,fac12,fac21,fac22,fac31,fac32,
     &                      fac41,fac42,fac51,fac52,fac61,fac62,
     &        wwA, wwD, wwP, wwe, Dstp,cfrq,hrms
# ifdef MUSTANG
     &       , wwU, fac71, fac72
# endif
      real    cff,cff1,cff2, kh,khd,kw, eps, wramp,m2p,isr2 
# ifdef BBL_OFFLINE
     &                                            , ab, orb_vel
# endif
      parameter (eps = 1.0d-10, m2p=1.2804973111, isr2=0.70710678)
!
# include "param.h"
# include "grid.h"
# include "forces.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ocean3d.h"
!
# include "compute_extended_bounds.h"
!
! Set coefficients for interpolation. Check that for the next time
! step [when time=time+dt] both weights will still be positive, and
! if not, set synchro_flag to signal that new data should be read
! from an appropriate netCDF input file (master thread only).
! After that, either load time-invariant data, or interpolate in time.
! Complain about error and signal to quit, if interpolation is
! needed, but not possible.
!
!
      it1=3-itww
      it2=itww
      fac=time+0.5*dt
      fac1=wwv_time(it2)-fac
      fac2=fac-wwv_time(it1)
!
! Load time-invariant wind-wave data.
! Time interpolation is not performed in this case.
!
      if (ww_cycle.lt.0.) then
        if (FIRST_TIME_STEP) then
        
            wwA=wwa_scale*wwap(itww)
            wwD=wwd_scale*wwdp(itww)
            wwP=wwp_scale*wwpp(itww)
# ifdef MUSTANG
            wwU=wwu_scale*wwup(itww)
# endif
# if defined MRL_WCI && defined WAVE_OFFLINE_BREAKING    
            wwEpb=wweb_scale*wwebp(itww)
# endif
# if defined MRL_WCI && defined WAVE_OFFLINE_FRICTION
            wwEpd=wwed_scale*wwedp(itww)
# endif
# if defined MRL_WCI && defined WAVE_OFFLINE_ROLLER   
            wwEpr=wwer_scale*wwerp(itww)
# endif
          
            do j=JstrR,JendR
              do i=IstrR,IendR
                Awave(i,j)=(wwagrd)*wwa_scale*wwag(i,j,itww) +
     &                     (1.-wwagrd)*wwA
                Dwave(i,j)=(wwdgrd)*wwd_scale*wwdg(i,j,itww) +
     &                     (1.-wwdgrd)*wwD 

                Pwave(i,j)=(wwpgrd)*wwp_scale*wwpg(i,j,itww) +
     &                     (1.-wwpgrd)*wwP
# ifdef MUSTANG
                Uwave(i,j)=(wwugrd)*wwu_scale*wwug(i,j,itww) +
     &                     (1.-wwugrd)*wwU
# endif
# if defined MRL_WCI && defined WAVE_OFFLINE_BREAKING
                Ebwave(i,j)=(wwebgrd)*wweb_scale*wweb(i,j,itww) +
     &                     (1.-wwebgrd)*wwEpb
# endif
# if defined MRL_WCI && defined WAVE_OFFLINE_FRICTION
                Edwave(i,j)=(wwedgrd)*wwed_scale*wwed(i,j,itww) +
     &                     (1.-wwedgrd)*wwEpd
# endif
# if defined MRL_WCI && defined WAVE_OFFLINE_ROLLER
                Erwave(i,j)=(wwergrd)*wwer_scale*wwer(i,j,itww) +
     &                     (1.-wwergrd)*wwEpr
# endif
              enddo
            enddo
        endif
!
! Time-interpolate wave data from gridded or point data.
! Make sure that for the next time step [when time=time+dt]
! time+dt is still between wwv_time(it1) and wwv_time(it2);
! and if not, set synchro_flag to signal that the new forcing data
! should be read from the netCDF input file (master thread only).
!
      elseif (fac1.ge.0. .and. fac2.ge.0.) then
       if (ZEROTH_TILE .and. fac1.lt.dt) synchro_flag=.TRUE.
      
        fac=wwa_scale/(fac1+fac2)
        fac11=fac*fac1
        fac12=fac*fac2
        wwA=fac11*wwap(it1)+fac12*wwap(it2)

        fac=wwd_scale/(fac1+fac2)
        fac21=fac*fac1
        fac22=fac*fac2
        wwD=fac21*wwdp(it1)+fac22*wwdp(it2)

        fac=wwp_scale/(fac1+fac2)
        fac31=fac*fac1
        fac32=fac*fac2
        wwP=fac31*wwpp(it1)+fac32*wwpp(it2)
# ifdef MUSTANG
        fac=wwu_scale/(fac1+fac2)
        fac71=fac*fac1
        fac72=fac*fac2
        wwU=fac71*wwup(it1)+fac72*wwup(it2)
# endif

# if defined MRL_WCI && defined WAVE_OFFLINE_BREAKING
        fac=wweb_scale/(fac1+fac2)
        fac41=fac*fac1
        fac42=fac*fac2
        wwEpb=fac41*wwebp(it1)+fac42*wwebp(it2)
# endif
# if defined MRL_WCI && defined WAVE_OFFLINE_FRICTION
        fac=wwed_scale/(fac1+fac2)
        fac51=fac*fac1
        fac52=fac*fac2
        wwEpd=fac51*wwedp(it1)+fac52*wwedp(it2)
# endif
# if defined MRL_WCI && defined WAVE_OFFLINE_ROLLER
        fac=wwer_scale/(fac1+fac2)
        fac61=fac*fac1
        fac62=fac*fac2
        wwEpr=fac61*wwerp(it1)+fac62*wwerp(it2)
# endif
        
        do j=JstrR,JendR
          do i=IstrR,IendR
            Awave(i,j)=wwagrd*(fac11*wwag(i,j,it1)+fac12*wwag(i,j,it2))
     &                +(1.-wwagrd)*wwA
            Dwave(i,j)=wwdgrd*(fac21*wwdg(i,j,it1)+fac22*wwdg(i,j,it2))
     &                +(1.-wwdgrd)*wwD 

            Pwave(i,j)=wwpgrd*(fac31*wwpg(i,j,it1)+fac32*wwpg(i,j,it2))
     &                +(1.-wwpgrd)*wwP
# ifdef MUSTANG
            Uwave(i,j)=wwugrd*(fac71*wwug(i,j,it1)+fac72*wwug(i,j,it2))
     &                +(1.-wwugrd)*wwU
# endif
# if defined MRL_WCI && defined WAVE_OFFLINE_BREAKING
            wepb(i,j)=wwebgrd*(fac41*wweb(i,j,it1)+fac42*wweb(i,j,it2))
     &                +(1.-wwebgrd)*wwEpb
# endif        
# if defined MRL_WCI && defined WAVE_OFFLINE_FRICTION
            wepd(i,j)=wwedgrd*(fac51*wwed(i,j,it1)+fac52*wwed(i,j,it2))
     &                +(1.-wwedgrd)*wwEpd
# endif
# if defined MRL_WCI && defined WAVE_OFFLINE_ROLLER
            wepr(i,j)=wwergrd*(fac61*wwer(i,j,it1)+fac62*wwer(i,j,it2))
     &                +(1.-wwergrd)*wwEpr
# endif
          enddo
        enddo
      endif
!
!--------------------------------------------------------------------
!  Convert wave direction Dwave to Cartesian convention
!--------------------------------------------------------------------
!
!  Wave direction conventions:
!  - Cartesian convention (CROCO): 
!          the direction to where the vector points, 
!          measured counterclockwise from the positive x-axis
!          of the model grid (must account for grid angle)
!  - Meteorological convention (i.e., ECMWF data ...): 
!          the direction where the waves come from, 
!          measured clockwise from geographic North.
!  - Oceanographic convention: 
!          the direction where the waves travel to, 
!          measured clockwise from geographic North.
!
# if defined CURVGRID || !defined SHOREFACE
      do j=jstrR,jendR
        do i=istrR,iendR
#  if !defined SHOREFACE && !defined MUSTANG 
!         Dwave(i,j)=0.5*pi+Dwave(i,j) ! Dwave=0 when waves to N (oceano)
          Dwave(i,j)=1.5*pi-Dwave(i,j) ! Dwave=0 when waves from N (meteo)
#  endif
#  ifdef CURVGRID
#   ifdef MUSTANG
          Dwave(i,j)=Dwave(i,j)-angler(i,j)*rad2deg
#   else
          Dwave(i,j)=Dwave(i,j)-angler(i,j)
#   endif
#  endif
        enddo
      enddo
# endif /* CURVGRID || !SHOREFACE */
!
!--------------------------------------------------------------------
! Simulate nearshore modifications for low resolution data
!
!   Mekong application provides an example that can be used
!   for other muddy areas (still needs work).
!--------------------------------------------------------------------
!
# if !defined MRL_WCI || defined MEKONG
      do j=jstrR,jendR
        do i=istrR,iendR
          Dstp=max(h(i,j)+z_w(i,j,N),1.e-3)
#  ifdef MEKONG
!
! First estimate wavenumber
!
          hrms=2*Awave(i,j)
          cfrq=2.0*pi/max(Pwave(i,j),0.1D0)
          khd=Dstp*cfrq*cfrq/g
          kh = sqrt(    khd*khd + khd/(1.0 + khd*(0.6666666666
     &                  +khd*(0.3555555555 + khd*(0.1608465608
     &                  +khd*(0.0632098765 + khd*(0.0217540484
     &                                  +khd*0.0065407983)))))) )
          kw=kh/Dstp
!
! Strong mudbank dissipation formulation: frequency and depth-dependent 
! dissipation rate is a function akin to that proposed by Bretschneider 
! and Reid (1954) and Grosskopf (1980): longer and higher waves are more 
! highly dissipated.
!
          Awave(i,j)=Awave(i,j)*
     &                (.5*(1-tanh(1.-2.*kw**.9*Dstp/hrms**0.2)))**2 
#  else
!
! Simulate surf-zone dissipation due to depth-induced wave breaking
!
          Awave(i,j)=min(0.2*Dstp, Awave(i,j))
#  endif
        enddo
      enddo
# endif /* !MRL_WCI || MEKONG */
!
!--------------------------------------------------------------------
! Compute variables for wave-current interactions (mrl_wci routine)
! --> store in same arrays as in OW_COUPLING case
!--------------------------------------------------------------------
!
# ifdef MRL_WCI
      do j=jstrR,jendR
        do i=istrR,iendR
          whrm(i,j)=2*Awave(i,j)                  ! Amp to Hrms conversion
          wfrq(i,j)=2.0*pi/max(Pwave(i,j),0.1D0)  ! Pwave peak period
          wdrx(i,j)=cos(Dwave(i,j))             
          wdre(i,j)=sin(Dwave(i,j))
#  ifdef MASKING
          wfrq(i,j)=wfrq(i,j)*rmask(i,j)
          whrm(i,j)=whrm(i,j)*rmask(i,j)
          wdrx(i,j)=wdrx(i,j)*rmask(i,j)
          wdre(i,j)=wdre(i,j)*rmask(i,j)
#   ifdef WAVE_OFFLINE_BREAKING
          wepb(i,j)=wepb(i,j)*rmask(i,j)
#   endif
#   ifdef WAVE_OFFLINE_FRICTION
          wepd(i,j)=wepd(i,j)*rmask(i,j)
#   endif
#   ifdef WAVE_OFFLINE_ROLLER
          wepr(i,j)=wepr(i,j)*rmask(i,j)
#   endif
#  endif
        enddo
      enddo
# endif /* MRL_WCI */
!
!--------------------------------------------------------------------
! Unable to set-up wave fields:
! Complain about the error and signal to quit (ONE THREAD ONLY).
!--------------------------------------------------------------------
!
      if ((ww_cycle.gt.0.) .and. (fac1.lt.0..or.fac2.lt.0.)) then
        if (ZEROTH_TILE) then
          write(stdout,1) 'wwv_time', tdays, wwv_time(it2)*sec2day
# ifdef USE_CALENDAR
     &  -ww_origin_date_in_sec*sec2day
# endif
  1       format(/,' SET_WWAVE_TILE - current model time exceeds',
     &             ' ending value for variable: ',a,/,14x,'TDAYS = ',
     &                                      g12.4,2x,'TEND = ',g12.4)
          may_day_flag=2
        endif
      endif
      return
      end
#else
      subroutine get_wwave_empty
      return
      end
#endif /* BBL && !ANA_WWAVE */
