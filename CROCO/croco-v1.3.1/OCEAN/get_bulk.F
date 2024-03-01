! $Id: get_bulk.F 1564 2014-06-24 17:39:21Z gcambon $
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
#if defined BULK_FLUX

      subroutine get_bulk
!
!  Read in wind speed and surface air temperature
!
!

      implicit none
# include "param.h"
# include "forces.h"
# include "scalars.h"
# include "netcdf.inc"
# include "ncscrum.h"
      real    cff
      integer i,ierr, lstr,lvar,lenstr, nf_fread, advance_cycle
!
! Initialization: Inquire about the contents of forcing NetCDF file:
!================  variables and dimensions. Check for consistency.
!
      if (may_day_flag.ne.0) return      !-->  EXIT
      if (itbulk.eq.0 .or. iic.eq.0) then
        lstr=lenstr(bulkname)
c*        call opencdf (bulkname,N)
c*        if (may_day_flag.ne.0) return    !-->  EXIT
!
! If not opened yet, open forcing NetCDF file for reading.
! Find and save IDs for relevant variables, determine whether
! SST is a field or scalar value.
!
        if (ncidbulk.eq.-1) then
          ierr=nf_open(bulkname(1:lstr), nf_nowrite, ncidbulk)
          if (ierr. ne. nf_noerr) goto 4                  !--> ERROR
        endif

        ierr=nf_inq_varid (ncidbulk, 'bulk_time', bulk_tid)
        if (ierr .ne. nf_noerr) then
          write(stdout,3) 'bulk_time', bulkname(1:lstr)
          goto 99                                         !--> ERROR
        endif

        lvar=lenstr(vname(1,indxTAIR))
        ierr=nf_inq_varid (ncidbulk, vname(1,indxTAIR)(1:lvar),tair_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidbulk, tair_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              ltairgrd=1
            else
              ltairgrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxTAIR)(1:lvar),bulkname(1:lstr)
          goto 99                                         !--> ERROR
        endif

        lvar=lenstr(vname(1,indxRHUM))
        ierr=nf_inq_varid (ncidbulk, vname(1,indxRHUM)(1:lvar),rhum_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidbulk, rhum_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lrhumgrd=1
            else
              lrhumgrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxRHUM)(1:lvar),bulkname(1:lstr)
          goto 99                                         !--> ERROR
        endif

        lvar=lenstr(vname(1,indxRADLW))
        ierr=nf_inq_varid (ncidbulk,vname(1,indxRADLW)(1:lvar),radlw_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidbulk, radlw_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lradlwgrd=1
            else
              lradlwgrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxRADLW)(1:lvar),bulkname(1:lstr)
          goto 99                                         !--> ERROR
        endif

        lvar=lenstr(vname(1,indxShflx_rsw))
        ierr=nf_inq_varid (ncidbulk,vname(1,indxShflx_rsw)
     &                                                (1:lvar),radsw_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidbulk, radsw_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lradswgrd=1
            else
              lradswgrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxShflx_rsw)(1:lvar),
     &                                                bulkname(1:lstr)
          goto 99                                        !--> ERROR
        endif

# ifdef READ_PATM
        lvar=lenstr(vname(1,indxPATM))
        ierr=nf_inq_varid (ncidbulk, vname(1,indxPATM)(1:lvar),patm_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidbulk, patm_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lpatmgrd=1
            else
              lpatmgrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxPATM)(1:lvar),bulkname(1:lstr)
          goto 99                                         !--> ERROR
        endif
# endif
# ifdef DIURNAL_INPUT_SRFLX
        lvar=lenstr(vname(1,indxShflx_rswbio))
        ierr=nf_inq_varid (ncidbulk,vname(1,indxShflx_rswbio)(1:lvar)
     &                                              ,radswbio_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidbulk, radswbio_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lradswbiogrd=1
            else
              lradswbiogrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxShflx_rswbio)(1:lvar),
     &                                              bulkname(1:lstr)
          goto 99                                         !--> ERROR
        endif
# endif /* DIURNAL_INPUT_SRFLX */

# ifdef SALINITY
        lvar=lenstr(vname(1,indxPRATE))
        ierr=nf_inq_varid (ncidbulk,vname(1,indxPRATE)(1:lvar),prate_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidbulk, prate_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lprategrd=1
            else
              lprategrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxPRATE)(1:lvar),bulkname(1:lstr)
          goto 99                                         !--> ERROR
        endif
# endif /* SALINITY */

        lvar=lenstr(vname(1,indxUWND))
        ierr=nf_inq_varid (ncidbulk,vname(1,indxUWND)(1:lvar),uwnd_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidbulk, uwnd_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              luwndgrd=1
            else
              luwndgrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxUWND)(1:lvar),bulkname(1:lstr)
          goto 99                                         !--> ERROR
        endif

        lvar=lenstr(vname(1,indxVWND))
        ierr=nf_inq_varid (ncidbulk,vname(1,indxVWND)(1:lvar),vwnd_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidbulk, vwnd_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lvwndgrd=1
            else
              lvwndgrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxVWND)(1:lvar),bulkname(1:lstr)
          goto 99                                         !--> ERROR
        endif
!
! Determine whether there is cycling to reuse the input data and
! find cycling period "bulk_cycle", set initial cycling index
! "wspd_ncycle" and record index "wspd_rec".
! Set initial value for time index "itbulk" and both time record
! bounds to large negative artificial values, so that it will
! trigger the logic in reading part below.
! Also set scale factor to convert input dQdSST from Watts/m2/Celsius
! to meter/second.
!
        call set_cycle (ncidbulk,   bulk_tid,    ntbulk,
     &                  bulk_cycle, bulk_ncycle, bulk_rec)
        if (may_day_flag.ne.0) return    !-->  EXIT
        itbulk=2
        bulk_time(1)=-1.E+20
        bulk_time(2)=-1.E+20
        srf_scale=1./(rho0*Cp)
# ifdef SALINITY
        stf_scale(isalt)=0.01/86400.
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
  1   i=3-itbulk
       cff=time+0.5*dt
        if (bulk_time(i).le.cff .and. cff.lt.bulk_time(itbulk))
     &            return
        ierr=advance_cycle (bulk_cycle,ntbulk,bulk_ncycle,bulk_rec)
        if (ierr .ne. 0) goto 7                           !--> ERROR

        ierr=nf_get_var1_FTYPE (ncidbulk, bulk_tid, bulk_rec, cff)
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'bulk_time', bulk_rec
          goto 99                                         !--> ERROR
        endif
# ifdef USE_CALENDAR
        call tool_origindate(ncidbulk,bulk_tid,
     &                       blk_origin_date_in_sec)
        cff=cff+blk_origin_date_in_sec*sec2day
# endif

        bulk_time(i)=cff*day2sec+bulk_cycle*bulk_ncycle
        if (bulk_time(itbulk).eq.-1.E+20)
     &          bulk_time(itbulk)=bulk_time(i)
!
! tair
!
        if (ltairgrd.eq.1) then
          ierr=nf_fread (tairg(START_2D_ARRAY,i), ncidbulk, tair_id,
     &                                            bulk_rec,  r2dvar)
        else
          ierr=nf_get_var1_FTYPE (ncidbulk,tair_id,bulk_rec,tairp(i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'TAIR', bulk_rec
          goto 99                                         !--> ERROR
        endif
!
! rhum
!
        if (lrhumgrd.eq.1) then
          ierr=nf_fread (rhumg(START_2D_ARRAY,i), ncidbulk, rhum_id,
     &                                            bulk_rec,  r2dvar)
        else
          ierr=nf_get_var1_FTYPE (ncidbulk,rhum_id,bulk_rec,rhump(i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'RHUM', bulk_rec
          goto 99                                         !--> ERROR
        endif
!
! radlw
!
        if (lradlwgrd.eq.1) then
          ierr=nf_fread (radlwg(START_2D_ARRAY,i), ncidbulk, radlw_id,
     &                                            bulk_rec,  r2dvar)
        else
          ierr=nf_get_var1_FTYPE (ncidbulk,radlw_id,bulk_rec,radlwp(i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'RADLW', bulk_rec
          goto 99                                         !--> ERROR
        endif
!
! radsw
!
        if (lradswgrd.eq.1) then
          ierr=nf_fread (radswg(START_2D_ARRAY,i), ncidbulk, radsw_id,
     &                                            bulk_rec,  r2dvar)
        else
          ierr=nf_get_var1_FTYPE (ncidbulk,radsw_id,bulk_rec,radswp(i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'RADSW', bulk_rec
          goto 99                                         !--> ERROR
        endif
!
! radswbio
!
# ifdef READ_PATM
! Patm
!
        if (lpatmgrd.eq.1) then
          ierr=nf_fread (patmg(START_2D_ARRAY,i), ncidbulk, patm_id,
     &                                            bulk_rec,  r2dvar)
        else
          ierr=nf_get_var1_FTYPE (ncidbulk,patm_id,bulk_rec,patmp(i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'PATM', bulk_rec
          goto 99                                         !--> ERROR
        endif

# endif
# ifdef DIURNAL_INPUT_SRFLX
        if (lradswbiogrd.eq.1) then
          ierr=nf_fread (radswbiog(START_2D_ARRAY,i), ncidbulk,
     &                        radswbio_id, bulk_rec,  r2dvar)
        else
          ierr=nf_get_var1_FTYPE (ncidbulk,radswbio_id,
     &                        bulk_rec, radswbiop(i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'RADSWBIO', bulk_rec
          goto 99                                         !--> ERROR
        endif
# endif /* DIURNAL_INPUT_SRFLX */
!
! prate
!
# ifdef SALINITY
        if (lprategrd.eq.1) then
          ierr=nf_fread (prateg(START_2D_ARRAY,i), ncidbulk, prate_id,
     &                                            bulk_rec,  r2dvar)
        else
          ierr=nf_get_var1_FTYPE (ncidbulk,prate_id,bulk_rec,pratep(i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'PRATE', bulk_rec
          goto 99                                         !--> ERROR
        endif
# endif
!
! uwnd
!
        if (luwndgrd.eq.1) then
          ierr=nf_fread(uwndg(START_2D_ARRAY,i), ncidbulk, uwnd_id,
     &                                              bulk_rec, u2dvar)
        else
          ierr=nf_get_var1_FTYPE(ncidbulk,uwnd_id,bulk_rec,uwndp(i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'uwnd', bulk_rec
          goto 99                                         !--> ERROR
        endif
!
! vwnd
!
        if (lvwndgrd.eq.1) then
          ierr=nf_fread(vwndg(START_2D_ARRAY,i), ncidbulk, vwnd_id,
     &                                              bulk_rec, v2dvar)
        else
          ierr=nf_get_var1_FTYPE(ncidbulk,vwnd_id,bulk_rec,vwndp(i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'vwnd', bulk_rec
          goto 99                                         !--> ERROR
        endif
       itbulk=i
       MPI_master_only write(stdout,'(6x,A,1x,A,1x,F10.4,1x,I4)')
     &              'GET_BULK   --',
     &              'Read fields for bulk formula   for time =', cff
# ifdef USE_CALENDAR
     &  -blk_origin_date_in_sec*sec2day
# endif
# ifdef MPI
     &                                                      , mynode
# endif
      if (ntbulk.gt.1) goto 1
      if (ntbulk.eq.1) return
!
! Sort out error messages: The following portion of the code is
!===== === ===== ========= not accessed unless something goes wrong.
!
  3   format(/,' GET_BULK - ERROR: unable to find forcing variable',
     &                   ': ',a,/,11x,'in forcing NetCDF file: ',a)
  4   write(stdout,5) bulkname(1:lstr)
  5   format(/,' GET_BULK - ERROR: unable to open forcing NetCDF ',
     &                                                  'file: ',a)
      goto 99
  6   format(/,' GET_BULK - ERROR while reading variable: ',a,2x,
     &                                      ' at TIME index = ',i4)
  7   write(stdout,8) bulk_rec, ntbulk, bulkname(1:lstr), tdays,
     &                                     bulk_time(itbulk)*sec2day
# ifdef USE_CALENDAR
     &  -blk_origin_date_in_sec*sec2day
# endif

  8   format(/,' GET_BULK - ERROR: requested time record ',I4,
     &       1x,'exceeds the last available',  /,  11x,'record ',I4,
     &       1x,'in forcing NetCDF file: ', a, /,   11x,'TDAYS = ',
     &               g12.4,2x,'last available BULK_TIME = ',g12.4)
  99  may_day_flag=2
      return
      end
!
!======================================================================
!
      subroutine set_bulk_tile (Istr,Iend,Jstr,Jend)
!
!======================================================================
!
! Set-up bulk data for current tile.
!
      implicit none
# include "param.h"
# include "forces.h"
# include "scalars.h"
# include "grid.h"
# ifdef CFB_WIND_TRA
#  include "ocean3d.h"
#  include "params_bulk.h"
# endif
      integer Istr,Iend,Jstr,Jend, i,j, it1,it2
      integer imax,jmax
      real cff,cff1,cff2, cff3,cff4
# ifdef SALINITY
      real cff5,cff6
# endif
      real val1,val2,val3,val4,val5,val6,val7,val8
# ifdef CFB_WIND_TRA
      real uwnd_r, vwnd_r
# endif
# ifdef DIURNAL_INPUT_SRFLX
      real val44
# endif /* DIURNAL_INPUT_SRFLX */
!
# include "compute_extended_bounds.h"
!
!
! Extended range is needed for KPP
!
# ifdef EW_PERIODIC
      imax=Iend+2
# else
      if (EASTERN_EDGE) then
        imax=Iend+1
      else
        imax=Iend+2
      endif
# endif
# ifdef NS_PERIODIC
      jmax=Jend+2
# else
      if (NORTHERN_EDGE) then
        jmax=Jend+1
      else
        jmax=Jend+2
      endif
# endif
!
      it1=3-itbulk
      it2=itbulk
      cff=time+0.5*dt
      cff1=bulk_time(it2)-cff
      cff2=cff-bulk_time(it1)
!
!----------------------------------------------------------------------
!  Load time invariant
!----------------------------------------------------------------------
!
      if (bulk_cycle.lt.0.) then

        if (FIRST_RST_TIME_STEP) then
          if (ltairgrd.eq.1) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                tair(i,j)=tairg(i,j,itbulk)
                rhum(i,j)=rhumg(i,j,itbulk)
                radlw(i,j)=srf_scale*radlwg(i,j,itbulk)
                radsw(i,j)=srf_scale*radswg(i,j,itbulk)
                srflx(i,j)=radsw(i,j)
# ifdef READ_PATM
                patm2d(i,j)=patmg(i,j,itbulk)
# endif
# ifdef DIURNAL_INPUT_SRFLX
                radswbio(i,j)=srf_scale*radswbiog(i,j,itbulk)
                srflxbio(i,j)=radswbio(i,j)
# endif /* DIURNAL_INPUT_SRFLX */
# ifdef SALINITY
                prate(i,j)=stf_scale(isalt)*prateg(i,j,itbulk)
# endif
              enddo
            enddo
            
            do j=JstrR,JendR
              do i=IstrR,imax
                uwnd(i,j)=uwndg(i,j,itbulk)
              enddo
            enddo
            do j=JstrR,jmax
              do i=IstrR,IendR
                vwnd(i,j)=vwndg(i,j,itbulk)
              enddo
            enddo
            do j=JstrR,min(JendR,Mm)
              do i=IstrR,min(IendR,Lm)
                wspd(i,j)=sqrt(0.25*(uwnd(i,j)+uwnd(i+1,j))**2+
     &                         0.25*(vwnd(i,j)+vwnd(i,j+1))**2)
# ifdef CFB_WIND_TRA
                cff      = 1.-swparam    ! current-wind coupling parameter: Ua => Ua-(1-sw)Uo
                uwnd_r   = 0.5*( uwnd(i+1,j)+uwnd(i,j)
     &                   - cff*( u(i+1,j,N,nrhs)+u(i,j,N,nrhs))
     &                   )
                vwnd_r   = 0.5*( vwnd(i,j+1)+vwnd(i,j)
     &                   - cff*( v(i,j+1,N,nrhs)+v(i,j,N,nrhs))
     &                   )
                wspd_cfb(i,j) = SQRT( uwnd_r*uwnd_r+vwnd_r*vwnd_r )
# endif     
              enddo
            enddo
# ifndef EW_PERIODIC
            if (WESTERN_EDGE) then
              do j=Jstr,Jend
                wspd(Istr-1,j)=wspd(Istr,j)
#  ifdef CFB_WIND_TRA
                wspd_cfb(Istr-1,j)=wspd_cfb(Istr,j)
#  endif                 
              enddo
            endif
            if (EASTERN_EDGE) then
              do j=Jstr,Jend
                wspd(Iend+1,j)=wspd(Iend,j)
#  ifdef CFB_WIND_TRA
                wspd_cfb(Iend+1,j)=wspd_cfb(Iend,j)       
#  endif            
              enddo
            endif
# endif
# ifndef NS_PERIODIC
            if (SOUTHERN_EDGE) then
              do i=Istr,Iend
                wspd(i,Jstr-1)=wspd(i,Jstr)
#  ifdef CFB_WIND_TRA 
                wspd_cfb(i,Jstr-1)=wspd_cfb(i,Jstr)
#  endif                  
              enddo
            endif
            if (NORTHERN_EDGE) then
              do i=Istr,Iend
                wspd(i,Jend+1)=wspd(i,Jend)
#  ifdef CFB_WIND_TRA 
                wspd_cfb(i,Jend+1)=wspd_cfb(i,Jend)
#  endif                  
              enddo
            endif
#  ifndef EW_PERIODIC
            if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
              wspd(Istr-1,Jstr-1)=wspd(Istr,Jstr)
#   ifdef CFB_WIND_TRA 
              wspd_cfb(Istr-1,Jstr-1)=wspd_cfb(Istr,Jstr)
#   endif                
            endif
            if (WESTERN_EDGE .and. NORTHERN_EDGE) then
              wspd(Istr-1,Jend+1)=wspd(Istr,Jend)
#   ifdef CFB_WIND_TRA  
              wspd_cfb(Istr-1,Jend+1)=wspd_cfb(Istr,Jend)
#   endif               
            endif
            if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
              wspd(Iend+1,Jstr-1)=wspd(Iend,Jstr)
#   ifdef CFB_WIND_TRA 
              wspd_cfb(Iend+1,Jstr-1)=wspd_cfb(Iend,Jstr)
#   endif                
            endif
            if (EASTERN_EDGE .and. NORTHERN_EDGE) then
              wspd(Iend+1,Jend+1)=wspd(Iend,Jend)
#   ifdef CFB_WIND_TRA 
              wspd_cfb(Iend+1,Jend+1)=wspd_cfb(Iend,Jend)
#   endif                
            endif
#  endif
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
            call exchange_r2d_tile(Istr,Iend,Jstr,Jend,
     &                             wspd(START_2D_ARRAY))
#  ifdef CFB_WIND_TRA 
            call exchange_r2d_tile(Istr,Iend,Jstr,Jend,
     &                             wspd_cfb(START_2D_ARRAY))
#  endif      
# endif
          else
            
            val1=tairp(itbulk)
            val2=rhump(itbulk)
            val3=srf_scale*radlwp(itbulk)
            val4=srf_scale*radswp(itbulk)
# ifdef DIURNAL_INPUT_SRFLX
            val44=srf_scale*radswbiop(itbulk)
# endif /* DIURNAL_INPUT_SRFLX */
# ifdef SALINITY
            val5=stf_scale(isalt)*pratep(itbulk)
# endif
            val6=uwndp(itbulk)
            val7=vwndp(itbulk)
# ifdef READ_PATM
            val8=patmp(itbulk)
# endif
            do j=JstrR,JendR
              do i=IstrR,IendR
                tair(i,j)=val1
                rhum(i,j)=val2
                radlw(i,j)=val3
                radsw(i,j)=val4
                srflx(i,j)=val4
# ifdef READ_PATM
                patm2d(i,j)=val8
# endif
# ifdef DIURNAL_INPUT_SRFLX
                radswbio(i,j)=val44
                srflxbio(i,j)=val44
# endif /* DIURNAL_INPUT_SRFLX */
# ifdef SALINITY
                prate(i,j)=val5
# endif
                uwnd(i,j)=val6
                vwnd(i,j)=val7
                wspd(i,j) = SQRT(val6**2.0+val7**2.0)
              enddo
            enddo
          endif
        endif
!     
!----------------------------------------------------------------------
!     Time interpolation
!
!     Time-interpolate SST and dQdSST from grided or point data.
!     Check that for the next time step [when time=time+dt] time+dt
!     is still between wspd_time(it1) and wspd_time(it2); and if not,
!     set synchro_flag top signal that the new forcing data should be
!     read from the netCDF input file (master thread only).
!----------------------------------------------------------------------
!
      elseif (cff1.ge.0. .and. cff2.ge.0.) then

        if (ZEROTH_TILE .and. cff1.lt.dt) synchro_flag=.TRUE.
        !note cff order maters
        cff=srf_scale/(cff1+cff2)
        cff3=cff1*cff
        cff4=cff2*cff
# ifdef SALINITY
        cff=stf_scale(isalt)/(cff1+cff2)
        cff5=cff1*cff
        cff6=cff2*cff
# endif
        cff=1./(cff1+cff2)
        cff1=cff1*cff
        cff2=cff2*cff

        if (ltairgrd.eq.1) then

          do j=JstrR,JendR
            do i=IstrR,IendR
              tair(i,j)=cff1*tairg(i,j,it1)+cff2*tairg(i,j,it2)
              rhum(i,j)=cff1*rhumg(i,j,it1)+cff2*rhumg(i,j,it2)
              radlw(i,j)=cff3*radlwg(i,j,it1)+cff4*radlwg(i,j,it2)
              radsw(i,j)=cff3*radswg(i,j,it1)+cff4*radswg(i,j,it2)
              srflx(i,j)=radsw(i,j)
# ifdef READ_PATM
              patm2d(i,j)=cff1*patmg(i,j,it1)+cff2*patmg(i,j,it2)
# endif
# ifdef DIURNAL_INPUT_SRFLX
              radswbio(i,j)=cff3*radswbiog(i,j,it1)
     &             +cff4*radswbiog(i,j,it2)
              srflxbio(i,j)=radswbio(i,j)
# endif /* DIURNAL_INPUT_SRFLX */
# ifdef SALINITY
              prate(i,j)=cff5*prateg(i,j,it1)+cff6*prateg(i,j,it2)
# endif
            enddo
          enddo

          do j=JstrR,JendR
            do i=IstrR,imax
              uwnd(i,j)=cff1*uwndg(i,j,it1)+cff2*uwndg(i,j,it2)
            enddo
          enddo
          do j=JstrR,jmax
            do i=IstrR,IendR
              vwnd(i,j)=cff1*vwndg(i,j,it1)+cff2*vwndg(i,j,it2)
            enddo
          enddo
          do j=JstrR,min(JendR,Mm)
            do i=IstrR,min(IendR,Lm)
              wspd(i,j)=sqrt(0.25*(uwnd(i,j)+uwnd(i+1,j))**2+
     &                       0.25*(vwnd(i,j)+vwnd(i,j+1))**2)
# ifdef CFB_WIND_TRA
              cff      = 1.-swparam    ! current-wind coupling parameter: Ua => Ua-(1-sw)Uo
              uwnd_r   = 0.5*( uwnd(i+1,j)+uwnd(i,j)
     &                 - cff*( u(i+1,j,N,nrhs)+u(i,j,N,nrhs))
     &                   )
              vwnd_r   = 0.5*( vwnd(i,j+1)+vwnd(i,j)
     &                 - cff*( v(i,j+1,N,nrhs)+v(i,j,N,nrhs))
     &                   )
              wspd_cfb(i,j) = SQRT( uwnd_r*uwnd_r+vwnd_r*vwnd_r )
# endif
            enddo
          enddo
# ifndef EW_PERIODIC
          if (WESTERN_EDGE) then
            do j=Jstr,Jend
              wspd(Istr-1,j)=wspd(Istr,j)
#  ifdef CFB_WIND_TRA
              wspd_cfb(Istr-1,j)=wspd_cfb(Istr,j)
#  endif
            enddo
          endif
          if (EASTERN_EDGE) then
            do j=Jstr,Jend
              wspd(Iend+1,j)=wspd(Iend,j)
#  ifdef CFB_WIND_TRA
              wspd_cfb(Iend+1,j)=wspd_cfb(Iend,j)
#  endif
            enddo
          endif
# endif
# ifndef NS_PERIODIC
          if (SOUTHERN_EDGE) then
            do i=Istr,Iend
              wspd(i,Jstr-1)=wspd(i,Jstr)
#  ifdef CFB_WIND_TRA
              wspd_cfb(i,Jstr-1)=wspd_cfb(i,Jstr)
#  endif
            enddo
          endif
          if (NORTHERN_EDGE) then
            do i=Istr,Iend
              wspd(i,Jend+1)=wspd(i,Jend)
#  ifdef CFB_WIND_TRA
              wspd_cfb(i,Jend+1)=wspd_cfb(i,Jend)
#  endif 
            enddo
          endif
#  ifndef EW_PERIODIC
          if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
            wspd(Istr-1,Jstr-1)=wspd(Istr,Jstr)
#   ifdef CFB_WIND_TRA
            wspd_cfb(Istr-1,Jstr-1)=wspd_cfb(Istr,Jstr)
#   endif
          endif
          if (WESTERN_EDGE .and. NORTHERN_EDGE) then
            wspd(Istr-1,Jend+1)=wspd(Istr,Jend)
#   ifdef CFB_WIND_TRA
            wspd_cfb(Istr-1,Jend+1)=wspd_cfb(Istr,Jend)   
#   endif          
          endif
          if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
            wspd(Iend+1,Jstr-1)=wspd(Iend,Jstr)
#   ifdef CFB_WIND_TRA 
            wspd_cfb(Iend+1,Jstr-1)=wspd_cfb(Iend,Jstr)
#   endif            
          endif
          if (EASTERN_EDGE .and. NORTHERN_EDGE) then
            wspd(Iend+1,Jend+1)=wspd(Iend,Jend)
#   ifdef CFB_WIND_TRA 
            wspd_cfb(Iend+1,Jend+1)=wspd_cfb(Iend,Jend)
#   endif           
          endif
#  endif
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
          call exchange_r2d_tile(Istr,Iend,Jstr,Jend,
     &                           wspd(START_2D_ARRAY))
#  ifdef CFB_WIND_TRA 
          call exchange_r2d_tile(Istr,Iend,Jstr,Jend,
     &                           wspd_cfb(START_2D_ARRAY)) 
#  endif    
# endif

        else
          
          val1=cff1*tairp(it1)+cff2*tairp(it2)
          val2=cff1*rhump(it1)+cff2*rhump(it2)
          val3=cff3*radlwp(it1)+cff4*radlwp(it2)
          val4=cff3*radswp(it1)+cff4*radswp(it2)
# ifdef DIURNAL_INPUT_SRFLX
          val44=cff3*radswbiop(it1)+cff4*radswbiop(it2)
# endif /* DIURNAL_INPUT_SRFLX */
# ifdef SALINITY
          val5=cff5*pratep(it1)+cff6*pratep(it2)
# endif
          val6=cff1*uwndp(it1)+cff2*uwndp(it2)
          val7=cff1*vwndp(it1)+cff2*vwndp(it2)
# ifdef READ_PATM
          val8=cff1*patmp(it1)+cff2*patmp(it2)
# endif

          do j=JstrR,JendR
            do i=IstrR,IendR
              tair(i,j)=val1
              rhum(i,j)=val2
              radlw(i,j)=val3
              radsw(i,j)=val4
              srflx(i,j)=val4
# ifdef DIURNAL_INPUT_SRFLX
              radswbio(i,j)=val44
              srflxbio(i,j)=val44
# endif /* DIURNAL_INPUT_SRFLX */
# ifdef SALINITY
              prate(i,j)=val5
# endif
              uwnd(i,j)=val6
              vwnd(i,j)=val7
# ifdef READ_PATM
              patm2d(i,j)=val8
# endif
              wspd(i,j)=sqrt(uwnd(i,j)**2.0+vwnd(i,j)**2.0)
            enddo
          enddo
        endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
          call exchange_r2d_tile(Istr,Iend,Jstr,Jend,
     &                           wspd(START_2D_ARRAY))
#  ifdef CFB_WIND_TRA 
          call exchange_r2d_tile(Istr,Iend,Jstr,Jend,
     &                           wspd_cfb(START_2D_ARRAY)) 
#  endif              
# endif
!     
!     Unable to set-up SST and dQdSST:
!     Complain about the error and signal to quit.
!
      else
        if (ZEROTH_TILE) then
          write(stdout,1) 'bulk_time',tdays,bulk_time(it2)*sec2day
# ifdef USE_CALENDAR
     &  -blk_origin_date_in_sec*sec2day
# endif
 1        format(/,' SET_BULK - current model time exceeds ending',
     &           1x,'value for variable: ',a,/,11x,'TDAYS = ',g12.4,
     &         2x,'TEND = ',g12.4)
          may_day_flag=2
        endif
      endif
      return
      end
#else
      subroutine get_bulk_empty
      return
      end
#endif /* BULK_FLUX */
