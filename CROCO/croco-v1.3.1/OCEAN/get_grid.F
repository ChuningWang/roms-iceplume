! $Id: get_grid.F 1458 2014-02-03 15:01:25Z gcambon $
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
#ifndef ANA_GRID
                             ! Read grid information
      subroutine get_grid    ! from GRID NetCDF file.
      implicit none
      character*1 char1
      integer ierr, ncid, varid, lstr, lenstr, checkdims, nf_fread
      real zob_temp
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "grid.h"
# include "netcdf.inc"
# ifdef MPI
      include 'mpif.h'
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif

!
! Open grid netCDF file for reading. Check that dimensions in that
! file are consistent with the model, then read all necessary
! variables.
!
      lstr=lenstr(grdname)
      ierr=nf_open(grdname(1:lstr), nf_nowrite, ncid)
      if (ierr .eq. nf_noerr) then
        ierr=checkdims (ncid, grdname, lstr, varid)
        if (ierr. ne. nf_noerr) goto 99
      else
        write(stdout,'(/3(1x,A)/)') 'GET_GRID ERROR: Cannot open',
     &                       'input NetCDF file:', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
!
! Logical switch for spherical grid configuration:
!
      ierr=nf_inq_varid (ncid, 'spherical', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_get_var1_text (ncid, varid, 1, char1)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'spherical', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'spherical', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif

      if (char1.eq.'t' .or. char1.eq.'T') then
# ifdef SPHERICAL
        MPI_master_only write(stdout,'(/1x,A/)')
     &                          'Spherical grid detected.'
# else
        write(stdout,'(/1x,A/12x,A/)')
     &      'GET_GRID - ERROR: Spherical grid detected, but',
     &               'model SPHERICAL CPP-switch is not set.'
        goto 99                                            !--> ERROR
# endif
      endif
!
! Physical dimensions of the basin in XI- and ETA-directions:
!
      ierr=nf_inq_varid (ncid, 'xl', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_get_var1_FTYPE (ncid, varid, 1, xl)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'xl', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'xl', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif

      ierr=nf_inq_varid (ncid, 'el', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_get_var1_FTYPE (ncid, varid, 1, el)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'el', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'el', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
!
! Read in grid arrays.
!===== == ==== =======
! Bathymetry:
!
      ierr=nf_inq_varid (ncid, 'h', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread(h (START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'h', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'h', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
! REMOVE WET DRY    
!     h(:,:)=max(h(:,:),10.0)
!
! Coriolis parameter.
!
      ierr=nf_inq_varid (ncid, 'f', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (f(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'f', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'f', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
!
! Coordinate transfomation metrics (m,n) associated with the
! differential distances in XI and ETA.
!
      ierr=nf_inq_varid (ncid, 'pm', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (pm(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'pm', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'pm', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif

      ierr=nf_inq_varid (ncid, 'pn', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (pn(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'pn', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'pn', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
!
! Coordinates (lon,lat [degrees]) or (x,y [meters]) at RHO-points.
!
# ifdef SPHERICAL
      ierr=nf_inq_varid (ncid, 'lon_rho', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (lonr(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'lon_rho', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'lon_rho', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif

      ierr=nf_inq_varid (ncid, 'lat_rho', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (latr(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'lat_rho', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'lat_rho', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
!--
      ierr=nf_inq_varid (ncid, 'lon_u', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (lonu(START_2D_ARRAY), ncid, varid, 0, u2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'lon_u', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'lon_u', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif

      ierr=nf_inq_varid (ncid, 'lat_u', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (latu(START_2D_ARRAY), ncid, varid, 0, u2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'lat_u', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'lat_u', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
!--
      ierr=nf_inq_varid (ncid, 'lon_v', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (lonv(START_2D_ARRAY), ncid, varid, 0, v2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'lon_v', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'lon_v', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif

      ierr=nf_inq_varid (ncid, 'lat_v', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (latv(START_2D_ARRAY), ncid, varid, 0, v2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'lat_v', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'lat_v', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
!--

# else
      ierr=nf_inq_varid (ncid, 'x_rho', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (xr(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'x_rho', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'x_rho', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif

      ierr=nf_inq_varid (ncid, 'y_rho', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (yr(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'y_rho', grdname(1:lstr)
          goto 99                                           !--> ERROR
        endif
      else
        write(stdout,1) 'y_rho', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
# endif
# ifdef CURVGRID
!
! Angle (radians) between XI-axis and EAST at RHO-points.
!
      ierr=nf_inq_varid (ncid, 'angle', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (angler(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'angle', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'angle', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
# endif
# ifdef MASKING
!
! Mask at RHO-points.
!
      ierr=nf_inq_varid (ncid, 'mask_rho', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread(rmask(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'mask_rho', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'mask_rho', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif

# endif /* MASKING */

# ifdef REDUC_SECTION
      ureduc(:,:)=1.0
      vreduc(:,:)=1.0
      ierr=nf_inq_varid (ncid, 'reduc_u', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (ureduc(START_2D_ARRAY), ncid, varid, 0, u2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'reduc_u', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'reduc_u', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
      ureduc(:,:)=min(ureduc(:,:),1.0)

      ierr=nf_inq_varid (ncid, 'reduc_v', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (vreduc(START_2D_ARRAY), ncid, varid, 0, v2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'reduc_v', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'reduc_v', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
      vreduc(:,:)=min(vreduc(:,:),1.0)
# endif
# ifdef Z0B_VAR
      zob(:,:)=0.0
      ierr=nf_inq_varid (ncid, 'z0b', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (zob(START_2D_ARRAY), ncid, varid, 0, v2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'zob', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'zob', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
      ! write(stdout,*)'zob',
     !&  minval(zob),maxval(zob),minloc(zob),maxloc(zob)
       zob(:,:)=max(zob(:,:),1e-7)
       zob_temp=maxval(zob)
#ifdef MPI
          if (NNODES.gt.1) then         ! Perform global max
            call MPI_allreduce(zob_temp,zobt,1,
     &              MPI_DOUBLE_PRECISION, MPI_MAX,MPI_COMM_WORLD,ierr)
          endif
# endif
# else 
       zob(:,:)=zobt 
# endif


      ierr=nf_close(ncid)
      return                                      !--> NORMAL RETURN

   1  format(/1x, 'GET_GRID - unable to find grid variable:',
     &           1x, A / 12x, 'in grid netCDF file:', 1x, A/)
   2  format(/1x, 'GET_GRID - error while reading variable:',
     &           1x, A / 12x, 'in grid netCDF file:', 1x, A/)
   3  format ('bornes h:',2(1x,F10.7),1x,4(i4,1x))
   4  format (A,2(1x,F10.7),1x,4(i4,1x))
  99  may_day_flag=2
      return                                              !--> ERROR
      end
#else
      subroutine get_grid_empty
      return
      end
#endif /* !ANA_GRID */
