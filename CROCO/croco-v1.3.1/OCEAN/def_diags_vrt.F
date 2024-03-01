! $Id: def_diags_vrt.F 1468 2014-02-10 12:58:35Z rblod $
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
#ifndef AVRH
# include "cppdefs.h"
#endif
#ifdef DIAGNOSTICS_VRT
#ifndef AVRH
      subroutine def_diags_vrt (ncid, total_rec, ierr)
#else
      subroutine def_diags_vrt_avg(ncid, total_rec, ierr)
#endif
!
!---------------------------------------------------------------
! Create diag data NetCDF file:
! Define its variables, dimensions, and attributes
!---------------------------------------------------------------
!
      implicit none
# include "param.h"
# include "mixing.h"
# include "ncscrum.h"
# include "scalars.h"
# include "strings.h"
# include "diags_vrt.h"
# include "netcdf.inc"
#ifdef NC4PAR
# include "mpi_cpl.h"
      include 'mpif.h'
#endif
      logical create_new_file, res
      integer ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim
     &      , r2dgrd(3),u2dgrd(3),v2dgrd(3),auxil(2),checkdims 
#ifdef SOLVE3D
     &      , r3dgrd(4),  u3dgrd(4), v3dgrd(4),  w3dgrd(4), itrc
     &      , p2dgrd(3), p3dgrd(4), pw3dgrd(4)
#endif
#ifdef NC4PAR
     &      , csize,cmode
#endif

      character*60 text

!
#ifndef AVRH

# define ncname diags_vrtname
# define rec_per_file nrpfdiags_vrt
# define wrt wrtdiags_vrt
# define vidTime diags_vrtTime
# define vidTime2 diags_vrtTime2
# define vidTstep diags_vrtTstep

#  define viddiags_vrtXadv diags_vrtXadv
#  define viddiags_vrtYadv diags_vrtYadv
#  define viddiags_vrtHdiff diags_vrtHdiff
#  define viddiags_vrtCor diags_vrtCor
#  define viddiags_vrtPrsgrd diags_vrtPrsgrd
#  define viddiags_vrtHmix diags_vrtHmix
#  define viddiags_vrtVmix diags_vrtVmix
#  define viddiags_vrtrate diags_vrtrate
#  define viddiags_vrtVmix2 diags_vrtVmix2
#  define viddiags_vrtWind diags_vrtWind
#  define viddiags_vrtDrag diags_vrtDrag
# if defined DIAGNOSTICS_BARO
#  define viddiags_vrtBaro diags_vrtBaro
# endif
# if defined M3FAST
#  define viddiags_vrtfast diags_vrtfast
# endif
#else

# define ncname diags_vrtname_avg
# define rec_per_file nrpfdiags_vrt_avg
# define wrt wrtdiags_vrt_avg
# define vidTime diags_vrtTime_avg
# define vidTime2 diags_vrtTime2_avg
# define vidTstep diags_vrtTstep_avg

#  define viddiags_vrtXadv diags_vrtXadv_avg
#  define viddiags_vrtYadv diags_vrtYadv_avg
#  define viddiags_vrtHdiff diags_vrtHdiff_avg
#  define viddiags_vrtCor diags_vrtCor_avg
#  define viddiags_vrtPrsgrd diags_vrtPrsgrd_avg
#  define viddiags_vrtHmix diags_vrtHmix_avg
#  define viddiags_vrtVmix diags_vrtVmix_avg
#  define viddiags_vrtrate diags_vrtrate_avg
#  define viddiags_vrtVmix2 diags_vrtVmix2_avg
#  define viddiags_vrtWind diags_vrtWind_avg
#  define viddiags_vrtDrag diags_vrtDrag_avg
# if defined DIAGNOSTICS_BARO
#  define viddiags_vrtBaro diags_vrtBaro_avg
# endif
# if defined M3FAST
#  define viddiags_vrtfast diags_vrtfast_avg
# endif
#endif

!
      if (may_day_flag.ne.0) return      !-->  EXIT
!
      ierr=0
      lstr=lenstr(ncname)
      if (rec_per_file.gt.0) then
        lvar=total_rec-(1+mod(total_rec-1, rec_per_file))
        call insert_time_index (ncname, lstr, lvar, ierr)
        if (ierr .ne. 0) goto 99
      endif
!
!  Create a new diagnostics data file.
!-------------------------------------
!
#ifdef DO_NOT_OVERWRITE
      if (lvar.eq.total_rec-1) then
# ifndef AVRH
        create_new_file=ldefdiags_vrt
# else
        create_new_file=ldefdiags_vrt_avg
# endif
      else
# ifndef AVRH
        create_new_file=.false.
# else
        create_new_file=.false.
# endif
      endif
#else
# ifndef AVRH
      create_new_file=ldefdiags_vrt
# else
      create_new_file=ldefdiags_vrt_avg
# endif
#endif

!
      if (ncid.ne.-1) create_new_file=.false.
#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      if (mynode.gt.0) create_new_file=.false.
#endif
!
 10   if (create_new_file) then
        lstr=lenstr(ncname)
#ifndef NC4PAR
        ierr=nf_create(ncname(1:lstr),NF_CLOBBER,ncid)
#else
        cmode = ior(nf_netcdf4,nf_classic_model)
        cmode = ior(cmode, nf_mpiio)
        csize=xi_rho*eta_rho/NNODES
        MPI_master_only write(stdout,*)'CREATE DIAG NC4 PARALLEL FILE'
        ierr=nf_create_par(ncname(1:lstr),cmode,
     &        MPI_COMM_WORLD,MPI_INFO_NULL,ncid)
#endif
        if (ierr.ne.nf_noerr) then
          write(stdout,11) ncname(1:lstr)
          may_day_flag=3
          return                         !-->  EXIT
        endif
!
! Put global attributes.
! --- ------ -----------
!
        call put_global_atts (ncid, ierr)
        if (ierr.ne.nf_noerr) then
          write(stdout,11) ncname(1:lstr)
          may_day_flag=3
          return                         !-->  EXIT
        endif
!
! Define dimensions of staggered fields.
! ------ ---------- -- --------- -------
!
        ierr=nf_def_dim (ncid, 'xi_rho',   xi_rho,   r2dgrd(1))
        ierr=nf_def_dim (ncid, 'xi_u',     xi_u,     u2dgrd(1))
        ierr=nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr=nf_def_dim (ncid, 'eta_v',    eta_v,    v2dgrd(2))
#ifdef SOLVE3D
        ierr=nf_def_dim (ncid, 's_rho',    N,        r3dgrd(3))
        ierr=nf_def_dim (ncid, 's_w',      N+1,      w3dgrd(3))
#endif
        ierr=nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr=nf_def_dim (ncid, 'auxil',    4,        auxil(1))
!
        auxil(2)=timedim
 
        r2dgrd(3)=timedim           ! Free surface
 
        u2dgrd(2)=r2dgrd(2)         ! 2D UBAR-type
        u2dgrd(3)=timedim
 
        v2dgrd(1)=r2dgrd(1)         ! 2D VBAR-type
        v2dgrd(3)=timedim

        p2dgrd(1)=u2dgrd(1)         ! 
        p2dgrd(2)=v2dgrd(2)         ! 2D vort-type
        p2dgrd(3)=timedim

#ifdef SOLVE3D

        r3dgrd(1)=r2dgrd(1)         !
        r3dgrd(2)=r2dgrd(2)         ! 3D RHO-type
        r3dgrd(4)=timedim           !
 
        u3dgrd(1)=u2dgrd(1)         !
        u3dgrd(2)=r2dgrd(2)         ! 3D U-type
        u3dgrd(3)=r3dgrd(3)         !
        u3dgrd(4)=timedim
 
        v3dgrd(1)=r2dgrd(1)         !
        v3dgrd(2)=v2dgrd(2)         ! 3D V-type
        v3dgrd(3)=r3dgrd(3)         !
        v3dgrd(4)=timedim
 
        w3dgrd(1)=r2dgrd(1)         !
        w3dgrd(2)=r2dgrd(2)         ! 3D W-type
        w3dgrd(4)=timedim           !

        p3dgrd(1)=u2dgrd(1)         !
        p3dgrd(2)=v2dgrd(2)         ! 3D vort-type
        p3dgrd(3)=r2dgrd(3)         ! 
        p3dgrd(4)=timedim           !

        pw3dgrd(1)=u2dgrd(1)         !
        pw3dgrd(2)=v2dgrd(2)         ! 3D vortW-type
        pw3dgrd(3)=w3dgrd(3)         ! 
        pw3dgrd(4)=timedim         

#endif
!
!
! Define running parameters :
!--------------------------
!
! Define variables (time step number and time record indices)
! and their attributes
!
        ierr=nf_def_var (ncid, 'time_step', nf_int, 2, auxil,
     &                                                 vidTstep)
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,vidTstep,nf_collective)
#endif
        ierr=nf_put_att_text (ncid, vidTstep, 'long_name', 48,
     &       'time step and record numbers from initialization')
!
! Time.
!
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, vidTime)
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,vidTime,nf_collective)
#endif
!
#ifndef AVRH
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', 
     &       lvar, vname(2,indxTime)(1:lvar))
#else
        text='avg'/ /vname(2,indxTime)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', 
     &       lvar, text(1:lvar))
#endif
!
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', lvar,
     &                                vname(2,indxTime)(1:lvar))
!
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))
!
        lvar=lenstr(vname(4,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'field',  lvar,
     &                                vname(4,indxTime)(1:lvar))

        call nf_add_attribute(ncid, vidTime, indxTime, 5, 
     &       NF_FOUT, ierr)
!
! Time2.
!
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_def_var (ncid, vname(1,indxTime2)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, vidTime2)
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,vidTime2,nf_collective)
#endif
!
#ifndef AVRH
        lvar=lenstr(vname(2,indxTime2))
        ierr=nf_put_att_text (ncid, vidTime2, 'long_name', 
     &       lvar, vname(2,indxTime)(1:lvar))
#else
        text='avg'/ /vname(2,indxTime2)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidTime2, 'long_name', 
     &       lvar, text(1:lvar))
#endif
!
        lvar=lenstr(vname(2,indxTime2))
        ierr=nf_put_att_text (ncid, vidTime2, 'long_name', lvar,
     &                                vname(2,indxTime2)(1:lvar))
!
        lvar=lenstr(vname(3,indxTime2))
        ierr=nf_put_att_text (ncid, vidTime2, 'units',  lvar,
     &                                vname(3,indxTime2)(1:lvar))
!
        lvar=lenstr(vname(4,indxTime2))
        ierr=nf_put_att_text (ncid, vidTime2, 'field',  lvar,
     &                                vname(4,indxTime2)(1:lvar))

        call nf_add_attribute(ncid, vidTime2, indxTime2, 5, 
     &       NF_FOUT, ierr)

! Momentum diagnostics variables.
! -------- ----------- --------

!
          itrc=1
!
          if (wrt(itrc)) then

!------------------------------
!
! indxvrtXadv
!          
          lvar=lenstr(vname(1,indxvrtXadv+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxvrtXadv+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, p2dgrd, viddiags_vrtXadv(itrc))
#ifdef NC4PAR
          ierr=nf_var_par_access(ncid,viddiags_vrtXadv(itrc),nf_collective)
#endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxvrtXadv+itrc-1)
# else
          text=vname(2,indxvrtXadv+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_vrtXadv(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxvrtXadv+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_vrtXadv(itrc), 'units', lvar,
     &                          vname(3,indxvrtXadv+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxvrtXadv+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_vrtXadv(itrc), 'field',
     &                      lvar, vname(4,indxvrtXadv+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, viddiags_vrtXadv(itrc), 
     &         indxvrtXadv+itrc-1, 5, NF_FOUT, ierr)
!------------------------------
!
! indxvrtYadv
!          
          lvar=lenstr(vname(1,indxvrtYadv+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxvrtYadv+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, p2dgrd, viddiags_vrtYadv(itrc))
#ifdef NC4PAR
          ierr=nf_var_par_access(ncid,viddiags_vrtYadv(itrc),nf_collective)
#endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxvrtYadv+itrc-1)
# else
          text=vname(2,indxvrtYadv+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_vrtYadv(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxvrtYadv+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_vrtYadv(itrc), 'units', lvar,
     &                          vname(3,indxvrtYadv+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxvrtYadv+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_vrtYadv(itrc), 'field',
     &                      lvar, vname(4,indxvrtYadv+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, viddiags_vrtYadv(itrc), 
     &        indxvrtYadv+itrc-1, 5, NF_FOUT, ierr)
!------------------------------
!
! indxvrtHdiff
!          
          lvar=lenstr(vname(1,indxvrtHdiff+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxvrtHdiff+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, p2dgrd, viddiags_vrtHdiff(itrc))
#ifdef NC4PAR
          ierr=nf_var_par_access(ncid,viddiags_vrtHdiff(itrc),nf_collective)
#endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxvrtHdiff+itrc-1)
# else
          text=vname(2,indxvrtHdiff+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_vrtHdiff(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxvrtHdiff+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_vrtHdiff(itrc), 'units', lvar,
     &                          vname(3,indxvrtHdiff+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxvrtHdiff+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_vrtHdiff(itrc), 'field',
     &                      lvar, vname(4,indxvrtHdiff+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, viddiags_vrtHdiff(itrc), 
     &            indxvrtHdiff+itrc-1, 5,NF_FOUT, ierr)
!------------------------------!
! indxvrtCor
!          
          lvar=lenstr(vname(1,indxvrtCor+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxvrtCor+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, p2dgrd, viddiags_vrtCor(itrc))
#ifdef NC4PAR
          ierr=nf_var_par_access(ncid,viddiags_vrtCor(itrc),nf_collective)
#endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxvrtCor+itrc-1)
# else
          text=vname(2,indxvrtCor+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_vrtCor(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxvrtCor+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_vrtCor(itrc), 'units', lvar,
     &                          vname(3,indxvrtCor+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxvrtCor+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_vrtCor(itrc), 'field',
     &                      lvar, vname(4,indxvrtCor+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, viddiags_vrtCor(itrc), 
     &            indxvrtCor+itrc-1, 5, NF_FOUT, ierr)

!------------------------------
!
! indxvrtPrsgrd
!          
          lvar=lenstr(vname(1,indxvrtPrsgrd+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxvrtPrsgrd+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, p2dgrd, viddiags_vrtPrsgrd(itrc))
#ifdef NC4PAR
          ierr=nf_var_par_access(ncid,viddiags_vrtPrsgrd(itrc),nf_collective)
#endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxvrtPrsgrd+itrc-1)
# else
          text=vname(2,indxvrtPrsgrd+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_vrtPrsgrd(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxvrtPrsgrd+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_vrtPrsgrd(itrc), 'units', lvar,
     &                          vname(3,indxvrtPrsgrd+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxvrtPrsgrd+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_vrtPrsgrd(itrc), 'field',
     &                      lvar, vname(4,indxvrtPrsgrd+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, viddiags_vrtPrsgrd(itrc), 
     &           indxvrtPrsgrd+itrc-1, 5,NF_FOUT, ierr)
!------------------------------
!
! indxvrtHmix
!          
          lvar=lenstr(vname(1,indxvrtHmix+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxvrtHmix+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, p2dgrd, viddiags_vrtHmix(itrc))
#ifdef NC4PAR
          ierr=nf_var_par_access(ncid,viddiags_vrtHmix(itrc),nf_collective)
#endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxvrtHmix+itrc-1)
# else
          text=vname(2,indxvrtHmix+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_vrtHmix(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxvrtHmix+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_vrtHmix(itrc), 'units', lvar,
     &                          vname(3,indxvrtHmix+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxvrtHmix+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_vrtHmix(itrc), 'field',
     &                      lvar, vname(4,indxvrtHmix+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, viddiags_vrtHmix(itrc),
     &            indxvrtHmix+itrc-1,5, NF_FOUT, ierr)
!------------------------------
!
! indxvrtVmix
!          
          lvar=lenstr(vname(1,indxvrtVmix+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxvrtVmix+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, p2dgrd, viddiags_vrtVmix(itrc))
#ifdef NC4PAR
          ierr=nf_var_par_access(ncid,viddiags_vrtVmix(itrc),nf_collective)
#endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxvrtVmix+itrc-1)
# else
          text=vname(2,indxvrtVmix+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_vrtVmix(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxvrtVmix+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_vrtVmix(itrc), 'units', lvar,
     &                          vname(3,indxvrtVmix+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxvrtVmix+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_vrtVmix(itrc), 'field',
     &                      lvar, vname(4,indxvrtVmix+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, viddiags_vrtVmix(itrc), 
     &            indxvrtVmix+itrc-1, 5,NF_FOUT, ierr)
!------------------------------
!
! indxvrtrate
!          
          lvar=lenstr(vname(1,indxvrtrate+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxvrtrate+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, p2dgrd, viddiags_vrtrate(itrc))
#ifdef NC4PAR
          ierr=nf_var_par_access(ncid,viddiags_vrtrate(itrc),nf_collective)
#endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxvrtrate+itrc-1)
# else
          text=vname(2,indxvrtrate+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_vrtrate(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxvrtrate+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_vrtrate(itrc), 'units', lvar,
     &                          vname(3,indxvrtrate+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxvrtrate+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_vrtrate(itrc), 'field',
     &                      lvar, vname(4,indxvrtrate+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, viddiags_vrtrate(itrc), 
     &             indxvrtrate+itrc-1, 5,NF_FOUT, ierr)

!------------------------------
!
! indxvrtVmix2
!          
          lvar=lenstr(vname(1,indxvrtVmix2+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxvrtVmix2+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, p2dgrd, viddiags_vrtVmix2(itrc))
#ifdef NC4PAR
          ierr=nf_var_par_access(ncid,viddiags_vrtVmix2(itrc),nf_collective)
#endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxvrtVmix2+itrc-1)
# else
          text=vname(2,indxvrtVmix2+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_vrtVmix2(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxvrtVmix2+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_vrtVmix2(itrc), 'units', lvar,
     &                          vname(3,indxvrtVmix2+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxvrtVmix2+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_vrtVmix2(itrc),'field',
     &                    lvar, vname(4,indxvrtVmix2+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, viddiags_vrtVmix2(itrc), 
     &               indxvrtVmix2+itrc-1, 5,NF_FOUT, ierr)
!


!------------------------------
!
! indxvrtWind
!          
          lvar=lenstr(vname(1,indxvrtWind+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxvrtWind+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, p2dgrd, viddiags_vrtWind(itrc))
#ifdef NC4PAR
          ierr=nf_var_par_access(ncid,viddiags_vrtWind(itrc),nf_collective)
#endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxvrtWind+itrc-1)
# else
          text=vname(2,indxvrtWind+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_vrtWind(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxvrtWind+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_vrtWind(itrc), 'units', lvar,
     &                          vname(3,indxvrtWind+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxvrtWind+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_vrtWind(itrc),'field',
     &                    lvar, vname(4,indxvrtWind+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, viddiags_vrtWind(itrc), 
     &           indxvrtWind+itrc-1, 5,NF_FOUT, ierr)

!------------------------------
!
! indxvrtDrag
!          
          lvar=lenstr(vname(1,indxvrtDrag+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxvrtDrag+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, p2dgrd, viddiags_vrtDrag(itrc))
#ifdef NC4PAR
          ierr=nf_var_par_access(ncid,viddiags_vrtDrag(itrc),nf_collective)
#endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxvrtDrag+itrc-1)
# else
          text=vname(2,indxvrtDrag+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_vrtDrag(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxvrtDrag+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_vrtDrag(itrc), 'units', lvar,
     &                          vname(3,indxvrtDrag+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxvrtDrag+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_vrtWind(itrc),'field',
     &                    lvar, vname(4,indxvrtDrag+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, viddiags_vrtWind(itrc), 
     &              indxvrtDrag+itrc-1, 5, NF_FOUT, ierr)
        endif
# if defined DIAGNOSTICS_BARO
!------------------------------
!
! indxvrtBaro
!
          lvar=lenstr(vname(1,indxvrtBaro+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxvrtBaro+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, p2dgrd, viddiags_vrtBaro(itrc))
#ifdef NC4PAR
           ierr=nf_var_par_access(ncid,viddiags_vrtBaro(itrc),nf_collective)
#endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxvrtBaro+itrc-1)
# else
          text=vname(2,indxvrtBaro+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_vrtBaro(itrc), 'long_name',
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxvrtBaro+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_vrtBaro(itrc), 'units', lvar,
     &                          vname(3,indxvrtBaro+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxvrtBaro+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_vrtBaro(itrc), 'field',
     &                      lvar, vname(4,indxvrtBaro+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, viddiags_vrtBaro(itrc),
     &           indxvrtBaro+itrc-1, 5,NF_FOUT, ierr)
# endif

# if defined M3FAST
!------------------------------
!
! indxvrtfast
!
          lvar=lenstr(vname(1,indxvrtfast+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxvrtfast+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, p2dgrd, viddiags_vrtfast(itrc))
#ifdef NC4PAR
           ierr=nf_var_par_access(ncid,viddiags_vrtfast(itrc),nf_collective)
#endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxvrtfast+itrc-1)
# else
          text=vname(2,indxvrtfast+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_vrtfast(itrc), 'long_name',
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxvrtfast+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_vrtfast(itrc), 'units', lvar,
     &                          vname(3,indxvrtfast+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxvrtfast+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_vrtfast(itrc), 'field',
     &                      lvar, vname(4,indxvrtfast+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, viddiags_vrtfast(itrc),
     &           indxvrtfast+itrc-1, 5,NF_FOUT, ierr)
# endif

!
! Leave definition mode.
! ----- ---------- -----
!
        ierr=nf_enddef(ncid)
        MPI_master_only write(stdout,'(6x,4A,1x,A,i4)')
#ifdef AVRH
     &        'DEF_DIAGS_VRT_AVG - Created ',
#else      
     &        'DEF_DIAGS_VRT - Created ',
#endif /*AVRH*/
     &                'new netCDF file ''', 
     &                 ncname(1:lstr), '''.'
     &                 MYID
!
! Open an existing file and prepare for appending data.
! ==== == ======== ==== === ======= === ========= =====
! Inquire about the dimensions and variables. Check for
! consistency with model dimensions. In the case when file 
! is rejected (whether it cannot be opened, or something
! is wrong with its dimensions) create a new file.2
!
! After that verify that all necessary variables are already
! defined there and find their netCDF IDs.
!
      elseif (ncid.eq.-1) then
#ifndef NC4PAR
        ierr=nf_open (ncname(1:lstr), nf_write, ncid)
#else
        MPI_master_only write(stdout,*) "Open file in parallel"
        ierr=nf_open_par (ncname(1:lstr), IOR(nf_write, nf_mpiio),
     &     MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
#endif
        if (ierr.eq. nf_noerr) then
          ierr=checkdims (ncid, ncname, lstr, rec) 
          if (ierr .eq. nf_noerr) then
            if (rec_per_file.eq.0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, rec_per_file))
            endif
            if (ierr.gt.0) then
              MPI_master_only write( stdout,
     &                 '(/1x,A,I5,1x,A/8x,3A,I5,/8x,A,I5,1x,A/)'
     &           ) 'WARNING: def_diags_vrt: Actual number of records',
     &               rec,  'in netCDF file',  '''',  ncname(1:lstr),
     &             ''' exceeds the record number from restart data',
     &             rec+1-ierr,'/', total_rec,', restart is assumed.'
              rec=rec-ierr
            elseif (rec_per_file.eq.0) then
              total_rec=rec+1           ! <-- set to the next record
#if defined MPI & !defined PARALLEL_FILES
              if (mynode.gt.0) total_rec=total_rec-1
#endif
            endif
            ierr=nf_noerr
          endif
        endif
!
        if (ierr. ne. nf_noerr) then
#if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
          if (mynode.eq.0) then
            create_new_file=.true.
            goto 10
          else
        write(stdout,'(/1x,4A,2x,A,I4/)') 'def_his/avg ERROR: ',
     &         'Cannot open file ''', ncname(1:lstr), '''.'
     &                   MYID
            goto 99                                       !--> ERROR
          endif
#else
          create_new_file=.true.
          goto 10
#endif
        endif
!
! Find netCDF IDs of evolving model variables:
! ---- ------ --- -- -------- ----- ----------
!
! Time step indices:
!         
        ierr=nf_inq_varid (ncid, 'time_step', vidTstep)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) 'time_step', ncname(1:lstr)
          goto 99                                         !--> ERROR
        endif 
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,vidTstep,nf_collective)
#endif
!
! Time.
!
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),vidTime)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime)(1:lvar), ncname(1:lstr)
          goto 99                                         !--> ERROR 
        endif
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,vidTime,nf_collective)
#endif
!
! Time2.
!
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_inq_varid (ncid,vname(1,indxTime2)(1:lvar),vidTime2)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime2)(1:lvar), ncname(1:lstr)
          goto 99                                         !--> ERROR 
        endif
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,vidTime2,nf_collective)
#endif

! momentum flux diagnostics variables :
! -------- ---- ----------- ----------
!
!
          itrc=1

          if (wrt(itrc)) then
!          
! indxvrtXadv
!

         lvar=lenstr(vname(1,indxvrtXadv+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxvrtXadv+itrc-1)(1:lvar),
     &                      viddiags_vrtXadv(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxvrtXadv+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
#ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_vrtXadv(itrc),
     &                          nf_collective)
#endif
!
! indxvrtYadv
!

         lvar=lenstr(vname(1,indxvrtYadv+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxvrtYadv+itrc-1)(1:lvar),
     &                      viddiags_vrtYadv(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxvrtYadv+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif      
#ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_vrtYadv(itrc),
     &                          nf_collective)
#endif
!
! indxvrtHdiff
!
         lvar=lenstr(vname(1,indxvrtHdiff+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxvrtHdiff+itrc-1)(1:lvar),
     &                      viddiags_vrtHdiff(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxvrtHdiff+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif
#ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_vrtHdiff(itrc),
     &                          nf_collective)
#endif
!
! indxvrtCor
!
         lvar=lenstr(vname(1,indxvrtCor+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxvrtCor+itrc-1)(1:lvar),
     &                      viddiags_vrtCor(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxvrtCor+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif
#ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_vrtCor(itrc),
     &                          nf_collective)
#endif
!
! indxvrtPrsgrd
!
         lvar=lenstr(vname(1,indxvrtPrsgrd+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxvrtPrsgrd+itrc-1)(1:lvar),
     &                      viddiags_vrtPrsgrd(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxvrtPrsgrd+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif
#ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_vrtPrsgrd(itrc),
     &                          nf_collective)
#endif
!
! indxvrtHmix
!
         lvar=lenstr(vname(1,indxvrtHmix+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxvrtHmix+itrc-1)(1:lvar),
     &                      viddiags_vrtHmix(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxvrtHmix+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif
#ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_vrtHmix(itrc),
     &                          nf_collective)
#endif
!
! indxvrtVmix
!
         lvar=lenstr(vname(1,indxvrtVmix+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxvrtVmix+itrc-1)(1:lvar),
     &                      viddiags_vrtVmix(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxvrtVmix+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
#ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_vrtVmix(itrc),
     &                          nf_collective)
#endif
!
! indxvrtrate
!
         lvar=lenstr(vname(1,indxvrtrate+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxvrtrate+itrc-1)(1:lvar),
     &                      viddiags_vrtrate(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxvrtrate+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
#ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_vrtrate(itrc),
     &                          nf_collective)
#endif

!
! indxvrtVmix2
!
         lvar=lenstr(vname(1,indxvrtVmix2+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxvrtVmix2+itrc-1)(1:lvar),
     &                      viddiags_vrtVmix2(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxvrtVmix2+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
#ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_vrtVmix2(itrc),
     &                          nf_collective)
#endif




!
! indxvrtWind
!
         lvar=lenstr(vname(1,indxvrtWind+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxvrtWind+itrc-1)(1:lvar),
     &                      viddiags_vrtWind(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxvrtWind+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
#ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_vrtWind(itrc),
     &                          nf_collective)
#endif



!
! indxvrtDrag
!
         lvar=lenstr(vname(1,indxvrtDrag+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxvrtDrag+itrc-1)(1:lvar),
     &                      viddiags_vrtDrag(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxvrtDrag+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
#ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_vrtDrag(itrc),
     &                          nf_collective)
#endif
# if defined DIAGNOSTICS_BARO
!
! indxvrtBaro
!
         lvar=lenstr(vname(1,indxvrtBaro+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxvrtBaro+itrc-1)(1:lvar),
     &                      viddiags_vrtBaro(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxvrtBaro+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif
#  ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_vrtBaro(itrc),
     &                          nf_collective)
#  endif
# endif

# if defined M3FAST
!
! indxvrtfast
!
         lvar=lenstr(vname(1,indxvrtfast+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxvrtfast+itrc-1)(1:lvar),
     &                      viddiags_vrtfast(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxvrtfast+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif
#  ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_vrtfast(itrc),
     &                          nf_collective)
#  endif
# endif

!
       endif

!
        MPI_master_only write(*,'(6x,2A,i4,1x,A,i4)') 
     &                     'def_diags_vrt: -- Opened ',
     &                     'existing file  from record =', rec
     &                      MYID
!
        MPI_master_only write(*,'(6x,2A,i4,1x,A,i4)') 
     &                     'def_diags_vrt: -- Opened ',
     &                     'existing file  from record =', rec
     &                      MYID
!
#if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      else
        ierr=nf_open (ncname(1:lstr), nf_write, ncid)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,'(/1x,4A,2x,A,I4/)')
     &                'def_diags_vrt: ERROR: ',
     &                'Cannot open file ''', ncname(1:lstr), '''.'
     &                 MYID
          goto 99                                         !--> ERROR
        endif
#endif
!
      endif             !<--  create_new_file  
!
      ierr=nf_set_fill (ncid, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
        MPI_master_only write(*,'(6x,2A,i4,1x,A,i4)')
     &    'def_diags_vrt ERROR: Cannot ',
     &    'switch to ''nf_nofill'' more; netCDF error code =', ierr 
      endif
   1  format(/1x,'def_diags_vrt ERROR: Cannot find variable ''',
     &                   A, ''' in netCDF file ''', A, '''.'/)

!
  11  format(/' def_diags_vrt - unable to create diag file: ',a)
  20  format(/' def_diags_vrt - error while writing variable: ',a,
     &        /,15x,'into diag  file: ',a)
!
  99  return
      end
!
#undef ncname
#undef rec_per_file
#undef wrt
#undef vidTime
#undef vidTime2
#undef vidTstep
!
#undef viddiags_vrtXadv
#undef viddiags_vrtYadv
#undef viddiags_vrtHdiff
#undef viddiags_vrtCor
#undef viddiags_vrtPrsgrd
#undef viddiags_vrtHmix
#undef viddiags_vrtVmix
#undef viddiags_vrtrate
#undef viddiags_vrtVmix2
#undef viddiags_vrtWind
#undef viddiags_vrtDrag
# if defined DIAGNOSTICS_BARO
#undef viddiags_vrtBaro
# endif
# if defined M3FAST
#undef viddiags_vrtfast
# endif
!
# ifdef AVERAGES
#  ifndef AVRH
#   define AVRH
#   include "def_diags_vrt.F"
#   endif
#  endif /* AVERAGES */
!
#else /*DIAGNOSTICS_VRT */
      subroutine def_diags_vrt_empty()
      return
      end
#endif /* DIAGNOSTICS_VRT */

