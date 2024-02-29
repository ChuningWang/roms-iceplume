! $Id: def_diags_eddy.F 1468 2014-02-10 12:58:35Z rblod $
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

#if defined DIAGNOSTICS_EDDY && ! defined XIOS

# ifndef AVRH
      subroutine def_diags_eddy (ncid, total_rec, ierr)
# else
      subroutine def_diags_eddy_avg(ncid, total_rec, ierr)
# endif
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
# include "diags_eddy.h"
# include "netcdf.inc"
      logical create_new_file, res
      integer ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim
     &      , r2dgrd(3),u2dgrd(3),v2dgrd(3),auxil(2),checkdims 
# ifdef SOLVE3D
     &      , r3dgrd(4),  u3dgrd(4), v3dgrd(4),  w3dgrd(4), itrc
     &      , p2dgrd(3), p3dgrd(4), pw3dgrd(4)
# endif
# ifdef NC4PAR
     &      , csize,cmode
#  include "mpi_cpl.h"
      include 'mpif.h'
# endif
      character*60 text
!
# ifndef AVRH

#  define ncname diags_eddyname
#  define rec_per_file nrpfdiags_eddy
#  define wrt wrtdiags_eddy
#  define vidTime diags_eddyTime
#  define vidTime2 diags_eddyTime2
#  define vidTstep diags_eddyTstep
#  define viddiags_zz diags_eddyzz
#  define viddiags_uu diags_eddyuu
#  define viddiags_vv diags_eddyvv
#  define viddiags_uv diags_eddyuv
#  define viddiags_ub diags_eddyub
#  define viddiags_vb diags_eddyvb
#  define viddiags_uw diags_eddyuw
#  define viddiags_vw diags_eddyvw
#  define viddiags_wb diags_eddywb
#  define viddiags_ubu diags_eddyubu
#  define viddiags_vbv diags_eddyvbv
#  define viddiags_usu diags_eddyusu
#  define viddiags_vsv diags_eddyvsv
#  define viddiags_ugsu diags_eddyugsu
#  define viddiags_vgsv diags_eddyvgsv

# else

#  define ncname diags_eddyname_avg
#  define rec_per_file nrpfdiags_eddy_avg
#  define wrt wrtdiags_eddy_avg
#  define vidTime diags_eddyTime_avg
#  define vidTime2 diags_eddyTime2_avg
#  define vidTstep diags_eddyTstep_avg
#  define viddiags_zz diags_eddyzz_avg
#  define viddiags_uu diags_eddyuu_avg
#  define viddiags_vv diags_eddyvv_avg
#  define viddiags_uv diags_eddyuv_avg
#  define viddiags_ub diags_eddyub_avg
#  define viddiags_vb diags_eddyvb_avg
#  define viddiags_uw diags_eddyuw_avg
#  define viddiags_vw diags_eddyvw_avg
#  define viddiags_wb diags_eddywb_avg
#  define viddiags_ubu diags_eddyubu_avg
#  define viddiags_vbv diags_eddyvbv_avg
#  define viddiags_usu diags_eddyusu_avg
#  define viddiags_vsv diags_eddyvsv_avg
#  define viddiags_ugsu diags_eddyugsu_avg
#  define viddiags_vgsv diags_eddyvgsv_avg
 
# endif
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
# ifdef DO_NOT_OVERWRITE
      if (lvar.eq.total_rec-1) then
#  ifndef AVRH
        create_new_file=ldefdiags_eddy
#  else
        create_new_file=ldefdiags_eddy_avg
#  endif
      else
        create_new_file=.false.
      endif
# else
#  ifndef AVRH
      create_new_file=ldefdiags_eddy
#  else
      create_new_file=ldefdiags_eddy_avg
#  endif
# endif
!
      if (ncid.ne.-1) create_new_file=.false.
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode.gt.0) create_new_file=.false.
# endif
!
 10   if (create_new_file) then
        lstr=lenstr(ncname)
# ifndef NC4PAR
        ierr=nf_create(ncname(1:lstr), NF_CLOBBER, ncid)
# else
        cmode = ior(nf_netcdf4, nf_classic_model)
        cmode = ior(cmode, nf_mpiio)
        csize=xi_rho*eta_rho/NNODES
        MPI_master_only write(stdout,*)
     &       'CREATE DIAGS_EDDY HIS/AVG NC4 PARALLEL FILE'
        ierr=nf_create_par(ncname(1:lstr), cmode,
     &                     MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
# endif
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
# ifdef SOLVE3D
        ierr=nf_def_dim (ncid, 's_rho',    N,        r3dgrd(3))
        ierr=nf_def_dim (ncid, 's_w',      N+1,      w3dgrd(3))
# endif
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

# ifdef SOLVE3D

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

# endif
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
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid, vidTstep, nf_collective)
# endif
        ierr=nf_put_att_text (ncid, vidTstep, 'long_name', 48,
     &       'time step and record numbers from initialization')
!
! Time.
!
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, vidTime)
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid, vidTime, nf_collective)
# endif
# ifndef AVRH
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', 
     &       lvar, vname(2,indxTime)(1:lvar))
# else
        text='avg'/ /vname(2,indxTime)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', 
     &       lvar, text(1:lvar))
# endif
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
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid, vidTime2, nf_collective)
# endif
# ifndef AVRH
        lvar=lenstr(vname(2,indxTime2))
        ierr=nf_put_att_text (ncid, vidTime2, 'long_name', 
     &       lvar, vname(2,indxTime)(1:lvar))
# else
        text='avg'/ /vname(2,indxTime2)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidTime2, 'long_name', 
     &       lvar, text(1:lvar))
# endif
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

! eddy diagnostics variables.
! -------- ----------- --------

        itrc=1

        if (wrt(itrc)) then

!------------------------------
! indxeddyzz
!          
           lvar=lenstr(vname(1,indxeddyzz+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyzz+itrc-1)(1:lvar),
     &                      NF_FOUT, 3, r2dgrd, viddiags_zz(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid,viddiags_zz(itrc),nf_collective)
# endif
# ifndef AVRH
           lvar=lenstr(vname(2,indxeddyzz+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_zz(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyzz+itrc-1)(1:lvar))
# else
           text='averaged '/ /vname(2,indxeddyzz+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, viddiags_zz(itrc), 'long_name',
     &                                            lvar, text(1:lvar))
# endif
           lvar=lenstr(vname(3,indxeddyzz+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_zz(itrc), 'units',
     &                          lvar,vname(3,indxeddyzz+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyzz+itrc-1))
           ierr=nf_put_att_text (ncid,viddiags_zz(itrc), 'field',
     &                       lvar, vname(4,indxeddyzz+itrc-1)(1:lvar))

           call nf_add_attribute(ncid,viddiags_zz(itrc)
     &                       ,indxeddyzz+itrc-1,5,NF_FOUT, ierr)

!------------------------------
! indxeddyuu
!          
           lvar=lenstr(vname(1,indxeddyuu+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyuu+itrc-1)(1:lvar),
     &                      NF_FOUT, 4, r3dgrd, viddiags_uu(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid,viddiags_uu(itrc),nf_collective)
# endif
# ifndef AVRH
           lvar=lenstr(vname(2,indxeddyuu+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_uu(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyuu+itrc-1)(1:lvar))
# else
           text='averaged '/ /vname(2,indxeddyuu+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, viddiags_uu(itrc), 'long_name',
     &                                            lvar, text(1:lvar))
# endif
           lvar=lenstr(vname(3,indxeddyuu+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_uu(itrc), 'units',
     &                          lvar,vname(3,indxeddyuu+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyuu+itrc-1))
           ierr=nf_put_att_text (ncid,viddiags_uu(itrc), 'field',
     &                       lvar, vname(4,indxeddyuu+itrc-1)(1:lvar))

           call nf_add_attribute(ncid,viddiags_uu(itrc)
     &                       ,indxeddyuu+itrc-1,5,NF_FOUT, ierr)
!
!------------------------------
! indxeddyvv
!          
           lvar=lenstr(vname(1,indxeddyvv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvv+itrc-1)(1:lvar),
     &                      NF_FOUT, 4, r3dgrd, viddiags_vv(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid,viddiags_vv(itrc),nf_collective)
# endif
# ifndef AVRH
           lvar=lenstr(vname(2,indxeddyvv+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_vv(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyvv+itrc-1)(1:lvar))
# else
           text='averaged '/ /vname(2,indxeddyvv+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, viddiags_vv(itrc), 'long_name',
     &                                            lvar, text(1:lvar))
# endif
           lvar=lenstr(vname(3,indxeddyvv+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_vv(itrc), 'units',
     &                          lvar,vname(3,indxeddyvv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvv+itrc-1))
           ierr=nf_put_att_text (ncid,viddiags_vv(itrc), 'field',
     &                       lvar, vname(4,indxeddyvv+itrc-1)(1:lvar))

           call nf_add_attribute(ncid,viddiags_vv(itrc)
     &                       ,indxeddyvv+itrc-1,5,NF_FOUT, ierr)
!
!------------------------------
! indxeddyuv
!          
           lvar=lenstr(vname(1,indxeddyuv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyuv+itrc-1)(1:lvar),
     &                      NF_FOUT, 4, r3dgrd, viddiags_uv(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid,viddiags_uv(itrc),nf_collective)
# endif
# ifndef AVRH
           lvar=lenstr(vname(2,indxeddyuv+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_uv(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyuv+itrc-1)(1:lvar))
# else
           text='averaged '/ /vname(2,indxeddyuv+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, viddiags_uv(itrc), 'long_name',
     &                                            lvar, text(1:lvar))
# endif
           lvar=lenstr(vname(3,indxeddyuv+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_uv(itrc), 'units',
     &                          lvar,vname(3,indxeddyuv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyuv+itrc-1))
           ierr=nf_put_att_text (ncid,viddiags_uv(itrc), 'field',
     &                       lvar, vname(4,indxeddyuv+itrc-1)(1:lvar))

           call nf_add_attribute(ncid,viddiags_uv(itrc)
     &                       ,indxeddyuv+itrc-1,5,NF_FOUT, ierr)
!
!------------------------------
! indxeddyuw
!          
           lvar=lenstr(vname(1,indxeddyuw+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyuw+itrc-1)(1:lvar),
     &                      NF_FOUT, 4, r3dgrd, viddiags_uw(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid,viddiags_uw(itrc),nf_collective)
# endif
# ifndef AVRH
           lvar=lenstr(vname(2,indxeddyuw+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_uw(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyuw+itrc-1)(1:lvar))
# else
           text='averaged '/ /vname(2,indxeddyuw+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, viddiags_uw(itrc), 'long_name',
     &                                            lvar, text(1:lvar))
# endif
           lvar=lenstr(vname(3,indxeddyuw+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_uw(itrc), 'units',
     &                          lvar,vname(3,indxeddyuw+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyuw+itrc-1))
           ierr=nf_put_att_text (ncid,viddiags_uw(itrc), 'field',
     &                       lvar, vname(4,indxeddyuw+itrc-1)(1:lvar))

           call nf_add_attribute(ncid,viddiags_uw(itrc)
     &                       ,indxeddyuw+itrc-1,5,NF_FOUT, ierr)
!
!------------------------------
! indxeddyvw
!          
           lvar=lenstr(vname(1,indxeddyvw+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvw+itrc-1)(1:lvar),
     &                      NF_FOUT, 4, r3dgrd, viddiags_vw(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid,viddiags_vw(itrc),nf_collective)
# endif
# ifndef AVRH
           lvar=lenstr(vname(2,indxeddyvw+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_vw(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyvw+itrc-1)(1:lvar))
# else
           text='averaged '/ /vname(2,indxeddyvw+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, viddiags_vw(itrc), 'long_name',
     &                                            lvar, text(1:lvar))
# endif
           lvar=lenstr(vname(3,indxeddyvw+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_vw(itrc), 'units',
     &                          lvar,vname(3,indxeddyvw+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvw+itrc-1))
           ierr=nf_put_att_text (ncid,viddiags_vw(itrc), 'field',
     &                       lvar, vname(4,indxeddyvw+itrc-1)(1:lvar))

           call nf_add_attribute(ncid,viddiags_vw(itrc)
     &                       ,indxeddyvw+itrc-1,5,NF_FOUT, ierr)
!
!------------------------------
! indxeddyub
!          
           lvar=lenstr(vname(1,indxeddyub+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyub+itrc-1)(1:lvar),
     &                      NF_FOUT, 4, r3dgrd, viddiags_ub(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid,viddiags_ub(itrc),nf_collective)
# endif
# ifndef AVRH
           lvar=lenstr(vname(2,indxeddyub+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_ub(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyub+itrc-1)(1:lvar))
# else
           text='averaged '/ /vname(2,indxeddyub+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, viddiags_ub(itrc), 'long_name',
     &                                            lvar, text(1:lvar))
# endif
           lvar=lenstr(vname(3,indxeddyub+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_ub(itrc), 'units',
     &                          lvar,vname(3,indxeddyub+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyub+itrc-1))
           ierr=nf_put_att_text (ncid,viddiags_ub(itrc), 'field',
     &                       lvar, vname(4,indxeddyub+itrc-1)(1:lvar))

           call nf_add_attribute(ncid,viddiags_ub(itrc)
     &                       ,indxeddyub+itrc-1,5,NF_FOUT, ierr)
!
!------------------------------
! indxeddyvb
!          
           lvar=lenstr(vname(1,indxeddyvb+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvb+itrc-1)(1:lvar),
     &                      NF_FOUT, 4, r3dgrd, viddiags_vb(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid,viddiags_vb(itrc),nf_collective)
# endif
# ifndef AVRH
           lvar=lenstr(vname(2,indxeddyvb+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_vb(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyvb+itrc-1)(1:lvar))
# else
           text='averaged '/ /vname(2,indxeddyvb+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, viddiags_vb(itrc), 'long_name',
     &                                            lvar, text(1:lvar))
# endif
           lvar=lenstr(vname(3,indxeddyvb+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_vb(itrc), 'units',
     &                          lvar,vname(3,indxeddyvb+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvb+itrc-1))
           ierr=nf_put_att_text (ncid,viddiags_vb(itrc), 'field',
     &                       lvar, vname(4,indxeddyvb+itrc-1)(1:lvar))

           call nf_add_attribute(ncid,viddiags_vb(itrc)
     &                       ,indxeddyvb+itrc-1,5,NF_FOUT, ierr)
!
!------------------------------
! indxeddywb
!          
           lvar=lenstr(vname(1,indxeddywb+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddywb+itrc-1)(1:lvar),
     &                      NF_FOUT, 4, r3dgrd, viddiags_wb(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid,viddiags_wb(itrc),nf_collective)
# endif
# ifndef AVRH
           lvar=lenstr(vname(2,indxeddywb+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_wb(itrc), 'long_name',
     &                       lvar, vname(2,indxeddywb+itrc-1)(1:lvar))
# else
           text='averaged '/ /vname(2,indxeddywb+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, viddiags_wb(itrc), 'long_name',
     &                                            lvar, text(1:lvar))
# endif
           lvar=lenstr(vname(3,indxeddywb+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_wb(itrc), 'units',
     &                          lvar,vname(3,indxeddywb+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddywb+itrc-1))
           ierr=nf_put_att_text (ncid,viddiags_wb(itrc), 'field',
     &                       lvar, vname(4,indxeddywb+itrc-1)(1:lvar))

           call nf_add_attribute(ncid,viddiags_wb(itrc)
     &                       ,indxeddywb+itrc-1,5,NF_FOUT, ierr)
!
!------------------------------
! indxeddyubu
!
           lvar=lenstr(vname(1,indxeddyubu+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyubu+itrc-1)(1:lvar),
     &                      NF_FOUT, 3, u2dgrd, viddiags_ubu(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid,viddiags_ubu(itrc),nf_collective)
# endif
# ifndef AVRH
           lvar=lenstr(vname(2,indxeddyubu+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_ubu(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyubu+itrc-1)(1:lvar))
# else
           text='averaged '/ /vname(2,indxeddyubu+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, viddiags_ubu(itrc), 'long_name',
     &                                            lvar, text(1:lvar))
# endif
           lvar=lenstr(vname(3,indxeddyubu+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_ubu(itrc), 'units',
     &                         lvar,vname(3,indxeddyubu+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyubu+itrc-1))
           ierr=nf_put_att_text (ncid,viddiags_ubu(itrc), 'field',
     &                       lvar, vname(4,indxeddyubu+itrc-1)(1:lvar))

           call nf_add_attribute(ncid,viddiags_ubu(itrc)
     &                       ,indxeddyubu+itrc-1,5,NF_FOUT, ierr)
!
!------------------------------
! indxeddyvbv
!
           lvar=lenstr(vname(1,indxeddyvbv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvbv+itrc-1)(1:lvar),
     &                      NF_FOUT, 3, v2dgrd, viddiags_vbv(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid,viddiags_vbv(itrc),nf_collective)
# endif
# ifndef AVRH
           lvar=lenstr(vname(2,indxeddyvbv+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_vbv(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyvbv+itrc-1)(1:lvar))
# else
           text='averaged '/ /vname(2,indxeddyvbv+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, viddiags_vbv(itrc), 'long_name',
     &                                            lvar, text(1:lvar))
# endif
           lvar=lenstr(vname(3,indxeddyvbv+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_vbv(itrc), 'units',
     &                         lvar,vname(3,indxeddyvbv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvbv+itrc-1))
           ierr=nf_put_att_text (ncid,viddiags_vbv(itrc), 'field',
     &                       lvar, vname(4,indxeddyvbv+itrc-1)(1:lvar))

           call nf_add_attribute(ncid,viddiags_vbv(itrc)
     &                       ,indxeddyvbv+itrc-1,5,NF_FOUT, ierr)    
!
!------------------------------
! indxeddyusu
!
           lvar=lenstr(vname(1,indxeddyusu+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyusu+itrc-1)(1:lvar),
     &                      NF_FOUT, 3, u2dgrd, viddiags_usu(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid,viddiags_usu(itrc),nf_collective)
# endif
# ifndef AVRH
           lvar=lenstr(vname(2,indxeddyusu+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_usu(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyusu+itrc-1)(1:lvar))
# else
           text='averaged '/ /vname(2,indxeddyusu+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, viddiags_usu(itrc), 'long_name',
     &                                            lvar, text(1:lvar))
# endif
           lvar=lenstr(vname(3,indxeddyusu+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_usu(itrc), 'units',
     &                         lvar,vname(3,indxeddyusu+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyusu+itrc-1))
           ierr=nf_put_att_text (ncid,viddiags_usu(itrc), 'field',
     &                       lvar, vname(4,indxeddyusu+itrc-1)(1:lvar))

           call nf_add_attribute(ncid,viddiags_usu(itrc)
     &                       ,indxeddyusu+itrc-1,5,NF_FOUT, ierr)
!
!------------------------------
! indxeddyvsv
!
           lvar=lenstr(vname(1,indxeddyvsv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvsv+itrc-1)(1:lvar),
     &                      NF_FOUT, 3, v2dgrd, viddiags_vsv(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid,viddiags_vsv(itrc),nf_collective)
# endif
# ifndef AVRH
           lvar=lenstr(vname(2,indxeddyvsv+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_vsv(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyvsv+itrc-1)(1:lvar))
# else
           text='averaged '/ /vname(2,indxeddyvsv+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, viddiags_vsv(itrc), 'long_name',
     &                                            lvar, text(1:lvar))
# endif
           lvar=lenstr(vname(3,indxeddyvsv+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_vsv(itrc), 'units',
     &                         lvar,vname(3,indxeddyvsv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvsv+itrc-1))
           ierr=nf_put_att_text (ncid,viddiags_vsv(itrc), 'field',
     &                       lvar, vname(4,indxeddyvsv+itrc-1)(1:lvar))

           call nf_add_attribute(ncid,viddiags_vsv(itrc)
     &                       ,indxeddyvsv+itrc-1,5,NF_FOUT, ierr)
     

!------------------------------
!
! indxeddyugsu
!
           lvar=lenstr(vname(1,indxeddyugsu+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyugsu+itrc-1)(1:lvar),
     &                      NF_FOUT, 3, u2dgrd, viddiags_ugsu(itrc))
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,viddiags_ugsu(itrc),nf_collective)
# endif
# ifndef AVRH
           lvar=lenstr(vname(2,indxeddyugsu+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_ugsu(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyugsu+itrc-1)(1:lvar))
# else
           text='averaged '/ /vname(2,indxeddyugsu+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, viddiags_ugsu(itrc), 'long_name',
     &                                            lvar, text(1:lvar))
# endif
           lvar=lenstr(vname(3,indxeddyugsu+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_ugsu(itrc), 'units',
     &            lvar,vname(3,indxeddyugsu+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyugsu+itrc-1))
           ierr=nf_put_att_text (ncid,viddiags_ugsu(itrc), 'field',
     &                       lvar, vname(4,indxeddyugsu+itrc-1)(1:lvar))

       call nf_add_attribute(ncid,viddiags_ugsu(itrc)
     &                       ,indxeddyugsu+itrc-1,5,NF_FOUT, ierr)

!
! indxeddyvgsv
!

           lvar=lenstr(vname(1,indxeddyvgsv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvgsv+itrc-1)(1:lvar),
     &                      NF_FOUT, 3, v2dgrd, viddiags_vgsv(itrc))
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,viddiags_vgsv(itrc),nf_collective)
# endif
# ifndef AVRH
           lvar=lenstr(vname(2,indxeddyvgsv+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_vgsv(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyvgsv+itrc-1)(1:lvar))
# else
           text='averaged '/ /vname(2,indxeddyvgsv+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, viddiags_vgsv(itrc), 'long_name',
     &                                            lvar, text(1:lvar))
# endif
           lvar=lenstr(vname(3,indxeddyvgsv+itrc-1))
           ierr=nf_put_att_text (ncid, viddiags_vgsv(itrc), 'units',
     &               lvar,vname(3,indxeddyvgsv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvgsv+itrc-1))
           ierr=nf_put_att_text (ncid,viddiags_vgsv(itrc), 'field',
     &                       lvar, vname(4,indxeddyvgsv+itrc-1)(1:lvar))

       call nf_add_attribute(ncid,viddiags_vgsv(itrc)
     &                       ,indxeddyvgsv+itrc-1,5,NF_FOUT, ierr)


!
        endif
!
! Leave definition mode.
! ----- ---------- -----
!
        ierr=nf_enddef(ncid)
        MPI_master_only write(stdout,'(6x,4A,1x,A,i4)')
# ifdef AVRH
     &        'DEF_DIAGS_EDDY_AVG - Created ',
# else      
     &        'DEF_DIAGS_EDDY - Created ',
# endif /*AVRH*/
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
# ifndef NC4PAR
        ierr=nf_open (ncname(1:lstr), nf_write, ncid)
# else
        ierr=nf_open_par (ncname(1:lstr), IOR(nf_write, nf_mpiio),
     &                    MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
# endif
        if (ierr .eq. nf_noerr) then
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
     &           ) 'WARNING: def_diags_eddy: Actual number of records',
     &               rec,  'in netCDF file',  '''',  ncname(1:lstr),
     &             ''' exceeds the record number from restart data',
     &             rec+1-ierr,'/', total_rec,', restart is assumed.'
              rec=rec-ierr
            elseif (rec_per_file.eq.0) then
              total_rec=rec+1           ! <-- set to the next record
# if defined MPI & !defined PARALLEL_FILES
              if (mynode.gt.0) total_rec=total_rec-1
# endif
            endif
            ierr=nf_noerr
          endif
        endif
!
        if (ierr. ne. nf_noerr) then
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
          if (mynode.eq.0) then
            create_new_file=.true.
            goto 10
          else
            write(stdout,'(/1x,4A,2x,A,I4/)') 'def_his/avg ERROR: ',
     &         'Cannot open file ''', ncname(1:lstr), '''.'
     &                   MYID
            goto 99                                       !--> ERROR
          endif
# else
          create_new_file=.true.
          goto 10
# endif
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
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid, vidTstep, nf_collective)
# endif
 
!
! Time.
!
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),vidTime)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime)(1:lvar), ncname(1:lstr)
          goto 99                                         !--> ERROR 
        endif
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid, vidTime, nf_collective)
# endif

!
! Time2.
!
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_inq_varid (ncid,vname(1,indxTime2)(1:lvar),vidTime2)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime2)(1:lvar), ncname(1:lstr)
          goto 99                                         !--> ERROR 
        endif
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid, vidTime2, nf_collective)
# endif


! momentum flux diagnostics variables :
! -------- ---- ----------- ----------

        itrc=1

        if (wrt(itrc)) then

! ---------------------------------
! indxeddyzz
!
          lvar=lenstr(vname(1,indxeddyzz+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyzz+itrc-1)(1:lvar),
     &                                            viddiags_zz(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyzz+itrc-1)(1:lvar),
     &                                         ncname(1:lstr) 
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_zz(itrc), nf_collective)
# endif
!
! ---------------------------------
! indxeddyuu
!
          lvar=lenstr(vname(1,indxeddyuu+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyuu+itrc-1)(1:lvar),
     &                                            viddiags_uu(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyuu+itrc-1)(1:lvar),
     &                                         ncname(1:lstr) 
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_uu(itrc), nf_collective)
# endif
!
! ---------------------------------
! indxeddyvv
!
          lvar=lenstr(vname(1,indxeddyvv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvv+itrc-1)(1:lvar),
     &                                            viddiags_vv(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvv+itrc-1)(1:lvar),
     &                                         ncname(1:lstr) 
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_vv(itrc), nf_collective)
# endif
!
! ---------------------------------
! indxeddyuv
!
          lvar=lenstr(vname(1,indxeddyuv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyuv+itrc-1)(1:lvar),
     &                                            viddiags_uv(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyuv+itrc-1)(1:lvar),
     &                                         ncname(1:lstr) 
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_uv(itrc), nf_collective)
# endif
!
! ---------------------------------
! indxeddyuw
!
          lvar=lenstr(vname(1,indxeddyuw+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyuw+itrc-1)(1:lvar),
     &                                            viddiags_uw(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyuw+itrc-1)(1:lvar),
     &                                         ncname(1:lstr) 
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_uw(itrc), nf_collective)
# endif
!
! ---------------------------------
! indxeddyvw
!
          lvar=lenstr(vname(1,indxeddyvw+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvw+itrc-1)(1:lvar),
     &                                            viddiags_vw(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvw+itrc-1)(1:lvar),
     &                                         ncname(1:lstr) 
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_vw(itrc), nf_collective)
# endif
!
! ---------------------------------
! indxeddyub
!
          lvar=lenstr(vname(1,indxeddyub+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyub+itrc-1)(1:lvar),
     &                                            viddiags_ub(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyub+itrc-1)(1:lvar),
     &                                         ncname(1:lstr) 
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_ub(itrc), nf_collective)
# endif
!
! ---------------------------------
! indxeddyvb
!
          lvar=lenstr(vname(1,indxeddyvb+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvb+itrc-1)(1:lvar),
     &                                            viddiags_vb(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvb+itrc-1)(1:lvar),
     &                                         ncname(1:lstr) 
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_vb(itrc), nf_collective)
# endif
!
! ---------------------------------
! indxeddywb
!
          lvar=lenstr(vname(1,indxeddywb+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddywb+itrc-1)(1:lvar),
     &                                            viddiags_wb(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddywb+itrc-1)(1:lvar),
     &                                         ncname(1:lstr) 
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_wb(itrc), nf_collective)
# endif
!
! ---------------------------------
! indxeddyubu
!
          lvar=lenstr(vname(1,indxeddyubu+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyubu+itrc-1)(1:lvar),
     &                                            viddiags_ubu(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyubu+itrc-1)(1:lvar),
     &                                         ncname(1:lstr)
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_ubu(itrc),nf_collective)
# endif
!
! ---------------------------------
! indxeddyvbv
!
          lvar=lenstr(vname(1,indxeddyvbv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvbv+itrc-1)(1:lvar),
     &                                            viddiags_vbv(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvbv+itrc-1)(1:lvar),
     &                                         ncname(1:lstr)
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_vbv(itrc),nf_collective)
# endif
!
! ---------------------------------
! indxeddyusu
!
          lvar=lenstr(vname(1,indxeddyusu+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyusu+itrc-1)(1:lvar),
     &                                            viddiags_usu(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyusu+itrc-1)(1:lvar),
     &                                         ncname(1:lstr)
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_usu(itrc),nf_collective)
# endif
!
! ---------------------------------
! indxeddyvsv
!
          lvar=lenstr(vname(1,indxeddyvsv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvsv+itrc-1)(1:lvar),
     &                                            viddiags_vsv(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvsv+itrc-1)(1:lvar),
     &                                         ncname(1:lstr)
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_vsv(itrc),nf_collective)
# endif
!
! ---------------------------------
! indxeddyugsu
!
          lvar=lenstr(vname(1,indxeddyugsu+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyugsu+itrc-1)(1:lvar),
     &                                            viddiags_ugsu(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyugsu+itrc-1)(1:lvar),
     &                                         ncname(1:lstr)
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_ugsu(itrc),nf_collective)
# endif
!
! ---------------------------------
! indxeddyvgsv
!
          lvar=lenstr(vname(1,indxeddyvgsv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvgsv+itrc-1)(1:lvar),
     &                                            viddiags_vgsv(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvgsv+itrc-1)(1:lvar),
     &                                         ncname(1:lstr)
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_vgsv(itrc),nf_collective)
# endif
!
        endif ! wrt(itrc)
!
! ---------------------------------
!
        MPI_master_only write(*,'(6x,2A,i4,1x,A,i4)') 
     &                     'def_diags_eddy: -- Opened ',
     &                     'existing file  from record =', rec
     &                      MYID
!
        MPI_master_only write(*,'(6x,2A,i4,1x,A,i4)') 
     &                     'def_diags_eddy: -- Opened ',
     &                     'existing file  from record =', rec
     &                      MYID
!
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      else
        ierr=nf_open (ncname(1:lstr), nf_write, ncid)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,'(/1x,4A,2x,A,I4/)')
     &                'def_diags_eddy: ERROR: ',
     &                'Cannot open file ''', ncname(1:lstr), '''.'
     &                 MYID
          goto 99                                         !--> ERROR
        endif
# endif
!
      endif             !<--  create_new_file  
!
      ierr=nf_set_fill (ncid, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
        MPI_master_only write(*,'(6x,2A,i4,1x,A,i4)')
     &    'def_diags_eddy ERROR: Cannot ',
     &    'switch to ''nf_nofill'' more; netCDF error code =', ierr 
      endif
   1  format(/1x,'def_diags_eddy ERROR: Cannot find variable ''',
     &                   A, ''' in netCDF file ''', A, '''.'/)

!
  11  format(/' def_diags_eddy - unable to create diag file: ',a)
  20  format(/' def_diags_eddy - error while writing variable: ',a,
     &        /,15x,'into diag  file: ',a)
!
  99  return
      end
!
# undef ncname
# undef rec_per_file
# undef wrt
# undef vidTime
# undef vidTime2
# undef vidTstep
!
# undef viddiags_zz
# undef viddiags_uu
# undef viddiags_vv
# undef viddiags_uv
# undef viddiags_ub
# undef viddiags_vb
# undef viddiags_wb
# undef viddiags_uw
# undef viddiags_vw
# undef viddiags_ubu
# undef viddiags_vbv
# undef viddiags_usu
# undef viddiags_vsv
# undef viddiags_ugsu
# undef viddiags_vgsv
!
# ifdef AVERAGES
#  ifndef AVRH
#   define AVRH
#   include "def_diags_eddy.F"
#  endif
# endif /* AVERAGES */
!
#else /*DIAGNOSTICS_EDDY */
      subroutine def_diags_eddy_empty()
      return
      end
#endif /* DIAGNOSTICS_EDDY */

