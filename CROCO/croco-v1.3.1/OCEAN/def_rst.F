! $Id: def_rst.F 1571 2014-07-01 12:38:05Z gcambon $
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
                                                  ! Create/open
      subroutine def_rst (ncid, total_rec, ierr)  ! restart netCDF
                                                  ! file.

#if defined FLOATS && defined AGRIF
        USE Agrif_Util
#endif
      implicit none  
      logical create_new_file
      integer ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim
     &      , r2dgrd(3),  u2dgrd(3), v2dgrd(3),  auxil(2),  checkdims
     &      , b3dgrd(4)
#ifdef NC4PAR
     &      , csize,cmode
#endif
#if defined SEDIMENT || defined BBL
     &      ,  indxWrk
#endif
#ifdef MUSTANG
     &      ,  indx 
#endif
#ifdef SOLVE3D
     &      , r3dgrd(4),  u3dgrd(4), v3dgrd(4),  w3dgrd(4), itrc
#endif
#ifdef USE_CALENDAR
      CHARACTER (len=19)    :: cdate,tool_sectodat
#endif
#ifdef FLOATS
      integer fltdim, NFVdim, Tinfovardim, NFTdim
     & , dim_Tinfo(3), dim_track(4), dim_fltgrd(2)
      character*65 vinfo(4)
#endif

#include "param.h"
#include "scalars.h"
#include "ncscrum.h"
#include "netcdf.inc"
#ifdef FLOATS
# include "floats.h"
# include "ncscrum_floats.h"
#endif
#ifdef NC4PAR
# include "mpi_cpl.h"
      include 'mpif.h'
#endif
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
      ierr=0
      lstr=lenstr(rstname)
      if (nrpfrst.gt.0) then
        lvar=total_rec - (1+mod(total_rec-1, nrpfrst))
        call insert_time_index (rstname, lstr, lvar, ierr)
#ifdef USE_CALENDAR
        if (nrpfrst.eq.1) then
          cdate = tool_sectodat(time)
          rstname=TRIM(rstname(1:lstr-9))//'.'//cdate(7:10)//cdate(4:5)
          rstname=TRIM(rstname)//cdate(1:2)//cdate(12:13)//cdate(15:16)
          rstname=TRIM(rstname)//cdate(18:19)//'.nc'
          lstr=lenstr(rstname)
        end if 
#endif
        if (ierr .ne. 0) goto 99
      endif
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
      create_new_file=ldefhis
      if (ncid.ne.-1) create_new_file=.false.
#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      if (mynode.gt.0) create_new_file=.false.
#endif
!
! Create new restart file:    Put global attributes
!======= === ======= =====    and define all variables.
!
  10  if (create_new_file) then

#ifndef NC4PAR
        ierr=nf_create(rstname(1:lstr),NF_CLOBBER, ncid)
#else
        cmode = ior(nf_netcdf4,nf_classic_model)
        cmode = ior(cmode, nf_mpiio)
        csize=xi_rho*eta_rho/NNODES
        MPI_master_only
     &       write(stdout,*)'CREATE RST NC4 PARALLEL FILE'
        ierr=nf_create_par(rstname(1:lstr),cmode,
     &        MPI_COMM_WORLD,MPI_INFO_NULL,ncid)
#endif

        if (ierr.ne.nf_noerr) then
          write(stdout,'(/3(1x,A)/)') 'ERROR in DEF_RST: Cannot',
     &             'create restart NetCDF file:', rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
        if (nrpfrst.eq.0) total_rec=0           
!
! Put global attributes.
! --- ------ -----------
!
        call put_global_atts (ncid, ierr)
!
! Define dimensions of staggered fields.
! ------ ---------- -- --------- -------
!
        ierr=nf_def_dim (ncid, 'xi_rho',   xi_rho,  r2dgrd(1))
        ierr=nf_def_dim (ncid, 'xi_u',     xi_u,     u2dgrd(1))
        ierr=nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr=nf_def_dim (ncid, 'eta_v',    eta_v,    v2dgrd(2))
#ifdef SEDIMENT
        ierr=nf_def_dim (ncid, 's_b',      NLAY,     b3dgrd(3))
#endif
#ifdef MUSTANG
        ierr=nf_def_dim (ncid, 's_b',      ksdmax,     b3dgrd(3))
#endif
#ifdef SOLVE3D
        ierr=nf_def_dim (ncid, 's_rho',    N,        r3dgrd(3))
        ierr=nf_def_dim (ncid, 's_w',      N+1,      w3dgrd(3))
#endif
        ierr=nf_def_dim (ncid, 'time', nf_unlimited, timedim)
#if defined OUTPUTS_SURFACE && ! defined XIOS
        ierr=nf_def_dim (ncid, 'auxil',    6,        auxil(1))
#else
        ierr=nf_def_dim (ncid, 'auxil',    4,        auxil(1))
#endif
        auxil(2)=timedim

        r2dgrd(3)=timedim           ! Free surface

        u2dgrd(2)=r2dgrd(2)         ! 2D UBAR-type
        u2dgrd(3)=timedim

        v2dgrd(1)=r2dgrd(1)         ! 2D VBAR-type
        v2dgrd(3)=timedim
        
        b3dgrd(1)=r2dgrd(1)         !
        b3dgrd(2)=r2dgrd(2)         ! 3D BED-type
        b3dgrd(4)=timedim           !
        
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
#endif
#ifdef FLOATS
# ifdef AGRIF
       if (Agrif_Root()) then
# endif
          ierr=nf_def_dim(ncid,'drifter',nfloats,fltdim)
          ierr=nf_def_dim(ncid,'Tinfovar',5,Tinfovardim)
          ierr=nf_def_dim(ncid,'NFV',6,NFVdim)
          ierr=nf_def_dim(ncid,'NFT',NFT+1,NFTdim)

          dim_Tinfo(1)=Tinfovardim
          dim_Tinfo(2)=fltdim
          dim_Tinfo(3)=timedim

          dim_track(1)=NFVdim
          dim_track(2)=NFTdim
          dim_track(3)=fltdim
          dim_track(4)=timedim

          dim_fltgrd(1)=fltdim
          dim_fltgrd(2)=timedim

# ifdef AGRIF
       endif
# endif
#endif /* FLOATS */

#ifdef PUT_GRID_INTO_RESTART
!
! Define grid variables.
! ------ ---- ----------
!        
        if (total_rec.le.1) then
# ifdef SOLVE3D
         call def_grid_3d(ncid, r2dgrd, u2dgrd, v2dgrd 
     &                   ,r3dgrd, w3dgrd)
# else
         call def_grid_2d(ncid, r2dgrd, u2dgrd, v2dgrd)
# endif
      endif
#endif
!
! Define evolving model variables:
! ------ -------- ----- ----------
!
!
! Time step number and time record numbers:
!
        ierr=nf_def_var (ncid, 'time_step', nf_int, 2, auxil,
     &       rstTstep)
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstTstep,nf_collective)
#endif
        ierr=nf_put_att_text (ncid, rstTstep, 'long_name', 48,
     &       'time step and record numbers from initialization')
!
! Time.
!
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar),
     &                              NF_DOUBLE, 1, timedim, rstTime)
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstTime,nf_collective)
#endif
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, rstTime, 'long_name', lvar,
     &                                  vname(2,indxTime)(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, rstTime, 'units',     lvar,
     &                                  vname(3,indxTime)(1:lvar))
        lvar=lenstr (vname(4,indxTime))
        ierr=nf_put_att_text(ncid, rstTime, 'field',     lvar,
     &                                  vname(4,indxTime)(1:lvar))

      call nf_add_attribute(ncid, rstTime, indxTime, 5, NF_DOUBLE, ierr)
!
! Time2.
!
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_def_var (ncid, vname(1,indxTime2)(1:lvar),
     &                              NF_DOUBLE, 1, timedim, rstTime2)
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstTime2,nf_collective)
#endif
        lvar=lenstr(vname(2,indxTime2))
        ierr=nf_put_att_text (ncid, rstTime2, 'long_name', lvar,
     &                                  vname(2,indxTime2)(1:lvar))
        lvar=lenstr(vname(3,indxTime2))
        ierr=nf_put_att_text (ncid, rstTime2, 'units',     lvar,
     &                                  vname(3,indxTime2)(1:lvar))
        lvar=lenstr (vname(4,indxTime2))
        ierr=nf_put_att_text(ncid, rstTime2, 'field',     lvar,
     &                                  vname(4,indxTime2)(1:lvar))

        call nf_add_attribute(ncid, rstTime2, indxTime2, 5,  
     &       NF_DOUBLE, ierr)
!     
!     Free-surface.
!
        lvar=lenstr(vname(1,indxZ))
        ierr=nf_def_var (ncid, vname(1,indxZ)(1:lvar),
     &                                  NF_DOUBLE, 3, r2dgrd,rstZ)
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstZ,nf_collective)
#endif
        lvar=lenstr(vname(2,indxZ))
        ierr=nf_put_att_text (ncid, rstZ, 'long_name', lvar,
     &                                     vname(2,indxZ)(1:lvar))
        lvar=lenstr(vname(3,indxZ))
        ierr=nf_put_att_text (ncid, rstZ, 'units', lvar,
     &                                     vname(3,indxZ)(1:lvar))
        lvar=lenstr(vname(4,indxZ))
        ierr=nf_put_att_text (ncid, rstZ, 'field', lvar,
     &                                     vname(4,indxZ)(1:lvar))

      call nf_add_attribute(ncid, rstZ, indxZ, 5, NF_DOUBLE, ierr)

!     
!   2D momenta in XI- and ETA-directions.
!
        lvar=lenstr(vname(1,indxUb))
        ierr=nf_def_var (ncid, vname(1,indxUb)(1:lvar),
     &                               NF_DOUBLE, 3, u2dgrd, rstUb)
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstUb,nf_collective)
#endif
        lvar=lenstr(vname(2,indxUb))
        ierr=nf_put_att_text (ncid, rstUb, 'long_name', lvar,
     &                                   vname(2,indxUb)(1:lvar))
        lvar=lenstr(vname(3,indxUb))
        ierr=nf_put_att_text (ncid, rstUb, 'units', lvar,
     &                                   vname(3,indxUb)(1:lvar))
        lvar=lenstr(vname(4,indxUb))
        ierr=nf_put_att_text (ncid, rstUb, 'field', lvar,
     &                                   vname(4,indxUb)(1:lvar))  
      call nf_add_attribute(ncid, rstUb, indxUb, 5, NF_DOUBLE, ierr)

        lvar=lenstr(vname(1,indxVb))
        ierr=nf_def_var (ncid, vname(1,indxVb)(1:lvar),
     &                               NF_DOUBLE, 3, v2dgrd, rstVb)
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstVb,nf_collective)
#endif
        lvar=lenstr(vname(2,indxVb))
        ierr=nf_put_att_text (ncid, rstVb, 'long_name', lvar,
     &                                   vname(2,indxVb)(1:lvar))
        lvar=lenstr(vname(3,indxVb))
        ierr=nf_put_att_text (ncid, rstVb, 'units',     lvar,
     &                                   vname(3,indxVb)(1:lvar))
        lvar=lenstr(vname(4,indxVb))
        ierr=nf_put_att_text (ncid, rstVb, 'field',     lvar,
     &                                  vname(4,indxVb)(1:lvar))
        
       call nf_add_attribute(ncid, rstVb, indxVb, 5, NF_DOUBLE, ierr)


#ifdef SOLVE3D
!
! 3D momenta in XI- and ETA-directions.
!
        lvar=lenstr(vname(1,indxU))
        ierr=nf_def_var (ncid, vname(1,indxU)(1:lvar),
     &                              NF_DOUBLE, 4, u3dgrd, rstU)
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstU,nf_collective)
# endif
        lvar=lenstr(vname(2,indxU))
        ierr=nf_put_att_text (ncid, rstU, 'long_name', lvar,
     &                                   vname(2,indxU)(1:lvar))
        lvar=lenstr(vname(3,indxU))
        ierr=nf_put_att_text (ncid, rstU, 'units',     lvar,
     &                                   vname(3,indxU)(1:lvar))
        lvar=lenstr(vname(4,indxU))
        ierr=nf_put_att_text (ncid, rstU, 'field',     lvar,
     &                                  vname(4,indxU)(1:lvar))

        call nf_add_attribute(ncid, rstU, indxU, 5, NF_DOUBLE, ierr)

        lvar=lenstr(vname(1,indxV))
        ierr=nf_def_var (ncid, vname(1,indxV)(1:lvar),
     &                              NF_DOUBLE, 4, v3dgrd, rstV)
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstV,nf_collective)
# endif
        lvar=lenstr(vname(2,indxV))
        ierr=nf_put_att_text (ncid, rstV, 'long_name', lvar,
     &                                   vname(2,indxV)(1:lvar))
        lvar=lenstr(vname(3,indxV))
        ierr=nf_put_att_text (ncid, rstV, 'units',     lvar,
     &                                   vname(3,indxV)(1:lvar))
        lvar=lenstr(vname(4,indxV))
        ierr=nf_put_att_text (ncid, rstV, 'field',     lvar,
     &                                  vname(4,indxV)(1:lvar))

        call nf_add_attribute(ncid, rstV, indxV, 5, NF_DOUBLE, ierr)

!
! Tracer variables.
!
# if defined TRACERS
        do itrc=1,NT
          lvar=lenstr(vname(1,indxV+itrc))
          ierr=nf_def_var (ncid, vname(1,indxV+itrc)(1:lvar),
     &         NF_DOUBLE, 4, r3dgrd, rstT(itrc))
#  ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstT(itrc),nf_collective)
#  endif
          lvar=lenstr(vname(2,indxV+itrc))
          ierr=nf_put_att_text (ncid, rstT(itrc), 'long_name',
     &                     lvar,   vname(2,indxV+itrc)(1:lvar))
          lvar=lenstr(vname(3,indxV+itrc))
          ierr=nf_put_att_text (ncid, rstT(itrc), 'units', lvar,
     &         vname(3,indxV+itrc)(1:lvar))
          lvar=lenstr(vname(4,indxV+itrc))
          ierr=nf_put_att_text (ncid, rstT(itrc), 'field', lvar,
     &         vname(4,indxV+itrc)(1:lvar))

          call nf_add_attribute(ncid, rstT(itrc), indxV+itrc, 5, 
     &         NF_DOUBLE,ierr)          
        enddo
# endif
        
# if defined LMD_MIXING 
!
! Depth of planetary boundary layer.
!     
#  if defined LMD_SKPP
        lvar=lenstr(vname(1,indxHbl))
        ierr=nf_def_var (ncid, vname(1,indxHbl)(1:lvar),
     &       NF_DOUBLE, 3, r2dgrd, rstHbl)
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstHbl,nf_collective)
#   endif
        lvar=lenstr(vname(2,indxHbl))
        ierr=nf_put_att_text (ncid, rstHbl, 'long_name',
     &       lvar,   vname(2,indxHbl)(1:lvar))
        lvar=lenstr(vname(3,indxHbl))
        ierr=nf_put_att_text (ncid, rstHbl, 'units', lvar,
     &       vname(3,indxHbl)(1:lvar))
        lvar=lenstr(vname(4,indxHbl))
        ierr=nf_put_att_text (ncid, rstHbl, 'field', lvar,
     &       vname(4,indxHbl)(1:lvar))
        
        call nf_add_attribute(ncid, rstHbl, indxHbl, 5,
     &       NF_DOUBLE,ierr)
#  endif        

#  ifdef LMD_BKPP
!
! Depth of bottom planetary boundary layer.
!        
        lvar=lenstr(vname(1,indxHbbl))
        ierr=nf_def_var (ncid, vname(1,indxHbbl)(1:lvar),
     &       NF_DOUBLE, 3, r2dgrd, rstHbbl)
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstHbbl,nf_collective)
#   endif
        lvar=lenstr(vname(2,indxHbbl))
        ierr=nf_put_att_text (ncid, rstHbbl, 'long_name',
     &       lvar,   vname(2,indxHbbl)(1:lvar))
        lvar=lenstr(vname(3,indxHbbl))
        ierr=nf_put_att_text (ncid, rstHbbl, 'units', lvar,
     &       vname(3,indxHbbl)(1:lvar))
        lvar=lenstr(vname(4,indxHbbl))
        ierr=nf_put_att_text (ncid, rstHbbl, 'field', lvar,
     &       vname(4,indxHbbl)(1:lvar))
        
        call nf_add_attribute(ncid, rstHbbl, indxHbbl, 5,
     &       NF_DOUBLE,ierr)
        
#  endif /*  LMD_BKPP  */ 
# endif /* LMD_MIXING */

# if defined GLS_MIXING
!
! Turbulent kinetic energy.
!
          lvar=lenstr(vname(1,indxTke))
          ierr=nf_def_var (ncid, vname(1,indxTke)(1:lvar),
     &         NF_DOUBLE, 4, w3dgrd, rstTke)
#  ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstTke,nf_collective)
#  endif
          lvar=lenstr(vname(2,indxTke))
          ierr=nf_put_att_text (ncid, rstTke, 'long_name',
     &                     lvar,   vname(2,indxTke)(1:lvar))
          lvar=lenstr(vname(3,indxTke))
          ierr=nf_put_att_text (ncid, rstTke, 'units', lvar,
     &         vname(3,indxTke)(1:lvar))
          lvar=lenstr(vname(4,indxTke))
          ierr=nf_put_att_text (ncid, rstTke, 'field', lvar,
     &         vname(4,indxTke)(1:lvar))

          call nf_add_attribute(ncid, rstTke, indxTke, 5,
     &         NF_DOUBLE,ierr)

!
! Generic length scale.
!
          lvar=lenstr(vname(1,indxGls))
          ierr=nf_def_var (ncid, vname(1,indxGls)(1:lvar),
     &         NF_DOUBLE, 4, w3dgrd, rstGls)
#  ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstGls,nf_collective)
#  endif
          lvar=lenstr(vname(2,indxGls))
          ierr=nf_put_att_text (ncid, rstGls, 'long_name',
     &                     lvar,   vname(2,indxGls)(1:lvar))
          lvar=lenstr(vname(3,indxGls))
          ierr=nf_put_att_text (ncid, rstGls, 'units', lvar,
     &         vname(3,indxGls)(1:lvar))
          lvar=lenstr(vname(4,indxGls))
          ierr=nf_put_att_text (ncid, rstGls, 'field', lvar,
     &         vname(4,indxGls)(1:lvar))

          call nf_add_attribute(ncid, rstGls, indxGls, 5,
     &         NF_DOUBLE,ierr)

!
! Vertical viscosity coefficient.
!

        lvar=lenstr(vname(1,indxAkv))
        ierr=nf_def_var (ncid, vname(1,indxAkv)(1:lvar),
     &       NF_DOUBLE, 4, w3dgrd, rstAkv)
#  ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstAkv,nf_collective)
#  endif
        lvar=lenstr(vname(2,indxAkv))
        ierr=nf_put_att_text (ncid, rstAkv, 'long_name',
     &       lvar,   vname(2,indxAkv)(1:lvar))
        lvar=lenstr(vname(3,indxAkv))
        ierr=nf_put_att_text (ncid, rstAkv, 'units', lvar,
     &       vname(3,indxAkv)(1:lvar))
        lvar=lenstr(vname(4,indxAkv))
        ierr=nf_put_att_text (ncid, rstAkv, 'field', lvar,
     &       vname(4,indxAkv)(1:lvar))
        
        call nf_add_attribute(ncid, rstAkv, indxAkv, 5,
     &       NF_DOUBLE,ierr)

#  ifdef TEMPERATURE        
!     
! Vertical diffusion coefficient for potential temperature.
!     
        
        lvar=lenstr(vname(1,indxAkt))
        ierr=nf_def_var (ncid, vname(1,indxAkt)(1:lvar),
     &       NF_DOUBLE, 4, w3dgrd, rstAkt)
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstAkt,nf_collective)
#   endif
        lvar=lenstr(vname(2,indxAkt))
        ierr=nf_put_att_text (ncid, rstAkt, 'long_name',
     &       lvar,   vname(2,indxAkt)(1:lvar))
        lvar=lenstr(vname(3,indxAkt))
        ierr=nf_put_att_text (ncid, rstAkt, 'units', lvar,
     &       vname(3,indxAkt)(1:lvar))
        lvar=lenstr(vname(4,indxAkt))
        ierr=nf_put_att_text (ncid, rstAkt, 'field', lvar,
     &       vname(4,indxAkt)(1:lvar))
        
        call nf_add_attribute(ncid, rstAkt, indxAkt, 5,
     &       NF_DOUBLE,ierr)
#  endif /* TEMPERATURE  */       

#  ifdef SALINITY
!
! Vertical diffusion coefficient for salinity.
!

          lvar=lenstr(vname(1,indxAks))
          ierr=nf_def_var (ncid, vname(1,indxAks)(1:lvar),
     &         NF_DOUBLE, 4, w3dgrd, rstAks)
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstAks,nf_collective)
#   endif
          lvar=lenstr(vname(2,indxAks))
          ierr=nf_put_att_text (ncid, rstAks, 'long_name',
     &                     lvar,   vname(2,indxAks)(1:lvar))
          lvar=lenstr(vname(3,indxAks))
          ierr=nf_put_att_text (ncid, rstAks, 'units', lvar,
     &         vname(3,indxAks)(1:lvar))
          lvar=lenstr(vname(4,indxAks))
          ierr=nf_put_att_text (ncid, rstAks, 'field', lvar,
     &         vname(4,indxAks)(1:lvar))

          call nf_add_attribute(ncid, rstAks, indxAks, 5,
     &         NF_DOUBLE,ierr)
#  endif /* SALINITY  */
# endif /* GLS_MIXING */
          
# if defined LMD_MIXING
#  ifdef M3FAST
!
!  Bottom stress u component 
!        
        lvar=lenstr(vname(1,indxBustr))
        ierr=nf_def_var (ncid, vname(1,indxBustr)(1:lvar),
     &       NF_DOUBLE, 3, u2dgrd, rstBustr)
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstBustr,nf_collective)
#   endif
        lvar=lenstr(vname(2,indxBustr))
        ierr=nf_put_att_text (ncid, rstBustr, 'long_name',
     &       lvar,   vname(2,indxBustr)(1:lvar))
        lvar=lenstr(vname(3,indxBustr))
        ierr=nf_put_att_text (ncid, rstBustr, 'units', lvar,
     &       vname(3,indxBustr)(1:lvar))
        lvar=lenstr(vname(4,indxBustr))
        ierr=nf_put_att_text (ncid, rstBustr, 'field', lvar,
     &       vname(4,indxBustr)(1:lvar))
        
        call nf_add_attribute(ncid, rstBustr, indxBustr, 5,
     &       NF_DOUBLE,ierr)
!
!  Bottom stress v component 
!        
        lvar=lenstr(vname(1,indxBvstr))
        ierr=nf_def_var (ncid, vname(1,indxBvstr)(1:lvar),
     &       NF_DOUBLE, 3, v2dgrd, rstBvstr)
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstBvstr,nf_collective)
#   endif
        lvar=lenstr(vname(2,indxBvstr))
        ierr=nf_put_att_text (ncid, rstBvstr, 'long_name',
     &       lvar,   vname(2,indxBvstr)(1:lvar))
        lvar=lenstr(vname(3,indxBvstr))
        ierr=nf_put_att_text (ncid, rstBvstr, 'units', lvar,
     &       vname(3,indxBvstr)(1:lvar))
        lvar=lenstr(vname(4,indxBvstr))
        ierr=nf_put_att_text (ncid, rstBvstr, 'field', lvar,
     &       vname(4,indxBvstr)(1:lvar))
        
        call nf_add_attribute(ncid, rstBvstr, indxBvstr, 5,
     &       NF_DOUBLE,ierr)
#  endif /* M3FAST */
# endif /* LMD_MIXING */

# ifdef EXACT_RESTART
!
! Forcing for barotropic equations in XI-direction.
!
        lvar=lenstr(vname(1,indxrufrc))
        ierr=nf_def_var (ncid, vname(1,indxrufrc)(1:lvar),
     &                         NF_DOUBLE, 3, u2dgrd, rstrufrc)
#  ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstrufrc,nf_collective)
#  endif
        lvar=lenstr(vname(2,indxrufrc))
        ierr=nf_put_att_text (ncid, rstrufrc, 'long_name', lvar,
     &                         vname(2,indxrufrc)(1:lvar))
        lvar=lenstr(vname(3,indxrufrc))
        ierr=nf_put_att_text (ncid, rstrufrc, 'units', lvar,
     &                         vname(3,indxrufrc)(1:lvar))
        lvar=lenstr(vname(4,indxrufrc))
        ierr=nf_put_att_text (ncid, rstrufrc, 'field', lvar,
     &                         vname(4,indxrufrc)(1:lvar))
      call nf_add_attribute(ncid, rstrufrc,
     &                         indxrufrc, 5, NF_DOUBLE, ierr)
!
! Forcing for barotropic equations in ETA-direction.
!
        lvar=lenstr(vname(1,indxrvfrc))
        ierr=nf_def_var (ncid, vname(1,indxrvfrc)(1:lvar),
     &                         NF_DOUBLE, 3, v2dgrd, rstrvfrc)
#  ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstrvfrc,nf_collective)
#  endif
        lvar=lenstr(vname(2,indxrvfrc))
        ierr=nf_put_att_text (ncid, rstrvfrc, 'long_name', lvar,
     &                         vname(2,indxrvfrc)(1:lvar))
        lvar=lenstr(vname(3,indxrvfrc))
        ierr=nf_put_att_text (ncid, rstrvfrc, 'units', lvar,
     &                         vname(3,indxrvfrc)(1:lvar))
        lvar=lenstr(vname(4,indxrvfrc))
        ierr=nf_put_att_text (ncid, rstrvfrc, 'field', lvar,
     &                         vname(4,indxrvfrc)(1:lvar))
      call nf_add_attribute(ncid, rstrvfrc,
     &                         indxrvfrc, 5, NF_DOUBLE, ierr)

#  ifdef M3FAST
!
! 3D rhs for M3FAST in XI- and ETA-directions.
!
        lvar=lenstr(vname(1,indxru_nbq))
        ierr=nf_def_var (ncid, vname(1,indxru_nbq)(1:lvar),
     &                              NF_DOUBLE, 4, u3dgrd, rstru_nbq)
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstru_nbq,nf_collective)
#   endif
        lvar=lenstr(vname(2,indxru_nbq))
        ierr=nf_put_att_text (ncid, rstru_nbq, 'long_name', lvar,
     &                                  vname(2,indxru_nbq)(1:lvar))
        lvar=lenstr(vname(3,indxru_nbq))
        ierr=nf_put_att_text (ncid, rstru_nbq, 'units',     lvar,
     &                                  vname(3,indxru_nbq)(1:lvar))
        lvar=lenstr(vname(4,indxru_nbq))
        ierr=nf_put_att_text (ncid, rstru_nbq, 'field',     lvar,
     &                                  vname(4,indxru_nbq)(1:lvar))

        call nf_add_attribute(ncid, rstru_nbq, indxru_nbq,
     &                                   5, NF_DOUBLE, ierr)

        lvar=lenstr(vname(1,indxrv_nbq))
        ierr=nf_def_var (ncid, vname(1,indxrv_nbq)(1:lvar),
     &                              NF_DOUBLE, 4, v3dgrd, rstrv_nbq)
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstrv_nbq,nf_collective)
#   endif
        lvar=lenstr(vname(2,indxrv_nbq))
        ierr=nf_put_att_text (ncid, rstrv_nbq, 'long_name', lvar,
     &                                  vname(2,indxrv_nbq)(1:lvar))
        lvar=lenstr(vname(3,indxrv_nbq))
        ierr=nf_put_att_text (ncid, rstrv_nbq, 'units',     lvar,
     &                                  vname(3,indxrv_nbq)(1:lvar))
        lvar=lenstr(vname(4,indxrv_nbq))
        ierr=nf_put_att_text (ncid, rstrv_nbq, 'field',     lvar,
     &                                  vname(4,indxrv_nbq)(1:lvar))

        call nf_add_attribute(ncid, rstrv_nbq, indxrv_nbq,
     &                                   5, NF_DOUBLE, ierr)

!
! 3D rhs for M3FAST in XI- and ETA-directions.
!
        lvar=lenstr(vname(1,indxru_nbq_avg2))
        ierr=nf_def_var (ncid, vname(1,indxru_nbq_avg2)(1:lvar),
     &                    NF_DOUBLE, 4, u3dgrd, rstru_nbq_avg2)
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstru_nbq_avg2,nf_collective)
#   endif
        lvar=lenstr(vname(2,indxru_nbq_avg2))
        ierr=nf_put_att_text (ncid, rstru_nbq_avg2, 'long_name',
     &                      lvar, vname(2,indxru_nbq_avg2)(1:lvar))
        lvar=lenstr(vname(3,indxru_nbq_avg2))
        ierr=nf_put_att_text (ncid, rstru_nbq_avg2, 'units', lvar,
     &                       vname(3,indxru_nbq_avg2)(1:lvar))
        lvar=lenstr(vname(4,indxru_nbq_avg2))
        ierr=nf_put_att_text (ncid, rstru_nbq_avg2, 'field', lvar,
     &                       vname(4,indxru_nbq_avg2)(1:lvar))

        call nf_add_attribute(ncid, rstru_nbq_avg2, indxru_nbq_avg2,
     &                                   5, NF_DOUBLE, ierr)

        lvar=lenstr(vname(1,indxrv_nbq_avg2))
        ierr=nf_def_var (ncid, vname(1,indxrv_nbq_avg2)(1:lvar),
     &                    NF_DOUBLE, 4, v3dgrd, rstrv_nbq_avg2)
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstrv_nbq_avg2,nf_collective)
#   endif
        lvar=lenstr(vname(2,indxrv_nbq_avg2))
        ierr=nf_put_att_text (ncid, rstrv_nbq_avg2, 'long_name',
     &                      lvar, vname(2,indxrv_nbq_avg2)(1:lvar))
        lvar=lenstr(vname(3,indxrv_nbq_avg2))
        ierr=nf_put_att_text (ncid, rstrv_nbq_avg2, 'units',  lvar,
     &                       vname(3,indxrv_nbq_avg2)(1:lvar))
        lvar=lenstr(vname(4,indxrv_nbq_avg2))
        ierr=nf_put_att_text (ncid, rstrv_nbq_avg2, 'field',  lvar,
     &                       vname(4,indxrv_nbq_avg2)(1:lvar))

        call nf_add_attribute(ncid, rstrv_nbq_avg2, indxrv_nbq_avg2,
     &                                   5, NF_DOUBLE, ierr)

!
! 3D rhs for M3FAST in XI- and ETA-directions.
!
        lvar=lenstr(vname(1,indxqdmu_nbq))
        ierr=nf_def_var (ncid, vname(1,indxqdmu_nbq)(1:lvar),
     &                    NF_DOUBLE, 4, u3dgrd, rstqdmu_nbq)
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstqdmu_nbq,nf_collective)
#   endif
        lvar=lenstr(vname(2,indxqdmu_nbq))
        ierr=nf_put_att_text (ncid, rstqdmu_nbq, 'long_name',
     &                      lvar, vname(2,indxqdmu_nbq)(1:lvar))
        lvar=lenstr(vname(3,indxqdmu_nbq))
        ierr=nf_put_att_text (ncid, rstqdmu_nbq, 'units', lvar,
     &                       vname(3,indxqdmu_nbq)(1:lvar))
        lvar=lenstr(vname(4,indxqdmu_nbq))
        ierr=nf_put_att_text (ncid, rstqdmu_nbq, 'field', lvar,
     &                       vname(4,indxqdmu_nbq)(1:lvar))

        call nf_add_attribute(ncid, rstqdmu_nbq, indxqdmu_nbq,
     &                                   5, NF_DOUBLE, ierr)

        lvar=lenstr(vname(1,indxqdmv_nbq))
        ierr=nf_def_var (ncid, vname(1,indxqdmv_nbq)(1:lvar),
     &                    NF_DOUBLE, 4, v3dgrd, rstqdmv_nbq)
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstqdmv_nbq,nf_collective)
#   endif
        lvar=lenstr(vname(2,indxqdmv_nbq))
        ierr=nf_put_att_text (ncid, rstqdmv_nbq, 'long_name',
     &                      lvar, vname(2,indxqdmv_nbq)(1:lvar))
        lvar=lenstr(vname(3,indxqdmv_nbq))
        ierr=nf_put_att_text (ncid, rstqdmv_nbq, 'units',  lvar,
     &                       vname(3,indxqdmv_nbq)(1:lvar))
        lvar=lenstr(vname(4,indxqdmv_nbq))
        ierr=nf_put_att_text (ncid, rstqdmv_nbq, 'field',  lvar,
     &                       vname(4,indxqdmv_nbq)(1:lvar))

        call nf_add_attribute(ncid, rstqdmv_nbq, indxqdmv_nbq,
     &                                   5, NF_DOUBLE, ierr)
#  endif  /* M3FAST */

#  ifdef TRACERS
#   ifdef TS_MIX_ISO_FILT
!
! density gradients for use in t3dmix
!
        lvar=lenstr(vname(1,indxdRdx))
        ierr=nf_def_var (ncid, vname(1,indxdRdx)(1:lvar),
     &                    NF_DOUBLE, 4, u3dgrd, rstdRdx)
#    ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstdRdx,nf_collective)
#    endif
        lvar=lenstr(vname(2,indxdRdx))
        ierr=nf_put_att_text (ncid, rstdRdx, 'long_name',
     &                      lvar, vname(2,indxdRdx)(1:lvar))
        lvar=lenstr(vname(3,indxdRdx))
        ierr=nf_put_att_text (ncid, rstdRdx, 'units', lvar,
     &                       vname(3,indxdRdx)(1:lvar))
        lvar=lenstr(vname(4,indxdRdx))
        ierr=nf_put_att_text (ncid, rstdRdx, 'field', lvar,
     &                       vname(4,indxdRdx)(1:lvar))

        call nf_add_attribute(ncid, rstdRdx, indxdRdx,
     &                                   5, NF_DOUBLE, ierr)

        lvar=lenstr(vname(1,indxdRde))
        ierr=nf_def_var (ncid, vname(1,indxdRde)(1:lvar),
     &                    NF_DOUBLE, 4, v3dgrd, rstdRde)
#    ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstdRde,nf_collective)
#    endif
        lvar=lenstr(vname(2,indxdRde))
        ierr=nf_put_att_text (ncid, rstdRde, 'long_name',
     &                      lvar, vname(2,indxdRde)(1:lvar))
        lvar=lenstr(vname(3,indxdRde))
        ierr=nf_put_att_text (ncid, rstdRde, 'units',  lvar,
     &                       vname(3,indxdRde)(1:lvar))
        lvar=lenstr(vname(4,indxdRde))
        ierr=nf_put_att_text (ncid, rstdRde, 'field',  lvar,
     &                       vname(4,indxdRde)(1:lvar))

        call nf_add_attribute(ncid, rstdRde, indxdRde,
     &                                   5, NF_DOUBLE, ierr)

#   endif /* TS_MIX_ISO_FILT */
#  endif /* TRACERS */
# endif  /* EXACT_RESTART */
#endif /* SOLVE3D */

#ifdef SEDIMENT
!
! Sediment bed thickness
!
        lvar=lenstr(vname(1,indxBTHK))
        ierr=nf_def_var (ncid, vname(1,indxBTHK)(1:lvar),
     &                              NF_DOUBLE, 4, b3dgrd, rstSed(1))
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstSed(1),nf_collective)
# endif
        lvar=lenstr(vname(2,indxBTHK))
        ierr=nf_put_att_text (ncid, rstSed(1), 'long_name', lvar,
     &                                   vname(2,indxBTHK)(1:lvar))
        lvar=lenstr(vname(3,indxBTHK))
        ierr=nf_put_att_text (ncid, rstSed(1), 'units',     lvar,
     &                                   vname(3,indxBTHK)(1:lvar))
        lvar=lenstr(vname(4,indxBTHK))
        ierr=nf_put_att_text (ncid, rstSed(1), 'field',     lvar,
     &                                  vname(4,indxBTHK)(1:lvar))

      call nf_add_attribute(ncid, rstSed(1), indxBTHK, 5, 
     &       NF_DOUBLE,ierr)

!
! Sediment bed porosity
!
        lvar=lenstr(vname(1,indxBPOR))
        ierr=nf_def_var (ncid, vname(1,indxBPOR)(1:lvar),
     &                              NF_DOUBLE, 4, b3dgrd, rstSed(2))
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstSed(2),nf_collective)
# endif
        lvar=lenstr(vname(2,indxBPOR))
        ierr=nf_put_att_text (ncid, rstSed(2), 'long_name', lvar,
     &                                   vname(2,indxBPOR)(1:lvar))
        lvar=lenstr(vname(3,indxBPOR))
        ierr=nf_put_att_text (ncid, rstSed(2), 'units',     lvar,
     &                                   vname(3,indxBPOR)(1:lvar))
        lvar=lenstr(vname(4,indxBPOR))
        ierr=nf_put_att_text (ncid, rstSed(2), 'field',     lvar,
     &                                  vname(4,indxBPOR)(1:lvar))
        
      call nf_add_attribute(ncid, rstSed(2), indxBPOR, 5, 
     &       NF_DOUBLE, ierr)
   
!
! Bed size class fractions.
!
        do itrc=1,NST
          indxWrk=indxBFRA(1)+itrc-1
          lvar=lenstr(vname(1,indxWrk))
          ierr=nf_def_var (ncid, vname(1,indxWrk)(1:lvar),
     &                          NF_DOUBLE, 4, b3dgrd, rstSed(itrc+2))
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,rstSed(itrc+2),nf_collective)
# endif
          lvar=lenstr(vname(2,indxWrk))
          ierr=nf_put_att_text (ncid, rstSed(itrc+2), 'long_name',
     &                     lvar,   vname(2,indxWrk)(1:lvar))
          lvar=lenstr(vname(3,indxWrk))
          ierr=nf_put_att_text (ncid, rstSed(itrc+2), 'units', lvar,
     &                             vname(3,indxWrk)(1:lvar))
          lvar=lenstr(vname(4,indxWrk))
          ierr=nf_put_att_text (ncid, rstSed(itrc+2), 'field', lvar,
     &                             vname(4,indxWrk)(1:lvar))

      call nf_add_attribute(ncid, rstSed(itrc+2), indxWrk, 5, 
     &         NF_DOUBLE, ierr)

        enddo
#endif /* SEDIMENT */

#ifdef MUSTANG
! ksmi, ksma
       do itrc=1,2
       indx=indxT+ntrc_salt+ntrc_substot+ntrc_subs+itrc+6
       lvar=lenstr(vname(1,indx))
       ierr=nf_def_var (ncid, vname(1,indx)(1:lvar),
     &            NF_INT, 3, r2dgrd, rstMUS(itrc)) ! ksmi,ksma

# ifdef NC4PAR
       ierr=nf_var_par_access(ncid,rstMUS(itrc),nf_collective)
# endif

        call nf_add_attribute(ncid, rstMUS(itrc), indx, 5,
     &         NF_INT, ierr)
       enddo

! dzs, tempsed, salsed
      do itrc=1,3
       indx=indxT+ntrc_salt+ntrc_substot+itrc+3
       ierr=nf_def_var (ncid, vname(1,indx),
     &            NF_DOUBLE, 4, b3dgrd, rstMUS(itrc+2)) ! 3D: dzs,tempsed,salsed

#  ifdef NC4PAR
       ierr=nf_var_par_access(ncid,rstMUS(itrc+2),nf_collective)
#  endif

        call nf_add_attribute(ncid, rstMUS(itrc+2), indx, 5,
     &         NF_DOUBLE, ierr)

      enddo

!cvsed
       do itrc=1,ntrc_subs
        indx=indxT+ntrc_salt+ntrc_substot+itrc+6           ! 3D: traceurSed
        ierr=nf_def_var (ncid, vname(1,indx),
     &            NF_DOUBLE, 4, b3dgrd, rstMUS(itrc+5)) ! 3D: traceurSed
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstMUS(itrc+5),nf_collective)
# endif

        call nf_add_attribute(ncid, rstMUS(itrc+5), indx, 5,
     &         NF_DOUBLE, ierr)

        enddo
#endif


#ifdef MORPHODYN 
        lvar=lenstr(vname(1,indxHm))
        ierr=nf_def_var (ncid, vname(1,indxHm)(1:lvar),
     &                            NF_DOUBLE, 3, r2dgrd, rstHm)
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstHm,nf_collective)
# endif
        lvar=lenstr(vname(2,indxHm))
        ierr=nf_put_att_text (ncid, rstHm, 'long_name', lvar,
     &                                vname(2,indxHm)(1:lvar))
        lvar=lenstr(vname(3,indxHm))
        ierr=nf_put_att_text (ncid, rstHm, 'units',     lvar,
     &                                vname(3,indxHm)(1:lvar))
        lvar=lenstr(vname(4,indxHm))
        ierr=nf_put_att_text (ncid, rstHm, 'field',     lvar,
     &                                vname(4,indxHm)(1:lvar))

        call nf_add_attribute(ncid, rstHm, indxHm, 5, NF_DOUBLE, ierr)

#endif

#ifdef BBL
!
! Bed ripple height
!
        lvar=lenstr(vname(1,indxHrip))
        ierr=nf_def_var (ncid, vname(1,indxHrip)(1:lvar),
     &                              NF_DOUBLE, 3, r2dgrd, rstBBL(1))
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstBBL(1),nf_collective)
# endif
        lvar=lenstr(vname(2,indxHrip))
        ierr=nf_put_att_text (ncid, rstBBL(1), 'long_name', lvar,
     &                                   vname(2,indxHrip)(1:lvar))
        lvar=lenstr(vname(3,indxHrip))
        ierr=nf_put_att_text (ncid, rstBBL(1), 'units',     lvar,
     &                                   vname(3,indxHrip)(1:lvar))
        lvar=lenstr(vname(4,indxHrip))
        ierr=nf_put_att_text (ncid, rstBBL(1), 'field',     lvar,
     &                                  vname(4,indxHrip)(1:lvar))

        call nf_add_attribute(ncid, rstBBL(1), indxHrip, 5, 
     &    NF_DOUBLE, ierr) 

 
!
! Bed ripple length
!
        lvar=lenstr(vname(1,indxLrip))
        ierr=nf_def_var (ncid, vname(1,indxLrip)(1:lvar),
     &                              NF_DOUBLE, 3, r2dgrd, rstBBL(2))
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstBBL(2),nf_collective)
# endif
        lvar=lenstr(vname(2,indxLrip))
        ierr=nf_put_att_text (ncid, rstBBL(2), 'long_name', lvar,
     &                                   vname(2,indxLrip)(1:lvar))
        lvar=lenstr(vname(3,indxLrip))
        ierr=nf_put_att_text (ncid, rstBBL(2), 'units',     lvar,
     &                                   vname(3,indxLrip)(1:lvar))
        lvar=lenstr(vname(4,indxLrip))
        ierr=nf_put_att_text (ncid, rstBBL(2), 'field',     lvar,
     &                                  vname(4,indxLrip)(1:lvar))

        call nf_add_attribute(ncid, rstBBL, indxLrip, 5, 
     &       NF_DOUBLE, ierr)
#endif /* BBL */


#ifdef FLOATS
# ifdef AGRIF
      if (Agrif_Root()) then
# endif

           vinfo(1)='nfloats'
           vinfo(2)='total number of floats'

           lvar=lenstr(vinfo(1))
           ierr=nf_def_var(ncid,vinfo(1)(1:lvar),nf_int,
     &                    1,timedim,rstnfloats)
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid,rstnfloats,nf_collective)
# endif
           lvar=lenstr(vinfo(2))
           ierr=nf_put_att_text(ncid,rstnfloats,'long_name',lvar,
     &                         vinfo(2)(1:lvar))

! define Tinfo

           vinfo(1)='Tinfo'
           vinfo(2)='initial float release information array'

           lvar=lenstr(vinfo(1))
           ierr=nf_def_var(ncid,vinfo(1)(1:lvar),NF_DOUBLE,
     &                    3,dim_Tinfo,rstTinfo)
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid,rstinfo,nf_collective)
# endif
           lvar=lenstr(vinfo(2))
           ierr=nf_put_att_text(ncid,rstTinfo,'long_name',lvar,
     &                         vinfo(2)(1:lvar))

! define grid level

           vinfo(1)='grid_level'
           vinfo(2)='float position in nested grids hierarchy'

           lvar=lenstr(vinfo(1))
           ierr=nf_def_var(ncid,vinfo(1)(1:lvar),nf_int,
     &                    2,dim_fltgrd,rstfltgrd)
           lvar=lenstr(vinfo(2))
           ierr=nf_put_att_text(ncid,rstfltgrd,'long_name',lvar,
     &                         vinfo(2)(1:lvar))

! define track

           vinfo(1)='track'
           vinfo(2)='Position and velocity array for released floats'
           vinfo(3)='x,y,z,u,v,w in this order'
           lvar=lenstr(vinfo(1))
           ierr=nf_def_var(ncid,vinfo(1)(1:lvar),NF_DOUBLE,
     &                    4,dim_track,rsttrack)
           lvar=lenstr(vinfo(2))
           ierr=nf_put_att_text(ncid,rsttrack,'long_name',lvar,
     &                         vinfo(2)(1:lvar))
           lvar=lenstr(vinfo(3))
           ierr=nf_put_att_text(ncid,rsttrack,'Content',lvar,
     &                         vinfo(3)(1:lvar))

# ifdef AGRIF
      endif
# endif
#endif /* FLOATS */

!
! Leave definition mode.                  Also initialize record
! ----- ---------- -----                  dimension size to zero.
!
        ierr=nf_enddef(ncid)
        MPI_master_only write(*,'(6x,4A,1x,A,i4)')
     &              'DEF_RST - Created new ',
     &              'netCDF file ''', rstname(1:lstr), '''.'
     &               MYID
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
      elseif (ncid.eq.-1) then
#ifndef NC4PAR
        ierr=nf_open (rstname(1:lstr), nf_write, ncid)
#else
        ierr=nf_open_par (rstname(1:lstr), IOR(nf_write, nf_mpiio),
     &     MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
#endif
        if (ierr. eq. nf_noerr) then
          ierr=checkdims (ncid, rstname, lstr, rec)
          if (ierr .eq. nf_noerr) then
            if (nrpfrst.eq.0) then
              ierr=rec+1 - nrecrst
            else
              ierr=rec+1 - (1+mod(nrecrst-1, abs(nrpfrst)))
            endif
            if (ierr.gt.0) then
              MPI_master_only write( stdout,
     &                 '(/1x,A,I5,1x,A/8x,3A,I5/8x,A,I5,1x,A/)'
     &           ) 'DEF_RST WARNING: Actual number of records', rec,
     &             'in netCDF file',  '''',  rstname(1:lstr),
     &             ''' exceeds the record number from restart data',
     &             rec+1-ierr,'/', total_rec,', restart is assumed.'
              rec=rec-ierr
            elseif (nrpfrst.eq.0) then
              total_rec=rec+1           ! <-- set to the next record
#if defined MPI & !defined PARALLEL_FILES
              if (mynode.gt.0) total_rec=total_rec-1
#endif
            endif
            ierr=nf_noerr
          endif
        endif

        if (ierr. ne. nf_noerr) then
#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
          if (mynode.eq.0) then
            create_new_file=.true.
            goto 10
          else
            write(stdout,'(/1x,4A, 1x,A,I4/)')     'DEF_RST ERROR: ',
     &     'Cannot open restart netCDF file ''',rstname(1:lstr),'''.'
     &               MYID
            goto 99                                     !--> ERROR 
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
        ierr=nf_inq_varid (ncid, 'time_step', rstTstep)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) 'time_step', rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstTstep,nf_collective)
#endif
!
! Time.
!
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid, vname(1,indxTime)(1:lvar), rstTime)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime)(1:lvar), rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstTime,nf_collective)
#endif
!
! Time2.
!
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_inq_varid (ncid, vname(1,indxTime2)(1:lvar), rstTime2)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime2)(1:lvar), rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstTime2,nf_collective)
#endif

! Free-surface.
!
        lvar=lenstr(vname(1,indxZ))
        ierr=nf_inq_varid (ncid, vname(1,indxZ)(1:lvar), rstZ)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxZ)(1:lvar),  rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstZ,nf_collective)
#endif
!
! 2D momenta in XI- and ETA-directions.
!
        lvar=lenstr(vname(1,indxUb))
        ierr=nf_inq_varid (ncid, vname(1,indxUb)(1:lvar), rstUb)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxUb)(1:lvar), rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstUb,nf_collective)
#endif

        lvar=lenstr(vname(1,indxVb))
        ierr=nf_inq_varid (ncid, vname(1,indxVb)(1:lvar), rstVb)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxVb)(1:lvar), rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstVb,nf_collective)
#endif

#ifdef SOLVE3D
!
! 3D momenta n XI- and ETA-directions.
!
        lvar=lenstr(vname(1,indxU))
        ierr=nf_inq_varid (ncid, vname(1,indxU)(1:lvar), rstU)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxU)(1:lvar), rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstU,nf_collective)
#endif

        lvar=lenstr(vname(1,indxV))
        ierr=nf_inq_varid (ncid, vname(1,indxV)(1:lvar), rstV)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxV)(1:lvar), rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstV,nf_collective)
#endif

!
! Tracer variables.
!
# if defined TRACERS
        do itrc=1,NT
          lvar=lenstr(vname(1,indxV+itrc))
          ierr=nf_inq_varid (ncid, vname(1,indxV+itrc)(1:lvar),
     &                                                rstT(itrc))
          if (ierr. ne. nf_noerr) then
            write(stdout,1) vname(1,indxV+itrc)(1:lvar),
     &                                            rstname(1:lstr)
            goto 99                                       !--> ERROR
          endif
#  ifdef NC4PAR
          ierr=nf_var_par_access(ncid,rstT(itrc),nf_collective)
#  endif
        enddo 
# endif

# if defined LMD_MIXING 
!
! Depth of planetary boundary layer.
!
#  ifdef LMD_SKPP
          lvar=lenstr(vname(1,indxHbl))
          ierr=nf_inq_varid (ncid,vname(1,indxHbl)(1:lvar), rstHbl)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxHbl)(1:lvar), rstname(1:lstr)
            goto 99                                       !--> ERROR
          endif
#   ifdef NC4PAR
          ierr=nf_var_par_access(ncid,rstHbl,nf_collective)
#   endif
#  endif
          
#  ifdef LMD_BKPP
!
! Depth of bottom planetary boundary layer.
!
          lvar=lenstr(vname(1,indxHbbl))
          ierr=nf_inq_varid (ncid,vname(1,indxHbbl)(1:lvar), rstHbbl)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxHbbl)(1:lvar), rstname(1:lstr)
            goto 99                                       !--> ERROR
          endif
#   ifdef NC4PAR
          ierr=nf_var_par_access(ncid,rstHbbl,nf_collective)
#   endif

#  endif /*  LMD_BKPP  */
# endif /* LMD_MIXING */

# if defined GLS_MIXING
!
! Turbulent kinetic energy.
!
          lvar=lenstr(vname(1,indxTke))
          ierr=nf_inq_varid (ncid,vname(1,indxTke)(1:lvar), rstTke)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTke)(1:lvar), rstname(1:lstr)
            goto 99                                       !--> ERROR
          endif
#  ifdef NC4PAR
          ierr=nf_var_par_access(ncid,rstTke,nf_collective)
#  endif

!
! Generic length scale.
!
          lvar=lenstr(vname(1,indxGls))
          ierr=nf_inq_varid (ncid,vname(1,indxGls)(1:lvar), rstGls)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxGls)(1:lvar), rstname(1:lstr)
            goto 99                                       !--> ERROR
          endif
#  ifdef NC4PAR
          ierr=nf_var_par_access(ncid,rstGls,nf_collective)
#  endif

!
! Vertical viscosity coefficient.
!
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_inq_varid (ncid,vname(1,indxAkv)(1:lvar), rstAkv)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxAkv)(1:lvar), rstname(1:lstr)
            goto 99                                       !--> ERROR
          endif
#  ifdef NC4PAR
          ierr=nf_var_par_access(ncid,rstAkv,nf_collective)
#   endif

#  ifdef TEMPERATURE

!
! Vertical diffusion coefficient for potential temperature.
!
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_inq_varid (ncid,vname(1,indxAkt)(1:lvar), rstAkt)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxAkt)(1:lvar), rstname(1:lstr)
            goto 99                                       !--> ERROR
          endif
#   ifdef NC4PAR
          ierr=nf_var_par_access(ncid,rstAkt,nf_collective)
#   endif

#  endif
          
#  ifdef SALINITY
!
! Vertical diffusion coefficient for salinity.
!
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_inq_varid (ncid,vname(1,indxAks)(1:lvar), rstAks)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxAks)(1:lvar), rstname(1:lstr)
            goto 99                                       !--> ERROR
          endif
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstAks,nf_collective)
#   endif
#  endif /* SALINITY  */
# endif /* GLS_MIXING */


# ifdef M3FAST
#  if defined LMD_MIXING
!
! Bottom stress u component
!
          lvar=lenstr(vname(1,indxBustr))
          ierr=nf_inq_varid (ncid,vname(1,indxBustr)(1:lvar), rstBustr)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxBustr)(1:lvar), rstname(1:lstr)
            goto 99                                       !--> ERROR
          endif
#   ifdef NC4PAR
          ierr=nf_var_par_access(ncid,rstBustr,nf_collective)
#   endif

!
! Bottom stress v component
!
          lvar=lenstr(vname(1,indxBvstr))
          ierr=nf_inq_varid (ncid,vname(1,indxBvstr)(1:lvar), rstBvstr)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxBvstr)(1:lvar), rstname(1:lstr)
            goto 99                                       !--> ERROR
          endif
#   ifdef NC4PAR
          ierr=nf_var_par_access(ncid,rstBvstr,nf_collective)
#   endif
#  endif
# endif
                    
# ifdef EXACT_RESTART
!
! Forcing for barotropic equations in XI-direction.
!
        lvar=lenstr(vname(1,indxrufrc))
        ierr=nf_inq_varid (ncid, vname(1,indxrufrc)(1:lvar), rstrufrc)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxrufrc)(1:lvar), rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#  ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstrufrc,nf_collective)
#  endif

!
! Forcing for barotropic equations in XI-direction.
!
        lvar=lenstr(vname(1,indxrvfrc))
        ierr=nf_inq_varid (ncid, vname(1,indxrvfrc)(1:lvar), rstrvfrc)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxrvfrc)(1:lvar), rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#  ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstrvfrc,nf_collective)
#  endif

#  ifdef M3FAST
        
!
! 3D rhs for M3FAST in XI- and ETA-directions.
!

        lvar=lenstr(vname(1,indxru_nbq))
        ierr=nf_inq_varid (ncid, vname(1,indxru_nbq)(1:lvar),
     &                                                 rstru_nbq)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxru_nbq)(1:lvar),
     &                                                 rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstru_nbq,nf_collective)
#   endif

        lvar=lenstr(vname(1,indxrv_nbq))
        ierr=nf_inq_varid (ncid, vname(1,indxrv_nbq)(1:lvar),
     &                                                 rstrv_nbq)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxrv_nbq)(1:lvar),
     &                                                 rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstrv_nbq,nf_collective)
#   endif

!
! 3D rhs for M3FAST in XI- and ETA-directions.
!

        lvar=lenstr(vname(1,indxru_nbq_avg2))
        ierr=nf_inq_varid (ncid, vname(1,indxru_nbq_avg2)(1:lvar),
     &                                            rstru_nbq_avg2)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxru_nbq_avg2)(1:lvar),
     &                                                 rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstru_nbq_avg2,nf_collective)
#   endif


        lvar=lenstr(vname(1,indxrv_nbq_avg2))
        ierr=nf_inq_varid (ncid, vname(1,indxrv_nbq_avg2)(1:lvar),
     &                                            rstrv_nbq_avg2)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxrv_nbq_avg2)(1:lvar),
     &                                                 rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstrv_nbq_avg2,nf_collective)
#   endif

!
! 3D rhs for M3FAST in XI- and ETA-directions.
!

        lvar=lenstr(vname(1,indxqdmu_nbq))
        ierr=nf_inq_varid (ncid, vname(1,indxqdmu_nbq)(1:lvar),
     &                                            rstqdmu_nbq)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxqdmu_nbq)(1:lvar),
     &                                                 rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstqdmu_nbq,nf_collective)
#   endif


        lvar=lenstr(vname(1,indxqdmv_nbq))
        ierr=nf_inq_varid (ncid, vname(1,indxqdmv_nbq)(1:lvar),
     &                                            rstqdmv_nbq)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxqdmv_nbq)(1:lvar),
     &                                                 rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#   ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstqdmv_nbq,nf_collective)
#   endif

#  endif  /* M3FAST */

#  ifdef TRACERS
#   ifdef TS_MIX_ISO_FILT
!
! density gradients for use in t3dmix
!

        lvar=lenstr(vname(1,indxdRdx))
        ierr=nf_inq_varid (ncid, vname(1,indxdRdx)(1:lvar),
     &                                            rstdRdx)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxdRdx)(1:lvar),
     &                                                 rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#    ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstdRdx,nf_collective)
#    endif

        lvar=lenstr(vname(1,indxdRde))
        ierr=nf_inq_varid (ncid, vname(1,indxdRde)(1:lvar),
     &                                            rstdRde)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxdRde)(1:lvar),
     &                                                 rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
#    ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstdRde,nf_collective)
#    endif

#   endif /* TS_MIX_ISO_FILT */
#  endif /* TRACERS */
# endif /* EXACT_RESTART */
#endif /* SOLVE3D */

#ifdef SEDIMENT
!
! Bed thickness
!
        lvar=lenstr(vname(1,indxBTHK))
        ierr=nf_inq_varid (ncid, vname(1,indxBTHK)(1:lvar), rstSed(1))
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxBTHK)(1:lvar), rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstSed(1),nf_collective)
# endif

!
! Bed porosity
!
        lvar=lenstr(vname(1,indxBPOR))
        ierr=nf_inq_varid (ncid, vname(1,indxBPOR)(1:lvar), rstSed(2))
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxBPOR)(1:lvar), rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstSed(2),nf_collective)
# endif

!
! Bed size class fractions.
!
        do itrc=1,NST
          indxWrk=indxBFRA(1)+itrc-1
          lvar=lenstr(vname(1,indxWrk))
          ierr=nf_inq_varid (ncid, vname(1,indxWrk)(1:lvar),
     &                                                 rstSed(itrc+2))
          if (ierr. ne. nf_noerr) then
            write(stdout,1) vname(1,indxWrk)(1:lvar),
     &                                            rstname(1:lstr)
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,rstSed(itrc+2),nf_collective)
# endif
        enddo
#endif /* SEDIMENT */

#ifdef MUSTANG
! ksmi, ksma
        do itrc=1,2
          indx=indxT+ntrc_salt+ntrc_substot+ntrc_substot+itrc+6
          lvar=lenstr(vname(1,indx))
          ierr=nf_inq_varid (ncid, vname(1,indx)(1:lvar),
     &                                                 rstMUS(itrc))
          if (ierr. ne. nf_noerr) then
            write(stdout,1) vname(1,indx)(1:lvar),
     &                                            rstMUS(1:lstr)
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,rstMUS(itrc),nf_collective)
# endif
        enddo

! dzs, tempsed, salsed
       do itrc=1,3
         indx=indxT+ntrc_salt+ntrc_substot+itrc+3
          lvar=lenstr(vname(1,indx))
          ierr=nf_inq_varid (ncid, vname(1,indx)(1:lvar),
     &                                                 rstMUS(itrc+2))
          if (ierr. ne. nf_noerr) then
            write(stdout,1) vname(1,indx)(1:lvar),
     &                                            rstMUS(1:lstr)
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,rstMUS(itrc+2),nf_collective)
# endif
        enddo


! cvsed
       do itrc=1,ntrc_subs
          indx=indxT+ntrc_salt+ntrc_substot+itrc+6           ! 3D: traceurSed
          lvar=lenstr(vname(1,indx))
          ierr=nf_inq_varid (ncid, vname(1,indx)(1:lvar),
     &                                                 rstMUS(itrc+5))
          if (ierr. ne. nf_noerr) then
            write(stdout,1) vname(1,indx)(1:lvar),
     &                                            rstMUS(1:lstr)
            goto 99                                       !--> ERROR
          endif
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,rstMUS(itrc+5),nf_collective)
# endif
        enddo

#endif


#ifdef MORPHODYN
!
! Time evolving bathymetry
!
        lvar=lenstr(vname(1,indxHm))
        ierr=nf_inq_varid (ncid, vname(1,indxHm)(1:lvar), rstHm)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxHm)(1:lvar), rstname(1:lstr)
          goto 99                                       !--> ERROR
        endif
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstHm,nf_collective)
# endif
#endif

#ifdef BBL
!
! Bed ripple height
!
        lvar=lenstr(vname(1,indxHrip))
        ierr=nf_inq_varid (ncid, vname(1,indxHrip)(1:lvar), rstBBL(1))
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxHrip)(1:lvar), rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstBBL(1),nf_collective)
# endif
!
! Bed ripple lengths
!
        lvar=lenstr(vname(1,indxLrip))
        ierr=nf_inq_varid (ncid, vname(1,indxLrip)(1:lvar), rstBBL(2))
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxLrip)(1:lvar), rstname(1:lstr)
          goto 99                                         !--> ERROR
        endif
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,rstBBL(2),nf_collective)
# endif
#endif /* BBL */
!
      MPI_master_only write(*,'(6x,2A,i4,1x,A,i4)') 
     &             'DEF_RST -- Opened ',
     &             'existing restart file,  record =', rec 
     &              MYID

#if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      else
        ierr=nf_open (rstname(1:lstr), nf_write, ncid)
        if (ierr .ne. nf_noerr) then
         MPI_master_only write(stdout,'(/1x,4A, 1x,A,I4/)') 
     &          'DEF_RST ERROR: Cannot',
     &          'open restart netCDF file ''', rstname(1:lstr), '''.'
     &           MYID
          goto 99                                         !--> ERROR
        endif
#endif
      endif              !<-- create_new_file
   1  format(/1x,'DEF_RST ERROR: Cannot find variable ''',
     &               A, ''' in netCDF file ''', A, '''.'/)
#ifdef PUT_GRID_INTO_RESTART
!
! Write grid variables.
! ----- ---- ----------
!
      if (total_rec.le.1) call wrt_grid (ncid, rstname, lstr)
#endif
  99  return                                              !--> ERROR
      end
