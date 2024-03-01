! $Id: def_surf.F 1468 2014-02-10 12:58:35Z rblod $
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
# if defined OUTPUTS_SURFACE && ! defined XIOS

#ifndef AVRH
      subroutine def_surf (ncid, total_rec, ierr)
#else
      subroutine def_surf_avg(ncid, total_rec, ierr)
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
# include "surf.h"
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
#ifndef AVRH

# define ncname surfname
# define rec_per_file nrpfsurf
# define wrt wrtsurf
# define vidTime surfTime
# define vidTime2 surfTime2
# define vidTstep surfTstep

# define vidsurft surf_surft
# define vidsurfs surf_surfs
# define vidsurfz surf_surfz
# define vidsurfu surf_surfu
# define vidsurfv surf_surfv


#else

# define ncname surfname_avg
# define rec_per_file nrpfsurf_avg
# define wrt wrtsurf_avg
# define vidTime surfTime_avg
# define vidTime2 surfTime2_avg
# define vidTstep surfTstep_avg

# define vidsurft surf_surft_avg
# define vidsurfs surf_surfs_avg
# define vidsurfz surf_surfz_avg
# define vidsurfu surf_surfu_avg
# define vidsurfv surf_surfv_avg

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
# ifdef DO_NOT_OVERWRITE
      if (lvar.eq.total_rec-1) then
#  ifndef AVRH
        create_new_file=ldefsurf
#  else
        create_new_file=ldefsurf_avg
#  endif
      else
#  ifndef AVRH
        create_new_file=.false.
#  else
        create_new_file=.false.
#  endif
      endif
# else
#  ifndef AVRH
      create_new_file=ldefsurf
#  else
      create_new_file=ldefsurf_avg
#  endif
# endif
!
      if (ncid.ne.-1) create_new_file=.false.
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode.gt.0) create_new_file=.false.
# endif
!
 10   if (create_new_file) then
# ifndef NC4PAR
        ierr=nf_create(ncname(1:lstr), NF_CLOBBER, ncid)
# else
        cmode = ior(nf_netcdf4, nf_classic_model)
        cmode = ior(cmode, nf_mpiio)
        csize=xi_rho*eta_rho/NNODES
        MPI_master_only write(stdout,*)
     &       'CREATE SURF HIS/AVG NC4 PARALLEL FILE'
        ierr=nf_create_par(ncname(1:lstr), cmode,
     &        MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
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
     &                   NF_DOUBLE, 1, timedim, vidTime)
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid, vidTime, nf_collective)
# endif
# ifndef AVRH
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', 
     &                        lvar, vname(2,indxTime)(1:lvar))
# else
        text='avg'/ /vname(2,indxTime)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid,  vidTime, 'long_name', 
     &                        lvar, text(1:lvar))
# endif
!
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', lvar,
     &                                vname(2,indxTime)(1:lvar))
!
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid,  vidTime, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))
!
        lvar=lenstr(vname(4,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'field',  lvar,
     &                                vname(4,indxTime)(1:lvar))

        call nf_add_attribute(ncid, vidTime, indxTime, 5, 
     &                        NF_FOUT, ierr)
!
! Time2.
!
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_def_var (ncid, vname(1,indxTime2)(1:lvar),
     &                          NF_DOUBLE, 1, timedim, vidTime2)
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid, vidTime2, nf_collective)
# endif
# ifndef AVRH
        lvar=lenstr(vname(2,indxTime2))
        ierr=nf_put_att_text (ncid, vidTime2, 'long_name', 
     &                        lvar, vname(2,indxTime)(1:lvar))
# else
        text='avg'/ /vname(2,indxTime2)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidTime2, 'long_name', 
     &                        lvar, text(1:lvar))
# endif
        lvar=lenstr(vname(2,indxTime2))
        ierr=nf_put_att_text (ncid, vidTime2, 'long_name', lvar,
     &                                vname(2,indxTime2)(1:lvar))
        lvar=lenstr(vname(3,indxTime2))
        ierr=nf_put_att_text (ncid, vidTime2, 'units',  lvar,
     &                                vname(3,indxTime2)(1:lvar))
        lvar=lenstr(vname(4,indxTime2))
        ierr=nf_put_att_text (ncid, vidTime2, 'field',  lvar,
     &                                vname(4,indxTime2)(1:lvar))

        call nf_add_attribute(ncid, vidTime2, indxTime2, 5, 
     &                        NF_FOUT, ierr)

! surface variables.
! -------- ----------- --------
!
          itrc=1
!
          if (wrt(itrc)) then


!------------------------------
!
! indxsurft
!          
          lvar=lenstr(vname(1,indxsurft+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxsurft+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, vidsurft(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid, vidsurft(itrc),
     &                                   nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxsurft+itrc-1)
# else
          text=vname(2,indxsurft+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidsurft(itrc), 'long_name',
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxsurft+itrc-1))
          ierr=nf_put_att_text (ncid, vidsurft(itrc), 'units', lvar,
     &                          vname(3,indxsurft+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxsurft+itrc-1))
          ierr=nf_put_att_text (ncid,vidsurft(itrc), 'field',
     &                      lvar, vname(4,indxsurft+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, vidsurft(itrc),
     &                       indxsurft+itrc-1, 5,
     &                        NF_FOUT, ierr)


!------------------------------
!
! indxsurfs
!          
          lvar=lenstr(vname(1,indxsurfs+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxsurfs+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, vidsurfs(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid, vidsurfs(itrc),
     &                                   nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxsurfs+itrc-1)
# else
          text=vname(2,indxsurfs+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidsurfs(itrc), 'long_name',
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxsurfs+itrc-1))
          ierr=nf_put_att_text (ncid, vidsurfs(itrc), 'units', lvar,
     &                          vname(3,indxsurfs+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxsurfs+itrc-1))
          ierr=nf_put_att_text (ncid,vidsurfs(itrc), 'field',
     &                      lvar, vname(4,indxsurfs+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, vidsurfs(itrc),
     &                       indxsurfs+itrc-1, 5,
     &                        NF_FOUT, ierr)

!------------------------------
!
! indxsurfz
!          
          lvar=lenstr(vname(1,indxsurfz+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxsurfz+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, vidsurfz(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid, vidsurfz(itrc),
     &                                   nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxsurfz+itrc-1)
# else
          text=vname(2,indxsurfz+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidsurfz(itrc), 'long_name',
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxsurfz+itrc-1))
          ierr=nf_put_att_text (ncid, vidsurfz(itrc), 'units', lvar,
     &                          vname(3,indxsurfz+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxsurfz+itrc-1))
          ierr=nf_put_att_text (ncid,vidsurfz(itrc), 'field',
     &                      lvar, vname(4,indxsurfz+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, vidsurfz(itrc),
     &                       indxsurfz+itrc-1, 5,
     &                        NF_FOUT, ierr)

!------------------------------
!          !
! indxsurfu
!          
          lvar=lenstr(vname(1,indxsurfu+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxsurfu+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, u2dgrd, vidsurfu(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid, vidsurfu(itrc),
     &                                   nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxsurfu+itrc-1)
# else
          text=vname(2,indxsurfu+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidsurfu(itrc), 'long_name',
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxsurfu+itrc-1))
          ierr=nf_put_att_text (ncid, vidsurfu(itrc), 'units', lvar,
     &                          vname(3,indxsurfu+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxsurfu+itrc-1))
          ierr=nf_put_att_text (ncid,vidsurfu(itrc), 'field',
     &                      lvar, vname(4,indxsurfu+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, vidsurfu(itrc),
     &                       indxsurfu+itrc-1, 5,
     &                        NF_FOUT, ierr)

!------------------------------
!          !
! indxsurfv
!          
          lvar=lenstr(vname(1,indxsurfv+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxsurfv+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, v2dgrd, vidsurfv(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid, vidsurfv(itrc),
     &                                   nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxsurfv+itrc-1)
# else
          text=vname(2,indxsurfv+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidsurfv(itrc), 'long_name',
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxsurfv+itrc-1))
          ierr=nf_put_att_text (ncid, vidsurfv(itrc), 'units', lvar,
     &                          vname(3,indxsurfv+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxsurfv+itrc-1))
          ierr=nf_put_att_text (ncid,vidsurfv(itrc), 'field',
     &                      lvar, vname(4,indxsurfv+itrc-1)(1:lvar))

        call nf_add_attribute(ncid, vidsurfv(itrc),
     &                       indxsurfv+itrc-1, 5,
     &                        NF_FOUT, ierr)

!------------------------------


        endif

!
! Leave definition mode.
! ----- ---------- -----
!
        ierr=nf_enddef(ncid)
        MPI_master_only write(stdout,'(6x,4A,1x,A,i4)') 
# ifdef AVRH
     &        'DEF_SURF_AVG - Created ',
# else      
     &        'DEF_SURF - Created ',
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
        if (ierr. eq. nf_noerr) then
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
     &           ) 'WARNING: def_surf: Actual number of records',
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
!          
! indxsurft
!
         
         lvar=lenstr(vname(1,indxsurft+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxsurft+itrc-1)(1:lvar),
     &                      vidsurft(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxsurft+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
# ifdef NC4PAR
         ierr=nf_var_par_access(ncid, vidsurft(itrc), nf_collective)
# endif

!          
! indxsurfs
!
         
         lvar=lenstr(vname(1,indxsurfs+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxsurfs+itrc-1)(1:lvar),
     &                      vidsurfs(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxsurfs+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
# ifdef NC4PAR
         ierr=nf_var_par_access(ncid, vidsurfs(itrc), nf_collective)
# endif

!          
! indxsurfz
!
         
         lvar=lenstr(vname(1,indxsurfz+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxsurfz+itrc-1)(1:lvar),
     &                      vidsurfz(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxsurfz+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
# ifdef NC4PAR
         ierr=nf_var_par_access(ncid, vidsurfz(itrc), nf_collective)
# endif
!          
! indxsurfu
!

         lvar=lenstr(vname(1,indxsurfu+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxsurfu+itrc-1)(1:lvar),
     &                      vidsurfu(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxsurfu+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
# ifdef NC4PAR
         ierr=nf_var_par_access(ncid, vidsurfu(itrc), nf_collective)
# endif
!          
! indxsurfv
!
         
         lvar=lenstr(vname(1,indxsurfv+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxsurfv+itrc-1)(1:lvar),
     &                      vidsurfv(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxsurfv+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
# ifdef NC4PAR
         ierr=nf_var_par_access(ncid, vidsurfv(itrc), nf_collective)
# endif
!
       endif

!

!
        MPI_master_only write(*,'(6x,2A,i4,1x,A,i4)') 
     &                     'def_surf: -- Opened ',
     &                     'existing file  from record =', rec
     &                      MYID
!
        MPI_master_only write(*,'(6x,2A,i4,1x,A,i4)') 
     &                     'def_surf: -- Opened ',
     &                     'existing file  from record =', rec
     &                      MYID
!
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      else
        ierr=nf_open (ncname(1:lstr), nf_write, ncid)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,'(/1x,4A,2x,A,I4/)')
     &                'def_surf: ERROR: ',
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
     &     'def_surf ERROR: Cannot ',
     &    'switch to ''nf_nofill'' more; netCDF error code =', ierr 
      endif
   1  format(/1x,'def_surf ERROR: Cannot find variable ''',
     &                   A, ''' in netCDF file ''', A, '''.'/)


  11  format(/' def_surf - unable to create diag file: ',a)
  20  format(/' def_surf - error while writing variable: ',a,
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

# undef vidsurft
# undef vidsurfs
# undef vidsurfz
# undef vidsurfu
# undef vidsurfv
!

!
# ifdef AVERAGES
#  ifndef AVRH
#   define AVRH
#   include "def_surf.F"
#   endif
#  endif /* AVERAGES */
!
#else /* OUTPUTS_SURFACE */
      subroutine def_surf_empty()
      return
      end
#endif /* OUTPUTS_SURFACE */

