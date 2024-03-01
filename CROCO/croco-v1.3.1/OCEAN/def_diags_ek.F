! $Id: def_diags_ek.F 1468 2014-02-10 12:58:35Z rblod $
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
#ifdef DIAGNOSTICS_EK

# ifndef AVRH
      subroutine def_diags_ek (ncid, total_rec, ierr)
# else
      subroutine def_diags_ek_avg(ncid, total_rec, ierr)
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
# include "diags_ek.h"
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

#  define ncname diags_ekname
#  define rec_per_file nrpfdiags_ek
#  define wrt wrtdiags_ek
#  define vidTime diags_ekTime
#  define vidTime2 diags_ekTime2
#  define vidTstep diags_ekTstep
#  define viddiags_ekHadv diags_ekHadv
#  define viddiags_ekHdiff diags_ekHdiff
#  define viddiags_ekVadv diags_ekVadv
#  define viddiags_ekCor diags_ekCor
#  define viddiags_ekPrsgrd diags_ekPrsgrd
#  define viddiags_ekHmix diags_ekHmix
#  define viddiags_ekVmix diags_ekVmix
#  define viddiags_ekrate diags_ekrate
#  define viddiags_ekvol diags_ekvol
#  define viddiags_ekVmix2 diags_ekVmix2
#  define viddiags_ekWind diags_ekWind
#  define viddiags_ekDrag diags_ekDrag
#  ifdef DIAGNOSTICS_BARO
#   define viddiags_ekBaro diags_ekBaro
#  endif
#  ifdef M3FAST
#   define viddiags_ekfast diags_ekfast
#  endif
#  ifdef DIAGNOSTICS_EK_MLD
#   define viddiags_ekHadv_mld diags_ekHadv_mld
#   define viddiags_ekHdiff_mld diags_ekHdiff_mld
#   define viddiags_ekVadv_mld diags_ekVadv_mld
#   define viddiags_ekCor_mld diags_ekCor_mld
#   define viddiags_ekPrsgrd_mld diags_ekPrsgrd_mld
#   define viddiags_ekHmix_mld diags_ekHmix_mld
#   define viddiags_ekVmix_mld diags_ekVmix_mld
#   define viddiags_ekrate_mld diags_ekrate_mld
#   define viddiags_ekvol_mld diags_ekvol_mld
#   define viddiags_ekVmix2_mld diags_ekVmix2_mld
#   ifdef DIAGNOSTICS_BARO
#    define viddiags_ekBaro_mld diags_ekBaro_mld
#   endif
#  endif

# else

#  define ncname diags_ekname_avg
#  define rec_per_file nrpfdiags_ek_avg
#  define wrt wrtdiags_ek_avg
#  define vidTime diags_ekTime_avg
#  define vidTime2 diags_ekTime2_avg
#  define vidTstep diags_ekTstep_avg
#  define viddiags_ekHadv diags_ekHadv_avg
#  define viddiags_ekHdiff diags_ekHdiff_avg
#  define viddiags_ekVadv diags_ekVadv_avg
#  define viddiags_ekCor diags_ekCor_avg
#  define viddiags_ekPrsgrd diags_ekPrsgrd_avg
#  define viddiags_ekHmix diags_ekHmix_avg
#  define viddiags_ekVmix diags_ekVmix_avg
#  define viddiags_ekrate diags_ekrate_avg
#  define viddiags_ekvol diags_ekvol_avg
#  define viddiags_ekVmix2 diags_ekVmix2_avg
#  define viddiags_ekWind diags_ekWind_avg
#  define viddiags_ekDrag diags_ekDrag_avg
#  ifdef DIAGNOSTICS_BARO
#   define viddiags_ekBaro diags_ekBaro_avg
#  endif
#  ifdef M3FAST
#   define viddiags_ekfast diags_ekfast_avg
#  endif
#  ifdef DIAGNOSTICS_EK_MLD
#   define viddiags_ekHadv_mld diags_ekHadv_mld_avg
#   define viddiags_ekHdiff_mld diags_ekHdiff_mld_avg
#   define viddiags_ekVadv_mld diags_ekVadv_mld_avg
#   define viddiags_ekCor_mld diags_ekCor_mld_avg
#   define viddiags_ekPrsgrd_mld diags_ekPrsgrd_mld_avg
#   define viddiags_ekHmix_mld diags_ekHmix_mld_avg
#   define viddiags_ekVmix_mld diags_ekVmix_mld_avg
#   define viddiags_ekrate_mld diags_ekrate_mld_avg
#   define viddiags_ekvol_mld diags_ekvol_mld_avg
#   define viddiags_ekVmix2_mld diags_ekVmix2_mld_avg
#   ifdef DIAGNOSTICS_BARO
#    define viddiags_ekBaro_mld diags_ekBaro_mld_avg
#   endif
#  endif

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
        create_new_file=ldefdiags_ek
#  else
        create_new_file=ldefdiags_ek_avg
#  endif
      else
        create_new_file=.false.
      endif
# else
#  ifndef AVRH
      create_new_file=ldefdiags_ek
#  else
      create_new_file=ldefdiags_ek_avg
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
     &       'CREATE DIAGS_EK HIS/AVG NC4 PARALLEL FILE'
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

! Energy diagnostics variables.
! -------- ----------- --------
!
          itrc=1
!
          if (wrt(itrc)) then
!
!------------------------------
!
! indxekHadv
!          
          lvar=lenstr(vname(1,indxekHadv+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxekHadv+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_ekHadv(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid, viddiags_ekHadv(itrc),
     &                                          nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxekHadv+itrc-1)
# else
          text=vname(2,indxekHadv+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid,  viddiags_ekHadv(itrc), 
     &                          'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekHadv+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekHadv(itrc), 
     &                          'units', lvar,
     &                          vname(3,indxekHadv+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekHadv+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekHadv(itrc), 'field',
     &                      lvar, vname(4,indxekHadv+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekHadv(itrc),
     &                          indxekHadv+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekHdiff
!          
          lvar=lenstr(vname(1,indxekHdiff+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxekHdiff+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_ekHdiff(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid, viddiags_ekHdiff(itrc),
     &                                          nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxekHdiff+itrc-1)
# else
          text=vname(2,indxekHdiff+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekHdiff(itrc), 
     &                          'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekHdiff+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekHdiff(itrc), 
     &                          'units', lvar,
     &                          vname(3,indxekHdiff+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekHdiff+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekHdiff(itrc), 'field',
     &                      lvar, vname(4,indxekHdiff+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekHdiff(itrc),
     &                          indxekHdiff+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekVadv
!          
          lvar=lenstr(vname(1,indxekVadv+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxekVadv+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_ekVadv(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid, viddiags_ekVadv(itrc),
     &                                          nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxekVadv+itrc-1)
# else
          text=vname(2,indxekVadv+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekVadv(itrc), 
     &                          'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekVadv+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekVadv(itrc),
     &                       'units', lvar,
     &                        vname(3,indxekVadv+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekVadv+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_ekVadv(itrc), 'field',
     &                      lvar, vname(4,indxekVadv+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekVadv(itrc),
     &                          indxekVadv+itrc-1, 5,
     &                           NF_FOUT, ierr)
!
!------------------------------
!
! indxekCor
!          
          lvar=lenstr(vname(1,indxekCor+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxekCor+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_ekCor(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid, viddiags_ekCor(itrc),
     &                                          nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxekCor+itrc-1)
# else
          text=vname(2,indxekCor+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekCor(itrc), 
     &                          'long_name', lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxekCor+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekCor(itrc),
     &                          'units', lvar,
     &                          vname(3,indxekCor+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekCor+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_ekCor(itrc), 'field',
     &                      lvar, vname(4,indxekCor+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekCor(itrc),
     &                          indxekCor+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekPrsgrd
!          
          lvar=lenstr(vname(1,indxekPrsgrd+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxekPrsgrd+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_ekPrsgrd(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid, viddiags_ekPrsgrd(itrc),
     &                                          nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxekPrsgrd+itrc-1)
# else
          text=vname(2,indxekPrsgrd+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekPrsgrd(itrc),
     &                         'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekPrsgrd+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekPrsgrd(itrc),
     &                         'units', lvar,
     &                          vname(3,indxekPrsgrd+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekPrsgrd+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_ekPrsgrd(itrc), 'field',
     &                      lvar, vname(4,indxekPrsgrd+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekPrsgrd(itrc),
     &                          indxekPrsgrd+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekHmix
!          
          lvar=lenstr(vname(1,indxekHmix+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxekHmix+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_ekHmix(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid, viddiags_ekHmix(itrc),
     &                                          nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxekHmix+itrc-1)
# else
          text=vname(2,indxekHmix+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekHmix(itrc),
     &                         'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekHmix+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekHmix(itrc),
     &                         'units', lvar,
     &                          vname(3,indxekHmix+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekHmix+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_ekHmix(itrc), 'field',
     &                      lvar, vname(4,indxekHmix+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekHmix(itrc),
     &                          indxekHmix+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekVmix
!          
          lvar=lenstr(vname(1,indxekVmix+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxekVmix+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_ekVmix(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid, viddiags_ekVmix(itrc),
     &                                          nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxekVmix+itrc-1)
# else
          text=vname(2,indxekVmix+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekVmix(itrc),
     &                         'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekVmix+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekVmix(itrc),
     &                         'units', lvar,
     &                          vname(3,indxekVmix+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekVmix+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_ekVmix(itrc), 'field',
     &                      lvar, vname(4,indxekVmix+itrc-1)(1:lvar))
          call nf_add_attribute(ncid, viddiags_ekVmix(itrc),
     &                          indxekVmix+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekrate
!          
          lvar=lenstr(vname(1,indxekrate+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxekrate+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_ekrate(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid, viddiags_ekrate(itrc),
     &                                          nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxekrate+itrc-1)
# else
          text=vname(2,indxekrate+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekrate(itrc),'long_name',
     &                          lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekrate+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekrate(itrc),
     &                          'units', lvar,
     &                          vname(3,indxekrate+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekrate+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_ekrate(itrc), 'field',
     &                      lvar, vname(4,indxekrate+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekrate(itrc),
     &                          indxekrate+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekvol
!          
          lvar=lenstr(vname(1,indxekvol+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxekvol+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_ekvol(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid, viddiags_ekvol(itrc),
     &                                         nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxekvol+itrc-1)
# else
          text=vname(2,indxekvol+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekvol(itrc), 'long_name',
     &                          lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekvol+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekvol(itrc),
     &                       'units', lvar,
     &                          vname(3,indxekvol+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekvol+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_ekvol(itrc), 'field',
     &                      lvar, vname(4,indxekvol+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekvol(itrc),
     &                          indxekvol+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekVmix2
!          
          lvar=lenstr(vname(1,indxekVmix2+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxekVmix2+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_ekVmix2(itrc))
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_ekVmix2(itrc),
     &                                          nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxekVmix2+itrc-1)
# else
          text=vname(2,indxekVmix2+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekVmix2(itrc),
     &                         'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekVmix2+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekVmix2(itrc),
     &                         'units', lvar,
     &                          vname(3,indxekVmix2+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekVmix2+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekVmix2(itrc),'field',
     &                    lvar, vname(4,indxekVmix2+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekVmix2(itrc),
     &                          indxekVmix2+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekWind
!          
          lvar=lenstr(vname(1,indxekWind+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxekWind+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_ekWind(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid, viddiags_ekWind(itrc),
     &                                          nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxekWind+itrc-1)
# else
          text=vname(2,indxekWind+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekWind(itrc),
     &                         'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekWind+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekWind(itrc),
     &                          'units', lvar,
     &                          vname(3,indxekWind+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekWind+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekWind(itrc),'field',
     &                    lvar, vname(4,indxekWind+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekWind(itrc),
     &                          indxekWind+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekDrag
!          
          lvar=lenstr(vname(1,indxekDrag+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxekDrag+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_ekDrag(itrc))
# ifdef NC4PAR
           ierr=nf_var_par_access(ncid, viddiags_ekDrag(itrc),
     &                                          nf_collective)
# endif
# ifdef AVRH
          text='averaged '/ /vname(2,indxekDrag+itrc-1)
# else
          text=vname(2,indxekDrag+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekDrag(itrc),
     &                         'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekDrag+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekDrag(itrc),
     &                          'units', lvar,
     &                          vname(3,indxekDrag+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekDrag+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekWind(itrc),'field',
     &                    lvar, vname(4,indxekDrag+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekWind(itrc),
     &                          indxekDrag+itrc-1, 5,
     &                          NF_FOUT, ierr)

# if defined DIAGNOSTICS_BARO
!------------------------------
!
! indxekBaro
!
          lvar=lenstr(vname(1,indxekBaro+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxekBaro+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_ekBaro(itrc))
#  ifdef NC4PAR
           ierr=nf_var_par_access(ncid, viddiags_ekBaro(itrc),
     &                                          nf_collective)
#  endif
#  ifdef AVRH
          text='averaged '/ /vname(2,indxekBaro+itrc-1)
#  else
          text=vname(2,indxekBaro+itrc-1)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekBaro(itrc),
     &                         'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekBaro+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekBaro(itrc),
     &                          'units', lvar,
     &                          vname(3,indxekBaro+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekBaro+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_ekBaro(itrc), 'field',
     &                      lvar, vname(4,indxekBaro+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekBaro(itrc),
     &                          indxekBaro+itrc-1, 5,
     &                          NF_FOUT, ierr)
# endif

# if defined M3FAST
!------------------------------
!
! indxekfast
!
          lvar=lenstr(vname(1,indxekfast+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxekfast+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_ekfast(itrc))
#  ifdef NC4PAR
           ierr=nf_var_par_access(ncid, viddiags_ekfast(itrc),
     &                                          nf_collective)
#  endif
#  ifdef AVRH
          text='averaged '/ /vname(2,indxekfast+itrc-1)
#  else
          text=vname(2,indxekfast+itrc-1)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekfast(itrc),
     &                         'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekfast+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekfast(itrc),
     &                         'units', lvar,
     &                          vname(3,indxekfast+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekfast+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_ekfast(itrc), 'field',
     &                      lvar, vname(4,indxekfast+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekfast(itrc),
     &                          indxekfast+itrc-1, 5,
     &                          NF_FOUT, ierr)
# endif

# ifdef DIAGNOSTICS_EK_MLD
!------------------------------
!
! indxekHadv_mld
!          
          lvar=lenstr(vname(1,indxekHadv_mld+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxekHadv_mld+itrc-1)(1:lvar),
     &                    NF_FOUT, 3, r2dgrd, viddiags_ekHadv_mld(itrc))
#  ifdef NC4PAR
           ierr=nf_var_par_access(ncid, viddiags_ekHadv_mld(itrc),
     &                                              nf_collective)
#  endif
#  ifdef AVRH
          text='averaged '/ /vname(2,indxekHadv_mld+itrc-1)
#  else
          text=vname(2,indxekHadv_mld+itrc-1)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekHadv_mld(itrc), 
     &                          'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekHadv_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekHadv_mld(itrc),
     &                          'units', lvar,
     &                          vname(3,indxekHadv_mld+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekHadv_mld+itrc-1))
          ierr=nf_put_att_text (ncid,viddiags_ekHadv_mld(itrc), 'field',
     &                     lvar, vname(4,indxekHadv_mld+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekHadv_mld(itrc),
     &                          indxekHadv_mld+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekHdiff_mld
!          
          lvar=lenstr(vname(1,indxekHdiff_mld+itrc-1))
          ierr=nf_def_var (ncid, 
     &                     vname(1,indxekHdiff_mld+itrc-1)(1:lvar),
     &                   NF_FOUT, 3, r2dgrd, viddiags_ekHdiff_mld(itrc))
#  ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_ekHdiff_mld(itrc),
     &                   i                          nf_collective)
#  endif
#  ifdef AVRH
          text='averaged '/ /vname(2,indxekHdiff_mld+itrc-1)
#  else
          text=vname(2,indxekHdiff_mld+itrc-1)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekHdiff_mld(itrc),
     &                          'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekHdiff_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekHdiff_mld(itrc),
     &                          'units', lvar,
     &                          vname(3,indxekHdiff_mld+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekHdiff_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekHdiff_mld(itrc),
     &                          'field',lvar, 
     &                          vname(4,indxekHdiff_mld+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekHdiff_mld(itrc),
     &                          indxekHdiff_mld+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekVadv_mld
!          
          lvar=lenstr(vname(1,indxekVadv_mld+itrc-1))
          ierr=nf_def_var (ncid, 
     &                     vname(1,indxekVadv_mld+itrc-1)(1:lvar),
     &                    NF_FOUT, 3, r2dgrd, viddiags_ekVadv_mld(itrc))
#  ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_ekVadv_mld(itrc),
     &                                          nf_collective)
#  endif
#  ifdef AVRH
          text='averaged '/ /vname(2,indxekVadv_mld+itrc-1)
#  else
          text=vname(2,indxekVadv_mld+itrc-1)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekVadv_mld(itrc),
     &                          'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekVadv_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekVadv_mld(itrc),
     &                          'units', lvar,
     &                          vname(3,indxekVadv_mld+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekVadv_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekVadv_mld(itrc), 
     &                          'field', lvar, 
     &                          vname(4,indxekVadv_mld+itrc-1)(1:lvar))
          call nf_add_attribute(ncid, viddiags_ekVadv_mld(itrc),
     &                          indxekVadv_mld+itrc-1, 5,
     &                          NF_FOUT, ierr)
!------------------------------
!
! indxekCor_mld
!          
          lvar=lenstr(vname(1,indxekCor_mld+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxekCor_mld+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_ekCor_mld(itrc))
#  ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_ekCor_mld(itrc),
     &                                            nf_collective)
#  endif
#  ifdef AVRH
          text='averaged '/ /vname(2,indxekCor_mld+itrc-1)
#  else
          text=vname(2,indxekCor_mld+itrc-1)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekCor_mld(itrc),
     &                          'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekCor_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekCor_mld(itrc),
     &                          'units', lvar,
     &                          vname(3,indxekCor_mld+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekCor_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekCor_mld(itrc), 'field',
     &                      lvar, vname(4,indxekCor_mld+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekCor_mld(itrc),
     &                          indxekCor_mld+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekPrsgrd_mld
!          
          lvar=lenstr(vname(1,indxekPrsgrd_mld+itrc-1))
          ierr=nf_def_var (ncid, 
     &                     vname(1,indxekPrsgrd_mld+itrc-1)(1:lvar),
     &                  NF_FOUT, 3, r2dgrd, viddiags_ekPrsgrd_mld(itrc))
#  ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_ekPrsgrd_mld(itrc),
     &                                               nf_collective)
#  endif
#  ifdef AVRH
          text='averaged '/ /vname(2,indxekPrsgrd_mld+itrc-1)
#  else
          text=vname(2,indxekPrsgrd_mld+itrc-1)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekPrsgrd_mld(itrc),
     &                          'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekPrsgrd_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekPrsgrd_mld(itrc),
     &                          'units', lvar,
     &                         vname(3,indxekPrsgrd_mld+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekPrsgrd_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekPrsgrd_mld(itrc),
     &                          'field', lvar, 
     &                         vname(4,indxekPrsgrd_mld+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekPrsgrd_mld(itrc),
     &                          indxekPrsgrd_mld+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekHmix_mld
!          
          lvar=lenstr(vname(1,indxekHmix_mld+itrc-1))
          ierr=nf_def_var (ncid, 
     &                     vname(1,indxekHmix_mld+itrc-1)(1:lvar),
     &                    NF_FOUT, 3, r2dgrd, viddiags_ekHmix_mld(itrc))
#  ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_ekHmix_mld(itrc),
     &                                             nf_collective)
#  endif
#  ifdef AVRH
          text='averaged '/ /vname(2,indxekHmix_mld+itrc-1)
#  else
          text=vname(2,indxekHmix_mld+itrc-1)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekHmix_mld(itrc),
     &                          'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekHmix_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekHmix_mld(itrc),
     &                          'units', lvar,
     &                          vname(3,indxekHmix_mld+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekHmix_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekHmix_mld(itrc),
     &                          'field', lvar,
     &                          vname(4,indxekHmix_mld+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekHmix_mld(itrc),
     &                          indxekHmix_mld+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekVmix_mld
!          
          lvar=lenstr(vname(1,indxekVmix_mld+itrc-1))
          ierr=nf_def_var (ncid, 
     &                     vname(1,indxekVmix_mld+itrc-1)(1:lvar),
     &                    NF_FOUT, 3, r2dgrd, viddiags_ekVmix_mld(itrc))
#  ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_ekVmix_mld(itrc),
     &                                             nf_collective)
#  endif
#  ifdef AVRH
          text='averaged '/ /vname(2,indxekVmix_mld+itrc-1)
#  else
          text=vname(2,indxekVmix_mld+itrc-1)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekVmix_mld(itrc),
     &                          'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekVmix_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekVmix_mld(itrc),
     &                          'units', lvar,
     &                          vname(3,indxekVmix_mld+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekVmix_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekVmix_mld(itrc),
     &                          'field', lvar, 
     &                          vname(4,indxekVmix_mld+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekVmix_mld(itrc),
     &                          indxekVmix_mld+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekrate_mld
!          
          lvar=lenstr(vname(1,indxekrate_mld+itrc-1))
          ierr=nf_def_var (ncid, 
     &                     vname(1,indxekrate_mld+itrc-1)(1:lvar),
     &                    NF_FOUT, 3, r2dgrd, viddiags_ekrate_mld(itrc))
#  ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_ekrate_mld(itrc),
     &                                             nf_collective)
#  endif
#  ifdef AVRH
          text='averaged '/ /vname(2,indxekrate_mld+itrc-1)
#  else
          text=vname(2,indxekrate_mld+itrc-1)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekrate_mld(itrc),
     &                          'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekrate_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekrate_mld(itrc),
     &                          'units', lvar,
     &                          vname(3,indxekrate_mld+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekrate_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekrate_mld(itrc),
     &                          'field', lvar,
     &                          vname(4,indxekrate_mld+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekrate_mld(itrc),
     &                          indxekrate_mld+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekvol_mld
!          
          lvar=lenstr(vname(1,indxekvol_mld+itrc-1))
          ierr=nf_def_var (ncid, 
     &                     vname(1,indxekvol_mld+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_ekvol_mld(itrc))
#  ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_ekvol_mld(itrc),
     &                                            nf_collective)
#  endif
#  ifdef AVRH
          text='averaged '/ /vname(2,indxekvol_mld+itrc-1)
#  else
          text=vname(2,indxekvol_mld+itrc-1)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekvol_mld(itrc),
     &                          'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekvol_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekvol_mld(itrc),
     &                          'units', lvar,
     &                          vname(3,indxekvol_mld+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekvol_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekvol_mld(itrc), 'field',
     &                      lvar, vname(4,indxekvol_mld+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekvol_mld(itrc),
     &                          indxekvol_mld+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
!------------------------------
!
! indxekVmix2_mld
!          
          lvar=lenstr(vname(1,indxekVmix2_mld+itrc-1))
          ierr=nf_def_var (ncid, 
     &                     vname(1,indxekVmix2_mld+itrc-1)(1:lvar),
     &                   NF_FOUT, 3, r2dgrd, viddiags_ekVmix2_mld(itrc))
#  ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_ekVmix2_mld(itrc),
     &                                              nf_collective)
#  endif
#  ifdef AVRH
          text='averaged '/ /vname(2,indxekVmix2_mld+itrc-1)
#  else
          text=vname(2,indxekVmix2_mld+itrc-1)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekVmix2_mld(itrc),
     &                          'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekVmix2_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekVmix2_mld(itrc),
     &                          'units', lvar,
     &                          vname(3,indxekVmix2_mld+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekVmix2_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekVmix2_mld(itrc),
     &                          'field', lvar,
     &                          vname(4,indxekVmix2_mld+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekVmix2_mld(itrc),
     &                          indxekVmix2_mld+itrc-1, 5,
     &                          NF_FOUT, ierr)
!
#  if defined DIAGNOSTICS_BARO
!------------------------------
!
! indxekBaro_mld
!
          lvar=lenstr(vname(1,indxekBaro_mld+itrc-1))
          ierr=nf_def_var (ncid, 
     &                     vname(1,indxekBaro_mld+itrc-1)(1:lvar),
     &                    NF_FOUT, 3, r2dgrd, viddiags_ekBaro_mld(itrc))
#   ifdef NC4PAR
          ierr=nf_var_par_access(ncid, viddiags_ekBaro_mld(itrc),
     &                                             nf_collective)
#   endif
#   ifdef AVRH
          text='averaged '/ /vname(2,indxekBaro_mld+itrc-1)
#   else
          text=vname(2,indxekBaro_mld+itrc-1)
#   endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_ekBaro_mld(itrc),
     &                         'long_name', lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxekBaro_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekBaro_mld(itrc),
     &                         'units', lvar,
     &                          vname(3,indxekBaro_mld+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxekBaro_mld+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_ekBaro_mld(itrc),
     &                          'field', lvar,
     &                          vname(4,indxekBaro_mld+itrc-1)(1:lvar))

          call nf_add_attribute(ncid, viddiags_ekBaro_mld(itrc),
     &                          indxekBaro_mld+itrc-1, 5,
     &                          NF_FOUT, ierr)
#  endif
# endif /* DIAGNOSTICS_EK_MLD */

        endif

!
! Leave definition mode.
! ----- ---------- -----
!
        ierr=nf_enddef(ncid)
        MPI_master_only write(stdout,'(6x,4A,1x,A,i4)')
# ifdef AVRH
     &        'DEF_DIAGS_EK_AVG - Created ',
# else      
     &        'DEF_DIAGS_EK - Created ',
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
     &           ) 'WARNING: def_diags_ek: Actual number of records',
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
! indxekHadv
!
         lvar=lenstr(vname(1,indxekHadv+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxekHadv+itrc-1)(1:lvar),
     &                      viddiags_ekHadv(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekHadv+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,viddiags_ekHadv(itrc),nf_collective)
# endif
!
! indxekHdiff
!
         lvar=lenstr(vname(1,indxekHdiff+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxekHdiff+itrc-1)(1:lvar),
     &                      viddiags_ekHdiff(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekHdiff+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif      
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,viddiags_ekHdiff(itrc),
     &                        nf_collective)
# endif
!
! indxekVadv
!
         lvar=lenstr(vname(1,indxekVadv+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxekVadv+itrc-1)(1:lvar),
     &                      viddiags_ekVadv(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekVadv+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif      
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,viddiags_ekVadv(itrc),
     &                        nf_collective)
# endif
!
! indxekCor
!
         lvar=lenstr(vname(1,indxekCor+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxekCor+itrc-1)(1:lvar),
     &                      viddiags_ekCor(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekCor+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,viddiags_ekCor(itrc),
     &                        nf_collective)
# endif
!
! indxekPrsgrd
!
         lvar=lenstr(vname(1,indxekPrsgrd+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxekPrsgrd+itrc-1)(1:lvar),
     &                      viddiags_ekPrsgrd(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekPrsgrd+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,viddiags_ekPrsgrd(itrc),
     &                        nf_collective)
# endif
!
! indxekHmix
!
         lvar=lenstr(vname(1,indxekHmix+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxekHmix+itrc-1)(1:lvar),
     &                      viddiags_ekHmix(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekHmix+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,viddiags_ekHmix(itrc),
     &                        nf_collective)
# endif
!
! indxekVmix
!
         lvar=lenstr(vname(1,indxekVmix+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxekVmix+itrc-1)(1:lvar),
     &                      viddiags_ekVmix(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekVmix+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
# ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekVmix(itrc),
     &                        nf_collective)
# endif
!
! indxekrate
!
         lvar=lenstr(vname(1,indxekrate+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxekrate+itrc-1)(1:lvar),
     &                      viddiags_ekrate(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekrate+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
# ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekrate(itrc),
     &                        nf_collective)
# endif
!
! indxekvol
!
         lvar=lenstr(vname(1,indxekvol+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxekvol+itrc-1)(1:lvar),
     &                      viddiags_ekvol(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekvol+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
# ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekvol(itrc),
     &                        nf_collective)
# endif

!
! indxekVmix2
!
         lvar=lenstr(vname(1,indxekVmix2+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxekVmix2+itrc-1)(1:lvar),
     &                      viddiags_ekVmix2(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekVmix2+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
# ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekVmix2(itrc),
     &                        nf_collective)
# endif

!
! indxekWind
!
         lvar=lenstr(vname(1,indxekWind+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxekWind+itrc-1)(1:lvar),
     &                      viddiags_ekWind(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekWind+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
# ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekWind(itrc),
     &                        nf_collective)
# endif

!
! indxekDrag
!
         lvar=lenstr(vname(1,indxekDrag+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxekDrag+itrc-1)(1:lvar),
     &                      viddiags_ekDrag(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekDrag+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
# ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekDrag(itrc),
     &                        nf_collective)
# endif

# if defined DIAGNOSTICS_BARO
!
! indxekBaro
!
         lvar=lenstr(vname(1,indxekBaro+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxekBaro+itrc-1)(1:lvar),
     &                      viddiags_ekBaro(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekBaro+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif
#  ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekBaro(itrc),
     &                        nf_collective)
#  endif
# endif

# if defined M3FAST
!
! indxekfast
!
         lvar=lenstr(vname(1,indxekfast+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxekfast+itrc-1)(1:lvar),
     &                      viddiags_ekfast(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekfast+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif
#  ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekfast(itrc),
     &                        nf_collective)
#  endif
# endif

# ifdef DIAGNOSTICS_EK_MLD
!          
! indxekHadv_mld
!
         lvar=lenstr(vname(1,indxekHadv_mld+itrc-1))
         ierr=nf_inq_varid (ncid, 
     &                      vname(1,indxekHadv_mld+itrc-1)(1:lvar),
     &                      viddiags_ekHadv_mld(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekHadv_mld+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
#  ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekHadv_mld(itrc),
     &                        nf_collective)
#  endif
!
! indxekHdiff_mld
!
         lvar=lenstr(vname(1,indxekHdiff_mld+itrc-1))
         ierr=nf_inq_varid (ncid, 
     &                      vname(1,indxekHdiff_mld+itrc-1)(1:lvar),
     &                      viddiags_ekHdiff_mld(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekHdiff_mld+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif      
#  ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekHdiff_mld(itrc),
     &                        nf_collective)
#  endif
!
! indxekVadv_mld
!
         lvar=lenstr(vname(1,indxekVadv_mld+itrc-1))
         ierr=nf_inq_varid (ncid, 
     &                      vname(1,indxekVadv_mld+itrc-1)(1:lvar),
     &                      viddiags_ekVadv_mld(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekVadv_mld+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif      
#  ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekVadv_mld(itrc),
     &                        nf_collective)
#  endif
!
! indxekCor_mld
!
         lvar=lenstr(vname(1,indxekCor_mld+itrc-1))
         ierr=nf_inq_varid (ncid, 
     &                      vname(1,indxekCor_mld+itrc-1)(1:lvar),
     &                      viddiags_ekCor_mld(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekCor_mld+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif
#  ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekCor_mld(itrc),
     &                        nf_collective)
#  endif
!
! indxekPrsgrd_mld
!
         lvar=lenstr(vname(1,indxekPrsgrd_mld+itrc-1))
         ierr=nf_inq_varid (ncid, 
     &                      vname(1,indxekPrsgrd_mld+itrc-1)(1:lvar),
     &                      viddiags_ekPrsgrd_mld(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekPrsgrd_mld+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif
#  ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekPrsgrd_mld(itrc),
     &                        nf_collective)
#  endif
!
! indxekHmix_mld
!
         lvar=lenstr(vname(1,indxekHmix_mld+itrc-1))
         ierr=nf_inq_varid (ncid, 
     &                      vname(1,indxekHmix_mld+itrc-1)(1:lvar),
     &                      viddiags_ekHmix_mld(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekHmix_mld+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif
#  ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekHmix_mld(itrc),
     &                        nf_collective)
#  endif
!
! indxekVmix_mld
!
         lvar=lenstr(vname(1,indxekVmix_mld+itrc-1))
         ierr=nf_inq_varid (ncid, 
     &                      vname(1,indxekVmix_mld+itrc-1)(1:lvar),
     &                      viddiags_ekVmix_mld(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekVmix_mld+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
#  ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekVmix_mld(itrc),
     &                        nf_collective)
#  endif
!
! indxekrate_mld
!
         lvar=lenstr(vname(1,indxekrate_mld+itrc-1))
         ierr=nf_inq_varid (ncid, 
     &                      vname(1,indxekrate_mld+itrc-1)(1:lvar),
     &                      viddiags_ekrate_mld(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekrate_mld+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
#  ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekrate_mld(itrc),
     &                        nf_collective)
#  endif
!
! indxekvol_mld
!
         lvar=lenstr(vname(1,indxekvol_mld+itrc-1))
         ierr=nf_inq_varid (ncid, 
     &                      vname(1,indxekvol_mld+itrc-1)(1:lvar),
     &                      viddiags_ekvol_mld(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekvol_mld+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
#  ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekvol_mld(itrc),
     &                        nf_collective)
#  endif
!
! indxekVmix2_mld
!
         lvar=lenstr(vname(1,indxekVmix2_mld+itrc-1))
         ierr=nf_inq_varid (ncid, 
     &                      vname(1,indxekVmix2_mld+itrc-1)(1:lvar),
     &                      viddiags_ekVmix2_mld(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekVmix2_mld+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
#  ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekVmix2_mld(itrc),
     &                        nf_collective)
#  endif
#  if defined DIAGNOSTICS_BARO
!
! indxekBaro_mld
!
         lvar=lenstr(vname(1,indxekBaro_mld+itrc-1))
         ierr=nf_inq_varid (ncid, 
     &                      vname(1,indxekBaro_mld+itrc-1)(1:lvar),
     &                      viddiags_ekBaro_mld(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxekBaro_mld+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif
#   ifdef NC4PAR
         ierr=nf_var_par_access(ncid,viddiags_ekBaro_mld(itrc),
     &                        nf_collective)
#   endif
#  endif
# endif
!
       endif
!
        MPI_master_only write(*,'(6x,2A,i4,1x,A,i4)') 
     &                     'def_diags_ek: -- Opened ',
     &                     'existing file  from record =', rec
     &                      MYID
!
        MPI_master_only write(*,'(6x,2A,i4,1x,A,i4)') 
     &                     'def_diags_ek: -- Opened ',
     &                     'existing file  from record =', rec
     &                      MYID
!
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      else
        ierr=nf_open (ncname(1:lstr), nf_write, ncid)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,'(/1x,4A,2x,A,I4/)')
     &                'def_diags_ek: ERROR: ',
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
     &     'def_diags_ek ERROR: Cannot ',
     &    'switch to ''nf_nofill'' more; netCDF error code =', ierr 
      endif
   1  format(/1x,'def_diags_ek ERROR: Cannot find variable ''',
     &                   A, ''' in netCDF file ''', A, '''.'/)

!
  11  format(/' def_diags_ek - unable to create diag file: ',a)
  20  format(/' def_diags_ek - error while writing variable: ',a,
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
# undef viddiags_ekHadv
# undef viddiags_ekHdiff
# undef viddiags_ekVadv
# undef viddiags_ekCor
# undef viddiags_ekPrsgrd
# undef viddiags_ekHmix
# undef viddiags_ekVmix
# undef viddiags_ekrate
# undef viddiags_ekvol
# undef viddiags_ekVmix2
# undef viddiags_ekWind
# undef viddiags_ekDrag
# if defined DIAGNOSTICS_BARO
#  undef viddiags_ekBaro
# endif
# if defined M3FAST
#  undef viddiags_ekfast
# endif
# ifdef DIAGNOSTICS_EK_MLD
#  undef viddiags_ekHadv_mld
#  undef viddiags_ekHdiff_mld
#  undef viddiags_ekVadv_mld
#  undef viddiags_ekCor_mld
#  undef viddiags_ekPrsgrd_mld
#  undef viddiags_ekHmix_mld
#  undef viddiags_ekVmix_mld
#  undef viddiags_ekrate_mld
#  undef viddiags_ekvol_mld
#  undef viddiags_ekVmix2_mld
#  if defined DIAGNOSTICS_BARO
#   undef viddiags_ekBaro_mld
#  endif
# endif
!
# ifdef AVERAGES
#  ifndef AVRH
#   define AVRH
#   include "def_diags_ek.F"
#   endif
#  endif /* AVERAGES */
!
#else /*DIAGNOSTICS_EK */
      subroutine def_diags_ek_empty()
      return
      end
#endif /* DIAGNOSTICS_EK */
