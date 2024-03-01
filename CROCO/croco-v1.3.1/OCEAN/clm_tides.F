! $Id: clm_tides.F 1476 2014-02-17 08:55:42Z rblod $
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
#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES 
      subroutine clm_tides (tile)
!
!================================================== Robert Hetland ===
!  Copyright (c) 2000 Rutgers/UCLA                                   !
!================================================ Hernan G. Arango ===
!                                                                    !
!  This routine adds tidal elevation (m) and tidal currents (m/s) to !
!  sea surface height and 2D momentum climatologies, respectively.   !
!                                                                    !
!=====================================================================
!
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
!
      trd=omp_get_thread_num()
      call clm_tides_tile (Istr,Iend,Jstr,Jend,
     &                     A2d(1,1,trd),A2d(1,2,trd),A2d(1,3,trd),
     &                     A2d(1,4,trd),A2d(1,5,trd),A2d(1,6,trd),
     &                                               A2d(1,7,trd))
      return
      end
!
!*********************************************************************
      subroutine clm_tides_tile (Istr,Iend,Jstr,Jend,Cangle,Sangle,
     &                           Cphase,Sphase,Etide,Utide,Vtide)
!*********************************************************************
!
#ifdef TIDES_MAS
      USE MOD_TIDES_MAS
#endif
      implicit none
#ifndef TIDES_MAS
# include "param.h"
#endif
# include "climat.h"
# include "grid.h"
# include "scalars.h"
# include "tides.h"
# include "boundary.h"
!
#ifdef TIDES_MAS
      CHARACTER (len=19)    :: cdate,tool_sectodat
      logical :: flag=.false.
      INTEGER               :: ian,imois,ijour,iheure,iminu,isec,jjulien
      integer :: nn,l,k,kk
      INTEGER, dimension(0:10) :: esp_corres=(/1,2,3,4,5,0,6,0,7,0,8/)
#endif
      integer
     &     Iend, Istr, Jend, Jstr, i, j, itide
      real
     &     angle, phase, omega_loc, ramp, cff
      real
     &        Cangle(PRIVATE_2D_SCRATCH_ARRAY),
     &        Cphase(PRIVATE_2D_SCRATCH_ARRAY),
     &        Sangle(PRIVATE_2D_SCRATCH_ARRAY),
     &        Sphase(PRIVATE_2D_SCRATCH_ARRAY),
     &        Etide(PRIVATE_2D_SCRATCH_ARRAY),
     &        Utide(PRIVATE_2D_SCRATCH_ARRAY),
     &        Vtide(PRIVATE_2D_SCRATCH_ARRAY)
!
# include "compute_auxiliary_bounds.h"
# ifdef MPI
#  define IRANGE Istr-1,Iend+1
#  define JRANGE Jstr-1,Jend+1
# else
#  define IRANGE IstrR,IendR
#  define JRANGE JstrR,JendR
# endif
!
      angle=0.
      phase=0.
      ramp=1.
# ifdef TIDERAMP
      ramp=TANH(dt*sec2day*float(iic-1))
# endif
      cff=2.*pi*time
!
# if defined SSH_TIDES
#  ifdef TIDES_MAS
!-------------
! MAS INIT
!-------------
      IF (FIRST_RST_TIME_STEP) THEN


! lecture des conditions aux limites
!-----------------------------------
! a faire une bonne fois pour toute
!----------------------------------
!
       IF (WESTERN_EDGE) THEN
        do j=JstrR,JendR
         num(:)=0
         rouest(:,:,j,:)=0.0
         gouest(:,:,j,:)=0.0
         nomd(:,:)=0
         fr1(:,:)=0.0
         fr2(:,:)=0.0
         do i=1,tide_nbharmcp
# ifdef TIDES_MAS_DBG
          if(j+jminmpi-1.eq.149) then 
             write(stdout,*) 'i=',nu(i),l_presence(nu(i))
          end if
# endif
          if (l_presence(nu(i))) then
           nol=ide(nu(i))
           ir0=rr0(nu(i))
           ir1=rr1(nu(i))
           ir2=rr2(nu(i))
           k=esp_corres(int(nol/100000))
           num(k)=num(k)+1
           nn=num(k)
# ifdef TIDES_MAS_DBG
           if(j+jminmpi-1.eq.149) then 
             write(stdout,*) 'ik=',nu(i),ide(nu(i)),i,k,nn
           end if
# endif
           do l=1,2 
            rouest(nn,k,j,l)=SSH_Tamp(Istr+l-2,j,i)*100.0
            gouest(nn,k,j,l)=SSH_Tphase(Istr+l-2,j,i)
            nomd(nn,k)=nol
            if(ir0 < 0) gouest(nn,k,j,l)=gouest(nn,k,j,l)+180.0
            fr1(nn,k)=ir1/float(ir0)
            fr2(nn,k)=ir2/float(ir0)
            hnouest(129,j,l)=-999.
           end do
          endif
         end do
# ifdef TIDES_MAS_DBG
        if(j+jminmpi-1.eq.149) then 
            do k=1,8
             write(stdout,*) 'k=',k,num(k)
               do kk=1,num(k)
               write(stdout,*) 'kk=',kk,rouest(kk,k,j,1)
               end do
            end do
        endif
# endif
        end do
       ENDIF
!
       IF (EASTERN_EDGE) THEN
        do j=JstrR,JendR
         num(:)=0
         rest(:,:,j,:)=0.0
         gest(:,:,j,:)=0.0
         nomd(:,:)=0
         fr1(:,:)=0.0
         fr2(:,:)=0.0
         do i=1,tide_nbharmcp
          if (l_presence(nu(i))) then
           nol=ide(nu(i))
           ir0=rr0(nu(i))
           ir1=rr1(nu(i))
           ir2=rr2(nu(i))
           k=esp_corres(int(nol/100000))
           num(k)=num(k)+1
           nn=num(k)
           do l=1,2 
            rest(nn,k,j,l)=SSH_Tamp(Iend+l-1,j,i)*100.0
            gest(nn,k,j,l)=SSH_Tphase(Iend+l-1,j,i)
            nomd(nn,k)=nol
            if(ir0 < 0) gest(nn,k,j,l)=gest(nn,k,j,l)+180.0
            fr1(nn,k)=ir1/float(ir0)
            fr2(nn,k)=ir2/float(ir0)
            hnest(129,j,l)=-999.
           end do
          endif
         end do
        end do
       ENDIF
!
       IF (NORTHERN_EDGE) THEN
        do i=IstrR,IendR
         num(:)=0
         rnord(:,:,i,:)=0.0
         gnord(:,:,i,:)=0.0
         fr1(:,:)=0.0
         fr2(:,:)=0.0
         do j=1,tide_nbharmcp
          if (l_presence(nu(j))) then
           nol=ide(nu(j))
           ir0=rr0(nu(j))
           ir1=rr1(nu(j))
           ir2=rr2(nu(j))
           k=esp_corres(int(nol/100000))
           num(k)=num(k)+1
           nn=num(k)
           do l=1,2 
            rnord(nn,k,i,l)=SSH_Tamp(i,Jend+l-1,j)*100.0
            gnord(nn,k,i,l)=SSH_Tphase(i,Jend+l-1,j)
            if(ir0 < 0) gnord(nn,k,i,l)=gnord(nn,k,i,l)+180.0
            nomd(nn,k)=nol
            fr1(nn,k)=ir1/float(ir0)
            fr2(nn,k)=ir2/float(ir0)
            hnnord(129,i,l)=-999.
           end do 
          end if
         end do
        end do
       ENDIF
!
       IF (SOUTHERN_EDGE) THEN
        do i=IstrR,IendR
         num(:)=0
         rsud(:,:,i,:)=0.0
         gsud(:,:,i,:)=0.0
         fr1(:,:)=0.0
         fr2(:,:)=0.0
         do j=1,tide_nbharmcp
          if (l_presence(nu(j))) then
           nol=ide(nu(j))
           ir0=rr0(nu(j))
           ir1=rr1(nu(j))
           ir2=rr2(nu(j))
           k=esp_corres(int(nol/100000))
           num(k)=num(k)+1
           nn=num(k)
           do l=1,2 
            rsud(nn,k,i,l)=SSH_Tamp(i,Jstr+l-2,j)*100.0
            gsud(nn,k,i,l)=SSH_Tphase(i,Jstr+l-2,j)
            if(ir0 < 0) gsud(nn,k,i,l)=gsud(nn,k,i,l)+180.0
            nomd(nn,k)=nol
            fr1(nn,k)=ir1/float(ir0)
            fr2(nn,k)=ir2/float(ir0)
            hnsud(129,i,l)=-999.
           end do 
          end if
         end do
        end do
       ENDIF
     
      ENDIF ! FIRST TIME STEP
#  endif /* MAS INIT */ 

#ifndef TIDES_MAS
!---------------------------------------------------------------------
!  Add tidal elevation (m) to sea surface height climatology.
!---------------------------------------------------------------------
!
      do j=JRANGE
        do i=IRANGE
          Etide(i,j)=0.
        enddo
      enddo
      do itide=1,Ntides
        if (Tperiod(itide).gt.0.) then
          omega_loc=cff/Tperiod(itide)
          do j=JRANGE
            do i=IRANGE
              Etide(i,j)=Etide(i,j)+ramp*
     &                   SSH_Tamp(i,j,itide)*
     &                   COS(omega_loc-SSH_Tphase(i,j,itide))
#  ifdef MASKING
              Etide(i,j)=Etide(i,j)*rmask(i,j)
#  endif
            enddo
          enddo
        endif
      enddo
#else
!--------------------------------------------
! USE MAS COMPOSITION FOR HARMONIC PREDICTION
!--------------------------------------------
      cdate = tool_sectodat(time)
      CALL tool_decompdate(cdate,ijour,imois,ian,iheure,iminu,isec)
      tempis =iheure*3600.0+iminu*60.0+isec
!
! calcul des 128 hauteurs par jour
!---------------------------------
!
      t_0=0.0
      k=0
      if(ijour /= jourmem) then

       if (SOUTHERN_EDGE) then
        do i=IstrR,IendR
         do l=1,2
           if(rmask(i,Jstr+l-2)/=0)
     &      call ma1(flag,mynode,i,l,
     &    t_0,ijour,imois,ian,hnsud(1,i,l),rsud(1,1,i,l),
     &                 gsud(1,1,i,l),nomd,fr1,fr2,num)
         end do
        end do
       end if

       if (NORTHERN_EDGE) then
        do i=IstrR,IendR
         do l=1,2
           if(rmask(i,Jend+l-1)/=0)  then 
             call ma1(flag,mynode,i,l,t_0,
     &                ijour,imois,ian,hnnord(1,i,l),rnord(1,1,i,l),
     &                gnord(1,1,i,l),nomd,fr1,fr2,num)
           end if
         end do
        end do
       end if
       
       if (WESTERN_EDGE) then
        do j=JstrR,JendR
         do l=1,2
           if(rmask(Istr+l-2,j)/=0)   then
             call ma1(flag,mynode,j,l,
     &            t_0,ijour,imois,ian,hnouest(1,j,l),rouest(1,1,j,l),
     &            gouest(1,1,j,l),nomd,fr1,fr2,num)
           endif
         end do
        end do
       end if

       if (EASTERN_EDGE) then
        do j=JstrR,JendR
         do l=1,2
           if(rmask(Iend+l-1,j)/=0)   
     &     call ma1(flag,mynode,j,l
     &         ,t_0,ijour,imois,ian,hnest(1,j,l),rest(1,1,j,l),
     &         gest(1,1,j,l),nomd,fr1,fr2,num)
         end do
        end do
       end if

      end if

      jourmem=ijour

! interpolation parmi les 128 hauteurs par jour
!----------------------------------------------

       if (SOUTHERN_EDGE) then
        do i=IstrR,IendR
         do l=1,2
           if(rmask(i,Jstr+l-2)/=0)
     &        call ma7(flag,tempis,hnsud(1,i,l),Etide(i,Jstr+l-2),k)
         end do
        end do
       end if
         
       if (NORTHERN_EDGE) then
        do i=IstrR,IendR
         do l=1,2
           if(rmask(i,Jend+l-1)/=0) then
           call ma7(flag,tempis,hnnord(1,i,l),Etide(i,Jend+l-1),k)
           endif
         end do
        end do
       end if

       if (WESTERN_EDGE) then
        do j=JstrR,JendR
         do l=1,2
           if(rmask(Istr+l-2,j)/=0) then
# ifdef TIDES_MAS_DBG
            if (j+jminmpi-1 .eq. 149 .and. l .eq. 2) then
              flag=.true.
            endif
# endif
            call ma7(flag,tempis,hnouest(1,j,l),Etide(Istr+l-2,j),k)
              flag=.false.
           endif
         end do
        end do
# ifdef TIDES_MAS_DBG
        do j=JstrR,JendR
          do i=IstrR,IendR
            if (i+iminmpi-1 .eq. 1 .and. j+jminmpi-1 .eq. 149) then
               write(stdout,*)'etide',Etide(i,j),
     &            tempis,ian,ijour,imois,
     &           iheure,iminu,isec,time,cdate,jourmem
            endif 
         enddo
        enddo
# endif
!
       end if

       if (EASTERN_EDGE) then
        do j=JstrR,JendR
         do l=1,2
           if(rmask(Iend+l-1,j)/=0)   
     &       call ma7(flag,tempis,hnest(1,j,l),Etide(Iend+l-1,j),k)
         end do
        end do
       end if
#  ifdef MASKING
       do j=JRANGE
         do i=IRANGE
           Etide(i,j)=Etide(i,j)*ramp*rmask(i,j)
         end do
       end do 
#  endif
#endif /* MAS COMPOSITION */

!----------------------------------------------------------------
!  Add sub-tidal forcing and adjust climatology to include tides.
!----------------------------------------------------------------
#  if defined ZCLIMATOLOGY
      do j=JRANGE
        do i=IRANGE
          ssh(i,j)=ssh(i,j)+Etide(i,j)
#   ifdef WET_DRY
          if (ssh(i,j) .lt. Dcrit(i,j)-h(i,j)) then
            ssh(i,j)=Dcrit(i,j)-h(i,j)
          endif
#   endif
        enddo
      enddo
#  endif /* ZCLIMATOLOGY */
!----------------------------------------------------------------
!  If appropriate, load tidal forcing into boundary arrays.
!----------------------------------------------------------------
#  ifdef Z_FRC_BRY
#   if defined OBC_WEST
      if (WESTERN_EDGE) then
        do j=JstrR,JendR
          zetabry_west(j)=zetabry_west(j)+Etide(Istr-1,j)
#ifdef TIDES_MAS_DBG
         if (j+jminmpi-1 .eq. 149) 
     &     write(102,*)'Etide',tool_sectodat(time),
     &     h(Istr-1,j),Etide(Istr-1,j),zetabry_west(j)
#endif
#    ifdef WET_DRY
          if (zetabry_west(j) .lt. Dcrit(Istr-1,j)-h(Istr-1,j)) then
            zetabry_west(j)=Dcrit(Istr-1,j)-h(Istr-1,j)
          endif
#    endif
        enddo
      endif
#   endif
#   if defined OBC_EAST
      if (EASTERN_EDGE) then
        do j=JstrR,JendR
          zetabry_east(j)=zetabry_east(j)+Etide(Iend+1,j)
#    ifdef WET_DRY
          if (zetabry_east(j) .lt. Dcrit(Iend+1,j)-h(Iend+1,j)) then
            zetabry_east(j)=Dcrit(Iend+1,j)-h(Iend+1,j)
          endif
#    endif
        enddo
      endif
#   endif
#   if defined OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do i=IstrR,IendR
          zetabry_south(i)=zetabry_south(i)+Etide(i,Jstr-1)
#    ifdef WET_DRY
          if (zetabry_south(i) .lt. Dcrit(i,Jstr-1)-h(i,Jstr-1)) then
            zetabry_south(i)=Dcrit(i,Jstr-1)-h(i,Jstr-1)
          endif
#    endif
        enddo
      endif
#   endif
#   if defined OBC_NORTH
      if (NORTHERN_EDGE) then
        do i=IstrR,IendR
          zetabry_north(i)=zetabry_north(i)+Etide(i,Jend+1)
#    ifdef WET_DRY
          if (zetabry_north(i) .lt. Dcrit(i,Jend+1)-h(i,Jend+1)) then
            zetabry_north(i)=Dcrit(i,Jend+1)-h(i,Jend+1)
          endif
#    endif
        enddo
      endif
#   endif
#  endif /* Z_FRC_BRY */
# endif /* SSH_TIDES */
!
# if defined UV_TIDES
!
!---------------------------------------------------------------------
!  Add tidal currents (m/s) to 2D momentum climatologies.
!---------------------------------------------------------------------
!
      do j=JstrR,JendR
        do i=Istr,IendR
          Utide(i,j)=0.
        enddo
      enddo
      do j=Jstr,JendR
        do i=IstrR,IendR
          Vtide(i,j)=0.
        enddo
      enddo
      do itide=1,Ntides
        if (Tperiod(itide).gt.0.) then
          omega_loc=cff/Tperiod(itide)
          do j=Jstr-1,JendR
            do i=Istr-1,IendR
              angle=UV_Tangle(i,j,itide)
#  ifdef CURVGRID
     &                      -angler(i,j)
#  endif
              phase=omega_loc-UV_Tphase(i,j,itide)
              Cangle(i,j)=COS(angle)
              Cphase(i,j)=COS(phase)
              Sangle(i,j)=SIN(angle)
              Sphase(i,j)=SIN(phase)
            enddo
          enddo
          do j=JstrR,JendR
            do i=Istr,IendR
              Utide(i,j)=Utide(i,j)+ramp*
     &                        0.125*((UV_Tmajor(i-1,j,itide)+
     &                                UV_Tmajor(i  ,j,itide))*
     &                               (Cangle(i-1,j)+Cangle(i,j))*
     &                               (Cphase(i-1,j)+Cphase(i,j))-
     &                               (UV_Tminor(i-1,j,itide)+
     &                                UV_Tminor(i  ,j,itide))*
     &                               (Sangle(i-1,j)+Sangle(i,j))*
     &                               (Sphase(i-1,j)+Sphase(i,j)))
#  ifdef MASKING
             Utide(i,j)=Utide(i,j)*umask(i,j)
#  endif
            enddo
          enddo
          do j=Jstr,JendR
            do i=IstrR,IendR
              Vtide(i,j)=Vtide(i,j)+ramp*
     &                        0.125*((UV_Tmajor(i,j-1,itide)+
     &                                UV_Tmajor(i,j  ,itide))*
     &                               (Sangle(i,j-1)+Sangle(i,j))*
     &                               (Cphase(i,j-1)+Cphase(i,j))+
     &                               (UV_Tminor(i,j-1,itide)+
     &                                UV_Tminor(i,j  ,itide))*
     &                               (Cangle(i,j-1)+Cangle(i,j))*
     &                               (Sphase(i,j-1)+Sphase(i,j)))
#  ifdef MASKING
              Vtide(i,j)=Vtide(i,j)*vmask(i,j)
#  endif
            enddo
          enddo
        endif
      enddo
#  if defined M2CLIMATOLOGY
!
!  Add sub-tidal forcing and adjust climatology to include tides.
!
      do j=JstrR,JendR
        do i=Istr,IendR
          ubclm(i,j)=ubclm(i,j)+Utide(i,j)
        enddo
      enddo
      do j=Jstr,JendR
        do i=IstrR,IendR
          vbclm(i,j)=vbclm(i,j)+Vtide(i,j)
        enddo
      enddo
#  endif
!
!  If appropriate, load tidal forcing into boundary arrays.
!
#  ifdef M2_FRC_BRY
#   ifdef OBC_WEST
      if (WESTERN_EDGE) THEN
        do j=JstrR,JendR
          ubarbry_west(j)=ubarbry_west(j)+Utide(Istr,j)
        enddo
        do j=Jstr,JendR
          vbarbry_west(j)=vbarbry_west(j)+Vtide(Istr-1,j)
        enddo
      endif
#   endif
#   ifdef OBC_EAST
      if (EASTERN_EDGE) THEN
        do j=JstrR,JendR
          ubarbry_east(j)=ubarbry_east(j)+Utide(Iend+1,j)
        enddo
        do j=Jstr,JendR
          vbarbry_east(j)=vbarbry_east(j)+Vtide(Iend+1,j)
        enddo
      endif
#   endif
#   ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) THEN
        do i=Istr,IendR
          ubarbry_south(i)=ubarbry_south(i)+Utide(i,Jstr-1)
        enddo
        do i=IstrR,IendR
          vbarbry_south(i)=vbarbry_south(i)+Vtide(i,Jstr)
        enddo
      endif
#   endif
#   ifdef OBC_NORTH
      if (NORTHERN_EDGE) THEN
        do i=Istr,IendR
          ubarbry_north(i)=ubarbry_north(i)+Utide(i,Jend+1)
        enddo
        do i=IstrR,IendR
          vbarbry_north(i)=vbarbry_north(i)+Vtide(i,Jend+1)
        enddo
      endif
#   endif
# endif /* M2_FRC_BRY */
# endif /* UV_TIDES */

# if defined POT_TIDES
!
!-----------------------------------------------------------------------
!  Compute tidal potential (m) 
!  -- > to be applied as a pressure-gradient force
!       in prsgrd (if SOLVE3D) or step2d (if not SOLVE3D)
!-----------------------------------------------------------------------
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          Ptide(i,j)=0.
        enddo
      enddo
      do itide=1,Ntides
        if (Tperiod(itide).gt.0.) then
          omega_loc=cff/Tperiod(itide)
          do j=JstrR,JendR
            do i=IstrR,IendR
              Ptide(i,j)=Ptide(i,j)+ 
     &                   ramp*POT_Tamp(i,j,itide)*
     &                   COS(omega_loc-POT_Tphase(i,j,itide))
#  ifdef MASKING
              Ptide(i,j)=Ptide(i,j)*rmask(i,j)
#  endif
            enddo
          enddo
        endif
      enddo

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        Ptide(START_2D_ARRAY))
#  endif

# endif /* POT_TIDES */

#else
      subroutine clm_tides_empty
#endif /* SSH_TIDES || UV_TIDES */
      return
      end
