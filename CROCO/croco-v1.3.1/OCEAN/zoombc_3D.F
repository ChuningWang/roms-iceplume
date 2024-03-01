! $Id: zoombc_3D.F 1458 2014-02-03 15:01:25Z gcambon $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! This routine belongs to the specific CROCO package.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"

#if defined AGRIF && defined SOLVE3D

!====================================================================
!               subroutine Agrif_u3dbc_interp_tile
!====================================================================
!
      subroutine u3dbc_interp_tile(Istr,Iend,Jstr,Jend)

      use AGRIF_Util
      implicit none
# include "param.h"
# include "boundary.h"  
# include "climat.h"
# include "grid.h"
# include "scalars.h"
# include "zoom.h"
# include "ocean3d.h"
      integer Istr,Iend,Jstr,Jend, i,j,k
      real     :: tinterp
      INTEGER  :: nbstep3dparent
      External :: u3dinterp
      real     :: ainterp(6)
      integer  :: irhot
      real     :: rrhot
      integer  :: nsource, ndest
      integer  :: parentnbstep
!
# include "compute_auxiliary_bounds.h"
!
# ifdef MPI
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif
!       return
      irhot=Agrif_Irhot()
      rrhot=real(irhot)

C$OMP BARRIER
C$OMP MASTER

      parentnbstep=Agrif_Parent_Nb_Step()

       IF ((((nbcoarse == 1)    .AND.(nnew == 3)).OR.
     &      ((nbcoarse == irhot).AND.(nnew /= 3))).AND.
     &      (UTimeindex.NE.parentnbstep)) THEN

        tinterp=1.
# ifdef MASKING
        Agrif_UseSpecialvalue=.true.
# endif        
        Agrif_Specialvalue=0.
        Call Agrif_Set_bc(uid,(/0,0/),InterpolationShouldbemade=.TRUE.)
        Call Agrif_Bc_variable(uid,calledweight=tinterp,
     &                         procname=u3dinterp)
        Agrif_UseSpecialValue=.false.
        UTimeindex = parentnbstep

      endif
C$OMP END MASTER
C$OMP BARRIER
      
       nbstep3dparent=Agrif_Parent(nbstep3d)
                    
       tinterp = real(nbcoarse)/rrhot
       IF (nnew == 3) tinterp = tinterp - 1./(2.*rrhot)

       ainterp = 0.
       
       IF (nnew == 3) then  
         ainterp(1) = 0.5-tinterp
         ainterp(2) = 0.5+tinterp
       else 
         ainterp(3) = -tinterp
         ainterp(4) = 1.+tinterp
       endif
       
       IF (nbstep3dparent .LE. 1) THEN
       ainterp = 0.
       ainterp(2) = 1.
       ENDIF
       
       if (nbstep3d .LE. (Agrif_irhot()-1)) then
        ainterp = 0.
        ainterp(2) = 1.
       endif

       IF ((nnew /= 3).AND.(nbcoarse == irhot)) THEN
         ainterp    = 0.
         ainterp(4) = 1.
       ENDIF
              
# ifdef AGRIF_OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do k=1,N      
          do i=Istr,IendR
#  if defined M3_FRC_BRY     
            ubry_south(i,k)=
#  else
            uclm(i,Jstr-1,k)=
#  endif
     &  (ainterp(1)*U_south(i,Jstr-1,k,1)+
     &   ainterp(2)*U_south(i,Jstr-1,k,2)+
     &   ainterp(3)*U_south(i,Jstr-1,k,3)+
     &   ainterp(4)*U_south(i,Jstr-1,k,4)+
     &   ainterp(5)*u(i,Jstr-1,k,3-nstp)+
     &   ainterp(6)*u(i,Jstr-1,k,nstp))
#  ifdef MASKING     
     &    *umask(i,Jstr-1)
#  endif
          enddo
        enddo
      endif
# endif

# ifdef AGRIF_OBC_NORTH
      if (NORTHERN_EDGE) then      
        do k=1,N      
          do i=Istr,IendR
#  if defined M3_FRC_BRY     
            ubry_north(i,k)=
#  else   
            uclm(i,Jend+1,k)=
#  endif
     &  (ainterp(1)*U_north(i,Jend+1,k,1)+
     &   ainterp(2)*U_north(i,Jend+1,k,2)+
     &   ainterp(3)*U_north(i,Jend+1,k,3)+
     &   ainterp(4)*U_north(i,Jend+1,k,4)+
     &   ainterp(5)*u(i,Jend+1,k,3-nstp)+
     &   ainterp(6)*u(i,Jend+1,k,nstp))     
#  ifdef MASKING     
     &    *umask(i,Jend+1)
#  endif
          enddo
        enddo
      endif      
# endif
 
# ifdef AGRIF_OBC_WEST
      if (WESTERN_EDGE) then      
        do k=1,N      
          do j=JstrR,JendR
#  if defined M3_FRC_BRY     
            ubry_west(j,k)=
#  else     
            uclm(Istr,j,k)=
#  endif
     &  (ainterp(1)*U_west(Istr,j,k,1)+
     &   ainterp(2)*U_west(Istr,j,k,2)+
     &   ainterp(3)*U_west(Istr,j,k,3)+
     &   ainterp(4)*U_west(Istr,j,k,4)+
     &   ainterp(5)*u(Istr,j,k,3-nstp)+
     &   ainterp(6)*u(Istr,j,k,nstp))
#  ifdef MASKING     
     &    *umask(Istr,j)
#  endif
          enddo
        enddo
      endif      
# endif
      
# ifdef AGRIF_OBC_EAST
      if (EASTERN_EDGE) then      
        do k=1,N      
          do j=JstrR,JendR
#  if defined M3_FRC_BRY  
            ubry_east(j,k)=
#  else
            uclm(Iend+1,j,k)=
#  endif
     &  (ainterp(1)*U_east(Iend+1,j,k,1)+
     &   ainterp(2)*U_east(Iend+1,j,k,2)+
     &   ainterp(3)*U_east(Iend+1,j,k,3)+
     &   ainterp(4)*U_east(Iend+1,j,k,4)+
     &   ainterp(5)*u(Iend+1,j,k,3-nstp)+
     &   ainterp(6)*u(Iend+1,j,k,nstp))         
#  ifdef MASKING     
     &    *umask(Iend+1,j)
#  endif
          enddo
        enddo
      endif
# endif
    
      return
      end
!
!====================================================================

      subroutine u3Dinterp(tabres,i1,i2,j1,j2,k1,k2,before,nb,ndir)

      implicit none
# include "param.h"
# include "ocean3d.h"
# include "scalars.h"
# include "zoom.h"
!
# ifdef MPI
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif
!
      integer i1,i2,j1,j2,k1,k2
      real    :: tabres(i1:i2,j1:j2,k1:k2)
      logical :: before
      integer :: nb,ndir
      integer :: nparent, nsource, ndest
      logical :: western_side, eastern_side
      logical :: northern_side,southern_side

      if (before) then

        IF ((iif == 0)) THEN
          nparent = 3
        ELSE
          nparent = nnew
        ENDIF

        tabres(i1:i2,j1:j2,k1:k2) = u(i1:i2,j1:j2,k1:k2,nparent)

      else

         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)

         IF (nnew == 3) THEN
           nsource = 2
         ELSE
           nsource = 4
         ENDIF
         ndest = nsource - 1

# ifdef AGRIF_OBC_SOUTH
        if (southern_side) then
#  ifdef MPI
            if (.not.SOUTH_INTER) then
#  endif
                U_south(i1:i2,j1:j2,1:N,ndest)=
     &          U_south(i1:i2,j1:j2,1:N,nsource)
                U_south(i1:i2,j1:j2,1:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif
# ifdef AGRIF_OBC_NORTH
        if (northern_side) then
#  ifdef MPI
            if (.not.NORTH_INTER) then
#  endif
                U_north(i1:i2,j1:j2,1:N,ndest)=
     &          U_north(i1:i2,j1:j2,1:N,nsource)
                U_north(i1:i2,j1:j2,1:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif
# ifdef AGRIF_OBC_WEST
        if (western_side) then
#  ifdef MPI
            if (.not.WEST_INTER) then
#  endif
                U_west(i1:i2,j1:j2,1:N,ndest)=
     &          U_west(i1:i2,j1:j2,1:N,nsource)
                U_west(i1:i2,j1:j2,1:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif
# ifdef AGRIF_OBC_EAST
        if (eastern_side) then
#  ifdef MPI
            if (.not.EAST_INTER) then
#  endif
                U_east(i1:i2,j1:j2,1:N,ndest)=
     &          U_east(i1:i2,j1:j2,1:N,nsource)
                U_east(i1:i2,j1:j2,1:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif

      endif

      return
      end
      
!====================================================================
!                   subroutine Agrif_v3dbc_interp_tile
!====================================================================
!
      subroutine v3dbc_interp_tile(Istr,Iend,Jstr,Jend)

      use AGRIF_Util
      implicit none
# include "param.h"
# include "boundary.h"  
# include "climat.h"
# include "grid.h"
# include "scalars.h"
# include "zoom.h"
# include "ocean3d.h"
      integer Istr,Iend,Jstr,Jend, i,j,k
      real     :: tinterp
      INTEGER  :: nbstep3dparent
      external :: v3dinterp
      real     :: ainterp(6)
      integer  :: irhot
      real     :: rrhot
      integer  :: nsource, ndest
      integer  :: parentnbstep
!
# include "compute_auxiliary_bounds.h"
!
# ifdef MPI
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif 

      irhot=Agrif_Irhot()
      rrhot=real(irhot) 

C$OMP BARRIER
C$OMP MASTER
      
      parentnbstep=Agrif_Parent_Nb_Step()
      
       IF ((((nbcoarse == 1)    .AND.(nnew == 3)).OR.
     &      ((nbcoarse == irhot).AND.(nnew /= 3))).AND.
     &      (VTimeindex.NE.parentnbstep)) THEN 

        tinterp=1.
# ifdef MASKING
        AGRIF_UseSpecialvalue=.true.
# endif        
        AGRIF_Specialvalue=0.
        Call Agrif_Set_bc(vid,(/0,0/),InterpolationShouldbemade=.TRUE.)
        Call Agrif_Bc_variable(vid,calledweight=tinterp,
     &                         procname=v3dinterp)
        AGRIF_UseSpecialValue=.false.
        VTimeindex=parentnbstep

      endif
C$OMP END MASTER
C$OMP BARRIER

       nbstep3dparent=Agrif_Parent(nbstep3d)
               
       tinterp = real(nbcoarse)/rrhot
       IF (nnew == 3) tinterp = tinterp - 1./(2.*rrhot)

       ainterp = 0.
       
       IF (nnew == 3) then  
         ainterp(1) = 0.5-tinterp
         ainterp(2) = 0.5+tinterp
       else 
         ainterp(3) = -tinterp
         ainterp(4) = 1.+tinterp
       endif

       IF (nbstep3dparent .LE. 1) THEN
         ainterp = 0.
         ainterp(2) = 1.
       ENDIF

       if (nbstep3d .LE. (Agrif_irhot()-1)) then
        ainterp = 0.
        ainterp(2) = 1.
       endif
       
       IF ((nnew /= 3).AND.(nbcoarse == irhot)) THEN
         ainterp    = 0.
         ainterp(4) = 1.
       ENDIF
       
# ifdef AGRIF_OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do k=1,N      
          do i=IstrR,IendR
#  if defined M3_FRC_BRY   
            vbry_south(i,k)=
#  else
            vclm(i,Jstr,k)=
#  endif
     &    (ainterp(1)*V_south(i,Jstr,k,1)+
     &     ainterp(2)*V_south(i,Jstr,k,2)+
     &     ainterp(3)*V_south(i,Jstr,k,3)+
     &     ainterp(4)*V_south(i,Jstr,k,4)+
     &     ainterp(5)*v(i,Jstr,k,3-nstp)+
     &     ainterp(6)*v(i,Jstr,k,nstp))
#  ifdef MASKING     
     &    *vmask(i,Jstr)
#  endif
          enddo
        enddo
      endif
# endif
# ifdef AGRIF_OBC_NORTH
      if (NORTHERN_EDGE) then      
        do k=1,N      
          do i=IstrR,IendR
#  if defined M3_FRC_BRY
            vbry_north(i,k)=
#  else
            vclm(i,Jend+1,k)=
#  endif
     &    (ainterp(1)*V_north(i,Jend+1,k,1)+
     &     ainterp(2)*V_north(i,Jend+1,k,2)+
     &     ainterp(3)*V_north(i,Jend+1,k,3)+
     &     ainterp(4)*V_north(i,Jend+1,k,4)+
     &     ainterp(5)*v(i,Jend+1,k,3-nstp)+
     &     ainterp(6)*v(i,Jend+1,k,nstp))     
#  ifdef MASKING     
     &    *vmask(i,Jend+1)
#  endif 
          enddo
        enddo
      endif      
# endif
# ifdef AGRIF_OBC_WEST
      if (WESTERN_EDGE) then      
        do k=1,N      
          do j=Jstr,JendR
#  if defined M3_FRC_BRY
            vbry_west(j,k)=
#  else
            vclm(Istr-1,j,k)=
#  endif
     &    (ainterp(1)*V_west(Istr-1,j,k,1)+
     &     ainterp(2)*V_west(Istr-1,j,k,2)+
     &     ainterp(3)*V_west(Istr-1,j,k,3)+
     &     ainterp(4)*V_west(Istr-1,j,k,4)+
     &     ainterp(5)*v(Istr-1,j,k,3-nstp)+
     &     ainterp(6)*v(Istr-1,j,k,nstp))     
#  ifdef MASKING     
     &    *vmask(Istr-1,j)
#  endif
          enddo
        enddo
      endif      
# endif
# ifdef AGRIF_OBC_EAST
      if (EASTERN_EDGE) then      
        do k=1,N      
          do j=Jstr,JendR
#  if defined M3_FRC_BRY  
            vbry_east(j,k)=
#  else
            vclm(Iend+1,j,k)=
#  endif
     &    (ainterp(1)*V_east(Iend+1,j,k,1)+
     &     ainterp(2)*V_east(Iend+1,j,k,2)+
     &     ainterp(3)*V_east(Iend+1,j,k,3)+
     &     ainterp(4)*V_east(Iend+1,j,k,4)+
     &     ainterp(5)*v(Iend+1,j,k,3-nstp)+
     &     ainterp(6)*v(Iend+1,j,k,nstp))      
#  ifdef MASKING     
     &    *vmask(Iend+1,j)
#  endif
          enddo
        enddo
      endif
# endif    

      return
      end

!====================================================================

      subroutine v3Dinterp(tabres,i1,i2,j1,j2,k1,k2,before,nb,ndir)

      implicit none
# include "param.h"
# include "ocean3d.h"
# include "scalars.h"
# include "zoom.h"
!
# ifdef MPI
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif
!
       integer i1,i2,j1,j2,k1,k2
       real    :: tabres(i1:i2,j1:j2,k1:k2)
       logical :: before
       integer :: nb,ndir
       integer :: nparent, nsource, ndest
       logical :: western_side, eastern_side
       logical :: northern_side,southern_side

       if (before) then

         IF ((iif == 0)) THEN
           nparent = 3
         ELSE
           nparent = nnew
         ENDIF

         tabres(i1:i2,j1:j2,k1:k2) = v(i1:i2,j1:j2,k1:k2,nparent)

      else

         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)

         IF (nnew == 3) THEN
           nsource = 2
         ELSE
           nsource = 4
         ENDIF
         ndest = nsource - 1

# ifdef AGRIF_OBC_SOUTH
        if (southern_side) then
#  ifdef MPI
            if (.not.SOUTH_INTER) then
#  endif
                V_south(i1:i2,j1:j2,1:N,ndest)=
     &          V_south(i1:i2,j1:j2,1:N,nsource)
                V_south(i1:i2,j1:j2,1:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif
# ifdef AGRIF_OBC_NORTH
        if (northern_side) then
#  ifdef MPI
            if (.not.NORTH_INTER) then
#  endif
                V_north(i1:i2,j1:j2,1:N,ndest)=
     &          V_north(i1:i2,j1:j2,1:N,nsource)
                V_north(i1:i2,j1:j2,1:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif
# ifdef AGRIF_OBC_WEST
        if (western_side) then
#  ifdef MPI
            if (.not.WEST_INTER) then
#  endif
                V_west(i1:i2,j1:j2,1:N,ndest)=
     &          V_west(i1:i2,j1:j2,1:N,nsource)
                V_west(i1:i2,j1:j2,1:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif
# ifdef AGRIF_OBC_EAST
        if (eastern_side) then
#  ifdef MPI
            if (.not.EAST_INTER) then
#  endif
                V_east(i1:i2,j1:j2,1:N,ndest)=
     &          V_east(i1:i2,j1:j2,1:N,nsource)
                V_east(i1:i2,j1:j2,1:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif

      endif

      return
      end
!
!====================================================================
!                   subroutine Agrif_t3dbc_interp_tile
!====================================================================
!
      subroutine t3dbc_interp_tile(Istr,Iend,Jstr,Jend,indx,itrc)

      use AGRIF_Util
      implicit none
# include "param.h"
# include "boundary.h"  
# include "climat.h"
# include "grid.h"
# include "scalars.h"
# include "zoom.h"
# include "ocean3d.h"
      integer Istr,Iend,Jstr,Jend,indx,itrc,i,j,k
      real     :: tinterp
      INTEGER  :: nbstep3dparent
      external :: t3dinterp
      real     :: ainterp(4)
      integer  :: irhot
      real     :: rrhot
      integer  :: nsource, ndest
      integer  :: parentnbstep
!
# include "compute_auxiliary_bounds.h"
!
# ifdef MPI
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif   

      irhot = Agrif_Irhot()
      rrhot = real(irhot) 

       nbstep3dparent=Agrif_Parent(nbstep3d)

C$OMP BARRIER
C$OMP MASTER

      parentnbstep=Agrif_Parent_Nb_Step()

       IF ((((nbcoarse == 1)    .AND.(nnew == 3)).OR.
     &      ((nbcoarse == irhot).AND.(nnew /= 3))).AND.
     &      (TTimeindex.NE.parentnbstep)) THEN 

        tinterp=1.
# ifdef MASKING        
        Agrif_UseSpecialvalue=.true.
# endif        
        Agrif_Specialvalue=0.
        Call Agrif_Set_bc(tid,(/-1,0/),
     &     InterpolationShouldbemade=.TRUE.)
        Call Agrif_Bc_variable(tid,calledweight=tinterp,
     &                         procname = t3dinterp)
      
        Agrif_UseSpecialValue=.false.

      endif
C$OMP END MASTER
C$OMP BARRIER
                    
       tinterp = real(nbcoarse)/rrhot
       IF (nnew == 3) tinterp = tinterp - 1./(2.*rrhot)

       ainterp = 0.
       
       IF (nnew == 3) then  
         ainterp(1) = 0.5-tinterp
         ainterp(2) = 0.5+tinterp
       else 
         ainterp(3) = -tinterp
         ainterp(4) = 1.+tinterp
       endif
       
       IF (nbstep3dparent .LE. 1) THEN
       ainterp = 0.
       ainterp(2) = 1.
       ENDIF

       if (nbstep3d .LE. (Agrif_irhot()-1)) then
        ainterp = 0.
        ainterp(2) = 1.
       endif
                           
       IF ((nnew /= 3).AND.(nbcoarse == irhot)) THEN
         ainterp    = 0.
         ainterp(4) = 1.
       ENDIF
!      
# ifdef AGRIF_OBC_SOUTH
       if (SOUTHERN_EDGE) then
         do k=1,N      
           do i=IstrR,IendR
#  ifdef T_FRC_BRY
            tbry_south(i,k,itrc)=
#  else
            tclm(i,Jstr-1,k,itrc)=
#  endif
     &     (ainterp(1)*T_south(i,Jstr-1,k,1,itrc)+
     &      ainterp(2)*T_south(i,Jstr-1,k,2,itrc)+
     &      ainterp(3)*T_south(i,Jstr-1,k,3,itrc)+
     &      ainterp(4)*T_south(i,Jstr-1,k,4,itrc))
#  ifdef MASKING     
     &        *rmask(i,Jstr-1)
#  endif
          enddo
        enddo
      endif
# endif
# ifdef AGRIF_OBC_NORTH
      if (NORTHERN_EDGE) then      
        do k=1,N      
          do i=IstrR,IendR
#  ifdef T_FRC_BRY
            tbry_north(i,k,itrc)=
#  else
            tclm(i,Jend+1,k,itrc)=
#  endif
     &     (ainterp(1)*T_north(i,Jend+1,k,1,itrc)+
     &      ainterp(2)*T_north(i,Jend+1,k,2,itrc)+
     &      ainterp(3)*T_north(i,Jend+1,k,3,itrc)+
     &      ainterp(4)*T_north(i,Jend+1,k,4,itrc))     
#  ifdef MASKING     
     &          *rmask(i,Jend+1)
#  endif
          enddo
        enddo
      endif      
# endif
# ifdef AGRIF_OBC_WEST
      if (WESTERN_EDGE) then      
        do k=1,N      
          do j=JstrR,JendR
#  ifdef T_FRC_BRY
            tbry_west(j,k,itrc)=
#  else
            tclm(Istr-1,j,k,itrc)=
#  endif
     &     (ainterp(1)*T_west(Istr-1,j,k,1,itrc)+
     &      ainterp(2)*T_west(Istr-1,j,k,2,itrc)+
     &      ainterp(3)*T_west(Istr-1,j,k,3,itrc)+
     &      ainterp(4)*T_west(Istr-1,j,k,4,itrc))      
#  ifdef MASKING     
     &         *rmask(Istr-1,j)
#  endif
          enddo
        enddo
      endif      
# endif
# ifdef AGRIF_OBC_EAST
      if (EASTERN_EDGE) then      
        do k=1,N      
          do j=JstrR,JendR
#  ifdef T_FRC_BRY
            tbry_east(j,k,itrc)= 
#  else
            tclm(Iend+1,j,k,itrc)=
#  endif
     &     (ainterp(1)*T_east(Iend+1,j,k,1,itrc)+
     &      ainterp(2)*T_east(Iend+1,j,k,2,itrc)+
     &      ainterp(3)*T_east(Iend+1,j,k,3,itrc)+
     &      ainterp(4)*T_east(Iend+1,j,k,4,itrc))       
#  ifdef MASKING     
     &        *rmask(Iend+1,j)
#  endif
           enddo
        enddo
      endif
# endif     
      return
      end

!======================================================================

      subroutine t3Dinterp(tabres,i1,i2,j1,j2,k1,k2,m1,m2,before,
     &                     nb,ndir)

      implicit none
# include "param.h"
# include "ocean3d.h"
# include "scalars.h"
# include "zoom.h"
# ifdef MPI
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif

       integer i1,i2,j1,j2,k1,k2,m1,m2
       real tabres(i1:i2,j1:j2,k1:k2,m1:m2)
       logical before
       integer nb,ndir
       integer nparent
       integer nsource, ndest
       logical :: western_side, eastern_side
       logical :: northern_side,southern_side

       if (before) then
         IF ((iif == 0)) THEN
           nparent = 3
         ELSE
           nparent = nnew
         ENDIF

         tabres(i1:i2,j1:j2,k1:k2,m1:m2) =
     &               t(i1:i2,j1:j2,k1:k2,nparent,m1:m2)

       else

          western_side  = (nb == 1).AND.(ndir == 1)
          eastern_side  = (nb == 1).AND.(ndir == 2)
          southern_side = (nb == 2).AND.(ndir == 1)
          northern_side = (nb == 2).AND.(ndir == 2)

          IF (nnew == 3) THEN
            nsource = 2
          ELSE
            nsource = 4
          ENDIF
          ndest = nsource - 1

# ifdef AGRIF_OBC_SOUTH
        if (southern_side) then
#  ifdef MPI
              if (.not.SOUTH_INTER) Then
#  endif
                T_south(i1:i2,j1:j2,1:N,ndest,1:NT)=
     &          T_south(i1:i2,j1:j2,1:N,nsource,1:NT)
                T_south(i1:i2,j1:j2,1:N,nsource,1:NT) = tabres
#  ifdef MPI
              endif
#  endif
        endif
# endif
# ifdef AGRIF_OBC_NORTH
        if (northern_side) then
#  ifdef MPI
              if (.not.NORTH_INTER) then
#  endif
                T_north(i1:i2,j1:j2,1:N,ndest,1:NT)=
     &          T_north(i1:i2,j1:j2,1:N,nsource,1:NT)
                T_north(i1:i2,j1:j2,1:N,nsource,1:NT) = tabres
#  ifdef MPI
              endif
#  endif
        endif
# endif
# ifdef AGRIF_OBC_WEST
        if (western_side) then
#  ifdef MPI
              if (.not.WEST_INTER) then
#  endif
                  T_west(i1:i2,j1:j2,1:N,ndest,1:NT)=
     &            T_west(i1:i2,j1:j2,1:N,nsource,1:NT)
                  T_west(i1:i2,j1:j2,1:N,nsource,1:NT) = tabres
#  ifdef MPI
              endif
#  endif
        endif
# endif
# ifdef AGRIF_OBC_EAST
        if (eastern_side) then
#  ifdef MPI
              if (.not.EAST_INTER) then
#  endif
                  T_east(i1:i2,j1:j2,1:N,ndest,1:NT)=
     &            T_east(i1:i2,j1:j2,1:N,nsource,1:NT)
                  T_east(i1:i2,j1:j2,1:N,nsource,1:NT) = tabres
#  ifdef MPI
              endif
#  endif
        endif
# endif

        endif
      return
      end
!
!====================================================================
!               subroutine Agrif_w3dbc_interp_tile
!====================================================================
!
# ifdef NBQ
      subroutine w3dbc_interp_tile(Istr,Iend,Jstr,Jend)

      use AGRIF_Util
      implicit none
#  include "param.h"
#  include "boundary.h"  
#  include "climat.h"
#  include "grid.h"
#  include "scalars.h"
#  include "zoom.h"
#  include "ocean3d.h"
      integer Istr,Iend,Jstr,Jend, i,j,k
      real     :: tinterp
      integer  :: nbstep3dparent
      External :: w3dinterp
      real     :: ainterp(6)
      integer  :: irhot,irhox,irhoy
      real     :: rrhot,cffx,cffy
      integer  :: nsource, ndest
      integer  :: parentnbstep
      
#  include "compute_auxiliary_bounds.h"
!
#  ifdef MPI
#   define LOCALLM Lmmpi
#   define LOCALMM Mmmpi
#  else
#   define LOCALLM Lm
#   define LOCALMM Mm
#  endif
!
      irhot=Agrif_Irhot()
      rrhot=real(irhot)
      irhox=Agrif_Irhox()
      cffx=2./(real(irhox)+1)
      irhoy=Agrif_Irhoy()
      cffy=2./(real(irhoy)+1)

C$OMP BARRIER
C$OMP MASTER

      parentnbstep=Agrif_Parent_Nb_Step()

      IF ((((nbcoarse == 1)     .AND. (nnew == 3)).OR.
     &     ((nbcoarse == irhot) .AND. (nnew /= 3))).AND.
     &     (WTimeindex.NE.parentnbstep)) THEN

        tinterp=1.
#  ifdef MASKING
        Agrif_UseSpecialvalue=.true.
#  endif        
        Agrif_Specialvalue=0.
        Call Agrif_Set_bc(wzid,(/0,0/),InterpolationShouldbemade=.TRUE.)
        Call Agrif_Bc_variable(wzid,calledweight=tinterp,
     &                         procname=w3dinterp)
        Agrif_UseSpecialValue=.false.
        WTimeindex = parentnbstep

      endif
C$OMP END MASTER
C$OMP BARRIER
      
       nbstep3dparent=Agrif_Parent(nbstep3d)
                    
       tinterp = real(nbcoarse)/rrhot
       IF (nnew == 3) tinterp = tinterp - 1./(2.*rrhot)

       ainterp = 0.
 
       IF (nnew == 3) then  
         ainterp(1) = 0.5-tinterp
         ainterp(2) = 0.5+tinterp
       else 
         ainterp(3) = -tinterp
         ainterp(4) = 1.+tinterp
       endif
 
       IF (nbstep3dparent .LE. 1) THEN
         ainterp = 0.
         ainterp(2) = 1.
       ENDIF
 
       if (nbstep3d .LE. (Agrif_irhot()-1)) then
         ainterp = 0.
         ainterp(2) = 1.
       endif
 
       IF ((nnew /= 3).AND.(nbcoarse == irhot)) THEN
         ainterp    = 0.
         ainterp(4) = 1.
       ENDIF
 
#  ifdef AGRIF_OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do k=0,N      
          do i=IstrR,IendR
#   ifdef W_FRC_BRY     
            wbry_south(i,k)=
#   else
            wclm(i,Jstr-1,k)=
#   endif
     &  (cffy*
     &  (ainterp(1)*W_south(i,Jstr-1,k,1)+
     &   ainterp(2)*W_south(i,Jstr-1,k,2)+
     &   ainterp(3)*W_south(i,Jstr-1,k,3)+
     &   ainterp(4)*W_south(i,Jstr-1,k,4)+
     &   ainterp(5)*wz(i,Jstr-1,k,3-nstp)+
     &   ainterp(6)*wz(i,Jstr-1,k,nstp))
     &   +(1.-cffy)*wz(i,Jstr  ,k,nnew))
#   ifdef MASKING     
     &    *rmask(i,Jstr-1)
#   endif
          enddo
        enddo
      endif
#  endif

#  ifdef AGRIF_OBC_NORTH
      if (NORTHERN_EDGE) then      
        do k=0,N      
          do i=IstrR,IendR
#   ifdef W_FRC_BRY     
            wbry_north(i,k)=
#   else   
            wclm(i,Jend+1,k)=
#   endif
     &  (cffy*
     &  (ainterp(1)*W_north(i,Jend+1,k,1)+
     &   ainterp(2)*W_north(i,Jend+1,k,2)+
     &   ainterp(3)*W_north(i,Jend+1,k,3)+
     &   ainterp(4)*W_north(i,Jend+1,k,4)+
     &   ainterp(5)*wz(i,Jend+1,k,3-nstp)+
     &   ainterp(6)*wz(i,Jend+1,k,nstp))
     &   +(1.-cffy)*wz(i,Jend  ,k,nnew))
#   ifdef MASKING     
     &    *rmask(i,Jend+1)
#   endif
          enddo
        enddo
      endif      
#  endif
 
#  ifdef AGRIF_OBC_WEST
      if (WESTERN_EDGE) then      
        do k=0,N      
          do j=JstrR,JendR
#   ifdef W_FRC_BRY     
            wbry_west(j,k)=
#   else     
            wclm(Istr-1,j,k)=
#   endif
     &  (cffx*
     &  (ainterp(1)*W_west(Istr-1,j,k,1)+
     &   ainterp(2)*W_west(Istr-1,j,k,2)+
     &   ainterp(3)*W_west(Istr-1,j,k,3)+
     &   ainterp(4)*W_west(Istr-1,j,k,4)+
     &   ainterp(5)*wz(Istr-1,j,k,3-nstp)+
     &   ainterp(6)*wz(Istr-1,j,k,nstp))
     &   +(1.-cffx)*wz(Istr  ,j,k,nnew))
#   ifdef MASKING     
     &    *rmask(Istr-1,j)
#   endif
          enddo
        enddo
      endif      
#  endif
      
#  ifdef AGRIF_OBC_EAST
      if (EASTERN_EDGE) then      
        do k=0,N      
          do j=JstrR,JendR
#   ifdef W_FRC_BRY  
            wbry_east(j,k)=
#   else
            wclm(Iend+1,j,k)=
#   endif
     &  (cffx*
     &  (ainterp(1)*W_east(Iend+1,j,k,1)+
     &   ainterp(2)*W_east(Iend+1,j,k,2)+
     &   ainterp(3)*W_east(Iend+1,j,k,3)+
     &   ainterp(4)*W_east(Iend+1,j,k,4)+
     &   ainterp(5)*wz(Iend+1,j,k,3-nstp)+
     &   ainterp(6)*wz(Iend+1,j,k,nstp))
     &   +(1.-cffx)*wz(Iend  ,j,k,nnew))        
#   ifdef MASKING     
     &    *rmask(Iend+1,j)
#   endif
          enddo
        enddo
      endif
#  endif
    
      return
      end

!======================================================================

      subroutine w3Dinterp(tabres,i1,i2,j1,j2,k1,k2,before,nb,ndir)

      implicit none
# include "param.h"
# include "ocean3d.h"
# include "scalars.h"
# include "zoom.h"
# ifdef MPI
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif
       integer i1,i2,j1,j2,k1,k2
       real tabres(i1:i2,j1:j2,k1:k2)
       logical before
       integer nb,ndir
       integer nparent
       integer nsource, ndest
       logical :: western_side, eastern_side
       logical :: northern_side,southern_side

       if (before) then

         IF ((iif == 0)) THEN
           nparent = 3
         ELSE
           nparent = nnew
         ENDIF

         tabres(i1:i2,j1:j2,k1:k2) = wz(i1:i2,j1:j2,k1:k2,nparent)

       else

         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)

         IF (nnew == 3) THEN
           nsource = 2
         ELSE
           nsource = 4
         ENDIF
         ndest = nsource - 1

#  ifdef AGRIF_OBC_SOUTH
        if (southern_side) then
#   ifdef MPI
            if (.not.SOUTH_INTER) then
#   endif
                W_south(i1:i2,j1:j2,0:N,ndest)=
     &          W_south(i1:i2,j1:j2,0:N,nsource)
                W_south(i1:i2,j1:j2,0:N,nsource) = tabres
#   ifdef MPI
            endif
#   endif
        endif
#  endif

#  ifdef AGRIF_OBC_NORTH
        if (northern_side) then
#   ifdef MPI
            if (.not.NORTH_INTER) then
#   endif
                W_north(i1:i2,j1:j2,0:N,ndest)=
     &          W_north(i1:i2,j1:j2,0:N,nsource)
                W_north(i1:i2,j1:j2,0:N,nsource) = tabres
#   ifdef MPI
            endif
#   endif
        endif
#  endif

#  ifdef AGRIF_OBC_WEST
        if (western_side) then
#   ifdef MPI
            if (.not.WEST_INTER) then
#   endif
                W_west(i1:i2,j1:j2,0:N,ndest)=
     &          W_west(i1:i2,j1:j2,0:N,nsource)
                W_west(i1:i2,j1:j2,0:N,nsource) = tabres
#   ifdef MPI
            endif
#   endif
        endif
#  endif

#  ifdef AGRIF_OBC_EAST
        if (eastern_side) then
#   ifdef MPI
            if (.not.EAST_INTER) then
#   endif
                W_east(i1:i2,j1:j2,0:N,ndest)=
     &          W_east(i1:i2,j1:j2,0:N,nsource)
                W_east(i1:i2,j1:j2,0:N,nsource) = tabres
#   ifdef MPI
            endif
#   endif
        endif
#  endif

       endif
      return
      end
# endif 

#else
      subroutine zoombc_3D_empty()
      return
      end
#endif /* AGRIF && SOLVE3D */
