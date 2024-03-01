! $Id: nbq_bry_store 1619 2015-01-07 13:53:03Z marchesiello $
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
!
!======================================================================
!                   subroutine nbq_bry_store
!======================================================================
!
#ifndef CHILD
!
# include "cppdefs.h"
# if defined OBC_NBQ && defined OBC_NBQORLANSKI
      subroutine nbq_bry_store_tile (Istr,Iend,Jstr,Jend)
#  ifdef AGRIF      
      use AGRIF_Util
      integer Istr,Iend,Jstr,Jend
      if (AGRIF_Root()) then
        call nbq_bry_store_parent_tile (Istr,Iend,Jstr,Jend)
      else
        call nbq_bry_store_child_tile (Istr,Iend,Jstr,Jend)
      endif
      return
      end
!
! PARENT
!
      subroutine nbq_bry_store_parent_tile (Istr,Iend,Jstr,Jend)
#  endif
!
! Set boundary arrays for U-component momentum qdmu_nbq for the parent grid.
!
# endif /* defined OBC_NBQ && defined OBC_NBQORLANSKI */
#else
# if defined OBC_NBQ && defined OBC_NBQORLANSKI
!
! CHILD
!
      subroutine nbq_bry_store_child_tile (Istr,Iend,Jstr,Jend)
!
! Set boundary arrays for U-component momentum
! qdmu_nbq for the child grid.
!
# endif /* defined OBC_NBQ && defined OBC_NBQORLANSKI */
#endif /* CHILD */
!
#if defined OBC_NBQ && defined OBC_NBQORLANSKI
!
! Common Code
!
# include "set_obc_definitions.h"
!
!----------------------------------------------------------------------
!  Store boundary values of nbq variables 
!  for use in radiative conditions
!----------------------------------------------------------------------
!
      implicit none          
      integer Istr,Iend,Jstr,Jend, i,j,k, itrc
# include "param.h"
# include "scalars.h"
# include "nbq.h"
!
# include "compute_auxiliary_bounds.h"
!
# ifdef OBC_COM_WEST
      if (WESTERN_EDGE) then
        do k=1,N
          do j=Jstr-1,Jend+1
            qdmu_nbq_west(j,k,1)=qdmu_nbq(Istr  ,j,k)
            qdmu_nbq_west(j,k,2)=qdmu_nbq(Istr+1,j,k)
            qdmv_nbq_west(j,k,1)=qdmv_nbq(Istr-1,j,k)
            qdmv_nbq_west(j,k,2)=qdmv_nbq(Istr  ,j,k)
#  ifdef NBQ
            qdmw_nbq_west(j,k,1)=qdmw_nbq(Istr-1,j,k)
            qdmw_nbq_west(j,k,2)=qdmw_nbq(Istr  ,j,k)
             rho_nbq_west(j,k,1)= rho_nbq(Istr-1,j,k)
             rho_nbq_west(j,k,2)= rho_nbq(Istr  ,j,k)
#  endif
          enddo
        enddo
#  ifdef NBQ
        do j=Jstr-1,Jend+1
          qdmw_nbq_west(j,0,1)=qdmw_nbq(Istr-1,j,0)
          qdmw_nbq_west(j,0,2)=qdmw_nbq(Istr  ,j,0)
        enddo
#  endif
      endif
# endif /* OBC_COM_WEST */
!
# ifdef OBC_COM_EAST
      if (EASTERN_EDGE) then
        do k=1,N
          do j=Jstr-1,Jend+1
            qdmu_nbq_east(j,k,1)=qdmu_nbq(Iend+1,j,k)
            qdmu_nbq_east(j,k,2)=qdmu_nbq(Iend  ,j,k)
            qdmv_nbq_east(j,k,1)=qdmv_nbq(Iend+1,j,k)
            qdmv_nbq_east(j,k,2)=qdmv_nbq(Iend  ,j,k)
#  ifdef NBQ
            qdmw_nbq_east(j,k,1)=qdmw_nbq(Iend+1,j,k)
            qdmw_nbq_east(j,k,2)=qdmw_nbq(Iend  ,j,k)
             rho_nbq_east(j,k,1)= rho_nbq(Iend+1,j,k)
             rho_nbq_east(j,k,2)= rho_nbq(Iend  ,j,k)
#  endif
          enddo
        enddo
#  ifdef NBQ
        do j=Jstr-1,Jend+1
          qdmw_nbq_east(j,0,1)=qdmw_nbq(Iend+1,j,0)
          qdmw_nbq_east(j,0,2)=qdmw_nbq(Iend  ,j,0)
        enddo
#  endif
      endif
# endif /* OBC_COM_EAST */
!
# ifdef OBC_COM_SOUTH
      if (SOUTHERN_EDGE) then
        do k=1,N
          do i=Istr-1,Iend+1
            qdmu_nbq_south(i,k,1)=qdmu_nbq(i,Jstr-1,k)
            qdmu_nbq_south(i,k,2)=qdmu_nbq(i,Jstr  ,k)
            qdmv_nbq_south(i,k,1)=qdmv_nbq(i,Jstr  ,k)
            qdmv_nbq_south(i,k,2)=qdmv_nbq(i,Jstr+1,k)
#  ifdef NBQ
            qdmw_nbq_south(i,k,1)=qdmw_nbq(i,Jstr-1,k)
            qdmw_nbq_south(i,k,2)=qdmw_nbq(i,Jstr  ,k)
             rho_nbq_south(i,k,1)= rho_nbq(i,Jstr-1,k)
             rho_nbq_south(i,k,2)= rho_nbq(i,Jstr  ,k)
#  endif
          enddo
        enddo
#  ifdef NBQ
        do i=Istr-1,Iend+1
          qdmw_nbq_south(i,0,1)=qdmw_nbq(i,Jstr-1,0)
          qdmw_nbq_south(i,0,2)=qdmw_nbq(i,Jstr  ,0)
        enddo
#  endif
      endif
# endif /* OBC_COM_SOUTH */
!
# ifdef OBC_COM_NORTH
      if (NORTHERN_EDGE) then
        do k=1,N
          do i=Istr-1,Iend+1
            qdmu_nbq_north(i,k,1)=qdmu_nbq(i,Jend+1,k)
            qdmu_nbq_north(i,k,2)=qdmu_nbq(i,Jend  ,k)
            qdmv_nbq_north(i,k,1)=qdmv_nbq(i,Jend+1,k)
            qdmv_nbq_north(i,k,2)=qdmv_nbq(i,Jend  ,k)
#  ifdef NBQ
            qdmw_nbq_north(i,k,1)=qdmw_nbq(i,Jend+1,k)
            qdmw_nbq_north(i,k,2)=qdmw_nbq(i,Jend  ,k)
             rho_nbq_north(i,k,1)= rho_nbq(i,Jend+1,k)
             rho_nbq_north(i,k,2)= rho_nbq(i,Jend  ,k)
#  endif
          enddo
        enddo
#  ifdef NBQ
        do i=Istr-1,Iend+1
          qdmw_nbq_north(i,0,1)=qdmw_nbq(i,Jend+1,0)
          qdmw_nbq_north(i,0,2)=qdmw_nbq(i,Jend  ,0)
        enddo
#  endif
      endif
# endif /* OBC_COM_NORTH */
!
      return
      end
#else
# ifndef CHILD
      subroutine nbq_bry_store_parent_empty
      end
# else
      subroutine nbq_bry_store_child_empty
      end
# endif
#endif /* OBC_NBQ && OBC_NBQORLANSKI */
#ifndef CHILD
# define CHILD
# ifdef AGRIF
#  include "nbq_bry_store.F"
# endif
# undef CHILD
#endif  /* !CHILD */

