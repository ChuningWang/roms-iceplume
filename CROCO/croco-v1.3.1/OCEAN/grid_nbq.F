! $Id: cppdefs.h 1628 2015-01-10 13:53:00Z marchesiello $
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
#ifdef NBQ
!
      subroutine grid_nbq_tile(Istr,Iend,Jstr,Jend,
     &                         Hzw_half_nbq_inv,Hzr_half_nbq_inv,
     &                         Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v,
     &                         Hzu_half_qdmu, Hzv_half_qdmv )     
!
!**********************************************************************
!
!                 Pre-computations for NH / NBQ modes
!
!**********************************************************************
!
      implicit none

      integer Istr,Iend,Jstr,Jend
      integer imin,imax,jmin,jmax
# include "param.h"
# include "scalars.h"
# include "private_scratch.h"
# include "grid.h"
# include "ocean3d.h"
# include "nbq.h"

      real Hzw_half_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,0:N)
      real Hzr_half_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,N)
      real Hzw_half_nbq_inv_u(PRIVATE_2D_SCRATCH_ARRAY,0:N)
      real Hzw_half_nbq_inv_v(PRIVATE_2D_SCRATCH_ARRAY,0:N)
      real Hzu_half_qdmu(PRIVATE_2D_SCRATCH_ARRAY,0:N)
      real Hzv_half_qdmv(PRIVATE_2D_SCRATCH_ARRAY,0:N)

      integer i,j,k,it
      double precision val1, val2

# include "compute_auxiliary_bounds.h"
!
# ifndef NBQ_MASS
#  define Hzr_half_nbq Hz
# endif

!
!----------------------------------------------------------------------
! Sets indices
!----------------------------------------------------------------------
!
#  ifdef EW_PERIODIC
      imin=Istr-2
      imax=Iend+2
#  else
      if (WESTERN_EDGE) then
        imin=Istr-1
      else
        imin=Istr-2
      endif
      if (EASTERN_EDGE) then
        imax=Iend+1
      else
        imax=Iend+2
      endif
#  endif
#  ifdef NS_PERIODIC
      jmin=Jstr-2
      jmax=Jend+2
#  else
      if (SOUTHERN_EDGE) then
        jmin=Jstr-1
      else
        jmin=Jstr-2
      endif
      if (NORTHERN_EDGE) then
        jmax=Jend+1
      else
        jmax=Jend+2
      endif
#  endif    
!
!----------------------------------------------------------------------
!  Compute other vertical grid variables
!          at m
!----------------------------------------------------------------------
!	
      do k=1,N-1
        do j=jmin,jmax
          do i=imin,imax
            Hzw_half_nbq(i,j,k)=z_r(i,j,k+1)-z_r(i,j,k)
          enddo
        enddo
      enddo
      
      do j=jmin,jmax
        do i=imin,imax
          Hzw_half_nbq(i,j,0)=z_r(i,j,1)-z_w(i,j,0)
          Hzw_half_nbq(i,j,N)=z_w(i,j,N)-z_r(i,j,N)
        enddo
      enddo 

      do k=1,N
        do j=JstrV-2,Jend+1
          do i=IstrU-2,Iend+1
            Hzr_half_nbq_inv(i,j,k)=1./max(1.e-30,Hzr(i,j,k))
# ifdef MASKING
     &                                            *rmask(i,j)
# endif
          enddo
        enddo
      enddo

      do k=1,N
        do j=Jstr,Jend
          do i=IstrU,Iend+1
            Hzu_half_qdmu(i,j,k)=0.5*(Hzr(i-1,j,k)+Hzr(i,j,k))*pm_u(i,j)
# ifdef MASKING
     &                                                       *umask(i,j)
# endif  
          enddo
        enddo
      enddo 
      do k=1,N
        do j=JstrV,Jend+1
          do i=Istr,Iend
            Hzv_half_qdmv(i,j,k)=0.5*(Hzr(i,j-1,k)+Hzr(i,j,k))*pn_v(i,j)
# ifdef MASKING
     &                                                       *vmask(i,j)
# endif  
          enddo
        enddo
      enddo 

      do k=0,N
        do j=JstrV-2,Jend+1
          do i=IstrU-2,Iend+1
            Hzw_half_nbq_inv(i,j,k)=1./max(1.e-30,Hzw_half_nbq(i,j,k))
# ifdef MASKING
     &                                                      *rmask(i,j)
# endif
          enddo
        enddo
      enddo

      do k=0,N
        do j=JstrV-2,Jend+1
          do i=IstrU-1,Iend+1
            Hzw_half_nbq_inv_u(i,j,k)=0.5/max(1.e-30,
     &                                        Hzw_half_nbq(i  ,j,k)+
     &                                        Hzw_half_nbq(i-1,j,k)) 
# ifdef MASKING
     &                                                   *umask(i,j)
# endif      
          enddo
        enddo
      enddo
              
      do k=0,N
        do j=JstrV-1,Jend+1
          do i=IstrU-2,Iend+1
            Hzw_half_nbq_inv_v(i,j,k)=0.5/max(1.e-30,
     &                                        Hzw_half_nbq(i,j  ,k)+
     &                                        Hzw_half_nbq(i,j-1,k)) 
# ifdef MASKING	
     &                                                   *vmask(i,j)
# endif      	  
          enddo
        enddo
      enddo          
        
      return
      end subroutine grid_nbq_tile

#else
      subroutine grid_nbq_tile_empty
      return
      end 
#endif


