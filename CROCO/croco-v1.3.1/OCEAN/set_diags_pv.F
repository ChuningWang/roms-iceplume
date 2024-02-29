#include "cppdefs.h"
#if defined DIAGNOSTICS_PV


      subroutine set_diags_pv (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call set_diags_pv_tile (istr,iend,jstr,jend,tile)
      return
      end
 
      subroutine set_diags_pv_tile (istr,iend,jstr,jend,tile)
!
!
      implicit none
# include "param.h"
      integer istr,iend,jstr,jend, i,j, ilc, iflux,
     &     imin,imax,jmin,jmax,tile
# ifdef SOLVE3D
     &                       , itrc, k
# endif
      real cff,cff1,cff2,
     &     dH(N), jstri(2), istri(2),
     &     jendi(2), iendi(2)
# include "scalars.h"
# include "ncscrum.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"

# ifdef DIAGNOSTICS_UV
#  include "diagnostics.h"
# endif
# ifdef DIAGNOSTICS_EK
#  include "diags_ek.h"
# endif

#include "diags_pv.h"
! diagnostic dissipation NJAL August 2017
# if defined DIAGNOSTICS_DISS       
      real alfa(PRIVATE_2D_SCRATCH_ARRAY), 
     &     beta(PRIVATE_2D_SCRATCH_ARRAY)
#endif
#include "compute_auxiliary_bounds.h"

        do itrc=1,NTA
          do k=1,N
            do j=JstrR,JendR
              do i=IstrR,IendR 
                Trhs(i,j,k,itrc) = TForc(i,j,k,itrc)
     &                           + THmix(i,j,k,itrc)
     &                           + TVmix(i,j,k,itrc)
              enddo
            enddo
          enddo
        enddo

# if defined DIAGNOSTICS_DISS
      do k=1,N
          call alfabeta_k_tile(Istr,Iend,Jstr,Jend,k,alfa,beta)
          do j=JstrR,JendR
              do i=IstrR,IendR
                  Trhs(i,j,k,1) = -alfa(i,j)*Trhs(i,j,k,1)
                  Trhs(i,j,k,2) = beta(i,j)*Trhs(i,j,k,2)
              enddo
          enddo
      enddo
# endif /* DIAGNOSTICS_DISS */




      return
      end
#else /* DIAGNOSTICS_PV*/
      subroutine set_diags_pv_empty
      end
#endif /* DIAGNOSTICS_PV */










