! $Id: ana_vmix.F 1458 2014-02-03 15:01:25Z gcambon $
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
#if defined ANA_VMIX && defined SOLVE3D
      subroutine ana_vmix (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call ana_vmix_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine ana_vmix_tile (Istr,Iend,Jstr,Jend)
!
! Set vertical mixing coefficients for momentum "Akv"
! and tracers "Akt" [m^2/s] using analytical expressions.
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
      real cff, kappustr
# ifdef SED_TOY
      real ustar0
# endif
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "mixing.h"
# include "scalars.h"
# ifdef MRL_WCI
#  include "ocean2d.h"
#  include "forces.h"
# endif
# ifdef SED_TOY
#  include "ocean2d.h"
#  include "forces.h"
# endif

# include "compute_auxiliary_bounds.h"

# ifdef CANYON_STRAT
      cff=1./50.                    ! Setup both top and bottom
      do k=1,N-1                    ! boundary layers.
        do j=Jstr,Jend
          do i=Istr,Iend
            Akv(i,j,k)=0.001+0.0095*(
     &                      exp(cff*(z_w(i,j,k)-z_w(i,j,N)))
     &                     +exp(cff*(z_w(i,j,0)-z_w(i,j,k)))
     &                                                     )
#  ifdef TEMPERATURE
            Akt(i,j,k,itemp)=Akt_bak(itemp)
#  endif
          enddo
        enddo
      enddo
# elif defined UPWELLING
      cff=1./150.
      do k=1,N-1
        do j=Jstr,Jend
          do i=Istr,Iend
            Akv(i,j,k)=0.002+0.008*exp(cff*(z_w(i,j,k)-z_w(i,j,N)))
#  ifdef TEMPERATURE
            Akt(i,j,k,itemp)=Akt_bak(itemp)
#  endif
#  ifdef SALINITY
            Akt(i,j,k,isalt)=Akt_bak(isalt)
#  endif
          enddo
        enddo
      enddo
# elif defined SHOREFACE
      do k=1,N-1
        do j=Jstr,Jend
          do i=Istr,Iend
            Akv(i,j,k)=0.011*(h(i,j)+z_w(i,j,k))*
     &                 (1.0-(h(i,j)+z_w(i,j,k))/
     &                  (h(i,j)+zeta(i,j,knew)))
#  ifdef TEMPERATURE
            Akt(i,j,k,itemp)=Akv(i,j,k)
#  endif
#  ifdef SALINITY
            Akt(i,j,k,isalt)=Akv(i,j,k)
#  endif
          enddo
        enddo
      enddo
# elif defined SED_TOY &&(defined SED_TOY_ROUSE || defined SED_TOY_FLOC)
      do k=1,N-1
        do j=Jstr,Jend
          do i=Istr,Iend
              ustar0=sqrt(
     &                  sqrt( (0.5*(bustr(i,j)+bustr(i+1,j)))**2
     &                       +(0.5*(bvstr(i,j)+bvstr(i,j+1)))**2  )
     &                 )
              Akv(i,j,k)=0.41*ustar0 *(h(i,j)+z_w(i,j,k))*
     &                (1.0-(h(i,j)+z_w(i,j,k))/
     &                (h(i,j)+zeta(i,j,knew)))
              Akt(i,j,k,itemp)=Akv(i,j,k)
#  ifdef SALINITY
              Akt(i,j,k,isalt)=Akv(i,j,k)
#  endif
          enddo
        enddo
      enddo
# else
      do k=1,N-1
        do j=Jstr,Jend
          do i=Istr,Iend
            Akv(i,j,k)=Akv_bak
#  ifdef TEMPERATURE
            Akt(i,j,k,itemp)=Akt_bak(itemp)
#  endif
#  ifdef SALINITY
            Akt(i,j,k,isalt)=Akt_bak(isalt)
#  endif
          enddo
        enddo
      enddo
# endif /* CANYON_STRAT */

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend, Akv)
#  ifdef TEMPERATURE
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                     Akt(START_2D_ARRAY,0,itemp))
#  endif
#  ifdef SALINITY
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                     Akt(START_2D_ARRAY,0,isalt))
#  endif
# endif
#else
      subroutine ana_vmix_empty
#endif /* ANA_VMIX */
      return
      end

