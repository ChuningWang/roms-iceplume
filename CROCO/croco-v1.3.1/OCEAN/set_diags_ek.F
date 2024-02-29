#include "cppdefs.h"
#if defined DIAGNOSTICS_EK


      subroutine set_diags_ek (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call set_diags_ek_tile (istr,iend,jstr,jend,tile)
      return
      end
 
      subroutine set_diags_ek_tile (istr,iend,jstr,jend,tile)
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


#include "diags_ek.h"
#include "compute_auxiliary_bounds.h"




! ------- ------------- ------ ------ - -----

        do j=JstrR,JendR
          do i=IstrR,IendR


          ekHadv(i,j) = 0.5 * (ekwrkHadv(i,j,1) + ekwrkHadv(i+1,j,1))
     &                + 0.5 * (ekwrkHadv(i,j,2) + ekwrkHadv(i,j+1,2))


          ekHdiff(i,j) = 0.5*(ekwrkHdiff(i,j,1)+ekwrkHdiff(i+1,j,1))
     &                + 0.5*(ekwrkHdiff(i,j,2) + ekwrkHdiff(i,j+1,2))


# if !defined DIAGNOSTICS_UV && !defined DIAGNOSTICS_PV && (defined UV_HADV_UP3 || defined UV_HADV_UP5)
          ekHdiff(i,j) = ekHadv(i,j) - ekHdiff(i,j)
# endif       
                 
          ekVadv(i,j) = 0.5 * (ekwrkVadv(i,j,1) + ekwrkVadv(i+1,j,1))
     &                + 0.5 * (ekwrkVadv(i,j,2) + ekwrkVadv(i,j+1,2))


          ekCor(i,j) = 0.5 * (ekwrkCor(i,j,1) + ekwrkCor(i+1,j,1))
     &                + 0.5 * (ekwrkCor(i,j,2) + ekwrkCor(i,j+1,2))


          ekPrsgrd(i,j)=0.5*(ekwrkPrsgrd(i,j,1)+ekwrkPrsgrd(i+1,j,1))
     &              + 0.5*(ekwrkPrsgrd(i,j,2) + ekwrkPrsgrd(i,j+1,2))


          ekHmix(i,j) =
     &         0.5 * (ekwrkHmix(i,j,1,nstp) + ekwrkHmix(i+1,j,1,nstp))
     &       + 0.5 * (ekwrkHmix(i,j,2,nstp) + ekwrkHmix(i,j+1,2,nstp))


          ekVmix(i,j) = 0.5 * (ekwrkVmix(i,j,1) + ekwrkVmix(i+1,j,1))
     &                + 0.5 * (ekwrkVmix(i,j,2) + ekwrkVmix(i,j+1,2))


          ekrate(i,j) = 0.5 * (ekwrkrate(i,j,1) + ekwrkrate(i+1,j,1))
     &                + 0.5 * (ekwrkrate(i,j,2) + ekwrkrate(i,j+1,2))


          ekvol(i,j) = 0.5 * (ekwrkvol(i,j,1) + ekwrkvol(i+1,j,1))
     &                + 0.5 * (ekwrkvol(i,j,2) + ekwrkvol(i,j+1,2))


          ekVmix2(i,j)=0.5*(ekwrkVmix2(i,j,1) + ekwrkVmix2(i+1,j,1))
     &                +0.5*(ekwrkVmix2(i,j,2) + ekwrkVmix2(i,j+1,2))


          ekWind(i,j) = 0.5 * (ekwrkWind(i,j,1) + ekwrkWind(i+1,j,1))
     &                + 0.5 * (ekwrkWind(i,j,2) + ekwrkWind(i,j+1,2))


          ekDrag(i,j) = 0.5 * (ekwrkDrag(i,j,1) + ekwrkDrag(i+1,j,1))
     &                + 0.5 * (ekwrkDrag(i,j,2) + ekwrkDrag(i,j+1,2))
# if defined DIAGNOSTICS_BARO
          ekBaro(i,j)=0.5*(ekwrkBaro(i,j,1)+ekwrkBaro(i+1,j,1))
     &              + 0.5*(ekwrkBaro(i,j,2) + ekwrkBaro(i,j+1,2))
# endif
# if defined M3FAST
          ekfast(i,j)=0.5*(ekwrkfast(i,j,1)+ekwrkfast(i+1,j,1))
     &              + 0.5*(ekwrkfast(i,j,2) + ekwrkfast(i,j+1,2))
# endif

# if defined DIAGNOSTICS_EK_MLD


          ekHadv_mld(i,j) = 0.5 * (ekwrkHadv_mld(i,j,1) + ekwrkHadv(i+1,j,1))
     &                + 0.5 * (ekwrkHadv_mld(i,j,2) + ekwrkHadv_mld(i,j+1,2))


          ekHdiff_mld(i,j) = 0.5*(ekwrkHdiff_mld(i,j,1)+ekwrkHdiff(i+1,j,1))
     &                + 0.5*(ekwrkHdiff_mld(i,j,2) + ekwrkHdiff_mld(i,j+1,2))


# if !defined DIAGNOSTICS_UV && !defined DIAGNOSTICS_PV && (defined UV_HADV_UP3 || defined UV_HADV_UP5)
          ekHdiff_mld(i,j) = ekHadv_mld(i,j) - ekHdiff_mld(i,j)
# endif       
                 
          ekVadv_mld(i,j) =
     &         0.5 * (ekwrkVadv_mld(i,j,1) + ekwrkVadv(i+1,j,1))
     &       + 0.5 * (ekwrkVadv_mld(i,j,2) + ekwrkVadv_mld(i,j+1,2))


          ekCor_mld(i,j) =
     &         0.5 * (ekwrkCor_mld(i,j,1) + ekwrkCor(i+1,j,1))
     &       + 0.5 * (ekwrkCor_mld(i,j,2) + ekwrkCor_mld(i,j+1,2))


          ekPrsgrd_mld(i,j)=
     &         0.5*(ekwrkPrsgrd_mld(i,j,1)+ekwrkPrsgrd(i+1,j,1))
     &       + 0.5*(ekwrkPrsgrd_mld(i,j,2) + ekwrkPrsgrd_mld(i,j+1,2))


          ekHmix_mld(i,j) =
     &         0.5 * (ekwrkHmix_mld(i,j,1,nstp) + ekwrkHmix(i+1,j,1,nstp))
     &       + 0.5 * (ekwrkHmix_mld(i,j,2,nstp) + ekwrkHmix_mld(i,j+1,2,nstp))


          ekVmix_mld(i,j) =
     &         0.5 * (ekwrkVmix_mld(i,j,1) + ekwrkVmix(i+1,j,1))
     &       + 0.5 * (ekwrkVmix_mld(i,j,2) + ekwrkVmix_mld(i,j+1,2))


          ekrate_mld(i,j) =
     &         0.5 * (ekwrkrate_mld(i,j,1) + ekwrkrate(i+1,j,1))
     &       + 0.5 * (ekwrkrate_mld(i,j,2) + ekwrkrate_mld(i,j+1,2))


          ekvol_mld(i,j) =
     &         0.5 * (ekwrkvol_mld(i,j,1) + ekwrkvol(i+1,j,1))
     &       + 0.5 * (ekwrkvol_mld(i,j,2) + ekwrkvol_mld(i,j+1,2))


          ekVmix2_mld(i,j)=
     &         0.5*(ekwrkVmix2_mld(i,j,1) + ekwrkVmix2(i+1,j,1))
     &        +0.5*(ekwrkVmix2_mld(i,j,2) + ekwrkVmix2_mld(i,j+1,2))

# if defined DIAGNOSTICS_BARO
          ekBaro_mld(i,j)=
     &         0.5*(ekwrkBaro_mld(i,j,1)+ekwrkBaro(i+1,j,1))
     &       + 0.5*(ekwrkBaro_mld(i,j,2) + ekwrkBaro_mld(i,j+1,2))
# endif
# endif

          enddo
         enddo

! ------- ------------- ------ ------ - -----
      return
      end
#else /* DIAGNOSTICS_EK*/
      subroutine set_diags_ek_empty
      end
#endif /* DIAGNOSTICS_EK */










