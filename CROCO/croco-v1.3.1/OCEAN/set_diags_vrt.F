#include "cppdefs.h"
#if defined DIAGNOSTICS_VRT


      subroutine set_diags_vrt (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call set_diags_vrt_tile (istr,iend,jstr,jend,tile)
      return
      end
 
      subroutine set_diags_vrt_tile (istr,iend,jstr,jend,tile)
!
      implicit none
# include "param.h"
      integer istr,iend,jstr,jend, i,j, ilc, iflux,
     &     imin,imax,jmin,jmax,tile
# ifdef SOLVE3D
     &                       , itrc, k
# endif
      real cff,cff1,cff2,
#ifdef DIAGNOSTICS_UV
     &     wrkXadv(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     wrkYadv(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     wrkHdiff(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     wrkCor(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     wrkPrsgrd(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     wrkHmix(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     wrkVmix(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     wrkVmix2(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     wrkrate(PRIVATE_2D_SCRATCH_ARRAY,2),
# if defined DIAGNOSTICS_BARO
     &     wrkBaro(PRIVATE_2D_SCRATCH_ARRAY,2),
# endif
# if defined M3FAST
     &     wrkfast(PRIVATE_2D_SCRATCH_ARRAY,2),
# endif
# endif
     &     dH(N), jstri(2), istri(2),
     &     jendi(2), iendi(2)
# include "scalars.h"
# include "ncscrum.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"

#if defined DIAGNOSTICS_UV
# include "diagnostics.h"
# endif

#include "diags_vrt.h"

#include "compute_auxiliary_bounds.h"


! ------- ------------- ------ ------ - -----

        jstri(1) = jstr-1
        jstri(2) = jstr

        istri(1) = istr
        istri(2) = istr-1

        jendi(1) = jend
        jendi(2) = jend

        iendi(1) = iend
        iendi(2) = iend

        !jstri(1) = jstrR
        !jstri(2) = jstrR

        !istri(1) = istrR
        !istri(2) = istrR

        !jendi(1) = jendR
        !jendi(2) = jendR

        !iendi(1) = iendR
        !iendi(2) = iendR


#ifdef DIAGNOSTICS_UV

        do itrc=1,2
         do j=jstri(itrc),jendi(itrc)
          do i=istri(itrc),iendi(itrc)


            if (itrc.eq.1) then
              do k=1,N
               dH(k) = (Hz(i,j,k)+Hz(i-1,j,k))
     &                 /(pm(i,j)+pm(i-1,j))
              enddo
            else          
              do k=1,N
               dH(k) = (Hz(i,j,k)+Hz(i,j-1,k))
     &                 /(pn(i,j)+pn(i,j-1))

              enddo
            endif

           ! ------- ------------- ------ ------ - -----
           wrkXadv(i,j,itrc)= MXadv(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkXadv(i,j,itrc)= wrkXadv(i,j,itrc)
     &         +MXadv(i,j,k,itrc)*dH(k)
           enddo
           ! ------- ------------- ------ ------ - -----
           wrkYadv(i,j,itrc)= MYadv(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkYadv(i,j,itrc)= wrkYadv(i,j,itrc)
     &         +MYadv(i,j,k,itrc)*dH(k)
           enddo
           ! ------- ------------- ------ ------ - -----
           wrkHdiff(i,j,itrc)= MHdiff(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkHdiff(i,j,itrc)= wrkHdiff(i,j,itrc)
     &         +MHdiff(i,j,k,itrc)*dH(k)
           enddo
           ! ------- ------------- ------ ------ - -----
           wrkCor(i,j,itrc)= MCor(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkCor(i,j,itrc)= wrkCor(i,j,itrc)
     &         +MCor(i,j,k,itrc)*dH(k)
           enddo
           ! ------- ------------- ------ ------ - -----
           wrkPrsgrd(i,j,itrc)= MPrsgrd(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkPrsgrd(i,j,itrc)= wrkPrsgrd(i,j,itrc)
     &         + MPrsgrd(i,j,k,itrc)*dH(k)
           enddo
           ! ------- ------------- ------ ------ - -----
           wrkHmix(i,j,itrc)= MHmix(i,j,1,itrc,nstp)*dH(1)

           do k=2,N
             wrkHmix(i,j,itrc)= wrkHmix(i,j,itrc)
     &         +MHmix(i,j,k,itrc,nstp)*dH(k)
           enddo
           ! ------- ------------- ------ ------ - -----
           wrkVmix(i,j,itrc)= MVmix(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkVmix(i,j,itrc)= wrkVmix(i,j,itrc)
     &         +MVmix(i,j,k,itrc)*dH(k)
           enddo
           ! ------- ------------- ------ ------ - -----
           wrkVmix2(i,j,itrc)= MVmix2(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkVmix2(i,j,itrc)= wrkVmix2(i,j,itrc)
     &         +MVmix2(i,j,k,itrc)*dH(k)
           enddo
           ! ------- ------------- ------ ------ - -----
           wrkrate(i,j,itrc)= Mrate(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkrate(i,j,itrc)= wrkrate(i,j,itrc)
     &         +Mrate(i,j,k,itrc)*dH(k)
           enddo
           ! ------- ------------- ------ ------ - -----
# if defined DIAGNOSTICS_BARO
           wrkBaro(i,j,itrc)= MBaro(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkBaro(i,j,itrc)= wrkBaro(i,j,itrc)
     &         + MBaro(i,j,k,itrc)*dH(k)
           enddo
# endif
           ! ------- ------------- ------ ------ - -----
# if defined M3FAST
           wrkfast(i,j,itrc)= Mfast(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkfast(i,j,itrc)= wrkfast(i,j,itrc)
     &         + Mfast(i,j,k,itrc)*dH(k)
           enddo
# endif
          enddo
         enddo
        enddo



#endif


      if (WESTERN_EDGE) then     ! Restrict extended ranges one
        imin=istr+1              ! point inward near the physical
      else                       ! boundary. Note that this version
        imin=istr                ! of code is suitable for MPI
      endif                      ! configuration too.

      if (SOUTHERN_EDGE) then
          jmin=jstr+1
      else
          jmin=jstr
      endif



! ------- ------------- ------ ------ - -----
! Compute Vorticity
! ------- ------------- ------ ------ - -----




         do j=jmin,jend
          do i=imin,iend      


          cff = 0.25*(pm(i,j) + pm(i-1,j) + pm(i,j-1) + pm(i-1,j-1))
     &         * 0.25*(pn(i,j) + pn(i-1,j) + pn(i,j-1) + pn(i-1,j-1))

          cff1 = cff
          cff2 = cff


          vrtXadv(i,j) = (wrkXadv(i,j,2) - wrkXadv(i-1,j,2)) * cff1
     &           - (wrkXadv(i,j,1) - wrkXadv(i,j-1,1)) * cff2


          vrtYadv(i,j) = (wrkYadv(i,j,2) - wrkYadv(i-1,j,2)) * cff1
     &           - (wrkYadv(i,j,1) - wrkYadv(i,j-1,1)) * cff2


          vrtHdiff(i,j) = (wrkHdiff(i,j,2)-wrkHdiff(i-1,j,2))* cff1
     &           - (wrkHdiff(i,j,1) - wrkHdiff(i,j-1,1)) * cff2


# if !defined DIAGNOSTICS_UV && (defined UV_HADV_UP3 || defined UV_HADV_UP5)
          vrtHdiff(i,j) = vrtXadv(i,j) + vrtYadv(i,j) - vrtHdiff(i,j)
# endif

          vrtCor(i,j) = (wrkCor(i,j,2) - wrkCor(i-1,j,2)) * cff1
     &           - (wrkCor(i,j,1) - wrkCor(i,j-1,1)) * cff2


          vrtPrsgrd(i,j)=(wrkPrsgrd(i,j,2)-wrkPrsgrd(i-1,j,2))*cff1
     &           - (wrkPrsgrd(i,j,1) - wrkPrsgrd(i,j-1,1)) * cff2

#ifdef DIAGNOSTICS_UV
          vrtHmix(i,j) =
     &           (wrkHmix(i,j,2) - wrkHmix(i-1,j,2)) * cff1
     &         - (wrkHmix(i,j,1) - wrkHmix(i,j-1,1)) * cff2
#else
          vrtHmix(i,j) =
     &           (wrkHmix(i,j,2,nstp) - wrkHmix(i-1,j,2,nstp)) * cff1
     &         - (wrkHmix(i,j,1,nstp) - wrkHmix(i,j-1,1,nstp)) * cff2
#endif

          vrtVmix(i,j) = (wrkVmix(i,j,2) - wrkVmix(i-1,j,2)) * cff1
     &           - (wrkVmix(i,j,1) - wrkVmix(i,j-1,1)) * cff2


          vrtVmix2(i,j) = (wrkVmix2(i,j,2)-wrkVmix2(i-1,j,2))*cff1
     &           - (wrkVmix2(i,j,1) - wrkVmix2(i,j-1,1)) * cff2


          vrtrate(i,j) = (wrkrate(i,j,2) - wrkrate(i-1,j,2)) * cff1
     &           - (wrkrate(i,j,1) - wrkrate(i,j-1,1)) * cff2

          vrtWind(i,j) = (wrkWind(i,j,2) - wrkWind(i-1,j,2)) * cff1
     &           - (wrkWind(i,j,1) - wrkWind(i,j-1,1)) * cff2

          vrtDrag(i,j) = (wrkDrag(i,j,2) - wrkDrag(i-1,j,2)) * cff1
     &           - (wrkDrag(i,j,1) - wrkDrag(i,j-1,1)) * cff2

# if defined DIAGNOSTICS_BARO
          vrtBaro(i,j)=(wrkBaro(i,j,2)-wrkBaro(i-1,j,2))*cff1
     &           - (wrkBaro(i,j,1) - wrkBaro(i,j-1,1)) * cff2
# endif

# if defined M3FAST
          vrtfast(i,j)=(wrkfast(i,j,2)-wrkfast(i-1,j,2))*cff1
     &           - (wrkfast(i,j,1) - wrkfast(i,j-1,1)) * cff2
# endif
          enddo
         enddo

! ------- ------------- ------ ------ - -----
      return
      end
#else /* DIAGNOSTICS_VRT*/
      subroutine set_diags_vrt_empty
      end
#endif /* DIAGNOSTICS_VRT */










