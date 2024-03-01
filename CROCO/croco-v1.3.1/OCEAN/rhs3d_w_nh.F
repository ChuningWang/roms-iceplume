#include "cppdefs.h"

#define W_HADV
#define W_VADV

#ifdef NBQ_FREESLIP
# define NBQ_WFREESLIP
#else
# undef  NBQ_WFREESLIP
#endif

#if defined SOLVE3D && defined NBQ

      SUBROUTINE rhs3d_w_nh(tile)

!======================================================================
!                   ***  Subroutine RHS3D  ***
! NBQ mode :  compute right-hand-side for the vertical velocity wz 
!======================================================================
! History :   2016-11  (F. Lemarié)  Original code
!----------------------------------------------------------------------

      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
!
!   A3d(1,1,trd), A3d(1,2,trd) contain ru,rv 
!   A3d(1,3,trd) contains rw 
!   A3d(1,3,trd) and A3d(1,4,trd) is available at this point
!
      CALL rhs3d_w_tile (Istr,Iend,Jstr,Jend,  A3d(1,3,trd),
     &             A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd),
     &             A2d(1,4,trd), A2d(1,5,trd), A2d(1,6,trd),
     &             A2d(1,7,trd), A3d(1,4,trd), A3d(1,5,trd))
      return
      END
      
!======================================================================
      
      SUBROUTINE rhs3d_w_tile (istr, iend,   jstr,   jend, rw ,
     &                           CF,   FC,  We_r ,    lap, Hzw, 
     &                          WFx,  WFe, HUon_w, HVom_w       )      
      IMPLICIT NONE
      !!
      INTEGER  :: Istr,Iend,Jstr,Jend
      INTEGER  :: i,j,k,kp
      INTEGER  :: imin,imax,jmin,jmax
      real        rr,rrm,rrp,limiter2,Cr,cdt,cu
      integer     kp2,kp1,km1
      !!
# include "param.h"
      REAL     :: rw     (PRIVATE_2D_SCRATCH_ARRAY,0:N  )
      REAL     :: CF     (PRIVATE_1D_SCRATCH_ARRAY,0:N  )
      REAL     :: FC     (PRIVATE_1D_SCRATCH_ARRAY,0:N+1)  
      REAL     :: We_r   (PRIVATE_1D_SCRATCH_ARRAY,0:N+1)        
      REAL     :: lap    (PRIVATE_2D_SCRATCH_ARRAY      )
      REAL     :: Hzw    (PRIVATE_1D_SCRATCH_ARRAY,0:N  )      
      REAL     :: WFx    (PRIVATE_2D_SCRATCH_ARRAY      ) 
      REAL     :: WFe    (PRIVATE_2D_SCRATCH_ARRAY      )
      REAL     :: HUon_w (PRIVATE_2D_SCRATCH_ARRAY,0:N  )
      REAL     :: HVom_w (PRIVATE_2D_SCRATCH_ARRAY,0:N  )
      !!    
      REAL     :: gamma, epsil, cff, cff1, cff2, cff3, Omeg_r
      REAL     :: wz_bot,wz_sfc
      PARAMETER (gamma=0.25  )
      PARAMETER (epsil=1.E-16)
      !!
# include "grid.h"
# include "ocean3d.h"
# include "coupling.h"
# include "forces.h"
# include "scalars.h"  
# include "nbq.h"
!
!--------------------------------------------------------------------
! Definition of flux operators: 1st, 2nd, 3rd, 4th, 5th or 6th order,
! used in UP5 and C6 advection schemes (and order degradation near 
! land masks). cdiff is part of laplacian diffusion in flux1 (used 
! near mask): 
!    0 --> flux1=flux2 (second order C2 advection scheme)
!    1 --> flux1 gives 1st order monotonic UP1 advection scheme
!--------------------------------------------------------------------
!
# if defined W_HADV_WENO5 || defined W_HADV_C6 || defined W_HADV_UP5 \
  || defined W_VADV_WENO5 || defined W_VADV_C6 || defined W_VADV_UP5
      REAL    :: q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2
      REAL    :: ua, vel, cdiff, cdif
      REAL    :: flux1, flux2, flux3, flux4, flux5, flux6
      REAL    :: flx2, flx3, flx4, flx5
      REAL    :: mask0, mask1, mask2, mask3

      flux2(q_im1, q_i, ua, cdiff) = 0.5*( q_i + q_im1 )

      flux1(q_im1, q_i, ua, cdiff) = flux2(q_im1, q_i, ua, cdiff) -
     &      0.5*cdiff*sign(1.,ua)*(q_i-q_im1)

      flux4(q_im2, q_im1, q_i, q_ip1, ua) = 
     &      ( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0

      flux3(q_im2, q_im1, q_i, q_ip1, ua) =
     &      flux4(q_im2, q_im1, q_i, q_ip1, ua) +
     &      sign(1.,ua)*((q_ip1 -
     &      q_im2)-3.*(q_i-q_im1))/12.0

      flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =
     &      ( 37.*(q_i+q_im1) - 8.*(q_ip1+q_im2)
     &      +(q_ip2+q_im3) )/60.0

      flux5(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =
     &      flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua)
     &      -sign(1.,ua)*( 
     &      (q_ip2-q_im3)-5.*(q_ip1-q_im2)+10.*(q_i-q_im1) )/60.0
# endif
# if  defined W_HADV_WENO5 || defined W_VADV_WENO5
      REAL    :: flux3_weno, flux5_weno
# endif

# include "compute_auxiliary_bounds.h"
 
# ifdef EW_PERIODIC 
#   define IU_RANGE Istr,Iend
#   define IV_RANGE Istr,Iend
# else
#   define IU_RANGE Istr,IendR
#   define IV_RANGE IstrR,IendR
# endif

# ifdef NS_PERIODIC 
#   define JU_RANGE Jstr,Jend
#   define JV_RANGE Jstr,Jend
# else
#   define JU_RANGE JstrR,JendR
#   define JV_RANGE Jstr,JendR
# endif

      DO j=Jstr,Jend
        DO i=Istr,Iend+1
          HUon_w(i,j,0)=0.5*HUon(i,j,1)                
          HUon_w(i,j,N)=0.5*HUon(i,j,N) 
        ENDDO
      ENDDO
      DO j=Jstr,Jend+1
        DO i=Istr,Iend
          HVom_w(i,j,0)=0.5*HVom(i,j,1)             
          HVom_w(i,j,N)=0.5*HVom(i,j,N) 
        ENDDO
      ENDDO 
      DO k=1,N-1
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            HUon_w(i,j,k)=0.5*( HUon(i,j,k)+HUon(i,j,k+1) )
          ENDDO
        ENDDO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            HVom_w(i,j,k)=0.5*( HVom(i,j,k)+HVom(i,j,k+1) )
          ENDDO
        ENDDO         
      ENDDO

      DO k=0,N
         DO j=Jstr,Jend   
            DO i=Istr,Iend        
              rw(i,j,k)=0.            
            ENDDO
         ENDDO
      ENDDO 

# if defined UV_ADV && defined W_HADV
!
!=======================================================================
!
!  Horizontal advection
!
!=======================================================================
!
#  if defined W_HADV_TVD || defined W_VADV_TVD
      if(nrhs.eq.3) then
        cdt = dt                     !<-- Corrector
      elseif(FIRST_TIME_STEP) then
        cdt = 0.5*dt
      else
        cdt = (1.-1./6.)*dt          !<-- Predictor
      endif
#  endif

#  if defined W_HADV_WENO5 || defined W_HADV_C6 || defined W_HADV_UP5
#   ifdef NS_PERIODIC
      jmin=1
      jmax=LOCALMM+1
#   else
#    ifdef MPI
      if (SOUTH_INTER) then
        jmin=1
      else
        jmin=3
      endif
      if (NORTH_INTER) then
        jmax=Mmmpi+1
      else
        jmax=Mmmpi-1
      endif
#    else
      jmin=3
      jmax=Mm-1
#    endif
#   endif
#   ifdef EW_PERIODIC
      imin=1
      imax=LOCALLM+1
#   else
#    ifdef MPI
      if (WEST_INTER) then
        imin=1
      else
        imin=3
      endif
      if (EAST_INTER) then
        imax=Lmmpi+1
      else
        imax=Lmmpi-1
      endif
#    else
      imin=3
      imax=Lm-1
#    endif
#   endif

#  else /* !W_HADV_WENO5 */
#   ifndef EW_PERIODIC
      IF (WESTERN_EDGE) THEN       
         imin=istr                  
      ELSE                        
         imin=istr-1             
      ENDIF                      
      IF (EASTERN_EDGE) THEN      
         imax=iend                 
      ELSE                  
         imax=iend+1           
      ENDIF
#   else
      imin=istr-1
      imax=iend+1
#   endif
#   ifndef NS_PERIODIC
      IF (SOUTHERN_EDGE) THEN
         jmin=jstr
      ELSE
         jmin=jstr-1
      ENDIF
      IF (NORTHERN_EDGE) then
         jmax=jend
      ELSE
         jmax=jend+1
      ENDIF
#   else
      jmin=jstr-1
      jmax=jend+1
#   endif
#  endif /* W_HADV_WENO5 */

 
      DO k=0,N             

#  ifdef W_HADV_C2
!
! === C2 horizontal advection scheme ===
!
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            WFx(i,j)=0.5*( wz(i,j,k,nrhs)+wz(i-1,j,k,nrhs)
     &                                     )*Huon_w(i,j,k)
          ENDDO    
        ENDDO
        DO j=jstr,jend+1
          DO i=istr,iend
            WFe(i,j)=0.5*( wz(i,j,k,nrhs)+wz(i,j-1,k,nrhs)
     &                                     )*Hvom_w(i,j,k)
          ENDDO
        ENDDO

#  elif defined W_HADV_TVD
!
! === TVD horizontal advection scheme ===
!
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            WFx(i,j)=0.5*( wz(i,j,k,nrhs)+wz(i-1,j,k,nrhs)
     &                                     )*Huon_w(i,j,k)
          ENDDO    
        ENDDO

        DO j=Jstr,Jend
          ! wz(i-2,..) out of bounds
          i = Istr
            if (k==0) then
              Hzw (i,k) = 0.25 * (HZR(i-1,j,k+1) + HZR(i  ,j,k+1))
            elseif (k==N) then
              Hzw (i,k) = 0.25 * (HZR(i-1,j,k  ) + HZR(i  ,j,k  ))
            else
              Hzw (i,k) = 0.25 * (HZR(i  ,j,k  ) + HZR(i  ,j,k+1)
     &                           +HZR(i-1,j,k  ) + HZR(i-1,j,k+1))
            endif
            cff = Huon_w(i,j,k)
            cu = pn_u(i,j)*cdt*cff
     &           /Hzw(i,k)*pm(i,j)

            rrp= (wz(i+1,j,k,nstp)-wz(i  ,j,k,nstp))
#   ifdef MASKING
     &                                 *umask(i+1,j)
#   endif
            rr = (wz(i  ,j,k,nstp)-wz(i-1,j,k,nstp))
#   ifdef MASKING
     &                                   *umask(i,j)
#   endif
            rrm= wz(i-1,j,k,nstp)
#   ifdef MASKING
     &                                 *umask(i-1,j)
#   endif
            cff1=0.5*(cff*(wz(i-1,j,k,nstp)+wz(i,j,k,nstp))
     &                                        -ABS(cff)*rr)
            Cr=limiter2(cu,WFx(i,j),cff1,rrm,rr,rrp)
            WFx(i,j) = (1-Cr)*cff1 + Cr* WFx(i,j)

          DO i=Istr+1,Iend
            if (k==0) then
              Hzw (i,k) = 0.25 * (HZR(i-1,j,k+1) + HZR(i  ,j,k+1))
            elseif (k==N) then
              Hzw (i,k) = 0.25 * (HZR(i-1,j,k  ) + HZR(i  ,j,k  ))
            else
              Hzw (i,k) = 0.25 * (HZR(i  ,j,k  ) + HZR(i  ,j,k+1)
     &                           +HZR(i-1,j,k  ) + HZR(i-1,j,k+1))
            endif
            cff = Huon_w(i,j,k)
            cu = pn_u(i,j)*cdt*cff
     &           /Hzw(i,k)*pm(i,j)

            rrp= (wz(i+1,j,k,nstp)-wz(i  ,j,k,nstp))
#   ifdef MASKING
     &                                 *umask(i+1,j)
#   endif
            rr = (wz(i  ,j,k,nstp)-wz(i-1,j,k,nstp))
#   ifdef MASKING
     &                                   *umask(i,j)
#   endif
            rrm= (wz(i-1,j,k,nstp)-wz(i-2,j,k,nstp))
#   ifdef MASKING
     &                                 *umask(i-1,j)
#   endif
            cff1=0.5*(cff*(wz(i-1,j,k,nstp)+wz(i,j,k,nstp))
     &                                        -ABS(cff)*rr)
            Cr=limiter2(cu,WFx(i,j),cff1,rrm,rr,rrp)
            WFx(i,j) = (1-Cr)*cff1 + Cr* WFx(i,j)
          ENDDO
          ! wz(i+1,.) out of bounds
          i=Iend+1
            if (k==0) then
              Hzw (i,k) = 0.25 * (HZR(i-1,j,k+1) + HZR(i  ,j,k+1))
            elseif (k==N) then
              Hzw (i,k) = 0.25 * (HZR(i-1,j,k  ) + HZR(i  ,j,k  ))
            else
              Hzw (i,k) = 0.25 * (HZR(i  ,j,k  ) + HZR(i  ,j,k+1)
     &                           +HZR(i-1,j,k  ) + HZR(i-1,j,k+1))
            endif
            cff = Huon_w(i,j,k)
            cu = pn_u(i,j)*cdt*cff
     &           /Hzw(i,k)*pm(i,j)

            rrp= wz(i  ,j,k,nstp)
#   ifdef MASKING
     &                                 *0.0 ! *umask(i+1,j) out of bounds
#   endif
            rr = (wz(i  ,j,k,nstp)-wz(i-1,j,k,nstp))
#   ifdef MASKING
     &                                   *umask(i,j)
#   endif
            rrm= (wz(i-1,j,k,nstp)-wz(i-2,j,k,nstp))
#   ifdef MASKING
     &                                 *umask(i-1,j)
#   endif
            cff1=0.5*(cff*(wz(i-1,j,k,nstp)+wz(i,j,k,nstp))
     &                                        -ABS(cff)*rr)
            Cr=limiter2(cu,WFx(i,j),cff1,rrm,rr,rrp)
            WFx(i,j) = (1-Cr)*cff1 + Cr* WFx(i,j)

        ENDDO
!---
        DO j=jstr,jend+1
          DO i=istr,iend
            WFe(i,j)=0.5*( Wz(i,j,k,nrhs)+Wz(i,j-1,k,nrhs)
     &                                     )*Hvom_w(i,j,k)
          ENDDO
        ENDDO 

        ! (wz(:,j-2)outof bounds
        j=JstrV-1
          DO i=Istr,Iend
            if (k==0) then
              Hzw(i,k) = 0.25 * (HZR(i,j,k+1) + HZR(i,j-1,k+1))
            elseif (k==N) then
              Hzw(i,k) = 0.25 * (HZR(i,j,k  ) + HZR(i,j-1,k  ))
            else
              Hzw(i,k) = 0.25 * (HZR(i,j  ,k) + HZR(i,j  ,k+1)
     &                         + HZR(i,j-1,k) + HZR(i,j-1,k+1))
            endif
            cff = Hvom_w(i,j,k)
            cu = pm_v(i,j)*cdt*cff
     &           /Hzw(i,k)*pn(i,j)

            rrp= (wz(i,j+1,k,nstp)-wz(i,j  ,k,nstp))
#   ifdef MASKING
     &                                 *vmask(i,j+1)
#   endif
            rr = (wz(i,j  ,k,nstp)-wz(i,j-1,k,nstp))
#   ifdef MASKING
     &                                   *vmask(i,j)
#   endif
            rrm= wz(i,j-1,k,nstp)
#   ifdef MASKING
     &                                 *vmask(i,j-1)
#   endif
            cff1=0.5*(cff*(wz(i,j,k,nstp)+wz(i,j-1,k,nstp))
     &                                        -ABS(cff)*rr)
            Cr=limiter2(cu,WFe(i,j),cff1,rrm,rr,rrp)
            WFe(i,j) = (1-Cr)*cff1  + Cr* WFe(i,j)
          ENDDO

        DO j=JstrV,Jend+1
          DO i=Istr,Iend
            if (k==0) then
              Hzw(i,k) = 0.25 * (HZR(i,j,k+1) + HZR(i,j-1,k+1))
            elseif (k==N) then
              Hzw(i,k) = 0.25 * (HZR(i,j,k  ) + HZR(i,j-1,k  ))
            else
              Hzw(i,k) = 0.25 * (HZR(i,j  ,k) + HZR(i,j  ,k+1)
     &                         + HZR(i,j-1,k) + HZR(i,j-1,k+1))
            endif
            cff = Hvom_w(i,j,k)
            cu = pm_v(i,j)*cdt*cff
     &           /Hzw(i,k)*pn(i,j)

            rrp= (wz(i,j+1,k,nstp)-wz(i,j  ,k,nstp))
#   ifdef MASKING
     &                                 *vmask(i,j+1)
#   endif
            rr = (wz(i,j  ,k,nstp)-wz(i,j-1,k,nstp))
#   ifdef MASKING
     &                                   *vmask(i,j)
#   endif
            rrm= (wz(i,j-1,k,nstp)-wz(i,j-2,k,nstp))
#   ifdef MASKING
     &                                 *vmask(i,j-1)
#   endif
            cff1=0.5*(cff*(wz(i,j,k,nstp)+wz(i,j-1,k,nstp))
     &                                        -ABS(cff)*rr)
            Cr=limiter2(cu,WFe(i,j),cff1,rrm,rr,rrp)
            WFe(i,j) = (1-Cr)*cff1  + Cr* WFe(i,j)
          ENDDO
        ENDDO

#  elif defined W_HADV_WENO5 || defined W_HADV_C6 || defined W_HADV_UP5
!
! === WENO5 or 5th/6th order horizontal advection scheme ===
!
#   ifdef W_HADV_WENO5
#    define FLUX5 flux5_weno
#    define FLUX3 flux3_weno
#    define FLUX2 flux1
#    define UP5_MASKING
          cdif=1.
#    include "w_hadv_order5.h"
#    undef FLUX5
#    undef FLUX3
#    undef FLUX2
#    undef UP5_MASKING

#   else /* C6 or UP5 */

        if (nnew.eq.3) then  ! predictor
#    define FLUX5 flux6
#    define FLUX3 flux4
#    define FLUX2 flux2
#    undef  UP5_MASKING
          cdif=0. 
#    include "w_hadv_order5.h"
#    undef FLUX5
#    undef FLUX3
#    undef FLUX2
        else                 ! corrector
#    ifdef W_HADV_C6
#     define FLUX5 flux6
#     define FLUX3 flux4
#     define FLUX2 flux2
#     undef  UP5_MASKING
          cdif=0.
#    else
#     define FLUX5 flux5
#     define FLUX3 flux3
#     define FLUX2 flux1
#     define UP5_MASKING
          cdif=0.5
#    endif
#    include "w_hadv_order5.h"
#    undef FLUX5
#    undef FLUX3
#    undef FLUX2
#    undef UP5_MASKING
        endif
#   endif /* WENO5 */

#  else /* W_HADV_C4 */
!
! === Fourth order horizontal advection scheme ===
!
        DO j=jstr,jend
          DO i=imin,imax+1
            WFx(i,j) = ( wz(i,j,k,nrhs) - wz(i-1,j,k,nrhs) )
#   ifdef MASKING
     &                                           *umask(i,j)
#   endif
          ENDDO
        ENDDO  

#   ifndef EW_PERIODIC
        IF (WESTERN_EDGE) THEN
          DO j=jstr,jend
            WFx(istr-1,j) = WFx(istr  ,j)
          ENDDO
        ENDIF
        IF (EASTERN_EDGE) then
          DO j=jstr,jend
            WFx(iend+2,j) = WFx(iend+1,j)
          ENDDO
        ENDIF
#   endif

        DO j=Jstr,Jend
          DO i=Istr-1,Iend+1          
            lap(i,j) = WFx(i+1,j) - WFx(i,j)
          ENDDO
        ENDDO  
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            IF ( Huon_w(i,j,k).gt. 0. ) THEN
              cff = lap( i-1, j )
            ELSE
              cff = lap( i  , j )
            ENDIF            
            WFx(i,j)=0.5*( wz(i,j,k,nrhs)+wz(i-1,j,k,nrhs)
     &                           -0.25*cff )*Huon_w(i,j,k)
          ENDDO    
        ENDDO
        DO j=jmin,jmax+1
          DO i=istr,iend
            WFe(i,j) = ( wz(i,j,k,nrhs)-wz(i,j-1,k,nrhs) ) 
#   ifdef MASKING
     &                                         *vmask(i,j)
#   endif
          ENDDO 
        ENDDO         
#   ifndef NS_PERIODIC
        IF (SOUTHERN_EDGE) then
          DO i=istr,iend
            WFe(i,jstr-1) = WFe(i,jstr)
          ENDDO
        ENDIF
        IF (NORTHERN_EDGE) then
          DO i=istr,iend
            WFe(i,jend+2) = WFe(i,jend+1)
          ENDDO
        ENDIF
#   endif         
 
        DO j=jstr-1,jend+1
          DO i=istr,iend
            lap(i,j) = WFe(i,j+1) - WFe(i,j  )
          ENDDO
        ENDDO  
        DO j=jstr,jend+1
          DO i=istr,iend
            IF ( HVom_w(i,j,k).gt. 0. ) THEN
               cff = lap( i , j-1 )
            ELSE
               cff = lap( i , j   )
            ENDIF 
            WFe(i,j)=0.5*( wz(i,j,k,nrhs)+wz(i,j-1,k,nrhs)
     &                           -0.25*cff )*Hvom_w(i,j,k)
          ENDDO
        ENDDO 
#  endif /* W_HADV_C2 */

!
!-------------------------------------------------------
! Finalize horizontal advection: compute flux divergence
!-------------------------------------------------------
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            rw(i,j,k)= rw(i,j,k) - WFx(i+1,j) + WFx(i,j)
     &                           - WFe(i,j+1) + WFe(i,j)    
          ENDDO
        ENDDO

      ENDDO          !<-- outer loop k

# endif /* W_HADV */


# if defined UV_ADV && defined W_VADV
!
!=======================================================================
!
!  Vertical advection
!
!=======================================================================
!
      DO j=Jstr,Jend 
 
        DO i=Istr,Iend          
          Hzw (i,0  ) = 0.5 * HZR(i,j,1) 
          Hzw (i,N  ) = 0.5 * HZR(i,j,N)  
          We_r(i,0  ) =        We(i,j,0)
          We_r(i,N  ) = 0.5 * (We(i,j,N)+We(i,j,N-1))             
          We_r(i,N+1) =        We(i,j,N)        
        ENDDO
        DO k=1,N-1
          DO i=Istr,Iend
            Hzw (i,k) = 0.5 * ( HZR(i,j,k) + HZR(i,j,k+1) )
            We_r(i,k) = 0.5 * ( We (i,j,k) + We (i,j,k-1) ) 
          ENDDO        
        ENDDO

#  ifdef W_VADV_SPLINES          
!
! === SPLINES vertical advection scheme ===
!
        DO i=Istr,Iend
          cff     =   0.5 / ( Hzw (i,1) + Hzw (i,0) )       !<- 1 / b(1)
          CF(i,1) =   cff * Hzw(i,0)                        !<- q(1) = c(1) / b(1)
          FC(i,1) =   cff * 3.*( Hzw(i,0)*wz(i,j,1,nrhs)    !<- f(1) / b(1)
     &                         + Hzw(i,1)*wz(i,j,0,nrhs) )
        ENDDO          
        DO k=2,N-1
          DO i=Istr,Iend          
            cff     = 1./( 2.*Hzw(i,k-1)+Hzw(i,k)*(2.-CF(i,k-1)) ) !<- p = 1 / ( b(k)+a(k)*q(k-1) )
            CF(i,k) = cff * Hzw(i,k-1)                             !<- c(k) * p
            FC(i,k) = cff * ( 3.*( Hzw(i,k-1)*wz(i,j,k  ,nrhs)     !<- f(k)=( f(k)-a(k)*f(k-1) )*p
     &                           + Hzw(i,k  )*wz(i,j,k-1,nrhs) )
     &                           - Hzw(i,k  )*FC(i,k-1)  )
          ENDDO
        ENDDO                   
        DO i=Istr,Iend           
          cff     = 1./( 2.*Hzw(i,N-1)+Hzw(i,N)*(2.-CF(i,N-1)) )   !<- p = 1 / ( b(N)+a(N)*q(k-1) )
          FC(i,N) = cff*( 3.*( Hzw(i,N-1)*wz(i,j,N  ,nrhs)         !<- f(N)=( f(N)-a(N)*f(N-1) )*p
     &                       + Hzw(i,N  )*wz(i,j,N-1,nrhs) )
     &                       - Hzw(i,N  )*FC(i,N-1)  )    
        ENDDO          
        DO k=N-1,1,-1
          DO i=Istr,Iend 
            FC(i,k  )=FC(i,k)-CF(i,k)*FC(i,k+1)
            FC(i,k+1)=FC(i,k+1)*We_r(i,k+1)  
          ENDDO
        ENDDO
        DO i=Istr,Iend 
          FC(i,1  )=FC(i,1  )*We_r(i,1)
          FC(i,0  )=0. 
          FC(i,N+1)=0.
        ENDDO

#  elif defined W_VADV_WENO5 || defined W_VADV_C6 || defined W_VADV_UP5
!
! === WENO5 vertical advection scheme ===
!
#   ifdef W_VADV_WENO5
#    define FLUX5 flux5_weno
#    define FLUX3 flux3_weno
#    define FLUX2 flux1
          cdif=1.
#    include "w_vadv_order5.h"
#    undef FLUX5
#    undef FLUX3
#    undef FLUX2
#   else /* C6 or UP5 */
        if (nnew.eq.3) then  ! predictor
#    define FLUX5 flux6
#    define FLUX3 flux4
#    define FLUX2 flux2
          cdif=0.
#    include "w_vadv_order5.h"
#    undef FLUX5
#    undef FLUX3
#    undef FLUX2
        else                 ! corrector
#    ifdef W_VADV_C6
#     define FLUX5 flux6
#     define FLUX3 flux4
#     define FLUX2 flux2
          cdif=0.
#    else
#     define FLUX5 flux5
#     define FLUX3 flux3
#     define FLUX2 flux1
          cdif=0.5
#    endif
#    include "w_vadv_order5.h"
#    undef FLUX5
#    undef FLUX3
#    undef FLUX2
        endif
#   endif /* WENO5 */

#  else
!
! === C2 vertical advection scheme with TVD correction ===
!
        DO k=2,N-1
          DO i=Istr,Iend
            FC(i,k) = 0.5*(wz(i,j,k,nrhs) + wz(i,j,k-1,nrhs))*We_r(i,k)
#   ifdef W_VADV_TVD
            CF(i,k)=(We_r(i,k) * (wz(i,j,k,nrhs) + wz(i,j,k-1,nrhs) )
     &          -abs(We_r(i,k))* (wz(i,j,k,nrhs) - wz(i,j,k-1,nrhs) ) 
     &               ) *0.5
#   endif
          ENDDO
        ENDDO

        DO i=Istr,Iend        
          FC(i,1) = 0.5*( wz(i,j,1,nrhs) + wz(i,j,0,nrhs)  )*We_r(i,1)
#   ifdef W_VADV_TVD
          CF(i,1)=(We_r(i,1) * (wz(i,j,1,nstp) + wz(i,j,0,nstp) )
     &        -abs(We_r(i,1))* (wz(i,j,1,nstp) - wz(i,j,0,nstp) ) )*0.5
#   endif
          FC(i,N) = 0.5*( wz(i,j,N,nrhs) + wz(i,j,N-1,nrhs))*We_r(i,N)
#   ifdef W_VADV_TVD
          CF(i,N)=(We_r(i,N) * (wz(i,j,N,nstp) + wz(i,j,N-1,nstp) )
     &        -abs(We_r(i,N))* (wz(i,j,N,nstp) - wz(i,j,N-1,nstp) ) )*0.5
#   endif
        ENDDO          

        DO i=Istr,Iend
          FC(i,0  )=0.
          FC(i,N+1)=0.
        ENDDO          

#   ifdef W_VADV_TVD
        DO k=1,N
          kp2=min(N+1,k+2)
          kp1=min(N+1,k+1)
          km1=max(1,k-1)
          DO i=IstrU,Iend
            cff= We_r(i,k)
            cu = cff*pm(i,j)*pn(i,j)
     &            *cdt/Hzw(i,k)
            rrp= (wz(i,j,min(kp2,N),nstp) -wz(i,j,min(kp1,N),nstp))
            rr = (wz(i,j,min(kp1,N),nstp) -wz(i ,j,k,nstp))
            rrm= (wz(i,j,k  ,nstp) -wz(i ,j,km1,nstp))
            Cr = limiter2(cu,FC(i,k),cff1,rrm,rr,rrp)
            FC(i,k) = (1-Cr)*CF(i,k) + Cr*FC(i,k)
          ENDDO
        ENDDO
#   endif

#  endif /* W_VADV_SPLINES */
!
!-----------------------------------------------------
! Finalize vertical advection: compute flux divergence
!-----------------------------------------------------
! 
        DO k=1,N
          DO i=Istr,Iend
            rw(i,j,k) = rw(i,j,k)  - FC(i,k+1) + FC(i,k)
#  ifdef MASKING              
            rw(i,j,k) = rw(i,j,k) * rmask(i,j)
#  endif             
          ENDDO
        ENDDO

      ENDDO  !<-- j loop
                     
# endif /* UV_ADV && W_VADV */

      return
      end
!
!=======================================================================
!
# if defined W_VADV_TVD || defined W_HADV_TVD

      function limiter2(nu,hf,lf,rrm,rr,rrp)

      implicit none
      real    :: cff,limiter2,Rj,nu,hf,lf,rrm,rr,rrp
!
!  Upwind        Limiter(Rj)=0.
!  Lax-Wendroff  Limiter(Rj)=1.
!  Min-Mod       Limiter(Rj)=max(0.,min(1.,Rj))
!  Suberbee      Limiter(Rj)=max(0.,max(min(1.,2*Rj),min(2.,Rj)))
!  Van Leer      Limiter(Rj)=(Rj+abs(Rj))/(1+abs(Rj))  => default
!
      if (abs(rr).gt.1.e-20) THEN
        if (nu .gt. 1.e-20) THEN
          Rj=rrm/rr
        else
          Rj=rrp/rr
        endif
      else
        if (nu.gt.1.e-20) then
          Rj=Rrm*1.E20
        else
          Rj=Rrp*1.E20
        endif
      endif

#  ifdef SUPERBEE
       limiter2=max(0.,max(min(1.,2.*Rj),min(2.,Rj)))
#  elif defined MINMOD
       limiter2=max(0.D0,min(1.,Rj))
#  else
       limiter2=(Rj+abs(Rj))/(1+abs(Rj))  ! Van Leer
#  endif

      end function

# endif /* W_VADV_TVD || W_VADV_TVD */
!
!=======================================================================
!
#else
      subroutine rhs3d_w_empty
      end
#endif /* SOLVE3D && NBQ */       
