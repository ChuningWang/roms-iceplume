! $Id: step3d_uv2.F 1480 2014-02-18 14:56:20Z rblod $
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
#ifdef SOLVE3D

      subroutine step3d_uv2 (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()

      call step3d_uv2_tile (Istr,Iend,Jstr,Jend, A2d(1,1,trd),
     &                                  A2d(1,2,trd), A2d(1,3,trd),
     &                                                A2d(1,4,trd)
#ifdef VADV_ADAPT_IMP
     &                                               ,A2d(1,5,trd)
#endif
     &                                                            )
      return
      end

      subroutine step3d_uv2_tile (Istr,Iend,Jstr,Jend, BC,CF,FC,DC
#ifdef VADV_ADAPT_IMP
     &                                                         ,WC
#endif
     &                                                            )
!
!--------------------------------------------------------------------
!  This routine  performs the timestep for the horizontal momentum
! equations. Due to the implicit treatment of vertical viscosity, it
! solves the following tri-diagonal problem:
!
!          FC(k-1) * F(k-1) + BC(k) * F(k) + FC(k) * F(k+1) = DC(k)
!
! where
!          BC(k) = Hz(k) - FC(k) - FC(k-1)
! and
!          FC(k) = - dt * Ak(k) / Hz(k)
! and
!          F(k) is the unknown field: u,v at time step nnew.
!
!  As long as all diffusivity/viscosity coefficients Ak(k) are
! nonnegative, the tri-diagonal matrics is diagonally dominant,
! so that a simple Gaussian elimination algorithm is stable,
! (e.g., Richtmeyer annd  Morton,  1967).
! The implicit vertical viscosity/diffusion terms are discretized
! using implicit backward time step.
!--------------------------------------------------------------------
!
      implicit none
# include "param.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "coupling.h"
# include "forces.h"
# include "mixing.h"
# include "scalars.h"
# include "sources.h"
# if defined M3NUDGING && defined M3CLIMATOLOGY
#  include "climat.h"
# endif
# ifdef DIAGNOSTICS_UV
#  include "diagnostics.h"
# endif
# if defined DIAGNOSTICS_VRT
#  include "diags_vrt.h"
# endif
# ifdef DIAGNOSTICS_EK
#  include "diags_ek.h"
# endif
# ifdef DIAGNOSTICS_PV
#  include "diags_pv.h"
# endif
      integer Istr,Iend,Jstr,Jend, i,j,k
# ifdef PSOURCE
     &       ,is
# endif
      real BC(PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &     CF(PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &     FC(PRIVATE_1D_SCRATCH_ARRAY,0:N), cff,
     &     DC(PRIVATE_1D_SCRATCH_ARRAY,0:N)
# ifdef DIAGNOSTICS_PV
      real cff1 , cff2
# endif
# ifdef VADV_ADAPT_IMP
      real WC(PRIVATE_1D_SCRATCH_ARRAY,0:N)
# endif
      real grad(PRIVATE_2D_SCRATCH_ARRAY)
      real dpth,aa,cc
# ifdef TS_MIX_ISO
      real dRx,dRe, cfilt
# endif
!
# include "compute_auxiliary_bounds.h"
!
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
# ifdef TS_MIX_ISO_FILT
      if (FIRST_TIME_STEP) then
       cfilt=1.                 ! exponential
      else                      ! time smoothing
       cfilt=csmooth            ! coefficient
      endif
# endif
!
!--------------------------------------------------------------------
! Time step momentum equation in the XI-direction.
!
! Compute coefficients FC [dt*Akv/Hz] for the implicit
! vertical viscosity term at future time level at horizontal
! U-points and vertical W-points.
! Load the right-hand-side terms into matrix DC;
! Also compute the main diagonal matrix coeddicient BC.
!--------------------------------------------------------------------
!

!###########################
# if defined DIAGNOSTICS_DISS
        do k=1,N
          do j=JstrR,JendR
                Mrhs(Istr,j,k,1) =  0.
                Mrhs(IendR,j,k,1)=  0.
          enddo
          do i=Istr,IendR
                Mrhs(i,JstrR,k,1) =  0.
                Mrhs(i,JendR,k,1)=  0.
          enddo
          do j=Jstr,JendR
                Mrhs(IstrR,j,k,2) =  0.
                Mrhs(IendR,j,k,2)=  0.
          enddo
          do i=IstrR,IendR
                Mrhs(i,Jstr,k,2) =  0.
                Mrhs(i,JendR,k,2)=  0.
          enddo
        enddo
# endif
!###########################


      do j=Jstr,Jend
!
! The vertical viscosity term is not added as an explicit divergence term
! like the other terms due to the semi-implicit backward Euler scheme.
! We store first u(i,j,k,nnew) ...
!
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL || defined DIAGNOSTICS_PV
        do k=1,N
          do i=IstrU,Iend
             MVmix(i,j,k,1)=u(i,j,k,nnew)
          enddo
        enddo
# elif defined DIAGNOSTICS_EK
        do i=IstrU,Iend
          k=1
             ekwrkVmix(i,j,1) =
     &           - u(i,j,k,nnew)
     &           *u(i,j,k,nrhs)
          do k=2,N
             ekwrkVmix(i,j,1) = ekwrkVmix(i,j,1) -
     &           u(i,j,k,nnew)
     &           *u(i,j,k,nrhs)
          enddo
# if defined DIAGNOSTICS_EK_MLD
             ekwrkVmix_mld(i,j,1) = 0.
          do k=kbl(i,j),N
             ekwrkVmix_mld(i,j,1) = ekwrkVmix_mld(i,j,1) -
     &           u(i,j,k,nnew)
     &           *u(i,j,k,nrhs)
          enddo
# endif
        enddo
# endif
# if defined DIAGNOSTICS_VRT && ! defined DIAGNOSTICS_UV
        do i=IstrU,Iend
          k=1
             wrkVmix(i,j,1) =
     &           - u(i,j,k,nnew)
          do k=2,N
             wrkVmix(i,j,1) = wrkVmix(i,j,1) -
     &           u(i,j,k,nnew)
          enddo
        enddo
# endif

!--------------------------------------------------------------------
! Resolve the tri-diagonal system. Also perform coupling between the
! barotropic and baroclinic modes. For the purpose of the processor
! efficiency these two operations are partially overlapped. So that
! in the third and fourth loop below operations above the separator
! (!-) belong to the tri-diagonal solver, while operations below it
! belong to the coupling.
!
! The coupling procedure implies vertical integration of the newly
! computed XI-component velocity [which ends up being stored in the
! scratch array DC(:,1:N) after the tri-diagonal problem is resolved,
! the fourth loop below] in order to obtain its vertical mean [stored
! as DC(:,0)], and replace it with the more accurate mean obtained
! from the 2D barotropic sub-model ubar=DU_avg1/(D*on_u), where
! D=CF(:,0) is the total depth of water column.
!
!--------------------------------------------------------------------
        do i=IstrU,Iend
          FC(i,0)=0.
          FC(i,N)=0.
# ifdef VADV_ADAPT_IMP
          DC(i,0)=0.25*dt*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
          WC(i,0)=0.
          WC(i,N)=0.
# endif
        enddo

        do k=1,N-1
          do i=IstrU,Iend
            FC(i,k)=-dt*(Akv(i,j,k)+Akv(i-1,j,k))
     &                 /( z_r(i,j,k+1)+z_r(i-1,j,k+1)
     &                   -z_r(i,j,k  )-z_r(i-1,j,k  ))
# ifdef VADV_ADAPT_IMP
            WC(i,k)=DC(i,0)*0.5*(Wi(i,j,k)+Wi(i-1,j,k))
# endif
          enddo
        enddo

        do k=1,N
          do i=IstrU,Iend
            DC(i,k)=u(i,j,k,nnew)
            BC(i,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))-FC(i,k)-FC(i,k-1)
# ifdef VADV_ADAPT_IMP
     &              +max(WC(i,k  ),0.)-min(WC(i,k-1),0.)
# endif
          enddo
        enddo      !--> discard DC(:,0); keep FC
!
! Apply wind forcing.
!
# ifndef BODYFORCE
        do i=IstrU,Iend
          DC(i,N)=DC(i,N)+dt*sustr(i,j)
#  ifndef BSTRESS_FAST
          DC(i,1)=DC(i,1)-dt*bustr(i,j)
#  endif
#  ifdef DIAGNOSTICS_VRT
            wrkDrag(i,j,1)= -bustr(i,j)*om_u(i,j)
     &                    SWITCH umask(i,j)
            wrkWind(i,j,1)= sustr(i,j)*om_u(i,j)
     &                    SWITCH umask(i,j)
#  endif /* DIAGNOSTICS_VRT */
#  ifdef DIAGNOSTICS_EK
          ekwrkDrag(i,j,1)= -bustr(i,j)*om_u(i,j)*on_u(i,j)
     &                    * u(i,j,1,nrhs) SWITCH umask(i,j)
          ekwrkWind(i,j,1)= sustr(i,j)*om_u(i,j)*on_u(i,j)
     &                    * u(i,j,N,nrhs) SWITCH umask(i,j)
#  endif /* DIAGNOSTICS_EK */
        enddo
# endif /* BODYFORCE */
        do i=IstrU,Iend
          cff=1./BC(i,1)
          CF(i,1)=cff*( FC(i,1)                  !<-- q(1) = c(1)/b(1)
# ifdef VADV_ADAPT_IMP
     &                +min( WC(i,1),0.)
# endif
     &                                 )
          DC(i,1)=cff*DC(i,1)
        enddo
        do k=2,N-1
          do i=IstrU,Iend
# ifdef VADV_ADAPT_IMP
              aa = FC(i,k-1)-max(WC(i,k-1),0.)   !<-- a = sub   (lower) diagonal
              cc = FC(i,k  )+min(WC(i,k  ),0.)   !<-- c = super (upper) diagonal
# else
              aa = FC(i,k-1)
              cc = FC(i,k  )
# endif
            cff=1./(BC(i,k)-aa*CF(i,k-1))        !<-- p    = 1/(b(k)-a(k)*q(k))
            CF(i,k)=cff*cc
            DC(i,k)=cff*(DC(i,k)-aa*DC(i,k-1))
          enddo
        enddo
        do i=IstrU,Iend
# ifdef VADV_ADAPT_IMP
              aa = FC(i,N-1)-max(WC(i,N-1),0.)   !<-- a = sub   (lower) diagonal
# else
              aa = FC(i,N-1)
# endif
          DC(i,N)=(DC(i,N)-aa*DC(i,N-1))
     &           /(BC(i,N)-aa*CF(i,N-1))

          CF(i,0)=0.5*(Hz(i,j,N)+Hz(i-1,j,N))
          DC(i,0)=CF(i,0)*DC(i,N)
        enddo
        do k=N-1,1,-1
          do i=IstrU,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)

            cff=0.5*(Hz(i,j,k)+Hz(i-1,j,k))
            CF(i,0)=CF(i,0)+cff
            DC(i,0)=DC(i,0)+cff*DC(i,k)
          enddo
        enddo      !--> discard FC,BC; keep DC,CF(:,0)

        do i=IstrU,Iend
          DC(i,0)=(DC(i,0)*on_u(i,j)-DU_avg1(i,j,nnew))
     &                             /(CF(i,0)*on_u(i,j))
        enddo




!###########################
        do k=1,N
          do i=IstrU,Iend
            u(i,j,k,nnew)=(DC(i,k)-DC(i,0)) SWITCH umask(i,j)
          enddo
        enddo        !--> discard DC,CF(:,0)
!###########################

!###########################
# if defined DIAGNOSTICS_DISS
! Compute the transport part of ke_vmix : d/dz[ Akv u du/dz]
! ------- ------------- ------ ------ - -----
      do k=1,N-1
        do i=IstrU,Iend

            FC(i,k)= 0.5*(Akv(i,j,k)+Akv(i-1,j,k))

            cff1= 2.* (2.*MVmix(i,j,k+1,1) / (Hz(i,j,k+1)+Hz(i-1,j,k+1))
     &               - 2.*MVmix(i,j,k,1) / (Hz(i,j,k  )+Hz(i-1,j,k  )) )
     &                 /( z_r(i,j,k+1)+z_r(i-1,j,k+1)
     &                   -z_r(i,j,k  )-z_r(i-1,j,k  ))

            BC(i,k)= 0.25 * (u(i,j,k+1,nrhs) + u(i,j,k+1,nnew)
     &                     + u(i,j,k  ,nrhs) + u(i,j,k  ,nnew)) * cff1

         enddo
       enddo

      do i=IstrU,Iend
          FC(i,0)= 0.
          BC(i,0)= 0.
          FC(i,N)= 1.
          BC(i,N)= 0.
!          BC(i,N)= 0.5 * (u(i,j,k,nrhs)+u(i,j,k,nnew)) * sustr(i,j)
      enddo

      do k=1,N
        do i=IstrU,Iend
            Mrhs(i,j,k,1) =  ( FC(i,k  ) * BC(i,k  )
     &                       - FC(i,k-1) * BC(i,k-1) )
     &                *  om_u(i,j)*on_u(i,j) SWITCH umask(i,j)
         enddo
       enddo
! ------- ------------- ------ ------ - -----

# endif



!###########################
        do k=1,N
          do i=IstrU,Iend
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL || defined DIAGNOSTICS_PV
            MVmix(i,j,k,1)=-( MVmix(i,j,k,1) - u(i,j,k,nnew)*
     &                           0.5*(Hz(i,j,k)+Hz(i-1,j,k)) )
     &                    *om_u(i,j)*on_u(i,j)/dt SWITCH umask(i,j)
# ifdef BODYFORCE
            MVmix(i,j,k,1)= MVmix(i,j,k,1) + Mbody(i,j,k,1)
# endif
# elif defined DIAGNOSTICS_EK
             ekwrkVmix(i,j,1) = ekwrkVmix(i,j,1) +
     &           0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nnew)
     &           * u(i,j,k,nrhs)
# if defined DIAGNOSTICS_EK_MLD
            if (k.ge.kbl(i,j)) then
             ekwrkVmix_mld(i,j,1) = ekwrkVmix_mld(i,j,1) +
     &           0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nnew)
     &           * u(i,j,k,nrhs)
            endif
# endif
# endif
# if defined DIAGNOSTICS_VRT && ! defined DIAGNOSTICS_UV
            wrkVmix(i,j,1) = wrkVmix(i,j,1) +
     &           0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nnew)
# endif
          enddo
        enddo        !--> discard DC,CF(:,0)


# if defined DIAGNOSTICS_VRT && !defined DIAGNOSTICS_UV
        do i=IstrU,Iend
          wrkVmix(i,j,1) = wrkVmix(i,j,1)
     &               *om_u(i,j)*on_u(i,j)/dt SWITCH umask(i,j)
# ifdef BODYFORCE
          wrkVmix(i,j,1) = wrkVmix(i,j,1)
     &               + (wrkDrag(i,j,1) + wrkWind(i,j,1))
     &               * on_u(i,j)
# endif
        enddo
# endif
# if defined DIAGNOSTICS_EK && ! defined DIAGNOSTICS_EK_FULL
        do i=IstrU,Iend
             ekwrkVmix(i,j,1) = ekwrkVmix(i,j,1)
     &               *om_u(i,j)*on_u(i,j)/dt SWITCH umask(i,j)
# ifdef BODYFORCE
             ekwrkVmix(i,j,1) = ekwrkVmix(i,j,1)
     &               + ekwrkDrag(i,j,1) +ekwrkWind(i,j,1)
# endif
# if defined DIAGNOSTICS_EK_MLD
             ekwrkVmix_mld(i,j,1) = ekwrkVmix_mld(i,j,1)
     &               *om_u(i,j)*on_u(i,j)/dt SWITCH umask(i,j)
# endif
        enddo
# endif
!###########################

# if defined DIAGNOSTICS_BARO
!
! Now add only coupling contribution to the pressure gradient diagnostics
!
        do k=1,N
          do i=IstrU,Iend
!            u(i,j,k,nnew)=(DC(i,k)-DC(i,0)) SWITCH umask(i,j)
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL || defined DIAGNOSTICS_PV
            MBaro(i,j,k,1)= -DC(i,0)
     &              *0.5*(Hz(i,j,k)+Hz(i-1,j,k))
     &              *om_u(i,j)*on_u(i,j)/dt SWITCH umask(i,j)
# elif defined DIAGNOSTICS_EK
             ekwrkBaro(i,j,1) = ekwrkBaro(i,j,1) -
     &           0.5*(Hz(i,j,k)+Hz(i-1,j,k))*DC(i,0)
     &           * u(i,j,k,nrhs)
     &           * om_u(i,j)*on_u(i,j)/dt SWITCH umask(i,j)
# if defined DIAGNOSTICS_EK_MLD
            if (k.ge.kbl(i,j)) then
             ekwrkBaro_mld(i,j,1) = ekwrkBaro_mld(i,j,1)
     &           - 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*DC(i,0)
     &           * u(i,j,k,nrhs)
     &           *om_u(i,j)*on_u(i,j)/dt SWITCH umask(i,j)
            endif
# endif
# endif
# if defined DIAGNOSTICS_VRT && ! defined DIAGNOSTICS_UV
            wrkBaro(i,j,1) = wrkBaro(i,j,1)
     &           - 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*DC(i,0)
     &              *om_u(i,j)*on_u(i,j)/dt SWITCH umask(i,j)
# endif
          enddo
        enddo        !--> discard DC,CF(:,0)
# endif




!
!--------------------------------------------------------------------
! Time step momentum equation in the ETA-direction.

! Compute off-diagonal matrix coefficients FC [dt*Akv/Hz]
! for the implicit vertical viscosity term at the new time step,
! at horizontal V-points and vertical W-points.
! Set vertical boundary conditions.
! After that load the right-hand-side terms into matrix DC;
! Also compute the main diagonal matrix coefficient BC.
!--------------------------------------------------------------------
!
        if (j.ge.JstrV) then
!
! The vertical viscosity term is not added as an explicit divergence term
! like the other terms due to the semi-implicit backward Euler scheme.
! We store first v(i,j,k,nnew) ...
!
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL || defined DIAGNOSTICS_PV
        do k=1,N
          do i=Istr,Iend
             MVmix(i,j,k,2)=v(i,j,k,nnew)
          enddo
        enddo
# elif defined DIAGNOSTICS_EK
        do i=Istr,Iend
          k=1
             ekwrkVmix(i,j,2) =
     &           - v(i,j,k,nnew)
     &           *v(i,j,k,nrhs)
          do k=2,N
             ekwrkVmix(i,j,2) = ekwrkVmix(i,j,2) -
     &           v(i,j,k,nnew)
     &           *v(i,j,k,nrhs)
          enddo
# if defined DIAGNOSTICS_EK_MLD
             ekwrkVmix_mld(i,j,2) = 0.
          do k=kbl(i,j),N
             ekwrkVmix_mld(i,j,2) = ekwrkVmix_mld(i,j,2) -
     &           v(i,j,k,nnew)
     &           *v(i,j,k,nrhs)
          enddo
# endif
        enddo
# endif
# if defined DIAGNOSTICS_VRT && ! defined DIAGNOSTICS_UV
        do i=Istr,Iend
          k=1
             wrkVmix(i,j,2) =
     &           - v(i,j,k,nnew)
          do k=2,N
             wrkVmix(i,j,2) = wrkVmix(i,j,2) -
     &           v(i,j,k,nnew)
          enddo
        enddo
# endif
          do i=Istr,Iend
            FC(i,0)=0.
            FC(i,N)=0.
# ifdef VADV_ADAPT_IMP
            DC(i,0)=0.25*dt*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
            WC(i,0)=0.
            WC(i,N)=0.
# endif
          enddo

          do k=1,N-1
            do i=Istr,Iend
              FC(i,k)=-dt*(Akv(i,j,k)+Akv(i,j-1,k))
     &                   /( z_r(i,j,k+1)+z_r(i,j-1,k+1)
     &                     -z_r(i,j,k  )-z_r(i,j-1,k  ))
# ifdef VADV_ADAPT_IMP
              WC(i,k)=DC(i,0)*0.5*(Wi(i,j,k)+Wi(i,j-1,k))
# endif
            enddo
          enddo

          do k=1,N
            do i=Istr,Iend
              DC(i,k)=v(i,j,k,nnew)
              BC(i,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))-FC(i,k)-FC(i,k-1)
# ifdef VADV_ADAPT_IMP
     &              +max(WC(i,k  ),0.)-min(WC(i,k-1),0.)
# endif
            enddo
          enddo
!
! Apply wind forcing
!
# ifndef BODYFORCE
          do i=Istr,Iend
            DC(i,N)=DC(i,N)+dt*svstr(i,j)
#  ifndef BSTRESS_FAST
            DC(i,1)=DC(i,1)-dt*bvstr(i,j)
#  endif
#  ifdef DIAGNOSTICS_VRT
            wrkDrag(i,j,2)= -bvstr(i,j)*on_v(i,j)
     &                    SWITCH vmask(i,j)
            wrkWind(i,j,2)= svstr(i,j)*on_v(i,j)
     &                    SWITCH vmask(i,j)
#  endif /* DIAGNOSTICS_VRT */
#  ifdef DIAGNOSTICS_EK
          ekwrkDrag(i,j,2)= -bvstr(i,j)*om_v(i,j)*on_v(i,j)
     &                    * v(i,j,1,nrhs) SWITCH vmask(i,j)
          ekwrkWind(i,j,2)= svstr(i,j)*om_v(i,j)*on_v(i,j)
     &                    * v(i,j,N,nrhs) SWITCH vmask(i,j)
#  endif /* DIAGNOSTICS_EK */
          enddo
# endif /* BODYFORCE */
!
! Resolve the tri-diagonal system. Also perform coupling between the
! barotropic and baroclinic modes (similarly to the XI-component, see
! comment above).
!
          do i=Istr,Iend
            cff=1./BC(i,1)
            CF(i,1)=cff*( FC(i,1)                  !<-- q(1) = c(1)/b(1)
# ifdef VADV_ADAPT_IMP
     &                +min( WC(i,1),0.)
# endif
     &                                 )
            DC(i,1)=cff*DC(i,1)
          enddo
          do k=2,N-1
            do i=Istr,Iend
# ifdef VADV_ADAPT_IMP
              aa = FC(i,k-1)-max(WC(i,k-1),0.)     !<-- a = sub   (lower) diagonal
              cc = FC(i,k  )+min(WC(i,k  ),0.)     !<-- c = super (upper) diagonal
# else
              aa = FC(i,k-1)
              cc = FC(i,k  )
# endif
              cff=1./(BC(i,k)-aa*CF(i,k-1))
              CF(i,k)=cff*cc
              DC(i,k)=cff*(DC(i,k)-aa*DC(i,k-1))
            enddo
          enddo
          do i=Istr,Iend
# ifdef VADV_ADAPT_IMP
              aa = FC(i,N-1)-max(WC(i,N-1),0.)     !<-- a = sub   (lower) diagonal
# else
              aa = FC(i,N-1)
# endif
            DC(i,N)=(DC(i,N)-aa*DC(i,N-1))
     &             /(BC(i,N)-aa*CF(i,N-1))

            CF(i,0)=0.5*(Hz(i,j,N)+Hz(i,j-1,N))
            DC(i,0)=CF(i,0)*DC(i,N)
          enddo
          do k=N-1,1,-1
            do i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)

              cff=0.5*(Hz(i,j,k)+Hz(i,j-1,k))
              CF(i,0)=CF(i,0)+cff
              DC(i,0)=DC(i,0)+cff*DC(i,k)
            enddo
          enddo      !--> discard FC,BC; keep DC,CF(:,0)

          do i=Istr,Iend
            DC(i,0)=(DC(i,0)*om_v(i,j)-DV_avg1(i,j,nnew))
     &                               /(CF(i,0)*om_v(i,j))
          enddo



!###########################
          do k=1,N
            do i=Istr,Iend
              v(i,j,k,nnew)=(DC(i,k)-DC(i,0)) SWITCH vmask(i,j)
          enddo
        enddo        !--> discard DC,CF(:,0)
!###########################

!###########################
# if defined DIAGNOSTICS_DISS
! Compute the transport part of ke_vmix : d/dz[ Akv v dv/dz]

! ------- ------------- ------ ------ - -----
      do k=1,N-1
        do i=Istr,Iend

            FC(i,k)= 0.5*(Akv(i,j,k)+Akv(i,j-1,k))

            cff1= 2.* (2.*MVmix(i,j,k+1,2)/(Hz(i,j,k+1)+Hz(i,j-1,k+1))
     &               - 2.*MVmix(i,j,k  ,2)/(Hz(i,j,k  )+Hz(i,j-1,k  )) )
     &                 /( z_r(i,j,k+1)+z_r(i,j-1,k+1)
     &                   -z_r(i,j,k  )-z_r(i,j-1,k  ))

            BC(i,k)= 0.25 * (v(i,j,k+1,nrhs) + v(i,j,k+1,nnew)
     &                     + v(i,j,k  ,nrhs) + v(i,j,k  ,nnew)) * cff1

         enddo
       enddo

      do i=Istr,Iend
          FC(i,0)= 0.
          BC(i,0)= 0.
          FC(i,N)= 1.
          BC(i,N)= 0.
!          BC(i,N)= 0.5 * (v(i,j,k,nrhs)+v(i,j,k,nnew)) * svstr(i,j)
      enddo

      do k=1,N
        do i=Istr,Iend
            Mrhs(i,j,k,2) =  ( FC(i,k  ) * BC(i,k  )
     &                       - FC(i,k-1) * BC(i,k-1) )
     &                *  om_v(i,j)*on_v(i,j) SWITCH vmask(i,j)
         enddo
       enddo
! ------- ------------- ------ ------ - -----

# endif



!###########################
          do k=1,N
            do i=Istr,Iend
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL || defined DIAGNOSTICS_PV
            MVmix(i,j,k,2)=-( MVmix(i,j,k,2) - v(i,j,k,nnew)*
     &                           0.5*(Hz(i,j,k)+Hz(i,j-1,k)) )
     &                    *om_v(i,j)*on_v(i,j)/dt SWITCH vmask(i,j)
# ifdef BODYFORCE
            MVmix(i,j,k,2)= MVmix(i,j,k,2) + Mbody(i,j,k,2)
# endif
# elif defined DIAGNOSTICS_EK
             ekwrkVmix(i,j,2) = ekwrkVmix(i,j,2) +
     &           0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nnew)
     &           * v(i,j,k,nrhs)
# if defined DIAGNOSTICS_EK_MLD
            if (k.ge.kbl(i,j)) then
             ekwrkVmix_mld(i,j,2) = ekwrkVmix_mld(i,j,2) +
     &           0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nnew)
     &           * v(i,j,k,nrhs)
            endif
# endif
# endif
# if defined DIAGNOSTICS_VRT && ! defined DIAGNOSTICS_UV
            wrkVmix(i,j,2) = wrkVmix(i,j,2) +
     &           0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nnew)
# endif
          enddo
        enddo        !--> discard DC,CF(:,0)
# if (defined DIAGNOSTICS_VRT && !defined DIAGNOSTICS_UV)
        do i=Istr,Iend
          wrkVmix(i,j,2) = wrkVmix(i,j,2)
     &               * om_v(i,j)*on_v(i,j)/dt SWITCH vmask(i,j)
# ifdef BODYFORCE
          wrkVmix(i,j,2) = wrkVmix(i,j,2)
     &               + (wrkDrag(i,j,2) +wrkWind(i,j,2))
     &               * om_v(i,j)
# endif
        enddo
# endif
# if defined DIAGNOSTICS_EK && ! defined DIAGNOSTICS_EK_FULL
        do i=Istr,Iend
             ekwrkVmix(i,j,2) = ekwrkVmix(i,j,2)
     &               *om_v(i,j)*on_v(i,j)/dt SWITCH vmask(i,j)
# ifdef BODYFORCE
             ekwrkVmix(i,j,2) = ekwrkVmix(i,j,2)
     &               + ekwrkDrag(i,j,2) +ekwrkWind(i,j,2)
# endif
# if defined DIAGNOSTICS_EK_MLD
             ekwrkVmix_mld(i,j,2) = ekwrkVmix_mld(i,j,2)
     &               *om_v(i,j)*on_v(i,j)/dt SWITCH vmask(i,j)
# endif
        enddo
# endif
!###########################
# if defined DIAGNOSTICS_BARO
!
! Now add only coupling contribution to the pressure gradient diagnostics
!
          do k=1,N
            do i=Istr,Iend
!              v(i,j,k,nnew)=(DC(i,k)-DC(i,0)) SWITCH vmask(i,j)
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL || defined DIAGNOSTICS_PV
            MBaro(i,j,k,2)= -DC(i,0)*
     &                           0.5*(Hz(i,j,k)+Hz(i,j-1,k))
     &                    *om_v(i,j)*on_v(i,j)/dt SWITCH vmask(i,j)
# elif defined DIAGNOSTICS_EK
             ekwrkBaro(i,j,2) = ekwrkBaro(i,j,2) -
     &           0.5*(Hz(i,j,k)+Hz(i,j-1,k))*DC(i,0)
     &           * v(i,j,k,nrhs)
     &           *om_v(i,j)*on_v(i,j)/dt SWITCH vmask(i,j)
# if defined DIAGNOSTICS_EK_MLD
            if (k.ge.kbl(i,j)) then
             ekwrkBaro_mld(i,j,2) = ekwrkBaro_mld(i,j,2)
     &           - 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*DC(i,0)
     &           * v(i,j,k,nrhs)
     &           *om_v(i,j)*on_v(i,j)/dt SWITCH vmask(i,j)
            endif
# endif
# endif
# if defined DIAGNOSTICS_VRT && ! defined DIAGNOSTICS_UV
            wrkBaro(i,j,2) = wrkBaro(i,j,2)
     &           - 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*DC(i,0)
     &           *om_v(i,j)*on_v(i,j)/dt SWITCH vmask(i,j)
# endif
          enddo
        enddo        !--> discard DC,CF(:,0)
# endif
        endif
      enddo      ! <-- j

!-------------------------------------------------------------------
! DIAGNOSTICS
!-------------------------------------------------------------------
# if defined DIAGNOSTICS_UV  || defined DIAGNOSTICS_PV
      do k=1,N
        do j=Jstr,Jend
          do i=IstrU,Iend
            Mrate(i,j,k,1)=( 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nnew)
     &              -0.5*(Hz_bak(i,j,k)+Hz_bak(i-1,j,k))*u(i,j,k,nstp) )
     &                                           *om_u(i,j)*on_u(i,j)/dt

          enddo
        enddo

        do j=JstrV,Jend
          do i=Istr,Iend
            Mrate(i,j,k,2)=( 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nnew)
     &              -0.5*(Hz_bak(i,j,k)+Hz_bak(i,j-1,k))*v(i,j,k,nstp) )
     &                                           *om_v(i,j)*on_v(i,j)/dt
          enddo
        enddo
      enddo
# endif
# if defined DIAGNOSTICS_EK
        do j=JstrR,JendR
          do i=IstrR,IendR
            ekwrkVmix2(i,j,1)=0
            ekwrkVmix2(i,j,2)=0
# if defined DIAGNOSTICS_EK_MLD
            ekwrkVmix2_mld(i,j,1)=0
            ekwrkVmix2_mld(i,j,2)=0
# endif
          enddo
        enddo
        do j=Jstr,Jend
          do i=IstrU,Iend
            do k=1,N
                ekwrkVmix2(i,j,1) = ekwrkVmix2(i,j,1) + 0.5*
     &            ( 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nnew)**2
     &         -0.5*(Hz_bak(i,j,k)+Hz_bak(i-1,j,k))*u(i,j,k,nstp)**2)
     &                                        *om_u(i,j)*on_u(i,j)/dt
#  ifdef MASKING
     &                     *umask(i,j)
#  endif
            enddo
# if defined DIAGNOSTICS_EK_MLD
            do k=kbl(i,j),N
                ekwrkVmix2_mld(i,j,1) = ekwrkVmix2_mld(i,j,1) + 0.5*
     &            ( 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nnew)**2
     &         -0.5*(Hz_bak(i,j,k)+Hz_bak(i-1,j,k))*u(i,j,k,nstp)**2)
     &                                        *om_u(i,j)*on_u(i,j)/dt
#  ifdef MASKING
     &                     *umask(i,j)
#  endif
            enddo
# endif
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            do k=1,N
              ekwrkVmix2(i,j,2) = ekwrkVmix2(i,j,2) + 0.5*
     &           ( 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nnew)**2
     &         -0.5*(Hz_bak(i,j,k)+Hz_bak(i,j-1,k))*v(i,j,k,nstp)**2)
     &                                        *om_v(i,j)*on_v(i,j)/dt
#  ifdef MASKING
     &                     *vmask(i,j)
#  endif
            enddo
# if defined DIAGNOSTICS_EK_MLD
            do k=kbl(i,j),N
              ekwrkVmix2_mld(i,j,2) = ekwrkVmix2_mld(i,j,2) + 0.5*
     &           ( 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nnew)**2
     &         -0.5*(Hz_bak(i,j,k)+Hz_bak(i,j-1,k))*v(i,j,k,nstp)**2)
     &                                        *om_v(i,j)*on_v(i,j)/dt
#  ifdef MASKING
     &                     *vmask(i,j)
#  endif
            enddo
# endif
          enddo
        enddo
# endif /* DIAGNOSTICS_EK */


!
!--------------------------------------------------------------------
!  Compute wet/dry masks for total velocity
!--------------------------------------------------------------------
!
# ifdef WET_DRY
      CALL wetdry_tile (Istr,Iend,Jstr,Jend)
# endif
!
!--------------------------------------------------------------------
! Set PHYSICAL lateral boundary conditions.
!--------------------------------------------------------------------
!
      call v3dbc_tile (Istr,Iend,Jstr,Jend, grad)
      call u3dbc_tile (Istr,Iend,Jstr,Jend, grad)


!
!
!--------------------------------------------------------------------
! Coupling 2D and 3D momentum equations:
!
! Compute inaccurate vertical mean of the three-dimensional
! velocity field, then subtract it and replace it with the vertically
! integrated (barotropic) velocity field computed from the two-
! dimensional submodel. After that compute mass fluxes through grid
! box faces.
!
! Meaning of scratch variables below:
!
! DC(i,k) [where k=1,N]  height of grid box for U- or V-cell.
! DC(i,0) total depth of water column at horizontal U- or V-points.
! CF(i,0) vertically integrated mass flux/error/correction for
!                       the new time step velocity u,v(:,:,:,nnew)
! FC(i,0) vertically integrated mass flux/error/correction for the
!             intermediate time step [n+1/2] mass fluxes Huon,Hvom.
!
! This procedure also replaces instantaneous (with respect to the
! fast time step) vertically integrated (barotropic) velocities with
! their values based on fast-time-averaged mass fluxes. These are to
! be used as initial conditions for the barotropic mode at the new
! time step.
!
! Explanation of horizontal loop indices: in the case of periodic
! boundaries (in either direction) the coupling and computation of
! mass fluxes Huon,Hvom is performed within the internal range of
! indices (excluding ghost zones), after that the ghost points for
! the newly computed arrays are exchanged; in the case of nonperiodic
! boundaries the coupling is done over the extendeed range of indices
! (that is including boundary points).
!--------------------------------------------------------------------
!
# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif
!
# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
!
!--------------------------------------------------
!  Couple XI-component of velocity
!--------------------------------------------------
!
      do j=JU_RANGE
        do i=IU_RANGE
          DC(i,0)=0.
          CF(i,0)=0.
          FC(i,0)=0.
        enddo
        do k=1,N,+1
          do i=IU_RANGE
            DC(i,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
            CF(i,0)=CF(i,0)+DC(i,k)*u(i,j,k,nnew)
          enddo
        enddo
        do i=IU_RANGE
          DC(i,0)=1./DC(i,0)
          CF(i,0)=DC(i,0)*(CF(i,0)-DU_avg1(i,j,nnew))
# ifdef MRL_WCI
     &                               +ust2d(i,j)
# endif
          ubar(i,j,knew)=DC(i,0)*DU_avg1(i,j,nnew)
# ifdef MRL_WCI
     &                               -ust2d(i,j)
#  ifdef MASKING
          ubar(i,j,knew)=ubar(i,j,knew)*umask(i,j)
     &                +ust2d(i,j)*(umask(i,j)-1.0)
#  endif
# endif
# ifdef WET_DRY
          ubar(i,j,knew)=ubar(i,j,knew)*umask_wet(i,j)
# endif
        enddo
        do k=N,1,-1
          do i=IU_RANGE
            u(i,j,k,nnew)=(u(i,j,k,nnew)-CF(i,0))
# ifdef MASKING
     &                                *umask(i,j)
#  ifdef MRL_WCI
     &               +ust(i,j,k)*(umask(i,j)-1.0)
#  endif
# endif
# ifdef WET_DRY
            u(i,j,k,nnew)=u(i,j,k,nnew)*umask_wet(i,j)
#  ifdef MRL_WCI
            ust(i,j,k)=ust(i,j,k)*umask_wet(i,j)
#  endif
# endif

!--# define TR
# ifdef TR
#  ifdef MRL_WCI
            FC(i,k)=DC(i,k)*(0.5*(u(i,j,k,nstp)+u(i,j,k,nnew))
     &                                            +ust(i,j,k))
#  else
            FC(i,k)=0.5*(u(i,j,k,nstp)+u(i,j,k,nnew))*DC(i,k)
#  endif
# else
#  define EPSIL 0.125
#  define DELTA 0.75
#  ifdef MRL_WCI
            FC(i,k)=DELTA*Huon(i,j,k) +
     &              EPSIL*DC(i,k)*(u(i,j,k,nstp)+u(i,j,k,nnew)
     &                                        +2.0*ust(i,j,k))
#  else
            FC(i,k)=DELTA*Huon(i,j,k) +
     &              EPSIL*DC(i,k)*(u(i,j,k,nstp)+u(i,j,k,nnew))
#  endif
# endif /* TR */

            FC(i,0)=FC(i,0)+FC(i,k)
          enddo
        enddo
        do i=IU_RANGE
          FC(i,0)=DC(i,0)*(FC(i,0)-DU_avg2(i,j))
        enddo
        do k=1,N,+1
          do i=IU_RANGE
            Huon(i,j,k)=FC(i,k)-DC(i,k)*FC(i,0)
          enddo
        enddo
!
!-------------------------------------------------
! Compute density XI gradients for use in t3dmix
!-------------------------------------------------
!
# if !defined TS_MIX_S && \
      (defined TS_DIF4 || defined TS_DIF2 || defined SPONGE_DIF2)
!
        do k=1,N-1,+1
          do i=IU_RANGE
            cff =0.5*(pm(i,j)+pm(i-1,j)) SWITCH umask(i,j)
#  ifdef TS_MIX_ISO
#   ifdef SPLIT_EOS
            dpth=0.5*( z_w(i,j,N)+z_w(i-1,j,N)
     &                -z_r(i,j,k)-z_r(i-1,j,k))
            dRx=cff*( rho1(i,j,k)-rho1(i-1,j,k)
     &             +(  qp1(i,j,k)- qp1(i-1,j,k))
     &                     *dpth*(1.-qp2*dpth) )
#   else
            dRx=cff*(rho1(i,j,k)-rho1(i-1,j,k))
#   endif
#   ifdef TS_MIX_ISO_FILT
            dRx  = cfilt*dRx + (1.-cfilt)*dRdx(i,j,k)  ! exp. time smoothing
#   endif
            dRdx(i,j,k)=dRx
#  elif defined TS_MIX_GEO
            dRdx(i,j,k)=cff*(z_r(i,j,k)-z_r(i-1,j,k))  ! dRdx get dZdx values
#  else
            dRdx(i,j,k)=0.
#  endif
          enddo
        enddo
        do i=IU_RANGE
          dRdx(i,j,N)=0.  ! no rotation at surface level
        enddo

!
!  Smooth vertically at interior points.
!
#  ifdef TS_MIX_ISO_FILT
# ifdef EXACT_RESTART
      if (mod(iic-1,ismooth).eq.0) then
# else
      if (mod(iic-ntstart,ismooth).eq.0) then
# endif
        do k=N-1,2,-1
          do i=IU_RANGE
            dRdx(i,j,k)=0.05*dRdx(i,j,k-1)+
     &                  0.90*dRdx(i,j,k  )+
     &                  0.05*dRdx(i,j,k+1)
          enddo
        enddo
        do i=IU_RANGE
          dRdx(i,j,1)=dRdx(i,j,2)
        enddo
      endif
#  endif /* TS_MIX_ISO_FILT */


!
# endif /* TS_DIF2 || TS_DIF4 || SPONGE_DIF2 */
!
!--------------------------------------------------
!  Couple ETA-component of velocity
!--------------------------------------------------
!
        if (j.ge.Jstr) then
          do i=IV_RANGE
            DC(i,0)=0.
            CF(i,0)=0.
            FC(i,0)=0.
          enddo
          do k=1,N,+1
            do i=IV_RANGE
              DC(i,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*v(i,j,k,nnew)
            enddo
          enddo
          do i=IV_RANGE
            DC(i,0)=1./DC(i,0)
            CF(i,0)=DC(i,0)*(CF(i,0)-DV_avg1(i,j,nnew))
# ifdef MRL_WCI
     &                                 +vst2d(i,j)
# endif
            vbar(i,j,knew)=DC(i,0)*DV_avg1(i,j,nnew)
# ifdef MRL_WCI
     &                                 -vst2d(i,j)
#  ifdef MASKING
            vbar(i,j,knew)=vbar(i,j,knew)*vmask(i,j)
     &                  +vst2d(i,j)*(vmask(i,j)-1.0)
#  endif
# endif
# ifdef WET_DRY
          vbar(i,j,knew)=vbar(i,j,knew)*vmask_wet(i,j)
# endif
          enddo
          do k=N,1,-1
            do i=IV_RANGE
              v(i,j,k,nnew)=(v(i,j,k,nnew)-CF(i,0))
# ifdef MASKING
     &                                  *vmask(i,j)
#  ifdef MRL_WCI
     &                 +vst(i,j,k)*(vmask(i,j)-1.0)
#  endif
# endif
# ifdef WET_DRY
              v(i,j,k,nnew)=v(i,j,k,nnew)*vmask_wet(i,j)
#  ifdef MRL_WCI
              vst(i,j,k)=vst(i,j,k)*vmask_wet(i,j)
#  endif
# endif

# ifdef TR
#  ifdef MRL_WCI
              FC(i,k)=DC(i,k)*(0.5*(v(i,j,k,nstp)+v(i,j,k,nnew))
     &                                              +vst(i,j,k))
#  else
              FC(i,k)=0.5*(v(i,j,k,nstp)+v(i,j,k,nnew))*DC(i,k)
#  endif
# else
#  ifdef MRL_WCI
              FC(i,k)=DELTA*Hvom(i,j,k)+
     &                EPSIL*DC(i,k)*(v(i,j,k,nstp)+v(i,j,k,nnew)
     &                                          +2.0*vst(i,j,k))
#  else
              FC(i,k)=DELTA*Hvom(i,j,k)+
     &                EPSIL*DC(i,k)*(v(i,j,k,nstp)+v(i,j,k,nnew))
#  endif
#  undef DELTA
#  undef EPSIL
# endif /* TR */

              FC(i,0)=FC(i,0)+FC(i,k)
            enddo
          enddo
          do i=IV_RANGE
            FC(i,0)=DC(i,0)*(FC(i,0)-DV_avg2(i,j))
          enddo
          do k=1,N,+1
            do i=IV_RANGE
              Hvom(i,j,k)=FC(i,k)-DC(i,k)*FC(i,0)
            enddo
          enddo
!
!-------------------------------------------------
! Compute density ETA gradients for use in t3dmix
!-------------------------------------------------
!
# if !defined TS_MIX_S && (defined TS_DIF4 || defined TS_DIF2)
          do k=1,N-1,+1
            do i=IV_RANGE
              cff=0.5*(pn(i,j)+pn(i,j-1)) SWITCH vmask(i,j)
#  ifdef  TS_MIX_ISO
#   ifdef SPLIT_EOS
              dpth=0.5*( z_w(i,j,N)+z_w(i,j-1,N)
     &                -  z_r(i,j,k)-z_r(i,j-1,k))
              dRe=cff*( rho1(i,j,k)-rho1(i,j-1,k)
     &               +(  qp1(i,j,k)- qp1(i,j-1,k))
     &                       *dpth*(1.-qp2*dpth) )
#   else
              dRe=cff*(rho1(i,j,k)-rho1(i,j-1,k))
#   endif
#   ifdef TS_MIX_ISO_FILT
              dRe  = cfilt*dRe + (1.-cfilt)*dRde(i,j,k)  ! exp. time smoothing
#   endif
              dRde(i,j,k)=dRe
#  elif defined TS_MIX_GEO
              dRde(i,j,k)=cff*(z_r(i,j,k)-z_r(i,j-1,k))  ! dRde gets dZde values
#  else
              dRde(i,j,k)=0.
#  endif
            enddo
          enddo
          do i=IV_RANGE
            dRde(i,j,N)=0.    ! no rotation at surface level
          enddo
!
!  Smooth vertically at interior points.
!
#  ifdef TS_MIX_ISO_FILT
# ifdef EXACT_RESTART
      if (mod(iic-1,ismooth).eq.0) then
# else
      if (mod(iic-ntstart,ismooth).eq.0) then
# endif
        do k=N-1,2,-1
          do i=IV_RANGE
            dRde(i,j,k)=0.05*dRde(i,j,k-1)+
     &                  0.90*dRde(i,j,k  )+
     &                  0.05*dRde(i,j,k+1)
          enddo
        enddo
        do i=IV_RANGE
          dRde(i,j,1)=dRde(i,j,2)
        enddo
      endif
#  endif /* TS_MIX_ISO_FILT */
!
# endif /* TS_DIF2 || TS_DIF4 */
!
!
        endif  !<-- j.ge.Jstr
      enddo    !<-- j loop

!
!--------------------------------------------------------------------
!
# if defined TS_HADV_TEST
#  if defined DIAGONAL_ADV || defined SOLID_BODY_ROT
      do k=1,N
        do j=jstrR,jendR
          do i=IstrR,iendR
            u(i,j,k,nnew) = u(i,j,k,nstp)
          enddo
        enddo
        do j=jstrR,jendR
          do i=IstrR,iendR
            v(i,j,k,nnew) = v(i,j,k,nstp)
          enddo
        enddo
      enddo
#  endif
#  if defined SOLID_BODY_PER
      cff = cos(2.*pi*(time+dt)/(float(ntimes)*dt))
      do k=1,N
        do j=jstrR,jendR
          do i=IstrR,iendR
             u(i,j,k,nnew) = u(i,j,k,nstp) * cff
             u(i,j,k,nrhs) = u(i,j,k,nstp)
          enddo
        enddo
        do j=JstrR,jendR
          do i=istrR,iendR
            v(i,j,k,nnew) = v(i,j,k,nstp) * cff
            v(i,j,k,nrhs) = v(i,j,k,nstp)
          enddo
        enddo
      enddo
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
#   ifdef THREE_GHOST_POINTS_UV
      call exchange_u3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             u(START_2D_ARRAY,1,nrhs))
      call exchange_v3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             v(START_2D_ARRAY,1,nrhs))
#   else
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        u(START_2D_ARRAY,1,nrhs))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        v(START_2D_ARRAY,1,nrhs))
#   endif
#  endif

      do j=Jstr,Jend
        do k=1,N
          do i=Istr,Iend
            Huon(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
     &                                         *u(i,j,k,nnew)
          enddo
          do i=Istr,Iend
            if (j.ge.JstrV) then
              Hvom(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
     &                                           *v(i,j,k,nnew)
            end if
          enddo
        enddo
      enddo
# endif
!
!--------------------------------------------------------------------
! Apply mass point sources
!--------------------------------------------------------------------
!
! This is done after the coupling is performed (for now) to avoid
! having (u,v) and (Huon,Hvom) affected by masked values or
! variables not containing mass point sources info during the
! coupling process (P. Marchesiello, 17 Jun 2004).
!
# ifdef PSOURCE
      do is=1,Nsrc
#  ifdef MPI
           i=Isrc_mpi(is,mynode)
           j=Jsrc_mpi(is,mynode)
#  else
            i=Isrc(is)
            j=Jsrc(is)
#  endif
        if (IstrR.le.i .and. i.le.IendR .and.
     &      JstrR.le.j .and. j.le.JendR) then
          if (Dsrc(is).eq.0) then
            do k=1,N
              u(i,j,k,nnew)=2.*Qsrc(is,k)/( on_u(i,j)*(
     &                       z_w(i-1,j,k)-z_w(i-1,j,k-1)
     &                      +z_w(i  ,j,k)-z_w(i  ,j,k-1)
     &                                                ))
              Huon(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)*
     &                                            u(i,j,k,nnew)
            enddo
          elseif (Dsrc(is).eq.1) then
            do k=1,N
              v(i,j,k,nnew)=2.*Qsrc(is,k)/( om_v(i,j)*(
     &                       z_w(i,j-1,k)-z_w(i,j-1,k-1)
     &                      +z_w(i,j  ,k)-z_w(i,j  ,k-1)
     &                                                ))
              Hvom(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)*
     &                                            v(i,j,k,nnew)
            enddo
          endif
        endif
      enddo
# endif
!
!--------------------------------------------------------------------
!  Set 3D Momemtum nudging
!--------------------------------------------------------------------
!
# if defined M3NUDGING && defined M3CLIMATOLOGY

#  ifdef ZONAL_NUDGING
      if (iic.eq.ntstart .or. mod(iic,10).eq.0) then
        call zonavg_3d(istr,iend,jstr,jend,
     &                        u(START_2D_ARRAY,1,nnew),uzon)
        call zonavg_3d(istr,iend,jstr,jend,
     &                        v(START_2D_ARRAY,1,nnew),vzon)
      endif
      if (iic.eq.ntstart) then
        call zonavg_3d(istr,iend,jstr,jend,
     &                        uclm(START_2D_ARRAY,1),uclmzon)
        call zonavg_3d(istr,iend,jstr,jend,
     &                        vclm(START_2D_ARRAY,1),vclmzon)
      endif
#  endif /* ZONAL_NUDGING */

      do k=1,N
        do j=JU_RANGE
          do i=IU_RANGE
            u(i,j,k,nnew)=u(i,j,k,nnew)+dt*M3nudgcof(i,j)*
#  ifdef ZONAL_NUDGING
     &                           (uclmzon(j,k)-uzon(j,k))
#  else
     &                           (uclm(i,j,k)-u(i,j,k,nnew))
#  endif /* ZONAL_NUDGING */
#  ifdef MASKING
     &                                           *umask(i,j)
#  endif
# ifdef WET_DRY
     &                                       *umask_wet(i,j)
# endif
          enddo
        enddo
      enddo
      do k=1,N
        do j=JV_RANGE
          do i=IV_RANGE
            v(i,j,k,nnew)=v(i,j,k,nnew)+dt*M3nudgcof(i,j)*
#  ifdef ZONAL_NUDGING
     &                            (vclmzon(j,k)-vzon(j,k))
#  else
     &                            (vclm(i,j,k)-v(i,j,k,nnew))
#  endif
#  ifdef MASKING
     &                                            *vmask(i,j)
#  endif
# ifdef WET_DRY
     &                                        *vmask_wet(i,j)
# endif
          enddo
        enddo
      enddo
# endif
!
!--------------------------------------------------------------------
!  Compute UV & PV diagnostics
!--------------------------------------------------------------------
!
#  if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_PV
      do k=1,N
        do j=JU_RANGE
          do i=IU_RANGE
            MVmix2(i,j,k,1)=( 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nnew)
     &              -0.5*(Hz_bak(i,j,k)+Hz_bak(i-1,j,k))*u(i,j,k,nstp) )
     &                                           *om_u(i,j)*on_u(i,j)/dt
     &                                           - Mrate(i,j,k,1)
            Mrate(i,j,k,1)= MVmix2(i,j,k,1) + Mrate(i,j,k,1)
#  if defined UV_HADV_UP3  || defined UV_HADV_UP5
            MHdiff(i,j,k,1)= MXadv(i,j,k,1) + MYadv(i,j,k,1)
     &                       - MHdiff(i,j,k,1)
# endif
!
! Divide all diagnostic terms by the cell volume Hz/(pm*pn).
! There after the unit of diag terms are :
! (unit of velocity) * s-1  =  m * s-2
!
            cff=2.*pm_u(i,j)*pn_u(i,j)/(Hz(i,j,k)+Hz(i-1,j,k))
     &                                           SWITCH umask(i,j)

            MXadv(i,j,k,1)   = MXadv(i,j,k,1)*cff
            MYadv(i,j,k,1)   = MYadv(i,j,k,1)*cff
# ifdef DIAGNOSTICS_UV
            MVadv(i,j,k,1)   = MVadv(i,j,k,1)*cff
            MCor(i,j,k,1)    = MCor(i,j,k,1)*cff
            MPrsgrd(i,j,k,1) = MPrsgrd(i,j,k,1)*cff
# if defined DIAGNOSTICS_BARO
            MBaro(i,j,k,1) = MBaro(i,j,k,1)*cff
# endif
# if defined M3FAST
            Mfast(i,j,k,1) = Mfast(i,j,k,1)*cff
# endif
# endif
            MHmix(i,j,k,1,3-nstp) = MHmix(i,j,k,1,3-nstp)*cff
            MVmix(i,j,k,1)   = MVmix(i,j,k,1)*cff
            MVmix2(i,j,k,1)   = MVmix2(i,j,k,1)*cff
            Mrate(i,j,k,1)   = Mrate(i,j,k,1)*cff
            MHdiff(i,j,k,1)   = MHdiff(i,j,k,1)*cff

#  ifdef MRL_WCI
            Mvf(i,j,k,1)     = Mvf(i,j,k,1)*cff
            Mbrk(i,j,k,1)    = Mbrk(i,j,k,1)*cff
            MStCo(i,j,k,1)   = MStCo(i,j,k,1)*cff
            MVvf(i,j,k,1)    = MVvf(i,j,k,1)*cff
            MPrscrt(i,j,k,1) = MPrscrt(i,j,k,1)*cff
            Msbk(i,j,k,1)    = Msbk(i,j,k,1)*cff
            Mbwf(i,j,k,1)    = Mbwf(i,j,k,1)*cff
            Mfrc(i,j,k,1)    = Mfrc(i,j,k,1)*cff
#  endif
          enddo
        enddo

        do j=JV_RANGE
          do i=IV_RANGE
            MVmix2(i,j,k,2)=( 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nnew)
     &              -0.5*(Hz_bak(i,j,k)+Hz_bak(i,j-1,k))*v(i,j,k,nstp) )
     &                                           *om_v(i,j)*on_v(i,j)/dt
     &                                           - Mrate(i,j,k,2)
            Mrate(i,j,k,2)= MVmix2(i,j,k,2) + Mrate(i,j,k,2)
#  if defined UV_HADV_UP3 || defined UV_HADV_UP5
            MHdiff(i,j,k,2)= MXadv(i,j,k,2) + MYadv(i,j,k,2)
     &                       - MHdiff(i,j,k,2)
#  endif
!
! Divide all diagnostic terms by the cell volume Hz/(pm*pn).
! There after the unit of diag terms are :
! (unit of velocity) * s-1  =  m * s-2

            cff=2.*pm_v(i,j)*pn_v(i,j)/(Hz(i,j,k)+Hz(i,j-1,k))
     &                                           SWITCH vmask(i,j)

            MXadv(i,j,k,2)   = MXadv(i,j,k,2)*cff
            MYadv(i,j,k,2)   = MYadv(i,j,k,2)*cff
# ifdef DIAGNOSTICS_UV
            MVadv(i,j,k,2)   = MVadv(i,j,k,2)*cff
            MCor(i,j,k,2)    = MCor(i,j,k,2)*cff
            MPrsgrd(i,j,k,2) = MPrsgrd(i,j,k,2)*cff
# if defined DIAGNOSTICS_BARO
            MBaro(i,j,k,2) = MBaro(i,j,k,2)*cff
# endif
# if defined M3FAST
            Mfast(i,j,k,2) = Mfast(i,j,k,2)*cff
# endif
#  endif
            MHmix(i,j,k,2,3-nstp) = MHmix(i,j,k,2,3-nstp)*cff
            MVmix(i,j,k,2)   = MVmix(i,j,k,2)*cff
            MVmix2(i,j,k,2)   = MVmix2(i,j,k,2)*cff
            Mrate(i,j,k,2)   = Mrate(i,j,k,2)*cff
            MHdiff(i,j,k,2)   = MHdiff(i,j,k,2)*cff

#  ifdef MRL_WCI
            Mvf(i,j,k,2)     = Mvf(i,j,k,2)*cff
            Mbrk(i,j,k,2)    = Mbrk(i,j,k,2)*cff
            MStCo(i,j,k,2)   = MStCo(i,j,k,2)*cff
            MVvf(i,j,k,2)    = MVvf(i,j,k,2)*cff
            MPrscrt(i,j,k,2) = MPrscrt(i,j,k,2)*cff
            Msbk(i,j,k,2)    = Msbk(i,j,k,2)*cff
            Mbwf(i,j,k,2)    = Mbwf(i,j,k,2)*cff
            Mfrc(i,j,k,2)    = Mfrc(i,j,k,2)*cff
#  endif
          enddo
        enddo
      enddo

# endif

# if defined DIAGNOSTICS_VRT && ! defined DIAGNOSTICS_UV

        do j=JU_RANGE
          do i=IU_RANGE
           k=1
           wrkrate(i,j,1)=(0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nnew)
     &           -0.5*(Hz_bak(i,j,k)+Hz_bak(i-1,j,k))*u(i,j,k,nstp))
     &                                       *om_u(i,j)*on_u(i,j)/dt
            do k=2,N
                wrkrate(i,j,1) = wrkrate(i,j,1) +
     &           ( 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nnew)
     &           -0.5*(Hz_bak(i,j,k)+Hz_bak(i-1,j,k))*u(i,j,k,nstp))
     &                                       *om_u(i,j)*on_u(i,j)/dt
            enddo

! Divide u-diagnostic terms by (1/pn).
! There after the unit of these terms are :
! m 3 s-2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            cff=0.25*(pm(i,j)+pm(i-1,j))
!     &             *(pn(i,j)+pn(i-1,j))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            cff=0.5*(pn(i,j)+pn(i-1,j))

            wrkXadv(i,j,1)=wrkXadv(i,j,1)*cff
            wrkYadv(i,j,1)=wrkYadv(i,j,1)*cff
            wrkHdiff(i,j,1)=wrkHdiff(i,j,1)*cff
            wrkCor(i,j,1)=wrkCor(i,j,1)*cff
            wrkPrsgrd(i,j,1)=wrkPrsgrd(i,j,1)*cff
            wrkHmix(i,j,1,3-nstp)=wrkHmix(i,j,1,3-nstp)*cff
            wrkVmix(i,j,1)=wrkVmix(i,j,1)*cff
            wrkVmix2(i,j,1)=wrkVmix2(i,j,1)*cff
            wrkrate(i,j,1)=wrkrate(i,j,1)*cff
# if defined DIAGNOSTICS_BARO
            wrkBaro(i,j,1)=wrkBaro(i,j,1)*cff
# endif
# if defined M3FAST
            wrkfast(i,j,1)=wrkfast(i,j,1)*cff
# endif
          enddo
        enddo
        do j=JV_RANGE
          do i=IV_RANGE
           k=1
           wrkrate(i,j,2)=(0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nnew)
     &           -0.5*(Hz_bak(i,j,k)+Hz_bak(i,j-1,k))*v(i,j,k,nstp))
     &                                       *om_v(i,j)*on_v(i,j)/dt
            do k=2,N
           wrkrate(i,j,2) = wrkrate(i,j,2) +
     &           ( 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nnew)
     &           -0.5*(Hz_bak(i,j,k)+Hz_bak(i,j-1,k))*v(i,j,k,nstp))
     &                                       *om_v(i,j)*on_v(i,j)/dt
            enddo

! Divide v-diagnostic terms by (1/pm).
! There after the unit of these terms are :
! m 3 s-2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            cff=0.25*(pm(i,j)+pm(i,j-1))
!     &             *(pn(i,j)+pn(i,j-1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            cff=0.5*(pm(i,j)+pm(i,j-1))

            wrkXadv(i,j,2)=wrkXadv(i,j,2)*cff
            wrkYadv(i,j,2)=wrkYadv(i,j,2)*cff
            wrkHdiff(i,j,2)=wrkHdiff(i,j,2)*cff
            wrkCor(i,j,2)=wrkCor(i,j,2)*cff
            wrkPrsgrd(i,j,2)=wrkPrsgrd(i,j,2)*cff
            wrkHmix(i,j,2,3-nstp)=wrkHmix(i,j,2,3-nstp)*cff
            wrkVmix(i,j,2)=wrkVmix(i,j,2)*cff
            wrkVmix2(i,j,2)=wrkVmix2(i,j,2)*cff
            wrkrate(i,j,2)=wrkrate(i,j,2)*cff
# if defined DIAGNOSTICS_BARO
            wrkBaro(i,j,2)=wrkBaro(i,j,2)*cff
# endif
# if defined M3FAST
            wrkfast(i,j,2)=wrkfast(i,j,2)*cff
# endif
          enddo
        enddo
# endif /* DIAGNOSTICS_VRT */

# if defined DIAGNOSTICS_EK_FULL
        do j=JstrR,JendR
          do i=IstrR,IendR
            ekwrkHadv(i,j,1) = 0.
            ekwrkVadv(i,j,1) = 0.
            ekwrkCor(i,j,1) = 0.
            ekwrkPrsgrd(i,j,1) = 0.
            ekwrkVmix(i,j,1) = 0.
            ekwrkHmix(i,j,1,3-nstp) = 0.
!            ekwrkVmix2(i,j,1) = 0.
            ekwrkHdiff(i,j,1) = 0.
# if defined DIAGNOSTICS_BARO
            ekwrkBaro(i,j,1) = 0.
# endif
# if defined M3FAST
            ekwrkfast(i,j,1) = 0.
# endif
            ekwrkHadv(i,j,2) = 0.
            ekwrkVadv(i,j,2) = 0.
            ekwrkCor(i,j,2) = 0.
            ekwrkPrsgrd(i,j,2) = 0.
            ekwrkVmix(i,j,2) = 0.
            ekwrkHmix(i,j,2,3-nstp) = 0.
!            ekwrkVmix2(i,j,2) = 0.
            ekwrkHdiff(i,j,2) = 0.
# if defined DIAGNOSTICS_BARO
            ekwrkBaro(i,j,2) = 0.
# endif
# if defined M3FAST
            ekwrkfast(i,j,2) = 0.
# endif

# if defined DIAGNOSTICS_EK_MLD
            ekwrkHadv_mld(i,j,1) = 0.
            ekwrkVadv_mld(i,j,1) = 0.
            ekwrkCor_mld(i,j,1) = 0.
            ekwrkPrsgrd_mld(i,j,1) = 0.
            ekwrkVmix_mld(i,j,1) = 0.
            ekwrkHmix_mld(i,j,1,3-nstp) = 0.
!            ekwrkVmix2_mld(i,j,1) = 0.
            ekwrkHdiff_mld(i,j,1) = 0.
# if defined DIAGNOSTICS_BARO
            ekwrkBaro_mld(i,j,1) = 0.
# endif
            ekwrkHadv_mld(i,j,2) = 0.
            ekwrkVadv_mld(i,j,2) = 0.
            ekwrkCor_mld(i,j,2) = 0.
            ekwrkPrsgrd_mld(i,j,2) = 0.
            ekwrkVmix_mld(i,j,2) = 0.
            ekwrkHmix_mld(i,j,2,3-nstp) = 0.
!            ekwrkVmix2_mld(i,j,2) = 0.
            ekwrkHdiff_mld(i,j,2) = 0.
# if defined DIAGNOSTICS_BARO
            ekwrkBaro_mld(i,j,2) = 0.
# endif
# endif
        enddo
      enddo

        do j=JU_RANGE
          do i=IU_RANGE
            do k=1,N
# if defined DIAGNOSTICS_UV  || defined DIAGNOSTICS_PV
            cff=(Hz(i,j,k)+Hz(i-1,j,k))/(2.*pm_u(i,j)*pn_u(i,j))
     &          * 0.5*(u(i,j,k,nstp)+u(i,j,k,nnew))
# else
            cff= 0.5*(u(i,j,k,nstp)+u(i,j,k,nnew))
# endif
#  ifdef MASKING
     &                     *umask(i,j)
#  endif
            ekwrkHadv(i,j,1) = ekwrkHadv(i,j,1)
     &        + (MXadv(i,j,k,1) + MYadv(i,j,k,1)) * cff
            ekwrkVmix(i,j,1) = ekwrkVmix(i,j,1)
     &        + MVmix(i,j,k,1)  * cff
            ekwrkHmix(i,j,1,3-nstp) = ekwrkHmix(i,j,1,3-nstp)
     &        + MHmix(i,j,k,1,3-nstp)  * cff
!            ekwrkVmix2(i,j,1) = ekwrkVmix2(i,j,1)
!     &        + MVmix2(i,j,k,1)  * cff
            ekwrkHdiff(i,j,1) = ekwrkHdiff(i,j,1)
     &        + MHdiff(i,j,k,1) * cff

# ifdef DIAGNOSTICS_UV
            cff=(Hz(i,j,k)+Hz(i-1,j,k))/(2.*pm_u(i,j)*pn_u(i,j))
     &          * 0.5*(u(i,j,k,nstp)+u(i,j,k,nnew))
# else
            cff= 0.5*(u(i,j,k,nstp)+u(i,j,k,nnew))
# endif
#  ifdef MASKING
     &                     *umask(i,j)
#  endif
            ekwrkVadv(i,j,1) = ekwrkVadv(i,j,1)
     &        + MVadv(i,j,k,1)  * cff
            ekwrkCor(i,j,1) = ekwrkCor(i,j,1)
     &        + MCor(i,j,k,1)  * cff
            ekwrkPrsgrd(i,j,1) = ekwrkPrsgrd(i,j,1)
     &        + MPrsgrd(i,j,k,1)  * cff
# if defined DIAGNOSTICS_BARO
            ekwrkBaro(i,j,1) = ekwrkBaro(i,j,1)
     &        + MBaro(i,j,k,1)  * cff
# endif
# if defined M3FAST
            ekwrkfast(i,j,1) = ekwrkfast(i,j,1)
     &        + Mfast(i,j,k,1)  * cff
# endif
          enddo

# if defined DIAGNOSTICS_EK_MLD

          do k=kbl(i,j),N
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_PV
            cff=(Hz(i,j,k)+Hz(i-1,j,k))/(2.*pm_u(i,j)*pn_u(i,j))
     &          * 0.5*(u(i,j,k,nstp)+u(i,j,k,nnew))
# else
            cff= 0.5*(u(i,j,k,nstp)+u(i,j,k,nnew))
# endif
#  ifdef MASKING
     &                     *umask(i,j)
#  endif
            ekwrkHadv_mld(i,j,1) = ekwrkHadv_mld(i,j,1)
     &        + (MXadv(i,j,k,1) + MYadv(i,j,k,1)) * cff

            ekwrkVmix_mld(i,j,1) = ekwrkVmix_mld(i,j,1)
     &        + MVmix(i,j,k,1)  * cff
            ekwrkHmix_mld(i,j,1,3-nstp) = ekwrkHmix_mld(i,j,1,3-nstp)
     &        + MHmix(i,j,k,1,3-nstp)  * cff
!            ekwrkVmix2_mld(i,j,1) = ekwrkVmix2_mld(i,j,1)
!     &        + MVmix2(i,j,k,1)  * cff
            ekwrkHdiff_mld(i,j,1) = ekwrkHdiff_mld(i,j,1)
     &        + MHdiff(i,j,k,1) * cff
# ifdef DIAGNOSTICS_UV
            cff=(Hz(i,j,k)+Hz(i-1,j,k))/(2.*pm_u(i,j)*pn_u(i,j))
     &          * 0.5*(u(i,j,k,nstp)+u(i,j,k,nnew))
# else
            cff= 0.5*(u(i,j,k,nstp)+u(i,j,k,nnew))
# endif
#  ifdef MASKING
     &                     *umask(i,j)
#  endif
            ekwrkVadv_mld(i,j,1) = ekwrkVadv_mld(i,j,1)
     &        + MVadv(i,j,k,1)  * cff
            ekwrkCor_mld(i,j,1) = ekwrkCor_mld(i,j,1)
     &        + MCor(i,j,k,1)  * cff
            ekwrkPrsgrd_mld(i,j,1) = ekwrkPrsgrd_mld(i,j,1)
     &        + MPrsgrd(i,j,k,1)  * cff
# if defined DIAGNOSTICS_BARO
            ekwrkBaro_mld(i,j,1) = ekwrkBaro_mld(i,j,1)
     &        + MBaro(i,j,k,1)  * cff
# endif
          enddo
# endif
        enddo
      enddo

      do j=JV_RANGE
        do i=IV_RANGE
            do k=1,N
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_PV
            cff=(Hz(i,j,k)+Hz(i,j-1,k))/(2.*pm_v(i,j)*pn_v(i,j))
     &          * 0.5*(v(i,j,k,nstp)+v(i,j,k,nnew))
# else
            cff= 0.5*(v(i,j,k,nstp)+v(i,j,k,nnew))
# endif
#  ifdef MASKING
     &                     *vmask(i,j)
#  endif
            ekwrkHadv(i,j,2) = ekwrkHadv(i,j,2)
     &        + (MXadv(i,j,k,2) + MYadv(i,j,k,2)) * cff
            ekwrkVmix(i,j,2) = ekwrkVmix(i,j,2)
     &        + MVmix(i,j,k,2)  * cff
            ekwrkHmix(i,j,2,3-nstp) = ekwrkHmix(i,j,2,3-nstp)
     &        + MHmix(i,j,k,2,3-nstp)  * cff
!            ekwrkVmix2(i,j,2) = ekwrkVmix2(i,j,2)
!     &        + MVmix2(i,j,k,2)  * cff
            ekwrkHdiff(i,j,2) = ekwrkHdiff(i,j,2)
     &        + MHdiff(i,j,k,2) * cff
# ifdef DIAGNOSTICS_UV
            cff=(Hz(i,j,k)+Hz(i,j-1,k))/(2.*pm_v(i,j)*pn_v(i,j))
     &          * 0.5*(v(i,j,k,nstp)+v(i,j,k,nnew))
# else
            cff= 0.5*(v(i,j,k,nstp)+v(i,j,k,nnew))
# endif
#  ifdef MASKING
     &                     *vmask(i,j)
#  endif
            ekwrkVadv(i,j,2) = ekwrkVadv(i,j,2)
     &        + MVadv(i,j,k,2)  * cff
            ekwrkCor(i,j,2) = ekwrkCor(i,j,2)
     &        + MCor(i,j,k,2)  * cff
            ekwrkPrsgrd(i,j,2) = ekwrkPrsgrd(i,j,2)
     &        + MPrsgrd(i,j,k,2)  * cff
# if defined DIAGNOSTICS_BARO
            ekwrkBaro(i,j,2) = ekwrkBaro(i,j,2)
     &        + MBaro(i,j,k,2)  * cff
# endif
# if defined M3FAST
            ekwrkfast(i,j,2) = ekwrkfast(i,j,2)
     &        + Mfast(i,j,k,2)  * cff
# endif
            enddo

# if defined DIAGNOSTICS_EK_MLD
          do k=kbl(i,j),N
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_PV
            cff=(Hz(i,j,k)+Hz(i,j-1,k))/(2.*pm_v(i,j)*pn_v(i,j))
     &          * 0.5*(v(i,j,k,nstp)+v(i,j,k,nnew))
# else
            cff= 0.5*(v(i,j,k,nstp)+v(i,j,k,nnew))
# endif
#  ifdef MASKING
     &                     *vmask(i,j)
#  endif
            ekwrkHadv_mld(i,j,2) = ekwrkHadv_mld(i,j,2)
     &        + (MXadv(i,j,k,2) + MYadv(i,j,k,2)) * cff
            ekwrkVmix_mld(i,j,2) = ekwrkVmix_mld(i,j,2)
     &        + MVmix(i,j,k,2)  * cff
            ekwrkHmix_mld(i,j,2,3-nstp) = ekwrkHmix_mld(i,j,2,3-nstp)
     &        + MHmix(i,j,k,2,3-nstp)  * cff
!            ekwrkVmix2_mld(i,j,2) = ekwrkVmix2_mld(i,j,2)
!     &        + MVmix2(i,j,k,2)  * cff
            ekwrkHdiff_mld(i,j,2) = ekwrkHdiff_mld(i,j,2)
     &        + MHdiff(i,j,k,2) * cff
# ifdef DIAGNOSTICS_UV
            cff=(Hz(i,j,k)+Hz(i,j-1,k))/(2.*pm_v(i,j)*pn_v(i,j))
     &          * 0.5*(v(i,j,k,nstp)+v(i,j,k,nnew))
# else
            cff= 0.5*(v(i,j,k,nstp)+v(i,j,k,nnew))
# endif
#  ifdef MASKING
     &                     *vmask(i,j)
#  endif
            ekwrkVadv_mld(i,j,2) = ekwrkVadv_mld(i,j,2)
     &        + MVadv(i,j,k,2)  * cff
            ekwrkCor_mld(i,j,2) = ekwrkCor_mld(i,j,2)
     &        + MCor(i,j,k,2)  * cff
            ekwrkPrsgrd_mld(i,j,2) = ekwrkPrsgrd_mld(i,j,2)
     &        + MPrsgrd(i,j,k,2)  * cff
# if defined DIAGNOSTICS_BARO
            ekwrkBaro_mld(i,j,2) = ekwrkBaro_mld(i,j,2)
     &        + MBaro(i,j,k,2)  * cff
# endif
          enddo
# endif
          enddo
        enddo
# endif /* DIAGNOSTICS_EK */

# if defined DIAGNOSTICS_EK

        do j=JU_RANGE
          do i=IU_RANGE

            ekwrkrate(i,j,1) = 0.
            ekwrkvol(i,j,1) = 0.

            do k=1,N

                ekwrkrate(i,j,1) = ekwrkrate(i,j,1) + 0.5*
     &            ( 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nnew)**2
     &         -0.5*(Hz_bak(i,j,k)+Hz_bak(i-1,j,k))*u(i,j,k,nstp)**2)
     &                                        *om_u(i,j)*on_u(i,j)/dt
#  ifdef MASKING
     &                     *umask(i,j)
#  endif
            ekwrkvol(i,j,1) = ekwrkvol(i,j,1)
     &                  + (0.5*(Hz_bak(i,j,k)+Hz_bak(i-1,j,k))
     &                  - 0.5*(Hz(i,j,k)+Hz(i-1,j,k))      )
     &                  * ((Hz_bak(i,j,k)+Hz_bak(i-1,j,k))
     &                  / (Hz(i,j,k)+Hz(i-1,j,k))
     &                  * 0.5 * u(i,j,k,nstp)**2
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL
     &                  + u(i,j,k,nstp) * (
     &                    MXadv(i,j,k,1)
     &                  + MYadv(i,j,k,1)
     &                  + MVadv(i,j,k,1)
     &                  + MCor(i,j,k,1)
     &                  + MPrsgrd(i,j,k,1)
     &                  + MHmix(i,j,k,1,3-nstp)
     &                  + MVmix(i,j,k,1)
# if defined M3FAST
     &                  + Mfast(i,j,k,1)
# endif
# if defined DIAGNOSTICS_UV
     &                  + MVmix2(i,j,k,1))
# else
     &                  ) / (Hz(i,j,k)+Hz(i-1,j,k))
     &                    * (2.*pm_u(i,j)*pn_u(i,j))
# endif
#  endif
     &                                     )  *om_u(i,j)*on_u(i,j)/dt
#  ifdef MASKING
     &                     *umask(i,j)
#  endif
            enddo

            ekwrkVmix2(i,j,1) = ekwrkrate(i,j,1) - ekwrkVmix2(i,j,1)


# if defined DIAGNOSTICS_EK_MLD

            ekwrkrate_mld(i,j,1) = 0.
            ekwrkvol_mld(i,j,1) = 0.

            do k=kbl(i,j),N

                ekwrkrate_mld(i,j,1) = ekwrkrate_mld(i,j,1) + 0.5*
     &            ( 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nnew)**2
     &         -0.5*(Hz_bak(i,j,k)+Hz_bak(i-1,j,k))*u(i,j,k,nstp)**2)
     &                                        *om_u(i,j)*on_u(i,j)/dt
#  ifdef MASKING
     &                     *umask(i,j)
#  endif

            ekwrkvol_mld(i,j,1) = ekwrkvol_mld(i,j,1)
     &                  + (0.5*(Hz_bak(i,j,k)+Hz_bak(i-1,j,k))
     &                  - 0.5*(Hz(i,j,k)+Hz(i-1,j,k))      )
     &                  * ((Hz_bak(i,j,k)+Hz_bak(i-1,j,k))
     &                  / (Hz(i,j,k)+Hz(i-1,j,k))
     &                  * 0.5 * u(i,j,k,nstp)**2
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL
     &                  + u(i,j,k,nstp) * (
     &                    MXadv(i,j,k,1)
     &                  + MYadv(i,j,k,1)
     &                  + MVadv(i,j,k,1)
     &                  + MCor(i,j,k,1)
     &                  + MPrsgrd(i,j,k,1)
     &                  + MHmix(i,j,k,1,3-nstp)
     &                  + MVmix(i,j,k,1)
# if defined DIAGNOSTICS_UV
     &                  + MVmix2(i,j,k,1))
# else
     &                  ) / (Hz(i,j,k)+Hz(i-1,j,k))
     &                    * (2.*pm_u(i,j)*pn_u(i,j))
# endif
#  endif
     &                                     )  *om_u(i,j)*on_u(i,j)/dt
#  ifdef MASKING
     &                     *umask(i,j)
#  endif

            enddo

            ekwrkVmix2_mld(i,j,1) = ekwrkrate_mld(i,j,1) - ekwrkVmix2_mld(i,j,1)
#  endif

          enddo
        enddo

        do j=JV_RANGE
          do i=IV_RANGE

            ekwrkrate(i,j,2)= 0.
            ekwrkvol(i,j,2) = 0.

            do k=1,N

              ekwrkrate(i,j,2) = ekwrkrate(i,j,2) + 0.5*
     &           ( 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nnew)**2
     &         -0.5*(Hz_bak(i,j,k)+Hz_bak(i,j-1,k))*v(i,j,k,nstp)**2)
     &                                        *om_v(i,j)*on_v(i,j)/dt
#  ifdef MASKING
     &                     *vmask(i,j)
#  endif
            ekwrkvol(i,j,2) = ekwrkvol(i,j,2)
     &                  + (0.5*(Hz_bak(i,j,k)+Hz_bak(i,j-1,k))
     &                  - 0.5*(Hz(i,j,k)+Hz(i,j-1,k))      )
     &                  * ((Hz_bak(i,j,k)+Hz_bak(i,j-1,k))
     &                  / (Hz(i,j,k)+Hz(i,j-1,k))
     &                  * 0.5 * v(i,j,k,nstp)**2
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL
     &                  + v(i,j,k,nstp) * (
     &                    MXadv(i,j,k,2)
     &                  + MYadv(i,j,k,2)
     &                  + MVadv(i,j,k,2)
     &                  + MCor(i,j,k,2)
     &                  + MPrsgrd(i,j,k,2)
     &                  + MHmix(i,j,k,2,3-nstp)
     &                  + MVmix(i,j,k,2)
# if defined M3FAST
     &                  + Mfast(i,j,k,2)
# endif
# if defined DIAGNOSTICS_UV
     &                  + MVmix2(i,j,k,2))
# else
     &                  ) / (Hz(i,j,k)+Hz(i,j-1,k))
     &                    * (2.*pm_v(i,j)*pn_v(i,j))
# endif
#  endif
     &                                     )   *om_v(i,j)*on_v(i,j)/dt
#  ifdef MASKING
     &                     *vmask(i,j)
#  endif
            enddo

            ekwrkVmix2(i,j,2) = ekwrkrate(i,j,2) - ekwrkVmix2(i,j,2)

# if defined DIAGNOSTICS_EK_MLD

            ekwrkrate_mld(i,j,2)= 0.
            ekwrkvol_mld(i,j,2) = 0.

            do k=kbl(i,j),N

              ekwrkrate_mld(i,j,2) = ekwrkrate_mld(i,j,2) + 0.5*
     &           ( 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nnew)**2
     &         -0.5*(Hz_bak(i,j,k)+Hz_bak(i,j-1,k))*v(i,j,k,nstp)**2)
     &                                        *om_v(i,j)*on_v(i,j)/dt
#  ifdef MASKING
     &                     *vmask(i,j)
#  endif

            ekwrkvol_mld(i,j,2) = ekwrkvol_mld(i,j,2)
     &                  + (0.5*(Hz_bak(i,j,k)+Hz_bak(i,j-1,k))
     &                  - 0.5*(Hz(i,j,k)+Hz(i,j-1,k))      )
     &                  * ((Hz_bak(i,j,k)+Hz_bak(i,j-1,k))
     &                  / (Hz(i,j,k)+Hz(i,j-1,k))
     &                  * 0.5 * v(i,j,k,nstp)**2
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL
     &                  + v(i,j,k,nstp) * (
     &                    MXadv(i,j,k,2)
     &                  + MYadv(i,j,k,2)
     &                  + MVadv(i,j,k,2)
     &                  + MCor(i,j,k,2)
     &                  + MPrsgrd(i,j,k,2)
     &                  + MHmix(i,j,k,2,3-nstp)
     &                  + MVmix(i,j,k,2)
# if defined DIAGNOSTICS_UV
     &                  + MVmix2(i,j,k,2))
# else
     &                  ) / (Hz(i,j,k)+Hz(i,j-1,k))
     &                    * (2.*pm_v(i,j)*pn_v(i,j))
# endif
#  endif
     &                                     )   *om_v(i,j)*on_v(i,j)/dt
#  ifdef MASKING
     &                     *vmask(i,j)
#  endif
            enddo

            ekwrkVmix2_mld(i,j,2) = ekwrkrate_mld(i,j,2) - ekwrkVmix2_mld(i,j,2)
#  endif
          enddo
        enddo
# endif /* DIAGNOSTICS_EK */




#if defined DIAGNOSTICS_PV
! Add non conservative terms to rhs.

          do k=1,N
            do j=JU_RANGE
              do i=IU_RANGE
# if defined DIAGNOSTICS_PV_FULL
                pv(i,j,k) = Mrhs(i,j,k,1)
# endif
                Mrhs(i,j,k,1) = (MVmix2(i,j,k,1)
     &                            + MVmix(i,j,k,1)
# if defined DIAGNOSTICS_BARO
     &                            - MBaro(i,j,k,1)
# endif
# if defined M3FAST
     &                            + Mfast(i,j,k,1)
# endif
     &                            + MHmix(i,j,k,1,3-nstp)
     &                            + MHdiff(i,j,k,1) )
# if defined DIAGNOSTICS_DISS
     &        * (Hz(i,j,k)+Hz(i-1,j,k))/(2.*pm_u(i,j)*pn_u(i,j))
     &        * 0.5*(u(i,j,k,nstp)+u(i,j,k,nnew))
# if ! defined DIAGNOSTICS_PV_FULL
     &                            - Mrhs(i,j,k,1)
# endif
# endif
              enddo
            enddo
          enddo

          do k=1,N
            do j=JV_RANGE
              do i=IV_RANGE
# if defined DIAGNOSTICS_PV_FULL
                pvd(i,j,k) = Mrhs(i,j,k,2)
# endif
                Mrhs(i,j,k,2) = (MVmix2(i,j,k,2)
     &                            + MVmix(i,j,k,2)
# if defined DIAGNOSTICS_BARO
     &                            - MBaro(i,j,k,2)
# endif
# if defined M3FAST
     &                            + Mfast(i,j,k,2)
# endif
     &                            + MHmix(i,j,k,2,3-nstp)
     &                            + MHdiff(i,j,k,2))
# if defined DIAGNOSTICS_DISS
     &        * (Hz(i,j,k)+Hz(i,j-1,k))/(2.*pm_v(i,j)*pn_v(i,j))
     &        * 0.5*(v(i,j,k,nstp)+v(i,j,k,nnew))
# if ! defined DIAGNOSTICS_PV_FULL
     &                            - Mrhs(i,j,k,2)
# endif
# endif
              enddo
            enddo
          enddo

#endif /* DIAGNOSTICS_PV */

# undef IU_RANGE
# undef JU_RANGE
# undef IV_RANGE
# undef JV_RANGE
!
!--------------------------------------------------------------------
! Exchange periodic boundaries and computational margins.
!--------------------------------------------------------------------
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
#  ifdef THREE_GHOST_POINTS_UV
      call exchange_u3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             u(START_2D_ARRAY,1,nnew))
      call exchange_v3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             v(START_2D_ARRAY,1,nnew))
#  else
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        u(START_2D_ARRAY,1,nnew))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        v(START_2D_ARRAY,1,nnew))
#  endif
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        Huon(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hvom(START_2D_ARRAY,1))

      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        ubar(START_2D_ARRAY,knew))

      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        vbar(START_2D_ARRAY,knew))
#  if defined TS_MIX_ISO || defined TS_MIX_GEO
      call exchange_u3d_tile (istr,iend,jstr,jend, dRdx  )
      call exchange_v3d_tile (istr,iend,jstr,jend, dRde  )
#  endif



! Need to exchange momentum terms for computation of vorticity terms
# if defined DIAGNOSTICS_UV && defined DIAGNOSTICS_VRT

      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MXadv(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MYadv(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MVadv(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MPrsgrd(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MCor(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                            MHmix(START_2D_ARRAY,1,1,3-nstp))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MHdiff(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MVmix(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 Mrate(START_2D_ARRAY,1,1))
# if defined DIAGNOSTICS_BARO
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MBaro(START_2D_ARRAY,1,1))
# endif
# if defined M3FAST
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 Mfast(START_2D_ARRAY,1,1))
# endif
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkWind(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkDrag(START_2D_ARRAY,1))



      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MXadv(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MYadv(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MVadv(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MPrsgrd(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MCor(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                            MHmix(START_2D_ARRAY,1,2,3-nstp))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MHdiff(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MVmix(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 Mrate(START_2D_ARRAY,1,2))
# if defined DIAGNOSTICS_BARO
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MBaro(START_2D_ARRAY,1,2))
# endif
# if defined M3FAST
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 Mfast(START_2D_ARRAY,1,2))
# endif
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkWind(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkDrag(START_2D_ARRAY,2))


# elif defined DIAGNOSTICS_VRT && ! defined DIAGNOSTICS_UV


      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkXadv(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkYadv(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkHdiff(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkPrsgrd(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkCor(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                              wrkHmix(START_2D_ARRAY,1,3-nstp))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkVmix(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkVmix2(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkrate(START_2D_ARRAY,1))
# if defined DIAGNOSTICS_BARO
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkBaro(START_2D_ARRAY,1))
# endif
# if defined M3FAST
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkfast(START_2D_ARRAY,1))
# endif
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkWind(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkDrag(START_2D_ARRAY,1))



      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkXadv(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkYadv(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkHdiff(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkPrsgrd(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkCor(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                              wrkHmix(START_2D_ARRAY,2,3-nstp))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkVmix(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkVmix2(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkrate(START_2D_ARRAY,2))
# if defined DIAGNOSTICS_BARO
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkBaro(START_2D_ARRAY,2))
# endif
# if defined M3FAST
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkfast(START_2D_ARRAY,2))
# endif
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkWind(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkDrag(START_2D_ARRAY,2))


# endif


# if defined TENDENCY



      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MXadv(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MYadv(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MVadv(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MPrsgrd(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MCor(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                            MHmix(START_2D_ARRAY,1,1,3-nstp))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MHdiff(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MVmix(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 Mrate(START_2D_ARRAY,1,1))
# if defined DIAGNOSTICS_BARO
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MBaro(START_2D_ARRAY,1,1))
# endif
# if defined M3FAST
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 Mfast(START_2D_ARRAY,1,1))
# endif
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MXadv(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MYadv(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MVadv(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MPrsgrd(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MCor(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                            MHmix(START_2D_ARRAY,1,2,3-nstp))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MHdiff(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MVmix(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 Mrate(START_2D_ARRAY,1,2))
# if defined DIAGNOSTICS_BARO
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MBaro(START_2D_ARRAY,1,2))
# endif
# if defined M3FAST
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 Mfast(START_2D_ARRAY,1,2))
# endif
# endif


# if defined DIAGNOSTICS_EK

      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkHadv(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkHdiff(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkVadv(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkPrsgrd(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkCor(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                          ekwrkHmix(START_2D_ARRAY,1,3-nstp))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkVmix(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkVmix2(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkrate(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkvol(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkWind(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkDrag(START_2D_ARRAY,1))
# if defined DIAGNOSTICS_BARO
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkBaro(START_2D_ARRAY,1))
# endif
# if defined M3FAST
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkfast(START_2D_ARRAY,1))
# endif
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkHadv(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkHdiff(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkVadv(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkPrsgrd(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkCor(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                          ekwrkHmix(START_2D_ARRAY,2,3-nstp))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkVmix(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkVmix2(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkrate(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkvol(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkWind(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkDrag(START_2D_ARRAY,2))
# if defined DIAGNOSTICS_BARO
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkBaro(START_2D_ARRAY,2))
# endif
# if defined M3FAST
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkfast(START_2D_ARRAY,2))
# endif
# if defined DIAGNOSTICS_EK_MLD

      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkHadv_mld(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkHdiff_mld(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkVadv_mld(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkPrsgrd_mld(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkCor_mld(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                          ekwrkHmix_mld(START_2D_ARRAY,1,3-nstp))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkVmix_mld(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkVmix2_mld(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkrate_mld(START_2D_ARRAY,1))
# if defined DIAGNOSTICS_BARO
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkBaro_mld(START_2D_ARRAY,1))
# endif
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkHadv_mld(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkHdiff_mld(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkVadv_mld(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkPrsgrd_mld(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkCor_mld(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                          ekwrkHmix_mld(START_2D_ARRAY,2,3-nstp))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkVmix_mld(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkVmix2_mld(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkrate_mld(START_2D_ARRAY,2))
# if defined DIAGNOSTICS_BARO
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 ekwrkBaro_mld(START_2D_ARRAY,2))
# endif
# endif

# endif



# endif
      return
      end
#else
      subroutine step3d_uv_empty
      return
      end
#endif /* SOLVE3D */
