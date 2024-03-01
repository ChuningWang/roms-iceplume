! $Id:$
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
#ifndef CHILD
!
# include "cppdefs.h"
# ifdef M3FAST
      subroutine unbq_bc_tile (Istr,Iend,Jstr,Jend,grad)
#  ifdef AGRIF      
      use AGRIF_Util
      integer Istr,Iend,Jstr,Jend
      real grad(PRIVATE_2D_SCRATCH_ARRAY)
      if (AGRIF_Root()) then
        call unbq_bc_parent_tile (Istr,Iend,Jstr,Jend,grad)
      else
        call unbq_bc_child_tile (Istr,Iend,Jstr,Jend,grad)
      endif
      return
      end
!
! PARENT
!
      subroutine unbq_bc_parent_tile (Istr,Iend,Jstr,Jend,grad)
#  endif
!
! Set lateral boundary conditions for U-component momentum
! qdmu_nbq for the parent grid.
!
# endif /* M3FAST */
#else
# ifdef M3FAST
!
! CHILD
!
      subroutine unbq_bc_child_tile (Istr,Iend,Jstr,Jend,grad)
!
! Set lateral boundary conditions for U-component momentum
! qdmu_nbq for the child grid.
!
# endif /* M3FAST */
#endif /* CHILD */
#ifdef M3FAST
!
! Common Code
!
# include "set_obc_definitions.h"
!
      implicit none
# include "param.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "climat.h"
# include "scalars.h"
# include "boundary.h"  
# ifdef MRL_WCI
#  include "forces.h"
# endif
#include "nbq.h"
      integer Istr,Iend,Jstr,Jend, i,j,k
      real    grad(PRIVATE_2D_SCRATCH_ARRAY)
      real    cff,eps,cff1,cff2,
     &        cx,cy, dft,dfx,dfy, tau,tau_in,tau_out
      parameter (eps=1.E-20)
!
# include "compute_auxiliary_bounds.h"
!
! Interpolations of the parent values to get unbqbry_east or unbqclm
!
# ifdef CHILD
      call unbq_bc_interp_tile(Istr,Iend,Jstr,Jend)
# endif
!
# if defined NBQ_FRC_BRY || defined NBQ_NUDGING
      tau_in =dtfast*tauM_in
      tau_out=dtfast*tauM_out
# endif

# ifndef EW_COM_PERIODIC
!
!====================================================================
!                            WESTERN BC
!====================================================================
      if (WESTERN_EDGE) then
!         write(6,*) 'WEST'
#  ifdef OBC_COM_WEST
#   if defined OBC_COM_NBQSPECIFIED || defined OBC_COM_NBQSPECIFIED_WEST
!                                            Western edge Specified BC
!                                            ======= ==== ========= ==
        do k=1,N
          do j=Jstr,Jend
#    ifdef NBQ_FRC_BRY
            qdmu_nbq(Istr,j,k)=unbqbry_west(j,k)         ! specified
#    else
            qdmu_nbq(Istr,j,k)=unbqclm(Istr,j,k)
#    endif
#    ifdef MASKING
     &                            *umask(Istr,j)
#    endif
#    ifdef WET_DRY
     &                        *umask_wet(Istr,j)
#    endif
          enddo
        enddo

#   elif defined OBC_COM_NBQORLANSKI
        do k=1,N                             ! Western edge radiation
          do j=Jstr,Jend+1                   ! ======= ==== =========
            grad(Istr  ,j)=(qdmu_nbq_west(j  ,k,1)
     &                     -qdmu_nbq_west(j-1,k,1))
#    ifdef MASKING
     &                             *pmask(Istr  ,j)
#    endif
            grad(Istr+1,j)=(qdmu_nbq_west(j  ,k,2)
     &                     -qdmu_nbq_west(j-1,k,2))
#    ifdef MASKING
     &                             *pmask(Istr+1,j)
#    endif
          enddo
          do j=Jstr,Jend
            dft=qdmu_nbq_west(j,k,2)-qdmu_nbq(Istr+1,j,k)
            dfx=qdmu_nbq(Istr+1,j,k)-qdmu_nbq(Istr+2,j,k)
            if (dfx*dft .lt. 0.) then
              dft=0.                       ! <-- cancel cx, if inflow
#    if defined NBQ_FRC_BRY || defined NBQ_NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
            if (dft*(grad(Istr+1,j)+grad(Istr+1,j+1)) .gt. 0.) then
              dfy=grad(Istr+1,j)
            else
              dfy=grad(Istr+1,j+1)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
            cy=min(cff,max(dft*dfy,-cff))
            qdmu_nbq(Istr,j,k)=(  cff*qdmu_nbq_west(j,k,1)
     &                            +cx*qdmu_nbq(Istr+1,j,k)
     &                          -max(cy,0.)*grad(Istr,j  )
     &                          -min(cy,0.)*grad(Istr,j+1)
     &                                          )/(cff+cx)
#    if defined NBQ_FRC_BRY  || defined NBQ_NUDGING 
            qdmu_nbq(Istr,j,k)=(1.-tau)*qdmu_nbq(Istr,j,k)
#     ifdef NBQ_FRC_BRY
     &                              +tau*unbqbry_west(j,k)
#     else     
     &                              +tau*unbqclm(Istr,j,k)
#     endif
#    endif
#    ifdef MASKING
     &                                      *umask(Istr,j)
#    endif
#    ifdef WET_DRY
     &                                  *umask_wet(Istr,j)
#    endif
          enddo
        enddo
!
#   else
!                                            Western edge gradient BC
!                                            ======= ==== ======== ==
        do k=1,N
          do j=Jstr,Jend
              qdmu_nbq(Istr,j,k)=qdmu_nbq(Istr+1,j,k)
#    ifdef MASKING
     &                                 *umask(Istr,j)
#    endif
#    ifdef WET_DRY
     &                             *umask_wet(Istr,j)
#    endif 
          enddo
        enddo
#   endif /* OBC_NBQORLANSKI */
!
!  Replace external mode
!
#   ifdef NBQ
#    define usum grad
        usum=0.
        do k=1,N
          do j=Jstr,Jend
            usum(Istr,j)=usum(Istr,j)+qdmu_nbq(Istr,j,k)
          enddo
        enddo
        do k=1,N
          do j=Jstr,Jend
            cff1=0.5*(Hz(Istr,j,k)+Hz(Istr-1,j,k))
            cff2=0.5*(h(Istr  ,j)+zeta(Istr  ,j,knew)
     &               +h(Istr-1,j)+zeta(Istr-1,j,knew))
            qdmu_nbq(Istr,j,k)=qdmu_nbq(Istr,j,k)
     &                         +(DU_nbq(Istr,j)
     &                            -usum(Istr,j))*cff1/cff2
#    ifdef MASKING
     &                                 *umask(Istr,j)
#    endif
#    ifdef WET_DRY
     &                             *umask_wet(Istr,j)
#    endif
          enddo
        enddo
#    undef usum
#   endif /* NBQ */
!
#  else           /* alternative to open */
        do k=1,N                                      ! Western edge closed
          do j=jstr,jend                              ! ======= ==== ======
#   ifdef MRL_WCI
            qdmu_nbq(Istr,j,k)=-ust(Istr,j,k)         ! no Lagrangian flux
     &                          *0.5*(Hz(Istr  ,j,k)
     &                              +Hz(Istr-1,j,k))
#    ifdef MASKING
     &                       *umask(Istr,j)
#    endif
#   else
            qdmu_nbq(Istr,j,k)=0.                     !  no-flux: default
#   endif
          enddo
        enddo

#  endif              /* OBC_COM_WEST */
!
      endif         !<-- WESTERN_EDGE
!
!====================================================================
!                            EASTERN BC
!====================================================================
      if (EASTERN_EDGE) then
#  ifdef OBC_COM_EAST
#   if defined OBC_COM_NBQSPECIFIED || defined OBC_COM_NBQSPECIFIED_EAST
!                                            Eastern edge Specified BC
!                                            ======= ==== ========= ==
        do k=1,N
          do j=Jstr,Jend
#    ifdef NBQ_FRC_BRY
            qdmu_nbq(Iend+1,j,k)=unbqbry_east(j,k)       ! specified
#    else
            qdmu_nbq(Iend+1,j,k)=unbqclm(Iend+1,j,k)
#    endif
#    ifdef MASKING
     &                              *umask(Iend+1,j)
#    endif
#    ifdef WET_DRY
     &                          *umask_wet(Iend+1,j)
#    endif
          enddo
        enddo
!
#   elif defined OBC_COM_NBQORLANSKI
        do k=1,N                             ! Eastern edge radiation
          do j=Jstr,Jend+1                   ! ======= ==== =========
            grad(Iend  ,j)=(qdmu_nbq_east(j  ,k,2)
     &                     -qdmu_nbq_east(j-1,k,2))
#    ifdef MASKING
     &                             *pmask(Iend  ,j)
#    endif
            grad(Iend+1,j)=(qdmu_nbq_east(j  ,k,1)
     &                     -qdmu_nbq_east(j-1,k,1))
#    ifdef MASKING
     &                             *pmask(Iend+1,j)
#    endif
          enddo
          do j=Jstr,Jend
            dft=qdmu_nbq_east(j,k,2)-qdmu_nbq(Iend  ,j,k)
            dfx=qdmu_nbq(Iend  ,j,k)-qdmu_nbq(Iend-1,j,k)
            if (dfx*dft .lt. 0.) then
              dft=0.                      ! <-- cancel cx, if inflow
#    if defined NBQ_FRC_BRY || defined NBQ_NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
            if (dft*(grad(Iend,j)+grad(Iend,j+1)) .gt. 0.) then
              dfy=grad(Iend,j)
            else
              dfy=grad(Iend,j+1)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
            cy=min(cff,max(dft*dfy,-cff))
            qdmu_nbq(Iend+1,j,k)=(    cff*qdmu_nbq_east(j,k,1)
     &                                  +cx*qdmu_nbq(Iend,j,k)
     &                            -max(cy,0.)*grad(Iend+1,j  )
     &                            -min(cy,0.)*grad(Iend+1,j+1)
     &                                              )/(cff+cx)
#    if defined NBQ_FRC_BRY  || defined NBQ_NUDGING 
            qdmu_nbq(Iend+1,j,k)=(1.-tau)*qdmu_nbq(Iend+1,j,k)
#     ifdef NBQ_FRC_BRY
     &                                  +tau*unbqbry_east(j,k)
#     else     
     &                                +tau*unbqclm(Iend+1,j,k)
#     endif
#    endif
#    ifdef MASKING
     &                                        *umask(Iend+1,j)
#    endif
#    ifdef WET_DRY
     &                                    *umask_wet(Iend+1,j)
#    endif
          enddo
        enddo
!
#   else
!                                            Eastern edge gradient BC
!                                            ======= ==== ======== ==
        do k=1,N
          do j=Jstr,Jend
            qdmu_nbq(Iend+1,j,k) = qdmu_nbq(Iend,j,k)
#    ifdef MASKING
     &                               *umask(Iend+1,j)
#    endif
#    ifdef WET_DRY
     &                           *umask_wet(Iend+1,j)
#    endif
          enddo
        enddo
#   endif /* OBC_COM_NBQORLANSKI */
!
!  Replace external mode
!
#   ifdef NBQ
#    define usum grad
        usum=0.
        do k=1,N
          do j=Jstr,Jend
            usum(Iend+1,j)=usum(Iend+1,j)+qdmu_nbq(Iend+1,j,k)
          enddo
        enddo
        do k=1,N
          do j=Jstr,Jend
            cff1=0.5*(Hz(Iend+1,j,k)+Hz(Iend,j,k))
            cff2=0.5*(h(Iend+1,j)+zeta(Iend+1,j,knew)
     &               +h(Iend  ,j)+zeta(Iend  ,j,knew))
            qdmu_nbq(Iend+1,j,k)=qdmu_nbq(Iend+1,j,k)
     &                           +(DU_nbq(Iend+1,j)
     &                              -usum(Iend+1,j))*cff1/cff2
#    ifdef MASKING
     &                               *umask(Iend+1,j)
#    endif
#    ifdef WET_DRY
     &                           *umask_wet(Iend+1,j)
#    endif
          enddo
        enddo
#    undef vsum
#   endif /* NBQ */
!
#  else
        do k=1,N                                     ! Eastern edge closed
          do j=jstr,jend                             ! ======= ==== ======
#   ifdef MRL_WCI
             qdmu_nbq(Iend+1,j,k)=-ust(Iend+1,j,k)   ! no Lagrangian flux
     &                            *0.5*(Hz(Iend+1,j,k)
     &                                 +Hz(Iend  ,j,k))
#    ifdef MASKING
     &                          *umask(Iend+1,j)
#    endif
#   else
             qdmu_nbq(Iend+1,j,k)=0.                 ! no-flux: default
#   endif
          enddo
        enddo

#  endif            /* OBC_COM_EAST */

      endif         !<-- EASTERN_EDGE

# endif             /* !EW_COM_PERIODIC */

# ifndef NS_COM_PERIODIC
!
!====================================================================
!                            SOUTHERN BC
!====================================================================
      if (SOUTHERN_EDGE) then
#  ifdef OBC_COM_SOUTH
#   if defined OBC_COM_NBQSPECIFIED || defined OBC_COM_NBQSPECIFIED_SOUTH
!                                           Southern edge Specified BC
!                                           ======== ==== ========= ==
        do k=1,N
          do i=IstrU,Iend
#    ifdef NBQ_FRC_BRY
            qdmu_nbq(i,Jstr-1,k)=unbqbry_south(i,k)      ! specified
#    else
            qdmu_nbq(i,Jstr-1,k)=unbqclm(i,Jstr-1,k)
#    endif
#    ifdef MASKING
     &                              *umask(i,Jstr-1)
#    endif
#    ifdef WET_DRY
     &                          *umask_wet(i,Jstr-1)
#    endif
          enddo
        enddo
!
#   elif defined OBC_COM_NBQORLANSKI
        do k=1,N                            ! Southern edge radiation
          do i=IstrU-1,Iend                 ! ======== ==== =========
            grad(i,Jstr-1)=qdmu_nbq_south(i+1,k,1)
     &                    -qdmu_nbq_south(i  ,k,1)
            grad(i,Jstr  )=qdmu_nbq_south(i+1,k,2)
     &                    -qdmu_nbq_south(i  ,k,2)
          enddo
          do i=IstrU,Iend
            dft=qdmu_nbq_south(i,k,2)-qdmu_nbq(i,Jstr  ,k)
            dfx=qdmu_nbq(i,Jstr  ,k) -qdmu_nbq(i,Jstr+1,k)
            if (dfx*dft .lt. 0.) then
              dft=0.                       ! <-- cancel cx, if inflow
#    if defined NBQ_FRC_BRY || defined NBQ_NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
            if (dft*(grad(i-1,Jstr)+grad(i,Jstr)) .gt. 0.) then
              dfy=grad(i-1,Jstr)
            else
              dfy=grad(i  ,Jstr)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
            cy=min(cff,max(dft*dfy,-cff))
            qdmu_nbq(i,Jstr-1,k)=(   cff*qdmu_nbq_south(i,k,1)
     &                                  +cx*qdmu_nbq(i,Jstr,k)
     &                            -max(cy,0.)*grad(i-1,Jstr-1)
     &                            -min(cy,0.)*grad(i  ,Jstr-1)
     &                                              )/(cff+cx)
#    if defined NBQ_FRC_BRY  || defined NBQ_NUDGING 
            qdmu_nbq(i,Jstr-1,k)=(1.-tau)*qdmu_nbq(i,Jstr-1,k)
#     ifdef NBQ_FRC_BRY
     &                                 +tau*unbqbry_south(i,k)
#     else     
     &                                +tau*unbqclm(i,Jstr-1,k)
#     endif
#    endif
#    ifdef MASKING
     &                                        *umask(i,Jstr-1)
#    endif
#    ifdef WET_DRY
     &                                    *umask_wet(i,Jstr-1)
#    endif
          enddo
        enddo
!
#   else
!                                           Southern edge gradient BC
!                                           ======== ==== ======== ==
        do k=1,N
          do i=IstrU,Iend
            qdmu_nbq(i,Jstr-1,k)=qdmu_nbq(i,Jstr,k)
#    ifdef MASKING
     &                             *umask(i,Jstr-1)
#    endif
#    ifdef WET_DRY
     &                         *umask_wet(i,Jstr-1)
#    endif
          enddo
        enddo

#   endif  /* OBC_COM_NBQORLANSKI */

!
!  Replace external mode
!
#   if defined NBQ && defined QDM_OBC_TANG_CORRECT
#    define usum grad
        usum=0.
        do k=1,N
          do i=IstrU,Iend
            usum(i,Jstr-1)=usum(i,Jstr-1)+qdmu_nbq(i,Jstr-1,k)
          enddo
        enddo
        do k=1,N
          do i=IstrU,Iend
            cff1=0.5*(Hz(i,Jstr-1,k)+Hz(i-1,Jstr-1,k))
            cff2=0.5*(h(i  ,Jstr-1)+zeta(i  ,Jstr-1,knew)
     &               +h(i-1,Jstr-1)+zeta(i-1,Jstr-1,knew))
            qdmu_nbq(i,Jstr-1,k)=qdmu_nbq(i,Jstr-1,k)
     &                           +(DU_nbq(i,Jstr-1)
     &                              -usum(i,Jstr-1))*cff1/cff2
#    ifdef MASKING
     &                             *umask(i,Jstr-1)
#    endif
#    ifdef WET_DRY
     &                         *umask_wet(i,Jstr-1)
#    endif
          enddo
        enddo
#    undef usum
#   endif

#  else
#   ifdef EW_COM_PERIODIC
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        do k=1,N                        ! Wall: free-slip (gamma2=+1)
          do i=I_RANGE                  ! =====   no-slip (gamma2=-1)
            qdmu_nbq(i,Jstr-1,k)=gamma2*qdmu_nbq(i,Jstr,k)
#   ifdef MASKING
     &                                    *umask(i,Jstr-1)
#   endif
          enddo
        enddo
#   undef I_RANGE

#  endif /* OBC_COM_SOUTH */

      endif              !<-- SOUTHERN_EDGE
!
!====================================================================
!                            NORTHERN BC
!====================================================================
      if (NORTHERN_EDGE) then
#  ifdef OBC_COM_NORTH
#   if defined OBC_COM_NBQSPECIFIED || defined OBC_COM_NBQSPECIFIED_NORTH
!                                           Northern edge Specified BC
!                                           ======== ==== ========= ==
        do k=1,N
          do i=IstrU,Iend
#    ifdef NBQ_FRC_BRY
            qdmu_nbq(i,Jend+1,k)=unbqbry_north(i,k)      ! specified
#    else
            qdmu_nbq(i,Jend+1,k)=unbqclm(i,Jend+1,k)
#    endif
#    ifdef MASKING
     &                          *umask(i,Jend+1)
#    endif
#    ifdef WET_DRY
     &                      *umask_wet(i,Jend+1)
#    endif
          enddo
        enddo
!
#   elif defined OBC_COM_NBQORLANSKI
        do k=1,N                            ! Northern edge radiation
          do i=IstrU-1,Iend                 ! ======== ==== =========
            grad(i,Jend+1)=qdmu_nbq_north(i+1,k,2)
     &                    -qdmu_nbq_north(i  ,k,2)
            grad(i,Jend  )=qdmu_nbq_north(i+1,k,1)
     &                    -qdmu_nbq_north(i  ,k,1)
          enddo
          do i=IstrU,Iend
            dft=qdmu_nbq_north(i,k,2)-qdmu_nbq(i,Jend  ,k)
            dfx=qdmu_nbq(i,Jend  ,k) -qdmu_nbq(i,Jend-1,k)
            if (dfx*dft .lt. 0.) then
              dft=0.                                 ! <-- INFLOW
#    if defined NBQ_FRC_BRY || defined NBQ_NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
            if (dft*(grad(i-1,Jend)+grad(i,Jend)) .gt. 0.) then
              dfy=grad(i-1,Jend)
            else
              dfy=grad(i  ,Jend)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
            cy=min(cff,max(dft*dfy,-cff))
            qdmu_nbq(i,Jend+1,k)=(   cff*qdmu_nbq_north(i,k,1)
     &                                  +cx*qdmu_nbq(i,Jend,k)
     &                            -max(cy,0.)*grad(i-1,Jend+1)
     &                            -min(cy,0.)*grad(i  ,Jend+1)
     &                                              )/(cff+cx)
#    if defined NBQ_FRC_BRY  || defined NBQ_NUDGING 
            qdmu_nbq(i,Jend+1,k)=(1.-tau)*qdmu_nbq(i,Jend+1,k)
#     ifdef NBQ_FRC_BRY
     &                                 +tau*unbqbry_north(i,k)
#     else     
     &                                +tau*unbqclm(i,Jend+1,k)
#     endif
#    endif
#    ifdef MASKING
     &                                       *umask(i,Jend+1)
#    endif
#    ifdef WET_DRY
     &                                   *umask_wet(i,Jend+1)
#    endif
          enddo
        enddo
!
#   else
!                                           Northern edge gradient BC
!                                           ======== ==== ======== ==
        do k=1,N
          do i=IstrU,Iend
            qdmu_nbq(i,Jend+1,k)=qdmu_nbq(i,Jend,k)
#    ifdef MASKING
     &                             *umask(i,Jend+1)
#    endif
#    ifdef WET_DRY
     &                         *umask_wet(i,Jend+1)
#    endif
          enddo
        enddo

#   endif /* OBC_COM_NBQORLANSKI */

!
!  Replace external mode
!
#  if defined NBQ && defined QDM_OBC_TANG_CORRECT
#   define usum grad
        usum=0.
        do k=1,N
          do i=IstrU,Iend
            usum(i,Jend+1)=usum(i,Jend+1)+qdmu_nbq(i,Jend+1,k)
          enddo
        enddo
        do k=1,N
          do i=IstrU,Iend
            cff1=0.5*(Hz(i,Jend+1,k)+Hz(i-1,Jend+1,k))
            cff2=0.5*(h(i  ,Jend+1)+zeta(i  ,Jend+1,knew)
     &               +h(i-1,Jstr-1)+zeta(i-1,Jend+1,knew))
            qdmu_nbq(i,Jend+1,k)=qdmu_nbq(i,Jend+1,k)
     &                           +(DU_nbq(i,Jend+1)
     &                              -usum(i,Jend+1))*cff1/cff2
#    ifdef MASKING
     &                             *umask(i,Jend+1)
#    endif
#    ifdef WET_DRY
     &                         *umask_wet(i,Jend+1)
#    endif
          enddo
        enddo
#   undef usum
#  endif

#  else
#   ifdef EW_COM_PERIODIC
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        do k=1,N                        ! Wall: free-slip (gamma2=+1)
          do i=I_RANGE                  ! =====   no-slip (gamma2=-1)
            qdmu_nbq(i,Jend+1,k)=gamma2*qdmu_nbq(i,Jend,k)
#   ifdef MASKING
     &                                    *umask(i,Jend+1)
#   endif
          enddo
        enddo
#   undef I_RANGE

#  endif  /* OBC_COM_NORTH */

      endif   !<-- NORTHERN_EDGE

# endif        /* !NS_COM_PERIODIC */

!====================================================================
                           ! Corners between adjacent open boundaries
                           ! ======= ======= ======== ==== ==========
 
# if defined OBC_COM_SOUTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          qdmu_nbq(Istr,Jstr-1,k)=0.5*(qdmu_nbq(Istr,Jstr,k)
     &                                +qdmu_nbq(Istr+1,Jstr-1,k))
#  ifdef MASKING
     &                                  *umask(Istr,Jstr-1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_SOUTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          qdmu_nbq(Iend+1,Jstr-1,k)=0.5*(qdmu_nbq(Iend+1,Jstr,k)
     &                                  +qdmu_nbq(Iend,Jstr-1,k))
#  ifdef MASKING
     &                                  *umask(Iend+1,Jstr-1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          qdmu_nbq(Istr,Jend+1,k)=0.5*(qdmu_nbq(Istr+1,Jend+1,k)
     &                                  +qdmu_nbq(Istr,Jend,k))
#  ifdef MASKING
     &                                  *umask(Istr,Jend+1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          qdmu_nbq(Iend+1,Jend+1,k)=0.5*(qdmu_nbq(Iend,Jend+1,k)
     &                                  +qdmu_nbq(Iend+1,Jend,k)) 
#  ifdef MASKING
     &                                  *umask(Iend+1,Jend+1)
#  endif
        enddo
      endif
# endif

      return
      end
#else
# ifndef CHILD
      subroutine unbq_bc_parent_empty
      end
# else
      subroutine unbq_bc_child_empty
      end
# endif
#endif /* M3FAST */
#ifndef CHILD
# define CHILD
# ifdef AGRIF
#  include "unbq_bc.F"
# endif
# undef CHILD
#endif  /* !CHILD */
