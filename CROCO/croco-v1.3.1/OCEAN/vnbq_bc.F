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
      subroutine vnbq_bc_tile (Istr,Iend,Jstr,Jend,grad)
#  ifdef AGRIF      
      use AGRIF_Util
      integer Istr,Iend,Jstr,Jend
      real grad(PRIVATE_2D_SCRATCH_ARRAY)
      if (AGRIF_Root()) then
        call vnbq_bc_parent_tile (Istr,Iend,Jstr,Jend,grad)
      else
        call vnbq_bc_child_tile (Istr,Iend,Jstr,Jend,grad)
!        call vnbq_bc_interp_tile(Istr,Iend,Jstr,Jend)
      endif
      return
      end
!
! PARENT
!
      subroutine vnbq_bc_parent_tile (Istr,Iend,Jstr,Jend,grad)
#  endif
!
! Set lateral boundary conditions for V-component momentum
! qdmv_nbq for the parent grid.
!
# endif /* M3FAST */
#else
# ifdef M3FAST
!
! CHILD
!
      subroutine vnbq_bc_child_tile (Istr,Iend,Jstr,Jend,grad)
!
! Set lateral boundary conditions for V-component momentum
! qdmv_nbq for the child grid.
!
# endif /* M3FAST */
#endif /* CHILD */
# ifdef M3FAST
!
! Common Code
!
# include "set_obc_definitions.h"
!
!
!
      implicit none
# include "param.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "nbq.h"
# include "climat.h"
# include "scalars.h"
# include "boundary.h"
# ifdef MRL_WCI
#  include "forces.h"
# endif
      integer Istr,Iend,Jstr,Jend, i,j,k
      real    grad(PRIVATE_2D_SCRATCH_ARRAY), cff,eps,cff1,cff2,
     &        cx,cy, dft,dfx,dfy, tau,tau_in,tau_out
      parameter (eps=1.E-20)
!
# include "compute_auxiliary_bounds.h"
!
! Interpolations of the parent values to get vnbqbry_east or vnbqclm
!
# ifdef CHILD
      call vnbq_bc_interp_tile(Istr,Iend,Jstr,Jend)
# endif
!
# if defined NBQ_FRC_BRY || defined NBQ_NUDGING
      tau_in =dtfast*tauM_in
      tau_out=dtfast*tauM_out
# endif
!
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
          do i=Istr,Iend
#    ifdef NBQ_FRC_BRY
            qdmv_nbq(i,Jstr,k)=vnbqbry_south(i,k)        ! specified
#    else
            qdmv_nbq(i,Jstr,k)=vnbqclm(i,Jstr,k)
#    endif
#    ifdef MASKING
     &                            *vmask(i,Jstr)
#    endif
#    ifdef WET_DRY
     &                        *vmask_wet(i,Jstr)
#    endif
          enddo
        enddo
!
#   elif defined OBC_COM_NBQORLANSKI
        do k=1,N                            ! Southern edge radiation
          do i=Istr,Iend+1                  ! ======== ==== =========
            grad(i,Jstr  )=(qdmv_nbq_south(i  ,k,1)
     &                     -qdmv_nbq_south(i-1,k,1))
#    ifdef MASKING
     &                              *pmask(i  ,Jstr)
#    endif
            grad(i,Jstr+1)=(qdmv_nbq_south(i  ,k,2)
     &                     -qdmv_nbq_south(i-1,k,2))
#    ifdef MASKING
     &                              *pmask(i,Jstr+1)
#    endif
          enddo
          do i=Istr,Iend
            dft=qdmv_nbq_south(i,k,2)-qdmv_nbq(i,Jstr+1,k)
            dfx=qdmv_nbq(i,Jstr+1,k) -qdmv_nbq(i,Jstr+2,k)
            if (dfx*dft .lt. 0.) then
              dft=0.                      ! <-- cancel cx, if inflow
#    if defined NBQ_FRC_BRY || defined NBQ_NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
            if (dft*(grad(i,Jstr+1)+grad(i+1,Jstr+1)) .gt. 0.) then
              dfy=grad(i  ,Jstr+1)
            else
              dfy=grad(i+1,Jstr+1)
            endif
            dfy=0.
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
            cy=min(cff,max(dft*dfy,-cff))
            qdmv_nbq(i,Jstr,k)=(      cff*qdmv_nbq_south(i,k,1)
     &                                 +cx*qdmv_nbq(i,Jstr+1,k)
     &                               -max(cy,0.)*grad(i  ,Jstr)
     &                               -min(cy,0.)*grad(i+1,Jstr)
     &                                               )/(cff+cx)
#    if defined NBQ_FRC_BRY || defined NBQ_NUDGING 
            qdmv_nbq(i,Jstr,k)=(1.-tau)*qdmv_nbq(i,Jstr,k)
#     ifdef NBQ_FRC_BRY
     &                                  +tau*vnbqbry_south(i,k)
#     else     
     &                                   +tau*vnbqclm(i,Jstr,k)
#     endif
#    endif
#    ifdef MASKING
     &                                           *vmask(i,Jstr)
#    endif
#    ifdef WET_DRY
     &                                       *vmask_wet(i,Jstr)
#    endif
          enddo
        enddo
!
#   else
!                                           Southern edge gradient BC
!                                           ======== ==== ======== ==
        do k=1,N
          do i=Istr,Iend
            qdmv_nbq(i,Jstr,k)=qdmv_nbq(i,Jstr+1,k)
#    ifdef MASKING
     &                               *vmask(i,Jstr)
#    endif
#    ifdef WET_DRY
     &                           *vmask_wet(i,Jstr)
#    endif
          enddo
        enddo
#   endif /* OBC_COM_NBQORLANSKI */

!
!  Replace external mode
!
#   ifdef NBQ
#    define vsum grad
        vsum=0.
        do k=1,N
          do i=Istr,Iend
            vsum(i,Jstr)=vsum(i,Jstr)+qdmv_nbq(i,Jstr,k)
          enddo
        enddo
        do k=1,N
          do i=Istr,Iend
            cff1=0.5*(Hz(i,Jstr,k)+Hz(i,Jstr-1,k))
            cff2=0.5*(h(i,Jstr  )+zeta(i,Jstr  ,knew)
     &              + h(i,Jstr-1)+zeta(i,Jstr-1,knew))
            qdmv_nbq(i,Jstr,k)=qdmv_nbq(i,Jstr,k)
     &                         +(DV_nbq(i,Jstr)
     &                            -vsum(i,Jstr))*cff1/cff2
#    ifdef MASKING
     &                               *vmask(i,Jstr)
#    endif
#    ifdef WET_DRY
     &                           *vmask_wet(i,Jstr)
#    endif
          enddo
        enddo
#    undef vsum
#   endif /* NBQ */

#  else
        do k=1,N                                     ! Southern edge closed
          do i=Istr,Iend                             ! ======== ==== ======
#   ifdef MRL_WCI
            qdmv_nbq(i,Jstr,k)=-vst(i,Jstr,k)        ! no Lagrangian flux
     &                          *0.5*(Hz(i,Jstr  ,k)
     &                              +Hz(i,Jstr-1,k))
#    ifdef MASKING
     &                         *vmask(i,Jstr)
#    endif
#   else
            qdmv_nbq(i,Jstr,k)=0.                    ! no-flux: default
#   endif

          enddo
        enddo
#  endif              /* OBC_COM_SOUTH */

      endif         !<-- SOUTHERN_EDGE
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
          do i=Istr,Iend
#    ifdef NBQ_FRC_BRY
            qdmv_nbq(i,Jend+1,k)=vnbqbry_north(i,k)      ! specified
#    else
            qdmv_nbq(i,Jend+1,k)=vnbqclm(i,Jend+1,k)
#    endif
#    ifdef MASKING
     &                              *vmask(i,Jend+1)
#    endif
#    ifdef WET_DRY
     &                          *vmask_wet(i,Jend+1)
#    endif
          enddo
        enddo
!
#   elif defined OBC_COM_NBQORLANSKI
        do k=1,N                            ! Northern edge radiation
          do i=Istr,Iend+1                  ! ======== ==== =========
            grad(i,Jend  )=(qdmv_nbq_north(i  ,k,2)
     &                     -qdmv_nbq_north(i-1,k,2))
#    ifdef MASKING
     &                                *pmask(i,Jend)
#    endif
            grad(i,Jend+1)=(qdmv_nbq_north(i  ,k,1)
     &                     -qdmv_nbq_north(i-1,k,1))
#    ifdef MASKING
     &                              *pmask(i,Jend+1)
#    endif
          enddo
          do i=Istr,Iend
            dft=qdmv_nbq_north(i,k,2)-qdmv_nbq(i,Jend  ,k)
            dfx=qdmv_nbq(i,Jend  ,k) -qdmv_nbq(i,Jend-1,k)
            if (dfx*dft .lt. 0.) then
              dft=0.                       ! <-- cancel cx, if inflow
#    if defined NBQ_FRC_BRY || defined NBQ_NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
            if (dft*(grad(i,Jend)+grad(i+1,Jend)) .gt. 0.) then
              dfy=grad(i,Jend)
            else
              dfy=grad(i+1,Jend)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
            cy=min(cff,max(dft*dfy,-cff))
            qdmv_nbq(i,Jend+1,k)=(   cff*qdmv_nbq_north(i,k,1)
     &                                  +cx*qdmv_nbq(i,Jend,k)
     &                            -max(cy,0.)*grad(i  ,Jend+1)
     &                            -min(cy,0.)*grad(i+1,Jend+1)
     &                                              )/(cff+cx)

#    if defined NBQ_FRC_BRY  || defined NBQ_NUDGING 
            qdmv_nbq(i,Jend+1,k)=(1.-tau)*qdmv_nbq(i,Jend+1,k)
#     ifdef NBQ_FRC_BRY
     &                                 +tau*vnbqbry_north(i,k)
#     else     
     &                                +tau*vnbqclm(i,Jend+1,k)
#     endif
#    endif
#    ifdef MASKING
     &                                        *vmask(i,Jend+1)
#    endif
#    ifdef WET_DRY
     &                                    *vmask_wet(i,Jend+1)
#    endif
          enddo
        enddo
!
#   else
        do k=1,N
          do i=Istr,Iend
!                                           Northern edge gradient BC
!                                           ======== ==== ======== ==
            qdmv_nbq(i,Jend+1,k)=qdmv_nbq(i,Jend,k)
#    ifdef MASKING
     &                             *vmask(i,Jend+1)
#    endif
#    ifdef WET_DRY
     &                         *vmask_wet(i,Jend+1)
#    endif
          enddo
        enddo
#   endif /* OBC_COM_NBQORLANSKI */

!
!  Replace external mode
!
#   ifdef NBQ
#    define vsum grad
        vsum=0.
        do k=1,N
          do i=Istr,Iend
            vsum(i,Jend+1)=vsum(i,Jend+1)+qdmv_nbq(i,Jend+1,k)
          enddo
        enddo
        do k=1,N
          do i=Istr,Iend
            cff1=0.5*(Hz(i,Jend+1,k)+ Hz(i,Jend,k))
            cff2=0.5*(h(i,Jend+1)+zeta(i,Jend+1,knew)
     &              + h(i,Jend  )+zeta(i,Jend  ,knew))
            qdmv_nbq(i,Jend+1,k)=qdmv_nbq(i,Jend+1,k)
     &                           +(DV_nbq(i,Jend+1)
     &                              -vsum(i,Jend+1))*cff1/cff2
#    ifdef MASKING
     &                             *vmask(i,Jend+1)
#    endif
#    ifdef WET_DRY
     &                         *vmask_wet(i,Jend+1)
#    endif
          enddo
        enddo
#    undef vsum
#   endif /* NBQ */
!
#  else
        do k=1,N                                    ! Northern edge closed
          do i=Istr,Iend                            ! ======== ==== ======
#   ifdef MRL_WCI
            qdmv_nbq(i,Jend+1,k)=-vst(i,Jend+1,k)   ! no Lagrangian flux
     &                           *0.5*(Hz(i,Jend+1,k)
     &                                +Hz(i,Jend ,k))
#    ifdef MASKING
     &                         *vmask(i,Jend+1)
#    endif
#   else
            qdmv_nbq(i,Jend+1,k)=0.                 ! no-flux: default
#   endif
          enddo
        enddo
#  endif

      endif     !<--  NORTHERN_EDGE
# endif          /* !NS_COM_PERIODIC */
 
# ifndef EW_COM_PERIODIC
!
!====================================================================
!                            WESTERN BC
!====================================================================
      if (WESTERN_EDGE) then
#  ifdef OBC_COM_WEST
#   if defined OBC_COM_NBQSPECIFIED || defined OBC_COM_NBQSPECIFIED_WEST
!                                            Western edge Specified BC
!                                            ======= ==== ========= ==
        do k=1,N
          do j=JstrV,Jend
#    ifdef NBQ_FRC_BRY
            qdmv_nbq(Istr-1,j,k)=vnbqbry_west(j,k)       ! specified
#    else
            qdmv_nbq(Istr-1,j,k)=vnbqclm(Istr-1,j,k)
#    endif
#    ifdef MASKING
     &                              *vmask(Istr-1,j)
#    endif
#    ifdef WET_DRY
     &                          *vmask_wet(Istr-1,j)
#    endif
          enddo
        enddo
!
#   elif defined OBC_COM_NBQORLANSKI
        do k=1,N                             ! Western edge radiation
          do j=JstrV-1,Jend                  ! ======= ==== =========
            grad(Istr-1,j)=qdmv_nbq_west(j+1,k,1)
     &                    -qdmv_nbq_west(j  ,k,1)
            grad(Istr  ,j)=qdmv_nbq_west(j+1,k,2)
     &                    -qdmv_nbq_west(j  ,k,2)
          enddo
          do j=JstrV,Jend
            dft=qdmv_nbq_west(j,k,2)-qdmv_nbq(Istr  ,j,k)
            dfx=qdmv_nbq(Istr  ,j,k)-qdmv_nbq(Istr+1,j,k)
            if (dfx*dft .lt. 0.) then
              dft=0.                       ! <-- cancel cx, if inflow
#    if defined NBQ_FRC_BRY || defined NBQ_NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
            if (dft*(grad(Istr,j-1)+grad(Istr,j)) .gt. 0.) then
              dfy=grad(Istr,j-1)
            else
              dfy=grad(Istr,j  )
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
            cy=min(cff,max(dft*dfy,-cff))
            qdmv_nbq(Istr-1,j,k)=(    cff*qdmv_nbq_west(j,k,1)
     &                                  +cx*qdmv_nbq(Istr,j,k)
     &                            -max(cy,0.)*grad(Istr-1,j-1)
     &                            -min(cy,0.)*grad(Istr-1,j  )
     &                                              )/(cff+cx)
#    if defined NBQ_FRC_BRY || defined NBQ_NUDGING 
            qdmv_nbq(Istr-1,j,k)=(1.-tau)*qdmv_nbq(Istr-1,j,k)
#     ifdef NBQ_FRC_BRY
     &                                  +tau*vnbqbry_west(j,k)
#     else
     &                                +tau*vnbqclm(Istr-1,j,k)
#     endif
#    endif
#    ifdef MASKING
     &                                        *vmask(Istr-1,j)
#    endif
#    ifdef WET_DRY
     &                                    *vmask_wet(Istr-1,j)
#    endif
          enddo
        enddo
!
#   else
!                                            Western edge gradient BC
!                                            ======= ==== ======== ==
        do k=1,N
          do j=JstrV,Jend
            qdmv_nbq(Istr-1,j,k)=qdmv_nbq(Istr,j,k)
#    ifdef MASKING
     &                             *vmask(Istr-1,j)
#    endif
#    ifdef WET_DRY
     &                         *vmask_wet(Istr-1,j)
#    endif
          enddo
        enddo    

#   endif /* OBC_COM_NBQORLANSKI */

!
!  Replace external mode
!
#  if defined NBQ && defined QDM_OBC_TANG_CORRECT
#   define vsum grad
        vsum=0.
        do k=1,N
          do j=JstrV,Jend
            vsum(Istr-1,j)=vsum(Istr-1,j)+qdmv_nbq(Istr-1,j,k)
          enddo
        enddo
        do k=1,N
          do j=JstrV,Jend
            cff1=0.5*(Hz(Istr-1,j,k)+Hz(Istr-1,j-1,k))
            cff2=0.5*(h(Istr-1,j  )+zeta(Istr-1,j  ,knew)
     &               +h(Istr-1,j-1)+zeta(Istr-1,j-1,knew))
            qdmv_nbq(Istr-1,j,k)=qdmv_nbq(Istr-1,j,k)
     &                           +(DV_nbq(Istr-1,j)
     &                              -vsum(Istr-1,j))*cff1/cff2
#    ifdef MASKING
     &                             *vmask(Istr-1,j)
#    endif
#    ifdef WET_DRY
     &                         *vmask_wet(Istr-1,j)
#    endif
          enddo
        enddo
#   undef vsum
#  endif

#  else
#   ifdef NS_COM_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        do k=1,N                        ! Wall: free-slip (gamma2=+1)
          do j=J_RANGE                  ! =====   no-slip (gamma2=-1)
            qdmv_nbq(Istr-1,j,k)=gamma2*qdmv_nbq(Istr,j,k)
#   ifdef MASKING
     &                                    *vmask(Istr-1,j)
#   endif
          enddo
        enddo
#   undef J_RANGE

#  endif /* OBC_WEST */

      endif          !<-- WESTERN_EDGE
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
          do j=JstrV,Jend
#    ifdef NBQ_FRC_BRY
            qdmv_nbq(Iend+1,j,k)=vnbqbry_east(j,k)       ! specified
#    else
            qdmv_nbq(Iend+1,j,k)=vnbqclm(Iend+1,j,k)
#    endif
#    ifdef MASKING
     &                              *vmask(Iend+1,j)
#    endif
#    ifdef WET_DRY
     &                          *vmask_wet(Iend+1,j)
#    endif
          enddo
        enddo
!
#   elif defined OBC_COM_NBQORLANSKI
        do k=1,N                             ! Eastern edge radiation
          do j=JstrV-1,Jend                  ! ======= ==== =========
            grad(Iend+1,j)=qdmv_nbq_east(j+1,k,2)
     &                    -qdmv_nbq_east(j  ,k,2)
            grad(Iend  ,j)=qdmv_nbq_east(j+1,k,1)
     &                    -qdmv_nbq_east(j  ,k,1)
          enddo
          do j=JstrV,Jend
            dft=qdmv_nbq_east(j,k,2)-qdmv_nbq(Iend  ,j,k)
            dfx=qdmv_nbq(Iend  ,j,k)-qdmv_nbq(Iend-1,j,k)
            if (dfx*dft .lt. 0.) then
              dft=0.                       ! <-- cancel cx, if inflow
#    if defined NBQ_FRC_BRY || defined NBQ_NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
            if (dft*(grad(Iend,j-1)+grad(Iend,j)) .gt. 0.) then
              dfy=grad(Iend,j-1)
            else
              dfy=grad(Iend,j  )
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
            cy=min(cff,max(dft*dfy,-cff))
            qdmv_nbq(Iend+1,j,k)=(    cff*qdmv_nbq_east(j,k,1)
     &                                  +cx*qdmv_nbq(Iend,j,k)
     &                            -max(cy,0.)*grad(Iend+1,j-1)
     &                            -min(cy,0.)*grad(Iend+1,j  )
     &                                              )/(cff+cx)
#    if defined NBQ_FRC_BRY  || defined NBQ_NUDGING 
            qdmv_nbq(Iend+1,j,k)=(1.-tau)*qdmv_nbq(Iend+1,j,k)
#     ifdef NBQ_FRC_BRY
     &                                  +tau*vnbqbry_east(j,k)
#     else     
     &                                +tau*vnbqclm(Iend+1,j,k)
#     endif
#    endif
#    ifdef MASKING
     &                                        *vmask(Iend+1,j)
#    endif
#    ifdef WET_DRY
     &                                    *vmask_wet(Iend+1,j)
#    endif
          enddo
        enddo
!
#   else
!                                            Eastern edge gradient BC
!                                            ======= ==== ======== ==
        do k=1,N
          do j=JstrV,Jend
            qdmv_nbq(Iend+1,j,k)=qdmv_nbq(Iend,j,k)
#    ifdef MASKING
     &                                  *vmask(Iend+1,j)
#    endif
#    ifdef WET_DRY
     &                              *vmask_wet(Iend+1,j)
#    endif
          enddo
        enddo      
#   endif /* OBC_COM_NBQORLANSKI */

!
!  Replace external mode
!
#   if defined NBQ && defined QDM_OBC_TANG_CORRECT
#    define vsum grad
        vsum=0.
        do k=1,N
          do j=JstrV,Jend
            vsum(Iend+1,j)=vsum(Iend+1,j)+qdmv_nbq(Iend+1,j,k)
          enddo
        enddo
        do k=1,N
          do j=JstrV,Jend
            cff1=0.5*(Hz(Iend+1,j,k)+Hz(Iend+1,j-1,k))
            cff2=0.5*(h(Iend+1,j  )+zeta(Iend+1,j  ,knew)
     &               +h(Iend+1,j-1)+zeta(Iend+1,j-1,knew))
            qdmv_nbq(Iend+1,j,k)=qdmv_nbq(Iend+1,j,k)
     &                           +(DV_nbq(Iend+1,j)
     &                              -vsum(Iend+1,j))*cff1/cff2
#    ifdef MASKING
     &                                  *vmask(Iend+1,j)
#    endif
#    ifdef WET_DRY
     &                              *vmask_wet(Iend+1,j)
#    endif
          enddo
        enddo
#   undef vsum
#  endif

#  else
#   ifdef NS_COM_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        do k=1,N                        ! Wall: free-slip (gamma2=+1)
          do j=J_RANGE                  ! ====    no-slip (gamma2=-1)
            qdmv_nbq(Iend+1,j,k)=gamma2*qdmv_nbq(Iend,j,k)
#   ifdef MASKING
     &                                    *vmask(Iend+1,j)
#   endif
          enddo
        enddo
#   undef J_RANGE

#  endif /* OBC_EAST */

      endif     !<-- EASTERN_EDGE

# endif          /* !EW_COM_PERIODIC */

!====================================================================
                           ! Corners between adjacent open boundaries
                           ! ======= ======= ======== ==== ==========
      
# if defined OBC_COM_SOUTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          qdmv_nbq(Istr-1,Jstr,k)=0.5*(qdmv_nbq(Istr-1,Jstr+1,k)
     &                                +qdmv_nbq(Istr,Jstr,k))
#  ifdef MASKING
     &                                  *vmask(Istr-1,Jstr)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_SOUTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          qdmv_nbq(Iend+1,Jstr,k)=0.5*(qdmv_nbq(Iend,Jstr+1,k)
     &                                +qdmv_nbq(Iend,Jstr,k))
#  ifdef MASKING
     &                                  *vmask(Iend+1,Jstr)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          qdmv_nbq(Istr-1,Jend+1,k)=0.5*(qdmv_nbq(Istr-1,Jend,k)
     &                                  +qdmv_nbq(Istr  ,Jend+1,k))     
#  ifdef MASKING
     &                                *vmask(Istr-1,Jend+1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          qdmv_nbq(Iend+1,Jend+1,k)=0.5*(qdmv_nbq(Iend+1,Jend,k)
     &                                  +qdmv_nbq(Iend  ,Jend+1,k))
#  ifdef MASKING
     &                                *vmask(Iend+1,Jend+1)
#  endif
        enddo
      endif
# endif

      return
      end
#else
# ifndef CHILD
      subroutine vnbq_bc_parent_empty
      end
# else
      subroutine vnbq_bc_child_empty
      end
# endif
#endif /* M3FAST */
#ifndef CHILD
# define CHILD
# ifdef AGRIF
#  include "vnbq_bc.F"
# endif
# undef CHILD
#endif  /* !CHILD */
 
