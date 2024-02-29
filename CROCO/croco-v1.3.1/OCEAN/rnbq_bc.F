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
# ifdef NBQ
      subroutine rnbq_bc_tile (Istr,Iend,Jstr,Jend,grad)
#  ifdef AGRIF      
      use AGRIF_Util
      integer Istr,Iend,Jstr,Jend
      real grad(PRIVATE_2D_SCRATCH_ARRAY)
      if (AGRIF_Root()) then
        call rnbq_bc_parent_tile (Istr,Iend,Jstr,Jend,grad)
      else
        call rnbq_bc_child_tile (Istr,Iend,Jstr,Jend,grad)
      endif
      return
      end
!
! PARENT
!
      subroutine rnbq_bc_parent_tile (Istr,Iend,Jstr,Jend,grad)
#  endif
!
! Set lateral boundary conditions for nbq density rho_nbq
! for the parent grid.
!
# endif /* NBQ */
#else
# ifdef NBQ
!
! CHILD
!
      subroutine rnbq_bc_child_tile (Istr,Iend,Jstr,Jend,grad)
!
! Set lateral boundary conditions for nbq density
! rho_nbq for the child grid.
!
# endif /* NBQ */
#endif /* CHILD */
# ifdef NBQ
!
! Common Code
!
# include "set_obc_definitions.h"
!
!
      implicit none
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "nbq.h"
# include "climat.h"
# include "scalars.h"
# include "boundary.h"  
# ifdef AGRIF
#  include "zoom.h"
# endif
      integer Istr,Iend,Jstr,Jend, i,j,k
      real    grad(PRIVATE_2D_SCRATCH_ARRAY)
      real    cff,eps,
     &        cx,cy, dft,dfx,dfy, tau,tau_in,tau_out
      parameter (eps=1.E-20)
!
# include "compute_auxiliary_bounds.h"
!
! Interpolations of the parent values to get rnbqbry or rnbqclm
!
# ifdef CHILD
      call rnbq_bc_interp_tile(Istr,Iend,Jstr,Jend)
# endif
!
# if defined NBQ_FRC_BRY || defined NBQ_NUDGING
      tau_in =dtfast*tauM_in
      tau_out=dtfast*tauM_out
# endif

# undef OBC_COM_NBQORLANSKI

# ifndef EW_COM_PERIODIC
!
!====================================================================
!                            WESTERN BC
!====================================================================
      if (WESTERN_EDGE) then
#  ifdef OBC_COM_WEST
#   ifdef OBC_COM_NBQORLANSKI
        do k=1,N                             ! Western edge radiation
          do j=Jstr,Jend+1                   ! ======= ==== =========
            grad(Istr-1,j)=(rho_nbq_west(j  ,k,1)
     &                     -rho_nbq_west(j-1,k,1))
#    ifdef MASKING
     &                            *vmask(Istr-1,j)
#    endif
            grad(Istr  ,j)=(rho_nbq_west(j  ,k,2)
     &                     -rho_nbq_west(j-1,k,2))
#    ifdef MASKING
     &                              *vmask(Istr,j)
#    endif
          enddo
          do j=Jstr,Jend
            dft=rho_nbq_west(j,k,2)-rho_nbq(Istr  ,j,k)
            dfx=rho_nbq(Istr  ,j,k)-rho_nbq(Istr+1,j,k)

            if (dfx*dft .lt. 0.) then
              dft=0.                       ! <-- cancel cx, if inflow
#    if defined NBQ_FRC_BRY || defined NBQ_NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
 
            if (dft*(grad(Istr,j)+grad(Istr,j+1)) .gt. 0.) then
              dfy=grad(Istr,j)
            else
              dfy=grad(Istr,j+1)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
            cy=min(cff,max(dft*dfy,-cff))
 
            rho_nbq(Istr-1,j,k)=(    cff*rho_nbq_west(j,k,1)
     &                               +cx*rho_nbq(Istr  ,j,k)
     &                          -max(cy,0.)*grad(Istr-1,j  )
     &                          -min(cy,0.)*grad(Istr-1,j+1)
     &                                            )/(cff+cx)

#    if defined NBQ_FRC_BRY  || defined NBQ_NUDGING 
            rho_nbq(Istr-1,j,k)=(1.-tau)*rho_nbq(Istr-1,j,k)
#     ifdef NBQ_FRC_BRY
     &                                +tau*rnbqbry_west(j,k)
#     else     
     &                              +tau*rnbqclm(Istr-1,j,k)
#     endif
#    endif
#    ifdef MASKING
     &                                      *rmask(Istr-1,j)
#    endif
#    ifdef WET_DRY
     &                                  *rmask_wet(Istr-1,j)
#    endif
          enddo
        enddo
!
#   elif defined OBC_COM_NBQSPECIFIED
!                                            Western edge Specified BC
!                                            ======= ==== ========= ==
        do k=1,N
          do j=Jstr,Jend
#    ifdef NBQ_FRC_BRY
            rho_nbq(Istr-1,j,k)=0.5*(rnbqbry_west(j,k)
     &                              +rho_nbq(Istr,j,k)) ! specified
#    else
            rho_nbq(Istr-1,j,k)=rnbqclm(Istr-1,j,k)
#    endif
#    ifdef MASKING
     &                             *rmask(Istr-1,j)
#    endif
#    ifdef WET_DRY
     &                         *rmask_wet(Istr-1,j)
#    endif
          enddo
        enddo
#   else
!                                            Western edge gradient BC
!                                            ======= ==== ======== ==
        do k=1,N
          do j=Jstr,Jend
           rho_nbq(Istr-1,j,k)=rho_nbq(Istr,j,k)
#    ifdef MASKING
     &                          *rmask(Istr-1,j)
#    endif
#    ifdef WET_DRY
     &                      *rmask_wet(Istr-1,j)
#    endif
          enddo
        enddo

#   endif /* OBC_COM_NBQORLANSKI */
!
#  else           /* alternative to open */
        do k=1,N                                     ! Western edge closed
          do j=jstr,jend                             ! ======= ==== ======
                                                     !  (no-slip: default)
            rho_nbq(Istr-1,j,k)=0.           
          enddo
        enddo

#  endif              /* OBC_COM_WEST */

      endif         !<-- WESTERN_EDGE
!
!====================================================================
!                            EASTERN BC
!====================================================================
      if (EASTERN_EDGE) then
#  ifdef OBC_COM_EAST
#   ifdef OBC_COM_NBQORLANSKI
        do k=1,N                             ! Eastern edge radiation
          do j=Jstr,Jend+1                   ! ======= ==== =========
            grad(Iend  ,j)=(rho_nbq_east(j  ,k,2)
     &                     -rho_nbq_east(j-1,k,2))
#    ifdef MASKING
     &                            *vmask(Iend  ,j)
#    endif
            grad(Iend+1,j)=(rho_nbq_east(j  ,k,1)
     &                     -rho_nbq_east(j-1,k,1))
#    ifdef MASKING
     &                            *vmask(Iend+1,j)
#    endif
          enddo
          do j=Jstr,Jend
            dft=rho_nbq_east(j,k,2)-rho_nbq(Iend  ,j,k)
            dfx=rho_nbq(Iend  ,j,k)-rho_nbq(Iend-1,j,k)
 
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
 
            rho_nbq(Iend+1,j,k)=(    cff*rho_nbq_east(j,k,1)
     &                               +cx*rho_nbq(Iend  ,j,k)
     &                          -max(cy,0.)*grad(Iend+1,j  )
     &                          -min(cy,0.)*grad(Iend+1,j+1)
     &                                            )/(cff+cx)
#    if defined NBQ_FRC_BRY  || defined NBQ_NUDGING 
            rho_nbq(Iend+1,j,k)=(1.-tau)*rho_nbq(Iend+1,j,k)
#     ifdef NBQ_FRC_BRY
     &                                +tau*rnbqbry_east(j,k)
#     else     
     &                              +tau*rnbqclm(Iend+1,j,k)
#     endif
#    endif
#    ifdef MASKING
     &                                      *rmask(Iend+1,j)
#    endif
#    ifdef WET_DRY
     &                                  *rmask_wet(Iend+1,j)
#    endif
          enddo
        enddo
!
#   elif defined OBC_COM_NBQSPECIFIED
!                                            Eastern edge Specified BC
!                                            ======= ==== ========= ==
        do k=1,N
          do j=Jstr,Jend
#    ifdef NBQ_FRC_BRY
             rho_nbq(Iend+1,j,k)=0.5*(rnbqbry_east(j,k) 
     &                               +rho_nbq(Iend,j,k)) ! specified
#    else
             rho_nbq(Iend+1,j,k)=rnbqclm(Iend+1,j,k)
#    endif
#    ifdef MASKING
     &                              *rmask(Iend+1,j)
#    endif
#    ifdef WET_DRY
     &                          *rmask_wet(Iend+1,j)
#    endif
          enddo
        enddo
#   else
!                                            Eastern edge gradient BC
!                                            ======= ==== ======== ==
        do k=1,N
          do j=Jstr,Jend
           rho_nbq(Iend+1,j,k)=rho_nbq(Iend,j,k)
#    ifdef MASKING
     &                          *rmask(Iend+1,j)
#    endif
#    ifdef WET_DRY
     &                      *rmask_wet(Iend+1,j)
#    endif
          enddo
        enddo

#   endif /* OBC_COM_NBQORLANSKI */
!
#  else
        do k=1,N                                    ! Eastern edge closed
          do j=jstr,jend                            ! ======= ==== ======
                                                    !  (no-slip: default)
             rho_nbq(Iend+1,j,k)=0.
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
#   ifdef OBC_COM_NBQORLANSKI
        do k=1,N                            ! Southern edge radiation
          do i=Istr,Iend+1                  ! ======== ==== =========
            grad(i,Jstr-1)=(rho_nbq_south(i  ,k,2)
     &                     -rho_nbq_south(i-1,k,2))
#    ifdef MASKING
     &                             *pmask(i  ,Jstr)
#    endif
            grad(i,Jstr  )=(rho_nbq_south(i  ,k,1)
     &                     -rho_nbq_south(i-1,k,1))
#    ifdef MASKING
     &                             *pmask(i,Jstr+1)
#    endif
          enddo
          do i=Istr,Iend
            dft=rho_nbq_south(i,k,2)-rho_nbq(i,Jstr  ,k)
            dfx=rho_nbq(i,Jstr  ,k) -rho_nbq(i,Jstr+1,k)

            if (dfx*dft .lt. 0.) then
              dft=0.                      ! <-- cancel cx, if inflow
#    if defined NBQ_FRC_BRY || defined NBQ_NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
 
            if (dft*(grad(i,Jstr)+grad(i+1,Jstr)) .gt. 0.) then
              dfy=grad(i,Jstr)
            else
              dfy=grad(i+1,Jstr)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
            cy=min(cff,max(dft*dfy,-cff))
 
            rho_nbq(i,Jstr-1,k)=(  cff*rho_nbq_south(i,k,1)
     &                                +cx*rho_nbq(i,Jstr,k)
     &                         -max(cy,0.)*grad(i  ,Jstr-1)
     &                         -min(cy,0.)*grad(i+1,Jstr-1)
     &                                           )/(cff+cx)

#    if defined NBQ_FRC_BRY  || defined NBQ_NUDGING 
            rho_nbq(i,Jstr-1,k)=(1.-tau)*rho_nbq(i,Jstr-1,k)
#     ifdef NBQ_FRC_BRY
     &                               +tau*rnbqbry_south(i,k)
#     else     
     &                              +tau*rnbqclm(i,Jstr-1,k)
#     endif
#    endif
#    ifdef MASKING
     &                                      *rmask(i,Jstr-1)
#    endif
#    ifdef WET_DRY
     &                                  *rmask_wet(i,Jstr-1)
#    endif
          enddo
        enddo
#   elif defined OBC_COM_NBQSPECIFIED
!                                           Southern edge Specified BC
!                                           ======== ==== ========= ==
        do k=1,N
          do i=Istr,Iend
#    ifdef NBQ_FRC_BRY
            rho_nbq(i,Jstr-1,k)=0.5*(rnbqbry_south(i,k)
     &                              +rho_nbq(i,Jstr,k))! specified
#    else
            rho_nbq(i,Jstr-1,k)=rnbqclm(i,Jstr-1,k)
#    endif
#    ifdef MASKING
     &                             *rmask(i,Jstr-1)
#    endif
#    ifdef WET_DRY
     &                         *rmask_wet(i,Jstr-1)
#    endif
          enddo
        enddo
#   else
!                                           Southern edge gradient BC
!                                           ======== ==== ======== ==
        do k=1,N
          do i=Istr,Iend
            rho_nbq(i,Jstr-1,k)=rho_nbq(i,Jstr,k)
#    ifdef MASKING
     &                           *rmask(i,Jstr-1)
#    endif
#    ifdef WET_DRY
     &                       *rmask_wet(i,Jstr-1)
#    endif
          enddo
        enddo

#   endif /* OBC_COM_NBQORLANSKI */

#  else
        do k=1,N                                     ! Southern edge closed
          do i=Istr,Iend                             ! ======== ==== ======
                                                     !  (no-slip: default)
            rho_nbq(i,Jstr-1,k)=0.
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
#   ifdef OBC_COM_NBQORLANSKI
        do k=1,N                            ! Northern edge radiation
          do i=Istr,Iend+1                  ! ======== ==== =========
            grad(i,Jend  )=(rho_nbq_north(i  ,k,2)
     &                     -rho_nbq_north(i-1,k,2))
#    ifdef MASKING
     &                               *umask(i,Jend)
#    endif
            grad(i,Jend+1)=(rho_nbq_north(i  ,k,1)
     &                     -rho_nbq_north(i-1,k,1))
#    ifdef MASKING
     &                             *umask(i,Jend+1)
#    endif
          enddo
          do i=Istr,Iend
            dft=rho_nbq_north(i,k,2)-rho_nbq(i,Jend  ,k)
            dfx=rho_nbq(i,Jend  ,k) -rho_nbq(i,Jend+1,k)

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
 
            rho_nbq(i,Jend+1,k)=(   cff*rho_nbq_north(i,k,1)
     &                                  +cx*rho_nbq(i,Jend,k)
     &                            -max(cy,0.)*grad(i  ,Jend+1)
     &                            -min(cy,0.)*grad(i+1,Jend+1)
     &                                              )/(cff+cx)

#    if defined NBQ_FRC_BRY  || defined NBQ_NUDGING 
            rho_nbq(i,Jend+1,k)=(1.-tau)*rho_nbq(i,Jend+1,k)
#     ifdef NBQ_FRC_BRY
     &                                 +tau*rnbqbry_north(i,k)
#     else     
     &                                +tau*rnbqclm(i,Jend+1,k)
#     endif
#    endif
#    ifdef MASKING
     &                                       *rmask(i,Jend+1)
#    endif
#    ifdef WET_DRY
     &                                   *rmask_wet(i,Jend+1)
#    endif
          enddo
        enddo
!
#   elif defined OBC_COM_NBQSPECIFIED
!                                           Northern edge Specified BC
!                                           ======== ==== ========= ==
        do k=1,N
          do i=Istr,Iend
#    ifdef NBQ_FRC_BRY
            rho_nbq(i,Jend+1,k)=0.5*(rnbqbry_north(i,k)
     &                        +      rho_nbq(i,Jend,k)) ! specified
#    else
            rho_nbq(i,Jend+1,k)=rnbqclm(i,Jend+1,k)
#    endif
#    ifdef MASKING
     &                             *rmask(i,Jend+1)
#    endif
#    ifdef WET_DRY
     &                         *rmask_wet(i,Jend+1)
#    endif
          enddo
        enddo
#   else
        do k=1,N
          do i=Istr,Iend
!                                           Northern edge gradient BC
!                                           ======== ==== ======== ==
            rho_nbq(i,Jend+1,k)=rho_nbq(i,Jend,k)
#    ifdef MASKING
     &                           *rmask(i,Jend+1)
#    endif
#    ifdef WET_DRY
     &                       *rmask_wet(i,Jend+1)
#    endif
          enddo
        enddo

#   endif /* OBC_COM_NBQORLANSKI */

#  else
        do k=1,N                                      ! Northern edge closed
          do i=Istr,Iend                              ! ======== ==== ======
                                                      !   (no-slip: default)
            rho_nbq(i,Jend+1,k)=0.
          enddo
        enddo
#  endif

      endif     !<--  NORTHERN_EDGE

# endif          /* !NS_COM_PERIODIC */

!====================================================================
                           ! Corners between adjacent open boundaries
                           ! ======= ======= ======== ==== ==========
 
# if defined OBC_COM_SOUTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          rho_nbq(Istr-1,Jstr-1,k)=0.5*(rho_nbq(Istr,Jstr-1,k)
     &                                 +rho_nbq(Istr-1,Jstr,k))
#  ifdef MASKING
     &                                    *rmask(Istr-1,Jstr-1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_SOUTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          rho_nbq(Iend+1,Jstr-1,k)=0.5*(rho_nbq(Iend,Jstr-1,k)
     &                                 +rho_nbq(Iend+1,Jstr,k))     
#  ifdef MASKING
     &                                    *rmask(Iend+1,Jstr-1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          rho_nbq(Istr-1,Jend+1,k)=0.5*(rho_nbq(Istr  ,Jend+1,k)
     &                                 +rho_nbq(Istr-1,Jend,k)) 
#  ifdef MASKING
     &                                     *rmask(Istr-1,Jend+1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          rho_nbq(Iend+1,Jend+1,k)=0.5*(rho_nbq(Iend  ,Jend+1,k)
     &                                 +rho_nbq(Iend+1,Jend,k))     
#  ifdef MASKING
     &                                     *rmask(Iend+1,Jend+1)
#  endif
        enddo
      endif
# endif

      return
      end
#else
# ifndef CHILD
      subroutine rnbq_bc_parent_empty
      end
# else
      subroutine rnbq_bc_child_empty
      end
# endif
#endif /* NBQ */
#ifndef CHILD
# define CHILD
# ifdef AGRIF
#  include "rnbq_bc.F"
# endif
# undef CHILD
#endif  /* !CHILD */
