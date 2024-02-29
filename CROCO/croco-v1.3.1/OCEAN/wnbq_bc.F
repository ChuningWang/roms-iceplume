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
      subroutine wnbq_bc_tile (Istr,Iend,Jstr,Jend,grad)
#  ifdef AGRIF      
      use AGRIF_Util
      integer Istr,Iend,Jstr,Jend
      real grad(PRIVATE_2D_SCRATCH_ARRAY)
      if (AGRIF_Root()) then
        call wnbq_bc_parent_tile (Istr,Iend,Jstr,Jend,grad)
      else
        call wnbq_bc_child_tile (Istr,Iend,Jstr,Jend,grad)
      endif
      return
      end
!
! PARENT
!
      subroutine wnbq_bc_parent_tile (Istr,Iend,Jstr,Jend,grad)
#  endif
!
! Set lateral boundary conditions for W-component momentum
! qdmw_nbq for the parent grid.
!
# endif /* NBQ */
#else
# ifdef NBQ
!
! CHILD
!
      subroutine wnbq_bc_child_tile (Istr,Iend,Jstr,Jend,grad)
!
! Set lateral boundary conditions for W-component momentum
! qdmw_nbq for the child grid.
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
! Interpolations of the parent values to get wnbqbry_east or wnbqclm
!
# ifdef CHILD
      call wnbq_bc_interp_tile(Istr,Iend,Jstr,Jend)
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
        do k=0,N                             ! Western edge radiation
          do j=Jstr,Jend+1                   ! ======= ==== =========
            grad(Istr-1,j)=(qdmw_nbq_west(j  ,k,1)
     &                     -qdmw_nbq_west(j-1,k,1))
#    ifdef MASKING
     &                             *vmask(Istr-1,j)
#    endif
            grad(Istr  ,j)=(qdmw_nbq_west(j  ,k,2)
     &                     -qdmw_nbq_west(j-1,k,2))
#    ifdef MASKING
     &                               *vmask(Istr,j)
#    endif
          enddo
          do j=Jstr,Jend
            dft=qdmw_nbq_west(j,k,2)-qdmw_nbq(Istr  ,j,k)
            dfx=qdmw_nbq(Istr  ,j,k)-qdmw_nbq(Istr+1,j,k)

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
 
            qdmw_nbq(Istr-1,j,k)=(    cff*qdmw_nbq_west(j,k,1)
     &                                +cx*qdmw_nbq(Istr  ,j,k)
     &                            -max(cy,0.)*grad(Istr-1,j  )
     &                            -min(cy,0.)*grad(Istr-1,j+1)
     &                                              )/(cff+cx)

#    if defined NBQ_FRC_BRY  || defined NBQ_NUDGING 
            qdmw_nbq(Istr-1,j,k)=(1.-tau)*qdmw_nbq(Istr-1,j,k)
#     ifdef NBQ_FRC_BRY
     &                                  +tau*wnbqbry_west(j,k)
#     else     
     &                                +tau*wnbqclm(Istr-1,j,k)
#     endif
#    endif
#    ifdef MASKING
     &                                        *rmask(Istr-1,j)
#    endif
#    ifdef WET_DRY
     &                                    *rmask_wet(Istr-1,j)
#    endif
          enddo
        enddo
!
#   elif defined OBC_COM_NBQSPECIFIED
!                                            Western edge Specified BC
!                                            ======= ==== ========= ==
        do k=0,N
          do j=Jstr,Jend
#    ifdef NBQ_FRC_BRY
            qdmw_nbq(Istr-1,j,k)=wnbqbry_west(j,k)         ! specified
#    else
            qdmw_nbq(Istr-1,j,k)=wnbqclm(Istr-1,j,k)
#    endif
#    ifdef MASKING
     &                              *rmask(Istr-1,j)
#    endif
#    ifdef WET_DRY
     &                          *rmask_wet(Istr-1,j)
#    endif
          enddo
        enddo
#   else
!                                            Western edge gradient BC
!                                            ======= ==== ======== ==
        do k=0,N
          do j=Jstr,Jend
           qdmw_nbq(Istr-1,j,k)=qdmw_nbq(Istr,j,k)
#    ifdef MASKING
     &                            *rmask(Istr-1,j)
#    endif
#    ifdef WET_DRY
     &                        *rmask_wet(Istr-1,j)
#    endif
          enddo
        enddo

#   endif /* OBC_COM_NBQORLANSKI */
!
#  else           /* alternative to open */
        do k=1,N                                     ! Western edge closed
          do j=jstr,jend                             ! ======= ==== ======
                                                     !  (no flux)
            qdmw_nbq(Istr-1,j,k)=qdmw_nbq(Istr,j,k)
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
        do k=0,N                             ! Eastern edge radiation
          do j=Jstr,Jend+1                   ! ======= ==== =========
            grad(Iend  ,j)=(qdmw_nbq_east(j  ,k,2)
     &                     -qdmw_nbq_east(j-1,k,2))
#    ifdef MASKING
     &                             *vmask(Iend  ,j)
#    endif
            grad(Iend+1,j)=(qdmw_nbq_east(j  ,k,1)
     &                     -qdmw_nbq_east(j-1,k,1))
#    ifdef MASKING
     &                             *vmask(Iend+1,j)
#    endif
          enddo
          do j=Jstr,Jend
            dft=qdmw_nbq_east(j,k,2)-qdmw_nbq(Iend  ,j,k)
            dfx=qdmw_nbq(Iend  ,j,k)-qdmw_nbq(Iend-1,j,k)
 
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
 
            qdmw_nbq(Iend+1,j,k)=(     cff*qdmw_nbq_east(j,k,1)
     &                                 +cx*qdmw_nbq(Iend  ,j,k)
     &                             -max(cy,0.)*grad(Iend+1,j  )
     &                             -min(cy,0.)*grad(Iend+1,j+1)
     &                                               )/(cff+cx)
#    if defined NBQ_FRC_BRY  || defined NBQ_NUDGING 
            qdmw_nbq(Iend+1,j,k)=(1.-tau)*qdmw_nbq(Iend+1,j,k)
#     ifdef NBQ_FRC_BRY
     &                                  +tau*wnbqbry_east(j,k)
#     else     
     &                                +tau*wnbqclm(Iend+1,j,k)
#     endif
#    endif
#    ifdef MASKING
     &                                        *rmask(Iend+1,j)
#    endif
#    ifdef WET_DRY
     &                                    *rmask_wet(Iend+1,j)
#    endif
          enddo
        enddo
!
#   elif defined OBC_COM_NBQSPECIFIED
!                                            Eastern edge Specified BC
!                                            ======= ==== ========= ==
        do k=0,N
          do j=Jstr,Jend
#    ifdef NBQ_FRC_BRY
             qdmw_nbq(Iend+1,j,k)=wnbqbry_east(j,k)       ! specified
#    else
             qdmw_nbq(Iend+1,j,k)=wnbqclm(Iend+1,j,k)
#    endif
#    ifdef MASKING
     &                               *rmask(Iend+1,j)
#    endif
#    ifdef WET_DRY
     &                           *rmask_wet(Iend+1,j)
#    endif
          enddo
        enddo
#   else
!                                            Eastern edge gradient BC
!                                            ======= ==== ======== ==
        do k=0,N
          do j=Jstr,Jend
           qdmw_nbq(Iend+1,j,k)=qdmw_nbq(Iend,j,k)
#    ifdef MASKING
     &                             *rmask(Iend+1,j)
#    endif
#    ifdef WET_DRY
     &                         *rmask_wet(Iend+1,j)
#    endif
          enddo
        enddo

#   endif /* OBC_COM_NBQORLANSKI */
!
#  else
        do k=1,N                                    ! Eastern edge closed
          do j=jstr,jend                            ! ======= ==== ======
                                                    !  (no flux)
             qdmw_nbq(Iend+1,j,k)=qdmw_nbq(Iend,j,k)
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
        do k=0,N                            ! Southern edge radiation
          do i=Istr,Iend+1                  ! ======== ==== =========
            grad(i,Jstr-1)=(qdmw_nbq_south(i  ,k,2)
     &                     -qdmw_nbq_south(i-1,k,2))
#    ifdef MASKING
     &                              *pmask(i  ,Jstr)
#    endif
            grad(i,Jstr  )=(qdmw_nbq_south(i  ,k,1)
     &                     -qdmw_nbq_south(i-1,k,1))
#    ifdef MASKING
     &                              *pmask(i,Jstr+1)
#    endif
          enddo
          do i=Istr,Iend
            dft=qdmw_nbq_south(i,k,2)-qdmw_nbq(i,Jstr  ,k)
            dfx=qdmw_nbq(i,Jstr  ,k) -qdmw_nbq(i,Jstr+1,k)

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
 
            qdmw_nbq(i,Jstr-1,k)=(   cff*qdmw_nbq_south(i,k,1)
     &                                  +cx*qdmw_nbq(i,Jstr,k)
     &                            -max(cy,0.)*grad(i  ,Jstr-1)
     &                            -min(cy,0.)*grad(i+1,Jstr-1)
     &                                              )/(cff+cx)

#    if defined NBQ_FRC_BRY  || defined NBQ_NUDGING 
            qdmw_nbq(i,Jstr-1,k)=(1.-tau)*qdmw_nbq(i,Jstr-1,k)
#     ifdef NBQ_FRC_BRY
     &                                 +tau*wnbqbry_south(i,k)
#     else     
     &                                +tau*wnbqclm(i,Jstr-1,k)
#     endif
#    endif
#    ifdef MASKING
     &                                        *rmask(i,Jstr-1)
#    endif
#    ifdef WET_DRY
     &                                    *rmask_wet(i,Jstr-1)
#    endif
          enddo
        enddo
#   elif defined OBC_COM_NBQSPECIFIED
!                                           Southern edge Specified BC
!                                           ======== ==== ========= ==
        do k=0,N
          do i=Istr,Iend
#    ifdef NBQ_FRC_BRY
            qdmw_nbq(i,Jstr-1,k)=wnbqbry_south(i,k)        ! specified
#    else
            qdmw_nbq(i,Jstr-1,k)=wnbqclm(i,Jstr-1,k)
#    endif
#    ifdef MASKING
     &                              *rmask(i,Jstr-1)
#    endif
#    ifdef WET_DRY
     &                          *rmask_wet(i,Jstr-1)
#    endif
          enddo
        enddo
#   else
!                                           Southern edge gradient BC
!                                           ======== ==== ======== ==
        do k=0,N
          do i=Istr,Iend
            qdmw_nbq(i,Jstr-1,k)=qdmw_nbq(i,Jstr,k)
#    ifdef MASKING
     &                             *rmask(i,Jstr-1)
#    endif
#    ifdef WET_DRY
     &                         *rmask_wet(i,Jstr-1)
#    endif
          enddo
        enddo

#   endif /* OBC_COM_NBQORLANSKI */

#  else
        do k=1,N                                     ! Southern edge closed
          do i=Istr,Iend                             ! ======== ==== ======
                                                     !  (no flux)
            qdmw_nbq(i,Jstr-1,k)=qdmw_nbq(i,Jstr,k)
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
        do k=0,N                            ! Northern edge radiation
          do i=Istr,Iend+1                  ! ======== ==== =========
            grad(i,Jend  )=(qdmw_nbq_north(i  ,k,2)
     &                     -qdmw_nbq_north(i-1,k,2))
#    ifdef MASKING
     &                                *umask(i,Jend)
#    endif
            grad(i,Jend+1)=(qdmw_nbq_north(i  ,k,1)
     &                     -qdmw_nbq_north(i-1,k,1))
#    ifdef MASKING
     &                              *umask(i,Jend+1)
#    endif
          enddo
          do i=Istr,Iend
            dft=qdmw_nbq_north(i,k,2)-qdmw_nbq(i,Jend  ,k)
            dfx=qdmw_nbq(i,Jend  ,k) -qdmw_nbq(i,Jend+1,k)

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
 
            qdmw_nbq(i,Jend+1,k)=(   cff*qdmw_nbq_north(i,k,1)
     &                                  +cx*qdmw_nbq(i,Jend,k)
     &                            -max(cy,0.)*grad(i  ,Jend+1)
     &                            -min(cy,0.)*grad(i+1,Jend+1)
     &                                              )/(cff+cx)

#    if defined NBQ_FRC_BRY  || defined NBQ_NUDGING 
            qdmw_nbq(i,Jend+1,k)=(1.-tau)*qdmw_nbq(i,Jend+1,k)
#     ifdef NBQ_FRC_BRY
     &                                 +tau*wnbqbry_north(i,k)
#     else     
     &                                +tau*wnbqclm(i,Jend+1,k)
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
        do k=0,N
          do i=Istr,Iend
#    ifdef NBQ_FRC_BRY
            qdmw_nbq(i,Jend+1,k)=wnbqbry_north(i,k)      ! specified
#    else
            qdmw_nbq(i,Jend+1,k)=wnbqclm(i,Jend+1,k)
#    endif
#    ifdef MASKING
     &                              *rmask(i,Jend+1)
#    endif
#    ifdef WET_DRY
     &                          *rmask_wet(i,Jend+1)
#    endif
          enddo
        enddo
#   else
        do k=0,N
          do i=Istr,Iend
!                                           Northern edge gradient BC
!                                           ======== ==== ======== ==
            qdmw_nbq(i,Jend+1,k)=qdmw_nbq(i,Jend,k)
#    ifdef MASKING
     &                             *rmask(i,Jend+1)
#    endif
#    ifdef WET_DRY
     &                         *rmask_wet(i,Jend+1)
#    endif
          enddo
        enddo

#   endif /* OBC_COM_NBQORLANSKI */

#  else
        do k=1,N                                      ! Northern edge closed
          do i=Istr,Iend                              ! ======== ==== ======
                                                      !   (no flux)
            qdmw_nbq(i,Jend+1,k)=qdmw_nbq(i,Jend,k)
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
        do k=0,N
          qdmw_nbq(Istr-1,Jstr-1,k)=0.5*(qdmw_nbq(Istr,Jstr-1,k)
     &                                  +qdmw_nbq(Istr-1,Jstr,k))
#  ifdef MASKING
     &                                *rmask(Istr-1,Jstr-1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_SOUTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=0,N
          qdmw_nbq(Iend+1,Jstr-1,k)=0.5*(qdmw_nbq(Iend,Jstr-1,k)
     &                                  +qdmw_nbq(Iend+1,Jstr,k))     
#  ifdef MASKING
     &                                *rmask(Iend+1,Jstr-1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=0,N
          qdmw_nbq(Istr-1,Jend+1,k)=0.5*(qdmw_nbq(Istr  ,Jend+1,k)
     &                                  +qdmw_nbq(Istr-1,Jend,k)) 
#  ifdef MASKING
     &                                *rmask(Istr-1,Jend+1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=0,N
          qdmw_nbq(Iend+1,Jend+1,k)=0.5*(qdmw_nbq(Iend  ,Jend+1,k)
     &                                  +qdmw_nbq(Iend+1,Jend,k))     
#  ifdef MASKING
     &                                *rmask(Iend+1,Jend+1)
#  endif
        enddo
      endif
# endif

      return
      end
#else
# ifndef CHILD
      subroutine wnbq_bc_parent_empty
      end
# else
      subroutine wnbq_bc_child_empty
      end
# endif
#endif /* NBQ */
#ifndef CHILD
# define CHILD
# ifdef AGRIF
#  include "wnbq_bc.F"
# endif
# undef CHILD
#endif  /* !CHILD */
