! $Id: v2dbc.F 1419 2014-01-13 11:16:42Z marchesiello $
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
      subroutine v2dbc_tile(Istr,Iend,Jstr,Jend,grad)
# ifdef AGRIF
      use AGRIF_Util
      integer Istr,Iend,Jstr,Jend
      real grad(PRIVATE_2D_SCRATCH_ARRAY)
      if (AGRIF_Root()) then
        call v2dbc_parent_tile(Istr,Iend,Jstr,Jend,grad)
      else
        call v2dbc_child_tile(Istr,Iend,Jstr,Jend,grad)
      endif
      return
      end
!
! PARENT
!
      subroutine v2dbc_parent_tile(Istr,Iend,Jstr,Jend,grad)
# endif
!
! Set lateral boundary conditions for the barotropic (i.e.
! vertically integrated) ETA-component velocity vbar(:,:,knew)
! for the parent grid.
!
#else
!
! CHILD
!
      subroutine v2dbc_child_tile (Istr,Iend,Jstr,Jend,grad)
!
! Set lateral boundary conditions for the barotropic (i.e.
! vertically integrated) ETA-component velocity vbar(:,:,knew)
! for the child grid.
!
#endif /* CHILD */
!
! Common Code
!
# include "set_obc_definitions.h"
!
      implicit none
#include "param.h"
      integer Istr,Iend,Jstr,Jend, i,j
      real    grad(PRIVATE_2D_SCRATCH_ARRAY)
#ifdef AGRIF
#include "private_scratch.h"
      integer trd
C$    integer omp_get_thread_num
#endif
      real    eps,cff, cx,cy,
     &        dft,dfx,dfy, tau,tau_in,tau_out,hx,zx
#ifdef OBC_REDUCED_PHYSICS
     &        ,bry_val,bry_pgr,bry_str,bry_cor
#endif
      parameter (eps=1.D-20)
#if defined WET_DRY || defined OBC_REDUCED_PHYSICS
      real cff1,cff2
#endif
#include "boundary.h"
#include "climat.h"
#include "grid.h"
#include "ocean2d.h"
#include "scalars.h"
#include "forces.h"
#include "zoom.h"
#ifdef NBQ
# include "nbq.h"
#endif
!
#include "compute_auxiliary_bounds.h"
!
!  If child grid, interpolate vbar from parent
!  AGRIF does not work yet with #undef SOLVE3D
!
#ifdef CHILD
#  ifdef SOLVE3D
#   ifdef NBQ
      if (M2bc_nbq_flag) then
#   endif
      trd=0
C$    trd=omp_get_thread_num()
      call v2dbc_interp_tile(Istr,Iend,Jstr,Jend
     &                    ,A1dETA(1,1)
     &                    ,A1dETA(1,2+NWEIGHT)
     &                    ,A1dETA(1,3+2*NWEIGHT)
     &                    ,A1dETA(1,4+2*NWEIGHT)
     &                    ,A1dXI(1,1)
     &                    ,A1dXI(1,2+NWEIGHT)
     &                    ,A1dXI(1,3+2*NWEIGHT)
     &                    ,A1dXI(1,4+2*NWEIGHT))
#   ifdef NBQ
      endif
#   endif
#  else
      Some coding is needed to have AGRIF without SOLVE3D.....
#  endif
#endif
!
#if defined M2_FRC_BRY || defined M2CLIMATOLOGY
      tau_in=dtfast*tauM_in
      tau_out=dtfast*tauM_out
#endif
!
#if defined OBC_COM_M2CHARACT || defined OBC_COM_M2CHARACT_WEST || defined OBC_COM_M2CHARACT_EAST
      grad = 1.
# if defined OBC_COM_SOUTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        grad(Istr,Jstr) = 0.5
      endif
# endif
# if defined OBC_COM_SOUTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        grad(Iend,Jstr) = 0.5
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        grad(Istr,Jend+1) = 0.5
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        grad(Iend,Jend+1) = 0.5
      endif
# endif
#endif

#ifndef NS_COM_PERIODIC
!
!====================================================================
!                            SOUTHERN BC
!====================================================================
      if (SOUTHERN_EDGE) then
# ifdef OBC_COM_SOUTH
#  ifdef OBC_COM_M2ORLANSKI
!                                          Southern edge radiation BC
!                                          ======== ==== ========= ==
        do i=Istr,Iend+1
          grad(i,Jstr  )=(vbar(i,Jstr  ,kstp)-vbar(i-1,Jstr  ,kstp))
#   ifdef MASKING
     &                                                *pmask(i,Jstr)
#   endif
          grad(i,Jstr+1)=(vbar(i,Jstr+1,kstp)-vbar(i-1,Jstr+1,kstp))
#   ifdef MASKING
     &                                              *pmask(i,Jstr+1)
#   endif
        enddo
        do i=Istr,Iend
          dft=vbar(i,Jstr+1,kstp)-vbar(i,Jstr+1,knew)
          dfx=vbar(i,Jstr+1,knew)-vbar(i,Jstr+2,knew)

          if (dfx*dft .lt. 0.) then
            dft=0.                       ! <-- cancel cx, if inflow
#   if defined M2_FRC_BRY || defined M2CLIMATOLOGY
            tau=tau_in
          else
            tau=tau_out
#   endif
          endif

          if (dft*(grad(i,Jstr+1)+grad(i+1,Jstr+1)) .gt. 0.) then
            dfy=grad(i  ,Jstr+1)
          else
            dfy=grad(i+1,Jstr+1)
          endif

#   ifdef OBC_COM_RAD_NORMAL
          dfy=0.
#   endif
          cff=max(dfx*dfx+dfy*dfy, eps)
          cx=dft*dfx
#   ifdef OBC_COM_RAD_NPO
          cy=0.
#   else
          cy=min(cff,max(dft*dfy,-cff))
#   endif

          vbar(i,Jstr,knew)=( cff*vbar(i,Jstr  ,kstp)
     &                        +cx*vbar(i,Jstr+1,knew)
     &                     -max(cy,0.)*grad(i  ,Jstr)
     &                     -min(cy,0.)*grad(i+1,Jstr)
     &                                     )/(cff+cx)
#   if defined M2_FRC_BRY || defined M2CLIMATOLOGY
          vbar(i,Jstr,knew)=(1.-tau)*vbar(i,Jstr,knew)
#    ifdef M2_FRC_BRY
     &                           +tau*vbarbry_south(i)
#    else
     &                           +tau*vbclm(i,Jstr)
#    endif
#   endif
#   ifdef MASKING
          vbar(i,Jstr,knew)=vbar(i,Jstr,knew)*vmask(i,Jstr)
#   endif
        enddo
#  elif defined OBC_COM_M2CHARACT
!                                            Southern edge Characteristic BC
!                                            ======= ==== ======= ==
        do i=Istr,Iend
          cff=0.5*(h(i,Jstr-1)+zeta(i,Jstr-1,kstp)+
     &             h(i,Jstr  )+zeta(i,Jstr  ,kstp))
          hx=sqrt(g/cff)
          cx=dtfast*cff*hx*0.5*(pn(i,Jstr-1)+pn(i,Jstr))
          zx=(0.5+cx)*zeta(i,jstr,kstp)+(0.5-cx)*zeta(i,jstr-1,kstp)
          if (cx .gt. 0.292893218813452) then
            zx=zx + ( zeta(i,jstr,knew) +cx*zeta(i,jstr-1,kstp)
     &                               -(1.+cx)*zeta(i,jstr,kstp)
     &                           )*(1.-0.292893218813452/cx)**2
          endif
#   ifdef OBC_REDUCED_PHYSICS
          cff1=1.0/(0.5*(h(i,Jstr-1)+zeta(i,Jstr-1,kstp)+
     &                   h(i,Jstr  )+zeta(i,Jstr  ,kstp)))
          bry_str=cff1*(svstr(i,Jstr)-bvstr(i,Jstr))
#    ifdef UV_COR
          bry_cor=-0.125*(ubar(i  ,Jstr-1,kstp)+
     &                    ubar(i+1,Jstr-1,kstp)+
     &                    ubar(i  ,Jstr  ,kstp)+
     &                    ubar(i+1,Jstr  ,kstp))*
     &                    (f(i,Jstr-1)+f(i,Jstr))
#    else
          bry_cor=0.
#    endif
          bry_pgr=-g*(zeta(i,Jstr,kstp)-zeta(i,Jstr-1,kstp))*
     &                          0.5*(pn(i,Jstr-1)+pn(i,Jstr))
          Cy=1.0/SQRT(g*0.5*(h(i,Jstr-1)+zeta(i,Jstr-1,kstp)+
     &                       h(i,Jstr  )+zeta(i,Jstr  ,kstp)))
          cff2=on_v(i,Jstr)*Cy
          bry_val=vbar(i,Jstr+1,kstp)+
     &                cff2*(bry_pgr+
     &                      bry_cor+
     &                      bry_str)
#   endif
          vbar(i,Jstr,knew)= 0.5*( (1.-cx)*vbar(i,Jstr,  kstp)
     &                                 +cx*vbar(i,Jstr+1,kstp)
#   ifdef M2_FRC_BRY
#    ifdef OBC_REDUCED_PHYSICS
     &                                     +bry_val
#    else
     &                                     +vbarbry_south(i)
#    endif
#   elif defined M2CLIMATOLOGY
     &                                     +vbclm(i,Jstr)
#    endif
     &                                     -hx*( zx
#   ifdef Z_FRC_BRY
     &                                        -zetabry_south(i)
#   elif defined ZCLIMATOLOGY
     &                                        -ssh(i,Jstr-1)
#   endif
     &                           ))
#   ifdef MASKING
          vbar(i,Jstr,knew)=vbar(i,Jstr,knew)*vmask(i,Jstr)
#   endif
        enddo
#  elif defined OBC_COM_M2SPECIFIED
!                                           Southern edge Specified BC
!                                           ======== ==== ========= ==
        do i=Istr,Iend
#   ifdef M2_FRC_BRY
          vbar(i,Jstr,knew)=vbarbry_south(i)
#   else
          vbar(i,Jstr,knew)=vbclm(i,Jstr)
#   endif
#   ifdef MASKING
     &                      *vmask(i,Jstr)
#   endif
        enddo
#  else
!                                           Southern edge gradient BC
!                                           ======== ==== ======== ==
        do i=Istr,Iend
          vbar(i,Jstr,knew)=vbar(i,Jstr+1,knew)
#   ifdef MASKING
     &                      *vmask(i,Jstr)
#   endif
        enddo

#  endif /* OBC_COM_ORLANSKI */
# else
                                              ! Southern edge closed
        do i=Istr,Iend                        ! ======== ==== ======
#  ifdef MRL_WCI
          vbar(i,Jstr,knew)=-vst2d(i,Jstr)    ! no Lagrangian flux
#   ifdef MASKING
     &                      *vmask(i,Jstr)
#   endif
#  else
          vbar(i,Jstr,knew)=0.                ! no-flux, default
#  endif
        enddo
# endif          /* OBC_COM_SOUTH */

      endif   !<-- SOUTHERN_EDGE
!
!====================================================================
!                            NORTHERN BC
!====================================================================
      if (NORTHERN_EDGE) then
# ifdef OBC_COM_NORTH
#  ifdef OBC_COM_M2ORLANSKI
!                                          Northern edge radiation BC
!                                          ======== ==== ========= ==
        do i=Istr,Iend+1
          grad(i,Jend  )=(vbar(i,Jend  ,kstp)-vbar(i-1,Jend  ,kstp))
#   ifdef MASKING
     &                                                *pmask(i,Jend)
#   endif
          grad(i,Jend+1)=(vbar(i,Jend+1,kstp)-vbar(i-1,Jend+1,kstp))
#   ifdef MASKING
     &                                              *pmask(i,Jend+1)
#   endif
        enddo
        do i=Istr,Iend
          dft=vbar(i,Jend,kstp)-vbar(i,Jend  ,knew)
          dfx=vbar(i,Jend,knew)-vbar(i,Jend-1,knew)

          if (dfx*dft .lt. 0.) then
            dft=0.                         ! <-- cancel cx, if inflow
#   if defined M2_FRC_BRY || defined M2CLIMATOLOGY
            tau=tau_in
          else
            tau=tau_out
#   endif
          endif

          if (dft*(grad(i,Jend)+grad(i+1,Jend)) .gt. 0.) then
            dfy=grad(i  ,Jend)
          else
            dfy=grad(i+1,Jend)
          endif

#   ifdef OBC_COM_RAD_NORMAL
          dfy=0.
#   endif
          cff=max(dfx*dfx+dfy*dfy, eps)
          cx=dft*dfx
#   ifdef OBC_COM_RAD_NPO
          cy=0.
#   else
          cy=min(cff,max(dft*dfy,-cff))
#   endif

          vbar(i,Jend+1,knew)=( cff*vbar(i,Jend+1,kstp)
     &                          +cx*vbar(i,Jend  ,knew)
     &                     -max(cy,0.)*grad(i  ,Jend+1)
     &                     -min(cy,0.)*grad(i+1,Jend+1)
     &                                       )/(cff+cx)
#   if defined M2_FRC_BRY || defined M2CLIMATOLOGY
          vbar(i,Jend+1,knew)=(1.-tau)*vbar(i,Jend+1,knew)
#    ifdef M2_FRC_BRY
     &                               +tau*vbarbry_north(i)
#    else
     &                               +tau*vbclm(i,Jend+1)
#    endif
#   endif
#   ifdef MASKING
          vbar(i,Jend+1,knew)=vbar(i,Jend+1,knew)*vmask(i,Jend+1)
#   endif
        enddo
#  elif defined OBC_COM_M2CHARACT
!                                            Northern edge Characteristic BC
!                                            ======= ==== ======= ==
        do i=Istr,Iend
          cff=0.5*(h(i,Jend  )+zeta(i,Jend  ,kstp)+
     &             h(i,Jend+1)+zeta(i,Jend+1,kstp))
          hx=sqrt(g/cff)
          cx=dtfast*cff*hx*0.5*(pn(i,Jend)+pn(i,Jend+1))
          zx=(0.5+cx)*zeta(i,jend,kstp)+(0.5-cx)*zeta(i,jend+1,kstp)
          if (cx .gt. 0.292893218813452) then
            zx=zx + ( zeta(i,jend,knew) +cx*zeta(i,jend+1,kstp)
     &                               -(1.+cx)*zeta(i,jend,kstp)
     &                           )*(1.-0.292893218813452/cx)**2
          endif
#   ifdef OBC_REDUCED_PHYSICS
          cff1=1.0/(0.5*(h(i,Jend  )+zeta(i,Jend  ,kstp)+
     &                   h(i,Jend+1)+zeta(i,Jend+1,kstp)))
          bry_str=cff1*(svstr(i,Jend+1)-bvstr(i,Jend+1))
#    ifdef UV_COR
          bry_cor=-0.125*(ubar(i  ,Jend,  kstp)+
     &                    ubar(i+1,Jend,  kstp)+
     &                    ubar(i  ,Jend+1,kstp)+
     &                    ubar(i+1,Jend+1,kstp))*
     &                    (f(i,Jend)+f(i,Jend+1))
#    else
          bry_cor=0.
#    endif
          bry_pgr=-g*(zeta(i,Jend+1,kstp)-zeta(i,Jend,kstp))*
     &                          0.5*(pn(i,Jend)+pn(i,Jend+1))
          Cy=1.0/SQRT(g*0.5*(h(i,Jend  )+zeta(i,Jend,  kstp)+
     &                       h(i,Jend+1)+zeta(i,Jend+1,kstp)))
          cff2=on_v(i,Jend+1)*Cy
          bry_val=vbar(i,Jend,kstp)+
     &                cff2*(bry_pgr+
     &                      bry_cor+
     &                      bry_str)
#   endif
          vbar(i,Jend+1,knew)= 0.5*( (1.-cx)*vbar(i,Jend+1,kstp)
     &                                   +cx*vbar(i,Jend  ,kstp)
#   ifdef M2_FRC_BRY
#    ifdef OBC_REDUCED_PHYSICS
     &                                       +bry_val
#    else
     &                                       +vbarbry_north(i)
#    endif
#   elif defined M2CLIMATOLOGY
     &                                       +vbclm(i,Jend+1)
#    endif
     &                                       +hx*( zx
#   ifdef Z_FRC_BRY
     &                                         -zetabry_north(i)
#   elif defined ZCLIMATOLOGY
     &                                         -ssh(i,Jend)
#   endif
     &                            ))
#   ifdef MASKING
          vbar(i,Jend+1,knew)=vbar(i,Jend+1,knew)*vmask(i,Jend+1)
#   endif
        enddo
#  elif defined OBC_COM_M2SPECIFIED
!                                           Northern edge Specified BC
!                                           ======== ==== ========= ==
        do i=Istr,Iend
#   ifdef M2_FRC_BRY
          vbar(i,Jend+1,knew)=vbarbry_north(i)
#   else
          vbar(i,Jend+1,knew)=vbclm(i,Jend+1)
#   endif
#   ifdef MASKING
     &                        *vmask(i,Jend+1)
#   endif
        enddo
#  else
!                                           Northern edge gradient BC
!                                           ======== ==== ======== ==
        do i=Istr,Iend
          vbar(i,Jend+1,knew)=vbar(i,Jend,knew)
#   ifdef MASKING
     &                        *vmask(i,Jend+1)
#   endif
        enddo

#  endif /* OBC_COM_ORLANSKI */
# else
                                              ! Northern edge closed
        do i=Istr,Iend                        ! ======== ==== ======
#  ifdef MRL_WCI
          vbar(i,Jend+1,knew)=-vst2d(i,Jend+1)  ! no Lagrangian flux
#   ifdef MASKING
     &                        *vmask(i,Jend+1)
#   endif
#  else
          vbar(i,Jend+1,knew)=0.                ! no-flux: default
#  endif
        enddo

# endif /* OBC_COM_NORTH */

      endif          !<--  NORTHERN_EDGE

#endif                  /* !NS_COM_PERIODIC */


#ifndef EW_COM_PERIODIC
!
!====================================================================
!                            WESTERN BC
!====================================================================
      if (WESTERN_EDGE) then
# ifdef OBC_COM_WEST
#  if defined OBC_COM_M2SPECIFIED || defined OBC_COM_M2SPECIFIED_WEST
!                                           Western edge Specified BC
!                                           ======= ==== ========= ==
       do j=JstrV,Jend
#   ifdef M2_FRC_BRY
          vbar(Istr-1,j,knew)=vbarbry_west(j)
#   else
          vbar(Istr-1,j,knew)=vbclm(Istr-1,j)
#   endif
#   ifdef MASKING
     &                        *vmask(Istr-1,j)
#   endif
        enddo
#  elif defined OBC_COM_M2ORLANSKI
!                                           Western edge radiation BC
!                                           ======= ==== ========= ==
        do j=JstrV-1,Jend
          grad(Istr-1,j)=vbar(Istr-1,j+1,kstp)-vbar(Istr-1,j,kstp)
          grad(Istr  ,j)=vbar(Istr  ,j+1,kstp)-vbar(Istr  ,j,kstp)
        enddo
        do j=JstrV,Jend
          dft=vbar(Istr,j,kstp)-vbar(Istr  ,j,knew)
          dfx=vbar(Istr,j,knew)-vbar(Istr+1,j,knew)

          if (dfx*dft .lt. 0.) then
            dft=0.                        ! <-- cancel cx, if inflow
#   if defined M2_FRC_BRY || defined M2CLIMATOLOGY
            tau=tau_in
          else
            tau=tau_out
#   endif
          endif

          if (dft*(grad(Istr,j-1)+grad(Istr,j)) .gt. 0.) then
            dfy=grad(Istr,j-1)
          else
            dfy=grad(Istr,j  )
          endif

#   ifdef OBC_COM_RAD_NORMAL
          dfy=0.
#   endif
          cff=max(dfx*dfx+dfy*dfy, eps)
          cx=dft*dfx
#   ifdef OBC_COM_RAD_NPO
          cy=0.
#   else
          cy=min(cff,max(dft*dfy,-cff))
#   endif

          vbar(Istr-1,j,knew)=( cff*vbar(Istr-1,j,kstp)
     &                          +cx*vbar(Istr  ,j,knew)
     &                     -max(cy,0.)*grad(Istr-1,j-1)
     &                     -min(cy,0.)*grad(Istr-1,j  )
     &                                       )/(cff+cx)
#   if defined M2_FRC_BRY || defined M2CLIMATOLOGY
          vbar(Istr-1,j,knew)=(1.-tau)*vbar(Istr-1,j,knew)
#    ifdef M2_FRC_BRY
     &                                +tau*vbarbry_west(j)
#    else
     &                                +tau*vbclm(Istr-1,j)
#    endif
#   endif
#   ifdef MASKING
          vbar(Istr-1,j,knew)=vbar(Istr-1,j,knew)*vmask(Istr-1,j)
#   endif
        enddo
#  elif defined OBC_COM_M2CHARACT
!                                             Western edge Chapman BC
!                                             ======= ==== ======= ==
        do j=JstrV,Jend
          cff=sqrt(0.5*g*(h(Istr-1,j-1)+zeta(Istr-1,j-1,kstp)+
     &                    h(Istr-1,j  )+zeta(Istr-1,j  ,kstp)))
          cx=dtfast*cff*0.5*(pm(Istr-1,j-1)+pm(Istr-1,j))
          vbar(Istr-1,j,knew)=( vbar(Istr-1,j,kstp)
     &                      +cx*vbar(Istr  ,j,knew) )/(1.+cx)
#   ifdef MASKING
     &                        *vmask(Istr-1,j)
#   endif
        enddo
#  elif defined OBC_COM_M2CHARACT0
!                                      Western edge Characteristic BC
!                                      ======= ==== ============== ==
        do j=JstrV,Jend
          if (ubar(Istr+1,j,knew).lt.0.) then
            vbar(Istr-1,j,knew)=vbar(Istr,j,knew)
          else
#   ifdef M2_FRC_BRY
            vbar(Istr-1,j,knew)=vbarbry_west(j)
#   else
            vbar(Istr-1,j,knew)=vbclm(Istr-1,j)
#   endif
          endif
#   ifdef MASKING
          vbar(Istr-1,j,knew)=vbar(Istr-1,j,knew)*vmask(Istr-1,j)
#   endif
        enddo
#  else
!                                           Western edge gradient BC
!                                           ======= ==== ======== ==
       do j=JstrV,Jend
          vbar(Istr-1,j,knew)=vbar(Istr,j,knew)
#   ifdef MASKING
     &                        *vmask(Istr-1,j)
#   endif
       enddo
#  endif
# else
#  ifdef NS_COM_PERIODIC
#   define J_RANGE JstrV,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif
                                   ! Closed BC: free-slip (gamma2=+1)
        do j=J_RANGE               ! ====== ===   no-slip (gamma2=-1)
          vbar(Istr-1,j,knew)=gamma2*vbar(Istr,j,knew)
#  ifdef MASKING
     &                                *vmask(Istr-1,j)
#  endif
        enddo
#  undef J_RANGE
# endif              /* OBC_COM_WEST */
      endif       !<-- WESTERN_EDGE
!
!====================================================================
!                            EASTERN BC
!====================================================================
      if (EASTERN_EDGE) then
# ifdef OBC_COM_EAST
#  if defined OBC_COM_M2SPECIFIED || defined OBC_COM_M2SPECIFIED_EAST
!                                           Eastern edge Specified BC
!                                           ======= ==== ========= ==
        do j=JstrV,Jend
#   ifdef M2_FRC_BRY
          vbar(Iend+1,j,knew)=vbarbry_east(j)
#   else
          vbar(Iend+1,j,knew)=vbclm(Iend+1,j)
#   endif
#   ifdef MASKING
     &                        *vmask(Iend+1,j)
#   endif
        enddo
#  elif defined OBC_COM_M2ORLANSKI
!                                           Eastern edge radiation BC
!                                           ======= ==== ========= ==
        do j=JstrV-1,Jend
          grad(Iend  ,j)=vbar(Iend  ,j+1,kstp)-vbar(Iend  ,j,kstp)
          grad(Iend+1,j)=vbar(Iend+1,j+1,kstp)-vbar(Iend+1,j,kstp)
        enddo
        do j=JstrV,Jend
          dft=vbar(Iend,j,kstp)-vbar(Iend  ,j,knew)
          dfx=vbar(Iend,j,knew)-vbar(Iend-1,j,knew)

          if (dfx*dft .lt. 0.) then
            dft=0.                                 ! <-- INFLOW
#   if defined M2_FRC_BRY || defined M2CLIMATOLOGY
            tau=tau_in
          else
            tau=tau_out
#   endif
          endif

          if (dft*(grad(Iend,j-1)+grad(Iend,j)) .gt. 0.) then
            dfy=grad(Iend,j-1)
          else
            dfy=grad(Iend,j  )
          endif

#   ifdef OBC_COM_RAD_NORMAL
          dfy=0.
#   endif
          cff=max(dfx*dfx+dfy*dfy, eps)
          cx=dft*dfx
#   ifdef OBC_COM_RAD_NPO
          cy=0.
#   else
          cy=min(cff,max(dft*dfy,-cff))
#   endif

          vbar(Iend+1,j,knew)=( cff*vbar(Iend+1,j,kstp)
     &                          +cx*vbar(Iend  ,j,knew)
     &                     -max(cy,0.)*grad(Iend+1,j-1)
     &                     -min(cy,0.)*grad(Iend+1,j  )
     &                                       )/(cff+cx)
#   if defined M2_FRC_BRY || defined M2CLIMATOLOGY
          vbar(Iend+1,j,knew)=(1.-tau)*vbar(Iend+1,j,knew)
#    ifdef M2_FRC_BRY
     &                                +tau*vbarbry_east(j)
#    else
     &                                +tau*vbclm(Iend+1,j)
#    endif
#   endif
#   ifdef MASKING
          vbar(Iend+1,j,knew)=vbar(Iend+1,j,knew)*vmask(Iend+1,j)
#   endif
        enddo
#  elif defined OBC_COM_M2CHARACT || defined OBC_COM_M2CHARACT_EAST
!                                           Eastern edge Chapman BC
!                                           ======= ==== ======= ==
        do j=JstrV,Jend
          cff=sqrt(0.5*g*(h(Iend+1,j-1)+zeta(Iend+1,j-1,kstp)+
     &                    h(Iend+1,j  )+zeta(Iend+1,j  ,kstp)))
          cx=dtfast*cff*0.5*(pm(Iend+1,j-1)+pm(Iend+1,j))
          vbar(Iend+1,j,knew)=(vbar(Iend+1,j,kstp)
     &                     +cx*vbar(Iend  ,j,knew) )/(1.+cx)
#   ifdef MASKING
     &                        *vmask(Iend+1,j)
#   endif
        enddo
#  elif defined OBC_COM_M2CHARACT0
!                                      Eastern edge Characteristic BC
!                                      ======= ==== ============== ==
        do j=JstrV,Jend
          if (ubar(Iend-1,j,knew).gt.0.) then
            vbar(Iend+1,j,knew)=vbar(Iend,j,knew)
          else
#   ifdef M2_FRC_BRY
            vbar(Iend+1,j,knew)=vbarbry_east(j)
#   else
            vbar(Iend+1,j,knew)=vbclm(Iend+1,j)
#   endif
          endif
#   ifdef MASKING
          vbar(Iend+1,j,knew)=vbar(Iend+1,j,knew)*vmask(Iend+1,j)
#   endif
        enddo
#  else
!                                            Eastern edge gradient BC
!                                            ======= ==== ======== ==
        do j=JstrV,Jend
          vbar(Iend+1,j,knew)=vbar(Iend,j,knew)
#   ifdef MASKING
     &                        *vmask(Iend+1,j)
#   endif
        enddo
#  endif
# else
#  ifdef NS_COM_PERIODIC
#   define J_RANGE JstrV,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif
                                        ! Wall: free-slip (gamma2=+1)
        do j=J_RANGE                    ! =====   no-slip (gamma2=-1)
          vbar(Iend+1,j,knew)=gamma2*vbar(Iend,j,knew)
#  ifdef MASKING
     &                                *vmask(Iend+1,j)
#  endif
        enddo
#  undef J_RANGE
# endif            /* OBC_COM_EAST */
      endif      !<-- EASTERN_EDGE
#endif             /* !EW_COM_PERIODIC */

!====================================================================
                           ! Corners between adjacent open boundaries
                           ! ======= ======= ======== ==== ==========

#if defined OBC_COM_SOUTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        vbar(Istr-1,Jstr,knew)=0.5*( vbar(Istr-1,Jstr+1,knew)
     &                              +vbar(Istr  ,Jstr  ,knew))
#  ifdef MASKING
     &                         *vmask(Istr-1,Jstr)
#  endif
      endif
#endif
#if defined OBC_COM_SOUTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        vbar(Iend+1,Jstr,knew)=0.5*( vbar(Iend+1,Jstr+1,knew)
     &                              +vbar(Iend  ,Jstr  ,knew))
#  ifdef MASKING
     &                         *vmask(Iend+1,Jstr)
#  endif
      endif
#endif
#if defined OBC_COM_NORTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        vbar(Istr-1,Jend+1,knew)=0.5*( vbar(Istr-1,Jend,knew)
     &                                +vbar(Istr,Jend+1,knew))
#  ifdef MASKING
     &                         *vmask(Istr-1,Jend+1)
#  endif
      endif
#endif
#if defined OBC_COM_NORTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        vbar(Iend+1,Jend+1,knew)=0.5*( vbar(Iend+1,Jend,knew)
     &                                +vbar(Iend,Jend+1,knew))
#  ifdef MASKING
     &                         *vmask(Iend+1,Jend+1)
#  endif
      endif
#endif

!------------------------------------------------------------
#ifdef NBQ
       if (M2bc_nbq_flag) then
#endif
!------------------------------------------------------------

#ifdef WET_DRY
!=======================================================================
!  Impose wetting and drying conditions.
!=======================================================================
!
# ifndef EW_COM_PERIODIC
      if (WESTERN_EDGE) then
        DO j=JstrV,Jend
          cff1=ABS(ABS(vmask_wet(Istr-1,j))-1.)
          cff2=0.5+SIGN(0.5,vbar(Istr-1,j,knew))*vmask_wet(Istr-1,j)
          vmask_wet(Istr-1,j)=0.5*vmask_wet(Istr-1,j)*cff1
     &                                               +cff2*(1.-cff1)
          vbar(Istr-1,j,knew)=vbar(Istr-1,j,knew)*vmask_wet(Istr-1,j)
        END DO
      END IF
      if (EASTERN_EDGE) then
        DO j=JstrV,Jend
          cff1=ABS(ABS(vmask_wet(Iend+1,j))-1.)
          cff2=0.5+SIGN(0.5,vbar(Iend+1,j,knew))*vmask_wet(Iend+1,j)
          vmask_wet(Iend+1,j)=0.5*vmask_wet(Iend+1,j)*cff1
     &                                               +cff2*(1.-cff1)
          vbar(Iend+1,j,knew)=vbar(Iend+1,j,knew)*vmask_wet(Iend+1,j)
        END DO
      END IF
# endif
# ifndef NS_COM_PERIODIC
      if (SOUTHERN_EDGE) then
        DO i=Istr,Iend
          cff1=ABS(ABS(vmask_wet(i,Jstr))-1.)
          cff2=0.5+SIGN(0.5,vbar(i,Jstr,knew))*vmask_wet(i,Jstr)
          vmask_wet(i,Jstr)=0.5*vmask_wet(i,Jstr)*cff1
     &                                           +cff2*(1.-cff1)
          vbar(i,Jstr,knew)=vbar(i,Jstr,knew)*vmask_wet(i,Jstr)
        END DO
      END IF
      if (NORTHERN_EDGE) then
        DO i=Istr,Iend
          cff1=ABS(ABS(vmask_wet(i,Jend+1))-1.)
          cff2=0.5+SIGN(0.5,vbar(i,Jend+1,knew))*vmask_wet(i,Jend+1)
          vmask_wet(i,Jend+1)=0.5*vmask_wet(i,Jend+1)*cff1
     &                                               +cff2*(1.-cff1)
          vbar(i,Jend+1,knew)=vbar(i,Jend+1,knew)*vmask_wet(i,Jend+1)
        END DO
      END IF
# endif
# if !defined EW_COM_PERIODIC && !defined NS_COM_PERIODIC
      if (SOUTHERN_EDGE .and. WESTERN_EDGE) then
        cff1=ABS(ABS(vmask_wet(Istr-1,Jstr))-1.)
        cff2=0.5+SIGN(0.5,vbar(Istr-1,Jstr,knew))*vmask_wet(Istr-1,Jstr)
        vmask_wet(Istr-1,Jstr)=0.5*vmask_wet(Istr-1,Jstr)*cff1
     &                                                   +cff2*(1.-cff1)
        vbar(Istr-1,Jstr,knew)=vbar(Istr-1,Jstr,knew)
     &                                        *vmask_wet(Istr-1,Jstr)
      END IF
      if (SOUTHERN_EDGE .and. EASTERN_EDGE) then
        cff1=ABS(ABS(vmask_wet(Iend+1,Jstr))-1.)
        cff2=0.5+SIGN(0.5,vbar(Iend+1,Jstr,knew))*vmask_wet(Iend+1,Jstr)
        vmask_wet(Iend+1,Jstr)=0.5*vmask_wet(Iend+1,Jstr)*cff1
     &                                                   +cff2*(1.-cff1)
        vbar(Iend+1,Jstr,knew)=vbar(Iend+1,Jstr,knew)
     &                                        *vmask_wet(Iend+1,Jstr)
      END IF
      if (NORTHERN_EDGE .and. WESTERN_EDGE) then
        cff1=ABS(ABS(vmask_wet(Istr-1,Jend+1))-1.)
        cff2=0.5+SIGN(0.5,vbar(Istr-1,Jend+1,knew))
     &                                       *vmask_wet(Istr-1,Jend+1)
        vmask_wet(Istr-1,Jend+1)=0.5*vmask_wet(Istr-1,Jend+1)*cff1
     &                                              +cff2*(1.-cff1)
        vbar(Istr-1,Jend+1,knew)=vbar(Istr-1,Jend+1,knew)
     &                                        *vmask_wet(Istr-1,Jend+1)
      END IF
      if (NORTHERN_EDGE .and. EASTERN_EDGE) then
        cff1=ABS(ABS(vmask_wet(Iend+1,Jend+1))-1.)
        cff2=0.5+SIGN(0.5,vbar(Iend+1,Jend+1,knew))
     &                                       *vmask_wet(Iend+1,Jend+1)
        vmask_wet(Iend+1,Jend+1)=0.5*vmask_wet(Iend+1,Jend+1)*cff1
     &                                              +cff2*(1.-cff1)
        vbar(Iend+1,Jend+1,knew)=vbar(Iend+1,Jend+1,knew)
     &                                        *vmask_wet(Iend+1,Jend+1)
      END IF
# endif
#endif  /* WET_DRY */

#ifdef NBQ
!=======================================================================
!  Computes normal boundary values of DU_nbq for NBQ OBCs
!=======================================================================
!
# ifndef NS_COM_PERIODIC
      if (SOUTHERN_EDGE) then
#  if defined OBC_COM_SOUTH || defined MRL_WCI
        do i=Istr,Iend
          DV_nbq(i,Jstr)=vbar(i,Jstr,knew)
     &                   *0.5*(h(i,Jstr  )+zeta(i,Jstr  ,knew)
     &                       + h(i,Jstr-1)+zeta(i,Jstr-1,knew))
#   ifdef NBQ_MASS
     &                         *0.5*(rhobar_nbq(i,Jstr  ,knew)
     &                              +rhobar_nbq(i,Jstr-1,knew))
#   endif
        enddo
#  else
        do i=Istr,Iend
          DV_nbq(i,Jstr)=0.  ! closed boundary
        enddo
#  endif
      endif
!
      if (NORTHERN_EDGE) then
#  if defined OBC_COM_NORTH || defined MRL_WCI
        do i=Istr,Iend
          DV_nbq(i,Jend+1)=vbar(i,Jend+1,knew)
     &                   *0.5*(h(i,Jend+1)+zeta(i,Jend+1,knew)
     &                       + h(i,Jend  )+zeta(i,Jend  ,knew))
#   ifdef NBQ_MASS
     &                           *0.5*(rhobar_nbq(i,Jend+1,knew)
     &                                +rhobar_nbq(i,Jend  ,knew))
#   endif
        enddo
#  else
        do i=Istr,Iend
          DV_nbq(i,Jend+1)=0.
        enddo
#  endif
      endif
# endif /* NS_COM_PERIODIC */
!
# ifndef EW_COM_PERIODIC
      if (WESTERN_EDGE) then
#  if defined OBC_COM_WEST || defined MRL_WCI
        do j=JstrV,Jend
          DV_nbq(Istr-1,j)=vbar(Istr-1,j,knew)
     &                 *0.5*(h(Istr-1,j  )+zeta(Istr-1,j,knew)
     &                      +h(Istr-1,j-1)+zeta(Istr-1,j,knew))
#   ifdef NBQ_MASS
     &                         *0.5*(rhobar_nbq(Istr-1,j,knew)
     &                              +rhobar_nbq(Istr-1,j,knew))
#   endif

        enddo
#  else
        do j=JstrV,Jend
          DV_nbq(Istr-1,j)=gamma2*DV_nbq(Istr,j)  ! slip/noslip boundary
        enddo
#  endif
      endif
!
      if (EASTERN_EDGE) then
#  if defined OBC_COM_EAST  || defined MRL_WCI
        do j=JstrV,Jend
          DV_nbq(Iend+1,j)=vbar(Iend+1,j,knew)
     &                 *0.5*(h(Iend+1,j  )+zeta(Iend+1,j,knew)
     &                      +h(Iend+1,j-1)+zeta(Iend+1,j,knew))
#   ifdef NBQ_MASS
     &                         *0.5*(rhobar_nbq(Iend+1,j,knew)
     &                              +rhobar_nbq(Iend+1,j,knew))
#   endif
        enddo
#  else
        do j=JstrV,Jend
          DV_nbq(Iend+1,j)=gamma2*DV_nbq(Iend,j)  ! slip/noslip boundary
        enddo
#  endif
      endif
# endif /* EW_COM_PERIODIC */
!
      endif ! M2bc_nbq_flag
!
#endif /* NBQ */


      return
      end
#ifndef CHILD
# define CHILD
# ifdef AGRIF
#  include "v2dbc.F"
# endif
# undef CHILD
#endif  /* !CHILD */
