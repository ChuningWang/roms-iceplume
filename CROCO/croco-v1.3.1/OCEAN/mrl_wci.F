! $Id$
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
#if defined MRL_WCI
      subroutine mrl_wci (tile)
      implicit none
      integer  tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
# ifdef WAVE_DRY
      call wavedry_tile(istr,iend,jstr,jend)  
# endif      
      call mrl_wci_tile (istr,iend,jstr,jend,
# ifdef SOLVE3D
     &            A3d(1,1,trd),A3d(1,2,trd),
# endif
     &            A2d(1,1,trd) ,A2d(1,2,trd) ,A2d(1,3,trd) ,
     &            A2d(1,4,trd) ,A2d(1,5,trd) ,A2d(1,6,trd) ,
     &            A2d(1,7,trd) ,A2d(1,8,trd) ,A2d(1,9,trd) ,
     &            A2d(1,10,trd),A2d(1,11,trd),A2d(1,12,trd),
     &            A2d(1,13,trd)
# if (defined WAVE_FRICTION && !defined WKB_WWAVE) || defined WAVE_STREAMING
     &           ,A2d(1,14,trd)
# endif
     &                  )
      return
      end

      subroutine mrl_wci_tile (istr,iend,jstr,jend, 
# ifdef SOLVE3D
     &            wrk1,wrk2,
# endif
     &            wh,fr,kw,brk,stk,Dstp,act,kD,inv_d,inv_f,frc,ebrk,erol
# if (defined WAVE_FRICTION && !defined WKB_WWAVE) || defined WAVE_STREAMING
     &           ,zob_my
# endif
     &     )
!
!----------------------------------------------------------------------
!   Evaluating wave-averaged terms and Stokes drift based on MRL04
!     (McWilliams, Restrepo & Lane, 2004, JFM, 511, pp.135-178)
!    
!   inputs: R.M.S. wave height (m);
!           peak wave frequency (rad/s);
!           mean wave direction: wdrx & wdre (non dimensional);
!           breaking dissipation   (epsilon_b / rho, m3/s3)
!           roller dissipation     (epsilon_r / rho, m3/s3) and
!           frictional dissipation (epsilon_d / rho, m3/s3)
!
!   Note that if wkb_wwave.F is used, all the dissipation terms are
!   divided by wave frequency, sigma (2pi/T).
!
!   REFERENCE
!   ---------
!   Uchiyama, Y., McWilliams, J.C. and Shchepetkin, A.F. (2010): 
!      Wave-current interaction in an oceanic circulation model with a 
!      vortex force formalism: Application to the surf zone.
!      Ocean Modelling Vol. 34:1-2, pp.16-35.
!
!----------------------------------------------------------------------
!
      implicit none
# include "param.h"
# include "forces.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# ifdef SOLVE3D
#  include "ocean3d.h"
#  ifdef WAVE_STREAMING
#   include "mixing.h"
#  endif
# endif
# ifdef ANA_BRY
#  include "boundary.h"
# endif
# ifdef WKB_WWAVE
#  include "wkb_wwave.h"
# endif
# ifdef BBL
#  include "bbl.h"
# endif
      integer istr,iend,jstr,jend,i,j,imin,jmin 
      real    cff,cff1,cff2,cff3,cff4,inv_g,khd,kh,
     &   eps,wave_ramp,wramp2,kbrk,fb,fb0,fb1,fb2,inv_fbs,
     &   bconst,tauc,tauw,ka_f00,hz0,z_tide,
     &   dtinv,inv_zb,a_brk,a_kv,c1o3,c4o3,inv_k,dd,ust_ker,intfb,
     &   fbsrf,fn1,fn2,khmax,bz1,bz2,beta,abot,delta,ks,
     &   a_frc,kfrc,Hrms,urms,fw,
     &             wh(PRIVATE_2D_SCRATCH_ARRAY),
     &             fr(PRIVATE_2D_SCRATCH_ARRAY),
     &             kw(PRIVATE_2D_SCRATCH_ARRAY),
     &            brk(PRIVATE_2D_SCRATCH_ARRAY),
     &            stk(PRIVATE_2D_SCRATCH_ARRAY),
     &           Dstp(PRIVATE_2D_SCRATCH_ARRAY),
     &            act(PRIVATE_2D_SCRATCH_ARRAY),
     &             kD(PRIVATE_2D_SCRATCH_ARRAY),
     &          inv_d(PRIVATE_2D_SCRATCH_ARRAY),
     &          inv_f(PRIVATE_2D_SCRATCH_ARRAY),
     &            frc(PRIVATE_2D_SCRATCH_ARRAY),
     &           ebrk(PRIVATE_2D_SCRATCH_ARRAY),
     &           erol(PRIVATE_2D_SCRATCH_ARRAY)
# if (defined WAVE_FRICTION && !defined WKB_WWAVE) || defined WAVE_STREAMING
     &        ,zob_my(PRIVATE_2D_SCRATCH_ARRAY)
# endif
# if defined WAVE_BREAK_BJ78
     &        ,Q0,Qb,HmaxBJ     
# endif      
# ifdef SOLVE3D
      integer  k,kk 
      real     wrk1(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &         wrk2(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &         wrk4(0:N),wrk5(0:N),wrk6(0:N), kvsurf
# endif
      parameter ( eps    = 1.e-10, 
     &            bconst = 0.03,              ! breaking contribution to KPP
     &            a_kv   = 1.2,               ! breaking scale for eddy visc.
     &            a_brk  = 0.2,               ! breaking scale for body force
     &            a_frc  = 3.,                ! friction scale for body force
     &            khmax  = 20.,               ! deep-water limit for k*H
     &            c1o3   = 0.3333333333333333,! 1/3
     &            c4o3   = 1.3333333333333333 ! 4/3
     &           )
!
# ifndef WKB_WWAVE
      real wave_gam,wave_btg,wave_roller      ! parameters for offline
#  ifdef WAVE_BREAK_BJ78
      parameter (wave_gam=0.73, wave_btg=1.0)
#  else      
      parameter (wave_gam=0.4, wave_btg=0.8)  ! wave forcing
#  endif      
#  ifdef WAVE_ROLLER
      parameter (wave_roller=1.0)
#  else
      parameter (wave_roller=0.0)
#  endif
# endif /* !WKB_WWAVE */
!
#include "compute_auxiliary_bounds.h"

# if (defined WAVE_FRICTION && !defined WKB_WWAVE) || defined WAVE_STREAMING
      if (Zobt.gt.0.) then                    ! default Zob for wave
        do j=jstr-1,jend+1                    ! dissipation or streaming
          do i=istr-1,iend+1      
#  ifdef BBL
            Zob_my(i,j)=Zbapp(i,j)
#  else
            Zob_my(i,j)=Zob(i,j)
#  endif
          enddo
        enddo                         
      else
        Zob_my=1.e-3
      endif
#  define Zob Zob_my
# endif

# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
# ifdef WET_DRY
#  define SWITCH_WET *
# else
#  define SWITCH_WET !
# endif
# ifdef WAVE_DRY
#  define SWITCH_WAVEDRY *
# else
#  define SWITCH_WAVEDRY !
# endif      

! explicit wavenumber estimator

# undef  KH_SOULSBY
# define KH_HUNT

! vertical distribution function F_Kv(z) for eddy viscosity

# undef  FKV_FUNC1
# define FKV_FUNC2
# undef  FKV_FUNC3

! vertical distribution function F_B(z) for breaking accerelation

# undef  FB_FUNC0
# undef  FB_FUNC1
# undef  FB_FUNC2
# define FB_FUNC3
# undef  FB_WSCALE

! initial ramping coefficients

# if defined WAVE_RAMP
#  if defined RIP_TOPO_2D
      wave_ramp = tanh(192.*dt*sec2day*float(iic-ntstart))
#  else
      wave_ramp = tanh(12.*dt*sec2day*float(iic-ntstart))
#  endif
# else
      wave_ramp = 1.
# endif
      if (iic.eq.0) wave_ramp=1.
      wramp2= wave_ramp**2
!
! Evaluate 2DH wave-current interaction variables.
! ================================================
!
! 2DH Stokes velocities, breaking, roller and bottom-friction 
! dissipation terms, defined at horizontal rho-points.
!
      inv_g = 1./g
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
# ifdef SOLVE3D
          Dstp(i,j)=h(i,j)+z_w(i,j,N)
# else
          Dstp(i,j)=h(i,j)+zeta(i,j,knew)
# endif          
          inv_d(i,j)=1./Dstp(i,j)

# ifdef WKB_WWAVE
!----------------------------------------------------------------------
!  Wave parameters from WKB model
!----------------------------------------------------------------------
!
          fr(i,j)=frq(i,j,wnew)
          wh(i,j)=wave_ramp*hrm(i,j,wnew)
          kw(i,j)=max(wvn(i,j,wnew),eps)
#  ifndef WAVE_ROLLER
          cff=wramp2*wsb(i,j,wnew)
          ebrk(i,j)=cff*fr(i,j)                  ! epsilon_b (m3/s3)
#  else
          ebrk(i,j)=wramp2*wsb(i,j,wnew)*fr(i,j) ! epsilon_b (m3/s3)
          erol(i,j)=wramp2*wsr(i,j,wnew)*fr(i,j) ! epsilon_r (m3/s3)
          cff =wramp2*( wsr(i,j,wnew)+           ! for breaking term
     &           (1.-wkb_roller)*wsb(i,j,wnew) ) ! primary + roller
#  endif
          brk(i,j)=cff*kw(i,j)                   ! 2DH breaking term
          frc(i,j)=wramp2*wfc(i,j,wnew)*kw(i,j)  ! 2DH friction term
#  ifdef MRL_CEW
          kD(i,j)=kw(i,j)*Dstp(i,j)
#  else
!          kD(i,j)=kw(i,j)*(h(i,j)+wkb_tide)
          kD(i,j)=kw(i,j)*Dstp(i,j)
#  endif
#  ifndef WAVE_ROLLER
          act(i,j)=wramp2*max(wac(i,j,wnew),0.)  ! wave action density
#  else
          act(i,j)=wramp2*(max(wac(i,j,wnew),0.)
     &                    +max(war(i,j,wnew),0.))
#  endif

# else    /* ifndef WKB_WWAVE */
!
!----------------------------------------------------------------------
!  Wave parameters from input file (read in get_wwave)
!----------------------------------------------------------------------
!
          wh(i,j)=wave_ramp*whrm(i,j)
          fr(i,j)=wfrq(i,j)
          inv_f(i,j)=1./max(fr(i,j),eps)
          khd=Dstp(i,j)*(fr(i,j)**2)*inv_g
#  ifdef BBL
! --- copy into arrays for BBL Soulsby parametrization
          Pwave(i,j) = 2.0*pi*inv_f(i,j)
          Awave(i,j) = wh(i,j)/2.0
          Dwave(i,j) = atan2(wdre(i,j),wdrx(i,j))
#  endif  

#  ifdef OW_COUPLING_FULL
          kw(i,j) =2.0*pi/max(wlm(i,j),0.1)
          !kw(i,j) =2.0*pi/max(wlm(i,j),eps)
          kD(i,j) =kw(i,j)*Dstp(i,j)
#  else          

#   ifdef KH_SOULSBY
          if (khd.ge.1.) then ! explicit wavenumber estimator
           kh=khd
          else
           kh=sqrt(khd)
          endif
          do k=1,3
            cff = tanh(kh)
            kh = kh-(kh*cff-khd)/max(cff+kh*(1.0-cff**2),eps)
          enddo
#   elif defined KH_HUNT
          kh = sqrt( khd*khd + khd/(1.0 + khd*(0.6666666666
     &               +khd*(0.3555555555 + khd*(0.1608465608
     &               +khd*(0.0632098765 + khd*(0.0217540484
     &                            +khd*0.0065407983)))))) )
#   endif
          kD(i,j) =kh
          kw(i,j) =kh*inv_d(i,j)
#  endif          
!
! --- Wave breaking ---
!   use wave_btg, wave_gam local parameters (as in WKB model)
!
#  ifdef WAVE_OFFLINE_BREAKING
          ebrk(i,j)=wramp2*wepb(i,j)                ! ep_b/rho (m3/s3)
#  else
#   ifdef OW_COUPLING_FULL
          ebrk(i,j)=wramp2*foc(i,j)/rho0 
#   else          
#    ifdef WAVE_BREAK_TG86
           cff=3./16.*sqrt(pi)*g*(wave_btg**3)/(wave_gam**4)/(2.*pi)
           ebrk(i,j)=cff/(Dstp(i,j)**5)*(wh(i,j)**7) ! ep_b/rho/sigma (m3/s2)
#    elif defined WAVE_BREAK_TG86A
           cff=3./16.*sqrt(pi)*g*(wave_btg**3)/(wave_gam**2)/(2.*pi)
           cff1=wh(i,j)/(wave_gam*Dstp(i,j))
           ebrk(i,j)=cff/(Dstp(i,j)**3)*(wh(i,j)**5)*
     &                           (1.-(1.+cff1**2)**(-2.5))
#    elif defined WAVE_BREAK_CT93
           cff=3./16.*sqrt(pi)*g*(wave_btg**3)/(2.*pi)
           cff1=wh(i,j)/(wave_gam*Dstp(i,j))
           ebrk(i,j)=cff/Dstp(i,j)*(wh(i,j)**3)
     &              *( 1.+tanh(8.0*(cff1-1.)) )*
     &               ( 1.-(1.0+cff1**2)**(-2.5) )
#    elif defined WAVE_BREAK_R93
           cff=0.25*g*(wave_btg**3)/(2.*pi)
           cff1=wh(i,j)/(wave_gam*Dstp(i,j))
           ebrk(i,j)=cff/Dstp(i,j)*wh(i,j)**3*cff1**4
     &                          *(1-exp (-cff1**5 ))
#    elif defined WAVE_BREAK_BJ78
           HmaxBJ=wave_gam/max(kw(i,j),eps)*tanh(kw(i,j)*Dstp(i,j))
           if (HmaxBJ.gt.eps) then
             cff1=wh(i,j)/HmaxBJ !=beta_BJ
           else
             cff1=0.0
           endif
           if (cff1.le.0.5) then
             Q0 = 0.0
           elseif(cff1.le.1.0) then
             Q0 = (2.0*cff1-1.0)**2
           endif
           if(cff1.le.0.2) then
             Qb = 0.0
           elseif(cff1.lt.1.0) then
             cff = exp((Q0-1.0)/cff1**2)
             Qb = Q0-cff1**2*(Q0-cff)/(cff1**2-cff)
             do k=1,3
               Qb = exp((Qb-1.0)/cff1**2)
             enddo
           else
             Qb = 1.0-eps
           endif
           if( (cff1**2.gt.eps) .and. (abs(cff1**2-Qb).gt.eps) ) then
             if(cff1**2.lt.1.0) then
               ebrk(i,j) = 0.25*wave_btg*g*Qb*HmaxBJ**2/(2.0*pi)
             else
               ebrk(i,j) = 0.25*wave_btg*g*cff1**4*HmaxBJ**2/(2.0*pi)
             endif  
           else
             ebrk(i,j) = 0.0
           endif  
#    endif  
           ebrk(i,j)=wramp2*ebrk(i,j)*fr(i,j)        ! ep_b/rho (m3/s3)
#   endif       
#  endif
!
#  ifdef WAVE_OFFLINE_ROLLER
          erol(i,j=wramp2*0.5762*Dstp(i,j)*         ! ep_r (m3/s3)
     &             wh(i,j)*wepr(i,j)*               ! 0.5762=0.06*sinb*g^2
     &             kw(i,j)*inv_f(i,j)               ! Nairn et al (1991)
          cff=(1.-wave_roller)*ebrk(i,j)+erol(i,j)
#  else
          cff=ebrk(i,j)
#  endif
          brk(i,j)=cff*kw(i,j)*inv_f(i,j)           ! 2DH breaking term
!
! --- Wave friction ---
!
#  ifdef WAVE_OFFLINE_FRICTION
          frc(i,j)=wramp2*wepd(i,j)*kw(i,j)
     &                             *inv_f(i,j)      ! 2DH friction term
#  elif defined WAVE_FRICTION
#   if defined OW_COUPLING_FULL
          urms=ubr(i,j)
          fw=1.39*(fr(i,j)*Zob(i,j)/max(urms,eps))**0.52
#   else          
          abot=wh(i,j)/max(2.*sinh(min(kD(i,j),khmax)),eps)
          urms=fr(i,j)*abot
          fw=1.39*(Zob(i,j)/max(abot,eps))**0.52    ! Soulsby (1995)
#   endif          
          cff=0.5/sqrt(pi)*fw*urms**3               ! ep_d/rho
          frc(i,j)=wramp2*cff*kw(i,j)*inv_f(i,j)    ! 2DH friction term
#  else
          frc(i,j)=0.
#  endif
!
! --- Wave action ---
!
          act(i,j)=inv_f(i,j)*(0.125*g*(wh(i,j)**2) ! wave action density
#  ifdef WAVE_OFFLINE_ROLLER
     &            +wramp2*0.03*g*Dstp(i,j)*wh(i,j)*wepr(i,j)
#  endif
     &                         )

# endif   /* WKB_WWAVE */
!
!----------------------------------------------------------------------
!  Finalize wave parameter computation
!----------------------------------------------------------------------
!
# ifdef MASKING
           kw(i,j) = kw(i,j)*rmask(i,j)
           kD(i,j) = kD(i,j)*rmask(i,j)
          act(i,j) =act(i,j)*rmask(i,j)
          brk(i,j) =brk(i,j)*rmask(i,j)
          frc(i,j) =frc(i,j)*rmask(i,j)
# endif /* MASKING */
# ifdef WET_DRY
           kw(i,j) = kw(i,j)*rmask_wet(i,j)
           kD(i,j) = kD(i,j)*rmask_wet(i,j)
          act(i,j) =act(i,j)*rmask_wet(i,j)
          brk(i,j) =brk(i,j)*rmask_wet(i,j)
          frc(i,j) =frc(i,j)*rmask_wet(i,j)
# endif /* WET_DRY */
# ifdef WAVE_DRY
           kw(i,j) = kw(i,j)*rmask_wavewet(i,j)
           kD(i,j) = kD(i,j)*rmask_wavewet(i,j)
          act(i,j) =act(i,j)*rmask_wavewet(i,j)
          brk(i,j) =brk(i,j)*rmask_wavewet(i,j)
          frc(i,j) =frc(i,j)*rmask_wavewet(i,j)
# endif          
        enddo
      enddo              ! <-- discard inv_d

!
! sup : quasi-static sea-level response, set-up. no interaction.
! ==============================================================
!
      do j=jstrR,jendR
        do i=istrR,iendR
          sup(i,j)=-0.125*(wh(i,j)**2)*kw(i,j)
     &              /max(sinh(2.*min(kD(i,j),khmax)),eps)
     &                                 SWITCH rmask(i,j)
     &                                 SWITCH_WET rmask_wet(i,j)
     &                                 SWITCH_WAVEDRY rmask_wavewet(i,j)
        enddo
      enddo
!
! Dissipation terms : copy into shared arrays
! ===========================================
!
      do j=jstrR,jendR
        do i=istrR,iendR
          wepb(i,j)=ebrk(i,j)         ! epsilon_b (m3/s3)
# ifdef WAVE_ROLLER
          wepr(i,j)=erol(i,j)         ! epsilon_r (m3/s3)
# endif
# ifdef WKB_WWAVE
          cff = wramp2*wfc(i,j,wnew)
# else
          cff = wramp2*frc(i,j)/kw(i,j)
# endif
          wepd(i,j)=cff*fr(i,j)       ! epsilon_d (m3/s3) 
# ifdef MASKING
          wepb(i,j)=wepb(i,j)*rmask(i,j)
          wepd(i,j)=wepd(i,j)*rmask(i,j)
#  ifdef WAVE_ROLLER
          wepr(i,j)=wepr(i,j)*rmask(i,j)
#  endif
# endif
# ifdef WET_DRY
          wepb(i,j)=wepb(i,j)*rmask_wet(i,j)
          wepd(i,j)=wepd(i,j)*rmask_wet(i,j)
#  ifdef WAVE_ROLLER
          wepr(i,j)=wepr(i,j)*rmask_wet(i,j)
#  endif
# endif
# ifdef WAVE_DRY
          wepb(i,j)=wepb(i,j)*rmask_wavewet(i,j)
          wepd(i,j)=wepd(i,j)*rmask_wavewet(i,j)
#  ifdef WAVE_ROLLER
          wepr(i,j)=wepr(i,j)*rmask_wavewet(i,j)
#  endif
# endif
        enddo
      enddo
# ifndef WKB_WWAVE
#  undef inv_f
# endif
!
! Wavenumber : copy into shared arrays for output
! ===============================================
!
      do j=jstrR,jendR
        do i=istrR,iendR
          wwkx(i,j)=kw(i,j)*wdrx(i,j)
          wwke(i,j)=kw(i,j)*wdre(i,j)
        enddo
      enddo
!
! 2DH depth-averaged breaking acceleration
! =========================================
!
      do j=jstrR,jendR
        do i=istr,iendR
          brk2dx(i,j) =0.5*(brk(i-1,j)*wdrx(i-1,j)+brk(i,j)*wdrx(i,j))
          frc2dx(i,j) =0.5*(frc(i-1,j)*wdrx(i-1,j)+frc(i,j)*wdrx(i,j))
# ifdef MASKING
          brk2dx(i,j) =brk2dx(i,j)*umask(i,j)
          frc2dx(i,j) =frc2dx(i,j)*umask(i,j)
# endif
# ifdef WET_DRY
          brk2dx(i,j) =brk2dx(i,j)*umask_wet(i,j)
          frc2dx(i,j) =frc2dx(i,j)*umask_wet(i,j)
# endif
# ifdef WAVE_DRY
          brk2dx(i,j) =brk2dx(i,j)*umask_wavewet(i,j)
          frc2dx(i,j) =frc2dx(i,j)*umask_wavewet(i,j)
# endif          
        enddo
      enddo
      do j=jstr,jendR
        do i=istrR,iendR
          brk2de(i,j) =0.5*(brk(i,j-1)*wdre(i,j-1)+brk(i,j)*wdre(i,j))
          frc2de(i,j) =0.5*(frc(i,j-1)*wdre(i,j-1)+frc(i,j)*wdre(i,j))
# ifdef MASKING
          brk2de(i,j) =brk2de(i,j)*vmask(i,j)
          frc2de(i,j) =frc2de(i,j)*vmask(i,j)
# endif
# ifdef WET_DRY
          brk2de(i,j) =brk2de(i,j)*vmask_wet(i,j)
          frc2de(i,j) =frc2de(i,j)*vmask_wet(i,j)
# endif
# ifdef WAVE_DRY
          brk2de(i,j) =brk2de(i,j)*vmask_wavewet(i,j)
          frc2de(i,j) =frc2de(i,j)*vmask_wavewet(i,j)
# endif
        enddo
      enddo

# ifdef SOLVE3D
!
! Conservative 3D wave-current interaction variables.
! ===================================================
!
#  define ust_r wrk1
#  define inv_ex inv_f
!
! 3D depth-dependent Stokes drift velocities
! ==========================================
!
#  ifdef STOKES_DRIFT
!
! Cell-averaged 3D Stokes drift velocity (ust_r) is estimated 
! assuming in surfzone that roller Stokes drift has the same 
! z-dependency as the primary waves:  
!
!                               cosh[2k(z+h)]
!          ust = 1/8 Hrms^2 f k -----------     (Ursell)
!                                sinh[kD]^2
!
!                           cosh[2k(z+h)]
!              = 2 act k^2  ----------- 
!                            sinh[2kD]
!
! with act=0.5gA^2/frq the wave action density (frq^2=gk*tanh(kD))
!
! Here, sinh and cosh functions are expanded to exponential function 
! so as not to be singular when kD gets very large (offshore). 
! In addition, analytical volume-averaging operation at z_r levels
! (between surrounding z_w levels: z_dw and z_up) gives:
!
!          ust = act*k/(1-exp(-4kD))/dz *
!                ( exp(2k(z_up+h-D)) - exp(-2k(z_up+h+D))
!                 -exp(2k(z_dw+h-D)) + exp(-2k(z_dw+h+D)) )
!
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
#   ifdef OW_COUPLING_FULL
          inv_ex(i,j)=1.0/max(2.0*kw(i,j)*(1.0+exp(-4.0*kD(i,j))),eps)
          ust_ker = ust_ext(i,j)*inv_ex(i,j)
#   else  
          inv_ex(i,j)=1./max(1.-exp(-4.*kD(i,j)),eps)        
          ust_ker=act(i,j)*kw(i,j)*inv_ex(i,j)
#   endif
          cff1=0.
          do k=1,N,+1              ! <-- irreversible
            cff2=exp( 2.*kw(i,j)*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &          -exp(-2.*kw(i,j)*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))
            ust_r(i,j,k)=ust_ker*(cff2-cff1)/Hz(i,j,k)
            cff1=cff2
          enddo
        enddo
      enddo
!
! Horizontal interpolation of ust_r to u- & v- points.
! (mask not applied).
!
      do k=1,N
        do j=jstrR,jendR
          do i=istr,iendR
            ust(i,j,k)=0.5*( ust_r(i-1,j,k)*wdrx(i-1,j)
     &                      +ust_r(i  ,j,k)*wdrx(i  ,j) )
          enddo
        enddo
        do j=jstr,jendR
          do i=istrR,iendR
            vst(i,j,k)=0.5*( ust_r(i,j-1,k)*wdre(i,j-1)
     &                      +ust_r(i,j  ,k)*wdre(i,j  ) )
          enddo
        enddo
      enddo
!
!  2DH depth-averaged Stokes drift 
!
      do j=JstrR,JendR
        do i=Istr,IendR
          cff1=0.5*(Hz(i-1,j,1)+Hz(i,j,1))
          ust2d(i,j)=cff1*ust(i,j,1)
          DO k=2,N
            cff1=0.5*(Hz(i-1,j,k)+Hz(i,j,k))
            ust2d(i,j)=ust2d(i,j)+cff1*ust(i,j,k)
          enddo
          cff2=2./(Dstp(i-1,j)+Dstp(i,j))
          ust2d(i,j)=ust2d(i,j)*cff2
#   ifdef MASKING
          ust2d(i,j)=ust2d(i,j)*umask(i,j)
#   endif
        enddo
      enddo
      do j=Jstr,JendR
        do i=IstrR,IendR
          cff1=0.5*(Hz(i,j-1,1)+Hz(i,j,1))
          vst2d(i,j)=cff1*vst(i,j,1)
          do k=2,N
            cff1=0.5*(Hz(i,j-1,k)+Hz(i,j,k))
            vst2d(i,j)=vst2d(i,j)+cff1*vst(i,j,k)
          enddo
          cff2=2./(Dstp(i,j-1)+Dstp(i,j))
          vst2d(i,j)=vst2d(i,j)*cff2
#   ifdef MASKING
          vst2d(i,j)=vst2d(i,j)*vmask(i,j)
#   endif
        enddo
      enddo
#  endif /* STOKES_DRIFT */
!
! calP & Kapsrf: surface pressure & Bernoulli head at rho-point
! =============================================================
!
#  define kv wrk4
#  define dkvdz wrk5
#  define d2kv wrk6
      do j=jstr,jend        ! kv (k dot v) at rho-point
        do i=istr,iend      ! and its 1st & 2nd derivertives
          do k=1,N          ! at rho-point
            kv(k) =0.5*kw(i,j)*(
     &         wdrx(i,j)* ( u(i,j,k,nstp)+u(i+1,j,k,nstp) )
     &       + wdre(i,j)* ( v(i,j,k,nstp)+v(i,j+1,k,nstp) ) )
          enddo
          kvsurf =1.5*kv(N)-0.5*kv(N-1) ! extrapolate to surface
          do k=1,N-1
            dkvdz(k) =2.0*(kv(k+1)-kv(k))/(Hz(i,j,k+1)+Hz(i,j,k))
          enddo
          dkvdz(0) = dkvdz(1)    !2.*dkvdz(1)-dkvdz(2) ! severe!
          dkvdz(N) = dkvdz(N-1)  !2.*dkvdz(N-1)-dkvdz(N-2)
          do k=1,N
            d2kv(k) =dkvdz(k)-dkvdz(k-1) ! d^2kv/dz^2 x Hz
          enddo
          cff3 = 0.0
          do k=1,N
            dd   = z_r(i,j,k)-z_w(i,j,N)
            cff3 = cff3 + d2kv(k)*(
     &             exp( 2.*kw(i,j)*(dd-Dstp(i,j)))
     &           + exp(-2.*kw(i,j)*(dd+Dstp(i,j))) )
          enddo
#  ifndef STOKES_DRIFT
          inv_ex(i,j)=1./max(1.-exp(-4.*kD(i,j)),eps)
#  endif
          cff1 =-2.*exp(-2.*kD(i,j))*inv_ex(i,j)*dkvdz(N)
          cff2 = dkvdz(0)/max(tanh(2.*kD(i,j)),eps)
          cff3 = cff3*inv_ex(i,j)
          cff4 =-2.*kw(i,j)*kvsurf
          calP(i,j) = act(i,j)*inv_g*tanh(kD(i,j))
     &                    *( cff1+cff2+cff3+cff4 )
     &                            SWITCH rmask(i,j)
     &                            SWITCH_WET rmask_wet(i,j)
     &                            SWITCH_WAVEDRY rmask_wavewet(i,j)
          cff=0.
          do k=1,N
            cff = cff +  d2kv(k)*
     &         ( exp( 2.*kw(i,j)*(z_w(i,j,N)-z_r(i,j,k)-Dstp(i,j)))
     &          -exp(-2.*kw(i,j)*(z_w(i,j,N)-z_r(i,j,k)+Dstp(i,j))) )
          enddo
#  define KAPSRF_SAFE
#  ifdef KAPSRF_SAFE
          Kapsrf(i,j) = ( 0.5*Kapsrf(i,j)  ! time-filter to suppress noises
#   ifdef OW_COUPLING_FULL
     &                   +0.5*bhd(i,j) )
#   else          
     &                   +0.5*cff*act(i,j)*inv_ex(i,j) )
#   endif          
#  else
#   ifdef OW_COUPLING_FULL
          Kapsrf(i,j) = bhd(i,j)
#   else          
          Kapsrf(i,j) = cff*act(i,j)*inv_ex(i,j)
#   endif          
#  endif
     &                            SWITCH rmask(i,j)
     &                            SWITCH_WET rmask_wet(i,j)
     &                            SWITCH_WAVEDRY rmask_wavewet(i,j)
        enddo
      enddo      ! <-- discard kv, dkvdz, d2kv, & inv_ex
#  undef kv
#  undef dkvdz
#  undef d2kv
#  undef inv_ex
!
! kvf : vertical vortex force term (K term) at rho-point
! ======================================================
!
      if (WESTERN_EDGE) then
        imin=IstrU-1
      else
        imin=IstrU
      endif
      if (SOUTHERN_EDGE) then
        jmin=JstrV-1
      else
        jmin=JstrV
      endif

#  ifdef STOKES_DRIFT

#   define kvr wrk2
#   define kvfx stk
#   define kvfy inv_f

      do k=1,N-1
        do j=JstrV-1,jend+1             ! du/dz & dv/dz at horz u- &
          do i=IstrU-1,iend+1           ! v- & vert w-point
            kvfx(i,j) = 2.0*(u(i,j,k+1,nstp)-u(i,j,k,nstp))
     &         /(z_r(i-1,j,k+1)+z_r(i,j,k+1)-z_r(i-1,j,k)-z_r(i,j,k))
            kvfy(i,j) = 2.0*(v(i,j,k+1,nstp)-v(i,j,k,nstp))
     &         /(z_r(i,j-1,k+1)+z_r(i,j,k+1)-z_r(i,j-1,k)-z_r(i,j,k))
          enddo
        enddo
        do j=jmin,jend
          do i=imin,iend     ! K term at horz rho- & vert w-point
            kvr(i,j,k) = 0.25*(ust_r(i,j,k)+ust_r(i,j,k+1))*
     &                   ( wdrx(i,j)*(kvfx(i,j)+kvfx(i+1,j))
     &                    +wdre(i,j)*(kvfy(i,j)+kvfy(i,j+1)) )
          enddo
        enddo
      enddo
      do j=jmin,jend                 ! bottom & top B.Cs.
        do i=imin,iend
          kvr(i,j,0)=0.D0
          kvr(i,j,N)=2.*kvr(i,j,N-1)-kvr(i,j,N-2)
        enddo
      enddo
      do k=1,N                       ! move to vert rho-point
        do j=jmin,jend               ! copy into a shared array
          do i=imin,iend
            kvf(i,j,k)=0.5*( kvr(i,j,k)+kvr(i,j,k-1) )
     &                            SWITCH rmask(i,j)
     &                            SWITCH_WET rmask_wet(i,j)
     &                            SWITCH_WAVEDRY rmask_wavewet(i,j)
          enddo
        enddo
      enddo           ! <-- discard ust_r, kvfx, kvfy, kvr

#   undef kvr
#   undef kvfx
#   undef kvfy

#  else
      do k=1,N
        do j=jmin,jend
          do i=imin,iend
            kvf(i,j,k)=0.
          enddo
        enddo
      enddo
#  endif /* STOKES_DRIFT */

#  undef ust_r

!
! Akw: primary wave-induced additional diffusivity for tracers at w-point
! =======================================================================
!
#  ifdef LMD_MIXING
#   define WAVE_NONBRK_VMIX_QIAO
!
#   ifdef WAVE_NONBRK_VMIX_U2010 /* weak and unstable in SHOREFACE case */
!
!  In Uchiyama et al 2010:
!  Akw = 0.5 d/dt[[A * sinh(k(z+d))/sinh(k*d)]^2]
!
      cff=wave_ramp/dt
      do j=jstrR,jendR
        do i=istrR,iendR
          cff1 = 0.25*wh(i,j)/max(1.-exp(-2.*kD(i,j)),eps)
          do k=0,N
            cff2 = cff1*( exp( kw(i,j)*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &                   -exp(-kw(i,j)*(z_w(i,j,k)+h(i,j)+Dstp(i,j))) )
            cff3 = cff2**2
            Akw(i,j,k) = cff*(cff3-E_pre(i,j,k))
     &                            SWITCH rmask(i,j)
     &                            SWITCH_WET rmask_wet(i,j)
     &                            SWITCH_WAVEDRY rmask_wavewet(i,j)
            E_pre(i,j,k) = cff3   ! store this for the next time step.
          enddo
        enddo
      enddo
#   elif defined WAVE_NONBRK_VMIX_QIAO
!
!  Qiao et al, GRL 2004 and OD 2010 (eq. 41b):
!    Akw = alpha * kw * sig * [A * sinh(k(z+d))/sinh(k*d)]^3
!  Wang et al. JGR 2010:
!    alpha=0.1 for KPP applications
!
      do j=jstrR,jendR
        do i=istrR,iendR
          cff1 = 0.5*wh(i,j)/max(1.-exp(-2.*kD(i,j)),eps)
          do k=0,N
            cff2 = cff1*( exp( kw(i,j)*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &                   -exp(-kw(i,j)*(z_w(i,j,k)+h(i,j)+Dstp(i,j))) )
            Akw(i,j,k) = 0.1*kw(i,j)*fr(i,j)*cff2**3
     &                            SWITCH rmask(i,j)
     &                            SWITCH_WET rmask_wet(i,j)
     &                            SWITCH_WAVEDRY rmask_wavewet(i,j)
          enddo
        enddo
      enddo
#   endif
#  endif /* LMD_MIXING */
!
! Non conservative wave-current interaction variables.
! ====================================================
!
! We estimate cell-averaged 3D breaking acceleration (brk_r),
! breaking-enhaced vertical eddy viscosity (Akb), and bottom
! wave streaming-induced acceleration (frc_r).  brk_r, Akb, and
! frc_r are assumed to follow an arbitrary vertical shape function,
! fb (or fkv, fwd for Akb & frc_r) with a surface/bottom scale,
! where kb^-1 = a_brk Hrms (for fb), or kb^-1 = a_kv Hrms (for fkv),
! or kb^-1 = a_frc x delta (for fwd), where Hrms is rms wave height
! and delta is turbulent wave boundary layer thickness.  a_brk,
! a_kv and a_frc are considered O(1) parameters which determine
! vertical penetration of breaking-driven momenta, additional eddy
! kinetic energy by breaking, and streaming-induced momenta near
! the bed.  fb (fkv and fwd) is designed with three (plus one 
! special case for fb) different ways:
!
!  FUNC1: gb~1-tanh(kb(\zeta-z))^4 : Warner et al (2008)
!  FUNC2: gb~1-tanh(kb(\zeta-z))^2 : modified from the above
!  FUNC3: gb~cosh(kb(z+h))         : analogous to primary wave
!  FUNC0: gb~1                     : vertically uniform (Fb only)
!
! where fb(z) = gb(z) / [\int^D gb(z') dz'] = 1, thus vertical
! integral of fb is normalized.  Then volume-averaging operation is 
! applied.  Exactly the same procedure is taken for fkv (not show 
! here), while slightly different consideration is taken for fwd 
! (see below).  In practice, vertical integrals of gb should have 
! the following forms:
!
!   1) FUNC1
!
!                      4      3e^{4kbz}+3e^{2kbz}+2
!   \int gb(z') dz' = ---*-------------------------------
!                     3kb  e^{6kbz}+3e^{4kbz}+3e^{2kbz}+1
!
!                      4  3e^{-2kbz}+3e^{-4kbz}+2e^{-6kbz}
!                   = ---*---------------------------------
!                     3kb 1+3e^{-2kbz}+3e^{-4kbz}+e^{-6kbz}
!
!   2) FUNC2
!
!                      2     1         2  e^{-2kbz}
!   \int gb(z') dz' = --*---------- = --*-----------
!                     kb e^{2kbz}+1   kb 1+e^{-2kbz}
!
!   3) FUNC3
!
!   \int gb(z') dz' = [sinh(kb*z)]/kb
!
!   and then for FUNC3
!
!           kb       1     z2
!   fb = ----------*--*\int  cosh[kb(z+h)] dz
!        sinh(kb*D) dz     z1
!
!                1                                       z2
!      = -----------------*[e^{kb(z+h-D)}-e^{-kb(z+h+D)}]
!        dz*(1-e^{-2kb*D})                               z1
!
! in which hyperbolic functions (sinh, cosh, and tanh) in fb
! are expanded to exponential function so as not to be singular
! when kD is very large.
!
! Additionally, we can introduce gb=1 (vertically uniform function)
! with FB_FUNC0 for experimental purposes.  This option can be
! chosen only for fb since it may be unappropriate to fkv and fwd.
!
! Note that the above expressions are introduced in order to
! avoid the integral being Inf or NaN when kb^-1 -> 0.
!
#  ifdef LMD_MIXING
#   define kb0 wrk4
#   define vs  cff1
#   define kvb cff2
!
! Akb: Eddy viscosity due to depth-induced wave breaking
! ======================================================
!
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
          inv_k =min(max(a_kv*wh(i,j),eps),Dstp(i,j))
          kbrk  =1./inv_k          ! vertical scale for breaking
#   ifdef FKV_FUNC1
          fn1   =3.*(exp(-2.*kbrk*Dstp(i,j))+exp(-4.*kbrk*Dstp(i,j)))
          fn2   =exp(-6.*kbrk*Dstp(i,j))
          fb1   =c4o3*inv_k*(fn1+2.*fn2)/(1.+fn1+fn2)
          intfb =c4o3*inv_k - fb1  ! intfb = \int^D gb(z') dz'
          fb0   =1./intfb          ! 1/[\int^D gb(z') dz']
#   elif defined FKV_FUNC2
          fn1   =exp(-2.*kbrk*Dstp(i,j))
          fb1   =2.*inv_k*fn1/(1.+fn1)
          intfb =inv_k - fb1
          fb0   =1./intfb
#   elif defined FKV_FUNC3
          fb1   =0.                ! fb1 = \int gb dz at z = -h
          fb0   =1./(1.-exp(-2.*kbrk*Dstp(i,j)))
#   endif
#   ifndef WAVE_ROLLER
          vs    =ebrk(i,j)**c1o3   ! representative velocity scale
#   else
          vs    =((1.-wkb_roller)*ebrk(i,j)+erol(i,j))**c1o3
#   endif
          kvb   =bconst*vs*wh(i,j) ! depth-averaged Akb
          do k=1,N,+1              ! <-- irreversible
#   ifdef FKV_FUNC1
            dd  =z_w(i,j,N)-z_w(i,j,k)
            fn1 =3.*(exp(-2.*kbrk*dd)+exp(-4.*kbrk*dd))
            fn2 =exp(-6.*kbrk*dd)
            fb2 =c4o3*inv_k*(fn1+2.*fn2)/(1.+fn1+fn2)
#   elif defined FKV_FUNC2
            fn1 =exp(-2.*kbrk*(z_w(i,j,N)-z_w(i,j,k)))
            fb2 =2.*inv_k*fn1/(1.+fn1)
#   elif defined FKV_FUNC3
            fb2 =exp( kbrk*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &          -exp(-kbrk*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))
#   endif
            fb =fb0*(fb2-fb1)/Hz(i,j,k)
            kb0(k) = kvb*fb*Dstp(i,j)
            fb1=fb2     ! recursive procedure
          enddo
          do k=1,N-1
            Akb(i,j,k)=max(0.5*(kb0(k)+kb0(k+1)),0.)
          enddo
          Akb(i,j,0)=0.
          Akb(i,j,N)=max(1.5*kb0(N)-0.5*kb0(N-1),0.)
#   ifdef MASKING
          do k=0,N
            Akb(i,j,k)=Akb(i,j,k)*rmask(i,j)
          enddo
#   endif
#   ifdef WET_DRY
          do k=0,N
            Akb(i,j,k)=Akb(i,j,k)*rmask_wet(i,j)
          enddo          
#   endif
#   ifdef WAVE_DRY
          do k=0,N
            Akb(i,j,k)=Akb(i,j,k)*rmask_wavewet(i,j)
          enddo          
#   endif
        enddo
      enddo        ! <-- discard kb0
#   undef kb0
#  endif    /* LMD_MIXING */

#  ifndef WAVE_SFC_BREAK
#   define brk_r wrk2
!
! 3D breaking acceleration term defined as a body force.
! =====================================================
!
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
#   ifdef FB_FUNC0
          do k=1,N                 ! vertically uniform case
            brk_r(i,j,k) = brk(i,j)*inv_d(i,j)
          enddo                    ! <-- discard inv_d
#   else
#    ifdef FB_WSCALE
          inv_k =1./max(2.*kw(i,j),eps)  ! Stokes scale
#    else
          inv_k =min(max(a_brk*wh(i,j),eps),Dstp(i,j))
#    endif
          kbrk  =1./inv_k          ! vertical scale for breaking
#    ifdef FB_FUNC1
          fn1   =3.*(exp(-2.*kbrk*Dstp(i,j))+exp(-4.*kbrk*Dstp(i,j)))
          fn2   =exp(-6.*kbrk*Dstp(i,j))
          fb1   =c4o3*inv_k*(fn1+2.*fn2)/(1.+fn1+fn2)
          intfb =c4o3*inv_k - fb1  ! intfb = \int^D gb(z') dz'
          fb0   =1./intfb          ! 1/[\int^D gb(z') dz']
#    elif defined FB_FUNC2
          fn1   =exp(-2.*kbrk*Dstp(i,j))
          fb1   =2.*inv_k*fn1/(1.+fn1)
          intfb =inv_k - fb1
          fb0   =1./intfb
#    elif defined FB_FUNC3
          fb1   =0.                ! fb1 = \int gb dz at z = -h
          fb0   =1./(1.-exp(-2.*kbrk*Dstp(i,j)))
#    endif

          do k=1,N,+1              ! <-- irreversible
#    ifdef FB_FUNC1
            dd  =z_w(i,j,N)-z_w(i,j,k)
            fn1 =3.*(exp(-2.*kbrk*dd)+exp(-4.*kbrk*dd))
            fn2 =exp(-6.*kbrk*dd)
            fb2 =c4o3*inv_k*(fn1+2.*fn2)/(1.+fn1+fn2)
#    elif defined FB_FUNC2
            fn1 =exp(-2.*kbrk*(z_w(i,j,N)-z_w(i,j,k)))
            fb2 =2.*inv_k*fn1/(1.+fn1)
#    elif defined FB_FUNC3
            fb2 =exp( kbrk*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &          -exp(-kbrk*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))
#    endif
            fb =fb0*(fb2-fb1)/Hz(i,j,k)
            brk_r(i,j,k) = fb*brk(i,j)
            fb1=fb2     ! recursive procedure
          enddo
#   endif       /* ifdef FB_FUNC0 */
        enddo
      enddo

      do k=1,N
        do j=jstrR,jendR
          do i=istr,iendR
            brk3dx(i,j,k)=0.5*( brk_r(i-1,j,k)*wdrx(i-1,j)
     &                         +brk_r(i  ,j,k)*wdrx(i  ,j) )
     &                            SWITCH umask(i,j)
     &                            SWITCH_WET umask_wet(i,j)
     &                            SWITCH_WAVEDRY umask_wavewet(i,j)
          enddo
        enddo

        do j=jstr,jendR
          do i=istrR,iendR
            brk3de(i,j,k)=0.5*( brk_r(i,j-1,k)*wdre(i,j-1)
     &                         +brk_r(i,j  ,k)*wdre(i,j  ) )
     &                            SWITCH vmask(i,j)
     &                            SWITCH_WET vmask_wet(i,j)
     &                            SWITCH_WAVEDRY vmask_wavewet(i,j)
          enddo
        enddo
      enddo       ! <-- discard brk_r, but keep ust_r

#   undef brk_r
#  endif        /* !WAVE_SFC_BREAK */
#  undef inv_d
!
! 3D bottom streaming acceleration term defined as a body force.
! ==============================================================
! The effect of bottom streaming in momentum balance is accounted for 
! by using wave dissipation due to bottom friction frc_r (depending
! on bed roughness) with an upward decaying vertical distribution fb:
!
!   fb ~ cosh(kfrc*(zeta-z)) (normalized so that integral is 1)
!
! kfrc is the inverse of decaying scale varying with the wave boundary
! layer thickness delta, estimated based on bed roughness.
!
! Hyperbolic functions (sinh, cosh, and tanh) in fb
! are expanded to exponential function so as not to be singular
! when kD is very large.
!
#  ifdef WAVE_STREAMING
#   define dfrac wrk1
#   define frc_r wrk2
#   define frc2d stk
!
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
          ks=30.*Zob_my(i,j)
#   ifdef OW_COUPLING_FULL
          abot=ubr(i,j)/max(fr(i,j),eps)
#   else          
          abot=wh(i,j)/max(2.*sinh(min(kD(i,j),khmax)),eps)
#   endif          
          delta=0.09*ks*(abot/ks)**0.82                ! WBL thickness
          kfrc=1./min(max(a_frc*delta,eps),Dstp(i,j))  ! friction scale
          frc2d(i,j)=0.
          fb1=1. -exp(-2.*kfrc*Dstp(i,j))
          fb0=1./(exp(-2.*kfrc*Dstp(i,j))-1.)
          do k=1,N
            fb2 =exp( kfrc*(z_w(i,j,N)-z_w(i,j,k)-Dstp(i,j)))
     &          -exp(-kfrc*(z_w(i,j,N)-z_w(i,j,k)+Dstp(i,j)))
            fb =fb0*(fb2-fb1)/Hz(i,j,k)
            fb1=fb2
            frc_r(i,j,k)=fb*frc(i,j)
            frc2d(i,j)=frc2d(i,j) + frc_r(i,j,k)*Hz(i,j,k)
          enddo
        enddo
      enddo
!
      do j=jstrR,jendR
        do i=istr,iendR
          frc2dx(i,j)=0.5*( frc2d(i-1,j)*wdrx(i-1,j)
     &                     +frc2d(i  ,j)*wdrx(i  ,j) )
     &                            SWITCH umask(i,j)
     &                            SWITCH_WET umask_wet(i,j)
     &                            SWITCH_WAVEDRY umask_wavewet(i,j)
          do k=1,N
            frc3dx(i,j,k)=0.5*( frc_r(i-1,j,k)*wdrx(i-1,j)
     &                         +frc_r(i  ,j,k)*wdrx(i  ,j) )
     &                            SWITCH umask(i,j)
     &                            SWITCH_WET umask_wet(i,j)
     &                            SWITCH_WAVEDRY umask_wavewet(i,j)
          enddo
        enddo
      enddo
      do j=jstr,jendR
        do i=istrR,iendR
          frc2de(i,j)=0.5*( frc2d(i,j-1)*wdre(i,j-1)
     &                     +frc2d(i,j  )*wdre(i,j  ) )
     &                            SWITCH vmask(i,j)
     &                            SWITCH_WET vmask_wet(i,j)
     &                            SWITCH_WAVEDRY vmask_wavewet(i,j)
          do k=1,N
            frc3de(i,j,k)=0.5*( frc_r(i,j-1,k)*wdre(i,j-1)
     &                         +frc_r(i,j  ,k)*wdre(i,j  ) )
     &                            SWITCH vmask(i,j)
     &                            SWITCH_WET vmask_wet(i,j)
     &                            SWITCH_WAVEDRY vmask_wavewet(i,j)
          enddo
        enddo
      enddo       ! <-- discard frc_r, frc2d
#   undef dfrac
#   undef frc_r
#   undef frc2d
#  endif   /* WAVE_STREAMING */

# endif       /* SOLVE3D */

# if (defined WAVE_FRICTION && !defined WKB_WWAVE) \
                        || defined WAVE_STREAMING
#   undef Zob
# endif
!
! Open boundary condition (valid only for idealized situations)
! =============================================================
!
# if defined ANA_BRY && !defined REGIONAL && !defined COASTAL
#  ifndef SHOREFACE
      z_tide=0.0 !wkb_tide
#  else
      z_tide=0.0
#  endif
#  if defined OBC_WEST && (defined Z_FRC_BRY || defined M2_FRC_BRY \
                                             || defined M3_FRC_BRY)
      if (WESTERN_EDGE) then
        do j=jstrR,jendR
#   ifdef Z_FRC_BRY
          zetabry_west(j)=z_tide + sup(istrR,j)
#   endif
#   ifdef M2_FRC_BRY
#    ifdef SHOREFACE
          ubarbry_west(j)=ubar(istr,j,knew)
#    else
          ubarbry_west(j)=-ust2d(istrR,j)    ! Stokes drift (xi)
#    endif
          vbarbry_west(j)=vbar(istr,j,knew) ! Stokes drift (eta)
#   endif
#   if defined SOLVE3D && defined M3_FRC_BRY
          do k=1,N
            ubry_west(j,k)=u(istr,j,k,nnew)
c            ubry_west(j,k)=-ust(istrR,j,k)
            vbry_west(j,k)=v(istr,j,k,nnew)
          enddo
#   endif
        enddo
      endif
#  endif  /* OBC_WEST */
#  if defined OBC_EAST && (defined Z_FRC_BRY || defined M2_FRC_BRY \
                                             || defined M3_FRC_BRY)
      if (EASTERN_EDGE) then
        do j=jstrR,jendR
#   ifdef Z_FRC_BRY
          zetabry_east(j)=z_tide + sup(iendR,j)
#   endif
#   ifdef M2_FRC_BRY
c         ubarbry_east(j)=ubar(iend,j,knew) ! Neumann
          ubarbry_east(j)=-ust2d(iend,j)  ! Stokes drift (xi)
          vbarbry_east(j)=vbar(iend,j,knew) !-vst2d(iendR,j)  ! Stokes drift (eta)
#   endif
#   if defined SOLVE3D && defined M3_FRC_BRY
          do k=1,N
            ubry_east(j,k)=u(iend,j,k,nnew)
c           ubry_east(j,k)=-ust(iendR,j,k)
            vbry_east(j,k)=v(iend,j,k,nnew)
          enddo
#   endif
        enddo
      endif
#  endif  /* OBC_EAST */
# endif  /* ANA_BRY */

# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile(istr,iend,jstr,jend,sup)
      call exchange_u2d_tile(istr,iend,jstr,jend,ust2d)
      call exchange_v2d_tile(istr,iend,jstr,jend,vst2d)
      call exchange_u2d_tile(istr,iend,jstr,jend,brk2dx)
      call exchange_v2d_tile(istr,iend,jstr,jend,brk2de)
      call exchange_u2d_tile(istr,iend,jstr,jend,frc2dx)
      call exchange_v2d_tile(istr,iend,jstr,jend,frc2de)
      call exchange_r2d_tile(istr,iend,jstr,jend,wepb)
      call exchange_r2d_tile(istr,iend,jstr,jend,wepd)
#  ifdef WAVE_ROLLER
      call exchange_r2d_tile(istr,iend,jstr,jend,wepr)
#  endif
#  ifdef SOLVE3D
      call exchange_r2d_tile(istr,iend,jstr,jend,Kapsrf)
      call exchange_r2d_tile(istr,iend,jstr,jend,calP)
      call exchange_u3d_tile(istr,iend,jstr,jend,ust)
      call exchange_v3d_tile(istr,iend,jstr,jend,vst)
      call exchange_r3d_tile(istr,iend,jstr,jend,wst)
      call exchange_r3d_tile(istr,iend,jstr,jend,kvf) 
#   ifdef LMD_MIXING
      call exchange_w3d_tile(istr,iend,jstr,jend,Akb)
      call exchange_w3d_tile(istr,iend,jstr,jend,Akw)
      call exchange_w3d_tile(istr,iend,jstr,jend,E_pre)
#   endif
#   ifndef WAVE_SFC_BREAK
      call exchange_u3d_tile (istr,iend,jstr,jend,brk3dx)
      call exchange_v3d_tile (istr,iend,jstr,jend,brk3de)
#   endif
#   ifdef WAVE_STREAMING
      call exchange_u3d_tile (istr,iend,jstr,jend,frc3dx)
      call exchange_v3d_tile (istr,iend,jstr,jend,frc3de)
#   endif
#  endif
# endif

      return
      end
!
!======================================================================
!
!                         WSTOKES
!
!======================================================================
!
# ifdef SOLVE3D

! Vertical Stokes drift velocity at rho-point
! ===========================================
! Note "wst" is not used in prognostic computation, but only 
! for output. 
!
      subroutine wstokes (tile)
      implicit none
#  include "param.h"
      integer tile, trd, omp_get_thread_num
#  include "private_scratch.h"
#  include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call wstokes_tile (istr,iend,jstr,jend, A2d(1,1,trd))
      return
      end

      subroutine wstokes_tile (istr,iend,jstr,jend,Wrk)
      implicit none
#  include "param.h"
      integer istr,iend,jstr,jend, i,j,k
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: Wrk

#  include "grid.h"
#  include "ocean3d.h"
#  include "forces.h"
#  include "scalars.h"

#  include "compute_auxiliary_bounds.h"

#  ifdef STOKES_DRIFT

      do j=jstr,jend
        do i=istr,iend
          Wrk(i,0)=0.
        enddo
        do k=1,N,+1        !--> recursive
          do i=istr,iend
            Wrk(i,k)=Wrk(i,k-1)-0.5*pm(i,j)*pn(i,j)*(
     &       (Hz(i+1,j  ,k)+Hz(i  ,j  ,k))*on_u(i+1,j  )*ust(i+1,j  ,k)
     &      -(Hz(i  ,j  ,k)+Hz(i-1,j  ,k))*on_u(i  ,j  )*ust(i  ,j  ,k)
     &      +(Hz(i  ,j+1,k)+Hz(i  ,j  ,k))*om_v(i  ,j+1)*vst(i  ,j+1,k)
     &      -(Hz(i  ,j  ,k)+Hz(i  ,j-1,k))*om_v(i  ,j  )*vst(i  ,j  ,k))
          enddo
        enddo
        do i=istr,iend
          wst(i,j,N)=+0.375*Wrk(i,N) +0.75*Wrk(i,N-1)
     &                              -0.125*Wrk(i,N-2)
        enddo
        do k=N-1,2,-1
          do i=istr,iend
            wst(i,j,k)=+0.5625*(Wrk(i,k  )+Wrk(i,k-1))
     &                 -0.0625*(Wrk(i,k+1)+Wrk(i,k-2))
          enddo
        enddo
        do i=istr,iend
          wst(i,j,  1)= -0.125*Wrk(i,2) +0.75*Wrk(i,1)
     &                                 +0.375*Wrk(i,0)
        enddo
      enddo

#   ifndef EW_PERIODIC
      if (WESTERN_EDGE) then                       ! Set lateral
        do k=1,N                                   ! boundary
          do j=jstr,jend                           ! conditions
            wst(istr-1,j,k)=wst(istr,j,k)
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do k=1,N
          do j=jstr,jend
            wst(iend+1,j,k)=wst(iend,j,k)
          enddo
        enddo
      endif
#   endif
#   ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do k=1,N
          do i=istr,iend
            wst(i,jstr-1,k)=wst(i,jstr,k)
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do k=1,N
          do i=istr,iend
            wst(i,jend+1,k)=wst(i,jend,k)
          enddo
        enddo
      endif
#    ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          wst(istr-1,jstr-1,k)=wst(istr,jstr,k)
        enddo
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          wst(istr-1,jend+1,k)=wst(istr,jend,k)
        enddo
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          wst(iend+1, jstr-1,k)=wst(iend,jstr,k)
        enddo
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          wst(iend+1,jend+1,k)=wst(iend,jend,k)
        enddo
      endif
#    endif
#   endif

#  else

      do k=1,N
        do j=jstr-1,jend+1
          do i=istr-1,iend+1
            wst(i,j,k)=0.
          enddo
        enddo
      enddo

#  endif  /* STOKES_DRIFT */
      return
      end
# endif  /* SOLVE3D */

#else
      subroutine mrl_wci_empty
      end
#endif

