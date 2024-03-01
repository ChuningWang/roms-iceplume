#include "cppdefs.h"

      module parameters

      integer,parameter :: rsh=4,rlg=8
      real,parameter    :: pi=3.141592654
      real(kind=rsh),parameter    :: valmanq=999.


      end module parameters 

!******************************************************
!
! Interpolations entres les 128 hauteurs precalculees
!
      subroutine ma7(flag,seci,hn,haut,k)
!
!******************************************************
!
      use parameters
      implicit none
     !
      logical :: flag
      integer                      :: idk,l,k
      real               :: d1,haut,fk,d,hh,d2 
      real,dimension(130):: hn
      real               :: seci 
!
      fk=seci/60.0/11.25+2.0
      k=fk+0.5
      d=fk-k
      if(k >= 130)then
       idk=k-129
       k=k-idk
       d=d+idk
      endif
      if (flag) write(6,*)"ma7 seci",seci
      if (flag) write(6,*)"ma7 k",k
      if (flag) write(6,*)"ma7 hn",hn
!
      l=k
      if(l == 130)l=129
      hh=hn(l)+hn(l)
      d1=hn(l+1)-hn(l-1)
      d2=hn(l+1)+hn(l-1)-hh
!
      haut=(hh+d*(d*d2+d1))/800.0
      if (flag) write(6,*)"ma7 haut",haut,l,hh,hn(l+1),hn(l-1),d1,d2
!
      end subroutine ma7
!
!***********************************************************
!
!  calcul de hn :128 hauteurs de maree par jour
!
      subroutine ma1(flag,my,ii,jj,t0,
     &               jour,mois,ia,hn,rr,gg,nno,fr1,fr2,nam)
!
!***********************************************************
!
      use parameters
      implicit none
      logical :: flag
      integer :: l,nj,i,julien2,nd0,n,m,jour,mois,ia,kg,ng,j,k
      integer :: ll,kk,ii,jj,my
      integer,dimension(8)              :: nam
      integer,dimension(30)             :: noa
      integer,dimension(30,8)           :: nno
      real                    :: t0,fq,co,so
      double precision        :: fnjd
      real,dimension(30)      :: f,q,q0,r1,r2,v0
      real,dimension(31)      :: g,r
      real,dimension(130)     :: hn
      real,dimension(362)     :: c
      real,dimension(91)      :: cs=(/ 
     &
     &  1.000000,0.999848,0.999391,0.998630,0.997564,0.996195,0.994522, 
     &  0.992546,0.990268,0.987688,0.984808,0.981627,0.978148,0.974370,
     &  0.970296,0.965926,0.961262,0.956305,0.951057,0.945519,0.939693,
     &  0.933580,0.927184,0.920505,0.913545,0.906308,0.898794,0.891007,
     &  0.882948,0.874620,0.866025,0.857167,0.848048,0.838671,0.829038,
     &  0.819152,0.809017,0.798636,0.788011,0.777146,0.766044,0.754710,
     &  0.743145,0.731354,0.719340,0.707107,0.694658,0.681998,0.669131,
     &  0.656059,0.642788,0.629320,0.615661,0.601815,0.587785,0.573576,
     &  0.559193,0.544639,0.529919,0.515038,0.500000,0.484810,0.469472,
     &  0.453991,0.438371,0.422618,0.406737,0.390731,0.374607,0.358368,
     &  0.342020,0.325568,0.309017,0.292372,0.275637,0.258819,0.241922,
     &  0.224951,0.207912,0.190809,0.173648,0.156434,0.139173,0.121869,
     &  0.104528,0.087156,0.069756,0.052336,0.034899,0.017452,0.000000/)
       real,dimension(30,8)    :: rr,gg,fr1,fr2
       real,dimension(11,3)    :: x,y
!
      equivalence (c(1),cs(1))
!
! executable statements
!-------------------------
!
      l=1
!
      x(:,:)=0.0
      y(:,:)=0.0
!
      c(271)=0.0
      do i=1,90
       c(182-i)=-c(i)
       c(362-i)=c(i)
       c(180+i)=-c(i)
      end do

      if (flag) write(6,*)"ma1 c",c(:)
!
      nj=julien2(ia,jour,mois)-2415021
      fnjd=real(nj+l,8)-DBLE(1.4921875)
      nj=mod(nj,7)+1
      if(nj <= 0)nj=nj+7
      nd0=-1
      n=0
      m=1
!
      if (flag) write(6,*)"ma1 nj",nj
      if (flag) write(6,*)"ma1 fnjd",fnjd

      do kg=1,8
       n=nam(kg)
       if(n /= 0) then
        do i=1,n
         r(i)=rr(i,kg)
         g(i)=gg(i,kg)
         noa(i)=nno(i,kg)
         r1(i)=fr1(i,kg)
         r2(i)=fr2(i,kg)
         if (flag)
     &  write(6,*)"ma1 before ma3",i,r(i),g(i),r1(i),r2(i),noa(i)
        end do
!


        call ma3(f,v0,q,n,r,g,fnjd,ng,t0,r1,r2,noa)
        
!
        do i=1,n
         r(i)=r(i)*f(i)
         q0(i)=v0(i)-g(i)+360.0
         q(i)=mod(q(i),360.0)
         q0(i)=mod(q0(i),360.0)
         if(q0(i) < 0.0)q0(i)=q0(i)+360.0
         if (flag) write(6,*)"ma1 after ma3",i,f(i),q0(i),v0(i),q(i),
     &                      r(i),g(i),ng,r1(i),r2(i),noa(i)
        end do
!
        do i=l,3
         do j=1,n
          fq=q0(j)
          k=fq+1
          fq=fq-k+1
          co=((c(k+1)-c(k))*fq+c(k))*r(j)
          if (flag) write(6,*)"x,y",j,k,fq,co
          k=k-90
          if(k <= 0) k=k+360
          so=((c(k+1)-c(k))*fq+c(k))*r(j)
          x(ng,i)=x(ng,i)+co
          y(ng,i)=y(ng,i)+so
          q0(j)=q0(j)+q(j)
          if(q0(j) > 360.0)q0(j)=q0(j)-360.0
          if (flag) write(6,*)"x,y2",j,so,x(ng,i),y(ng,i),q0
         end do
        end do
       end if
      end do
!
      if (flag) write(6,*)'ma1 before ma4',x,y
      call ma4(hn,l,x,y)
      if (flag) write(6,*)'ma1 after ma4',hn
!
      end subroutine ma1
!
!**********************************************************************
!
!  Calcul des hauteurs reduites h pour les jours jour-1,jour et jour+1
!  et appel de la subroutine de calcul des 128 hauteurs par jour (ma5)
!
      subroutine ma4(hn,l,hm1,hm2)
!
!**********************************************************************
!
      use parameters
      implicit none
      integer                        :: k,l,i
      real                 :: fnm
      real,dimension(130)  :: hn
      real,dimension(11,3) :: hm1,hm2
      real,dimension(128,3):: h
      complex,dimension(128)         :: x
!
! executable statements
!-------------------------
!
      do  k=l,3
       h(:,k)=0.0
       x(:)=0.0
       x(1)=hm1(1,k)
       fnm=x(1)*2.0
       hm2(1,k)=0.0
       hm1(1,k)=0.0
       do i=2,11
        x(i)=cmplx(hm1(i,k),hm2(i,k))
        hm1(i,k)=0.0
        hm2(i,k)=0.0
        x(130-i)=conjg(x(i))
       end do
!
       call fft(x)
!
       h(:,k)=x(:)
      end do
!
      call ma5(h,hn)
!
      end subroutine ma4
!
!**********************************************************************
!
!  Calcul des 128 hauteurs par interpolation entre les hauteurs reduites
!
      subroutine ma5(h,x)
!
!**********************************************************************
!
      use parameters
      implicit none
      integer                        :: k
      real                 :: b,hh
      real,dimension(128,3):: h
      real,dimension(130)  :: x
!
! executable statements
!-------------------------
!
      b=-0.5078125
      x(1)=x(129)
      x(2)=x(130)
      do k=1,128
       hh=h(k,2)+h(k,2)
       b=b+0.0078125
       x(k+2)=((h(k,1)+h(k,3)-hh)*b+h(k,3)-h(k,1))*b+hh
      end do
      if(x(1) == -valmanq)then
       x(2)=3.0*(x(3)-x(4))+x(5)
       x(1)=6.0*x(3)-8.0*x(4)+3.0*x(5)
      endif
!
      end subroutine ma5
!
!********************************************
!
!  Transformee de Fourier rapide
!
      subroutine fft(a)
!
!********************************************
! 
      use parameters
      implicit none
      integer                :: m,nv2,nm1,n,i,j,ip,l,k,le,le1
      real         :: fij
      complex,dimension(128) :: a
      complex                :: u,w,t
!
! executable statements
!-------------------------
!
      m=7
      nv2=64
      nm1=127
      n=128
      j=1
!
      do  i=1,nm1
      if(i < j) then
       t=a(j)
       a(j)=a(i)
       a(i)=t
      end if
      k=nv2
    6 continue
      if(k < j) then
       j=j-k
       k=k/2
       goto 6
      end if
      j=j+k
      end do
!
      do l=1,m
      le=2**l
      le1=le/2
      u=cmplx(1.0,0.0)
      fij=pi/le1
      w=cexp(cmplx(0.0,fij))
      do j=1,le1
       do i=j,n,le
        ip=i+le1
        t=a(ip)*u
        a(ip)=a(i)-t
        a(i)=a(i)+t
       end do
       u=u*w
      end do
      end do
!
      end subroutine fft
!
!**************************************************************
!
!  Calcul des amplitudes et phases des composantes
!
      subroutine ma3(f,v0,q,n,r,g,fnjd,ng,t0,r1,r2,noa)
!
!**************************************************************
!
      use parameters
       implicit none
       integer                     :: i,ng,n
       integer,dimension(6)        :: nd
       integer,dimension(30)       :: noa
       real              :: fa1,fb1,fa2,fb2,a1,a2,asfov,t0
       double precision            :: fnjd
       real,dimension(30):: q,v0,f,r,g,r1,r2
       double precision            :: dq
       complex                     :: fc
!
!----------------------------------------------
!
      do i=1,n
        nd(6)=noa(i)
        nd(5)=noa(i)/10
        nd(4)=noa(i)/100
        nd(3)=noa(i)/1000
        nd(2)=noa(i)/10000
        nd(1)=noa(i)/100000
        nd(6)=nd(6)-nd(5)*10
        nd(5)=nd(5)-nd(4)*10
        nd(4)=nd(4)-nd(3)*10
        nd(3)=nd(3)-nd(2)*10
        nd(2)=nd(2)-nd(1)*10
        ng=nd(1)+1
        if(nd(1) == 0) r(i)=r(i)*2.0
        call masfo(nd,dq)
        q(i)=dq
        g(i)=g(i)+q(i)*t0/24.0
        g(i)=mod(g(i),360.0)
        if(g(i) < 0.0) g(i)=g(i)+360.0
        v0(i)=asfov(nd)+mod(REAL(dq*fnjd,8),DBLE(360.))
        fc=1.0
        if( r1(i) /= 0.0 .or. r2(i) /= 0.0 ) then
          fb2=9.242202e-4
          fb1=-fb2
          fa2=1.760045
          fa1=4.523139
          if(noa(i) == 275555)then
           fb1=2*fb2
           fa1=2*fa2
          endif
          a1=fa1+fb1*fnjd
          a2=fa2+fb2*fnjd
          fc=1.0+r1(i)*cmplx(cos(a1),sin(a1))+
     &    r2(i)*cmplx(cos(a2),sin(a2))
        end if
!
        v0(i)=v0(i)+atan2(aimag(fc),real(fc))*57.29578
        f(i)=cabs(fc)
        if(noa(i) == 355555.or.noa(i) == 382555 .or.noa(i) == 164555)
     &          v0(i)=v0(i)-90.0
      end do
!
      end subroutine ma3
!
!***************************************************
!
! calcul des vitesses angulaires des composantes
!
      subroutine masfo(nd,asfo)
!
!***************************************************
!
      use parameters
      implicit none
      integer,dimension(6) :: nd
      double precision       :: asfo
!
      asfo =(DBLE(360.)-DBLE(12.19074939))*(nd(1))
     &     +DBLE(13.17639673)*(nd(2)-5)+DBLE(0.98564734)*(nd(3)-5) 
     &     +DBLE(0.11140408)*(nd(4)-5)+DBLE(0.05295392)*(nd(5)-5)+
     &      DBLE(0.00004707)*(nd(6)-5)
!
      end subroutine masfo
!
!***************************************************
!
! Calcul des arguments astronomiques
!
      function asfov(nd)
!
!***************************************************
!
      use parameters
      implicit none
      integer,dimension(6) :: nd
      real       :: fov,asfov
!
      fov=280.1895*(nd(1)+nd(3)-5)+
     &  277.0248*(nd(2)-nd(1)-5)+334.3853*(nd(4)-5)
     &  +100.8432*(nd(5)-5)+281.2209*(nd(6)-5)+(mod(nd(1),2))*90.0
      asfov=(mod(fov,360.0))
!
      end function asfov
!
!**************************************************************
!
! Calcul du nombre de jours ecoules depuis le 1er Janvier 1900
!
      function julien2(ia,jou,moi)
!
!*********************************************************
!
      use parameters
      implicit none 
      integer                 :: i,julien2,ia,jou,moi,ibs,jour,mois,iy,m
      integer,dimension(12,2),PARAMETER :: jo=
     &     reshape((/0,31,59,90,120,151,181,212,243,273,304,334,
     &     0,31,60,91,121,152,182,213,244,274,305,335/),(/12,2/))
      real          :: a,b
!
      jour=jou
      mois=moi
!
      if(mois==1) then
      ibs=1
      if( ia >= 1582 .and.(ia /= 1582 .or. jou >= 277) ) then
       ibs=mod(ia,4)+2
       if(ibs /= 2)ibs=1
      end if
       do i=1,12
        mois=12-i+1
        jour=jou-jo(mois,ibs)
        if(jour > 0) goto 20
       end do
   20 continue
      end if
!
      a=ia+mois/100.0+jour/10000.0
      b=0.0
      iy=ia
      m=mois
!
      if(mois <= 2)then
        iy=iy-1
        m=m+12
      endif
!
      if(a >= 1582.1015)then
       b=2-int(iy/100)+int(iy/100)/4
      endif
!
      julien2=int(365.25*iy)+int(30.6001*(m+1))+jour+1720995+b
!
      end function julien2
