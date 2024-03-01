!
!====================================================================
!                        Internal wave maker
!
! Mode-1 internal wave (e.g., Hall et al., JPO 2013)
!====================================================================
!
!--------------------------------------------------------------------
!  Configurations
!--------------------------------------------------------------------
!
        wp=12.*3600             ! period
        wa=0.08                 ! amplitude m/s
!
!--------------------------------------------------------------------
!  Initialisation
!--------------------------------------------------------------------
!
        ramp=tanh(dt/wp*float(iic-ntstart))
        wa=wa*ramp
        wf=2*pi/wp
!
!--------------------------------------------------------------------
!  Sea level zetabry
!--------------------------------------------------------------------
!
#   ifdef Z_FRC_BRY
        do j=JstrR,JendR
          zetabry_west(j)=0.
        enddo
#   endif /* Z_FRC_BRY */
!
!--------------------------------------------------------------------
!  XI velocity components ubry and ubarbry
!--------------------------------------------------------------------
!
#   ifdef M3_FRC_BRY
        do j=JstrR,JendR
          h0=0.5*(h(0,j)+h(1,j))
          do k=1,N
            Zu=0.5*(z_r(0,j,k)+z_r(1,j,k))
            ubry_west(j,k)=wa*cos(pi*Zu/h0)*sin(wf*time)
          enddo
        enddo
#   endif /* M3_FRC_BRY */

#   ifdef M2_FRC_BRY
        do j=JstrR,JendR
          cff4=0.
          cff5=0.
          do k=1,N
            cff4=cff4+ubry_west(j,k)*(Hz(0,j,k)+Hz(1,j,k))
            cff5=cff5+(Hz(0,j,k)+Hz(1,j,k))
          enddo
          ubarbry_west(j)=cff4/cff5
        enddo
#   endif /* M2_FRC_BRY */
!
!--------------------------------------------------------------------
!  ETA velocity components vbry and vbarbry
!--------------------------------------------------------------------
!
#   ifdef M3_FRC_BRY
        do j=JstrV,JendR
          do k=1,N
#    ifdef UV_COR
            Zv=0.5*(z_r(0,j,k)+z_r(0,j-1,k))
            vbry_west(j,k)=wa*cos(pi*Zv/h0)*cos(wf*time)
     &                       *0.5*(f(0,j)+f(0,j-1))/wf
#    else
            vbry_west(j,k)=0.
#    endif
          enddo
        enddo
#   endif /* M3_FRC_BRY */
#   ifdef M2_FRC_BRY
#    ifdef UV_COR
        do j=JstrR,JendR
          cff4=0.
          cff5=0.
          do k=1,N
            cff4=cff4+vbry_west(j,k)*(Hz(0,j,k)+Hz(0,j-1,k))
            cff5=cff5+(Hz(0,j,k)+Hz(0,j-1,k))
          enddo
          vbarbry_west(j)=cff4/cff5
        enddo
#    else
        do j=JstrV,JendR
          vbarbry_west(j)=0.
        enddo
#    endif
#   endif /* M2_FRC_BRY */
!
!--------------------------------------------------------------------
!  Z velocity component wbry
!--------------------------------------------------------------------
!
#   ifdef W_FRC_BRY
        if (FIRST_TIME_STEP) then
          do k=1,N-1
            do j=JstrR,JendR
              bvf0bry_west(j,k)=bvf(0,j,k) ! init. stratif.
            enddo
          enddo
          do j=JstrR,JendR
            bvf0bry_west(j,N)=bvf(0,j,N-1)
          enddo
        endif
        cff2=wf**2
        do j=JstrR,JendR
          h0=h(0,j)
          do k=1,N
            Zr=z_w(0,j,k)
            cff3=f(0,j)**2                 !         f**2
            cff1=bvf0bry_west(j,k)         ! initial N**2
#    ifdef UV_COR
            wbry_west(j,k)=wa*sin(pi*Zr/h0)*cos(wf*time)
     &                       *sqrt((cff2-cff3)/(cff1-cff2))
#    else
            wbry_west(j,k)=wa*sin(pi*Zr/h0)*cos(wf*time)
     &                       *sqrt(cff2/(cff1-cff2))
#    endif
          enddo
        enddo
#   endif /* W_FRC_BRY */
!
!--------------------------------------------------------------------
!  TRACERS tbry
!--------------------------------------------------------------------
!
#   ifdef T_FRC_BRY 
        if (FIRST_TIME_STEP) then
          do k=1,N
            do j=JstrR,JendR
              do itrc=1,NT
                tbry_west(j,k,itrc)=t(0,j,k,1,itrc)
              enddo
            enddo
          enddo
        endif
#    if !defined NONLIN_EOS && !defined SALINITY \
                            &&  defined TEMPERATURE
        if (FIRST_TIME_STEP) then
          do k=1,N
            do j=JstrR,JendR
              t0bry_west(j,k)=t(0,j,k,1,itemp) ! save initial
              bvf0bry_west(j,k)=bvf(0,j,k)     ! stratification
            enddo
          enddo
          do j=JstrR,JendR
            bvf0bry_west(j,0)=bvf(0,j,1)
            bvf0bry_west(j,N)=bvf(0,j,N-1)
          enddo
        endif
        cff2=(1./wf)**2  ! inverse frequency squared
        cff3=Tcoef/rho0  ! thermal expansion coefficient °C-1
        do j=JstrR,JendR
          h0=h(0,j)
          do k=1,N
            Zr=z_r(0,j,k)
            cff1=sqrt(0.5*(bvf0bry_west(j,k-1)+
     &                     bvf0bry_west(j,k  )))     ! init N
            tbry_west(j,k,itemp)= t0bry_west(j,k)
     &                            -wa*sin(pi*Zr/h0)*sin(wf*time)
     &                               *cff1/(g*cff3)
#     ifdef UV_COR
     &                               *sqrt(1.-cff2*f(0,j)**2)
#     endif
          enddo
        enddo
#    endif
#   endif
!
!--------------------------------------------------------------------
!  NBQ variables: unqbry, vnbqbru and wnbqbry
!--------------------------------------------------------------------
!
#   ifdef NBQ_FRC_BRY
        do k=1,N
          do j=JstrR,JendR
            unbqbry_west(j,k)=0.5*(Hz(0,j,k)+Hz(1,j,k))
     &                                  *ubry_west(j,k)
          enddo
          do j=JstrV,JendR
            vnbqbry_west(j,k)=0.5*(Hz(0,j,k)+Hz(0,j-1,k))
     &                                    *vbry_west(j,k)
          enddo
        enddo
        do k=1,N-1
          do j=JstrR,JendR
            wnbqbry_west(j,k)=0.5*(Hz(0,j,k)+Hz(0,j,k+1))
     &                                    *wbry_west(j,k)
          enddo
        enddo
        k=N
        do j=JstrR,JendR
          wnbqbry_west(j,k)=0.5*Hz(0,j,k)*wbry_west(j,k)
        enddo
#   endif

