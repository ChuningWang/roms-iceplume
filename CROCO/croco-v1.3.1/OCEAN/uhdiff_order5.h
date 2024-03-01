!
!===============================================================
!
! Compute diffusive part of UP5
!
!===============================================================
!
#  ifdef NS_PERIODIC
          jmin=1
          jmax=LOCALMM+1
#  else
#   ifdef MPI
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
#   else
          jmin=3
          jmax=Mm-1
#   endif
#  endif
#  ifdef EW_PERIODIC
          imin=0
          imax=LOCALLM+1
#  else
#   ifdef MPI
          if (WEST_INTER) then
            imin=0
          else
            imin=3
          endif
          if (EAST_INTER) then
            imax=Lmmpi+1
          else
            imax=Lmmpi-2
          endif
#   else
          imin=3
          imax=Lm-2
#   endif
#  endif
!
!----------------------------------------------------------------------
!  j loop: UFe4
!----------------------------------------------------------------------
!
          DO j = Jstr,Jend+1  !j_loop_y_flux_5
                                                  !
            IF ( j.ge.jmin .and. j.le.jmax ) THEN ! use full stencil
                                                  !
              DO i = IstrU,Iend
                if ( i.ge.imin .and. i.le.imax ) then
                  vel = flux6(Hvom(i-3,j,k),Hvom(i-2,j,k),Hvom(i-1,j,k),
     &                        Hvom(i  ,j,k),Hvom(i+1,j,k),Hvom(i+2,j,k),1.)
                else
                  vel = 0.5*(Hvom(i-1,j,k)+Hvom(i,j,k))
                endif
                flx5 = vel*FLUX5(
     &             u(i,j-3,k,nrhs), u(i,j-2,k,nrhs), 
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs),
     &             u(i,j+1,k,nrhs), u(i,j+2,k,nrhs), vel )
#  ifdef MASKING
!                vel = flux4(Hvom(i-2,j,k),Hvom(i-1,j,k),
!     &                      Hvom(i  ,j,k),Hvom(i+1,j,k),1.)
                vel = 0.5*(Hvom(i-1,j,k)+Hvom(i,j,k))
                flx3 = vel*FLUX3(
     &             u(i,j-2,k,nrhs), u(i,j-1,k,nrhs),
     &             u(i,j  ,k,nrhs), u(i,j+1,k,nrhs), vel ) 
!                vel = 0.5*(Hvom(i-1,j,k)+Hvom(i,j,k))
                flx2 = vel*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), vel, cdif)
#   ifdef UP5_MASKING
                mask2=umask(i,j-2)*umask(i,j+1)
                IF (vel.gt.0) THEN
                  mask1=umask(i,j-2)
                  mask0=umask(i,j-3)*mask2          
                ELSE
                  mask1=umask(i,j+1)
                  mask0=umask(i,j+2)*mask2
                ENDIF
#   else
                mask1=umask(i,j-2)*umask(i,j+1)
                mask2=umask(i,j-3)*umask(i,j+2)
                mask0=mask1*mask2
#   endif
                UFe4(i,j)=mask0*flx5+(1-mask0)*mask1*flx3+
     &                              (1-mask0)*(1-mask1)*flx2

#  else
                UFe4(i,j)=flx5
#  endif /* MASKING */
              ENDDO
                                           !
            ELSE IF ( j.eq.jmin-2 ) THEN   ! 2nd order flux next to south
                                           ! boundary
              DO i = IstrU,Iend
                vel = 0.5*(Hvom(i-1,j,k)+Hvom(i,j,k))
                UFe4(i,j) = vel*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), vel, cdif)
              ENDDO
                                                             !
            ELSE IF ( j.eq.jmin-1 .and. jmax.ge.jmin ) THEN  ! 3rd of 4th order flux 2 in
                                                             ! from south boundary
              DO i = IstrU,Iend
               vel = 0.5*(Hvom(i-1,j,k)+ Hvom(i,j,k))
!                vel = flux4(Hvom(i-2,j,k),Hvom(i-1,j,k),
!     &                      Hvom(i  ,j,k),Hvom(i+1,j,k),1.)
                flx3 = vel*FLUX3(
     &             u(i,j-2,k,nrhs), u(i,j-1,k,nrhs),
     &             u(i,j  ,k,nrhs), u(i,j+1,k,nrhs), vel )
#  ifdef MASKING
!                vel = 0.5*(Hvom(i-1,j,k)+ Hvom(i,j,k))
                flx2 = vel*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), vel, cdif)
                mask1=umask(i,j-2)*umask(i,j+1)
                UFe4(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                UFe4(i,j)=flx3
#  endif
              ENDDO

            ELSE IF ( j.eq.jmax+2 ) THEN  ! 2nd order flux next to north
                                          ! boundary
              DO i = IstrU,Iend
                vel = 0.5*(Hvom(i-1,j,k)+Hvom(i,j,k))
                UFe4(i,j) = vel*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), vel, cdif)
              ENDDO
                                          !
            ELSE IF ( j.eq.jmax+1 ) THEN  ! 3rd or 4th order flux 2 in from
                                          ! north boundary
              DO i = IstrU,Iend
                vel = 0.5*(Hvom(i-1,j,k)+ Hvom(i,j,k))
!               vel = flux4(Hvom(i-2,j,k),Hvom(i-1,j,k),
!     &                     Hvom(i  ,j,k),Hvom(i+1,j,k),1.)
                flx3 = vel*FLUX3(
     &             u(i,j-2,k,nrhs), u(i,j-1,k,nrhs),
     &             u(i,j  ,k,nrhs), u(i,j+1,k,nrhs), vel )
#  ifdef MASKING
!                vel = 0.5*(Hvom(i-1,j,k)+ Hvom(i,j,k))
                flx2 = vel*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), vel, cdif)
                mask1=umask(i,j-2)*umask(i,j+1)
                UFe4(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                UFe4(i,j)=flx3
#  endif
              ENDDO
            ENDIF
          ENDDO ! j_loop_y_flux_5
!
!----------------------------------------------------------------------
!  i loop: UFx4
!----------------------------------------------------------------------
!
          DO i = IstrU-1,Iend  !i_loop_x_flux_5
                                                  !
            IF ( i.ge.imin .and. i.le.imax ) THEN ! use full stencil
                                                  !
              DO j = Jstr,Jend
                vel = flux6(Huon(i-2,j,k),Huon(i-1,j,k),Huon(i  ,j,k),
     &                      Huon(i+1,j,k),Huon(i+2,j,k),Huon(i+3,j,k),1.)
!               vel = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                flx5 = vel*FLUX5(
     &             u(i-2,j,k,nrhs), u(i-1,j,k,nrhs),
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs),
     &             u(i+2,j,k,nrhs), u(i+3,j,k,nrhs), vel )
#  ifdef MASKING
                vel = flux4(Huon(i-1,j,k),Huon(i  ,j,k),
     &                      Huon(i+1,j,k),Huon(i+2,j,k),1.)
                flx3 = vel*FLUX3(
     &             u(i-1,j,k,nrhs), u(i  ,j,k,nrhs),
     &             u(i+1,j,k,nrhs), u(i+2,j,k,nrhs), vel )       
                vel = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                flx2 = vel*FLUX2(
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs), vel, cdif)
#   ifdef UP5_MASKING
                mask2=umask(i-1,j)*umask(i+2,j)
                IF (vel.gt.0) THEN
                  mask1=umask(i-1,j)
                  mask0=umask(i-2,j)*mask2          
                ELSE
                  mask1=umask(i+2,j)
                  mask0=umask(i+3,j)*mask2
                ENDIF
#   else
                mask1=umask(i-1,j)*umask(i+2,j)
                mask2=umask(i-2,j)*umask(i+3,j)
                mask0=mask1*mask2
#   endif
                UFx4(i,j)=mask0*flx5+(1-mask0)*mask1*flx3+
     &                              (1-mask0)*(1-mask1)*flx2
#  else
                UFx4(i,j)=flx5
#  endif /* MASKING */
              ENDDO
                                           !
            ELSE IF ( i.eq.imin-2 ) THEN   ! 2nd order flux next to south
                                           ! boundary
              DO j = Jstr,Jend
                vel = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                UFx4(i,j) = vel*FLUX2(
     &             u(i,j,k,nrhs), u(i+1,j,k,nrhs), vel, cdif)
              ENDDO
                                                             !
            ELSE IF ( i.eq.imin-1 .and. imax.ge.imin ) THEN  ! 3rd of 4th order flux 2 in
                                                             ! from south boundary
              DO j = Jstr,Jend
!               vel = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                vel = flux4(Huon(i-1,j,k),Huon(i  ,j,k),
     &                      Huon(i+1,j,k),Huon(i+2,j,k),1.)    
                flx3 = vel*FLUX3(
     &             u(i-1,j,k,nrhs), u(i  ,j,k,nrhs),
     &             u(i+1,j,k,nrhs), u(i+2,j,k,nrhs), vel )
#  ifdef MASKING
                vel = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                flx2 = vel*FLUX2(
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs), vel, cdif)
                mask1=umask(i-1,j)*umask(i+2,j)
                UFx4(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                UFx4(i,j)=flx3
#  endif
              ENDDO
                                          !
            ELSE IF ( i.eq.imax+2 ) THEN  ! 2nd order flux next to north
                                          ! boundary
              DO j = Jstr,Jend
                vel = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                UFx4(i,j) = vel*FLUX2(
     &             u(i,j,k,nrhs), u(i+1,j,k,nrhs), vel, cdif)
              ENDDO
                                          !
            ELSE IF ( i.eq.imax+1 ) THEN  ! 3rd or 4th order flux 2 in from
                                          ! north boundary
              DO j = Jstr,Jend
!               vel = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                vel = flux4(Huon(i-1,j,k),Huon(i  ,j,k),
     &                      Huon(i+1,j,k),Huon(i+2,j,k),1.)  
                flx3 = vel*FLUX3(
     &             u(i-1,j,k,nrhs), u(i  ,j,k,nrhs),
     &             u(i+1,j,k,nrhs), u(i+2,j,k,nrhs),  vel )
#  ifdef MASKING
                vel = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                flx2 = vel*FLUX2(
     &             u(i,j,k,nrhs), u(i+1,j,k,nrhs), vel, cdif)
                mask1=umask(i-1,j)*umask(i+2,j)
                UFx4(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                UFx4(i,j)=flx3
#  endif
              ENDDO
            ENDIF
          ENDDO ! i_loop_x_flux_5

