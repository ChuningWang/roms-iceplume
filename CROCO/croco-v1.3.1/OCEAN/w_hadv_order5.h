!
!===============================================================
!
! Compute wz 5th order horizontal advection
!
!===============================================================
!
!----------------------------------------------------------------------
!  j loop: WFe
!----------------------------------------------------------------------
!
          DO j = Jstr,Jend+1  !j_loop_y_flux_5
                                                  !
            IF ( j.ge.jmin .and. j.le.jmax ) THEN ! use full stencil
                                                  !
              DO i = Istr,Iend
                vel = Hvom_w(i,j,k)
                flx5 = vel*FLUX5(
     &             wz(i,j-3,k,nrhs), wz(i,j-2,k,nrhs), 
     &             wz(i,j-1,k,nrhs), wz(i,j  ,k,nrhs),
     &             wz(i,j+1,k,nrhs), wz(i,j+2,k,nrhs),  vel )
#  ifdef MASKING 
                flx3 = vel*FLUX3(
     &             wz(i,j-2,k,nrhs), wz(i,j-1,k,nrhs),
     &             wz(i,j  ,k,nrhs), wz(i,j+1,k,nrhs),  vel ) 
                flx2 = vel*FLUX2(
     &             wz(i,j-1,k,nrhs), wz(i,j,k,nrhs), vel, cdif)
#   ifdef UP5_MASKING
                mask0=rmask(i,j-1)*rmask(i,j)
                mask2=rmask(i,j-2)*mask0*rmask(i,j+1)
                IF (vel.gt.0) THEN
                  mask1=rmask(i,j-2)*mask0
                  mask3=rmask(i,j-3)*mask2          
                ELSE
                  mask1=rmask(i,j+1)*mask0
                  mask3=rmask(i,j+2)*mask2
                ENDIF
                WFe(i,j)=mask3*flx5+(1-mask3)*mask1*flx3+
     &                             (1-mask3)*(1-mask1)*mask0*flx2
#   else
                mask1=rmask(i,j-2)*rmask(i,j+1)
                mask2=rmask(i,j-3)*rmask(i,j+2)
                mask0=mask1*mask2
                WFe(i,j)=mask0*flx5+(1-mask0)*mask1*flx3+
     &                         (1-mask0)*(1-mask1)*flx2
#   endif /* UP5_MASKING */
#  else
                WFe(i,j)=flx5
#  endif /* MASKING */
              ENDDO
                                           !
            ELSE IF ( j.eq.jmin-2 ) THEN   ! 2nd order flux next to south
                                           ! boundary
              DO i = Istr,Iend
                vel = Hvom_w(i,j,k)
                WFe(i,j) = vel*FLUX2(
     &             wz(i,j-1,k,nrhs), wz(i,j,k,nrhs), vel, cdif)
              ENDDO
                                                             !
            ELSE IF ( j.eq.jmin-1 .and. jmax.ge.jmin ) THEN  ! 3rd of 4th order flux 2 in
                                                             ! from south boundary
              DO i = Istr,Iend
                vel = Hvom_w(i,j,k)
                flx3 = vel*FLUX3(
     &             wz(i,j-2,k,nrhs), wz(i,j-1,k,nrhs),
     &             wz(i,j  ,k,nrhs), wz(i,j+1,k,nrhs),  vel )
#  ifdef MASKING
                flx2 = vel*FLUX2(
     &             wz(i,j-1,k,nrhs), wz(i,j,k,nrhs), vel, cdif)
                mask1=rmask(i,j-2)*rmask(i,j+1)
                WFe(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                WFe(i,j)=flx3
#  endif
              ENDDO
                                          !
            ELSE IF ( j.eq.jmax+2 ) THEN  ! 2nd order flux next to north
                                          ! boundary
              DO i = Istr,Iend
                vel = Hvom_w(i,j,k)
                WFe(i,j) = vel*FLUX2(
     &             wz(i,j-1,k,nrhs), wz(i,j,k,nrhs), vel, cdif)
              ENDDO
                                          !
            ELSE IF ( j.eq.jmax+1 ) THEN  ! 3rd or 4th order flux 2 in from
                                          ! north boundary
              DO i = Istr,Iend
                vel = Hvom_w(i,j,k)
                flx3 = vel*FLUX3(
     &             wz(i,j-2,k,nrhs), wz(i,j-1,k,nrhs),
     &             wz(i,j  ,k,nrhs), wz(i,j+1,k,nrhs),  vel )
#  ifdef MASKING
                flx2 = vel*FLUX2(
     &             wz(i,j-1,k,nrhs), wz(i,j,k,nrhs), vel, cdif)
                mask1=rmask(i,j-2)*rmask(i,j+1)
                WFe(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                WFe(i,j)=flx3
#  endif
              ENDDO
            ENDIF
          ENDDO ! j_loop_y_flux_5
!
!----------------------------------------------------------------------
!  i loop: WFx
!----------------------------------------------------------------------
!
          DO i = Istr,Iend+1  !i_loop_x_flux_5
                                                  !
            IF ( i.ge.imin .and. i.le.imax ) THEN ! use full stencil
                                                  !
              DO j = Jstr,Jend
                vel = Huon_w(i,j,k)
                flx5 = vel*FLUX5(
     &             wz(i-3,j,k,nrhs), wz(i-2,j,k,nrhs),
     &             wz(i-1,j,k,nrhs), wz(i  ,j,k,nrhs),
     &             wz(i+1,j,k,nrhs), wz(i+2,j,k,nrhs),  vel )
#  ifdef MASKING
                flx3 = vel*FLUX3(
     &             wz(i-2,j,k,nrhs), wz(i-1,j,k,nrhs),
     &             wz(i  ,j,k,nrhs), wz(i+1,j,k,nrhs),  vel )
                flx2 = vel*FLUX2(
     &             wz(i-1,j,k,nrhs), wz(i,j,k,nrhs), vel, cdif)
#   ifdef UP5_MASKING
                mask0=rmask(i-1,j)*rmask(i,j)
                mask2=rmask(i-2,j)*mask0*rmask(i+1,j)
                IF (vel.gt.0) THEN
                  mask1=rmask(i-2,j)*mask0
                  mask3=rmask(i-3,j)*mask2          
                ELSE
                  mask1=rmask(i+1,j)*mask0
                  mask3=rmask(i+2,j)*mask2
                ENDIF
                WFx(i,j)=mask3*flx5+(1-mask3)*mask1*flx3+
     &                             (1-mask3)*(1-mask1)*mask0*flx2
#   else
                mask1=rmask(i-2,j)*rmask(i+1,j)
                mask2=rmask(i-3,j)*rmask(i+2,j)
                mask0=mask1*mask2
                WFx(i,j)=mask0*flx5+(1-mask0)*mask1*flx3+
     &                         (1-mask0)*(1-mask1)*flx2
#   endif /* UP5_MASKING */
#  else
                WFx(i,j)=flx5
#  endif /* MASKING */
              ENDDO
                                           !
            ELSE IF ( i.eq.imin-2 ) THEN   ! 2nd order flux next to south
                                           ! boundary
              DO j = Jstr,Jend
                vel = Huon_w(i,j,k)
                WFx(i,j) = vel*FLUX2(
     &             wz(i-1,j,k,nrhs), wz(i,j,k,nrhs), vel, cdif)
              ENDDO
                                                             !
            ELSE IF ( i.eq.imin-1 .and. imax.ge.imin ) THEN  ! 3rd of 4th order flux 2 in
                                                             ! from south boundary
              DO j = Jstr,Jend
                vel = Huon_w(i,j,k)
                flx3 = vel*FLUX3(
     &             wz(i-2,j,k,nrhs), wz(i-1,j,k,nrhs),
     &             wz(i  ,j,k,nrhs), wz(i+1,j,k,nrhs),  vel )
#  ifdef MASKING
                flx2 = vel*FLUX2(
     &             wz(i-1,j,k,nrhs), wz(i,j,k,nrhs), vel, cdif)
                mask1=rmask(i-2,j)*rmask(i+1,j)
                WFx(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                WFx(i,j)=flx3
#  endif
              ENDDO
                                          !
            ELSE IF ( i.eq.imax+2 ) THEN  ! 2nd order flux next to north
                                          ! boundary
              DO j = Jstr,Jend
                vel = Huon_w(i,j,k)
                WFx(i,j) = vel*FLUX2(
     &             wz(i-1,j,k,nrhs), wz(i,j,k,nrhs), vel, cdif)
              ENDDO
                                          !
            ELSE IF ( i.eq.imax+1 ) THEN  ! 3rd or 4th order flux 2 in from
                                          ! north boundary
              DO j = Jstr,Jend
                vel = Huon_w(i,j,k)
                flx3 = vel*FLUX3(
     &             wz(i-2,j,k,nrhs), wz(i-1,j,k,nrhs),
     &             wz(i  ,j,k,nrhs), wz(i+1,j,k,nrhs),  vel )
#  ifdef MASKING
                flx2 = vel*FLUX2(
     &             wz(i-1,j,k,nrhs), wz(i,j,k,nrhs), vel, cdif)
                mask1=rmask(i-2,j)*rmask(i+1,j)
                WFx(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                WFx(i,j)=flx3
#  endif
              ENDDO
            ENDIF
          ENDDO ! i_loop_x_flux_5


