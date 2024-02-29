!
!----------------------------------------------------------
! Compute wz vertical advective fluxes 
! using 5th-order WENO scheme
!----------------------------------------------------------
!
          do k=3,N-2
            do i=Istr,Iend
              FC(i,k)=We_r(i,k)*FLUX5(
     &             wz(i,j,k-3,nrhs), wz(i,j,k-2,nrhs), 
     &             wz(i,j,k-1,nrhs), wz(i,j,k  ,nrhs),
     &             wz(i,j,k+1,nrhs), wz(i,j,k+2,nrhs), We_r(i,k))
            enddo
          enddo
          do i=Istr,Iend
            FC(i,2)=We_r(i,2)*FLUX3(
     &           wz(i,j,0,nrhs), wz(i,j,1,nrhs), 
     &           wz(i,j,2,nrhs), wz(i,j,3,nrhs), We_r(i,2))

            FC(i,N-1)=We_r(i,N-1)*FLUX3(
     &           wz(i,j,N-3,nrhs), wz(i,j,N-2,nrhs), 
     &           wz(i,j,N-1,nrhs), wz(i,j,N  ,nrhs), We_r(i,N-1))

            FC(i,1)=We_r(i,1)*FLUX2( wz(i,j,0  ,nrhs),
     &                               wz(i,j,1  ,nrhs), We_r(i,1), 1.)

            FC(i,N)=We_r(i,N)*FLUX2( wz(i,j,N-1,nrhs),
     &                               wz(i,j,N  ,nrhs), We_r(i,N), 1.)

            FC(i,0)=0.
            FC(i,N+1)=0.
          enddo

