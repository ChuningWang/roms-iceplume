!
!----------------------------------------------------------
! Compute vertical advective fluxes 
! using 5th-order WENO scheme
!----------------------------------------------------------
!
#  ifdef EW_PERIODIC
          imin=1
          imax=LOCALLM+1
#  else
#   ifdef MPI
          if (WEST_INTER) then
            imin=1
          else
            imin=3
          endif
          if (EAST_INTER) then
            imax=Lmmpi+1
          else
            imax=Lmmpi-1
          endif
#   else
          imin=3
          imax=Lm-1
#   endif
#  endif
!
!----------------------------------------------------------------------
!  k loop: FC
!----------------------------------------------------------------------
!
        do k=3,N-3
          do i=IstrU,Iend
            if ( i.ge.imin .and. i.le.imax ) then
              vel=flux6(We(i-3,j,k),We(i-2,j,k),We(i-1,j,k),
     &                  We(i  ,j,k),We(i+1,j,k),We(i+2,j,k),1.)
            else
              vel=0.5*(We(i-1,j,k)+We(i,j,k))
            endif
            FC(i,k)=vel*FLUX5(
     &           u(i,j,k-2,nrhs), u(i,j,k-1,nrhs), 
     &           u(i,j,k  ,nrhs), u(i,j,k+1,nrhs),
     &           u(i,j,k+2,nrhs), u(i,j,k+3,nrhs), vel)
          enddo
        enddo

        do i=IstrU,Iend
          if ( i.ge.imin .and. i.le.imax ) then
            vel=flux6(We(i-3,j,2),We(i-2,j,2),We(i-1,j,2),
     &                We(i  ,j,2),We(i+1,j,2),We(i+2,j,2),1.)
          else
            vel=0.5*(We(i-1,j,2)+We(i,j,2))
          endif
          FC(i,2)=vel*FLUX3(
     &         u(i,j,1,nrhs), u(i,j,2,nrhs), 
     &         u(i,j,3,nrhs), u(i,j,4,nrhs), vel)

          if ( i.ge.imin .and. i.le.imax ) then
            vel=flux6(We(i-3,j,N-2),We(i-2,j,N-2),We(i-1,j,N-2),
     &                We(i  ,j,N-2),We(i+1,j,N-2),We(i+2,j,N-2),1.)
          else   
            vel=0.5*(We(i-1,j,N-2)+We(i,j,N-2))
          endif
          FC(i,N-2)=vel*FLUX3(
     &         u(i,j,N-3,nrhs), u(i,j,N-2,nrhs), 
     &         u(i,j,N-1,nrhs), u(i,j,N  ,nrhs), vel)

          if ( i.ge.imin .and. i.le.imax ) then
            vel=flux6(We(i-3,j,1),We(i-2,j,1),We(i-1,j,1),
     &                We(i  ,j,1),We(i+1,j,1),We(i+2,j,1),1.)
          else
            vel=0.5*(We(i-1,j,1)+We(i,j,1))
          endif
          FC(i,1)=vel*FLUX2(
     &         u(i,j,1,nrhs), u(i,j,2,nrhs), vel, cdif)

          if ( i.ge.imin .and. i.le.imax ) then
            vel=flux6(We(i-3,j,N-1),We(i-2,j,N-1),We(i-1,j,N-1),
     &                We(i  ,j,N-1),We(i+1,j,N-1),We(i+2,j,N-1),1.)
          else
            vel=0.5*(We(i-1,j,N-1)+We(i,j,N-1))
          endif
          FC(i,N-1)=vel*FLUX2(
     &         u(i,j,N-1,nrhs), u(i,j,N,nrhs), vel, cdif)
	    
          FC(i,0)=0.
          FC(i,N)=0.
        enddo
