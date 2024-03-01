!
!----------------------------------------------------------
! Compute vertical advective fluxes 
! using 5th-order WENO scheme
!----------------------------------------------------------
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
!
!----------------------------------------------------------------------
!  k loop: FC
!----------------------------------------------------------------------
!
        do k=3,N-3
          do i=Istr,Iend
            if ( j.ge.jmin .and. j.le.jmax ) then
              vel=flux6(We(i,j-3,k),We(i,j-2,k),We(i,j-1,k),
     &                  We(i,j  ,k),We(i,j+1,k),We(i,j+2,k),1.)
            else
              vel=0.5*(We(i,j-1,k)+We(i,j,k))
            endif
            FC(i,k)=vel*FLUX5(
     &           v(i,j,k-2,nrhs), v(i,j,k-1,nrhs), 
     &           v(i,j,k  ,nrhs), v(i,j,k+1,nrhs),
     &           v(i,j,k+2,nrhs), v(i,j,k+3,nrhs), vel)
          enddo
        enddo

        do i=Istr,Iend
          if ( j.ge.jmin .and. j.le.jmax ) then
            vel=flux6(We(i,j-3,2),We(i,j-2,2),We(i,j-1,2),
     &                We(i,j  ,2),We(i,j+1,2),We(i,j+2,2),1.)
          else
            vel=0.5*(We(i,j-1,2)+We(i,j,2))
          endif
          FC(i,2)=vel*FLUX3(
     &         v(i,j,1,nrhs), v(i,j,2,nrhs), 
     &         v(i,j,3,nrhs), v(i,j,4,nrhs), vel)

          if ( j.ge.jmin .and. j.le.jmax ) then
            vel=flux6(We(i,j-3,N-2),We(i,j-2,N-2),We(i,j-1,N-2),
     &                We(i,j  ,N-2),We(i,j+1,N-2),We(i,j+2,N-2),1.)
          else   
            vel=0.5*(We(i,j-1,N-2)+We(i,j,N-2))
          endif
          FC(i,N-2)=vel*FLUX3(
     &         v(i,j,N-3,nrhs), v(i,j,N-2,nrhs), 
     &         v(i,j,N-1,nrhs), v(i,j,N  ,nrhs), vel)

          if ( j.ge.jmin .and. j.le.jmax ) then
            vel=flux6(We(i,j-3,1),We(i,j-2,1),We(i,j-1,1),
     &                We(i,j  ,1),We(i,j+1,1),We(i,j+2,1),1.)
          else
            vel=0.5*(We(i,j-1,1)+We(i,j,1))
          endif
          FC(i,1)=vel*FLUX2(
     &         v(i,j,1,nrhs), v(i,j,2,nrhs), vel, cdif)

          if ( j.ge.jmin .and. j.le.jmax ) then
            vel=flux6(We(i,j-3,N-1),We(i,j-2,N-1),We(i,j-1,N-1),
     &                We(i,j  ,N-1),We(i,j+1,N-1),We(i,j+2,N-1),1.)
          else
            vel=0.5*(We(i,j-1,N-1)+We(i,j,N-1))
          endif
          FC(i,N-1)=vel*FLUX2(
     &         v(i,j,N-1,nrhs), v(i,j,N,nrhs), vel, cdif)
	    
          FC(i,0)=0.
          FC(i,N)=0.
        enddo

