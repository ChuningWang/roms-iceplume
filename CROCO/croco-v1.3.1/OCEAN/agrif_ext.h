 
!
!-----------------------------------------------------------------------
! Load rhs values into additional AGRIF shared array for nesting
!-----------------------------------------------------------------------
!  

#ifdef AGRIF
      if (FIRST_2D_STEP) then
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            Zt_avg3(i,j,0)=zeta(i,j,kstp)       
          enddo
        enddo 
        do j=JstrR,JendR
          do i=Istr,IendR
          du_avg3(i,j,0)  = DUon(i,j)
          enddo
        enddo 
        do j=Jstr,JendR
          do i=IstrR,IendR
          dv_avg3(i,j,0)  = DVom(i,j)
          enddo
        enddo 
      endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   Zt_avg3(START_2D_ARRAY,0))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   du_avg3(START_2D_ARRAY,0))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   dv_avg3(START_2D_ARRAY,0))
# endif

# ifdef RVTK_DEBUG_ADVANCED
       if (.not.agrif_Root()) then
C$OMP BARRIER
C$OMP MASTER
       call check_tab2d(Zt_avg3(:,:,0),'Zt_avg3 (index 0) step2d','r')
       call check_tab2d(DU_avg3(:,:,0),'DU_avg3 (index 0) step2d','u')
       call check_tab2d(DV_avg3(:,:,0),'DV_avg3 (index 0) step2d','v')
C$OMP END MASTER  
       endif
# endif  
#endif /* AGRIF */        
