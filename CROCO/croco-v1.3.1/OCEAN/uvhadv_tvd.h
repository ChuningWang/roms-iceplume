        do j=Jstr,Jend
          do i=IstrU-1,Iend
            cff = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
            cu = cff*pn_u(i,j)*cdt/Hz(i,j,k) * pm(i,j)
            rrp= (u(i+2,j,k,nstp)-u(i+1,j,k,nstp)) 
#   ifdef MASKING
     &         *rmask(i+1,j)
#   endif             
            rr = (u(i+1,j,k,nstp)-u(i  ,j,k,nstp)) 
#   ifdef MASKING
     &         *rmask(i  ,j)
#   endif
            rrm= (u(i  ,j,k,nstp)-u(i-1,j,k,nstp))
#   ifdef MASKING
     &         *rmask(i-1,j)
#   endif
            cff1=( cff*(u(i+1,j,k,nstp)+u(i,j,k,nstp))
     &       -abs(cff)*rr )*0.5         
            Cr=limiter(cu,UFx(i,j),cff1,rrm,rr,rrp)
            UFx(i,j) = cff1*(1-Cr) + Cr*UFx(i,j)
          enddo
        enddo

        do j=JstrV-1,Jend
          do i=Istr,Iend
            cff = 0.5*(Hvom(i,j,k)+Hvom(i,j+1,k))
            cu = cff*pm_v(i,j)*cdt/Hz(i,j,k)*pn(i,j)
            rrp= (v(i,j+2,k,nstp)-v(i,j+1,k,nstp)) 
#   ifdef MASKING
     &         *rmask(i,j+1)
#   endif
            rr = (v(i,j+1,k,nstp)-v(i,j,k,nstp)) 
#   ifdef MASKING
     &         *rmask(i,j  )
#   endif
            rrm= (v(i,j  ,k,nstp)-v(i,j-1,k,nstp))
#   ifdef MASKING
     &         *rmask(i,j-1)
#   endif                    
            cff1=( cff*(v(i,j,k,nstp) + v(i,j+1,k,nstp))
     &        -abs(cff)*rr )*0.5
            Cr=limiter(cu,VFe(i,j),cff1,rrm,rr,rrp)
            VFe(i,j) = cff1*(1-Cr) + Cr*VFe(i,j)
          enddo
        enddo

        do j=Jstr,Jend+1
          do i=IstrU,Iend
            cff = 0.5*(Hvom(i,j,k)+Hvom(i-1,j,k))
            cu = cff*pm_u(i,j) /on_p(i,j)*cdt
     &    /( 0.25*( Hz(i  ,j,k)+Hz(i  ,j-1,k)
     &           +  Hz(i-1,j,k)+Hz(i-1,j-1,k) ) )
            rrp= (u(i,j+1,k,nstp)-u(i,j  ,k,nstp)) 
#   ifdef MASKING
     &         *pmask(i,j+1)
#   endif
            rr = (u(i,j  ,k,nstp)-u(i,j-1,k,nstp)) 
#   ifdef MASKING
     &         *pmask(i,j  )
#   endif
            rrm= (u(i,j-1,k,nstp)-u(i,j-2,k,nstp))
#   ifdef MASKING
     &         *pmask(i,j-1)
#   endif
            cff1=( cff*(u(i,j,k,nstp)+u(i,j-1,k,nstp))
     &        -abs(cff)*rr )*0.5      
            Cr=limiter(cu,UFe(i,j),cff1,rrm,rr,rrp)       
            UFe(i,j) =cff1*(1-Cr) + Cr*UFe(i,j)
          enddo
        enddo

        do j=JstrV,Jend
          do i=Istr,Iend+1
            cff = 0.5*(Huon(i,j,k)+Huon(i,j-1,k))
            cu = cff*pn_v(i,j)/om_p(i,j)*cdt
     &    /( 0.25*( Hz(i  ,j,k)+Hz(i  ,j-1,k)
     &           +  Hz(i-1,j,k)+Hz(i-1,j-1,k) ) )
            rrp=(v(i+1,j,k,nstp)-v(i  ,j,k,nstp))
#   ifdef MASKING
     &         *pmask(i+1,j)
#   endif
            rr =(v(i  ,j,k,nstp)-v(i-1,j,k,nstp))
#   ifdef MASKING
     &         *pmask(i  ,j)
#   endif
            rrm=(v(i-1,j,k,nstp)-v(i-2,j,k,nstp))
#   ifdef MASKING
     &         *pmask(i-1,j)
#   endif
            cff1=( cff*(v(i,j,k,nstp)+v(i-1,j,k,nstp))
     &       -abs(cff)*rr )*0.5
            Cr=limiter(cu,VFx(i,j),cff1,rrm,rr,rrp)           
            VFx(i,j) =cff1*(1-Cr) + Cr*VFx(i,j)
          enddo
        enddo
