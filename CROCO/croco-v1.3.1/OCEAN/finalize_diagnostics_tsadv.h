! $Id: finalize_diagnostics_tsadv ???? 2020-11-14 15:01:25Z jgula $
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
! This block is inside  itrc,k,j,i loops

              cff1= Hz(i,j,k) / (pm(i,j)*pn(i,j))
     &              * 0.5*(t(i,j,k,nstp,itrc)+t(i,j,k,nnew,itrc))

              TXadv(i,j,k,itrc)=TXadv(i,j,k,itrc)*cff1
              TYadv(i,j,k,itrc)=TYadv(i,j,k,itrc)*cff1
              TVadv(i,j,k,itrc)=TVadv(i,j,k,itrc)*cff1
              TForc(i,j,k,itrc)=TForc(i,j,k,itrc)*cff1
              Trate(i,j,k,itrc)=Trate(i,j,k,itrc)*cff1
              TVmix(i,j,k,itrc)=TVmix(i,j,k,itrc)*cff1
              THmix(i,j,k,itrc)=THmix(i,j,k,itrc)*cff1
!              TVmixt(i,j,k,itrc)=TVmixt(i,j,k,itrc)*cff

