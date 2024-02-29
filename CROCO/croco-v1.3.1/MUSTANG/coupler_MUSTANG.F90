#include "cppdefs.h"
!----------------------------------------------------------------------------
MODULE coupler_MUSTANG
!---------------------------------------------------------------------------- 

#if defined MUSTANG 

!&E==========================================================================
!&E                   ***  MODULE  coupler_MUSTANG  ***
!&E
!&E ** Purpose : concerns coupling MUSTANG with hydro code
!&E              
!&E ** Description :
!&E   subroutine coupl_conv2MUSTANG ! calculation of some variables needed 
!&E                                   by MUSTANG
!&E   subroutine coupl_MUSTANG2hydro ! transfert from MUSTANG to hydro code
!&E
!&E==========================================================================

#include "coupler_define_MUSTANG.h"

    USE comMUSTANG
    USE comsubstance
    USE module_MUSTANG
    USE module_substance

    IMPLICIT NONE

    ! functions & routines of this module, called outside :
    PUBLIC coupl_conv2MUSTANG, coupl_MUSTANG2hydro

    PRIVATE

    CONTAINS
    !!=======================================================================
    SUBROUTINE coupl_conv2MUSTANG(ifirst, ilast, jfirst, jlast, iappel,     &
                                  BATHY_H0, ssh, WATER_CONCENTRATION)
    !&E----------------------------------------------------------------------
    !&E                 ***  ROUTINE coupl_conv2MUSTANG  ***
    !&E
    !&E ** Purpose : transfer cv_Wat and htot et alt_cw1 computation
    !&E               only inside the domain, not at boundaries meshes
    !&E
    !&E ** Description :  
    !&E  arguments IN : BATHY_H0, ssh, WATER_CONCENTRATION
    !&E  arguments OUT: no (all variables  in comMUSTANG)
    !&E     
    !&E   variables OUT :   
    !&E       htot (total water height )
    !&E       epn_bottom (thickness of water the bottom layer)
    !&E       sal_bottom_MUSTANG,temp_bottom_MUSTANG (salinity, temperature 
    !&E           in the water bottom layer)
    !&E       cw_bottom_MUSTANG ( concentrations in the water bottom layer)
    !&E       roswat_bot ( water density  in the water bottom layer)
    !&E     
    !&E   initial call  (iappel=0 for initialization):  
    !&E       extraction of thickness, salinity, temperature and water 
    !&E           concentrations and densities in the bottom of the water 
    !&E           column
    !&E
    !&E   first call (iappel=1) :  
    !&E       extraction of thickness, salinity, temperature and water 
    !&E           concentrations and densities in the bottom of the water
    !&E           column
    !&E       calculation of total water height    
    !&E       calculation of  alt_cw1  : altitude of the computation point 
    !&E           of Cw in the  bottom layer
    !&E     
    !&E   second call (iappel=2) :  
    !&E       extraction of thickness, salinity, temperature and water 
    !&E           concentrations and densities in the bottom of the water 
    !&E           column
    !&E       calculation of total water height    
    !&E       if not MARS : conversion to transmit to MUSTANG the hydro 
    !&E           variables: SETTL_FLUXSUM_w2s: effective deposit flux  
    !&E           of the particle variables during transport
    !&E    
    !&E ** Called by :  MUSTANG_init (iappel=0)
    !&E                 sed_MUSTANG_update (iappel=1)
    !&E                 sed_MUSTANG_deposition (iappel=2)
    !&E
    !&E----------------------------------------------------------------------
   !! * Modules used
#include "scalars_F90.h"

   !! * Arguments 
   INTEGER, INTENT(IN)  :: ifirst, ilast, jfirst, jlast, iappel
   REAL(KIND=rsh),DIMENSION(ARRAY_BATHY_H0),INTENT(IN)  :: BATHY_H0                         
   REAL(KIND=rsh),DIMENSION(ARRAY_WATER_ELEVATION),INTENT(IN):: ssh                         
   REAL(KIND=rsh),DIMENSION(ARRAY_WATER_CONC), INTENT(IN) :: WATER_CONCENTRATION   
   !! * Local declarations
   INTEGER  :: iv, i, j

   !! * Executable part
   DO j=jfirst,jlast
      DO i=ifirst,ilast
      ! WARNING : not need to calculate at boundaries meshes where MUSTANG is not applied

           ! extraction of  concentrations in the bottom of the water column
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CROCO vecteur au temps 1, 2 ou 3 ????
# ifdef SALINITY
            sal_bottom_MUSTANG(i,j)=WATER_CONCENTRATION(i,j,1,1,itemp+1)
# else
            sal_bottom_MUSTANG(i,j)=35.
# endif
# ifdef TEMPERATURE
            temp_bottom_MUSTANG(i,j)=WATER_CONCENTRATION(i,j,1,1,itemp)
# else
            temp_bottom_MUSTANG(i,j)=15.
# endif
            cw_bottom_MUSTANG(1:nv_adv,i,j)=WATER_CONCENTRATION(i,j,1,1,itsubs1:itsubs2)

            ! thickness of the bottom water layer or altitude at the top of the bottom layer
            ! + water density in the bottom water layer
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! CROCO z_w connu via scalars.h
               epn_bottom_MUSTANG(i,j)=z_w(i,j,1)-z_w(i,j,0)
               roswat_bot(i,j)= rho(i,j,1)+rho0  

       ENDDO
   ENDDO

   IF(iappel > 0 ) THEN
   ! WARNING : need to calculate htot at all meshes (imin+1:imax, jmin+1: jmax)
   ! at interior meshes : for ljmin-1 and ljmax+1 , BATHY_H0 and ssh are known if MPI exchange was done after change
    DO j=jfirst-1,jlast+1
      DO i=ifirst-1,ilast+1
          !  htot : total water height
          htot(i,j)=z_w(i,j,N)+h(i,j)
       ENDDO
    ENDDO

   ENDIF
   IF (iappel == 1) THEN
   ! first call before evaluation of settling velocities, erosion, consolidation, diffusion
   
     DO j=jfirst,jlast
       DO i=ifirst,ilast
          ! alt_cw1 : altitude of the computation point of Cw in the  bottom layer
          ! could be eliminated if the concentration calculation point is always in the middle of the layer
          ! but if not true, you have to calculate this altitude differently
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          alt_cw1(i,j) = z_r(i,j,1)-z_w(i,j,0) 
       ENDDO
     ENDDO

    ELSE IF (iappel==2) THEN

   ! second call before evaluation of sediment deposition :  cumulated settling flux 

     DO j=jfirst,jlast
      DO i=ifirst,ilast
#if defined key_noTSdiss_insed
         DO iv=1,nv_adv
#else
         DO iv=-1,nv_adv
#endif
                flx_w2s_sum(iv,i,j)=SETTL_FLUXSUM_w2s(i,j,IV_HOSTMODEL) ! flux de depot cumule effectif apres transport
         ENDDO
     ENDDO
     ENDDO

    ENDIF

  END SUBROUTINE coupl_conv2MUSTANG      
!!===========================================================================
    SUBROUTINE coupl_MUSTANG2hydro(ifirst, ilast, jfirst, jlast)                                         
    !&E----------------------------------------------------------------------
    !&E                 ***  ROUTINE coupl_MUSTANG2flx ***
    !&E
    !&E ** Purpose : transfer flx_w2s , ws3 ,flx_s2w for hydro code 
    !&E
    !&E ** Description : conversion for hydro code 
    !&E  arguments OUT: no because stored in comMUSTANG
    !&E     SETTL_FLUX_w2s : deposit trends 
    !&E     EROS_FLUX_s2w : erosion flux 
    !&E     EROS_FLUX_TEMP_s2w et eros_flix_SAL : erosion flux for 
    !&E                                           temperature, salinity
    !&E     
    !&E ** Called by :  sed_MUSTANG_update
    !&E
    !&E----------------------------------------------------------------------
    !! * Arguments 
    INTEGER, INTENT(IN) :: ifirst, ilast, jfirst, jlast                     
    !! * Local declarations
    INTEGER :: iv, i, j

    !! * Executable part
    ! exchange erosion  fluxes
    DO j = jfirst, jlast
        DO i = ifirst, ilast
            DO iv = 1, nvp
                SETTL_FLUX_w2s(i, j, IV_HOSTMODEL) = flx_w2s(iv, i, j)
            ENDDO
            DO iv = 1, nv_adv
                EROS_FLUX_s2w(i, j, IV_HOSTMODEL) = flx_s2w(iv, i, j)
            ENDDO
            ! temperature
            EROS_FLUX_s2w(i, j, ITEMP_HOSTMODEL) = flx_s2w(-1, i, j)
            ! salinity
            EROS_FLUX_s2w(i, j, ISAL_HOSTMODEL) = flx_s2w(0, i, j)
        ! no transfer of SETTL_FLUX_w2s_TEMP et SAL and for dissolved subst. 
        ! because they are merged in EROS_FLUX_s2w for dissolved variables
        ! (EROS_FLUX_s2w=erosion-settling+consolidation-diffusion) 
        ENDDO
    ENDDO

END SUBROUTINE coupl_MUSTANG2hydro    
!!===========================================================================
#endif /* ifdef MUSTANG */

END MODULE coupler_MUSTANG
