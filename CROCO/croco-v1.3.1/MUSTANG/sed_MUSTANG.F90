!---------------------------------------------------------------------------
MODULE sed_MUSTANG
!---------------------------------------------------------------------------
   
#include "cppdefs.h"

#if defined MUSTANG 

   !&E==========================================================================
   !&E                   ***  MODULE  sed_MUSTANG  ***
   !&E
   !&E ** Purpose : concerns subroutines related to sediment dynamics 
   !&E 
   !&E ** Description :
   !&E     subroutine MUSTANG_update      ! updates forcings and remoddeling sediment  
   !&E     subroutine MUSTANG_deposition     ! sedimentation  sediment
   !&E     subroutine MUSTANG_morpho         ! update bathymetry with sediment thickness (morphodynamic)
   !&E 
   !&E     subroutine sed_MUSTANG_comp_z0sed     ! computation of sediment roughness length (based on mean diameter)
   !&E     subroutine sed_MUSTANG_comp_z0hydro   ! computation of hydrodynamic roughness length (based on mean diameter)
   !&E     subroutine sed_MUSTANG_depflx         ! computation of deposition tendancy
   !&E     subroutine sed_MUSTANG_sandconcextrap ! extrapolation of sand concentration
   !&E     subroutine sed_MUSTANG_erosion        ! computes erosion thickness
   !&E     subroutine sed_MUSTANG_comp_tocr_mixsed    ! critical shear stress for erosion (Pa)
   !&E     subroutine sed_MUSTANG_comp_eros_flx ! comp_eros_flx : erosion flux (kg/m2/s)
   !&E     subroutine sed_MUSTANG_effdep         ! effective deposition
   !&E     subroutine sed_MUSTANG_fusion         ! fusion of two first bottom layers
   !&E     subroutine sed_MUSTANG_slipdepo       ! slip lateral deposit
   !&E     subroutine sed_MUSTANG_split_surflayer ! splitting surface layer when too thick si key_splitlayersurf
   !&E     subroutine sed_MUSTANG_consol_diff_bioturb   
   !&E     subroutine sed_MUSTANG_coefbioturb_part        ! estimating bioturbation  coefficients in the sediment
   !&E     subroutine sed_MUSTANG_coefbioturb_diss        ! estimating  biodiffusion coefficients in the sediment
   !&E     subroutine sed_MUSTANG_constitutivrel  !
   !&E     subroutine sed_MUSTANG_outres        ! preparation of the tables of results to be written according to the choice of the user
   !&E     subroutine sed_MUSTANG_interpout_dzs ! interpolation of  output layer thickness 
   !&E     subroutine sed_MUSTANG_interpout_cvs ! interpolation of sediment concentrations in output layers
   !&E
   !&E  ---for new version MUSTANG V2 ----
   !&E     subroutine MUSTANGV2_manage_active_layer  ! thickness et composition of active layer 
   !&E     subroutine MUSTANGV2_fusion_with_poro    ! fusion of two surficial sediment layers + poro 
   !&E     subroutine MUSTANGV2_comp_eros_flx_indep  ! erosion flux but independant between sediment classes
   !&E     subroutine MUSTANGV2_borne_and_apply_erosion_tot
   !&E     subroutine MUSTANGV2_manage_small_mass_in_ksmax
   !&E     subroutine MUSTANGV2_comp_poro_mixsed  ! porosity estimation in case of non-cohesive/cohesive sediment mixtures
   !&E     subroutine MUSTANGV2_eval_bedload            ! bedload transport
   !&E     subroutine MUSTANGV2_eval_dissvar_IWSflux            ! dissolved variables in interstitial waters and fluxes at interface water-sed
   !&E     FUNCTION isitcohesive                ! function which return TRUE if sediment is cohesif (.F. if not)
   !&E     FUNCTION dzsminvar                   ! function which return dzsmin
   !&E     FUNCTION MUSTANG_E0sand              ! function to compute E0 sand

   !&E  **********************************************************************
   !&E  *** requires to have advected substances in a specific order       ***
   !&E  ***  in water : cw (iv,k) k=1,kmax=NB_LAYER_WAT bottom to surface  ***
   !&E  ***  in sediment : cs (iv,k) k=ksdmin,ksdmax                       ***
   !&E  ***                                                                ***
   !&E  ***       First nvp particulate substances                         ***
   !&E  ***    first GRAVEL (indices igrav1 to igrav2) : coarser to finer  ***
   !&E  ***    then  SAND (isand1=igrav1+1 to isand2) : coarser to finer   ***
   !&E  ***    then MUDS (imud1=isand2+1 to imud2                          ***
   !&E  ***                                                                ***
   !&E  ***    then nvd Dissolved substances                               *** 
   !&E  ***                                                                ***
   !&E  **********************************************************************
   !&E
   !&E     Used with specific date issued for hydrodynamic module : conversion in coupleur_MUSTANG.F90
   !&E       kmax : number of layers in water column
   !&E       h0fond : residual water thickness (in m). This thickness is added to the computed water level (xe)
   !&E                to avoid possible mass loss during drying (due to the fact that the scheme for water continuity is not positive)
   !&E       RHOREF : reference sea water density (in kg/m3)
   !&E       dt : time step 
   !&E
   !&E     Used with specific date issued for substances module : 
   !&E       nv_adv : number of advected (transported) state variables
   !&E       nv_tot : total number of  variables (state+intermediate+driving..)
   !&E       nv_mud : total number of  variables type MUDS
   !&E       nvp :  number of particulate advected variables
   !&E       nvpc :  number of constitutive particulate advected variables
   !&E       igrav1,igrav2 : first and last indices of gravel
   !&E       isand1,isand2 : first and last indices of sand (first : corsaer sand - last : finer sand)
   !&E       imud1,imud2 :first and last indices of mud
   !&E       typart(1:nv_state) : =0 for non-constitutive part. variable, typart=1 for constitutive particulate variable
   !&E       ws_free_opt,ws_free_para,ws_free_min,ws_free_max,ws_hind_opt,ws_hind_para : parameters defininf settling velocities
   !&E       water concentrations [cw(iv,k)] iv=1:nv_tot & k=1:kmax  (Unity : Masse/m3 of water)
   !&E       sediment concentrations [cs(iv,k)] iv=-1:nv_tot & k=ksdmin:ksdmax 
   !&E                                iv=-1 : temperature  - iv=0 : salinity 
   !&E                                iv=1,nvpc => Unity : kg/m3 of total dry sediment      
   !&E                                iv=nvpc+1,nvp => Unity : Masse/m3 of total dry sediment     
   !&E                                iv=nvp+1,nv_state => Unity : Masse/m3 of interstitial water     
   !&E       tocd : critical stress of deposition
   !&E       irkm_var_assoc : number of the constitutive particulate variable associated with a non constitutive particulate variable (both have same settling velocity)
   !&E                        If = 0, particulate variable with its own settling velocity (or dissolved variable)
   !&E
   !&E
   !&E===================================================================================================================
#include "coupler_define_MUSTANG.h"

   !! variables  SUBSTANCE and variable from croco known via 
   USE comsubstance
   USE module_substance

   USE comMUSTANG 
   USE coupler_MUSTANG 

   IMPLICIT NONE
   
   !! * Accessibility
   ! functions & routines of this module, called outside :
   PUBLIC MUSTANG_update, MUSTANG_deposition, sed_MUSTANG_outres
   PUBLIC sed_MUSTANG_comp_z0hydro, MUSTANG_E0sand
#if defined MORPHODYN
   PUBLIC MUSTANG_morpho
#endif
#if defined key_MUSTANG_V2
   PUBLIC MUSTANGV2_comp_poro_mixsed
#endif
#ifdef key_MUSTANG_splitlayersurf
   PUBLIC sed_MUSTANG_split_surflayer
#endif

   PRIVATE
   
 CONTAINS
 
!!==============================================================================
 
  SUBROUTINE MUSTANG_update(ifirst, ilast, jfirst, jlast,     &
               WATER_CONCENTRATION, z0hydro,                  &
               WATER_ELEVATION,                               &
#if defined key_MUSTANG_lateralerosion || defined key_MUSTANG_bedload
               BAROTROP_VELOCITY_U, BAROTROP_VELOCITY_V,             &
#endif
               saliref_lin, temperef_lin, dt_true)

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE MUSTANG_update ***
   !&E
   !&E ** Purpose : update forcings (bottom stress, u* , settling velocities..)
   !&E               + remodelling (erosion, consolidation, diffusion, bioturbation
   !&E               (only inside the domain, not at boundaries)
   !&E
   !&E  arguments IN : 
   !&E         loops  :ifirst,ilast,jfirst,jlast
   !&E         parametres ref  :RHOREF, saliref_lin,temperef_lin
   !&E         time  :dt_true,t (DOUBLE PRECISION)
   !&E         hydro  :WATER_ELEVATION,BAROTROP_VELOCITY_U,BAROTROP_VELOCITY_V
   !&E         concentrations  : WATER_CONCENTRATION,SALINITY_MOD,TEMPERATURE_MOD
   !&E         [settling velocities (transmitted as argument or by USE as in MARS or in CROCO)]
   !&E
   !&E  arguments OUT:
   !&E         [settling velocities ]
   !&E         and the sediment has been remodelled
   !&E          cv_sed, dzs, ksma, poro have changed
   !&E                
   !&E ** Called by :  step before mouvement and transport equations
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
    USE sed_MUSTANG_HOST,    ONLY :  sed_MUSTANG_settlveloc
    USE sed_MUSTANG_HOST,    ONLY :  sed_skinstress
    USE sed_MUSTANG_HOST,    ONLY :  sed_gradvit
#ifdef key_MUSTANG_bedload
    USE sed_MUSTANG_HOST,    ONLY :  sed_bottom_slope
#if defined MPI 
      USE sed_MUSTANG_HOST,    ONLY :  sed_exchange_flxbedload
      USE sed_MUSTANG_HOST,    ONLY :  sed_exchange_maskbedload
#endif
#endif
#if defined MPI  && defined key_MUSTANG_lateralerosion
    USE sed_MUSTANG_HOST,    ONLY :  sed_exchange_s2w
#endif
#if defined MUSTANG_CORFLUX
    USE sed_MUSTANG_HOST,    ONLY :  sed_obc_corflu
    USE sed_MUSTANG_HOST,    ONLY :  sed_meshedges_corflu
#if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      USE sed_MUSTANG_HOST,    ONLY :  sed_exchange_corflu
#endif
#endif
#if defined key_BLOOM_insed && defined key_oxygen && ! defined key_biolo_opt2
   USE reactionsinsed,  ONLY : reactions_in_sed
   USE bioloinit,     ONLY : p_txfiltbenthmax
#endif
#if defined key_MUSTANG_flocmod
   USE flocmod,  ONLY : flocmod_main
#endif

   !! * Arguments
   INTEGER, INTENT(IN)                                       :: ifirst, ilast, jfirst, jlast                           
   REAL(KIND=rsh),INTENT(IN)                                 :: saliref_lin, temperef_lin 
   REAL(KIND=rlg),INTENT(IN)                                 :: dt_true  ! !  (dt_true=halfdt in MARS)
   REAL(KIND=rsh),DIMENSION(ARRAY_Z0HYDRO),INTENT(INOUT)          :: z0hydro                         
   REAL(KIND=rsh),DIMENSION(ARRAY_WATER_ELEVATION),INTENT(INOUT)  :: WATER_ELEVATION                         
#if defined key_MUSTANG_lateralerosion || defined key_MUSTANG_bedload                        
   REAL(KIND=rsh),DIMENSION(ARRAY_VELOCITY_U),INTENT(IN)          :: BAROTROP_VELOCITY_U                        
   REAL(KIND=rsh),DIMENSION(ARRAY_VELOCITY_V),INTENT(IN)          :: BAROTROP_VELOCITY_V   
#endif                      
#if defined key_MUSTANG_flocmod || defined key_BLOOM_insed
   REAL(KIND=rsh),DIMENSION(ARRAY_WATER_CONC), INTENT(INOUT) :: WATER_CONCENTRATION         
#else
   REAL(KIND=rsh),DIMENSION(ARRAY_WATER_CONC), INTENT(IN) :: WATER_CONCENTRATION         
#endif 

   !! * Local declarations
   INTEGER        ::  i,j,k,iv,ivp,ksmin,ksmax,iappel
   REAL(KIND=rsh) ::  dtinv
   REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY) :: workexch
   !!---------------------------------------------------------------------------
   !! * Executable part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!     COUPLING Z0 hydro : update z0hydro  if l_z0_coupl      !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(l_z0hydro_coupl)CALL sed_MUSTANG_comp_z0hydro(ifirst,ilast,jfirst,jlast,z0hydro)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! PRELIMINARY EVALUATIONS FOR COUPLING !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! coupler for some calculations specific to each model ==> htot, alt_cw1, epn_bottom
    ! + computation of bottom concentrations..
    iappel=1
    CALL coupl_conv2MUSTANG(ifirst,ilast,jfirst,jlast,iappel,BATHY_H0,WATER_ELEVATION,   &
                       WATER_CONCENTRATION)
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! FORCING :    bottom shear stress (tauskin)!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Update sediment roughness length
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (.not. l_z0seduni) then
      call sed_MUSTANG_comp_z0sed(ifirst, ilast, jfirst, jlast, BATHY_H0)
    endif

    call sed_skinstress(ifirst, ilast, jfirst, jlast)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! FORCING : G dU/dz (taux de cisaillement : shear rate)   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! computing  gradvit - depending on hydro code ==> in MUSTANG_HYDROCODE.F90
      ! out : gradvit in common comMUSTANG
    call sed_gradvit(ifirst, ilast, jfirst, jlast)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! vitesse de chute : settling rate     !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call sed_MUSTANG_settlveloc(ifirst, ilast, jfirst, jlast,   &
                           WATER_CONCENTRATION)

#ifdef key_MUSTANG_flocmod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! FLOCMOD :    compute aggregation /fragmentation processes  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j=jfirst,jlast
        DO i=ifirst,ilast   
            IF(htot(i,j) > h0fond) THEN
                DO k=1,NB_LAYER_WAT
                    CALL flocmod_main( dt_true, &
                        t(i,j,k,nstp,itsubs1-1+imud1:itsubs1-1+nvpc),  &
                        gradvit(k,i,j) )
                ENDDO
            ENDIF
        ENDDO
    ENDDO

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! tendance au depot : deposit tendency   flx_w2s (m.s-1) !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
    CALL sed_MUSTANG_depflx(ifirst, ilast, jfirst, jlast)
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Correction for sand horizontal transport  corflux,corfluy ,flx_w2s(sand)  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL sed_MUSTANG_sandconcextrap(ifirst, ilast, jfirst, jlast)

   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! interpolation of corflux,corfluy at mesh edges  and                               !!
!!!!!! treatment of corflux,corfluy at grid corners and exchange between MPI processors  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined MUSTANG_CORFLUX
    IF (isand2.GE.1) THEN
        CALL sed_obc_corflu(ifirst, ilast, jfirst, jlast)
#if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
        CALL sed_exchange_corflu(ifirst, ilast, jfirst, jlast, 0)
#endif
        ! corflux are interpolated on mesh edges (in u & v) 
        ! depends on model mesh (ARAKAWA grid), coded in sed_MUSTANG_HOST
        CALL sed_meshedges_corflu(ifirst, ilast, jfirst, jlast)
        CALL sed_obc_corflu(ifirst, ilast, jfirst, jlast)
#if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
        CALL sed_exchange_corflu(ifirst, ilast, jfirst, jlast, 1)
#endif
        !! for substances which are sorbed or associated with sand variables
         DO ivp = nvpc+1, nvp
            IF(irkm_var_assoc(ivp) .NE. 0 .AND. irkm_var_assoc(ivp) .LT. imud1 ) THEN
                corflux(ivp, :, :) = corflux(irkm_var_assoc(ivp), :, :)
                corfluy(ivp, :, :) = corfluy(irkm_var_assoc(ivp), :, :)
            ENDIF
         ENDDO             
    ENDIF
#endif

   dtinv=1.0_rsh/REAL(dt_true,rsh)
  
#if ! defined key_noTSdiss_insed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TEMERATURE in SEDIMENT                             !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   CALL sed_MUSTANG_Temperatur_in_sed(ifirst, ilast, jfirst, jlast,  &
                                                dt_true, dtinv)
#endif
                           
#if defined key_BLOOM_insed && defined key_oxygen && ! defined key_biolo_opt2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BIO PROCESSES in SEDIMENT                                        !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL reactions_in_sed(ifirst, ilast, jfirst, jlast, BATHY_H0, dt_true, dtinv)

    !**TODO** : create sed_exchange_cvwat in sed_MUSTANG_CROCO
    !IF(p_txfiltbenthmax .NE. 0.0_rsh) THEN
    ! CALL sed_exchange_cvwat_MARS(WATER_CONCENTRATION)
    !ENDIF

#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CONSOLIDATION & DIFFUSION   !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   IF(l_dyn_insed)CALL sed_MUSTANG_consol_diff_bioturb(ifirst, ilast, jfirst, jlast,  &
                saliref_lin, temperef_lin, dt_true)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          EROSION             !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined key_MUSTANG_lateralerosion
! initialization LATERAL EROSION 
   IF(coef_erolat .NE. 0.0_rsh) THEN
      flx_s2w_corip1(:,:,:) = 0.0_rsh
      flx_s2w_corim1(:,:,:) = 0.0_rsh
      flx_s2w_corjp1(:,:,:) = 0.0_rsh
      flx_s2w_corjm1(:,:,:) = 0.0_rsh
#if ! defined key_nofluxwat_IWS
      phieau_s2w_corip1(:,:) = 0.0_rsh
      phieau_s2w_corim1(:,:) = 0.0_rsh
      phieau_s2w_corjp1(:,:) = 0.0_rsh
      phieau_s2w_corjm1(:,:) = 0.0_rsh
#endif
   ENDIF  ! end if coef_erolat
#endif
   
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
! initialization BEDLOAD fluxes and masks 
   flx_bx(:,:,:)=0.0_rsh
   flx_by(:,:,:)=0.0_rsh

   sedimask_h0plusxe(:,:) = 0.0_rsh
   DO j = jfirst-1, jlast+1
     DO i = ifirst-1, ilast+1
         IF (htot(i,j) .GT. hmin_bedload) THEN
           sedimask_h0plusxe(i,j) = 1.0_rsh
         ENDIF
     ENDDO
   ENDDO
#if defined MPI
   CALL sed_exchange_maskbedload(ifirst, ilast, jfirst, jlast)
#endif

#ifdef key_MUSTANG_specif_outputs
  ! flx_bx and by
  varspecif3Dnv_save(5:6,:,:,:)=0.0_rsh
#endif

#if defined MORPHODYN  
     IF (l_slope_effect_bedload .AND. it_morphoYes==1 ) THEN
        CALL sed_bottom_slope(ifirst, ilast, jfirst, jlast, BATHY_H0)
        it_morphoYes = 0
     ENDIF
#endif   

#endif  /*end key_MUSTANG_bedload (version V2)*/

   CALL sed_MUSTANG_erosion(ifirst, ilast, jfirst, jlast, dtinv,     &
#if defined key_MUSTANG_lateralerosion || defined key_MUSTANG_bedload
                           BAROTROP_VELOCITY_U, BAROTROP_VELOCITY_V,     &
#endif
                             dt_true)

#if defined MPI && defined key_MUSTANG_bedload
       if (float(ifirst+ii*Lm) .EQ. IMIN_GRID) then
        flx_bx(:,ifirst-1,:)=flx_bx(:,ifirst,:)
        flx_by(:,ifirst-1,:)=flx_by(:,ifirst,:)
       endif
# if (!defined DUNE    || (defined DUNE    && defined DUNE3D))
       if (float(jfirst+jj*Mm) .EQ. JMIN_GRID) then
        flx_bx(:,:,jfirst-1)=flx_bx(:,:,jfirst)
        flx_by(:,:,jfirst-1)=flx_by(:,:,jfirst)
       endif
# endif
#endif
#if (!defined MPI && defined key_MUSTANG_bedload)
        flx_bx(:,ifirst-1,:)=flx_bx(:,ifirst,:)
        flx_by(:,ifirst-1,:)=flx_by(:,ifirst,:)
# if (!defined DUNE    || (defined DUNE    && defined DUNE3D))
        flx_bx(:,:,jfirst-1)=flx_bx(:,:,jfirst)
        flx_by(:,:,jfirst-1)=flx_by(:,:,jfirst)
# endif
#endif

#if defined key_MUSTANG_bedload && defined MPI 
    call sed_exchange_flxbedload(ifirst, ilast, jfirst, jlast)
#endif                             
                           

#if defined key_MUSTANG_lateralerosion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Fluxes correction if LATERAL EROSION !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   IF(coef_erolat .NE. 0.0_rsh) THEN
     ! lateral erosion of dry cell
#if defined MPI
    call sed_exchange_s2w(ifirst, ilast, jfirst, jlast)
#endif

      ! correction : neighboring cells of eroded laterally  dry cell receive one fraction of eroded sediment 
      DO j=jfirst,jlast
        DO i=ifirst,ilast
      ! warning it may be different for complex grid
             DO iv=-1,nv_adv
               flx_s2w(iv,i,j)=flx_s2w(iv,i,j)+ dtinv*(                            &
                       +flx_s2w_corip1(iv,i-1,j)+flx_s2w_corim1(iv,i+1,j)   &
                       +flx_s2w_corjp1(iv,i,j-1)+flx_s2w_corjm1(iv,i,j+1))
             ENDDO
#if ! defined key_nofluxwat_IWS
             phieau_s2w(i,j)=phieau_s2w(i,j)+phieau_s2w_corip1(i-1,j)+  &
                    phieau_s2w_corim1(i+1,j)+phieau_s2w_corjp1(i,j-1)+phieau_s2w_corjm1(i,j+1)
#endif
        END DO
      END DO
                          
   ENDIF  ! end if coef_erolat
#endif

#if defined key_MUSTANG_specif_outputs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! save cumulated erosion Fluxes  of constitutive particulate variables !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef key_MUSTANG_V2
   IF(l_outsed_flx_WS_all .OR. l_outsed_flx_WS_int) THEN
#else
   IF(l_outsed_flx_WS_all) THEN
#endif
      DO j=jfirst,jlast
        DO i=ifirst,ilast
          DO iv=1,nvpc
             varspecif3Dnv_save(2,iv,i,j)=varspecif3Dnv_save(2,iv,i,j)+  &
                                           flx_s2w(iv,i,j)*REAL(dt_true,rsh) 
             ! cumulate flx_s2w (in kg/m2 integrated on time step dt_true as for deposit flux )
          ENDDO
#ifdef key_MUSTANG_V2
          DO iv=isand1,isand2
             !flx_s2w_noncoh
             varspecif2D_save(11,i,j)=varspecif2D_save(11,i,j)+  &
                                        flx_s2w(iv,i,j)*REAL(dt_true,rsh)
          END DO
          DO iv=imud1,imud2
             !flx_s2w_coh
            varspecif2D_save(13,i,j)=varspecif2D_save(13,i,j)+   &
                                        flx_s2w(iv,i,j)*REAL(dt_true,rsh)
          END DO
#endif
       END DO
     END DO
   ENDIF
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   conversion of deposit flux for  hydro code                                                        !!!!!
! + conversion of erosion flux                                                                         !!!!
! + water flux at interface sediment/water resulting from deposit at previous time step and of erosion  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL coupl_MUSTANG2hydro(ifirst,ilast,jfirst,jlast)

#if defined key_nofluxwat_IWS || defined key_noTSdiss_insed
!!  not taking into account water fluxes threw sediment-water interface
#else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    taking into account the water flux from sediment to water in the bottom layer                   !!!!!!
!    due to erosion and consolidation                                                                !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! WATER_FLUX_INPUT_BOTCELL = WATER_FLUX_INPUTS (k,i,j) in bottom cell = phieau(1,:,:)
    WATER_FLUX_INPUT_BOTCELL=WATER_FLUX_INPUT_BOTCELL+phieau_s2w(:,:)/dt_true
    phieau_s2w(:,:)=0.0_rlg
#endif

  END SUBROUTINE MUSTANG_update
  
  !!==============================================================================
  
  SUBROUTINE MUSTANG_deposition(ifirst, ilast, jfirst, jlast,  &
                        WATER_ELEVATION,                       &
                        WATER_CONCENTRATION)
              
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE MUSTANG_deposition  ***
   !&E
   !&E ** Purpose : settling from water to sediment, sliding processes 
   !&E                    and sediment sedimentation 
   !&E               (only inside the domain, not at boundaries)
   !&E
   !&E         for MARS : ifirst=imin+2, ilast=imax-1, jfirst=jmin+2,  jlast=jmax-1
   !&E                    for interior processors : ifirst=limin, ilast=limax, jfirst=ljmin,  jlast=ljmax
   !&E
   !&E ** Description : 
   !&E  arguments IN : 
   !&E         loops  :ifirst,ilast,jfirst,jlast
   !&E         hydro  :WATER_ELEVATION
   !&E         concentrations  : WATER_CONCENTRATION,SALINITY_MOD,TEMPERATURE_MOD
   !&E
   !&E  arguments OUT:
   !&E
   !&E         the sediment has been remodelled
   !&E          cv_sed, dzs, ksma, poro have changed
   !&E         
   !&E
   !&E ** Called by :  step after transport equations
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
#if defined MPI  && defined key_MUSTANG_slipdeposit
    USE sed_MUSTANG_HOST,    ONLY :  sed_exchange_w2s
#endif
   !! * Arguments
   INTEGER, INTENT(IN)  :: ifirst, ilast, jfirst, jlast 
   REAL(KIND=rsh),DIMENSION(ARRAY_WATER_ELEVATION),INTENT(INOUT) :: WATER_ELEVATION
#if defined key_BLOOM_insed
   REAL(KIND=rsh),DIMENSION(ARRAY_WATER_CONC), INTENT(INOUT)  :: WATER_CONCENTRATION   
#else
   REAL(KIND=rsh),DIMENSION(ARRAY_WATER_CONC), INTENT(IN)  :: WATER_CONCENTRATION   
#endif


   !! * Local declarations
   INTEGER        :: iappel,iexchge_MPI_cvwat
   !! * Executable part 
  
   !deposit if strong bottom slope   
     flx_w2s_corim1(:,:,:)=0.0_rsh   
     flx_w2s_corip1(:,:,:)=0.0_rsh
     flx_w2s_corjm1(:,:,:)=0.0_rsh
     flx_w2s_corjp1(:,:,:)=0.0_rsh
     flx_w2s_corin(:,:,:)=0.0_rsh 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! coupler for some calculations specific to each model                  !!!!!!
!!!!! ==> htot, alt_cw1, epn_bottom and new bottom concentrations           !!!!!!
!!!!! and conversion of effectif deposit flux from hydro code to MUSTANG    !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    iappel=2
    CALL coupl_conv2MUSTANG(ifirst,ilast,jfirst,jlast,iappel,    &
                            BATHY_H0,WATER_ELEVATION,            &
                            WATER_CONCENTRATION )        
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if ! defined key_noTSdiss_insed
   flx_w2s(nvp+1:nv_adv,:,:)=0.0_rsh
   flx_w2s(-1:0,:,:)=0.0_rsh
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! deposit slip if steep slope (slidepo)        !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef key_MUSTANG_slipdeposit
   IF(slopefac .NE. 0.0_rsh) THEN
     CALL sed_MUSTANG_slipdepo(ifirst, ilast, jfirst, jlast)
#if defined MPI
    CALL sed_exchange_w2s(ifirst, ilast, jfirst, jlast)
#endif
   ENDIF 
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! effective deposit after variables transport and settling =>  sedimentation   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call sed_MUSTANG_effdep(ifirst, ilast, jfirst, jlast, iexchge_MPI_cvwat)

    !**TODO** code sed_exchange_cvwat CALL sed_exchange_cvwat_MARS(WATER_CONCENTRATION,iexchge_MPI_cvwat)

  END SUBROUTINE MUSTANG_deposition
 
!!===========================================================================
#if defined MORPHODYN
    SUBROUTINE MUSTANG_morpho(ifirst, ilast, jfirst, jlast, dhsed )
    !&E--------------------------------------------------------------------------
    !&E                 ***  ROUTINE MUSTANG_morpho  ***
    !&E
    !&E ** Purpose : compute bathymetry variation dhsed (morphodynamic)
    !&E
    !&E ** Called by :  step, at the end
    !&E
    !&E--------------------------------------------------------------------------

    !! * Arguments
    INTEGER, INTENT(IN) :: ifirst, ilast, jfirst, jlast
    REAL(KIND=rsh),DIMENSION(ARRAY_DHSED),INTENT(INOUT) :: dhsed                        
    !! * Local declarations
    INTEGER                  :: i,j,k
    !!--------------------------------------------------------------------------
    !! * Executable part

! **TODO** keep these lines commented?? create a specific subroutine for zero gradient boundaries??
 ! if key_MUSTANG_bedload : choice of zero gradient at boundaries or no flux 
 ! if zero gradient at one open boundary : remove comment at this boundary
   ! south boundary
   !    IF (jfirst == JMIN_GRID) hsed(:,jfirst)=hsed(:,jfirst+1)
   ! north boundary
   !    IF (jlast == JMAX_GRID) hsed(:,jlast)=hsed(:,jlast-1)
   ! West boundary
   !    IF (ifirst == IMIN_GRID) hsed(ifirst,:)=hsed(ifirst+1,:)
   ! East boundary
   !    IF (ilast == IMAX_GRID) hsed(ilast,:)=hsed(ilast-1,:)

    DO j=jfirst,jlast
        DO i=ifirst,ilast
            hsed(i,j)=0.0_rsh
            DO k=ksmi(i,j),ksma(i,j)
                hsed(i,j)=hsed(i,j)+dzs(k,i,j)
            ENDDO
            IF (l_MF_dhsed) THEN
                dhsed(i,j)=(hsed_previous(i,j)-hsed(i,j))*MF_dhsed
            ELSE
                dhsed(i,j)=hsed_previous(i,j)-hsed(i,j)
            ENDIF
            hsed_previous(i,j)=hsed(i,j)
        ENDDO
    ENDDO

    t_morpho = t_morpho + dt_morpho

#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
!   bottom slope must be updated for bedload
    it_morphoYes = 1
#endif

  END SUBROUTINE MUSTANG_morpho      
#endif /* if defined MORPHODYN */
!=========================================================================== 
  
  SUBROUTINE sed_MUSTANG_outres(ifirst,ilast,jfirst,jlast,nv_out,h0_out,mask_h0,  &       
            var2D_ksma,var2D_tauskin,var2D_hsed)
                    
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_outres  ***
   !&E
   !&E ** Purpose : preparation of output sediment variables before writing
   !&E
   !&E
   !&E ** Description : 
   !&E  arguments IN : 
   !&E         loops  :ifirst,ilast,jfirst,jlast
   !&E         nv_out : number of variables that will be written
   !&E         h0_out : writing depths
   !&E
   !&E  arguments OUT:
   !&E         mask_h0 : mask for writing
   !&E         var2D.. : 2D variables to write 
   !&E
   !&E ** Called by :  sed_outMUSTANG_MARS
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
#ifdef key_BLOOM_insed
   USE bioloinit,  ONLY :  ndiag_tot,ndiag_3d_sed,diag_3d_sed,diag_2d_sed,ndiag_1d,ndiag_2d,ndiag_2d_sed
#endif

   !! * Arguments

   INTEGER, INTENT(IN)                                     :: ifirst,ilast,jfirst,jlast,nv_out
   REAL(KIND=riosh),DIMENSION(PROC_IN_ARRAY), INTENT(IN)   :: h0_out
   LOGICAL, DIMENSION(PROC_IN_ARRAY), INTENT(INOUT)          :: mask_h0
   REAL(KIND=riosh),  DIMENSION(PROC_IN_ARRAY),INTENT(OUT)   :: var2D_ksma,var2D_hsed,var2D_tauskin

   !! * Local declarations
   INTEGER        :: i,j,k,l
   REAL(KIND=rsh) :: eptmp,eptmp1,cstmp,epintbottom,dzs_estim,unitmudbinv,sumdzs
   
   !! * Executable part 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! calculation of masks at time t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     mask_h0(:,:) = .TRUE.
     DO j = jfirst, jlast
       DO i = ifirst, ilast
         IF ( h0_out(i,j) .NE. -valmanq .AND. ksma(i,j) > 0 ) mask_h0(i,j) = .FALSE.
       END DO
     END DO  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! preparation of var2D_ksma and var2D_hsed and var2D_tauskin arrays 
!!!!         number of sediment layer & total thickness of sediment & skinstress
!!!!         + additional variables if key_wave
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     var2D_ksma(PROC_IN_ARRAY) = -rg_valmanq_io
     var2D_hsed(PROC_IN_ARRAY) = -rg_valmanq_io
     var2D_tauskin(PROC_IN_ARRAY) = -rg_valmanq_io
#ifdef key_MUSTANG_specif_outputs
     varspecif2D_out(:,PROC_IN_ARRAY) = -rg_valmanq_io
     varspecif3Dnv_out(:,:,PROC_IN_ARRAY) = -rg_valmanq_io
#endif
#ifdef key_BLOOM_insed
     IF (l_out_subs_diag_sed) THEN
       var2D_diagsed(PROC_IN_ARRAY,:) = -rg_valmanq_io
     ENDIF
#endif

     DO j = jfirst, jlast
     DO i = ifirst, ilast
        DO k = ksma(i,j)+1, ksdmax 
         dzs(k,i,j) = -rg_valmanq_io
        END DO
        var2D_hsed(i,j) = 0.0_rsh
        IF(mask_h0(i,j)) THEN
           var2D_ksma(i,j) = -rg_valmanq_io
           var2D_hsed(i,j) = -rg_valmanq_io
           var2D_tauskin(i,j) = -rg_valmanq_io

#ifdef key_MUSTANG_specif_outputs
           varspecif2D_out(:,i,j) = -rg_valmanq_io
           DO k = 1, nv_out 
            varspecif3Dnv_out(:,k,i,j) = -rg_valmanq_io 
           END DO
#endif
#ifdef key_BLOOM_insed
           IF (l_out_subs_diag_sed) THEN
             var2D_diagsed(i,j,:) = -rg_valmanq_io
           ENDIF
#endif
        ELSE 
           var2D_ksma(i,j) = ksma(i,j)
           var2D_tauskin(i,j) = tauskin(i,j)

#ifdef key_MUSTANG_specif_outputs
           varspecif2D_out(:,i,j) = varspecif2D_save(:,i,j)
           DO k = 1, nv_out
            varspecif3Dnv_out(:,k,i,j) = varspecif3Dnv_save(:,k,i,j) 
            ENDDO
#endif
           sumdzs = 0.0_rsh
           DO k = ksdmin, ksma(i,j)
              sumdzs = sumdzs + dzs(k,i,j)
           END DO
           var2D_hsed(i,j) = sumdzs
#ifdef key_BLOOM_insed
           IF (l_out_subs_diag_sed) THEN
             var2D_diagsed(i,j,ndiag_1d+ndiag_2d-ndiag_2d_sed+1:ndiag_1d+ndiag_2d) = diag_2D_sed(ndiag_1d+ndiag_2d-ndiag_2d_sed+1:ndiag_1d+ndiag_2d,i,j)
           ENDIF
#endif
        ENDIF

     END DO
     END DO

        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! preparation of 3D arrays - interpolation according to the user's choice
!!!! var3D_dzs : thickness of sediment layers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     var3D_dzs(:,:,:) = 0.0_rsh
     CALL sed_MUSTANG_interpout_dzs(ifirst,ilast,jfirst,jlast,mask_h0,var3D_dzs,var2D_hsed)
     
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! preparation of 3D arrays - interpolation according to the user's choice
!!!! var3D_TEMP, var3D_SAL : temperature and salinity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(l_outsed_saltemp) THEN
       var3D_TEMP(:,:,:) = 0.0_rsh
       var3D_SAL(:,:,:) = 0.0_rsh
       unitmudbinv = 1._rsh
       CALL sed_MUSTANG_interpout_cvs(ifirst,ilast,jfirst,jlast,0,unitmudbinv, &
                                      mask_h0,cv_sed(-1,:,:,:),var3D_TEMP)
       CALL sed_MUSTANG_interpout_cvs(ifirst,ilast,jfirst,jlast,0,unitmudbinv,  &
                                       mask_h0,cv_sed(0,:,:,:),var3D_SAL)
     ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! preparation of 3D arrays - interpolation according to the user's choice
!!!! var3D_cvsed : substance concentrations in sediment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     DO l = 1, nv_out
         unitmudbinv = 1._rsh/unit_modif_mudbio_N2dw(irk_fil(l))
         var3D_cvsed(:,:,:,l) = 0.0_rsh
         CALL sed_MUSTANG_interpout_cvs(ifirst,ilast,jfirst,jlast,l,unitmudbinv,mask_h0,  &
                                        cv_sed(l,:,:,:),var3D_cvsed(:,:,:,l))
     ENDDO

#ifdef key_BLOOM_insed
     IF (l_out_subs_diag_sed) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! preparation of 3D arrays - interpolation according to the user's choice
!!!! var3D_diagsed : diagnostics variables 3D in sediment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       DO l = ndiag_tot-ndiag_3d_sed+1, ndiag_tot
         var3D_diagsed(:,:,:,l) = 0.0_rsh
         unitmudbinv = 1._rsh

         CALL sed_MUSTANG_interpout_cvs(ifirst,ilast,jfirst,jlast,0,unitmudbinv,mask_h0,  &
                                        diag_3d_sed(l,:,:,:),var3D_diagsed(:,:,:,l))

       ENDDO
     ENDIF
#endif
  
#ifdef key_MUSTANG_specif_outputs
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! preparation of 3D arrays - interpolation according to the user's choice
!!!! var3D_specifout : new output variables in sediment
!!!!   poro_save (and  poro_mud_save)
!!!! if key_MUSTANG_add_consol_outputs (version 2) : 
!!!!   loadograv_save, permeab_save, sigmapsg_save, dtsdzs_save, hinder_save
!!!!   sed_rate_save, sigmadjge_save, stateconsol_save

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DO l = 1, nv_out3Dk_specif
         unitmudbinv = 1._rsh
         var3D_specifout(:,:,:,l) = 0.0_rsh
         CALL sed_MUSTANG_interpout_cvs(ifirst,ilast,jfirst,jlast,0,unitmudbinv,mask_h0,  &
                                        varspecif3Dk_save(l,:,:,:),var3D_specifout(:,:,:,l))
     ENDDO
#endif

  END SUBROUTINE sed_MUSTANG_outres
  
!=========================================================================== 
  
  SUBROUTINE sed_MUSTANG_interpout_dzs(ifirst,ilast,jfirst,jlast,mask_h0,var3D_dzs,var2D_hsed)
              
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_interpout_dzs  ***
   !&E
   !&E ** Purpose : interpolation sediment thickness for output
   !&E
   !&E
   !&E ** Description : 
   !&E  arguments IN : 
   !&E         loops  :ifirst,ilast,jfirst,jlast
   !&E         mask_h0 : mask for writing
   !&E         var2D_hsed : total sediment thickness
   !&E
   !&E  arguments OUT:
   !&E         var3D_dzs : layers sediment thickness
   !&E
   !&E
   !&E ** Called by :  sed_MUSTANG_outres
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used


   !! * Arguments
   INTEGER, INTENT(IN)                                                 :: ifirst,ilast,jfirst,jlast
   LOGICAL,DIMENSION(PROC_IN_ARRAY), INTENT(IN)                        :: mask_h0
   REAL(KIND=riosh),DIMENSION(PROC_IN_ARRAY), INTENT(IN)                 :: var2D_hsed   
   REAL(KIND=riosh),DIMENSION(nk_nivsed_out,PROC_IN_ARRAY), INTENT(OUT)  :: var3D_dzs   


   !! * Local declarations
   INTEGER              :: i,j,k,kniv,finkniv
   REAL(KIND=rsh)       :: eptmp

   
   !! * Executable part 
     DO j = jfirst,jlast
       DO i = ifirst,ilast
         IF(mask_h0(i,j)) THEN
          var3D_dzs(:,i,j)=-rg_valmanq_io
         ELSE 
          IF(choice_nivsed_out == 1) THEN
              ! all layers are saved - k=1 bottom layer
              IF(ksma(i,j) >= 1 .AND. ksma(i,j) /= INT(rg_valmanq_io)) THEN
                var3D_dzs(ksdmin:ksma(i,j),i,j)=dzs(ksdmin:ksma(i,j),i,j)
              ELSE
                var3D_dzs(:,i,j)=-rg_valmanq_io
              END IF
              
          ELSE IF ( choice_nivsed_out == 2) THEN
            ! nk_nivsed_out layers are saved (from surface) - k=1 most deep layer
              IF(ksma(i,j) >= 1 .AND. ksma(i,j) /= INT(rg_valmanq_io)) THEN
                  DO kniv=1,nk_nivsed_out
                    k=ksma(i,j)-nk_nivsed_out+kniv
                    IF(k >= 1) THEN
                       var3D_dzs(kniv,i,j)=dzs(k,i,j)
                    ELSE
                       var3D_dzs(kniv,i,j)=-rg_valmanq_io
                    ENDIF
                  ENDDO
              ELSE
                  var3D_dzs(:,i,j)=-rg_valmanq_io
              END IF
              
          ELSE IF ( choice_nivsed_out == 3) THEN
              ! all layers till a max thickness given in paraMUSTANG.txt - k=1 : surface layer
              ! last layer is the layer whoch the bottom level exceeds the maximum thickness 
              IF(ksma(i,j) >= 1 .AND. ksma(i,j) /= INT(rg_valmanq_io)) THEN
                eptmp=0.0_rsh
                finkniv=0
                var3D_dzs(:,i,j)=0.
                kniv=1
                DO k = ksma(i,j),ksmi(i,j),-1
                  eptmp=eptmp+dzs(k,i,j)
                  IF(eptmp < epmax_nivsed_out) THEN
                    var3D_dzs(kniv,i,j)=dzs(k,i,j)
                    kniv=kniv+1
                  ELSE
                    finkniv=1
                    IF(kniv <= nk_nivsed_out) THEN
                      var3D_dzs(kniv,i,j)=var3D_dzs(kniv,i,j)+dzs(k,i,j)
                    ENDIF
                  ENDIF
                END DO
              ELSE
                var3D_dzs(:,i,j)=-rg_valmanq_io
              END IF
              
            ELSE IF ( choice_nivsed_out == 4) THEN
              ! a fixed number af layers is saved whith constant and fixed thickness given in paraMUSTANG.txt - k=1 : surface layer
              ! the last layer in addition integrates the results till the bottom sediment
              IF(ksma(i,j) >= 1 .AND. ksma(i,j) /= INT(rg_valmanq_io)) THEN
                eptmp=0.0_rsh
                kniv=1
                finkniv=0
                var3D_dzs(:,i,j)=-rg_valmanq_io
                k=ksma(i,j)
                DO WHILE (k >= ksdmin .AND. finkniv==0)
                  eptmp=eptmp+dzs(k,i,j)
                  IF(eptmp < ep_nivsed_outp1(kniv) .AND. kniv==nk_nivsed_out ) THEN
                    var3D_dzs(kniv,i,j)=var2D_hsed(i,j)-SUM(var3D_dzs(1:kniv-1,i,j))
                    finkniv=1
                  ELSE IF (eptmp < ep_nivsed_outp1(kniv) )THEN
                    var3D_dzs(kniv,i,j)=ep_nivsed_outp1(kniv)
                    eptmp=eptmp-ep_nivsed_outp1(kniv)
                    kniv=kniv+1  
                    IF(eptmp > 0 ) k=k+1
                  ELSE
                    DO WHILE (eptmp > ep_nivsed_outp1(kniv).AND. finkniv==0 ) 
                      eptmp=eptmp-ep_nivsed_outp1(kniv)
                      IF(kniv < nk_nivsed_out) THEN
                         var3D_dzs(kniv,i,j)=ep_nivsed_outp1(kniv)
                         kniv=kniv+1
                      ELSE                                ! kniv=nk_nivsed_out
                         var3D_dzs(kniv,i,j)=var2D_hsed(i,j)-SUM(var3D_dzs(1:kniv-1,i,j))
                         finkniv=1
                      ENDIF
                    ENDDO
                  ENDIF
                  k=k-1
                END DO
                IF(k==0 .AND. finkniv==0 ) THEN
                  IF(eptmp > 0.0_rsh) THEN
                    var3D_dzs(kniv,i,j)=eptmp
                    finkniv=1
                  ELSE
                    var3D_dzs(kniv-1,i,j)=var3D_dzs(kniv-1,i,j)+eptmp
                    finkniv=1
                  ENDIF
                ENDIF
              ELSE
                var3D_dzs(:,i,j)=-rg_valmanq_io
              END IF
            ENDIF
           ENDIF
         END DO 
       END DO
   
  END SUBROUTINE sed_MUSTANG_interpout_dzs
  
!=========================================================================== 
  
  SUBROUTINE sed_MUSTANG_interpout_cvs(ifirst,ilast,jfirst,jlast,iconvCV, &
                                        unitmudbinv,mask_h0,var3D,var3D_cvs)
              
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_interpout_cvs  ***
   !&E
   !&E ** Purpose : interpolation sediment concentration for output
   !&E
   !&E
   !&E ** Description : 
   !&E  arguments IN : 
   !&E         loops  :ifirst,ilast,jfirst,jlast
   !&E         mask_h0 : mask for writing
   !&E         var3D   : 3D variables to interpolate
   !&E         unitmudbinv   : 1/unit_modif_mudbio_N2dw
   !&E
   !&E  arguments OUT:
   !&E         var3D_cvs   : interpolated 3D variables to write
   !&E
   !&E         
   !&E
   !&E ** Called by  :  sed_MUSTANG_outres
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used

   !! * Arguments
   INTEGER, INTENT(IN)                                                 :: ifirst,ilast,jfirst,jlast,iconvCV
   LOGICAL,DIMENSION(PROC_IN_ARRAY), INTENT(IN)                        :: mask_h0
   REAL(KIND=rsh),DIMENSION(ksdmin:ksdmax,PROC_IN_ARRAY), INTENT(IN)   ::  var3D   
   REAL(KIND=riosh),DIMENSION(nk_nivsed_out,PROC_IN_ARRAY), INTENT(OUT)  ::  var3D_cvs  
   REAL(KIND=rsh),INTENT(IN)                                           :: unitmudbinv

   !! * Local declarations
   INTEGER              :: i,j,k,kniv,finkniv
   REAL(KIND=rsh)       :: eptmp,eptmp1,cstmp,epintbottom

   !! * Executable part
    
   DO j = jfirst,jlast
     DO i = ifirst,ilast

         IF(mask_h0(i,j)) THEN
               var3D_cvs(:,i,j)=-rg_valmanq_io
         ELSE 
         
              IF(choice_nivsed_out == 1) THEN
              ! all layers are saved - k=1 bottom layer
                IF(ksma(i,j) >= 1 .AND. ksma(i,j) /= INT(rg_valmanq_io)) THEN
                  var3D_cvs(ksdmin:ksma(i,j),i,j)=var3D(ksdmin:ksma(i,j),i,j)*unitmudbinv
!#if defined key_BLOOM_insed
!                  IF(iconvCV .GE. nvpc+1 .AND. iconvCV .LE. nvp) THEN
!                  ! repassage en micromol/gSedSec des variables NoCP (only concentrations, not diagnostics)
!                    var3D_cvs(ksdmin:ksma(i,j),i,j)=var3D_cvs(ksdmin:ksma(i,j),i,j)/c_sedtot(ksdmin:ksma(i,j),i,j)
!                  ENDIF
!#endif
                ELSE
                  var3D_cvs(:,i,j)=-rg_valmanq_io
                END IF
                
              ELSE IF (choice_nivsed_out == 2) THEN
             ! nk_nivsed_out layers are saved (from surface) - k=1 most deep layer
                IF(ksma(i,j) >= 1 .AND. ksma(i,j) /= INT(rg_valmanq_io)) THEN
                   DO kniv=1,nk_nivsed_out
                    k=ksma(i,j)-nk_nivsed_out+kniv
                    IF (k >= 1) THEN
                      var3D_cvs(kniv,i,j)=var3D(k,i,j)*unitmudbinv
                    ELSE
                      var3D_cvs(kniv,i,j)=-rg_valmanq_io
                    ENDIF
                  ENDDO
                ELSE
                  var3D_cvs(:,i,j)=-rg_valmanq_io
                END IF
              ELSE IF (choice_nivsed_out == 3) THEN
              ! all layers till a max thickness given in paraMUSTANG.txt
              ! last layer is the layer which the bottom level exceeds the maximum thickness 
                IF(ksma(i,j) >= 1 .AND. ksma(i,j) /= INT(rg_valmanq_io)) THEN
                  eptmp=0.0_rsh
                  finkniv=0
                  var3D_cvs(:,i,j)=0.0_rsh
                  kniv=1
                  epintbottom=0.0_rsh
                  DO k = ksma(i,j),ksmi(i,j),-1
                    eptmp=eptmp+dzs(k,i,j)
                    IF(eptmp < epmax_nivsed_out) THEN
                      var3D_cvs(kniv,i,j)=var3D(k,i,j)*unitmudbinv
                      kniv=kniv+1
                    ELSE
                      finkniv=1
                      epintbottom=epintbottom+dzs(k,i,j)
                      IF(kniv <= nk_nivsed_out) THEN
                        var3D_cvs(kniv,i,j)=var3D_cvs(kniv,i,j)+var3D(k,i,j)*dzs(k,i,j)*unitmudbinv
                      ENDIF
                    ENDIF
                  END DO
                  IF(finkniv==1) var3D_cvs(kniv,i,j)=var3D_cvs(kniv,i,j)/epintbottom
                ELSE
                  var3D_cvs(:,i,j)=-rg_valmanq_io
                END IF
              ELSE IF (choice_nivsed_out == 4) THEN
              ! a fixed number af layers is saved whith constant et fixed thickness given in paraMUSTANG.txt
              ! the last layer in addition integrates the results till the bottom sediment
                IF(ksma(i,j) >= 1 .AND. ksma(i,j) /= INT(rg_valmanq_io)) THEN
                  eptmp1=0.0_rsh
                  eptmp=0.0_rsh
                  cstmp=0.0_rsh
                  kniv=1
                  finkniv=0
                  var3D_cvs(:,i,j)=-rg_valmanq_io
                  k=ksma(i,j)
                  DO WHILE (k >= ksdmin .and. finkniv ==0)
                    eptmp=eptmp+dzs(k,i,j)
                    IF(eptmp < ep_nivsed_outp1(kniv) .AND. kniv==nk_nivsed_out) THEN
                      var3D_cvs(kniv,i,j)=(cstmp+SUM(var3D(ksdmin:k,i,j)*dzs(ksdmin:k,i,j)*unitmudbinv))  &
                                            /(SUM(dzs(ksdmin:k,i,j))+eptmp-dzs(k,i,j))
                      finkniv=1
                    ELSE
                     IF(eptmp <= ep_nivsed_outp1(kniv)) THEN
                         cstmp=cstmp+dzs(k,i,j)*var3D(k,i,j)*unitmudbinv 
                         eptmp1=eptmp
                     ENDIF
                     DO WHILE (eptmp > ep_nivsed_outp1(kniv) .and. finkniv ==0 ) 
                        cstmp=cstmp+(ep_nivsed_outp1(kniv)-eptmp1)*var3D(k,i,j)*unitmudbinv
                        eptmp=eptmp-ep_nivsed_outp1(kniv)
                        IF(kniv < nk_nivsed_out) THEN
                          var3D_cvs(kniv,i,j)=cstmp/(ep_nivsed_outp1(kniv))
                          kniv=kniv+1
                          IF(eptmp > ep_nivsed_outp1(kniv)) THEN
                             cstmp=0.0_rsh
                             eptmp1=0.0_rsh
                          ELSE
                             eptmp1=eptmp
                             cstmp=eptmp1*var3D(k,i,j)*unitmudbinv
                          ENDIF
                        ELSE                                ! kniv=nk_nivsed_out
                          IF(k==1) THEN
                            var3D_cvs(kniv,i,j)=var3D(k,i,j)*unitmudbinv
                          ELSE
                            var3D_cvs(kniv,i,j)=(cstmp+SUM(var3D(ksdmin:k-1,i,j)*dzs(ksdmin:k-1,i,j)*unitmudbinv))  &
                                                 /(eptmp+SUM(dzs(ksdmin:k-1,i,j)))
                          ENDIF
                          finkniv=1
                        ENDIF
                      ENDDO 
                    ENDIF
                    k=k-1
                  END DO
                  IF(k==0 .AND. finkniv==0 .AND. eptmp > 0.0_rsh) THEN
                    var3D_cvs(kniv,i,j)=var3D(1,i,j)*unitmudbinv
                    finkniv=1
                  ENDIF
                ELSE
                  var3D_cvs(:,i,j)=-rg_valmanq_io
                END IF
              ENDIF
             ENDIF
           END DO 
         END DO
  
  END SUBROUTINE sed_MUSTANG_interpout_cvs
  
   !!==============================================================================

   SUBROUTINE sed_MUSTANG_comp_z0sed(ifirst,ilast,jfirst,jlast,BATHY_H0)

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_comp_z0sed ***
   !&E
   !&E ** Purpose : computes sediment roughness length z0sed in meter
   !&E               depending on sediment diameter
   !&E
   !&E ** Called by :  MUSTANG_update
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used

   !! * Arguments
   INTEGER, INTENT(IN)                                  :: ifirst,ilast,jfirst,jlast
   REAL(KIND=rsh),DIMENSION(ARRAY_BATHY_H0),INTENT(IN)  :: BATHY_H0

   !! * Local declarations
   INTEGER        :: i,j,ksmax
   REAL(KIND=rsh) :: diamgravsan
   !!---------------------------------------------------------------------------
   !! * Executable part

     DO j=jfirst,jlast
       DO i=ifirst,ilast
         IF(BATHY_H0(i,j).EQ.-valmanq) THEN
           ! Here value not set to -valmanq, since in skinstress
           ! computation, there are calls in i+1, i-1, j+1, j-1 whihtout test on
           ! land : in case of land velocities =0 ; and velocities are averaged
           z0sed(i,j) = z0sedmud
         ELSE
           ksmax = ksma(i,j)
           IF((ksmax.NE.0).AND. c_sedtot(ksmax,i,j).NE.0.0_rsh)THEN !There is sediment
             diamgravsan = MAX(diam_sed(isand2),  &
                          SUM( (cv_sed(igrav1:isand2,ksmax,i,j)/c_sedtot(ksmax,i,j))*diam_sed(igrav1:isand2) ))
             ! Application of Nikurasde z0 = (d50/12)
             ! If only mud diamgravsan = 0., so z0seduni is used
             z0sed(i,j) = MAX(z0seduni,diamgravsan/12.0_rsh)
           ELSE ! There is no sediment
             z0sed(i,j) = z0sedbedrock
           ENDIF
         ENDIF 
#if defined key_MUSTANG_specif_outputs        
         varspecif2D_save(10,i,j)=z0sed(i,j)
#endif
       ENDDO
     ENDDO

   END SUBROUTINE sed_MUSTANG_comp_z0sed

   !!==============================================================================

   SUBROUTINE sed_MUSTANG_comp_z0hydro(ifirst,ilast,jfirst,jlast,z0hydro)

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_comp_z0hydro ***
   !&E
   !&E ** Purpose : computes hydrodynamical roughness length depending on sediment diameter
   !&E
   !&E ** Description : use arguments and common variable 
   !&E
   !&E  arguments IN : 
   !&E         loops  :ifirst,ilast,jfirst,jlast
   !&E          
   !&E  variables OUT:
   !&E         z0hydro : roughness length (m)
   !&E
   !&E
   !&E ** Called by :  initMUSTANG and MUSTANG_update
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used

   !! * Arguments
   INTEGER, INTENT(IN)                                     :: ifirst,ilast,jfirst,jlast
   REAL(KIND=rsh),DIMENSION(ARRAY_Z0HYDRO),INTENT(INOUT)   :: z0hydro                         

   !! * Local declarations
   INTEGER        :: i,j,k,iv
   REAL(KIND=rsh) :: diamsan,somsan,hsedz0,sommud
   !!---------------------------------------------------------------------------
   !! * Executable part

     DO j=jfirst,jlast
       DO i=ifirst,ilast
           IF(ksma(i,j) .GT. ksmi(i,j))THEN
              diamsan=0.0_rsh
              somsan=0.0_rsh
              hsedz0=0.0_rsh
              sommud=0.0_rsh
              k=ksma(i,j)
              !  thickness of 1 cm for sand diameter evaluation( 1 cm superficial) 
              DO WHILE (hsedz0 .LE. 0.01_rsh .AND. k .GT. (ksmi(i,j)+1))
                DO iv=isand1,isand2
                      somsan=somsan+cv_sed(iv,k,i,j)*dzs(k,i,j)
                      diamsan=diamsan+diam_sed(iv)*cv_sed(iv,k,i,j)*dzs(k,i,j)
                END DO
                DO iv=imud1,imud2
                      sommud=sommud+cv_sed(iv,k,i,j)*dzs(k,i,j)
                ENDDO
                hsedz0=hsedz0+dzs(k,i,j)
                k=k-1
              END DO
              IF (sommud/(somsan+sommud+epsilon_MUSTANG) .LT. 0.3_rsh)  THEN
                 diamsan=MAX(diamsan/(somsan+epsilon_MUSTANG),diam_sed(isand2))
                 z0hydro(i,j)=coef_z0_coupl*diamsan  
              ELSE
                 z0hydro(i,j)=z0_hydro_mud
              ENDIF
           ELSE
              z0hydro(i,j)=z0_hydro_bed
           ENDIF
 
#if defined key_MUSTANG_specif_outputs        
           varspecif2D_save(15,i,j)=z0hydro(i,j)
#endif
       ENDDO
     ENDDO

   END SUBROUTINE sed_MUSTANG_comp_z0hydro

   !!==============================================================================

  SUBROUTINE sed_MUSTANG_depflx(ifirst, ilast, jfirst, jlast)
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_depflx  ***
   !&E
   !&E ** Purpose : computes deposition flux trends  before advection computations
   !&E
   !&E ** Description : use arguments and common variable 
   !&E
   !&E  arguments IN : 
   !&E         loops  :ifirst,ilast,jfirst,jlast
   !&E          
   !&E  variables OUT:
   !&E         flx_w2s : deposit trends (m/s) ranged in comMUSTANG
   !&E         
   !&E  need to be know by hydrodynamic code:
   !&E         kmax=NB_LAYER_WAT  (known from coupleur_dimhydro_MUSTANG.h)
   !&E
   !&E  need to be know by code treated substance 
   !&E         igrav2,isand2, nvpc, nvp, nv_adv : 
   !&E         irkm_var_assoc (ivp)
   !&E         tocd(iv)  
   !&E
   !&E ** Called by :  MUSTANG_update
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
 

   !! * Arguments
   INTEGER, INTENT(IN)  :: ifirst, ilast, jfirst, jlast
   
   !! * Local declarations
   INTEGER                          :: ivp, i, j
   REAL(KIND=rlg)                   :: tocdpe

   !!---------------------------------------------------------------------------
   !! * Executable part

      DO j = jfirst, jlast
      DO i = ifirst, ilast

       ! Separation of the loops because not necessarily the same unit for non constitutive variables
       ! no calculation for gravels now
       DO ivp = igrav2+1, nvpc
         IF(ivp .LE. isand2) THEN 
            ! Sanford et Halka (1993), Winterwerp (2004), Dufois (2008)
            flx_w2s(ivp, i, j) = ws3_bottom_MUSTANG(ivp, i, j) * fwet(i, j)
         ELSE
            tocdpe = tocd(ivp) + 0.00001_rsh
            flx_w2s(ivp, i, j) = ws3_bottom_MUSTANG(ivp, i, j) * fwet(i, j)  &
                          * MAX(0.0_rsh, 1.0_rsh - (tauskin(i, j) / tocdpe)) &  
                          * tocd(ivp) / tocdpe
         ENDIF
       ENDDO
#ifdef key_Pconstitonly_insed
       DO ivp = nvpc+1, nvp
          flx_w2s(ivp, i, j) = 0.0_rsh
       ENDDO
#else
       DO ivp = nvpc+1, nvp
         IF(irkm_var_assoc(ivp) > 0 .AND. irkm_var_assoc(ivp) < imud1) THEN
           ! substance sorbed on sand
           flx_w2s(ivp, i, j) = ws3_bottom_MUSTANG(ivp, i, j) * fwet(i, j)
         ELSE
           tocdpe = tocd(ivp) + 0.00001_rsh
           flx_w2s(ivp, i, j) = ws3_bottom_MUSTANG(ivp, i, j) * fwet(i, j)  &
                          * MAX(0.0_rsh,1.0_rsh - (tauskin(i, j) / tocdpe)) &  
                          * tocd(ivp) / tocdpe
         ENDIF
         IF(irkm_var_assoc(ivp) > 0 ) THEN
          IF (flx_w2s(irkm_var_assoc(ivp),i,j) == 0.0_rsh) flx_w2s(ivp, i, j) = 0.0_rsh
         END IF
       ENDDO
#endif

       ENDDO
       ENDDO
         
  END SUBROUTINE sed_MUSTANG_depflx

!!==============================================================================

  SUBROUTINE sed_MUSTANG_sandconcextrap(ifirst, ilast, jfirst, jlast)

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_sandconcextrap  ***
   !&E
   !&E ** Purpose : interpolates sand concentration in the first layer
   !&E              assuming a Rouse profile 
   !&E
   !&E ** Description : use arguments and common variable 
   !&E
   !&E  arguments IN : 
   !&E         loops  :ifirst,ilast,jfirst,jlast
   !&E
   !&E  variables OUT : (in comMUSTANG) 
   !&E     flx_w2s : tendances aux depots corriges pour sables 
   !&E     corflux, corfluy : corrections des flux horizontaux pour les sables
   !&E     
   !&E     
   !&E  need to be know by hydrodynamic code:
   !&E         kmax=NB_LAYER_WAT (known from coupleur_dimhydro_MUSTANG.h)
   !&E         alt_cw1 , htot: evaluated in coupleur_conv2MUSTANG
   !&E         
   !&E  need to be know by code treated substance:
   !&E         isand1,isand2,nv_adv 
   !&E     
   !&E  use module MUSTANG variables  :
   !&E         aref_sand
   !&E         ustarbot (evaluated in sed_gravit)
   !&E
   !&E ** Called by :  MUSTANG_update
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used

   !! * Arguments
   INTEGER, INTENT(IN)                        :: ifirst, ilast, jfirst, jlast
   !! * Local declarations
   INTEGER                              :: izz,ivp,i,j,ivpp
   REAL(KIND=rsh)                       :: altc1,extrap,rouse,som,dzcche,hzed,     &
                                           einstein,alogaltc1sz0,ashmapwrouse,rouse1
   INTEGER,PARAMETER                    :: nbcche=20
   REAL(KIND=rsh),DIMENSION(1:nbcche+1) :: hzi
#ifdef key_sand2D
   REAL(KIND=rsh),DIMENSION(1:nbcche+1) :: hzisdbot
#endif


   !!---------------------------------------------------------------------------
   !! * Executable part


    corflux(:,:,:)=1.0_rsh
    corfluy(:,:,:)=1.0_rsh

      DO j=jfirst,jlast
      DO i=ifirst,ilast

          altc1=alt_cw1(i,j)
          IF(altc1 .LE. aref_sand .OR. htot(i,j) .LE. h0fond)THEN
            extrap=1.0_rsh
            ! in this case, flx_w2s(ivp,i,j) is not changed
            DO ivp=isand1,isand2
              corflux(ivp,i,j)=1.0_rsh
              corfluy(ivp,i,j)=1.0_rsh
            ENDDO
          ELSE
            alogaltc1sz0=LOG(altc1/z0sed(i,j))
            hzi(1)=epn_bottom_MUSTANG(i,j)
            DO izz=1,nbcche
              hzi(izz+1)=hzi(izz)-(hzi(1)-aref_sand)/2.0_rsh**izz
            ENDDO     
            hzi(nbcche+1)=aref_sand
#ifdef key_sand2D
            hzisdbot(1)=htot(i,j)
            DO izz=1,nbcche
              hzisdbot(izz+1)=hzisdbot(izz)-(hzisdbot(1)-aref_sand)/2.0_rsh**izz
            ENDDO     
            hzisdbot(nbcche+1)=aref_sand
#endif
            DO ivp=isand1,isand2
            ! Rouse number depends on settling velocity, then on the sand variable
            ! modif Julie 2012 to taking into account flat bottom
              IF(ws3_bottom_MUSTANG(ivp,i,j)>1.34_rsh*ustarbot(i,j)) THEN
                 rouse=1.2_rsh
              ELSE
                rouse1=ws3_bottom_MUSTANG(ivp,i,j)/ustarbot(i,j)
                IF(rouse1>=0.75_rsh) THEN
                   rouse=0.353_rsh*rouse1+0.727
                ELSEIF(rouse1>0.1_rsh .AND. rouse1<0.75_rsh) THEN
                   rouse=rouse1/(0.4_rsh+0.8_rsh*rouse1*rouse1)   ! rouse=rouse1/(0.4*beta) avec beta=1+2*rouse1*rouse1
                ELSE
                  rouse=rouse1/0.4_rsh
                ENDIF
              ENDIF
            ! integral concentration and concentration*vitesse log
            ! in first layer 
              som=0.0_rsh
              einstein=0.0_rsh
              DO izz=1,nbcche
#ifdef key_sand2D
               IF(l_subs2D(ivp)) THEN
                hzed=0.5_rsh*(hzisdbot(izz)+hzisdbot(izz+1))
                dzcche=hzisdbot(izz)-hzisdbot(izz+1) 
               ELSE
#endif
                hzed=0.5_rsh*(hzi(izz)+hzi(izz+1))
                dzcche=hzi(izz)-hzi(izz+1) 
#ifdef key_sand2D
               ENDIF
#endif
                som=som+dzcche*((htot(i,j)-hzed)/hzed)**rouse
                einstein=einstein+dzcche*((htot(i,j)-hzed)/hzed)**rouse      &  
                                  *LOG(hzed/z0sed(i,j))
              ENDDO
#ifdef key_sand2D
              IF(l_subs2D(ivp)) THEN
                rouse2D(ivp,i,j)=rouse
                sum_tmp(ivp,i,j)=som
              ENDIF
#endif 
!!
              ashmapwrouse=(aref_sand/(htot(i,j)-aref_sand))**rouse
              som=som*ashmapwrouse

              einstein=einstein*ashmapwrouse/som
              ! the computation of einstein and USE of corflux should be checked
              !??? modif P.LeHir              if(hex(i,j).le.hm)then
              !??? modif P.LeHir              if(hey(i,j).le.hm)then

#ifdef key_sand2D
              IF(l_subs2D(ivp)) THEN
! with SAN2D in CROCO, transport is done in the bottom layer 
! (in MARS, we consider the full depth although advection is done using the bottom velocity)
! i.e. in CROCO, C is representative of the bottom layer, 
! in MARS C is representative of the depth averaged concentration
                extrap=(hzi(1)-aref_sand)/som
                corflux(ivp,i,j)=einstein/alogaltc1sz0
                corfluy(ivp,i,j)=einstein/alogaltc1sz0
              ELSE
#endif
                extrap=(hzi(1)-aref_sand)/som

                corflux(ivp,i,j)=einstein/alogaltc1sz0
                corfluy(ivp,i,j)=einstein/alogaltc1sz0
                
#ifdef key_sand2D
              ENDIF
#endif 
              flx_w2s(ivp,i,j)=flx_w2s(ivp,i,j)*extrap 
              ! for substances which are sorbed on sand 
              DO ivpp=nvpc+1,nvp
               IF(irkm_var_assoc(ivpp) .EQ. ivp) THEN
                flx_w2s(ivpp,i,j)=flx_w2s(ivpp,i,j)*extrap
               ENDIF
              ENDDO             
            ENDDO              
          ENDIF


       ENDDO
       ENDDO
         
  END SUBROUTINE sed_MUSTANG_sandconcextrap
      
!!==============================================================================
      
#ifdef key_MUSTANG_V2
  SUBROUTINE sed_MUSTANG_erosion(ifirst, ilast, jfirst, jlast, dtinv, &
#if defined key_MUSTANG_lateralerosion || defined key_MUSTANG_bedload
                                    BAROTROP_VELOCITY_U, BAROTROP_VELOCITY_V, &
#endif
                                    dt_true) 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_erosion version V2  ***
   !&E
   !&E ** Purpose : computes erosion fluxes (version B. Mengual et P. Le Hir)
   !&E
   !&E ** Description : 
   !&E       arguments IN :
   !&E          loops  :ifirst,ilast,jfirst,jlast
   !&E          dtinv, dt_true : 1/dt  and dt
   !&E          u,v :   BAROTROP_VELOCITY_U, BAROTROP_VELOCITY_V
   !&E      
   !&E     variables OUT : 
   !&E          flx_s2w : concentration flux sediment to water
   !&E          phieau_s2w : water flux  sediment to water
   !&E          the sediment has been remodelled
   !&E          cv_sed, dzs, ksma, poro have changed
   !&E          flx_bx, flx_by : bedload fluxes
   !&E
   !&E
   !&E ** Called by :  MUSTANG_update
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used


   INTEGER, INTENT(IN)                        :: ifirst, ilast, jfirst, jlast
   REAL(KIND=rsh),INTENT(IN)                  :: dtinv
   REAL(KIND=rlg),INTENT(IN)                  :: dt_true  ! =halfdt in MARS
#if defined key_MUSTANG_lateralerosion || defined key_MUSTANG_bedload
   REAL(KIND=rsh),DIMENSION(ARRAY_VELOCITY_U),INTENT(IN)   :: BAROTROP_VELOCITY_U
   REAL(KIND=rsh),DIMENSION(ARRAY_VELOCITY_V),INTENT(IN)   :: BAROTROP_VELOCITY_V 
#endif

   !! * Local declarations
   INTEGER        ::  i,j,k,iv,ksmax,ksup,ksmaxa,isplit
   REAL(KIND=rsh) ::  dt1,toce,cvolgrv,csanmud,phieau_ero_ij,             &
                      erodab,ero,erosi,ddzs,dzsa,dzsisu,dzsam1,dflx_s2w,  &
                      cvolp,volerod,poroa,poroam1,dflusve,        &
                      xeros, sed_eros_flx,excespowr, &
                      heauw,heaue,heaun,heaus,heau_milieu,eroe,erow,eros,eron
   REAL(KIND=rsh) :: diamgravsan,somgravsan,frmudcr1,cv_sed_tot,dzs_activelayer_ij, &
                     dt_ero_max,dts2,mass_tot,cvolgrvsan,sommud,  &
                     niter_ero_noncoh,niter_ero_coh,isthere_erosion
   REAL(KIND=rsh) :: dzs_ini,poro_ini
   REAL(KIND=rsh),DIMENSION(1:nvpc) :: frac_sed
   REAL(KIND=rsh),DIMENSION(1:nvp) :: mass_sed,sed_eros_flx_class_by_class,dt_ero
   REAL(KIND=rsh),DIMENSION(1:nvp) :: flx_bxij,flx_byij
   REAL(KIND=rsh),DIMENSION(-1:nv_tot) :: cv_sed_ini
   REAL(KIND=rsh),DIMENSION(-1:nv_adv)  ::  flx_s2w_eroij
#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
   REAL(KIND=rsh) :: porowater1,porowater_new
#endif

   !!---------------------------------------------------------------------------
   !! * Executable part

    DO j=jfirst,jlast
      DO i=ifirst,ilast

        flx_s2w(-1,i,j)=0.0_rsh
        flx_s2w( 0,i,j)=0.0_rsh
        DO iv=1,nv_adv
          flx_s2w(iv,i,j)=0.0_rsh
        ENDDO
        phieau_ero_ij=0.0_rsh
        flx_s2w_eroij(:)=0.0_rsh
 
        DO k=ksmi(i,j),ksma(i,j)
          sommud=0.0_rsh
          DO iv=imud1,imud2
            sommud=sommud+cv_sed(iv,k,i,j)
          ENDDO
          cvolgrvsan=0.0_rsh
          DO iv=igrav1,isand2
            cvolgrvsan=cvolgrvsan+cv_sed(iv,k,i,j)/ros(iv)
          ENDDO
          crel_mud(k,i,j)=sommud/(1.0_rsh-cvolgrvsan)
        END DO

#ifdef key_MUSTANG_specif_outputs
       ! Outputs >--
        varspecif2D_save(1,i,j)=sommud/(c_sedtot(ksma(i,j),i,j)+epsi30_MUSTANG) ! frmudsup last sommud computed in k=ksma(i,j)
        varspecif3Dnv_save(1,:,i,j)=0.0_rsh  ! toce_save             
        varspecif2D_save(2:4,i,j)=0.0_rsh !(dzs_ksmax_save,dzs_aclay_comp_save,dzs_aclay_kept_save
        varspecif3Dnv_save(4,:,i,j)=0.0_rsh  ! pephm_fcor_save
#ifdef key_MUSTANG_bedload
        varspecif3Dnv_save(8,:,i,j)=0.0_rsh  ! fsusp_save
#endif
       ! --<
#endif       
        ksmax=ksma(i,j)
#ifdef key_MUSTANG_debug
               IF (l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
                 print *,'  > deb erosion',i,j
                 print *,'  t=',CURRENT_TIME, ' ksmax=',ksmax
                 print *,'  dzs(ksmax-3:ksmax,i,j)=',dzs(ksmax-3:ksmax,i,j)
                 print *,'  cv_sed(:,ksmax-3,i,j)=',cv_sed(:,ksmax-3,i,j)
                 print *,'  cv_sed(:,ksmax-2,i,j)=',cv_sed(:,ksmax-2,i,j)
                 print *,'  cv_sed(:,ksmax-1,i,j)=',cv_sed(:,ksmax-1,i,j)
                 print *,'  cv_sed(:,ksmax,i,j)=',cv_sed(:,ksmax,i,j)
               END IF
#endif 
        ero=0.0_rsh
        ! niter put in real beacause of ratio to estimate statistics 
        niter_ero_noncoh=0.0_rsh
        niter_ero_coh=0.0_rsh

        IF(ksmax > ksmi(i,j))THEN  ! this means that the ksmi layer is never eroded

            dt1=REAL(dt_true,rsh)
            k=ksmax
                       
      2     CONTINUE

          IF (.NOT. l_eroindep_noncoh) THEN 

            ! sediment always eroded as a mixture
            frmudcr1=0.0_rsh ! sediment always eroded as a mixture
            diamgravsan=0.0_rsh

          ELSE
            !!! Test if sediment et ksmax is cohesif or not
            !!! ===========================================
            diamgravsan=0.0_rsh
            somgravsan=0.0_rsh
            DO iv=igrav1,isand2
              diamgravsan=diamgravsan+diam_sed(iv)*cv_sed(iv,ksmax,i,j)
              somgravsan=somgravsan+cv_sed(iv,ksmax,i,j)
            ENDDO
            IF (isand2>0 .AND. somgravsan>0.0_rsh) THEN
              diamgravsan=MAX(diamgravsan/(somgravsan+epsi30_MUSTANG),diam_sed(isand2))
            END IF
            frmudcr1=MIN(coef_frmudcr1*diamgravsan,frmudcr2)
            
          END IF
          l_isitcohesive(i,j)=isitcohesive(cv_sed(:,ksmax,i,j),frmudcr1)
!
#ifdef key_MUSTANG_debug
          IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
              print *,''
              print *,' ************************'
              print *,' ENTER sed_erosion_V2    '
              print *,' ************************'
              print *,''
              print *,'    diamgravsan / frmudcr1 / l_isitcohesive(i,j) = ',diamgravsan,frmudcr1,l_isitcohesive(i,j)              
              print *,'    dt1=',dt1
              print *,'    ksmax et dzs=',ksmax,dzs(ksmax,i,j)
              print *,''
           END IF
#endif
           IF  (.NOT. l_isitcohesive(i,j)) THEN

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !!! CASE 1: NON COHESIVE SEDIMENT --> EROSION CLASS BY CLASS (+ BEDLOAD, not operational) !!!
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           

                 !print *,''
                 !print *,'CASE 1: NON COHESIVE SEDIMENT --> EROSION CLASS BY CLASS'

               niter_ero_noncoh=niter_ero_noncoh+1.0_rsh

              CALL sed_MUSTANG_comp_tocr_mixsed(ksmax, i, j, xeros, excespowr, toce)

              IF(tauskin(i, j) .GT. toce)THEN

              
#ifdef key_MUSTANG_debug
               IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
                 print *,'    TAUSKIN =',tauskin(i,j),' > TOCE=',TOCE
                 print *,'    Carac of ksmax layer before managing active layer : '
                 print *,'      ksmax=',ksmax
                 print *,'      dzs(ksmax,i,j)=',dzs(ksmax,i,j)
                 print *,'      cv_sed(:,ksmax,i,j)=',cv_sed(:,ksmax,i,j)
                 print *,'      c_sedtot(ksmax,i,j)=',c_sedtot(ksmax,i,j)
                 print *,'      poro(ksmax,i,j)=',poro(ksmax,i,j),'poro_mud(ksmax,i,j)=',poro_mud(ksmax,i,j) 
                 print *,'      crel_mud(ksmax,i,j)=',crel_mud(ksmax,i,j)
               END IF
#endif
               ! IN : i,j,ksmax / OUT : active layer: updates of ksmax, dzs, cv_sed, c_sedtot, poro
               CALL MUSTANGV2_manage_active_layer(i,j,ksmax  &
#if ! defined key_noTSdiss_insed 
                        ,flx_s2w_eroij                                &
#endif
#if ! defined key_nofluxwat_IWS
                        ,phieau_ero_ij                                &
#endif
                        )

              END IF  


#ifdef key_MUSTANG_debug
              IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
                 print *,'    Carac of ksmax layer after managing active layer : '
                 print *,'      ksmax=',ksmax
                 print *,'      dzs(ksmax,i,j)=',dzs(ksmax,i,j)
                 print *,'      cv_sed(:,ksmax,i,j)=',cv_sed(:,ksmax,i,j)
                 print *,'      c_sedtot(ksmax,i,j)=',c_sedtot(ksmax,i,j)
                 print *,'      poro(ksmax,i,j)=',poro(ksmax,i,j),'poro_mud(ksmax,i,j)=',poro_mud(ksmax,i,j) 
                 print *,'      crel_mud(ksmax,i,j)=',crel_mud(ksmax,i,j)
               END IF
#endif
              ksmaxa=ksmax ! needs to be memorised
              poro_ini=poro(ksmax,i,j)
              dzs_ini=dzs(ksmax,i,j)
              cv_sed_ini(:)=cv_sed(:,ksmax,i,j)

#ifdef key_MUSTANG_bedload 
              ! IN : i,j,ksmax / OUT : flx_bxij,flx_byij (bedload Flux in kg/m/s)

              CALL MUSTANGV2_eval_bedload(i, j, ksmax, flx_bxij, flx_byij)  
#else
              flx_bxij(:) = 0.0_rsh
              flx_byij(:) = 0.0_rsh
#endif

              ! On calcule un flux derosion pour chaque classe (en kg/m2/s)
              ! IN : i,j,ksmax / OUT : sed_eros_flx_class_by_class(iv)

              CALL MUSTANGV2_comp_eros_flx_indep(i,j,ksmax,        &
#ifdef key_MUSTANG_bedload
                                        CELL_DX,CELL_DY,flx_bxij,flx_byij,       &
                                        BAROTROP_VELOCITY_U,BAROTROP_VELOCITY_V, &
#endif
                                        sed_eros_flx_class_by_class)


              ! IN : i,j,dt1
              ! INOUT : ksmax,flx_bxij,flx_byij,sed_eros_flx_class_by_class 
                        !!! Attention [flx_bxij,flx_byij] IN : kg/m/s --> OUT : kg
              !             [sed_eros_flx_class_by_class] IN : kg/m2/s --> OUT : kg
              ! OUT : dt_ero
              !       update dzs, cv_sed, c_sedtot, and poro in the ksmax (potentially changed) layer after erosion

#ifdef key_MUSTANG_debug
                 IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
                    print *,'    SED BORNE AND APPLY EROSION TOT'
                 END IF
#endif

              ! on ne fait rien si aucune variable particulaire constitutive ne bouge
              !  mais probleme pour les variables particulaires non constitutives non associees 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!?????????????????????????????????
              isthere_erosion=sum(flx_bxij(1:nvpc))+sum(flx_byij(1:nvpc))+sum(sed_eros_flx_class_by_class(1:nvpc))
              IF(isthere_erosion .NE. 0._rsh) THEN
                CALL MUSTANGV2_borne_and_apply_erosion_tot(i,j,ksmax,flx_bxij,flx_byij, &
#if ! defined key_noTSdiss_insed 
                        flx_s2w_eroij,                                &
#endif
#if ! defined key_nofluxwat_IWS
                        phieau_ero_ij,                                &
#endif
                        sed_eros_flx_class_by_class,dt1,dt_ero)

                !!  ==> erosion of one layer or elimination of the entire layer

#ifdef key_MUSTANG_debug
                IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
                  print *,'    INTEGRATION of FLX_BX/Y et FLX_S2W'
                END IF
#endif

                DO iv=1,nvp

                  flx_s2w(iv,i,j)=flx_s2w(iv,i,j)+(sed_eros_flx_class_by_class(iv)/CELL_SURF(i,j))/MF 
                        ! in kg.m-2 (will be multiplied by dtinv at the end of halfdt)
#ifdef key_MUSTANG_bedload
                  flx_bx(iv,i,j)=flx_bx(iv,i,j)+flx_bxij(iv)/MF ! in kg
                  flx_by(iv,i,j)=flx_by(iv,i,j)+flx_byij(iv)/MF
#ifdef key_MUSTANG_specif_outputs
                  varspecif3Dnv_save(5,iv,i,j)=flx_bx(iv,i,j)
                  varspecif3Dnv_save(6,iv,i,j)=flx_by(iv,i,j)
#endif
#endif

#ifdef key_MUSTANG_debug
                 IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
                    print *,'    iv=',iv
                    print *,'    flx_s2w(iv,i,j)=flx_s2w(iv,i,j)+(sed_eros_flx_class_by_class(iv)/CELL_SURF(i,j))=',flx_s2w(iv,i,j)
#ifdef key_MUSTANG_bedload
                    print *,'    flx_bx(iv,i,j)=flx_bx(iv,i,j)+flx_bxij(iv)=',flx_bx(iv,i,j)
                    print *,'    flx_by(iv,i,j)=flx_by(iv,i,j)+flx_byij(iv)=',flx_by(iv,i,j)
#endif
                  END IF
#endif

                END DO

                dt_ero_max=maxval(dt_ero)

#ifdef key_MUSTANG_specif_outputs
                ! Output: cumulated time (in hours) elapsed in non cohesive regime              
                varspecif2D_save(5,i,j)=varspecif2D_save(5,i,j)+(dt_ero_max/3600.0_rsh)  !tero_non_coh
#endif
                dt1=dt1-dt_ero_max

#ifdef key_MUSTANG_debug
                IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
                  print *,'    TEMPS RESTANT ?'
                  print *,'      > dt_ero_max=',dt_ero_max
                  print *,'      > new dt1=',dt1
                END IF
#endif

                IF (dt1 .GT. 0.0_rsh .AND. dt_ero_max .GT. 0.0_rsh .AND. ksmax .GT. ksmi(i,j)) THEN
                  !print *,' '
                  !print *,' !!!!!! =======> Time is not consumed dt1=',dt1,' ==> CONTINUE EROSION'

#ifdef key_MUSTANG_debug
                  IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
                    print *,'    ==> CONTINUE EROSION !'
                    print *,''
                  END IF
#endif

                  GOTO 2
                END IF

#if ! defined key_nofluxwat_IWS
                phieau_s2w(i,j)=phieau_s2w(i,j)+phieau_ero_ij
#endif
#if ! defined key_noTSdiss_insed 
                flx_s2w(-1:0,i,j)=flx_s2w(-1:0,i,j)+flx_s2w_eroij(-1:0)
                flx_s2w(nvp+1:nv_adv,i,j)=flx_s2w(nvp+1:nv_adv,i,j)+flx_s2w_eroij(nvp+1:nv_adv)
#endif
                !END IF ! IF(tauskin(i,j).GT.toce)

#ifdef key_MUSTANG_splitlayersurf
                !! Splitting surface layers if too thick
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                isplit=0
                DO k=ksmax,ksmax-nlayer_surf_sed+1,-1
                  IF(k > ksmi(i,j)) THEN
                    IF(dzs(k,i,j) > dzsmax(i,j) + 5.0_rsh* dzsmin) isplit=1
                  ENDIF
                ENDDO
                IF(isplit==1 ) then
                 CALL sed_MUSTANG_split_surflayer(i,j,ksmax)
                ENDIF
#else
                ! to avoid increasing the thickness of the surface layer 
                IF(ksmax .LT. ksdmax .AND. ksmax > ksmi(i,j)) THEN
                    IF(dzs(ksmax,i,j) > dzsmax(i,j) + 5.0_rsh* dzsmin) THEN
#ifdef key_MUSTANG_debug
                       IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
                         print *,'    SPLIT SURFACE LAYER BECAUSE dzs > dzsmax '
                         print *,ksmax,'  layers become', ksmax+1, 'layers'
                       END IF
#endif
                       dzs(ksmax+1,i,j)=MIN(dzs(ksmax,i,j)-dzsmax(i,j),dzsmax(i,j))
                       dzs(ksmax,i,j)=dzs(ksmax,i,j)-dzs(ksmax+1,i,j)
                       poro(ksmax+1,i,j)=poro(ksmax,i,j)
                       poro_mud(ksmax+1,i,j)=poro_mud(ksmax,i,j)
                       crel_mud(ksmax+1,i,j)=crel_mud(ksmax,i,j)
                       cv_sed(:,ksmax+1,i,j)=cv_sed(:,ksmax,i,j)
                       c_sedtot(ksmax+1,i,j)=c_sedtot(ksmax,i,j)
                       ksmax=ksmax+1
                    ENDIF
                ENDIF
#endif

              ENDIF  ! no erosion (non cohesive sediment)

            ELSE              ! l_isitcohesive(i,j)==.TRUE.

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !!!              CASE 2: COHESIVE SEDIMENT --> SEDIMENT ERODED AS A MIXTURE               !!!
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              !print *,''
              !print *,'CASE 2: COHESIVE SEDIMENT --> SEDIMENT ERODED AS A MIXTURE'

              k=ksmax
                    
              CALL sed_MUSTANG_comp_tocr_mixsed(k, i, j, xeros, excespowr, toce)

#ifdef key_MUSTANG_specif_outputs
              varspecif3Dnv_save(1,:,i,j)=toce  ! toce_save
#endif
              cvolgrv=0.0_rsh
              DO iv=igrav1,igrav2
                cvolgrv=cvolgrv+cv_sed(iv,k,i,j)/ros(iv)
              ENDDO
              csanmud=0.0_rsh
              DO iv=igrav2+1,nvpc
                csanmud=csanmud+cv_sed(iv,k,i,j)
              ENDDO
              erodab=csanmud*dzs(k,i,j)
              ero=0.0_rsh

              IF(htot(i,j) .GT. h0fond) THEN
                IF(tauskin(i, j) .GT. toce)THEN
                  CALL sed_MUSTANG_comp_eros_flx(tauskin(i, j), toce, excespowr, xeros, sed_eros_flx)
                  ero = sed_eros_flx * fwet(i, j)    

#if defined key_MUSTANG_lateralerosion
                    ! lateral erosion :  wet cell (cellule mouillee)
                    IF (coef_erolat .NE. 0.0_rsh .AND. l_erolat_wet_cell) THEN  
                                                            
                        heaue = HTOT_NEAR_E - htncrit_eros
                        heauw = HTOT_NEAR_W - htncrit_eros
                        heaus = HTOT_NEAR_S - htncrit_eros
                        heaun = HTOT_NEAR_N - htncrit_eros
                        heau_milieu = htot(i, j) - h0fond
                        eroe = max(0.0_rsh, coef_tauskin_lat / 4.0_rsh * V_NEAR_E**2 - toce) &
                            * max(0.0_rsh, heaue - heau_milieu)
                        erow = max(0.0_rsh, coef_tauskin_lat / 4.0_rsh * V_NEAR_W**2 - toce) &
                            * max(0.0_rsh, heauw - heau_milieu)
                        eros = max(0.0_rsh, coef_tauskin_lat / 4.0_rsh * U_NEAR_S**2 - toce) &
                            * max(0.0_rsh, heaus - heau_milieu)
                        eron = max(0.0_rsh, coef_tauskin_lat / 4.0_rsh * U_NEAR_N**2 - toce) &
                            * max(0.0_rsh, heaun - heau_milieu)
                        ero = ero + coef_erolat * (eroe + erow + eros + eron)
                                                                        
                    ENDIF
#endif
                ENDIF
              ELSE
#if defined key_MUSTANG_lateralerosion
                ! lateral erosion :  dry cell 
                IF (coef_erolat .NE. 0.0_rsh) THEN

                    heaue = max(0., HTOT_NEAR_E - htncrit_eros)
                    heauw = max(0., HTOT_NEAR_W - htncrit_eros)
                    heaus = max(0., HTOT_NEAR_S - htncrit_eros)
                    heaun = max(0., HTOT_NEAR_N - htncrit_eros)
                    eroe = coef_erolat * max(0.0_rsh, coef_tauskin_lat / &
                        4.0_rsh * V_NEAR_E**2 - toce) * heaue
                    erow = coef_erolat * max(0.0_rsh, coef_tauskin_lat / &
                        4.0_rsh * V_NEAR_W**2 - toce) * heauw
                    eros = coef_erolat * max(0.0_rsh, coef_tauskin_lat / &
                        4.0_rsh * U_NEAR_S**2 - toce) * heaus
                    eron = coef_erolat * max(0.0_rsh, coef_tauskin_lat / &
                        4.0_rsh * U_NEAR_N**2 - toce) * heaun
                    ero = eroe + erow + eros + eron

                ENDIF
#endif
              ENDIF

              erosi=MF*dt1*ero

              IF(erosi >  0.0_rsh)THEN
              
                 niter_ero_coh=niter_ero_coh+1.0_rsh
                 dzsa=dzs(k,i,j)
                 poroa=poro(k,i,j)
                           
                 IF(erosi.LE.erodab)THEN
       
                          
#ifdef key_MUSTANG_specif_outputs
                      ! Output: cumulated time (in hours) elapsed in cohesive regime
                      ! --> in this case, all sand/mud sediments can be eroded within dt1
                   varspecif2D_save(6,i,j)=varspecif2D_save(6,i,j)+(dt1/3600.0_rsh)  ! tero_coh
#endif
                   ddzs=erosi/csanmud
                
                   DO iv=igrav1,igrav2
                     mass_sed(iv)=cv_sed(iv,k,i,j)*dzsa
                   END DO

                   DO iv=igrav2+1,nv_use ! nv_use=nvp (or nvpc if key_Pconstitonly_insed)
                     dflx_s2w=ddzs*cv_sed(iv,k,i,j)
                     !flx_s2w(iv,i,j)=flx_s2w(iv,i,j)+dflx_s2w/MF 
                     flx_s2w_eroij(iv)=flx_s2w_eroij(iv)+dflx_s2w/MF
                     !cv_sed(iv,k,i,j)=(cv_sed(iv,k,i,j)*dzsa-dflx_s2w)*dzsisu
                     mass_sed(iv)=(cv_sed(iv,k,i,j)*dzsa)-dflx_s2w
                   ENDDO
             
                   mass_tot=0.0_rsh
                   DO iv=1,nvpc
                     mass_tot=mass_tot+mass_sed(iv)
                   END DO

                   IF (mass_tot .GT. 0.0_rsh) THEN
 
                     DO iv=1,imud2
                      frac_sed(iv)=mass_sed(iv)/mass_tot
                     END DO

                     CALL MUSTANGV2_comp_poro_mixsed(frac_sed, poro_mud(k,i,j), &
                                crel_mud(k,i,j), poro(k,i,j)) 
                             ! poro_mud is unchanged here

                     !IF (crel_mud(k,i,j) .GT. 1500.0_rsh) THEN
                     !  print *,'in sed_erosion_mixsed 2'
                     !  print *,' > crel_mud(k,i,j) = ',crel_mud(k,i,j)
                     !  print *,' > sommud = ',sommud
                     !  print *,' > cvolgrvsan = ',cvolgrvsan
                     !  print *,' > cv_sed(:,k,i,j) = ',cv_sed(:,k,i,j)
                     !END IF


                     dzs(k,i,j)=mass_tot/((1.0_rsh-poro(k,i,j))*ros(1))

#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
                     ! dissolved variable in pore waters and water fluxes at the interface
                     !--------------------------------------------------------------------
                     !  ==> flx_s2w and cv_sed and phieau_s2w
                     ! layer modified by erosion  ==> code 3
                     ! dissolved concentrations in this layer unchanged by erosion 
                       porowater_new=poro(ksmax,i,j)*dzs(ksmax,i,j)
                       porowater1=dzsa*poroa
                       CALL MUSTANGV2_eval_dissvar_IWSflux(i,j,ksmax,3,  &
                                    phieau_ero_ij=phieau_ero_ij,           &
                                    flx_s2w_eroij=flx_s2w_eroij,           &
                                    porowater1=porowater1,                 &
                                    cv_sed1=cv_sed(-1:nv_adv,ksmax,i,j),   &
                                    porowater_new=porowater_new)
#endif

                     !dzsmin=(1.0_rsh-coeff_dzsmin)*dzsminuni + coeff_dzsmin*SUM( (cv_sed(1:nvpc,k,i,j)/c_sedtot(k,i,j))*diam_sed(1:nvpc) )
                     !dzsmin=dzsminuni + coeff_dzsmin*SUM( (cv_sed(1:isand2,k,i,j)/c_sedtot(k,i,j))*diam_sed(1:isand2) ) ??
                     !dzsmin=(1.0_rsh-coeff_dzsmin)*dzsminuni + coeff_dzsmin*SUM( frac_sed(1:nvpc)*diam_sed(1:nvpc) )
                     dzsmin=dzsminvar(frac_sed)


                     IF (dzs(k,i,j) .GT. dzsmin) THEN

                       DO iv=1,nvp
                         cv_sed(iv,k,i,j)=mass_sed(iv)/dzs(k,i,j)
                       END DO
                       c_sedtot(k,i,j)=mass_tot/dzs(k,i,j)

                     ELSE  ! dzs<=dzsmin

                      CALL MUSTANGV2_manage_small_mass_in_ksmax(i,j,ksmax, &
#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
                                  phieau_ero_ij,flx_s2w_eroij,   &
#endif
                                                        mass_sed)
                      ! dissolved concentrations and fluxes are updated during fusion in routine MUSTANGV2_manage_small_mass_in_ksmax

                     END IF

                   ELSE  ! mass_tot<=0
 
#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
                      ! dissolved variable in pore waters and water fluxes at the interface
                      !--------------------------------------------------------------------
                      !  ==> flx_s2w and phieau_s2w
                      ! layer elimination ==> code 1
                     porowater1=dzsa*poroa
                     CALL MUSTANGV2_eval_dissvar_IWSflux(i,j,ksmax,1,  &
                            phieau_ero_ij=phieau_ero_ij,           &
                            flx_s2w_eroij=flx_s2w_eroij,           &
                            porowater1=porowater1,                   &
                            cv_sed1=cv_sed(-1:nv_adv,ksmax,i,j))
#endif
                     dzs(ksmax,i,j)=0.0_rsh
                     cv_sed(:,ksmax,i,j)=0.0_rsh
                     c_sedtot(ksmax,i,j)=0.0_rsh
                     poro(ksmax,i,j)=0.0_rsh
                     poro_mud(ksmax,i,j)=0.0_rsh
                     crel_mud(ksmax,i,j)=0.0_rsh

                     ksmax=ksmax-1

                   END IF

#if defined key_MUSTANG_lateralerosion
                   ! memorisation of lateral erosion for dry cell
                   IF(htot(i,j).LE.h0fond .AND. coef_erolat .NE. 0.0_rsh) THEN
                     DO iv=-1,nv_adv
                       flx_s2w_corip1(iv,i,j)=flx_s2w_corip1(iv,i,j)+  &
                                              flx_s2w_eroij(iv)*eroe/ero*CELL_SURF(i,j)/SURF_NEAR_E
                       flx_s2w_corim1(iv,i,j)=flx_s2w_corim1(iv,i,j)+  &
                                              flx_s2w_eroij(iv)*erow/ero*CELL_SURF(i,j)/SURF_NEAR_W
                       flx_s2w_corjm1(iv,i,j)=flx_s2w_corjm1(iv,i,j)+  &
                                              flx_s2w_eroij(iv)*eros/ero*CELL_SURF(i,j)/SURF_NEAR_S
                       flx_s2w_corjp1(iv,i,j)=flx_s2w_corjp1(iv,i,j)+  &
                                              flx_s2w_eroij(iv)*eron/ero*CELL_SURF(i,j)/SURF_NEAR_N
                     ENDDO
#if ! defined key_nofluxwat_IWS
                     phieau_s2w_corip1(i,j)=phieau_s2w_corip1(i,j)+phieau_ero_ij*eroe/ero
                     phieau_s2w_corim1(i,j)=phieau_s2w_corim1(i,j)+phieau_ero_ij*erow/ero
                     phieau_s2w_corjm1(i,j)=phieau_s2w_corjm1(i,j)+phieau_ero_ij*eros/ero
                     phieau_s2w_corjp1(i,j)=phieau_s2w_corjp1(i,j)+phieau_ero_ij*eron/ero
#endif
                   ENDIF
#endif

                   phieau_s2w(i,j)=phieau_s2w(i,j)+phieau_ero_ij
                   flx_s2w(:,i,j)=flx_s2w(:,i,j)+flx_s2w_eroij(:)

                 ELSE IF(csanmud /= 0.0_rsh) THEN  !erosi/erodab
                  ! ELSE ! --> suffisant car probablement inutile etant donne que l_isit_cohesive = TRUE

#ifdef key_MUSTANG_specif_outputs
                    ! Output: cumulated time (in hours) elapsed in cohesive regime
                    varspecif2D_save(6,i,j)=varspecif2D_save(6,i,j)+(dt1*(erodab/erosi)/3600.0_rsh)   ! tero_coh
#endif                    
                    ! correction PLH Juillet 2015
                    !volerod=erodab*fwet(i,j)/csanmud
                    !volerod=erodab/csanmud
                
                    DO iv=igrav2+1,nv_use ! nv_use=nvp (or nvpc if key_Pconstitonly_insed)
                       !flx_s2w(iv,i,j)=flx_s2w(iv,i,j)+volerod*cv_sed(iv,k,i,j)/MF
                       !flx_s2w(iv,i,j)=flx_s2w(iv,i,j)+dzs(k,i,j)*cv_sed(iv,k,i,j)/MF
                       flx_s2w_eroij(iv)=flx_s2w_eroij(iv)+dzs(k,i,j)*cv_sed(iv,k,i,j)/MF
                       cv_sed(iv,k,i,j)=0.0_rsh
                    ENDDO

                
                    !IF(cvolgrv.LE.0.0_rsh)THEN
                    IF(cvolgrv.LE.epsilon_MUSTANG)THEN ! Attention, a la base cvolgrv.LE.0.0_rsh, conservativite ?
                    
#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
                       ! dissolved variable in pore waters and water fluxes at the interface
                       !--------------------------------------------------------------------
                       !  ==> flx_s2w and phieau_s2w
                       ! layer elimination ==> code 1
                       porowater1=dzsa*poroa
                       CALL MUSTANGV2_eval_dissvar_IWSflux(i,j,ksmax,1, &
                            phieau_ero_ij=phieau_ero_ij,           &
                            flx_s2w_eroij=flx_s2w_eroij,           &
                            porowater1=porowater1, &
                            cv_sed1=cv_sed(-1:nv_adv,ksmax,i,j))
#endif

                       ! reset  dzs and poro
                       dzs(ksmax,i,j)=0.0_rsh
                       poro(ksmax,i,j)=0.0_rsh
                       poro_mud(ksmax,i,j)=0.0_rsh
                       cv_sed(:,ksmax,i,j)=0.0_rsh
                       c_sedtot(ksmax,i,j)=0.0_rsh
                       crel_mud(ksmax,i,j)=0.0_rsh
                  
                       ksmax=ksmax-1
                  
                    ELSE

                       mass_sed(:)=0.0_rsh
                       mass_tot=0.0_rsh
                       DO iv=igrav1,igrav2
                         mass_sed(iv)=cv_sed(iv,k,i,j)*dzsa
                         mass_tot=mass_tot+mass_sed(iv)
                       END DO

                       frac_sed(:)=0.0_rsh
                       DO iv=igrav1,igrav2
                         frac_sed(iv)=mass_sed(iv)/mass_tot
                       END DO

                       
                       !  crel_mud(k,i,j) = 0.0_rsh ????
                       CALL MUSTANGV2_comp_poro_mixsed(frac_sed, poro_mud(k,i,j), &
                                crel_mud(k,i,j),poro(k,i,j)) 


                       !IF (crel_mud(k,i,j) .GT. 1500.0_rsh) THEN
                       !  print *,'in sed_erosion_mixsed 2'
                       !  print *,' > crel_mud(k,i,j) = ',crel_mud(k,i,j)
                       !  print *,' > sommud = ',sommud
                       !  print *,' > cvolgrvsan = ',cvolgrvsan
                       !  print *,' > cv_sed(:,k,i,j) = ',cv_sed(:,k,i,j)
                       !END IF

                       dzs(k,i,j)=mass_tot/((1.0_rsh-poro(k,i,j))*ros(1))

#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
                        ! dissolved variable in pore waters and water fluxes at the interface
                        !--------------------------------------------------------------------
                        !  ==> flx_s2w and cv_sed and phieau_s2w
                        ! dissolved concentrations in this layer unchanged by erosion  
                        porowater_new=poro(ksmax,i,j)*dzs(ksmax,i,j)
                        porowater1=dzsa*poroa
                        CALL MUSTANGV2_eval_dissvar_IWSflux(i,j,ksmax,3,  &
                                    phieau_ero_ij=phieau_ero_ij,           &
                                    flx_s2w_eroij=flx_s2w_eroij,           &
                                    porowater1=porowater1,                &
                                    cv_sed1=cv_sed(-1:nv_adv,ksmax,i,j),  &
                                    porowater_new=porowater_new)
#endif

                       ! BM: In case of very small remaining masses of gravel(s), the new dzs(k,i,j) can be 0
                       ! To prevent this, a test is done on dzs to determine if the layer characteristics
                       ! can be updated or if these small sediment quantities have to be put in the underlying
                       ! layer. The criterion is dzsmin here

                       !dzsmin=(1.0_rsh-coeff_dzsmin)*dzsminuni + coeff_dzsmin*SUM( (cv_sed(igrav1:igrav2,k,i,j)/c_sedtot(k,i,j))*diam_sed(igrav1:igrav2) ) ! NON 
                       !dzsmin=(1.0_rsh-coeff_dzsmin)*dzsminuni + coeff_dzsmin* &
                       !        SUM(frac_sed(igrav1:igrav2)*diam_sed(igrav1:igrav2))                    
                       dzsmin=dzsminvar(frac_sed)

                       IF (dzs(k,i,j) .GT. dzsmin) THEN

                         c_sedtot(k,i,j)=0.0_rsh
                         DO iv=igrav1,igrav2
                           cv_sed(iv,k,i,j)=mass_sed(iv)/dzs(k,i,j)
                           c_sedtot(k,i,j)=c_sedtot(k,i,j)+cv_sed(iv,k,i,j)
                         END DO
                         

                       ELSE   ! dzs<=dzsmin

                        CALL MUSTANGV2_manage_small_mass_in_ksmax(i,j,ksmax,   &
#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
                                       phieau_ero_ij,flx_s2w_eroij,  &
#endif
                                                   mass_sed)
                         ! dissolved concentrations IWS fluxes updated in MUSTANGV2_manage_small_mass_in_ksmax
                   
                       ENDIF
                   ENDIF  ! en loop on cvolgrv
              
#if defined key_MUSTANG_lateralerosion
                   ! memorisation of lateral erosion for dry cell
                   IF(htot(i,j).LE.h0fond .AND. coef_erolat .NE. 0.0_rsh) THEN
                      DO iv=-1,nv_adv
                        flx_s2w_corip1(iv,i,j)=flx_s2w_corip1(iv,i,j)+   &
                                               flx_s2w_eroij(iv)*eroe/ero*CELL_SURF(i,j)/SURF_NEAR_E
                        flx_s2w_corim1(iv,i,j)=flx_s2w_corim1(iv,i,j)+   &
                                               flx_s2w_eroij(iv)*erow/ero*CELL_SURF(i,j)/SURF_NEAR_W
                        flx_s2w_corjm1(iv,i,j)=flx_s2w_corjm1(iv,i,j)+   &
                                               flx_s2w_eroij(iv)*eros/ero*CELL_SURF(i,j)/SURF_NEAR_S
                        flx_s2w_corjp1(iv,i,j)=flx_s2w_corjp1(iv,i,j)+   &
                                               flx_s2w_eroij(iv)*eron/ero*CELL_SURF(i,j)/SURF_NEAR_N
                      ENDDO
#if ! defined key_nofluxwat_IWS
                      phieau_s2w_corip1(i,j)=phieau_s2w_corip1(i,j)+phieau_ero_ij*eroe/ero  ! *SURF_NEAR_E/CELL_SURF(i,j)
                      phieau_s2w_corim1(i,j)=phieau_s2w_corim1(i,j)+phieau_ero_ij*erow/ero  !*SURF_NEAR_W/CELL_SURF(i,j)
                      phieau_s2w_corjm1(i,j)=phieau_s2w_corjm1(i,j)+phieau_ero_ij*eros/ero  !*SURF_NEAR_S/CELL_SURF(i,j)
                      phieau_s2w_corjp1(i,j)=phieau_s2w_corjp1(i,j)+phieau_ero_ij*eron/ero  !*SURF_NEAR_N/CELL_SURF(i,j)
#endif
                   ENDIF
#endif                

                   IF (ksmax .GT. ksmi(i,j)) THEN
                    dt1=dt1*(1.0_rsh-erodab/erosi) !  erosi=MF*dt1*ero
                    !print *,' !!!!!! =======> Il reste du temps dt1=dt1-erodab/ero=',dt1,' ==> CONTINUE EROSION'
                    GOTO 2
                  ENDIF 
                  phieau_s2w(i,j)=phieau_s2w(i,j)+phieau_ero_ij
                  flx_s2w(:,i,j)=flx_s2w(:,i,j)+flx_s2w_eroij(:)


              ENDIF   ! fin test erosi > erodab
          
              ! in order to avoid increasing thickness of the surface layer
              ! if after manage_small dzs(ksmax) > dzsmax+2*dzsmin ==> split surface layer
              IF(ksmax .LT. ksdmax .AND. dzs(ksmax,i,j) > dzsmax(i,j) + 2* dzsmin) THEN
#ifdef key_MUSTANG_debug
                  IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
                    print *,'    SPLIT SURFACE LAYER BECAUSE dzs > dzsmax '
                    print *,ksmax,'  layers become', ksmax+1, 'layers'
                  END IF
#endif
                dzs(ksmax+1,i,j)=MIN(dzs(ksmax,i,j)-dzsmax(i,j),dzsmax(i,j))   
                dzs(ksmax,i,j)=dzs(ksmax,i,j)-dzs(ksmax+1,i,j)
                poro(ksmax+1,i,j)=poro(ksmax,i,j)
                poro_mud(ksmax+1,i,j)=poro_mud(ksmax,i,j)
                cv_sed(:,ksmax+1,i,j)=cv_sed(:,ksmax,i,j)
                crel_mud(ksmax+1,i,j)=crel_mud(ksmax,i,j)
                c_sedtot(ksmax+1,i,j)=c_sedtot(ksmax,i,j)
                ksmax=ksmax+1
              ENDIF
          
             ENDIF    ! fin test erosi=0
            ENDIF     ! fin test l_isit_cohesive
            
        ENDIF       ! fin test ksmax > ksmi(i,j) 

          ! to find an erosion flux in .../m2/s (for particulates only):
        DO iv=1,nv_use ! nv_use=nvp (or nvpc if key_Pconstitonly_insed)
            flx_s2w(iv,i,j)=flx_s2w(iv,i,j)*dtinv
        ENDDO
        
          !==============================================
          
#if ! defined key_noTSdiss_insed
          ! for dissolved variable in pore waters
          ! add water-to-sediment flux explicitly treated with the sediment-to-water flux also explicitly treated
#if ! defined key_Pconstitonly_insed
          !   for dissolved variables, we add the diffusive flux to the interface
          !    and the flux which results from the expulsion of interstitial water by consolidation:
        DO iv=nvp+1,nv_adv
            flx_s2w(iv,i,j)=(flx_s2w(iv,i,j)-flx_w2s(iv,i,j)+fluconsol(iv,i,j)-fludif(iv,i,j))*dtinv
        ENDDO
#endif

        DO iv=-1,0
            flx_s2w(iv,i,j)=(flx_s2w(iv,i,j)-flx_w2s(iv,i,j)+fluconsol(iv,i,j)-fludif(iv,i,j))*dtinv
        ENDDO
#endif

        ksma(i,j)=ksmax

          ! no correction of water column height (could be !)
          
#ifdef key_MUSTANG_specif_outputs
          ! Stats on non-coh/coh erosion iterations during halfdt
        varspecif2D_save(9,i,j)=niter_ero_noncoh+niter_ero_coh           ! niter_ero
        varspecif2D_save(7,i,j)=niter_ero_noncoh/(varspecif2D_save(9,i,j)+epsilon_MUSTANG) ! pct_iter_noncoh
        varspecif2D_save(8,i,j)=niter_ero_coh/(varspecif2D_save(9,i,j)+epsilon_MUSTANG)    !pct_iter_coh
#endif

#ifdef key_MUSTANG_debug
               IF (l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
                 print *,'  > fin erosion'
                 print *,'  t=',CURRENT_TIME,' ksmax=',ksmax
                 print *,'  dzs(ksmax-3:ksmax,i,j)=',dzs(ksmax-3:ksmax,i,j)
                 print *,'  cv_sed(:,ksmax-3,i,j)=',cv_sed(:,ksmax-3,i,j)
                 print *,'  cv_sed(:,ksmax-2,i,j)=',cv_sed(:,ksmax-2,i,j)
                 print *,'  cv_sed(:,ksmax-1,i,j)=',cv_sed(:,ksmax-1,i,j)
                 print *,'  cv_sed(:,ksmax,i,j)=',cv_sed(:,ksmax,i,j)
               END IF
#endif 

     END DO
   END DO

   ! version V2
  END SUBROUTINE sed_MUSTANG_erosion


#else
!  version V1
   !!==============================================================================
  SUBROUTINE sed_MUSTANG_erosion(ifirst, ilast, jfirst, jlast, dtinv,  &
#if defined key_MUSTANG_lateralerosion 
                                 BAROTROP_VELOCITY_U, BAROTROP_VELOCITY_V, &
#endif
                                  dt_true) 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_erosion  version V1 ***
   !&E
   !&E ** Purpose : computes erosion fluxes
   !&E
   !&E ** Description : 
   !&E       arguments IN :
   !&E          loops  :ifirst,ilast,jfirst,jlast
   !&E          dtinv, dt_true : 1/dt  and dt
   !&E          CELL_SURF : cells surface
   !&E          u,v :   BAROTROP_VELOCITY_U, BAROTROP_VELOCITY_V
   !&E      
   !&E     variables OUT : 
   !&E          flx_s2w : concentration flux sediment to water
   !&E          phieau_s2w : water flux  sediment to water
   !&E          the sediment has been remodelled
   !&E          cv_sed, dzs, ksma, poro have changed
   !&E
   !&E
   !&E ** Called by :  MUSTANG_update
  !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used

   INTEGER, INTENT(IN)                        :: ifirst, ilast, jfirst, jlast
   REAL(KIND=rsh),INTENT(IN)                  :: dtinv
   REAL(KIND=rlg),INTENT(IN)                  :: dt_true  ! =halfdt in MARS
#if defined key_MUSTANG_lateralerosion 
   REAL(KIND=rsh),DIMENSION(ARRAY_VELOCITY_U),INTENT(IN)   :: BAROTROP_VELOCITY_U
   REAL(KIND=rsh),DIMENSION(ARRAY_VELOCITY_V),INTENT(IN)   :: BAROTROP_VELOCITY_V 
#endif

   !! * Local declarations
   INTEGER        ::  i,j,k,iv,ksmax,ksup,ksmaxa,isplit
   LOGICAL        ::  fusion
   REAL(KIND=rsh) ::  dt1,toce,cvolgrv,csanmud,phieau_ero_ij,                  &
                      erodab,ero,erosi,ddzs,dzsa,dzsisu,dzsam1,dflx_s2w,       &
                      cvolp,volerod,poroa,poroam1,dflusve,      &
                      xeros, sed_eros_flx,excespowr,dzpoi,      &
                      heauw,heaue,heaun,heaus,heau_milieu,eroe,erow,eros,eron
! because of lateral erosion, it is necessary to keep in memory the accumulated eroded fluxes
   REAL(KIND=rsh),DIMENSION(-1:nv_adv)  ::  flx_s2w_eroij

   !!---------------------------------------------------------------------------
   !! * Executable part

      DO j=jfirst,jlast
        DO i=ifirst,ilast

          flx_s2w(-1,i,j)=0.0_rsh
          flx_s2w( 0,i,j)=0.0_rsh
          DO iv=1,nv_adv
            flx_s2w(iv,i,j)=0.0_rsh
          ENDDO
          phieau_ero_ij=0.0_rsh
          flx_s2w_eroij(:)=0.0_rsh
        
          ksmax=ksma(i,j)
          ero=0.0_rsh

          IF(ksmax.GE.ksmi(i,j))THEN

            dt1=REAL(dt_true,rsh)
            k=ksmax
                       
        2   CONTINUE

            CALL sed_MUSTANG_comp_tocr_mixsed(k, i, j, xeros, excespowr, toce)
            
#ifdef key_MUSTANG_specif_outputs
            varspecif3Dnv_save(1,:,i,j)=toce  ! toce_save
#endif
            cvolgrv=0.0_rsh
            DO iv=igrav1,igrav2
              cvolgrv=cvolgrv+cv_sed(iv,k,i,j)/ros(iv)
            ENDDO
            csanmud=0.0_rsh
            DO iv=igrav2+1,nvpc
              csanmud=csanmud+cv_sed(iv,k,i,j)
            ENDDO
            erodab=csanmud*dzs(k,i,j)
            ero=0.0_rsh

            IF(htot(i,j) .GT. h0fond) THEN
              IF(tauskin(i,j) .GT. toce)THEN
                CALL sed_MUSTANG_comp_eros_flx(tauskin(i,j), toce, excespowr, xeros, sed_eros_flx)
                ero=sed_eros_flx*fwet(i,j)    

#if defined key_MUSTANG_lateralerosion
                ! lateral erosion :  wet cell (cellule mouillee)
                IF (coef_erolat .NE. 0.0_rsh .AND. l_erolat_wet_cell) THEN                                                              
                    heaue = HTOT_NEAR_E - htncrit_eros
                    heauw = HTOT_NEAR_W - htncrit_eros
                    heaus = HTOT_NEAR_S - htncrit_eros
                    heaun = HTOT_NEAR_N - htncrit_eros
                    heau_milieu = htot(i,j) - h0fond
                    eroe = max(0.0_rsh, coef_tauskin_lat / 4.0_rsh * V_NEAR_E**2 - toce) &
                        * max(0.0_rsh, heaue - heau_milieu)
                    erow = max(0.0_rsh, coef_tauskin_lat / 4.0_rsh * V_NEAR_W**2 - toce) &
                        * max(0.0_rsh, heauw - heau_milieu)
                    eros = max(0.0_rsh, coef_tauskin_lat / 4.0_rsh * U_NEAR_S**2 - toce) &
                        * max(0.0_rsh, heaus - heau_milieu)
                    eron = max(0.0_rsh, coef_tauskin_lat / 4.0_rsh * U_NEAR_N**2 - toce) &
                        * max(0.0_rsh, heaun - heau_milieu)
                    ero = ero + coef_erolat * (eroe + erow + eros + eron)                                                                       
                ENDIF
#endif
              ENDIF
            ELSE
#if defined key_MUSTANG_lateralerosion
              ! lateral erosion :  dry cell 
                IF (coef_erolat .NE. 0.0_rsh) THEN
                    !write(*,*)'l_erolat_dry_cell'
                    heaue = max(0., HTOT_NEAR_E - htncrit_eros)
                    heauw = max(0., HTOT_NEAR_W - htncrit_eros)
                    heaus = max(0., HTOT_NEAR_S - htncrit_eros)
                    heaun = max(0., HTOT_NEAR_N - htncrit_eros)
                    eroe = coef_erolat * max(0.0_rsh, coef_tauskin_lat / &
                        4.0_rsh * V_NEAR_E**2 - toce) * heaue
                    erow = coef_erolat * max(0.0_rsh, coef_tauskin_lat / &
                        4.0_rsh * V_NEAR_W**2 - toce) * heauw
                    eros = coef_erolat * max(0.0_rsh, coef_tauskin_lat / &
                        4.0_rsh * U_NEAR_S**2 - toce) * heaus
                    eron = coef_erolat * max(0.0_rsh, coef_tauskin_lat / &
                        4.0_rsh * U_NEAR_N**2 - toce) * heaun
                    ero = eroe + erow + eros + eron
                    !write(*,*)'l_erolat_dry_cell',ero
                ENDIF
#endif
            ENDIF

            erosi = MF * dt1 * ero

            IF(erosi >  0.0_rsh)THEN
              IF(erosi.LE.erodab)THEN
                ddzs=erosi/csanmud
                dzsa=dzs(k,i,j)
                dzs(k,i,j)=MAX(dzsa-ddzs,dzsa*cvolgrv/cvolmaxsort)
                IF(dzs(k,i,j).LE.0.0_rsh)THEN
                  !WRITE(*,*)'en kij:',k,i,j,dzs(k,i,j),dzsa,ddzs,csanmud,cvolgrv
                  dzs(k,i,j)=epsilon_MUSTANG
                ENDIF
                dzsisu=1.0_rsh/dzs(k,i,j)
                DO iv=igrav1,igrav2
                  cv_sed(iv,k,i,j)=cv_sed(iv,k,i,j)*dzsa*dzsisu
                ENDDO
                DO iv=igrav2+1,nv_use ! nv_use=nvp (or nvpc if key_Pconstitonly_insed)
                  dflx_s2w=ddzs*cv_sed(iv,k,i,j)
                  !flx_s2w(iv,i,j)=flx_s2w(iv,i,j)+dflx_s2w/MF 
                  flx_s2w_eroij(iv)=flx_s2w_eroij(iv)+dflx_s2w/MF
                  cv_sed(iv,k,i,j)=(cv_sed(iv,k,i,j)*dzsa-dflx_s2w)*dzsisu
                ENDDO
                c_sedtot(k,i,j)=0.0_rsh
                cvolp=0.0_rsh
                DO iv=1,nvpc
                  c_sedtot(k,i,j)=c_sedtot(k,i,j)+cv_sed(iv,k,i,j)
                  cvolp=cvolp+cv_sed(iv,k,i,j)/ros(iv)
                ENDDO
                poroa=poro(k,i,j)
                poro(k,i,j)=1.0_rsh-cvolp

#if ! defined key_noTSdiss_insed
                ! dissolved variable in pore waters
#if ! defined key_Pconstitonly_insed
                DO iv=nvp+1,nv_adv
                  dflusve=ddzs*poroa*cv_sed(iv,k,i,j)
                  flx_s2w_eroij(iv)=flx_s2w_eroij(iv)+dflusve/MF
                  ! Fev2017 BT : qd il y a erosion, pourquoi conc dissouste changerait elle dans les eaux intertertitielles ?
                  ! cv_sed(iv,k,i,j)=(cv_sed(iv,k,i,j)*dzsa*poroa-dflusve)*dzsisu/poro(k,i,j) 
                ENDDO
#endif

                ! Temperature and salinity
                DO iv=-1,0
                  dflusve=ddzs*poroa*cv_sed(iv,k,i,j)
                  flx_s2w_eroij(iv)=flx_s2w_eroij(iv)+dflusve/MF
                  ! Fev2017 BT : qd il y a erosion, pourquoi conc dissouste changerait elle dans les eaux intertertitielles ?
                  ! cv_sed(iv,k,i,j)=(cv_sed(iv,k,i,j)*dzsa*poroa-dflusve)*dzsisu/poro(k,i,j)
                ENDDO
#endif

#if ! defined key_nofluxwat_IWS
                ! flux d eau du a l erosion
                phieau_ero_ij=phieau_ero_ij+REAL(ddzs*poroa*CELL_SURF(i,j),rlg)/MF
#endif

                ! memorisation of lateral erosion for dry cell
                IF(htot(i,j).LE.h0fond .AND. coef_erolat .NE. 0.0_rsh) THEN
                   DO iv=-1,nv_adv
                     flx_s2w_corip1(iv,i,j)=flx_s2w_corip1(iv,i,j)+    &
                                            flx_s2w_eroij(iv)*eroe/ero*CELL_SURF(i,j)/SURF_NEAR_E
                     flx_s2w_corim1(iv,i,j)=flx_s2w_corim1(iv,i,j)+    &
                                            flx_s2w_eroij(iv)*erow/ero*CELL_SURF(i,j)/SURF_NEAR_W
                     flx_s2w_corjm1(iv,i,j)=flx_s2w_corjm1(iv,i,j)+    &
                                            flx_s2w_eroij(iv)*eros/ero*CELL_SURF(i,j)/SURF_NEAR_S
                     flx_s2w_corjp1(iv,i,j)=flx_s2w_corjp1(iv,i,j)+    &
                                            flx_s2w_eroij(iv)*eron/ero*CELL_SURF(i,j)/SURF_NEAR_N
                   ENDDO
#if ! defined key_nofluxwat_IWS
                   phieau_s2w_corip1(i,j)=phieau_s2w_corip1(i,j)+phieau_ero_ij*eroe/ero
                   phieau_s2w_corim1(i,j)=phieau_s2w_corim1(i,j)+phieau_ero_ij*erow/ero
                   phieau_s2w_corjm1(i,j)=phieau_s2w_corjm1(i,j)+phieau_ero_ij*eros/ero
                   phieau_s2w_corjp1(i,j)=phieau_s2w_corjp1(i,j)+phieau_ero_ij*eron/ero
#endif
                ENDIF
                phieau_s2w(i,j)=phieau_s2w(i,j)+phieau_ero_ij
                flx_s2w(:,i,j)=flx_s2w(:,i,j)+flx_s2w_eroij(:)

              ELSE IF(csanmud /= 0.0_rsh) THEN  !erosi/erodab
                ! correction PLH Juillet 2015
                !volerod=erodab*fwet(i,j)/csanmud
                volerod=erodab/csanmud
                DO iv=igrav2+1,nv_use ! nv_use=nvp (or nvpc if key_Pconstitonly_insed)
                  !flx_s2w(iv,i,j)=flx_s2w(iv,i,j)+volerod*cv_sed(iv,k,i,j)/MF
                  flx_s2w_eroij(iv)=flx_s2w_eroij(iv)+volerod*cv_sed(iv,k,i,j)/MF
                  cv_sed(iv,k,i,j)=0.0_rsh
                ENDDO
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
                ! dissolved variable in pore waters
                DO iv=nvp+1,nv_adv
                  flx_s2w_eroij(iv)=flx_s2w_eroij(iv)+volerod*poro(k,i,j)*cv_sed(iv,k,i,j)/MF
                  cv_sed(iv,k,i,j)=0.0_rsh
                ENDDO
#endif
                ! Temperature and salinity
                DO iv=-1,0
                  flx_s2w_eroij(iv)=flx_s2w_eroij(iv)+volerod*poro(k,i,j)*cv_sed(iv,k,i,j)/MF
                  cv_sed(iv,k,i,j)=0.0_rsh
                ENDDO
#endif
#if ! defined key_nofluxwat_IWS
               ! water flux due to erosion
                phieau_ero_ij=phieau_ero_ij+ REAL(volerod*poro(k,i,j)*CELL_SURF(i,j),rlg)/MF
#endif
                IF(cvolgrv.LE.0.0_rsh)THEN
                  ! reset  dzs and poro
                  dzs(ksmax,i,j)=0.0_rsh
                  poro(ksmax,i,j)=0.0_rsh
                  ksmax=ksmax-1
                ELSE
                  dzsa=dzs(k,i,j)
                  dzs(k,i,j)=dzsa*cvolgrv/cvolmaxsort
                  !??? dzs(k,i,j)=dzsa*cvolgrv/cvolmaxmel
                  c_sedtot(k,i,j)=0.0_rsh
                  DO iv=igrav1,igrav2
                      cv_sed(iv,k,i,j)=cv_sed(iv,k,i,j)*dzsa/dzs(k,i,j)
                      c_sedtot(k,i,j)=c_sedtot(k,i,j)+cv_sed(iv,k,i,j)
                  ENDDO
                  poroa=poro(k,i,j)
                  poro(k,i,j)=1.0_rsh-cvolmaxsort
                  !??? poro(k,i,j)=1.0_rsh-cvolmaxmel
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
                  ! dissolved variable in pore waters
                  DO iv=nvp+1,nv_adv
                    !flx_s2w(iv,i,j)=flx_s2w(iv,i,j)+cv_sed(iv,k,i,j)*(dzsa*poroa   &
                    flx_s2w_eroij(iv)=flx_s2w_eroij(iv)+cv_sed(iv,k,i,j)*(dzsa*poroa   &
                                   -dzs(k,i,j)*poro(k,i,j))
                  ENDDO
#endif

                  ! Temperature and salinity
                  DO iv=-1,0
                    !flx_s2w(iv,i,j)=flx_s2w(iv,i,j)+cv_sed(iv,k,i,j)*(dzsa*poroa   &
                    flx_s2w_eroij(iv)=flx_s2w_eroij(iv)+cv_sed(iv,k,i,j)*(dzsa*poroa   &
                                   -dzs(k,i,j)*poro(k,i,j))
                  ENDDO
#endif
#if ! defined key_nofluxwat_IWS
               ! water flux due to erosion
                  !phieau_s2w(i,j)=phieau_s2w(i,j)+                         &
                  phieau_ero_ij=phieau_ero_ij + &
                                 REAL((dzsa*poroa-dzs(k,i,j)*poro(k,i,j))*CELL_SURF(i,j),rlg)
#endif                                 
                ENDIF  ! end loop on cvolgrv
              
#if defined key_MUSTANG_lateralerosion
                ! memorisation of lateral erosion for dry cell
                IF(htot(i,j).LE.h0fond .AND. coef_erolat .NE. 0.0_rsh) THEN
                   DO iv=-1,nv_adv
                     flx_s2w_corip1(iv,i,j)=flx_s2w_corip1(iv,i,j)+   &
                                            flx_s2w_eroij(iv)*eroe/ero*CELL_SURF(i,j)/SURF_NEAR_E
                     flx_s2w_corim1(iv,i,j)=flx_s2w_corim1(iv,i,j)+   &
                                            flx_s2w_eroij(iv)*erow/ero*CELL_SURF(i,j)/SURF_NEAR_W
                     flx_s2w_corjm1(iv,i,j)=flx_s2w_corjm1(iv,i,j)+   &
                                            flx_s2w_eroij(iv)*eros/ero*CELL_SURF(i,j)/SURF_NEAR_S
                     flx_s2w_corjp1(iv,i,j)=flx_s2w_corjp1(iv,i,j)+   &
                                            flx_s2w_eroij(iv)*eron/ero*CELL_SURF(i,j)/SURF_NEAR_N
                   ENDDO
#if ! defined key_nofluxwat_IWS
                   phieau_s2w_corip1(i,j)=phieau_s2w_corip1(i,j)+phieau_ero_ij*eroe/ero  ! *SURF_NEAR_E/CELL_SURF(i,j)
                   phieau_s2w_corim1(i,j)=phieau_s2w_corim1(i,j)+phieau_ero_ij*erow/ero  !*SURF_NEAR_W/CELL_SURF(i,j)
                   phieau_s2w_corjm1(i,j)=phieau_s2w_corjm1(i,j)+phieau_ero_ij*eros/ero  !*SURF_NEAR_S/CELL_SURF(i,j)
                   phieau_s2w_corjp1(i,j)=phieau_s2w_corjp1(i,j)+phieau_ero_ij*eron/ero  !*SURF_NEAR_N/CELL_SURF(i,j)
#endif
                ENDIF
#endif


                k=k-1
                IF(k.GE.ksmi(i,j))THEN
                  dt1=dt1-erodab/(MF*ero)
                  GOTO 2
                ENDIF
                phieau_s2w(i,j)=phieau_s2w(i,j)+phieau_ero_ij
                flx_s2w(:,i,j)=flx_s2w(:,i,j)+flx_s2w_eroij(:)

              ENDIF   ! end test erosi > erodab
            ENDIF    ! end test erosi=0
          ENDIF       ! end test ksmax > ksmi(i,j) 
 
          ! fusion of the two superficial layers, as long as dzs (ksmax) <dzsmin
          !        et si l_consolid=TRUE et couche diminuee par erosion
          !        dans ce cas, il faudrait rajouter un test 
          !        sur la difference entre les concentrations des 2 couches superficielles
          !        pas encore fait
          ! ----------------------------------------------------------
          IF(ero > 0.0_rsh) THEN
            fusion = .true.
            DO WHILE (fusion)
                IF (ksmax.GT.ksmi(i,j)) THEN
                    IF (dzs(ksmax,i,j).LE.dzsmin .AND. (dzs(ksmax-1,i,j)< dzsmax(i,j)  .OR. l_consolid)) THEN
                        ksmaxa=ksmax
                        dzsa=dzs(ksmax,i,j)
                        dzsam1=dzs(ksmax-1,i,j)
                        poroa=poro(ksmax,i,j)
                        poroam1=poro(ksmax-1,i,j)
                        ksmax=ksmax-1
                        dzs(ksmax,i,j)=dzsa+dzsam1
                        dzsisu=1.0_rsh/dzs(ksmax,i,j)
                        c_sedtot(ksmax,i,j)=0.0_rsh
                        cvolp=0.0_rsh
                        DO iv=1,nv_use ! nv_use=nvp (or nvpc if key_Pconstitonly_insed)
                        cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmaxa,i,j)*dzsa             &
                                            +cv_sed(iv,ksmax,i,j)*dzsam1)*dzsisu
                        c_sedtot(ksmax,i,j)=c_sedtot(ksmax,i,j)                      &
                                            +cv_sed(iv,ksmax,i,j)*typart(iv)
                        cvolp=cvolp+cv_sed(iv,ksmax,i,j)*typart(iv)/ros(iv)
                        ENDDO
                        poro(ksmax,i,j)=1.0_rsh-cvolp
                        dzpoi=dzsisu/poro(ksmax,i,j)
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
                        ! dissolved variable in pore waters
                        DO iv=nvp+1,nv_adv
                            ! in .../m3 of interstial waters
                        cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmaxa,i,j)*poroa*dzsa+        &
                                            cv_sed(iv,ksmax,i,j)*poroam1*dzsam1)*dzpoi
                        ENDDO
#endif

                        ! Temperature and salinity
                        DO iv=-1,0
                        cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmaxa,i,j)*poroa*dzsa+        &
                                            cv_sed(iv,ksmax,i,j)*poroam1*dzsam1)*dzpoi
                        ENDDO
#endif
                    ELSE
                        fusion = .false.
                    ENDIF ! dzs(ksmax,i,j).LE.dzsmin .AND. (dzs(ksmax-1,i,j)< dzsmax(i,j)  .OR. l_consolid)
                ELSE
                    fusion = .false.
                ENDIF ! ksmax.GT.ksmi(i,j)
            ENDDO ! end of do while (ksmax.GT.ksmi(i,j))

           IF(ksmax .LT. ksdmax .AND. ksmax > ksmi(i,j)) THEN
              IF(dzs(ksmax,i,j) > dzsmax(i,j) + 5.0_rsh* dzsmin) THEN
                dzs(ksmax+1,i,j)=MIN(dzs(ksmax,i,j)-dzsmax(i,j),dzsmax(i,j))
                dzs(ksmax,i,j)=dzs(ksmax,i,j)-dzs(ksmax+1,i,j)
                poro(ksmax+1,i,j)=poro(ksmax,i,j)
                cv_sed(:,ksmax+1,i,j)=cv_sed(:,ksmax,i,j)
                c_sedtot(ksmax+1,i,j)=c_sedtot(ksmax,i,j)
                ksmax=ksmax+1
               ENDIF
            ENDIF
#ifdef key_MUSTANG_splitlayersurf
           !! Splitting surface layers if too thick
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           isplit=0
           DO k=ksmax,ksmax-nlayer_surf_sed+1,-1
             IF(k > ksmi(i,j)) THEN
               IF(dzs(k,i,j) > dzsmax(i,j) + 5.0_rsh* dzsmin) isplit=1
             ENDIF
           ENDDO
           IF(isplit==1 ) then
             CALL sed_MUSTANG_split_surflayer(i,j,ksmax)
           ENDIF
#endif
          ENDIF  ! ero>0

          ! to find an erosion flux in .../m2/s (for particulates only):
          DO iv=1,nv_use ! nv_use=nvp (or nvpc if key_Pconstitonly_insed)
            flx_s2w(iv,i,j)=flx_s2w(iv,i,j)*dtinv
          ENDDO
          !==============================================
#if ! defined key_noTSdiss_insed
          ! for dissolved variable in pore waters
          ! add water-to-sediment flux explicitly treated with the sediment-to-water flux also explicitly treated
#if ! defined key_Pconstitonly_insed
          !   for dissolved variables, we add the diffusive flux to the interface
          !    and the flux which results from the expulsion of interstitial water by consolidation:
          DO iv=nvp+1,nv_adv
            flx_s2w(iv,i,j)=(flx_s2w(iv,i,j)-flx_w2s(iv,i,j)+fluconsol(iv,i,j)-fludif(iv,i,j))*dtinv
          ENDDO
#endif

          DO iv=-1,0
            flx_s2w(iv,i,j)=(flx_s2w(iv,i,j)-flx_w2s(iv,i,j)+fluconsol(iv,i,j)-fludif(iv,i,j))*dtinv
          ENDDO
#endif
          ! updating ksma (as ksmax can be modified in the routine)
          ksma(i,j)=ksmax
          dzs(ksma(i,j)+1:ksdmax,i,j)=0.0_rsh
          poro(ksma(i,j)+1:ksdmax,i,j)=0.0_rsh

          ! no correction of water column height (could be !)

      END DO
   END DO

  !!  version V1
  END SUBROUTINE sed_MUSTANG_erosion
#endif
   !!==============================================================================

  SUBROUTINE sed_MUSTANG_comp_tocr_mixsed(k, i, j, xeros, excespowr, taucr)

  !&E--------------------------------------------------------------------------
  !&E                 ***  SUBROUTINE sed_MUSTANG_comp_tocr_mixsed  ***
  !&E
  !&E ** Purpose : Compute the erosion parameters in sediment layer k : 
  !&E   - xeros : erosion constant (kg.m-2.s-1)
  !&E   - excespowr : power of the excess shear stress
  !&E   - taucr : critical shear stress (N.m-2)
  !&E
  !&E ** Description : 0D in a cell, variables OUT : xeros,excespowr,taucr
  !&E          Variables used from module
  !&E                         cv_sed, ros_sand_homogen, ros
  !&E                         diam_sed, stresscri0
  !&E                         E0_sand, n_eros_sand, coef_frmudcr1, frmudcr2
  !&E                         E0_sand_option, E0_sand_para, E0_sand_Cst
  !&E                         E0_mud, x1toce_mud, x2toce_mud, n_eros_mud
  !&E                         l_xexp_ero_cst
  !&E                         in V2 : crel_mud
  !&E  
  !&E  !**************************************************************!
  !&E  ! WARNING : the proposed formulas for erosion of               !
  !&E  ! sandy-muddy mixtures are questionable                        !
  !&E  ! and must be tested, improved, changed ...                    !
  !&E  ! Each user can write his law                                  !  
  !&E  !**************************************************************!
  !&E
  !&E ** Called by :  sed_MUSTANG_erosion
  !&E 
  !&E ** External calls : none
  !&E
  !&E--------------------------------------------------------------------------

  !! * Arguments
  INTEGER,INTENT(in) :: k, i, j           ! cell location (i,j) and layer location (k)
  REAL(KIND=rsh),INTENT(out) :: xeros     ! erosion constant (kg.m-2.s-1)
  REAL(KIND=rsh),INTENT(out) :: excespowr ! power of the excess shear stress
  REAL(KIND=rsh),INTENT(out) :: taucr      ! ! critical shear stress (N.m-2)

  !! * Local declarations
  INTEGER            :: iv  ! sediment class
  REAL(KIND=rsh)     :: diamsan ! diameter characteristic of gravel and sand present in the layer (balanced sum of diam of gravels and sands) in m
  REAL(KIND=rsh)     :: taucr_sand ! critical shear stress characteristic of gravel and sand present in the layer(balanced sum of tauc_susp of gravels and sands) in N.m-2 
  REAL(KIND=rsh)     :: cmudr ! relative mud concentration in kg.m-3
  REAL(KIND=rsh)     :: frmudsup ! mud fraction 
  REAL(KIND=rsh)     :: frmudcr1 ! critical mud fraction under which the behaviour is purely sandy
  REAL(KIND=rsh)     :: sommud, somsan, somgrav
  REAL(KIND=rsh)     :: frvolsan, frvolgrv,  frvolsangrv
  REAL(KIND=rsh)     :: pinterp, taucr_mud
  REAL(KIND=rsh)     :: diamsanstar, wssand
  ! use of the van Rijn relation for erosion flux (only valid if sand only)
  REAL(KIND=rsh)     :: rossan, xeromud, coef_tmp, rapexpcoef
  REAL(KIND=rsh)     :: E0_sand_loc

  !!---------------------------------------------------------------------------
  !! * Executable part

  ! compute each type (mud/sand/gravel) quantities and fraction
  sommud = 0.0_rsh     
  somsan = 0.0_rsh
  somgrav = 0.0_rsh
  frvolsan = 0.0_rsh   
  frvolgrv = 0.0_rsh
  ! calculation of diamsan and taucr_sand from gravels and sands properties
  diamsan = 0.0_rsh
  taucr_sand = 0.0_rsh
  E0_sand_loc = 0.0_rsh
  rossan = ros_sand_homogen !!! WARNING: EVEN IF SEVERAL SANDS, WE ASSUME THAT THEY HAVE THE SAME DENSITY

#ifdef key_MUSTANG_V2 /* gravel & sands */
  DO iv = imud1, imud2
    sommud = sommud + cv_sed(iv,k,i,j)
  ENDDO
  DO iv = isand1, isand2
    somsan = somsan + cv_sed(iv,k,i,j)
    frvolsan = frvolsan + cv_sed(iv,k,i,j) / ros(iv)
  ENDDO
  DO iv = igrav1, igrav2
    somgrav = somgrav + cv_sed(iv,k,i,j)
    frvolgrv = frvolgrv + cv_sed(iv,k,i,j) / ros(iv)
  ENDDO
  frvolsangrv = frvolgrv + frvolsan
  frmudsup = sommud / (somgrav + sommud + somsan)

  DO iv = igrav1, isand2   
    diamsan = diamsan + diam_sed(iv) * cv_sed(iv,k,i,j)
    taucr_sand = taucr_sand + stresscri0(iv) * cv_sed(iv,k,i,j)
    E0_sand_loc = E0_sand_loc + E0_sand(iv) * cv_sed(iv,k,i,j)
  ENDDO
  IF (isand2 > 0 .AND. (somsan + somgrav) > 0.0_rsh) THEN ! if there is at least one sand or gravel defined and in the layer k 
    diamsan = MAX(diamsan / (somsan + somgrav + epsilon_MUSTANG), diam_sed(isand2))
    taucr_sand = MAX(taucr_sand / (somsan + somgrav + epsilon_MUSTANG), stresscri0(isand2))
    E0_sand_loc = MAX(E0_sand_loc / (somsan + somgrav + epsilon_MUSTANG), E0_sand(isand2))
  ENDIF
  cmudr = crel_mud(k,i,j)

#else  /*version V1 :  gravel are not taking into account in V1*/     
  DO iv = imud1, imud2
    sommud = sommud + cv_sed(iv,k,i,j)
  ENDDO
  DO iv = isand1, isand2
    IF (diam_sed(iv) .LT. 0.002) THEN   ! to remove gravels that are declared as sand 
      ! in our configuration before computing mean parameters 
      !(i.e. gravels are not working and are declared as sands in Mustang V1)
        somsan = somsan + cv_sed(iv,k,i,j)
        frvolsan = frvolsan + cv_sed(iv,k,i,j) / ros(iv)
        diamsan = diamsan + diam_sed(iv) * cv_sed(iv,k,i,j)
        taucr_sand = taucr_sand + stresscri0(iv) * cv_sed(iv,k,i,j)
    ENDIF
  ENDDO
  DO iv = igrav1, igrav2
    somgrav = somgrav + cv_sed(iv,k,i,j)
    frvolgrv = frvolgrv + cv_sed(iv,k,i,j) / ros(iv)
  ENDDO
  frvolsangrv = frvolgrv + frvolsan
  frmudsup = sommud / (somgrav + sommud + somsan)
  
   IF (isand2 > 0 .AND. somsan > 0.0_rsh) THEN
     diamsan = MAX(diamsan / (somsan + epsilon_MUSTANG), diam_sed(isand2))
     taucr_sand = MAX(taucr_sand / (somsan + epsilon_MUSTANG), stresscri0(isand2))
     diamsanstar = diamsan * 10000.0_rsh * (GRAVITY * (rossan / RHOREF - 1.0_rsh))**0.33_rsh
     ! according to Soulsby, 1997, and if viscosity = 10-6 m/s :
     wssand = .000001_rsh * ((107.33_rsh + 1.049_rsh * diamsanstar**3)**0.5_rsh - 10.36_rsh) / diamsan
     E0_sand_loc = MUSTANG_E0sand(diamsan, taucr_sand, rossan, wssand) 
   ENDIF
   cmudr = sommud / (1.0_rsh - frvolsangrv)
#endif /*version V1/V2 */


  IF (sommud + somsan .LE. 0.0_rsh) THEN ! sediment could not pass in suspension : only diam > 0.002, gravels
#ifdef key_MUSTANG_V2 /*taucr needed for activelayer calculation*/
    IF (ero_option == 0) THEN ! keep mud constant values
      xeros = E0_mud
      taucr = taucr_mud
      excespowr = n_eros_mud
    ELSE
      xeros = E0_sand_loc
      taucr = taucr_sand
      excespowr = n_eros_sand
    ENDIF
#else /*version V1, gravel not in suspension*/
    xeros = 0.0_rsh
    taucr = 1000.0_rsh
    excespowr = 0.0_rsh
#endif
  ELSE ! there is mud and/or sand

   frmudcr1 = MIN(coef_frmudcr1 * diamsan, frmudcr2)
   taucr_mud = x1toce_mud * cmudr**x2toce_mud
   
    IF (ero_option == 0) THEN ! keep mud constant values

      xeros = E0_mud
      taucr = taucr_mud
      excespowr = n_eros_mud

    ELSE IF (ero_option == 1) THEN

      IF(frmudsup  .LE. frmudcr1) THEN ! I) Sandy behavior of the mixture
        xeros = E0_sand_loc
        taucr = taucr_sand
        excespowr = n_eros_sand
      ELSE IF(frmudsup .LE. frmudcr2) THEN ! II) Intermediate sand / mud 
        pinterp = (frmudcr2 - frmudsup) / (frmudcr2 - frmudcr1)
        !!! Initial formulation with linear transition in the case of an intermediate mixture
        !!! Le Hir et al (2011) 
        xeros = pinterp * E0_sand_loc + (1.0_rsh - pinterp) * E0_mud
        excespowr = pinterp * n_eros_sand + (1.0_rsh - pinterp) * n_eros_mud
        taucr = taucr_sand * pinterp + (1.0_rsh - pinterp) * taucr_mud
      ELSE ! III) Mud Behavior 
        xeros = E0_mud
        taucr = taucr_mud
        excespowr = n_eros_mud
      ENDIF ! frmudsusp

    ELSE IF (ero_option .EQ. 2) THEN

      IF(frmudsup  .LE. frmudcr1) THEN ! I) Sandy behavior of the mixture
        xeros = E0_sand_loc
        taucr = taucr_sand
        excespowr = n_eros_sand
      ELSE IF(frmudsup .LE. frmudcr2) THEN ! II) Intermediate sand / mud 
        pinterp = (frmudcr2 - frmudsup) / (frmudcr2 - frmudcr1)
        !!! Formulation of Julie Vareilles (2013) 
        !!! according to Van Ledden (2001) and Carniello et al (2012)
        xeros = ( (E0_sand_loc / E0_mud)**pinterp ) * E0_mud
        taucr = ( (taucr_sand / taucr_mud)**pinterp ) * taucr_mud
        excespowr = ( (n_eros_sand / n_eros_mud)**pinterp ) * n_eros_mud
      ELSE ! III) Mud Behavior 
        xeros = E0_mud
        taucr = taucr_mud
        excespowr = n_eros_mud
      ENDIF ! frmudsusp

    ELSE IF (ero_option .EQ. 3) THEN  

        !!! Formulation of Baptiste Mengual (2017)
        !!! Mengual, B.; Le Hir, P.; Cayocca, F.; Garlan, T. 
        !!! Modelling Fine Sediment Dynamics: Towards a Common Erosion Law for Fine Sand, Mud and Mixtures. Water 2017, 9, 564.
        !!! Exponential Formulations of Parameters for a Sand/mud Mixture            
#ifdef key_MUSTANG_V2
      IF(frmudsup  .LE. frmudcr1) THEN ! I) Sandy behavior of the mixture
        xeros = E0_sand_loc
        taucr = taucr_sand
        excespowr = n_eros_sand
      ELSE ! II) Intermediate sand / mud, no frmudcr2 in V2
        IF (.NOT. l_xexp_ero_cst) xexp_ero = 78.48_rsh * frmudcr1 - 2.037_rsh !xexp_ero=f(diamsan) with min/max values based on frmudcr1
                                                                  !Mainly based on Fig 2 in Le Hir et al. CSR (2011) but
                                                                  !consitent with other experiments from the literature
        coef_tmp = xexp_ero * (frmudcr1 - frmudsup) ! Here, frmudcr2 is just a max prescribed for frmudcr1 
                                                ! but does not constitute a critical mud fraction 
        ! the correction of F.Ganthy is not done here for version V2. 
        !  Should we add it to do as in version V1 ?
        rapexpcoef = EXP(coef_tmp)
        taucr = (taucr_sand - taucr_mud) * rapexpcoef + taucr_mud
        excespowr = (n_eros_sand - n_eros_mud) * rapexpcoef + n_eros_mud
        xeros = (E0_sand_loc - E0_mud) * rapexpcoef + E0_mud 
      END IF !frmudsusp
#else /* V1 */
      IF(frmudsup  .LE. frmudcr1) THEN ! I) Sandy behavior of the mixture
        xeros = E0_sand_loc
        taucr = taucr_sand
        excespowr = n_eros_sand
      ELSE IF(frmudsup .LE. frmudcr2) THEN ! II) Intermediate sand / mud 
        coef_tmp = xexp_ero * (frmudcr1 - frmudsup) / (frmudcr2 - frmudcr1) 
        ! correction of  F.Ganthy which allows to avoid the shift when one approaches frmudcr2, 
        !   and allows to have a linear relation when xexp_ero tends towards 0 (but must remain different from 0)
        rapexpcoef = (EXP(coef_tmp) - 1.0_rsh) / (EXP(xexp_ero) - 1.0_rsh)
        taucr = (taucr_sand - taucr_mud) * rapexpcoef + taucr_mud
        excespowr = (n_eros_sand - n_eros_mud) * rapexpcoef + n_eros_mud
        xeros = (E0_sand_loc - E0_mud) * rapexpcoef + E0_mud 
      ELSE ! III) Mud Behavior 
        xeros = E0_mud
        taucr = taucr_mud
        excespowr = n_eros_mud
      ENDIF ! frmudsusp
#endif    

    ENDIF ! ero_option = 0/1/2/3 
  ENDIF ! sommud + somsand 

  taucr = MAX(0.00005_rsh, taucr)

#ifdef key_MUSTANG_specif_outputs
  varspecif2D_save(1,i,j)=frmudsup           
#endif
  END SUBROUTINE sed_MUSTANG_comp_tocr_mixsed
   !!==============================================================================

  SUBROUTINE sed_MUSTANG_comp_eros_flx(tenfo, toce, excespowr, xeros, sed_eros_flx)
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  SUBROUTINE sed_erosion_flux  ***
   !&E
   !&E ** Purpose : computes erosion flux for given stress, critical stress and
   !&E              surficial sediment concentration
   !&E              uses toce, excespowr and xeros computed by sed_tocr
   !&E
   !&E ** Description : Partheniades or adapted Partheniades
   !&E
   !&E ** Called by :  sed_erosion
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
   !! * Arguments
   REAL(KIND=rsh),INTENT(in) :: tenfo, toce
   REAL(KIND=rsh),INTENT(in) :: excespowr, xeros
   REAL(KIND=rsh),INTENT(out) :: sed_eros_flx

   !!---------------------------------------------------------------------------
   !! * Executable part
   
   ! Partheniades :
   sed_eros_flx = xeros * (tenfo/toce - 1.0_rsh)**excespowr
!   or choice one other formula 
!  sed_eros_flx=xeros*corfluero*(tenfo-toce)**excespowr
    
  END SUBROUTINE sed_MUSTANG_comp_eros_flx

   !!==============================================================================
#ifdef key_MUSTANG_V2
 SUBROUTINE sed_MUSTANG_effdep(ifirst, ilast, jfirst, jlast, iexchge_MPI_cvwat)
    
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_effdep version V2 ***
   !&E
   !&E ** Purpose : computes effective deposition
   !&E
   !&E ** Description : 
   !&E
   !&E       arguments IN :
   !&E          loops  :ifirst,ilast,jfirst,jlast
   !&E          CELL_SURF
   !&E          
   !&E       variables OUT :
   !&E           phieau_s2w
   !&E           flx_w2s
   !&E
   !&E ** Called by :  MUSTANG_deposition
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used

   !! * Arguments
   INTEGER, INTENT(IN)  :: ifirst,ilast,jfirst,jlast
   INTEGER, INTENT(INOUT)  :: iexchge_MPI_cvwat

   !! * Local declarations
   REAL(KIND=rsh),DIMENSION(nvp)  :: flx_w2s_loc
   REAL(KIND=rsh),DIMENSION(nvpc) :: frdep,flx_w2s_loca          
   INTEGER                        :: i,j,k,iv,ksmin,ksmax,kl,ivp_assoc,isplit
   REAL(KIND=rsh) :: fludep,ddzs,ddzs1,ddzs2,dzsa,dzsisu,dzsi,dzsaici,      &
                     dzsgrv,dzssan,dzsmud,dzsnew,dzsa2,dzsa3,ddzs3,ddzs4,ddzs5,ddzs6,    &
                     poroa,somalp,somala,cdepiv,porosi,porodep,fluxdissous,              &
                     ddzsici,voldepgrv,voldepsan,masdepmud,dzskmanew,porewatera,         &
                     sommud,cvolinigrv,cvolinisan,cvolp,cmudr,porewater,porewaterdep,    &
                     dvolgrv,dvolsan,dmasmud,dmasmudgrav,dmasmudsand,cordepflu,          &
                     cordepfluw,cordepflue,cordepflus,cordepflun,qdep,vitdepo
   LOGICAL        :: l_createnewlayer,l_increase_dep,l_isitcohesive_dep
   REAL(KIND=rsh) :: mass_tot_dep, poro_mud_dep, poro_dep, dzs_dep, &
                     mass_tot, poro_muda, poro_mud_new,             &
                     dzsmin_dep,dzsmina,crel_mud_new,crel_mud_dep,  &
                     diamgravsan_dep,frmud_dep,frmudcr_dep,mass_tot_depa
                 
   REAL(KIND=rsh),DIMENSION(nvp) :: mass_sed
   REAL(KIND=rsh),DIMENSION(nvpc):: frac_sed_dep, cv_sed_dep, frac_sed, frac_seda,frac_sed_depa
#ifdef key_MUSTANG_bedload
   REAL(KIND=rsh),DIMENSION(nvp):: flx_bedload_in
#endif
   !!----------------------------------------------------------------------
   !! * Executable part
#if defined key_MUSTANG_debug
   IF (l_debug_effdep .AND. CURRENT_TIME> t_start_debug .AND. htot(i_MUSTANG_debug,j_MUSTANG_debug) > h0fond ) THEN
     print *,''
     print *,' ************************'
     print *,' ENTER sed_effdep_mixsed'
     print *,' ************************'
   ENDIF        
#endif

   ddzsici=-1000.0_rsh
   iexchge_MPI_cvwat=0       

      DO j=jfirst,jlast
        DO i=ifirst,ilast
          ! test 0 (water or not)
          IF(htot(i,j) > h0fond) THEN

            ksmax=ksma(i,j)
            ksmin=ksmi(i,j)
            flx_w2s_loc(:)=0.0_rsh
            flx_w2s_loca(:)=0.0_rsh
            frdep(:)=0.0_rsh
            frac_sed_depa(:)=0.0_rsh

#ifdef key_MUSTANG_debug
               IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
                 print *,'  > deb effdep',i,j
                 print *,'  t=',CURRENT_TIME, 'ksmax=',ksmax
                 print *,'  dzs(ksmax-3:ksmax,i,j)=',dzs(ksmax-3:ksmax,i,j)
                 print *,'  cv_sed(:,ksmax-3,i,j)=',cv_sed(:,ksmax-3,i,j)
                 print *,'  cv_sed(:,ksmax-2,i,j)=',cv_sed(:,ksmax-2,i,j)
                 print *,'  cv_sed(:,ksmax-1,i,j)=',cv_sed(:,ksmax-1,i,j)
                 print *,'  cv_sed(:,ksmax,i,j)=',cv_sed(:,ksmax,i,j)
               END IF
#endif 

#ifdef key_MUSTANG_bedload
            flx_bedload_in(:)=0.0_rsh
            ! bedload fluxes
            !   ATTENTION : need to know fls_bx in i+1,J+1,i-1,j-1
            !               have been exchanged with the neighboring processors  in sedim_MUSTANG_update after erosion
            DO iv=ibedload1,ibedload2

                ! en kg
                flx_bedload_in(iv)=-MIN(flx_bx(iv,i+1,j),0.0_rsh)+MAX(flx_bx(iv,i-1,j),0.0_rsh) &
                                   -MIN(flx_by(iv,i,j+1),0.0_rsh)+MAX(flx_by(iv,i,j-1),0.0_rsh)

#ifdef key_MUSTANG_specif_outputs
                ! bil_bedload(iv,i,j) a mettre a 0 en debut de run --> cumule tout au long de la simu
                varspecif3Dnv_save(7,iv,i,j) = varspecif3Dnv_save(7,iv,i,j) + ( (flx_bedload_in(iv)  &
                                - ABS(flx_bx(iv,i,j)) - ABS(flx_by(iv,i,j)))/CELL_SURF(i,j) ) ! cumul des bilans en kg/m2 
#endif

                ! in kg/m2
                flx_bedload_in(iv)=flx_bedload_in(iv)/CELL_SURF(i,j)

                flx_w2s_sum(iv,i,j)=flx_w2s_sum(iv,i,j)+flx_bedload_in(iv) ! en kg/m2

            END DO
            DO iv=imud2+1,nvpc

                flx_bedload_in(iv)=-MIN(flx_bx(iv,i+1,j),0.0_rsh)+MAX(flx_bx(iv,i-1,j),0.0_rsh) &
                                  -MIN(flx_by(iv,i,j+1),0.0_rsh)+MAX(flx_by(iv,i,j-1),0.0_rsh)

                ! in kg/m2
                flx_bedload_in(iv)=flx_bedload_in(iv)/CELL_SURF(i,j)

                flx_w2s_sum(iv,i,j)=flx_w2s_sum(iv,i,j)+flx_bedload_in(iv) 

            END DO
            
#ifdef key_MUSTANG_specif_outputs
            !bil_bedload_int
            varspecif2D_save(18,i,j)=SUM(varspecif3Dnv_save(7,ibedload1:ibedload2,i,j))
#endif
#endif
           
            ! updating effective deposition :  (flx_w2s) is implicit in the vertical advection scheme
            !      flux exprime en Masse/m2 integre sur le  pas de temps vrai (demi pas de temps dans MARS)
            ! loop including gravels as there may now be gravel coming from adjacent meshes
            ! In addition the algorithm was written to take into account a deposit of gravel

            DO iv=1,nv_use   ! nv_use=nvp (or nvpc if key_Pconstitonly_insed)
              flx_w2s_loc(iv)=MAX(0.0_rsh,flx_w2s_sum(iv,i,j)+flx_w2s_corin(iv,i,j)      &
                           +flx_w2s_corip1(iv,i-1,j)+flx_w2s_corim1(iv,i+1,j)         &
                           +flx_w2s_corjp1(iv,i,j-1)+flx_w2s_corjm1(iv,i,j+1))

              flx_w2s_sum(iv,i,j)=0.0_rsh

              ! MF /= 1 only if l_morphocoupl
              flx_w2s_loc(iv) = MF * flx_w2s_loc(iv)
              
#if defined key_MUSTANG_debug
               IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. &
                       ( CURRENT_TIME> t_start_debug)) THEN
                 print *,'flx_w2s_loc(',iv,')=',flx_w2s_loc(iv)
              END IF             
#endif
            ENDDO

            fludep=0.0_rsh
            DO iv=1,nvpc
             fludep=fludep+flx_w2s_loc(iv)
#ifdef key_MUSTANG_specif_outputs
              !flx_w2s_save
#ifdef key_MUSTANG_bedload
             varspecif3Dnv_save(3,iv,i,j)=varspecif3Dnv_save(3,iv,i,j)  &
                                       +flx_w2s_loc(iv) -flx_bedload_in(iv)
#else
             varspecif3Dnv_save(3,iv,i,j)=varspecif3Dnv_save(3,iv,i,j)+flx_w2s_loc(iv)
#endif
#endif
            ENDDO

#ifdef key_MUSTANG_specif_outputs
            DO iv=isand1,isand2
#ifdef key_MUSTANG_bedload
              varspecif2D_save(12,i,j)=varspecif2D_save(12,i,j)  &
                                       +flx_w2s_loc(iv)-flx_bedload_in(iv)  ! flx_w2s_noncoh
#else
              varspecif2D_save(12,i,j)=varspecif2D_save(12,i,j)+flx_w2s_loc(iv)
#endif
            END DO
            DO iv=imud1,imud2
              varspecif2D_save(14,i,j)=varspecif2D_save(14,i,j)+flx_w2s_loc(iv) ! flx_w2s_coh
            END DO
#endif
            voldepgrv=0.0_rsh
            voldepsan=0.0_rsh
            masdepmud=0.0_rsh
              
            !test 1 (deposit)
            IF(fludep.GT.0.0_rsh)THEN
               l_increase_dep=.FALSE.

#ifdef key_MUSTANG_debug
               IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND.   &
                       (CURRENT_TIME> t_start_debug)) THEN
                 print *,'fludep=',fludep,' > 0 --> there is deposition'
                 print *,''
                 print *,'  > Charac depsition'
              END IF
#endif
              ! case 1:  there is deposition
              ! ****************************

              ! characterising deposition:
              DO iv=igrav1,igrav2
                voldepgrv=voldepgrv+flx_w2s_loc(iv)/ros(iv) ! in m
              ENDDO
              DO iv=igrav1,igrav2
                frdep(iv)=flx_w2s_loc(iv)/(voldepgrv+epsi30_MUSTANG) ! in kg/m3
              ENDDO
              DO iv=isand1,isand2
                voldepsan=voldepsan+flx_w2s_loc(iv)/ros(iv) ! in m
              ENDDO
              DO iv=isand1,isand2
                frdep(iv)=flx_w2s_loc(iv)/(voldepsan+epsi30_MUSTANG) ! in kg/m3
              ENDDO
              DO iv=imud1,imud2
                masdepmud=masdepmud+flx_w2s_loc(iv) ! in kg/m2
              ENDDO
              DO iv=imud1,imud2
                frdep(iv)=flx_w2s_loc(iv)/(masdepmud+epsi30_MUSTANG) ! SU
              ENDDO 

              mass_tot_dep=0.0_rsh
              DO iv=igrav1,imud2
                mass_tot_dep=mass_tot_dep+flx_w2s_loc(iv)
              END DO 

              DO iv=igrav1,imud2
                frac_sed_dep(iv)=flx_w2s_loc(iv)/mass_tot_dep
              END DO
   
              poro_mud_dep=1.0_rsh-(cfreshmud/ros(1))

              crel_mud_dep=cfreshmud

              CALL MUSTANGV2_comp_poro_mixsed(frac_sed_dep, poro_mud_dep, &
                                 crel_mud_dep, poro_dep)

              dzs_dep=mass_tot_dep/((1.0_rsh-poro_dep)*ros(1))
             
              !dzsmin_dep=(1.0_rsh-coeff_dzsmin)*dzsminuni +   &
              !            coeff_dzsmin*SUM( frac_sed_dep(1:nvpc)*diam_sed(1:nvpc) )
              dzsmin_dep=dzsminvar(frac_sed_dep)

              l_createnewlayer=.TRUE. ! According to different tests, the booleen will be re-refreshed

#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
                ! to keep in mind the amount of water coming from the water column by the deposit (used in case of l_increase_dep=T)
                 porewaterdep=poro_dep*dzs_dep
#endif


#ifdef key_MUSTANG_debug
               IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND.   &
                       (CURRENT_TIME> t_start_debug)) THEN
                print *,'  poro_mud_dep=',poro_mud_dep
                print *,'  poro_dep=',poro_dep
                print *,'  frac_sed_dep=',frac_sed_dep
                print *,'  mass_tot_dep=',mass_tot_dep
                print *,'  dzs_dep=',dzs_dep
                print *,'  dzsmin_dep=',dzsmin_dep
              END IF
#endif

              !test  (sediment exist)
              IF(ksmax.GE.ksmi(i,j))THEN
          
               ! case 1.1: sediment already exists
               ! ---------------------------------
               ! characterising previous surficial sediment:

               sommud=0.0_rsh      
               DO iv=imud1,imud2
                 sommud=sommud+cv_sed(iv,ksmax,i,j)
               ENDDO
               cvolinigrv=0.0_rsh      
               DO iv=igrav1,igrav2
                 cvolinigrv=cvolinigrv+cv_sed(iv,ksmax,i,j)/ros(iv)
               ENDDO
               cvolinisan=0.0_rsh
               DO iv=isand1,isand2
                 cvolinisan=cvolinisan+cv_sed(iv,ksmax,i,j)/ros(iv)
               ENDDO
               cmudr=crel_mud(ksmax,i,j)

               poro_muda=poro_mud(ksmax,i,j)
               poroa=poro(ksmax,i,j)
               dzsa=dzs(ksmax,i,j)

#ifdef key_MUSTANG_debug
               IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND.   &
                       (CURRENT_TIME> t_start_debug)) THEN
                 print *,''
                 print *,'ksmax=',ksmax,' >= ',ksmi(i,j)
                 print *,''
                 print *,'  > Charac existing sed in ksmax'
               END IF
#endif 
               frac_seda(:)=0.0_rsh
               DO iv=1,nvpc
                 frac_seda(iv)=cv_sed(iv,ksmax,i,j)/c_sedtot(ksmax,i,j)
               END DO
               !dzsmina=(1.0_rsh-coeff_dzsmin)*dzsminuni +  &
                !                coeff_dzsmin*SUM( frac_seda(1:nvpc)*diam_sed(1:nvpc) )
               dzsmina=dzsminvar(frac_seda)

#ifdef key_MUSTANG_debug
               IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND.   &
                       (CURRENT_TIME> t_start_debug)) THEN
                 print *,'   sommud=',sommud
                 print *,'   cvolinigrv=',cvolinigrv
                 print *,'   cvolinisan=',cvolinisan
                 print *,'   cmudr (crel_mud)=',cmudr
                 print *,'   poro_muda=',poro_muda
                 print *,'   poroa=',poroa
                 print *,'   dzsa=',dzsa
                   print *,'   frac_seda=',frac_seda
                 print *,'   dzsmina=',dzsmina
                 print *,''
                 print *,'MIXING OR NOT ?'
                 print *,''
               END IF
#endif                          
               !test  (mixing layer)
               IF (cmudr .LE. cmudcr) THEN
               ! The surface sediment in ksmax is not consolidated, 
               ! the mixture is therefore the solution envisaged, while taking care to respect the dzsmax criterion.
               ! Note: a layer thickness may exceed dzsmax if dzs (ksmax, i, j) <dzsmax (i, j) at the beginning of 
               !       the time step and it is decided to mix with a significant deposit.
               !       It is only at the next time step that we take into account that we have exceeded dzsmax.
                 IF( (dzsa .LT. MAX(dzsmax(i,j), 2.0_rsh*dzsmina))) THEN
                   ! If we want to mix and we can do it while respecting the dzsmax criterion
                   ! Important note: we take the MAX (dzsmax (i, j), 2.0_rsh * dzsmina) for cases where we would prescribe 
                   ! a dzsmax smaller than the diameter of the coarsest class ... 
                   ! (typically 5 mm in the Seine configuration then we have a gravel of 10 mm)
                   !  In this case -> it will be a problem if we increase a deposition with dzsmina (l_increase_dep = True) 
                   !   ... because there will be a risk of negative thickness of the existing layer
                   l_createnewlayer=.FALSE.
                   l_increase_dep=.FALSE.   
                 ELSE
                   ! If we want to mix but we do not because the existing layer has already reached the max 
                   ! between dzsmax and 2 * dzsmina
                   ! So we create a layer by asking the question of the physical meaning of the quantities 
                   ! that we are preparing to deposit
                   l_createnewlayer=.TRUE.
                   IF (dzs_dep .LT. dzsmin_dep) THEN
                     ! In this case the quantities to be deposited are very small and considered not sufficient to justify 
                     ! the creation of a new layer.
                     ! So we decide to take some of the existing ksmax layer to increase this deposit
                     l_increase_dep=.TRUE.
                   ELSE
                     ! In this case, the quantities to be deposited are sufficient,
                     ! a new layer is created with the amounts to be deposited.
                     l_increase_dep=.FALSE.
                   END IF
                 END IF
               ELSE
                 ! The sediment in ksmax is consolidated and we want to preserve the layer
                 ! We decide to create a new layer by asking the question
                 ! 1- the physical meaning of the quantities that we have to deposit
                 ! 2- the nature of the deposit, cohesive or not
                 ! 3- the thickness of the existing ksmax layer
                 IF (dzs_dep .LT. dzsmin_dep) THEN
                   ! In this case, the quantities to be deposited are very small and are not
                   ! sufficient to justify the creation of a new layer

                   ! Test whether the sediment to be deposited is of the non-cohesive or cohesive type
                   diamgravsan_dep=0.0_rsh
                   DO iv=igrav1,isand2
                     diamgravsan_dep=diamgravsan_dep+diam_sed(iv)*(frac_sed_dep(iv)  &
                                       /(SUM(frac_sed_dep(igrav1:isand2))+epsilon_MUSTANG))
                   ENDDO
                   IF (isand2>0) THEN
                    diamgravsan_dep=MAX(diamgravsan_dep,diam_sed(isand2))
                   END IF

                   frmud_dep=SUM(frac_sed_dep(imud1:imud2))
                   frmudcr_dep=MIN(coef_frmudcr1*diamgravsan_dep,frmudcr2)

                   IF (frmud_dep .GE. frmudcr_dep) THEN
                     l_isitcohesive_dep=.TRUE.
                   ELSE
                     l_isitcohesive_dep=.FALSE.
                   END IF

                   IF ( (.NOT. l_isitcohesive_dep) .AND. (dzsa .LT. dzsmax(i,j)) ) THEN
                      ! The sediments to be deposited are essentially non-cohesive
                      ! In this case, in order to avoid overly artificially increasing our
                      ! erodable vase stock in ksmax, we decided to mix these non-cohesive deposits
                      ! with the existing layer until reaching dzsmax                 
                       l_createnewlayer=.FALSE.
                       l_increase_dep=.FALSE.
                   ELSE
                     IF (dzsa .GT. 2.0_rsh*dzsmina) THEN
                       ! So we decide to take some of the existing ksmax layer to increase
                       ! this deposit because it is thick enough to be reduced by half                    
                       l_createnewlayer=.TRUE.
                       l_increase_dep=.TRUE.
                     ELSE
                       ! The ksmax layer is too thin to increase the deposit
                       ! -> we mix                    
                       l_createnewlayer=.FALSE.
                       l_increase_dep=.FALSE.
                     END IF
                  END IF
                 ELSE
                  ! In this case, the quantities to be deposited are sufficient and a new layer is created.
                  l_createnewlayer=.TRUE.
                  l_increase_dep=.FALSE.
                 END IF
               END IF
               !!! 
               !IF (l_createnewlayer==.FALSE. .AND. l_increase_dep==.TRUE.) &
               !   print *,' !!! l_createnewlayer==.FALSE. .AND. l_increase_dep==.TRUE. !!! --> contradictory'
               
            
               IF (.NOT. l_createnewlayer) THEN

                 !!!!!!!!!!!!!!!!!!!!!!!
                 !!!!!! MIXING !!!!!!!!!
                 !!!!!!!!!!!!!!!!!!!!!!!
           
           
                 !print *,' ==> l_createnewlayer=',l_createnewlayer

#ifdef key_MUSTANG_debug
                 IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND.   &
                       (CURRENT_TIME> t_start_debug)) THEN
                    print *,'  > MIXING of deposits with sed in ksmax, because :'
                    IF ((cmudr .LE. cmudcr) .AND. (dzsa .LT. dzsmax(i,j))) &
                      print *,cmudr,' (cmudr) <= ',cmudcr,' (cmudcr) AND ',dzsa,' (dzsa) < ',dzsmax(i,j),' (dzsmax(i,j))'
                    IF ((dzs_dep .LT. dzsmin_dep) .AND. (dzsa .LT. 2.0_rsh*dzsmina)) &
                       print *,dzs_dep,' (dzs_dep) < ',dzsmin_dep,' (dzsmin_dep) AND ',dzsa, &
                                 ' (dzsa) < ',2.0_rsh*dzsmina,' (2*dzsmina)'
                    print *,'  l_createnewlayer=',l_createnewlayer, '*******MIXING********'
                  END IF
#endif            

                  ! mixing deposits with upper layer with no restriction
                  ! if surficial sediment is not consolidated and if the upper
                  ! layer thickness at the beginning of the true dt (1/2 dt in MARS) does not exceed dzsmax.
                  ! Simplification: Henceforth, dzsmax /dzsmin criteria do not intervene within a time step 

                  IF (masdepmud+sommud*dzsa .GT. 0.0_rsh) THEN
                    poro_mud_new=(masdepmud/(masdepmud+sommud*dzsa))*poro_mud_dep + &
                                 ((sommud*dzsa)/(masdepmud+sommud*dzsa))*poro_muda
                    crel_mud_new=(masdepmud/(masdepmud+sommud*dzsa))*crel_mud_dep + &
                                 ((sommud*dzsa)/(masdepmud+sommud*dzsa))*crel_mud(ksmax,i,j)
                  ELSE
                    poro_mud_new=0.0_rsh
                    crel_mud_new=0.0_rsh
                  END IF

                  mass_tot = 0.0_rsh
                  DO iv=igrav1,imud2
                     mass_sed(iv)=cv_sed(iv,ksmax,i,j)*dzsa + flx_w2s_loc(iv)
                     mass_tot=mass_tot+mass_sed(iv)
                  END DO

                  DO iv=igrav1,imud2
                    frac_sed(iv)=mass_sed(iv)/mass_tot
                  END DO

#ifdef key_MUSTANG_debug
                 IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND.   &
                       (CURRENT_TIME> t_start_debug)) THEN
                    print *,'  > MIXING of deposits with sed in ksmax, : avant comp_poro_mixsed'
                    print *,'frac_sed = ',frac_sed
                    print *,'poro_mud_new = ',poro_mud_new
                    print *,'crel_mud_new = ',crel_mud_new
                    print *,'poro = ',poro(ksmax,i,j)
                  END IF
#endif            

                  CALL MUSTANGV2_comp_poro_mixsed(frac_sed, poro_mud_new,   &
                                crel_mud_new, poro(ksmax,i,j))

#ifdef key_MUSTANG_debug
                 IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND.   &
                       (CURRENT_TIME> t_start_debug)) THEN
                    print *,'  > MIXING of deposits with sed in ksmax, : apres comp_poro_mixsed'
                    print *,'frac_sed = ',frac_sed
                    print *,'poro_mud_new = ',poro_mud_new
                    print *,'crel_mud_new = ',crel_mud_new
                    print *,'poro = ',poro(ksmax,i,j)
                  END IF
#endif            
                  dzs(ksmax,i,j)=mass_tot/((1.0_rsh-poro(ksmax,i,j))*ros(1))
                  dzsi=1.0_rsh/dzs(ksmax,i,j)

                  DO iv=igrav1,imud2
                    cv_sed(iv,ksmax,i,j)=mass_sed(iv)*dzsi
                  END DO
  
                  IF (dzs(ksmax,i,j) .LT. 0.0_rsh) THEN
                    print *,'  something goes wrong in case of l_createnewlayer==.FALSE.'
                    print *,'   dzsa=',dzsa,'  new dzs(ksmax,i,j)=',dzs(ksmax,i,j),' new poro=',poro(k,i,j)
                    print *,'   flx_w2s_loc(:)=',flx_w2s_loc(:)
                    print *,'   poro_dep=',poro_dep,' dzs_dep=',dzs_dep
                    print *,'   crel_mud(ksmax,i,j)=',crel_mud(ksmax,i,j),' poro(ksmax,i,j)=',poro(ksmax,i,j)
                  END IF
 
                  DO iv=nvpc+1,nvp
#ifdef key_Pconstitonly_insed
                    cv_sed(iv,ksmax,i,j)=0.0_rsh
#else
                    ivp_assoc=irkm_var_assoc(iv)
                    IF (ivp_assoc == 0) THEN
                      ! its own settling velocity , totaly integrated in this layer
                      cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa+flx_w2s_loc(iv))*dzsi

                    ELSE IF (ivp_assoc < imud1) THEN
                      cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa+        &
                                     flx_w2s_loc(iv)/(flx_w2s_loc(ivp_assoc)+epsi30_MUSTANG)*    &
                                     mass_tot_dep*frac_sed_dep(ivp_assoc))*dzsi  

                    ELSE 
                      cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa             &
                                     +flx_w2s_loc(iv)/(flx_w2s_loc(ivp_assoc)+epsi30_MUSTANG)*    &
                                      masdepmud*frdep(ivp_assoc))*dzsi 

                    END IF

#endif
                  ENDDO
                  k=ksmax
                  c_sedtot(k,i,j)=0.0_rsh
                  DO iv=1,nvpc
                    c_sedtot(k,i,j)=c_sedtot(k,i,j)+cv_sed(iv,k,i,j)
                  ENDDO    

                  crel_mud(k,i,j)=crel_mud_new
                  poro_mud(k,i,j)=poro_mud_new ! Actualisation de la porosite de la vase liee au melange

                  !print *,''
                  !print *,' CHARACTERISTICS OF THE MIXED SEDIMENT LAYER'
                  !print *,'ksmax=',ksmax,' check: must be equal to ',k
                  !print *,'poro(ksmax,i,j)=',poro(ksmax,i,j),' (suite au CALL MUSTANGV2_comp_poro_mixsed)'
                  !print *,'poro_mud(ksmax,i,j)=poro_mud_new=',poro_mud(ksmax,i,j)
                  !print *,'dzs(ksmax,i,j)=',dzs(k,i,j)
                  !print *,'cv_sed(:,ksmax,i,j)=',cv_sed(:,k,i,j)
                  !print *,''              

#ifdef key_MUSTANG_debug
                 IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND.   &
                       (CURRENT_TIME> t_start_debug)) THEN
                    IF (masdepmud+sommud*dzsa .GT. 0.0_rsh) THEN
                      print *,'  poro_mud(ksmax,i,j)=',(masdepmud/(masdepmud+sommud*dzsa)),' * ',&
                                poro_mud_dep,' + ',((sommud*dzsa)/(masdepmud+sommud*dzsa)),' * ',&
                                poro_muda,'=',poro_mud(k,i,j)
                    ELSE
                      print *,'  poro_mud(ksmax,i,j)=',poro_mud(k,i,j)
                    END IF
                    print *,'  mass_sed (cv_sed(iv,ksmax,i,j)*dzsa + flx_w2s_loc(iv)) = ',mass_sed
                    print *,'  mass_tot=',mass_tot
                    print *,'  frac_sed=',frac_sed
                    print *,'  poro(ksmax,i,j)=',poro(k,i,j)
                    print *,'  dzs(ksmax,i,j)=',dzs(k,i,j)
                    print *,'  verif dimi1 cv_sed iv,nvpc+1,nvp =',iv,nvpc+1,nvp
                    print *,'  verif dimi2 cv_sed ksmax,ksdmin,ksdmax =',ksmax,ksdmin,ksdmax
                    print *,'  cv_sed(:,ksmax,i,j)=',cv_sed(:,k,i,j)
                    !print *,'  cv_sed(iv,ksmax,i,j)=',cv_sed(iv,k,i,j)
                  END IF
#endif

#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
                 ! dissolved variable in pore waters and water fluxes at the interface
                 !--------------------------------------------------------------------
                 !  ==> flx_s2w and cv_sed and phieau_s2w
                 ! deposit mixing in surface layer ==> code 5
                 porewatera=poroa*dzsa
                 porewater=poro(ksmax,i,j)*dzs(ksmax,i,j)
                 CALL MUSTANGV2_eval_dissvar_IWSflux(i,j,ksmax,5, &
                         porowater1=porewatera,                   &
                         cv_sed1=cv_sed(-1:nv_adv,ksmax,i,j),     &
                         porowater_new=porewater)
#endif 
                              
               ENDIF ! test on l_createnewlayer=F (mixing)
               
               IF (l_increase_dep) THEN 

                !!!!!!!!!!!!!!!!!!!!!!!
                !!!!!! NEW LAYER !!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!
 
#ifdef key_MUSTANG_debug
                 IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND.   &
                       (CURRENT_TIME> t_start_debug)) THEN
                   print *,'  > A part of ksmax layer is added to small deposits to become suficient'
                   print *,'cmudr=',cmudr,' / cmudcr=',cmudcr
                   print *,'dzsa=',dzsa,' / dzsmax(i,j) = ',dzsmax(i,j)
                   !print *,'l_isitcohesive_dep=',l_isitcohesive_dep,' frmud_dep = ',frmud_dep,' frmudcr_dep', frmudcr_dep
                   print *,'  l_createnewlayer=',l_createnewlayer,  '***** NEW LAYER ***'
                 END IF
#endif

                 !!! Modif of deposits charac
                 IF (masdepmud+sommud*dzsmina .GT. 0.0_rsh) THEN
                   poro_mud_dep=(masdepmud/(masdepmud+sommud*dzsmina))*(1.0_rsh-(cfreshmud/ros(1))) + &
                                   ((sommud*dzsmina)/(masdepmud+sommud*dzsmina))*poro_muda
                   crel_mud_dep=(masdepmud/(masdepmud+sommud*dzsmina))*cfreshmud + &
                                    ((sommud*dzsmina)/(masdepmud+sommud*dzsmina))*crel_mud(ksmax,i,j)
                 ELSE
                   poro_mud_dep=0.0_rsh
                   crel_mud_dep=0.0_rsh
                 END IF

#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
                ! to keep in mind the amount of water coming from the ksmax layer to go into the new layer
                 porewatera=poroa*dzsmina
#endif
#if ! defined key_Pconstitonly_insed
                 ! keep in memory flx_w2s_loc before increase_dep for evaluating associated non constitutive variables
                 flx_w2s_loca(igrav1:imud2)=flx_w2s_loc(igrav1:imud2)
                 mass_tot_depa= mass_tot_dep
                 frac_sed_depa(:)=frac_sed_dep(:)
#endif
                 mass_tot_dep = 0.0_rsh
                 DO iv=igrav1,imud2
                   flx_w2s_loc(iv)=flx_w2s_loc(iv)+cv_sed(iv,ksmax,i,j)*dzsmina
                   mass_tot_dep=mass_tot_dep+flx_w2s_loc(iv)
                 END DO

                 DO iv=igrav1,imud2
                   frac_sed_dep(iv)=flx_w2s_loc(iv)/mass_tot_dep
                 END DO

                 CALL MUSTANGV2_comp_poro_mixsed(frac_sed_dep, poro_mud_dep,   &
                                            crel_mud_dep, poro_dep)
              
                 dzs_dep=mass_tot_dep/((1.0_rsh-poro_dep)*ros(1))


                 !!! Modif of ksmax layer 
                 dzs(ksmax,i,j)=dzsa-dzsmina

                 IF (dzs(ksmax,i,j) .LT. dzsmina) THEN 
                   print *,' in case of l_increase_dep=',l_increase_dep
                   print *,'  something goes wrong :'
                   print *,'  dzsa=',dzsa,' dzsmina=',dzsmina,' new dzs(ksmax,i,j)=',dzs(ksmax,i,j)
                 END IF

#ifdef key_MUSTANG_debug
                 IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND.   &
                       (CURRENT_TIME> t_start_debug)) THEN
                   IF (masdepmud+sommud*dzsa .GT. 0.0_rsh) THEN
                     print *,'  poro_mud_dep=',(masdepmud/(masdepmud+sommud*dzsmina)),' * ',  &
                             (1.0_rsh-(cfreshmud/ros(1))),' + ',((sommud*dzsmina)/(masdepmud+sommud*dzsmina)), &
                              ' * ',poro_muda,'=',poro_mud_dep
                   ELSE
                     print *,'  poro_mud_dep=',poro_mud_dep
                   END IF
                   print *,'  flx_w2s_loc(iv)=flx_w2s_loc(iv)+cv_sed(iv,ksmax,i,j)*dzsmina = ',flx_w2s_loc(:)
                   print *,'  mass_tot_dep=',mass_tot_dep
                   print *,'  frac_sed_dep=',frac_sed_dep
                   print *,'  poro_dep=',poro_dep
                   print *,'  dzs_dep=',dzs_dep
                   print *,'  dzs(ksmax,i,j)=dzsa-dzsmina=',dzs(ksmax,i,j)
                 END IF
#endif
               END IF ! test on l_increase_dep
              ENDIF !  (sedment exist) test on ksmax < or > than ksmi(i,j)
  
              ! ** if l_createnewlayer==.False. **
              !            no existing sediment or no mixing with underlying layer
              !            because of constrasting features (existing sediment already consolidated)
              !            --> both conditions lead to the creation of a new layer 
              ! -----------------------------------------------------------


              ! test  (new layer)
              IF (l_createnewlayer) THEN 

#ifdef key_MUSTANG_debug
               IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND.   &
                       (CURRENT_TIME> t_start_debug)) THEN
                 print *,''
                 print *,'  > Creation of a new layer, l_createnewlayer=',l_createnewlayer, ' because : '
                 IF (cmudr .GT. cmudcr) print *,cmudr,' (cmudr) > ',cmudcr,' (cmudcr)'
                 IF (dzsa .GE. dzsmax(i,j)) print *,dzsa,' (dzsa) >= ',dzsmax(i,j),' (dzsmax(i,j))'
                 IF (dzs_dep .GE. dzsmin) print *,dzs_dep,' (dzs_dep) >= ',dzsmin,' (dzsmin)'
               END IF
#endif

               ! actual constitution of the new layer:
               ! -------------------------------------

               IF(ksmax.EQ.ksdmax) THEN
#ifdef key_MUSTANG_debug
               IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND.   &
                       (CURRENT_TIME> t_start_debug)) THEN
                 print *,''
                 print *,'  > fusion because ksmaxx ==ksdmax and l_createnewlayer=',l_createnewlayer
                 print *,'dzs before fusion  ',dzs(ksmi(i,j):ksmax,i,j)
                 porewater=0._rsh
                 do k=ksmi(i,j),ksmax
                    porewater=porewater+dzs(k,i,j)*poro(k,i,j)
                 enddo
                 print *,'water vol tot before fusion',porewater
               END IF
#endif
                 CALL sed_MUSTANG_fusion(i,j,ksmax)
                 !print *,'CALL sed_MUSTANG_fusion(i,j,ksmax) at ',i,j
#ifdef key_MUSTANG_debug
               IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND.   &
                       (CURRENT_TIME> t_start_debug)) THEN
                 print *,'dzs after fusion ',dzs(ksmi(i,j):ksmax,i,j)
                 porewater=0._rsh
                 do k=ksmi(i,j),ksmax
                    porewater=porewater+dzs(k,i,j)*poro(k,i,j)
                 enddo
                 print *,'water vol tot after fusion',porewater
               END IF
#endif

               END IF

               ksmax=ksmax+1
               dzs(ksmax,i,j)=dzs_dep
               dzsi=1.0_rsh/dzs(ksmax,i,j)

               DO iv=igrav1,imud2
                 !cv_sed(iv,ksmax,i,j)=cv_sed_dep(iv)
                 cv_sed(iv,ksmax,i,j)=flx_w2s_loc(iv)*dzsi
               END DO

#if ! defined key_Pconstitonly_insed

               IF(l_increase_dep) THEN
                 DO iv=nvpc+1,nvp
                   ivp_assoc=irkm_var_assoc(iv)
                    IF (ivp_assoc == 0) THEN
                      cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax-1,i,j)*dzsmina+flx_w2s_loc(iv))*dzsi
                    ELSE IF (ivp_assoc < imud1) THEN
                     cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax-1,i,j)*dzsmina+        &
                                     flx_w2s_loc(iv)/(flx_w2s_loca(ivp_assoc)+epsi30_MUSTANG)*    &
                                     mass_tot_depa*frac_sed_depa(ivp_assoc))*dzsi  

                    ELSE 
                      cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax-1,i,j)*dzsmina             &
                                     +flx_w2s_loc(iv)/(flx_w2s_loca(ivp_assoc)+epsi30_MUSTANG)*    &
                                      masdepmud*frdep(ivp_assoc))*dzsi 

                    END IF
                 ENDDO
               ELSE
                 DO iv=nvpc+1,nvp
                   ivp_assoc=irkm_var_assoc(iv)
                   IF (ivp_assoc .NE. 0) THEN
                     cv_sed(iv,ksmax,i,j)= flx_w2s_loc(iv)/(flx_w2s_loc(ivp_assoc)+epsi30_MUSTANG)*    &
                                        cv_sed(ivp_assoc,ksmax,i,j)
                   ELSE
                     cv_sed(iv,ksmax,i,j)=flx_w2s_loc(iv)/dzs_dep
                   END IF 
                 ENDDO
               ENDIF
#endif

               ! computing total concentration and porosity of the surf.layer:
               ! -------------------------------------------------------------
               k=ksmax
               c_sedtot(k,i,j)=0.0_rsh
               DO iv=1,nvpc
                 c_sedtot(k,i,j)=c_sedtot(k,i,j)+cv_sed(iv,k,i,j)
               ENDDO
               crel_mud(k,i,j)=crel_mud_dep
               poro(k,i,j)=poro_dep
               poro_mud(k,i,j)=poro_mud_dep

#ifdef key_MUSTANG_debug
               IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND.   &
                       (CURRENT_TIME> t_start_debug)) THEN
                 print *,'  > Charac of the new sediment layer'
                 print *,'  ksmax=',ksmax
                 print *,'  poro(ksmax,i,j)=',poro(k,i,j)
                 print *,'  poro_mud(ksmax,i,j)=',poro_mud(k,i,j)
                 print *,'  dzs(ksmax,i,j)=',dzs(k,i,j)
                 print *,'  cv_sed(:,ksmax,i,j)=',cv_sed(:,k,i,j)
                 print *,'  c_sedtot(ksmax,i,j)=',c_sedtot(k,i,j)
               END IF
#endif 

#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
                 ! dissolved variable in pore waters and water fluxes at the interface
                 !--------------------------------------------------------------------
                 !  ==> flx_s2w and cv_sed and phieau_s2w
               porewater=poro(ksmax,i,j)*dzs(ksmax,i,j)
               IF (l_increase_dep) THEN
                 ! create a new surface layer with deposit and the old ksmax layer==>  7  (input of pore water from water column)
                 CALL MUSTANGV2_eval_dissvar_IWSflux(i,j,ksmax,7, &
                         porowater1=porewatera,                    &
                         cv_sed1=cv_sed(-1:nv_adv,ksmax-1,i,j),   &
                         porowater2=porewaterdep,                 &
                         porowater_new=porewater)

               ELSE
                 ! create a new surface layer with deposit ==> code 6 (input of pore water from water column)
                 CALL MUSTANGV2_eval_dissvar_IWSflux(i,j,ksmax,6, &
                         porowater_new=porewater)
               ENDIF
#endif         

                ! to avoid increasing the thickness of the surface layer 
                IF(ksmax .LT. ksdmax .AND. ksmax > ksmi(i,j)) THEN
                    IF(dzs(ksmax,i,j) > dzsmax(i,j) + 5.0_rsh* dzsmin) THEN
#ifdef key_MUSTANG_debug
                       IF ( l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
                         print *,'    SPLIT SURFACE LAYER BECAUSE dzs > dzsmax '
                         print *,ksmax,'  layers become', ksmax+1, 'layers'
                       END IF
#endif
                     dzs(ksmax+1,i,j)=MIN(dzs(ksmax,i,j)-dzsmax(i,j),dzsmax(i,j))
                     dzs(ksmax,i,j)=dzs(ksmax,i,j)-dzs(ksmax+1,i,j)
                     poro(ksmax+1,i,j)=poro(ksmax,i,j)
                     poro_mud(ksmax+1,i,j)=poro_mud(ksmax,i,j)
                     crel_mud(ksmax+1,i,j)=crel_mud(ksmax,i,j)
                     cv_sed(:,ksmax+1,i,j)=cv_sed(:,ksmax,i,j)
                     c_sedtot(ksmax+1,i,j)=c_sedtot(ksmax,i,j)
                     ksmax=ksmax+1
                   ENDIF
                ENDIF

              ENDIF  !  test  (new layer)
 
            ELSE   ! test  (fludep < 0 no deposition)

             ! no deposition
             ! ****************
 
               ! computing flx_w2s for dissolved var., to be used in the next time step:
               ! -----------------------------------------------------------------------
#if ! defined key_noTSdiss_insed
               DO iv=nvp+1,nv_adv
                 flx_w2s(iv,i,j)=0.0_rsh
               ENDDO
               flx_w2s(-1,i,j)=0.0_rsh
               flx_w2s(0,i,j)=0.0_rsh
#endif

               ! accounting for possible deposition of non constitutive particulate variable:
               ! ----------------------------------------------------------------------------
               IF(ksmax.GE.ksmi(i,j))THEN
                 DO iv=nvpc+1,nvp
                   cv_sed(iv,ksmax,i,j)=cv_sed(iv,ksmax,i,j)+flx_w2s_loc(iv)/dzs(ksmax,i,j)
                 ENDDO
               ELSE
                 ! no sediment
                 DO iv=nvpc+1,nvp
#ifdef key_MARS
#if defined key_siggen || defined key_gencoord
                  cw_bottom_MUSTANG(iv,i,j)=cw_bottom_MUSTANG(iv,i,j)                                        &
                               +flx_w2s_loc(iv)/(epn_bottom_MUSTANG(i,j))
#else
                  cw_bottom_MUSTANG(iv,i,j)=cw_bottom_MUSTANG(iv,i,j)                                        &
                               +flx_w2s_loc(iv)/(epn_bottom_MUSTANG(i,j))
#endif
#else
                ! To Program  **TODO**
#endif
                  iexchge_MPI_cvwat=1
                 ENDDO
               ENDIF


            ENDIF  ! test  (fludep < 0 no deposition)

            ! updating ksma ( ksmax can be modified in the routine)
            ksma(i,j)=ksmax
            dzs(ksma(i,j)+1:ksdmax,i,j)=0.0_rsh
            poro(ksma(i,j)+1:ksdmax,i,j)=0.0_rsh

#ifdef key_MUSTANG_specif_outputs
            varspecif3Dk_save(1,:,i,j)=0.0_rsh       ! poro_save
            DO k=ksmi(i,j),ksmax
              !!! BM: check  porosity (output)
              varspecif3Dk_save(1,k,i,j)=poro(k,i,j)  ! poro_save
            ENDDO
            varspecif2D_save(2,i,j)=dzs(ksmax,i,j)  ! dzs at sediment surface
#endif
#ifdef key_MUSTANG_debug
               IF (l_debug_effdep .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
                 print *,'  > fin effdep',i,j
                 print *,'  t=',CURRENT_TIME, 'ksmax=',ksmax
                 print *,'  dzs(ksmax-3:ksmax,i,j)=',dzs(ksmax-3:ksmax,i,j)
                 print *,'  cv_sed(:,ksmax-3,i,j)=',cv_sed(:,ksmax-3,i,j)
                 print *,'  cv_sed(:,ksmax-2,i,j)=',cv_sed(:,ksmax-2,i,j)
                 print *,'  cv_sed(:,ksmax-1,i,j)=',cv_sed(:,ksmax-1,i,j)
                 print *,'  cv_sed(:,ksmax,i,j)=',cv_sed(:,ksmax,i,j)
               END IF
#endif 
        
          END IF ! test on htot
        END DO  ! loop on i
    END DO    ! loop on j

  END SUBROUTINE sed_MUSTANG_effdep
   ! end version V2

   !!==============================================================================

#else
!! version V1
  SUBROUTINE sed_MUSTANG_effdep(ifirst, ilast, jfirst, jlast, iexchge_MPI_cvwat)
    
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_effdep version V1 ***
   !&E
   !&E ** Purpose : computes effective deposition
   !&E
   !&E ** Description : 
   !&E
   !&E       arguments IN :
   !&E          loops  :ifirst,ilast,jfirst,jlast
   !&E          
   !&E       variables OUT :
   !&E           phieau_s2w
   !&E           flx_w2s
   !&E
   !&E ** Called by :  MUSTANG_deposition
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used

   !! * Arguments
   INTEGER, INTENT(IN)  :: ifirst,ilast,jfirst,jlast
   INTEGER, INTENT(INOUT)  :: iexchge_MPI_cvwat

   !! * Local declarations
   REAL(KIND=rsh),DIMENSION(nvp)  :: flx_w2s_loc
   REAL(KIND=rsh),DIMENSION(nvpc) :: frdep          
   INTEGER                        :: i,j,k,iv,ksmin,ksmax,kl,ivp_assoc,isplit
   REAL(KIND=rsh) :: fludep,ddzs,ddzs1,ddzs2,dzsa,dzsisu,dzsi,dzsaici,      &
                     dzsgrv,dzssan,dzsmud,dzsnew,dzsa2,dzsa3,ddzs3,ddzs4,ddzs5,ddzs6,    &
                     poroa,somalp,somala,cdepiv,porosi,dzpoi,porodep,fluxdissous,        &
                     ddzsici,voldepgrv,voldepsan,masdepmud,dzskmanew,porewatera,         &
                     sommud,cvolinigrv,cvolinisan,cvolp,cmudr,porewater,                 &
                     dvolgrv,dvolsan,dmasmud,dmasmudgrav,dmasmudsand,cordepflu,          &
                     cordepfluw,cordepflue,cordepflus,cordepflun,                        &
                     dmasmud1,dmasmud2,dmasmud3,dmasmud4

   !!----------------------------------------------------------------------
   !! * Executable part

   ddzsici=-1000.0_rsh
   flx_w2s_loc(:)=0.0_rsh
   frdep(:)=0.0_rsh        
   iexchge_MPI_cvwat=0         

      DO j=jfirst,jlast
        DO i=ifirst,ilast
          IF(htot(i,j) > h0fond) THEN

            ksmax=ksma(i,j)
            ksmin=ksmi(i,j)

            ! updating effective deposition :  (flx_w2s) is implicit in the vertical advection scheme
            ! hence multiplied here by the newly computed concentration at each sub time step and cumulated
            DO iv=igrav2+1,nv_use   ! nv_use=nvp (or nvpc if key_Pconstitonly_insed)
              flx_w2s_loc(iv)=MAX(0.0_rsh,flx_w2s_sum(iv,i,j)+flx_w2s_corin(iv,i,j)      &
                           +flx_w2s_corip1(iv,i-1,j)+flx_w2s_corim1(iv,i+1,j)         &
                           +flx_w2s_corjp1(iv,i,j-1)+flx_w2s_corjm1(iv,i,j+1))

              flx_w2s_sum(iv,i,j)=0.0_rsh

              ! MF /= 1 only if l_morphocoupl
              flx_w2s_loc(iv) = MF * flx_w2s_loc(iv)
#ifdef key_MUSTANG_specif_outputs
              !flx_w2s_save
             varspecif3Dnv_save(3,iv,i,j)=varspecif3Dnv_save(3,iv,i,j)+flx_w2s_loc(iv)
#endif
            ENDDO

            fludep=0.0_rsh
            DO iv=1,nvpc
             fludep=fludep+flx_w2s_loc(iv)
            ENDDO

            voldepgrv=0.0_rsh
            voldepsan=0.0_rsh
            masdepmud=0.0_rsh
            dmasmud1 = 0.0_rsh
            dmasmud2 = 0.0_rsh
            dmasmud3 = 0.0_rsh
            dmasmud4 = 0.0_rsh
              
            IF(fludep.GT.0.0_rsh)THEN

             ! case 1:  there is deposition
             ! ****************************

             ! characterising deposition:
             DO iv=igrav1,igrav2
              voldepgrv=voldepgrv+flx_w2s_loc(iv)/ros(iv)
             ENDDO
             DO iv=igrav1,igrav2
              frdep(iv)=flx_w2s_loc(iv)/(voldepgrv+epsi30_MUSTANG)
             ENDDO
             DO iv=isand1,isand2
              voldepsan=voldepsan+flx_w2s_loc(iv)/ros(iv)
             ENDDO
             DO iv=isand1,isand2
              frdep(iv)=flx_w2s_loc(iv)/(voldepsan+epsi30_MUSTANG)
             ENDDO
             DO iv=imud1,imud2
              masdepmud=masdepmud+flx_w2s_loc(iv)
             ENDDO
             DO iv=imud1,imud2
              frdep(iv)=flx_w2s_loc(iv)/(masdepmud+epsi30_MUSTANG)
             ENDDO

             IF(ksmax.GE.ksmi(i,j))THEN
          
              ! case 1.1: sediment already exists
              ! ---------------------------------
              ! characterising previous surficial sediment:

              sommud=0.0_rsh      
              DO iv=imud1,imud2
                sommud=sommud+cv_sed(iv,ksmax,i,j)
              ENDDO
              cvolinigrv=0.0_rsh      
              DO iv=igrav1,igrav2
                cvolinigrv=cvolinigrv+cv_sed(iv,ksmax,i,j)/ros(iv)
              ENDDO
              cvolinisan=0.0_rsh
              DO iv=isand1,isand2
                cvolinisan=cvolinisan+cv_sed(iv,ksmax,i,j)/ros(iv)
              ENDDO
              cmudr=sommud/(1.0_rsh-cvolinisan-cvolinigrv)
            
              IF(cmudr.LE.cmudcr)THEN
            
                ! case 1.1.1: surficial sediment is not consolidated
                ! --------------------------------------------------
                ! mixing deposits with upper layer, until completion:

                dzsa=dzs(ksmax,i,j)

                ! first, gravels deposition :
                ! ---------------------------
                ! dvolgrv represents the volume of gravel that is inserted into
                ! the existing sediment, by densifying it
                ! it will remain to be deposited at the volume concentration
                ! cvolmaxsort, whose ddzs2 in the existing layer, ddzs1 in
                ! a new layer
                dvolgrv=MIN(voldepgrv,dzsa*(cvolmaxsort-cvolinigrv),dzsa*(cvolmaxmel-cvolinigrv-cvolinisan))
                dvolgrv=MAX(dvolgrv,0.0_rsh)
                ddzs=(voldepgrv-dvolgrv)/cvolmaxsort
                ddzsici=ddzs
                dzsaici=dzsa
                IF((dzsa+ddzs).LE.dzsmax(i,j))THEN
                  ddzs1=0.0_rsh
                ELSEIF(dzsa.GE.dzsmax(i,j))THEN
                  ddzs1=ddzs
                ELSE
                  ddzs1=dzsa+ddzs-dzsmax(i,j)
                ENDIF              
                !ddzs1=max(0.0_rsh,dzsa+ddzs-dzsmax(i,j))t
                ddzs2=ddzs-ddzs1
                dzsa2=dzsa+ddzs2
                voldepgrv=ddzs1*cvolmaxsort

                ! sand deposition :
                ! -----------------
                ! modification du calcul de dvolsan (9 janvier 2012)
                !dvolsan=MIN(voldepsan,dzsa2*cvolmaxsort-dzsa*cvolinisan, &
                !                      dzsa2*cvolmaxmel-dzsa*(cvolinigrv+cvolinisan)-ddzs2*cvolmaxsort)
                ! ajout dvolgrv
                dvolsan=MIN(voldepsan,dzsa2*cvolmaxsort-dzsa*cvolinisan, &
                                      dzsa2*cvolmaxmel-dzsa*(cvolinigrv+cvolinisan)-ddzs2*cvolmaxsort)
                dvolsan=MAX(dvolsan,0.0_rsh)
                ! une partie des sables dvolsab (avec repartition frdep) 
                ! vient densifier l epaisseur dzsa2,l autre partie vient creer 
                ! une surepaisseur ddzs4 (repartition frdep) et
                ! contribuer a une nouvelle couche pour ddzs3 (repartition frdep)
                 ! some of the sands dvolsab (with distribution frdep)
                 ! just densify the thickness dzsa2, the other part just create
                 ! a large amount of ddzs4 (distribution frdep) and
                 ! contribute to a new layer for ddzs3 (distribution frdep)
                ddzs=(voldepsan-dvolsan)/cvolmaxsort
                ! suppression de la ligne suivante (9 janvier 2012)
                !voldepsan=0.0_rsh
                IF((dzsa2+ddzs).LE.dzsmax(i,j))THEN
                  ddzs3=0.0_rsh
                ELSEIF(dzsa2.GE.dzsmax(i,j))THEN
                  ddzs3=ddzs
                ELSE
                  ddzs3=dzsa2+ddzs-dzsmax(i,j)
                ENDIF              
                !ddzs3=max(0.0_rsh,dzsa2+ddzs-dzsmax(i,j))
                ddzs4=ddzs-ddzs3
                dzsa3=dzsa2+ddzs4
                voldepsan=ddzs3*cvolmaxsort

                ! mud deposition :
                ! ----------------
                ! une partie des vases dmasmud (avec repartition frdep) 
                ! vient densifier l epaisseur dzsa3, l autre partie vient creer 
                ! une surepaisseur ddzs6 (repartition frdep) et
                ! contribuer a une nouvelle couche pour ddzs5 (repartition frdep)
                 ! some of the dmasmud muds (with frdep distribution)
                 ! just densify the thickness dzsa3, the other part just create
                 ! ddzs6 overdimension (distribution frdep) and
                 ! contribute to a new layer for ddzs5 (frdep distribution)

!!! B.Mengual (17/09/2015): 
!!! Le depot de la vase se fait en commencant a melanger depuis la surface
!!! dans le cas ou un excedant de sables ou de graviers conduira forcement
!!! a la creation dune nouvelle couche
!!! BUT : ne pas pieger de la vase par melange dans la couche ksmax-1 suite
!!!       a la creation dune nouvelle couche


                IF ((voldepgrv+voldepsan) .GT. 0.0_rsh) THEN

                  !!! Part de vase qui va etre melangee avec la future couche
                  ! On enleve les cas a cvolmaxmel pour etre coherent avec les ddzs1 et ddzs3
                  ! calcules avec cvolmaxsort
                   !!! Part of the mud that will be mixed with the future layer
                   ! We remove the cases with cvolmaxmel to be coherent with the ddzs1 and ddzs3
                   ! calculated with cvolmaxsort
                  dmasmud1 = MIN(masdepmud,(ddzs1+ddzs3)*cfreshmud*(1-cvolmaxsort))
                  !!! Part de vase melangee avec les depots en sables / graviers effectues 
                  !!! au cours du pas de temps dans la couche ksmax
                  ! Meme commentaire
                   !!! Part of mud mixed with sand / gravel deposits
                   !!! during the time step in the ksmax layer
                   ! Even comment
                  dmasmud2 = MIN(masdepmud-dmasmud1,(dzsa3-dzsa)*cfreshmud*(1-cvolmaxsort))
                  !!! Part de vase melangee avec le sediment present initialement dans la couche ksmax
                   !!! Part of mud mixed with the sediment initially present in the ksmax layer
                  IF((cvolinigrv+cvolinisan) .GT. 0.0_rsh) THEN
                    dmasmud3 = MIN(masdepmud-dmasmud1-dmasmud2,dzsa*(cfreshmud-cmudr)*(1-cvolinigrv-cvolinisan))
                  END IF
                  dmasmud3 = MAX(0.0_rsh,dmasmud3)
                  dmasmud = dmasmud2 + dmasmud3
                  !!! Calcul de l excedant de vase a ajouter dans la futur couche une fois que tout a ete melange
                   !!! Calculation of the excess of mud to add in the future layer once everything has been mixed
                  dmasmud4 = masdepmud - dmasmud1 - dmasmud
                  !!! Calcul de l epaisseur totale de vase qui sera attribuee a la future couche
                  !!! Attention a ne pas confondre ddzs5 et dzsmud car une partie de ddzs5 sera melangee a lexcedant 
                  !!! de sable gravier avant daccroitre lepaisseur de la couche
                   !!! Calculation of the total mud thickness to be attributed to the future layer
                   !!! Be careful not to confuse ddzs5 and dzsmud because a part of ddzs5 will be mixed with the excess
                   !!! of sand gravel before increasing the thickness of the layer
                  ddzs5 = (dmasmud4/cfreshmud) + (dmasmud1/cfreshmud)
                  ddzs6 = 0.0_rsh

                  !IF (dmasmud1+dmasmud2+dmasmud3+dmasmud4 .LT. masdepmud) THEN
                  !  print*,'*** BUG ***'
                  !  print*,'La somme des dmasmud est inferieure a masdepmud au point (i,j) : ',i,' ',j
                  !  print*,'dmasmud1,dmasmud2,dmasmud3,dmasmud4,masdepmud : ',dmasmud1,dmasmud2,dmasmud3,dmasmud4,masdepmud
                  !END IF

                  !IF (dmasmud1+dmasmud2+dmasmud3+dmasmud4 .GT. masdepmud) THEN
                  !  print*,'*** BUG ***'
                  !  print*,'La somme des dmasmud est superieure a masdepmud au point (i,j) : ',i,' ',j
                  !  print*,'dmasmud1,dmasmud2,dmasmud3,dmasmud4,masdepmud : ',dmasmud1,dmasmud2,dmasmud3,dmasmud4,masdepmud
                  !END IF


                ELSE
                  dmasmud2 = MIN(masdepmud,(dzsa3-dzsa)*cfreshmud*(1-cvolmaxsort))
                  IF((cvolinigrv+cvolinisan) .GT. 0.0_rsh) THEN
                    dmasmud3 = MIN(masdepmud-dmasmud2,dzsa*(cfreshmud-cmudr)*(1-cvolinigrv-cvolinisan))
                  END IF
                  dmasmud3 = MAX(0.0_rsh,dmasmud3)
                  dmasmud = dmasmud2 + dmasmud3

!                 ! une partie des vases dmasmud (avec repartition frdep) 
!                 ! vient densifier l epaisseur dzsa3, l autre partie vient creer 
!                 ! une surepaisseur ddzs6 (repartition frdep) et
!                 ! contribuer a une nouvelle couche pour ddzs5 (repartition frdep)
!                   ! some of the  mud dmasmud (with frdep distribution)
!                   ! just densify the thickness dzsa3, the other part just create
!                   ! ddzs6 overdimension (distribution frdep) and
!                   ! contribute to a new layer for ddzs5 (frdep distribution)
!                 IF(cvolinigrv >0.0_rsh .and. cvolinisan > 0.0_rsh) THEN
!                   ! mixture sand+gravel  
!                   dmasmud=MIN(masdepmud,(dzsa3-dzsa)*cfreshmud*(1.0_rsh-cvolmaxmel))  ! masse de vase qui s insere dans un melange
!                 ELSE IF ((cvolinisan > 0.0_rsh .AND. isand2>isand1) .OR. (cvolinigrv > 0.0_rsh .AND. igrav2>igrav1) ) THEN
!                   ! mixture sand or gravel
!                  dmasmud=MIN(masdepmud,(dzsa3-dzsa)*cfreshmud*(1.0_rsh-cvolmaxmel))  ! masse de vase qui s insere dans un melange
!                 ELSE
!                   ! one type of sediment 
!                   dmasmud=MIN(masdepmud,(dzsa3-dzsa)*cfreshmud*(1.0_rsh-cvolmaxsort))  ! masse de vase qui s insere dans le gravier ou sable (1-cvolmaxsort)
!                 ENDIF

                  ddzs=(masdepmud-dmasmud)/cfreshmud
                  IF((dzsa3+ddzs).LE.dzsmax(i,j))THEN
                    ddzs5=0.0_rsh
                  ELSEIF(dzsa3.GE.dzsmax(i,j))THEN
                    ddzs5=ddzs
                  ELSE
                    ddzs5=dzsa3+ddzs-dzsmax(i,j)
                  ENDIF
                  !ddzs5=max(0.0_rsh,dzsa3+ddzs-dzsmax(i,j))
                  ddzs6=ddzs-ddzs5

                END IF

                masdepmud=ddzs5*cfreshmud

                dzs(ksmax,i,j)=dzsa3+ddzs6
                dzsi=1.0_rsh/dzs(ksmax,i,j)
                DO iv=igrav1,igrav2
                  cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa+(dvolgrv    &
                                     +ddzs2*cvolmaxsort)*frdep(iv))*dzsi
                ENDDO
                DO iv=isand1,isand2
                  cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa+(dvolsan    &
                                     +ddzs4*cvolmaxsort)*frdep(iv))*dzsi
                ENDDO
                DO iv=imud1,imud2
                  cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa+(dmasmud    &
                                     +ddzs6*cfreshmud)*frdep(iv))*dzsi
                ENDDO

                DO iv=nvpc+1,nvp
#ifdef key_Pconstitonly_insed
                  cv_sed(iv,ksmax,i,j)=0.0_rsh
#else
                  ivp_assoc=irkm_var_assoc(iv)
                  IF (ivp_assoc == 0) THEN
                    ! its own settling velocity , totaly integrated in this layer
                    cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa+flx_w2s_loc(iv))/max(dzsmin,dzs(ksmax,i,j))
                    ! cancel the flx_w2s_loc (iv) so as not to add it to a new hypothetical layer above
                    flx_w2s_loc(iv)=0.0_rsh
                  ELSE IF (ivp_assoc < isand1) THEN
                    cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa+        &
                                     flx_w2s_loc(iv)/(flx_w2s_loc(ivp_assoc)+epsi30_MUSTANG)*    &
                                     (dvolgrv+ddzs2*cvolmaxsort)*frdep(ivp_assoc))*dzsi
                  ELSE IF (ivp_assoc < imud1) THEN
                    cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa+         &
                                      flx_w2s_loc(iv)/(flx_w2s_loc(ivp_assoc)+epsi30_MUSTANG)*   &
                                    (dvolsan+ddzs4*cvolmaxsort)*frdep(ivp_assoc))*dzsi
                  ELSE 
                    cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa             &
                                     +flx_w2s_loc(iv)/(flx_w2s_loc(ivp_assoc)+epsi30_MUSTANG)*    &
                                     (dmasmud+ddzs6*cfreshmud)*frdep(ivp_assoc))*dzsi
                  END IF
#endif
                ENDDO
                k=ksmax
                c_sedtot(k,i,j)=0.0_rsh
                cvolp=0.0_rsh
                DO iv=1,nvpc
                  c_sedtot(k,i,j)=c_sedtot(k,i,j)+cv_sed(iv,k,i,j)
                  cvolp=cvolp+cv_sed(iv,k,i,j)/ros(iv)
                ENDDO    
                poroa=poro(k,i,j)
                poro(k,i,j)=1.0_rsh-cvolp
          
                !  dissolved variable in pore waters
                !------------------------------------------
                porewatera=poroa*dzsa
                porewater=poro(k,i,j)*dzs(ksmax,i,j)
                IF(porewater > porewatera) THEN
                   !input of pore water from water column
                   dzpoi=dzsi/poro(ksmax,i,j)
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
                   DO iv=nvp+1,nv_adv
                     fluxdissous=(porewater-porewatera)*cw_bottom_MUSTANG(iv,i,j)  
                     ! flux dissous ou flx_w2s  sera divise par dt qd pris en compte dans erosion, pour etre en quantite/m2/s
                      ! dissolved flux or flx_w2s will be divided by dt when taken into account in erosion, to be in quantity/m2/s                     cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*porewatera +fluxdissous)*dzpoi
                     cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*porewatera +fluxdissous)*dzpoi
                     flx_w2s(iv,i,j)=flx_w2s(iv,i,j)+fluxdissous
                   ENDDO
#endif
                   ! temperature and salinity
                   fluxdissous=(porewater-porewatera)*temp_bottom_MUSTANG(i,j)  
                   cv_sed(-1,ksmax,i,j)=(cv_sed(-1,ksmax,i,j)*porewatera +fluxdissous)*dzpoi
                   flx_w2s(-1,i,j)=flx_w2s(-1,i,j)+fluxdissous
                   fluxdissous=(porewater-porewatera)*sal_bottom_MUSTANG(i,j)  
                   cv_sed(0,ksmax,i,j)=(cv_sed(0,ksmax,i,j)*porewatera +fluxdissous)*dzpoi
                   flx_w2s(0,i,j)=flx_w2s(0,i,j)+fluxdissous
#endif
#if ! defined key_nofluxwat_IWS
                   ! water flux 
                   phieau_s2w(i,j)=phieau_s2w(i,j)+REAL((porewatera-porewater)*CELL_SURF(i,j),rlg)
#endif
                ELSE
                   !output of pore-water to water column ; cvsed are inchanged
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
                   DO iv=nvp+1,nv_adv
                     flx_w2s(iv,i,j)=flx_w2s(iv,i,j)+(porewater-porewatera)*cv_sed(iv,k,i,j)  
                     ! flx_w2s will be divided by dt when taken into account in erosion, to be in quantity/m2/s                     
                     !cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*porewatera +fluxdissous)*dzpoi
                   ENDDO
#endif
                   ! temperature and salinity
                   flx_w2s(-1,i,j)=flx_w2s(-1,i,j)+(porewater-porewatera)*cv_sed(-1,k,i,j)  
                   flx_w2s(0,i,j)=flx_w2s(0,i,j)+(porewater-porewatera)*cv_sed(0,k,i,j)
#endif
#if ! defined key_nofluxwat_IWS
                   ! water flux 
                   phieau_s2w(i,j)=phieau_s2w(i,j)+REAL((porewatera-porewater)*CELL_SURF(i,j),rlg)
#endif
                ENDIF              
              ENDIF
             ENDIF
  
             ! case 1.2 or following 1.1: computation of remaining deposits
             ! -----------------------------------------------------------

             dzsgrv=0.0_rsh
             dzssan=0.0_rsh
             dzsmud=0.0_rsh
             ! gravel deposition:
             ! ------------------
             IF(voldepgrv.GE.0.0_rsh)THEN
              dzsgrv=voldepgrv/cvolmaxsort
             ENDIF

             ! sand deposition:
             ! ----------------
             IF(voldepsan.GE.0.0_rsh)THEN
              dvolsan=MIN(voldepsan,dzsgrv*(cvolmaxmel-cvolmaxsort))
              dzssan=(voldepsan-dvolsan)/cvolmaxsort
             ENDIF

             ! mud deposition:
             ! ---------------
             IF(masdepmud.GE.0.0_rsh)THEN
               dmasmud=MIN(masdepmud,(dzsgrv+dzssan)*cfreshmud*(1.0_rsh-cvolmaxmel))
               dzsmud=(masdepmud-dmasmud)/cfreshmud
             ENDIF
          
             dzsnew=dzsgrv+dzssan+dzsmud

             IF(dzsnew.GT.0.0_rsh)THEN
          
              ! case 1.2 or following 1.1: possible creation of a new layer
              ! ----------------------------------------------------------

              !we choose to do the same test with or whithout consolidation. 
              !We will manage very small dzs in consolidation module
              IF(ksmax == 0 .OR. l_consolid) THEN
                dzskmanew=0.0_rsh
              ELSE
                dzskmanew=dzs(ksmax,i,j)+dzsnew
              ENDIF

#ifdef key_BLOOM_insed
              ! change test because dzsnew ca not be to small if there are biological variables
              ! if deposit of biological variables and no contitituve variables ==> /dzs => cvsed very big
              dzskmanew=0.0_rsh
#endif 

              ! creation only if dzsnew>dzsmin or if overtake dzsmax
              IF(((dzsnew >= dzsmin .OR. dzskmanew > dzsmax(i,j)) .AND. ksdmax > 1) .OR. ksmax == 0)THEN

                ! actual constitution of the new layer:
                ! -------------------------------------
                IF(ksmax.EQ.ksdmax)CALL sed_MUSTANG_fusion(i,j,ksmax)
                dzsa=0.0_rsh
                ksmax=ksmax+1
                dzs(ksmax,i,j)=dzsnew
                dzsi=1.0_rsh/dzs(ksmax,i,j)
                DO iv=igrav1,igrav2
                  cv_sed(iv,ksmax,i,j)=voldepgrv*frdep(iv)*dzsi
                ENDDO
                DO iv=isand1,isand2
                  cv_sed(iv,ksmax,i,j)=voldepsan*frdep(iv)*dzsi
                ENDDO
                DO iv=imud1,imud2
                  cv_sed(iv,ksmax,i,j)=masdepmud*frdep(iv)*dzsi
                ENDDO

#if ! defined key_Pconstitonly_insed
                DO iv=nvpc+1,nvp
                  ivp_assoc=irkm_var_assoc(iv)
                  IF (ivp_assoc .NE. 0) THEN
                    cv_sed(iv,ksmax,i,j)= flx_w2s_loc(iv)/(flx_w2s_loc(ivp_assoc)+epsi30_MUSTANG)*    &
                                        cv_sed(ivp_assoc,ksmax,i,j)
                  ELSE IF(dzsnew > dzsmin) THEN
                   cv_sed(iv,ksmax,i,j)=flx_w2s_loc(iv)*dzsi
                  END IF 
                ENDDO
#endif
                ! computing total concentration and porosity of the surf.layer:
                ! -------------------------------------------------------------
                k=ksmax
                c_sedtot(k,i,j)=0.0_rsh
                cvolp=0.0_rsh
                DO iv=1,nvpc
                  c_sedtot(k,i,j)=c_sedtot(k,i,j)+cv_sed(iv,k,i,j)
                  cvolp=cvolp+cv_sed(iv,k,i,j)/ros(iv)
                ENDDO

                poro(k,i,j)=1.0_rsh-cvolp
              
                !  dissolved variable in pore waters
                !------------------------------------------
                porewater=poro(k,i,j)*dzs(ksmax,i,j)
                !input of pore water from water column
                dzpoi=dzsi/poro(ksmax,i,j)
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
                DO iv=nvp+1,nv_adv
                  cv_sed(iv,ksmax,i,j)=cw_bottom_MUSTANG(iv,i,j)
                  flx_w2s(iv,i,j)=flx_w2s(iv,i,j)+porewater*cw_bottom_MUSTANG(iv,i,j)   
                  ! flx_w2s will be divided by dt when taken into account in erosion, to be in quantity/m2/s                    
                  ! cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*porewatera +fluxdissous)*dzpoi
                ENDDO
#endif
                cv_sed(-1,ksmax,i,j)=temp_bottom_MUSTANG(i,j)
                flx_w2s(-1,i,j)=flx_w2s(-1,i,j)+porewater*temp_bottom_MUSTANG(i,j) 
                cv_sed(0,ksmax,i,j)=sal_bottom_MUSTANG(i,j)
                flx_w2s(0,i,j)=flx_w2s(0,i,j)+porewater*sal_bottom_MUSTANG(i,j) 
#endif
#if ! defined key_nofluxwat_IWS
                ! water flux  
                phieau_s2w(i,j)=phieau_s2w(i,j)-REAL(porewater*CELL_SURF(i,j),rlg)
#endif                
              ELSE
            
                ! as the new layer is too small, back to an increase of the surf. layer:
                ! ----------------------------------------------------------------------
                dzsa=dzs(ksmax,i,j)
                dzs(ksmax,i,j)=dzs(ksmax,i,j)+dzsnew
                dzsi=1.0_rsh/dzs(ksmax,i,j)
                DO iv=igrav1,igrav2
                  cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa                 &
                                  +voldepgrv*frdep(iv))*dzsi
                ENDDO
                DO iv=isand1,isand2
                  cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa                 &
                                 +voldepsan*frdep(iv))*dzsi
                ENDDO
                DO iv=imud1,imud2
                  cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa                 &
                                  +masdepmud*frdep(iv))*dzsi
                ENDDO

#if ! defined key_Pconstitonly_insed
                DO iv=nvpc+1,nvp
                  ivp_assoc=irkm_var_assoc(iv)
                  IF (ivp_assoc == 0) THEN
                   ! its own settling velocity 
                   cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa+                    &
                             flx_w2s_loc(iv))*dzsi

                  ELSE IF (ivp_assoc < isand1) THEN
                   cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa+                    &
                                       flx_w2s_loc(iv)/(flx_w2s_loc(ivp_assoc)+epsi30_MUSTANG)*  &
                                       voldepgrv*frdep(ivp_assoc))*dzsi
                  ELSE IF (ivp_assoc < imud1) THEN
                   cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa+                    &
                                     flx_w2s_loc(iv)/(flx_w2s_loc(ivp_assoc)+epsi30_MUSTANG)*   &
                                     voldepsan*frdep(ivp_assoc))*dzsi
                  ELSE 
                   cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*dzsa+                    &
                                      flx_w2s_loc(iv)/(flx_w2s_loc(ivp_assoc)+epsi30_MUSTANG)*   &
                                      masdepmud*frdep(ivp_assoc))*dzsi
                  END IF 
                ENDDO
#endif
            
                ! computing total concentration and porosity of the surf.layer:
                ! -------------------------------------------------------------
                k=ksmax
                c_sedtot(k,i,j)=0.0_rsh
                cvolp=0.0_rsh
                DO iv=1,nvpc
                  c_sedtot(k,i,j)=c_sedtot(k,i,j)+cv_sed(iv,k,i,j)
                  cvolp=cvolp+cv_sed(iv,k,i,j)/ros(iv)
                ENDDO
                poroa=poro(k,i,j)
                poro(k,i,j)=1.0_rsh-cvolp
 
                !  dissolved variable in pore waters
                !------------------------------------------
                porewatera=poroa*dzsa
                porewater=poro(k,i,j)*dzs(ksmax,i,j)
                IF(porewater > porewatera) THEN
                  !input of pore water from water column
                  dzpoi=dzsi/poro(ksmax,i,j)
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
                  DO iv=nvp+1,nv_adv
                    fluxdissous=(porewater-porewatera)*cw_bottom_MUSTANG(iv,i,j)  
                    ! will be divided by dt when taken into account in erosion, to be in quantity/m2/s
                    cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*porewatera +fluxdissous)*dzpoi
                    flx_w2s(iv,i,j)=flx_w2s(iv,i,j)+fluxdissous
                  ENDDO
#endif
                  ! salinity and temperature
                  fluxdissous=(porewater-porewatera)*temp_bottom_MUSTANG(i,j)  
                  cv_sed(-1,ksmax,i,j)=(cv_sed(-1,ksmax,i,j)*porewatera +fluxdissous)*dzpoi
                  flx_w2s(-1,i,j)=flx_w2s(-1,i,j)+fluxdissous
                  fluxdissous=(porewater-porewatera)*sal_bottom_MUSTANG(i,j)  
                  cv_sed(0,ksmax,i,j)=(cv_sed(0,ksmax,i,j)*porewatera +fluxdissous)*dzpoi
                  flx_w2s(0,i,j)=flx_w2s(0,i,j)+fluxdissous                          
#endif
#if ! defined key_nofluxwat_IWS
                  ! water flux  
                  phieau_s2w(i,j)=phieau_s2w(i,j)+REAL((porewatera-porewater)*CELL_SURF(i,j),rlg)
#endif
                ELSE
                !output of pore-water to water column ; cvsed are unchanged
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
                  DO iv=nvp+1,nv_adv
                    flx_w2s(iv,i,j)=flx_w2s(iv,i,j)+(porewater-porewatera)*cv_sed(iv,k,i,j)  
                    ! will be divided by dt when taken into account in erosion, to be in quantity/m2/s
                    ! cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmax,i,j)*porewatera +fluxdissous)*dzpoi
                  ENDDO
#endif
                ! salinity and temperature
                  DO iv=-1,0
                    flx_w2s(iv,i,j)=flx_w2s(iv,i,j)+(porewater-porewatera)*cv_sed(iv,k,i,j)  
                  ENDDO
#endif
#if ! defined key_nofluxwat_IWS
                ! water flux 
                  phieau_s2w(i,j)=phieau_s2w(i,j)+REAL((porewatera-porewater)*CELL_SURF(i,j),rlg)
#endif
                ENDIF  !  porewater > porewatera
              ENDIF  ! dwsnew > dzsmin
 
              ! control writing :            
              !if(poro(k,i,j).ge.0.9999_rsh)then
              !  write(*,*)'dans effdep, a',t,'  en kij:',k,i,j
              !  write(*,555)k,dzs(k,i,j),c_sedtot(k,i,j),poro(k,i,j)
              !  write(*,*)dzsa,dzsgrv,dzssan,dzsmud,voldepgrv,voldepsan, &
              !            masdepmud,frdep,(cv_sed(iv,k,i,j),iv=1,4)
              !  write(*,*)igrav1,igrav2,isand1,isand2,imud1,imud2
              !  write(*,*)(fluevs(iv,i,j),iv=1,5)
              !  write(*,*)cvasr,cvascr,ddzs,ddzs1,ddzs2,ddzs3,ddzs4,ddzs5
              !  write(*,*)ddzsici,dvolgrv,dvolsan,cvolinigrv,cvolinisan,dzsaici
              !endif

           ! pour eviter l augmentation de l epaisseur de la couche de surface 
              IF(ksmax .LT. ksdmax .AND. ksmax > ksmi(i,j)) THEN
                 IF(dzs(ksmax,i,j) > dzsmax(i,j) + 5.0_rsh* dzsmin) THEN
                   dzs(ksmax+1,i,j)=MIN(dzs(ksmax,i,j)-dzsmax(i,j),dzsmax(i,j))
                   dzs(ksmax,i,j)=dzs(ksmax,i,j)-dzs(ksmax+1,i,j)
                   poro(ksmax+1,i,j)=poro(ksmax,i,j)
                   cv_sed(:,ksmax+1,i,j)=cv_sed(:,ksmax,i,j)
                   c_sedtot(ksmax+1,i,j)=c_sedtot(ksmax,i,j)
                   ksmax=ksmax+1
                 ENDIF
              ENDIF

             ENDIF  ! dzsnew > 0

            ELSE   !

             ! case 2. : no deposition
             ! ***********************  
               !      write(*,*)'no deposition ',i,j,fludep

               ! computing flx_w2s for dissolved var., to be used in the next time step:
               ! -----------------------------------------------------------------------
#if ! defined key_noTSdiss_insed
               DO iv=nvp+1,nv_adv
                 flx_w2s(iv,i,j)=0.0_rsh
               ENDDO
               flx_w2s(-1,i,j)=0.0_rsh
               flx_w2s(0,i,j)=0.0_rsh
#endif

               ! accounting for possible deposition of non constitutive particulate variable:
               ! ----------------------------------------------------------------------------
               IF(ksmax.GE.ksmi(i,j))THEN
                 DO iv=nvpc+1,nvp
                   cv_sed(iv,ksmax,i,j)=cv_sed(iv,ksmax,i,j)+flx_w2s_loc(iv)/MAX(dzsmin,dzs(ksmax,i,j))
                 ENDDO
               ELSE
                 ! no sediment
                 DO iv=nvpc+1,nvp

                  
                  !**TODO** Check if this is to code in CROCO ??
                  !cw_bottom_MUSTANG(iv,i,j)=cw_bottom_MUSTANG(iv,i,j)                                        &
                  !             +flx_w2s_loc(iv)/(epn_bottom_MUSTANG(i,j))
                  !iexchge_MPI_cvwat=1

                 ENDDO
               ENDIF


            ENDIF

              !!! : check  porosity (output)
#ifdef key_MUSTANG_specif_outputs
            DO k=ksmi(i,j),ksmax
              varspecif3Dk_save(1,k,i,j)=poro(k,i,j)  
            ENDDO
#endif

            ! updating ksma ( ksmax can be modified in the routine)
            ksma(i,j)=ksmax

#ifdef key_MUSTANG_specif_outputs
            varspecif3Dk_save(1,:,i,j)=0.0_rsh       ! poro_save
            DO k=ksmi(i,j),ksmax
              !!! BM: check  porosity (output)
              varspecif3Dk_save(1,k,i,j)=poro(k,i,j)  ! poro_save
            ENDDO
            varspecif2D_save(2,i,j)=dzs(ksmax,i,j)  ! dzs at sediment surface
#endif
       
          END IF ! test on htot
        END DO  ! loop on i
    END DO    ! loop on j

  END SUBROUTINE sed_MUSTANG_effdep
! fin effdep version V1
#endif

  !!==============================================================================
  SUBROUTINE sed_MUSTANG_fusion(i,j,ksmaxhere)
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_fusion  ***
   !&E
   !&E ** Purpose : makes fusion of 2  layers of the sediment (at the bottom)
   !&E              when the max. number of layers is reached
   !&E              and updating of the sediment column
   !&E
   !&E ** Description : 1DV 
   !&E  arguments IN : 
   !&E         i,j : cell 
   !&E         ksmaxhere : number of sediment layer in this cell i,j 
   !&E
   !&E  arguments OUT:
   !&E         ksmaxhere : new number of sediment layer after fusion 
   !&E              
   !&E ** Called by :  sed_effdep
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used

   !! * Arguments
   INTEGER,INTENT(IN)     ::  i,j
   INTEGER,INTENT(INOUT)  ::  ksmaxhere

   !! * Local declarations
   INTEGER        ::  k,iv,ksmin,kk
   REAL(KIND=rsh) ::  dzsa,poroa,cvolp,dzsmax_bottloc
#ifdef key_MUSTANG_V2 
   REAL(KIND=rsh) ::  sommud,cvolgrvsan,mass_tot,poro_mud_new,crel_mud_new,poro_kij
   REAL(KIND=rsh),DIMENSION(nvpc) :: frac_sed
   REAL(KIND=rsh),DIMENSION(1:nvp)  :: mass_sed
#endif            
#if ! defined key_noTSdiss_insed 
!#ifdef key_MUSTANG_V2 
!   REAL(KIND=rsh) ::  porowater_new,porowater1,porowater2
!#else
   REAL(KIND=rsh) ::  dzspwi
!#endif   
#endif   

   !!----------------------------------------------------------------------
   !! * Executable part

#ifdef key_MUSTANG_V2 
   ! fusion is done on ksmin+1 et ksmin+2 layers because ksmin layer cannot be
   ! eroded, and it enables to prevent the case in which there is a thick layer
   ! not erodable
   ksmin=ksmi(i,j)+1
#else 
   ksmin=ksmi(i,j) 
#endif
   kk=ksmin  
   dzsmax_bottloc=dzsmax_bottom
   DO WHILE (dzs(kk,i,j) .GE. dzsmax_bottloc  .AND. kk <ksmaxhere-nlayer_surf_sed)
      kk=kk+1
      IF(kk==ksmaxhere-nlayer_surf_sed) THEN
        kk=ksmin
        dzsmax_bottloc=dzsmax_bottloc+dzsmax_bottom
      ENDIF
   ENDDO
#ifdef key_MUSTANG_V2
   ! local evaluation of porosity as V1 in order to stay conservatif
   cvolp=0.0_rsh
   DO iv=1,nvp
       cvolp=cvolp+cv_sed(iv,kk,i,j)*typart(iv)/ros(iv)
   ENDDO
   poroa=1.0_rsh-cvolp
#else
   poroa=poro(kk,i,j)
#endif
   dzsa=dzs(kk,i,j)
   dzs(kk,i,j)=dzs(kk,i,j)+dzs(kk+1,i,j)
   cvolp=0.0_rsh
   c_sedtot(kk,i,j)=0.0_rsh
   DO iv=1,nvp
       cv_sed(iv,kk,i,j)=(cv_sed(iv,kk,i,j)*dzsa+cv_sed(iv,kk+1,i,j)*dzs(kk+1,i,j))/dzs(kk,i,j)
       c_sedtot(kk,i,j)=c_sedtot(kk,i,j)+cv_sed(iv,kk,i,j)*typart(iv)
       cvolp=cvolp+cv_sed(iv,kk,i,j)*typart(iv)/ros(iv)
   ENDDO
   poro(kk,i,j)=1.0_rsh-cvolp
   
#ifdef key_MUSTANG_V2    
   sommud=0.0_rsh
   DO iv=imud1,imud2
     sommud=sommud+cv_sed(iv,kk,i,j)
   ENDDO
   cvolgrvsan=0.0_rsh
   DO iv=igrav1,isand2
     cvolgrvsan=cvolgrvsan+cv_sed(iv,kk,i,j)/ros(iv)
   ENDDO
   crel_mud(kk,i,j)=sommud/(1.0_rsh-cvolgrvsan)
   crel_mud_new=crel_mud(kk,i,j)
   if(crel_mud(kk,i,j) ==0.0_rsh) then
     poro_mud_new=0.0_rsh
   else
     poro_mud_new=1.0_rsh-crel_mud(kk,i,j)/ros(1)
   endif
   
#endif
#if ! defined key_noTSdiss_insed 
     ! dissolved variable in pore waters and water fluxes at the interface
     ! --------------------------------------------------------------------
     !  ==> only cv_sed 
     !  ici on fait comme pour V1 pas d appel a MUSTANGV2_eval_dissvar_IWSflux
!#ifdef key_MUSTANG_V2
     ! fusion of 2 layers ==> code 2
     !porowater_new=poro(kk,i,j)*dzs(kk,i,j)
     !porowater1=dzsa*poroa
     !porowater2=dzs(kk+1,i,j)*poro(kk+1,i,j)
     !CALL MUSTANGV2_eval_dissvar_IWSflux(i,j,ksmax,2,  &
     !      porowater1=porowater1,            &
     !      cv_sed1=cv_sed(iv,kk,i,j),        &
     !      porowater2=porowater2,            &
     !      cv_sed2=cv_sed(iv,kk+1,i,j),      &
     !      porowater_new=porowater_new)
!#else
  ! MUSTANG_V1 et ici V2
   dzspwi=1.0_rsh/(dzs(kk,i,j)*poro(kk,i,j))
   DO iv=-1,0
      cv_sed(iv,kk,i,j)=dzspwi*(cv_sed(iv,kk,i,j)*dzsa*poroa &
                           +cv_sed(iv,kk+1,i,j)*dzs(kk+1,i,j)*poro(kk+1,i,j))
   ENDDO
#if ! defined key_Pconstitonly_insed
   DO iv=nvp+1,nv_adv
      cv_sed(iv,kk,i,j)=dzspwi*(cv_sed(iv,kk,i,j)*dzsa*poroa &
                                +cv_sed(iv,kk+1,i,j)*dzs(kk+1,i,j)*poro(kk+1,i,j))
   ENDDO
#endif
#endif

#ifdef key_MUSTANG_V2
    ! re-evaluation of poro as V2
     mass_tot = 0.0_rsh
     DO iv=igrav1,imud2
       mass_sed(iv)=cv_sed(iv,kk,i,j)*dzs(kk,i,j)
       mass_tot=mass_tot+mass_sed(iv)
     END DO

     DO iv=igrav1,imud2
       frac_sed(iv)=mass_sed(iv)/mass_tot
     END DO
     CALL MUSTANGV2_comp_poro_mixsed(frac_sed, poro_mud_new,   &
                                crel_mud_new, poro_kij) 
                             !   pi, crel_mud_new, poro_kij) 
     poro(kk,i,j)=poro_kij
     poro_mud(kk,i,j)=poro_mud_new
#endif

   DO k=kk+1,ksmaxhere-1
       dzs(k,i,j)=dzs(k+1,i,j)
       poro(k,i,j)=poro(k+1,i,j)
#ifdef key_MUSTANG_V2    
       poro_mud(k,i,j)=poro_mud(k+1,i,j)
       crel_mud(k,i,j)=crel_mud(k+1,i,j)
#endif
       c_sedtot(k,i,j)=c_sedtot(k+1,i,j)
       DO iv=-1,nv_adv
           cv_sed(iv,k,i,j)=cv_sed(iv,k+1,i,j)
       ENDDO
   ENDDO
   ksmaxhere=ksmaxhere-1
  

  END SUBROUTINE sed_MUSTANG_fusion
      

 
!!==============================================================================
#if defined key_MUSTANG_slipdeposit
  SUBROUTINE sed_MUSTANG_slipdepo(ifirst, ilast, jfirst, jlast)
   
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_slipdepo  ***
   !&E
   !&E ** Purpose : computes slipping settling if slope bathy
   !&E
   !&E ** Description :  
   !&E  arguments IN : 
   !&E          loops : ifirst,ilast,jfirst,jlast
   !&E          
   !&E   variables in  (comMUSTANG)      flx_w2s_sum 
   !&E          
   !&E   variables OUT: (comMUSTANG)
   !&E         flx_w2s_corip1,flx_w2s_corim1,flx_w2s_corjp1,flx_w2s_corjm1
   !&E         flx_w2s_corin
   !&E          
   !&E ** Called by :  MUSTANG_deposition
   !&E
   !&E--------------------------------------------------------------------------
   !! * Arguments
   INTEGER, INTENT(IN)    :: ifirst, ilast, jfirst, jlast

   !! * Local declarations
   INTEGER                :: iv, i, j
   REAL(KIND=rsh)         :: cordepfluw, cordepflue, cordepflus, cordepflun, cordepflu

   !! * Executable part 
   DO j = jfirst, jlast 
     DO i = ifirst, ilast
        IF(htot(i,j) > h0fond) THEN

         cordepfluw = max(0.0_rsh, slopefac * SLOPE_W)
         cordepflue = max(0.0_rsh, slopefac * SLOPE_E)
         cordepflus = max(0.0_rsh, slopefac * SLOPE_S)
         cordepflun = max(0.0_rsh, slopefac * SLOPE_N)
         cordepflu = cordepfluw + cordepflue + cordepflus + cordepflun
         IF (cordepflu.gt.1.0_rsh)THEN
           cordepfluw = cordepfluw / cordepflu
           cordepflue = cordepflue / cordepflu
           cordepflus = cordepflus / cordepflu
           cordepflun = cordepflun / cordepflu
           cordepflu = 1.0_rsh
         ENDIF
         DO iv = isand2+1, nvp
           flx_w2s_corim1(iv,i,j) = cordepfluw * flx_w2s_sum(iv,i,j) * CELL_SURF(i,j) / SURF_NEAR_W
           flx_w2s_corip1(iv,i,j) = cordepflue * flx_w2s_sum(iv,i,j) * CELL_SURF(i,j) / SURF_NEAR_E
           flx_w2s_corjm1(iv,i,j) = cordepflus * flx_w2s_sum(iv,i,j) * CELL_SURF(i,j) / SURF_NEAR_S
           flx_w2s_corjp1(iv,i,j) = cordepflun * flx_w2s_sum(iv,i,j) * CELL_SURF(i,j) / SURF_NEAR_N
           flx_w2s_corin(iv,i,j) = -cordepflu * flx_w2s_sum(iv,i,j)
         ENDDO
        ENDIF

      ENDDO
   ENDDO

  END SUBROUTINE sed_MUSTANG_slipdepo
#endif
    !!==============================================================================
#if ! defined key_noTSdiss_insed
   SUBROUTINE sed_MUSTANG_Temperatur_in_sed(ifirst, ilast, jfirst, jlast, dt_true, dtinv)
! 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_Temperatur_in_sed ***
   !&E
   !&E ** Purpose : dynamic in sediment : processes of Tempertur diffusion in sediment 
   !&E
   !&E ** Description :
   !&E        arguments IN :
   !&E            dt_true : time step
   !&E            parameters : 
   !&E
   !&E        variales OUT :
   !&E            fludiff : substance flux de temperature at the  interface water/sediment due to diffusion
   !&E
   !&E ** Called by :  MUSTANG_update
   !&E
   !&E ** External calls : 
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
    USE cometemp,  ONLY : chp


   !! * Arguments
   INTEGER, INTENT(IN)            :: ifirst,ilast,jfirst,jlast                           
   REAL(KIND=rsh),INTENT(IN)      :: dtinv
   REAL(KIND=rlg),INTENT(IN)      :: dt_true

   
   !! * Local declarations
   INTEGER                                   :: i,j,k,iv,ksmin,ksmax
   REAL(KIND=rsh)                            :: somcmud,phi_surfsed,phi_bottsed
   REAL(KIND=rsh)                            :: mu0_tempsed,frmud,hsedloc
   !REAL(KIND=rsh),DIMENSION(ksdmin:ksdmax,-1:nv_adv)      :: difbio                   
   REAL(KIND=rsh),DIMENSION(ksdmin  :ksdmax)              :: aa,bb,cc,dd
   REAL(KIND=rsh)                                         :: hcrit,dzssdt
   REAL(KIND=rsh),DIMENSION(ksdmin :ksdmax)               :: disvi,dzsiin




   !!----------------------------------------------------------------------
   !! * Executable part

     DO j=jfirst,jlast
      DO i=ifirst,ilast

       IF (ksma(i,j) .GE. ksmi(i,j) .AND. ksma(i,j) > 0)THEN
           
         ! memorization of the porosities, the volume of interstitial water
         ! the masses of dissolved substances
         ! and the thickness of the surface layer
         ! ----------------------------------------------------------
         DO k=ksmi(i,j),ksma(i,j)-1
            dzsiin(k)=2.0_rsh/(dzs(k,i,j)+dzs(k+1,i,j))
         ENDDO
              
         ksmax=ksma(i,j)
         ksmin=ksmi(i,j)
        
         !  IF(l_biodiffs ) THEN
         !     ! calculation of biodiffusion coef for temperature
         !     CALL sed_MUSTANG_coefbioturb(i,j,difbio) !  To review to differentiate particulate and dissolved
         !  ENDIF

         hcrit=MAX(0.0_rsh,htot(i,j)-h0fond)
         disvi(:)=0.0_rsh
         phi_surfsed=0.0_rsh
         phi_bottsed=0.0_rsh

         ! calculation of diffusion / biodiffusion rates in the sediment
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!FG(28/06/2018) mu_tempsed = f(frmud) dans la couche de surface
         k=ksmax
         !frmud=0.0_rsh
         !somcmud=0.0_rsh
         !DO iv=imud1,imud2
         !    somcmud=somcmud+cv_sed(iv,k,i,j)
         !ENDDO
         !frmud=somcmud/c_sedtot(k,i,j)
         !mu0_tempsed=mu_tempsed1*frmud**2+mu_tempsed2*frmud+mu_tempsed3

            ! utiliser diffusivite thermique dans l eau et non dans le sediment
            ! eta : conductivite thermique de l eau = 0.8
            ! chp = capacite calorifique de l eau 
         disvi(ksmax)=0.8_rsh/(min(epdifi,htot(i,j))+epsilon_MUSTANG)/(chp*roswat_bot(i,j))
         disvi(ksmax)=disvi(ksmax)*hcrit/(hcrit+epsilon_MUSTANG)
         hsedloc=0.0_rsh
         DO k=ksmin,ksmax-1
               hsedloc=hsedloc+dzs(k,i,j)
               ! calcul de la diffusivite thermique dans le sediment
               !!!FG(28/06/2018) mu_tempsed = f(frmud) dans chaque couche
               frmud=0.0_rsh
               somcmud=0.0_rsh
               DO iv=imud1,imud2
                 somcmud=somcmud+cv_sed(iv,k,i,j)
               ENDDO
               frmud=somcmud/c_sedtot(k,i,j)
               mu0_tempsed=mu_tempsed1*frmud**2+mu_tempsed2*frmud+mu_tempsed3
               !mu0_tempsed=mu_tempsed3
               !!!FG(29/06/2018)mu_tempsedk=mu0_tempsed
               ! formulation en fonction de la tortuosite Dsed=Dpure/Tortuosite^2
               ! tortuosite fontion de la porosite (eq 4.120)
               ! in Boudreau 1997 p 132 ! totuosity^2=1-log(poro^2)
               !!!FG(29/06/2018)disvi(k,-1)=(mu_tempsedk/(1.0_rsh-2.0_rsh*LOG(poro(k,i,j)))+difbio(k,-1))*dzsiin(k)
               !disvi(k)=(mu0_tempsed/(1.0_rsh-2.0_rsh*LOG(poro(k,i,j))))*dzsiin(k)
               disvi(k)=mu0_tempsed*dzsiin(k)
         ENDDO


         hsedloc=hsedloc+dzs(ksmax,i,j)
         phi_surfsed=phitemp_s(i,j)*dtinv
         IF(hsedloc > 0.0_rsh)phi_bottsed=MAX(0.0_rsh,MIN(phi_surfsed,phi_surfsed*(epsedmax_tempsed-hsedloc)*(epsedmin_tempsed/hsedloc)))
         ! end of cumul : reset phitemp_s (si pas le meme pas de temps , mais ici on a le meme dt_true)
         phitemp_s(i,j)=0.0_rsh
       
         iv=-1
         IF(ksmax > ksmin)THEN
!
                !cas ou il y a plusieurs couches:
                !--------------------------------
                !coefficients de la matrice
  
           k=ksmin
           dzssdt=dzs(k,i,j)*dtinv
           aa(k)=0.0_rsh
           bb(k)=dzssdt+disvi(k) 
           cc(k)=-disvi(k) 
           dd(k)=cv_sed(iv,k,i,j)*dzs(k,i,j)*dtinv-phi_bottsed

           DO k=ksmin+1,ksmax-1
                 dzssdt=dzs(k,i,j)*dtinv
                 aa(k)=-disvi(k-1) 
                 bb(k)=dzssdt+disvi(k-1)+disvi(k) 
                 cc(k)=-disvi(k) 
                 dd(k)=cv_sed(iv,k,i,j)*dzs(k,i,j)*dtinv
           ENDDO

           k=ksmax
           aa(k)=-disvi(k-1) 
           bb(k)=dzs(k,i,j)*dtinv+disvi(k-1)+disvi(k)  
           cc(k)=0.0_rsh
           fludif(iv,i,j)=disvi(k)*temp_bottom_MUSTANG(i,j)
           dd(k)=cv_sed(iv,k,i,j)*dzs(k,i,j)*dtinv+fludif(iv,i,j)+phi_surfsed
!
           !resolution du systeme tridiagonal
           dd(1)=dd(1)/bb(1)
           DO k=ksmin+1,ksmax
                  bb(k)=bb(k)-aa(k)*cc(k-1)/bb(k-1)
                  dd(k)=(dd(k)-aa(k)*dd(k-1))/bb(k)
           ENDDO
           cv_sed(iv,ksmax,i,j)=dd(ksmax)
           DO k=ksmax-1,ksmin,-1
                  cv_sed(iv,k,i,j)=dd(k)-cc(k)*cv_sed(iv,k+1,i,j)/bb(k)
           ENDDO

          !  flux d echange a l interface eau/sediment: TEMPERATURE
          !  -------------------------------------------------------
          !  fludif est le flux par diffusion a l interface, >0 vers le bas  ==> M/m2
           fludif(iv,i,j)=  &
                    (fludif(iv,i,j)-disvi(ksmax)*cv_sed(iv,ksmax,i,j))*REAL(dt_true,rsh)

!
         ELSE
!
               !cas ou il n y a qu une couche (ksmax=ksmin): TEMPERATURE
               !--------------------------------------------       

           k=ksmax  
           IF(dzs(k,i,j) > dzsmin) THEN
             fludif(iv,i,j)=disvi(k)*temp_bottom_MUSTANG(i,j)   
             cv_sed(iv,k,i,j)=  &
                    (cv_sed(iv,k,i,j)*dzs(k,i,j)+fludif(iv,i,j)*REAL(dt_true,rsh)    &
                     +phi_surfsed-phi_bottsed)           &
                   /(dzs(k,i,j)+(disvi(k))*REAL(dt_true,rsh))
             fludif(iv,i,j)=  &
                    (fludif(iv,i,j)-disvi(ksmax)*cv_sed(iv,ksmax,i,j))*REAL(dt_true,rsh)

           ELSE
             fludif(iv,i,j)=0.0_rsh
             cv_sed(iv,k,i,j)=temp_bottom_MUSTANG(i,j)
           ENDIF

         ENDIF
!
        ELSE
          phitemp_s(i,j)=0.0_rsh
        ENDIF    ! test if existence sediment
    
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!   END OF DYNAMIC PROCESS IN SEDIMENT                !!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ENDDO
    ENDDO
   
  END SUBROUTINE sed_MUSTANG_Temperatur_in_sed
#endif
    !!==============================================================================

   SUBROUTINE sed_MUSTANG_consol_diff_bioturb(ifirst, ilast, jfirst, jlast,  &
                saliref_lin, temperef_lin, dt_true)
! 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_consol_diff_bioturb ***
   !&E
   !&E ** Purpose : dynamic in sediment : processes of consolidation 
   !&E                   and/or diffusion and/or bioturbation
   !&E
   !&E ** Description :
   !&E        arguments IN :
   !&E            dt_true : time step
   !&E            parameters : saliref_lin,temperef_lin
   !&E
   !&E        variales OUT :
   !&E            phieau_s2w : water flux at the interface water/sediment due to consolidation
   !&E                        (in m3 water exchanged during consolidation time step)
   !&E            fluconsol : substance flux at the  interface water/sediment due to consolidation
   !&E            fludiff : substance flux de substance at the  interface water/sediment due to diffusion
   !&E            +  sediment changed  (cv_sed, poro, dzs, ksma.. modified)
   !&E
   !&E ** Called by :  MUSTANG_update
   !&E
   !&E ** External calls : 
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used


   !! * Arguments
   INTEGER, INTENT(IN)            :: ifirst, ilast, jfirst, jlast                           
   REAL(KIND=rsh),INTENT(IN)      :: saliref_lin, temperef_lin
   REAL(KIND=rlg),INTENT(IN)      :: dt_true
   
   !! * Local declarations
   INTEGER                                   :: i,j,k,iv,ksmin,ksmax,celldry,ksmaxa,kp,kl,nv_use2,ivv
   REAL(KIND=rsh)                            :: volpwa_tot,volpwnv_tot,volpwk,roro,dt_sed_inv,dtsdzs,   &
                                                sigmadjge,somdeltaro,halfmaslayer,somcsan,somcgrav,   &
                                                somcmud,cvolp,rapcsed,csednew,poroam1,poroa,     &
                                                dzsa,dzsam1,dzsisu,rowinv,dtsdzsmin,  &
                                                volak,volakm1,vola,vola_inv,ksmax_turb_part,ksmax_turb_diss
   !REAL(KIND=rsh)                            :: dvol,dvol_fusion 
   REAL(KIND=rsh),DIMENSION(ksdmin  :ksdmax) :: volpwak,stateconsol,permeab,sigmapsg,  &
                                                poroin,dzsiin,loadograv,sed_rate,hinder,cseda,poroanc
   REAL(KIND=rsh),DIMENSION(ksdmin-1:ksdmax) :: winter
   REAL(KIND=rlg)                            :: dt_sed_cor,dt_sed_eff,dtrest,dtiter                                           
   REAL(KIND=rsh),DIMENSION(ksdmin:ksdmax,-1:nv_adv)      :: difbio                   
   REAL(KIND=rsh),DIMENSION(nvp+1,ksdmin:ksdmax)           :: cvsednew
#if ! defined key_noTSdiss_insed
 !  REAL(KIND=rsh),DIMENSION(-1:nv_adv)                    :: flu_fusion              
   REAL(KIND=rsh),DIMENSION(ksdmin  :ksdmax)              :: aa,bb,cc,dd
   REAL(KIND=rsh)                                         :: hcrit,dzssdt,xdifs1b
   REAL(KIND=rsh)                                         :: Sc,p,mu,nu,D0
   REAL(KIND=rsh),DIMENSION(ksdmin :ksdmax,0:nv_adv-nvp) :: disvi
   REAL(KIND=rsh),DIMENSION(0:nv_adv-nvp)                :: conc_bottom
#endif
#ifdef key_MUSTANG_V2
   CHARACTER(len=19)                :: tool_sectodat
   REAL(KIND=rsh)                   :: cv_sed_tot,mass_tot,sommud,cvolgrvsan
   REAL(KIND=rsh),DIMENSION(1:nvpc) :: frac_sed
   REAL(KIND=rsh),DIMENSION(1:nvp)  :: mass_sed
#endif



   !!----------------------------------------------------------------------
   !! * Executable part
#if ! defined key_Pconstitonly_insed
   nv_use2=nv_adv
#else
   nv_use2=0
#endif

   IF(CURRENT_TIME .GE. tstart_dyninsed .AND. CURRENT_TIME .GE. t_dyninsed) THEN
   
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!! INTERNAL DYNAMICS INTO THE SEDIMENT   !!!!!!
      !!!!  CONSOLIDATION/BIOTURBATION/DIFFUSION  !!!!!!
      !!!!  UPDATE ALL dt_dyninsed                !!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     dt_sed_cor=MAX(dt_dyninsed,dt_true)
     dt_sed_eff=dt_sed_cor+(CURRENT_TIME-t_dyninsed)
     dt_sed_inv=1.0_rsh/REAL(dt_sed_eff,rsh)
     phieau_s2w_consol(:,:)=0.0_rsh

     DO j=jfirst,jlast
      DO i=ifirst,ilast
       difbio(ksdmin:ksdmax,-1:nv_adv)=0.0_rsh
       poroin(ksdmin:ksdmax)=0.0_rsh
       volpwak(:)=0.0_rsh
       IF (ksma(i,j) .GE. ksmi(i,j) .AND. ksma(i,j) > 0)THEN
         celldry=0
         IF(htot(i,j) .LE. h0fond)celldry=1
           
         ! memorization of the porosities, the volume of interstitial water
         ! the masses of dissolved substances
         ! and the thickness of the surface layer
         ! ----------------------------------------------------------
         volpwak(ksma(i,j))=poro(ksma(i,j),i,j)*dzs(ksma(i,j),i,j)
         volpwa_tot=volpwak(ksma(i,j))
         poroin(ksma(i,j))=poro(ksma(i,j),i,j)
         DO k=ksmi(i,j),ksma(i,j)-1
            dzsiin(k)=2.0_rsh/(dzs(k,i,j)+dzs(k+1,i,j))
            poroin(k)=0.5_rsh*(poro(k+1,i,j)+poro(k,i,j))
            volpwak(k)=poro(k,i,j)*dzs(k,i,j)
            volpwa_tot=volpwa_tot+volpwak(k)
         ENDDO

         winter(:)=0.0_rsh
                      
         ksmax=ksma(i,j)
         ksmin=ksmi(i,j)

         !!! calculating the density of interstitial water
          ! it is assumed to be equal to the density of the water at the bottom of the water column
          ! if there is no water, the salinity and temperature are not known
          ! if defined key_noTSdiss_insed, temp and sal are not calculated in the sediment
          ! and are false (= valmanq) in the new layers
         IF(celldry ==1) THEN
         ! no water
         !  density estimated from the bottom sediment layer or surface sediment layer
#ifdef key_noTSdiss_insed
           !roro=(RHOREF*(1.0_rsh+0.0008_rsh*(cv_sed(0,1,i,j)-saliref_lin)  &
           !                            -0.00016_rsh*(cv_sed(-1,1,i,j)-temperef_lin)))
           roro=RHOREF
#else
           roro=(RHOREF*(1.0_rsh+0.0008_rsh*(cv_sed(0,ksmax,i,j)-saliref_lin)  &
                                      -0.00016_rsh*(cv_sed(-1,ksmax,i,j)-temperef_lin)))
#endif
         ELSE
               roro=roswat_bot(i,j)
         ENDIF
         rowinv=1.0_rsh/roro

         sed_rate(ksdmin:ksdmax)=0.0_rsh
         hinder(ksdmin:ksdmax)=0.0_rsh
         !dvol=0.0_rsh
         !dvol_fusion=0.0_rsh
         volak=0.0_rsh
         volakm1=0.0_rsh
!#if ! defined key_noTSdiss_insed
!         flu_fusion(:)=0.0_rsh
!#endif
         
         dtrest=dt_sed_eff
         
         ! test : to avoid to much computation, we do a single call for estimate particulate and dissolved mixing coef 
         ! even if there is consolidation processes with varying concentrations 
         IF(l_bioturb) THEN
           ksmax_turb_part=ksmax 
          ! if(j==5 .OR. j==20)write(*,*)'avt 1',t,j,difbio(ksma(i,j)-5:ksma(i,j),-1)
           CALL sed_MUSTANG_coefbioturb_part(i,j,difbio) !  To review to differentiate particulate mixing coef
          ! if(j==5 .OR. j==20)write(*,*)'aprs 1',t,j,difbio(ksma(i,j)-5:ksma(i,j),-1)
         ENDIF
#if ! defined key_noTSdiss_insed
         IF(l_biodiffs) THEN
           ksmax_turb_diss=ksmax 
           CALL sed_MUSTANG_coefbioturb_diss(i,j,difbio) !  To review to differentiate dissolved mixing coef
         ENDIF
#endif

  100    CONTINUE
                   
         cseda(:)=0.0_rsh     
         DO k=ksmin,ksmax
              cseda(k)=c_sedtot(k,i,j)
         ENDDO
            
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!  CONSOLIDATION OR BIOTURBATION  !!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         dtiter=MIN(dtrest,subdt_consol)

         IF(l_consolid .OR. l_bioturb) THEN
             
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!  CONSOLIDATION   !!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          IF(l_consolid .AND. ksma(i,j) > ksmi(i,j)) THEN
          
          ! computation of permeability, effective stress and hydraulic load
          ! ----------------------------------------------------------------
            CALL sed_MUSTANG_constitutivrel(i,j,stateconsol,permeab,sigmapsg)        

            sigmadjge=0.0_rsh
            DO k=ksmax,ksmin,-1
              somdeltaro=0.0_rsh
              DO iv=1,nvpc
              ! rosmrowsros=(ros-row)/ros
                somdeltaro=somdeltaro+cv_sed(iv,k,i,j)*rosmrowsros(iv)
              ENDDO
              halfmaslayer=.5*dzs(k,i,j)*somdeltaro
              sigmadjge=sigmadjge+halfmaslayer
              ! loadograv : exces de pression d eau interstitielle au milieu de la couche
              ! sigmadjge : sigma dejauge (sans la part de l eau)
              ! sigmapsg : contrainte effective (transmise de grains a grains)
               ! loadograv: excess of interstitial water pressure in the middle of the layer
               ! sigmadjge: sigma unseparated (without the share of water)
               ! sigmapsg: effective stress (transmitted from grain to grain)
              loadograv(k)=MAX(0.0_rsh,sigmadjge-sigmapsg(k))
#if defined key_MUSTANG_specif_outputs && defined key_MUSTANG_add_consol_outputs
              varspecif3D_save(3,k,i,j)=loadograv(k)
              varspecif3D_save(4,k,i,j)=permeab(k)
              varspecif3D_save(5,k,i,j)=sigmapsg(k)
              varspecif3D_save(9,k,i,j)=sigmadjge
              varspecif3D_save(10,k,i,j)=stateconsol(k)
#endif

              ! total weight of all layers above layer k
              sigmadjge=sigmadjge+halfmaslayer
            ENDDO
          
            ! computing new concentrations from mass conservation
            !  only particulate variables are modified during consolidation
            ! ---------------------------------------------------       
            k=ksmin
            !sed_rate : advection speed of mud particles
            !  Darcy's law
            ! stateconsol : state of consolidation indicator
            sed_rate(k)=-rowinv*(permeab(k+1)+permeab(k))                                 &
                *MIN(0.0_rsh,loadograv(k+1)-loadograv(k))*stateconsol(k)         &
                /(dzs(k+1,i,j)+dzs(k,i,j))
            somcgrav=0.0_rsh
            DO iv=igrav1,igrav2
              somcgrav=somcgrav+cv_sed(iv,k,i,j)
            ENDDO
            somcsan=0.0_rsh
            DO iv=isand1,isand2
               somcsan=somcsan+cv_sed(iv,k,i,j)
            ENDDO
            somcmud=0.0_rsh
            DO iv=imud1,imud2
             somcmud=somcmud+cv_sed(iv,k,i,j)
            ENDDO
            !hinder : entravement du sable/gravier (sans dimension) entre 0 et 1
            ! ATTENTION : hypothese de ros homogene quelque soit le type de sediment
            ! entrave par le trop de vase qui empeche la segregation sable/vase (1er terme)
            ! entrave par le trop de sablequi ne permet pas le passage entre les grains (2ieme terme)
            ! hinder: shackling sand / gravel (dimensionless) between 0 and 1
            ! WARNING: homogeneous ros hypothesis whatever the type of sediment
            ! obstruction by the too much vase which prevents the segregation sand / vase (1st term)
            ! obstructed by too much sand that does not allow the passage between the grains (2nd term)
            hinder(k)=(MAX(0.0_rsh,MIN(1.0_rsh-somcmud/(1-(somcsan+somcgrav)/ros(1))/csegreg,   &
                   1.0_rsh-(somcsan+somcgrav)/csandseg)))**4.65_rsh
#if defined key_MUSTANG_specif_outputs && defined key_MUSTANG_add_consol_outputs
            varspecif3D_save(7,k,i,j)=hinder(k)
            varspecif3D_save(8,k,i,j)=sed_rate(k)
#endif
                   
            DO k=ksmin+1,ksmax-1

              sed_rate(k)=-rowinv*(permeab(k+1)+permeab(k))                               &
                  *MIN(0.0_rsh,loadograv(k+1)-loadograv(k))*stateconsol(k)       &
                  /(dzs(k+1,i,j)+dzs(k,i,j))
              somcgrav=0.0_rsh
              DO iv=igrav1,igrav2
                somcgrav=somcgrav+cv_sed(iv,k,i,j)
              ENDDO
              somcsan=0.0_rsh
              DO iv=isand1,isand2
                somcsan=somcsan+cv_sed(iv,k,i,j)
              ENDDO
              somcmud=0.0_rsh
              DO iv=imud1,imud2
                somcmud=somcmud+cv_sed(iv,k,i,j)
              ENDDO
              hinder(k)=(MAX(0.0_rsh,MIN(1.0_rsh-somcmud/(1-(somcsan+somcgrav)/ros(1))/csegreg,   &
                   1.0_rsh-(somcsan+somcgrav)/csandseg)))**4.65_rsh
            ENDDO
        
            k=ksmax
            dtsdzs=REAL(dtiter,rsh)/dzs(k,i,j)
#if defined key_MUSTANG_specif_outputs && defined key_MUSTANG_add_consol_outputs
            varspecif3D_save(6,k,i,j)=dtsdzs
#endif
           ! implicit scheme decentred upstream -cv_sed(k+1) known just before, so implicit
            DO iv=igrav1,isand2
              cv_sed(iv,ksmax,i,j)=cv_sed(iv,ksmax,i,j)                               &
                                 /(1.0_rsh+dtsdzs*MAX(sed_rate(ksmax-1),hinder(ksmax-1)*ws_sand(iv)))
            ENDDO
            DO iv=imud1,nvp
              cv_sed(iv,ksmax,i,j)=cv_sed(iv,ksmax,i,j)/(1.0_rsh+dtsdzs*sed_rate(ksmax-1))
            ENDDO
#if ! defined key_noTSdiss_insed
              cv_sed(-1,ksmax,i,j)=cv_sed(-1,ksmax,i,j)/(1.0_rsh+dtsdzs*sed_rate(ksmax-1))
#endif
            DO k=ksmax-1,ksmin+1,-1
              dtsdzs=REAL(dtiter,rsh)/dzs(k,i,j)
#if defined key_MUSTANG_specif_outputs && defined key_MUSTANG_add_consol_outputs
              varspecif3D_save(6,k,i,j)=dtsdzs
#endif
              DO iv=igrav1,isand2
                cv_sed(iv,k,i,j)=(cv_sed(iv,k,i,j)+dtsdzs*MAX(sed_rate(k),ws_sand(iv)*               &
                               hinder(k))*cv_sed(iv,k+1,i,j))                        &
                               /(1.0_rsh+dtsdzs*MAX(sed_rate(k-1),hinder(k-1)*ws_sand(iv)))
              ENDDO
              DO iv=imud1,nvp
                cv_sed(iv,k,i,j)=(cv_sed(iv,k,i,j)+dtsdzs*sed_rate(k)*cv_sed(iv,k+1,i,j))/(1.0_rsh+dtsdzs*sed_rate(k-1))
              ENDDO
#if ! defined key_noTSdiss_insed
              cv_sed(-1,k,i,j)=(cv_sed(-1,k,i,j)+dtsdzs*sed_rate(k)*cv_sed(-1,k+1,i,j))/(1.0_rsh+dtsdzs*sed_rate(k-1))
#endif
            ENDDO
            k=ksmin
            dtsdzs=REAL(dtiter,rsh)/dzs(k,i,j)
#if defined key_MUSTANG_specif_outputs && defined key_MUSTANG_add_consol_outputs
            varspecif3D_save(6,k,i,j)=dtsdzs
#endif
            DO iv=igrav1,isand2
              cv_sed(iv,k,i,j)=(cv_sed(iv,k,i,j)+dtsdzs*MAX(sed_rate(k),ws_sand(iv)*hinder(k))*cv_sed(iv,k+1,i,j))
            ENDDO
            DO iv=imud1,nvp
              cv_sed(iv,k,i,j)=(cv_sed(iv,k,i,j)+dtsdzs*sed_rate(k)*cv_sed(iv,k+1,i,j))
            ENDDO
#if ! defined key_noTSdiss_insed
              cv_sed(-1,k,i,j)=(cv_sed(-1,k,i,j)+dtsdzs*sed_rate(k)*cv_sed(-1,k+1,i,j))
#endif

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!  UPDATE OF PARTICULATE CONCENTRATIONS           !!!!!!!!!!!
            !!!!!!!!!!!       AFTER CONSOLIDATION                       !!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! computing new concentrations/porosities :
            DO k=ksmin,ksmax
              c_sedtot(k,i,j)=0.0_rsh
              cvolp=0.0_rsh
              DO iv=1,nvpc
               c_sedtot(k,i,j)=c_sedtot(k,i,j)+cv_sed(iv,k,i,j)
               cvolp=cvolp+cv_sed(iv,k,i,j)/ros(iv)
              ENDDO
              poro(k,i,j)=1.0_rsh-cvolp

#ifdef key_MUSTANG_V2
!!! BM :  on verifiera a la fin du dt consol les caracteristiques de notre couche
!        --> pourra aboutir a un calcul de nouvelle porosite et donc a une actualisation de lepaisseur
!           de la couche en cas d un resultat juge  non coherent 
              IF(poro(k,i,j).LT.poro_min)THEN
                poro(k,i,j)=poro_min
                !!! CAUTION: we assume that ros is identical regardless of the particles
                csednew=ros(1)*(1.0_rsh-poro_min)
#else           
               ! if(porosity too small (<32%), it is increased until (1.0_rsh-cvolmaxmel):
               !  dissolved variables concentrations are inchanged
               !------------------------------------------------------------------------
              IF(poro(k,i,j).LT.1.0_rsh-cvolmaxmel-0.01_rsh)THEN
                poro(k,i,j)=1.0_rsh-cvolmaxmel
                !!! CAUTION: we assume that ros is identical regardless of the particles
                csednew=ros(1)*cvolmaxmel
#endif
            
                rapcsed=csednew/c_sedtot(k,i,j)
                DO iv=1,nvp
                  cv_sed(iv,k,i,j)=cv_sed(iv,k,i,j)*rapcsed
                ENDDO
#if ! defined key_noTSdiss_insed
                cv_sed(-1,k,i,j)=cv_sed(-1,k,i,j)*rapcsed
#endif
                c_sedtot(k,i,j)=csednew
                dzs(k,i,j)=dzs(k,i,j)/rapcsed
              ENDIF 
            
#ifdef key_MUSTANG_V2
              sommud=0.0_rsh
              DO iv=imud1,imud2
                sommud=sommud+cv_sed(iv,k,i,j)
              ENDDO
              cvolgrvsan=0.0_rsh
              DO iv=igrav1,isand2
                cvolgrvsan=cvolgrvsan+cv_sed(iv,k,i,j)/ros(iv)
              ENDDO
              crel_mud(k,i,j)=sommud/(1.0_rsh-cvolgrvsan)
              poro_mud(k,i,j)=1.0_rsh-sommud/ros(1)
#endif
            
            ENDDO
        
            ! again calculation for ksmax because we vary dzs in the surface layer
            ! we lost material but c_sedtot remains constant => dzs change
            k=ksmax
#if ! defined key_MUSTANG_V2   
            IF(c_sedtot(k,i,j) .GT. 1.e-20) THEN
#endif
               dzsa=dzs(k,i,j)
               dzs(k,i,j)=c_sedtot(k,i,j)*dzsa/cseda(k)
               cvolp=0.0_rsh
               c_sedtot(k,i,j)=0.0_rsh
               DO iv=1,nvp
                 cv_sed(iv,k,i,j)=cv_sed(iv,k,i,j)*dzsa/dzs(k,i,j)
                 c_sedtot(k,i,j)=c_sedtot(k,i,j)+typart(iv)*cv_sed(iv,k,i,j)
                 cvolp=cvolp+typart(iv)*cv_sed(iv,k,i,j)/ros(iv)
               ENDDO
#if ! defined key_noTSdiss_insed
               cv_sed(-1,k,i,j)=cv_sed(-1,k,i,j)*dzsa/dzs(k,i,j)
#endif
               ! in the case where ros is homogeneous whatever the type of sediment,
               ! poro should remain unchanged too
                poro(k,i,j)=1.0_rsh-cvolp
#if ! defined key_MUSTANG_V2   
             ELSE
                ksmax=MAX(1,ksmax-1)
             ENDIF 
#endif
#ifdef key_MUSTANG_V2
             sommud=0.0_rsh
             DO iv=imud1,imud2
               sommud=sommud+cv_sed(iv,k,i,j)
             ENDDO
             cvolgrvsan=0.0_rsh
             DO iv=igrav1,isand2
              cvolgrvsan=cvolgrvsan+cv_sed(iv,k,i,j)/ros(iv)
             ENDDO
             crel_mud(k,i,j)=sommud/(1.0_rsh-cvolgrvsan)
             poro_mud(k,i,j)=1.0_rsh-sommud/ros(1)
#endif
          
             ! put in comment because useless
             ! if(porosity too small (<32%), it is increased until (1.0_rsh-cvolmaxmel):
             !IF(poro(k,i,j) < 1.0_rsh-cvolmaxmel-0.01_rsh)THEN
             !   poro(k,i,j)=1.0_rsh-cvolmaxmel
             !   !!! ATTENTION : on suppose ros identique quelque soit les particules
             !   csednew=ros(1)*cvolmaxmel
             !   rapcsed=csednew/c_sedtot(k,i,j)
             !   DO iv=1,nvp
             !     cv_sed(iv,k,i,j)=cv_sed(iv,k,i,j)*rapcsed
             !   ENDDO
             !   c_sedtot(k,i,j)=csednew
             !   dzs(k,i,j)=dzs(k,i,j)/rapcsed
             !ENDIF 
 
             ksmaxa=ksmax
             DO k=ksmaxa,ksmin+1,-1
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   !!!  FUSION OF 2 LAYERS IF TOO FINE OR TOO LITTLE MATERIAL        !!!!!
                   !!! (HYP : csed_tot of bottom layer is unchanged )                !!!!!
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef key_MUSTANG_V2
                !dzsmin=(1.0_rsh-coeff_dzsmin)*dzsminuni    &
                !               +coeff_dzsmin*SUM( (cv_sed(1:nvpc,k,i,j)/c_sedtot(k,i,j))*diam_sed(1:nvpc) )
                frac_sed(1:nvpc)=cv_sed(1:nvpc,k,i,j)/c_sedtot(k,i,j)
                dzsmin=dzsminvar(frac_sed)

#endif

                IF(c_sedtot(k,i,j) .LE. csedmin .OR. dzs(k,i,j) .LE. dzsmin)THEN

                  ! state of the layer before fusion
                  dzsa=dzs(k,i,j)
                  dzsam1=dzs(k-1,i,j)
                  poroam1=poro(k-1,i,j)
                  poroa=poro(k,i,j)
              
                  ! csed_tot(k-1) Cste : l epaisseur de la couche k-1 n est donc pas la somme des deux epaisseurs.
                  !dzs(k-1,i,j)=(dzsam1*c_sedtot(k-1,i,j)+c_sedtot(k,i,j)*dzsa)   &
                  !             /c_sedtot(k-1,i,j)

                  ! ou autre hypothese pour conserver eau et variables dissoutes
                  ! csed_tot non cste mais volume d eau et quantite de particules restent constantes
                  dzs(k-1,i,j)=dzsam1+dzsa
                  dzsisu=1.0_rsh/dzs(k-1,i,j)
              
                  c_sedtot(k-1,i,j)=0.0_rsh
                  cvolp=0.0_rsh
                  DO iv=1,nvp
                    cv_sed(iv,k-1,i,j)=(cv_sed(iv,k,i,j)*dzsa+                           &
                                    cv_sed(iv,k-1,i,j)*dzsam1)*dzsisu
                    c_sedtot(k-1,i,j)=c_sedtot(k-1,i,j)+cv_sed(iv,k-1,i,j)*typart(iv)
                    cvolp=cvolp+cv_sed(iv,k-1,i,j)*typart(iv)/ros(iv)
                  ENDDO
                  poro(k-1,i,j)=1.0_rsh-cvolp
#if ! defined key_noTSdiss_insed
                  cv_sed(-1,k-1,i,j)=(cv_sed(-1,k,i,j)*dzsa+cv_sed(-1,k-1,i,j)*dzsam1)*dzsisu
#endif
#ifdef key_MUSTANG_V2
                  sommud=0.0_rsh
                  DO iv=imud1,imud2
                    sommud=sommud+cv_sed(iv,k-1,i,j)
                  ENDDO
                  cvolgrvsan=0.0_rsh
                  DO iv=igrav1,isand2
                    cvolgrvsan=cvolgrvsan+cv_sed(iv,k-1,i,j)/ros(iv)
                  ENDDO
                  crel_mud(k-1,i,j)=sommud/(1.0_rsh-cvolgrvsan)
                  poro_mud(k-1,i,j)=1.0_rsh-sommud/ros(1)
#endif
                  volak=poroa*dzsa
                  volakm1=poroam1*dzsam1
                  vola=volak+volakm1
              
                 ! put in comment because useless
                  ! if(porosity too small (<32%), it is increased until (1.0_rsh-cvolmaxmel):
                  !IF(poro(k-1,i,j) < 1.0_rsh-cvolmaxmel-0.01_rsh)THEN
                  !  poro(k-1,i,j)=1.0_rsh-cvolmaxmel
                !!!! ATTENTION : on suppose ros identique quelque soit les particules
                  !  csednew=ros(1)*cvolmaxmel
                  !  rapcsed=csednew/c_sedtot(k-1,i,j)
                  !  DO iv=1,nvp
                  !    cv_sed(iv,k-1,i,j)=cv_sed(iv,k-1,i,j)*rapcsed
                  !  ENDDO
                  !  c_sedtot(k-1,i,j)=csednew
                  !  dzs(k-1,i,j)=dzs(k-1,i,j)/rapcsed
                  !ENDIF   
             
                 ! calcul des flux d eau et de matiere dissoute dus a la fusion des deux couches
                 ! dvol_fusion  en m3/m2; 
                  !dvol=vola-poro(k-1,i,j)*dzs(k-1,i,j)
                  !dvol_fusion=dvol +dvol_fusion            

#if ! defined key_noTSdiss_insed
                  !! ici pas d appel a routine MUSTANGV2_eval_dissvar_IWSflux 
                  !! on fait pareil qu en V1
                  ! calculation of new concentrations of dissolved variables  in pore water
                  ! in old water volume => mixed in layer k-1
                  vola_inv=1._rsh/vola

                  ! iv=0 : salinity (temperature is treated elsewhere
                  iv=0
                  cv_sed(iv,k-1,i,j)=vola_inv*(cv_sed(iv,k-1,i,j)*volakm1 + &
                                       cv_sed(iv,k,i,j)*volak)
                  !  flu_fusion(iv)=flu_fusion(iv)+cv_sed(iv,k-1,i,j)*dvol
#if ! defined key_Pconstitonly_insed
                  DO iv=nvp+1,nv_use2
                    cv_sed(iv,k-1,i,j)=vola_inv*(cv_sed(iv,k-1,i,j)*volakm1 + &
                                   cv_sed(iv,k,i,j)*volak)
                   ! flu_fusion(iv)=flu_fusion(iv)+cv_sed(iv,k-1,i,j)*dvol
                  ENDDO
#endif
#endif

                  volpwak(k-1)=volpwak(k-1)+volpwak(k)
              
                  DO kp=k,ksmax-1
                    dzs(kp,i,j)=dzs(kp+1,i,j)
                    volpwak(kp)=volpwak(kp+1)
                    DO iv=0,nv_use2
                      cv_sed(iv,kp,i,j)=cv_sed(iv,kp+1,i,j)
                    ENDDO
                    c_sedtot(kp,i,j)=c_sedtot(kp+1,i,j)
                    poro(kp,i,j)=poro(kp+1,i,j)
#ifdef key_MUSTANG_V2
                    crel_mud(kp,i,j)=crel_mud(kp+1,i,j)
                    poro_mud(kp,i,j)=poro_mud(kp+1,i,j)
#endif
                   ENDDO
                  dzs(ksmax:ksdmax,i,j)=0.0_rsh
                  c_sedtot(ksmax:ksdmax,i,j)=-valmanq
                  cv_sed(0:nv_use2,ksmax:ksdmax,i,j)=-valmanq
                  poro(ksmax:ksdmax,i,j)=0.0_rsh
#ifdef key_MUSTANG_V2
                  crel_mud(ksmax:ksdmax,i,j)=0.0_rsh
                  poro_mud(ksmax:ksdmax,i,j)=0.0_rsh
#endif
                  ksmax=ksmax-1
                ENDIF  
                ! end fusion

             ENDDO        
             
              ! updating ksma (ksmax can be modified in the routine)
             ksmaxa=ksma(i,j)
             ksma(i,j)=ksmax
          ENDIF            
          ! end of l_consolid

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!  PARTICULATE BIOTURBATION     !!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                                      
          ! explicit scheme by fractional step :  bioturbation treated as a diffusion mixing
          IF(l_bioturb .AND. ksma(i,j) > ksmi(i,j)) THEN
         
            ! again estimation of bioturb and biodiff coefs if sed layers number has changed during diter
            IF(ksmax .NE. ksmax_turb_part) THEN
                   CALL sed_MUSTANG_coefbioturb_part(i,j,difbio) 
                ksmax_turb_part=ksmax
            ENDIF
            cvsednew(:,:)=0.0_rsh
            k=ksmax
            dtsdzs=REAL(dtiter,rsh)/dzs(k,i,j)
            dtsdzsmin=REAL(dtiter,rsh)/MIN(dzs(k,i,j),dzs(k-1,i,j))
            IF(dtsdzsmin>0.5_rsh*dzs(k,i,j)/(difbio(k-1,imud1)+1.e-10_rsh) .AND. difbio(k-1,imud1) > 1.e-13) THEN
            ! instant mixing
              DO iv=isand1,nvp
                cv_sed(iv,k,i,j)= (cv_sed(iv,k,i,j)*dzs(k,i,j)+cv_sed(iv,k-1,i,j)*dzs(k-1,i,j)) &
                                   /(dzs(k,i,j)+dzs(k-1,i,j))
                cv_sed(iv,k-1,i,j)= cv_sed(iv,k,i,j) 
                cvsednew(iv,k)= cv_sed(iv,k,i,j)                                                     
                cvsednew(iv,k)= cv_sed(iv,k-1,i,j)                                                     
              ENDDO
#if ! defined key_noTSdiss_insed
              cv_sed(-1,k,i,j)= (cv_sed(-1,k,i,j)*dzs(k,i,j)+cv_sed(-1,k-1,i,j)*dzs(k-1,i,j)) &
                                 /(dzs(k,i,j)+dzs(k-1,i,j))
              cv_sed(-1,k-1,i,j)= cv_sed(-1,k,i,j) 
              cvsednew(nvp+1,k)= cv_sed(-1,k,i,j)                                                     
              cvsednew(nvp+1,k)= cv_sed(-1,k-1,i,j)
#endif                                                     
            ELSE
              DO iv=isand1,nvp
                cvsednew(iv,k)= cv_sed(iv,k,i,j) + &
                                dtsdzs*difbio(k-1,iv)*dzsiin(k-1)*(cv_sed(iv,k-1,i,j)-cv_sed(iv,k,i,j))                                     
              ENDDO
#if ! defined key_noTSdiss_insed
              cvsednew(nvp+1,k)= cv_sed(-1,k,i,j) + &
                             dtsdzs*difbio(k-1,iv)*dzsiin(k-1)*(cv_sed(-1,k-1,i,j)-cv_sed(-1,k,i,j))
#endif                                     
            ENDIF
            DO k=ksmax-1,ksmin+1,-1
              dtsdzs=REAL(dtiter,rsh)/dzs(k,i,j)
              dtsdzsmin=REAL(dtiter,rsh)/MIN(dzs(k,i,j),dzs(k-1,i,j),dzs(k+1,i,j))
              IF(dtsdzsmin>0.5_rsh*dzs(k,i,j)/(difbio(k,imud1)+1.e-10_rsh).AND. difbio(k,imud1) > 1.e-13) THEN
            ! instant mixing
               DO iv=isand1,nvp
                  cv_sed(iv,k,i,j)= (cv_sed(iv,k,i,j)*dzs(k,i,j)+cv_sed(iv,k-1,i,j)*dzs(k-1,i,j)+ &
                          cv_sed(iv,k+1,i,j)*dzs(k+1,i,j))/(dzs(k,i,j)+dzs(k-1,i,j)+dzs(k+1,i,j))
                  cv_sed(iv,k-1,i,j)= cv_sed(iv,k,i,j) 
                  cv_sed(iv,k+1,i,j)= cv_sed(iv,k,i,j) 
                  cvsednew(iv,k)= cv_sed(iv,k,i,j)                                                     
                  cvsednew(iv,k-1)= cv_sed(iv,k-1,i,j)                                                     
                  cvsednew(iv,k+1)= cv_sed(iv,k+1,i,j)                                                     
               ENDDO
#if ! defined key_noTSdiss_insed
               cv_sed(-1,k,i,j)= (cv_sed(-1,k,i,j)*dzs(k,i,j)+cv_sed(-1,k-1,i,j)*dzs(k-1,i,j)+ &
                         cv_sed(-1,k+1,i,j)*dzs(k+1,i,j))/(dzs(k,i,j)+dzs(k-1,i,j)+dzs(k+1,i,j))
               cv_sed(-1,k-1,i,j)= cv_sed(-1,k,i,j) 
               cv_sed(-1,k+1,i,j)= cv_sed(-1,k,i,j) 
               cvsednew(nvp+1,k)= cv_sed(-1,k,i,j)                                                     
               cvsednew(nvp+1,k-1)= cv_sed(-1,k-1,i,j)                                                     
               cvsednew(nvp+1,k+1)= cv_sed(-1,k+1,i,j)  
#endif                                                   
              ELSE
                DO iv=isand1,nvp
                  cvsednew(iv,k)= cv_sed(iv,k,i,j) + dtsdzs*(   &
                                difbio(k,iv)*dzsiin(k)*cv_sed(iv,k+1,i,j)  &
                              -(difbio(k,iv)*dzsiin(k)+difbio(k-1,iv)*dzsiin(k-1))*cv_sed(iv,k,i,j) &
                              + difbio(k-1,iv)*dzsiin(k-1)*cv_sed(iv,k-1,i,j))                                     
                ENDDO 
#if ! defined key_noTSdiss_insed             
                cvsednew(nvp+1,k)= cv_sed(-1,k,i,j) + dtsdzs*(   &
                                difbio(k,iv)*dzsiin(k)*cv_sed(-1,k+1,i,j)  &
                              -(difbio(k,iv)*dzsiin(k)+difbio(k-1,iv)*dzsiin(k-1))*cv_sed(-1,k,i,j) &
                              + difbio(k-1,iv)*dzsiin(k-1)*cv_sed(-1,k-1,i,j))                                     
#endif
              ENDIF
            ENDDO              
            k=ksmin
            dtsdzs=REAL(dtiter,rsh)/dzs(k,i,j)    
            dtsdzsmin=REAL(dtiter,rsh)/MIN(dzs(k,i,j),dzs(k+1,i,j))
            IF(dtsdzsmin>0.5_rsh*dzs(k,i,j)/(difbio(k,imud1)+1.e-10_rsh).AND. difbio(k,imud1) > 1.e-13) THEN
            ! instant mixing
              DO iv=isand1,nvp
                cv_sed(iv,k,i,j)= (cv_sed(iv,k,i,j)*dzs(k,i,j)+cv_sed(iv,k+1,i,j)*dzs(k+1,i,j))&
                                   /(dzs(k,i,j)+dzs(k+1,i,j))
                cv_sed(iv,k+1,i,j)= cv_sed(iv,k,i,j) 
                cvsednew(iv,k)= cv_sed(iv,k,i,j)                                                     
                cvsednew(iv,k+1)= cv_sed(iv,k+1,i,j)                                                     
              ENDDO
              cv_sed(-1,k,i,j)= (cv_sed(-1,k,i,j)*dzs(k,i,j)+cv_sed(-1,k+1,i,j)*dzs(k+1,i,j))&
                                   /(dzs(k,i,j)+dzs(k+1,i,j))
              cv_sed(-1,k+1,i,j)= cv_sed(-1,k,i,j) 
              cvsednew(nvp+1,k)= cv_sed(-1,k,i,j)                                                     
              cvsednew(nvp+1,k+1)= cv_sed(-1,k+1,i,j)                                                     
            ELSE
              DO iv=isand1,nvp
               cvsednew(iv,k)= cv_sed(iv,k,i,j) + &
                                dtsdzs*difbio(k,iv)*dzsiin(k)*(cv_sed(iv,k+1,i,j)-cv_sed(iv,k,i,j))                                     
              ENDDO              
              cvsednew(nvp+1,k)= cv_sed(-1,k,i,j) + &
                                dtsdzs*difbio(k,iv)*dzsiin(k)*(cv_sed(-1,k+1,i,j)-cv_sed(-1,k,i,j))                                     
            ENDIF
            DO k=ksmin,ksmax
              cv_sed(isand1:nvp,k,i,j)=cvsednew(isand1:nvp,k)
              cv_sed(-1,k,i,j)=cvsednew(nvp+1,k)
            ENDDO              
            ! computing new concentrations/porosities :
            DO k=ksmin,ksmax
              c_sedtot(k,i,j)=0.0_rsh
              cvolp=0.0_rsh
              DO iv=1,nvpc
                c_sedtot(k,i,j)=c_sedtot(k,i,j)+cv_sed(iv,k,i,j)
                cvolp=cvolp+cv_sed(iv,k,i,j)/ros(iv)
              ENDDO
              poro(k,i,j)=1.0_rsh-cvolp
#ifdef key_MUSTANG_V2
              sommud=0.0_rsh
              DO iv=imud1,imud2
                sommud=sommud+cv_sed(iv,k,i,j)
              ENDDO
              cvolgrvsan=0.0_rsh
              DO iv=igrav1,isand2
                cvolgrvsan=cvolgrvsan+cv_sed(iv,k,i,j)/ros(iv)
              ENDDO
              crel_mud(k,i,j)=sommud/(1.0_rsh-cvolgrvsan)
              poro_mud(k,i,j)=1.0_rsh-sommud/ros(1)
#endif
            ENDDO
             ! mixing of particulate leads to dissolved mixing (taking into account further with winter)
          ENDIF   
          ! end of l_bioturb

         ELSE
            dtiter=dtrest
         ENDIF
         ! several possible steps of consolidation or bioturbation
         dtrest=dtrest-dtiter
         IF(dtrest.GE.1.0_rlg)GOTO 100

#if ! defined key_noTSdiss_insed && ! defined key_nofluxwat_IWS
         IF((l_consolid .OR. l_bioturb)) THEN

          ! Exchange flux at the water/sediment interface
          ! and intertitial water vertical rates
          ! ---------------------------------------------
          !  fluconsol is flux resulting from consolidation, >0 upwards, centered at the upper interface
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
          !!!! interstitial water vertical rates  : winter (m3/m2/s) k :between k and k+1
          !!!! total flux of dissolved substance due to fusion (melting) : M/m2/s
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
          winter(ksmin-1)=0.0_rsh
          volpwnv_tot=0.0_rsh
          DO k=ksmin,ksmax
           volpwk=poro(k,i,j)*dzs(k,i,j)
           volpwnv_tot=volpwnv_tot+volpwk
           winter(k)=winter(k-1)+(volpwak(k)-volpwk)*dt_sed_inv
          ENDDO
         
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
          !!!! water fluxes between water and sediment due to consolidation and fusion (m3/s)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
#if ! defined key_nofluxwat_IWS
          phieau_s2w_consol(i,j)=(volpwa_tot-volpwnv_tot)*CELL_SURF(i,j)
#endif
#if ! defined key_noTSdiss_insed
          !  flu_dyninsed is the flux due to consolidation at the interface, >0 upwards (M/m2)
          flu_dyninsed(0,i,j)=winter(ksmax)*cv_sed(0,ksmax,i,j)*REAL(dt_sed_eff,rsh) 
          flu_dyninsed(nvp+1:nv_use2,i,j)= &
                        winter(ksmax)*cv_sed(nvp+1:nv_use2,ksmax,i,j)*REAL(dt_sed_eff,rsh) 
                
#endif
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                      
          !!!!!!!  END consolidation -  bioturbation                 !!!
          !!!!!!!    update thicknesses and intermediate porosities  !!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                      
           
          poroin(ksma(i,j))=poro(ksma(i,j),i,j)
          DO k=ksmi(i,j),ksma(i,j)-1
              dzsiin(k)=2.0_rsh/(dzs(k,i,j)+dzs(k+1,i,j))
              poroin(k)=0.5_rsh*(poro(k+1,i,j)+poro(k,i,j))          
          ENDDO
         ENDIF  ! end test l_consolid or l_bioturb
#endif
        
        !print *,' --> fin consolidation : poro(ksmax,i,j)=',poro(ksmax,i,j)
        !print *,' --> fin consolidation : poro(ksmax-1,i,j)=',poro(ksmax-1,i,j)
#if defined key_MUSTANG_specif_outputs 
        DO k=ksmi(i,j),ksma(i,j)
          varspecif3Dk_save(1,k,i,j)=poro(k,i,j)       ! poro_save
        END DO
#endif
        
#if ! defined key_noTSdiss_insed
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!  DIFFUSION -BIOTURBATION DES VARIABLES DISSOUTES       !!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF(l_biodiffs .AND. ksmax .NE. ksmax_turb_diss) THEN
            ! calculation of biodiffusion coef to be updated if the number of layers has changed during diter 
            CALL sed_MUSTANG_coefbioturb_diss(i,j,difbio) 
            ksmax_turb_diss=ksmax
         ENDIF

         hcrit=MAX(0.0_rsh,htot(i,j)-h0fond)
         disvi(:,:)=0.0_rsh
         conc_bottom(:)=0.0_rsh
         
              ! calculation of diffusion / biodiffusion rates in the sediment
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF(l_diffused .OR. l_biodiffs) THEN

           DO ivv=0,nv_use2-nvp
              iv=ivdiss(ivv)

              k=ksmax
              IF(choice_flxdiss_diffsed == 1) THEN
                disvi(k,ivv)=2.0_rsh*xdifsi1/(dzs(k,i,j)+epn_bottom_MUSTANG(i,j))
              ELSE IF(choice_flxdiss_diffsed == 2) THEN
                disvi(k,ivv)=xdifsi1/(dzs(k,i,j)*0.5_rsh+epdifi)
              ELSE 
               ! in the boundary layer : Boudreau formulation  (1997)
               ! diffusion rate= D/deltae=beta (page 182) eq 5.55 a 5.56)
               ! and beta =0.0889 u* Sc^(-0.704)
               ! formula of Shaw and Hanratty (1977) eq 5.58 in Boudreau 
               ! Sc  Schmidt number = kinematic viscosity / mass diffusivity
               ! of the order of 1000 but may vary with the substance  !!!!!!! TO REVIEW!!!!!!!!!

               ! disvi(k,ivv)=0.0889*ustarbot(i,j)*(1000._rsh)**(-0.704)
                IF(iv<1) THEN     ! Sc kept = 1000 for  Temp and Sal
                  Sc = 1000.0_rsh
                ELSE              ! Modif Martin, Schmidt number for all subs excepted T and S : Sc= nu/D (Sideman & Pinczewski,1975)
                  p = 0.0_rsh     ! No impact of pressure (for now)
                  ! mu = dynamic viscosity in g/cm/s (Kulkula et al. 1987, in Boudreau p.94)
                  mu = 0.01*(1.791_rsh - 6.144e-02_rsh*temp_bottom_MUSTANG(i,j) + 1.451e-03_rsh*temp_bottom_MUSTANG(i,j)**2 &
                       - 1.6826e-05_rsh*temp_bottom_MUSTANG(i,j)**3 - 1.529e-04_rsh*p + 8.3885e-08_rsh*p*p &
                       + 2.4727e-03_rsh*sal_bottom_MUSTANG(i,j) + temp_bottom_MUSTANG(i,j)*(6.0574e-06_rsh*p - 2.676e-09_rsh*p*p) &
                       + sal_bottom_MUSTANG(i,j)*(4.8429e-05_rsh*temp_bottom_MUSTANG(i,j) - 4.7172e-06_rsh*temp_bottom_MUSTANG(i,j)**2 &
                       + 7.5986e-08_rsh*temp_bottom_MUSTANG(i,j)**3))
                  nu = mu*rowinv*1000                                ! 1/roro or rowinv in cm3/g and nu = cinematic viscosity in cm2/s
                  IF(D0_funcT_opt(iv) == 1) D0=(D0_m0(iv)+D0_m1(iv)*temp_bottom_MUSTANG(i,j))/1000000_rsh     ! D0 in cm2/s    
                  IF(D0_funcT_opt(iv) == 2) D0=(D0_m0(iv)+D0_m1(iv)*(temp_bottom_MUSTANG(i,j)+273.15_rsh)/(mu*100.0_rsh))/100000_rsh
                  Sc = nu/D0
                ENDIF
                disvi(k,ivv)=0.0889*ustarbot(i,j)*Sc**(-0.704)
                !write(*,*)'Beta mass transfert coef at the interface =',disvi(k)
              ENDIF
              disvi(k,ivv)=disvi(k,ivv)*hcrit/(hcrit+epsilon_MUSTANG)
              DO k=ksmin,ksmax-1
                IF(iv<1) THEN     ! xdifd1 unchanged for Temp and Sal for now
               ! formulation according to the tortuosity Dsed=Dpure/Tortuosite^2
               ! tortuosity function of porosity (eq 4.120)
               ! in Boudreau 1997 p 132 ! tortuosity^2=1-log(poro^2)
               !disvi(k,ivv)=xdifs1*dzsiin(k,i,j)/(1-LOG(poro(k,i,j)**2))
                  disvi(k,ivv)=(xdifs1/(1.0_rsh-2.0_rsh*LOG(poro(k,i,j)))+difbio(k,iv))*dzsiin(k)
                ELSE              ! Modif Martin, diffusion for all substances excepted T,S
                  p = 0.0_rsh     ! No impact of the pressure (for now)  
                  ! mu = dynamic viscosivity in g/cm/s (Kulkula et al. 1987, in Boudreau p.94)
                  mu = 0.01*(1.791_rsh - 6.144e-02_rsh*cv_sed(-1,k,i,j) + 1.451e-03_rsh*cv_sed(-1,k,i,j)**2 &
                       - 1.6826e-05_rsh*cv_sed(-1,k,i,j)**3 - 1.529e-04_rsh*p + 8.3885e-08_rsh*p*p &
                       + 2.4727e-03_rsh*cv_sed(0,k,i,j) + cv_sed(-1,k,i,j)*(6.0574e-06_rsh*p - 2.676e-09_rsh*p*p) &
                       + cv_sed(0,k,i,j)*(4.8429e-05_rsh*cv_sed(-1,k,i,j) - 4.7172e-06_rsh*cv_sed(-1,k,i,j)**2 &
                       + 7.5986e-08_rsh*cv_sed(-1,k,i,j)**3))
                  IF(D0_funcT_opt(iv) == 1) D0=(D0_m0(iv)+D0_m1(iv)*cv_sed(-1,k,i,j))/1000000_rsh     ! D0 in cm2/s    
                  IF(D0_funcT_opt(iv) == 2) D0=(D0_m0(iv)+D0_m1(iv)*(cv_sed(-1,k,i,j)+273.15_rsh)/(mu*100.0_rsh))/100000_rsh
                  xdifs1b = D0*0.94_rsh  ! convertion from 'infinite-dilution' diff into 'porewater' diff (Li & Gregory, 1974; in Boudreau p.125)
                  xdifs1b=xdifs1b/10000.0_rsh      ! from cm2/s to m2/s
                  disvi(k,ivv)=(xdifs1b/(1.0_rsh-2.0_rsh*LOG(poro(k,i,j)))+difbio(k,iv))*dzsiin(k)
                ENDIF
              ENDDO
           ENDDO

           conc_bottom(0)=sal_bottom_MUSTANG(i,j)
           DO ivv=1,nv_use2-nvp
               iv=ivv+nvp
               conc_bottom(ivv)=cw_bottom_MUSTANG(iv,i,j)
           ENDDO 

         ENDIF  !(diffused)

!
              !resolution of the advection / diffusion equation (advection by consolidation)
!
         IF(ksmax > ksmin)THEN
!
                !where there are several layers:
                !--------------------------------
                !coefficients of the matrix
!
           DO ivv=0,nv_use2-nvp
               iv=ivdiss(ivv)
                k=ksmin
                dzssdt=dzs(k,i,j)*dt_sed_inv
                aa(k)=0.0_rsh
                bb(k)=poro(k,i,j)*dzssdt+disvi(k,ivv)*poroin(k)+winter(k)*fexcs
                cc(k)=-disvi(k,ivv)*poroin(k)+winter(k)*cexcs
                dd(k)=cv_sed(iv,k,i,j)*volpwak(k)*dt_sed_inv
                DO k=ksmin+1,ksmax-1
                 dzssdt=dzs(k,i,j)*dt_sed_inv
                 aa(k)=-disvi(k-1,ivv)*poroin(k-1)-winter(k-1)*fexcs
                 bb(k)=poro(k,i,j)*dzssdt+disvi(k-1,ivv)*poroin(k-1)+disvi(k,ivv)*poroin(k) &
                                +winter(k)*fexcs-winter(k-1)*cexcs
                 cc(k)=-disvi(k,ivv)*poroin(k)+winter(k)*cexcs
                 dd(k)=cv_sed(iv,k,i,j)*volpwak(k)*dt_sed_inv
                ENDDO
                k=ksmax
                aa(k)=-disvi(k-1,ivv)*poroin(k-1)-winter(k-1)*fexcs
                bb(k)=dzs(k,i,j)*poro(k,i,j)*dt_sed_inv+disvi(k-1,ivv)*poroin(k-1)+disvi(k,ivv)+winter(k)-winter(k-1)*cexcs
                cc(k)=0.0_rsh
                fludif(iv,i,j)=disvi(k,ivv)*conc_bottom(ivv)
                dd(k)=cv_sed(iv,k,i,j)*volpwak(k)*dt_sed_inv+fludif(iv,i,j)
!
                ! resolution of the tridiagonal system
!
                dd(1)=dd(1)/bb(1)
                DO k=ksmin+1,ksmax
                  bb(k)=bb(k)-aa(k)*cc(k-1)/bb(k-1)
                  dd(k)=(dd(k)-aa(k)*dd(k-1))/bb(k)
                ENDDO
                cv_sed(iv,ksmax,i,j)=dd(ksmax)
                DO k=ksmax-1,ksmin,-1
                  cv_sed(iv,k,i,j)=dd(k)-cc(k)*cv_sed(iv,k+1,i,j)/bb(k)
                ENDDO
           ENDDO           

!
         ELSE
!
               !where there is only one layer (ksmax=ksmin):
               !--------------------------------------------       
           DO ivv=0,nv_use2-nvp
               iv=ivdiss(ivv)
                k=ksmax        
                fludif(iv,i,j)=disvi(k,ivv)*conc_bottom(ivv)   
                cv_sed(iv,k,i,j)=  &
                    (cv_sed(iv,k,i,j)*volpwak(k)+fludif(iv,i,j)*REAL(dt_sed_eff,rsh))    &
                   /(dzs(k,i,j)*poro(k,i,j)+(disvi(k,ivv)+winter(k))*REAL(dt_sed_eff,rsh))
           ENDDO
         ENDIF
!
              !  exchange flux at the water / sediment interface:
              !  ------------------------------------------
              !  fludif is the diffusion flux at the interface, >0 down  ==> M/m2
         IF(l_diffused .OR. l_biodiffs) THEN
           DO ivv=0,nv_use2-nvp
               iv=ivdiss(ivv)
               fludif(iv,i,j)=  &
                   (fludif(iv,i,j)-disvi(ksmax,ivv)*cv_sed(iv,ksmax,i,j))*REAL(dt_sed_eff,rsh)
           ENDDO                         ! end loop on  iv

         ENDIF                     ! test on iv ( var dissoutes seules)

#endif
        ENDIF    ! test if existence sediment
    
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!   END OF DYNAMIC PROCESS IN SEDIMENT                !!
       !!!!!!!!!!!!  CALCULATION OF WATER  FLUX                         !!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! si la cellule est mise a sec  ==> on cumule les apports d  eau et les flux de 
       !  substances dissoutes pour ne pas perdre de la matiere
       ! if the cell is dry ==> one cumulates the contributions of water and the flows of 
        ! dissolved substances so as not to lose matter
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF(htot(i,j) .LE. h0fond) THEN  
            ! no water
#if ! defined key_nofluxwat_IWS
            ! phieau_s2w_drycell in m3  
            phieau_s2w_drycell(i,j)=phieau_s2w_drycell(i,j)+REAL(phieau_s2w_consol(i,j),rlg)
#endif
#if ! defined key_noTSdiss_insed
            !fluconsol_drycell : cumul => M/m2
            fluconsol_drycell(0,i,j)=fluconsol_drycell(0,i,j) +flu_dyninsed(0,i,j)
            fluconsol_drycell(nvp+1:nv_use2,i,j)=fluconsol_drycell(nvp+1:nv_use2,i,j)  &
                                               +flu_dyninsed(nvp+1:nv_use2,i,j)
            fluconsol(0,i,j)=0.0_rsh
            fluconsol(nvp+1:nv_use2,i,j)=0.0_rsh                                    
#endif 
        ELSE
#if ! defined key_nofluxwat_IWS
           ! phieau_s2w in m3  
            phieau_s2w(i,j)=phieau_s2w(i,j)+phieau_s2w_drycell(i,j) + &
                                          REAL(phieau_s2w_consol(i,j),rlg)
            phieau_s2w_drycell(i,j)=0.0_rlg
#endif
#if ! defined key_noTSdiss_insed
            ! fluconsol  => M/m2 will be distributed in flx_s2w (M/m2/s)
            fluconsol(0,i,j)=flu_dyninsed(0,i,j) +fluconsol_drycell(0,i,j)
            fluconsol(nvp+1:nv_use2,i,j)=flu_dyninsed(nvp+1:nv_use2,i,j)  &
                                       +fluconsol_drycell(nvp+1:nv_use2,i,j)
            fluconsol_drycell(0,i,j)=0.0_rsh
            fluconsol_drycell(nvp+1:nv_use2,i,j)=0.0_rsh
#endif 
        ENDIF
  
      ENDDO
    ENDDO

     t_dyninsed=CURRENT_TIME+dt_sed_cor

    ELSE

            fluconsol(:,:,:)=0.0_rsh
            fludif(0:nv_adv,:,:)=0.0_rsh

    ENDIF
   
  END SUBROUTINE sed_MUSTANG_consol_diff_bioturb

!!==============================================================================
#ifdef key_MUSTANG_splitlayersurf
  SUBROUTINE sed_MUSTANG_split_surflayer(i,j,ksmax)
! 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_split_surflayer  ***
   !&E
   !&E ** Purpose : splitting surface layer if too thick 
   !&E
   !&E ** Description : can be used when one wants to always keep thin layers 
   !&E                  on the surface of the sediment, for biology for example.
   !&E                  nlayer_surf_sed : numbre of thin layers at sediment surface
   !&E                  dzsmax : max thickness at surface sediment (dzsmax_bottom below) 
   !&E
   !&E ** Called by :  MUSTANG_update
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used


   !! * Arguments
   INTEGER, INTENT(IN)            :: i,j
   INTEGER, INTENT(INOUT)         :: ksmax   
   
  !! * Local declarations
   INTEGER                          :: k,kk,ks_surf,nlayer,nn,ks_split_min,iv,nk
   INTEGER                          :: nlayer_a_fusion, nlayer_splitt
   INTEGER ,DIMENSION(ksdmax)       :: nlayer_splitk
   REAL(KIND=rsh)                   :: dzsa,poroa,cvolp,zz
   REAL(KIND=rsh)                   :: thick_surf_new ,thick_surf_old,stcvsed_surf_new,stcvsed_surf_old
   REAL(KIND=rsh)                   :: stcsedtot_new,stcsedtot_old
#if ! defined key_noTSdiss_insed
!#ifdef key_MUSTANG_V2
!   REAL(KIND=rsh)                   :: porowater_new,porowater1,porowater2
!#else
   REAL(KIND=rsh)                   :: dzspwi
!#endif
#endif


   !! * executable statement
   !!   --------------------
             ! splitting : evaluating number of layers to be splitting (ks_split_min to ks_surf)
            !                    and number of new layers in each layer to be splitting (nayler_split(k))
            !                    and number of new layers at the end of splitting (nayler)
         ks_surf=ksmax
         IF(ks_surf > ksmi(i,j))THEN
           zz=0.0_rsh
           k=ks_surf
           nlayer_splitk(:)=0
           nlayer_splitt=0
           DO WHILE (zz < (nlayer_surf_sed*dzsmax(i,j)-dzsmin)  &
                     .AND. ks_surf-k < nlayer_surf_sed .AND. k> ksmi(i,j))
             zz=zz+dzs(k,i,j)
             !  nombre de couches a prevoir  (number of layers to predict)
             nlayer_splitt=nlayer_splitt+INT(dzs(k,i,j)/(dzsmax(i,j)+ 5.0_rsh*dzsmin))
             k=k-1
           ENDDO
           IF(nlayer_splitt > 0)THEN
            !split yes
             ks_split_min=MAX(2,k)
             ! number of layers limited by the max depth chosen for the well-discrete surface layer
             nlayer_splitt=nlayer_splitt-INT((zz-nlayer_surf_sed*dzsmax(i,j))/(dzsmax(i,j)+ 5.0_rsh*dzsmin))+1
             IF(nlayer_splitt > ksdmax-ks_surf)THEN
             !  melting at the bottom because not enough layers for split
                nlayer_a_fusion=nlayer_splitt-(ksdmax-ks_surf)
                DO k=1,nlayer_a_fusion
                  CALL sed_MUSTANG_fusion(i,j,ks_surf)
                ENDDO
                ! everything is shifted: ks_surf and ks_split_min go down if merge at the bottom
                ks_split_min=ks_split_min-nlayer_a_fusion
             ENDIF
               
#ifdef key_test_conservativity_splitsurf
             ! to test conservativity
             thick_surf_old=0.0_rsh
             !stcvsed_surf_old=0.0_rsh
             stcsedtot_old=0.0_rsh
             DO kk=MAX(1,ks_split_min-2),ks_surf
               thick_surf_old=thick_surf_old+dzs(kk,i,j)
              ! stcvsed_surf_old=stcvsed_surf_old+cv_sed(3,kk,i,j)*dzs(kk,i,j)
               stcsedtot_old=stcsedtot_old+c_sedtot(kk,i,j)*dzs(kk,i,j)
             ENDDO
#endif
             
             ! splitting surface layer first above ks_surf
             k=ks_surf  
             ! first layer at surface , split only if dzs>dzsmax+10 dzsmin
             ! number of layers to split (not counting the ks_surf layer itself)
             nlayer_splitk(k)=MIN(nlayer_surf_sed,INT(dzs(ks_surf,i,j)/(dzsmax(i,j)+5._rsh*dzsmin)))
             ! split of layers ks_surf+1 to ks_surf+nlayer_split 
             DO kk=k+1,k+nlayer_splitk(k)
                 dzs(kk,i,j)=dzsmax(i,j)
                 cv_sed(:,kk,i,j)=cv_sed(:,k,i,j)
                 c_sedtot(kk,i,j)=c_sedtot(k,i,j)
                 poro(kk,i,j)=poro(k,i,j)
             ENDDO
             ! calculating the residual thickness of the previously surface layer  
             dzs(ks_surf,i,j)=dzs(ks_surf,i,j)-nlayer_splitk(k)*dzsmax(i,j)

             ! number of added layers 
             nlayer_splitt=nlayer_splitk(k)
             ! new ksmax
             ksmax=ks_surf+nlayer_splitk(k)
             ! loop on the extra layers below that are to split 
             DO k=ks_surf-1,ks_split_min,-1

              ! if the layer above has not been splited because too thin, we go
              IF(nlayer_splitk(k+1) > 0 ) THEN
                 ! we test the thickness of the layer above to merge it if it is too thin
                 IF( dzs(k+1,i,j) < dzsmin ) THEN
                    ! fusion of above layer (k+1) which was splitting before and layer to be split (k)
                    dzsa=dzs(k,i,j)
                    poroa=poro(k,i,j)
                    dzs(k,i,j)=dzs(k,i,j)+dzs(k+1,i,j)
                    cvolp=0.0_rsh
                    c_sedtot(k,i,j)=0.0_rsh
                    DO iv=1,nvp
                     cv_sed(iv,k,i,j)=(cv_sed(iv,k,i,j)*dzsa+cv_sed(iv,k+1,i,j)*dzs(k+1,i,j))/dzs(k,i,j)
                     c_sedtot(k,i,j)=c_sedtot(k,i,j)+cv_sed(iv,k,i,j)*typart(iv)
                     cvolp=cvolp+cv_sed(iv,k,i,j)*typart(iv)/ros(iv)
                    ENDDO
                    poro(k,i,j)=1.0_rsh-cvolp
                    
#if ! defined key_noTSdiss_insed
                    ! dissolved variable in pore waters and water fluxes at the interface
     ! --------------------------------------------------------------------
!#ifdef key_MUSTANG_V2
  !!! enleve pour eviter d apeler une routine alors que peu de calcul
                    ! fusion of 2 layers ==> code 2
 !                   porowater_new=poro(k,i,j)*dzs(k,i,j)
 !                   porowater1=dzsa*poroa
 !                   porowater2=dzs(k+1,i,j)*poro(k+1,i,j)
 !                   CALL MUSTANGV2_eval_dissvar_IWSflux(i,j,ksmax,2,  &
 !                         porowater1=porowater1,                    &
 !                         cv_sed1=cv_sed(-1:nv_adv,k,i,j),          &
 !                         porowater2=porowater2,                    &
 !                         cv_sed2=cv_sed(-1:nv_adv,k+1,i,j),        &
 !                         porowater_new=porowater_new)
!#else
  ! MUSTANG_V1  et ici V2
                    dzspwi=1.0_rsh/(dzs(k,i,j)*poro(k,i,j))
                    DO iv=-1,0
                     cv_sed(iv,k,i,j)=dzspwi*(cv_sed(iv,k,i,j)*dzsa*poroa &
                           +cv_sed(iv,k+1,i,j)*dzs(k+1,i,j)*poro(k+1,i,j))
                    ENDDO
#if ! defined key_Pconstitonly_insed
                    DO iv=nvp+1,nv_adv
                     cv_sed(iv,k,i,j)=dzspwi*(cv_sed(iv,k,i,j)*dzsa*poroa &
                                +cv_sed(iv,k+1,i,j)*dzs(k+1,i,j)*poro(k+1,i,j))
                    ENDDO
#endif
#endif

                    DO kk=k+1,ksmax-1
                         dzs(kk,i,j)=dzs(kk+1,i,j)
                         poro(kk,i,j)=poro(kk+1,i,j)
                         c_sedtot(kk,i,j)=c_sedtot(kk+1,i,j)
                         DO iv=-1,nv_adv
                          cv_sed(iv,kk,i,j)=cv_sed(iv,kk+1,i,j)
                         ENDDO
#ifdef key_MUSTANG_V2
                         crel_mud(kk,i,j)=crel_mud(kk+1,i,j)
                         poro_mud(kk,i,j)=poro_mud(kk+1,i,j)
#endif
                    ENDDO
                    ksmax=ksmax-1 
                    nlayer_splitk(k+1)=nlayer_splitk(k+1)-1
                  ENDIF
              ENDIF
              ! number of  layers to split in layer k, limited if the number of layers nlayer_surf_sed is reached 
              nlayer_splitk(k)=MIN(nlayer_surf_sed-(ksmax-k-1),INT(dzs(k,i,j)/(dzsmax(i,j)+5.0_rsh*dzsmin)))
              nlayer_splitk(k)=MAX(0,nlayer_splitk(k))
              IF(nlayer_splitk(k) > 0) THEN
                ! moving all surface layers above for addition new splitting layers issued from k layer
                ! toutes les couches au dessus de la couche k sont decalees vers le haut pour faire la place aux nouvelles couches
                dzs(k+nlayer_splitk(k)+1:ksmax+nlayer_splitk(k),i,j)=dzs(k+1:ksmax,i,j)
                cv_sed(:,k+nlayer_splitk(k)+1:ksmax+nlayer_splitk(k),i,j)=cv_sed(:,k+1:ksmax,i,j)
                c_sedtot(k+nlayer_splitk(k)+1:ksmax+nlayer_splitk(k),i,j)=c_sedtot(k+1:ksmax,i,j)
                poro(k+nlayer_splitk(k)+1:ksmax+nlayer_splitk(k),i,j)=poro(k+1:ksmax,i,j)
                nlayer_splitt=nlayer_splitt+nlayer_splitk(k)
  
                ! splitting layer k
                DO kk=k+1,k+nlayer_splitk(k)
                  dzs(kk,i,j)=dzsmax(i,j)
                  cv_sed(:,kk,i,j)=cv_sed(:,k,i,j)
                  c_sedtot(kk,i,j)=c_sedtot(k,i,j)
                  poro(kk,i,j)=poro(k,i,j)
#ifdef key_MUSTANG_V2
                  crel_mud(kk,i,j)=crel_mud(k,i,j)
                  poro_mud(kk,i,j)=poro_mud(k,i,j)
#endif
                ENDDO
                ksmax=ksmax+nlayer_splitk(k)
              
                ! rest of k layer thickness 
                dzs(k,i,j)=dzs(k,i,j)-nlayer_splitk(k)*dzsmax(i,j)
              ENDIF

             ENDDO  
            ks_split_min=ksmax-nlayer_surf_sed

            ! fusion of the last split layer with the layer underneath (elimination of the remainder of the last split layer)
             IF(dzs(ks_split_min,i,j) < dzsmax(i,j)*0.5_rsh .AND. ks_split_min > ksmi(i,j)) THEN
               k=ks_split_min-1
               dzsa=dzs(k,i,j)
               poroa=poro(k,i,j)
               dzs(k,i,j)=dzs(k,i,j)+dzs(k+1,i,j)
               cvolp=0.0_rsh
               c_sedtot(k,i,j)=0.0_rsh
               DO iv=1,nvp
                 cv_sed(iv,k,i,j)=(cv_sed(iv,k,i,j)*dzsa+cv_sed(iv,k+1,i,j)*dzs(k+1,i,j))/dzs(k,i,j)
                 c_sedtot(k,i,j)=c_sedtot(k,i,j)+cv_sed(iv,k,i,j)*typart(iv)
                 cvolp=cvolp+cv_sed(iv,k,i,j)*typart(iv)/ros(iv)
               ENDDO
               poro(k,i,j)=1.0_rsh-cvolp
               !!  ici on fait pareil en V2, pas d appel a routine MUSTANGV2_eval_dissvar_IWSflux
#if ! defined key_noTSdiss_insed
               dzspwi=1.0_rsh/(dzs(k,i,j)*poro(k,i,j))
               DO iv=-1,0
                 cv_sed(iv,k,i,j)=dzspwi*(cv_sed(iv,k,i,j)*dzsa*poroa &
                           +cv_sed(iv,k+1,i,j)*dzs(k+1,i,j)*poro(k+1,i,j))
               ENDDO
#if ! defined key_Pconstitonly_insed
               DO iv=nvp+1,nv_adv
                 cv_sed(iv,k,i,j)=dzspwi*(cv_sed(iv,k,i,j)*dzsa*poroa &
                                +cv_sed(iv,k+1,i,j)*dzs(k+1,i,j)*poro(k+1,i,j))
               ENDDO
#endif
#endif
              ! if fusion, elimination of the ks_split_min layer and moving the others
               dzs(ks_split_min:ksmax-1,i,j)=dzs(ks_split_min+1:ksmax,i,j)
               cv_sed(:,ks_split_min:ksmax-1,i,j)=cv_sed(:,ks_split_min+1:ksmax,i,j)
               c_sedtot(ks_split_min:ksmax-1,i,j)=c_sedtot(ks_split_min+1:ksmax,i,j)
               poro(ks_split_min:ksmax-1,i,j)=poro(ks_split_min+1:ksmax,i,j)
               ksmax=ksmax-1
               dzs(ksmax+1:ksdmax,i,j)=0.0_rsh
               cv_sed(:,ksmax+1:ksdmax,i,j)=-valmanq
               c_sedtot(ksmax+1:ksdmax,i,j)=-valmanq
               poro(ksmax+1:ksdmax,i,j)=0.0_rsh
#ifdef key_MUSTANG_V2
               crel_mud(ksmax+1:ksdmax,i,j)=0.0_rsh
               poro_mud(ksmax+1:ksdmax,i,j)=0.0_rsh
#endif

             ENDIF
             
             ! decoupage de la  couche dessous si trop epaisse
             IF(dzs(ks_split_min,i,j) > dzsmax(i,j)+dzsmax(i,j)*0.1_rsh) THEN
                DO k=ksmax+1,ks_split_min+2,-1
                   dzs(k,i,j)=dzs(k-1,i,j)
                   cv_sed(:,k,i,j)=cv_sed(:,k-1,i,j)
                   c_sedtot(k,i,j)=c_sedtot(k-1,i,j)
                   poro(k,i,j)=poro(k-1,i,j)
                ENDDO
                ksmax=ksmax+1
                dzs(ks_split_min+1,i,j)=dzs(ks_split_min,i,j)*0.5_rsh
                cv_sed(:,ks_split_min+1,i,j)=cv_sed(:,ks_split_min,i,j)
                c_sedtot(ks_split_min+1,i,j)=c_sedtot(ks_split_min,i,j)
                poro(ks_split_min+1,i,j)=poro(ks_split_min,i,j)
                dzs(ks_split_min,i,j)=dzs(ks_split_min,i,j)*0.5_rsh                
#ifdef key_MUSTANG_V2
                crel_mud(ks_split_min+1,i,j)=crel_mud(ks_split_min,i,j)
                poro_mud(ks_split_min+1,i,j)=poro_mud(ks_split_min,i,j)
#endif
             ENDIF
             
             ksma(i,j)=ksmax
            
#ifdef key_test_conservativity_splitsurf
             ! Verification of the conservativity of thicknesses and concentrations
             thick_surf_new=0.0_rsh
            ! stcvsed_surf_new=0.0_rsh
             stcsedtot_new=0.0_rsh
             DO k= MAX(1,ks_split_min-2), ksma(i,j)
               thick_surf_new=thick_surf_new+dzs(k,i,j)
               !stcvsed_surf_new=stcvsed_surf_new+cv_sed(3,k,i,j)*dzs(k,i,j)
               stcsedtot_new=stcsedtot_new+c_sedtot(k,i,j)*dzs(k,i,j)
             ENDDO          
             IF((thick_surf_new-thick_surf_old)/thick_surf_old*100 > 1.e-3 ) THEN
                write(*,*)'Probleme thick ',t,i,j,ksma(i,j),thick_surf_new,thick_surf_old, &
                          'diff %=',(thick_surf_new-thick_surf_old)/thick_surf_old*100
              ! write(*,*)'dzs old',j,ks_split_min,ks_surf,dzs_old(ks_split_min-2:ks_surf)
               write(*,*)'dzs new',j,dzs(ks_split_min-2:ksma(i,j),i,j)
             ENDIF
             !IF((stcvsed_surf_new-stcvsed_surf_old)/stcvsed_surf_old*100 > 1.e-3) THEN
             !   write(*,*)'Probleme cvsed ',i,j,ksma(i,j),stcvsed_surf_new,stcvsed_surf_old, &
             !              'diff %=',(stcvsed_surf_new-stcvsed_surf_old)/stcvsed_surf_old*100
             !ENDIF
             IF((stcsedtot_new-stcsedtot_old)/stcsedtot_old*100 > 1.e-3) THEN
                write(*,*)'Probleme csedtot ',i,j,ksma(i,j),stcsedtot_new,stcsedtot_old, &
                            'diff %=',(stcsedtot_new-stcsedtot_old)/stcsedtot_old*100
             ENDIF
#endif 
          
           ENDIF   ! if nlayer_splitt > 0                
         ENDIF   ! if ks_surf > 0 

  END SUBROUTINE sed_MUSTANG_split_surflayer
#endif
!!===========================================================================================
!
      SUBROUTINE sed_MUSTANG_coefbioturb_part(i,j,difbio)
!
!  !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_coefbioturb_part  ***
   !&E
   !&E ** Purpose :     Bioturbation coefficient in sediment (particulate )
   !&E
   !&E ** Description : 1DV 
   !&E  arguments IN : i,j         
   !&E  arguments OUT: difbio bioturbation coefficient (treated as a diffusion coefficient)
   !&E 
   !&E ** Called by :  sed_MUSTANG_consol_diff_bioturb
   !&E
! **********************************************************************
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used

   !! * Arguments
   INTEGER, INTENT(IN)            :: i,j   
   REAL(KIND=rsh),DIMENSION(ksdmin:ksdmax,-1:nv_adv),INTENT(OUT)  ::difbio                   
   
   !! * Local declarations
   INTEGER        :: iv,k,ksmib,ksmin,ksmax               
   REAL(KIND=rsh) :: dsurf,somcgravsand,somcmud,coefdifbio,frmud,hsedloc

   !!---------------------------------------------------------------------------
   !! * Executable part

        ksmax=ksma(i,j)
        ksmin=ksmi(i,j)
        difbio(ksdmin:ksdmax,1:nvp)=0.0_rsh
        IF (ksmax.GE.ksmin)THEN
          hsedloc=0.0_rsh
          DO k=ksmin,ksmax
            hsedloc=hsedloc+dzs(k,i,j)
            difbio(k,1:nvp)=0.0_rsh
          ENDDO   
          IF(hsedloc < dbiotu0_part) THEN
            ksmib=ksmax
          ELSE
            dsurf=-dzs(ksmax,i,j) ! evaluation of Db (k) at the upper limit
            DO k=ksmax,ksmin,-1
              dsurf=dsurf+dzs(k,i,j)
              IF(dsurf .GE. dbiotu0_part) THEN
                ksmib=k
                GOTO 11
              ENDIF
              coefdifbio=1.0_rsh
              IF(dsurf .GE. dbiotum_part) THEN
               coefdifbio=coefdifbio*EXP(-xbioturbk_part*(dsurf-dbiotum_part)/(dbiotu0_part-dbiotum_part))
              ENDIF
              difbio(k,isand1:nvp)=coefdifbio*xbioturbmax_part
#if ! defined key_noTSdiss_insed
              difbio(k,-1)=coefdifbio*xbioturbmax_part
#endif
            ENDDO    
            ksmib=1
          ENDIF
 11       DO k=ksmin,ksmib
            difbio(k,1:nvp)=0.0_rsh
#if ! defined key_noTSdiss_insed
            difbio(k,-1)=0.0_rsh
#endif

          ENDDO
        ENDIF
    
  END SUBROUTINE sed_MUSTANG_coefbioturb_part
!!===========================================================================================
!
#if ! defined key_noTSdiss_insed
  SUBROUTINE sed_MUSTANG_coefbioturb_diss(i,j,difbio)
!
!  !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_coefbioturb_diss  ***
   !&E
   !&E ** Purpose :     Bioturbation coefficient in sediment dissolved
   !&E
   !&E ** Description : 1DV 
   !&E  arguments IN : i,j                   
   !&E  arguments OUT: difbio bioturbation coefficient (treated as a diffusion coefficient)
   !&E              
   !&E ** Called by :  sed_MUSTANG_consol_diff_bioturb
   !&E
! **********************************************************************
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used

   !! * Arguments
   INTEGER, INTENT(IN)            :: i,j   
   REAL(KIND=rsh),DIMENSION(ksdmin:ksdmax,-1:nv_adv),INTENT(OUT)  ::difbio                   
   
   !! * Local declarations
   INTEGER        :: iv,k,ksmib,ksmin,ksmax               
   REAL(KIND=rsh) :: dsurf,somcgravsand,somcmud,coefdifbio,frmud,hsedloc

   !!---------------------------------------------------------------------------
   !! * Executable part

        ksmax=ksma(i,j)
        ksmin=ksmi(i,j)
        difbio(ksdmin:ksdmax,nvp+1:nv_adv)=0.0_rsh
        difbio(ksdmin:ksdmax,0)=0.0_rsh

        IF (ksmax.GE.ksmin)THEN
          hsedloc=0.0_rsh
          DO k=ksmin,ksmax
            hsedloc=hsedloc+dzs(k,i,j)
            difbio(k,nvp+1:nv_adv)=0.0_rsh
          ENDDO   
          IF(hsedloc < dbiotu0_diss) THEN
            ksmib=ksmax
          ELSE
            dsurf=-dzs(ksmax,i,j) ! evaluation of Db (k) at the upper limit
            DO k=ksmax,ksmin,-1
              dsurf=dsurf+dzs(k,i,j)
              IF(dsurf .GE. dbiotu0_diss) THEN
                ksmib=k
                GOTO 11
              ENDIF
              !difbio dependent on mud fraction ??
              somcgravsand=0.0_rsh
              DO iv=igrav1,isand2
                somcgravsand=somcgravsand+cv_sed(iv,k,i,j)
              ENDDO
              somcmud=0.0_rsh
              DO iv=imud1,imud2
                somcmud=somcmud+cv_sed(iv,k,i,j)
              ENDDO
              frmud=somcmud/(somcgravsand+somcmud)
              IF(frmud < frmud_db_min) THEN
                coefdifbio=0.0_rsh
              ELSE IF(frmud > frmud_db_max) THEN
                coefdifbio=1.0_rsh
              ELSE
                coefdifbio=(frmud-frmud_db_min)/(frmud_db_max-frmud_db_min)
              ENDIF
              IF(dsurf .GE. dbiotum_diss) THEN
               ! difbio=0 for gravels

               coefdifbio=coefdifbio*EXP(-xbioturbk_diss*(dsurf-dbiotum_diss)/(dbiotu0_diss-dbiotum_diss))
              ENDIF
               ! difbio same for sand, muds and dissolved variables here - to be changed by user
                difbio(k,nvp+1:nv_adv)=coefdifbio*xbioturbmax_diss
                difbio(k,0)=coefdifbio*xbioturbmax_diss
            ENDDO    
            ksmib=1
          ENDIF
 11       DO k=ksmin,ksmib
            difbio(k,:)=0.0_rsh
            difbio(k,0)=0.0_rsh
          ENDDO
        ENDIF
 
    
  END SUBROUTINE sed_MUSTANG_coefbioturb_diss
#endif   
!!==============================================================================
  SUBROUTINE sed_MUSTANG_constitutivrel(i,j,stateconsol,permeab,sigmapsg)
! 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_constitutivrel  ***
   !&E
   !&E ** Purpose : applies constitutive relationships for permeability, 
   !&E              effective stess
   !&E
   !&E ** Description :
   !&E
   !&E ** Note : GRAVITY must be known as a parameters transmtted by coupleur 
   !&E           in MARS : coupleur_dimhydro.h (USE ..)
   !&E           in CROCO : module_MUSTANG.F (include..)
   !&E
   !&E ** Called by :  sed_MUSTANG_consol_diff_bioturb
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used

   !! * Arguments
   INTEGER,INTENT(IN) :: i,j
   REAL(KIND=rsh),DIMENSION(ksdmin:ksdmax),INTENT(out) :: stateconsol,permeab,sigmapsg
   
  !! * Local declarations
   INTEGER                        :: k,iv,icond
   REAL(KIND=rsh),DIMENSION(nvpc) :: cvol
   REAL(KIND=rsh)                 :: cvolsed,cvolsangrv,fracvolcum,diarepres,       &
                                     voidratio,diamstarrep,wsoulsby,somcmud,crelmud,cvolrelmud, &
                                     permeab_e


   !! * executable statement
   !!   --------------------
   sigmapsg(:)=0.0_rsh
   stateconsol(:)=0.0_rsh

   DO k=ksmi(i,j),ksma(i,j)
     cvol(:)=0.0_rsh
     cvolsed=0.0_rsh
     DO iv=igrav1,isand2
       cvol(iv)=cv_sed(iv,k,i,j)/ros(iv)
       cvolsed=cvolsed+cvol(iv)
     ENDDO
     cvolsangrv=cvolsed
     somcmud=0.0_rsh
     DO iv=isand2+1,nvpc
       cvol(iv)=cv_sed(iv,k,i,j)/ros(iv)
       cvolsed=cvolsed+cvol(iv)
       somcmud=somcmud+cv_sed(iv,k,i,j)
     ENDDO
       
     !  permeability with voidratio:
     voidratio=poro(k,i,j)/(1.0000001_rsh-poro(k,i,j))
     voidratio=MAX(voidratio,0.0_rsh)
     permeab_e=2.0e-9_rsh*voidratio**3.7_rsh   ! best param for Co<120 kg/m3 (Grasso et al., Ocean Dynamics 2014)

     !  permeability with cvolrelmud (relative volumic concentration of mud, cvolrelmud=cvolmud/(1-cvolsand), cf. Merckelbach & Kranenburg, 2004)
     crelmud=somcmud/(1.0_rsh-cvolsangrv)
     cvolrelmud=crelmud/ros(1)
     !!!FG(06/03/2015) to correct the case where cvolrelmud is too small (Cvolrelmud = 0.01 -> Crelmud = 26 g/l -> permeab = 4 m/s)
     IF(cvolrelmud .LE. 0.01) THEN
       permeab(k)=permeab_e
     ELSE
       permeab(k)=xperm1*cvolrelmud**xperm2
     ENDIF

!     !  minimalisation of permeability
     permeab(k)=MIN(permeab(k),permeab_e,(1.0_rsh+voidratio)*0.0001_rsh/1.65_rsh)   ! equivalent permeability for ws=0.1 mm/s

     !  effective stress :
     icond=1
     DO iv=igrav1,isand2
       icond=icond*MAX(0.0_rsh,100.0_rsh*(cvolmaxsort-cvol(iv)-.02_rsh))
     ENDDO
     IF(cvolsed.LE.cvolmaxmel.AND.icond.NE.0)THEN
       sigmapsg(k)=xsigma1/GRAVITY*cvolrelmud**xsigma2 !FG(04/09/2013) to use Merckelback & Kranenburg s (2004) formulation
       stateconsol(k)=1.0_rsh
     ELSE
       !    in this case, the effective stress is maximum, 
       !    but porewater overpressure is zero :
       stateconsol(k)=0.0_rsh
     ENDIF
   ENDDO

  END SUBROUTINE sed_MUSTANG_constitutivrel

 !!==============================================================================
#if defined key_MUSTANG_V2
 !!==============================================================================

   SUBROUTINE MUSTANGV2_manage_active_layer(i,j,ksmax &
#if ! defined key_noTSdiss_insed 
                        ,flx_s2w_eroij                                &
#endif
#if ! defined key_nofluxwat_IWS
                        ,phieau_ero_ij                                &
#endif
                        )

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE MUSTANGV2_manage_active_layer  ***
   !&E
   !&E ** Purpose : compute the active layer thickness and update its composition
   !&E
   !&E ** Description : 0D
   !&E                  variables OUT : ksmax,dzs_activelayer_comp, dzs, poro
   !&E                                  c_sedtot,poro_mud,crel_mud,cv_sed
   !&E                  variables IN : ksmi,ksmax,diam_sed, cv_sed,stresscri0
   !&E                         k1HW97,tauskin(i,j),k2HW97,coef_frmudcr1,frmudcr2
   !&E                         isitcohesive,dzs
   !&E                         c_sedtot,poro,poro_mud,crel_mud,
   !&E 
   !&E ** Called by :  sed_MUSTANG__erosion
   !&E
   !&E--------------------------------------------------------------------------

   !! * Modules used

   !! * Arguments
   INTEGER, INTENT(IN) :: i,j
   INTEGER, INTENT(INOUT) :: ksmax

#if ! defined key_noTSdiss_insed 
   REAL(KIND=rsh),DIMENSION(-1:nv_adv),INTENT(INOUT)  ::  flx_s2w_eroij
#endif
#if ! defined key_nofluxwat_IWS
   REAL(KIND=rsh),INTENT(INOUT) ::   phieau_ero_ij
#endif


   !! * Local declaration
   INTEGER :: iv,ksmin,ksmaxa,nb_iter_aclay
   REAL(KIND=rsh) :: diamsan,somsan,frmudcr_fusion,&
                     dzs_activelayer_comp,dzs_excess,dzsa,&
                     dzsr,c_sedtotr,poror,poro_mudr,crel_mudr,&
                     somgravsan,diamgravsan,critstresgravsan
   REAL(KIND=rsh),DIMENSION(nvpc):: frac_sed
   REAL(KIND=rsh),DIMENSION(-1:nv_tot) :: cv_sedr
   LOGICAL :: l_ksmaxm1_cohesive,l_stop_fusion_in_activelayer


   !!---------------------------------------------------------------------------
   !! * Executable part

   !print *,''
   !print *,'************************************'
   !print *,'   > Enter MUSTANGV2_manage_active_layer'
   !print *,''

   ksmin=ksmi(i,j)

   ! A boolean l_stop_fusion_in_activelayer (initially at false) is introduced to manage the fusion 
   ! within the active layer and to prevent the final number of iteration (potentially time consuming)

   l_stop_fusion_in_activelayer=.FALSE.
   nb_iter_aclay=0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!! INITIALISATION !!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! First, this boolean needs to be initialiased according to different criteria which are:
   ! (1) the ksmax layer thickness (initial or post fusion) does not exceed the theorical one 
   !                                       computed from Harris and Wiberg (1997) formulation
   ! (2) the ksmax-1 layer is not cohesive
   ! (3) ksmax-1 > ksmin (because the ksmin layer is never eroded)

   ! If one of the three criteria is not respected, l_stop_fusion_in_activelayer=.TRUE. and 
   !                   the ksmax layer will be our active layer and no fusion process occurs

     ! --> Criterion (1)

   somgravsan=0.0_rsh
   diamgravsan=0.0_rsh
   critstresgravsan=0.0_rsh
   DO iv=igrav1,isand2
     somgravsan=somgravsan+cv_sed(iv,ksmax,i,j)
     diamgravsan=diamgravsan+diam_sed(iv)*cv_sed(iv,ksmax,i,j)
     critstresgravsan=critstresgravsan+stresscri0(iv)*cv_sed(iv,ksmax,i,j)
   END DO

   diamgravsan=MAX(diamgravsan/(somgravsan+epsilon_MUSTANG),diam_sed(isand2))
   critstresgravsan=MAX(critstresgravsan/(somgravsan+epsilon_MUSTANG),stresscri0(isand2))

   ! With this formulation, the active layer thickness cannot be lower than 6 times the diameter 
   !       of the finest sand (0.6 mm in the SEINE configuration)
         ! tauskin in dy/cm2 and diam in cm in this formula
   dzs_activelayer_comp = (k1HW97 * 10.0_rsh * MAX(tauskin(i,j) - critstresgravsan, 0.0_rsh)) + (k2HW97 * diamgravsan * 100.0_rsh) 
   dzs_activelayer_comp = dzs_activelayer_comp / 100.0_rsh ! .. now dzs in m

   IF (dzs(ksmax,i,j) .GE. dzs_activelayer_comp) l_stop_fusion_in_activelayer=.TRUE.

     ! --> Criterion (2) : ne devrait on pas prendre en compte le gravier ?? 
   ! FD+PLH 202112 Attention: critère frmudcr_fusion calculé à partir de la couche active (i.e. ksmax)
   !   
   frmudcr_fusion=MIN(fusion_para_activlayer*coef_frmudcr1*diamgravsan,frmudcr2)
   l_ksmaxm1_cohesive=isitcohesive(cv_sed(:,ksmax-1,i,j),frmudcr_fusion) 
            ! return l_ksmaxm1_cohesive=.TRUE. if frmud(ksmax-1) >= frmudcr_fusion
   IF (l_ksmaxm1_cohesive) l_stop_fusion_in_activelayer=.TRUE.

     ! --> Criterion 3

   IF (ksmax-1 .EQ. ksmin) l_stop_fusion_in_activelayer=.TRUE.


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!! FUSION IN ACTIVE LAYER !!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! While l_stop_fusion_in_activelayer=.FALSE. according to the 3 criteria, fusion in surficial sediment continues

   ! Updates at each iteration : 
   !   > the ksmax, composition, and thickness of the fusioned layer by accounting of the new porosity computation
   !   > Active layer thickness from the new ksmax layer composition 
   !   > Boolean l_stop_fusion_in_activelayer

   DO WHILE (.NOT. l_stop_fusion_in_activelayer)

     nb_iter_aclay=nb_iter_aclay+1

     dzs_excess=dzs(ksmax,i,j)+dzs(ksmax-1,i,j)-dzs_activelayer_comp
     frac_sed(1:nvpc)=cv_sed(1:nvpc,ksmax-1,i,j)/c_sedtot(ksmax-1,i,j)
     !dzsmin=(1.0_rsh-coeff_dzsmin)*dzsminuni + &
     !        coeff_dzsmin*SUM( (cv_sed(1:nvpc,ksmax-1,i,j)/c_sedtot(ksmax-1,i,j))*diam_sed(1:nvpc) )
      dzsmin=dzsminvar(frac_sed)

     ! In this case, the fusion process is stopped to prevent the final number of iteration (potentially time consuming)
     IF (dzs(ksmax,i,j)+dzs(ksmax-1,i,j) .GT. dzs_activelayer_comp) THEN
         l_stop_fusion_in_activelayer=.TRUE.

        IF (dzs_excess .GT. dzsmin) THEN

          ! We keep in memory the excess of the layer ksmax-1 that will be added after fusion below the new layer ksmax
          dzsr=dzs_excess
          cv_sedr(:)=cv_sed(:,ksmax-1,i,j)
          c_sedtotr=c_sedtot(ksmax-1,i,j)
          poror=poro(ksmax-1,i,j)
          poro_mudr=poro_mud(ksmax-1,i,j)
          crel_mudr=crel_mud(ksmax-1,i,j)

          ! The thickness of the ksmax-1 layer is reduced to that necessary to satisfy the active layer criterion
          ! (concentrations, porosite, crelmud .. remain unchanged)
          dzs(ksmax-1,i,j)=dzs_activelayer_comp-dzs(ksmax,i,j)

          ! Ksmax layer is fused with this part of ksmax-1 to satisfy the active layer criterion
           !  computing new cv_sed, including  dissolved variables conc, and IWS fluxes and water fluxes
          CALL MUSTANGV2_fusion_with_poro(i,j,ksmax   &
#if ! defined key_noTSdiss_insed 
                        ,flx_s2w_eroij                                &
#endif
#if ! defined key_nofluxwat_IWS
                        ,phieau_ero_ij                                &
#endif
                        )
     
          ! We add dzs_excess corresponding to the excedent in the initial ksmax-1 layer below the ksmax layer
          ksmaxa=ksmax
          ksmax=ksmax+1
    
          ! the layer ksmax is set to +1
          dzs(ksmax,i,j)=dzs(ksmaxa,i,j)
          cv_sed(:,ksmax,i,j)=cv_sed(:,ksmaxa,i,j)
          c_sedtot(ksmax,i,j)=c_sedtot(ksmaxa,i,j)
          poro(ksmax,i,j)=poro(ksmaxa,i,j)
          poro_mud(ksmax,i,j)=poro_mud(ksmaxa,i,j)
          crel_mud(ksmax,i,j)=crel_mud(ksmaxa,i,j)
    
          ! We add the excedent in the new layer at ksmax-1  
          dzs(ksmaxa,i,j)=dzsr
          cv_sed(:,ksmaxa,i,j)=cv_sedr(:)
          c_sedtot(ksmaxa,i,j)=c_sedtotr
          poro(ksmaxa,i,j)=poror
          poro_mud(ksmaxa,i,j)=poro_mudr
          crel_mud(ksmaxa,i,j)=crel_mudr

        ELSE

           !!! Fusion of the two upper layers
            !  computing new cv_sed, including  dissolved variables conc, and IWS fluxes and water fluxes
         CALL MUSTANGV2_fusion_with_poro(i,j,ksmax  &
#if ! defined key_noTSdiss_insed 
                        ,flx_s2w_eroij                                &
#endif
#if ! defined key_nofluxwat_IWS
                        ,phieau_ero_ij                                &
#endif
                        )

        END IF
     ELSE

       !!! Fusion of the two upper layers
           !  computing new cv_sed, including  dissolved variables conc, and IWS fluxes and water fluxes
       CALL MUSTANGV2_fusion_with_poro(i,j,ksmax   &
#if ! defined key_noTSdiss_insed 
                        ,flx_s2w_eroij                                &
#endif
#if ! defined key_nofluxwat_IWS
                        ,phieau_ero_ij                                &
#endif
                         )

     END IF

     !!! Update of l_stop_fusion_in_activelayer according to the 3 criteria (see above initialisation section)
       ! If one of the three criteria is not respected, l_stop_fusion_in_activelayer=.TRUE. and 
       !                   the ksmax layer will be our active layer and no fusion process occurs

     ! --> Criterion (1)

     somgravsan=0.0_rsh
     diamgravsan=0.0_rsh
     critstresgravsan=0.0_rsh
     DO iv=igrav1,isand2
       somgravsan=somgravsan+cv_sed(iv,ksmax,i,j)
       diamgravsan=diamgravsan+diam_sed(iv)*cv_sed(iv,ksmax,i,j)
       critstresgravsan=critstresgravsan+stresscri0(iv)*cv_sed(iv,ksmax,i,j)
     END DO

     diamgravsan=MAX(diamgravsan/(somgravsan+epsilon_MUSTANG),diam_sed(isand2))
     critstresgravsan=MAX(critstresgravsan/(somgravsan+epsilon_MUSTANG),stresscri0(isand2))

     ! With this formulation, the active layer thickness cannot be lower than 6 times the diameter of the finnest sand 
     !                          (0.6 mm in the SEINE configuration)
     dzs_activelayer_comp = (k1HW97*10.0_rsh*MAX(tauskin(i,j)-critstresgravsan,0.0_rsh)) + (k2HW97*diamgravsan*100.0_rsh) 
                     ! tauskin in dy/cm2 and diam in cm
     dzs_activelayer_comp = dzs_activelayer_comp/100.0_rsh ! .. now dzs in m

     IF (dzs(ksmax,i,j) .GE. dzs_activelayer_comp) l_stop_fusion_in_activelayer=.TRUE.

     ! --> Criterion (2)

     frmudcr_fusion=MIN(fusion_para_activlayer*coef_frmudcr1*diamgravsan,frmudcr2)
     l_ksmaxm1_cohesive=isitcohesive(cv_sed(:,ksmax-1,i,j),frmudcr_fusion)
               ! return l_ksmaxm1_cohesive=.TRUE. if frmud(ksmax-1) >= frmudcr_fusion
     IF (l_ksmaxm1_cohesive) l_stop_fusion_in_activelayer=.TRUE.

     ! --> Criterion 3

     IF (ksmax-1 .EQ. ksmin) l_stop_fusion_in_activelayer=.TRUE.


   END DO ! while l_stop_fusion_in_activelayer==.FALSE.

#if defined key_MUSTANG_specif_outputs 
   varspecif2D_save(3,i,j)=dzs_activelayer_comp ! dzs_aclay_comp_save
   varspecif2D_save(4,i,j)=dzs(ksmax,i,j)       ! dzs_aclay_kept_save
#endif

   !print *,'Exit MUSTANGV2_manage_active_layer'
   !print *,'************************************'
   !print *,''


END SUBROUTINE MUSTANGV2_manage_active_layer
  !!==============================================================================
SUBROUTINE MUSTANGV2_fusion_with_poro(i,j,ksmax   &
#if ! defined key_noTSdiss_insed 
                        ,flx_s2w_eroij                                &
#endif
#if ! defined key_nofluxwat_IWS
                        ,phieau_ero_ij                                &
#endif
                        )

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE MUSTANGV2_fusion_with_poro  ***
   !&E
   !&E ** Purpose : fusion of two surficial sediment layers accounting for new 
   !&E              computation of porosity
   !&E
   !&E ** Description : 
   !&E                  variables OUT : ksmax, dzs, poro,poro_mud,crel_mud, cv_sed
   !&E                  variables IN : ksmi,ksmax,dzs, poro,poro_mud,crel_mud,
   !&E
   !&E ** Called by :  MUSTANGV2_manage_active_layer
   !&E
   !&E--------------------------------------------------------------------------

   !! * Modules used

   !! * Arguments
   INTEGER, INTENT(IN)        :: i,j
   INTEGER, INTENT(INOUT)     :: ksmax

#if ! defined key_noTSdiss_insed 
   REAL(KIND=rsh),DIMENSION(-1:nv_adv),INTENT(INOUT)  ::  flx_s2w_eroij
#endif
#if ! defined key_nofluxwat_IWS
   REAL(KIND=rsh),INTENT(INOUT) ::   phieau_ero_ij
#endif


   !! * Local declaration
   INTEGER :: iv,ksmaxa,ivp_assoc
   REAL(KIND=rsh) :: sommud_k,sommud_km1,cvolinigrv_k,cvolinigrv_km1,cvolinisan_k,cvolinisan_km1,&
                     cmudr_k,cmudr_km1,poro_mud_k,poro_mud_km1,dzsa,dzsam1,&
                     massmud_k,massmud_km1,massmud_tot,poro_kij,&
                     poroa,poroam1,poro_mud_new,mass_tot,&
                     dzsi,crel_mud_new

   REAL(KIND=rsh),DIMENSION(nvp) :: mass_sed
   REAL(KIND=rsh),DIMENSION(nvpc):: frac_sed
#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
   REAL(KIND=rsh) :: porowater1,porowater2, porowater_new
#endif

   !!---------------------------------------------------------------------------
   !! * Executable part

   !print *,''
   !print *,'************************************'
   !print *,'   > Enter MUSTANGV2_fusion_with_poro'
   !print *,'************************************'
   !print *,''
   
     ! 1- We look at the composition and the porosity of the mud in layers ksmax and ksmax-1
     ! 2- A porosity of the representative mud is calculated (pro rata of the masses)
     ! 3- The sediment masses / fractions of the fused layer are calculated
     ! 4- We calculate the porosity -> ksmax = ksmax-1 -> MAJ dzs, cv_sed ksmax
          
     ! 1- On regarde la composition et la porosite de la vase des couches ksmax et ksmax-1
     ! 2- On calcule une porosite de la vase representative (au prorata des masses) 
     ! 3- On calcule les masses/fractions sedimentaires de la couche fusionnee
     ! 4- On calcule la porosite --> ksmax=ksmax-1 --> MAJ dzs, cv_sed ksmax

     !1-
     sommud_k=0.0_rsh
     sommud_km1=0.0_rsh
     DO iv=imud1,imud2
       sommud_k=sommud_k+cv_sed(iv,ksmax,i,j)
       sommud_km1=sommud_km1+cv_sed(iv,ksmax-1,i,j)
     ENDDO

     dzsa=dzs(ksmax,i,j)
     dzsam1=dzs(ksmax-1,i,j)

     massmud_k=sommud_k*dzsa
     massmud_km1=sommud_km1*dzsam1
     massmud_tot=massmud_k+massmud_km1

     !poro_mud_k=1.0_rsh-(sommud_k/ros(1))
     !poro_mud_km1=1.0_rsh-(sommud_km1/ros(1))
     
     poro_mud_k=poro_mud(ksmax,i,j)
     poro_mud_km1=poro_mud(ksmax-1,i,j)

     poroa=poro(ksmax,i,j)
     poroam1=poro(ksmax-1,i,j)

#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
     ! for  dissolved concentrations et fluxes estimation
     porowater1=poroa*dzsa
     porowater2=poroam1*dzsam1
#endif

     !print *,'**Couche ksmax**'
     !print *,'ksmax=',ksmax
     !print *,'dzs(ksmax,i,j)=',dzs(ksmax,i,j)
     !print *,'cv_sed(:,ksmax,i,j)=',cv_sed(:,ksmax,i,j)
     !print *,'poro_mud en ksmax=',poro_mud_k
     !print *,'poro(ksmax,i,j)=',poro(ksmax,i,j)
     !print *,''
     !print *,'**Couche ksmax-1**'
     !print *,'ksmax-1=',ksmax-1
     !print *,'dzs(ksmax-1,i,j)=',dzs(ksmax-1,i,j)
     !print *,'cv_sed(:,ksmax-1,i,j)=',cv_sed(:,ksmax-1,i,j)
     !print *,'poro_mud en ksmax-1=',poro_mud_km1
     !print *,'poro(ksmax-1,i,j)=',poro(ksmax-1,i,j)
     !print *,''

     !2-
     IF (massmud_tot .GT. 0.0_rsh) THEN
       poro_mud_new=(massmud_k/massmud_tot)*poro_mud_k + (massmud_km1/massmud_tot)*poro_mud_km1
       crel_mud_new=(massmud_k/massmud_tot)*crel_mud(ksmax,i,j) + (massmud_km1/massmud_tot)*crel_mud(ksmax-1,i,j)
     ELSE
       poro_mud_new=0.0_rsh
       crel_mud_new=0.0_rsh
     END IF

     !3-
     mass_tot = 0.0_rsh
     DO iv=igrav1,imud2
       mass_sed(iv)=cv_sed(iv,ksmax,i,j)*dzsa + cv_sed(iv,ksmax-1,i,j)*dzsam1
       mass_tot=mass_tot+mass_sed(iv)
     END DO

     DO iv=igrav1,imud2
       frac_sed(iv)=mass_sed(iv)/mass_tot
       !frac_sed(iv)=MAX(mass_sed(iv),0.00000000001_rlg)/mass_tot ! Vraiment utile ?! a priori non suite aux dernieres modifs ==> je commente
     END DO

     !4-
     CALL MUSTANGV2_comp_poro_mixsed(frac_sed, poro_mud_new,   &
                                 crel_mud_new, poro_kij)

     ! Question : faut il prescrire a la toute fin que dzs et cv_sed en ksmax avant fusion
     ! sont desormais = 0.0 (ou -999.0)

     ksmaxa=ksmax

     ksmax=ksmax-1

     dzs(ksmax,i,j)=mass_tot/((1.0_rsh-poro_kij)*ros(1))
     dzsi=1.0_rsh/dzs(ksmax,i,j)

     crel_mud(ksmax,i,j)=crel_mud_new
     poro(ksmax,i,j)=poro_kij     
     poro_mud(ksmax,i,j)=poro_mud_new ! Actualisation de poro_mud liee au melange !

     DO iv=1,nvpc
       cv_sed(iv,ksmax,i,j)=mass_sed(iv)*dzsi
     END DO

     !IF (crel_mud(ksmax,i,j) .GT. 1500.0_rsh) THEN
     !  print *,'in sed_fusion_with_poro'
     !  print *,' > crel_mud(ksmax,i,j) = ',crel_mud(ksmax,i,j)
     !  print *,' > cv_sed(:,ksmax,i,j) = ',cv_sed(:,ksmax,i,j)
     !END IF


      ! non constitutive particulate variables
     DO iv=nvpc+1,nvp
#ifdef key_Pconstitonly_insed
       cv_sed(iv,ksmax,i,j)=0.0_rsh
#else
       cv_sed(iv,ksmax,i,j)=(cv_sed(iv,ksmaxa,i,j)*dzsa+cv_sed(iv,ksmax,i,j)*dzsam1)*dzsi
#endif
     ENDDO

#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
     ! dissolved variable in pore waters and water fluxes at the interface
     ! --------------------------------------------------------------------
     !  ==> flx_s2w and cv_sed and phieau_s2w
     ! fusion in active layer ==> mixing with water conscentrations ==> cv_sed=cw_bottom 
     ! code=4
       porowater_new=poro(ksmax,i,j)*dzs(ksmax,i,j)
       CALL MUSTANGV2_eval_dissvar_IWSflux(i,j,ksmax,4,    &
           phieau_ero_ij=phieau_ero_ij,          &
           flx_s2w_eroij=flx_s2w_eroij,          &
           porowater1=porowater1, &
           cv_sed1=cv_sed(-1:nv_adv,ksmaxa,i,j),     &
           porowater2=porowater2, &
           cv_sed2=cv_sed(-1:nv_adv,ksmax,i,j),       &
           porowater_new=porowater_new)
#endif

     c_sedtot(ksmax,i,j)=0.0_rsh
     DO iv=1,nv_use ! nv_use=nvp (or nvpc if key_Pconstitonly_insed)
       c_sedtot(ksmax,i,j)=c_sedtot(ksmax,i,j)+cv_sed(iv,ksmax,i,j)*typart(iv)
     END DO

     !print *,'Apres fusion'
     !print *,'ksmax=',ksmax
     !print *,'dzs(ksmax,i,j)=',dzs(ksmax,i,j)
     !print *,'cv_sed(:,ksmax,i,j)=',cv_sed(:,ksmax,i,j)
     !print *,'poro_mud_new en ksmax=',poro_mud_new
     !print *,'poro(ksmax,i,j)=',poro(ksmax,i,j)
     !print *,''

   !print *,'Exit sed_fusion_with_poro'
   !print *,'************************************'
   !print *,''


END SUBROUTINE MUSTANGV2_fusion_with_poro
  !!==============================================================================
  
  SUBROUTINE MUSTANGV2_comp_eros_flx_indep(i, j, ksmax,                            &
#ifdef key_MUSTANG_bedload
                                      CELL_DX, CELL_DY, flx_bxij, flx_byij,         &
                                      BAROTROP_VELOCITY_U, BAROTROP_VELOCITY_V,   &
#endif
                                      sed_eros_flx_class_by_class)

   !&E--------------------------------------------------------------------------
   !&E                 *** MUSTANGV2_comp_eros_flx_indep  ***
   !&E
   !&E ** Purpose : Non-cohesive class independent erosion
   !&E
   !&E ** Description : 0D
   !&E         variables IN : BAROTROP_VELOCITY_U,BAROTROP_VELOCITY_V,
   !&E                        E0_sand,ksmax,cv_sed,c_sedtot,diam_sed
   !&E                        l_peph_suspension, stresscri0, tauskin
   !&E                        l_fsusp,ws_sand, n_eros_sand,n_eros_mud
   !&E                        crel_mud, tau_cri_mud_option_eroindep
   !&E                        x1toce_mud,x2toce_mud,E0_mud_para_indep,E0_mud
   !&E
   !&E         variables OUT : sed_eros_flx_class_by_class
   !&E
   !&E ** Called by :  sed_erosion
   !&E
   !&E--------------------------------------------------------------------------

   !! * Modules used

   !! * Arguments
   INTEGER, INTENT(IN) :: i, j, ksmax
   REAL(KIND=rsh),DIMENSION(1:nvp),INTENT(OUT) :: sed_eros_flx_class_by_class 
#ifdef key_MUSTANG_bedload
   REAL(KIND=rsh),DIMENSION(1:nvp),INTENT(IN)              :: flx_bxij 
   REAL(KIND=rsh),DIMENSION(1:nvp),INTENT(IN)              :: flx_byij
   REAL(KIND=rsh),DIMENSION(ARRAY_CELL_DX),INTENT(IN)      :: CELL_DX
   REAL(KIND=rsh),DIMENSION(ARRAY_CELL_DY),INTENT(IN)      :: CELL_DY
   REAL(KIND=rsh),DIMENSION(ARRAY_VELOCITY_U),INTENT(IN)   :: BAROTROP_VELOCITY_U                       
   REAL(KIND=rsh),DIMENSION(ARRAY_VELOCITY_V),INTENT(IN)   :: BAROTROP_VELOCITY_V                         
#endif

   !! * Local declaration
   INTEGER                  :: iv, jiv
   REAL(KIND=rsh),PARAMETER :: m = 0.6_rsh ! for hinding / exposure processes
   REAL(KIND=rsh)           :: ph, pe, pephm_fcor, sed_eros_flxsand, frac_sand
   REAL(KIND=rsh)           :: c_sed_tot, sommud, cvolinigrv, cvolinisan, cmudr, tauc_mud, speed, fsusp 
   REAL(KIND=rsh)           :: somsan, toce_meansan, toce_minsan, E0_mud_loc
   REAL(KIND=rsh),DIMENSION(isand1:isand2) :: toce_sand
   REAL(KIND=rsh),DIMENSION(nvpc)          :: frac_sed, toce_loc
   REAL(KIND=rsh),DIMENSION(nvp)           :: E0_sand_loc

   !!---------------------------------------------------------------------------
   !! * Executable part

   !print *,''
   !print *,'************************************'
   !print *,'   > Enter sed_comp_eros_flx_indep'
   !print *,''

   ! Needs to be initialised here
   sed_eros_flx_class_by_class(:)=0.0_rsh

   ! Critical shear stress for each class based on masking / exposure processes
   toce_loc(:)=0.0_rsh

   E0_sand_loc(:)=E0_sand(:)

   DO iv=igrav1,imud2
     frac_sed(iv)=cv_sed(iv,ksmax,i,j)/c_sedtot(ksmax,i,j)
   END DO
 
#ifdef key_MUSTANG_debug
    IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
     print *,'    SED COMP EROS FLX INDEP'
     print *,''
     print *,'      > Ini sed_eros_flx_class_by_class(:)=', sed_eros_flx_class_by_class(:)
     print *,'      > E0_sand_loc(:)=', E0_sand_loc(:)
     print *,'      > frac_sed(:)=', frac_sed(:)
     print *,'      > tauskin(i,j)=', tauskin(i,j)
     print *,''
     print *,'      ** NON COHESIVE SEDIMENTS **'
     print *,''
   END IF
#endif
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!! SANDY SEDIMENTS        !!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   DO iv=isand1,isand2
#ifdef key_MUSTANG_debug
     IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
       print *,'      iv=',iv
     END IF
#endif

     !!! Critical shear stress toce in N/m2
     IF (l_peph_suspension) THEN
       ph = 0.0_rsh
       pe = 0.0_rsh
       ! Hindering / exposure coefficients for sediment class iv, ph and pe, respectively
       DO jiv=1,isand2 ! imud2 ou isand2 --> reflechir --> preference pour isand2
         ph = ph + (diam_sed(jiv)/(diam_sed(jiv)+diam_sed(iv)))*  &
                     ( max(0.0_rsh,cv_sed(jiv,ksmax,i,j))/(max(epsilon_MUSTANG,c_sedtot(ksmax,i,j))) )
         pe = pe + (diam_sed(iv) /(diam_sed(jiv)+diam_sed(iv)))*  &
                     ( max(0.0_rsh,cv_sed(jiv,ksmax,i,j))/(max(epsilon_MUSTANG,c_sedtot(ksmax,i,j))) )
       END DO
       pephm_fcor = (pe / ph)**(-m)
       toce_loc(iv) = stresscri0(iv) * pephm_fcor
#if defined key_MUSTANG_specif_outputs        
       varspecif3Dnv_save(4, iv, i, j)=pephm_fcor
#endif
     ELSE
       toce_loc(iv) = stresscri0(iv)
     END IF

#if defined key_MUSTANG_specif_outputs        
     varspecif3Dnv_save(1,iv,i,j) = toce_loc(iv)  ! toce_save
#endif

     IF (tauskin(i,j) .GT. toce_loc(iv)) THEN

#ifdef key_MUSTANG_bedload
       ! fsusp is the suspension part of the whole iv transport (suspension + bedload) according to Wu and Lin (2014) 
       ! It is applied to E0_sand parameter to prevent any overestimation of sediment transport in the event that
       !      bedload & suspension are accounted for
       ! Warning tauskin suspension = tauskin bedload while different in Wu and Lin

       IF (l_fsusp .and. iv.le.ibedload2) THEN
         speed = SQRT( CURRENTU_ij**2+ CURRENTV_ij**2) 
         fsusp= (0.0000262_rsh*((speed/ws_sand(iv))**1.74_rsh)) /  &
                ( (0.0000262_rsh*((speed/ws_sand(iv))**1.74_rsh))  &
                + (0.0053_rsh*(tauskin(i,j)/toce_loc(iv)-1.0_rsh)**0.46_rsh) )
         E0_sand_loc(iv)=fsusp*E0_sand_loc(iv)
#if defined key_MUSTANG_specif_outputs        
         varspecif3Dnv_save(8,iv,i,j)=fsusp ! check output
#endif
       END IF
#endif

       !!! Compute of erosion flux sed_eros_flx_class_by_class in kg/m2/s
       ! As already mentioned in paraMUSTANG, if E0_sand_option == 3, 
       !  n_eros_sand should be equal to 1.7 (to be consistent with Wu and Lin formulation)

       sed_eros_flx_class_by_class(iv)=MF*fwet(i,j)*frac_sed(iv)*E0_sand_loc(iv) &
                       *((tauskin(i,j)/toce_loc(iv))-1.0_rsh)**n_eros_sand

#ifdef key_MUSTANG_debug
       IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
         print *,'       EROSION of iv !'
         IF (l_peph_suspension) print *,'       - pe/ph=',pe,' / ',ph
         print *,'       - toce(iv)=',toce_loc(iv)
         print *,'       - sed_eros_flx_class_by_class(iv)=',sed_eros_flx_class_by_class(iv)
       END IF
#endif

     END IF

   END DO

      
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! MUDDY SEDIMENTS    !!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef key_MUSTANG_debug
   IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
     print *,''
     print *,'      ** COHESIVE SEDIMENTS **'
     print *,''
   END IF
#endif

   IF( .NOT. l_eroindep_mud) THEN
     !! mud erosion is proportional to total sand erosion
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      sed_eros_flxsand=0.0_rsh 
      frac_sand=0.0_rsh  
      DO iv=isand1,isand2
#ifdef key_MUSTANG_bedload
         sed_eros_flxsand=sed_eros_flxsand+sed_eros_flx_class_by_class(iv)       &
                      +ABS(flx_bxij(iv)/CELL_DX(i,j))+ABS(flx_byij(iv)/CELL_DY(i,j))
#else
         sed_eros_flxsand=sed_eros_flxsand+sed_eros_flx_class_by_class(iv)
#endif
         frac_sand=frac_sand+frac_sed(iv)
      ENDDO

      DO iv=imud1,imud2
        sed_eros_flx_class_by_class(iv)=(frac_sed(iv)/(frac_sand+0.01_rsh))*sed_eros_flxsand
#if defined key_MUSTANG_specif_outputs        
        varspecif3Dnv_save(1,iv,i,j)=0.0_rsh  ! toce_save
#endif
      END DO
      
    ELSE
     !! mud erosion is independent of sand erosion
     !! depend of tau_cri and choised options
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      sommud=0.0_rsh
!      DO iv=imud1,imud2
!        sommud=sommud+cv_sed(iv,ksmax,i,j)
!      ENDDO
!      cvolinigrv=0.0_rsh
!      DO iv=igrav1,igrav2
!        cvolinigrv=cvolinigrv+cv_sed(iv,ksmax,i,j)/ros(iv)
!      ENDDO
!      cvolinisan=0.0_rsh
!      DO iv=isand1,isand2
!        cvolinisan=cvolinisan+cv_sed(iv,ksmax,i,j)/ros(iv)
!      ENDDO

!      cmudr=sommud/(1.0_rsh-cvolinisan-cvolinigrv)
       cmudr=crel_mud(ksmax,i,j)

        !!!! *** Critical shear stress *** !!!!

        ! 4 options to compute critical shear stress for mud erosion
        ! 0- Default --> =f(cmudr)
        IF (tau_cri_mud_option_eroindep .EQ. 0) THEN
          tauc_mud=(x1toce_mud*cmudr**x2toce_mud)
          tauc_mud=max(tauc_mud,0.00005_rsh)
        END IF

        ! 1- 
        IF (tau_cri_mud_option_eroindep .EQ. 1) THEN

           somsan=0.0_rsh
           toce_meansan=0.0_rsh
           DO iv=isand1,isand2
             somsan=somsan+cv_sed(iv,ksmax,i,j)
             toce_meansan=toce_meansan+toce_loc(iv)*cv_sed(iv,ksmax,i,j)
           END DO

           IF (somsan .GT. epsilon_MUSTANG) THEN
             tauc_mud=toce_meansan/somsan
           ELSE
             tauc_mud=(x1toce_mud*cmudr**x2toce_mud) + 0.00005_rsh
           END IF

        END IF
   
        ! 2-
        IF (tau_cri_mud_option_eroindep .EQ. 2) THEN
          toce_sand(:)=0.0_rsh
          DO iv=isand1,isand2
           toce_sand(iv)=(toce_loc(iv)*cv_sed(iv,ksmax,i,j))/(cv_sed(iv,ksmax,i,j)+epsilon_MUSTANG)
          END DO
          toce_minsan=minval(toce_sand)
          IF (toce_minsan .NE. 0.0_rsh) THEN
           tauc_mud=toce_minsan
          ELSE
           tauc_mud=(x1toce_mud*cmudr**x2toce_mud) + 0.00005_rsh
          END IF
        END IF

        ! 3-  
        IF (tau_cri_mud_option_eroindep .EQ. 3) THEN
           tauc_mud=min( (x1toce_mud*cmudr**x2toce_mud) + 0.00005_rsh, stresscri0(isand2))
        END IF

        !!!! *** Erodibility parameter *** !!!!

        E0_mud_loc=E0_mud_para_indep*E0_mud

        DO iv=imud1,imud2
         sed_eros_flx_class_by_class(iv)=MF*fwet(i,j)*frac_sed(iv)*E0_mud_loc*  &
                     MAX((tauskin(i,j)/tauc_mud)-1.0_rsh,0.0_rsh)**n_eros_mud
     
#if defined key_MUSTANG_specif_outputs        
         varspecif3Dnv_save(1,iv,i,j)=tauc_mud  ! toce_save
#endif

#ifdef key_MUSTANG_debug
         IF ( l_debug_erosion .AND. CURRENT_TIME> t_start_debug  &
               .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug ) THEN
             print *,'       - tauc_mud=',tauc_mud
             print *,'       - sed_eros_flx_class_by_class(iv)=',sed_eros_flx_class_by_class(iv)
         END IF
#endif
        END DO

    ENDIF
    
   !print *,'Exit MUSTANGV2_comp_eros_flx_indep'
   !print *,'************************************'
   !print *,''

END SUBROUTINE MUSTANGV2_comp_eros_flx_indep
  !!==============================================================================
  
  SUBROUTINE  MUSTANGV2_borne_and_apply_erosion_tot(i,j,ksmax,flx_bxij,flx_byij, &
#if ! defined key_noTSdiss_insed 
                        flx_s2w_eroij,                                &
#endif
#if ! defined key_nofluxwat_IWS
                        phieau_ero_ij,                                &
#endif
                       sed_eros_flx_class_by_class,dt1,dt_ero)

   !&E--------------------------------------------------------------------------
   !&E                 *** MUSTANGV2_borne_and_apply_erosion_tot  ***
   !&E
   !&E ** Purpose : limitation of erosion and bedload extraction fluxes 
   !&E              that can not exceed the mass available in the surface layer
   !&E
   !&E ** Description : 0D
   !&E         variables IN : ksmax,dt1, dzs, cv_sed, crel_mud
   !&E                        ros 
   !&E                        
   !&E                        
   !&E                        
   !&E
   !&E         variables OUT : flx_bxij,flx_byij,sed_eros_flx_class_by_class
   !&E                          dt_ero, crel_mud, poro, dzs, cv_sed, dzsmin,
   !&E                          c_sedtot, poro_mud, ksmax
   !&E
   !&E ** Called by :  sed_erosion
   !&E
   !&E--------------------------------------------------------------------------

      !! * Arguments
      INTEGER,INTENT(IN) :: i,j
      INTEGER,INTENT(INOUT) :: ksmax
      REAL(KIND=rsh),INTENT(IN) :: dt1
      REAL(KIND=rsh),DIMENSION(1:nvp),INTENT(INOUT) :: flx_bxij,flx_byij,sed_eros_flx_class_by_class
      REAL(KIND=rsh),DIMENSION(1:nvp),INTENT(OUT)          :: dt_ero
#if ! defined key_noTSdiss_insed 
      REAL(KIND=rsh),DIMENSION(-1:nv_adv),INTENT(INOUT)   :: flx_s2w_eroij
#endif
#if ! defined key_nofluxwat_IWS
      REAL(KIND=rsh),INTENT(INOUT)                         :: phieau_ero_ij
#endif




      !! * Local declarations
      INTEGER:: iv,ivp_assoc
      REAL(KIND=rsh) :: flx_tot,dzsa,mass_tot
      REAL(KIND=rsh),DIMENSION(1:nvp) :: mass_avail,ero_tot,massinactivlayer,massinactivlayer_ini
      REAL(KIND=rsh),DIMENSION(1:nvpc) :: frac_sed
      LOGICAL,DIMENSION(1:nvp) :: l_empty
#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
      REAL(KIND=rsh) :: poroa,porowater_new,porowater1
#endif

   !!----------------------------------------------------------------------


     !print *,''
     !print *,'************************************'
     !print *,'   > Enter MUSTANGV2_borne_and_apply_erosion_tot'
     !print *,''

     dzsa=dzs(ksmax,i,j)
#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
     poroa=poro(ksmax,i,j)
#endif
     ! Initial mass in active layer for sediment class iv
     DO iv=1,nvp
       massinactivlayer(iv)=cv_sed(iv,ksmax,i,j)*dzsa
     END DO

     massinactivlayer_ini(:)=massinactivlayer(:)

#ifdef key_MUSTANG_debug
     IF ( l_debug_erosion .AND. CURRENT_TIME> t_start_debug   &
               .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug ) THEN
       print *,'    Carac of ksmax layer before borning/applying erosion : '
       print *,'      ksmax=',ksmax
       print *,'      dzs(ksmax,i,j)=',dzs(ksmax,i,j)
       print *,'      cv_sed(:,ksmax,i,j)=',cv_sed(:,ksmax,i,j)
       print *,'      c_sedtot(ksmax,i,j)=',c_sedtot(ksmax,i,j)
       print *,'      poro(ksmax,i,j)=',poro(ksmax,i,j),'poro_mud(ksmax,i,j)=',poro_mud(ksmax,i,j)
       print *,'      massinactivlayer_ini(:)=',massinactivlayer_ini(:)
       print *,'      crel_mud(ksmax,i,j)=',crel_mud(ksmax,i,j)
      ! print *,'        > CELL_SURF(i,j)=',CELL_SURF(i,j)
      ! print *,'        > CELL_DX(i,j)=',CELL_DX(i,j)
      ! print *,'        > CELL_DY(i,j)=',CELL_DY(i,j)
     END IF
#endif

     l_empty(:)=.FALSE.
     ero_tot(:)=0.0_rsh ! needs to be initialised

     DO iv=1,nvpc

#ifdef key_MUSTANG_debug
         IF ( l_debug_erosion .AND. CURRENT_TIME> t_start_debug   &
               .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug ) THEN
         print *,'      iv=',iv
       END IF
#endif

       IF (cv_sed(iv,ksmax,i,j) .LT. 10e-3) THEN ! Bof ?

         flx_bxij(iv)=0.0_rsh
         flx_byij(iv)=0.0_rsh
         sed_eros_flx_class_by_class(iv)=0.0_rsh

         dt_ero(iv)=0.0_rsh
         !dt_ero(iv)=dt1 

#ifdef key_MUSTANG_debug
         IF ( l_debug_erosion .AND. CURRENT_TIME> t_start_debug   &
               .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug ) THEN
           print *,'        > cv_sed(iv,ksmax,i,j) .LT. 10e-3 ==> flx_bx(:), flx_by(:), sed_eros(:)= 0. pour tt iv / dt_ero(:)=dt1'
         END IF
#endif

       ELSE
 
         flx_tot=   ABS(flx_bxij(iv)*CELL_DY(i,j))   &
                  + ABS(flx_byij(iv)*CELL_DX(i,j))   &
                  + sed_eros_flx_class_by_class(iv)*CELL_SURF(i,j) ! kg/s
         mass_avail(iv)=cv_sed(iv,ksmax,i,j)*dzs(ksmax,i,j)*CELL_SURF(i,j)

         IF  (flx_tot*dt1 .GE. mass_avail(iv)) THEN
           l_empty(iv)=.TRUE.
           dt_ero(iv)=mass_avail(iv)/flx_tot
         ELSE
           !l_empty(iv)=.FALSE.
           dt_ero(iv)=dt1
         END IF

#ifdef key_MUSTANG_debug
         IF ( l_debug_erosion .AND. CURRENT_TIME> t_start_debug   &
               .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug ) THEN
           print *,'        > flx_tot=',flx_tot,' mass_avail=',mass_avail(iv)
           print *,'        > l_empty(iv)=',l_empty(iv),' / dt_ero(iv)=',dt_ero(iv),' (dt1=',dt1,')'
           print *,'        Initial fluxes : '
           print *,'        > flx_bxij(iv)=',flx_bxij(iv)
           print *,'        > flx_byij(iv)=',flx_byij(iv)
           print *,'        > sed_eros_flx_class_by_class(iv)=',sed_eros_flx_class_by_class(iv)
         END IF
#endif

         ! Updating erosion/bedload fluxes according to available sediment masses in active layer !!
         ! Unit changes (kg/m2/s or kg/m/s --> kg)

         sed_eros_flx_class_by_class(iv)=sed_eros_flx_class_by_class(iv)*CELL_SURF(i,j)*dt_ero(iv) ! in kg

#ifdef key_MUSTANG_bedload
         flx_bxij(iv)=dt_ero(iv)*flx_bxij(iv)*CELL_DY(i,j) ! in kg
         flx_byij(iv)=dt_ero(iv)*flx_byij(iv)*CELL_DX(i,j) ! in kg
#endif

#ifdef key_MUSTANG_debug
         IF ( l_debug_erosion .AND. CURRENT_TIME> t_start_debug   &
               .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug ) THEN
           print *,'        > flx_tot=',flx_tot,' mass_avail=',mass_avail(iv)
           print *,'        > l_empty(iv)=',l_empty(iv),' / dt_ero(iv)=',dt_ero(iv),' (dt1=',dt1,')'
           print *,'        Updates of fluxes in kg (*dx, dy, or surf and *dt_ero(iv): '
           print *,'        > flx_bxij(iv)=',flx_bxij(iv)
           print *,'        > flx_byij(iv)=',flx_byij(iv)
           print *,'        > sed_eros_flx_class_by_class(iv)=',sed_eros_flx_class_by_class(iv)
         END IF
#endif

         ! Updating masses in active layer according to divergence of actual (limited or not) erosion/bedload fluxes 

           ero_tot(iv)=sed_eros_flx_class_by_class(iv)
#ifdef key_MUSTANG_bedload
           ero_tot(iv)=ero_tot(iv)+ABS(flx_bxij(iv))+ABS(flx_byij(iv))
#endif
         IF (.NOT. l_empty(iv)) THEN
           massinactivlayer(iv)=(cv_sed(iv,ksmax,i,j)*dzsa) - (ero_tot(iv)/CELL_SURF(i,j))
         ELSE
           massinactivlayer(iv)=0.0_rsh
         END IF

#ifdef key_MUSTANG_debug
         IF ( l_debug_erosion .AND. CURRENT_TIME> t_start_debug   &
               .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug ) THEN
           print *,'        Updating masses in active layer according to divergence of actual (limited or not) erosion/bedload fluxes'
           print *,'        > massinactivlayer(iv)=',massinactivlayer(iv),' (l_empty(iv)=',l_empty(iv),')'
         END IF
#endif        
       ENDIF

     END DO

     DO iv=nvpc+1,nvp
#ifdef key_Pconstitonly_insed
       sed_eros_flx_class_by_class(iv)=0.0_rsh
       flx_bxij(iv)=0.0_rsh
       flx_byij(iv)=0.0_rsh
       massinactivlayer(iv)=0.0_rsh
#else
       ivp_assoc=irkm_var_assoc(iv)

       IF (ivp_assoc == 0) THEN
         sed_eros_flx_class_by_class(iv)=massinactivlayer_ini(iv)*CELL_SURF(i,j)*  &
            (1.0_rsh-( SUM(massinactivlayer(1:nvpc))/SUM(massinactivlayer_ini(:))))
         flx_bxij(iv)=0.0_rsh
         flx_byij(iv)=0.0_rsh
         massinactivlayer(iv)=massinactivlayer_ini(iv)-(sed_eros_flx_class_by_class(iv)/CELL_SURF(i,j))
       ELSE
         sed_eros_flx_class_by_class(iv)=massinactivlayer_ini(iv)*CELL_SURF(i,j)  &
             *(1.0_rsh - (massinactivlayer(ivp_assoc)/(massinactivlayer_ini(ivp_assoc)+epsi30_MUSTANG))) &
                *(sed_eros_flx_class_by_class(ivp_assoc)/(ero_tot(ivp_assoc)+epsi30_MUSTANG))
         flx_bxij(iv)=massinactivlayer_ini(iv)*CELL_SURF(i,j)*(1.0_rsh - &
                 (massinactivlayer(ivp_assoc)/(massinactivlayer_ini(ivp_assoc)+epsi30_MUSTANG))) &
                  *(flx_bxij(ivp_assoc)/(ero_tot(ivp_assoc)+epsi30_MUSTANG))
         flx_byij(iv)=massinactivlayer_ini(iv)*CELL_SURF(i,j)*(1.0_rsh - &
                 (massinactivlayer(ivp_assoc)/(massinactivlayer_ini(ivp_assoc)+epsi30_MUSTANG))) &
                  *(flx_byij(ivp_assoc)/(ero_tot(ivp_assoc)+epsi30_MUSTANG))
         massinactivlayer(iv)=massinactivlayer_ini(iv)*(massinactivlayer(ivp_assoc) &
                               /(massinactivlayer_ini(ivp_assoc)+epsi30_MUSTANG))
       END IF
#endif
     END DO


     !print *,'Masses remaining in activelayer after erosion = ',massinactivlayer(:)

     mass_tot=0.0_rsh
     DO iv=1,imud2
       mass_tot=mass_tot+massinactivlayer(iv)
     END DO

#ifdef key_MUSTANG_debug
         IF ( l_debug_erosion .AND. CURRENT_TIME> t_start_debug   &
               .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug ) THEN
       print *,'    Update of the characteristics of the new ksmax layer after erosion : '
       print *,'      > massinactivlayer(:)=',massinactivlayer(:)
       print *,'      > mass_tot=',mass_tot
     END IF
#endif
     !print *,'massinactivlayer(:)=',massinactivlayer(:)


     IF (mass_tot .GT. 0.0_rsh) THEN
       !print *,'Remaining masses, update of the ksmax layer composition and thickness'

       DO iv=1,imud2
         frac_sed(iv)=massinactivlayer(iv)/mass_tot
       END DO

       CALL MUSTANGV2_comp_poro_mixsed(frac_sed, poro_mud(ksmax,i,j),  &
                                crel_mud(ksmax,i,j), poro(ksmax,i,j))

       !IF (crel_mud(ksmax,i,j) .GT. 1500.0_rsh) THEN
       !  print *,'in sed_erosion_mixsed 2'
       !  print *,' > crel_mud(ksmax,i,j) = ',crel_mud(ksmax,i,j)
       !  print *,' > massinactivlayer(iv) = ',massinactivlayer(:)
       !  print *,' > frac_sed(iv) = ',frac_sed(:)
       !END IF


       dzs(ksmax,i,j)=mass_tot/((1.0_rsh-poro(ksmax,i,j))*ros(1))

!       dzsmin=(1.0_rsh-coeff_dzsmin)*dzsminuni + coeff_dzsmin*SUM( (cv_sed(1:nvpc,ksmax,i,j)/c_sedtot(ksmax,i,j))*diam_sed(1:nvpc) ) ! NON
       !dzsmin=dzsminuni + coeff_dzsmin*SUM( (cv_sed(1:isand2,k,i,j)/c_sedtot(k,i,j))*diam_sed(1:isand2) ) ??
       !dzsmin=(1.0_rsh-coeff_dzsmin)*dzsminuni + coeff_dzsmin*SUM( frac_sed(1:nvpc)*diam_sed(1:nvpc) )
       dzsmin=dzsminvar(frac_sed)

#ifdef key_MUSTANG_debug
         IF ( l_debug_erosion .AND. CURRENT_TIME> t_start_debug   &
               .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug ) THEN
         print *,'    mass_tot > 0'
         print *,'      > CALL MUSTANGV2_comp_poro_mixsed(frac_sed, poro_mud(ksmax,i,j), poro(ksmax,i,j))'
         print *,'           - frac_sed(:)=',frac_sed(:)
         print *,'           - poro_mud(ksmax,i,j)=',poro_mud(ksmax,i,j)
         print *,'           - OUT ==> poro(ksmax,i,j)=',poro(ksmax,i,j)
         print *,'           - OUT ==> crel_mud(ksmax,i,j)=',crel_mud(ksmax,i,j)
         print *,'      > dzs(ksmax,i,j)=',dzs(ksmax,i,j)
         print *,'      > dzsmin=',dzsmin
       END IF
#endif

#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
         ! dissolved variable in pore waters and water fluxes at the interface
         ! --------------------------------------------------------------------
         !  ==> flx_s2w and cv_sed and phieau_s2w
         ! in the active layer which has been modified by erosion and bedload ==> code 3
         ! dissolved concentrations in this layer unchanged by erosion and bedload 
         !  because concentrations in activelayer have been already mixed and = cw_bottom during fusion
          porowater_new=poro(ksmax,i,j)*dzs(ksmax,i,j)
          porowater1=dzsa*poroa
          CALL MUSTANGV2_eval_dissvar_IWSflux(i,j,ksmax,3,   &
              phieau_ero_ij=phieau_ero_ij,          &
              flx_s2w_eroij=flx_s2w_eroij,          &
              porowater1=porowater1,                &
              cv_sed1=cv_sed(-1:nv_adv,ksmax,i,j),  &
              porowater_new=porowater_new)
#endif

       IF (dzs(ksmax,i,j) .GT. dzsmin) THEN

         DO iv=1,nvp
           cv_sed(iv,ksmax,i,j)=massinactivlayer(iv)/dzs(ksmax,i,j)
         END DO

         c_sedtot(ksmax,i,j)=mass_tot/dzs(ksmax,i,j)

       ELSE

#ifdef key_MUSTANG_debug
         IF ( l_debug_erosion .AND. CURRENT_TIME> t_start_debug   &
               .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug ) THEN
           print *,'      > CALL MUSTANGV2_manage_small_mass_in_ksmax(i,j,ksmax,massinactivlayer)'
         END IF
#endif

         CALL MUSTANGV2_manage_small_mass_in_ksmax(i,j,ksmax,   &
#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
                                phieau_ero_ij,flx_s2w_eroij,   &
#endif
                                massinactivlayer)


       END IF

     ELSE
       !print *,'Total erosion of the ksmax layer'
       ! et que deviennent les particulaires non constitutives ?

#ifdef key_MUSTANG_debug
         IF ( l_debug_erosion .AND. CURRENT_TIME> t_start_debug   &
               .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug ) THEN
         print *,'    mass_tot = 0 ==> total erosion of the layer'
       END IF
#endif
#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
         ! dissolved variable in pore waters and water fluxes at the interface
         ! --------------------------------------------------------------------
         !  ==> flx_s2w and phieau_s2w
         ! layer elimination ==> code 1
          porowater1=dzsa*poroa
          CALL MUSTANGV2_eval_dissvar_IWSflux(i,j,ksmax,1,   &
                           phieau_ero_ij=phieau_ero_ij,          &
                           flx_s2w_eroij=flx_s2w_eroij,          &
                           porowater1=porowater1,                 &
                           cv_sed1=cv_sed(-1:nv_adv,ksmax,i,j))
#endif
       dzs(ksmax,i,j)=0.0_rsh
       cv_sed(:,ksmax,i,j)=0.0_rsh
       c_sedtot(ksmax,i,j)=0.0_rsh
       poro(ksmax,i,j)=0.0_rsh
       poro_mud(ksmax,i,j)=0.0_rsh
       crel_mud(ksmax,i,j)=0.0_rsh

       ksmax=ksmax-1


     END IF

#ifdef key_MUSTANG_debug
         IF ( l_debug_erosion .AND. CURRENT_TIME> t_start_debug   &
               .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug ) THEN
       print *,'    Carac of ksmax layer after borning/applying erosion : '
       print *,'      ksmax=',ksmax
       print *,'      dzs(ksmax,i,j)=',dzs(ksmax,i,j)
       print *,'      cv_sed(:,ksmax,i,j)=',cv_sed(:,ksmax,i,j)
       print *,'      c_sedtot(ksmax,i,j)=',c_sedtot(ksmax,i,j)
       print *,'      poro(ksmax,i,j)=',poro(ksmax,i,j),'poro_mud(ksmax,i,j)=',poro_mud(ksmax,i,j)
       print *,'      crel_mud(ksmax,i,j)=',crel_mud(ksmax,i,j)
     END IF
#endif

! To see later
!     DO iv=imud2+1,nvp
!       ivp_assoc=irkm_var_assoc(iv)
!#ifdef key_MUSTANG_bedload
!       flx_bxij(iv)=flx_bxij(ivp_assoc)*cv_sed(iv,ksmax,i,j)/cv_sed(ivp_assoc,ksmax,i,j)
!       flx_byij(iv)=flx_byij(ivp_assoc)*cv_sed(iv,ksmax,i,j)/cv_sed(ivp_assoc,ksmax,i,j)
!#endif
!       sed_eros_flx_class_by_class(iv)=sed_eros_flx_class_by_class(ivp_assoc)*cv_sed(iv,ksmax,i,j)/cv_sed(ivp_assoc,ksmax,i,j)
!     END DO


    !print *,''
    !print *,'Exit sed_borne_and_apply_erosion_tot'
    !print *,'************************************'
    !print *,''

END SUBROUTINE MUSTANGV2_borne_and_apply_erosion_tot

  !!==============================================================================
SUBROUTINE MUSTANGV2_manage_small_mass_in_ksmax(i,j,ksmax,      &
#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
                                phieau_ero_ij,flx_s2w_eroij,   &
#endif
                                mass_sed_ijksmax)

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE MUSTANGV2_manage_small_mass_in_ksmax  ***
   !&E
   !&E ** Purpose : when sediment masses remaining in ksmax become too small,
   !&E              they are integrated in ksmax-1 layer
   !&E
   !&E ** Description : 0D
   !&E              variables IN : ksmax, dzs, cv_sed, poro, poro_mud, crel_mud
   !&E                             c_sedtot
   !&E
   !&E              variables OUT : ksmax,dzs, cv_sed, poro, poro_mud, crel_mud
   !&E                             c_sedtot
   !&E
   !&E ** Called by :  MUSTANGV2_comp_poro_mixsed
   !&E
   !&E--------------------------------------------------------------------------

   !! * Modules used

   !! * Arguments
   INTEGER, INTENT(IN) :: i,j
   INTEGER, INTENT(INOUT) :: ksmax
   REAL(KIND=rsh), DIMENSION(nvp), INTENT(IN) :: mass_sed_ijksmax

#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
   REAL(KIND=rsh),DIMENSION(-1:nv_adv),INTENT(INOUT)           :: flx_s2w_eroij                       
   REAL(KIND=rsh),INTENT(INOUT)                                :: phieau_ero_ij                    
#endif

   !! * Local declaration
   INTEGER :: iv,ivp_assoc
   REAL(KIND=rsh) :: poro_mud_k,poro_mud_km1,dzsam1, &
                     massmud_k,massmud_km1,massmud_tot,   &
                     poroa,poroam1,poro_mud_new,mass_tot, &
                     dzsi,crel_mud_new, poro_kij

   REAL(KIND=rsh),DIMENSION(nvp) :: mass_sed
   REAL(KIND=rsh),DIMENSION(nvpc):: frac_sed
#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
   REAL(KIND=rsh) :: porowater_new,porowater1,porowater2,dzsa
#endif
   !!---------------------------------------------------------------------------
   !! * Executable part

   !print *,''
   !print *,'*****************************************************'
   !print *,'   > Enter MUSTANGV2_manage_small_mass_in_ksmax'
   !print *,'*****************************************************'
   !print *,''

     ! 1- We look at the composition and the porosity of the mud layers ksmax and ksmax-1
     ! 2- A porosity of the representative mud is calculated (pro rata of the masses)
     ! 3- The sediment masses / fractions of the fused layer are calculated
     ! 4- We calculate the porosity -> ksmax = ksmax-1 -> MAJ dzs, cv_sed ksmax
     
     ! 1- On regarde la composition et la porosite de la vase des couches ksmax et ksmax-1
     ! 2- On calcule une porosite de la vase representative (au prorata des masses) 
     ! 3- On calcule les masses/fractions sedimentaires de la couche fusionnee
     ! 4- On calcule la porosite --> ksmax=ksmax-1 --> MAJ dzs, cv_sed ksmax

     !1-
     massmud_k=0.0_rsh
     massmud_km1=0.0_rsh
     DO iv=imud1,imud2
       massmud_k=massmud_k+mass_sed_ijksmax(iv)
       massmud_km1=massmud_km1+(cv_sed(iv,ksmax-1,i,j)*dzs(ksmax-1,i,j))
     ENDDO

     massmud_tot=massmud_k+massmud_km1

     poro_mud_k=poro_mud(ksmax,i,j)
     poro_mud_km1=poro_mud(ksmax-1,i,j)

     poroa=poro(ksmax,i,j) ! needs to be updated according to remaining sediments before calling the subroutine
     poroam1=poro(ksmax-1,i,j)

     dzsam1=dzs(ksmax-1,i,j)
#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
     dzsa=dzs(ksmax,i,j)
     porowater1=dzsa*poroa
     porowater2=dzsam1*poroam1
#endif
     !print *,'**Couche ksmax**'
     !print *,'ksmax=',ksmax
     !print *,'mass_sed_ijksmax(:)=',mass_sed_ijksmax(:)
     !print *,'poro_mud en ksmax=',poro_mud_k
     !print *,'poro(ksmax,i,j)=',poroa
     !print *,''
     !print *,'**Couche ksmax-1**'
     !print *,'ksmax-1=',ksmax-1
     !print *,'dzs(ksmax-1,i,j)=',dzs(ksmax-1,i,j)
     !print *,'cv_sed(:,ksmax-1,i,j)=',cv_sed(:,ksmax-1,i,j)
     !print *,'poro_mud en ksmax-1=',poro_mud_km1
     !print *,'poro(ksmax-1,i,j)=',poroam1
     !print *,''

     !2-
     IF (massmud_tot .GT. 0.0_rsh) THEN
       poro_mud_new=(massmud_k*poro_mud_k + massmud_km1*poro_mud_km1)/massmud_tot
       crel_mud_new=(massmud_k*crel_mud(ksmax,i,j) + massmud_km1*crel_mud(ksmax-1,i,j))/massmud_tot
     ELSE
       poro_mud_new=0.0_rsh
       crel_mud_new=0.0_rsh
     END IF

     !3-
     mass_tot = 0.0_rsh
     DO iv=igrav1,imud2
       mass_sed(iv)=mass_sed_ijksmax(iv) + (cv_sed(iv,ksmax-1,i,j)*dzsam1)
       mass_tot=mass_tot+mass_sed(iv)
     END DO

     DO iv=igrav1,imud2
       frac_sed(iv)=mass_sed(iv)/mass_tot
     END DO

     !4-
     CALL MUSTANGV2_comp_poro_mixsed(frac_sed, poro_mud_new,    &
                                   crel_mud_new, poro_kij)

     ! Question : faut il prescrire a la toute fin que dzs et cv_sed en ksmax avant fusion
     ! sont desormais = 0.0 (ou -999.0)

     ksmax=ksmax-1

     dzs(ksmax,i,j)=mass_tot/((1.0_rsh-poro_kij)*ros(1))
     dzsi=1.0_rsh/dzs(ksmax,i,j)

     crel_mud(ksmax,i,j)=crel_mud_new
     poro(ksmax,i,j)=poro_kij
     poro_mud(ksmax,i,j)=poro_mud_new ! Actualisation de poro_mud liee au melange !

     c_sedtot(ksmax,i,j)=0.0_rsh
     DO iv=1,nvpc
       cv_sed(iv,ksmax,i,j)=mass_sed(iv)*dzsi
       c_sedtot(ksmax,i,j)=c_sedtot(ksmax,i,j)+cv_sed(iv,ksmax,i,j)
     END DO

     !IF (crel_mud(ksmax,i,j) .GT. 1500.0_rsh) THEN
     !  print *,'in sed_fusion_with_poro'
     !  print *,' > crel_mud(ksmax,i,j) = ',crel_mud(ksmax,i,j)
     !  print *,' > cv_sed(:,ksmax,i,j) = ',cv_sed(:,ksmax,i,j)
     !END IF


     DO iv=nvpc+1,nvp
#ifdef key_Pconstitonly_insed
       cv_sed(iv,ksmax,i,j)=0.0_rsh
#else
       cv_sed(iv,ksmax,i,j)=(mass_sed_ijksmax(iv)+cv_sed(iv,ksmax,i,j)*dzsam1)*dzsi
#endif
     ENDDO


#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
     ! dissolved variable in pore waters and water fluxes at the interface
     ! --------------------------------------------------------------------
     ! fusion of 2 layers at surface ==> code 8  new poro could introduce water inputs or output
     porowater_new=poro(ksmax,i,j)*dzs(ksmax,i,j)
     CALL MUSTANGV2_eval_dissvar_IWSflux(i,j,ksmax,8,  &
           phieau_ero_ij=phieau_ero_ij,           &
           flx_s2w_eroij=flx_s2w_eroij,           &
           porowater1=porowater1,                      &
           cv_sed1=cv_sed(-1:nv_adv,ksmax+1,i,j),        &
           porowater2=porowater2,                      &
           cv_sed2=cv_sed(-1:nv_adv,ksmax,i,j),      &
           porowater_new=porowater_new)
#endif

     ! old layer ksmax =0
     poro(ksmax+1,i,j)=0.0_rsh
     poro_mud(ksmax+1,i,j)=0.0_rsh
     dzs(ksmax+1,i,j)=0.0_rsh
     cv_sed(:,ksmax+1,i,j)=0.0_rsh
     
     !print *,'Apres fusion'
     !print *,'ksmax=',ksmax
     !print *,'dzs(ksmax,i,j)=',dzs(ksmax,i,j)
     !print *,'cv_sed(:,ksmax,i,j)=',cv_sed(:,ksmax,i,j)
     !print *,'poro_mud_new en ksmax=',poro_mud_new
     !print *,'poro(ksmax,i,j)=',poro(ksmax,i,j)
     !print *,''

   !print *,'Exit MUSTANGV2_manage_small_mass_in_ksmax'
   !print *,'************************************'
   !print *,''


END SUBROUTINE MUSTANGV2_manage_small_mass_in_ksmax

   !!==============================================================================
  SUBROUTINE MUSTANGV2_comp_poro_mixsed(frac_sed_kij, poro_mud_kij, &
                                crel_mud_kij, poro_kij)

   !&E--------------------------------------------------------------------------
   !&E                 ***  SUBROUTINE MUSTANGV2_comp_poro_mixsed  ***
   !&E
   !&E ** Purpose : compute porosity in case of non-cohesive/cohesive sediment mixtures
   !&E
   !&E ** Description : Wooster et al. (2008) and Wu and Li (2017)
   !&E                  variables OUT : poro,crel_mud
   !&E                  variables IN : frac_sed,psi_sed,Awooster,Bwooster,Bmax_wu      
   !&E                                  poro_option,diam_sed,poro_mud,crel_mud,
   !&E
   !&E ** Note : NUMBER_PI must be known as a parameters transmtted by coupleur 
   !&E           in MARS : coupleur_dimhydro.h (USE ..)
   !&E           in CROCO : module_MUSTANG.F (include..)
   !&E
   !&E
   !&E ** Called by :  fusion_with_poro, manage_small_mass__in_ksmax
   !&E                 borne_and_apply_erosion_tot, effdep, erosion
   !&E
   !&E--------------------------------------------------------------------------

    !! * Arguments
    REAL(KIND=rsh),DIMENSION(1:imud2), INTENT(IN)    :: frac_sed_kij
    REAL(KIND=rsh),                    INTENT(INOUT)    :: poro_mud_kij
    REAL(KIND=rsh),                    INTENT(INOUT) :: crel_mud_kij
    REAL(KIND=rsh),                    INTENT(OUT)   :: poro_kij

    !! * Local declarations
    INTEGER         ::  iv
    REAL(KIND=rsh)  ::  psim,siggeo2,siggeo,stdgeo, &
                        frac_gravsan,diam_gravsan,poro_gravsan, &
                        frac_mud,beta,Nc,n,Rn,P2,Bi, &
                        Term1,Term2,Term3,SumTerm,frac_mud_cr, &
                        f1b,f1c,poro_icp,f2b,f2c,poro_ifp,diam_mud

    REAL(KIND=rsh),DIMENSION(igrav1:isand2) :: frac_noncoh

   !!---------------------------------------------------------------------------
   !! * Executable part


    !print *,''
    !print *,' --> enter MUSTANGV2_comp_poro_mixsed'
    !print *,''


    !!!!!! Step I - Computation of the porosity associated to the non-cohesive sediment
    !               fraction (i.e. igrav1 --> isand2), using the experimental relation derived
    !               by Wooster et al. (2008) between geometric standard deviation and porosity

    !               REF : Sediment supply and relative size distribution effects on fine sediment
    !                     infiltration into immobile gravels
    !
    !               Remarks : - see Cui et al. (1996) for more details concerning the geometric standard deviation 
    !                           computation (hereafter referred to as stdgeo)
    !                         - Awooster = 0.621 and Bwooster = -0.659 in the original paper of Wooster et al. (2008)
    !                         - Representative diameter of muddy particles is set at 31.5 microns (considering a %mud = %<63 microns)

    frac_gravsan=0.0_rsh
    DO iv=igrav1,isand2
      frac_gravsan=frac_gravsan+frac_sed_kij(iv)
    END DO 

    IF (frac_gravsan .GT. epsilon_MUSTANG) THEN

      frac_noncoh(:)=0.0_rsh
      DO iv=igrav1,isand2
        IF (frac_sed_kij(iv) .GT. 0.0001_rsh*epsilon_MUSTANG) THEN
          frac_noncoh(iv)=frac_sed_kij(iv)/frac_gravsan
        END IF
      END DO

      psim=0.0_rsh
      DO iv=igrav1,isand2
        psim=psim+(psi_sed(iv)*frac_noncoh(iv))
      END DO

      siggeo2=0.0_rsh
      DO iv=igrav1,isand2
        siggeo2=siggeo2+( frac_noncoh(iv)*(psi_sed(iv)-psim)**2.0_rsh )
      END DO

      siggeo  = siggeo2**(0.5_rsh)
      stdgeo  = MAX(1.0_rsh,MIN(2.0_rsh**(siggeo),4.0_rsh)) ! to remain in ranges of values described by Wooster et al. (2008)

      diam_gravsan = (2.0_rsh**(psim))/1000.0_rsh !psi_sed is computed with diam_sed in mm --> /1e3
      poro_gravsan = Awooster*(stdgeo**(Bwooster))

      !!!!!!

      frac_mud=0.0_rsh
      DO iv=imud1,imud2
        frac_mud=frac_mud+frac_sed_kij(iv)
      END DO

     ! modif plehir nov 2019
     ! IF (frac_mud .GT. 0.001_rsh .AND. crel_mud_kij .GT. 0.001_rsh) THEN
      frac_mud_cr=crel_mud_kij/ros(1)/(crel_mud_kij/ros(1)+(1-poro_gravsan)/poro_gravsan)
      IF(frac_mud .GT. frac_mud_cr) THEN

        IF (poro_option .EQ. 1) THEN

          !!!!!! Step II - Computation of the porosity representative of the non-cohesive / cohesive mixture
    
          !                REF : Wu and Li (2017); Porosity of bimodal sediment mixture with particle filling

          ! calcul de diam_mud
          diam_mud=0.0_rsh
          DO iv=imud1,imud2
            diam_mud=diam_mud+diam_sed(iv)*frac_sed_kij(iv)
          ENDDO
  
          ! beta : In the configuration where fines particles cover the surface
          !        of a coarse particle, beta represents the arc angle linking the two 
          !        centers of adjacent fine particles
          beta = 2.0_rsh*ASIN(diam_mud/(diam_gravsan+diam_mud))

          ! Nc : Number of fines particles needed to cover the surface of a coarse
          !      particle
          Nc = 4.0_rsh*NUMBER_PI/( (beta**2.0_rsh) * COS(NUMBER_PI/6.0_rsh) )

          ! n : Minimum number of layers of fine particles required to fill up 
          !     the voids of coarse particles
          n = 0.1124_rsh*(diam_gravsan/diam_mud)
          !print *,'n',n

          ! Rn : The number of fine particles available to cover a coarse
          !      particle, which is equal to the ratio of the numbers of
          !      fine and coarse particles in the sediment mixture
          Rn = (frac_mud*diam_gravsan**3.)/(frac_gravsan*diam_mud**3.)
          !print *,'Rn',Rn
  
          ! P2 : The probability of a coarse particle contacting with a
          !      fine particle
          P2 = (frac_mud/diam_mud)/((frac_gravsan/diam_gravsan) + (frac_mud/diam_mud))

          ! Bi : The portion of the coarse sediment class participating 
          !      in filling
          !!! Envisager de revoir formulation de Bi dans cas ou sable heterometrique, en faisant intervenir ecart type geometrique !!!
          Bi = MIN( (P2*Rn)/(n*Nc) , Bmax_wu )

          Term1 = (frac_gravsan/(1.0_rsh-poro_gravsan))*(1.0_rsh-Bi) ! Volume of the coarse size class without filling
          Term2 = frac_gravsan*Bi                                    ! Volume of the coarse size class with filling
          Term3 = (frac_mud/(1.0_rsh-poro_mud_kij))                  ! Volume of the fine size class

          SumTerm = Term1+Term2+Term3
          poro_kij = 1.0_rsh - (1.0_rsh/SumTerm)


        ELSE IF (poro_option .EQ. 2) THEN

          ! modif plehir automn 2019
          !crel_mud_kij=MIN( crel_mud_kij, ros(1)* ( (1.0_rsh-poro_gravsan)/poro_gravsan )* ( (1.0_rsh/frac_gravsan)-1.0_rsh ) )
          !poro_kij = (ros(1)-crel_mud_kij)/(ros(1)+crel_mud_kij*((1.0_rsh/frac_mud) -1.0_rsh))
          poro_kij=(ros(1)-crel_mud_kij)/(ros(1)+crel_mud_kij*((1.0_rsh/frac_mud)-1.0_rsh))

        END IF

      ELSE

        ! if no mud, only gravel and sand ==> poro_kij=poro_gravsan
        !poro_kij=poro_gravsan
        !crel_mud_kij=0.0_rsh

        ! modif plehir automn 2019
        ! faible fraction vaseuse, on reste en ideal coarse packing
        crel_mud_kij=ros(1)*((1.0_rsh- poro_gravsan)/poro_gravsan)*frac_mud/(1-frac_mud)
        poro_kij=poro_gravsan-(1-poro_gravsan)*frac_mud/(1-frac_mud)
       ! poro_mud_kij=1.0_rsh/ros(1)
# if defined key_ANA_bedload ||  defined ANA_DUNE  ||  defined DUNE 
        poro_kij=0.4
# endif
      END IF

    ELSE 

      ! if mud only, poro=poro_mud
      IF (poro_option .EQ. 2) THEN
        poro_kij = MIN(1.0_rsh-(crel_mud_kij/ros(1)), 1.0_rsh-(cfreshmud/ros(1)))
      ELSE
        poro_kij = poro_mud_kij
      END IF

    END IF

    IF (poro_option .EQ. 2) poro_kij = MAX(poro_kij, poro_min)

    IF (poro_option .EQ. 2 .AND. poro_kij == 1.0_rsh) THEN
      print *,'in MUSTANGV2_comp_poro_mixsed poro=1'
      print *,'  > frac_gravsan=',frac_gravsan
      print *,'  > poro_gravsan=',poro_gravsan
      print *,'  > frac_mud=',frac_mud
      print *,'  > poro_mud_kij=',poro_mud_kij
      print *,'  > crel_mud_kij=',crel_mud_kij
      print *,'  > poro_icp=',poro_icp
      print *,'  > poro_ifp=',poro_ifp
      print *,'  > poro_kij=',poro_kij      
    END IF

    !print *,''
    !print *,' --> exit MUSTANGV2_comp_poro_mixsed'
    !print *,''


  END SUBROUTINE MUSTANGV2_comp_poro_mixsed

  !!============================================================================== 
 
#if ! defined key_noTSdiss_insed || ! defined key_nofluxwat_IWS
  SUBROUTINE MUSTANGV2_eval_dissvar_IWSflux(i,j,k,code,phieau_ero_ij,    &
                                            flx_s2w_eroij,porowater1,cv_sed1, &
                                            porowater2,cv_sed2,porowater_new)

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE MUSTANGV2_eval_dissvar_IWSflux  ***
   !&E
   !&E ** Purpose : Compute concentrations of dissolved variables in interstitial waters
   !&E              with water fluxes and variables fluxes at the interface
   !&E
   !&E ** Description : 0D
   !&E                variables IN :  k,code,porowater1,cv_sed1, 
   !&E                                            porowater2,cv_sed2,
   !&E                                            porowater_new 
   !&E                variables OUT : cv_sed, flx_s2w,flx_w2s,phieau_s2w
   !&E
   !&E  ATTENTION : MF is not taking into account - 
   !&E              morphodynamic with MF different to 1 is not compatible with dissolved variables and interstitial water
   !&E
   !&E ** Called by :  sed_erosion,..
   !&E
   !&E--------------------------------------------------------------------------


   !! * Modules used  

   !! * Arguments
   INTEGER,INTENT(IN)                                         :: i,j,k,code
   REAL(KIND=rsh),INTENT(INOUT), OPTIONAL                     :: phieau_ero_ij
   REAL(KIND=rsh),DIMENSION(-1:nv_adv),INTENT(INOUT), OPTIONAL:: flx_s2w_eroij
   REAL(KIND=rsh),INTENT(IN), OPTIONAL                        :: porowater1,porowater2,porowater_new
   REAL(KIND=rsh),DIMENSION(-1:nv_adv),INTENT(IN), OPTIONAL   :: cv_sed1,cv_sed2

   !! * Local declarations
   INTEGER        ::  iv
   REAL(KIND=rsh) :: dporow
   !!---------------------------------------------------------------------------
   !! * Executable part

   IF(code==1) THEN
     !! case where suppression of layer k
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
        !flx_s2w(nvp+1:nv_adv,i,j)=flx_s2w(nvp+1:nv_adv,i,j)+  &
        flx_s2w_eroij(nvp+1:nv_adv)=flx_s2w_eroij(nvp+1:nv_adv)+  &
                                   cv_sed1(nvp+1:nv_adv)*porowater1
#endif
        !flx_s2w(-1:0,i,j)=flx_s2w(-1:0,i,j)+  &
        flx_s2w_eroij(-1:0)=flx_s2w_eroij(-1:0)+  &
                                   cv_sed1(-1:0)*porowater1
#endif
#if ! defined key_nofluxwat_IWS
        phieau_ero_ij=phieau_ero_ij+REAL(porowater1*CELL_SURF(i,j),rlg)
#endif


   ELSE IF (code ==2 ) THEN
      !! case where 2 layers are melted (no fluxes) (not at surface or in the case of active layer)
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
        cv_sed(nvp+1:nv_adv,k,i,j)=(cv_sed1(nvp+1:nv_adv)*porowater1+ &
                                    cv_sed2(nvp+1:nv_adv)*porowater2) &
                                      /porowater_new
#endif
        cv_sed(-1:0,k,i,j)=(cv_sed1(-1:0)*porowater1+cv_sed2(-1:0)*porowater2) &
                                      /porowater_new
#endif

   ELSE IF (code==3) THEN
      !! case where a part of the surface layer is eroded 
      dporow=porowater1-porowater_new
      IF( dporow > 0.0_rsh) THEN
          ! output of porewater in sediment

#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
         !flx_s2w(nvp+1:nv_adv,i,j)=flx_s2w(nvp+1:nv_adv,i,j)+  &
         flx_s2w_eroij(nvp+1:nv_adv)=flx_s2w_eroij(nvp+1:nv_adv)+  &
                                cv_sed1(nvp+1:nv_adv)*dporow
#endif
         !flx_s2w(-1:0,i,j)=flx_s2w(-1:0,i,j)+ cv_sed1(-1:0)*dporow
         flx_s2w_eroij(-1:0)=flx_s2w_eroij(-1:0)+ cv_sed1(-1:0)*dporow
#endif
#if ! defined key_nofluxwat_IWS
         !phieau_s2w(i,j)=phieau_s2w(i,j)+REAL(dporow*CELL_SURF(i,j),rlg)
         phieau_ero_ij=phieau_ero_ij+REAL(dporow*CELL_SURF(i,j),rlg)
#endif
      ELSE
          ! input of porewater in sediment dporow <0
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
        cv_sed(nvp+1:nv_adv,k,i,j)=(cv_sed1(nvp+1:nv_adv)*porowater1  &
                                   -cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)*dporow) /porowater_new
        !flx_s2w(nvp+1:nv_adv,i,j)=flx_s2w(nvp+1:nv_adv,i,j)+ dporow*cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)
        flx_s2w_eroij(nvp+1:nv_adv)=flx_s2w_eroij(nvp+1:nv_adv)+ dporow*cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)

#endif
        cv_sed(-1,k,i,j)=(cv_sed1(-1)*porowater1  &
                                   -temp_bottom_MUSTANG(i,j)*dporow) /porowater_new
        cv_sed(0,k,i,j)=(cv_sed1(0)*porowater1  &
                                   -sal_bottom_MUSTANG(i,j)*dporow) /porowater_new
        !flx_s2w(-1,i,j)=flx_s2w(-1,i,j)+ temp_bottom_MUSTANG(i,j)*dporow
        flx_s2w_eroij(-1)=flx_s2w_eroij(-1)+ temp_bottom_MUSTANG(i,j)*dporow
        !flx_s2w(0,i,j)=flx_s2w(0,i,j)+ sal_bottom_MUSTANG(i,j)*dporow
        flx_s2w_eroij(0)=flx_s2w_eroij(0)+ sal_bottom_MUSTANG(i,j)*dporow
#endif
#if ! defined key_nofluxwat_IWS
        !phieau_s2w(i,j)=phieau_s2w(i,j)+REAL(dporow*CELL_SURF(i,j),rlg)
        phieau_ero_ij=phieau_ero_ij+REAL(dporow*CELL_SURF(i,j),rlg)
#endif
      ENDIF

   ELSE IF (code==4) THEN
      !! case where 2 layers are melted in the case of active layer formation
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
        cv_sed(nvp+1:nv_adv,k,i,j)=cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)
        !flx_s2w(nvp+1:nv_adv,i,j)=flx_s2w(nvp+1:nv_adv,i,j)+  &
        flx_s2w_eroij(nvp+1:nv_adv)=flx_s2w_eroij(nvp+1:nv_adv)+  &
             cv_sed1(nvp+1:nv_adv)*porowater1+cv_sed2(nvp+1:nv_adv)*porowater2  &
             -cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)*porowater_new
#endif
        cv_sed(0,k,i,j)=sal_bottom_MUSTANG(i,j)
        cv_sed(-1,k,i,j)=temp_bottom_MUSTANG(i,j)
        !flx_s2w(0,i,j)=flx_s2w(0,i,j)+cv_sed1(0)*porowater1+cv_sed2(0)*porowater2 &
        flx_s2w_eroij(0)=flx_s2w_eroij(0)+cv_sed1(0)*porowater1+cv_sed2(0)*porowater2 &
                            -sal_bottom_MUSTANG(i,j)*porowater_new
        !flx_s2w(-1,i,j)=flx_s2w(-1,i,j)+cv_sed1(-1)*porowater1+cv_sed2(-1)*porowater2  &
        flx_s2w_eroij(-1)=flx_s2w_eroij(-1)+cv_sed1(-1)*porowater1+cv_sed2(-1)*porowater2  &
                             -temp_bottom_MUSTANG(i,j)*porowater_new
#endif
#if ! defined key_nofluxwat_IWS
        !phieau_s2w(i,j)=phieau_s2w(i,j)+REAL((porowater1+porowater2-porowater_new)*CELL_SURF(i,j),rlg)
        phieau_ero_ij=phieau_ero_ij+REAL((porowater1+porowater2-porowater_new)*CELL_SURF(i,j),rlg)
#endif
      
    ELSE IF (code==5) THEN
      !! case where there is a deposit and the surface layer is increased
      dporow=porowater_new-porowater1
      IF( dporow > 0.0_rsh) THEN
          ! input of porewater in sediment
          
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
        cv_sed(nvp+1:nv_adv,k,i,j)=(cv_sed1(nvp+1:nv_adv)*porowater1  &
                                   +cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)*dporow) /porowater_new
        flx_w2s(nvp+1:nv_adv,i,j)=flx_w2s(nvp+1:nv_adv,i,j)+ dporow*cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)
#endif
        cv_sed(0,k,i,j)=(porowater1*cv_sed1(0)+dporow*sal_bottom_MUSTANG(i,j))/porowater_new
        flx_w2s(0,i,j)=flx_w2s(0,i,j)+ dporow*sal_bottom_MUSTANG(i,j)
        cv_sed(-1,k,i,j)=(porowater1*cv_sed1(-1)+dporow*temp_bottom_MUSTANG(i,j))/porowater_new
        flx_w2s(-1,i,j)=flx_w2s(-1,i,j)+ dporow*temp_bottom_MUSTANG(i,j)
#endif
#if ! defined key_nofluxwat_IWS
        phieau_s2w(i,j)=phieau_s2w(i,j)-REAL(dporow*CELL_SURF(i,j),rlg)
#endif

      ELSE 
          ! output of porewater in sediment (cv_sed unchanged)
          
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
        flx_w2s(nvp+1:nv_adv,i,j)=flx_w2s(nvp+1:nv_adv,i,j)+ dporow*cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)
#endif
        flx_w2s(0,i,j)=flx_w2s(0,i,j)+ dporow*sal_bottom_MUSTANG(i,j)
        flx_w2s(-1,i,j)=flx_w2s(-1,i,j)+ dporow*temp_bottom_MUSTANG(i,j)
#endif
#if ! defined key_nofluxwat_IWS
        phieau_s2w(i,j)=phieau_s2w(i,j)-REAL(dporow*CELL_SURF(i,j),rlg)
#endif
      ENDIF 

    ELSE IF (code==6) THEN
      !! case where there is a deposit and a new surface layer is created
      !! l_increase_dep=False
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
        cv_sed(nvp+1:nv_adv,k,i,j)=cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)
        flx_w2s(nvp+1:nv_adv,i,j)=flx_w2s(nvp+1:nv_adv,i,j)+  &
             cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)*porowater_new
#endif
        cv_sed(0,k,i,j)=sal_bottom_MUSTANG(i,j)
        cv_sed(-1,k,i,j)=temp_bottom_MUSTANG(i,j)
        flx_w2s(0,i,j)=flx_w2s(0,i,j)+sal_bottom_MUSTANG(i,j)*porowater_new
        flx_w2s(-1,i,j)=flx_w2s(-1,i,j)+temp_bottom_MUSTANG(i,j)*porowater_new
#endif
#if ! defined key_nofluxwat_IWS
        phieau_s2w(i,j)=phieau_s2w(i,j)-REAL(porowater_new*CELL_SURF(i,j),rlg)
#endif

     !  if(cv_sed(nv_adv,k,i,j) < 9.999999999_rsh)then
     !     write(*,*)t,k,i,j,code,cw_bottom_MUSTANG(nv_adv,i,j),porowater_new
     !  endif

    ELSE IF (code==7) THEN
      !! case where there is a deposit and a new surface layer is created
      !!l_icrease_dep=true  
      ! 
      ! water quantity in the new layer coming from deposit  (porowater2)
      dporow=porowater1+porowater2-porowater_new
      IF( dporow > 0.0_rsh) THEN
          ! output of porewater from sediment to water

#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
        cv_sed(nvp+1:nv_adv,k,i,j)=(cv_sed1(nvp+1:nv_adv)*porowater1  &
                                   +cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)*porowater2) /(porowater_new+dporow)
        flx_w2s(nvp+1:nv_adv,i,j)=flx_w2s(nvp+1:nv_adv,i,j)- cv_sed(nvp+1:nv_adv,k,i,j)*dporow
#endif
        cv_sed(0,k,i,j)=(porowater1*cv_sed1(0)+porowater2*sal_bottom_MUSTANG(i,j))/(porowater_new+dporow)
        !flx_w2s(0,i,j)=flx_w2s(0,i,j)+ porowater2*sal_bottom_MUSTANG(i,j)
        flx_w2s(0,i,j)=flx_w2s(0,i,j)- cv_sed(0,k,i,j)*dporow
        cv_sed(-1,k,i,j)=(porowater1*cv_sed1(-1)+porowater2*temp_bottom_MUSTANG(i,j))/(porowater_new+dporow)
        !flx_w2s(-1,i,j)=flx_w2s(-1,i,j)+ porowater2*temp_bottom_MUSTANG(i,j)
        flx_w2s(-1,i,j)=flx_w2s(-1,i,j)+ - cv_sed(-1,k,i,j)*dporow
#endif
#if ! defined key_nofluxwat_IWS
        phieau_s2w(i,j)=phieau_s2w(i,j)+REAL(dporow*CELL_SURF(i,j),rlg)
#endif

      ELSE
          ! input of porewater in sediment  dporow <0
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
        cv_sed(nvp+1:nv_adv,k,i,j)=(cv_sed1(nvp+1:nv_adv)*porowater1+ &
                                    cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)*(porowater2-dporow)) /porowater_new
        flx_w2s(nvp+1:nv_adv,i,j)=flx_w2s(nvp+1:nv_adv,i,j)- cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)*(porowater2-dporow)
        !flx_w2s(nvp+1:nv_adv,i,j)=flx_w2s(nvp+1:nv_adv,i,j)- cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)*dporow
#endif
        cv_sed(-1,k,i,j)=(cv_sed1(-1)*porowater1+   &
                                   temp_bottom_MUSTANG(i,j)*(porowater2-dporow)) /porowater_new
        cv_sed(0,k,i,j)=(cv_sed1(0)*porowater1     &
                                   -sal_bottom_MUSTANG(i,j)*(porowater2-dporow)) /porowater_new
       ! flx_w2s(-1,i,j)=flx_w2s(-1,i,j)- temp_bottom_MUSTANG(i,j)*dporow
       ! flx_w2s(0,i,j)=flx_w2s(0,i,j)- sal_bottom_MUSTANG(i,j)*dporow
        flx_w2s(-1,i,j)=flx_w2s(-1,i,j)- temp_bottom_MUSTANG(i,j)*(porowater2-dporow)
        flx_w2s(0,i,j)=flx_w2s(0,i,j)- sal_bottom_MUSTANG(i,j)*(porowater2-dporow)
#endif
#if ! defined key_nofluxwat_IWS
         !phieau_s2w(i,j)=phieau_s2w(i,j)+REAL(dporow*CELL_SURF(i,j),rlg)
         phieau_s2w(i,j)=phieau_s2w(i,j)+REAL((porowater2-dporow)*CELL_SURF(i,j),rlg)
#endif
      ENDIF


    ELSE IF (code ==8 ) THEN
      !! case where 2 layers are melted at surface but not in the case of active layer
      dporow=porowater1+porowater2-porowater_new
      IF( dporow > 0.0_rsh) THEN
          ! output of porewater in sediment

#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
         cv_sed(nvp+1:nv_adv,k,i,j)=(cv_sed1(nvp+1:nv_adv)*porowater1+ &
                                    cv_sed2(nvp+1:nv_adv)*porowater2) &
                                      /(porowater_new+dporow)
         !flx_s2w(nvp+1:nv_adv,i,j)=flx_s2w(nvp+1:nv_adv,i,j)+  &
         flx_s2w_eroij(nvp+1:nv_adv)=flx_s2w_eroij(nvp+1:nv_adv)+  &
                                cv_sed(nvp+1:nv_adv,k,i,j)*dporow
#endif
         cv_sed(-1:0,k,i,j)=(cv_sed1(-1:0)*porowater1+ &
                                    cv_sed2(-1:0)*porowater2) &
                                      /(porowater_new+dporow)
         !flx_s2w(-1:0,i,j)=flx_s2w(-1:0,i,j)+ cv_sed(-1:0,k,i,j)*dporow
         flx_s2w_eroij(-1:0)=flx_s2w_eroij(-1:0)+ cv_sed(-1:0,k,i,j)*dporow
#endif
#if ! defined key_nofluxwat_IWS
         !phieau_s2w(i,j)=phieau_s2w(i,j)+REAL(dporow*CELL_SURF(i,j),rlg)
         phieau_ero_ij=phieau_ero_ij+REAL(dporow*CELL_SURF(i,j),rlg)
#endif
      ELSE
          ! input of porewater in sediment  dporow <0
#if ! defined key_noTSdiss_insed
#if ! defined key_Pconstitonly_insed
        cv_sed(nvp+1:nv_adv,k,i,j)=(cv_sed1(nvp+1:nv_adv)*porowater1+ &
                                    cv_sed2(nvp+1:nv_adv)*porowater2  &
                                   -cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)*dporow) /porowater_new
        !flx_s2w(nvp+1:nv_adv,i,j)=flx_s2w(nvp+1:nv_adv,i,j)-dporow*cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)
        flx_s2w_eroij(nvp+1:nv_adv)=flx_s2w_eroij(nvp+1:nv_adv)-dporow*cw_bottom_MUSTANG(nvp+1:nv_adv,i,j)

#endif
        cv_sed(-1,k,i,j)=(cv_sed1(-1)*porowater1+ cv_sed2(-1)*porowater2  &
                                   -temp_bottom_MUSTANG(i,j)*dporow) /porowater_new
        cv_sed(0,k,i,j)=(cv_sed1(0)*porowater1+ cv_sed2(0)*porowater2     &
                                   -sal_bottom_MUSTANG(i,j)*dporow) /porowater_new
        !flx_s2w(-1,i,j)=flx_s2w(-1,i,j)+ temp_bottom_MUSTANG(i,j)*dporow
        flx_s2w_eroij(-1)=flx_s2w_eroij(-1)+ temp_bottom_MUSTANG(i,j)*dporow
        !flx_s2w(0,i,j)=flx_s2w(0,i,j)+ sal_bottom_MUSTANG(i,j)*dporow
        flx_s2w_eroij(0)=flx_s2w_eroij(0)+ sal_bottom_MUSTANG(i,j)*dporow
#endif
#if ! defined key_nofluxwat_IWS
         !phieau_s2w(i,j)=phieau_s2w(i,j)+REAL(dporow*CELL_SURF(i,j),rlg)
         phieau_ero_ij=phieau_ero_ij+REAL(dporow*CELL_SURF(i,j),rlg)
#endif
      ENDIF

 
     ENDIF
  
END SUBROUTINE MUSTANGV2_eval_dissvar_IWSflux
#endif
 
   !!==============================================================================
 
! end MUSTANG_V2
#endif
                                
   !!===========================================================================
 

#if defined key_MUSTANG_bedload
  !!============================================================================== 
SUBROUTINE MUSTANGV2_eval_bedload(i, j, ksmax, flx_bxij, flx_byij) 

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE MUSTANGV2_eval_bedload  ***
   !&E
   !&E ** Purpose : Compute bedload transport from the formulation of Wu and Lin (2014)
   !&E              with hinding/exposure processes
   !&E
   !&E ** Description : 0D
   !&E                variables IN :  ksmax,CELL_DX,CELL_DY, ibedload1,ibedload2
   !&E                              diam_sed,cv_sed,c_sedtot,roswat_bot
   !&E                              stresscri0, tauskin,ros,raphbx,raphby,
   !&E                              tauskin_c_u,tauskin_c_v,
   !&E                              l_peph_bedload,l_slope_effect_bedload
   !&E                variables OUT : flx_bx,flx_by
   !&E
   !&E ** Note : NUMBER_PI and GRAVITY must be known as a parameters transmtted by coupleur 
   !&E           in MARS : coupleur_dimhydro.h (USE ..)
   !&E           in CROCO : module_MUSTANG.F (include..)
   !&E
   !&E ** Called by :  sed_erosion
   !&E
   !&E--------------------------------------------------------------------------


   !! * Modules used  
# if defined key_ANA_bedload || defined ANA_DUNE
#  include "ocean2d.h"
# endif

   !! * Arguments
   INTEGER,INTENT(IN)                                      :: i, j, ksmax
   REAL(KIND=rsh),DIMENSION(1:nvp),INTENT(out)             :: flx_bxij, flx_byij 

   !! * Local declarations
   INTEGER        ::  iv, jiv
   REAL(KIND=rsh) ::  ph, pe, pephm_fcor
   REAL(KIND=rsh) ::  phi_bed, qb
   REAL(KIND=rsh) ::  dhdx, dhdy, betas, betan, alphas, alphan, flx_bxij_star, flx_byij_star
   REAL(KIND=rsh),DIMENSION(nvpc) :: toce_loc
   REAL(KIND=rsh),PARAMETER :: m = 0.6_rsh ! for hinding / exposure processes
   REAL(KIND=rsh), PARAMETER :: tand30 = 0.577350269189626  ! TAND(30.0_rsh)=TAN(30*NUMBER_PI/180)=TAN(0.523599)
   !!---------------------------------------------------------------------------
   !! * Executable part

   !nu=0.00000136_rsh !a utiliser si capacite de transport pour suspension

    ! Bedload Flux  calculated in center of the mesh i,j (will intervene in the computation of "mass_out")
    ! see Rivier et al. (2017) Coastal Dynamics)   
    flx_bxij(:) = 0.0_rsh
    flx_byij(:) = 0.0_rsh
    ! Critical shear stress for each class based on masking / exposure processes   
    toce_loc(:) = 0.0_rsh

#ifdef key_MUSTANG_debug
    IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
        print *,'    EVAL BEDLOAD'
        print *,''
        print *,'      > Ini flx_bxij(:)=',flx_bxij(:)
        print *,'            flx_byij(:)=',flx_byij(:)
        print *,'            toce(:)=',toce_loc(:)
        print *,''
      END IF
#endif

   DO iv=ibedload1,ibedload2

     !!! Critical shear stress toce_loc in N/m2
     IF (l_peph_bedload) THEN
       ph=0.0_rsh
       pe=0.0_rsh
       ! Hindering / exposure coefficients for sediment class iv, ph and pe, respectively
       DO jiv=1,isand2 ! imud2 ou isand2 --> reflechir --> preference pour isand2
         ph=ph+(diam_sed(jiv)/(diam_sed(jiv)+diam_sed(iv)))  &
                    *( max(0.0_rsh,cv_sed(jiv,ksmax,i,j))/(max(epsilon_MUSTANG,c_sedtot(ksmax,i,j))) )
         pe=pe+(diam_sed(iv) /(diam_sed(jiv)+diam_sed(iv)))  &
                    *( max(0.0_rsh,cv_sed(jiv,ksmax,i,j))/(max(epsilon_MUSTANG,c_sedtot(ksmax,i,j))) )
       END DO
       pephm_fcor=(pe/ph)**(-m)
#if defined key_MUSTANG_specif_outputs        
       varspecif3Dnv_save(4,iv,i,j)=pephm_fcor
#endif
       toce_loc(iv)=stresscri0(iv)*pephm_fcor
       !IF ((ph .LT. 0.0_rsh) .OR. (pe.LT. 0.0_rsh) .OR. (ph .GT. 1.0_rsh) .OR. (pe .GT. 1.0_rsh)) &
                        ! write(*,*) 'peph', i,j,pe, ph
     ELSE
       toce_loc(iv)=stresscri0(iv)
     END IF

# if defined key_ANA_bedload || defined ANA_DUNE
     phi_bed=0.001*ubar(i,j,nrhs)**3.0_rsh
     qb=phi_bed & !  m2/s
                  *ros(iv)*cv_sed(iv,ksmax,i,j)/(c_sedtot(ksmax,i,j)+epsilon_MUSTANG) ! kg/m/s
# else
     phi_bed=0.0053_rsh*(max((tauskin(i,j)/toce_loc(iv))-1.0_rsh,0.0_rsh))**2.2_rsh

     ! Calculation of the rate of transport by bedload for class iv.
     ! We multiply by ros (iv) so we choose to have qb in kg/m/ s and not in m2/s     
     qb=phi_bed*sqrt((ros(iv)/RHOREF-1.0_rsh)*GRAVITY*diam_sed(iv)**3.0_rsh)  &
                    *ros(iv)*cv_sed(iv,ksmax,i,j)/(c_sedtot(ksmax,i,j)+epsilon_MUSTANG) ! kg/m/s
# endif
     !qb_ini(iv,i,j)=qb !pour ecriture en sortie

     !============Projection sur x et y en fonction de la direction de la tension sur le fond ==============

# if defined key_ANA_bedload || defined ANA_DUNE
     flx_byij(iv)=0.
#endif

     flx_bxij(iv) = qb * tauskin_x(i, j) / (tauskin_c(i, j) + epsilon_MUSTANG)
     flx_byij(iv) = qb * tauskin_y(i, j) / (tauskin_c(i, j) + epsilon_MUSTANG)


#ifdef key_MUSTANG_debug
    IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
          print *,'      iv=',iv
          IF (l_peph_bedload) print *,'        - pe/ph = ',pe,' / ',ph
          print *,'        - toce(iv)=',toce_loc(iv)
          print *,'        - phi_bed=',phi_bed
          print *,'        - qb=',qb
          print *,'        - flx_bxij(iv)=',flx_bxij(iv),' in kg/m/s'
          print *,'        - flx_byij(iv)=',flx_byij(iv)
        END IF
#endif


     IF (l_slope_effect_bedload) THEN
            
           !============Accounting for slope effects in computations of bedload fluxes============

           ! parametres Lesser : alphabs and alphabn are defined in paraMUSTANGV2.txt
           !============Preparation parametres Lesser===
           IF(htot(i+1,j) .GT. h0fond .AND. htot(i-1,j) .GT. h0fond ) THEN
              dhdx=slope_dhdx(i,j)
           ELSE
              dhdx=0.0_rsh
           ENDIF
           IF(htot(i,j+1) .GT. h0fond .AND. htot(i,j-1).GT. h0fond ) THEN
              dhdy=slope_dhdy(i,j)
           ELSE
              dhdy=0.0_rsh
           ENDIF
           betas=MIN(ATAN(dhdx*flx_bxij(iv)/(qb+epsilon_MUSTANG)           &
                 +dhdy*flx_byij(iv)/(qb+epsilon_MUSTANG))          &
                               ,0.9_rsh*30.0_rsh/180.0_rsh*NUMBER_PI)
           betan=ATAN(-dhdx*flx_byij(iv)/(qb+epsilon_MUSTANG)+dhdy*flx_bxij(iv)/(qb+epsilon_MUSTANG))

           alphas=1.0_rsh+alphabs*(tand30/(COS(betas)*(tand30-TAN(betas))+epsilon_MUSTANG)-1)
           alphan=alphabn*(toce_loc(iv)/(tauskin(i,j)+epsilon_MUSTANG))**0.5_rsh*TAN(betan)

           !============Application Lesser==============
           flx_bxij_star= alphas*flx_bxij(iv)
           flx_byij_star= alphas*flx_byij(iv)

           flx_bxij(iv)=flx_bxij_star-alphan*flx_byij_star
           flx_byij(iv)=flx_byij_star+alphan*flx_bxij_star

           !!==============================================================================
     ENDIF                                   

     ! sedimask_h0plusxe : = 1 si BATHY_H0(i,j)+WATER_ELEVATION(i,j) .GT. 1    = 0 sinon
     ! ==> Le flux charrie en X et Y est mis a 0 si la maille voisine est a terre
     !TODO : put this in subroutine in sed_MUSTANG_HOST because it is host dependant (raphbx&raphby)

     flx_bxij(iv) = (flx_bxij(iv) + abs(flx_bxij(iv))) * 0.5_rsh * raphbx(i+1, j) * sedimask_h0plusxe(i+1, j)+ &
                    (flx_bxij(iv) - abs(flx_bxij(iv))) * 0.5_rsh * raphbx(  i, j) * sedimask_h0plusxe(i-1, j)

     flx_byij(iv) = (flx_byij(iv) + abs(flx_byij(iv))) * 0.5_rsh * raphby(i, j+1) * sedimask_h0plusxe(i, j+1)+ &
                    (flx_byij(iv) - abs(flx_byij(iv))) * 0.5_rsh * raphby(i,   j) * sedimask_h0plusxe(i, j-1)

    ! Morphological factor

     flx_bxij(iv) = MF * fwet(i, j) * flx_bxij(iv)
     flx_byij(iv) = MF * fwet(i, j) * flx_byij(iv)

#if defined key_MUSTANG_specif_outputs        
     ! flx_bx_int and flx_by_int
     varspecif2D_save(16,i,j)=varspecif2D_save(16,i,j)+flx_bxij(iv) !pour ecriture en sortie
     varspecif2D_save(17,i,j)=varspecif2D_save(17,i,j)+flx_byij(iv) !pour ecriture en sortie
#endif

#ifdef key_MUSTANG_debug
        IF ( l_debug_erosion .AND. i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug .AND. CURRENT_TIME> t_start_debug) THEN
          print *,'      apres application des masques'
          print *,'        - flx_bxij(iv)=',flx_bxij(iv),' in kg/m/s'
          print *,'        - flx_byij(iv)=',flx_byij(iv)
        END IF
#endif
     ! So we have at the output the bedload fluxes coming out of the mesh i,j en kg/m/s


   END DO


END SUBROUTINE MUSTANGV2_eval_bedload

! end key_MUSTANG_bedload
#endif
!!==============================================================================



!!===========================================================================
real function MUSTANG_E0sand(diamsan, taucr, rossan, ws_sand)
!&E--------------------------------------------------------------------------
!&E                 *** FUNCTION MUSTANG_E0sand  ***
!&E Compute the sand erosion constant in kg.m-2.s-1 given e0_sand_option. 
!&E * If e0_sand_option = 0 : use of  e0_sand read in namelist
!&E * If e0_sand_option = 1 : sand erosion constant computed from Van Rijn (1984) 
!&E                           formulation
!&E * If e0_sand_option = 2 : sand erosion constant computed by an erodimetry 
!&E                           formulation (function of mean sand diameter)
!&E * If e0_sand_option = 3 : sand erosion constant computed in order to use a 
!&E                           pick-up function for the erosion fluxe. The volumetric 
!&E                           concentration at reference level aref is from  
!&E                           Wu and Lin (2014) eq.34.
!&E--------------------------------------------------------------------------
implicit none
real, intent(in) :: diamsan ! sediment diameter (m)
real, intent(in) :: taucr   ! critical shear stress (N.m-2)
real, intent(in) :: rossan  ! sediment density (kg/m3)
real, intent(in) :: ws_sand ! sediment settling velocity (m/s)
! SOURCE
real             :: diamsan_star, ws_diamsan

if (E0_sand_option==0) then 
  ! use of E0_sand read in namelist (E0_sand_para read in paraMUSTANG)
  MUSTANG_E0sand = E0_sand_Cst
elseif (E0_sand_option == 1 .and. rossan .gt. 1000_rsh) then
  ! E0_sand computed from Van Rijn (1984) formulation
  MUSTANG_E0sand = 0.00033_rsh * rossan*((rossan / 1000.0_rsh - 1.0_rsh) * 9.81 * diamsan)**0.5_rsh &
    * (diamsan * ((rossan / 1000.0_rsh - 1) * 9.81_rsh/(0.000001_rsh)**2)**(1.0_rsh/3.0_rsh))**(0.3_rsh)
elseif (E0_sand_option == 2) then
  ! E0_sand computed by an erodimetry formulation (function of mean sand diameter)
  MUSTANG_E0sand = taucr**n_eros_sand * min(0.27, 1000. * diamsan - 0.01)
elseif (E0_sand_option == 3) then
  ! E0_sand computed from WU and Lin (2014) concentration at reference level and Van Rijn pick-up function
  MUSTANG_E0sand = ws_sand * rossan * 0.0032 * diamsan / aref_sand 
endif

MUSTANG_E0sand = E0_sand_para * MUSTANG_E0sand

end function ! MUSTANG_E0sand

#if defined key_MUSTANG_V2
!!==============================================================================
   FUNCTION isitcohesive(cvsed_sup_i_j, criterion_cohesive)
   !&E--------------------------------------------------------------------------
   !&E                 *** FUNCTION isitcohesive  ***
   !&E      returns true if the superficial sediment is cohesive, false otherwise
   !&E--------------------------------------------------------------------------

      LOGICAL :: isitcohesive
      REAL(KIND=rsh), DIMENSION(-1:nv_tot):: cvsed_sup_i_j
      REAL(KIND=rsh) :: sommud, somsan, somgrav, criterion_cohesive
      INTEGER :: iv
      REAL(KIND=rsh), PARAMETER  ::  epsilon_isit = 0.000000001_rsh

      sommud = 0.0_rsh
      somsan = 0.0_rsh
      somgrav = 0.0_rsh

      DO iv = imud1, imud2
        sommud = sommud + cvsed_sup_i_j(iv)
      ENDDO
      DO iv = isand1, isand2
        somsan = somsan + cvsed_sup_i_j(iv)
      ENDDO
      DO iv = igrav1,igrav2
        somgrav = somgrav + cvsed_sup_i_j(iv)
      ENDDO

      IF (sommud / (sommud + somsan + somgrav + epsilon_isit) .GE. criterion_cohesive) THEN
       isitcohesive = .TRUE.
      ELSE
       isitcohesive = .FALSE.
      END IF

   END FUNCTION isitcohesive
!!==============================================================================
   FUNCTION dzsminvar(frac_sediv)
   !&E--------------------------------------------------------------------------
   !&E                 *** FUNCTION dzsminvar  ***
   !&E           returns dzsmin from sediment composition
   !&E--------------------------------------------------------------------------

      REAL(KIND=rsh), DIMENSION(1:nvpc)    :: frac_sediv
      REAL(KIND=rsh)                      :: dzsminvar

       dzsminvar = (1.0_rsh - coeff_dzsmin) * dzsminuni + &
                  coeff_dzsmin * SUM( frac_sediv(1:nvpc) * diam_sed(1:nvpc) )

   END FUNCTION dzsminvar
!!==============================================================================
#endif /* end if define key_MUSTANG_V2 */
#endif  /* end if define MUSTANG */

END MODULE sed_MUSTANG

