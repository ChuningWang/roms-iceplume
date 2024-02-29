#include "cppdefs.h"

MODULE par_pisces
   !!======================================================================
   !!                        ***  par_pisces  ***
   !! TOP :   set the PISCES parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   IMPLICIT NONE
   PUBLIC

#if defined key_pisces
   !!---------------------------------------------------------------------
   !!   'key_pisces'   :                         standard PISCES bio-model
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_pisces     = .TRUE.  !: PISCES flag 
#if defined key_pisces_quota
#   if defined key_ligand
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 40      !: number of PISCES passive tracers
#   else
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 39      !: number of PISCES passive tracers
#   endif
#else
#   if defined key_ligand
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 25      !: number of PISCES passive tracers
#   else
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 24      !: number of PISCES passive tracers
#   endif
#endif
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_2d  = 11      !: additional 2d output ('key_trc_diaadd')
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_3d  = 16      !: additional 3d output ('key_trc_diaadd')

   ! assign an index in trc arrays for each LOBSTER prognostic variables
   !    WARNING: be carefull about the order when reading the restart
        !   !!gm  this warning should be obsolet with IOM
   INTEGER, PUBLIC, PARAMETER ::   jpdic =  1    !: dissolved inoganic carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jptal =  2    !: total alkalinity 
   INTEGER, PUBLIC, PARAMETER ::   jpoxy =  3    !: oxygen carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpcal =  4    !: calcite  concentration 
   INTEGER, PUBLIC, PARAMETER ::   jppo4 =  5    !: phosphate concentration 
   INTEGER, PUBLIC, PARAMETER ::   jppoc =  6    !: small particulate organic phosphate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpsil =  7    !: silicate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpphy =  8    !: phytoplancton concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpzoo =  9    !: zooplancton concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdoc = 10    !: dissolved organic carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpdia = 11    !: Diatoms Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpmes = 12    !: Mesozooplankton Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpgsi = 13    !: (big) Silicate Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpfer = 14    !: Iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpbfe = 15    !: Big iron particles Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpgoc = 16    !: big particulate organic phosphate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpsfe = 17    !: Small iron particles Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdfe = 18    !: Diatoms iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdsi = 19    !: Diatoms Silicate Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnfe = 20    !: Nano iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnch = 21    !: Nano Chlorophyll Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdch = 22    !: Diatoms Chlorophyll Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpno3 = 23    !: Nitrates Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnh4 = 24    !: Ammonium Concentration
#if defined key_ligand
   INTEGER, PUBLIC, PARAMETER ::   jplgw = 25    !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpp4z = 25    !: Number of P4Z tracers
#else
   INTEGER, PUBLIC ::   jplgw    !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpp4z = 24    !: Number of p4z tracers
#endif
#if defined key_pisces_quota
   INTEGER, PUBLIC, PARAMETER ::   jpdon = jpp4z + 1   !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdop = jpp4z + 2   !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jppon = jpp4z + 3    !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jppop = jpp4z + 4   !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnph = jpp4z + 5   !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jppph = jpp4z + 6   !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpndi = jpp4z + 7   !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jppdi = jpp4z + 8   !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jppic = jpp4z + 9   !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnpi = jpp4z + 10  !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpppi = jpp4z + 11  !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jppfe = jpp4z + 12  !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jppch = jpp4z + 13  !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpgon = jpp4z + 14  !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpgop = jpp4z + 15  !: Ammonium Concentration
#else
   INTEGER, PUBLIC ::   jpdon    !: DON concentration 
   INTEGER, PUBLIC ::   jpdop    !: DOP concentration 
   INTEGER, PUBLIC ::   jppon    !: PON concentration
   INTEGER, PUBLIC ::   jppop    !: POP concentration
   INTEGER, PUBLIC ::   jpnph     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppph     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpndi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppdi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppic     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpnpi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpppi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppfe     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppch     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpgon    !: GON concentration
   INTEGER, PUBLIC ::   jpgop    !: GOP concentration
#endif

   INTEGER, PUBLIC ::   jp_flxco2  
   INTEGER, PUBLIC ::   jp_flxo2   
   INTEGER, PUBLIC ::   jp_kgco2   
   INTEGER, PUBLIC ::   jp_dpco2   
   INTEGER, PUBLIC ::   jp_sinkco2 
   INTEGER, PUBLIC ::   jp_sinkfer 
   INTEGER, PUBLIC ::   jp_sinksil 
   INTEGER, PUBLIC ::   jp_sinkcal 
   INTEGER, PUBLIC ::   jp_heup    
   INTEGER, PUBLIC ::   jp_sildep   
   INTEGER, PUBLIC ::   jp_po4dep   
   INTEGER, PUBLIC ::   jp_no3dep   
   INTEGER, PUBLIC ::   jp_nh4dep   
   INTEGER, PUBLIC ::   jp_nitrpot 

   INTEGER, PUBLIC ::   jp_hi      
   INTEGER, PUBLIC ::   jp_co3     
   INTEGER, PUBLIC ::   jp_co3sat  
   INTEGER, PUBLIC ::   jp_etot    
   INTEGER, PUBLIC ::   jp_pphy    
   INTEGER, PUBLIC ::   jp_pphy2   
   INTEGER, PUBLIC ::   jp_pnew    
   INTEGER, PUBLIC ::   jp_pnew2   
   INTEGER, PUBLIC ::   jp_pbsi    
   INTEGER, PUBLIC ::   jp_pfed    
   INTEGER, PUBLIC ::   jp_pfen    
   INTEGER, PUBLIC ::   jp_pnewo2  
   INTEGER, PUBLIC ::   jp_prego2  
   INTEGER, PUBLIC ::   jp_grapoc   
   INTEGER, PUBLIC ::   jp_grapoc2   
   INTEGER, PUBLIC ::   jp_mico2  
   INTEGER, PUBLIC ::   jp_meso2  
   INTEGER, PUBLIC ::   jp_nitrifo2 
   INTEGER, PUBLIC ::   jp_remino2 
   INTEGER, PUBLIC ::   jp_nfixo2  
   INTEGER, PUBLIC ::   jp_irondep  
   INTEGER, PUBLIC ::   jp_ironsed  
#else
   !!---------------------------------------------------------------------
   !!   Default                                   No CFC geochemical model
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_pisces     = .FALSE.  !: CFC flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     =  0       !: No CFC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_2d  =  0       !: No CFC additional 2d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_3d  =  0       !: No CFC additional 3d output arrays 
#endif

   ! Starting/ending PISCES do-loop indices (N.B. no PISCES : jpl_pcs < jpf_pcs the do-loop are never done)
   INTEGER, PUBLIC, PARAMETER ::   jptra       = jp_pisces                  !: First index of PISCES tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0     = 1                  !: First index of PISCES tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1     = jp_pisces          !: Last  index of PISCES tracers

   !!======================================================================
END MODULE par_pisces
