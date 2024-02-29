#include "cppdefs.h"

MODULE oce_sed
   !!======================================================================
   !!                        ***  sed  ***
   !! Sediment :   set sediment global variables
   !!======================================================================
   !! History :
   !!        !  06-12  (C. Ethe)  Orignal
   !!----------------------------------------------------------------------
#if defined key_pisces
   USE par_sed

   USE sms_pisces, ONLY : tmask     =>   tmask           !: Mask
   USE sms_pisces, ONLY : wsbio4    =>   wsbio4          !: sinking flux for POC
   USE sms_pisces, ONLY : wsbio3    =>   wsbio3          !: sinking flux for GOC
   USE sms_pisces, ONLY : wsbio2    =>   wsbio2           !: sinking flux for calcite
   USE sms_pisces, ONLY : wsbio     =>   wsbio           !: sinking flux for calcite
   USE sms_pisces, ONLY : ln_p5z    =>   ln_p5z          !: PISCES-QUOTA flag
   USE p4zche, ONLY     : akb3      =>   akb3            !: Chemical constants  
   USE sms_pisces, ONLY : ak13      =>   ak13            !: Chemical constants  
   USE sms_pisces, ONLY : ak23      =>   ak23            !: Chemical constants  
   USE p4zche, ONLY     : akw3      =>   akw3            !: Chemical constants  
   USE sms_pisces, ONLY : aksp      =>   aksp            !: Chemical constants  
   USE p4zche, ONLY     : borat     =>   borat           !: Chemical constants ( borat ) 
   USE p4zche, ONLY     : ak1p3     =>   ak1p3           !: Chemical constants  
   USE p4zche, ONLY     : ak2p3     =>   ak2p3           !: Chemical constants  
   USE p4zche, ONLY     : ak3p3     =>   ak3p3           !: Chemical constants  
   USE p4zche, ONLY     : aksi3     =>   aksi3           !: Chemical constants  
   USE p4zche, ONLY     : aks3      =>   aks3            !: Chemical constants  
   USE p4zche, ONLY     : akf3      =>   akf3            !: Chemical constants  
   USE p4zche, ONLY     : fluorid   =>   fluorid         !: Chemical constants  
   USE p4zche, ONLY     : sulfat    =>   sulfat          !: Chemical constants  
   USE p4zche, ONLY     : sio3eq    =>   sio3eq          !: Chemical constants  
   USE p4zsbc, ONLY     : dust      =>   dust

#endif

END MODULE oce_sed


