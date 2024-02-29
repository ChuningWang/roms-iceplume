#include "cppdefs.h"

MODULE sedinitrc
   !!======================================================================
   !!              ***  MODULE  sedinitrc  ***
   !! Sediment : define sediment variables
   !!=====================================================================
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   sed_init    : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   !! * Modules used
   USE sed     ! sediment global variable
   USE sed_oce
   USE sedini
   USE seddta
   USE sedrst
   USE sedco3
   USE sedchem
   USE sedarr

   IMPLICIT NONE
   PRIVATE

   !!* Substitution
#  include "ocean2pisces.h90"

   REAL(wp)    ::  &
      ryear = 365. * 24. * 3600. !:  1 year converted in second

   !! *  Routine accessibility
   PUBLIC sed_initrc          ! routine called by opa.F90

   !! $Id: sedini.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
CONTAINS


   SUBROUTINE sed_initrc
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_init  ***
      !!
      !! ** Purpose :  Initialization of sediment module
      !!               - Reading namelist
      !!               - Read the deepest water layer thickness
      !!                 ( using as mask ) in Netcdf file
      !!               - Convert unity if necessary
      !!               - sets initial sediment composition
      !!                 ( only clay or reading restart file )
      !!               - sets sediment grid, porosity and others constants
      !!
      !!   History :
      !!        !  04-10  (N. Emprin, M. Gehlen )  Original code
      !!        !  06-07  (C. Ethe)  Re-organization
      !!----------------------------------------------------------------------
      INTEGER :: ji, jj, ikt
      !!----------------------------------------------------------------------


      ! Initialize the sediment tracers concentrations
      !------------------------------------------------

      IF(lwp) WRITE(numsed,*) ' sed_initrc : Initialization of sediment concentration '
      IF(lwp) WRITE(numsed,*) ' '

      ! Determination of sediments number of points and allocate global variables

      ! sets initial sediment composition
      ! ( only clay or reading restart file )
      !---------------------------------------
      CALL sed_init_data


      CALL sed_init_wri


   END SUBROUTINE sed_initrc


   SUBROUTINE sed_init_data
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_init_data  ***
      !!
      !! ** Purpose :  Initialization of sediment module
      !!               - sets initial sediment composition
      !!                 ( only clay or reading restart file )
      !!
      !!   History :
      !!        !  06-07  (C. Ethe)  original
      !!----------------------------------------------------------------------
 
      ! local variables
      INTEGER :: &
         ji, jk, zhipor

      !--------------------------------------------------------------------
 

      IF( .NOT. ln_rst_sed ) THEN

         IF (lwp) WRITE(numsed,*) ' Initilization of default values of sediment components'

         ! default values for initial pore water concentrations [mol/l]
         pwcp(:,:,:) = 0.
         ! default value for initial solid component (fraction of dry weight dim=[0])
         ! clay
         solcp(:,:,:) = 0.
         solcp(:,2:jpksed,jsclay) = 1.0 * 0.965
         solcp(:,2:jpksed,jsfeo)  = 1.0 * 0.035

         ! Initialization of [h+] and [co3--]

         zhipor = 8.0
         ! Initialization of [h+] in mol/kg
         DO jk = 1, jpksed
            DO ji = 1, jpoce
               hipor (ji,jk) = 10.**( -1. * zhipor )
            ENDDO
         ENDDO

         co3por(:,:) = 1E-6

      ELSE   
  
         IF (lwp) WRITE(numsed,*) ' Initilization of Sediment components from restart'
         CALL sed_rst_read

      ENDIF


      ! Load initial Pisces Data for bot. wat. Chem and fluxes
      CALL sed_dta ( nitsed000 ) 

      ! Initialization of chemical constants
      CALL sed_chem ( nitsed000 )

      ! Stores initial sediment data for mass balance calculation
      pwcp0 (1:jpoce,1:jpksed,1:jpwat ) = pwcp (1:jpoce,1:jpksed,1:jpwat ) 
      solcp0(1:jpoce,1:jpksed,1:jpsol ) = solcp(1:jpoce,1:jpksed,1:jpsol) 

      ! Conversion of [h+] in mol/Kg to get it in mol/l ( multiplication by density)
      DO jk = 1, jpksed
         hipor(1:jpoce,jk) = hipor(1:jpoce,jk) * densSW(1:jpoce)
      ENDDO


      ! In default case - no restart - sedco3 is run to initiate [h+] and [co32-]
      ! Otherwise initiate values of pH and co3 read in restart
      IF( .NOT. ln_rst_sed ) THEN
         ! sedco3 is run to initiate[h+] [co32-] in mol/l of solution
         CALL sed_co3 ( nitsed000 )

      ENDIF
            
   END SUBROUTINE sed_init_data

   SUBROUTINE sed_init_wri

      INTEGER :: jpij, jk

      IF (lwp) THEN
         jpij = jpi*jpj
         WRITE(numsed,*)' '
         WRITE(numsed,*)'======== Write summary of sediment char.  ============'
         WRITE(numsed,*)' '
         WRITE(numsed,*)' '
         WRITE(numsed,*)'-------------------------------------------------------------------'
         WRITE(numsed,*)' Initial Conditions '
         WRITE(numsed,*)'-------------------------------------------------------------------'
         WRITE(numsed,*)'dzm = dzkbot minimum to calculate ', 0.
         WRITE(numsed,*)'Local zone : jpi, jpj, jpksed : ',jpi, jpj, jpksed
         WRITE(numsed,*)'jpoce = ',jpoce,' nbtot pts = ',jpij,' nb earth pts = ',jpij - jpoce
         WRITE(numsed,*)'sublayer thickness dz(1) [cm] : ', dz(1)
         WRITE(numsed,*)'Vertical domain of the sediment'
         WRITE(numsed,*)'-------------------------------'
         WRITE(numsed,*)' Indice, profsed, dz'
         DO jk = 2, jpksed
            WRITE(numsed,*) jk,profsed(jk),dz(jk) 
         END DO
         WRITE(numsed,*)' nb solid comp : ',jpsol
         WRITE(numsed,*)'(1=opal,2=clay,3=POC,4=CaCO3), 5=POS, 6=POR, 7=FEO, 8=FeS'
         WRITE(numsed,*)'weight mol 1,2,3,4,5,6,7'
         WRITE(numsed,'(8(F0.2,3X))')mol_wgt(jsopal),mol_wgt(jsclay),mol_wgt(jspoc),mol_wgt(jscal),   &
         &                           mol_wgt(jspos),mol_wgt(jspor),mol_wgt(jsfeo),mol_wgt(jsfes)
         WRITE(numsed,*)'nb dissolved comp',jpwat
         WRITE(numsed,*)'1=silicic acid,,2=O2,3=DIC,4=NO3,5=PO4,6=Alk,7=NH4,8=ODU'
         WRITE(numsed,*)'redfield coef C,O,N P Dit '
         WRITE(numsed,'(5(F0.2,3X))')1./spo4r,so2ut/spo4r,srno3/spo4r,spo4r/spo4r,srDnit/spo4r
         WRITE(numsed,*) ' '
         WRITE(numsed,*) ' End Of Initialization '
         WRITE(numsed,*) ' '
      ENDIF
!
   END SUBROUTINE sed_init_wri

#endif

END MODULE sedinitrc
