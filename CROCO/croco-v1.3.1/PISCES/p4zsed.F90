#include "cppdefs.h"

MODULE p4zsed
   !!======================================================================
   !!                         ***  MODULE p4sed  ***
   !! TOP :   PISCES Compute loss of organic matter in the sediments
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12 (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06 (C. Ethe) USE of fldread
   !!             3.5  !  2012-07 (O. Aumont) improvment of river input of nutrients 
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!   p4z_sed        :  Compute loss of organic matter in the sediments
   !!----------------------------------------------------------------------
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zlim          !  Co-limitations of differents nutrients
   USE p4zsbc          !  External source of nutrients 
   USE p4zint          !  interpolation and computation of various fields
   USE sed             !  Sediment module
!   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sed  
   PUBLIC   p4z_sed_alloc

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"
 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: nitrpot    !: Nitrogen fixation 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:  ) :: sdenit     !: Nitrate reduction in the sediments
   REAL(wp) :: r1_rday                  !: inverse of rday
   LOGICAL, SAVE :: lk_sed

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zsed.F90 10780 2019-03-20 17:53:44Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_sed( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed  ***
      !!
      !! ** Purpose :   Compute loss of organic matter in the sediments. This
      !!              is by no way a sediment model. The loss is simply 
      !!              computed to balance the inout from rivers and dust
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  zrivalk, zrivsil, zrivno3
      REAL(wp) ::  zwflux, zlim, zfact, zfactcal
      REAL(wp) ::  zo2, zno3, zflx, zpdenit, z1pdenit, zolimit
      REAL(wp) ::  zsiloss, zcaloss, zws3, zws4, zwsc, zdep
      REAL(wp) ::  zwstpoc, zwstpon, zwstpop, zmsk
      REAL(wp) ::  ztrfer, ztrpo4s, ztrdp, zwdust, zmudia, ztemp
      REAL(wp) ::  xdiano3, xdianh4, zrfact2
      !
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zsoufer, zlight
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY) :: zdenit2d, zbureff, zwork
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY) :: zwsbio3, zwsbio4
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY) :: zsedcal, zsedsi, zsedc
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zno3dep, znh4dep
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ztrpo4, ztrdop, zirondep, zpdep
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zsidep
      !!---------------------------------------------------------------------
      !
      IF( kt == nittrc000 .AND. knt == 1 )   THEN
          r1_rday  = 1. / rday
          IF (ln_sediment .AND. ln_sed_2way) THEN
             lk_sed = .TRUE.
          ELSE
             lk_sed = .FALSE.
          ENDIF
      ENDIF
      !
      ! Allocate temporary workspace
      ALLOCATE( ztrpo4(PRIV_3D_BIOARRAY) )
      IF( ln_p5z )    ALLOCATE( ztrdop(PRIV_3D_BIOARRAY) )

      zdenit2d(:,:) = 0.e0
      zbureff (:,:) = 0.e0
      zwork   (:,:) = 0.e0
      zsedsi  (:,:) = 0.e0
      zsedcal (:,:) = 0.e0
      zsedc   (:,:) = 0.e0

      ! Add the external input of nutrients from dust deposition
      ! ----------------------------------------------------------
      IF( ln_dust ) THEN
         !                                              
         ALLOCATE( zsidep(PRIV_2D_BIOARRAY), zpdep(PRIV_3D_BIOARRAY), zirondep(PRIV_3D_BIOARRAY) )
         !                                              ! Iron and Si deposition at the surface
         DO jj = JRANGE
            DO ji = IRANGE
               zirondep(ji,jj,1) = dustsolub  * dust(ji,jj) * mfrac * rfact2 / e3t_n(ji,jj,KSURF) / 55.85
               zsidep(ji,jj)   = 8.8 * 0.075 * dust(ji,jj) * mfrac * rfact2 / e3t_n(ji,jj,KSURF) / 28.1 
               zpdep (ji,jj,1) = 0.1 * 0.021 * dust(ji,jj) * mfrac * rfact2 / e3t_n(ji,jj,KSURF) / 31. / po4r 
            END DO
         END DO
         !                                              ! Iron solubilization of particles in the water column
         !                                              ! dust in kg/m2/s ---> 1/55.85 to put in mol/Fe ;  wdust in m/j
         zwdust = 0.03 * rday / ( wdust * 55.85 ) / ( 270. * rday )
         DO jk = KRANGEL
            DO jj = JRANGE
               DO ji = IRANGE 
                  zirondep(ji,jj,jk) = dust(ji,jj) * mfrac * zwdust * rfact2   &
                     &               * EXP( -gdept_n(ji,jj,K) / 540. )
                  zpdep   (ji,jj,jk) = zirondep(ji,jj,jk) * 0.023
               END DO
            END DO
         END DO
         !                                              ! Iron solubilization of particles in the water column
         tra(:,:,1,jpsil) = tra(:,:,1,jpsil) + zsidep  (:,:)
         DO jk = KRANGE
            tra(:,:,jk,jppo4) = tra(:,:,jk,jppo4) + zpdep   (:,:,jk)
            tra(:,:,jk,jpfer) = tra(:,:,jk,jpfer) + zirondep(:,:,jk) 
         ENDDO
         ! 
#if defined key_iomput
         IF( lk_iomput ) THEN
            IF( knt == nrdttrc ) THEN
                IF( iom_use( "Irondep" ) )   &
                &  CALL iom_put( "Irondep", zirondep(:,:,1) * 1.e+3 * rfact2r * e3t_n(:,:,1) * tmask(:,:,1) ) ! surface downward dust depo of iron
                IF( iom_use( "pdust" ) )   &
                &  CALL iom_put( "pdust"  , dust(:,:) / ( wdust * rday )  * tmask(:,:,1) ) ! dust concentration at surface
            ENDIF
         ENDIF
#endif
         !
#if defined key_trc_diaadd
         zrfact2 = 1.e+3 * rfact2r
         DO jj = JRANGE
            DO ji = IRANGE
               zmsk = zrfact2 * tmask(ji,jj,1)
               trc2d(ji,jj,jp_sildep)   = zsidep(ji,jj)  * zmsk                ! iron deposition
               trc2d(ji,jj,jp_po4dep)   = zpdep(ji,jj,1) * po4r * zmsk                ! iron deposition
            END DO
         END DO

         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  zmsk = zrfact2 * tmask(ji,jj,jk)
                  trc3d(ji,jj,K,jp_irondep)  = zirondep(ji,jj,jk) * zmsk                ! iron flux from dust
               END DO
            END DO
         ENDDO
         DEALLOCATE( zsidep, zpdep, zirondep )
# endif
         !                                              
      ENDIF
     
      ! Add the external input of nutrients from river
      ! ----------------------------------------------------------
      IF( ln_river ) THEN
         DO jj = JRANGE
            DO ji = IRANGE
               tra(ji,jj,1,jppo4) = tra(ji,jj,1,jppo4) +  rivdip(ji,jj) * rfact2
               tra(ji,jj,1,jpno3) = tra(ji,jj,1,jpno3) +  rivdin(ji,jj) * rfact2
               tra(ji,jj,1,jpfer) = tra(ji,jj,1,jpfer) +  rivdic(ji,jj) * 5.e-5 * rfact2
               tra(ji,jj,1,jpsil) = tra(ji,jj,1,jpsil) +  rivdsi(ji,jj) * rfact2
               tra(ji,jj,1,jpdic) = tra(ji,jj,1,jpdic) +  rivdic(ji,jj) * rfact2
               tra(ji,jj,1,jptal) = tra(ji,jj,1,jptal) +  ( rivalk(ji,jj) - rno3 * rivdin(ji,jj) ) * rfact2
               tra(ji,jj,1,jpdoc) = tra(ji,jj,1,jpdoc) +  rivdoc(ji,jj) * rfact2
            ENDDO
         ENDDO
         IF (ln_ligand) THEN
            DO jj = JRANGE
               DO ji = IRANGE
                  tra(ji,jj,1,jplgw) = tra(ji,jj,1,jplgw) +  rivdic(ji,jj) * 5.e-5 * rfact2
               ENDDO
            ENDDO
         ENDIF
         IF( ln_p5z ) THEN
            DO jj = JRANGE
               DO ji = IRANGE
                  tra(ji,jj,1,jpdop) = tra(ji,jj,1,jpdop) + rivdop(ji,jj) * rfact2
                  tra(ji,jj,1,jpdon) = tra(ji,jj,1,jpdon) + rivdon(ji,jj) * rfact2
               ENDDO
            ENDDO
         ENDIF
      ENDIF
      
      ! Add the external input of nutrients from nitrogen deposition
      ! ----------------------------------------------------------
      IF( ln_ndepo ) THEN
         ALLOCATE( zno3dep(PRIV_2D_BIOARRAY), znh4dep(PRIV_2D_BIOARRAY) )
         DO jj = JRANGE
            DO ji = IRANGE
               ! conversion from KgN/m2/s to molC/L/s
               zfact = rfact2 / rno3 / 14. / e3t_n(ji,jj,KSURF)
               zno3dep(ji,jj) =  zfact * no3dep(ji,jj)
               znh4dep(ji,jj) =  zfact * nh4dep(ji,jj)
               !
               tra(ji,jj,1,jpno3) = tra(ji,jj,1,jpno3) + zno3dep(ji,jj)
               tra(ji,jj,1,jpnh4) = tra(ji,jj,1,jpnh4) + znh4dep(ji,jj)
               tra(ji,jj,1,jptal) = tra(ji,jj,1,jptal) + rno3 * ( znh4dep(ji,jj) - zno3dep(ji,jj) )
            END DO
         END DO
#if defined key_trc_diaadd
         zrfact2 = 1.e+3 * rfact2r
         DO jj = JRANGE
            DO ji = IRANGE
               zmsk = zrfact2 * tmask(ji,jj,1)
               trc2d(ji,jj,jp_no3dep )  = zno3dep(ji,jj) * rno3 * zmsk                ! iron deposition
               trc2d(ji,jj,jp_nh4dep )  = znh4dep(ji,jj) * rno3 * zmsk                ! iron deposition
            END DO
         END DO
# endif
      !
      DEALLOCATE( zno3dep, znh4dep )
      !
      ENDIF

      ! OA: Warning, the following part is necessary to avoid CFL problems above the sediments
      ! --------------------------------------------------------------------
      DO jj = JRANGE
         DO ji = IRANGE
            zdep = e3t_n(ji,jj,KSED) / xstep
            zwsbio4(ji,jj) = MIN( 0.99 * zdep, wsbio4(ji,jj,ikt) )
            zwsbio3(ji,jj) = MIN( 0.99 * zdep, wsbio3(ji,jj,ikt) )
         END DO
      END DO
      !
      IF( .NOT.lk_sed ) THEN
!
         ! Add the external input of iron from sediment mobilization
         ! ------------------------------------------------------
         IF( ln_ironsed ) THEN
            tra(:,:,:,jpfer) = tra(:,:,:,jpfer) + ironsed(:,:,:) * rfact2
            !
#if defined key_iomput
            IF( lk_iomput .AND. knt == nrdttrc .AND. iom_use( "Ironsed" ) )   &
               &   CALL iom_put( "Ironsed", ironsed(:,:,:) * 1.e+3 * tmask(:,:,:) ) ! iron inputs from sediments
#endif
#if defined key_trc_diaadd
        DO jk = KRANGE
           DO jj = JRANGE
              DO ji = IRANGE
                 trc3d(ji,jj,K,jp_ironsed ) = ironsed(ji,jj,jk) * 1e+3 * tmask(ji,jj,K)  ! iron from  sediment
             END DO
           END DO
        ENDDO
#endif
         !
         ENDIF

         ! Computation of the sediment denitrification proportion: The metamodel from midlleburg (2006) is being used
         ! Computation of the fraction of organic matter that is permanently buried from Dunne's model
         ! -------------------------------------------------------
         DO jj = JRANGE
            DO ji = IRANGE
              IF( tmask(ji,jj,1) == 1 ) THEN
                 zflx = (  trb(ji,jj,KSED,jpgoc) * zwsbio4(ji,jj)   &
                   &     + trb(ji,jj,KSED,jppoc) * zwsbio3(ji,jj) )  * 1E3 * 1E6 / 1E4
                 zflx  = LOG10( MAX( 1E-3, zflx ) )
                 zo2   = LOG10( MAX( 10. , trb(ji,jj,KSED,jpoxy) * 1E6 ) )
                 zno3  = LOG10( MAX( 1.  , trb(ji,jj,KSED,jpno3) * 1E6 * rno3 ) )
                 zdep  = LOG10( gdepw_n(ji,jj,ikt+1) )
                 zdenit2d(ji,jj) = -2.2567 - 1.185 * zflx - 0.221 * zflx**2 - 0.3995 * zno3 * zo2 + 1.25 * zno3    &
                   &                + 0.4721 * zo2 - 0.0996 * zdep + 0.4256 * zflx * zo2
                 zdenit2d(ji,jj) = 10.0**( zdenit2d(ji,jj) )
                   !
                 zflx = (  trb(ji,jj,KSED,jpgoc) * zwsbio4(ji,jj)   &
                   &     + trb(ji,jj,KSED,jppoc) * zwsbio3(ji,jj) ) * 1E6
                 zbureff(ji,jj) = 0.013 + 0.53 * zflx**2 / ( 7.0 + zflx )**2
              ENDIF
            END DO
         END DO 
         !
      ENDIF

      ! This loss is scaled at each bottom grid cell for equilibrating the total budget of silica in the ocean.
      ! Thus, the amount of silica lost in the sediments equal the supply at the surface (dust+rivers)
      ! ------------------------------------------------------
      IF( .NOT.lk_sed )  zrivsil = 1.0 - sedsilfrac

      DO jj = JRANGE
         DO ji = IRANGE
            zdep = xstep / e3t_n(ji,jj,KSED) 
            zwsc = zwsbio4(ji,jj) * zdep
            zsiloss = trb(ji,jj,KSED,jpgsi) * zwsc
            zcaloss = trb(ji,jj,KSED,jpcal) * zwsc
            !
            tra(ji,jj,ikt,jpgsi) = tra(ji,jj,ikt,jpgsi) - zsiloss
            tra(ji,jj,ikt,jpcal) = tra(ji,jj,ikt,jpcal) - zcaloss
         END DO
      END DO
      !
      IF( .NOT.lk_sed ) THEN
         DO jj = JRANGE
            DO ji = IRANGE
               zdep = xstep / e3t_n(ji,jj,KSED) 
               zwsc = zwsbio4(ji,jj) * zdep
               zsiloss = trb(ji,jj,KSED,jpgsi) * zwsc
               zcaloss = trb(ji,jj,KSED,jpcal) * zwsc
               tra(ji,jj,ikt,jpsil) = tra(ji,jj,ikt,jpsil) + zsiloss * zrivsil 
               !
               zfactcal = MIN( excess(ji,jj,ikt), 0.2 )
               zfactcal = MIN( 1., 1.3 * ( 0.2 - zfactcal ) / ( 0.4 - zfactcal ) )
               zrivalk  = sedcalfrac * zfactcal
               tra(ji,jj,ikt,jptal) =  tra(ji,jj,ikt,jptal) + zcaloss * zrivalk * 2.0
               tra(ji,jj,ikt,jpdic) =  tra(ji,jj,ikt,jpdic) + zcaloss * zrivalk
               zsedcal(ji,jj) = (1.0 - zrivalk) * zcaloss * e3t_n(ji,jj,KSED) 
               zsedsi (ji,jj) = (1.0 - zrivsil) * zsiloss * e3t_n(ji,jj,KSED) 
            END DO
         END DO
      ENDIF
      !
      DO jj = JRANGE
         DO ji = IRANGE
            zdep = xstep / e3t_n(ji,jj,KSED) 
            zws4 = zwsbio4(ji,jj) * zdep
            zws3 = zwsbio3(ji,jj) * zdep
            tra(ji,jj,ikt,jpgoc) = tra(ji,jj,ikt,jpgoc) - trb(ji,jj,KSED,jpgoc) * zws4 
            tra(ji,jj,ikt,jppoc) = tra(ji,jj,ikt,jppoc) - trb(ji,jj,KSED,jppoc) * zws3
            tra(ji,jj,ikt,jpbfe) = tra(ji,jj,ikt,jpbfe) - trb(ji,jj,KSED,jpbfe) * zws4
            tra(ji,jj,ikt,jpsfe) = tra(ji,jj,ikt,jpsfe) - trb(ji,jj,KSED,jpsfe) * zws3
         END DO
      END DO
      !
      IF( ln_p5z ) THEN
         DO jj = JRANGE
            DO ji = IRANGE
               zdep = xstep / e3t_n(ji,jj,KSED) 
               zws4 = zwsbio4(ji,jj) * zdep
               zws3 = zwsbio3(ji,jj) * zdep
               tra(ji,jj,ikt,jpgon) = tra(ji,jj,ikt,jpgon) - trb(ji,jj,KSED,jpgon) * zws4
               tra(ji,jj,ikt,jppon) = tra(ji,jj,ikt,jppon) - trb(ji,jj,KSED,jppon) * zws3
               tra(ji,jj,ikt,jpgop) = tra(ji,jj,ikt,jpgop) - trb(ji,jj,KSED,jpgop) * zws4
               tra(ji,jj,ikt,jppop) = tra(ji,jj,ikt,jppop) - trb(ji,jj,KSED,jppop) * zws3
            END DO
         END DO
      ENDIF

      IF( .NOT.lk_sed ) THEN
         ! The 0.5 factor in zpdenit is to avoid negative NO3 concentration after
         ! denitrification in the sediments. Not very clever, but simpliest option.
         DO jj = JRANGE
            DO ji = IRANGE
               zdep = xstep / e3t_n(ji,jj,KSED) 
               zws4 = zwsbio4(ji,jj) * zdep
               zws3 = zwsbio3(ji,jj) * zdep
               zrivno3 = 1. - zbureff(ji,jj)
               zwstpoc = trb(ji,jj,KSED,jpgoc) * zws4 + trb(ji,jj,KSED,jppoc) * zws3
               zpdenit  = MIN( 0.5 * ( trb(ji,jj,KSED,jpno3) - rtrn )   &
                  &     / rdenit, zdenit2d(ji,jj) * zwstpoc * zrivno3 )
               z1pdenit = zwstpoc * zrivno3 - zpdenit
               zolimit = MIN( ( trb(ji,jj,KSED,jpoxy) - rtrn ) / o2ut, z1pdenit * ( 1.- nitrfac(ji,jj,ikt) ) )
               tra(ji,jj,ikt,jpdoc) = tra(ji,jj,ikt,jpdoc) + z1pdenit - zolimit
               tra(ji,jj,ikt,jppo4) = tra(ji,jj,ikt,jppo4) + zpdenit + zolimit
               tra(ji,jj,ikt,jpnh4) = tra(ji,jj,ikt,jpnh4) + zpdenit + zolimit
               tra(ji,jj,ikt,jpno3) = tra(ji,jj,ikt,jpno3) - rdenit * zpdenit
               tra(ji,jj,ikt,jpoxy) = tra(ji,jj,ikt,jpoxy) - zolimit * o2ut
               tra(ji,jj,ikt,jptal) = tra(ji,jj,ikt,jptal) + rno3 * (zolimit + (1.+rdenit) * zpdenit )
               tra(ji,jj,ikt,jpdic) = tra(ji,jj,ikt,jpdic) + zpdenit + zolimit 
               sdenit(ji,jj) = rdenit * zpdenit * e3t_n(ji,jj,KSED)
               zsedc(ji,jj)   = (1. - zrivno3) * zwstpoc * e3t_n(ji,jj,KSED)
               IF( ln_p5z ) THEN
                  zwstpop              = trb(ji,jj,KSED,jpgop) * zws4 + trb(ji,jj,KSED,jppop) * zws3
                  zwstpon              = trb(ji,jj,KSED,jpgon) * zws4 + trb(ji,jj,KSED,jppon) * zws3
                  tra(ji,jj,ikt,jpdon) = tra(ji,jj,ikt,jpdon) + ( z1pdenit - zolimit ) * zwstpon / (zwstpoc + rtrn)
                  tra(ji,jj,ikt,jpdop) = tra(ji,jj,ikt,jpdop) + ( z1pdenit - zolimit ) * zwstpop / (zwstpoc + rtrn)
               ENDIF
            END DO
         END DO
      ENDIF

      ! Nitrogen fixation process
      ! Small source iron from particulate inorganic iron
      !-----------------------------------
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zlight (ji,jj,jk) =  ( 1.- EXP( -etot_ndcy(ji,jj,jk) / diazolight ) ) * ( 1. - fr_i(ji,jj) ) 
               zsoufer(ji,jj,jk) = zlight(ji,jj,jk) * 2E-11 / ( 2E-11 + biron(ji,jj,jk) )
            END DO
         END DO
      END DO
      IF( ln_p4z ) THEN
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  !                      ! Potential nitrogen fixation dependant on temperature and iron
                  ztemp = tsn(ji,jj,K,jp_tem)
                  zmudia = MAX( 0.,-0.001096*ztemp**2 + 0.057*ztemp -0.637 ) / rno3
                  !       Potential nitrogen fixation dependant on temperature and iron
                  xdianh4 = trb(ji,jj,K,jpnh4) / ( concnnh4 + trb(ji,jj,K,jpnh4) )
                  xdiano3 = trb(ji,jj,K,jpno3) / ( concnno3   &
                     &    + trb(ji,jj,K,jpno3) ) * (1. - xdianh4)
                  zlim = ( 1.- xdiano3 - xdianh4 )
                  IF( zlim <= 0.1 )   zlim = 0.01
                  zfact = zlim * rfact2
                  ztrfer = biron(ji,jj,jk) / ( concfediaz + biron(ji,jj,jk) )
                  ztrpo4(ji,jj,jk) = trb(ji,jj,K,jppo4)   &
                    &              / ( 1E-6 + trb(ji,jj,K,jppo4) )
                  ztrdp = ztrpo4(ji,jj,jk)
                  nitrpot(ji,jj,jk) =  zmudia * r1_rday * zfact * MIN( ztrfer, ztrdp ) * zlight(ji,jj,jk)
               END DO
            END DO
         END DO
      ELSE       ! p5z
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  !                      ! Potential nitrogen fixation dependant on temperature and iron
                  ztemp = tsn(ji,jj,K,jp_tem)
                  zmudia = MAX( 0.,-0.001096*ztemp**2 + 0.057*ztemp -0.637 )  / rno3
                  !       Potential nitrogen fixation dependant on temperature and iron
                  xdianh4 = trb(ji,jj,K,jpnh4) / ( concnnh4 + trb(ji,jj,K,jpnh4) )
                  xdiano3 = trb(ji,jj,K,jpno3) / ( concnno3   &
                     &    + trb(ji,jj,K,jpno3) ) * (1. - xdianh4)
                  zlim = ( 1.- xdiano3 - xdianh4 )
                  IF( zlim <= 0.1 )   zlim = 0.01
                  zfact = zlim * rfact2
                  ztrfer = biron(ji,jj,jk) / ( concfediaz + biron(ji,jj,jk) )
                  ztrpo4(ji,jj,jk) = trb(ji,jj,K,jppo4)   &
                     &             / ( 1E-6 + trb(ji,jj,K,jppo4) )
                  ztrdop(ji,jj,jk) = trb(ji,jj,K,jpdop)   &
                     &             / ( 1E-6 + trb(ji,jj,K,jpdop) ) * (1. - ztrpo4(ji,jj,jk))
                  ztrdp = ztrpo4(ji,jj,jk) + ztrdop(ji,jj,jk)
                  nitrpot(ji,jj,jk) =  zmudia * r1_rday * zfact * MIN( ztrfer, ztrdp ) * zlight(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF

      ! Nitrogen change due to nitrogen fixation
      ! ----------------------------------------
      IF( ln_p4z ) THEN
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  zfact = nitrpot(ji,jj,jk) * nitrfix
                  tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zfact / 3.0
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zfact / 3.0
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) - zfact * 2.0 / 3.0
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zfact * 1.0 / 3.0
                  tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zfact * 1.0 / 3.0 * 2.0 / 3.0
                  tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zfact * 1.0 / 3.0 * 1.0 / 3.0
                  tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + ( o2ut + o2nit ) * zfact * 2.0 / 3.0 + o2nit * zfact / 3.0
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - 30E-6 * zfact * 1.0 / 3.0
                  tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + 30E-6 * zfact * 1.0 / 3.0 * 2.0 / 3.0
                  tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + 30E-6 * zfact * 1.0 / 3.0 * 1.0 / 3.0
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + 0.002 * 4E-10 * zsoufer(ji,jj,jk) * rfact2 / rday
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + concdnh4 / ( concdnh4 + trb(ji,jj,K,jppo4) ) &
                  &                     * 0.001 * trb(ji,jj,K,jpdoc) * xstep
              END DO
            END DO 
         END DO
      ELSE    ! p5z
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  zfact = nitrpot(ji,jj,jk) * nitrfix
                  tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zfact / 3.0
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zfact / 3.0
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) - 16.0 / 46.0 * zfact * ( 1.0 - 1.0 / 3.0 ) &
                  &                     * ztrpo4(ji,jj,jk) / (ztrpo4(ji,jj,jk) + ztrdop(ji,jj,jk) + rtrn)
                  tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + zfact * 1.0 / 3.0
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zfact * 1.0 / 3.0
                  tra(ji,jj,jk,jpdop) = tra(ji,jj,jk,jpdop) + 16.0 / 46.0 * zfact / 3.0  &
                  &                     - 16.0 / 46.0 * zfact * ztrdop(ji,jj,jk)   &
                  &                     / (ztrpo4(ji,jj,jk) + ztrdop(ji,jj,jk) + rtrn)
                  tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zfact * 1.0 / 3.0 * 2.0 / 3.0
                  tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + zfact * 1.0 / 3.0 * 2.0 /3.0
                  tra(ji,jj,jk,jppop) = tra(ji,jj,jk,jppop) + 16.0 / 46.0 * zfact * 1.0 / 3.0 * 2.0 /3.0
                  tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zfact * 1.0 / 3.0 * 1.0 / 3.0
                  tra(ji,jj,jk,jpgon) = tra(ji,jj,jk,jpgon) + zfact * 1.0 / 3.0 * 1.0 /3.0
                  tra(ji,jj,jk,jpgop) = tra(ji,jj,jk,jpgop) + 16.0 / 46.0 * zfact * 1.0 / 3.0 * 1.0 /3.0
                  tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + ( o2ut + o2nit ) * zfact * 2.0 / 3.0 + o2nit * zfact / 3.0
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - 30E-6 * zfact * 1.0 / 3.0 
                  tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + 30E-6 * zfact * 1.0 / 3.0 * 2.0 / 3.0
                  tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + 30E-6 * zfact * 1.0 / 3.0 * 1.0 / 3.0
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + 0.002 * 4E-10 * zsoufer(ji,jj,jk) * rfact2 / rday
              END DO
            END DO 
         END DO
         !
      ENDIF

#if defined key_iomput
      IF( lk_iomput ) THEN
         IF( knt == nrdttrc ) THEN
            zfact = 1.e+3 * rfact2r !  conversion from molC/l/kt  to molN/m3/s
            IF( iom_use("Nfix"   ) ) CALL iom_put( "Nfix", nitrpot(:,:,:) * nitrfix * rno3 * zfact * tmask(:,:,:) )  ! nitrogen fixation 
            IF( iom_use("INTNFIX") ) THEN   ! nitrogen fixation rate in ocean ( vertically integrated )
               zwork(:,:) = 0.
               DO jk = KRANGE
                  DO jj = JRANGE
                     DO ji = IRANGE
                        zwork(ji,jj) = zwork(ji,jj) + nitrpot(ji,jj,jk) * nitrfix * rno3    &
                        &              * zfact * e3t_n(ji,jj,K) * tmask(ji,jj,jk)
                     END DO
                  END DO
               END DO
               CALL iom_put( "INTNFIX" , zwork ) 
            ENDIF
            IF( iom_use("SedCal" ) ) CALL iom_put( "SedCal", zsedcal(:,:) * zfact )
            IF( iom_use("SedSi" ) )  CALL iom_put( "SedSi",  zsedsi (:,:) * zfact )
            IF( iom_use("SedC" ) )   CALL iom_put( "SedC",   zsedc  (:,:) * zfact )
            IF( iom_use("Sdenit" ) ) CALL iom_put( "Sdenit", sdenit (:,:) * zfact * rno3 )
         ENDIF
      ENDIF
#endif
      !
#if defined key_trc_diaadd
        zfact = 1.e+3 * rfact2r
        zwork(:,:) = 0.
        DO jk = KRANGE
           DO jj = JRANGE
              DO ji = IRANGE
                 zwork(ji,jj) = zwork(ji,jj) + nitrpot(ji,jj,jk) * nitrfix * rno3    &
                 &              * zfact * e3t_n(ji,jj,K) * tmask(ji,jj,jk)
              END DO
           END DO
        END DO

        DO jj = JRANGE
           DO ji = IRANGE
              trc2d(ji,jj,jp_nfix   )  = zwork(ji,jj)       ! nitrogen fixation at surface
           END DO
        END DO

        zrfact2 = 1.e+3 * rfact2r

        DO jk = KRANGE
           DO jj = JRANGE
              DO ji = IRANGE
                 zmsk = zrfact2 * tmask(ji,jj,K)
                 trc3d(ji,jj,K,jp_nfixo2 )  = nitrpot(ji,jj,jk) * rno3 * zmsk * nitrfix * o2nit  ! O2 production by Nfix
             END DO
           END DO
        ENDDO
# endif
      !
      IF(ln_ctl) THEN  ! print mean trends (USEd for debugging)
         WRITE(charout, fmt="('sed ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( charout, ltra='tra')
      ENDIF
      !
      IF( ln_p5z )    DEALLOCATE( ztrpo4, ztrdop )
      !
   END SUBROUTINE p4z_sed


   INTEGER FUNCTION p4z_sed_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( nitrpot(PRIV_3D_BIOARRAY), sdenit(PRIV_2D_BIOARRAY), STAT=p4z_sed_alloc )
      !
      IF( p4z_sed_alloc /= 0 )   CALL ctl_warn( 'p4z_sed_alloc: failed to allocate arrays' )
      !
   END FUNCTION p4z_sed_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sed                         ! Empty routine
   END SUBROUTINE p4z_sed
#endif 

   !!======================================================================
END MODULE p4zsed
