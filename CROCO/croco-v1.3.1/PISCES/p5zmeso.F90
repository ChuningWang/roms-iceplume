#include "cppdefs.h"

MODULE p5zmeso
   !!======================================================================
   !!                         ***  MODULE p5zmeso  ***
   !! TOP :   PISCES Compute the sources/sinks for mesozooplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!   p5z_meso       :   Compute the sources/sinks for mesozooplankton
   !!   p5z_meso_init  :   Initialization of the parameters for mesozooplankton
   !!----------------------------------------------------------------------
   USE sms_pisces      !  PISCES Source Minus Sink variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p5z_meso              ! called in p5zbio.F90
   PUBLIC   p5z_meso_init         ! called in trcsms_pisces.F90

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"

   !! * Shared module variables
   REAL(wp), PUBLIC ::  part2        !: part of calcite not dissolved in mesozoo guts
   REAL(wp), PUBLIC ::  xpref2c      !: mesozoo preference for POC 
   REAL(wp), PUBLIC ::  xpref2n      !: mesozoo preference for nanophyto
   REAL(wp), PUBLIC ::  xpref2z      !: mesozoo preference for zooplankton
   REAL(wp), PUBLIC ::  xpref2d      !: mesozoo preference for Diatoms 
   REAL(wp), PUBLIC ::  xpref2m      !: mesozoo preference for mesozoo
   REAL(wp), PUBLIC ::  xthresh2zoo  !: zoo feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2dia  !: diatoms feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2phy  !: nanophyto feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2poc  !: poc feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2mes  !: mesozoo feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2     !: feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  resrat2      !: exsudation rate of mesozooplankton
   REAL(wp), PUBLIC ::  mzrat2       !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::  grazrat2     !: maximal mesozoo grazing rate
   REAL(wp), PUBLIC ::  xkgraz2      !: Half-saturation constant of assimilation
   REAL(wp), PUBLIC ::  unass2c      !: Non-assimilated fraction of food
   REAL(wp), PUBLIC ::  unass2n      !: Non-assimilated fraction of food
   REAL(wp), PUBLIC ::  unass2p      !: Non-assimilated fraction of food
   REAL(wp), PUBLIC ::  epsher2      !: Growth efficiency of mesozoo
   REAL(wp), PUBLIC ::  epsher2min   !: Minimum growth efficiency of mesozoo
   REAL(wp), PUBLIC ::  ssigma2      !: Fraction excreted as semi-labile DOM
   REAL(wp), PUBLIC ::  srespir2     !: Active respiration
   REAL(wp), PUBLIC ::  grazflux     !: mesozoo flux feeding rate
   LOGICAL,  PUBLIC ::  bmetexc2     !: Use of excess carbon for respiration

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p5zmeso.F90 10362 2018-11-30 15:38:17Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p5z_meso( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_meso  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for mesozooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompadi, zcompaph, zcompapoc, zcompaz, zcompam, zcompames
      REAL(wp) :: zgraze2, zdenom, zfact, zfood, zfoodlim, zproport
      REAL(wp) :: zmortzgoc, zfracc, zfracn, zfracp, zfracfe, zratio, zratio2
      REAL(wp) :: zepsherf, zepshert, zepsherv, zrespirc, zrespirn, zrespirp, zbasresb, zbasresi
      REAL(wp) :: zgraztotc, zgraztotn, zgraztotp, zgraztotf, zbasresn, zbasresp, zbasresf
      REAL(wp) :: zgradoc, zgradon, zgradop, zgratmp, zgradoct, zgradont, zgrareft, zgradopt
      REAL(wp) :: zgrapoc, zgrapon, zgrapop, zgrapof, zprcaca, zmortz
      REAL(wp) :: zexcess, zgrarem, zgraren, zgrarep, zgraref
      REAL(wp) :: zbeta, zrespz, ztortz, zgrasratp, zgrasratn, zgrasratf
      REAL(wp) :: ztmp1, ztmp2, ztmp3, ztmp4, ztmp5, ztmptot
      REAL(wp) :: zgrazdc, zgrazz, zgrazm, zgrazpof, zgrazcal, zfracal
      REAL(wp) :: zgraznc, zgrazpoc, zgrazpon, zgrazpop, zgraznf, zgrazdf
      REAL(wp) :: zgraznp, zgraznn, zgrazdn, zgrazdp
      REAL(wp) :: zgrazfffp, zgrazfffg, zgrazffep, zgrazffeg
      REAL(wp) :: zgrazffnp, zgrazffng, zgrazffpp, zgrazffpg
      CHARACTER (len=25) :: charout
      REAL(wp) :: zrfact2, zmetexcess
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY) :: zgrazing, zfezoo2
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d, zz2ligprod

      !!---------------------------------------------------------------------
      !
      zgrazing(:,:,:) = 0.
      zfezoo2 (:,:,:) = 0.
      !
      IF (ln_ligand) THEN
         ALLOCATE( zz2ligprod(PRIV_3D_BIOARRAY) )
         zz2ligprod(:,:,:) = 0.
      ENDIF

      zmetexcess = 0.0
      IF ( bmetexc2 ) zmetexcess = 1.0

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zcompam   = MAX( ( trb(ji,jj,K,jpmes) - 1.e-9 ), 0.e0 )
               zfact     = xstep * tgfunc2(ji,jj,jk) * zcompam

               !   Michaelis-Menten mortality rates of mesozooplankton
               !   ---------------------------------------------------
               zrespz   = resrat2 * zfact * ( trb(ji,jj,K,jpmes)   &
               &          / ( xkmort + trb(ji,jj,K,jpmes) )  &
               &          + 3. * nitrfac(ji,jj,jk) )

               !   Zooplankton mortality. A square function has been selected with
               !   no real reason except that it seems to be more stable and may mimic predation
               !   ---------------------------------------------------------------
               ztortz   = mzrat2 * 1.e6 * zfact * trb(ji,jj,K,jpmes) * (1. - nitrfac(ji,jj,jk))

               !   Computation of the abundance of the preys
               !   A threshold can be specified in the namelist
               !   --------------------------------------------
               zcompadi  = MAX( ( trb(ji,jj,K,jpdia) - xthresh2dia ), 0.e0 )
               zcompaz   = MAX( ( trb(ji,jj,K,jpzoo) - xthresh2zoo ), 0.e0 )
               zcompaph  = MAX( ( trb(ji,jj,K,jpphy) - xthresh2phy ), 0.e0 )
               zcompapoc = MAX( ( trb(ji,jj,K,jppoc) - xthresh2poc ), 0.e0 )
               zcompames = MAX( ( trb(ji,jj,K,jpmes) - xthresh2mes ), 0.e0 )

               !   Mesozooplankton grazing
               !   ------------------------
               zfood     = xpref2d * zcompadi + xpref2z * zcompaz + xpref2n * zcompaph + xpref2c * zcompapoc   &
               &           + xpref2m * zcompames 
               zfoodlim  = MAX( 0., zfood - MIN( 0.5 * zfood, xthresh2 ) )
               zdenom    = zfoodlim / ( xkgraz2 + zfoodlim )
               zgraze2   = grazrat2 * xstep * tgfunc2(ji,jj,jk)   &
               &          * trb(ji,jj,K,jpmes) * (1. - nitrfac(ji,jj,jk)) 

               !   An active switching parameterization is used here.
               !   We don't use the KTW parameterization proposed by 
               !   Vallina et al. because it tends to produce to steady biomass
               !   composition and the variance of Chl is too low as it grazes
               !   too strongly on winning organisms. Thus, instead of a square
               !   a 1.5 power value is used which decreases the pressure on the
               !   most abundant species
               !   ------------------------------------------------------------  
               ztmp1 = xpref2n * zcompaph**1.5
               ztmp2 = xpref2m * zcompames**1.5
               ztmp3 = xpref2c * zcompapoc**1.5
               ztmp4 = xpref2d * zcompadi**1.5
               ztmp5 = xpref2z * zcompaz**1.5
               ztmptot = ztmp1 + ztmp2 + ztmp3 + ztmp4 + ztmp5 + rtrn
               ztmp1 = ztmp1 / ztmptot
               ztmp2 = ztmp2 / ztmptot
               ztmp3 = ztmp3 / ztmptot
               ztmp4 = ztmp4 / ztmptot
               ztmp5 = ztmp5 / ztmptot

               !   Mesozooplankton regular grazing on the different preys
               !   ------------------------------------------------------
               zgrazdc   = zgraze2 * ztmp4 * zdenom
               zgrazdn   = zgrazdc * trb(ji,jj,K,jpndi)   &
               &          / ( trb(ji,jj,K,jpdia) + rtrn)
               zgrazdp   = zgrazdc * trb(ji,jj,K,jppdi)   &
               &          / ( trb(ji,jj,K,jpdia) + rtrn)
               zgrazdf   = zgrazdc * trb(ji,jj,K,jpdfe)   &
               &          / ( trb(ji,jj,K,jpdia) + rtrn)
               zgrazz    = zgraze2 * ztmp5 * zdenom
               zgrazm    = zgraze2 * ztmp2 * zdenom
               zgraznc   = zgraze2 * ztmp1 * zdenom
               zgraznn   = zgraznc * trb(ji,jj,K,jpnph)   &
               &          / ( trb(ji,jj,K,jpphy) + rtrn)
               zgraznp   = zgraznc * trb(ji,jj,K,jppph)   &
               &          / ( trb(ji,jj,K,jpphy) + rtrn)
               zgraznf   = zgraznc * trb(ji,jj,K,jpnfe)   &
               &          / ( trb(ji,jj,K,jpphy) + rtrn)
               zgrazpoc  = zgraze2 * ztmp3 * zdenom
               zgrazpon  = zgrazpoc * trb(ji,jj,K,jppon)  &
               &          / ( trb(ji,jj,K,jppoc) + rtrn)
               zgrazpop  = zgrazpoc * trb(ji,jj,K,jppop)  &
               &          / ( trb(ji,jj,K,jppoc) + rtrn)
               zgrazpof  = zgrazpoc * trb(ji,jj,K,jpsfe)  &
               &          / ( trb(ji,jj,K,jppoc) + rtrn)

               !   Mesozooplankton flux feeding on GOC
               !   ----------------------------------
               zgrazffeg = grazflux  * xstep * wsbio4(ji,jj,jk)      &
               &           * tgfunc2(ji,jj,jk) * trb(ji,jj,K,jpgoc)  &
               &           * trb(ji,jj,K,jpmes)  &
               &           * (1. - nitrfac(ji,jj,jk))
               zgrazfffg = zgrazffeg * trb(ji,jj,K,jpbfe)   &
               &          / (trb(ji,jj,K,jpgoc) + rtrn)
               zgrazffng = zgrazffeg * trb(ji,jj,K,jpgon)   &
               &          / (trb(ji,jj,K,jpgoc) + rtrn)
               zgrazffpg = zgrazffeg * trb(ji,jj,K,jpgop)   &
               &          / (trb(ji,jj,K,jpgoc) + rtrn)
               zgrazffep = grazflux  * xstep *  wsbio3(ji,jj,jk)     &
               &           * tgfunc2(ji,jj,jk) * trb(ji,jj,K,jppoc)  &
               &           * trb(ji,jj,K,jpmes)   &
               &           * (1. - nitrfac(ji,jj,jk))
               zgrazfffp = zgrazffep * trb(ji,jj,K,jpsfe)   &
               &           / (trb(ji,jj,K,jppoc) + rtrn)
               zgrazffnp = zgrazffep * trb(ji,jj,K,jppon)   &
               &           / (trb(ji,jj,K,jppoc) + rtrn)
               zgrazffpp = zgrazffep * trb(ji,jj,K,jppop)   &
               &           / (trb(ji,jj,K,jppoc) + rtrn)
               !
               zgraztotc  = zgrazdc + zgrazz + zgraznc + zgrazm + zgrazpoc + zgrazffep + zgrazffeg

               !   Compute the proportion of filter feeders
               !   ----------------------------------------  
               zproport  = (zgrazffep + zgrazffeg)/(rtrn + zgraztotc)

               !   Compute fractionation of aggregates. It is assumed that 
               !   diatoms based aggregates are more prone to fractionation
               !   since they are more porous (marine snow instead of fecal pellets)
               !   ----------------------------------------------------------------
               zratio    = trb(ji,jj,K,jpgsi) / ( trb(ji,jj,K,jpgoc) + rtrn )
               zratio2   = zratio * zratio
               zfracc    = zproport * grazflux  * xstep * wsbio4(ji,jj,jk)      &
               &          * trb(ji,jj,K,jpgoc) * trb(ji,jj,K,jpmes)          &
               &          * ( 0.2 + 3.8 * zratio2 / ( 1.**2 + zratio2 ) )
               zfracfe   = zfracc * trb(ji,jj,K,jpbfe) / (trb(ji,jj,K,jpgoc) + rtrn)
               zfracn    = zfracc * trb(ji,jj,K,jpgon) / (trb(ji,jj,K,jpgoc) + rtrn)
               zfracp    = zfracc * trb(ji,jj,K,jpgop) / (trb(ji,jj,K,jpgoc) + rtrn)

               zgrazffep = zproport * zgrazffep   ;   zgrazffeg = zproport * zgrazffeg
               zgrazfffp = zproport * zgrazfffp   ;   zgrazfffg = zproport * zgrazfffg
               zgrazffnp = zproport * zgrazffnp   ;   zgrazffng = zproport * zgrazffng
               zgrazffpp = zproport * zgrazffpp   ;   zgrazffpg = zproport * zgrazffpg

               zgraztotc  = zgrazdc + zgrazz + zgraznc + zgrazm + zgrazpoc + zgrazffep + zgrazffeg
               zgraztotf  = zgrazdf + zgraznf + ( zgrazz + zgrazm ) * ferat3 + zgrazpof &
               &            + zgrazfffp + zgrazfffg
               zgraztotn  = zgrazdn + (zgrazm + zgrazz) * no3rat3 + zgraznn + zgrazpon  &
               &            + zgrazffnp + zgrazffng
               zgraztotp  = zgrazdp + (zgrazz + zgrazm) * po4rat3 + zgraznp + zgrazpop  &
               &            + zgrazffpp + zgrazffpg


               ! Total grazing ( grazing by microzoo is already computed in p5zmicro )
               zgrazing(ji,jj,jk) = zgraztotc

               !   Stoichiometruc ratios of the food ingested by zooplanton 
               !   --------------------------------------------------------
               zgrasratf  =  (zgraztotf + rtrn) / ( zgraztotc + rtrn )
               zgrasratn  =  (zgraztotn + rtrn) / ( zgraztotc + rtrn )
               zgrasratp  =  (zgraztotp + rtrn) / ( zgraztotc + rtrn )

               !   Growth efficiency is made a function of the quality 
               !   and the quantity of the preys
               !   ---------------------------------------------------
               zepshert  = MIN( 1., zgrasratn/ no3rat3, zgrasratp/ po4rat3, zgrasratf / ferat3)
               zbeta     = MAX(0., (epsher2 - epsher2min) )
               zepsherf  = epsher2min + zbeta / ( 1.0 + 0.04E6 * 12. * zfood * zbeta )
               zepsherv  = zepsherf * zepshert

               !   Respiration of mesozooplankton
               !   Excess carbon in the food is used preferentially
               !   ----------------  ------------------------------
               zexcess  = zgraztotc * zepsherf * (1.0 - zepshert) * zmetexcess 
               zbasresb = MAX(0., zrespz - zexcess)
               zbasresi = zexcess + MIN(0., zrespz - zexcess)
               zrespirc = srespir2 * zepsherv * zgraztotc + zbasresb

               !   When excess carbon is used, the other elements in excess
               !   are also used proportionally to their abundance
               !   --------------------------------------------------------
               zexcess  = ( zgrasratn/ no3rat3 - zepshert ) / ( 1.0 - zepshert + rtrn)
               zbasresn = zbasresi * zexcess * zgrasratn
               zexcess  = ( zgrasratp/ po4rat3 - zepshert ) / ( 1.0 - zepshert + rtrn)
               zbasresp = zbasresi * zexcess * zgrasratp
               zexcess  = ( zgrasratf/ ferat3 - zepshert ) / ( 1.0 - zepshert + rtrn)
               zbasresf = zbasresi * zexcess * zgrasratf

               !   Voiding of the excessive elements as organic matter
               !   --------------------------------------------------------
               zgradoct = (1. - unass2c - zepsherv) * zgraztotc - zbasresi
               zgradont = (1. - unass2n) * zgraztotn - zepsherv * no3rat3 * zgraztotc - zbasresn
               zgradopt = (1. - unass2p) * zgraztotp - zepsherv * po4rat3 * zgraztotc - zbasresp
               zgrareft = (1. - unass2c) * zgraztotf - zepsherv * ferat3 * zgraztotc - zbasresf
               ztmp1   = ( 1. - epsher2 - unass2c ) /( 1. - 0.8 * epsher2 ) * ztortz
               zgradoc = (zgradoct + ztmp1) * ssigma2
               zgradon = (zgradont + no3rat3 * ztmp1) * ssigma2
               zgradop = (zgradopt + po4rat3 * ztmp1) * ssigma2
               zgratmp = 0.2 * epsher2 /( 1. - 0.8 * epsher2 ) * ztortz

               !  Since only semilabile DOM is represented in PISCES
               !  part of DOM is in fact labile and is then released
               !  as dissolved inorganic compounds (ssigma2)
               !  --------------------------------------------------
               zgrarem = zgratmp + ( zgradoct + ztmp1 ) * (1.0 - ssigma2)
               zgraren = no3rat3 * zgratmp + ( zgradont + no3rat3 * ztmp1 ) * (1.0 - ssigma2)
               zgrarep = po4rat3 * zgratmp + ( zgradopt + po4rat3 * ztmp1 ) * (1.0 - ssigma2)
               zgraref = zgrareft + ferat3 * ( ztmp1 + zgratmp )

               !   Defecation as a result of non assimilated products
               !   --------------------------------------------------
               zgrapoc  = zgraztotc * unass2c + unass2c / ( 1. - 0.8 * epsher2 ) * ztortz
               zgrapon  = zgraztotn * unass2n + no3rat3 * unass2n / ( 1. - 0.8 * epsher2 ) * ztortz
               zgrapop  = zgraztotp * unass2p + po4rat3 * unass2p / ( 1. - 0.8 * epsher2 ) * ztortz
               zgrapof  = zgraztotf * unass2c + ferat3  * unass2c / ( 1. - 0.8 * epsher2 ) * ztortz

               !  Addition of respiration to the release of inorganic nutrients
               !  -------------------------------------------------------------
               zgrarem = zgrarem + zbasresi + zrespirc
               zgraren = zgraren + zbasresn + zrespirc * no3rat3
               zgrarep = zgrarep + zbasresp + zrespirc * po4rat3
               zgraref = zgraref + zbasresf + zrespirc * ferat3

               !   Update the arrays TRA which contain the biological sources and
               !   sinks
               !   --------------------------------------------------------------
               tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zgrarep 
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zgraren
               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zgradoc
               !
               IF( ln_ligand ) THEN
                  tra(ji,jj,jk,jplgw)  = tra(ji,jj,jk,jplgw) + zgradoc * ldocz
                  zz2ligprod(ji,jj,jk) = zgradoc * ldocz
               ENDIF
               !
               tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + zgradon
               tra(ji,jj,jk,jpdop) = tra(ji,jj,jk,jpdop) + zgradop
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2ut * zgrarem
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zgraref
               zfezoo2(ji,jj,jk)   = zgraref
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrarem
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zgraren
               tra(ji,jj,jk,jpmes) = tra(ji,jj,jk,jpmes) + zepsherv * zgraztotc - zrespirc   &
               &                     - ztortz - zgrazm
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zgrazdc
               tra(ji,jj,jk,jpndi) = tra(ji,jj,jk,jpndi) - zgrazdn
               tra(ji,jj,jk,jppdi) = tra(ji,jj,jk,jppdi) - zgrazdp
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zgrazdf
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - zgrazz
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zgraznc
               tra(ji,jj,jk,jpnph) = tra(ji,jj,jk,jpnph) - zgraznn
               tra(ji,jj,jk,jppph) = tra(ji,jj,jk,jppph) - zgraznp
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zgraznf
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zgraznc * trb(ji,jj,K,jpnch)   &
               &                    / ( trb(ji,jj,K,jpphy) + rtrn )
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - zgrazdc * trb(ji,jj,K,jpdch)   &
               &                    / ( trb(ji,jj,K,jpdia) + rtrn )
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) - zgrazdc * trb(ji,jj,K,jpdsi)   &
               &                    / ( trb(ji,jj,K,jpdia) + rtrn )
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) + zgrazdc * trb(ji,jj,K,jpdsi)   &
               &                    / ( trb(ji,jj,K,jpdia) + rtrn )

               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zgrazpoc - zgrazffep + zfracc
               prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zfracc
               conspoc(ji,jj,jk) = conspoc(ji,jj,jk) - zgrazpoc - zgrazffep
               tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) - zgrazpon - zgrazffnp + zfracn
               tra(ji,jj,jk,jppop) = tra(ji,jj,jk,jppop) - zgrazpop - zgrazffpp + zfracp
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) - zgrazffeg + zgrapoc - zfracc
               prodgoc(ji,jj,jk) = prodgoc(ji,jj,jk) + zgrapoc
               consgoc(ji,jj,jk) = consgoc(ji,jj,jk) - zgrazffeg - zfracc
               tra(ji,jj,jk,jpgon) = tra(ji,jj,jk,jpgon) - zgrazffng + zgrapon - zfracn
               tra(ji,jj,jk,jpgop) = tra(ji,jj,jk,jpgop) - zgrazffpg + zgrapop - zfracp
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) - zgrazpof - zgrazfffp + zfracfe
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) - zgrazfffg + zgrapof - zfracfe
               zfracal = trb(ji,jj,K,jpcal) / ( trb(ji,jj,K,jpgoc) + rtrn )
               zgrazcal = zgrazffeg * (1. - part2) * zfracal

               !  calcite production
               !  ------------------
               zprcaca = xfracal(ji,jj,jk) * zgraznc
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
               zprcaca = part2 * zprcaca
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrazcal - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + 2. * ( zgrazcal - zprcaca )
               tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) - zgrazcal + zprcaca
            END DO
         END DO
      END DO
      !
#if defined key_iomput
      IF( lk_iomput .AND. knt == nrdttrc ) THEN
         ALLOCATE( zw3d(PRIV_3D_BIOARRAY) )
         IF( iom_use( "GRAZ2" ) ) THEN
            zw3d(:,:,:) = zgrazing(:,:,:) * 1.e+3 * rfact2r * tmask(:,:,:)  !   Total grazing of phyto by zooplankton
            CALL iom_put( "GRAZ2", zw3d )
         ENDIF
         IF( iom_use( "PCAL" ) ) THEN
            zw3d(:,:,:) = prodcal(:,:,:) * 1.e+3 * rfact2r * tmask(:,:,:)   !  Calcite production
            CALL iom_put( "PCAL", zw3d )
         ENDIF
         IF( iom_use( "FEZOO2" ) ) THEN
            zw3d(:,:,:) = zfezoo2(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:)   !
            CALL iom_put( "FEZOO2", zw3d )
         ENDIF
         IF( iom_use( "LPRODZ2" ) .AND. ln_ligand )  THEN
            zw3d(:,:,:) = zz2ligprod(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:)
            CALL iom_put( "LPRODZ2"  , zw3d )
         ENDIF
         DEALLOCATE( zw3d )
      ENDIF
#endif
      !
#if defined key_trc_diaadd
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               trc3d(ji,jj,K,jp_grapoc2) = zgrazing(ji,jj,jk) * 1.e+3 * rfact2r * tmask(ji,jj,jk) !  grazing of phyto by mesozoo
            END DO
         END DO
      END DO

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               trc3d(ji,jj,K,jp_meso2) = zgrazing(ji,jj,jk) * ( 1. - epsher2 - unass2c ) &
                     &                   * (-o2ut) * ssigma2 * 1.e+3 * rfact2r * tmask(ji,jj,jk) ! o2 consumption by Mesozoo
            END DO
         END DO
      END DO
#endif
      !
      IF (ln_ligand)  DEALLOCATE( zz2ligprod )
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('meso')")
        CALL prt_ctl_trc_info(charout)
        CALL prt_ctl_trc( charout, ltra='tra')
      ENDIF
      !
   END SUBROUTINE p5z_meso


   SUBROUTINE p5z_meso_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p5z_meso_init  ***
      !!
      !! ** Purpose :   Initialization of mesozooplankton parameters
      !!
      !! ** Method  :   Read the nampismes namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampismes
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!
      NAMELIST/namp5zmes/part2, bmetexc2, grazrat2, resrat2, mzrat2, xpref2c, xpref2n, xpref2z, &
         &                xpref2m, xpref2d, xthresh2dia, xthresh2phy, xthresh2zoo, xthresh2poc, &
         &                xthresh2mes, xthresh2, xkgraz2, epsher2, epsher2min, ssigma2, unass2c, &
         &                unass2n, unass2p, srespir2, grazflux
      !!----------------------------------------------------------------------
      !
      REWIND( numnatp_ref )              ! Namelist nampismes in reference namelist : Pisces mesozooplankton
      READ  ( numnatp_ref, namp5zmes, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismes in reference namelist', lwp )
      !
      REWIND( numnatp_cfg )              ! Namelist nampismes in configuration namelist : Pisces mesozooplankton
      READ  ( numnatp_cfg, namp5zmes, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'nampismes in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, namp5zmes )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' ' 
         WRITE(numout,*) ' Namelist parameters for mesozooplankton, namp5zmes'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    part of calcite not dissolved in mesozoo guts  part2       = ', part2
         WRITE(numout,*) '    mesozoo preference for nano.                   xpref2n     = ', xpref2n
         WRITE(numout,*) '    mesozoo preference for diatoms                 xpref2d     = ', xpref2d
         WRITE(numout,*) '    mesozoo preference for zoo                     xpref2z     = ', xpref2z
         WRITE(numout,*) '    mesozoo preference for mesozoo                 xpref2m     = ', xpref2m
         WRITE(numout,*) '    mesozoo preference for poc                     xpref2c     = ', xpref2c
         WRITE(numout,*) '    microzoo feeding threshold  for mesozoo        xthresh2zoo = ', xthresh2zoo
         WRITE(numout,*) '    diatoms feeding threshold  for mesozoo         xthresh2dia = ', xthresh2dia
         WRITE(numout,*) '    nanophyto feeding threshold for mesozoo        xthresh2phy = ', xthresh2phy
         WRITE(numout,*) '    poc feeding threshold for mesozoo              xthresh2poc = ', xthresh2poc
         WRITE(numout,*) '    mesozoo feeding threshold for mesozoo          xthresh2mes = ', xthresh2mes
         WRITE(numout,*) '    feeding threshold for mesozooplankton          xthresh2    = ', xthresh2
         WRITE(numout,*) '    exsudation rate of mesozooplankton             resrat2     = ', resrat2
         WRITE(numout,*) '    mesozooplankton mortality rate                 mzrat2      = ', mzrat2
         WRITE(numout,*) '    maximal mesozoo grazing rate                   grazrat2    = ', grazrat2
         WRITE(numout,*) '    mesozoo flux feeding rate                      grazflux    = ', grazflux
         WRITE(numout,*) '    C egested fraction of food by mesozoo          unass2c     = ', unass2c
         WRITE(numout,*) '    N egested fraction of food by mesozoo          unass2n     = ', unass2n
         WRITE(numout,*) '    P egested fraction of food by mesozoo          unass2p     = ', unass2p
         WRITE(numout,*) '    Efficicency of Mesozoo growth                  epsher2     = ', epsher2
         WRITE(numout,*) '    Minimum Efficiency of Mesozoo growth           epsher2min  =', epsher2min
         WRITE(numout,*) '    Fraction excreted as semi-labile DOM           ssigma2     = ', ssigma2
         WRITE(numout,*) '    Active respiration                             srespir2    = ', srespir2
         WRITE(numout,*) '    half sturation constant for grazing 2          xkgraz2     = ', xkgraz2
         WRITE(numout,*) '    Use excess carbon for respiration              bmetexc2    = ', bmetexc2
      ENDIF
      !
   END SUBROUTINE p5z_meso_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p5z_meso                    ! Empty routine
   END SUBROUTINE p5z_meso
#endif 

   !!======================================================================
END MODULE p5zmeso
