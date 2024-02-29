#include "cppdefs.h"

MODULE p5zlim
   !!======================================================================
   !!                         ***  MODULE p5zlim  ***
   !! TOP :   PISCES with variable stoichiometry 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-04  (O. Aumont, C. Ethe) Limitation for iron modelled in quota 
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   p5z_lim        :   Compute the nutrients limitation terms 
   !!   p5z_lim_init   :   Read the namelist 
   !!----------------------------------------------------------------------
   USE p4zlim
   USE sms_pisces      ! PISCES variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC p5z_lim    
   PUBLIC p5z_lim_init    
   PUBLIC p5z_lim_alloc

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"

   !! * Shared module variables
   REAL(wp), PUBLIC ::  concpno3    !:  NO3, PO4 half saturation   
   REAL(wp), PUBLIC ::  concpnh4    !:  NH4 half saturation for phyto  
   REAL(wp), PUBLIC ::  concnpo4    !:  NH4 half saturation for diatoms
   REAL(wp), PUBLIC ::  concppo4    !:  NH4 half saturation for diatoms
   REAL(wp), PUBLIC ::  concdpo4    !:  NH4 half saturation for diatoms
   REAL(wp), PUBLIC ::  concpfer    !:  Iron half saturation for nanophyto 
   REAL(wp), PUBLIC ::  concbpo4    !:  PO4 half saturation for bacteria
   REAL(wp), PUBLIC ::  xsizepic    !:  Minimum size criteria for diatoms
   REAL(wp), PUBLIC ::  xsizerp     !:  Size ratio for nanophytoplankton
   REAL(wp), PUBLIC ::  qfnopt      !:  optimal Fe quota for nanophyto
   REAL(wp), PUBLIC ::  qfpopt      !:  optimal Fe quota for nanophyto
   REAL(wp), PUBLIC ::  qfdopt      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qnnmin      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qnnmax      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qpnmin      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qpnmax      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qnpmin      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qnpmax      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qppmin      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qppmax      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qndmin      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qndmax      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qpdmin      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qpdmax      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qfnmax      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qfpmax      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qfdmax      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  zpsinh4
   REAL(wp), PUBLIC ::  zpsino3
   REAL(wp), PUBLIC ::  zpsiuptk

   !!*  Allometric variations of the quotas
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqnnmin    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqnnmax    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpnmin    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpnmax    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqnpmin    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqnpmax    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqppmin    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqppmax    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqndmin    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqndmax    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpdmin    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpdmax    !: ???

   !!* Phytoplankton limitation terms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpicono3   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpiconh4   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpicopo4   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanodop   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpicodop   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatdop   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanofer   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpicofer   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatfer   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimpic    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimpfe    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   fvnuptk
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   fvpuptk
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   fvduptk

   ! Coefficient for iron limitation
   REAL(wp) ::  xcoef1   = 0.00167  / 55.85
   REAL(wp) ::  xcoef2   = 1.21E-5 * 14. / 55.85 / 7.625 * 0.5 * 1.5
   REAL(wp) ::  xcoef3   = 1.15E-4 * 14. / 55.85 / 7.625 * 0.5 
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p5zlim.F90 10070 2018-08-28 14:30:54Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p5z_lim( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_lim  ***
      !!
      !! ** Purpose :   Compute the co-limitations by the various nutrients
      !!              for the various phytoplankton species
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in)  :: kt, knt
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zlim1, zlim2, zlim3, zlim4, zno3, zferlim
      REAL(wp) ::   z1_trndia, z1_trnpic, z1_trnphy, ztem1, ztem2, zetot1
      REAL(wp) ::   zratio, zration, zratiof, znutlim, zfalim
      REAL(wp) ::   zconc1d, zconc1dnh4, zconc0n, zconc0nnh4, zconc0npo4, zconc0dpo4
      REAL(wp) ::   zconc0p, zconc0pnh4, zconc0ppo4, zconcpfe, zconcnfe, zconcdfe
      REAL(wp) ::   fanano, fananop, fananof, fadiat, fadiatp, fadiatf
      REAL(wp) ::   fapico, fapicop, fapicof
      REAL(wp) ::   zrpho, zrass, zcoef, zfuptk, zratchl
      REAL(wp) ::   zfvn, zfvp, zfvf, zsizen, zsizep, zsized, znanochl, zpicochl, zdiatchl
      REAL(wp) ::   zqfemn, zqfemp, zqfemd, zbactno3, zbactnh4
      !!---------------------------------------------------------------------
      !
      zratchl = 6.0
      !
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               ! 
               ! Tuning of the iron concentration to a minimum level that is set to the detection limit
               !-------------------------------------
               zno3    = trb(ji,jj,K,jpno3) / 40.e-6
               zferlim = MAX( 3e-11 * zno3 * zno3, 5e-12 )
               zferlim = MIN( zferlim, 7e-11 )
               trb(ji,jj,K,jpfer) = MAX( trb(ji,jj,K,jpfer), zferlim )

               ! Computation of the mean relative size of each community
               ! -------------------------------------------------------
               z1_trnphy   = 1. / ( trb(ji,jj,K,jpphy) + rtrn )
               z1_trnpic   = 1. / ( trb(ji,jj,K,jppic) + rtrn )
               z1_trndia   = 1. / ( trb(ji,jj,K,jpdia) + rtrn )
               znanochl = trb(ji,jj,K,jpnch) * z1_trnphy
               zpicochl = trb(ji,jj,K,jppch) * z1_trnpic
               zdiatchl = trb(ji,jj,K,jpdch) * z1_trndia

               ! Computation of a variable Ks for iron on diatoms taking into account
               ! that increasing biomass is made of generally bigger cells
               !------------------------------------------------
               zsized            = sized(ji,jj,jk)**0.81
               zconcdfe          = concdfer * zsized
               zconc1d           = concdno3 * zsized
               zconc1dnh4        = concdnh4 * zsized
               zconc0dpo4        = concdpo4 * zsized

               zsizep            = 1.
               zconcpfe          = concpfer * zsizep
               zconc0p           = concpno3 * zsizep
               zconc0pnh4        = concpnh4 * zsizep
               zconc0ppo4        = concppo4 * zsizep

               zsizen            = 1.
               zconcnfe          = concnfer * zsizen
               zconc0n           = concnno3 * zsizen
               zconc0nnh4        = concnnh4 * zsizen
               zconc0npo4        = concnpo4 * zsizen

               ! Allometric variations of the minimum and maximum quotas
               ! From Talmy et al. (2014) and Maranon et al. (2013)
               ! -------------------------------------------------------
               xqnnmin(ji,jj,jk) = qnnmin
               xqnnmax(ji,jj,jk) = qnnmax
               xqndmin(ji,jj,jk) = qndmin * sized(ji,jj,jk)**(-0.27) 
               xqndmax(ji,jj,jk) = qndmax
               xqnpmin(ji,jj,jk) = qnpmin
               xqnpmax(ji,jj,jk) = qnpmax

               ! Computation of the optimal allocation parameters
               ! Based on the different papers by Pahlow et al., and Smith et al.
               ! -----------------------------------------------------------------
               znutlim = MAX( trb(ji,jj,K,jpnh4) / zconc0nnh4,    &
                 &         trb(ji,jj,K,jpno3) / zconc0n)
               fanano = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               znutlim = trb(ji,jj,K,jppo4) / zconc0npo4
               fananop = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               znutlim = biron(ji,jj,jk) / zconcnfe
               fananof = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               znutlim = MAX( trb(ji,jj,K,jpnh4) / zconc0pnh4,    &
                 &         trb(ji,jj,K,jpno3) / zconc0p)
               fapico = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               znutlim = trb(ji,jj,K,jppo4) / zconc0ppo4
               fapicop = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               znutlim = biron(ji,jj,jk) / zconcpfe
               fapicof = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               znutlim = MAX( trb(ji,jj,K,jpnh4) / zconc1dnh4,    &
                 &         trb(ji,jj,K,jpno3) / zconc1d )
               fadiat = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               znutlim = trb(ji,jj,K,jppo4) / zconc0dpo4
               fadiatp = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               znutlim = biron(ji,jj,jk) / zconcdfe
               fadiatf = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               !
               ! Michaelis-Menten Limitation term for nutrients Small bacteria
               ! -------------------------------------------------------------
               zbactnh4 = trb(ji,jj,K,jpnh4) / ( concbnh4 + trb(ji,jj,K,jpnh4) )
               zbactno3 = trb(ji,jj,K,jpno3) / ( concbno3   &
                  &     + trb(ji,jj,K,jpno3) ) * (1. - zbactnh4)
               !
               zlim1    = zbactno3 + zbactnh4
               zlim2    = trb(ji,jj,K,jppo4) / ( trb(ji,jj,K,jppo4) + concbpo4)
               zlim3    = biron(ji,jj,jk) / ( concbfe + biron(ji,jj,jk) )
               zlim4    = trb(ji,jj,K,jpdoc) / ( xkdoc   + trb(ji,jj,K,jpdoc) )
               xlimbacl(ji,jj,jk) = MIN( zlim1, zlim2, zlim3 )
               xlimbac (ji,jj,jk) = xlimbacl(ji,jj,jk) * zlim4
               !
               ! Michaelis-Menten Limitation term for nutrients Small flagellates
               ! -----------------------------------------------
               zfalim = (1.-fanano) / fanano
               xnanonh4(ji,jj,jk) = (1. - fanano) * trb(ji,jj,K,jpnh4)   &
                  &                / ( zfalim * zconc0nnh4 + trb(ji,jj,K,jpnh4) )
               xnanono3(ji,jj,jk) = (1. - fanano) * trb(ji,jj,K,jpno3)   &
               &                    / ( zfalim * zconc0n + trb(ji,jj,K,jpno3) )  &
               &                    * (1. - xnanonh4(ji,jj,jk))
               !
               zfalim = (1.-fananop) / fananop
               xnanopo4(ji,jj,jk) = (1. - fananop) * trb(ji,jj,K,jppo4)  &
               &                   / ( trb(ji,jj,K,jppo4) + zfalim * zconc0npo4 )
               xnanodop(ji,jj,jk) = trb(ji,jj,K,jpdop)   &
               &                    / ( trb(ji,jj,K,jpdop) + xkdoc )   &
               &                    * ( 1.0 - xnanopo4(ji,jj,jk) )
               xnanodop(ji,jj,jk) = 0.
               !
               zfalim = (1.-fananof) / fananof
               xnanofer(ji,jj,jk) = (1. - fananof) * biron(ji,jj,jk) / ( biron(ji,jj,jk) + zfalim * zconcnfe )
               !
               zratiof   = trb(ji,jj,K,jpnfe) * z1_trnphy
               zqfemn = xcoef1 * znanochl + xcoef2 + xcoef3 * xnanono3(ji,jj,jk)
               !
               zration = trb(ji,jj,K,jpnph) * z1_trnphy
               zration = MIN(xqnnmax(ji,jj,jk), MAX( 2. * xqnnmin(ji,jj,jk), zration ))
               fvnuptk(ji,jj,jk) = 1. / zpsiuptk * rno3 * 2. * xqnnmin(ji,jj,jk) / (zration + rtrn)  &
               &                   * MAX(0., (1. - zratchl * znanochl / 12. ) )
               !
               zlim1    = max(0., (zration - 2. * xqnnmin(ji,jj,jk) )  &
               &          / (xqnnmax(ji,jj,jk) - 2. * xqnnmin(ji,jj,jk) ) ) * xqnnmax(ji,jj,jk)  &
               &          / (zration + rtrn)
               zlim3    = MAX( 0.,( zratiof - zqfemn ) / qfnopt ) 
               xlimnfe(ji,jj,jk) = MIN( 1., zlim3 )
               xlimphy(ji,jj,jk) = MIN( 1., zlim1, zlim3 )
               !
               ! Michaelis-Menten Limitation term for nutrients picophytoplankton
               ! ----------------------------------------------------------------
               zfalim = (1.-fapico) / fapico 
               xpiconh4(ji,jj,jk) = (1. - fapico) * trb(ji,jj,K,jpnh4)   &
               &                   / ( zfalim * zconc0pnh4 + trb(ji,jj,K,jpnh4) )
               xpicono3(ji,jj,jk) = (1. - fapico) * trb(ji,jj,K,jpno3)            &
               &                    / ( zfalim * zconc0p + trb(ji,jj,K,jpno3) )   &
               &                    * (1. - xpiconh4(ji,jj,jk))
               !
               zfalim = (1.-fapicop) / fapicop 
               xpicopo4(ji,jj,jk) = (1. - fapicop) * trb(ji,jj,K,jppo4)   &
               &                   / ( trb(ji,jj,K,jppo4) + zfalim * zconc0ppo4 )
               xpicodop(ji,jj,jk) = trb(ji,jj,K,jpdop)   &
               &                    / ( trb(ji,jj,K,jpdop) + xkdoc )   &
               &                    * ( 1.0 - xpicopo4(ji,jj,jk) )
               xpicodop(ji,jj,jk) = 0.
               !
               zfalim = (1.-fapicof) / fapicof
               xpicofer(ji,jj,jk) = (1. - fapicof) * biron(ji,jj,jk) / ( biron(ji,jj,jk) + zfalim * zconcpfe )
               !
               zratiof   = trb(ji,jj,K,jppfe) * z1_trnpic
               zqfemp = xcoef1 * zpicochl + xcoef2 + xcoef3 * xpicono3(ji,jj,jk)
               !
               zration   = trb(ji,jj,K,jpnpi) * z1_trnpic
               zration = MIN(xqnpmax(ji,jj,jk), MAX( 2. * xqnpmin(ji,jj,jk), zration ))
               fvpuptk(ji,jj,jk) = 1. / zpsiuptk * rno3 * 2. * xqnpmin(ji,jj,jk) / (zration + rtrn)  &
               &                   * MAX(0., (1. - zratchl * zpicochl / 12. ) ) 
               !
               zlim1    = max(0., (zration - 2. * xqnpmin(ji,jj,jk) )  &
               &          / (xqnpmax(ji,jj,jk) - 2. * xqnpmin(ji,jj,jk) ) ) * xqnpmax(ji,jj,jk)  &
               &          / (zration + rtrn)
               zlim3    = MAX( 0.,( zratiof - zqfemp ) / qfpopt )
               xlimpfe(ji,jj,jk) = MIN( 1., zlim3 )
               xlimpic(ji,jj,jk) = MIN( 1., zlim1, zlim3 )
               !
               !   Michaelis-Menten Limitation term for nutrients Diatoms
               !   ------------------------------------------------------
               zfalim = (1.-fadiat) / fadiat 
               xdiatnh4(ji,jj,jk) = (1. - fadiat) * trb(ji,jj,K,jpnh4)   &
               &                   / ( zfalim * zconc1dnh4 + trb(ji,jj,K,jpnh4) )
               xdiatno3(ji,jj,jk) = (1. - fadiat) * trb(ji,jj,K,jpno3)   &
               &                    / ( zfalim * zconc1d + trb(ji,jj,K,jpno3) )  &
               &                    * (1. - xdiatnh4(ji,jj,jk))
               !
               zfalim = (1.-fadiatp) / fadiatp
               xdiatpo4(ji,jj,jk) = (1. - fadiatp) * trb(ji,jj,K,jppo4)   &
               &                    / ( trb(ji,jj,K,jppo4) + zfalim * zconc0dpo4 )
               xdiatdop(ji,jj,jk) = trb(ji,jj,K,jpdop)  &
               &                    / ( trb(ji,jj,K,jpdop) + xkdoc )  &
               &                    * ( 1.0 - xdiatpo4(ji,jj,jk) )
               xdiatdop(ji,jj,jk) = 0.
               !
               zfalim = (1.-fadiatf) / fadiatf
               xdiatfer(ji,jj,jk) = (1. - fadiatf) * biron(ji,jj,jk) / ( biron(ji,jj,jk) + zfalim * zconcdfe )
               !
               zratiof   = trb(ji,jj,K,jpdfe) * z1_trndia
               zqfemd = xcoef1 * zdiatchl + xcoef2 + xcoef3 * xdiatno3(ji,jj,jk)
               !
               zration   = trb(ji,jj,K,jpndi) * z1_trndia
               zration = MIN(xqndmax(ji,jj,jk), MAX( 2. * xqndmin(ji,jj,jk), zration ))
               fvduptk(ji,jj,jk) = 1. / zpsiuptk * rno3 * 2. * xqndmin(ji,jj,jk) / (zration + rtrn)   &
               &                   * MAX(0., (1. - zratchl * zdiatchl / 12. ) ) 
               !
               zlim1    = max(0., (zration - 2. * xqndmin(ji,jj,jk) )    &
               &          / (xqndmax(ji,jj,jk) - 2. * xqndmin(ji,jj,jk) ) )   &
               &          * xqndmax(ji,jj,jk) / (zration + rtrn)
               zlim3    = trb(ji,jj,K,jpsil) / ( trb(ji,jj,K,jpsil) + xksi(ji,jj) )
               zlim4    = MAX( 0., ( zratiof - zqfemd ) / qfdopt )
               xlimdfe(ji,jj,jk) = MIN( 1., zlim4 )
               xlimdia(ji,jj,jk) = MIN( 1., zlim1, zlim3, zlim4 )
               xlimsi(ji,jj,jk)  = MIN( zlim1, zlim4 )
            END DO
         END DO
      END DO
      !
      ! Compute the phosphorus quota values. It is based on Litchmann et al., 2004 and Daines et al, 2013.
      ! The relative contribution of three fonctional pools are computed: light harvesting apparatus, 
      ! nutrient uptake pool and assembly machinery. DNA is assumed to represent 1% of the dry mass of 
      ! phytoplankton (see Daines et al., 2013). 
      ! --------------------------------------------------------------------------------------------------
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               ! Size estimation of nanophytoplankton
               ! ------------------------------------
               zfvn = 2. * fvnuptk(ji,jj,jk)
               sizen(ji,jj,jk) = MAX(1., MIN(xsizern, 1.0 / ( MAX(rtrn, zfvn) ) ) )

               ! N/P ratio of nanophytoplankton
               ! ------------------------------
               zfuptk = 0.23 * zfvn
               zrpho = 2.24 * trb(ji,jj,K,jpnch)   &
                  &  / ( trb(ji,jj,K,jpnph) * rno3 * 15. + rtrn )
               zrass = 1. - 0.2 - zrpho - zfuptk
               xqpnmax(ji,jj,jk) = ( zfuptk + zrpho ) * 0.0128 * 16. + zrass * 1./ 7.2 * 16.
               xqpnmax(ji,jj,jk) = xqpnmax(ji,jj,jk) * trb(ji,jj,K,jpnph)   &
                 &               / ( trb(ji,jj,K,jpphy) + rtrn ) + 0.13
               xqpnmin(ji,jj,jk) = 0.13 + 0.23 * 0.0128 * 16.

               ! Size estimation of picophytoplankton
               ! ------------------------------------
               zfvn = 2. * fvpuptk(ji,jj,jk)
               sizep(ji,jj,jk) = MAX(1., MIN(xsizerp, 1.0 / ( MAX(rtrn, zfvn) ) ) )

               ! N/P ratio of picophytoplankton
               ! ------------------------------
               zfuptk = 0.35 * zfvn
               zrpho = 2.24 * trb(ji,jj,K,jppch)   &
               &      / ( trb(ji,jj,K,jpnpi) * rno3 * 15. + rtrn )
               zrass = 1. - 0.4 - zrpho - zfuptk
               xqppmax(ji,jj,jk) =  (zrpho + zfuptk) * 0.0128 * 16. + zrass * 1./ 9. * 16.
               xqppmax(ji,jj,jk) = xqppmax(ji,jj,jk) * trb(ji,jj,K,jpnpi)   &
                  &              / ( trb(ji,jj,K,jppic) + rtrn ) + 0.13
               xqppmin(ji,jj,jk) = 0.13

               ! Size estimation of diatoms
               ! --------------------------
               zfvn = 2. * fvduptk(ji,jj,jk)
               sized(ji,jj,jk) = MAX(1., MIN(xsizerd, 1.0 / ( MAX(rtrn, zfvn) ) ) )
               zcoef = trb(ji,jj,K,jpdia) - MIN(xsizedia, trb(ji,jj,K,jpdia) )
               sized(ji,jj,jk) = 1. + xsizerd * zcoef *1E6 / ( 1. + zcoef * 1E6 )

               ! N/P ratio of diatoms
               ! --------------------
               zfuptk = 0.2 * zfvn
               zrpho = 2.24 * trb(ji,jj,K,jpdch) / ( trb(ji,jj,K,jpndi)   &
               &      * rno3 * 15. + rtrn )
               zrass = 1. - 0.2 - zrpho - zfuptk
               xqpdmax(ji,jj,jk) = ( zfuptk + zrpho ) * 0.0128 * 16. + zrass * 1./ 7.2 * 16.
               xqpdmax(ji,jj,jk) = xqpdmax(ji,jj,jk) * trb(ji,jj,K,jpndi)   &
                  &              / ( trb(ji,jj,K,jpdia) + rtrn ) + 0.13
               xqpdmin(ji,jj,jk) = 0.13 + 0.2 * 0.0128 * 16.

            END DO
         END DO
      END DO

      ! Compute the fraction of nanophytoplankton that is made of calcifiers
      ! --------------------------------------------------------------------
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zlim1 =  trb(ji,jj,K,jpnh4) / ( trb(ji,jj,K,jpnh4) + concnnh4 )   &
               &        + trb(ji,jj,K,jpno3)    &
               &        / ( trb(ji,jj,K,jpno3) + concnno3 )   &
               &        * ( 1.0 - trb(ji,jj,K,jpnh4)   &
               &        / ( trb(ji,jj,K,jpnh4) + concnnh4 ) )
               zlim2  = trb(ji,jj,K,jppo4) / ( trb(ji,jj,K,jppo4) + concnpo4 )
               zlim3  = trb(ji,jj,K,jpfer) / ( trb(ji,jj,K,jpfer) +  5.E-11 ) 
               ztem1  = MAX( 0., tsn(ji,jj,K,jp_tem) )
               ztem2  = tsn(ji,jj,K,jp_tem) - 10.
               zetot1 = MAX( 0., etot(ji,jj,jk) - 1.) / ( 4. + etot(ji,jj,jk) ) * 20. / ( 20. + etot(ji,jj,jk) ) 

!               xfracal(ji,jj,jk) = caco3r * MIN( zlim1, zlim2, zlim3 )                  &
               xfracal(ji,jj,jk) = caco3r                 &
               &                   * ztem1 / ( 1. + ztem1 ) * MAX( 1., trb(ji,jj,K,jpphy)*1E6 )   &
                  &                * ( 1. + EXP(-ztem2 * ztem2 / 25. ) )         &
                  &                * zetot1 * MIN( 1., 50. / ( hmld(ji,jj) + rtrn ) )
               xfracal(ji,jj,jk) = MAX( 0.02, MIN( 0.8 , xfracal(ji,jj,jk) ) )
            END DO
         END DO
      END DO
      !
      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               ! denitrification factor computed from O2 levels
               nitrfac(ji,jj,jk) = MAX(  0.e0, 0.4 * ( 6.e-6  - trb(ji,jj,K,jpoxy) )    &
                  &                                / ( oxymin + trb(ji,jj,K,jpoxy) )  )
               nitrfac(ji,jj,jk) = MIN( 1., nitrfac(ji,jj,jk) )
            END DO
         END DO
      END DO
      !
#if defined key_iomput
      IF( lk_iomput .AND. knt == nrdttrc ) THEN        ! save output diagnostics
        IF( iom_use( "xfracal" ) ) CALL iom_put( "xfracal", xfracal(:,:,:) * tmask(:,:,:) )  ! euphotic layer deptht
        IF( iom_use( "LNnut"   ) ) CALL iom_put( "LNnut"  , xlimphy(:,:,:) * tmask(:,:,:) )  ! Nutrient limitation term
        IF( iom_use( "LPnut"   ) ) CALL iom_put( "LPnut"  , xlimpic(:,:,:) * tmask(:,:,:) )  ! Nutrient limitation term
        IF( iom_use( "LDnut"   ) ) CALL iom_put( "LDnut"  , xlimdia(:,:,:) * tmask(:,:,:) )  ! Nutrient limitation term
        IF( iom_use( "LNFe"    ) ) CALL iom_put( "LNFe"   , xlimnfe(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "LPFe"    ) ) CALL iom_put( "LPFe"   , xlimpfe(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "LDFe"    ) ) CALL iom_put( "LDFe"   , xlimdfe(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "SIZEN"   ) ) CALL iom_put( "SIZEN"  , sizen(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "SIZEP"   ) ) CALL iom_put( "SIZEP"  , sizep(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "SIZED"   ) ) CALL iom_put( "SIZED"  , sized(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
      ENDIF
#endif
      !
   END SUBROUTINE p5z_lim


   SUBROUTINE p5z_lim_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p5z_lim_init  ***
      !!
      !! ** Purpose :   Initialization of nutrient limitation parameters
      !!
      !! ** Method  :   Read the nampislim and nampisquota namelists and check
      !!      the parameters called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampislim
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!
      NAMELIST/namp5zlim/ concnno3, concpno3, concdno3, concnnh4, concpnh4, concdnh4,  &
         &                concnfer, concpfer, concdfer, concbfe, concnpo4, concppo4,   &
         &                concdpo4, concbno3, concbnh4, concbpo4, xsizedia, xsizepic,  &
         &                xsizephy, xsizern, xsizerp, xsizerd, xksi1, xksi2, xkdoc,    &
         &                caco3r, oxymin
         !
      NAMELIST/namp5zquota/ qnnmin, qnnmax, qpnmin, qpnmax, qnpmin, qnpmax, qppmin,      &
         &                  qppmax, qndmin, qndmax, qpdmin, qpdmax, qfnmax, qfpmax, qfdmax,  &
         &                  qfnopt, qfpopt, qfdopt
      !!----------------------------------------------------------------------
      !
      REWIND( numnatp_ref )              ! Namelist nampislim in reference namelist : Pisces nutrient limitation parameters
      READ  ( numnatp_ref, namp5zlim, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampislim in reference namelist', lwp )
      !
      REWIND( numnatp_cfg )              ! Namelist nampislim in configuration namelist : Pisces nutrient limitation parameters 
      READ  ( numnatp_cfg, namp5zlim, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'nampislim in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, namp5zlim )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for nutrient limitations, namp5zlim'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    mean rainratio                           caco3r    = ', caco3r
         WRITE(numout,*) '    NO3 half saturation of nanophyto         concnno3  = ', concnno3
         WRITE(numout,*) '    NO3 half saturation of picophyto         concpno3  = ', concpno3
         WRITE(numout,*) '    NO3 half saturation of diatoms           concdno3  = ', concdno3
         WRITE(numout,*) '    NH4 half saturation for phyto            concnnh4  = ', concnnh4
         WRITE(numout,*) '    NH4 half saturation for pico             concpnh4  = ', concpnh4
         WRITE(numout,*) '    NH4 half saturation for diatoms          concdnh4  = ', concdnh4
         WRITE(numout,*) '    PO4 half saturation for phyto            concnpo4  = ', concnpo4
         WRITE(numout,*) '    PO4 half saturation for pico             concppo4  = ', concppo4
         WRITE(numout,*) '    PO4 half saturation for diatoms          concdpo4  = ', concdpo4
         WRITE(numout,*) '    half saturation constant for Si uptake   xksi1     = ', xksi1
         WRITE(numout,*) '    half saturation constant for Si/C        xksi2     = ', xksi2
         WRITE(numout,*) '    half-sat. of DOC remineralization        xkdoc     = ', xkdoc
         WRITE(numout,*) '    Iron half saturation for nanophyto       concnfer  = ', concnfer
         WRITE(numout,*) '    Iron half saturation for picophyto       concpfer  = ', concpfer
         WRITE(numout,*) '    Iron half saturation for diatoms         concdfer  = ', concdfer
         WRITE(numout,*) '    size ratio for nanophytoplankton         xsizern   = ', xsizern
         WRITE(numout,*) '    size ratio for picophytoplankton         xsizerp   = ', xsizerp
         WRITE(numout,*) '    size ratio for diatoms                   xsizerd   = ', xsizerd
         WRITE(numout,*) '    NO3 half saturation of bacteria          concbno3  = ', concbno3
         WRITE(numout,*) '    NH4 half saturation for bacteria         concbnh4  = ', concbnh4
         WRITE(numout,*) '    Minimum size criteria for diatoms        xsizedia  = ', xsizedia
         WRITE(numout,*) '    Minimum size criteria for picophyto      xsizepic  = ', xsizepic
         WRITE(numout,*) '    Minimum size criteria for nanophyto      xsizephy  = ', xsizephy
         WRITE(numout,*) '    Fe half saturation for bacteria          concbfe   = ', concbfe
         WRITE(numout,*) '    halk saturation constant for anoxia       oxymin   =' , oxymin
      ENDIF

      REWIND( numnatp_ref )              ! Namelist nampislim in reference namelist : Pisces nutrient limitation parameters
      READ  ( numnatp_ref, namp5zquota, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisquota in reference namelist', lwp )
      !
      REWIND( numnatp_cfg )              ! Namelist nampislim in configuration namelist : Pisces nutrient limitation parameters 
      READ  ( numnatp_cfg, namp5zquota, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 ) CALL ctl_nam ( ios , 'nampisquota in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, namp5zquota )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for nutrient limitations, namp5zquota'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    optimal Fe quota for nano.               qfnopt    = ', qfnopt
         WRITE(numout,*) '    optimal Fe quota for pico.               qfpopt    = ', qfpopt
         WRITE(numout,*) '    Optimal Fe quota for diatoms             qfdopt    = ', qfdopt
         WRITE(numout,*) '    Minimal N quota for nano                 qnnmin    = ', qnnmin
         WRITE(numout,*) '    Maximal N quota for nano                 qnnmax    = ', qnnmax
         WRITE(numout,*) '    Minimal P quota for nano                 qpnmin    = ', qpnmin
         WRITE(numout,*) '    Maximal P quota for nano                 qpnmax    = ', qpnmax
         WRITE(numout,*) '    Minimal N quota for pico                 qnpmin    = ', qnpmin
         WRITE(numout,*) '    Maximal N quota for pico                 qnpmax    = ', qnpmax
         WRITE(numout,*) '    Minimal P quota for pico                 qppmin    = ', qppmin
         WRITE(numout,*) '    Maximal P quota for pico                 qppmax    = ', qppmax
         WRITE(numout,*) '    Minimal N quota for diatoms              qndmin    = ', qndmin
         WRITE(numout,*) '    Maximal N quota for diatoms              qndmax    = ', qndmax
         WRITE(numout,*) '    Minimal P quota for diatoms              qpdmin    = ', qpdmin
         WRITE(numout,*) '    Maximal P quota for diatoms              qpdmax    = ', qpdmax
         WRITE(numout,*) '    Maximal Fe quota for nanophyto.          qfnmax    = ', qfnmax
         WRITE(numout,*) '    Maximal Fe quota for picophyto.          qfpmax    = ', qfpmax
         WRITE(numout,*) '    Maximal Fe quota for diatoms             qfdmax    = ', qfdmax
      ENDIF
      !
      zpsino3 = 2.3 * rno3
      zpsinh4 = 1.8 * rno3
      zpsiuptk = 2.3 * rno3
      !
      nitrfac (:,:,:) = 0.
      !
   END SUBROUTINE p5z_lim_init


   INTEGER FUNCTION p5z_lim_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_lim_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER ::   ierr(2)        ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      !*  Biological arrays for phytoplankton growth
      ALLOCATE( xpicono3(PRIV_3D_BIOARRAY), xpiconh4(PRIV_3D_BIOARRAY),       &
         &      xpicopo4(PRIV_3D_BIOARRAY), xpicodop(PRIV_3D_BIOARRAY),       &
         &      xnanodop(PRIV_3D_BIOARRAY), xdiatdop(PRIV_3D_BIOARRAY),       &
         &      xnanofer(PRIV_3D_BIOARRAY), xdiatfer(PRIV_3D_BIOARRAY),       &
         &      xpicofer(PRIV_3D_BIOARRAY), xlimpfe (PRIV_3D_BIOARRAY),       &
         &      fvnuptk (PRIV_3D_BIOARRAY), fvduptk (PRIV_3D_BIOARRAY),       &
         &      fvpuptk (PRIV_3D_BIOARRAY), xlimpic (PRIV_3D_BIOARRAY),    STAT=ierr(1) )
         !
      !*  Minimum/maximum quotas of phytoplankton
      ALLOCATE( xqnnmin (PRIV_3D_BIOARRAY), xqnnmax(PRIV_3D_BIOARRAY),       &
         &      xqpnmin (PRIV_3D_BIOARRAY), xqpnmax(PRIV_3D_BIOARRAY),       &
         &      xqnpmin (PRIV_3D_BIOARRAY), xqnpmax(PRIV_3D_BIOARRAY),       &
         &      xqppmin (PRIV_3D_BIOARRAY), xqppmax(PRIV_3D_BIOARRAY),       &
         &      xqndmin (PRIV_3D_BIOARRAY), xqndmax(PRIV_3D_BIOARRAY),       &
         &      xqpdmin (PRIV_3D_BIOARRAY), xqpdmax(PRIV_3D_BIOARRAY),     STAT=ierr(2) )
         !
      p5z_lim_alloc = MAXVAL( ierr )
      !
      IF( p5z_lim_alloc /= 0 ) CALL ctl_warn( 'p5z_lim_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p5z_lim_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p5z_lim                   ! Empty routine
   END SUBROUTINE p5z_lim
#endif

   !!======================================================================
END MODULE p5zlim
