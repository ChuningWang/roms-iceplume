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
#include "cppdefs.h"
#ifdef BULK_FLUX
!======================================================================
      subroutine bulk_flux (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call bulk_flux_tile (Istr,Iend,Jstr,Jend,
     &                               A2d(1,1,trd),A2d(1,2,trd))
      return
      end

      subroutine bulk_flux_tile (Istr,Iend,Jstr,Jend, aer,cer)
      !
      !======================================================================
      !                 ***  SUBROUTINE  bulk_flux  ***
      !  This routine computes the turbulent and radiative components of
      !  air-sea fluxes for the specification of surface boundary conditions
      !=======================================================================
      !
      implicit none
#include "param.h"
#include "grid.h"
#include "ocean3d.h"
#include "forces.h"
#include "scalars.h"
#include "params_bulk.h"
!++ Local integers
      INTEGER              :: i, j, m, iter
      INTEGER              :: Istr, Iend, Jstr, Jend
      INTEGER              :: imin, imax, jmin, jmax
      INTEGER, PARAMETER   :: IterFl = 3
#ifdef CFB_WIND_TRA
      INTEGER, PARAMETER :: mb = 2
#else
      INTEGER, PARAMETER :: mb = 1
#endif

!++ Local scalars
      REAL           :: cff  , cff2
      REAL           :: rho0i, cpi  , RH   , Hlv     , wspd0
      REAL           :: TseaC, TseaK, Qsea , spec_hum, Tvstar
      REAL           :: TairC, TairK, Q    , rhoAir  , Bf
      REAL           :: delW(2) , delT , delQ , ddelW(2) , ddelT , ddelQ
      REAL           :: Cd   , Ch   , Ce   , qsat    , iexns , iexna
      REAL           :: Wstar(2), Tstar, Qstar, Wstar0  , Tstar0, Qstar0
      REAL           :: psi_u(2), ZoLu(2) , psi_t, ZoLt    , patm
      REAL           :: hfsen, hflat, hflw , upvel   , evap
# ifdef CFB_WIND_TRA
      REAL           :: wspd0_cfb
# endif
# ifdef BULK_ECUMEV0
      REAL            :: Cdn(2), Chn, Cen
# endif
# ifdef BULK_ECUMEV6
      REAL            :: partn, parqn
# endif
# ifdef BULK_WASP
      INTEGER         :: idx
      REAL            :: chi,aa,bb
# endif
!++ Parameter values (blk_ZW and blk_ZT should not be parameters with ABL)
      REAL, PARAMETER ::  blk_ZW =  10.0
      REAL, PARAMETER ::  blk_ZT =   2.0
      REAL, PARAMETER ::  blk_ZToZW = blk_ZT/blk_ZW
      REAL            ::  logus10,logts10
# if defined BULK_ECUMEV0 || defined BULK_ECUMEV6
      PARAMETER ( logus10   = LOG(blk_ZW/10.0) )
      PARAMETER ( logts10   = LOG(blk_ZT/10.0) )
      REAL           :: bulk_psiu_ecume, bulk_psit_ecume
      REAL           :: deltau10n(2), deltat10n, deltaq10n
# else
      REAL,PARAMETER :: Log10oLogZw = LOG(  10.0*10000.)/
     &                                LOG(blk_ZW*10000.)
      REAL           :: bulk_psiu_coare, bulk_psit_coare, air_visc
      REAL           :: iZoW , iZoT  , CC, Ch10, Ribcu , charn
      REAL           :: iZo10, iZoT10, Ri, Rr  , VisAir
# endif
!++ Local arrays
      REAL           ::  aer(PRIVATE_2D_SCRATCH_ARRAY)
      REAL           ::  cer(PRIVATE_2D_SCRATCH_ARRAY)
#ifdef CFB_STRESS
      REAL           ::  stau(PRIVATE_2D_SCRATCH_ARRAY)
#endif
!======================================================================
! Extended range  (needed for subgrid scale closures)
!======================================================================
# ifdef EW_PERIODIC
        imin=Istr-2
        imax=Iend+2
# else
      if (WESTERN_EDGE) then
        imin=Istr-1
      else
        imin=Istr-2
      endif
      if (EASTERN_EDGE) then
        imax=Iend+1
      else
        imax=Iend+2
      endif
# endif
# ifdef NS_PERIODIC
        jmin=Jstr-2
        jmax=Jend+2
# else
      if (SOUTHERN_EDGE) then
        jmin=Jstr-1
      else
        jmin=Jstr-2
      endif
      if (NORTHERN_EDGE) then
        jmax=Jend+1
      else
        jmax=Jend+2
      endif
# endif
!======================================================================
! Initialization of various constants
!======================================================================
!++  Inverse seawater density
      rho0i=1.0/rho0
!++  Inverse of specific heat for seawater (kg-degC/Joule) [cp is defined in scalars.h]
      cpi=1.0/cp
!++
#ifndef READ_PATM
!++  Inverse of exner function at the surface
      iexns  = (psurf*ip00)**(-rdocpd)
#endif
!----------------------
      DO j=jmin,jmax
        DO i=imin,imax
!----------------------
!======================================================================
! Initialization of bulk quantities (common to all bulk algorithms)
!======================================================================
          wspd0     = wspd(i,j)
          wspd0     = MAX ( wspd0 , 0.1 * MIN(10., blk_ZW) )
#ifdef CFB_WIND_TRA
          wspd0_cfb = wspd_cfb(i,j)
          wspd0_cfb = MAX ( wspd0_cfb , 0.1 * MIN(10., blk_ZW) )
#endif
          TairC = tair(i,j)
          TairK = TairC + CtoK
# ifdef SST_SKIN
          TseaC = sst_skin(i,j)
# else
          TseaC = t(i,j,N,nrhs,itemp)
# endif
          TseaK = TseaC + CtoK
!
# ifdef READ_PATM
!++  Get pressure from ATM file in Pa
          psurf  = patm2d(i,j)
!++  Inverse of exner function at the surface
          iexns  = (psurf*ip00)**(-rdocpd)
# endif
!
          RH     = rhum(i,j)
          Q      = spec_hum (RH,psurf,TairC)
          CALL exner_patm_from_tairabs(iexna,patm,Q,TairK,blk_ZT,psurf)
!++  Air density
          rhoAir = patm*(1.+Q) / ( blk_Rgas*TairK*(1.+MvoMa*Q) )
!++  Specific humidity at saturation
          Qsea   = qsat(TseaK,psurf,0.98)
!++  Air-sea gradients delW, delT, delQ through the interface
#ifdef BULK_GUSTINESS
          delW(1) = SQRT(wspd0*wspd0+0.25)
#  ifdef CFB_WIND_TRA
          delW(2) = SQRT(wspd0_cfb*wspd0_cfb+0.25)
#  endif
#else
          delW(1) = wspd0
#  ifdef CFB_WIND_TRA
          delW(2) = wspd0_cfb
#  endif
#endif
          delQ = Q-Qsea
          cff  = CtoK*(iexna-iexns)
          delT = TairC*iexna - TseaC*iexns + cff

# if defined BULK_ECUMEV0 || defined BULK_ECUMEV6
!
!======================================================================
!   @@@@@@@@@@  @@@@@@@@@@  @@      @@  @@@    @@@  @@@@@@@@@@
!   @@          @@          @@      @@  @@ @  @ @@  @@
!   @@@@@@@     @@          @@      @@  @@  @@  @@  @@@@@@@
!   @@          @@          @@      @@  @@      @@  @@
!   @@@@@@@@@@  @@@@@@@@@@   @@@@@@@@   @@      @@  @@@@@@@@@@
!======================================================================
!
!++ Initial guess for ECUME algorithms
          ddelW(1:mb)     = SIGN(MAX(ABS(delW(1:mb)),10.0*dWstar0),delW(1:mb))
          ddelT           = SIGN(MAX(ABS(delT),10.0*dTstar0),delT)
          ddelQ           = SIGN(MAX(ABS(delQ),10.0*dQstar0),delQ)
          deltau10n(1:mb) = ddelW(1:mb)
          deltat10n       = ddelT
          deltaq10n       = ddelQ
!======================================================================
!  ITERATIVE LOOP TO COMPUTE U*, T*, Q*.
!======================================================================
!----------------------
          DO Iter=1,IterFl
!----------------------
#ifdef BULK_ECUMEV0
            DO m=1,mb
              cff    = deltau10n(m)*deltau10n(m)
!++  Neutral coefficient for wind speed cdn (ECUME_v0 formulation)
              IF ( deltau10n(m) <= utu1 ) THEN
                Cdn(m) = coefu10+coefu11*deltau10n(m) +
     &                   (coefu12+coefu13*deltau10n(m))*cff
              ELSEIF ( deltau10n(m) <= utu2 ) THEN
                Cdn(m) = coefu20+coefu21*deltau10n(m)
     &                 + (coefu22+coefu23*deltau10n(m))*cff
     &                 + coefu24 * cff * cff
              ELSE
                Cdn(m) = Cdn0
              ENDIF
            ENDDO
!++ Neutral coefficient for temperature chn (ECUME_v0 formulation)
            IF ( deltau10n(mb) <= utt ) THEN
              Chn =  coeft0+coeft1*deltau10n(mb)
     &            +  (coeft2+coeft3*deltau10n(mb))*cff
     &            +  (coeft4+coeft5*deltau10n(mb)) * cff * cff
            ELSE
              Chn = Chn0
            ENDIF
!++ Neutral coefficient for humidity cen (ECUME_v0 formulation)
            IF ( deltau10n(mb) <= utq1 ) THEN
              Cen = coefq10+coefq11*deltau10n(mb)
     &            + (coefq12+coefq13*deltau10n(mb))*cff
     &            + (coefq14                      )* cff* cff
            ELSEIF ( deltau10n(mb) <= utq2 ) then
              Cen = coefq20 + coefq21 * deltau10n(mb) +  ( coefq22 )*cff
            ELSE
              Cen = Cen0
            ENDIF
!++ Estimate Monin-Obukhov similarity scales (ECUME_v0 formulation)
            Wstar(1:mb) = SQRT(Cdn(1:mb))*deltau10n(1:mb)
            Tstar       = Chn/SQRT(Cdn(mb))*deltat10n
            Qstar       = Cen/SQRT(Cdn(mb))*deltaq10n
#else
            DO m=1,mb
              cff    = deltau10n(m)*deltau10n(m)
              cff2   = 1./deltau10n(m)
!++  Neutral parameter for wind speed (ECUME_V6 formulation)
              IF (deltau10n(m) <= utu) THEN
                Wstar(m) = (coefu0+coefu1*deltau10n(m))
     &                   + (coefu2+coefu3*deltau10n(m))*cff
     &                   + (coefu4+coefu5*deltau10n(m))*cff*cff
              ELSE
                Wstar(m) = cdiru*(deltau10n(m)-utu) + ordou
              ENDIF
            ENDDO
!++  Neutral parameter for temperature (ECUME_V6 formulation)
            IF (deltau10n(mb) <= utt) then
              partn = (coeft0+coeft1*deltau10n(mb))
     &              + (coeft2+coeft3*deltau10n(mb))*cff
     &              +  coeft4*cff*cff
            ELSE
              partn = cdirt*(deltau10n(mb)-utt) + ordot
            ENDIF
            Tstar = partn*deltat10n*cff2
!++  Neutral parameter for humidity (ECUME_V6 formulation)
            IF (deltau10n(mb) <= utq) then
              parqn =   coefq0+coefq1*deltau10n(mb)+coefq2*cff
            ELSE
              parqn = cdirq*(deltau10n(mb)-utq) + ordoq
            ENDIF
            Qstar = parqn*deltaq10n*cff2
#endif
!++ Obukhovs stability param. z/l
            DO m=1,mb
              ZoLu(m)  = blk_ZW*g*vonKar*(Tstar*(1.+cpvir*Q)
     &                   +cpvir*TairK*Qstar)/
     &                   (TairK*Wstar(m)*Wstar(m)*(1.+cpvir*Q)+eps)
              ZoLu(m)  = MAX(MIN(ZoLu(m),LMOmax),LMOmin)
              psi_u(m) = bulk_psiu_ecume ( ZoLu(m) )
            ENDDO
            ZoLt = ZoLu(mb)*blk_ZToZW
            ZoLt = MAX(MIN(ZoLt,LMOmax),LMOmin)
!++ Stability function psi (see Liu et al, 1979 ; Dyer and Hicks, 1970)
!++  Modified to include convective form following Fairall (unpublished)
            psi_t = bulk_psit_ecume ( ZoLt )
!
#ifdef BULK_GUSTINESS
!++ Add gustiness contribution to wind speed in unstable conditions
            Tvstar   = Tstar*(1.0+cpvir*Q)+cpvir*TairK*Qstar
            Bf       = MAX(0.0,-g/TairK*Wstar(1)*Tvstar)
            cff2     = blk_beta*blk_beta*(Bf*blk_Zabl)**(2./3.)
            delW (1) = SQRT(wspd0*wspd0+cff2)
            ddelW(1) = SIGN(MAX(ABS(delW(1)),10.0*dWstar0),delW(1))
#  ifdef CFB_WIND_TRA
            Bf       = MAX(0.0,-g/TairK*Wstar(2)*Tvstar)
            cff2     = blk_beta*blk_beta*(Bf*blk_Zabl)**(2./3.)
            delW (2) = SQRT(wspd0_cfb*wspd0_cfb+cff2)
            ddelW(2) = SIGN(MAX(ABS(delW(2)),10.0*dWstar0),delW(2))
#  endif
#endif
!++ Update 10 meter Neutral quantities
            cff          = 1./vonKar
            cff2         = ddelW(1)-Wstar(1)*(logus10-psi_u(1))*cff
            deltau10n(1) = SIGN(MAX(ABS(cff2),10.0*dWstar0),cff2)
#ifdef CFB_WIND_TRA
            cff2         = ddelW(2)-Wstar(2)*(logus10-psi_u(2))*cff
            deltau10n(2) = SIGN(MAX(ABS(cff2),10.0*dWstar0),cff2)
#endif
            cff2         = ddelT-Tstar*(logts10-psi_t)*cff
            deltat10n    = SIGN(MAX(ABS(cff2),10.0*dTstar0),cff2)
            cff2         = ddelQ-Qstar*(logts10-psi_t)*cff
            deltaq10n    = SIGN(MAX(ABS(cff2),10.0*dQstar0),cff2)
!----------------------
      ENDDO  !<-- terminate iterations-loop
!----------------------
#else
!
!======================================================================
!   @@@@@@@@@@   @@@@@@@@@@   @@@@@@@@@@  @@@@@@@@@@  @@@@@@@@@@
!   @@          @@        @@  @@      @@  @@      @@  @@
!   @@          @@        @@  @@@@@@@@@@  @@@@@@@@@@  @@@@@@@
!   @@          @@        @@  @@      @@  @@    @@    @@
!   @@@@@@@@@@   @@@@@@@@@@   @@      @@  @@      @@  @@@@@@@@@@
!
!                                     @@
!                                   @@
!                                 @@
!                               @@
!                             @@
!
!           @@      @@  @@@@@@@@@@  @@@@@@@@@@  @@@@@@@@@@
!           @@      @@  @@      @@  @@          @@      @@
!           @@  @@  @@  @@@@@@@@@@  @@@@@@@@@@  @@@@@@@@@@
!            @@ @@ @@   @@      @@          @@  @@
!             @@  @@    @@      @@  @@@@@@@@@@  @@
!======================================================================
!
!++ Initial guess for COARE algorithm
          Wstar(1:mb)  = 0.035*delW(1:mb)*Log10oLogZw   ! Friction velocity
          VisAir = air_visc ( TairC )       ! Molecular viscosity
          charn  = 0.011                    ! Charnock parameter
          Ch10   = 0.00115                                                ! 10m exchange coefficient for heat
          Ribcu  = -    blk_ZW / ( blk_Zabl * 0.004 * blk_beta**3 )
          DO m=1,mb
            iZo10  = g*Wstar(m) /
     &               (charn*Wstar(m)*Wstar(m)*Wstar(m)+0.11*g*VisAir)    ! inverse of roughness length
            iZoT10 = 0.1 * exp(vonKar*vonKar /
     &                     ( Ch10*LOG( 10.0*iZo10 ) ) ) ! inverse of thermal roughness length
!++ Obukhovs stability param. z/l
            CC     = LOG( blk_ZW*iZo10 )*LOG( blk_ZW*iZo10 )/
     &               LOG( blk_ZT*iZoT10 )
            Ri     =  g * blk_ZW * ( delT+cpvir*TairK*delQ )/
     &                             ( TairK*delW(m)*delW(m) )
            IF ( Ri < 0.0 ) THEN
              ZoLu(m)=CC*Ri/(1.0+Ri/Ribcu)       ! Unstable
            ELSE
              ZoLu(m)=CC*Ri/(1.0+3.0*Ri/CC)      ! Stable
            ENDIF
            psi_u(m)     = bulk_psiu_coare ( ZoLu(m) )
            logus10   = LOG(blk_ZW* iZo10)
            Wstar(m)  = delW(m)*vonKar/(logus10-psi_u(m))
          ENDDO
          ZoLt      = ZoLu(mb)*blk_ZToZW        ! rescale in case blk_ZW != blk_ZT
!++ Stability functions
          psi_t     = bulk_psit_coare ( ZoLt    )
!++ Initial guess for Monin-Obukhov similarity scales
          logts10   = LOG(blk_ZT*iZoT10)
          cff       =      vonKar/(logts10-psi_t)
          Tstar     = delT*cff
          Qstar     = delQ*cff
!++  Compute Charnock coefficient
#  ifdef BULK_WASP
          IF ( wspd0 >= 7.0 .AND. wspd0 < 23.0 ) then
            idx = 2
          ELSEIF ( wspd0 >= 23.0 .AND. wspd0 < 25.0 ) then
            idx = 3
          ELSEIF ( wspd0 >= 25.0) THEN
            idx = 4
          ELSE
            idx = 1
          ENDIF
          chi   = CWage * ( wspd0 / Wstar(1) )   ! Wave age
          IF(idx==1) chi = wspd0
          cff   = wspd0 * wspd0
          aa    =     Awasp(0,idx) + Awasp(1,idx)*wspd0
     &        + cff*( Awasp(2,idx) + Awasp(3,idx)*wspd0 )
          bb    =     Bwasp(0,idx) + Bwasp(1,idx)*wspd0
     &        + cff*( Bwasp(2,idx) + Bwasp(3,idx)*wspd0 )
!
          charn = aa *(chi**bb)             ! Wind-dependent Charnock parameter
!++  Limit Charnock parameter
          IF(idx==3 .AND. charn < Charn0 ) charn = Charn0
          IF(idx==4 .AND. charn < Charn2 ) charn = Charn2
          IF(charn > Charn1) charn = Charn1
#  else
          charn = 0.011
          IF     ( delW(1) > 18.0 ) then
            charn = 0.018
          ELSEIF ( delW(1) > 10.0 ) then
            charn = 0.011+0.125*(0.018-0.011)*(delW(1)-10.)
          ENDIF
#  endif
!
!======================================================================
!  ITERATIVE LOOP TO COMPUTE U*, T*, Q*.
!======================================================================
!----------------------
          DO Iter=1,IterFl
!----------------------
            DO m=1,mb
!++  Inverse of roughness length
              iZoW    = g*Wstar(m) / ( charn*Wstar(m)*Wstar(m)*Wstar(m)
     &                                 +0.11*g*VisAir )
!++  Inverse of thermal roughness length
#  ifdef BULK_WASP
              iZoT    = EXP(vonKar*vonKar/(Ch10*LOG(blk_ZW*iZoW)) )/
     &                      blk_ZW
#  else
              Rr      = Wstar(m)/(iZow*VisAir)
              iZoT    = MAX(8695.65,18181.8*(Rr**0.6))
#  endif
!++  Obukhovs stability param. z/l
              ZoLu(m)    = vonKar*g*blk_ZW*
     &             (Tstar*(1.0+cpvir*Q)+cpvir*TairK*Qstar)/
     &             (TairK*Wstar(m)*Wstar(m)*(1.0+cpvir*Q)+eps)
!++ Stability functions
              psi_u(m)     = bulk_psiu_coare ( ZoLu(m) )
              logus10 = LOG(blk_ZW*iZoW)
              Wstar(m)   = delW(m)*vonKar/(logus10-psi_u(m))
            ENDDO

            ZoLt      = ZoLu(mb)*blk_ZToZW
            psi_t     = bulk_psit_coare ( ZoLt )

!++  Compute Monin-Obukhov similarity scales

            logts10   = LOG(blk_ZT*iZoT)
            cff       = vonKar/(logts10-psi_t)
            Tstar     = delT*cff
            Qstar     = delQ*cff
!++  Add BULK_GUSTINESS
# ifdef BULK_GUSTINESS
              Bf=-g/TairK*Wstar(1)*(Tstar+cpvir*TairK*Qstar)
              if (Bf.gt.0.0) then
                cff=blk_beta*(Bf*blk_Zabl)**r3
              else
                cff=0.2
              endif
              delW(1)  = SQRT(wspd0*wspd0+cff*cff)
#  ifdef CFB_WIND_TRA
              Bf=-g/TairK*Wstar(2)*(Tstar+cpvir*TairK*Qstar)
              if (Bf.gt.0.0) then
                cff=blk_beta*(Bf*blk_Zabl)**r3
              else
                cff=0.2
              endif
              delW(2)  = SQRT(wspd0_cfb*wspd0_cfb+cff*cff)
#  endif

# else
            delW(1)  = wspd0
#  ifdef CFB_WIND_TRA
            delW(2)  = wspd0_cfb
#  endif
# endif
            ddelW = delW
!----------------------
          ENDDO !<-- terminate iterations-loop
! ---------------------
#endif

!===============================================================================
!   @@@@@@@@   @@@@@@@@@@      @@@@@@@@@@ @@          @@        @@  @@     @@
!  @@      @@  @@      @@      @@         @@          @@        @@    @@  @@
!  @@      @@  @@@@@@@@@@  &&  @@@@@@@    @@          @@        @@      @@
!  @@      @@  @@      @@      @@         @@          @@        @@    @@  @@
!   @@@@@@@@   @@      @@      @@         @@@@@@@@@@   @@@@@@@@@@   @@      @@
!===============================================================================

!++  Compute turbulent sensible heat flux (W/m2).
          hfsen = - blk_Cpa*rhoAir*Wstar(mb)*Tstar
!++ Latent heat of vaporisation (J/kg)
          Hlv   = (2.5008 - 0.0023719*TseaC)*1.0E+6
!++  Compute turbulent latent heat flux (W/m2).
          hflat = - Hlv*rhoAir*Wstar(mb)*Qstar
!++ Longwave Radiation.
# ifndef BULK_LW
          hflw   = -radlw(i,j) ! positive downward, i.e., warming the ocean
# else
!++  Treat input longwave data as downwelling radiation only and add
!++  outgoing IR from model sea surface temperature.
          hflw   = radlw(i,j)  ! positive downward
     &           - emiss_lw*rho0i*cpi*SigmaSB*TseaK*TseaK*TseaK*TseaK
# endif
!++  Compute Webb correction (Webb effect) to latent heat flux, Hlw.
          upvel=-1.61*Wstar(mb)*Qstar-(1.0+1.61*Q)*Wstar(mb)
     &                                           *     Tstar/TairK
!++  Compute turbulent latent heat flux (W/m2).
          hflat=hflat+rhoAir*Hlv*upvel*Q
!
!======================================================================
!  Latent and sensible heat fluxes
!======================================================================
!
          hflat=-hflat*rho0i*cpi
          hfsen=-hfsen*rho0i*cpi
!++  Total surface heat flux
          stflx(i,j,itemp)=srflx(i,j)+hflw+hflat+hfsen
# ifdef SALINITY
!++  Total surface salt flux
          evap=-cp*hflat/Hlv
          stflx(i,j,isalt)=(evap-prate(i,j))*t(i,j,N,nrhs,isalt)
# endif
!
# ifdef MASKING
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask(i,j)
# endif
# ifdef WET_DRY
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask_wet(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask_wet(i,j)
# endif
!++  Save for the calculation of wind stress
          aer(i,j)  = rhoAir*delW(1)
          Cd        = (Wstar(1)/ddelW(1))**2
          cer(i,j)  = Cd
# ifdef CFB_STRESS
          stau(i,j) = cfb_slope * delW(1) + cfb_offset
# endif
!++  Store fluxes for diagnostics, when storing in netCDF file, multiplied by rho0.Cp to get W/m2
          shflx_rsw(i,j)=srflx(i,j)
          shflx_lat(i,j)=hflat
          shflx_sen(i,j)=hfsen
          shflx_rlw(i,j)=hflw
!----------------------
        ENDDO !<-- terminate i-loop
      ENDDO   !<-- terminate j-loop
!----------------------




!
!======================================================================
!  SWITCH OFF shortwave radiation over "sea-ice"
!======================================================================
#ifdef SEA_ICE_NOFLUX
      do j=jmin,jmax
        do i=imin,imax
              if( t(i,j,N,nrhs,itemp) .le. -1.8 ) then
                   stflx(i,j,itemp)=0.
#   if defined LMD_SKPP || defined LMD_BKPP || defined GLS_MIXING
                   srflx(i,j)=0.
#    endif
              endif
        enddo
      enddo
#endif
!
!======================================================================
!  Compute kinematic, surface wind stress (m2/s2).
!======================================================================
!
      do j=jmin,jmax
        do i=imin+1,imax
          cff =0.5*(cer(i-1,j)+cer(i,j))
          cff2=0.5*(aer(i-1,j)+aer(i,j))
          sustr(i,j)= ( cff*cff2*uwnd(i,j)
#ifdef CFB_STRESS
     &     + 0.5*(stau(i,j)+stau(i-1,j))*u(i,j,N,nrhs)
#endif
     &                             )*rho0i
#  ifdef MASKING
          sustr(i,j)=sustr(i,j)*umask(i,j)
#  endif
#  ifdef WET_DRY
          sustr(i,j)=sustr(i,j)*umask_wet(i,j)
#  endif
        enddo
      enddo
!
      do j=jmin+1,jmax
        do i=imin,imax
          cff =0.5*(cer(i,j-1)+cer(i,j))
          cff2=0.5*(aer(i,j-1)+aer(i,j))
          svstr(i,j)=( cff*cff2*vwnd(i,j)
#ifdef CFB_STRESS
     &   + 0.5*(stau(i,j)+stau(i,j-1))*v(i,j,N,nrhs)
#endif
     &                                        )*rho0i
#  ifdef MASKING
          svstr(i,j)=svstr(i,j)*vmask(i,j)
#  endif
#  ifdef WET_DRY
          svstr(i,j)=svstr(i,j)*vmask_wet(i,j)
#  endif
        enddo
      enddo

      return
      end


!======================================================================
!                     EXTERNAL FUNCTIONS
!======================================================================


!=======================================================================
!   Compute specific humidity from relative humidity
!=======================================================================
      FUNCTION spec_hum (RH,psfc,TairC)
      IMPLICIT NONE
#include "params_bulk.h"
      REAL       ::  RH   , spec_hum , cff
      REAL       ::  psfc, TairC
!++  Compute air saturation vapor pressure (mb), using Teten formula.
      cff=(1.0007+3.46e-6*0.01*psfc)*6.1121*
     &        exp(17.502*TairC/(240.97+TairC))
!++  Compute specific humidity, Q (kg/kg).
      IF (RH.lt.2.0) then                             ! RH fraction
         cff=cff*RH                                   ! Vapor pres (mb)
         spec_hum=MvoMa*(cff/(psfc*0.01-0.378*cff))  ! Spec hum (kg/kg)
      ELSE             !RH input was actually specific humidity in g/kg
         spec_hum=0.001*RH                            ! Spec Hum (kg/kg)
      ENDIF
      END FUNCTION spec_hum
!=======================================================================
!    Compute Exner function from absolute air temperature
!=======================================================================
      SUBROUTINE exner_patm_from_tairabs (iexn,pair,q,tairabs,z,psfc)
      IMPLICIT NONE
#include "params_bulk.h"
      REAL,INTENT(  out)   :: iexn,pair
      REAL,INTENT(in   )   :: q, tairabs, z, psfc
!++  Local variable declarations.
      REAL                 :: xm,q_sat,qsat
      REAL, PARAMETER      ::  g    = 9.80665
      INTEGER              :: iter
      INTEGER, PARAMETER   :: Niter = 3
!
      pair = psfc
      DO Iter = 1, Niter
        q_sat = qsat(tairabs, pair, 1.)
        xm    =  mm_dryair + (q/q_sat) * ( mm_water - mm_dryair )
        pair  = psfc * EXP( -g * xm * z / ( r_gas * tairabs ) )
      ENDDO
      iexn =  (pair*ip00)**(-rdocpd)
!
      return
      END SUBROUTINE exner_patm_from_tairabs
!=======================================================================
!    Compute Saturation humidity
!=======================================================================
      FUNCTION qsat (TairK, patm, coeff)
      IMPLICIT NONE
#include "params_bulk.h"
      REAL                 ::  qsat
      REAL                 ::  TairK, patm, coeff
!++  Local variables declarations.
      REAL                 ::  psat
      REAL, PARAMETER      ::  alpw    = 60.2227554
      REAL, PARAMETER      ::  betaw   = 6822.40088
      REAL, PARAMETER      ::  gamw    = 5.13926744
      REAL, PARAMETER      ::  alpi    = 32.62117980819471
      REAL, PARAMETER      ::  betai   = 6295.421338904806
      REAL, PARAMETER      ::  gami    = 0.5631331575423155
!++  Compute Saturation Vapor Pressure
      IF (TairK .LE. CtoK) then
        psat = EXP( alpi - betai/TairK - gami*LOG(TairK) )
      ELSE
        psat = EXP( alpw - betaw/TairK - gamw*LOG(TairK) )
      ENDIF
      psat = coeff * psat
!++ Compute Saturation Humidity
      qsat = (MvoMa*psat)/(patm+(MvoMa-1.0)*psat)
      return
      END FUNCTION qsat

# if defined BULK_ECUMEV0 || defined BULK_ECUMEV6
!=======================================================================
!    Ecume stability function for velocity
!=======================================================================
      FUNCTION bulk_psiu_ecume (ZoL)
      IMPLICIT NONE
#include "params_bulk.h"
      REAL                 ::  bulk_psiu_ecume
      REAL                 ::  ZoL
!++  Local variables declarations.
      REAL                 ::  chik, psik
      REAL                 ::  chic, psic
!
      IF (ZoL >= 0.0) THEN
        bulk_psiu_ecume = -7.0*ZoL
      ELSE
        chik  = (1.0-16.0*ZoL)**0.25
        psik  = 2.0*LOG(0.5*(1.0+chik))+LOG(0.5*(1.0+chik**2))
     &              -2.0*ATAN(chik)+pis2
        chic  = (1.0-12.87*ZoL)**r3    !for very unstable conditions
        psic  = 1.5*LOG(r3*(chic**2+chic+1.0))
     &               - sqr3*ATAN((2.0*chic+1.0)/sqr3)+2.*pis2osqr3
        bulk_psiu_ecume = psic+(psik-psic)/(1.0+ZoL**2) !match Kansas & free-conv. forms
      ENDIF
      return
      END FUNCTION bulk_psiu_ecume
!=======================================================================
!    Ecume stability function for tracers
!=======================================================================
      FUNCTION bulk_psit_ecume (ZoL)
      IMPLICIT NONE
#include "params_bulk.h"
      REAL                 ::  bulk_psit_ecume
      REAL                 ::  ZoL
      !  Local variables declarations.
      REAL                 ::  chik, psik
      REAL                 ::  chic, psic
      !
      IF (ZoL >= 0.0) THEN
        bulk_psit_ecume = -7.0*ZoL
      ELSE
        chik  = (1.0-16.0*ZoL)**0.25
        psik  = 2.0*LOG(0.5*(1.0+chik**2))
        chic  = (1.0-12.87*ZoL)**r3     !for very unstable conditions
        psic  = 1.5*LOG((chic**2+chic+1.0)*r3)
     &               -sqr3*ATAN((2.0*chic+1.0)/sqr3)
     &               +2.*pis2osqr3
        bulk_psit_ecume = psic+(psik-psic)/(1.0+ZoL**2)    !match Kansas & free-conv. forms
      ENDIF
      return
      END FUNCTION bulk_psit_ecume
!
# else
!=======================================================================
!    Compute molecular viscosity as a function of air temperature
!=======================================================================
      FUNCTION air_visc(TairC)
      REAL                 :: air_visc,cff
      REAL, PARAMETER      :: c0 = 1.326E-5
      REAL, PARAMETER      :: c1 = 6.542E-3
      REAL, PARAMETER      :: c2 = 8.301E-6
      REAL, PARAMETER      :: c3 = 4.84E-9
      cff      = TairC*TairC
      air_visc = c0*(1.+c1*TairC+c2*cff-c3*cff*TairC)
      return
      END FUNCTION air_visc
!=======================================================================
!    Coare stability function for velocity
!=======================================================================
      FUNCTION bulk_psiu_coare (ZoL)
      IMPLICIT NONE
#include "params_bulk.h"
      REAL                 ::  bulk_psiu_coare
      REAL                 ::  ZoL
!  Local variables declarations.
      REAL                 ::  chik, psik
      REAL                 ::  chic, psic
!
      IF (ZoL <= 0.0) then  !  Unstable conditions.
        chik = (1.0-15.0*ZoL)**0.25
        psik = 2.0*LOG(0.5*(1.0+chik))+LOG(0.5*(1.0+chik**2))
     &                                  -2.0*ATAN(chik)+pis2
        chic = (1.0-10.15*ZoL)**r3
        psic  = 1.5*LOG(r3*(chic**2+chic+1.0))
     &          - sqr3*ATAN((2.0*chic+1.0)/sqr3)+2.*pis2osqr3
        bulk_psiu_coare=psic+(psik-psic)/(1.0+ZoL**2)
      ELSE                    !  Stable conditions
        chic=-MIN(50.0,0.35*ZoL)
        bulk_psiu_coare=-((1.0+ZoL)+0.6667*(ZoL-14.28)*EXP(chic)+8.525)
      ENDIF
      return
      END FUNCTION bulk_psiu_coare
!=======================================================================
!    Coare stability function for tracers
!=======================================================================
      FUNCTION bulk_psit_coare (ZoL)
      IMPLICIT NONE
#include "params_bulk.h"
      REAL                 ::  bulk_psit_coare
      REAL                 ::  ZoL
!  Local variables declarations.
      REAL                 ::  chik, psik
      REAL                 ::  chic, psic
!
      IF (ZoL < 0.0) THEN
        chik  = (1.0-15.0*ZoL)**0.25
        psik  = 2.0*LOG(0.5*(1.0+chik**2))
        chic  = (1.0-34.15*ZoL)**r3     !for very unstable conditions
        psic  = 1.5*LOG((chic**2+chic+1.0)*r3)
     &               -sqr3*ATAN((2.0*chic+1.0)/sqr3)
     &               +2.*pis2osqr3
        bulk_psit_coare = psic+(psik-psic)/(1.0+ZoL**2) !match Kansas & free-conv. forms
      ELSE
        chic=-MIN(50.0,0.35*ZoL)
        bulk_psit_coare = -((1.0+2.0*ZoL/3.0)**1.5+
     &            0.6667*(ZoL-14.28)*EXP(chic)+8.525)
      ENDIF
      return
      END FUNCTION bulk_psit_coare
!=======================================================================
# endif

#else
      subroutine bulk_flux_empty
      return
      end
#endif
