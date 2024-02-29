      REAL, PARAMETER ::  CtoK      =  273.16   ! conversion factor for [C] to [K]
      REAL, PARAMETER ::  blk_Rgas  =  287.0596736665907  ! gas constant for dry air  [J/(kg K)]
      REAL, PARAMETER ::  blk_Rvap  =  461.5249933083879  ! gas constant for water vapor [J/(kg K)]
      REAL, PARAMETER ::  blk_Cpa   = 1004.708857833067
      REAL, PARAMETER ::  ip00      =    1.E-5
      REAL, PARAMETER ::  p00       =    1.E+5
      REAL, PARAMETER ::  rdocpd    = blk_Rgas/blk_Cpa
      REAL, PARAMETER ::  cpdord    = blk_Cpa/blk_Rgas
      REAL, PARAMETER ::  r_gas     = 8.314510
      REAL, PARAMETER ::  mm_dryair = 28.9644E-3
      REAL, PARAMETER ::  mm_water  = 18.0153E-3
      REAL, PARAMETER ::  cpvir     = blk_Rvap/blk_Rgas - 1.
      REAL, PARAMETER ::  MvoMa     = mm_water/mm_dryair
      REAL, PARAMETER ::  eps       = 1.d-8
      REAL, PARAMETER ::  r3        = 1./3.
      REAL, PARAMETER ::  pis2      = 2.*ATAN(1.)
      REAL, PARAMETER ::  sqr3      = SQRT(3.)
      REAL, PARAMETER ::  pis2osqr3 = pis2/sqr3
      REAL, PARAMETER ::  LMOmin    = -200.
      REAL, PARAMETER ::  LMOmax    = 0.25
      REAL, PARAMETER ::  blk_beta  = 1.2
      REAL, PARAMETER ::  blk_Zabl  = 600.
      REAL, PARAMETER ::  dWstar0   = 1.E-6
      REAL, PARAMETER ::  dTstar0   = 1.E-6
      REAL, PARAMETER ::  dQstar0   = 1.E-9
#ifdef BULK_LW
      REAL, PARAMETER ::  emiss_lw  =  0.985
      REAL, PARAMETER ::  SigmaSB   =  5.6697E-8
#endif
      REAL            :: psurf
#ifndef READ_PATM
      PARAMETER (psurf=100000.0)
#endif
#ifdef BULK_WASP
      REAL            ::  Awasp(0:3,4)
      REAL            ::  Bwasp(0:3,4)
      PARAMETER( Awasp = reshape((/  0.7    , 0.      , 0.      , 0.        ,
     &                              -9.202  , 2.265   ,-0.134   , 2.35e-3   ,
     &                               2.27   ,-6.67e-2 , 0.      , 0.        ,
     &                               0.0981 ,-4.13e-3 , 4.34e-5 , 1.16e-8  /)
     &                                                      , shape(Awasp)) )
      PARAMETER( Bwasp = reshape((/ -2.52   , 0.      , 0.      , 0.      ,
     &                              -0.4124 ,-0.2225  , 0.01178 ,-1.616e-4,
     &                              -2.41   , 4.30e-2 , 0.      , 0.      ,
     &                               0.     , 0.      , 0.      , 0.      /)
     &                                                      , shape(Bwasp)) )
      REAL, PARAMETER :: CWage    = 9.80665 / (16.*ATAN(1.))  ! g / (4 Pi)
      REAL, PARAMETER :: Charn0   = 0.018
      REAL, PARAMETER :: Charn1   = 0.1
      REAL, PARAMETER :: Charn2   = 0.002
#endif
#ifdef BULK_ECUMEV0
      real utu1,utu2,utt,utq1,utq2
      parameter(utu1=16.8,utu2=50.0)
      parameter(utt =33.0)
      parameter(utq1=29.0,utq2=33.0)
      real coefu10,coefu11,coefu12,coefu13
      parameter (coefu10= 1.3013E-03);parameter (coefu11=-1.2719E-04)
      parameter (coefu12= 1.3067E-05);parameter (coefu13=-2.2261E-07)
      real coefu20,coefu21,coefu22,coefu23,coefu24,Cdn0
      parameter (coefu20= 1.3633E-03);parameter (coefu21=-1.3056E-04)
      parameter (coefu22= 1.6212E-05);parameter (coefu23=-4.8208E-07)
      parameter (coefu24= 4.2684E-09);parameter (Cdn0   = 1.7828E-03)
      real coeft0,coeft1,coeft2,coeft3,coeft4,coeft5
      parameter (coeft0= 1.2536E-03);      parameter (coeft3=-4.3701E-07)
      parameter (coeft1=-1.2455E-04);      parameter (coeft4= 3.4517E-09)
      parameter (coeft2= 1.6038E-05);      parameter (coeft5= 3.5763E-12)
      real Chn0
      parameter (Chn0  = 3.1374E-03)
      real coefq10,coefq11,coefq12,coefq13,coefq14
      parameter (coefq10= 1.2687E-03);parameter (coefq11=-1.1384E-04)
      parameter (coefq12= 1.1467E-05);parameter (coefq13=-3.9144E-07)
      parameter (coefq14= 5.0864E-09)
      real coefq20,coefq21,coefq22
      parameter (coefq20= -1.3526E-03);parameter (coefq21=1.8229E-04)
      parameter (coefq22= -2.6995E-06)
      real Cen0
      parameter (Cen0  = 1.7232E-03)
#endif
#ifdef BULK_ECUMEV6
      real coefu0,coefu1,coefu2,coefu3,coefu4,coefu5
      parameter (coefu0= 1.00E-03);      parameter (coefu3= 2.32E-04)
      parameter (coefu1= 3.66E-02);      parameter (coefu4=-7.02E-06)
      parameter (coefu2=-1.92E-03);      parameter (coefu5= 6.40E-08)
      real coeft0,coeft1,coeft2,coeft3,coeft4,coeft5
      parameter (coeft0= 5.36E-03);      parameter (coeft3= 4.50E-04)
      parameter (coeft1= 2.90E-02);      parameter (coeft4=-2.06E-05)
      parameter (coeft2=-1.24E-03);      parameter (coeft5= 0.0     )
      real coefq0,coefq1,coefq2,coefq3,coefq4,coefq5
      parameter (coefq0= 1.00E-03);      parameter (coefq3= 0.0)
      parameter (coefq1= 3.59E-02);      parameter (coefq4= 0.0)
      parameter (coefq2=-2.87E-04);      parameter (coefq5= 0.0)
      real utu,utt,utq
      parameter(utu=40.0)
      parameter(utt=14.4)
      parameter(utq=10.0)
      real cdiru,cdirt,cdirq
      parameter(cdiru = coefu1+2.*coefu2*utu+3.*coefu3*utu**2
     &                  + 4.*coefu4*utu**3 + 5.*coefu5*utu**4 )
      parameter(cdirt = coeft1+2.*coeft2*utt+3.*coeft3*utt**2
     &                  + 4.*coeft4*utt**3                    )
      parameter(cdirq = coefq1 + 2.0*coefq2*utq               )
      !
      real ordou,ordot,ordoq
      parameter(ordou = coefu0 + coefu1*utu + coefu2*utu**2
     &                 + coefu3*utu**3
     &                 + coefu4*utu**4 + coefu5*utu**5      )
      parameter(ordot = coeft0 + coeft1*utt + coeft2*utt**2
     &                 + coeft3*utt**3
     &                 + coeft4*utt**4                      )
      parameter(ordoq = coefq0 + coefq1*utq + coefq2*utq**2 )
#endif
#ifdef SFLUX_CFB
# ifdef CFB_WIND_TRA
      ! wind correction: Ua-(1-sw)*Uo
      ! this is only used to correct heat flux (bulk_flux)
      real swparam
      parameter (swparam=0.3)
# endif
# ifdef CFB_STRESS
      ! wind-stress correction using wind speed:  rho0*sustr + s_tau*Uo
      !   s_tau = cfb_slope * wspd + cfb_offset [N.m^-3.s]
      !  (recommendended and default if BULK_FLUX - needs wspd data)
      real cfb_slope, cfb_offset
      parameter (cfb_slope=-0.0029)
      parameter (cfb_offset=0.008)
# endif
#endif
