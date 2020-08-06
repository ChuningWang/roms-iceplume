/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2020 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for ROMS-ICEPLUME Test.
**
** Application flag:   ICEPLUME_TEST
** Input script:       roms_iceplume_test.in
*/

/* general */
#define MASKING
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS

/* iceplume */
#define ICEPLUME
#ifdef ICEPLUME
# define ICEPLUME_TRACER
# define ICEPLUME_MELT
# define ICEPLUME_DET_AVERAGE
# define ICEPLUME_MELT_TRACER
# undef ICEPLUME_MIX
#endif

/* advection, dissipation, pressure grad, etc. */
#define UV_ADV
#define UV_COR

#define UV_VIS2
#define MIX_S_UV
#define TS_DIF2
#define MIX_GEO_TS

#ifdef SOLVE3D
# define DJ_GRADPS
#endif

/* vertical mixing */
#ifdef SOLVE3D
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
#endif

#ifdef SOLVE3D
# define GLS_MIXING

# define LIMIT_VDIFF
# define LIMIT_VVISC

# if defined GLS_MIXING || defined MY25_MIXING
#  define N2S2_HORAVG
#  define CRAIG_BANNER
#  define KANTHA_CLAYSON
#  define CHARNOK
# endif
#endif

#define ANA_DRAG
#define UV_DRAG_GRID
#define UV_LDRAG

/* tracers */
#define T_PASSIVE
#define ANA_PASSIVE
#define ANA_PSOURCE

/* analytical functionals */
#define ANA_GRID
#define ANA_MASK
#define ANA_INITIAL
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_M3OBC
#define ANA_TOBC
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_SPFLUX
#define ANA_BPFLUX
