/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2017 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Estuary with Sediment Transport Test.
**
** Application flag:   ESTUARY_TEST
** Input script:       ocean_estuary_test.in
**                     sediment_estuary_test.in
*/

/* general */
#define MASKING
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS

/* iceplume */
#define ICEPLUME
#ifdef ICEPLUME
# define ICEPLUME_VIRTUAL_MIX
# define ICEPLUME_DETRAIN_AVERAGE
# define ICEPLUME_TRACER
# define ICEPLUME_MELT
# undef ICEPLUME_MELT_TRACER
#endif

/* advection, dissipation, pressure grad, etc. */
#define UV_ADV
#define UV_COR

#define TS_MPDATA
#ifndef TS_MPDATA
# define TS_U3HADVECTION
# ifdef SOLVE3D
#  undef TS_C4VADVECTION
# endif
#endif

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
