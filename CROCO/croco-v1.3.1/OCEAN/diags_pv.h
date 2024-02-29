! This is include file "diags_ek.h"
!  ==== == ======= ==== ==========
!

#ifdef DIAGNOSTICS_PV
# if defined DIAGNOSTICS_PV_FULL
      real pv(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE pv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real pvd(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE pvd(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
      real Mrhs(GLOBAL_2D_ARRAY,N,2)
!CSDISTRIBUTE_RESHAPE Mrhs(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real Trhs(GLOBAL_2D_ARRAY,N,NTA)
!CSDISTRIBUTE_RESHAPE Trhs(BLOCK_PATTERN,*) BLOCK_CLAUSE


# ifdef AVERAGES
      real timediags_pv_avg
# if defined DIAGNOSTICS_PV_FULL
      real pv_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE pv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real pvd_avg(GLOBAL_2D_ARRAY,0:N)
!CSDISTRIBUTE_RESHAPE pvd_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
      real Mrhs_avg(GLOBAL_2D_ARRAY,N,2)
!CSDISTRIBUTE_RESHAPE Mrhs_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real Trhs_avg(GLOBAL_2D_ARRAY,N,NTA)
!CSDISTRIBUTE_RESHAPE Trhs_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE

# endif

# if defined DIAGNOSTICS_PV_FULL
      common /diag_pv/pv
     &       /diag_pvd/pvd
# endif
      common /diag_Mrhs/Mrhs
     &       /diag_Trhs/Trhs

# ifdef AVERAGES
      common /diag_timediags_pv_avg/timediags_pv_avg
# if defined DIAGNOSTICS_PV_FULL
      common /diag_pv_avg/pv_avg
     &       /diag_pvd_avg/pvd_avg
# endif
      common /diag_Mrhs_avg/Mrhs_avg
     &       /diag_Trhs_avg/Trhs_avg
# endif      


# if defined DIAGNOSTICS_PV && ! defined DIAGNOSTICS_UV && ! defined DIAGNOSTICS_EK_FULL
      real MXadv(GLOBAL_2D_ARRAY,N,2)
      real MYadv(GLOBAL_2D_ARRAY,N,2)
      real MHdiff(GLOBAL_2D_ARRAY,N,2)
      real MHmix(GLOBAL_2D_ARRAY,N,2,2)
      real MVmix(GLOBAL_2D_ARRAY,N,2)

      common /diag_MXadv/MXadv
     &       /diag_MYadv/MYadv
     &       /diag_MHdiff/MHdiff
     &       /diag_MHmix/MHmix
     &       /diag_MVmix/MVmix
     
# endif

# if defined DIAGNOSTICS_PV && ! defined DIAGNOSTICS_UV
      real MVmix2(GLOBAL_2D_ARRAY,N,2)
      real Mrate(GLOBAL_2D_ARRAY,N,2)
      common /diag_MVmix2/MVmix2
     &       /diag_Mrate/Mrate
# endif

# if defined DIAGNOSTICS_PV && ! defined DIAGNOSTICS_TS

      real THmix(GLOBAL_2D_ARRAY,N,NT)
      real TVmix(GLOBAL_2D_ARRAY,N,NT)
      real TForc(GLOBAL_2D_ARRAY,N,NT)
      real Trate(GLOBAL_2D_ARRAY,N,NT)

      common /diag_TForc/TForc
     &       /diag_THmix/THmix
     &       /diag_TVmix/TVmix
     &       /diag_Trate/Trate

# endif

# endif /* DIAGNOSTICS_PV */



 
