#include "cppdefs.h"

MODULE sedarr
   !!======================================================================
   !!                       ***  MODULE sedarr   ***
   !!              transform 1D (2D) array to a 2D (1D) table
   !!======================================================================
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   arr_2d_1d  : 2-D to 1-D
   !!   arr_1d_2d  : 1-D to 2-D
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_sed
   USE sed

   IMPLICIT NONE
   PRIVATE

   INTERFACE pack_arr
      MODULE PROCEDURE pack_arr_2d_1d , pack_arr_3d_2d 
   END INTERFACE

   INTERFACE unpack_arr
      MODULE PROCEDURE unpack_arr_1d_2d , unpack_arr_2d_3d 
   END INTERFACE

   !! * Routine accessibility
   PUBLIC pack_arr
   PUBLIC unpack_arr

   !!* Substitution
#  include "ocean2pisces.h90"

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: sedarr.F90 10222 2018-10-25 09:42:23Z aumont $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE pack_arr_2d_1d ( ndim1d, tab1d, tab2d, tab_ind )
#ifdef AGRIF
      USE ocean2pisces
#endif
      INTEGER, INTENT(in) ::  ndim1d
      REAL(wp), DIMENSION (PRIV_2D_BIOARRAY), INTENT(in) :: tab2d
      INTEGER, DIMENSION (ndim1d), INTENT (in) :: tab_ind
      REAL(wp), DIMENSION(ndim1d), INTENT (out) ::  tab1d

      INTEGER ::  jn, jid, jjd

      DO jn = 1, ndim1d
         jid        = MOD( tab_ind(jn) - 1, jpi ) + 1
         jjd        = ( tab_ind(jn) - 1 ) / jpi + 1
         tab1d(jn)  = tab2d(jid, jjd)
      END DO 

   END SUBROUTINE pack_arr_2d_1d

   SUBROUTINE unpack_arr_1d_2d ( ndim1d, tab2d, tab_ind, tab1d )

      INTEGER, INTENT ( in) ::  ndim1d
      INTEGER, DIMENSION (ndim1d) , INTENT (in) ::   tab_ind
      REAL(wp), DIMENSION(ndim1d), INTENT (in) ::   tab1d  
      REAL(wp), DIMENSION (PRIV_2D_BIOARRAY), INTENT ( out) ::  tab2d
      INTEGER ::  jn, jid, jjd

      DO jn = 1, ndim1d
         jid             = MOD( tab_ind(jn) - 1, jpi) + 1
         jjd             =    ( tab_ind(jn) - 1 ) / jpi  + 1
         tab2d(jid, jjd) = tab1d(jn)
      END DO

   END SUBROUTINE unpack_arr_1d_2d

   SUBROUTINE pack_arr_3d_2d ( ndim1d, tab2d, tab3d, tab_ind )

      INTEGER, INTENT(in) :: ndim1d      
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY,jpksed), INTENT(in) ::   tab3d      
      INTEGER, DIMENSION(ndim1d), INTENT (in) ::   tab_ind      
      REAL(wp), DIMENSION(ndim1d,jpksed), INTENT (out) ::    tab2d 
      INTEGER, DIMENSION(ndim1d) ::  jid, jjd        
      INTEGER ::    jk, jn , ji, jj

      DO jn = 1, ndim1d
         jid(jn) = MOD( tab_ind(jn) - 1, jpi ) + 1
         jjd(jn) = ( tab_ind(jn) - 1 ) / jpi + 1
      END DO
 
      DO jk = 1, jpksed
         DO jn = 1, ndim1d
            ji = jid(jn)
            jj = jjd(jn)
            tab2d(jn,jk)  = tab3d(ji,jj,jk) 
         ENDDO
      ENDDO

   END SUBROUTINE pack_arr_3d_2d


   SUBROUTINE unpack_arr_2d_3d ( ndim1d, tab3d, tab_ind, tab2d )

      INTEGER, INTENT(in) :: ndim1d      
      REAL(wp), DIMENSION(ndim1d,jpksed), INTENT(in) ::   tab2d      
      INTEGER, DIMENSION(ndim1d), INTENT (in) ::   tab_ind      
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY,jpksed), INTENT (out) ::    tab3d 
      INTEGER, DIMENSION(ndim1d) ::  jid, jjd        
      INTEGER ::   jk, jn , ji, jj
      !
      DO jn = 1, ndim1d
         jid(jn) = MOD( tab_ind(jn) - 1, jpi ) + 1
         jjd(jn) = ( tab_ind(jn) - 1 ) / jpi + 1
      END DO
 
      DO jk = 1, jpksed
         DO jn = 1, ndim1d
            ji = jid(jn)
            jj = jjd(jn)
            tab3d(ji, jj, jk) = tab2d(jn,jk)
         ENDDO
      ENDDO

   END SUBROUTINE unpack_arr_2d_3d

#endif

END MODULE sedarr
