#include "cppdefs.h"

#if defined MUSTANG

      module plug_MUSTANG_CROCO

      USE module_MUSTANG
      USE initMUSTANG, ONLY : MUSTANG_init
      USE sed_MUSTANG, ONLY : MUSTANG_update
      USE sed_MUSTANG, ONLY : MUSTANG_deposition
# ifdef MORPHODYN
      USE sed_MUSTANG, ONLY : MUSTANG_morpho
# endif

# include "coupler_define_MUSTANG.h"

      IMPLICIT NONE

      PRIVATE

      PUBLIC   mustang_update_main
      PUBLIC   mustang_deposition_main
      PUBLIC   mustang_init_main
# ifdef MORPHODYN
      PUBLIC   mustang_morpho_main
# endif

CONTAINS
!
!-----------------------------------------------------------------------
!
      subroutine mustang_update_main (tile)

      INTEGER :: tile
# include "ocean2d.h"
# include "compute_tile_bounds.h"
      CALL MUSTANG_update (Istr, Iend, Jstr, Jend,  & 
                   WATER_CONCENTRATION, Z0HYDRO,    &
                   WATER_ELEVATION,                 &
# if defined key_MUSTANG_lateralerosion || defined key_MUSTANG_bedload
                   BAROTROP_VELOCITY_U,             &
                   BAROTROP_VELOCITY_V,             &
# endif
                   SALREF_LIN, TEMPREF_LIN,         &
                   TRANSPORT_TIME_STEP)
      end subroutine
!
!-----------------------------------------------------------------------
!
      subroutine mustang_deposition_main (tile)

      integer :: tile
# include "ocean2d.h"
# include "compute_tile_bounds.h"
      CALL MUSTANG_deposition (Istr, Iend, Jstr, Jend, &
                   WATER_ELEVATION,                    &
                   WATER_CONCENTRATION)
      end subroutine
!
!-----------------------------------------------------------------------
!
      subroutine mustang_init_main (tile)

      integer :: tile
# include "ocean2d.h"
# include "compute_tile_bounds.h"

      CALL MUSTANG_init (Istr, Iend, Jstr, Jend, &
                    WATER_ELEVATION,                       &
# if defined MORPHODYN
                    DHSED,                                 &
# endif
                    RESIDUAL_THICKNESS_WAT, Z0HYDRO,       &
                    WATER_CONCENTRATION)
      end subroutine
!
!-----------------------------------------------------------------------
# ifdef MORPHODYN
      subroutine mustang_morpho_main (tile)

      integer :: tile
# include "ocean2d.h"
# include "compute_tile_bounds.h"

      CALL MUSTANG_morpho (Istr, Iend, Jstr, Jend, DHSED)

      end subroutine
# endif
!-----------------------------------------------------------------------
!
      end module plug_MUSTANG_CROCO

#else

      module plug_MUSTANG_CROCO_empty
      end module plug_MUSTANG_CROCO_empty
      
#endif /* MUSTANG */
