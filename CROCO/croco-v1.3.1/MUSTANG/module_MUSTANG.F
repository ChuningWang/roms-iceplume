#include "cppdefs.h"

!!! give access to dimensions, variables, parameters from hydro model to MUSTANG

#if defined MUSTANG

      module module_MUSTANG

   ! pass from fix (F) to free format (F90)

      use module_substance

#include "ocean3d.h"

      end module module_MUSTANG
#else
      module module_MUSTANG_empty
      end module module_MUSTANG_empty
      
#endif /* MUSTANG */
