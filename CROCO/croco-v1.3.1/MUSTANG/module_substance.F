#include "cppdefs.h"

!!! give access to dimensions, variables, parameters from hydro model to substance

#if defined SUBSTANCE

      module module_substance

   ! pass from fix (F) to free format (F90)

#include "param.h"
#include "ncscrum.h"
#include "scalars.h"
#include "forces.h"
#include "grid.h"

      end module module_substance
#else
      module module_substance_empty
      end module module_substance_empty
      
#endif /* SUBSTANCE */
