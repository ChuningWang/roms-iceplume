!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
!$AGRIF_DO_NOT_TREAT
      INTEGER :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT 

#if defined OA_COUPLING || defined OW_COUPLING 
      INTEGER :: comp_id                       ! component identification	
      CHARACTER(len=6)   :: comp_name = 'crocox'

      INTEGER :: comp_ierror
      INTEGER :: oasis_part_id      
      INTEGER :: oasis_var_nodims(2) 
      INTEGER :: oasis_var_shape(4) 
      INTEGER :: oasis_var_type
      INTEGER , dimension(5) :: oasis_ig_paral ! Box partiton
 
      INTEGER, PARAMETER ::   nmaxfld = 60 ! Maximum number of coupling fields
      INTEGER, PARAMETER ::   nmaxatm =  5 ! Maximum number of atmospheric models

      CHARACTER(len = 64), DIMENSION(nmaxfld) :: srcv_clname, ssnd_clname   ! Coupling fields 
      INTEGER, DIMENSION(0:nmaxatm, nmaxfld) :: srcv_nid   , ssnd_nid
      common /exchange_fields_oasis3/ srcv_clname,ssnd_clname,
     &                                srcv_nid   ,ssnd_nid

      INTEGER :: oasis_time, oasis_runtime
      common /exchange_times_oasis3/ oasis_time, oasis_runtime

      REAL cplmsk(GLOBAL_2D_ARRAY,0:nmaxatm)
      common /coupling_mask/cplmsk
#endif /* OA_COUPLING */

