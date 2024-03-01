module flocmod
!
! This module manage the flocculation distribution (in one cell "i,j,k").
!
! No use, no commons, no include this module is independant.
!
! Authors : R.Verney et all, 2011
! adaptation from MARS version to an independant module S.Le Gac, 2022
!

implicit none
private
! variables
public :: f_ws, f_diam, f_vol, f_rho, f_mass
! subroutines
public :: flocmod_alloc
public :: flocmod_init
public :: flocmod_main
public :: flocmod_comp_g

! Declaration

! parameters
real, parameter    :: pi = 3.14159265358979323846 ! pi value in the module
real, parameter    :: gravity = 9.81              ! gravity value in the module
integer, parameter :: rlg = 8                     ! double precision

! variables of the module
integer :: nv_mud ! number of mud class
real :: rhoref ! water density
logical :: l_0Dcase ! true if it is a run for 0Dcase comparison (need to stop settling process)
integer :: ierrorlog ! to write in the appropriate file unit

!TODO add units !!
real, dimension(:), allocatable :: f_ws   ! flocs settling velocities
real, dimension(:), allocatable :: f_diam ! floc diameter
real, dimension(:), allocatable :: f_vol  ! floc volume
real, dimension(:), allocatable :: f_rho  ! floc density
real, dimension(:), allocatable :: f_mass ! floc mass

real, dimension(:,:), allocatable :: f_coll_prob_sh ! shear agregation collision probability
real, dimension(:,:), allocatable :: f_coll_prob_ds ! differential settling collision probability

real, dimension(:,:,:), allocatable :: f_g1_sh ! shear agregation gain term
real, dimension(:,:,:), allocatable :: f_g1_ds ! differential settling agregation gain term
real, dimension(:,:), allocatable :: f_l1_sh ! shear agregation loss term
real, dimension(:,:), allocatable :: f_l1_ds ! differential settling agregation loss term  
real, dimension(:,:), allocatable :: f_g3 ! fragmentation gain term
real, dimension(:), allocatable :: f_l3 ! fragmentation loss term

logical :: l_ASH ! set to .true. if aggregation by shear
logical :: l_ADS ! set to .true. if aggregation by differential settling
logical :: l_COLLFRAG ! set to .true. if fragmentation by collision
integer :: f_ero_iv ! fragment class (mud variable index corresponding to 
    ! the eroded particle size - typically 1)
real :: f_ater ! ternary fragmentation factor : proportion of 
    ! flocs fragmented as half the size of the initial binary fragments 
    ! (0.0 if full binary fragmentation, 0.5 if ternary fragmentation)
real :: f_dmin_frag ! minimum diameter for fragmentation 
    ! (default 10e-6 microns)
real :: f_ero_frac ! floc erosion (% of floc mass eroded) 
    ! (default 0.05)
real :: f_ero_nbfrag ! nb of fragments produced by erosion 
    ! (default 2.0)
real :: f_mneg_param ! negative mass after 
    ! flocculation/fragmentation allowed before redistribution 
    ! (default 0.001 g/l)
real :: f_collfragparam ! fraction of shear aggregation leading 
    ! to fragmentation by collision (default 0.0, must be less than 1.0)
real :: f_cfcst ! fraction of mass lost by flocs if fragmentation
    ! by collision .. (default : =3./16.)
real :: f_fp ! relative depth of inter particle penetration  
    ! (default =0.1) (McAnally, 1999)
real :: f_fy ! floc yield strength  (default= 1.0e-10) 
    ! (Winterwerp, 2002)
real :: f_dp0 ! primary particle size (default 4.e-6 m)
real :: f_alpha ! flocculation efficiency parameter 
    ! (default 0.15)
real :: f_beta ! floc break up parameter (default 0.1)
real :: f_nb_frag ! nb of fragments of equal size by shear 
    ! fragmentation (default 2.0 as binary fragmentation)
real :: f_nf ! fractal dimension (default 2.0, usual range from 
    ! 1.6 to 2.8)
real :: f_clim ! min concentration below which flocculation 
    !processes are not calculated


contains

!==============================================================================
subroutine flocmod_alloc(n)
! Initialize the dimensions of number of mud class and allocate the global 
! arrays of the module.
implicit none

integer, intent(in) :: n  ! number of mud classes

nv_mud = n
  
! floc characteristics
allocate(  f_ws(1:nv_mud))       ! flocs settling velocities 
allocate(f_diam(1:nv_mud))       ! floc diameter
allocate( f_vol(1:nv_mud))       ! floc volume
allocate( f_rho(1:nv_mud))       ! floc density
allocate(f_mass(0:nv_mud + 1))   ! floc mass

! agregation kernels
allocate(f_coll_prob_sh(1:nv_mud, 1:nv_mud)) ! shear agregation collision probability
allocate(f_coll_prob_ds(1:nv_mud, 1:nv_mud)) ! differential settling collision probability

allocate(f_g1_sh(1:nv_mud, 1:nv_mud, 1:nv_mud)) ! shear agregation gain term
allocate(f_g1_ds(1:nv_mud, 1:nv_mud, 1:nv_mud)) ! differential settling agregation gain term
allocate(f_l1_sh(1:nv_mud, 1:nv_mud)) ! shear agregation loss term
allocate(f_l1_ds(1:nv_mud, 1:nv_mud)) ! differential settling agregation loss term 
allocate(f_g3(1:nv_mud, 1:nv_mud)) ! fragmentation gain term
allocate(f_l3(1:nv_mud)) ! fragmentation loss term

f_ws(1:nv_mud) = 0.0
f_diam(1:nv_mud) = 0.
f_vol(1:nv_mud) = 0.
f_rho(1:nv_mud) = 0.
f_mass(0:nv_mud+1) = 0.

f_coll_prob_sh(1:nv_mud,1:nv_mud) = 0.
f_coll_prob_ds(1:nv_mud,1:nv_mud) = 0.

f_g1_sh(1:nv_mud,1:nv_mud,1:nv_mud) = 0.
f_g1_ds(1:nv_mud,1:nv_mud,1:nv_mud) = 0.
f_l1_sh(1:nv_mud,1:nv_mud) = 0.
f_l1_ds(1:nv_mud,1:nv_mud) = 0.
f_g3(1:nv_mud,1:nv_mud) = 0.
f_l3(1:nv_mud) = 0.

return
end subroutine ! flocmod_alloc

!==================================================================================================
subroutine flocmod_init(l_ADS_in, l_ASH_in, l_COLLFRAG_in,       &
    f_dp0_in, f_nf_in, f_nb_frag_in, f_alpha_in, f_beta_in, f_ater_in,    &
    f_ero_frac_in, f_ero_nbfrag_in, f_ero_iv_in, f_mneg_param_in,   &
    f_collfragparam_in, f_dmin_frag_in, f_cfcst_in, f_fp_in, f_fy_in,  &
    f_clim_in, f_diam_in, ros_in, rhoref_in, l_0Dcase_in, ierrorlog_in)
! Initalization of flocs characteristics from input values
implicit none

logical, intent(in)   :: l_ASH_in ! set to .true. if aggregation by shear
logical, intent(in)   :: l_ADS_in ! set to .true. if aggregation by differential settling
logical, intent(in)   :: l_COLLFRAG_in ! set to .true. if fragmentation by collision
integer, intent(in)   :: f_ero_iv_in ! fragment class (mud variable index corresponding to 
    ! the eroded particle size - typically 1)
real, intent(in)   :: f_ater_in ! ternary fragmentation factor : proportion of 
    ! flocs fragmented as half the size of the initial binary fragments 
    ! (0.0 if full binary fragmentation, 0.5 if ternary fragmentation)
real, intent(in)   :: f_dmin_frag_in ! minimum diameter for fragmentation 
    ! (default 10e-6 microns)
real, intent(in)   :: f_ero_frac_in ! floc erosion (% of floc mass eroded) 
    ! (default 0.05)
real, intent(in)   :: f_ero_nbfrag_in ! nb of fragments produced by erosion 
    ! (default 2.0)
real, intent(in)   :: f_mneg_param_in ! negative mass after 
    ! flocculation/fragmentation allowed before redistribution 
    ! (default 0.001 g/l)
real, intent(in)   :: f_collfragparam_in ! fraction of shear aggregation leading 
    ! to fragmentation by collision (default 0.0, must be less than 1.0)
real, intent(in)   :: f_cfcst_in ! fraction of mass lost by flocs if fragmentation
    ! by collision .. (default : =3./16.)
real, intent(in)   :: f_fp_in ! relative depth of inter particle penetration  
    ! (default =0.1) (McAnally, 1999)
real, intent(in)   :: f_fy_in ! floc yield strength  (default= 1.0e-10) 
    ! (Winterwerp, 2002)
real, intent(in)   :: f_dp0_in ! primary particle size (default 4.e-6 m)
real, intent(in)   :: f_alpha_in ! flocculation efficiency parameter 
    ! (default 0.15)
real, intent(in)   :: f_beta_in ! floc break up parameter (default 0.1)
real, intent(in)   :: f_nb_frag_in ! nb of fragments of equal size by shear 
    ! fragmentation (default 2.0 as binary fragmentation)
real, intent(in)   :: f_nf_in ! fractal dimension (default 2.0, usual range from 
    ! 1.6 to 2.8)
real, intent(in)   :: f_clim_in ! min concentration below which flocculation 
    !processes are not calculated
real, dimension(nv_mud), intent(in)   :: f_diam_in
real, dimension(nv_mud), intent(in)   :: ros_in
real, intent(in)   :: rhoref_in
logical, intent(in) :: l_0Dcase_in
integer, intent(in) :: ierrorlog_in

l_ADS = l_ADS_in
l_ASH = l_ASH_in
l_COLLFRAG = l_COLLFRAG_in 
f_dp0 = f_dp0_in
f_nf = f_nf_in
f_nb_frag = f_nb_frag_in
f_alpha = f_alpha_in
f_beta = f_beta_in
f_ater = f_ater_in
f_ero_frac = f_ero_frac_in
f_ero_nbfrag = f_ero_nbfrag_in
f_ero_iv = f_ero_iv_in
f_mneg_param = f_mneg_param_in
f_collfragparam = f_collfragparam_in
f_dmin_frag = f_dmin_frag_in
f_cfcst = f_cfcst_in
f_fp = f_fp_in
f_fp = f_fp_in
f_clim = f_clim_in
rhoref = rhoref_in
l_0Dcase = l_0Dcase_in
ierrorlog = ierrorlog_in

f_diam(:) = f_diam_in(:)
f_vol(:) = pi / 6. * (f_diam(:))**3
f_rho(:) = rhoref + (ros_in(:) - rhoref) * (f_dp0 / f_diam(:))**(3. - f_nf)
f_mass(1:nv_mud) = f_vol(1:nv_mud) * (f_rho(1:nv_mud) - rhoref) / (1 - rhoref / ros_in(1:nv_mud))
f_mass(nv_mud+1) = f_mass(nv_mud) * 2. + 1.  
if (f_diam(1) .eq. f_dp0) then
  f_mass(1) = f_vol(1) * ros_in(1)
endif
if (l_0Dcase) then
    ! for 0D test case, we need to suppress all settling process
    f_ws(1:nv_mud) = 0.
else
    f_ws(1:nv_mud) = gravity * (f_rho(1:nv_mud) - rhoref) * f_diam(1:nv_mud)**2. / (18. * 0.001)
endif

! kernels computation
call flocmod_kernels()    

return
end subroutine ! flocmod_init


!==================================================================================================
subroutine flocmod_kernels()
! Purpose : computations of agregation/fragmentation kernels for FLOCMOD
! Called by : flocmod_init
implicit none

real ::  f_weight, mult
integer  ::  iv1, iv2, iv3

! compute collision probability
call flocmod_agregation_statistics

!********************************************************************************
! agregation : GAIN : f_g1_sh and f_g1_ds
!********************************************************************************

do iv1=1,nv_mud-1
    do iv2=1,nv_mud
        do iv3=iv2,nv_mud
            if(        (f_mass(iv2)+f_mass(iv3)) .gt. f_mass(iv1-1) &
                .and. ((f_mass(iv2)+f_mass(iv3)) .LE. f_mass(iv1))) then

                f_weight=(f_mass(iv2)+f_mass(iv3)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))

            else if (    (f_mass(iv2)+f_mass(iv3)) .gt. f_mass(iv1) &
                .and.   ((f_mass(iv2)+f_mass(iv3)) .LT. f_mass(iv1+1))) then
            
                f_weight=1.-(f_mass(iv2)+f_mass(iv3)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1))

            else

                f_weight=0.0

            endif

            f_g1_sh(iv2,iv3,iv1)=f_weight*f_alpha*f_coll_prob_sh(iv2,iv3)*(f_mass(iv2)+f_mass(iv3))/f_mass(iv1)
            f_g1_ds(iv2,iv3,iv1)=f_weight*f_alpha*f_coll_prob_ds(iv2,iv3)*(f_mass(iv2)+f_mass(iv3))/f_mass(iv1)

        enddo
    enddo
enddo
! difference for coarsest mud class >> no transfer towards coarser class because it is the last one
iv1=nv_mud
do iv2=1,nv_mud
    do iv3=iv2,nv_mud
        if(        (f_mass(iv2)+f_mass(iv3)) .gt. f_mass(iv1-1) &
            .and. ((f_mass(iv2)+f_mass(iv3)) .LE. f_mass(iv1))) then

            f_weight=(f_mass(iv2)+f_mass(iv3)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))

        else if (  (f_mass(iv2)+f_mass(iv3)) .gt. f_mass(iv1) &
            .and. ((f_mass(iv2)+f_mass(iv3)) .LT. f_mass(iv1+1))) then

            f_weight=1.

        else
            f_weight=0.0
        endif

        f_g1_sh(iv2,iv3,iv1)=f_weight*f_alpha*f_coll_prob_sh(iv2,iv3)*(f_mass(iv2)+f_mass(iv3))/f_mass(iv1)
        f_g1_ds(iv2,iv3,iv1)=f_weight*f_alpha*f_coll_prob_ds(iv2,iv3)*(f_mass(iv2)+f_mass(iv3))/f_mass(iv1)

    enddo
enddo

!********************************************************************************
! Shear fragmentation : GAIN : f_g3
!********************************************************************************    
do iv1=1,nv_mud
    do iv2=iv1,nv_mud     
        if (f_diam(iv2) > f_dmin_frag) then
        ! binary fragmentation
            if (f_mass(iv2)/f_nb_frag .gt. f_mass(iv1-1) &
            .and. f_mass(iv2)/f_nb_frag .LE. f_mass(iv1)) then
                if (iv1 == 1) then 
                    f_weight=1.
                else
                    f_weight=(f_mass(iv2)/f_nb_frag-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))
                endif
            elseif (f_mass(iv2)/f_nb_frag .gt. f_mass(iv1) &
                    .and. f_mass(iv2)/f_nb_frag .LT. f_mass(iv1+1)) then

                f_weight=1.-(f_mass(iv2)/f_nb_frag-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1))
            else      
                f_weight=0.
            endif
        else
            f_weight=0.0
        endif
    
        f_g3(iv2,iv1)=f_g3(iv2,iv1)+(1.-f_ero_frac)*(1.-f_ater)*f_weight*f_beta &
                    *f_diam(iv2)*((f_diam(iv2)-f_dp0)/f_dp0)**(3.-f_nf)           &
                    *f_mass(iv2)/f_mass(iv1)

        ! ternary fragmentation
        if (f_diam(iv2) .gt. f_dmin_frag) then
            if (f_mass(iv2)/(2.*f_nb_frag) .gt. f_mass(iv1-1) &
            .and. f_mass(iv2)/(2.*f_nb_frag) .LE. f_mass(iv1)) then

                if (iv1 == 1) then 
                    f_weight=1.
                else
                    f_weight=(f_mass(iv2)/(2.*f_nb_frag)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))
                endif

            else if (f_mass(iv2)/(2.*f_nb_frag) .gt. f_mass(iv1) &
            .and. f_mass(iv2)/(2.*f_nb_frag) .LT. f_mass(iv1+1)) then

                f_weight=1.-(f_mass(iv2)/(2.*f_nb_frag)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1))

            else
                f_weight=0.
        
            endif
            ! update for ternary fragments
            f_g3(iv2,iv1)=f_g3(iv2,iv1)+(1.-f_ero_frac)*(f_ater)*f_weight*f_beta &
                    *f_diam(iv2)*((f_diam(iv2)-f_dp0)/f_dp0)**(3.-f_nf)           &
                    *f_mass(iv2)/f_mass(iv1)   

            ! Floc erosion

            if ((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .gt. f_mass(f_ero_iv)) then

                if (((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .gt. f_mass(iv1-1)) &
                .and. (f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .LE. f_mass(iv1)) then
            
                    if (iv1 == 1) then
                        f_weight=1.
                    else
                        f_weight=(f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))
                    endif
        
                else if ((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .gt. f_mass(iv1) &
                .and. (f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .LT. f_mass(iv1+1)) then
            
                    f_weight=1.-(f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1))
        
                else
                    f_weight=0.
                endif
    
                ! update for eroded floc masses 
    
                f_g3(iv2,iv1)=f_g3(iv2,iv1)+f_ero_frac*f_weight*f_beta                    &
                        *f_diam(iv2)*((f_diam(iv2)-f_dp0)/f_dp0)**(3.-f_nf)           &
                        *(f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag)/f_mass(iv1)
    
                if (iv1 == f_ero_iv) then
            
                    f_g3(iv2,iv1)=f_g3(iv2,iv1)+f_ero_frac*f_beta                           &
                        *f_diam(iv2)*((f_diam(iv2)-f_dp0)/f_dp0)**(3.-f_nf)           &
                        *f_ero_nbfrag*f_mass(f_ero_iv)/f_mass(iv1)
                endif
            endif
        endif ! condition on f_dmin_frag
    enddo
enddo  
 
!********************************************************************************
!  Shear agregation : LOSS : f_l1
!********************************************************************************
do iv1 = 1, nv_mud
    do iv2 = 1, nv_mud
        if(iv2 == iv1) then
            mult=2.
        else
            mult=1.
        endif 
        f_l1_sh(iv2, iv1) = mult * f_alpha * f_coll_prob_sh(iv2, iv1) 
        f_l1_ds(iv2, iv1) = mult * f_alpha * f_coll_prob_ds(iv2, iv1) 
    enddo
enddo
 
!********************************************************************************
!  Shear fragmentation : LOSS : f_l2
!********************************************************************************
do iv1 = 1, nv_mud
    if (f_diam(iv1) > f_dmin_frag) then
        ! shear fragmentation
        f_l3(iv1) = f_l3(iv1)+(1. - f_ero_frac)*f_beta*f_diam(iv1)*((f_diam(iv1)-f_dp0)/f_dp0)**(3. - f_nf)
        ! shear erosion
        if ((f_mass(iv1) - f_mass(f_ero_iv)*f_ero_nbfrag) > f_mass(f_ero_iv)) then
            f_l3(iv1)=f_l3(iv1)+f_ero_frac*f_beta*f_diam(iv1)*((f_diam(iv1)-f_dp0)/f_dp0)**(3. - f_nf)
        endif
    endif    
enddo
return
end subroutine flocmod_kernels

!==================================================================================================
subroutine flocmod_agregation_statistics
! Purpose : computation of shear / differential settling statistics
! Called by : flocmod_kernels
implicit none

integer :: iv1, iv2
    
do iv1 = 1, nv_mud
    do iv2 = 1, nv_mud
        f_coll_prob_sh(iv1, iv2) = 1. / 6. * (f_diam(iv1) + f_diam(iv2))**3.
        
        f_coll_prob_ds(iv1, iv2) = 0.25 * pi * (f_diam(iv1) + f_diam(iv2))**2. &
                                  * abs(f_ws(iv1) - f_ws(iv2))
    enddo
enddo
return
end subroutine flocmod_agregation_statistics

!==================================================================================================
subroutine flocmod_comp_g(Gval, time)
! Purpose : compute shear rate to estimate shear aggregation and erosion  
! reproducing flocculation experiment Verney et al., 2011
implicit none 
    
real,intent(out)    :: Gval
real(kind=rlg),intent(in)     :: time ! time in seconds

Gval=0.0
if (time .lt. 7201.0) then
    Gval=1.0
elseif (time .lt. 8401.0) then
    Gval=2.0
elseif (time .lt. 9601.0) then
    Gval=3.0  
elseif (time .lt. 10801.0) then
    Gval=4.0
elseif (time .lt. 12601.0) then
    Gval=12.0
elseif (time .lt. 13801.0) then
    Gval=4.0
elseif (time .lt. 15001.0) then
    Gval=3.0
elseif (time .lt. 16201.0) then
    Gval=2.0
elseif (time .lt. 21601.0) then
    Gval=1.0
elseif (time .lt. 25201.0) then
    Gval=0.0
elseif (time .lt. 30601.0) then
    Gval=1.0
elseif (time .lt. 31801.0) then
    Gval=2.0                     
elseif (time .lt. 33001.0) then
    Gval=3.0       
elseif (time .lt. 34201.0) then
    Gval=4.0
elseif (time .lt. 36001.0) then
    Gval=12.0
elseif (time .lt. 37201.0) then
    Gval=4.0
elseif (time .lt. 38401.0) then
    Gval=3.0
elseif (time .lt. 39601.0) then
    Gval=2.0
elseif (time .lt. 45001.0) then
    Gval=1.0
elseif (time .lt. 48601.0) then
    Gval=0.0
elseif (time .lt. 54001.0) then
    Gval=1.0
elseif (time .lt. 55201.0) then
    Gval=2.0                     
elseif (time .lt. 56401.0) then
    Gval=3.0       
elseif (time .lt. 57601.0) then
    Gval=4.0
else 
    Gval=12.0
endif

return
end subroutine flocmod_comp_g

!==================================================================================================
subroutine flocmod_main(dt, cwater, gradvit)
!Purpose : main routine used to compute aggregation /fragmentation processes in a cell
implicit none
                        
real(kind=rlg), intent(in) :: dt 
real, intent(in) :: gradvit
real,dimension(nv_mud), intent(inout) :: cwater         

! Local declarations
real                      :: cvtotmud, cvtotmudref, Gval, mneg, sum_flocs
real(kind=rlg)            :: dttemp, f_dt
real, dimension(1:nv_mud) :: cv_tmp, NNin, NNout

!Gval : shear rate value calculated or estimated from the hydrodynmic model
Gval = gradvit

f_dt = dt
dttemp = 0.0_rlg

cv_tmp(1:nv_mud) = cwater   ! concentration of all mud classes in one grid cell
cvtotmudref = sum(cv_tmp(1:nv_mud))
cv_tmp(1:nv_mud) = MAX(0.0, cv_tmp(1:nv_mud))
cvtotmud = sum(cv_tmp(1:nv_mud))

NNin(1:nv_mud) = cv_tmp(1:nv_mud) / f_mass(1:nv_mud)

if (cvtotmud .gt. f_clim) then

    do while (dttemp .LE. dt)

        call flocmod_comp_fsd(NNin,NNout,Gval,f_dt)
        call flocmod_mass_control(NNout,mneg)
        if (mneg .gt. f_mneg_param) then
            do while (mneg .gt. f_mneg_param)
                f_dt=MIN(f_dt/2._rlg,dt-dttemp)
                call flocmod_comp_fsd(NNin,NNout,Gval,f_dt)
                call flocmod_mass_control(NNout,mneg)  
            enddo
        endif
        dttemp=dttemp+f_dt
        NNin(:)=NNout(:) ! update new Floc size distribution

        ! redistribute negative masses if any on positive classes, depends on f_mneg_param
        call flocmod_mass_redistribute(NNin) 

        if (dttemp == dt) exit

    enddo ! loop on full dt

    sum_flocs = sum(NNin(1:nv_mud) * f_mass(1:nv_mud))
    NNin(:) = NNin(:) * cvtotmudref / sum_flocs

    if (abs(sum_flocs-cvtotmud).le. 0.01 * cvtotmud) then
        NNin(:)=NNin(:)*cvtotmudref/sum_flocs
    else
        write(ierrorlog,*) 'CAUTION flocculation routine not conservative!!!'
        write(ierrorlog,*) 'before : cvtotmud= ',cvtotmud
        write(ierrorlog,*) 'after  : cvtotmud= ',sum_flocs
        write(ierrorlog,*) 'absolute difference  : cvtotmud= ',abs(cvtotmud-sum_flocs)
        write(ierrorlog,*) 'Simultation stopped'

        STOP
    endif

endif ! only if cvtotmud > f_clim

! update mass concentration for all mud classes
cwater = NNin(1:nv_mud) * f_mass(1:nv_mud)

return
end subroutine flocmod_main

!==================================================================================================
subroutine flocmod_comp_fsd(NNin, NNout, Gval, f_dt)
! Purpose : computation of floc size distribution
implicit none

real, intent(in) :: Gval
real(kind=rlg), intent(in) :: f_dt
real, dimension(1:nv_mud), intent(in)  :: NNin
real, dimension(1:nv_mud), intent(out) :: NNout


integer  :: iv1, iv2, iv3
real :: mu, tmp_g1, tmp_g3, tmp_l1, tmp_l3, tmp_l4, tmp_g4
real, dimension(1:nv_mud, 1:nv_mud, 1:nv_mud) :: f_g1_tmp
real, dimension(1:nv_mud, 1:nv_mud) :: f_l1_tmp
real, dimension(1:nv_mud,1:nv_mud,1:nv_mud) :: f_g4 ! Collision fragmentation gain term
real, dimension(1:nv_mud,1:nv_mud) :: f_l4

tmp_g1 = 0.0
tmp_g3 = 0.0
tmp_g4 = 0.0
tmp_l1 = 0.0
tmp_l3 = 0.0
tmp_l4 = 0.0    
f_g1_tmp(1:nv_mud, 1:nv_mud, 1:nv_mud) = 0.0
f_l1_tmp(1:nv_mud, 1:nv_mud) = 0.0

if (l_COLLFRAG) call flocmod_collfrag(Gval, f_g4, f_l4)

do iv1=1,nv_mud
    do iv2=1,nv_mud
        do iv3=1,nv_mud
            if (l_ASH) then
                f_g1_tmp(iv2,iv3,iv1)=f_g1_tmp(iv2,iv3,iv1)+f_g1_sh(iv2,iv3,iv1)*Gval   
            endif
            if (l_ADS) then
                f_g1_tmp(iv2,iv3,iv1)=f_g1_tmp(iv2,iv3,iv1)+f_g1_ds(iv2,iv3,iv1)   
            endif

            tmp_g1=tmp_g1+(NNin(iv3)*(f_g1_tmp(iv2,iv3,iv1))*NNin(iv2))

            if (l_COLLFRAG) then
                tmp_g4=tmp_g4+(NNin(iv3)*(f_g4(iv2,iv3,iv1)*Gval)*NNin(iv2))
            endif
        enddo

        tmp_g3=tmp_g3+f_g3(iv2,iv1)*NNin(iv2)*Gval**1.5

        if (l_ASH) then
            f_l1_tmp(iv2,iv1)=f_l1_tmp(iv2,iv1)+f_l1_sh(iv2,iv1)*Gval
        endif
        if (l_ADS) then
            f_l1_tmp(iv2,iv1)=f_l1_tmp(iv2,iv1)+f_l1_ds(iv2,iv1)
        endif

        tmp_l1=tmp_l1+(f_l1_tmp(iv2,iv1))*NNin(iv2)

        if (l_COLLFRAG) then
            tmp_l4=tmp_l4+(f_l4(iv2,iv1)*Gval)*NNin(iv2)
        endif
    enddo

    tmp_l1=tmp_l1*NNin(iv1)
    tmp_l4=tmp_l4*NNin(iv1)

    tmp_l3=f_l3(iv1)*Gval**1.5*NNin(iv1)

    NNout(iv1)=NNin(iv1)+f_dt*(tmp_g1+tmp_g3+tmp_g4-(tmp_l1+tmp_l3+tmp_l4))

    tmp_g1=0.0
    tmp_g3=0.0
    tmp_g4=0.0
    tmp_l1=0.0
    tmp_l3=0.0
    tmp_l4=0.0    
enddo  

end subroutine flocmod_comp_fsd

!==================================================================================================
subroutine flocmod_collfrag(Gval, f_g4, f_l4) 
! Purpose : computation of collision fragmentation term, based on McAnally and Mehta, 2001
!
! Called by : flocmod_comp_fsd
implicit none

real, intent(in) :: Gval
real, dimension(1:nv_mud,1:nv_mud,1:nv_mud), intent(out) :: f_g4   ! Collision fragmentation gain term
real, dimension(1:nv_mud,1:nv_mud), intent(out)          :: f_l4   ! Collision fragmentation loss term  


integer :: iv1, iv2, iv3
real :: gcolfragmin, gcolfragmax, gcolfragiv1, gcolfragiv2, f_weight, mult

f_g4(1:nv_mud,1:nv_mud,1:nv_mud) = 0.0
f_l4(1:nv_mud,1:nv_mud) = 0.0

do iv1 = 1, nv_mud
do iv2 = 1, nv_mud
do iv3 = iv2, nv_mud
! fragmentation after collision probability based on Gval for particles iv2 and iv3
! gcolfrag=(collision induced shear) / (floc strength)

gcolfragmin=2.*(Gval*(f_diam(iv2)+f_diam(iv3)))**2.*f_mass(iv2)*f_mass(iv3)  &
  /(pi*f_fy*f_fp*f_diam(iv3)**2.*(f_mass(iv2)+f_mass(iv3))         &
 *((f_rho(iv3)-rhoref)/rhoref)**(2./(3.-f_nf)))

gcolfragmax=2.*(Gval*(f_diam(iv2)+f_diam(iv3)))**2.*f_mass(iv2)*f_mass(iv3)  &
 /(pi*f_fy*f_fp*f_diam(iv2)**2.*(f_mass(iv2)+f_mass(iv3))         &
*((f_rho(iv2)-rhoref)/rhoref)**(2./(3.-f_nf)))


! first case : iv3 not eroded, iv2 eroded forming 2 particles : iv3+f_cfcst*iv2 / iv2-f_cfcst*iv2
if (gcolfragmin.lt.1. .and. gcolfragmax.ge.1) then  

    if (((f_mass(iv3)+f_cfcst*f_mass(iv2)) .gt. f_mass(iv1-1)) .and.  &
    ((f_mass(iv3)+f_cfcst*f_mass(iv2)).le.f_mass(iv1))) then
        f_weight=((f_mass(iv3)+f_cfcst*f_mass(iv2)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1)))
    elseif (f_mass(iv3)+f_cfcst*f_mass(iv2).gt.f_mass(iv1)  .and. &
    f_mass(iv3)+f_cfcst*f_mass(iv2).lt.f_mass(iv1+1)) then
        if (iv1 == nv_mud) then 
            f_weight=1.
        else
            f_weight=1.-((f_mass(iv3)+f_cfcst*f_mass(iv2)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)))
        endif
    else
        f_weight=0.0
    endif 

    f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   &
        *(f_mass(iv3)+f_cfcst*f_mass(iv2))/f_mass(iv1)

    if (f_mass(iv2)-f_cfcst*f_mass(iv2).gt.f_mass(iv1-1)   .and. &
    f_mass(iv2)-f_cfcst*f_mass(iv2).le.f_mass(iv1)) then
        f_weight=((f_mass(iv2)-f_cfcst*f_mass(iv2)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1)))
    elseif (f_mass(iv2)-f_cfcst*f_mass(iv2).gt.f_mass(iv1)  .and.  &
    f_mass(iv2)-f_cfcst*f_mass(iv2).lt.f_mass(iv1+1)) then
        if (iv1.eq.nv_mud) then 
            f_weight=1.
        else
            f_weight=1.-((f_mass(iv2)-f_cfcst*f_mass(iv2)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)))
        endif
    else
        f_weight=0.0
    endif 

    f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   &
       *(f_mass(iv2)-f_cfcst*f_mass(iv2))/f_mass(iv1)


! second case : iv3 eroded and iv2 eroded forming 3 particles : iv3-f_cfcst*iv3 / iv2-f_cfcst*iv2 / f_cfcst*iv3+f_cfcst*iv2
elseif (gcolfragmin.ge.1. .and. gcolfragmax.ge.1) then  ! iv2 and iv3 eroded forming new (third) particle

    if (f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3).gt.f_mass(iv1-1) .and.  &
    f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3).le.f_mass(iv1)) then
        f_weight=((f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1)))
    elseif (f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3).gt.f_mass(iv1) .and.  &
    f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3).lt.f_mass(iv1+1)) then
        if (iv1.eq.nv_mud) then 
        f_weight=1.
        else
        f_weight=1.-((f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)))
        endif
    else
        f_weight=0.0
    endif 

    f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   &
      *(f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3))/f_mass(iv1)

    if ((1.-f_cfcst)*f_mass(iv2).gt.f_mass(iv1-1) .and.  &
    (1.-f_cfcst)*f_mass(iv2).le.f_mass(iv1)) then
        f_weight=((1.-f_cfcst)*f_mass(iv2)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))
    elseif ((1.-f_cfcst)*f_mass(iv2).gt.f_mass(iv1) .and.  &
    (1.-f_cfcst)*f_mass(iv2).lt.f_mass(iv1+1)) then
        if (iv1.eq.nv_mud) then 
        f_weight=1.
        else
        f_weight=1.-(((1.-f_cfcst)*f_mass(iv2)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)))
        endif
    else
        f_weight=0.0
    endif 

    f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   &
      *((1.-f_cfcst)*f_mass(iv2))/f_mass(iv1)

    if ((1.-f_cfcst)*f_mass(iv3).gt.f_mass(iv1-1) .and.  &
    (1.-f_cfcst)*f_mass(iv3).le.f_mass(iv1)) then
        f_weight=((1.-f_cfcst)*f_mass(iv3)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))
    elseif ((1.-f_cfcst)*f_mass(iv3).gt.f_mass(iv1) .and.  &
    (1.-f_cfcst)*f_mass(iv3).lt.f_mass(iv1+1)) then
        if (iv1.eq.nv_mud) then 
        f_weight=1.
        else
        f_weight=1.-(((1.-f_cfcst)*f_mass(iv3)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)))
        endif
    else
        f_weight=0.0
    endif 

    f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   &
       *((1.-f_cfcst)*f_mass(iv3))/f_mass(iv1)

endif ! end collision test case
enddo   
enddo
enddo

do iv1 = 1, nv_mud
    do iv2 = 1, nv_mud

        gcolfragiv1 = 2.*(Gval*(f_diam(iv1)+f_diam(iv2)))**2.*f_mass(iv1)*f_mass(iv2)  &
            /(pi*f_fy*f_fp*f_diam(iv1)**2.*(f_mass(iv1)+f_mass(iv2))         &
            *((f_rho(iv1)-rhoref)/rhoref)**(2./(3.-f_nf)))

        gcolfragiv2 = 2.*(Gval*(f_diam(iv1)+f_diam(iv2)))**2.*f_mass(iv1)*f_mass(iv2)  &
            /(pi*f_fy*f_fp*f_diam(iv2)**2.*(f_mass(iv1)+f_mass(iv2))         &
            *((f_rho(iv2)-rhoref)/rhoref)**(2./(3.-f_nf)))    

        mult = 1.
        if (iv1.eq.iv2) mult = 2.

        if (iv1 .eq. MAX(iv1,iv2) .and. gcolfragiv1 .ge. 1.) then
            f_l4(iv2, iv1) = f_l4(iv2, iv1)+ mult * (f_coll_prob_sh(iv1, iv2))
        elseif (iv1 .eq. MIN(iv1,iv2) .and. gcolfragiv2 .ge. 1.) then
            f_l4(iv2, iv1) = f_l4(iv2, iv1)+ mult * (f_coll_prob_sh(iv1, iv2))
        endif

    enddo
enddo

f_g4(1:nv_mud, 1:nv_mud, 1:nv_mud) = f_g4(1:nv_mud, 1:nv_mud, 1:nv_mud) * f_collfragparam
f_l4(1:nv_mud, 1:nv_mud) = f_l4(1:nv_mud, 1:nv_mud) * f_collfragparam

end subroutine flocmod_collfrag

!==================================================================================================
subroutine flocmod_mass_control(NN,mneg)
! Purpose : Compute mass in every class after flocculation and 
! returns negative mass if any
!
! Called by : flocmod_main
implicit none

real, dimension(1:nv_mud), intent(in)     :: NN
real, intent(out)     :: mneg

integer :: iv1

mneg = 0.0

do iv1 = 1, nv_mud
    if (NN(iv1) .lt. 0.0) then
        mneg = mneg - NN(iv1) * f_mass(iv1)
    endif
enddo

end subroutine flocmod_mass_control

!==================================================================================================
subroutine flocmod_mass_redistribute(NN)
! Purpose : based on a tolerated negative mass parameter, negative masses  
! are redistributed linearly towards remaining postive masses 
! and negative masses are set to 0
!
! Called by : flocmod_main
implicit none

real, dimension(1:nv_mud), intent(inout) :: NN

integer  :: iv
real     :: npos
real     :: mneg
real,dimension(1:nv_mud) :: NNtmp

mneg = 0.0
npos = 0.0
NNtmp(:) = NN(:)

do iv = 1, nv_mud
    if (NN(iv) .lt. 0.0) then
        mneg = mneg - NN(iv) * f_mass(iv)
        NNtmp(iv) = 0.0
    else
        npos = npos + 1.
    endif
enddo

if (mneg.gt.0.0) then 
    if (npos .eq. 0.) then
        write(ierrorlog,*) 'CAUTION : all floc sizes have negative mass!'
        write(ierrorlog,*) 'SIMULATION STOPPED'
        stop    
    else
        do iv=1,nv_mud
            if (NN(iv) .gt. 0.0) then
                NN(iv) = NN(iv) - mneg / sum(NNtmp) * NN(iv) / f_mass(iv)
            else
                NN(iv) = 0.0
            endif
        enddo
    endif
endif

return
end subroutine flocmod_mass_redistribute

!==================================================================================================  
end module flocmod
