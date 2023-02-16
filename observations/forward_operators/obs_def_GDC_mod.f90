! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! Fortran has a limit of 32 characters for variable names. Hence,
! each column can be at most 32 characters wide.
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
! BEGIN DART PREPROCESS TYPE DEFINITIONS
! GDC_TEMPERATURE,                 QTY_TEMPERATURE,            COMMON_CODE
! GDC_TEMPERATURE_ION,             QTY_TEMPERATURE_ION,        COMMON_CODE
! GDC_VELOCITY_U,                  QTY_VELOCITY_U,             COMMON_CODE
! GDC_VELOCITY_V,                  QTY_VELOCITY_V,             COMMON_CODE
! GDC_DENSITY_ION_OP,              QTY_DENSITY_ION_OP,         COMMON_CODE
! GDC_DENSITY_NEUTRAL_O,           QTY_DENSITY_NEUTRAL_O
! GDC_DENSITY_NEUTRAL_O2,          QTY_DENSITY_NEUTRAL_O2
! GDC_DENSITY_NEUTRAL_N2,          QTY_DENSITY_NEUTRAL_N2
! GDC_DENSITY,                     QTY_DENSITY
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_GDC_mod, only : get_expected_GDC_upper_atm_density
!  use obs_def_GDC_mod, only : get_expected_GDC_atomic_oxygen_number_density
!  use obs_def_GDC_mod, only : get_expected_GDC_molecular_oxygen_number_density
!  use obs_def_GDC_mod, only : get_expected_GDC_molecular_nitrogen_number_density
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!       case(GDC_DENSITY) 
!            call get_expected_GDC_upper_atm_density(state_handle, ens_size, location, expected_obs, istatus)
!       case(GDC_DENSITY_NEUTRAL_O)
!            call get_expected_GDC_atomic_oxygen_number_density(state_handle, ens_size, location, expected_obs, istatus)
!       case(GDC_DENSITY_NEUTRAL_O2)
!            call get_expected_GDC_molecular_oxygen_number_density(state_handle, ens_size, location, expected_obs, istatus)
!       case(GDC_DENSITY_NEUTRAL_N2)
!            call get_expected_GDC_molecular_nitrogen_number_density(state_handle, ens_size, location, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF



! BEGIN DART PREPROCESS READ_OBS_DEF
! case(GDC_DENSITY) 
!      continue
! case(GDC_DENSITY_NEUTRAL_O) 
!      continue
! case(GDC_DENSITY_NEUTRAL_O2)
!      continue
! case(GDC_DENSITY_NEUTRAL_N2)
!      continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
! case(GDC_DENSITY) 
!      continue
! case(GDC_DENSITY_NEUTRAL_O) 
!      continue
! case(GDC_DENSITY_NEUTRAL_O2)
!      continue
! case(GDC_DENSITY_NEUTRAL_N2)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
! case(GDC_DENSITY) 
!      continue
! case(GDC_DENSITY_NEUTRAL_O) 
!      continue
! case(GDC_DENSITY_NEUTRAL_O2)
!      continue
! case(GDC_DENSITY_NEUTRAL_N2)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_GDC_mod

use        types_mod, only : r8, MISSING_R8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, get_location, set_location, &
                             VERTISHEIGHT, VERTISLEVEL
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : QTY_DENSITY_NEUTRAL_O, &
                             QTY_DENSITY_NEUTRAL_O2, &
                             QTY_DENSITY_NEUTRAL_N2, &
                             QTY_ATOMIC_OXYGEN_MIXING_RATIO, &
                             QTY_ATOMIC_H_MIXING_RATIO, &
                             QTY_ION_O_MIXING_RATIO, &
                             QTY_MOLEC_OXYGEN_MIXING_RATIO, &
                             QTY_TEMPERATURE, &
                             QTY_PRESSURE, &
                             QTY_DENSITY, &
                             QTY_DENSITY_ION_E, &
                             QTY_ELECTRON_DENSITY, &
                             QTY_GEOPOTENTIAL_HEIGHT, &
                             QTY_GEOMETRIC_HEIGHT, &
                             QTY_O_N2_COLUMN_DENSITY_RATIO
use  ensemble_manager_mod, only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none
private
public :: get_expected_GDC_upper_atm_density, &
          get_expected_GDC_molecular_nitrogen_number_density, &
          get_expected_GDC_molecular_oxygen_number_density, &
          get_expected_GDC_atomic_oxygen_number_density
! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

real(r8), parameter :: N2_molar_mass = 28.0_r8
real(r8), parameter :: O_molar_mass  = 16.0_r8
real(r8), parameter :: O2_molar_mass = 32.0_r8
real(r8), parameter :: H_molar_mass  =  1.0_r8

! WACCM-X; put into common/types_mod.f90?
real(r8), parameter :: kboltz = 1.380648E-23_r8    ! [N*m/K]
real(r8), parameter :: universal_gas_constant = 8314.0_r8 ! [J/K/kmol]
real(r8), parameter :: molar_mass_dry_air = 28.9644_r8
integer,  parameter :: MAXLEVELS = 300 ! more than max levels expected in the model (waccm-x has 126)
character(len=512) :: string1, string2, string3

contains

!-----------------------------------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

!-----------------------------------------------------------------------------
!>@todo Test RMA.
! Given DART state vector and a location, 
! it computes thermospheric neutral density [Kg/m3] 
! The istatus variable should be returned as 0 unless there is a problem

subroutine get_expected_GDC_upper_atm_density(state_handle, ens_size, location, obs_val, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
real(r8),           intent(out) :: obs_val(ens_size)
integer,            intent(out) :: istatus(ens_size)

real(r8) :: mmro1(ens_size), mmro2(ens_size), mmrh1(ens_size) ! mass mixing ratio 
real(r8) :: mass_reciprocal(ens_size), pressure(ens_size), temperature(ens_size)
integer  :: this_istatus(ens_size)
logical  :: return_now

if ( .not. module_initialized ) call initialize_module

istatus = 0

! Some models (i.e. GITM) have density as part of the state.
! If it is available, just return it. If density is not state,
! then we need to create it from its constituents.

call interpolate(state_handle, ens_size, location, QTY_DENSITY, obs_val, istatus)
if(any(istatus == 0)) return ! density is part of the state


! This part was implemented for TIEGCM. Check the units for use with
! other models.
istatus(:) = 0
call interpolate(state_handle, ens_size, location, QTY_ATOMIC_OXYGEN_MIXING_RATIO, mmro1, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_MOLEC_OXYGEN_MIXING_RATIO, mmro2, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_ATOMIC_H_MIXING_RATIO, mmrh1, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_PRESSURE, pressure, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_TEMPERATURE, temperature, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return


! density [Kg/m3] =  pressure [N/m2] * M [g/mol] / temperature [K] / R [N*m/K/kmol] 
! where M is the mean molar mass 
! 1/M = sum(wi/Mi) where wi are mass mixing fractions and Mi are individual molar masses

where (istatus == 0) 
   mass_reciprocal = mmro1/O_molar_mass + mmro2/O2_molar_mass + mmrh1/H_molar_mass +&
                    (1.0_r8-mmro1-mmro2-mmrh1)/N2_molar_mass

   obs_val = pressure / mass_reciprocal / temperature / universal_gas_constant 
endwhere

end subroutine get_expected_GDC_upper_atm_density

!-----------------------------------------------------------------------------
!> Given DART state vector and a location, it computes O number density density [1/cm^3].

subroutine get_expected_GDC_atomic_oxygen_number_density(state_handle, ens_size, location, obs_val, istatus)
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(out) :: istatus(ens_size)
real(r8),            intent(out) :: obs_val(ens_size)

real(r8), dimension(ens_size)  :: mmr_o1, mmr_o2, mmr_n2, mmr_h1, mmr_op   ! mass mixing ratio 
real(r8), dimension(ens_size)  :: mbar, pressure, temperature 
integer,  dimension(ens_size)  :: this_istatus
real(r8), dimension(3)  :: loc_vals
logical :: return_now

istatus = 0 ! Need to have istatus = 0 for track_status()
call interpolate(state_handle, ens_size, location, QTY_DENSITY_NEUTRAL_O, obs_val, istatus)
if(any(istatus == 0)) return ! density is part of the state

! cam-fv returns volume mixing ratio, not mass mixing ratio. undo for computation below.
call interpolate(state_handle, ens_size, location, QTY_ATOMIC_OXYGEN_MIXING_RATIO, mmr_o1, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return
mmr_o1 = mmr_o1 / (molar_mass_dry_air/O_molar_mass)

call interpolate(state_handle, ens_size, location, QTY_MOLEC_OXYGEN_MIXING_RATIO, mmr_o2, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return
mmr_o2 = mmr_o2 / (molar_mass_dry_air/O2_molar_mass)

call interpolate(state_handle, ens_size, location, QTY_ATOMIC_H_MIXING_RATIO, mmr_h1, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return
mmr_h1 = mmr_h1 / (molar_mass_dry_air/H_molar_mass)

call interpolate(state_handle, ens_size, location, QTY_PRESSURE, pressure, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_TEMPERATURE, temperature, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return

!------------------------------------------------------------------------------------------------------
!  Need to get number density (cgs units) from mass mixing ratio (kg/kg).  
!  mbar is g/mole, same as rMass units
!       kg/kg * (g/mole)/(g/mole) * (Pa = N/m^2)/((Joules/K = N*m/K) * (K)) = m-3 * 1E-06 = cm-3
!------------------------------------------------------------------------------------------------------

loc_vals = get_location(location)

where (istatus == 0) 
   mmr_n2 = 1.0_r8 - (mmr_o1 + mmr_o2 + mmr_h1)
   mbar   = 1.0_r8/( mmr_o1/O_molar_mass   &
                   + mmr_o2/O2_molar_mass  &
                   + mmr_h1/H_molar_mass   &
                   + mmr_n2/N2_molar_mass)
   obs_val = mmr_o1 * mbar/O_molar_mass * pressure/(kboltz * temperature) * 1.E-06_r8
end where

end subroutine get_expected_GDC_atomic_oxygen_number_density

!------------------------------------------------------------------------------------------------------
subroutine get_expected_GDC_molecular_oxygen_number_density(state_handle, ens_size, location, obs_val, istatus)
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(out) :: istatus(ens_size)
real(r8),            intent(out) :: obs_val(ens_size)

real(r8), dimension(ens_size)  :: mmr_o1, mmr_o2, mmr_n2, mmr_h1, mmr_op   ! mass mixing ratio 
real(r8), dimension(ens_size)  :: mbar, pressure, temperature 
integer,  dimension(ens_size)  :: this_istatus
real(r8), dimension(3)  :: loc_vals
logical :: return_now

istatus = 0 ! Need to have istatus = 0 for track_status()
call interpolate(state_handle, ens_size, location, QTY_DENSITY_NEUTRAL_O2, obs_val, istatus)
if(any(istatus == 0)) return ! density is part of the state

! cam-fv returns volume mixing ratio, not mass mixing ratio. undo for computation below.
call interpolate(state_handle, ens_size, location, QTY_ATOMIC_OXYGEN_MIXING_RATIO, mmr_o1, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return
mmr_o1 = mmr_o1 / (molar_mass_dry_air/O_molar_mass)

call interpolate(state_handle, ens_size, location, QTY_MOLEC_OXYGEN_MIXING_RATIO, mmr_o2, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return
mmr_o2 = mmr_o2 / (molar_mass_dry_air/O2_molar_mass)

call interpolate(state_handle, ens_size, location, QTY_ATOMIC_H_MIXING_RATIO, mmr_h1, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return
mmr_h1 = mmr_h1 / (molar_mass_dry_air/H_molar_mass)

call interpolate(state_handle, ens_size, location, QTY_PRESSURE, pressure, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_TEMPERATURE, temperature, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return

!------------------------------------------------------------------------------------------------------
!  Need to get number density (cgs units) from mass mixing ratio (kg/kg).  
!  mbar is g/mole, same as rMass units
!       kg/kg * (g/mole)/(g/mole) * (Pa = N/m^2)/((Joules/K = N*m/K) * (K)) = m-3 * 1E-06 = cm-3
!------------------------------------------------------------------------------------------------------

loc_vals = get_location(location)

where (istatus == 0) 
   mmr_n2 = 1.0_r8 - (mmr_o1 + mmr_o2 + mmr_h1)
   mbar   = 1.0_r8/( mmr_o1/O_molar_mass   &
                   + mmr_o2/O2_molar_mass  &
                   + mmr_h1/H_molar_mass   &
                   + mmr_n2/N2_molar_mass)
   obs_val = mmr_o2 * mbar/O2_molar_mass * pressure/(kboltz * temperature) * 1.E-06_r8
end where

end subroutine get_expected_GDC_molecular_oxygen_number_density

!------------------------------------------------------------------------------------------------------
subroutine get_expected_GDC_molecular_nitrogen_number_density(state_handle, ens_size, location, obs_val, istatus)
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(out) :: istatus(ens_size)
real(r8),            intent(out) :: obs_val(ens_size)

real(r8), dimension(ens_size)  :: mmr_o1, mmr_o2, mmr_n2, mmr_h1, mmr_op   ! mass mixing ratio 
real(r8), dimension(ens_size)  :: mbar, pressure, temperature 
integer,  dimension(ens_size)  :: this_istatus
real(r8), dimension(3)  :: loc_vals
logical :: return_now

istatus = 0 ! Need to have istatus = 0 for track_status()
call interpolate(state_handle, ens_size, location, QTY_DENSITY_NEUTRAL_N2, obs_val, istatus)
if(any(istatus == 0)) return ! density is part of the state

! cam-fv returns volume mixing ratio, not mass mixing ratio. undo for computation below.
call interpolate(state_handle, ens_size, location, QTY_ATOMIC_OXYGEN_MIXING_RATIO, mmr_o1, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return
mmr_o1 = mmr_o1 / (molar_mass_dry_air/O_molar_mass)

call interpolate(state_handle, ens_size, location, QTY_MOLEC_OXYGEN_MIXING_RATIO, mmr_o2, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return
mmr_o2 = mmr_o2 / (molar_mass_dry_air/O2_molar_mass)

call interpolate(state_handle, ens_size, location, QTY_ATOMIC_H_MIXING_RATIO, mmr_h1, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return
mmr_h1 = mmr_h1 / (molar_mass_dry_air/H_molar_mass)

call interpolate(state_handle, ens_size, location, QTY_PRESSURE, pressure, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_TEMPERATURE, temperature, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return

!------------------------------------------------------------------------------------------------------
!  Need to get number density (cgs units) from mass mixing ratio (kg/kg).  
!  mbar is g/mole, same as rMass units
!       kg/kg * (g/mole)/(g/mole) * (Pa = N/m^2)/((Joules/K = N*m/K) * (K)) = m-3 * 1E-06 = cm-3
!------------------------------------------------------------------------------------------------------

loc_vals = get_location(location)

where (istatus == 0) 
   mmr_n2 = 1.0_r8 - (mmr_o1 + mmr_o2 + mmr_h1)
   mbar   = 1.0_r8/( mmr_o1/O_molar_mass   &
                   + mmr_o2/O2_molar_mass  &
                   + mmr_h1/H_molar_mass   &
                   + mmr_n2/N2_molar_mass)
   obs_val = mmr_n2 * mbar/N2_molar_mass * pressure/(kboltz * temperature) * 1.E-06_r8
end where

end subroutine get_expected_GDC_molecular_nitrogen_number_density

end module obs_def_GDC_mod
! END DART PREPROCESS MODULE CODE      

