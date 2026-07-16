!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS) Coupler.
!*
!* FMS Coupler is free software: you can redistribute it and/or modify
!* it under the terms of the GNU Lesser General Public License as
!* published by the Free Software Foundation, either version 3 of the
!* License, or (at your option) any later version.
!*
!* FMS Coupler is distributed in the hope that it will be useful, bu
!* WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!* General Public License for more details.
 !*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS Coupler.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!! @file
!> @parblock
!! Module coupler_mod contains the public coupler driver routines for the
!! fully coupled GFDL climate model (atmosphere, land, sea ice, and ocean).
!!
!! Each routine in this module is a wrapper around a "science" routine.
!! The wrappers add:
!! - fms_mpp_clock_begin and end to measure performance
!! - Optional checksum computation if do_chksum namelist flag is .true.
!! - Optional memory-usage reporting if do_debug namelist flag is .true.
!! - mpp_set_current_pelist calls for MPI synchronization.
!! @endparblock
module full_coupler_mod

  use omp_lib

  use FMS
  use FMSconstants, only: fmsconstants_ini

#ifdef use_deprecated_io
  use fms_io_mod, only: fms_io_exi
#endif

  use atmos_model_mod, only: atmos_model_init, atmos_model_end
  use atmos_model_mod, only: update_atmos_model_dynamics
  use atmos_model_mod, only: update_atmos_model_down
  use atmos_model_mod, only: update_atmos_model_up
  use atmos_model_mod, only: atmos_data_type
  use atmos_model_mod, only: land_ice_atmos_boundary_type
  use atmos_model_mod, only: atmos_data_type_chksum
  use atmos_model_mod, only: lnd_ice_atm_bnd_type_chksum
  use atmos_model_mod, only: lnd_atm_bnd_type_chksum
  use atmos_model_mod, only: ice_atm_bnd_type_chksum
  use atmos_model_mod, only: atmos_model_restar
  use atmos_model_mod, only: update_atmos_model_radiation
  use atmos_model_mod, only: update_atmos_model_state

  use land_model_mod, only: land_model_init, land_model_end
  use land_model_mod, only: land_data_type, atmos_land_boundary_type
  use land_model_mod, only: update_land_model_fast, update_land_model_slow
  use land_model_mod, only: atm_lnd_bnd_type_chksum
  use land_model_mod, only: land_data_type_chksum
  use land_model_mod, only: land_model_restar

  use ice_model_mod, only: ice_model_init, share_ice_domains, ice_model_end, ice_model_restar
  use ice_model_mod, only: update_ice_model_fast, set_ice_surface_fields
  use ice_model_mod, only: ice_data_type, land_ice_boundary_type
  use ice_model_mod, only: ocean_ice_boundary_type, atmos_ice_boundary_type
  use ice_model_mod, only: ice_data_type_chksum, ocn_ice_bnd_type_chksum
  use ice_model_mod, only: atm_ice_bnd_type_chksum, lnd_ice_bnd_type_chksum
  use ice_model_mod, only: unpack_ocean_ice_boundary, exchange_slow_to_fast_ice
  use ice_model_mod, only: ice_model_fast_cleanup, unpack_land_ice_boundary
  use ice_model_mod, only: exchange_fast_to_slow_ice, update_ice_model_slow

  use ocean_model_mod, only: update_ocean_model, ocean_model_init,  ocean_model_end
  use ocean_model_mod, only: ocean_public_type, ocean_state_type, ice_ocean_boundary_type
  use ocean_model_mod, only: ocean_model_restar
  use ocean_model_mod, only: ocean_public_type_chksum, ice_ocn_bnd_type_chksum

  use combined_ice_ocean_driver, only: update_slow_ice_and_ocean, ice_ocean_driver_type
  use combined_ice_ocean_driver, only: ice_ocean_driver_init, ice_ocean_driver_end
!
! flux_ calls translate information between model grids - see flux_exchange.f90
!

  use flux_exchange_mod, only: flux_exchange_init, gas_exchange_init, sfc_boundary_layer
  use flux_exchange_mod, only: generate_sfc_xgrid, send_ice_mask_sic
  use flux_exchange_mod, only: flux_down_from_atmos, flux_up_to_atmos
  use flux_exchange_mod, only: flux_land_to_ice, flux_ice_to_ocean, flux_ocean_to_ice
  use flux_exchange_mod, only: flux_ice_to_ocean_finish, flux_ocean_to_ice_finish
  use flux_exchange_mod, only: flux_check_stocks, flux_init_stocks
  use flux_exchange_mod, only: flux_ocean_from_ice_stocks, flux_ice_to_ocean_stocks
  use flux_exchange_mod, only: flux_atmos_to_ocean, flux_ex_arrays_dealloc

  use atmos_tracer_driver_mod, only: atmos_tracer_driver_gather_data
  use gex_mod, only: gex_ini

  use iso_fortran_env

  implicit none
  private

  public :: atmos_data_type
  public :: atmos_ice_boundary_type
  public :: atmos_land_boundary_type
  public :: ice_data_type
  public :: ice_ocean_boundary_type
  public :: ice_ocean_driver_type
  public :: land_data_type
  public :: land_ice_atmos_boundary_type
  public :: land_ice_boundary_type
  public :: ocean_ice_boundary_type
  public :: ocean_public_type
  public :: ocean_state_type

  public :: fmsconstants_ini

  ! need to be made public in order to call from coupler_main.F90
  public :: flux_ice_to_ocean_finish
  public :: flux_ice_to_ocean_stocks
  public :: flux_ocean_from_ice_stocks
  public :: send_ice_mask_sic
  public :: update_slow_ice_and_ocean

  public :: atmos_model_restar
  public :: ice_model_restar
  public :: land_model_restar
  public :: ocean_model_restar

  public :: atm_ice_bnd_type_chksum
  public :: atm_lnd_bnd_type_chksum
  public :: atmos_data_type_chksum
  public :: ice_atm_bnd_type_chksum
  public :: ice_data_type_chksum
  public :: ice_ocn_bnd_type_chksum
  public :: land_data_type_chksum
  public :: lnd_atm_bnd_type_chksum
  public :: lnd_ice_atm_bnd_type_chksum
  public :: lnd_ice_bnd_type_chksum
  public :: ocean_public_type_chksum
  public :: ocn_ice_bnd_type_chksum

  public :: coupler_end
  public :: coupler_ini
  public :: coupler_intermediate_restar
  public :: coupler_restar
  public :: coupler_summarize_timestep

  public :: coupler_atmos_tracer_driver_gather_data
  public :: coupler_exchange_fast_to_slow_ice
  public :: coupler_exchange_slow_to_fast_ice
  public :: coupler_flux_atmos_to_ocean
  public :: coupler_flux_check_stocks
  public :: coupler_flux_down_from_atmos
  public :: coupler_flux_ice_to_ocean
  public :: coupler_flux_init_finish_stocks
  public :: coupler_flux_land_to_ice
  public :: coupler_flux_ocean_to_ice
  public :: coupler_flux_up_to_atmos
  public :: coupler_generate_sfc_xgrid
  public :: coupler_set_ice_surface_fields
  public :: coupler_sfc_boundary_layer
  public :: coupler_unpack_land_ice_boundary
  public :: coupler_unpack_ocean_ice_boundary
  public :: coupler_update_atmos_model_down
  public :: coupler_update_atmos_model_dynamics
  public :: coupler_update_atmos_model_radiation
  public :: coupler_update_atmos_model_state
  public :: coupler_update_atmos_model_up
  public :: coupler_update_ice_model_fas
  public :: coupler_update_ice_model_slow_and_stocks
  public :: coupler_update_land_model_fas
  public :: coupler_update_land_model_slow
  public :: coupler_update_ocean_model

  public :: coupler_clock_type
  public :: coupler_chksum_type
  public :: coupler_components_type

#include <file_version.fh>

  integer, dimension(6), public :: restart_interval = (/ 0, 0, 0, 0, 0, 0/)
  !< is a variable in the coupler namelist (coupler_nml) with format (yr, mo, day, hr, min, sec)
  !! to set the time interval for writing the intermediate restarts.
  !! Intermediate restarts are not written if restart_interval = [0,0,0,0,0,0]
  !! Example nml to write restarts every 6 hours:
  !! &coupler_nml
  !!   restart_interval = 0, 0, 0, 6, 0, 0
  !! /

  integer, dimension(6) :: current_date = (/ 0, 0, 0, 0, 0, 0 /)
  !> is a variable in the coupler namelist (coupler_nml) with forma
  !! (yr, mo, day, hr, min, sec) to set the model start date.
  !! Example — start a run on 1 January 2000:
  !! &coupler_nml
  !!   current_date = 2000, 1, 1, 0, 0, 0
  !! /

  character(len=17) :: calendar = '                 '
  !< is a variable in the coupler namelist (coupler_nml) to set the
  !! calendar type.  Valid values are those from FMS/time_manager_mod:
  !! 'gregorian', 'julian', 'noleap', or 'thirty_day'.

  logical :: force_date_from_namelist = .false.
  !> is a flag in the coupler namelist (coupler_nml) where a .true. value enforces
  !! starting date to be current_date from the namelist.

  integer, public :: months=0
    !< is a namelist variable to set the number of additional months to simulate
  integer, public :: days=0
    !< is a namelist variable to set the number of additional days to simulate
  integer, public :: hours=0
    !< is a namelist variable to set the number of additional hours to simulate
  integer, public :: minutes=0
    !< is a namelist variable to set the number of additional minutes to simulate
  integer, public :: seconds=0
    !< is a namelist variable to set the number of additional seconds to simulate
  integer, public :: dt_atmos = 0
    !< is a namelist variable to set the the time step [s] for the atmospheric model dynamics and
    !! fast coupling with land and sea ice
  integer, public :: dt_cpld  = 0
    !< is a namelist variable to set the time step [s] for coupling between ocean and atmosphere.  This mus
    !! be an integral multiple of dt_atmos and dt_ocean.  This is the "slow" timestep.
  integer, public :: atmos_npes=0
    !< is a namelist variable to set the number of MPI ranks (processing element) for the atmosphere
  integer, public :: ocean_npes=0
    !< is a namelist variable to set the number of MPI ranks (processing element) for the ocean
  integer, public :: ice_npes=0
    !< is a namelist variable to set the number of MPI ranks (processing element) for the ice
  integer, public :: land_npes=0
    !< is a namelist variable to set the number of MPI ranks (processing element) for the land
  integer, public :: atmos_nthreads=1
    !< is a namelist variable to set the number of OpenMP threads to use for the atmosphere
  integer, public :: ocean_nthreads=1
    !< is a namelist variable to set the number of OpenMP threads to use for the ocean
  integer, public :: radiation_nthreads=1
    !< is a namelist variable to set the number of threads to use for radiation;
    !! is set to atmos_nthreads if do_concurrent_radiation is .false.

  logical, public :: do_atmos =.true.
    !< is a namelist flag where if .false., skip atmosphere update
  logical, public :: do_land =.true.
    !< is a namelist flag where if .false., skip land model update
  logical, public :: do_ice =.true.
    !< is a namelist flag where if .false., skip sea-ice model update
  logical, public :: do_ocean=.true.
    !< is a namelist flag where if .false., skip ocean model update
  logical, public :: do_flux =.true.
    !< is a namelist flag where if .false., skip all flux exchanges between components.

  logical, public :: concurrent=.FALSE.
    !< is a namelist flag where if .TRUE., the ocean model updates concurrently with the atmosphere-land-ice
    !! on a separate set of PEs.  Concurrent should be .TRUE. if concurrent_ice is .TRUE.
    !! If .FALSE., the execution is serial: call atmos... followed by call ocean...

  logical, public :: do_concurrent_radiation=.FALSE.
    !< is a namelist flag where if .TRUE., radiation is updated concurrently with
    !! atmospheric physics with OpenMP threading (radiation_nthreads)

  logical, public :: use_lag_fluxes=.TRUE.
  !< is a namelist flag where if .TRUE., the ocean is forced with surface boundary conditions (SBCs)
  !! from one coupling timestep ago.  If .FALSE., the ocean is forced with most recent SBCs.
  !! For an old leapfrog MOM4 coupling with dt_cpld=dt_ocean, lag fluxes can be shown to be stable
  !! and current fluxes to be unconditionally unstable.  For dt_cpld>dt_ocean there
  !! is probably sufficient damping for MOM4.  For more modern ocean models (such as
  !! MOM5, GOLD or MOM6) that do not use leapfrog timestepping, use_lag_fluxes=.False.
  !! should be much more stable.  Controls whether ice-to-ocean flux exchange are called
  !! at the beginning or at the end of the time loop

  logical, public :: concurrent_ice=.FALSE.
  !< is a namelist flag where if .TRUE., the slow sea-ice is forced with the fluxes that were used for the
  !! fast ice processes one timestep before.  When used in conjuction with setting
  !! slow_ice_with_ocean=.TRUE., this approach allows the atmosphere and
  !! ocean to run concurrently even if use_lag_fluxes=.FALSE., and it can
  !! be shown to ameliorate or eliminate several ice-ocean coupled instabilities.

  logical, public :: slow_ice_with_ocean=.FALSE.
  !< is a namelist flag where if true, the slow sea-ice is advanced on the ocean PEs.  Otherwise
  !! the slow sea-ice processes are on the same PEs as the fast sea-ice.

  logical, public :: combined_ice_and_ocean=.FALSE.
  !< is a namelist flag where if true, there is a single call from the coupler to advance
  !! both the slow sea-ice and the ocean. slow_ice_with_ocean and
  !! concurrent_ice must both be true if combined_ice_and_ocean is true.

  logical, public :: do_chksum=.FALSE.
    !< is a namelist flag where if .TRUE., compute checksums throughout the model simulation
  logical, public :: do_endpoint_chksum=.TRUE.
    !< is a namelist flag where if .TRUE., do checksums of the initial and final states
    !! regardless of do_chksum value.
  logical, public :: do_debug=.FALSE.
    !< is a namelist flag where if .TRUE., print additional debugging messages.
  integer, public :: check_stocks = 0
    !< is a namelist flag where value of -1: don't compute stocks;
    !! value of 0: compute stocks at the end of the run; value = n > 0:  compute stocks every n coupled steps
  logical, public :: use_hyper_thread = .false.
    !< is a namelist flag where if .TRUE., enable use of hyperthreading on supported hardware

  namelist /coupler_nml/ current_date, calendar, force_date_from_namelist,         &
                         months, days, hours, minutes, seconds, dt_cpld, dt_atmos, &
                         do_atmos, do_land, do_ice, do_ocean, do_flux,             &
                         atmos_npes, ocean_npes, ice_npes, land_npes,              &
                         atmos_nthreads, ocean_nthreads, radiation_nthreads,       &
                         concurrent, do_concurrent_radiation, use_lag_fluxes,      &
                         check_stocks, restart_interval, do_debug, do_chksum,      &
                         use_hyper_thread, concurrent_ice, slow_ice_with_ocean,    &
                         do_endpoint_chksum, combined_ice_and_ocean

  !> coupler_clock_type is a derived type that holds all FMS performance-clock IDs
  !! used for profiling the coupled model time loop.  Each integer member is an
  !! clock id returned by fms_mpp_clock_id that is used in fms_mpp_clock_begin/end
  type coupler_clock_type
    integer :: atm !< Outer clock enclosing all atmosphere-side work each coupled step
    integer :: atmos_loop !< Clock for the inner atmospheric sub-step loop
    integer :: atmos_model_init !< Clock for atmospheric model initialization
    integer :: atmos_tracer_driver_gather_data !< Clock for gathering bottom-level tracers before flux exchange
    integer :: concurrent_atmos !< Clock for the concurrent-radiation atmospheric physics thread
    integer :: final_flux_check_stocks !< Clock for the end-of-run stock conservation check
    integer :: flux_check_stocks !< Clock for periodic mid-run stock conservation checks
    integer :: flux_down_from_atmos !< Clock for mapping atmosphere→land/ice fluxes via the exchange grid
    integer :: flux_exchange_init !< Clock for flux-exchange infrastructure initialization
    integer :: flux_ice_to_ocean !< Clock for packaging ice fluxes for the ocean
    integer :: flux_ice_to_ocean_stocks !< Clock for bookkeeping water/heat/salt stocks transferred ice→ocean
    integer :: flux_land_to_ice !< Clock for transferring land runoff and calving to the ice grid
    integer :: flux_ocean_to_ice !< Clock for transferring ocean state to the ice model
    integer :: flux_up_to_atmos !< Clock for mapping updated surface state back to the atmosphere
    integer :: generate_sfc_xgrid !< Clock for rebuilding the surface exchange grid
    integer :: ice_model_init !< Clock for sea-ice model initialization
    integer :: initialization !< Clock for the full coupler_init phase
    integer :: intermediate_restart !< Clock for writing intermediate restart files
    integer :: land_model_init !< Clock for land model initialization
    integer :: main !< Clock for the main coupled time loop
    integer :: ocean !< Clock for the ocean model update
    integer :: ocean_model_init !< Clock for ocean model initialization
    integer :: radiation !< Clock for the atmospheric radiation calculation (serial or concurrent)
    integer :: set_ice_surface_exchange !< Clock for the fast↔slow ice state exchange
    integer :: set_ice_surface_fast !< Clock for computing fast-ice surface fields (albedo, roughness, etc.)
    integer :: set_ice_surface_slow !< Clock for computing slow-ice surface fields after ocean-ice exchange
    integer :: sfc_boundary_layer !< Clock for the surface boundary-layer turbulent-flux computation
    integer :: termination !< Clock for the coupler_end finalization phase
    integer :: update_atmos_model_down !< Clock for the downward atmospheric physics sweep
    integer :: update_atmos_model_dynamics  !< Clock for the atmospheric dynamical core update
    integer :: update_atmos_model_state !< Clock for committing tendencies and advancing the atmospheric clock
    integer :: update_atmos_model_up  !< Clock for the upward atmospheric physics sweep
    integer :: update_ice_model_fast !< Clock for the fast (atmospheric timestep) sea-ice update
    integer :: update_ice_model_slow_exchange !< Clock for exchange_fast_to_slow_ice (concurrent-ice mode)
    integer :: update_ice_model_slow_fast  !< Clock for fast-PE portion of the slow-ice update (cleanup + unpack)
    integer :: update_ice_model_slow_slow !< Clock for the slow sea-ice thermodynamics/dynamics update
    integer :: update_land_model_fast !< Clock for the fast (atmospheric timestep) land model update
    integer :: update_land_model_slow !< Clock for the slow (coupled timestep) land model update
  end type coupler_clock_type

  !> coupler_components_type is a convenient object that holds pointers to all
  !! model component derived types and inter-component boundary derived types.
  !! Its primary purpose is to reduce the length of the argument list when calling
  !! checksum subroutines that only need read-only access to all components.
  type coupler_components_type
    private
    type(atmos_data_type), pointer :: Atm !< Pointer to the atmosphere component derived type
    type(land_data_type), pointer :: Land !< Pointer to the land component derived type
    type(ice_data_type), pointer :: Ice !< Pointer to the sea-ice component derived type
    type(ocean_public_type), pointer :: Ocean !< Pointer to the ocean component public derived type
    type(land_ice_atmos_boundary_type), pointer :: Land_ice_atmos_boundary !< Pointer to the
                                                                            !! land-ice to atm boundary fluxes
    type(atmos_land_boundary_type), pointer :: Atmos_land_boundary !< Pointer to the atm to land boundary fluxes
    type(atmos_ice_boundary_type), pointer :: Atmos_ice_boundary !< Pointer to the atm to ice boundary fluxes
    type(land_ice_boundary_type), pointer :: Land_ice_boundary !< Pointer to the land to ice boundary
    type(ice_ocean_boundary_type), pointer :: Ice_ocean_boundary !< Pointer to the ice to ocean boundary fluxes
    type(ocean_ice_boundary_type), pointer :: Ocean_ice_boundary !< Pointer to the ocean to ice boundary state

  contains
    procedure, public :: initialize_coupler_components_obj !< Associates all pointer members to the components
    procedure, public :: get_component  !< Retrieves a pointer to a named component from this objec
  end type coupler_components_type

  !> coupler_chksum_type contains all checksum-related operations for the coupler.
  !! It holds a pointer to a coupler_components_type object that can be retrieved
  !! during checksum computation.  This object was created to avoid passing a long
  !! list of component derived types to the checksum subroutines
  type coupler_chksum_type
    private
    type(coupler_components_type), pointer :: components !< Pointer to the container of all component derived types
  contains
    procedure, public :: initialize_coupler_chksum_obj  !< Associates components pointer to an initialized
                                                        !! coupler_components_type
    procedure, public :: get_components_obj !< Returns a copy of the coupler_components_type pointed to by this objec
    procedure, public :: get_atmos_ice_land_ocean_chksums !< Computes checksums for all four components
                                                          !! (atm, land, ice, ocean)
    procedure, public :: get_atmos_ice_land_chksums !< Computes checksums for atmosphere, land, and fast-ice fields
    procedure, public :: get_slow_ice_chksums !< Computes checksums for slow-ice and ocean-ice boundary fields
    procedure, public :: get_ocean_chksums !< Computes checksums for ocean and ice-ocean boundary fields
    procedure, public :: get_coupler_chksums !< Computes checksums for selected atmosphere, land, and ice state fields
  end type coupler_chksum_type

  character(len=80) :: tex
    !< is a temporary string variable used for printing messages
  character(len=48), parameter :: mod_name = 'coupler_main_mod'
    !< is a string variable for the module name, used in error messages

  integer :: calendar_type = INVALID_CALENDAR
    !< is the calendar_type initialized to INVALID_CALENDAR

  integer :: date_init(6) = (/ 0, 0, 0, 0, 0, 0 /)
    !< is the initial date for the model run in (yr, mo, day, hr, min, sec) forma

  contains

  !> @parblock
  !! Subroutine coupler_init initializes all component models and the flux-exchange infrastructure,
  !! and sets all runtime configurations.  Coupler_init must be called before the time stepping loops.
  !! @endparblock
  subroutine coupler_init(Atm, Ocean, Land, Ice, Ocean_state, Atmos_land_boundary, Atmos_ice_boundary, &
      Ocean_ice_boundary, Ice_ocean_boundary, Land_ice_atmos_boundary, Land_ice_boundary,              &
      Ice_ocean_driver_CS, Ice_bc_restart, Ocn_bc_restart, ensemble_pelist, slow_ice_ocean_pelist, conc_nthreads, &
      coupler_clocks, coupler_components_obj, coupler_chksum_obj, Time_step_cpld, Time_step_atmos, Time_atmos, &
      Time_ocean, num_cpld_calls, num_atmos_calls, Time, Time_start, Time_end, Time_restart, Time_restart_current)

    implicit none

    type(atmos_data_type), intent(inout) :: Atm
      !< is the atmosphere derived type
    type(land_data_type), intent(inout) :: Land
      !< is the land derived type
    type(ice_data_type), intent(inout) :: Ice
      !< is the sea-ice derived type
    type(ocean_public_type), intent(inout) :: Ocean
      !< is the ocean public derived type
    type(ocean_state_type), pointer, intent(inout) :: Ocean_state
      !< is the ocean internal state derived type
    type(atmos_land_boundary_type), intent(inout) :: Atmos_land_boundary
      !< is the derived type holding atmos to land fluxes and properies
    type(atmos_ice_boundary_type), intent(inout) :: Atmos_ice_boundary
      !< is the derived type holding atmos to ice fluxes and properties
    type(ice_ocean_boundary_type), intent(inout) :: Ice_ocean_boundary
      !< is the derived type holding ice to ocean fluxes and properties
    type(ocean_ice_boundary_type), intent(inout) :: Ocean_ice_boundary
      !< is the derived type holding ocean to ice fluxes and properties
    type(land_ice_boundary_type), intent(inout) :: Land_ice_boundary
      !< is the derived type holding land to ice fluxes and properties
    type(ice_ocean_driver_type), pointer, intent(inout) :: Ice_ocean_driver_CS
      !< is the control structure for combined ice-ocean driver
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_ice_atmos_boundary
      !< is the derived type holding atmos to/from land and ice fluxes and properties
    type(FmsNetcdfDomainFile_t), pointer, dimension(:), intent(inout) :: Ice_bc_restar
      !< is the fms2_io file object to write ice boundary condition restarts
    type(FmsNetcdfDomainFile_t), pointer, dimension(:), intent(inout) :: Ocn_bc_restar
      !< is the fms2_io file object to write ocean boundary condition restarts

    integer, intent(inout) :: conc_nthreads
      !< is the number of concurrent OpenMP threads; set to 2 when do_concurrent_radiation=.true.
    integer, allocatable, dimension(:,:), intent(inout) :: ensemble_pelis
      !< is the PE list for each ensemble member (ensemble_size, npes) for an ensemble run
    integer, allocatable, dimension(:),   intent(inout) :: slow_ice_ocean_pelis
      !< is the union of slow-ice and ocean PE lists

    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< is the collection of FMS clock IDs for profiling
    type(coupler_components_type), intent(inout) :: coupler_components_obj
      !< is the object holding pointers to all component model structures
    type(coupler_chksum_type), intent(inout) :: coupler_chksum_obj
     !< is the object used for computing and printing checksums

    type(FMSTime_type), intent(inout) :: Time_step_cpld
      !< is the coupled (slow) time step
    type(FMSTime_type), intent(inout) :: Time_step_atmos
      !< is the atmosphere (fast) time step
    type(FMSTime_type), intent(inout) :: Time_atmos
      !< is the current atmosphere model time
    type(FMSTime_type), intent(inout) :: Time_ocean
      !< is the current ocean model time
    type(FMSTime_type), intent(inout) :: Time
      !< is the current model time
    type(FMSTime_type), intent(inout) :: Time_star
      !< is the model start time
    type(FMSTime_type), intent(inout) :: Time_end
      !< is the model end time
    type(FMSTime_type), intent(inout) :: Time_restar
      !< is the time of the next intermediate restart write
    type(FMSTime_type), intent(inout) :: Time_restart_curren
      !< is the time of the most recent intermediate restart write

    integer, intent(inout) :: num_cpld_calls
      !< is the total number of coupled (slow) time steps in this run
    integer, intent(inout) :: num_atmos_calls
      !< is the number of atmosphere time steps per coupled time step

    character(len=64), parameter :: sub_name = 'coupler_init'
    character(len=256), parameter:: error_header = &
         '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
    character(len=256), parameter :: note_header = &
         '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

    integer :: ierr, io, m, i, outunit, logunit, erruni
    integer :: date(6) ! is the current model date as (year, month, day, hour, minute, second)
    type (FmsTime_type) :: Run_length ! is the total run length as a FMSTime_type
    character(len=9) :: month  ! is the name of the current month (for log messages)
    integer :: pe, npes ! is the current PE rank and total number of PEs in the ensemble member

    integer :: ens_siz(6), ensemble_size ! is the ensemble sizing array and scalar ensemble size
    integer :: ensemble_id = 1  ! is the ID (index) of this ensemble member (1-based)

    integer :: atmos_pe_start=0, atmos_pe_end=0, & ! First and last PE indices of the atmosphere PE range
               ocean_pe_start=0, ocean_pe_end=0 ! First and last PE indices of the ocean PE range
    integer :: n ! General-purpose loop or count index
    integer :: diag_model_subset=DIAG_ALL ! Diagnostic subset flag passed to fms_diag_ini
    logical :: other_fields_exist ! Scratch flag used when checking for optional restart fields
    character(len=256) :: err_msg ! Error message string returned from FMS routines
    integer :: date_restart(6) ! Date of the most recent intermediate restart as (yr,mo,day,hr,min,sec)
    character(len=64)  :: filename, fieldname ! Scratch strings for restart file and field names
    integer :: id_restart, l ! Loop indices for restart file processing
    character(len=8)  :: walldate ! Wall-clock date string from DATE_AND_TIME
    character(len=10) :: walltime ! Wall-clock time string from DATE_AND_TIME
    character(len=5)  :: wallzone ! Wall-clock time-zone string from DATE_AND_TIME
    integer :: wallvalues(8) ! Wall-clock integer values from DATE_AND_TIME
    character(len=:), dimension(:), allocatable :: restart_file ! Restart file saved as a string
    integer :: time_stamp_unit ! Unit of the time_stamp file
    integer :: ascii_unit  ! Unit of a dummy ascii file

    type(FmsTime_type) :: Time_ini

    type(FmsCoupler1dBC_type), pointer :: gas_fields_atm => NULL()
      ! A pointer to the type containing the atmospheric gas fields
    type(FmsCoupler1dBC_type), pointer :: gas_fields_ocn => NULL()
     ! A pointer to the type containing the ocean and ice surface gas fields
    type(FmsCouplerGasFluxes_type), pointer :: gas_fluxes => NULL()
      ! A pointer to the type containing the atmosphere-ocean gas and tracer fluxes.

    integer :: num_ice_bc_restart, num_ocn_bc_restar

    !> @parblock
    !! INITIALIZE STDOUT, STDERR, STDLOG.
    !! @endparblock
    outunit = fms_mpp_stdout()
    errunit = fms_mpp_stderr()
    logunit = fms_mpp_stdlog()

    !> @parblock
    !! WRITE WALLDATE AND WALLTIME TO STDERR.
    !! @endparblock
    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Entering coupler_init at '//trim(walldate)//' '//trim(walltime)
    endif

    !> @parblock
    !! WRITE COUPLER VERSION
    !! @endparblock
    call fms_write_version_number('FULL_COUPLER_MOD', version)

    !> @parblock
    !! READ COUPLER_NML.
    !! @endparblock
    read (fms_mpp_input_nml_file, coupler_nml, iostat=io)
    ierr = fms_check_nml_error (io, 'coupler_nml')

    !> @parblock
    !! IF INPUT/COUPLER.RES EXISTS, READ CALENDAR_TYPE, DATE_INIT, AND DATE
    !! IF FILE DOES NOT EXIST, SET FORCE_DATE_FROM_NAMELIST = .TRUE.
    !! @endparblock
    if (fms2_io_file_exists('INPUT/coupler.res')) then
       call fms2_io_ascii_read('INPUT/coupler.res', restart_file)
       read(restart_file(1), *) calendar_type
       read(restart_file(2), *) date_ini
       read(restart_file(3), *) date
       deallocate(restart_file)
    else
      force_date_from_namelist = .true.
    endif

    !> @parblock
    !! IF FORCE_DATE_FROM_NAMELIST = .TRUE., SET DATE TO CURRENT_DATE AND CALENDAR_TYPE FROM NAMELIST.
    !! @endparblock
    if ( force_date_from_namelist ) then
      if ( sum(current_date) <= 0 ) then
        call fms_error_mesg ('program coupler',  &
             'no namelist value for base_date or current_date', FATAL)
      else
        date = current_date
      endif
      select case( fms_mpp_uppercase(trim(calendar)) )
      case( 'GREGORIAN' )
        calendar_type = GREGORIAN
      case( 'JULIAN' )
        calendar_type = JULIAN
      case( 'NOLEAP' )
        calendar_type = NOLEAP
      case( 'THIRTY_DAY' )
        calendar_type = THIRTY_DAY_MONTHS
      case( 'NO_CALENDAR' )
        calendar_type = NO_CALENDAR
      end selec
    endif
    call fms_time_manager_set_calendar_type (calendar_type, err_msg)
    if (err_msg /= '') then
      call fms_mpp_error(FATAL, 'ERROR in coupler_init: '//trim(err_msg))
    endif

    if (concurrent .AND. .NOT.(use_lag_fluxes .OR. concurrent_ice) ) call fms_mpp_error( WARNING, &
            & 'coupler_init: you have set concurrent=TRUE, &
            & use_lag_fluxes=FALSE, and concurrent_ice=FALSE &
            & in coupler_nml. When not using lag fluxes, components &
            & will synchronize at two points, and thus run serially.' )
    if (concurrent_ice .AND. .NOT.slow_ice_with_ocean ) call fms_mpp_error(WARNING, &
           &'coupler_init: concurrent_ice is true, but slow ice_with_ocean is &
           & false in coupler_nml.  These two flags should both be true to avoid &
           & effectively serializing the run.' )
    if (use_lag_fluxes .AND. concurrent_ice ) call fms_mpp_error(WARNING, &
           &'coupler_init: use_lag_fluxes and concurrent_ice are both true. &
           & These two coupling options are intended to be exclusive.' )

    !> @parblock
    !! INITIALIZE FMS ENSEMBLE_MANAGER
    !! FOR AN EMSEMBLE RUN, ENSEMBLE_MANAGER_INIT WILL RENAME ALL THE
    !! RESTART AND DIAGNOSTIC FILES TO CONTAIN THE NUMBER OF THE ENSEMBLE MEMBER.
    !! TO RESTART AN ENSEMBLE RUN, RESTART FILES MUST EXISTS FOR EACH ENSEMBLE MEMBER.
    !! @endparblock
    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Starting initializing ensemble_manager at '//trim(walldate)//' '//trim(walltime)
    endif
    call fms_ensemble_manager_init() ! init pelists for ensembles
    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Finished initializing ensemble_manager at '//trim(walldate)//' '//trim(walltime)
    endif
    ens_siz = fms_ensemble_manager_get_ensemble_size()
    ensemble_size = ens_siz(1)
    npes = ens_siz(2)

    !> @parblock
    !! CHECK PE ALLOCATION
    !! @endparblock
    if (concurrent) then
      !atmos_npes + ocean_npes must equal npes
      if (atmos_npes.EQ.0 ) atmos_npes = npes - ocean_npes
      if (ocean_npes.EQ.0 ) ocean_npes = npes - atmos_npes
      if (atmos_npes.EQ.0 .OR. ocean_npes.EQ.0 ) &
        call fms_mpp_error( FATAL, 'coupler_init: atmos_npes or ocean_npes must be specified for concurrent coupling.' )
      if (atmos_npes+ocean_npes.NE.npes ) &
        call fms_mpp_error( FATAL, 'coupler_init: atmos_npes+ocean_npes must equal npes for concurrent coupling.' )
    else !serial timestepping
      if ((atmos_npes.EQ.0) .and. (do_atmos .or. do_land .or. do_ice)) atmos_npes = npes
      if ((ocean_npes.EQ.0) .and. (do_ocean)) ocean_npes = npes
      if (max(atmos_npes,ocean_npes).EQ.npes) then !overlapping pelists
        ! do nothing
      else !disjoint pelists
        if (atmos_npes+ocean_npes.NE.npes ) call fms_mpp_error( FATAL,  &
             'coupler_init: atmos_npes+ocean_npes must equal npes for serial coupling on disjoint pelists.' )
      endif
    endif

    !> @parblock
    !! IF LAND_NPES IS 0, SET LAND_NPES = ATMOS_NPES.
    !! LAND_NPES SHOULD BE LESS THAN OR EQUAL TO ATMOS_NPES.
    !! @endparblock
    if (land_npes == 0 ) land_npes = atmos_npes
    if (land_npes > atmos_npes) call fms_mpp_error(FATAL, 'coupler_init: land_npes > atmos_npes')

    !> @parblock
    !! IF ICE_NPES IS 0, SET ICE_NPES = ATMOS_NPES.
    !! ICE_NPES SHOULD BE LESS THAN OR EQUAL TO ATMOS_NPES
    !! @endparblock
    if (ice_npes  == 0 ) ice_npes  = atmos_npes
    if (ice_npes  > atmos_npes) call fms_mpp_error(FATAL, 'coupler_init: ice_npes > atmos_npes')

    !> @parblock
    !! SETUP ATM%PELIST, OCEAN%PELIST, LAND%PELIST, AND ICE%FAST_PELIST
    !! IF CONCURRENT = .TRUE. ATMOS AND OCEAN GET DISTINCT PELIST.  ELSE
    !! ATMOS AND OCEAN HAVE OVERLAPPING PELISTS OR SAME PELIST.
    !! LAND%PELIST AND ICE%FAST_PELIST ARE SUBSETS OF ATM%PELIST
    !! @endparblock
    allocate( Atm%pelist  (atmos_npes) )
    allocate( Ocean%pelist(ocean_npes) )
    allocate( Land%pelist (land_npes) )
    allocate( Ice%fast_pelist(ice_npes) )
    call fms_ensemble_manager_ensemble_pelist_setup(concurrent, atmos_npes, ocean_npes, land_npes, ice_npes, &
                               Atm%pelist, Ocean%pelist, Land%pelist, Ice%fast_pelist)
    ensemble_id = fms_ensemble_manager_get_ensemble_id()
    if(allocated(ensemble_pelist)) call fms_mpp_error(FATAL, 'ensemble_pelist unexpectedly has already been allocated')
    allocate(ensemble_pelist(1:ensemble_size,1:npes))
    call fms_ensemble_manager_get_ensemble_pelist(ensemble_pelist)

    !> @parblock
    !! SET PE IDENTITY FOR ATM, OCEAN, AND LAND,
    !! I.E., ATM%PE = .TRUE. FOR PE IN ATM%PELIST
    !! @endparblock
    Atm%pe = ANY(Atm%pelist .EQ. fms_mpp_pe())
    Ocean%is_ocean_pe = ANY(Ocean%pelist .EQ. fms_mpp_pe())
    Land%pe = ANY(Land%pelist .EQ. fms_mpp_pe())

    !> @parblock
    !! SET ICE%PELISTS.  OUTLINE BELOW IS FOR WHEN DO_ATMOS IS TRUE:
    !! IF SLOW_ICE_WITH_OEAN = .FALSE. THEN ICE%SLOW_PELIST = ICE%FAST_PELIST.
    !! IF SLOW_ICE_WITH_OEAN = .TRUE., THEN ICE%SLOW_PELIST = OCEAN%PELIST, AND ICE%FAST_PELIST = ATM%PELIST.
    !! @endparblock
    Ice%shared_slow_fast_PEs = .not.slow_ice_with_ocean
    ! However, if using a data atmosphere and slow_ice_with_ocean then shared_slow_fast_PEs
    ! will be true. In this case, all procesors do the ocean, slow ice, and fast ice.
    if (slow_ice_with_ocean.and.(.not.do_atmos)) Ice%shared_slow_fast_PEs = .true.
    ! This is where different settings would be applied if the fast and slow ice occurred on different PEs.
    if (do_atmos) then
      if (Ice%shared_slow_fast_PEs) then
        ! Fast and slow ice processes occur on the same PEs.
        allocate( Ice%pelist  (ice_npes) )
        Ice%pelist(:) = Ice%fast_pelist(:)
        allocate( Ice%slow_pelist(ice_npes) )
        Ice%slow_pelist(:) = Ice%fast_pelist(:)
        if(concurrent) then
          if(.not.allocated(slow_ice_ocean_pelist)) then
            allocate(slow_ice_ocean_pelist(ocean_npes+ice_npes))
          else
            call fms_mpp_error(FATAL, 'allocation of slow_ice_ocean_pelist unexpectedly has already been allocated')
          end if
          slow_ice_ocean_pelist(1:ice_npes) = Ice%slow_pelist(:)
          slow_ice_ocean_pelist(ice_npes+1:ice_npes+ocean_npes) = Ocean%pelist(:)
        else
          if(ice_npes .GE. ocean_npes) then
             allocate(slow_ice_ocean_pelist(ice_npes))
             slow_ice_ocean_pelist(:) = Ice%slow_pelist(:)
          else
             allocate(slow_ice_ocean_pelist(ocean_npes))
             slow_ice_ocean_pelist(:) = Ocean%pelist(:)
          endif
        endif
      else
        ! Fast ice processes occur a subset of the atmospheric PEs, while
        ! slow ice processes occur on the ocean PEs.
        allocate( Ice%slow_pelist(ocean_npes) )
        Ice%slow_pelist(:) = Ocean%pelist(:)
        allocate( Ice%pelist  (ice_npes+ocean_npes) )
        ! Set Ice%pelist() to be the union of Ice%fast_pelist and Ice%slow_pelist.
        Ice%pelist(1:ice_npes) = Ice%fast_pelist(:)
        Ice%pelist(ice_npes+1:ice_npes+ocean_npes) = Ocean%pelist(:)
        allocate(slow_ice_ocean_pelist(ocean_npes))
        slow_ice_ocean_pelist(:) = Ocean%pelist(:)
      endif
    elseif (.not.do_atmos) then
      ! In the no atmos cases, shared_slow_fast_PEs is not enough to distinguish
      ! the slow and fast ice procesor layout; slow_ice_with_ocean should be used instead.
      if (slow_ice_with_ocean) then
        ! data atmos, using combined ice-ocean driver
        ! Both fast ice and slow ice processes occur on the same PEs,
        ! since the Atmos and Ocean PEs are shared
        allocate( Ice%slow_pelist(ocean_npes) )
        Ice%slow_pelist(:) = Ocean%pelist(:)
        allocate( Ice%pelist  (ice_npes) )
        Ice%pelist(1:ice_npes) = Ice%fast_pelist(:)
        allocate(slow_ice_ocean_pelist(ocean_npes))
        slow_ice_ocean_pelist(:) = Ocean%pelist(:)
      else
        ! data atmos, not using combined ice-ocean driver
        allocate( Ice%pelist  (ice_npes) )
        Ice%pelist(:) = Ice%fast_pelist(:)
        allocate( Ice%slow_pelist(ice_npes) )
        Ice%slow_pelist(:) = Ice%fast_pelist(:)
        if(ice_npes .GE. ocean_npes) then
           allocate(slow_ice_ocean_pelist(ice_npes))
           slow_ice_ocean_pelist(:) = Ice%slow_pelist(:)
        else
           allocate(slow_ice_ocean_pelist(ocean_npes))
           slow_ice_ocean_pelist(:) = Ocean%pelist(:)
        endif
      endif
    endif
    Ice%fast_ice_pe = ANY(Ice%fast_pelist(:) .EQ. fms_mpp_pe())
    Ice%slow_ice_pe = ANY(Ice%slow_pelist(:) .EQ. fms_mpp_pe())
    Ice%pe = Ice%fast_ice_pe .OR. Ice%slow_ice_pe
    call fms_mpp_declare_pelist(slow_ice_ocean_pelist)

    !> @parblock
    !! SET OMP FOR WHEN DO_CONCURRENT_RADIATION = .TRUE.
    !! @endparblock
    !--- dynamic threading turned off when affinity placement is in use
!$  call omp_set_dynamic(.FALSE.)
    !--- nested OpenMP enabled for OpenMP concurrent components
!$  call omp_set_max_active_levels(3)

    if (Atm%pe) then
      call fms_mpp_set_current_pelist( Atm%pelist )
!$    if (.not.do_concurrent_radiation) radiation_nthreads=atmos_nthreads
!$    if (do_concurrent_radiation) conc_nthreads=2
      !--- setting affinity
      if (do_concurrent_radiation) then
!$      call fms_affinity_set('ATMOS', use_hyper_thread, atmos_nthreads + radiation_nthreads)
!$      call omp_set_num_threads(atmos_nthreads+radiation_nthreads)
      else
!$      call fms_affinity_set('ATMOS', use_hyper_thread, atmos_nthreads)
!$      call omp_set_num_threads(atmos_nthreads)
      endif
    endif

    !> @parblock
    !! INITIALIZE CLOCKS FOR PROFILING.
    !! @endparblock
    ! The pelists need to be set before initializing the clocks
    call coupler_set_clock_ids(coupler_clocks, Atm, Land, Ice, Ocean, ensemble_pelist, &
                               slow_ice_ocean_pelist, ensemble_id)

    !> @parblock
    !! WRITE PELISTS TO LOG.
    !! @endparblock
    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      write( text,'(a,2i6,a,i2.2)' ) &
        'Atmos PE range: ', Atm%pelist(1), Atm%pelist(atmos_npes), ' ens_', ensemble_id
      call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
      if (ocean_npes .gt. 0) then
        write( text,'(a,2i6,a,i2.2)' )&
          'Ocean PE range: ', Ocean%pelist(1),Ocean%pelist(ocean_npes),' ens_', ensemble_id
        call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
      else
        write( text,'(a,i2.2)' ) &
          'Ocean PE range is not set (do_ocean=.false. and concurrent=.false.) for ens_', ensemble_id
        call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
      endif
      write( text,'(a,2i6,a,i2.2)' ) &
        'Land PE range: ', Land%pelist(1)  , Land%pelist(land_npes),' ens_', ensemble_id
      call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
      if (.not.concurrent_ice) then
        write( text,'(a,2i6,a,i2.2)' ) &
          'Ice PE range: ', Ice%pelist(1), Ice%pelist(ice_npes), ' ens_', ensemble_id
        call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
      elseif (concurrent_ice) then
        if (do_atmos) then
          write( text,'(a,2i6,a,i2.2)' ) &
            'Ice PE range: ', Ice%pelist(1), Ice%pelist(ice_npes+ocean_npes), ' ens_', ensemble_id
          call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
        elseif ((.not.do_atmos)) then
          write( text,'(a,2i6,a,i2.2)' ) &
            'Ice PE range: ', Ice%pelist(1), Ice%pelist(ice_npes), ' ens_', ensemble_id
          call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
        endif
        call fms_mpp_error( NOTE, 'coupler_init: Running with CONCURRENT ICE coupling.' )
        write( text,'(a,2i6,a,i2.2)' ) &
          'slow Ice PE range: ', Ice%slow_pelist(1), Ice%slow_pelist(ocean_npes), ' ens_', ensemble_id
        call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
        write( text,'(a,2i6,a,i2.2)' ) &
          'fast Ice PE range: ', Ice%fast_pelist(1), Ice%fast_pelist(ice_npes), ' ens_', ensemble_id
        call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
      endif
      if (concurrent) then
        call fms_mpp_error( NOTE, 'coupler_init: Running with CONCURRENT coupling.' )

        write( logunit,'(a)' )'Using concurrent coupling...'
        write( logunit,'(a,4i6)' ) &
              'atmos_pe_start, atmos_pe_end, ocean_pe_start, ocean_pe_end=', &
              Atm%pelist(1)  , Atm%pelist(atmos_npes), Ocean%pelist(1), Ocean%pelist(ocean_npes)
      else
        call fms_mpp_error( NOTE, 'coupler_init: Running with SERIAL coupling.' )
      endif
      if (use_lag_fluxes) then
        call fms_mpp_error( NOTE, 'coupler_init: Sending LAG fluxes to ocean.' )
      else
        call fms_mpp_error( NOTE, 'coupler_init: Sending most recent fluxes to ocean.' )
      endif
      if (concurrent_ice) call fms_mpp_error( NOTE, &
        'coupler_init: using lagged slow-ice coupling mode.')
      if (combined_ice_and_ocean) call fms_mpp_error( NOTE, &
        'coupler_init: advancing the ocean and slow-ice in a single call.')
      if (combined_ice_and_ocean .and. .not.concurrent_ice) call fms_mpp_error( FATAL, &
        'coupler_init: concurrent_ice must be true if combined_ice_and_ocean is true.')
      if (combined_ice_and_ocean .and. .not.slow_ice_with_ocean) call fms_mpp_error( FATAL, &
        'coupler_init: slow_ice_with_ocean must be true if combined_ice_and_ocean is true.')
    endif

    !> @parblock
    !! WRITE NAMELIST TO LOG.
    !! @endparblock
    if (fms_mpp_pe() == fms_mpp_root_pe() )write( logunit, nml=coupler_nml )

    !> @parblock
    !! WRITE MODEL INITIAL DATE TO LOGFILE.
    !! @endparblock
    if ( fms_mpp_pe().EQ.fms_mpp_root_pe() ) &
      write( logunit, 16 )date(1),trim(fms_time_manager_month_name(date(2))),date(3:6)
16  format ('  current date used = ',i4,1x,a,2i3,2(':',i2.2),' gmt')

!jwd Fork here is somewhat dangerous. It relies on "no side effects" from
!    diag_manager_init. diag_manager_init or this section should be
!    re-architected to guarantee this or remove this assumption.
!    For instance, what follows assumes that get_base_date has the same
!    time for both Atm and Ocean pes. While this should be the case, the
!    possible error condition needs to be checked

    !> @parblock
    !! INITIALIZE DIAG_MANAGER AND READ DIAG_TABLE.
    !! @endparblock
    diag_model_subset=DIAG_ALL
    if (Atm%pe) then
      call fms_mpp_set_current_pelist(Atm%pelist)
      if (atmos_npes /= npes) diag_model_subset = DIAG_OTHER  ! change diag_model_subset from DIAG_ALL
    elseif (Ocean%is_ocean_pe) then  ! Error check above for disjoint pelists should catch any problem
      call fms_mpp_set_current_pelist(Ocean%pelist)
      ! The FMS diag manager has a convention that segregates files with "ocean"
      ! in their names from the other files to handle long diag tables.  This
      ! does not work if the ice is on the ocean PEs.
      if ((ocean_npes /= npes) .and. .not.slow_ice_with_ocean) &
        diag_model_subset = DIAG_OCEAN  ! change diag_model_subset from DIAG_ALL
    endif
    if ( fms_mpp_pe() == fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Starting to initialize diag_manager at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
    ! initialize diag_manager for processor subset outpu
    call fms_diag_init(DIAG_MODEL_SUBSET=diag_model_subset, TIME_INIT=date)
    call fms_memutils_print_memuse_stats( 'diag_manager_init' )
    if ( fms_mpp_pe() == fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Finished initializing diag_manager at '&
                       //trim(walldate)//' '//trim(walltime)
    endif

    call fms_mpp_set_current_pelist()

    !> @parblock
    !! OVERRIDE DATE_INIT WITH BASE DATE FROM DIAG_MANAGER IF BASE DATE EXISTS IN DIAG_TABLE.
    !! @endparblock
    call fms_diag_get_base_date ( date_init(1), date_init(2), date_init(3), &
         date_init(4), date_init(5), date_init(6)  )

    if ( date_init(1) == 0 ) date_init = date

    !> @parblock
    !! SET TIME_INIT, TIME, AND TIME_START FROM DATE_INIT.
    !! @endparblock
    Time_init = fms_time_manager_set_date (date_init(1), date_init(2), date_init(3), &
         date_init(4), date_init(5), date_init(6))

    !> @parblock
    !! SET TIME FROM DATE.
    !! @endparblock
    Time  = fms_time_manager_set_date (date(1), date(2), date(3), date(4), date(5), date(6))

    !> @parblock
    !! SET TIME_START = TIME (from date)
    Time_start = Time

    !> @parblock
    !! COMPUTE TIME_END  FROM MONTHS, DAYS, HOURS, MINUTES, AND SECONDS.
    !! @endparblock
    Time_end = Time
    do m=1,months
       Time_end = Time_end + fms_time_manager_set_time(0,fms_time_manager_days_in_month(Time_end))
    enddo
    Time_end   = Time_end + fms_time_manager_set_time(hours*3600+minutes*60+seconds, days)

    !> @parblock
    !! CALL FMS_DIAG_SET_TIME_END WITH TIME_END.
    !! @endparblock
    !Need to pass Time_end into diag_manager for multiple thread case.
    call fms_diag_set_time_end(Time_end)

    !> @parblock
    !! GET RUN_LENGTH = TIME_END - TIME.
    !! @endparblock
    Run_length = Time_end - Time

    !> @parblock
    !! IF INPUT/COUPLER.INTERMEDIATE.RES EXISTS, READ DATE_RESTART FROM THIS FILE.
    !! ELSE SET DATE_RESTART = DATE.
    !! @endparblock
    if (fms2_io_file_exists('INPUT/coupler.intermediate.res')) then
       call fms2_io_ascii_read('INPUT/coupler.intermediate.res', restart_file)
       read(restart_file(1), *) date_restar
       deallocate(restart_file)
    else
       date_restart = date
    endif

    !> @parblock
    !! SET TIME_RESTART.
    !! @endparblock
    Time_restart_current = Time
    if (ALL(restart_interval ==0)) then
       Time_restart = fms_time_manager_increment_date(Time_end, 0, 0, 10, 0, 0, 0)   ! no intermediate restar
    else
       Time_restart = fms_time_manager_set_date(date_restart(1), date_restart(2), date_restart(3),  &
                               date_restart(4), date_restart(5), date_restart(6) )
       Time_restart = fms_time_manager_increment_date(Time_restart, restart_interval(1), restart_interval(2), &
            restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )
       if (Time_restart <= Time) call fms_mpp_error(FATAL, &
            '==>Error from program coupler: The first intermediate restart time is no larger than the start time')
    endif

    !> @parblock
    !! WRITE STARTING AND ENDING TIME TO LOG.
    !! @endparblock
    if ( fms_mpp_pe().EQ.fms_mpp_root_pe() ) &
      open(newunit = time_stamp_unit, file='time_stamp.out', status='replace', form='formatted')
    month = fms_time_manager_month_name(date(2))
    if ( fms_mpp_pe().EQ.fms_mpp_root_pe() ) write (time_stamp_unit,20) date, month(1:3)
    call fms_time_manager_get_date (Time_end, date(1), date(2), date(3), date(4), date(5), date(6))
    month = fms_time_manager_month_name(date(2))
    if ( fms_mpp_pe().EQ.fms_mpp_root_pe() ) write (time_stamp_unit,20) date, month(1:3)
    if ( fms_mpp_pe().EQ.fms_mpp_root_pe() ) close(time_stamp_unit)
20  format (i6,5i4,2x,a3)

    !> @parblock
    !! SET TIME STEPS AND TOTAL NUMBER OF LOOP ITERATIONS
    !! CHECK FOR CONSISTENCY IN TIME STEPS AND RUN LENGTH.
    !! @endparblock
    Time_step_cpld  = fms_time_manager_set_time (dt_cpld ,0)
    Time_step_atmos = fms_time_manager_set_time (dt_atmos,0)
    num_cpld_calls  = Run_length      / Time_step_cpld
    num_atmos_calls = Time_step_cpld  / Time_step_atmos
    if ( Time_init > Time ) &
      call fms_error_mesg ('program coupler', 'initial time is greater than current time', FATAL)
    if ( num_cpld_calls * Time_step_cpld  /= Run_length )  &
      call fms_error_mesg ('program coupler', 'run length must be multiple of coupled time step', FATAL)
    if ( num_atmos_calls * Time_step_atmos /= Time_step_cpld )  &
      call fms_error_mesg ('program coupler', 'cpld time step is not a multiple of the atmos time step', FATAL)
!
!       Initialize the tracer manager. This needs to be done on all PEs,
!       before the individual models are initialized.
!
    !> @parblock
    !! INITIALIZE TRACER MANAGER AND GAS EXCHANGE FLUXES.
    !! @endparblock
    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Starting to initialize tracer_manager at '//trim(walldate)//' '//trim(walltime)
    endif
    call fms_tracer_manager_init()
    ! Initialize the gas-exchange fluxes so this information can be made
    ! available to the individual components.
    call gas_exchange_init(gas_fields_atm, gas_fields_ocn, gas_fluxes)
    call fms_coupler_types_init()
    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Finished initializing tracer_manager at '//trim(walldate)//' '//trim(walltime)
    endif
    ! Initialize atm/land exchange (not for tracers)
    call gex_init()

    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Beginning to initialize component models at '//trim(walldate)//' '//trim(walltime)
    endif

    !> @parblock
    !! INITIALIZE ATM MODEL ON ATM%PES INCLUDING DATA_OVERRIDE_INIT FOR ATM.
    !! @endparblock
    if (Atm%pe) then
        call fms_mpp_set_current_pelist(Atm%pelist)
        !---- atmosphere ----
        if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Starting to initialize atmospheric model at '//trim(walldate)//' '//trim(walltime)
        endif

        call fms_mpp_clock_begin(coupler_clocks%atmos_model_init)
        call atmos_model_init( Atm, Time_init, Time, Time_step_atmos, do_concurrent_radiation)
        call fms_mpp_clock_end(coupler_clocks%atmos_model_init)

        if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Finished initializing atmospheric model at '//trim(walldate)//' '//trim(walltime)
        endif
        call fms_memutils_print_memuse_stats( 'atmos_model_init' )
        call fms_data_override_init(Atm_domain_in = Atm%domain)
    endif

    !> @parblock
    !! INITIALIZE LAND MODEL ON LAND%PES INCLUDING DATA_OVERRIDE_INIT FOR LAND.
    !! @endparblock
    if (Land%pe) then
      call fms_mpp_set_current_pelist(Land%pelist)
      if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
        write(errunit,*) 'Starting to initialize land model at '//trim(walldate)//' '//trim(walltime)
      endif

      call fms_mpp_clock_begin(coupler_clocks%land_model_init)
      call land_model_init( Atmos_land_boundary, Land, Time_init, Time, Time_step_atmos, Time_step_cpld )
      call fms_mpp_clock_end(coupler_clocks%land_model_init)

      if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
        write(errunit,*) 'Finished initializing land model at '//trim(walldate)//' '//trim(walltime)
      endif
      call fms_memutils_print_memuse_stats( 'land_model_init' )
      call fms_data_override_init(Land_domain_in = Land%domain)
#ifndef _USE_LEGACY_LAND_
      call fms_data_override_init(Land_domainUG_in = Land%ug_domain)
#endif
    endif

    !> @parblock
    !! INITIALIZE ICE MODEL ON BOTH FAST AND SLOW ICE%PES INCLUDING DATA_OVERRIDE_INIT FOR ICE.
    !! @endparblock
    if (Ice%pe) then
      if (Ice%fast_ice_pe) then
        call fms_mpp_set_current_pelist(Ice%fast_pelist)
      elseif (Ice%slow_ice_pe) then
        call fms_mpp_set_current_pelist(Ice%slow_pelist)
      else
        call fms_mpp_error(FATAL, "All Ice%pes must be a part of Ice%fast_ice_pe or Ice%slow_ice_pe")
      endif
      if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
        write(errunit,*) 'Starting to initialize ice model at '//trim(walldate)//' '//trim(walltime)
      endif
      call fms_mpp_clock_begin(coupler_clocks%ice_model_init)
      call ice_model_init(Ice, Time_init, Time, Time_step_atmos, Time_step_cpld, Verona_coupler=.false., &
                          concurrent_ice=concurrent_ice, gas_fluxes=gas_fluxes, gas_fields_ocn=gas_fields_ocn )
      call fms_mpp_clock_end(coupler_clocks%ice_model_init)
      ! This must be called using the union of the ice PE_lists.
      call fms_mpp_set_current_pelist(Ice%pelist)
      call share_ice_domains(Ice)

      if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
        write(errunit,*) 'Finished initializing ice model at '//trim(walldate)//' '//trim(walltime)
      endif
      call fms_memutils_print_memuse_stats( 'ice_model_init' )
      if (Ice%fast_ice_pe) then
        call fms_mpp_set_current_pelist(Ice%fast_pelist)
        call fms_data_override_init(Ice_domain_in = Ice%domain)
      endif
    endif

    !> @parblock
    !! INITIALIZE OCEAN MODEL ON OCEAN%PES INCLUDING DATA_OVERRIDE_INIT FOR OCEAN AND OPENMP THREADS
    !! @endparblock
    if (Ocean%is_ocean_pe) then
      call fms_mpp_set_current_pelist(Ocean%pelist)
      if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
        write(errunit,*) 'Starting to initialize ocean model at '//trim(walldate)//' '//trim(walltime)
      endif

      call fms_mpp_clock_begin(coupler_clocks%ocean_model_init)
      call ocean_model_init( Ocean, Ocean_state, Time_init, Time, &
                             gas_fields_ocn=gas_fields_ocn  )
      call fms_mpp_clock_end(coupler_clocks%ocean_model_init)

      if (concurrent) then
        call fms_mpp_set_current_pelist( Ocean%pelist )
!$      call fms_affinity_set('OCEAN', use_hyper_thread, ocean_nthreads)
!$      call omp_set_num_threads(ocean_nthreads)
      else
        ocean_nthreads = atmos_nthreads
        !--- omp_num_threads has already been set by the Atmos-pes, but set again to ensure
!$      call omp_set_num_threads(ocean_nthreads)
      endif

      if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
        write(errunit,*) 'Finished initializing ocean model at '//trim(walldate)//' '//trim(walltime)
      endif
      call fms_memutils_print_memuse_stats( 'ocean_model_init' )
      if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
        write(errunit,*) 'Starting to initialize data_override at '//trim(walldate)//' '//trim(walltime)
      endif
      call fms_data_override_init(Ocean_domain_in = Ocean%domain )
      if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
        write(errunit,*) 'Finished initializing data_override at '//trim(walldate)//' '//trim(walltime)
      endif

      if (combined_ice_and_ocean) &
        call ice_ocean_driver_init(ice_ocean_driver_CS, Time_init, Time)
    endif

    !> @parblock
    !! CALL MPP_DOMAINS_BROADCAST_DOMAIN FOR ICE AND OCEAN TO SHARE DOMAIN INFORMATION.
    !! @endparblock
    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Finished initializing component models at '//trim(walldate)//' '//trim(walltime)
    endif
    call fms_mpp_set_current_pelist(ensemble_pelist(ensemble_id,:))
    call fms_mpp_domains_broadcast_domain(Ice%domain)
    call fms_mpp_domains_broadcast_domain(Ice%slow_domain_NH)
    call fms_mpp_domains_broadcast_domain(Ocean%domain)

    !> @parblock
    !! INITIALIZE FLUX EXCHANGE.
    !! @endparblock
    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Starting to initialize flux_exchange at '//trim(walldate)//' '//trim(walltime)
    endif

    call fms_mpp_clock_begin(coupler_clocks%flux_exchange_init)
    if(do_flux) call flux_exchange_init ( Time, Atm, Land, Ice, Ocean, Ocean_state,&
         atmos_ice_boundary, land_ice_atmos_boundary, land_ice_boundary, ice_ocean_boundary, ocean_ice_boundary, &
         do_ocean, slow_ice_ocean_pelist, dt_atmos=dt_atmos, dt_cpld=dt_cpld)
    call fms_mpp_set_current_pelist(ensemble_pelist(ensemble_id,:))
    call fms_mpp_clock_end(coupler_clocks%flux_exchange_init)

    call fms_mpp_set_current_pelist()
    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Finished initializing flux_exchange at '//trim(walldate)//' '//trim(walltime)
    endif

    !> @parblock
    !! SET TIME_ATMOS = TIME_OCEAN = TIME.
    !! @endparblock
    Time_atmos = Time
    Time_ocean = Time

    !> @parblock
    !! READ ICE FROM RESTART.
    !! @endparblock
    if ( Ice%slow_ice_pe ) then
      call fms_mpp_set_current_pelist(Ice%slow_pelist)
      call fms_coupler_type_register_restarts(Ice%ocean_fluxes, Ice_bc_restart, &
             num_ice_bc_restart, Ice%slow_domain_NH, to_read=.true., ocean_restart=.false., directory="INPUT/")

      ! Restore the fields from the restart files
      do l = 1, num_ice_bc_restar
         if(fms2_io_check_if_open(Ice_bc_restart(l))) call fms2_io_read_restart(Ice_bc_restart(l))
      enddo

      ! Check whether the restarts were read successfully.
      call fms_coupler_type_restore_state(Ice%ocean_fluxes, use_fms2_io=.true., test_by_field=.true.)

      do l = 1, num_ice_bc_restar
        if(fms2_io_check_if_open(Ice_bc_restart(l))) call fms2_io_close_file(Ice_bc_restart(l))
      enddo
    endif

    !> @parblock
    !! READ OCEAN FROM RESTART
    !! @endparblock
    if ( Ocean%is_ocean_pe ) then
      call fms_mpp_set_current_pelist(Ocean%pelist)

      call fms_coupler_type_register_restarts(Ocean%fields, Ocn_bc_restart, &
               num_ocn_bc_restart, Ocean%domain, to_read=.true., ocean_restart=.true., directory="INPUT/")

      ! Restore the fields from the restart files
      do l = 1, num_ocn_bc_restar
         if(fms2_io_check_if_open(Ocn_bc_restart(l))) call fms2_io_read_restart(Ocn_bc_restart(l))
      enddo

      ! Check whether the restarts were read successfully.
      call fms_coupler_type_restore_state(Ocean%fields, use_fms2_io=.true., test_by_field=.true.)

      do l = 1, num_ocn_bc_restar
         if(fms2_io_check_if_open(Ocn_bc_restart(l))) call fms2_io_close_file(Ocn_bc_restart(l))
      enddo
    endif

    call fms_mpp_set_current_pelist()

    !> @parblock
    !! MISCELLANEOUS INCLUDING CALLING DIAG_GRID_END TO FREE UP MEMORY USED
    !! DURING REGIONAL OUTPUT SETUP.
    !! @endparblock
    if ( fms_mpp_pe().EQ.fms_mpp_root_pe() ) then
       open(newunit = ascii_unit, file='RESTART/file', status='replace', form='formatted')
       close(ascii_unit,status="delete")
    endif
    call fms_diag_grid_end()

    !> @parblock
    !! INITIALIZE COUPLER_COMPONENT_OBJ
    !! @endparblock
    call coupler_components_obj%initialize_coupler_components_obj(Atm, Land, Ice, Ocean, Land_ice_atmos_boundary,&
        Atmos_land_boundary, Atmos_ice_boundary, Land_ice_boundary, Ice_ocean_boundary, Ocean_ice_boundary)

    !> @parblock
    !! INITIALIZE COUPLER CHECKSUM OBJECT.
    !! @endparblock
    call coupler_chksum_obj%initialize_coupler_chksum_obj(coupler_components_obj)

    !> @parblock
    !! IF DO_ENDPOINT_CHKSUM IS TRUE, COMPUTE CHECKSUM.
    !! @endparblock
    if ( do_endpoint_chksum ) then
      call coupler_chksum_obj%get_atmos_ice_land_ocean_chksums('coupler_init+', 0)
      if (Ice%slow_ice_PE) then
        call fms_mpp_set_current_pelist(Ice%slow_pelist)
        call coupler_chksum_obj%get_slow_ice_chksums('coupler_init+', 0)
      end if
    end if

    call fms_mpp_set_current_pelist()

    !> @parblock
    !! LOG.
    !! @endparblock
    call fms_memutils_print_memuse_stats('coupler_init')
    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Exiting coupler_init at '//trim(walldate)//' '//trim(walltime)
    endif

  end subroutine coupler_ini

  !> @parblock
  !! Subroutine initialize_coupler_components_obj is a typed-bound procedure to the coupler_components_type.
  !! This subroutine associates each pointer member of the coupler_components_type object to the corresponding
  !! model component derived type.
  !! @endparblock
  subroutine initialize_coupler_components_obj(this, Atm, Land, Ice, Ocean, Land_ice_atmos_boundary, &
      Atmos_land_boundary, Atmos_ice_boundary, Land_ice_boundary, Ice_ocean_boundary, Ocean_ice_boundary)

    implicit none
    class(coupler_components_type), intent(inout) :: this
      !< is the reference to self (coupler_components_type)
    type(atmos_data_type), target, intent(in) :: Atm
      !< is the Atm derived type containing atmospheric model data and metadata
    type(land_data_type), target, intent(in) :: Land
      !< is the Land derived type containing land model data and metadata
    type(ice_data_type), target, intent(in) :: Ice
      !< is the Ice derived type containing ice model data and metadata
    type(ocean_public_type), target, intent(in) :: Ocean
      !< is the Ocean derived type containing ocean model data and metadata
    type(land_ice_atmos_boundary_type), target, intent(in) :: Land_ice_atmos_boundary
      !< is the Land_ice_atmos_boundary derived type containing data and metadata for the land-ice-atmosphere boundary
    type(atmos_land_boundary_type), target, intent(in) :: Atmos_land_boundary
      !< is the Atmos_land_boundary derived type containing data and metadata for the atmosphere-land boundary
    type(atmos_ice_boundary_type), target, intent(in) :: Atmos_ice_boundary
      !< is the Atmos_ice_boundary derived type containing data and metadata for the atmosphere-ice boundary
    type(land_ice_boundary_type), target, intent(in) :: Land_ice_boundary
      !< is the Land_ice_boundary derived type containing data and metadata for the land-ice boundary
    type(ice_ocean_boundary_type), target, intent(in) :: Ice_ocean_boundary
      !< is the Ice_ocean_boundary derived type containing data and metadata for the ice-ocean boundary
    type(ocean_ice_boundary_type), target, intent(in) :: Ocean_ice_boundary
      !< is the Ocean_ice_boundary derived type containing data and metadata for the ocean-ice boundary

    !> @parblock
    !! POINTER ASSOCIATION
    !! @endparblock
    this%Atm => Atm
    this%Land => Land
    this%Ice => Ice
    this%Ocean => Ocean
    this%Land_ice_atmos_boundary => Land_ice_atmos_boundary
    this%Atmos_land_boundary => Atmos_land_boundary
    this%Atmos_ice_boundary => Atmos_ice_boundary
    this%Land_ice_boundary => Land_ice_boundary
    this%Ice_ocean_boundary => Ice_ocean_boundary
    this%Ocean_ice_boundary => Ocean_ice_boundary

  end subroutine initialize_coupler_components_obj

  !> @parblock
  !! Subroutine get_component is a type-bound procedure to coupler_commponents_type and
  !! retrieves the requested component.  For example,
  !! coupler_components_obj%get_component(Atm) retrieves coupler_components_obj%Atm,
  !! which is a pointer to the atmospheric component derived type and a private member
  !! of coupler_components_obj.
  !! @endparblock
  subroutine get_component(this, retrieve_component )

    implicit none
    class(coupler_components_type), intent(in) :: this
      !< is the reference to self (coupler_components_type object)
    class(*), intent(out) :: retrieve_componen
      !< is the requested component to be retrieve.
      !! retrieve_component can be of type atmos_data_type, land_data_type, ice_data_type,
      !! ocean_public_type, land_ice_atmos_boundary_type, atmos_land_boundary_type,
      !! atmos_ice_boundary_type, land_ice_boundary_type, ice_ocean_boundary_type,
      !! ocean_ice_boundary_type

    select type(retrieve_component)
    type is(atmos_data_type) ; retrieve_component = this%Atm
    type is(land_data_type)  ; retrieve_component = this%Land
    type is(ice_data_type)   ; retrieve_component = this%Ice
    type is(ocean_public_type) ; retrieve_component = this%Ocean
    type is(land_ice_atmos_boundary_type) ; retrieve_component = this%Land_ice_atmos_boundary
    type is(atmos_land_boundary_type) ; retrieve_component = this%Atmos_land_boundary
    type is(atmos_ice_boundary_type)  ; retrieve_component = this%Atmos_ice_boundary
    type is(land_ice_boundary_type)   ; retrieve_component = this%Land_ice_boundary
    type is(ice_ocean_boundary_type)  ; retrieve_component = this%Ice_ocean_boundary
    type is(ocean_ice_boundary_type)  ; retrieve_component = this%Ocean_ice_boundary
    class defaul
      call fms_mpp_error(FATAL, "failure retrieving component in coupler_components_type object, &
                         cannot recognize the type of requested component")
    end selec

  end subroutine get_componen

  !> @parblock
  !! Subroutine initialize_coupler_chksum_obj is a type-bound procedure to coupler_chksum_type and
  !! associates coupler_chksum_obj%components => components_obj.   After this call, the chksum objec
  !! can access all component model derived types through the pointer.
  !! @endparblock
  subroutine initialize_coupler_chksum_obj(this, components_obj)

    implicit none
    class(coupler_chksum_type), intent(inout) :: this
      !< The coupler_chksum_type object being initialized
    type(coupler_components_type), intent(in), target :: components_obj
      !< The components object whose address will be stored

    this%components => components_obj

  end subroutine initialize_coupler_chksum_obj

  !> @parblock
  !! Subroutine get_components_obj is a type-bound procedure to coupler_chksum_type and
  !! retrieves the coupler_components_type object stored inside a coupler_chksum_type object.
  !! For example, coupler_chksum_obj%get_components_obj(components_obj) retrieves the components_obj,
  !! which is a pointer to the coupler_components_type object and a private member of coupler_chksum_obj.
  !! @endparblock
  subroutine get_components_obj(this, components_obj)

    implicit none

    class(coupler_chksum_type), intent(in) :: this
      !< is a reference to self (coupler_chksum_type)
    type(coupler_components_type), intent(out) :: components_obj
      !< is the coupler_components_type to be returned

    components_obj = this%components

  end subroutine get_components_obj

  !> @parblock
  !! Subroutine coupler_end finalizes all component models (such as deallocating arrays),
  !! writes restart files, and calls fms_diag_end to flush and close all diagnostic output files.
  !! Checksums are computed when do_chksum or do_endpoint_chksum is .true.
  !! @endparblock
  subroutine coupler_end(Atm, Land, Ice, Ocean, Ocean_state, Land_ice_atmos_boundary, Atmos_ice_boundary,&
                         Atmos_land_boundary, Ice_ocean_boundary, Ocean_ice_boundary, Ocn_bc_restart,    &
                         Ice_bc_restart, current_timestep, Time_current, Time_start, Time_end, Time_restart_current,&
                         coupler_chksum_obj, coupler_clocks)

    implicit none

    type(atmos_data_type), intent(inout) :: Atm
      !< is the atmospheric derived type
    type(land_data_type), intent(inout) :: Land
      !< is the land derived type
    type(ice_data_type), intent(inout) :: Ice
      !< is the ice derived type
    type(ocean_public_type), intent(inout) :: Ocean
      !< is the ocean derived type
    type(ocean_state_type), pointer, intent(inout) :: Ocean_state
      !< is the ocean state derived type
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_ice_atmos_boundary
      !< is the land-ice-atmosphere boundary derived type
    type(atmos_ice_boundary_type), intent(inout) :: Atmos_ice_boundary
      !< is the atmosphere-ice boundary derived type
    type(atmos_land_boundary_type), intent(inout) :: Atmos_land_boundary
      !< is the atmosphere-land boundary derived type
    type(ice_ocean_boundary_type), intent(inout) :: Ice_ocean_boundary
      !< is the ice-ocean boundary derived type
    type(ocean_ice_boundary_type), intent(inout) :: Ocean_ice_boundary
      !< is the ocean-ice boundary derived type
    type(FmsNetcdfDomainFile_t), dimension(:), pointer, intent(inout) :: Ocn_bc_restar
      !< is required to write restart files
    type(FmsNetcdfDomainFile_t), dimension(:), pointer, intent(inout) :: Ice_bc_restar
      !< is required to write restart files
    integer, intent(in) :: current_timestep
      !< is the current timestep (nc)
    type(coupler_clock_type), intent(in)  :: coupler_clocks
      !< are the coupler clocks
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj
      !< is required for chksum computations

    type(FmsTime_type), intent(in) :: Time_curren
      !< is the current timestep
    type(FmsTime_type), intent(in) :: Time_star
      !< is the model starting time
    type(FmsTime_type), intent(in) :: Time_end
      !< is the model ending time
    type(FmsTime_type), intent(in) :: Time_restart_curren
      !< is the time corresponding to last restart time

    call fms_mpp_clock_begin(coupler_clocks%termination)

    !> @parblock
    !! IF DO_CHKSUM AND/OR DO_ENDPOINT_CHKSUM IS TRUE, COMPUTE CHECKSUMS
    !! @endparblock
    if (do_chksum) call coupler_chksum_obj%get_coupler_chksums('coupler_end-', current_timestep)
    if ( do_endpoint_chksum ) then
      call coupler_chksum_obj%get_atmos_ice_land_ocean_chksums('coupler_end', 0)
      if (Ice%slow_ice_PE) then
        call fms_mpp_set_current_pelist(Ice%slow_pelist)
        call coupler_chksum_obj%get_slow_ice_chksums('coupler_end', 0)
      end if
    endif
    call fms_mpp_set_current_pelist()

    !> @parblock
    !! CHECK TIME_CURRENT == TIME_END
    !! @endparblock
    if (Time_current /= Time_end) call fms_error_mesg ('program coupler',  &
         'final time does not match expected ending time', WARNING)

    !> @parblock
    !! CALL OCEAN_MODEL_END
    !! @endparblock
    if (Ocean%is_ocean_pe) then
      call fms_mpp_set_current_pelist(Ocean%pelist)
      call ocean_model_end (Ocean, Ocean_state, Time_current)
    endif

    !> @parblock
    !! CALL ATMOS_MODEL_END
    !! @endparblock
    if (Atm%pe) then
      call fms_mpp_set_current_pelist(Atm%pelist)
      call atmos_model_end ( Atm )
    endif

    !> @parblock
    !! CALL LAND_MODEL_END
    !! @endparblock
    if (Land%pe) then
      call fms_mpp_set_current_pelist(Land%pelist)
      call land_model_end (Atmos_land_boundary, Land)
    endif

    !> @parblock
    !! CALL ICE_MODEL_END
    !! @endparblock
    if (Ice%pe) then
      if (Ice%slow_ice_PE) then
        call fms_mpp_set_current_pelist(Ice%slow_pelist)
      else ! This must be a fast ice PE.
        call fms_mpp_set_current_pelist(Ice%fast_pelist)
      endif
      call ice_model_end (Ice)
    endif

    !> @parblock
    !! WRITE RESTART FILE
    !! @endparblock
    call coupler_restart(Atm, Ice, Ocean, Ocn_bc_restart, Ice_bc_restart, &
                         Time_current, Time_restart_current, Time_start)

    !> @parblock
    !! FINALIZE FMS DIAGNOSTICS MANAGER
    !! @endparblock
    call fms_diag_end (Time_current)
#ifdef use_deprecated_io
    call fms_io_exi
#endif

    !> @parblock
    !! END CLOCKS
    !! @endparblock
    call fms_mpp_set_current_pelist()
    call fms_mpp_clock_end(coupler_clocks%termination)

  end subroutine coupler_end

  !> @parblock
  !! Subroutine add_domain_dimension_data writes indices for the x and y dimensions
  !! into domain-decomposed fms2_io restart files.
  !! This is required so that the FMS tile-combining tool can reconstruct the global
  !! field correctly when the I/O layout is not (1,1).  Without this call the combiner
  !! cannot determine the global position of each tile's data.
  !! @endparblock
  subroutine add_domain_dimension_data(fileobj)
    type(FmsNetcdfDomainFile_t) :: fileobj
      !< is the fms2io domain decomposed fileobj
    integer, dimension(:), allocatable :: buffer
      !< is a buffer array with axis data
    integer :: is, ie ! Starting and Ending indices for data

    call fms2_io_get_global_io_domain_indices(fileobj, "xaxis_1", is, ie, indices=buffer)
    call fms2_io_write_data(fileobj, "xaxis_1", buffer)
    deallocate(buffer)

    call fms2_io_get_global_io_domain_indices(fileobj, "yaxis_1", is, ie, indices=buffer)
    call fms2_io_write_data(fileobj, "yaxis_1", buffer)
    deallocate(buffer)

  end subroutine add_domain_dimension_data


  !> @parblock
  !! Subroutine coupler_restart writes all coupler-owned restart files.
  !!
  !! Files written:
  !! - RESTART/coupler.res or RESTART/time_stamp.coupler.res (if time_stamp
  !!   is present): ASCII file containing calendar type integer, model
  !!   start date (yr,mo,day,hr,min,sec), and current model date.
  !! - RESTART/coupler.intermediate.res or RESTART/time_stamp.coupler.intermediate.res (if time_stamp
  !!   is present: ASCII file with the time of the most recent intermediate restart.
  !!   Written only if Time_restart_current > Time_start.
  !! - Ocean boundary-condition fields: registered via fms_coupler_type_register_restarts
  !! - Ice boundary-condition fields (Ice%ocean_fluxes): registered via fms_coupler_type_register_restarts
  !!
  !! The optional time_stamp argument, when present, prefixes all file names so tha
  !! multiple intermediate restart sets can coexist in RESTART/.
  !! @endparblock
  subroutine coupler_restart(Atm, Ice, Ocean, Ocn_bc_restart, Ice_bc_restart, &
                            Time_current, Time_restart_current, Time_start, time_stamp)

    implicit none

    type(atmos_data_type), intent(inout) :: Atm
      !< is the atmospheric component derived type
    type(ice_data_type), intent(inout) :: Ice
      !< is the ice component derived type
    type(ocean_public_type), intent(inout) :: Ocean
      !< is the ocean component derived type

    type(FmsNetcdfDomainFile_t), dimension(:), pointer, intent(inout) :: Ocn_bc_restar
      !< is the ocean boundary condition restart fms2_io fileobj
    type(FmsNetcdfDomainFile_t), dimension(:), pointer, intent(inout) :: Ice_bc_restar
      !< is the ice boundary condition restart fms2_io fileobj
    type(FmsTime_type), intent(in) :: Time_curren
      !< is the current model runtime (Time)
    type(FmsTime_type), intent(in) :: Time_restart_curren
      !< is the current restart time
    type(FmsTime_type), intent(in) :: Time_star
      !< is the model start time
    character(len=*), intent(in),  optional :: time_stamp
      !< is used to determine the restart file as 'RESTART/time_stamp/coupler.res'
      !! and 'RESTART/time_stamp/coupler.intermediate.res'.

    character(len=128) :: file_run, file_res

    integer :: yr, mon, day, hr, min, sec, date(6), n
    integer ::  num_ice_bc_restart, num_ocn_bc_restar
    integer :: restart_unit ! Unit for the coupler restart file

    call fms_mpp_set_current_pelist()

    !> @parblock
    !! SET COUPLER.RES AND COUPLER.INTERMEDIATE.RES FILE NAMES
    !! @endparblock
    if (present(time_stamp)) then
      file_run = 'RESTART/'//trim(time_stamp)//'.coupler.res'
      file_res = 'RESTART/'//trim(time_stamp)//'.coupler.intermediate.res'
    else
      file_run = 'RESTART/coupler.res'
      file_res = 'RESTART/coupler.intermediate.res'
    endif

    !> @parblock
    !! WRITE TIME TO COUPLER.RES
    !! @endparblock
    call fms_time_manager_get_date (Time_current, date(1), date(2), date(3), date(4), date(5), date(6))
    if ( fms_mpp_pe().EQ.fms_mpp_root_pe()) then
       open(newunit = restart_unit, file=file_run, status='replace', form='formatted')
       write(restart_unit, '(i6,8x,a)' ) calendar_type, &
            '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'
       write(restart_unit, '(6i6,8x,a)' )date_init, 'Model start time:   year, month, day, hour, minute, second'
       write(restart_unit, '(6i6,8x,a)' )date, 'Current model time:  year, month, day, hour, minute, second'
       close(restart_unit)
    endif

    !> @parblock
    !! WRITE DATE TO COUPLER.INTERMEDIATE.RES
    !! @endparblock
    if (Time_restart_current > Time_start) then
      if ( fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        open(newunit = restart_unit, file=file_res, status='replace', form='formatted')
        call fms_time_manager_get_date(Time_restart_current, yr,mon,day,hr,min,sec)
        write(restart_unit, '(6i6,8x,a)' )yr,mon,day,hr,min,sec, &
             'Current intermediate restart time:  year, month, day, hour, minute, second'
        close(restart_unit)
      endif
    endif

    !> @parblock
    !! WRITE OCEAN%FIELDS RESTARTS
    !! @endparblock
    if (Ocean%is_ocean_pe) then
      call fms_mpp_set_current_pelist(Ocean%pelist)
      if (associated(Ocn_bc_restart)) deallocate(Ocn_bc_restart)
      call fms_coupler_type_register_restarts(Ocean%fields, Ocn_bc_restart, &
               num_ocn_bc_restart, Ocean%domain, to_read=.false., ocean_restart=.true., directory="RESTART/")
      do n = 1, num_ocn_bc_restar
         if (fms2_io_check_if_open(Ocn_bc_restart(n))) then
             call fms2_io_write_restart(Ocn_bc_restart(n))
             call add_domain_dimension_data(Ocn_bc_restart(n))
             call fms2_io_close_file(Ocn_bc_restart(n))
          endif
       enddo
    endif

    !> @parblock
    !! WRITE ICE%OCEAN_FLUXES RESTARTS
    !! @endparblock
    if (Atm%pe) then
      call fms_mpp_set_current_pelist(Atm%pelist)
      if (associated(Ice_bc_restart)) deallocate(Ice_bc_restart)
      call fms_coupler_type_register_restarts(Ice%ocean_fluxes, Ice_bc_restart, &
           num_ice_bc_restart, Ice%slow_domain_NH, to_read=.false., ocean_restart=.false., directory="RESTART/")
      do n = 1, num_ice_bc_restar
        if (fms2_io_check_if_open(Ice_bc_restart(n))) then
          call fms2_io_write_restart(Ice_bc_restart(n))
          call add_domain_dimension_data(Ice_bc_restart(n))
          call fms2_io_close_file(Ice_bc_restart(n))
        endif
      enddo
    endif

  end subroutine coupler_restar

  !> @parblock
  !! Subroutine get_coupler_chksums computes chksums with fms_mpp_chksums
  !! for the following:
  !! Atmosphere (Atm) fields:
  !! - atm%t_bot (temperature at bottom)
  !! - atm%z_bot (height at bottom)
  !! - atm%p_bot (pressure at bottom)
  !! - atm%u_bot (u-wind at bottom)
  !! - atm%v_bot (v-wind at bottom)
  !! - atm%p_surf (surface pressure)
  !! - atm%gust (gustiness)
  !! - atm%tr_bot (atmospheric tracers - dynamically included)
  !!
  !! Land fields:
  !! - land%t_surf (surface temperature)
  !! - land%t_ca (canopy air temperature)
  !! - land%rough_mom (momentum roughness)
  !! - land%rough_heat (heat roughness)
  !! - land%rough_scale (roughness scale)
  !! - land%tr (land tracers - dynamically included)
  !!
  !! Ice fields:
  !! - ice%t_surf (surface temperature)
  !! - ice%rough_mom (momentum roughness)
  !! - ice%rough_heat (heat roughness)
  !! - ice%rough_moist (moisture roughness)
  !! - ice%ocean_fields (ocean-ice boundary fields)
  !! @endparblock
  subroutine get_coupler_chksums(this, id, timestep)

    implicit none

    class(coupler_chksum_type), intent(in) :: this
      !< is a reference to self (coupler_chksum_type object)
    character(len=*), intent(in) :: id
      !< id to label CHECKSUMS in stdout, e.g., 'coupler_init+', 'top_of_coupled_loop+', 'coupler_end-', etc
    integer, intent(in) :: timestep
      !< timestep to label CHECKSUMS in stdou

    type :: tracer_ind_type
      integer :: atm, ice, lnd ! indices of the tracer in the respective models
    end type tracer_ind_type

    integer :: n_atm_tr, n_lnd_tr, n_exch_tr
    integer :: n_atm_tr_tot, n_lnd_tr_to
    integer :: i, tr, n, m, outuni
    type(tracer_ind_type), allocatable :: tr_table(:)
    character(32) :: tr_name

    call fms_tracer_manager_get_number_tracers (MODEL_ATMOS, num_tracers=n_atm_tr_tot, num_prog=n_atm_tr)
    call fms_tracer_manager_get_number_tracers (MODEL_LAND, num_tracers=n_lnd_tr_tot, num_prog=n_lnd_tr)

    ! Assemble the table of tracer number translation by matching names of
    ! prognostic tracers in the atmosphere and surface models; skip all atmos.
    ! tracers that have no corresponding surface tracers.
    allocate(tr_table(n_atm_tr))
    n = 1
    do i = 1,n_atm_tr
      call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, i, tr_name )
      tr_table(n)%atm = i
      tr_table(n)%ice = fms_tracer_manager_get_tracer_index ( MODEL_ICE,  tr_name )
      tr_table(n)%lnd = fms_tracer_manager_get_tracer_index ( MODEL_LAND, tr_name )
      if (tr_table(n)%ice/=NO_TRACER .or. tr_table(n)%lnd/=NO_TRACER) n = n+1
    enddo
    n_exch_tr = n-1

100 FORMAT("CHECKSUM::",A32," = ",Z20)
101 FORMAT("CHECKSUM::",A16,a,'%',a," = ",Z20)

    if (this%components%Atm%pe) then
      call fms_mpp_set_current_pelist(this%components%Atm%pelist)

      outunit = fms_mpp_stdout()
      write(outunit,*) 'BEGIN CHECKSUM(Atm):: ', id, timestep
      write(outunit,100) 'atm%t_bot',  fms_mpp_chksum(this%components%Atm%t_bot)
      write(outunit,100) 'atm%z_bot',  fms_mpp_chksum(this%components%Atm%z_bot)
      write(outunit,100) 'atm%p_bot',  fms_mpp_chksum(this%components%Atm%p_bot)
      write(outunit,100) 'atm%u_bot',  fms_mpp_chksum(this%components%Atm%u_bot)
      write(outunit,100) 'atm%v_bot',  fms_mpp_chksum(this%components%Atm%v_bot)
      write(outunit,100) 'atm%p_surf', fms_mpp_chksum(this%components%Atm%p_surf)
      write(outunit,100) 'atm%gust',   fms_mpp_chksum(this%components%Atm%gust)
      do tr = 1,n_exch_tr
         n = tr_table(tr)%atm
         if (n /= NO_TRACER) then
            call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
            write(outunit,100) 'atm%'//trim(tr_name), fms_mpp_chksum(this%components%Atm%tr_bot(:,:,n))
          endif
       enddo

      write(outunit,100) 'land%t_surf', fms_mpp_chksum(this%components%Land%t_surf)
      write(outunit,100) 'land%t_ca',   fms_mpp_chksum(this%components%Land%t_ca)
      write(outunit,100) 'land%rough_mom',   fms_mpp_chksum(this%components%Land%rough_mom)
      write(outunit,100) 'land%rough_heat',  fms_mpp_chksum(this%components%Land%rough_heat)
      write(outunit,100) 'land%rough_scale', fms_mpp_chksum(this%components%Land%rough_scale)
      do tr = 1,n_exch_tr
        n = tr_table(tr)%lnd
        if (n /= NO_TRACER) then
          call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
#ifndef _USE_LEGACY_LAND_
          write(outunit,100) 'land%'//trim(tr_name), fms_mpp_chksum(this%components%Land%tr(:,:,n))
#else
          write(outunit,100) 'land%'//trim(tr_name), fms_mpp_chksum(this%components%Land%tr(:,:,:,n))
#endif
        endif
      enddo

      write(outunit,100) 'ice%t_surf', fms_mpp_chksum(this%components%Ice%t_surf)
      write(outunit,100) 'ice%rough_mom', fms_mpp_chksum(this%components%Ice%rough_mom)
      write(outunit,100) 'ice%rough_heat', fms_mpp_chksum(this%components%Ice%rough_heat)
      write(outunit,100) 'ice%rough_moist', fms_mpp_chksum(this%components%Ice%rough_moist)
      write(outunit,*) 'STOP CHECKSUM(Atm):: ', id, timestep

    !if (Ocean%is_ocean_pe) call mpp_set_current_pelist(Ocean%pelist)

      write(outunit,*) 'BEGIN CHECKSUM(Ice):: ', id, timestep
      call fms_coupler_type_write_chksums(this%components%Ice%ocean_fields, outunit, 'ice%')
      write(outunit,*) 'STOP CHECKSUM(Ice):: ', id, timestep

    endif

    deallocate(tr_table)

    call fms_mpp_set_current_pelist()

  end subroutine get_coupler_chksums

  !#######################################################################

  !> @parblock
  !! Subroutine get_atmos_ice_land_ocean_chksums calls get_atmos_ice_land_chksums
  !! and get_ocean_chksums
  !! @endparblock
  subroutine get_atmos_ice_land_ocean_chksums(this, id, timestep)

    implicit none

    class(coupler_chksum_type), intent(in) :: this
      !< is a reference to self (coupler_chksum_type object)
    character(len=*), intent(in) :: id
      !< is the id labelling the set of checksums in the logfile
    integer, intent(in) :: timestep
      !< is the timestep

    if (this%components%Atm%pe) then
      call fms_mpp_set_current_pelist(this%components%Atm%pelist)
      call this%get_atmos_ice_land_chksums(trim(id), timestep)
    endif
    if (this%components%Ocean%is_ocean_pe) then
      call fms_mpp_set_current_pelist(this%components%Ocean%pelist)
      call this%get_ocean_chksums(trim(id), timestep)
    endif

    call fms_mpp_set_current_pelist()

  end subroutine get_atmos_ice_land_ocean_chksums

  !> @parblock
  !! Subroutine get_atmos_ice_land_chksums computes and prints checksums for
  !! atmosphere, fast-ice, and land fields.
  !!
  !! The pelist must be set (synchronize) before calling this subroutine:
  !! if (Atm%pe) then
  !!    call fms_mpp_set_current_pelist(Atm%pelist)
  !!    call coupler_chksum_obj%get_atmos_ice_land_chksums('MAIN_LOOP-', nc)
  !! endif
  !! @endparblock
  subroutine get_atmos_ice_land_chksums(this, id, timestep)

    class(coupler_chksum_type), intent(in) :: this !< self
    character(len=*), intent(in) :: id       !< id to label CHECKSUMS in stdou
    integer         , intent(in) :: timestep !< timestep

    call atmos_data_type_chksum(     id, timestep, this%components%Atm)
    call lnd_ice_atm_bnd_type_chksum(id, timestep, this%components%Land_ice_atmos_boundary)

    if (this%components%Ice%fast_ice_pe) then
      call fms_mpp_set_current_pelist(this%components%Ice%fast_pelist)
      call ice_data_type_chksum(   id, timestep, this%components%Ice)
      call atm_ice_bnd_type_chksum(id, timestep, this%components%Atmos_ice_boundary)
    endif
    if (this%components%Land%pe) then
      call fms_mpp_set_current_pelist(this%components%Land%pelist)
      call land_data_type_chksum(  id, timestep, this%components%Land)
      call atm_lnd_bnd_type_chksum(id, timestep, this%components%Atmos_land_boundary)
    endif

    call fms_mpp_set_current_pelist(this%components%Atm%pelist)

  end subroutine get_atmos_ice_land_chksums

  !> @parblock
  !! Subroutine get_slow_ice_chksums calls subroutine that will print ou
  !! checksums for slow ice and ocean-ice boundary fields.
  !! The pelist must be set (synchronize) before calling this subroutine:
  !! if (Ice%slow_ice_pe) then
  !!    call mpp_set_current_pelist(Ice%slow_pelist)
  !!    call slow_ice_chksum('MAIN_LOOP-', nc)
  !! endif
  !! @endparblock
  subroutine get_slow_ice_chksums(this, id, timestep)

    class(coupler_chksum_type), intent(in) :: this !< self
    character(len=*), intent(in) :: id       !<id to label CHECKSUMS in stdou
    integer         , intent(in) :: timestep !< timestep

    call ice_data_type_chksum(    id, timestep, this%components%Ice)
    call ocn_ice_bnd_type_chksum( id, timestep, this%components%Ocean_ice_boundary)

  end subroutine get_slow_ice_chksums

  !> @parblock
  !! Subroutine get_ocean_chksums calls subroutine that will print ou
  !! checksums for ocean and ice-ocean boundary fields.
  !! The pelist must be set (synchronize) before calling this subroutine:
  !! if (Ocean%is_ocean_pe) then
  !!    call mpp_set_current_pelist(Ocean%pelist)
  !!    call ocean_chksum('MAIN_LOOP-', nc)
  !! endif
  !! @endparblock
  subroutine get_ocean_chksums(this, id, timestep)

    class(coupler_chksum_type), intent(in) :: this !< self
    character(len=*), intent(in) :: id       !< ID labelling the set of CHECKSUMS
    integer         , intent(in) :: timestep !< Timestep

    call ocean_public_type_chksum(id, timestep, this%components%Ocean)
    call ice_ocn_bnd_type_chksum( id, timestep, this%components%Ice_ocean_boundary)

  end subroutine get_ocean_chksums

  !> @parblock
  !! Subroutine coupler_set_clock_ids registers all FMS performance-clock IDs for the
  !! coupled model and stores them in the coupler_clocks struct.
  !!
  !! Clocks are registered on the PE list most appropriate for each phase: atmosphere
  !! clocks on Atm%pelist, ocean clocks on Ocean%pelist, ice clocks on Ice%fast_pelis
  !! or Ice%slow_pelist, and ocean-ice flux clocks slow_ice_ocean_pelist.
  !! Global clocks (main loop, termination, flux_check_stocks) are registered on all PEs.
  !! This routine must be called after PE lists have been set up but before any clock is started.
  !! @endparblock
  subroutine coupler_set_clock_ids(coupler_clocks, Atm, Land, Ice, Ocean, ensemble_pelist,&
                                   slow_ice_ocean_pelist, ensemble_id)

    implicit none
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< is a derived type containing clocks for profiling
    type(atmos_data_type), intent(in) :: Atm
      !< is the atm derived type, required for atm pelis
    type(land_data_type),  intent(in) :: Land
      !< is the land derived type, required for land pelis
    type(ocean_public_type), intent(in) :: Ocean
      !< is the ocean derived type, required for ocean pelis
    type(ice_data_type), intent(in) :: Ice
      !< is the ice derived type, required for ice pelis
    integer, dimension(:), intent(in) :: slow_ice_ocean_pelis
      !< is the slow_ice_ocean_pelist, required for slow_ice_ocean pelis
    integer, dimension(:,:), intent(in) :: ensemble_pelis
      !< is the ensemble_pelist, to register clocks for ensemble members
    integer, intent(in) :: ensemble_id
      !< is the ensemble_id used as index in ensemble_pelis

    if (Atm%pe) then
      call fms_mpp_set_current_pelist(Atm%pelist)
      coupler_clocks%atmos_model_init = fms_mpp_clock_id( '  Init: atmos_model_init ' )
    endif
    if (Land%pe) then
      call fms_mpp_set_current_pelist(Land%pelist)
      coupler_clocks%land_model_init  = fms_mpp_clock_id( '  Init: land_model_init ' )
    endif
    if (Ice%pe) then
      if (Ice%shared_slow_fast_PEs) then ; call fms_mpp_set_current_pelist(Ice%pelist)
      elseif (Ice%fast_ice_pe)      then ; call fms_mpp_set_current_pelist(Ice%fast_pelist)
      elseif (Ice%slow_ice_pe)      then ; call fms_mpp_set_current_pelist(Ice%slow_pelist)
      else ; call fms_mpp_error(FATAL, "All Ice%pes must be a part of Ice%fast_ice_pe or Ice%slow_ice_pe")
      endif
      coupler_clocks%ice_model_init   = fms_mpp_clock_id( '  Init: ice_model_init ' )
    endif
    if (Ocean%is_ocean_pe) then
      call fms_mpp_set_current_pelist(Ocean%pelist)
      coupler_clocks%ocean_model_init = fms_mpp_clock_id( '  Init: ocean_model_init ' )
    endif
    call fms_mpp_set_current_pelist(ensemble_pelist(ensemble_id,:))
    coupler_clocks%flux_exchange_init = fms_mpp_clock_id( '  Init: flux_exchange_init' )

    call fms_mpp_set_current_pelist()
    coupler_clocks%main = fms_mpp_clock_id( 'Main loop' )
    coupler_clocks%termination = fms_mpp_clock_id( 'Termination' )

    If(Atm%pe) then
      call fms_mpp_set_current_pelist(Atm%pelist)
      coupler_clocks%generate_sfc_xgrid = fms_mpp_clock_id( 'generate_sfc_xgrid' )
    end if
    if (Ice%slow_ice_PE .or. Ocean%is_ocean_pe) then
      call fms_mpp_set_current_pelist(slow_ice_ocean_pelist)
      coupler_clocks%flux_ocean_to_ice = fms_mpp_clock_id( 'flux_ocean_to_ice' )
      coupler_clocks%flux_ice_to_ocean = fms_mpp_clock_id( 'flux_ice_to_ocean' )
    endif
    if (Atm%pe) then
      call fms_mpp_set_current_pelist(Atm%pelist)
      coupler_clocks%atm         = fms_mpp_clock_id( 'ATM' )
      coupler_clocks%atmos_loop  = fms_mpp_clock_id( ' ATM: atmos loop' )
      coupler_clocks%atmos_tracer_driver_gather_data  &
          = fms_mpp_clock_id( '  A-L: atmos_tracer_driver_gather_data' )
      coupler_clocks%sfc_boundary_layer           = fms_mpp_clock_id( '  A-L: sfc_boundary_layer' )
      coupler_clocks%update_atmos_model_dynamics  = fms_mpp_clock_id( '  A-L: update_atmos_model_dynamics')
      if (.not. do_concurrent_radiation) &
          coupler_clocks%radiation            = fms_mpp_clock_id( '  A-L: serial radiation' )
      coupler_clocks%update_atmos_model_down  = fms_mpp_clock_id( '  A-L: update_atmos_model_down' )
      coupler_clocks%flux_down_from_atmos     = fms_mpp_clock_id( '  A-L: flux_down_from_atmos' )
      coupler_clocks%update_land_model_fast   = fms_mpp_clock_id( '  A-L: update_land_model_fast' )
      coupler_clocks%update_ice_model_fast    = fms_mpp_clock_id( '  A-L: update_ice_model_fast' )
      coupler_clocks%flux_up_to_atmos         = fms_mpp_clock_id( '  A-L: flux_up_to_atmos' )
      coupler_clocks%update_atmos_model_up    = fms_mpp_clock_id( '  A-L: update_atmos_model_up' )
      if (do_concurrent_radiation) then
        coupler_clocks%radiation             = fms_mpp_clock_id( '  A-L: concurrent radiation' )
        coupler_clocks%concurrent_atmos      = fms_mpp_clock_id( '  A-L: concurrent atmos' )
      endif
      coupler_clocks%update_atmos_model_state  = fms_mpp_clock_id( '  A-L: update_atmos_model_state')
      coupler_clocks%update_land_model_slow    = fms_mpp_clock_id( ' ATM: update_land_model_slow' )
      coupler_clocks%flux_land_to_ice          = fms_mpp_clock_id( ' ATM: flux_land_to_ice' )
    endif
    if (Ice%pe) then
      if (Ice%fast_ice_pe) call fms_mpp_set_current_pelist(Ice%fast_pelist)
      coupler_clocks%set_ice_surface_fast       = fms_mpp_clock_id( ' Ice: set_ice_surface fast' )
      coupler_clocks%update_ice_model_slow_fast = fms_mpp_clock_id( ' Ice: update_ice_model_slow fast' )

      if (Ice%slow_ice_pe) call fms_mpp_set_current_pelist(Ice%slow_pelist)
      coupler_clocks%set_ice_surface_slow       = fms_mpp_clock_id( ' Ice: set_ice_surface slow' )
      coupler_clocks%update_ice_model_slow_slow = fms_mpp_clock_id( ' Ice: update_ice_model_slow slow' )
      coupler_clocks%flux_ice_to_ocean_stocks   = fms_mpp_clock_id( ' Ice: flux_ice_to_ocean_stocks' )

      call fms_mpp_set_current_pelist(Ice%pelist)
      coupler_clocks%set_ice_surface_exchange       = fms_mpp_clock_id( ' Ice: set_ice_surface exchange' )
      coupler_clocks%update_ice_model_slow_exchange = fms_mpp_clock_id( ' Ice: update_ice_model_slow exchange' )

    endif
    if (Ocean%is_ocean_pe) then
      call fms_mpp_set_current_pelist(Ocean%pelist)
      coupler_clocks%ocean = fms_mpp_clock_id( 'OCN' )
    endif

    call fms_mpp_set_current_pelist()
    coupler_clocks%flux_check_stocks       = fms_mpp_clock_id( 'flux_check_stocks' )
    coupler_clocks%intermediate_restart    = fms_mpp_clock_id( 'intermediate restart' )
    coupler_clocks%final_flux_check_stocks = fms_mpp_clock_id( 'final flux_check_stocks' )

  end subroutine coupler_set_clock_ids

  !> @parblock
  !! Subroutine coupler_flux_init_finish_stocks initializes or finalizes stock computation
  !! to check for water, heat, and salt conservation.
  !!
  !! - When init_stocks=.true., calls flux_init_stocks to establish the baseline
  !!   globally integrated water, heat, and salt stocks (q_start) for all four
  !!   component models at the start of the run.
  !! - When finish_stocks=.true., calls flux_check_stocks (if check_stocks >= 0) to
  !!   compute final stocks, compares them to q_start, and reports conservation errors
  !!   to the stocks output file.
  !! @endparblock
  subroutine coupler_flux_init_finish_stocks(Time, Atm, Land, Ice, Ocean_state, &
                                             coupler_clocks, init_stocks, finish_stocks)

    implicit none

    type(FmsTime_type), intent(in) :: Time
      !< is the current model time
    type(atmos_data_type), intent(inout) :: Atm
      !< is the atmosphere derived type
    type(land_data_type), intent(inout) :: Land
      !< is the land derived type
    type(ice_data_type), intent(inout) :: Ice
      !< is the ice derived type
    type(ocean_state_type), pointer, intent(inout) :: Ocean_state
      !< is the ocean state derived type
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< contains the clocks for profiling
    logical, optional, intent(in) :: init_stocks
      !< is a flag where if true, call flux_init_stocks
    logical, optional, intent(in) :: finish_stocks
      !< is a flag where if true, call final flux_check_stocks

    logical :: init, finish
      ! local control flags: init defaults to false unless init_stocks is provided;
      ! finish defaults to false unless finish_stocks is provided

    init=.False.   ; if(present(init_stocks)) init=init_stocks
    finish=.False. ; if(present(finish_stocks)) finish=finish_stocks

    if(init) then
      call fms_mpp_set_current_pelist()
      call flux_init_stocks(Time, Atm, Land, Ice, Ocean_state)
    else if(finish) then
      call fms_mpp_set_current_pelist()
      call fms_mpp_clock_begin(coupler_clocks%final_flux_check_stocks)
      if (check_stocks >= 0) then
        call fms_mpp_set_current_pelist()
        call flux_check_stocks(Time=Time, Atm=Atm, Lnd=Land, Ice=Ice, Ocn_state=Ocean_state)
      endif
      call fms_mpp_clock_end(coupler_clocks%final_flux_check_stocks)
    else
      call fms_mpp_error(FATAL, 'coupler_flux_init_finish_stocks: either init or finish needs to be .True.')
    end if

  end subroutine coupler_flux_init_finish_stocks

  !> @parblock
  !! Subroutine coupler_flux_check_stocks periodically computes
  !! water, heat, and salt stocks for all four components.
  !! @endparblock
  subroutine coupler_flux_check_stocks(nc, Time, Atm, Land, Ice, Ocean_state, coupler_clocks)

    implicit none

    integer, intent(in) :: nc
      !< is the current outer-loop timestep
    type(FmsTime_type), intent(in) :: Time
      !< is the current model time
    type(atmos_data_type), intent(inout) :: Atm
      !< is the atmosphere componen
    type(land_data_type), intent(inout) :: Land
      !< is the land componen
    type(ice_data_type), intent(inout) :: Ice
      !< is the ice componen
    type(ocean_state_type), pointer, intent(inout) :: Ocean_state
      !< is the ocean state componen
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< are the coupler clocks

    call fms_mpp_clock_begin(coupler_clocks%flux_check_stocks)
    if (check_stocks*((nc-1)/check_stocks) == nc-1 .AND. nc > 1) then
      call fms_mpp_set_current_pelist()
      call flux_check_stocks(Time=Time, Atm=Atm, Lnd=Land, Ice=Ice, Ocn_state=Ocean_state)
    endif
    call fms_mpp_clock_end(coupler_clocks%flux_check_stocks)

  end subroutine coupler_flux_check_stocks

  !> @parblock
  !! Subroutine coupler_flux_ocean_to_ice calls flux_ocean_to_ice to
  !! transfers the current ocean state (SST, surface currents, salinity, sea-surface height)
  !! into the Ocean_ice_boundary derived type in preparation for the slow-ice update.
  !! The call occurs on slow_ice_ocean_pelist (the union of slow-ice and ocean PEs)
  !! and is profiled with coupler_clocks%flux_ocean_to_ice.
  !! @endparblock
  subroutine coupler_flux_ocean_to_ice(Ocean, Ice, Ocean_ice_boundary, coupler_clocks, slow_ice_ocean_pelist)

    implicit none

    type(ocean_public_type), intent(inout) :: Ocean
      !< is the ocean componen
    type(ice_data_type), intent(in) :: Ice
      !< is the ice componen
    type(ocean_ice_boundary_type), intent(inout) :: Ocean_ice_boundary
      !< is the ocean-ice boundary componen
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< are the coupler clocks
    integer, dimension(:), intent(in) :: slow_ice_ocean_pelis
      !< is the slow ice-ocean PE lis

    call fms_mpp_set_current_pelist(slow_ice_ocean_pelist)
    call fms_mpp_clock_begin(coupler_clocks%flux_ocean_to_ice)

    call flux_ocean_to_ice(Ocean, Ice, Ocean_ice_boundary)

    call fms_mpp_clock_end(coupler_clocks%flux_ocean_to_ice)

  end subroutine coupler_flux_ocean_to_ice

  !> @parblock
  !! Subroutine coupler_flux_ice_to_ocean updates the accumulated ice-to-ocean
  !! forcing fluxes (heat, freshwater, salt, momentum, shortwave) in Ice_ocean_boundary
  !! in preparation for the ocean model update.
  !!
  !! The optional set_current_slow_ice_ocean_pelist flag controls whether
  !! fms_mpp_set_current_pelist(slow_ice_ocean_pelist) is called.  It defaults to
  !! .false. because when this routine follows coupler_flux_ocean_to_ice, the PE
  !! list is already set to slow_ice_ocean_pelist by that routine.  Pass .true.
  !! when calling coupler_flux_ice_to_ocean independently (e.g., in lag-flux mode).
  !! @endparblock
  subroutine coupler_flux_ice_to_ocean(Ice, Ocean, Ice_ocean_boundary, coupler_clocks, &
                                       slow_ice_ocean_pelist, set_current_slow_ice_ocean_pelist)

    implicit none

    type(ice_data_type), intent(inout)  :: Ice
      !< is the Ice componen
    type(ocean_public_type), intent(inout)  :: Ocean
      !< is the Ocean componen
    type(ice_ocean_boundary_type), intent(inout) :: Ice_ocean_boundary
      !< is the Ice_ocean_boundary componen
    type(coupler_clock_type), intent(inout) :: coupler_clocks
     !< are the coupler_clocks
    integer, dimension(:), optional, intent(in) :: slow_ice_ocean_pelis
      !< is the slow_ice_ocean_pelis
    logical, optional, intent(in) :: set_current_slow_ice_ocean_pelis
      !< is a flag where if true, call mpp_set_current_pelist(slow_ice_ocean_pelist)

    logical :: set_current_slow_ice_ocean_pelist_in ! .F. by default; set to equal set_current_slow_ice_ocean_pelis

    ! mpp_set_current_pelist(slow_ice_ocean_pelist) is not required if coupler_flux_ice_to_ocean is called after
    ! coupler_flux_ocean_to_ice since mpp_set_current_pelist(slow_ice_ocean_pelist) is called
    ! in coupler_flux_ocean_to_ice
    set_current_slow_ice_ocean_pelist_in=.False.
    if(present(set_current_slow_ice_ocean_pelist)) &
        set_current_slow_ice_ocean_pelist_in = set_current_slow_ice_ocean_pelis

    ! Update Ice_ocean_boundary; the first iteration is supplied by restarts

    if(set_current_slow_ice_ocean_pelist_in) call fms_mpp_set_current_pelist(slow_ice_ocean_pelist)

    call fms_mpp_clock_begin(coupler_clocks%flux_ice_to_ocean)
    call flux_ice_to_ocean(Ice, Ocean, Ice_ocean_boundary)
    call fms_mpp_clock_end(coupler_clocks%flux_ice_to_ocean)

  end subroutine coupler_flux_ice_to_ocean

  !> @parblock
  !! Subroutine coupler_unpack_ocean_ice_boundary, called after coupler_flux_ocean_to_ice,
  !! first calls flux_ocean_to_ice_finish to override data (if field exists in data_table)
  !! and then calls unpack_ocean_ice_boundary to unpack the ocean-ice boundary data into the ice model state.
  !! slow_ice_chksums are computed if do_chksum is true.
  !! @endparblock
  subroutine coupler_unpack_ocean_ice_boundary(nc, Time_flux_ocean_to_ice, Ice, Ocean_ice_boundary, coupler_clocks, &
                                               coupler_chksum_obj)

    implicit none

    integer, intent(in) :: nc
      !< is the current outer loop timestep
    type(FmsTime_type),  intent(inout) :: Time_flux_ocean_to_ice
      !< is the time for flux_ocean_to_ice
    type(ice_data_type), intent(inout) :: Ice
      !< is the Ice componen
    type(ocean_ice_boundary_type), intent(inout) :: Ocean_ice_boundary
      !< is the Ocean_ice_boundary
    type(coupler_clock_type),      intent(inout) :: coupler_clocks
      !< are the coupler_clocks
    type(coupler_chksum_type),     intent(in)  :: coupler_chksum_obj
      !< is used for computing slow-ice checksums when do_chksum=.true.

    call fms_mpp_set_current_pelist(Ice%slow_pelist)
    call fms_mpp_clock_begin(coupler_clocks%set_ice_surface_slow)

    call flux_ocean_to_ice_finish( Time_flux_ocean_to_ice, Ice, Ocean_Ice_Boundary )
    call unpack_ocean_ice_boundary( Ocean_ice_boundary, Ice )
    if (do_chksum) call coupler_chksum_obj%get_slow_ice_chksums('update_ice_slow+', nc)

    call fms_mpp_clock_end(coupler_clocks%set_ice_surface_slow)

  end subroutine coupler_unpack_ocean_ice_boundary

  !> @parblock
  !! Subroutine coupler_exchange_slow_to_fast_ice transfers updated ocean boundary state
  !! to fast ice procesess by calling exchange_slow_to_fast_ice from ice_model_mod.
  !! This subroutine is called after coupler_flux_ocean_to_ice and coupler_unpack_ocean_ice_boundary
  !! @endparblock
  subroutine coupler_exchange_slow_to_fast_ice(Ice, coupler_clocks)

    implicit none
    type(ice_data_type), intent(inout) :: Ice
      !< is the Ice componen
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< are the coupler_clocks

    ! This could be a point where the model is serialized if the fast and
    ! slow ice are on different PEs.
    if (.not.Ice%shared_slow_fast_PEs) call fms_mpp_set_current_pelist(Ice%pelist)
    call fms_mpp_clock_begin(coupler_clocks%set_ice_surface_exchange)
    call exchange_slow_to_fast_ice(Ice)
    call fms_mpp_clock_end(coupler_clocks%set_ice_surface_exchange)

  end subroutine coupler_exchange_slow_to_fast_ice

  !> @parblock
  !! Subroutine coupler_exchange_fast_to_slow_ice calls exchange_fast_to_slow_ice form
  !! ice_model_mod to copy fast part of sea-ice to slow part of sea-ice.
  !!
  !! The optional set_ice_current_pelist flag, when .true., calls
  !! fms_mpp_set_current_pelist(Ice%pelist) to set and synchronize the pes in the pelist.
  !! @endparblock
  subroutine coupler_exchange_fast_to_slow_ice(Ice, coupler_clocks, set_ice_current_pelist)

    implicit none
    type(ice_data_type), intent(inout) :: Ice
      !< is the ice componen
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< are the coupler_clocks
    logical, optional, intent(in) :: set_ice_current_pelis
     !< is a flag where if true, call fms_mpp_set_current_pelist(Ice%pelist)

    logical :: set_ice_current_pelist_in

    set_ice_current_pelist_in = .False.
    if(present(set_ice_current_pelist)) set_ice_current_pelist_in = set_ice_current_pelis

    if(set_ice_current_pelist_in .and. .not.Ice%shared_slow_fast_PEs) call fms_mpp_set_current_pelist(Ice%pelist)
    call fms_mpp_clock_begin(coupler_clocks%update_ice_model_slow_exchange)
    call exchange_fast_to_slow_ice(Ice)
    call fms_mpp_clock_end(coupler_clocks%update_ice_model_slow_exchange)

  end subroutine coupler_exchange_fast_to_slow_ice

  !> @parblock
  !! Subroutine coupler_set_ice_surface_fields calls set_ice_surface_fields from
  !! ice_model_mod to prepare the ice surface state for atmosphere fast physics,
  !! as well as pre-calculate ice radiative properties.
  !! @endparblock
  subroutine coupler_set_ice_surface_fields(Ice, coupler_clocks)

    implicit none
    type(ice_data_type), intent(inout) :: Ice
      !< is the Ice componen
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< are the coupler_clocks

    if (.not.Ice%shared_slow_fast_PEs) call fms_mpp_set_current_pelist(Ice%fast_pelist)
    call fms_mpp_clock_begin(coupler_clocks%set_ice_surface_fast)
    call set_ice_surface_fields(Ice)
    call fms_mpp_clock_end(coupler_clocks%set_ice_surface_fast)

  end subroutine coupler_set_ice_surface_fields

  !> @parblock
  !! Subroutine coupler_generate_sfc_xgrid calls generate_sfc_xgrid to rebuild
  !! the atmosphere-surface exchange grid (xmap_sfc) from the current land mask
  !! and ice concentration.
  !! @endparblock
  subroutine coupler_generate_sfc_xgrid(Land, Ice, coupler_clocks)

    implicit none
    type(land_data_type), intent(inout) :: Land
      !< is the Land componen
    type(ice_data_type),  intent(inout) :: Ice
      !< is the Ice componen
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< are the coupler_clocks

    call fms_mpp_clock_begin(coupler_clocks%generate_sfc_xgrid)
    call generate_sfc_xgrid( Land, Ice )
    call fms_mpp_clock_end(coupler_clocks%generate_sfc_xgrid)

  end subroutine coupler_generate_sfc_xgrid

  !> @parblock
  !! Subroutine coupler_atmos_tracer_driver_gather_data calls atmos_tracer_driver_gather_data from
  !! atmos_tracer_driver_mod to gather CO2, NH3, and tagged/isotopic NH3 at the bottom atm layer
  !! in preparation for flux exchange.
  !! @endparblock
  subroutine coupler_atmos_tracer_driver_gather_data(Atm, coupler_clocks)

    implicit none

    type(atmos_data_type), intent(inout) :: Atm
      !< is the atm componen
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< are the coupler_clocks
    call fms_mpp_clock_begin(coupler_clocks%atmos_tracer_driver_gather_data)
    call atmos_tracer_driver_gather_data(Atm%fields, Atm%tr_bot)
    call fms_mpp_clock_end(coupler_clocks%atmos_tracer_driver_gather_data)

  end subroutine coupler_atmos_tracer_driver_gather_data

  !> @parblock
  !! Subroutine coupler_sfc_boundary_layer sets the clock and calls sfc_boundary_layer
  !! to compute fluxes at the surface.  Chksum is computed if do_chksum is true.
  !! @endparblock
  subroutine coupler_sfc_boundary_layer(Atm, Land, Ice, Land_ice_atmos_boundary, &
                                        Time_atmos, current_timestep, coupler_chksum_obj, coupler_clocks)

    implicit none
    type(atmos_data_type), intent(inout) :: Atm
      !< is the atm componen
    type(land_data_type), intent(inout) :: Land
      !< is the land componen
    type(ice_data_type), intent(inout) :: Ice
      !< is the Ice componen
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_ice_atmos_boundary
      !< is the Land_ice_atmos_boundary componen
    type(FmsTime_type), intent(in) :: Time_atmos
      !< is the Atmos time
    integer, intent(in) :: current_timestep
       !< is the timestep (nc-1)*num_atmos_cal + na
    type(coupler_chksum_type), intent(in)   :: coupler_chksum_obj
      !< is the coupler_chksum_obj
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< are the coupler_clocks

    call fms_mpp_clock_begin(coupler_clocks%sfc_boundary_layer)

    call sfc_boundary_layer( real(dt_atmos), Time_atmos, Atm, Land, Ice, Land_ice_atmos_boundary )
    if(do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('sfc+', current_timestep)

    call fms_mpp_clock_end(coupler_clocks%sfc_boundary_layer)

  end subroutine coupler_sfc_boundary_layer

  !> @parblock
  !! Subroutine coupler_update_atmos_model_dynamics calls update_atmos_model_dynamics from atmos_driver
  !! to advance the atmospheric dynamical core by one timestep.  Checksums are computed when do_chksum=.true.,
  !! and memory usage is printed when do_debug=.true.
  !! @endparblock
  subroutine coupler_update_atmos_model_dynamics(Atm, current_timestep, coupler_chksum_obj, coupler_clocks)

    implicit none
    type(atmos_data_type), intent(inout) :: Atm
      !< is the Atm componen
    integer, intent(in) :: current_timestep
      !< is the current timestep
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj
      !< is the coupler_chksum_obj for computing chksums
    type(coupler_clock_type),  intent(inout) :: coupler_clocks
      !< are the coupler_clocks

    call fms_mpp_clock_begin(coupler_clocks%update_atmos_model_dynamics)
    call update_atmos_model_dynamics(Atm)
    call fms_mpp_clock_end(coupler_clocks%update_atmos_model_dynamics)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('update_atmos_model_dynamics', current_timestep)
    if (do_debug)  call fms_memutils_print_memuse_stats( 'update dyn')

  end subroutine coupler_update_atmos_model_dynamics

  !> @parblock
  !! Subroutine coupler_update_atmos_model_radiation calls update_atmos_model_radiation in
  !! atmos_driver to update the radiative heating rates, boundary radiative fluxes, and other
  !! properties. Checksums are computed if do_chksum is true and do_concurrent_radiation = .false.
  !! (due to threading restrictions in mpp_chksum). Memory usage is printed when do_debug=.true.
  !! @endparblock
  subroutine coupler_update_atmos_model_radiation(Atm, Land_ice_atmos_boundary, coupler_clocks, &
                                                  current_timestep, coupler_chksum_obj)

    implicit none

    type(atmos_data_type), intent(inout) :: Atm
      !< is the Atm componen
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_ice_atmos_boundary
      !< is the Land_ice_atmos_boundary componen
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< are the coupler_clocks
    integer,  optional, intent(in) :: current_timestep
      !< is the current timestep
    type(coupler_chksum_type), optional, intent(in) :: coupler_chksum_obj
      !< is the coupler_chksum_obj for computing chksums

    character(128) :: memuse_stats_id = 'update serial rad' !< used to label mem usage

    call fms_mpp_clock_begin(coupler_clocks%radiation)
    call update_atmos_model_radiation( Land_ice_atmos_boundary, Atm )
    call fms_mpp_clock_end(coupler_clocks%radiation)

    if(do_chksum) then
      ! cannot put mpp_chksum for concurrent_radiation as it requires the ability to have two different OpenMP threads
      ! inside of MPI at the same time which is not currently allowed
      if(.not.do_concurrent_radiation) &
          call coupler_chksum_obj%get_atmos_ice_land_chksums('update_atmos_model_radiation(ser)',current_timestep)
    end if

    if (do_debug) then
      if(do_concurrent_radiation) memuse_stats_id = 'update concurrent rad'
      call fms_memutils_print_memuse_stats(trim(memuse_stats_id))
    end if

  end subroutine coupler_update_atmos_model_radiation

  !> @parblock
  !! Subroutine coupler_update_atmos_model_down calls update_atmos_model_down from
  !! atmos_driver to execute the downward atmospheric physics sweep for heat/moisture.
  !! Checksums are computed when do_chksum=.true., and memory usage is printed when do_debug=.true.
  !! @endparblock
  subroutine coupler_update_atmos_model_down(Atm, Land_ice_atmos_boundary, current_timestep, &
                                             coupler_chksum_obj, coupler_clocks)

    implicit none
    type(atmos_data_type), intent(inout) :: Atm
      !< is the atmosphere model derived type
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_ice_atmos_boundary
      !< is the derived type containing quantities going from land and ice to atmos
    integer, intent(in) :: current_timestep
      !< is the current coupled timestep index used for checksum labelling
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj
      !< is the coupler checksum object used to compute checksums
    type(coupler_clock_type),  intent(inout) :: coupler_clocks
      !< are the coupler clocks used to measure runtime of the downward physics sweep

    call fms_mpp_clock_begin(coupler_clocks%update_atmos_model_down)
    call update_atmos_model_down( Land_ice_atmos_boundary, Atm )
    call fms_mpp_clock_end(coupler_clocks%update_atmos_model_down)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('update_atmos_down+', current_timestep)
    if (do_debug)  call fms_memutils_print_memuse_stats( 'update down')

  end subroutine coupler_update_atmos_model_down

  !> @parblock
  !! Subroutine coupler_flux_down_from_atmos calls flux_down_from_atmos to map fluxes from
  !! atmosphere to land and ice components. Runtime is measured by the clock for flux_down_from_atmos,
  !! and checksums are computed when do_chksum=.true.
  !! @endparblock
  subroutine coupler_flux_down_from_atmos(Atm, Land, Ice, Land_ice_atmos_boundary, Atmos_land_boundary, &
              Atmos_ice_boundary, Time_atmos, current_timestep, coupler_clocks, coupler_chksum_obj)

    implicit none
    type(atmos_data_type), intent(inout) :: Atm
      !< is the Atm derived type
    type(land_data_type), intent(inout) :: Land
      !< is the Land derived type
    type(ice_data_type), intent(inout) :: Ice
      !< is the Ice derived type
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_ice_atmos_boundary
      !< is the Land_ice_atmos_boundary derived type
    type(atmos_land_boundary_type), intent(inout) :: Atmos_land_boundary
      !< is the Atmos_land_boundary derived type
    type(atmos_ice_boundary_type), intent(inout) :: Atmos_ice_boundary
      !< is the Atmos_ice_boundary derived type
    type(FmsTime_type), intent(in) :: Time_atmos
      !< is the Time_atmos FmsTime_type containing time in seconds
    integer, intent(in) :: current_timestep
      !< is the current timestep
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< is the coupler_clocks
    type(coupler_chksum_type), intent(in)   :: coupler_chksum_obj
      !< is the coupler_chksum_obj for computing chksums

    call fms_mpp_clock_begin(coupler_clocks%flux_down_from_atmos)
    call flux_down_from_atmos(Time_atmos, Atm, Land, Ice, Land_ice_atmos_boundary, &
                              Atmos_land_boundary, Atmos_ice_boundary )
    call fms_mpp_clock_end(coupler_clocks%flux_down_from_atmos)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('flux_down_from_atmos+', current_timestep)

  end subroutine coupler_flux_down_from_atmos

  !> @parblock
  !! Subroutine coupler_update_land_model_fast calls update_land_model_fast from land_model_mod
  !! to advance fast land processes by one atmospheric timestep. Clocks are initialized to
  !! measure runtime, pelist is set and synchronized before and after fast land model update, and
  !! checksums and memory usages are computed if do_chksum and do_debug are
  !! true respectively.
  !! @endparblock
  subroutine coupler_update_land_model_fast(Land, Atmos_land_boundary, atm_pelist, current_timestep, &
                                            coupler_chksum_obj, coupler_clocks)

    implicit none
    type(land_data_type), intent(inout) :: Land
      !< is the land model derived type
    type(atmos_land_boundary_type), intent(inout) :: Atmos_land_boundary
      !< is the atmosphere-to-land boundary derived type containing atm to land fluxes
    integer, dimension(:), intent(in) :: atm_pelis
      !< is the atmosphere PE list used to reset the current PE list after the land update
    integer, intent(in) :: current_timestep
      !< is the current coupled timestep index used for checksum labelling
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj
      !< is the coupler checksum object used to compute checksums
    type(coupler_clock_type),  intent(inout) :: coupler_clocks
      !< are the coupler clocks used to measure runtime of the land update

    call fms_mpp_clock_begin(coupler_clocks%update_land_model_fast) !< current pelist=Atm%pelis
    if (land_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(Land%pelist)

    call update_land_model_fast( Atmos_land_boundary, Land )

    if (land_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(atm_pelist)
    call fms_mpp_clock_end(coupler_clocks%update_land_model_fast)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('update_land_fast+', current_timestep)
    if (do_debug)  call fms_memutils_print_memuse_stats( 'update land')

  end subroutine coupler_update_land_model_fas

  !> @parblock
  !! Subroutine coupler_update_ice_model_fast calls update_ice_model_fast from ice_model_mod
  !! to advance fast sea ice by one atmospheric timestep.  Pelists are set and synchronized
  !! before and after the fast ice model update.  Runtime is measured by update_ice_model_fast.
  !! Checksums and memory usage reporting are controlled by do_chksum and do_debug.
  !! @endparblock
  subroutine coupler_update_ice_model_fast(Ice, Atmos_ice_boundary, atm_pelist, current_timestep, &
                                           coupler_chksum_obj, coupler_clocks)

    implicit none
    type(ice_data_type), intent(inout) :: Ice
      !< is the ice model derived type
    type(Atmos_ice_boundary_type), intent(inout) :: Atmos_ice_boundary
      !< is the atmosphere-to-ice boundary derived type containing fluxes passed down from the atmosphere
    integer, dimension(:), intent(in) :: atm_pelis
      !< is the atmosphere PE list used to reset the current PE list after the fast ice update
    integer, intent(in) :: current_timestep
      !< is the current coupled timestep index used for checksum labelling
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj
      !< is the coupler checksum object used to compute checksums
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< are the coupler clocks used to measure runtime of the fast ice update

    call fms_mpp_clock_begin(coupler_clocks%update_ice_model_fast)  !< current pelist = Atm%pelis
    if (ice_npes .NE. atmos_npes)call fms_mpp_set_current_pelist(Ice%fast_pelist)

    call update_ice_model_fast( Atmos_ice_boundary, Ice )

    if (ice_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(atm_pelist)
    call fms_mpp_clock_end(coupler_clocks%update_ice_model_fast)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('update_ice_fast+', current_timestep)
    if (do_debug)  call fms_memutils_print_memuse_stats( 'update ice')

  end subroutine coupler_update_ice_model_fas

  !> @parblock
  !! Subroutine coupler_flux_up_to_atmos calls flux_up_to_atmos
  !! to transfer updated surface states from land and ice to atmosphere.
  !! Runtime is measured by coupler_clocks%flux_up_to_atmos, and checksums are computed
  !! if do_chksum=.true.
  !! @endparblock
  subroutine coupler_flux_up_to_atmos(Land, Ice, Land_ice_atmos_boundary, Atmos_land_boundary, Atmos_ice_boundary,&
                                      Time_atmos, current_timestep, coupler_chksum_obj, coupler_clocks)

    implicit none
    type(land_data_type), intent(inout) :: Land
      !< is the land model derived type
    type(ice_data_type),  intent(inout) :: Ice
      !< is the ice model derived type
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_ice_atmos_boundary
      !< is the land-ice-to-atmosphere boundary derived type accumulating surface fluxes to the atmosphere
    type(atmos_land_boundary_type), intent(inout) :: Atmos_land_boundary
      !< is the atmosphere-to-land boundary derived type used to get dimensions
    type(atmos_ice_boundary_type), intent(inout) :: Atmos_ice_boundary
      !< is the atmosphere-to-ice boundary derived type used to get dimensions
    type(FmsTime_type), intent(in) :: Time_atmos
      !< is the current atmospheric model time in seconds
    integer, intent(in) :: current_timestep
      !< is the current coupled timestep index used for checksum labelling
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj
      !< is the coupler checksum object used to compute checksums
    type(coupler_clock_type), intent(in) :: coupler_clocks
      !< are the coupler clocks used to measure runtime

    call fms_mpp_clock_begin(coupler_clocks%flux_up_to_atmos)
    call flux_up_to_atmos(Time_atmos, Land, Ice, Land_ice_atmos_boundary, Atmos_land_boundary, Atmos_ice_boundary)
    call fms_mpp_clock_end(coupler_clocks%flux_up_to_atmos)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('flux_up2atmos+', current_timestep)

  end subroutine coupler_flux_up_to_atmos

  !> @parblock
  !! Subroutine coupler_update_atmos_model_up calls update_atmos_model_up from atmos_driver
  !! to finish the upward sweep of the tridiagonal eliminiation for heat/moisture and to compute
  !! the convective and large-scale tendencies.  Runtime is measured by coupler_clocks%update_atmos_model_up.
  !! Checksums are computed when do_chksum=.true., and memory usage is printed when do_debug=.true.
  !! @endparblock
  subroutine coupler_update_atmos_model_up(Atm, Land_ice_atmos_boundary, current_timestep, &
                                           coupler_chksum_obj, coupler_clocks)

    implicit none
    type(atmos_data_type),  intent(inout) :: Atm
      !< is the atmosphere model derived type
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_ice_atmos_boundary
      !< is the land-ice-to-atmosphere boundary derived type containing surface fluxes returned to the atmosphere
    integer, intent(in) :: current_timestep
      !< is the current coupled timestep index used for checksum labelling
    type(coupler_chksum_type),intent(in) :: coupler_chksum_obj
      !< is the coupler checksum object used to compute checksums
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< are the coupler clocks used to measure runtime of the upward atmospheric physics sweep

    call fms_mpp_clock_begin(coupler_clocks%update_atmos_model_up)
    call update_atmos_model_up(Land_ice_atmos_boundary, Atm)
    call fms_mpp_clock_end(coupler_clocks%update_atmos_model_up)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('update_atmos_up+', current_timestep)
    if (do_debug) call fms_memutils_print_memuse_stats( 'update up')

  end subroutine coupler_update_atmos_model_up

  !> @parblock
  !! Subroutine coupler_flux_atmos_to_ocean calls flux_atmos_to_ocean to compute
  !! atmosphere-to-ocean/ice gas deposition fluxes.
  !! @endparblock
  subroutine coupler_flux_atmos_to_ocean(Atm, Atmos_ice_boundary, Ice, Time_atmos)

    implicit none
    type(atmos_data_type), intent(inout) :: Atm
      !< is the atmosphere model derived type
    type(atmos_ice_boundary_type), intent(inout) :: Atmos_ice_boundary
      !< is the atmosphere-to-ice boundary derived type used to pass gas and deposition fluxes to the ocean
    type(ice_data_type), intent(inout) :: Ice
      !< is the ice model derived type
    type(FmsTime_type),  intent(in)    :: Time_atmos
      !< is the current atmospheric model time in seconds

    call flux_atmos_to_ocean(Time_atmos, Atm, Atmos_ice_boundary, Ice)
    call flux_ex_arrays_dealloc()

  end subroutine coupler_flux_atmos_to_ocean

  !> @parblock
  !! Subroutine coupler_update_atmos_model_state calls update_atmos_model_state in atmos_model
  !! to update atmospheric state and diagnostic fields in Atm at the end of the atmospheric timestep.
  !! Runtime is measured by coupler_clocks%update_atmos_model_state. Checksums are computed when
  !! do_chksum=.true., and memory usage is printed when do_debug=.true.
  !! @endparblock
  subroutine coupler_update_atmos_model_state(Atm, current_timestep, coupler_chksum_obj, coupler_clocks)

    implicit none
    type(atmos_data_type), intent(inout)  :: Atm
      !< is the atmosphere model derived type
    integer, intent(in)     :: current_timestep
      !< is the current coupled timestep index used for checksum labelling
    type(coupler_chksum_type), intent(in)    :: coupler_chksum_obj
      !< is the coupler checksum object used to compute checksums
    type(coupler_clock_type),  intent(inout) :: coupler_clocks
      !< are the coupler clocks used to measure runtime of the atmospheric state update

    call fms_mpp_clock_begin(coupler_clocks%update_atmos_model_state)
    call update_atmos_model_state( Atm )
    call fms_mpp_clock_end(coupler_clocks%update_atmos_model_state)

    if (do_chksum) &
        call coupler_chksum_obj%get_atmos_ice_land_chksums('update_atmos_model_state+', current_timestep)
    if (do_debug)  call fms_memutils_print_memuse_stats( 'update state')

  end subroutine coupler_update_atmos_model_state

  !> @parblock
  !! Subroutine coupler_update_land_model_slow calls update_land_model_slow in land_model_mod
  !! to advance the land model on the slow (coupled) timestep.  Pelist are set and synchronized
  !! before and after the call.  Runtime is measured with coupler_clocks%update_land_model_slow.
  !! Checksums are computed when do_chksum=.true.
  !! @endparblock
  subroutine coupler_update_land_model_slow(Land, Atmos_land_boundary, atm_pelist, current_timestep, &
                                            coupler_chksum_obj, coupler_clocks)

    implicit none
    type(land_data_type), intent(inout) :: Land
      !< is the land model derived type
    type(atmos_land_boundary_type), intent(inout) :: Atmos_land_boundary
      !< is the atmosphere-to-land boundary derived type containing fluxes passed down from the atmosphere
    integer, dimension(:), intent(in) :: atm_pelis
      !< is the atmosphere PE list used to reset the current PE list after the slow land update
    integer, intent(in) :: current_timestep
      !< is the current coupled timestep index used for checksum labelling
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj
      !< is the coupler checksum object used to compute checksums
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< are the coupler clocks used to measure runtime of the slow land update

    call fms_mpp_clock_begin(coupler_clocks%update_land_model_slow)

    if (Land%pe) then
      if (land_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(Land%pelist)
      call update_land_model_slow(Atmos_land_boundary,Land)
    endif

    if (land_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(atm_pelist)
    call fms_mpp_clock_end(coupler_clocks%update_land_model_slow)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('update_land_slow+', current_timestep)

  end subroutine coupler_update_land_model_slow

  !> @parblock
  !! Subroutine coupler_flux_land_to_ice calls flux_land_to_ice to transfer
  !! freshwater discharge from the land model to the ice/ocean grid.
  !! Runtime is measured with coupler_clocks%flux_land_to_ice, and checksums are
  !! computed when do_chksum=.true.
  !! @endparblock
  subroutine coupler_flux_land_to_ice(Land, Ice, Land_ice_boundary, Time, current_timestep, &
                                      coupler_chksum_obj, coupler_clocks)

    implicit none
    type(land_data_type), intent(inout) :: Land
      !< is the land model derived type
    type(ice_data_type),  intent(inout) :: Ice
      !< is the ice model derived type
    type(land_ice_boundary_type), intent(inout) :: Land_ice_boundary
      !< is the land-to-ice boundary derived type receiving runoff and other land fluxes
    type(FmsTime_type), intent(in) :: Time
      !< is the current model time in seconds passed to flux_land_to_ice
    integer, intent(in) :: current_timestep
      !< is the current coupled timestep index used for checksum labelling
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj
      !< is the coupler checksum object used to compute checksums
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< are the coupler clocks used to measure runtime of the land-to-ice flux transfer

    call fms_mpp_clock_begin(coupler_clocks%flux_land_to_ice)
    call flux_land_to_ice( Time, Land, Ice, Land_ice_boundary )
    call fms_mpp_clock_end(coupler_clocks%flux_land_to_ice)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('fluxlnd2ice+', current_timestep)

  end subroutine coupler_flux_land_to_ice

  !> @parblock
  !! Subroutine coupler_unpack_land_ice_boundary calls ice_model_fast_cleanup
  !! and unpack_land_ice_boundary from ice_model_mod to prepare the fast-ice model to receive
  !! the new land discharge fields and then copies them into its internal state.
  !! Ice_model_fast_cleanup resets the fast-ice accumulation buffers so tha
  !! the incoming runoff/calving values replace, rather than accumulate on top of,
  !! values from previous steps.
  !! Unpack_land_ice_boundary(Ice, Land_ice_boundary copies runoff, calving, runoff_hflx, calving_hflx
  !! from Land_ice_boundary that was populated by coupler_flux_land_to_ice into the ice model's internal
  !! fast-ice derived types.
  !! @endparblock
  subroutine coupler_unpack_land_ice_boundary(Ice, Land_ice_boundary, coupler_clocks)

    implicit none
    type(ice_data_type), intent(inout) :: Ice
      !< is the ice model derived type
    type(land_ice_boundary_type), intent(inout) :: Land_ice_boundary
      !< is the land-to-ice boundary derived type whose fields are unpacked into the ice model
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< are the coupler clocks used to measure runtime of unpacking

    if (ice_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(Ice%fast_pelist)
    call fms_mpp_clock_begin(coupler_clocks%update_ice_model_slow_fast)

    ! These two calls occur on whichever PEs handle the fast ice processess.
    call ice_model_fast_cleanup(Ice)
    call unpack_land_ice_boundary(Ice, Land_ice_boundary)

    call fms_mpp_clock_end(coupler_clocks%update_ice_model_slow_fast)

  end subroutine coupler_unpack_land_ice_boundary

  !> @parblock
  !! Subroutine coupler_update_ice_model_slow_and_stocks calls update_ice_model_slow
  !! from ice_model_mod to advance the slow sea-ice model and then calls flux_ice_to_ocean_stocks
  !! to compute stocks.  Update_ice_model_slow runs slow-timescale sea-ice processes including
  !! dynamics, freezing and melting, precipitation, and transport.  Flux_ice_to_ocean_stocks
  !! updates for stocks transferred from ice to ocean.  Runtime is measured by
  !! coupler_clocks%update_ice_model_slow_slow and coupler_clocks%flux_ice_to_ocean_stocks (inner clock).
  !! @endparblock
  subroutine coupler_update_ice_model_slow_and_stocks(Ice, coupler_clocks)

    implicit none
    type(ice_data_type), intent(inout) :: Ice
      !< is the ice model derived type
    type(coupler_clock_type), intent(inout) :: coupler_clocks
      !< is the coupler timing clock set used to measure runtime of the slow ice update and stock flux steps

    if (slow_ice_with_ocean) call fms_mpp_set_current_pelist(Ice%slow_pelist)
    call fms_mpp_clock_begin(coupler_clocks%update_ice_model_slow_slow)

    call update_ice_model_slow(Ice)

    call fms_mpp_clock_begin(coupler_clocks%flux_ice_to_ocean_stocks)
    call flux_ice_to_ocean_stocks(Ice)
    call fms_mpp_clock_end(coupler_clocks%flux_ice_to_ocean_stocks)

    call fms_mpp_clock_end(coupler_clocks%update_ice_model_slow_slow)

  end subroutine coupler_update_ice_model_slow_and_stocks

  !> @parblock
  !! Subroutine coupler_update_ocean_model calls update_ocean_model from ocean_model_mod
  !! to advance the ocean model by one coupled timestep with ice-ocean boundary forcing.
  !! Time_ocean is advanced by Time_step_cpld, and checksums are computed when
  !! do_chksum=.true.
  !! @endparblock
  subroutine coupler_update_ocean_model(Ocean, Ocean_state, Ice_ocean_boundary, &
                                        Time_ocean, Time_step_cpld, current_timestep, coupler_chksum_obj)

    implicit none
    type(ocean_public_type), intent(inout) :: Ocean
      !< is the ocean model public derived type
    type(ocean_state_type), pointer, intent(inout) :: Ocean_state
      !< is the pointer to the internal ocean model state
    type(Ice_ocean_boundary_type),   intent(inout) :: Ice_ocean_boundary
      !< is the ice-to-ocean boundary derived type containing forcing fluxes passed to the ocean
    type(FmsTime_type), intent(inout) :: Time_ocean
      !< is the current ocean model time; advanced by Time_step_cpld on outpu
    type(FmsTime_type), intent(in) :: Time_step_cpld
      !< is the duration of one coupled (slow) timestep passed to update_ocean_model
    integer, intent(in) :: current_timestep
      !< is the current coupled timestep index used for checksum labelling
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj
      !< is the coupler checksum object used to compute and report field checksums

    call update_ocean_model(Ice_ocean_boundary, Ocean_state,  Ocean, Time_ocean, Time_step_cpld)
    if (do_chksum) call coupler_chksum_obj%get_ocean_chksums('update_ocean_model+', current_timestep)

  end subroutine coupler_update_ocean_model

  !> @parblock
  !! Subroutine coupler_intermediate_restart writes mid-run restart files for all
  !! component models and the coupler.
  !!
  !! Component restarts are written on their respective PE sets: atmosphere, land,
  !! and ice restarts are written on PEs where Atm%pe is true, and the ocean
  !! restart is written on PEs where Ocean%is_ocean_pe is true. Coupler-specific
  !! boundary-condition restart data (Ocn_bc_restart and Ice_bc_restart) are
  !! written by coupler_restart in FMS.
  !!
  !! After all files are written, Time_restart is advanced by restart_interval
  !! to set the next scheduled intermediate restart write.
  !! @endparblock
  subroutine coupler_intermediate_restart(Atm, Ice, Ocean, Ocean_state, Ocn_bc_restart, Ice_bc_restart,&
                                          Time_current, Time_restart, Time_restart_current, Time_start)

    implicit none
    type(atmos_data_type),   intent(inout) :: Atm
      !< is the atmosphere model derived type; restart is written by atmos_model_restar
    type(ice_data_type), intent(inout) :: Ice
      !< is the ice model derived type; restart is written by ice_model_restar
    type(ocean_public_type), intent(inout) :: Ocean
      !< is the ocean model public derived type; restart is written by ocean_model_restar
    type(ocean_state_type), pointer, intent(inout) :: Ocean_state
      !< is the pointer to the internal ocean model state passed to ocean_model_restar
    type(FmsNetcdfDomainFile_t), pointer, intent(inout) :: Ocn_bc_restart(:)
      !< is the array of fms2_io fileobjs used to write ocean boundary-condition coupler restart data
    type(FmsNetcdfDomainFile_t), pointer, intent(inout) :: Ice_bc_restart(:)
      !< is the array of fms2_io fileobjs used to write ice boundary-condition coupler restart data
    type(FmsTime_type), intent(in) :: Time_curren
      !< is the current model time stamped on the intermediate restart files
    type(FmsTime_type), intent(in) :: Time_star
      !< is the model start time passed to coupler_restar
    type(FmsTime_type), intent(inout) :: Time_restar
      !< is the next scheduled intermediate restart time; updated by this subroutine after writing restarts
    type(FmsTime_type), intent(inout) :: Time_restart_curren
      !< is the current intermediate restart time; set to Time_current at the start of this subroutine

    character(len=32) :: timestamp ! Time in string
    integer :: outunit             ! stdou

    Time_restart_current = Time_curren

    timestamp = fms_time_manager_date_to_string(Time_restart_current)
    outunit= fms_mpp_stdout()
    write(outunit,*) '=> NOTE from program coupler: intermediate restart file is written and ', &
                      trim(timestamp),' is appended as prefix to each restart file name'
    if (Atm%pe) then
      call atmos_model_restart(Atm, timestamp)
      call land_model_restart(timestamp)
      call ice_model_restart(Ice, timestamp)
    endif
    if (Ocean%is_ocean_pe) call ocean_model_restart(Ocean_state, timestamp)

    call coupler_restart(Atm, Ice, Ocean, Ocn_bc_restart, Ice_bc_restart, &
                         Time_current, Time_restart_current, Time_start, timestamp)

    Time_restart = fms_time_manager_increment_date(Time_current, restart_interval(1), restart_interval(2), &
                   restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )

  end subroutine coupler_intermediate_restar

  !> @parblock
  !! Subroutine coupler_summarize_timestep reports coupled-timestep progress,
  !! memory usage, and optional concurrent-radiation timing diagnostics.
  !! Checksums are computed when do_chksum=.true. and summary text is written to
  !! stdout each timestep.
  !! @endparblock
  subroutine coupler_summarize_timestep(current_timestep, num_cpld_calls, coupler_chksum_obj, &
                                        is_atmos_pe, omp_sec, imb_sec)

    implicit none
    integer, intent(in) :: current_timestep
      !< is the current coupled timestep index (nc) used for checksum labelling and progress reporting
    integer, intent(in) :: num_cpld_calls
      !< is the total number of coupled (outer-loop) timesteps in the run, used for progress reporting
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj
      !< is the coupler checksum object used to compute and report end-of-timestep field checksums
    logical, intent(in)  :: is_atmos_pe
      !< is Atm%pe; true if this PE belongs to the atmosphere PE list, required for concurrent-radiation timing outpu
    real, dimension(:), intent(inout) :: omp_sec
      !< is the elapsed wall-clock seconds for each concurrent OpenMP section (atmosphere, radiation)
    real, dimension(:), intent(inout) :: imb_sec
      !< is the OpenMP load-imbalance seconds for each concurrent OpenMP section

    integer :: outunit        ! stdou
    character(len=80) :: text ! text to be written out to stdou

    if (do_chksum) call coupler_chksum_obj%get_coupler_chksums('MAIN_LOOP+', current_timestep)
    write( text,'(a,i6)' )'Main loop at coupling timestep=', current_timestep
    call fms_memutils_print_memuse_stats(text)
    outunit= fms_mpp_stdout()

    if (fms_mpp_pe() == fms_mpp_root_pe() .and. is_atmos_pe .and. do_concurrent_radiation) &
        write(outunit,102) 'At coupling step ', current_timestep,' of ',num_cpld_calls, ' Atm & Rad (imbalance): ', &
                            omp_sec(1),' (',imb_sec(1),')  ',omp_sec(2),' (',imb_sec(2),')'

    call flush(outunit)

102 format(A17,i5,A4,i5,A24,f10.4,A2,f10.4,A3,f10.4,A2,f10.4,A1)

  end subroutine coupler_summarize_timestep

end module full_coupler_mod
