#Coupler_main - top-level program for the full FMSCoupler

## Introduction
Program coupler_main contains the main time loops to call the time-stepping 
dynamics to advance the coupled model. Coupler_main also calls flux exchange between 
atmosphere, ocean, land, and sea-ice.

## Time integration
There are two nested time integration loops:
  - Slow (coupled) loop (coupled_timestep_loop, nc = 1 … num_cpld_calls):
    advances the ocean and slow sea-ice by one coupled timestep dt_cpld.
    Ocean–ice fluxes are exchanged once per iteration.
  - Fast (atmospheric) loop (fast_integration_loop, na = 1 … num_atmos_calls):
    advances the atmosphere, land surface, and fast sea-ice by one atmospheric
    timestep dt_atmos.  The fast loop is nested inside the slow loop.

## Fast loop
Heat and moisture are exchanged between the atmosphere and the surface
 (land/ice) using an tridiagonal scheme for implicit vertical diffusion in the fast loop:
 1. coupler_update_atmos_model_down — forward (downward) sweep from atmospheric top to surface.
 2. coupler_flux_down_from_atmos — transfers forward-elimination coefficients to land/ice.
 3. Land and ice fast updates (coupler_update_land_model_fast,
    coupler_update_ice_model_fast) — compute new surface temperatures.
 4. coupler_flux_up_to_atmos — applies implicit surface flux corrections using
    updated surface temperatures.
 5. coupler_update_atmos_model_up — back-substitution (upward sweep), convection,
    and large-scale condensation.

 When `do_concurrent_radiation = .true.`, the atmospheric dynamics/physics and
 radiation updates run simultaneously with OpenMP threading (see `atmos_nthreads`
 and `radiation_nthreads` namelist variables).

## Slow loop
Ocean fluxes are exchanged explicitly (one coupling step behind).  All fluxes
reaching the ocean — including atmospheric fluxes and land runoff — are passed
through the sea-ice model via Ice_ocean_boundary.

The below section outlines what is not used:
Sea-ice physics is split into two timescales that can run on different MPI PE sets:
  - Fast ice on Ice%fast_ice_pe: thermodynamics and surface-flux coupling at the
    atmospheric timestep.  Fast ice always runs on a subset of the atmosphere PEs
  - Slow ice on Ice%slow_ice_pe: dynamics, freezing/melting, and transport at the
    coupled (ocean) timestep.  The placement of the slow ice PEs depends on
    the slow_ice_with_ocean namelist variable:
      - slow_ice_with_ocean = .false. (default): slow and fast ice share the same PEs
      - slow_ice_with_ocean = .true.: slow ice runs on the ocean PEs.  In this case,
        Ice%pelist is the union of the fast (atmos) and slow (ocean) PE sets.
The flag concurrent_ice = .true. runs the fast ice and slow ice processes concurrently
and requires slow_ice_with_ocean = .true.  
The flag combine_ice_and_ocean = .true. advances the slow ice and ocean processes together 
on the ocean PEs.  The flags concurent_ice and slow_ice_with_ocean must be .true. to use combine_ice_and_ocean.

## Pseudocode
do nc = 1, num_cpld_calls

  ! Redistribute ocean surface state onto ice grid
  call flux_ocean_to_ice(...)

  ! If slow_ice_pe: override ocean-ice BCs, unpack into Ice type
  call flux_ocean_to_ice_finish(...)
  call unpack_ocean_ice_boundary(...)

  ! Exchange slow-ice state to fast-ice data structures
  call exchange_slow_to_fast_ice(...)

  ! Prepare ice surface fields (albedo, T, etc.) for atmos surface flux calc
  call set_ice_surface_fields(...)

  do na = 1, num_atmos_calls

    !Copy Atm%tr_bot → Atm%fields for gas-exchange tracers
    call atmos_tracer_driver_gather_data(Atm%fields, Atm%tr_bot)

    ! Compute surface exchange coefficients and turbulent fluxes on the
    ! atm-land-ice exchange grid
    call sfc_boundary_layer(...)

    ! Atmosphere dynamical core (FV3)
    call update_atmos_model_dynamics(...)

    ! Radiation (sequential, or concurrent on a separate OMP team)
    call update_atmos_model_radiation(...)

    ! Forward (downward) sweep of the implicit tridiagonal diffusion
    call update_atmos_model_down(...)

    ! Apply implicit atm diffusion correction; pass updated surface fluxes
    ! (heat, moisture, momentum) to land and ice boundary types
    call flux_down_from_atmos(...)

    ! Fast land physics (hydrology, canopy, soil temperature)
    call update_land_model_fast(...)

    ! Fast ice thermodynamics (surface energy balance, melt ponds)
    call update_ice_model_fast(...)

    ! Recompute surface fluxes using updated land/ice surface temperatures
    call flux_up_to_atmos(...)

    ! Back-substitution (upward) sweep; convection; large-scale condensation
    call update_atmos_model_up(...)

    ! Remap atmosphere gas/tracer fields onto exchange grid;
    ! compute air-sea deposition fluxes; deallocate exchange grid arrays
    call flux_atmos_to_ocean(...)
    call flux_ex_arrays_dealloc()

    ! Advance atmos diagnostics and tracer state
    call update_atmos_model_state(...)

  end do  ! fast loop

  ! Slow land physics (routing, carbon, DGVM)
  call update_land_model_slow(...)

  ! Interpolate land runoff and calving onto ice grid
  call flux_land_to_ice(...)

  ! Reset fast-ice accumulators; copy Land_ice_boundary into ice internals
  call ice_model_fast_cleanup(...)
  call unpack_land_ice_boundary(...)

  ! Exchange fast-ice averages to slow-ice side
  call exchange_fast_to_slow_ice(...)

  ! Slow ice physics (dynamics, freezing/melting, transport)
  call update_ice_model_slow(...)

  if(Ocean%is_ocean_pe) then
    ! Interpolate ice-bottom fluxes onto ocean grid -> Ice_ocean_boundary
    call flux_ice_to_ocean(...)
    ! Override/diagnose Ice_ocean_boundary fields; send to diag_manager
    call flux_ice_to_ocean_finish(...)
    ! Advance ocean state by dt_cpld using Ice_ocean_boundary forcing
    call update_ocean_model(...)
  endif

  ! Bookkeep ocean stocks from ice-ocean flux transfer
  call flux_ocean_from_ice_stocks(...)

  ! Flush diagnostic send-data buffer for this coupled step
  call fms_diag_send_complete(Time_step_cpld)

end do

call coupler_restart(...) ! write coupler.res and component restart files
call fms_diag_end(...)    ! flush and close diagnostic output 

## MPI Parallelization
Users can specify the number of processing elements (MPI ranks here on abbreviated as 'pes') for each 
model component in the coupler namelist as shown below:
```
&coupler_nml
atmos_npes = 10
ocean_npes = 10
ice_npes = 10
land_npes = 10
\
```
The number of pes for each component must meet the following:
  * At least atmos_npes or ocean_npes must be specified.  
  * land_npes <= atmos_npes
  * ice_npes <= atmos_npes 
  * atmos_npes + ocean_npes = npes (total number of pes determined with FMS)

When concurrent = .true., concurrent_ice = .false, and slow_ice_with_ocean = .false.
  * atm and ocean will have distinct set of pelists
  * land%pelist will be a subset of atm%pelist
  * ice%pelist = ice%slow_pelist = ice%fast_pelist = subset of atm%pelist

## OpenMP Parallelization
Users can also specify the number of OpenMP threads as below:
```
&coupler_nml
do_concurrent_radiaton = .true.
use_hyper_thread = .true.
conc_nthreads = 2
atmos_nthreads = 1
radiation_nthreads = 1
ocean_nthreads = 1
\
```
Note, the model must be compiled with OpenMP enabled (this can be achieved 
with fre by specifying targets with "-openmp" such as "prod-openmp")

When do_concurrent_radiation is true, conc_nthreads will be set to 2:
thread 0 on the atm%pes will run atmosphere dynamics and physics while thread 1
will run the radiation dynamics.  Else, radiation will run sequentuaally after 
the atmosphere update.  Note, atmos_nthreads will affect the number of threads 
within the atmosphere dynamics.

## Model Component State Types
The following are derived types holding data for each component 
  - Atm (atmos_data_type):  Instantaneous atm model state at current timestep
  - Land (land_data_type):  Instantaneous land model state at current timestep
  - Ice (ice_data_type):  Instantaneous ice model state at current timestep
  - Ocean (ocean_public_type): Ocean model state at current timestep; target for the Ocean_state pointer
  - Ocean_state (ocean_state_type pointer): Alias pointer to Ocean datatype

## Boundary Exchange Types
  - Atmos_land_boundary (atmos_land_boundary_type):  Data to exchange between atmos and land
  - Atmos_ice_boundary (atmos_ice_boundary_type):  Data to exchange between atmos and sea ice
  - Land_ice_atmos_boundary (land_ice_atmos_boundary_type):  Data to exchange between land and ice to atmos
  - Land_ice_boundary (land_ice_boundary_type):  Data to exchange between land and ice
  - Ice_ocean_boundary (ice_ocean_boundary_type):  Data to exchange from ice and ocean
  - Ocean_ice_boundary (ocean_ice_boundary_type):  Data to exchange from ocean to ice
  - ice_ocean_driver_CS (ice_ocean_driver_type (pointer)):  Pointer alias containing control
    parameters for the combined ice-ocean driver

## Time Variables
  - Time_step_atmos (FmsTime_type):  Timestep in the fast loop
  - Time_step_cpld (FmsTime_type):  Timestep in the slow loop
  - Time_atmos (FmsTime_type):  Time tracked in the fast loop
  - Time_ocean (FmsTime_type):  Time tracked for the ocean model
  - Time_flux_ice_to_ocean (FmsTime_type):  Time tracked for lag fluxes from ice to ocean
  - Time_flux_ocean_to_ice (FmsTime_type):  Time tracked for flux exchange from ocean to ice
  - Time_restart (FmsTime_type):  Next timepoint to write intermediate restarts
  - Time_restart_current (FmsTime_type):  Last timepoint when intermediate restarts were written
  - Time_start (FmsTime_type):  Model start time
  - Time_end (FmsTime_type):  Model end time

## Loop Counters
  - num_atmos_calls (integer):  Number of timesteps in the fast-integration loop
  - na (integer):  Do-loop counter in the fast-integration loop
  - num_cpld_calls (integer):  Number of timesteps in the slow-integration loop
  - nc (integer):  Do-loop counter in the slow-integration loop
  - current_timestep (integer):  Total number of fast loop iteration, equal to 
    (nc-1)*num_atmos_calls + na.  
  

