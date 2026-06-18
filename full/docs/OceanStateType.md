# ocean_state_type
ocean_state_type contains information about the state of the ocean.
Its members are private and are not accessible outside ocean_model_MOM.F90.  
The coupler holds only a pointer to it (type(ocean_state_type), pointer :: Ocean_state).

## Ocean_state%is_ocean_PE
Ocean_state%is_ocean_PE, logical, is .true. if this PE is part of the ocean pelist; initialised to .false. so non-ocean PEs that hold the pointer do nothing.

## Ocean_state%Time
Ocean_state%Time, type(time_type), is The ocean model's master clock; advanced each call to update_ocean_model.
## Ocean_state%Time_dyn
Ocean_state%Time_dyn, type(time_type), is The ocean model's dynamics time; equals Time after a complete timestep but may lag during a split dynamics/thermodynamics step.

## Ocean_state%nstep
Ocean_state%nstep, integer, is Counter for the number of update_ocean calls that advanced the dynamics.
## Ocean_state%nstep_thermo
Ocean_state%nstep_thermo, integer, is Counter for the number of update_ocean calls that advanced the thermodynamics.
## Ocean_state%single_step_call
Ocean_state%single_step_call, logical, is If .true. (default), dynamics and thermodynamics are advanced together in a single call. If .false., separate calls are made and dt/dt_therm below are used.
## Ocean_state%dt
Ocean_state%dt, real, is Baroclinic dynamics timestep [T ~> s]; only used when single_step_call=.false..
## Ocean_state%dt_therm
Ocean_state%dt_therm, real, is Thermodynamics timestep [T ~> s]; only used when single_step_call=.false..
## Ocean_state%thermo_spans_coupling
Ocean_state%thermo_spans_coupling, logical, is If .true., thermodynamic and tracer timesteps can span multiple coupled timesteps.
## Ocean_state%diabatic_first
Ocean_state%diabatic_first, logical, is If .true., diabatic and thermodynamic processes are applied before the dynamics step.

## Ocean_state%Restart_control
Ocean_state%Restart_control, integer, is Bit-field controlling restart file writing: bit 0 (+1) saves generic restart files; bit 1 (+2) saves time-stamped files. A negative value suppresses restart writing at run end.

## Ocean_state%use_ice_shelf
Ocean_state%use_ice_shelf, logical, is If .true., the MOM6 ice shelf model (MOM_ice_shelf) is enabled.
## Ocean_state%use_waves
Ocean_state%use_waves, logical, is If .true., surface wave coupling is active.
## Ocean_state%icebergs_alter_ocean
Ocean_state%icebergs_alter_ocean, logical, is If .true., icebergs can modify ocean dynamics and forcing fluxes.
## Ocean_state%calve_ice_shelf_bergs
Ocean_state%calve_ice_shelf_bergs, logical, is If .true., icebergs are seeded from ice-shelf flux through the ice front.
## Ocean_state%offline_tracer_mode
Ocean_state%offline_tracer_mode, logical, is If .false. (default), full prognostic dynamics and thermodynamics are integrated. If .true., only tracer advection/diffusion is integrated using velocity fields read from a previous run.

## Ocean_state%C_p
Ocean_state%C_p, real, is Specific heat capacity of seawater [J degC⁻¹ kg⁻¹].
## Ocean_state%press_to_z
Ocean_state%press_to_z, real, is Conversion factor from pressure to ocean depth, typically 1/(ρ₀g) [Z T² R⁻¹ L⁻² ~> m Pa⁻¹].

## Ocean_state%forces
Ocean_state%forces, type(mech_forcing), is Mechanical surface forcing: wind stress, sea-level pressure, and wave-related fields.
## Ocean_state%fluxes
Ocean_state%fluxes, type(forcing), is Primary thermodynamic ocean forcing: heat flux, freshwater flux, salt flux, shortwave penetration, etc..
## Ocean_state%flux_tmp
Ocean_state%flux_tmp, type(forcing), is Secondary forcing structure used when multiple coupled timesteps are taken per thermodynamic step; accumulates forcing between thermodynamic updates.
## Ocean_state%sfc_state
Ocean_state%sfc_state, type(surface), is Ocean surface state fields (SST, SSS, surface currents, boundary layer depth) returned to the coupler after each update.

## Ocean_state%grid
Ocean_state%grid, type(ocean_grid_type), pointer, is MOM6 horizontal grid structure: cell areas, distances, coordinates, metric terms, and land masks.
## Ocean_state%GV
Ocean_state%GV, type(verticalGrid_type), pointer, is MOM6 vertical grid structure: layer thicknesses, target densities (for isopycnal coordinates), and ALE remapping parameters.
## Ocean_state%US
Ocean_state%US, type(unit_scale_type), pointer, is MOM6 dimensional unit-scaling factors used to convert between external MKS units and MOM6's internal non-dimensionalised units.

## Ocean_state%MOM_CSp
Ocean_state%MOM_CSp, type(MOM_control_struct), is MOM6 master control structure; holds all module-level control structures, parameter settings, and diagnostic handles for the ocean model.
## Ocean_state%Ice_shelf_CSp
Ocean_state%Ice_shelf_CSp, type(ice_shelf_CS), pointer, is Control structure for the MOM6 ice shelf model; null if use_ice_shelf=.false..
## Ocean_state%marine_ice_CSp
Ocean_state%marine_ice_CSp, type(marine_ice_CS), pointer, is Control structure for the marine ice effects module (e.g. iceberg melt parameterisation).
## Ocean_state%Waves
Ocean_state%Waves, type(wave_parameters_cs), pointer, is Control structure for surface wave coupling; null if use_waves=.false..
## Ocean_state%forcing_CSp
Ocean_state%forcing_CSp, type(surface_forcing_CS), pointer, is MOM6 surface forcing control structure; handles the translation from coupler boundary conditions to internal MOM6 forcing arrays.
## Ocean_state%diag
Ocean_state%diag, type(diag_ctrl), pointer, is MOM6 diagnostic control structure; manages registration and posting of all ocean diagnostics to the FMS diag_manager.

## Ocean_state%dirs
Ocean_state%dirs, type(directories), is Structure containing relevant directory paths for input, output, and restart files.
