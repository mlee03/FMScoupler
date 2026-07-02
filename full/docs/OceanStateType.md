# `ocean_state_type` — MOM6 Interior Ocean State

## Overview

`ocean_state_type` contains information about the full interior state of the MOM6 ocean model. All members are **private** and are not accessible outside `ocean_model_MOM.F90`. The coupler holds only a pointer to it (`type(ocean_state_type), pointer :: Ocean_state`); the public surface fields are exposed separately through `ocean_public_type`.

**Defined in:** `ocean_model_MOM.F90`  
**Related types:** `ocean_public_type`  
**Key subroutines:** `update_ocean_model`, `ocean_model_init`, `ocean_model_end`

---

## PE and Time Fields

| Field | Type | Description |
|---|---|---|
| `Ocean_state%is_ocean_PE` | logical | `.true.` if this PE is part of the ocean pelist; initialised to `.false.` so non-ocean PEs that hold the pointer do nothing. |
| `Ocean_state%Time` | `type(time_type)` | The ocean model's master clock; advanced each call to `update_ocean_model`. |
| `Ocean_state%Time_dyn` | `type(time_type)` | The ocean model's dynamics time; equals `Time` after a complete timestep but may lag during a split dynamics/thermodynamics step. |

---

## Timestep and Step-Count Control

| Field | Type | Default | Description |
|---|---|---|---|
| `Ocean_state%nstep` | integer | 0 | Counter for the number of `update_ocean` calls that advanced the dynamics. |
| `Ocean_state%nstep_thermo` | integer | 0 | Counter for the number of `update_ocean` calls that advanced the thermodynamics. |
| `Ocean_state%single_step_call` | logical | `.true.` | If `.true.`, dynamics and thermodynamics are advanced together in a single call. If `.false.`, separate calls are made using `dt` and `dt_therm`. |
| `Ocean_state%dt` | real | — | Baroclinic dynamics timestep [T ~> s]; only used when `single_step_call=.false.`. |
| `Ocean_state%dt_therm` | real | — | Thermodynamics timestep [T ~> s]; only used when `single_step_call=.false.`. |
| `Ocean_state%thermo_spans_coupling` | logical | — | If `.true.`, thermodynamic and tracer timesteps can span multiple coupled timesteps. |
| `Ocean_state%diabatic_first` | logical | — | If `.true.`, diabatic and thermodynamic processes are applied before the dynamics step. |
| `Ocean_state%Restart_control` | integer | — | Bit-field controlling restart file writing: bit 0 (+1) saves generic restart files; bit 1 (+2) saves time-stamped files. A negative value suppresses restart writing at run end. |

---

## Feature Flags

| Field | Type | Description |
|---|---|---|
| `Ocean_state%use_ice_shelf` | logical | If `.true.`, the MOM6 ice shelf model (`MOM_ice_shelf`) is enabled. |
| `Ocean_state%use_waves` | logical | If `.true.`, surface wave coupling is active. |
| `Ocean_state%icebergs_alter_ocean` | logical | If `.true.`, icebergs can modify ocean dynamics and forcing fluxes. |
| `Ocean_state%calve_ice_shelf_bergs` | logical | If `.true.`, icebergs are seeded from ice-shelf flux through the ice front. |
| `Ocean_state%offline_tracer_mode` | logical | If `.false.` (default), full prognostic dynamics and thermodynamics are integrated. If `.true.`, only tracer advection/diffusion is integrated using velocity fields read from a previous run. |

---

## Physical Constants

| Field | Type | Units | Description |
|---|---|---|---|
| `Ocean_state%C_p` | real | J degC⁻¹ kg⁻¹ | Specific heat capacity of seawater. |
| `Ocean_state%press_to_z` | real | Z T² R⁻¹ L⁻² ~> m Pa⁻¹ | Conversion factor from pressure to ocean depth, typically `1/(ρ₀g)`. |

---

## Forcing Structures

| Field | Type | Description |
|---|---|---|
| `Ocean_state%forces` | `type(mech_forcing)` | Mechanical surface forcing: wind stress, sea-level pressure, and wave-related fields. |
| `Ocean_state%fluxes` | `type(forcing)` | Primary thermodynamic ocean forcing: heat flux, freshwater flux, salt flux, shortwave penetration, etc. |
| `Ocean_state%flux_tmp` | `type(forcing)` | Secondary forcing structure used when multiple coupled timesteps are taken per thermodynamic step; accumulates forcing between thermodynamic updates. |
| `Ocean_state%sfc_state` | `type(surface)` | Ocean surface state fields (SST, SSS, surface currents, boundary layer depth) returned to the coupler after each update. |

---

## Grid Structures

| Field | Type | Description |
|---|---|---|
| `Ocean_state%grid` | `type(ocean_grid_type)` (pointer) | MOM6 horizontal grid structure: cell areas, distances, coordinates, metric terms, and land masks. |
| `Ocean_state%GV` | `type(verticalGrid_type)` (pointer) | MOM6 vertical grid structure: layer thicknesses, target densities (for isopycnal coordinates), and ALE remapping parameters. |
| `Ocean_state%US` | `type(unit_scale_type)` (pointer) | MOM6 dimensional unit-scaling factors; converts between external MKS units and MOM6's internal non-dimensionalised units. |

---

## MOM6 Control Structures

| Field | Type | Description |
|---|---|---|
| `Ocean_state%MOM_CSp` | `type(MOM_control_struct)` | MOM6 master control structure; holds all module-level control structures, parameter settings, and diagnostic handles for the ocean model. |
| `Ocean_state%Ice_shelf_CSp` | `type(ice_shelf_CS)` (pointer) | Control structure for the MOM6 ice shelf model; null if `use_ice_shelf=.false.`. |
| `Ocean_state%marine_ice_CSp` | `type(marine_ice_CS)` (pointer) | Control structure for the marine ice effects module (e.g., iceberg melt parameterisation). |
| `Ocean_state%Waves` | `type(wave_parameters_cs)` (pointer) | Control structure for surface wave coupling; null if `use_waves=.false.`. |
| `Ocean_state%forcing_CSp` | `type(surface_forcing_CS)` (pointer) | MOM6 surface forcing control structure; handles the translation from coupler boundary conditions to internal MOM6 forcing arrays. |
| `Ocean_state%diag` | `type(diag_ctrl)` (pointer) | MOM6 diagnostic control structure; manages registration and posting of all ocean diagnostics to the FMS `diag_manager`. |
| `Ocean_state%dirs` | `type(directories)` | Structure containing relevant directory paths for input, output, and restart files. |
