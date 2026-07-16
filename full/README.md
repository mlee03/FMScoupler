# Coupler_main — FMSCoupler Full Coupler Top-Level Program

## Overview

`coupler_main` (defined in `coupler_main.F90`) is the top-level program for the **FMSCoupler full coupler**, which drives the coupled GFDL Earth-system model. It orchestrates time-stepping and flux exchange between four model components: **atmosphere**, **ocean**, **land**, and **sea ice**. This document describes the data types, time-integration loops, parallelization strategies, and namelist configuration options used by the full coupler.

---

## Model Component State Types

The full coupler uses the following Fortran derived types to hold the instantaneous state for each component. Each type is declared in `coupler_main.F90`.

| Variable Name | Type | Description |
|---|---|---|
| `Atm` | `atmos_data_type` | Holds the atmosphere model state. |
| `Land` | `land_data_type` | Holds the land model state. |
| `Ice` | `ice_data_type` | Holds the sea-ice model state. |
| `Ocean` | `ocean_public_type` | Contains public fields for the ocean model; is the target of the `Ocean_state` pointer. |
| `Ocean_state` | `ocean_state_type` (pointer) | Points to the full Ocean private state. |

---

## Boundary Exchange Types

The full coupler uses the following Fortran derived types to hold the fields exchanged at each component interface. Each type is declared in `coupler_main.F90`.

| Variable Name | Type | Description |
|---|---|---|
| `Atmos_land_boundary` | `atmos_land_boundary_type` | Fields exchanged between atmosphere and land. |
| `Atmos_ice_boundary` | `atmos_ice_boundary_type` | Fields exchanged between atmosphere and sea ice. |
| `Land_Ice_Atmos_Boundary` | `land_ice_atmos_boundary_type` | Aggregated surface state returned from land and ice to the atmosphere. |
| `Land_ice_boundary` | `land_ice_boundary_type` | Runoff and calving fields passed from land to ice. |
| `Ice_ocean_boundary` | `ice_ocean_boundary_type` | Fluxes passed from ice to the ocean. |
| `Ocean_ice_boundary` | `ocean_ice_boundary_type` | Ocean surface state passed to the ice model. |
| `Ice_ocean_driver_CS` | `ice_ocean_driver_type` (pointer) | Control parameters for the combined ice–ocean driver. |

---

## Time Integration

The full coupler has two nested time-integration loops:

- **Slow (coupled) loop**: Advances the ocean by one coupled timestep `dt_cpld`.
- **Fast (atmospheric) loop**: Advances the atmosphere, land surface, and fast sea ice by one atmospheric timestep `dt_atmos`.

The fast loop is nested inside the slow loop.

### Fast Loop — Atmosphere/Land/Ice Implicit Coupling

In the fast loop, heat and moisture are exchanged between the atmosphere and the surface (land/ice) using a **tridiagonal scheme for implicit vertical diffusion**. The calls execute in the following order:

1. **`coupler_update_atmos_model_down`** — Forward (downward) sweep from the atmospheric top to the surface.
2. **`coupler_update_atmos_model_radiation`** — Calls the radiation driver.
3. **`coupler_flux_down_from_atmos`** — Transfers forward-elimination coefficients to the land and ice models.
4. **`coupler_update_land_model_fast`** — Updates land hydrology, canopy physics, and surface temperature.
5. **`coupler_update_ice_model_fast`** — Updates fast sea-ice thermodynamics and surface temperature.
6. **`coupler_flux_up_to_atmos`** — Applies implicit surface-flux corrections using the updated surface temperatures.
7. **`coupler_update_atmos_model_up`** — Back-substitution (upward sweep); updates convection and large-scale condensation.

> **Concurrent radiation option**: When `do_concurrent_radiation = .true.`, radiation is updated concurrently with the remaining atmosphere, fast ice, and fast land updates.

### Slow Loop — Ocean Coupling

In the slow loop, when `concurrent = .true.`, the ocean model is updated **concurrently with the atmosphere** and is one coupling step behind.

- **MOM6 (recommended)**: Set `use_lag_fluxes = .false.` where `flux_ice_to_ocean` is called **after** `update_ocean_model`.
- **Older ocean models**: Set `use_lag_fluxes = .true.` for numerical stability — `flux_ice_to_ocean` is called **before** the ocean model update.

---

## Sea-Ice Physics: Fast and Slow Timescales

In the full coupler, sea-ice physics can be split into two timescales that run on different MPI PE (processing element) sets:

| Timescale | Processes | PE Assignment |
|---|---|---|
| **Fast ice** | Thermodynamics and surface-flux coupling at the atmospheric timestep | Always runs on a subset of the atmosphere PEs (`Ice%fast_ice_pe`) |
| **Slow ice** | Dynamics, freezing/melting, and transport at the coupled (ocean) timestep | Controlled by `slow_ice_with_ocean` namelist variable (`Ice%slow_ice_pe`) |

**`slow_ice_with_ocean` behavior:**

- `slow_ice_with_ocean = .false.` *(default)*: Slow and fast ice share the same PEs.
- `slow_ice_with_ocean = .true.`: Slow ice runs on the ocean PEs. `Ice%pelist` becomes the union of fast (atmosphere) and slow (ocean) PE sets.

**Uncommonly used options:**

- `concurrent_ice = .true.`: Fast and slow-ice processes run concurrently. Requires `slow_ice_with_ocean = .true.`.
- `combine_ice_and_ocean = .true.`: Slow ice and ocean are advanced together on the ocean PEs. Requires both `concurrent_ice = .true.` and `slow_ice_with_ocean = .true.`.

---

## Pseudocode — Main Science Call Sequence

The following pseudocode outlines the main science calls for the full coupler with these settings:

```
concurrent = .true.
use_lag_forces = .false.
concurrent_ice = .false.
do_concurrent_radiation = .false.
do_atm/land/ice/ocean = .true.
```

Only processes running on ocean PEs are explicitly shown. The rest runs on subsets or on atmosphere PEs.

```fortran
do nc = 1, num_cpld_calls

  ! Redistribute ocean surface states to Ocean_ice_boundary
  if (Ocean%is_ocean_pe) call flux_ocean_to_ice

  ! Map Ocean_ice_boundary fields to Ice state
  call unpack_ocean_ice_boundary

  ! Map slow-ice control structure (Ice%sCS) to fast-ice (Ice%fCS)
  call exchange_slow_to_fast_ice

  ! Prepare ice surface fields for atmosphere coupling
  call set_ice_surface_fields

  ! Generate surface exchange grid between land, ice, and atmosphere
  call generate_sfc_xgrid

  do na = 1, num_atmos_calls

    ! Copy atmosphere tracer bottom fields to atmosphere fields
    call atmos_tracer_driver_gather_data

    ! Compute surface exchange fluxes at the atmosphere–surface boundary
    call sfc_boundary_layer

    call update_atmos_model_dynamics

    call update_atmos_model_radiation

    ! Downward sweep of the implicit tridiagonal diffusion
    call update_atmos_model_down(...)

    ! Apply implicit atmosphere diffusion correction;
    ! pass updated surface fluxes (heat, moisture, momentum)
    ! to land and ice boundary types
    call flux_down_from_atmos

    call update_land_model_fast

    call update_ice_model_fast

    ! Recompute surface fluxes using updated land/ice surface temperatures
    call flux_up_to_atmos

    ! Back-substitution (upward) sweep
    call update_atmos_model_up

    ! Compute air–sea deposition fluxes
    call flux_atmos_to_ocean

    ! Update atmosphere and tracer state
    call update_atmos_model_state

  end do  ! fast (atmospheric) loop

  call update_land_model_slow(...)

  ! Interpolate land runoff, calving, and heat fluxes onto the ice grid
  call flux_land_to_ice(...)

  ! Reset fast-ice accumulators; copy Land_ice_boundary into Ice
  call ice_model_fast_cleanup(...)

  ! Exchange fast-ice averages to slow-ice
  call exchange_fast_to_slow_ice()

  call update_ice_model_slow()

  call flux_ice_to_ocean(...)
  if (Ocean%is_ocean_pe) call update_ocean_model(...)

end do  ! slow (coupled) loop
```

---

# MPI Parallelization

In the full coupler, the number of MPI processing elements (PEs) for each model component is specified in the `coupler_nml` namelist:

```fortran
&coupler_nml
  atmos_npes = 10
  ocean_npes = 10
  ice_npes   = 10
  land_npes  = 10
/
```

**PE layout rules:**

- atmosphere and ocean have distinct PE lists when `concurent = .true.`
- `land%pelist` is a subset of `atm%pelist` with `land_npes` less than or equal to `atmos_npes`
- `ice%pelist = ice%slow_pelist = ice%fast_pelist` are subsets of `atm%pelist` with `ice_npes` less than or equal to `atmos_npes`
- `atmos_npes` + `ocean_npes` = total number of PEs (`npes`)

---

# OpenMP Parallelization

OpenMP thread counts are configured in `coupler_nml`:

```fortran
&coupler_nml
  do_concurrent_radiation = .true.
  use_hyper_thread        = .true.
  conc_nthreads           = 2
  atmos_nthreads          = 1
  radiation_nthreads      = 1
  ocean_nthreads          = 1
/
```

* `do_concurrent_radiation`: When `.true.`, radiation runs concurrently with atmosphere dynamics/physics. Thread 0 runs atmosphere; thread 1 runs radiation. When `.false.`, radiation runs sequentially after the atmosphere update.
* `conc_nthreads`: Number of concurrent threads when `do_concurrent_radiation = .true.`; set to 2. 
* `atmos_nthreads`: Number of OpenMP threads used within atmosphere dynamics.
* `radiation_nthreads`: Number of OpenMP threads used for radiation calculations. 
* `ocean_nthreads`: Number of OpenMP threads used for ocean calculations. 
* `use_hyper_thread`: Enables use of hardware hyper-threading. 


The model must be compiled with OpenMP. When using the FRE build system, ensure the target name contains the `-openmp` suffix (e.g., `prod-openmp`).

**Thread behavior:**

| Namelist Variable | Description |
|---|---|

---

## Time Variables

The full coupler uses the following `FmsTime_type` variables to track model time:

| Variable | Description |
|---|---|
| `Time_step_atmos` | Timestep for the fast loop (atmospheric), `dt_atmos`. |
| `Time_step_cpld` | Timestep for the slow loop (coupled), `dt_cpld`. |
| `Time_atmos` | Current model time tracked in the fast loop. |
| `Time_ocean` | Current model time tracked for the ocean. |
| `Time_flux_ice_to_ocean` | Time when flux was last exchanged from ice to ocean. |
| `Time_flux_ocean_to_ice` | Time when flux was last exchanged from ocean to ice. |
| `Time_restart` | Next scheduled time for writing intermediate restarts. |
| `Time_restart_current` | Time at which the most recent intermediate restart was written. |
| `Time_start` | Model start time. |
| `Time_end` | Model end time. |

---

## Compiling a Model

Coming soon.

---

## Glossary

| Term | Definition |
|---|---|
| **FMSCoupler** | GFDL Flexible Modeling System Coupler; orchestrates component models in a coupled Earth-system simulation. |
| **Full coupler** | The `coupler_main` program in `coupler_main.F90`; drives the complete coupled model with atmosphere, land, ice, and ocean. |
| **MOM6** | Modular Ocean Model version 6; the ocean component. |
| **SIS2** | Sea Ice Simulator version 2; the sea-ice component. |
| **LM4** | Land Model version 4; the land surface component. |
| **dt_atmos** | Atmospheric (fast loop) timestep. |
| **dt_cpld** | Coupled (slow loop) timestep. |
| **PE** | MPI Processing Element; a parallel compute task. |
| **Tridiagonal scheme** | An implicit numerical method for vertical diffusion that solves a tridiagonal matrix system in a downward (forward elimination) and upward (back-substitution) sweep. |
| **concurrent** | Namelist flag; when `.true.`, the ocean runs simultaneously with the atmosphere on separate PEs, one coupling step behind. |
| **use_lag_fluxes** | Namelist flag; when `.true.`, `flux_ice_to_ocean` is called before `update_ocean_model` for numerical stability. |
| **slow_ice_with_ocean** | Namelist flag; when `.true.`, slow sea-ice dynamics run on the ocean PEs. |
| **concurrent_ice** | Namelist flag; when `.true.`, fast and slow ice run concurrently (requires `slow_ice_with_ocean = .true.`). |
| **combine_ice_and_ocean** | Namelist flag; when `.true.`, slow ice and ocean advance together on ocean PEs. |
| **do_concurrent_radiation** | Namelist flag; when `.true.`, radiation is computed concurrently with atmosphere dynamics using a separate OpenMP thread. |
| **FRE** | Flexible Runtime Environment; GFDL's build and workflow system. |
| **coupler_nml** | Fortran namelist used to configure the full coupler (PE counts, thread counts, flags). |
| **diag_table** | FMS diagnostics configuration file; may contain `date_init` to set the model start time. |
| **coupler.res** | Restart file at `INPUT/coupler.res`; stores start date and calendar type for restarts. |
