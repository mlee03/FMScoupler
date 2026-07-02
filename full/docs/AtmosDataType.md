# `atmos_data_type` — Atmosphere Model State

## Overview

`atmos_data_type` is the main Fortran derived type that holds all fields and states of the atmosphere model. An instance named `Atm` is declared in `coupler_main.F90` and passed to flux-exchange subroutines throughout the full coupler. The type carries the lowest-level atmospheric state, radiative and precipitation fluxes, grid/domain metadata, implicit-coupling coefficients (`Surf_diff`), and grid geometry (`grid`).

**Related types:** `land_data_type`, `ice_data_type`, `atmos_land_boundary_type`, `atmos_ice_boundary_type`, `land_ice_atmos_boundary_type`

**Key subroutines that read or write `Atm` fields:** `sfc_boundary_layer`, `flux_down_from_atmos`, `flux_up_to_atmos`, `update_atmos_model_down`, `update_atmos_model_up`

---

## Grid and Domain Fields

| Field | Type / Dimensions | Description |
|---|---|---|
| `Atm%domain` | `type(domain2d)` | FMS domain decomposition for the atmosphere; defines the MPI tile layout and halo widths. |
| `Atm%axes` | `integer(4)` | Diag-manager axis indices for x, y, pfull, and phalf; used when registering and sending diagnostic fields. |
| `Atm%lon_bnd` | `real 2D` | Longitude of grid-box corners on the local compute domain [radians]. |
| `Atm%lat_bnd` | `real 2D` | Latitude of grid-box corners on the local compute domain [radians]. |
| `Atm%lon` | `real 2D` | Longitude of grid-box centres on the local compute domain [radians]. |
| `Atm%lat` | `real 2D` | Latitude of grid-box centres on the local compute domain [radians]. |
| `Atm%grid` | `type(grid_box_type)` | Grid geometry needed for second-order conservative remapping on the cubic-sphere exchange grid. See [Grid Geometry Fields](#grid-geometry-fields-atmgrid) below. |
| `Atm%maskmap` | `logical(:,:)` (pointer) | Mask indicating which logical processors are active for ocean code; processors covering all-land points may not be assigned to physical PEs. Dummy field — must be present for compilation but need not be set. |

---

## Lowest Atmospheric Level State Fields

These fields carry the atmospheric state at the bottom model level and are the primary inputs to `sfc_boundary_layer`. All are 2D arrays on the local compute domain.

| Field | Units | Description |
|---|---|---|
| `Atm%t_bot` | K | Temperature at the lowest model level. |
| `Atm%tr_bot` | — | Tracer mixing ratios at the lowest model level; 3D array, third dimension indexes the tracer table. Specific humidity (`sphum`) is always present. |
| `Atm%z_bot` | m | Height of the lowest model level above the surface. |
| `Atm%p_bot` | Pa | Pressure at the lowest model level. |
| `Atm%u_bot` | m/s | Zonal wind component at the lowest model level. |
| `Atm%v_bot` | m/s | Meridional wind component at the lowest model level. |
| `Atm%p_surf` | Pa | Surface pressure. |
| `Atm%slp` | Pa | Sea-level pressure. |
| `Atm%gust` | m/s | Gustiness factor — a minimum wind speed added in quadrature to the resolved wind to account for sub-grid convective gusts in surface flux calculations. |
| `Atm%coszen` | dimensionless | Cosine of the solar zenith angle; used to weight shortwave fluxes and partition direct vs. diffuse radiation. |

---

## Radiative Flux Fields

All fields are real 2D arrays. These fluxes are computed by the atmosphere and passed to land and ice by `flux_down_from_atmos`.

| Field | Units | Description |
|---|---|---|
| `Atm%flux_sw` | W/m² | Total net shortwave flux at the surface (absorbed by the surface). |
| `Atm%flux_sw_dir` | W/m² | Direct-beam component of the net shortwave flux. |
| `Atm%flux_sw_dif` | W/m² | Diffuse component of the net shortwave flux. |
| `Atm%flux_sw_down_vis_dir` | W/m² | Downward direct-beam flux in the visible band (0.2–0.7 µm). |
| `Atm%flux_sw_down_vis_dif` | W/m² | Downward diffuse flux in the visible band. |
| `Atm%flux_sw_down_total_dir` | W/m² | Downward direct-beam broadband shortwave flux. |
| `Atm%flux_sw_down_total_dif` | W/m² | Downward diffuse broadband shortwave flux. |
| `Atm%flux_sw_vis` | W/m² | Net (downward minus reflected) visible-band shortwave flux at the surface. |
| `Atm%flux_sw_vis_dir` | W/m² | Direct-beam component of the net visible shortwave flux. |
| `Atm%flux_sw_vis_dif` | W/m² | Diffuse component of the net visible shortwave flux. |
| `Atm%flux_lw` | W/m² | Net downward longwave flux at the surface. |

---

## Precipitation Fields

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Atm%lprec` | real 2D | kg/m² | Mass of liquid precipitation accumulated since the last time step; equivalent to a rate in kg/m²/s when divided by `dt_atmos`. |
| `Atm%fprec` | real 2D | kg/m² | Mass of frozen (solid) precipitation accumulated since the last time step. |

---

## Generic Exchange and Tracer Boundary Condition Fields

| Field | Type / Dimensions | Description |
|---|---|---|
| `Atm%gex_atm2lnd` | real 3D | Generic exchange fields sent from the atmosphere to the land model (e.g., CO₂, aerosol deposition); third dimension indexes the exchange field list as defined in `field_table`. |
| `Atm%gex_lnd2atm` | real 3D | Generic exchange fields returned from the land model to the atmosphere (e.g., surface emission fluxes); third dimension indexes the exchange field list. |
| `Atm%fields` | `type(coupler_2d_bc_type)` | Array of additional tracer boundary-condition fields for atmosphere–ocean gas exchange (CO₂, O₂, CFCs, etc.); registered and populated by `atmos_tracer_flux_init`. |

---

## Time and PE Metadata Fields

| Field | Type | Description |
|---|---|---|
| `Atm%Time` | `type(time_type)` | Current model time; passed to `diag_manager` `send_data` calls and to `fms_data_override`. |
| `Atm%Time_step` | `type(time_type)` | Atmospheric model timestep duration (`dt_atmos`). |
| `Atm%Time_init` | `type(time_type)` | Reference (initial) time for the model run. |
| `Atm%pelist` | `integer 1D` | List of MPI PE numbers on which the atmosphere is running. |
| `Atm%pe` | `logical` | `.true.` on PEs that are part of the atmosphere pelist; used to gate atmosphere-only code blocks. |

---

## Implicit Vertical Diffusion Coefficients (`Atm%Surf_diff`)

`Atm%Surf_diff` is of type `surf_diff_type`, defined in `atmos_phys/atmos_param/vert_diff/vert_diff.F90`. It carries the forward-elimination coefficients from the implicit vertical diffusion scheme that couples the atmosphere to the surface models. Fields are accessed as, e.g., `Atm%Surf_diff%dtmass`.

These fields support the **tridiagonal implicit surface coupling** between the atmosphere and land/ice. They are populated during `update_atmos_model_down` (forward sweep) and consumed during `update_atmos_model_up` (back-substitution) and `flux_down_from_atmos` / `flux_up_to_atmos`.

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Atm%Surf_diff%dtmass` | real 2D | s·m²/kg | `dt/mass` — ratio of the atmospheric timestep to the surface-layer air mass; scales flux tendencies to temperature/tracer tendencies. |
| `Atm%Surf_diff%dflux_t` | real 2D | W/m²/K | `d(sensible heat flux)/d(T_surf)` — linearisation of the surface heat flux with respect to surface temperature; used to form the implicit coupling term. |
| `Atm%Surf_diff%delta_t` | real 2D | K | Forward-elimination coefficient for temperature from the implicit tridiagonal scheme; represents the accumulated atmospheric temperature forcing at the bottom level waiting for the surface response. |
| `Atm%Surf_diff%delta_u` | real 2D | m/s | Forward-elimination coefficient for zonal wind from the implicit scheme. |
| `Atm%Surf_diff%delta_v` | real 2D | m/s | Forward-elimination coefficient for meridional wind from the implicit scheme. |
| `Atm%Surf_diff%dflux_tr` | real 3D | varies | `d(tracer flux)/d(tracer_surf)` — linearisation of tracer surface fluxes with respect to surface tracer concentration; third dimension indexes tracers. |
| `Atm%Surf_diff%delta_tr` | real 3D | varies | Forward-elimination coefficient for each tracer from the implicit scheme; third dimension indexes tracers. |
| `Atm%Surf_diff%tdt_dyn` | real 3D | K/s | Temperature tendency from dynamics (advection, etc.) passed through the diffusion scheme. |
| `Atm%Surf_diff%qdt_dyn` | real 3D | kg/kg/s | Moisture tendency from dynamics. |
| `Atm%Surf_diff%dgz_dyn` | real 3D | m²/s³ | Geopotential height tendency from dynamics. |
| `Atm%Surf_diff%ddp_dyn` | real 3D | Pa/s | Pressure-thickness tendency from dynamics. |
| `Atm%Surf_diff%tdt_rad` | real 3D | K/s | Temperature tendency from radiation; used in MIZ (marginal ice zone) forecast mode. |

---

## Grid Geometry Fields (`Atm%grid`)

`Atm%grid` is of type `grid_box_type`, defined in `FMS/exchange/xgrid`. It holds the geometric quantities needed for **second-order conservative flux remapping** between the atmosphere and surface component grids on a cubic-sphere mesh. Fields are accessed as, e.g., `Atm%grid%dx`.

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Atm%grid%dx` | real 2D | m | Grid-box width in the x-direction. |
| `Atm%grid%dy` | real 2D | m | Grid-box width in the y-direction. |
| `Atm%grid%area` | real 2D | m² | Grid-box area. |
| `Atm%grid%edge_w` | real 1D | m | Western edge lengths of grid boxes along the boundary. |
| `Atm%grid%edge_e` | real 1D | m | Eastern edge lengths. |
| `Atm%grid%edge_s` | real 1D | m | Southern edge lengths. |
| `Atm%grid%edge_n` | real 1D | m | Northern edge lengths. |
| `Atm%grid%en1` | real 3D | — | First unit normal vector at grid-box edges; used to project vector fields (winds, stresses) during remapping. |
| `Atm%grid%en2` | real 3D | — | Second unit normal vector at grid-box edges. |
| `Atm%grid%vlon` | real 3D | — | Unit vector in the local longitude direction at each grid point; used to rotate between geographic and local coordinate frames during exchange. |
| `Atm%grid%vlat` | real 3D | — | Unit vector in the local latitude direction at each grid point. |
