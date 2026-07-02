# `atmos_ice_boundary_type` — Atmosphere-to-Ice Boundary Fields

## Overview

`atmos_ice_boundary_type` holds all data passed from the coupler to the sea-ice model (SIS2) at each atmospheric timestep. An instance named `Atmos_ice_boundary` is declared in `coupler_main.F90`. Fields are 3D arrays dimensioned `(:, :, n_categories)` where the third dimension indexes the sea-ice thickness categories; category 1 corresponds to open ocean.

**Populated by:** `flux_down_from_atmos`, `sfc_boundary_layer`  
**Consumed by:** `update_ice_model_fast`, `flux_ice_to_ocean`  
**Related types:** `atmos_land_boundary_type`, `land_ice_atmos_boundary_type`, `ice_ocean_boundary_type`

---

## Wind Stress Fields

All fields are real 3D arrays on an A-grid (not rotated to the model grid).

| Field | Units | Description |
|---|---|---|
| `Atmos_ice_boundary%u_flux` | Pa | True-eastward wind stress from the atmosphere to the ocean or ice in each thickness category. |
| `Atmos_ice_boundary%v_flux` | Pa | True-northward wind stress from the atmosphere to the ocean or ice in each thickness category. |
| `Atmos_ice_boundary%u_star` | Pa | Atmospheric friction velocity on an A-grid. |

---

## Sensible Heat and Moisture Flux Fields

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Atmos_ice_boundary%t_flux` | real 3D | W/m² | Net sensible heat flux from the ocean or ice surface into the atmosphere. |
| `Atmos_ice_boundary%q_flux` | real 3D | kg/m²/s | Moisture flux from the ice or ocean to the atmosphere due to evaporation or sublimation. |

---

## Implicit Coupling Derivative Fields

These linearisation (derivative) terms are required to close the **implicit tridiagonal surface diffusion scheme** between the atmosphere and ice. They are populated during `update_atmos_model_down` (forward sweep) and used by `update_ice_model_fast` to update the ice surface temperature.

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Atmos_ice_boundary%dhdt` | real 3D | W/m²/°C | `d(upward sensible heat flux)/d(T_surf)` — derivative of sensible heat flux with respect to surface temperature. |
| `Atmos_ice_boundary%dedt` | real 3D | kg/m²/s/°C | `d(sublimation+evaporation rate)/d(T_surf)` — derivative of the moisture flux with respect to surface temperature. |
| `Atmos_ice_boundary%drdt` | real 3D | W/m²/°C | `d(net upward longwave flux)/d(T_surf)` — derivative of the net upward longwave flux with respect to surface temperature. |

---

## Radiation Flux Fields

All fields are real 3D arrays in W/m².

| Field | Units | Description |
|---|---|---|
| `Atmos_ice_boundary%lw_flux` | W/m² | Net downward longwave radiation flux from the atmosphere into the ice or ocean. |
| `Atmos_ice_boundary%sw_flux_vis_dir` | W/m² | Net direct visible shortwave radiation flux into the ice or ocean. |
| `Atmos_ice_boundary%sw_flux_vis_dif` | W/m² | Net diffuse visible shortwave radiation flux into the ice or ocean. |
| `Atmos_ice_boundary%sw_flux_nir_dir` | W/m² | Net direct near-infrared shortwave radiation flux into the ice or ocean. |
| `Atmos_ice_boundary%sw_flux_nir_dif` | W/m² | Net diffuse near-infrared shortwave radiation flux into the ice or ocean. |
| `Atmos_ice_boundary%sw_down_vis_dir` | W/m² | Downward direct visible shortwave radiation flux from the atmosphere. |
| `Atmos_ice_boundary%sw_down_vis_dif` | W/m² | Downward diffuse visible shortwave radiation flux from the atmosphere. |
| `Atmos_ice_boundary%sw_down_nir_dir` | W/m² | Downward direct near-infrared shortwave radiation flux from the atmosphere. |
| `Atmos_ice_boundary%sw_down_nir_dif` | W/m² | Downward diffuse near-infrared shortwave radiation flux from the atmosphere. |
| `Atmos_ice_boundary%coszen` | dimensionless (≤ 1) | Cosine of the solar zenith angle averaged over the next radiation timestep (not the timestep used to compute the `sw_flux` fields). |

---

## Precipitation Fields

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Atmos_ice_boundary%lprec` | real 3D | kg/m²/s | Liquid precipitation (rain) from the atmosphere onto the ice or ocean in each thickness category. Rain falling on snow is currently assumed to drain directly through the ice into the ocean. |
| `Atmos_ice_boundary%fprec` | real 3D | kg/m²/s | Frozen precipitation (snowfall, sleet, hail, graupel) from the atmosphere to the ice or ocean. All forms of frozen precipitation are treated as snow in SIS2. |

---

## Metadata and Transfer Fields

| Field | Type | Description |
|---|---|---|
| `Atmos_ice_boundary%p` | real 3D | Atmospheric surface pressure [Pa]; typically ~1×10⁵ Pa. |
| `Atmos_ice_boundary%xtype` | integer | Transfer mode for the atmosphere-to-ice exchange: `REGRID` (1), `REDIST` (2), or `DIRECT` (3). |
| `Atmos_ice_boundary%fluxes` | `type(coupler_3d_bc_type)` | Array of additional per-tracer gas and deposition fluxes from the atmosphere to the ice. |
