# `land_ice_atmos_boundary_type` — Aggregated Surface State Returned to the Atmosphere

## Overview

`land_ice_atmos_boundary_type` contains surface quantities passed from the land and ice models to the atmosphere. All quantities are on the exchange grid.

---

## Surface Temperature and Albedo Fields

All fields are real 2D arrays on the exchange grid.

| Field | Units | Description |
|---|---|---|
| `Land_ice_atmos_boundary%t` | K | Area-weighted surface temperature seen by the atmosphere for radiation calculations; weighted over land and ice fractions. |
| `Land_ice_atmos_boundary%t_ocean` | K | Ocean surface temperature for radiation calculations; sourced from `Ice%t_surf` through the exchange grid. |
| `Land_ice_atmos_boundary%albedo` | dimensionless | Broadband surface albedo. |
| `Land_ice_atmos_boundary%albedo_vis_dir` | dimensionless | Direct-beam visible-band surface albedo. |
| `Land_ice_atmos_boundary%albedo_nir_dir` | dimensionless | Direct-beam near-infrared surface albedo. |
| `Land_ice_atmos_boundary%albedo_vis_dif` | dimensionless | Diffuse visible-band surface albedo. |
| `Land_ice_atmos_boundary%albedo_nir_dif` | dimensionless | Diffuse near-infrared surface albedo. |
| `Land_ice_atmos_boundary%land_frac` | dimensionless | Fraction of the atmospheric grid cell covered by land. |
| `Land_ice_atmos_boundary%frac_open_sea` | dimensionless | Non-sea-ice fraction of the grid cell; complement of the sea-ice concentration. |
| `Land_ice_atmos_boundary%rough_mom` | m | Area-weighted surface roughness length for momentum. |
| `Land_ice_atmos_boundary%rough_heat` | m | Area-weighted surface roughness length for heat. |

---

## Reference-Height Diagnostic Fields

All are real 2D arrays.

| Field | Units | Description |
|---|---|---|
| `Land_ice_atmos_boundary%u_ref` | m/s | Zonal wind at the momentum reference height (`z_ref_mom`). |
| `Land_ice_atmos_boundary%v_ref` | m/s | Meridional wind at the momentum reference height (`z_ref_mom`). |
| `Land_ice_atmos_boundary%t_ref` | K | Air temperature at the heat reference height (`z_ref_heat`). |
| `Land_ice_atmos_boundary%q_ref` | kg/kg | Specific humidity at the heat reference height (`z_ref_heat`). |
| `Land_ice_atmos_boundary%wind` | m/s | Absolute wind speed at the lowest atmospheric model level including gust corrections. |
| `Land_ice_atmos_boundary%thv_atm` | K | Virtual potential temperature at the lowest atmospheric model level. |
| `Land_ice_atmos_boundary%thv_surf` | K | Virtual potential temperature at the surface. |

---

## Implicit Coupling Output Fields

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Land_ice_atmos_boundary%dt_t` | real 2D | K/s | Temperature tendency correction at the lowest atmospheric level from the implicit surface flux scheme. |
| `Land_ice_atmos_boundary%dt_tr` | real 3D | tracer units/s | Tracer mixing-ratio tendency correction at the lowest level; third dimension indexes tracers. |

---

## Wind Stress and Turbulence Fields

All fields are real 2D arrays.

| Field | Units | Description |
|---|---|---|
| `Land_ice_atmos_boundary%u_flux` | Pa | Zonal wind stress on the atmosphere. |
| `Land_ice_atmos_boundary%v_flux` | Pa | Meridional wind stress on the atmosphere. |
| `Land_ice_atmos_boundary%dtaudu` | Pa·s/m | `d(zonal wind stress)/d(u)` — implicit coupling coefficient for zonal momentum. |
| `Land_ice_atmos_boundary%dtaudv` | Pa·s/m | `d(meridional wind stress)/d(v)` — implicit coupling coefficient for meridional momentum. |
| `Land_ice_atmos_boundary%u_star` | m/s | Friction velocity (surface turbulent velocity scale). |
| `Land_ice_atmos_boundary%b_star` | m/s² | Buoyancy scale used in Monin-Obukhov similarity theory. |
| `Land_ice_atmos_boundary%q_star` | kg/kg | Moisture scale used in Monin-Obukhov similarity theory. |

---

## Surface Heat Flux Fields

> **Note:** `shflx` and `lhflx` are not compiled when `use_AM3_physics` is defined.

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Land_ice_atmos_boundary%shflx` | real 2D | W/m² | Sensible heat flux at the surface. |
| `Land_ice_atmos_boundary%lhflx` | real 2D | W/m² | Latent heat flux at the surface. |

---

## Internal and Generic Exchange Fields

| Field | Type | Description |
|---|---|---|
| `Land_ice_atmos_boundary%data` | real 3D | Collective array providing named access to the scalar fields above; used internally for data-override and exchange-grid operations. |
| `Land_ice_atmos_boundary%gex_lnd2atm` | real 3D | Generic exchange fields returned from the land model to the atmosphere (e.g., surface emission fluxes); third dimension indexes the exchange field list. |
| `Land_ice_atmos_boundary%xtype` | integer | Transfer mode for the exchange-grid-to-atmosphere remap: `REGRID` (1), `REDIST` (2), or `DIRECT` (3). |
