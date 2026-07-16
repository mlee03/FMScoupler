# `atmos_land_boundary_type` — Atmosphere-to-Land Boundary Fields

## Overview

`atmos_land_boundary_type` carries all data passed from atmos to the land model (LM4).  All fields are pointers dimensioned `(grid_index, tile_number)` unless noted otherwise; the tile dimension supports LM4's unstructured multi-tile land representation.

--- 

## Radiation Flux Fields

These fields are real 2D arrays in W/m², dimensioned `(grid_index, tile_number)`.

| Field | Units | Description |
|---|---|---|
| `Atmos_land_boundary%t_flux` | W/m² | Sensible heat flux into the land surface. |
| `Atmos_land_boundary%lw_flux` | W/m² | Net longwave radiation flux at the land surface. |
| `Atmos_land_boundary%lwdn_flux` | W/m² | Downward longwave radiation flux at the land surface. |
| `Atmos_land_boundary%sw_flux` | W/m² | Net shortwave radiation flux at the land surface. |
| `Atmos_land_boundary%swdn_flux` | W/m² | Downward shortwave radiation flux at the land surface. |
| `Atmos_land_boundary%sw_flux_down_vis_dir` | W/m² | Downward direct-beam visible shortwave flux. |
| `Atmos_land_boundary%sw_flux_down_total_dir` | W/m² | Downward direct-beam total (broadband) shortwave flux. |
| `Atmos_land_boundary%sw_flux_down_vis_dif` | W/m² | Downward diffuse visible shortwave flux. |
| `Atmos_land_boundary%sw_flux_down_total_dif` | W/m² | Downward diffuse total (broadband) shortwave flux. |

---

## Precipitation Fields

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Atmos_land_boundary%lprec` | real 2D | kg/m²/s | Liquid precipitation rate. |
| `Atmos_land_boundary%fprec` | real 2D | kg/m²/s | Frozen precipitation rate. |
| `Atmos_land_boundary%tprec` | real 2D | K | Temperature of precipitation. |

---

## Derive terms

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Atmos_land_boundary%dhdt` | real 2D | W/m²/K | `d(sensible heat flux)/d(T_surf)` — derivative of sensible heat flux with respect to surface temperature. |
| `Atmos_land_boundary%dhdq` | real 2D | W/m²/(kg/kg) | `d(sensible heat flux)/d(q_surf)` — derivative of sensible heat flux with respect to surface specific humidity. |
| `Atmos_land_boundary%drdt` | real 2D | W/m²/K | `d(longwave flux)/d(T_surf)` — derivative of longwave flux with respect to surface radiative temperature. |

---

## Turbulence and Surface Layer Fields

All fields are real 2D arrays dimensioned `(grid_index, tile_number)`.

| Field | Units | Description |
|---|---|---|
| `Atmos_land_boundary%cd_m` | dimensionless | Drag coefficient for momentum. |
| `Atmos_land_boundary%cd_t` | dimensionless | Drag coefficient for tracers (heat and moisture). |
| `Atmos_land_boundary%ustar` | m/s | Turbulent wind scale (friction velocity). |
| `Atmos_land_boundary%bstar` | m/s | Turbulent buoyancy scale. |
| `Atmos_land_boundary%wind` | m/s | Absolute wind speed at the bottom of the atmospheric layer. |
| `Atmos_land_boundary%z_bot` | m | Height of the bottom atmospheric layer above the land surface. |
| `Atmos_land_boundary%drag_q` | m/s | Product of the moisture drag coefficient and wind speed (`cd_q × wind`); used in land surface moisture flux calculations. |
| `Atmos_land_boundary%p_surf` | Pa | Surface pressure. |

---

## Tracer Flux Fields

Dimensioned `(grid_index, tile_number, tracer_index)`.

| Field | Units | Description |
|---|---|---|
| `Atmos_land_boundary%tr_flux` | tracer units · kg air / (m²·s) | Flux of each tracer into the land surface, including water vapor flux. |
| `Atmos_land_boundary%dfdtr` | varies | `d(tracer flux)/d(tracer_surf)` — derivative of the tracer flux with respect to the surface tracer value, including evaporation over surface specific humidity. |

---

## Metadata Field

| Field | Type | Description |
|---|---|---|
| `Atmos_land_boundary%xtype` | integer | Transfer mode for the atmosphere-to-land exchange: `REGRID` (1), `REDIST` (2), or `DIRECT` (3). |
