# `land_data_type` — Land Model State (LM4)

## Overview

`land_data_type` carries the state of the land model. 

---

## Surface State Fields

Dimensioned `(grid_index, tile_number)`.

| Field | Units | Description |
|---|---|---|
| `Land%tile_size` | dimensionless, 0–1 | Fractional coverage of the atmospheric grid cell by this tile; used to area-weight tile quantities back onto the atmosphere grid. |
| `Land%t_surf` | K | Ground (radiative) surface temperature; used in longwave radiation and sensible heat flux calculations. |
| `Land%t_ca` | K | Canopy air temperature — near-surface air temperature within the plant canopy layer; differs from `t_surf` over vegetated tiles. |
| `Land%tr` | — | Surface tracer mixing ratios on each tile, including canopy air specific humidity as the first tracer; additional tracers (e.g., CO₂) follow the tracer table order. Dimensioned `(grid_index, tile_number, tracer_index)`. |

---

## Surface Albedo Fields

Dimensioned `(grid_index, tile_number)`.

| Field | Units | Description |
|---|---|---|
| `Land%albedo` | dimensionless | Broadband surface albedo; legacy field — per-band albedos below are preferred. |
| `Land%albedo_vis_dir` | dimensionless | Surface albedo for direct-beam visible radiation (0.2–0.7 µm). |
| `Land%albedo_nir_dir` | dimensionless | Surface albedo for direct-beam near-infrared radiation. |
| `Land%albedo_vis_dif` | dimensionless | Surface albedo for diffuse visible radiation. |
| `Land%albedo_nir_dif` | dimensionless | Surface albedo for diffuse near-infrared radiation. |

---

## Surface Roughness Fields

Dimensioned `(grid_index, tile_number)`.

| Field | Units | Description |
|---|---|---|
| `Land%rough_mom` | m | Surface roughness length for momentum; used in Monin-Obukhov flux calculations. |
| `Land%rough_heat` | m | Surface roughness length for heat and tracers. |
| `Land%rough_scale` | dimensionless | Topographic form-drag scaling factor for momentum; accounts for sub-grid orographic drag. |

---

## Land Discharge and Runoff Fields

Dimensioned `(lon, lat)`.

| Field | Units | Description |
|---|---|---|
| `Land%discharge` | kg/m²/s | Liquid water discharge (river runoff) from land to ocean. |
| `Land%discharge_heat` | W/m² | Sensible heat carried by liquid discharge, using 0 °C as datum. |
| `Land%discharge_snow` | kg/m²/s | Solid water (snow/ice) discharge from land to ocean. |
| `Land%discharge_snow_heat` | W/m² | Sensible heat carried by solid discharge, using 0 °C as datum. |

---

## Domain and PE Metadata

| Field | Type | Description |
|---|---|---|
| `Land%mask` | `logical 2D` | `.true.` where the grid cell contains land; used to gate land-only computations. |
| `Land%axes(1)` | integer | Diag-manager axis ID for the unstructured land grid; used when registering tiled land diagnostics. |
| `Land%domain` | `type(domain2D)` | FMS structured-grid domain for the land model; used for halo exchanges and exchange-grid setup. |
| `Land%ug_domain` | `type(domainUG)` | FMS unstructured-grid domain for the land model; carries the tile-based decomposition used by LM4. |
| `Land%pelist` | `integer 1D` | List of MPI PE numbers on which the land model is running. |
| `Land%pe` | logical | `.true.` on PEs that are part of the land pelist; used to gate land stock calculations. |
