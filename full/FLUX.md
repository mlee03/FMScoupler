# Flux Exchange — FMSCoupler Component Coupling

*Original authors: Bruce Wyman, V. Balaji, Sergey Malyshev*

## Overview

Six modules couple the atmosphere, ocean, land, and sea-ice components through flux exchange in the FMSCoupler full coupler.

| Module | Role |
|---|---|
| `atm_land_ice_flux_exchange` | Top-level: exchanges fluxes between atmosphere, land, and ice; registers all diagnostic fields. |
| `flux_exchange` | Top-level: initializes all flux-exchange modules; contains subroutines for stock computation. |
| `atmos_ocean_fluxes_calc` | Computes non-deposition gas fluxes between atmosphere and ocean. |
| `atmos_ocean_dep_fluxes_calc` | Computes deposition gas fluxes between atmosphere and ocean. |
| `ice_ocean_flux_exchange` | Exchanges fluxes between ice and ocean. |
| `land_ice_flux_exchange` | Exchanges fluxes between land and ice. |

---

## Design Principles

### Grid Layout

The flux exchange supports **physically independent atmosphere, land, and sea-ice grids**. Ice and ocean must share the same physical grid, though their MPI domain decompositions may differ. Constraints:
- The masked region of the land grid and the ice-ocean grid must tile each other such that every atmosphere grid cell is covered by either land or ice-ocean, but not both.
- The masked regions of the ice and ocean grids must be identical.

The three component grids tile the sphere. `|xxx|` marks a masked (inactive) grid point:

```
ATMOSPHERE  |----|----|----|----|----|----|----|----|
      LAND  |---|---|---|---|xxx|xxx|xxx|xxx|xxx|xxx|
       ICE  |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
     OCEAN  |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
```

### Exchange Methods

| Interface | Exchange Method | Notes |
|---|---|---|
| Atmosphere ↔ land ↔ ice | `xmap_sfc` exchange grid using conservative interpolation (REGRID) | — |
| Land ↔ ice (runoff) | `xmap_runoff` exchange grid | — |
| Ice ↔ ocean | `mpp_redistribute` (REDIST) or direct copy (DIRECT) | No exchange grid needed; grids are physically identical. |
| Atmosphere → ocean | Via ice model: `flux_down_from_atmos` → `flux_ice_to_ocean` | Atmospheric fluxes reach the ocean through the ice model. |

### Time Propagation

Sensible heat flux and surface evaporation can depend **implicitly** on surface temperature. Therefore:
- Land and sea-ice temperature updates must run on the **atmospheric timestep**.
- `update_land_model_fast` and `update_ice_model_fast` must update the surface temperature each atmospheric timestep for the implicit diffusion scheme to be correct.
- Surface fluxes for all other tracers and for momentum are treated as explicit functions of the surface state.
- The module supports simultaneous implicit time integration on both sides of the surface interface.

### Configuring Additional Fields

- Additional tracer and gas-exchange fluxes are configured through `field_table`.
- Additional named boundary-condition fields are configured in the coupler boundary types.
- Any field exchanged between components can be replaced at runtime by a constant or file-based value using FMS `data_override` and `data_table`.
- Each component model must expose a public data type containing the boundary fields needed by the coupler.

---

## Transfer Modes

There are three modes of flux exchange within FMSCoupler, selected by the `xtype` integer field on each boundary type:

| Mode | `xtype` Value | When Used | Mechanism |
|---|---|---|---|
| `REGRID` | 1 | Grids are physically distinct (e.g., atmosphere ↔ land ↔ ice) | Maps data through the exchange grid using conservative interpolation. |
| `REDIST` | 2 | Grids share the same physical grid but have different MPI decompositions (e.g., ice-ocean when `slow_ice_with_ocean=.true.`) | Moves data between PE layouts using `mpp_redistribute`. |
| `DIRECT` | 3 | Physical grid and MPI decomposition are identical (e.g., ice-ocean when `slow_ice_with_ocean` is not used) | Data are copied directly. |

For `REGRID`: data on one component grid is first mapped onto the exchange grid, computations are carried out on the exchange grid, and the result is then mapped to the receiving component grid. Computed fields and fluxes can be overwritten by `data_override`, but the override is applied only if the tracer is specified in the `tracer_table`.

---

## Namelist Parameters (`flux_exchange_nml`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `z_ref_heat` | real | 2.0 m | Reference height for temperature and relative-humidity diagnostics (`t_ref`, `rh_ref`, `del_h`, `del_q`). |
| `z_ref_mom` | real | 10.0 m | Reference height for momentum diagnostics (`u_ref`, `v_ref`, `del_m`). |
| `do_area_weighted_flux` | logical | `.false.` | When `.true.`, fluxes passed to the ocean are multiplied by the ice area fraction before redistribution so the ocean receives the grid-cell-mean flux rather than the per-ice-area flux. |
| `debug_stocks` | logical | `.false.` | Enables additional stock-conservation debug output when `.true.`. |
| `divert_stocks_report` | logical | `.false.` | Redirects stock reporting to `stocks.out` rather than the standard log when `.true.`. |
| `do_runoff` | logical | `.true.` | Enables interpolation of land runoff to the ocean. |
| `do_forecast` | logical | `.false.` | Enables forecast-mode behavior in the flux coupler when `.true.`. |
| `nblocks` | integer | 1 | Number of blocks used to divide `n_xgrid_sfc` for OpenMP parallelism; often set to match `coupler_nml%atmos_nthreads`. If left at 1 when threading is active, the model will emit a warning and reset it automatically. |
| `partition_fprec_from_lprec` | logical | `.false.` | For atmosphere-override experiments where liquid and frozen precipitation are combined: converts liquid precipitation to snow when `t_ref < tfreeze`. |
| `scale_precip_2d` | logical | `.false.` | Rescales `lprec` by a 2-D field read from `data_table`. |

---

## Data Override Capabilities

> **Warning:** The original authors strongly advise against using data override capabilities for non-experts.

Any field in the lists below can be replaced at runtime with a constant or file-based value by adding a matching entry to `data_table`. A data override is applied only when a matching entry exists; otherwise the model-computed value is used unchanged.

### Overridable Fields in `sfc_boundary_layer`

**Atmosphere boundary → exchange grid:**

`t_bot`, `z_bot`, `p_bot`, `u_bot`, `v_bot`, `p_surf`, `slp`, `gust`, and fields in the coupler bc type.

**Ice boundary → exchange grid:**

`t_surf`, `rough_mom`, `rough_heat`, `rough_moist`, `albedo`, `albedo_vis_dir`, `albedo_nir_dir`, `albedo_vis_dif`, `albedo_nir_dif`, `u_surf`, `v_surf`

**Land boundary → exchange grid:**

`t_surf`, `t_ca`, `rough_mom`, `rough_heat`, `albedo`, `tracers`, `albedo_vis_dir`, `albedo_nir_dir`, `albedo_vis_dif`, `albedo_nir_dif`

**Exchange grid → `Land_ice_atmos_boundary`:**

`t`, `albedo`, `albedo_vis_dir`, `albedo_nir_dir`, `albedo_vis_dif`, `albedo_nir_dif`, `land_frac`, `dt_t`, `dt_tr`, `u_flux`, `v_flux`, `dtaudu`, `dtaudv`, `u_star`, `b_star`, `rough_mom`

### Overridable Fields in `flux_down_from_atmos`

**Atmosphere boundary → exchange grid:**

`flux_sw`, `flux_sw_dir`, `flux_sw_dif`, `flux_sw_down_vis_dir`, `flux_sw_down_vis_dif`, `flux_sw_down_total_dir`, `flux_sw_down_total_dif`, `flux_sw_vis`, `flux_sw_vis_dir`, `flux_sw_vis_dif`, `flux_lw`, `lprec`, `frac_precip`, `fprec`, `coszen`, `dtmass`, `delta_t`, `dflux_t`, `delta_tr`, `dflux_tr`

**Exchange grid → land boundary:**

`drag_q`, `lwdn_flux`, `cd_m`, `cd_t`, `bstar`, `ustar`, `wind`, `z_bot`, `t_flux`, `lw_flux`, `sw_flux`, `sw_flux_down_vis_dir`, `sw_flux_down_total_dir`, `sw_flux_down_vis_dif`, `sw_flux_down_total_dif`, `lprec`, `fprec`, `dhdt`, `drdt`, `p_surf`, `tr_flux`, `dfdtr`

**Exchange grid → ice boundary:**

`u_flux`, `v_flux`, `t_flux`, `q_flux`, `lw_flux`, `sw_flux_nir_dir`, `sw_flux_vis_dir`, `sw_flux_nir_dif`, `sw_flux_vis_dif`, `sw_down_vis_dir`, `sw_down_vis_dif`, `sw_down_nir_dir`, `sw_down_nir_dif`, `lprec`, `fprec`, `dhdt`, `dedt`, `drdt`, `coszen`, `p`

### Overridable Fields in `flux_up_to_atmos`

**Ice boundary → atmosphere boundary:** `t_surf`

**Land boundary → atmosphere boundary:** `t_ca`, `t_surf`, `tr`

### Overridable Fields in `flux_land_to_ice`

**Land boundary → ice boundary:** `runoff`, `calving`, `runoff_hflx`, `calving_hflx`

> `do_runoff` namelist flag must be `.true.` for this exchange to occur.

### Overridable Fields in `flux_ice_to_ocean`

**Ice boundary → ocean boundary:**

`u_flux`, `v_flux`, `t_flux`, `q_flux`, `salt_flux`, `lw_flux`, `sw_flux_nir_dir`, `sw_flux_nir_dif`, `sw_flux_vis_dir`, `sw_flux_vis_dif`, `lprec`, `fprec`, `runoff`, `calving`, `runoff_hflx`, `calving_hflx`, `p`, `mi`, `ustar_berg`, `area_berg`, `mass_berg`

### Overridable Fields in `flux_ocean_to_ice`

**Ocean boundary → ice boundary:** `u`, `v`, `t`, `s`, `frazil`, `sea_level`

---

## Diagnostic Fields

All fields below are registered in `atm_land_ice_flux_exchange.F90` inside `diag_field_init`.

**Static fields:** `land_mask`, `height2m`, `height10m`, `sftlf`

**Atmosphere surface fields:**

`ice_mask`, `wind`, `drag_moist`, `drag_heat`, `drag_mom`, `rough_moist`, `rough_heat`, `rough_mom`, `u_star`, `b_star`, `q_star`, `thv_atm`, `thv_surf`, `tau_x`, `tau_y`, `t_ocean`, `t_surf`, `t_ca`, `z_atm`, `p_atm`, `slp`, `gust`, `shflx`, `lwflx`, `t_atm`, `u_atm`, `v_atm`, `t_ref`, `rh_ref`, `rh_ref_cmip`, `u_ref`, `v_ref`, `wind_ref`, `del_h`, `del_m`, `del_q`, `q_ref`, `rough_scale`, `evap`, `co2_bot`

**Atmosphere tracer fields** (per-tracer, name-prefixed):

`*_tot_con_atm`, `*_tot_con_ref`, `*_atm`, `*_surf`, `*_flux`, `*_ref`, `*_mol_flux`, `*_atm_dvmr`, `*_surf_dvmr`, `*_mol_flux_atm0`

**CMIP fields** (registered with `register_cmip_diag_field_2d` or `fms_diag_register_diag_field` with `use_AM3_physics`):

`tas`, `uas`, `vas`, `sfcWind`, `huss`, `hurs`, `rhs`, `ts`, `psl`, `tauu`, `tauv`, `hfss`, `hfls`, `evspsbl`, `tslsi`, `tos`, `sic`

**Global scalar time-series fields** (registered with `register_global_diag_field`, only without `use_AM3_physics`):

`evspsbl`, `ts`, `tas`, `tasl`, `hfss`, `hfls`, `rls`

**Land axes fields** (registered with `register_tiled_diag_field` or `fms_diag_register_diag_field` with `_USE_LEGACY_LAND_`):

`t_ref`, `q_ref`, `rh_ref`, `u_ref`, `v_ref`, `evap`, `shflx`, `tasLut`, `hussLut`, `*_tot_con_atm`, `*_tot_con_ref`, `*_flux`, `*_mol_flux`, `*_ref`

---

## Required Variables in Component Data Types

The following fields must be defined in each component's public data type for the flux exchange to function correctly.

### Atmosphere (`atmos_data_type`)

**Grid fields** (must be grid-box corner coordinates in radians, monotonic order):
`lon_bnd`, `lat_bnd`

**Required bottom-level and surface fields** (primary inputs to `sfc_boundary_layer`):
`t_bot`, `q_bot`, `z_bot`, `p_bot`, `u_bot`, `v_bot`, `p_surf`, `slp`, `gust`

**Radiative and precipitation fields** (passed to land and ice by `flux_down_from_atmos`):
`flux_sw`, `flux_lw`, `lprec`, `fprec`, `coszen`

**Diagnostic axis IDs** (required for FMS diagnostic registration):
`axes`

**Implicit time-stepping fields** (support implicit coupling between atmosphere and surface models):
`dtmass`, `delta_t`, `delta_q`, `dflux_t`, `dflux_q`

### Land (`land_data_type`)

**Grid fields** (grid-box corner coordinates in radians, monotonic):
`lon_bnd`, `lat_bnd`

| Field | Description |
|---|---|
| `mask` | Land-sea mask; `.true.` over land points. |
| `glacier` | Glacier mask; `.true.` over glacier points. |
| `tile_size` | Fractional area of each land tile within the atmospheric grid cell [0–1]. |
| `t_surf` | Surface temperature; used for turbulent flux and radiation calculations. |
| `albedo`, `rough_mom`, `rough_heat` | Surface state fields for turbulent flux and radiation. |
| `t_ca`, `q_ca` | Canopy air temperature and specific humidity; returned to the atmosphere by `flux_up_to_atmos`. |
| `stomatal`, `snow`, `water`, `max_water` | Additional surface properties used in flux parameterisations. |

### Ice (`ice_data_type`)

**Grid fields** (all boundary arrays in radians, monotonic):
`lon_bnd`, `lat_bnd`, `lon_bnd_uv`, `lat_bnd_uv`

**Mask fields:**
`mask` (ocean-land mask for tracer points), `mask_uv` (ocean-land mask for momentum points), `ice_mask` (optional explicit sea-ice mask)

**Coverage fields:**
`part_size` (fractional area of each ice thickness category), `part_size_uv`

**Atmosphere–ice interface fields** (provided to the atmosphere each fast timestep):
`t_surf`, `albedo`, `rough_mom`, `rough_heat`, `rough_moist`, `u_srf`, `v_surf`

**Ice–ocean interface fields** (populated by `flux_down_from_atmos` and `flux_land_to_ice`, then passed to the ocean by `flux_ice_to_ocean`):
`flux_u`, `flux_v`, `flux_t`, `flux_q`, `flux_salt`, `flux_lw`, `flux_sw_vis_dir`, `flux_sw_vis_dif`, `flux_sw_nir_dir`, `flux_sw_nir_dif`, `lprec`, `fprec`, `runoff`, `calving`, `runoff_hflx`, `calving_hflx`, `p_surf`

**Optional iceberg fields** (allocated only when the iceberg module is active):
`ustar_berg`, `area_berg`, `mass_berg`

### Ocean (`ocean_public_type`)

Required fields: `t_surf`, `s_surf`, `u_surf`, `v_surf`, `frazil`, `sea_lev`, `Data%mask`, `Data%mask_uv`, `Ocean%mask`, `Ocean%mask_uv`
