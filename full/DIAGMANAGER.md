# Diagnostic Manager Variables in FMSCoupler/full

This document lists every diagnostic field registered via `fms_diag_register_diag_field`,
`fms_diag_register_static_field`, `register_tiled_diag_field`, `register_cmip_diag_field_2d`,
or `register_global_diag_field` in the top-level `.F90` files of `FMSCoupler/full`, organised
by source file and enclosing subroutine.

"send_data?" indicates whether the returned ID variable is subsequently passed to
`fms_diag_send_data`, `send_tile_data`, or `send_global_diag` somewhere in the same file.

---

## `atm_land_ice_flux_exchange.F90`

All register calls in this file appear inside `diag_field_init` (lines 4432–5149).
`send_data` calls are scattered across `sfc_boundary_layer`, `flux_down_from_atmos`,
`flux_up_to_atmos`, `update_ice_model_slow_finish`, `flux_atmos_to_ocean`, and the
helper `diag_sic` — the relevant subroutine is noted in parentheses in the
"send_data?" column where it differs from an obvious single location.

### `diag_field_init` — static fields

| Field name | ID variable | Description | Units | send_data? |
|---|---|---|---|---|
| `'land_mask'` | `id_land_mask` | fractional amount of land | `none` | Yes (`sfc_boundary_layer`) |
| `'height2m'` | `id_height2m` | Height (scalar axis, 2 m reference level) | `m` | Yes (`sfc_boundary_layer`) |
| `'height10m'` | `id_height10m` | Height (scalar axis, 10 m reference level) | `m` | Yes (`sfc_boundary_layer`) |
| `'sftlf'` | `id_sftlf` | Fraction of the Grid Cell Occupied by Land | `1.0` | Yes (`sfc_boundary_layer`) |

### `diag_field_init` — time-varying fields on the atmosphere grid (`mod_name` / `atmos_axes`)

| Field name | ID variable | Description | Units | send_data? |
|---|---|---|---|---|
| `'ice_mask'` | `id_ice_mask` | fractional amount of sea ice | `none` | Yes (`diag_sic`) |
| `'wind'` | `id_wind` | wind speed for flux calculations | `m/s` | Yes (`sfc_boundary_layer`) |
| `'drag_moist'` | `id_drag_moist` | drag coeff for moisture | `none` | Yes (`sfc_boundary_layer`) |
| `'drag_heat'` | `id_drag_heat` | drag coeff for heat | `none` | Yes (`sfc_boundary_layer`) |
| `'drag_mom'` | `id_drag_mom` | drag coeff for momentum | `none` | Yes (`sfc_boundary_layer`) |
| `'rough_moist'` | `id_rough_moist` | surface roughness for moisture | `m` | Yes (`sfc_boundary_layer`) |
| `'rough_heat'` | `id_rough_heat` | surface roughness for heat | `m` | Yes (`sfc_boundary_layer`) |
| `'rough_mom'` | `id_rough_mom` | surface roughness for momentum | `m` | Yes (`sfc_boundary_layer`) |
| `'u_star'` | `id_u_star` | friction velocity | `m/s` | Yes (`sfc_boundary_layer`) |
| `'b_star'` | `id_b_star` | buoyancy scale | `m/s2` | Yes (`sfc_boundary_layer`) |
| `'q_star'` | `id_q_star` | moisture scale | `kg water/kg air` | Yes (`sfc_boundary_layer`) |
| `'thv_atm'` | `id_thv_atm` | surface air virtual potential temperature | `K` | Yes (`sfc_boundary_layer`) |
| `'thv_surf'` | `id_thv_surf` | surface virtual potential temperature | `K` | Yes (`sfc_boundary_layer`) |
| `'tau_x'` | `id_u_flux` | zonal wind stress | `pa` | Yes (`flux_down_from_atmos`) |
| `'tau_y'` | `id_v_flux` | meridional wind stress | `pa` | Yes (`flux_down_from_atmos`) |
| `'t_ocean'` | `id_t_ocean` | surface temperature from ocean output | `deg_k` | Yes (`flux_up_to_atmos`) |
| `'t_surf'` | `id_t_surf` | surface temperature | `deg_k` | Yes (`flux_up_to_atmos`) |
| `'t_ca'` | `id_t_ca` | canopy air temperature | `deg_k` | Yes (`flux_up_to_atmos`) |
| `'z_atm'` | `id_z_atm` | height of btm level | `m` | Yes (`sfc_boundary_layer`) |
| `'p_atm'` | `id_p_atm` | pressure at btm level | `pa` | Yes (`sfc_boundary_layer`) |
| `'slp'` | `id_slp` | sea level pressure | `pa` | Yes (`sfc_boundary_layer`) |
| `'gust'` | `id_gust` | gust scale | `m/s` | Yes (`sfc_boundary_layer`) |
| `'shflx'` | `id_t_flux` | sensible heat flux | `w/m2` | Yes (`flux_up_to_atmos`) |
| `'lwflx'` | `id_r_flux` | net (down-up) longwave flux | `w/m2` | Yes (`flux_up_to_atmos`) |
| `'t_atm'` | `id_t_atm` | temperature at btm level | `deg_k` | Yes (`sfc_boundary_layer`) |
| `'u_atm'` | `id_u_atm` | u wind component at btm level | `m/s` | Yes (`sfc_boundary_layer`) |
| `'v_atm'` | `id_v_atm` | v wind component at btm level | `m/s` | Yes (`sfc_boundary_layer`) |
| `'t_ref'` | `id_t_ref` | temperature at reference height (label_zh) | `deg_k` | Yes (`sfc_boundary_layer`) |
| `'rh_ref'` | `id_rh_ref` | relative humidity at reference height (label_zh) | `percent` | Yes (`sfc_boundary_layer`) |
| `'rh_ref_cmip'` | `id_rh_ref_cmip` | relative humidity at reference height (label_zh) | `percent` | Yes (`sfc_boundary_layer`) |
| `'u_ref'` | `id_u_ref` | zonal wind component at reference height (label_zm) | `m/s` | Yes (`sfc_boundary_layer`) |
| `'v_ref'` | `id_v_ref` | meridional wind component at reference height (label_zm) | `m/s` | Yes (`sfc_boundary_layer`) |
| `'wind_ref'` | `id_wind_ref` | absolute value of wind at reference height (label_zm) | `m/s` | Yes (`sfc_boundary_layer`) |
| `'del_h'` | `id_del_h` | ref height interp factor for heat | `none` | Yes (`sfc_boundary_layer`) |
| `'del_m'` | `id_del_m` | ref height interp factor for momentum | `none` | Yes (`sfc_boundary_layer`) |
| `'del_q'` | `id_del_q` | ref height interp factor for moisture | `none` | Yes (`sfc_boundary_layer`) |
| `'q_ref'` | `id_q_ref` | specific humidity at reference height (label_zh) | `kg/kg` | Yes (`sfc_boundary_layer`) |
| `'rough_scale'` | `id_rough_scale` | topographic scaling factor for momentum drag | `1` | Yes (`sfc_boundary_layer`) |
| `'evap'` | `id_q_flux` | evaporation rate | `kg/m2/s` | Yes (`flux_up_to_atmos`) |
| `'co2_bot'` | `id_co2_bot` | co2_bot from data_override | `ppmv` | Yes (`sfc_boundary_layer`) |

#### Per-tracer fields on the atmosphere grid (loop over `n_exch_tr`; `name` / `longname` / `units` from tracer table)

| Field name | ID variable | Description | Units | send_data? |
|---|---|---|---|---|
| `trim(name)//'_tot_con_atm'` | `id_tr_con_atm(tr)` | vd of `longname` | `m/s` | Yes (`flux_up_to_atmos`) |
| `trim(name)//'_tot_con_ref'` | `id_tr_con_ref(tr)` | vd of `longname` at reference height | `m/s` | Yes (`flux_up_to_atmos`) |
| `trim(name)//'_atm'` | `id_tr_atm(tr)` | `longname` at btm level | `units` | Yes (`sfc_boundary_layer`) |
| `trim(name)//'_surf'` | `id_tr_surf(tr)` | `longname` at the surface | `units` | Yes (`flux_up_to_atmos`) |
| `trim(name)//'_flux'` | `id_tr_flux(tr)` | flux of `longname` | `units kg air/(m2 s)` | Yes (`flux_up_to_atmos`) |
| `trim(name)//'_ref'` | `id_tr_ref(tr)` | `longname` at reference height *(skipped for sphum)* | `units` | Yes (`sfc_boundary_layer`) |
| `trim(name)//'_mol_flux'` | `id_tr_mol_flux(tr)` | flux of `longname` | `mol CO2/(m2 s)` (CO2) or `mol/(m2 s)` (other) | Yes (`flux_up_to_atmos`) |
| `trim(name)//'_atm_dvmr'` | `id_co2_atm_dvmr` | `longname` at btm level *(CO2 only)* | `mol CO2 /mol air` | Yes (`sfc_boundary_layer`) |
| `trim(name)//'_surf_dvmr'` | `id_co2_surf_dvmr` | `longname` at the surface *(CO2 only)* | `mol CO2 /mol air` | Yes (`flux_up_to_atmos`) |
| `trim(name)//'_mol_flux_atm0'` | `id_tr_mol_flux0(tr)` | gross flux of `longname` | `mol/(m2 s)` | Yes (`sfc_boundary_layer`) |

### `diag_field_init` — CMIP / use_AM3_physics fields on the atmosphere grid

These fields are registered under two mutually exclusive code paths (the `use_AM3_physics`
preprocessor macro and the default `register_cmip_diag_field_2d` path). In both cases the
field names, descriptions, and units are identical; only the registration function differs.

| Field name | ID variable | Description | Units | Registration function | send_data? |
|---|---|---|---|---|---|
| `'tas'` | `id_tas` | Near-Surface Air Temperature | `K` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`sfc_boundary_layer`) |
| `'uas'` | `id_uas` | Eastward Near-Surface Wind | `m s-1` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`sfc_boundary_layer`) |
| `'vas'` | `id_vas` | Northward Near-Surface Wind | `m s-1` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`sfc_boundary_layer`) |
| `'sfcWind'` | `id_sfcWind` | Near-Surface Wind Speed | `m s-1` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`sfc_boundary_layer`) |
| `'huss'` | `id_huss` | Near-Surface Specific Humidity | `1.0` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`sfc_boundary_layer`) |
| `'hurs'` | `id_hurs` | Near-Surface Relative Humidity | `%` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`sfc_boundary_layer`) |
| `'rhs'` | `id_rhs` | Near-Surface Relative Humidity | `%` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`sfc_boundary_layer`) |
| `'ts'` | `id_ts` | Surface Temperature | `K` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`flux_up_to_atmos`) |
| `'psl'` | `id_psl` | Sea Level Pressure | `Pa` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`sfc_boundary_layer`) |
| `'tauu'` | `id_tauu` | Surface Downward Eastward Wind Stress | `Pa` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`flux_down_from_atmos`) |
| `'tauv'` | `id_tauv` | Surface Downward Northward Wind Stress | `Pa` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`flux_down_from_atmos`) |
| `'hfss'` | `id_hfss` | Surface Upward Sensible Heat Flux | `W m-2` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`flux_up_to_atmos`) |
| `'hfls'` | `id_hfls` | Surface Upward Latent Heat Flux | `W m-2` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`flux_up_to_atmos`) |
| `'evspsbl'` | `id_evspsbl` | Evaporation | `kg m-2 s-1` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`flux_up_to_atmos`) |
| `'tslsi'` | `id_tslsi` | Surface Temperature Where Land or Sea Ice | `K` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`flux_up_to_atmos`) |
| `'tos'` | `id_tos` | Sea Surface Temperature | `K` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`flux_up_to_atmos`) |
| `'sic'` | `id_sic` | Sea Ice Area Fraction | `1.0` | `fms_diag_register_diag_field` / `register_cmip_diag_field_2d` | Yes (`diag_sic`) |

### `diag_field_init` — global-integral fields (`register_global_diag_field`)

These fields produce globally averaged scalar time-series output (separate from the
2-D grid diagnostics above). The `register_global_diag_field` call is only compiled in when
the `use_AM3_physics` macro is **not** defined.

| Field name | ID variable | Description | Units | send_data? |
|---|---|---|---|---|
| `'evspsbl'` | `id_evspsbl_g` | Evaporation | `mm d-1` | Yes (`flux_up_to_atmos`) |
| `'ts'` | `id_ts_g` | Surface Temperature | `K` | Yes (`flux_up_to_atmos`) |
| `'tas'` | `id_tas_g` | Near-Surface Air Temperature | `K` | Yes (`sfc_boundary_layer`) |
| `'tasl'` | `id_tasl_g` | Near-Surface Air Temperature (Land Only) | `K` | Yes (`sfc_boundary_layer`) |
| `'hfss'` | `id_hfss_g` | Surface Upward Sensible Heat Flux | `W m-2` | Yes (`flux_up_to_atmos`) |
| `'hfls'` | `id_hfls_g` | Surface Upward Latent Heat Flux | `W m-2` | Yes (`flux_up_to_atmos`) |
| `'rls'` | `id_rls_g` | Net Longwave Surface Radiation | `W m-2` | Yes (`flux_up_to_atmos`) |

### `diag_field_init` — tiled land fields (`register_tiled_diag_field` / `fms_diag_register_diag_field` on land axes)

These fields are registered only when `land_pe` is `.true.`. Under the `_USE_LEGACY_LAND_`
preprocessor branch the calls use `fms_diag_register_diag_field`; otherwise they use
`register_tiled_diag_field`. Field names, descriptions, and units are the same in both paths.
The module name is `'flux_land'` except where noted.

| Field name | ID variable | Description | Units | send_data? |
|---|---|---|---|---|
| `'t_ref'` (mod `flux_land`) | `id_t_ref_land` | temperature at reference height over land | `deg_k` | Yes (`sfc_boundary_layer`) |
| `'q_ref'` (mod `flux_land`) | `id_q_ref_land` | specific humidity at reference height over land | `kg/kg` | Yes (`sfc_boundary_layer`) |
| `'rh_ref'` (mod `flux_land`) | `id_rh_ref_land` | relative humidity at reference height over land | `percent` | Yes (`sfc_boundary_layer`) |
| `'u_ref'` (mod `flux_land`) | `id_u_ref_land` | zonal wind component at reference height over land | `m/s` | Yes (`sfc_boundary_layer`) |
| `'v_ref'` (mod `flux_land`) | `id_v_ref_land` | meridional wind component at reference height over land | `m/s` | Yes (`sfc_boundary_layer`) |
| `'evap'` (mod `flux_land`) | `id_q_flux_land` | evaporation rate over land | `kg/m2/s` | Yes (`flux_up_to_atmos`) |
| `'shflx'` (mod `flux_land`) | `id_t_flux_land` | sensible heat flux | `W/m2` | Yes (`flux_up_to_atmos`) |
| `'tasLut'` (mod `cmor_land`) | `id_tasLut_land` | Near-Surface Air Temperature (reference height above displacement height) on Land Use Tile | `K` | Yes (`sfc_boundary_layer`) |
| `'hussLut'` (mod `cmor_land`) | `id_hussLut_land` | Near-Surface Specific Humidity on Land Use Tile | `1.0` | Yes (`sfc_boundary_layer`) |

#### Per-tracer tiled land fields (loop over `n_exch_tr`)

| Field name | ID variable | Description | Units | send_data? |
|---|---|---|---|---|
| `trim(name)//'_tot_con_atm'` (mod `flux_land`) | `id_tr_con_atm_land(tr)` | vd of `longname` *(new-land path only)* | `m/s` | Yes (`flux_up_to_atmos`) |
| `trim(name)//'_tot_con_ref'` (mod `flux_land`) | `id_tr_con_ref_land(tr)` | vd of `longname` at reference height *(new-land path only)* | `m/s` | Yes (`flux_up_to_atmos`) |
| `trim(name)//'_flux'` (mod `flux_land`) | `id_tr_flux_land(tr)` | flux of `longname` | `units kg air/(m2 s)` | Yes (`flux_up_to_atmos`) |
| `trim(name)//'_mol_flux'` (mod `flux_land`) | `id_tr_mol_flux_land(tr)` | flux of `longname` | `mol CO2/(m2 s)` (CO2) or `mol/(m2 s)` (other) | Yes (`flux_up_to_atmos`) |
| `trim(name)//'_ref'` (mod `flux_land`) | `id_tr_ref_land(tr)` | `longname` at reference height over land *(skipped for sphum; new-land path only)* | `units` | Yes (`sfc_boundary_layer`) |

---

## `ice_ocean_flux_exchange.F90`

No `register*diag*` calls are present in this file. The file does contain two calls to
`fms_coupler_type_send_data` (lines 452 and 572), which send data for coupler boundary-condition
type fields that are registered elsewhere (in the coupler infrastructure), not in this file
directly.

---

## `land_ice_flux_exchange.F90`

No `register*diag*` calls are present in this file.

---

## `flux_exchange.F90`

No `register*diag*` calls are present in this file.
