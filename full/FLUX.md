# Flux Exchange

## Configuration

`flux_exchange_mod` is configured through the `flux_exchange_nml` namelist in `input.nml`.

| Variable Name | Type | Default Value | Description |
| --- | --- | --- | --- |
| `z_ref_heat` | `real` | `2.0` | Reference height in meters for temperature and relative humidity diagnostics (`t_ref`, `rh_ref`, `del_h`, `del_q`). |
| `z_ref_mom` | `real` | `10.0` | Reference height in meters for momentum diagnostics (`u_ref`, `v_ref`, `del_m`). |
| `do_area_weighted_flux` | `logical` | `.FALSE.` | Enables area-weighted flux handling. |
| `debug_stocks` | `logical` | `.FALSE.` | Enables additional stock-debug output. |
| `divert_stocks_report` | `logical` | `.FALSE.` | Diverts stock reporting output. |
| `do_runoff` | `logical` | `.TRUE.` | Turns land runoff interpolation to the ocean on or off. |
| `do_forecast` | `logical` | `.FALSE.` | Enables forecast-mode behavior in the flux coupler. |
| `nblocks` | `integer` | `1` | Number of blocks used to divide `n_xgrid_sfc`. This primarily supports OpenMP execution. In practice this is often set to match `coupler_nml%atmos_nthreads`. |
| `partition_fprec_from_lprec` | `logical` | `.FALSE.` | For atmosphere override experiments where liquid and frozen precipitation are combined, convert liquid precipitation to snow when `t_ref < tfreeze`. |
| `scale_precip_2d` | `logical` | `.FALSE.` | Rescale `Atm%lprec` using a field read from `data_table`. |

## Module Overview

Authors:
- Bruce Wyman <Bruce.Wyman@noaa.gov>
- V. Balaji <V.Balaji@noaa.gov>
- Sergey Malyshev <Sergey.Malyshev@noaa.gov>

The `flux_exchange_mod` module provides the interfaces used to couple atmosphere, ocean, land, and ice components. Interpolation between physically distinct model grids is handled by the exchange grid (`xgrid_mod`) with conservation of the interpolated quantities.

### Coupling assumptions and behavior

1. `flux_exchange_mod` supports physically independent atmosphere, land, and sea-ice grids. Ice and ocean must share the same physical grid, although their domain decompositions may differ.
2. Grid information is read from the grid specification file. The masked region of the land grid and the ice-ocean grid must tile each other, and the masked regions of the ice and ocean grids must be identical.
3. The atmosphere, land, and ice grids exchange information using the surface exchange grid `xmap_sfc`.
4. The land and ice grids exchange runoff data using the exchange grid `xmap_runoff`.
5. Ice-bottom to ocean transfer does not require an exchange grid because those grids are physically identical. The flux routines automatically redistribute data when decompositions differ.
6. Information from the atmosphere reaches the ocean through the ice model: first atmosphere to ice, then ice to ocean.
7. Each component model must expose a public data type containing the boundary fields needed by the coupler.
8. Sensible heat flux and surface evaporation can depend implicitly on surface temperature, so land and sea-ice temperature updates must run on the atmospheric time step.
9. Surface fluxes for all other tracers and for momentum are treated as explicit functions of the surface state.
10. The module is designed to support simultaneous implicit time integration on both sides of the surface interface.
11. Because of that implicit coupling, the diffusion part of the land and ice models must also run on the atmospheric time step.
12. Additional tracer and gas-exchange fluxes are configured through `field_table` and named boundary-condition fields in the coupler boundary types.
13. Any field exchanged between components can be replaced by a constant or file-based value using the FMS `data_override` facility configured through `data_table`.

The original documentation strongly advises against using the data override capabilities until the model configuration is well understood.

### Grid layout

```text
        ATMOSPHERE  |----|----|----|----|----|----|----|----|

              LAND  |---|---|---|---|xxx|xxx|xxx|xxx|xxx|xxx|

               ICE  |xxx|xxx|xxx|xxx|---|---|---|---|---|---|

              OCEAN |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
```

Here `|xxx|` marks a masked grid point.

## Data Override Capabilities

The module supports runtime data override in the following paths.

### Atmosphere boundary to exchange grid in `sfc_boundary_layer`

`t_bot`, `q_bot`, `z_bot`, `p_bot`, `u_bot`, `v_bot`, `p_surf`, `slp`, `gust`

### Ice boundary to exchange grid in `sfc_boundary_layer`

`t_surf`, `rough_mom`, `rough_heat`, `rough_moist`, `albedo`, `u_surf`, `v_surf`

### Land boundary to exchange grid in `sfc_boundary_layer`

`t_surf`, `t_ca`, `q_ca`, `rough_mom`, `rough_heat`, `albedo`

### Exchange grid to `land_ice_atmos_boundary` in `sfc_boundary_layer`

`t`, `albedo`, `land_frac`, `dt_t`, `dt_q`, `u_flux`, `v_flux`, `dtaudu`, `dtaudv`, `u_star`, `b_star`, `rough_mom`

### Atmosphere boundary to exchange grid in `flux_down_from_atmos`

`flux_sw`, `flux_lw`, `lprec`, `fprec`, `coszen`, `dtmass`, `delta_t`, `delta_q`, `dflux_t`, `dflux_q`

### Exchange grid to land boundary in `flux_down_from_atmos`

`t_flux`, `q_flux`, `lw_flux`, `sw_flux`, `lprec`, `fprec`, `dhdt`, `dedt`, `dedq`, `drdt`, `drag_q`, `p_surf`

### Exchange grid to ice boundary in `flux_down_from_atmos`

`u_flux`, `v_flux`, `t_flux`, `q_flux`, `lw_flux`, `lw_flux_dn`, `sw_flux`, `sw_flux_dn`, `lprec`, `fprec`, `dhdt`, `dedt`, `drdt`, `coszen`, `p`

### Land boundary to ice boundary in `flux_land_to_ice`

`runoff`, `calving`

### Ice boundary to ocean boundary in `flux_ice_to_ocean`

`u_flux`, `v_flux`, `t_flux`, `q_flux`, `salt_flux`, `lw_flux`, `sw_flux`, `lprec`, `fprec`, `runoff`, `calving`, `p`, `ustar_berg`, `area_berg`, `mass_berg`

### Ocean boundary to ice boundary in `flux_ocean_to_ice`

`u`, `v`, `t`, `s`, `frazil`, `sea_level`

### Ice boundary to atmosphere boundary in `flux_up_to_atmos`

`t_surf`

### Land boundary to atmosphere boundary in `flux_up_to_atmos`

`t_ca`, `t_surf`, `q_ca`

## Diagnostic Fields

The `flux` diagnostic module provides the following fields.

| Field Name | Units | Description |
| --- | --- | --- |
| `land_mask` | `none` | Fractional amount of land |
| `wind` | `m/s` | Wind speed for flux calculations |
| `drag_moist` | `none` | Drag coefficient for moisture |
| `drag_heat` | `none` | Drag coefficient for heat |
| `drag_mom` | `none` | Drag coefficient for momentum |
| `rough_moist` | `m` | Surface roughness for moisture |
| `rough_heat` | `m` | Surface roughness for heat |
| `rough_mom` | `m` | Surface roughness for momentum |
| `u_star` | `m/s` | Friction velocity |
| `b_star` | `m/s` | Buoyancy scale |
| `q_star` | `kg water/kg air` | Moisture scale |
| `t_atm` | `deg_k` | Temperature at bottom level |
| `u_atm` | `m/s` | Zonal wind component at bottom level |
| `v_atm` | `m/s` | Meridional wind component at bottom level |
| `q_atm` | `kg/kg` | Specific humidity at bottom level |
| `p_atm` | `pa` | Pressure at bottom level |
| `z_atm` | `m` | Height of bottom level |
| `gust` | `m/s` | Gust scale |
| `rh_ref` | `percent` | Relative humidity at reference height |
| `t_ref` | `deg_k` | Temperature at reference height |
| `u_ref` | `m/s` | Zonal wind component at reference height |
| `v_ref` | `m/s` | Meridional wind component at reference height |
| `del_h` | `none` | Reference-height interpolation factor for heat |
| `del_m` | `none` | Reference-height interpolation factor for momentum |
| `del_q` | `none` | Reference-height interpolation factor for moisture |
| `tau_x` | `pa` | Zonal wind stress |
| `tau_y` | `pa` | Meridional wind stress |
| `ice_mask` | `none` | Fractional amount of sea ice |
| `t_surf` | `deg_k` | Surface temperature |
| `t_ca` | `deg_k` | Canopy air temperature |
| `q_surf` | `kg/kg` | Surface specific humidity |
| `shflx` | `w/m2` | Sensible heat flux |
| `evap` | `kg/m2/s` | Evaporation rate |
| `lwflx` | `w/m2` | Net downward-minus-upward longwave flux |

## Main Program Example

The main coupling loop follows this pattern.

```fortran
DO slow time steps (ocean)
   call flux_ocean_to_ice

   call ICE_SLOW_UP

   DO fast time steps (atmos)
      call sfc_boundary_layer

      call ATMOS_DOWN

      call flux_down_from_atmos

      call LAND_FAST

      call ICE_FAST

      call flux_up_to_atmos

      call ATMOS_UP
   END DO

   call ICE_SLOW_DN

   call flux_ice_to_ocean

   call OCEAN
END DO
```

`LAND_FAST` and `ICE_FAST` must update the surface temperature.

## Required Variables in Defined Data Types for Component Models

### Atmosphere

```fortran
type (atmos_boundary_data_type) :: Atm

real, dimension(:) :: Atm%lon_bnd &
                       Atm%lat_bnd
real, dimension(:,:) :: Atm%t_bot   &
                        Atm%q_bot   &
                        Atm%z_bot   &
                        Atm%p_bot   &
                        Atm%u_bot   &
                        Atm%v_bot   &
                        Atm%p_surf  &
                        Atm%slp     &
                        Atm%gust    &
                        Atm%flux_sw &
                        Atm%flux_lw &
                        Atm%lprec   &
                        Atm%fprec   &
                        Atm%coszen
integer, dimension(4) :: Atm%axes
```

- `Atm%lon_bnd`, `Atm%lat_bnd`: Grid-box boundaries in radians and must be monotonic.
- `Atm%t_bot`, `Atm%q_bot`, `Atm%z_bot`, `Atm%p_bot`, `Atm%u_bot`, `Atm%v_bot`: State at the lowest model level.
- `Atm%p_surf`, `Atm%slp`, `Atm%gust`: Surface pressure, sea-level pressure, and gustiness factor.
- `Atm%flux_sw`, `Atm%flux_lw`, `Atm%lprec`, `Atm%fprec`, `Atm%coszen`: Surface radiative and precipitation inputs.
- `Atm%axes`: Axis identifiers returned by `diag_axis_init` for the atmospheric X, Y, `Z_full`, and `Z_half` axes.

The following fields are gathered for convenience in simultaneous implicit time stepping between the atmosphere and surface models. The original documentation points to `flux_exchange.tech.ps` and `vert_diff_mod` for more detail.

```fortran
type (surf_diff_type) :: Atm%Surf_Diff

real, dimension(:,:) :: Atm%Surf_Diff%dtmass  &
                        Atm%Surf_Diff%delta_t &
                        Atm%Surf_Diff%delta_q &
                        Atm%Surf_Diff%dflux_t &
                        Atm%Surf_Diff%dflux_q
```

- `dtmass`: `dt / mass`, where `dt` is the atmospheric time step.
- `delta_t`, `delta_q`: Leapfrog increments for lowest-layer temperature and specific humidity.
- `dflux_t`, `dflux_q`: Derivatives of the implicit downward flux terms with respect to lowest-layer temperature and humidity.

### Land

```fortran
type (land_boundary_data_type) :: Land

real, dimension(:) :: Land%lon_bnd &
                      Land%lat_bnd

logical, dimension(:,:,:) :: Land%mask    &
                             Land%glacier

real, dimension(:,:,:) :: Land%tile_size  &
                          Land%t_surf     &
                          Land%albedo     &
                          Land%rough_mom  &
                          Land%rough_heat &
                          Land%stomatal   &
                          Land%snow       &
                          Land%water      &
                          Land%max_water
```

- `Land%lon_bnd`, `Land%lat_bnd`: Grid-box boundaries in radians and must be monotonic.
- `Land%mask`: Land-sea mask, true over land.
- `Land%glacier`: Glacier mask, true over glacier.
- `Land%tile_size`: Fractional area of each land tile.
- `Land%t_surf`, `Land%albedo`, `Land%rough_mom`, `Land%rough_heat`: Surface state needed by the coupler.
- `Land%stomatal`, `Land%snow`, `Land%water`, `Land%max_water`: Additional land-surface properties.

### Ice

```fortran
type (ice_boundary_data_type) :: Ice

real, dimension(:) :: Ice%lon_bnd    &
                      Ice%lat_bnd    &
                      Ice%lon_bnd_uv &
                      Ice%lat_bnd_uv

logical, dimension(:,:,:) :: Ice%mask    &
                             Ice%mask_uv &
                             Ice%ice_mask

real, dimension(:,:,:) :: Ice%part_size    &
                          Ice%part_size_uv
```

- Boundary arrays are in radians and must be monotonic.
- `Ice%mask` and `Ice%mask_uv` are ocean-land masks for temperature and momentum points.
- `Ice%ice_mask` is an optional explicit ice mask.
- `Ice%part_size` and `Ice%part_size_uv` are fractional partition areas.

Fields on the ice top grid:

```fortran
real, dimension(:,:,:) :: Ice%t_surf     &
                          Ice%albedo     &
                          Ice%rough_mom  &
                          Ice%rough_heat &
                          Ice%u_surf     &
                          Ice%v_surf
```

Fields on the ice bottom grid:

```fortran
real, dimension(:,:,:) :: Ice%flux_u  &
                          Ice%flux_v  &
                          Ice%flux_t  &
                          Ice%flux_q  &
                          Ice%flux_sw &
                          Ice%flux_lw &
                          Ice%lprec   &
                          Ice%fprec   &
                          Ice%runoff
```

### Ocean

```fortran
type (ocean_boundary_data_type) :: Ocean

real, dimension(:) :: Ocean%Data%lon_bnd     &
                      Ocean%Data%lat_bnd     &
                      Ocean%Data%lon_bnd_uv  &
                      Ocean%Data%lat_bnd_uv  &
                      Ocean%Ocean%lon_bnd    &
                      Ocean%Ocean%lat_bnd    &
                      Ocean%Ocean%lon_bnd_uv &
                      Ocean%Ocean%lat_bnd_uv
```

All longitude and latitude boundary arrays must be monotonic.

```fortran
logical, dimension(:,:) :: Ocean%Data%mask    &
                           Ocean%Data%mask_uv &
                           Ocean%Ocean%mask   &
                           Ocean%Ocean%mask_uv
real, dimension(:,:) :: Ocean%t_surf_data &
                        Ocean%t_surf      &
                        Ocean%u_surf      &
                        Ocean%v_surf      &
                        Ocean%frazil
```

- `Ocean%Data%mask`, `Ocean%Data%mask_uv`: Ocean-land masks on the ocean data grid.
- `Ocean%Ocean%mask`, `Ocean%Ocean%mask_uv`: Ocean-land masks on the ocean model grid.
- `Ocean%t_surf_data`, `Ocean%t_surf`, `Ocean%u_surf`, `Ocean%v_surf`, `Ocean%frazil`: Surface state on the ocean data and model grids.
