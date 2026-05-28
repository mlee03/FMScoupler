# Flux Exchange

Authors:
- Bruce Wyman <Bruce.Wyman@noaa.gov>
- V. Balaji <V.Balaji@noaa.gov>
- Sergey Malyshev <Sergey.Malyshev@noaa.gov>

## Overview
There are six modules to couple the atmosphere, ocean, land, and ice components
through flux exchange:
* atm_land_ice_flux_exchange:  exchange fluxes between atm, land, and ice 
* atmos_ocean_fluxes_calc:  compute non-deposition gas fluxes between atm and ocean 
* atmos_ocean_dep_fluxes_calc:  compute deposition gas fluxes between atm and ocean
* ice_ocean_flux_exchange:  exchange fluxes between ice and ocean 
* land_ice_flux_exchange: exchange fluxes between land and ice
* flux_exchange:  top level module that initializes the various flux_exchange module; 
                  also contains subroutines for stock computation

Note the following:
1. Flux exchange in coupler supports physically independent atmosphere, land, and sea-ice grids. 
   However, ice and ocean must share the same physical grid, although their domain decompositions may differ.
2. The masked region of the land grid and the ice-ocean grid must tile each other.
3. The masked regions of the ice and ocean grids must be identical.
4. The atmosphere, land, and ice grids exchange information using the surface exchange grid `xmap_sfc`
   with conservative interpolation
5. The land and ice grids exchange runoff data using the exchange grid `xmap_runoff`.
6. Ice-bottom to ocean transfer does not require an exchange grid because those grids are physically identical. 
   The flux data are automatically redistributed when decompositions differ (REDIST=2)
7. Information from the atmosphere reaches the ocean through the ice model: first atmosphere to ice, then ice to ocean.
8. Each component model must expose a public data type containing the boundary fields needed by the coupler.
9. Sensible heat flux and surface evaporation can depend implicitly on surface temperature, 
   so land and sea-ice temperature updates must run on the atmospheric time step.
10. Surface fluxes for all other tracers and for momentum are treated as explicit functions of the surface state.
11. The module is designed to support simultaneous implicit time integration on both sides of the surface interface.
12. Because of that implicit coupling, the diffusion part of the land and ice models must also run on the 
    atmospheric time step.
13. Additional tracer and gas-exchange fluxes are configured through `field_table` and named boundary-condition 
    fields in the coupler boundary types.
14. Any field exchanged between components can be replaced by a constant or file-based value using the FMS `data_override` 
    facility configured through `data_table`.


## Configuration
The below can be configured with the flux_exchange_nml in input.nml
* `z_ref_heat` (real, default = 2.0):  reference height in meters for temperature and relative humidity diagnostics (t_ref, rh_ref, del_h, del_q)
* `z_ref_mom` (real, default 10.0):  reference height in meters for momentum diagnostics (u_ref, v_ref, del_m)
* `do_area_weighted_flux` (logical, default = .FALSE.): enables area-weighted flux handling. When .TRUE., fluxes 
  passed to the ocean are multiplied by the ice area fraction before redistribution, so the ocean receives the
  grid-cell-mean flux rather than the ice-covered-area flux.
* `debug_stocks` (logical, default = .FALSE.): enables additional stock-debug output.
* `divert_stocks_report` (logical, default = .FALSE.): diverts stock reporting output 
  to `stocks.out` rather than the standard log.
* `do_runoff` (logical, default = .TRUE.): turns land runoff interpolation to the ocean on or off.
* `do_forecast` (logical, default = .FALSE.):  enables forecast-mode behavior in the flux coupler. 
* `nblocks` (integer, default = 1) Number of blocks used to divide n_xgrid_sfc for OpenMP parallelism. 
  In practice this is often set to match `coupler_nml%atmos_nthreads`
* `partition_fprec_from_lprec` (logical): default = .FALSE,  For atmosphere override experiments where 
  liquid and frozen precipitation are combined, convert liquid precipitation to snow when t_ref < tfreeze
* `scale_precip_2d` (logical, default = .FALSE.):  rescale `Atm%lprec` using a field read from data_table




### Grid layout
```
        ATMOSPHERE  |----|----|----|----|----|----|----|----|
              LAND  |---|---|---|---|xxx|xxx|xxx|xxx|xxx|xxx|
               ICE  |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
             OCEAN  |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
```
Here `|xxx|` marks a masked grid point.

## Data Override Capabilities
The original documentation strongly advises against using the data override capabilities 
until the model configuration is well understood.  The module supports runtime data override in the following paths.
* Atmosphere boundary to exchange grid in `sfc_boundary_layer`:
  - `t_bot`, `q_bot`, `z_bot`, `p_bot`, `u_bot`, `v_bot`, `p_surf`, `slp`, `gust`

* Ice boundary to exchange grid in `sfc_boundary_layer`:
  - `t_surf`, `rough_mom`, `rough_heat`, `rough_moist`, `albedo`, `u_surf`, `v_surf`

* Land boundary to exchange grid in `sfc_boundary_layer`:
  - `t_surf`, `t_ca`, `q_ca`, `rough_mom`, `rough_heat`, `albedo`

* Exchange grid to `land_ice_atmos_boundary` in `sfc_boundary_layer`:
  - `t`, `t_ocean`, `frac_open_sea`, `albedo`, `albedo_vis_dir`, `albedo_nir_dir`, `albedo_vis_dif`, `albedo_nir_dif`, `land_frac`, `rough_mom`, `rough_heat`, `u_flux`, `v_flux`, `dtaudu`, `dtaudv`, `u_star`, `b_star`, `q_star`, `u_ref`, `v_ref`, `wind`

* Atmosphere boundary to exchange grid in `flux_down_from_atmos`:
  - `flux_sw`, `flux_lw`, `lprec`, `fprec`, `coszen`, `dtmass`, `delta_t`, `delta_q`, `dflux_t`, `dflux_q`

* Exchange grid to land boundary in `flux_down_from_atmos`:
  - `t_flux`, `lw_flux`, `lwdn_flux`, `sw_flux`, `sw_flux_down_vis_dir`, `sw_flux_down_total_dir`, `sw_flux_down_vis_dif`, `sw_flux_down_total_dif`, `lprec`, `fprec`, `tprec`, `dhdt`, `dedt`, `dedq`, `drdt`, `drag_q`, `p_surf`, `tr_flux`, `dfdtr`

* Exchange grid to ice boundary in `flux_down_from_atmos`:
  - `u_flux`, `v_flux`, `t_flux`, `q_flux`, `lw_flux`, `sw_flux_vis_dir`, `sw_flux_nir_dir`, `sw_flux_vis_dif`, `sw_flux_nir_dif`, `sw_down_vis_dir`, `sw_down_nir_dir`, `sw_down_vis_dif`, `sw_down_nir_dif`, `lprec`, `fprec`, `dhdt`, `dedt`, `drdt`, `u_star`, `coszen`, `p`, `fluxes`

* Land boundary to ice boundary in `flux_land_to_ice`:
  - `runoff`, `calving`, `runoff_hflx`, `calving_hflx`

* Ice boundary to ocean boundary in `flux_ice_to_ocean`:
  - `u_flux`, `v_flux`, `t_flux`, `q_flux`, `salt_flux`, `lw_flux`, `sw_flux_vis_dir`, `sw_flux_vis_dif`, `sw_flux_nir_dir`, `sw_flux_nir_dif`, `lprec`, `fprec`, `runoff`, `calving`, `runoff_hflx`, `calving_hflx`, `p`, `mi`, `ustar_berg`\*, `area_berg`\*, `mass_berg`\*`, `fluxes`

\* Allocated and passed only when the corresponding pointer is associated in the `Ice` data type (i.e., when the iceberg model is active).

### Ocean boundary to ice boundary in `flux_ocean_to_ice`

`u`, `v`, `t`, `s`, `frazil`, `sea_level`, `fields`

### Ice boundary to atmosphere boundary in `flux_up_to_atmos`

`dt_t`, `dt_tr`

From ice surface: `t_surf`

From land surface: `t_ca`, `t_surf`, `q_ca`

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

The main coupling loop follows this pattern. See `coupler_main.F90` for the full annotated pseudocode.

```fortran
! Initialization
call atmos_model_init / land_model_init / ice_model_init / ocean_model_init
call flux_exchange_init       ! build atm-land-ice exchange grids

do nc = 1, num_cpld_calls     ! slow (ocean/ice) loop

   call flux_ocean_to_ice     ! ocean SST, currents, frazil → ice

   call exchange_slow_to_fast_ice
   call set_ice_surface_fields

   do na = 1, num_atmos_calls ! fast (atmos/land/ice) loop

      call sfc_boundary_layer         ! turbulent fluxes, exchange coefficients

      call update_atmos_model_dynamics
      call update_atmos_model_radiation
      call update_atmos_model_down    ! downward tridiagonal sweep

      call flux_down_from_atmos       ! SW/LW/precip/implicit-diff → land, ice

      call update_land_model_fast     ! land surface temperature update
      call update_ice_model_fast      ! ice surface temperature update

      call flux_up_to_atmos           ! updated t_surf → atmosphere

      call update_atmos_model_up      ! back-substitution, convection

      call flux_atmos_to_ocean        ! deposition gas fluxes → ice/ocean

   end do

   call update_land_model_slow
   call flux_land_to_ice             ! runoff, calving, heat fluxes → ice

   call exchange_fast_to_slow_ice
   call update_ice_model_slow        ! ice dynamics, freeze/melt

   call flux_ice_to_ocean            ! all ice-bottom fluxes → ocean
   call update_ocean_model

end do
```

`update_land_model_fast` and `update_ice_model_fast` must update the surface temperature each atmospheric time step for the implicit diffusion scheme to be correct.

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

- `Atm%lon_bnd`, `Atm%lat_bnd`: Grid-box boundaries in radians; must be monotonic.
- `Atm%t_bot`, `Atm%q_bot`, `Atm%z_bot`, `Atm%p_bot`, `Atm%u_bot`, `Atm%v_bot`: State at the lowest model level.
- `Atm%p_surf`, `Atm%slp`, `Atm%gust`: Surface pressure, sea-level pressure, and gustiness factor.
- `Atm%flux_sw`, `Atm%flux_lw`, `Atm%lprec`, `Atm%fprec`, `Atm%coszen`: Surface radiative and precipitation inputs.
- `Atm%axes`: Axis identifiers returned by `diag_axis_init` for the atmospheric X, Y, `Z_full`, and `Z_half` axes.

The following fields support the simultaneous implicit time-stepping between the atmosphere and surface models. See `flux_exchange.tech.ps` and `vert_diff_mod` for details.

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
                          Land%t_ca       &
                          Land%q_ca       &
                          Land%albedo     &
                          Land%rough_mom  &
                          Land%rough_heat &
                          Land%stomatal   &
                          Land%snow       &
                          Land%water      &
                          Land%max_water
```

- `Land%lon_bnd`, `Land%lat_bnd`: Grid-box boundaries in radians; must be monotonic.
- `Land%mask`: Land-sea mask, true over land.
- `Land%glacier`: Glacier mask, true over glacier.
- `Land%tile_size`: Fractional area of each land tile.
- `Land%t_surf`, `Land%albedo`, `Land%rough_mom`, `Land%rough_heat`: Surface state used by the coupler for turbulent flux and radiation calculations.
- `Land%t_ca`, `Land%q_ca`: Canopy air temperature and specific humidity; returned to the atmosphere by `flux_up_to_atmos`.
- `Land%stomatal`, `Land%snow`, `Land%water`, `Land%max_water`: Additional land-surface properties used in surface flux parameterizations.

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

Fields on the ice **top** grid (atmosphere–ice interface):

```fortran
real, dimension(:,:,:) :: Ice%t_surf      &
                          Ice%albedo      &
                          Ice%rough_mom   &
                          Ice%rough_heat  &
                          Ice%rough_moist &
                          Ice%u_surf      &
                          Ice%v_surf
```

Fields on the ice **bottom** grid (ice–ocean interface), populated by `flux_down_from_atmos` and `flux_land_to_ice`, then passed to the ocean by `flux_ice_to_ocean`:

```fortran
real, dimension(:,:,:) :: Ice%flux_u          &  ! zonal wind stress [Pa]
                          Ice%flux_v          &  ! meridional wind stress [Pa]
                          Ice%flux_t          &  ! sensible heat flux [W/m2]
                          Ice%flux_q          &  ! moisture flux [kg/m2/s]
                          Ice%flux_salt       &  ! salt flux [kg/m2/s]
                          Ice%flux_lw         &  ! net longwave flux [W/m2]
                          Ice%flux_sw_vis_dir &  ! direct visible SW [W/m2]
                          Ice%flux_sw_vis_dif &  ! diffuse visible SW [W/m2]
                          Ice%flux_sw_nir_dir &  ! direct near-IR SW [W/m2]
                          Ice%flux_sw_nir_dif &  ! diffuse near-IR SW [W/m2]
                          Ice%lprec           &  ! liquid precipitation [kg/m2/s]
                          Ice%fprec           &  ! frozen precipitation [kg/m2/s]
                          Ice%runoff          &  ! liquid runoff from land [kg/m2/s]
                          Ice%calving         &  ! solid discharge from land [kg/m2/s]
                          Ice%runoff_hflx     &  ! heat flux with runoff [W/m2]
                          Ice%calving_hflx    &  ! heat flux with calving [W/m2]
                          Ice%p_surf             ! atmospheric surface pressure [Pa]
```

Optional iceberg fields (allocated only when the iceberg model is active):

```fortran
real, dimension(:,:) :: Ice%ustar_berg  &  ! iceberg friction velocity [m/s]
                        Ice%area_berg   &  ! iceberg area fraction
                        Ice%mass_berg      ! iceberg mass [kg/m2]
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
                        Ocean%frazil      &
                        Ocean%s
```

- `Ocean%Data%mask`, `Ocean%Data%mask_uv`: Ocean-land masks on the ocean data grid.
- `Ocean%Ocean%mask`, `Ocean%Ocean%mask_uv`: Ocean-land masks on the ocean model grid.
- `Ocean%t_surf_data`, `Ocean%t_surf`: SST on the data and model grids; passed to ice as `Ocean_ice_boundary%t`.
- `Ocean%u_surf`, `Ocean%v_surf`: Surface ocean currents; passed to ice as `Ocean_ice_boundary%u`, `%v`.
- `Ocean%frazil`: Frazil ice heat flux [J/m²]; passed to ice as `Ocean_ice_boundary%frazil`.
- `Ocean%s`: Surface salinity [psu]; passed to ice as `Ocean_ice_boundary%s`.
