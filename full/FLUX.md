# Flux Exchange

Authors:
- Bruce Wyman <Bruce.Wyman@noaa.gov>
- V. Balaji <V.Balaji@noaa.gov>
- Sergey Malyshev <Sergey.Malyshev@noaa.gov>

## Overview
There are six modules to couple the atmosphere, ocean, land, and ice components
through flux exchange:
* `atm_land_ice_flux_exchange`:  exchange fluxes between atm, land, and ice 
* `atmos_ocean_fluxes_calc`:  compute non-deposition gas fluxes between atm and ocean 
* `atmos_ocean_dep_fluxes_calc`:  compute deposition gas fluxes between atm and ocean
* `ice_ocean_flux_exchange`:  exchange fluxes between ice and ocean 
* `land_ice_flux_exchange`: exchange fluxes between land and ice
* `flux_exchange`:  top level module that initializes the various flux_exchange module; 
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
15. `update_land_model_fast` and `update_ice_model_fast` must update the surface temperature each atmospheric 
    time step for the implicit diffusion scheme to be correct.

## Configuration
The below can be configured with the flux_exchange_nml in input.nml
* `z_ref_heat` (real, default = 2.0):  reference height in meters for temperature and relative humidity diagnostics 
  (t_ref, rh_ref, del_h, del_q)
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
|---|---|
| `Atm%t_bot` | Temperature at the lowest atmospheric level [K] |
| `Atm%z_bot` | Height of the lowest atmospheric level [m] |
| `Atm%p_bot` | Pressure at the lowest atmospheric level [Pa] |
| `Atm%u_bot` | Zonal wind at the lowest atmospheric level [m/s] |
| `Atm%v_bot` | Meridional wind at the lowest atmospheric level [m/s] |
| `Atm%p_surf` | Surface pressure [Pa] |
| `Atm%slp` | Sea-level pressure [Pa] |
| `Atm%gust` | Gustiness velocity used to augment surface wind speed in flux calculations [m/s] |
| `atm%fields%bc(n)%field(m)%values` | Per-tracer atmospheric surface fields (e.g. tracer concentrations at the lowest model level) |

* Ice boundary to exchange grid in `sfc_boundary_layer`:
|---|---|
| `Ice%t_surf` | Ice/ocean surface skin temperature [K] |
| `Ice%rough_mom` | Surface roughness length for momentum over ice [m] |
| `Ice%rough_heat` | Surface roughness length for heat over ice [m] |
| `Ice%rough_moist` | Surface roughness length for moisture over ice [m] |
| `Ice%albedo` | Broadband surface albedo over ice [dimensionless] |
| `Ice%albedo_vis_dir` | Direct-beam visible-band albedo over ice [dimensionless] |
| `Ice%albedo_nir_dir` | Direct-beam near-infrared albedo over ice [dimensionless] |
| `Ice%albedo_vis_dif` | Diffuse visible-band albedo over ice [dimensionless] |
| `Ice%albedo_nir_dif` | Diffuse near-infrared albedo over ice [dimensionless] |
| `Ice%u_surf` | Zonal surface current velocity of ice/ocean [m/s] |
| `Ice%v_surf` | Meridional surface current velocity of ice/ocean [m/s] |
| `Ice%ocean_fields` | Coupler boundary-condition type holding ocean/ice-top gas and tracer fields used in atmosphere-ocean flux calculations |

* Land boundary to exchange grid in `sfc_boundary_layer`:
|---|---|
| `Land%t_surf` | Land surface (radiative) temperature [K] |
| `Land%t_ca` | Canopy air temperature — near-surface air temperature within the plant canopy [K] |
| `Land%rough_mom` | Surface roughness length for momentum over land [m] |
| `Land%rough_heat` | Surface roughness length for heat over land [m] |
| `Land%albedo` | Broadband surface albedo over land [dimensionless] |
| `Land%tr(:,:,:,tr)` | Surface tracer mixing ratio over land; one entry per exchanged land tracer |
| `Land%albedo_vis_dir` | Direct-beam visible-band albedo over land [dimensionless] |
| `Land%albedo_nir_dir` | Direct-beam near-infrared albedo over land [dimensionless] |
| `Land%albedo_vis_dif` | Diffuse visible-band albedo over land [dimensionless] |
| `Land%albedo_nir_dif` | Diffuse near-infrared albedo over land [dimensionless] |

* Exchange grid to `land_ice_atmos_boundary` in `sfc_boundary_layer`:
|---|---|
| `Land_Ice_Atmos_Boundary%t` | Surface temperature (area-weighted over land and ice fractions) seen by the atmosphere [K] |
| `Land_Ice_Atmos_Boundary%albedo` | Broadband surface albedo [dimensionless] |
| `Land_Ice_Atmos_Boundary%albedo_vis_dir` | Direct-beam visible-band surface albedo [dimensionless] |
| `Land_Ice_Atmos_Boundary%albedo_nir_dir` | Direct-beam near-infrared surface albedo [dimensionless] |
| `Land_Ice_Atmos_Boundary%albedo_vis_dif` | Diffuse visible-band surface albedo [dimensionless] |
| `Land_Ice_Atmos_Boundary%albedo_nir_dif` | Diffuse near-infrared surface albedo [dimensionless] |
| `Land_Ice_Atmos_Boundary%land_frac` | Fraction of atmospheric grid cell covered by land [dimensionless] |
| `Land_Ice_Atmos_Boundary%dt_t` | Implicit correction to atmospheric temperature from surface flux scheme [K] |
| `Land_Ice_Atmos_Boundary%dt_tr(:,:,tr)` | Implicit correction to each atmospheric tracer mixing ratio from surface flux scheme; one entry per exchanged tracer |
| `Land_Ice_Atmos_Boundary%u_flux` | Zonal wind stress on the atmosphere [Pa] |
| `Land_Ice_Atmos_Boundary%v_flux` | Meridional wind stress on the atmosphere [Pa] |
| `Land_Ice_Atmos_Boundary%dtaudu` | d(wind stress)/d(u) — implicit coupling coefficient for zonal momentum [Pa·s/m] |
| `Land_Ice_Atmos_Boundary%dtaudv` | d(wind stress)/d(v) — implicit coupling coefficient for meridional momentum [Pa·s/m] |
| `Land_Ice_Atmos_Boundary%u_star` | Friction velocity (surface turbulent velocity scale) [m/s] |
| `Land_Ice_Atmos_Boundary%b_star` | Buoyancy scale used in Monin-Obukhov similarity theory [m/s²] |
| `Land_Ice_Atmos_Boundary%rough_mom` | Roughness length for momentum [m] |

* Atmosphere boundary to exchange grid in `flux_down_from_atmos`:
|---|---|
| `Atm%flux_sw` | Total net shortwave flux at the surface [W/m²] |
| `Atm%flux_sw_dir` | Direct-beam component of net shortwave flux [W/m²] |
| `Atm%flux_sw_dif` | Diffuse component of net shortwave flux [W/m²] |
| `Atm%flux_sw_down_vis_dir` | Downward direct-beam visible shortwave flux [W/m²] |
| `Atm%flux_sw_down_vis_dif` | Downward diffuse visible shortwave flux [W/m²] |
| `Atm%flux_sw_down_total_dir` | Downward direct-beam total (broadband) shortwave flux [W/m²] |
| `Atm%flux_sw_down_total_dif` | Downward diffuse total (broadband) shortwave flux [W/m²] |
| `Atm%flux_sw_vis` | Net visible-band shortwave flux at the surface [W/m²] |
| `Atm%flux_sw_vis_dir` | Direct-beam component of net visible shortwave flux [W/m²] |
| `Atm%flux_sw_vis_dif` | Diffuse component of net visible shortwave flux [W/m²] |
| `Atm%flux_lw` | Net downward longwave flux at the surface [W/m²] |
| `Atm%lprec` | Liquid precipitation rate [kg/m²/s] |
| `frac_precip` | 2-D scaling field for liquid precipitation; applied when `scale_precip_2d=.true.` [dimensionless] |
| `Atm%fprec` | Frozen (solid) precipitation rate [kg/m²/s] |
| `Atm%coszen` | Cosine of the solar zenith angle [dimensionless] |
| `Atm%Surf_Diff%dtmass` | dt/mass — ratio of timestep to surface layer mass used in the implicit diffusion scheme [s·m²/kg] |
| `Atm%Surf_Diff%delta_t` | Forward-elimination temperature coefficient from the implicit vertical diffusion scheme [K] |
| `Atm%Surf_Diff%dflux_t` | d(sensible heat flux)/d(T_surf) — linearisation coefficient for the implicit heat flux scheme [W/m²/K] |
| `Atm%Surf_Diff%delta_tr(:,:,tr)` | Forward-elimination coefficient for each tracer from the implicit diffusion scheme |
| `Atm%Surf_Diff%dflux_tr(:,:,tr)` | d(tracer flux)/d(tracer_surf) — linearisation coefficient for the implicit tracer flux scheme |

* Exchange grid to land boundary in `flux_down_from_atmos`:
 |---|---|
| `Land_boundary%drag_q` | Drag coefficient for moisture used in land surface flux calculations [dimensionless] |
| `Land_boundary%lwdn_flux` | Downward longwave radiation flux at the land surface [W/m²] |
| `Land_boundary%cd_m` | Drag coefficient for momentum over land [dimensionless] |
| `Land_boundary%cd_t` | Drag coefficient for heat over land [dimensionless] |
| `Land_boundary%bstar` | Buoyancy scale passed to the land model [m/s²] |
| `Land_boundary%ustar` | Friction velocity passed to the land model [m/s] |
| `Land_boundary%wind` | Wind speed at the lowest atmospheric level for land surface calculations [m/s] |
| `Land_boundary%z_bot` | Height of the lowest atmospheric level above the land surface [m] |
| `Land_boundary%t_flux` | Sensible heat flux into the land surface [W/m²] |
| `Land_boundary%lw_flux` | Net longwave flux at the land surface [W/m²] |
| `Land_boundary%sw_flux` | Net shortwave flux at the land surface [W/m²] |
| `Land_boundary%sw_flux_down_vis_dir` | Downward direct-beam visible shortwave flux over land [W/m²] |
| `Land_boundary%sw_flux_down_total_dir` | Downward direct-beam total shortwave flux over land [W/m²] |
| `Land_boundary%sw_flux_down_vis_dif` | Downward diffuse visible shortwave flux over land [W/m²] |
| `Land_boundary%sw_flux_down_total_dif` | Downward diffuse total shortwave flux over land [W/m²] |
| `Land_boundary%lprec` | Liquid precipitation rate over land [kg/m²/s] |
| `Land_boundary%fprec` | Frozen precipitation rate over land [kg/m²/s] |
| `Land_boundary%dhdt` | d(sensible heat flux)/d(T_surf) — implicit coupling coefficient for land heat flux [W/m²/K] |
| `Land_boundary%drdt` | d(longwave flux)/d(T_surf) — implicit coupling coefficient for land longwave flux [W/m²/K] |
| `Land_boundary%p_surf` | Surface pressure over land [Pa] |
| `Land_boundary%tr_flux(:,:,:,tr)` | Flux of each exchanged tracer into the land surface [kg/m²/s] |
| `Land_boundary%dfdtr(:,:,:,tr)` | d(tracer flux)/d(tracer_surf) — implicit coupling coefficient for each land tracer flux |

* Exchange grid to ice boundary in `flux_down_from_atmos`:
|---|---|
| `Ice_boundary%u_flux` | Zonal wind stress on the ice surface [Pa] |
| `Ice_boundary%v_flux` | Meridional wind stress on the ice surface [Pa] |
| `Ice_boundary%t_flux` | Sensible heat flux into the ice surface [W/m²] |
| `Ice_boundary%q_flux` | Latent heat (moisture) flux into the ice surface [W/m²] |
| `Ice_boundary%lw_flux` | Net longwave flux at the ice surface [W/m²] |
| `Ice_boundary%lw_flux` | Downward longwave flux at the ice surface; alternate override target for the same field [W/m²] |
| `Ice_boundary%sw_flux_nir_dir` | Direct-beam near-infrared shortwave flux over ice [W/m²] |
| `Ice_boundary%sw_flux_vis_dir` | Direct-beam visible shortwave flux over ice [W/m²] |
| `Ice_boundary%sw_flux_nir_dif` | Diffuse near-infrared shortwave flux over ice [W/m²] |
| `Ice_boundary%sw_flux_vis_dif` | Diffuse visible shortwave flux over ice [W/m²] |
| `Ice_boundary%sw_down_vis_dir` | Downward direct-beam visible shortwave flux over ice [W/m²] |
| `Ice_boundary%sw_down_vis_dif` | Downward diffuse visible shortwave flux over ice [W/m²] |
| `Ice_boundary%sw_down_nir_dir` | Downward direct-beam near-infrared shortwave flux over ice [W/m²] |
| `Ice_boundary%sw_down_nir_dif` | Downward diffuse near-infrared shortwave flux over ice [W/m²] |
| `Ice_boundary%lprec` | Liquid precipitation rate over ice [kg/m²/s] |
| `Ice_boundary%fprec` | Frozen precipitation rate over ice [kg/m²/s] |
| `Ice_boundary%dhdt` | d(sensible heat flux)/d(T_surf) — implicit coupling coefficient for ice heat flux [W/m²/K] |
| `Ice_boundary%dedt` | d(latent heat flux)/d(T_surf) — implicit coupling coefficient for ice moisture flux [W/m²/K] |
| `Ice_boundary%drdt` | d(longwave flux)/d(T_surf) — implicit coupling coefficient for ice longwave flux [W/m²/K] |
| `Ice_boundary%coszen` | Cosine of the solar zenith angle over ice [dimensionless] |
| `Ice_boundary%p` | Surface pressure over ice [Pa] |
| `Ice_boundary%fluxes` | Coupler boundary-condition type holding all per-tracer gas and deposition fluxes from atmosphere to ice |

* Ice boundary to atmosphere boundary in `flux_up_to_atmos`
| `Ice%t_surf` | Updated ice surface temperature after the ice model step [K] |

* Land boundary to atmosphere boundary in `flux_up_to_atmos`
| `Land%t_ca` | Updated canopy air temperature after the land model step [K] |
| `Land%t_surf` | Updated land surface temperature after the land model step [K] |
| `Land%tr(:,:,:,tr)` | Updated surface tracer mixing ratio over land after the land model step; one entry per exchanged tracer |

* Land boundary to ice boundary in `flux_land_to_ice`
|---|---|
| `Land_Ice_Boundary%runoff` | Liquid runoff (river discharge) from land to ocean/ice [kg/m²/s] |
| `Land_Ice_Boundary%calving` | Solid calving flux (iceberg / glacier discharge) from land to ocean/ice [kg/m²/s] |
| `Land_Ice_Boundary%runoff_hflx` | Heat flux carried by liquid runoff [W/m²] |
| `Land_Ice_Boundary%calving_hflx` | Heat flux carried by calving discharge [W/m²] |

* Ice boundary to ocean boundary in `flux_ice_to_ocean`
|---|---|
| `Ice_Ocean_Boundary%u_flux` | Zonal wind/ice stress on the ocean surface [Pa] |
| `Ice_Ocean_Boundary%v_flux` | Meridional wind/ice stress on the ocean surface [Pa] |
| `Ice_Ocean_Boundary%t_flux` | Sensible heat flux into the ocean [W/m²] |
| `Ice_Ocean_Boundary%q_flux` | Latent heat (freshwater) flux into the ocean [W/m²] |
| `Ice_Ocean_Boundary%salt_flux` | Salt flux from sea-ice to ocean (brine rejection / melting) [kg/m²/s] |
| `Ice_Ocean_Boundary%lw_flux` | Net longwave flux at the ocean surface [W/m²] |
| `Ice_Ocean_Boundary%sw_flux_nir_dir` | Direct-beam near-infrared shortwave flux into the ocean [W/m²] |
| `Ice_Ocean_Boundary%sw_flux_nir_dif` | Diffuse near-infrared shortwave flux into the ocean [W/m²] |
| `Ice_Ocean_Boundary%sw_flux_vis_dir` | Direct-beam visible shortwave flux into the ocean [W/m²] |
| `Ice_Ocean_Boundary%sw_flux_vis_dif` | Diffuse visible shortwave flux into the ocean [W/m²] |
| `Ice_Ocean_Boundary%lprec` | Liquid precipitation flux reaching the ocean surface [kg/m²/s] |
| `Ice_Ocean_Boundary%fprec` | Frozen precipitation flux reaching the ocean surface [kg/m²/s] |
| `Ice_Ocean_Boundary%runoff` | Liquid runoff from land routed to the ocean [kg/m²/s] |
| `Ice_Ocean_Boundary%calving` | Solid calving flux routed to the ocean [kg/m²/s] |
| `Ice_Ocean_Boundary%runoff_hflx` | Heat flux carried by liquid runoff into the ocean [W/m²] |
| `Ice_Ocean_Boundary%calving_hflx` | Heat flux carried by calving into the ocean [W/m²] |
| `Ice_Ocean_Boundary%p` | Surface atmospheric pressure at the ocean surface [Pa] |
| `Ice_Ocean_Boundary%mi` | Sea-ice mass per unit area (used for pressure loading on the ocean) [kg/m²] |
| `Ice_Ocean_Boundary%ustar_berg` | Friction velocity beneath icebergs [m/s]; only present when iceberg module is active |
| `Ice_Ocean_Boundary%area_berg` | Fractional area of ocean cell covered by icebergs [dimensionless]; only present when iceberg module is active |
| `Ice_Ocean_Boundary%mass_berg` | Iceberg mass per unit area [kg/m²]; only present when iceberg module is active |
| `Ice_Ocean_Boundary%fluxes` | Coupler boundary-condition type holding all per-tracer gas fluxes from ice/atmosphere to ocean |

* Ocean boundary to ice boundary in `flux_ocean_to_ice`
|---|---|
| `Ocean_Ice_Boundary%u` | Zonal ocean surface current velocity [m/s] |
| `Ocean_Ice_Boundary%v` | Meridional ocean surface current velocity [m/s] |
| `Ocean_Ice_Boundary%t` | Sea-surface temperature seen by the ice model [K] |
| `Ocean_Ice_Boundary%s` | Sea-surface salinity seen by the ice model [psu] |
| `Ocean_Ice_Boundary%frazil` | Frazil ice heat flux from the ocean to the ice [W/m²] |
| `Ocean_Ice_Boundary%sea_level` | Sea-surface height / sea level [m] |
| `Ocean_Ice_Boundary%fields` | Coupler boundary-condition type holding per-tracer ocean surface fields passed to the ice model |


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


## Required Variables in Component Datatypes

### Atmosphere
type (atmos_boundary_data_type) :: Atm
real, dimension(:) :: Atm%lon_bnd, Atm%lat_bnd
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
- `Atm%t_bot`, `Atm%q_bot`, `Atm%z_bot`, `Atm%p_bot`, `Atm%u_bot`, `Atm%v_bot`: State at the lowest atmosphere level.
- `Atm%p_surf`, `Atm%slp`, `Atm%gust`: Surface pressure, sea-level pressure, and gustiness factor.
- `Atm%flux_sw`, `Atm%flux_lw`, `Atm%lprec`, `Atm%fprec`, `Atm%coszen`: Surface radiative and precipitation inputs.
- `Atm%axes`: Axis identifiers returned by `diag_axis_init` for the atmospheric X, Y, `Z_full`, and `Z_half` axes.

The following fields support the simultaneous implicit time-stepping between the atmosphere and surface models:
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
```
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
```
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

### Fields on the ice **top** grid (atmosphere–ice interface):

```
real, dimension(:,:,:) :: Ice%t_surf      &
                          Ice%albedo      &
                          Ice%rough_mom   &
                          Ice%rough_heat  &
                          Ice%rough_moist &
                          Ice%u_surf      &
                          Ice%v_surf
```

### Fields on the ice **bottom** grid (ice–ocean interface), populated by flux_down_from_atmos` 
### and `flux_land_to_ice`, then passed to the ocean by `flux_ice_to_ocean`:

```
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
```
real, dimension(:,:) :: Ice%ustar_berg  &  ! iceberg friction velocity [m/s]
                        Ice%area_berg   &  ! iceberg area fraction
                        Ice%mass_berg      ! iceberg mass [kg/m2]
```

### Ocean
```
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
