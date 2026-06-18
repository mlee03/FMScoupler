# Flux Exchange

Original authors for this documentation: Bruce Wyman, V. Balaji, Sergey Malyshev

## Overview
Six modules couple the atmosphere, ocean, land, and sea-ice components through flux exchange. atm_land_ice_flux_exchange exchanges fluxes between atmosphere, land, and ice. atmos_ocean_fluxes_calc computes non-deposition gas fluxes between atmosphere and ocean. atmos_ocean_dep_fluxes_calc computes deposition gas fluxes between atmosphere and ocean. ice_ocean_flux_exchange exchanges fluxes between ice and ocean. land_ice_flux_exchange exchanges fluxes between land and ice. flux_exchange is the top-level module that initializes all flux-exchange modules and also contains subroutines for stock computation.

## Design Principles

### Grids 
The flux exchange supports physically independent atmosphere, land, and sea-ice grids; ice and ocean must share the same physical grid though their MPI domain decompositions may differ. The masked region of the land grid and the ice-ocean grid must tile each other such that every atmosphere grid cell is covered by either land or ice-ocean, but not both. The masked regions of the ice and ocean grids must be identical.

The atmosphere, land, and ice grids exchange information through the surface exchange grid xmap_sfc using conservative interpolation (REGRID). The land and ice grids exchange runoff data using the exchange grid xmap_runoff. Ice-bottom to ocean transfer does not require an exchange grid because those grids are physically identical; flux data are moved between PE layouts using mpp_redistribute when decompositions differ (REDIST), or copied directly when the decomposition is the same (DIRECT). Atmospheric fluxes reach the ocean through the ice model: first atmosphere to ice via flux_down_from_atmos, then ice to ocean via flux_ice_to_ocean.

### Time propagation
Sensible heat flux and surface evaporation can depend implicitly on surface temperature.  Therefore, land and sea-ice temperature updates must run on the atmospheric timestep. Surface fluxes for all other tracers and for momentum are treated as explicit functions of the surface state. The module is designed to support simultaneous implicit time integration on both sides of the surface interface, which requires that the diffusion part of the land and ice models also run on the atmospheric timestep. update_land_model_fast and update_ice_model_fast must update the surface temperature each atmospheric timestep for the implicit diffusion scheme to be correct.

### Additional fields
Additional tracer and gas-exchange fluxes are configured through field_table.  
Additional named boundary-condition fields are configured in the coupler boundary types. 
Any field exchanged between components can be replaced by a constant or file-based value using FMS data_override and data_table.
Each component model must expose a public data type containing the boundary fields needed by the coupler.

### Grid Layout
The three component grids tile the sphere. |xxx| marks a masked (inactive) grid point:

    ATMOSPHERE  |----|----|----|----|----|----|----|----|
          LAND  |---|---|---|---|xxx|xxx|xxx|xxx|xxx|xxx|
           ICE  |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
         OCEAN  |xxx|xxx|xxx|xxx|---|---|---|---|---|---|


## Transfer Types
There are three modes of flux exchange within FMSCoupler. REGRID (boundary_type%xtype=1) is used when the grids are physically distinct; it maps data through the exchange grid using conservative interpolation. REDIST (boundary_type%xtype=2) applies when component grids share the same physical grid but have different MPI decompositions, as is the case for ice and ocean when slow_ice_with_ocean=.true.; data are moved between PE layouts using mpp_redistribute. DIRECT (boundary_type%xtype=3) applies when the physical grid and MPI decomposition are identical, as for ice and ocean when slow_ice_with_ocean is not used; data are copied directly.

For REGRID, data on one component grid is first mapped onto the exchange grid, computations are carried out on the exchange grid, and the result is then mapped to the receiving component grid. Computed fields and fluxes can be overwritten by calls to data_override, but the override is applied only if the tracer is specified in the tracer_table.

## Namelist Parameters (flux_exchange_nml)

### namelist parameter z_ref_heat
z_ref_heat (namelist parameter in flux_exchange_nml) is a real parameter with a default of 2.0 that sets the reference height in meters for temperature and relative-humidity diagnostics (t_ref, rh_ref, del_h, del_q). 

### namelist parameter z_ref_mom
z_ref_mom (namelist parameter in flux_exchange_nml) is a real parameter with a default of 10.0 that sets the reference height in meters for momentum diagnostics (u_ref, v_ref, del_m).

### namelist parameter do_area_weighted_flux
do_area_weighted_flux (namelist parameter in flux_exchange_nml) is a logical parameter that defaults to .false.; when .true., fluxes passed to the ocean are multiplied by the ice area fraction before redistribution so the ocean receives the grid-cell-mean flux rather than the per-ice-area flux. 

### namelist parameter debug_stocks
debug_stocks (namelist parameter in flux_exchange_nml) is a logical parameter that defaults to .false. and enables additional stock-conservation debug output when .true.. divert_stocks_report is a logical parameter that defaults to .false. and redirects stock reporting to stocks.out rather than the standard log when .true.. 

### namelist parameter do_runoff
do_runoff (namelist parameter in flux_exchange_nml) is a logical parameter that defaults to .true. and enables interpolation of land runoff to the ocean. 

### namelist parameter do_forecast
do_forecast (namelist parameter in flux_exchange_nml) is a logical parameter that defaults to .false. and enables forecast-mode behavior in the flux coupler when .true..

### namelist parameter nblocks
nblocks (namelist parameter in flux_exchange_nml) is an integer parameter with a default of 1 that sets the number of blocks used to divide n_xgrid_sfc for OpenMP parallelism; in practice it is often set to match coupler_nml%atmos_nthreads, and if left at 1 when threading is active the model will emit a warning and reset it automatically. 

### namelist parameter partition_fprec_from_lprec
partition_fprec_from_lprec is a logical parameter that defaults to .false.; for atmosphere-override experiments where liquid and frozen precipitation are combined, it converts liquid precipitation to snow when t_ref < tfreeze. 

### namelist parameter scale_precip_2d
scale_precip_2d is a logical parameter that defaults to .false. and rescales lprec by a 2-D field read from data_table.

## Data Override Capabilities

### Context
Warning, the original authors strongly advise against using data override capabilities for "non-experts".
Any field in the lists below can be replaced at runtime with a constant or file-based value by adding a matching entry to data_table.
A data override is applied only when a matching entry exists; otherwise the model-computed value is used unchanged. 

### override in sfc_boundary_layer
In sfc_boundary_layer, these fields from the atmosphere boundary to the exchange grid can be overwritten:
t_bot, z_bot, p_bot, u_bot, v_bot, p_surf, slp, and gust, and fields in the coupler bc type.

These fields from the ice boundary to the exchange grid can be ovewritten:
t_surf, rough_mom, rough_heat, rough_moist, albedo, albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif, 
u_surf, v_surf.

These fields from the land boundary to the exchange grid can be overwritten:
t_surf, t_ca, rough_mom, rough_heat, albedo, tracers, albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif

These fields from the exchange grid to Land_ice_atmos_boundary can be overwritten:
t, albedo, albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif, land_frac, dt_t, dt_tr, u_flux, 
v_flux, dtaudu, dtaudv, u_star, b_star, rough_mom

### override in flux_down_from_atmos
In flux_down_from_atmos, these fields from atmosphere boundary to the exchange grid can be overwritten:
flux_sw, flux_sw_dir, flux_sw_dif, flux_sw_down_vis_dir, flux_sw_down_vis_dif, flux_sw_down_total_dir, 
flux_sw_down_total_dif, flux_sw_vis, flux_sw_vis_dir, flux_sw_vis_dif, flux_lw, lprec, frac_precip, 
fprec, coszen, dtmass, delta_t, dflux_t, delta_tr, dflux_tr

These fields from the exchange grid to the land boundary can be overwritten:
drag_q, lwdn_flux, cd_m, cd_t, bstar, ustar, wind, z_bot, t_flux, lw_flux, sw_flux, sw_flux_down_vis_dir, 
sw_flux_down_total_dir, sw_flux_down_vis_dif, sw_flux_down_total_dif, lprec, fprec, dhdt, drdt, p_surf, tr_flux, dfdtr

These fields from the exchange grid to the ice boundary can be overwritten:
u_flux, v_flux, t_flux, q_flux, lw_flux, sw_flux_nir_dir, sw_flux_vis_dir, sw_flux_nir_dif, sw_flux_vis_dif,
sw_down_vis_dir, sw_down_vis_dif, sw_down_nir_dir, sw_down_nir_dif, lprec, fprec, dhdt, dedt, drdt, coszen, p

### override in flux_up_to_atmos
In flux_up_to_atmos, this field from the ice boundary to atmosphere boundary can be overwritten:  t_surf

These fields from the land boundary to atmosphere boundary can be overwritten:
t_ca, t_surf, tr

### override in flux_land_to_ice
In flux_land_to_ice, these fields from the land boundary to ice boundary can be overwritten:  
runoff, calving, runoff_hflx, calving_hflx.  Note, do_runoff namelist flag must be .true. for this exchange to occur.

### override in flux_ice_to_ocean
In flux_ice_to_ocean, these fields from the ice boundary to the ocean boundary can be overwritten:
u_flux, v_flux, t_flux, q_flux, salt_flux, lw_flux, sw_flux_nir_dir, sw_flux_nir_dif, sw_flux_vis_dir, 
sw_flux_vis_dif, lprec, fprec, runoff, calving, runoff_hflx, calving_hflx, p, mi, ustar_berg, area_berg, mass_berg

### override in flux_ocean_to_ice
In flux_ocean_to_ice, these fields from the ocean boundary to ice boundary can be overwritten:
u, v, t, s, frazil, sea_level

## Diagnostic Fields

All fields below are registered in atm_land_ice_flux_exchange.F90 inside diag_field_init. 

static fields:
land_mask, height2m, height10m, sftlf

atmosphere fields:
ice_mask, wind, drag_moist, drag_heat, drag_mom, rough_moist, rough_heat, rough_mom, u_star, b_star, q_star, thv_atm, thv_surf, tau_x, tau_y, t_ocean, t_surf, t_ca, z_atm, p_atm, slp, gust, shflx, lwflx, t_atm, u_atm, v_atm, t_ref, rh_ref, rh_ref_cmip, u_ref, v_ref, wind_ref, del_h, del_m, del_q, q_ref, rough_scale, evap, co2_bot

atmosphere tracer fields:
*_tot_con_atm, *_tot_con_ref, *_atm, *_surf, *_flux, *_ref, *_mol_flux, *_atm_dvmr, *_surf_dvmr, *_mol_flux_atm0

CMIP fields that are registerd with  register_cmip_diag_field_2d or fms_diag_register_diag_field with use_AM3_physics:
tas, uas, vas, sfcWind, huss, hurs, rhs, ts, psl, tauu, tauv, hfss, hfls, evspsbl, tslsi, tos, sic

These fields (registered with register_global_diag_field) to produce globally averaged scalar time-series output when not using use_AM3_physics:
evspsbl, ts, tas, tasl, hfss, hfls, rls

These fields are registered on the land axes using register_tiled_diag_field or fms_diag_register_diag_field with _USE_LEGACY_LAND_:
t_ref, q_ref, rh_ref, u_ref, v_ref, evap, shflx, tasLut, hussLut, *_tot_con_atm, *_tot_con_ref, *_flux, *_mol_flux, *_ref


## Required Variables in Component Data Types

### Atmosphere (atmos_data_type)
These atmos_data_type fields (for example Atm%lon_bnd) must be defined:

lon_bnd, lat_bnd   
t_bot, q_bot, z_bot, p_bot, u_bot, v_bot, p_surf, slp, gust, flux_sw, flux_lw, 
lprec, fprec, coszen
axes  

As additional notes, lon_bnd and lat_bnd are grid-box corner coordinates in radians and in monotonic order.
Fields with "_bot" are states at the lowest atmospheric model level and are the primary inputs to sfc_boundary_layer.
p_surf, slp, gust are the surface pressure, sea-level pressure, and gustiness factor, respectively.
flux_sw, flux_lw, lprec, fprec, coszen are surface radiative and precipitation forcing passed to land and ice by flux_down_from_atmos.
axes are axis IDs returned by diag_axis_init and are required for diagnostic registration with FMS.

The following fields support implicit time-stepping between the atmosphere and the surface models:
dtmass, delta_t, delta_q, dflux_t, dflux_q     


### Land (land_data_type)
These land_data_type (for example, Land%albedo) must be defined:

lon_bnd, lat_bnd
mask, glacier
tile_size, t_surf, t_ca, q_ca, albedo, rough_mom, rough_heat, stomatal, snow, water, max_water

As additional notes, lon_bnd and lat_bnd are grid-box corner coordinates in radians and must be monotonic. 
Mask is the land-sea mask and is .true. over land points. 
Glacier is the glacier mask and is .true. over glacier points. 
Tile_siz, valued between 0 and 1, is the fractional area of each land tile within the atmospheric grid cell. 
t_surf, albedo, rough_mom, and rough_heat are the surface state fields used for turbulent flux and radiation calculations. 
t_ca and q_ca are the canopy air temperature and specific humidity, which are returned to the atmosphere by flux_up_to_atmos. 
stomatal, snow, water, and max_water are additional surface properties used in flux parameterizations.

### Ice (ice_data_type)

These ice_data_type fields (for example, Ice%lon_bnd) must be defined:

lon_bnd, lat_bnd, lon_bnd_uv, lat_bnd_uv
mask, mask_uv, ice_mask
part_size, part_size_uv

All boundary arrays are in radians and must be monotonic. 
mask and mask_uv are the ocean-land masks for tracer and momentum points, respectively. 
ice_mask is an optional explicit sea-ice mask. 
part_size and part_size_uv are the fractional area of each ice thickness category.

These fields on ice top (atmosphere–ice interface), and are provided to the atmosphere each fast timestep:
t_surf, albedo, rough_mom, rough_heat, rough_moist, u_srf, and v_surf 

These fields on the ice bottom (ice–ocean interface) are populated by flux_down_from_atmos 
and flux_land_to_ice, then passed to the ocean by flux_ice_to_ocean:
flux_u, flux_v, flux_t, flux_q, flux_salt, flux_lw,
flux_sw_vis_dir, flux_sw_vis_dif, flux_sw_nir_dir, flux_sw_nir_dif,
lprec, fprec, runoff, calving, runoff_hflx, calving_hflx, p_surf

These optional iceberg fields, allocated only when the iceberg module is active:
ustar_berg, area_berg, mass_berg

### Ocean (ocean_public_type)

These ocean_public_type fields must be defined:
t_surf, s_surf, u_surf, v_surf, frazil, sea_lev, Data%mask, Data%mask_uv, Ocean%mask, Ocean%mask_uv
