# land_data_type
land_data_type carries land surface states passed to the couper.
The land model uses an unstructured tile representation where
multiple land-cover types (soil, vegetation,
lakes, etc.) can co-exist within a single atmospheric grid cell.


## Land%tile_size
Land%tile_size, is Fractional coverage of the atmospheric grid cell by this tile [dimensionless, 0–1]; used to area-weight tile quantities back onto the atmosphere grid.  The variable has dimension of (grid_index, tile_number).
## Land%t_surf
Land%t_surf, is Ground (radiative) surface temperature [K]; used in longwave radiation and sensible heat flux calculations.
The variable has dimension of (grid_index, tile_number).
## Land%t_ca
Land%t_ca, is Canopy air temperature — near-surface air temperature within the plant canopy layer [K]; differs from t_surf over vegetated tiles.
The variable has dimension of (grid_index, tile_number).
## Land%albedo
Land%albedo, is Broadband surface albedo [dimensionless]; legacy field, per-band albedos below are preferred.
The variable has dimension of (grid_index, tile_number).
## Land%albedo_vis_dir
Land%albedo_vis_dir, is Surface albedo for direct-beam visible radiation (0.2–0.7 µm) [dimensionless].
The variable has dimension of(grid_index, tile_number).
## Land%albedo_nir_dir
Land%albedo_nir_dir, is Surface albedo for direct-beam near-infrared radiation [dimensionless].
The variable has dimension of (grid_index, tile_number).
## Land%albedo_vis_dif
Land%albedo_vis_dif, is Surface albedo for diffuse visible radiation [dimensionless].
The variable has dimension of (grid_index, tile_number).
## Land%albedo_nir_dif
Land%albedo_nir_dif, is Surface albedo for diffuse near-infrared radiation [dimensionless].
The variable has dimension of (grid_index, tile_number).
## Land%rough_mom
Land%rough_mom, is Surface roughness length for momentum [m]; used in Monin-Obukhov flux calculations.
The variable has dimension of (grid_index, tile_number).
## Land%rough_heat
Land%rough_heat, is Surface roughness length for heat and tracers [m].
The variable has dimension of (grid_index, tile_number).
## Land%rough_scale
Land%rough_scale, is Topographic form-drag scaling factor for momentum [dimensionless]; accounts for sub-grid orographic drag.
The variable has dimension of (grid_index, tile_number).

## Land%tr
Land%tr, is Surface tracer mixing ratios on each tile, including canopy air specific humidity as the first tracer; additional tracers (e.g. CO₂) follow the tracer table order.  The variable has dimension of (grid_index, tile_number, tracer_index)

## Land%discharge
Land%discharge, is Liquid water discharge (river runoff) from land to ocean [kg/m**2/s].
This field carries freshwater and heat leaving the land surface and routed to the ocean/ice via flux_land_to_ice.
The field has dimension of (lon, lat).
## Land%discharge_heat
Land%discharge_heat, is Sensible heat carried by liquid discharge, using 0 °C as datum W/m**2.
This field carries freshwater and heat leaving the land surface and routed to the ocean/ice via flux_land_to_ice.
The field has dimension of (lon, lat).
## Land%discharge_snow
Land%discharge_snow, is Solid water (snow/ice) discharge from land to ocean [kg/m**2/s].
This field carries freshwater and heat leaving the land surface and routed to the ocean/ice via flux_land_to_ice.
The field has dimension of (lon, lat).
## Land%discharge_snow_heat
Land%discharge_snow_heat, is Sensible heat carried by solid discharge, using 0 °C as datum W/m**2.
This field carries freshwater and heat leaving the land surface and routed to the ocean/ice via flux_land_to_ice.
The field has dimension of (lon, lat).

## Land%mask
Land%mask, a logical 2D array, is .true. where the grid cell contains land; used to gate land-only computations.
## Land%axes(1)
Land%axes(1), integer, is Diag-manager axis ID for the unstructured land grid; used when registering tiled land diagnostics.
## Land%domain
Land%domain, type(domain2D), is FMS structured-grid domain for the land model; used for halo exchanges and exchange-grid setup.
## Land%ug_domain
Land%ug_domain, type(domainUG), is FMS unstructured-grid domain for the land model; carries the tile-based decomposition used by LM4.
## Land%pelist
Land%pelist, integer 1D array, is List of MPI PE numbers on which the land model is running.
## Land%pe
Land%pe, logical, is .true. on PEs that are part of the land pelist; used to gate land stock calculations.

