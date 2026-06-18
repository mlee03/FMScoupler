# atmos_land_boundary_type
atmos_land_boundary_type carries all data passed from the coupler to the land model (LM4).
All fields are pointers with dimension (grid_idex, tile_number) 

## Atmos_land_boundary%t_flux
Atmos_land_boundary%t_flux, a real 2D array, is Sensible heat flux into the land surface W/m**2.
## Atmos_land_boundary%lw_flux
Atmos_land_boundary%lw_flux, a real 2D array, is Net longwave radiation flux at the land surface W/m**2.
## Atmos_land_boundary%lwdn_flux
Atmos_land_boundary%lwdn_flux, a real 2D array, is Downward longwave radiation flux at the land surface W/m**2.
## Atmos_land_boundary%sw_flux
Atmos_land_boundary%sw_flux, a real 2D array, is Net shortwave radiation flux at the land surface W/m**2.
## Atmos_land_boundary%swdn_flux
Atmos_land_boundary%swdn_flux, a real 2D array, is Downward shortwave radiation flux at the land surface W/m**2.
## Atmos_land_boundary%sw_flux_down_vis_dir
Atmos_land_boundary%sw_flux_down_vis_dir, a real 2D array, is Downward direct-beam visible shortwave flux W/m**2.
## Atmos_land_boundary%sw_flux_down_total_dir
Atmos_land_boundary%sw_flux_down_total_dir, a real 2D array, is Downward direct-beam total (broadband) shortwave flux W/m**2.
## Atmos_land_boundary%sw_flux_down_vis_dif
Atmos_land_boundary%sw_flux_down_vis_dif, a real 2D array, is Downward diffuse visible shortwave flux W/m**2.
## Atmos_land_boundary%sw_flux_down_total_dif
Atmos_land_boundary%sw_flux_down_total_dif, a real 2D array, is Downward diffuse total (broadband) shortwave flux W/m**2.
## Atmos_land_boundary%lprec
Atmos_land_boundary%lprec, a real 2D array, is Liquid precipitation rate [kg/m**2/s].
## Atmos_land_boundary%fprec
Atmos_land_boundary%fprec, a real 2D array, is Frozen precipitation rate [kg/m**2/s].
## Atmos_land_boundary%tprec
Atmos_land_boundary%tprec, a real 2D array, is Temperature of precipitation [K].

## Atmos_land_boundary%dhdt
Atmos_land_boundary%dhdt, a real 2D array, is d(sensible heat flux)/d(T_surf) — derivative of sensible heat flux with respect to surface temperature [W/m**2/K].  It is a derivative quantity needed to close the implicit tridiagonal surface diffusion scheme between the atmosphere and land.

## Atmos_land_boundary%dhdq
Atmos_land_boundary%dhdq, a real 2D array, is d(sensible heat flux)/d(q_surf) — derivative of sensible heat flux with respect to surface specific humidity [W/m**2/(kg/kg)].  It is a derivative quantity needed to close the implicit tridiagonal surface diffusion scheme between the atmosphere and land.

## Atmos_land_boundary%drdt
Atmos_land_boundary%drdt, a real 2D array, is d(longwave flux)/d(T_surf) — derivative of longwave flux with respect to surface radiative temperature [W/m**2/K].  It is a derivative quantity needed to close the implicit tridiagonal surface diffusion scheme between the atmosphere and land.

## Atmos_land_boundary%cd_m
Atmos_land_boundary%cd_m, a real 2D array, is Drag coefficient for momentum [dimensionless].
## Atmos_land_boundary%cd_t
Atmos_land_boundary%cd_t, a real 2D array, is Drag coefficient for tracers (heat and moisture) [dimensionless].
## Atmos_land_boundary%ustar
Atmos_land_boundary%ustar, a real 2D array, is Turbulent wind scale (friction velocity) [m/s].
## Atmos_land_boundary%bstar
Atmos_land_boundary%bstar, a real 2D array, is Turbulent buoyancy scale [m/s].
## Atmos_land_boundary%wind
Atmos_land_boundary%wind, a real 2D array, is Absolute wind speed at the bottom of the atmospheric layer [m/s].
## Atmos_land_boundary%z_bot
Atmos_land_boundary%z_bot, a real 2D array, is Height of the bottom atmospheric layer above the land surface [m].
## Atmos_land_boundary%drag_q
Atmos_land_boundary%drag_q, a real 2D array, is Product of the moisture drag coefficient and wind speed (cd_q × wind); used in land surface moisture flux calculations [m/s].
## Atmos_land_boundary%p_surf
Atmos_land_boundary%p_surf, a real 2D array, is Surface pressure [Pa].

## Tracer fluxes
dimension (grid_index, tile_number, tracer_index)

## Atmos_land_boundary%tr_flux
Atmos_land_boundary%tr_flux, a real 3D array, is Flux of each tracer into the land surface, including water vapor flux; dimensioned (grid_index, tile, tracer) [tracer units · kg air / (m**2·s)].
## Atmos_land_boundary%dfdtr
Atmos_land_boundary%dfdtr, a real 3D array, is d(tracer flux)/d(tracer_surf) — derivative of the tracer flux with respect to the surface tracer value, including evaporation over surface specific humidity; dimensioned (grid_index, tile, tracer).

## Atmos_land_boundary%xtype
Atmos_land_boundary%xtype, integer, is Transfer mode for the atmosphere-to-land exchange: REGRID (1), REDIST (2), or DIRECT (3).
