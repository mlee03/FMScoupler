# land_ice_atmos_boundary_type
land_ice_atmos_boundary_type contains surface quantities going from land and ice to the atmosphere.
All quantities are on the exchange grid. 

## Land_ice_atmos_boundary%t
Land_ice_atmos_boundary%t, a real 2D array, is Area-weighted surface temperature seen by the atmosphere for radiation calculations [K]; weighted over land and ice fractions.
## Land_ice_atmos_boundary%t_ocean
Land_ice_atmos_boundary%t_ocean, a real 2D array, is Ocean surface temperature for radiation calculations; sourced from Ice%t_surf through the exchange grid [K].
## Land_ice_atmos_boundary%albedo
Land_ice_atmos_boundary%albedo, a real 2D array, is Broadband surface albedo [dimensionless].
## Land_ice_atmos_boundary%albedo_vis_dir
Land_ice_atmos_boundary%albedo_vis_dir, a real 2D array, is Direct-beam visible-band surface albedo [dimensionless].
## Land_ice_atmos_boundary%albedo_nir_dir
Land_ice_atmos_boundary%albedo_nir_dir, a real 2D array, is Direct-beam near-infrared surface albedo [dimensionless].
## Land_ice_atmos_boundary%albedo_vis_dif
Land_ice_atmos_boundary%albedo_vis_dif, a real 2D array, is Diffuse visible-band surface albedo [dimensionless].
## Land_ice_atmos_boundary%albedo_nir_dif
Land_ice_atmos_boundary%albedo_nir_dif, a real 2D array, is Diffuse near-infrared surface albedo [dimensionless].
## Land_ice_atmos_boundary%land_frac
Land_ice_atmos_boundary%land_frac, a real 2D array, is Fraction of the atmospheric grid cell covered by land [dimensionless].
## Land_ice_atmos_boundary%frac_open_sea
Land_ice_atmos_boundary%frac_open_sea, a real 2D array, is Non-sea-ice fraction of the grid cell [dimensionless]; complement of the sea-ice concentration.
## Land_ice_atmos_boundary%rough_mom
Land_ice_atmos_boundary%rough_mom, a real 2D array, is Area-weighted surface roughness length for momentum [m].
## Land_ice_atmos_boundary%rough_heat
Land_ice_atmos_boundary%rough_heat, a real 2D array, is Area-weighted surface roughness length for heat [m].

## Land_ice_atmos_boundary%u_ref
Land_ice_atmos_boundary%u_ref, a real 2D array, is Zonal wind at the momentum reference height (z_ref_mom) [m/s].
This field is interpolated from the surface-layer profile to the reference heights z_ref_heat and z_ref_mom configured in flux_exchange_nml.

## Land_ice_atmos_boundary%v_ref
Land_ice_atmos_boundary%v_ref, a real 2D array, is Meridional wind at the momentum reference height (z_ref_mom) [m/s].
This field is interpolated from the surface-layer profile to the reference heights z_ref_heat and z_ref_mom configured in flux_exchange_nml.

## Land_ice_atmos_boundary%t_ref
Land_ice_atmos_boundary%t_ref, a real 2D array, is Air temperature at the heat reference height (z_ref_heat) [K].
This field is interpolated from the surface-layer profile to the reference heights z_ref_heat and z_ref_mom configured in flux_exchange_nml.

## Land_ice_atmos_boundary%q_ref
Land_ice_atmos_boundary%q_ref, a real 2D array, is Specific humidity at the heat reference height (z_ref_heat) [kg/kg].
This field is interpolated from the surface-layer profile to the reference heights z_ref_heat and z_ref_mom configured in flux_exchange_nml.

## Land_ice_atmos_boundary%wind
Land_ice_atmos_boundary%wind, a real 2D array, is Absolute wind speed at the lowest atmospheric model level including gust corrections [m/s].
This field is interpolated from the surface-layer profile to the reference heights z_ref_heat and z_ref_mom configured in flux_exchange_nml.

## Land_ice_atmos_boundary%thv_atm
Land_ice_atmos_boundary%thv_atm, a real 2D array, is Virtual potential temperature at the lowest atmospheric model level [K].
This field is interpolated from the surface-layer profile to the reference heights z_ref_heat and z_ref_mom configured in flux_exchange_nml.

## Land_ice_atmos_boundary%thv_surf
Land_ice_atmos_boundary%thv_surf, a real 2D array, is Virtual potential temperature at the surface [K].
This field is interpolated from the surface-layer profile to the reference heights z_ref_heat and z_ref_mom configured in flux_exchange_nml.

## Land_ice_atmos_boundary%dt_t
Land_ice_atmos_boundary%dt_t, a real 2D array, is Temperature tendency correction at the lowest atmospheric level from the implicit surface flux scheme [K/s].  This field is an output of the tridiagonal back-substitution step. They correct the atmospheric temperature and tracer tendencies for the implicit surface coupling.
## Land_ice_atmos_boundary%dt_tr
Land_ice_atmos_boundary%dt_tr, a real 3D array, is Tracer mixing-ratio tendency correction at the lowest level; third dimension indexes tracers [tracer units/s].  This field is an output of the tridiagonal back-substitution step. They correct the atmospheric temperature and tracer tendencies for the implicit surface coupling.

## Land_ice_atmos_boundary%u_flux
Land_ice_atmos_boundary%u_flux, a real 2D array, is Zonal wind stress on the atmosphere [Pa].
## Land_ice_atmos_boundary%v_flux
Land_ice_atmos_boundary%v_flux, a real 2D array, is Meridional wind stress on the atmosphere [Pa].
## Land_ice_atmos_boundary%dtaudu
Land_ice_atmos_boundary%dtaudu, a real 2D array, is d(zonal wind stress)/d(u) — implicit coupling coefficient for zonal momentum [Pa·s/m].
## Land_ice_atmos_boundary%dtaudv
Land_ice_atmos_boundary%dtaudv, a real 2D array, is d(meridional wind stress)/d(v) — implicit coupling coefficient for meridional momentum [Pa·s/m].
## Land_ice_atmos_boundary%u_star
Land_ice_atmos_boundary%u_star, a real 2D array, is Friction velocity (surface turbulent velocity scale) [m/s].
## Land_ice_atmos_boundary%b_star
Land_ice_atmos_boundary%b_star, a real 2D array, is Buoyancy scale used in Monin-Obukhov similarity theory [m/s**2].
## Land_ice_atmos_boundary%q_star
Land_ice_atmos_boundary%q_star, a real 2D array, is Moisture scale used in Monin-Obukhov similarity theory [kg/kg].
## Land_ice_atmos_boundary%shflx
Land_ice_atmos_boundary%shflx, a real 2D array, is Sensible heat flux at the surface W/m**2; not compiled when use_AM3_physics is defined.
## Land_ice_atmos_boundary%lhflx
Land_ice_atmos_boundary%lhflx, a real 2D array, is Latent heat flux at the surface W/m**2; not compiled when use_AM3_physics is defined.

## Land_ice_atmos_boundary%data
Land_ice_atmos_boundary%data, a real 3D array, is Collective array providing named access to the scalar fields above; used internally for data-override and exchange-grid operations.
## Land_ice_atmos_boundary%gex_lnd2atm
Land_ice_atmos_boundary%gex_lnd2atm, a real 3D array, is Generic exchange fields returned from the land model to the atmosphere (e.g., surface emission fluxes); third dimension indexes the exchange field list.
## Land_ice_atmos_boundary%xtype
Land_ice_atmos_boundary%xtype, integer, is Transfer mode for the exchange-grid-to-atmosphere remap: REGRID (1), REDIST (2), or DIRECT (3).
