# atmos_ice_boundary_type
atmos_ice_boundary_type holds data passed between atmosphere and the sea-ice model (SIS2).

---

## Wind stress

## Atmos_ice_boundary_type%u_flux
Atmos_ice_boundary_type%u_flux, a real 3D array, is True-eastward wind stress from the atmosphere to the ocean or ice in each thickness category, on an A-grid and not rotated to the model grid [Pa].
## Atmos_ice_boundary_type%v_flux
Atmos_ice_boundary_type%v_flux, a real 3D array, is True-northward wind stress from the atmosphere to the ocean or ice in each thickness category, on an A-grid and not rotated to the model grid [Pa].
## Atmos_ice_boundary_type%u_star
Atmos_ice_boundary_type%u_star, a real 3D array, is Atmospheric friction velocity on an A-grid [Pa].

## Atmos_ice_boundary_type%t_flux
Atmos_ice_boundary_type%t_flux, a real 3D array, is Net sensible heat flux from the ocean or ice surface into the atmosphere W/m**2.
## Atmos_ice_boundary_type%q_flux
Atmos_ice_boundary_type%q_flux, a real 3D array, is Moisture flux from the ice or ocean to the atmosphere due to evaporation or sublimation [kg/m²/s].

## Atmos_ice_boundary_type%dhdt
Atmos_ice_boundary_type%dhdt, a real 3D array, is d(upward sensible heat flux)/d(T_surf) — derivative of sensible heat flux with respect to surface temperature [W/m²/°C].  It's a linearization (derivative) term needed to close the implicit tridiagonal surface diffusion scheme between the atmosphere and ice.
## Atmos_ice_boundary_type%dedt
Atmos_ice_boundary_type%dedt, a real 3D array, is d(sublimation+evaporation rate)/d(T_surf) — derivative of the moisture flux with respect to surface temperature [kg/m²/s/°C].  It's a linearization (derivative) term needed to close the implicit tridiagonal surface diffusion scheme between the atmosphere and ice.
## Atmos_ice_boundary_type%drdt
Atmos_ice_boundary_type%drdt, a real 3D array, is d(net upward longwave flux)/d(T_surf) — derivative of the net upward longwave flux (i.e. -lw_flux) with respect to surface temperature [W/m²/°C].  It's a linearization (derivative) term needed to close the implicit tridiagonal surface diffusion scheme between the atmosphere and ice.

## Atmos_ice_boundary_type%lw_flux
Atmos_ice_boundary_type%lw_flux, a real 3D array, is Net downward longwave radiation flux from the atmosphere into the ice or ocean W/m**2.
## Atmos_ice_boundary_type%sw_flux_vis_dir
Atmos_ice_boundary_type%sw_flux_vis_dir, a real 3D array, is Net direct visible shortwave radiation flux into the ice or ocean W/m**2.
## Atmos_ice_boundary_type%sw_flux_vis_dif
Atmos_ice_boundary_type%sw_flux_vis_dif, a real 3D array, is Net diffuse visible shortwave radiation flux into the ice or ocean W/m**2.
## Atmos_ice_boundary_type%sw_flux_nir_dir
Atmos_ice_boundary_type%sw_flux_nir_dir, a real 3D array, is Net direct near-infrared shortwave radiation flux into the ice or ocean W/m**2.
## Atmos_ice_boundary_type%sw_flux_nir_dif
Atmos_ice_boundary_type%sw_flux_nir_dif, a real 3D array, is Net diffuse near-infrared shortwave radiation flux into the ice or ocean W/m**2.
## Atmos_ice_boundary_type%sw_down_vis_dir
Atmos_ice_boundary_type%sw_down_vis_dir, a real 3D array, is Downward direct visible shortwave radiation flux from the atmosphere W/m**2.
## Atmos_ice_boundary_type%sw_down_vis_dif
Atmos_ice_boundary_type%sw_down_vis_dif, a real 3D array, is Downward diffuse visible shortwave radiation flux from the atmosphere W/m**2.
## Atmos_ice_boundary_type%sw_down_nir_dir
Atmos_ice_boundary_type%sw_down_nir_dir, a real 3D array, is Downward direct near-infrared shortwave radiation flux from the atmosphere W/m**2.
## Atmos_ice_boundary_type%sw_down_nir_dif
Atmos_ice_boundary_type%sw_down_nir_dif, a real 3D array, is Downward diffuse near-infrared shortwave radiation flux from the atmosphere W/m**2.
## Atmos_ice_boundary_type%coszen
Atmos_ice_boundary_type%coszen, a real 3D array, is Cosine of the solar zenith angle averaged over the next radiation timestep (not the timestep used to compute the sw_flux fields) [dimensionless, ≤ 1].

## Atmos_ice_boundary_type%lprec
Atmos_ice_boundary_type%lprec, a real 3D array, is Liquid precipitation (rain) from the atmosphere onto the ice or ocean in each thickness category [kg/m²/s]; rain falling on snow is currently assumed to drain directly through the ice into the ocean.
## Atmos_ice_boundary_type%fprec
Atmos_ice_boundary_type%fprec, a real 3D array, is Frozen precipitation (snowfall, sleet, hail, graupel) from the atmosphere to the ice or ocean [kg/m²/s]; all forms of frozen precipitation are treated as snow in SIS2.

## Atmos_ice_boundary_type%p
Atmos_ice_boundary_type%p, a real 3D array, is Atmospheric surface pressure [Pa]; typically ~1×10⁵ Pa.

## Atmos_ice_boundary_type%xtype
Atmos_ice_boundary_type%xtype, integer, is Transfer mode for the atmosphere-to-ice exchange: REGRID (1), REDIST (2), or DIRECT (3).
## Atmos_ice_boundary_type%fluxes
Atmos_ice_boundary_type%fluxes, type(coupler_3d_bc_type), is Array of additional per-tracer gas and deposition fluxes from the atmosphere to the ice.
