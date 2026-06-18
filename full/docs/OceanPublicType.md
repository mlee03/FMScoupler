# ocean_public_type
ocean_public_type is the publicly visible face of the MOM6 ocean model.
It contains only the surface fields and domain metadata for the coupler.

## Ocean%Domain
Ocean%Domain, type(domain2d), is FMS domain decomposition for the ocean surface fields; defines the MPI tile layout used for coupler remapping.
## Ocean%is_ocean_pe
Ocean%is_ocean_pe, logical, is .true. on processors that run the ocean model; used throughout the coupler to gate ocean-only code paths.
## Ocean%pelist
Ocean%pelist, integer(:), is List of MPI PE numbers assigned to the ocean model.
## Ocean%maskmap
Ocean%maskmap, logical(:,:), is Pointer to a mask indicating which logical processors are active for ocean computation; logical processors covering all-land points may not be mapped to physical PEs. Need not be set if all processors are used.
## Ocean%instance_name
Ocean%instance_name, character(32), is Optional name identifying this ocean model instance; used in ensemble runs to disambiguate log messages.

## Ocean%t_surf
Ocean%t_surf, is Sea-surface temperature (SST) on tracer (T) cells [K].
This field is set by the ocean model after each update and read by the
coupler to construct the ocean_ice_boundary_type passed to the sea-ice
model via flux_ocean_to_ice.

## Ocean%s_surf
Ocean%s_surf, is Sea-surface salinity (SSS) on T cells [ppt].
This field is set by the ocean model after each update and read by the
coupler to construct the ocean_ice_boundary_type passed to the sea-ice
model via flux_ocean_to_ice.

## Ocean%u_surf
Ocean%u_surf, is Surface current i-velocity at the locations indicated by stagger [m/s].
This field is set by the ocean model after each update and read by the
coupler to construct the ocean_ice_boundary_type passed to the sea-ice
model via flux_ocean_to_ice.

## Ocean%v_surf
Ocean%v_surf, is Surface current j-velocity at the locations indicated by stagger [m/s].
This field is set by the ocean model after each update and read by the
coupler to construct the ocean_ice_boundary_type passed to the sea-ice
model via flux_ocean_to_ice.

## Ocean%sea_lev
Ocean%sea_lev, is Sea level corrected for surface pressure: dzt(1) + η + p_atm/(ρ₀g) [m]; passed to the ice model as sea_level.
This field is set by the ocean model after each update and read by the
coupler to construct the ocean_ice_boundary_type passed to the sea-ice
model via flux_ocean_to_ice.

## Ocean%frazil
Ocean%frazil, is Accumulated heating from frazil ice formation in the ocean since the last coupling step [J/m**2]; delivered to the ice model so it can account for ocean-side freezing.
This field is set by the ocean model after each update and read by the
coupler to construct the ocean_ice_boundary_type passed to the sea-ice
model via flux_ocean_to_ice.

## Ocean%melt_potential
Ocean%melt_potential, is Instantaneous heat available to melt sea ice from below [J/m**2]; computed when the ocean boundary layer depth exceeds HFrz.
This field is set by the ocean model after each update and read by the
coupler to construct the ocean_ice_boundary_type passed to the sea-ice
model via flux_ocean_to_ice.

## Ocean%OBLD
Ocean%OBLD, is Ocean boundary layer depth [m]; used to determine the depth over which melt potential is computed.
This field is set by the ocean model after each update and read by the
coupler to construct the ocean_ice_boundary_type passed to the sea-ice
model via flux_ocean_to_ice.

## Ocean%area
Ocean%area, is Grid-cell area of each ocean surface cell [m**2]; used for conservative flux remapping.
This field is set by the ocean model after each update and read by the
coupler to construct the ocean_ice_boundary_type passed to the sea-ice
model via flux_ocean_to_ice.

## Ocean%calving
Ocean%calving, is Mass per unit area of ice-shelf flux to be converted to icebergs [kg/m**2]; passed to the iceberg module.
This field is set by the ocean model after each update and read by the
coupler to construct the ocean_ice_boundary_type passed to the sea-ice
model via flux_ocean_to_ice.

## Ocean%calving_hflx
Ocean%calving_hflx, is Heat flux associated with calving W/m**2.
This field is set by the ocean model after each update and read by the
coupler to construct the ocean_ice_boundary_type passed to the sea-ice
model via flux_ocean_to_ice.

## Ocean%stagger
Ocean%stagger, integer, is Arakawa staggering of the surface velocity components (u_surf, v_surf) relative to tracer points. Valid values: AGRID, BGRID_NE, CGRID_NE, BGRID_SW, CGRID_SW. Set to -999 before initialization so a global max can propagate the value to non-ocean PEs.

## Ocean%fields
Ocean%fields, type(coupler_2d_bc_type), is Named arrays of tracer-related ocean surface fields (e.g. pCO₂, O₂ saturation) used in atmosphere-ocean gas flux calculations; populated by atmos_ocean_fluxes_calc.
## Ocean%avg_kount
Ocean%avg_kount, integer, is Counter tracking the number of contributions accumulated in the running averages stored in this type; used externally by FMSCoupler to manage time averaging.
## Ocean%axes(2)
Ocean%axes(2), integer(2), is Diag-manager axis IDs available for I/O using this surface data.
