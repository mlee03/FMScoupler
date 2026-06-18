# `atmos_data_type`
is the main derived type holding fields and states of the atmosphere.  Below are all the fields in atmos_data_type:


## Atm%domain
Atm%domain, type(domain2d), is FMS domain decomposition object for the atmosphere; defines the MPI tile layout and halo widths.
## Atm%axes
Atm%axes, integer(4), is Diag-manager axis indices for x, y, pfull, and phalf; used when registering and sending diagnostic fields.
## Atm%lon_bnd
Atm%lon_bnd, a real 2D array, is Longitude of grid-box corners on the local compute domain [radians].
## Atm%lat_bnd
Atm%lat_bnd, a real 2D array, is Latitude of grid-box corners on the local compute domain [radians].
## Atm%lon
Atm%lon, a real 2D array, is Longitude of grid-box centres on the local compute domain [radians].
## Atm%lat
Atm%lat, a real 2D array, is Latitude of grid-box centres on the local compute domain [radians].
## Atm%grid
Atm%grid, type(grid_box_type), is Grid geometry needed for second-order conservative remapping on the cubic-sphere exchange grid (see below).
## Atm%maskmap
Atm%maskmap, logical(:,:), is Pointer to a mask indicating which logical processors are active for ocean code; processors covering all-land points may not be assigned to physical PEs. Dummy field — must be present for compilation but need not be set.
## Atm%t_bot
Atm%t_bot, a real 2D array, is Temperature at the lowest model level [K].
## Atm%tr_bot
Atm%tr_bot, a real 3D array, is Tracer mixing ratios at the lowest model level; third dimension indexes the tracer table (specific humidity sphum is always present).
## Atm%z_bot
Atm%z_bot, a real 2D array, is Height of the lowest model level above the surface [m].
## Atm%p_bot
Atm%p_bot, a real 2D array, is Pressure at the lowest model level [Pa].
## Atm%u_bot
Atm%u_bot, a real 2D array, is Zonal wind component at the lowest model level [m/s].
## Atm%v_bot
Atm%v_bot, a real 2D array, is Meridional wind component at the lowest model level [m/s].
## Atm%p_surf
Atm%p_surf, a real 2D array, is Surface pressure [Pa].
## Atm%slp
Atm%slp, a real 2D array, is Sea-level pressure [Pa].
## Atm%gust
Atm%gust, a real 2D array, is Gustiness factor — a minimum wind speed added in quadrature to the resolved wind to account for sub-grid convective gusts in surface flux calculations [m/s].
## Atm%coszen
Atm%coszen, a real 2D array, is Cosine of the solar zenith angle; used to weight shortwave fluxes and partition direct vs. diffuse radiation [dimensionless].

## Atm%flux_sw
Atm%flux_sw, a real 2D array, is Total net shortwave flux at the surface (absorbed by the surface) [W/m²].
## Atm%flux_sw_dir
Atm%flux_sw_dir, a real 2D array, is Direct-beam component of the net shortwave flux [W/m²].
## Atm%flux_sw_dif
Atm%flux_sw_dif, a real 2D array, is Diffuse component of the net shortwave flux [W/m²].
## Atm%flux_sw_down_vis_dir
Atm%flux_sw_down_vis_dir, a real 2D array, is Downward direct-beam flux in the visible band (0.2–0.7 µm) [W/m²].
## Atm%flux_sw_down_vis_dif
Atm%flux_sw_down_vis_dif, a real 2D array, is Downward diffuse flux in the visible band [W/m²].
## Atm%flux_sw_down_total_dir
Atm%flux_sw_down_total_dir, a real 2D array, is Downward direct-beam broadband shortwave flux [W/m²].
## Atm%flux_sw_down_total_dif
Atm%flux_sw_down_total_dif, a real 2D array, is Downward diffuse broadband shortwave flux [W/m²].
## Atm%flux_sw_vis
Atm%flux_sw_vis, a real 2D array, is Net (downward minus reflected) visible-band shortwave flux at the surface [W/m²].
## Atm%flux_sw_vis_dir
Atm%flux_sw_vis_dir, a real 2D array, is Direct-beam component of the net visible shortwave flux [W/m²].
## Atm%flux_sw_vis_dif
Atm%flux_sw_vis_dif, a real 2D array, is Diffuse component of the net visible shortwave flux [W/m²].
## Atm%flux_lw
Atm%flux_lw, a real 2D array, is Net downward longwave flux at the surface [W/m²].
## Atm%lprec
Atm%lprec, a real 2D array, is Mass of liquid precipitation accumulated since the last time step [kg/m²]; equivalent to a rate in kg/m²/s when divided by dt_atm.
## Atm%fprec
Atm%fprec, a real 2D array, is Mass of frozen (solid) precipitation accumulated since the last time step [kg/m²].
## Atm%gex_atm2lnd
Atm%gex_atm2lnd, a real 3D array, is Generic exchange fields sent from the atmosphere to the land model; the field table defines which quantities are exchanged (e.g. CO₂, aerosol deposition); third dimension indexes the exchange field list.
## Atm%gex_lnd2atm
Atm%gex_lnd2atm, a real 3D array, is Generic exchange fields returned from the land model to the atmosphere (e.g. surface emission fluxes); third dimension indexes the exchange field list.
## Atm%fields
Atm%fields, type(coupler_2d_bc_type), is Array of additional tracer boundary-condition fields used for atmosphere-ocean gas exchange (CO₂, O₂, CFCs, etc.); registered and populated by atmos_tracer_flux_init.

## Atm%Time
Atm%Time, type(time_type), is Current model time; passed to diag_manager send_data calls and to fms_data_override.
## Atm%Time_step
Atm%Time_step, type(time_type), is Atmospheric model timestep duration.
## Atm%Time_init
Atm%Time_init, type(time_type), is Reference (initial) time for the model run.
## Atm%pelist
Atm%pelist, an integer 1D array, is List of MPI PE numbers on which the atmosphere is running.
## Atm%pe
Atm%pe, logical, is .true. on PEs that are part of the atmosphere pelist; used to gate atmosphere-only code blocks.

## Surf_diff
`Surf_diff` is of type `surf_diff_type` defined in atmos_phys/atmos_phys/atmos_param/vert_diff/vert_diff.F90. 
It carries the forward-elimination coefficients from the implicit vertical diffusion scheme that couples the atmosphere to the surface models.Surf_diff is a component of Atm (for example Atm%surf_diff%dtmass).

## Atm%Surf_diff%dtmass
Atm%Surf_diff%dtmass, a real 2D array, is dt/mass — ratio of the atmospheric timestep to the surface-layer air mass [s·m²/kg]; scales flux tendencies to temperature/tracer tendencies.
## Atm%Surf_diff%dflux_t
Atm%Surf_diff%dflux_t, a real 2D array, is d(sensible heat flux)/d(T_surf) — linearisation of the surface heat flux with respect to surface temperature; used to form the implicit coupling term [W/m²/K].
## Atm%Surf_diff%delta_t
Atm%Surf_diff%delta_t, a real 2D array, is Forward-elimination coefficient for temperature from the implicit tridiagonal scheme; represents the accumulated atmospheric temperature forcing at the bottom level waiting for the surface response [K].
## Atm%Surf_diff%delta_u
Atm%Surf_diff%delta_u, a real 2D array, is Forward-elimination coefficient for zonal wind from the implicit scheme [m/s].
## Atm%Surf_diff%delta_v
Atm%Surf_diff%delta_v, a real 2D array, is Forward-elimination coefficient for meridional wind from the implicit scheme [m/s].
## Atm%Surf_diff%dflux_tr
Atm%Surf_diff%dflux_tr, a real 3D array, is d(tracer flux)/d(tracer_surf) — linearisation of tracer surface fluxes with respect to surface tracer concentration; third dimension indexes tracers.
## Atm%Surf_diff%delta_tr
Atm%Surf_diff%delta_tr, a real 3D array, is Forward-elimination coefficient for each tracer from the implicit scheme; third dimension indexes tracers.
## Atm%Surf_diff%tdt_dyn
Atm%Surf_diff%tdt_dyn, a real 3D array, is Temperature tendency from dynamics (advection, etc.) passed through the diffusion scheme.
## Atm%Surf_diff%qdt_dyn
Atm%Surf_diff%qdt_dyn, a real 3D array, is Moisture tendency from dynamics.
## Atm%Surf_diff%dgz_dyn
Atm%Surf_diff%dgz_dyn, a real 3D array, is Geopotential height tendency from dynamics.
## Atm%Surf_diff%ddp_dyn
Atm%Surf_diff%ddp_dyn, a real 3D array, is Pressure-thickness tendency from dynamics.
## Atm%Surf_diff%tdt_rad
Atm%Surf_diff%tdt_rad, a real 3D array, is Temperature tendency from radiation; used in the MIZ (marginal ice zone) forecast mode.

## Grid geometry subfields
Atm%grid is of type grid_box_type defined in FMS/exchange/xgrid.
It holds the geometric quantities needed for second-order conservative flux
remapping between the atmosphere and surface component grids on a
cubic-sphere mesh.  grid is a component of Atm(for example, Atm%grid%dx)

## Atm%grid%dx
Atm%grid%dx, a real 2D array, is Grid-box width in the x-direction [m].
## Atm%grid%dy
Atm%grid%dy, a real 2D array, is Grid-box width in the y-direction [m].
## Atm%grid%area
Atm%grid%area, a real 2D array, is Grid-box area [m²].
## Atm%grid%edge_w
Atm%grid%edge_w, real 1D array, is Western edge lengths of grid boxes along the boundary [m].
## Atm%grid%edge_e
Atm%grid%edge_e, real 1D array, is Eastern edge lengths [m].
## Atm%grid%edge_s
Atm%grid%edge_s, real 1D array, is Southern edge lengths [m].
## Atm%grid%edge_n
Atm%grid%edge_n, real 1D array, is Northern edge lengths [m].
## Atm%grid%en1
Atm%grid%en1, a real 3D array, is First unit normal vector at grid-box edges; used to project vector fields (winds, stresses) during remapping.
## Atm%grid%en2
Atm%grid%en2, a real 3D array, is Second unit normal vector at grid-box edges.
## Atm%grid%vlon
Atm%grid%vlon, a real 3D array, is Unit vector in the local longitude direction at each grid point; used to rotate between geographic and local coordinate frames during exchange.
## Atm%grid%vlat
vlat, a real 3D array, is Unit vector in the local latitude direction at each grid point.

