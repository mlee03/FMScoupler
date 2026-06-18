# ice_data_type
ice_data_type is the publicly visible face of the SIS2 sea-ice model. 
It carries the surface state, inter-component fluxes, and PE/domain metadata that the coupler needs.
 All internal SIS2 state is accessed through the private control structures fCS (fast ice) and sCS (slow ice), 
 which are opaque to the coupler. 

The type supports a split fast-ice / slow-ice architecture. Fast processes (surface thermodynamics, atmosphere–ice flux coupling) run on the atmospheric timestep on atmosphere PEs; slow processes (ice dynamics, freezing/melting, transport) run on the coupled timestep and may run on ocean PEs when slow_ice_with_ocean = .true.. See [README-new.md](README-new.md) for a description of how this affects the PE layout.

## Ice%pe
Ice%pe, logical, is .true. on any PE that participates in ice computation (fast or slow).
## Ice%fast_ice_pe
Ice%fast_ice_pe, logical, is .true. on PEs in the fast-ice pelist; these PEs handle atmosphere-timestep ice processes.
## Ice%slow_ice_pe
Ice%slow_ice_pe, logical, is .true. on PEs in the slow-ice pelist; these PEs handle coupled-timestep ice dynamics and ice-ocean exchange.
## Ice%shared_slow_fast_PEs
Ice%shared_slow_fast_PEs, logical, is .true. when fast and slow ice use the same PE set and domain decomposition (i.e. slow_ice_with_ocean=.false.); .false. when slow ice runs on the ocean PEs.
## Ice%pelist
Ice%pelist, integer(:), is Combined list of all ice PEs (union of fast and slow pelists); used for flux exchange.
## Ice%fast_pelist
Ice%fast_pelist, integer(:), is MPI PE numbers for fast-ice processes.
## Ice%slow_pelist
Ice%slow_pelist, integer(:), is MPI PE numbers for slow-ice processes.
## Ice%Domain
Ice%Domain, type(domain2D), is Copy of the fast-ice FMS domain without halos; used for exchange-grid setup.
## Ice%slow_Domain_NH
Ice%slow_Domain_NH, type(domain2D), is Copy of the slow-ice FMS domain without halos; used for ice-ocean flux redistribution.
## Ice%fast_domain
Ice%fast_domain, type(domain2D) pointer, is Pointer to the fast-ice MPI domain (or an allocated copy on slow-ice PEs).
## Ice%slow_domain
Ice%slow_domain, type(domain2D) pointer, is Pointer to the slow-ice MPI domain (or an allocated copy on fast-ice PEs).
## Ice%ocean_pt
Ice%ocean_pt, logical(:,:), is Mask array; .true. at ocean (non-land) points.
## Ice%xtype
Ice%xtype, integer, is Transfer mode for ice-ocean flux exchange: DIRECT (3) when ice and ocean share the same decomposition, REDIST (2) when they differ.
## Ice%axes
Ice%axes, integer(3), is Diag-manager axis IDs for the ice surface grid.
## Ice%Time
Ice%Time, type(time_type), is The sea-ice model's current clock time.

## Ice%part_size
Ice%part_size, a real 3D array, is Fractional coverage of the grid cell by each ice thickness category [dimensionless, 0–1]; category 1 is open ocean.
Dimensioned (:, :, n_categories) and provide per-ice-thickness-category information to the atmosphere each fast timestep. Category 1 is open ocean; the sum of part_size over all categories equals 1.
## Ice%t_surf
Ice%t_surf, a real 3D array, is Surface temperature of the ocean (category 1) or each ice thickness category [K].
Dimensioned (:, :, n_categories) and provide per-ice-thickness-category information to the atmosphere each fast timestep. Category 1 is open ocean; the sum of part_size over all categories equals 1.
## Ice%albedo
Ice%albedo, a real 3D array, is Broadband surface albedo averaged across all wavelength and orientation bands within each ice category [dimensionless, 0–1].
Dimensioned (:, :, n_categories) and provide per-ice-thickness-category information to the atmosphere each fast timestep. Category 1 is open ocean; the sum of part_size over all categories equals 1.
## Ice%albedo_vis_dir
Ice%albedo_vis_dir, a real 3D array, is Surface albedo for direct visible shortwave radiation in each ice category [dimensionless].
Dimensioned (:, :, n_categories) and provide per-ice-thickness-category information to the atmosphere each fast timestep. Category 1 is open ocean; the sum of part_size over all categories equals 1.
## Ice%albedo_nir_dir
Ice%albedo_nir_dir, a real 3D array, is Surface albedo for diffuse visible shortwave radiation in each ice category [dimensionless].
Dimensioned (:, :, n_categories) and provide per-ice-thickness-category information to the atmosphere each fast timestep. Category 1 is open ocean; the sum of part_size over all categories equals 1.
## Ice%albedo_vis_dif
Ice%albedo_vis_dif, a real 3D array, is Surface albedo for direct near-infrared shortwave radiation in each ice category [dimensionless].
Dimensioned (:, :, n_categories) and provide per-ice-thickness-category information to the atmosphere each fast timestep. Category 1 is open ocean; the sum of part_size over all categories equals 1.
## Ice%albedo_nir_dif
Ice%albedo_nir_dif, a real 3D array, is Surface albedo for diffuse near-infrared shortwave radiation in each ice category [dimensionless].
Dimensioned (:, :, n_categories) and provide per-ice-thickness-category information to the atmosphere each fast timestep. Category 1 is open ocean; the sum of part_size over all categories equals 1.
## Ice%rough_mom
Ice%rough_mom, a real 3D array, is Surface roughness length for momentum at the ocean/ice surface, as provided by ocean_rough_mod [m].
Dimensioned (:, :, n_categories) and provide per-ice-thickness-category information to the atmosphere each fast timestep. Category 1 is open ocean; the sum of part_size over all categories equals 1.
## Ice%rough_heat
Ice%rough_heat, a real 3D array, is Surface roughness length for heat at the ocean/ice surface [m].
Dimensioned (:, :, n_categories) and provide per-ice-thickness-category information to the atmosphere each fast timestep. Category 1 is open ocean; the sum of part_size over all categories equals 1.
## Ice%rough_moist
Ice%rough_moist, a real 3D array, is Surface roughness length for moisture at the ocean/ice surface [m].
Dimensioned (:, :, n_categories) and provide per-ice-thickness-category information to the atmosphere each fast timestep. Category 1 is open ocean; the sum of part_size over all categories equals 1.
## Ice%u_surf
Ice%u_surf, a real 3D array, is Eastward surface velocity of the ocean (category 1, :,:,1) or sea ice [m/s]; used as a lower boundary condition for wind stress.  Dimensioned (:, :, n_categories) and provide per-ice-thickness-category information to the atmosphere each fast timestep. Category 1 is open ocean; the sum of part_size over all categories equals 1.
## Ice%v_surf
Ice%v_surf, a real 3D array, is Northward surface velocity of the ocean (category 1, :,:,1) or sea ice [m/s].
Dimensioned (:, :, n_categories) and provide per-ice-thickness-category information to the atmosphere each fast timestep. Category 1 is open ocean; the sum of part_size over all categories equals 1.
## Ice%flux_uv_stagger
Ice%flux_uv_stagger, integer, is Staggering of u_surf/v_surf relative to tracer points; valid values are AGRID, BGRID_NE, CGRID_NE, BGRID_SW, CGRID_SW (Arakawa notation). Initialized to -999 so that a global max across all PEs propagates the correct value to non-ice PEs.  Dimensioned (:, :, n_categories) and provide per-ice-thickness-category information to the atmosphere each fast timestep. Category 1 is open ocean; the sum of part_size over all categories equals 1.

## Ice%s_surf
Ice%s_surf, a real 2D array, is Ocean surface salinity [g salt / kg seawater]; populated by flux_ocean_to_ice.
## Ice%SST_C
Ice%SST_C, a real 2D array, is Ocean surface temperature [°C]; used in forecast mode.


## Ice%area
Ice%area, a real 2D array, is Area of each ocean cell [m**2]; land cells have area = 0 and this field can double as a mask.
## Ice%mi
Ice%mi, a real 2D array, is Total ice + snow mass per unit area [kg/m**2]; passed to the ocean for pressure loading and to the wave model.

## Ice%flux_u
Ice%flux_u, a real 2D array, is Flux of x-momentum into the ocean (zonal wind/ice stress) [Pa].  
It is computed by the slow-ice model.
## Ice%flux_v
Ice%flux_v, a real 2D array, is Flux of y-momentum into the ocean (meridional wind/ice stress) [Pa].
It is computed by the slow-ice model.
## Ice%flux_t
Ice%flux_t, a real 2D array, is Sensible heat flux out of the ocean W/m**2.
It is computed by the slow-ice model.
## Ice%flux_q
Ice%flux_q, a real 2D array, is Evaporative moisture flux out of the ocean [kg/m**2/s].
It is computed by the slow-ice model.
## Ice%flux_lh
Ice%flux_lh, a real 2D array, is Latent heat flux out of the ocean W/m**2.
It is computed by the slow-ice model.
## Ice%flux_lw
Ice%flux_lw, a real 2D array, is Net longwave flux out of the ocean W/m**2.
It is computed by the slow-ice model.
## Ice%flux_sw_vis_dir
Ice%flux_sw_vis_dir, a real 2D array, is Direct visible shortwave heat flux into the ocean W/m**2.
It is computed by the slow-ice model.
## Ice%flux_sw_vis_dif
Ice%flux_sw_vis_dif, a real 2D array, is Diffuse visible shortwave heat flux into the ocean W/m**2.
It is computed by the slow-ice model.
## Ice%flux_sw_nir_dir
Ice%flux_sw_nir_dir, a real 2D array, is Direct near-infrared shortwave heat flux into the ocean W/m**2.
It is computed by the slow-ice model.
## Ice%flux_sw_nir_dif
Ice%flux_sw_nir_dif, a real 2D array, is Diffuse near-infrared shortwave heat flux into the ocean W/m**2.
It is computed by the slow-ice model.
## Ice%flux_salt
Ice%flux_salt, a real 2D array, is Salt flux out of the ocean (brine rejection / melting) [kg/m**2/s].
It is computed by the slow-ice model.
## Ice%salt_left_behind
Ice%salt_left_behind, a real 2D array, is Salt remaining in the ocean during ice growth (not ejected as brine) [kg/m**2/s].
It is computed by the slow-ice model.
## Ice%lprec
Ice%lprec, a real 2D array, is Liquid precipitation flux into the ocean [kg/m**2].
It is computed by the slow-ice model.
## Ice%fprec
Ice%fprec, a real 2D array, is Frozen precipitation flux into the ocean [kg/m**2].
It is computed by the slow-ice model.
## Ice%runoff
Ice%runoff, a real 2D array, is Liquid runoff from land into the ocean [kg/m**2].
It is computed by the slow-ice model.
## Ice%runoff_hflx
Ice%runoff_hflx, a real 2D array, is Heat flux associated with liquid runoff, relative to a reference temperature W/m**2.
It is computed by the slow-ice model.
## Ice%runoff_carbon
Ice%runoff_carbon, a real 2D array, is Carbon content of liquid runoff into the ocean [kg/m**2].
It is computed by the slow-ice model.
## Ice%calving
Ice%calving, a real 2D array, is Calving of ice or frozen freshwater runoff into the ocean [kg/m**2].
It is computed by the slow-ice model.
## Ice%calving_hflx
Ice%calving_hflx, a real 2D array, is Heat flux associated with calving, relative to a reference temperature W/m**2.
It is computed by the slow-ice model.
## Ice%p_surf
Ice%p_surf, a real 2D array, is Pressure at the ocean surface [Pa]; may or may not include atmospheric pressure depending on configuration.
It is computed by the slow-ice model.
## Ice%stress_mag
Ice%stress_mag, a real 2D array, is Time-mean magnitude of the ice-ocean stress [Pa]; passed to the ocean when pass_stress_mag = .true. in SIS_slow_CS.
It is computed by the slow-ice model.

## Ice%ustar_berg
Ice%ustar_berg, a real 2D array, is Friction velocity contribution beneath icebergs [m/s]; used for iceberg-ocean drag.
## Ice%area_berg
Ice%area_berg, a real 2D array, is Fraction of the grid cell covered by icebergs [m**2/m**2].
## Ice%mass_berg
Ice%mass_berg, a real 2D array, is Mass of icebergs per unit area [kg/m**2].

## Ice%ocean_fields
Ice%ocean_fields, type(coupler_3d_bc_type), is Named surface-state fields shared with the atmosphere (SST, SSS, piston velocities, etc.) for atmosphere-ocean gas exchange; populated by flux_ocean_to_ice.
## Ice%ocean_fluxes
Ice%ocean_fluxes, type(coupler_2d_bc_type), is Computed gas fluxes from the ice to the ocean for additional tracers.
## Ice%ocean_fluxes_top
Ice%ocean_fluxes_top, type(coupler_3d_bc_type), is Gas flux boundary conditions at the top of the ice (atmosphere side); archaic and flagged for eventual removal.

## Ice%fCS
Ice%fCS, type(SIS_fast_CS) pointer, is Control structure for the SIS2 fast ice thermodynamics; lives on atmosphere PEs; contains the fast-ice grid (fCS%G), diagnostics, and all fast-timestep state.  This is a private field.
## Ice%sCS
Ice%sCS, type(SIS_slow_CS) pointer, is Control structure for the SIS2 slow ice dynamics and thermodynamics; may live on ocean PEs when slow_ice_with_ocean=.true.; contains the slow-ice grid (sCS%G), dynamics, tracer registry, and all slow-timestep state.  This is a private field.
## Ice%icebergs
Ice%icebergs, type(icebergs) pointer, is Control structure for the Lagrangian iceberg module; null when do_icebergs=.false..
This is a private field.
## Ice%US
Ice%US, type(unit_scale_type) pointer, is SIS2 dimensional unit-scaling factors; used to convert between external MKS values and SIS2 internal non-dimensionalized units.  This is a private field.
## Ice%Ice_restart
Ice%Ice_restart, type(SIS_restart_CS) pointer, is Control structure for writing and reading slow-ice restart files.
This is a private field.
## Ice%Ice_fast_restart
Ice%Ice_fast_restart, type(SIS_restart_CS) pointer, is Control structure for writing and reading fast-ice restart files.
This is a private field.
## Ice%OBC
Ice%OBC, type(ice_OBC_type) pointer, is Control structure for ice open boundary conditions; null when OBCs are not configured.
This is a private field.
## Ice%restart_output_dir
Ice%restart_output_dir, character(240), is Directory path for restart file output; default is './RESTART/'.
This is a private field.
