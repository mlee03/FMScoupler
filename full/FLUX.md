# Flux Exchange

*Original authors: Bruce Wyman, V. Balaji, Sergey Malyshev*

## Background information

Six modules couple the atmosphere, ocean, land, and sea-ice components with flux exchange in the FMSCoupler full coupler.

| Module | Role |
|---|---|
| `atm_land_ice_flux_exchange` | Exchanges fluxes between atmosphere, land, and ice; registers all diagnostic fields. |
| `flux_exchange` | Top-level module that initializes all flux-exchange modules; contains subroutines for stock computation. |
| `atmos_ocean_fluxes_calc` | Computes non-deposition gas fluxes between atmosphere and ocean. |
| `atmos_ocean_dep_fluxes_calc` | Computes deposition gas fluxes between atmosphere and ocean. |
| `ice_ocean_flux_exchange` | Exchanges fluxes between ice and ocean. |
| `land_ice_flux_exchange` | Exchanges fluxes between land and ice. |

GFDL coupled models represent atmosphere and land on the same cubed-sphere grid; the land grid, however, is masked 
(for cells containing ice or water) and data is stored as arrays of rank 1 while atm data is stored as arrays of rank 2.  
Ice and ocean must share the same physical grid, though their MPI domain decompositions may differ.  
The masked region of the land grid and the ice-ocean grid must tile each other such that every atmosphere grid 
cell is covered by either land or ice-ocean, but not both.  The masked regions of the ice and ocean grids must be identical.

The three component grids tile the sphere. `|xxx|` marks a masked (inactive) grid point:

```
ATMOSPHERE  |----|----|----|----|----|----|----|----|
      LAND  |---|---|---|---|xxx|xxx|xxx|xxx|xxx|xxx|
       ICE  |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
     OCEAN  |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
```

## Useful information

* Each component model must have a public defined data type containing specific
  boundary fields. A list of these quantities is located in the NOTES of this document.

* The surface flux of sensible heat and surface evaporation can be implicit functions
  of surface temperature. As a consequence, the parts of the land and sea-ice models
  that update the surface temperature must be called on the atmospheric time step

* The surface fluxes of all other tracers and momentum are assumed to be explicit
  functions of all surface parameters.

* While no explicit reference is made within this module to the implicit treatment
  of vertical diffusion in the atmosphere and in the land or sea-ice models, the
  module is designed to allow for simultaneous implicit time integration on both
  sides of the surface interface.

* Due to the previous point, the diffusion part of the land and ice models must be called on the
  atmospheric time step, although in the case of concurrent-ice coupling, this
  version of the sea-ice that is called by the atmosphere may later be replaced
  by a version of the ice that is tightly coupled with the ocean.

* The fluxes of additional tracers related to biological quantities or the
  air-sea exchange of gases are accomplished by specifying fields that will
  be passed between components via the "field_table" and the use of named
  fields in the coupler_..._bc_types.

* Any field passed from one component to another may be "faked" to a constant
  value, or to data acquired from a file, using the data_override feature of FMS.
  The fields to override are runtime configurable, using the text file
  data_table for input. See the data_override_mod documentation for more
  details.  It is NOT RECOMMENDED to exercise the data override capabilities 
  of the FMS coupler until the user has acquired considerable sophistication in running FMS.

* model1_model2_boundary_type (e.g., atmos_land_boundary_type) contains fields that model2 (land) gets
  from model1 (atm), may also include fluxes. These are declared by flux_exchange_mod and have private components. 


## Exchange Modes

* The atmosphere, land, and ice grids exchange information via xmap_sfc exchange grid:  first, data is mapped
  onto the exchange grid, then computation is carried out on the exchange grid before data is interpolated onto the 
  receiving component grid.

* The land and ice grids exchange runoff data using the exchange grid xmap_runoff for conservative interpolation

* Transfer of data between the ice bottom and ocean does not require an exchange
  grid as the grids are physically identical. The flux routines will automatically
  detect and redistribute data if their domain decompositions are different, or copy data if 
  the domain decomposition is identical

* To get information from the atmosphere to the ocean it must pass through the
  ice model, first by interpolating from the atmospheric grid to the ice grid,
  and then transferring from the ice grid to the ocean grid.


## Namelist Parameters (`flux_exchange_nml`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `z_ref_heat` | real | 2.0 m | Reference height for temperature and relative-humidity diagnostics (`t_ref`, `rh_ref`, `del_h`, `del_q`). |
| `z_ref_mom` | real | 10.0 m | Reference height for momentum diagnostics (`u_ref`, `v_ref`, `del_m`). |
| `do_area_weighted_flux` | logical | `.false.` | When `.true.`, fluxes passed to the ocean are multiplied by the ice area fraction before redistribution. |
| `debug_stocks` | logical | `.false.` | Enables additional stock-conservation debug output when `.true.`. |
| `divert_stocks_report` | logical | `.false.` | If `.true.`, output stock reporting to `stocks.out` instead to stdout. |
| `do_runoff` | logical | `.true.` | Enables interpolation of land runoff to the ocean. |
| `do_forecast` | logical | `.false.` | Do forecast if `.true.`|
| `nblocks` | integer | 1 | Number of OpenMP blocks for parallel computation on `sfc_xgrid` |
| `partition_fprec_from_lprec` | logical | `.false.` | If `.true.`, converts liquid precipitation to snow when `t_ref < tfreeze`. |
| `scale_precip_2d` | logical | `.false.` | Rescales `lprec` by a 2-D field read from `data_table`. |


## Data Override Capabilities
Any field passed from one component to another may be "faked" to a constant value, or to data acquired from a file, using the
data_override feature of FMS. The fields to override are runtime configurable, using data_table for input.
See the data_override_mod documentation for more details.

The original authors DO NOT RECOMMEND exercising the data override capabilities of the FMS coupler until the user has acquired considerable
sophistication in running FMS.

| Subroutine | Transfer | Overridable fields |
|---|---|---|
| `sfc_boundary_layer` | Atm boundary -> exchange grid | `t_bot`, `z_bot`, `p_bot`, `u_bot`, `v_bot`, `p_surf`, `slp`, `gust`, fields in the coupler bc type |
| `sfc_boundary_layer` | Ice boundary -> exchange grid | `t_surf`, `rough_mom`, `rough_heat`, `rough_moist`, `albedo`, `albedo_vis_dir`, `albedo_nir_dir`, `albedo_vis_dif`, `albedo_nir_dif`, `u_surf`, `v_surf` |
| `sfc_boundary_layer` | Land boundary -> exchange grid | `t_surf`, `t_ca`, `rough_mom`, `rough_heat`, `albedo`, `tracers`, `albedo_vis_dir`, `albedo_nir_dir`, `albedo_vis_dif`, `albedo_nir_dif` |
| `sfc_boundary_layer` | Exchange grid -> `Land_ice_atmos_boundary` | `t`, `albedo`, `albedo_vis_dir`, `albedo_nir_dir`, `albedo_vis_dif`, `albedo_nir_dif`, `land_frac`, `dt_t`, `dt_tr`, `u_flux`, `v_flux`, `dtaudu`, `dtaudv`, `u_star`, `b_star`, `rough_mom` |
| `flux_down_from_atmos` | Atmosphere boundary -> exchange grid | `flux_sw`, `flux_sw_dir`, `flux_sw_dif`, `flux_sw_down_vis_dir`, `flux_sw_down_vis_dif`, `flux_sw_down_total_dir`, `flux_sw_down_total_dif`, `flux_sw_vis`, `flux_sw_vis_dir`, `flux_sw_vis_dif`, `flux_lw`, `lprec`, `frac_precip`, `fprec`, `coszen`, `dtmass`, `delta_t`, `dflux_t`, `delta_tr`, `dflux_tr` |
| `flux_down_from_atmos` | Exchange grid -> land boundary | `drag_q`, `lwdn_flux`, `cd_m`, `cd_t`, `bstar`, `ustar`, `wind`, `z_bot`, `t_flux`, `lw_flux`, `sw_flux`, `sw_flux_down_vis_dir`, `sw_flux_down_total_dir`, `sw_flux_down_vis_dif`, `sw_flux_down_total_dif`, `lprec`, `fprec`, `dhdt`, `drdt`, `p_surf`, `tr_flux`, `dfdtr` |
| `flux_down_from_atmos` | Exchange grid -> ice boundary | `u_flux`, `v_flux`, `t_flux`, `q_flux`, `lw_flux`, `sw_flux_nir_dir`, `sw_flux_vis_dir`, `sw_flux_nir_dif`, `sw_flux_vis_dif`, `sw_down_vis_dir`, `sw_down_vis_dif`, `sw_down_nir_dir`, `sw_down_nir_dif`, `lprec`, `fprec`, `dhdt`, `dedt`, `drdt`, `coszen`, `p` |
| `flux_up_to_atmos` | Ice boundary -> atmosphere boundary | `t_surf` |
| `flux_up_to_atmos` | Land boundary -> atmosphere boundary | `t_ca`, `t_surf`, `tr` |
| `flux_land_to_ice` | Land boundary -> ice boundary | `runoff`, `calving`, `runoff_hflx`, `calving_hflx` |
| `flux_ice_to_ocean` | Ice boundary -> ocean boundary | `u_flux`, `v_flux`, `t_flux`, `q_flux`, `salt_flux`, `lw_flux`, `sw_flux_nir_dir`, `sw_flux_nir_dif`, `sw_flux_vis_dir`, `sw_flux_vis_dif`, `lprec`, `fprec`, `runoff`, `calving`, `runoff_hflx`, `calving_hflx`, `p`, `mi`, `ustar_berg`, `area_berg`, `mass_berg` |
| `flux_ocean_to_ice` | Ocean boundary -> ice boundary | `u`, `v`, `t`, `s`, `frazil`, `sea_level` |

---

## Diagnostic Fields

All fields below are registered in `atm_land_ice_flux_exchange.F90` inside `diag_field_init`.

| Category | Registration/Notes | Fields |
|---|---|---|
| Static fields | Registered in `diag_field_init` | `land_mask`, `height2m`, `height10m`, `sftlf` |
| Atmosphere surface fields | Registered in `diag_field_init` | `ice_mask`, `wind`, `drag_moist`, `drag_heat`, `drag_mom`, `rough_moist`, `rough_heat`, `rough_mom`, `u_star`, `b_star`, `q_star`, `thv_atm`, `thv_surf`, `tau_x`, `tau_y`, `t_ocean`, `t_surf`, `t_ca`, `z_atm`, `p_atm`, `slp`, `gust`, `shflx`, `lwflx`, `t_atm`, `u_atm`, `v_atm`, `t_ref`, `rh_ref`, `rh_ref_cmip`, `u_ref`, `v_ref`, `wind_ref`, `del_h`, `del_m`, `del_q`, `q_ref`, `rough_scale`, `evap`, `co2_bot` |
| Atmosphere tracer fields (per-tracer, name-prefixed) | Pattern-based tracer diagnostics | `*_tot_con_atm`, `*_tot_con_ref`, `*_atm`, `*_surf`, `*_flux`, `*_ref`, `*_mol_flux`, `*_atm_dvmr`, `*_surf_dvmr`, `*_mol_flux_atm0` |
| CMIP fields | Registered with `register_cmip_diag_field_2d` or `fms_diag_register_diag_field` with `use_AM3_physics` | `tas`, `uas`, `vas`, `sfcWind`, `huss`, `hurs`, `rhs`, `ts`, `psl`, `tauu`, `tauv`, `hfss`, `hfls`, `evspsbl`, `tslsi`, `tos`, `sic` |
| Global scalar time-series fields | Registered with `register_global_diag_field`, only without `use_AM3_physics` | `evspsbl`, `ts`, `tas`, `tasl`, `hfss`, `hfls`, `rls` |
| Land axes fields | Registered with `register_tiled_diag_field` or `fms_diag_register_diag_field` with `_USE_LEGACY_LAND_` | `t_ref`, `q_ref`, `rh_ref`, `u_ref`, `v_ref`, `evap`, `shflx`, `tasLut`, `hussLut`, `*_tot_con_atm`, `*_tot_con_ref`, `*_flux`, `*_mol_flux`, `*_ref` |

---

## Required Variables in Component Data Types

The following fields must be defined in each component's public data type for flux exchange.

!! type (atmos_boundary_data_type) :: Atm
!!
!! real, dimension(:) :: Atm%lon_bnd & ! longitude axis grid box boundaries in radians
!!                                     ! must be monotonic
!!                       Atm%lat_bnd   ! latitude axis grid box boundaries in radians
!!                                     ! must be monotonic
!! real, dimension(:,:) :: Atm%t_bot   & ! temperature at lowest model level
!!                         Atm%q_bot   & ! specific humidity at lowest model level
!!                         Atm%z_bot   & !    height above the surface for the lowest model level (m)
!!                         Atm%p_bot   & !    pressure at lowest model level (pa)
!!                         Atm%u_bot   & !    zonal wind component at lowest model level (m/s)
!!                         Atm%v_bot   & !    meridional wind component at lowest model level (m/s)
!!                         Atm%p_surf  & !   surface pressure (pa)
!!                         Atm%slp     & !   sea level pressure (pa)
!!                         Atm%gust    & !   gustiness factor (m/s)
!!                         Atm%flux_sw & !  net shortwave flux at the surface
!!                         Atm%flux_lw & !  downward longwave flux at the surface
!!                         Atm%lprec   & !  liquid precipitation (kg/m2)
!!                         Atm%fprec   & !  water equivalent frozen precipitation (kg/m2)
!!                         Atm%coszen  & !  cosine of the zenith angle
!! integer, dimension(4) :: Atm%axes ! Axis identifiers returned by diag_axis_init for the
!!                                   ! atmospheric model axes: X, Y, Z_full, Z_half.

!! The following five fields are gathered into a data type for convenience in passing
!! this information through the different levels of the atmospheric model --
!! these fields are rlated to the simultaneous implicit time steps in the
!! atmosphere and surface models -- they are described more fully in
!! flux_exchange.tech.ps and in the documntation for vert_diff_mod
!!
!! ~~~~~~~~~~{.f90}
!! type (surf_diff_type) :: Atm%Surf_Diff
!!
!! real, dimension(:,:) :: Atm%Surf_Diff%dtmass  & !dt/mass where dt=atmospheric time step ((i+1)=(i-1) for leapfrog)(s)
!!                                                 ! mass = mass per unit area of lowest atmosphehic layer  (Kg/m2))
!!                         Atm%Surf_Diff%delta_t & ! increment ((i+1) = (i-1) for leapfrog) in temperature of
!!                                                 ! lowest atmospheric layer  (K)
!!                         Atm%Surf_Diff%delta_q & ! increment ((i+1) = (i-1) for leapfrog) in specific humidity of
!!                                                 ! lowest atmospheric layer (nondimensional -- Kg/Kg)
!!                         Atm%Surf_Diff%dflux_t & ! derivative of implicit part of downward temperature flux at top of
!!                                                 ! lowest atmospheric layer with respect to temperature
!!                                                 ! of lowest atmospheric layer (Kg/(m2 s))
!!                         Atm%Surf_Diff%dflux_q   ! derivative of implicit part of downward moisture flux at top of
!!                                                 ! lowest atmospheric layer with respect to specific humidity of
!!                                                 ! of lowest atmospheric layer (Kg/(m2 s))

!! type (land_boundary_data_type) :: Land
!!
!! real, dimension(:) :: Land%lon_bnd & ! longitude axis grid box boundaries in radians
!!                                      ! must be monotonic
!!                       Land%lat_bnd   ! latitude axis grid box boundaries in radians
!!                                      ! must be monotonic
!!
!! logical, dimension(:,:,:) :: Land%mask & ! land/sea mask (true for land)
!!                              Land%glacier ! glacier mask  (true for glacier)
!!
!! real, dimension(:,:,:) :: Land%tile_size  & !  fractional area of each tile (partition)
!!                           Land%t_surf     & ! surface temperature (deg k)
!!                           Land%albedo     & ! surface albedo (fraction)
!!                           Land%rough_mom  & ! surface roughness for momentum (m)
!!                           Land%rough_heat & ! surface roughness for heat/moisture (m)
!!                           Land%stomatal   & ! stomatal resistance
!!                           Land%snow       & ! snow depth (water equivalent) (kg/m2)
!!                           Land%water      & ! water depth of the uppermost bucket (kg/m2)
!!                           Land%max_water    ! maximum water depth allowed in the uppermost bucket (kg/m2)
!! ~~~~~~~~~~
!! type (ice_boundary_data_type) :: Ice
!!
!! real, dimension(:) :: Ice%lon_bnd    & ! longitude axis grid box boundaries for temperature points
!!                                        ! in radians (must be monotonic)
!!                       Ice%lat_bnd    & ! latitude axis grid box boundaries for temperature points
!!                                        ! in radians (must be monotonic)
!!                       Ice%lon_bnd_uv & ! longitude axis grid box boundaries for momentum points
!!                                        ! in radians (must be monotonic)
!!                       Ice%lat_bnd_uv   ! latitude axis grid box boundaries for momentum points
!!                                        ! in radians (must be monotonic)
!!
!! logical, dimension(:,:,:) :: Ice%mask    & ! ocean/land mask for temperature points
!!                                            ! (true for ocean, with or without ice)
!!                              Ice%mask_uv & ! ocean/land mask for momentum points
!!                                            ! (true for ocean, with or without ice)
!!                              Ice%ice_mask  ! optional ice mask (true for ice)
!!
!! real, dimension(:,:,:) :: Ice%part_size  & ! fractional area of each partition of a temperature grid box
!!                           Ice%part_size_uv ! fractional area of each partition of a momentum grid box
!! type (ice_boundary_data_type) :: Ice
!!
!! real, dimension(:) :: Ice%lon_bnd    & ! longitude axis grid box boundaries for temperature points
!!                                        ! in radians (must be monotonic)
!!                       Ice%lat_bnd    & ! latitude axis grid box boundaries for temperature points
!!                                        ! in radians (must be monotonic)
!!                       Ice%lon_bnd_uv & ! longitude axis grid box boundaries for momentum points
!!                                        ! in radians (must be monotonic)
!!                       Ice%lat_bnd_uv   ! latitude axis grid box boundaries for momentum points
!!                                        ! in radians (must be monotonic)
!!
!! logical, dimension(:,:,:) :: Ice%mask    & ! ocean/land mask for temperature points
!!                                            ! (true for ocean, with or without ice)
!!                              Ice%mask_uv & ! ocean/land mask for momentum points
!!                                            ! (true for ocean, with or without ice)
!!                              Ice%ice_mask  ! optional ice mask (true for ice)
!!
!! real, dimension(:,:,:) :: Ice%part_size  & ! fractional area of each partition of a temperature grid box
!!                           Ice%part_size_uv ! fractional area of each partition of a momentum grid box
!!
!! The following fields are located on the ice top grid
!!
!! ~~~~~~~~~~{.f90}
!! real, dimension(:,:,:) :: Ice%t_surf     & ! surface temperature (deg k)
!!                           Ice%albedo     & ! surface albedo (fraction)
!!                           Ice%rough_mom  & ! surface roughness for momentum (m)
!!                           Ice%rough_heat & ! surface roughness for heat/moisture (m)
!!                           Ice%u_surf     & ! zonal (ocean/ice) current at the surface (m/s)
!!                           Ice%v_surf       ! meridional (ocean/ice) current at the surface (m/s)
!! ~~~~~~~~~~
!!
!! The following fields are located on the ice bottom grid
!!
!! ~~~~~~~~~~{.f90}
!! real, dimension(:,:,:) :: Ice%flux_u  & ! zonal wind stress (Pa)
!!                           Ice%flux_v  & ! meridional wind stress (Pa)
!!                           Ice%flux_t  & ! sensible heat flux (w/m2)
!!                           Ice%flux_q  & ! specific humidity flux (kg/m2/s)
!!                           Ice%flux_sw & ! net (down-up) shortwave flux (w/m2)
!!                           Ice%flux_lw & ! net (down-up) longwave flux (w/m2)
!!                           Ice%lprec   & ! mass of liquid precipitation since last time step (Kg/m2)
!!                           Ice%fprec   & ! mass of frozen precipitation since last time step (Kg/m2)
!!                           Ice%runoff    ! mass of runoff water since last time step (Kg/m2)
!! type (ocean_boundary_data_type) :: Ocean
!!
!! real, dimension(:) :: Ocean%Data%lon_bnd     & ! longitude axis grid box boundaries for temperature
!!                                                ! points on the ocean DATA GRID (radians)
!!                       Ocean%Data%lat_bnd     & ! latitude axis grid box boundaries for temperature
!!                                                ! points on the ocean DATA GRID (radians)
!!                       Ocean%Data%lon_bnd_uv  & ! longitude axis grid box boundaries for momentum
!!                                                ! points on the ocean DATA GRID (radians)
!!                       Ocean%Data%lat_bnd_uv  & ! latitude axis grid box boundaries for momentum
!!                                                ! points on the ocean DATA GRID (radians)
!!                       Ocean%Ocean%lon_bnd    & ! longitude axis grid box boundaries for temperature
!!                                                ! points on the ocean MODEL GRID (radians)
!!                       Ocean%Ocean%lat_bnd    & ! latitude axis grid box boundaries for temperature
!!                                                ! points on the ocean MODEL GRID (radians)
!!                       Ocean%Ocean%lon_bnd_uv & ! longitude axis grid box boundaries for momentum
!!                                                ! points on the ocean MODEL GRID (radians)
!!                       Ocean%Ocean%lat_bnd_uv & ! latitude axis grid box boundaries for momentum
!!                                                ! points on the ocean MODEL GRID (radians)
!! ~~~~~~~~~~
!!
!! \note The data values in all longitude and latitude grid box boundary
!!       array must be monotonic.
!!
!! ~~~~~~~~~~{.f90}
!! logical, dimension(:,:) :: Ocean%Data%mask    & ! ocean/land mask for temperature points on the ocean
!!                                                 ! DATA GRID (true for ocean)
!!                            Ocean%Data%mask_uv & ! ocean/land mask for momentum points on the ocean
!!                                                 ! DATA GRID (true for ocean)
!!                            Ocean%Ocean%mask   & ! ocean/land mask for temperature points on the ocean
!!                                                 ! MODEL GRID (true for ocean)
!!                            Ocean%Ocean%mask_uv  ! ocean/land mask for momentum points on the ocean
!!                                                 ! MODEL GRID (true for ocean)
!! real, dimension(:,:) :: Ocean%t_surf_data & ! surface temperature on the ocean DATA GRID (deg k)
!!                         Ocean%t_surf      & ! surface temperature on the ocean MODEL GRID (deg k)
!!                         Ocean%u_surf      & ! zonal ocean current at the surface on the ocean
!!                                             ! MODEL GRID (m/s)
!!                         Ocean%v_surf      & ! meridional ocean current at the surface on the
!!                                             ! ocean MODEL GRID (m/s)
!!                         Ocean%frazil        ! frazil at temperature points on the ocean MODEL GRID
