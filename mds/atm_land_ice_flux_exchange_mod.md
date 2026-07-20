# atm_land_ice_flux_exchange_mod

Module atm_land_ice_flux_exchange_mod handles flux exchange between atmosphere to land and ice.

## variables

| Name | Type | Definition |
|------|------|------------|
| version | character(len=128) | is a program version number that is set automatically during compile time |
| tag | character(len=128) | is a Github tag that is set automatically during compile time |
| xmap_sfc | type(fmsxgridxmap_type), save | is a FmsXgridXmap_type that holds the exchange grid between different components |
| n_xgrid_sfc | integer | is the total number of exchange grid cells |
| mod_name | character(len=4), parameter | is the module name used when registering variable for diag_manager |
| id_drag_mois | integer | is a diag_manager register field id for 'drag coefficient for moisture' |
| id_drag_hea | integer | is a diag_manager register field id for 'drag coefficient for heat' |
| id_drag_mom | integer | is a diag_manager register field id for 'drag coefficient for momentum' |
| id_rough_mois | integer | is a diag_manager register field id for 'surface roughness for moisture' |
| id_rough_hea | integer | is a diag_manager register field id for 'surface roughness for heat' |
| id_rough_mom | integer | is a diag_manager register field id for 'surface roughness for momentum' |
| id_land_mask | integer | is a diag_manager register field id for 'fractional amount of sea ice' |
| id_ice_mask | integer | is a diag_manager register field id for 'fractional amount of land' |
| id_u_star | integer | is a diag_manager register field id for 'friction velocity' |
| id_b_star | integer | is a diag_manager register field id for 'bouyancy scale' |
| id_q_star | integer | is a diag_manager register field id for 'moisture scale' |
| id_u_flux | integer | is a diag_manager register field id for 'zonal wind stress' |
| id_v_flux | integer | is a diag_manager register field id for 'meridional wind stress' |
| id_t_surf | integer | is a diag_manager register field id for 'surface temperature' |
| id_t_ocean | integer | is a diag_manager register field id for 'surface temperature from ocean output' |
| id_t_flux | integer | is a diag_manager register field id for 'sensible heat flux' |
| id_r_flux | integer | is a diag_manager register field id for 'net (down-up) longwave flux' |
| id_q_flux | integer | is a diag_manager register field id for 'evaporation rate' |
| id_slp | integer | is a diag_manager register field id for 'sea level pressure' |
| id_t_atm | integer | is a diag_manager register field id for 'temperature at lowest atmospheric level' |
| id_u_atm | integer | is a diag_manager register field id for 'u wind component at lowest atmospheric level' |
| id_v_atm | integer | is a diag_manager register field id for 'v wind component at lowest atmospheric level' |
| id_wind | integer | is a diag_manager register field id for 'wind speed for flux calculations' |
| id_thv_atm | integer | is a diag_manager register field id for 'surface air virtual potential temperature' |
| id_thv_surf | integer | is a diag_manager register field id for 'surface virtual potential temperature' |
| id_t_ref | integer | is a diag_manager register field id for 'temperature at z_ref_heat' |
| id_rh_ref | integer | is a diag_manager register field id for 'relative humidity at z_ref_heat' |
| id_u_ref | integer | is a diag_manager register field id for 'zonal wind component at z_ref_mom' |
| id_v_ref | integer | is a diag_manager register field id for 'meridional wind component at z_ref_mom' |
| id_wind_ref | integer | is a diag_manager register field id for 'absolute value of wind at z_ref_mom' |
| id_del_h | integer | is a diag_manager register field id for 'ref height interp factor for for heat' |
| id_del_m | integer | is a diag_manager register field id for 'ref height interp factor for momentum' |
| id_del_q | integer | is a diag_manager register field id for 'ref height interp factor for moisture' |
| id_rough_scale | integer | is a diag_manager register field id for 'topographic scaling fractor for momentum drag' |
| id_t_ca | integer | is a diag_manager register field id for 'canopy air temperature' |
| id_z_atm | integer | is a diag_manager register field id for 'height of lowest atmospheric level' |
| id_p_atm | integer | is a diag_manager register field id for 'pressure at lowest atmospheric level' |
| id_gus | integer | is a diag_manager register field id for 'gust scale' |
| id_t_ref_land | integer | is a diag_manager register field id for 'temperature at z_ref_heat over land' |
| id_rh_ref_land | integer | is a diag_manager register field id for 'relative humidity at z_ref_heat over land' |
| id_u_ref_land | integer | is a diag_manager register field id for 'zonal wind component at z_ref_mom over land' |
| id_v_ref_land | integer | is a diag_manager register field id for 'meridional wind component at z_ref_mom over land' |
| id_q_ref | integer | is a diag_manager register field id for 'specific humidity at z_ref_heat' |
| id_q_ref_land | integer | is a diag_manager register field id for 'specific humidity at z_ref_heat over land' |
| id_q_flux_land | integer | is a diag_manager register field id for 'evaporation rate over land' |
| id_rh_ref_cmip | integer | is a diag_manager register field id for 'relative humidity at z_ref_heat' |
| id_husslut_land | integer | is a diag_manager register field id for 'near-surface specific humidity on land use tile' |
| id_taslut_land | integer | is a diag_manager register field id for 'near-surface air temperature a z_ref_heat above displacement height on land-use tile' |
| id_t_flux_land | integer | is a diag_manager register field id for 'sensible heat flux over land' |
| id_co2_atm_dvmr | integer | is a diag_manager register field id for 'co2 dry volume mixing ratio at lowest atmospheric level' |
| id_co2_surf_dvmr | integer | is a diag_manager register field id for 'c02 dry volume mixing ratio at surface' |
| id_co2_bo | integer | is a diag_manager register field id for 'concentration of co2 to be passed to land/photosynthesis' |
| id_co2_flux_pcair_atm | integer | is a diag_manager register field id for 'concentration of co2 to be passed to ocean' |
| id_o2_flux_pcair_atm | integer | is a diag_manager register field id for 'concentration of o2 to be passed to ocean' |
| id_tr_atm | integer, dimension(:), allocatable | is an array of diag_manager register field ids for 'tracers at lowest atmospheric level' |
| id_tr_surf | integer, dimension(:), allocatable | is an array of diag_manager register field ids for 'tracers at surface' |
| id_tr_flux | integer, dimension(:), allocatable | is an array of diag_manager register field ids for 'tracers fluxes' |
| id_tr_mol_flux | integer, dimension(:), allocatable | is an array of diag_manager register field ids for 'flux of co2 concentration in [mol/m2*s]' |
| id_tr_ref | integer, dimension(:), allocatable | is an array of diag_manager register field ids for 'tracers at z_ref_heat' |
| id_tr_ref_land | integer, dimension(:), allocatable | is an array of diag_manager register field ids for 'tracer flux at z_ref_heat over land' |
| id_tr_mol_flux0 | integer, dimension(:), allocatable | is an array of diag_manager register field ids for 'gross flux of tracer concentration over land in [mol/m2*s]' |
| id_tr_flux_land | integer, dimension(:), allocatable | is an array of diag_manager register field ids for 'flux of tracer concentration over land in [kg/m2*s]' |
| id_tr_mol_flux_land | integer, dimension(:), allocatable | is an array of diag_manager register field ids for 'flux of tracer concentration over land in [mol/m2*s]' |
| id_tr_con_atm_land | integer, dimension(:), allocatable | is an array of diag_manager register field ids for 'deposition velocity at lowest atmospheric level over land' Used only when USE_LEGACY_LAND macro is set at compile time |
| id_tr_con_ref_land | integer, dimension(:), allocatable | is an array of diag_manager register field id for 'deposition velocity at reference height over land' |
| id_tr_con_atm | integer, dimension(:), allocatable | is an array of diag_manager register field ids for 'deposition velocity at lowest atmospheric level (atm)'. Used only when USE_LEGACY_LAND macro is set at compile time |
| id_tr_con_ref | integer, dimension(:), allocatable | is an array of diag_manager register field ids for 'deposition velocity at reference height (atm)' |
| id_tas | integer | is a diag_manager register field id for 'near-surface air temperature' (cmip) |
| id_uas | integer | is a diag_manager register field id for 'eastward near-surface wind' (cmip) |
| id_vas | integer | is a diag_manager register field id for 'northward near-surface wind' (cmip) |
| id_ts | integer | is a diag_manager register field id for 'surface temperature' (cmip) |
| id_psl | integer | is a diag_manager register field id for 'air pressure at sea level' (cmip) |
| id_sfcwind | integer | is a diag_manager register field id for 'near-surface wind speed' (cmip) |
| id_tauu | integer | is a diag_manager register field id for 'surface downward eastward wind stress' (cmip) |
| id_tauv | integer | is a diag_manager register field id for 'surface downward northward wind stress' (cmip) |
| id_hurs | integer | is a diag_manager register field id for 'near-surface relative humidty' (cmip) |
| id_huss | integer | is a diag_manager register field id for 'near-surface specific humidity' (cmip) |
| id_evspsbl | integer | is a diag_manager register field id for 'water evaporation flux' (cmip) |
| id_hfls | integer | is a diag_manager register field id for 'surface upward latent heat flux' (cmip) |
| id_hfss | integer | is a diag_manager register field id for 'surface upward sensible heat flux' (cmip) |
| id_rhs | integer | is a diag_manager register field id for 'near-surface relative humidty' (cmip) |
| id_sftlf | integer | is a diag_manager register field id for 'fraction of the grid cell occupied by land' (cmip) |
| id_tos | integer | is a diag_manager register field id for 'sea surface temperature' (cmip) |
| id_sic | integer | is a diag_manager register field id for 'sea ice area fraction' (cmip) |
| id_tslsi | integer | is a diag_manager register field id for 'surface temperature on land or sea ice' (cmip) |
| id_height2m | integer | is a diag_manager register field id for 'near surface height' (cmip) |
| id_height10m | integer | is a diag_manager register field id for 'near surface height' (cmip) |
| id_evspsbl_g | integer | is a diag_manager register field id for 'global integral of water evaporation flux' |
| id_ts_g | integer | is a diag_manager register field id for 'global integral of surface temperature' |
| id_tas_g | integer | is a diag_manager register field id for 'global integral of near-surface air temperature' |
| id_tasl_g | integer | is a diag_manager register field id for 'global integral of near-surface air temperature on land' |
| id_hfss_g | integer | is a diag_manager register field id for 'global integral of surface upward sensible heat flux' |
| id_hfls_g | integer | is a diag_manager register field id for 'global integral of surface upward latent heat flux' |
| id_rls_g | integer | is a diag_manager register field id for 'global integral of near-surface relative humidty' |
| first_static | logical | is a flag where if true, land_mask, sftlf, height2m, and height10m are saved once per file at first call to sf_boundary_layer |
| do_init | logical | is a flag where if true, if atm_land_ice_flux_exchnge_init has been called |
| remap_method | integer | is a flag value to indicate first or second order conservative remapping onto exchange grid |
| d622 | real, parameter | is the value rdgas/rvgas |
| d378 | real, parameter | is the value 1.0-d622 |
| d608 | real, parameter | is the value d378/d622 |
| tfreeze | real, parameter | is the freezing point of water at 1 atm [K] |
| frac_precip | real, dimension(:,:), allocatable | is an array holding the scale values for precipitation |
| z_ref_heat | real | is the reference height (meters) for temperature and relative humidity diagnostics (t_ref, rh_ref, del_h, del_q) |
| z_ref_mom | real | is the reference height (meters) for momentum diagnostics (u_ref, v_ref, del_m) |
| do_forecast | logical | is a flag |
| nblocks | integer | is the OpenMP number of thread. Do loops on the exchange grid are parallelized into noblocks |
| partition_fprec_from_lprec | logical | is a flag where if true, liquid precip is converted to snow when t_ref < tfreeze. Used for atm override experiments where liquid and frozen precip are combined |
| scale_precip_2d | logical | is a flag where if true, scale mass of liqud preciptation |
| my_nblocks | integer | is the number of blocks in OpenMP parallelization, defaults to 1 |
| block_start | integer, dimension(:), allocatable | is the starting do loop indices for OpenMP thread |
| block_end | integer, dimension(:), allocatable | is the ending do loop indices for OpenMP thread |
| ex_t_surf | real, dimension(:), allocatable | is the surface temperature for radiation calc on exchange grid [K]. Note, T canopy is only differet from t_surf over vegetated land |
| ex_t_surf_miz | real, dimension(:), allocatable | is the surface temperature in the marginal ice zone (MIZ) on the exchange grid [K]; Used in forecast mode (do_forecast=.true.) when #ifdef AM3_physics is defined |
| ex_t_ca | real, dimension(:), allocatable | is the near-surface (canopy) air temperature on exchange grid [K] |
| ex_p_surf | real, dimension(:), allocatable | is the surface pressure on exchange grid on the exchange grid |
| ex_slp | real, dimension(:), allocatable | is the surface pressure on exchange grid |
| ex_flux_ | real, dimension(:), allocatable | is the sens heat flux on the exchange grid |
| ex_flux_lw | real, dimension(:), allocatable | is the longwave radiation flux on the exchange grid |
| ex_dhdt_surf | real, dimension(:), allocatable | is d(sens.heat.flux)/d(T canopy) on the exchange grid |
| ex_dedt_surf | real, dimension(:), allocatable | is d(water.vap.flux)/d(T canopy) on the exchange grid |
| ex_dqsatdt_surf | real, dimension(:), allocatable | is the d(water.vap.flux)/d(q canopy) on the exchange grid |
| ex_e_q_n | real, dimension(:), allocatable | is dt/mass * dedet_surf * gamma on the exchange grid |
| ex_drdt_surf | real, dimension(:), allocatable | is d(LW flux)/d(T surf) on the exchange grid |
| ex_dhdt_atm | real, dimension(:), allocatable | is d(sens.heat.flux)/d(T atm) on the exchange grid |
| ex_flux_u | real, dimension(:), allocatable | is the u stress on atmosphere on the exchange grid |
| ex_flux_v | real, dimension(:), allocatable | is the v stress on atmosphere on the exchange grid |
| ex_dtaudu_atm | real, dimension(:), allocatable | is d(stress)/d(u) on the exchange grid |
| ex_dtaudv_atm | real, dimension(:), allocatable | is d(stress)/d(v) on the exchange grid |
| ex_seawater | real, dimension(:), allocatable | is a mask array of seaice fractions on the exchange grid. Takes value of 1 when there is any open water in the OCN grid cell. Takes value of 0 when there is no open water in the OCN grid cell (i.e, totally covered with ice or land). Note, ex_seawater should not be mistaken with ex_avail where ex_avail is 1 for all OCN grid cells regardless of sea-ice coverage. |
| ex_albedo_vis_dir_fix | real, dimension(:), allocatable | is the fixed (overridden) direct visible-band surface albedo on the exchange grid [dimensionless]; used to apply data_override corrections to the albedo before flux calculations |
| ex_albedo_nir_dir_fix | real, dimension(:), allocatable | is the fixed (overridden) direct near-infrared surface albedo on the exchange grid [dimensionless] |
| ex_albedo_vis_dif_fix | real, dimension(:), allocatable | is the fixed (overridden) diffuse visible-band surface albedo on the exchange grid [dimensionless] |
| ex_albedo_nir_dif_fix | real, dimension(:), allocatable | is the fixed (overridden) diffuse near-infrared surface albedo on the exchange grid [dimensionless] |
| ex_drag_q | real, dimension(:), allocatable | is the q drag coefficient on the exchange grid |
| ex_cd_ | real, dimension(:), allocatable | is the drag coefficient for heat on the exchange grid |
| ex_cd_m | real, dimension(:), allocatable | is the drag coefficient for momentum on the exchange grid |
| ex_b_star | real, dimension(:), allocatable | is the boyuancy scale on the exchange grid |
| ex_u_star | real, dimension(:), allocatable | is the friction velocity on exchange grid |
| ex_wind | real, dimension(:), allocatable | is the wind speed on exchange grid |
| ex_z_atm | real, dimension(:), allocatable | is the height of lowest atmospheric level on exchange grid |
| ex_con_atm | real, dimension(:), allocatable | is the deposition velocity at lowest atmospheric level on the exchange grid |
| ex_dhdt_surf_forland | real, dimension(:), allocatable |  |
| ex_dedt_surf_forland | real, dimension(:), allocatable |  |
| ex_dedq_surf_forland | real, dimension(:), allocatable |  |
| ex_tr_surf | real, dimension(:,:), allocatable | is the surface temperature for radiation calc on exchange grid [K] |
| ex_flux_tr | real, dimension(:,:), allocatable | is the tracer fluxes on the exchange grid |
| ex_dfdtr_surf | real, dimension(:,:), allocatable | is the d(tracer flux)/d(surf tracer) on the exchange grid |
| ex_dfdtr_atm | real, dimension(:,:), allocatable | is the d(tracer flux)/d(atm tracer) on the exchange grid |
| ex_e_tr_n | real, dimension(:,:), allocatable | is the coefficient in implicit scheme on the exchange grid |
| ex_f_tr_delt_n | real, dimension(:,:), allocatable | is the coefficient in implicit scheme on the exchange grid |
| ex_tr_con_ref | real, dimension(:,:), allocatable | is the deposition velocity at reference height on the exchange grid |
| ex_tr_con_atm | real, dimension(:,:), allocatable | is the deposition velocity at atmospheric height on the exchange grid |
| ex_avail | logical, dimension(:), allocatable | if a mask array where if true, the exchange grid cell is over ocean and/or seaice |
| ex_land | logical, dimension(:), allocatable | is a mask array where if true, the exchange grid cell is over land |
| ex_e_t_n | real, dimension(:), allocatable | is the implicit coupling coefficient e_t^n for sensible heat on the exchange grid; used in the semi-implicit surface flux scheme to couple atmospheric temperature to surface temperature |
| ex_f_t_delt_n | real, dimension(:), allocatable | is the implicit coupling coefficient f_t^{delta n} for sensible heat on the exchange grid; represents the accumulated forcing term in the semi-implicit time integration of surface heat flux |
| n_atm_tr | integer | is the number of prognostic tracers in the atmos model |
| n_atm_tr_to | integer | is the number of prognostic tracers in the atmos model |
| n_lnd_tr | integer | is the number of prognostic tracers in the land model |
| n_lnd_tr_to | integer | is the number of prognostic tracers in the land model |
| n_exch_tr | integer | is the number of tracers exchanged between models |
| n_gex_atm2lnd | integer | is the number of gex fields exchanged between land and atmosphere |
| n_gex_lnd2atm | integer | is the number of gex fields exchanged between atmosphere and land |
| tr_table | type(tracer_ind_type), dimension(:), allocatable | is the table of tracers passed through flux exchange |
| tr_table_map | type(tracer_exch_ind_type), dimension(:), allocatable | is a map to map atm tracers to exchange, ice and land variables |
| isphum | integer | is the specific humidity index. Initialized as NO_TRACER |
| ico2 | integer | is the co2 tracer index. Initialized as NO_TRACER |
| inh3 | integer | is the nh3 tracer index. Initialized as NO_TRACER |
| ex_gas_fields_atm | type(fmscoupler1dbc_type), pointer | contains atmospheric gas fields used for atm-ocean flux exchange |
| ex_gas_fields_ice | type(fmscoupler1dbc_type), pointer | contains ice-top and ocean_surface gas fields |
| ex_gas_fluxes | type(fmscoupler1dbc_type), pointer | contains gas fluxes between atmosphere and ocean |
| ni_atm | integer | is the number of x gridpoints in the atm compute domain |
| nj_atm | integer | is the number of y gridpoints in the atm compute domain |
| regrid | integer, parameter | is the boundary_typextype value when grids are physically different and data between model components needs to be exchanged via the exchange grid |
| redist | integer, parameter | is the boundary_typextype value when grids are physically same, but differ in domain decomposition. |
| direct | integer, parameter | is the boundary_typextype value when grids and domains are identical and data can be copied directly beteween components |
| cplclock | integer | is a FMS clock id for profiling general processes |
| sfcclock | integer | is a FMS clock id for profiling sfc_boundary_layer |
| fluxatmdnclock | integer | is a FMS clock id for profiling flux down from atmosphere |
| regenclock | integer | is a FMS clock for profiling exchange grid generation |
| fluxatmupclock | integer | is a FMS clock for profiling flux up to atmosphere |
| x1_grid_atm | integer | is the exchange grid index for xgrid_stock_move. Set to value of 1 |
| x1_grid_ice | integer | is the exchange grid index for xgrid_stock_move. Set to value of 2 |
| x1_grid_lnd | integer | is the exchange grid index for xgrid_stock_move. Set to value of 3 |
| dt_atm | real | is the atmospheric timestep [s] |
| dt_cpl | real | is the coupled timestep [s] |
| nxc_ice | integer | is the number of x gridpoints in ice compute domain |
| nyc_ice | integer | is the number of y gridpoints in ice compute domain |
| nk_ice | integer | is the number of vertical levels in ice |
| nxc_lnd | integer | is the number of x gridpoints in land compute domain |
| nyc_lnd | integer | is the number of y gridpoints in land compute domain |


## atm_land_ice_flux_exchange_init
### intro
Module atm_land_ice_flux_exchange_mod handles flux exchange between atmosphere to land and ice.  atm_land_ice_flux_exchange_init is a subroutine in atm_land_ice_flux_exchange_mod.
### description
  Subroutine atm_land_ice_flux_exchange_init initializes atm_land_ice_flux_exchange_mod by allocating and seting default values for module level variable; and calling initialization routines in FMS modules. This subroutine must be called before calling any other public procedures in this module.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| time | intent(inout) | atm_land_ice_flux_exchange_init | is the current model time |
| atm | intent(inout) | atm_land_ice_flux_exchange_init | is a derived data type holding atmosphere boundary data |
| land | intent(inout) | atm_land_ice_flux_exchange_init | is a derived data type holding land boundary data |
| ice | intent(inout) | atm_land_ice_flux_exchange_init | is a derived data type holding ice boundary data |

### flowchart
atm_land_ice_flux_exchange_init does the following:  
Step 1: INITIALIZE MODULE-LEVEL VARIABLES.
Step 2: GET FILE UNIT FOR STDOUT AND STDLOG FOR INTERNAL LOGGING PURPOSES
Step 3: FROM THE TRACER TABLE, GET THE TOTAL NUMBER TRACERS, TOTAL NUMBER OF SPECIFIC HUMIDITY TRACER, AND TOTAL NUMBER PROGNOSTIC TRACERS IN ATMOSPHERE AND LAND MODELS
Step 4: CONSTRUCT THE TRACER TABLE (TR_TABLE): FOR EACH TRACER, RECORD THE TRACER_INDEX IN THE ATM MODEL, ICE MO! DEL, AND LAND MODEL. SKIP ALL ATMOS TRACERS THAT DO NOT HAVE CORRESPONDING SURFACE TRA! CERS.
Step 5: GET THE TOTAL NUMBER OF GENERIC EXCHANGE FIELDS BETWEEN ATMOSPHER! E AND LAND
Step 6: GET THE TRACER INDEX OF SPECIFIC HUMIDITY.
Step 7: INITIALIZE FRAC_PRECIP IF SCALE_PRECIP_2D IS TRUE
Step 8: SET UP THE EXCHANGE GRID AND SET X1_GRID_ATM = 1, X1_GRID_ICE = 2! , AND X1_GRID_LAND = 3. SETS XMAP_SFC(1)GRIDS FOR ATM, XMAP_SFC(2)GRIDS FOR ICE, XMAP_S! FC(3)GRIDS FOR LAND
Step 9: INITIALIZE SURFACE_FLUX MODULE
Step 10: INITIALLIZE FMS DIAG_INTEGRAL FIELDS FOR EVAP, T_SURF, T_REF GLOB! AL INTEGRAL QUANTITIES.call diag_integral_field_init ('prec', 'f6.3')
Step 11: REGISTER FMS DIAGNOSTIC FIELDS IN DIAG_MANAGER.
Step 12: GET THE SIZE OF THE ATM COMPUTE DOMAIN.
Step 13: ALLOCATE ATMOS_ICE_BOUNDARY AND SET FIELDS EQUAL TO ZERO.
Step 14: ALLOCATE LAND_ICE_ATMOS_BOUNDARY AND SET FIELDS EQUAL TO ZERO EXC! EPT FOR T_OCEAN WHICH IS SET TO 200 K, T_REF TO 273 K, AND ROUGHNESS LENG! THS TO 0.01 m.
Step 15: COPY EX_GAS_FIELDS_ATM TO ATMFIELDS.
Step 16: GET THE SIZE OF ICE COMPUTE DOMAIN.
Step 17: GET THE SIZE OF LAND COMPUTE DOMAIN.
Step 18: DECLARE CLOCKS FOR PROFILING.
Step 19: SET DO_INIT = .FALSE. IN ORDER TO AVOID RE-INITIALIZATION THE MOD! ULE IF THIS SUBROUTINE IS CALLED AGAIN.


## sfc_boundary_layer
### intro
Module atm_land_ice_flux_exchange_mod handles flux exchange between atmosphere to land and ice.  sfc_boundary_layer is a subroutine in atm_land_ice_flux_exchange_mod.
### description
  Subroutine sfc_boundary_layer computes and exchanges the following ! fluxes:t_surf_atm, albedo_atm, rough_mom_atm, land_frac_atm, dtaudu_atm, ! dtaudv_atm, flux_u_atm, flux_v_atm u_star_atm, b_star_atmNote, u_star and b_star are defined so that u_star**2 is the magnit! ude of surface stress divided by density of air at the surface, and u_star*b_star is the ! buoyancy flux at the surface.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| time | intent(inout) | sfc_boundary_layer | is the current model time |
| atm | intent(inout) | sfc_boundary_layer | is a derived type holding atmosphere boundary data |
| land | intent(inout) | sfc_boundary_layer | is a derived type holding land boundary data |
| ice | intent(inout) | sfc_boundary_layer | is a derived type holding ice boundary data |

### flowchart
sfc_boundary_layer does the following:  
Step 1: CHECK MODULE INITIALIZATION.
Step 2: INITIALIZE CLOCKS FOR PROFILING.
Step 3: ALLOCATE ARRAY FOR EXCHANGE FIELDS. THE ARRAYS ARE DEALLOCATED I! N FLUX_UP_TO_ATMOS.
Step 4: ALLOCATE EX_GAS_FIELDS_ICE ARRAYS FOR OCEAN_ICE_BOUNDARY EXCHANGE ! FIELDS.
Step 5: ALLOCATE EX_GAS_FIELDS_ATM ARRAYS FOR ATMOSPHERE EXCHANGE FIELDS.
Step 6: ALLOCATE EX_GAS_FLUXES FOR ADDITIONAL EXCHANGE FIELDS.
Step 7: ON THE EXCHANGE GRID, SET INITIAL VALUES FOR ALBEDO, DRAG COEFFICIENTS, AND OPEN WATER MASK.
Step 8: OVERRIDE ATM FIELDS IF FIELD EXISTS IN THE DATA_OVERRIDE DATA_TABLE.
Step 9: CONVERT CO2 TRACER UNITS TO WET_MMR UNITS.
Step 10: OVERRIDE CO2 VALUES IF THE FIELD EXISTS IN DATA_OVERRIDE DATA_TABLE AND SEND DATA TO THE DIAG_MANAGER BUFFER.
Step 11: OVERRIDE ICE AND LAND FIELD IF THE FIELD EXISTS IN DATA_OVERRIDE DATA_TABLE.
Step 12: MAP ATM FIELDS ONTO THE EXCHANGE GRID.
Step 13: INITIALIZE EX_TR_SURF TO BE THE AMOUNT OF TRACERS AT THE BOTTOM-M! OST ATMOSPHERE LAYER.
Step 14: MAP ICE FIELDS ONTO THE EXCHANGE GRID.
Step 15: ON THE EXCHANGE GRID, GENERATE DYNAMIC WET MASK ARRAY WITH VALUE ! OF 1.O FOR OPEN WATER.
Step 16: ON THE EXCHANGE GRID, INITIALIZE CANOPY TEMPERATURE EQUAL TO SURFACE TEMPERATURE.
Step 17: MAP LAND EXCHANGE FIELDS ONTO THE EXCHANGE GRID.
Step 18: ON THE EXCHANGE GRID, COMPUTE EXPLICIT FLUXES AND TENDENCIES.
Step 19: CALL MONIN_OBUKHOV_MO_PROFILE IN FMS. ON THE EXCHANGE GRID, COMPUTE ZONAL AND MERIDIONAL WINDS AT THE BOUNDARY LAYER AND AT REFERENCE HEIGHTS.
Step 20: ON THE EXCHANGE GRID, CALCULATE ATMOSPHERIC CONDUCTANCE.
Step 21: ON THE EXCHANGE GRID, COMPUTE DERIVATIVES OF TRACER FLUXES.
Step 22: ON THE EXCHANGE GRID, COMPUTE EXPLICIT FLUXES BETWEEN ATM AND OCEAN.
Step 23: OVERRIDE LAND AND ICE TRACER FLUXES IF FIELD EXISTS IN DATA_OVERRIDE DATA_TABLE.
Step 24: ON THE EXCHANGE GRID, COMPUTE T_SURF**4. NOTE, TO COMPUTE FLUXES, T_SURF**4 (T_SURF TO THE FOURTH POWER) IS SENT TO THE EXCHANGE GRID AND NOT T_SURF DUE TO NONLINEARITY IN THE STEFAN-BOLTZMANN LAW WHERE LOGWAVE_FLUX = STEFAN_BOLTZMANN_CONSTANT * T**4. ON THE EXCHANGE GRID, AS QUANTITIES ARE REMAPPED, FIELDS ARE AREA-WEIGHTED (AVERAGED) SUCH THAT OUTPUT_TEMPERATURE = SUM(INPUT_TEMPERATURE * (XGRID_AREA)/(INPUT_GRID_AREA)) WHERE THE SUM IS OVER ALL XGRID CELLS THAT OVERLAP WITH THE OUTPUT CELL. BECAUSE OF THIS WEIGHTING, THE COMPUTED FLUX WOULD DIFFER FROM USING <T**4> VS <T>**4.
Step 25: MAP FIELDS FROM THE EXCHANGE GRID TO THE ATM GRID.
Step 26: ON THE ATM GRID, COMPUTE T**0.25.
Step 27: DATA OVERRIDE ATMOSPHERIC QUANTITIES IF FOUND IN DATA_TABLE.
Step 28: ON THE EXCHANGE GRID, INITIALIZE ARRAYS FOR FIXING THE ALBEDO.
Step 29: MAP ALBEDO FIELDS TO THE EXCHANGE GRID
Step 30: ON THE EXCHANGE GRID, COMPUTE THE ALBEDO FIXING FACTORS.
Step 31: SAVE STATIC FIELDS.
Step 32: MAP ATMOSPHERIC DATA FROM THE EXCHANGE TO ATM GRID AND SEND DATA ! TO DIAG_MANAGER BUFFER.
Step 33: COMPUTE DIAGNOSTIC FIELDS AT REFERENCE LEVELS WITH FMS_MONIN_OBUH! KOV_MO_PROFILE.
Step 34: MAP LAND AND ADDITIONAL ATMOSPHERIC FIELDS FROM THE EXCHANGE GRID ! TO THE COMPONENT GRID AND SEND DATA TO DIAG_MANAGER BUFFER FOR DIAGNOSTIC OUTPUT NOTE, DATA WILL ONLY BE OUTPUTTED IF VARIABLE SPECIFICATION IS FO! UND IN THE DIAG_TABLE
Step 35: END CLOCKS FOR PROFILING


## flux_down_from_atmos
### intro
Module atm_land_ice_flux_exchange_mod handles flux exchange between atmosphere to land and ice.  flux_down_from_atmos is a subroutine in atm_land_ice_flux_exchange_mod.
### description
  Subroutine flux_down_from_atmos corrects for the implicit treatment ! of atmospheric diffisuve fluxes in flux exchange from atm to land and ice.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| time | intent(inout) | flux_down_from_atmos | is the current model time |
| atm | intent(inout) | flux_down_from_atmos | is a derived type holding atmosphere boundary data |
| land | intent(inout) | flux_down_from_atmos | is a derived type holding land boundary data |
| ice | intent(inout) | flux_down_from_atmos | is a derived type holding ice boundary data |
| atmos_boundary | intent(inout) | flux_down_from_atmos | is a derived type holding properties and fluxes passed from excha! nge grid to atmosphere land and ice |

### flowchart
flux_down_from_atmos does the following:  
Step 1: START CLOCKS FOR PROFILING
Step 2: INITIALIZE REUSABLE FLAG. DATA_OVERRIDE WILL RETURN OV=.TRUE. IF ! DATA WAS OVERWRITTEN
Step 3: OVERRIDE ATM SHORTWAVE/LONGWAVE, DIRECT/DOWNWARD DIFFUSIVE FLUXES NOTE, DATA_OVERRIDE WILL ONLY OVERWRITE IF THE FIELD IS SPECIFIED ! IN THE DATA_TABLE.
Step 4: SCALE LIQUID PRECIPITATION BY FRAC_PRECIP IF SCALE_PRECIP_2D IS T! RUE SCALE_PRECIP_2D IS SET DURING MODULE INITIALIZATION CALL TO ATM_L! AND_ICE_FLUX_EXCHANGE_INIT FRAC_PRECIP VALUES ARE SET WITH DATA_OVERRIDE
Step 5: PARTITION PRECIPTATION TO LIQUID PRECIPITATION AND FROZEN PRECIPI! TATION IF PARTITION_FPREC_FROM_LPREC = .TRUE. PARTIION_FPREC_FROM_LPREC IS SET AS PART OF MODULE INITIALIZATION ! CALL IN ATM_LAND_ICE_FLUX_EXCHANGE
Step 6: OVERRIDE ATM FPREC, COZEN, AND SURF_DIFF FIELDS. NOTE, DATA_OVERRIDE WILL ONLY OVERWRITE ARRAY IF THE FIELD IS SPE! CIFIED IN THE DATA_TABLE
Step 7: MAP ATMOSPHERE QUANTITIES ONTO THE EXCHANGE GRID
Step 8: ON THE EXCHANGE GRID, UPDATE U AND V STRESS
Step 9: ON THE EXCHANGE GRID, TAKE INTO ACCOUNT FOR ALBEDO VARIATION IN SHORTWAVE RADIATION FLUX OF VISIBLE LIGHT
Step 10: ON THE EXCHANGE GRID, ADJUST FLUXES FOR IMPLICIT DEPENDENCE
Step 11: MAP FLUXES FROM THE EXCHANGE GRID TO THE LAND GRID AND OVERRIDE FIELDS WITH DATA_OVERRIDE WHERE DATA WILL BE OVERWRITTEN IF THE FIELD IS SPECIFIED IN DATA_TABLE
Step 12: OVERRIDE LAND FLUXES. NOTE, DATA_OVERRIDE WILL ONLY OVERWRITE ARRAY IF THE FIELD IS SPECIFIED IN THE DATA_TABLE
Step 13: MAP ICE FIELDS FROM THE EXCHANGE GRID TO THE ICE GRID
Step 14: OVERRIDE ICE FIELDS. NOTE, DATA_OVERRIDE WILL ONLY OVERWRITE ARRAY IF THE FIELD IS SPECIFIED IN THE DATA_TABLE
Step 15: COMPUTE STOCK EXCHANGES BETWEEN COMPONENTS
Step 16: SEND U_FLUX AND V_FLUX TO THE DIAG_MANAGER BUFFER NOTE, DATA WILL ONLY BE OUTPUTTED IF VARIABLE SPECIFICATION IS FO! UND IN DIAG_TABLE.YAML
Step 17: END CLOCK FOR PROFILING


## generate_sfc_xgrid
### intro
Module atm_land_ice_flux_exchange_mod handles flux exchange between atmosphere to land and ice.  generate_sfc_xgrid is a subroutine in atm_land_ice_flux_exchange_mod.
### description
  Subroutine generate_sfc_xgrid updates fractional areas on the surfa! ce exchange grid and recompute the number of active exchange-grid cells. The fractional ! area measures the portion of each exchange-grid cell that corresponds to land or ice. fms_xg! rid_set_frac_area is called for both the OCN (ice) and LND (land) grids to reflect the current ! sea-ice concentration and land-tile coverage.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| land | intent(inout) | generate_sfc_xgrid | is a derived data type to specify land boundary data |
| ice | intent(inout) | generate_sfc_xgrid | is a derived data type to specify ice boundary da |

### flowchart
generate_sfc_xgrid does the following:  
Step 1: INITIALIZE CLOCK FOR PROFILING
Step 2: GET ICE COMPUTE DOMAIN INDICES
Step 3: UPDATE FRACTIONAL AREAS OF THE EXCHANGE GRID THAT ARE ICE AND LAN! D
Step 4: UPDATE THE NUMBER OF EXCHANGE GRID CELLS SAVED IN THE MODULE
Step 5: END CLOCK FOR PROFILING


## flux_up_to_atmos
### intro
Module atm_land_ice_flux_exchange_mod handles flux exchange between atmosphere to land and ice.  flux_up_to_atmos is a subroutine in atm_land_ice_flux_exchange_mod.
### description
  Subroutine flux_up_to_atmos corrects the fluxes to take into accoun the new surface temperatures in land and ice models.The following elements of the land_ice_atmos_boundary_type are computed: dt_t = temperature change at the lowest atmospheric level [K] dt_q = specific humidity change at the lowest atmospheric level [kg/kg].
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| time | intent(inout) | flux_up_to_atmos | Current time is the current model time |
| land | intent(inout) | flux_up_to_atmos | is a derived type holding land boundary data |
| ice | intent(inout) | flux_up_to_atmos | is a derived type holding ice boundary data |

### flowchart
flux_up_to_atmos does the following:  
Step 1: START CLOCK FOR PROFILING
Step 2: OVERRIDE ICET_SURF, LANDT_CA, LANDT_SURF AND LAND SURFACE TRACERS NOTE, DATA_OVERRIDE WILL ONLY OVERWRITE DATA IF THE FIELD IS SPECIFIED IN THE DATA_TABLE
Step 3: INITIALIZE EX_T_SURF_NEW = 200.0
Step 4: MAP ICET_SURF, LANDT_CA AND LANDT_SURF ONTO THE EXCHANGE GRID.
Step 5: ON THE EXCHANGE GRID, COMPUTE CHANGES IN SURFACE TEMPERATURE AND RADIATIVE TEMPERATURE.
Step 6: ON THE EXCHANGE GRID, UPDATE FLUXES AND ATMOSPHERIC INCREMENTS FOR IMPLICIT DEPENDENCE ON SURFACE TEMPERATURE.
Step 7: ON THE EXCHANGE GRID, UPDATE TRACER TENDENCIES IN THE ATMOSPHERE
Step 8: MAP DT_T, SHFLX, and LHFLX FIELDS IN LAND_ICE_ATMOS_BOUNDARY FROM THE EXCHANGE GRID TO THE ATMOSPERE GRID.
Step 9: MAP DATA FROM THE EXCHANGE GRID TO OCN/ATM/LND GRID AND SEND DATA TO THE DIAG_MANAGER BUFFER. NOTE, DATA WILL ONLY BE OUTPUTTED IF VARIABLE SPECIFICATION IS FOUND IN DIAG_TABLE.YAML.
Step 10: COMPUTE STOCK EXCHANGE BETWEEN MODEL COMPONENTS.
Step 11: END CLOCK FOR PROFILING


## flux_ex_arrays_dealloc
### intro
Module atm_land_ice_flux_exchange_mod handles flux exchange between atmosphere to land and ice.  flux_ex_arrays_dealloc is a subroutine in atm_land_ice_flux_exchange_mod.
### description
  Subroutine flux_ex_arrays_dealloc deallocates the model-level excha! nge grid related arrays that were allocated in sfc_boundary_layer.
### arguments
None
### flowchart
flux_ex_arrays_dealloc does the following:  


## flux_atmos_to_ocean
### intro
Module atm_land_ice_flux_exchange_mod handles flux exchange between atmosphere to land and ice.  flux_atmos_to_ocean is a subroutine in atm_land_ice_flux_exchange_mod.
### description
  Subroutine flux_atmos_to_ocean computes deposition gas fluxes betwe! en atmosphere and ocean This subroutine is called only if the do_flux namelist variable is ! set to .True.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| time | intent(inout) | flux_atmos_to_ocean | is the current time |
| atm | intent(inout) | flux_atmos_to_ocean | is a derived data type holding atmosphere boundary data |
| ice_boundary | intent(inout) | flux_atmos_to_ocean | is a derived data type holding properties and fluxes passed from ! atmosphere to ice |
| ice | intent(inout) | flux_atmos_to_ocean | is a derived type holding ice boundary tdata |

### flowchart
flux_atmos_to_ocean does the following:  
Step 1: MAP ATMOSPHERE FIELDS TO THE EXCHANGE MAP FOR FLUX EXCHANGE WITH ! OCEAN.
Step 2: ON THE EXCHANGE GRID, CALCULATE OCEAN EXPLICIT FLUX BY CALLING AT! MOS_OCEAN_DEP_FLUXES_CALC.
Step 3: MAP AIR_SEA_DEPOSITION FLUX FROM THE EXCHANGE GRID TO THE ICE GRI! D FOLLOWED BY CALL DATA_OVERRIDE WHERE DATA WILL BE OVERWRITTEN IF THE FLUX FIE! LDS ARE SPECIFIED IN THE DATA_TABLE. THEN SEND_DATA TO THE DIAG_MANA! GER BUFFER.
Step 4: UPDATE ICE FIELDS THAT ARE LABELED AS AIR_SEA_DEOOSITION FLUXES B! Y CALLING UPDATE_ICE_ATM_DEPOSITION_FLUX.


## put_logical_to_real_sg
### intro
Module atm_land_ice_flux_exchange_mod handles flux exchange between atmosphere to land and ice.  put_logical_to_real_sg is a subroutine in atm_land_ice_flux_exchange_mod.
### description
  Subroutine put_logical_to_real_sg maps 2D logical mask arrays to re! al arrays where .true. is equal to 1.0 and .false. is equal to 0.0. The real ! array is then mapped onto the exchange grid. This subroutine is used internally to con! vert Landmask on structured grid for example, when #ifndef USE_LEGACY_LAND is f! alse.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| mask | intent(inout) | put_logical_to_real_sg | is the land/ice mask |
| id | intent(inout) | put_logical_to_real_sg | is the component id |
| ex_mask | intent(inout) | put_logical_to_real_sg | is the mapped mask on exchange grid |
| xmap | intent(inout) | put_logical_to_real_sg | is the xmap |

### flowchart
put_logical_to_real_sg does the following:  
Step 1: MAP (LOGICAL) MASK TO RMASK
Step 2: MAP RMASK TO THE EXCHANGE GRID


## put_logical_to_real_ug
### intro
Module atm_land_ice_flux_exchange_mod handles flux exchange between atmosphere to land and ice.  put_logical_to_real_ug is a subroutine in atm_land_ice_flux_exchange_mod.
### description
  Subroutine put_logical_to_real_ug maps 2D logical mask arrays to re! al arrays where .true. is equal to 1.0 and .false. is equal to 0.0. The real ! array is then mapped onto the exchange grid. This subroutine is used internally to con! vert Landmask on unstructured grid (when #ifndef USE_LEGACY_LAND is true).
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| mask | intent(inout) | put_logical_to_real_ug | is the mask on component grid |
| id | intent(inout) | put_logical_to_real_ug | is the component id |
| ex_mask | intent(inout) | put_logical_to_real_ug | is the mapped mask on exchange grid |
| xmap | intent(inout) | put_logical_to_real_ug | is the xmap |

### flowchart
put_logical_to_real_ug does the following:  
Step 1: MAP (LOGICAL) MASK TO RMASK
Step 2: MAP RMASK TO THE EXCHANGE GRID


## diag_field_init
### intro
Module atm_land_ice_flux_exchange_mod handles flux exchange between atmosphere to land and ice.  diag_field_init is a subroutine in atm_land_ice_flux_exchange_mod.
### description
  Subroutine diag_field_init registers the diagnostic fields in this ! module to the diag_manager. Note, diagnostic fields must be registered in diag_manager and all ! diagnostics fields must be specified in diag_table.yaml in order for the data to be ou! tputted to a NetCDF file at the end of the model run. This subroutine is c! alled during module initialization in subroutine atm_land_ice_flux_exchan! ge_init.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| time | intent(inout) | diag_field_init | is the curent model time |
| atmos_axes | intent(inout) | diag_field_init | is the array size for atmospheric diagnostic fields |
| land_axes | intent(inout) | diag_field_init | is the array size for land diagnostic fields |
| land_pe | intent(inout) | diag_field_init | is the land pe number |

### flowchart
diag_field_init does the following:  
Step 1: CONVERT DIAGNOSTIC LABELS FROM INTEGERS TO STRINGS
Step 2: CALL FMS_DIAG_REGISTER_DIAG_FIELD


## divide_by_area
### intro
Module atm_land_ice_flux_exchange_mod handles flux exchange between atmosphere to land and ice.  divide_by_area is a subroutine in atm_land_ice_flux_exchange_mod.
### description
  Subroutine divide_by_area divides data on a grid by the grid cell area only for cells with non-zero area. This subroutine iscurrently not used. Note, a similar subroutine also exists in ice_ocean_flux_exchange_mod.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| data | intent(inout) | divide_by_area | is the data to be divided |
| area | intent(inout) | divide_by_area | is the area used as denominator |

### flowchart
divide_by_area does the following:  
Step 1: CHECK TO ENSURE SHAPE OF DATA IS THE SAME AS SHAPE OF AREA IF SHAPES MISMATCH, RETURN
Step 2: DIVIDE DATA BY GRID CELL AREA WHERE AREA /= 0.0


## send_ice_mask_sic
### intro
Module atm_land_ice_flux_exchange_mod handles flux exchange between atmosphere to land and ice.  send_ice_mask_sic is a subroutine in atm_land_ice_flux_exchange_mod.
### description
  If the variables ice_mask or sic have been registered with diag_manager, this subroutine, send_ice_mask_sic, maps the fractional amount of sea ice from the OCN grid to the ATM grid and sends the data to the diag_manager buffer.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| time | intent(inout) | send_ice_mask_sic | is the current model time |

### flowchart
send_ice_mask_sic does the following:  
Step 1: IF ID_ICE_MASK > 0 OR ID_SIC > 0
Step 2: INITIALIZE ICE_FRAC.
Step 3: MAP ICE_MASK FROM THE OCN GRID TO THE EXCHANGE GRID.
Step 4: MAP ICE_MASK FROM THE EXCHANGE GRID TO THE ATM GRID.
Step 5: IF ID_ICE_MASK > 0, SEND ICE_MASK TO THE DIAG_MANAGER BUFFER.
Step 6: FOR CMIP, IF ID_SIC > 0, COMPUTE SEA ICE FRACTIONAL AREA FOR ATM GRID CELLS THAT ARE OVER THE OCEAN AND NORMALIZE AREA BY THE FRACTION OF ATMOS GRID CELL THAT IS OCEAN.


## atm_stock_integrate
### intro
Module atm_land_ice_flux_exchange_mod handles flux exchange between atmosphere to land and ice.  atm_stock_integrate is a subroutine in atm_land_ice_flux_exchange_mod.
### description
  Subroutine atm_stock_integrate integrates over the total precipitation (liquid and frozen) in the atmosphere and multiply the integrated value by the timestep dt. This subroutine is called in flux_exchange_mod/flux_check_stocks.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| atm | intent(inout) | atm_stock_integrate | is the derived type holding the atmosphere boundary data |
| res | intent(inout) | atm_stock_integrate | is the integrated value |

### flowchart
atm_stock_integrate does the following:  
Step 1: CALL FMS_XGRID_STOCK_INTEGRATE

