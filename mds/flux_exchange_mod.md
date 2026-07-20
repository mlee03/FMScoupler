# flux_exchange_mod

Flux_exchange_mod is the top level module for flux exchange between components.

## variables

| Name | Type | Definition |
|------|------|------------|
| version | character(len=128) | is the program version string set automatically at compile time. |
| tag | character(len=128) | is a string set automatically at compile time. |
| do_init | logical | is a flag where if .TRUE., initialize module |
| bound_tol | real, parameter | is the tolerance value used when checking grid-boundary coordinate consistency. |
| d622 | real, parameter | is the ratio of dry-air and water-vapor gas constants. |
| d378 | real, parameter | is the complement of d622, used in humidity conversions. |
| z_ref_heat | real | is the reference height [m] for temperature and relative humidity diagnostics (t_ref, rh_ref, del_h, del_q). |
| z_ref_mom | real | is the reference height [m] for momentum diagnostics (u_ref, v_ref, del_m). |
| do_area_weighted_flux | logical | is a namelist flag where if .TRUE., normalize exchanged fluxes by the area; used in ice_ocean_flux_exchange. |
| debug_stocks | logical | is a namelist flag where if .TRUE., enable extra stock-conservation output for debugging. |
| divert_stocks_report | logical | is a namelist flag where if .TRUE., write stock reports 'stocks.out'; else write to stdout. |
| do_runoff | logical | is a namelist flag where if .TRUE., turn on the land runoff interpolation to the ocean |
| do_forecast | logical | is a namelist flag. |
| nblocks | integer | is a namelist variable for number of OpenMP blocks, defaults to 1. |
| partition_fprec_from_lprec | logical | is a namelist flag where if .TRUE., convert liquid precip to snow when t_ref is less than tfreeze parameter |
| tfreeze | real, parameter | is the freezing point of water at one atmosphere in [K]. |
| scale_precip_2d | logical | is a namelist flag where if .TRUE., rescale liquid precipitation using a 2-D field from data override. |
| gas_fluxes_initialized | logical | is a flag to indicate component fluxes have been initialized. |
| ex_gas_fields_atm | type(fmscoupler1dbc_type), target | is a derived type containing atmospheric surface variables that are used in calculating atmosphere-ocean gas fluxes. |
| ex_gas_fields_ice | type(fmscoupler1dbc_type), target | is a derived type containing ice-top and ocean surface variables that are used in calculating atmosphere-ocean gas fluxes. |
| ex_gas_fluxes | type(fmscoupler1dbc_type), target | is a derived type for exchanging gas or tracer fluxes between the atmosphere and ocean, defined by the field table. Also a place holder of intermediate calculations. |
| ni_atm | integer | is the number of x gridpoints in the atm compute domain |
| nj_atm | integer | is the number of y gridpoints in the atm compute domain |
| ccc | real, dimension(3) | is a temporary array used for conservation-check summaries; not used. |
| cplclock | integer | is the FMS clock id to profile land-ice-atmos coupler. |
| dt_atm | real | is the atmosphere timestep [s] |
| dt_cpl | real | is the coupled timesteps in [s]. |
| atm_precip_new | real | is used to take into account implicit evaporation in stock computation |
| mod_name | character(len=14), parameter |  |
| id_drag_moist | integer |  |
| id_drag_heat | integer |  |
| id_drag_mom | integer |  |
| id_rough_moist | integer |  |
| id_rough_heat | integer |  |
| id_rough_mom | integer |  |
| id_u_star | integer |  |
| id_b_star | integer |  |
| id_q_star | integer |  |
| id_u_flux | integer |  |
| id_v_flux | integer |  |
| id_t_surf | integer |  |
| id_t_flux | integer |  |
| id_q_flux | integer |  |
| id_r_flux | integer |  |
| id_t_atm | integer |  |
| id_u_atm | integer |  |
| id_v_atm | integer |  |
| id_wind | integer |  |
| id_thv_atm | integer |  |
| id_thv_surf | integer |  |
| id_t_ref | integer |  |
| id_rh_ref | integer |  |
| id_u_ref | integer |  |
| id_v_ref | integer |  |
| id_q_ref | integer |  |
| id_del_h | integer |  |
| id_del_m | integer |  |
| id_del_q | integer |  |
| id_albedo | integer |  |
| id_gust | integer |  |
| id_t_ca | integer |  |
| id_q_surf | integer |  |
| id_q_atm | integer |  |
| id_z_atm | integer |  |
| id_p_atm | integer |  |
| id_land_mask | integer |  |
| id_ice_mask | integer |  |
| id_rough_scale | integer |  |
| id_albedo_vis_dir | integer |  |
| id_albedo_nir_dir | integer |  |
| id_albedo_vis_dif | integer |  |
| id_albedo_nir_dif | integer |  |
| id_tas | integer |  |
| id_uas | integer |  |
| id_vas | integer |  |
| id_ts | integer |  |
| id_psl | integer |  |
| id_sfcwind | integer |  |
| id_tauu | integer |  |
| id_tauv | integer |  |
| id_hurs | integer |  |
| id_huss | integer |  |
| id_evspsbl | integer |  |
| id_hfls | integer |  |
| id_hfss | integer |  |
| id_height2m | integer |  |
| id_height10m | integer |  |
| first_static | logical |  |
| do_read_nml | logical |  |
| isphum | integer |  |
| n_atm_tr_tot | integer |  |
| n_atm_tr | integer |  |
| use_existing_grid_spec | logical |  |
| all_ocean | logical |  |
| all_land | logical |  |
| is | integer |  |
| ie | integer |  |
| js | integer |  |
| je | integer |  |
| t_surf | real, dimension(:,:), allocatable |  |
| t_ca | real, dimension(:,:), allocatable |  |
| q_surf | real, dimension(:,:), allocatable |  |
| p_surf | real, dimension(:,:), allocatable |  |
| e_t_n | real, dimension(:,:), allocatable |  |
| f_t_delt_n | real, dimension(:,:), allocatable |  |
| e_q_n | real, dimension(:,:), allocatable |  |
| f_q_delt_n | real, dimension(:,:), allocatable |  |
| dhdt_surf | real, dimension(:,:), allocatable |  |
| dedt_surf | real, dimension(:,:), allocatable |  |
| dedq_surf | real, dimension(:,:), allocatable |  |
| drdt_surf | real, dimension(:,:), allocatable |  |
| dhdt_atm | real, dimension(:,:), allocatable |  |
| dedq_atm | real, dimension(:,:), allocatable |  |
| flux_t | real, dimension(:,:), allocatable |  |
| flux_q | real, dimension(:,:), allocatable |  |
| flux_lw | real, dimension(:,:), allocatable |  |
| flux_u | real, dimension(:,:), allocatable |  |
| flux_v | real, dimension(:,:), allocatable |  |
| drag_q | real, dimension(:,:), allocatable |  |
| dtaudu_atm | real, dimension(:,:), allocatable |  |
| dtaudv_atm | real, dimension(:,:), allocatable |  |
| cd_t | real, dimension(:,:), allocatable |  |
| cd_m | real, dimension(:,:), allocatable |  |
| b_star | real, dimension(:,:), allocatable |  |
| u_star | real, dimension(:,:), allocatable |  |
| wind | real, dimension(:,:), allocatable |  |
| used | logical |  |


## gas_exchange_init
### intro
Flux_exchange_mod is the top level module for flux exchange between components.  gas_exchange_init is a subroutine in flux_exchange_mod.
### description
  Subroutine gas_exchange_init initializes the fms_atmos_ocean_type_fluxes, ocean_model_fluxes, and atmos_tracer_flux. The subroutine also calls fms_atmos_ocean_fluxes to initialize ex_gas_fluxes and fields.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| gas_fields_atm | intent(inout) | gas_exchange_init | is a derived type containing atmospheric surface variables tha are used in computing atmosphere-ocean gas fluxes. |
| gas_fields_ice | intent(inout) | gas_exchange_init | is a derived type containing ice-top and ocean surface variables that are used in computing atmosphere-ocean gas fluxes. |
| gas_fluxes | intent(inout) | gas_exchange_init | is a derived type for exchanging gas or tracer fluxes between the atmosphere and ocean, defined by the field table, as well as a place holder of intermediate calculations, such as piston velocities, and parameters that impact the fluxes. |

### flowchart
gas_exchange_init does the following:  
Step 1: CALL ATMOS_TRACER_FLUX_INIT(), OCEAN_MODEL_FLUX_INIT(), ATMOS_TRACER_FLUX_INIT(). ALSO CALLS FMS_ATMOS_OCEAN_FLUXES_INIT() TO ALLOCATE DERIVED TYPES.
Step 2: SET MODULE LEVEL GAS_FIELDS_ATM, GAS_FIELDS_ICE, AND GAS_FLUXES.


## flux_exchange_init
### intro
Flux_exchange_mod is the top level module for flux exchange between components.  flux_exchange_init is a subroutine in flux_exchange_mod.
### description
  Subroutine flux_exchange_init setups derived types and variables, and initializes modules that will be used in flux exchange. Ocean_tracer_flux_init is called first to get restart filenames for tracer fluxes for restart model runs. Atmos_tracer_flux_init is called last in order to use tracer values set in ocean_tracer_flux_init.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| time | intent(inout) | flux_exchange_init | is the model current time |
| atm | intent(inout) | flux_exchange_init | is a derived type to specify atm boundary data |
| land | intent(inout) | flux_exchange_init | is a derived type to specify land boundary data |
| ice | intent(inout) | flux_exchange_init | is a derived type to specify ice boundary data |
| ocean | intent(inout) | flux_exchange_init | is a derived type to specify ocean boundary data |
| ocean_state | intent(inout) | flux_exchange_init | is a pointer to the ocean model's internal state |
| atmos_ice_boundary | intent(inout) | flux_exchange_init | is a derived type holding properties and fluxes passed from atmosphere to ice |
| land_ice_atmos_boundary | intent(inout) | flux_exchange_init | is a derived type holding properties and fluxes passed from land and ice to atm |
| land_ice_boundary | intent(inout) | flux_exchange_init | is a derived type holding properties and fluxes passed from land to ice |
| ice_ocean_boundary | intent(inout) | flux_exchange_init | is a derived type holding properties and fluxes passed from ice to ocean |
| ocean_ice_boundary | intent(inout) | flux_exchange_init | is a derived type holding properties and fluxes passed from ocean to ice |
| do_ocean | intent(inout) | flux_exchange_init | is a flag indicating whether the ocean component is active |
| dt_atmos | intent(inout) | flux_exchange_init | is the atmosphere time step in [s] |
| dt_cpld | intent(inout) | flux_exchange_init | is the coupled time step in [s] |

### flowchart
flux_exchange_init does the following:  
Step 1: CALL FMS_SAT_VAPOR_PRES_INIT.
Step 2: SETUP OPENMP PARAMETERS.
Step 3: SET LOGFILE.
Step 4: READ FLUX_EXCHANGE_NML.
Step 5: WRITE NAMELIST TO LOGFILE.
Step 6: SET MODULE LEVEL DT_ATM AND DT_CPL TIMESTEPS.
Step 7: GET OCEAN MODEL GRID CELL AREAS FROM GRID_SPEC.
Step 8: IF ATMPE, CALL ATM_LAND_ICE_FLUX_EXCHANGE_INIT() AND LAND_ICE_FLUX_EXCHANGE_INIT() ALSO CHECK ATM_GRID CONSISTENCY WITH PROVIDED GRID_SPEC.
Step 9: CALL ICE_OCEAN_FLUX_EXCHANGE_INIT().
Step 10: SET DO_INIT TO .FALSE. TO SKIP INITIALIZATION IF FLUX_EXCHANGE_INIT IS CALLED AGAIN.


## flux_check_stocks
### intro
Flux_exchange_mod is the top level module for flux exchange between components.  flux_check_stocks is a subroutine in flux_exchange_mod.
### description
  Subroutine flux_check_stocks computes the current stock values for atm, land, ice, and ocean; and outputs the stock differences with respect to the initial values in the logfile.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| time | intent(inout) | flux_check_stocks | is the model's current time |
| atm | intent(inout) | flux_check_stocks | is the atmosphere boundary data type used to compute atmosphere stocks |
| lnd | intent(inout) | flux_check_stocks | is the land boundary data type used to compute land stocks |
| ice | intent(inout) | flux_check_stocks | is the ice boundary data type used to compute ice stocks |
| ocn_state | intent(inout) | flux_check_stocks | is a pointer to the ocean model's internal state used to compute ocean stocks |

### flowchart
flux_check_stocks does the following:  
Step 1: FOR WATER, HEAT, AND SALT STOCKS FOR EACH COMPONENT, GET CURRENT STOCK VALUE AND COMPARE WITH INTEGRATED FLUXES FOR ATM WATER STOCK. FOR ATM, INTEGRATE ATM_PRECIP_NEW FOR IMPLICIT EVAPORATION.
Step 2: PRINT FOR EACH ELEMENT, S(t): TOTAL STOCK, S(t)-S(0): CHANGE IN STOCK WITH RESPECT TO INITIAL VALUE, F(t): CUMULATIVE FLUX INTO COMPONENT FROM OTHER COMPONENTS F(t) - [S(t)-S(0)]: DIFFERENCE BETWEEN THE FLUXES AND STOCK CHANGE (S(t)-S(0))/F(t): RELATIVE ERROR


## flux_init_stocks
### intro
Flux_exchange_mod is the top level module for flux exchange between components.  flux_init_stocks is a subroutine in flux_exchange_mod.
### description
  Subroutine flux_init_stocks initializes the stock values for the atmosphere, land, ice, and ocean. Stocks are the globally integrated total amoun of conserved quantities such as mass and energy and is used to check conservation.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| time | intent(inout) | flux_init_stocks | is the model's current time |
| atm | intent(inout) | flux_init_stocks | is a derived type holding atmosphere boundary data |
| lnd | intent(inout) | flux_init_stocks | is a derived type holding land boundary data |
| ice | intent(inout) | flux_init_stocks | is a derived type holding ice boundary data |
| ocn_state | intent(inout) | flux_init_stocks | is a pointer to ocean model's internal state |

### flowchart
flux_init_stocks does the following:  
Step 1: IF DIVERT_STOCKS_REPORT IS FALSE, OPEN STOCKS OUTPUT FILE TO STDOUT. IF DIVERT_STOCKS_REPORT IS TRUE, OPEN STOCKS OUTPUT FILE TO "stocks.out". ONLY THE ROOT PE WILL WRITE TO THE FILE.
Step 2: INITIALIZE WATER, HEAT, AND SALT STOCK VALUES FOR EACH COMPONENT. FOR ATMOSPHERE, INTEGRATE ATM_PRECIP_NEW TO GET THE INITIAL ISTOCK_WATER.
Step 3: INITIALIZE STOCKS IN FMS.


## check_atm_grid
### intro
Flux_exchange_mod is the top level module for flux exchange between components.  check_atm_grid is a subroutine in flux_exchange_mod.
### description
  Subroutine check_atm_grid checks the consistency of the atmosphere grid specified in the model with the grid specified in the grid_file.

atm


is a derived type holding atmosphere boundary and grid data



grid_file


is the path to the grid specification file.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| atm | intent(inout) | check_atm_grid | is a derived type holding atmosphere boundary and grid data |
| grid_file | intent(inout) | check_atm_grid | is the path to the grid specification file |

### flowchart
check_atm_grid does the following:  
Step 1: GET GLOBAL, COMPUTE, AND DATA DOMAIN INDICES AND SIZES FOR THE ATMOSPHERE COMPONENT.
Step 2: OPEN GRID_FILE.
Step 3: CHECK GRID SIZES ARE CONSISTENT.
Step 4: CHECK LON, LAT, AND GRID CELL AREAS ARE CONSISTENT.


## sfc_boundary_layer
### intro
Flux_exchange_mod is the top level module for flux exchange between components.  sfc_boundary_layer is a subroutine in flux_exchange_mod.
### description
  
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| dt | intent(inout) | sfc_boundary_layer | Time step |
| time | intent(inout) | sfc_boundary_layer | Current time |
| atm | intent(inout) | sfc_boundary_layer | A derived data type to specify atmospheric boundary data |
| land | intent(inout) | sfc_boundary_layer | A derived data type to specify land boundary data |
| ice | intent(inout) | sfc_boundary_layer | A derived data type to specify ice boundary data |
| boundary | intent(inout) | sfc_boundary_layer | A derived data type to specify properties and fluxes passed from exchange grid to the atmosphere, |

### flowchart
sfc_boundary_layer does the following:  


## flux_down_from_atmos
### intro
Flux_exchange_mod is the top level module for flux exchange between components.  flux_down_from_atmos is a subroutine in flux_exchange_mod.
### description
  
### arguments
None
### flowchart
flux_down_from_atmos does the following:  


## flux_up_to_atmos
### intro
Flux_exchange_mod is the top level module for flux exchange between components.  flux_up_to_atmos is a subroutine in flux_exchange_mod.
### description
  
### arguments
None
### flowchart
flux_up_to_atmos does the following:  


## flux_exchange_init
### intro
Flux_exchange_mod is the top level module for flux exchange between components.  flux_exchange_init is a subroutine in flux_exchange_mod.
### description
  
### arguments
None
### flowchart
flux_exchange_init does the following:  


## read_namelist
### intro
Flux_exchange_mod is the top level module for flux exchange between components.  read_namelist is a subroutine in flux_exchange_mod.
### description
  
### arguments
None
### flowchart
read_namelist does the following:  


## diag_field_init
### intro
Flux_exchange_mod is the top level module for flux exchange between components.  diag_field_init is a subroutine in flux_exchange_mod.
### description
  
### arguments
None
### flowchart
diag_field_init does the following:  


## flux_exchange_end
### intro
Flux_exchange_mod is the top level module for flux exchange between components.  flux_exchange_end is a subroutine in flux_exchange_mod.
### description
  
### arguments
None
### flowchart
flux_exchange_end does the following:  


## surface_flux_2d
### intro
Flux_exchange_mod is the top level module for flux exchange between components.  surface_flux_2d is a subroutine in flux_exchange_mod.
### description
  
### arguments
None
### flowchart
surface_flux_2d does the following:  

