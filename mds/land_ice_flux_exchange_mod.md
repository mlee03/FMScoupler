# land_ice_flux_exchange_mod

Module land_ice_flux_exchange_mod handles freshwater discharge (runoff and calving) and associated heat exchanges via the exchange grid between land and ice grids.

## variables

| Name | Type | Definition |
|------|------|------------|
| xmap_runoff | type(fmsxgridxmap_type), save | is the exchange grid map between land and ice/ocean. |
| n_xgrid_runoff | integer | is the number of exchange grid cells in xmap_runoff. |
| x2_grid_lnd | integer | is the land index for xmap_runoff; used to identify the source side when calling fms_xgrid_stock_move. set to 1 |
| x2_grid_ice | integer | is the ice/ocean index for xmap_runoff; used to identify the destination side when calling fms_xgrid_stock_move. set to 2 |
| cplclock | integer | is the clock ID for timing flux_land_to_ice calls. |
| fluxlandiceclock | integer | is the clock ID for timing the flux_land_to_ice transfer |
| do_runoff | logical | is a flag where if .TRUE., land discharge is transferred to ice. If .FALSE., all runoff/calving fields are zeroed. |
| dt_cpl | real | is the coupled (slow) timestep in seconds; used in stock computation. |


## land_ice_flux_exchange_init
### intro
Module land_ice_flux_exchange_mod handles freshwater discharge (runoff and calving) and associated heat exchanges via the exchange grid between land and ice grids.  land_ice_flux_exchange_init is a subroutine in land_ice_flux_exchange_mod.
### description
  Subrutine land_ice_flux_exchange_init initializes the land-ice flux exchange module for flux exchange between land and ice.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| land | intent(inout) | land_ice_flux_exchange_init | is a derived type holding land boundary data |
| ice | intent(inout) | land_ice_flux_exchange_init | is a derived type holding ice boundary data |
| land_ice_boundary | intent(inout) | land_ice_flux_exchange_init | is a derived type holding properties and fluxes passed from land to ice |
| dt_cpl_in | intent(inout) | land_ice_flux_exchange_init | is the coupled (slow) timestep in seconds to set module level dt_cpl |
| do_runoff_in | intent(inout) | land_ice_flux_exchange_init | is a flag to set module level do_runoff. |
| cplclock_in | intent(inout) | land_ice_flux_exchange_init | is the FMS MPP clock id for the top-level coupler profiling clock; stored module-wide so tha flux_land_to_ice can bracket its work. |

### flowchart
land_ice_flux_exchange_init does the following:  
Step 1: SET DO_RUNOFF, CPLCLOCK, DT_CPL, AND FLUXLANDICECLOCK.
Step 2: ALLOCATE LAND_ICE_BOUNDARYRUNOFF, CALVING, RUNOFF_HFLX, AND CALVING_HFLX, AND INITIALIZE THEM TO ZERO.


## flux_land_to_ice
### intro
Module land_ice_flux_exchange_mod handles freshwater discharge (runoff and calving) and associated heat exchanges via the exchange grid between land and ice grids.  flux_land_to_ice is a subroutine in land_ice_flux_exchange_mod.
### description
  Subroutine flux_land_to_ice handles conservative transfer of water and snow discharge from land to sea ice/ocean. The following elements are transferred from the Land to the Land_ice_boundary: discharge to runoff (kg/m2). discharge_snow to calving (kg/m2).
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| time | intent(inout) | flux_land_to_ice | is the current time |
| land | intent(inout) | flux_land_to_ice | is a derived type holding land boundary data |
| ice | intent(inout) | flux_land_to_ice | is a derived type holding ice boundary data |
| land_ice_boundary | intent(inout) | flux_land_to_ice | is a derived type holding properties and fluxes passed from land to ice |

### flowchart
flux_land_to_ice does the following:  
Step 1: INITIALIZE CLOCK.
Step 2: IF DO_RUNOFF, TRANSFER DATA, ELSE SET LAND_ICE_BOUNDARYRUNOFF, CALVING RUNOFF_HFLX, AND CALVING_HFLX TO ZERO.
Step 3: TRANSFER DISCHARGE* FIELDS FROM THE LAND TO ICE VIA THE EXCHANGE GRID.
Step 4: OVERRIDE TRANSFERRED DATA WITH DATA_OVERRIDE IF FIELD EXISTS IN DATA_TABLE.
Step 5: COMPUTE WATER STOCK ON THE EXCHANGE GRID TO MEASURE STOCK BEING TRANSFERRED FROM LAND TO ICE.
Step 6: END CLOCK.

