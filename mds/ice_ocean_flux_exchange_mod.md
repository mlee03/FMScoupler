# ice_ocean_flux_exchange_mod

Module ice_ocean_flux_exchange_mod handles data transfer between ice and ocean as well as stock computation.

## variables

| Name | Type | Definition |
|------|------|------------|
| regrid | integer, parameter | is a flag used to indicate ice and ocean are on physically different grids. Data will be transferred via the exchange grid |
| redist | integer, parameter | is a flag used to indicate grids for ocean and ice are same but with differen domain decomposition. Data will be transferred with fms_mpp_redistribute. |
| direct | integer, parameter | is a flag used to indicate grids for ocean and ice are same. Data can be copied directly. |
| debug_stocks | logical | is a flag where if .TRUE., call check_flux_conservation at module initialization |
| do_area_weighted_flux | logical | is a flag where if .TRUE., scale fluxes by source cell area and divide by the destination cell area to preserve the global area-weighted integral if redistributing. |
| cplocnclock | integer | is a FMS clock ID to time flux_ice_to_ocean and flux_ocean_to_ice |
| fluxoceaniceclock | integer | is a FMS clock ID to time flux_ocean_to_ice transfer |
| fluxiceoceanclock | integer | is a FMS clock ID to time flux_ice_to_ocean transfer |
| dt_cpl | real | is the Coupled (slow) timestep in seconds; used in stock computation |
| slow_ice_ocean_pelist | integer, dimension(:), allocatable | is the Combined MPI pelist of the slow-ice and ocean pes; set during initialization |


## ice_ocean_flux_exchange_init
### intro
Module ice_ocean_flux_exchange_mod handles data transfer between ice and ocean as well as stock computation.  ice_ocean_flux_exchange_init is a subroutine in ice_ocean_flux_exchange_mod.
### description
  Subroutine ice_ocean_flux_exchange_init initializes the module for flux exchange between ice and ocean.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| time | intent(inout) | ice_ocean_flux_exchange_init | is the model's current time |
| ice | intent(inout) | ice_ocean_flux_exchange_init | is a derived data type holding ice boundary data |
| ocean | intent(inout) | ice_ocean_flux_exchange_init | is a derived data type holding ocean boundary data |
| ocean_state | intent(inout) | ice_ocean_flux_exchange_init | is a pointer pointing to the ocean model's internal state |
| ice_ocean_boundary | intent(inout) | ice_ocean_flux_exchange_init | is a derived data type holding properties and fluxes passed from ice to ocean |
| ocean_ice_boundary | intent(inout) | ice_ocean_flux_exchange_init | is a derived data type holding properties and fluxes passed from ocean to ice |
| dt_cpl_in | intent(inout) | ice_ocean_flux_exchange_init | is the coupled (slow) timestep in seconds to set module level Dt_cpl that's used in stock computation. |
| debug_stocks_in | intent(inout) | ice_ocean_flux_exchange_init | is used to set module level variable debug_stocks. If TRUE, stocks will be computed for flux exchange consistency |
| do_area_weighted_flux_in | intent(inout) | ice_ocean_flux_exchange_init | is used to set module level do_area_weighted_flux. If TRUE, flux between ice and ocean will be area-weighted |
| ex_gas_fluxes | intent(inout) | ice_ocean_flux_exchange_init | is used to spawn matching arrays in ice_ocean_boundary and Iceocean_fluxes. |
| do_ocean | intent(inout) | ice_ocean_flux_exchange_init | is a flag where if .TRUE., ocean_ice_boundarystagger = oceanstagger, else defaults to AGRID |
| slow_ice_ocean_pelist_in | intent(inout) | ice_ocean_flux_exchange_init | is the combined MPI pelist of the slow-ice and ocean processing element used to set module level slow_ice_ocean_pelis |

### flowchart
ice_ocean_flux_exchange_init does the following:  
Step 1: INITIALIZE OCEAN_ICE_BOUNDARY FIELDS TO ZERO. INITIALIZE T TO 273.0 [K].
Step 2: SPAWN GAS FIELDS TO OCEAN_ICE_BOUNDARYFIELDS.
Step 3: SPAWN GAS FLUXES TO ICEOCEAN_FLUXES.
Step 4: ALLOCATE ICE_OCEAN_BOUNDARY FIELDS AND INITIALIZE TO ZERO. IF ICEBERG FIELDS ARE ASSOCIATED IN ICE, ALLOCATE ICE_BERGS FIELDS IN ICE_OCEAN_BOUNDARY AND INITIALIZE TO ZERO.
Step 5: SPAWN GAS FIELDS AND FLUXES TO ICE_OCEAN_BOUNDARYFLUXES.
Step 6: SPAWN GAS FIELDS TO OCEANFIELDS.
Step 7: INITIALIZE BOUNDARY VALUES OCEAN_ICE_BOUNDARYXTYPE TO DIRECT IF THE ICE AND OCEAN DOMAINS ARE THE SAME, OTHERWISE REDIST. (USED IN DATA_OVERRIDE)
Step 8: CALL OCEAN_MODEL_INIT_SFC TO COMPLETE OCEAN SURFACE FIELD INITIALIZATION.
Step 9: CHECK FLUX CONSERVATION IF DEBUG_STOCKS IS TRUE.
Step 10: ALLOCATE SLOW_ICE_OCEAN_PELIST. INITIALIZE CLOCKS TO MEASURE PERFORMANCE.


## flux_ice_to_ocean
### intro
Module ice_ocean_flux_exchange_mod handles data transfer between ice and ocean as well as stock computation.  flux_ice_to_ocean is a subroutine in ice_ocean_flux_exchange_mod.
### description
  Subroutine flux_ice_to_ocean interpolates data from Ice to Ice_Ocean_Boundary in order to exchange fluxes at the bottom of ice to the ocean model. The following quantities are transferred from the Ice to Ice_Ocean_Boundary: flux_u = zonal wind stress [Pa] flux_v = meridional wind stress [Pa] flux_t = sensible heat flux [W/m2] flux_q = specific humidity flux [Kg/m2/s] flux_salt = salt flux [Kg/m2/s] flux_sw = net (down-up) shortwave flux [W/m2] flux_lw = net (down-up) longwave flux [W/m2] lprec = mass of liquid precipitation since last time step [Kg/m2] fprec = mass of frozen precipitation since last time step [Kg/m2] runoff = mass of runoff since last time step [Kg/m2] calving = mass of calving since last time step [Kg/m2] p_surf = surface pressure [Pa].
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| ice | intent(inout) | flux_ice_to_ocean | is a derived data type containg ice boundary data |
| ocean | intent(inout) | flux_ice_to_ocean | is a derived data type to containing ocean boundary data |
| ice_ocean_boundary | intent(inout) | flux_ice_to_ocean | is a derived data type to specify properties and fluxes passed from ice to ocean |

### flowchart
flux_ice_to_ocean does the following:  


## flux_ice_to_ocean_finish
### intro
Module ice_ocean_flux_exchange_mod handles data transfer between ice and ocean as well as stock computation.  flux_ice_to_ocean_finish is a subroutine in ice_ocean_flux_exchange_mod.
### description
  Subroutine flux_ice_to_ocean_finish mainly calls fms_data_override to override fluxes in Ice_Ocean_Boundary before transferring flux from Ice to Ocean. NOTE, fms_data_override will only override data if field entry is found in the data_table. This subroutine is only called by the ocean pe.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| time | intent(inout) | flux_ice_to_ocean_finish | is the current time |
| ice_ocean_boundary | intent(inout) | flux_ice_to_ocean_finish | is a derived data type containing fluxes and properties passed from ice to ocean |

### flowchart
flux_ice_to_ocean_finish does the following:  


## flux_ocean_to_ice
### intro
Module ice_ocean_flux_exchange_mod handles data transfer between ice and ocean as well as stock computation.  flux_ocean_to_ice is a subroutine in ice_ocean_flux_exchange_mod.
### description
  Subroutine flux_ocean_to_ice interpolates data from Ocean to Ocean_Ice_Boundary in order to exchange fluxes from ocean to bottom of ice. The following quantities are remapped from the Ocean to Ocean_Ice_Boundary: t_surf = surface temperature [deg K] frazil = frazil fluxes since the last coupling step [J/m2] u_surf = zonal ocean current/ice motion [m/s] v_surf = meridional ocean current/ice motion [m/s] sea_lev = sea level used to drive ice accelerations [m].
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| ocean | intent(inout) | flux_ocean_to_ice | is a derived data type holding ocean boundary data |
| ice | intent(inout) | flux_ocean_to_ice | is a derived data type holding ice boundary data |
| ocean_ice_boundary | intent(inout) | flux_ocean_to_ice | is a derived data type holding properties and fluxes passed from ocean to ice |

### flowchart
flux_ocean_to_ice does the following:  


## flux_ocean_to_ice_finish
### intro
Module ice_ocean_flux_exchange_mod handles data transfer between ice and ocean as well as stock computation.  flux_ocean_to_ice_finish is a subroutine in ice_ocean_flux_exchange_mod.
### description
  Subroutine flux_ocean_to_ice_finish carrries out a final set of tasks that should only occur on the slow-ice processors, including data override and perhaps saving diagnostics.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| time | intent(inout) | flux_ocean_to_ice_finish | is the current time |
| ice | intent(inout) | flux_ocean_to_ice_finish | is a derived type holding ice boundary data |
| ocean_ice_boundary | intent(inout) | flux_ocean_to_ice_finish | is a derived type holding properties and fluxes passed from ocean to ice |

### flowchart
flux_ocean_to_ice_finish does the following:  


## flux_ice_to_ocean_stocks
### intro
Module ice_ocean_flux_exchange_mod handles data transfer between ice and ocean as well as stock computation.  flux_ice_to_ocean_stocks is a subroutine in ice_ocean_flux_exchange_mod.
### description
  Subroutine flux_ice_to_ocean_stocks integrates the fluxes from ice to ocean over the surface and in time. Ice stocks are decremented at the base of the ice and incremented to the ocean stocks at the ocean surface.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| ice | intent(inout) | flux_ice_to_ocean_stocks | A derived data type to specify ice boundary data |

### flowchart
flux_ice_to_ocean_stocks does the following:  
Step 1: COMPUTE STOCKS CHANGE FOR QUANTITY (PRECIP - EVAP).
Step 2: COMPUTE STOCKS FOR RIVER.
Step 3: COMPUTE STOCKS FOR HEAT (SENSIBLE + SHORTWAVE + LONGWAVE + LATENT).
Step 4: COMPUTE STOCKS FOR HEAT FROM RADIATIVE AND TURBLENT FLUXES AND HEAT CARRIED BY RIVER AND PME (assuming reference temperature of 0 degC and river/pme temp = surface temp). Note: it does not matter what the ref temperature is but it must be consistent with that in OCN and ICE. PME = preciptation minus evaporation
Step 5: COMPUTE STOCKS FOR FLUX_SALT.


## flux_ocean_from_ice_stocks
### intro
Module ice_ocean_flux_exchange_mod handles data transfer between ice and ocean as well as stock computation.  flux_ocean_from_ice_stocks is a subroutine in ice_ocean_flux_exchange_mod.
### description
  Subroutine flux_ocean_from_ice_stocks updates stocks in Ocean after flux transfer from Ice. Unlike subroutine flux_ice_to_ocean_stocks() that uses Icefluxes to update the stocks, this subroutine uses Ice_Ocean_boundaryfluxes to calculate the amount of input to Ocean. These fluxes are the ones that Ocean model uses internally to calculate its budgets. Hence there should be no difference between this input and what Ocean model internal diagnostics uses. This bypasses the possible mismatch in cell areas between Ice and Ocean in diagnosing the stocks of Ocean and should report a conserving Ocean component regardless of the glitches in fluxes. The use of this subroutine in conjunction with subroutine flux_ice_to_ocean_stocks() will also allow to directly diagnose the amount "stocks lost in exchange" between Ice and Ocean.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| ocean_state | intent(inout) | flux_ocean_from_ice_stocks | is a pointer to the ocean model's internal state; used to retrieve ocean-side grid and flux data via ocean_model_data_get. |
| ocean | intent(inout) | flux_ocean_from_ice_stocks | is a derived type containing the Ocean public boundary data type; provides the MPI domain and ocean pe information to get compute domain. |
| ice_ocean_boundary | intent(inout) | flux_ocean_from_ice_stocks | is a derived type containing fluxes passed from ice to ocean. |

### flowchart
flux_ocean_from_ice_stocks does the following:  
Step 1: USE THE RETRIEVER FROM OCEAN_MODEL_MOD TO GET AREA, MASK, SURFACE TEMPERATURE, PME TEMPERATURE, CALVING TEMPERATURE, RUNOFF TEMPERATURE, BOTTOM HEAT FLUX, AND SPECIFIC HEAT CAPACITY FIELDS FROM THE OCEAN MODEL.
Step 2: COMPUTE STOCK TRANSFER OF (PRECIP - EVAP) TO OCEAN SURFACE AND FROM LATERALLY.
Step 3: COMPUTE STOCK TRANSFER OF (SENSIBLE HEAT + SHORTWAVE + LONGWAVE + LATENT HEAT) TO OCEAN LATERALLY.
Step 4: COMPUTE STOCK TRANSFER OF HEAT CARRIED BY RIVER + PME (ASSUMING REFERENCE TEMPERATURE OF 0 DEGC AND RIVER/PME TEMP = SURFACE TEMP).
Step 5: COMPUTE STOCK TRANSFER OF BOTTOM HEAT FLUX.
Step 6: COMPUTE STOCK TRANSFER OF FRAZIL HEAT.
Step 7: COMPUTE STOCK TRANSFER OF SALT FLUX.


## flux_ice_to_ocean_redistribute
### intro
Module ice_ocean_flux_exchange_mod handles data transfer between ice and ocean as well as stock computation.  flux_ice_to_ocean_redistribute is a subroutine in ice_ocean_flux_exchange_mod.
### description
  Subroutine flux_ice_to_ocean_redistribute performs a globally conservative flux redistribution across ICE/OCN. If domain decomposition is identical for ocean and ice, data is copied from Ice to ICE_OCEAN_BOUNDARY If domain decomposition differs, data is copied from Ice to ICE_OCEAN_BOUNDARY with fms_mpp_domains_redistribute. (Assumes Ice and Ocean are on the same grid.).
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| ice | intent(inout) | flux_ice_to_ocean_redistribute | is the ice boundary data type; provides the slow-ice MPI domain (slow_Domain_NH) and cell areas (Icearea) used in redistribution. |
| ocean | intent(inout) | flux_ice_to_ocean_redistribute | is the ocean public boundary data type; provides the ocean MPI domain (OceanDomain) and cell areas (Oceanarea) used in redistribution. |
| ice_data | intent(inout) | flux_ice_to_ocean_redistribute | is the flux field on the ice domain to be transferred to the ocean boundary. |
| ocn_bnd_data | intent(inout) | flux_ice_to_ocean_redistribute | is the flux field on the ocean domain; filled with the redistributed (and optionally area-weighted) values from ice_data. |
| type | intent(inout) | flux_ice_to_ocean_redistribute | is the transfer type: DIRECT (same MPI decomposition, copy directly) or REDIST (same grid, different MPI decomposition, use mpp_redistribute). |
| do_area_weighted | intent(inout) | flux_ice_to_ocean_redistribute | is a flag where if .TRUE., scale flux by ice cell area before redistribution and divide by ocean cell area after, preserving the global area-weighted integral. |

### flowchart
flux_ice_to_ocean_redistribute does the following:  


## divide_by_area
### intro
Module ice_ocean_flux_exchange_mod handles data transfer between ice and ocean as well as stock computation.  divide_by_area is a subroutine in ice_ocean_flux_exchange_mod.
### description
  Subroutine divide_by_area divides data by area for area > 0.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| data | intent(inout) | divide_by_area | is the data to divide by area; modified in-place |
| area | intent(inout) | divide_by_area | is the area field to divide by |

### flowchart
divide_by_area does the following:  
Step 1: IF DATA AND AREA DIFFER IN SIZE, RETURN WITHOUT MODIFYING DATA
Step 2: WHERE(AREA /= 0.0) DATA = DATA / AREA


## check_flux_conservation
### intro
Module ice_ocean_flux_exchange_mod handles data transfer between ice and ocean as well as stock computation.  check_flux_conservation is a subroutine in ice_ocean_flux_exchange_mod.
### description
  Subroutine check_flux_conservation checks for flux conservation after flux_ice_to_ocean_redistribute.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| ice | intent(inout) | check_flux_conservation | is the Ice boundary data type; provides the ice MPI domain, cell areas (Icearea), and flux array sizes used to allocate test data. |
| ocean | intent(inout) | check_flux_conservation | is the Ocean public boundary data type; provides the ocean MPI domain and cell areas (Oceanarea) used to compute redistributed sums. |
| ice_ocean_boundary | intent(inout) | check_flux_conservation | is the Ice-to-ocean boundary type; provides xtype (DIRECT or REDIST) and q_flux array size used to allocate the ocean-side test buffer. |

### flowchart
check_flux_conservation does the following:  
Step 1: SET OUTUNIT TO STDOUT.
Step 2: ALLOCATE ICE_DATA AND OCN_DATA FOR TESTING.
Step 3: INITIALIZE ICE_DATA WITH RANDOM NUMBERS.
Step 4: CALL FLUX_ICE_TO_OCEAN_DISTRIBUTE WITH AREA_WEIGHTED_SUM = .FALSE. AND GET GLOBAL SUM.
Step 5: CALL FLUX_ICE_TO_OCEAN_DISTRIBUTE WITH AREA_WEIGHTED_SUM = .TRUE. AND GET GLOBAL SUM.
Step 6: WRITE REPORT TO OUTUNIT.

