# atmos_ocean_dep_fluxes_calc_mod

Module atmos_ocean_dep_fluxes_calc_mod handles computation of ocean and atmosphere deposition gas fluxes.

## variables

| Name | Type | Definition |
|------|------|------------|
| mod_name | character(len= *), parameter | mod_name used when printing error messages |


## atmos_ocean_dep_fluxes_calc
### intro
Module atmos_ocean_dep_fluxes_calc_mod handles computation of ocean and atmosphere deposition gas fluxes.  atmos_ocean_dep_fluxes_calc is a subroutine in atmos_ocean_dep_fluxes_calc_mod.
### description
  Subroutine atmos_ocean_dep_fluxes_calc calculates atmosphere to ocean wet and dry deposition fluxes.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| gas_fields_atm | intent(inout) | atmos_ocean_dep_fluxes_calc | is a derived type containing atmospheric surface variables |
| gas_fields_ice | intent(inout) | atmos_ocean_dep_fluxes_calc | is a derived type containing ice-top and ocean surface variables |
| gas_fluxes | intent(inout) | atmos_ocean_dep_fluxes_calc | is a derived type containing gas fluxes between the atmosphere and the ocean, and related parameters |
| seawater | intent(inout) | atmos_ocean_dep_fluxes_calc | is a mask with value of 1 for the open water, 0 if ice or land. |

### flowchart
atmos_ocean_dep_fluxes_calc does the following:  
Step 1: RETURN IF THE NUMBER OF GAS FLUXES AT BOUNDARY IS ZERO.
Step 2: ERROR IF GAS FLUXES BC ARRAY IS NOT ASSOCIATED.
Step 3: COMPUTE DEPOSITION FLUXES IF FLUX WAS NOT OVERRIDDEN BY DATA_OVERRIDE AND IF FLUX TYPE IS AIR-SEA-DEPOSITION.
Step 4: CALCULATE DEPOSITION FLUXES FOR OPEN WATER CELLS. SET FLUXES TO ZERO FOR ICE AND LAND CELLS.

