# atmos_ocean_fluxes_calc_mod

Module atmos_ocean_fluxes_calc_mod calculates gas fluxes between atmosphere and ocean.

## variables

| Name | Type | Definition |
|------|------|------------|
| mod_name | character(len= *), parameter | is the module name used when printing error messages |
| epsln | real, parameter | is a really small number used to prevent divide-by-zero |


## atmos_ocean_fluxes_calc
### intro
Module atmos_ocean_fluxes_calc_mod calculates gas fluxes between atmosphere and ocean.  atmos_ocean_fluxes_calc is a subroutine in atmos_ocean_fluxes_calc_mod.
### description
  Subroutine atmos_ocean_fluxes_calc calculates atmos-ocean gas fluxes. All fluxes are in units of [mol/m^2/s] with values > 0 for upward flux. Deposition fluxes are computed in atmos_ocean_dep_fluxes_calc. All calculations are done on the exchange grid.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| gas_fields_atm | intent(inout) | atmos_ocean_fluxes_calc | is a derived type containing atmospheric surface variables |
| gas_fields_ice | intent(inout) | atmos_ocean_fluxes_calc | is a derived type containing ice-top and ocean surface variables |
| gas_fluxes | intent(inout) | atmos_ocean_fluxes_calc | is a derived type containing gas fluxes between the atmosphere and the ocean, and parameters |
| seawater | intent(inout) | atmos_ocean_fluxes_calc | is a mask with value of 1 for the open water, 0 if ice or land. |
| tsurf | intent(inout) | atmos_ocean_fluxes_calc | is the sea-surface temperature [K]; used to compute gas-phase and liquid-phase transfer velocities. |
| ustar | intent(inout) | atmos_ocean_fluxes_calc | is the friction velocity [m/s]. When provided, overrides the internally computed value u_{10}{C_D} used inside calc_kw. |
| cd_m | intent(inout) | atmos_ocean_fluxes_calc | is the drag coefficient [dimensionless]. Only used when ustar is provided; otherwise calc_kw uses the bulk formula C_D = 6.1*10^{-4}+0.63*10^{-4}u_{10}. |

### flowchart
atmos_ocean_fluxes_calc does the following:  
Step 1: RETURN IF THE NUMBER OF GAS FLUXES AT ATMOSPHERE AND OCEAN BOUNDARY IS ZERO.
Step 2: COMPUTE FLUXES AT BOUNDARY AS FOLLOWS: OCMIP2, DUCE, OR JOHNSON IMPLEMENTATIONS FOR AIR_SEA_GAS_FLUX_GENERIC FLUXES; OCMIP2, OCMIP2_DATA, OR LINEAR IMPLEMENTATIONS FOR AIR_SEA_GAS_FLUX FLUXES; RIVER IMPLEMENTATION FOR LAND_SEA_RUNOFF FLUXES. AIR_SEA_DEPOSITION FLUXES ARE COMPUTED ELSEWHERE IN ATMOS_OCEAN_DEP_FLUXES_MOD.


## calc_kw
### intro
Module atmos_ocean_fluxes_calc_mod calculates gas fluxes between atmosphere and ocean.  calc_kw is a real function in atmos_ocean_fluxes_calc_mod.
### description
  Function calc_kw calculates the total transfer velocities from the "point of view" of liquid following Johnson Implementation from Johnson, Ocean Science, 2010. (http://doi.org/10.5194/os-6-913-2010) Uses equations defined in Liss[1974], F = K_g(c_g - H C_l) = K_l(c_g/H - C_l) where, F is the flux of gas across air-water interface, c_g and C_l are the bulk gas and liquid concentrations, H is the Henry's law constant (H = c_{sg}/C_{sl}), C_{sg} is the equilibrium concentration in gas phase [g/cm^3 of air] and C_{sl} is the equilibrium concentration of un-ionized dissolved gas in liquid phase [g/cm^3of water]), and K_g and K_l are the gas-phase and liquid-phase exchange constants, respectively. 1/K_g = 1/k_g + H/k_l 1/K_l = 1/k_l + 1/(H*k_g).
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| tk | intent(inout) | calc_kw | is the temperature at surface [K] |
| p | intent(inout) | calc_kw | is the pressure at surface [Pa] |
| u10 | intent(inout) | calc_kw | is the wind speed at 10m above the surface [m/s] |
| h | intent(inout) | calc_kw | is the Henry's law constant (H = c_sg/C_sl) (unitless) |
| vb | intent(inout) | calc_kw | is the Molar volume [m^3/mol] |
| mw | intent(inout) | calc_kw | is the molecular weight [g/mol] |
| sc_w | intent(inout) | calc_kw | is the Schmidt number [dimensionless] used to scale the liquid-phase piston velocity k_l relative to the reference Schmidt number of 660 (CO2 at 20 °C). |
| ustar | intent(inout) | calc_kw | is the Friction velocity [m/s]. If not provided, ustar = u_{10}*sqrt{C_D}. |
| cd_m | intent(inout) | calc_kw | is the Drag coefficient ($C_D). Used only if ustar is provided. If ustar is not provided, cd_m = 6.1x10^{-4} + 0.63x10^{-4} * u_10 |

### flowchart
calc_kw does the following:  


## calc_ka
### intro
Module atmos_ocean_fluxes_calc_mod calculates gas fluxes between atmosphere and ocean.  calc_ka is a real function in atmos_ocean_fluxes_calc_mod.
### description
Function calc_ka calculates total transfer velocities from the "point of view" of gas following Johnson Implementation from Johnson, Ocean Science, 2010. (http://doi.org/10.5194/os-6-913-2010) Uses equations defined in Liss[1974], F = K_g(c_g - H C_l) = K_l(c_g/H - C_l) where, F is the flux of gas across air-water interface, c_g and C_l are the bulk gas and liquid concentrations, H is the Henry's law constant (H = c_{sg}/C_{sl}), C_{sg} is the equilibrium concentration in gas phase [g/cm^3 of air] and C_{sl} is the equilibrium concentration of unionised dissolved gas in liquid phase [g/cm^3of water]), and K_g and K_l are the gas-phase and liquid-phase exchange constants, respectively. 1/K_g = 1/k_g + H/k_l 1/K_l = 1/k_l + 1/(H*k_g).  
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| t | intent(inout) | calc_ka | is the temperature at surface in [C] |
| p | intent(inout) | calc_ka | is the pressure at surface in [Pa] |
| mw | intent(inout) | calc_ka | is the molecular weight [g/mol] |
| vb | intent(inout) | calc_ka | is the molar volume [m^3/mol] |
| u10 | intent(inout) | calc_ka | is the wind speed at 10m above the surface in [m/s] |
| ustar | intent(inout) | calc_ka | is the Friction velocity [m/s]. If not provided, ustar = u_{10}*sqrt{C_D}. |
| cd_m | intent(inout) | calc_ka | is the Drag coefficient C_D. Used only if ustar is provided. If ustar is not provided, cd_m = 6.1x10^{-4} + 0.63x10^{-4} * u_10 |

### flowchart
calc_ka does the following:  


## calc_kl
### intro
Module atmos_ocean_fluxes_calc_mod calculates gas fluxes between atmosphere and ocean.  calc_kl is a real function in atmos_ocean_fluxes_calc_mod.
### description
  Function calc_kl computes k_l, the liquid-side transfer velocity. See Johnson, Ocean Science, 2010. (http://doi.org/10.5194/os-6-913-2010) and Nightingale, Global Biogeochemical Cycles, 2000 (https://doi.org/10.1029/1999GB900091).
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| t | intent(inout) | calc_kl | is the temperature at surface [C] |
| v | intent(inout) | calc_kl | is the wind speed at surface [m/s] |
| sc | intent(inout) | calc_kl | is the Schmidt number [dimensionless] |

### flowchart
calc_kl does the following:  


## schmidt_g
### intro
Module atmos_ocean_fluxes_calc_mod calculates gas fluxes between atmosphere and ocean.  schmidt_g is a real function in atmos_ocean_fluxes_calc_mod.
### description
  Function schmidt_g computes the Schmidt number of gas in air.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| t | intent(inout) | schmidt_g | is the temperature at surface [C] |
| p | intent(inout) | schmidt_g | is the pressure at surface [Pa] |
| mw | intent(inout) | schmidt_g | is the molecular weight [g/mol] |
| vb | intent(inout) | schmidt_g | is the molar volume [cm^3/mol] |

### flowchart
schmidt_g does the following:  


## d_air
### intro
Module atmos_ocean_fluxes_calc_mod calculates gas fluxes between atmosphere and ocean.  d_air is a real function in atmos_ocean_fluxes_calc_mod.
### description
  Function d_air computes the diffusion coefficient of gas in air [m^2/s] following Fuller, Industrial & Engineering Chemistry (https://doi.org/10.1021/ie50677a007).
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| t | intent(inout) | d_air | is the temperature [C] |
| p | intent(inout) | d_air | is the pressure [Pa] |
| mw | intent(inout) | d_air | is the molecular weight [g/mol] |
| vb | intent(inout) | d_air | is the diffusion coefficient [cm^3/mol] |

### flowchart
d_air does the following:  


## p_air
### intro
Module atmos_ocean_fluxes_calc_mod calculates gas fluxes between atmosphere and ocean.  p_air is a real function in atmos_ocean_fluxes_calc_mod.
### description
  Function p_air computes the density of air [kg/m^3] as a cubic polynomial in temperature. Coefficients sd_0, ..., sd_3 approximate the dry-air density at standard pressure.

t


is the temperature at surface [C].
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| t | intent(inout) | p_air | is the temperature at surface [C] |

### flowchart
p_air does the following:  


## v_air
### intro
Module atmos_ocean_fluxes_calc_mod calculates gas fluxes between atmosphere and ocean.  v_air is a real function in atmos_ocean_fluxes_calc_mod.
### description
  Function v_air computes the kinematic viscosity in air [m^2/s].
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| t | intent(inout) | v_air | is the temperature at surface [C] |

### flowchart
v_air does the following:  


## n_air
### intro
Module atmos_ocean_fluxes_calc_mod calculates gas fluxes between atmosphere and ocean.  n_air is a real function in atmos_ocean_fluxes_calc_mod.
### description
  Function n_air computes the dynamic viscosity in air [Pa*s].
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| t | intent(inout) | n_air | is the temperature at surface [C] |

### flowchart
n_air does the following:  

