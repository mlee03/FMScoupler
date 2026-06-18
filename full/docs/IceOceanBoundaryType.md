# `ice_ocean_boundary_type`

## Ice_ocean_boundary_type
`ice_ocean_boundary_type` holds all surface forcing passed from the sea-ice model to MOM6 each coupled timestep.

## Ice_ocean_boundary_type%u_flux
Ice_ocean_boundary_type%u_flux, a real 2D array, is i-direction wind/ice stress on the ocean surface [Pa].
## Ice_ocean_boundary_type%v_flux
Ice_ocean_boundary_type%v_flux, a real 2D array, is j-direction wind/ice stress on the ocean surface [Pa].
## Ice_ocean_boundary_type%stress_mag
Ice_ocean_boundary_type%stress_mag, a real 2D array, is Time-mean magnitude of the stress on the ocean [Pa]; present when pass_stress_mag=.true. in SIS_slow_CS. FMS cap only..
## Ice_ocean_boundary_type%wind_stagger
Ice_ocean_boundary_type%wind_stagger, integer, is Spatial discretization of the wind stresses; may be set by the flux-exchange code based on what the sea-ice model provides, otherwise taken from the surface forcing control structure.
## Ice_ocean_boundary_type%u10_sqr
Ice_ocean_boundary_type%u10_sqr, a real 2D array, is Wind speed squared at 10 m height [m**2/s**2]. NUOPC cap only..

## Ice_ocean_boundary_type%t_flux
Ice_ocean_boundary_type%t_flux, a real 2D array, is Sensible heat flux into the ocean W/m**2.
## Ice_ocean_boundary_type%lw_flux
Ice_ocean_boundary_type%lw_flux, a real 2D array, is Net longwave radiation flux into the ocean W/m**2.
## Ice_ocean_boundary_type%sw_flux_vis_dir
Ice_ocean_boundary_type%sw_flux_vis_dir, a real 2D array, is Direct visible shortwave radiation into theIceOceanBoundaryType.md ocean W/m**2.
## Ice_ocean_boundary_type%sw_flux_vis_dif
Ice_ocean_boundary_type%sw_flux_vis_dif, a real 2D array, is Diffuse visible shortwave radiation into the ocean W/m**2.
## Ice_ocean_boundary_type%sw_flux_nir_dir
Ice_ocean_boundary_type%sw_flux_nir_dir, a real 2D array, is Direct near-infrared shortwave radiation into the ocean W/m**2.
## Ice_ocean_boundary_type%sw_flux_nir_dif
Ice_ocean_boundary_type%sw_flux_nir_dif, a real 2D array, is Diffuse near-infrared shortwave radiation into the ocean W/m**2.
## Ice_ocean_boundary_type%seaice_melt_heat
Ice_ocean_boundary_type%seaice_melt_heat, a real 2D array, is Heat flux from sea ice and snow melting W/m**2. NUOPC cap only..
## Ice_ocean_boundary_type%swnet_afracr
Ice_ocean_boundary_type%swnet_afracr, a real 2D array, is Net shortwave radiation multiplied by the atmosphere fraction, positive into the ocean W/m**2. NUOPC cap only..
## Ice_ocean_boundary_type%swpen_ifrac_n
Ice_ocean_boundary_type%swpen_ifrac_n, a real 3D array, is Net shortwave radiation penetrating into ice and ocean, multiplied by ice fraction per thickness category; positive into the ocean W/m**2; third dimension indexes ice categories. NUOPC cap only..

## Ice_ocean_boundary_type%hrofl
Ice_ocean_boundary_type%hrofl, a real 2D array, is Heat content from liquid runoff W/m**2.
## Ice_ocean_boundary_type%hrofi
Ice_ocean_boundary_type%hrofi, a real 2D array, is Heat content from frozen runoff (calving) W/m**2.
## Ice_ocean_boundary_type%hrofl_glc
Ice_ocean_boundary_type%hrofl_glc, a real 2D array, is Heat content from liquid glacier runoff via the river-routing model W/m**2.
## Ice_ocean_boundary_type%hrofi_glc
Ice_ocean_boundary_type%hrofi_glc, a real 2D array, is Heat content from frozen glacier runoff via the river-routing model W/m**2.
## Ice_ocean_boundary_type%hrain
Ice_ocean_boundary_type%hrain, a real 2D array, is Heat content from liquid precipitation W/m**2.
## Ice_ocean_boundary_type%hsnow
Ice_ocean_boundary_type%hsnow, a real 2D array, is Heat content from frozen precipitation W/m**2.
## Ice_ocean_boundary_type%hevap
Ice_ocean_boundary_type%hevap, a real 2D array, is Heat content from evaporation W/m**2.
## Ice_ocean_boundary_type%hcond
Ice_ocean_boundary_type%hcond, a real 2D array, is Heat content from condensation W/m**2.

## Ice_ocean_boundary_type%q_flux
Ice_ocean_boundary_type%q_flux, a real 2D array, is Specific humidity (freshwater) flux into the ocean [kg/m**2/s].
## Ice_ocean_boundary_type%salt_flux
Ice_ocean_boundary_type%salt_flux, a real 2D array, is Salt flux from sea ice into the ocean (brine rejection / melting) [kg/m**2/s].
## Ice_ocean_boundary_type%excess_salt
Ice_ocean_boundary_type%excess_salt, a real 2D array, is Salt left behind in the ocean by brine rejection rather than ejected as a salt flux [kg/m**2/s]. FMS cap only..
## Ice_ocean_boundary_type%seaice_melt
Ice_ocean_boundary_type%seaice_melt, a real 2D array, is Water flux due to sea ice and snow melting [kg/m**2/s]. NUOPC cap only..
## Ice_ocean_boundary_type%lprec
Ice_ocean_boundary_type%lprec, a real 2D array, is Mass flux of liquid precipitation into the ocean [kg/m**2/s].
## Ice_ocean_boundary_type%fprec
Ice_ocean_boundary_type%fprec, a real 2D array, is Mass flux of frozen precipitation into the ocean [kg/m**2/s].

## Ice_ocean_boundary_type%runoff
Ice_ocean_boundary_type%runoff, a real 2D array, is Mass flux of liquid runoff from land into the ocean [kg/m**2/s]. (FMS cap only)
## Ice_ocean_boundary_type%runoff_carbon
Ice_ocean_boundary_type%runoff_carbon, a real 2D array, is Mass flux of carbon carried by liquid runoff [kg/m**2/s]. (FMS cap only)
## Ice_ocean_boundary_type%runoff_hflx
Ice_ocean_boundary_type%runoff_hflx, a real 2D array, is Heat content of liquid runoff relative to 0 °C W/m**2. (FMS cap only)
## Ice_ocean_boundary_type%calving
Ice_ocean_boundary_type%calving, a real 2D array, is Mass flux of frozen runoff (calving) into the ocean [kg/m**2/s]; offered first to icebergs if active. (FMS cap only)
## Ice_ocean_boundary_type%calving_hflx
Ice_ocean_boundary_type%calving_hflx, a real 2D array, is Heat content of frozen runoff relative to 0 °C W/m**2. (FMS cap only)
## Ice_ocean_boundary_type%lrunoff
Ice_ocean_boundary_type%lrunoff, a real 2D array, is Liquid runoff [kg/m**2/s]. (NUOPC cap only)
## Ice_ocean_boundary_type%frunoff
Ice_ocean_boundary_type%frunoff, a real 2D array, is Frozen (ice) runoff [kg/m**2/s]. (NUOPC cap only)
## Ice_ocean_boundary_type%lrunoff_glc
Ice_ocean_boundary_type%lrunoff_glc, a real 2D array, is Liquid glacier runoff delivered via the river-routing model [kg/m**2/s]. (NUOPC cap only)
## Ice_ocean_boundary_type%frunoff_glc
Ice_ocean_boundary_type%frunoff_glc, a real 2D array, is Frozen glacier runoff delivered via the river-routing model [kg/m**2/s]. (NUOPC cap only)

## Ice_ocean_boundary_type%p
Ice_ocean_boundary_type%p, a real 2D array, is Pressure of overlying ice and atmosphere on the ocean surface [Pa].
## Ice_ocean_boundary_type%mi
Ice_ocean_boundary_type%mi, a real 2D array, is Mass of sea ice per unit ocean area [kg/m**2]; used for ice-pressure loading.
## Ice_ocean_boundary_type%ice_rigidity
Ice_ocean_boundary_type%ice_rigidity, a real 2D array, is Rigidity of sea ice and ice shelves expressed as a divergence-damping coefficient [m³/s]; determined outside the ocean model.
## Ice_ocean_boundary_type%ice_fraction
Ice_ocean_boundary_type%ice_fraction, a real 2D array, is Fractional ice area [dimensionless]. NUOPC cap only..
## Ice_ocean_boundary_type%ifrac_n
Ice_ocean_boundary_type%ifrac_n, a real 3D array, is Ice fraction per ice thickness category [dimensionless]; third dimension indexes categories. NUOPC cap only..
## Ice_ocean_boundary_type%ice_ncat
Ice_ocean_boundary_type%ice_ncat, integer, is Number of ice categories provided by the coupler; 1 means per-category data is not used. NUOPC cap only..
## Ice_ocean_boundary_type%afracr
Ice_ocean_boundary_type%afracr, a real 2D array, is Fractional atmosphere coverage relative to the ocean grid cell [dimensionless]. NUOPC cap only..

## Ice_ocean_boundary_type%ustar_berg
Ice_ocean_boundary_type%ustar_berg, a real 2D array, is Frictional velocity beneath icebergs [m/s].
## Ice_ocean_boundary_type%area_berg
Ice_ocean_boundary_type%area_berg, a real 2D array, is Fractional area of the ocean cell covered by icebergs [m**2/m**2].
## Ice_ocean_boundary_type%mass_berg
Ice_ocean_boundary_type%mass_berg, a real 2D array, is Mass of icebergs per unit ocean area [kg/m**2].

## Ice_ocean_boundary_type%shelf_sfc_mass_flux
Ice_ocean_boundary_type%shelf_sfc_mass_flux, a real 2D array, is Mass flux to the surface of the ice sheet [kg/m**2/s].

## Ice_ocean_boundary_type%nhx_dep
Ice_ocean_boundary_type%nhx_dep, a real 2D array, is Reduced nitrogen (NHx) deposition flux [kg/m**2/s].
This field support ocean biogeochemistry modules that require atmospheric deposition forcing.
## Ice_ocean_boundary_type%noy_dep
Ice_ocean_boundary_type%noy_dep, a real 2D array, is Oxidized nitrogen (NOy) deposition flux [kg/m**2/s].
This field support ocean biogeochemistry modules that require atmospheric deposition forcing.

## Ice_ocean_boundary_type%atm_co2_prog
Ice_ocean_boundary_type%atm_co2_prog, a real 2D array, is Prognostic atmospheric CO₂ concentration [ppm].
This field support ocean biogeochemistry modules that require atmospheric deposition forcing.

## Ice_ocean_boundary_type%atm_co2_diag
Ice_ocean_boundary_type%atm_co2_diag, a real 2D array, is Diagnostic atmospheric CO₂ concentration [ppm].
This field support ocean biogeochemistry modules that require atmospheric deposition forcing.

## Ice_ocean_boundary_type%atm_fine_dust_flux
Ice_ocean_boundary_type%atm_fine_dust_flux, a real 2D array, is Fine dust deposition flux from the atmosphere [kg/m**2/s].
This field support ocean biogeochemistry modules that require atmospheric deposition forcing.

## Ice_ocean_boundary_type%atm_coarse_dust_flux
Ice_ocean_boundary_type%atm_coarse_dust_flux, a real 2D array, is Coarse dust deposition flux from the atmosphere [kg/m**2/s].
This field support ocean biogeochemistry modules that require atmospheric deposition forcing.

## Ice_ocean_boundary_type%seaice_dust_flux
Ice_ocean_boundary_type%seaice_dust_flux, a real 2D array, is Dust flux released from sea ice [kg/m**2/s].
This field support ocean biogeochemistry modules that require atmospheric deposition forcing.

## Ice_ocean_boundary_type%atm_bc_flux
Ice_ocean_boundary_type%atm_bc_flux, a real 2D array, is Black carbon deposition flux from the atmosphere [kg/m**2/s].
This field support ocean biogeochemistry modules that require atmospheric deposition forcing.

## Ice_ocean_boundary_type%seaice_bc_flux
Ice_ocean_boundary_type%seaice_bc_flux, a real 2D array, is Black carbon flux released from sea ice [kg/m**2/s].
This field support ocean biogeochemistry modules that require atmospheric deposition forcing.

## Ice_ocean_boundary_type%lamult
Ice_ocean_boundary_type%lamult, a real 2D array, is Langmuir turbulence enhancement factor [dimensionless].
This field support Langmuir turbulence and wave-driven mixing parameterizations in MOM6.
## Ice_ocean_boundary_type%stk_wavenumbers
Ice_ocean_boundary_type%stk_wavenumbers, a real 1D array, is Central wavenumber of each Stokes drift band [rad/m]; dimensioned (num_stk_bands).
This field support Langmuir turbulence and wave-driven mixing parameterizations in MOM6.
## Ice_ocean_boundary_type%ustkb
Ice_ocean_boundary_type%ustkb, a real 3D array, is Stokes drift spectrum, zonal component, at u-points [m/s]; third dimension indexes wavenumber bands.
This field support Langmuir turbulence and wave-driven mixing parameterizations in MOM6.
## Ice_ocean_boundary_type%vstkb
Ice_ocean_boundary_type%vstkb, a real 3D array, is Stokes drift spectrum, meridional component, at v-points [m/s]; third dimension indexes wavenumber bands.
This field support Langmuir turbulence and wave-driven mixing parameterizations in MOM6.
## Ice_ocean_boundary_type%num_stk_bands
Ice_ocean_boundary_type%num_stk_bands, integer, is Number of Stokes drift wavenumber bands passed through the coupler.
This field support Langmuir turbulence and wave-driven mixing parameterizations in MOM6.

## Ice_ocean_boundary_type%xtype
Ice_ocean_boundary_type%xtype, integer, is Transfer mode for the ice-to-ocean exchange: REGRID (1), REDIST (2), or DIRECT (3).
## Ice_ocean_boundary_type%fluxes
Ice_ocean_boundary_type%fluxes, type(coupler_2d_bc_type), is Named array of additional per-tracer passive tracer fluxes from ice/atmosphere to ocean.
