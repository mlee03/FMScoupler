# `ice_ocean_boundary_type` — Ice-to-Ocean Boundary Fields

## Overview

`ice_ocean_boundary_type` holds all surface forcing passed from the sea-ice model (SIS2) to the ocean model (MOM6) each coupled timestep. An instance named `Ice_ocean_boundary` is declared in `coupler_main.F90`. Fields are 2D arrays unless otherwise noted.

**Populated by:** `flux_ice_to_ocean`  
**Consumed by:** `update_ocean_model`  
**Related types:** `ocean_ice_boundary_type`, `ice_data_type`

> **Cap-specific fields:** Some fields are only present in the FMS coupler cap; others only in the NUOPC cap. These are labelled `(FMS cap only)` or `(NUOPC cap only)` respectively.

---

## Wind Stress Fields

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Ice_ocean_boundary%u_flux` | real 2D | Pa | i-direction wind/ice stress on the ocean surface. |
| `Ice_ocean_boundary%v_flux` | real 2D | Pa | j-direction wind/ice stress on the ocean surface. |
| `Ice_ocean_boundary%stress_mag` | real 2D | Pa | Time-mean magnitude of the stress on the ocean; present when `pass_stress_mag=.true.` in `SIS_slow_CS`. (FMS cap only) |
| `Ice_ocean_boundary%wind_stagger` | integer | — | Spatial discretization of the wind stresses; may be set by the flux-exchange code based on what the sea-ice model provides, otherwise taken from the surface forcing control structure. |
| `Ice_ocean_boundary%u10_sqr` | real 2D | m²/s² | Wind speed squared at 10 m height. (NUOPC cap only) |

---

## Heat Flux Fields

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Ice_ocean_boundary%t_flux` | real 2D | W/m² | Sensible heat flux into the ocean. |
| `Ice_ocean_boundary%lw_flux` | real 2D | W/m² | Net longwave radiation flux into the ocean. |
| `Ice_ocean_boundary%sw_flux_vis_dir` | real 2D | W/m² | Direct visible shortwave radiation into the ocean. |
| `Ice_ocean_boundary%sw_flux_vis_dif` | real 2D | W/m² | Diffuse visible shortwave radiation into the ocean. |
| `Ice_ocean_boundary%sw_flux_nir_dir` | real 2D | W/m² | Direct near-infrared shortwave radiation into the ocean. |
| `Ice_ocean_boundary%sw_flux_nir_dif` | real 2D | W/m² | Diffuse near-infrared shortwave radiation into the ocean. |
| `Ice_ocean_boundary%seaice_melt_heat` | real 2D | W/m² | Heat flux from sea ice and snow melting. (NUOPC cap only) |
| `Ice_ocean_boundary%swnet_afracr` | real 2D | W/m² | Net shortwave radiation multiplied by the atmosphere fraction, positive into the ocean. (NUOPC cap only) |
| `Ice_ocean_boundary%swpen_ifrac_n` | real 3D | W/m² | Net shortwave radiation penetrating into ice and ocean, multiplied by ice fraction per thickness category; third dimension indexes ice categories. (NUOPC cap only) |
| `Ice_ocean_boundary%hrofl` | real 2D | W/m² | Heat content from liquid runoff. |
| `Ice_ocean_boundary%hrofi` | real 2D | W/m² | Heat content from frozen runoff (calving). |
| `Ice_ocean_boundary%hrofl_glc` | real 2D | W/m² | Heat content from liquid glacier runoff via the river-routing model. |
| `Ice_ocean_boundary%hrofi_glc` | real 2D | W/m² | Heat content from frozen glacier runoff via the river-routing model. |
| `Ice_ocean_boundary%hrain` | real 2D | W/m² | Heat content from liquid precipitation. |
| `Ice_ocean_boundary%hsnow` | real 2D | W/m² | Heat content from frozen precipitation. |
| `Ice_ocean_boundary%hevap` | real 2D | W/m² | Heat content from evaporation. |
| `Ice_ocean_boundary%hcond` | real 2D | W/m² | Heat content from condensation. |

---

## Freshwater and Salt Flux Fields

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Ice_ocean_boundary%q_flux` | real 2D | kg/m²/s | Specific humidity (freshwater) flux into the ocean. |
| `Ice_ocean_boundary%salt_flux` | real 2D | kg/m²/s | Salt flux from sea ice into the ocean (brine rejection / melting). |
| `Ice_ocean_boundary%excess_salt` | real 2D | kg/m²/s | Salt left behind in the ocean by brine rejection rather than ejected as a salt flux. (FMS cap only) |
| `Ice_ocean_boundary%seaice_melt` | real 2D | kg/m²/s | Water flux due to sea ice and snow melting. (NUOPC cap only) |
| `Ice_ocean_boundary%lprec` | real 2D | kg/m²/s | Mass flux of liquid precipitation into the ocean. |
| `Ice_ocean_boundary%fprec` | real 2D | kg/m²/s | Mass flux of frozen precipitation into the ocean. |

---

## Land Runoff and Calving Fields

### FMS cap

| Field | Units | Description |
|---|---|---|
| `Ice_ocean_boundary%runoff` | kg/m²/s | Mass flux of liquid runoff from land into the ocean. |
| `Ice_ocean_boundary%runoff_carbon` | kg/m²/s | Mass flux of carbon carried by liquid runoff. |
| `Ice_ocean_boundary%runoff_hflx` | W/m² | Heat content of liquid runoff relative to 0 °C. |
| `Ice_ocean_boundary%calving` | kg/m²/s | Mass flux of frozen runoff (calving) into the ocean; offered first to icebergs if active. |
| `Ice_ocean_boundary%calving_hflx` | W/m² | Heat content of frozen runoff relative to 0 °C. |

### NUOPC cap

| Field | Units | Description |
|---|---|---|
| `Ice_ocean_boundary%lrunoff` | kg/m²/s | Liquid runoff. |
| `Ice_ocean_boundary%frunoff` | kg/m²/s | Frozen (ice) runoff. |
| `Ice_ocean_boundary%lrunoff_glc` | kg/m²/s | Liquid glacier runoff delivered via the river-routing model. |
| `Ice_ocean_boundary%frunoff_glc` | kg/m²/s | Frozen glacier runoff delivered via the river-routing model. |

---

## Pressure, Mass Loading, and Sea-Ice State

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Ice_ocean_boundary%p` | real 2D | Pa | Pressure of overlying ice and atmosphere on the ocean surface. |
| `Ice_ocean_boundary%mi` | real 2D | kg/m² | Mass of sea ice per unit ocean area; used for ice-pressure loading. |
| `Ice_ocean_boundary%ice_rigidity` | real 2D | m³/s | Rigidity of sea ice and ice shelves expressed as a divergence-damping coefficient; determined outside the ocean model. |
| `Ice_ocean_boundary%ice_fraction` | real 2D | dimensionless | Fractional ice area. (NUOPC cap only) |
| `Ice_ocean_boundary%ifrac_n` | real 3D | dimensionless | Ice fraction per ice thickness category; third dimension indexes categories. (NUOPC cap only) |
| `Ice_ocean_boundary%ice_ncat` | integer | — | Number of ice categories provided by the coupler; 1 means per-category data is not used. (NUOPC cap only) |
| `Ice_ocean_boundary%afracr` | real 2D | dimensionless | Fractional atmosphere coverage relative to the ocean grid cell. (NUOPC cap only) |
| `Ice_ocean_boundary%shelf_sfc_mass_flux` | real 2D | kg/m²/s | Mass flux to the surface of the ice sheet. |

---

## Iceberg Fields

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Ice_ocean_boundary%ustar_berg` | real 2D | m/s | Frictional velocity beneath icebergs. |
| `Ice_ocean_boundary%area_berg` | real 2D | m²/m² | Fractional area of the ocean cell covered by icebergs. |
| `Ice_ocean_boundary%mass_berg` | real 2D | kg/m² | Mass of icebergs per unit ocean area. |

---

## Biogeochemistry Deposition Fields

These fields support ocean biogeochemistry modules that require atmospheric deposition forcing.

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Ice_ocean_boundary%nhx_dep` | real 2D | kg/m²/s | Reduced nitrogen (NHx) deposition flux. |
| `Ice_ocean_boundary%noy_dep` | real 2D | kg/m²/s | Oxidized nitrogen (NOy) deposition flux. |
| `Ice_ocean_boundary%atm_co2_prog` | real 2D | ppm | Prognostic atmospheric CO₂ concentration. |
| `Ice_ocean_boundary%atm_co2_diag` | real 2D | ppm | Diagnostic atmospheric CO₂ concentration. |
| `Ice_ocean_boundary%atm_fine_dust_flux` | real 2D | kg/m²/s | Fine dust deposition flux from the atmosphere. |
| `Ice_ocean_boundary%atm_coarse_dust_flux` | real 2D | kg/m²/s | Coarse dust deposition flux from the atmosphere. |
| `Ice_ocean_boundary%seaice_dust_flux` | real 2D | kg/m²/s | Dust flux released from sea ice. |
| `Ice_ocean_boundary%atm_bc_flux` | real 2D | kg/m²/s | Black carbon deposition flux from the atmosphere. |
| `Ice_ocean_boundary%seaice_bc_flux` | real 2D | kg/m²/s | Black carbon flux released from sea ice. |

---

## Langmuir Turbulence and Wave Fields

These fields support Langmuir turbulence and wave-driven mixing parameterizations in MOM6.

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Ice_ocean_boundary%lamult` | real 2D | dimensionless | Langmuir turbulence enhancement factor. |
| `Ice_ocean_boundary%stk_wavenumbers` | real 1D | rad/m | Central wavenumber of each Stokes drift band; dimensioned `(num_stk_bands)`. |
| `Ice_ocean_boundary%ustkb` | real 3D | m/s | Stokes drift spectrum, zonal component, at u-points; third dimension indexes wavenumber bands. |
| `Ice_ocean_boundary%vstkb` | real 3D | m/s | Stokes drift spectrum, meridional component, at v-points; third dimension indexes wavenumber bands. |
| `Ice_ocean_boundary%num_stk_bands` | integer | — | Number of Stokes drift wavenumber bands passed through the coupler. |

---

## Transfer and Tracer Metadata

| Field | Type | Description |
|---|---|---|
| `Ice_ocean_boundary%xtype` | integer | Transfer mode for the ice-to-ocean exchange: `REGRID` (1), `REDIST` (2), or `DIRECT` (3). |
| `Ice_ocean_boundary%fluxes` | `type(coupler_2d_bc_type)` | Named array of additional per-tracer passive tracer fluxes from ice/atmosphere to ocean. |
