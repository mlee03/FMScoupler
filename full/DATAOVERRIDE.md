# Data Override Variables in FMSCoupler/full

This document lists all variables passed to `fms_data_override` (or
`fms_coupler_type_data_override`) in the top-level `.F90` files of
`FMSCoupler/full`, organized by the subroutine in which the call appears.
A data override replaces the model-computed value with externally supplied
data only when a matching entry exists in the `data_table`; otherwise the
model value is left unchanged.

---

## `atm_land_ice_flux_exchange.F90`

### `sfc_boundary_layer`

This subroutine computes turbulent surface fluxes between the atmosphere,
land, and sea ice. Data overrides here allow replacement of atmospheric
boundary-layer state and surface properties before flux calculations.

#### Component: ATM

| Fortran field | Description |
|---|---|
| `Atm%t_bot` | Temperature at the lowest atmospheric level [K] |
| `Atm%z_bot` | Height of the lowest atmospheric level [m] |
| `Atm%p_bot` | Pressure at the lowest atmospheric level [Pa] |
| `Atm%u_bot` | Zonal wind at the lowest atmospheric level [m/s] |
| `Atm%v_bot` | Meridional wind at the lowest atmospheric level [m/s] |
| `Atm%p_surf` | Surface pressure [Pa] |
| `Atm%slp` | Sea-level pressure [Pa] |
| `Atm%gust` | Gustiness velocity used to augment surface wind speed in flux calculations [m/s] |
| `atm%fields%bc(n)%field(m)%values` | Per-tracer atmospheric surface fields (e.g. tracer concentrations at the lowest model level); iterated over all boundary-condition fields in `Atm%fields` |

#### Component: ATM — `Land_Ice_Atmos_Boundary` fields
*(Fields returned from the surface back to the atmosphere)*

| Fortran field | Description |
|---|---|
| `Land_Ice_Atmos_Boundary%t` | Surface temperature (area-weighted over land and ice fractions) seen by the atmosphere [K] |
| `Land_Ice_Atmos_Boundary%albedo` | Broadband surface albedo [dimensionless] |
| `Land_Ice_Atmos_Boundary%albedo_vis_dir` | Direct-beam visible-band surface albedo [dimensionless] |
| `Land_Ice_Atmos_Boundary%albedo_nir_dir` | Direct-beam near-infrared surface albedo [dimensionless] |
| `Land_Ice_Atmos_Boundary%albedo_vis_dif` | Diffuse visible-band surface albedo [dimensionless] |
| `Land_Ice_Atmos_Boundary%albedo_nir_dif` | Diffuse near-infrared surface albedo [dimensionless] |
| `Land_Ice_Atmos_Boundary%land_frac` | Fraction of atmospheric grid cell covered by land [dimensionless] |
| `Land_Ice_Atmos_Boundary%dt_t` | Implicit correction to atmospheric temperature from surface flux scheme [K] |
| `Land_Ice_Atmos_Boundary%dt_tr(:,:,tr)` | Implicit correction to each atmospheric tracer mixing ratio from surface flux scheme; one entry per exchanged tracer |
| `Land_Ice_Atmos_Boundary%u_flux` | Zonal wind stress on the atmosphere [Pa] |
| `Land_Ice_Atmos_Boundary%v_flux` | Meridional wind stress on the atmosphere [Pa] |
| `Land_Ice_Atmos_Boundary%dtaudu` | d(wind stress)/d(u) — implicit coupling coefficient for zonal momentum [Pa·s/m] |
| `Land_Ice_Atmos_Boundary%dtaudv` | d(wind stress)/d(v) — implicit coupling coefficient for meridional momentum [Pa·s/m] |
| `Land_Ice_Atmos_Boundary%u_star` | Friction velocity (surface turbulent velocity scale) [m/s] |
| `Land_Ice_Atmos_Boundary%b_star` | Buoyancy scale used in Monin-Obukhov similarity theory [m/s²] |
| `Land_Ice_Atmos_Boundary%rough_mom` | Roughness length for momentum [m] |

#### Component: ICE

| Fortran field | Description |
|---|---|
| `Ice%t_surf` | Ice/ocean surface skin temperature [K] |
| `Ice%rough_mom` | Surface roughness length for momentum over ice [m] |
| `Ice%rough_heat` | Surface roughness length for heat over ice [m] |
| `Ice%rough_moist` | Surface roughness length for moisture over ice [m] |
| `Ice%albedo` | Broadband surface albedo over ice [dimensionless] |
| `Ice%albedo_vis_dir` | Direct-beam visible-band albedo over ice [dimensionless] |
| `Ice%albedo_nir_dir` | Direct-beam near-infrared albedo over ice [dimensionless] |
| `Ice%albedo_vis_dif` | Diffuse visible-band albedo over ice [dimensionless] |
| `Ice%albedo_nir_dif` | Diffuse near-infrared albedo over ice [dimensionless] |
| `Ice%u_surf` | Zonal surface current velocity of ice/ocean [m/s] |
| `Ice%v_surf` | Meridional surface current velocity of ice/ocean [m/s] |
| `Ice%ocean_fields` | Coupler boundary-condition type holding ocean/ice-top gas and tracer fields used in atmosphere-ocean flux calculations |

#### Component: LND

| Fortran field | Description |
|---|---|
| `Land%t_surf` | Land surface (radiative) temperature [K] |
| `Land%t_ca` | Canopy air temperature — near-surface air temperature within the plant canopy [K] |
| `Land%rough_mom` | Surface roughness length for momentum over land [m] |
| `Land%rough_heat` | Surface roughness length for heat over land [m] |
| `Land%albedo` | Broadband surface albedo over land [dimensionless] |
| `Land%tr(:,:,:,tr)` | Surface tracer mixing ratio over land; one entry per exchanged land tracer |
| `Land%albedo_vis_dir` | Direct-beam visible-band albedo over land [dimensionless] |
| `Land%albedo_nir_dir` | Direct-beam near-infrared albedo over land [dimensionless] |
| `Land%albedo_vis_dif` | Diffuse visible-band albedo over land [dimensionless] |
| `Land%albedo_nir_dif` | Diffuse near-infrared albedo over land [dimensionless] |

---

### `flux_down_from_atmos`

This subroutine transfers downward atmospheric forcing to land and ice.
Data overrides here replace computed atmospheric fluxes and implicit
diffusion coefficients before they are handed off to the surface models.

#### Component: ATM — downward forcing fields

| Fortran field | Description |
|---|---|
| `Atm%flux_sw` | Total net shortwave flux at the surface [W/m²] |
| `Atm%flux_sw_dir` | Direct-beam component of net shortwave flux [W/m²] |
| `Atm%flux_sw_dif` | Diffuse component of net shortwave flux [W/m²] |
| `Atm%flux_sw_down_vis_dir` | Downward direct-beam visible shortwave flux [W/m²] |
| `Atm%flux_sw_down_vis_dif` | Downward diffuse visible shortwave flux [W/m²] |
| `Atm%flux_sw_down_total_dir` | Downward direct-beam total (broadband) shortwave flux [W/m²] |
| `Atm%flux_sw_down_total_dif` | Downward diffuse total (broadband) shortwave flux [W/m²] |
| `Atm%flux_sw_vis` | Net visible-band shortwave flux at the surface [W/m²] |
| `Atm%flux_sw_vis_dir` | Direct-beam component of net visible shortwave flux [W/m²] |
| `Atm%flux_sw_vis_dif` | Diffuse component of net visible shortwave flux [W/m²] |
| `Atm%flux_lw` | Net downward longwave flux at the surface [W/m²] |
| `Atm%lprec` | Liquid precipitation rate [kg/m²/s] |
| `frac_precip` | 2-D scaling field for liquid precipitation; applied when `scale_precip_2d=.true.` [dimensionless] |
| `Atm%fprec` | Frozen (solid) precipitation rate [kg/m²/s] |
| `Atm%coszen` | Cosine of the solar zenith angle [dimensionless] |
| `Atm%Surf_Diff%dtmass` | dt/mass — ratio of timestep to surface layer mass used in the implicit diffusion scheme [s·m²/kg] |
| `Atm%Surf_Diff%delta_t` | Forward-elimination temperature coefficient from the implicit vertical diffusion scheme [K] |
| `Atm%Surf_Diff%dflux_t` | d(sensible heat flux)/d(T_surf) — linearisation coefficient for the implicit heat flux scheme [W/m²/K] |
| `Atm%Surf_Diff%delta_tr(:,:,tr)` | Forward-elimination coefficient for each tracer from the implicit diffusion scheme |
| `Atm%Surf_Diff%dflux_tr(:,:,tr)` | d(tracer flux)/d(tracer_surf) — linearisation coefficient for the implicit tracer flux scheme |

#### Component: LND — fields sent to land

| Fortran field | Description |
|---|---|
| `Land_boundary%drag_q` | Drag coefficient for moisture used in land surface flux calculations [dimensionless] |
| `Land_boundary%lwdn_flux` | Downward longwave radiation flux at the land surface [W/m²] |
| `Land_boundary%cd_m` | Drag coefficient for momentum over land [dimensionless] |
| `Land_boundary%cd_t` | Drag coefficient for heat over land [dimensionless] |
| `Land_boundary%bstar` | Buoyancy scale passed to the land model [m/s²] |
| `Land_boundary%ustar` | Friction velocity passed to the land model [m/s] |
| `Land_boundary%wind` | Wind speed at the lowest atmospheric level for land surface calculations [m/s] |
| `Land_boundary%z_bot` | Height of the lowest atmospheric level above the land surface [m] |
| `Land_boundary%t_flux` | Sensible heat flux into the land surface [W/m²] |
| `Land_boundary%lw_flux` | Net longwave flux at the land surface [W/m²] |
| `Land_boundary%sw_flux` | Net shortwave flux at the land surface [W/m²] |
| `Land_boundary%sw_flux_down_vis_dir` | Downward direct-beam visible shortwave flux over land [W/m²] |
| `Land_boundary%sw_flux_down_total_dir` | Downward direct-beam total shortwave flux over land [W/m²] |
| `Land_boundary%sw_flux_down_vis_dif` | Downward diffuse visible shortwave flux over land [W/m²] |
| `Land_boundary%sw_flux_down_total_dif` | Downward diffuse total shortwave flux over land [W/m²] |
| `Land_boundary%lprec` | Liquid precipitation rate over land [kg/m²/s] |
| `Land_boundary%fprec` | Frozen precipitation rate over land [kg/m²/s] |
| `Land_boundary%dhdt` | d(sensible heat flux)/d(T_surf) — implicit coupling coefficient for land heat flux [W/m²/K] |
| `Land_boundary%drdt` | d(longwave flux)/d(T_surf) — implicit coupling coefficient for land longwave flux [W/m²/K] |
| `Land_boundary%p_surf` | Surface pressure over land [Pa] |
| `Land_boundary%tr_flux(:,:,:,tr)` | Flux of each exchanged tracer into the land surface [kg/m²/s] |
| `Land_boundary%dfdtr(:,:,:,tr)` | d(tracer flux)/d(tracer_surf) — implicit coupling coefficient for each land tracer flux |

#### Component: ICE — fields sent to sea ice

| Fortran field | Description |
|---|---|
| `Ice_boundary%u_flux` | Zonal wind stress on the ice surface [Pa] |
| `Ice_boundary%v_flux` | Meridional wind stress on the ice surface [Pa] |
| `Ice_boundary%t_flux` | Sensible heat flux into the ice surface [W/m²] |
| `Ice_boundary%q_flux` | Latent heat (moisture) flux into the ice surface [W/m²] |
| `Ice_boundary%lw_flux` | Net longwave flux at the ice surface [W/m²] |
| `Ice_boundary%lw_flux` | Downward longwave flux at the ice surface; alternate override target for the same field [W/m²] |
| `Ice_boundary%sw_flux_nir_dir` | Direct-beam near-infrared shortwave flux over ice [W/m²] |
| `Ice_boundary%sw_flux_vis_dir` | Direct-beam visible shortwave flux over ice [W/m²] |
| `Ice_boundary%sw_flux_nir_dif` | Diffuse near-infrared shortwave flux over ice [W/m²] |
| `Ice_boundary%sw_flux_vis_dif` | Diffuse visible shortwave flux over ice [W/m²] |
| `Ice_boundary%sw_down_vis_dir` | Downward direct-beam visible shortwave flux over ice [W/m²] |
| `Ice_boundary%sw_down_vis_dif` | Downward diffuse visible shortwave flux over ice [W/m²] |
| `Ice_boundary%sw_down_nir_dir` | Downward direct-beam near-infrared shortwave flux over ice [W/m²] |
| `Ice_boundary%sw_down_nir_dif` | Downward diffuse near-infrared shortwave flux over ice [W/m²] |
| `Ice_boundary%lprec` | Liquid precipitation rate over ice [kg/m²/s] |
| `Ice_boundary%fprec` | Frozen precipitation rate over ice [kg/m²/s] |
| `Ice_boundary%dhdt` | d(sensible heat flux)/d(T_surf) — implicit coupling coefficient for ice heat flux [W/m²/K] |
| `Ice_boundary%dedt` | d(latent heat flux)/d(T_surf) — implicit coupling coefficient for ice moisture flux [W/m²/K] |
| `Ice_boundary%drdt` | d(longwave flux)/d(T_surf) — implicit coupling coefficient for ice longwave flux [W/m²/K] |
| `Ice_boundary%coszen` | Cosine of the solar zenith angle over ice [dimensionless] |
| `Ice_boundary%p` | Surface pressure over ice [Pa] |
| `Ice_boundary%fluxes` | Coupler boundary-condition type holding all per-tracer gas and deposition fluxes from atmosphere to ice |

---

### `flux_up_to_atmos`

This subroutine collects updated surface states from land and ice and
transfers them back to the atmosphere after each surface model update.

#### Component: ICE

| Fortran field | Description |
|---|---|
| `Ice%t_surf` | Updated ice surface temperature after the ice model step [K] |

#### Component: LND

| Fortran field | Description |
|---|---|
| `Land%t_ca` | Updated canopy air temperature after the land model step [K] |
| `Land%t_surf` | Updated land surface temperature after the land model step [K] |
| `Land%tr(:,:,:,tr)` | Updated surface tracer mixing ratio over land after the land model step; one entry per exchanged tracer |

---

### `flux_atmos_to_ocean`

This subroutine computes atmosphere-to-ocean/ice deposition gas fluxes and
passes them to the ice boundary type for downstream ocean use.

#### Component: ICE

| Fortran field | Description |
|---|---|
| `Ice_boundary%fluxes%bc(n)%field(m)%values` | Per-gas, per-field entries in the atmosphere-ice flux coupler boundary-condition type; iterated over all boundary conditions and their sub-fields (e.g. gas partial pressures, piston velocities, fluxes) |

---

## `land_ice_flux_exchange.F90`

### `flux_land_to_ice`

This subroutine transfers runoff and calving fluxes from the land model to
the sea-ice/ocean model.

#### Component: ICE

| Fortran field | Description |
|---|---|
| `Land_Ice_Boundary%runoff` | Liquid runoff (river discharge) from land to ocean/ice [kg/m²/s] |
| `Land_Ice_Boundary%calving` | Solid calving flux (iceberg / glacier discharge) from land to ocean/ice [kg/m²/s] |
| `Land_Ice_Boundary%runoff_hflx` | Heat flux carried by liquid runoff [W/m²] |
| `Land_Ice_Boundary%calving_hflx` | Heat flux carried by calving discharge [W/m²] |

---

## `ice_ocean_flux_exchange.F90`

### `flux_ice_to_ocean_finish`

This subroutine finalises the ice-to-ocean boundary by applying any data
overrides to the fluxes that will be passed to the ocean model.

#### Component: OCN

| Fortran field | Description |
|---|---|
| `Ice_Ocean_Boundary%u_flux` | Zonal wind/ice stress on the ocean surface [Pa] |
| `Ice_Ocean_Boundary%v_flux` | Meridional wind/ice stress on the ocean surface [Pa] |
| `Ice_Ocean_Boundary%t_flux` | Sensible heat flux into the ocean [W/m²] |
| `Ice_Ocean_Boundary%q_flux` | Latent heat (freshwater) flux into the ocean [W/m²] |
| `Ice_Ocean_Boundary%salt_flux` | Salt flux from sea-ice to ocean (brine rejection / melting) [kg/m²/s] |
| `Ice_Ocean_Boundary%lw_flux` | Net longwave flux at the ocean surface [W/m²] |
| `Ice_Ocean_Boundary%sw_flux_nir_dir` | Direct-beam near-infrared shortwave flux into the ocean [W/m²] |
| `Ice_Ocean_Boundary%sw_flux_nir_dif` | Diffuse near-infrared shortwave flux into the ocean [W/m²] |
| `Ice_Ocean_Boundary%sw_flux_vis_dir` | Direct-beam visible shortwave flux into the ocean [W/m²] |
| `Ice_Ocean_Boundary%sw_flux_vis_dif` | Diffuse visible shortwave flux into the ocean [W/m²] |
| `Ice_Ocean_Boundary%lprec` | Liquid precipitation flux reaching the ocean surface [kg/m²/s] |
| `Ice_Ocean_Boundary%fprec` | Frozen precipitation flux reaching the ocean surface [kg/m²/s] |
| `Ice_Ocean_Boundary%runoff` | Liquid runoff from land routed to the ocean [kg/m²/s] |
| `Ice_Ocean_Boundary%calving` | Solid calving flux routed to the ocean [kg/m²/s] |
| `Ice_Ocean_Boundary%runoff_hflx` | Heat flux carried by liquid runoff into the ocean [W/m²] |
| `Ice_Ocean_Boundary%calving_hflx` | Heat flux carried by calving into the ocean [W/m²] |
| `Ice_Ocean_Boundary%p` | Surface atmospheric pressure at the ocean surface [Pa] |
| `Ice_Ocean_Boundary%mi` | Sea-ice mass per unit area (used for pressure loading on the ocean) [kg/m²] |
| `Ice_Ocean_Boundary%ustar_berg` | Friction velocity beneath icebergs [m/s]; only present when iceberg module is active |
| `Ice_Ocean_Boundary%area_berg` | Fractional area of ocean cell covered by icebergs [dimensionless]; only present when iceberg module is active |
| `Ice_Ocean_Boundary%mass_berg` | Iceberg mass per unit area [kg/m²]; only present when iceberg module is active |
| `Ice_Ocean_Boundary%fluxes` | Coupler boundary-condition type holding all per-tracer gas fluxes from ice/atmosphere to ocean |

---

### `flux_ocean_to_ice_finish`

This subroutine finalises the ocean-to-ice boundary by applying data
overrides to the ocean state fields passed up to the sea-ice model.

#### Component: ICE

| Fortran field | Description |
|---|---|
| `Ocean_Ice_Boundary%u` | Zonal ocean surface current velocity [m/s] |
| `Ocean_Ice_Boundary%v` | Meridional ocean surface current velocity [m/s] |
| `Ocean_Ice_Boundary%t` | Sea-surface temperature seen by the ice model [K] |
| `Ocean_Ice_Boundary%s` | Sea-surface salinity seen by the ice model [psu] |
| `Ocean_Ice_Boundary%frazil` | Frazil ice heat flux from the ocean to the ice [W/m²] |
| `Ocean_Ice_Boundary%sea_level` | Sea-surface height / sea level [m] |
| `Ocean_Ice_Boundary%fields` | Coupler boundary-condition type holding per-tracer ocean surface fields passed to the ice model |
