# `ocean_public_type` — Ocean Model Public Surface Fields (MOM6)

## Overview

`ocean_public_type` is the publicly visible face of the MOM6 ocean model. It contains only the **surface fields and domain metadata** needed by the coupler; the full MOM6 interior state is held in the private `ocean_state_type`. An instance named `Ocean` is declared in `coupler_main.F90`, with `Ocean_state` as a separate pointer to the interior state.

**Defined in:** `ocean_model_MOM.F90`  
**Populated by:** `update_ocean_model`  
**Read by:** `flux_ocean_to_ice` to construct `ocean_ice_boundary_type`  
**Related types:** `ocean_state_type`, `ocean_ice_boundary_type`

---

## Domain and PE Metadata

| Field | Type | Description |
|---|---|---|
| `Ocean%Domain` | `type(domain2d)` | FMS domain decomposition for the ocean surface fields; defines the MPI tile layout used for coupler remapping. |
| `Ocean%is_ocean_pe` | logical | `.true.` on processors that run the ocean model; used throughout the coupler to gate ocean-only code paths. |
| `Ocean%pelist` | `integer(:)` | List of MPI PE numbers assigned to the ocean model. |
| `Ocean%maskmap` | `logical(:,:)` (pointer) | Mask indicating which logical processors are active for ocean computation; logical processors covering all-land points may not be mapped to physical PEs. Need not be set if all processors are used. |
| `Ocean%instance_name` | `character(32)` | Optional name identifying this ocean model instance; used in ensemble runs to disambiguate log messages. |
| `Ocean%stagger` | integer | Arakawa staggering of the surface velocity components (`u_surf`, `v_surf`) relative to tracer points. Valid values: `AGRID`, `BGRID_NE`, `CGRID_NE`, `BGRID_SW`, `CGRID_SW`. Set to -999 before initialization so a global max can propagate the value to non-ocean PEs. |

---

## Ocean Surface State Fields

All fields are set by the ocean model after each `update_ocean_model` call and read by the coupler to construct the `ocean_ice_boundary_type` passed to the sea-ice model via `flux_ocean_to_ice`.

| Field | Units | Description |
|---|---|---|
| `Ocean%t_surf` | K | Sea-surface temperature (SST) on tracer (T) cells. |
| `Ocean%s_surf` | ppt | Sea-surface salinity (SSS) on T cells. |
| `Ocean%u_surf` | m/s | Surface current i-velocity at the locations indicated by `stagger`. |
| `Ocean%v_surf` | m/s | Surface current j-velocity at the locations indicated by `stagger`. |
| `Ocean%sea_lev` | m | Sea level corrected for surface pressure: `dzt(1) + η + p_atm/(ρ₀g)`; passed to the ice model as `sea_level`. |
| `Ocean%frazil` | J/m² | Accumulated heating from frazil ice formation in the ocean since the last coupling step; delivered to the ice model so it can account for ocean-side freezing. |
| `Ocean%melt_potential` | J/m² | Instantaneous heat available to melt sea ice from below; computed when the ocean boundary layer depth exceeds `HFrz`. |
| `Ocean%OBLD` | m | Ocean boundary layer depth; used to determine the depth over which melt potential is computed. |
| `Ocean%area` | m² | Grid-cell area of each ocean surface cell; used for conservative flux remapping. |

---

## Calving Fields

| Field | Units | Description |
|---|---|---|
| `Ocean%calving` | kg/m² | Mass per unit area of ice-shelf flux to be converted to icebergs; passed to the iceberg module. |
| `Ocean%calving_hflx` | W/m² | Heat flux associated with calving. |

---

## Tracer and Diagnostic Fields

| Field | Type | Description |
|---|---|---|
| `Ocean%fields` | `type(coupler_2d_bc_type)` | Named arrays of tracer-related ocean surface fields (e.g., pCO₂, O₂ saturation) used in atmosphere-ocean gas flux calculations; populated by `atmos_ocean_fluxes_calc`. |
| `Ocean%avg_kount` | integer | Counter tracking the number of contributions accumulated in the running averages stored in this type; used externally by FMSCoupler to manage time averaging. |
| `Ocean%axes(2)` | `integer(2)` | Diag-manager axis IDs available for I/O using this surface data. |
