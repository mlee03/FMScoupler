# `ice_data_type` — Sea-Ice Model State (SIS2)

## Overview

`ice_data_type` is the publicly visible face of the sea-ice model.

---

## PE and Domain Metadata

| Field | Type | Description |
|---|---|---|
| `Ice%pe` | logical | `.true.` on any PE that participates in ice computation (fast or slow). |
| `Ice%fast_ice_pe` | logical | `.true.` on PEs in the fast-ice pelist; these PEs handle atmosphere-timestep ice processes. |
| `Ice%slow_ice_pe` | logical | `.true.` on PEs in the slow-ice pelist; these PEs handle coupled-timestep ice dynamics and ice-ocean exchange. |
| `Ice%shared_slow_fast_PEs` | logical | `.true.` when fast and slow ice use the same PE set and domain decomposition (`slow_ice_with_ocean=.false.`); `.false.` when slow ice runs on the ocean PEs. |
| `Ice%pelist` | `integer(:)` | Combined list of all ice PEs (union of fast and slow pelists); used for flux exchange. |
| `Ice%fast_pelist` | `integer(:)` | MPI PE numbers for fast-ice processes. |
| `Ice%slow_pelist` | `integer(:)` | MPI PE numbers for slow-ice processes. |
| `Ice%Domain` | `type(domain2D)` | Copy of the fast-ice FMS domain without halos; used for exchange-grid setup. |
| `Ice%slow_Domain_NH` | `type(domain2D)` | Copy of the slow-ice FMS domain without halos; used for ice-ocean flux redistribution. |
| `Ice%fast_domain` | `type(domain2D)` (pointer) | Pointer to the fast-ice MPI domain (or an allocated copy on slow-ice PEs). |
| `Ice%slow_domain` | `type(domain2D)` (pointer) | Pointer to the slow-ice MPI domain (or an allocated copy on fast-ice PEs). |
| `Ice%ocean_pt` | `logical(:,:)` | Mask array; `.true.` at ocean (non-land) points. |
| `Ice%xtype` | integer | Transfer mode for ice-ocean flux exchange: `DIRECT` (3) when ice and ocean share the same decomposition, `REDIST` (2) when they differ. |
| `Ice%axes` | `integer(3)` | Diag-manager axis IDs for the ice surface grid. |
| `Ice%Time` | `type(time_type)` | The sea-ice model's current clock time. |

---

## Per-Category Surface Fields (Atmosphere–Ice Interface)

These fields are 3D arrays where the third index holds the paritition (e.g., 1 for open water)

| Field | Units | Description |
|---|---|---|
| `Ice%part_size` | dimensionless, 0–1 | Fractional coverage of the grid cell by each ice thickness category. |
| `Ice%t_surf` | K | Surface temperature of the ocean (category 1) or each ice thickness category. |
| `Ice%albedo` | dimensionless, 0–1 | Broadband surface albedo averaged across all wavelength and orientation bands within each ice category. |
| `Ice%albedo_vis_dir` | dimensionless | Surface albedo for direct visible shortwave radiation in each ice category. |
| `Ice%albedo_nir_dir` | dimensionless | Surface albedo for diffuse visible shortwave radiation in each ice category. |
| `Ice%albedo_vis_dif` | dimensionless | Surface albedo for direct near-infrared shortwave radiation in each ice category. |
| `Ice%albedo_nir_dif` | dimensionless | Surface albedo for diffuse near-infrared shortwave radiation in each ice category. |
| `Ice%rough_mom` | m | Surface roughness length for momentum at the ocean/ice surface, as provided by `ocean_rough_mod`. |
| `Ice%rough_heat` | m | Surface roughness length for heat at the ocean/ice surface. |
| `Ice%rough_moist` | m | Surface roughness length for moisture at the ocean/ice surface. |
| `Ice%u_surf` | m/s | Eastward surface velocity of the ocean (category 1) or sea ice; used as a lower boundary condition for wind stress. |
| `Ice%v_surf` | m/s | Northward surface velocity of the ocean (category 1) or sea ice. |
| `Ice%flux_uv_stagger` | integer | Staggering of `u_surf`/`v_surf` relative to tracer points; valid values are `AGRID`, `BGRID_NE`, `CGRID_NE`, `BGRID_SW`, `CGRID_SW` (Arakawa notation). Initialized to -999 so that a global max across all PEs propagates the correct value to non-ice PEs. |

---

## Scalar Ocean-Side Fields

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Ice%s_surf` | real 2D | g salt/kg seawater | Ocean surface salinity. |
| `Ice%SST_C` | real 2D | °C | Ocean surface temperature; used in forecast mode. |
| `Ice%area` | real 2D | m² | Area of each ocean cell; land cells have area = 0 and this field can double as a mask. |
| `Ice%mi` | real 2D | kg/m² | Total ice + snow mass per unit area; passed to the ocean for pressure loading and to the wave model. |

---

## Ice-to-Ocean Flux Fields (Computed by the Slow-Ice Model)

| Field | Units | Description |
|---|---|---|
| `Ice%flux_u` | Pa | Flux of x-momentum into the ocean (zonal wind/ice stress). |
| `Ice%flux_v` | Pa | Flux of y-momentum into the ocean (meridional wind/ice stress). |
| `Ice%flux_t` | W/m² | Sensible heat flux out of the ocean. |
| `Ice%flux_q` | kg/m²/s | Evaporative moisture flux out of the ocean. |
| `Ice%flux_lh` | W/m² | Latent heat flux out of the ocean. |
| `Ice%flux_lw` | W/m² | Net longwave flux out of the ocean. |
| `Ice%flux_sw_vis_dir` | W/m² | Direct visible shortwave heat flux into the ocean. |
| `Ice%flux_sw_vis_dif` | W/m² | Diffuse visible shortwave heat flux into the ocean. |
| `Ice%flux_sw_nir_dir` | W/m² | Direct near-infrared shortwave heat flux into the ocean. |
| `Ice%flux_sw_nir_dif` | W/m² | Diffuse near-infrared shortwave heat flux into the ocean. |
| `Ice%flux_salt` | kg/m²/s | Salt flux out of the ocean (brine rejection / melting). |
| `Ice%salt_left_behind` | kg/m²/s | Salt remaining in the ocean during ice growth (not ejected as brine). |
| `Ice%lprec` | kg/m² | Liquid precipitation flux into the ocean. |
| `Ice%fprec` | kg/m² | Frozen precipitation flux into the ocean. |
| `Ice%runoff` | kg/m² | Liquid runoff from land into the ocean. |
| `Ice%runoff_hflx` | W/m² | Heat flux associated with liquid runoff, relative to a reference temperature. |
| `Ice%runoff_carbon` | kg/m² | Carbon content of liquid runoff into the ocean. |
| `Ice%calving` | kg/m² | Calving of ice or frozen freshwater runoff into the ocean. |
| `Ice%calving_hflx` | W/m² | Heat flux associated with calving, relative to a reference temperature. |
| `Ice%p_surf` | Pa | Pressure at the ocean surface; may or may not include atmospheric pressure depending on configuration. |
| `Ice%stress_mag` | Pa | Time-mean magnitude of the ice-ocean stress; passed to the ocean when `pass_stress_mag = .true.` in `SIS_slow_CS`. |

---

## Iceberg Fields

Allocated only when the iceberg module is active (`do_icebergs = .true.`).

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Ice%ustar_berg` | real 2D | m/s | Friction velocity contribution beneath icebergs; used for iceberg-ocean drag. |
| `Ice%area_berg` | real 2D | m²/m² | Fraction of the grid cell covered by icebergs. |
| `Ice%mass_berg` | real 2D | kg/m² | Mass of icebergs per unit area. |

---

## Gas Exchange and Tracer Fields

| Field | Type | Description |
|---|---|---|
| `Ice%ocean_fields` | `type(coupler_3d_bc_type)` | Named surface-state fields shared with the atmosphere (SST, SSS, piston velocities, etc.) for atmosphere-ocean gas exchange; populated by `flux_ocean_to_ice`. |
| `Ice%ocean_fluxes` | `type(coupler_2d_bc_type)` | Computed gas fluxes from the ice to the ocean for additional tracers. |
| `Ice%ocean_fluxes_top` | `type(coupler_3d_bc_type)` | Gas flux boundary conditions at the top of the ice (atmosphere side); flagged for eventual removal. |

---

## Private Control Structures

These fields are private to the ice model and opaque to the coupler. They are listed here for reference only.

| Field | Type | Description |
|---|---|---|
| `Ice%fCS` | `type(SIS_fast_CS)` (pointer) | Control structure for the SIS2 fast ice thermodynamics; lives on atmosphere PEs. Contains the fast-ice grid (`fCS%G`), diagnostics, and all fast-timestep state. |
| `Ice%sCS` | `type(SIS_slow_CS)` (pointer) | Control structure for the SIS2 slow ice dynamics and thermodynamics; may live on ocean PEs when `slow_ice_with_ocean=.true.`. Contains the slow-ice grid (`sCS%G`), dynamics, tracer registry, and all slow-timestep state. |
| `Ice%icebergs` | `type(icebergs)` (pointer) | Control structure for the Lagrangian iceberg module; null when `do_icebergs=.false.`. |
| `Ice%US` | `type(unit_scale_type)` (pointer) | SIS2 dimensional unit-scaling factors; converts between external MKS values and SIS2 internal non-dimensionalised units. |
| `Ice%Ice_restart` | `type(SIS_restart_CS)` (pointer) | Control structure for writing and reading slow-ice restart files. |
| `Ice%Ice_fast_restart` | `type(SIS_restart_CS)` (pointer) | Control structure for writing and reading fast-ice restart files. |
| `Ice%OBC` | `type(ice_OBC_type)` (pointer) | Control structure for ice open boundary conditions; null when OBCs are not configured. |
| `Ice%restart_output_dir` | `character(240)` | Directory path for restart file output; default is `./RESTART/`. |
