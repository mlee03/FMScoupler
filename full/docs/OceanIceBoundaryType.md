# `ocean_ice_boundary_type` — Ocean-to-Ice Boundary Fields

## Overview

`ocean_ice_boundary_type` carries ocean surface state passed from the ocean model to the sea-ice model. 

---

## Surface Velocity Fields

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Ocean_ice_boundary%u` | real 2D | m/s | x-direction ocean surface velocity at a position determined by `stagger`. |
| `Ocean_ice_boundary%v` | real 2D | m/s | y-direction ocean surface velocity at a position determined by `stagger`. |
| `Ocean_ice_boundary%stagger` | integer | — | Spatial staggering of `u` and `v` relative to tracer points; default is `BGRID_NE`. |

---

## Ocean Surface Thermodynamic Fields

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Ocean_ice_boundary%t` | real 2D | K | Ocean surface temperature. |
| `Ocean_ice_boundary%s` | real 2D | g salt/kg seawater | Ocean surface salinity. |
| `Ocean_ice_boundary%frazil` | real 2D | J/m² | Frazil heat rejected by the ocean since the last coupling step; delivered to SIS2 so it can account for ocean-side freezing. |
| `Ocean_ice_boundary%sea_level` | real 2D | m | Sea level after adjustment for any surface pressure that the ocean allows to be expressed. |

---

## Calving Fields

| Field | Type / Dimensions | Units | Description |
|---|---|---|---|
| `Ocean_ice_boundary%calving` | real 2D | kg/m²/s | Mass flux per unit area of ice-shelf flux to be converted to icebergs. |
| `Ocean_ice_boundary%calving_hflx` | real 2D | W/m² | Heat flux associated with calving. |

---

## Transfer and Metadata Fields

| Field | Type | Description |
|---|---|---|
| `Ocean_ice_boundary%data` | real 3D | Collective array providing named access to the scalar fields above; used internally for data-override and exchange-grid operations. |
| `Ocean_ice_boundary%xtype` | integer | Transfer mode for the ocean-to-ice exchange: `REGRID` (1), `REDIST` (2), or `DIRECT` (3). |
| `Ocean_ice_boundary%fields` | `type(coupler_2d_bc_type)` | Named array of additional per-tracer ocean surface fields (e.g., pCO₂, SSS for gas exchange) passed to the ice model. |
