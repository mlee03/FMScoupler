# `ice_ocean_driver_type` — Combined Ice–Ocean Driver Control Structure

## Overview

`ice_ocean_driver_type` is the control structure for the combined ice-ocean driver module (`combined_ice_ocean_driver`). It is used only when `combine_ice_and_ocean`, along with `concurrent_ice` and `slow_ice_with_ocean`, are true in `coupler_nml`.

---

## Control Structure Fields

| Field | Type | Default | Description |
|---|---|---|---|
| `Ice_ocean_driver_CS%CS_is_initialized` | logical | `.false.` | `.true.` once `ice_ocean_driver_init` has completed successfully; guards against use before initialization. |
| `Ice_ocean_driver_CS%single_MOM_call` | logical | `.true.` | If `.true.`, MOM6 dynamics and thermodynamics are advanced together in a single call to `update_ocean_model`. If `.false.`, separate calls are made for the two phases. |
| `Ice_ocean_driver_CS%intersperse_ice_ocn` | logical | `.false.` | If `.true.`, ice and ocean thermodynamic and dynamic updates are interspersed rather than sequential. Requires `single_MOM_call=.false.`. |
| `Ice_ocean_driver_CS%use_intersperse_bug` | logical | `.false.` | If `.true.`, retains a bug in the intersperse option where the ocean state was not being passed back to the sea ice between interspersed steps; exists for backward compatibility. |
| `Ice_ocean_driver_CS%dt_coupled_dyn` | real | -1 (s) | Timestep for coupling ice and ocean dynamics when `intersperse_ice_ocn=.true.`; set to < 0 to use the standard coupled timestep `dt_cpld`. |
