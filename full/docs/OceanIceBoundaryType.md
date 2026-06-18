# ocean_ice_boundary_type
ocean_ice_boundary_type carries ocean surface state passed from MOM6 to the sea-ice model (SIS2)


## Ocean_ice_boundary_type%u
Ocean_ice_boundary_type%u, a real 2D array, is x-direction ocean surface velocity at a position determined by stagger [m/s].
## Ocean_ice_boundary_type%v
Ocean_ice_boundary_type%v, a real 2D array, is y-direction ocean surface velocity at a position determined by stagger [m/s].
## Ocean_ice_boundary_type%stagger
Ocean_ice_boundary_type%stagger, integer, is Spatial staggering of u and v relative to tracer points; default is BGRID_NE.

## Ocean_ice_boundary_type%t
Ocean_ice_boundary_type%t, a real 2D array, is Ocean surface temperature [K].
## Ocean_ice_boundary_type%s
Ocean_ice_boundary_type%s, a real 2D array, is Ocean surface salinity [g salt / kg seawater].
## Ocean_ice_boundary_type%frazil
Ocean_ice_boundary_type%frazil, a real 2D array, is Frazil heat rejected by the ocean since the last coupling step [J/m**2]; delivered to SIS2 so it can account for ocean-side freezing.
## Ocean_ice_boundary_type%sea_level
Ocean_ice_boundary_type%sea_level, a real 2D array, is Sea level after adjustment for any surface pressure that the ocean allows to be expressed [m].

## Ocean_ice_boundary_type%calving
Ocean_ice_boundary_type%calving, a real 2D array, is Mass flux per unit area of ice-shelf flux to be converted to icebergs [kg/m**2/s].
## Ocean_ice_boundary_type%calving_hflx
Ocean_ice_boundary_type%calving_hflx, a real 2D array, is Heat flux associated with calving W/m**2.

## Ocean_ice_boundary_type%data
Ocean_ice_boundary_type%data, a real 3D array, is Collective array providing named access to the scalar fields above; used internally for data-override and exchange-grid operations.
## Ocean_ice_boundary_type%xtype
Ocean_ice_boundary_type%xtype, integer, is Transfer mode for the ocean-to-ice exchange: REGRID (1), REDIST (2), or DIRECT (3).
## Ocean_ice_boundary_type%fields
Ocean_ice_boundary_type%fields, type(coupler_2d_bc_type), is Named array of additional per-tracer ocean surface fields (e.g., pCO₂, SSS for gas exchange) passed to the ice model.
