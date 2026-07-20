!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS) Coupler.
!*
!* FMS Coupler is free software: you can redistribute it and/or modify
!* it under the terms of the GNU Lesser General Public License as
!* published by the Free Software Foundation, either version 3 of the
!* License, or (at your option) any later version.
!*
!* FMS Coupler is distributed in the hope that it will be useful, bu
!* WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!* General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS Coupler.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> @file
!> @parblock
!! Module land_ice_flux_exchange_mod handles freshwater discharge (runoff and calving) and
!! associated heat exchanges via the exchange grid between land and ice grids.
 !! @endparblock
module land_ice_flux_exchange_mod

!! FMS
  use FMS
  use FMSconstants, only: RADIUS
  use land_model_mod,      only: land_data_type
  use ice_model_mod,       only: ice_data_type, land_ice_boundary_type

  implicit none
  private

  type(FmsXgridXmap_type), save :: xmap_runoff
    !< is the exchange grid map between land and ice/ocean.
  integer :: n_xgrid_runoff=0
    !< is the number of exchange grid cells in xmap_runoff.

  integer :: X2_GRID_LND
    !< is the land index for xmap_runoff; used to identify
    !! the source side when calling fms_xgrid_stock_move.  set to 1
  integer :: X2_GRID_ICE
    !< is the ice/ocean index for xmap_runoff; used to identify
    !! the destination side when calling fms_xgrid_stock_move.  set to 2

  public :: flux_land_to_ice, land_ice_flux_exchange_init

  integer :: cplClock
    !< is the clock ID for timing flux_land_to_ice calls.
  integer :: fluxLandIceClock
    !< is the clock ID for timing the flux_land_to_ice transfer
  logical :: do_runoff
    !< is a flag where if .TRUE., land discharge is transferred to ice.
    !! If .FALSE., all runoff/calving fields are zeroed.
  real :: Dt_cpl
    !< is the coupled (slow) timestep in seconds; used in stock computation.

contains

  !> @parblock
  !! Subrutine land_ice_flux_exchange_init initializes the land-ice flux exchange module
  !! for flux exchange between land and ice.
   !! @endparblock
  subroutine land_ice_flux_exchange_init(Land, Ice, land_ice_boundary, Dt_cpl_in, do_runoff_in, cplClock_in)
    type(land_data_type), intent(in)    :: Land
      !< is a derived type holding land boundary data
    type(ice_data_type), intent(inout) :: Ice
      !< is a derived type holding ice boundary data
    type(land_ice_boundary_type), intent(inout) :: land_ice_boundary
      !< is a derived type holding properties and fluxes passed from land to ice
    real, intent(in)  :: Dt_cpl_in
      !< is the coupled (slow) timestep in seconds to set module level dt_cpl
    logical, intent(in) :: do_runoff_in
      !< is a flag to set module level do_runoff.
    integer, intent(in) :: cplClock_in
      !< is the FMS MPP clock id for the top-level coupler profiling clock; stored module-wide so tha
      !! flux_land_to_ice can bracket its work.

    integer :: is, ie, js, je

    !> @parblock
    !! SET DO_RUNOFF, CPLCLOCK, DT_CPL, AND FLUXLANDICECLOCK.
    !! @endparblock
    do_runoff = do_runoff_in
    cplClock = cplClock_in
    Dt_cpl   = Dt_cpl_in
    fluxLandIceClock = fms_mpp_clock_id( 'Flux land to ice', flags=fms_clock_flag_default, grain=CLOCK_ROUTINE )

    if (do_runoff) then
       call fms_xgrid_setup_xmap(xmap_runoff, (/ 'LND', 'OCN' /),       &
            (/ Land%Domain, Ice%Domain /),                    &
            "INPUT/grid_spec.nc"             )
       ! exchange grid indices
       X2_GRID_LND = 1; X2_GRID_ICE = 2;
       n_xgrid_runoff = max(fms_xgrid_count(xmap_runoff),1)
       if (n_xgrid_runoff.eq.1) write (*,'(a,i6,6x,a)') 'PE = ', fms_mpp_pe(), 'Runoff  exchange size equals one.'
    endif

    call fms_mpp_domains_get_compute_domain( Ice%domain, is, ie, js, je )

    !> @parblock
    !! ALLOCATE LAND_ICE_BOUNDARY%RUNOFF, CALVING, RUNOFF_HFLX, AND CALVING_HFLX, AND INITIALIZE THEM TO ZERO.
    !! @endparblock
    allocate( land_ice_boundary%runoff(is:ie,js:je) )
    allocate( land_ice_boundary%calving(is:ie,js:je) )
    allocate( land_ice_boundary%runoff_hflx(is:ie,js:je) )
    allocate( land_ice_boundary%calving_hflx(is:ie,js:je) )
    land_ice_boundary%runoff=0.0
    land_ice_boundary%calving=0.0
    land_ice_boundary%runoff_hflx=0.0
    land_ice_boundary%calving_hflx=0.0

  end subroutine land_ice_flux_exchange_init

  !#######################################################################
  !> @parblock
  !! Subroutine flux_land_to_ice handles conservative transfer of water and snow discharge
  !! from land to sea ice/ocean. The following elements are transferred from the Land to the Land_ice_boundary:
  !!        discharge to runoff (kg/m2).
  !!        discharge_snow to calving (kg/m2).
   !! @endparblock
  subroutine flux_land_to_ice( Time, Land, Ice, Land_Ice_Boundary )
    type(FmsTime_type),  intent(in) :: Time
      !< is the current time
    type(land_data_type), intent(in) :: Land
      !< is a derived type holding land boundary data
    type(ice_data_type), intent(in) :: Ice
      !< is a derived type holding ice boundary data
    type(land_ice_boundary_type), intent(inout):: Land_Ice_Boundary
      !< is a derived type holding properties and fluxes passed from land to ice

    integer :: ier
      ! Error code returned by fms_xgrid_stock_move; non-zero indicates a stock accounting error.
    real, dimension(n_xgrid_runoff) :: ex_runoff
      ! Liquid runoff (kg/m2) on the exchange grid, gathered from the land domain.
    real, dimension(n_xgrid_runoff) :: ex_calving
      ! Snow discharge / calving (kg/m2) on the exchange grid, gathered from the land domain.
    real, dimension(n_xgrid_runoff) :: ex_runoff_hflx
      ! Heat flux associated with liquid runoff (W/m2) on the exchange grid.
    real, dimension(n_xgrid_runoff) :: ex_calving_hflx
      ! Heat flux associated with snow discharge (W/m2) on the exchange grid.
    real, dimension(size(Land_Ice_Boundary%runoff,1),size(Land_Ice_Boundary%runoff,2),1) :: ice_buf
      ! Temporary 3-D buffer used to receive exchange-grid data and copy it into the 2-D ice-domain fields.

    !> @parblock
    !! INITIALIZE CLOCK.
    !! @endparblock
    call fms_mpp_clock_begin(cplClock)
    call fms_mpp_clock_begin(fluxLandIceClock)

    !> @parblock
    !! IF DO_RUNOFF, TRANSFER DATA, ELSE SET LAND_ICE_BOUNDARY%RUNOFF, CALVING
    !! RUNOFF_HFLX, AND CALVING_HFLX TO ZERO.
    !! @endparblock
    if (do_runoff) then
        !> @parblock
        !! TRANSFER DISCHARGE* FIELDS FROM THE LAND TO ICE VIA THE EXCHANGE GRID.
        !! @endparblock
       call fms_xgrid_put_to_xgrid ( Land%discharge,      'LND', ex_runoff,  xmap_runoff)
       call fms_xgrid_put_to_xgrid ( Land%discharge_snow, 'LND', ex_calving, xmap_runoff)
       call fms_xgrid_put_to_xgrid ( Land%discharge_heat,      'LND', ex_runoff_hflx,  xmap_runoff)
       call fms_xgrid_put_to_xgrid ( Land%discharge_snow_heat, 'LND', ex_calving_hflx, xmap_runoff)
       call fms_xgrid_get_from_xgrid (ice_buf, 'OCN', ex_runoff,  xmap_runoff)
       Land_Ice_Boundary%runoff = ice_buf(:,:,1);
       call fms_xgrid_get_from_xgrid (ice_buf, 'OCN', ex_calving, xmap_runoff)
       Land_Ice_Boundary%calving = ice_buf(:,:,1);
       call fms_xgrid_get_from_xgrid (ice_buf, 'OCN', ex_runoff_hflx,  xmap_runoff)
       Land_Ice_Boundary%runoff_hflx = ice_buf(:,:,1);
       call fms_xgrid_get_from_xgrid (ice_buf, 'OCN', ex_calving_hflx, xmap_runoff)
       Land_Ice_Boundary%calving_hflx = ice_buf(:,:,1);

       !> @parblock
       !! OVERRIDE TRANSFERRED DATA WITH DATA_OVERRIDE IF FIELD EXISTS IN DATA_TABLE.
       !! @endparblock
       call fms_data_override('ICE', 'runoff' , Land_Ice_Boundary%runoff , Time)
       call fms_data_override('ICE', 'calving', Land_Ice_Boundary%calving, Time)
       call fms_data_override('ICE', 'runoff_hflx' , Land_Ice_Boundary%runoff_hflx , Time)
       call fms_data_override('ICE', 'calving_hflx', Land_Ice_Boundary%calving_hflx, Time)

       !> @parblock
       !! COMPUTE WATER STOCK ON THE EXCHANGE GRID TO MEASURE STOCK BEING
       !! TRANSFERRED FROM LAND TO ICE.
       !! @endparblock
       ice_buf(:,:,1) = Land_Ice_Boundary%runoff + Land_Ice_Boundary%calving
       call fms_xgrid_stock_move(from=fms_stock_constants_lnd_stock(ISTOCK_WATER), &
            & to=fms_stock_constants_ice_stock(ISTOCK_WATER), &
            & grid_index=X2_GRID_ICE, &
            & stock_data3d=ice_buf, &
            & xmap=xmap_runoff, &
            & delta_t=Dt_cpl, &
            & from_side=ISTOCK_SIDE, to_side=ISTOCK_SIDE, &
            & radius=Radius, ier=ier, verbose='stock move RUNOFF+CALVING (Lnd->Ice) ')
    else
       Land_Ice_Boundary%runoff = 0.0
       Land_Ice_Boundary%calving = 0.0
       Land_Ice_Boundary%runoff_hflx = 0.0
       Land_Ice_Boundary%calving_hflx = 0.0
    endif

    !> @parblock
    !! END CLOCK.
    !! @endparblock
    call fms_mpp_clock_end(fluxLandIceClock)
    call fms_mpp_clock_end(cplClock)

  end subroutine flux_land_to_ice


!#######################################################################

end module land_ice_flux_exchange_mod
