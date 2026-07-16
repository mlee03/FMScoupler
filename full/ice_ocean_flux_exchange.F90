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
!! Module ice_ocean_flux_exchange_mod handles data transfer between ice and ocean
!! as well as stock computation
!! @endparblock
module ice_ocean_flux_exchange_mod

  use FMS
  use FMSconstants, only: HLF, HLV, CP_OCEAN
  use ice_model_mod,       only: ice_data_type, ocean_ice_boundary_type
  use ocean_model_mod,     only: ocean_public_type, ice_ocean_boundary_type
  use ocean_model_mod,     only: ocean_state_type, ocean_model_data_ge
  use ocean_model_mod,     only: ocean_model_init_sfc

  implicit none ; private


  public :: ice_ocean_flux_exchange_init, &
       flux_ice_to_ocean, flux_ice_to_ocean_finish, &
       flux_ocean_to_ice, flux_ocean_to_ice_finish, &
       flux_ice_to_ocean_stocks,&
       flux_ocean_from_ice_stocks

  !Balaji, sets boundary_type%xtype
  !  REGRID: grids are physically different, pass via exchange grid
  !  REDIST: same physical grid, different decomposition, must move data around
  !  DIRECT: same physical grid, same domain decomposition, can directly copy data
  integer, parameter :: REGRID=1
    !< is a flag used to indicate ice and ocean are on physically different grids.
    !! Data will be transferred via the exchange grid
  integer, parameter :: REDIST=2
    !< is a flag used to indicate grids for ocean and ice are same but with differen
    !! domain decomposition.  Data will be transferred with fms_mpp_redistribute.
  integer, parameter :: DIRECT=3
    !< is a flag used to indicate grids for ocean and ice are same.
    !! Data can be copied directly.

  logical :: debug_stocks = .false.
    !< is a flag where if .TRUE., call check_flux_conservation at module initialization
  logical :: do_area_weighted_flux = .false.
    !< is a flag where if .TRUE., scale fluxes by source cell area and divide by the destination cell
    !! area to preserve the global area-weighted integral if redistributing.

  integer :: cplOcnClock
    !< is a FMS clock ID to time flux_ice_to_ocean and flux_ocean_to_ice
  integer :: fluxOceanIceClock
    !< is a FMS clock ID to time flux_ocean_to_ice transfer
  integer :: fluxIceOceanClock
    !< is a FMS clock ID to time flux_ice_to_ocean transfer

  real :: Dt_cpl
    !< is the Coupled (slow) timestep in seconds; used in stock computation

  integer, allocatable :: slow_ice_ocean_pelist(:)
    !< is the Combined MPI pelist of the slow-ice and ocean pes;
    !! set during initialization

contains

  !> @parblock
  !! Subroutine ice_ocean_flux_exchange_init initializes the module for flux exchange between ice and ocean.
  !! @endparblock
  subroutine ice_ocean_flux_exchange_init(Time, Ice, Ocean, Ocean_state, ice_ocean_boundary, &
                                          ocean_ice_boundary, Dt_cpl_in, debug_stocks_in,    &
                                          do_area_weighted_flux_in, ex_gas_fields_ice, ex_gas_fluxes, &
                                          do_ocean, slow_ice_ocean_pelist_in )

    type(FmsTime_type), intent(in) :: Time
      !< is the model's current time
    type(ice_data_type), intent(inout) :: Ice
      !< is a derived data type holding ice boundary data
    type(ocean_public_type), intent(inout) :: Ocean
      !< is a derived data type holding ocean boundary data
    type(ocean_state_type), pointer :: Ocean_state
      !< is a pointer pointing to the ocean model's internal state
    type(ice_ocean_boundary_type), intent(inout) :: ice_ocean_boundary
      !< is a derived data type holding properties and fluxes passed from ice to ocean
    type(ocean_ice_boundary_type), intent(inout) :: ocean_ice_boundary
      !< is a derived data type holding properties and fluxes passed from ocean to ice
    real, intent(in) :: Dt_cpl_in
      !< is the coupled (slow) timestep in seconds to set module level Dt_cpl that's used in stock computation.
    logical, intent(in) :: debug_stocks_in
      !< is used to set module level variable debug_stocks.  If TRUE, stocks will be
      !! computed for flux exchange consistency
    logical, intent(in) :: do_area_weighted_flux_in
      !< is used to set module level do_area_weighted_flux.  If TRUE,
      !! flux between ice and ocean will be area-weighted
    type(FmsCoupler1dBC_type), intent(in) :: ex_gas_fields_ice
      !! is used to spawn matching arrays in ocean_ice_boundary and Ocean%fields.
    type(FmsCoupler1dBC_type), intent(in) :: ex_gas_fluxes
      !< is used to spawn matching arrays in ice_ocean_boundary and Ice%ocean_fluxes.
    logical, intent(in) :: do_ocean
      !< is a flag where if .TRUE., ocean_ice_boundary%stagger = ocean%stagger, else defaults to AGRID
    integer, dimension(:), intent(in) :: slow_ice_ocean_pelist_in
      !< is the combined MPI pelist of the slow-ice and ocean processing element used
      !! to set module level slow_ice_ocean_pelis

    integer :: is, ie, js, je

    Dt_cpl = Dt_cpl_in
    debug_stocks = debug_stocks_in
    do_area_weighted_flux = do_area_weighted_flux_in

    !> @parblock
    !! INITIALIZE OCEAN_ICE_BOUNDARY FIELDS TO ZERO.  INITIALIZE T TO 273.0 [K].
    !! @endparblock
    !ocean_ice_boundary and ice_ocean_boundary must be done on all PES
    !domain boundaries will assure no space is allocated on non-relevant PEs.
    call fms_mpp_domains_get_compute_domain( Ice%slow_Domain_NH, is, ie, js, je )
    !allocate ocean_ice_boundary
    allocate( ocean_ice_boundary%u(is:ie,js:je) )
    allocate( ocean_ice_boundary%v(is:ie,js:je) )
    allocate( ocean_ice_boundary%t(is:ie,js:je) )
    allocate( ocean_ice_boundary%s(is:ie,js:je) )
    !frazil and sea_level are optional, if not present they should be nullified
    allocate( ocean_ice_boundary%frazil(is:ie,js:je) )
    allocate( ocean_ice_boundary%sea_level(is:ie,js:je) )
    ! initialize boundary fields for override experiments
    ocean_ice_boundary%u=0.0
    ocean_ice_boundary%v=0.0
    ocean_ice_boundary%t=273.0
    ocean_ice_boundary%s=0.0
    ocean_ice_boundary%frazil=0.0
    ocean_ice_boundary%sea_level=0.0

    !> @parblock
    !! SPAWN GAS FIELDS TO OCEAN_ICE_BOUNDARY%FIELDS.
    !! @endparblock
    if (.not.fms_coupler_type_initialized(ocean_ice_boundary%fields)) &
      call fms_coupler_type_spawn(ex_gas_fields_ice, ocean_ice_boundary%fields, (/is,is,ie,ie/), &
                              (/js,js,je,je/), suffix='_ocn_ice')
    if (Ice%pe) &
      call fms_coupler_type_set_diags(ocean_ice_boundary%fields, "ice_flux", Ice%axes(1:2), Time)

    !> @parblock
    !! SPAWN GAS FLUXES TO ICE%OCEAN_FLUXES.
    !! @endparblock
    if (.not.fms_coupler_type_initialized(Ice%ocean_fluxes)) &
      call fms_coupler_type_spawn(ex_gas_fluxes, Ice%ocean_fluxes, (/is,is,ie,ie/), &
                              (/js,js,je,je/),  suffix = '_ice')

    ! This was never being sent, so comment it out for now.
    ! if (Ice%pe) &
    !   call coupler_type_set_diags(Ice%ocean_fluxes, "ice_flux", Ice%axes(1:2), Time)


    !> @parblock
    !! ALLOCATE ICE_OCEAN_BOUNDARY FIELDS AND INITIALIZE TO ZERO.
    !! IF ICEBERG FIELDS ARE ASSOCIATED IN ICE, ALLOCATE ICE_BERGS FIELDS IN
    !! ICE_OCEAN_BOUNDARY AND INITIALIZE TO ZERO.
    !! @endparblock
    call fms_mpp_domains_get_compute_domain( Ocean%domain, is, ie, js, je )
    !ML ocean only requires t, q, lw, sw, fprec, calving
    !AMIP ocean needs no input fields
    !choice of fields will eventually be done at runtime
    !via field_manager
    allocate( ice_ocean_boundary%u_flux(is:ie,js:je) ); ice_ocean_boundary%u_flux = 0.0
    allocate( ice_ocean_boundary%v_flux(is:ie,js:je) ); ice_ocean_boundary%v_flux = 0.0
    allocate( ice_ocean_boundary%t_flux(is:ie,js:je) ); ice_ocean_boundary%t_flux = 0.0
    allocate( ice_ocean_boundary%q_flux(is:ie,js:je) ); ice_ocean_boundary%q_flux = 0.0
    allocate( ice_ocean_boundary%salt_flux(is:ie,js:je) ); ice_ocean_boundary%salt_flux = 0.0
    allocate( ice_ocean_boundary%lw_flux(is:ie,js:je) ); ice_ocean_boundary%lw_flux = 0.0
    allocate( ice_ocean_boundary%sw_flux_vis_dir(is:ie,js:je) ); ice_ocean_boundary%sw_flux_vis_dir = 0.0
    allocate( ice_ocean_boundary%sw_flux_vis_dif(is:ie,js:je) ); ice_ocean_boundary%sw_flux_vis_dif = 0.0
    allocate( ice_ocean_boundary%sw_flux_nir_dir(is:ie,js:je) ); ice_ocean_boundary%sw_flux_nir_dir = 0.0
    allocate( ice_ocean_boundary%sw_flux_nir_dif(is:ie,js:je) ); ice_ocean_boundary%sw_flux_nir_dif = 0.0
    allocate( ice_ocean_boundary%lprec(is:ie,js:je) ); ice_ocean_boundary%lprec = 0.0
    allocate( ice_ocean_boundary%fprec(is:ie,js:je) ); ice_ocean_boundary%fprec = 0.0
    allocate( ice_ocean_boundary%runoff(is:ie,js:je) ); ice_ocean_boundary%runoff = 0.0
    allocate( ice_ocean_boundary%calving(is:ie,js:je) ); ice_ocean_boundary%calving = 0.0
    allocate( ice_ocean_boundary%runoff_hflx(is:ie,js:je) ); ice_ocean_boundary%runoff_hflx = 0.0
    allocate( ice_ocean_boundary%calving_hflx(is:ie,js:je) ); ice_ocean_boundary%calving_hflx = 0.0
    allocate( ice_ocean_boundary%p(is:ie,js:je) ); ice_ocean_boundary%p = 0.0
    allocate( ice_ocean_boundary%mi(is:ie,js:je) ); ice_ocean_boundary%mi = 0.0
    !Allocating iceberg fields, if the corresponding fields are associated in the sea ice model(s)
    if (associated(Ice%ustar_berg)) then
      allocate( ice_ocean_boundary%ustar_berg (is:ie,js:je) ); ice_ocean_boundary%ustar_berg = 0.0
    endif
    if (associated(Ice%area_berg)) then
      allocate( ice_ocean_boundary%area_berg  (is:ie,js:je) ); ice_ocean_boundary%area_berg = 0.0
    endif
    if (associated(Ice%mass_berg)) then
      allocate( ice_ocean_boundary%mass_berg  (is:ie,js:je) ); ice_ocean_boundary%mass_berg = 0.0
    endif
    ! Copy the stagger indication variables from the ice processors the ocean
    ! PEs and vice versa.  The defaults are large negative numbers, so the
    ! global max here picks out only values that have been set on active PEs.
    call fms_mpp_max(Ice%flux_uv_stagger)
    call fms_mpp_max(Ocean%stagger)
    ice_ocean_boundary%wind_stagger = Ice%flux_uv_stagger
    if(do_ocean) then
       ocean_ice_boundary%stagger = Ocean%stagger
    else
       ocean_ice_boundary%stagger = AGRID
    endif

    !> @parblock
    !! SPAWN GAS FIELDS AND FLUXES TO ICE_OCEAN_BOUNDARY%FLUXES.
    !! @endparblock
    if (.not.fms_coupler_type_initialized(ice_ocean_boundary%fluxes)) &
      call fms_coupler_type_spawn(ex_gas_fluxes, ice_ocean_boundary%fluxes, (/is,is,ie,ie/), &
                              (/js,js,je,je/), suffix='_ice_ocn')
    if (Ocean%is_ocean_pe) &
      call fms_coupler_type_set_diags(ice_ocean_boundary%fluxes, "ocean_flux", Ocean%axes(1:2), Time)

    !> @parblock
    !! SPAWN GAS FIELDS TO OCEAN%FIELDS.
    !! @endparblock
    if (.not.fms_coupler_type_initialized(Ocean%fields)) &
      call fms_coupler_type_spawn(ex_gas_fields_ice, Ocean%fields, (/is,is,ie,ie/), &
                              (/js,js,je,je/), suffix = '_ocn')

    !> @parblock
    !! INITIALIZE BOUNDARY VALUES OCEAN_ICE_BOUNDARY%XTYPE TO DIRECT IF
    !! THE ICE AND OCEAN DOMAINS ARE THE SAME, OTHERWISE REDIST.
    !! (USED IN DATA_OVERRIDE)
    !! @endparblock
    ocean_ice_boundary%xtype = REDIST
    if( Ocean%domain.EQ.Ice%slow_Domain_NH )ocean_ice_boundary%xtype = DIRECT
    ice_ocean_boundary%xtype = ocean_ice_boundary%xtype

    !       initialize the Ocean type for extra fields for surface fluxes
    ! Same allocation of arrays and stuff
    !       (this must be done after the Ocean fields are allocated as the fields on the Ocean%fields
    !       are read in in this subroutine)
    !

    !> @parblock
    !! CALL OCEAN_MODEL_INIT_SFC TO COMPLETE OCEAN SURFACE FIELD INITIALIZATION.
    !! @endparblock
    if ( Ocean%is_ocean_pe ) then
       call fms_mpp_set_current_pelist(Ocean%pelist)
       call ocean_model_init_sfc(Ocean_state, Ocean)
    endif
    call fms_mpp_set_current_pelist()

    !> @parblock
    !! CHECK FLUX CONSERVATION IF DEBUG_STOCKS IS TRUE.
    !! @endparblock
    if(debug_stocks) call check_flux_conservation(Ice, Ocean, Ice_Ocean_Boundary)

    !> @parblock
    !! ALLOCATE SLOW_ICE_OCEAN_PELIST.  INITIALIZE CLOCKS TO MEASURE PERFORMANCE.
    !! @endparblock
    if (Ice%slow_ice_PE .or. Ocean%is_ocean_pe) then
      allocate(slow_ice_ocean_pelist(size(slow_ice_ocean_pelist_in(:))))
      slow_ice_ocean_pelist = slow_ice_ocean_pelist_in
      call fms_mpp_set_current_pelist(slow_ice_ocean_pelist)
      cplOcnClock = fms_mpp_clock_id( 'Ice-ocean coupler', flags=fms_clock_flag_default, grain=CLOCK_COMPONENT )
      fluxIceOceanClock = fms_mpp_clock_id( 'Flux ice to ocean', flags=fms_clock_flag_default, grain=CLOCK_ROUTINE )
      fluxOceanIceClock = fms_mpp_clock_id( 'Flux ocean to ice', flags=fms_clock_flag_default, grain=CLOCK_ROUTINE )
    endif

  end subroutine ice_ocean_flux_exchange_ini


  !> @parblock
  !! Subroutine flux_ice_to_ocean interpolates data from Ice to Ice_Ocean_Boundary in order to
  !! exchange fluxes at the bottom of ice to the ocean model.
  !! The following quantities are transferred from the Ice to Ice_Ocean_Boundary:
  !!       flux_u = zonal wind stress [Pa]
  !!       flux_v = meridional wind stress [Pa]
  !!       flux_t = sensible heat flux [W/m2]
  !!       flux_q = specific humidity flux [Kg/m2/s]
  !!    flux_salt = salt flux [Kg/m2/s]
  !!      flux_sw = net (down-up) shortwave flux [W/m2]
  !!      flux_lw = net (down-up) longwave flux [W/m2]
  !!        lprec = mass of liquid precipitation since last time step [Kg/m2]
  !!        fprec = mass of frozen precipitation since last time step [Kg/m2]
  !!       runoff = mass of runoff since last time step [Kg/m2]
  !!       calving = mass of calving since last time step [Kg/m2]
  !!       p_surf = surface pressure [Pa]
  !! @endparblock
  subroutine flux_ice_to_ocean ( Ice, Ocean, Ice_Ocean_Boundary )

    type(ice_data_type), intent(in) :: Ice
      !< is a derived data type containg ice boundary data
    type(ocean_public_type), intent(in) :: Ocean
      !< is a derived data type to containing ocean boundary data
    type(ice_ocean_boundary_type), intent(inout) :: Ice_Ocean_Boundary
      !< is a derived data type to specify properties and fluxes passed from ice to ocean

    integer :: m
    integer :: n
    logical :: used

    call fms_mpp_clock_begin(cplOcnClock)
    call fms_mpp_clock_begin(fluxIceOceanClock)

    if(ASSOCIATED(Ice_Ocean_Boundary%u_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_u, Ice_Ocean_Boundary%u_flux, Ice_Ocean_Boundary%xtype, .FALSE. )

    if(ASSOCIATED(Ice_Ocean_Boundary%v_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_v, Ice_Ocean_Boundary%v_flux, Ice_Ocean_Boundary%xtype, .FALSE. )

    if(ASSOCIATED(Ice_Ocean_Boundary%p) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%p_surf, Ice_Ocean_Boundary%p     , Ice_Ocean_Boundary%xtype, .FALSE. )

    if(ASSOCIATED(Ice_Ocean_Boundary%mi) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%mi,     Ice_Ocean_Boundary%mi    , Ice_Ocean_Boundary%xtype, .FALSE. )
    ! Extra fluxes
    if (Ice_Ocean_Boundary%xtype == DIRECT) then
       call fms_coupler_type_copy_data(Ice%ocean_fluxes, Ice_Ocean_Boundary%fluxes)
    else
       call fms_coupler_type_redistribute_data(Ice%ocean_fluxes, Ice%slow_Domain_NH, &
                     Ice_Ocean_Boundary%fluxes, ocean%Domain, complete=.true.)
    endif

    !! @endparblock
    !--- The following variables may require conserved flux exchange from ice to ocean because the
    !--- ice area maybe different from ocean area.
    if(ASSOCIATED(Ice_Ocean_Boundary%t_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_t, Ice_Ocean_Boundary%t_flux, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%salt_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_salt, Ice_Ocean_Boundary%salt_flux, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%sw_flux_nir_dir) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_sw_nir_dir, Ice_Ocean_Boundary%sw_flux_nir_dir, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%sw_flux_nir_dif) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_sw_nir_dif, Ice_Ocean_Boundary%sw_flux_nir_dif, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%sw_flux_vis_dir) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_sw_vis_dir, Ice_Ocean_Boundary%sw_flux_vis_dir, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%sw_flux_vis_dif) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_sw_vis_dif, Ice_Ocean_Boundary%sw_flux_vis_dif, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%lw_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_lw, Ice_Ocean_Boundary%lw_flux, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%lprec) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%lprec, Ice_Ocean_Boundary%lprec, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%fprec) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%fprec, Ice_Ocean_Boundary%fprec, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%runoff) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%runoff, Ice_Ocean_Boundary%runoff, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%calving) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%calving, Ice_Ocean_Boundary%calving, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%ustar_berg) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
       Ice%ustar_berg, Ice_Ocean_Boundary%ustar_berg, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%area_berg) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%area_berg, Ice_Ocean_Boundary%area_berg, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%mass_berg) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%mass_berg, Ice_Ocean_Boundary%mass_berg, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%runoff_hflx) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%runoff_hflx, Ice_Ocean_Boundary%runoff_hflx, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%calving_hflx) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%calving_hflx, Ice_Ocean_Boundary%calving_hflx, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%q_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_q, Ice_Ocean_Boundary%q_flux, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    call fms_mpp_clock_end(fluxIceOceanClock)
    call fms_mpp_clock_end(cplOcnClock)

  end subroutine flux_ice_to_ocean

  !> @parblock
  !! Subroutine flux_ice_to_ocean_finish mainly calls fms_data_override to override fluxes in Ice_Ocean_Boundary
  !! before transferring flux from Ice to Ocean. NOTE, fms_data_override will only override data if field entry
  !! is found in the data_table.  This subroutine is only called by the ocean pe.
  !! @endparblock
  subroutine flux_ice_to_ocean_finish ( Time, Ice_Ocean_Boundary )

    type(FmsTime_type), intent(in) :: Time
      !< is the current time
    type(ice_ocean_boundary_type), intent(inout) :: Ice_Ocean_Boundary
      !< is a derived data type containing fluxes and properties passed from ice to ocean

    call fms_data_override('OCN', 'u_flux', Ice_Ocean_Boundary%u_flux, Time )
    call fms_data_override('OCN', 'v_flux', Ice_Ocean_Boundary%v_flux, Time )
    call fms_data_override('OCN', 't_flux', Ice_Ocean_Boundary%t_flux, Time )
    call fms_data_override('OCN', 'q_flux', Ice_Ocean_Boundary%q_flux, Time )
    call fms_data_override('OCN', 'salt_flux', Ice_Ocean_Boundary%salt_flux, Time )
    call fms_data_override('OCN', 'lw_flux', Ice_Ocean_Boundary%lw_flux, Time )
    call fms_data_override('OCN', 'sw_flux_nir_dir', Ice_Ocean_Boundary%sw_flux_nir_dir, Time )
    call fms_data_override('OCN', 'sw_flux_nir_dif', Ice_Ocean_Boundary%sw_flux_nir_dif, Time )
    call fms_data_override('OCN', 'sw_flux_vis_dir', Ice_Ocean_Boundary%sw_flux_vis_dir, Time )
    call fms_data_override('OCN', 'sw_flux_vis_dif', Ice_Ocean_Boundary%sw_flux_vis_dif, Time )
    call fms_data_override('OCN', 'lprec', Ice_Ocean_Boundary%lprec, Time )
    call fms_data_override('OCN', 'fprec', Ice_Ocean_Boundary%fprec, Time )
    call fms_data_override('OCN', 'runoff', Ice_Ocean_Boundary%runoff, Time )
    call fms_data_override('OCN', 'calving', Ice_Ocean_Boundary%calving, Time )
    call fms_data_override('OCN', 'runoff_hflx', Ice_Ocean_Boundary%runoff_hflx, Time )
    call fms_data_override('OCN', 'calving_hflx', Ice_Ocean_Boundary%calving_hflx, Time )
    call fms_data_override('OCN', 'p', Ice_Ocean_Boundary%p, Time )
    call fms_data_override('OCN', 'mi', Ice_Ocean_Boundary%mi, Time )

    if (ASSOCIATED(Ice_Ocean_Boundary%ustar_berg) ) &
      call fms_data_override('OCN', 'ustar_berg', Ice_Ocean_Boundary%ustar_berg, Time )
    if (ASSOCIATED(Ice_Ocean_Boundary%area_berg)  ) &
      call fms_data_override('OCN', 'area_berg',  Ice_Ocean_Boundary%area_berg , Time )
    if (ASSOCIATED(Ice_Ocean_Boundary%mass_berg)  ) &
      call fms_data_override('OCN', 'mass_berg',  Ice_Ocean_Boundary%mass_berg , Time )

    call fms_coupler_type_data_override('OCN', Ice_Ocean_Boundary%fluxes, Time )

    call fms_coupler_type_send_data(Ice_Ocean_Boundary%fluxes, Time )

  end subroutine flux_ice_to_ocean_finish

  !#######################################################################
  !> @parblock
  !! Subroutine flux_ocean_to_ice interpolates data from Ocean to Ocean_Ice_Boundary in order to exchange fluxes
  !! from ocean to bottom of ice.  The following quantities are remapped from the Ocean to Ocean_Ice_Boundary:
  !!        t_surf = surface temperature [deg K]
  !!        frazil = frazil fluxes since the last coupling step [J/m2]
  !!        u_surf = zonal ocean current/ice motion [m/s]
  !!        v_surf = meridional ocean current/ice motion [m/s]
  !!       sea_lev = sea level used to drive ice accelerations [m]
  !! @endparblock
  subroutine flux_ocean_to_ice ( Ocean, Ice, Ocean_Ice_Boundary )

    type(ocean_public_type), intent(in) :: Ocean
      !< is a derived data type holding ocean boundary data
    type(ice_data_type), intent(in) :: Ice
      !< is a derived data type holding ice boundary data
    type(ocean_ice_boundary_type), intent(inout) :: Ocean_Ice_Boundary
      !< is a derived data type holding properties and fluxes passed from ocean to ice

    real, allocatable, dimension(:,:) :: tmp
    integer :: m
    integer :: n
    logical :: used

    call fms_mpp_clock_begin(cplOcnClock)
    call fms_mpp_clock_begin(fluxOceanIceClock)

    select case (Ocean_Ice_Boundary%xtype)
    case(DIRECT)
       !same grid and domain decomp for ocean and ice
       if( ASSOCIATED(Ocean_Ice_Boundary%u) )Ocean_Ice_Boundary%u = Ocean%u_surf
       if( ASSOCIATED(Ocean_Ice_Boundary%v) )Ocean_Ice_Boundary%v = Ocean%v_surf
       if( ASSOCIATED(Ocean_Ice_Boundary%t) )Ocean_Ice_Boundary%t = Ocean%t_surf
       if( ASSOCIATED(Ocean_Ice_Boundary%s) )Ocean_Ice_Boundary%s = Ocean%s_surf
       if( ASSOCIATED(Ocean_Ice_Boundary%sea_level) )Ocean_Ice_Boundary%sea_level = Ocean%sea_lev
       if( ASSOCIATED(Ocean_Ice_Boundary%frazil) ) then
          if(do_area_weighted_flux) then
             Ocean_Ice_Boundary%frazil = Ocean%frazil * Ocean%area
             call divide_by_area(data=Ocean_Ice_Boundary%frazil, area=Ice%area)
          else
             Ocean_Ice_Boundary%frazil = Ocean%frazil
          endif
       endif

       ! Extra fluxes
       call fms_coupler_type_copy_data(Ocean%fields, Ocean_Ice_Boundary%fields)

    case(REDIST)
       !same grid, different domain decomp for ocean and ice
       if( ASSOCIATED(Ocean_Ice_Boundary%u) )                     &
            call fms_mpp_domains_redistribute(Ocean%Domain, Ocean%u_surf, Ice%slow_Domain_NH, Ocean_Ice_Boundary%u)
       if( ASSOCIATED(Ocean_Ice_Boundary%v) )                     &
            call fms_mpp_domains_redistribute(Ocean%Domain, Ocean%v_surf, Ice%slow_Domain_NH, Ocean_Ice_Boundary%v)
       if( ASSOCIATED(Ocean_Ice_Boundary%t) )                     &
            call fms_mpp_domains_redistribute(Ocean%Domain, Ocean%t_surf, Ice%slow_Domain_NH, Ocean_Ice_Boundary%t)
       if( ASSOCIATED(Ocean_Ice_Boundary%s) )                     &
            call fms_mpp_domains_redistribute(Ocean%Domain, Ocean%s_surf, Ice%slow_Domain_NH, Ocean_Ice_Boundary%s)

       if( ASSOCIATED(Ocean_Ice_Boundary%sea_level) )             &
            call fms_mpp_domains_redistribute(Ocean%Domain, Ocean%sea_lev, Ice%slow_Domain_NH, &
                                              Ocean_Ice_Boundary%sea_level)
       if( ASSOCIATED(Ocean_Ice_Boundary%frazil) ) then
          if(do_area_weighted_flux) then
             if (Ocean%is_ocean_pe) then
               allocate(tmp(size(Ocean%area,1), size(Ocean%area,2)))
               tmp(:,:) = Ocean%frazil(:,:) * Ocean%area(:,:)
             endif
             call fms_mpp_domains_redistribute( Ocean%Domain, tmp, Ice%slow_Domain_NH, Ocean_Ice_Boundary%frazil)
             if (Ice%slow_ice_pe) &
               call divide_by_area(data=Ocean_Ice_Boundary%frazil, area=Ice%area)
             if (Ocean%is_ocean_pe) deallocate(tmp)
          else
             call fms_mpp_domains_redistribute(Ocean%Domain,Ocean%frazil, Ice%slow_Domain_NH, Ocean_Ice_Boundary%frazil)
          endif
       endif

       ! Extra fluxes
       call fms_coupler_type_redistribute_data(Ocean%fields, Ocean%Domain, &
                     Ocean_Ice_Boundary%fields, Ice%slow_Domain_NH)
    case DEFAULT
       call fms_mpp_error( FATAL, 'flux_ocean_to_ice: Ocean_Ice_Boundary%xtype must be DIRECT or REDIST.' )
    end selec

    call fms_mpp_clock_end(fluxOceanIceClock)
    call fms_mpp_clock_end(cplOcnClock)

  end subroutine flux_ocean_to_ice

  !> @parblock
  !! Subroutine flux_ocean_to_ice_finish carrries out a final set of tasks that should only occur on
  !! the slow-ice processors, including data override and perhaps saving diagnostics.
  !! @endparblock
  subroutine flux_ocean_to_ice_finish( Time, Ice, Ocean_Ice_Boundary )

    type(FmsTime_type), intent(in) :: Time
      !< is the current time
    type(ice_data_type), intent(in) :: Ice
      !< is a derived type holding ice boundary data
    type(ocean_ice_boundary_type), intent(inout) :: Ocean_Ice_Boundary
      !< is a derived type holding properties and fluxes passed from ocean to ice
    real :: from_dq

    call fms_data_override('ICE', 'u',         Ocean_Ice_Boundary%u,         Time)
    call fms_data_override('ICE', 'v',         Ocean_Ice_Boundary%v,         Time)
    call fms_data_override('ICE', 't',         Ocean_Ice_Boundary%t,         Time)
    call fms_data_override('ICE', 's',         Ocean_Ice_Boundary%s,         Time)
    call fms_data_override('ICE', 'frazil',    Ocean_Ice_Boundary%frazil,    Time)
    call fms_data_override('ICE', 'sea_level', Ocean_Ice_Boundary%sea_level, Time)
    call fms_coupler_type_data_override('ICE', Ocean_Ice_Boundary%fields, Time)

    !  Perform diagnostic output for the ocean_ice_boundary fields
    call fms_coupler_type_send_data( Ocean_Ice_Boundary%fields, Time)

    ! frazil (already in J/m^2 so no need to multiply by Dt_cpl)
    from_dq = SUM( Ice%area * Ocean_Ice_Boundary%frazil )
    fms_stock_constants_ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) = &
            fms_stock_constants_ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) - from_dq
    fms_stock_constants_ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) = &
            fms_stock_constants_ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) + from_dq

  end subroutine flux_ocean_to_ice_finish


  !#######################################################################
  !> @parblock
  !> Subroutine flux_ice_to_ocean_stocks integrates the fluxes from ice to ocean over the surface and in time.
  !! Ice stocks are decremented at the base of the ice and incremented to the ocean stocks at the ocean surface.
  !! @endparblock
  subroutine flux_ice_to_ocean_stocks(Ice)

    type(ice_data_type), intent(in) :: Ice
      !< A derived data type to specify ice boundary data

    real :: from_dq

    !> @parblock
    !! COMPUTE STOCKS CHANGE FOR QUANTITY (PRECIP - EVAP).
    !! @endparblock
    from_dq = Dt_cpl * SUM( Ice%area * (Ice%lprec+Ice%fprec-Ice%flux_q) )
    fms_stock_constants_ice_stock(ISTOCK_WATER)%dq(ISTOCK_BOTTOM) = &
            fms_stock_constants_ice_stock(ISTOCK_WATER)%dq(ISTOCK_BOTTOM) - from_dq
    fms_stock_constants_ocn_stock(ISTOCK_WATER)%dq(ISTOCK_TOP   ) = &
            fms_stock_constants_ocn_stock(ISTOCK_WATER)%dq(ISTOCK_TOP   ) + from_dq

    !> @parblock
    !! COMPUTE STOCKS FOR RIVER.
    !! @endparblock
    from_dq = Dt_cpl * SUM( Ice%area * (Ice%runoff + Ice%calving) )
    fms_stock_constants_ice_stock(ISTOCK_WATER)%dq(ISTOCK_BOTTOM) = &
            fms_stock_constants_ice_stock(ISTOCK_WATER)%dq(ISTOCK_BOTTOM) - from_dq
    fms_stock_constants_ocn_stock(ISTOCK_WATER)%dq(ISTOCK_SIDE  ) = &
            fms_stock_constants_ocn_stock(ISTOCK_WATER)%dq(ISTOCK_SIDE  ) + from_dq

    !> @parblock
    !! COMPUTE STOCKS FOR HEAT (SENSIBLE + SHORTWAVE + LONGWAVE + LATENT).
    !! @endparblock
    from_dq = Dt_cpl * SUM( Ice%area * ( &
         &   Ice%flux_sw_vis_dir+Ice%flux_sw_vis_dif &
         & + Ice%flux_sw_nir_dir+Ice%flux_sw_nir_dif + Ice%flux_lw &
         & - (Ice%fprec + Ice%calving)*HLF - Ice%flux_t - Ice%flux_q*HLV) )
    fms_stock_constants_ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) = &
            fms_stock_constants_ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) - from_dq
    fms_stock_constants_ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) = &
            fms_stock_constants_ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) + from_dq

    !> @parblock
    !! COMPUTE STOCKS FOR HEAT FROM RADIATIVE AND TURBLENT FLUXES AND HEAT CARRIED BY
    !! RIVER AND PME (assuming reference temperature of 0 degC and river/pme temp = surface temp).
    !! Note: it does not matter what the ref temperature is but it must be consistent with that in OCN and ICE.
    !! PME = preciptation minus evaporation
    !! @endparblock
    from_dq = Dt_cpl * SUM( Ice%area * ( &
         & (Ice%lprec+Ice%fprec-Ice%flux_q + Ice%runoff+Ice%calving)*CP_OCEAN*Ice%SST_C(:,:)) )
    fms_stock_constants_ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) = &
            fms_stock_constants_ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) - from_dq
    fms_stock_constants_ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) = &
            fms_stock_constants_ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) + from_dq

    !> @parblock
    !! COMPUTE STOCKS FOR FLUX_SALT.
    !! @endparblock
    from_dq = Dt_cpl* SUM( Ice%area * ( -Ice%flux_salt ))
    fms_stock_constants_ice_stock(ISTOCK_SALT)%dq(ISTOCK_BOTTOM) = &
            fms_stock_constants_ice_stock(ISTOCK_SALT)%dq(ISTOCK_BOTTOM) - from_dq
    fms_stock_constants_ocn_stock(ISTOCK_SALT)%dq(ISTOCK_TOP   ) = &
            fms_stock_constants_ocn_stock(ISTOCK_SALT)%dq(ISTOCK_TOP   ) + from_dq


  end subroutine flux_ice_to_ocean_stocks


  !> @parblock
  !! Subroutine flux_ocean_from_ice_stocks updates stocks in Ocean after flux transfer from Ice.
  !! Unlike subroutine flux_ice_to_ocean_stocks() that uses Ice%fluxes to update the stocks, this
  !! subroutine uses Ice_Ocean_boundary%fluxes to calculate the amount of input to Ocean. These fluxes
  !! are the ones that Ocean model uses internally to calculate its budgets. Hence there should be no difference between
  !! this input and what Ocean model internal diagnostics uses.
  !! This bypasses the possible mismatch in cell areas between Ice and Ocean in diagnosing the stocks of Ocean
  !! and should report a conserving Ocean component regardless of the glitches in fluxes.
  !! The use of this subroutine in conjunction with  subroutine flux_ice_to_ocean_stocks() will also allow to directly
  !! diagnose the amount "stocks lost in exchange" between Ice and Ocean.
  !! @endparblock
  subroutine flux_ocean_from_ice_stocks(ocean_state,Ocean,Ice_Ocean_boundary)
    type(ocean_state_type), pointer :: ocean_state
      !< is a pointer to the ocean model's internal state; used to retrieve
      !! ocean-side grid and flux data via ocean_model_data_get.
    type(ocean_public_type), intent(in) :: Ocean
      !< is a derived type containing the Ocean public boundary data type; provides the MPI domain and
      !! ocean pe information to get compute domain.
    type(ice_ocean_boundary_type), intent(in) :: Ice_Ocean_Boundary
      !< is a derived type containing fluxes passed from ice to ocean.

    real :: from_dq, cp_ocn
    real, dimension(size(Ice_Ocean_Boundary%lprec,1), size(Ice_Ocean_Boundary%lprec,2)) :: &
      ocean_cell_area, wet, t_surf, t_pme, t_calving, t_runoff, btfHea
    integer :: isc, iec, jsc, jec

    !> @parblock
    !! USE THE RETRIEVER FROM OCEAN_MODEL_MOD TO GET AREA, MASK, SURFACE TEMPERATURE,
    !! PME TEMPERATURE, CALVING TEMPERATURE, RUNOFF TEMPERATURE, BOTTOM HEAT FLUX,
    !! AND SPECIFIC HEAT CAPACITY FIELDS FROM THE OCEAN MODEL.
    !! @endparblock
    call fms_mpp_domains_get_compute_domain(Ocean%Domain, isc, iec, jsc, jec)
    call ocean_model_data_get(ocean_state,Ocean,'area'  , ocean_cell_area,isc,jsc)
    call ocean_model_data_get(ocean_state,Ocean,'mask', wet,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'t_surf', t_surf,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'t_runoff', t_runoff,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'t_pme', t_pme,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'t_calving', t_calving,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'btfHeat', btfHeat,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'c_p', cp_ocn )

    ! fluxes from ice -> ocean, integrate over surface and in time

    !> @parblock
    !! COMPUTE STOCK TRANSFER OF (PRECIP - EVAP) TO OCEAN SURFACE AND FROM LATERALLY.
    !! @endparblock
    !precip - evap
    from_dq = SUM(ocean_cell_area * wet * (Ice_Ocean_Boundary%lprec+Ice_Ocean_Boundary%fprec-Ice_Ocean_Boundary%q_flux))
    fms_stock_constants_ocn_stock(ISTOCK_WATER)%dq_IN(ISTOCK_TOP   ) = &
            fms_stock_constants_ocn_stock(ISTOCK_WATER)%dq_IN(ISTOCK_TOP   ) + from_dq * Dt_cpl

    from_dq = SUM( ocean_cell_area * wet * (Ice_Ocean_Boundary%runoff+Ice_Ocean_Boundary%calving) )
    fms_stock_constants_ocn_stock(ISTOCK_WATER)%dq_IN(ISTOCK_SIDE  ) = &
            fms_stock_constants_ocn_stock(ISTOCK_WATER)%dq_IN(ISTOCK_SIDE  ) + from_dq * Dt_cpl

    !> @parblock
    !! COMPUTE STOCK TRANSFER OF (SENSIBLE HEAT + SHORTWAVE + LONGWAVE + LATENT HEAT) TO OCEAN LATERALLY.
    !! @endparblock
    from_dq = SUM( ocean_cell_area * wet *( Ice_Ocean_Boundary%sw_flux_vis_dir + Ice_Ocean_Boundary%sw_flux_vis_dif &
         +Ice_Ocean_Boundary%sw_flux_nir_dir + Ice_Ocean_Boundary%sw_flux_nir_dif &
         +Ice_Ocean_Boundary%lw_flux &
         - (Ice_Ocean_Boundary%fprec + Ice_Ocean_Boundary%calving)*HLF &
         - Ice_Ocean_Boundary%t_flux - Ice_Ocean_Boundary%q_flux*HLV ))
    fms_stock_constants_ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE  ) = &
            fms_stock_constants_ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE  ) + from_dq * Dt_cpl

    !> @parblock
    !! COMPUTE STOCK TRANSFER OF HEAT CARRIED BY RIVER + PME (ASSUMING REFERENCE TEMPERATURE OF 0 DEGC
    !! AND RIVER/PME TEMP = SURFACE TEMP).
    !! @endparblock
    ! Note: it does not matter what the ref temperature is but it must be consistent with that in OCN and ICE
    from_dq = SUM( ocean_cell_area * wet * cp_ocn *&
         ((Ice_Ocean_Boundary%lprec+Ice_Ocean_Boundary%fprec-Ice_Ocean_Boundary%q_flux)*t_pme &
         +Ice_Ocean_Boundary%calving * t_calving &
         +Ice_Ocean_Boundary%runoff  * t_runoff  ))
    fms_stock_constants_ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE  ) = &
            fms_stock_constants_ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE ) + from_dq * Dt_cpl

    !> @parblock
    !! COMPUTE STOCK TRANSFER OF BOTTOM HEAT FLUX.
    !! @endparblock
    from_dq = - SUM( ocean_cell_area * wet * btfHeat)
    fms_stock_constants_ocn_stock(ISTOCK_HEAT)%dq_IN( ISTOCK_BOTTOM ) = &
            fms_stock_constants_ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_BOTTOM ) + from_dq * Dt_cpl

    !> @parblock
    !! COMPUTE STOCK TRANSFER OF FRAZIL HEAT.
    !! @endparblock
    from_dq =  SUM( ocean_cell_area *wet * Ocean%frazil )
    fms_stock_constants_ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE  ) = &
            fms_stock_constants_ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE ) + from_dq

    !> @parblock
    !! COMPUTE STOCK TRANSFER OF SALT FLUX.
    !! @endparblock
    from_dq = SUM( ocean_cell_area * wet * ( -Ice_Ocean_Boundary%salt_flux))
    fms_stock_constants_ocn_stock(ISTOCK_SALT)%dq_IN(ISTOCK_TOP  ) = &
            fms_stock_constants_ocn_stock(ISTOCK_SALT)%dq_IN(ISTOCK_TOP   ) + from_dq  * Dt_cpl


  end subroutine flux_ocean_from_ice_stocks

  !#######################################################################
  !> @parblock
  !! Subroutine flux_ice_to_ocean_redistribute performs a globally conservative flux redistribution across ICE/OCN.
  !! If domain decomposition is identical for ocean and ice, data is copied from Ice to ICE_OCEAN_BOUNDARY
  !! If domain decomposition differs, data is copied from Ice to ICE_OCEAN_BOUNDARY with
  !! fms_mpp_domains_redistribute.  (Assumes Ice and Ocean are on the same grid.)
  !! @endparblock
  subroutine flux_ice_to_ocean_redistribute(Ice, Ocean, ice_data, ocn_bnd_data, type, do_area_weighted )

    ! Performs a globally conservative flux redistribution across ICE/OCN.
    ! Assumes that the ice/ocn grids are the same. If ocean is present,
    ! then assume different mpp domans and redistribute

    ! should be invoked by all PEs

    type(ice_data_type), intent(in) :: Ice
      !< is the ice boundary data type; provides the slow-ice MPI domain
      !! (slow_Domain_NH) and cell areas (Ice%area) used in redistribution.
    type(ocean_public_type), intent(in) :: Ocean
      !< is the ocean public boundary data type; provides the ocean MPI domain
      !! (Ocean%Domain) and cell areas (Ocean%area) used in redistribution.
    real, dimension(:,:), intent(in) :: ice_data
      !< is the flux field on the ice domain to be transferred to the ocean boundary.
    real, dimension(:,:), intent(out) :: ocn_bnd_data
      !< is the flux field on the ocean domain; filled with the redistributed
      !! (and optionally area-weighted) values from ice_data.
    integer, intent(in) :: type
      !< is the transfer type: DIRECT (same MPI decomposition, copy directly) or
      !! REDIST (same grid, different MPI decomposition, use mpp_redistribute).
    logical, intent(in) :: do_area_weighted
      !< is a flag where if .TRUE., scale flux by ice cell area before redistribution and divide by
      !! ocean cell area after, preserving the global area-weighted integral.

    real, allocatable, dimension(:,:) :: tmp
      ! Temporary work array used on ice PEs to hold ice_data * ice%area
      ! before MPI redistribution when do_area_weighted is .TRUE.

    select case(type)
    case (DIRECT)
       if(do_area_weighted) then
          ocn_bnd_data = ice_data * ice%area
          call divide_by_area(data=ocn_bnd_data, area=ocean%area)
       else
          ocn_bnd_data = ice_data
       endif
    case (REDIST)
       if (do_area_weighted) then
         if ( Ice%slow_ice_pe ) then
            allocate(tmp(size(ice%area,1), size(ice%area,2)))
            tmp(:,:) = ice_data(:,:) * ice%area(:,:)
          endif
          call fms_mpp_domains_redistribute(Ice%slow_Domain_NH, tmp, ocean%Domain, ocn_bnd_data)
          if (ocean%is_ocean_pe) call divide_by_area(ocn_bnd_data, area=ocean%area)
          if (Ice%slow_ice_pe) deallocate(tmp)
       else
          call fms_mpp_domains_redistribute(Ice%slow_Domain_NH, ice_data, ocean%Domain, ocn_bnd_data)
       endif
    case DEFAULT
       call fms_mpp_error( FATAL, 'FLUX_ICE_TO_OCEAN: Ice_Ocean_Boundary%xtype must be DIRECT or REDIST.' )
    end selec

  end subroutine flux_ice_to_ocean_redistribute

  !> @parblock
  !! Subroutine divide_by_area divides data by area for area > 0.
  !! @endparblock
  subroutine divide_by_area(data, area)
    real, intent(inout) :: data(:,:)
      !< is the data to divide by area; modified in-place
    real, intent(in) :: area(:,:)
      !< is the area field to divide by

    !> @parblock
    !! IF DATA AND AREA DIFFER IN SIZE, RETURN WITHOUT MODIFYING DATA
    !! @endparblock
    if(size(data, dim=1) /= size(area, dim=1) .or. size(data, dim=2) /= size(area, dim=2)) then
       return
    endif

    !> @parblock
    !! WHERE(AREA /= 0.0) DATA = DATA / AREA
    !! @endparblock
    where(area /= 0.0)
       data = data / area
    end where

  end subroutine divide_by_area


  !> @parblock
  !! Subroutine check_flux_conservation checks for flux conservation
  !! after flux_ice_to_ocean_redistribute.
  !! @endparblock
  subroutine check_flux_conservation(Ice, Ocean, Ice_Ocean_Boundary)
    type(ice_data_type), intent(inout) :: Ice
      !< is the Ice boundary data type; provides the ice MPI domain, cell areas
      !! (Ice%area), and flux array sizes used to allocate test data.
    type(ocean_public_type), intent(inout) :: Ocean
      !< is the Ocean public boundary data type; provides the ocean MPI domain
      !! and cell areas (Ocean%area) used to compute redistributed sums.
    type(ice_ocean_boundary_type), intent(inout) :: Ice_Ocean_Boundary
      !< is the Ice-to-ocean boundary type; provides xtype (DIRECT or REDIST)
      !! and q_flux array size used to allocate the ocean-side test buffer.

    real, allocatable, dimension(:,:) :: ice_data
      ! Random ice_data that will be filled with random numbers for testing
    real, allocatable, dimension(:,:) :: ocn_data
      ! Random test data to receive from ice_data
    real :: ice_sum
      ! Global area-weighted sum of ice_data on the ice domain
      ! [ice_data units × m2]; serves as the conservation reference value.
    real :: area_weighted_sum
      ! Global area-weighted sum of ocn_data after redistribution with do_area_weighted=.true.
    real :: non_area_weighted_sum
      ! Global area-weighted sum of ocn_data after redistribution with do_area_weighted=.false.
    integer :: outuni
      ! Fortran unit number for stdout; used to write the diagnostic report.

    !> @parblock
    !! SET OUTUNIT TO STDOUT.
    !! @endparblock
    outunit = fms_mpp_stdout()

    !> @parblock
    !! ALLOCATE ICE_DATA AND OCN_DATA FOR TESTING.
    !! @endparblock
    allocate(ice_data(size(Ice%flux_q,1), size(Ice%flux_q,2) ) )
    allocate(ocn_data(size(Ice_Ocean_Boundary%q_flux,1), size(Ice_Ocean_Boundary%q_flux,2) ) )

    !> @parblock
    !! INITIALIZE ICE_DATA WITH RANDOM NUMBERS.
    !! @endparblock
    call random_number(ice_data)
    ice_sum = sum(ice_data*Ice%area)
    call fms_mpp_sum(ice_sum)

    !> @parblock
    !! CALL FLUX_ICE_TO_OCEAN_DISTRIBUTE WITH AREA_WEIGHTED_SUM = .FALSE. AND GET GLOBAL SUM.
    !! @endparblock
    ocn_data = 0.0
    call flux_ice_to_ocean_redistribute( Ice, Ocean, ice_data, ocn_data, Ice_Ocean_Boundary%xtype, .false.)
    non_area_weighted_sum = sum(ocn_data*Ocean%area)
    call fms_mpp_sum(non_area_weighted_sum)

    !> @parblock
    !! CALL FLUX_ICE_TO_OCEAN_DISTRIBUTE WITH AREA_WEIGHTED_SUM = .TRUE. AND GET GLOBAL SUM.
    !! @endparblock
    ocn_data = 0.0
    call flux_ice_to_ocean_redistribute( Ice, Ocean, ice_data, ocn_data, Ice_Ocean_Boundary%xtype, .true.)
    area_weighted_sum = sum(ocn_data*Ocean%area)
    call fms_mpp_sum(area_weighted_sum)

    !> @parblock
    !! WRITE REPORT TO OUTUNIT.
    !! @endparblock
    write(outunit,*)"NOTE from flux_exchange_mod: check for flux conservation for flux_ice_to_ocean"
    write(outunit,*)"***** The global area sum of random number on ice domain (input data) is ", ice_sum
    write(outunit,*)"***** The global area sum of data after flux_ice_to_ocean_redistribute with "// &
         "do_area_weighted_flux = false is ", non_area_weighted_sum, &
         " and the difference from global input area sum = ", ice_sum - non_area_weighted_sum
    write(outunit,*)"***** The global area sum of data after flux_ice_to_ocean_redistribute with "// &
         "do_area_weighted_flux = true is ", area_weighted_sum, &
         " and the difference from global input area sum = ", ice_sum - area_weighted_sum

  end subroutine check_flux_conservation

end module ice_ocean_flux_exchange_mod
