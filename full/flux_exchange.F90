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
!* FMS Coupler is distributed in the hope that it will be useful, but
!* WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!* General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS Coupler.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> @defgroup flux_exchange_mod flux_exchange_mod
!! @{
!! @parblock
!! Flux_exchange_mod is the top level module for flux exchange between components
!! @endparblock
module flux_exchange_mod

  use FMS
  use FMSconstants, only: rdgas, rvgas, cp_air, stefan, WTMAIR, &
                          HLV, HLF, Radius, PI, CP_OCEAN, WTMCO2, WTMC

  !! Components
  use land_model_mod,             only: Lnd_stock_pe
  use ocean_model_mod,            only: Ocean_stock_pe
  use atmos_model_mod,            only: Atm_stock_pe
  use atm_land_ice_flux_exchange_mod, only: atm_land_ice_flux_exchange_init, sfc_boundary_layer
  use atm_land_ice_flux_exchange_mod, only: generate_sfc_xgrid, flux_down_from_atmos
  use atm_land_ice_flux_exchange_mod, only: flux_up_to_atmos, atm_stock_integrate, send_ice_mask_sic
  use atm_land_ice_flux_exchange_mod, only: flux_atmos_to_ocean, flux_ex_arrays_dealloc
  use land_ice_flux_exchange_mod,     only: flux_land_to_ice, land_ice_flux_exchange_init
  use ice_ocean_flux_exchange_mod,    only: ice_ocean_flux_exchange_init
  use ice_ocean_flux_exchange_mod,    only: flux_ocean_to_ice, flux_ocean_to_ice_finish
  use ice_ocean_flux_exchange_mod,    only: flux_ice_to_ocean, flux_ice_to_ocean_finish
  use ice_ocean_flux_exchange_mod,    only: flux_ice_to_ocean_stocks, flux_ocean_from_ice_stocks
  use atmos_model_mod,    only: atmos_data_type, land_ice_atmos_boundary_type
  use ocean_model_mod,    only: ocean_public_type, ice_ocean_boundary_type
  use ocean_model_mod,    only: ocean_state_type
  use ice_model_mod,      only: ice_data_type, land_ice_boundary_type, &
                                ocean_ice_boundary_type, atmos_ice_boundary_type, Ice_stock_pe
  use land_model_mod,     only: land_data_type, atmos_land_boundary_type
  use atmos_ocean_fluxes_mod,     only: atmos_ocean_fluxes_init, atmos_ocean_type_fluxes_init
  use atmos_ocean_fluxes_calc_mod, only: atmos_ocean_fluxes_calc
  use ocean_model_mod,            only: ocean_model_init_sfc, ocean_model_flux_init
  use atmos_tracer_driver_mod,    only: atmos_tracer_flux_init

  implicit none ; private

  public :: flux_exchange_init, gas_exchange_init, &
     sfc_boundary_layer,   &
     generate_sfc_xgrid,   &
     flux_down_from_atmos, &
     flux_up_to_atmos,     &
     flux_land_to_ice,     &
     flux_atmos_to_ocean,  &
     flux_ex_arrays_dealloc,&
     flux_ice_to_ocean,    &
     flux_ice_to_ocean_finish, &
     flux_ocean_to_ice,    &
     flux_ocean_to_ice_finish, &
     flux_check_stocks,    &
     flux_init_stocks,     &
     flux_ice_to_ocean_stocks,&
     flux_ocean_from_ice_stocks,&
     send_ice_mask_sic

  !-----------------------------------------------------------------------
  character(len=128) :: version = '$Id$'
   !< is the program version string set automatically at compile time.
  character(len=128) :: tag = '$Name$'
   !< is a string set automatically at compile time.

  logical :: do_init = .true.
   !< is a flag where if .TRUE., initialize module

  real, parameter :: bound_tol = 1e-7
   !< is the tolerance value used when checking grid-boundary coordinate consistency.

  real, parameter :: d622 = rdgas/rvgas
   !< is the ratio of dry-air and water-vapor gas constants.

  real, parameter :: d378 = 1.0-d622
   !< is the complement of d622, used in humidity conversions.

  real :: z_ref_heat =  2.
   !< is the reference height [m] for temperature and relative humidity diagnostics
   !! (t_ref, rh_ref, del_h, del_q).

  real :: z_ref_mom  = 10.
   !< is the reference height [m] for momentum diagnostics (u_ref, v_ref, del_m).

  logical :: do_area_weighted_flux = .FALSE.
   !< is a namelist flag where if .TRUE., normalize exchanged fluxes by the area;
   !! used in ice_ocean_flux_exchange.

  logical :: debug_stocks = .FALSE.
   !< is a namelist flag where if .TRUE., enable extra stock-conservation output for debugging.

  logical :: divert_stocks_report = .FALSE.
   !< is a namelist flag where if .TRUE., write stock reports 'stocks.out'; else write to stdout.

  logical :: do_runoff = .TRUE.
   !< is a namelist flag where if .TRUE., turn on the land runoff interpolation to the ocean

  logical :: do_forecast = .false.
   !< is a namelist flag.

  integer :: nblocks = 1
   !< is a namelist variable for number of OpenMP blocks, defaults to 1.

  logical :: partition_fprec_from_lprec = .FALSE.
   !< is a namelist flag where if .TRUE., convert liquid precip to snow when t_ref is less than
   !! tfreeze parameter

  real, parameter :: tfreeze = 273.15
   !< is the freezing point of water at one atmosphere in [K].

  logical :: scale_precip_2d = .false.
   !< is a namelist flag where if .TRUE., rescale liquid precipitation using a 2-D field from data override.

  namelist /flux_exchange_nml/ z_ref_heat, z_ref_mom,&
       & do_area_weighted_flux, debug_stocks, divert_stocks_report, do_runoff, do_forecast, nblocks,&
       & partition_fprec_from_lprec, scale_precip_2d

   logical :: gas_fluxes_initialized = .false.
   !< is a flag to indicate component fluxes have been initialized.

   type(FmsCoupler1dBC_type), target :: ex_gas_fields_atm
   !< is a derived type containing atmospheric surface variables that are used in
   !! calculating atmosphere-ocean gas fluxes.

   type(FmsCoupler1dBC_type), target :: ex_gas_fields_ice
   !< is a derived type containing ice-top and ocean surface variables that are used
   !! in calculating atmosphere-ocean gas fluxes.

  type(FmsCoupler1dBC_type), target :: ex_gas_fluxes
   !< is a derived type for exchanging gas or tracer fluxes between the atmosphere and ocean,
   !! defined by the field table.  Also a place holder of intermediate calculations.

   integer :: ni_atm
      !< is the number of x gridpoints in the atm compute domain
   integer :: nj_atm
      !< is the number of y gridpoints in the atm compute domain
   real, dimension(3) :: ccc
      !< is a temporary array used for conservation-check summaries; not used.

  integer :: cplClock
   !< is the FMS clock id to profile land-ice-atmos coupler.

  real :: Dt_atm
  !< is the atmosphere timestep [s]

  real :: Dt_cpl
   !< is the coupled timesteps in [s].

  real :: ATM_PRECIP_NEW
   !< is used to take into account implicit evaporation in stock computation

contains

  !#######################################################################
  !> @parblock
  !! Subroutine gas_exchange_init initializes the fms_atmos_ocean_type_fluxes,
  !! ocean_model_fluxes, and atmos_tracer_flux.  The subroutine also calls
  !! fms_atmos_ocean_fluxes to initialize ex_gas_fluxes and fields
   !! @endparblock
  subroutine gas_exchange_init (gas_fields_atm, gas_fields_ice, gas_fluxes)
    type(FmsCoupler1dBC_type), optional, pointer :: gas_fields_atm
      !< is a derived type containing atmospheric surface variables that
      !! are used in computing atmosphere-ocean gas fluxes.
    type(FmsCoupler1dBC_type), optional, pointer :: gas_fields_ice
      !< is a derived type containing ice-top and ocean surface variables
      !! that are used in computing atmosphere-ocean gas fluxes.
    type(FmsCoupler1dBC_type), optional, pointer :: gas_fluxes
      !< is a derived type for exchanging gas or tracer fluxes between the
      !! atmosphere and ocean, defined by the field table, as well as a place holder
      !! of intermediate calculations, such as piston velocities, and parameters
      !! that impact the fluxes.

    !> @parblock
    !! CALL ATMOS_TRACER_FLUX_INIT(), OCEAN_MODEL_FLUX_INIT(), ATMOS_TRACER_FLUX_INIT().
    !! ALSO CALLS FMS_ATMOS_OCEAN_FLUXES_INIT() TO ALLOCATE DERIVED TYPES.
    !! @endparblock
    if (.not.gas_fluxes_initialized) then
      call fms_atmos_ocean_type_fluxes_init( )
      call ocean_model_flux_init( )
      call atmos_tracer_flux_init( )
      call fms_atmos_ocean_fluxes_init(ex_gas_fluxes, ex_gas_fields_atm, ex_gas_fields_ice)
      gas_fluxes_initialized = .true.
    endif

    !> @parblock
    !! SET MODULE LEVEL GAS_FIELDS_ATM, GAS_FIELDS_ICE, AND GAS_FLUXES.
    !! @endparblock
    if (present(gas_fields_atm)) gas_fields_atm => ex_gas_fields_atm
    if (present(gas_fields_ice)) gas_fields_ice => ex_gas_fields_ice
    if (present(gas_fluxes)) gas_fluxes => ex_gas_fluxes

  end subroutine gas_exchange_init

  !#######################################################################
  !> @parblock
  !! Subroutine flux_exchange_init setups derived types and variables, and
  !! initializes modules that will be used in flux exchange.
  !! Ocean_tracer_flux_init is called first to get restart filenames for tracer fluxes
  !! for restart model runs.  Atmos_tracer_flux_init is called last in order
  !! to use tracer values set in ocean_tracer_flux_init.
   !! @endparblock
  subroutine flux_exchange_init ( Time, Atm, Land, Ice, Ocean, Ocean_state,&
       atmos_ice_boundary, land_ice_atmos_boundary, &
       land_ice_boundary, ice_ocean_boundary, ocean_ice_boundary, &
       do_ocean, slow_ice_ocean_pelist, dt_atmos, dt_cpld )

    type(FmsTime_type), intent(in) :: Time
      !< is the model current time
    type(atmos_data_type), intent(inout) :: Atm
      !< is a derived type to specify atm boundary data
    type(land_data_type), intent(in) :: Land
      !< is a derived type to specify land boundary data
    type(ice_data_type), intent(inout) :: Ice
      !< is a derived type to specify ice boundary data
    type(ocean_public_type), intent(inout) :: Ocean
      !< is a derived type to specify ocean boundary data
    type(ocean_state_type), pointer :: Ocean_state
      !< is a pointer to the ocean model's internal state
    type(atmos_ice_boundary_type), intent(inout) :: atmos_ice_boundary
      !< is a derived type holding properties and fluxes passed from atmosphere to ice
    type(land_ice_atmos_boundary_type),intent(inout) :: land_ice_atmos_boundary
      !< is a derived type holding properties and fluxes passed from land and ice to atm
    type(land_ice_boundary_type),  intent(inout) :: land_ice_boundary
      !< is a derived type holding properties and fluxes passed from land to ice
    type(ice_ocean_boundary_type), intent(inout) :: ice_ocean_boundary
      !< is a derived type holding properties and fluxes passed from ice to ocean
    type(ocean_ice_boundary_type), intent(inout) :: ocean_ice_boundary
      !< is a derived type holding properties and fluxes passed from ocean to ice
    logical, intent(in)  :: do_ocean
      !< is a flag indicating whether the ocean component is active
    integer, dimension(:), intent(in) :: slow_ice_ocean_pelist
      !< is an array holding pes for slow ice-ocean exchange
    integer, optional,  intent(in)  :: dt_atmos
      !< is the atmosphere time step in [s]
    integer, optional, intent(in) :: dt_cpld
      !< is the coupled time step in [s]

    character(len=64),  parameter :: grid_file = 'INPUT/grid_spec.nc'
    integer :: ierr, io
    integer :: logunit, uni
    character(len=256) :: errmsg
    integer :: omp_get_num_threads, nthreads

    !> @parblock
    !! CALL FMS_SAT_VAPOR_PRES_INIT.
    !! @endparblock
    call fms_sat_vapor_pres_init()

    !> @parblock
    !! SETUP OPENMP PARAMETERS.
    !! @endparblock
    nthreads = 1
    ! assign nblocks to number of threads.
    !$OMP PARALLEL
    !$  nthreads = omp_get_num_threads()
    !$OMP END PARALLEL
    nblocks = nthreads

    !> @parblock
    !! SET LOGFILE.
    !! @endparblock
    logunit = fms_mpp_stdlog()

    !> @parblock
    !! READ FLUX_EXCHANGE_NML.
    !! @endparblock
    read (fms_mpp_input_nml_file, flux_exchange_nml, iostat=io)
    ierr = fms_check_nml_error (io, 'flux_exchange_nml')

    !> @parblock
    !! WRITE NAMELIST TO LOGFILE.
    !! @endparblock
    call fms_write_version_number (version, tag)
    if( fms_mpp_pe() == fms_mpp_root_pe() )write( logunit, nml=flux_exchange_nml )
    if(nblocks<1) call fms_error_mesg ('flux_exchange_mod',  &
         'flux_exchange_nml nblocks must be positive', FATAL)
    if(nblocks .NE. nthreads) then
       write(errmsg, '(a,i3,a,i3)')'flux_exchange_nml nblocks is set to ', nblocks, &
            ' is different from the default value (number of threads) = ', nthreads
       call fms_error_mesg ('flux_exchange_mod', errmsg, NOTE)
    endif

    !> @parblock
    !! SET MODULE LEVEL DT_ATM AND DT_CPL TIMESTEPS.
    !! @endparblock
    ! required by stock_move, all fluxes used to update stocks will be zero if dt_atmos,
    ! and dt_cpld are absent
    Dt_atm = 0.0
    Dt_cpl = 0.0
    if(present(dt_atmos)) Dt_atm = real(dt_atmos)
    if(present(dt_cpld )) Dt_cpl = real(dt_cpld)

    !> @parblock
    !! GET OCEAN MODEL GRID CELL AREAS FROM GRID_SPEC.
    !! @endparblock
    call fms_xgrid_get_ocean_model_area_elements(Ocean%domain, grid_file)

    !> @parblock
    !! IF ATM%PE, CALL ATM_LAND_ICE_FLUX_EXCHANGE_INIT() AND LAND_ICE_FLUX_EXCHANGE_INIT()
    !! ALSO CHECK ATM_GRID CONSISTENCY WITH PROVIDED GRID_SPEC.
    !! @endparblock
    if( Atm%pe )then
       call fms_mpp_set_current_pelist(Atm%pelist)
       cplClock = fms_mpp_clock_id( 'Land-ice-atm coupler', flags=fms_clock_flag_default, grain=CLOCK_COMPONENT )
       call check_atm_grid(Atm, grid_file)
       call atm_land_ice_flux_exchange_init(Time, Atm, Land, Ice, atmos_ice_boundary, land_ice_atmos_boundary, &
            Dt_atm, Dt_cpl, z_ref_heat, z_ref_mom,  &
            do_area_weighted_flux, do_forecast,  &
            partition_fprec_from_lprec, scale_precip_2d, nblocks, cplClock, &
            ex_gas_fields_atm, ex_gas_fields_ice, ex_gas_fluxes)
       call land_ice_flux_exchange_init(Land, Ice, land_ice_boundary, Dt_cpl, do_runoff, cplClock)
    end if

    !> @parblock
    !! CALL ICE_OCEAN_FLUX_EXCHANGE_INIT().
    !! @endparblock
    call fms_mpp_set_current_pelist()
    call ice_ocean_flux_exchange_init(Time, Ice, Ocean, Ocean_state,ice_ocean_boundary, ocean_ice_boundary, &
         Dt_cpl, debug_stocks, do_area_weighted_flux, ex_gas_fields_ice, ex_gas_fluxes, do_ocean, slow_ice_ocean_pelist)

    !> @parblock
    !! SET DO_INIT TO .FALSE. TO SKIP INITIALIZATION IF FLUX_EXCHANGE_INIT IS CALLED AGAIN.
    !! @endparblock
    do_init = .false.

  end subroutine flux_exchange_init

  !> @parblock
  !! Subroutine flux_check_stocks computes the current stock values for atm, land, ice, and ocean; and
  !! outputs the stock differences with respect to the initial values in the logfile.
   !! @endparblock
  subroutine flux_check_stocks(Time, Atm, Lnd, Ice, Ocn_state)

    type(FmsTime_type), intent(in) :: Time
      !< is the model's current time
    type(atmos_data_type), intent(inout), optional :: Atm
      !< is the atmosphere boundary data type used to compute atmosphere stocks
    type(land_data_type), intent(inout), optional :: Lnd
      !< is the land boundary data type used to compute land stocks
    type(ice_data_type), intent(inout), optional :: Ice
      !< is the ice boundary data type used to compute ice stocks
    type(ocean_state_type), intent(inout), optional, pointer :: Ocn_state
      !< is a pointer to the ocean model's internal state used to compute ocean stocks

    real :: ref_value
    integer :: i

    !> @parblock
    !! FOR WATER, HEAT, AND SALT STOCKS FOR EACH COMPONENT,
    !! GET CURRENT STOCK VALUE AND COMPARE WITH INTEGRATED FLUXES
    !! FOR ATM WATER STOCK.  FOR ATM, INTEGRATE ATM_PRECIP_NEW FOR IMPLICIT EVAPORATION.
    !! @endparblock
    do i = 1, NELEMS !< constant from fms/stock_constants_mod

       if(present(Atm)) then
          ref_value = 0.0
          call Atm_stock_pe(Atm, index=i, value=ref_value)
          if(i==ISTOCK_WATER .and. Atm%pe ) then
             ! decrease the Atm stock by the precip adjustment to reflect the fact that
             ! after an update_atmos_up call, the precip will be that of the future time step.
             ! Thus, the stock call will represent the (explicit ) precip at
             ! the beginning of the preceding time step, and the (implicit) evap at the
             ! end of the preceding time step
             call atm_stock_integrate(Atm, ATM_PRECIP_NEW)
             ref_value = ref_value + ATM_PRECIP_NEW
          endif

          fms_stock_constants_atm_stock(i)%q_now = ref_value
       endif

       if(present(Lnd)) then
          ref_value = 0.0
          call Lnd_stock_pe(Lnd, index=i, value=ref_value)
          fms_stock_constants_lnd_stock(i)%q_now = ref_value
       endif

       if(present(Ice)) then
          ref_value = 0.0
          call Ice_stock_pe(Ice, index=i, value=ref_value)
          fms_stock_constants_ice_stock(i)%q_now = ref_value
       endif

       if(present(Ocn_state)) then
          ref_value = 0.0
          call Ocean_stock_pe(Ocn_state, index=i, value=ref_value)
          fms_stock_constants_ocn_stock(i)%q_now = ref_value
       endif
    enddo

    !> @parblock
    !! PRINT FOR EACH ELEMENT,
    !! S(t): TOTAL STOCK,
    !! S(t)-S(0): CHANGE IN STOCK WITH RESPECT TO INITIAL VALUE,
    !! F(t): CUMULATIVE FLUX INTO COMPONENT FROM OTHER COMPONENTS
    !! F(t) - [S(t)-S(0)]: DIFFERENCE BETWEEN THE FLUXES AND STOCK CHANGE
    !! (S(t)-S(0))/F(t): RELATIVE ERROR
    !! @endparblock
    call fms_stock_constants_stocks_report(Time)


  end subroutine flux_check_stocks

  !#######################################################################
  !> @parblock
  !! Subroutine flux_init_stocks initializes the stock values for the atmosphere,
  !! land, ice, and ocean.  Stocks are the globally integrated total amoun
  !! of conserved quantities such as mass and energy and is used to check conservation.
   !! @endparblock
  subroutine flux_init_stocks(Time, Atm, Lnd, Ice, Ocn_state)
    type(FmsTime_type) , intent(in) :: Time
      !< is the model's current time
    type(atmos_data_type) :: Atm
      !< is a derived type holding atmosphere boundary data
    type(land_data_type) :: Lnd
      !< is a derived type holding land boundary data
    type(ice_data_type) :: Ice
      !< is a derived type holding ice boundary data
    type(ocean_state_type), pointer :: Ocn_state
      !< is a pointer to ocean model's internal state

    integer :: i

    !> @parblock
    !! IF DIVERT_STOCKS_REPORT IS FALSE, OPEN STOCKS OUTPUT FILE TO STDOUT.
    !! IF DIVERT_STOCKS_REPORT IS TRUE, OPEN STOCKS OUTPUT FILE TO "stocks.out".
    !! ONLY THE ROOT PE WILL WRITE TO THE FILE.
    !! @endparblock
    fms_stock_constants_stocks_file=fms_mpp_stdout()
    if(fms_mpp_pe()==fms_mpp_root_pe() .and. divert_stocks_report) then
       open(newunit = fms_stock_constants_stocks_file, file='stocks.out', status='replace', form='formatted')
    endif

    !> @parblock
    !! INITIALIZE WATER, HEAT, AND SALT STOCK VALUES FOR EACH COMPONENT.
    !! FOR ATMOSPHERE, INTEGRATE ATM_PRECIP_NEW TO GET THE INITIAL ISTOCK_WATER.
    !! @endparblock
    do i = 1, NELEMS !from fms/stock_constants_mod
       call Atm_stock_pe(   Atm , index=i, value=fms_stock_constants_atm_stock(i)%q_start)

       if(i==ISTOCK_WATER .and. Atm%pe ) then
          call atm_stock_integrate(Atm, ATM_PRECIP_NEW)
          fms_stock_constants_atm_stock(i)%q_start = fms_stock_constants_atm_stock(i)%q_start + ATM_PRECIP_NEW
       endif

       call Lnd_stock_pe(   Lnd , index=i, value=fms_stock_constants_lnd_stock(i)%q_start)
       call Ice_stock_pe(   Ice , index=i, value=fms_stock_constants_ice_stock(i)%q_start)
       call Ocean_stock_pe( Ocn_state , index=i, value=fms_stock_constants_ocn_stock(i)%q_start)
    enddo

    !> @parblock
    !! INITIALIZE STOCKS IN FMS.
    !! @endparblock
    call fms_stocks_report_init(Time)


  end subroutine flux_init_stocks

  !> @parblock
  !! Subroutine check_atm_grid checks the consistency of the atmosphere grid specified in the model
  !! with the grid specified in the grid_file.
  !! @endparblock
  subroutine check_atm_grid(Atm, grid_file)
    type(atmos_data_type), intent(in) :: Atm
      !< is a derived type holding atmosphere boundary and grid data
    character(len=*), intent(in) :: grid_file
      !< is the path to the grid specification file

    integer :: isg, ieg, jsg, jeg
    integer :: isc, iec, jsc, jec
    integer :: isd, ied, jsd, jed
    integer :: isc2, iec2, jsc2, jec2
    integer :: nxg, nyg, ioff, joff
    integer :: nlon, nlat, siz(4)
    integer :: i, j
    type(FmsMppDomain2D) :: domain2
    real, dimension(:,:), allocatable :: tmpx, tmpy
    real, dimension(:), allocatable :: atmlonb, atmlatb
    character(len=256) :: atm_mosaic_file, tile_file, buffer

    integer, dimension(:), allocatable :: pes
      ! are the current process IDs in the pelist
    type(FmsNetcdfFile_t) :: grid_file_obj, atm_mosaic_file_obj
      ! are the fms2 I/O file objects for the grid specification and atmosphere mosaic files
    type(FmsNetcdfDomainFile_t) :: tile_file_obj
      ! is the fms2 I/O domain file object for the atmosphere mosaic tile
    character(len=20) :: dim_names(2)
      ! are the dimension names for variables in the atmosphere mosaic tile file
    integer :: ppos

    !> @parblock
    !! GET GLOBAL, COMPUTE, AND DATA DOMAIN INDICES AND SIZES FOR THE ATMOSPHERE COMPONENT.
    !! @endparblock
    call fms_mpp_domains_get_global_domain(Atm%domain, isg, ieg, jsg, jeg, xsize=nxg, ysize=nyg)
    call fms_mpp_domains_get_compute_domain(Atm%domain, isc, iec, jsc, jec)
    call fms_mpp_domains_get_data_domain(Atm%domain, isd, ied, jsd, jed)

    !> @parblock
    !! OPEN GRID_FILE.
    !! @endparblock
    allocate(pes(fms_mpp_npes()))
    call fms_mpp_get_current_pelist(pes)

    if ( .not. fms2_io_open_file(grid_file_obj, grid_file, "read", pelist=pes)) then
         call fms_error_mesg ('atm_land_ice_flux_exchange_mod',  &
              & 'Error opening '//trim(grid_file), FATAL)
    endif

    !> @parblock
    !! CHECK GRID SIZES ARE CONSISTENT.
    !! @endparblock
    if(size(Atm%lon_bnd,1) .NE. iec-isc+2 .OR. size(Atm%lon_bnd,2) .NE. jec-jsc+2) then
       call fms_error_mesg ('atm_land_ice_flux_exchange_mod',  &
            'size of Atm%lon_bnd does not match the Atm computational domain', FATAL)
    endif

    ioff = lbound(Atm%lon_bnd,1) - isc
    joff = lbound(Atm%lon_bnd,2) - jsc

    !> @parblock
    !! CHECK LON, LAT, AND GRID CELL AREAS ARE CONSISTENT.
    !! @endparblock
    if(fms2_io_variable_exists(grid_file_obj, "AREA_ATM" ) ) then  ! old grid
       call fms2_io_get_variable_size(grid_file_obj, "AREA_ATM", siz(1:2))
       nlon = siz(1)
       nlat = siz(2)

       if (nlon /= nxg .or. nlat /= nyg) then
          if (fms_mpp_pe()==fms_mpp_root_pe()) then
             print *, 'grid_spec.nc has', nlon, 'longitudes,', nlat, 'latitudes; ', &
                  'atmosphere has', nxg, 'longitudes,', &
                  nyg, 'latitudes (see xba.dat and yba.dat)'
          end if
          call fms_error_mesg ('atm_land_ice_flux_exchange_mod',  &
               'grid_spec.nc incompatible with atmosphere resolution', FATAL)
       end if
       allocate( atmlonb(isg:ieg+1) )
       allocate( atmlatb(jsg:jeg+1) )
       call fms2_io_read_data(grid_file_obj, 'xba', atmlonb)
       call fms2_io_read_data(grid_file_obj, 'yba', atmlatb)

       do i=isc, iec+1
          if(abs(atmlonb(i)-Atm%lon_bnd(i+ioff,jsc+joff)*45.0/atan(1.0))>bound_tol) then
             print *, 'GRID_SPEC/ATMOS LONGITUDE INCONSISTENCY at i= ',i, ': ', &
                  atmlonb(i),  Atm%lon_bnd(i+ioff,jsc+joff)*45.0/atan(1.0)
             call fms_error_mesg ('atm_land_ice_flux_exchange_mod', &
                  'grid_spec.nc incompatible with atmosphere longitudes (see xba.dat and yba.dat)'&
                  , FATAL)
          endif
       enddo
       do j=jsc, jec+1
          if(abs(atmlatb(j)-Atm%lat_bnd(isc+ioff,j+joff)*45.0/atan(1.0))>bound_tol) then
             print *, 'GRID_SPEC/ATMOS LATITUDE INCONSISTENCY at j= ',j, ': ', &
                  atmlatb(j),  Atm%lat_bnd(isc+ioff, j+joff)*45.0/atan(1.0)
             call fms_error_mesg ('atm_land_ice_flux_exchange_mod', &
                  'grid_spec.nc incompatible with atmosphere latitudes (see xba.dat and yba.dat)'&
                  , FATAL)
          endif
       enddo
       deallocate(atmlonb, atmlatb)
    else if(fms2_io_variable_exists(grid_file_obj, "atm_mosaic_file" ) ) then  ! mosaic grid file.
       call fms2_io_read_data(grid_file_obj, 'atm_mosaic_file', atm_mosaic_file)

       if ( .not. fms2_io_open_file(atm_mosaic_file_obj, "INPUT/"//trim(atm_mosaic_file)//"", "read", pelist=pes)) then
           call fms_error_mesg ('atm_land_ice_flux_exchange_mod',  &
              & 'Error opening '//trim(atm_mosaic_file), FATAL)
       endif

       call fms2_io_read_data(atm_mosaic_file_obj, "gridfiles", buffer, corner=1)

       !< Remove the .tile from the filename to get basename
       ppos = index(trim(buffer),".tile")
       if ( ppos > 0 ) then
          tile_file = buffer(1:ppos-1)//".nc"
       else
          tile_file = buffer
       endif

       call fms2_io_close_file(atm_mosaic_file_obj)

       call fms_mpp_domains_copy_domain(Atm%domain, domain2)
       call fms_mpp_domains_create_super_grid_domain(domain2)
       call fms_mpp_domains_define_io_domain  (domain2, (/1,1/))

       call fms_mpp_domains_get_compute_domain(domain2, isc2, iec2, jsc2, jec2)

       if(isc2 .NE. 2*isc-1 .OR. iec2 .NE. 2*iec+1 .OR. jsc2 .NE. 2*jsc-1 .OR. jec2 .NE. 2*jec+1) then
          call fms_mpp_error(FATAL, 'atm_land_ice_flux_exchange_mod: supergrid domain is not set properly')
       endif

       !< This is will open the correct atm_mosaic_file for the current tile, i.e "C96_grid.tile1.nc"
       if ( .not. fms2_io_open_file(tile_file_obj, "INPUT/"//trim(tile_file)//"", "read", domain2)) then
          call fms_error_mesg ('atm_land_ice_flux_exchange_mod',  &
              & 'Error opening '//trim(tile_file), FATAL)
       endif

       call fms2_io_get_variable_size(tile_file_obj, 'area', siz(1:2))
       nlon = siz(1); nlat = siz(2)
       if( mod(nlon,2) .NE. 0) call fms_mpp_error(FATAL,  &
            'atm_land_ice_flux_exchange_mod: atmos supergrid longitude size can not be divided by 2')
       if( mod(nlat,2) .NE. 0) call fms_mpp_error(FATAL,  &
            'atm_land_ice_flux_exchange_mod: atmos supergrid latitude size can not be divided by 2')
       nlon = nlon/2
       nlat = nlat/2
       if (nlon /= nxg .or. nlat /= nyg) then
          if (fms_mpp_pe()==fms_mpp_root_pe()) then
             print *, 'atmosphere mosaic tile has', nlon, 'longitudes,', nlat, 'latitudes; ', &
                  'atmosphere has', nxg, 'longitudes,', nyg, 'latitudes'
          end if
          call fms_error_mesg ('atm_land_ice_flux_exchange_mod',  &
               'atmosphere mosaic tile grid file incompatible with atmosphere resolution', FATAL)
       end if

       allocate(tmpx(isc2:iec2,jsc2:jec2), tmpy(isc2:iec2,jsc2:jec2) )

       !< Register the dimension of the variables "x" and "y" in the atm_mosaic_file
       call fms2_io_get_variable_dimension_names(tile_file_obj, "x", dim_names)
       call fms2_io_register_axis(tile_file_obj, dim_names(1), "x")
       call fms2_io_register_axis(tile_file_obj, dim_names(2), "y")
       call fms2_io_register_field(tile_file_obj, "x", "double", dim_names)
       call fms2_io_register_field(tile_file_obj, "y", "double", dim_names)

       !< Read the variables "x" and "y" as domain decomposed variables from the atm_moasic_file
       call fms2_io_read_data( tile_file_obj, 'x', tmpx)
       call fms2_io_read_data( tile_file_obj, 'y', tmpy)

       call fms2_io_close_file(tile_file_obj)

       call fms_mpp_domains_deallocate_domain(domain2)

       do j = jsc, jec+1
          do i = isc, iec+1
             if (abs(tmpx(2*i-1,2*j-1)-Atm%lon_bnd(i+ioff,j+joff)*45.0/atan(1.0))>bound_tol) then
                print *, 'GRID_SPEC/ATMOS LONGITUDE INCONSISTENCY at i= ',i, ', j= ', j, ': ', &
                     tmpx(2*i-1,2*j-1),  Atm%lon_bnd(i+ioff,j+joff)*45.0/atan(1.0)
                call fms_error_mesg ('atm_land_ice_flux_exchange_mod', &
                     'grid_spec.nc incompatible with atmosphere longitudes (see '//trim(tile_file)//')'&
                     ,FATAL)
             end if
             if (abs(tmpy(2*i-1,2*j-1)-Atm%lat_bnd(i+ioff,j+joff)*45.0/atan(1.0))>bound_tol) then
                print *, 'GRID_SPEC/ATMOS LATITUDE INCONSISTENCY at i= ',i, ', j= ', j, ': ', &
                     tmpy(2*i-1,2*j-1),  Atm%lat_bnd(i+ioff,j+joff)*45.0/atan(1.0)
                call fms_error_mesg ('atm_land_ice_flux_exchange_mod', &
                     'grid_spec.nc incompatible with atmosphere latitudes (see '//trim(tile_file)//')'&
                     ,FATAL)
             end if
          end do
       end do
       deallocate(tmpx, tmpy)
    else
       call fms_mpp_error(FATAL, &
            'atm_land_ice_flux_exchange_mod: both AREA_ATMxOCN and ocn_mosaic_file does not exist in '//trim(grid_file))
    end if

    call fms2_io_close_file(grid_file_obj)
  end subroutine check_atm_grid

end module flux_exchange_mod
!> @}
