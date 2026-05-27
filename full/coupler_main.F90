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
!> @defgroup coupler_main coupler_main
!! @ingroup FMSCoupler
!!
!! @brief Main driver program for the fully coupled (atmosphere, land, sea ice,
!! and ocean) GFDL climate models.
!!
!! Please see the [**main page**](index.html) for additional information.
!!
!! @author Bruce Wyman <Bruce.Wyman@noaa.gov>
!! @author V. Balaji <V.Balaji@noaa.gov>

!> @file
!! @brief Main driver program for the fully coupled GFDL climate model.
!!
!! \parblock
!! The program coupler_main couples the atmosphere, ocean, land, and sea-ice components,
!! each modeled on independent grids; and also calls the drives for each component.
!!
!! The model start time (Time_start) is determined in coupler_init by the following:
!! If date_init exists in the diag_table, it is used to set the model start time.
!! If date_init is not found in the diag_table, the model start time is set as below:
!! If INPUT/coupler.res exists, the start date and the calendar values are
!! read in to set the model start date and the calendar type.  These values can be 
!! overwritten if force_date_from_namelist = .true. and current_date with calendar_type 
!! is defined in coupler_nml.  If date_init is not found in the diag_table and INPUT/coupler.res
!! does not exist, the start date is taken from current_date and calendar in coupler_nml.  
!!
!! There are two nested time integration loops:
!! - Slow (coupled) loop (coupled_timestep_loop, nc = 1 … num_cpld_calls):
!!   advances the ocean and slow sea-ice by one coupled timestep dt_cpld.
!!   Ocean–ice fluxes are exchanged once per iteration.
!! - Fast (atmospheric) loop** (fast_integration_loop, na = 1 … num_atmos_calls):
!!   advances the atmosphere, land surface, and fast sea-ice by one atmospheric
!!   timestep dt_atmos.  The fast loop is nested inside the slow loop.
!!
!! Heat and moisture are exchanged between the atmosphere and the surface
!! (land/ice) using an tridiagonal scheme for implicit vertical diffusion in the fast loop:
!! 1. coupler_update_atmos_model_down — forward (downward) sweep from atmospheric top to surface.
!! 2. coupler_flux_down_from_atmos — transfers forward-elimination coefficients to land/ice.
!! 3. Land and ice fast updates (coupler_update_land_model_fast,
!!    coupler_update_ice_model_fast) — compute new surface temperatures.
!! 4. coupler_flux_up_to_atmos — applies implicit surface flux corrections using
!!    updated surface temperatures.
!! 5. coupler_update_atmos_model_up — back-substitution (upward sweep), convection,
!!    and large-scale condensation.
!!
!! When `do_concurrent_radiation = .true.`, the atmospheric dynamics/physics and
!! radiation updates run simultaneously with OpenMP threading (see `atmos_nthreads`
!! and `radiation_nthreads` namelist variables).

!! Ocean fluxes are exchanged explicitly (one coupling step behind).  All fluxes
!! reaching the ocean — including atmospheric fluxes and land runoff — are passed
!! through the sea-ice model via Ice_ocean_boundary.
!!
!! Sea-ice physics is split into two timescales that can run on different MPI PE sets:
!! - Fast ice on Ice%fast_ice_pe: thermodynamics and surface-flux coupling at the
!!   atmospheric timestep.  Fast ice always runs on a subset of the atmosphere PEs
!! - Slow ice on Ice%slow_ice_pe: dynamics, freezing/melting, and transport at the
!!   coupled (ocean) timestep.  The placement of the slow ice PEs depends on
!!   the slow_ice_with_ocean namelist variable:
!!   - slow_ice_with_ocean = .false. (default): slow and fast ice share the same PEs
!!   - slow_ice_with_ocean = .true.: slow ice runs on the ocean PEs.  In this case,
!!     `Ice%pelist` is the union of the fast (atmos) and slow (ocean) PE sets.
!! The flag concurrent_ice = .true. runs the fast ice and slow ice processes concurrently
!! and requires slow_ice_with_ocean = .true.  
!! The flag combine_ice_and_ocean = .true. advances the slow ice and ocean processes together 
!! on the ocean PEs.  The flags concurent_ice and slow_ice_with_ocean must be .true. to use combine_ice_and_ocean.
!!
!! Full coupling is configured through three namelists:
!! - @ref coupler_config "coupler_nml"
!! - @ref flux_exchange_conf "flux_exchange_nml"
!! - @ref surface_flux_config "surface_flux_nml"
!!
!! Pseudocode:
!! call fms_diag_init(...)        ! open diagnostic output
!! call fms_tracer_manager_init() ! register tracers
!! call gas_exchange_init(...)    ! register air-sea gas/tracer BCs
!! call flux_exchange_init(...)   ! build atm-land-ice exchange grids
!! call atmos_model_init(...)     ! initialize atmosphere
!! call land_model_init(...)      ! initialize land
!! call ice_model_init(...)       ! initialize sea ice (fast + slow)
!! call ocean_model_init(...)     ! initialize ocean
!!
!! do nc = 1, num_cpld_calls
!!
!!   ! Redistribute ocean surface state onto ice grid
!!   call flux_ocean_to_ice(...)

!!   ! If slow_ice_pe: override ocean-ice BCs, unpack into Ice type
!!   call flux_ocean_to_ice_finish(...)
!!   call unpack_ocean_ice_boundary(...)
!!
!!   ! Exchange slow-ice state to fast-ice data structures
!!   call exchange_slow_to_fast_ice(...)
!!
!!   ! Prepare ice surface fields (albedo, T, etc.) for atmos surface flux calc
!!   call set_ice_surface_fields(...)
!!
!!   do na = 1, num_atmos_calls
!!
!!     ! Copy Atm%tr_bot → Atm%fields for gas-exchange tracers
!!     call atmos_tracer_driver_gather_data(Atm%fields, Atm%tr_bot)
!!
!!     ! Compute surface exchange coefficients and turbulent fluxes on the
!!     ! atm-land-ice exchange grid
!!     call sfc_boundary_layer(...)
!!
!!     ! Atmosphere dynamical core (FV3)
!!     call update_atmos_model_dynamics(...)
!!
!!     ! Radiation (sequential, or concurrent on a separate OMP team)
!!     call update_atmos_model_radiation(...)
!!
!!     ! Forward (downward) sweep of the implicit tridiagonal diffusion
!!     call update_atmos_model_down(...)
!!
!!     ! Apply implicit atm diffusion correction; pass updated surface fluxes
!!     ! (heat, moisture, momentum) to land and ice boundary types
!!     call flux_down_from_atmos(...)
!!
!!     ! Fast land physics (hydrology, canopy, soil temperature)
!!     call update_land_model_fast(...)
!!
!!     ! Fast ice thermodynamics (surface energy balance, melt ponds)
!!     call update_ice_model_fast(...)
!!
!!     ! Recompute surface fluxes using updated land/ice surface temperatures
!!     call flux_up_to_atmos(...)
!!
!!     ! Back-substitution (upward) sweep; convection; large-scale condensation
!!     call update_atmos_model_up(...)
!!
!!     ! Remap atmosphere gas/tracer fields onto exchange grid;
!!     ! compute air-sea deposition fluxes; deallocate exchange grid arrays
!!     call flux_atmos_to_ocean(...)
!!     call flux_ex_arrays_dealloc()
!!
!!     ! Advance atmos diagnostics and tracer state
!!     call update_atmos_model_state(...)
!!
!!   end do  ! fast loop
!!
!!   ! Slow land physics (routing, carbon, DGVM)
!!   call update_land_model_slow(...)
!!
!!   ! Interpolate land runoff and calving onto ice grid
!!   call flux_land_to_ice(...)
!!
!!   ! Reset fast-ice accumulators; copy Land_ice_boundary into ice internals
!!   call ice_model_fast_cleanup(...)
!!   call unpack_land_ice_boundary(...)
!!
!!   ! Exchange fast-ice averages to slow-ice side
!!   call exchange_fast_to_slow_ice(...)
!!
!!   ! Slow ice physics (dynamics, freezing/melting, transport)
!!   call update_ice_model_slow(...)
!!
!!   ! Bookkeep ice-to-ocean flux stocks (water, heat, salt)
!!   call flux_ice_to_ocean_stocks(...)
!!
!!   ! Interpolate ice-bottom fluxes onto ocean grid -> Ice_ocean_boundary
!!   call flux_ice_to_ocean(...)
!!   ! Override/diagnose Ice_ocean_boundary fields; send to diag_manager
!!   call flux_ice_to_ocean_finish(...)
!!
!!   ! Advance ocean state by dt_cpld using Ice_ocean_boundary forcing
!!   call update_ocean_model(...)
!!   ! (or: call update_slow_ice_and_ocean(...) if combined_ice_and_ocean)
!!
!!   ! Bookkeep ocean stocks from ice-ocean flux transfer
!!   call flux_ocean_from_ice_stocks(...)
!!
!!   ! Flush diagnostic send-data buffer for this coupled step
!!   call fms_diag_send_complete(Time_step_cpld)
!!
!! end do
!!
!! call coupler_restart(...) ! write coupler.res and component restart files
!! call fms_diag_end(...)    ! flush and close diagnostic output
!! \endparblock


!> @ingroup coupler_main
program coupler_main
  !--- F90 module for OpenMP
  use omp_lib
  use FMS
  use full_coupler_mod

  implicit none

  !> model defined types.
  !! Targets to pointers in coupler_components_obj
 
  !> Datatype holding instantaneous atm model state at current timestep
  type (atmos_data_type), target :: Atm
  !> datatype holding instantaneous land model state at current timestep
  type (land_data_type), target :: Land 
  !> datatype holding instantaneous ice model state at current timestep
  type (ice_data_type), target :: Ice 
  ! allow members of ocean type to be aliased (ap)
  !> datatype holding ocean model state at current timestep.  Target for the Ocean_state pointer
  type (ocean_public_type), target  :: Ocean
  !> Alias pointer to Ocean datatype
  type (ocean_state_type),  pointer :: Ocean_state => NULL()

  !> datatype holding data to exchange between atmos and land
  type(atmos_land_boundary_type), target :: Atmos_land_boundary
  !> datatype holding data to exchange between atmos and sea ice  
  type(atmos_ice_boundary_type), target  :: Atmos_ice_boundary
  !> datatype holding data to exchange between land and ice to atmos
  type(land_ice_atmos_boundary_type), target  :: Land_ice_atmos_boundary
  !> datatype holding data to exchange between land and ice
  type(land_ice_boundary_type), target  :: Land_ice_boundary
  !> datatype holding data to exchange from ice and ocean
  type(ice_ocean_boundary_type), target :: Ice_ocean_boundary
  !> datatype holding data to exchange from ocean to ice
  type(ocean_ice_boundary_type), target :: Ocean_ice_boundary
  !> Pointer alias to ice_ocean_driver_type containing control parameters to combined ice-ocean-driver
  type(ice_ocean_driver_type), pointer  :: ice_ocean_driver_CS => NULL()

  !> current model time
  type(FmsTime_type) :: Time
  !> timestep used in the fast timestepping loop 
  type(FmsTime_type) :: Time_step_atmos 
  !> timestep used in the slow timestepping loop
  type(FmsTime_type) :: Time_step_cpld 
  !> time tracked in the fast timestepping loop
  type(FmsTime_type) :: Time_atmos
  !> time tracked for the ocean model
  type(FmsTime_type) :: Time_ocean
  !> time tracked for lag_fluxes from ice to ocean
  type(FmsTime_type) :: Time_flux_ice_to_ocean 
  !> time tracked for flux exchange from ocean to ice
  type(FmsTime_type) :: Time_flux_ocean_to_ice 

  !> number of timesteps in the fast-integration loop
  integer :: num_atmos_calls 
  !> do loop counter in the fast-integration loop
  integer :: na
  !> number of timesteps in the slow-integration loop 
  integer :: num_cpld_calls 
  !> do loop counter in the slow-integration loop
  integer :: nc
  !> Accumulated count of fast (atmospheric) timestep iterations across the
  !! entire run so far; equals `(nc-1)*num_atmos_calls + na` and is used as a
  !! unique step index for checksum labels and diagnostic timestamps.
  integer :: current_timestep

  !> fms2_io file type to read/write data on a decomposed domain
  type(FmsNetcdfDomainFile_t), dimension(:), pointer :: Ice_bc_restart => NULL()
  !> fms2_io file type to read/write data on a decomposed domain  
  type(FmsNetcdfDomainFile_t), dimension(:), pointer :: Ocn_bc_restart => NULL()

  !> the next timepoint to write intermediate restarts
  type(FmsTime_type) :: Time_restart 
  !> model start time
  type(FmsTime_type) :: Time_start 
  !> model end time
  type(FmsTime_type) :: Time_end
  !> last timepoint when intermediate restarts were written
  type(FmsTime_type) :: Time_restart_current 

  !> derived type holding clock ids for clocks used in runtime profiling and debugging  
  type(coupler_clock_type) :: coupler_clocks 
  !> object containing pointers to all the model component derived types.  Used primarily in coupler_chksum_type
  type(coupler_components_type), target :: coupler_components_obj
  !> convenient object holding pointers to model component derived types.  Used when generating CHECKSUMS
  type(coupler_chksum_type) :: coupler_chksum_obj 

  !> MPI PE list for ensemble runs; shape `(npes_per_member, num_ensemble_members)`
  integer, allocatable :: ensemble_pelist(:, :)
  !> Combined MPI PE list of the slow sea-ice and ocean PEs; used to set the
  !! current pelist when both components need to communicate together
  integer, allocatable :: slow_ice_ocean_pelist(:)
  !> Total number of concurrent OpenMP thread teams used within the fast
  !! integration loop.  Equals `atmos_nthreads + radiation_nthreads` when
  !! `do_concurrent_radiation = .true.`; otherwise 1 (serial radiation).
  integer :: conc_nthreads = 1
  !> Scratch variable that records the OMP wall-clock start time (via
  !! `omp_get_wtime`) at the beginning of each thread team's work, used to
  !! accumulate elapsed time into `omp_sec`.
  real :: dsec
  !> Accumulated OMP wall-clock time (seconds) within the current coupled
  !! timestep for each concurrent thread team:
  !!   - `omp_sec(1)` — atmosphere dynamics/physics team
  !!   - `omp_sec(2)` — concurrent radiation team
  !! Reset to zero after each call to `coupler_summarize_timestep`.
  real :: omp_sec(2)=0.0
  !> Accumulated OMP load-imbalance time (seconds) within the current coupled
  !! timestep for each thread team:
  !!   - `imb_sec(1)` — atmosphere team idle time waiting for radiation
  !!   - `imb_sec(2)` — radiation team idle time waiting for atmosphere
  !! Reset to zero after each call to `coupler_summarize_timestep`.
  real :: imb_sec(2)=0.0

  
  !> INITIALIZE FMS MPP_MOD.  MPP_INIT MUST BE INITIALIZED FIRST
  !! BEFORE CREATING A CLOCK WITH MPP_CLOCK_ID
  call fms_mpp_init()

  
  !> START CLOCK TO MEASURE INITIALIZATION ROUTINE
  coupler_clocks%initialization = fms_mpp_clock_id( 'Initialization' )
  call fms_mpp_clock_begin(coupler_clocks%initialization)

  
  !> INITIALIZE FMS
  !{
  call fms_init
  call fmsconstants_init
  call fms_affinity_init
  !}

  
  !> INITIALIZE COUPLER PROGRAM VARIABLES
  call coupler_init(Atm, Ocean, Land, Ice, Ocean_state, Atmos_land_boundary, Atmos_ice_boundary, &
    Ocean_ice_boundary, Ice_ocean_boundary, Land_ice_atmos_boundary, Land_ice_boundary,          &
    Ice_ocean_driver_CS, Ice_bc_restart, Ocn_bc_restart, ensemble_pelist, slow_ice_ocean_pelist, &
    conc_nthreads, coupler_clocks, coupler_components_obj, coupler_chksum_obj, &
    Time_step_cpld, Time_step_atmos, Time_atmos, Time_ocean, num_cpld_calls,   &
    num_atmos_calls, Time, Time_start, Time_end, Time_restart, Time_restart_current)

  
  !> IF DO_CHKSUM, COMPUTE CHECKSUM OF ATM, LAND, AND ICE FIELDS
  if (do_chksum) call coupler_chksum_obj%get_coupler_chksums('coupler_init+', 0)

  
  !> SYNCHRONIZE ALL PES
  call fms_mpp_set_current_pelist()

  
  !> END CLOCK TO MEASURE INITIALIZATION ROUTINE
  call fms_mpp_clock_end(coupler_clocks%initialization)

  
  !> START CLOCK TO MEASURE MAIN LOOP
  call fms_mpp_clock_begin(coupler_clocks%main)


  !> IF CHECK_STOCKS >= 0 AND DO_FLUX, COMPUTE FLUX STOCKS
  if (check_stocks >= 0 .and. do_flux) &
       call coupler_flux_init_finish_stocks(Time, Atm, Land, Ice, Ocean_state, coupler_clocks, init_stocks=.True.)

  
  !> START OCEAN/SLOW-ICE INTEGRATION LOOP
  coupled_timestep_loop : do nc = 1, num_cpld_calls

     
    !> IF DO_CHKSUM, COMPUTE CHECKSUMS OF ATM, LAND, ICE, AND OCEAN FIELDS
    if (do_chksum) then
      call coupler_chksum_obj%get_coupler_chksums('top_of_coupled_loop+', nc)
      call coupler_chksum_obj%get_atmos_ice_land_ocean_chksums('MAIN_LOOP-', nc)
    end if


    
    !> FOR SLOW_ICE_PES AND OCEAN PES, CALL FLUX_OCEAN_TO_ICE, REDISTRIBUTE FLUXES FROM OCEAN TO ICE
    !! AND STORE IN OCEAN_ICE_BOUNDARY TYPE.  IF USE_LAG_FLUXES IS TRUE, REDISTRIBUTE FLUXES AT THE
    !! BOTTOM OF THE ICE TO THE OCEAN MODEL GRID
    ! Calls to flux_ocean_to_ice and flux_ice_to_ocean are all PE communication
    ! points when running concurrently. The calls are placed next to each other in
    ! concurrent mode to avoid multiple synchronizations within the main loop.
    ! With concurrent_ice, these only occur on the ocean PEs.
    !{
    if (Ice%slow_ice_PE .or. Ocean%is_ocean_pe) then
      ! Redistribute quantities from Ocean to Ocean_ice_boundary
      call coupler_flux_ocean_to_ice(Ocean, Ice, Ocean_ice_boundary, coupler_clocks, slow_ice_ocean_pelist)
      Time_flux_ocean_to_ice = Time
      ! Update Ice_ocean_boundary; the first iteration is supplied by restarts
      if(use_lag_fluxes) then
        call coupler_flux_ice_to_ocean(Ice, Ocean, Ice_ocean_boundary, coupler_clocks)
        Time_flux_ice_to_ocean = Time
      end if
   end if
   !}


   !> IF DO_CHKSUM, COMPUTE CHECKSUMS OF ATM, LAND, ICE, AND OCEAN FIELDS
   !{
    if (do_chksum) then
      call coupler_chksum_obj%get_coupler_chksums('flux_ocn2ice+', nc)
      call coupler_chksum_obj%get_atmos_ice_land_ocean_chksums('flux_ocn2ice+', nc)
   end if
   !}
   

   !> IF CHECK_STOCKS > 0 AND DO_FLUX IS TRUE, COMPUTE FLUX STOCKS
   !{
   ! needs to sit here rather than at the end of the coupler loop.
   if (check_stocks > 0 .and. do_flux) &
        call coupler_flux_check_stocks(nc, Time, Atm, Land, Ice, Ocean_state, coupler_clocks)
   !}


   !> IF DO_ICE, AND CURRENT PE IS ICE PE,
   !! (1) COPY INFORMATION FROM OCEAN_ICE_BOUNDARY INTO SLOW PART OF ICE DATATYPE
   !! (2) COPY INFORMATION FROM THE FAST PART OF SEA-ICE TO SLOW PART OF SEA ICE
   !! (3) PREPARE ICE SURFACE STATE FOR ATMOSPHERE FAST PHYSICS
   !{
   if (do_ice .and. Ice%pe) then
      if (Ice%slow_ice_pe) call coupler_unpack_ocean_ice_boundary(nc, Time_flux_ocean_to_ice, Ice, Ocean_ice_boundary,&
           coupler_clocks, coupler_chksum_obj)

      ! This could be a point where the model is serialized if the fast and
      ! slow ice are on different PEs.  call fms_mpp_set_current_pelist(Ice%pelist)
      ! is called if(.not.Ice%shared_slow_fast_PEs)
      call coupler_exchange_slow_to_fast_ice(Ice, coupler_clocks)
      ! This call occurs all ice PEs.
      if (concurrent_ice) call coupler_exchange_fast_to_slow_ice(Ice, coupler_clocks)
      ! call fms_mpp_set_current_pelist(Ice%pelist) is called if(.not.Ice%shared_slow_fast_PEs)
      if (Ice%fast_ice_pe) call coupler_set_ice_surface_fields(Ice, coupler_clocks)
   endif
   !}
   

   !> IF PE IS ATMP%PE
    atm_pe_block : if (Atm%pe) then

      !> SYNCHRONIZE ATM PES
      if (.NOT.(do_ice.and.Ice%pe) .OR. (ice_npes.NE.atmos_npes)) call fms_mpp_set_current_pelist(Atm%pelist)

      !> CALL CHECKSUM FOR ATMOS ICE LAND FIELDS
      if(do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('set_ice_surface+', nc)

      !> START CLOCK FOR PROFILING ATM
      call fms_mpp_clock_begin(coupler_clocks%atm)

      !> IF DO_FLUX, GENERATE THE SURFACE EXCHANGE GRID TO EXCHANGE FLUXES BETWEEN LAND AND ICE
      if (do_flux) call coupler_generate_sfc_xgrid(Land, Ice, coupler_clocks)

      !> SEND ICE MASK TO DIAG_MANAGER BUFFER
      call send_ice_mask_sic(Time)

      !> START CLOCK TO PROFILE FAST LOOP INTEGRATION
      call fms_mpp_clock_begin(coupler_clocks%atmos_loop)


      !> START ATMOS/FAST-LAND/FAST-ICE INTEGRATION LOOP
      fast_integration_loop : do na = 1, num_atmos_calls

        !> INCREMENT TIME_ATMOS BY TIME_STEP_ATMOS
        Time_atmos = Time_atmos + Time_step_atmos

        !> INCREMENT CURRENT_TIMESTEP
        current_timestep = (nc-1)*num_atmos_calls+na

        !> IF DO_CHECKSUM, COMPUTE CHECKSUMS FOR ATMOS, ICE, LAND FIELDS
        if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('top_of_atmos_loop-', current_timestep)

        !> IF DO_ATMOS, COPY %ATM%TR_BOT TO ATM%FIELDS FOR GASES
        if (do_atmos) call coupler_atmos_tracer_driver_gather_data(Atm, coupler_clocks)

        !> IF DO_FLUX, COMPUTE THE FLUXES BETWEN MODEL COMPONENTS
        if (do_flux) call coupler_sfc_boundary_layer(Atm, Land, Ice, Land_ice_atmos_boundary, &
             Time_atmos, current_timestep, coupler_chksum_obj, coupler_clocks)


!$OMP   PARALLEL  &
!$OMP&    NUM_THREADS(conc_nthreads)  &
!$OMP&    DEFAULT(NONE)  &
!$OMP&    PRIVATE(conc_nthreads) &
!$OMP&    SHARED(atmos_nthreads, radiation_nthreads, nc, na, num_atmos_calls, atmos_npes, land_npes, ice_npes) &
!$OMP&    SHARED(Time_atmos, Atm, Land, Ice, Land_ice_atmos_boundary, Atmos_land_boundary, Atmos_ice_boundary) &
!$OMP&    SHARED(Ocean_ice_boundary) &
!$OMP&    SHARED(do_debug, do_flux, do_chksum, do_atmos, do_land, do_ice, do_concurrent_radiation, omp_sec, imb_sec) &
!$OMP&    SHARED(coupler_clocks, current_timestep, coupler_chksum_obj)
!$      if (omp_get_thread_num() == 0) then
!$OMP     PARALLEL &
!$OMP&      NUM_THREADS(1) &
!$OMP&      DEFAULT(NONE) &
!$OMP&      PRIVATE(dsec) &
!$OMP&      SHARED(atmos_nthreads, radiation_nthreads, nc, na, num_atmos_calls, atmos_npes, land_npes, ice_npes) &
!$OMP&      SHARED(Time_atmos, Atm, Land, Ice, Land_ice_atmos_boundary, Atmos_land_boundary, Atmos_ice_boundary) &
!$OMP&      SHARED(Ocean_ice_boundary) &
!$OMP&      SHARED(do_debug, do_flux, do_chksum, do_atmos, do_land, do_ice, do_concurrent_radiation, omp_sec, imb_sec) &
!$OMP&      SHARED(coupler_clocks, current_timestep, coupler_chksum_obj)
!$        call omp_set_num_threads(atmos_nthreads)
!$        dsec=omp_get_wtime()

        
        !>  START CLOCK TO MEASURE DO_CONCURRENT_RADIATION
        if (do_concurrent_radiation) call fms_mpp_clock_begin(coupler_clocks%concurrent_atmos)
        
        !> IF DO_ATMOS, CALL ATMOSPHERE DRIVER, CALL FV DYNAMICAL CORE DRIVER
        if (do_atmos) &
             call coupler_update_atmos_model_dynamics(Atm, current_timestep, coupler_chksum_obj, coupler_clocks)
        
        !> IF NOT DO_CONCURRENT_RADIATION, CALL THE RADIATION_DRIVER
        if (.not.do_concurrent_radiation) call coupler_update_atmos_model_radiation(Atm, Land_ice_atmos_boundary, &
             coupler_clocks, current_timestep, coupler_chksum_obj)
        
        !> IF DO_ATMOS, CALL PHYSICS_DRIVER_DOWN TO COMPUTE ATMOSPEHRIC TENDENCIES FOR DYNAMICS,
        !! RADIATION, VERTICAL DIFFUSION OF MOMENTUM, TRACERS, AND HEAT/MOISTURE.
        !! FOR HEAT/MOISTURE, ONLY THE DOWNWARD SWEEP OF THE TRIDONAL ELIMINATION IS PERFORMED 
        if (do_atmos) call coupler_update_atmos_model_down(Atm, Land_ice_atmos_boundary, &
             current_timestep, coupler_chksum_obj, coupler_clocks)
        
        !> IF DO_FLUX, CORRECT FOR IMPLICIT TREATMENT OF ATMOSPHERIC DIFFUSIVE FLUXES IN FLUX EXCHANGE
        !! FROM ATM TO LAND AND ICE. CHECKSUMS ARE COMPUTED IF DO_CHKSUM=.TRUE.
        if (do_flux) call coupler_flux_down_from_atmos(Atm, Land, Ice, Land_ice_atmos_boundary, Atmos_land_boundary, &
             Atmos_ice_boundary, Time_atmos, current_timestep, coupler_clocks, coupler_chksum_obj)
        
        
        !> IF DO_LAND, CALL LAND DYNAMICS DRIVER FOR PROCESSES OCCURING AT FAST TIMESCALE 
        if (do_land .AND. land%pe) call coupler_update_land_model_fast(Land, Atmos_land_boundary, Atm%pelist, &
             current_timestep, coupler_chksum_obj, coupler_clocks)
        
        !> IF DO_ICE, RECORD FLUXES IN ICE TYPE AND CALCULATE ICE TEMPERATURE 
        if (do_ice .AND. Ice%fast_ice_pe) call coupler_update_ice_model_fast(Ice, Atmos_ice_boundary, Atm%pelist, &
             current_timestep, coupler_chksum_obj, coupler_clocks)
        
        
        !> IF DO_FLUX, CORRECT FOR FLUXES TO TAKE INTO ACCOUNT THE NEW SURFACE TEMPERATURES IN LAND AND ICE MODELS
        if (do_flux) call coupler_flux_up_to_atmos(Land, Ice, Land_ice_atmos_boundary, Atmos_land_boundary, &
             Atmos_ice_boundary, Time_atmos, current_timestep, coupler_chksum_obj, coupler_clocks)

        
        !> IF DO_ATMOS, COMPUTE UPWARD VERTICAL DIFFUSION OF HEAT/MOISTURE AND MOISTURE PROCESSES:
        !! COMPUTE UPWARD SWEEP OF THE TRIDONAL ELIMINATION FOR HEAT/MOISTURE AND COMPUTE
        !! THE CONVECTIVE AND LARGE-SCALE TENDENCIES
        if (do_atmos) call coupler_update_atmos_model_up(Atm, Land_ice_atmos_boundary, current_timestep, &
             coupler_chksum_obj, coupler_clocks)

        
        !> IF DO_FLUX, COMPUTES DEPOSITION GAS FLUXES BETWEEN ATMOSPHERE AND OCEAN
        if (do_flux) call coupler_flux_atmos_to_ocean(Atm, Atmos_ice_boundary, Ice, Time_atmos)

        
        !> IF DO_CONCURRENT_RADIATION, END CLOCK TO MEASURE CONCURRENT_ATMOS
        if (do_concurrent_radiation) call fms_mpp_clock_end(coupler_clocks%concurrent_atmos)

        
        !> IF DO_CONCURRENT_RADIATION, CALL THE RADIATION DRIVER
        !{
!$        omp_sec(1) = omp_sec(1) + (omp_get_wtime() - dsec)
!$OMP END PARALLEL
!$      endif
!$      if (omp_get_thread_num() == max(0,omp_get_num_threads()-1)) then
          !      ---- atmosphere radiation ----
        if (do_concurrent_radiation) then
!$OMP PARALLEL &
!$OMP&      NUM_THREADS(1) &
!$OMP&      DEFAULT(NONE) &
!$OMP&      PRIVATE(dsec) &
!$OMP&      SHARED(Atm, Land, Ice, Land_ice_atmos_boundary, Atmos_ice_boundary, Ocean_ice_boundary,Atmos_land_boundary)&
!$OMP&      SHARED(do_chksum, do_debug, omp_sec, num_atmos_calls, na, radiation_nthreads) &
!$OMP&      SHARED(coupler_clocks)
!$          call omp_set_num_threads(radiation_nthreads)
!$          dsec=omp_get_wtime()
           call coupler_update_atmos_model_radiation(Atm, Land_ice_atmos_boundary, coupler_clocks)
!$          omp_sec(2) = omp_sec(2) + (omp_get_wtime() - dsec)
!$OMP END PARALLEL
        endif
!$      endif
!$      imb_sec(omp_get_thread_num()+1) = imb_sec(omp_get_thread_num()+1) - omp_get_wtime()
!$OMP END PARALLEL
!$      imb_sec(1) = imb_sec(1) + omp_get_wtime()
!$      if (do_concurrent_radiation) imb_sec(2) = imb_sec(2) + omp_get_wtime()
!$      call omp_set_num_threads(atmos_nthreads+(conc_nthreads-1)*radiation_nthreads)
        !}


        !> UPDATE STATE OF THE ATMOS MODEL 
        call coupler_update_atmos_model_state(Atm, current_timestep, coupler_chksum_obj, coupler_clocks )

        
      enddo fast_integration_loop ! end of na (fast loop)


      !> END CLOCK TO MEASURE ATMOS_LOOP
      call fms_mpp_clock_end(coupler_clocks%atmos_loop)
      

      !> IF DO_LAND, CALL LAND DYNAMICS DRIVER FOR PROCESSES OCCURING AT SLOW TIMESCALE 
      if (do_land) call coupler_update_land_model_slow(Land, Atmos_land_boundary, &
                   Atm%pelist, current_timestep, coupler_chksum_obj, coupler_clocks)

      
      !> TRANSLATE RUNOFF FROM LAND TO ICE GRIDS
      call coupler_flux_land_to_ice(Land, Ice, Land_ice_boundary, Time, current_timestep, &
           coupler_chksum_obj, coupler_clocks)

      
      !> SET ATMOSPHERIC SURFACE PRESSURE TO 0 IN ATMOS_ICE_BOUNDARY TYPE
      ! call flux_atmos_to_ice_slow ?
      Atmos_ice_boundary%p = 0.0

      
      !> UPDATE CURRENT TIME
      Time = Time_atmos

      
      !> END CLOCK FOR MEASURING ATM PROCESSES
      call fms_mpp_clock_end(coupler_clocks%atm)

    endif atm_pe_block

    !> Ice is still using ATM pelist and need to be included in ATM clock
    !> ATM clock is used for load-balancing the coupled models
    start_atm_clock2: if(Atm%pe) then
      call fms_mpp_clock_begin(coupler_clocks%atm)
    end if start_atm_clock2

    !> IF DO_ICE AND ICE%PE
    if (do_ice .and. Ice%pe) then
       
       !> IF ICE_FAST_ICE_PE, CONVERT FIELDS IN LAND_ICE_BOUNDARY TO PRIVATE, FAST_ICE_AVG_TYPE FIELDS IN ICE
       if (Ice%fast_ice_PE) call coupler_unpack_land_ice_boundary(Ice, Land_ice_boundary, coupler_clocks)
       
       
       !> IF NOT CONCURRENT_ICE, COPIES INFORMATION FROM FAST PART OF SEA_ICE TO SLOW PART OF SEA_ICE
       ! This could be a point where the model is serialized; This calls on all ice PEs
       if (.not.concurrent_ice) &
            call coupler_exchange_fast_to_slow_ice(Ice, coupler_clocks, set_ice_current_pelist=.True.)
       
       
       !> IF SLOW_ICE_PE AND NOT COMBINED_ICE_AND_OCEAN, UPDATE SEA-ICE STATE OCCURING AT
       !! SLOWRE TIMESCALE WHICH INCLUDES DYNAMICS, FREEZING AND MELTING, PRECIPITATION,
       !! AND TRANSPORT PROCESSES.  COMPUTE STOCKS AFTERWARDS
       ! This call occurs on whichever PEs handle the slow ice processess.
       if (Ice%slow_ice_PE .and. .not.combined_ice_and_ocean) &
            call coupler_update_ice_model_slow_and_stocks(Ice, coupler_clocks)

       !> IF DO_CHECKSUM, GET CHECKSUM
       if (do_chksum) call coupler_chksum_obj%get_slow_ice_chksums('update_ice_slow+', nc)
    endif  ! End of Ice%pe block

    
    !>  SYNCHRONIZE AND END CLOCK
    end_atm_clock2: if(Atm%pe) then
      call fms_mpp_set_current_pelist(Atm%pelist)
      call fms_mpp_clock_end(coupler_clocks%atm)
    endif end_atm_clock2

    
    !> IF (CONCURRENT_ICE OR NOT USE LAG FLUXES) AND NOT COMBINED_ICE_AND_OCEAN
    if ((concurrent_ice .or. .not.use_lag_fluxes) .and. .not.combined_ice_and_ocean) then
      !this could serialize unless slow_ice_with_ocean is true.
      if ((.not.do_ice) .or. (.not.slow_ice_with_ocean)) call fms_mpp_set_current_pelist()

      !> IF SLOW_ICE_PE OR OCEAN PE, INTERPOLATE ICE MODEL STATE (FLUXES AT THE BOTTOM OF ICE)
      !! AND INTERPOLATE TO THE OCEAN MODEL GRID      
      if (Ice%slow_ice_PE .or. Ocean%is_ocean_pe) &
           call coupler_flux_icell_to_ocean(Ice, Ocean, ice_ocean_boundary, coupler_clocks, &
           slow_ice_ocean_pelist=slow_ice_ocean_pelist, set_current_slow_ice_ocean_pelist=.True.)
      
      !> UPDATE TIME FOR FLUX_ICE_TO_OCEAN
      Time_flux_ice_to_ocean = Time
    endif

    !> IF OCEAN PE
    if (Ocean%is_ocean_pe) then


      !> START CLOCK TO MEASURE OCEAN PROCESSES
      call fms_mpp_set_current_pelist(Ocean%pelist)
      call fms_mpp_clock_begin(coupler_clocks%ocean)

      !> DATA_OVERRIDE AND SEND DATA TO DIAG_MANAGER FOR FIELDS IN ICE_OCEAN_BOUNDARY
      ! This may do data override or diagnostics on Ice_ocean_boundary.
      call flux_ice_to_ocean_finish(Time_flux_ice_to_ocean, Ice_ocean_boundary)


      if (combined_ice_and_ocean) then

         !> IF COMBINED ICE AND OCEAN, GET STOCKS FOR FLUXES FROM ICE TO OCEAN
         call flux_ice_to_ocean_stocks(Ice)

         !> IF COMBINED ICE AND OCEAN,
         !! USES FORCING STORED IN ICE DATATYPE TO ADVANCE SEA-ICE AND ICEBERGS, AND OCEAN STATES
         call update_slow_ice_and_ocean(ice_ocean_driver_CS, Ice, Ocean_state, Ocean, &
              Ice_ocean_boundary, Time_ocean, Time_step_cpld )
      else

         !> IF NOT COMBINED ICE AND OCEAN AND IF DO_CHKSUM, CALL 
         if (do_chksum) call coupler_chksum_obj%get_ocean_chksums('update_ocean_model-', nc)

         !> IF NOT COMBINED ICE AND OCEAN, 
         !! USE FORCINGS IN ICE_OCEAN_BOUNDARY TYPE TO ADVANCE THE OCEAN MODEL'S STATE
         ! update_ocean_model since fluxes don't change here
         if (do_ocean) call coupler_update_ocean_model(Ocean, Ocean_state, Ice_ocean_boundary,&
              Time_ocean, Time_step_cpld, nc, coupler_chksum_obj)
      end if

      !> COMPUTE STOCK 
      ! Get stocks from "Ice_ocean_boundary" and add them to Ocean stocks.
      ! This call is just for record keeping of stocks transfer and
      ! does not modify either Ocean or Ice_ocean_boundary
      call flux_ocean_from_ice_stocks(Ocean_state, Ocean, Ice_ocean_boundary)

      !> CALL FMS DIAG_MANAGER SEND_COMPLETE TO COMPLETE THE SEND_DATA CALLS FOR CURRENT TIMESTEP
      call fms_diag_send_complete(Time_step_cpld)

      !> UPDATE TIME_OCEAN
      Time_ocean = Time_ocean +  Time_step_cpld

      !> UPDATE CURRENT TIME
      Time = Time_ocean

      !> END CLOCK 
      call fms_mpp_clock_end(coupler_clocks%ocean)
    endif

    !> WRITE OUT INTERMEDIATE RESTART FILE WHEN NEEDED.
    if (Time >= Time_restart) &
        call coupler_intermediate_restart(Atm, Ice, Ocean, Ocean_state, Ocn_bc_restart, Ice_bc_restart, &
                                          Time, Time_restart, Time_restart_current, Time_start)

    call coupler_summarize_timestep(nc, num_cpld_calls, coupler_chksum_obj, Atm%pe, omp_sec, imb_sec)

    omp_sec(:)=0.
    imb_sec(:)=0.

  enddo coupled_timestep_loop

  
  if(check_stocks >=0 .and. do_flux) call coupler_flux_init_finish_stocks(Time, Atm, Land, Ice, Ocean_state, &
       coupler_clocks, finish_stocks=.True.)

  call fms_mpp_set_current_pelist()
  call fms_mpp_clock_end(coupler_clocks%main)

  !> CALL COUPLER_END
  call coupler_end(Atm, Land, Ice, Ocean, Ocean_state, Land_ice_atmos_boundary, Atmos_ice_boundary,&
      Atmos_land_boundary, Ice_ocean_boundary, Ocean_ice_boundary, Ocn_bc_restart, Ice_bc_restart, &
      nc, Time, Time_start, Time_end, Time_restart_current, coupler_chksum_obj, coupler_clocks)

  call fms_memutils_print_memuse_stats( 'Memory HiWaterMark', always=.True. )
  call fms_end


end program coupler_main
