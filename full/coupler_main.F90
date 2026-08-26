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
!! @author Bruce Wyman <Bruce.Wyman@noaa.gov>
!! @author V. Balaji <V.Balaji@noaa.gov>

!> @file
!! @brief Main driver program for the fully coupled (atmosphere, land, sea ice, and ocean) GFDL climate model

!> @ingroup coupler_main
program coupler_main

  use omp_lib
  use FMS
  use full_coupler_mod

  implicit none

  ! model defined types. Targets to pointers in coupler_components_obj

  type (atmos_data_type), target :: Atm
  type (land_data_type), target :: Land
  type (ice_data_type), target :: Ice
  ! allow members of ocean type to be aliased (ap)
  type (ocean_public_type), target  :: Ocean
  type (ocean_state_type),  pointer :: Ocean_state => NULL()

  type(atmos_land_boundary_type), target :: Atmos_land_boundary
  type(atmos_ice_boundary_type), target  :: Atmos_ice_boundary
  type(land_ice_atmos_boundary_type), target  :: Land_ice_atmos_boundary
  type(land_ice_boundary_type), target  :: Land_ice_boundary
  type(ice_ocean_boundary_type), target :: Ice_ocean_boundary
  type(ocean_ice_boundary_type), target :: Ocean_ice_boundary
  type(ice_ocean_driver_type), pointer  :: ice_ocean_driver_CS => NULL()

  type(FmsTime_type) :: Time
  type(FmsTime_type) :: Time_step_atmos
  type(FmsTime_type) :: Time_step_cpld
  type(FmsTime_type) :: Time_atmos
  type(FmsTime_type) :: Time_ocean
  type(FmsTime_type) :: Time_flux_ice_to_ocean
  type(FmsTime_type) :: Time_flux_ocean_to_ice

  integer :: num_atmos_calls
  integer :: na
  integer :: num_cpld_calls
  integer :: nc
  integer :: current_timestep

  type(FmsNetcdfDomainFile_t), dimension(:), pointer :: Ice_bc_restart => NULL()
  type(FmsNetcdfDomainFile_t), dimension(:), pointer :: Ocn_bc_restart => NULL()

  type(FmsTime_type) :: Time_restart
  type(FmsTime_type) :: Time_start
  type(FmsTime_type) :: Time_end
  type(FmsTime_type) :: Time_restart_current

  type(coupler_clock_type) :: coupler_clocks
  type(coupler_components_type), target :: coupler_components_obj
  type(coupler_chksum_type) :: coupler_chksum_obj

  integer, allocatable :: ensemble_pelist(:, :)
  integer, allocatable :: slow_ice_ocean_pelist(:)
  integer :: conc_nthreads = 1
  real :: dsec
  real :: omp_sec(2)=0.0
  real :: imb_sec(2)=0.0


  !> @parblock
  !! INITIALIZE FMS MPP_MOD.  MPP_INIT MUST BE INITIALIZED FIRST
  !! BEFORE CREATING A CLOCK WITH MPP_CLOCK_ID
  !! @endparblock
  call fms_mpp_init()


  !> @parblock
  !! START CLOCK TO MEASURE INITIALIZATION ROUTINE
  !! @endparblock
  coupler_clocks%initialization = fms_mpp_clock_id( 'Initialization' )
  call fms_mpp_clock_begin(coupler_clocks%initialization)


  !> @parblock
  !! INITIALIZE FMS
  !! @endparblock
  call fms_init
  call fmsconstants_init
  call fms_affinity_init


  !> @parblock
  !! INITIALIZE COUPLER PROGRAM VARIABLES
  !! @endparblock
  call coupler_init(Atm, Ocean, Land, Ice, Ocean_state, Atmos_land_boundary, Atmos_ice_boundary, &
    Ocean_ice_boundary, Ice_ocean_boundary, Land_ice_atmos_boundary, Land_ice_boundary,          &
    Ice_ocean_driver_CS, Ice_bc_restart, Ocn_bc_restart, ensemble_pelist, slow_ice_ocean_pelist, &
    conc_nthreads, coupler_clocks, coupler_components_obj, coupler_chksum_obj, &
    Time_step_cpld, Time_step_atmos, Time_atmos, Time_ocean, num_cpld_calls,   &
    num_atmos_calls, Time, Time_start, Time_end, Time_restart, Time_restart_current)


  !> @parblock
  !! IF DO_CHKSUM, COMPUTE CHECKSUM OF ATM, LAND, AND ICE FIELDS
  !! @endparblock
  if (do_chksum) call coupler_chksum_obj%get_coupler_chksums('coupler_init+', 0)


  !> @parblock
  !! SYNCHRONIZE ALL PES
  !! @endparblock
  call fms_mpp_set_current_pelist()


  !> @parblock
  !! END CLOCK TO MEASURE INITIALIZATION ROUTINE
  !! @endparblock
  call fms_mpp_clock_end(coupler_clocks%initialization)


  !> @parblock
  !! START CLOCK TO MEASURE MAIN LOOP
  !! @endparblock
  call fms_mpp_clock_begin(coupler_clocks%main)


  !> @parblock
  !! IF CHECK_STOCKS >= 0 AND DO_FLUX, COMPUTE FLUX STOCKS
  !! @endparblock
  if (check_stocks >= 0 .and. do_flux) &
       call coupler_flux_init_finish_stocks(Time, Atm, Land, Ice, Ocean_state, coupler_clocks, init_stocks=.True.)


  !> @parblock
  !! START OCEAN/SLOW-ICE INTEGRATION LOOP
  !! @endparblock
  coupled_timestep_loop : do nc = 1, num_cpld_calls


    !> @parblock
    !! IF DO_CHKSUM, COMPUTE CHECKSUMS OF ATM, LAND, ICE, AND OCEAN FIELDS
    !! @endparblock
    if (do_chksum) then
      call coupler_chksum_obj%get_coupler_chksums('top_of_coupled_loop+', nc)
      call coupler_chksum_obj%get_atmos_ice_land_ocean_chksums('MAIN_LOOP-', nc)
    end if



    !> @parblock
    !! FOR SLOW_ICE_PES AND OCEAN PES, CALL FLUX_OCEAN_TO_ICE, REDISTRIBUTE FLUXES FROM OCEAN TO ICE
    !! AND STORE IN OCEAN_ICE_BOUNDARY TYPE.  IF USE_LAG_FLUXES IS TRUE, REDISTRIBUTE FLUXES AT THE
    !! BOTTOM OF THE ICE TO THE OCEAN MODEL GRID
    !! @endparblock
    ! Calls to flux_ocean_to_ice and flux_ice_to_ocean are all PE communication
    ! points when running concurrently. The calls are placed next to each other in
    ! concurrent mode to avoid multiple synchronizations within the main loop.
    ! With concurrent_ice, these only occur on the ocean PEs.
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


   !> @parblock
   !! IF DO_CHKSUM, COMPUTE CHECKSUMS OF ATM, LAND, ICE, AND OCEAN FIELDS
   !! @endparblock
    if (do_chksum) then
      call coupler_chksum_obj%get_coupler_chksums('flux_ocn2ice+', nc)
      call coupler_chksum_obj%get_atmos_ice_land_ocean_chksums('flux_ocn2ice+', nc)
   end if


   !> @parblock
   !! IF CHECK_STOCKS > 0 AND DO_FLUX IS TRUE, COMPUTE FLUX STOCKS
   !! @endparblock
   ! needs to sit here rather than at the end of the coupler loop.
   if (check_stocks > 0 .and. do_flux) &
        call coupler_flux_check_stocks(nc, Time, Atm, Land, Ice, Ocean_state, coupler_clocks)


   !> @parblock
   !! IF DO_ICE, AND CURRENT PE IS ICE PE,
   !! (1) COPY INFORMATION FROM OCEAN_ICE_BOUNDARY INTO SLOW PART OF ICE DATATYPE
   !! (2) COPY INFORMATION FROM THE FAST PART OF SEA-ICE TO SLOW PART OF SEA ICE
   !! (3) PREPARE ICE SURFACE STATE FOR ATMOSPHERE FAST PHYSICS
   !! @endparblock
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


   !> @parblock
   !! IF PE IS ATMP%PE
   !! @endparblock
    atm_pe_block : if (Atm%pe) then

      !> @parblock
      !! SYNCHRONIZE ATM PES
      !! @endparblock
      if (.NOT.(do_ice.and.Ice%pe) .OR. (ice_npes.NE.atmos_npes)) call fms_mpp_set_current_pelist(Atm%pelist)

      !> @parblock
      !! CALL CHECKSUM FOR ATMOS ICE LAND FIELDS
      !! @endparblock
      if(do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('set_ice_surface+', nc)

      !> @parblock
      !! START CLOCK FOR PROFILING ATM
      !! @endparblock
      call fms_mpp_clock_begin(coupler_clocks%atm)

      !> @parblock
      !! IF DO_FLUX, GENERATE THE SURFACE EXCHANGE GRID TO EXCHANGE FLUXES BETWEEN LAND AND ICE
      !! @endparblock
      if (do_flux) call coupler_generate_sfc_xgrid(Land, Ice, coupler_clocks)

      !> @parblock
      !! SEND ICE MASK TO DIAG_MANAGER BUFFER
      !! @endparblock
      call send_ice_mask_sic(Time)

      !> @parblock
      !! START CLOCK TO PROFILE FAST LOOP INTEGRATION
      !! @endparblock
      call fms_mpp_clock_begin(coupler_clocks%atmos_loop)


      !> @parblock
      !! START ATMOS/FAST-LAND/FAST-ICE INTEGRATION LOOP
      !! @endparblock
      fast_integration_loop : do na = 1, num_atmos_calls

        !> @parblock
        !! INCREMENT TIME_ATMOS BY TIME_STEP_ATMOS
        !! @endparblock
        Time_atmos = Time_atmos + Time_step_atmos

        !> @parblock
        !! INCREMENT CURRENT_TIMESTEP
        !! @endparblock
        current_timestep = (nc-1)*num_atmos_calls+na

        !> @parblock
        !! IF DO_CHECKSUM, COMPUTE CHECKSUMS FOR ATMOS, ICE, LAND FIELDS
        !! @endparblock
        if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('top_of_atmos_loop-', current_timestep)

        !> @parblock
        !! IF DO_ATMOS, COPY %ATM%TR_BOT TO ATM%FIELDS FOR GASES
        !! @endparblock
        if (do_atmos) call coupler_atmos_tracer_driver_gather_data(Atm, coupler_clocks)

        !> @parblock
        !! IF DO_FLUX, COMPUTE THE FLUXES BETWEEN MODEL COMPONENTS
        !! @endparblock
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


        !> @parblock
        !!  START CLOCK TO MEASURE DO_CONCURRENT_RADIATION
        !! @endparblock
        if (do_concurrent_radiation) call fms_mpp_clock_begin(coupler_clocks%concurrent_atmos)

        !> @parblock
        !! IF DO_ATMOS, CALL ATMOSPHERE DRIVER, CALL FV DYNAMICAL CORE DRIVER
        !! @endparblock
        if (do_atmos) &
             call coupler_update_atmos_model_dynamics(Atm, current_timestep, coupler_chksum_obj, coupler_clocks)

        !> @parblock
        !! IF NOT DO_CONCURRENT_RADIATION, CALL THE RADIATION_DRIVER
        !! @endparblock
        if (.not.do_concurrent_radiation) call coupler_update_atmos_model_radiation(Atm, Land_ice_atmos_boundary, &
             coupler_clocks, current_timestep, coupler_chksum_obj)

        !> @parblock
        !! IF DO_ATMOS, CALL PHYSICS_DRIVER_DOWN TO COMPUTE ATMOSPEHRIC TENDENCIES FOR DYNAMICS,
        !! RADIATION, VERTICAL DIFFUSION OF MOMENTUM, TRACERS, AND HEAT/MOISTURE.
        !! FOR HEAT/MOISTURE, ONLY THE DOWNWARD SWEEP OF THE TRIDONAL ELIMINATION IS PERFORMED
        !! @endparblock
        if (do_atmos) call coupler_update_atmos_model_down(Atm, Land_ice_atmos_boundary, &
             current_timestep, coupler_chksum_obj, coupler_clocks)

        !> @parblock
        !! IF DO_FLUX, CORRECT FOR IMPLICIT TREATMENT OF ATMOSPHERIC DIFFUSIVE FLUXES IN FLUX EXCHANGE
        !! FROM ATM TO LAND AND ICE. CHECKSUMS ARE COMPUTED IF DO_CHKSUM=.TRUE.
        !! @endparblock
        if (do_flux) call coupler_flux_down_from_atmos(Atm, Land, Ice, Land_ice_atmos_boundary, Atmos_land_boundary, &
             Atmos_ice_boundary, Time_atmos, current_timestep, coupler_clocks, coupler_chksum_obj)


        !> @parblock
        !! IF DO_LAND, CALL LAND DYNAMICS DRIVER FOR PROCESSES OCCURRING AT FAST TIMESCALE
        !! @endparblock
        if (do_land .AND. land%pe) call coupler_update_land_model_fast(Land, Atmos_land_boundary, Atm%pelist, &
             current_timestep, coupler_chksum_obj, coupler_clocks)

        !> @parblock
        !! IF DO_ICE, RECORD FLUXES IN ICE TYPE AND CALCULATE ICE TEMPERATURE
        !! @endparblock
        if (do_ice .AND. Ice%fast_ice_pe) call coupler_update_ice_model_fast(Ice, Atmos_ice_boundary, Atm%pelist, &
             current_timestep, coupler_chksum_obj, coupler_clocks)


        !> @parblock
        !! IF DO_FLUX, CORRECT FOR FLUXES TO TAKE INTO ACCOUNT THE NEW SURFACE TEMPERATURES IN LAND AND ICE MODELS
        !! @endparblock
        if (do_flux) call coupler_flux_up_to_atmos(Land, Ice, Land_ice_atmos_boundary, Atmos_land_boundary, &
             Atmos_ice_boundary, Time_atmos, current_timestep, coupler_chksum_obj, coupler_clocks)


        !> @parblock
        !! IF DO_ATMOS, COMPUTE UPWARD VERTICAL DIFFUSION OF HEAT/MOISTURE AND MOISTURE PROCESSES:
        !! COMPUTE UPWARD SWEEP OF THE TRIDONAL ELIMINATION FOR HEAT/MOISTURE AND COMPUTE
        !! THE CONVECTIVE AND LARGE-SCALE TENDENCIES
        !! @endparblock
        if (do_atmos) call coupler_update_atmos_model_up(Atm, Land_ice_atmos_boundary, current_timestep, &
             coupler_chksum_obj, coupler_clocks)


        !> @parblock
        !! IF DO_FLUX, COMPUTES DEPOSITION GAS FLUXES BETWEEN ATMOSPHERE AND OCEAN
        !! @endparblock
        if (do_flux) call coupler_flux_atmos_to_ocean(Atm, Atmos_ice_boundary, Ice, Time_atmos)


        !> @parblock
        !! IF DO_CONCURRENT_RADIATION, END CLOCK TO MEASURE CONCURRENT_ATMOS
        !! @endparblock
        if (do_concurrent_radiation) call fms_mpp_clock_end(coupler_clocks%concurrent_atmos)


        !> @parblock
        !! IF DO_CONCURRENT_RADIATION, CALL THE RADIATION DRIVER
        !! @endparblock
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


        !> @parblock
        !! UPDATE STATE OF THE ATMOS MODEL
        !! @endparblock
        call coupler_update_atmos_model_state(Atm, current_timestep, coupler_chksum_obj, coupler_clocks )


      enddo fast_integration_loop ! end of na (fast loop)


      !> @parblock
      !! END CLOCK TO MEASURE ATMOS_LOOP
      !! @endparblock
      call fms_mpp_clock_end(coupler_clocks%atmos_loop)


      !> @parblock
      !! IF DO_LAND, CALL LAND DYNAMICS DRIVER FOR PROCESSES OCCURRING AT SLOW TIMESCALE
      !! @endparblock
      if (do_land) call coupler_update_land_model_slow(Land, Atmos_land_boundary, &
                   Atm%pelist, current_timestep, coupler_chksum_obj, coupler_clocks)


      !> @parblock
      !! TRANSLATE RUNOFF FROM LAND TO ICE GRIDS
      !! @endparblock
      call coupler_flux_land_to_ice(Land, Ice, Land_ice_boundary, Time, current_timestep, &
           coupler_chksum_obj, coupler_clocks)


      !> @parblock
      !! SET ATMOSPHERIC SURFACE PRESSURE TO 0 IN ATMOS_ICE_BOUNDARY TYPE
      !! @endparblock
      ! call flux_atmos_to_ice_slow ?
      Atmos_ice_boundary%p = 0.0


      !> @parblock
      !! UPDATE CURRENT TIME
      !! @endparblock
      Time = Time_atmos


      !> @parblock
      !! END CLOCK FOR MEASURING ATM PROCESSES
      !! @endparblock
      call fms_mpp_clock_end(coupler_clocks%atm)

    endif atm_pe_block

    !> @parblock
    !! Ice is still using ATM pelist and need to be included in ATM clock
    !! ATM clock is used for load-balancing the coupled models
    !! @endparblock
    start_atm_clock2: if(Atm%pe) then
      call fms_mpp_clock_begin(coupler_clocks%atm)
    end if start_atm_clock2

    !> @parblock
    !! IF DO_ICE AND ICE%PE
    !! @endparblock
    if (do_ice .and. Ice%pe) then

       !> @parblock
       !! IF ICE_FAST_ICE_PE, CONVERT FIELDS IN LAND_ICE_BOUNDARY TO PRIVATE, FAST_ICE_AVG_TYPE FIELDS IN ICE
       !! @endparblock
       if (Ice%fast_ice_PE) call coupler_unpack_land_ice_boundary(Ice, Land_ice_boundary, coupler_clocks)


       !> @parblock
       !! IF NOT CONCURRENT_ICE, COPIES INFORMATION FROM FAST PART OF SEA_ICE TO SLOW PART OF SEA_ICE
       !! @endparblock
       ! This could be a point where the model is serialized; This calls on all ice PEs
       if (.not.concurrent_ice) &
            call coupler_exchange_fast_to_slow_ice(Ice, coupler_clocks, set_ice_current_pelist=.True.)


       !> @parblock
       !! IF SLOW_ICE_PE AND NOT COMBINED_ICE_AND_OCEAN, UPDATE SEA-ICE STATE OCCURRING AT
       !! SLOWRE TIMESCALE WHICH INCLUDES DYNAMICS, FREEZING AND MELTING, PRECIPITATION,
       !! AND TRANSPORT PROCESSES.  COMPUTE STOCKS AFTERWARDS
       !! @endparblock
       ! This call occurs on whichever PEs handle the slow ice processes.
       if (Ice%slow_ice_PE .and. .not.combined_ice_and_ocean) &
            call coupler_update_ice_model_slow_and_stocks(Ice, coupler_clocks)

       !> @parblock
       !! IF DO_CHECKSUM, GET CHECKSUM
       !! @endparblock
       if (do_chksum) call coupler_chksum_obj%get_slow_ice_chksums('update_ice_slow+', nc)
    endif  ! End of Ice%pe block


    !> @parblock
    !!  SYNCHRONIZE AND END CLOCK
    !! @endparblock
    end_atm_clock2: if(Atm%pe) then
      call fms_mpp_set_current_pelist(Atm%pelist)
      call fms_mpp_clock_end(coupler_clocks%atm)
    endif end_atm_clock2


    !> @parblock
    !! IF (CONCURRENT_ICE OR NOT USE LAG FLUXES) AND NOT COMBINED_ICE_AND_OCEAN
    !! @endparblock
    if ((concurrent_ice .or. .not.use_lag_fluxes) .and. .not.combined_ice_and_ocean) then
      !this could serialize unless slow_ice_with_ocean is true.
      if ((.not.do_ice) .or. (.not.slow_ice_with_ocean)) call fms_mpp_set_current_pelist()

      !> @parblock
      !! IF SLOW_ICE_PE OR OCEAN PE, INTERPOLATE ICE MODEL STATE (FLUXES AT THE BOTTOM OF ICE)
      !! AND INTERPOLATE TO THE OCEAN MODEL GRID
      !! @endparblock
      if (Ice%slow_ice_PE .or. Ocean%is_ocean_pe) &
           call coupler_flux_ice_to_ocean(Ice, Ocean, ice_ocean_boundary, coupler_clocks, &
           slow_ice_ocean_pelist=slow_ice_ocean_pelist, set_current_slow_ice_ocean_pelist=.True.)

      !> @parblock
      !! UPDATE TIME FOR FLUX_ICE_TO_OCEAN
      !! @endparblock
      Time_flux_ice_to_ocean = Time
    endif

    !> @parblock
    !! IF OCEAN PE
    !! @endparblock
    if (Ocean%is_ocean_pe) then


      !> @parblock
      !! START CLOCK TO MEASURE OCEAN PROCESSES
      !! @endparblock
      call fms_mpp_set_current_pelist(Ocean%pelist)
      call fms_mpp_clock_begin(coupler_clocks%ocean)

      !> @parblock
      !! DATA_OVERRIDE AND SEND DATA TO DIAG_MANAGER FOR FIELDS IN ICE_OCEAN_BOUNDARY
      !! @endparblock
      ! This may do data override or diagnostics on Ice_ocean_boundary.
      call flux_ice_to_ocean_finish(Time_flux_ice_to_ocean, Ice_ocean_boundary)


      if (combined_ice_and_ocean) then

         !> @parblock
         !! IF COMBINED ICE AND OCEAN, GET STOCKS FOR FLUXES FROM ICE TO OCEAN
         !! @endparblock
         call flux_ice_to_ocean_stocks(Ice)

         !> @parblock
         !! IF COMBINED ICE AND OCEAN,
         !! USES FORCING STORED IN ICE DATATYPE TO ADVANCE SEA-ICE AND ICEBERGS, AND OCEAN STATES
         !! @endparblock
         call update_slow_ice_and_ocean(ice_ocean_driver_CS, Ice, Ocean_state, Ocean, &
              Ice_ocean_boundary, Time_ocean, Time_step_cpld )
      else

         !> @parblock
         !! IF NOT COMBINED ICE AND OCEAN AND IF DO_CHKSUM, CALL
         !! @endparblock
         if (do_chksum) call coupler_chksum_obj%get_ocean_chksums('update_ocean_model-', nc)

         !> @parblock
         !! IF NOT COMBINED ICE AND OCEAN,
         !! USE FORCINGS IN ICE_OCEAN_BOUNDARY TYPE TO ADVANCE THE OCEAN MODEL'S STATE
         !! @endparblock
         ! update_ocean_model since fluxes don't change here
         if (do_ocean) call coupler_update_ocean_model(Ocean, Ocean_state, Ice_ocean_boundary,&
              Time_ocean, Time_step_cpld, nc, coupler_chksum_obj)
      end if

      !> @parblock
      !! COMPUTE STOCK
      !! @endparblock
      ! Get stocks from "Ice_ocean_boundary" and add them to Ocean stocks.
      ! This call is just for record keeping of stocks transfer and
      ! does not modify either Ocean or Ice_ocean_boundary
      call flux_ocean_from_ice_stocks(Ocean_state, Ocean, Ice_ocean_boundary)

      !> @parblock
      !! CALL FMS DIAG_MANAGER SEND_COMPLETE TO COMPLETE THE SEND_DATA CALLS FOR CURRENT TIMESTEP
      !! @endparblock
      call fms_diag_send_complete(Time_step_cpld)

      !> @parblock
      !! UPDATE TIME_OCEAN
      !! @endparblock
      Time_ocean = Time_ocean +  Time_step_cpld

      !> @parblock
      !! UPDATE CURRENT TIME
      !! @endparblock
      Time = Time_ocean

      !> @parblock
      !! END CLOCK
      !! @endparblock
      call fms_mpp_clock_end(coupler_clocks%ocean)
    endif

    !> @parblock
    !! WRITE OUT INTERMEDIATE RESTART FILE WHEN NEEDED.
    !! @endparblock
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

  !> @parblock
  !! CALL COUPLER_END
  !! @endparblock
  call coupler_end(Atm, Land, Ice, Ocean, Ocean_state, Land_ice_atmos_boundary, Atmos_ice_boundary,&
      Atmos_land_boundary, Ice_ocean_boundary, Ocean_ice_boundary, Ocn_bc_restart, Ice_bc_restart, &
      nc, Time, Time_start, Time_end, Time_restart_current, coupler_chksum_obj, coupler_clocks)

  call fms_memutils_print_memuse_stats( 'Memory HiWaterMark', always=.True. )
  call fms_end


end program coupler_main
