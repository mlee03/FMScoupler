# full_coupler_mod

## variables

| Name | Type | Definition |
|------|------|------------|
| restart_interval | integer, dimension(6), public | is a variable in the coupler namelist (coupler_nml) with format (yr, mo, day, hr, min, sec) to set the time interval for writing the intermediate restarts. Intermediate restarts are not written if restart_interval = [0,0,0,0,0,0] Example nml to write restarts every 6 hours: &coupler_nml restart_interval = 0, 0, 0, 6, 0, 0 / |
| current_date | integer, dimension(6) |  |
| calendar | character(len=17) | is a variable in the coupler namelist (coupler_nml) with forma (yr, mo, day, hr, min, sec) to set the model start date. Example — start a run on 1 January 2000: &coupler_nml current_date = 2000, 1, 1, 0, 0, 0 / |
| force_date_from_namelist | logical |  |
| months | integer, public | is a flag in the coupler namelist (coupler_nml) where a .true. value enforces starting date to be current_date from the namelist. |
| days | integer, public | is a namelist variable to set the number of additional days to simulate |
| hours | integer, public | is a namelist variable to set the number of additional hours to simulate |
| minutes | integer, public | is a namelist variable to set the number of additional minutes to simulate |
| seconds | integer, public | is a namelist variable to set the number of additional seconds to simulate |
| dt_atmos | integer, public | is a namelist variable to set the the time step [s] for the atmospheric model dynamics and fast coupling with land and sea ice |
| dt_cpld | integer, public | is a namelist variable to set the time step [s] for coupling between ocean and atmosphere. This mus be an integral multiple of dt_atmos and dt_ocean. This is the "slow" timestep. |
| atmos_npes | integer, public | is a namelist variable to set the number of MPI ranks (processing element) for the atmosphere |
| ocean_npes | integer, public | is a namelist variable to set the number of MPI ranks (processing element) for the ocean |
| ice_npes | integer, public | is a namelist variable to set the number of MPI ranks (processing element) for the ice |
| land_npes | integer, public | is a namelist variable to set the number of MPI ranks (processing element) for the land |
| atmos_nthreads | integer, public | is a namelist variable to set the number of OpenMP threads to use for the atmosphere |
| ocean_nthreads | integer, public | is a namelist variable to set the number of OpenMP threads to use for the ocean |
| radiation_nthreads | integer, public | is a namelist variable to set the number of threads to use for radiation; is set to atmos_nthreads if do_concurrent_radiation is .false. |
| do_atmos | logical, public | is a namelist flag where if .false., skip atmosphere update |
| do_land | logical, public | is a namelist flag where if .false., skip land model update |
| do_ice | logical, public | is a namelist flag where if .false., skip sea-ice model update |
| do_ocean | logical, public | is a namelist flag where if .false., skip ocean model update |
| do_flux | logical, public | is a namelist flag where if .false., skip all flux exchanges between components. |
| concurrent | logical, public | is a namelist flag where if .TRUE., the ocean model updates concurrently with the atmosphere-land-ice on a separate set of PEs. Concurrent should be .TRUE. if concurrent_ice is .TRUE. If .FALSE., the execution is serial: call atmos... followed by call ocean... |
| do_concurrent_radiation | logical, public | is a namelist flag where if .TRUE., radiation is updated concurrently with atmospheric physics with OpenMP threading (radiation_nthreads) |
| use_lag_fluxes | logical, public | is a namelist flag where if .TRUE., the ocean is forced with surface boundary conditions (SBCs) from one coupling timestep ago. If .FALSE., the ocean is forced with most recent SBCs. For an old leapfrog MOM4 coupling with dt_cpld=dt_ocean, lag fluxes can be shown to be stable and current fluxes to be unconditionally unstable. For dt_cpld>dt_ocean there is probably sufficient damping for MOM4. For more modern ocean models (such as MOM5, GOLD or MOM6) that do not use leapfrog timestepping, use_lag_fluxes=.False. should be much more stable. Controls whether ice-to-ocean flux exchange are called at the beginning or at the end of the time loop |
| concurrent_ice | logical, public | is a namelist flag where if .TRUE., the slow sea-ice is forced with the fluxes that were used for the fast ice processes one timestep before. When used in conjuction with setting slow_ice_with_ocean=.TRUE., this approach allows the atmosphere and ocean to run concurrently even if use_lag_fluxes=.FALSE., and it can be shown to ameliorate or eliminate several ice-ocean coupled instabilities. |
| slow_ice_with_ocean | logical, public | is a namelist flag where if true, the slow sea-ice is advanced on the ocean PEs. Otherwise the slow sea-ice processes are on the same PEs as the fast sea-ice. |
| combined_ice_and_ocean | logical, public | is a namelist flag where if true, there is a single call from the coupler to advance both the slow sea-ice and the ocean. slow_ice_with_ocean and concurrent_ice must both be true if combined_ice_and_ocean is true. |
| do_chksum | logical, public | is a namelist flag where if .TRUE., compute checksums throughout the model simulation |
| do_endpoint_chksum | logical, public | is a namelist flag where if .TRUE., do checksums of the initial and final states regardless of do_chksum value. |
| do_debug | logical, public | is a namelist flag where if .TRUE., print additional debugging messages. |
| check_stocks | integer, public | is a namelist flag where value of -1: don't compute stocks; value of 0: compute stocks at the end of the run; value = n > 0: compute stocks every n coupled steps |
| use_hyper_thread | logical, public | is a namelist flag where if .TRUE., enable use of hyperthreading on supported hardware |
| tex | character(len=80) | is a temporary string variable used for printing messages |
| mod_name | character(len=48), parameter | is a string variable for the module name, used in error messages |
| calendar_type | integer | is the calendar_type initialized to INVALID_CALENDAR |
| date_init | integer, dimension(6) | is the initial date for the model run in (yr, mo, day, hr, min, sec) forma |


## coupler_init
### intro
coupler_init is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_init initializes all component models and the flux-exchange infrastructure, and sets all runtime configurations. Coupler_init must be called before the time stepping loops.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| atm | intent(inout) | coupler_init | is the atmosphere derived type |
| land | intent(inout) | coupler_init | is the land derived type |
| ice | intent(inout) | coupler_init | is the sea-ice derived type |
| ocean | intent(inout) | coupler_init | is the ocean public derived type |
| ocean_state | intent(inout) | coupler_init | is the ocean internal state derived type |
| atmos_land_boundary | intent(inout) | coupler_init | is the derived type holding atmos to land fluxes and properies |
| atmos_ice_boundary | intent(inout) | coupler_init | is the derived type holding atmos to ice fluxes and properties |
| ice_ocean_boundary | intent(inout) | coupler_init | is the derived type holding ice to ocean fluxes and properties |
| ocean_ice_boundary | intent(inout) | coupler_init | is the derived type holding ocean to ice fluxes and properties |
| land_ice_boundary | intent(inout) | coupler_init | is the derived type holding land to ice fluxes and properties |
| ice_ocean_driver_cs | intent(inout) | coupler_init | is the control structure for combined ice-ocean driver |
| land_ice_atmos_boundary | intent(inout) | coupler_init | is the derived type holding atmos to/from land and ice fluxes and properties |
| conc_nthreads | intent(inout) | coupler_init | is the number of concurrent OpenMP threads; set to 2 when do_concurrent_radiation=.true. |
| coupler_clocks | intent(inout) | coupler_init | is the collection of FMS clock IDs for profiling |
| coupler_components_obj | intent(inout) | coupler_init | is the object holding pointers to all component model structures |
| coupler_chksum_obj | intent(inout) | coupler_init | is the object used for computing and printing checksums |
| time_step_cpld | intent(inout) | coupler_init | is the coupled (slow) time step |
| time_step_atmos | intent(inout) | coupler_init | is the atmosphere (fast) time step |
| time_atmos | intent(inout) | coupler_init | is the current atmosphere model time |
| time_ocean | intent(inout) | coupler_init | is the current ocean model time |
| time | intent(inout) | coupler_init | is the current model time |
| time_end | intent(inout) | coupler_init | is the model end time |
| num_cpld_calls | intent(inout) | coupler_init | is the total number of coupled (slow) time steps in this run |
| num_atmos_calls | intent(inout) | coupler_init | is the number of atmosphere time steps per coupled time step |

### flowchart
coupler_init does the following:  
Step 1: INITIALIZE STDOUT, STDERR, STDLOG.
Step 2: WRITE WALLDATE AND WALLTIME TO STDERR.
Step 3: WRITE COUPLER VERSION
Step 4: READ COUPLER_NML.
Step 5: IF INPUT/COUPLER.RES EXISTS, READ CALENDAR_TYPE, DATE_INIT, AND DATE IF FILE DOES NOT EXIST, SET FORCE_DATE_FROM_NAMELIST = .TRUE.
Step 6: IF FORCE_DATE_FROM_NAMELIST = .TRUE., SET DATE TO CURRENT_DATE AND CALENDAR_TYPE FROM NAMELIST.
Step 7: INITIALIZE FMS ENSEMBLE_MANAGER FOR AN EMSEMBLE RUN, ENSEMBLE_MANAGER_INIT WILL RENAME ALL THE RESTART AND DIAGNOSTIC FILES TO CONTAIN THE NUMBER OF THE ENSEMBLE MEMBER. TO RESTART AN ENSEMBLE RUN, RESTART FILES MUST EXISTS FOR EACH ENSEMBLE MEMBER.
Step 8: CHECK PE ALLOCATION
Step 9: IF LAND_NPES IS 0, SET LAND_NPES = ATMOS_NPES. LAND_NPES SHOULD BE LESS THAN OR EQUAL TO ATMOS_NPES.
Step 10: IF ICE_NPES IS 0, SET ICE_NPES = ATMOS_NPES. ICE_NPES SHOULD BE LESS THAN OR EQUAL TO ATMOS_NPES
Step 11: SETUP ATMPELIST, OCEANPELIST, LANDPELIST, AND ICEFAST_PELIST IF CONCURRENT = .TRUE. ATMOS AND OCEAN GET DISTINCT PELIST. ELSE ATMOS AND OCEAN HAVE OVERLAPPING PELISTS OR SAME PELIST. LANDPELIST AND ICEFAST_PELIST ARE SUBSETS OF ATMPELIST
Step 12: SET PE IDENTITY FOR ATM, OCEAN, AND LAND, I.E., ATMPE = .TRUE. FOR PE IN ATMPELIST
Step 13: SET ICEPELISTS. OUTLINE BELOW IS FOR WHEN DO_ATMOS IS TRUE: IF SLOW_ICE_WITH_OEAN = .FALSE. THEN ICESLOW_PELIST = ICEFAST_PELIST. IF SLOW_ICE_WITH_OEAN = .TRUE., THEN ICESLOW_PELIST = OCEANPELIST, AND ICEFAST_PELIST = ATMPELIST.
Step 14: SET OMP FOR WHEN DO_CONCURRENT_RADIATION = .TRUE.
Step 15: INITIALIZE CLOCKS FOR PROFILING.
Step 16: WRITE PELISTS TO LOG.
Step 17: WRITE NAMELIST TO LOG.
Step 18: WRITE MODEL INITIAL DATE TO LOGFILE.
Step 19: INITIALIZE DIAG_MANAGER AND READ DIAG_TABLE.
Step 20: OVERRIDE DATE_INIT WITH BASE DATE FROM DIAG_MANAGER IF BASE DATE EXISTS IN DIAG_TABLE.
Step 21: SET TIME_INIT, TIME, AND TIME_START FROM DATE_INIT.
Step 22: SET TIME FROM DATE.
Step 23: SET TIME_START = TIME (from date)COMPUTE TIME_END FROM MONTHS, DAYS, HOURS, MINUTES, AND SECONDS.
Step 24: COMPUTE TIME_END FROM MONTHS, DAYS, HOURS, MINUTES, AND SECONDS.
Step 25: CALL FMS_DIAG_SET_TIME_END WITH TIME_END.
Step 26: GET RUN_LENGTH = TIME_END - TIME.
Step 27: IF INPUT/COUPLER.INTERMEDIATE.RES EXISTS, READ DATE_RESTART FROM THIS FILE. ELSE SET DATE_RESTART = DATE.
Step 28: SET TIME_RESTART.
Step 29: WRITE STARTING AND ENDING TIME TO LOG.
Step 30: SET TIME STEPS AND TOTAL NUMBER OF LOOP ITERATIONS CHECK FOR CONSISTENCY IN TIME STEPS AND RUN LENGTH.
Step 31: INITIALIZE TRACER MANAGER AND GAS EXCHANGE FLUXES.
Step 32: INITIALIZE ATM MODEL ON ATMPES INCLUDING DATA_OVERRIDE_INIT FOR ATM.
Step 33: INITIALIZE LAND MODEL ON LANDPES INCLUDING DATA_OVERRIDE_INIT FOR LAND.
Step 34: INITIALIZE ICE MODEL ON BOTH FAST AND SLOW ICEPES INCLUDING DATA_OVERRIDE_INIT FOR ICE.
Step 35: INITIALIZE OCEAN MODEL ON OCEANPES INCLUDING DATA_OVERRIDE_INIT FOR OCEAN AND OPENMP THREADS
Step 36: CALL MPP_DOMAINS_BROADCAST_DOMAIN FOR ICE AND OCEAN TO SHARE DOMAIN INFORMATION.
Step 37: INITIALIZE FLUX EXCHANGE.
Step 38: SET TIME_ATMOS = TIME_OCEAN = TIME.
Step 39: READ ICE FROM RESTART.
Step 40: READ OCEAN FROM RESTART
Step 41: MISCELLANEOUS INCLUDING CALLING DIAG_GRID_END TO FREE UP MEMORY USED DURING REGIONAL OUTPUT SETUP.
Step 42: INITIALIZE COUPLER_COMPONENT_OBJ
Step 43: INITIALIZE COUPLER CHECKSUM OBJECT.
Step 44: IF DO_ENDPOINT_CHKSUM IS TRUE, COMPUTE CHECKSUM.
Step 45: LOG.


## initialize_coupler_components_obj
### intro
initialize_coupler_components_obj is a subroutine in full_coupler_mod.
### description
  Subroutine initialize_coupler_components_obj is a typed-bound procedure to the coupler_components_type. This subroutine associates each pointer member of the coupler_components_type object to the corresponding model component derived type.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| this | intent(inout) | initialize_coupler_components_obj | is the reference to self (coupler_components_type) |
| atm | intent(inout) | initialize_coupler_components_obj | is the Atm derived type containing atmospheric model data and metadata |
| land | intent(inout) | initialize_coupler_components_obj | is the Land derived type containing land model data and metadata |
| ice | intent(inout) | initialize_coupler_components_obj | is the Ice derived type containing ice model data and metadata |
| ocean | intent(inout) | initialize_coupler_components_obj | is the Ocean derived type containing ocean model data and metadata |
| land_ice_atmos_boundary | intent(inout) | initialize_coupler_components_obj | is the Land_ice_atmos_boundary derived type containing data and metadata for the land-ice-atmosphere boundary |
| atmos_land_boundary | intent(inout) | initialize_coupler_components_obj | is the Atmos_land_boundary derived type containing data and metadata for the atmosphere-land boundary |
| atmos_ice_boundary | intent(inout) | initialize_coupler_components_obj | is the Atmos_ice_boundary derived type containing data and metadata for the atmosphere-ice boundary |
| land_ice_boundary | intent(inout) | initialize_coupler_components_obj | is the Land_ice_boundary derived type containing data and metadata for the land-ice boundary |
| ice_ocean_boundary | intent(inout) | initialize_coupler_components_obj | is the Ice_ocean_boundary derived type containing data and metadata for the ice-ocean boundary |
| ocean_ice_boundary | intent(inout) | initialize_coupler_components_obj | is the Ocean_ice_boundary derived type containing data and metadata for the ocean-ice boundary |

### flowchart
initialize_coupler_components_obj does the following:  
Step 1: POINTER ASSOCIATION


## get_component
### intro
get_component is a subroutine in full_coupler_mod.
### description
  Subroutine get_component is a type-bound procedure to coupler_commponents_type and retrieves the requested component. For example, coupler_components_objget_component(Atm) retrieves coupler_components_objAtm, which is a pointer to the atmospheric component derived type and a private member of coupler_components_obj.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| this | intent(inout) | get_component | is the reference to self (coupler_components_type object) |

### flowchart
get_component does the following:  


## initialize_coupler_chksum_obj
### intro
initialize_coupler_chksum_obj is a subroutine in full_coupler_mod.
### description
  Subroutine initialize_coupler_chksum_obj is a type-bound procedure to coupler_chksum_type and associates coupler_chksum_objcomponents => components_obj. After this call, the chksum objec can access all component model derived types through the pointer.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| this | intent(inout) | initialize_coupler_chksum_obj | The coupler_chksum_type object being initialized |
| components_obj | intent(inout) | initialize_coupler_chksum_obj | The components object whose address will be stored |

### flowchart
initialize_coupler_chksum_obj does the following:  


## get_components_obj
### intro
get_components_obj is a subroutine in full_coupler_mod.
### description
  Subroutine get_components_obj is a type-bound procedure to coupler_chksum_type and retrieves the coupler_components_type object stored inside a coupler_chksum_type object. For example, coupler_chksum_objget_components_obj(components_obj) retrieves the components_obj, which is a pointer to the coupler_components_type object and a private member of coupler_chksum_obj.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| this | intent(inout) | get_components_obj | is a reference to self (coupler_chksum_type) |
| components_obj | intent(inout) | get_components_obj | is the coupler_components_type to be returned |

### flowchart
get_components_obj does the following:  


## coupler_end
### intro
coupler_end is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_end finalizes all component models (such as deallocating arrays), writes restart files, and calls fms_diag_end to flush and close all diagnostic output files. Checksums are computed when do_chksum or do_endpoint_chksum is .true.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| atm | intent(inout) | coupler_end | is the atmospheric derived type |
| land | intent(inout) | coupler_end | is the land derived type |
| ice | intent(inout) | coupler_end | is the ice derived type |
| ocean | intent(inout) | coupler_end | is the ocean derived type |
| ocean_state | intent(inout) | coupler_end | is the ocean state derived type |
| land_ice_atmos_boundary | intent(inout) | coupler_end | is the land-ice-atmosphere boundary derived type |
| atmos_ice_boundary | intent(inout) | coupler_end | is the atmosphere-ice boundary derived type |
| atmos_land_boundary | intent(inout) | coupler_end | is the atmosphere-land boundary derived type |
| ice_ocean_boundary | intent(inout) | coupler_end | is the ice-ocean boundary derived type |
| ocean_ice_boundary | intent(inout) | coupler_end | is the ocean-ice boundary derived type |
| current_timestep | intent(inout) | coupler_end | is the current timestep (nc) |
| coupler_clocks | intent(inout) | coupler_end | are the coupler clocks |
| coupler_chksum_obj | intent(inout) | coupler_end | is required for chksum computations |
| time_end | intent(inout) | coupler_end | is the model ending time |

### flowchart
coupler_end does the following:  
Step 1: IF DO_CHKSUM AND/OR DO_ENDPOINT_CHKSUM IS TRUE, COMPUTE CHECKSUMS
Step 2: CHECK TIME_CURRENT == TIME_END
Step 3: CALL OCEAN_MODEL_END
Step 4: CALL ATMOS_MODEL_END
Step 5: CALL LAND_MODEL_END
Step 6: CALL ICE_MODEL_END
Step 7: WRITE RESTART FILE
Step 8: FINALIZE FMS DIAGNOSTICS MANAGER
Step 9: END CLOCKS


## add_domain_dimension_data
### intro
add_domain_dimension_data is a subroutine in full_coupler_mod.
### description
  Subroutine add_domain_dimension_data writes indices for the x and y dimensions into domain-decomposed fms2_io restart files. This is required so that the FMS tile-combining tool can reconstruct the global field correctly when the I/O layout is not (1,1). Without this call the combiner cannot determine the global position of each tile's data.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| fileobj | intent(inout) | add_domain_dimension_data | is the fms2io domain decomposed fileobj |

### flowchart
add_domain_dimension_data does the following:  


## coupler_restart
### intro
coupler_restart is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_restart writes all coupler-owned restart files.Files written:
RESTART/coupler.res or RESTART/time_stamp.coupler.res (if time_stamp is present): ASCII file containing calendar type integer, model start date (yr,mo,day,hr,min,sec), and current model date.RESTART/coupler.intermediate.res or RESTART/time_stamp.coupler.intermediate.res (if time_stamp is present: ASCII file with the time of the most recent intermediate restart. Written only if Time_restart_current > Time_start.Ocean boundary-condition fields: registered via fms_coupler_type_register_restartsIce boundary-condition fields (Iceocean_fluxes): registered via fms_coupler_type_register_restarts
The optional time_stamp argument, when present, prefixes all file names so tha multiple intermediate restart sets can coexist in RESTART/.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| atm | intent(inout) | coupler_restart | is the atmospheric component derived type |
| ice | intent(inout) | coupler_restart | is the ice component derived type |
| ocean | intent(inout) | coupler_restart | is the ocean component derived type |
| time_stamp | intent(inout) | coupler_restart | is used to determine the restart file as 'RESTART/time_stamp/coupler.res' and 'RESTART/time_stamp/coupler.intermediate.res'. |

### flowchart
coupler_restart does the following:  
Step 1: SET COUPLER.RES AND COUPLER.INTERMEDIATE.RES FILE NAMES
Step 2: WRITE TIME TO COUPLER.RES
Step 3: WRITE DATE TO COUPLER.INTERMEDIATE.RES
Step 4: WRITE OCEANFIELDS RESTARTS
Step 5: WRITE ICEOCEAN_FLUXES RESTARTS


## get_coupler_chksums
### intro
get_coupler_chksums is a subroutine in full_coupler_mod.
### description
  Subroutine get_coupler_chksums computes chksums with fms_mpp_chksums for the following: Atmosphere (Atm) fields:
atmt_bot (temperature at bottom)atmz_bot (height at bottom)atmp_bot (pressure at bottom)atmu_bot (u-wind at bottom)atmv_bot (v-wind at bottom)atmp_surf (surface pressure)atmgust (gustiness)atmtr_bot (atmospheric tracers - dynamically included)
Land fields:
landt_surf (surface temperature)landt_ca (canopy air temperature)landrough_mom (momentum roughness)landrough_heat (heat roughness)landrough_scale (roughness scale)landtr (land tracers - dynamically included)
Ice fields:
icet_surf (surface temperature)icerough_mom (momentum roughness)icerough_heat (heat roughness)icerough_moist (moisture roughness)iceocean_fields (ocean-ice boundary fields).
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| this | intent(inout) | get_coupler_chksums | is a reference to self (coupler_chksum_type object) |
| id | intent(inout) | get_coupler_chksums | id to label CHECKSUMS in stdout, e.g., 'coupler_init+', 'top_of_coupled_loop+', 'coupler_end-', etc |
| timestep | intent(inout) | get_coupler_chksums | timestep to label CHECKSUMS in stdou |

### flowchart
get_coupler_chksums does the following:  


## get_atmos_ice_land_ocean_chksums
### intro
get_atmos_ice_land_ocean_chksums is a subroutine in full_coupler_mod.
### description
  Subroutine get_atmos_ice_land_ocean_chksums calls get_atmos_ice_land_chksums and get_ocean_chksums.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| this | intent(inout) | get_atmos_ice_land_ocean_chksums | is a reference to self (coupler_chksum_type object) |
| id | intent(inout) | get_atmos_ice_land_ocean_chksums | is the id labelling the set of checksums in the logfile |
| timestep | intent(inout) | get_atmos_ice_land_ocean_chksums | is the timestep |

### flowchart
get_atmos_ice_land_ocean_chksums does the following:  


## get_atmos_ice_land_chksums
### intro
get_atmos_ice_land_chksums is a subroutine in full_coupler_mod.
### description
  Subroutine get_atmos_ice_land_chksums computes and prints checksums for atmosphere, fast-ice, and land fields.The pelist must be set (synchronize) before calling this subroutine: if (Atmpe) then call fms_mpp_set_current_pelist(Atmpelist) call coupler_chksum_objget_atmos_ice_land_chksums('MAIN_LOOP-', nc) endif.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| this | intent(inout) | get_atmos_ice_land_chksums | self |
| id | intent(inout) | get_atmos_ice_land_chksums | id to label CHECKSUMS in stdou |

### flowchart
get_atmos_ice_land_chksums does the following:  


## get_slow_ice_chksums
### intro
get_slow_ice_chksums is a subroutine in full_coupler_mod.
### description
  Subroutine get_slow_ice_chksums calls subroutine that will print ou checksums for slow ice and ocean-ice boundary fields. The pelist must be set (synchronize) before calling this subroutine: if (Iceslow_ice_pe) then call mpp_set_current_pelist(Iceslow_pelist) call slow_ice_chksum('MAIN_LOOP-', nc) endif.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| this | intent(inout) | get_slow_ice_chksums | self |
| id | intent(inout) | get_slow_ice_chksums | id to label CHECKSUMS in stdou |

### flowchart
get_slow_ice_chksums does the following:  


## get_ocean_chksums
### intro
get_ocean_chksums is a subroutine in full_coupler_mod.
### description
  Subroutine get_ocean_chksums calls subroutine that will print ou checksums for ocean and ice-ocean boundary fields. The pelist must be set (synchronize) before calling this subroutine: if (Oceanis_ocean_pe) then call mpp_set_current_pelist(Oceanpelist) call ocean_chksum('MAIN_LOOP-', nc) endif.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| this | intent(inout) | get_ocean_chksums | self |
| id | intent(inout) | get_ocean_chksums | ID labelling the set of CHECKSUMS |

### flowchart
get_ocean_chksums does the following:  


## coupler_set_clock_ids
### intro
coupler_set_clock_ids is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_set_clock_ids registers all FMS performance-clock IDs for the coupled model and stores them in the coupler_clocks struct.Clocks are registered on the PE list most appropriate for each phase: atmosphere clocks on Atmpelist, ocean clocks on Oceanpelist, ice clocks on Icefast_pelis or Iceslow_pelist, and ocean-ice flux clocks slow_ice_ocean_pelist. Global clocks (main loop, termination, flux_check_stocks) are registered on all PEs. This routine must be called after PE lists have been set up but before any clock is started.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| coupler_clocks | intent(inout) | coupler_set_clock_ids | is a derived type containing clocks for profiling |
| atm | intent(inout) | coupler_set_clock_ids | is the atm derived type, required for atm pelis |
| land | intent(inout) | coupler_set_clock_ids | is the land derived type, required for land pelis |
| ocean | intent(inout) | coupler_set_clock_ids | is the ocean derived type, required for ocean pelis |
| ice | intent(inout) | coupler_set_clock_ids | is the ice derived type, required for ice pelis |
| ensemble_id | intent(inout) | coupler_set_clock_ids | is the ensemble_id used as index in ensemble_pelis |

### flowchart
coupler_set_clock_ids does the following:  


## coupler_flux_init_finish_stocks
### intro
coupler_flux_init_finish_stocks is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_flux_init_finish_stocks initializes or finalizes stock computation to check for water, heat, and salt conservation.
When init_stocks=.true., calls flux_init_stocks to establish the baseline globally integrated water, heat, and salt stocks (q_start) for all four component models at the start of the run.When finish_stocks=.true., calls flux_check_stocks (if check_stocks >= 0) to compute final stocks, compares them to q_start, and reports conservation errors to the stocks output file.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| time | intent(inout) | coupler_flux_init_finish_stocks | is the current model time |
| atm | intent(inout) | coupler_flux_init_finish_stocks | is the atmosphere derived type |
| land | intent(inout) | coupler_flux_init_finish_stocks | is the land derived type |
| ice | intent(inout) | coupler_flux_init_finish_stocks | is the ice derived type |
| ocean_state | intent(inout) | coupler_flux_init_finish_stocks | is the ocean state derived type |
| coupler_clocks | intent(inout) | coupler_flux_init_finish_stocks | contains the clocks for profiling |
| init_stocks | intent(inout) | coupler_flux_init_finish_stocks | is a flag where if true, call flux_init_stocks |
| finish_stocks | intent(inout) | coupler_flux_init_finish_stocks | is a flag where if true, call final flux_check_stocks |

### flowchart
coupler_flux_init_finish_stocks does the following:  


## coupler_flux_check_stocks
### intro
coupler_flux_check_stocks is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_flux_check_stocks periodically computes water, heat, and salt stocks for all four components.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| nc | intent(inout) | coupler_flux_check_stocks | is the current outer-loop timestep |
| time | intent(inout) | coupler_flux_check_stocks | is the current model time |
| atm | intent(inout) | coupler_flux_check_stocks | is the atmosphere componen |
| land | intent(inout) | coupler_flux_check_stocks | is the land componen |
| ice | intent(inout) | coupler_flux_check_stocks | is the ice componen |
| ocean_state | intent(inout) | coupler_flux_check_stocks | is the ocean state componen |
| coupler_clocks | intent(inout) | coupler_flux_check_stocks | are the coupler clocks |

### flowchart
coupler_flux_check_stocks does the following:  


## coupler_flux_ocean_to_ice
### intro
coupler_flux_ocean_to_ice is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_flux_ocean_to_ice calls flux_ocean_to_ice to transfers the current ocean state (SST, surface currents, salinity, sea-surface height) into the Ocean_ice_boundary derived type in preparation for the slow-ice update. The call occurs on slow_ice_ocean_pelist (the union of slow-ice and ocean PEs) and is profiled with coupler_clocksflux_ocean_to_ice.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| ocean | intent(inout) | coupler_flux_ocean_to_ice | is the ocean componen |
| ice | intent(inout) | coupler_flux_ocean_to_ice | is the ice componen |
| ocean_ice_boundary | intent(inout) | coupler_flux_ocean_to_ice | is the ocean-ice boundary componen |
| coupler_clocks | intent(inout) | coupler_flux_ocean_to_ice | are the coupler clocks |

### flowchart
coupler_flux_ocean_to_ice does the following:  


## coupler_flux_ice_to_ocean
### intro
coupler_flux_ice_to_ocean is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_flux_ice_to_ocean updates the accumulated ice-to-ocean forcing fluxes (heat, freshwater, salt, momentum, shortwave) in Ice_ocean_boundary in preparation for the ocean model update.The optional set_current_slow_ice_ocean_pelist flag controls whether fms_mpp_set_current_pelist(slow_ice_ocean_pelist) is called. It defaults to .false. because when this routine follows coupler_flux_ocean_to_ice, the PE list is already set to slow_ice_ocean_pelist by that routine. Pass .true. when calling coupler_flux_ice_to_ocean independently (e.g., in lag-flux mode).
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| ice | intent(inout) | coupler_flux_ice_to_ocean | is the Ice componen |
| ocean | intent(inout) | coupler_flux_ice_to_ocean | is the Ocean componen |
| ice_ocean_boundary | intent(inout) | coupler_flux_ice_to_ocean | is the Ice_ocean_boundary componen |
| coupler_clocks | intent(inout) | coupler_flux_ice_to_ocean | are the coupler_clocks |

### flowchart
coupler_flux_ice_to_ocean does the following:  


## coupler_unpack_ocean_ice_boundary
### intro
coupler_unpack_ocean_ice_boundary is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_unpack_ocean_ice_boundary, called after coupler_flux_ocean_to_ice, first calls flux_ocean_to_ice_finish to override data (if field exists in data_table) and then calls unpack_ocean_ice_boundary to unpack the ocean-ice boundary data into the ice model state. slow_ice_chksums are computed if do_chksum is true.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| nc | intent(inout) | coupler_unpack_ocean_ice_boundary | is the current outer loop timestep |
| time_flux_ocean_to_ice | intent(inout) | coupler_unpack_ocean_ice_boundary | is the time for flux_ocean_to_ice |
| ice | intent(inout) | coupler_unpack_ocean_ice_boundary | is the Ice componen |
| ocean_ice_boundary | intent(inout) | coupler_unpack_ocean_ice_boundary | is the Ocean_ice_boundary |
| coupler_clocks | intent(inout) | coupler_unpack_ocean_ice_boundary | are the coupler_clocks |
| coupler_chksum_obj | intent(inout) | coupler_unpack_ocean_ice_boundary | is used for computing slow-ice checksums when do_chksum=.true. |

### flowchart
coupler_unpack_ocean_ice_boundary does the following:  


## coupler_exchange_slow_to_fast_ice
### intro
coupler_exchange_slow_to_fast_ice is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_exchange_slow_to_fast_ice transfers updated ocean boundary state to fast ice procesess by calling exchange_slow_to_fast_ice from ice_model_mod. This subroutine is called after coupler_flux_ocean_to_ice and coupler_unpack_ocean_ice_boundary.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| ice | intent(inout) | coupler_exchange_slow_to_fast_ice | is the Ice componen |
| coupler_clocks | intent(inout) | coupler_exchange_slow_to_fast_ice | are the coupler_clocks |

### flowchart
coupler_exchange_slow_to_fast_ice does the following:  


## coupler_exchange_fast_to_slow_ice
### intro
coupler_exchange_fast_to_slow_ice is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_exchange_fast_to_slow_ice calls exchange_fast_to_slow_ice form ice_model_mod to copy fast part of sea-ice to slow part of sea-ice.The optional set_ice_current_pelist flag, when .true., calls fms_mpp_set_current_pelist(Icepelist) to set and synchronize the pes in the pelist.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| ice | intent(inout) | coupler_exchange_fast_to_slow_ice | is the ice componen |
| coupler_clocks | intent(inout) | coupler_exchange_fast_to_slow_ice | are the coupler_clocks |

### flowchart
coupler_exchange_fast_to_slow_ice does the following:  


## coupler_set_ice_surface_fields
### intro
coupler_set_ice_surface_fields is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_set_ice_surface_fields calls set_ice_surface_fields from ice_model_mod to prepare the ice surface state for atmosphere fast physics, as well as pre-calculate ice radiative properties.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| ice | intent(inout) | coupler_set_ice_surface_fields | is the Ice componen |
| coupler_clocks | intent(inout) | coupler_set_ice_surface_fields | are the coupler_clocks |

### flowchart
coupler_set_ice_surface_fields does the following:  


## coupler_generate_sfc_xgrid
### intro
coupler_generate_sfc_xgrid is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_generate_sfc_xgrid calls generate_sfc_xgrid to rebuild the atmosphere-surface exchange grid (xmap_sfc) from the current land mask and ice concentration.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| land | intent(inout) | coupler_generate_sfc_xgrid | is the Land componen |
| ice | intent(inout) | coupler_generate_sfc_xgrid | is the Ice componen |
| coupler_clocks | intent(inout) | coupler_generate_sfc_xgrid | are the coupler_clocks |

### flowchart
coupler_generate_sfc_xgrid does the following:  


## coupler_atmos_tracer_driver_gather_data
### intro
coupler_atmos_tracer_driver_gather_data is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_atmos_tracer_driver_gather_data calls atmos_tracer_driver_gather_data from atmos_tracer_driver_mod to gather CO2, NH3, and tagged/isotopic NH3 at the bottom atm layer in preparation for flux exchange.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| atm | intent(inout) | coupler_atmos_tracer_driver_gather_data | is the atm componen |
| coupler_clocks | intent(inout) | coupler_atmos_tracer_driver_gather_data | are the coupler_clocks |

### flowchart
coupler_atmos_tracer_driver_gather_data does the following:  


## coupler_sfc_boundary_layer
### intro
coupler_sfc_boundary_layer is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_sfc_boundary_layer sets the clock and calls sfc_boundary_layer to compute fluxes at the surface. Chksum is computed if do_chksum is true.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| atm | intent(inout) | coupler_sfc_boundary_layer | is the atm componen |
| land | intent(inout) | coupler_sfc_boundary_layer | is the land componen |
| ice | intent(inout) | coupler_sfc_boundary_layer | is the Ice componen |
| land_ice_atmos_boundary | intent(inout) | coupler_sfc_boundary_layer | is the Land_ice_atmos_boundary componen |
| time_atmos | intent(inout) | coupler_sfc_boundary_layer | is the Atmos time |
| current_timestep | intent(inout) | coupler_sfc_boundary_layer | is the timestep (nc-1)*num_atmos_cal + na |
| coupler_chksum_obj | intent(inout) | coupler_sfc_boundary_layer | is the coupler_chksum_obj |
| coupler_clocks | intent(inout) | coupler_sfc_boundary_layer | are the coupler_clocks |

### flowchart
coupler_sfc_boundary_layer does the following:  


## coupler_update_atmos_model_dynamics
### intro
coupler_update_atmos_model_dynamics is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_update_atmos_model_dynamics calls update_atmos_model_dynamics from atmos_driver to advance the atmospheric dynamical core by one timestep. Checksums are computed when do_chksum=.true., and memory usage is printed when do_debug=.true.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| atm | intent(inout) | coupler_update_atmos_model_dynamics | is the Atm componen |
| current_timestep | intent(inout) | coupler_update_atmos_model_dynamics | is the current timestep |
| coupler_chksum_obj | intent(inout) | coupler_update_atmos_model_dynamics | is the coupler_chksum_obj for computing chksums |
| coupler_clocks | intent(inout) | coupler_update_atmos_model_dynamics | are the coupler_clocks |

### flowchart
coupler_update_atmos_model_dynamics does the following:  


## coupler_update_atmos_model_radiation
### intro
coupler_update_atmos_model_radiation is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_update_atmos_model_radiation calls update_atmos_model_radiation in atmos_driver to update the radiative heating rates, boundary radiative fluxes, and other properties. Checksums are computed if do_chksum is true and do_concurrent_radiation = .false. (due to threading restrictions in mpp_chksum). Memory usage is printed when do_debug=.true.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| atm | intent(inout) | coupler_update_atmos_model_radiation | is the Atm componen |
| land_ice_atmos_boundary | intent(inout) | coupler_update_atmos_model_radiation | is the Land_ice_atmos_boundary componen |
| coupler_clocks | intent(inout) | coupler_update_atmos_model_radiation | are the coupler_clocks |
| current_timestep | intent(inout) | coupler_update_atmos_model_radiation | is the current timestep |
| coupler_chksum_obj | intent(inout) | coupler_update_atmos_model_radiation | is the coupler_chksum_obj for computing chksums |

### flowchart
coupler_update_atmos_model_radiation does the following:  


## coupler_update_atmos_model_down
### intro
coupler_update_atmos_model_down is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_update_atmos_model_down calls update_atmos_model_down from atmos_driver to execute the downward atmospheric physics sweep for heat/moisture. Checksums are computed when do_chksum=.true., and memory usage is printed when do_debug=.true.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| atm | intent(inout) | coupler_update_atmos_model_down | is the atmosphere model derived type |
| land_ice_atmos_boundary | intent(inout) | coupler_update_atmos_model_down | is the derived type containing quantities going from land and ice to atmos |
| current_timestep | intent(inout) | coupler_update_atmos_model_down | is the current coupled timestep index used for checksum labelling |
| coupler_chksum_obj | intent(inout) | coupler_update_atmos_model_down | is the coupler checksum object used to compute checksums |
| coupler_clocks | intent(inout) | coupler_update_atmos_model_down | are the coupler clocks used to measure runtime of the downward physics sweep |

### flowchart
coupler_update_atmos_model_down does the following:  


## coupler_flux_down_from_atmos
### intro
coupler_flux_down_from_atmos is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_flux_down_from_atmos calls flux_down_from_atmos to map fluxes from atmosphere to land and ice components. Runtime is measured by the clock for flux_down_from_atmos, and checksums are computed when do_chksum=.true.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| atm | intent(inout) | coupler_flux_down_from_atmos | is the Atm derived type |
| land | intent(inout) | coupler_flux_down_from_atmos | is the Land derived type |
| ice | intent(inout) | coupler_flux_down_from_atmos | is the Ice derived type |
| land_ice_atmos_boundary | intent(inout) | coupler_flux_down_from_atmos | is the Land_ice_atmos_boundary derived type |
| atmos_land_boundary | intent(inout) | coupler_flux_down_from_atmos | is the Atmos_land_boundary derived type |
| atmos_ice_boundary | intent(inout) | coupler_flux_down_from_atmos | is the Atmos_ice_boundary derived type |
| time_atmos | intent(inout) | coupler_flux_down_from_atmos | is the Time_atmos FmsTime_type containing time in seconds |
| current_timestep | intent(inout) | coupler_flux_down_from_atmos | is the current timestep |
| coupler_clocks | intent(inout) | coupler_flux_down_from_atmos | is the coupler_clocks |
| coupler_chksum_obj | intent(inout) | coupler_flux_down_from_atmos | is the coupler_chksum_obj for computing chksums |

### flowchart
coupler_flux_down_from_atmos does the following:  


## coupler_update_land_model_fast
### intro
coupler_update_land_model_fast is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_update_land_model_fast calls update_land_model_fast from land_model_mod to advance fast land processes by one atmospheric timestep. Clocks are initialized to measure runtime, pelist is set and synchronized before and after fast land model update, and checksums and memory usages are computed if do_chksum and do_debug are true respectively.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| land | intent(inout) | coupler_update_land_model_fast | is the land model derived type |
| atmos_land_boundary | intent(inout) | coupler_update_land_model_fast | is the atmosphere-to-land boundary derived type containing atm to land fluxes |
| current_timestep | intent(inout) | coupler_update_land_model_fast | is the current coupled timestep index used for checksum labelling |
| coupler_chksum_obj | intent(inout) | coupler_update_land_model_fast | is the coupler checksum object used to compute checksums |
| coupler_clocks | intent(inout) | coupler_update_land_model_fast | are the coupler clocks used to measure runtime of the land update |
| coupler_clocks | intent(inout) | coupler_update_land_model_fast | current pelist=Atmpelis |

### flowchart
coupler_update_land_model_fast does the following:  


## coupler_update_ice_model_fast
### intro
coupler_update_ice_model_fast is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_update_ice_model_fast calls update_ice_model_fast from ice_model_mod to advance fast sea ice by one atmospheric timestep. Pelists are set and synchronized before and after the fast ice model update. Runtime is measured by update_ice_model_fast. Checksums and memory usage reporting are controlled by do_chksum and do_debug.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| ice | intent(inout) | coupler_update_ice_model_fast | is the ice model derived type |
| atmos_ice_boundary | intent(inout) | coupler_update_ice_model_fast | is the atmosphere-to-ice boundary derived type containing fluxes passed down from the atmosphere |
| current_timestep | intent(inout) | coupler_update_ice_model_fast | is the current coupled timestep index used for checksum labelling |
| coupler_chksum_obj | intent(inout) | coupler_update_ice_model_fast | is the coupler checksum object used to compute checksums |
| coupler_clocks | intent(inout) | coupler_update_ice_model_fast | are the coupler clocks used to measure runtime of the fast ice update |
| coupler_clocks | intent(inout) | coupler_update_ice_model_fast | current pelist = Atmpelis |

### flowchart
coupler_update_ice_model_fast does the following:  


## coupler_flux_up_to_atmos
### intro
coupler_flux_up_to_atmos is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_flux_up_to_atmos calls flux_up_to_atmos to transfer updated surface states from land and ice to atmosphere. Runtime is measured by coupler_clocksflux_up_to_atmos, and checksums are computed if do_chksum=.true.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| land | intent(inout) | coupler_flux_up_to_atmos | is the land model derived type |
| ice | intent(inout) | coupler_flux_up_to_atmos | is the ice model derived type |
| land_ice_atmos_boundary | intent(inout) | coupler_flux_up_to_atmos | is the land-ice-to-atmosphere boundary derived type accumulating surface fluxes to the atmosphere |
| atmos_land_boundary | intent(inout) | coupler_flux_up_to_atmos | is the atmosphere-to-land boundary derived type used to get dimensions |
| atmos_ice_boundary | intent(inout) | coupler_flux_up_to_atmos | is the atmosphere-to-ice boundary derived type used to get dimensions |
| time_atmos | intent(inout) | coupler_flux_up_to_atmos | is the current atmospheric model time in seconds |
| current_timestep | intent(inout) | coupler_flux_up_to_atmos | is the current coupled timestep index used for checksum labelling |
| coupler_chksum_obj | intent(inout) | coupler_flux_up_to_atmos | is the coupler checksum object used to compute checksums |
| coupler_clocks | intent(inout) | coupler_flux_up_to_atmos | are the coupler clocks used to measure runtime |

### flowchart
coupler_flux_up_to_atmos does the following:  


## coupler_update_atmos_model_up
### intro
coupler_update_atmos_model_up is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_update_atmos_model_up calls update_atmos_model_up from atmos_driver to finish the upward sweep of the tridiagonal eliminiation for heat/moisture and to compute the convective and large-scale tendencies. Runtime is measured by coupler_clocksupdate_atmos_model_up. Checksums are computed when do_chksum=.true., and memory usage is printed when do_debug=.true.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| atm | intent(inout) | coupler_update_atmos_model_up | is the atmosphere model derived type |
| land_ice_atmos_boundary | intent(inout) | coupler_update_atmos_model_up | is the land-ice-to-atmosphere boundary derived type containing surface fluxes returned to the atmosphere |
| current_timestep | intent(inout) | coupler_update_atmos_model_up | is the current coupled timestep index used for checksum labelling |
| coupler_chksum_obj | intent(inout) | coupler_update_atmos_model_up | is the coupler checksum object used to compute checksums |
| coupler_clocks | intent(inout) | coupler_update_atmos_model_up | are the coupler clocks used to measure runtime of the upward atmospheric physics sweep |

### flowchart
coupler_update_atmos_model_up does the following:  


## coupler_flux_atmos_to_ocean
### intro
coupler_flux_atmos_to_ocean is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_flux_atmos_to_ocean calls flux_atmos_to_ocean to compute atmosphere-to-ocean/ice gas deposition fluxes.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| atm | intent(inout) | coupler_flux_atmos_to_ocean | is the atmosphere model derived type |
| atmos_ice_boundary | intent(inout) | coupler_flux_atmos_to_ocean | is the atmosphere-to-ice boundary derived type used to pass gas and deposition fluxes to the ocean |
| ice | intent(inout) | coupler_flux_atmos_to_ocean | is the ice model derived type |
| time_atmos | intent(inout) | coupler_flux_atmos_to_ocean | is the current atmospheric model time in seconds |

### flowchart
coupler_flux_atmos_to_ocean does the following:  


## coupler_update_atmos_model_state
### intro
coupler_update_atmos_model_state is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_update_atmos_model_state calls update_atmos_model_state in atmos_model to update atmospheric state and diagnostic fields in Atm at the end of the atmospheric timestep. Runtime is measured by coupler_clocksupdate_atmos_model_state. Checksums are computed when do_chksum=.true., and memory usage is printed when do_debug=.true.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| atm | intent(inout) | coupler_update_atmos_model_state | is the atmosphere model derived type |
| current_timestep | intent(inout) | coupler_update_atmos_model_state | is the current coupled timestep index used for checksum labelling |
| coupler_chksum_obj | intent(inout) | coupler_update_atmos_model_state | is the coupler checksum object used to compute checksums |
| coupler_clocks | intent(inout) | coupler_update_atmos_model_state | are the coupler clocks used to measure runtime of the atmospheric state update |

### flowchart
coupler_update_atmos_model_state does the following:  


## coupler_update_land_model_slow
### intro
coupler_update_land_model_slow is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_update_land_model_slow calls update_land_model_slow in land_model_mod to advance the land model on the slow (coupled) timestep. Pelist are set and synchronized before and after the call. Runtime is measured with coupler_clocksupdate_land_model_slow. Checksums are computed when do_chksum=.true.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| land | intent(inout) | coupler_update_land_model_slow | is the land model derived type |
| atmos_land_boundary | intent(inout) | coupler_update_land_model_slow | is the atmosphere-to-land boundary derived type containing fluxes passed down from the atmosphere |
| current_timestep | intent(inout) | coupler_update_land_model_slow | is the current coupled timestep index used for checksum labelling |
| coupler_chksum_obj | intent(inout) | coupler_update_land_model_slow | is the coupler checksum object used to compute checksums |
| coupler_clocks | intent(inout) | coupler_update_land_model_slow | are the coupler clocks used to measure runtime of the slow land update |

### flowchart
coupler_update_land_model_slow does the following:  


## coupler_flux_land_to_ice
### intro
coupler_flux_land_to_ice is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_flux_land_to_ice calls flux_land_to_ice to transfer freshwater discharge from the land model to the ice/ocean grid. Runtime is measured with coupler_clocksflux_land_to_ice, and checksums are computed when do_chksum=.true.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| land | intent(inout) | coupler_flux_land_to_ice | is the land model derived type |
| ice | intent(inout) | coupler_flux_land_to_ice | is the ice model derived type |
| land_ice_boundary | intent(inout) | coupler_flux_land_to_ice | is the land-to-ice boundary derived type receiving runoff and other land fluxes |
| time | intent(inout) | coupler_flux_land_to_ice | is the current model time in seconds passed to flux_land_to_ice |
| current_timestep | intent(inout) | coupler_flux_land_to_ice | is the current coupled timestep index used for checksum labelling |
| coupler_chksum_obj | intent(inout) | coupler_flux_land_to_ice | is the coupler checksum object used to compute checksums |
| coupler_clocks | intent(inout) | coupler_flux_land_to_ice | are the coupler clocks used to measure runtime of the land-to-ice flux transfer |

### flowchart
coupler_flux_land_to_ice does the following:  


## coupler_unpack_land_ice_boundary
### intro
coupler_unpack_land_ice_boundary is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_unpack_land_ice_boundary calls ice_model_fast_cleanup and unpack_land_ice_boundary from ice_model_mod to prepare the fast-ice model to receive the new land discharge fields and then copies them into its internal state. Ice_model_fast_cleanup resets the fast-ice accumulation buffers so tha the incoming runoff/calving values replace, rather than accumulate on top of, values from previous steps. Unpack_land_ice_boundary(Ice, Land_ice_boundary copies runoff, calving, runoff_hflx, calving_hflx from Land_ice_boundary that was populated by coupler_flux_land_to_ice into the ice model's internal fast-ice derived types.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| ice | intent(inout) | coupler_unpack_land_ice_boundary | is the ice model derived type |
| land_ice_boundary | intent(inout) | coupler_unpack_land_ice_boundary | is the land-to-ice boundary derived type whose fields are unpacked into the ice model |
| coupler_clocks | intent(inout) | coupler_unpack_land_ice_boundary | are the coupler clocks used to measure runtime of unpacking |

### flowchart
coupler_unpack_land_ice_boundary does the following:  


## coupler_update_ice_model_slow_and_stocks
### intro
coupler_update_ice_model_slow_and_stocks is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_update_ice_model_slow_and_stocks calls update_ice_model_slow from ice_model_mod to advance the slow sea-ice model and then calls flux_ice_to_ocean_stocks to compute stocks. Update_ice_model_slow runs slow-timescale sea-ice processes including dynamics, freezing and melting, precipitation, and transport. Flux_ice_to_ocean_stocks updates for stocks transferred from ice to ocean. Runtime is measured by coupler_clocksupdate_ice_model_slow_slow and coupler_clocksflux_ice_to_ocean_stocks (inner clock).
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| ice | intent(inout) | coupler_update_ice_model_slow_and_stocks | is the ice model derived type |
| coupler_clocks | intent(inout) | coupler_update_ice_model_slow_and_stocks | is the coupler timing clock set used to measure runtime of the slow ice update and stock flux steps |

### flowchart
coupler_update_ice_model_slow_and_stocks does the following:  


## coupler_update_ocean_model
### intro
coupler_update_ocean_model is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_update_ocean_model calls update_ocean_model from ocean_model_mod to advance the ocean model by one coupled timestep with ice-ocean boundary forcing. Time_ocean is advanced by Time_step_cpld, and checksums are computed when do_chksum=.true.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| ocean | intent(inout) | coupler_update_ocean_model | is the ocean model public derived type |
| ocean_state | intent(inout) | coupler_update_ocean_model | is the pointer to the internal ocean model state |
| ice_ocean_boundary | intent(inout) | coupler_update_ocean_model | is the ice-to-ocean boundary derived type containing forcing fluxes passed to the ocean |
| time_ocean | intent(inout) | coupler_update_ocean_model | is the current ocean model time; advanced by Time_step_cpld on outpu |
| time_step_cpld | intent(inout) | coupler_update_ocean_model | is the duration of one coupled (slow) timestep passed to update_ocean_model |
| current_timestep | intent(inout) | coupler_update_ocean_model | is the current coupled timestep index used for checksum labelling |
| coupler_chksum_obj | intent(inout) | coupler_update_ocean_model | is the coupler checksum object used to compute and report field checksums |

### flowchart
coupler_update_ocean_model does the following:  


## coupler_intermediate_restart
### intro
coupler_intermediate_restart is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_intermediate_restart writes mid-run restart files for all component models and the coupler.Component restarts are written on their respective PE sets: atmosphere, land, and ice restarts are written on PEs where Atmpe is true, and the ocean restart is written on PEs where Oceanis_ocean_pe is true. Coupler-specific boundary-condition restart data (Ocn_bc_restart and Ice_bc_restart) are written by coupler_restart in FMS.After all files are written, Time_restart is advanced by restart_interval to set the next scheduled intermediate restart write.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| atm | intent(inout) | coupler_intermediate_restart | is the atmosphere model derived type; restart is written by atmos_model_restar |
| ice | intent(inout) | coupler_intermediate_restart | is the ice model derived type; restart is written by ice_model_restar |
| ocean | intent(inout) | coupler_intermediate_restart | is the ocean model public derived type; restart is written by ocean_model_restar |
| ocean_state | intent(inout) | coupler_intermediate_restart | is the pointer to the internal ocean model state passed to ocean_model_restar |
| ocn_bc_restart | intent(inout) | coupler_intermediate_restart | is the array of fms2_io fileobjs used to write ocean boundary-condition coupler restart data |
| ice_bc_restart | intent(inout) | coupler_intermediate_restart | is the array of fms2_io fileobjs used to write ice boundary-condition coupler restart data |

### flowchart
coupler_intermediate_restart does the following:  


## coupler_summarize_timestep
### intro
coupler_summarize_timestep is a subroutine in full_coupler_mod.
### description
  Subroutine coupler_summarize_timestep reports coupled-timestep progress, memory usage, and optional concurrent-radiation timing diagnostics. Checksums are computed when do_chksum=.true. and summary text is written to stdout each timestep.
### arguments
| Name | Type | Subroutine | Definition |
|------|------|------------|------------|
| current_timestep | intent(inout) | coupler_summarize_timestep | is the current coupled timestep index (nc) used for checksum labelling and progress reporting |
| num_cpld_calls | intent(inout) | coupler_summarize_timestep | is the total number of coupled (outer-loop) timesteps in the run, used for progress reporting |
| coupler_chksum_obj | intent(inout) | coupler_summarize_timestep | is the coupler checksum object used to compute and report end-of-timestep field checksums |
| is_atmos_pe | intent(inout) | coupler_summarize_timestep | is Atmpe; true if this PE belongs to the atmosphere PE list, required for concurrent-radiation timing outpu |
| omp_sec | intent(inout) | coupler_summarize_timestep | is the elapsed wall-clock seconds for each concurrent OpenMP section (atmosphere, radiation) |
| imb_sec | intent(inout) | coupler_summarize_timestep | is the OpenMP load-imbalance seconds for each concurrent OpenMP section |

### flowchart
coupler_summarize_timestep does the following:  

