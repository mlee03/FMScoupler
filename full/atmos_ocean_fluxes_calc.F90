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
!> \file
!> \brief Calculates gas fluxes for atmosphere and ocean
module atmos_ocean_fluxes_calc_mod

  use FMS
  use FMSconstants, only: wtmair, rdgas, vonkarm

  implicit none
  private

  public atmos_ocean_fluxes_calc

  character(len=*), parameter :: mod_name = "cdwfe"
  !< is the module name used when printing error messages

  real, parameter :: epsln=1.0e-30
  !< is a really small number used to prevent division by zero

contains
  !> \parblock
  !! atmos_ocean_fluxes_calc calculates atmos-ocean gas fluxes.
  !! All fluxes are in units of [mol/m^2/s] with values > 0 for upward flux.
  !! Deposition fluxes are computed in atmos_ocean_dep_flluxes_calc
  !! \endparblock
  subroutine atmos_ocean_fluxes_calc(gas_fields_atm, gas_fields_ice,&
      & gas_fluxes, seawater, tsurf, ustar, cd_m)
    type(FmsCoupler1dBC_type), intent(in) :: gas_fields_atm
      !< is a derived type containing atmospheric surface variables
    type(FmsCoupler1dBC_type), intent(in) :: gas_fields_ice
      !< is a derived type containing ice-top and ocean surface variables
    type(FmsCoupler1dBC_type), intent(inout) :: gas_fluxes
      !< is a derived type containing the gas fluxes between the atmosphere and the ocean and parameters
    real, dimension(:), intent(in) :: seawater
      !< is a mask with value of 1 for the open water category, 0 if ice or land.
    real, dimension(:), intent(in) :: tsurf !
      !< is the sea-surface temperature [K]; used by the Johnson implementation to compute gas-phase and liquid-phase
      !! transfer velocities.
    real, dimension(:), intent(in), optional :: ustar
      !< is the friction velocity [m/s].  When provided, overrides the internally computed value $u_{10}\sqrt{C_D}
      !! used inside calc_kw.
    real, dimension(:), intent(in), optional :: cd_m
      !< is the drag coefficient (dimensionless).  Only used when ustar is also provided;
     !! otherwise calc_kw uses the bulk formula C_D = 6.1\times10^{-4}+0.63\times10^{-4}u_{10}.

    character(len=*), parameter   :: sub_name = 'atmos_ocean_fluxes_calc'
    character(len=*), parameter   :: error_header =&
        & '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

    integer                                 :: n
    integer                                 :: i
    integer                                 :: length
    real, dimension(:), allocatable         :: kw
    real, dimension(:), allocatable         :: cair
    character(len=128)                      :: error_string

    real, parameter :: permeg=1.0e-6

    !> RETURN IF NUMBER OF GAS FLUXES AT ATMOSPHERE AND OCEAN BOUNDARY IS ZERO
    if (gas_fluxes%num_bcs .le. 0) return

    if (.not. associated(gas_fluxes%bc)) then
      if (gas_fluxes%num_bcs .ne. 0) then
        call fms_mpp_error(FATAL, trim(error_header) // ' Number of gas fluxes not zero')
      else
        return
      endif
    endif

    !> COMPUTE FLUXES AT BOUNDARY FOLLOWING
    !! OCMIP2, DUCE, OR JOHNSON IMPLEMENTATIONS FOR AIR_SEA_GAS_FLUX_GENERIC FLUXES;
    !! OCMIP2, OCMIP2_DATA, OR LINEAR IMPLEMENTATIONS FOR AIR_SEA_GAS_FLUX FLUXES;
    !! RIVER IMPLEMENTATION FOR LAND_SEA_RUNOFF FLUXES.
    !! AIR_SEA_DEPOSITION FLUXES ARE COMPUTED ELSEWHERE IN ATMOS_OCEAN_DEP_FLUXES_MOD
    do n = 1, gas_fluxes%num_bcs
      ! only do calculations if the flux has not been overridden
      if ( .not. gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%override) then
        if (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_gas_flux_generic') then
          length = size(gas_fluxes%bc(n)%field(1)%values(:))

          if (.not. allocated(kw)) then
            allocate( kw(length) )
            allocate ( cair(length) )
          elseif (size(kw(:)) .ne. length) then
            call fms_mpp_error(FATAL, trim(error_header) // ' Lengths of flux fields do not match')
          endif

          if (gas_fluxes%bc(n)%implementation .eq. 'ocmip2') then
            do i = 1, length
              if (seawater(i) == 1.) then
                gas_fluxes%bc(n)%field(fms_coupler_ind_kw)%values(i) =&
                    & gas_fluxes%bc(n)%param(1) * gas_fields_atm%bc(n)%field(fms_coupler_ind_u10)%values(i)**2
                cair(i) = &
                    gas_fields_ice%bc(n)%field(fms_coupler_ind_alpha)%values(i) * &
                    gas_fields_atm%bc(n)%field(fms_coupler_ind_pCair)%values(i) * &
                    gas_fields_atm%bc(n)%field(fms_coupler_ind_psurf)%values(i) * gas_fluxes%bc(n)%param(2)
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) =&
                    & gas_fluxes%bc(n)%field(fms_coupler_ind_kw)%values(i) *&
                    & sqrt(660. / (gas_fields_ice%bc(n)%field(fms_coupler_ind_sc_no)%values(i) + epsln)) *&
                    & (gas_fields_ice%bc(n)%field(fms_coupler_ind_csurf)%values(i) - cair(i))
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux0)%values(i) =&
                    & gas_fluxes%bc(n)%field(fms_coupler_ind_kw)%values(i) *&
                    & sqrt(660. / (gas_fields_ice%bc(n)%field(fms_coupler_ind_sc_no)%values(i) + epsln)) *&
                    & gas_fields_ice%bc(n)%field(fms_coupler_ind_csurf)%values(i)
                gas_fluxes%bc(n)%field(fms_coupler_ind_deltap)%values(i) =&
                    & (gas_fields_ice%bc(n)%field(fms_coupler_ind_csurf)%values(i) - cair(i)) / &
                   (gas_fields_ice%bc(n)%field(fms_coupler_ind_alpha)%values(i) * permeg + epsln)
              else
                gas_fluxes%bc(n)%field(fms_coupler_ind_kw)%values(i) = 0.0
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) = 0.0
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux0)%values(i) = 0.0
                gas_fluxes%bc(n)%field(fms_coupler_ind_deltap)%values(i) = 0.0
                cair(i) = 0.0
              endif
            enddo
          elseif (gas_fluxes%bc(n)%implementation .eq. 'duce') then
            do i = 1, length
              if (seawater(i) == 1.) then
                gas_fluxes%bc(n)%field(fms_coupler_ind_kw)%values(i) = &
                    & gas_fields_atm%bc(n)%field(fms_coupler_ind_u10)%values(i) /&
                    & (770.+45.*gas_fluxes%bc(n)%param(1)**(1./3.)) *&
                    & 101325./(rdgas*wtmair*1e-3*tsurf(i) *&
                    & max(gas_fields_ice%bc(n)%field(fms_coupler_ind_alpha)%values(i),epsln))
                !alpha: mol/m3/atm
                cair(i) = &
                    gas_fields_ice%bc(n)%field(fms_coupler_ind_alpha)%values(i) * &
                    gas_fields_atm%bc(n)%field(fms_coupler_ind_pCair)%values(i) * &
                    gas_fields_atm%bc(n)%field(fms_coupler_ind_psurf)%values(i) * 9.86923e-6
                cair(i) = max(cair(i),0.)
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) =&
                    & gas_fluxes%bc(n)%field(fms_coupler_ind_kw)%values(i) *&
                    & (max(gas_fields_ice%bc(n)%field(fms_coupler_ind_csurf)%values(i),0.) - cair(i))
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux0)%values(i) =&
                    & gas_fluxes%bc(n)%field(fms_coupler_ind_kw)%values(i) *&
                    & max(gas_fields_ice%bc(n)%field(fms_coupler_ind_csurf)%values(i),0.)
                gas_fluxes%bc(n)%field(fms_coupler_ind_deltap)%values(i) =&
                    & (max(gas_fields_ice%bc(n)%field(fms_coupler_ind_csurf)%values(i),0.) - cair(i)) /&
                    & (gas_fields_ice%bc(n)%field(fms_coupler_ind_alpha)%values(i) * permeg + epsln)
              else
                gas_fluxes%bc(n)%field(fms_coupler_ind_kw)%values(i) = 0.0
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) = 0.0
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux0)%values(i) = 0.0
                gas_fluxes%bc(n)%field(fms_coupler_ind_deltap)%values(i) = 0.0
                cair(i) = 0.0
              endif
            enddo
          elseif (gas_fluxes%bc(n)%implementation .eq. 'johnson') then
            !f1p: not sure how to pass salinity. For now, just force at 35.
            do i = 1, length
              if (seawater(i) == 1.) then
                !calc_kw(tk,p,u10,h,vb,mw,sc_w,ustar,cd_m)
                gas_fluxes%bc(n)%field(fms_coupler_ind_kw)%values(i) =&
                    & calc_kw(tsurf(i),&
                    & gas_fields_atm%bc(n)%field(fms_coupler_ind_psurf)%values(i),&
                    & gas_fields_atm%bc(n)%field(fms_coupler_ind_u10)%values(i),&
                    & 101325./(rdgas*wtmair*1e-3*tsurf(i)* &
                               max(gas_fields_ice%bc(n)%field(fms_coupler_ind_alpha)%values(i),epsln)),&
                    & gas_fluxes%bc(n)%param(2),&
                    & gas_fluxes%bc(n)%param(1),&
                    & gas_fields_ice%bc(n)%field(fms_coupler_ind_sc_no)%values(i))
                cair(i) =&
                    & gas_fields_ice%bc(n)%field(fms_coupler_ind_alpha)%values(i) *&
                    & gas_fields_atm%bc(n)%field(fms_coupler_ind_pCair)%values(i) *&
                    & gas_fields_atm%bc(n)%field(fms_coupler_ind_psurf)%values(i) * 9.86923e-6
                cair(i) = max(cair(i),0.)
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) =&
                    & gas_fluxes%bc(n)%field(fms_coupler_ind_kw)%values(i) *&
                    & (max(gas_fields_ice%bc(n)%field(fms_coupler_ind_csurf)%values(i),0.) - cair(i))
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux0)%values(i) =&
                    & gas_fluxes%bc(n)%field(fms_coupler_ind_kw)%values(i) *&
                    & max(gas_fields_ice%bc(n)%field(fms_coupler_ind_csurf)%values(i),0.)
                gas_fluxes%bc(n)%field(fms_coupler_ind_deltap)%values(i) =&
                    & (max(gas_fields_ice%bc(n)%field(fms_coupler_ind_csurf)%values(i),0.) - cair(i)) /&
                    & (gas_fields_ice%bc(n)%field(fms_coupler_ind_alpha)%values(i) * permeg + epsln)
              else
                gas_fluxes%bc(n)%field(fms_coupler_ind_kw)%values(i) = 0.0
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) = 0.0
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux0)%values(i) = 0.0
                gas_fluxes%bc(n)%field(fms_coupler_ind_deltap)%values(i) = 0.0
                cair(i) = 0.0
              endif
            enddo
          else
            call fms_mpp_error(FATAL, ' Unknown implementation (' //&
                & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
          endif
        elseif (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_gas_flux') then
          length = size(gas_fluxes%bc(n)%field(1)%values(:))

          if (.not. allocated(kw)) then
            allocate( kw(length) )
            allocate ( cair(length) )
          elseif (size(kw(:)) .ne. length) then
            call fms_mpp_error(FATAL, trim(error_header) // ' Lengths of flux fields do not match')
          endif

          if (gas_fluxes%bc(n)%implementation .eq. 'ocmip2_data') then
            do i = 1, length
              if (seawater(i) == 1.) then
                kw(i) = gas_fluxes%bc(n)%param(1) * gas_fields_atm%bc(n)%field(fms_coupler_ind_u10)%values(i)
                cair(i) =&
                    & gas_fields_ice%bc(n)%field(fms_coupler_ind_alpha)%values(i) *&
                    & gas_fields_atm%bc(n)%field(fms_coupler_ind_pCair)%values(i) *&
                    & gas_fields_atm%bc(n)%field(fms_coupler_ind_psurf)%values(i) * gas_fluxes%bc(n)%param(2)
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) = kw(i) *&
                    & (gas_fields_ice%bc(n)%field(fms_coupler_ind_csurf)%values(i) - cair(i))
              else
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) = 0.0
                cair(i) = 0.0
                kw(i) = 0.0
              endif
            enddo
          elseif (gas_fluxes%bc(n)%implementation .eq. 'ocmip2') then
            do i = 1, length
              if (seawater(i) == 1.) then
                kw(i) = gas_fluxes%bc(n)%param(1) * gas_fields_atm%bc(n)%field(fms_coupler_ind_u10)%values(i)**2
                cair(i) =&
                    & gas_fields_ice%bc(n)%field(fms_coupler_ind_alpha)%values(i) *&
                    & gas_fields_atm%bc(n)%field(fms_coupler_ind_pCair)%values(i) *&
                    & gas_fields_atm%bc(n)%field(fms_coupler_ind_psurf)%values(i) * gas_fluxes%bc(n)%param(2)
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) = kw(i) *&
                    & (gas_fields_ice%bc(n)%field(fms_coupler_ind_csurf)%values(i) - cair(i))
              else
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) = 0.0
                cair(i) = 0.0
                kw(i) = 0.0
              endif
            enddo
          elseif (gas_fluxes%bc(n)%implementation .eq. 'linear') then
            do i = 1, length
              if (seawater(i) == 1.) then
                kw(i) = gas_fluxes%bc(n)%param(1) *&
                    & max(0.0, gas_fields_atm%bc(n)%field(fms_coupler_ind_u10)%values(i) - gas_fluxes%bc(n)%param(2))
                cair(i) =&
                    & gas_fields_ice%bc(n)%field(fms_coupler_ind_alpha)%values(i) *&
                    & gas_fields_atm%bc(n)%field(fms_coupler_ind_pCair)%values(i) *&
                    & gas_fields_atm%bc(n)%field(fms_coupler_ind_psurf)%values(i) * gas_fluxes%bc(n)%param(3)
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) = kw(i) *&
                    & (gas_fields_ice%bc(n)%field(fms_coupler_ind_csurf)%values(i) - cair(i))
              else
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) = 0.0
                cair(i) = 0.0
                kw(i) = 0.0
              endif
            enddo
          else
            call fms_mpp_error(FATAL, ' Unknown implementation (' //&
                & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
          endif
        elseif (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_deposition') then
          cycle !air_sea_deposition is done in another subroutine
        elseif (gas_fluxes%bc(n)%flux_type .eq. 'land_sea_runoff') then
          if (gas_fluxes%bc(n)%param(1) .le. 0.0) then
            write (error_string, '(1pe10.3)') gas_fluxes%bc(n)%param(1)
            call fms_mpp_error(FATAL, ' Bad parameter (' // trim(error_string) //&
                & ') for land_sea_runoff for ' // trim(gas_fluxes%bc(n)%name))
          endif

          length = size(gas_fluxes%bc(n)%field(1)%values(:))

          if (gas_fluxes%bc(n)%implementation .eq. 'river') then
            do i = 1, length
              if (seawater(i) == 1.) then
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) =&
                    & gas_fields_atm%bc(n)%field(fms_coupler_ind_deposition)%values(i) /&
                    & gas_fluxes%bc(n)%param(1)
              else
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) = 0.0
              endif
            enddo
          else
            call fms_mpp_error(FATAL, ' Unknown implementation (' //&
                & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
          endif
        else
          call fms_mpp_error(FATAL, ' Unknown flux_type (' // trim(gas_fluxes%bc(n)%flux_type) //&
              & ') for ' // trim(gas_fluxes%bc(n)%name))
        endif
      endif
    enddo

    if (allocated(kw)) then
      deallocate(kw)
      deallocate(cair)
    endif
  end subroutine  atmos_ocean_fluxes_calc

  !> \parblock
  !! Calculate total transfer velocities from the "point of view" of liquid
  !! following Johnson Implementation from Johnson, Ocean Science, 2010.
  !! (http://doi.org/10.5194/os-6-913-2010)
  !! Uses equations defined in Liss[1974],
  !!  F = K_g(c_g - H C_l) = K_l(c_g/H - C_l)
  !! where,
  !! F is the flux of gas across air-water interface,
  !! c_g and C_l are the bulk gas and liquid concentrations,
  !! H is the Henry's law constant (H = c_{sg}/C_{sl}),
  !! C_{sg} is the equilibrium concentration in gas phase [g/cm^3 of air] and
  !! C_{sl} is the equilibrium concentration of unionised dissolved gas in liquid phase
  !! [g/cm^3of water]), and
  !! K_g and K_l are the gas-phase and liquid-phase exchange constants, respectively.
  !!    1/K_g = 1/k_g + H/k_l
  !!    1/K_l = 1/k_l + 1/(H*k_g)
  !! \endparblock
  real function calc_kw(tk, p, u10, h, vb, mw, sc_w, ustar, cd_m)
    real, intent(in) :: tk !< is the temperature at surface [K]
    real, intent(in) :: p !< is the pressure at surface [Pa]
    real, intent(in) :: u10 !< is the wind speed at 10m above the surface [m/s]
    real, intent(in) :: h !< is the Henry's law constant (H = c_sg/C_sl) (unitless)
    real, intent(in) :: vb !< is the Molar volume [m^3/mol]
    real, intent(in) :: mw !< is the molecular weight [g/mol]
    real, intent(in) :: sc_w
      !< is the Schmidt number of the gas in seawater (dimensionless); used
      !! to scale the liquid-phase piston velocity k_l relative to the
      !! reference Schmidt number of 660 (CO2 at 20 °C).
    real, intent(in), optional :: ustar
      !< is the Friction velocity [m/s].  If not provided, ustar = u_{10}*sqrt{C_D}.
    real, intent(in), optional :: cd_m
      !< is the Drag coefficient ($C_D).  Used only if ustar is provided.
      !! If ustar is not provided, cd_m = 6.1x10^{-4} + 0.63x10^{-4} * u_10

    real :: ra,rl,tc

    tc = tk-273.15
    ra = 1./max(h*calc_ka(tc,p,mw,vb,u10,ustar,cd_m),epsln)
    rl = 1./max(calc_kl(tc,u10,sc_w),epsln)
    calc_kw = 1./max(ra+rl,epsln)
  end function calc_kw

  !> Calculate total transfer velocities from the "point of view" of gas
  !! following Johnson Implementation from Johnson, Ocean Science, 2010.
  !! (http://doi.org/10.5194/os-6-913-2010)
  !! Uses equations defined in Liss[1974],
  !!  F = K_g(c_g - H C_l) = K_l(c_g/H - C_l)
  !! where,
  !! F is the flux of gas across air-water interface,
  !! c_g and C_l are the bulk gas and liquid concentrations,
  !! H is the Henry's law constant (H = c_{sg}/C_{sl}),
  !! C_{sg} is the equilibrium concentration in gas phase [g/cm^3 of air] and
  !! C_{sl} is the equilibrium concentration of unionised dissolved gas in liquid phase
  !! [g/cm^3of water]), and
  !! K_g and K_l are the gas-phase and liquid-phase exchange constants, respectively.
  !!    1/K_g = 1/k_g + H/k_l
  !!    1/K_l = 1/k_l + 1/(H*k_g)
  !! \endparblock
  real function calc_ka(t, p, mw, vb, u10, ustar, cd_m)
    real, intent(in) :: t !< is the temperature at surface in [C]
    real, intent(in) :: p !< is the pressure at surface in [Pa]
    real, intent(in) :: mw !< is the molecular weight [g/mol]
    real, intent(in) :: vb !< is the molar volume [m^3/mol]
    real, intent(in) :: u10 !< is the wind speed at 10m above the surface in [m/s]
    real, intent(in), optional :: ustar
      !< is the Friction velocity [m/s].  If not provided, ustar = u_{10}*sqrt{C_D}.
    real, intent(in), optional :: cd_m
      !< is the Drag coefficient C_D.  Used only if ustar is provided.
      !! If ustar is not provided, cd_m = 6.1x10^{-4} + 0.63x10^{-4} * u_10

    real             :: sc
    real             :: ustar_t, cd_m_t

    if (.not. present(ustar)) then
      !drag coefficient
      cd_m_t = 6.1e-4 +0.63e-4*u10
      !friction velocity
      ustar_t = u10*sqrt(cd_m_t)
    else
      cd_m_t = cd_m
      ustar_t = ustar
    end if
    sc = schmidt_g(t,p,mw,vb)
    calc_ka = 1e-3+ustar_t/(13.3*sqrt(sc)+1/sqrt(cd_m_t)-5.+log(sc)/(2.*vonkarm))
  end function calc_ka

  !> \parblock
  !! Compute k_l, the liquid-side transfer velocity.  See Johnson, Ocean Science, 2010.
  !! (http://doi.org/10.5194/os-6-913-2010) and Nightingale, Global Biogeochemical Cycles, 2000
  !! (https://doi.org/10.1029/1999GB900091)
  !! \endparblock
  real function calc_kl(t, v, sc)
    real, intent(in) :: t !< is the temperature at surface in C
    real, intent(in) :: v !< is the wind speed at surface in m/s
    real, intent(in) :: sc !< is the Schmidt number of the gas in seawater (dimensionless)

    calc_kl = (((0.222*v**2)+0.333*v)*(max(sc,epsln)/600.)**(-0.5))/(100.*3600.)
  end function calc_kl

  !> \parblock
  !! Compute Schmidt number of the gas in air
  !! \endparblock
  real function schmidt_g(t, p, mw, vb)
    real, intent(in) :: t !< is the temperature at surface in C
    real, intent(in) :: p !< is the pressure at surface in pa
    real, intent(in) :: mw !< is the molecular weight (g/mol)
    real, intent(in) :: vb !< is the molar volume

    real :: d,v

    d = d_air(t,p,mw,vb)
    v = v_air(t)
    schmidt_g = v / d
  end function schmidt_g

  !> \parblock
  !! Compute the diffusion coefficient of the gas in air [m^2/s] following
  !! Fuller, Industrial & Engineering Chemistry (https://doi.org/10.1021/ie50677a007)
  !! \endparblock
  real function d_air(t, p, mw, vb)
    real, intent(in) :: t  !< is the temperature [C]
    real, intent(in) :: p  !< is the pressure [Pa]
    real, intent(in) :: mw !< is the molecular weight [g/mol]
    real, intent(in) :: vb !< is the diffusion coefficient [cm^3/mol]

    real, parameter :: ma = 28.97d0 !< is the molecular weight of air [g/mol]
    real, parameter :: va = 20.1d0  !< is the diffusion volume for air [cm^3/mol]

    real            :: pa

    ! convert p to atm
    pa = 9.8692d-6*p
    d_air = 1d-3 *&
        & (t+273.15d0)**(1.75d0)*sqrt(1d0/ma + 1d0/mw)/(pa*(va**(1d0/3d0)+vb**(1d0/3d0))**2d0)
    ! d_air is in cm2/s convert to m2/s
    d_air = d_air * 1d-4
  end function d_air

  !> \parblock
  !! Compute the density of air [kg/m^3] as a cubic polynomial in temperature.
  !! Coefficients sd_0..sd_3 approximate the dry-air density at standard pressure.
  real function p_air(t)
    real, intent(in) :: t !< is the temperature at surface [C]

    real, parameter :: sd_0 = 1.293393662d0,&
        & sd_1 = -5.538444326d-3,&
        & sd_2 = 3.860201577d-5,&
        & sd_3 = -5.2536065d-7
    p_air = sd_0+(sd_1*t)+(sd_2*t**2)+(sd_3*t**3)
  end function p_air

  !> \parblock
  !! Compute the kinematic viscosity in air [m^2/s]
  !! \endparblock
  real function v_air(t)
    real, intent(in) :: t !< is the temperature at surface [C]
    v_air = n_air(t)/p_air(t)
  end function v_air

  !> \parblock
  !! Compute the dynamic viscosity in air [Pa.s]
  !! \endparblock
  real function n_air(t)
    real, intent(in) :: t !< is the temperature at surface [C]

    real, parameter :: sv_0 = 1.715747771d-5,&
        & sv_1 = 4.722402075d-8,&
        & sv_2 = -3.663027156d-10,&
        & sv_3 = 1.873236686d-12,&
        & sv_4 = -8.050218737d-14
    ! in n.s/m^2 (pa.s)
    n_air = sv_0+(sv_1*t)+(sv_2*t**2)+(sv_3*t**3)+(sv_4*t**4)
  end function n_air
end module atmos_ocean_fluxes_calc_mod
