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
!> @parblock
!! Module atmos_ocean_dep_fluxes_calc_mod handles computation of ocean and atmosphere deposition gas fluxes
 !! @endparblock
module atmos_ocean_dep_fluxes_calc_mod

  use FMS

  implicit none

  character(len=*), parameter :: mod_name = "aodfc"
  !< mod_name used when printing error messages

contains

  !> @parblock
  !! Subroutine atmos_ocean_dep_fluxes_calc calculates atmosphere-to-ocean wet and dry deposition fluxes.
   !! @endparblock
  subroutine atmos_ocean_dep_fluxes_calc(gas_fields_atm, gas_fields_ice, gas_fluxes, seawater)
    type(FmsCoupler1dBC_type), intent(in) :: gas_fields_atm
      !< is a derived type containing atmospheric surface variables that are used in the calculation
      !! of the atmosphere-ocean gas fluxes.
    type(FmsCoupler1dBC_type), intent(in) :: gas_fields_ice
      !< is a derived type containing ice-top and ocean surface variables that are
      !! used in the calculation of the atmosphere-ocean gas fluxes.
    type(FmsCoupler1dBC_type), intent(inout) :: gas_fluxes
      !< is a derived type containing the gas fluxes between the atmosphere and the ocean and related parameters
    real, dimension(:), intent(in)    :: seawater
      !< is a mask with value of 1 for the open water category, 0 if ice or land.

    character(len=64), parameter    :: sub_name = 'atmos_ocean_dep_fluxes_calc'
    character(len=256), parameter   :: error_header = &
        &'==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

    integer :: n ! Loop index over gas flux boundary conditions
    integer :: i ! Loop index over grid cells
    integer :: length ! Number of grid cells in the current BC flux field
    character(len=128) :: error_string ! Scratch string for formatted error messages

    real, parameter :: permeg=1.0e-6 ! Conversion factor: parts-per-million to fraction (1e-6)

    !> @parblock
    !! RETURN IF NUMBER OF GAS FLUXES AT BOUNDARY IS ZERO.
    !! @endparblock
    if (gas_fluxes%num_bcs .le. 0) return

    !> @parblock
    !! ERROR IF GAS FLUXES BC ARRAY NOT ASSOCIATED.
    !! @endparblock
    if (.not. associated(gas_fluxes%bc)) then
      if (gas_fluxes%num_bcs .ne. 0) then
        call fms_mpp_error(FATAL, trim(error_header) // ' Number of gas fluxes not zero')
      else
        return
      endif
    endif

    !> @parblock
    !! COMPUTE DEPOSITION FLUXES IF FLUX. WAS NOT OVERRIDDEN BY DATA_OVERRIDE AND 
    !! IF FLUX TYPE IS AIR-SEA-DEPOSITION.
    !! @endparblock
    do n = 1, gas_fluxes%num_bcs
      if ( .not. gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%override) then
        if (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_deposition') then
          if (gas_fluxes%bc(n)%param(1) .le. 0.0) then
            write (error_string, '(1pe10.3)') gas_fluxes%bc(n)%param(1)
            call fms_mpp_error(FATAL, 'atmos_ocean_dep_fluxes_calc: Bad parameter (' //&
                & trim(error_string) // ') for air_sea_deposition for ' //&
                & trim(gas_fluxes%bc(n)%name))
          endif

          length = size(gas_fluxes%bc(n)%field(1)%values(:))

          !> @parblock
          !! CALCULATE DEPOSITION FLUXES FOR OPEN WATER CELLS. SET FLUXES TO ZERO FOR ICE AND LAND CELLS.
          !! @endparblock
          if (gas_fluxes%bc(n)%implementation .eq. 'dry') then
            do i = 1, length
              if (seawater(i) == 1.) then
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) = &
                    gas_fields_atm%bc(n)%field(fms_coupler_ind_deposition)%values(i) / gas_fluxes%bc(n)%param(1)
              else
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) = 0.0
              endif
            enddo
          elseif (gas_fluxes%bc(n)%implementation .eq. 'wet') then
            do i = 1, length
              if (seawater(i) == 1.) then
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) = &
                    gas_fields_atm%bc(n)%field(fms_coupler_ind_deposition)%values(i) / gas_fluxes%bc(n)%param(1)
              else
                gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) = 0.0
              endif
            enddo
          else
            call fms_mpp_error(FATAL, 'atmos_ocean_dep_fluxes_calc: Unknown implementation ('&
                & // trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
          endif
        else
          cycle
        endif
      endif
    enddo
  end subroutine  atmos_ocean_dep_fluxes_calc
end module atmos_ocean_dep_fluxes_calc_mod
