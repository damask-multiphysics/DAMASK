! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief Dummy homogenization scheme for 1 constituent per material point
!--------------------------------------------------------------------------------------------------
submodule(homogenization:thermal) thermal_pass

contains

module subroutine pass_init()

  integer :: &
    ho


  print'(/,1x,a)', '<<<+-  homogenization:thermal:pass init  -+>>>'

  if (count(thermal_type == THERMAL_PASS_ID) == 0) return

  print'(/,a,i0)', ' # homogenizations: ',count(thermal_type == THERMAL_PASS_ID)

  do ho = 1, size(material_name_homogenization)

    if (thermal_type(ho) /= THERMAL_PASS_ID) cycle

    print'(/,1x,a,1x,i0,a)', 'homogenization',ho,': '//material_name_homogenization(ho)

    if (homogenization_Nconstituents(ho) /= 1) &
      call IO_error(211,ext_msg='(pass) with N_constituents !=1')

  end do

end subroutine pass_init

end submodule thermal_pass
