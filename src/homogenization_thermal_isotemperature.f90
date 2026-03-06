! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief Isotemperature homogenization
!--------------------------------------------------------------------------------------------------
submodule(homogenization:thermal) isotemperature

contains

module subroutine isotemperature_init()

  integer :: ho


  print'(/,1x,a)', '<<<+-  homogenization:thermal:isotemperature init  -+>>>'

  if (count(thermal_type == THERMAL_ISOTEMPERATURE_ID) == 0) return

  print'(/,a,i0)', ' # homogenizations: ',count(thermal_type == THERMAL_ISOTEMPERATURE_ID)

  do ho = 1, size(material_name_homogenization)

    if (thermal_type(ho) /= THERMAL_ISOTEMPERATURE_ID) cycle

    print'(/,1x,a,1x,i0,a)', 'homogenization',ho,': '//material_name_homogenization(ho)

  end do

end subroutine isotemperature_init

end submodule isotemperature
