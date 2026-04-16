! SPDX-License-Identifier: AGPL-3.0-or-later
submodule(homogenization:chemical) chemical_pass

contains

module subroutine pass_init()

  integer :: &
    ho


  print'(/,a)', ' <<<+-  homogenization:chemical:pass init  -+>>>'

  if (count(chemical_type == CHEMICAL_PASS_ID) == 0) return

  print'(/,a,i0)', ' # homogenizations: ',count(chemical_type == CHEMICAL_PASS_ID)

  do ho = 1, size(material_name_homogenization)

    if (chemical_type(ho) /= CHEMICAL_PASS_ID) cycle

    print'(/,1x,a,1x,i0,a)', 'homogenization',ho,': '//material_name_homogenization(ho)

  end do
end subroutine pass_init

end submodule chemical_pass
