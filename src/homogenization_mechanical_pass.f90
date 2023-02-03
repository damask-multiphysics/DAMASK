!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief dummy homogenization homogenization scheme for 1 constituent per material point
!--------------------------------------------------------------------------------------------------
submodule(homogenization:mechanical) mechanical_pass

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all necessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
module subroutine pass_init()

  integer :: &
    ho, &
    Nmembers

  print'(/,1x,a)', '<<<+-  homogenization:mechanical:pass init  -+>>>'

  print'(/,a,i0)', ' # homogenizations: ',count(mechanical_type == MECHANICAL_PASS_ID)
  flush(IO_STDOUT)

  do ho = 1, size(mechanical_type)
    if (mechanical_type(ho) /= MECHANICAL_PASS_ID) cycle

    if (homogenization_Nconstituents(ho) /= 1) &
      call IO_error(211,ext_msg='(pass) with N_constituents !=1')

    Nmembers = count(material_ID_homogenization == ho)
    homogState(ho)%sizeState = 0
    allocate(homogState(ho)%state0(0,Nmembers))
    allocate(homogState(ho)%state (0,Nmembers))

  end do

end subroutine pass_init

end submodule mechanical_pass
