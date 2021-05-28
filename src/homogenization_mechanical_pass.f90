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
module subroutine pass_init

  integer :: &
    ho, &
    Nmembers

  print'(/,a)', ' <<<+-  homogenization:mechanical:pass init  -+>>>'

  print'(a,i0)', ' # homogenizations: ',count(homogenization_type == HOMOGENIZATION_NONE_ID)
  flush(IO_STDOUT)

  do ho = 1, size(homogenization_type)
    if(homogenization_type(ho) /= HOMOGENIZATION_NONE_ID) cycle

    if(homogenization_Nconstituents(ho) /= 1) &
      call IO_error(211,ext_msg='N_constituents (pass)')

    Nmembers = count(material_homogenizationID == ho)
    homogState(ho)%sizeState = 0
    allocate(homogState(ho)%state0(0,Nmembers))
    allocate(homogState(ho)%state (0,Nmembers))

  enddo

end subroutine pass_init

end submodule mechanical_pass
