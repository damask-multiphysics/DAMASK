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
    Ninstances, &
    h, &
    Nmaterialpoints

  print'(/,a)', ' <<<+-  homogenization:mechanical:pass init  -+>>>'

  Ninstances = count(homogenization_type == HOMOGENIZATION_NONE_ID)
  print'(a,i2)', ' # instances: ',Ninstances; flush(IO_STDOUT)

  do h = 1, size(homogenization_type)
    if(homogenization_type(h) /= HOMOGENIZATION_NONE_ID) cycle

    if(homogenization_Nconstituents(h) /= 1) &
      call IO_error(211,ext_msg='N_constituents (pass)')

    Nmaterialpoints = count(material_homogenizationAt == h)
    homogState(h)%sizeState = 0
    allocate(homogState(h)%state0   (0,Nmaterialpoints))
    allocate(homogState(h)%state    (0,Nmaterialpoints))

  enddo

end subroutine pass_init

end submodule mechanical_pass
