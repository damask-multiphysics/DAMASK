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

  do ho = 1, size(thermal_type)

    if (thermal_type(ho) /= THERMAL_PASS_ID) cycle

    if (homogenization_Nconstituents(ho) /= 1) &
      call IO_error(211,ext_msg='(pass) with N_constituents !=1')

  end do

end subroutine pass_init

end submodule thermal_pass
