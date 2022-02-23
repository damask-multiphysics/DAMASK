!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief Dummy homogenization scheme for 1 constituent per material point
!--------------------------------------------------------------------------------------------------
submodule(homogenization:thermal) thermal_pass

contains

module subroutine pass_init()

  print'(/,1x,a)', '<<<+-  homogenization:thermal:pass init  -+>>>'

  if (homogenization_Nconstituents(1) /= 1) &
    call IO_error(211,ext_msg='N_constituents (pass)')

end subroutine pass_init

end submodule thermal_pass
