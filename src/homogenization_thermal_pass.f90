!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief Dummy homogenization scheme for 1 constituent per material point
!--------------------------------------------------------------------------------------------------
submodule(homogenization:thermal) thermal_pass

contains

module subroutine pass_init()
  
  print'(/,a)', ' <<<+-  homogenization:thermal:pass init  -+>>>'

end subroutine pass_init

end submodule thermal_pass
