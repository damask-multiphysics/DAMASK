!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief Dummy homogenization scheme for 1 constituent per material point
!--------------------------------------------------------------------------------------------------
submodule(homogenization:damage) damage_pass

contains

module subroutine pass_init()
  
  print'(/,a)', ' <<<+-  homogenization:damage:pass init  -+>>>'

end subroutine pass_init

end submodule damage_pass
