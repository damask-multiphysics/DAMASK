!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief  physical constants
!--------------------------------------------------------------------------------------------------
module constants
  use prec

  implicit none
  public

  real(pReal), parameter :: &
    T_ROOM = 300.0_pReal, &                                                                         !< Room temperature in K
    K_B = 1.38e-23_pReal                                                                             !< Boltzmann constant in J/Kelvin

end module constants
