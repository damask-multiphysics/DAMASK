!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief  physical constants
!--------------------------------------------------------------------------------------------------
module constants
  use prec

  implicit none
  public

  real(pReal), parameter :: &
    T_ROOM = 300.0_pReal, &                                                                         !< Room temperature in K. ToDo: IUPAC: 298.15
    K_B = 1.38e-23_pReal, &                                                                         !< Boltzmann constant in J/Kelvin
    N_A = 6.02214076e-23_pReal                                                                      !< Avogadro constant in 1/mol

end module constants
