!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief  physical constants
!--------------------------------------------------------------------------------------------------
module constants
  use prec

  implicit none(type,external)
  public

  real(pReal), parameter :: &
    T_ROOM = 293.15_pReal, &                                                                        !< Room temperature (20°C) in K (https://en.wikipedia.org/wiki/ISO_1)
    K_B = 1.380649e-23_pReal, &                                                                     !< Boltzmann constant in J/Kelvin (https://doi.org/10.1351/goldbook)
    N_A = 6.02214076e23_pReal                                                                       !< Avogadro constant in 1/mol (https://doi.org/10.1351/goldbook)

end module constants
