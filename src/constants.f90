!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief Constants.
!--------------------------------------------------------------------------------------------------
module constants
  use prec

  implicit none(type,external)
  public

  real(pREAL), parameter :: &
    T_ROOM = 293.15_pREAL, &                                                                        !< Room temperature (20°C) in K (https://en.wikipedia.org/wiki/ISO_1)
    K_B = 1.380649e-23_pREAL, &                                                                     !< Boltzmann constant in J/Kelvin (https://doi.org/10.1351/goldbook)
    N_A = 6.02214076e23_pREAL                                                                       !< Avogadro constant in 1/mol (https://doi.org/10.1351/goldbook)

  character, parameter :: &
    CR = achar(13), &
    LF = new_line('DAMASK')

  character(len=*),          parameter :: LOWER = 'abcdefghijklmnopqrstuvwxyz'
  character(len=len(LOWER)), parameter :: UPPER = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  enum, bind(c); enumerator :: &
    STATUS_OK, &
    STATUS_ITERATING, &
    STATUS_FAIL_PHASE_MECHANICAL, &
    STATUS_FAIL_PHASE_MECHANICAL_STATE, &
    STATUS_FAIL_PHASE_MECHANICAL_DELTASTATE, &
    STATUS_FAIL_PHASE_MECHANICAL_STRESS, &
    STATUS_FAIL_PHASE_MECHANICAL_STRESS_LP_MAXIT, &
    STATUS_FAIL_PHASE_MECHANICAL_STRESS_LP_NAN, &
    STATUS_FAIL_PHASE_MECHANICAL_STRESS_LP_DGESV, &
    STATUS_FAIL_PHASE_DAMAGE, &
    STATUS_FAIL_PHASE_DAMAGE_STATE, &
    STATUS_FAIL_PHASE_DAMAGE_DELTASTATE, &
    STATUS_FAIL_PHASE_THERMAL, &
    STATUS_FAIL_PHASE_THERMAL_DOTSTATE
  end enum

end module constants
