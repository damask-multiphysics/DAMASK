!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief global variables for flow control
!--------------------------------------------------------------------------------------------------
module FEsolving

  implicit none
  public

  integer, dimension(2) :: &
    FEsolving_execElem, &                                                                           !< for ping-pong scheme always whole range, otherwise one specific element
    FEsolving_execIP                                                                                !< for ping-pong scheme always range to max IP, otherwise one specific IP

end module FEsolving
