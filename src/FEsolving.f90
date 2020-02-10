!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief global variables for flow control
!--------------------------------------------------------------------------------------------------
module FEsolving
  use prec
   
  implicit none
 
  logical :: &
    terminallyIll = .false.                                                                         !< at least one material point is terminally ill

  integer, dimension(2) :: &
    FEsolving_execElem, &                                                                           !< for ping-pong scheme always whole range, otherwise one specific element
    FEsolving_execIP                                                                                !< for ping-pong scheme always range to max IP, otherwise one specific IP
    
#if defined(Marc4DAMASK)
  logical, dimension(:,:), allocatable :: &
    calcMode                                                                                        !< do calculation or simply collect when using ping pong scheme
#endif

end module FEsolving
