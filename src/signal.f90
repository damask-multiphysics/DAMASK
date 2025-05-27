!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Handling of UNIX signals.
!--------------------------------------------------------------------------------------------------
module signal
  use prec

  implicit none(type,external)
  private

  logical, volatile, public, protected :: &
    signal_SIGINT  = .false., &                                                                    !< interrupt signal
    signal_SIGUSR1 = .false., &                                                                    !< 1. user-defined signal
    signal_SIGUSR2 = .false.                                                                       !< 2. user-defined signal

  public :: &
    signal_init, &
    signal_setSIGINT, &
    signal_setSIGUSR1, &
    signal_setSIGUSR2

  interface

    subroutine init_signal_C() bind(C)
    end subroutine init_signal_C

  end interface

contains


!--------------------------------------------------------------------------------------------------
!> @brief Init C signal handler.
!--------------------------------------------------------------------------------------------------
subroutine signal_init()

  call init_signal_C()

end subroutine signal_init


!--------------------------------------------------------------------------------------------------
!> @brief Set global variable signal_SIGINT.
!--------------------------------------------------------------------------------------------------
subroutine signal_setSIGINT(state)

  logical, intent(in) :: state


  signal_SIGINT = state
  print '(/a)', 'SIGINT ' // trim(merge('received', 'cleared ', state))

end subroutine signal_setSIGINT

!--------------------------------------------------------------------------------------------------
!> @brief Set global variable signal_SIGINT to .true.
!> @details To be called from the signal handler in C.
!--------------------------------------------------------------------------------------------------
subroutine signal_setSIGINT_true_Fortran() bind(C,name='signal_setSIGINT_true_Fortran')

  call signal_setSIGINT(.true.)

end subroutine signal_setSIGINT_true_Fortran


!--------------------------------------------------------------------------------------------------
!> @brief Set global variable signal_SIGUSR.
!--------------------------------------------------------------------------------------------------
subroutine signal_setSIGUSR1(state)

  logical, intent(in) :: state


  signal_SIGUSR1 = state
  print '(/a)', 'SIGUSR1 ' // trim(merge('received', 'cleared ', state))

end subroutine signal_setSIGUSR1

!--------------------------------------------------------------------------------------------------
!> @brief Set global variable signal_SIGUSR1 to .true.
!> @details To be called from the signal handler in C.
!--------------------------------------------------------------------------------------------------
subroutine signal_setSIGUSR1_true_Fortran() bind(C,name='signal_setSIGUSR1_true_Fortran')

  call signal_setSIGUSR1(.true.)

end subroutine signal_setSIGUSR1_true_Fortran

!--------------------------------------------------------------------------------------------------
!> @brief Set global variable signal_SIGUSR2.
!--------------------------------------------------------------------------------------------------
subroutine signal_setSIGUSR2(state)

  logical, intent(in) :: state


  signal_SIGUSR2 = state
  print '(/a)', 'SIGUSR2 ' // trim(merge('received', 'cleared ', state))

end subroutine signal_setSIGUSR2

!--------------------------------------------------------------------------------------------------
!> @brief Set global variable signal_SIGUSR2 to .true.
!> @details To be called from the signal handler in C.
!--------------------------------------------------------------------------------------------------
subroutine signal_setSIGUSR2_true_Fortran() bind(C,name='signal_setSIGUSR2_true_Fortran')

  call signal_setSIGUSR2(.true.)

end subroutine signal_setSIGUSR2_true_Fortran


end module signal
