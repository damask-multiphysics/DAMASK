! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Handling of UNIX signals.
!--------------------------------------------------------------------------------------------------
module signal
  use prec
  use iso_c_binding

  implicit none(type,external)
  private

  logical(C_BOOL), public, protected, volatile, bind(C, name='f_sigint') :: &
    signal_SIGINT  = .false._C_BOOL
  logical(C_BOOL), public, protected, volatile, bind(C, name='f_sigusr1') :: &
    signal_SIGUSR1 = .false._C_BOOL
  logical(C_BOOL), public, protected, volatile, bind(C, name='f_sigusr2') :: &
    signal_SIGUSR2 = .false._C_BOOL
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


  signal_SIGINT = logical(state,C_BOOL)
  print '(/a)', 'SIGINT ' // trim(merge('received', 'cleared ', state))

end subroutine signal_setSIGINT

!--------------------------------------------------------------------------------------------------
!> @brief Set global variable signal_SIGUSR.
!--------------------------------------------------------------------------------------------------
subroutine signal_setSIGUSR1(state)

  logical, intent(in) :: state


  signal_SIGUSR1 = logical(state,C_BOOL)
  print '(/a)', 'SIGUSR1 ' // trim(merge('received', 'cleared ', state))

end subroutine signal_setSIGUSR1

!--------------------------------------------------------------------------------------------------
!> @brief Set global variable signal_SIGUSR2.
!--------------------------------------------------------------------------------------------------
subroutine signal_setSIGUSR2(state)

  logical, intent(in) :: state


  signal_SIGUSR2 = logical(state,C_BOOL)
  print '(/a)', 'SIGUSR2 ' // trim(merge('received', 'cleared ', state))

end subroutine signal_setSIGUSR2

end module signal
