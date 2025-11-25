! SPDX-License-Identifier: AGPL-3.0-or-later
module test_OS
  use OS

  implicit none(type,external)

  private
  public :: test_OS_run

  contains

subroutine test_OS_run()

  call OS_selfTest()

end subroutine test_OS_run

end module test_OS
