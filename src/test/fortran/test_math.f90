! SPDX-License-Identifier: AGPL-3.0-or-later
module test_math
  use math

  implicit none(type,external)

  private
  public :: test_math_run

  contains

subroutine test_math_run()

  call math_selfTest()

end subroutine test_math_run

end module test_math
