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
