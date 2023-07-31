module test_misc
  use misc

  implicit none(type,external)

  private
  public :: test_misc_run

  contains

subroutine test_misc_run()

  call misc_selfTest()

end subroutine test_misc_run

end module test_misc
