module test_types
  use types

  implicit none(type,external)

  private
  public :: test_types_run

  contains

subroutine test_types_run()

  call types_selfTest()

end subroutine test_types_run

end module test_types
