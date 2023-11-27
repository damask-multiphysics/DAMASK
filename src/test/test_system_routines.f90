module test_system_routines
  use system_routines

  implicit none(type,external)

  private
  public :: test_system_routines_run

  contains

subroutine test_system_routines_run()

  call system_routines_selfTest()

end subroutine test_system_routines_run

end module test_system_routines
