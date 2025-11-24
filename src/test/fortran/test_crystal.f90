module test_crystal
  use crystal

  implicit none(type,external)

  private
  public :: test_crystal_run

  contains

subroutine test_crystal_run()

  call crystal_selfTest()

end subroutine test_crystal_run

end module test_crystal
