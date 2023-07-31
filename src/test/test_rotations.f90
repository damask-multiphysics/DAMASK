module test_rotations
  use rotations

  implicit none(type,external)

  private
  public :: test_rotations_run

  contains

subroutine test_rotations_run()

  call rotations_selfTest()

end subroutine test_rotations_run

end module test_rotations
