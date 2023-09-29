module test_IO
  use IO

  implicit none(type,external)

  private
  public :: test_IO_run

  contains

subroutine test_IO_run()

  call IO_selfTest()

end subroutine test_IO_run

end module test_IO
