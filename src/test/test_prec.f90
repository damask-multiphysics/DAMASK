module test_prec
  use prec

  implicit none(type,external)

  private
  public :: test_prec_run

  contains

subroutine test_prec_run()

  call prec_selfTest()

end subroutine test_prec_run

end module test_prec
