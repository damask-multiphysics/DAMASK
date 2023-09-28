module test_polynomials
  use polynomials

  implicit none(type,external)

  private
  public :: test_polynomials_run

  contains

subroutine test_polynomials_run()

  call polynomials_selfTest()

end subroutine test_polynomials_run

end module test_polynomials
