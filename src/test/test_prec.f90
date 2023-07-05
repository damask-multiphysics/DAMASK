module test_prec
  use prec

  implicit none(type,external)

  private
  public :: prec_test

  contains

subroutine prec_test()

  print*, 'begin test prec'
  call prec_selfTest()
  print*, 'end test prec'

end subroutine prec_test

end module test_prec
