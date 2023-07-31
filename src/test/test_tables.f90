module test_tables
  use tables

  implicit none(type,external)

  private
  public :: test_tables_run

  contains

subroutine test_tables_run()

  call tables_selfTest()

end subroutine test_tables_run

end module test_tables
