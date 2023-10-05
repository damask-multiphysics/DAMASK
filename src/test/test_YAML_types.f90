module test_YAML_types
  use YAML_types

  implicit none(type,external)

  private
  public :: test_YAML_types_run

  contains

subroutine test_YAML_types_run()

  call YAML_types_selfTest()

end subroutine test_YAML_types_run

end module test_YAML_types
