!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Reads in the material, numerics & debug configuration from their respective file
!> @details Reads the material configuration file, where solverJobName.yaml takes
!! precedence over material.yaml.
!--------------------------------------------------------------------------------------------------
module config
  use prec
  use DAMASK_interface
  use IO
  use YAML_parse
  use YAML_types

#ifdef PETSc
#include <petsc/finclude/petscsys.h>
   use petscsys
#endif

  implicit none
  private

  class(tNode), pointer, public :: &
    config_material, &
    config_numerics, &
    config_debug

  public :: &
    config_init, &
    config_deallocate

contains

!--------------------------------------------------------------------------------------------------
!> @brief calls subroutines that reads material, numerics and debug configuration files
!--------------------------------------------------------------------------------------------------
subroutine config_init

  write(6,'(/,a)') ' <<<+-  config init  -+>>>'; flush(6)

  call parse_material
  call parse_numerics
  call parse_debug

end subroutine config_init


!--------------------------------------------------------------------------------------------------
!> @brief reads material.yaml
!--------------------------------------------------------------------------------------------------
subroutine parse_material

  logical :: fileExists
  character(len=:), allocatable :: fname

  fname = getSolverJobName()//'.yaml'
  inquire(file=fname,exist=fileExists)
  if(.not. fileExists) then
    fname = 'material.yaml'
    inquire(file=fname,exist=fileExists)
    if(.not. fileExists) call IO_error(100,ext_msg=fname)
  endif
  print*, 'reading '//fname; flush(6)
  config_material => YAML_parse_file(fname)

end subroutine parse_material


!--------------------------------------------------------------------------------------------------
!> @brief reads in parameters from numerics.yaml and sets openMP related parameters. Also does
! a sanity check
!--------------------------------------------------------------------------------------------------
subroutine parse_numerics

  logical :: fexist

  config_numerics => emptyDict
  inquire(file='numerics.yaml', exist=fexist)
  if (fexist) then
    print*, 'reading numerics.yaml'; flush(6)
    config_numerics => YAML_parse_file('numerics.yaml')
  endif

end subroutine parse_numerics


!--------------------------------------------------------------------------------------------------
!> @brief reads in parameters from debug.yaml
!--------------------------------------------------------------------------------------------------
subroutine parse_debug

  logical :: fexist

  config_debug => emptyDict
  inquire(file='debug.yaml', exist=fexist)
  fileExists: if (fexist) then
    print*, 'reading debug.yaml'; flush(6)
    config_debug => YAML_parse_file('debug.yaml')
  endif fileExists

end subroutine parse_debug


!--------------------------------------------------------------------------------------------------
!> @brief deallocates material.yaml structure
!ToDo: deallocation of numerics debug (optional)
!--------------------------------------------------------------------------------------------------
subroutine config_deallocate

  deallocate(config_material)

end subroutine config_deallocate

end module config
