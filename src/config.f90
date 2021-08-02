!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Reads in the material, numerics & debug configuration from their respective file
!> @details Reads the material configuration file, where solverJobName.yaml takes
!! precedence over material.yaml.
!--------------------------------------------------------------------------------------------------
module config
  use IO
  use YAML_parse
  use YAML_types
  use results
  use parallelization

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
!> @brief Real *.yaml configuration files.
!--------------------------------------------------------------------------------------------------
subroutine config_init

  print'(/,a)', ' <<<+-  config init  -+>>>'; flush(IO_STDOUT)


  call parse_material
  call parse_numerics
  call parse_debug

end subroutine config_init


!--------------------------------------------------------------------------------------------------
!> @brief Read material.yaml or <jobname>.yaml.
!--------------------------------------------------------------------------------------------------
subroutine parse_material()

  logical :: fileExists
  character(len=:), allocatable :: fileContent


  inquire(file='material.yaml',exist=fileExists)
  if(.not. fileExists) call IO_error(100,ext_msg='material.yaml')

  if (worldrank == 0) then
    print*, 'reading material.yaml'; flush(IO_STDOUT)
    fileContent = IO_read('material.yaml')
    call results_openJobFile(parallel=.false.)
    call results_writeDataset_str(fileContent,'setup','material.yaml','main configuration')
    call results_closeJobFile
  endif
  call parallelization_bcast_str(fileContent)

  config_material => YAML_parse_str(fileContent)

end subroutine parse_material


!--------------------------------------------------------------------------------------------------
!> @brief Read numerics.yaml.
!--------------------------------------------------------------------------------------------------
subroutine parse_numerics()

  logical :: fileExists
  character(len=:), allocatable :: fileContent


  config_numerics => emptyDict

  inquire(file='numerics.yaml', exist=fileExists)
  if (fileExists) then

    if (worldrank == 0) then
      print*, 'reading numerics.yaml'; flush(IO_STDOUT)
      fileContent = IO_read('numerics.yaml')
      call results_openJobFile(parallel=.false.)
      call results_writeDataset_str(fileContent,'setup','numerics.yaml','numerics configuration')
      call results_closeJobFile
    endif
    call parallelization_bcast_str(fileContent)

    config_numerics => YAML_parse_str(fileContent)

  endif

end subroutine parse_numerics


!--------------------------------------------------------------------------------------------------
!> @brief Read debug.yaml.
!--------------------------------------------------------------------------------------------------
subroutine parse_debug()

  logical :: fileExists
  character(len=:), allocatable :: fileContent


  config_debug => emptyDict

  inquire(file='debug.yaml', exist=fileExists)
  if (fileExists) then

    if (worldrank == 0) then
      print*, 'reading debug.yaml'; flush(IO_STDOUT)
      fileContent = IO_read('debug.yaml')
      call results_openJobFile(parallel=.false.)
      call results_writeDataset_str(fileContent,'setup','debug.yaml','debug configuration')
      call results_closeJobFile
    endif
    call parallelization_bcast_str(fileContent)

    config_debug => YAML_parse_str(fileContent)

  endif

end subroutine parse_debug


!--------------------------------------------------------------------------------------------------
!> @brief Deallocate config_material.
!ToDo: deallocation of numerics and debug (optional)
!--------------------------------------------------------------------------------------------------
subroutine config_deallocate

  deallocate(config_material)

end subroutine config_deallocate

end module config
