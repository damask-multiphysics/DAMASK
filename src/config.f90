!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Read in the material and numerics configuration from their respective file.
!--------------------------------------------------------------------------------------------------
module config
  use IO
  use CLI
  use misc
  use YAML_parse
  use YAML_types
  use result
  use parallelization

  implicit none(type,external)
  private

  type(tDict), pointer, public :: &
    config_material, &
    config_numerics

  public :: &
    config_init, &
    config_material_deallocate, &
    config_numerics_deallocate, &
    config_listReferences

contains

!--------------------------------------------------------------------------------------------------
!> @brief Read *.yaml configuration files.
!--------------------------------------------------------------------------------------------------
subroutine config_init()

  print'(/,1x,a)', '<<<+-  config init  -+>>>'; flush(IO_STDOUT)

  call parse_material()
  call parse_numerics()

end subroutine config_init


!--------------------------------------------------------------------------------------------------
!> @brief Deallocate config_material.
!--------------------------------------------------------------------------------------------------
subroutine config_material_deallocate()

  print'(/,1x,a)', 'deallocating material configuration'; flush(IO_STDOUT)
  deallocate(config_material)

end subroutine config_material_deallocate

!--------------------------------------------------------------------------------------------------
!> @brief Deallocate config_numerics if present.
!--------------------------------------------------------------------------------------------------
subroutine config_numerics_deallocate()

  if (.not. associated(config_numerics, emptyDict)) then
    print'(/,1x,a)', 'deallocating numerics configuration'; flush(IO_STDOUT)
    deallocate(config_numerics)
  end if

end subroutine config_numerics_deallocate


!--------------------------------------------------------------------------------------------------
!> @brief Return string with references from dict.
!--------------------------------------------------------------------------------------------------
function config_listReferences(config,indent) result(references)

  type(tDict) :: config
  integer, optional :: indent
  character(len=:), allocatable :: references


  type(tList), pointer :: ref
  character(len=:), allocatable :: filler
  integer :: r


  filler = repeat(' ',misc_optional(indent,0))
  ref => config%get_list('references',emptyList)
  if (ref%length == 0) then
    references = ''
  else
    references = 'references:'
    do r = 1, ref%length
      references = references//IO_EOL//filler//'- '//IO_wrapLines(ref%get_asString(r),filler=filler//'  ')
    end do
  end if

end function config_listReferences


!--------------------------------------------------------------------------------------------------
!> @brief Read material.yaml.
!--------------------------------------------------------------------------------------------------
subroutine parse_material()

  logical :: fileExists
  character(len=:), allocatable :: &
    fileContent, fname

  if (worldrank == 0) then
    print'(/,1x,a)', 'reading material.yaml'; flush(IO_STDOUT)
    fileContent = IO_read(CLI_materialFile)
    fname = CLI_materialFile
    if (scan(fname,'/') /= 0) fname = fname(scan(fname,'/',.true.)+1:)
    call result_openJobFile(parallel=.false.)
    call result_writeDataset_str(fileContent,'setup',fname,'main configuration')
    call result_closeJobFile()
  end if
  call parallelization_bcast_str(fileContent)

  config_material => YAML_parse_str_asDict(fileContent)

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
      print'(1x,a)', 'reading numerics.yaml'; flush(IO_STDOUT)
      fileContent = IO_read('numerics.yaml')
      if (len(fileContent) > 0) then
        call result_openJobFile(parallel=.false.)
        call result_writeDataset_str(fileContent,'setup','numerics.yaml','numerics configuration')
        call result_closeJobFile()
      end if
    end if
    call parallelization_bcast_str(fileContent)

    config_numerics => YAML_parse_str_asDict(fileContent)

  end if

end subroutine parse_numerics

end module config
