!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Read in the material and numerics configuration from their respective file.
!--------------------------------------------------------------------------------------------------
module config
  use IO
  use misc
  use YAML
  use types
  use result
  use parallelization
#if   defined(MESH) || defined(GRID)
  use CLI
#endif
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

#if defined(MESH) || defined(GRID)
  config_material => parse(CLI_materialFile,'material configuration')
#else
  config_material => parse('material.yaml','material configuration')
#endif

  config_numerics => emptyDict
#if defined(MESH) || defined(GRID)
  if (allocated(CLI_numericsFile)) &
    config_numerics => parse(CLI_numericsFile,'numerics configuration')
#else
  MSCMarc: block
    logical :: exists
    inquire(file='numerics.yaml',exist=exists)
    if (exists) config_numerics => parse('numerics.yaml','numerics configuration')
  end block MSCMarc
#endif

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

  type(tDict), intent(in) :: config
  integer, intent(in), optional :: indent
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
      references = references//IO_EOL//filler//'- '//IO_wrapLines(ref%get_asStr(r),filler=filler//'  ')
    end do
  end if

end function config_listReferences


!--------------------------------------------------------------------------------------------------
!> @brief Read configuration, spread over all processes, and add to DADF5.
!--------------------------------------------------------------------------------------------------
function parse(fname,description)

  character(len=*), intent(in) :: fname, description
  type(tDict), pointer :: parse

  character(len=:), allocatable :: fileContent


  if (worldrank == 0) then
    print'(/,1x,a)', 'reading '//description; flush(IO_STDOUT)
    fileContent = IO_read(fname)
    call result_openJobFile(parallel=.false.)
    call result_addSetupFile(fileContent,fname,description)
    call result_closeJobFile()
  end if
  call parallelization_bcast_str(fileContent)

  parse => YAML_str_asDict(fileContent)

end function parse

end module config
