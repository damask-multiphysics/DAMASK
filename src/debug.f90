!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Reading in and interpretating the debugging settings for the various modules
!--------------------------------------------------------------------------------------------------
module debug
  use prec
  use IO
  use YAML_types
  use YAML_parse

  implicit none
  private

  class(tNode), pointer, public :: &
    debug_root

  public :: debug_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief reads in parameters from debug.config and allocates arrays
!--------------------------------------------------------------------------------------------------
subroutine debug_init

  character(len=:), allocatable :: &
    debug_input, &
    debug_inFlow
  logical :: fexist 

  write(6,'(/,a)')   ' <<<+-  debug init  -+>>>'
#ifdef DEBUG
  write(6,'(a)') achar(27)//'[31m <<<+-  DEBUG version  -+>>>'//achar(27)//'[0m'
#endif

  debug_root => emptyDict
  inquire(file='debug.yaml', exist=fexist)
  fileExists: if (fexist) then
    debug_input  = IO_read('debug.yaml') 
    debug_inFlow = to_flow(debug_input)
    debug_root   => parse_flow(debug_inFlow,defaultVal=emptyDict)
  endif fileExists

end subroutine debug_init

end module debug
