!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief needs a good name and description
!--------------------------------------------------------------------------------------------------
module materialpoint
  use parallelization
  use CLI
  use OS
  use signal
  use prec
  use misc
  use IO
  use types
  use YAML
  use HDF5
  use HDF5_utilities
  use result
  use config
  use math
  use rotations
  use polynomials
  use tables
  use crystal
  use material
  use phase
  use homogenization
  use discretization
#if   defined(MESH)
  use discretization_mesh
#elif defined(GRID)
  use base64
  use zlib
  use discretization_grid
#endif

  implicit none(type,external)
  public

contains


!--------------------------------------------------------------------------------------------------
!> @brief Initialize all modules.
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_initAll()

  call parallelization_init()
  call prec_init()
  call OS_init()
  call misc_init()
  call IO_init()
  call CLI_init()                                                                                   ! grid and mesh commandline interface
  call signal_init()
#if defined(GRID)
   call zlib_init()
   call base64_init()
#endif
  call types_init()
  call YAML_init()
  call HDF5_utilities_init()
  call result_init(restart=CLI_restartInc>0)
  call config_init()
  call math_init()
  call rotations_init()
  call polynomials_init()
  call tables_init()
  call crystal_init()
#if   defined(MESH)
  call discretization_mesh_init()
#elif defined(GRID)
  call discretization_grid_init()
#endif
  call material_init()
  call phase_init()
  call homogenization_init()
  call materialpoint_init()
  call config_material_deallocate()

end subroutine materialpoint_initAll


!--------------------------------------------------------------------------------------------------
!> @brief Read restart information if needed.
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_init()

  integer(HID_T) :: fileHandle


  print'(/,1x,a)', '<<<+-  materialpoint init  -+>>>'; flush(IO_STDOUT)

  if (CLI_restartInc > 0) then
    print'(/,1x,a,1x,i0)', 'loading restart information of increment',CLI_restartInc; flush(IO_STDOUT)

    fileHandle = HDF5_openFile(CLI_jobName//'_restart.hdf5','r')

    call homogenization_restartRead(fileHandle)
    call phase_restartRead(fileHandle)

    call HDF5_closeFile(fileHandle)
  end if

end subroutine materialpoint_init


!--------------------------------------------------------------------------------------------------
!> @brief Write restart information.
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_restartWrite()

  integer(HID_T) :: fileHandle


  print'(1x,a)', 'saving field and constitutive data required for restart';flush(IO_STDOUT)

  fileHandle = HDF5_openFile(CLI_jobName//'_restart.hdf5','a')

  call homogenization_restartWrite(fileHandle)
  call phase_restartWrite(fileHandle)

  call HDF5_closeFile(fileHandle)

end subroutine materialpoint_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief Forward data for new time increment.
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_forward()

  call homogenization_forward()
  call phase_forward()

end subroutine materialpoint_forward


!--------------------------------------------------------------------------------------------------
!> @brief Trigger writing of results.
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_result(inc,time)

  integer,     intent(in) :: inc
  real(pREAL), intent(in) :: time

  call result_openJobFile()
  call result_addIncrement(inc,time)
  call phase_result()
  call homogenization_result()
  call discretization_result()
  call result_finalizeIncrement()
  call result_closeJobFile()

end subroutine materialpoint_result

end module materialpoint
