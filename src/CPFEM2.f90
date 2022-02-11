!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief needs a good name and description
!--------------------------------------------------------------------------------------------------
module CPFEM2
  use parallelization
  use DAMASK_interface
  use prec
  use IO
  use YAML_types
  use YAML_parse
  use HDF5
  use HDF5_utilities
  use results
  use config
  use math
  use rotations
  use polynomials
  use lattice
  use material
  use phase
  use homogenization

  use discretization
#if   defined(MESH)
  use FEM_quadrature
  use discretization_mesh
#elif defined(GRID)
  use base64
  use discretization_grid
#endif

  implicit none
  public

contains


!--------------------------------------------------------------------------------------------------
!> @brief Initialize all modules.
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_initAll

  call parallelization_init
  call DAMASK_interface_init                                                                        ! Spectral and FEM interface to commandline
  call prec_init
  call IO_init
#if   defined(MESH)
  call FEM_quadrature_init
#elif defined(GRID)
   call base64_init
#endif
  call YAML_types_init
  call YAML_parse_init
  call HDF5_utilities_init
  call results_init(restart=interface_restartInc>0)
  call config_init
  call math_init
  call rotations_init
  call polynomials_init
  call lattice_init
#if   defined(MESH)
  call discretization_mesh_init(restart=interface_restartInc>0)
#elif defined(GRID)
  call discretization_grid_init(restart=interface_restartInc>0)
#endif
  call material_init(restart=interface_restartInc>0)
  call phase_init
  call homogenization_init
  call CPFEM_init
  call config_deallocate

end subroutine CPFEM_initAll


!--------------------------------------------------------------------------------------------------
!> @brief Read restart information if needed.
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_init

  integer(HID_T) :: fileHandle


  print'(/,1x,a)', '<<<+-  CPFEM init  -+>>>'; flush(IO_STDOUT)


  if (interface_restartInc > 0) then
    print'(/,a,i0,a)', ' reading restart information of increment from file'; flush(IO_STDOUT)

    fileHandle = HDF5_openFile(getSolverJobName()//'_restart.hdf5','r')

    call homogenization_restartRead(fileHandle)
    call phase_restartRead(fileHandle)

    call HDF5_closeFile(fileHandle)
  endif

end subroutine CPFEM_init


!--------------------------------------------------------------------------------------------------
!> @brief Write restart information.
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_restartWrite

  integer(HID_T) :: fileHandle


  print*, ' writing field and constitutive data required for restart to file';flush(IO_STDOUT)

  fileHandle = HDF5_openFile(getSolverJobName()//'_restart.hdf5','a')

  call homogenization_restartWrite(fileHandle)
  call phase_restartWrite(fileHandle)

  call HDF5_closeFile(fileHandle)

end subroutine CPFEM_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief Forward data for new time increment.
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_forward

  call homogenization_forward
  call phase_forward

end subroutine CPFEM_forward


!--------------------------------------------------------------------------------------------------
!> @brief Trigger writing of results.
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_results(inc,time)

  integer,     intent(in) :: inc
  real(pReal), intent(in) :: time

  call results_openJobFile
  call results_addIncrement(inc,time)
  call phase_results
  call homogenization_results
  call discretization_results
  call results_finalizeIncrement
  call results_closeJobFile

end subroutine CPFEM_results

end module CPFEM2
