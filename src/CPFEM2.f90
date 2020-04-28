!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief needs a good name and description
!--------------------------------------------------------------------------------------------------
module CPFEM2
  use prec
  use numerics
  use debug
  use config
  use FEsolving
  use math
  use rotations
  use material
  use lattice
  use IO
  use DAMASK_interface
  use results
  use discretization
  use HDF5_utilities
  use homogenization
  use constitutive
  use crystallite
#if    defined(Mesh)
  use FEM_quadrature
  use discretization_mesh
#elif defined(Grid)
  use discretization_grid
#endif

  implicit none
  public

contains


!--------------------------------------------------------------------------------------------------
!> @brief call all module initializations
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_initAll

  call DAMASK_interface_init                                                                        ! Spectral and FEM interface to commandline
  call prec_init
  call IO_init
#ifdef Mesh
  call FEM_quadrature_init
#endif
  call numerics_init
  call debug_init
  call config_init
  call math_init
  call rotations_init
  call lattice_init
  call HDF5_utilities_init
  call results_init
#if    defined(Mesh)
  call discretization_mesh_init
#elif defined(Grid)
  call discretization_grid_init
#endif
  call material_init
  call constitutive_init
  call crystallite_init
  call homogenization_init
  call CPFEM_init

end subroutine CPFEM_initAll


!--------------------------------------------------------------------------------------------------
!> @brief Read restart information if needed.
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_init

  write(6,'(/,a)') ' <<<+-  CPFEM init  -+>>>'; flush(6)

  if (interface_restartInc > 0) call crystallite_restartRead

end subroutine CPFEM_init


!--------------------------------------------------------------------------------------------------
!> @brief Write restart information.
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_restartWrite

  call crystallite_restartWrite

end subroutine CPFEM_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief Forward data for new time increment.
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_forward

  call crystallite_forward

end subroutine CPFEM_forward


!--------------------------------------------------------------------------------------------------
!> @brief Trigger writing of results.
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_results(inc,time)

  integer,     intent(in) :: inc
  real(pReal), intent(in) :: time

  call results_openJobFile
  call results_addIncrement(inc,time)
  call constitutive_results
  call crystallite_results
  call homogenization_results
  call discretization_results
  call results_finalizeIncrement
  call results_closeJobFile

end subroutine CPFEM_results

end module CPFEM2
