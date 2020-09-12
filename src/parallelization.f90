!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Inquires variables related to parallelization (openMP, MPI)
!--------------------------------------------------------------------------------------------------
module parallelization
  use prec

#ifdef PETSc
#include <petsc/finclude/petscsys.h>
   use petscsys
#endif
!$ use OMP_LIB

  implicit none
  private

  public :: &
    parallelization_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief calls subroutines that reads material, numerics and debug configuration files
!--------------------------------------------------------------------------------------------------
subroutine parallelization_init

  write(6,'(/,a)') ' <<<+-  parallelization init  -+>>>'; flush(6)
  
end subroutine parallelization_init

end module parallelization
