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

!$ integer :: got_env, DAMASK_NUM_THREADS
!$ character(len=6) NumThreadsString

  write(6,'(/,a)') ' <<<+-  parallelization init  -+>>>'; flush(6)

!$ call get_environment_variable(name='DAMASK_NUM_THREADS',value=NumThreadsString,STATUS=got_env)   ! get environment variable DAMASK_NUM_THREADS...
!$ if(got_env /= 0) then                                                                            ! could not get number of threads, set it to 1
!$   write(6,*) 'Could not determine value of $DAMASK_NUM_THREADS'
!$   DAMASK_NUM_THREADS = 1_pI32
!$ else
!$   read(NumThreadsString,'(i6)') DAMASK_NUM_THREADS
!$   if (DAMASK_NUM_THREADS < 1_pI32) then
!$     write(6,*) 'Invalid DAMASK_NUM_THREADS: '//trim(NumThreadsString)
!$     DAMASK_NUM_THREADS = 1_pI32
!$   endif
!$ endif
!$ write(6,'(a,i8,/)')   ' DAMASK_NUM_THREADS: ',DAMASK_NUM_THREADS
!$ call omp_set_num_threads(DAMASK_NUM_THREADS)
  
end subroutine parallelization_init

end module parallelization
