!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Inquires variables related to parallelization (openMP, MPI)
!--------------------------------------------------------------------------------------------------
module parallelization
  use, intrinsic :: ISO_fortran_env, only: &
    OUTPUT_UNIT

#ifdef PETSc
#include <petsc/finclude/petscsys.h>
   use petscsys
#endif
!$ use OMP_LIB

  use prec

  implicit none
  private

  integer, protected, public :: &
    worldrank = 0, &                                                                                !< MPI worldrank (/=0 for MPI simulations only)
    worldsize = 1                                                                                   !< MPI worldsize (/=1 for MPI simulations only)

  public :: &
    parallelization_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief calls subroutines that reads material, numerics and debug configuration files
!--------------------------------------------------------------------------------------------------
subroutine parallelization_init

  integer :: err, typeSize
!$ integer :: got_env, DAMASK_NUM_THREADS, threadLevel
!$ character(len=6) NumThreadsString
#ifdef PETSc
  PetscErrorCode :: petsc_err

#else
  print'(/,a)', ' <<<+-  parallelization init  -+>>>'; flush(OUTPUT_UNIT)
#endif

#ifdef PETSc
#ifdef _OPENMP
  ! If openMP is enabled, check if the MPI libary supports it and initialize accordingly.
  ! Otherwise, the first call to PETSc will do the initialization.
  call MPI_Init_Thread(MPI_THREAD_FUNNELED,threadLevel,err)
  if (err /= 0)                              error stop 'MPI init failed'
  if (threadLevel<MPI_THREAD_FUNNELED)       error stop 'MPI library does not support OpenMP'
#endif

  call PetscInitializeNoArguments(petsc_err)                                                        ! first line in the code according to PETSc manual
  CHKERRQ(petsc_err)

#if defined(DEBUG) && defined(__INTEL_COMPILER)
  call PetscSetFPTrap(PETSC_FP_TRAP_ON,petsc_err)
#else
  call PetscSetFPTrap(PETSC_FP_TRAP_OFF,petsc_err)
#endif
  CHKERRQ(petsc_err)

call MPI_Comm_rank(PETSC_COMM_WORLD,worldrank,err)
  if (err /= 0)                              error stop 'Could not determine worldrank'

  if (worldrank == 0) print'(/,a)',  ' <<<+-  parallelization init  -+>>>'

  call MPI_Comm_size(PETSC_COMM_WORLD,worldsize,err)
  if (err /= 0)                              error stop 'Could not determine worldsize'
  if (worldrank == 0) print'(a,i3)', ' MPI processes: ',worldsize

  call MPI_Type_size(MPI_INTEGER,typeSize,err)
  if (err /= 0)                              error stop 'Could not determine MPI integer size'
  if (typeSize*8 /= bit_size(0))             error stop 'Mismatch between MPI and DAMASK integer'

  call MPI_Type_size(MPI_DOUBLE,typeSize,err)
  if (err /= 0)                              error stop 'Could not determine MPI real size'
  if (typeSize*8 /= storage_size(0.0_pReal)) error stop 'Mismatch between MPI and DAMASK real'
#endif

  if (worldrank /= 0) then
    close(OUTPUT_UNIT)                                                                              ! disable output
    open(OUTPUT_UNIT,file='/dev/null',status='replace')                                             ! close() alone will leave some temp files in cwd
  endif

!$ call get_environment_variable(name='DAMASK_NUM_THREADS',value=NumThreadsString,STATUS=got_env)
!$ if(got_env /= 0) then
!$   print*, 'Could not determine value of $DAMASK_NUM_THREADS'
!$   DAMASK_NUM_THREADS = 1_pI32
!$ else
!$   read(NumThreadsString,'(i6)') DAMASK_NUM_THREADS
!$   if (DAMASK_NUM_THREADS < 1_pI32) then
!$     print*, 'Invalid DAMASK_NUM_THREADS: '//trim(NumThreadsString)
!$     DAMASK_NUM_THREADS = 1_pI32
!$   endif
!$ endif
!$ print'(a,i2)',   ' DAMASK_NUM_THREADS: ',DAMASK_NUM_THREADS
!$ call omp_set_num_threads(DAMASK_NUM_THREADS)

end subroutine parallelization_init

end module parallelization
