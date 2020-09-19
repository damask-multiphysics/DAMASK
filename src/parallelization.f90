!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Inquires variables related to parallelization (openMP, MPI)
!--------------------------------------------------------------------------------------------------
module parallelization
  use prec
  use, intrinsic :: iso_fortran_env

#ifdef PETSc
#include <petsc/finclude/petscsys.h>
   use petscsys
#endif
!$ use OMP_LIB

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

  integer        :: err, typeSize
!$ integer :: got_env, DAMASK_NUM_THREADS, threadLevel
!$ character(len=6) NumThreadsString
#ifdef PETSc
  PetscErrorCode :: petsc_err

#else
  print'(/,a)', ' <<<+-  parallelization init  -+>>>'; flush(6)
#endif

#ifdef PETSc
#ifdef _OPENMP
  ! If openMP is enabled, check if the MPI libary supports it and initialize accordingly.
  ! Otherwise, the first call to PETSc will do the initialization.
  call MPI_Init_Thread(MPI_THREAD_FUNNELED,threadLevel,err)
  if (err /= 0)                        error stop 'MPI init failed'
  if (threadLevel<MPI_THREAD_FUNNELED) error stop 'MPI library does not support OpenMP'
#endif

  call PETScInitializeNoArguments(petsc_err)                                                        ! first line in the code according to PETSc manual
  CHKERRQ(petsc_err)

  call MPI_Comm_rank(PETSC_COMM_WORLD,worldrank,err)
  if (err /= 0) error stop 'Could not determine worldrank'

  if (worldrank == 0) print'(/,a)',  ' <<<+-  parallelization init  -+>>>'
  if (worldrank == 0) print'(a,i3)', ' MPI processes: ',worldsize

  call MPI_Comm_size(PETSC_COMM_WORLD,worldsize,err)
  if (err /= 0) error stop 'Could not determine worldsize'

  call MPI_Type_size(MPI_INTEGER,typeSize,err)
  if (err /= 0)                              error stop 'Could not determine MPI integer size'
  if (typeSize*8 /= bit_size(0))             error stop 'Mismatch between MPI and DAMASK integer'

  call MPI_Type_size(MPI_DOUBLE,typeSize,err)
  if (err /= 0)                              error stop 'Could not determine MPI real size'
  if (typeSize*8 /= storage_size(0.0_pReal)) error stop 'Mismatch between MPI and DAMASK real'
#endif

  mainProcess: if (worldrank == 0) then
    if (output_unit /= 6) error stop 'STDOUT != 6'
    if (error_unit /= 0)  error stop 'STDERR != 0'
  else mainProcess
    close(6)                                                                                        ! disable output for non-master processes (open 6 to rank specific file for debug)
    open(6,file='/dev/null',status='replace')                                                       ! close(6) alone will leave some temp files in cwd
  endif mainProcess


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
