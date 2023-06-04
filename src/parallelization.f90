!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Inquires variables related to parallelization (openMP, MPI)
!--------------------------------------------------------------------------------------------------
module parallelization
  use, intrinsic :: ISO_fortran_env, only: &
    OUTPUT_UNIT, &
    ERROR_UNIT

#ifdef PETSC
#include <petsc/finclude/petscsys.h>
  use PETScSys
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  use MPI_f08
#endif
!$ use OMP_LIB
#endif

  use prec

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  implicit none(type,external)
#else
  implicit none
#endif
  private

#ifndef PETSC
  integer, parameter, public :: &
    MPI_INTEGER_KIND = pI64                                                                         !< needed for MSC.Marc
  integer(MPI_INTEGER_KIND), parameter, public :: &
#else
  integer(MPI_INTEGER_KIND), protected, public :: &
#endif
    worldrank = 0_MPI_INTEGER_KIND, &                                                               !< MPI worldrank (/=0 for MPI simulations only)
    worldsize = 1_MPI_INTEGER_KIND                                                                  !< MPI worldsize (/=1 for MPI simulations only)

#ifndef PETSC
public :: parallelization_bcast_str

contains
subroutine parallelization_bcast_str(str)
  character(len=:), allocatable, intent(inout) :: str
end subroutine parallelization_bcast_str

#else
  public :: &
    parallelization_init, &
    parallelization_chkerr, &
    parallelization_bcast_str

contains

!--------------------------------------------------------------------------------------------------
!> @brief Initialize shared memory (openMP) and distributed memory (MPI) parallelization.
!--------------------------------------------------------------------------------------------------
subroutine parallelization_init()

  integer(MPI_INTEGER_KIND) :: err_MPI, typeSize, version, subversion, devNull
  character(len=4) :: rank_str
  character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: MPI_library_version
!$ integer :: got_env, threadLevel
!$ integer(pI32) :: OMP_NUM_THREADS
!$ character(len=6) NumThreadsString

  PetscErrorCode :: err_PETSc
#ifdef _OPENMP
  ! If openMP is enabled, check if the MPI libary supports it and initialize accordingly.
  call MPI_Init_Thread(MPI_THREAD_FUNNELED,threadLevel,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND)   error stop 'MPI init failed'
  if (threadLevel<MPI_THREAD_FUNNELED) error stop 'MPI library does not support OpenMP'
#else
  call MPI_Init(err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND)   error stop 'MPI init failed'
#endif

#if defined(DEBUG)
  call PetscInitialize(PETSC_NULL_CHARACTER,err_PETSc)
#else
  call PetscInitializeNoArguments(err_PETSc)
#endif
  CHKERRQ(err_PETSc)

#if defined(DEBUG) && defined(__INTEL_COMPILER)
  call PetscSetFPTrap(PETSC_FP_TRAP_ON,err_PETSc)
#else
  call PetscSetFPTrap(PETSC_FP_TRAP_OFF,err_PETSc)
#endif
  CHKERRQ(err_PETSc)

  call MPI_Comm_rank(MPI_COMM_WORLD,worldrank,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) &
    error stop 'Could not determine worldrank'

#ifdef LOGFILE
  write(rank_str,'(i4.4)') worldrank
  open(OUTPUT_UNIT,file='out.'//rank_str,status='replace',encoding='UTF-8')
  open(ERROR_UNIT,file='err.'//rank_str,status='replace',encoding='UTF-8')
#else
  if (worldrank /= 0) then
    close(OUTPUT_UNIT)                                                                              ! disable output
    open(OUTPUT_UNIT,file='/dev/null',status='replace')                                             ! close() alone will leave some temp files in cwd
  else
    open(OUTPUT_UNIT,encoding='UTF-8')                                                              ! for special characters in output
  end if
#endif

  print'(/,1x,a)', '<<<+-  parallelization init  -+>>>'

  call MPI_Get_library_version(MPI_library_version,devNull,err_MPI)
  print'(/,1x,a)', trim(MPI_library_version)
  call MPI_Get_version(version,subversion,err_MPI)
  print'(1x,a,i0,a,i0)', 'MPI standard: ',version,'.',subversion
#ifdef _OPENMP
  print'(1x,a,i0)',      'OpenMP version: ',openmp_version
#endif

  call MPI_Comm_size(MPI_COMM_WORLD,worldsize,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) &
    error stop 'Could not determine worldsize'
  if (worldrank == 0) print'(/,1x,a,i0)', 'MPI processes: ',worldsize

  call MPI_Type_size(MPI_INTEGER,typeSize,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) &
    error stop 'Could not determine size of MPI_INTEGER'
  if (typeSize*8_MPI_INTEGER_KIND /= int(bit_size(0),MPI_INTEGER_KIND)) &
    error stop 'Mismatch between MPI_INTEGER and DAMASK default integer'

  call MPI_Type_size(MPI_INTEGER8,typeSize,err_MPI)
  if (err_MPI /= 0) &
    error stop 'Could not determine size of MPI_INTEGER8'
  if (typeSize*8_MPI_INTEGER_KIND /= int(bit_size(0_pI64),MPI_INTEGER_KIND)) &
    error stop 'Mismatch between MPI_INTEGER8 and DAMASK pI64'

  call MPI_Type_size(MPI_DOUBLE,typeSize,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) &
    error stop 'Could not determine size of MPI_DOUBLE'
  if (typeSize*8_MPI_INTEGER_KIND /= int(storage_size(0.0_pREAL),MPI_INTEGER_KIND)) &
    error stop 'Mismatch between MPI_DOUBLE and DAMASK pREAL'

!$ call get_environment_variable(name='OMP_NUM_THREADS',value=NumThreadsString,STATUS=got_env)
!$ if (got_env /= 0) then
!$   print'(1x,a)', 'Could not get $OMP_NUM_THREADS, using default'
!$   OMP_NUM_THREADS = 4_pI32
!$ else
!$   read(NumThreadsString,'(i6)') OMP_NUM_THREADS
!$   if (OMP_NUM_THREADS < 1_pI32) then
!$     print'(1x,a)', 'Invalid OMP_NUM_THREADS: "'//trim(NumThreadsString)//'", using default'
!$     OMP_NUM_THREADS = 4_pI32
!$   end if
!$ end if
!$ print'(1x,a,i0)',   'OMP_NUM_THREADS: ',OMP_NUM_THREADS
!$ call omp_set_num_threads(OMP_NUM_THREADS)

end subroutine parallelization_init


!--------------------------------------------------------------------------------------------------
!> @brief Check for MPI error.
!--------------------------------------------------------------------------------------------------
subroutine parallelization_chkerr(e)

  integer(MPI_INTEGER_KIND), intent(in) :: e


  if (e/=0_MPI_INTEGER_KIND) error stop 'MPI error'

end subroutine parallelization_chkerr


!--------------------------------------------------------------------------------------------------
!> @brief Broadcast a string from process 0.
!--------------------------------------------------------------------------------------------------
subroutine parallelization_bcast_str(str)

  character(len=:), allocatable, intent(inout) :: str

  integer(MPI_INTEGER_KIND) :: strlen, err_MPI


  if (worldrank == 0) strlen = len(str,MPI_INTEGER_KIND)
  call MPI_Bcast(strlen,1_MPI_INTEGER_KIND,MPI_INTEGER,0_MPI_INTEGER_KIND,MPI_COMM_WORLD, err_MPI)
  if (worldrank /= 0) allocate(character(len=strlen)::str)

  call MPI_Bcast(str,strlen,MPI_CHARACTER,0_MPI_INTEGER_KIND,MPI_COMM_WORLD, err_MPI)


end subroutine parallelization_bcast_str

#endif

end module parallelization
