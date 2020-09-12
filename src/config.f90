!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Reads in the material, numerics & debug configuration from their respective file
!> @details Reads the material configuration file, where solverJobName.yaml takes
!! precedence over material.yaml. 
!--------------------------------------------------------------------------------------------------
module config
  use prec
  use DAMASK_interface
  use IO
  use YAML_parse
  use YAML_types

#ifdef PETSc
#include <petsc/finclude/petscsys.h>
   use petscsys
#endif
!$ use OMP_LIB

  implicit none
  private

  class(tNode), pointer, public :: &
    material_root, &
    numerics_root, &
    debug_root

  integer, protected, public :: &
    worldrank                  = 0, &                                                               !< MPI worldrank (/=0 for MPI simulations only)
    worldsize                  = 1                                                                  !< MPI worldsize (/=1 for MPI simulations only)
  integer(4), protected, public :: &
    DAMASK_NumThreadsInt       = 0                                                                  !< value stored in environment variable DAMASK_NUM_THREADS, set to zero if no OpenMP directive

  public :: &
    config_init, &
    config_deallocate

contains

!--------------------------------------------------------------------------------------------------
!> @brief calls subroutines that reads material, numerics and debug configuration files
!--------------------------------------------------------------------------------------------------
subroutine config_init

  write(6,'(/,a)') ' <<<+-  config init  -+>>>'; flush(6)
  
  call parse_material
  call parse_numerics
  call parse_debug

end subroutine config_init


!--------------------------------------------------------------------------------------------------
!> @brief reads material.yaml
!--------------------------------------------------------------------------------------------------
subroutine parse_material

  logical :: fileExists
  character(len=:), allocatable :: fname

  fname = getSolverJobName()//'.yaml'
  inquire(file=fname,exist=fileExists)
  if(.not. fileExists) then
    fname = 'material.yaml'
    inquire(file=fname,exist=fileExists)
    if(.not. fileExists) call IO_error(100,ext_msg=fname)
  endif
  write(6,*) 'reading '//fname; flush(6)
  material_root => parse_flow(to_flow(IO_read(fname)))

end subroutine parse_material


!--------------------------------------------------------------------------------------------------
!> @brief reads in parameters from numerics.yaml and sets openMP related parameters. Also does
! a sanity check
!--------------------------------------------------------------------------------------------------
subroutine parse_numerics

  logical :: fexist
  integer :: ierr
!$ integer :: gotDAMASK_NUM_THREADS = 1
!$ character(len=6) DAMASK_NumThreadsString                                                         ! environment variable DAMASK_NUM_THREADS

#ifdef PETSc
  call MPI_Comm_rank(PETSC_COMM_WORLD,worldrank,ierr);CHKERRQ(ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,worldsize,ierr);CHKERRQ(ierr)
#endif

!$ call GET_ENVIRONMENT_VARIABLE(NAME='DAMASK_NUM_THREADS',VALUE=DAMASK_NumThreadsString,STATUS=gotDAMASK_NUM_THREADS)   ! get environment variable DAMASK_NUM_THREADS...
!$ if(gotDAMASK_NUM_THREADS /= 0) then                                                              ! could not get number of threads, set it to 1
!$   call IO_warning(35,ext_msg='BEGIN:'//DAMASK_NumThreadsString//':END')
!$   DAMASK_NumThreadsInt = 1_4
!$ else
!$   read(DAMASK_NumThreadsString,'(i6)') DAMASK_NumThreadsInt                                      ! read as integer
!$   if (DAMASK_NumThreadsInt < 1_4) DAMASK_NumThreadsInt = 1_4                                     ! in case of string conversion fails, set it to one
!$ endif
!$ call omp_set_num_threads(DAMASK_NumThreadsInt)                                                   ! set number of threads for parallel execution

  numerics_root => emptyDict
  inquire(file='numerics.yaml', exist=fexist)
  if (fexist) then
    write(6,*) 'reading numerics.yaml'; flush(6)
    numerics_root =>  parse_flow(to_flow(IO_read('numerics.yaml')))
  endif

!--------------------------------------------------------------------------------------------------
! openMP parameter
 !$  write(6,'(a24,1x,i8,/)')   ' number of threads:      ',DAMASK_NumThreadsInt

end subroutine parse_numerics


!--------------------------------------------------------------------------------------------------
!> @brief reads in parameters from debug.yaml
!--------------------------------------------------------------------------------------------------
subroutine parse_debug

  logical :: fexist 

  debug_root => emptyDict
  inquire(file='debug.yaml', exist=fexist)
  fileExists: if (fexist) then
    write(6,*) 'reading debug.yaml'; flush(6)
    debug_root  => parse_flow(to_flow(IO_read('debug.yaml')))
  endif fileExists

end subroutine parse_debug


!--------------------------------------------------------------------------------------------------
!> @brief deallocates material.yaml structure
!--------------------------------------------------------------------------------------------------
subroutine config_deallocate

  deallocate(material_root)                                                                         !ToDo: deallocation of numerics debug (optional)

end subroutine config_deallocate

end module config
