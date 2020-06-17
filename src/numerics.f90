!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Managing of parameters related to numerics
!--------------------------------------------------------------------------------------------------
module numerics
 use prec
 use IO
 use YAML_types
 use YAML_parse

#ifdef PETSc
#include <petsc/finclude/petscsys.h>
   use petscsys
#endif
!$ use OMP_LIB

 implicit none
 private

 class(tNode), pointer, public :: &
   numerics_root
 integer, protected, public    :: &
   worldrank                  =  0, &                                                               !< MPI worldrank (/=0 for MPI simulations only)
   worldsize                  =  1                                                                  !< MPI worldsize (/=1 for MPI simulations only)
 integer(4), protected, public :: &
   DAMASK_NumThreadsInt       =  0                                                                  !< value stored in environment variable DAMASK_NUM_THREADS, set to zero if no OpenMP directive

!--------------------------------------------------------------------------------------------------
! field parameters:
 integer, protected, public :: &
   itmax                      =  250, &                                                             !< maximum number of iterations
   itmin                      =  1                                                                  !< minimum number of iterations

!--------------------------------------------------------------------------------------------------
! spectral parameters:
#ifdef Grid
 character(len=pStringLen), protected, public :: &
   petsc_options              = ''
#endif

!--------------------------------------------------------------------------------------------------
! Mesh parameters:
#ifdef Mesh
 character(len=pStringLen), protected, public :: &
   petsc_options           = ''
#endif

 public :: numerics_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief reads in parameters from numerics.config and sets openMP related parameters. Also does
! a sanity check
!--------------------------------------------------------------------------------------------------
subroutine numerics_init
!$ integer ::                                gotDAMASK_NUM_THREADS = 1
 integer :: i, ierr
 character(len=:), allocatable :: &
   numerics_input, &
   numerics_inFlow, &
   key
 class (tNode), pointer :: &
   num_grid
 logical :: fexist
!$ character(len=6) DAMASK_NumThreadsString                                                         ! environment variable DAMASK_NUM_THREADS

#ifdef PETSc
 call MPI_Comm_rank(PETSC_COMM_WORLD,worldrank,ierr);CHKERRQ(ierr)
 call MPI_Comm_size(PETSC_COMM_WORLD,worldsize,ierr);CHKERRQ(ierr)
#endif
 write(6,'(/,a)') ' <<<+-  numerics init  -+>>>'

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
 
 fileExists: if (fexist) then
   write(6,'(a,/)') ' using values from config file'
   flush(6)
   numerics_input =  IO_read('numerics.yaml')
   numerics_inFlow = to_flow(numerics_input)
   numerics_root =>  parse_flow(numerics_inFlow,defaultVal=emptyDict)
!--------------------------------------------------------------------------------------------------
! spectral parameters
   num_grid => numerics_root%get('grid',defaultVal=emptyDict)
   do i=1,num_grid%length
     key = num_grid%getKey(i)
     select case(key) 
       case ('itmax')
         itmax = num_grid%get_asInt(key)
       case ('itmin')
         itmin = num_grid%get_asInt(key)
#ifdef PETSC
       case ('petsc_options')
         petsc_options = num_grid%get_asString(key)
#endif 
     endselect
   enddo
 else fileExists
   write(6,'(a,/)') ' using standard values'
   flush(6)
 endif fileExists

!--------------------------------------------------------------------------------------------------
! openMP parameter
 !$  write(6,'(a24,1x,i8,/)')   ' number of threads:      ',DAMASK_NumThreadsInt

!--------------------------------------------------------------------------------------------------
! field parameters
 write(6,'(a24,1x,i8)')      ' itmax:                  ',itmax
 write(6,'(a24,1x,i8)')      ' itmin:                  ',itmin

!--------------------------------------------------------------------------------------------------
#ifdef PETSC
 write(6,'(a24,1x,a)')       ' PETSc_options:          ',trim(petsc_options)
#endif

!--------------------------------------------------------------------------------------------------
! sanity checks
 if (itmax <= 1)                           call IO_error(301,ext_msg='itmax')
 if (itmin > itmax .or. itmin < 1)         call IO_error(301,ext_msg='itmin')

end subroutine numerics_init

end module numerics
