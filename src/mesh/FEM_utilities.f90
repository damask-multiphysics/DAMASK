!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Utilities used by the FEM solver
!--------------------------------------------------------------------------------------------------
module FEM_utilities
#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscis.h>

  use PETScdmplex
  use PETScdmda
  use PETScis
  
  use prec
  use FEsolving
  use homogenization
  use numerics
  use debug
  use math
  use mesh

  implicit none
  private

!--------------------------------------------------------------------------------------------------
  logical,       public       :: cutBack = .false.                                                  !< cut back of BVP solver in case convergence is not achieved or a material point is terminally ill
  integer, public, parameter  :: maxFields = 6
  integer, public             :: nActiveFields = 0
 
!--------------------------------------------------------------------------------------------------
! grid related information information
  real(pReal),   public       :: wgt                                                                !< weighting factor 1/Nelems
  

!--------------------------------------------------------------------------------------------------
! field labels information
  character(len=*),                         parameter,            public :: &
    FIELD_MECH_label     = 'mechanical'
 
  enum, bind(c)
    enumerator :: FIELD_UNDEFINED_ID, &
                  FIELD_MECH_ID
  end enum
  enum, bind(c)
    enumerator :: COMPONENT_UNDEFINED_ID, &
                  COMPONENT_MECH_X_ID, &
                  COMPONENT_MECH_Y_ID, &
                  COMPONENT_MECH_Z_ID
  end enum
 
!--------------------------------------------------------------------------------------------------
! variables controlling debugging
 logical :: &
   debugPETSc                                                                                       !< use some in debug defined options for more verbose PETSc solution

!--------------------------------------------------------------------------------------------------
! derived types
  type, public :: tSolutionState                                                                    !< return type of solution from FEM solver variants
    logical       :: converged         = .true.   
    logical       :: stagConverged     = .true.   
    integer :: iterationsNeeded  = 0
  end type tSolutionState
 
  type, public :: tComponentBC
    integer(kind(COMPONENT_UNDEFINED_ID))          :: ID
    real(pReal),                       allocatable :: Value(:)
    logical,                           allocatable :: Mask(:)   
  end type tComponentBC
 
  type, public :: tFieldBC
    integer(kind(FIELD_UNDEFINED_ID))  :: ID
    integer                            :: nComponents = 0
    type(tComponentBC),    allocatable :: componentBC(:)
  end type tFieldBC
 
  type, public :: tLoadCase
    real(pReal)  :: time                   = 0.0_pReal                                              !< length of increment
    integer      :: incs                   = 0, &                                                   !< number of increments
                    outputfrequency        = 1, &                                                   !< frequency of result writes
                    restartfrequency       = 0, &                                                   !< frequency of restart writes
                    logscale               = 0                                                      !< linear/logarithmic time inc flag
    logical      :: followFormerTrajectory = .true.                                                 !< follow trajectory of former loadcase
    integer,        allocatable :: faceID(:)
    type(tFieldBC),       allocatable :: fieldBC(:)
  end type tLoadCase
  
  public :: &
    utilities_init, &
    utilities_constitutiveResponse, &
    utilities_projectBCValues, &
    FIELD_MECH_ID, &
    COMPONENT_MECH_X_ID, &
    COMPONENT_MECH_Y_ID, &
    COMPONENT_MECH_Z_ID

contains 

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, sets debug flags
!--------------------------------------------------------------------------------------------------
subroutine utilities_init
   
  character(len=1024)                :: petsc_optionsPhysics
  PetscErrorCode                     :: ierr

  write(6,'(/,a)')   ' <<<+-  DAMASK_FEM_utilities init  -+>>>'
 
!--------------------------------------------------------------------------------------------------
! set debugging parameters
  debugPETSc      = iand(debug_level(debug_SPECTRAL),debug_SPECTRALPETSC)      /= 0
  if(debugPETSc) write(6,'(3(/,a),/)') &
                 ' Initializing PETSc with debug options: ', &
                 trim(PETScDebug), &
                 ' add more using the PETSc_Options keyword in numerics.config '
  flush(6)
  call PetscOptionsClear(PETSC_NULL_OPTIONS,ierr)
  CHKERRQ(ierr)
  if(debugPETSc) call PetscOptionsInsertString(PETSC_NULL_OPTIONS,trim(PETSCDEBUG),ierr)
  CHKERRQ(ierr)
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,trim(petsc_defaultOptions),ierr)
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,trim(petsc_options),ierr)
  CHKERRQ(ierr)
  write(petsc_optionsPhysics,'(a,i0)') '-mechFE_petscspace_degree '   , structOrder
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,trim(petsc_optionsPhysics),ierr)
  CHKERRQ(ierr)
  
  wgt = 1.0/real(mesh_maxNips*mesh_NcpElemsGlobal,pReal)


end subroutine utilities_init

!--------------------------------------------------------------------------------------------------
!> @brief calculates constitutive response
!--------------------------------------------------------------------------------------------------
subroutine utilities_constitutiveResponse(timeinc,P_av,forwardData)
 
  real(pReal), intent(in)                 :: timeinc                                                !< loading time
  logical,     intent(in)                 :: forwardData                                            !< age results
  
  real(pReal),intent(out), dimension(3,3) :: P_av                                                   !< average PK stress
  
  PetscErrorCode :: ierr

  write(6,'(/,a)') ' ... evaluating constitutive response ......................................'

  call materialpoint_stressAndItsTangent(.true.,timeinc)                                            ! calculate P field

  restartWrite = .false.                                                                            ! reset restartWrite status
  cutBack = .false.                                                                                 ! reset cutBack status
  
  P_av = sum(sum(materialpoint_P,dim=4),dim=3) * wgt                                                ! average of P 
  call MPI_Allreduce(MPI_IN_PLACE,P_av,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)

end subroutine utilities_constitutiveResponse


!--------------------------------------------------------------------------------------------------
!> @brief Project BC values to local vector
!--------------------------------------------------------------------------------------------------
subroutine utilities_projectBCValues(localVec,section,field,comp,bcPointsIS,BCValue,BCDotValue,timeinc)

  Vec                  :: localVec
  PetscInt             :: field, comp, nBcPoints, point, dof, numDof, numComp, offset
  PetscSection         :: section
  IS                   :: bcPointsIS
  PetscInt,    pointer :: bcPoints(:)
  PetscScalar, pointer :: localArray(:)
  PetscScalar          :: BCValue,BCDotValue,timeinc
  PetscErrorCode       :: ierr

  call PetscSectionGetFieldComponents(section,field,numComp,ierr); CHKERRQ(ierr)
  call ISGetSize(bcPointsIS,nBcPoints,ierr); CHKERRQ(ierr)
  if (nBcPoints > 0) call ISGetIndicesF90(bcPointsIS,bcPoints,ierr)
  call VecGetArrayF90(localVec,localArray,ierr); CHKERRQ(ierr)
  do point = 1, nBcPoints
    call PetscSectionGetFieldDof(section,bcPoints(point),field,numDof,ierr)
    CHKERRQ(ierr)
    call PetscSectionGetFieldOffset(section,bcPoints(point),field,offset,ierr)
    CHKERRQ(ierr)
    do dof = offset+comp+1, offset+numDof, numComp
      localArray(dof) = localArray(dof) + BCValue + BCDotValue*timeinc
    enddo
  enddo    
  call VecRestoreArrayF90(localVec,localArray,ierr); CHKERRQ(ierr)
  call VecAssemblyBegin(localVec, ierr); CHKERRQ(ierr)
  call VecAssemblyEnd  (localVec, ierr); CHKERRQ(ierr)
  if (nBcPoints > 0) call ISRestoreIndicesF90(bcPointsIS,bcPoints,ierr)
 
end subroutine utilities_projectBCValues

end module FEM_utilities
