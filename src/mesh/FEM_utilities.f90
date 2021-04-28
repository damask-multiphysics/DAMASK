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
  use config
  use math
  use IO
  use discretization_mesh
  use homogenization

  implicit none
  private

!--------------------------------------------------------------------------------------------------
  logical, public             :: cutBack = .false.                                                  !< cut back of BVP solver in case convergence is not achieved or a material point is terminally ill
  integer, public, parameter  :: maxFields = 6
  integer, public             :: nActiveFields = 0
 
!--------------------------------------------------------------------------------------------------
! grid related information information
  real(pReal),   public       :: wgt                                                                !< weighting factor 1/Nelems
  

!--------------------------------------------------------------------------------------------------
! field labels information
  character(len=*),                         parameter,            public :: &
    FIELD_MECH_label     = 'mechanical'
 
  enum, bind(c); enumerator :: &
    FIELD_UNDEFINED_ID, &
    FIELD_MECH_ID
  end enum
  enum, bind(c); enumerator :: &
    COMPONENT_UNDEFINED_ID, &
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
    logical :: converged        = .true.   
    logical :: stagConverged    = .true.   
    integer :: iterationsNeeded = 0
  end type tSolutionState
 
  type, public :: tComponentBC
    integer(kind(COMPONENT_UNDEFINED_ID))          :: ID
    real(pReal),                       allocatable, dimension(:) :: Value
    logical,                           allocatable, dimension(:) :: Mask 
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
                    logscale               = 0                                                      !< linear/logarithmic time inc flag
    logical      :: followFormerTrajectory = .true.                                                 !< follow trajectory of former loadcase
    integer,        allocatable, dimension(:) :: faceID
    type(tFieldBC), allocatable, dimension(:) :: fieldBC
  end type tLoadCase
  
  public :: &
    FEM_utilities_init, &
    utilities_constitutiveResponse, &
    utilities_projectBCValues, &
    FIELD_MECH_ID, &
    COMPONENT_UNDEFINED_ID, &
    COMPONENT_MECH_X_ID, &
    COMPONENT_MECH_Y_ID, &
    COMPONENT_MECH_Z_ID

contains 

!ToDo: use functions in variable call
!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, sets debug flags
!--------------------------------------------------------------------------------------------------
subroutine FEM_utilities_init
   
  character(len=pStringLen) :: petsc_optionsOrder
  class(tNode), pointer :: &
    num_mesh, &
    debug_mesh                                                                                      ! pointer to mesh debug options
  integer :: structOrder                                                                            !< order of displacement shape functions
  character(len=*), parameter :: &
    PETSCDEBUG = ' -snes_view -snes_monitor '

  PetscErrorCode            :: ierr

  print'(/,a)',   ' <<<+-  FEM_utilities init  -+>>>'
 
  num_mesh    => config_numerics%get('mesh',defaultVal=emptyDict)
  structOrder =  num_mesh%get_asInt('structOrder', defaultVal = 2)

  debug_mesh  => config_debug%get('mesh',defaultVal=emptyList)
  debugPETSc  =  debug_mesh%contains('PETSc')

  if(debugPETSc) print'(3(/,a),/)', &
                 ' Initializing PETSc with debug options: ', &
                 trim(PETScDebug), &
                 ' add more using the "PETSc_options" keyword in numerics.yaml'
  flush(IO_STDOUT)
  call PetscOptionsClear(PETSC_NULL_OPTIONS,ierr)
  CHKERRQ(ierr)
  if(debugPETSc) call PetscOptionsInsertString(PETSC_NULL_OPTIONS,trim(PETSCDEBUG),ierr)
  CHKERRQ(ierr)
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,'-mechanical_snes_type newtonls &
                               &-mechanical_snes_linesearch_type cp -mechanical_snes_ksp_ew &
                               &-mechanical_snes_ksp_ew_rtol0 0.01 -mechanical_snes_ksp_ew_rtolmax 0.01 &
                               &-mechanical_ksp_type fgmres -mechanical_ksp_max_it 25 &
                               &-mechanical_pc_type ml -mechanical_mg_levels_ksp_type chebyshev &
                               &-mechanical_mg_levels_pc_type sor -mechanical_pc_ml_nullspace user',ierr)
  CHKERRQ(ierr)
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,num_mesh%get_asString('PETSc_options',defaultVal=''),ierr)
  CHKERRQ(ierr)
  write(petsc_optionsOrder,'(a,i0)') '-mechFE_petscspace_degree ', structOrder
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,trim(petsc_optionsOrder),ierr)
  CHKERRQ(ierr)
  
  wgt = 1.0/real(mesh_maxNips*mesh_NcpElemsGlobal,pReal)


end subroutine FEM_utilities_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates constitutive response
!--------------------------------------------------------------------------------------------------
subroutine utilities_constitutiveResponse(timeinc,P_av,forwardData)
 
  real(pReal), intent(in)                 :: timeinc                                                !< loading time
  logical,     intent(in)                 :: forwardData                                            !< age results
  
  real(pReal),intent(out), dimension(3,3) :: P_av                                                   !< average PK stress
  
  PetscErrorCode :: ierr

  print'(/,a)', ' ... evaluating constitutive response ......................................'

  call materialpoint_stressAndItsTangent(timeinc,[1,mesh_maxNips],[1,mesh_NcpElems])                ! calculate P field

  cutBack = .false.                                                                                 ! reset cutBack status
  
  P_av = sum(homogenization_P,dim=3) * wgt
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
