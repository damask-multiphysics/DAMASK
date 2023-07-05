!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Utilities used by the FEM solver
!--------------------------------------------------------------------------------------------------
module FEM_utilities
#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscis.h>
  use PETScDMplex
  use PETScDMDA
  use PETScIS
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  use MPI_f08
#endif

  use prec
  use config
  use math
  use IO
  use discretization_mesh
  use homogenization
  use FEM_quadrature

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  implicit none(type,external)
#else
  implicit none
#endif
  private

  logical,     public             :: cutBack = .false.                                              !< cut back of BVP solver in case convergence is not achieved or a material point is terminally ill
  real(pREAL), public, protected  :: wgt                                                            !< weighting factor 1/Nelems


!--------------------------------------------------------------------------------------------------
! field labels information
  character(len=*), parameter, public :: &
    FIELD_MECH_label = 'mechanical'

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
! derived types
  type, public :: tSolutionState                                                                    !< return type of solution from FEM solver variants
    logical :: converged        = .true.
    logical :: stagConverged    = .true.
    PetscInt :: iterationsNeeded = 0_pPETSCINT
  end type tSolutionState

  type, public :: tComponentBC
    integer(kind(COMPONENT_UNDEFINED_ID)) :: ID
    real(pREAL), allocatable, dimension(:) :: Value
    logical,     allocatable, dimension(:) :: Mask
  end type tComponentBC

  type, public :: tFieldBC
    integer(kind(FIELD_UNDEFINED_ID))  :: ID
    integer                            :: nComponents = 0
    type(tComponentBC), allocatable, dimension(:) :: componentBC
  end type tFieldBC

  external :: &                                                                                     ! ToDo: write interfaces
    PetscSectionGetFieldComponents, &
    PetscSectionGetFieldDof, &
    PetscSectionGetFieldOffset

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
!> @brief Allocate all neccessary fields.
!--------------------------------------------------------------------------------------------------
subroutine FEM_utilities_init

  character(len=pSTRLEN) :: petsc_optionsOrder
  type(tDict), pointer :: &
    num_mesh
  integer :: &
    p_s, &                                                                                          !< order of shape functions
    p_i                                                                                             !< integration order (quadrature rule)
  PetscErrorCode :: err_PETSc


  print'(/,1x,a)',   '<<<+-  FEM_utilities init  -+>>>'

  num_mesh => config_numerics%get_dict('mesh',defaultVal=emptyDict)

  p_s = num_mesh%get_asInt('p_s',defaultVal = 2)
  p_i = num_mesh%get_asInt('p_i',defaultVal = p_s)

  if (p_s < 1 .or. p_s > size(FEM_nQuadrature,2)) &
    call IO_error(821,ext_msg='shape function order (p_s) out of bounds')
  if (p_i < max(1,p_s-1) .or. p_i > p_s) &
    call IO_error(821,ext_msg='integration order (p_i) out of bounds')

  flush(IO_STDOUT)
  call PetscOptionsClear(PETSC_NULL_OPTIONS,err_PETSc)
  CHKERRQ(err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,'-mechanical_snes_type newtonls &
                               &-mechanical_snes_linesearch_type cp -mechanical_snes_ksp_ew &
                               &-mechanical_snes_ksp_ew_rtol0 0.01 -mechanical_snes_ksp_ew_rtolmax 0.01 &
                               &-mechanical_ksp_type fgmres -mechanical_ksp_max_it 25', err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,num_mesh%get_asStr('PETSc_options',defaultVal=''),err_PETSc)
  CHKERRQ(err_PETSc)
  write(petsc_optionsOrder,'(a,i0)') '-mechFE_petscspace_degree ', p_s
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,trim(petsc_optionsOrder),err_PETSc)
  CHKERRQ(err_PETSc)

  wgt = real(mesh_maxNips*mesh_NcpElemsGlobal,pREAL)**(-1)


end subroutine FEM_utilities_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates constitutive response
!--------------------------------------------------------------------------------------------------
subroutine utilities_constitutiveResponse(timeinc,P_av,forwardData)

  real(pREAL), intent(in)                 :: timeinc                                                !< loading time
  logical,     intent(in)                 :: forwardData                                            !< age results
  real(pREAL),intent(out), dimension(3,3) :: P_av                                                   !< average PK stress

  integer(MPI_INTEGER_KIND) :: err_MPI

  print'(/,1x,a)', '... evaluating constitutive response ......................................'

  call homogenization_mechanical_response(timeinc,1,mesh_maxNips*mesh_NcpElems)                     ! calculate P field
  if (.not. terminallyIll) &
    call homogenization_mechanical_response2(timeinc,[1,mesh_maxNips],[1,mesh_NcpElems])
  cutBack = .false.

  P_av = sum(homogenization_P,dim=3) * wgt
  call MPI_Allreduce(MPI_IN_PLACE,P_av,9_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'


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
  real(pREAL), pointer :: localArray(:)
  real(pREAL)          :: BCValue,BCDotValue,timeinc
  PetscErrorCode       :: err_PETSc


  call PetscSectionGetFieldComponents(section,field,numComp,err_PETSc)
  CHKERRQ(err_PETSc)
  call ISGetSize(bcPointsIS,nBcPoints,err_PETSc)
  CHKERRQ(err_PETSc)
  if (nBcPoints > 0) call ISGetIndicesF90(bcPointsIS,bcPoints,err_PETSc)
  call VecGetArrayF90(localVec,localArray,err_PETSc)
  CHKERRQ(err_PETSc)
  do point = 1, nBcPoints
    call PetscSectionGetFieldDof(section,bcPoints(point),field,numDof,err_PETSc)
    CHKERRQ(err_PETSc)
    call PetscSectionGetFieldOffset(section,bcPoints(point),field,offset,err_PETSc)
    CHKERRQ(err_PETSc)
    do dof = offset+comp+1, offset+numDof, numComp
      localArray(dof) = localArray(dof) + BCValue + BCDotValue*timeinc
    end do
  end do
  call VecRestoreArrayF90(localVec,localArray,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecAssemblyBegin(localVec, err_PETSc)
  CHKERRQ(err_PETSc)
  call VecAssemblyEnd  (localVec, err_PETSc)
  CHKERRQ(err_PETSc)
  if (nBcPoints > 0) call ISRestoreIndicesF90(bcPointsIS,bcPoints,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine utilities_projectBCValues

end module FEM_utilities
