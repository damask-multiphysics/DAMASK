!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Utilities used by the FEM solver
!--------------------------------------------------------------------------------------------------
module FEM_utilities
#include <petsc/finclude/petsc.h>
  use PETSc
#ifndef PETSC_HAVE_MPI_F90MODULE_VISIBILITY
  use MPI_f08
#endif

  use prec
  use config
  use math
  use misc
  use IO
  use parallelization
  use discretization_mesh
  use homogenization
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<18)
  use FEM_quadrature
#endif
  use constants

#ifndef PETSC_HAVE_MPI_F90MODULE_VISIBILITY
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

!--------------------------------------------------------------------------------------------------
! derived types
  type, public :: tSolutionState                                                                    !< return type of solution from FEM solver variants
    logical :: converged        = .true.
    logical :: stagConverged    = .true.
    PetscInt :: iterationsNeeded = 0_pPETSCINT
  end type tSolutionState

  type, public :: tMechBC
    integer :: nComponents = 0
    real(pREAL), allocatable, dimension(:) :: dot_u
    logical,     allocatable, dimension(:) :: active
  end type tMechBC

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<23)
  external :: &
    PetscSectionGetFieldComponents, &
    PetscSectionGetFieldDof, &
    PetscSectionGetFieldOffset
#endif

  public :: &
    FEM_utilities_init, &
    utilities_constitutiveResponse, &
    utilities_projectBCValues

contains

!ToDo: use functions in variable call
!--------------------------------------------------------------------------------------------------
!> @brief Allocate all neccessary fields.
!--------------------------------------------------------------------------------------------------
subroutine FEM_utilities_init(num_mesh)

  type(tDict), pointer, intent(in) :: &
    num_mesh

  type(tDict), pointer :: &
    num_mech
  character(len=:), allocatable :: &
    petsc_options
  integer :: &
    p_s, &                                                                                          !< order of shape functions
    p_i                                                                                             !< integration order (quadrature rule)
  PetscErrorCode  :: err_PETSc


  print'(/,1x,a)',   '<<<+-  FEM_utilities init  -+>>>'

  num_mech => num_mesh%get_dict('mechanical', defaultVal=emptyDict)

  p_s = num_mesh%get_asInt('p_s',defaultVal = 2)
  p_i = num_mesh%get_asInt('p_i',defaultVal = p_s)

#if (PETSC_VERSION_MINOR>=18)
  if (p_s < 1) &
#else
  if (p_s < 1 .or. p_s > size(FEM_nQuadrature,2)) &
#endif
    call IO_error(301,ext_msg='shape function order (p_s) out of bounds')
  if (p_i < max(1,p_s-1) .or. p_i > p_s) &
    call IO_error(301,ext_msg='integration order (p_i) out of bounds')

  flush(IO_STDOUT)

  petsc_options = misc_prefixOptions('-snes_type newtonls &
                                     &-ksp_type gmres -ksp_max_it 25 -pc_type eisenstat &
                                     &-snes_ksp_ew -snes_ksp_ew_rtol0 0.01 -snes_ksp_ew_rtolmax 0.01 &
                                     &-petscspace_degree ' // IO_intAsStr(p_s) // ' &
                                     &-petscdualspace_lagrange_node_type equispaced &
                                     &-petscdualspace_lagrange_node_endpoints 1 '// &
                                     num_mech%get_asStr('PETSc_options',defaultVal=''),&
                                     'mechanical_')
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,petsc_options,err_PETSc)
  CHKERRQ(err_PETSc)

  call PetscOptionsSetValue(PETSC_NULL_OPTIONS,'-petscds_force_quad','0',err_PETSc)
  CHKERRQ(err_PETSc)

  wgt = real(mesh_maxNips*mesh_nElems,pREAL)**(-1)

end subroutine FEM_utilities_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates constitutive response
!--------------------------------------------------------------------------------------------------
subroutine utilities_constitutiveResponse(status, Delta_t,P_av,forwardData)

  integer(kind(STATUS_OK)),  intent(out)  :: status
  real(pREAL), intent(in)                 :: Delta_t                                                !< loading time
  logical,     intent(in)                 :: forwardData                                            !< age results
  real(pREAL),intent(out), dimension(3,3) :: P_av                                                   !< average PK stress

  integer(MPI_INTEGER_KIND) :: err_MPI


  print'(/,1x,a)', '... evaluating constitutive response ......................................'

  call homogenization_mechanical_response(status,Delta_t,1,int(mesh_maxNips*mesh_nElems))         ! calculate P field
  cutBack = .false.

  P_av = sum(homogenization_P,dim=3) * wgt
  call MPI_Allreduce(MPI_IN_PLACE,P_av,9_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)


end subroutine utilities_constitutiveResponse


!--------------------------------------------------------------------------------------------------
!> @brief Project BC values to local vector
!--------------------------------------------------------------------------------------------------
subroutine utilities_projectBCValues(dm_local,solution_local,section,mechBC,Delta_t,dimPlex)
  DM  :: dm_local
  Vec :: solution_local
  PetscSection   :: section
  type(tMechBC),  dimension(:), intent(in) :: &
    mechBC
  real(pREAL),    intent(in) :: Delta_t
  PetscInt, intent(in) :: dimPlex

  PetscInt       :: component, boundary, bcSize, &
                    nBcPoints, point, dof, numDof, numComp, offset
  IS             :: bcPointsIS
  PetscErrorCode :: err_PETSc
  PetscInt,    pointer :: bcPoints(:)
  real(pREAL), pointer :: localArray(:)

  character(len=11) :: setLabel


  do boundary = 1_pPETSCINT, mesh_Nboundaries; do component = 1_pPETSCINT, dimPlex
   if (mechBC(boundary)%active(component)) then
     setLabel = mesh_BCTypeLabel(mesh_boundariesIdx(boundary))
     call DMGetStratumSize(dm_local,setLabel,mesh_boundariesIS(boundary),bcSize,err_PETSc)
     if (bcSize > 0_pPETSCINT) then
       call DMGetStratumIS(dm_local,setLabel,mesh_boundariesIS(boundary),bcPointsIS,err_PETSc)
       CHKERRQ(err_PETSc)
       call PetscSectionGetFieldComponents(section,0_pPETSCINT,numComp,err_PETSc)
       CHKERRQ(err_PETSc)
       call ISGetSize(bcPointsIS,nBcPoints,err_PETSc)
       CHKERRQ(err_PETSc)
       if (nBcPoints > 0_pPETSCINT) call ISGetIndices(bcPointsIS,bcPoints,err_PETSc)
       call VecGetArray(solution_local,localArray,err_PETSc)
       CHKERRQ(err_PETSc)
       do point = 1_pPETSCINT, nBcPoints
         call PetscSectionGetFieldDof(section,bcPoints(point),0_pPETSCINT,numDof,err_PETSc)
         CHKERRQ(err_PETSc)
         call PetscSectionGetFieldOffset(section,bcPoints(point),0_pPETSCINT,offset,err_PETSc)
         CHKERRQ(err_PETSc)
         do dof = offset+component, offset+numDof, numComp
           localArray(dof) = localArray(dof) + mechBC(boundary)%dot_u(component)*Delta_t
         end do
       end do
       call VecRestoreArray(solution_local,localArray,err_PETSc)
       CHKERRQ(err_PETSc)
       call VecAssemblyBegin(solution_local, err_PETSc)
       CHKERRQ(err_PETSc)
       call VecAssemblyEnd(solution_local, err_PETSc)
       CHKERRQ(err_PETSc)
       if (nBcPoints > 0_pPETSCINT) call ISRestoreIndices(bcPointsIS,bcPoints,err_PETSc)
       CHKERRQ(err_PETSc)
       call ISDestroy(bcPointsIS,err_PETSc)
       CHKERRQ(err_PETSc)
     end if
   end if
  end do; end do

end subroutine utilities_projectBCValues

end module FEM_utilities
