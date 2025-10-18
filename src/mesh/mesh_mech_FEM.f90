!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief FEM PETSc solver
!--------------------------------------------------------------------------------------------------
module mesh_mechanical_FEM
#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petsc.h>
  use PETScSNES
  use PETScDM
  use PETScDMplex
#if (PETSC_VERSION_MAJOR==3 && (PETSC_VERSION_MINOR>18 && PETSC_VERSION_MINOR<23))
  use PETScDT
#endif
#ifndef PETSC_HAVE_MPI_F90MODULE_VISIBILITY
  use MPI_f08
#endif

  use prec
  use FEM_utilities
  use discretization
  use discretization_mesh
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<18)
  use FEM_quadrature
#endif
  use homogenization
  use math
  use constants

#ifndef PETSC_HAVE_MPI_F90MODULE_VISIBILITY
  implicit none(type,external)
#else
  implicit none
#endif
  private

!--------------------------------------------------------------------------------------------------
! derived types
  type tSolutionParams
    type(tMechBC),allocatable, dimension(:)  :: mechBC
    real(pREAL)    :: Delta_t
  end type tSolutionParams

  type(tSolutionParams)  :: params

  type, private :: tNumerics
    PetscInt :: &
      p_i, &                                                                                        !< integration order (quadrature rule)
      itmax
    logical :: &
      BBarStabilization
    real(pREAL) :: &
      eps_struct_atol, &                                                                            !< absolute tolerance for mechanical equilibrium
      eps_struct_rtol                                                                               !< relative tolerance for mechanical equilibrium
  end type tNumerics

  type(tNumerics), private :: num
!--------------------------------------------------------------------------------------------------
! PETSc data
  SNES                           :: mechanical_snes
  Vec                            :: solution, solution_rate, solution_local
  PetscInt                       :: dimPlex, cellDof, nBasis
  PetscInt                       :: nQuadrature
  PetscReal, allocatable, target :: qWeights(:)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<18)
  PetscReal, allocatable, target :: qPoints(:)
#endif
  MatNullSpace                   :: matnull

!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
  character(len=pSTRLEN) :: incInfo
  real(pREAL), dimension(3,3) :: &
    P_av = 0.0_pREAL
  logical :: ForwardData
  integer(kind(STATUS_OK)) :: status
  real(pREAL), parameter :: eps = 1.0e-18_pREAL

#if (PETSC_VERSION_MAJOR==3 && (PETSC_VERSION_MINOR>14 && PETSC_VERSION_MINOR<23))
  external :: &                                                                                     ! ToDo: write interfaces
#if (PETSC_VERSION_MINOR<16)
    ISDestroy, &
#endif
#if (PETSC_VERSION_MINOR>18 && PETSC_VERSION_MINOR<22)
    DMAddField, &
#endif
    PetscSectionGetNumFields, &
    PetscFESetQuadrature, &
    PetscFEGetDimension, &
    PetscFEDestroy, &
    PetscSectionGetDof, &
    PetscFEGetDualSpace, &
    PetscDualSpaceGetFunctional
#endif

  public :: &
    FEM_mechanical_init, &
    FEM_mechanical_solution, &
    FEM_mechanical_forward, &
    FEM_mechanical_updateCoords

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data
!--------------------------------------------------------------------------------------------------
subroutine FEM_mechanical_init(mechBC,num_mesh)

  type(tMechBC), dimension(:), intent(in) :: mechBC
  type(tDict),   pointer,      intent(in) :: num_mesh

  DM       :: mechanical_mesh
  IS       :: bcPoint
  DMLabel  :: BCLabel
  PetscFE  :: mechFE
  PetscDS  :: mechDS
  PetscInt :: numActiveBC, bcSize, nc, &
              component, boundary, topologDim, &
              cellStart, cellEnd, &
              nCoords
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
  PetscBool       :: isSimplex
#endif
  PetscSection    :: section
  PetscQuadrature :: mechQuad
  PetscDualSpace  :: mechDualSpace
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
  PetscObject     :: obj
#endif
  PetscErrorCode  :: err_PETSc
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<23)
  DMLabel,   dimension(:), pointer     :: PETSC_NULL_DMLABEL_ARRAY => NULL()
#endif
  IS,        dimension(:), pointer     :: pBcComps, pBcPoints
  PetscInt,  dimension(:), pointer     :: pNumComp, pNumDof, pBcField, pBcPoint
  PetscInt,  dimension(:), allocatable :: idx
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>17)
  PetscReal, dimension(:), pointer     :: qWeightsP
#else
  PetscReal, dimension(:), pointer     :: qPointsP, qWeightsP
#endif
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<23)
  PetscReal, dimension(:), pointer     :: PETSC_NULL_REAL_POINTER => NULL()
#endif
  real(pREAL), dimension(:), allocatable :: nodeCoords

  character(len=*), parameter :: prefix = 'mechanical_'
  character(len=11)           :: setLabel

  real(pREAL), dimension(3,3) :: devNull
  type(tDict), pointer        :: num_mech

  print'(/,1x,a)', '<<<+-  FEM_mech init  -+>>>'; flush(IO_STDOUT)

!-----------------------------------------------------------------------------
! read numerical parametes and do sanity checks
  num_mech => num_mesh%get_dict('mechanical', defaultVal=emptyDict)

  num%p_i               = int(num_mesh%get_asInt('p_i',defaultVal=2),pPETSCINT)
  num%BBarStabilization = num_mesh%get_asBool('bbarstabilization',defaultVal=.false.)

  num%itmax           = int(num_mech%get_asInt('N_iter_max',defaultVal=250),pPETSCINT)
  num%eps_struct_atol = num_mech%get_asReal('eps_abs_div(P)', defaultVal=1.0e-10_pREAL)
  num%eps_struct_rtol = num_mech%get_asReal('eps_rel_div(P)', defaultVal=1.0e-4_pREAL)

  if (num%itmax <= 1_pPETSCINT)         call IO_error(301,ext_msg='N_iter_max')
  if (num%eps_struct_rtol <= 0.0_pREAL) call IO_error(301,ext_msg='eps_rel_div(P)')
  if (num%eps_struct_atol <= 0.0_pREAL) call IO_error(301,ext_msg='eps_abs_div(P)')

!--------------------------------------------------------------------------------------------------
! Setup FEM mech mesh
  call DMClone(geomMesh,mechanical_mesh,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetDimension(mechanical_mesh,dimPlex,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMSetFromOptions(mechanical_mesh,err_PETSc)
  CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
  call DMPlexIsSimplex(mechanical_mesh,isSimplex,err_PETSc)
  CHKERRQ(err_PETSc)
  if (.not. isSimplex) num%p_i = num%p_i + 1_pPETSCINT                                              ! adjust for quad/hex (non-simplex)
#endif

!--------------------------------------------------------------------------------------------------
! Setup FEM mech discretization
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>17)
#if (PETSC_VERSION_MINOR>=24)
  if (isSimplex) then
    call PetscDTSimplexQuadrature(dimPlex,num%p_i,PETSCDTSIMPLEXQUAD_DEFAULT, &
                                  mechQuad,err_PETSc)
  else
    call PetscDTGaussTensorQuadrature(dimPlex,dimPlex,num%p_i,-1.0_pREAL,1.0_pREAL, &
                                      mechQuad,err_PETSc)
  end if
#elif (PETSC_VERSION_MINOR==23)
  call PetscDTSimplexQuadrature(dimPlex,num%p_i,PETSCDTSIMPLEXQUAD_DEFAULT, &
                                mechQuad,err_PETSc)
#else
  call PetscDTSimplexQuadrature(dimPlex,num%p_i,-1,mechQuad,err_PETSc)
#endif
  CHKERRQ(err_PETSc)
#else
  qPoints  = FEM_quadrature_points( dimPlex,num%p_i)%p
  qWeights = FEM_quadrature_weights(dimPlex,num%p_i)%p
  nQuadrature = FEM_nQuadrature(    dimPlex,num%p_i)
  qPointsP  => qPoints
  qWeightsP => qWeights
  call PetscQuadratureCreate(PETSC_COMM_SELF,mechQuad,err_PETSc)
  CHKERRQ(err_PETSc)
  nc = dimPlex
  call PetscQuadratureSetData(mechQuad,dimPlex,nc,int(nQuadrature,pPETSCINT), &
                              qPointsP,qWeightsP,err_PETSc)
  CHKERRQ(err_PETSc)
#endif

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>21)
  call PetscQuadratureGetData(mechQuad,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                              nQuadrature,PETSC_NULL_REAL_POINTER,qWeightsP,err_PETSc)
  CHKERRQ(err_PETSc)
  qWeights = qWeightsP
  call PetscQuadratureRestoreData(mechQuad,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                                  PETSC_NULL_INTEGER,PETSC_NULL_REAL_POINTER,qWeightsP, &
                                  err_PETSc)
#else
  call PetscQuadratureGetData(mechQuad,PETSC_NULL_INTEGER(1),PETSC_NULL_INTEGER(1), &
                              nQuadrature,PETSC_NULL_REAL_POINTER,qWeightsP,err_PETSc)
  CHKERRQ(err_PETSc)
  qWeights = qWeightsP
  call PetscQuadratureRestoreData(mechQuad,PETSC_NULL_INTEGER(1),PETSC_NULL_INTEGER(1), &
                                  PETSC_NULL_INTEGER(1),PETSC_NULL_REAL_POINTER,qWeightsP, &
                                  err_PETSc)
#endif
  CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
  if (.not. isSimplex) qWeights = [(qWeights(nc), nc=1,size(qWeights),int(dimPlex))]
#endif
  nc = dimPlex
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
  call PetscFECreateDefault(PETSC_COMM_SELF,dimPlex,nc,isSimplex,prefix, &
#else
  call PetscFECreateDefault(PETSC_COMM_SELF,dimPlex,nc,PETSC_TRUE,prefix, &
#endif
                            num%p_i,mechFE,err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscFESetQuadrature(mechFE,mechQuad,err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscFEGetDimension(mechFE,nBasis,err_PETSc)
  CHKERRQ(err_PETSc)
  nBasis = nBasis/nc
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
  call DMAddField(mechanical_mesh,PETSC_NULL_DMLABEL,PetscObjectCast(mechFE),err_PETSc)
#else
  call DMAddField(mechanical_mesh,PETSC_NULL_DMLABEL,mechFE,err_PETSc)
#endif
  CHKERRQ(err_PETSc)
  call DMCreateDS(mechanical_mesh,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetDS(mechanical_mesh,mechDS,err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscDSGetTotalDimension(mechDS,cellDof,err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscFEDestroy(mechFE,err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscQuadratureDestroy(mechQuad,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! Setup FEM mech boundary conditions
  call DMGetLabel(mechanical_mesh,'Face Sets',BCLabel,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMPlexLabelComplete(mechanical_mesh,BCLabel,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetLocalSection(mechanical_mesh,section,err_PETSc)
  CHKERRQ(err_PETSc)
  allocate(pnumComp(1), source=dimPlex)
  allocate(pnumDof(0:dimPlex), source = 0_pPETSCINT)
  do topologDim = 0, dimPlex
    call DMPlexGetDepthStratum(mechanical_mesh,topologDim,cellStart,cellEnd,err_PETSc)
    CHKERRQ(err_PETSc)
    call PetscSectionGetDof(section,cellStart,pnumDof(topologDim),err_PETSc)
    CHKERRQ(err_PETSc)
  end do
  numActiveBC = sum([(count(mechBC(boundary)%active), boundary = 1, size(mechBC))])                ! Number of active DOF in BC
  allocate(pbcField(numActiveBC), source = 0_pPETSCINT)
  allocate(pbcComps(numActiveBC))
  allocate(pbcPoints(numActiveBC))
  numActiveBC = 0_pPETSCINT
  do boundary = 1_pPETSCINT, mesh_Nboundaries; do component = 1_pPETSCINT, dimPlex
    if (mechBC(boundary)%active(component)) then
      numActiveBC = numActiveBC + 1_pPETSCINT
      setLabel = mesh_BCTypeLabel(mesh_boundariesIdx(boundary))
      call ISCreateGeneral(PETSC_COMM_WORLD,1_pPETSCINT,[component-1_pPETSCINT],PETSC_COPY_VALUES, &
                           pbcComps(numActiveBC),err_PETSc)
      CHKERRQ(err_PETSc)
      call DMGetStratumSize(mechanical_mesh,setLabel,mesh_boundariesIS(boundary),bcSize,err_PETSc)
      CHKERRQ(err_PETSc)
      if (bcSize > 0) then
        call DMGetStratumIS(mechanical_mesh,setLabel,mesh_boundariesIS(boundary),bcPoint,err_PETSc)
        CHKERRQ(err_PETSc)
        call ISGetIndices(bcPoint,pBcPoint,err_PETSc)
        CHKERRQ(err_PETSc)
        call ISCreateGeneral(PETSC_COMM_WORLD,bcSize,pBcPoint,PETSC_COPY_VALUES,pbcPoints(numActiveBC),err_PETSc)
        CHKERRQ(err_PETSc)
        call ISRestoreIndices(bcPoint,pBcPoint,err_PETSc)
        CHKERRQ(err_PETSc)
        call ISDestroy(bcPoint,err_PETSc)
        CHKERRQ(err_PETSc)
      else
        call ISCreateGeneral(PETSC_COMM_WORLD,0_pPETSCINT,[0_pPETSCINT],PETSC_COPY_VALUES,pbcPoints(numActiveBC),err_PETSc)
        CHKERRQ(err_PETSc)
      end if
    end if
  end do; end do
  call DMPlexCreateSection(mechanical_mesh,PETSC_NULL_DMLABEL_ARRAY,pNumComp,pNumDof, &
                           numActiveBC,pBcField,pBcComps,pBcPoints,PETSC_NULL_IS,section,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMSetLocalSection(mechanical_mesh,section,err_PETSc)
  CHKERRQ(err_PETSc)
  do boundary = 1_pPETSCINT, numActiveBC
    call ISDestroy(pbcPoints(boundary),err_PETSc)
    CHKERRQ(err_PETSc)
  end do

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
  call SNESCreate(PETSC_COMM_WORLD,mechanical_snes,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetOptionsPrefix(mechanical_snes,'mechanical_',err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetDM(mechanical_snes,mechanical_mesh,err_PETSc)                                         ! set the mesh for non-linear solver
  CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(mechanical_mesh,solution, err_PETSc)                                    ! locally owned displacement Dofs
  CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(mechanical_mesh,solution_rate, err_PETSc)                               ! locally owned velocity Dofs to guess solution at next load step
  CHKERRQ(err_PETSc)
  call DMCreateLocalVector (mechanical_mesh,solution_local,err_PETSc)                               ! locally owned velocity Dofs to guess solution at next load step
  CHKERRQ(err_PETSc)
  call DMSNESSetFunctionLocal(mechanical_mesh,FEM_mechanical_formResidual,PETSC_NULL_VEC,err_PETSc) ! function to evaluate residual forces
  CHKERRQ(err_PETSc)
  call DMSNESSetJacobianLocal(mechanical_mesh,FEM_mechanical_formJacobian,PETSC_NULL_VEC,err_PETSc) ! function to evaluate stiffness matrix
  CHKERRQ(err_PETSc)
  call SNESSetMaxLinearSolveFailures(mechanical_snes, huge(1_pPETSCINT), err_PETSc)                 ! ignore linear solve failures
  CHKERRQ(err_PETSc)
  call SNESSetConvergenceTest(mechanical_snes,FEM_mechanical_converged,PETSC_NULL_VEC,PETSC_NULL_FUNCTION,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetTolerances(mechanical_snes,1.0_pREAL,0.0_pREAL,0.0_pREAL,num%itmax,num%itmax,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetFromOptions(mechanical_snes,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! init fields
  call VecSet(solution     ,0.0_pREAL,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecSet(solution_rate,0.0_pREAL,err_PETSc)
  CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
  call PetscDSGetDiscretization(mechDS,0_pPETSCINT,obj,err_PETSc)
  PetscObjectSpecificCast(mechFE,obj)
#else
  call PetscDSGetDiscretization(mechDS,0_pPETSCINT,mechFE,err_PETSc)
#endif
  CHKERRQ(err_PETSc)
  call PetscFEGetDualSpace(mechFE,mechDualSpace,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMPlexGetHeightStratum(mechanical_mesh,0_pPETSCINT,cellStart,cellEnd,err_PETSc)
  CHKERRQ(err_PETSc)

  nCoords = size(x_n(1:dimPlex,:),kind=pPETSCINT)
  nodeCoords = pack(x_n(1:dimPlex,:), .true.)
  idx = [(nc, nc = 0_pPETSCINT, nCoords - 1_pPETSCINT)]
  call VecSetValuesLocal(solution_local, nCoords, idx, nodeCoords, INSERT_VALUES, err_PETSc)        ! initial node coordinates (undeformed)

  call utilities_constitutiveResponse(status,0.0_pREAL,devNull,.true.)

end subroutine FEM_mechanical_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the FEM load step
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function FEM_mechanical_solution( &
             incInfoIn,Delta_t,Delta_t_prev,mechBC)

!--------------------------------------------------------------------------------------------------
! input data for solution
  real(pREAL), intent(in) :: &
    Delta_t, &                                                                                      !< increment in time for current solution
    Delta_t_prev                                                                                    !< increment in time of last increment
  type(tMechBC), dimension(:),intent(in) :: &
    mechBC
  character(len=*), intent(in) :: &
    incInfoIn

  PetscErrorCode :: err_PETSc
  SNESConvergedReason :: reason

  incInfo = incInfoIn
  FEM_mechanical_solution%converged = .false.
!--------------------------------------------------------------------------------------------------
! set module wide availabe data
  params%Delta_t = Delta_t
  params%mechBC = mechBC

  call SNESSolve(mechanical_snes,PETSC_NULL_VEC,solution,err_PETSc)                                 ! solve mechanical_snes based on solution guess (result in solution)
  CHKERRQ(err_PETSc)
  call SNESGetConvergedReason(mechanical_snes,reason,err_PETSc)                                     ! solution converged?
  CHKERRQ(err_PETSc)

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
  if (reason%v <= SNES_CONVERGED_ITERATING%v) then
#else
  if (reason < 1) then                                                                              ! 0: still iterating (will not occur), negative -> convergence error
#endif
    FEM_mechanical_solution%converged = .false.
    FEM_mechanical_solution%iterationsNeeded = num%itmax
  else                                                                                              ! >= 1 proper convergence (or broken)
    FEM_mechanical_solution%converged = .true.
    call SNESGetIterationNumber(mechanical_snes,FEM_mechanical_solution%iterationsNeeded,err_PETSc)
    CHKERRQ(err_PETSc)
  end if

  print'(/,1x,a)', '==========================================================================='
  flush(IO_STDOUT)

end function FEM_mechanical_solution


!--------------------------------------------------------------------------------------------------
!> @brief Form the FEM residual vector.
!--------------------------------------------------------------------------------------------------
subroutine FEM_mechanical_formResidual(dm_local,xx_local,f_local,dummy,err_PETSc)

  DM                                 :: dm_local
  PetscObject,intent(in)             :: dummy
  PetscErrorCode                     :: err_PETSc
  integer(MPI_INTEGER_KIND)          :: err_MPI
  Vec                                :: f_local, xx_local

  PetscDS                            :: prob
  Vec                                :: x_local
  PetscSection                       :: section
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
  PetscQuadrature                    :: quadrature
#endif
  real(pREAL), dimension(:), pointer      :: x_scal, pf_scal
  real(pREAL), dimension(cellDof), target :: f_scal
  PetscReal,   dimension(dimPlex,dimPlex) :: invCellJ

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
  PetscReal, dimension(:), pointer :: pCellJ, pInvCellJ, pDetJ
#else
  PetscReal, dimension(:), pointer :: pV0, pCellJ, pInvCellJ
  PetscReal                        :: detJ
#endif
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
  PetscTabulation, pointer :: tab(:)
#else
  PetscReal, dimension(:), pointer :: basisFieldDer, &
                                      dev_null
#endif
  PetscInt  :: cellStart, cellEnd, cell, &
               qPt, basis, comp, cidx, &
               numFields, m,i
  PetscReal :: detFAvg
  PetscReal, dimension(dimPlex*dimPlex,cellDof) :: BMat
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
  PetscBool :: isSimplex
#endif

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
  allocate(pCellJ(nQuadrature*dimPlex**2))
  allocate(pInvCellJ(nQuadrature*dimPlex**2))
  allocate(pDetJ(nQuadrature))
#else
  allocate(pV0(dimPlex))
  allocate(pCellJ(dimPlex**2))
  allocate(pInvCellJ(dimPlex**2))
#endif
  allocate(x_scal(cellDof))

  call DMGetLocalSection(dm_local,section,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMPlexGetHeightStratum(dm_local,0_pPETSCINT,cellStart,cellEnd,err_PETSc)
  CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
  call DMPlexIsSimplex(dm_local,isSimplex,err_PETSc)
  CHKERRQ(err_PETSc)
#endif

  call DMGetLocalVector(dm_local,x_local,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecWAXPY(x_local,1.0_pREAL,xx_local,solution_local,err_PETSc)
  CHKERRQ(err_PETSc)

  call utilities_projectBCValues(dm_local,x_local,section,params%mechBC,params%Delta_t,dimPlex)

!--------------------------------------------------------------------------------------------------
! evaluate field derivatives
  call DMGetDS(dm_local,prob,err_PETSc)
  CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
  call PetscDSGetTabulation(prob,tab,err_PETSc)
#else
  call PetscDSGetTabulation(prob,0_pPETSCINT,dev_null,basisFieldDer,err_PETSc)
#endif
  CHKERRQ(err_PETSc)

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
  if (isSimplex) then
    call PetscDTSimplexQuadrature(dimPlex,num%p_i,PETSCDTSIMPLEXQUAD_DEFAULT, &
                                  quadrature,err_PETSc)
  else
    call PetscDTGaussTensorQuadrature(dimPlex,dimPlex,num%p_i,-1.0_pREAL,1.0_pREAL, &
                                      quadrature,err_PETSc)
  end if
#endif
  CHKERRQ(err_PETSc)

  do cell = cellStart, cellEnd-1_pPETSCINT                                                          !< loop over all elements
    call PetscSectionGetNumFields(section,numFields,err_PETSc)
    CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
    call DMPlexVecGetClosure(dm_local,section,x_local,cell,PETSC_NULL_INTEGER,x_scal,err_PETSc)     ! get Dofs belonging to element
#else
    call DMPlexVecGetClosure(dm_local,section,x_local,cell,x_scal,err_PETSc)                        ! get Dofs belonging to element
#endif
    CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
    call DMPlexComputeCellGeometryFEM(dm_local,cell,quadrature,PETSC_NULL_REAL_ARRAY,pCellJ, &
                                      pInvCellJ,pDetJ,err_PETSc)
#else
    call DMPlexComputeCellGeometryAffineFEM(dm_local,cell,pV0,pCellJ,pInvCellJ,detJ,err_PETSc)
#endif
    CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<24)
    invCellJ = reshape(pInvCellJ, shape=[dimPlex,dimPlex])
#endif
    do qPt = 0_pPETSCINT, nQuadrature-1_pPETSCINT
      m = cell*nQuadrature + qPt+1_pPETSCINT
      BMat = 0.0_pREAL
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
      invCellJ = reshape(pInvCellJ(qPt*dimPlex**2+1_pPETSCINT:(qPt+1_pPETSCINT)*dimPlex**2), &
                         shape=[dimPlex,dimPlex])
#endif
      do basis = 0_pPETSCINT, nBasis-1_pPETSCINT
        do comp = 0_pPETSCINT, dimPlex-1_pPETSCINT
          cidx = basis*dimPlex+comp
          i = ((qPt*nBasis + basis)*dimPlex + comp)*dimPlex+comp
          BMat(comp*dimPlex+1_pPETSCINT:(comp+1_pPETSCINT)*dimPlex,basis*dimPlex+comp+1_pPETSCINT) = &
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
            matmul(invCellJ,tab(1)%ptr%T(2)%ptr(i*dimPlex+1_pPETSCINT:(i+1_pPETSCINT)*dimPlex))
#else
            matmul(invCellJ,basisFieldDer(i*dimPlex+1_pPETSCINT:(i+1_pPETSCINT)*dimPlex))
#endif

        end do
      end do
      homogenization_F(1:dimPlex,1:dimPlex,m) = reshape(matmul(BMat,x_scal),shape=[dimPlex,dimPlex], order=[2,1])
    end do
    if (num%BBarStabilization) then
      detFAvg = math_det33(sum(homogenization_F(1:3,1:3,cell*nQuadrature+1:(cell+1)*nQuadrature),dim=3)/real(nQuadrature,pREAL))
      do qPt = 0, nQuadrature-1
        m = cell*nQuadrature + qPt+1
        homogenization_F(1:dimPlex,1:dimPlex,m) = homogenization_F(1:dimPlex,1:dimPlex,m) &
                                                * (detFAvg/math_det33(homogenization_F(1:3,1:3,m)))**(1.0_pREAL/real(dimPlex,pREAL))

      end do
    end if
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
    call DMPlexVecRestoreClosure(dm_local,section,x_local,cell,PETSC_NULL_INTEGER,x_scal,err_PETSc)
#else
    call DMPlexVecRestoreClosure(dm_local,section,x_local,cell,x_scal,err_PETSc)
#endif
    CHKERRQ(err_PETSc)
  end do

!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
  call utilities_constitutiveResponse(status,params%Delta_t,P_av,ForwardData)
  call MPI_Allreduce(MPI_IN_PLACE,status,1_MPI_INTEGER_KIND,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  ForwardData = .false.

!--------------------------------------------------------------------------------------------------
! integrating residual
  do cell = cellStart, cellEnd-1_pPETSCINT                                                          ! loop over all elements
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
    call DMPlexVecGetClosure(dm_local,section,x_local,cell,PETSC_NULL_INTEGER,x_scal,err_PETSc)     ! get Dofs belonging to element
#else
    call DMPlexVecGetClosure(dm_local,section,x_local,cell,x_scal,err_PETSc)                        ! get Dofs belonging to element
#endif
    CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
    call DMPlexComputeCellGeometryFEM(dm_local,cell,quadrature,PETSC_NULL_REAL_ARRAY,pCellJ, &
                                      pInvCellJ,pDetJ,err_PETSc)
#else
    call DMPlexComputeCellGeometryAffineFEM(dm_local,cell,pV0,pCellJ,pInvCellJ,detJ,err_PETSc)
#endif
    CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<24)
    invCellJ = reshape(pInvCellJ, shape=[dimPlex,dimPlex])
#endif
    f_scal = 0.0_pREAL
    do qPt = 0_pPETSCINT, nQuadrature-1_pPETSCINT
      m = cell*nQuadrature + qPt+1_pPETSCINT
      BMat = 0.0_pREAL
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
      invCellJ = reshape(pInvCellJ(qPt*dimPlex**2+1_pPETSCINT:(qPt+1_pPETSCINT)*dimPlex**2), &
                         shape=[dimPlex,dimPlex])
#endif
      do basis = 0_pPETSCINT, nBasis-1_pPETSCINT
        do comp = 0_pPETSCINT, dimPlex-1_pPETSCINT
          cidx = basis*dimPlex+comp
          i = ((qPt*nBasis + basis)*dimPlex + comp)*dimPlex+comp
          BMat(comp*dimPlex+1_pPETSCINT:(comp+1_pPETSCINT)*dimPlex,basis*dimPlex+comp+1_pPETSCINT) = &
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
            matmul(invCellJ,tab(1)%ptr%T(2)%ptr(i*dimPlex+1_pPETSCINT:(i+1_pPETSCINT)*dimPlex))
#else
            matmul(invCellJ,basisFieldDer(i*dimPlex+1_pPETSCINT:(i+1_pPETSCINT)*dimPlex))
#endif
        end do
      end do
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
      f_scal = f_scal + pDetJ(qPt+1_pPETSCINT) * qWeights(qPt+1_pPETSCINT) &
#else
      f_scal = f_scal + abs(detJ) * qWeights(qPt+1_pPETSCINT) &
#endif
             * matmul(transpose(BMat), &
                      reshape(transpose(homogenization_P(1:dimPlex,1:dimPlex,m)), &
                              shape=[dimPlex*dimPlex]))
    end do
    pf_scal => f_scal
    call DMPlexVecSetClosure(dm_local,section,f_local,cell,pf_scal,ADD_VALUES,err_PETSc)
    CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
    call DMPlexVecRestoreClosure(dm_local,section,x_local,cell,PETSC_NULL_INTEGER,x_scal,err_PETSc)
#else
    call DMPlexVecRestoreClosure(dm_local,section,x_local,cell,x_scal,err_PETSc)
#endif
    CHKERRQ(err_PETSc)
  end do
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
  call PetscDSRestoreTabulation(prob,tab,err_PETSc)
#else
  call PetscDSRestoreTabulation(prob,0_pPETSCINT,dev_null,basisFieldDer,err_PETSc)
#endif
  CHKERRQ(err_PETSc)
  call DMRestoreLocalVector(dm_local,x_local,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine FEM_mechanical_formResidual


!--------------------------------------------------------------------------------------------------
!> @brief Form the FEM stiffness matrix.
!--------------------------------------------------------------------------------------------------
subroutine FEM_mechanical_formJacobian(dm_local,xx_local,J,Jp,dummy,err_PETSc)

  DM                      :: dm_local
  Mat                     :: J, Jp
  PetscObject, intent(in) :: dummy
  PetscErrorCode          :: err_PETSc

  PetscDS      :: prob
  Vec          :: x_local, xx_local
  PetscSection :: section, gSection
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
  PetscQuadrature :: quadrature
#endif

  PetscReal, dimension(1,         cellDof) :: MatB
  PetscReal, dimension(dimPlex**2,cellDof) :: BMat, BMatAvg, MatA
  PetscReal, dimension(3,3) :: F, FAvg, FInv

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
  real(pREAL), dimension(:), pointer :: pCellJ, pInvCellJ, pDetJ
#else
  real(pREAL), dimension(:), pointer :: pV0, pCellJ, pInvCellJ
  PetscReal :: detJ
#endif
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
  PetscTabulation, pointer :: tab(:)
#else
  PetscReal,   dimension(:), pointer :: basisFieldDer, &
                                        dev_null
#endif
  real(pREAL), dimension(:), pointer :: pK_e, x_scal

  real(pREAL), dimension(cellDOF,cellDOF), target :: K_e
  real(pREAL), dimension(cellDOF,cellDOF)         :: K_eA, K_eB
  real(pREAL), dimension(dimPlex,dimPlex)         :: invCellJ

  PetscInt :: cellStart, cellEnd, cell, &
              qPt, basis, comp, cidx, ce, i
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
  PetscBool :: isSimplex
#endif


#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
  allocate(pCellJ(nQuadrature*dimPlex**2))
  allocate(pInvcellJ(nQuadrature*dimPlex**2))
  allocate(pDetJ(nQuadrature))
#else
  allocate(pV0(dimPlex))
  allocate(pCellJ(dimPlex**2))
  allocate(pInvCellJ(dimPlex**2))
#endif

  call MatSetOption(Jp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatSetOption(Jp,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatZeroEntries(Jp,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetDS(dm_local,prob,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetLocalSection(dm_local,section,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetGlobalSection(dm_local,gSection,err_PETSc)
  CHKERRQ(err_PETSc)

  call DMGetLocalVector(dm_local,x_local,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecWAXPY(x_local,1.0_pREAL,xx_local,solution_local,err_PETSc)
  CHKERRQ(err_PETSc)

  call utilities_projectBCValues(dm_local,x_local,section,params%mechBC,params%Delta_t,dimPlex)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
  call PetscDSGetTabulation(prob,tab,err_PETSc)
#else
  call PetscDSGetTabulation(prob,0_pPETSCINT,dev_null,basisFieldDer,err_PETSc)
#endif
  CHKERRQ(err_PETSc)
  call DMPlexGetHeightStratum(dm_local,0_pPETSCINT,cellStart,cellEnd,err_PETSc)
  CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
  call DMPlexIsSimplex(dm_local,isSimplex,err_PETSc)
  CHKERRQ(err_PETSc)
  if (isSimplex) then
    call PetscDTSimplexQuadrature(dimPlex,num%p_i,PETSCDTSIMPLEXQUAD_DEFAULT, &
                                  quadrature,err_PETSc)
  else
    call PetscDTGaussTensorQuadrature(dimPlex,dimPlex,num%p_i,-1.0_pREAL,1.0_pREAL, &
                                      quadrature,err_PETSc)
  end if
  CHKERRQ(err_PETSc)
#endif
  do cell = cellStart, cellEnd-1                                                                    !< loop over all elements
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
    call DMPlexVecGetClosure(dm_local,section,x_local,cell,PETSC_NULL_INTEGER,x_scal,err_PETSc)     !< get Dofs belonging to el
#else
    call DMPlexVecGetClosure(dm_local,section,x_local,cell,x_scal,err_PETSc)                        !< get Dofs belonging to element
#endif
    CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
    call DMPlexComputeCellGeometryFEM(dm_local,cell,quadrature,PETSC_NULL_REAL_ARRAY,pCellJ, &
                                      pInvCellJ,pDetJ,err_PETSc)
#else
    call DMPlexComputeCellGeometryAffineFEM(dm_local,cell,pV0,pCellJ,pInvCellJ,detJ,err_PETSc)
#endif
    CHKERRQ(err_PETSc)
    K_eA = 0.0_pREAL
    K_eB = 0.0_pREAL
    MatB = 0.0_pREAL
    FAvg = 0.0_pREAL
    BMatAvg = 0.0_pREAL
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<24)
    invCellJ = reshape(pInvCellJ, shape=[dimPlex,dimPlex])
#endif
    do qPt = 0_pPETSCINT, nQuadrature-1_pPETSCINT
      ce = cell*nQuadrature + qPt + 1_pPETSCINT
      BMat = 0.0_pREAL
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
      invCellJ = reshape(pInvCellJ(qPt*dimPlex**2+1_pPETSCINT:(qPt+1_pPETSCINT)*dimPlex**2), &
                         shape=[dimPlex,dimPlex])
#endif
      do basis = 0_pPETSCINT, nBasis-1_pPETSCINT
        do comp = 0_pPETSCINT, dimPlex-1_pPETSCINT
          cidx = basis*dimPlex+comp
          i = ((qPt*nBasis + basis)*dimPlex + comp)*dimPlex+comp
          BMat(comp*dimPlex+1_pPETSCINT:(comp+1_pPETSCINT)*dimPlex,basis*dimPlex+comp+1_pPETSCINT) = &
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
            matmul(invCellJ,tab(1)%ptr%T(2)%ptr(i*dimPlex+1_pPETSCINT:(i+1_pPETSCINT)*dimPlex))
#else
            matmul(invCellJ,basisFieldDer(i*dimPlex+1_pPETSCINT:(i+1_pPETSCINT)*dimPlex))
#endif
        end do
      end do
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=24)
      MatA = qWeights(qPt+1_pPETSCINT) * pDetJ(qPt+1_pPETSCINT) &
#else
      MatA = qWeights(qPt+1_pPETSCINT) * abs(detJ) &
#endif
           * matmul(reshape(reshape(homogenization_dPdF(1:dimPlex,1:dimPlex,1:dimPlex,1:dimPlex,ce), &
                                    shape=[dimPlex,dimPlex,dimPlex,dimPlex], order=[2,1,4,3]), &
                            shape=[dimPlex*dimPlex,dimPlex*dimPlex]),BMat)

      if (num%BBarStabilization) then
        F(1:dimPlex,1:dimPlex) = reshape(matmul(BMat,x_scal),shape=[dimPlex,dimPlex])
        FInv = math_inv33(F)
        K_eA = K_eA + matmul(transpose(BMat),MatA)*math_det33(FInv)**(1.0_pREAL/real(dimPlex,pREAL))
        K_eB = K_eB - &
               matmul(transpose(matmul(reshape(homogenization_F(1:dimPlex,1:dimPlex,ce),shape=[dimPlex**2,1_pPETSCINT]), &
                                       matmul(reshape(FInv(1:dimPlex,1:dimPlex), &
                                                      shape=[1_pPETSCINT,dimPlex**2],order=[2,1]),BMat))),MatA)
        MatB = MatB &
             + matmul(reshape(homogenization_F(1:dimPlex,1:dimPlex,ce),shape=[1_pPETSCINT,dimPlex**2]),MatA)
        FAvg = FAvg + F
        BMatAvg = BMatAvg + BMat
      else
        K_eA = K_eA + matmul(transpose(BMat),MatA)
      end if
    end do
    if (num%BBarStabilization) then
      FInv = math_inv33(FAvg)
      K_e = K_eA*math_det33(FAvg/real(nQuadrature,pREAL))**(1.0_pREAL/real(dimPlex,pREAL)) + &
            (matmul(matmul(transpose(BMatAvg), &
                           reshape(FInv(1:dimPlex,1:dimPlex),shape=[dimPlex**2,1_pPETSCINT],order=[2,1])),MatB) + &
             K_eB)/real(dimPlex,pREAL)
    else
      K_e = K_eA
    end if
    K_e = (K_e + eps*math_eye(int(cellDof)))
    pK_e(1:cellDOF**2) => K_e
    call DMPlexMatSetClosure(dm_local,section,gSection,Jp,cell,pK_e,ADD_VALUES,err_PETSc)
    CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
    call DMPlexVecRestoreClosure(dm_local,section,x_local,cell,PETSC_NULL_INTEGER,x_scal,err_PETSc)
#else
    call DMPlexVecRestoreClosure(dm_local,section,x_local,cell,x_scal,err_PETSc)
#endif
    CHKERRQ(err_PETSc)
  end do
  call MatAssemblyBegin(Jp,MAT_FINAL_ASSEMBLY,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatAssemblyEnd(Jp,MAT_FINAL_ASSEMBLY,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMRestoreLocalVector(dm_local,x_local,err_PETSc)
  CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
  call PetscDSRestoreTabulation(prob,tab,err_PETSc)
#else
  call PetscDSRestoreTabulation(prob,0_pPETSCINT,dev_null,basisFieldDer,err_PETSc)
#endif
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! apply boundary conditions
  call DMPlexCreateRigidBody(dm_local,0_pPETSCINT,matnull,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatSetNullSpace(Jp,matnull,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatSetNearNullSpace(Jp,matnull,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatNullSpaceDestroy(matnull,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine FEM_mechanical_formJacobian


!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine FEM_mechanical_forward(guess,Delta_t,Delta_t_prev,mechBC)

  type(tMechBC),  dimension(:), intent(in) :: &
    mechBC
  real(pREAL),    intent(in) :: &
    Delta_t_prev, &
    Delta_t
  logical,        intent(in) :: &
    guess

  DM             :: dm_local
  Vec            :: x_local
  PetscSection   :: section
  PetscErrorCode :: err_PETSc

!--------------------------------------------------------------------------------------------------
! forward last inc
  if (guess .and. .not. cutBack) then
    ForwardData = .True.
    call SNESGetDM(mechanical_snes,dm_local,err_PETSc)                                              ! retrieve mesh info from mechanical_snes into dm_local
    CHKERRQ(err_PETSc)
    call DMGetLocalSection(dm_local,section,err_PETSc)
    CHKERRQ(err_PETSc)
    call DMGetLocalVector(dm_local,x_local,err_PETSc)
    CHKERRQ(err_PETSc)
    call VecSet(x_local,0.0_pREAL,err_PETSc)
    CHKERRQ(err_PETSc)
    call DMGlobalToLocalBegin(dm_local,solution,INSERT_VALUES,x_local,err_PETSc)                    ! retrieve my partition of global solution vector
    CHKERRQ(err_PETSc)
    call DMGlobalToLocalEnd(dm_local,solution,INSERT_VALUES,x_local,err_PETSc)
    CHKERRQ(err_PETSc)
    call VecAXPY(solution_local,1.0_pREAL,x_local,err_PETSc)
    CHKERRQ(err_PETSc)

    call utilities_projectBCValues(dm_local,solution_local,section,mechBC,Delta_t_prev,dimPlex)

    call DMRestoreLocalVector(dm_local,x_local,err_PETSc)
    CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! update rate and forward last inc
    call VecCopy(solution,solution_rate,err_PETSc)
    CHKERRQ(err_PETSc)
    call VecScale(solution_rate,Delta_t_prev**(-1),err_PETSc)
    CHKERRQ(err_PETSc)
  end if
  call VecCopy(solution_rate,solution,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecScale(solution,Delta_t,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine FEM_mechanical_forward


!--------------------------------------------------------------------------------------------------
!> @brief reporting
!--------------------------------------------------------------------------------------------------
subroutine FEM_mechanical_converged(snes_local,PETScIter,xnorm,snorm,fnorm,reason,dummy,err_PETSc)

  SNES :: snes_local
  PetscInt :: PETScIter
  PetscReal :: xnorm,snorm,fnorm,divTol
  SNESConvergedReason :: reason
  PetscObject :: dummy
  PetscErrorCode :: err_PETSc

!--------------------------------------------------------------------------------------------------
! report
  divTol = max(maxval(abs(P_av(1:dimPlex,1:dimPlex)))*num%eps_struct_rtol,num%eps_struct_atol)
  call SNESConvergedDefault(snes_local,PETScIter,xnorm,snorm,fnorm/divTol,reason,dummy,err_PETSc)
  CHKERRQ(err_PETSc)
  if (status /= STATUS_OK) reason = SNES_DIVERGED_FUNCTION_DOMAIN
  print'(/,1x,a,a,i0,a,f0.3)', trim(incInfo), &
                  ' @ Iteration ',PETScIter,' mechanical residual norm = ',fnorm/divTol
  print'(/,1x,a,/,2(3(2x,f12.4,1x)/),3(2x,f12.4,1x))', &
    'Piola--Kirchhoff stress / MPa =',transpose(P_av)*1.e-6_pREAL
  flush(IO_STDOUT)

end subroutine FEM_mechanical_converged


!--------------------------------------------------------------------------------------------------
!> @brief Calculate current coordinates (both nodal and ip coordinates)
!--------------------------------------------------------------------------------------------------
subroutine FEM_mechanical_updateCoords()

  PetscReal, pointer, dimension(:,:) :: &
    nodeCoords                                                                                      !< nodal coordinates (3,nNodes)
  real(pREAL), pointer, dimension(:,:,:) :: &
    ipCoords                                                                                        !< ip coordinates (3,nQuadrature,mesh_nElems)

  integer :: &
    qPt, &
    comp, &
    qOffset, &
    nOffset

  DM  :: dm_local
  Vec :: x_local
  PetscErrorCode :: err_PETSc
  PetscInt :: nNodes, cellStart, cellEnd, q, c, n
  PetscSection :: section
  PetscDS :: mechDS
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
  PetscTabulation, pointer :: tab(:)
#else
  PetscReal,   dimension(:), pointer :: basisField, &
                                        dev_null
#endif
  PetscReal,   dimension(:), pointer :: nodeCoordsDM                                                ! nodal coordinates read from DM (dimPlex*nNodes)
  real(pREAL), dimension(:), pointer :: x_scal

  call SNESGetDM(mechanical_snes,dm_local,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetDS(dm_local,mechDS,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetLocalSection(dm_local,section,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetLocalVector(dm_local,x_local,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetDimension(dm_local,dimPlex,err_PETSc)
  CHKERRQ(err_PETSc)

  ! write nodes displacements
  call VecGetArray(x_local,nodeCoordsDM,err_PETSc)
  CHKERRQ(err_PETSc)
  nNodes = size(nodeCoordsDM,kind=pPETSCINT)/dimPlex
  allocate(nodeCoords(3,nNodes),source=0.0_pREAL)
  nodeCoords(1:dimPlex,:) = reshape(nodeCoordsDM, [dimPlex, nNodes])
  call discretization_setNodeCoords(nodeCoords)
  call VecRestoreArray(x_local,nodeCoordsDM,err_PETSc)
  CHKERRQ(err_PETSc)

  ! write ip displacements
  call DMPlexGetHeightStratum(dm_local,0_pPETSCINT,cellStart,cellEnd,err_PETSc)
  CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
  call PetscDSGetTabulation(mechDS,tab,err_PETSc)
#else
  call PetscDSGetTabulation(mechDS,0_pPETSCINT,basisField,dev_null,err_PETSc)
#endif
  CHKERRQ(err_PETSc)
  allocate(ipCoords(3,nQuadrature,mesh_nElems),source=0.0_pREAL)
  do c = cellStart, cellEnd - 1_pPETSCINT
    qOffset=0
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
    call DMPlexVecGetClosure(dm_local,section,x_local,c,PETSC_NULL_INTEGER,x_scal,err_PETSc)        ! get nodal coordinates of each element
#else
    call DMPlexVecGetClosure(dm_local,section,x_local,c,x_scal,err_PETSc)                           ! get nodal coordinates of each element
#endif
    CHKERRQ(err_PETSc)
    do qPt=0,nQuadrature-1
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
      qOffset = qPt * (size(tab(1)%ptr%T(1)%ptr)/nQuadrature)
#else
      qOffset = qPt * (size(basisField)/nQuadrature)
#endif
      do comp=0,dimPlex-1                                                                           ! loop over components
        nOffset=0
        q = comp
        do n=0,nBasis-1
          ipCoords(comp+1,qPt+1,c+1)=ipCoords(comp+1,qPt+1,c+1)+&
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
                                     sum(tab(1)%ptr%T(1)%ptr(qOffset+(q*dimPlex)+1:qOffset+(q*dimPlex)+dimPlex)*&
#else
                                     sum(basisField(qOffset+(q*dimPlex)+1:qOffset+(q*dimPlex)+dimPlex)*&
#endif
                                     x_scal(nOffset+1:nOffset+dimPlex))
          q = q+dimPlex
          nOffset = nOffset+dimPlex
        end do
      end do
    end do
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
    call DMPlexVecRestoreClosure(dm_local,section,x_local,c,PETSC_NULL_INTEGER,x_scal,err_PETSc)
#else
    call DMPlexVecRestoreClosure(dm_local,section,x_local,c,x_scal,err_PETSc)
#endif
    CHKERRQ(err_PETSc)
  end do
  call discretization_setIPcoords(reshape(ipCoords,[3,int(mesh_nElems*nQuadrature)]))
  call DMRestoreLocalVector(dm_local,x_local,err_PETSc)
  CHKERRQ(err_PETSc)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>22)
  call PetscDSRestoreTabulation(mechDS,tab,err_PETSc)
#else
  call PetscDSRestoreTabulation(mechDS,0_pPETSCINT,basisField,dev_null,err_PETSc)
#endif
  CHKERRQ(err_PETSc)

end subroutine FEM_mechanical_updateCoords

end module mesh_mechanical_FEM
