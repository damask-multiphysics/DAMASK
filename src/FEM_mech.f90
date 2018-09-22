!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief FEM PETSc solver
!--------------------------------------------------------------------------------------------------
module FEM_mech
#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petsc.h>

 use PETScsnes
 use PETScDM
 use PETScDMplex
 use PETScDT
 use prec, only: & 
   pInt, &
   pReal
 use math, only: &
   math_I3
 use FEM_utilities, only: &
   tSolutionState, &
   tFieldBC, &
   tComponentBC
 use mesh, only: &
   mesh_Nboundaries, &
   mesh_boundaries

 implicit none
 private
   
!--------------------------------------------------------------------------------------------------
! derived types
 type tSolutionParams 
   type(tFieldBC)  :: fieldBC
   real(pReal)     :: timeinc
   real(pReal)     :: timeincOld
 end type tSolutionParams
 
 type(tSolutionParams),   private :: params

!--------------------------------------------------------------------------------------------------
! PETSc data
 SNES,                           private :: mech_snes
 Vec,                            private :: solution, solution_rate, solution_local
 PetscInt,                       private :: dimPlex, cellDof, nQuadrature, nBasis
 PetscReal, allocatable, target, private :: qPoints(:), qWeights(:)
 MatNullSpace,                   private :: matnull

!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
 character(len=1024), private :: incInfo   
 real(pReal), private, dimension(3,3) :: &
   P_av = 0.0_pReal
 logical, private :: ForwardData
 real(pReal), parameter, private :: eps = 1.0e-18_pReal

 public :: &
   FEM_mech_init, &
   FEM_mech_solution ,&
   FEM_mech_forward, &
   FEM_mech_destroy
contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
subroutine FEM_mech_init(fieldBC)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 use IO, only: &
   IO_timeStamp, &
   IO_error
 use DAMASK_interface, only: &
   getSolverJobName
 use mesh, only: &
   geomMesh
 use numerics, only: &
   itmax, &
   integrationOrder
 use FEM_Zoo, only: &
   FEM_Zoo_nQuadrature, &
   FEM_Zoo_QuadraturePoints, &
   FEM_Zoo_QuadratureWeights  

 implicit none
 type(tFieldBC),             intent(in) :: fieldBC
 DM                                     :: mech_mesh
 PetscFE                                :: mechFE
 PetscQuadrature                        :: mechQuad, functional
 PetscDS                                :: mechDS
 PetscDualSpace                         :: mechDualSpace
 DMLabel                                :: BCLabel
 PetscInt,          allocatable, target :: numComp(:), numDoF(:), bcField(:)
 PetscInt,                      pointer :: pNumComp(:), pNumDof(:), pBcField(:), pBcPoint(:)
 PetscInt                               :: numBC, bcSize, nc
 IS                                     :: bcPoint
 IS,                allocatable, target :: bcComps(:), bcPoints(:)
 IS,                            pointer :: pBcComps(:), pBcPoints(:)
 PetscSection                           :: section
 PetscInt                               :: field, faceSet, topologDim, nNodalPoints
 PetscReal,      dimension(:)     ,  pointer :: qPointsP, qWeightsP, &
                                           nodalPointsP, nodalWeightsP
 PetscReal,         allocatable, target :: nodalPoints(:), nodalWeights(:)
 PetscScalar,                   pointer :: px_scal(:)
 PetscScalar,       allocatable, target ::  x_scal(:)
 PetscReal                              :: detJ
 PetscReal,         allocatable, target :: v0(:), cellJ(:), invcellJ(:), cellJMat(:,:)
 PetscReal,                     pointer :: pV0(:), pCellJ(:), pInvcellJ(:)
 PetscInt                               :: cellStart, cellEnd, cell, basis 
 character(len=7)                       :: prefix = 'mechFE_'
 PetscErrorCode                         :: ierr

 write(6,'(/,a)') ' <<<+-  FEM_mech init  -+>>>'
 write(6,'(a15,a)')   ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

!--------------------------------------------------------------------------------------------------
! Setup FEM mech mesh
 call DMClone(geomMesh,mech_mesh,ierr); CHKERRQ(ierr)
 call DMGetDimension(mech_mesh,dimPlex,ierr); CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! Setup FEM mech discretization
 qPoints = FEM_Zoo_QuadraturePoints(dimPlex,integrationOrder)%p
 qWeights = FEM_Zoo_QuadratureWeights(dimPlex,integrationOrder)%p
 nQuadrature = FEM_Zoo_nQuadrature(dimPlex,integrationOrder)
 qPointsP => qPoints
 qWeightsP => qWeights
 call PetscQuadratureCreate(PETSC_COMM_SELF,mechQuad,ierr); CHKERRQ(ierr)
 CHKERRQ(ierr)
 ! what is the number of components, nc?
 nc = dimPlex
 call PetscQuadratureSetData(mechQuad,dimPlex,nc,nQuadrature,qPointsP,qWeightsP,ierr)
 CHKERRQ(ierr)
! call PetscFECreateDefault(mech_mesh,dimPlex,nc,PETSC_TRUE,prefix, &
!                           integrationOrder,mechFE,ierr); CHKERRQ(ierr)
 call PetscFECreateDefault(PETSC_COMM_SELF,dimPlex,nc,PETSC_TRUE,prefix, &
                           integrationOrder,mechFE,ierr); CHKERRQ(ierr)
!   Polar decomposition failed in run time
 call PetscFESetQuadrature(mechFE,mechQuad,ierr); CHKERRQ(ierr)
 call PetscFEGetDimension(mechFE,nBasis,ierr); CHKERRQ(ierr)
 nBasis = nBasis/nc
 call DMGetDS(mech_mesh,mechDS,ierr); CHKERRQ(ierr)
 call PetscDSAddDiscretization(mechDS,mechFE,ierr); CHKERRQ(ierr)
 call PetscDSGetTotalDimension(mechDS,cellDof,ierr); CHKERRQ(ierr)
 call PetscFEDestroy(mechFE,ierr); CHKERRQ(ierr)
 call PetscQuadratureDestroy(mechQuad,ierr); CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! Setup FEM mech boundary conditions
 write(6,*) 'starting to set up boundary conditions';flush(6)
 call DMGetLabel(mech_mesh,'Face Sets',BCLabel,ierr); CHKERRQ(ierr)
 call DMPlexLabelComplete(mech_mesh,BCLabel,ierr); CHKERRQ(ierr)
 call DMGetSection(mech_mesh,section,ierr); CHKERRQ(ierr)
 allocate(numComp(1), source=dimPlex); pNumComp => numComp
 allocate(numDof(dimPlex+1), source = 0); pNumDof  => numDof
 do topologDim = 0, dimPlex
   call DMPlexGetDepthStratum(mech_mesh,topologDim,cellStart,cellEnd,ierr)
   CHKERRQ(ierr)
   call PetscSectionGetDof(section,cellStart,numDof(topologDim+1),ierr)
   CHKERRQ(ierr)
 enddo
 numBC = 0
 do field = 1, dimPlex; do faceSet = 1, mesh_Nboundaries
   if (fieldBC%componentBC(field)%Mask(faceSet)) numBC = numBC + 1
 enddo; enddo
 allocate(bcField(numBC), source=0); pBcField => bcField 
 allocate(bcComps(numBC)); pBcComps => bcComps
 allocate(bcPoints(numBC)); pBcPoints => bcPoints
 numBC = 0
 do field = 1, dimPlex; do faceSet = 1, mesh_Nboundaries
   if (fieldBC%componentBC(field)%Mask(faceSet)) then
     numBC = numBC + 1
     write(6,*) 'adding boundary condition', numBC
     call ISCreateGeneral(PETSC_COMM_WORLD,1,[field-1],PETSC_COPY_VALUES,bcComps(numBC),ierr)
     CHKERRQ(ierr)
     call DMGetStratumSize(mech_mesh,'Face Sets',mesh_boundaries(faceSet),bcSize,ierr)
     CHKERRQ(ierr)
     if (bcSize > 0) then
       call DMGetStratumIS(mech_mesh,'Face Sets',mesh_boundaries(faceSet),bcPoint,ierr)
       CHKERRQ(ierr)
       call ISGetIndicesF90(bcPoint,pBcPoint,ierr); CHKERRQ(ierr)
       call ISCreateGeneral(PETSC_COMM_WORLD,bcSize,pBcPoint,PETSC_COPY_VALUES,bcPoints(numBC),ierr)
       CHKERRQ(ierr)
       call ISRestoreIndicesF90(bcPoint,pBcPoint,ierr); CHKERRQ(ierr)
       call ISDestroy(bcPoint,ierr); CHKERRQ(ierr)
     else
       call ISCreateGeneral(PETSC_COMM_WORLD,0,[0],PETSC_COPY_VALUES,bcPoints(numBC),ierr)
       CHKERRQ(ierr)
     endif  
   endif  
 enddo; enddo
 call DMPlexCreateSection(mech_mesh,dimPlex,1,pNumComp,pNumDof, &
                          numBC,pBcField,pBcComps,pBcPoints,PETSC_NULL_IS, &
                          section,ierr)
 CHKERRQ(ierr)
 call DMSetSection(mech_mesh,section,ierr); CHKERRQ(ierr)
 do faceSet = 1, numBC
   call ISDestroy(bcPoints(faceSet),ierr); CHKERRQ(ierr)
 enddo
 
!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
 call SNESCreate(PETSC_COMM_WORLD,mech_snes,ierr);CHKERRQ(ierr) 
 call SNESSetOptionsPrefix(mech_snes,'mech_',ierr);CHKERRQ(ierr) 
 call SNESSetDM(mech_snes,mech_mesh,ierr); CHKERRQ(ierr)                                             !< set the mesh for non-linear solver
 call DMCreateGlobalVector(mech_mesh,solution        ,ierr); CHKERRQ(ierr)                           !< locally owned displacement Dofs
 call DMCreateGlobalVector(mech_mesh,solution_rate   ,ierr); CHKERRQ(ierr)                           !< locally owned velocity Dofs to guess solution at next load step
 call DMCreateLocalVector (mech_mesh,solution_local  ,ierr); CHKERRQ(ierr)                           !< locally owned velocity Dofs to guess solution at next load step
 call DMSNESSetFunctionLocal(mech_mesh,FEM_mech_formResidual,PETSC_NULL_VEC,ierr)                 !< function to evaluate residual forces
 CHKERRQ(ierr)
 call DMSNESSetJacobianLocal(mech_mesh,FEM_mech_formJacobian,PETSC_NULL_VEC,ierr)                 !< function to evaluate stiffness matrix
 CHKERRQ(ierr)
 call SNESSetMaxLinearSolveFailures(mech_snes, huge(1), ierr); CHKERRQ(ierr)                         !< ignore linear solve failures 
 call SNESSetConvergenceTest(mech_snes,FEM_mech_converged,PETSC_NULL_VEC,PETSC_NULL_FUNCTION,ierr)
 CHKERRQ(ierr)
 call SNESSetTolerances(mech_snes,1.0,0.0,0.0,itmax,itmax,ierr)
 CHKERRQ(ierr)
 call SNESSetFromOptions(mech_snes,ierr); CHKERRQ(ierr)
  
!--------------------------------------------------------------------------------------------------
! init fields
 call VecSet(solution        ,0.0,ierr); CHKERRQ(ierr)
 call VecSet(solution_rate   ,0.0,ierr); CHKERRQ(ierr)
 allocate(x_scal(cellDof))
 allocate(nodalPoints (dimPlex))
 allocate(nodalWeights(1))
 nodalPointsP  => nodalPoints
 nodalWeightsP => nodalWeights
 allocate(v0(dimPlex))
 allocate(cellJ(dimPlex*dimPlex))
 allocate(invcellJ(dimPlex*dimPlex))
 allocate(cellJMat(dimPlex,dimPlex))
 pV0 => v0
 pCellJ => cellJ
 pInvcellJ => invcellJ
 call DMGetSection(mech_mesh,section,ierr); CHKERRQ(ierr)
 call DMGetDS(mech_mesh,mechDS,ierr); CHKERRQ(ierr)
 call PetscDSGetDiscretization(mechDS,0,mechFE,ierr)
 CHKERRQ(ierr)
 call PetscFEGetDualSpace(mechFE,mechDualSpace,ierr); CHKERRQ(ierr)
 call DMPlexGetHeightStratum(mech_mesh,0,cellStart,cellEnd,ierr)
 CHKERRQ(ierr)
 do cell = cellStart, cellEnd-1                                                                     !< loop over all elements 
   x_scal = 0.0
   call  DMPlexComputeCellGeometryAffineFEM(mech_mesh,cell,pV0,pCellJ,pInvcellJ,detJ,ierr) 
   CHKERRQ(ierr)
   cellJMat = reshape(pCellJ,shape=[dimPlex,dimPlex])
   do basis = 0, nBasis-1
     call PetscDualSpaceGetFunctional(mechDualSpace,basis,functional,ierr)
     CHKERRQ(ierr)
     call PetscQuadratureGetData(functional,dimPlex,nc,nNodalPoints,nodalPointsP,nodalWeightsP,ierr)
     CHKERRQ(ierr)
     x_scal(basis*dimPlex+1:(basis+1)*dimPlex) = pV0 + matmul(transpose(cellJMat),nodalPointsP + 1.0)
   enddo
   px_scal => x_scal
   call DMPlexVecSetClosure(mech_mesh,section,solution_local,cell,px_scal,INSERT_ALL_VALUES,ierr)
   CHKERRQ(ierr) 
 enddo 

end subroutine FEM_mech_init
  
!--------------------------------------------------------------------------------------------------
!> @brief solution for the FEM load step
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function FEM_mech_solution( &
             incInfoIn,timeinc,timeinc_old,fieldBC)
 use numerics, only: &
   itmax
 use FEsolving, only: &
   terminallyIll

 implicit none
!--------------------------------------------------------------------------------------------------
! input data for solution
 real(pReal), intent(in) :: &
   timeinc, &                                                                                       !< increment in time for current solution
   timeinc_old                                                                                      !< increment in time of last increment
 type(tFieldBC),      intent(in) :: &
   fieldBC
 character(len=*), intent(in) :: &
   incInfoIn
 
!--------------------------------------------------------------------------------------------------
! 
 PetscErrorCode :: ierr   
 SNESConvergedReason :: reason

 incInfo = incInfoIn
 FEM_mech_solution%converged =.false.
!--------------------------------------------------------------------------------------------------
! set module wide availabe data 
 params%timeinc = timeinc
 params%timeincOld = timeinc_old
 params%fieldBC = fieldBC

 call SNESSolve(mech_snes,PETSC_NULL_VEC,solution,ierr); CHKERRQ(ierr)                           ! solve mech_snes based on solution guess (result in solution)
 call SNESGetConvergedReason(mech_snes,reason,ierr); CHKERRQ(ierr)                                  ! solution converged?
 terminallyIll = .false.

 if (reason < 1) then                                                                               ! 0: still iterating (will not occur), negative -> convergence error
   FEM_mech_solution%converged = .false.
   FEM_mech_solution%iterationsNeeded = itmax
 else                                                                                               ! >= 1 proper convergence (or terminally ill)
   FEM_mech_solution%converged = .true.
   call SNESGetIterationNumber(mech_snes,FEM_mech_solution%iterationsNeeded,ierr)
   CHKERRQ(ierr)
 endif

 write(6,'(/,a)') ' ==========================================================================='
 flush(6) 

end function FEM_mech_solution


!--------------------------------------------------------------------------------------------------
!> @brief forms the FEM residual vector
!--------------------------------------------------------------------------------------------------
subroutine FEM_mech_formResidual(dm_local,xx_local,f_local,dummy,ierr)
 use numerics, only: &
   BBarStabilisation
 use FEM_utilities, only: &
   utilities_projectBCValues, &
   utilities_constitutiveResponse
 use homogenization, only: &
   materialpoint_F, &
   materialpoint_P
 use math, only: &
   math_I3, &
   math_det33, &
   math_inv33  
 use FEsolving, only: &
   terminallyIll

 implicit none
 DM                                 :: dm_local
 PetscDS                            :: prob
 Vec                                :: x_local, f_local, xx_local
 PetscSection                       :: section
 PetscScalar, dimension(:), pointer :: x_scal, pf_scal
 PetscScalar,                target :: f_scal(cellDof) 
 PetscReal                          :: detJ, IcellJMat(dimPlex,dimPlex)
 PetscReal,                  target :: v0(dimPlex), cellJ(dimPlex*dimPlex), &
                                                    invcellJ(dimPlex*dimPlex)
 PetscReal,                 pointer :: pV0(:), pCellJ(:), pInvcellJ(:)
 PetscReal,                 pointer :: basisField(:), basisFieldDer(:)
 PetscInt                           :: cellStart, cellEnd, cell, field, face, &
                                       qPt, basis, comp, cidx
 PetscReal                          :: detFAvg
 PetscReal                          :: BMat(dimPlex*dimPlex,cellDof)
 PetscObject                        :: dummy
 PetscInt                           :: bcSize
 IS                                 :: bcPoints
 PetscErrorCode                     :: ierr
 
 pV0 => v0
 pCellJ => cellJ
 pInvcellJ => invcellJ
 call DMGetSection(dm_local,section,ierr); CHKERRQ(ierr)
 call DMGetDS(dm_local,prob,ierr); CHKERRQ(ierr)
 call PetscDSGetTabulation(prob,0,basisField,basisFieldDer,ierr)
 CHKERRQ(ierr)
 call DMPlexGetHeightStratum(dm_local,0,cellStart,cellEnd,ierr); CHKERRQ(ierr)
 call DMGetLocalVector(dm_local,x_local,ierr); CHKERRQ(ierr)
 !call VecWAXPY(x_local,1.0,xx_local,solution_local,ierr); CHKERRQ(ierr)
 call VecCopy(xx_local,x_local,ierr); CHKERRQ(ierr)
 do field = 1, dimPlex; do face = 1, mesh_Nboundaries
   if (params%fieldBC%componentBC(field)%Mask(face)) then
     call DMGetStratumSize(dm_local,'Face Sets',mesh_boundaries(face),bcSize,ierr)
     if (bcSize > 0) then
       call DMGetStratumIS(dm_local,'Face Sets',mesh_boundaries(face),bcPoints,ierr)
       CHKERRQ(ierr)
       call utilities_projectBCValues(x_local,section,0,field-1,bcPoints, &
                                      0.0,params%fieldBC%componentBC(field)%Value(face),params%timeinc)
       call ISDestroy(bcPoints,ierr); CHKERRQ(ierr)
     endif
   endif
 enddo; enddo

!--------------------------------------------------------------------------------------------------
! evaluate field derivatives
 do cell = cellStart, cellEnd-1                                                                     !< loop over all elements 
   call DMPlexVecGetClosure(dm_local,section,x_local,cell,x_scal,ierr)                              !< get Dofs belonging to element
   CHKERRQ(ierr)
   call  DMPlexComputeCellGeometryAffineFEM(dm_local,cell,pV0,pCellJ,pInvcellJ,detJ,ierr) 
   CHKERRQ(ierr)
   IcellJMat = reshape(pInvcellJ,shape=[dimPlex,dimPlex])
   do qPt = 0, nQuadrature-1
     BMat = 0.0
     do basis = 0, nBasis-1
       do comp = 0, dimPlex-1
         cidx = basis*dimPlex+comp
         BMat(comp*dimPlex+1:(comp+1)*dimPlex,basis*dimPlex+comp+1) = &
           matmul(IcellJMat,basisFieldDer((qPt*nBasis*dimPlex+cidx  )*dimPlex+1: &
                                          (qPt*nBasis*dimPlex+cidx+1)*dimPlex   ))
       enddo
     enddo
     materialpoint_F(1:dimPlex,1:dimPlex,qPt+1,cell+1) = &
       math_I3 + &
       reshape(matmul(BMat,x_scal),shape=[dimPlex,dimPlex], order=[2,1])
   enddo    
   if (BBarStabilisation) then
     detFAvg = math_det33(sum(materialpoint_F(1:3,1:3,1:nQuadrature,cell+1),dim=3)/real(nQuadrature))
     do qPt = 1, nQuadrature
       materialpoint_F(1:dimPlex,1:dimPlex,qPt,cell+1) = &
         materialpoint_F(1:dimPlex,1:dimPlex,qPt,cell+1)* &
         (detFAvg/math_det33(materialpoint_F(1:3,1:3,qPt,cell+1)))**(1.0/real(dimPlex))
         
     enddo  
   endif
   call DMPlexVecRestoreClosure(dm_local,section,x_local,cell,x_scal,ierr)
   CHKERRQ(ierr) 
 enddo 

!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
 call Utilities_constitutiveResponse(params%timeinc,P_av,ForwardData)
 call MPI_Allreduce(MPI_IN_PLACE,terminallyIll,1,MPI_LOGICAL,MPI_LOR,PETSC_COMM_WORLD,ierr)
 ForwardData = .false.

!--------------------------------------------------------------------------------------------------
! integrating residual 
 do cell = cellStart, cellEnd-1                                                                     !< loop over all elements 
   call DMPlexVecGetClosure(dm_local,section,x_local,cell,x_scal,ierr)                              !< get Dofs belonging to element
   CHKERRQ(ierr)
   call  DMPlexComputeCellGeometryAffineFEM(dm_local,cell,pV0,pCellJ,pInvcellJ,detJ,ierr) 
   CHKERRQ(ierr)
   IcellJMat = reshape(pInvcellJ,shape=[dimPlex,dimPlex])
   f_scal = 0.0
   do qPt = 0, nQuadrature-1
     BMat = 0.0
     do basis = 0, nBasis-1
       do comp = 0, dimPlex-1
         cidx = basis*dimPlex+comp
         BMat(comp*dimPlex+1:(comp+1)*dimPlex,basis*dimPlex+comp+1) = &
           matmul(IcellJMat,basisFieldDer((qPt*nBasis*dimPlex+cidx  )*dimPlex+1: &
                                          (qPt*nBasis*dimPlex+cidx+1)*dimPlex   ))
       enddo
     enddo
     f_scal = f_scal + &
              matmul(transpose(BMat), &
                     reshape(transpose(materialpoint_P(1:dimPlex,1:dimPlex,qPt+1,cell+1)), &
                             shape=[dimPlex*dimPlex]))*qWeights(qPt+1)          
   enddo    
   f_scal = f_scal*abs(detJ)  
   pf_scal => f_scal
   call DMPlexVecSetClosure(dm_local,section,f_local,cell,pf_scal,ADD_VALUES,ierr)
   CHKERRQ(ierr) 
   call DMPlexVecRestoreClosure(dm_local,section,x_local,cell,x_scal,ierr)
   CHKERRQ(ierr) 
 enddo
 call DMRestoreLocalVector(dm_local,x_local,ierr); CHKERRQ(ierr)
 
end subroutine FEM_mech_formResidual


!--------------------------------------------------------------------------------------------------
!> @brief forms the FEM stiffness matrix
!--------------------------------------------------------------------------------------------------
subroutine FEM_mech_formJacobian(dm_local,xx_local,Jac_pre,Jac,dummy,ierr)
 use numerics, only: &
   BBarStabilisation
 use homogenization, only: &
   materialpoint_dPdF, &
   materialpoint_F
 use math, only: &
   math_inv33, &
   math_identity2nd, &
   math_det33
 use FEM_utilities, only: &
   utilities_projectBCValues

 implicit none

 DM                                   :: dm_local
 PetscDS                              :: prob
 Vec                                  :: x_local, xx_local
 Mat                                  :: Jac_pre, Jac
 PetscSection                         :: section, gSection
 PetscReal                            :: detJ, IcellJMat(dimPlex,dimPlex)
 PetscReal,                    target :: v0(dimPlex), cellJ(dimPlex*dimPlex), &
                                                      invcellJ(dimPlex*dimPlex)
 PetscReal,                   pointer :: pV0(:), pCellJ(:), pInvcellJ(:)
 PetscReal,   dimension(:),   pointer :: basisField, basisFieldDer
 PetscInt                             :: cellStart, cellEnd, cell, field, face, &
                                         qPt, basis, comp, cidx
 PetscScalar,                  target :: K_e   (cellDof,cellDof), &
                                         K_eA  (cellDof,cellDof), &  
                                         K_eB  (cellDof,cellDof), &  
                                         K_eVec(cellDof*cellDof)
 PetscReal                            :: BMat   (dimPlex*dimPlex,cellDof), &
                                         BMatAvg(dimPlex*dimPlex,cellDof), &
                                         MatA   (dimPlex*dimPlex,cellDof), &
                                         MatB   (1              ,cellDof)
 PetscScalar, dimension(:),   pointer :: pK_e, x_scal
 PetscReal,   dimension(3,3)          :: F = math_I3, FAvg, FInv
 PetscObject                          :: dummy
 PetscInt                             :: bcSize
 IS                                   :: bcPoints
 PetscErrorCode                       :: ierr
 
 pV0 => v0
 pCellJ => cellJ
 pInvcellJ => invcellJ
 call MatSetOption(Jac,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr); CHKERRQ(ierr)
 call MatSetOption(Jac,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr); CHKERRQ(ierr)
 call MatZeroEntries(Jac,ierr); CHKERRQ(ierr)
 call DMGetDS(dm_local,prob,ierr); CHKERRQ(ierr)
 call PetscDSGetTabulation(prob,0,basisField,basisFieldDer,ierr)
 call DMGetSection(dm_local,section,ierr); CHKERRQ(ierr)
 call DMGetGlobalSection(dm_local,gSection,ierr); CHKERRQ(ierr)
 
 call DMGetLocalVector(dm_local,x_local,ierr); CHKERRQ(ierr)
 call VecWAXPY(x_local,1.0,xx_local,solution_local,ierr); CHKERRQ(ierr)
 do field = 1, dimPlex; do face = 1, mesh_Nboundaries
   if (params%fieldBC%componentBC(field)%Mask(face)) then
     call DMGetStratumSize(dm_local,'Face Sets',mesh_boundaries(face),bcSize,ierr)
     if (bcSize > 0) then
       call DMGetStratumIS(dm_local,'Face Sets',mesh_boundaries(face),bcPoints,ierr)
       CHKERRQ(ierr)
       call utilities_projectBCValues(x_local,section,0,field-1,bcPoints, &
                                      0.0,params%fieldBC%componentBC(field)%Value(face),params%timeinc)
       call ISDestroy(bcPoints,ierr); CHKERRQ(ierr)
     endif
   endif
 enddo; enddo
 call DMPlexGetHeightStratum(dm_local,0,cellStart,cellEnd,ierr); CHKERRQ(ierr)
 do cell = cellStart, cellEnd-1                                                                     !< loop over all elements 
   call DMPlexVecGetClosure(dm_local,section,x_local,cell,x_scal,ierr)                              !< get Dofs belonging to element
   CHKERRQ(ierr)
   call  DMPlexComputeCellGeometryAffineFEM(dm_local,cell,pV0,pCellJ,pInvcellJ,detJ,ierr) 
   CHKERRQ(ierr)
   IcellJMat = reshape(pInvcellJ, shape = [dimPlex,dimPlex])
   K_eA = 0.0
   K_eB = 0.0
   MatB = 0.0
   FAvg = 0.0
   BMatAvg = 0.0
   do qPt = 0, nQuadrature-1
     BMat = 0.0
     do basis = 0, nBasis-1
       do comp = 0, dimPlex-1
         cidx = basis*dimPlex+comp
         BMat(comp*dimPlex+1:(comp+1)*dimPlex,basis*dimPlex+comp+1) = &
           matmul(IcellJMat,basisFieldDer((qPt*nBasis*dimPlex+cidx  )*dimPlex+1: &
                                          (qPt*nBasis*dimPlex+cidx+1)*dimPlex   ))
       enddo
     enddo
     MatA = matmul(reshape(reshape(materialpoint_dPdF(1:dimPlex,1:dimPlex,1:dimPlex,1:dimPlex,qPt+1,cell+1), &
                                   shape=[dimPlex,dimPlex,dimPlex,dimPlex], order=[2,1,4,3]), &
                           shape=[dimPlex*dimPlex,dimPlex*dimPlex]),BMat)*qWeights(qPt+1)        
     if (BBarStabilisation) then
       F(1:dimPlex,1:dimPlex) = reshape(matmul(BMat,x_scal),shape=[dimPlex,dimPlex])
       FInv = math_inv33(F)
       K_eA = K_eA + matmul(transpose(BMat),MatA)*math_det33(FInv)**(1.0/real(dimPlex))        
       K_eB = K_eB - &
              matmul(transpose(matmul(reshape(materialpoint_F(1:dimPlex,1:dimPlex,qPt+1,cell+1), &
                                              shape=[dimPlex*dimPlex,1]), &
                                      matmul(reshape(FInv(1:dimPlex,1:dimPlex), &
                                                     shape=[1,dimPlex*dimPlex],order=[2,1]),BMat))),MatA)
       MatB = MatB + &
              matmul(reshape(materialpoint_F(1:dimPlex,1:dimPlex,qPt+1,cell+1),shape=[1,dimPlex*dimPlex]),MatA)  
       FAvg = FAvg + F      
       BMatAvg = BMatAvg + BMat
     else
       K_eA = K_eA + matmul(transpose(BMat),MatA)
     endif  
   enddo  
   if (BBarStabilisation) then
     FInv = math_inv33(FAvg)
     K_e = K_eA*math_det33(FAvg/real(nQuadrature))**(1.0/real(dimPlex)) + &
           (matmul(matmul(transpose(BMatAvg), &
                          reshape(FInv(1:dimPlex,1:dimPlex),shape=[dimPlex*dimPlex,1],order=[2,1])),MatB) + &
            K_eB)/real(dimPlex)              
     
   else
     K_e = K_eA
   endif  
   K_e = K_e + eps*math_identity2nd(cellDof)
   K_eVec = reshape(K_e, [cellDof*cellDof])*abs(detJ) 
   pK_e => K_eVec      
   call DMPlexMatSetClosure(dm_local,section,gSection,Jac,cell,pK_e,ADD_VALUES,ierr)
   CHKERRQ(ierr)
   call DMPlexVecRestoreClosure(dm_local,section,x_local,cell,x_scal,ierr)
   CHKERRQ(ierr) 
 enddo 
 call MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyBegin(Jac_pre,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call MatAssemblyEnd(Jac_pre,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
 call DMRestoreLocalVector(dm_local,x_local,ierr); CHKERRQ(ierr)
 
!--------------------------------------------------------------------------------------------------
! apply boundary conditions 
 call DMPlexCreateRigidBody(dm_local,matnull,ierr); CHKERRQ(ierr) !MD: linker error
 call MatSetNullSpace(Jac,matnull,ierr); CHKERRQ(ierr)
 call MatSetNearNullSpace(Jac,matnull,ierr); CHKERRQ(ierr)
 call MatNullSpaceDestroy(matnull,ierr); CHKERRQ(ierr)

end subroutine FEM_mech_formJacobian

!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine FEM_mech_forward(guess,timeinc,timeinc_old,fieldBC)
 use FEM_utilities, only: &
   cutBack
 use homogenization, only: &
   materialpoint_F0, &
   materialpoint_F
 use FEM_utilities, only: &
   utilities_projectBCValues

 implicit none
 type(tFieldBC), intent(in) :: &
   fieldBC
 real(pReal), intent(in) :: &
   timeinc_old, &
   timeinc 
 logical, intent(in) :: &
   guess
 PetscInt                :: field, face
 DM                      :: dm_local
 Vec                     :: x_local
 PetscSection            :: section
 PetscInt                :: bcSize
 IS                      :: bcPoints
 PetscErrorCode          :: ierr   

!--------------------------------------------------------------------------------------------------
! forward last inc
 if (guess .and. .not. cutBack) then 
   ForwardData = .True.
   materialpoint_F0 = materialpoint_F
   call SNESGetDM(mech_snes,dm_local,ierr); CHKERRQ(ierr)                                                 !< retrieve mesh info from mech_snes into dm_local
   call DMGetSection(dm_local,section,ierr); CHKERRQ(ierr)
   call DMGetLocalVector(dm_local,x_local,ierr); CHKERRQ(ierr)
   call VecSet(x_local,0.0,ierr); CHKERRQ(ierr)
   call DMGlobalToLocalBegin(dm_local,solution,INSERT_VALUES,x_local,ierr)                       !< retrieve my partition of global solution vector
   CHKERRQ(ierr)
   call DMGlobalToLocalEnd(dm_local,solution,INSERT_VALUES,x_local,ierr)
   CHKERRQ(ierr)
   call VecAXPY(solution_local,1.0,x_local,ierr); CHKERRQ(ierr)
   do field = 1, dimPlex; do face = 1, mesh_Nboundaries
     if (fieldBC%componentBC(field)%Mask(face)) then
       call DMGetStratumSize(dm_local,'Face Sets',mesh_boundaries(face),bcSize,ierr)
       if (bcSize > 0) then
         call DMGetStratumIS(dm_local,'Face Sets',mesh_boundaries(face),bcPoints,ierr)
         CHKERRQ(ierr)
         call utilities_projectBCValues(solution_local,section,0,field-1,bcPoints, &
                                        0.0,fieldBC%componentBC(field)%Value(face),timeinc_old)
         call ISDestroy(bcPoints,ierr); CHKERRQ(ierr)
       endif
     endif
   enddo; enddo
   call DMRestoreLocalVector(dm_local,x_local,ierr); CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! update rate and forward last inc
   call VecCopy(solution,solution_rate,ierr); CHKERRQ(ierr)
   call VecScale(solution_rate,1.0/timeinc_old,ierr); CHKERRQ(ierr)
 endif
 call VecCopy(solution_rate,solution,ierr); CHKERRQ(ierr)
 call VecScale(solution,timeinc,ierr); CHKERRQ(ierr)
 
end subroutine FEM_mech_forward


!--------------------------------------------------------------------------------------------------
!> @brief reporting
!--------------------------------------------------------------------------------------------------
subroutine FEM_mech_converged(snes_local,PETScIter,xnorm,snorm,fnorm,reason,dummy,ierr)
 use numerics, only: &
   worldrank, &
   err_struct_tolAbs, &
   err_struct_tolRel
 use IO, only: &
   IO_intOut
 use FEsolving, only: &
   terminallyIll

 implicit none
 SNES :: snes_local
 PetscInt :: PETScIter
 PetscReal :: xnorm,snorm,fnorm,divTol
 SNESConvergedReason :: reason
 PetscObject :: dummy
 PetscErrorCode :: ierr

!--------------------------------------------------------------------------------------------------
! report
 divTol = max(maxval(abs(P_av(1:dimPlex,1:dimPlex)))*err_struct_tolRel,err_struct_tolAbs)
 call SNESConvergedDefault(snes_local,PETScIter,xnorm,snorm,fnorm/divTol,reason,dummy,ierr)
 CHKERRQ(ierr)
 if (terminallyIll) reason = SNES_DIVERGED_FUNCTION_DOMAIN
 if (worldrank == 0) then 
   write(6,'(1/,1x,a,a,i0,a,i0,f0.3)') trim(incInfo), &
                   ' @ Iteration ',PETScIter,' mechanical residual norm = ', &
                                                   int(fnorm/divTol),fnorm/divTol-int(fnorm/divTol)
   write(6,'(/,a,/,3(3(2x,f12.4,1x)/))',advance='no') ' Piola--Kirchhoff stress / MPa =',&
                                                       transpose(P_av)*1.e-6_pReal
   flush(6) 
 endif
 
end subroutine FEM_mech_converged


!--------------------------------------------------------------------------------------------------
!> @brief destroy routine
!--------------------------------------------------------------------------------------------------
subroutine FEM_mech_destroy()

 implicit none
 PetscErrorCode :: ierr

 call VecDestroy(solution,ierr); CHKERRQ(ierr)
 call VecDestroy(solution_rate,ierr); CHKERRQ(ierr)
 call SNESDestroy(mech_snes,ierr); CHKERRQ(ierr)

end subroutine FEM_mech_destroy

end module FEM_mech
