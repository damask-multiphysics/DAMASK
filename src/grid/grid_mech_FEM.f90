!--------------------------------------------------------------------------------------------------
!> @author Arko Jyoti Bhattacharjee, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Grid solver for mechanics: FEM
!--------------------------------------------------------------------------------------------------
module grid_mechanical_FEM
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscdmda.h>
  use PETScDMDA
  use PETScSNES
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  use MPI_f08
#endif

  use prec
  use parallelization
  use CLI
  use IO
  use HDF5
  use HDF5_utilities
  use math
  use rotations
  use spectral_utilities
  use config
  use homogenization
  use discretization
  use discretization_grid

  implicit none
  private

  type(tSolutionParams) :: params

  type :: tNumerics
    integer :: &
      itmin, &                                                                                      !< minimum number of iterations
      itmax                                                                                         !< maximum number of iterations
    real(pReal) :: &
      eps_div_atol, &                                                                               !< absolute tolerance for equilibrium
      eps_div_rtol, &                                                                               !< relative tolerance for equilibrium
      eps_stress_atol, &                                                                            !< absolute tolerance for fullfillment of stress BC
      eps_stress_rtol                                                                               !< relative tolerance for fullfillment of stress BC
  end type tNumerics

  type(tNumerics) :: num                                                                            ! numerics parameters. Better name?

  logical :: debugRotation

!--------------------------------------------------------------------------------------------------
! PETSc data
  DM   :: mechanical_grid
  SNES :: SNES_mechanical
  Vec  :: solution_current, solution_lastInc, solution_rate

!--------------------------------------------------------------------------------------------------
! common pointwise data
  real(pReal), dimension(:,:,:,:,:), allocatable :: F, P_current, F_lastInc
  real(pReal) :: detJ
  real(pReal), dimension(3)   :: delta
  real(pReal), dimension(3,8) :: BMat
  real(pReal), dimension(8,8) :: HGMat

!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
  real(pReal), dimension(3,3) :: &
    F_aimDot = 0.0_pReal, &                                                                         !< assumed rate of average deformation gradient
    F_aim = math_I3, &                                                                              !< current prescribed deformation gradient
    F_aim_lastInc = math_I3, &                                                                      !< previous average deformation gradient
    P_av = 0.0_pReal, &                                                                             !< average 1st Piola--Kirchhoff stress
    P_aim = 0.0_pReal
  character(len=:), allocatable :: incInfo                                                          !< time and increment information
  real(pReal), dimension(3,3,3,3) :: &
    C_volAvg = 0.0_pReal, &                                                                         !< current volume average stiffness
    C_volAvgLastInc = 0.0_pReal, &                                                                  !< previous volume average stiffness
    S = 0.0_pReal                                                                                   !< current compliance (filled up with zeros)

  real(pReal) :: &
    err_BC                                                                                          !< deviation from stress BC

  integer :: &
    totalIter = 0                                                                                   !< total iteration in current increment

  public :: &
    grid_mechanical_FEM_init, &
    grid_mechanical_FEM_solution, &
    grid_mechanical_FEM_forward, &
    grid_mechanical_FEM_updateCoords, &
    grid_mechanical_FEM_restartWrite

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all necessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_FEM_init

  real(pReal), parameter :: HGCoeff = 0.0e-2_pReal
  real(pReal), parameter, dimension(4,8) :: &
    HGcomp = reshape([ 1.0_pReal, 1.0_pReal, 1.0_pReal,-1.0_pReal, &
                       1.0_pReal,-1.0_pReal,-1.0_pReal, 1.0_pReal, &
                      -1.0_pReal, 1.0_pReal,-1.0_pReal, 1.0_pReal, &
                      -1.0_pReal,-1.0_pReal, 1.0_pReal,-1.0_pReal, &
                      -1.0_pReal,-1.0_pReal, 1.0_pReal, 1.0_pReal, &
                      -1.0_pReal, 1.0_pReal,-1.0_pReal,-1.0_pReal, &
                       1.0_pReal,-1.0_pReal,-1.0_pReal,-1.0_pReal, &
                       1.0_pReal, 1.0_pReal, 1.0_pReal, 1.0_pReal], [4,8])
  real(pReal), dimension(3,3,3,3) :: devNull
  PetscErrorCode :: err_PETSc
  integer(MPI_INTEGER_KIND) :: err_MPI
  PetscScalar, pointer, dimension(:,:,:,:) :: &
    u_current,u_lastInc
  PetscInt, dimension(0:worldsize-1) :: localK
  integer(HID_T) :: fileHandle, groupHandle
  class(tNode), pointer :: &
    num_grid, &
    debug_grid
  character(len=pStringLen) :: &
    extmsg = ''

  print'(/,1x,a)', '<<<+-  grid_mechanical_FEM init  -+>>>'; flush(IO_STDOUT)

!-------------------------------------------------------------------------------------------------
! debugging options
  debug_grid => config_debug%get('grid',defaultVal=emptyList)
  debugRotation = debug_grid%contains('rotation')

!-------------------------------------------------------------------------------------------------
! read numerical parameters and do sanity checks
  num_grid => config_numerics%get('grid',defaultVal=emptyDict)

  num%eps_div_atol    = num_grid%get_asFloat('eps_div_atol',   defaultVal=1.0e-4_pReal)
  num%eps_div_rtol    = num_grid%get_asFloat('eps_div_rtol',   defaultVal=5.0e-4_pReal)
  num%eps_stress_atol = num_grid%get_asFloat('eps_stress_atol',defaultVal=1.0e3_pReal)
  num%eps_stress_rtol = num_grid%get_asFloat('eps_stress_rtol',defaultVal=1.0e-3_pReal)
  num%itmin           = num_grid%get_asInt  ('itmin',defaultVal=1)
  num%itmax           = num_grid%get_asInt  ('itmax',defaultVal=250)

  if (num%eps_div_atol <= 0.0_pReal)             extmsg = trim(extmsg)//' eps_div_atol'
  if (num%eps_div_rtol < 0.0_pReal)              extmsg = trim(extmsg)//' eps_div_rtol'
  if (num%eps_stress_atol <= 0.0_pReal)          extmsg = trim(extmsg)//' eps_stress_atol'
  if (num%eps_stress_rtol < 0.0_pReal)           extmsg = trim(extmsg)//' eps_stress_rtol'
  if (num%itmax <= 1)                            extmsg = trim(extmsg)//' itmax'
  if (num%itmin > num%itmax .or. num%itmin < 1)  extmsg = trim(extmsg)//' itmin'

  if (extmsg /= '') call IO_error(301,ext_msg=trim(extmsg))

!--------------------------------------------------------------------------------------------------
! set default and user defined options for PETSc
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS, &
                                '-mechanical_snes_type newtonls -mechanical_ksp_type fgmres &
                                &-mechanical_ksp_max_it 25', &
                                err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,num_grid%get_asString('petsc_options',defaultVal=''),err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! allocate global fields
  allocate(F (3,3,cells(1),cells(2),cells3),source = 0.0_pReal)
  allocate(P_current (3,3,cells(1),cells(2),cells3),source = 0.0_pReal)
  allocate(F_lastInc (3,3,cells(1),cells(2),cells3),source = 0.0_pReal)

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
  call SNESCreate(PETSC_COMM_WORLD,SNES_mechanical,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetOptionsPrefix(SNES_mechanical,'mechanical_',err_PETSc)
  CHKERRQ(err_PETSc)
  localK            = 0_pPetscInt
  localK(worldrank) = int(cells3,pPetscInt)
  call MPI_Allreduce(MPI_IN_PLACE,localK,worldsize,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  if(err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  call DMDACreate3d(PETSC_COMM_WORLD, &
         DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, &
         DMDA_STENCIL_BOX, &
         int(cells(1),pPetscInt),int(cells(2),pPetscInt),int(cells(3),pPetscInt), &                 ! global cells
         1_pPetscInt, 1_pPetscInt, int(worldsize,pPetscInt), &
         3_pPetscInt, 1_pPetscInt, &                                                                ! #dof (u, vector), ghost boundary width (domain overlap)
         [int(cells(1),pPetscInt)],[int(cells(2),pPetscInt)],localK, &                              ! local cells
         mechanical_grid,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMsetFromOptions(mechanical_grid,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMsetUp(mechanical_grid,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDASetUniformCoordinates(mechanical_grid,0.0_pReal,geomSize(1),0.0_pReal,geomSize(2),0.0_pReal,geomSize(3),err_PETSc)
  CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(mechanical_grid,solution_current,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(mechanical_grid,solution_lastInc,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(mechanical_grid,solution_rate   ,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMSNESSetFunctionLocal(mechanical_grid,formResidual,PETSC_NULL_SNES,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMSNESSetJacobianLocal(mechanical_grid,formJacobian,PETSC_NULL_SNES,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetConvergenceTest(SNES_mechanical,converged,PETSC_NULL_SNES,PETSC_NULL_FUNCTION,err_PETSc) ! specify custom convergence check function "_converged"
  CHKERRQ(err_PETSc)
  call SNESSetMaxLinearSolveFailures(SNES_mechanical, huge(1_pPetscInt), err_PETSc)                 ! ignore linear solve failures
  CHKERRQ(err_PETSc)
  call SNESSetDM(SNES_mechanical,mechanical_grid,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetFromOptions(SNES_mechanical,err_PETSc)                                                ! pull it all together with additional cli arguments
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! init fields
  call VecSet(solution_current,0.0_pReal,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecSet(solution_lastInc,0.0_pReal,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecSet(solution_rate   ,0.0_pReal,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(mechanical_grid,solution_current,u_current,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(mechanical_grid,solution_lastInc,u_lastInc,err_PETSc)
  CHKERRQ(err_PETSc)

  delta = geomSize/real(cells,pReal)                                                                ! grid spacing
  detJ = product(delta)                                                                             ! cell volume

  BMat = reshape(real([-delta(1)**(-1),-delta(2)**(-1),-delta(3)**(-1), &
                        delta(1)**(-1),-delta(2)**(-1),-delta(3)**(-1), &
                       -delta(1)**(-1), delta(2)**(-1),-delta(3)**(-1), &
                        delta(1)**(-1), delta(2)**(-1),-delta(3)**(-1), &
                       -delta(1)**(-1),-delta(2)**(-1), delta(3)**(-1), &
                        delta(1)**(-1),-delta(2)**(-1), delta(3)**(-1), &
                       -delta(1)**(-1), delta(2)**(-1), delta(3)**(-1), &
                        delta(1)**(-1), delta(2)**(-1), delta(3)**(-1)],pReal), [3,8])/4.0_pReal    ! shape function derivative matrix

  HGMat = matmul(transpose(HGcomp),HGcomp) &
        * HGCoeff*(delta(1)*delta(2) + delta(2)*delta(3) + delta(3)*delta(1))/16.0_pReal            ! hourglass stabilization matrix

!--------------------------------------------------------------------------------------------------
! init fields
  restartRead: if (CLI_restartInc > 0) then
    print'(/,1x,a,i0,a)', 'reading restart data of increment ', CLI_restartInc, ' from file'

    fileHandle  = HDF5_openFile(getSolverJobName()//'_restart.hdf5','r')
    groupHandle = HDF5_openGroup(fileHandle,'solver')

    call HDF5_read(P_aim,groupHandle,'P_aim',.false.)
    call MPI_Bcast(P_aim,9_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
    if(err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
    call HDF5_read(F_aim,groupHandle,'F_aim',.false.)
    call MPI_Bcast(F_aim,9_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
    if(err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
    call HDF5_read(F_aim_lastInc,groupHandle,'F_aim_lastInc',.false.)
    call MPI_Bcast(F_aim_lastInc,9_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
    if(err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
    call HDF5_read(F_aimDot,groupHandle,'F_aimDot',.false.)
    call MPI_Bcast(F_aimDot,9_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
    if(err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
    call HDF5_read(F,groupHandle,'F')
    call HDF5_read(F_lastInc,groupHandle,'F_lastInc')
    call HDF5_read(u_current,groupHandle,'u')
    call HDF5_read(u_lastInc,groupHandle,'u_lastInc')

  elseif (CLI_restartInc == 0) then restartRead
    F_lastInc = spread(spread(spread(math_I3,3,cells(1)),4,cells(2)),5,cells3)                      ! initialize to identity
    F         = spread(spread(spread(math_I3,3,cells(1)),4,cells(2)),5,cells3)
  endif restartRead

  homogenization_F0 = reshape(F_lastInc, [3,3,product(cells(1:2))*cells3])                          ! set starting condition for homogenization_mechanical_response
  call utilities_updateCoords(F)
  call utilities_constitutiveResponse(P_current,P_av,C_volAvg,devNull, &                            ! stress field, stress avg, global average of stiffness and (min+max)/2
                                      F, &                                                          ! target F
                                      0.0_pReal)                                                    ! time increment
  call DMDAVecRestoreArrayF90(mechanical_grid,solution_current,u_current,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecRestoreArrayF90(mechanical_grid,solution_lastInc,u_lastInc,err_PETSc)
  CHKERRQ(err_PETSc)

  restartRead2: if (CLI_restartInc > 0) then
    print'(1x,a,i0,a)', 'reading more restart data of increment ', CLI_restartInc, ' from file'
    call HDF5_read(C_volAvg,groupHandle,'C_volAvg',.false.)
    call MPI_Bcast(C_volAvg,81_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
    if(err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
    call HDF5_read(C_volAvgLastInc,groupHandle,'C_volAvgLastInc',.false.)
    call MPI_Bcast(C_volAvgLastInc,81_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
    if(err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'

    call HDF5_closeGroup(groupHandle)
    call HDF5_closeFile(fileHandle)

  endif restartRead2

end subroutine grid_mechanical_FEM_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the FEM scheme with internal iterations
!--------------------------------------------------------------------------------------------------
function grid_mechanical_FEM_solution(incInfoIn) result(solution)

!--------------------------------------------------------------------------------------------------
! input data for solution
  character(len=*),            intent(in) :: &
    incInfoIn
  type(tSolutionState)                    :: &
    solution
!--------------------------------------------------------------------------------------------------
! PETSc Data
  PetscErrorCode :: err_PETSc
  SNESConvergedReason :: reason

  incInfo = incInfoIn

!--------------------------------------------------------------------------------------------------
! update stiffness (and gamma operator)
  S = utilities_maskedCompliance(params%rotation_BC,params%stress_mask,C_volAvg)

  call SNESsolve(SNES_mechanical,PETSC_NULL_VEC,solution_current,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESGetConvergedReason(SNES_mechanical,reason,err_PETSc)
  CHKERRQ(err_PETSc)

  solution%converged = reason > 0
  solution%iterationsNeeded = totalIter
  solution%termIll = terminallyIll
  terminallyIll = .false.
  P_aim = merge(P_av,P_aim,params%stress_mask)

end function grid_mechanical_FEM_solution


!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!> @details find new boundary conditions and best F estimate for end of current timestep
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_FEM_forward(cutBack,guess,Delta_t,Delta_t_old,t_remaining,&
                                 deformation_BC,stress_BC,rotation_BC)

  logical,                  intent(in) :: &
    cutBack, &
    guess
  real(pReal),              intent(in) :: &
    Delta_t_old, &
    Delta_t, &
    t_remaining                                                                                     !< remaining time of current load case
  type(tBoundaryCondition), intent(in) :: &
    stress_BC, &
    deformation_BC
  type(tRotation),          intent(in) :: &
    rotation_BC
  PetscErrorCode :: err_PETSc
  PetscScalar, pointer, dimension(:,:,:,:) :: &
    u_current,u_lastInc


  call DMDAVecGetArrayF90(mechanical_grid,solution_current,u_current,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(mechanical_grid,solution_lastInc,u_lastInc,err_PETSc)
  CHKERRQ(err_PETSc)

  if (cutBack) then
    C_volAvg = C_volAvgLastInc
  else
    C_volAvgLastInc    = C_volAvg

    F_aimDot = merge(merge(.0_pReal,(F_aim-F_aim_lastInc)/Delta_t_old,stress_BC%mask),.0_pReal,guess) ! estimate deformation rate for prescribed stress components
    F_aim_lastInc = F_aim

    !-----------------------------------------------------------------------------------------------
    ! calculate rate for aim
    if     (deformation_BC%myType=='L') then                                                        ! calculate F_aimDot from given L and current F
      F_aimDot = F_aimDot &
               + matmul(merge(.0_pReal,deformation_BC%values,deformation_BC%mask),F_aim_lastInc)
    elseif (deformation_BC%myType=='dot_F') then                                                    ! F_aimDot is prescribed
      F_aimDot = F_aimDot &
               + merge(.0_pReal,deformation_BC%values,deformation_BC%mask)
    elseif (deformation_BC%myType=='F') then                                                        ! aim at end of load case is prescribed
      F_aimDot = F_aimDot &
               + merge(.0_pReal,(deformation_BC%values - F_aim_lastInc)/t_remaining,deformation_BC%mask)
    endif

    if (guess) then
      call VecWAXPY(solution_rate,-1.0_pReal,solution_lastInc,solution_current,err_PETSc)
      CHKERRQ(err_PETSc)
      call VecScale(solution_rate,1.0_pReal/Delta_t_old,err_PETSc)
      CHKERRQ(err_PETSc)
    else
      call VecSet(solution_rate,0.0_pReal,err_PETSc)
      CHKERRQ(err_PETSc)
    endif
    call VecCopy(solution_current,solution_lastInc,err_PETSc)
    CHKERRQ(err_PETSc)

    F_lastInc = F

    homogenization_F0 = reshape(F, [3,3,product(cells(1:2))*cells3])
  endif

!--------------------------------------------------------------------------------------------------
! update average and local deformation gradients
  F_aim = F_aim_lastInc + F_aimDot * Delta_t
  if (stress_BC%myType=='P')     P_aim = P_aim &
                                       + merge(.0_pReal,(stress_BC%values - P_aim)/t_remaining,stress_BC%mask)*Delta_t
  if (stress_BC%myType=='dot_P') P_aim = P_aim &
                                       + merge(.0_pReal,stress_BC%values,stress_BC%mask)*Delta_t

  call VecAXPY(solution_current,Delta_t,solution_rate,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecRestoreArrayF90(mechanical_grid,solution_current,u_current,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecRestoreArrayF90(mechanical_grid,solution_lastInc,u_lastInc,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! set module wide available data
  params%stress_mask = stress_BC%mask
  params%rotation_BC = rotation_BC
  params%Delta_t     = Delta_t

end subroutine grid_mechanical_FEM_forward


!--------------------------------------------------------------------------------------------------
!> @brief Update coordinates
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_FEM_updateCoords

  call utilities_updateCoords(F)

end subroutine grid_mechanical_FEM_updateCoords


!--------------------------------------------------------------------------------------------------
!> @brief Write current solver and constitutive data for restart to file
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_FEM_restartWrite

  PetscErrorCode :: err_PETSc
  integer(HID_T) :: fileHandle, groupHandle
  PetscScalar, dimension(:,:,:,:), pointer :: u_current,u_lastInc


  call DMDAVecGetArrayF90(mechanical_grid,solution_current,u_current,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(mechanical_grid,solution_lastInc,u_lastInc,err_PETSc)
  CHKERRQ(err_PETSc)

  print'(1x,a)', 'writing solver data required for restart to file'; flush(IO_STDOUT)

  fileHandle  = HDF5_openFile(getSolverJobName()//'_restart.hdf5','w')
  groupHandle = HDF5_addGroup(fileHandle,'solver')
  call HDF5_write(F,groupHandle,'F')
  call HDF5_write(F_lastInc,groupHandle,'F_lastInc')
  call HDF5_write(u_current,groupHandle,'u')
  call HDF5_write(u_lastInc,groupHandle,'u_lastInc')
  call HDF5_closeGroup(groupHandle)
  call HDF5_closeFile(fileHandle)

  if (worldrank == 0) then
    fileHandle  = HDF5_openFile(getSolverJobName()//'_restart.hdf5','a',.false.)
    groupHandle = HDF5_openGroup(fileHandle,'solver')
    call HDF5_write(P_aim,groupHandle,'P_aim',.false.)
    call HDF5_write(F_aim,groupHandle,'F_aim',.false.)
    call HDF5_write(F_aim_lastInc,groupHandle,'F_aim_lastInc',.false.)
    call HDF5_write(F_aimDot,groupHandle,'F_aimDot',.false.)
    call HDF5_write(C_volAvg,groupHandle,'C_volAvg',.false.)
    call HDF5_write(C_volAvgLastInc,groupHandle,'C_volAvgLastInc',.false.)
    call HDF5_closeGroup(groupHandle)
    call HDF5_closeFile(fileHandle)
  endif

  call DMDAVecRestoreArrayF90(mechanical_grid,solution_current,u_current,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecRestoreArrayF90(mechanical_grid,solution_lastInc,u_lastInc,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine grid_mechanical_FEM_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief convergence check
!--------------------------------------------------------------------------------------------------
subroutine converged(snes_local,PETScIter,devNull1,devNull2,fnorm,reason,dummy,err_PETSc)

  SNES :: snes_local
  PetscInt,  intent(in) :: PETScIter
  PetscReal, intent(in) :: &
    devNull1, &
    devNull2, &
    fnorm
  SNESConvergedReason :: reason
  PetscObject :: dummy
  PetscErrorCode :: err_PETSc
  real(pReal) :: &
    err_div, &
    divTol, &
    BCTol

  err_div = fnorm*sqrt(wgt)*geomSize(1)/scaledGeomSize(1)/detJ
  divTol = max(maxval(abs(P_av))*num%eps_div_rtol, num%eps_div_atol)
  BCTol  = max(maxval(abs(P_av))*num%eps_stress_rtol, num%eps_stress_atol)

  if ((totalIter >= num%itmin .and. all([err_div/divTol, err_BC/BCTol] < 1.0_pReal)) &
       .or. terminallyIll) then
    reason = 1
  elseif (totalIter >= num%itmax) then
    reason = -1
  else
    reason = 0
  endif

  print'(/,1x,a)', '... reporting .............................................................'
  print'(/,1x,a,f12.2,a,es8.2,a,es9.2,a)', 'error divergence = ', &
          err_div/divTol,  ' (',err_div,' / m, tol = ',divTol,')'
  print'(1x,a,f12.2,a,es8.2,a,es9.2,a)',    'error stress BC  = ', &
          err_BC/BCTol,    ' (',err_BC, ' Pa,  tol = ',BCTol,')'
  print'(/,1x,a)', '==========================================================================='
  flush(IO_STDOUT)
  err_PETSc = 0

end subroutine converged


!--------------------------------------------------------------------------------------------------
!> @brief forms the residual vector
!--------------------------------------------------------------------------------------------------
subroutine formResidual(da_local,x_local, &
                        f_local,dummy,err_PETSc)

  DM                   :: da_local
  Vec                  :: x_local, f_local
  PetscScalar, pointer,dimension(:,:,:,:) :: x_scal, r
  PetscScalar, dimension(8,3) :: x_elem,  f_elem
  PetscInt             :: i, ii, j, jj, k, kk, ctr, ele
  PetscInt :: &
    PETScIter, &
    nfuncs
  PetscObject :: dummy
  PetscErrorCode :: err_PETSc
  integer(MPI_INTEGER_KIND) :: err_MPI
  real(pReal), dimension(3,3,3,3) :: devNull

  call SNESGetNumberFunctionEvals(SNES_mechanical,nfuncs,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESGetIterationNumber(SNES_mechanical,PETScIter,err_PETSc)
  CHKERRQ(err_PETSc)

  if (nfuncs == 0 .and. PETScIter == 0) totalIter = -1                                              ! new increment

!--------------------------------------------------------------------------------------------------
! begin of new iteration
  newIteration: if (totalIter <= PETScIter) then
    totalIter = totalIter + 1
    print'(1x,a,3(a,i0))', trim(incInfo), ' @ Iteration ', num%itmin, '≤',totalIter+1, '≤', num%itmax
    if (debugRotation) print'(/,1x,a,/,2(3(f12.7,1x)/),3(f12.7,1x))', &
      'deformation gradient aim (lab) =', transpose(params%rotation_BC%rotate(F_aim,active=.true.))
    print'(/,1x,a,/,2(3(f12.7,1x)/),3(f12.7,1x))', &
      'deformation gradient aim       =', transpose(F_aim)
    flush(IO_STDOUT)
  endif newIteration

!--------------------------------------------------------------------------------------------------
! get deformation gradient
  call DMDAVecGetArrayF90(da_local,x_local,x_scal,err_PETSc)
  CHKERRQ(err_PETSc)
  do k = cells3Offset+1, cells3Offset+cells3; do j = 1, cells(2); do i = 1, cells(1)
    ctr = 0
    do kk = -1, 0; do jj = -1, 0; do ii = -1, 0
      ctr = ctr + 1
      x_elem(ctr,1:3) = x_scal(0:2,i+ii,j+jj,k+kk)
    enddo; enddo; enddo
    F(1:3,1:3,i,j,k-cells3Offset) = params%rotation_BC%rotate(F_aim,active=.true.) + transpose(matmul(BMat,x_elem))
  enddo; enddo; enddo
  call DMDAVecRestoreArrayF90(da_local,x_local,x_scal,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
  call utilities_constitutiveResponse(P_current,&
                                      P_av,C_volAvg,devNull, &
                                      F,params%Delta_t,params%rotation_BC)
  call MPI_Allreduce(MPI_IN_PLACE,terminallyIll,1_MPI_INTEGER_KIND,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,err_MPI)
  if(err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'

!--------------------------------------------------------------------------------------------------
! stress BC handling
  F_aim = F_aim - math_mul3333xx33(S, P_av - P_aim)                                                 ! S = 0.0 for no bc
  err_BC = maxval(abs(merge(.0_pReal,P_av - P_aim,params%stress_mask)))

!--------------------------------------------------------------------------------------------------
! constructing residual
  call VecSet(f_local,0.0_pReal,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(da_local,f_local,r,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(da_local,x_local,x_scal,err_PETSc)
  CHKERRQ(err_PETSc)
  ele = 0
  do k = cells3Offset+1, cells3Offset+cells3; do j = 1, cells(2); do i = 1, cells(1)
    ctr = 0
    do kk = -1, 0; do jj = -1, 0; do ii = -1, 0
      ctr = ctr + 1
      x_elem(ctr,1:3) = x_scal(0:2,i+ii,j+jj,k+kk)
    enddo; enddo; enddo
    ele = ele + 1
    f_elem = matmul(transpose(BMat),transpose(P_current(1:3,1:3,i,j,k-cells3Offset)))*detJ + &
             matmul(HGMat,x_elem)*(homogenization_dPdF(1,1,1,1,ele) + &
                                   homogenization_dPdF(2,2,2,2,ele) + &
                                   homogenization_dPdF(3,3,3,3,ele))/3.0_pReal
    ctr = 0
    do kk = -1, 0; do jj = -1, 0; do ii = -1, 0
      ctr = ctr + 1
      r(0:2,i+ii,j+jj,k+kk) = r(0:2,i+ii,j+jj,k+kk) + f_elem(ctr,1:3)
    enddo; enddo; enddo
  enddo; enddo; enddo
  call DMDAVecRestoreArrayF90(da_local,x_local,x_scal,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecRestoreArrayF90(da_local,f_local,r,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! applying boundary conditions
  call DMDAVecGetArrayF90(da_local,f_local,r,err_PETSc)
  CHKERRQ(err_PETSc)
  if (cells3Offset == 0) then
    r(0:2,0,       0,       0) = 0.0_pReal
    r(0:2,cells(1),0,       0) = 0.0_pReal
    r(0:2,0,       cells(2),0) = 0.0_pReal
    r(0:2,cells(1),cells(2),0) = 0.0_pReal
  end if
  if (cells3+cells3Offset == cells(3)) then
    r(0:2,0,       0,       cells(3)) = 0.0_pReal
    r(0:2,cells(1),0,       cells(3)) = 0.0_pReal
    r(0:2,0,       cells(2),cells(3)) = 0.0_pReal
    r(0:2,cells(1),cells(2),cells(3)) = 0.0_pReal
  end if
  call DMDAVecRestoreArrayF90(da_local,f_local,r,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine formResidual


!--------------------------------------------------------------------------------------------------
!> @brief forms the FEM stiffness matrix
!--------------------------------------------------------------------------------------------------
subroutine formJacobian(da_local,x_local,Jac_pre,Jac,dummy,err_PETSc)

  DM                                   :: da_local
  Vec                                  :: x_local, coordinates
  Mat                                  :: Jac_pre, Jac
  MatStencil,dimension(4,24)           :: row, col
  PetscScalar,pointer,dimension(:,:,:,:) :: x_scal
  PetscScalar,dimension(24,24)         :: K_ele
  PetscScalar,dimension(9,24)          :: BMatFull
  PetscInt                             :: i, ii, j, jj, k, kk, ctr, ce
  PetscInt,dimension(3),parameter      :: rows = [0, 1, 2]
  PetscScalar                          :: diag
  PetscObject                          :: dummy
  MatNullSpace                         :: matnull
  PetscErrorCode                       :: err_PETSc

  BMatFull = 0.0_pReal
  BMatFull(1:3,1 :8 ) = BMat
  BMatFull(4:6,9 :16) = BMat
  BMatFull(7:9,17:24) = BMat
  call MatSetOption(Jac,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatSetOption(Jac,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatZeroEntries(Jac,err_PETSc)
  CHKERRQ(err_PETSc)
  ce = 0
  do k = cells3Offset+1, cells3Offset+cells3; do j = 1, cells(2); do i = 1, cells(1)
    ctr = 0
    do kk = -1, 0; do jj = -1, 0; do ii = -1, 0
      ctr = ctr + 1
      col(MatStencil_i,ctr   ) = i+ii
      col(MatStencil_j,ctr   ) = j+jj
      col(MatStencil_k,ctr   ) = k+kk
      col(MatStencil_c,ctr   ) = 0
      col(MatStencil_i,ctr+8 ) = i+ii
      col(MatStencil_j,ctr+8 ) = j+jj
      col(MatStencil_k,ctr+8 ) = k+kk
      col(MatStencil_c,ctr+8 ) = 1
      col(MatStencil_i,ctr+16) = i+ii
      col(MatStencil_j,ctr+16) = j+jj
      col(MatStencil_k,ctr+16) = k+kk
      col(MatStencil_c,ctr+16) = 2
    enddo; enddo; enddo
    row = col
    ce = ce + 1
    K_ele = 0.0_pReal
    K_ele(1 :8 ,1 :8 ) = HGMat*(homogenization_dPdF(1,1,1,1,ce) + &
                                homogenization_dPdF(2,2,2,2,ce) + &
                                homogenization_dPdF(3,3,3,3,ce))/3.0_pReal
    K_ele(9 :16,9 :16) = HGMat*(homogenization_dPdF(1,1,1,1,ce) + &
                                homogenization_dPdF(2,2,2,2,ce) + &
                                homogenization_dPdF(3,3,3,3,ce))/3.0_pReal
    K_ele(17:24,17:24) = HGMat*(homogenization_dPdF(1,1,1,1,ce) + &
                                homogenization_dPdF(2,2,2,2,ce) + &
                                homogenization_dPdF(3,3,3,3,ce))/3.0_pReal
    K_ele = K_ele + &
            matmul(transpose(BMatFull), &
                   matmul(reshape(reshape(homogenization_dPdF(1:3,1:3,1:3,1:3,ce), &
                                          shape=[3,3,3,3], order=[2,1,4,3]),shape=[9,9]),BMatFull))*detJ
    call MatSetValuesStencil(Jac,24_pPETScInt,row,24_pPetscInt,col,K_ele,ADD_VALUES,err_PETSc)
    CHKERRQ(err_PETSc)
  enddo; enddo; enddo
  call MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatAssemblyBegin(Jac_pre,MAT_FINAL_ASSEMBLY,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatAssemblyEnd(Jac_pre,MAT_FINAL_ASSEMBLY,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! applying boundary conditions
  diag = (C_volAvg(1,1,1,1)/delta(1)**2 + C_volAvg(2,2,2,2)/delta(2)**2 + C_volAvg(3,3,3,3)/delta(3)**2) &
       * detJ
  call MatZeroRowsColumns(Jac,size(rows,kind=pPetscInt),rows,diag,PETSC_NULL_VEC,PETSC_NULL_VEC,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGetGlobalVector(da_local,coordinates,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(da_local,coordinates,x_scal,err_PETSc)
  CHKERRQ(err_PETSc)
  ce = 0
  do k = cells3Offset+1, cells3Offset+cells3; do j = 1, cells(2); do i = 1, cells(1)
    ce = ce + 1
    x_scal(0:2,i-1,j-1,k-1) = discretization_IPcoords(1:3,ce)
  enddo; enddo; enddo
  call DMDAVecRestoreArrayF90(da_local,coordinates,x_scal,err_PETSc)
  CHKERRQ(err_PETSc)                                                                                ! initialize to undeformed coordinates (ToDo: use ip coordinates)
  call MatNullSpaceCreateRigidBody(coordinates,matnull,err_PETSc)
  CHKERRQ(err_PETSc)                                                                                ! get rigid body deformation modes
  call DMRestoreGlobalVector(da_local,coordinates,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatSetNullSpace(Jac,matnull,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatSetNearNullSpace(Jac,matnull,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatNullSpaceDestroy(matnull,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine formJacobian

end module grid_mechanical_FEM
