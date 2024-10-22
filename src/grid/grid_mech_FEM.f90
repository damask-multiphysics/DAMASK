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
  use IO
  use misc
  use CLI
  use HDF5
  use HDF5_utilities
  use math
  use rotations
  use spectral_utilities
  use grid_utilities
  use config
  use homogenization
  use discretization
  use discretization_grid
  use constants

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  implicit none(type,external)
#else
  implicit none
#endif
  private

  type(tSolutionParams) :: params

  type :: tNumerics
    integer :: &
      itmin, &                                                                                      !< minimum number of iterations
      itmax                                                                                         !< maximum number of iterations
    real(pREAL) :: &
      eps_div_atol, &                                                                               !< absolute tolerance for equilibrium
      eps_div_rtol, &                                                                               !< relative tolerance for equilibrium
      eps_stress_atol, &                                                                            !< absolute tolerance for fullfillment of stress BC
      eps_stress_rtol                                                                               !< relative tolerance for fullfillment of stress BC
  end type tNumerics

  type(tNumerics) :: num                                                                            ! numerics parameters. Better name?

!--------------------------------------------------------------------------------------------------
! PETSc data
  DM   :: DM_mech
  SNES :: SNES_mech
  Vec  :: u_PETSc, u_lastInc_PETSc, uDot_PETSc

!--------------------------------------------------------------------------------------------------
! common pointwise data
  real(pREAL), dimension(:,:,:,:,:), allocatable :: F, P_current, F_lastInc
  real(pREAL) :: detJ
  real(pREAL), dimension(3)   :: delta
  real(pREAL), dimension(3,8) :: BMat
  real(pREAL), dimension(8,8) :: HGMat

!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
  real(pREAL), dimension(3,3) :: &
    F_aimDot = 0.0_pREAL, &                                                                         !< assumed rate of average deformation gradient
    F_aim = math_I3, &                                                                              !< current prescribed deformation gradient
    F_aim_lastInc = math_I3, &                                                                      !< previous average deformation gradient
    P_av = 0.0_pREAL, &                                                                             !< average 1st Piola--Kirchhoff stress
    P_aim = 0.0_pREAL
  character(len=:), allocatable :: incInfo                                                          !< time and increment information
  real(pREAL), dimension(3,3,3,3) :: &
    C_volAvg = 0.0_pREAL, &                                                                         !< current volume average stiffness
    C_volAvgLastInc = 0.0_pREAL, &                                                                  !< previous volume average stiffness
    S = 0.0_pREAL                                                                                   !< current compliance (filled up with zeros)

  real(pREAL) :: &
    err_BC                                                                                          !< deviation from stress BC

  integer :: totalIter = 0                                                                          !< total iteration in current increment
  integer(kind(STATUS_OK)) :: status

  public :: &
    grid_mechanical_FEM_init, &
    grid_mechanical_FEM_solution, &
    grid_mechanical_FEM_forward, &
    grid_mechanical_FEM_updateCoords, &
    grid_mechanical_FEM_restartWrite

contains

!--------------------------------------------------------------------------------------------------
!> @brief Allocate all necessary fields and fill them with data, potentially from restart info.
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_FEM_init(num_grid)

  type(tDict), pointer, intent(in) :: num_grid

  real(pREAL), parameter :: HGCoeff = 0.0e-2_pREAL
  real(pREAL), parameter, dimension(4,8) :: &
    HGcomp = reshape([ 1.0_pREAL, 1.0_pREAL, 1.0_pREAL,-1.0_pREAL, &
                       1.0_pREAL,-1.0_pREAL,-1.0_pREAL, 1.0_pREAL, &
                      -1.0_pREAL, 1.0_pREAL,-1.0_pREAL, 1.0_pREAL, &
                      -1.0_pREAL,-1.0_pREAL, 1.0_pREAL,-1.0_pREAL, &
                      -1.0_pREAL,-1.0_pREAL, 1.0_pREAL, 1.0_pREAL, &
                      -1.0_pREAL, 1.0_pREAL,-1.0_pREAL,-1.0_pREAL, &
                       1.0_pREAL,-1.0_pREAL,-1.0_pREAL,-1.0_pREAL, &
                       1.0_pREAL, 1.0_pREAL, 1.0_pREAL, 1.0_pREAL], [4,8])
  real(pREAL), dimension(3,3,3,3) :: devNull
  real(pREAL), dimension(3,3,product(cells(1:2))*cells3) :: temp33n
  real(pREAL), dimension(3,product(cells(1:2))*cells3) :: temp3n
  PetscErrorCode :: err_PETSc
  integer(MPI_INTEGER_KIND) :: err_MPI
  PetscScalar, pointer, dimension(:,:,:,:) :: &
    u,u_lastInc
  integer(MPI_INTEGER_KIND), dimension(0:worldsize-1) :: cells3_global
  integer(HID_T) :: fileHandle, groupHandle
  type(tDict), pointer :: &
    num_grid_mech
  character(len=:), allocatable :: &
    extmsg, &
    petsc_options


  print'(/,1x,a)', '<<<+-  grid_mechanical_FEM init  -+>>>'; flush(IO_STDOUT)

!-------------------------------------------------------------------------------------------------
! read numerical parameters and do sanity checks
  num_grid_mech => num_grid%get_dict('mechanical',defaultVal=emptyDict)

  num%itmin           = num_grid_mech%get_asInt('N_iter_min',defaultVal=1)
  num%itmax           = num_grid_mech%get_asInt('N_iter_max',defaultVal=100)
  num%eps_div_atol    = num_grid_mech%get_asReal('eps_abs_div(P)',defaultVal=1.0e-4_pREAL)
  num%eps_div_rtol    = num_grid_mech%get_asReal('eps_rel_div(P)',defaultVal=5.0e-4_pREAL)
  num%eps_stress_atol = num_grid_mech%get_asReal('eps_abs_P',     defaultVal=1.0e3_pREAL)
  num%eps_stress_rtol = num_grid_mech%get_asReal('eps_rel_P',     defaultVal=1.0e-3_pREAL)

  extmsg = ''
  if (num%eps_div_atol <= 0.0_pREAL)             extmsg = trim(extmsg)//' eps_abs_div(P)'
  if (num%eps_div_rtol <= 0.0_pREAL)             extmsg = trim(extmsg)//' eps_rel_div(P)'
  if (num%eps_stress_atol <= 0.0_pREAL)          extmsg = trim(extmsg)//' eps_abs_P'
  if (num%eps_stress_rtol <= 0.0_pREAL)          extmsg = trim(extmsg)//' eps_rel_P'
  if (num%itmax < 1)                             extmsg = trim(extmsg)//' N_iter_max'
  if (num%itmin > num%itmax .or. num%itmin < 1)  extmsg = trim(extmsg)//' N_iter_min'

  if (extmsg /= '') call IO_error(301,ext_msg=trim(extmsg))

!--------------------------------------------------------------------------------------------------
! set default and user defined options for PETSc

  petsc_options = misc_prefixOptions('-snes_type newtonls -ksp_type fgmres -ksp_max_it 25 '// &
                                     num_grid_mech%get_asStr('PETSc_options',defaultVal='') ,'mechanical_')
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,petsc_options,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! allocate global fields
  allocate(F (3,3,cells(1),cells(2),cells3),source = 0.0_pREAL)
  allocate(P_current (3,3,cells(1),cells(2),cells3),source = 0.0_pREAL)
  allocate(F_lastInc (3,3,cells(1),cells(2),cells3),source = 0.0_pREAL)

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
  call SNESCreate(PETSC_COMM_WORLD,SNES_mech,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetOptionsPrefix(SNES_mech,'mechanical_',err_PETSc)
  CHKERRQ(err_PETSc)
  call MPI_Allgather(int(cells3,MPI_INTEGER_KIND),1_MPI_INTEGER_KIND,MPI_INTEGER,&
                     cells3_global,1_MPI_INTEGER_KIND,MPI_INTEGER,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call DMDACreate3d(PETSC_COMM_WORLD, &
         DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, &
         DMDA_STENCIL_BOX, &
         int(cells(1),pPETSCINT),int(cells(2),pPETSCINT),int(cells(3),pPETSCINT), &                 ! global cells
         1_pPETSCINT, 1_pPETSCINT, int(worldsize,pPETSCINT), &
         3_pPETSCINT, 1_pPETSCINT, &                                                                ! #dof (u, vector), ghost boundary width (domain overlap)
         [int(cells(1),pPETSCINT)],[int(cells(2),pPETSCINT)],int(cells3_global,pPETSCINT), &        ! local cells
         DM_mech,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMsetFromOptions(DM_mech,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMsetUp(DM_mech,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDASetUniformCoordinates(DM_mech,0.0_pREAL,geomSize(1),0.0_pREAL,geomSize(2),0.0_pREAL,geomSize(3),err_PETSc)
  CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(DM_mech,u_PETSc,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(DM_mech,u_lastInc_PETSc,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(DM_mech,uDot_PETSc,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMSNESSetFunctionLocal(DM_mech,formResidual,PETSC_NULL_SNES,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMSNESSetJacobianLocal(DM_mech,formJacobian,PETSC_NULL_SNES,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetConvergenceTest(SNES_mech,converged,PETSC_NULL_SNES,PETSC_NULL_FUNCTION,err_PETSc)    ! specify custom convergence check function "_converged"
  CHKERRQ(err_PETSc)
  call SNESSetMaxLinearSolveFailures(SNES_mech, huge(1_pPETSCINT), err_PETSc)                       ! ignore linear solve failures
  CHKERRQ(err_PETSc)
  call SNESSetDM(SNES_mech,DM_mech,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetFromOptions(SNES_mech,err_PETSc)                                                      ! pull it all together with additional cli arguments
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! init fields
  call VecSet(u_PETSc,0.0_pREAL,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecSet(u_lastInc_PETSc,0.0_pREAL,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecSet(uDot_PETSc   ,0.0_pREAL,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(DM_mech,u_PETSc,u,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(DM_mech,u_lastInc_PETSc,u_lastInc,err_PETSc)
  CHKERRQ(err_PETSc)

  delta = geomSize/real(cells,pREAL)                                                                ! grid spacing
  detJ = product(delta)                                                                             ! cell volume

  BMat = reshape(real([-delta(1)**(-1),-delta(2)**(-1),-delta(3)**(-1), &
                        delta(1)**(-1),-delta(2)**(-1),-delta(3)**(-1), &
                       -delta(1)**(-1), delta(2)**(-1),-delta(3)**(-1), &
                        delta(1)**(-1), delta(2)**(-1),-delta(3)**(-1), &
                       -delta(1)**(-1),-delta(2)**(-1), delta(3)**(-1), &
                        delta(1)**(-1),-delta(2)**(-1), delta(3)**(-1), &
                       -delta(1)**(-1), delta(2)**(-1), delta(3)**(-1), &
                        delta(1)**(-1), delta(2)**(-1), delta(3)**(-1)],pREAL), [3,8])/4.0_pREAL    ! shape function derivative matrix

  HGMat = matmul(transpose(HGcomp),HGcomp) &
        * HGCoeff*(delta(1)*delta(2) + delta(2)*delta(3) + delta(3)*delta(1))/16.0_pREAL            ! hourglass stabilization matrix

!--------------------------------------------------------------------------------------------------
! init fields
  restartRead: if (CLI_restartInc > 0) then
    print'(/,1x,a,1x,i0)', 'loading restart data of increment', CLI_restartInc

    fileHandle  = HDF5_openFile(getSolverJobName()//'_restart.hdf5','r')
    groupHandle = HDF5_openGroup(fileHandle,'solver')

    call HDF5_read(P_aim,groupHandle,'P_aim',.false.)
    call MPI_Bcast(P_aim,9_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
    call parallelization_chkerr(err_MPI)
    call HDF5_read(F_aim,groupHandle,'F_aim',.false.)
    call MPI_Bcast(F_aim,9_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
    call parallelization_chkerr(err_MPI)
    call HDF5_read(F_aim_lastInc,groupHandle,'F_aim_lastInc',.false.)
    call MPI_Bcast(F_aim_lastInc,9_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
    call parallelization_chkerr(err_MPI)
    call HDF5_read(F_aimDot,groupHandle,'F_aimDot',.false.)
    call MPI_Bcast(F_aimDot,9_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
    call parallelization_chkerr(err_MPI)
    call HDF5_read(temp33n,groupHandle,'F')
    F = reshape(temp33n,[3,3,cells(1),cells(2),cells3])
    call HDF5_read(temp33n,groupHandle,'F_lastInc')
    F_lastInc = reshape(temp33n,[3,3,cells(1),cells(2),cells3])
    call HDF5_read(temp3n,groupHandle,'u')
    u = reshape(temp3n,[3,cells(1),cells(2),cells3])
    call HDF5_read(temp3n,groupHandle,'u_lastInc')
    u_lastInc = reshape(temp3n,[3,cells(1),cells(2),cells3])

  elseif (CLI_restartInc == 0) then restartRead
    F_lastInc = spread(spread(spread(math_I3,3,cells(1)),4,cells(2)),5,cells3)                      ! initialize to identity
    F         = spread(spread(spread(math_I3,3,cells(1)),4,cells(2)),5,cells3)
  end if restartRead

  call utilities_updateCoords(F)
  call utilities_constitutiveResponse(status,P_current,P_av,C_volAvg,devNull, &                     ! stress field, stress avg, global average of stiffness and (min+max)/2
                                      F, &                                                          ! target F
                                      0.0_pREAL)                                                    ! time increment
  call DMDAVecRestoreArrayF90(DM_mech,u_PETSc,u,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecRestoreArrayF90(DM_mech,u_lastInc_PETSc,u_lastInc,err_PETSc)
  CHKERRQ(err_PETSc)

  restartRead2: if (CLI_restartInc > 0) then
    print'(1x,a,1x,i0)', 'loading additional restart data of increment', CLI_restartInc
    call HDF5_read(C_volAvg,groupHandle,'C_volAvg',.false.)
    call MPI_Bcast(C_volAvg,81_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
    call parallelization_chkerr(err_MPI)
    call HDF5_read(C_volAvgLastInc,groupHandle,'C_volAvgLastInc',.false.)
    call MPI_Bcast(C_volAvgLastInc,81_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
    call parallelization_chkerr(err_MPI)

    call HDF5_closeGroup(groupHandle)
    call HDF5_closeFile(fileHandle)

  end if restartRead2

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

  S = utilities_maskedCompliance(params%rotation_BC,params%stress_mask,C_volAvg)

  call SNESsolve(SNES_mech,PETSC_NULL_VEC,u_PETSc,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESGetConvergedReason(SNES_mech,reason,err_PETSc)
  CHKERRQ(err_PETSc)

  solution%converged = reason > 0
  solution%iterationsNeeded = totalIter
  P_aim = merge(P_av,P_aim,params%stress_mask)

end function grid_mechanical_FEM_solution


!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!> @details find new boundary conditions and best F estimate for end of current timestep
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_FEM_forward(cutBack,guess,Delta_t,Delta_t_old,t_remaining, &
                                       deformation_BC,stress_BC,rotation_BC)

  logical,                  intent(in) :: &
    cutBack, &
    guess
  real(pREAL),              intent(in) :: &
    Delta_t_old, &
    Delta_t, &
    t_remaining                                                                                     !< remaining time of current load case
  type(tBoundaryCondition), intent(in) :: &
    stress_BC, &
    deformation_BC
  type(tRotation),          intent(in) :: &
    rotation_BC

  PetscErrorCode :: err_PETSc


  if (cutBack) then
    C_volAvg = C_volAvgLastInc
  else
    C_volAvgLastInc    = C_volAvg

    F_aimDot = merge(merge(.0_pREAL,(F_aim-F_aim_lastInc)/Delta_t_old,stress_BC%mask),.0_pREAL,guess) ! estimate deformation rate for prescribed stress components
    F_aim_lastInc = F_aim

    !-----------------------------------------------------------------------------------------------
    ! calculate rate for aim
    if     (deformation_BC%myType=='L') then                                                        ! calculate F_aimDot from given L and current F
      F_aimDot = F_aimDot &
               + matmul(merge(.0_pREAL,deformation_BC%values,deformation_BC%mask),F_aim_lastInc)
    elseif (deformation_BC%myType=='dot_F') then                                                    ! F_aimDot is prescribed
      F_aimDot = F_aimDot &
               + merge(.0_pREAL,deformation_BC%values,deformation_BC%mask)
    elseif (deformation_BC%myType=='F') then                                                        ! aim at end of load case is prescribed
      F_aimDot = F_aimDot &
               + merge(.0_pREAL,(deformation_BC%values - F_aim_lastInc)/t_remaining,deformation_BC%mask)
    end if

    if (guess) then
      call VecWAXPY(uDot_PETSc,-1.0_pREAL,u_lastInc_PETSc,u_PETSc,err_PETSc)
      CHKERRQ(err_PETSc)
      call VecScale(uDot_PETSc,1.0_pREAL/Delta_t_old,err_PETSc)
      CHKERRQ(err_PETSc)
    else
      call VecSet(uDot_PETSc,0.0_pREAL,err_PETSc)
      CHKERRQ(err_PETSc)
    end if
    call VecCopy(u_PETSc,u_lastInc_PETSc,err_PETSc)
    CHKERRQ(err_PETSc)

    F_lastInc = F

  end if

!--------------------------------------------------------------------------------------------------
! update average and local deformation gradients
  F_aim = F_aim_lastInc + F_aimDot * Delta_t
  if (stress_BC%myType=='P')     P_aim = P_aim &
                                       + merge(.0_pREAL,(stress_BC%values - P_aim)/t_remaining,stress_BC%mask)*Delta_t
  if (stress_BC%myType=='dot_P') P_aim = P_aim &
                                       + merge(.0_pREAL,stress_BC%values,stress_BC%mask)*Delta_t

  call VecAXPY(u_PETSc,Delta_t,uDot_PETSc,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! set module wide available data
  params%stress_mask = stress_BC%mask
  params%rotation_BC = rotation_BC
  params%Delta_t     = Delta_t

end subroutine grid_mechanical_FEM_forward


!--------------------------------------------------------------------------------------------------
!> @brief Update coordinates.
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_FEM_updateCoords()

  call utilities_updateCoords(F)

end subroutine grid_mechanical_FEM_updateCoords


!--------------------------------------------------------------------------------------------------
!> @brief Write current solver and constitutive data for restart to file.
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_FEM_restartWrite()

  PetscErrorCode :: err_PETSc
  integer(HID_T) :: fileHandle, groupHandle
  PetscScalar, dimension(:,:,:,:), pointer :: u,u_lastInc


  call DMDAVecGetArrayReadF90(DM_mech,u_PETSc,u,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayReadF90(DM_mech,u_lastInc_PETSc,u_lastInc,err_PETSc)
  CHKERRQ(err_PETSc)

  print'(1x,a)', 'saving solver data required for restart'; flush(IO_STDOUT)

  fileHandle  = HDF5_openFile(getSolverJobName()//'_restart.hdf5','w')
  groupHandle = HDF5_addGroup(fileHandle,'solver')
  call HDF5_write(reshape(F,[3,3,product(cells(1:2))*cells3]),groupHandle,'F')
  call HDF5_write(reshape(F_lastInc,[3,3,product(cells(1:2))*cells3]),groupHandle,'F_lastInc')
  call HDF5_write(reshape(u,[3,product(cells(1:2))*cells3]),groupHandle,'u')
  call HDF5_write(reshape(u_lastInc,[3,product(cells(1:2))*cells3]),groupHandle,'u_lastInc')
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
  end if

  call DMDAVecRestoreArrayReadF90(DM_mech,u_PETSc,u,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecRestoreArrayReadF90(DM_mech,u_lastInc_PETSc,u_lastInc,err_PETSc)
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
  real(pREAL) :: &
    err_div, &
    divTol, &
    BCTol


  err_div = fnorm*sqrt(wgt)*geomSize(1)/scaledGeomSize(1)/detJ
  divTol = max(maxval(abs(P_av))*num%eps_div_rtol, num%eps_div_atol)
  BCTol  = max(maxval(abs(P_av))*num%eps_stress_rtol, num%eps_stress_atol)

  if (totalIter >= num%itmin .and. all([err_div/divTol, err_BC/BCTol] < 1.0_pREAL) &
       .and. status == STATUS_OK) then
    reason = 1
  elseif (totalIter >= num%itmax) then
    reason = -1
  else
    reason = 0
  end if

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
!> @brief Form the residual vector.
!--------------------------------------------------------------------------------------------------
subroutine formResidual(da_local,x_local, &
                        f_local,dummy,err_PETSc)

  DM          :: da_local
  Vec         :: x_local, f_local
  PetscObject :: dummy
  PetscErrorCode :: err_PETSc

  real(pREAL), pointer,dimension(:,:,:,:) :: x_scal, r
  real(pREAL), dimension(8,3) :: x_elem,  f_elem
  PetscInt :: i, ii, j, jj, k, kk, ctr, ce, &
    PETScIter, nfuncs
  integer(MPI_INTEGER_KIND) :: err_MPI
  real(pREAL), dimension(3,3,3,3) :: devNull


  call SNESGetNumberFunctionEvals(SNES_mech,nfuncs,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESGetIterationNumber(SNES_mech,PETScIter,err_PETSc)
  CHKERRQ(err_PETSc)


  if (nfuncs == 0 .and. PETScIter == 0) totalIter = 0                                              ! new increment

!--------------------------------------------------------------------------------------------------
! begin of new iteration
  newIteration: if (totalIter <= PETScIter) then
    totalIter = totalIter + 1
    print'(1x,a,3(a,i0))', trim(incInfo), ' @ Iteration ', num%itmin, '≤',totalIter, '≤', num%itmax
    if (any(dNeq(params%rotation_BC%asQuaternion(), real([1.0, 0.0, 0.0, 0.0],pREAL)))) &
      print'(/,1x,a,/,2(3(f12.7,1x)/),3(f12.7,1x))', &
      'deformation gradient aim (lab) =', transpose(params%rotation_BC%rotate(F_aim,active=.true.))
    print'(/,1x,a,/,2(3(f12.7,1x)/),3(f12.7,1x))', &
      'deformation gradient aim       =', transpose(F_aim)
    flush(IO_STDOUT)
  end if newIteration

!--------------------------------------------------------------------------------------------------
! get deformation gradient
  call DMDAVecGetArrayReadF90(da_local,x_local,x_scal,err_PETSc)
  CHKERRQ(err_PETSc)
  do k = cells3Offset+1, cells3Offset+cells3; do j = 1, cells(2); do i = 1, cells(1)
    ctr = 0
    do kk = -1, 0; do jj = -1, 0; do ii = -1, 0
      ctr = ctr + 1
      x_elem(ctr,1:3) = x_scal(0:2,i+ii,j+jj,k+kk)
    end do; end do; end do
    F(1:3,1:3,i,j,k-cells3Offset) = params%rotation_BC%rotate(F_aim,active=.true.) + transpose(matmul(BMat,x_elem))
  end do; end do; end do
  call DMDAVecRestoreArrayReadF90(da_local,x_local,x_scal,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
  call utilities_constitutiveResponse(status,P_current,&
                                      P_av,C_volAvg,devNull, &
                                      F,params%Delta_t,params%rotation_BC)
  call MPI_Allreduce(MPI_IN_PLACE,status,1_MPI_INTEGER_KIND,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

!--------------------------------------------------------------------------------------------------
! stress BC handling
  F_aim = F_aim - math_mul3333xx33(S, P_av - P_aim)                                                 ! S = 0.0 for no bc
  err_BC = maxval(abs(merge(.0_pREAL,P_av - P_aim,params%stress_mask)))

!--------------------------------------------------------------------------------------------------
! constructing residual
  call DMDAVecGetArrayF90(da_local,f_local,r,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayReadF90(da_local,x_local,x_scal,err_PETSc)
  CHKERRQ(err_PETSc)
  ce = 0
  r = 0.0_pREAL
  do k = cells3Offset+1, cells3Offset+cells3; do j = 1, cells(2); do i = 1, cells(1)
    ctr = 0
    do kk = -1, 0; do jj = -1, 0; do ii = -1, 0
      ctr = ctr + 1
      x_elem(ctr,1:3) = x_scal(0:2,i+ii,j+jj,k+kk)
    end do; end do; end do
    ce = ce + 1
    f_elem = matmul(transpose(BMat),transpose(P_current(1:3,1:3,i,j,k-cells3Offset)))*detJ + &
             matmul(HGMat,x_elem)*(homogenization_dPdF(1,1,1,1,ce) + &
                                   homogenization_dPdF(2,2,2,2,ce) + &
                                   homogenization_dPdF(3,3,3,3,ce))/3.0_pREAL
    ctr = 0
    do kk = -1, 0; do jj = -1, 0; do ii = -1, 0
      ctr = ctr + 1
      r(0:2,i+ii,j+jj,k+kk) = r(0:2,i+ii,j+jj,k+kk) + f_elem(ctr,1:3)
    end do; end do; end do
  end do; end do; end do
  call DMDAVecRestoreArrayReadF90(da_local,x_local,x_scal,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! applying boundary conditions
  if (cells3Offset == 0) then
    r(0:2,0,       0,       0) = 0.0_pREAL
    r(0:2,cells(1),0,       0) = 0.0_pREAL
    r(0:2,0,       cells(2),0) = 0.0_pREAL
    r(0:2,cells(1),cells(2),0) = 0.0_pREAL
  end if
  if (cells3+cells3Offset == cells(3)) then
    r(0:2,0,       0,       cells(3)) = 0.0_pREAL
    r(0:2,cells(1),0,       cells(3)) = 0.0_pREAL
    r(0:2,0,       cells(2),cells(3)) = 0.0_pREAL
    r(0:2,cells(1),cells(2),cells(3)) = 0.0_pREAL
  end if
  call DMDAVecRestoreArrayF90(da_local,f_local,r,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine formResidual


!--------------------------------------------------------------------------------------------------
!> @brief Form the FEM stiffness matrix.
!--------------------------------------------------------------------------------------------------
subroutine formJacobian(da_local,x_local,Jac_pre,Jac,dummy,err_PETSc)

  DM                                   :: da_local
  Vec                                  :: x_local
  Mat                                  :: Jac_pre, Jac
  PetscObject                          :: dummy
  PetscErrorCode                       :: err_PETSc

  MatStencil,dimension(4,24)           :: row, col
  real(pREAL),pointer,dimension(:,:,:,:) :: x_scal
  real(pREAL),dimension(24,24)         :: K_ele
  real(pREAL),dimension(9,24)          :: BMatFull
  PetscInt                             :: i, ii, j, jj, k, kk, ctr, ce
  PetscInt,dimension(3),parameter      :: rows = [0, 1, 2]
  real(pREAL)                          :: diag
  MatNullSpace                         :: matnull
  Vec                                  :: coordinates


  BMatFull = 0.0_pREAL
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
    end do; end do; end do
    row = col
    ce = ce + 1
    K_ele = 0.0_pREAL
    K_ele(1 :8 ,1 :8 ) = HGMat*(homogenization_dPdF(1,1,1,1,ce) + &
                                homogenization_dPdF(2,2,2,2,ce) + &
                                homogenization_dPdF(3,3,3,3,ce))/3.0_pREAL
    K_ele(9 :16,9 :16) = HGMat*(homogenization_dPdF(1,1,1,1,ce) + &
                                homogenization_dPdF(2,2,2,2,ce) + &
                                homogenization_dPdF(3,3,3,3,ce))/3.0_pREAL
    K_ele(17:24,17:24) = HGMat*(homogenization_dPdF(1,1,1,1,ce) + &
                                homogenization_dPdF(2,2,2,2,ce) + &
                                homogenization_dPdF(3,3,3,3,ce))/3.0_pREAL
    K_ele = K_ele + &
            matmul(transpose(BMatFull), &
                   matmul(reshape(reshape(homogenization_dPdF(1:3,1:3,1:3,1:3,ce), &
                                          shape=[3,3,3,3], order=[2,1,4,3]),shape=[9,9]),BMatFull))*detJ
    call MatSetValuesStencil(Jac,24_pPETScInt,row,24_pPetscInt,col,K_ele,ADD_VALUES,err_PETSc)
    CHKERRQ(err_PETSc)
  end do; end do; end do
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
  x_scal = reshape(discretization_IPcoords,[3,cells(1),cells(2),cells3])
  call DMDAVecRestoreArrayF90(da_local,coordinates,x_scal,err_PETSc)                                ! ToDo: use undeformed or deformed configuration?
  CHKERRQ(err_PETSc)
  call MatNullSpaceCreateRigidBody(coordinates,matnull,err_PETSc)                                   ! get rigid body deformation modes
  CHKERRQ(err_PETSc)
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
