!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Yi Hu, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Grid solver for mechanics: Spectral Galerkin
!--------------------------------------------------------------------------------------------------
module grid_mechanical_spectral_Galerkin
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
  use misc
  use IO
  use HDF5
  use HDF5_utilities
  use math
  use rotations
  use spectral_utilities
  use grid_utilities
  use homogenization
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
  Vec  :: F_PETSc
  Mat  :: Jac_PETSc

!--------------------------------------------------------------------------------------------------
! common pointwise data
  real(pREAL), dimension(:,:,:,:,:), allocatable ::  &
    F_lastInc, &                                                                                    !< field of previous compatible deformation gradients
    Fdot                                                                                            !< field of assumed rate of compatible deformation gradient

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
    C_minMaxAvg = 0.0_pREAL, &                                                                      !< current (min+max)/2 stiffness
    C_minMaxAvgLastInc = 0.0_pREAL, &                                                               !< previous (min+max)/2 stiffness
    C_minMaxAvgRestart = 0.0_pREAL, &                                                               !< (min+max)/2 stiffnes (restart)
    S = 0.0_pREAL                                                                                   !< current compliance (filled up with zeros)

  real(pREAL) :: &
    err_BC, &                                                                                       !< deviation from stress BC
    err_div                                                                                         !< RMS of div of P

  integer :: totalIter = 0                                                                          !< total iteration in current increment
  integer(kind(STATUS_OK)) :: status

  public :: &
    grid_mechanical_spectral_Galerkin_init, &
    grid_mechanical_spectral_Galerkin_solution, &
    grid_mechanical_spectral_Galerkin_forward, &
    grid_mechanical_spectral_Galerkin_updateCoords, &
    grid_mechanical_spectral_Galerkin_restartWrite

  ! Missing interfaces for some PETSc versions
  interface MatCreateShell
    subroutine MatCreateShell(comm,mloc,nloc,m,n,ctx,mat,ierr)
      use petscmat
      MPI_Comm :: comm
      PetscInt :: mloc,nloc,m,n
      Vec :: ctx
      Mat :: mat
      PetscErrorCode :: ierr
    end subroutine MatCreateShell
  end interface MatCreateShell

  interface MatShellSetContext
    subroutine MatShellSetContext(mat,ctx,ierr)
      use petscmat
      Mat :: mat
      Vec :: ctx
      PetscErrorCode :: ierr
    end subroutine MatShellSetContext
  end interface MatShellSetContext

  interface MatShellGetContext
    subroutine MatShellGetContext(mat,ctx,ierr)
      use petscmat
      Mat :: mat
      Vec, Pointer :: ctx
      PetscErrorCode :: ierr
    end subroutine MatShellGetContext
  end interface MatShellGetContext

  interface MatShellSetOperation
    subroutine MatShellSetOperation(mat,op_num,op_callback,ierr)
      use petscmat
      Mat :: mat
      PetscEnum :: op_num
      external :: op_callback
      PetscErrorCode :: ierr
    end subroutine MatShellSetOperation
  end interface MatShellSetOperation

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<22)
  interface SNESSetJacobian
    subroutine SNESSetJacobian(snes_mech,A,P,jac_callback,ctx,ierr)
      use petscsnes
      SNES :: snes_mech
      Mat :: A, P
      external :: jac_callback
      integer :: ctx
      PetscErrorCode :: ierr
    end subroutine SNESSetJacobian
  end interface SNESSetJacobian
#endif

  interface SNESSetUpdate
    subroutine SNESSetUpdate(snes_mech,upd_callback,ierr)
      use petscsnes
      SNES :: snes_mech
      external :: upd_callback
      PetscErrorCode :: ierr
    end subroutine SNESSetUpdate
  end interface SNESSetUpdate

contains

!--------------------------------------------------------------------------------------------------
!> @brief Allocate all necessary fields and fill them with data, potentially from restart info.
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_spectral_Galerkin_init(num_grid)

  type(tDict), pointer, intent(in) :: num_grid

  real(pREAL), dimension(3,3,cells(1),cells(2),cells3) :: P
  PetscErrorCode :: err_PETSc
  integer(MPI_INTEGER_KIND) :: err_MPI
  real(pREAL), pointer, dimension(:,:,:,:) :: &
    F                                                                                               ! pointer to solution data
  integer(MPI_INTEGER_KIND), dimension(0:worldsize-1) :: cells3_global
  real(pREAL), dimension(3,3,product(cells(1:2))*cells3) :: temp33n
  integer(HID_T) :: fileHandle, groupHandle
  type(tDict), pointer :: &
    num_grid_fft, &
    num_grid_mech
  character(len=:), allocatable :: &
    extmsg, &
    petsc_options


  print'(/,1x,a)', '<<<+-  grid_mechanical_spectral_Galerkin init  -+>>>'; flush(IO_STDOUT)

  print'(/,1x,a)', 'J. Vondřejc et al., Computers & Mathematics with Applications 68(3):156–173, 2014'
  print'(  1x,a)', 'https://doi.org/10.1016/j.camwa.2014.05.014'//IO_EOL

  print'(  1x,a)', 'T.W.J. de Geus et al., Computer Methods in Applied Mechanics and Engineering 318:412–430, 2017'
  print'(  1x,a)', 'https://doi.org/10.1016/j.cma.2016.12.032'//IO_EOL

  print'(  1x,a)', 'S Lucarini et al., Modelling and Simulation in Materials Science and Engineering 30(2):023002, 2021'
  print'(  1x,a)', 'https://doi.org/10.1088/1361-651X/ac34e1'//IO_EOL

  print'(  1x,a)', 'µSpectre'
  print'(  1x,a)', 'https://gitlab.com/muspectre'

!-------------------------------------------------------------------------------------------------
! read numerical parameters and do sanity checks
  num_grid_fft =>  num_grid%get_dict('FFT',defaultVal=emptyDict)
  num_grid_mech => num_grid%get_dict('mechanical',defaultVal=emptyDict)

  num%itmin           = num_grid_mech%get_asInt('N_iter_min',defaultVal=1)
  num%itmax           = num_grid_mech%get_asInt('N_iter_max',defaultVal=100)
  num%eps_div_atol    = num_grid_mech%get_asReal('eps_abs_div(P)', defaultVal=1.0e-4_pREAL)
  num%eps_div_rtol    = num_grid_mech%get_asReal('eps_rel_div(P)', defaultVal=5.0e-4_pREAL)
  num%eps_stress_atol = num_grid_mech%get_asReal('eps_abs_P',      defaultVal=1.0e3_pREAL)
  num%eps_stress_rtol = num_grid_mech%get_asReal('eps_rel_P',      defaultVal=1.0e-3_pREAL)

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
  petsc_options = &
    misc_prefixOptions('-snes_type newtonls -ksp_type gmres -snes_linesearch_type bt '// &
                      num_grid_mech%get_asStr('PETSc_options',defaultVal=''), 'mechanical_')
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,petsc_options,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! allocate global fields
  allocate(F_lastInc(3,3,cells(1),cells(2),cells3),source = 0.0_pREAL)
  allocate(Fdot     (3,3,cells(1),cells(2),cells3),source = 0.0_pREAL)

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
         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &                                    ! cut off stencil at boundary
         DMDA_STENCIL_BOX, &                                                                        ! Moore (26) neighborhood around central point
         int(cells(1),pPETSCINT),int(cells(2),pPETSCINT),int(cells(3),pPETSCINT), &                 ! global cells
         1_pPETSCINT, 1_pPETSCINT, int(worldsize,pPETSCINT), &
         9_pPETSCINT, 0_pPETSCINT, &                                                                ! #dof (F, tensor), ghost boundary width (domain overlap)
         [int(cells(1),pPETSCINT)],[int(cells(2),pPETSCINT)],int(cells3_global,pPETSCINT), &        ! local cells
         DM_mech,err_PETSc)                                                                         ! handle, error
  CHKERRQ(err_PETSc)
  call DMsetFromOptions(DM_mech,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMsetUp(DM_mech,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMcreateGlobalVector(DM_mech,F_PETSc,err_PETSc)                                              ! global solution vector (cells x 9, i.e. every def grad tensor)
  CHKERRQ(err_PETSc)
  call DMDASNESsetFunctionLocal(DM_mech,INSERT_VALUES,formResidual,PETSC_NULL_SNES,err_PETSc)       ! residual vector of same shape as solution vector
  CHKERRQ(err_PETSc)

  ! Set shell Jacobian
  call MatCreateShell(PETSC_COMM_WORLD,&
                      int(9*product(cells(1:2))*cells3,pPETSCINT),&
                      int(9*product(cells(1:2))*cells3,pPETSCINT),&
                      int(9*product(cells(1:2))*sum(cells3_global),pPETSCINT),&
                      int(9*product(cells(1:2))*sum(cells3_global),pPETSCINT),&
                      F_PETSc,Jac_PETSc,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatShellSetOperation(Jac_PETSc,MATOP_MULT,GK_op,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetDM(SNES_mech,DM_mech,err_PETSc)                                                       ! associate snes with dm first, otherwise jac_shell not written
  CHKERRQ(err_PETSc)
  call SNESSetJacobian(SNES_mech,Jac_PETSc,Jac_PETSc,PETSC_NULL_FUNCTION,0,err_PETSc)
  CHKERRQ(err_PETSc)


  call SNESSetUpdate(SNES_mech,set_F_aim,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESsetConvergenceTest(SNES_mech,converged,PETSC_NULL_SNES,PETSC_NULL_FUNCTION,err_PETSc)    ! specify custom convergence check function "converged"
  CHKERRQ(err_PETSc)
  call SNESsetFromOptions(SNES_mech,err_PETSc)                                                      ! pull it all together with additional CLI arguments
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! init fields
  call DMDAVecGetArrayF90(DM_mech,F_PETSc,F,err_PETSc)                                              ! places pointer on PETSc data
  CHKERRQ(err_PETSc)

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
    F = reshape(temp33n,[9,cells(1),cells(2),cells3])
    call HDF5_read(temp33n,groupHandle,'F_lastInc')
    F_lastInc = reshape(temp33n,[3,3,cells(1),cells(2),cells3])

  elseif (CLI_restartInc == 0) then restartRead
    F_lastInc = spread(spread(spread(math_I3,3,cells(1)),4,cells(2)),5,cells3)                      ! initialize to identity
    F = reshape(F_lastInc,[9,cells(1),cells(2),cells3])
  end if restartRead

  call utilities_updateCoords(reshape(F,shape(F_lastInc)))
  call utilities_constitutiveResponse(status,P,P_av,C_volAvg,C_minMaxAvg, &                         ! stress field, stress avg, global average of stiffness and (min+max)/2
                                      reshape(F,shape(F_lastInc)), &                                ! target F
                                      0.0_pREAL)                                                    ! time increment
  call DMDAVecRestoreArrayF90(DM_mech,F_PETSc,F,err_PETSc)                                          ! deassociate pointer
  CHKERRQ(err_PETSc)

  restartRead2: if (CLI_restartInc > 0) then
    print'(1x,a,1x,i0)', 'loading additional restart data of increment', CLI_restartInc
    call HDF5_read(C_volAvg,groupHandle,'C_volAvg',.false.)
    call MPI_Bcast(C_volAvg,81_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
    call parallelization_chkerr(err_MPI)
    call HDF5_read(C_volAvgLastInc,groupHandle,'C_volAvgLastInc',.false.)
    call MPI_Bcast(C_volAvgLastInc,81_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
    call parallelization_chkerr(err_MPI)
    call HDF5_read(C_minMaxAvg,groupHandle,'C_minMaxAvg',.false.)
    call MPI_Bcast(C_minMaxAvg,81_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
    call parallelization_chkerr(err_MPI)

    call HDF5_closeGroup(groupHandle)
    call HDF5_closeFile(fileHandle)

  end if restartRead2

end subroutine grid_mechanical_spectral_Galerkin_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the Galerkin scheme with internal iterations
!--------------------------------------------------------------------------------------------------
function grid_mechanical_spectral_Galerkin_solution(incInfoIn) result(solution)

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
! update stiffness
  S = utilities_maskedCompliance(params%rotation_BC,params%stress_mask,C_volAvg)

  call SNESsolve(SNES_mech,PETSC_NULL_VEC,F_PETSc,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESGetConvergedReason(SNES_mech,reason,err_PETSc)
  CHKERRQ(err_PETSc)

  solution%converged = reason > 0
  solution%iterationsNeeded = totalIter
  P_aim = merge(P_av,P_aim,params%stress_mask)

end function grid_mechanical_spectral_Galerkin_solution


!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!> @details find new boundary conditions and best F estimate for end of current timestep
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_spectral_Galerkin_forward(cutBack,guess,Delta_t,Delta_t_old,t_remaining,&
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
  type(tRotation),           intent(in) :: &
    rotation_BC
  PetscErrorCode :: err_PETSc
  real(pREAL), pointer, dimension(:,:,:,:) :: F


  call DMDAVecGetArrayF90(DM_mech,F_PETSc,F,err_PETSc)
  CHKERRQ(err_PETSc)

  if (cutBack) then
    C_volAvg    = C_volAvgLastInc
    C_minMaxAvg = C_minMaxAvgLastInc
  else
    C_volAvgLastInc    = C_volAvg
    C_minMaxAvgLastInc = C_minMaxAvg

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

    Fdot = utilities_calculateRate(guess, &
                                   F_lastInc,reshape(F,[3,3,cells(1),cells(2),cells3]),Delta_t_old, &
                                   rotation_BC%rotate(F_aimDot,active=.true.))
    F_lastInc = reshape(F,[3,3,cells(1),cells(2),cells3])
  end if

!--------------------------------------------------------------------------------------------------
! update average and local deformation gradients
  F_aim = F_aim_lastInc + F_aimDot * Delta_t

  if (stress_BC%myType=='P')     P_aim = P_aim &
                                       + merge(.0_pREAL,(stress_BC%values - P_aim)/t_remaining,stress_BC%mask)*Delta_t
  if (stress_BC%myType=='dot_P') P_aim = P_aim &
                                       + merge(.0_pREAL,stress_BC%values,stress_BC%mask)*Delta_t

  F = reshape(utilities_forwardTensorField(Delta_t,F_lastInc,Fdot, &                                ! estimate of F at end of time+Delta_t that matches rotated F_aim on average
              rotation_BC%rotate(F_aim,active=.true.)),[9,cells(1),cells(2),cells3])
  call DMDAVecRestoreArrayF90(DM_mech,F_PETSc,F,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! set module wide available data
  params%stress_mask = stress_BC%mask
  params%rotation_BC = rotation_BC
  params%Delta_t     = Delta_t

end subroutine grid_mechanical_spectral_Galerkin_forward


!--------------------------------------------------------------------------------------------------
!> @brief Update coordinates.
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_spectral_Galerkin_updateCoords()

  PetscErrorCode :: err_PETSc
  real(pREAL), dimension(:,:,:,:), pointer :: F

  call DMDAVecGetArrayReadF90(DM_mech,F_PETSc,F,err_PETSc)
  CHKERRQ(err_PETSc)
  call utilities_updateCoords(reshape(F,[3,3,size(F,2),size(F,3),size(F,4)]))
  call DMDAVecRestoreArrayReadF90(DM_mech,F_PETSc,F,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine grid_mechanical_spectral_Galerkin_updateCoords


!--------------------------------------------------------------------------------------------------
!> @brief Write current solver and constitutive data for restart to file.
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_spectral_Galerkin_restartWrite()

  PetscErrorCode :: err_PETSc
  integer(HID_T) :: fileHandle, groupHandle
  real(pREAL), dimension(:,:,:,:), pointer :: F

  call DMDAVecGetArrayReadF90(DM_mech,F_PETSc,F,err_PETSc)
  CHKERRQ(err_PETSc)

  print'(1x,a)', 'saving solver data required for restart'; flush(IO_STDOUT)

  fileHandle  = HDF5_openFile(getSolverJobName()//'_restart.hdf5','w')
  groupHandle = HDF5_addGroup(fileHandle,'solver')
  call HDF5_write(reshape(F,[3,3,product(cells(1:2))*cells3]),groupHandle,'F')
  call HDF5_write(reshape(F_lastInc,[3,3,product(cells(1:2))*cells3]),groupHandle,'F_lastInc')
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
    call HDF5_write(C_minMaxAvgRestart,groupHandle,'C_minMaxAvg',.false.)
    call HDF5_closeGroup(groupHandle)
    call HDF5_closeFile(fileHandle)
  end if

  call DMDAVecRestoreArrayReadF90(DM_mech,F_PETSc,F,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine grid_mechanical_spectral_Galerkin_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief convergence check
!--------------------------------------------------------------------------------------------------
subroutine converged(snes_local,PETScIter,devNull1,devNull2,rhs_norm,reason,dummy,err_PETSc)

  SNES :: snes_local
  PetscInt,  intent(in) :: PETScIter
  PetscReal, intent(in) :: &
    devNull1, &
    devNull2, &
    rhs_norm
  SNESConvergedReason :: reason
  PetscObject :: dummy
  PetscErrorCode :: err_PETSc
  real(pREAL) :: &
    divTol, &
    BCTol

  divTol = max(maxval(abs(P_av))*num%eps_div_rtol, num%eps_div_atol)
  BCTol = max(maxval(abs(P_av))*num%eps_stress_rtol, num%eps_stress_atol)

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
!> @brief Construct the residual vector.
!--------------------------------------------------------------------------------------------------
subroutine formResidual(residual_subdomain, F, &
                        r, dummy, err_PETSc)

  DMDALocalInfo, dimension(DMDA_LOCAL_INFO_SIZE) :: &
    residual_subdomain                                                                              !< DMDA info (needs to be named "in" for macros like XRANGE to work)
  real(pREAL), dimension(3,3,cells(1),cells(2),cells3), intent(in) :: &
    F                                                                                               !< deformation gradient field
  real(pREAL), dimension(3,3,cells(1),cells(2),cells3), intent(out) :: &
    r                                                                                               !< residuum field
  PetscObject :: dummy
  PetscErrorCode :: err_PETSc

  real(pREAL), dimension(3,3) :: dP_with_BC

  PetscInt :: &
    PETScIter, &
    nfuncs
  integer(MPI_INTEGER_KIND) :: err_MPI


  call SNESGetNumberFunctionEvals(SNES_mech,nfuncs,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESGetIterationNumber(SNES_mech,PETScIter,err_PETSc)
  CHKERRQ(err_PETSc)

  if (nfuncs == 0 .and. PETScIter == 0) totalIter = -1                                              ! new increment

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

  associate (P => r)
    call utilities_constitutiveResponse(status,P, &
                                        P_av,C_volAvg,C_minMaxAvg, &
                                        F,params%Delta_t,params%rotation_BC)
    call MPI_Allreduce(MPI_IN_PLACE,status,1_MPI_INTEGER_KIND,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,err_MPI)
    call parallelization_chkerr(err_MPI)
    err_div = utilities_divergenceRMS(P)
  end associate

  dP_with_BC = merge(.0_pREAL,P_av-P_aim,params%stress_mask)
  err_BC     = maxval(abs(dP_with_BC))

  ! Yi: TODO: rotation
  r = utilities_G_Convolution(r,params%stress_mask,dP_with_BC)

end subroutine formResidual


!--------------------------------------------------------------------------------------------------
!> @brief Matrix-free operation GK_op -> GK_op(dF) = Fourier_inv( G_hat : Fourier(K:dF) )
!--------------------------------------------------------------------------------------------------
subroutine GK_op(Jac,dF_global,output_global,err_PETSc)

  DM                                   :: dm_local
  Vec                                  :: dF_global
  Vec                                  :: dF_local
  Vec                                  :: output_global
  Mat                                  :: Jac
  PetscErrorCode                       :: err_PETSc

  real(pREAL), pointer,dimension(:,:,:,:) :: dF_scal, output_scal

  real(pREAL), dimension(3,3,cells(1),cells(2),cells3) :: &
   dF
  real(pREAL), dimension(3,3,cells(1),cells(2),cells3) :: &
   output

  integer :: i, j, k, ce


  call SNESGetDM(SNES_mech,dm_local,err_PETSc)
  CHKERRQ(err_PETSc)

  call DMGetLocalVector(dm_local,dF_local,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGlobalToLocalBegin(dm_local,dF_global,INSERT_VALUES,dF_local,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMGlobalToLocalEnd(dm_local,dF_global,INSERT_VALUES,dF_local,err_PETSc)
  CHKERRQ(err_PETSc)

  call DMDAVecGetArrayReadF90(dm_local,dF_local,dF_scal,err_PETSc)
  CHKERRQ(err_PETSc)
  dF = reshape(dF_scal, [3,3,cells(1),cells(2),cells3])

  ! ===== K:dF operartor, i.e. dP = K:dF =====
  ce = 0
  do k = 1, cells3; do j = 1, cells(2); do i = 1, cells(1)
    ce = ce + 1
    output(1:3,1:3,i,j,k) = &
      math_mul3333xx33(homogenization_dPdF(1:3,1:3,1:3,1:3,ce), dF(1:3,1:3,i,j,k))
  end do; end do; end do

  ! ===== G* operator =====
  output = utilities_G_Convolution(output,params%stress_mask)

  call DMDAVecGetArrayF90(dm_local,output_global,output_scal,err_PETSc)
  CHKERRQ(err_PETSc)
  output_scal = reshape(output, [9,cells(1),cells(2),cells3])
  call DMDAVecRestoreArrayF90(dm_local,output_global,output_scal,err_PETSc)
  CHKERRQ(err_PETSc)

  call DMDAVecRestoreArrayF90(dm_local,dF_local,dF_scal,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMRestoreLocalVector(dm_local,dF_local,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine GK_op


!--------------------------------------------------------------------------------------------------
!> @brief Update F_aim only in newton step (not in line search step)
!--------------------------------------------------------------------------------------------------
subroutine set_F_aim(snes, step, ierr)
  SNES      :: snes
  PetscInt  :: step ! curr completed petsc iter
  PetscErrorCode, intent(out) :: ierr


  real(pREAL), dimension(3,3) :: &
    deltaF_aim

  deltaF_aim = math_mul3333xx33(S, P_av - P_aim)
  F_aim = F_aim - deltaF_aim

  ierr = 0_pPETSCERRORCODE

end subroutine set_F_aim

end module grid_mechanical_spectral_Galerkin
