!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Grid solver for mechanics: Spectral Polarisation
!--------------------------------------------------------------------------------------------------
module grid_mechanical_spectral_polarisation
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscdmda.h>
  use PETScdmda
  use PETScsnes

  use prec
  use parallelization
  use DAMASK_interface
  use IO
  use HDF5_utilities
  use math
  use rotations
  use spectral_utilities
  use config
  use homogenization
  use discretization_grid

  implicit none
  private

  type(tSolutionParams) :: params

  type :: tNumerics
    logical :: update_gamma                                                                         !< update gamma operator with current stiffness
    integer :: &
      itmin, &                                                                                      !< minimum number of iterations
      itmax                                                                                         !< maximum number of iterations
    real(pReal) :: &
      eps_div_atol, &                                                                               !< absolute tolerance for equilibrium
      eps_div_rtol, &                                                                               !< relative tolerance for equilibrium
      eps_curl_atol, &                                                                              !< absolute tolerance for compatibility
      eps_curl_rtol, &                                                                              !< relative tolerance for compatibility
      eps_stress_atol, &                                                                            !< absolute tolerance for fullfillment of stress BC
      eps_stress_rtol                                                                               !< relative tolerance for fullfillment of stress BC
    real(pReal) :: &
      alpha, &                                                                                      !< polarization scheme parameter 0.0 < alpha < 2.0. alpha = 1.0 ==> AL scheme, alpha = 2.0 ==> accelerated scheme
      beta                                                                                          !< polarization scheme parameter 0.0 < beta < 2.0. beta = 1.0 ==> AL scheme, beta = 2.0 ==> accelerated scheme
  end type tNumerics

  type(tNumerics) :: num                                                                            ! numerics parameters. Better name?

  logical :: debugRotation

!--------------------------------------------------------------------------------------------------
! PETSc data
  DM   :: da
  SNES :: snes
  Vec  :: solution_vec

!--------------------------------------------------------------------------------------------------
! common pointwise data
  real(pReal), dimension(:,:,:,:,:), allocatable :: &
    F_lastInc, &                                                                                    !< field of previous compatible deformation gradients
    F_tau_lastInc, &                                                                                !< field of previous incompatible deformation gradient
    Fdot, &                                                                                         !< field of assumed rate of compatible deformation gradient
    F_tauDot                                                                                        !< field of assumed rate of incopatible deformation gradient

!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
  real(pReal), dimension(3,3) :: &
    F_aimDot = 0.0_pReal, &                                                                         !< assumed rate of average deformation gradient
    F_aim = math_I3, &                                                                              !< current prescribed deformation gradient
    F_aim_lastInc = math_I3, &                                                                      !< previous average deformation gradient
    F_av = 0.0_pReal, &                                                                             !< average incompatible def grad field
    P_av = 0.0_pReal, &                                                                             !< average 1st Piola--Kirchhoff stress
    P_aim = 0.0_pReal
  character(len=:), allocatable :: incInfo                                                          !< time and increment information
  real(pReal), dimension(3,3,3,3) :: &
    C_volAvg = 0.0_pReal, &                                                                         !< current volume average stiffness
    C_volAvgLastInc = 0.0_pReal, &                                                                  !< previous volume average stiffness
    C_minMaxAvg = 0.0_pReal, &                                                                      !< current (min+max)/2 stiffness
    C_minMaxAvgLastInc = 0.0_pReal, &                                                               !< previous (min+max)/2 stiffness
    S = 0.0_pReal, &                                                                                !< current compliance (filled up with zeros)
    C_scale = 0.0_pReal, &
    S_scale = 0.0_pReal

  real(pReal) :: &
    err_BC, &                                                                                       !< deviation from stress BC
    err_curl, &                                                                                     !< RMS of curl of F
    err_div                                                                                         !< RMS of div of P

  integer :: &
    totalIter = 0                                                                                   !< total iteration in current increment

  public :: &
    grid_mechanical_spectral_polarisation_init, &
    grid_mechanical_spectral_polarisation_solution, &
    grid_mechanical_spectral_polarisation_forward, &
    grid_mechanical_spectral_polarisation_updateCoords, &
    grid_mechanical_spectral_polarisation_restartWrite

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all necessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_spectral_polarisation_init

  real(pReal), dimension(3,3,grid(1),grid(2),grid3) :: P
  PetscErrorCode :: ierr
  PetscScalar, pointer, dimension(:,:,:,:) :: &
    FandF_tau, &                                                                                    ! overall pointer to solution data
    F, &                                                                                            ! specific (sub)pointer
    F_tau                                                                                           ! specific (sub)pointer
  PetscInt, dimension(0:worldsize-1) :: localK
  integer(HID_T) :: fileHandle, groupHandle
  integer        :: fileUnit
  class (tNode), pointer :: &
    num_grid, &
    debug_grid

  print'(/,a)', ' <<<+-  grid_mechanical_spectral_polarization init  -+>>>'; flush(IO_STDOUT)

  print*, 'P. Shanthraj et al., International Journal of Plasticity 66:31–45, 2015'
  print*, 'https://doi.org/10.1016/j.ijplas.2014.02.006'

!-------------------------------------------------------------------------------------------------
! debugging options
  debug_grid => config_debug%get('grid',defaultVal=emptyList)
  debugRotation = debug_grid%contains('rotation')

!-------------------------------------------------------------------------------------------------
! read numerical parameters and do sanity checks
  num_grid => config_numerics%get('grid',defaultVal=emptyDict)

  num%update_gamma    = num_grid%get_asBool ('update_gamma',   defaultVal=.false.)
  num%eps_div_atol    = num_grid%get_asFloat('eps_div_atol',   defaultVal=1.0e-4_pReal)
  num%eps_div_rtol    = num_grid%get_asFloat('eps_div_rtol',   defaultVal=5.0e-4_pReal)
  num%eps_curl_atol   = num_grid%get_asFloat('eps_curl_atol',  defaultVal=1.0e-10_pReal)
  num%eps_curl_rtol   = num_grid%get_asFloat('eps_curl_rtol',  defaultVal=5.0e-4_pReal)
  num%eps_stress_atol = num_grid%get_asFloat('eps_stress_atol',defaultVal=1.0e3_pReal)
  num%eps_stress_rtol = num_grid%get_asFloat('eps_stress_rtol',defaultVal=1.0e-3_pReal)
  num%itmin           = num_grid%get_asInt  ('itmin',          defaultVal=1)
  num%itmax           = num_grid%get_asInt  ('itmax',          defaultVal=250)
  num%alpha           = num_grid%get_asFloat('alpha',          defaultVal=1.0_pReal)
  num%beta            = num_grid%get_asFloat('beta',           defaultVal=1.0_pReal)

  if (num%eps_div_atol <= 0.0_pReal)                      call IO_error(301,ext_msg='eps_div_atol')
  if (num%eps_div_rtol < 0.0_pReal)                       call IO_error(301,ext_msg='eps_div_rtol')
  if (num%eps_curl_atol <= 0.0_pReal)                     call IO_error(301,ext_msg='eps_curl_atol')
  if (num%eps_curl_rtol < 0.0_pReal)                      call IO_error(301,ext_msg='eps_curl_rtol')
  if (num%eps_stress_atol <= 0.0_pReal)                   call IO_error(301,ext_msg='eps_stress_atol')
  if (num%eps_stress_rtol < 0.0_pReal)                    call IO_error(301,ext_msg='eps_stress_rtol')
  if (num%itmax <= 1)                                     call IO_error(301,ext_msg='itmax')
  if (num%itmin > num%itmax .or. num%itmin < 1)           call IO_error(301,ext_msg='itmin')
  if (num%alpha <= 0.0_pReal .or. num%alpha >  2.0_pReal) call IO_error(301,ext_msg='alpha')
  if (num%beta < 0.0_pReal .or. num%beta > 2.0_pReal)     call IO_error(301,ext_msg='beta')

!--------------------------------------------------------------------------------------------------
! set default and user defined options for PETSc
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,'-mechanical_snes_type ngmres',ierr)
  CHKERRQ(ierr)
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,num_grid%get_asString('petsc_options',defaultVal=''),ierr)
  CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! allocate global fields
  allocate(F_lastInc    (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
  allocate(Fdot         (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
  allocate(F_tau_lastInc(3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
  allocate(F_tauDot     (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
  call SNESCreate(PETSC_COMM_WORLD,snes,ierr); CHKERRQ(ierr)
  call SNESSetOptionsPrefix(snes,'mechanical_',ierr);CHKERRQ(ierr)
  localK            = 0
  localK(worldrank) = grid3
  call MPI_Allreduce(MPI_IN_PLACE,localK,worldsize,MPI_INTEGER,MPI_SUM,PETSC_COMM_WORLD,ierr)
  call DMDACreate3d(PETSC_COMM_WORLD, &
         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &                                    ! cut off stencil at boundary
         DMDA_STENCIL_BOX, &                                                                        ! Moore (26) neighborhood around central point
         grid(1),grid(2),grid(3), &                                                                 ! global grid
         1 , 1, worldsize, &
         18, 0, &                                                                                   ! #dof (F tensor), ghost boundary width (domain overlap)
         [grid(1)],[grid(2)],localK, &                                                              ! local grid
         da,ierr)                                                                                   ! handle, error
  CHKERRQ(ierr)
  call SNESSetDM(snes,da,ierr); CHKERRQ(ierr)                                                       ! connect snes to da
  call DMsetFromOptions(da,ierr); CHKERRQ(ierr)
  call DMsetUp(da,ierr); CHKERRQ(ierr)
  call DMcreateGlobalVector(da,solution_vec,ierr); CHKERRQ(ierr)                                    ! global solution vector (grid x 18, i.e. every def grad tensor)
  call DMDASNESsetFunctionLocal(da,INSERT_VALUES,formResidual,PETSC_NULL_SNES,ierr)                 ! residual vector of same shape as solution vector
  CHKERRQ(ierr)
  call SNESsetConvergenceTest(snes,converged,PETSC_NULL_SNES,PETSC_NULL_FUNCTION,ierr)              ! specify custom convergence check function "converged"
  CHKERRQ(ierr)
  call SNESsetFromOptions(snes,ierr); CHKERRQ(ierr)                                                 ! pull it all together with additional CLI arguments

!--------------------------------------------------------------------------------------------------
! init fields
  call DMDAVecGetArrayF90(da,solution_vec,FandF_tau,ierr); CHKERRQ(ierr)                             ! places pointer on PETSc data
  F     => FandF_tau(0: 8,:,:,:)
  F_tau => FandF_tau(9:17,:,:,:)

  restartRead: if (interface_restartInc > 0) then
    print'(/,a,i0,a)', ' reading restart data of increment ', interface_restartInc, ' from file'

    fileHandle  = HDF5_openFile(getSolverJobName()//'_restart.hdf5','r')
    groupHandle = HDF5_openGroup(fileHandle,'solver')

    call HDF5_read(P_aim,groupHandle,'P_aim',.false.)
    call MPI_Bcast(P_aim,9,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)
    if(ierr /=0) error stop 'MPI error'
    call HDF5_read(F_aim,groupHandle,'F_aim',.false.)
    call MPI_Bcast(F_aim,9,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)
    if(ierr /=0) error stop 'MPI error'
    call HDF5_read(F_aim_lastInc,groupHandle,'F_aim_lastInc',.false.)
    call MPI_Bcast(F_aim_lastInc,9,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)
    if(ierr /=0) error stop 'MPI error'
    call HDF5_read(F_aimDot,groupHandle,'F_aimDot',.false.)
    call MPI_Bcast(F_aimDot,9,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)
    if(ierr /=0) error stop 'MPI error'
    call HDF5_read(F,groupHandle,'F')
    call HDF5_read(F_lastInc,groupHandle,'F_lastInc')
    call HDF5_read(F_tau,groupHandle,'F_tau')
    call HDF5_read(F_tau_lastInc,groupHandle,'F_tau_lastInc')

  elseif (interface_restartInc == 0) then restartRead
    F_lastInc = spread(spread(spread(math_I3,3,grid(1)),4,grid(2)),5,grid3)                         ! initialize to identity
    F = reshape(F_lastInc,[9,grid(1),grid(2),grid3])
    F_tau = 2.0_pReal*F
    F_tau_lastInc = 2.0_pReal*F_lastInc
  endif restartRead

  homogenization_F0 = reshape(F_lastInc, [3,3,product(grid(1:2))*grid3])                            ! set starting condition for materialpoint_stressAndItsTangent
  call utilities_updateCoords(reshape(F,shape(F_lastInc)))
  call utilities_constitutiveResponse(P,P_av,C_volAvg,C_minMaxAvg, &                                ! stress field, stress avg, global average of stiffness and (min+max)/2
                                      reshape(F,shape(F_lastInc)), &                                ! target F
                                      0.0_pReal)                                                    ! time increment
  call DMDAVecRestoreArrayF90(da,solution_vec,FandF_tau,ierr); CHKERRQ(ierr)                        ! deassociate pointer

  restartRead2: if (interface_restartInc > 0) then
    print'(a,i0,a)', ' reading more restart data of increment ', interface_restartInc, ' from file'
    call HDF5_read(C_volAvg,groupHandle,'C_volAvg',.false.)
    call MPI_Bcast(C_volAvg,81,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)
    if(ierr /=0) error stop 'MPI error'
    call HDF5_read(C_volAvgLastInc,groupHandle,'C_volAvgLastInc',.false.)
    call MPI_Bcast(C_volAvgLastInc,81,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)
    if(ierr /=0) error stop 'MPI error'

    call HDF5_closeGroup(groupHandle)
    call HDF5_closeFile(fileHandle)

    call MPI_File_open(PETSC_COMM_WORLD, trim(getSolverJobName())//'.C_ref', &
                       MPI_MODE_RDONLY,MPI_INFO_NULL,fileUnit,ierr)
    call MPI_File_read(fileUnit,C_minMaxAvg,81,MPI_DOUBLE,MPI_STATUS_IGNORE,ierr)
    call MPI_File_close(fileUnit,ierr)
  endif restartRead2

  call utilities_updateGamma(C_minMaxAvg)
  call utilities_saveReferenceStiffness
  C_scale = C_minMaxAvg
  S_scale = math_invSym3333(C_minMaxAvg)

end subroutine grid_mechanical_spectral_polarisation_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the Polarisation scheme with internal iterations
!--------------------------------------------------------------------------------------------------
function grid_mechanical_spectral_polarisation_solution(incInfoIn) result(solution)

!--------------------------------------------------------------------------------------------------
! input data for solution
  character(len=*), intent(in) :: &
    incInfoIn
  type(tSolutionState)                    :: &
    solution
!--------------------------------------------------------------------------------------------------
! PETSc Data
  PetscErrorCode :: ierr
  SNESConvergedReason :: reason

  incInfo = incInfoIn

!--------------------------------------------------------------------------------------------------
! update stiffness (and gamma operator)
  S = utilities_maskedCompliance(params%rotation_BC,params%stress_mask,C_volAvg)
  if(num%update_gamma) then
    call utilities_updateGamma(C_minMaxAvg)
    C_scale = C_minMaxAvg
    S_scale = math_invSym3333(C_minMaxAvg)
  endif

!--------------------------------------------------------------------------------------------------
! solve BVP
  call SNESsolve(snes,PETSC_NULL_VEC,solution_vec,ierr); CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! check convergence
  call SNESGetConvergedReason(snes,reason,ierr); CHKERRQ(ierr)

  solution%converged = reason > 0
  solution%iterationsNeeded = totalIter
  solution%termIll = terminallyIll
  terminallyIll = .false.
  P_aim = merge(P_aim,P_av,params%stress_mask)

end function grid_mechanical_spectral_polarisation_solution


!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!> @details find new boundary conditions and best F estimate for end of current timestep
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_spectral_polarisation_forward(cutBack,guess,Delta_t,Delta_t_old,t_remaining,&
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
  type(rotation),           intent(in) :: &
    rotation_BC
  PetscErrorCode :: ierr
  PetscScalar, pointer, dimension(:,:,:,:) :: FandF_tau, F, F_tau
  integer :: i, j, k
  real(pReal), dimension(3,3) :: F_lambda33


  call DMDAVecGetArrayF90(da,solution_vec,FandF_tau,ierr); CHKERRQ(ierr)
  F     => FandF_tau(0: 8,:,:,:)
  F_tau => FandF_tau(9:17,:,:,:)

  if (cutBack) then
    C_volAvg    = C_volAvgLastInc
    C_minMaxAvg = C_minMaxAvgLastInc
  else
    C_volAvgLastInc    = C_volAvg
    C_minMaxAvgLastInc = C_minMaxAvg

    F_aimDot = merge(merge((F_aim-F_aim_lastInc)/Delta_t_old,0.0_pReal,stress_BC%mask), 0.0_pReal, guess)  ! estimate deformation rate for prescribed stress components
    F_aim_lastInc = F_aim

    !-----------------------------------------------------------------------------------------------
    ! calculate rate for aim
    if     (deformation_BC%myType=='L') then                                                        ! calculate F_aimDot from given L and current F
      F_aimDot = F_aimDot &
               + merge(matmul(deformation_BC%values, F_aim_lastInc),.0_pReal,deformation_BC%mask)
    elseif (deformation_BC%myType=='dot_F') then                                                    ! F_aimDot is prescribed
      F_aimDot = F_aimDot &
               + merge(deformation_BC%values,.0_pReal,deformation_BC%mask)
    elseif (deformation_BC%myType=='F') then                                                        ! aim at end of load case is prescribed
      F_aimDot = F_aimDot &
               + merge((deformation_BC%values - F_aim_lastInc)/t_remaining,.0_pReal,deformation_BC%mask)
    endif

    Fdot     = utilities_calculateRate(guess, &
                                       F_lastInc,reshape(F,[3,3,grid(1),grid(2),grid3]),Delta_t_old, &
                                       rotation_BC%rotate(F_aimDot,active=.true.))
    F_tauDot = utilities_calculateRate(guess, &
                                       F_tau_lastInc,reshape(F_tau,[3,3,grid(1),grid(2),grid3]), Delta_t_old, &
                                       rotation_BC%rotate(F_aimDot,active=.true.))
    F_lastInc     = reshape(F,    [3,3,grid(1),grid(2),grid3])
    F_tau_lastInc = reshape(F_tau,[3,3,grid(1),grid(2),grid3])

    homogenization_F0 = reshape(F,[3,3,product(grid(1:2))*grid3])
  endif

!--------------------------------------------------------------------------------------------------
! update average and local deformation gradients
  F_aim = F_aim_lastInc + F_aimDot * Delta_t
  if(stress_BC%myType=='P')     P_aim = P_aim &
                                      + merge((stress_BC%values - P_aim)/t_remaining,0.0_pReal,stress_BC%mask)*Delta_t
  if(stress_BC%myType=='dot_P') P_aim = P_aim &
                                      + merge(stress_BC%values,0.0_pReal,stress_BC%mask)*Delta_t

  F = reshape(utilities_forwardField(Delta_t,F_lastInc,Fdot, &                                      ! estimate of F at end of time+Delta_t that matches rotated F_aim on average
                                     rotation_BC%rotate(F_aim,active=.true.)),&
              [9,grid(1),grid(2),grid3])
  if (guess) then
     F_tau = reshape(Utilities_forwardField(Delta_t,F_tau_lastInc,F_taudot), &
                     [9,grid(1),grid(2),grid3])                                                     ! does not have any average value as boundary condition
   else
    do k = 1, grid3; do j = 1, grid(2); do i = 1, grid(1)
       F_lambda33 = reshape(F_tau(1:9,i,j,k)-F(1:9,i,j,k),[3,3])
       F_lambda33 = math_I3 &
                  + math_mul3333xx33(S_scale,0.5_pReal*matmul(F_lambda33, &
                    math_mul3333xx33(C_scale,matmul(transpose(F_lambda33),F_lambda33)-math_I3)))
       F_tau(1:9,i,j,k) = reshape(F_lambda33,[9])+F(1:9,i,j,k)
    enddo; enddo; enddo
  endif

  call DMDAVecRestoreArrayF90(da,solution_vec,FandF_tau,ierr); CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! set module wide available data
  params%stress_mask = stress_BC%mask
  params%rotation_BC = rotation_BC
  params%timeinc     = Delta_t

end subroutine grid_mechanical_spectral_polarisation_forward


!--------------------------------------------------------------------------------------------------
!> @brief Update coordinates
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_spectral_polarisation_updateCoords

  PetscErrorCode :: ierr
  PetscScalar, dimension(:,:,:,:), pointer :: FandF_tau

  call DMDAVecGetArrayF90(da,solution_vec,FandF_tau,ierr); CHKERRQ(ierr)
  call utilities_updateCoords(FandF_tau(0:8,:,:,:))
  call DMDAVecRestoreArrayF90(da,solution_vec,FandF_tau,ierr); CHKERRQ(ierr)

end subroutine grid_mechanical_spectral_polarisation_updateCoords


!--------------------------------------------------------------------------------------------------
!> @brief Write current solver and constitutive data for restart to file
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_spectral_polarisation_restartWrite

  PetscErrorCode :: ierr
  integer(HID_T) :: fileHandle, groupHandle
  PetscScalar, dimension(:,:,:,:), pointer :: FandF_tau, F, F_tau

  call DMDAVecGetArrayF90(da,solution_vec,FandF_tau,ierr); CHKERRQ(ierr)
  F     => FandF_tau(0: 8,:,:,:)
  F_tau => FandF_tau(9:17,:,:,:)

  print*, 'writing solver data required for restart to file'; flush(IO_STDOUT)

  fileHandle  = HDF5_openFile(getSolverJobName()//'_restart.hdf5','w')
  groupHandle = HDF5_addGroup(fileHandle,'solver')

  call HDF5_write(F_aim,groupHandle,'P_aim',.false.)
  call HDF5_write(F_aim,groupHandle,'F_aim',.false.)
  call HDF5_write(F_aim_lastInc,groupHandle,'F_aim_lastInc',.false.)
  call HDF5_write(F_aimDot,groupHandle,'F_aimDot',.false.)
  call HDF5_write(F,groupHandle,'F')
  call HDF5_write(F_lastInc,groupHandle,'F_lastInc')
  call HDF5_write(F_tau,groupHandle,'F_tau')
  call HDF5_write(F_tau_lastInc,groupHandle,'F_tau_lastInc')

  call HDF5_write(C_volAvg,groupHandle,'C_volAvg',.false.)
  call HDF5_write(C_volAvgLastInc,groupHandle,'C_volAvgLastInc',.false.)

  call HDF5_closeGroup(groupHandle)
  call HDF5_closeFile(fileHandle)

  if(num%update_gamma) call utilities_saveReferenceStiffness

  call DMDAVecRestoreArrayF90(da,solution_vec,FandF_tau,ierr); CHKERRQ(ierr)

end subroutine grid_mechanical_spectral_polarisation_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief convergence check
!--------------------------------------------------------------------------------------------------
subroutine converged(snes_local,PETScIter,devNull1,devNull2,devNull3,reason,dummy,ierr)

  SNES :: snes_local
  PetscInt,  intent(in) :: PETScIter
  PetscReal, intent(in) :: &
    devNull1, &
    devNull2, &
    devNull3
  SNESConvergedReason :: reason
  PetscObject :: dummy
  PetscErrorCode :: ierr
  real(pReal) :: &
    curlTol, &
    divTol, &
    BCTol

  curlTol    = max(maxval(abs(F_aim-math_I3))*num%eps_curl_rtol  ,num%eps_curl_atol)
  divTol     = max(maxval(abs(P_av))         *num%eps_div_rtol   ,num%eps_div_atol)
  BCTol      = max(maxval(abs(P_av))         *num%eps_stress_rtol,num%eps_stress_atol)

  if ((totalIter >= num%itmin .and. all([err_div/divTol, err_curl/curlTol, err_BC/BCTol] < 1.0_pReal)) &
       .or. terminallyIll) then
    reason = 1
  elseif (totalIter >= num%itmax) then
    reason = -1
  else
    reason = 0
  endif

!--------------------------------------------------------------------------------------------------
! report
  print'(1/,a)', ' ... reporting .............................................................'
  print'(1/,a,f12.2,a,es8.2,a,es9.2,a)', ' error divergence = ', &
            err_div/divTol,  ' (',err_div, ' / m, tol = ',divTol,')'
  print  '(a,f12.2,a,es8.2,a,es9.2,a)', ' error curl       = ', &
            err_curl/curlTol,' (',err_curl,' -,   tol = ',curlTol,')'
  print  '(a,f12.2,a,es8.2,a,es9.2,a)', ' error stress BC  = ', &
            err_BC/BCTol,    ' (',err_BC,  ' Pa,  tol = ',BCTol,')'
  print'(/,a)', ' ==========================================================================='
  flush(IO_STDOUT)

end subroutine converged


!--------------------------------------------------------------------------------------------------
!> @brief forms the residual vector
!--------------------------------------------------------------------------------------------------
subroutine formResidual(in, FandF_tau, &
                        residuum, dummy,ierr)

  DMDALocalInfo, dimension(DMDA_LOCAL_INFO_SIZE) :: in                                              !< DMDA info (needs to be named "in" for macros like XRANGE to work)
  PetscScalar, dimension(3,3,2,XG_RANGE,YG_RANGE,ZG_RANGE), &
    target, intent(in) :: FandF_tau
  PetscScalar, dimension(3,3,2,X_RANGE,Y_RANGE,Z_RANGE),&
    target,  intent(out) :: residuum                                                                !< residuum field
  PetscScalar, pointer, dimension(:,:,:,:,:) :: &
    F, &
    F_tau, &
    residual_F, &
    residual_F_tau
  PetscInt :: &
    PETScIter, &
    nfuncs
  PetscObject :: dummy
  PetscErrorCode :: ierr
 integer :: &
    i, j, k, e

!---------------------------------------------------------------------------------------------------

  F              => FandF_tau(1:3,1:3,1,&
                                 XG_RANGE,YG_RANGE,ZG_RANGE)
  F_tau          => FandF_tau(1:3,1:3,2,&
                                 XG_RANGE,YG_RANGE,ZG_RANGE)
  residual_F     => residuum(1:3,1:3,1,&
                                 X_RANGE, Y_RANGE, Z_RANGE)
  residual_F_tau => residuum(1:3,1:3,2,&
                                 X_RANGE, Y_RANGE, Z_RANGE)

  F_av = sum(sum(sum(F,dim=5),dim=4),dim=3) * wgt
  call MPI_Allreduce(MPI_IN_PLACE,F_av,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)

  call SNESGetNumberFunctionEvals(snes,nfuncs,ierr); CHKERRQ(ierr)
  call SNESGetIterationNumber(snes,PETScIter,ierr);  CHKERRQ(ierr)

  if (nfuncs == 0 .and. PETScIter == 0) totalIter = -1                                              ! new increment

!--------------------------------------------------------------------------------------------------
! begin of new iteration
  newIteration: if (totalIter <= PETScIter) then
    totalIter = totalIter + 1
    print'(1x,a,3(a,i0))', trim(incInfo), ' @ Iteration ', num%itmin, '≤',totalIter, '≤', num%itmax
    if (debugRotation) print'(/,a,/,2(3(f12.7,1x)/),3(f12.7,1x))', &
      ' deformation gradient aim (lab) =', transpose(params%rotation_BC%rotate(F_aim,active=.true.))
    print'(/,a,/,2(3(f12.7,1x)/),3(f12.7,1x))', &
      ' deformation gradient aim       =', transpose(F_aim)
    flush(IO_STDOUT)
  endif newIteration

!--------------------------------------------------------------------------------------------------
!
  tensorField_real = 0.0_pReal
  do k = 1, grid3; do j = 1, grid(2); do i = 1, grid(1)
    tensorField_real(1:3,1:3,i,j,k) = &
      num%beta*math_mul3333xx33(C_scale,F(1:3,1:3,i,j,k) - math_I3) -&
      num%alpha*matmul(F(1:3,1:3,i,j,k), &
                         math_mul3333xx33(C_scale,F_tau(1:3,1:3,i,j,k) - F(1:3,1:3,i,j,k) - math_I3))
  enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! doing convolution in Fourier space
  call utilities_FFTtensorForward
  call utilities_fourierGammaConvolution(params%rotation_BC%rotate(num%beta*F_aim,active=.true.))
  call utilities_FFTtensorBackward

!--------------------------------------------------------------------------------------------------
! constructing residual
  residual_F_tau = num%beta*F - tensorField_real(1:3,1:3,1:grid(1),1:grid(2),1:grid3)

!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
  call utilities_constitutiveResponse(residual_F, &                                                 ! "residuum" gets field of first PK stress (to save memory)
                                      P_av,C_volAvg,C_minMaxAvg, &
                                      F - residual_F_tau/num%beta,params%timeinc,params%rotation_BC)
  call MPI_Allreduce(MPI_IN_PLACE,terminallyIll,1,MPI_LOGICAL,MPI_LOR,PETSC_COMM_WORLD,ierr)

!--------------------------------------------------------------------------------------------------
! stress BC handling
  F_aim = F_aim - math_mul3333xx33(S, P_av - P_aim)                                                 ! S = 0.0 for no bc
  err_BC = maxval(abs(merge(P_av-P_aim, &
                            math_mul3333xx33(C_scale,F_aim-params%rotation_BC%rotate(F_av)),&
                            params%stress_mask)))
! calculate divergence
  tensorField_real = 0.0_pReal
  tensorField_real(1:3,1:3,1:grid(1),1:grid(2),1:grid3) = residual_F                                !< stress field in disguise
  call utilities_FFTtensorForward
  err_div = utilities_divergenceRMS()                                                               !< root mean squared error in divergence of stress

!--------------------------------------------------------------------------------------------------
! constructing residual
  e = 0
  do k = 1, grid3; do j = 1, grid(2); do i = 1, grid(1)
    e = e + 1
    residual_F(1:3,1:3,i,j,k) = &
      math_mul3333xx33(math_invSym3333(homogenization_dPdF(1:3,1:3,1:3,1:3,e) + C_scale), &
                       residual_F(1:3,1:3,i,j,k) - matmul(F(1:3,1:3,i,j,k), &
                       math_mul3333xx33(C_scale,F_tau(1:3,1:3,i,j,k) - F(1:3,1:3,i,j,k) - math_I3))) &
                       + residual_F_tau(1:3,1:3,i,j,k)
  enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! calculating curl
  tensorField_real = 0.0_pReal
  tensorField_real(1:3,1:3,1:grid(1),1:grid(2),1:grid3) = F
  call utilities_FFTtensorForward
  err_curl = utilities_curlRMS()

end subroutine formResidual

end module grid_mechanical_spectral_polarisation
