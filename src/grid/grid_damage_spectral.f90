!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Shaokang Zhang, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Spectral solver for nonlocal damage
!--------------------------------------------------------------------------------------------------
module grid_damage_spectral
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscdmda.h>
  use PETScdmda
  use PETScsnes

  use prec
  use parallelization
  use IO
  use spectral_utilities
  use discretization_grid
  use homogenization
  use YAML_types
  use config

  implicit none
  private

  type :: tNumerics
    integer :: &
      itmax                                                                                         !< maximum number of iterations
    real(pReal) :: &
      residualStiffness, &                                                                          !< non-zero residual damage
      eps_damage_atol, &                                                                            !< absolute tolerance for damage evolution
      eps_damage_rtol                                                                               !< relative tolerance for damage evolution
  end type tNumerics

  type(tNumerics) :: num

  type(tSolutionParams) :: params
!--------------------------------------------------------------------------------------------------
! PETSc data
  SNES :: damage_snes
  Vec  :: solution_vec
  PetscInt :: xstart, xend, ystart, yend, zstart, zend
  real(pReal), dimension(:,:,:), allocatable :: &
    phi_current, &                                                                                  !< field of current damage
    phi_lastInc, &                                                                                  !< field of previous damage
    phi_stagInc                                                                                     !< field of staggered damage

!--------------------------------------------------------------------------------------------------
! reference diffusion tensor, mobility etc.
  integer                     :: totalIter = 0                                                      !< total iteration in current increment
  real(pReal), dimension(3,3) :: K_ref
  real(pReal)                 :: mu_ref

  public :: &
    grid_damage_spectral_init, &
    grid_damage_spectral_solution, &
    grid_damage_spectral_forward

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data
! ToDo: Restart not implemented
!--------------------------------------------------------------------------------------------------
subroutine grid_damage_spectral_init()

  PetscInt, dimension(0:worldsize-1) :: localK
  DM :: damage_grid
  Vec :: uBound, lBound
  PetscErrorCode :: ierr
  class(tNode), pointer :: &
    num_grid, &
    num_generic
  character(len=pStringLen) :: &
    snes_type

  print'(/,a)', ' <<<+-  grid_spectral_damage init  -+>>>'

  print*, 'P. Shanthraj et al., Handbook of Mechanics of Materials, 2019'
  print*, 'https://doi.org/10.1007/978-981-10-6855-3_80'

!-------------------------------------------------------------------------------------------------
! read numerical parameters and do sanity checks
  num_grid => config_numerics%get('grid',defaultVal=emptyDict)
  num%itmax           = num_grid%get_asInt   ('itmax',defaultVal=250)
  num%eps_damage_atol = num_grid%get_asFloat ('eps_damage_atol',defaultVal=1.0e-2_pReal)
  num%eps_damage_rtol = num_grid%get_asFloat ('eps_damage_rtol',defaultVal=1.0e-6_pReal)

  num_generic => config_numerics%get('generic',defaultVal=emptyDict)
  num%residualStiffness = num_generic%get_asFloat('residualStiffness', defaultVal=1.0e-6_pReal)

  if (num%residualStiffness < 0.0_pReal) call IO_error(301,ext_msg='residualStiffness')
  if (num%itmax <= 1)                    call IO_error(301,ext_msg='itmax')
  if (num%eps_damage_atol <= 0.0_pReal)  call IO_error(301,ext_msg='eps_damage_atol')
  if (num%eps_damage_rtol <= 0.0_pReal)  call IO_error(301,ext_msg='eps_damage_rtol')

!--------------------------------------------------------------------------------------------------
! set default and user defined options for PETSc
 call PetscOptionsInsertString(PETSC_NULL_OPTIONS,'-damage_snes_type newtonls -damage_snes_mf &
                               &-damage_snes_ksp_ew -damage_ksp_type fgmres',ierr)
 CHKERRQ(ierr)
 call PetscOptionsInsertString(PETSC_NULL_OPTIONS,num_grid%get_asString('petsc_options',defaultVal=''),ierr)
 CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
  call SNESCreate(PETSC_COMM_WORLD,damage_snes,ierr); CHKERRQ(ierr)
  call SNESSetOptionsPrefix(damage_snes,'damage_',ierr);CHKERRQ(ierr)
  localK            = 0
  localK(worldrank) = grid3
  call MPI_Allreduce(MPI_IN_PLACE,localK,worldsize,MPI_INTEGER,MPI_SUM,PETSC_COMM_WORLD,ierr)
  call DMDACreate3D(PETSC_COMM_WORLD, &
         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &                                    ! cut off stencil at boundary
         DMDA_STENCIL_BOX, &                                                                        ! Moore (26) neighborhood around central point
         grid(1),grid(2),grid(3), &                                                                 ! global grid
         1, 1, worldsize, &
         1, 0, &                                                                                    ! #dof (damage phase field), ghost boundary width (domain overlap)
         [grid(1)],[grid(2)],localK, &                                                              ! local grid
         damage_grid,ierr)                                                                          ! handle, error
  CHKERRQ(ierr)
  call SNESSetDM(damage_snes,damage_grid,ierr); CHKERRQ(ierr)                                       ! connect snes to da
  call DMsetFromOptions(damage_grid,ierr); CHKERRQ(ierr)
  call DMsetUp(damage_grid,ierr); CHKERRQ(ierr)
  call DMCreateGlobalVector(damage_grid,solution_vec,ierr); CHKERRQ(ierr)                           ! global solution vector (grid x 1, i.e. every def grad tensor)
  call DMDASNESSetFunctionLocal(damage_grid,INSERT_VALUES,formResidual,PETSC_NULL_SNES,ierr)        ! residual vector of same shape as solution vector
  CHKERRQ(ierr)
  call SNESSetFromOptions(damage_snes,ierr); CHKERRQ(ierr)                                          ! pull it all together with additional CLI arguments
  call SNESGetType(damage_snes,snes_type,ierr); CHKERRQ(ierr)
  if (trim(snes_type) == 'vinewtonrsls' .or. &
      trim(snes_type) == 'vinewtonssls') then
    call DMGetGlobalVector(damage_grid,lBound,ierr); CHKERRQ(ierr)
    call DMGetGlobalVector(damage_grid,uBound,ierr); CHKERRQ(ierr)
    call VecSet(lBound,0.0_pReal,ierr); CHKERRQ(ierr)
    call VecSet(uBound,1.0_pReal,ierr); CHKERRQ(ierr)
    call SNESVISetVariableBounds(damage_snes,lBound,uBound,ierr)                                    ! variable bounds for variational inequalities like contact mechanics, damage etc.
    call DMRestoreGlobalVector(damage_grid,lBound,ierr); CHKERRQ(ierr)
    call DMRestoreGlobalVector(damage_grid,uBound,ierr); CHKERRQ(ierr)
  endif

!--------------------------------------------------------------------------------------------------
! init fields
  call DMDAGetCorners(damage_grid,xstart,ystart,zstart,xend,yend,zend,ierr)
  CHKERRQ(ierr)
  xend = xstart + xend - 1
  yend = ystart + yend - 1
  zend = zstart + zend - 1
  allocate(phi_current(grid(1),grid(2),grid3), source=1.0_pReal)
  allocate(phi_lastInc(grid(1),grid(2),grid3), source=1.0_pReal)
  allocate(phi_stagInc(grid(1),grid(2),grid3), source=1.0_pReal)

  call VecSet(solution_vec,1.0_pReal,ierr); CHKERRQ(ierr)

  call updateReference

end subroutine grid_damage_spectral_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the spectral damage scheme with internal iterations
!--------------------------------------------------------------------------------------------------
function grid_damage_spectral_solution(timeinc) result(solution)

  real(pReal), intent(in) :: &
    timeinc                                                                                         !< increment in time for current solution
  integer :: i, j, k, ce
  type(tSolutionState) :: solution
  PetscInt  :: devNull
  PetscReal :: phi_min, phi_max, stagNorm, solnNorm

  PetscErrorCode :: ierr
  SNESConvergedReason :: reason

  solution%converged =.false.

!--------------------------------------------------------------------------------------------------
! set module wide availabe data
  params%timeinc = timeinc

  call SNESSolve(damage_snes,PETSC_NULL_VEC,solution_vec,ierr); CHKERRQ(ierr)
  call SNESGetConvergedReason(damage_snes,reason,ierr); CHKERRQ(ierr)

  if (reason < 1) then
    solution%converged = .false.
    solution%iterationsNeeded = num%itmax
  else
    solution%converged = .true.
    solution%iterationsNeeded = totalIter
  endif
  stagNorm = maxval(abs(phi_current - phi_stagInc))
  solnNorm = maxval(abs(phi_current))
  call MPI_Allreduce(MPI_IN_PLACE,stagNorm,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
  call MPI_Allreduce(MPI_IN_PLACE,solnNorm,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
  phi_stagInc = phi_current
  solution%stagConverged = stagNorm < max(num%eps_damage_atol, num%eps_damage_rtol*solnNorm)

!--------------------------------------------------------------------------------------------------
! updating damage state
  ce = 0
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    ce = ce + 1
    call homogenization_set_phi(phi_current(i,j,k),ce)
  enddo; enddo; enddo

  call VecMin(solution_vec,devNull,phi_min,ierr); CHKERRQ(ierr)
  call VecMax(solution_vec,devNull,phi_max,ierr); CHKERRQ(ierr)
  if (solution%converged) &
    print'(/,a)', ' ... nonlocal damage converged .....................................'
  print'(/,a,f8.6,2x,f8.6,2x,e11.4)', ' Minimum|Maximum|Delta Damage      = ', phi_min, phi_max, stagNorm
  print'(/,a)', ' ==========================================================================='
  flush(IO_STDOUT)

end function grid_damage_spectral_solution


!--------------------------------------------------------------------------------------------------
!> @brief spectral damage forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine grid_damage_spectral_forward(cutBack)

  logical, intent(in) :: cutBack
  integer :: i, j, k, ce
  DM :: dm_local
  PetscScalar,  dimension(:,:,:), pointer :: x_scal
  PetscErrorCode :: ierr

  if (cutBack) then
    phi_current = phi_lastInc
    phi_stagInc = phi_lastInc
!--------------------------------------------------------------------------------------------------
! reverting damage field state
    call SNESGetDM(damage_snes,dm_local,ierr); CHKERRQ(ierr)
    call DMDAVecGetArrayF90(dm_local,solution_vec,x_scal,ierr); CHKERRQ(ierr)                       !< get the data out of PETSc to work with
    x_scal(xstart:xend,ystart:yend,zstart:zend) = phi_current
    call DMDAVecRestoreArrayF90(dm_local,solution_vec,x_scal,ierr); CHKERRQ(ierr)
    ce = 0
    do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
      ce = ce + 1
      call homogenization_set_phi(phi_current(i,j,k),ce)
    enddo; enddo; enddo
  else
    phi_lastInc = phi_current
    call updateReference
  endif

end subroutine grid_damage_spectral_forward


!--------------------------------------------------------------------------------------------------
!> @brief forms the spectral damage residual vector
!--------------------------------------------------------------------------------------------------
subroutine formResidual(in,x_scal,f_scal,dummy,ierr)

  DMDALocalInfo, dimension(DMDA_LOCAL_INFO_SIZE) :: &
    in
  PetscScalar, dimension( &
    XG_RANGE,YG_RANGE,ZG_RANGE), intent(in) :: &
    x_scal
  PetscScalar, dimension( &
    X_RANGE,Y_RANGE,Z_RANGE), intent(out) :: &
    f_scal
  PetscObject :: dummy
  PetscErrorCode :: ierr
  integer :: i, j, k, ce

  phi_current = x_scal
!--------------------------------------------------------------------------------------------------
! evaluate polarization field
  scalarField_real = 0.0_pReal
  scalarField_real(1:grid(1),1:grid(2),1:grid3) = phi_current
  call utilities_FFTscalarForward
  call utilities_fourierScalarGradient                                                              !< calculate gradient of damage field
  call utilities_FFTvectorBackward
  ce = 0
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    ce = ce + 1
    vectorField_real(1:3,i,j,k) = matmul(homogenization_K_phi(ce) - K_ref, vectorField_real(1:3,i,j,k))
  enddo; enddo; enddo
  call utilities_FFTvectorForward
  call utilities_fourierVectorDivergence                                                            !< calculate damage divergence in fourier field
  call utilities_FFTscalarBackward
  ce = 0
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    ce = ce + 1
    scalarField_real(i,j,k) = params%timeinc*(scalarField_real(i,j,k) + homogenization_f_phi(phi_current(i,j,k),ce)) &
                            + homogenization_mu_phi(ce)*(phi_lastInc(i,j,k) - phi_current(i,j,k)) &
                            + mu_ref*phi_current(i,j,k)
  enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! convolution of damage field with green operator
  call utilities_FFTscalarForward
  call utilities_fourierGreenConvolution(K_ref, mu_ref, params%timeinc)
  call utilities_FFTscalarBackward

  where(scalarField_real(1:grid(1),1:grid(2),1:grid3) > phi_lastInc) &
        scalarField_real(1:grid(1),1:grid(2),1:grid3) = phi_lastInc
  where(scalarField_real(1:grid(1),1:grid(2),1:grid3) < num%residualStiffness) &
        scalarField_real(1:grid(1),1:grid(2),1:grid3) = num%residualStiffness

!--------------------------------------------------------------------------------------------------
! constructing residual
  f_scal = scalarField_real(1:grid(1),1:grid(2),1:grid3) - phi_current

end subroutine formResidual


!--------------------------------------------------------------------------------------------------
!> @brief update reference viscosity and conductivity
!--------------------------------------------------------------------------------------------------
subroutine updateReference()

  integer :: ce,ierr


  K_ref = 0.0_pReal
  mu_ref = 0.0_pReal
  do ce = 1, product(grid(1:2))*grid3
    K_ref  = K_ref  + homogenization_K_phi(ce)
    mu_ref = mu_ref + homogenization_mu_phi(ce)
  enddo

  K_ref = K_ref*wgt
  call MPI_Allreduce(MPI_IN_PLACE,K_ref,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
  mu_ref = mu_ref*wgt
  call MPI_Allreduce(MPI_IN_PLACE,mu_ref,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)

end subroutine updateReference


end module grid_damage_spectral
