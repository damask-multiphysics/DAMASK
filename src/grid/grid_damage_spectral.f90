!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Shaokang Zhang, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Spectral solver for nonlocal damage
!--------------------------------------------------------------------------------------------------
module grid_damage_spectral
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
  use CLI
  use HDF5_utilities
  use HDF5
  use spectral_utilities
  use discretization_grid
  use homogenization
  use YAML_types
  use config

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  implicit none(type,external)
#else
  implicit none
#endif
  private

  type :: tNumerics
    integer :: &
      itmax                                                                                         !< maximum number of iterations
    real(pReal) :: &
      phi_min, &                                                                                    !< non-zero residual damage
      eps_damage_atol, &                                                                            !< absolute tolerance for damage evolution
      eps_damage_rtol                                                                               !< relative tolerance for damage evolution
  end type tNumerics

  type(tNumerics) :: num

  type(tSolutionParams) :: params
!--------------------------------------------------------------------------------------------------
! PETSc data
  SNES :: SNES_damage
  Vec  :: solution_vec
  real(pReal), dimension(:,:,:), allocatable :: &
    phi, &                                                                                          !< field of current damage
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
    grid_damage_spectral_restartWrite, &
    grid_damage_spectral_forward

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data
!--------------------------------------------------------------------------------------------------
subroutine grid_damage_spectral_init()

  PetscInt, dimension(0:worldsize-1) :: localK
  integer :: i, j, k, ce
  DM :: damage_grid
  real(pReal), dimension(:,:,:), pointer :: phi_PETSc
  Vec :: uBound, lBound
  integer(MPI_INTEGER_KIND) :: err_MPI
  PetscErrorCode :: err_PETSc
  integer(HID_T) :: fileHandle, groupHandle
  real(pReal), dimension(1,product(cells(1:2))*cells3) :: tempN
  type(tDict), pointer :: &
    num_grid, &
    num_generic
  character(len=pStringLen) :: &
    snes_type

  print'(/,1x,a)', '<<<+-  grid_spectral_damage init  -+>>>'

  print'(/,1x,a)', 'P. Shanthraj et al., Handbook of Mechanics of Materials, 2019'
  print'(  1x,a)', 'https://doi.org/10.1007/978-981-10-6855-3_80'

!-------------------------------------------------------------------------------------------------
! read numerical parameters and do sanity checks
  num_grid => config_numerics%get_dict('grid',defaultVal=emptyDict)
  num%itmax           = num_grid%get_asInt   ('itmax',defaultVal=250)
  num%eps_damage_atol = num_grid%get_asFloat ('eps_damage_atol',defaultVal=1.0e-2_pReal)
  num%eps_damage_rtol = num_grid%get_asFloat ('eps_damage_rtol',defaultVal=1.0e-6_pReal)

  num_generic => config_numerics%get_dict('generic',defaultVal=emptyDict)
  num%phi_min = num_generic%get_asFloat('phi_min', defaultVal=1.0e-6_pReal)

  if (num%phi_min < 0.0_pReal) call IO_error(301,ext_msg='phi_min')
  if (num%itmax <= 1)                    call IO_error(301,ext_msg='itmax')
  if (num%eps_damage_atol <= 0.0_pReal)  call IO_error(301,ext_msg='eps_damage_atol')
  if (num%eps_damage_rtol <= 0.0_pReal)  call IO_error(301,ext_msg='eps_damage_rtol')

!--------------------------------------------------------------------------------------------------
! set default and user defined options for PETSc
 call PetscOptionsInsertString(PETSC_NULL_OPTIONS,'-damage_snes_type newtonls -damage_snes_mf &
                               &-damage_snes_ksp_ew -damage_ksp_type fgmres',err_PETSc)
 CHKERRQ(err_PETSc)
 call PetscOptionsInsertString(PETSC_NULL_OPTIONS,num_grid%get_asString('petsc_options',defaultVal=''),err_PETSc)
 CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! init fields
  phi = discretization_grid_getInitialCondition('phi')
  phi_lastInc = phi
  phi_stagInc = phi

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
  call SNESCreate(PETSC_COMM_WORLD,SNES_damage,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetOptionsPrefix(SNES_damage,'damage_',err_PETSc)
  CHKERRQ(err_PETSc)
  localK            = 0_pPetscInt
  localK(worldrank) = int(cells3,pPetscInt)
  call MPI_Allreduce(MPI_IN_PLACE,localK,worldsize,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  call DMDACreate3D(PETSC_COMM_WORLD, &
         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &                                    ! cut off stencil at boundary
         DMDA_STENCIL_BOX, &                                                                        ! Moore (26) neighborhood around central point
         int(cells(1),pPetscInt),int(cells(2),pPetscInt),int(cells(3),pPetscInt), &                 ! global cells
         1_pPetscInt, 1_pPetscInt, int(worldsize,pPetscInt), &
         1_pPetscInt, 0_pPetscInt, &                                                                ! #dof (phi, scalar), ghost boundary width (domain overlap)
         [int(cells(1),pPetscInt)],[int(cells(2),pPetscInt)],localK, &                              ! local cells
         damage_grid,err_PETSc)                                                                     ! handle, error
  CHKERRQ(err_PETSc)
  call DMsetFromOptions(damage_grid,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMsetUp(damage_grid,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(damage_grid,solution_vec,err_PETSc)                                     ! global solution vector (cells x 1, i.e. every def grad tensor)
  CHKERRQ(err_PETSc)
  call DMDASNESSetFunctionLocal(damage_grid,INSERT_VALUES,formResidual,PETSC_NULL_SNES,err_PETSc)   ! residual vector of same shape as solution vector
  CHKERRQ(err_PETSc)
  call SNESSetDM(SNES_damage,damage_grid,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetFromOptions(SNES_damage,err_PETSc)                                                    ! pull it all together with additional CLI arguments
  CHKERRQ(err_PETSc)
  call SNESGetType(SNES_damage,snes_type,err_PETSc)
  CHKERRQ(err_PETSc)
  if (trim(snes_type) == 'vinewtonrsls' .or. &
      trim(snes_type) == 'vinewtonssls') then
    call DMGetGlobalVector(damage_grid,lBound,err_PETSc)
    CHKERRQ(err_PETSc)
    call DMGetGlobalVector(damage_grid,uBound,err_PETSc)
    CHKERRQ(err_PETSc)
    call VecSet(lBound,0.0_pReal,err_PETSc)
    CHKERRQ(err_PETSc)
    call VecSet(uBound,1.0_pReal,err_PETSc)
    CHKERRQ(err_PETSc)
    call SNESVISetVariableBounds(SNES_damage,lBound,uBound,err_PETSc)                               ! variable bounds for variational inequalities
    CHKERRQ(err_PETSc)
    call DMRestoreGlobalVector(damage_grid,lBound,err_PETSc)
    CHKERRQ(err_PETSc)
    call DMRestoreGlobalVector(damage_grid,uBound,err_PETSc)
    CHKERRQ(err_PETSc)
  end if

  restartRead: if (CLI_restartInc > 0) then
    print'(/,1x,a,i0,a)', 'reading restart data of increment ', CLI_restartInc, ' from file'

    fileHandle  = HDF5_openFile(getSolverJobName()//'_restart.hdf5','r')
    groupHandle = HDF5_openGroup(fileHandle,'solver')

    call HDF5_read(tempN,groupHandle,'phi',.false.)
    phi = reshape(tempN,[cells(1),cells(2),cells3])
    call HDF5_read(tempN,groupHandle,'phi_lastInc',.false.)
    phi_lastInc = reshape(tempN,[cells(1),cells(2),cells3])
  end if restartRead

  ce = 0
  do k = 1, cells3; do j = 1, cells(2); do i = 1, cells(1)
    ce = ce + 1
    call homogenization_set_phi(phi(i,j,k),ce)
  end do; end do; end do

  call DMDAVecGetArrayF90(damage_grid,solution_vec,phi_PETSc,err_PETSc)
  CHKERRQ(err_PETSc)
  phi_PETSc = phi
  call DMDAVecRestoreArrayF90(damage_grid,solution_vec,phi_PETSc,err_PETSc)
  CHKERRQ(err_PETSc)

  call updateReference()

end subroutine grid_damage_spectral_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the spectral damage scheme with internal iterations
!--------------------------------------------------------------------------------------------------
function grid_damage_spectral_solution(Delta_t) result(solution)

  real(pReal), intent(in) :: &
    Delta_t                                                                                         !< increment in time for current solution
  integer :: i, j, k, ce
  type(tSolutionState) :: solution
  PetscInt  :: devNull
  PetscReal :: phi_min, phi_max, stagNorm

  integer(MPI_INTEGER_KIND) :: err_MPI
  PetscErrorCode :: err_PETSc
  SNESConvergedReason :: reason

  solution%converged =.false.

!--------------------------------------------------------------------------------------------------
! set module wide availabe data
  params%Delta_t = Delta_t

  call SNESSolve(SNES_damage,PETSC_NULL_VEC,solution_vec,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESGetConvergedReason(SNES_damage,reason,err_PETSc)
  CHKERRQ(err_PETSc)

  if (reason < 1) then
    solution%converged = .false.
    solution%iterationsNeeded = num%itmax
  else
    solution%converged = .true.
    solution%iterationsNeeded = totalIter
  end if
  stagNorm = maxval(abs(phi - phi_stagInc))
  call MPI_Allreduce(MPI_IN_PLACE,stagNorm,1_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  solution%stagConverged = stagNorm < max(num%eps_damage_atol, num%eps_damage_rtol*maxval(phi))
  call MPI_Allreduce(MPI_IN_PLACE,solution%stagConverged,1_MPI_INTEGER_KIND,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  phi_stagInc = phi

!--------------------------------------------------------------------------------------------------
! updating damage state
  ce = 0
  do k = 1, cells3;  do j = 1, cells(2);  do i = 1,cells(1)
    ce = ce + 1
    call homogenization_set_phi(phi(i,j,k),ce)
  end do; end do; end do

  call VecMin(solution_vec,devNull,phi_min,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecMax(solution_vec,devNull,phi_max,err_PETSc)
  CHKERRQ(err_PETSc)
  if (solution%converged) &
    print'(/,1x,a)', '... nonlocal damage converged .....................................'
  print'(/,1x,a,f8.6,2x,f8.6,2x,e11.4)', 'Minimum|Maximum|Delta Damage      = ', phi_min, phi_max, stagNorm
  print'(/,1x,a)', '==========================================================================='
  flush(IO_STDOUT)

end function grid_damage_spectral_solution


!--------------------------------------------------------------------------------------------------
!> @brief spectral damage forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine grid_damage_spectral_forward(cutBack)

  logical, intent(in) :: cutBack

  integer :: i, j, k, ce
  DM :: dm_local
  real(pReal),  dimension(:,:,:), pointer :: phi_PETSc
  PetscErrorCode :: err_PETSc


  if (cutBack) then
    phi = phi_lastInc
    phi_stagInc = phi_lastInc
!--------------------------------------------------------------------------------------------------
! reverting damage field state
    call SNESGetDM(SNES_damage,dm_local,err_PETSc)
    CHKERRQ(err_PETSc)
    call DMDAVecGetArrayF90(dm_local,solution_vec,phi_PETSc,err_PETSc)                              !< get the data out of PETSc to work with
    CHKERRQ(err_PETSc)
    phi_PETSc = phi
    call DMDAVecRestoreArrayF90(dm_local,solution_vec,phi_PETSc,err_PETSc)
    CHKERRQ(err_PETSc)
    ce = 0
    do k = 1, cells3;  do j = 1, cells(2);  do i = 1,cells(1)
      ce = ce + 1
      call homogenization_set_phi(phi(i,j,k),ce)
    end do; end do; end do
  else
    phi_lastInc = phi
    call updateReference
  end if

end subroutine grid_damage_spectral_forward


!--------------------------------------------------------------------------------------------------
!> @brief Write current solver and constitutive data for restart to file.
!--------------------------------------------------------------------------------------------------
subroutine grid_damage_spectral_restartWrite

  PetscErrorCode :: err_PETSc
  DM :: dm_local
  integer(HID_T) :: fileHandle, groupHandle
  PetscScalar, dimension(:,:,:), pointer :: phi

  call SNESGetDM(SNES_damage,dm_local,err_PETSc);
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(dm_local,solution_vec,phi,err_PETSc);
  CHKERRQ(err_PETSc)

  print'(1x,a)', 'writing damage solver data required for restart to file'; flush(IO_STDOUT)

  fileHandle  = HDF5_openFile(getSolverJobName()//'_restart.hdf5','a')
  groupHandle = HDF5_openGroup(fileHandle,'solver')
  call HDF5_write(reshape(phi,[1,product(cells(1:2))*cells3]),groupHandle,'phi')
  call HDF5_write(reshape(phi_lastInc,[1,product(cells(1:2))*cells3]),groupHandle,'phi_lastInc')
  call HDF5_closeGroup(groupHandle)
  call HDF5_closeFile(fileHandle)

  call DMDAVecRestoreArrayF90(dm_local,solution_vec,phi,err_PETSc);
  CHKERRQ(err_PETSc)

end subroutine grid_damage_spectral_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief Construct the residual vector.
!--------------------------------------------------------------------------------------------------
subroutine formResidual(residual_subdomain,x_scal,r,dummy,err_PETSc)

  DMDALocalInfo, dimension(DMDA_LOCAL_INFO_SIZE) :: &
    residual_subdomain
  real(pReal), dimension(cells(1),cells(2),cells3), intent(in) :: &
    x_scal
  real(pReal), dimension(cells(1),cells(2),cells3), intent(out) :: &
    r                                                                                               !< residual
  PetscObject :: dummy
  PetscErrorCode, intent(out) :: err_PETSc

  integer :: i, j, k, ce
  real(pReal), dimension(3,cells(1),cells(2),cells3) :: vectorField


  phi = x_scal
  vectorField = utilities_ScalarGradient(phi)
  ce = 0
  do k = 1, cells3;  do j = 1, cells(2);  do i = 1,cells(1)
    ce = ce + 1
    vectorField(1:3,i,j,k) = matmul(homogenization_K_phi(ce) - K_ref, vectorField(1:3,i,j,k))
  end do; end do; end do
  r = utilities_VectorDivergence(vectorField)
  ce = 0
  do k = 1, cells3;  do j = 1, cells(2);  do i = 1,cells(1)
    ce = ce + 1
    r(i,j,k) = params%Delta_t*(r(i,j,k) + homogenization_f_phi(phi(i,j,k),ce)) &
             + homogenization_mu_phi(ce)*(phi_lastInc(i,j,k) - phi(i,j,k)) &
             + mu_ref*phi(i,j,k)
  end do; end do; end do

  r = max(min(utilities_GreenConvolution(r, K_ref, mu_ref, params%Delta_t),phi_lastInc),num%phi_min) &
    - phi
  err_PETSc = 0

end subroutine formResidual


!--------------------------------------------------------------------------------------------------
!> @brief Update reference viscosity and conductivity.
!--------------------------------------------------------------------------------------------------
subroutine updateReference()

  integer :: ce
  integer(MPI_INTEGER_KIND) :: err_MPI


  K_ref = 0.0_pReal
  mu_ref = 0.0_pReal
  do ce = 1, product(cells(1:2))*cells3
    K_ref  = K_ref  + homogenization_K_phi(ce)
    mu_ref = mu_ref + homogenization_mu_phi(ce)
  end do

  K_ref = K_ref*wgt
  call MPI_Allreduce(MPI_IN_PLACE,K_ref,9_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  mu_ref = mu_ref*wgt
  call MPI_Allreduce(MPI_IN_PLACE,mu_ref,1_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'

end subroutine updateReference


end module grid_damage_spectral
