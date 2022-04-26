!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Shaokang Zhang, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Spectral solver for thermal conduction
!--------------------------------------------------------------------------------------------------
module grid_thermal_spectral
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

  implicit none
  private

  type :: tNumerics
    integer :: &
      itmax                                                                                         !< maximum number of iterations
    real(pReal) :: &
      eps_thermal_atol, &                                                                           !< absolute tolerance for thermal equilibrium
      eps_thermal_rtol                                                                              !< relative tolerance for thermal equilibrium
  end type tNumerics

  type(tNumerics) :: num

  type(tSolutionParams) :: params
!--------------------------------------------------------------------------------------------------
! PETSc data
  SNES :: SNES_thermal
  Vec :: solution_vec
  real(pReal), dimension(:,:,:), allocatable :: &
    T_current, &                                                                                    !< field of current temperature
    T_lastInc, &                                                                                    !< field of previous temperature
    T_stagInc                                                                                       !< field of staggered temperature

!--------------------------------------------------------------------------------------------------
! reference diffusion tensor, mobility etc.
  integer                     :: totalIter = 0                                                      !< total iteration in current increment
  real(pReal), dimension(3,3) :: K_ref
  real(pReal)                 :: mu_ref

  public :: &
    grid_thermal_spectral_init, &
    grid_thermal_spectral_solution, &
    grid_thermal_spectral_restartWrite, &
    grid_thermal_spectral_forward

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data
!--------------------------------------------------------------------------------------------------
subroutine grid_thermal_spectral_init()

  PetscInt, dimension(0:worldsize-1) :: localK
  integer :: i, j, k, ce
  DM :: thermal_grid
  PetscScalar, dimension(:,:,:), pointer :: T_PETSc
  integer(MPI_INTEGER_KIND) :: err_MPI
  PetscErrorCode :: err_PETSc
  integer(HID_T) :: fileHandle, groupHandle
  class(tNode), pointer :: &
    num_grid

  print'(/,1x,a)', '<<<+-  grid_thermal_spectral init  -+>>>'

  print'(/,1x,a)', 'P. Shanthraj et al., Handbook of Mechanics of Materials, 2019'
  print'(  1x,a)', 'https://doi.org/10.1007/978-981-10-6855-3_80'

!-------------------------------------------------------------------------------------------------
! read numerical parameters and do sanity checks
  num_grid => config_numerics%get('grid',defaultVal=emptyDict)
  num%itmax            = num_grid%get_asInt   ('itmax',           defaultVal=250)
  num%eps_thermal_atol = num_grid%get_asFloat ('eps_thermal_atol',defaultVal=1.0e-2_pReal)
  num%eps_thermal_rtol = num_grid%get_asFloat ('eps_thermal_rtol',defaultVal=1.0e-6_pReal)

  if (num%itmax <= 1)                    call IO_error(301,ext_msg='itmax')
  if (num%eps_thermal_atol <= 0.0_pReal) call IO_error(301,ext_msg='eps_thermal_atol')
  if (num%eps_thermal_rtol <= 0.0_pReal) call IO_error(301,ext_msg='eps_thermal_rtol')

!--------------------------------------------------------------------------------------------------
! set default and user defined options for PETSc
 call PetscOptionsInsertString(PETSC_NULL_OPTIONS,'-thermal_snes_type newtonls -thermal_snes_mf &
                               &-thermal_snes_ksp_ew -thermal_ksp_type fgmres',err_PETSc)
 CHKERRQ(err_PETSc)
 call PetscOptionsInsertString(PETSC_NULL_OPTIONS,num_grid%get_asString('petsc_options',defaultVal=''),err_PETSc)
 CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! init fields
  T_current = discretization_grid_getInitialCondition('T')
  T_lastInc = T_current
  T_stagInc = T_current

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
  call SNESCreate(PETSC_COMM_WORLD,SNES_thermal,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetOptionsPrefix(SNES_thermal,'thermal_',err_PETSc)
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
         1_pPetscInt, 0_pPetscInt, &                                                                ! #dof (T, scalar), ghost boundary width (domain overlap)
         [int(cells(1),pPetscInt)],[int(cells(2),pPetscInt)],localK, &                              ! local cells
         thermal_grid,err_PETSc)                                                                    ! handle, error
  CHKERRQ(err_PETSc)
  call DMsetFromOptions(thermal_grid,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMsetUp(thermal_grid,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(thermal_grid,solution_vec,err_PETSc)                                    ! global solution vector (cells x 1, i.e. every def grad tensor)
  CHKERRQ(err_PETSc)
  call DMDASNESSetFunctionLocal(thermal_grid,INSERT_VALUES,formResidual,PETSC_NULL_SNES,err_PETSc)  ! residual vector of same shape as solution vector
  CHKERRQ(err_PETSc)
  call SNESSetDM(SNES_thermal,thermal_grid,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetFromOptions(SNES_thermal,err_PETSc)                                                   ! pull it all together with additional CLI arguments
  CHKERRQ(err_PETSc)


  restartRead: if (CLI_restartInc > 0) then
    print'(/,1x,a,i0,a)', 'reading restart data of increment ', CLI_restartInc, ' from file'

    fileHandle  = HDF5_openFile(getSolverJobName()//'_restart.hdf5','r')
    groupHandle = HDF5_openGroup(fileHandle,'solver')

    call HDF5_read(T_current,groupHandle,'T',.false.)
    call HDF5_read(T_lastInc,groupHandle,'T_lastInc',.false.)
  end if restartRead

  ce = 0
  do k = 1, cells3; do j = 1, cells(2); do i = 1, cells(1)
    ce = ce + 1
    call homogenization_thermal_setField(T_current(i,j,k),0.0_pReal,ce)
  end do; end do; end do

  call DMDAVecGetArrayF90(thermal_grid,solution_vec,T_PETSc,err_PETSc)
  CHKERRQ(err_PETSc)
  T_PETSc = T_current
  call DMDAVecRestoreArrayF90(thermal_grid,solution_vec,T_PETSc,err_PETSc)
  CHKERRQ(err_PETSc)

  call updateReference()

end subroutine grid_thermal_spectral_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the spectral thermal scheme with internal iterations
!--------------------------------------------------------------------------------------------------
function grid_thermal_spectral_solution(Delta_t) result(solution)

  real(pReal), intent(in) :: &
    Delta_t                                                                                         !< increment in time for current solution
  integer :: i, j, k, ce
  type(tSolutionState) :: solution
  PetscInt  :: devNull
  PetscReal :: T_min, T_max, stagNorm

  integer(MPI_INTEGER_KIND) :: err_MPI
  PetscErrorCode :: err_PETSc
  SNESConvergedReason :: reason

  solution%converged =.false.

!--------------------------------------------------------------------------------------------------
! set module wide availabe data
  params%Delta_t = Delta_t

  call SNESSolve(SNES_thermal,PETSC_NULL_VEC,solution_vec,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESGetConvergedReason(SNES_thermal,reason,err_PETSc)
  CHKERRQ(err_PETSc)

  if (reason < 1) then
    solution%converged = .false.
    solution%iterationsNeeded = num%itmax
  else
    solution%converged = .true.
    solution%iterationsNeeded = totalIter
  end if
  stagNorm = maxval(abs(T_current - T_stagInc))
  call MPI_Allreduce(MPI_IN_PLACE,stagNorm,1_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  solution%stagConverged = stagNorm < max(num%eps_thermal_atol, num%eps_thermal_rtol*maxval(T_current))
  call MPI_Allreduce(MPI_IN_PLACE,solution%stagConverged,1_MPI_INTEGER_KIND,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  T_stagInc = T_current

!--------------------------------------------------------------------------------------------------
! updating thermal state
  ce = 0
  do k = 1, cells3;  do j = 1, cells(2);  do i = 1,cells(1)
    ce = ce + 1
    call homogenization_thermal_setField(T_current(i,j,k),(T_current(i,j,k)-T_lastInc(i,j,k))/params%Delta_t,ce)
  end do; end do; end do

  call VecMin(solution_vec,devNull,T_min,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecMax(solution_vec,devNull,T_max,err_PETSc)
  CHKERRQ(err_PETSc)
  if (solution%converged) &
    print'(/,1x,a)', '... thermal conduction converged ..................................'
  print'(/,1x,a,f8.4,2x,f8.4,2x,f8.4)', 'Minimum|Maximum|Delta Temperature / K = ', T_min, T_max, stagNorm
  print'(/,1x,a)', '==========================================================================='
  flush(IO_STDOUT)

end function grid_thermal_spectral_solution


!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine grid_thermal_spectral_forward(cutBack)

  logical, intent(in) :: cutBack
  integer :: i, j, k, ce
  DM :: dm_local
  PetscScalar,  dimension(:,:,:), pointer :: T_PETSc
  PetscErrorCode :: err_PETSc

  if (cutBack) then
    T_current = T_lastInc
    T_stagInc = T_lastInc

!--------------------------------------------------------------------------------------------------
! reverting thermal field state
    call SNESGetDM(SNES_thermal,dm_local,err_PETSc)
    CHKERRQ(err_PETSc)
    call DMDAVecGetArrayF90(dm_local,solution_vec,T_PETSc,err_PETSc)                                 !< get the data out of PETSc to work with
    CHKERRQ(err_PETSc)
    T_PETSc = T_current
    call DMDAVecRestoreArrayF90(dm_local,solution_vec,T_PETSc,err_PETSc)
    CHKERRQ(err_PETSc)
    ce = 0
    do k = 1, cells3;  do j = 1, cells(2);  do i = 1,cells(1)
      ce = ce + 1
      call homogenization_thermal_setField(T_current(i,j,k),(T_current(i,j,k)-T_lastInc(i,j,k))/params%Delta_t,ce)
    end do; end do; end do
  else
    T_lastInc = T_current
    call updateReference
  end if

end subroutine grid_thermal_spectral_forward


!--------------------------------------------------------------------------------------------------
!> @brief Write current solver and constitutive data for restart to file
!--------------------------------------------------------------------------------------------------
subroutine grid_thermal_spectral_restartWrite

  PetscErrorCode :: err_PETSc
  DM :: dm_local
  integer(HID_T) :: fileHandle, groupHandle
  PetscScalar, dimension(:,:,:), pointer :: T

  call SNESGetDM(SNES_thermal,dm_local,err_PETSc);
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(dm_local,solution_vec,T,err_PETSc);
  CHKERRQ(err_PETSc)

  print'(1x,a)', 'writing thermal solver data required for restart to file'; flush(IO_STDOUT)

  fileHandle  = HDF5_openFile(getSolverJobName()//'_restart.hdf5','a')
  groupHandle = HDF5_openGroup(fileHandle,'solver')
  call HDF5_write(T,groupHandle,'T')
  call HDF5_write(T_lastInc,groupHandle,'T_lastInc')
  call HDF5_closeGroup(groupHandle)
  call HDF5_closeFile(fileHandle)

  call DMDAVecRestoreArrayF90(dm_local,solution_vec,T,err_PETSc);
  CHKERRQ(err_PETSc)

end subroutine grid_thermal_spectral_restartWrite



!--------------------------------------------------------------------------------------------------
!> @brief forms the spectral thermal residual vector
!--------------------------------------------------------------------------------------------------
subroutine formResidual(in,x_scal,r,dummy,err_PETSc)

  DMDALocalInfo, dimension(DMDA_LOCAL_INFO_SIZE) :: &
    in
  PetscScalar, dimension( &
    XG_RANGE,YG_RANGE,ZG_RANGE), intent(in) :: &
    x_scal
  PetscScalar, dimension( &
    X_RANGE,Y_RANGE,Z_RANGE), intent(out) :: &
    r
  PetscObject :: dummy
  PetscErrorCode :: err_PETSc
  integer :: i, j, k, ce

  T_current = x_scal
!--------------------------------------------------------------------------------------------------
! evaluate polarization field
  scalarField_real = 0.0_pReal
  scalarField_real(1:cells(1),1:cells(2),1:cells3) = T_current
  call utilities_FFTscalarForward
  call utilities_fourierScalarGradient                                                              !< calculate gradient of temperature field
  call utilities_FFTvectorBackward
  ce = 0
  do k = 1, cells3;  do j = 1, cells(2);  do i = 1,cells(1)
    ce = ce + 1
    vectorField_real(1:3,i,j,k) = matmul(homogenization_K_T(ce) - K_ref, vectorField_real(1:3,i,j,k))
  end do; end do; end do
  call utilities_FFTvectorForward
  call utilities_fourierVectorDivergence                                                            !< calculate temperature divergence in fourier field
  call utilities_FFTscalarBackward
  ce = 0
  do k = 1, cells3;  do j = 1, cells(2);  do i = 1,cells(1)
    ce = ce + 1
    scalarField_real(i,j,k) = params%Delta_t*(scalarField_real(i,j,k) + homogenization_f_T(ce)) &
                            + homogenization_mu_T(ce) * (T_lastInc(i,j,k) - T_current(i,j,k)) &
                            + mu_ref*T_current(i,j,k)
  end do; end do; end do

!--------------------------------------------------------------------------------------------------
! convolution of temperature field with green operator
  call utilities_FFTscalarForward
  call utilities_fourierGreenConvolution(K_ref, mu_ref, params%Delta_t)
  call utilities_FFTscalarBackward

!--------------------------------------------------------------------------------------------------
! constructing residual
  r = T_current - scalarField_real(1:cells(1),1:cells(2),1:cells3)
  err_PETSc = 0

end subroutine formResidual


!--------------------------------------------------------------------------------------------------
!> @brief update reference viscosity and conductivity
!--------------------------------------------------------------------------------------------------
subroutine updateReference()

  integer :: ce
  integer(MPI_INTEGER_KIND) :: err_MPI


  K_ref = 0.0_pReal
  mu_ref = 0.0_pReal
  do ce = 1, product(cells(1:2))*cells3
    K_ref  = K_ref  + homogenization_K_T(ce)
    mu_ref = mu_ref + homogenization_mu_T(ce)
  end do

  K_ref = K_ref*wgt
  call MPI_Allreduce(MPI_IN_PLACE,K_ref,9_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  mu_ref = mu_ref*wgt
  call MPI_Allreduce(MPI_IN_PLACE,mu_ref,1_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'

end subroutine updateReference


end module grid_thermal_spectral
