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
  use misc
  use CLI
  use HDF5_utilities
  use HDF5
  use spectral_utilities
  use discretization_grid
  use homogenization
  use types
  use config
  use constants

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  implicit none(type,external)
#else
  implicit none
#endif
  private

  type :: tNumerics
    integer :: &
      itmax                                                                                         !< maximum number of iterations
    real(pREAL) :: &
      eps_thermal_atol, &                                                                           !< absolute tolerance for thermal equilibrium
      eps_thermal_rtol                                                                              !< relative tolerance for thermal equilibrium
  end type tNumerics

  type(tNumerics) :: num

!--------------------------------------------------------------------------------------------------
! PETSc data
  SNES :: SNES_thermal
  Vec :: T_PETSc
  real(pREAL), dimension(:,:,:), allocatable :: &
    T_lastInc, &                                                                                    !< field of previous temperature
    T_stagInc, &                                                                                    !< field of staggered temperature
    dotT_lastInc
!--------------------------------------------------------------------------------------------------
! reference diffusion tensor, mobility etc.
  integer                     :: totalIter = 0                                                      !< total iteration in current increment
  real(pREAL), dimension(3,3) :: K_ref
  real(pREAL)                 :: mu_ref, Delta_t_
  integer(kind(STATUS_OK))    :: status

  public :: &
    grid_thermal_spectral_init, &
    grid_thermal_spectral_solution, &
    grid_thermal_spectral_restartWrite, &
    grid_thermal_spectral_forward

contains

!--------------------------------------------------------------------------------------------------
!> @brief Allocate all necessary fields and fill them with data, potentially from restart info.
!--------------------------------------------------------------------------------------------------
subroutine grid_thermal_spectral_init(num_grid)

  type(tDict), pointer, intent(in) :: num_grid

  integer(MPI_INTEGER_KIND), dimension(0:worldsize-1) :: cells3_global
  integer :: ce
  DM :: DM_thermal
  real(pREAL), dimension(:,:,:), pointer :: T                                                       ! 0-indexed
  integer(MPI_INTEGER_KIND) :: err_MPI
  PetscErrorCode :: err_PETSc
  integer(HID_T) :: fileHandle, groupHandle
  real(pREAL), dimension(1,product(cells(1:2))*cells3) :: tempN
  type(tDict), pointer :: &
    num_grid_thermal
  character(len=:), allocatable :: &
    extmsg, &
    petsc_options


  print'(/,1x,a)', '<<<+-  grid_thermal_spectral init  -+>>>'

  print'(/,1x,a)', 'P. Shanthraj et al., Handbook of Mechanics of Materials, 2019'
  print'(  1x,a)', 'https://doi.org/10.1007/978-981-10-6855-3_80'

  if (.not. homogenization_thermal_active()) call IO_error(501,ext_msg='thermal')

!-------------------------------------------------------------------------------------------------
! read numerical parameters and do sanity checks
  num_grid_thermal => num_grid%get_dict('thermal',defaultVal=emptyDict)

  num%itmax            = num_grid_thermal%get_asInt('N_iter_max', defaultVal=100)
  num%eps_thermal_atol = num_grid_thermal%get_asReal('eps_abs_T', defaultVal=1.0e-2_pREAL)
  num%eps_thermal_rtol = num_grid_thermal%get_asReal('eps_rel_T', defaultVal=1.0e-6_pREAL)

  extmsg = ''
  if (num%eps_thermal_atol <= 0.0_pREAL)        extmsg = trim(extmsg)//' eps_abs_T'
  if (num%eps_thermal_rtol <= 0.0_pREAL)        extmsg = trim(extmsg)//' eps_rel_T'
  if (num%itmax < 1)                            extmsg = trim(extmsg)//' N_iter_max'

  if (extmsg /= '') call IO_error(301,ext_msg=trim(extmsg))
!--------------------------------------------------------------------------------------------------
! set default and user defined options for PETSc
  petsc_options = misc_prefixOptions('-snes_type newtonls -snes_mf -snes_ksp_ew -ksp_type fgmres '// &
                                     num_grid_thermal%get_asStr('PETSc_options',defaultVal=''), 'thermal_')
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,petsc_options,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
  call SNESCreate(PETSC_COMM_WORLD,SNES_thermal,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetOptionsPrefix(SNES_thermal,'thermal_',err_PETSc)
  CHKERRQ(err_PETSc)
  call MPI_Allgather(int(cells3,pPETSCINT),1_MPI_INTEGER_KIND,MPI_INTEGER,&
                     cells3_global,1_MPI_INTEGER_KIND,MPI_INTEGER,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call DMDACreate3D(PETSC_COMM_WORLD, &
         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &                                    ! cut off stencil at boundary
         DMDA_STENCIL_BOX, &                                                                        ! Moore (26) neighborhood around central point
         int(cells(1),pPETSCINT),int(cells(2),pPETSCINT),int(cells(3),pPETSCINT), &                 ! global cells
         1_pPETSCINT, 1_pPETSCINT, int(worldsize,pPETSCINT), &
         1_pPETSCINT, 0_pPETSCINT, &                                                                ! #dof (T, scalar), ghost boundary width (domain overlap)
         [int(cells(1),pPETSCINT)],[int(cells(2),pPETSCINT)],int(cells3_global,pPETSCINT), &        ! local cells
         DM_thermal,err_PETSc)                                                                      ! handle, error
  CHKERRQ(err_PETSc)
  call DMsetFromOptions(DM_thermal,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMsetUp(DM_thermal,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(DM_thermal,T_PETSc,err_PETSc)                                           ! global solution vector (cells x 1, i.e. every def grad tensor)
  CHKERRQ(err_PETSc)
  call DMDASNESSetFunctionLocal(DM_thermal,INSERT_VALUES,formResidual,PETSC_NULL_SNES,err_PETSc)    ! residual vector of same shape as solution vector
  CHKERRQ(err_PETSc)
  call SNESSetDM(SNES_thermal,DM_thermal,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetFromOptions(SNES_thermal,err_PETSc)                                                   ! pull it all together with additional CLI arguments
  CHKERRQ(err_PETSc)

  call DMDAVecGetArrayF90(DM_thermal,T_PETSc,T,err_PETSc)                                           ! returns 0-indexed T
  CHKERRQ(err_PETSc)

  restartRead: if (CLI_restartInc > 0) then
    print'(/,1x,a,1x,i0)', 'loading restart data of increment', CLI_restartInc

    fileHandle  = HDF5_openFile(getSolverJobName()//'_restart.hdf5','r')
    groupHandle = HDF5_openGroup(fileHandle,'solver')

    call HDF5_read(tempN,groupHandle,'T',.false.)
    T = reshape(tempN,[cells(1),cells(2),cells3])
    call HDF5_read(tempN,groupHandle,'T_lastInc',.false.)
    T_lastInc = reshape(tempN,[cells(1),cells(2),cells3])
    T_stagInc = T_lastInc
    call HDF5_read(tempN,groupHandle,'dotT_lastInc',.false.)
    dotT_lastInc = reshape(tempN,[cells(1),cells(2),cells3])
  else
    T = discretization_grid_getInitialCondition('T')
    T_lastInc = T(0:,0:,0:)
    T_stagInc = T_lastInc
    dotT_lastInc = 0.0_pREAL * T_lastInc
  end if restartRead

  call homogenization_thermal_setField(reshape(T,[product(cells(1:2))*cells3]), &
                                       [(0.0_pReal, ce = 1,product(cells(1:2))*cells3)])

  call DMDAVecRestoreArrayF90(DM_thermal,T_PETSc,T,err_PETSc)
  CHKERRQ(err_PETSc)

  call updateReference()

end subroutine grid_thermal_spectral_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the spectral thermal scheme with internal iterations
!--------------------------------------------------------------------------------------------------
function grid_thermal_spectral_solution(Delta_t) result(solution)

  real(pREAL), intent(in) ::  Delta_t                                                               !< increment in time for current solution

  type(tSolutionState) :: solution
  PetscInt  :: devNull
  PetscReal :: T_min, T_max, stagNorm
  DM :: DM_thermal
  real(pREAL), dimension(:,:,:), pointer :: T                                                       ! 0-indexed
  integer(MPI_INTEGER_KIND) :: err_MPI
  PetscErrorCode :: err_PETSc
  SNESConvergedReason :: reason


!--------------------------------------------------------------------------------------------------
! set module wide availabe data
  Delta_t_ = Delta_t

  call SNESSolve(SNES_thermal,PETSC_NULL_VEC,T_PETSc,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESGetConvergedReason(SNES_thermal,reason,err_PETSc)
  CHKERRQ(err_PETSc)

  solution%converged = reason > 0 .and. status == STATUS_OK
  solution%iterationsNeeded = merge(totalIter,num%itmax,solution%converged)

  call SNESGetDM(SNES_thermal,DM_thermal,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(DM_thermal,T_PETSc,T,err_PETSc)                                           ! returns 0-indexed T
  CHKERRQ(err_PETSc)

  T_min = minval(T)
  T_max = maxval(T)
  stagNorm = maxval(abs(T - T_stagInc))
  call MPI_Allreduce(MPI_IN_PLACE,stagNorm,1_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  solution%stagConverged = stagNorm < max(num%eps_thermal_atol, num%eps_thermal_rtol*T_max)
  call MPI_Allreduce(MPI_IN_PLACE,solution%stagConverged,1_MPI_INTEGER_KIND,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  T_stagInc = T

  call homogenization_thermal_setField(reshape(T,[product(cells(1:2))*cells3]), &
                                       reshape(T-T_lastInc,[product(cells(1:2))*cells3])/Delta_t_)

  call DMDAVecRestoreArrayF90(DM_thermal,T_PETSc,T,err_PETSc)
  CHKERRQ(err_PETSc)

  if (solution%converged) &
    print'(/,1x,a)', '... thermal conduction converged ..................................'
  print'(/,1x,a,f8.4,2x,f8.4,2x,f8.4)', 'Minimum|Maximum|Delta Temperature / K = ', T_min, T_max, stagNorm
  print'(/,1x,a)', '==========================================================================='
  flush(IO_STDOUT)

end function grid_thermal_spectral_solution


!--------------------------------------------------------------------------------------------------
!> @brief Set DAMASK data to current solver status.
!--------------------------------------------------------------------------------------------------
subroutine grid_thermal_spectral_forward(cutBack)

  logical, intent(in) :: cutBack

  DM :: DM_thermal
  real(pREAL),  dimension(:,:,:), pointer :: T                                                      ! 0-indexed
  PetscErrorCode :: err_PETSc


  call SNESGetDM(SNES_thermal,DM_thermal,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(DM_thermal,T_PETSc,T,err_PETSc)                                           ! returns 0-indexed T
  CHKERRQ(err_PETSc)

  if (cutBack) then
    call homogenization_thermal_setField(reshape(T_lastInc,[product(cells(1:2))*cells3]), &
                                         reshape(dotT_lastInc,[product(cells(1:2))*cells3]))
    T = T_lastInc
    T_stagInc = T_lastInc
  else
    dotT_lastInc = (T - T_lastInc)/Delta_t_
    T_lastInc = T
    call updateReference()
  end if

  call DMDAVecRestoreArrayF90(DM_thermal,T_PETSc,T,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine grid_thermal_spectral_forward


!--------------------------------------------------------------------------------------------------
!> @brief Write current solver and constitutive data for restart to file.
!--------------------------------------------------------------------------------------------------
subroutine grid_thermal_spectral_restartWrite()

  PetscErrorCode :: err_PETSc
  DM :: DM_thermal
  integer(HID_T) :: fileHandle, groupHandle
  real(pREAL), dimension(:,:,:), pointer :: T                                                       ! 0-indexed


  call SNESGetDM(SNES_thermal,DM_thermal,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayReadF90(DM_thermal,T_PETSc,T,err_PETSc)                                       ! returns 0-indexed T
  CHKERRQ(err_PETSc)

  print'(1x,a)', 'saving thermal solver data required for restart'; flush(IO_STDOUT)

  fileHandle  = HDF5_openFile(getSolverJobName()//'_restart.hdf5','a')
  groupHandle = HDF5_openGroup(fileHandle,'solver')
  call HDF5_write(reshape(T,[1,product(cells(1:2))*cells3]),groupHandle,'T')
  call HDF5_write(reshape(T_lastInc,[1,product(cells(1:2))*cells3]),groupHandle,'T_lastInc')
  call HDF5_write(reshape(dotT_lastInc,[1,product(cells(1:2))*cells3]),groupHandle,'dotT_lastInc')
  call HDF5_closeGroup(groupHandle)
  call HDF5_closeFile(fileHandle)

  call DMDAVecRestoreArrayReadF90(DM_thermal,T_PETSc,T,err_PETSc);
  CHKERRQ(err_PETSc)

end subroutine grid_thermal_spectral_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief Construct the residual vector.
!--------------------------------------------------------------------------------------------------
subroutine formResidual(residual_subdomain,x_scal,r,dummy,err_PETSc)

  DMDALocalInfo, dimension(DMDA_LOCAL_INFO_SIZE) :: &
    residual_subdomain
  real(pREAL), dimension(cells(1),cells(2),cells3), intent(in) :: &
    x_scal
  real(pREAL), dimension(cells(1),cells(2),cells3), intent(out) :: &
    r                                                                                               !< residual
  PetscObject :: dummy
  PetscErrorCode, intent(out) :: err_PETSc

  integer :: i, j, k, ce
  real(pREAL), dimension(3,cells(1),cells(2),cells3) :: vectorField


  call homogenization_thermal_response(status,Delta_t_,1,product(cells(1:2))*cells3)

  associate(T => x_scal)
    vectorField = utilities_ScalarGradient(T)
    ce = 0
    do k = 1, cells3;  do j = 1, cells(2);  do i = 1,cells(1)
      ce = ce + 1
      vectorField(1:3,i,j,k) = matmul(homogenization_K_T(ce) - K_ref, vectorField(1:3,i,j,k))
    end do; end do; end do
    r = utilities_VectorDivergence(vectorField)
    ce = 0
    do k = 1, cells3;  do j = 1, cells(2);  do i = 1,cells(1)
      ce = ce + 1
      r(i,j,k) = Delta_t_*(r(i,j,k) + homogenization_f_T(ce)) &
               + homogenization_mu_T(ce) * (T_lastInc(i,j,k) - T(i,j,k)) &
               + mu_ref*T(i,j,k)
    end do; end do; end do

    r = T &
      - utilities_GreenConvolution(r, K_ref, mu_ref, Delta_t_)
  end associate
  err_PETSc = 0

end subroutine formResidual


!--------------------------------------------------------------------------------------------------
!> @brief Update reference viscosity and conductivity.
!--------------------------------------------------------------------------------------------------
subroutine updateReference()

  integer :: ce
  integer(MPI_INTEGER_KIND) :: err_MPI


  K_ref = 0.0_pREAL
  mu_ref = 0.0_pREAL
  do ce = 1, product(cells(1:2))*cells3
    K_ref  = K_ref  + homogenization_K_T(ce)
    mu_ref = mu_ref + homogenization_mu_T(ce)
  end do

  K_ref = K_ref*wgt
  call MPI_Allreduce(MPI_IN_PLACE,K_ref,9_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  mu_ref = mu_ref*wgt
  call MPI_Allreduce(MPI_IN_PLACE,mu_ref,1_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

end subroutine updateReference


end module grid_thermal_spectral
