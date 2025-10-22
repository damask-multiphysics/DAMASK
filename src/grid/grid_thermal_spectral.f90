!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Shaokang Zhang, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Spectral solver for thermal conduction
!--------------------------------------------------------------------------------------------------
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscdmda.h>
module grid_thermal_spectral
  use PETScDMDA
  use PETScSNES
#ifndef PETSC_HAVE_MPI_F90MODULE_VISIBILITY
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

#ifndef PETSC_HAVE_MPI_F90MODULE_VISIBILITY
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
  Vec :: T_vec
  real(pREAL), dimension(:,:,:), allocatable :: &
    T_lastinc, &                                                                                    !< field of previous temperature
    T_staginc, &                                                                                    !< field of staggered temperature
    T_dot_lastinc
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
subroutine grid_thermal_spectral_init(num_grid_thermal)

  type(tDict), pointer, intent(in) :: num_grid_thermal

  integer(MPI_INTEGER_KIND), dimension(0:worldsize-1) :: cells3_global
  integer :: ce
  DM :: DM_thermal
  real(pREAL), dimension(:,:,:), pointer :: T                                                       ! 0-indexed
  integer(MPI_INTEGER_KIND) :: err_MPI
  PetscErrorCode :: err_PETSc
  integer(HID_T) :: fileHandle, groupHandle
  real(pREAL), dimension(1,product(cells(1:2))*cells3) :: tempN
  character(len=:), allocatable :: &
    extmsg, &
    petsc_options


  print'(/,1x,a)', '<<<+-  grid_thermal_spectral init  -+>>>'

  print'(/,1x,a)', 'P. Shanthraj et al., Handbook of Mechanics of Materials, 2019'
  print'(  1x,a)', 'https://doi.org/10.1007/978-981-10-6855-3_80'

  if (.not. homogenization_thermal_active()) call IO_error(501,ext_msg='thermal')

!-------------------------------------------------------------------------------------------------
! read numerical parameters and do sanity checks
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
  call DMCreateGlobalVector(DM_thermal,T_vec,err_PETSc)                                             ! global solution vector (cells x 1, i.e. every def grad tensor)
  CHKERRQ(err_PETSc)
  call DMDASNESSetFunctionLocal(DM_thermal,INSERT_VALUES,form_residual,PETSC_NULL_SNES,err_PETSc)   ! residual vector of same shape as solution vector
  CHKERRQ(err_PETSc)
  call SNESSetDM(SNES_thermal,DM_thermal,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetFromOptions(SNES_thermal,err_PETSc)                                                   ! pull it all together with additional CLI arguments
  CHKERRQ(err_PETSc)

  call DMDAVecGetArray(DM_thermal,T_vec,T,err_PETSc)                                                ! returns 0-indexed T
  CHKERRQ(err_PETSc)

  restartRead: if (CLI_restartInc > 0) then
    print'(/,1x,a,1x,i0)', 'loading restart data of increment', CLI_restartInc

    fileHandle  = HDF5_openFile(CLI_jobName//'_restart.hdf5','r')
    groupHandle = HDF5_openGroup(fileHandle,'solver')

    call HDF5_read(tempN,groupHandle,'T',.false.)
    T = reshape(tempN,[cells(1),cells(2),cells3])
    call HDF5_read(tempN,groupHandle,'T_lastinc',.false.)                                           ! ToDo 4.0: rename dataset name
    T_lastinc = reshape(tempN,[cells(1),cells(2),cells3])
    T_staginc = T_lastinc
    call HDF5_read(tempN,groupHandle,'dotT_lastinc',.false.)
    T_dot_lastinc = reshape(tempN,[cells(1),cells(2),cells3])                                       ! ToDo 4.0: rename dataset name
  else
    T = discretization_grid_getScalarInitialCondition('T')
    T_lastinc = T(0:,0:,lbound(T,3):)
    T_staginc = T_lastinc
    T_dot_lastinc = 0.0_pREAL * T_lastinc
  end if restartRead

  call homogenization_thermal_setField(reshape(T,[product(cells(1:2))*cells3]), &
                                       [(0.0_pREAL, ce = 1,product(cells(1:2))*cells3)])

  call DMDAVecRestoreArray(DM_thermal,T_vec,T,err_PETSc)
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

  call SNESSolve(SNES_thermal,PETSC_NULL_VEC,T_vec,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESGetConvergedReason(SNES_thermal,reason,err_PETSc)
  CHKERRQ(err_PETSc)

#if PETSC_VERSION_MINOR<23
  solution%converged = reason > SNES_CONVERGED_ITERATING .and. status == STATUS_OK
#else
  solution%converged = reason%v > SNES_CONVERGED_ITERATING%v .and. status == STATUS_OK
#endif
  solution%iterationsNeeded = merge(totalIter,num%itmax,solution%converged)

  call VecMin(T_vec,devNull,T_min,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecMax(T_vec,devNull,T_max,err_PETSc)
  CHKERRQ(err_PETSc)

  call SNESGetDM(SNES_thermal,DM_thermal,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayRead(DM_thermal,T_vec,T,err_PETSc)                                            ! returns 0-indexed T
  CHKERRQ(err_PETSc)

  stagNorm = maxval(abs(T - T_staginc))
  call MPI_Allreduce(MPI_IN_PLACE,stagNorm,1_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  solution%stagConverged = stagNorm < max(num%eps_thermal_atol, num%eps_thermal_rtol*T_max)
  call MPI_Allreduce(MPI_IN_PLACE,solution%stagConverged,1_MPI_INTEGER_KIND,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  T_staginc = T

  call homogenization_thermal_setField(reshape(T,[product(cells(1:2))*cells3]), &
                                       reshape(T-T_lastinc,[product(cells(1:2))*cells3])/Delta_t_)

  call DMDAVecRestoreArrayRead(DM_thermal,T_vec,T,err_PETSc)
  CHKERRQ(err_PETSc)

  if (solution%converged) &
    print'(/,1x,a)', '... thermal conduction converged ..................................'
  print'(/,1x,a,f8.3,2x,f8.3,2x,f8.4)', 'Minimum|Maximum|Delta Temperature / K = ', T_min, T_max, stagNorm
  print'(/,1x,a)', '==========================================================================='
  flush(IO_STDOUT)

end function grid_thermal_spectral_solution


!--------------------------------------------------------------------------------------------------
!> @brief Set DAMASK data to current solver status.
!--------------------------------------------------------------------------------------------------
subroutine grid_thermal_spectral_forward(cutBack, Delta_t)

  logical,     intent(in) :: cutBack
  real(pREAL), intent(in) :: Delta_t

  DM :: DM_thermal
  real(pREAL),  dimension(:,:,:), pointer :: T                                                      ! 0-indexed
  PetscErrorCode :: err_PETSc


  call SNESGetDM(SNES_thermal,DM_thermal,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArray(DM_thermal,T_vec,T,err_PETSc)                                                ! returns 0-indexed T
  CHKERRQ(err_PETSc)

  if (cutBack) then
    call homogenization_thermal_setField(reshape(T_lastinc,[product(cells(1:2))*cells3]), &
                                         reshape(T_dot_lastinc,[product(cells(1:2))*cells3]))
    T = T_lastinc
    T_staginc = T_lastinc
  else
    T_dot_lastinc = (T - T_lastinc)/Delta_t
    T_lastinc = T
    call updateReference()
  end if

  call DMDAVecRestoreArray(DM_thermal,T_vec,T,err_PETSc)
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
  call DMDAVecGetArrayRead(DM_thermal,T_vec,T,err_PETSc)                                            ! returns 0-indexed T
  CHKERRQ(err_PETSc)

  print'(1x,a)', 'saving thermal solver data required for restart'; flush(IO_STDOUT)

  fileHandle  = HDF5_openFile(CLI_jobName//'_restart.hdf5','a')
  groupHandle = HDF5_openGroup(fileHandle,'solver')
  call HDF5_write(reshape(T,[1,product(shape(T))]),groupHandle,'T')
  call HDF5_write(reshape(T_lastinc,[1,product(shape(T_lastinc))]),groupHandle,'T_lastinc')         ! ToDo 4.0: rename dataset name
  call HDF5_write(reshape(T_dot_lastinc,[1,product(shape(T_dot_lastinc))]),groupHandle,'dotT_lastinc') ! ToDo 4.0: rename dataset name
  call HDF5_closeGroup(groupHandle)
  call HDF5_closeFile(fileHandle)

  call DMDAVecRestoreArrayRead(DM_thermal,T_vec,T,err_PETSc);
  CHKERRQ(err_PETSc)

end subroutine grid_thermal_spectral_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief Construct the residual vector.
!--------------------------------------------------------------------------------------------------
subroutine form_residual(residual_subdomain,T,r,dummy,err_PETSc)

#if PETSC_VERSION_MINOR<22
  DMDALocalInfo, dimension(DMDA_LOCAL_INFO_SIZE) :: &
#else
  DMDALocalInfo :: &
#endif
    residual_subdomain
  real(pREAL), dimension(cells(1),cells(2),cells3), intent(in) :: &
    T                                                                                               !< temperature
  real(pREAL), dimension(cells(1),cells(2),cells3), intent(out) :: &
    r                                                                                               !< residual
  PetscObject :: dummy
  PetscErrorCode, intent(out) :: err_PETSc

  integer :: i, j, k, ce
  real(pREAL), dimension(3,cells(1),cells(2),cells3) :: vector_field


  call homogenization_thermal_response(status,Delta_t_,1,product(cells(1:2))*cells3)

  vector_field = utilities_ScalarGradient(T)
  ce = 0
  do k = 1, cells3;  do j = 1, cells(2);  do i = 1,cells(1)
    ce = ce + 1
    vector_field(1:3,i,j,k) = matmul(homogenization_K_T(ce) - K_ref, vector_field(1:3,i,j,k))
  end do; end do; end do
  r = utilities_VectorDivergence(vector_field)
  ce = 0
  do k = 1, cells3;  do j = 1, cells(2);  do i = 1,cells(1)
    ce = ce + 1
    r(i,j,k) = Delta_t_*(r(i,j,k) + homogenization_f_T(ce)) &
             + homogenization_mu_T(ce) * (T_lastinc(i,j,k) - T(i,j,k)) &
             + mu_ref*T(i,j,k)
  end do; end do; end do

  r = T &
    - utilities_GreenConvolution(r, K_ref, mu_ref, Delta_t_)
  err_PETSc = 0_pPETSCERRORCODE

end subroutine form_residual


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
