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
  use misc
  use CLI
  use HDF5_utilities
  use HDF5
  use spectral_utilities
  use discretization_grid
  use homogenization
  use types
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
    real(pREAL) :: &
      phi_min, &                                                                                    !< non-zero residual damage
      eps_damage_atol, &                                                                            !< absolute tolerance for damage evolution
      eps_damage_rtol                                                                               !< relative tolerance for damage evolution
  end type tNumerics

  type(tNumerics) :: num

!--------------------------------------------------------------------------------------------------
! PETSc data
  SNES :: SNES_damage
  Vec :: phi_PETSc
  real(pREAL), dimension(:,:,:), allocatable :: &
    phi_lastInc, &                                                                                  !< field of previous damage
    phi_stagInc                                                                                     !< field of staggered damage

!--------------------------------------------------------------------------------------------------
! reference diffusion tensor, mobility etc.
  integer                     :: totalIter = 0                                                      !< total iteration in current increment
  real(pREAL), dimension(3,3) :: K_ref
  real(pREAL)                 :: mu_ref, Delta_t_

  public :: &
    grid_damage_spectral_init, &
    grid_damage_spectral_solution, &
    grid_damage_spectral_restartWrite, &
    grid_damage_spectral_forward

contains

!--------------------------------------------------------------------------------------------------
!> @brief Allocate all necessary fields and fill them with data, potentially from restart file.
!--------------------------------------------------------------------------------------------------
subroutine grid_damage_spectral_init(num_grid)

  type(tDict), pointer, intent(in) :: num_grid

  integer(MPI_INTEGER_KIND), dimension(0:worldsize-1) :: cells3_global
  DM :: DM_damage
  real(pREAL), dimension(:,:,:), pointer :: phi                                                     ! 0-indexed
  Vec :: uBound, lBound
  integer(MPI_INTEGER_KIND) :: err_MPI
  PetscErrorCode :: err_PETSc
  integer(HID_T) :: fileHandle, groupHandle
  real(pREAL), dimension(1,product(cells(1:2))*cells3) :: tempN
  type(tDict), pointer :: &
    num_grid_damage
  character(len=pSTRLEN) :: &
    snes_type
  character(len=:), allocatable :: &
    extmsg, &
    petsc_options


  print'(/,1x,a)', '<<<+-  grid_spectral_damage init  -+>>>'

  print'(/,1x,a)', 'P. Shanthraj et al., Handbook of Mechanics of Materials, 2019'
  print'(  1x,a)', 'https://doi.org/10.1007/978-981-10-6855-3_80'

  if (.not. homogenization_damage_active()) call IO_error(501,ext_msg='damage')

!-------------------------------------------------------------------------------------------------
! read numerical parameters and do sanity checks
  num_grid_damage => num_grid%get_dict('damage',defaultVal=emptyDict)

  num%itmax           = num_grid_damage%get_asInt ('N_iter_max', defaultVal=100)
  num%eps_damage_atol = num_grid_damage%get_asReal('eps_abs_phi',defaultVal=1.0e-2_pREAL)
  num%eps_damage_rtol = num_grid_damage%get_asReal('eps_rel_phi',defaultVal=1.0e-6_pREAL)
  num%phi_min         = num_grid_damage%get_asReal('phi_min',    defaultVal=1.0e-6_pREAL)

  extmsg = ''
  if (num%eps_damage_atol <= 0.0_pREAL) extmsg = trim(extmsg)//' eps_abs_phi'
  if (num%eps_damage_rtol <= 0.0_pREAL) extmsg = trim(extmsg)//' eps_rel_phi'
  if (num%phi_min <= 0.0_pREAL)         extmsg = trim(extmsg)//' phi_min'
  if (num%itmax < 1)                    extmsg = trim(extmsg)//' N_iter_max'

  if (extmsg /= '') call IO_error(301,ext_msg=trim(extmsg))

!--------------------------------------------------------------------------------------------------
! set default and user defined options for PETSc
  petsc_options = misc_prefixOptions('-snes_type newtonls -snes_mf -snes_ksp_ew -ksp_type fgmres '// &
                                     num_grid_damage%get_asStr('PETSc_options',defaultVal=''),'damage_')
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,petsc_options,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
  call SNESCreate(PETSC_COMM_WORLD,SNES_damage,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetOptionsPrefix(SNES_damage,'damage_',err_PETSc)
  CHKERRQ(err_PETSc)
  call MPI_Allgather(int(cells3,MPI_INTEGER_KIND),1_MPI_INTEGER_KIND,MPI_INTEGER,&
                     cells3_global,1_MPI_INTEGER_KIND,MPI_INTEGER,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call DMDACreate3D(PETSC_COMM_WORLD, &
         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &                                    ! cut off stencil at boundary
         DMDA_STENCIL_BOX, &                                                                        ! Moore (26) neighborhood around central point
         int(cells(1),pPETSCINT),int(cells(2),pPETSCINT),int(cells(3),pPETSCINT), &                 ! global cells
         1_pPETSCINT, 1_pPETSCINT, int(worldsize,pPETSCINT), &
         1_pPETSCINT, 0_pPETSCINT, &                                                                ! #dof (phi, scalar), ghost boundary width (domain overlap)
         [int(cells(1),pPetscInt)],[int(cells(2),pPetscInt)],int(cells3_global,pPETSCINT), &        ! local cells
         DM_damage,err_PETSc)                                                                       ! handle, error
  CHKERRQ(err_PETSc)
  call DMsetFromOptions(DM_damage,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMsetUp(DM_damage,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(DM_damage,phi_PETSc,err_PETSc)                                          ! global solution vector (cells x 1, i.e. every def grad tensor)
  CHKERRQ(err_PETSc)
  call DMDASNESSetFunctionLocal(DM_damage,INSERT_VALUES,formResidual,PETSC_NULL_SNES,err_PETSc)     ! residual vector of same shape as solution vector
  CHKERRQ(err_PETSc)
  call SNESSetDM(SNES_damage,DM_damage,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetFromOptions(SNES_damage,err_PETSc)                                                    ! pull it all together with additional CLI arguments
  CHKERRQ(err_PETSc)
  call SNESGetType(SNES_damage,snes_type,err_PETSc)
  CHKERRQ(err_PETSc)
  if (trim(snes_type) == 'vinewtonrsls' .or. &
      trim(snes_type) == 'vinewtonssls') then
    call DMGetGlobalVector(DM_damage,lBound,err_PETSc)
    CHKERRQ(err_PETSc)
    call DMGetGlobalVector(DM_damage,uBound,err_PETSc)
    CHKERRQ(err_PETSc)
    call VecSet(lBound,0.0_pREAL,err_PETSc)
    CHKERRQ(err_PETSc)
    call VecSet(uBound,1.0_pREAL,err_PETSc)
    CHKERRQ(err_PETSc)
    call SNESVISetVariableBounds(SNES_damage,lBound,uBound,err_PETSc)                               ! variable bounds for variational inequalities
    CHKERRQ(err_PETSc)
    call DMRestoreGlobalVector(DM_damage,lBound,err_PETSc)
    CHKERRQ(err_PETSc)
    call DMRestoreGlobalVector(DM_damage,uBound,err_PETSc)
    CHKERRQ(err_PETSc)
  end if

  call DMDAVecGetArrayF90(DM_damage,phi_PETSc,phi,err_PETSc)                                        ! returns 0-indexed phi
  CHKERRQ(err_PETSc)

  restartRead: if (CLI_restartInc > 0) then
    print'(/,1x,a,1x,i0)', 'loading restart data of increment', CLI_restartInc

    fileHandle  = HDF5_openFile(getSolverJobName()//'_restart.hdf5','r')
    groupHandle = HDF5_openGroup(fileHandle,'solver')

    call HDF5_read(tempN,groupHandle,'phi',.false.)
    phi = reshape(tempN,[cells(1),cells(2),cells3])
    call HDF5_read(tempN,groupHandle,'phi_lastInc',.false.)
    phi_lastInc = reshape(tempN,[cells(1),cells(2),cells3])
    phi_stagInc = phi_lastInc
  else
    phi = discretization_grid_getInitialCondition('phi')
    phi_lastInc = phi(0:,0:,0:)
    phi_stagInc = phi_lastInc
  end if restartRead

  call homogenization_set_phi(reshape(phi,[product(cells(1:2))*cells3]))

  call DMDAVecRestoreArrayF90(DM_damage,phi_PETSc,phi,err_PETSc)
  CHKERRQ(err_PETSc)

  call updateReference()

end subroutine grid_damage_spectral_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the spectral damage scheme with internal iterations
!--------------------------------------------------------------------------------------------------
function grid_damage_spectral_solution(Delta_t) result(solution)

  real(pREAL), intent(in) :: Delta_t                                                                !< increment in time for current solution

  type(tSolutionState) :: solution
  PetscInt  :: devNull
  PetscReal :: phi_min, phi_max, stagNorm
  DM :: DM_damage
  real(pREAL), dimension(:,:,:), pointer :: phi                                                     ! 0-indexed
  integer(MPI_INTEGER_KIND) :: err_MPI
  PetscErrorCode :: err_PETSc
  SNESConvergedReason :: reason


!--------------------------------------------------------------------------------------------------
! set module wide availabe data
  Delta_t_ = Delta_t

  call SNESSolve(SNES_damage,PETSC_NULL_VEC,phi_PETSc,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESGetConvergedReason(SNES_damage,reason,err_PETSc)
  CHKERRQ(err_PETSc)

  solution%converged = reason > 0
  solution%iterationsNeeded = merge(totalIter,num%itmax,solution%converged)

  call SNESGetDM(SNES_damage,DM_damage,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(DM_damage,phi_PETSc,phi,err_PETSc)                                        ! returns 0-indexed phi
  CHKERRQ(err_PETSc)

  phi_min = minval(phi)
  phi_max = maxval(phi)
  stagNorm = maxval(abs(phi - phi_stagInc))
  call MPI_Allreduce(MPI_IN_PLACE,stagNorm,1_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  solution%stagConverged = stagNorm < max(num%eps_damage_atol, num%eps_damage_rtol*phi_max)
  call MPI_Allreduce(MPI_IN_PLACE,solution%stagConverged,1_MPI_INTEGER_KIND,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  phi_stagInc = phi

  call homogenization_set_phi(reshape(phi,[product(cells(1:2))*cells3]))

  call DMDAVecRestoreArrayF90(DM_damage,phi_PETSc,phi,err_PETSc)
  CHKERRQ(err_PETSc)

  if (solution%converged) &
    print'(/,1x,a)', '... nonlocal damage converged .....................................'
  print'(/,1x,a,f8.6,2x,f8.6,2x,e11.4)', 'Minimum|Maximum|Delta Damage      = ', phi_min, phi_max, stagNorm
  print'(/,1x,a)', '==========================================================================='
  flush(IO_STDOUT)

end function grid_damage_spectral_solution


!--------------------------------------------------------------------------------------------------
!> @brief Set DAMASK data to current solver status.
!--------------------------------------------------------------------------------------------------
subroutine grid_damage_spectral_forward(cutBack)

  logical, intent(in) :: cutBack

  DM :: DM_damage
  real(pREAL),  dimension(:,:,:), pointer :: phi                                                    ! 0-indexed
  PetscErrorCode :: err_PETSc


  call SNESGetDM(SNES_damage,DM_damage,err_PETSc)
    CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(DM_damage,phi_PETSc,phi,err_PETSc)                                        ! returns 0-indexed T
    CHKERRQ(err_PETSc)

  if (cutBack) then
    call homogenization_set_phi(reshape(phi_lastInc,[product(cells(1:2))*cells3]))
    phi = phi_lastInc
    phi_stagInc = phi_lastInc
  else
    phi_lastInc = phi
    call updateReference()
  end if

end subroutine grid_damage_spectral_forward


!--------------------------------------------------------------------------------------------------
!> @brief Write current solver and constitutive data for restart to file.
!--------------------------------------------------------------------------------------------------
subroutine grid_damage_spectral_restartWrite()

  PetscErrorCode :: err_PETSc
  DM :: DM_damage
  integer(HID_T) :: fileHandle, groupHandle
  real(pREAL), dimension(:,:,:), pointer :: phi                                                     ! 0-indexed


  call SNESGetDM(SNES_damage,DM_damage,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayReadF90(DM_damage,phi_PETSc,phi,err_PETSc)                                    ! returns 0-indexed T
  CHKERRQ(err_PETSc)

  print'(1x,a)', 'saving damage solver data required for restart'; flush(IO_STDOUT)

  fileHandle  = HDF5_openFile(getSolverJobName()//'_restart.hdf5','a')
  groupHandle = HDF5_openGroup(fileHandle,'solver')
  call HDF5_write(reshape(phi,[1,product(cells(1:2))*cells3]),groupHandle,'phi')
  call HDF5_write(reshape(phi_lastInc,[1,product(cells(1:2))*cells3]),groupHandle,'phi_lastInc')
  call HDF5_closeGroup(groupHandle)
  call HDF5_closeFile(fileHandle)

  call DMDAVecRestoreArrayReadF90(DM_damage,phi_PETSc,phi,err_PETSc);
  CHKERRQ(err_PETSc)

end subroutine grid_damage_spectral_restartWrite


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


  associate(phi => x_scal)
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
      r(i,j,k) = Delta_t_*(r(i,j,k) + homogenization_f_phi(phi(i,j,k),ce)) &
               + homogenization_mu_phi(ce)*(phi_lastInc(i,j,k) - phi(i,j,k)) &
               + mu_ref*phi(i,j,k)
    end do; end do; end do

    r = max(min(utilities_GreenConvolution(r, K_ref, mu_ref, Delta_t_),phi_lastInc),num%phi_min) &
      - phi
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
    K_ref  = K_ref  + homogenization_K_phi(ce)
    mu_ref = mu_ref + homogenization_mu_phi(ce)
  end do

  K_ref = K_ref*wgt
  call MPI_Allreduce(MPI_IN_PLACE,K_ref,9_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  mu_ref = mu_ref*wgt
  call MPI_Allreduce(MPI_IN_PLACE,mu_ref,1_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

end subroutine updateReference


end module grid_damage_spectral
