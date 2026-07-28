!--------------------------------------------------------------------------------------------------
!> @author Tilak Raj Pant, Indian Institute of Science
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Grid solver for electricals: Spectral Lippmann-Schwinger
!--------------------------------------------------------------------------------------------------
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscdmda.h>
module grid_electrical_spectral
  use PETScDMDA
  use PETScSNES
#ifndef PETSC_EXPOSES_MPI
  use MPI_f08
#endif

  use prec
  use math
  use parallelization
  use IO
  use misc
  use CLI
  use HDF5_utilities
  use HDF5
  use spectral_utilities
  use grid_utilities
  use discretization_grid
  use homogenization
  use types
  use config
  use constants

#ifndef PETSC_EXPOSES_MPIF90
  implicit none(type,external)
#else
  implicit none
#endif
  private

  type :: tElectricalSolutionParams
    logical, dimension(3) :: E_mask = .true.   !< .true. => E_i unknown/free (Neumann)
    logical, dimension(3) :: J_mask = .true.   !< .true. => J_i unknown/free (Dirichlet)
  end type tElectricalSolutionParams

  type(tElectricalSolutionParams) :: params

  type :: tNumerics
    logical     :: update_gamma
    integer     :: itmin
    integer     :: itmax
    real(pREAL) :: eps_abs_E
    real(pREAL) :: eps_rel_E
    real(pREAL) :: eps_abs_J
    real(pREAL) :: eps_rel_J
  end type tNumerics

  type(tNumerics) :: num

  DM   :: DM_electrical
  SNES :: SNES_electrical
  Vec  :: E_PETSc

  real(pREAL), dimension(:,:,:,:), allocatable :: &
    E_lastInc, &
    E_stagInc

  real(pREAL), dimension(3) :: &
    E_aim         = 0.0_pREAL, &
    E_aim_lastInc = 0.0_pREAL, &
    J_aim         = 0.0_pREAL, &
    J_av          = 0.0_pREAL

  real(pREAL), dimension(3,3) :: sigma_ref, sigma_ref_inv

  real(pREAL)                   :: err_J     = 0.0_pREAL
  integer                       :: totalIter = 0
  integer(kind(STATUS_OK))      :: status
  character(len=:), allocatable :: incInfo
  PetscInt :: xs, ys, zs, xm, ym, zm

  public :: &
    grid_electrical_spectral_init, &
    grid_electrical_spectral_solution, &
    grid_electrical_spectral_restartWrite, &
    grid_electrical_spectral_forward

contains

!--------------------------------------------------------------------------------------------------
!> @brief Allocate all necessary fields and fill them with data, potentially from restart info.
!--------------------------------------------------------------------------------------------------
subroutine grid_electrical_spectral_init(num_grid)

  type(tDict), pointer, intent(in) :: num_grid

  integer(MPI_INTEGER_KIND), dimension(0:worldsize-1) :: cells3_global
  real(pREAL), pointer, dimension(:,:,:,:) :: E
  real(pREAL), dimension(3,product(cells(1:2))*cells3) :: tempN
  integer(HID_T)  :: fileHandle, groupHandle
  integer(MPI_INTEGER_KIND) :: err_MPI
  PetscErrorCode  :: err_PETSc
  type(tDict), pointer :: num_grid_electrical
  character(len=:), allocatable :: extmsg, petsc_options
  real(pREAL), dimension(3,cells(1),cells(2),cells3) :: J_init


  print'(/,1x,a)', '<<<+-  grid_electrical_spectral init  -+>>>'

  if (.not. homogenization_electrical_active()) call IO_error(501,ext_msg='electrical')

  num_grid_electrical => num_grid%get_dict('electrical',defaultVal=emptyDict)

  num%itmin        = num_grid_electrical%get_asInt ('N_iter_min',   defaultVal=1)
  num%itmax        = num_grid_electrical%get_asInt ('N_iter_max',   defaultVal=100)
  num%update_gamma = num_grid_electrical%get_asBool('update_gamma', defaultVal=.false.)
  num%eps_abs_E    = num_grid_electrical%get_asReal('eps_abs_E',    defaultVal=1.0e-4_pREAL)
  num%eps_rel_E    = num_grid_electrical%get_asReal('eps_rel_E',    defaultVal=1.0e-4_pREAL)
  num%eps_abs_J    = num_grid_electrical%get_asReal('eps_abs_J',    defaultVal=1.0e-4_pREAL)
  num%eps_rel_J    = num_grid_electrical%get_asReal('eps_rel_J',    defaultVal=1.0e-3_pREAL)

  extmsg = ''
  if (num%itmax < 1)                            extmsg = trim(extmsg)//' N_iter_max'
  if (num%itmin > num%itmax .or. num%itmin < 1) extmsg = trim(extmsg)//' N_iter_min'
  if (num%eps_abs_E <= 0.0_pREAL)               extmsg = trim(extmsg)//' eps_abs_E'
  if (num%eps_rel_E <= 0.0_pREAL)               extmsg = trim(extmsg)//' eps_rel_E'
  if (num%eps_abs_J <= 0.0_pREAL)               extmsg = trim(extmsg)//' eps_abs_J'
  if (num%eps_rel_J <= 0.0_pREAL)               extmsg = trim(extmsg)//' eps_rel_J'

  if (extmsg /= '') call IO_error(301,ext_msg=trim(extmsg))

  petsc_options = misc_prefixOptions( &
      '-snes_type anderson -snes_anderson_m 10 '// &
      num_grid_electrical%get_asStr('PETSc_options',defaultVal=''), &
      'electrical_')
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,petsc_options,err_PETSc)
  CHKERRQ(err_PETSc)

  allocate(E_lastInc(3,cells(1),cells(2),cells3), source=0.0_pREAL)
  allocate(E_stagInc(3,cells(1),cells(2),cells3), source=0.0_pREAL)

  call SNESCreate(PETSC_COMM_WORLD,SNES_electrical,err_PETSc); CHKERRQ(err_PETSc)
  call SNESSetOptionsPrefix(SNES_electrical,'electrical_',err_PETSc); CHKERRQ(err_PETSc)

  call MPI_Allgather(int(cells3,MPI_INTEGER_KIND),1_MPI_INTEGER_KIND,MPI_INTEGER, &
                     cells3_global,1_MPI_INTEGER_KIND,MPI_INTEGER,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

  call DMDACreate3D(PETSC_COMM_WORLD, &
         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
         DMDA_STENCIL_BOX, &
         int(cells(1),pPETSCINT), int(cells(2),pPETSCINT), int(cells(3),pPETSCINT), &
         1_pPETSCINT, 1_pPETSCINT, int(worldsize,pPETSCINT), &
         3_pPETSCINT, 0_pPETSCINT, &
         [int(cells(1),pPETSCINT)], [int(cells(2),pPETSCINT)], int(cells3_global,pPETSCINT), &
         DM_electrical,err_PETSc); CHKERRQ(err_PETSc)
  call DMsetFromOptions(DM_electrical,err_PETSc); CHKERRQ(err_PETSc)
  call DMsetUp(DM_electrical,err_PETSc);          CHKERRQ(err_PETSc)
  call DMcreateGlobalVector(DM_electrical,E_PETSc,err_PETSc); CHKERRQ(err_PETSc)
  call DMDASNESSetFunctionLocal(DM_electrical,INSERT_VALUES,form_residual, &
                                PETSC_NULL_SNES,err_PETSc); CHKERRQ(err_PETSc)
  call SNESSetConvergenceTest(SNES_electrical,converged, &
                              PETSC_NULL_SNES,PETSC_NULL_FUNCTION,err_PETSc); CHKERRQ(err_PETSc)
  call SNESSetDM(SNES_electrical,DM_electrical,err_PETSc); CHKERRQ(err_PETSc)
  call SNESSetFromOptions(SNES_electrical,err_PETSc);         CHKERRQ(err_PETSc)

  call DMDAVecGetArray(DM_electrical,E_PETSc,E,err_PETSc); CHKERRQ(err_PETSc)

  restartRead: if (CLI_restartInc /= -1) then
    print'(/,1x,a,1x,i0)', 'loading restart data of increment', CLI_restartInc

    fileHandle  = HDF5_openFile(CLI_jobName//'_restart.hdf5','r')
    groupHandle = HDF5_openGroup(fileHandle,'solver')

    call HDF5_read(E_aim,groupHandle,'E_aim',.false.)
    call HDF5_read(E_aim_lastInc,groupHandle,'E_aim_lastInc',.false.)
    call HDF5_read(tempN,groupHandle,'E')
    E(0:2,:,:,:) = reshape(tempN,[3,cells(1),cells(2),cells3])
    call HDF5_read(tempN,groupHandle,'E_lastInc')
    E_lastInc = reshape(tempN,[3,cells(1),cells(2),cells3])
    E_stagInc = E_lastInc

    call HDF5_closeGroup(groupHandle)
    call HDF5_closeFile(fileHandle)

  else restartRead
    E             = 0.0_pREAL
    E_lastInc     = 0.0_pREAL
    E_stagInc     = 0.0_pREAL
    E_aim         = 0.0_pREAL
    E_aim_lastInc = 0.0_pREAL

  end if restartRead

  call utilities_electricalResponse(status, J_init, J_av,E)

  call DMDAVecRestoreArray(DM_electrical,E_PETSc,E,err_PETSc); CHKERRQ(err_PETSc)

  call updateReference()

end subroutine grid_electrical_spectral_init


!--------------------------------------------------------------------------------------------------
!> @brief Solution — one SNESSolve call, no outer loop.
!--------------------------------------------------------------------------------------------------
function grid_electrical_spectral_solution(incInfoIn) result(solution)

  character(len=*), intent(in) :: incInfoIn
  type(tSolutionState)         :: solution

  PetscErrorCode      :: err_PETSc
  SNESConvergedReason :: reason
  real(pREAL), pointer, dimension(:,:,:,:) :: E
  PetscReal   :: norm_E
  real(pREAL), dimension(    cells(1),cells(2),cells3) :: temp_norm
  real(pREAL), dimension(3,  cells(1),cells(2),cells3) :: J_publish
  real(pREAL), dimension(3,product(cells(1:2))*cells3) :: E_publish
  integer(MPI_INTEGER_KIND) :: err_MPI
  integer :: i, j, k, ce

  incInfo = incInfoIn

  if (num%update_gamma) call updateReference()

  call SNESSolve(SNES_electrical,PETSC_NULL_VEC,E_PETSc,err_PETSc); CHKERRQ(err_PETSc)
  call SNESGetConvergedReason(SNES_electrical,reason,err_PETSc);    CHKERRQ(err_PETSc)

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<23)
  solution%converged = reason > SNES_CONVERGED_ITERATING
#else
  solution%converged = reason%v > SNES_CONVERGED_ITERATING%v
#endif

  call DMDAVecGetArrayRead(DM_electrical,E_PETSc,E,err_PETSc); CHKERRQ(err_PETSc)

  if (solution%converged) then
    print'(/,1x,a)', '... electrical conduction converged, publishing fields.'

call DMDAGetCorners(DM_electrical, xs, ys, zs, xm, ym, zm, err_PETSc)
CHKERRQ(err_PETSc)

  do k = zs, zs + zm - 1
    do j = 0, cells(2)-1
      do i = 0, cells(1)-1
        ce = (k - zs) * cells(2) * cells(1) + j * cells(1) + i + 1
        E_publish(1:3,ce) = E(0:2,i,j,k)
      end do
    end do
  end do

    ! 2. Safe execution: Update constituent states ONLY now
    call homogenization_electrical_response(status,E_publish,1,product(cells(1:2))*cells3)

    ! 3. Standard results logging
    call utilities_electricalResponse(status,J_publish,J_av,E)
  else
    print'(/,1x,a)', '... electrical conduction FAILED to converge.'
  end if

  temp_norm = sqrt(  (E(0,:,:,:) - E_stagInc(1,:,:,:))**2 &
                   + (E(1,:,:,:) - E_stagInc(2,:,:,:))**2 &
                   + (E(2,:,:,:) - E_stagInc(3,:,:,:))**2 )

  call VecNorm(E_PETSc,NORM_2,norm_E,err_PETSc); CHKERRQ(err_PETSc)

  solution%stagConverged = maxval(temp_norm) < max(num%eps_abs_E, num%eps_rel_E * norm_E)
  call MPI_Allreduce(MPI_IN_PLACE,solution%stagConverged,1_MPI_INTEGER_KIND, &
                     MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

  E_stagInc(1:3,:,:,:) = E(0:2,:,:,:)

  call DMDAVecRestoreArrayRead(DM_electrical,E_PETSc,E,err_PETSc); CHKERRQ(err_PETSc)

  print'(/,1x,a,es12.4)', 'Norm of E-field / V/m = ', norm_E
  flush(IO_STDOUT)

end function grid_electrical_spectral_solution


!--------------------------------------------------------------------------------------------------
!> @brief Set solver state at start of new increment or after cutback.
!--------------------------------------------------------------------------------------------------
subroutine grid_electrical_spectral_forward(cutBack,Delta_t,t_remaining, &
                                            electric_field_BC,current_density_BC)

  logical,                  intent(in) :: cutBack
  real(pREAL),              intent(in) :: Delta_t, t_remaining
  type(tBCelectrical),      intent(in) :: electric_field_BC, current_density_BC

  real(pREAL), pointer, dimension(:,:,:,:) :: E
  real(pREAL), dimension(3,cells(1),cells(2),cells3) :: J_restore
  PetscErrorCode :: err_PETSc
  integer :: i


  call DMDAVecGetArray(DM_electrical,E_PETSc,E,err_PETSc); CHKERRQ(err_PETSc)

  if (cutBack) then
    E_aim        = E_aim_lastInc
    E(0:2,:,:,:) = E_lastInc(1:3,:,:,:)
    E_stagInc    = E_lastInc
    call utilities_electricalResponse(status,J_restore,J_av,E_lastInc)                              ! restore constitutive state via utilities — single entry point
  else
    E_aim_lastInc        = E_aim
    E_lastInc(1:3,:,:,:) = E(0:2,:,:,:)
    E_stagInc            = E_lastInc
    call updateReference()
  end if

  call DMDAVecRestoreArray(DM_electrical,E_PETSc,E,err_PETSc); CHKERRQ(err_PETSc)

  params%E_mask = electric_field_BC%mask
  params%J_mask = current_density_BC%mask

  do i = 1, 3
    if (.not. params%E_mask(i)) E_aim(i) = electric_field_BC%values(i)  ! E_aim for Dirichlet components (E_mask=F means E is prescribed)
    if (.not. params%J_mask(i)) then
      J_aim(i) = current_density_BC%values(i)
    else
      J_aim(i) = 0.0_pREAL
    end if
  end do                                                                ! J_aim for Neumann components (J_mask=F means J is prescribed)

end subroutine grid_electrical_spectral_forward


!--------------------------------------------------------------------------------------------------
!> @brief Write current solver data for restart.
!--------------------------------------------------------------------------------------------------
subroutine grid_electrical_spectral_restartWrite()

  PetscErrorCode :: err_PETSc
  integer(HID_T) :: fileHandle, groupHandle
  real(pREAL), pointer, dimension(:,:,:,:) :: E


  call DMDAVecGetArrayRead(DM_electrical,E_PETSc,E,err_PETSc); CHKERRQ(err_PETSc)

  if (num%update_gamma) call updateReference()

  print'(1x,a)', 'saving electrical solver data required for restart'; flush(IO_STDOUT)

  fileHandle  = HDF5_openFile(CLI_jobName//'_restart.hdf5','a')
  groupHandle = HDF5_addGroup(fileHandle,'solver')
  call HDF5_write(reshape(E(0:2,:,:,:),[3,product(cells(1:2))*cells3]),groupHandle,'E')
  call HDF5_write(reshape(E_lastInc,[3,product(cells(1:2))*cells3]),groupHandle,'E_lastInc')
  call HDF5_closeGroup(groupHandle)
  call HDF5_closeFile(fileHandle)

  if (worldrank == 0) then
    fileHandle  = HDF5_openFile(CLI_jobName//'_restart.hdf5','a',.false.)
    groupHandle = HDF5_openGroup(fileHandle,'solver')
    call HDF5_write(E_aim,        groupHandle,'E_aim',        .false.)
    call HDF5_write(E_aim_lastInc,groupHandle,'E_aim_lastInc',.false.)
    call HDF5_closeGroup(groupHandle)
    call HDF5_closeFile(fileHandle)
  end if

  call DMDAVecRestoreArrayRead(DM_electrical,E_PETSc,E,err_PETSc); CHKERRQ(err_PETSc)

end subroutine grid_electrical_spectral_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief Custom PETSc convergence callback
!--------------------------------------------------------------------------------------------------
subroutine converged(snes_local,PETScIter,devNull1,devNull2,devNull3,reason,dummy,err_PETSc)

  SNES                  :: snes_local
  PetscInt,  intent(in) :: PETScIter
  PetscReal, intent(in) :: devNull1, devNull2, devNull3
  SNESConvergedReason   :: reason
  PetscObject           :: dummy
  PetscErrorCode        :: err_PETSc

  real(pREAL)           :: J_tol


  J_tol = max(maxval(abs(J_av)) * num%eps_rel_J, num%eps_abs_J)                                     ! tolerance scaled with |J_av|

  if (totalIter >= num%itmin .and. err_J < J_tol .and. status == STATUS_OK) then
    reason = SNES_CONVERGED_USER
  elseif (totalIter >= num%itmax) then
    reason = SNES_DIVERGED_USER
  else
    reason = SNES_CONVERGED_ITERATING
  end if

  print'(/,1x,a)', '... electrical reporting .........................................'
  print'(/,1x,a,f12.2,a,es8.2,a,es9.2,a)', 'error current BC = ', &
        err_J/J_tol,' (',err_J,' A/m^2, tol = ',J_tol,')'
  print'(1x,a,/,3(es12.4,1x))', 'E_aim / V/m   = ', E_aim
  print'(1x,a,/,3(es12.4,1x))', 'J_av  / A/m^2 = ', J_av
  print'(/,1x,a)', '==========================================================================='
  flush(IO_STDOUT)
  err_PETSc = 0

end subroutine converged


!--------------------------------------------------------------------------------------------------
!> @brief Construct residual and update E_aim (Neumann components) cleanly.
!> @details The residual is defined as R(e) = e - [E_aim + Gamma * p(e)]
!--------------------------------------------------------------------------------------------------
subroutine form_residual(residual_subdomain, x_vec, r_vec, dummy, err_PETSc)

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<22)
  DMDALocalInfo, dimension(DMDA_LOCAL_INFO_SIZE) :: &
#else
  DMDALocalInfo :: &
#endif
    residual_subdomain
  real(pREAL), dimension(3,cells(1),cells(2),cells3), intent(in)  :: x_vec                          !< input E-field guess from PETSc
  real(pREAL), dimension(3,cells(1),cells(2),cells3), intent(out) :: r_vec                          !< output residual vector to PETSc
  PetscObject                                                     :: dummy
  PetscErrorCode,                                     intent(out) :: err_PETSc

  integer :: i, j, k, ce
  real(pREAL), dimension(3,cells(1),cells(2),cells3) :: J_field, polarization, e_correction
  real(pREAL), dimension(3)   :: deltaE_aim
  real(pREAL), dimension(3,3) :: Psig
  integer(MPI_INTEGER_KIND)   :: err_MPI
  PetscInt :: PETScIter, nfuncs


  call SNESGetNumberFunctionEvals(SNES_electrical,nfuncs,err_PETSc); CHKERRQ(err_PETSc)
  call SNESGetIterationNumber(SNES_electrical,PETScIter,err_PETSc);  CHKERRQ(err_PETSc)

  if (nfuncs == 0 .and. PETScIter == 0) totalIter = -1

! ------------------------------------------------------------------------------------------------
!True New Iteration Guard (Updates E_aim strictly once per physical Newton step)
  newIteration: if (totalIter <= PETScIter) then
    totalIter = totalIter + 1
    print'(1x,a,3(a,i0))', trim(incInfo), &
         ' @ Electrical Iteration ', num%itmin,'<=',totalIter,'<=',num%itmax
    print'(/,1x,a,/,3(es12.4,1x))', 'E_aim / V/m = ', E_aim
    flush(IO_STDOUT)

    call utilities_electricalResponse(status, J_field, J_av, x_vec)                                ! safe, read-only calculation of the macro J tracking parameters

    call MPI_Allreduce(MPI_IN_PLACE,status,1_MPI_INTEGER_KIND, &
                       MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err_MPI)
    call parallelization_chkerr(err_MPI)

    deltaE_aim = matmul(sigma_ref_inv, J_av - J_aim)                                                ! update macro driving force for Neumann/Current components
    where (.not. params%E_mask) deltaE_aim = 0.0_pREAL                                              ! locked on Dirichlet (E prescribed)

    E_aim = E_aim - deltaE_aim
    err_J = maxval(abs(merge(0.0_pREAL, J_av - J_aim, params%J_mask)))                              ! update the boundary condition error tracker
  end if newIteration

! ------------------------------------------------------------------------------------------------
! Build Polarization Field (Local algebraic perturbation vector)
  !$OMP PARALLEL DO PRIVATE(ce, Psig)
  do k = 1, cells3
    do j = 1, cells(2)
      do i = 1, cells(1)
        ce = (k-1) * cells(2) * cells(1) + (j-1) * cells(1) + i
        Psig = matmul(sigma_ref_inv, homogenization_sigma_E(ce)) - math_I3
        polarization(:,i,j,k) = matmul(Psig, x_vec(:,i,j,k))
      end do
    end do
  end do
  !$OMP END PARALLEL DO

! ------------------------------------------------------------------------------------------------
! Spectral Lippmann-Schwinger Green Operator Convolution via FFT
  e_correction = utilities_GreenConvolution_electrical(polarization)

! ------------------------------------------------------------------------------------------------
! Construct Residual: R(e) = e - [E_aim + Gamma*p(e)]
  !$OMP PARALLEL DO PRIVATE(k, j, i) SHARED(r_vec, x_vec, E_aim, e_correction)
  do k = 1, cells3
    do j = 1, cells(2)
      do i = 1, cells(1)
        r_vec(:,i,j,k) = x_vec(:,i,j,k) - (E_aim + e_correction(:,i,j,k))
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  err_PETSc = 0

end subroutine form_residual


!--------------------------------------------------------------------------------------------------
!> @brief Update volume-averaged reference conductivity and refresh the Green operator.
!--------------------------------------------------------------------------------------------------
subroutine updateReference()

  integer :: ce, ierr
  integer(MPI_INTEGER_KIND) :: err_MPI
  real(pREAL) :: det

  sigma_ref = 0.0_pREAL
  do ce = 1, product(cells(1:2))*cells3
    sigma_ref = sigma_ref + homogenization_sigma_E(ce)
  end do
  sigma_ref = sigma_ref * wgt

  call MPI_Allreduce(MPI_IN_PLACE,sigma_ref,9_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM, &
                     MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

  det = math_det33(sigma_ref)

  if (abs(det) < 1.0e-20_pREAL) then
    write(IO_STDOUT,'(a,es12.4)') &
      'FATAL ERROR in updateReference: sigma_ref is singular; det = ', det
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end if

  sigma_ref_inv = math_inv33(sigma_ref)

  call utilities_updateGamma_electrical()

end subroutine updateReference

end module grid_electrical_spectral
