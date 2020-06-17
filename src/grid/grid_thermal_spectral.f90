!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Shaokang Zhang, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Spectral solver for thermal conduction
!--------------------------------------------------------------------------------------------------
module grid_thermal_spectral
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscdmda.h>
  use PETScdmda
  use PETScsnes

  use prec
  use IO
  use spectral_utilities
  use discretization_grid
  use thermal_conduction
  use YAML_types
  use numerics
  use material
 
  implicit none
  private

!--------------------------------------------------------------------------------------------------
! derived types
  type(tSolutionParams), private :: params

!--------------------------------------------------------------------------------------------------
! PETSc data
  SNES,     private :: thermal_snes
  Vec,      private :: solution_vec
  PetscInt, private :: xstart, xend, ystart, yend, zstart, zend
  real(pReal), private, dimension(:,:,:), allocatable :: &
    T_current, &                                                                                    !< field of current temperature
    T_lastInc, &                                                                                    !< field of previous temperature
    T_stagInc                                                                                       !< field of staggered temperature

!--------------------------------------------------------------------------------------------------
! reference diffusion tensor, mobility etc. 
  integer,                     private :: totalIter = 0                                             !< total iteration in current increment
  real(pReal), dimension(3,3), private :: K_ref
  real(pReal), private                 :: mu_ref
  
  public :: &
    grid_thermal_spectral_init, &
    grid_thermal_spectral_solution, &
    grid_thermal_spectral_forward
  private :: &
    formResidual

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data
! ToDo: Restart not implemented
!--------------------------------------------------------------------------------------------------
subroutine grid_thermal_spectral_init

  PetscInt, dimension(0:worldsize-1) :: localK  
  integer :: i, j, k, cell
  DM :: thermal_grid
  PetscScalar,  dimension(:,:,:), pointer :: x_scal
  PetscErrorCode :: ierr
  class(tNode), pointer :: &
    num_generic
  character(len=pStringLen) :: &
    petsc_options

  write(6,'(/,a)') ' <<<+-  grid_thermal_spectral init  -+>>>'

  write(6,'(/,a)') ' Shanthraj et al., Handbook of Mechanics of Materials, 2019'
  write(6,'(a)')   ' https://doi.org/10.1007/978-981-10-6855-3_80'

!-------------------------------------------------------------------------------------------------
! read numerical parameter
  num_generic => numerics_root%get('generic',defaultVal=emptyDict)
  petsc_options = num_generic%get_asString('petsc_options',defaultVal='')

!--------------------------------------------------------------------------------------------------
! set default and user defined options for PETSc
 call PETScOptionsInsertString(PETSC_NULL_OPTIONS,'-thermal_snes_type ngmres',ierr)
 CHKERRQ(ierr)
 call PETScOptionsInsertString(PETSC_NULL_OPTIONS,trim(petsc_options),ierr)
 CHKERRQ(ierr)
 
!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
  call SNESCreate(PETSC_COMM_WORLD,thermal_snes,ierr); CHKERRQ(ierr)
  call SNESSetOptionsPrefix(thermal_snes,'thermal_',ierr);CHKERRQ(ierr) 
  localK            = 0
  localK(worldrank) = grid3
  call MPI_Allreduce(MPI_IN_PLACE,localK,worldsize,MPI_INTEGER,MPI_SUM,PETSC_COMM_WORLD,ierr)
  call DMDACreate3D(PETSC_COMM_WORLD, &
         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &                                    ! cut off stencil at boundary
         DMDA_STENCIL_BOX, &                                                                        ! Moore (26) neighborhood around central point
         grid(1),grid(2),grid(3), &                                                                 ! global grid
         1, 1, worldsize, &
         1, 0, &                                                                                    ! #dof (thermal phase field), ghost boundary width (domain overlap)
         [grid(1)],[grid(2)],localK, &                                                              ! local grid
         thermal_grid,ierr)                                                                         ! handle, error
  CHKERRQ(ierr)
  call SNESSetDM(thermal_snes,thermal_grid,ierr); CHKERRQ(ierr)                                     ! connect snes to da
  call DMsetFromOptions(thermal_grid,ierr); CHKERRQ(ierr)
  call DMsetUp(thermal_grid,ierr); CHKERRQ(ierr)
  call DMCreateGlobalVector(thermal_grid,solution_vec,ierr); CHKERRQ(ierr)                          ! global solution vector (grid x 1, i.e. every def grad tensor)
  call DMDASNESSetFunctionLocal(thermal_grid,INSERT_VALUES,formResidual,PETSC_NULL_SNES,ierr)       ! residual vector of same shape as solution vector
  CHKERRQ(ierr) 
  call SNESSetFromOptions(thermal_snes,ierr); CHKERRQ(ierr)                                         ! pull it all together with additional CLI arguments

!--------------------------------------------------------------------------------------------------
! init fields             
  call DMDAGetCorners(thermal_grid,xstart,ystart,zstart,xend,yend,zend,ierr)
  CHKERRQ(ierr)
  xend = xstart + xend - 1
  yend = ystart + yend - 1
  zend = zstart + zend - 1 
  allocate(T_current(grid(1),grid(2),grid3), source=0.0_pReal)
  allocate(T_lastInc(grid(1),grid(2),grid3), source=0.0_pReal)
  allocate(T_stagInc(grid(1),grid(2),grid3), source=0.0_pReal)
  cell = 0
  do k = 1, grid3; do j = 1, grid(2); do i = 1,grid(1)
    cell = cell + 1
    T_current(i,j,k) = temperature(material_homogenizationAt(cell))% &
                                   p(thermalMapping(material_homogenizationAt(cell))%p(1,cell))
    T_lastInc(i,j,k) = T_current(i,j,k)
    T_stagInc(i,j,k) = T_current(i,j,k)
  enddo; enddo; enddo
  call DMDAVecGetArrayF90(thermal_grid,solution_vec,x_scal,ierr); CHKERRQ(ierr)                     !< get the data out of PETSc to work with
  x_scal(xstart:xend,ystart:yend,zstart:zend) = T_current
  call DMDAVecRestoreArrayF90(thermal_grid,solution_vec,x_scal,ierr); CHKERRQ(ierr)

  call updateReference

end subroutine grid_thermal_spectral_init

  
!--------------------------------------------------------------------------------------------------
!> @brief solution for the spectral thermal scheme with internal iterations
!--------------------------------------------------------------------------------------------------
function grid_thermal_spectral_solution(timeinc,timeinc_old) result(solution)
 
  real(pReal), intent(in) :: &
    timeinc, &                                                                                      !< increment in time for current solution
    timeinc_old                                                                                     !< increment in time of last increment
  integer :: i, j, k, cell, &
    itmax                                                                                           !< maximum number of iterations
  type(tSolutionState) :: solution
  class(tNode), pointer :: &
    num_generic
  PetscInt  :: devNull
  PetscReal :: T_min, T_max, stagNorm, solnNorm

  PetscErrorCode :: ierr   
  SNESConvergedReason :: reason

!-------------------------------------------------------------------
! reading numerical parameter and do sanity check
  num_generic => numerics_root%get('generic',defaultVal=emptyDict)
  itmax = num_generic%get_asInt('itmax',defaultVal=250)
  if (itmax <= 1)   call IO_error(301,ext_msg='itmax')

  solution%converged =.false.
 
!--------------------------------------------------------------------------------------------------
! set module wide availabe data 
  params%timeinc = timeinc
  params%timeincOld = timeinc_old

  call SNESSolve(thermal_snes,PETSC_NULL_VEC,solution_vec,ierr); CHKERRQ(ierr)
  call SNESGetConvergedReason(thermal_snes,reason,ierr); CHKERRQ(ierr)

  if (reason < 1) then
    solution%converged = .false.
    solution%iterationsNeeded = itmax
  else
    solution%converged = .true.
    solution%iterationsNeeded = totalIter
  endif
  stagNorm = maxval(abs(T_current - T_stagInc))
  solnNorm = maxval(abs(T_current))
  call MPI_Allreduce(MPI_IN_PLACE,stagNorm,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
  call MPI_Allreduce(MPI_IN_PLACE,solnNorm,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
  T_stagInc = T_current
  solution%stagConverged = stagNorm < max(1.0e-2_pReal, 1.0e-6_pReal*solnNorm)

!--------------------------------------------------------------------------------------------------
! updating thermal state 
  cell = 0
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    cell = cell + 1
    call thermal_conduction_putTemperatureAndItsRate(T_current(i,j,k), &
                                                     (T_current(i,j,k)-T_lastInc(i,j,k))/params%timeinc, &
                                                     1,cell)
  enddo; enddo; enddo

  call VecMin(solution_vec,devNull,T_min,ierr); CHKERRQ(ierr)
  call VecMax(solution_vec,devNull,T_max,ierr); CHKERRQ(ierr)
  if (solution%converged) &
    write(6,'(/,a)') ' ... thermal conduction converged ..................................'
  write(6,'(/,a,f8.4,2x,f8.4,2x,f8.4,/)',advance='no') ' Minimum|Maximum|Delta Temperature / K = ',&
                                                        T_min, T_max, stagNorm
  write(6,'(/,a)') ' ==========================================================================='
  flush(6) 

end function grid_thermal_spectral_solution


!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine grid_thermal_spectral_forward(cutBack)
   
  logical, intent(in) :: cutBack
  integer :: i, j, k, cell
  DM :: dm_local
  PetscScalar,  dimension(:,:,:), pointer :: x_scal
  PetscErrorCode :: ierr
  
  if (cutBack) then 
    T_current = T_lastInc
    T_stagInc = T_lastInc

!--------------------------------------------------------------------------------------------------
! reverting thermal field state 
    cell = 0
    call SNESGetDM(thermal_snes,dm_local,ierr); CHKERRQ(ierr)
    call DMDAVecGetArrayF90(dm_local,solution_vec,x_scal,ierr); CHKERRQ(ierr)                       !< get the data out of PETSc to work with
    x_scal(xstart:xend,ystart:yend,zstart:zend) = T_current
    call DMDAVecRestoreArrayF90(dm_local,solution_vec,x_scal,ierr); CHKERRQ(ierr)
    do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
      cell = cell + 1
      call thermal_conduction_putTemperatureAndItsRate(T_current(i,j,k), &
                                                       (T_current(i,j,k) - &
                                                        T_lastInc(i,j,k))/params%timeinc, &
                                                       1,cell)
    enddo; enddo; enddo
  else
    T_lastInc = T_current
    call updateReference
  endif
 
end subroutine grid_thermal_spectral_forward


!--------------------------------------------------------------------------------------------------
!> @brief forms the spectral thermal residual vector
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
  integer :: i, j, k, cell
  real(pReal)   :: Tdot, dTdot_dT

  T_current = x_scal 
!--------------------------------------------------------------------------------------------------
! evaluate polarization field
  scalarField_real = 0.0_pReal
  scalarField_real(1:grid(1),1:grid(2),1:grid3) = T_current 
  call utilities_FFTscalarForward
  call utilities_fourierScalarGradient                                                              !< calculate gradient of temperature field
  call utilities_FFTvectorBackward
  cell = 0
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    cell = cell + 1
    vectorField_real(1:3,i,j,k) = matmul(thermal_conduction_getConductivity(1,cell) - K_ref, &
                                         vectorField_real(1:3,i,j,k))
  enddo; enddo; enddo
  call utilities_FFTvectorForward
  call utilities_fourierVectorDivergence                                                            !< calculate temperature divergence in fourier field
  call utilities_FFTscalarBackward
  cell = 0
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    cell = cell + 1
    call thermal_conduction_getSourceAndItsTangent(Tdot, dTdot_dT, T_current(i,j,k), 1, cell)
    scalarField_real(i,j,k) = params%timeinc*(scalarField_real(i,j,k) + Tdot) &
                            + thermal_conduction_getMassDensity (1,cell)* &
                              thermal_conduction_getSpecificHeat(1,cell)*(T_lastInc(i,j,k)  - &
                                                                          T_current(i,j,k))&
                            + mu_ref*T_current(i,j,k)
  enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! convolution of temperature field with green operator
  call utilities_FFTscalarForward
  call utilities_fourierGreenConvolution(K_ref, mu_ref, params%timeinc)
  call utilities_FFTscalarBackward
 
!--------------------------------------------------------------------------------------------------
! constructing residual
  f_scal = T_current - scalarField_real(1:grid(1),1:grid(2),1:grid3)

end subroutine formResidual


!--------------------------------------------------------------------------------------------------
!> @brief update reference viscosity and conductivity
!--------------------------------------------------------------------------------------------------
subroutine updateReference

  integer :: i,j,k,cell,ierr
  
  cell = 0
  K_ref = 0.0_pReal
  mu_ref = 0.0_pReal
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    cell = cell + 1
    K_ref  = K_ref  + thermal_conduction_getConductivity(1,cell)
    mu_ref = mu_ref + thermal_conduction_getMassDensity(1,cell)* thermal_conduction_getSpecificHeat(1,cell)
  enddo; enddo; enddo
  K_ref = K_ref*wgt
  call MPI_Allreduce(MPI_IN_PLACE,K_ref,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
  mu_ref = mu_ref*wgt
  call MPI_Allreduce(MPI_IN_PLACE,mu_ref,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)

end subroutine updateReference


end module grid_thermal_spectral
