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
  use prec, only: & 
    pReal
  use spectral_utilities, only: &
    tSolutionState, &
    tSolutionParams
 
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
    temperature_current, &                                                                          !< field of current temperature
    temperature_lastInc, &                                                                          !< field of previous temperature
    temperature_stagInc                                                                             !< field of staggered temperature

!--------------------------------------------------------------------------------------------------
! reference diffusion tensor, mobility etc. 
  integer,                     private :: totalIter = 0                                             !< total iteration in current increment
  real(pReal), dimension(3,3), private :: D_ref
  real(pReal), private                 :: mobility_ref
  
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
  use spectral_utilities, only: &
    wgt
  use mesh, only: &
    grid, &
    grid3
  use thermal_conduction, only: &
    thermal_conduction_getConductivity33, &
    thermal_conduction_getMassDensity, &
    thermal_conduction_getSpecificHeat
  use material, only: &
    material_homogenizationAt, &
    temperature, &
    thermalMapping
  use numerics, only: &
    worldrank, &
    worldsize, &
    petsc_options

  implicit none
  PetscInt, dimension(worldsize) :: localK  
  integer :: i, j, k, cell
  DM :: thermal_grid
  PetscScalar,  dimension(:,:,:), pointer :: x_scal
  PetscErrorCode :: ierr

  write(6,'(/,a)') ' <<<+-  grid_thermal_spectral init  -+>>>'

  write(6,'(/,a)') ' Shanthraj et al., Handbook of Mechanics of Materials, 2019'
  write(6,'(a)')   ' https://doi.org/10.1007/978-981-10-6855-3_80'

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
  localK              = 0
  localK(worldrank+1) = grid3
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
  allocate(temperature_current(grid(1),grid(2),grid3), source=0.0_pReal)
  allocate(temperature_lastInc(grid(1),grid(2),grid3), source=0.0_pReal)
  allocate(temperature_stagInc(grid(1),grid(2),grid3), source=0.0_pReal)
  cell = 0
  do k = 1, grid3; do j = 1, grid(2); do i = 1,grid(1)
    cell = cell + 1
    temperature_current(i,j,k) = temperature(material_homogenizationAt(cell))% &
                                   p(thermalMapping(material_homogenizationAt(cell))%p(1,cell))
    temperature_lastInc(i,j,k) = temperature_current(i,j,k)
    temperature_stagInc(i,j,k) = temperature_current(i,j,k)
  enddo; enddo; enddo
  call DMDAVecGetArrayF90(thermal_grid,solution_vec,x_scal,ierr); CHKERRQ(ierr)                     !< get the data out of PETSc to work with
  x_scal(xstart:xend,ystart:yend,zstart:zend) = temperature_current
  call DMDAVecRestoreArrayF90(thermal_grid,solution_vec,x_scal,ierr); CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! thermal reference diffusion update
  cell = 0
  D_ref = 0.0_pReal
  mobility_ref = 0.0_pReal
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    cell = cell + 1
    D_ref = D_ref + thermal_conduction_getConductivity33(1,cell)
    mobility_ref = mobility_ref + thermal_conduction_getMassDensity(1,cell)* &
                                  thermal_conduction_getSpecificHeat(1,cell)
  enddo; enddo; enddo
  D_ref = D_ref*wgt
  call MPI_Allreduce(MPI_IN_PLACE,D_ref,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
  mobility_ref = mobility_ref*wgt
  call MPI_Allreduce(MPI_IN_PLACE,mobility_ref,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)

end subroutine grid_thermal_spectral_init

  
!--------------------------------------------------------------------------------------------------
!> @brief solution for the spectral thermal scheme with internal iterations
!--------------------------------------------------------------------------------------------------
function grid_thermal_spectral_solution(timeinc,timeinc_old,loadCaseTime) result(solution)
  use numerics, only: &
    itmax, &
    err_thermal_tolAbs, &
    err_thermal_tolRel
  use mesh, only: &
    grid, &
    grid3
  use thermal_conduction, only: &
    thermal_conduction_putTemperatureAndItsRate
 
  implicit none
  real(pReal), intent(in) :: &
    timeinc, &                                                                                      !< increment in time for current solution
    timeinc_old, &                                                                                  !< increment in time of last increment
    loadCaseTime                                                                                    !< remaining time of current load case
  integer :: i, j, k, cell
  type(tSolutionState) :: solution
  PetscInt  :: position
  PetscReal :: minTemperature, maxTemperature, stagNorm, solnNorm

  PetscErrorCode :: ierr   
  SNESConvergedReason :: reason

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
  stagNorm = maxval(abs(temperature_current - temperature_stagInc))
  solnNorm = maxval(abs(temperature_current))
  call MPI_Allreduce(MPI_IN_PLACE,stagNorm,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
  call MPI_Allreduce(MPI_IN_PLACE,solnNorm,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
  temperature_stagInc = temperature_current
  solution%stagConverged = stagNorm < min(err_thermal_tolAbs, err_thermal_tolRel*solnNorm)

!--------------------------------------------------------------------------------------------------
! updating thermal state 
  cell = 0
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    cell = cell + 1
    call thermal_conduction_putTemperatureAndItsRate(temperature_current(i,j,k), &
                                                     (temperature_current(i,j,k)-temperature_lastInc(i,j,k))/params%timeinc, &
                                                     1,cell)
  enddo; enddo; enddo

  call VecMin(solution_vec,position,minTemperature,ierr); CHKERRQ(ierr)
  call VecMax(solution_vec,position,maxTemperature,ierr); CHKERRQ(ierr) 
  if (solution%converged) &
    write(6,'(/,a)') ' ... thermal conduction converged ..................................'
  write(6,'(/,a,f8.4,2x,f8.4,2x,f8.4,/)',advance='no') ' Minimum|Maximum|Delta Temperature / K = ',&
                                                        minTemperature, maxTemperature, stagNorm
  write(6,'(/,a)') ' ==========================================================================='
  flush(6) 

end function grid_thermal_spectral_solution


!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine grid_thermal_spectral_forward
  use mesh, only: &
    grid, &
    grid3
  use spectral_utilities, only: &
    cutBack, &
    wgt
  use thermal_conduction, only: &
    thermal_conduction_putTemperatureAndItsRate, &
    thermal_conduction_getConductivity33, &
    thermal_conduction_getMassDensity, &
    thermal_conduction_getSpecificHeat
    
  implicit none
  integer :: i, j, k, cell
  DM :: dm_local
  PetscScalar,  dimension(:,:,:), pointer :: x_scal
  PetscErrorCode :: ierr
  
  if (cutBack) then 
    temperature_current = temperature_lastInc
    temperature_stagInc = temperature_lastInc

!--------------------------------------------------------------------------------------------------
! reverting thermal field state 
    cell = 0
    call SNESGetDM(thermal_snes,dm_local,ierr); CHKERRQ(ierr)
    call DMDAVecGetArrayF90(dm_local,solution_vec,x_scal,ierr); CHKERRQ(ierr)                       !< get the data out of PETSc to work with
    x_scal(xstart:xend,ystart:yend,zstart:zend) = temperature_current
    call DMDAVecRestoreArrayF90(dm_local,solution_vec,x_scal,ierr); CHKERRQ(ierr)
    do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
      cell = cell + 1
      call thermal_conduction_putTemperatureAndItsRate(temperature_current(i,j,k), &
                                                       (temperature_current(i,j,k) - &
                                                        temperature_lastInc(i,j,k))/params%timeinc, &
                                                       1,cell)
    enddo; enddo; enddo
  else
!--------------------------------------------------------------------------------------------------
! update rate and forward last inc
    temperature_lastInc = temperature_current
    cell = 0
    D_ref = 0.0_pReal
    mobility_ref = 0.0_pReal
    do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
      cell = cell + 1
      D_ref = D_ref + thermal_conduction_getConductivity33(1,cell)
      mobility_ref = mobility_ref + thermal_conduction_getMassDensity(1,cell)* &
                                    thermal_conduction_getSpecificHeat(1,cell)
    enddo; enddo; enddo
    D_ref = D_ref*wgt
    call MPI_Allreduce(MPI_IN_PLACE,D_ref,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
    mobility_ref = mobility_ref*wgt
    call MPI_Allreduce(MPI_IN_PLACE,mobility_ref,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
  endif
 
end subroutine grid_thermal_spectral_forward


!--------------------------------------------------------------------------------------------------
!> @brief forms the spectral thermal residual vector
!--------------------------------------------------------------------------------------------------
subroutine formResidual(in,x_scal,f_scal,dummy,ierr)
  use mesh, only: &
    grid, &
    grid3
  use math, only: &
    math_mul33x3
  use spectral_utilities, only: &
    scalarField_real, &
    vectorField_real, &
    utilities_FFTvectorForward, &
    utilities_FFTvectorBackward, &
    utilities_FFTscalarForward, &
    utilities_FFTscalarBackward, &
    utilities_fourierGreenConvolution, &
    utilities_fourierScalarGradient, &
    utilities_fourierVectorDivergence
  use thermal_conduction, only: &
    thermal_conduction_getSourceAndItsTangent, &
    thermal_conduction_getConductivity33, &
    thermal_conduction_getMassDensity, &
    thermal_conduction_getSpecificHeat
 
  implicit none
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

  temperature_current = x_scal 
!--------------------------------------------------------------------------------------------------
! evaluate polarization field
  scalarField_real = 0.0_pReal
  scalarField_real(1:grid(1),1:grid(2),1:grid3) = temperature_current 
  call utilities_FFTscalarForward
  call utilities_fourierScalarGradient                                                              !< calculate gradient of damage field
  call utilities_FFTvectorBackward
  cell = 0
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    cell = cell + 1
    vectorField_real(1:3,i,j,k) = math_mul33x3(thermal_conduction_getConductivity33(1,cell) - D_ref, &
                                               vectorField_real(1:3,i,j,k))
  enddo; enddo; enddo
  call utilities_FFTvectorForward
  call utilities_fourierVectorDivergence                                                            !< calculate damage divergence in fourier field
  call utilities_FFTscalarBackward
  cell = 0
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    cell = cell + 1
    call thermal_conduction_getSourceAndItsTangent(Tdot, dTdot_dT, temperature_current(i,j,k), 1, cell)
    scalarField_real(i,j,k) = params%timeinc*scalarField_real(i,j,k) + &
                              params%timeinc*Tdot + &
                              thermal_conduction_getMassDensity (1,cell)* &
                              thermal_conduction_getSpecificHeat(1,cell)*(temperature_lastInc(i,j,k)  - &
                                                                          temperature_current(i,j,k)) + &
                              mobility_ref*temperature_current(i,j,k)
  enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! convolution of damage field with green operator
  call utilities_FFTscalarForward
  call utilities_fourierGreenConvolution(D_ref, mobility_ref, params%timeinc)
  call utilities_FFTscalarBackward
 
!--------------------------------------------------------------------------------------------------
! constructing residual
  f_scal = temperature_current - scalarField_real(1:grid(1),1:grid(2),1:grid3)

end subroutine formResidual

end module grid_thermal_spectral
