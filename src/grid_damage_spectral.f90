!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Shaokang Zhang, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Spectral solver for nonlocal damage
!--------------------------------------------------------------------------------------------------
module grid_damage_spectral
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
 SNES, private :: damage_snes
 Vec,  private :: solution_vec
 PetscInt, private :: xstart, xend, ystart, yend, zstart, zend
 real(pReal), private, dimension(:,:,:), allocatable :: &
   damage_current, &                                                                           !< field of current damage
   damage_lastInc, &                                                                           !< field of previous damage
   damage_stagInc                                                                              !< field of staggered damage

!--------------------------------------------------------------------------------------------------
! reference diffusion tensor, mobility etc. 
 integer,               private :: totalIter = 0                                              !< total iteration in current increment
 real(pReal), dimension(3,3), private :: D_ref
 real(pReal), private                 :: mobility_ref
 
 public :: &
   grid_damage_spectral_init, &
   grid_damage_spectral_solution, &
   grid_damage_spectral_forward
 private :: &
   formResidual

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data
! ToDo: Restart not implemented
!--------------------------------------------------------------------------------------------------
subroutine grid_damage_spectral_init
  use IO, only: &
    IO_intOut
  use spectral_utilities, only: &
    wgt
  use mesh, only: &
    grid, &
    grid3
  use damage_nonlocal, only: &
    damage_nonlocal_getDiffusion33, &
    damage_nonlocal_getMobility
  use numerics, only: &
    worldrank, &
    worldsize, &
    petsc_options
    
  implicit none
  PetscInt, dimension(worldsize) :: localK  
  integer :: i, j, k, cell
  DM :: damage_grid
  Vec :: uBound, lBound
  PetscErrorCode :: ierr
  character(len=100) :: snes_type
 
  write(6,'(/,a)') ' <<<+-  grid_spectral_damage init  -+>>>'

  write(6,'(/,a)') ' Shanthraj et al., Handbook of Mechanics of Materials, 2019'
  write(6,'(a)')   ' https://doi.org/10.1007/978-981-10-6855-3_80'
 
!--------------------------------------------------------------------------------------------------
! set default and user defined options for PETSc
 call PETScOptionsInsertString(PETSC_NULL_OPTIONS,'-damage_snes_type ngmres',ierr)
 CHKERRQ(ierr)
 call PETScOptionsInsertString(PETSC_NULL_OPTIONS,trim(petsc_options),ierr)
 CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
  call SNESCreate(PETSC_COMM_WORLD,damage_snes,ierr); CHKERRQ(ierr)
  call SNESSetOptionsPrefix(damage_snes,'damage_',ierr);CHKERRQ(ierr) 
  localK              = 0
  localK(worldrank+1) = grid3
  call MPI_Allreduce(MPI_IN_PLACE,localK,worldsize,MPI_INTEGER,MPI_SUM,PETSC_COMM_WORLD,ierr)
  call DMDACreate3D(PETSC_COMM_WORLD, &
         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &                                    !< cut off stencil at boundary
         DMDA_STENCIL_BOX, &                                                                        !< Moore (26) neighborhood around central point
         grid(1),grid(2),grid(3), &                                                                 !< global grid
         1, 1, worldsize, &
         1, 0, &                                                                                    !< #dof (damage phase field), ghost boundary width (domain overlap)
         [grid(1)],[grid(2)],localK, &                                                              !< local grid
         damage_grid,ierr)                                                                          !< handle, error
  CHKERRQ(ierr)
  call SNESSetDM(damage_snes,damage_grid,ierr); CHKERRQ(ierr)                                       !< connect snes to da
  call DMsetFromOptions(damage_grid,ierr); CHKERRQ(ierr)
  call DMsetUp(damage_grid,ierr); CHKERRQ(ierr)
  call DMCreateGlobalVector(damage_grid,solution_vec,ierr); CHKERRQ(ierr)                               !< global solution vector (grid x 1, i.e. every def grad tensor)
  call DMDASNESSetFunctionLocal(damage_grid,INSERT_VALUES,formResidual,&
                                                                             PETSC_NULL_SNES,ierr)  !< residual vector of same shape as solution vector
  CHKERRQ(ierr) 
  call SNESSetFromOptions(damage_snes,ierr); CHKERRQ(ierr)                                          !< pull it all together with additional CLI arguments
  call SNESGetType(damage_snes,snes_type,ierr); CHKERRQ(ierr)
  if (trim(snes_type) == 'vinewtonrsls' .or. &
      trim(snes_type) == 'vinewtonssls') then
    call DMGetGlobalVector(damage_grid,lBound,ierr); CHKERRQ(ierr)
    call DMGetGlobalVector(damage_grid,uBound,ierr); CHKERRQ(ierr)
    call VecSet(lBound,0.0_pReal,ierr); CHKERRQ(ierr)
    call VecSet(uBound,1.0_pReal,ierr); CHKERRQ(ierr)
    call SNESVISetVariableBounds(damage_snes,lBound,uBound,ierr)                                    !< variable bounds for variational inequalities like contact mechanics, damage etc.
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
  call VecSet(solution_vec,1.0_pReal,ierr); CHKERRQ(ierr)
  allocate(damage_current(grid(1),grid(2),grid3), source=1.0_pReal)
  allocate(damage_lastInc(grid(1),grid(2),grid3), source=1.0_pReal)
  allocate(damage_stagInc(grid(1),grid(2),grid3), source=1.0_pReal)

!--------------------------------------------------------------------------------------------------
! damage reference diffusion update
  cell = 0
  D_ref = 0.0_pReal
  mobility_ref = 0.0_pReal
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    cell = cell + 1
    D_ref = D_ref + damage_nonlocal_getDiffusion33(1,cell)
    mobility_ref = mobility_ref + damage_nonlocal_getMobility(1,cell)
  enddo; enddo; enddo
  D_ref = D_ref*wgt
  call MPI_Allreduce(MPI_IN_PLACE,D_ref,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
  mobility_ref = mobility_ref*wgt
  call MPI_Allreduce(MPI_IN_PLACE,mobility_ref,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)

end subroutine grid_damage_spectral_init
  
!--------------------------------------------------------------------------------------------------
!> @brief solution for the spectral damage scheme with internal iterations
!--------------------------------------------------------------------------------------------------
function grid_damage_spectral_solution(timeinc,timeinc_old,loadCaseTime) result(solution)
  use numerics, only: &
    itmax, &
    err_damage_tolAbs, &
    err_damage_tolRel
  use mesh, only: &
    grid, &
    grid3
  use damage_nonlocal, only: &
    damage_nonlocal_putNonLocalDamage
 
  implicit none
  real(pReal), intent(in) :: &
    timeinc, &                                                                                      !< increment in time for current solution
    timeinc_old, &                                                                                  !< increment in time of last increment
    loadCaseTime                                                                                    !< remaining time of current load case

  integer :: i, j, k, cell
  PetscInt  ::position
  PetscReal ::  minDamage, maxDamage, stagNorm, solnNorm
  PetscErrorCode :: ierr
 type(tSolutionState) :: &
   solution
  SNESConvergedReason :: reason

 solution%converged =.false.
 
!--------------------------------------------------------------------------------------------------
! set module wide availabe data 
  params%timeinc = timeinc
  params%timeincOld = timeinc_old
 
  call SNESSolve(damage_snes,PETSC_NULL_VEC,solution_vec,ierr); CHKERRQ(ierr)
  call SNESGetConvergedReason(damage_snes,reason,ierr); CHKERRQ(ierr)
 
  if (reason < 1) then
    solution%converged = .false.
    solution%iterationsNeeded = itmax
  else
    solution%converged = .true.
    solution%iterationsNeeded = totalIter
  endif
  stagNorm = maxval(abs(damage_current - damage_stagInc))
  solnNorm = maxval(abs(damage_current))
  call MPI_Allreduce(MPI_IN_PLACE,stagNorm,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
  call MPI_Allreduce(MPI_IN_PLACE,solnNorm,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
  damage_stagInc = damage_current
  solution%stagConverged = stagNorm < min(err_damage_tolAbs,err_damage_tolRel*solnNorm)

!--------------------------------------------------------------------------------------------------
! updating damage state 
  cell = 0
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    cell = cell + 1
    call damage_nonlocal_putNonLocalDamage(damage_current(i,j,k),1,cell)
  enddo; enddo; enddo
 
  call VecMin(solution_vec,position,minDamage,ierr); CHKERRQ(ierr)
  call VecMax(solution_vec,position,maxDamage,ierr); CHKERRQ(ierr)
  if (solution%converged) &
     write(6,'(/,a)') ' ... nonlocal damage converged .....................................'
   write(6,'(/,a,f8.6,2x,f8.6,2x,f8.6,/)',advance='no') ' Minimum|Maximum|Delta Damage      = ',&
                                                        minDamage, maxDamage, stagNorm
  write(6,'(/,a)') ' ==========================================================================='
  flush(6) 

end function grid_damage_spectral_solution


!--------------------------------------------------------------------------------------------------
!> @brief spectral damage forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine grid_damage_spectral_forward()
 use mesh, only: &
   grid, &
   grid3
 use spectral_utilities, only: &
   cutBack, &
   wgt
 use damage_nonlocal, only: &
   damage_nonlocal_putNonLocalDamage, &
   damage_nonlocal_getDiffusion33, &
   damage_nonlocal_getMobility
   
 implicit none
 integer                               :: i, j, k, cell
 DM :: dm_local
 PetscScalar,  dimension(:,:,:), pointer     :: x_scal
 PetscErrorCode                              :: ierr

 if (cutBack) then 
   damage_current = damage_lastInc
   damage_stagInc = damage_lastInc
!--------------------------------------------------------------------------------------------------
! reverting damage field state 
   cell = 0
   call SNESGetDM(damage_snes,dm_local,ierr); CHKERRQ(ierr)
   call DMDAVecGetArrayF90(dm_local,solution_vec,x_scal,ierr); CHKERRQ(ierr)                            !< get the data out of PETSc to work with
   x_scal(xstart:xend,ystart:yend,zstart:zend) = damage_current
   call DMDAVecRestoreArrayF90(dm_local,solution_vec,x_scal,ierr); CHKERRQ(ierr)
   do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
     cell = cell + 1                                                                           
     call damage_nonlocal_putNonLocalDamage(damage_current(i,j,k),1,cell)
   enddo; enddo; enddo
 else
!--------------------------------------------------------------------------------------------------
! update rate and forward last inc
   damage_lastInc = damage_current
   cell = 0
   D_ref = 0.0_pReal
   mobility_ref = 0.0_pReal
   do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
     cell = cell + 1
     D_ref = D_ref + damage_nonlocal_getDiffusion33(1,cell)
     mobility_ref = mobility_ref + damage_nonlocal_getMobility(1,cell)
   enddo; enddo; enddo
   D_ref = D_ref*wgt
   call MPI_Allreduce(MPI_IN_PLACE,D_ref,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
   mobility_ref = mobility_ref*wgt
   call MPI_Allreduce(MPI_IN_PLACE,mobility_ref,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 endif  

end subroutine grid_damage_spectral_forward


!--------------------------------------------------------------------------------------------------
!> @brief forms the spectral damage residual vector
!--------------------------------------------------------------------------------------------------
subroutine formResidual(in,x_scal,f_scal,dummy,ierr)
 use numerics, only: &
   residualStiffness
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
 use damage_nonlocal, only: &
   damage_nonlocal_getSourceAndItsTangent,&
   damage_nonlocal_getDiffusion33, &
   damage_nonlocal_getMobility

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
 real(pReal)   :: phiDot, dPhiDot_dPhi, mobility

 damage_current = x_scal 
!--------------------------------------------------------------------------------------------------
! evaluate polarization field
 scalarField_real = 0.0_pReal
 scalarField_real(1:grid(1),1:grid(2),1:grid3) = damage_current 
 call utilities_FFTscalarForward()
 call utilities_fourierScalarGradient()                                                             !< calculate gradient of damage field
 call utilities_FFTvectorBackward()
 cell = 0
 do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
   cell = cell + 1
   vectorField_real(1:3,i,j,k) = math_mul33x3(damage_nonlocal_getDiffusion33(1,cell) - D_ref, &
                                              vectorField_real(1:3,i,j,k))
 enddo; enddo; enddo
 call utilities_FFTvectorForward()
 call utilities_fourierVectorDivergence()                                                           !< calculate damage divergence in fourier field
 call utilities_FFTscalarBackward()
 cell = 0
 do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
   cell = cell + 1
   call damage_nonlocal_getSourceAndItsTangent(phiDot, dPhiDot_dPhi, damage_current(i,j,k), 1, cell)
   mobility = damage_nonlocal_getMobility(1,cell)
   scalarField_real(i,j,k) = params%timeinc*scalarField_real(i,j,k) + &
                             params%timeinc*phiDot + &
                             mobility*damage_lastInc(i,j,k) - &
                             mobility*damage_current(i,j,k) + &
                             mobility_ref*damage_current(i,j,k)
 enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! convolution of damage field with green operator
 call utilities_FFTscalarForward()
 call utilities_fourierGreenConvolution(D_ref, mobility_ref, params%timeinc)
 call utilities_FFTscalarBackward()
 where(scalarField_real(1:grid(1),1:grid(2),1:grid3) > damage_lastInc) &
       scalarField_real(1:grid(1),1:grid(2),1:grid3) = damage_lastInc
 where(scalarField_real(1:grid(1),1:grid(2),1:grid3) < residualStiffness) &
       scalarField_real(1:grid(1),1:grid(2),1:grid3) = residualStiffness
 
!--------------------------------------------------------------------------------------------------
! constructing residual
 f_scal = scalarField_real(1:grid(1),1:grid(2),1:grid3) - damage_current

end subroutine formResidual


end module grid_damage_spectral
