!--------------------------------------------------------------------------------------------------
! $Id: spectral_damage.f90 4082 2015-04-11 20:28:07Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Shaokang Zhang, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Spectral solver for nonlocal damage
!--------------------------------------------------------------------------------------------------
module spectral_damage
 use prec, only: & 
   pInt, &
   pReal
 use math, only: &
   math_I3
 use DAMASK_spectral_Utilities, only: &
   tSolutionState, &
   tSolutionParams
 use numerics, only: &
   worldrank, &
   worldsize

 implicit none
 private
#include <petsc/finclude/petsc.h90>

 character (len=*), parameter, public :: &
   spectral_damage_label = 'spectraldamage'
   
!--------------------------------------------------------------------------------------------------
! derived types
 type(tSolutionParams), private :: params

!--------------------------------------------------------------------------------------------------
! PETSc data
 SNES, private :: damage_snes
 Vec,  private :: solution
 PetscInt, private :: xstart, xend, ystart, yend, zstart, zend
 real(pReal), private, dimension(:,:,:), allocatable :: &
   damage_current, &                                                                           !< field of current damage
   damage_lastInc, &                                                                           !< field of previous damage
   damage_stagInc                                                                              !< field of staggered damage

!--------------------------------------------------------------------------------------------------
! reference diffusion tensor, mobility etc. 
 integer(pInt),               private :: totalIter = 0_pInt                                         !< total iteration in current increment
 real(pReal), dimension(3,3), private :: D_ref
 real(pReal), private                 :: mobility_ref
 character(len=1024),         private :: incInfo
 
 public :: &
   spectral_damage_init, &
   spectral_damage_solution, &
   spectral_damage_forward, &
   spectral_damage_destroy
 external :: &
   VecDestroy, &
   DMDestroy, &
   DMDACreate3D, &
   DMCreateGlobalVector, &
   DMDASNESSetFunctionLocal, &
   PETScFinalize, &
   SNESDestroy, &
   SNESGetNumberFunctionEvals, &
   SNESGetIterationNumber, &
   SNESSolve, &
   SNESSetDM, &
   SNESGetConvergedReason, &
   SNESSetConvergenceTest, &
   SNESSetFromOptions, &
   SNESCreate, &
   MPI_Abort, &
   MPI_Bcast, &
   MPI_Allreduce

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
subroutine spectral_damage_init()
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 use IO, only: &
   IO_intOut, &
   IO_read_realFile, &
   IO_timeStamp
 use DAMASK_spectral_Utilities, only: &
   wgt
 use mesh, only: &
   grid, &
   grid3
 use damage_nonlocal, only: &
   damage_nonlocal_getDiffusion33, &
   damage_nonlocal_getMobility
   
 implicit none
 DM :: damage_grid
 Vec :: uBound, lBound
 PetscErrorCode :: ierr
 PetscObject    :: dummy
 integer(pInt), dimension(:), allocatable :: localK  
 integer(pInt) :: proc
 integer(pInt) :: i, j, k, cell
 character(len=100) :: snes_type

 mainProcess: if (worldrank == 0_pInt) then
   write(6,'(/,a)') ' <<<+-  spectral_damage init  -+>>>'
   write(6,'(a)') ' $Id: spectral_damage.f90 4082 2015-04-11 20:28:07Z MPIE\m.diehl $'
   write(6,'(a15,a)')   ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
 call SNESCreate(PETSC_COMM_WORLD,damage_snes,ierr); CHKERRQ(ierr)
 call SNESSetOptionsPrefix(damage_snes,'damage_',ierr);CHKERRQ(ierr) 
 allocate(localK(worldsize), source = 0); localK(worldrank+1) = grid3
 do proc = 1, worldsize
   call MPI_Bcast(localK(proc),1,MPI_INTEGER,proc-1,PETSC_COMM_WORLD,ierr)
 enddo  
 call DMDACreate3d(PETSC_COMM_WORLD, &
        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &                                     !< cut off stencil at boundary
        DMDA_STENCIL_BOX, &                                                                         !< Moore (26) neighborhood around central point
        grid(1),grid(2),grid(3), &                                                                  !< global grid
        1, 1, worldsize, &
        1, 0, &                                                                                     !< #dof (damage phase field), ghost boundary width (domain overlap)
        grid(1),grid(2),localK, &                                                                   !< local grid
        damage_grid,ierr)                                                                           !< handle, error
 CHKERRQ(ierr)
 call SNESSetDM(damage_snes,damage_grid,ierr); CHKERRQ(ierr)                                        !< connect snes to da
 call DMCreateGlobalVector(damage_grid,solution,ierr); CHKERRQ(ierr)                                !< global solution vector (grid x 1, i.e. every def grad tensor)
 call DMDASNESSetFunctionLocal(damage_grid,INSERT_VALUES,spectral_damage_formResidual,dummy,ierr)   !< residual vector of same shape as solution vector
 CHKERRQ(ierr) 
 call SNESSetFromOptions(damage_snes,ierr); CHKERRQ(ierr)                                           !< pull it all together with additional cli arguments
 call SNESGetType(damage_snes,snes_type,ierr); CHKERRQ(ierr)
 if (trim(snes_type) == 'vinewtonrsls' .or. &
     trim(snes_type) == 'vinewtonssls') then
   call DMGetGlobalVector(damage_grid,lBound,ierr); CHKERRQ(ierr)
   call DMGetGlobalVector(damage_grid,uBound,ierr); CHKERRQ(ierr)
   call VecSet(lBound,0.0,ierr); CHKERRQ(ierr)
   call VecSet(uBound,1.0,ierr); CHKERRQ(ierr)
   call SNESVISetVariableBounds(damage_snes,lBound,uBound,ierr)                                        !< variable bounds for variational inequalities like contact mechanics, damage etc.
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
 call VecSet(solution,1.0,ierr); CHKERRQ(ierr)
 allocate(damage_current(grid(1),grid(2),grid3), source=1.0_pReal)
 allocate(damage_lastInc(grid(1),grid(2),grid3), source=1.0_pReal)
 allocate(damage_stagInc(grid(1),grid(2),grid3), source=1.0_pReal)

!--------------------------------------------------------------------------------------------------
! damage reference diffusion update
 cell = 0_pInt
 D_ref = 0.0_pReal
 mobility_ref = 0.0_pReal
 do k = 1_pInt, grid3;  do j = 1_pInt, grid(2);  do i = 1_pInt,grid(1)
   cell = cell + 1_pInt
   D_ref = D_ref + damage_nonlocal_getDiffusion33(1,cell)
   mobility_ref = mobility_ref + damage_nonlocal_getMobility(1,cell)
 enddo; enddo; enddo
 D_ref = D_ref*wgt
 call MPI_Allreduce(MPI_IN_PLACE,D_ref,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 mobility_ref = mobility_ref*wgt
 call MPI_Allreduce(MPI_IN_PLACE,mobility_ref,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)

end subroutine spectral_damage_init
  
!--------------------------------------------------------------------------------------------------
!> @brief solution for the spectral damage scheme with internal iterations
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function spectral_damage_solution(guess,timeinc,timeinc_old,loadCaseTime)
 use numerics, only: &
   itmax, &
   err_damage_tolAbs, &
   err_damage_tolRel
 use DAMASK_spectral_Utilities, only: &
   tBoundaryCondition, &
   Utilities_maskedCompliance, &
   Utilities_updateGamma
 use mesh, only: &
   grid, &
   grid3
 use damage_nonlocal, only: &
   damage_nonlocal_putNonLocalDamage

 implicit none

!--------------------------------------------------------------------------------------------------
! input data for solution
 real(pReal), intent(in) :: &
   timeinc, &                                                                                       !< increment in time for current solution
   timeinc_old, &                                                                                   !< increment in time of last increment
   loadCaseTime                                                                                     !< remaining time of current load case
 logical, intent(in) :: guess
 integer(pInt) :: i, j, k, cell
 PetscInt  ::position
 PetscReal ::  minDamage, maxDamage, stagNorm, solnNorm
 
!--------------------------------------------------------------------------------------------------
! PETSc Data
 PetscErrorCode :: ierr   
 SNESConvergedReason :: reason
 
 spectral_damage_solution%converged =.false.
 
!--------------------------------------------------------------------------------------------------
! set module wide availabe data 
 params%timeinc = timeinc
 params%timeincOld = timeinc_old

 call SNESSolve(damage_snes,PETSC_NULL_OBJECT,solution,ierr); CHKERRQ(ierr)
 call SNESGetConvergedReason(damage_snes,reason,ierr); CHKERRQ(ierr)

 if (reason < 1) then
   spectral_damage_solution%converged = .false.
   spectral_damage_solution%iterationsNeeded = itmax
 else
   spectral_damage_solution%converged = .true.
   spectral_damage_solution%iterationsNeeded = totalIter
 endif
 stagNorm = maxval(abs(damage_current - damage_stagInc))
 solnNorm = maxval(abs(damage_current))
 call MPI_Allreduce(MPI_IN_PLACE,stagNorm,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,solnNorm,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
 damage_stagInc = damage_current
 spectral_damage_solution%stagConverged =     stagNorm < err_damage_tolAbs &
                                         .or. stagNorm < err_damage_tolRel*solnNorm

!--------------------------------------------------------------------------------------------------
! updating damage state 
 cell = 0_pInt                                                                                      !< material point = 0
 do k = 1_pInt, grid3;  do j = 1_pInt, grid(2);  do i = 1_pInt,grid(1)
   cell = cell + 1_pInt                                                                             !< material point increase
   call damage_nonlocal_putNonLocalDamage(damage_current(i,j,k),1,cell)
 enddo; enddo; enddo

 call VecMin(solution,position,minDamage,ierr); CHKERRQ(ierr)
 call VecMax(solution,position,maxDamage,ierr); CHKERRQ(ierr)
 if (worldrank == 0) then 
   if (spectral_damage_solution%converged) &
     write(6,'(/,a)') ' ... nonlocal damage converged .....................................'
   write(6,'(/,a,f8.6,2x,f8.6,2x,f8.6,/)',advance='no') ' Minimum|Maximum|Delta Damage      = ',&
                                                       minDamage, maxDamage, stagNorm
   write(6,'(/,a)') ' ==========================================================================='
   flush(6) 
 endif 

end function spectral_damage_solution


!--------------------------------------------------------------------------------------------------
!> @brief forms the spectral damage residual vector
!--------------------------------------------------------------------------------------------------
subroutine spectral_damage_formResidual(in,x_scal,f_scal,dummy,ierr)
 use numerics, only: &
   residualStiffness
 use mesh, only: &
   grid, &
   grid3
 use math, only: &
   math_mul33x3
 use DAMASK_spectral_Utilities, only: &
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
   XG_RANGE,YG_RANGE,ZG_RANGE) :: &
   x_scal
 PetscScalar, dimension( &
   X_RANGE,Y_RANGE,Z_RANGE) :: &
   f_scal
 PetscObject :: dummy
 PetscErrorCode :: ierr
 integer(pInt) :: i, j, k, cell
 real(pReal)   :: phiDot, dPhiDot_dPhi, mobility

 damage_current = x_scal 
!--------------------------------------------------------------------------------------------------
! evaluate polarization field
 scalarField_real = 0.0_pReal
 scalarField_real(1:grid(1),1:grid(2),1:grid3) = damage_current 
 call utilities_FFTscalarForward()
 call utilities_fourierScalarGradient()                                                             !< calculate gradient of damage field
 call utilities_FFTvectorBackward()
 cell = 0_pInt
 do k = 1_pInt, grid3;  do j = 1_pInt, grid(2);  do i = 1_pInt,grid(1)
   cell = cell + 1_pInt
   vectorField_real(1:3,i,j,k) = math_mul33x3(damage_nonlocal_getDiffusion33(1,cell) - D_ref, &
                                              vectorField_real(1:3,i,j,k))
 enddo; enddo; enddo
 call utilities_FFTvectorForward()
 call utilities_fourierVectorDivergence()                                                           !< calculate damage divergence in fourier field
 call utilities_FFTscalarBackward()
 cell = 0_pInt
 do k = 1_pInt, grid3;  do j = 1_pInt, grid(2);  do i = 1_pInt,grid(1)
   cell = cell + 1_pInt
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

end subroutine spectral_damage_formResidual

!--------------------------------------------------------------------------------------------------
!> @brief spectral damage forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine spectral_damage_forward(guess,timeinc,timeinc_old,loadCaseTime)
 use mesh, only: &
   grid, &
   grid3
 use DAMASK_spectral_Utilities, only: &
   cutBack, &
   wgt
 use damage_nonlocal, only: &
   damage_nonlocal_putNonLocalDamage, &
   damage_nonlocal_getDiffusion33, &
   damage_nonlocal_getMobility
   
 implicit none
 real(pReal), intent(in) :: &
   timeinc_old, &
   timeinc, &
   loadCaseTime                                                                                     !< remaining time of current load case
 logical,     intent(in)                     :: guess
 PetscErrorCode                              :: ierr
 integer(pInt)                               :: i, j, k, cell
 DM :: dm_local
 PetscScalar,  dimension(:,:,:), pointer     :: x_scal

 if (cutBack) then 
   damage_current = damage_lastInc
   damage_stagInc = damage_lastInc
!--------------------------------------------------------------------------------------------------
! reverting damage field state 
   cell = 0_pInt
   call SNESGetDM(damage_snes,dm_local,ierr); CHKERRQ(ierr)
   call DMDAVecGetArrayF90(dm_local,solution,x_scal,ierr); CHKERRQ(ierr)                            !< get the data out of PETSc to work with
   x_scal(xstart:xend,ystart:yend,zstart:zend) = damage_current
   call DMDAVecRestoreArrayF90(dm_local,solution,x_scal,ierr); CHKERRQ(ierr)
   do k = 1_pInt, grid3;  do j = 1_pInt, grid(2);  do i = 1_pInt,grid(1)
     cell = cell + 1_pInt                                                                           
     call damage_nonlocal_putNonLocalDamage(damage_current(i,j,k),1,cell)
   enddo; enddo; enddo
 else
!--------------------------------------------------------------------------------------------------
! update rate and forward last inc
   damage_lastInc = damage_current
   cell = 0_pInt
   D_ref = 0.0_pReal
   mobility_ref = 0.0_pReal
   do k = 1_pInt, grid3;  do j = 1_pInt, grid(2);  do i = 1_pInt,grid(1)
     cell = cell + 1_pInt
     D_ref = D_ref + damage_nonlocal_getDiffusion33(1,cell)
     mobility_ref = mobility_ref + damage_nonlocal_getMobility(1,cell)
   enddo; enddo; enddo
   D_ref = D_ref*wgt
   call MPI_Allreduce(MPI_IN_PLACE,D_ref,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
   mobility_ref = mobility_ref*wgt
   call MPI_Allreduce(MPI_IN_PLACE,mobility_ref,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 endif  

 end subroutine spectral_damage_forward

!--------------------------------------------------------------------------------------------------
!> @brief destroy routine
!--------------------------------------------------------------------------------------------------
subroutine spectral_damage_destroy()

 implicit none
 PetscErrorCode :: ierr

 call VecDestroy(solution,ierr); CHKERRQ(ierr)
 call SNESDestroy(damage_snes,ierr); CHKERRQ(ierr)

end subroutine spectral_damage_destroy

end module spectral_damage
