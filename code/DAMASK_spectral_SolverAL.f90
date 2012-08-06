!--------------------------------------------------------------------------------------------------
! $Id: DAMASK_spectral_SolverAL.f90 1654 2012-08-03 09:25:48Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief AL scheme solver
!--------------------------------------------------------------------------------------------------
module DAMASK_spectral_SolverAL
 
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 
 use prec, only: & 
   pInt, &
   pReal
 
 use math, only: &
   math_I3
 
 use DAMASK_spectral_Utilities, only: &
   solutionState
 
 implicit none
 character (len=*), parameter, public :: &
   DAMASK_spectral_SolverAL_label = 'AL'
   
!--------------------------------------------------------------------------------------------------
! derived types
 type solutionParams 
   real(pReal), dimension(3,3) :: P_BC, rotation_BC
   real(pReal) :: timeinc
 end type solutionParams
 
 type(solutionParams), private :: params

!--------------------------------------------------------------------------------------------------
! PETSc data
  SNES, private :: snes
  DM, private :: da
  Vec, private :: x,r
  PetscMPIInt, private :: rank
  integer(pInt), private :: iter
  PetscInt, private :: xs,xm,gxs,gxm
  PetscInt, private :: ys,ym,gys,gym
  PetscInt, private :: zs,zm,gzs,gzm
  character(len=1024), private :: PetSc_options = '-snes_type ngmres -snes_ngmres_anderson -snes_monitor -snes_view'
 
!--------------------------------------------------------------------------------------------------
! common pointwise data
  real(pReal), private, dimension(:,:,:,:,:), allocatable ::  F, F_lastInc, F_lambda, F_lambda_lastInc, P
  real(pReal), private, dimension(:,:,:,:),   allocatable ::  coordinates
  real(pReal), private, dimension(:,:,:),     allocatable ::  temperature
 
!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
  real(pReal), private, dimension(3,3) :: &
    F_aim = math_I3, &
    F_aim_lastInc = math_I3, &
    P_av
   
  real(pReal), private, dimension(3,3,3,3) :: &
    C = 0.0_pReal, &
    S = 0.0_pReal, &
    C_scale = 0.0_pReal, &
    S_scale = 0.0_pReal
 
  real(pReal), private :: err_stress, err_f, err_p
  logical, private :: ForwardData
  real(pReal), private, dimension(3,3) :: &
    mask_stress = 0.0_pReal
 
  contains
 
!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
  subroutine AL_init()
      
    use IO, only: &
      IO_read_JobBinaryFile, &
      IO_write_JobBinaryFile
    
    use FEsolving, only: &
      restartInc
   
    use DAMASK_interface, only: &
      getSolverJobName
        
    use DAMASK_spectral_Utilities, only: &
      Utilities_init, &
      Utilities_constitutiveResponse, &
      Utilities_updateGamma, &
      debugrestart
         
    use mesh, only: &
      res, &
      geomdim
     
    implicit none
    integer(pInt) :: i,j,k
    real(pReal), dimension(3,3) :: temp33_Real
    PetscErrorCode  ierr_psc
   
    call Utilities_init()
    
    write(6,'(a)') ''
    write(6,'(a)') ' <<<+-  DAMASK_spectral_solverAL init  -+>>>'
    write(6,'(a)') ' $Id: DAMASK_spectral_SolverAL.f90 1654 2012-08-03 09:25:48Z MPIE\m.diehl $'
   #include "compilation_info.f90"
    write(6,'(a)') ''
   
    allocate (F          (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
    allocate (F_lastInc  (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
    allocate (F_lambda   (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
    allocate (F_lambda_lastInc(res(1),res(2),res(3),3,3),  source = 0.0_pReal)
    allocate (P          (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
    allocate (coordinates(  res(1),  res(2),res(3),3),    source = 0.0_pReal)
    allocate (temperature(  res(1),  res(2),res(3)),      source = 0.0_pReal)
   
 !--------------------------------------------------------------------------------------------------
 ! init fields                 
    if (restartInc == 1_pInt) then                                                                     ! no deformation (no restart)
      do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
        F(i,j,k,1:3,1:3) = math_I3
        F_lastInc(i,j,k,1:3,1:3) = math_I3
        F_lambda(i,j,k,1:3,1:3) = math_I3
        F_lambda_lastInc(i,j,k,1:3,1:3) = math_I3
        coordinates(i,j,k,1:3) = geomdim/real(res,pReal)*real([i,j,k],pReal) &
                               - geomdim/real(2_pInt*res,pReal)
      enddo; enddo; enddo
    elseif (restartInc > 1_pInt) then                                                                  ! using old values from file                                                      
      if (debugRestart) write(6,'(a,i6,a)') 'Reading values of increment ',&
                                                restartInc - 1_pInt,' from file' 
      call IO_read_jobBinaryFile(777,'convergedSpectralDefgrad',&
                                                   trim(getSolverJobName()),size(F))
      read (777,rec=1) F
      close (777)
      call IO_read_jobBinaryFile(777,'convergedSpectralDefgrad_lastInc',&
                                                   trim(getSolverJobName()),size(F_lastInc))
      read (777,rec=1) F_lastInc
      close (777)
      call IO_read_jobBinaryFile(777,'convergedSpectralDefgradLambda',&
                                                   trim(getSolverJobName()),size(F_lambda))
      read (777,rec=1) F
      close (777)
      call IO_read_jobBinaryFile(777,'convergedSpectralDefgradLambda_lastInc',&
                                                   trim(getSolverJobName()),size(F_lambda_lastInc))
      read (777,rec=1) F_lastInc
      close (777)
      call IO_read_jobBinaryFile(777,'F_aim',trim(getSolverJobName()),size(F_aim))
      read (777,rec=1) F_aim
      close (777)
      call IO_read_jobBinaryFile(777,'F_aim_lastInc',trim(getSolverJobName()),size(F_aim_lastInc))
      read (777,rec=1) F_aim_lastInc
      close (777)
  
      coordinates = 0.0 ! change it later!!!
    endif
   
    call constitutiveResponse(coordinates,F,F_lastInc,temperature,0.0_pReal,&
                                 P,C,P_av,.false.,math_I3)
   
 !--------------------------------------------------------------------------------------------------
 ! reference stiffness
    if (restartInc == 1_pInt) then
      call IO_write_jobBinaryFile(777,'C_ref',size(C))
      write (777,rec=1) C
      close(777)
    elseif (restartInc > 1_pInt) then
      call IO_read_jobBinaryFile(777,'C_ref',trim(getSolverJobName()),size(C))
      read (777,rec=1) C
      close (777)
    endif
   
    call Utilities_updateGamma(C)
    C_scale = C
    S_scale = math_invSym3333(C)
   
 !--------------------------------------------------------------------------------------------------
 ! PETSc Init
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr_psc)
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr_psc)
   
    call SNESCreate(PETSC_COMM_WORLD,snes,ierr_psc)
    call DMDACreate3d(PETSC_COMM_WORLD,                               &
              DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, &
              DMDA_STENCIL_BOX,res(1),res(2),res(3),PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, &
              18,1,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da,ierr_psc)
    call DMCreateGlobalVector(da,x,ierr_psc)
    call VecDuplicate(x,r,ierr_psc)
    call DMDASetLocalFunction(da,AL_FormRHS,ierr_psc)
  
    call SNESSetDM(snes,da,ierr_psc)
    call SNESSetFunction(snes,r,SNESDMDAComputeFunction,da,ierr_psc)
    call SNESSetConvergenceTest(snes,AL_converged,dummy,PETSC_NULL_FUNCTION,ierr_psc)
    call PetscOptionsInsertString(PetSc_options,ierr_psc)
    call SNESSetFromOptions(snes,ierr_psc)  
    call DMDAGetCorners(da,xs,ys,zs,xm,ym,zm,ierr_psc)
    call DMDAGetCorners(da,gxs,gys,gzs,gxm,gym,gzm,ierr_psc)
  
    xs = xs+1; gxs = gxs+1 
    xm = xm-1; gxm = gxm-1 
    ys = ys+1; gys = gys+1 
    ym = ym-1; gym = gym-1
    zs = zs+1; gzs = gzs+1 
    zm = zm-1; gzm = gzm-1
 
  end subroutine AL_init
  
!--------------------------------------------------------------------------------------------------
!> @brief solution for the AL scheme with internal iterations
!--------------------------------------------------------------------------------------------------
  type(solutionState) function AL_solution(guessmode,timeinc,timeinc_old,P_BC,F_BC,temperature_bc,rotation_BC)
   
   use numerics, only: &
     update_gamma
   use math, only: &
     math_mul33x33 ,&
     math_rotate_backward33, &
     deformed_fft
   use mesh, only: &
     res,&
     geomdim
   use IO, only: &
     IO_write_JobBinaryFile
     
   use DAMASK_spectral_Utilities, only: &
     boundaryCondition, &
     Utilities_forwardField, &
     Utilities_maskedCompliance, &
     Utilities_updateGamma
       
   use FEsolving, only: &
     restartWrite
   
   implicit none
!--------------------------------------------------------------------------------------------------
! input data for solution
   real(pReal), intent(in) :: timeinc, timeinc_old, temperature_bc, guessmode
   type(boundaryCondition),      intent(in) :: P_BC,F_BC
   real(pReal), dimension(3,3), intent(in) :: rotation_BC
   
  
   
   real(pReal), dimension(3,3)            :: deltaF_aim, &
                                             F_aim_lab, &
                                             F_aim_lab_lastIter
!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
   real(pReal)   :: err_div, err_stress       
   integer(pInt) :: iter, row, column, i, j, k
   real(pReal)   :: defgradDet, defgradDetMax, defgradDetMin
   real(pReal), dimension(3,3)            :: temp33_Real 
 
!--------------------------------------------------------------------------------------------------
! 
   PetscScalar, pointer :: xx_psc(:)
   PetscErrorCode ierr_psc
   
!--------------------------------------------------------------------------------------------------
! restart information for spectral solver
   if (restartWrite) then
     write(6,'(a)') 'writing converged results for restart'
     call IO_write_jobBinaryFile(777,'convergedSpectralDefgrad',size(F_lastInc))
     write (777,rec=1) F_LastInc
     close (777)
     call IO_write_jobBinaryFile(777,'C',size(C))
     write (777,rec=1) C
     close(777)
   endif 
  
!--------------------------------------------------------------------------------------------------
! winding forward of deformation aim in loadcase system
   if (F_BC%myType=='l') then                                                        ! calculate deltaF_aim from given L and current F
     deltaF_aim = timeinc * F_BC%maskFloat * math_mul33x33(F_BC%values, F_aim)
   elseif(F_BC%myType=='fdot')   then                                                                                      ! deltaF_aim = fDot *timeinc where applicable
     deltaF_aim = timeinc * F_BC%maskFloat * F_BC%values
   endif
   temp33_Real = F_aim                                            
   F_aim = F_aim &                                                                         
           + guessmode * P_BC%maskFloat * (F_aim - F_aim_lastInc)*timeinc/timeinc_old &      
           + deltaF_aim
   F_aim_lastInc = temp33_Real
   F_aim_lab = math_rotate_backward33(F_aim,rotation_BC)                            ! boundary conditions from load frame into lab (Fourier) frame
  
!--------------------------------------------------------------------------------------------------
! update local deformation gradient and coordinates
   deltaF_aim = math_rotate_backward33(deltaF_aim,rotation_BC)
   call Utilities_forwardField(deltaF_aim,timeinc,timeinc_old,guessmode,F_lastInc,F)
   call Utilities_forwardField(deltaF_aim,timeinc,timeinc_old,guessmode,F_lambda_lastInc,F_lambda)
   call deformed_fft(res,geomdim,math_rotate_backward33(F_aim,rotation_BC),1.0_pReal,F_lastInc,coordinates)
  
!--------------------------------------------------------------------------------------------------
! update stiffness (and gamma operator)
   S = Utilities_maskedCompliance(rotation_BC,P_BC%maskLogical,C)
   if (update_gamma) call Utilities_updateGamma(C)
   
   iter = 0_pInt
   ForwardData = .True.
   mask_stress = P_BC%maskFloat
   params%P_BC = P_BC%values
   params%rotation_BC = rotation_BC
   params%timeinc = timeinc
   
   call VecGetArrayF90(x,xx_psc,ierr_psc)
   call AL_InitialGuess(xx_psc)
   call VecRestoreArrayF90(x,xx_psc,ierr_psc) 
   call SNESSolve(snes,PETSC_NULL_OBJECT,x,ierr_psc)
 
 end function AL_solution

! ------------------------------------------------------------------- 

 subroutine AL_InitialGuess(xx_psc)

   implicit none
 #include <finclude/petsc.h>
  
!  Input/output variables:

   PetscScalar xx_psc(0:17,gxs:(gxs+gxm),gys:(gys+gym),gxs:(gzs+gzm))
   integer(pInt) :: i, j, k

!  Compute function over the locally owned part of the grid


   do k=gzs,gzs+gzm; do j=gys,gys+gym; do i=gxs,gxs+gxm     
     xx_psc(0,i,j,k) = F(i,j,k,1,1)
     xx_psc(1,i,j,k) = F(i,j,k,1,2)
     xx_psc(2,i,j,k) = F(i,j,k,1,3)
     xx_psc(3,i,j,k) = F(i,j,k,2,1)
     xx_psc(4,i,j,k) = F(i,j,k,2,2)
     xx_psc(5,i,j,k) = F(i,j,k,2,3)
     xx_psc(6,i,j,k) = F(i,j,k,3,1)
     xx_psc(7,i,j,k) = F(i,j,k,3,2)
     xx_psc(8,i,j,k) = F(i,j,k,3,3)
     xx_psc(9,i,j,k) =  F_lambda(i,j,k,1,1)
     xx_psc(10,i,j,k) = F_lambda(i,j,k,1,2)
     xx_psc(11,i,j,k) = F_lambda(i,j,k,1,3)
     xx_psc(12,i,j,k) = F_lambda(i,j,k,2,1)
     xx_psc(13,i,j,k) = F_lambda(i,j,k,2,2)
     xx_psc(14,i,j,k) = F_lambda(i,j,k,2,3)
     xx_psc(15,i,j,k) = F_lambda(i,j,k,3,1)
     xx_psc(16,i,j,k) = F_lambda(i,j,k,3,2)
     xx_psc(17,i,j,k) = F_lambda(i,j,k,3,3)
   enddo; enddo; enddo

   return
 end subroutine AL_InitialGuess

! ---------------------------------------------------------------------
!
!  Input Parameter:
!  x - local vector data
!
!  Output Parameters:
!  f - local vector data, f(x)
!  ierr - error code 
!
!  Notes:
!  This routine uses standard Fortran-style computations over a 3-dim array.
!
 subroutine AL_FormRHS(in,x_scal,f_scal,dummy,ierr_psc)
  
   use numerics, only: &
     itmax, &
     itmin
   use math, only: &
     math_rotate_backward33, &
     math_transpose33, &
     math_mul3333xx33
   use mesh, only: &
     res
   use DAMASK_spectral_Utilities, only: &
     field_real, &
     Utilities_forwardFFT, &
     Utilities_fourierConvolution, &
     Utilities_backwardFFT, &
     Utilities_constitutiveResponse
  
   implicit none
 #include <finclude/petsc.h>

   integer(pInt) :: i,j,k
   Input/output variables:
   DMDALocalInfo in(DMDA_LOCAL_INFO_SIZE)
   PetscScalar x_scal(0:17,XG_RANGE,YG_RANGE,ZG_RANGE)  
   PetscScalar f_scal(0:17,X_RANGE,Y_RANGE,Z_RANGE)
   real(pReal), dimension (3,3) ::  temp33_real
   PetscObject dummy
   PetscErrorCode ierr_psc

   iter = iter + 1_pInt
  
 !--------------------------------------------------------------------------------------------------
 ! report begin of new iteration
   write(6,'(a)') ''
   write(6,'(a)') '=================================================================='
   write(6,'(3(a,i6.6))') ' @ Iter. ',itmin,' < ',iter,' < ',itmax
   write(6,'(a,/,3(3(f12.7,1x)/))',advance='no') 'deformation gradient aim =',&
                                                             math_transpose33(F_aim)
  
   do k=gzs,gzs+gzm; do j=gys,gys+gym; do i=gxs,gxs+gxm
     F(i,j,k,1,1) = x_scal(0,i,j,k)
     F(i,j,k,1,2) = x_scal(1,i,j,k)
     F(i,j,k,1,3) = x_scal(2,i,j,k)
     F(i,j,k,2,1) = x_scal(3,i,j,k)
     F(i,j,k,2,2) = x_scal(4,i,j,k)
     F(i,j,k,2,3) = x_scal(5,i,j,k)
     F(i,j,k,3,1) = x_scal(6,i,j,k)
     F(i,j,k,3,2) = x_scal(7,i,j,k)
     F(i,j,k,3,3) = x_scal(8,i,j,k)
     F_lambda(i,j,k,1,1) = x_scal(9,i,j,k)
     F_lambda(i,j,k,1,2) = x_scal(10,i,j,k)
     F_lambda(i,j,k,1,3) = x_scal(11,i,j,k)
     F_lambda(i,j,k,2,1) = x_scal(12,i,j,k)
     F_lambda(i,j,k,2,2) = x_scal(13,i,j,k)
     F_lambda(i,j,k,2,3) = x_scal(14,i,j,k)
     F_lambda(i,j,k,3,1) = x_scal(15,i,j,k)
     F_lambda(i,j,k,3,2) = x_scal(16,i,j,k)
     F_lambda(i,j,k,3,3) = x_scal(17,i,j,k)
   enddo; enddo; enddo
  
 !--------------------------------------------------------------------------------------------------
 ! evaluate constitutive response
   call constitutiveResponse(coordinates,F,F_lastInc,temperature,params%timeinc,&
                                 P,C,P_av,ForwardData,params%rotation_BC)
   ForwardData = .False.
   
!--------------------------------------------------------------------------------------------------
! stress BC handling
   F_aim = F_aim - math_mul3333xx33(S, ((P_av - params%P_BC))) !S = 0.0 for no bc
   err_stress = maxval(mask_stress * (P_av - params%P_BC))     ! mask = 0.0 for no bc

                                
   F_aim_lab = math_rotate_backward33(F_aim,params%rotation_BC)                            ! boundary conditions from load frame into lab (Fourier) frame
   
 !--------------------------------------------------------------------------------------------------
 ! doing Fourier transform
   field_real = 0.0_pReal
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
     field_real(i,j,k,1:3,1:3) = math_mul3333xx33(C_scale,F_lambda(i,j,k,1:3,1:3)-F(i,j,k,1:3,1:3))
   enddo; enddo; enddo
  
   call Utilities_forwardFFT()
   call Utilities_fourierConvolution(F_aim_lab) 
   call Utilities_backwardFFT()

   err_f = 0.0_pReal
   err_p = 0.0_pReal
  
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
     temp33_real = field_real(i,j,k,1:3,1:3) - F(i,j,k,1:3,1:3)
     err_f = err_f + sum(temp33_real*temp33_real)
    
     temp33_real = F_lambda(i,j,k,1:3,1:3) - &
                   math_mul3333xx33(S_scale,P(i,j,k,1:3,1:3)) + math_I3
     err_p = err_p + sum(temp33_real*temp33_real)
   enddo; enddo; enddo
                                       
   err_f = wgt*sqrt(err_f)/sum((F_aim-math_I3)*(F_aim-math_I3)))
   err_p = wgt*sqrt(err_p)/sum((F_aim-math_I3)*(F_aim-math_I3)))
  
   do k=zs,ze; do j=ys,ye; do i=xs,xe
     temp33_real = math_mul3333xx33(S_scale,P(i,j,k,1:3,1:3)) + math_I3 - F_lambda(i,j,k,1:3,1:3) &
            + F(i,j,k,1:3,1:3) - field_real(i,j,k,1:3,1:3)
     f_scal(0,i,j,k) = temp33_real(1,1)
     f_scal(1,i,j,k) = temp33_real(1,2)
     f_scal(2,i,j,k) = temp33_real(1,3)
     f_scal(3,i,j,k) = temp33_real(2,1)
     f_scal(4,i,j,k) = temp33_real(2,2)
     f_scal(5,i,j,k) = temp33_real(2,3)
     f_scal(6,i,j,k) = temp33_real(3,1)
     f_scal(7,i,j,k) = temp33_real(3,2)
     f_scal(8,i,j,k) = temp33_real(3,3)
     f_scal(9,i,j,k)  = F(i,j,k,1,1) - field_real(i,j,k,1,1)
     f_scal(10,i,j,k) = F(i,j,k,1,2) - field_real(i,j,k,1,2)
     f_scal(11,i,j,k) = F(i,j,k,1,3) - field_real(i,j,k,1,3)
     f_scal(12,i,j,k) = F(i,j,k,2,1) - field_real(i,j,k,2,1)
     f_scal(13,i,j,k) = F(i,j,k,2,2) - field_real(i,j,k,2,2)
     f_scal(14,i,j,k) = F(i,j,k,2,3) - field_real(i,j,k,2,3)
     f_scal(15,i,j,k) = F(i,j,k,3,1) - field_real(i,j,k,3,1)
     f_scal(16,i,j,k) = F(i,j,k,3,2) - field_real(i,j,k,3,2)
     f_scal(17,i,j,k) = F(i,j,k,3,3) - field_real(i,j,k,3,3)
   enddo; enddo; enddo

   return
 end subroutine AL_FormRHS

 ! ---------------------------------------------------------------------
 ! User defined convergence check
 !
 subroutine AL_converged(snes,it,xnorm,snorm,fnorm,reason,dummy,ierr_psc)
  
   use numerics, only: &
    itmax, &
    itmin, &
    err_f_tol, &
    err_p_tol, &
    err_stress_tolrel, &
    err_stress_tolabs
   
   implicit none
 #include <finclude/petsc.h>

 ! Input/output variables:
   SNES snes
   PetscInt it
   PetscReal xnorm, snorm, fnorm
   SNESConvergedReason reason
   PetscObject dummy
   PetscErrorCode ierr_psc
   logical :: Converged
             
   Converged = (iter < itmax) .and. (iter > itmin) .and. &
                 all([ err_f/sqrt(sum((F_aim-math_I3)*(F_aim-math_I3)))/err_f_tol, &
                       err_p/sqrt(sum((F_aim-math_I3)*(F_aim-math_I3)))/err_p_tol, &
                       err_stress/min(maxval(abs(P_av))*err_stress_tolrel,err_stress_tolabs)] < 1.0_pReal)
   
   if (Converged) then
     reason = 1
   else  
     reason = 0
   endif 
   
   write(6,'(a,es14.7)') 'error stress BC = ', err_stress/min(maxval(abs(P_av))*err_stress_tolrel,err_stress_tolabs)  
   write(6,'(a,es14.7)') 'error F         = ', err_f/sqrt(sum((F_aim-math_I3)*(F_aim-math_I3)))/err_f_tol
   write(6,'(a,es14.7)') 'error P         = ', err_p/sqrt(sum((F_aim-math_I3)*(F_aim-math_I3)))/err_p_tol
   return
 end subroutine AL_converged

 subroutine AL_destroy()

 implicit none

 call VecDestroy(x,ierr_psc)
 call VecDestroy(r,ierr_psc)
 call SNESDestroy(snes,ierr_psc)
 call DMDestroy(da,ierr_psc)
 call PetscFinalize(ierr_psc)
 call Utilities_destroy()

 end subroutine AL_destroy

 end module DAMASK_spectral_SolverAL
