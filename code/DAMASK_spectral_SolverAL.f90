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
#ifdef PETSC
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscdmda.h>
#include <finclude/petscis.h>
#include <finclude/petscmat.h>
#include <finclude/petscksp.h>
#include <finclude/petscpc.h>
#include <finclude/petscsnes.h>
#include <finclude/petscvec.h90>
#include <finclude/petscsnes.h90>
#endif

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
  Vec, private :: solution_vec
 
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
      
    use numerics, only: &
      petsc_options  
         
    use mesh, only: &
      res, &
      geomdim
      
    use math, only: &
      math_invSym3333
      
    implicit none
    
    integer(pInt) :: i,j,k
    
    PetscErrorCode :: ierr_psc
    PetscObject :: dummy
    PetscMPIInt :: rank
    DM :: da
    
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
   
    call Utilities_constitutiveResponse(coordinates,F,F_lastInc,temperature,0.0_pReal,&
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
    call DMCreateGlobalVector(da,solution_vec,ierr_psc)
    call DMDASetLocalFunction(da,AL_formResidual,ierr_psc)
  
    call SNESSetDM(snes,da,ierr_psc)
    call SNESSetConvergenceTest(snes,AL_converged,dummy,PETSC_NULL_FUNCTION,ierr_psc)
    call PetscOptionsInsertString(petsc_options,ierr_psc)
    call SNESSetFromOptions(snes,ierr_psc)  
    call DMDestroy(da,ierr_psc)
 
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
   SNESConvergedReason reason
   real(pReal), dimension(3,3)            :: deltaF_aim, &
                                             F_aim_lab
!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
   real(pReal), dimension(3,3)            :: temp33_Real 
   integer(pInt) :: ctr, i, j, k
 
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
  AL_solution%converged =.false.
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
   
   ForwardData = .True.
   mask_stress = P_BC%maskFloat
   params%P_BC = P_BC%values
   params%rotation_BC = rotation_BC
   params%timeinc = timeinc
   
   ctr = 1_pInt
   call VecGetArrayF90(solution_vec,xx_psc,ierr_psc)
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)    
     xx_psc(ctr:ctr+8_pInt) = reshape(F(i,j,k,1:3,1:3),[9])
     xx_psc(ctr+9_pInt:ctr+17_pInt) = reshape(F_lambda(i,j,k,1:3,1:3),[9])
     ctr = ctr + 18_pInt
   enddo; enddo; enddo
   call VecRestoreArrayF90(solution_vec,xx_psc,ierr_psc) 
   call SNESSolve(snes,PETSC_NULL_OBJECT,solution_vec,ierr_psc)
   call SNESGetConvergedReason(snes,reason,ierr_psc)
   if (reason > 0 ) AL_solution%converged = .true.

 end function AL_solution

!--------------------------------------------------------------------------------------------------
!> @brief forms the AL residual vector
!--------------------------------------------------------------------------------------------------
 subroutine AL_formResidual(in,x_scal,f_scal,dummy,ierr_psc)
  
   use numerics, only: &
     itmax, &
     itmin
   use math, only: &
     math_rotate_backward33, &
     math_transpose33, &
     math_mul3333xx33
   use mesh, only: &
     res, &
     wgt
   use DAMASK_spectral_Utilities, only: &
     field_real, &
     Utilities_forwardFFT, &
     Utilities_fourierConvolution, &
     Utilities_backwardFFT, &
     Utilities_constitutiveResponse
  
   implicit none

   integer(pInt) :: i,j,k
   real(pReal), dimension (3,3) ::  temp33_real
   
   DMDALocalInfo :: in(DMDA_LOCAL_INFO_SIZE)
   PetscScalar :: x_scal(in(DMDA_LOCAL_INFO_DOF),XG_RANGE,YG_RANGE,ZG_RANGE)  
   PetscScalar :: f_scal(in(DMDA_LOCAL_INFO_DOF),X_RANGE,Y_RANGE,Z_RANGE) 
   PetscInt :: iter, nfuncs
   PetscObject :: dummy
   PetscErrorCode :: ierr_psc

   call SNESGetNumberFunctionEvals(snes,nfuncs,ierr_psc)
   call SNESGetIterationNumber(snes,iter,ierr_psc)
  
 !--------------------------------------------------------------------------------------------------
 ! report begin of new iteration
   write(6,'(a)') ''
   write(6,'(a)') '=================================================================='
   write(6,'(4(a,i6.6))') ' @ Iter. ',itmin,' < ',iter,' < ',itmax, ' | # Func. calls = ',nfuncs
   write(6,'(a,/,3(3(f12.7,1x)/))',advance='no') 'deformation gradient aim =',&
                                                             math_transpose33(F_aim)

   do k=in(DMDA_LOCAL_INFO_ZS)+1,in(DMDA_LOCAL_INFO_ZS)+in(DMDA_LOCAL_INFO_ZM) 
     do j=in(DMDA_LOCAL_INFO_YS)+1,in(DMDA_LOCAL_INFO_YS)+in(DMDA_LOCAL_INFO_YM) 
       do i=in(DMDA_LOCAL_INFO_XS)+1,in(DMDA_LOCAL_INFO_XS)+in(DMDA_LOCAL_INFO_XM)
         F(i,j,k,1:3,1:3) = reshape(x_scal(1:9,i,j,k),[3,3])
         F_lambda(i,j,k,1:3,1:3) = reshape(x_scal(10:18,i,j,k),[3,3])
   enddo; enddo; enddo

 !--------------------------------------------------------------------------------------------------
 ! evaluate constitutive response
   call Utilities_constitutiveResponse(coordinates,F_lastInc,F,temperature,params%timeinc,&
                                 P,C,P_av,ForwardData,params%rotation_BC)
   ForwardData = .False.
   
!--------------------------------------------------------------------------------------------------
! stress BC handling
   F_aim = F_aim - math_mul3333xx33(S, ((P_av - params%P_BC))) ! S = 0.0 for no bc
   err_stress = maxval(mask_stress * (P_av - params%P_BC))     ! mask = 0.0 for no bc

   
 !--------------------------------------------------------------------------------------------------
 ! doing Fourier transform
   field_real = 0.0_pReal
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
     field_real(i,j,k,1:3,1:3) = math_mul3333xx33(C_scale,F_lambda(i,j,k,1:3,1:3)-F(i,j,k,1:3,1:3))
   enddo; enddo; enddo
  
   call Utilities_forwardFFT()
   call Utilities_fourierConvolution(math_rotate_backward33(F_aim,params%rotation_BC)) 
   call Utilities_backwardFFT()

   err_f = 0.0_pReal
   err_p = 0.0_pReal
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
     temp33_real = field_real(i,j,k,1:3,1:3) - F(i,j,k,1:3,1:3)
     err_f = err_f + sum(temp33_real*temp33_real)
    
     temp33_real = F_lambda(i,j,k,1:3,1:3) - &
                   (math_mul3333xx33(S_scale,P(i,j,k,1:3,1:3)) + math_I3)
     err_p = err_p + sum(temp33_real*temp33_real)
   enddo; enddo; enddo                              
   err_f = wgt*sqrt(err_f)
   err_p = wgt*sqrt(err_p)

   do k=in(DMDA_LOCAL_INFO_ZS)+1,in(DMDA_LOCAL_INFO_ZS)+in(DMDA_LOCAL_INFO_ZM) 
     do j=in(DMDA_LOCAL_INFO_YS)+1,in(DMDA_LOCAL_INFO_YS)+in(DMDA_LOCAL_INFO_YM) 
       do i=in(DMDA_LOCAL_INFO_XS)+1,in(DMDA_LOCAL_INFO_XS)+in(DMDA_LOCAL_INFO_XM)
         temp33_real = math_mul3333xx33(S_scale,P(i,j,k,1:3,1:3)) + math_I3 - F_lambda(i,j,k,1:3,1:3) &
                + F(i,j,k,1:3,1:3) - field_real(i,j,k,1:3,1:3)
         f_scal(1:9,i,j,k) = reshape(temp33_real,[9])
         f_scal(10:18,i,j,k) = reshape(F(i,j,k,1:3,1:3) - field_real(i,j,k,1:3,1:3),[9])
   enddo; enddo; enddo

   return
 end subroutine AL_formResidual

!--------------------------------------------------------------------------------------------------
!> @brief convergence check
!--------------------------------------------------------------------------------------------------
 subroutine AL_converged(snes_local,it,xnorm,snorm,fnorm,reason,dummy,ierr_psc)
  
   use numerics, only: &
    itmax, &
    itmin, &
    err_f_tol, &
    err_p_tol, &
    err_stress_tolrel, &
    err_stress_tolabs
   
   implicit none

   SNES snes_local
   PetscInt it
   PetscReal xnorm, snorm, fnorm
   SNESConvergedReason reason
   PetscObject dummy
   PetscErrorCode ierr_psc
   logical :: Converged
             
   Converged = (it > itmin .and. &
                 all([ err_f/sqrt(sum((F_aim-math_I3)*(F_aim-math_I3)))/err_f_tol, &
                       err_p/sqrt(sum((F_aim-math_I3)*(F_aim-math_I3)))/err_p_tol, &
                       err_stress/min(maxval(abs(P_av))*err_stress_tolrel,err_stress_tolabs)] < 1.0_pReal))
   
   if (Converged) then
     reason = 1
   elseif (it > itmax) then
     reason = -1
   else  
     reason = 0
   endif 
   
   write(6,'(a,es14.7)') 'error stress BC = ', err_stress/min(maxval(abs(P_av))*err_stress_tolrel,err_stress_tolabs)  
   write(6,'(a,es14.7)') 'error F         = ', err_f/sqrt(sum((F_aim-math_I3)*(F_aim-math_I3)))/err_f_tol
   write(6,'(a,es14.7)') 'error P         = ', err_p/sqrt(sum((F_aim-math_I3)*(F_aim-math_I3)))/err_p_tol
   return

 end subroutine AL_converged

!--------------------------------------------------------------------------------------------------
!> @brief destroy routine
!--------------------------------------------------------------------------------------------------
 subroutine AL_destroy()
 
   use DAMASK_spectral_Utilities, only: &
     Utilities_destroy
   implicit none
   PetscErrorCode ierr_psc
  
   call VecDestroy(solution_vec,ierr_psc)
   call SNESDestroy(snes,ierr_psc)
   call PetscFinalize(ierr_psc)
   call Utilities_destroy()

 end subroutine AL_destroy

 end module DAMASK_spectral_SolverAL
