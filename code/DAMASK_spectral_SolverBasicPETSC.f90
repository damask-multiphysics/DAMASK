!--------------------------------------------------------------------------------------------------
! $Id: DAMASK_spectral_SolverAL.f90 1654 2012-08-03 09:25:48Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Basic scheme PETSc solver
!--------------------------------------------------------------------------------------------------
module DAMASK_spectral_SolverBasicPETSC
 use prec, only: & 
   pInt, &
   pReal
 
 use math, only: &
   math_I3
 
 use DAMASK_spectral_Utilities, only: &
   tSolutionState

 implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscdmda.h>
#include <finclude/petscis.h>
#include <finclude/petscmat.h>
#include <finclude/petscksp.h>
#include <finclude/petscpc.h>
#include <finclude/petscsnes.h>
#include <finclude/petscvec.h90>
#include <finclude/petscdmda.h90>
#include <finclude/petscsnes.h90>

 character (len=*), parameter, public :: &
   DAMASK_spectral_SolverBasicPETSC_label = 'basicpetsc'
   
!--------------------------------------------------------------------------------------------------
! derived types
 type tSolutionParams 
   real(pReal), dimension(3,3) :: P_BC, rotation_BC
   real(pReal) :: timeinc
 end type tSolutionParams
 
 type(tSolutionParams), private :: params

!--------------------------------------------------------------------------------------------------
! PETSc data
  DM, private :: da
  SNES, private :: snes
  Vec, private :: solution_vec
 
!--------------------------------------------------------------------------------------------------
! common pointwise data
  real(pReal), private, dimension(:,:,:,:,:), allocatable ::  F_lastInc, Fdot
  real(pReal), private, dimension(:,:,:,:),   allocatable ::  coordinates
  real(pReal), private, dimension(:,:,:),     allocatable ::  temperature
 
!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
  real(pReal), private, dimension(3,3) :: &
    F_aim = math_I3, &
    F_aim_lastInc = math_I3, &
    P_av
  character(len=1024), private :: incInfo   
  real(pReal), private, dimension(3,3,3,3) :: &
    C = 0.0_pReal, C_lastInc= 0.0_pReal, &
    S = 0.0_pReal
 
  real(pReal), private :: err_stress, err_div
  logical, private :: ForwardData
  real(pReal), private, dimension(3,3) :: mask_stress = 0.0_pReal
 
  contains
 
!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
subroutine BasicPETSC_init()
      
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 
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
   debugRestart
      
 use numerics, only: &
   petsc_options  
         
 use mesh, only: &
   res, &
   geomdim, &
   mesh_NcpElems
      
 use math, only: &
   math_invSym3333
      
 implicit none
 integer(pInt) :: i,j,k
 real(pReal), dimension(:,:,:,:,:), allocatable ::  P
    
 PetscErrorCode :: ierr_psc
 PetscObject :: dummy
 PetscMPIInt :: rank
 PetscScalar, pointer :: xx_psc(:,:,:,:), F(:,:,:,:)
    
    call Utilities_init()
    
    write(6,'(/,a)') ' <<<+-  DAMASK_spectral_solverBasicPETSC init  -+>>>'
    write(6,'(a)') ' $Id: DAMASK_spectral_SolverBasicPETSC.f90 1654 2012-08-03 09:25:48Z MPIE\m.diehl $'
#include "compilation_info.f90"
    write(6,'(a)') ''
   
    allocate (F_lastInc  (3,3,  res(1),  res(2),res(3)),  source = 0.0_pReal)
    allocate (Fdot  (3,3,  res(1),  res(2),res(3)),  source = 0.0_pReal)
    allocate (P          (3,3,  res(1),  res(2),res(3)),  source = 0.0_pReal)
    allocate (coordinates(  res(1),  res(2),res(3),3),    source = 0.0_pReal)
    allocate (temperature(  res(1),  res(2),res(3)),      source = 0.0_pReal)
    
 !--------------------------------------------------------------------------------------------------
 ! PETSc Init
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr_psc)
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr_psc)
    call SNESCreate(PETSC_COMM_WORLD,snes,ierr_psc)
    call DMDACreate3d(PETSC_COMM_WORLD,                               &
              DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, &
              DMDA_STENCIL_BOX,res(1),res(2),res(3),PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, &
              9,1,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da,ierr_psc)

    call DMCreateGlobalVector(da,solution_vec,ierr_psc)
    call DMDASetLocalFunction(da,BasicPETSC_formResidual,ierr_psc)
    call SNESSetDM(snes,da,ierr_psc)
    call SNESSetConvergenceTest(snes,BasicPETSC_converged,dummy,PETSC_NULL_FUNCTION,ierr_psc)
    call PetscOptionsInsertString(petsc_options,ierr_psc)
    call SNESSetFromOptions(snes,ierr_psc)  

 !--------------------------------------------------------------------------------------------------
 ! init fields                 
    call DMDAVecGetArrayF90(da,solution_vec,xx_psc,ierr_psc)
    F => xx_psc(0:8,:,:,:)
    if (restartInc == 1_pInt) then                                                                     ! no deformation (no restart)
      F_lastInc         = spread(spread(spread(math_I3,3,res(1)),4,res(2)),5,res(3))                   ! initialize to identity
      F = reshape(F_lastInc,[9,res(1),res(2),res(3)])
      do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
        coordinates(i,j,k,1:3) = geomdim/real(res,pReal)*real([i,j,k],pReal) &
                               - geomdim/real(2_pInt*res,pReal)
      enddo; enddo; enddo
    elseif (restartInc > 1_pInt) then                                                                  ! using old values from file                                                      
      if (debugRestart) write(6,'(a,i6,a)') 'Reading values of increment ',&
                                                restartInc - 1_pInt,' from file' 
      call IO_read_jobBinaryFile(777,'convergedSpectralDefgrad',&
                                                   trim(getSolverJobName()),size(F_lastInc))
      read (777,rec=1) F
      close (777)
      call IO_read_jobBinaryFile(777,'convergedSpectralDefgrad_lastInc',&
                                                   trim(getSolverJobName()),size(F_lastInc))
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
   
    call Utilities_constitutiveResponse(coordinates,F,F,temperature,0.0_pReal,P,C,P_av,.false.,math_I3)
    call DMDAVecRestoreArrayF90(da,solution_vec,xx_psc,ierr_psc)

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
 
  end subroutine BasicPETSC_init
  
!--------------------------------------------------------------------------------------------------
!> @brief solution for the Basic PETSC scheme with internal iterations
!--------------------------------------------------------------------------------------------------
  type(tSolutionState) function &
    basicPETSc_solution(incInfoIn,guessmode,timeinc,timeinc_old,P_BC,F_BC,temperature_bc,rotation_BC)
   use numerics, only: &
     update_gamma
   use math, only: &
     math_mul33x33 ,&
     math_rotate_backward33
   use mesh, only: &
     res,&
     geomdim,&
     deformed_fft
   use IO, only: &
     IO_write_JobBinaryFile
   use DAMASK_spectral_Utilities, only: &
     tBoundaryCondition, &
     Utilities_calculateRate, &
     Utilities_forwardField, &
     Utilities_maskedCompliance, &
     Utilities_updateGamma, &
     cutBack
   use FEsolving, only: &
     restartWrite, &
     terminallyIll
   implicit none
!--------------------------------------------------------------------------------------------------
! input data for solution
   real(pReal), intent(in) :: timeinc, timeinc_old, temperature_bc, guessmode
   type(tBoundaryCondition),      intent(in) :: P_BC,F_BC
   real(pReal), dimension(3,3), intent(in) :: rotation_BC
 character(len=*), intent(in) :: incInfoIn
   real(pReal), dimension(3,3),save       :: F_aimDot=0.0_pReal
   real(pReal), dimension(3,3)            :: &
                                             F_aim_lab
!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
   real(pReal), dimension(3,3)            :: temp33_Real 
 
!--------------------------------------------------------------------------------------------------
! 
   PetscScalar, pointer :: xx_psc(:,:,:,:), F(:,:,:,:)
   PetscErrorCode :: ierr_psc   
   SNESConvergedReason :: reason

!--------------------------------------------------------------------------------------------------
! restart information for spectral solver
   incInfo = incInfoIn
   if (restartWrite) then
     write(6,'(a)') 'writing converged results for restart'
     call IO_write_jobBinaryFile(777,'convergedSpectralDefgrad',size(F_lastInc))
     write (777,rec=1) F_LastInc
     close (777)
     call IO_write_jobBinaryFile(777,'C',size(C))
     write (777,rec=1) C
     close(777)
   endif 
  call DMDAVecGetArrayF90(da,solution_vec,xx_psc,ierr_psc)
  F => xx_psc(0:8,:,:,:)
 
if ( cutBack) then 
  F_aim = F_aim_lastInc

  F = reshape(F_lastInc,[9,res(1),res(2),res(3)]) 
  C = C_lastInc
else
 
  C_lastInc = C

!--------------------------------------------------------------------------------------------------
! calculate rate for aim
   if (F_BC%myType=='l') then                                                                       ! calculate f_aimDot from given L and current F
     f_aimDot = F_BC%maskFloat * math_mul33x33(F_BC%values, F_aim)
   elseif(F_BC%myType=='fdot')   then                                                               ! f_aimDot is prescribed
     f_aimDot = F_BC%maskFloat * F_BC%values
   endif
   f_aimDot  = f_aimDot &                                                                         
             + guessmode * P_BC%maskFloat * (F_aim - F_aim_lastInc)/timeinc_old
   F_aim_lastInc = F_aim

!--------------------------------------------------------------------------------------------------
! update coordinates and rate and forward last inc
   call deformed_fft(res,geomdim,math_rotate_backward33(F_aim_lastInc,rotation_BC), &
                                                                   1.0_pReal,F_lastInc,coordinates)
   Fdot =  Utilities_calculateRate(math_rotate_backward33(f_aimDot,rotation_BC), &
                                           timeinc,timeinc_old,guessmode,F_lastInc,reshape(F,[3,3,res(1),res(2),res(3)]))
   F_lastInc = reshape(F,[3,3,res(1),res(2),res(3)])
 endif
 F_aim = F_aim + f_aimDot * timeinc


 F = reshape(Utilities_forwardField(timeinc,F_aim,F_lastInc,Fdot),[9,res(1),res(2),res(3)])
 call DMDAVecRestoreArrayF90(da,solution_vec,xx_psc,ierr_psc)
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

   call SNESSolve(snes,PETSC_NULL_OBJECT,solution_vec,ierr_psc)
   call SNESGetConvergedReason(snes,reason,ierr_psc)
   basicPETSc_solution%termIll = terminallyIll
   terminallyIll = .false.
   BasicPETSC_solution%converged =.false.
   if (reason > 0 ) BasicPETSC_solution%converged = .true.

 end function BasicPETSC_solution

!--------------------------------------------------------------------------------------------------
!> @brief forms the AL residual vector
!--------------------------------------------------------------------------------------------------
 subroutine BasicPETSC_formResidual(in,x_scal,f_scal,dummy,ierr_psc)
  
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
     Utilities_FFTforward, &
     Utilities_FFTbackward, &
     Utilities_fourierConvolution, &
     Utilities_constitutiveResponse, &
     Utilities_divergenceRMS
   use IO, only : IO_intOut 
   implicit none

   real(pReal), dimension(3,3) :: F_aim_lab_lastIter, F_aim_lab
   
   DMDALocalInfo :: in(DMDA_LOCAL_INFO_SIZE)
   PetscScalar, target :: x_scal(3,3,XG_RANGE,YG_RANGE,ZG_RANGE)  
   PetscScalar, target :: f_scal(3,3,X_RANGE,Y_RANGE,Z_RANGE) 
   PetscScalar, pointer :: F(:,:,:,:,:)
   PetscScalar, pointer :: residual_F(:,:,:,:,:)
   PetscInt :: iter, nfuncs
   PetscObject :: dummy
   PetscErrorCode :: ierr_psc
   F => x_scal(:,:,:,:,:)
   residual_F => f_scal(:,:,:,:,:)
   
   call SNESGetNumberFunctionEvals(snes,nfuncs,ierr_psc)
   call SNESGetIterationNumber(snes,iter,ierr_psc)
  
 !--------------------------------------------------------------------------------------------------
 ! report begin of new iteration
   write(6,'(/,a,3(a,'//IO_intOut(itmax)//'))') trim(incInfo), &
                    ' @ Iter. ', itmin, '<',iter, '≤', itmax
   write(6,'(a,/,3(3(f12.7,1x)/))',advance='no') 'deformation gradient aim =',&
                                                             math_transpose33(F_aim)
   F_aim_lab_lastIter = math_rotate_backward33(F_aim,params%rotation_BC)
   
 !--------------------------------------------------------------------------------------------------
 ! evaluate constitutive response
   call Utilities_constitutiveResponse(coordinates,F_lastInc,F,temperature,params%timeinc, &
                                       residual_F,C,P_av,ForwardData,params%rotation_BC)
   ForwardData = .false.
   
!--------------------------------------------------------------------------------------------------
! stress BC handling
   F_aim = F_aim - math_mul3333xx33(S, ((P_av - params%P_BC))) ! S = 0.0 for no bc
   err_stress = maxval(abs(mask_stress * (P_av - params%P_BC)))     ! mask = 0.0 for no bc
   F_aim_lab = math_rotate_backward33(F_aim,params%rotation_BC) 
   
!--------------------------------------------------------------------------------------------------
! updated deformation gradient using fix point algorithm of basic scheme
   field_real = 0.0_pReal
   field_real(1:res(1),1:res(2),1:res(3),1:3,1:3) = reshape(residual_F,[res(1),res(2),res(3),3,3],&
                                                               order=[4,5,1,2,3]) ! field real has a different order
   call Utilities_FFTforward()
   err_div = Utilities_divergenceRMS()
   call Utilities_fourierConvolution(F_aim_lab_lastIter - F_aim_lab) 
   call Utilities_FFTbackward()
   
 !--------------------------------------------------------------------------------------------------
 ! constructing residual                         
   residual_F = reshape(field_real(1:res(1),1:res(2),1:res(3),1:3,1:3),shape(F),order=[3,4,5,1,2])  
   write(6,'(/,a)') '=========================================================================='
 end subroutine BasicPETSC_formResidual

!--------------------------------------------------------------------------------------------------
!> @brief convergence check
!--------------------------------------------------------------------------------------------------
 subroutine BasicPETSC_converged(snes_local,it,xnorm,snorm,fnorm,reason,dummy,ierr_psc)
  
   use numerics, only: &
    itmax, &
    itmin, &
    err_div_tol, &
    err_stress_tolrel, &
    err_stress_tolabs
  
   use math, only: &
     math_mul33x33, &
     math_eigenvalues33, &
     math_transpose33
   
   implicit none

   SNES snes_local
   PetscInt it
   PetscReal xnorm, snorm, fnorm
   SNESConvergedReason reason
   PetscObject dummy
   PetscErrorCode ierr_psc
   logical :: Converged
   real(pReal) :: pAvgDivL2
             
   pAvgDivL2 = sqrt(maxval(math_eigenvalues33(math_mul33x33(P_av,math_transpose33(P_av)))))
   Converged = (it > itmin .and. &
                 all([ err_div/pAvgDivL2/err_div_tol, &
                       err_stress/min(maxval(abs(P_av))*err_stress_tolrel,err_stress_tolabs)] < 1.0_pReal))
   
   if (Converged) then
     reason = 1
   elseif (it > itmax) then
     reason = -1
   else  
     reason = 0
   endif 
  
   write(6,'(a,f6.2,a,es11.4,a)') 'error divergence = ', err_div/pAvgDivL2/err_div_tol,&
                                                         ' (',err_div/pAvgDivL2,' N/m³)'
   write(6,'(a,f6.2,a,es11.4,a)') 'error stress =     ', err_stress/min(maxval(abs(P_av))*err_stress_tolrel,err_stress_tolabs), &
                                                         ' (',err_stress,' Pa)'  
   return

 end subroutine BasicPETSC_converged

!--------------------------------------------------------------------------------------------------
!> @brief destroy routine
!--------------------------------------------------------------------------------------------------
 subroutine BasicPETSC_destroy()
 
   use DAMASK_spectral_Utilities, only: &
     Utilities_destroy
   implicit none
   PetscErrorCode ierr_psc
  
   call VecDestroy(solution_vec,ierr_psc)
   call SNESDestroy(snes,ierr_psc)
   call DMDestroy(da,ierr_psc)
   call PetscFinalize(ierr_psc)
   call Utilities_destroy()

 end subroutine BasicPETSC_destroy

 end module DAMASK_spectral_SolverBasicPETSC
