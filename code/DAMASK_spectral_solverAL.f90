!--------------------------------------------------------------------------------------------------
! $Id: DAMASK_spectral_solverAL.f90 1654 2012-08-03 09:25:48Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief AL scheme solver
!--------------------------------------------------------------------------------------------------
module DAMASK_spectral_solverAL
 use prec, only: & 
   pInt, &
   pReal
 use math, only: &
   math_I3
 use DAMASK_spectral_utilities, only: &
   tSolutionState
 
 implicit none
#include <finclude/petscsys.h>
#include <finclude/petscdmda.h>
#include <finclude/petscsnes.h>

 character (len=*), parameter, public :: &
   DAMASK_spectral_SolverAL_label = 'al'
   
!--------------------------------------------------------------------------------------------------
! derived types ToDo: use here the type definition for a full loadcase including mask
 type tSolutionParams 
   real(pReal), dimension(3,3) :: P_BC, rotation_BC
   real(pReal) :: timeinc
 end type tSolutionParams
 
 type(tSolutionParams), private :: params
 real(pReal), private, dimension(3,3) :: mask_stress = 0.0_pReal

!--------------------------------------------------------------------------------------------------
! PETSc data
 DM,   private :: da
 SNES, private :: snes
 Vec,  private :: solution_vec
 
!--------------------------------------------------------------------------------------------------
! common pointwise data
 real(pReal), private, dimension(:,:,:,:,:), allocatable :: &
   F_lastInc, &                                                                                     !< field of previous compatible deformation gradients
   F_lambda_lastInc, &                                                                              !< field of previous incompatible deformation gradient 
   Fdot, &                                                                                          !< field of assumed rate of compatible deformation gradient
   F_lambdaDot                                                                                      !< field of assumed rate of incopatible deformation gradient
 real(pReal), private                                    ::  temperature                            !< temperature, no spatial quantity at the moment
 
!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
 real(pReal), private, dimension(3,3) :: &
   F_aimDot, &                                                                                      !< assumed rate of average deformation gradient
   F_aim = math_I3, &                                                                               !< current prescribed deformation gradient
   F_aim_lastInc = math_I3, &                                                                       !< previous average deformation gradient
   P_av = 0.0_pReal                                                                                 !< average 1st Piola--Kirchhoff stress
 character(len=1024), private :: incInfo                                                            !< time and increment information
 real(pReal), private, dimension(3,3,3,3) :: &
   C = 0.0_pReal, &                                                                                 !< current average stiffness
   C_lastInc = 0.0_pReal, &                                                                         !< previous average stiffness
   S = 0.0_pReal, &                                                                                 !< current compliance (filled up with zeros)
   C_scale = 0.0_pReal, &                             
   S_scale = 0.0_pReal
 
 real(pReal), private :: &
   err_stress, &                                                                                    !< deviation from stress BC
   err_f, &                                                                                         !< difference between compatible and incompatible deformation gradient
   err_p                                                                                            !< difference of stress resulting from compatible and incompatible F
 logical, private :: ForwardData
 integer(pInt) :: reportIter = 0_pInt
contains
 
!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
subroutine AL_init()
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
   wgt, &
   mesh_NcpElems, &
   mesh_ipCoordinates
 use math, only: &
   math_invSym3333
   
 implicit none
#include <finclude/petscdmda.h90>
#include <finclude/petscsnes.h90>
 integer(pInt) :: i,j,k
 real(pReal), dimension(3,3,  res(1),  res(2),res(3)) ::  P
 real(pReal), dimension(3,3) :: &
   temp33_Real = 0.0_pReal
 real(pReal), dimension(3,3,3,3) :: &
   temp3333_Real = 0.0_pReal, &
   temp3333_Real2 = 0.0_pReal

 PetscErrorCode :: ierr
 PetscObject :: dummy
 PetscScalar, pointer, dimension(:,:,:,:) :: xx_psc, F, F_lambda
 
 call Utilities_init()
 write(6,'(/,a)') ' <<<+-  DAMASK_spectral_solverAL init  -+>>>'
 write(6,'(a)') ' $Id: DAMASK_spectral_SolverAL.f90 1654 2012-08-03 09:25:48Z MPIE\m.diehl $'
#include "compilation_info.f90"
 write(6,'(a)') ''
   
 allocate (F_lastInc  (3,3,  res(1),  res(2),res(3)),  source = 0.0_pReal)
 allocate (Fdot  (3,3,  res(1),  res(2),res(3)),  source = 0.0_pReal)
 !   allocate (Fdot,source = F_lastInc) somethin like that should be possible
 allocate (F_lambda_lastInc(3,3,  res(1),  res(2),res(3)),  source = 0.0_pReal)
 allocate (F_lambdaDot(3,3,  res(1),  res(2),res(3)),  source = 0.0_pReal)
    
!--------------------------------------------------------------------------------------------------
! PETSc Init
 call SNESCreate(PETSC_COMM_WORLD,snes,ierr)
 CHKERRQ(ierr)
 call DMDACreate3d(PETSC_COMM_WORLD,                               &
           DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, &
           DMDA_STENCIL_BOX,res(1),res(2),res(3),PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, &
           18,1,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da,ierr)
 CHKERRQ(ierr)
 call DMCreateGlobalVector(da,solution_vec,ierr)
 CHKERRQ(ierr)
 call DMDASetLocalFunction(da,AL_formResidual,ierr)
 CHKERRQ(ierr)
 call SNESSetDM(snes,da,ierr)
 CHKERRQ(ierr)
 call SNESSetConvergenceTest(snes,AL_converged,dummy,PETSC_NULL_FUNCTION,ierr)
 CHKERRQ(ierr)
 call SNESSetFromOptions(snes,ierr)  
 CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! init fields                 
 call DMDAVecGetArrayF90(da,solution_vec,xx_psc,ierr)                                               ! places pointer xx_psc on PETSc data
 F => xx_psc(0:8,:,:,:)
 F_lambda => xx_psc(9:17,:,:,:)
 if (restartInc == 1_pInt) then                                                                     ! no deformation (no restart)
   F_lastInc         = spread(spread(spread(math_I3,3,res(1)),4,res(2)),5,res(3))                   ! initialize to identity
   F_lambda_lastInc = F_lastInc
   F = reshape(F_lastInc,[9,res(1),res(2),res(3)])
   F_lambda = F
 elseif (restartInc > 1_pInt) then                                                                  ! using old values from file                                                      
   if (debugRestart) write(6,'(a,i6,a)') 'Reading values of increment ',&
                                             restartInc - 1_pInt,' from file' 
   flush(6)
   call IO_read_jobBinaryFile(777,'F',&
                                                trim(getSolverJobName()),size(F))
   read (777,rec=1) F
   close (777)
   call IO_read_jobBinaryFile(777,'F_lastInc',&
                                                trim(getSolverJobName()),size(F_lastInc))
   read (777,rec=1) F_lastInc
   close (777)
   F_aim         = reshape(sum(sum(sum(F,dim=4),dim=3),dim=2) * wgt, [3,3])                         ! average of F
   F_aim_lastInc = sum(sum(sum(F_lastInc,dim=5),dim=4),dim=3) * wgt                                 ! average of F_lastInc 
   call IO_read_jobBinaryFile(777,'F_lambda',&
                                           trim(getSolverJobName()),size(F_lambda))
   read (777,rec=1) F_lambda
   close (777)
   call IO_read_jobBinaryFile(777,'F_lambda_lastInc',&
                                        trim(getSolverJobName()),size(F_lambda_lastInc))
   read (777,rec=1) F_lambda_lastInc
   close (777)
   call IO_read_jobBinaryFile(777,'C_lastInc',trim(getSolverJobName()),size(C_lastInc))
   read (777,rec=1) C_lastInc
   close (777)
   call IO_read_jobBinaryFile(777,'C',trim(getSolverJobName()),size(C))
   read (777,rec=1) C
   close (777)
   call IO_read_jobBinaryFile(777,'F_aimDot',trim(getSolverJobName()),size(f_aimDot))
   read (777,rec=1) f_aimDot
   close (777)
   call IO_read_jobBinaryFile(777,'C_ref',trim(getSolverJobName()),size(temp3333_Real))
   read (777,rec=1) temp3333_Real
   close (777)
 endif
 mesh_ipCoordinates = 0.0_pReal !reshape(mesh_deformedCoordsFFT(geomdim,&
                             !reshape(F,[3,3,res(1),res(2),res(3)])),[3,1,mesh_NcpElems])
 call Utilities_constitutiveResponse(F,F,temperature,0.0_pReal,P,temp3333_Real2,&
                                temp33_Real,.false.,math_I3)
 call DMDAVecRestoreArrayF90(da,solution_vec,xx_psc,ierr)
 CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! reference stiffness
 if (restartInc == 1_pInt) then                                                                     ! use initial stiffness as reference stiffness
   temp3333_Real = temp3333_Real2
   C = temp3333_Real2
 endif 

 call Utilities_updateGamma(temp3333_Real,.True.)
 C_scale = temp3333_Real
 S_scale = math_invSym3333(temp3333_Real)
 
end subroutine AL_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the AL scheme with internal iterations
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function &
  AL_solution(incInfoIn,guess,timeinc,timeinc_old,P_BC,F_BC,temperature_bc,rotation_BC)
 
 use numerics, only: &
   update_gamma
 use math, only: &
   math_mul33x33 ,&
   math_rotate_backward33
 use mesh, only: &
   res,&
   geomdim,&
   mesh_ipCoordinates
 use IO, only: &
   IO_write_JobBinaryFile
 use DAMASK_spectral_Utilities, only: &
   tBoundaryCondition, &
   Utilities_forwardField, &
   Utilities_calculateRate, &
   Utilities_maskedCompliance, &
   Utilities_updateGamma, &
   cutBack
 use FEsolving, only: &
   restartWrite, &
   terminallyIll
 
 implicit none
#include <finclude/petscdmda.h90>
#include <finclude/petscsnes.h90>
!--------------------------------------------------------------------------------------------------
! input data for solution
 real(pReal), intent(in) :: timeinc, timeinc_old, temperature_bc
 logical, intent(in) :: guess
 type(tBoundaryCondition),      intent(in) :: P_BC,F_BC
 character(len=*), intent(in) :: incInfoIn
 real(pReal), dimension(3,3), intent(in) :: rotation_BC

!--------------------------------------------------------------------------------------------------
! PETSc Data
 PetscScalar, dimension(:,:,:,:), pointer :: xx_psc, F, F_lambda
 PetscErrorCode :: ierr   
 SNESConvergedReason ::reason

 incInfo = incInfoIn
 
 call DMDAVecGetArrayF90(da,solution_vec,xx_psc,ierr)
 F => xx_psc(0:8,:,:,:)
 F_lambda => xx_psc(9:17,:,:,:)
 
!--------------------------------------------------------------------------------------------------
! restart information for spectral solver
 if (restartWrite) then
   write(6,'(a)') 'writing converged results for restart'
   call IO_write_jobBinaryFile(777,'F',size(F))                                                     ! writing deformation gradient field to file
   write (777,rec=1) F
   close (777)
   call IO_write_jobBinaryFile(777,'F_lastInc',size(F_lastInc))                                     ! writing F_lastInc field to file
   write (777,rec=1) F_lastInc
   close (777)
   call IO_write_jobBinaryFile(777,'F_lambda',size(F_lambda))                                                     ! writing deformation gradient field to file
   write (777,rec=1) F_lambda
   close (777)
   call IO_write_jobBinaryFile(777,'F_lambda_lastInc',size(F_lambda_lastInc))                                     ! writing F_lastInc field to file
   write (777,rec=1) F_lambda_lastInc
   close (777)
   call IO_write_jobBinaryFile(777,'C',size(C))
   write (777,rec=1) C
   close(777)
   call IO_write_jobBinaryFile(777,'C_lastInc',size(C_lastInc))
   write (777,rec=1) C_lastInc
   close(777)
   call IO_write_jobBinaryFile(777,'F_aimDot',size(F_aimDot))
   write (777,rec=1) F_aimDot
   close(777)
 endif 
 AL_solution%converged =.false.

 if ( cutBack) then 
   F_aim = F_aim_lastInc
   F_lambda= reshape(F_lambda_lastInc,[9,res(1),res(2),res(3)]) 
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
   if (guess) f_aimDot  = f_aimDot + P_BC%maskFloat * (F_aim - F_aim_lastInc)/timeinc_old
   F_aim_lastInc = F_aim

!--------------------------------------------------------------------------------------------------
! update coordinates and rate and forward last inc
   mesh_ipCoordinates = 0.0_pReal !reshape(mesh_deformedCoordsFFT(geomdim,&
                             !reshape(F,[3,3,res(1),res(2),res(3)])),[3,1,mesh_NcpElems])
   Fdot =  Utilities_calculateRate(math_rotate_backward33(f_aimDot,rotation_BC), &
                  timeinc_old,guess,F_lastInc,reshape(F,[3,3,res(1),res(2),res(3)]))
   F_lambdaDot =  Utilities_calculateRate(math_rotate_backward33(f_aimDot,rotation_BC), &
                  timeinc_old,guess,F_lambda_lastInc,reshape(F_lambda,[3,3,res(1),res(2),res(3)]))  
                
   F_lastInc = reshape(F,[3,3,res(1),res(2),res(3)])
   F_lambda_lastInc = reshape(F_lambda,[3,3,res(1),res(2),res(3)])
 endif
 F_aim = F_aim + f_aimDot * timeinc

!--------------------------------------------------------------------------------------------------
! update local deformation gradient and coordinates
!  deltaF_aim = math_rotate_backward33(deltaF_aim,rotation_BC)

 F = reshape(Utilities_forwardField(timeinc,F_aim,F_lastInc,Fdot),[9,res(1),res(2),res(3)])
 F_lambda = reshape(Utilities_forwardField(timeinc,F_aim,F_lambda_lastInc,F_lambdadot),[9,res(1),res(2),res(3)])

 call DMDAVecRestoreArrayF90(da,solution_vec,xx_psc,ierr)
 CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! update stiffness (and gamma operator)
 S = Utilities_maskedCompliance(rotation_BC,P_BC%maskLogical,C)
 if (update_gamma) call Utilities_updateGamma(C,restartWrite)
 
 ForwardData = .True.
 mask_stress = P_BC%maskFloat
 params%P_BC = P_BC%values
 params%rotation_BC = rotation_BC
 params%timeinc = timeinc

 call SNESSolve(snes,PETSC_NULL_OBJECT,solution_vec,ierr)
 CHKERRQ(ierr)
 call SNESGetConvergedReason(snes,reason,ierr)
 CHKERRQ(ierr)
 
 AL_solution%termIll = terminallyIll
 terminallyIll = .false.
 if (reason > 0 ) then 
   AL_solution%converged = .true.
   AL_solution%iterationsNeeded = reportIter
 endif

end function AL_solution


!--------------------------------------------------------------------------------------------------
!> @brief forms the AL residual vector
!--------------------------------------------------------------------------------------------------
subroutine AL_formResidual(in,x_scal,f_scal,dummy,ierr)

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
   Utilities_fourierConvolution, &
   Utilities_FFTbackward, &
   Utilities_constitutiveResponse
 use IO, only: IO_intOut

 implicit none
 integer(pInt) :: i,j,k
 integer(pInt), save :: callNo = 3_pInt
 real(pReal), dimension(3,3)            :: temp33_Real
 logical :: report
 
 DMDALocalInfo :: in(DMDA_LOCAL_INFO_SIZE)
 PetscScalar, target :: x_scal(3,3,2,XG_RANGE,YG_RANGE,ZG_RANGE)  
 PetscScalar, target :: f_scal(3,3,2,X_RANGE,Y_RANGE,Z_RANGE) 
 PetscScalar, pointer :: F(:,:,:,:,:), F_lambda(:,:,:,:,:) 
 PetscScalar, pointer :: residual_F(:,:,:,:,:), residual_F_lambda(:,:,:,:,:)
 PetscInt :: iter, nfuncs
 PetscObject :: dummy
 PetscErrorCode :: ierr

 F => x_scal(:,:,1,:,:,:)
 F_lambda => x_scal(:,:,2,:,:,:)
 residual_F => f_scal(:,:,1,:,:,:)
 residual_F_lambda => f_scal(:,:,2,:,:,:)
 
 call SNESGetNumberFunctionEvals(snes,nfuncs,ierr)
 CHKERRQ(ierr)
 call SNESGetIterationNumber(snes,iter,ierr)
 CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! report begin of new iteration
 if (iter == 0 .and. callNo>2) then
   callNo = 0_pInt
   reportIter = 0_pInt
 endif
 if (callNo == 0 .or. mod(callNo,2) == 1_pInt) then
   write(6,'(/,a,3(a,'//IO_intOut(itmax)//'))') trim(incInfo), &
                    ' @ Iter. ', itmin, '<',reportIter, '≤', itmax
   write(6,'(a,/,3(3(f12.7,1x)/))',advance='no') 'deformation gradient aim =',&
                                                             math_transpose33(F_aim)
   reportIter = reportIter + 1_pInt
 endif
 callNo = callNo +1_pInt

!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
 call Utilities_constitutiveResponse(F_lastInc,F,temperature,params%timeinc, &
                                     residual_F,C,P_av,ForwardData,params%rotation_BC)
 ForwardData = .False.
  
!--------------------------------------------------------------------------------------------------
! stress BC handling
 F_aim = F_aim - math_mul3333xx33(S, ((P_av - params%P_BC))) ! S = 0.0 for no bc
 err_stress = maxval(abs(mask_stress * (P_av - params%P_BC)))     ! mask = 0.0 for no bc
  
!--------------------------------------------------------------------------------------------------
! 
 field_real = 0.0_pReal
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   temp33_Real = math_mul3333xx33(S_scale,residual_F(1:3,1:3,i,j,k)) + math_I3
   residual_F(1:3,1:3,i,j,k) = temp33_Real
   field_real(i,j,k,1:3,1:3) = -math_mul3333xx33(C_scale,F_lambda(1:3,1:3,i,j,k)-F(1:3,1:3,i,j,k))
 enddo; enddo; enddo
 
!--------------------------------------------------------------------------------------------------
! doing convolution in Fourier space 
 call Utilities_FFTforward()
 call Utilities_fourierConvolution(math_rotate_backward33(F_aim,params%rotation_BC)) 
 call Utilities_FFTbackward()
 
!--------------------------------------------------------------------------------------------------
! constructing residual                         
 residual_F_lambda = F - reshape(field_real(1:res(1),1:res(2),1:res(3),1:3,1:3),&
                                 [3,3,res(1),res(2),res(3)],order=[3,4,5,1,2])
 residual_F = residual_F - F_lambda + residual_F_lambda
 
!--------------------------------------------------------------------------------------------------
! calculating errors  
 err_f = wgt*sqrt(sum(residual_F_lambda**2.0_pReal))
 err_p = wgt*sqrt(sum((residual_F - residual_F_lambda)**2.0_pReal))
   
end subroutine AL_formResidual


!--------------------------------------------------------------------------------------------------
!> @brief convergence check
!--------------------------------------------------------------------------------------------------
subroutine AL_converged(snes_local,it,xnorm,snorm,fnorm,reason,dummy,ierr)

 use numerics, only: &
  itmax, &
  itmin, &
  err_f_tol, &
  err_p_tol, &
  err_stress_tolrel, &
  err_stress_tolabs
 
 implicit none

 SNES :: snes_local
 PetscInt :: it
 PetscReal :: xnorm, snorm, fnorm
 SNESConvergedReason :: reason
 PetscObject :: dummy
 PetscErrorCode ::ierr
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
 
 write(6,'(a,f12.7,1x,1a,1x,es9.3)') 'error stress BC = ', &
   err_stress/min(maxval(abs(P_av))*err_stress_tolrel,err_stress_tolabs),&
   '@',min(maxval(abs(P_av))*err_stress_tolrel,err_stress_tolabs)
 write(6,'(a,f12.7,1x,1a,1x,es9.3)') 'error F         = ',&
   err_f/sqrt(sum((F_aim-math_I3)*(F_aim-math_I3)))/err_f_tol,&
   '@',err_f_tol
 write(6,'(a,f12.7,1x,1a,1x,es9.3)') 'error P         = ', &
   err_p/sqrt(sum((F_aim-math_I3)*(F_aim-math_I3)))/err_p_tol,&
   '@',err_p_tol
 write(6,'(/,a)') '=========================================================================='
end subroutine AL_converged

!--------------------------------------------------------------------------------------------------
!> @brief destroy routine
!--------------------------------------------------------------------------------------------------
subroutine AL_destroy()

 use DAMASK_spectral_Utilities, only: &
   Utilities_destroy
 implicit none
 PetscErrorCode :: ierr

 call VecDestroy(solution_vec,ierr)
 CHKERRQ(ierr)
 call SNESDestroy(snes,ierr)
 CHKERRQ(ierr)
 call DMDestroy(da,ierr)
 CHKERRQ(ierr)
 call PetscFinalize(ierr)
 CHKERRQ(ierr)
 call Utilities_destroy()

end subroutine AL_destroy

end module DAMASK_spectral_SolverAL
