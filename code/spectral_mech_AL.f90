!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief AL scheme solver
!--------------------------------------------------------------------------------------------------
module spectral_mech_AL
 use prec, only: & 
   pInt, &
   pReal
 use math, only: &
   math_I3
 use spectral_utilities, only: &
   tSolutionState, &
   tSolutionParams

 implicit none
 private
#include <petsc/finclude/petsc.h90>

 character (len=*), parameter, public :: &
   DAMASK_spectral_solverAL_label = 'al'
   
!--------------------------------------------------------------------------------------------------
! derived types 
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

!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
 real(pReal), private, dimension(3,3) :: &
   F_aimDot, &                                                                                      !< assumed rate of average deformation gradient
   F_aim = math_I3, &                                                                               !< current prescribed deformation gradient
   F_aim_lastInc = math_I3, &                                                                       !< previous average deformation gradient
   F_av = 0.0_pReal, &                                                                              !< average incompatible def grad field
   P_av = 0.0_pReal, &                                                                              !< average 1st Piola--Kirchhoff stress
   P_avLastEval = 0.0_pReal                                                                         !< average 1st Piola--Kirchhoff stress last call of CPFEM_general
 character(len=1024), private :: incInfo                                                            !< time and increment information
 real(pReal), private, dimension(3,3,3,3) :: &
   C_volAvg = 0.0_pReal, &                                                                          !< current volume average stiffness 
   C_volAvgLastInc = 0.0_pReal, &                                                                   !< previous volume average stiffness
   C_minMaxAvg = 0.0_pReal, &                                                                       !< current (min+max)/2 stiffness
   S = 0.0_pReal, &                                                                                 !< current compliance (filled up with zeros)
   C_scale = 0.0_pReal, &                             
   S_scale = 0.0_pReal
 
 real(pReal), private :: &
   err_BC, &                                                                                        !< deviation from stress BC
   err_curl, &                                                                                      !< RMS of curl of F
   err_div                                                                                          !< RMS of div of P
 logical, private :: ForwardData
 integer(pInt), private :: &
   totalIter = 0_pInt                                                                               !< total iteration in current increment
 
 public :: &
   AL_init, &
   AL_solution, &
   AL_forward, &
   AL_destroy
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
!> @todo use sourced allocation, e.g. allocate(Fdot,source = F_lastInc)
!--------------------------------------------------------------------------------------------------
subroutine AL_init
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 use IO, only: &
   IO_intOut, &
   IO_read_realFile, &
   IO_timeStamp
 use debug, only: &
  debug_level, &
  debug_spectral, &
  debug_spectralRestart
 use FEsolving, only: &
   restartInc
 use numerics, only: &
   worldrank, &
   worldsize
 use DAMASK_interface, only: &
   getSolverJobName
 use spectral_utilities, only: &
   Utilities_constitutiveResponse, &
   Utilities_updateGamma, &
   Utilities_updateIPcoords
 use mesh, only: &
   grid, &
   grid3
 use math, only: &
   math_invSym3333
   
 implicit none
 real(pReal), dimension(3,3,grid(1),grid(2),grid3) :: P
 real(pReal), dimension(3,3) :: &
   temp33_Real = 0.0_pReal

 PetscErrorCode :: ierr
 PetscObject :: dummy
 PetscScalar, pointer, dimension(:,:,:,:) :: xx_psc, F, F_lambda
 integer(pInt), dimension(:), allocatable :: localK  
 integer(pInt) :: proc
 character(len=1024) :: rankStr
 
 if (worldrank == 0_pInt) then
   write(6,'(/,a)') ' <<<+-  DAMASK_spectral_solverAL init  -+>>>'
   write(6,'(a15,a)')   ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif

!--------------------------------------------------------------------------------------------------
! allocate global fields
 allocate (F_lastInc       (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
 allocate (Fdot            (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
 allocate (F_lambda_lastInc(3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
 allocate (F_lambdaDot     (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
    
!--------------------------------------------------------------------------------------------------
! PETSc Init
 call SNESCreate(PETSC_COMM_WORLD,snes,ierr); CHKERRQ(ierr)
 call SNESSetOptionsPrefix(snes,'mech_',ierr);CHKERRQ(ierr) 
 allocate(localK(worldsize), source = 0); localK(worldrank+1) = grid3
 do proc = 1, worldsize
   call MPI_Bcast(localK(proc),1,MPI_INTEGER,proc-1,PETSC_COMM_WORLD,ierr)
 enddo  
 call DMDACreate3d(PETSC_COMM_WORLD, &
        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &                                     ! cut off stencil at boundary
        DMDA_STENCIL_BOX, &                                                                         ! Moore (26) neighborhood around central point
        grid(1),grid(2),grid(3), &                                                                  ! global grid
        1 , 1, worldsize, &
        18, 0, &                                                                                    ! #dof (F tensor), ghost boundary width (domain overlap)
        grid(1),grid(2),localK, &                                                                   ! local grid
        da,ierr)                                                                                    ! handle, error
 CHKERRQ(ierr)
 call SNESSetDM(snes,da,ierr); CHKERRQ(ierr)
 call DMCreateGlobalVector(da,solution_vec,ierr); CHKERRQ(ierr)
 call DMDASNESSetFunctionLocal(da,INSERT_VALUES,AL_formResidual,dummy,ierr)
 CHKERRQ(ierr)
 call SNESSetConvergenceTest(snes,AL_converged,dummy,PETSC_NULL_FUNCTION,ierr)
 CHKERRQ(ierr)
 call SNESSetFromOptions(snes,ierr); CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! init fields                 
 call DMDAVecGetArrayF90(da,solution_vec,xx_psc,ierr); CHKERRQ(ierr)                                ! places pointer xx_psc on PETSc data
 F => xx_psc(0:8,:,:,:)
 F_lambda => xx_psc(9:17,:,:,:)
 restart: if (restartInc > 1_pInt) then
   if (iand(debug_level(debug_spectral),debug_spectralRestart)/= 0 .and. worldrank == 0_pInt) &
     write(6,'(/,a,'//IO_intOut(restartInc-1_pInt)//',a)') &
     'reading values of increment ', restartInc - 1_pInt, ' from file'
   flush(6)
   write(rankStr,'(a1,i0)')'_',worldrank
   call IO_read_realFile(777,'F'//trim(rankStr), trim(getSolverJobName()),size(F))
   read (777,rec=1) F
   close (777)
   call IO_read_realFile(777,'F_lastInc'//trim(rankStr), trim(getSolverJobName()),size(F_lastInc))
   read (777,rec=1) F_lastInc
   close (777)
   call IO_read_realFile(777,'F_lambda'//trim(rankStr),trim(getSolverJobName()),size(F_lambda))
   read (777,rec=1) F_lambda
   close (777)
   call IO_read_realFile(777,'F_lambda_lastInc'//trim(rankStr),&
                                        trim(getSolverJobName()),size(F_lambda_lastInc))
   read (777,rec=1) F_lambda_lastInc
   close (777)
   call IO_read_realFile(777,'F_aim', trim(getSolverJobName()),size(F_aim))
   read (777,rec=1) F_aim
   close (777)
   call IO_read_realFile(777,'F_aim_lastInc', trim(getSolverJobName()),size(F_aim_lastInc))
   read (777,rec=1) F_aim_lastInc
   close (777)
   call IO_read_realFile(777,'F_aimDot',trim(getSolverJobName()),size(f_aimDot))
   read (777,rec=1) f_aimDot
   close (777)
 elseif (restartInc == 1_pInt) then restart
   F_lastInc = spread(spread(spread(math_I3,3,grid(1)),4,grid(2)),5,grid3)                          ! initialize to identity
   F = reshape(F_lastInc,[9,grid(1),grid(2),grid3])
   F_lambda = F
   F_lambda_lastInc = F_lastInc
 endif restart

 
 call Utilities_updateIPcoords(reshape(F,shape(F_lastInc)))
 call Utilities_constitutiveResponse(F_lastInc, reshape(F,shape(F_lastInc)), &
                   0.0_pReal,P,C_volAvg,C_minMaxAvg,temp33_Real,.false.,math_I3)
 nullify(F)
 nullify(F_lambda)
 call DMDAVecRestoreArrayF90(da,solution_vec,xx_psc,ierr); CHKERRQ(ierr)                            ! write data back to PETSc
                             
 readRestart: if (restartInc > 1_pInt) then
   if (iand(debug_level(debug_spectral),debug_spectralRestart)/= 0 .and. worldrank == 0_pInt) &
     write(6,'(/,a,'//IO_intOut(restartInc-1_pInt)//',a)') &
     'reading more values of increment', restartInc - 1_pInt, 'from file'
   flush(6)
   call IO_read_realFile(777,'C_volAvg',trim(getSolverJobName()),size(C_volAvg))
   read (777,rec=1) C_volAvg
   close (777)
   call IO_read_realFile(777,'C_volAvgLastInc',trim(getSolverJobName()),size(C_volAvgLastInc))
   read (777,rec=1) C_volAvgLastInc
   close (777)
   call IO_read_realFile(777,'C_ref',trim(getSolverJobName()),size(C_minMaxAvg))
   read (777,rec=1) C_minMaxAvg
   close (777)
 endif readRestart

 call Utilities_updateGamma(C_minMaxAvg,.True.)
 C_scale = C_minMaxAvg
 S_scale = math_invSym3333(C_minMaxAvg)
 
end subroutine AL_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the AL scheme with internal iterations
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function &
  AL_solution(incInfoIn,guess,timeinc,timeinc_old,loadCaseTime,P_BC,F_BC,rotation_BC)
 use IO, only: &
   IO_error
 use numerics, only: &
   update_gamma
 use math, only: &
   math_invSym3333
 use spectral_utilities, only: &
   tBoundaryCondition, &
   Utilities_maskedCompliance, &
   Utilities_updateGamma
 use FEsolving, only: &
   restartWrite, &
   terminallyIll
 
 implicit none

!--------------------------------------------------------------------------------------------------
! input data for solution
 real(pReal), intent(in) :: &
   timeinc, &                                                                                       !< increment in time for current solution
   timeinc_old, &                                                                                   !< increment in time of last increment
   loadCaseTime                                                                                     !< remaining time of current load case
 logical, intent(in) :: &
   guess
 type(tBoundaryCondition),      intent(in) :: &
   P_BC, &
   F_BC
 character(len=*), intent(in) :: &
   incInfoIn
 real(pReal), dimension(3,3), intent(in) :: rotation_BC
 
!--------------------------------------------------------------------------------------------------
! PETSc Data
 PetscErrorCode :: ierr   
 SNESConvergedReason :: reason

 incInfo = incInfoIn

!--------------------------------------------------------------------------------------------------
! update stiffness (and gamma operator)
 S = Utilities_maskedCompliance(rotation_BC,P_BC%maskLogical,C_volAvg)
 if (update_gamma) then
   call Utilities_updateGamma(C_minMaxAvg,restartWrite)
   C_scale = C_minMaxAvg
   S_scale = math_invSym3333(C_minMaxAvg)
 endif  

!--------------------------------------------------------------------------------------------------
! set module wide availabe data 
 mask_stress = P_BC%maskFloat
 params%P_BC = P_BC%values
 params%rotation_BC = rotation_BC
 params%timeinc = timeinc
 params%timeincOld = timeinc_old

!--------------------------------------------------------------------------------------------------
! solve BVP 
 call SNESSolve(snes,PETSC_NULL_OBJECT,solution_vec,ierr)
 CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! check convergence
 call SNESGetConvergedReason(snes,reason,ierr)
 CHKERRQ(ierr)
 AL_solution%termIll = terminallyIll
 terminallyIll = .false.
 if (reason == -4) call IO_error(893_pInt)
 if (reason < 1) AL_solution%converged = .false.
 AL_solution%iterationsNeeded = totalIter

end function AL_solution


!--------------------------------------------------------------------------------------------------
!> @brief forms the AL residual vector
!--------------------------------------------------------------------------------------------------
subroutine AL_formResidual(in,x_scal,f_scal,dummy,ierr)
 use numerics, only: &
   itmax, &
   itmin, &
   polarAlpha, &
   polarBeta, &
   worldrank
 use mesh, only: &
   grid3, &
   grid
 use IO, only: &
   IO_intOut
 use math, only: &
   math_rotate_backward33, &
   math_transpose33, &
   math_mul3333xx33, &
   math_invSym3333, &
   math_mul33x33
 use spectral_utilities, only: &
   wgt, &
   tensorField_real, &
   utilities_FFTtensorForward, &
   utilities_fourierGammaConvolution, &
   utilities_FFTtensorBackward, &
   Utilities_constitutiveResponse, &
   Utilities_divergenceRMS, &
   Utilities_curlRMS
 use debug, only: &
   debug_level, &
   debug_spectral, &
   debug_spectralRotation
 use homogenization, only: &
   materialpoint_dPdF
 use FEsolving, only: &
   terminallyIll

 implicit none
!--------------------------------------------------------------------------------------------------
! strange syntax in the next line because otherwise macros expand beyond 132 character limit 
 DMDALocalInfo,        dimension(&
   DMDA_LOCAL_INFO_SIZE) :: &
   in
 PetscScalar, target, dimension(3,3,2, &
   XG_RANGE,YG_RANGE,ZG_RANGE) :: &
   x_scal
 PetscScalar, target, dimension(3,3,2, &
   X_RANGE,Y_RANGE,Z_RANGE) :: &
   f_scal
 PetscScalar, pointer, dimension(:,:,:,:,:) :: &
   F, &
   F_lambda, &
   residual_F, &
   residual_F_lambda
 PetscInt :: &
   PETScIter, &
   nfuncs
 PetscObject :: dummy
 PetscErrorCode :: ierr
 integer(pInt) :: &
   i, j, k, e

 F                => x_scal(1:3,1:3,1,&
  XG_RANGE,YG_RANGE,ZG_RANGE)
 F_lambda         => x_scal(1:3,1:3,2,&
  XG_RANGE,YG_RANGE,ZG_RANGE)
 residual_F       => f_scal(1:3,1:3,1,&
  X_RANGE,Y_RANGE,Z_RANGE)
 residual_F_lambda => f_scal(1:3,1:3,2,&
  X_RANGE,Y_RANGE,Z_RANGE)
 
 call SNESGetNumberFunctionEvals(snes,nfuncs,ierr); CHKERRQ(ierr)
 call SNESGetIterationNumber(snes,PETScIter,ierr); CHKERRQ(ierr)

 F_av = sum(sum(sum(F,dim=5),dim=4),dim=3) * wgt
 call MPI_Allreduce(MPI_IN_PLACE,F_av,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 
 if(nfuncs== 0 .and. PETScIter == 0) totalIter = -1_pInt                                            ! new increment
 newIteration: if(totalIter <= PETScIter) then
!--------------------------------------------------------------------------------------------------
! report begin of new iteration
   totalIter = totalIter + 1_pInt
   if (worldrank == 0_pInt) then
     write(6,'(1x,a,3(a,'//IO_intOut(itmax)//'))') trim(incInfo), &
                    ' @ Iteration ', itmin, '≤',totalIter, '≤', itmax
     if (iand(debug_level(debug_spectral),debug_spectralRotation) /= 0) &
       write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') ' deformation gradient aim (lab) =', &
                             math_transpose33(math_rotate_backward33(F_aim,params%rotation_BC))
     write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') ' deformation gradient aim =', &
                                   math_transpose33(F_aim)
     flush(6)
   endif
 endif newIteration

!--------------------------------------------------------------------------------------------------
! 
 tensorField_real = 0.0_pReal
 do k = 1_pInt, grid3; do j = 1_pInt, grid(2); do i = 1_pInt, grid(1)
   tensorField_real(1:3,1:3,i,j,k) = &
     polarBeta*math_mul3333xx33(C_scale,F(1:3,1:3,i,j,k) - math_I3) -&
     polarAlpha*math_mul33x33(F(1:3,1:3,i,j,k), &
                              math_mul3333xx33(C_scale,F_lambda(1:3,1:3,i,j,k) - math_I3))

 enddo; enddo; enddo
 
!--------------------------------------------------------------------------------------------------
! doing convolution in Fourier space 
 call utilities_FFTtensorForward()
 call utilities_fourierGammaConvolution(math_rotate_backward33(polarBeta*F_aim,params%rotation_BC)) 
 call utilities_FFTtensorBackward()

!--------------------------------------------------------------------------------------------------
! constructing residual                         
 residual_F_lambda = polarBeta*F - tensorField_real(1:3,1:3,1:grid(1),1:grid(2),1:grid3)

!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
 P_avLastEval = P_av
 call Utilities_constitutiveResponse(F_lastInc,F - residual_F_lambda/polarBeta,params%timeinc, &
                                     residual_F,C_volAvg,C_minMaxAvg,P_av,ForwardData,params%rotation_BC)
 call MPI_Allreduce(MPI_IN_PLACE,terminallyIll,1,MPI_LOGICAL,MPI_LOR,PETSC_COMM_WORLD,ierr)
 ForwardData = .False.

!--------------------------------------------------------------------------------------------------
! calculate divergence
 tensorField_real = 0.0_pReal
 tensorField_real(1:3,1:3,1:grid(1),1:grid(2),1:grid3) = residual_F
 call utilities_FFTtensorForward()
 err_div = Utilities_divergenceRMS()
 call utilities_FFTtensorBackward()
 
!--------------------------------------------------------------------------------------------------
! constructing residual
 e = 0_pInt
 do k = 1_pInt, grid3; do j = 1_pInt, grid(2); do i = 1_pInt, grid(1)
   e = e + 1_pInt
   residual_F(1:3,1:3,i,j,k) = math_mul3333xx33(math_invSym3333(materialpoint_dPdF(1:3,1:3,1:3,1:3,1,e) + C_scale), &
                                                residual_F(1:3,1:3,i,j,k) - &
                                                math_mul33x33(F(1:3,1:3,i,j,k), &
                                                math_mul3333xx33(C_scale,F_lambda(1:3,1:3,i,j,k) - math_I3))) &
                                                + residual_F_lambda(1:3,1:3,i,j,k)
 enddo; enddo; enddo
 
!--------------------------------------------------------------------------------------------------
! calculating curl 
 tensorField_real = 0.0_pReal
 tensorField_real(1:3,1:3,1:grid(1),1:grid(2),1:grid3) = F
 call utilities_FFTtensorForward()
 err_curl = Utilities_curlRMS()
 call utilities_FFTtensorBackward()
 
end subroutine AL_formResidual


!--------------------------------------------------------------------------------------------------
!> @brief convergence check
!--------------------------------------------------------------------------------------------------
subroutine AL_converged(snes_local,PETScIter,xnorm,snorm,fnorm,reason,dummy,ierr)
 use numerics, only: &
   itmax, &
   itmin, &
   err_div_tolRel, &
   err_div_tolAbs, &
   err_curl_tolRel, &
   err_curl_tolAbs, &
   err_stress_tolAbs, &
   err_stress_tolRel, &
   worldrank
 use math, only: &
   math_mul3333xx33
 use FEsolving, only: &
   terminallyIll

 implicit none
 SNES :: snes_local
 PetscInt :: PETScIter
 PetscReal :: &
   xnorm, &
   snorm, &
   fnorm
 SNESConvergedReason :: reason
 PetscObject :: dummy
 PetscErrorCode ::ierr
 real(pReal) :: &
   curlTol, &
   divTol, &
   BC_tol
     
!--------------------------------------------------------------------------------------------------
! stress BC handling
 F_aim = F_aim - math_mul3333xx33(S, ((P_av - params%P_BC)))                                        ! S = 0.0 for no bc
 err_BC = maxval(abs((-mask_stress+1.0_pReal)*math_mul3333xx33(C_scale,F_aim-F_av) + &
                         mask_stress              *(P_av - params%P_BC)))                           ! mask = 0.0 for no bc

!--------------------------------------------------------------------------------------------------
! error calculation
 curlTol    = max(maxval(abs(F_aim-math_I3))*err_curl_tolRel,err_curl_tolAbs)
 divTol     = max(maxval(abs(P_av))         *err_div_tolRel,err_div_tolAbs)
 BC_tol     = max(maxval(abs(P_av))    *err_stress_tolrel,err_stress_tolabs)

 converged: if ((totalIter >= itmin .and. &
                           all([ err_div/divTol, &
                                 err_curl/curlTol, &
                                 err_BC/BC_tol       ] < 1.0_pReal)) &
             .or.    terminallyIll) then
   reason = 1
 elseif (totalIter >= itmax) then converged
   reason = -1
 else converged
   reason = 0
 endif converged

!--------------------------------------------------------------------------------------------------
! report
 if (worldrank == 0_pInt) then
   write(6,'(1/,a)') ' ... reporting .............................................................'
   write(6,'(/,a,f12.2,a,es8.2,a,es9.2,a)') ' error curl =       ', &
            err_curl/curlTol,' (',err_curl,' -,   tol =',curlTol,')'
   write(6,'  (a,f12.2,a,es8.2,a,es9.2,a)') ' error divergence = ', &
            err_div/divTol,  ' (',err_div, ' / m, tol =',divTol,')'
   write(6,'  (a,f12.2,a,es8.2,a,es9.2,a)') ' error BC =         ', &
            err_BC/BC_tol, ' (',err_BC,    ' Pa,  tol =',BC_tol,')' 
   write(6,'(/,a)') ' ==========================================================================='
   flush(6) 
 endif

end subroutine AL_converged

!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine AL_forward(guess,timeinc,timeinc_old,loadCaseTime,F_BC,P_BC,rotation_BC)
 use math, only: &
   math_mul33x33, &
   math_mul3333xx33, &
   math_transpose33, &
   math_rotate_backward33
 use numerics, only: &
   worldrank 
 use mesh, only: &
   grid3, &
   grid
 use spectral_utilities, only: &
   Utilities_calculateRate, &
   Utilities_forwardField, &
   Utilities_updateIPcoords, &
   tBoundaryCondition, &
   cutBack
 use IO, only: &
   IO_write_JobRealFile
 use FEsolving, only: &
   restartWrite

 implicit none
 real(pReal), intent(in) :: &
   timeinc_old, &
   timeinc, &
   loadCaseTime                                                                                     !< remaining time of current load case
 type(tBoundaryCondition),      intent(in) :: &
   P_BC, &
   F_BC
 real(pReal), dimension(3,3), intent(in) :: rotation_BC
 logical, intent(in) :: &
   guess
 PetscErrorCode :: ierr
 PetscScalar, dimension(:,:,:,:), pointer :: xx_psc, F, F_lambda
 integer(pInt) :: i, j, k
 real(pReal), dimension(3,3) :: F_lambda33
 character(len=1024) :: rankStr

!--------------------------------------------------------------------------------------------------
! update coordinates and rate and forward last inc
 call DMDAVecGetArrayF90(da,solution_vec,xx_psc,ierr)
 F => xx_psc(0:8,:,:,:)
 F_lambda => xx_psc(9:17,:,:,:)
 if (restartWrite) then
   if (worldrank == 0_pInt) then
     write(6,'(/,a)') ' writing converged results for restart'
     flush(6)
   endif
   write(rankStr,'(a1,i0)')'_',worldrank
   call IO_write_jobRealFile(777,'F'//trim(rankStr),size(F))                                                       ! writing deformation gradient field to file
   write (777,rec=1) F
   close (777)
   call IO_write_jobRealFile(777,'F_lastInc'//trim(rankStr),size(F_lastInc))                                       ! writing F_lastInc field to file
   write (777,rec=1) F_lastInc
   close (777)
   call IO_write_jobRealFile(777,'F_lambda'//trim(rankStr),size(F_lambda))                                         ! writing deformation gradient field to file
   write (777,rec=1) F_lambda
   close (777)
   call IO_write_jobRealFile(777,'F_lambda_lastInc'//trim(rankStr),size(F_lambda_lastInc))                         ! writing F_lastInc field to file
   write (777,rec=1) F_lambda_lastInc
   close (777)
   if (worldrank == 0_pInt) then
     call IO_write_jobRealFile(777,'F_aim',size(F_aim))
     write (777,rec=1) F_aim
     close(777)
     call IO_write_jobRealFile(777,'F_aim_lastInc',size(F_aim_lastInc))
     write (777,rec=1) F_aim_lastInc
     close(777)
     call IO_write_jobRealFile(777,'F_aimDot',size(F_aimDot))
     write (777,rec=1) F_aimDot
     close(777)
     call IO_write_jobRealFile(777,'C_volAvg',size(C_volAvg))
     write (777,rec=1) C_volAvg
     close(777)
     call IO_write_jobRealFile(777,'C_volAvgLastInc',size(C_volAvgLastInc))
     write (777,rec=1) C_volAvgLastInc
     close(777)
   endif
 endif

 call utilities_updateIPcoords(F)

 if (cutBack) then 
   F_aim    = F_aim_lastInc
   F_lambda = reshape(F_lambda_lastInc,[9,grid(1),grid(2),grid3]) 
   F        = reshape(F_lastInc,       [9,grid(1),grid(2),grid3]) 
   C_volAvg = C_volAvgLastInc
 else
   ForwardData = .True.
   C_volAvgLastInc = C_volAvg
!--------------------------------------------------------------------------------------------------
! calculate rate for aim
   if (F_BC%myType=='l') then                                                                       ! calculate f_aimDot from given L and current F
     f_aimDot = F_BC%maskFloat * math_mul33x33(F_BC%values, F_aim)
   elseif(F_BC%myType=='fdot') then                                                                 ! f_aimDot is prescribed
     f_aimDot = F_BC%maskFloat * F_BC%values
   elseif(F_BC%myType=='f') then                                                                    ! aim at end of load case is prescribed
     f_aimDot = F_BC%maskFloat * (F_BC%values -F_aim)/loadCaseTime
   endif
   if (guess) f_aimDot  = f_aimDot + P_BC%maskFloat * (F_aim - F_aim_lastInc)/timeinc_old
   F_aim_lastInc = F_aim

!--------------------------------------------------------------------------------------------------
! update coordinates and rate and forward last inc
   call utilities_updateIPcoords(F)
   Fdot =  Utilities_calculateRate(math_rotate_backward33(f_aimDot,rotation_BC), &
                  timeinc_old,guess,F_lastInc,reshape(F,[3,3,grid(1),grid(2),grid3]))
   F_lambdaDot =  Utilities_calculateRate(math_rotate_backward33(f_aimDot,rotation_BC), &
                  timeinc_old,guess,F_lambda_lastInc,reshape(F_lambda,[3,3,grid(1),grid(2),grid3]))  
   F_lastInc        = reshape(F,       [3,3,grid(1),grid(2),grid3])
   F_lambda_lastInc = reshape(F_lambda,[3,3,grid(1),grid(2),grid3])
 endif

 F_aim = F_aim + f_aimDot * timeinc

!--------------------------------------------------------------------------------------------------
! update local deformation gradient
 F = reshape(Utilities_forwardField(timeinc,F_lastInc,Fdot, &                                       ! ensure that it matches rotated F_aim
                                    math_rotate_backward33(F_aim,rotation_BC)), &
             [9,grid(1),grid(2),grid3])
 F_lambda = reshape(Utilities_forwardField(timeinc,F_lambda_lastInc,F_lambdadot), &
                    [9,grid(1),grid(2),grid3])                                                      ! does not have any average value as boundary condition
 if (.not. guess) then                                                                              ! large strain forwarding
   do k = 1_pInt, grid3; do j = 1_pInt, grid(2); do i = 1_pInt, grid(1)
      F_lambda33 = reshape(F_lambda(1:9,i,j,k),[3,3])
      F_lambda33 = math_mul3333xx33(S_scale,math_mul33x33(F_lambda33, &
                                  math_mul3333xx33(C_scale,&
                                                   math_mul33x33(math_transpose33(F_lambda33),&
                                                                 F_lambda33) -math_I3))*0.5_pReal)&
                              + math_I3
      F_lambda(1:9,i,j,k) = reshape(F_lambda33,[9])
   enddo; enddo; enddo
 endif
 call DMDAVecRestoreArrayF90(da,solution_vec,xx_psc,ierr); CHKERRQ(ierr)

end subroutine AL_forward

!--------------------------------------------------------------------------------------------------
!> @brief destroy routine
!--------------------------------------------------------------------------------------------------
subroutine AL_destroy()
 use spectral_utilities, only: &
   Utilities_destroy

 implicit none
 PetscErrorCode :: ierr

 call VecDestroy(solution_vec,ierr); CHKERRQ(ierr)
 call SNESDestroy(snes,ierr); CHKERRQ(ierr)
 call DMDestroy(da,ierr); CHKERRQ(ierr)

end subroutine AL_destroy

end module spectral_mech_AL
