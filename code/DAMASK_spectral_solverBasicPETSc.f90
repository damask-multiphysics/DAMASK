!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Basic scheme PETSc solver
!--------------------------------------------------------------------------------------------------
module DAMASK_spectral_SolverBasicPETSc
 use prec, only: & 
   pInt, &
   pReal
 use math, only: &
   math_I3
 use DAMASK_spectral_Utilities, only: &
   tSolutionState, &
   tSolutionParams

 implicit none
 private
#include <petsc/finclude/petsc.h90>

 character (len=*), parameter, public :: &
   DAMASK_spectral_SolverBasicPETSC_label = 'basicpetsc'
   
!--------------------------------------------------------------------------------------------------
! derived types
 type(tSolutionParams), private :: params

!--------------------------------------------------------------------------------------------------
! PETSc data
 DM,   private :: da
 SNES, private :: snes
 Vec,  private :: solution_vec

!--------------------------------------------------------------------------------------------------
! common pointwise data
 real(pReal), private, dimension(:,:,:,:,:), allocatable ::  F_lastInc, Fdot

!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
 real(pReal), private, dimension(3,3) :: &
   F_aim = math_I3, &
   F_aim_lastIter = math_I3, &
   F_aim_lastInc = math_I3, &
   P_av = 0.0_pReal, &
   F_aimDot=0.0_pReal
 character(len=1024), private :: incInfo   
 real(pReal), private, dimension(3,3,3,3) :: &
   C_volAvg = 0.0_pReal, &                                                                          !< current volume average stiffness 
   C_volAvgLastInc = 0.0_pReal, &                                                                   !< previous volume average stiffness
   C_minMaxAvg = 0.0_pReal, &                                                                       !< current (min+max)/2 stiffness
   S = 0.0_pReal                                                                                    !< current compliance (filled up with zeros)
 real(pReal), private :: err_stress, err_div
 logical, private :: ForwardData
 integer(pInt), private :: &
   totalIter = 0_pInt                                                                               !< total iteration in current increment
 real(pReal), private, dimension(3,3) :: mask_stress = 0.0_pReal

 public :: &
   basicPETSc_init, &
   basicPETSc_solution, &
   BasicPETSc_forward, &
   basicPETSc_destroy
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
subroutine basicPETSc_init
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
 use DAMASK_spectral_Utilities, only: &
   Utilities_constitutiveResponse, &
   Utilities_updateGamma, &
   utilities_updateIPcoords, &
   wgt
 use mesh, only: &
   grid, &
   grid3
 use math, only: &
   math_invSym3333
   
 implicit none
 real(pReal), dimension(:,:,:,:,:), allocatable :: P
 PetscScalar,  dimension(:,:,:,:), pointer     ::  F
 PetscErrorCode :: ierr
 PetscObject    :: dummy
 real(pReal), dimension(3,3) :: &
   temp33_Real = 0.0_pReal
 integer(pInt), dimension(:), allocatable :: localK  
 integer(pInt) :: proc
 character(len=1024) :: rankStr

 mainProcess: if (worldrank == 0_pInt) then
   write(6,'(/,a)') ' <<<+-  DAMASK_spectral_solverBasicPETSc init  -+>>>'
   write(6,'(a)') ' $Id$'
   write(6,'(a15,a)')   ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 allocate (P         (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
!--------------------------------------------------------------------------------------------------
! allocate global fields
 allocate (F_lastInc (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
 allocate (Fdot      (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
    
!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
 call SNESCreate(PETSC_COMM_WORLD,snes,ierr); CHKERRQ(ierr)
 call SNESSetOptionsPrefix(snes,'mech_',ierr);CHKERRQ(ierr) 
 allocate(localK(worldsize), source = 0); localK(worldrank+1) = grid3
 do proc = 1, worldsize
   call MPI_Bcast(localK(proc),1,MPI_INTEGER,proc-1,PETSC_COMM_WORLD,ierr)
 enddo  
 call DMDACreate3d(PETSC_COMM_WORLD, &
        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &                                     ! cut off stencil at boundary
        DMDA_STENCIL_BOX, &                                                                         ! Moore (26) neighborhood around central point
        grid(1),grid(2),grid(3), &                                                ! global grid
        1, 1, worldsize, &
        9, 0, &                                                                                     ! #dof (F tensor), ghost boundary width (domain overlap)
        grid (1),grid (2),localK, &                                                       ! local grid
        da,ierr)                                                                                    ! handle, error
 CHKERRQ(ierr)
 call SNESSetDM(snes,da,ierr); CHKERRQ(ierr)
 call DMCreateGlobalVector(da,solution_vec,ierr); CHKERRQ(ierr)                                     ! global solution vector (grid x 9, i.e. every def grad tensor)
 call DMDASNESSetFunctionLocal(da,INSERT_VALUES,BasicPETSC_formResidual,dummy,ierr)                 ! residual vector of same shape as solution vector
 CHKERRQ(ierr) 
 call SNESSetDM(snes,da,ierr); CHKERRQ(ierr)                                                        ! connect snes to da
 call SNESSetConvergenceTest(snes,BasicPETSC_converged,dummy,PETSC_NULL_FUNCTION,ierr)              ! specify custom convergence check function "_converged"
 CHKERRQ(ierr)
 call SNESSetFromOptions(snes,ierr); CHKERRQ(ierr)                                                  ! pull it all together with additional cli arguments

!--------------------------------------------------------------------------------------------------
! init fields                 
 call DMDAVecGetArrayF90(da,solution_vec,F,ierr); CHKERRQ(ierr)                                     ! get the data out of PETSc to work with

 restart: if (restartInc > 1_pInt) then                                                     
   if (iand(debug_level(debug_spectral),debug_spectralRestart)/= 0 .and. worldrank == 0_pInt) &
     write(6,'(/,a,'//IO_intOut(restartInc-1_pInt)//',a)') &
     'reading values of increment ', restartInc - 1_pInt, ' from file'
   flush(6)
   write(rankStr,'(a1,i0)')'_',worldrank
   call IO_read_realFile(777,'F'//trim(rankStr),trim(getSolverJobName()),size(F))
   read (777,rec=1) F
   close (777)
   call IO_read_realFile(777,'F_lastInc'//trim(rankStr),trim(getSolverJobName()),size(F_lastInc))
   read (777,rec=1) F_lastInc
   close (777)
   call IO_read_realFile(777,'F_aimDot',trim(getSolverJobName()),size(f_aimDot))
   read (777,rec=1) f_aimDot
   close (777)
   F_aim         = reshape(sum(sum(sum(F,dim=4),dim=3),dim=2) * wgt, [3,3])                         ! average of F
   F_aim_lastInc = sum(sum(sum(F_lastInc,dim=5),dim=4),dim=3) * wgt                                 ! average of F_lastInc 
 elseif (restartInc == 1_pInt) then restart 
   F_lastInc = spread(spread(spread(math_I3,3,grid(1)),4,grid(2)),5,grid3)                          ! initialize to identity
   F = reshape(F_lastInc,[9,grid(1),grid(2),grid3])
 endif restart

 call Utilities_updateIPcoords(F)
 call Utilities_constitutiveResponse(F_lastInc, F, &
    0.0_pReal, &
    P, &
    C_volAvg,C_minMaxAvg, &                                                                         ! global average of stiffness and (min+max)/2
    temp33_Real, &
    .false., &
    math_I3)

 call DMDAVecRestoreArrayF90(da,solution_vec,F,ierr); CHKERRQ(ierr)                                 ! write data back to PETSc

 restartRead: if (restartInc > 1_pInt) then                                                    
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
 endif restartRead
   
 call Utilities_updateGamma(C_minmaxAvg,.True.)

end subroutine basicPETSc_init
  
!--------------------------------------------------------------------------------------------------
!> @brief solution for the Basic PETSC scheme with internal iterations
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function basicPETSc_solution( &
     incInfoIn,guess,timeinc,timeinc_old,loadCaseTime,P_BC,F_BC,rotation_BC)
 use numerics, only: &
   update_gamma, &
   itmax
 use DAMASK_spectral_Utilities, only: &
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
 type(tBoundaryCondition),      intent(in) :: &
   P_BC, &
   F_BC
 character(len=*), intent(in) :: &
   incInfoIn
 real(pReal), dimension(3,3), intent(in) :: rotation_BC
 logical, intent(in) :: &
   guess
 
!--------------------------------------------------------------------------------------------------
! PETSc Data
 PetscErrorCode :: ierr   
 SNESConvergedReason :: reason
 incInfo = incInfoIn

!--------------------------------------------------------------------------------------------------
! update stiffness (and gamma operator)
 S = Utilities_maskedCompliance(rotation_BC,P_BC%maskLogical,C_volAvg)
 if (update_gamma) call Utilities_updateGamma(C_minmaxAvg,restartWrite)
 
 BasicPETSc_solution%converged =.false.
 
!--------------------------------------------------------------------------------------------------
! set module wide availabe data 
 mask_stress = P_BC%maskFloat
 params%P_BC = P_BC%values
 params%rotation_BC = rotation_BC
 params%timeinc = timeinc
 params%timeincOld = timeinc_old

 call SNESSolve(snes,PETSC_NULL_OBJECT,solution_vec,ierr); CHKERRQ(ierr)
 call SNESGetConvergedReason(snes,reason,ierr); CHKERRQ(ierr)
 basicPETSc_solution%termIll = terminallyIll
 terminallyIll = .false.

 if (reason < 1) then
   basicPETSC_solution%converged = .false.
   basicPETSC_solution%iterationsNeeded = itmax
 else
   basicPETSC_solution%converged = .true.
   basicPETSC_solution%iterationsNeeded = totalIter
 endif

end function BasicPETSc_solution


!--------------------------------------------------------------------------------------------------
!> @brief forms the AL residual vector
!--------------------------------------------------------------------------------------------------
subroutine BasicPETSC_formResidual(in,x_scal,f_scal,dummy,ierr)
 use numerics, only: &
   itmax, &
   itmin
 use numerics, only: &
   worldrank
 use mesh, only: &
   grid, &
   grid3
 use math, only: &
   math_rotate_backward33, &
   math_transpose33, &
   math_mul3333xx33
 use debug, only: &
   debug_level, &
   debug_spectral, &
   debug_spectralRotation
 use DAMASK_spectral_Utilities, only: &
   tensorField_real, &
   utilities_FFTtensorForward, &
   utilities_FFTtensorBackward, &
   utilities_fourierGammaConvolution, &
   Utilities_inverseLaplace, &
   Utilities_constitutiveResponse, &
   Utilities_divergenceRMS
 use IO, only: &
   IO_intOut 
 use FEsolving, only: &
   terminallyIll

 implicit none
 DMDALocalInfo, dimension(DMDA_LOCAL_INFO_SIZE) :: &
   in
 PetscScalar, dimension(3,3, &
   XG_RANGE,YG_RANGE,ZG_RANGE) :: &
   x_scal
 PetscScalar, dimension(3,3, &
   X_RANGE,Y_RANGE,Z_RANGE) :: &
   f_scal
 PetscInt :: &
   PETScIter, &
   nfuncs
 PetscObject :: dummy
 PetscErrorCode :: ierr

 call SNESGetNumberFunctionEvals(snes,nfuncs,ierr); CHKERRQ(ierr)
 call SNESGetIterationNumber(snes,PETScIter,ierr); CHKERRQ(ierr)

 if(nfuncs== 0 .and. PETScIter == 0) totalIter = -1_pInt                                            ! new increment
 newIteration: if (totalIter <= PETScIter) then
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
! evaluate constitutive response
 call Utilities_constitutiveResponse(F_lastInc,x_scal,params%timeinc, &
                                     f_scal,C_volAvg,C_minmaxAvg,P_av,ForwardData,params%rotation_BC)
 call MPI_Allreduce(MPI_IN_PLACE,terminallyIll,1,MPI_LOGICAL,MPI_LOR,PETSC_COMM_WORLD,ierr)
 ForwardData = .false.
  
!--------------------------------------------------------------------------------------------------
! stress BC handling
 F_aim_lastIter = F_aim
 F_aim = F_aim - math_mul3333xx33(S, ((P_av - params%P_BC)))                                        ! S = 0.0 for no bc
 err_stress = maxval(abs(mask_stress * (P_av - params%P_BC)))                                       ! mask = 0.0 for no bc
 
!--------------------------------------------------------------------------------------------------
! updated deformation gradient using fix point algorithm of basic scheme
 tensorField_real = 0.0_pReal
 tensorField_real(1:3,1:3,1:grid(1),1:grid(2),1:grid3) = f_scal
 call utilities_FFTtensorForward()
 err_div = Utilities_divergenceRMS()
 call utilities_fourierGammaConvolution(math_rotate_backward33(F_aim_lastIter-F_aim,params%rotation_BC))
 call utilities_FFTtensorBackward()
 
!--------------------------------------------------------------------------------------------------
! constructing residual
 f_scal = tensorField_real(1:3,1:3,1:grid(1),1:grid(2),1:grid3)

end subroutine BasicPETSc_formResidual


!--------------------------------------------------------------------------------------------------
!> @brief convergence check
!--------------------------------------------------------------------------------------------------
subroutine BasicPETSc_converged(snes_local,PETScIter,xnorm,snorm,fnorm,reason,dummy,ierr)
 use numerics, only: &
   itmax, &
   itmin, &
   err_div_tolRel, &
   err_div_tolAbs, &
   err_stress_tolRel, &
   err_stress_tolAbs, &
   worldrank
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
 PetscErrorCode :: ierr
 real(pReal) :: &
   divTol, &
   stressTol 
 
 divTol    = max(maxval(abs(P_av))*err_div_tolRel,err_div_tolAbs)
 stressTol = max(maxval(abs(P_av))*err_stress_tolrel,err_stress_tolabs)
 
 converged: if ((totalIter >= itmin .and. &
                           all([ err_div/divTol, &
                                 err_stress/stressTol       ] < 1.0_pReal)) &
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
   write(6,'(1/,a,f12.2,a,es8.2,a,es9.2,a)') ' error divergence = ', &
            err_div/divTol,  ' (',err_div,' / m, tol =',divTol,')'
   write(6,'(a,f12.2,a,es8.2,a,es9.2,a)')   ' error stress BC =  ', &
            err_stress/stressTol, ' (',err_stress, ' Pa,  tol =',stressTol,')' 
   write(6,'(/,a)') ' ==========================================================================='
   flush(6) 
 endif
 
end subroutine BasicPETSc_converged

!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine BasicPETSc_forward(guess,timeinc,timeinc_old,loadCaseTime,F_BC,P_BC,rotation_BC)
 use math, only: &
   math_mul33x33 ,&
   math_rotate_backward33
 use mesh, only: &
   grid, &
   grid3
 use DAMASK_spectral_Utilities, only: &
   Utilities_calculateRate, &
   Utilities_forwardField, &
   utilities_updateIPcoords, &
   tBoundaryCondition, &
   cutBack
 use IO, only: &
   IO_write_JobRealFile
 use FEsolving, only: &
   restartWrite
 use numerics, only: &
   worldrank

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
 PetscScalar, pointer :: F(:,:,:,:)
 PetscErrorCode :: ierr
 character(len=1024) :: rankStr

 call DMDAVecGetArrayF90(da,solution_vec,F,ierr)
!--------------------------------------------------------------------------------------------------
! restart information for spectral solver
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
   if (worldrank == 0_pInt) then
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
   F_aim = F_aim_lastInc
   F    = reshape(F_lastInc,    [9,grid(1),grid(2),grid3]) 
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
   F_lastInc     = reshape(F,       [3,3,grid(1),grid(2),grid3])
 endif
 F_aim = F_aim + f_aimDot * timeinc

!--------------------------------------------------------------------------------------------------
! update local deformation gradient
 F     = reshape(Utilities_forwardField(timeinc,F_lastInc,Fdot, &                                   ! ensure that it matches rotated F_aim
           math_rotate_backward33(F_aim,rotation_BC)),[9,grid(1),grid(2),grid3])
 call DMDAVecRestoreArrayF90(da,solution_vec,F,ierr); CHKERRQ(ierr)

end subroutine BasicPETSc_forward

!--------------------------------------------------------------------------------------------------
!> @brief destroy routine
!--------------------------------------------------------------------------------------------------
subroutine BasicPETSc_destroy()
 use DAMASK_spectral_Utilities, only: &
   Utilities_destroy

 implicit none
 PetscErrorCode :: ierr

 call VecDestroy(solution_vec,ierr); CHKERRQ(ierr)
 call SNESDestroy(snes,ierr); CHKERRQ(ierr)
 call DMDestroy(da,ierr); CHKERRQ(ierr)

end subroutine BasicPETSc_destroy

end module DAMASK_spectral_SolverBasicPETSc
