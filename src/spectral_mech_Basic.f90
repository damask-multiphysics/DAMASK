!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Basic scheme PETSc solver
!--------------------------------------------------------------------------------------------------
module spectral_mech_basic
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
   F_aim_lastInc = math_I3, &
   P_av = 0.0_pReal, &
   F_aimDot = 0.0_pReal
 character(len=1024), private :: incInfo   
 real(pReal), private, dimension(3,3,3,3) :: &
   C_volAvg = 0.0_pReal, &                                                                          !< current volume average stiffness 
   C_volAvgLastInc = 0.0_pReal, &                                                                   !< previous volume average stiffness
   C_minMaxAvg = 0.0_pReal, &                                                                       !< current (min+max)/2 stiffness
   C_minMaxAvgLastInc = 0.0_pReal, &                                                                !< previous (min+max)/2 stiffness
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
   PETScFinalize, &
   MPI_Abort, &
   MPI_Bcast, &
   MPI_Allreduce

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all necessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
subroutine basicPETSc_init
#ifdef __GFORTRAN__
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
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
 use homogenization, only: &
   materialpoint_F0
 use DAMASK_interface, only: &
   getSolverJobName
 use spectral_utilities, only: &
   Utilities_constitutiveResponse, &
   Utilities_updateGamma, &
   Utilities_updateIPcoords, &
   wgt
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
 PetscScalar, pointer, dimension(:,:,:,:)   ::  F

 integer(pInt), dimension(:), allocatable :: localK  
 integer(pInt) :: proc
 character(len=1024) :: rankStr
 
 external :: &
   SNESCreate, &
   SNESSetOptionsPrefix, &
   DMDACreate3D, &
   SNESSetDM, &
   DMCreateGlobalVector, &
   DMDASNESSetFunctionLocal, &
   SNESGetConvergedReason, &
   SNESSetConvergenceTest, &
   SNESSetFromOptions
   
 write(6,'(/,a)') ' <<<+-  DAMASK_spectral_solverBasicPETSc init  -+>>>'
 write(6,'(a15,a)')   ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

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
        grid(1),grid(2),grid(3), &                                                                  ! global grid
        1 , 1, worldsize, &
        9, 0, &                                                                                     ! #dof (F tensor), ghost boundary width (domain overlap)
        grid(1),grid(2),localK, &                                                                   ! local grid
        da,ierr)                                                                                    ! handle, error
 CHKERRQ(ierr)
 call SNESSetDM(snes,da,ierr); CHKERRQ(ierr)
 call DMCreateGlobalVector(da,solution_vec,ierr); CHKERRQ(ierr)                                     ! global solution vector (grid x 9, i.e. every def grad tensor)
 call DMDASNESSetFunctionLocal(da,INSERT_VALUES,BasicPETSC_formResidual,PETSC_NULL_OBJECT,ierr)     ! residual vector of same shape as solution vector
 CHKERRQ(ierr) 
 call SNESSetDM(snes,da,ierr); CHKERRQ(ierr)                                                        ! connect snes to da
 call SNESSetConvergenceTest(snes,BasicPETSC_converged,PETSC_NULL_OBJECT,PETSC_NULL_FUNCTION,ierr)  ! specify custom convergence check function "_converged"
 CHKERRQ(ierr)
 call SNESSetFromOptions(snes,ierr); CHKERRQ(ierr)                                                  ! pull it all together with additional cli arguments

!--------------------------------------------------------------------------------------------------
! init fields                 
 call DMDAVecGetArrayF90(da,solution_vec,F,ierr); CHKERRQ(ierr)                                     ! get the data out of PETSc to work with

 restart: if (restartInc > 1_pInt) then                                                     
   if (iand(debug_level(debug_spectral),debug_spectralRestart)/= 0) &
     write(6,'(/,a,'//IO_intOut(restartInc-1_pInt)//',a)') &
     'reading values of increment ', restartInc - 1_pInt, ' from file'
   flush(6)
   write(rankStr,'(a1,i0)')'_',worldrank
   call IO_read_realFile(777,'F'//trim(rankStr),trim(getSolverJobName()),size(F))
   read (777,rec=1) F; close (777)
   call IO_read_realFile(777,'F_lastInc'//trim(rankStr),trim(getSolverJobName()),size(F_lastInc))
   read (777,rec=1) F_lastInc; close (777)
   call IO_read_realFile(777,'F_aimDot',trim(getSolverJobName()),size(f_aimDot))
   read (777,rec=1) f_aimDot; close (777)
   F_aim         = reshape(sum(sum(sum(F,dim=4),dim=3),dim=2) * wgt, [3,3])                         ! average of F
   F_aim_lastInc = sum(sum(sum(F_lastInc,dim=5),dim=4),dim=3) * wgt                                 ! average of F_lastInc 
 elseif (restartInc == 1_pInt) then restart 
   F_lastInc = spread(spread(spread(math_I3,3,grid(1)),4,grid(2)),5,grid3)                          ! initialize to identity
   F = reshape(F_lastInc,[9,grid(1),grid(2),grid3])
 endif restart

 materialpoint_F0 = reshape(F_lastInc, [3,3,1,product(grid(1:2))*grid3])                            ! set starting condition for materialpoint_stressAndItsTangent
 call Utilities_updateIPcoords(reshape(F,shape(F_lastInc)))
 call Utilities_constitutiveResponse(P, temp33_Real, C_volAvg,C_minMaxAvg, &                        ! stress field, stress avg, global average of stiffness and (min+max)/2
                                     reshape(F,shape(F_lastInc)), &                                 ! target F
                                     0.0_pReal, &                                                   ! time increment
                                     math_I3)                                                       ! no rotation of boundary condition
 call DMDAVecRestoreArrayF90(da,solution_vec,F,ierr); CHKERRQ(ierr)                                 ! write data back to PETSc
                                                                                                    ! QUESTION: why not writing back right after reading (l.189)?

 restartRead: if (restartInc > 1_pInt) then                                                         ! QUESTION: are those values not calc'ed by constitutiveResponse? why reading from file?
   if (iand(debug_level(debug_spectral),debug_spectralRestart)/= 0 .and. worldrank == 0_pInt) &
     write(6,'(/,a,'//IO_intOut(restartInc-1_pInt)//',a)') &
     'reading more values of increment', restartInc-1_pInt, 'from file'
   flush(6)
   call IO_read_realFile(777,'C_volAvg',trim(getSolverJobName()),size(C_volAvg))
   read (777,rec=1) C_volAvg; close (777)
   call IO_read_realFile(777,'C_volAvgLastInc',trim(getSolverJobName()),size(C_volAvgLastInc))
   read (777,rec=1) C_volAvgLastInc; close (777)
   call IO_read_realFile(777,'C_ref',trim(getSolverJobName()),size(C_minMaxAvg))
   read (777,rec=1) C_minMaxAvg; close (777)
 endif restartRead
   
 call Utilities_updateGamma(C_minmaxAvg,.true.)

end subroutine basicPETSc_init
  
!--------------------------------------------------------------------------------------------------
!> @brief solution for the Basic PETSC scheme with internal iterations
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function basicPETSc_solution(incInfoIn,timeinc,timeinc_old,stress_BC,rotation_BC)
 use IO, only: &
   IO_error
 use numerics, only: &
   update_gamma
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
 character(len=*), intent(in) :: &
   incInfoIn
 real(pReal), intent(in) :: &
   timeinc, &                                                                                       !< increment time for current solution
   timeinc_old                                                                                      !< increment time of last successful increment
 type(tBoundaryCondition),      intent(in) :: &
   stress_BC
 real(pReal), dimension(3,3), intent(in) :: rotation_BC
 
!--------------------------------------------------------------------------------------------------
! PETSc Data
 PetscErrorCode :: ierr   
 SNESConvergedReason :: reason

 external :: &
   SNESSolve, &
   SNESGetConvergedReason

 incInfo = incInfoIn

!--------------------------------------------------------------------------------------------------
! update stiffness (and gamma operator)
 S = Utilities_maskedCompliance(rotation_BC,stress_BC%maskLogical,C_volAvg)
 if (update_gamma) call Utilities_updateGamma(C_minmaxAvg,restartWrite)
 

!--------------------------------------------------------------------------------------------------
! set module wide availabe data
 mask_stress        = stress_BC%maskFloat
 params%stress_BC   = stress_BC%values
 params%rotation_BC = rotation_BC
 params%timeinc     = timeinc
 params%timeincOld  = timeinc_old

!--------------------------------------------------------------------------------------------------
! solve BVP 
 call SNESSolve(snes,PETSC_NULL_OBJECT,solution_vec,ierr)
 CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! check convergence
 call SNESGetConvergedReason(snes,reason,ierr); CHKERRQ(ierr)
 
 BasicPETSc_solution%converged = reason > 0
 basicPETSC_solution%iterationsNeeded = totalIter
 basicPETSc_solution%termIll = terminallyIll
 terminallyIll = .false.
 if (reason == -4) call IO_error(893_pInt)                                                         ! MPI error

end function BasicPETSc_solution


!--------------------------------------------------------------------------------------------------
!> @brief forms the basic residual vector
!--------------------------------------------------------------------------------------------------
subroutine BasicPETSC_formResidual(in,x_scal,f_scal,dummy,ierr)
 use numerics, only: &
   itmax, &
   itmin
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
 use spectral_utilities, only: &
   tensorField_real, &
   utilities_FFTtensorForward, &
   utilities_fourierGammaConvolution, &
   utilities_FFTtensorBackward, &
   Utilities_constitutiveResponse, &
   Utilities_divergenceRMS
 use IO, only: &
   IO_intOut 
 use FEsolving, only: &
   terminallyIll

 implicit none
 DMDALocalInfo, dimension(DMDA_LOCAL_INFO_SIZE) :: in
 PetscScalar, &
   dimension(3,3, XG_RANGE,YG_RANGE,ZG_RANGE), intent(in) :: x_scal                         !< what is this?
 PetscScalar, &
   dimension(3,3, X_RANGE,Y_RANGE,Z_RANGE),   intent(out) :: f_scal                         !< what is this?
 PetscInt :: &
   PETScIter, &
   nfuncs
 PetscObject :: dummy
 PetscErrorCode :: ierr
 real(pReal), dimension(3,3) :: &
   deltaF_aim

 external :: &
   SNESGetNumberFunctionEvals, &
   SNESGetIterationNumber

 call SNESGetNumberFunctionEvals(snes,nfuncs,ierr); CHKERRQ(ierr)
 call SNESGetIterationNumber(snes,PETScIter,ierr); CHKERRQ(ierr)

 if (nfuncs == 0 .and. PETScIter == 0) totalIter = -1_pInt                                            ! new increment
!--------------------------------------------------------------------------------------------------
! begin of new iteration
 newIteration: if (totalIter <= PETScIter) then
   totalIter = totalIter + 1_pInt
   write(6,'(1x,a,3(a,'//IO_intOut(itmax)//'))') &
           trim(incInfo), ' @ Iteration ', itmin, '≤',totalIter, '≤', itmax
   if (iand(debug_level(debug_spectral),debug_spectralRotation) /= 0) &
     write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') &
             ' deformation gradient aim (lab) =', math_transpose33(math_rotate_backward33(F_aim,params%rotation_BC))
   write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') &
             ' deformation gradient aim       =', math_transpose33(F_aim)
   flush(6)
 endif newIteration

!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
 call Utilities_constitutiveResponse(f_scal,P_av,C_volAvg,C_minmaxAvg, &
                                     x_scal,params%timeinc, params%rotation_BC)
 call MPI_Allreduce(MPI_IN_PLACE,terminallyIll,1,MPI_LOGICAL,MPI_LOR,PETSC_COMM_WORLD,ierr)
  
!--------------------------------------------------------------------------------------------------
! stress BC handling
 deltaF_aim = math_mul3333xx33(S, P_av - params%stress_BC)
 F_aim = F_aim - deltaF_aim
 err_stress = maxval(abs(mask_stress * (P_av - params%stress_BC)))                                  ! mask = 0.0 when no stress bc

!--------------------------------------------------------------------------------------------------
! updated deformation gradient using fix point algorithm of basic scheme
 tensorField_real = 0.0_pReal
 tensorField_real(1:3,1:3,1:grid(1),1:grid(2),1:grid3) = f_scal
 call utilities_FFTtensorForward()                                                                  ! FFT forward of global "tensorField_real"
 err_div = Utilities_divergenceRMS()                                                                ! divRMS of tensorField_fourier
 call utilities_fourierGammaConvolution(math_rotate_backward33(deltaF_aim,params%rotation_BC))      ! convolution of Gamma and tensorField_fourier, with arg 
 call utilities_FFTtensorBackward()                                                                 ! FFT backward of global tensorField_fourier
 
!--------------------------------------------------------------------------------------------------
! constructing residual
 f_scal = tensorField_real(1:3,1:3,1:grid(1),1:grid(2),1:grid3)                                     ! Gamma*P gives correction towards div(P) = 0, so needs to be zero, too

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
   err_stress_tolAbs
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
 write(6,'(1/,a)') ' ... reporting .............................................................'
 write(6,'(1/,a,f12.2,a,es8.2,a,es9.2,a)') ' error divergence = ', &
         err_div/divTol,  ' (',err_div,' / m, tol =',divTol,')'
 write(6,'(a,f12.2,a,es8.2,a,es9.2,a)')   ' error stress BC =  ', &
         err_stress/stressTol, ' (',err_stress, ' Pa,  tol =',stressTol,')' 
 write(6,'(/,a)') ' ==========================================================================='
flush(6) 
 
end subroutine BasicPETSc_converged

!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!> @details find new boundary conditions and best F estimate for end of current timestep
!> possibly writing restart information, triggering of state increment in DAMASK, and updating of IPcoordinates
!--------------------------------------------------------------------------------------------------
subroutine BasicPETSc_forward(guess,timeinc,timeinc_old,loadCaseTime,deformation_BC,stress_BC,rotation_BC)
  use math, only: &
    math_mul33x33 ,&
    math_rotate_backward33
  use numerics, only: &
    worldrank 
 use homogenization, only: &
   materialpoint_F0
  use mesh, only: &
    grid, &
    grid3
  use CPFEM2, only: &
    CPFEM_age
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
  logical, intent(in) :: &
    guess
  real(pReal), intent(in) :: &
    timeinc_old, &
    timeinc, &
    loadCaseTime                                                                                     !< remaining time of current load case
  type(tBoundaryCondition),    intent(in) :: &
    stress_BC, &
    deformation_BC
  real(pReal), dimension(3,3), intent(in) ::&
    rotation_BC
  PetscErrorCode :: ierr 
  PetscScalar, pointer :: F(:,:,:,:)

  character(len=32) :: rankStr

  call DMDAVecGetArrayF90(da,solution_vec,F,ierr); CHKERRQ(ierr)
  
  if (cutBack) then
    C_volAvg    = C_volAvgLastInc                                                                  ! QUESTION: where is this required?
    C_minMaxAvg = C_minMaxAvgLastInc                                                               ! QUESTION: where is this required?
  else
  !--------------------------------------------------------------------------------------------------
    ! restart information for spectral solver
    if (restartWrite) then                                                                           ! QUESTION: where is this logical properly set?
      write(6,'(/,a)') ' writing converged results for restart'
      flush(6)

      if (worldrank == 0_pInt) then
        call IO_write_jobRealFile(777,'C_volAvg',size(C_volAvg))
        write (777,rec=1) C_volAvg; close(777)
        call IO_write_jobRealFile(777,'C_volAvgLastInc',size(C_volAvgLastInc))
        write (777,rec=1) C_volAvgLastInc; close(777)
        call IO_write_jobRealFile(777,'C_minMaxAvg',size(C_volAvg))
        write (777,rec=1) C_minMaxAvg; close(777)
        call IO_write_jobRealFile(777,'C_minMaxAvgLastInc',size(C_volAvgLastInc))
        write (777,rec=1) C_minMaxAvgLastInc; close(777)
      endif

      write(rankStr,'(a1,i0)')'_',worldrank
      call IO_write_jobRealFile(777,'F'//trim(rankStr),size(F))                                      ! writing deformation gradient field to file
      write (777,rec=1) F; close (777)
      call IO_write_jobRealFile(777,'F_lastInc'//trim(rankStr),size(F_lastInc))                      ! writing F_lastInc field to file
      write (777,rec=1) F_lastInc; close (777)
    endif

    call CPFEM_age()                                                                                 ! age state and kinematics
    call utilities_updateIPcoords(F)

    C_volAvgLastInc    = C_volAvg
    C_minMaxAvgLastInc = C_minMaxAvg

    if (guess) then                                                                                   ! QUESTION: better with a =  L ? x:y
      F_aimDot = stress_BC%maskFloat * (F_aim - F_aim_lastInc)/timeinc_old                            ! initialize with correction based on last inc
    else
      F_aimDot = 0.0_pReal
    endif
    F_aim_lastInc = F_aim
    !--------------------------------------------------------------------------------------------------
    ! calculate rate for aim
    if     (deformation_BC%myType=='l') then                                                          ! calculate f_aimDot from given L and current F
      F_aimDot = &
      F_aimDot + deformation_BC%maskFloat * math_mul33x33(deformation_BC%values, F_aim_lastInc)
    elseif(deformation_BC%myType=='fdot') then                                                        ! f_aimDot is prescribed
      F_aimDot = &
      F_aimDot + deformation_BC%maskFloat * deformation_BC%values
    elseif (deformation_BC%myType=='f') then                                                          ! aim at end of load case is prescribed
      F_aimDot = &
      F_aimDot + deformation_BC%maskFloat * (deformation_BC%values - F_aim_lastInc)/loadCaseTime
    endif


    Fdot =  Utilities_calculateRate(guess, &
                                    F_lastInc,reshape(F,[3,3,grid(1),grid(2),grid3]),timeinc_old, &
                                    math_rotate_backward33(f_aimDot,rotation_BC))
    F_lastInc        = reshape(F,         [3,3,grid(1),grid(2),grid3])                                ! winding F forward
    materialpoint_F0 = reshape(F_lastInc, [3,3,1,product(grid(1:2))*grid3])                           ! set starting condition for materialpoint_stressAndItsTangent
  endif

!--------------------------------------------------------------------------------------------------
! update average and local deformation gradients
  F_aim = F_aim_lastInc + f_aimDot * timeinc
  F = reshape(Utilities_forwardField(timeinc,F_lastInc,Fdot, &                                       ! estimate of F at end of time+timeinc that matches rotated F_aim on average
              math_rotate_backward33(F_aim,rotation_BC)),[9,grid(1),grid(2),grid3])
  call DMDAVecRestoreArrayF90(da,solution_vec,F,ierr); CHKERRQ(ierr)
  
end subroutine BasicPETSc_forward

!--------------------------------------------------------------------------------------------------
!> @brief destroy routine
!--------------------------------------------------------------------------------------------------
subroutine BasicPETSc_destroy()
 use spectral_utilities, only: &
   Utilities_destroy

 implicit none
 PetscErrorCode :: ierr

 external :: &
   VecDestroy, &
   SNESDestroy, &
   DMDestroy

 call VecDestroy(solution_vec,ierr); CHKERRQ(ierr)
 call SNESDestroy(snes,ierr); CHKERRQ(ierr)
 call DMDestroy(da,ierr); CHKERRQ(ierr)

end subroutine BasicPETSc_destroy

end module spectral_mech_basic
