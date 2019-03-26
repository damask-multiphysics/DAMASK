!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Grid solver for mechanics: Spectral Polarisation
!--------------------------------------------------------------------------------------------------
module grid_mech_spectral_polarisation
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscdmda.h>
 use PETScdmda
 use PETScsnes
 use prec, only: & 
   pReal
 use math, only: &
   math_I3
 use spectral_utilities, only: &
   tSolutionState, &
   tSolutionParams

 implicit none
 private
 
!--------------------------------------------------------------------------------------------------
! derived types
  type(tSolutionParams), private :: params
  
  type, private :: tNumerics
    logical :: &
      update_gamma !< update gamma operator with current stiffness
  end type tNumerics
  
  type(tNumerics) :: num                                                                            ! numerics parameters. Better name?
  
!--------------------------------------------------------------------------------------------------
! PETSc data
 DM,   private :: da
 SNES, private :: snes
 Vec,  private :: solution_vec

!--------------------------------------------------------------------------------------------------
! common pointwise data
 real(pReal), private, dimension(:,:,:,:,:), allocatable :: &
   F_lastInc, &                                                                                     !< field of previous compatible deformation gradients
   F_tau_lastInc, &                                                                                 !< field of previous incompatible deformation gradient 
   Fdot, &                                                                                          !< field of assumed rate of compatible deformation gradient
   F_tauDot                                                                                         !< field of assumed rate of incopatible deformation gradient

!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
 real(pReal), private, dimension(3,3) :: &
   F_aimDot = 0.0_pReal, &                                                                          !< assumed rate of average deformation gradient
   F_aim = math_I3, &                                                                               !< current prescribed deformation gradient
   F_aim_lastInc = math_I3, &                                                                       !< previous average deformation gradient
   F_av = 0.0_pReal, &                                                                              !< average incompatible def grad field
   P_av = 0.0_pReal                                                                                 !< average 1st Piola--Kirchhoff stress
 
 character(len=1024), private :: incInfo                                                            !< time and increment information
 real(pReal), private, dimension(3,3,3,3) :: &
   C_volAvg = 0.0_pReal, &                                                                          !< current volume average stiffness 
   C_volAvgLastInc = 0.0_pReal, &                                                                   !< previous volume average stiffness
   C_minMaxAvg = 0.0_pReal, &                                                                       !< current (min+max)/2 stiffness
   C_minMaxAvgLastInc = 0.0_pReal, &                                                                !< previous (min+max)/2 stiffness
   S = 0.0_pReal, &                                                                                 !< current compliance (filled up with zeros)
   C_scale = 0.0_pReal, &                             
   S_scale = 0.0_pReal

 real(pReal), private :: &
   err_BC, &                                                                                        !< deviation from stress BC
   err_curl, &                                                                                      !< RMS of curl of F
   err_div                                                                                          !< RMS of div of P
  
 integer, private :: &
   totalIter = 0                                                                                    !< total iteration in current increment
 
 public :: &
   grid_mech_spectral_polarisation_init, &
   grid_mech_spectral_polarisation_solution, &
   grid_mech_spectral_polarisation_forward
 private :: &
   converged, &
   formResidual

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all necessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
subroutine grid_mech_spectral_polarisation_init
 use IO, only: &
   IO_intOut, &
   IO_error, &
   IO_open_jobFile_binary
 use FEsolving, only: &
   restartInc
 use config, only :&
   config_numerics
 use numerics, only: &
   worldrank, &
   worldsize, &
   petsc_options
 use homogenization, only: &
   materialpoint_F0
 use DAMASK_interface, only: &
   getSolverJobName
 use spectral_utilities, only: &
   utilities_constitutiveResponse, &
   utilities_updateGamma, &
   utilities_updateIPcoords, &
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
 PetscScalar, pointer, dimension(:,:,:,:) :: &
   FandF_tau, &                                                                                     ! overall pointer to solution data
   F, &                                                                                             ! specific (sub)pointer
   F_tau                                                                                            ! specific (sub)pointer
 PetscInt, dimension(worldsize) :: localK 
 integer :: fileUnit
 character(len=1024) :: rankStr
 
 write(6,'(/,a)') ' <<<+-  grid_mech_spectral_polarisation init  -+>>>'

 write(6,'(/,a)') ' Shanthraj et al., International Journal of Plasticity 66:31–45, 2015'
 write(6,'(a)')   ' https://doi.org/10.1016/j.ijplas.2014.02.006'
 
 num%update_gamma = config_numerics%getInt('update_gamma',defaultVal=0) > 0

!--------------------------------------------------------------------------------------------------
! set default and user defined options for PETSc
 call PETScOptionsInsertString(PETSC_NULL_OPTIONS,'-mech_snes_type ngmres',ierr)
 CHKERRQ(ierr)
 call PETScOptionsInsertString(PETSC_NULL_OPTIONS,trim(petsc_options),ierr)
 CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! allocate global fields
 allocate (F_lastInc    (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
 allocate (Fdot         (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
 allocate (F_tau_lastInc(3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
 allocate (F_tauDot     (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
    
!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
 call SNESCreate(PETSC_COMM_WORLD,snes,ierr); CHKERRQ(ierr)
 call SNESSetOptionsPrefix(snes,'mech_',ierr);CHKERRQ(ierr) 
 localK              = 0
 localK(worldrank+1) = grid3
 call MPI_Allreduce(MPI_IN_PLACE,localK,worldsize,MPI_INTEGER,MPI_SUM,PETSC_COMM_WORLD,ierr)
 call DMDACreate3d(PETSC_COMM_WORLD, &
        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &                                     ! cut off stencil at boundary
        DMDA_STENCIL_BOX, &                                                                         ! Moore (26) neighborhood around central point
        grid(1),grid(2),grid(3), &                                                                  ! global grid
        1 , 1, worldsize, &
        18, 0, &                                                                                    ! #dof (F tensor), ghost boundary width (domain overlap)
        [grid(1)],[grid(2)],localK, &                                                               ! local grid
        da,ierr)                                                                                    ! handle, error
 CHKERRQ(ierr)
 call SNESSetDM(snes,da,ierr); CHKERRQ(ierr)                                                        ! connect snes to da
 call DMsetFromOptions(da,ierr); CHKERRQ(ierr)
 call DMsetUp(da,ierr); CHKERRQ(ierr)
 call DMcreateGlobalVector(da,solution_vec,ierr); CHKERRQ(ierr)                                     ! global solution vector (grid x 18, i.e. every def grad tensor)
 call DMDASNESsetFunctionLocal(da,INSERT_VALUES,formResidual,PETSC_NULL_SNES,ierr)     ! residual vector of same shape as solution vector
 CHKERRQ(ierr) 
 call SNESsetConvergenceTest(snes,converged,PETSC_NULL_SNES,PETSC_NULL_FUNCTION,ierr)  ! specify custom convergence check function "converged"
 CHKERRQ(ierr)
 call SNESsetFromOptions(snes,ierr); CHKERRQ(ierr)                                                  ! pull it all together with additional CLI arguments

!--------------------------------------------------------------------------------------------------
! init fields                 
 call DMDAVecGetArrayF90(da,solution_vec,FandF_tau,ierr); CHKERRQ(ierr)                             ! places pointer on PETSc data
 F        => FandF_tau( 0: 8,:,:,:)
 F_tau    => FandF_tau( 9:17,:,:,:)

 restart: if (restartInc > 0) then
   write(6,'(/,a,'//IO_intOut(restartInc)//',a)') ' reading values of increment ', restartInc, ' from file'

   fileUnit = IO_open_jobFile_binary('F_aim')
   read(fileUnit) F_aim; close(fileUnit)
   fileUnit = IO_open_jobFile_binary('F_aim_lastInc')
   read(fileUnit) F_aim_lastInc; close(fileUnit)
   fileUnit = IO_open_jobFile_binary('F_aimDot')
   read(fileUnit) F_aimDot; close(fileUnit)

   write(rankStr,'(a1,i0)')'_',worldrank

   fileUnit = IO_open_jobFile_binary('F'//trim(rankStr))
   read(fileUnit) F; close (fileUnit)
   fileUnit = IO_open_jobFile_binary('F_lastInc'//trim(rankStr))
   read(fileUnit) F_lastInc; close (fileUnit)
   fileUnit = IO_open_jobFile_binary('F_tau'//trim(rankStr))
   read(fileUnit) F_tau; close (fileUnit)
   fileUnit = IO_open_jobFile_binary('F_tau_lastInc'//trim(rankStr))
   read(fileUnit) F_tau_lastInc; close (fileUnit)

 elseif (restartInc == 0) then restart
   F_lastInc = spread(spread(spread(math_I3,3,grid(1)),4,grid(2)),5,grid3)                          ! initialize to identity
   F = reshape(F_lastInc,[9,grid(1),grid(2),grid3])
   F_tau = 2.0_pReal*F
   F_tau_lastInc = 2.0_pReal*F_lastInc
 endif restart

 materialpoint_F0 = reshape(F_lastInc, [3,3,1,product(grid(1:2))*grid3])                            ! set starting condition for materialpoint_stressAndItsTangent
 call Utilities_updateIPcoords(reshape(F,shape(F_lastInc)))
 call Utilities_constitutiveResponse(P,temp33_Real,C_volAvg,C_minMaxAvg, &                          ! stress field, stress avg, global average of stiffness and (min+max)/2
                                     reshape(F,shape(F_lastInc)), &                                 ! target F
                                     0.0_pReal, &                                                   ! time increment
                                     math_I3)                                                       ! no rotation of boundary condition
 call DMDAVecRestoreArrayF90(da,solution_vec,FandF_tau,ierr); CHKERRQ(ierr)                         ! deassociate pointer

 restartRead: if (restartInc > 0) then
   write(6,'(/,a,'//IO_intOut(restartInc)//',a)') ' reading more values of increment ', restartInc, ' from file'
   fileUnit = IO_open_jobFile_binary('C_volAvg')
   read(fileUnit) C_volAvg; close(fileUnit)
   fileUnit = IO_open_jobFile_binary('C_volAvgLastInv')
   read(fileUnit) C_volAvgLastInc; close(fileUnit)
   fileUnit = IO_open_jobFile_binary('C_ref')
   read(fileUnit) C_minMaxAvg; close(fileUnit)
 endif restartRead

 call Utilities_updateGamma(C_minMaxAvg,.true.)
 C_scale = C_minMaxAvg
 S_scale = math_invSym3333(C_minMaxAvg)
 
end subroutine grid_mech_spectral_polarisation_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the Polarisation scheme with internal iterations
!--------------------------------------------------------------------------------------------------
function grid_mech_spectral_polarisation_solution(incInfoIn,timeinc,timeinc_old,stress_BC,rotation_BC) result(solution)
 use math, only: &
   math_invSym3333
 use spectral_utilities, only: &
   tBoundaryCondition, &
    utilities_maskedCompliance, &
    utilities_updateGamma
 use FEsolving, only: &
   restartWrite, &
   terminallyIll

 implicit none

!--------------------------------------------------------------------------------------------------
! input data for solution
 character(len=*), intent(in) :: &
   incInfoIn
 real(pReal), intent(in) :: &
   timeinc, &                                                                                       !< time increment of current solution
   timeinc_old                                                                                      !< time increment of last successful increment
 type(tBoundaryCondition),    intent(in) :: &
   stress_BC
 real(pReal), dimension(3,3), intent(in) :: rotation_BC
 type(tSolutionState)                    :: &
   solution
!--------------------------------------------------------------------------------------------------
! PETSc Data
 PetscErrorCode :: ierr   
 SNESConvergedReason :: reason

 incInfo = incInfoIn

!--------------------------------------------------------------------------------------------------
! update stiffness (and gamma operator)
 S = Utilities_maskedCompliance(rotation_BC,stress_BC%maskLogical,C_volAvg)
 if (num%update_gamma) then
   call utilities_updateGamma(C_minMaxAvg,restartWrite)
   C_scale = C_minMaxAvg
   S_scale = math_invSym3333(C_minMaxAvg)
 endif  

!--------------------------------------------------------------------------------------------------
! set module wide available data 
 params%stress_mask = stress_BC%maskFloat
 params%stress_BC   = stress_BC%values
 params%rotation_BC = rotation_BC
 params%timeinc     = timeinc
 params%timeincOld  = timeinc_old

!--------------------------------------------------------------------------------------------------
! solve BVP 
 call SNESsolve(snes,PETSC_NULL_VEC,solution_vec,ierr); CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! check convergence
 call SNESGetConvergedReason(snes,reason,ierr); CHKERRQ(ierr)

 solution%converged = reason > 0
 solution%iterationsNeeded = totalIter
 solution%termIll = terminallyIll
 terminallyIll = .false.

end function grid_mech_spectral_polarisation_solution


!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!> @details find new boundary conditions and best F estimate for end of current timestep
!> possibly writing restart information, triggering of state increment in DAMASK, and updating of IPcoordinates
!--------------------------------------------------------------------------------------------------
subroutine grid_mech_spectral_polarisation_forward(guess,timeinc,timeinc_old,loadCaseTime,deformation_BC,stress_BC,rotation_BC)
 use math, only: &
   math_mul33x33, &
   math_mul3333xx33, &
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
    utilities_calculateRate, &
    utilities_forwardField, &
    utilities_updateIPcoords, &
    tBoundaryCondition, &
    cutBack
  use IO, only: &
    IO_open_jobFile_binary
  use FEsolving, only: &
    restartWrite

 implicit none
  logical, intent(in) :: &
    guess
  real(pReal), intent(in) :: &
    timeinc_old, &
    timeinc, &
    loadCaseTime                                                                                     !< remaining time of current load case
  type(tBoundaryCondition),      intent(in) :: &
    stress_BC, &
    deformation_BC
  real(pReal), dimension(3,3), intent(in) ::&
    rotation_BC
  PetscErrorCode :: ierr
  PetscScalar, dimension(:,:,:,:), pointer :: FandF_tau, F, F_tau
  integer :: i, j, k
  real(pReal), dimension(3,3) :: F_lambda33

  integer :: fileUnit
  character(len=32) :: rankStr

  call DMDAVecGetArrayF90(da,solution_vec,FandF_tau,ierr); CHKERRQ(ierr)
  F        => FandF_tau( 0: 8,:,:,:)
  F_tau    => FandF_tau( 9:17,:,:,:)

  if (cutBack) then
    C_volAvg    = C_volAvgLastInc                                                                  ! QUESTION: where is this required?
    C_minMaxAvg = C_minMaxAvgLastInc                                                               ! QUESTION: where is this required?
  else
  !--------------------------------------------------------------------------------------------------
    ! restart information for spectral solver
    if (restartWrite) then                                                                           ! QUESTION: where is this logical properly set?
      write(6,'(/,a)') ' writing converged results for restart'
      flush(6)

      if (worldrank == 0) then
        fileUnit = IO_open_jobFile_binary('C_volAvg','w')
        write(fileUnit) C_volAvg; close(fileUnit)
        fileUnit = IO_open_jobFile_binary('C_volAvgLastInv','w')
        write(fileUnit) C_volAvgLastInc; close(fileUnit)
        fileUnit = IO_open_jobFile_binary('F_aim','w')
        write(fileUnit) F_aim; close(fileUnit)
        fileUnit = IO_open_jobFile_binary('F_aim_lastInc','w')
        write(fileUnit) F_aim_lastInc; close(fileUnit)
        fileUnit = IO_open_jobFile_binary('F_aimDot','w')
        write(fileUnit) F_aimDot; close(fileUnit)
      endif

      write(rankStr,'(a1,i0)')'_',worldrank
      fileUnit = IO_open_jobFile_binary('F'//trim(rankStr),'w')
      write(fileUnit) F; close (fileUnit)
      fileUnit = IO_open_jobFile_binary('F_lastInc'//trim(rankStr),'w')
      write(fileUnit) F_lastInc; close (fileUnit)
      fileUnit = IO_open_jobFile_binary('F_tau'//trim(rankStr),'w')
      write(fileUnit) F_tau; close (fileUnit)
      fileUnit = IO_open_jobFile_binary('F_tau_lastInc'//trim(rankStr),'w')
      write(fileUnit) F_tau_lastInc; close (fileUnit)
    endif

    call CPFEM_age                                                                                  ! age state and kinematics
    call utilities_updateIPcoords(F)

    C_volAvgLastInc    = C_volAvg
    C_minMaxAvgLastInc = C_minMaxAvg

    F_aimDot = merge(stress_BC%maskFloat*(F_aim-F_aim_lastInc)/timeinc_old, 0.0_pReal, guess)
    F_aim_lastInc = F_aim

    !--------------------------------------------------------------------------------------------------
    ! calculate rate for aim
    if     (deformation_BC%myType=='l') then                                                          ! calculate F_aimDot from given L and current F
      F_aimDot = &
      F_aimDot + deformation_BC%maskFloat * math_mul33x33(deformation_BC%values, F_aim_lastInc)
    elseif(deformation_BC%myType=='fdot') then                                                        ! F_aimDot is prescribed
      F_aimDot = &
      F_aimDot + deformation_BC%maskFloat * deformation_BC%values
    elseif (deformation_BC%myType=='f') then                                                          ! aim at end of load case is prescribed
      F_aimDot = &
      F_aimDot + deformation_BC%maskFloat * (deformation_BC%values - F_aim_lastInc)/loadCaseTime
    endif


    Fdot        = utilities_calculateRate(guess, &
                                          F_lastInc,reshape(F,[3,3,grid(1),grid(2),grid3]),timeinc_old, &
                                          math_rotate_backward33(F_aimDot,rotation_BC))
    F_tauDot    = utilities_calculateRate(guess, &
                                          F_tau_lastInc,reshape(F_tau,[3,3,grid(1),grid(2),grid3]), timeinc_old, &
                                          math_rotate_backward33(F_aimDot,rotation_BC))
    F_lastInc        = reshape(F,         [3,3,grid(1),grid(2),grid3])                                ! winding F forward
    F_tau_lastInc    = reshape(F_tau,     [3,3,grid(1),grid(2),grid3])                                ! winding F_tau forward
    materialpoint_F0 = reshape(F_lastInc, [3,3,1,product(grid(1:2))*grid3])                           ! set starting condition for materialpoint_stressAndItsTangent
  endif

!--------------------------------------------------------------------------------------------------
! update average and local deformation gradients
  F_aim = F_aim_lastInc + F_aimDot * timeinc
  F = reshape(utilities_forwardField(timeinc,F_lastInc,Fdot, &                                        ! estimate of F at end of time+timeinc that matches rotated F_aim on average
                                     math_rotate_backward33(F_aim,rotation_BC)),&
              [9,grid(1),grid(2),grid3])
 if (guess) then
    F_tau = reshape(Utilities_forwardField(timeinc,F_tau_lastInc,F_taudot), &
                    [9,grid(1),grid(2),grid3])                                                        ! does not have any average value as boundary condition
  else
   do k = 1, grid3; do j = 1, grid(2); do i = 1, grid(1)
      F_lambda33 = reshape(F_tau(1:9,i,j,k)-F(1:9,i,j,k),[3,3])
      F_lambda33 = math_mul3333xx33(S_scale,math_mul33x33(F_lambda33, &
                                  math_mul3333xx33(C_scale,&
                                                   math_mul33x33(transpose(F_lambda33),&
                                                                 F_lambda33)-math_I3))*0.5_pReal)&
                              + math_I3
      F_tau(1:9,i,j,k) = reshape(F_lambda33,[9])+F(1:9,i,j,k)
   enddo; enddo; enddo
 endif

 call DMDAVecRestoreArrayF90(da,solution_vec,FandF_tau,ierr); CHKERRQ(ierr)

end subroutine grid_mech_spectral_polarisation_forward


!--------------------------------------------------------------------------------------------------
!> @brief convergence check
!--------------------------------------------------------------------------------------------------
subroutine converged(snes_local,PETScIter,xnorm,snorm,fnorm,reason,dummy,ierr)
 use numerics, only: &
   itmax, &
   itmin, &
   err_div_tolRel, &
   err_div_tolAbs, &
   err_curl_tolRel, &
   err_curl_tolAbs, &
   err_stress_tolRel, &
   err_stress_tolAbs
 use FEsolving, only: &
   terminallyIll

 implicit none
 SNES :: snes_local
 PetscInt :: PETScIter
 PetscReal :: &
   xnorm, &                                                                                        ! not used
   snorm, &                                                                                        ! not used
   fnorm                                                                                           ! not used
 SNESConvergedReason :: reason
 PetscObject :: dummy
 PetscErrorCode :: ierr
 real(pReal) :: &
   curlTol, &
   divTol, &
   BCTol

 curlTol    = max(maxval(abs(F_aim-math_I3))*err_curl_tolRel  ,err_curl_tolAbs)
 divTol     = max(maxval(abs(P_av))         *err_div_tolRel   ,err_div_tolAbs)
 BCTol      = max(maxval(abs(P_av))         *err_stress_tolRel,err_stress_tolAbs)

 if ((totalIter >= itmin .and. &
                           all([ err_div /divTol, &
                                 err_curl/curlTol, &
                                 err_BC  /BCTol       ] < 1.0_pReal)) &
             .or.    terminallyIll) then
   reason = 1
 elseif (totalIter >= itmax) then
   reason = -1
 else
   reason = 0
 endif

!--------------------------------------------------------------------------------------------------
! report
 write(6,'(1/,a)') ' ... reporting .............................................................'
 write(6,'(1/,a,f12.2,a,es8.2,a,es9.2,a)') ' error divergence = ', &
           err_div/divTol,  ' (',err_div, ' / m, tol = ',divTol,')'
 write(6,  '(a,f12.2,a,es8.2,a,es9.2,a)') ' error curl       = ', &
           err_curl/curlTol,' (',err_curl,' -,   tol = ',curlTol,')'
 write(6,  '(a,f12.2,a,es8.2,a,es9.2,a)') ' error BC         = ', &
           err_BC/BCTol,    ' (',err_BC,  ' Pa,  tol = ',BCTol,')' 
 write(6,'(/,a)') ' ==========================================================================='
 flush(6) 

end subroutine converged


!--------------------------------------------------------------------------------------------------
!> @brief forms the polarisation residual vector
!--------------------------------------------------------------------------------------------------
subroutine formResidual(in, FandF_tau, &
                        residuum, dummy,ierr)
 use numerics, only: &
   itmax, &
   itmin, &
   polarAlpha, &
   polarBeta
 use mesh, only: &
   grid, &
   grid3
 use math, only: &
   math_rotate_forward33, &
   math_rotate_backward33, &
   math_mul3333xx33, &
   math_invSym3333, &
   math_mul33x33
 use debug, only: &
   debug_level, &
   debug_spectral, &
   debug_spectralRotation
 use spectral_utilities, only: &
   wgt, &
   tensorField_real, &
   utilities_FFTtensorForward, &
   utilities_fourierGammaConvolution, &
   utilities_FFTtensorBackward, &
   utilities_constitutiveResponse, &
   utilities_divergenceRMS, &
   utilities_curlRMS
 use IO, only: &
    IO_intOut 
 use homogenization, only: &
   materialpoint_dPdF
 use FEsolving, only: &
   terminallyIll

 implicit none
 DMDALocalInfo, dimension(DMDA_LOCAL_INFO_SIZE) :: in                                               !< DMDA info (needs to be named "in" for macros like XRANGE to work)
 PetscScalar, dimension(3,3,2,XG_RANGE,YG_RANGE,ZG_RANGE), &
   target, intent(in) :: FandF_tau
 PetscScalar, dimension(3,3,2,X_RANGE,Y_RANGE,Z_RANGE),&
   target,  intent(out) :: residuum                                                                 !< residuum field
 PetscScalar, pointer, dimension(:,:,:,:,:) :: &
   F, &
   F_tau, &
   residual_F, &
   residual_F_tau
 PetscInt :: &
   PETScIter, &
   nfuncs
 PetscObject :: dummy
 PetscErrorCode :: ierr
 integer :: &
   i, j, k, e

 F                 => FandF_tau(1:3,1:3,1,&
                                XG_RANGE,YG_RANGE,ZG_RANGE)
 F_tau             => FandF_tau(1:3,1:3,2,&
                                XG_RANGE,YG_RANGE,ZG_RANGE)
 residual_F        => residuum(1:3,1:3,1,&
                                X_RANGE, Y_RANGE, Z_RANGE)
 residual_F_tau    => residuum(1:3,1:3,2,&
                                X_RANGE, Y_RANGE, Z_RANGE)

 F_av = sum(sum(sum(F,dim=5),dim=4),dim=3) * wgt
 call MPI_Allreduce(MPI_IN_PLACE,F_av,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 
 call SNESGetNumberFunctionEvals(snes,nfuncs,ierr); CHKERRQ(ierr)
 call SNESGetIterationNumber(snes,PETScIter,ierr);  CHKERRQ(ierr)

 if (nfuncs == 0 .and. PETScIter == 0) totalIter = -1                                                ! new increment
!--------------------------------------------------------------------------------------------------
! begin of new iteration
 newIteration: if (totalIter <= PETScIter) then
   totalIter = totalIter + 1
   write(6,'(1x,a,3(a,'//IO_intOut(itmax)//'))') &
           trim(incInfo), ' @ Iteration ', itmin, '≤',totalIter, '≤', itmax
   if (iand(debug_level(debug_spectral),debug_spectralRotation) /= 0) &
     write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') &
             ' deformation gradient aim (lab) =', transpose(math_rotate_backward33(F_aim,params%rotation_BC))
   write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') &
             ' deformation gradient aim       =', transpose(F_aim)
   flush(6)
 endif newIteration

!--------------------------------------------------------------------------------------------------
! 
 tensorField_real = 0.0_pReal
 do k = 1, grid3; do j = 1, grid(2); do i = 1, grid(1)
   tensorField_real(1:3,1:3,i,j,k) = &
     polarBeta*math_mul3333xx33(C_scale,F(1:3,1:3,i,j,k) - math_I3) -&
     polarAlpha*math_mul33x33(F(1:3,1:3,i,j,k), &
                        math_mul3333xx33(C_scale,F_tau(1:3,1:3,i,j,k) - F(1:3,1:3,i,j,k) - math_I3))
 enddo; enddo; enddo
 
!--------------------------------------------------------------------------------------------------
! doing convolution in Fourier space 
 call utilities_FFTtensorForward
 call utilities_fourierGammaConvolution(math_rotate_backward33(polarBeta*F_aim,params%rotation_BC)) 
 call utilities_FFTtensorBackward

!--------------------------------------------------------------------------------------------------
! constructing residual                         
 residual_F_tau = polarBeta*F - tensorField_real(1:3,1:3,1:grid(1),1:grid(2),1:grid3)

!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
 call utilities_constitutiveResponse(residual_F, &                                                    ! "residuum" gets field of first PK stress (to save memory)
                                     P_av,C_volAvg,C_minMaxAvg, &
                                     F - residual_F_tau/polarBeta,params%timeinc,params%rotation_BC)
 call MPI_Allreduce(MPI_IN_PLACE,terminallyIll,1,MPI_LOGICAL,MPI_LOR,PETSC_COMM_WORLD,ierr) 
 
!--------------------------------------------------------------------------------------------------
! stress BC handling
 F_aim = F_aim - math_mul3333xx33(S, ((P_av - params%stress_BC)))                                   ! S = 0.0 for no bc
 err_BC = maxval(abs((1.0_pReal-params%stress_mask) * math_mul3333xx33(C_scale,F_aim &
                                                -math_rotate_forward33(F_av,params%rotation_BC)) + &
                                params%stress_mask  * (P_av-params%stress_BC)))                     ! mask = 0.0 for no bc
! calculate divergence
 tensorField_real = 0.0_pReal
 tensorField_real(1:3,1:3,1:grid(1),1:grid(2),1:grid3) = residual_F                                   !< stress field in disguise
 call utilities_FFTtensorForward
 err_div = Utilities_divergenceRMS()                                                                  !< root mean squared error in divergence of stress
!--------------------------------------------------------------------------------------------------
! constructing residual
 e = 0
 do k = 1, grid3; do j = 1, grid(2); do i = 1, grid(1)
   e = e + 1
   residual_F(1:3,1:3,i,j,k) = &
     math_mul3333xx33(math_invSym3333(materialpoint_dPdF(1:3,1:3,1:3,1:3,1,e) + C_scale), &
                      residual_F(1:3,1:3,i,j,k) - math_mul33x33(F(1:3,1:3,i,j,k), &
                      math_mul3333xx33(C_scale,F_tau(1:3,1:3,i,j,k) - F(1:3,1:3,i,j,k) - math_I3))) &
                      + residual_F_tau(1:3,1:3,i,j,k)
 enddo; enddo; enddo
 
!--------------------------------------------------------------------------------------------------
! calculating curl
 tensorField_real = 0.0_pReal
 tensorField_real(1:3,1:3,1:grid(1),1:grid(2),1:grid3) = F
 call utilities_FFTtensorForward
 err_curl = Utilities_curlRMS()

end subroutine formResidual

end module grid_mech_spectral_polarisation
