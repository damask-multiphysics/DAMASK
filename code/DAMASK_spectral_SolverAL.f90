module DAMASK_spectral_SolverAL
 
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 
 use DAMASK_spectral_Utilities
 
 use math
   
 use mesh,  only : &
   mesh_spectral_getResolution, &
   mesh_spectral_getDimension
 
 implicit none
#include <finclude/petsc.h>
#include <finclude/petscvec.h90>
 
 character (len=*), parameter, public :: &
   DAMASK_spectral_SolverAL_label = 'AL'

!--------------------------------------------------------------------------------------------------
! PETSc data
 SNES     snes
 KSP      ksp
 DM       da
 Vec      x,r
 PetscErrorCode  ierr_psc
 PetscMPIInt rank
 PetscObject dummy
 PetscInt xs,xm,gxs,gxm
 PetscInt ys,ym,gys,gym
 PetscInt zs,zm,gzs,gzm
 character(len=1024) :: PetSc_options = '-snes_type ngmres -snes_ngmres_anderson -snes_monitor -snes_view'
 
 external FormFunctionLocal, SNESConverged_Interactive

!--------------------------------------------------------------------------------------------------
! common pointwise data
 real(pReal),    dimension(:,:,:,:,:), allocatable ::  F, F_lastInc, F_lambda, F_lambda_lastInc, P
 real(pReal),    dimension(:,:,:,:),   allocatable ::  coordinates
 real(pReal),    dimension(:,:,:),     allocatable ::  temperature
 
!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
 real(pReal), dimension(3,3) :: &
   F_aim = math_I3, &
   F_aim_lastInc = math_I3, &
   P_av
   
 real(pReal), dimension(3,3,3,3) :: &
   C_ref = 0.0_pReal, &
   C = 0.0_pReal
 
 integer(pInt) :: iter
 real(pReal)   :: err_div, err_stress 
 
 contains
 
 subroutine AL_init()
   
   use IO, only: &
     IO_read_JobBinaryFile, &
     IO_write_JobBinaryFile
 
   use FEsolving, only: &
     restartInc

   use DAMASK_interface, only: &
     getSolverJobName
     
   implicit none
   integer(pInt) :: i,j,k
   
   call Utilities_init()
   
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
   
   call Utilities_updateGamma(C_ref)
   
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
   call DMDASetLocalFunction(da,FormFunctionLocal,ierr_psc)
  
   call SNESSetDM(snes,da,ierr_psc)
   call SNESSetFunction(snes,r,SNESDMDAComputeFunction,da,ierr_psc)
   call SNESSetConvergenceTest(snes,SNESConverged_Interactive,dummy,PETSC_NULL_FUNCTION,ierr_psc)
   call PetscOptionsInsertString(PetSc_options,ierr_psc)
   call SNESSetFromOptions(snes,ierr_psc)  
   call DMDAGetCorners(da,xs,ys,zs,xm,ym,zm,ierr_psc)
   call DMDAGetCorners(da,gxs,gys,gzs,gxm,gym,gzm,ierr_psc)
  
   xs = xs+1; gxs = gxs+1; xm = xm-1; gxm = gxm-1 
   ys = ys+1; gys = gys+1; ym = ym-1; gym = gym-1
   zs = zs+1; gzs = gzs+1; zm = zm-1; gzm = gzm-1
 
 end subroutine AL_init
 
type(solutionState) function AL_solution(guessmode,timeinc,timeinc_old,P_BC,F_BC,mask_stressVector,velgrad,rotation_BC)
 
 use numerics, only: &
   itmax, &
   itmin, &
   update_gamma

 use IO, only: &
   IO_write_JobBinaryFile
   
 use FEsolving, only: &
   restartWrite
 
 implicit none

!--------------------------------------------------------------------------------------------------
! input data for solution

 real(pReal), intent(in) :: timeinc, timeinc_old
 real(pReal), intent(in) :: guessmode
 logical,     intent(in) :: velgrad
 real(pReal), dimension(3,3), intent(in) :: P_BC,F_BC,rotation_BC
 logical,     dimension(9),   intent(in) :: mask_stressVector
 
!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.

 real(pReal), dimension(3,3), parameter :: ones = 1.0_pReal, zeroes = 0.0_pReal    
 real(pReal), dimension(3,3)            :: temp33_Real  
 real(pReal), dimension(3,3,3,3)        :: S
 real(pReal), dimension(3,3)            :: mask_stress, &
                                           mask_defgrad, &
                                           deltaF_aim, &
                                           F_aim_lab, &
                                           F_aim_lab_lastIter
 integer(pInt) :: i, j, k
 logical       :: ForwardData
 real(pReal)   :: defgradDet
 real(pReal)   :: defgradDetMax, defgradDetMin
 
 PetscScalar, pointer :: xx_psc(:)

 mask_stress = merge(ones,zeroes,reshape(mask_stressVector,[3,3]))                                   
 mask_defgrad  = merge(zeroes,ones,reshape(mask_stressVector,[3,3]))
 
 if (restartWrite) then
   write(6,'(a)') 'writing converged results for restart'
   call IO_write_jobBinaryFile(777,'convergedSpectralDefgrad',size(F_lastInc))                        ! writing deformation gradient field to file
   write (777,rec=1) F_LastInc
   close (777)
   call IO_write_jobBinaryFile(777,'C',size(C))
   write (777,rec=1) C
   close(777)
 endif 

 ForwardData = .True.
 if (velgrad) then                                                        ! calculate deltaF_aim from given L and current F
   deltaF_aim = timeinc * mask_defgrad * math_mul33x33(F_BC, F_aim)
 else                                                                                         ! deltaF_aim = fDot *timeinc where applicable
   deltaF_aim = timeinc * mask_defgrad * F_BC
 endif

!--------------------------------------------------------------------------------------------------
! winding forward of deformation aim in loadcase system
 temp33_Real = F_aim                                            
 F_aim = F_aim &                                                                         
         + guessmode * mask_stress * (F_aim - F_aim_lastInc)*timeinc/timeinc_old &      
         + deltaF_aim
 F_aim_lastInc = temp33_Real
 F_aim_lab = math_rotate_backward33(F_aim,rotation_BC)                            ! boundary conditions from load frame into lab (Fourier) frame

!--------------------------------------------------------------------------------------------------
! update local deformation gradient and coordinates
 deltaF_aim = math_rotate_backward33(deltaF_aim,rotation_BC)
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   temp33_Real = F(i,j,k,1:3,1:3)
   F(i,j,k,1:3,1:3) = F(i,j,k,1:3,1:3) &                                                             ! decide if guessing along former trajectory or apply homogeneous addon
                   + guessmode * (F(i,j,k,1:3,1:3) - F_lastInc(i,j,k,1:3,1:3))* &
                   timeinc/timeinc_old + (1.0_pReal-guessmode) * deltaF_aim                                ! if not guessing, use prescribed average deformation where applicable
   F_lastInc(i,j,k,1:3,1:3) = temp33_Real 
   temp33_Real = F_lambda(i,j,k,1:3,1:3)
   F_lambda(i,j,k,1:3,1:3) = F_lambda(i,j,k,1:3,1:3) &                                                             ! decide if guessing along former trajectory or apply homogeneous addon
                   + guessmode * (F_lambda(i,j,k,1:3,1:3) - F_lambda_lastInc(i,j,k,1:3,1:3))* &
                   timeinc/timeinc_old + (1.0_pReal-guessmode) * deltaF_aim                                ! if not guessing, use prescribed average deformation where applicable
   F_lambda_lastInc(i,j,k,1:3,1:3) = temp33_Real 
 enddo; enddo; enddo
 call deformed_fft(res,geomdim,math_rotate_backward33(F_aim,rotation_BC),&          ! calculate current coordinates
                                                          1.0_pReal,F_lastInc,coordinates)

 iter = 0_pInt
 S = Utilities_stressBC(rotation_BC,mask_stressVector,C)
 if (update_gamma) call Utilities_updateGamma(C)
 
 call VecGetArrayF90(x,xx_psc,ierr_psc)
 call FormInitialGuessLocal(xx_psc)
 call VecRestoreArrayF90(x,xx_psc,ierr_psc)    
 call SNESSolve(snes,PETSC_NULL_OBJECT,x,ierr_psc)
   
 convergenceLoop: do while((iter < itmax .and. (any([err_div ,err_stress] > 1.0_pReal)))&
          .or. iter < itmin)
   
   iter = iter + 1_pInt
!--------------------------------------------------------------------------------------------------
! report begin of new iteration
   write(6,'(a)') ''
   write(6,'(a)') '=================================================================='
   write(6,'(3(a,i6.6))') ' @ Iter. ',itmin,' < ',iter,' < ',itmax
   write(6,'(a,/,3(3(f12.7,1x)/))',advance='no') 'deformation gradient aim =',&
                                                             math_transpose33(F_aim)
   F_aim_lab_lastIter = math_rotate_backward33(F_aim,rotation_BC)

!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
   call constitutiveResponse(coordinates,F,F_lastInc,temperature,timeinc,&
                                P,C,P_av,ForwardData,rotation_BC)
   ForwardData = .False.
   
!--------------------------------------------------------------------------------------------------
! stress BC handling
   if(any(mask_stressVector)) then                                                             ! calculate stress BC if applied
     F_aim = F_aim - math_mul3333xx33(S, ((P_av - P_BC)))
     err_stress = mask_stress * (P_av - P_BC)))
   else
     err_stress = 0.0_pReal
   endif
                                
  F_aim_lab = math_rotate_backward33(F_aim,rotation_BC)                           ! boundary conditions from load frame into lab (Fourier) frame
 
!--------------------------------------------------------------------------------------------------
! updated deformation gradient
  field_real(1:res(1),1:res(2),1:res(3),1:3,1:3) = P
  call FFT_forward()
  err_div = calcDivergence()
  call convolution_fourier(F_aim_lab_lastIter - F_aim_lab, C_ref) 
  call FFT_backward()

  do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
    F(i,j,k,1:3,1:3) = F(i,j,k,1:3,1:3) - field_real(i,j,k,1:3,1:3)                       ! F(x)^(n+1) = F(x)^(n) + correction;  *wgt: correcting for missing normalization
  enddo; enddo; enddo
  
!--------------------------------------------------------------------------------------------------
! calculate some additional output
  if(debugGeneral) then
    maxCorrectionSkew = 0.0_pReal
    maxCorrectionSym  = 0.0_pReal
    temp33_Real = 0.0_pReal
    do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
      maxCorrectionSym  = max(maxCorrectionSym,&
                              maxval(math_symmetric33(field_real(i,j,k,1:3,1:3))))
      maxCorrectionSkew = max(maxCorrectionSkew,&
                              maxval(math_skew33(field_real(i,j,k,1:3,1:3))))
      temp33_Real = temp33_Real + field_real(i,j,k,1:3,1:3)
    enddo; enddo; enddo
    write(6,'(a,1x,es11.4)') 'max symmetric correction of deformation =',&
                                  maxCorrectionSym*wgt
    write(6,'(a,1x,es11.4)') 'max skew      correction of deformation =',&
                                  maxCorrectionSkew*wgt
    write(6,'(a,1x,es11.4)') 'max sym/skew of avg correction =         ',&
                                  maxval(math_symmetric33(temp33_real))/&
                                  maxval(math_skew33(temp33_real))
  endif

!--------------------------------------------------------------------------------------------------
! calculate bounds of det(F) and report
  if(debugGeneral) then
    defgradDetMax = -huge(1.0_pReal)
    defgradDetMin = +huge(1.0_pReal)
    do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
      defgradDet = math_det33(F(i,j,k,1:3,1:3))
      defgradDetMax = max(defgradDetMax,defgradDet)
      defgradDetMin = min(defgradDetMin,defgradDet) 
    enddo; enddo; enddo

    write(6,'(a,1x,es11.4)') 'max determinant of deformation =', defgradDetMax
    write(6,'(a,1x,es11.4)') 'min determinant of deformation =', defgradDetMin
  endif
 enddo convergenceLoop

end function AL_solution

subroutine AL_destroy()

implicit none

call VecDestroy(x,ierr_psc)
call VecDestroy(r,ierr_psc)
call SNESDestroy(snes,ierr_psc)
call DMDestroy(da,ierr_psc)
call PetscFinalize(ierr_psc)
call Utilities_destroy()

end subroutine AL_destroy

! ------------------------------------------------------------------- 

subroutine FormInitialGuessLocal(xx_psc)

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
end subroutine FormInitialGuessLocal

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
subroutine FormFunctionLocal(in,x_scal,f_scal,dummy,ierr_psc)
  
  use numerics, only: &
   itmax, &
   itmin
  
  implicit none
#include <finclude/petsc.h>

!  Input/output variables:
  DMDALocalInfo in(DMDA_LOCAL_INFO_SIZE)
  PetscScalar x_scal(0:17,XG_RANGE,YG_RANGE,ZG_RANGE)  
  PetscScalar f_scal(0:17,X_RANGE,Y_RANGE,Z_RANGE)
  real(pReal), dimension (3,3) ::  temp
  PetscObject dummy
  
! Compute function over the locally owned part of the grid

  iter = iter + 1_pInt
  
!--------------------------------------------------------------------------------------------------
! report begin of new iteration
  write(6,'(a)') ''
  write(6,'(a)') '=================================================================='
  write(6,'(3(a,i6.6))') ' @ Iter. ',itmin,' < ',iter,' < ',itmax
  write(6,'(a,/,3(3(f12.7,1x)/))',advance='no') 'deformation gradient aim =',&
                                                            math_transpose33(F_aim)
  
  F_star_av = 0.0
  lambda_av = 0.0
  do k=gzs,gze; do j=gys,gye; do i=gxs,gxe
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
    F_star_av = F_star_av + F(i,j,k,1:3,1:3)
    lambda_av = lambda_av + F_lambda(i,j,k,1:3,1:3)
  enddo; enddo; enddo
  F_star_av = F_star_av *wgt
  lambda_av = math_mul3333xx33(C_inc0,lambda_av*wgt-math_I3)
  
!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
   call constitutiveResponse(coordinates,F,F_lastInc,temperature,timeinc,&
                                P,C,P_av,ForwardData,rotation_BC)
   ForwardData = .False.
   
!--------------------------------------------------------------------------------------------------
! stress BC handling
   if(any(mask_stressVector)) then                                                             ! calculate stress BC if applied
     F_aim = F_aim - math_mul3333xx33(S, ((P_av - P_BC)))
     err_stress = mask_stress * (P_av - P_BC)))
   else
     err_stress = 0.0_pReal
   endif
                                
  F_aim_lab = math_rotate_backward33(F_aim,rotation_BC)  
  
!--------------------------------------------------------------------------------------------------
! doing Fourier transform
  field_real = 0.0_pReal
  do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
    field_real(i,j,k,1:3,1:3) = math_mul3333xx33(C_ref,F_lambda(i,j,k,1:3,1:3)-F(i,j,k,1:3,1:3))
    
  enddo; enddo; enddo
  
  call Utilities_forwardFFT()
  call Utilities_fourierConvolution(F_aim_lab) 
  call Utilities_backwardFFT()

  err_f = 0.0_pReal
  err_f_point = 0.0_pReal
  err_p = 0.0_pReal
  err_p_point = 0.0_pReal
  
  do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
    temp33_real = field_real(i,j,k,1:3,1:3) - F(i,j,k,1:3,1:3)
    err_f_point = max(err_f_point, maxval(abs(temp33_real)))
    err_f = err_f + sum(temp33_real*temp33_real)
    
    temp33_real = F_lambda(i,j,k,1:3,1:3) - &
                  math_mul3333xx33(S_inc0,P(i,j,k,1:3,1:3)) + math_I3
    err_p_point = max(err_p_point, maxval(abs(temp33_real)))
    err_p = err_p + sum(temp33_real*temp33_real)
  enddo; enddo; enddo
                                       
  err_f = wgt*sqrt(err_f/sum((F_aim-math_I3)*(F_aim-math_I3)))
  err_p = wgt*sqrt(err_p/sum((F_aim-math_I3)*(F_aim-math_I3)))

  write(6,'(a,es14.7,es14.7)') 'error stress     = ',err_stress/err_stress_tol
  write(6,*) '  ' 
  write(6,'(a,es14.7)') 'max abs err F', err_f
  write(6,'(a,es14.7)') 'max abs err P', err_p
  
  do k=zs,ze; do j=ys,ye; do i=xs,xe
    temp = math_mul3333xx33(S_inc0,P(i,j,k,1:3,1:3)) + math_I3 - F_lambda(i,j,k,1:3,1:3) &
           + F(i,j,k,1:3,1:3) - field_real(i,j,k,1:3,1:3)
    f_scal(0,i,j,k) = temp(1,1)
    f_scal(1,i,j,k) = temp(1,2)
    f_scal(2,i,j,k) = temp(1,3)
    f_scal(3,i,j,k) = temp(2,1)
    f_scal(4,i,j,k) = temp(2,2)
    f_scal(5,i,j,k) = temp(2,3)
    f_scal(6,i,j,k) = temp(3,1)
    f_scal(7,i,j,k) = temp(3,2)
    f_scal(8,i,j,k) = temp(3,3)
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
end subroutine FormFunctionLocal

! ---------------------------------------------------------------------
! User defined convergence check
!
subroutine SNESConverged_Interactive(snes,it,xnorm,snorm,fnorm,reason,dummy,ierr_psc)
  
  implicit none
#include <finclude/petsc.h>

! Input/output variables:
  SNES snes
  PetscInt it
  PetscReal xnorm, snorm, fnorm
  SNESConvergedReason reason
  PetscObject dummy
  PetscErrorCode ierr_psc

  err_crit = max(err_stress/err_stress_tol, &
             err_f/1e-6, err_p/1e-5)
             !fnorm*wgt/sqrt(sum((F_star_av-math_I3)*(F_star_av-math_I3)))/err_div_tol)
             
  if ((err_crit > 1.0_pReal .or. it < itmin) .and. it < itmax) then
    reason = 0
  else  
    reason = 1
  endif    

  return
end subroutine SNESConverged_Interactive

end module DAMASK_spectral_SolverAL
