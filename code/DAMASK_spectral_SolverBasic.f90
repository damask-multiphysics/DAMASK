!--------------------------------------------------------------------------------------------------
!* $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Basic scheme solver
!--------------------------------------------------------------------------------------------------
module DAMASK_spectral_SolverBasic
 
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 
 use DAMASK_spectral_Utilities
 
 use math, only: &
   math_I3, &
   math_mul33x33,
   
 use mesh,  only : &
   mesh_spectral_getResolution, &
   mesh_spectral_getDimension, &
   math_rotate_backward33, &
   math_transpose33,&
   math_mul3333xx33, &
   math_eigenvalues33
 
 implicit none
 
 character (len=*), parameter, public :: &
   DAMASK_spectral_SolverBasic_label = 'basic'

!--------------------------------------------------------------------------------------------------
! common pointwise data
 real(pReal),    dimension(:,:,:,:,:), allocatable ::  F, F_lastInc, P
 real(pReal),    dimension(:,:,:,:),   allocatable ::  coordinates
 real(pReal),    dimension(:,:,:),     allocatable ::  temperature
 
!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
 real(pReal), dimension(3,3) :: &
   F_aim = math_I3, &
   F_aim_lastInc = math_I3
 real(pReal), dimension(3,3,3,3) :: &
   C = 0.0_pReal
 
 
 contains
 
 subroutine basic_init()
   
   use IO, only: &
     IO_read_JobBinaryFile, &
     IO_write_JobBinaryFile
 
   use FEsolving, only: &
     restartInc

   use DAMASK_interface, only: &
     getSolverJobName
     
   implicit none
   integer(pInt) :: i,j,k

   real(pReal), dimension(3,3) :: temp33_Real
 
   write(6,'(a)') ''
   write(6,'(a)') ' <<<+-  DAMASK_spectral_solverBasic init  -+>>>'
   write(6,'(a)') ' $Id$'
#include "compilation_info.f90"
   write(6,'(a)') ''
     call Utilities_Init()
   
   allocate (F          (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
   allocate (F_lastInc  (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
   allocate (P          (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
   allocate (coordinates(  res(1),  res(2),res(3),3),    source = 0.0_pReal)
   allocate (temperature(  res(1),  res(2),res(3)),      source = 0.0_pReal)
   
!--------------------------------------------------------------------------------------------------
! init fields                 
   if (restartInc == 1_pInt) then                                                                     ! no deformation (no restart)
     do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
       F(i,j,k,1:3,1:3) = math_I3
       F_lastInc(i,j,k,1:3,1:3) = math_I3
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
     call IO_read_jobBinaryFile(777,'F_aim',trim(getSolverJobName()),size(F_aim))
     read (777,rec=1) F_aim
     close (777)
     call IO_read_jobBinaryFile(777,'F_aim_lastInc',trim(getSolverJobName()),size(F_aim_lastInc))
     read (777,rec=1) F_aim_lastInc
     close (777)
  
     coordinates = 0.0 ! change it later!!!
   endif
   
   call Utilities_constitutiveResponse(coordinates,F,F_lastInc,temperature,0.0_pReal,&
                                P,C,temp33_Real,.false.,math_I3)
   
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
 
 end subroutine basic_init
 
type(solutionState) function basic_solution(guessmode,timeinc,timeinc_old,P_BC,F_BC,mask_stressVector,velgrad,rotation_BC)
 
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
                                           F_aim_lab_lastIter, &
                                           P_av
 real(pReal)   :: err_div, err_stress       
 integer(pInt) :: iter
 integer(pInt) :: i, j, k
 logical       :: ForwardData
 real(pReal)   :: defgradDet
 real(pReal)   :: defgradDetMax, defgradDetMin

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
                   + guessmode * (F(i,j,k,1:3,1:3) - F_lastInc(i,j,k,1:3,1:3))*timeinc/timeinc_old&  ! guessing... 
                   + (1.0_pReal-guessmode) * deltaF_aim                                ! if not guessing, use prescribed average deformation where applicable
   F_lastInc(i,j,k,1:3,1:3) = temp33_Real 
 enddo; enddo; enddo
 call deformed_fft(res,geomdim,math_rotate_backward33(F_aim,rotation_BC),&          ! calculate current coordinates
                                                          1.0_pReal,F_lastInc,coordinates)

 iter = 0_pInt
 S = Utilities_maskedCompliance(rotation_BC,mask_stressVector,C)
 if (update_gamma) call Utilities_updateGamma(C)
   
 convergenceLoop: do 
   
   iter = iter + 1_pInt
!--------------------------------------------------------------------------------------------------
! report begin of new iteration
   write(6,'(a)') ''
   write(6,'(a)') '=================================================================='
   write(6,'(3(a,i6.6))') ' Iter. ',itmin,' < ',iter,' < ',itmax + 1_pInt
   write(6,'(a,/,3(3(f12.7,1x)/))',advance='no') 'deformation gradient aim =',&
                                                             math_transpose33(F_aim)
   F_aim_lab_lastIter = math_rotate_backward33(F_aim,rotation_BC)

!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
print*, 'FLast 111', F_lastInc(1,1,1,1:3,1:3)
   call Utilities_constitutiveResponse(coordinates,F_lastInc,F,temperature,timeinc,&
                                P,C,P_av,ForwardData,rotation_BC)
   ForwardData = .False.
   
!--------------------------------------------------------------------------------------------------
! stress BC handling
   F_aim = F_aim - math_mul3333xx33(S, ((P_av - P_BC))) !S = 0.0 for no bc
   err_stress = maxval(mask_stress * (P_av - P_BC))     ! mask = 0.0 for no bc

                                
  F_aim_lab = math_rotate_backward33(F_aim,rotation_BC)                            ! boundary conditions from load frame into lab (Fourier) frame
 
!--------------------------------------------------------------------------------------------------
! updated deformation gradient
  field_real = 0.0_pReal
  field_real(1:res(1),1:res(2),1:res(3),1:3,1:3) = P
  call Utilities_forwardFFT()
  err_div = Utilities_divergenceRMS()
  call Utilities_fourierConvolution(F_aim_lab_lastIter - F_aim_lab) 
  call Utilities_backwardFFT()

  temp33_real =0.0_pReal

  do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
    F(i,j,k,1:3,1:3) = F(i,j,k,1:3,1:3) - field_real(i,j,k,1:3,1:3)                       ! F(x)^(n+1) = F(x)^(n) + correction;  *wgt: correcting for missing normalization
    temp33_real = temp33_real + field_real(i,j,k,1:3,1:3)
  enddo; enddo; enddo
  if (Convergenced(err_div,P_av,err_stress,P_av,iter)) exit
  
 enddo convergenceLoop

end function basic_solution

logical function Convergenced(err_div,P_av,err_stress,P_av2,iter)

 use numerics, only: &
  itmax, &
  itmin, &
  err_div_tol, &
  err_stress_tolrel, &
  err_stress_tolabs
    
 implicit none
  
 real(pReal), dimension(3,3) :: P_av, P_av2
 real(pReal) :: err_div, err_stress, field_av_L2
 integer(pInt) :: iter
  
 field_av_L2 = sqrt(maxval(math_eigenvalues33(math_mul33x33(P_av,&                    ! L_2 norm of average stress (http://mathworld.wolfram.com/SpectralNorm.html)
                                              math_transpose33(P_av)))))
 Convergenced = (iter < itmax) .and. (iter > itmin) .and. &
                 all([err_div/field_av_L2/err_div_tol,&
                       err_stress/min(maxval(abs(P_av2))*err_stress_tolrel,err_stress_tolabs)] < 1.0_pReal)
  
  
  write(6,'(a,f6.2,a,es11.4,a)') 'error stress = ', err_stress/min(maxval(abs(P_av2))*err_stress_tolrel,err_stress_tolabs), &
                                                                            ' (',err_stress,' Pa)'  
  write(6,'(a,f6.2,a,es11.4,a)') 'error divergence = ', err_div/field_av_L2/err_div_tol,&
                                                           ' (',err_div,' N/m³)'
end function Convergenced

subroutine basic_destroy()

implicit none
call Utilities_destroy()

end subroutine basic_destroy

end module DAMASK_spectral_SolverBasic


!--------------------------------------------------------------------------------------------------
! calculate some additional output
 ! if(debugGeneral) then
 !   maxCorrectionSkew = 0.0_pReal
 !   maxCorrectionSym  = 0.0_pReal
 !   temp33_Real = 0.0_pReal
 !   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
 !     maxCorrectionSym  = max(maxCorrectionSym,&
 !                             maxval(math_symmetric33(field_real(i,j,k,1:3,1:3))))
 !     maxCorrectionSkew = max(maxCorrectionSkew,&
 !                             maxval(math_skew33(field_real(i,j,k,1:3,1:3))))
 !     temp33_Real = temp33_Real + field_real(i,j,k,1:3,1:3)
 !   enddo; enddo; enddo
 !   write(6,'(a,1x,es11.4)') 'max symmetric correction of deformation =',&
 !                                 maxCorrectionSym*wgt
 !   write(6,'(a,1x,es11.4)') 'max skew      correction of deformation =',&
 !                                 maxCorrectionSkew*wgt
 !   write(6,'(a,1x,es11.4)') 'max sym/skew of avg correction =         ',&
 !                                 maxval(math_symmetric33(temp33_real))/&
 !                                 maxval(math_skew33(temp33_real))
 ! endif

!--------------------------------------------------------------------------------------------------
! calculate bounds of det(F) and report
  ! if(debugGeneral) then
    ! defgradDetMax = -huge(1.0_pReal)
    ! defgradDetMin = +huge(1.0_pReal)
    ! do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
      ! defgradDet = math_det33(F(i,j,k,1:3,1:3))
      ! defgradDetMax = max(defgradDetMax,defgradDet)
      ! defgradDetMin = min(defgradDetMin,defgradDet) 
    ! enddo; enddo; enddo

    ! write(6,'(a,1x,es11.4)') 'max determinant of deformation =', defgradDetMax
    ! write(6,'(a,1x,es11.4)') 'min determinant of deformation =', defgradDetMin
  ! endif
