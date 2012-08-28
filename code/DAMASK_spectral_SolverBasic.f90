!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Basic scheme solver
!> @details this solver follows closely the original large strain formulation presented by
!> Suquet. The iterative procedure is solved using a fix-point iteration
!--------------------------------------------------------------------------------------------------
module DAMASK_spectral_SolverBasic
 use prec, only: & 
   pInt, &
   pReal
 
 use math, only: &
   math_I3
 
 use DAMASK_spectral_Utilities, only: &
   solutionState
 
 implicit none
 character (len=*), parameter, public :: &
   DAMASK_spectral_SolverBasic_label = 'basic'

!--------------------------------------------------------------------------------------------------
! pointwise data
 real(pReal),  private,  dimension(:,:,:,:,:), allocatable ::  F, F_lastInc, P
 real(pReal),  private,  dimension(:,:,:,:),   allocatable ::  coordinates
 real(pReal),  private,  dimension(:,:,:),     allocatable ::  temperature
 
!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
 real(pReal), private, dimension(3,3) :: &
   F_aim = math_I3, &
   F_aim_lastInc = math_I3
 real(pReal), private,dimension(3,3,3,3) :: &
   C = 0.0_pReal
 
contains
 
!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
subroutine basic_init()
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
   debugrestart
      
 use mesh, only: &
   res, &
   geomdim

 implicit none

 integer(pInt) :: i,j,k
 real(pReal), dimension(3,3) :: temp33_Real
 
 
 call Utilities_Init()
 write(6,'(a)') ''
 write(6,'(a)') ' <<<+-  DAMASK_spectral_solverBasic init  -+>>>'
 write(6,'(a)') ' $Id$'
#include "compilation_info.f90"
 write(6,'(a)') ''
  

 allocate (F          (  3,3,res(1),  res(2),res(3)),  source = 0.0_pReal)
 allocate (F_lastInc  (  3,3,res(1),  res(2),res(3)),  source = 0.0_pReal)
 allocate (P          (  3,3,res(1),  res(2),res(3)),  source = 0.0_pReal)
 allocate (coordinates(  res(1),  res(2),res(3),3),    source = 0.0_pReal)
 allocate (temperature(  res(1),  res(2),res(3)),      source = 0.0_pReal)
   
!--------------------------------------------------------------------------------------------------
! init fields                 
 if (restartInc == 1_pInt) then                                                                     ! no deformation (no restart)
   F         = spread(spread(spread(math_I3,3,res(1)),4,res(2)),5,res(3))                           ! initialize to identity
   F_lastInc = F
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
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
 !no rotation bc call deformed_fft(res,geomdim,math_rotate_backward33(F_aim,rotation_BC),1.0_pReal,F_lastInc,coordinates) 
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


!--------------------------------------------------------------------------------------------------
!> @brief solution for the basic scheme with internal iterations
!--------------------------------------------------------------------------------------------------
type(solutionState) function basic_solution(guessmode,timeinc,timeinc_old,P_BC,F_BC,temperature_bc,rotation_BC)
 
 use numerics, only: &
   itmax, &
   itmin, &
   update_gamma
 use math, only: &
   math_mul33x33 ,&
   math_rotate_backward33, &
   math_transpose33, &
   math_mul3333xx33
 use mesh, only: &
   res,&
   geomdim, &
   deformed_fft
 use IO, only: &
   IO_write_JobBinaryFile
   
 use DAMASK_spectral_Utilities, only: &
   boundaryCondition, &
   field_real, &
   Utilities_forwardField, &
   Utilities_maskedCompliance, &
   Utilities_forwardFFT, &
   Utilities_divergenceRMS, &
   Utilities_fourierConvolution, &
   Utilities_backwardFFT, &
   Utilities_updateGamma, &
   Utilities_constitutiveResponse
     
 use FEsolving, only: &
   restartWrite
 
 implicit none
!--------------------------------------------------------------------------------------------------
! input data for solution
 real(pReal), intent(in) :: timeinc, timeinc_old, temperature_bc, guessmode
 type(boundaryCondition),      intent(in) :: P_BC,F_BC
 real(pReal), dimension(3,3), intent(in) :: rotation_BC
 

 
 real(pReal), dimension(3,3,3,3)        :: S
 real(pReal), dimension(3,3)            :: deltaF_aim, &
                                           F_aim_lab, &
                                           F_aim_lab_lastIter, &
                                           P_av
!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
 real(pReal)   :: err_div, err_stress       
 integer(pInt) :: iter, row, column, i, j, k
 logical       :: ForwardData
 real(pReal)   :: defgradDet, defgradDetMax, defgradDetMin
 real(pReal), dimension(3,3)            :: temp33_Real 

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
 call deformed_fft(res,geomdim,math_rotate_backward33(F_aim,rotation_BC),1.0_pReal,F_lastInc,coordinates)

!--------------------------------------------------------------------------------------------------
! update stiffness (and gamma operator)
 S = Utilities_maskedCompliance(rotation_BC,P_BC%maskLogical,C)
 if (update_gamma) call Utilities_updateGamma(C)
 
 iter = 0_pInt
 ForwardData = .True.
 convergenceLoop: do while(iter < itmax)
   
   iter = iter + 1_pInt
!--------------------------------------------------------------------------------------------------
! report begin of new iteration
   write(6,'(a)') ''
   write(6,'(a)') '=================================================================='
   write(6,'(3(a,i6.6))') ' Iter. ',itmin,' < ',iter,' < ',itmax + 1_pInt
   write(6,'(a,/,3(3(f12.7,1x)/))',advance='no') 'deformation gradient aim =', &
                                                                        math_transpose33(F_aim)
   F_aim_lab_lastIter = math_rotate_backward33(F_aim,rotation_BC)

!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
   call Utilities_constitutiveResponse(coordinates,F_lastInc,F,temperature,timeinc,&
                                 P,C,P_av,ForwardData,rotation_BC)
   ForwardData = .False.
   
!--------------------------------------------------------------------------------------------------
! stress BC handling
   F_aim = F_aim - math_mul3333xx33(S, ((P_av - P_BC%values))) !S = 0.0 for no bc
   err_stress = maxval(P_BC%maskFloat * (P_av - P_BC%values))     ! mask = 0.0 for no bc
   F_aim_lab = math_rotate_backward33(F_aim,rotation_BC)                            ! boundary conditions from load frame into lab (Fourier) frame
 
!--------------------------------------------------------------------------------------------------
! updated deformation gradient using fix point algorithm of basic scheme
   field_real = 0.0_pReal
   field_real(1:res(1),1:res(2),1:res(3),1:3,1:3) = reshape(P,[res(1),res(2),res(3),3,3],&
                                                               order=[4,5,1,2,3]) ! field real has a different order
   call Utilities_forwardFFT()
   err_div = Utilities_divergenceRMS()
   call Utilities_fourierConvolution(F_aim_lab_lastIter - F_aim_lab) 
   call Utilities_backwardFFT()
   F = F - reshape(field_real(1:res(1),1:res(2),1:res(3),1:3,1:3),shape(F),order=[3,4,5,1,2])                       ! F(x)^(n+1) = F(x)^(n) + correction;  *wgt: correcting for missing normalization
   basic_solution%converged = basic_Converged(err_div,P_av,err_stress,P_av)
   if (basic_solution%converged .and. iter > itmin) exit  
 enddo convergenceLoop

end function basic_solution


!--------------------------------------------------------------------------------------------------
!> @brief convergence check for basic scheme based on div of P and deviation from stress aim
!--------------------------------------------------------------------------------------------------
logical function basic_Converged(err_div,pAvgDiv,err_stress,pAvgStress)

 use numerics, only: &
   itmin, &
   err_div_tol, &
   err_stress_tolrel, &
   err_stress_tolabs
  
 use math, only: &
   math_mul33x33, &
   math_eigenvalues33, &
   math_transpose33
    
 implicit none
  
 real(pReal), dimension(3,3), intent(in) :: &
   pAvgDiv,&
   pAvgStress
 
 real(pReal), intent(in) :: &
   err_div, &
   err_stress
 
 real(pReal) :: &
   err_stress_tol, &
   pAvgDivL2
 

  
 pAvgDivL2 = sqrt(maxval(math_eigenvalues33(math_mul33x33(pAvgDiv,math_transpose33(pAvgDiv)))))                    ! L_2 norm of average stress (http://mathworld.wolfram.com/SpectralNorm.html)
 err_stress_tol = min(maxval(abs(pAvgStress))*err_stress_tolrel,err_stress_tolabs)
 
 basic_Converged = all([ err_div/pAvgDivL2/err_div_tol,&
                           err_stress/err_stress_tol    ]  < 1.0_pReal)
  
 write(6,'(a,f6.2,a,es11.4,a)') 'error divergence = ', err_div/pAvgDivL2/err_div_tol,&
                                                       ' (',err_div,' N/m³)'
 write(6,'(a,f6.2,a,es11.4,a)') 'error stress =     ', err_stress/err_stress_tol, &
                                                       ' (',err_stress,' Pa)'  

end function basic_Converged

subroutine basic_destroy()
 
 use DAMASK_spectral_Utilities, only: &
   Utilities_destroy
 
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
