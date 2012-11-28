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
   tSolutionState
 
 implicit none
 character (len=*), parameter, public :: &
   DAMASK_spectral_SolverBasic_label = 'basic'

!--------------------------------------------------------------------------------------------------
! pointwise global data
 real(pReal),  private,  dimension(:,:,:,:,:), allocatable ::  &
   F, &                                                                                             !< deformation gradient field
   F_lastInc, &                                                                                     !< deformation gradient field last increment
   Fdot                                                                                             !< assumed rate for F n to F n+1
 real(pReal),  private                                     ::  temperature                          !< not pointwise
 
!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
 real(pReal), private, dimension(3,3) :: &
   F_aim = math_I3, &                                                                               !< deformation gradient aim
   F_aim_lastInc = math_I3, &                                                                       !< deformation gradient aim last increment
   F_aimDot = 0.0_pReal                                                                             !< assumed rate
 real(pReal), private,dimension(3,3,3,3) :: &
   C = 0.0_pReal, &                                                                                 !< average stiffness
   C_lastInc = 0.0_pReal                                                                            !< average stiffness last increment
 
contains
 
!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
subroutine basic_init()
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 
 use IO, only: &
   IO_read_JobBinaryFile, &
   IO_write_JobBinaryFile, &
   IO_intOut
 use FEsolving, only: &
   restartInc
 use DAMASK_interface, only: &
   getSolverJobName
 use DAMASK_spectral_Utilities, only: &
   Utilities_init, &
   Utilities_constitutiveResponse, &
   Utilities_updateGamma, &
   debugRestart
 use mesh, only: &
   res, &
   geomdim, &
   mesh_ipCoordinates, &
   mesh_NcpElems, &
   mesh_deformedCoordsFFT

 implicit none
 real(pReal), dimension(3,3,res(1),res(2),res(3)) :: P
 integer(pInt) :: &
   i, j, k
 real(pReal), dimension(3,3) :: &
   temp33_Real
 real(pReal), dimension(3,3,3,3) :: &
   temp3333_Real
 

 call Utilities_Init()
 write(6,'(/,a)') ' <<<+-  DAMASK_spectral_solverBasic init  -+>>>'
 write(6,'(a)') ' $Id$'
#include "compilation_info.f90"
 write(6,'(a)') ''
  
!--------------------------------------------------------------------------------------------------
! allocate global fields
 allocate (F          (  3,3,res(1),  res(2),res(3)),  source = 0.0_pReal)
 allocate (F_lastInc  (  3,3,res(1),  res(2),res(3)),  source = 0.0_pReal)
 allocate (Fdot       (  3,3,res(1),  res(2),res(3)),  source = 0.0_pReal)
   
!--------------------------------------------------------------------------------------------------
! init fields and average quantities
 if (restartInc == 1_pInt) then                                                                     ! no deformation (no restart)
   F         = spread(spread(spread(math_I3,3,res(1)),4,res(2)),5,res(3))                           ! initialize to identity
   F_lastInc = F
 elseif (restartInc > 1_pInt) then                                                                  ! using old values from file                                                      
   if (debugRestart) write(6,'(a,'//IO_intOut(restartInc-1_pInt)//',a)') &
                             'Reading values of increment', restartInc - 1_pInt, 'from file' 
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
 call Utilities_constitutiveResponse(F,F,temperature,0.0_pReal,P,C,temp33_Real,.false.,math_I3) ! constitutive response with no deformation in no time to get reference stiffness
 if (restartInc == 1_pInt) then                                                                     ! use initial stiffness as reference stiffness
   temp3333_Real = C
 endif 
   
 call Utilities_updateGamma(temp3333_Real,.True.)
 
end subroutine basic_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the basic scheme with internal iterations
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function & 
  basic_solution(incInfo,guess,timeinc,timeinc_old,P_BC,F_BC,temperature_bc,rotation_BC)
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
   wgt, &
   mesh_ipCoordinates,&
   mesh_NcpElems, &
   mesh_deformedCoordsFFT
 use IO, only: &
   IO_write_JobBinaryFile, &
   IO_intOut
 use DAMASK_spectral_Utilities, only: &
   tBoundaryCondition, &
   field_real, &
   Utilities_forwardField, &
   Utilities_maskedCompliance, &
   Utilities_FFTforward, &
   Utilities_divergenceRMS, &
   Utilities_fourierConvolution, &
   Utilities_FFTbackward, &
   Utilities_updateGamma, &
   Utilities_constitutiveResponse, &
   Utilities_calculateRate, &
   debugRotation
 use FEsolving, only: &
   restartWrite, &
   restartRead, &
   terminallyIll
 use DAMASK_spectral_Utilities, only: &
   cutBack 
 
 implicit none
!--------------------------------------------------------------------------------------------------
! input data for solution
 real(pReal), intent(in) :: &
   timeinc, &                                                                                       !< increment in time for current solution
   timeinc_old, &                                                                                   !< increment in time of last increment
   temperature_bc                                                                                   !< temperature (I wonder, why it is intent(in)
 logical, intent(in) :: &
   guess                                                                                            !< if .false., assume homogeneous addon
 type(tBoundaryCondition),      intent(in) :: &
   P_BC,&                                                                                           !< stress boundary conditions
   F_BC                                                                                             !< deformation boundary conditions
 character(len=*), intent(in) :: &
   incInfo                                                                                          !< string with information of current increment for output to screen
 real(pReal), dimension(3,3), intent(in) :: &
   rotation_BC                                                                                      !< rotation load to lab
 
 real(pReal), dimension(3,3,3,3)        :: &
   S                                                                                                !< current average compliance 
 real(pReal), dimension(3,3)            :: &                            
   F_aim_lastIter, &                                                                                !< aim of last iteration
   P_av
 real(pReal), dimension(3,3,res(1),res(2),res(3)) :: P
!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
 real(pReal)   :: err_div, err_stress       
 integer(pInt) :: iter, row, column, i, j, k
 logical       :: ForwardData
 real(pReal)   :: & 
   defgradDet, &
   defgradDetMax, &
   defgradDetMin 
 real(pReal), dimension(3,3)            :: temp33_Real 

!--------------------------------------------------------------------------------------------------
! write restart information for spectral solver
 if (restartWrite) then
   write(6,'(a)') 'writing converged results for restart'
   call IO_write_jobBinaryFile(777,'convergedSpectralDefgrad',size(F))                        ! writing deformation gradient field to file
   write (777,rec=1) F
   close (777)
   call IO_write_jobBinaryFile(777,'convergedSpectralDefgrad_lastInc',size(F_lastInc))        ! writing F_lastInc field to file
   write (777,rec=1) F_lastInc
   close (777)
   call IO_write_jobBinaryFile(777,'F_aim',size(F_aim))
   write (777,rec=1) F_aim
   close(777)
   call IO_write_jobBinaryFile(777,'F_aim_lastInc',size(F_aim_lastInc))
   write (777,rec=1) F_aim_lastInc
   close(777)
   call IO_write_jobBinaryFile(777,'C',size(C))
   write (777,rec=1) C
   close(777)
   call IO_write_jobBinaryFile(777,'C_lastInc',size(C_lastInc))
   write (777,rec=1) C_lastInc
   close(777)
   call IO_write_jobBinaryFile(777,'F_aimDot',size(f_aimDot))
   write (777,rec=1) f_aimDot
   close(777)
 endif 

!--------------------------------------------------------------------------------------------------
! forward data 
 if (cutBack) then  
   F_aim = F_aim_lastInc
   F = F_lastInc
   C = C_lastInc
 else
   C_lastInc = C
   mesh_ipCoordinates = 0.0_pReal !reshape(mesh_deformedCoordsFFT(geomdim,&
                             !reshape(F,[3,3,res(1),res(2),res(3)])),[3,1,mesh_NcpElems])
   if (F_BC%myType=='l') then                                                                       ! calculate f_aimDot from given L and current F
     f_aimDot = F_BC%maskFloat * math_mul33x33(F_BC%values, F_aim)
   elseif(F_BC%myType=='fdot')   then                                                               ! f_aimDot is prescribed
     f_aimDot = F_BC%maskFloat * F_BC%values
   endif
   if (guess) f_aimDot  = f_aimDot + P_BC%maskFloat * (F_aim - F_aim_lastInc)/timeinc_old
   F_aim_lastInc = F_aim
   Fdot =  utilities_calculateRate(math_rotate_backward33(f_aimDot,rotation_BC), &
                                                        timeinc_old,guess,F_lastInc,F)
   F_lastInc = F
 endif
 F_aim = F_aim + f_aimDot * timeinc
 F = Utilities_forwardField(timeinc,math_rotate_backward33(F_aim,rotation_BC),F_lastInc,Fdot)                                           !I think F aim should be rotated here

!--------------------------------------------------------------------------------------------------
! update stiffness (and gamma operator)
 S = Utilities_maskedCompliance(rotation_BC,P_BC%maskLogical,C)
 if (update_gamma) call Utilities_updateGamma(C,restartWrite)
 
!--------------------------------------------------------------------------------------------------
! iteration till converged
 if (.not. restartRead) ForwardData = .True.
 iter = 0_pInt
 convergenceLoop: do while(iter < itmax)
   iter = iter + 1_pInt
!--------------------------------------------------------------------------------------------------
! report begin of new iteration
   write(6,'(/,a,3(a,'//IO_intOut(itmax)//'))') trim(incInfo), &
                    ' @ Iter. ', itmin, '≤',iter, '≤', itmax
   if (debugRotation) &
   write(6,'(a,/,3(3(f12.7,1x)/))',advance='no') 'deformation gradient aim (lab)=', &
                                        math_transpose33(math_rotate_backward33(F_aim,rotation_BC))
   write(6,'(a,/,3(3(f12.7,1x)/))',advance='no') 'deformation gradient aim =', &
                                        math_transpose33(F_aim)
!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
   F_aim_lastIter = F_aim
   basic_solution%termIll = .false.
   call Utilities_constitutiveResponse(F_lastInc,F,temperature,timeinc,&
                                 P,C,P_av,ForwardData,rotation_BC)
   basic_solution%termIll = terminallyIll
   terminallyIll = .false.
   ForwardData = .false.
   
!--------------------------------------------------------------------------------------------------
! stress BC handling
   F_aim = F_aim - math_mul3333xx33(S, ((P_av - P_BC%values)))                                      ! S = 0.0 for no bc
   err_stress = maxval(abs(P_BC%maskFloat * (P_av - P_BC%values)))                                  ! mask = 0.0 for no bc

!--------------------------------------------------------------------------------------------------
! updated deformation gradient using fix point algorithm of basic scheme
   field_real = 0.0_pReal
   field_real(1:res(1),1:res(2),1:res(3),1:3,1:3) = reshape(P,[res(1),res(2),res(3),3,3],&
                                                               order=[4,5,1,2,3]) ! field real has a different order
   call Utilities_FFTforward()
   err_div = Utilities_divergenceRMS()
   call Utilities_fourierConvolution(math_rotate_backward33(F_aim_lastIter-F_aim,rotation_BC))
   call Utilities_FFTbackward()
   F = F - reshape(field_real(1:res(1),1:res(2),1:res(3),1:3,1:3),shape(F),order=[3,4,5,1,2])                       ! F(x)^(n+1) = F(x)^(n) + correction;  *wgt: correcting for missing normalization
   basic_solution%converged = basic_Converged(err_div,P_av,err_stress,P_av)
   write(6,'(/,a)') '=========================================================================='
   if ((basic_solution%converged .and. iter >= itmin) .or. basic_solution%termIll) then
     basic_solution%iterationsNeeded = iter
     exit
   endif
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
