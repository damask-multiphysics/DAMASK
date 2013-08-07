! Copyright 2011-13 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced MAterial Simulation Kit.
!
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DAMASK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DAMASK. If not, see <http://www.gnu.org/licenses/>.
!
!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Basic scheme solver
!> @details this solver follows closely the original large strain formulation presented by
!> Suquet. The iterative procedure is solved using a fix-point iteration.
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
 private
 character (len=*), parameter, public :: &
   DAMASK_spectral_SolverBasic_label = 'basic'

!--------------------------------------------------------------------------------------------------
! pointwise global data
 real(pReal),  private,  dimension(:,:,:,:,:), allocatable ::  &
   F, &                                                                                             !< deformation gradient field
   F_lastInc, &                                                                                     !< deformation gradient field last increment
   Fdot                                                                                             !< assumed rate for F n to F n+1

!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
 real(pReal), private, dimension(3,3) :: &
   F_aim = math_I3, &                                                                               !< deformation gradient aim
   F_aim_lastInc = math_I3, &                                                                       !< deformation gradient aim last increment
   F_aimDot = 0.0_pReal                                                                             !< assumed rate
 real(pReal), private,dimension(3,3,3,3) :: &
   C = 0.0_pReal, C_minmaxAvg = 0.0_pReal, &                                                        !< average stiffness
   C_lastInc = 0.0_pReal                                                                            !< average stiffness last increment
 public :: &
   basic_init, &
   basic_solution, &
   basic_destroy
contains
 
!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
subroutine basic_init(temperature)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 
 use IO, only: &
   IO_read_JobBinaryFile, &
   IO_write_JobBinaryFile, &
   IO_intOut, &
   IO_timeStamp
 use debug, only: &
   debug_level, &
   debug_spectral, &
   debug_spectralRestart
 use FEsolving, only: &
   restartInc
 use DAMASK_interface, only: &
   getSolverJobName
 use DAMASK_spectral_Utilities, only: &
   Utilities_init, &
   Utilities_constitutiveResponse, &
   Utilities_updateGamma, &
   grid, &
   wgt, &
   geomSize
 use mesh, only: &
   mesh_ipCoordinates, &
   mesh_deformedCoordsFFT

 implicit none
 real(pReal), intent(inout) :: &
   temperature
 real(pReal), dimension(:,:,:,:,:), allocatable :: P
 real(pReal), dimension(3,3) :: &
   temp33_Real = 0.0_pReal
 real(pReal), dimension(3,3,3,3) :: &
   temp3333_Real
 
 call Utilities_Init()
 write(6,'(/,a)')   ' <<<+-  DAMASK_spectral_solverBasic init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 allocate (P         (3,3,grid(1),  grid(2),grid(3)),  source = 0.0_pReal)
!--------------------------------------------------------------------------------------------------
! allocate global fields
 allocate (F         (3,3,grid(1),  grid(2),grid(3)),  source = 0.0_pReal)
 allocate (F_lastInc (3,3,grid(1),  grid(2),grid(3)),  source = 0.0_pReal)
 allocate (Fdot      (3,3,grid(1),  grid(2),grid(3)),  source = 0.0_pReal)
   
!--------------------------------------------------------------------------------------------------
! init fields and average quantities
 if (restartInc == 1_pInt) then                                                                     ! no deformation (no restart)
   F         = spread(spread(spread(math_I3,3,grid(1)),4,grid(2)),5,grid(3))                        ! initialize to identity
   F_lastInc = F
 elseif (restartInc > 1_pInt) then                                                                  ! using old values from file                                                      
   if (iand(debug_level(debug_spectral),debug_spectralRestart)/= 0) &
     write(6,'(/,a,'//IO_intOut(restartInc-1_pInt)//',a)') &
     'reading values of increment', restartInc - 1_pInt, 'from file'
   flush(6)
   call IO_read_jobBinaryFile(777,'F',&
                                                  trim(getSolverJobName()),size(F))
   read (777,rec=1) F
   close (777)
   call IO_read_jobBinaryFile(777,'F_lastInc',&
                                                  trim(getSolverJobName()),size(F_lastInc))
   read (777,rec=1) F_lastInc
   close (777)
   
   F_aim         = sum(sum(sum(F,dim=5),dim=4),dim=3) * wgt                                         ! average of F
   F_aim_lastInc = sum(sum(sum(F_lastInc,dim=5),dim=4),dim=3) * wgt                                 ! average of F_lastInc 
   
   call IO_read_jobBinaryFile(777,'F_aimDot',trim(getSolverJobName()),size(f_aimDot))
   read (777,rec=1) f_aimDot
   close (777)
   call IO_read_jobBinaryFile(777,'C',trim(getSolverJobName()),size(C))
   read (777,rec=1) C
   close (777)
   call IO_read_jobBinaryFile(777,'C_lastInc',trim(getSolverJobName()),size(C_lastInc))
   read (777,rec=1) C_lastInc
   close (777)
   call IO_read_jobBinaryFile(777,'C_ref',trim(getSolverJobName()),size(temp3333_Real))
   read (777,rec=1) temp3333_Real
   close (777)
 endif
 mesh_ipCoordinates = reshape(mesh_deformedCoordsFFT(geomSize,F),[3,1,product(grid)])
 call Utilities_constitutiveResponse(F,F,temperature,0.0_pReal,P,C,C_minmaxAvg,&
                                     temp33_Real,.false.,math_I3)                                   ! constitutive response with no deformation in no time to get reference stiffness
 if (restartInc == 1_pInt) then                                                                     ! use initial stiffness as reference stiffness
   temp3333_Real = C_minmaxAvg
 endif 
   
 call Utilities_updateGamma(temp3333_Real,.True.)
 
end subroutine basic_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the basic scheme with internal iterations
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function basic_solution(&
             incInfo,guess,timeinc,timeinc_old,loadCaseTime,P_BC,F_BC,temperature_bc,rotation_BC)
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
   mesh_ipCoordinates,&
   mesh_deformedCoordsFFT
 use IO, only: &
   IO_write_JobBinaryFile, &
   IO_intOut
 use debug, only: &
   debug_level, &
   debug_spectral, &
   debug_spectralRotation
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
   grid,&
   geomSize, &
   wgt
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
   loadCaseTime                                                                                     !< remaining time of current load case
 logical, intent(in) :: &
   guess                                                                                            !< if .false., assume homogeneous addon
 type(tBoundaryCondition),      intent(in) :: &
   P_BC,&                                                                                           !< stress boundary conditions
   F_BC                                                                                             !< deformation boundary conditions
 character(len=*), intent(in) :: &
   incInfo                                                                                          !< string with information of current increment for output to screen
 real(pReal), dimension(3,3), intent(in) :: &
   rotation_BC                                                                                      !< rotation load to lab
 real(pReal), intent(inout) :: &
   temperature_bc  
 real(pReal), dimension(3,3,3,3)        :: &
   S                                                                                                !< current average compliance 
 real(pReal), dimension(3,3)            :: &                            
   F_aim_lastIter, &                                                                                !< aim of last iteration
   P_av
 real(pReal), dimension(3,3,grid(1),grid(2),grid(3)) :: P
!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
 real(pReal)   :: err_div, err_stress       
 integer(pInt) :: iter
 logical       :: ForwardData

!--------------------------------------------------------------------------------------------------
! write restart information for spectral solver
 if (restartWrite) then
   write(6,'(/,a)') ' writing converged results for restart'
   flush(6)
   call IO_write_jobBinaryFile(777,'F',size(F))                                                     ! writing deformation gradient field to file
   write (777,rec=1) F
   close (777)
   call IO_write_jobBinaryFile(777,'F_lastInc',size(F_lastInc))                                     ! writing F_lastInc field to file
   write (777,rec=1) F_lastInc
   close (777)
   call IO_write_jobBinaryFile(777,'F_aimDot',size(f_aimDot))
   write (777,rec=1) f_aimDot
   close(777)
   call IO_write_jobBinaryFile(777,'C',size(C))
   write (777,rec=1) C
   close(777)
   call IO_write_jobBinaryFile(777,'C_lastInc',size(C_lastInc))
   write (777,rec=1) C_lastInc
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
   mesh_ipCoordinates = reshape(mesh_deformedCoordsFFT(geomSize,F),[3,1,product(grid)])

!--------------------------------------------------------------------------------------------------
! calculate rate for aim
   if (F_BC%myType=='l') then                                                                       ! calculate f_aimDot from given L and current F
     f_aimDot = F_BC%maskFloat * math_mul33x33(F_BC%values, F_aim)
   elseif(F_BC%myType=='fdot')   then                                                               ! f_aimDot is prescribed
     f_aimDot = F_BC%maskFloat * F_BC%values
   elseif(F_BC%myType=='f') then                                                                    ! aim at end of load case is prescribed
     f_aimDot = F_BC%maskFloat * (F_BC%values - F_aim)/loadCaseTime
   endif
   if (guess) f_aimDot  = f_aimDot + P_BC%maskFloat * (F_aim - F_aim_lastInc)/timeinc_old
   F_aim_lastInc = F_aim

!--------------------------------------------------------------------------------------------------
! update rate and forward last inc
   Fdot =  utilities_calculateRate(math_rotate_backward33(f_aimDot,rotation_BC), &
                                                        timeinc_old,guess,F_lastInc,F)
   F_lastInc = F
 endif
 F_aim = F_aim + f_aimDot * timeinc
 F = Utilities_forwardField(timeinc,F_lastInc,Fdot,math_rotate_backward33(F_aim,rotation_BC))

!--------------------------------------------------------------------------------------------------
! update stiffness (and gamma operator)
 S = Utilities_maskedCompliance(rotation_BC,P_BC%maskLogical,C)
 if (update_gamma) call Utilities_updateGamma(C_minmaxAvg,restartWrite)
 
!--------------------------------------------------------------------------------------------------
! iteration till converged
 if (.not. restartRead) ForwardData = .True.
 iter = 0_pInt
 basic_solution%iterationsNeeded = itmax
 convergenceLoop: do while(iter < itmax)
   iter = iter + 1_pInt
!--------------------------------------------------------------------------------------------------
! report begin of new iteration
   write(6,'(1x,a,3(a,'//IO_intOut(itmax)//'))') trim(incInfo), &
                    ' @ Iteration ', itmin, '≤',iter, '≤', itmax
   if (iand(debug_level(debug_spectral),debug_spectralRotation) /= 0) &
   write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') ' deformation gradient aim (lab)=', &
                                        math_transpose33(math_rotate_backward33(F_aim,rotation_BC))
   write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') ' deformation gradient aim =', &
                                        math_transpose33(F_aim)
   flush(6)
!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
   F_aim_lastIter = F_aim
   basic_solution%termIll = .false.
   call Utilities_constitutiveResponse(F_lastInc,F,temperature_bc,timeinc,&
                                 P,C,C_minmaxAvg,P_av,ForwardData,rotation_BC)
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
   field_real(1:grid(1),1:grid(2),1:grid(3),1:3,1:3) = reshape(P,[grid(1),grid(2),grid(3),3,3],&
                                                               order=[4,5,1,2,3])                   ! field real has a different order
   call Utilities_FFTforward()
   err_div = Utilities_divergenceRMS()
   call Utilities_fourierConvolution(math_rotate_backward33(F_aim_lastIter-F_aim,rotation_BC))
   call Utilities_FFTbackward()
   F = F - reshape(field_real(1:grid(1),1:grid(2),1:grid(3),1:3,1:3),shape(F),order=[3,4,5,1,2])    ! F(x)^(n+1) = F(x)^(n) + correction;  *wgt: correcting for missing normalization
   basic_solution%converged = basic_Converged(err_div,P_av,err_stress)
   write(6,'(/,a)') ' ==========================================================================='
   flush(6)
   if ((basic_solution%converged .and. iter >= itmin) .or. basic_solution%termIll) then
     basic_solution%iterationsNeeded = iter
     exit
   endif
 enddo convergenceLoop

end function basic_solution


!--------------------------------------------------------------------------------------------------
!> @brief convergence check for basic scheme based on div of P and deviation from stress aim
!--------------------------------------------------------------------------------------------------
logical function basic_Converged(err_div,P_av,err_stress)
 use numerics, only: &
   itmin, &
   err_div_tolRel, &
   err_div_tolAbs, &
   err_stress_tolrel, &
   err_stress_tolabs
 use math, only: &
   math_mul33x33, &
   math_eigenvalues33, &
   math_transpose33
    
 implicit none
 real(pReal), dimension(3,3), intent(in) :: &
   P_av
 
 real(pReal), intent(in) :: &
   err_div, &
   err_stress
 
 real(pReal) :: &
   stressTol, &
   divTol
  
 divTol  = max(maxval(abs(P_av))*err_div_tolRel,err_div_tolAbs)
 stressTol = max(maxval(abs(P_av))*err_stress_tolrel,err_stress_tolabs)
 
 basic_Converged = all([ err_div/divTol,&
                         err_stress/stressTol    ]  < 1.0_pReal)

!--------------------------------------------------------------------------------------------------
! report
 write(6,'(1/,a)') ' ... reporting .............................................................'
 write(6,'(1/,a,f12.2,a,es8.2,a,es9.2,a)') ' error divergence = ', &
            err_div/divTol,  ' (',err_div,' / m, tol =',divTol,')'
 write(6,'(a,f12.2,a,es8.2,a,es9.2,a)')   ' error stress BC =  ', &
            err_stress/stressTol, ' (',err_stress, ' Pa,  tol =',stressTol,')' 
 flush(6) 

end function basic_Converged


!--------------------------------------------------------------------------------------------------
!> @brief does the cleaning up after the simulation has finished
!--------------------------------------------------------------------------------------------------
subroutine basic_destroy()
 
use DAMASK_spectral_Utilities, only: &
  Utilities_destroy
 
 implicit none
 call Utilities_destroy()

end subroutine basic_destroy


end module DAMASK_spectral_SolverBasic
