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
 private
#include <finclude/petscsys.h>
#include <finclude/petscdmda.h>
#include <finclude/petscsnes.h>

 character (len=*), parameter, public :: &
   DAMASK_spectral_SolverAL_label = 'al'
   
!--------------------------------------------------------------------------------------------------
! derived types 
 type tSolutionParams                                                                               !< @todo use here the type definition for a full loadcase including mask
   real(pReal), dimension(3,3) :: P_BC, rotation_BC
   real(pReal) :: timeinc
   real(pReal) :: temperature
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
   F_tau_lastInc, &                                                                                 !< field of previous incompatible deformation gradient 
   Fdot, &                                                                                          !< field of assumed rate of compatible deformation gradient
   F_tauDot                                                                                         !< field of assumed rate of incopatible deformation gradient

!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
 real(pReal), private, dimension(3,3) :: &
   F_aimDot, &                                                                                      !< assumed rate of average deformation gradient
   F_aim = math_I3, &                                                                               !< current prescribed deformation gradient
   F_aim_lastInc = math_I3, &                                                                       !< previous average deformation gradient
   P_av = 0.0_pReal                                                                                 !< average 1st Piola--Kirchhoff stress
 character(len=1024), private :: incInfo                                                            !< time and increment information
 real(pReal), private, dimension(3,3,3,3) :: &
   C = 0.0_pReal, C_minmaxAvg = 0.0_pReal, &                                                        !< current average stiffness
   C_lastInc = 0.0_pReal, &                                                                         !< previous average stiffness
   S = 0.0_pReal, &                                                                                 !< current compliance (filled up with zeros)
   C_scale = 0.0_pReal, &                             
   S_scale = 0.0_pReal
 
 real(pReal), private :: &
   err_stress, &                                                                                    !< deviation from stress BC
   err_f, &                                                                                         !< difference between compatible and incompatible deformation gradient
   err_p                                                                                            !< difference of stress resulting from compatible and incompatible F
 logical, private :: ForwardData
 integer(pInt), private :: reportIter = 0_pInt
 
 public :: &
   AL_init, &
   AL_solution, &
   AL_destroy
 external :: &
   VecDestroy, &
   DMDestroy, &
   DMDACreate3D, &
   DMCreateGlobalVector, &
   DMDASetLocalFunction, &
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
   MPI_Abort

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data, potentially from restart info
!> @todo use sourced allocation, e.g. allocate(Fdot,source = F_lastInc)
!--------------------------------------------------------------------------------------------------
subroutine AL_init(temperature)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 use IO, only: &
   IO_intOut, &
   IO_read_JobBinaryFile, &
   IO_write_JobBinaryFile, &
   IO_timeStamp
 use debug, only : &
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
   geomSize, &
   wgt
 use mesh, only: &
   mesh_ipCoordinates, &
   mesh_deformedCoordsFFT
 use math, only: &
   math_invSym3333
   
 implicit none
 real(pReal), intent(inout) :: &
   temperature
#include <finclude/petscdmda.h90>
#include <finclude/petscsnes.h90>
 real(pReal), dimension(:,:,:,:,:), allocatable :: P
 real(pReal), dimension(3,3) :: &
   temp33_Real = 0.0_pReal
 real(pReal), dimension(3,3,3,3) :: &
   temp3333_Real = 0.0_pReal, &
   temp3333_Real2 = 0.0_pReal

 PetscErrorCode :: ierr
 PetscObject :: dummy
 PetscScalar, pointer, dimension(:,:,:,:) :: xx_psc, F, F_tau
 
 call Utilities_init()
 write(6,'(/,a)') ' <<<+-  DAMASK_spectral_solverAL init  -+>>>'
 write(6,'(a)') ' $Id$'
 write(6,'(a16,a)')   ' Current time : ',IO_timeStamp()
#include "compilation_info.f90"

 allocate (P            (3,3,grid(1),grid(2),grid(3)),source = 0.0_pReal)
!--------------------------------------------------------------------------------------------------
! allocate global fields
 allocate (F_lastInc    (3,3,grid(1),grid(2),grid(3)),source = 0.0_pReal)
 allocate (Fdot         (3,3,grid(1),grid(2),grid(3)),source = 0.0_pReal)
 allocate (F_tau_lastInc(3,3,grid(1),grid(2),grid(3)),source = 0.0_pReal)
 allocate (F_tauDot     (3,3,grid(1),grid(2),grid(3)),source = 0.0_pReal)
    
!--------------------------------------------------------------------------------------------------
! PETSc Init
 call SNESCreate(PETSC_COMM_WORLD,snes,ierr); CHKERRQ(ierr)
 call DMDACreate3d(PETSC_COMM_WORLD,                               &
           DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, &
           DMDA_STENCIL_BOX,grid(1),grid(2),grid(3),PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, &
           18,1,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da,ierr)
 CHKERRQ(ierr)
 call DMCreateGlobalVector(da,solution_vec,ierr); CHKERRQ(ierr)
 call DMDASetLocalFunction(da,AL_formResidual,ierr); CHKERRQ(ierr)
 call SNESSetDM(snes,da,ierr); CHKERRQ(ierr)
 call SNESSetConvergenceTest(snes,AL_converged,dummy,PETSC_NULL_FUNCTION,ierr)
 CHKERRQ(ierr)
 call SNESSetFromOptions(snes,ierr); CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! init fields                 
 call DMDAVecGetArrayF90(da,solution_vec,xx_psc,ierr); CHKERRQ(ierr)                                ! places pointer xx_psc on PETSc data
 F => xx_psc(0:8,:,:,:)
 F_tau => xx_psc(9:17,:,:,:)
 if (restartInc == 1_pInt) then                                                                     ! no deformation (no restart)
   F_lastInc     = spread(spread(spread(math_I3,3,grid(1)),4,grid(2)),5,grid(3))                    ! initialize to identity
   F_tau_lastInc = F_lastInc
   F = reshape(F_lastInc,[9,grid(1),grid(2),grid(3)])
   F_tau = F
 elseif (restartInc > 1_pInt) then 
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
   F_aim         = reshape(sum(sum(sum(F,dim=4),dim=3),dim=2) * wgt, [3,3])                         ! average of F
   F_aim_lastInc = sum(sum(sum(F_lastInc,dim=5),dim=4),dim=3) * wgt                                 ! average of F_lastInc 
   call IO_read_jobBinaryFile(777,'F_tau',&
                                           trim(getSolverJobName()),size(F_tau))
   read (777,rec=1) F_tau
   close (777)
   call IO_read_jobBinaryFile(777,'F_tau_lastInc',&
                                        trim(getSolverJobName()),size(F_tau_lastInc))
   read (777,rec=1) F_tau_lastInc
   close (777)
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
   read (777,rec=1) C_minmaxAvg
   close (777)
 endif
 mesh_ipCoordinates = reshape(mesh_deformedCoordsFFT(geomSize,reshape(&
                                              F,[3,3,grid(1),grid(2),grid(3)])),[3,1,product(grid)])
 call Utilities_constitutiveResponse(F,F,temperature,0.0_pReal,P,temp3333_Real,temp3333_Real2,&
                                temp33_Real,.false.,math_I3)
 call DMDAVecRestoreArrayF90(da,solution_vec,xx_psc,ierr); CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! reference stiffness
 if (restartInc == 1_pInt) then                                                                     ! use initial stiffness as reference stiffness
   C_minmaxAvg = temp3333_Real2
   C = temp3333_Real
 endif 

 call Utilities_updateGamma(temp3333_Real2,.True.)
 C_scale = temp3333_Real2
 S_scale = math_invSym3333(temp3333_Real2)
 
end subroutine AL_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the AL scheme with internal iterations
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function &
  AL_solution(incInfoIn,guess,timeinc,timeinc_old,loadCaseTime,P_BC,F_BC,temperature_bc,rotation_BC)
 use numerics, only: &
   update_gamma, &
   itmax
 use math, only: &
   math_mul33x33 ,&
   math_rotate_backward33, &
   math_invSym3333
use mesh, only: &
   mesh_ipCoordinates, &
   mesh_deformedCoordsFFT
 use IO, only: &
   IO_write_JobBinaryFile
 use DAMASK_spectral_Utilities, only: &
   grid, &
   geomSize, &
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
  real(pReal), intent(in) :: &
   timeinc, &                                                                                       !< increment in time for current solution
   timeinc_old, &                                                                                   !< increment in time of last increment
   loadCaseTime, &                                                                                  !< remaining time of current load case
   temperature_bc
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
 PetscScalar, dimension(:,:,:,:), pointer :: xx_psc, F, F_tau
 PetscErrorCode :: ierr   
 SNESConvergedReason :: reason

 incInfo = incInfoIn
 
 call DMDAVecGetArrayF90(da,solution_vec,xx_psc,ierr)
 F => xx_psc(0:8,:,:,:)
 F_tau => xx_psc(9:17,:,:,:)
 
!--------------------------------------------------------------------------------------------------
! restart information for spectral solver
 if (restartWrite) then
   write(6,'(/,a)') ' writing converged results for restart'
   flush(6)
   call IO_write_jobBinaryFile(777,'F',size(F))                                                     ! writing deformation gradient field to file
   write (777,rec=1) F
   close (777)
   call IO_write_jobBinaryFile(777,'F_lastInc',size(F_lastInc))                                     ! writing F_lastInc field to file
   write (777,rec=1) F_lastInc
   close (777)
   call IO_write_jobBinaryFile(777,'F_tau',size(F_tau))                                             ! writing deformation gradient field to file
   write (777,rec=1) F_tau
   close (777)
   call IO_write_jobBinaryFile(777,'F_tau_lastInc',size(F_tau_lastInc))                             ! writing F_lastInc field to file
   write (777,rec=1) F_tau_lastInc
   close (777)
   call IO_write_jobBinaryFile(777,'F_aimDot',size(F_aimDot))
   write (777,rec=1) F_aimDot
   close(777)
   call IO_write_jobBinaryFile(777,'C',size(C))
   write (777,rec=1) C
   close(777)
   call IO_write_jobBinaryFile(777,'C_lastInc',size(C_lastInc))
   write (777,rec=1) C_lastInc
   close(777)
 endif 
 AL_solution%converged =.false.

 if (cutBack) then 
   F_aim = F_aim_lastInc
   F_tau= reshape(F_tau_lastInc,[9,grid(1),grid(2),grid(3)]) 
   F    = reshape(F_lastInc,    [9,grid(1),grid(2),grid(3)]) 
   C = C_lastInc
 else
   C_lastInc = C
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
   mesh_ipCoordinates = reshape(mesh_deformedCoordsFFT(geomSize,reshape(&
                                            F,[3,3,grid(1),grid(2),grid(3)])),[3,1,product(grid)])
   Fdot =  Utilities_calculateRate(math_rotate_backward33(f_aimDot,rotation_BC), &
                  timeinc_old,guess,F_lastInc,reshape(F,[3,3,grid(1),grid(2),grid(3)]))
   F_tauDot =  Utilities_calculateRate(math_rotate_backward33(f_aimDot,rotation_BC), &
                  timeinc_old,guess,F_tau_lastInc,reshape(F_tau,[3,3,grid(1),grid(2),grid(3)]))  
                
   F_lastInc     = reshape(F,       [3,3,grid(1),grid(2),grid(3)])
   F_tau_lastInc = reshape(F_tau,[3,3,grid(1),grid(2),grid(3)])
 endif
 F_aim = F_aim + f_aimDot * timeinc

!--------------------------------------------------------------------------------------------------
! update local deformation gradient
 F     = reshape(Utilities_forwardField(timeinc,F_lastInc,Fdot, &                                   ! ensure that it matches rotated F_aim
                               math_rotate_backward33(F_aim,rotation_BC)),[9,grid(1),grid(2),grid(3)])
 F_tau = reshape(Utilities_forwardField(timeinc,F_tau_lastInc,F_taudot),  [9,grid(1),grid(2),grid(3)]) ! does not have any average value as boundary condition
 call DMDAVecRestoreArrayF90(da,solution_vec,xx_psc,ierr)
 CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! update stiffness (and gamma operator)
 S = Utilities_maskedCompliance(rotation_BC,P_BC%maskLogical,C)
 if (update_gamma) then
   call Utilities_updateGamma(C_minmaxAvg,restartWrite)
   C_scale = C_minmaxAvg
   S_scale = math_invSym3333(C_minmaxAvg)
 endif  
 
 ForwardData = .True.
 mask_stress = P_BC%maskFloat
 params%P_BC = P_BC%values
 params%rotation_BC = rotation_BC
 params%timeinc = timeinc
 params%temperature = temperature_bc

 call SNESSolve(snes,PETSC_NULL_OBJECT,solution_vec,ierr)
 CHKERRQ(ierr)
 call SNESGetConvergedReason(snes,reason,ierr)
 CHKERRQ(ierr)
 
 AL_solution%termIll = terminallyIll
 terminallyIll = .false.
 if (reason < 1 ) then 
   AL_solution%converged = .false.
   AL_solution%iterationsNeeded = itmax
 else
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
   itmin, &
   polarAlpha, &
   polarBeta
 use IO, only: &
   IO_intOut
 use math, only: &
   math_rotate_backward33, &
   math_transpose33, &
   math_mul3333xx33, &
   math_invSym3333
 use DAMASK_spectral_Utilities, only: &
   grid, &
   wgt, &
   field_real, &
   Utilities_FFTforward, &
   Utilities_fourierConvolution, &
   Utilities_FFTbackward, &
   Utilities_constitutiveResponse
 use debug, only: &
   debug_level, &
   debug_spectral, &
   debug_spectralRotation
 use homogenization, only: &
   materialpoint_P, &
   materialpoint_dPdF

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
   F_tau, &
   residual_F, &
   residual_F_tau
 PetscInt :: &
   iter, &
   nfuncs
 PetscObject :: dummy
 PetscErrorCode :: ierr
 integer(pInt) :: &
   i, j, k, n_ele

 F              => x_scal(1:3,1:3,1,&
  XG_RANGE,YG_RANGE,ZG_RANGE)
 F_tau          => x_scal(1:3,1:3,2,&
  XG_RANGE,YG_RANGE,ZG_RANGE)
 residual_F     => f_scal(1:3,1:3,1,&
  X_RANGE,Y_RANGE,Z_RANGE)
 residual_F_tau => f_scal(1:3,1:3,2,&
  X_RANGE,Y_RANGE,Z_RANGE)
 
 call SNESGetNumberFunctionEvals(snes,nfuncs,ierr); CHKERRQ(ierr)
 call SNESGetIterationNumber(snes,iter,ierr); CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! report begin of new iteration
 if (iter == 0 .and. nfuncs == 0) then                                                             ! new increment
   reportIter = -1_pInt
 endif
 if (reportIter <= iter) then                                                                      ! new iteration
   reportIter = reportIter + 1_pInt
   write(6,'(1x,a,3(a,'//IO_intOut(itmax)//'))') trim(incInfo), &
                    ' @ Iteration ', itmin, '≤',reportIter, '≤', itmax
   if (iand(debug_level(debug_spectral),debug_spectralRotation) /= 0) &
   write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') ' deformation gradient aim (lab) =', &
                                 math_transpose33(math_rotate_backward33(F_aim,params%rotation_BC))
   write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') ' deformation gradient aim =', &
                                 math_transpose33(F_aim)
   flush(6)
 endif
  
!--------------------------------------------------------------------------------------------------
! 
 field_real = 0.0_pReal
 do k = 1_pInt, grid(3); do j = 1_pInt, grid(2); do i = 1_pInt, grid(1)
   field_real(i,j,k,1:3,1:3) = math_mul3333xx33(C_scale,(polarAlpha + polarBeta)*F(1:3,1:3,i,j,k) - &
                                                               (polarAlpha)*F_tau(1:3,1:3,i,j,k))
 enddo; enddo; enddo
 
!--------------------------------------------------------------------------------------------------
! doing convolution in Fourier space 
 call Utilities_FFTforward()
 call Utilities_fourierConvolution(math_rotate_backward33(polarBeta*F_aim,params%rotation_BC)) 
 call Utilities_FFTbackward()
 
!--------------------------------------------------------------------------------------------------
! constructing residual                         
 residual_F_tau = polarBeta*F - reshape(field_real(1:grid(1),1:grid(2),1:grid(3),1:3,1:3),&
                                      [3,3,grid(1),grid(2),grid(3)],order=[3,4,5,1,2])

!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
 call Utilities_constitutiveResponse(F_lastInc,F - residual_F_tau/polarBeta,params%temperature,params%timeinc, &
                                     residual_F,C,C_minmaxAvg,P_av,ForwardData,params%rotation_BC)
 ForwardData = .False.
  
!--------------------------------------------------------------------------------------------------
! stress BC handling
 F_aim = F_aim - math_mul3333xx33(S, ((P_av - params%P_BC)))                                        ! S = 0.0 for no bc
 err_stress = maxval(abs(mask_stress * (P_av - params%P_BC)))                                       ! mask = 0.0 for no bc
 
!--------------------------------------------------------------------------------------------------
! constructing residual
 n_ele = 0_pInt
 err_p = 0.0_pReal
 do k = 1_pInt, grid(3); do j = 1_pInt, grid(2); do i = 1_pInt, grid(1)
   n_ele = n_ele + 1_pInt
   err_p = err_p + sum((math_mul3333xx33(S_scale,residual_F(1:3,1:3,i,j,k)) - &
                        (F_tau(1:3,1:3,i,j,k) - &
                         F(1:3,1:3,i,j,k) + residual_F_tau(1:3,1:3,i,j,k)/polarBeta))**2.0_pReal)
   residual_F(1:3,1:3,i,j,k) = math_mul3333xx33(math_invSym3333(materialpoint_dPdF(:,:,:,:,1,n_ele) + C_scale), &
                                                residual_F(1:3,1:3,i,j,k) - &
                                                math_mul3333xx33(C_scale,F_tau(1:3,1:3,i,j,k) - F(1:3,1:3,i,j,k))) &
                               + residual_F_tau(1:3,1:3,i,j,k)
 enddo; enddo; enddo
 
!--------------------------------------------------------------------------------------------------
! calculating errors  
 err_f = wgt*sqrt(sum(residual_F_tau**2.0_pReal))/polarBeta
 err_p = wgt*sqrt(err_p) 
   
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
 use FEsolving, only: &
   terminallyIll

 implicit none
 SNES :: snes_local
 PetscInt :: it
 PetscReal :: &
   xnorm, &
   snorm, &
   fnorm
 SNESConvergedReason :: reason
 PetscObject :: dummy
 PetscErrorCode ::ierr
 logical :: Converged
 real(pReal) :: err_stress_tol
 
 err_stress_tol = max(maxval(abs(P_av))*err_stress_tolrel,err_stress_tolabs)
 Converged = (it >= itmin .and. &
               all([ err_f/err_f_tol, &
                     err_p/err_p_tol, &
                     err_stress/err_stress_tol] < 1.0_pReal)) &
             .or.    terminallyIll     
 
 if (Converged) then
   reason = 1
 elseif (it >= itmax) then
   reason = -1
 else  
   reason = 0
 endif 
 write(6,'(1/,a)') ' ... reporting ....................................................'
 write(6,'(/,a,f8.2,a,es11.5,a,es11.4,a)') ' mismatch F =      ', &
                     err_f/err_f_tol, &
                ' (',err_f,' -,  tol =',err_f_tol,')'
 write(6,'(a,f8.2,a,es11.5,a,es11.4,a)')   ' mismatch P =      ', &
                     err_p/err_p_tol, &
                ' (',err_p,' -,  tol =',err_p_tol,')'
 write(6,'(a,f8.2,a,es11.5,a,es11.4,a)')   ' error stress BC = ', &
                   err_stress/err_stress_tol, ' (',err_stress, ' Pa, tol =',err_stress_tol,')' 
 write(6,'(/,a)') ' =========================================================================='
 flush(6)

end subroutine AL_converged

!--------------------------------------------------------------------------------------------------
!> @brief destroy routine
!--------------------------------------------------------------------------------------------------
subroutine AL_destroy()
 use DAMASK_spectral_Utilities, only: &
   Utilities_destroy
 
 implicit none
 PetscErrorCode :: ierr

 call VecDestroy(solution_vec,ierr); CHKERRQ(ierr)
 call SNESDestroy(snes,ierr); CHKERRQ(ierr)
 call DMDestroy(da,ierr); CHKERRQ(ierr)
 call PetscFinalize(ierr); CHKERRQ(ierr)
 call Utilities_destroy()

end subroutine AL_destroy

end module DAMASK_spectral_SolverAL
