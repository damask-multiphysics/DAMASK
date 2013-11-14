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
   phaseFieldDataBin, &
   maxPhaseFields

 implicit none
 private
#include <finclude/petscsys.h>
#include <finclude/petscdmda.h>
#include <finclude/petscsnes.h>
 character (len=*), parameter, public :: &
   DAMASK_spectral_SolverBasicPETSC_label = 'basicpetsc'
   
!--------------------------------------------------------------------------------------------------
! derived types
 type tSolutionParams 
   real(pReal), dimension(3,3) :: P_BC, rotation_BC
   real(pReal) :: timeinc
   real(pReal) :: timeincOld
   real(pReal) :: temperature
   real(pReal) :: density
   integer(pInt) :: nActivePhaseFields
   type(phaseFieldDataBin) :: phaseFieldData(maxPhaseFields)
 end type tSolutionParams
 
 type(tSolutionParams), private :: params

!--------------------------------------------------------------------------------------------------
! PETSc data
 DM,   private :: da
 SNES, private :: snes
 Vec,  private :: solution_vec

!--------------------------------------------------------------------------------------------------
! common pointwise data
 real(pReal), private, dimension(:,:,:,:,:), allocatable ::  F_lastInc, Fdot, F_lastInc2
 real(pReal), private, dimension(:,:,:,:), allocatable :: &
   phaseFieldRHS_lastInc, &
   phaseField_lastInc, &
   phaseFieldRHS, &
   phaseFieldDot 
 complex(pReal), private, dimension(:,:,:,:,:), allocatable :: inertiaField_fourier

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
 real(pReal), private :: err_stress, err_div, err_divPrev, err_divDummy
 real(pReal), private, dimension(:), allocatable :: err_phaseField, phaseField_Avg
 logical, private :: ForwardData
 integer(pInt), private :: &
   totalIter = 0_pInt                                                                               !< total iteration in current increment
 real(pReal), private, dimension(3,3) :: mask_stress = 0.0_pReal

 public :: &
   basicPETSc_init, &
   basicPETSc_solution ,&
   basicPETSc_destroy
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
!--------------------------------------------------------------------------------------------------
subroutine basicPETSc_init(temperature,nActivePhaseFields,phaseFieldData)
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
 use DAMASK_interface, only: &
   getSolverJobName
 use DAMASK_spectral_Utilities, only: &
   Utilities_init, &
   Utilities_constitutiveResponse, &
   Utilities_updateGamma, &
   grid, &
   grid1Red, &
   wgt, &
   geomSize 
 use mesh, only: &
   mesh_ipCoordinates, &
   mesh_deformedCoordsFFT
 use math, only: &
   math_invSym3333
   
 implicit none
 integer(pInt), intent(in) :: nActivePhaseFields
 type(phaseFieldDataBin), intent(in) :: phaseFieldData(nActivePhaseFields)
 real(pReal), intent(inOut) :: temperature
#include <finclude/petscdmda.h90>
#include <finclude/petscsnes.h90>
#include <finclude/petscvec.h>
 real(pReal), dimension(:,:,:,:,:), allocatable :: P
 PetscScalar,  dimension(:,:,:,:), pointer     ::  xx_psc, F
 PetscErrorCode :: ierr
 PetscObject    :: dummy
 real(pReal), dimension(3,3) :: &
   temp33_Real = 0.0_pReal
 real(pReal), dimension(3,3,3,3) :: &
   temp3333_Real = 0.0_pReal
 KSP :: ksp
 integer(pInt) :: i

 call Utilities_init()
 write(6,'(/,a)') ' <<<+-  DAMASK_spectral_solverBasicPETSc init  -+>>>'
 write(6,'(a)') ' $Id$'
 write(6,'(a15,a)')   ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
!--------------------------------------------------------------------------------------------------
! allocate global fields
 allocate (P         (3,3,grid(1),grid(2),grid(3)),source = 0.0_pReal)
 allocate (F_lastInc (3,3,grid(1),grid(2),grid(3)),source = 0.0_pReal)
 allocate (F_lastInc2(3,3,grid(1),grid(2),grid(3)),source = 0.0_pReal)
 allocate (Fdot      (3,3,grid(1),grid(2),grid(3)),source = 0.0_pReal)
 allocate (inertiaField_fourier (grid1Red,grid(2),grid(3),3,3),source = cmplx(0.0_pReal,0.0_pReal,pReal))
 allocate (phaseFieldRHS_lastInc (nActivePhaseFields,grid(1),grid(2),grid(3)),source = 0.0_pReal)
 allocate (phaseField_lastInc    (nActivePhaseFields,grid(1),grid(2),grid(3)),source = 0.0_pReal)
 allocate (phaseFieldDot         (nActivePhaseFields,grid(1),grid(2),grid(3)),source = 0.0_pReal)
 allocate (phaseFieldRHS         (nActivePhaseFields,grid(1),grid(2),grid(3)),source = 0.0_pReal)
 allocate (err_phaseField(nActivePhaseFields), source = 0.0_pReal) 
 allocate (phaseField_Avg(nActivePhaseFields), source = 0.0_pReal)
    
!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
 call SNESCreate(PETSC_COMM_WORLD,snes,ierr); CHKERRQ(ierr)
 call DMDACreate3d(PETSC_COMM_WORLD, &
        DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, &
        DMDA_STENCIL_BOX,grid(1),grid(2),grid(3),PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, &
        9+nActivePhaseFields,1,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da,ierr)
   CHKERRQ(ierr)
 call DMCreateGlobalVector(da,solution_vec,ierr); CHKERRQ(ierr)
 !call DMDASNESSetFunctionLocal(da,INSERT_VALUES,BasicPETSC_formResidual,dummy,ierr); CHKERRQ(ierr)     ! needed for newer versions of petsc
 call DMDASetLocalFunction(da,BasicPETSC_formResidual,ierr); CHKERRQ(ierr)
 call SNESSetDM(snes,da,ierr); CHKERRQ(ierr)
 call SNESSetConvergenceTest(snes,BasicPETSC_converged,dummy,PETSC_NULL_FUNCTION,ierr); CHKERRQ(ierr)
 call SNESGetKSP(snes,ksp,ierr); CHKERRQ(ierr)
 call KSPSetConvergenceTest(ksp,BasicPETSC_convergedKSP,dummy,PETSC_NULL_FUNCTION,ierr); CHKERRQ(ierr)
 call SNESSetFromOptions(snes,ierr); CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! init fields                 
 call DMDAVecGetArrayF90(da,solution_vec,xx_psc,ierr); CHKERRQ(ierr)                                 ! get the data out of PETSc to work with
 F => xx_psc(0:8,:,:,:)
 if (restartInc == 1_pInt) then                                                                      ! no deformation (no restart)
   F_lastInc = spread(spread(spread(math_I3,3,grid(1)),4,grid(2)),5,grid(3))                         ! initialize to identity
   xx_psc(0:8,:,:,:) = reshape(F_lastInc,[9,grid(1),grid(2),grid(3)])
   F_lastInc2 = F_lastInc
   do i = 1, nActivePhaseFields
     xx_psc(8+i,:,:,:)           = phaseFieldData(i)%phaseField0
     phaseField_lastInc(i,:,:,:) = phaseFieldData(i)%phaseField0
   enddo
 elseif (restartInc > 1_pInt) then                                                                  ! using old values from file                                                      
   if (iand(debug_level(debug_spectral),debug_spectralRestart)/= 0) &
     write(6,'(/,a,'//IO_intOut(restartInc-1_pInt)//',a)') &
     'reading values of increment', restartInc - 1_pInt, 'from file'
   flush(6)
   call IO_read_realFile(777,'F',trim(getSolverJobName()),size(F))
   read (777,rec=1) F
   close (777)
   call IO_read_realFile(777,'F_lastInc',trim(getSolverJobName()),size(F_lastInc))
   read (777,rec=1) F_lastInc
   close (777)
   call IO_read_realFile(777,'F_lastInc2',trim(getSolverJobName()),size(F_lastInc2))
   read (777,rec=1) F_lastInc2
   close (777)
   F_aim         = reshape(sum(sum(sum(F,dim=4),dim=3),dim=2) * wgt, [3,3])                         ! average of F
   F_aim_lastInc = sum(sum(sum(F_lastInc,dim=5),dim=4),dim=3) * wgt                                 ! average of F_lastInc 

   call IO_read_realFile(777,'F_aimDot',trim(getSolverJobName()),size(f_aimDot))
   read (777,rec=1) f_aimDot
   close (777)
   call IO_read_realFile(777,'C_volAvg',trim(getSolverJobName()),size(C_volAvg))
   read (777,rec=1) C_volAvg
   close (777)
   call IO_read_realFile(777,'C_volAvgLastInc',trim(getSolverJobName()),size(C_volAvgLastInc))
   read (777,rec=1) C_volAvgLastInc
   close (777)
   call IO_read_realFile(777,'C_ref',trim(getSolverJobName()),size(temp3333_Real))
   read (777,rec=1) temp3333_Real
   close (777)
 endif
 mesh_ipCoordinates = reshape(mesh_deformedCoordsFFT(geomSize,reshape(&
                                              F,[3,3,grid(1),grid(2),grid(3)])),[3,1,product(grid)])
 call Utilities_constitutiveResponse(&
    reshape(F,[3,3,grid(1),grid(2),grid(3)]),&
    reshape(F,[3,3,grid(1),grid(2),grid(3)]),&
    temperature,0.0_pReal,P,C_volAvg,C_minmaxAvg,temp33_Real,.false.,math_I3)
 call DMDAVecRestoreArrayF90(da,solution_vec,xx_psc,ierr); CHKERRQ(ierr)                                 ! write data back into PETSc
 if (restartInc == 1_pInt) then                                                                     ! use initial stiffness as reference stiffness
   temp3333_Real = C_minMaxAvg
 endif 
   
 call Utilities_updateGamma(temp3333_Real,.True.)

end subroutine basicPETSc_init
  
!--------------------------------------------------------------------------------------------------
!> @brief solution for the Basic PETSC scheme with internal iterations
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function basicPETSc_solution( &
             incInfoIn,guess,timeinc,timeinc_old,loadCaseTime,P_BC,F_BC,temperature_bc,rotation_BC,density, &
             nActivePhaseFields,phaseFieldData)
 use numerics, only: &
   update_gamma, &
   itmax
 use math, only: &
   math_mul33x33 ,&
   math_rotate_backward33
 use mesh, only: &
   mesh_ipCoordinates,&
   mesh_deformedCoordsFFT
 use IO, only: &
   IO_write_JobRealFile
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
 use homogenization, only: &                        
   materialpoint_heat

 implicit none
#include <finclude/petscdmda.h90>
#include <finclude/petscsnes.h90>

!--------------------------------------------------------------------------------------------------
! input data for solution
 integer(pInt), intent(in) :: nActivePhaseFields
 type(phaseFieldDataBin), intent(in) :: phaseFieldData(nActivePhaseFields)
 real(pReal), intent(in) :: &
   timeinc, &                                                                                       !< increment in time for current solution
   timeinc_old, &                                                                                   !< increment in time of last increment
   loadCaseTime, &                                                                                  !< remaining time of current load case
   temperature_bc, &
   density
 logical, intent(in) :: &
   guess
 type(tBoundaryCondition),      intent(in) :: &
   P_BC, &
   F_BC
 character(len=*), intent(in) :: &
   incInfoIn
 real(pReal), dimension(3,3), intent(in) :: rotation_BC
 integer(pInt) :: i
 
!--------------------------------------------------------------------------------------------------
! PETSc Data
 PetscScalar, pointer :: xx_psc(:,:,:,:), F(:,:,:,:)
 PetscErrorCode :: ierr   
 SNESConvergedReason :: reason
 incInfo = incInfoIn

 call DMDAVecGetArrayF90(da,solution_vec,xx_psc,ierr); CHKERRQ(ierr)                                     ! get the data out of PETSc to work with
 F => xx_psc(0:8,:,:,:)
!--------------------------------------------------------------------------------------------------
! restart information for spectral solver
 if (restartWrite) then
   write(6,'(/,a)') ' writing converged results for restart'
   flush(6)
   call IO_write_jobRealFile(777,'F',size(F))                                                     ! writing deformation gradient field to file
   write (777,rec=1) F
   close (777)
   call IO_write_jobRealFile(777,'F_lastInc',size(F_lastInc))                                     ! writing F_lastInc field to file
   write (777,rec=1) F_lastInc
   close (777)
   call IO_write_jobRealFile(777,'F_lastInc2',size(F_lastInc2))                                   ! writing F_lastInc field to file
   write (777,rec=1) F_lastInc2
   close (777)
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
 mesh_ipCoordinates = reshape(mesh_deformedCoordsFFT(geomSize,reshape(&
                                             F,[3,3,grid(1),grid(2),grid(3)])),[3,1,product(grid)])
 if (cutBack) then 
   F_aim = F_aim_lastInc
   xx_psc(0:8,:,:,:) = reshape(F_lastInc,[9,grid(1),grid(2),grid(3)]) 
   C_volAvg = C_volAvgLastInc
   do i = 1, nActivePhaseFields
     xx_psc(8+i,:,:,:) = phaseField_lastInc(i,:,:,:)
   enddo
 else
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
   mesh_ipCoordinates = reshape(mesh_deformedCoordsFFT(geomSize,reshape(&
                                            F,[3,3,grid(1),grid(2),grid(3)])),[3,1,product(grid)])
   Fdot =  Utilities_calculateRate(math_rotate_backward33(f_aimDot,params%rotation_BC), &
                                 timeinc_old,guess,F_lastInc,reshape(F,[3,3,grid(1),grid(2),grid(3)]))
   do i = 1, nActivePhaseFields
     phaseFieldDot(i,:,:,:) = (xx_psc(8+i,:,:,:) - phaseField_lastInc(i,:,:,:))/timeinc_old
     phaseField_lastInc(i,:,:,:) = xx_psc(8+i,:,:,:)
     phaseFieldRHS_lastInc(i,:,:,:) = phaseFieldRHS(i,:,:,:)
   enddo
   F_lastInc2 = F_lastInc
   F_lastInc = reshape(F,[3,3,grid(1),grid(2),grid(3)])  
 endif
 F_aim = F_aim + f_aimDot * timeinc

 F = reshape(Utilities_forwardField(timeinc,F_lastInc,Fdot,math_rotate_backward33(F_aim, &
                                                            rotation_BC)),[9,grid(1),grid(2),grid(3)])
 do i = 1, nActivePhaseFields
   xx_psc(8+i,:,:,:) = phaseField_lastInc(i,:,:,:) + phaseFieldDot(i,:,:,:)*timeinc
 enddo
 call DMDAVecRestoreArrayF90(da,solution_vec,xx_psc,ierr); CHKERRQ(ierr)
  
!--------------------------------------------------------------------------------------------------
! update stiffness (and gamma operator)
 S = Utilities_maskedCompliance(rotation_BC,P_BC%maskLogical,C_volAvg)
 if (update_gamma) call Utilities_updateGamma(C_minmaxAvg,restartWrite)
 
 ForwardData = .True.

!--------------------------------------------------------------------------------------------------
! set module wide availabe data 
 mask_stress = P_BC%maskFloat
 params%P_BC = P_BC%values
 params%rotation_BC = rotation_BC
 params%timeinc = timeinc
 params%timeincOld = timeinc_old
 params%temperature = temperature_BC
 params%density = density
 params%nActivePhaseFields = nActivePhaseFields
 params%phaseFieldData(1:nActivePhaseFields) = phaseFieldData(1:nActivePhaseFields)

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
 use math, only: &
   math_rotate_backward33, &
   math_transpose33, &
   math_mul3333xx33
 use debug, only: &
   debug_level, &
   debug_spectral, &
   debug_spectralRotation
 use DAMASK_spectral_Utilities, only: &
   grid, &
   geomSize, &
   wgt, &
   field_real, &
   field_fourier, &
   phaseField_real, &
   phaseField_fourier, &
   Utilities_FFTforward, &
   Utilities_FFTbackward, &
   utilities_scalarFFTforward, &
   utilities_scalarFFTbackward, &
   Utilities_fourierConvolution, &
   Utilities_inverseLaplace, &
   Utilities_diffusion, &
   Utilities_constitutiveResponse, &
   Utilities_divergenceRMS
 use IO, only: &
   IO_intOut 
 use crystallite, only: &
   crystallite_temperature
 use homogenization, only: &                        
   materialpoint_heat, &
   materialpoint_P

 implicit none
 DMDALocalInfo, dimension(DMDA_LOCAL_INFO_SIZE) :: &
   in
 PetscScalar, target, dimension(9+params%nActivePhaseFields, &
   XG_RANGE,YG_RANGE,ZG_RANGE) :: &
   x_scal
 PetscScalar, target, dimension(9+params%nActivePhaseFields, &
   X_RANGE,Y_RANGE,Z_RANGE) :: &
   f_scal
 PetscScalar, pointer, dimension(:,:,:,:) :: &
   F, &
   residual_F
 PetscInt :: &
   PETScIter, &
   nfuncs
 PetscObject :: dummy
 PetscErrorCode :: ierr
 integer(pInt) :: i

 F                   => x_scal(1:9,1:grid(1),1:grid(2),1:grid(3))
 residual_F          => f_scal(1:9,1:grid(1),1:grid(2),1:grid(3))
 
 call SNESGetNumberFunctionEvals(snes,nfuncs,ierr); CHKERRQ(ierr)
 call SNESGetIterationNumber(snes,PETScIter,ierr); CHKERRQ(ierr)

 if(nfuncs== 0 .and. PETScIter == 0) totalIter = -1_pInt                                            ! new increment
 if (totalIter <= PETScIter) then                                                                   ! new iteration
!--------------------------------------------------------------------------------------------------
! report begin of new iteration
   totalIter = totalIter + 1_pInt
   write(6,'(1x,a,3(a,'//IO_intOut(itmax)//'))') trim(incInfo), &
                    ' @ Iteration ', itmin, '≤',totalIter, '≤', itmax
   if (iand(debug_level(debug_spectral),debug_spectralRotation) /= 0) &
   write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') ' deformation gradient aim (lab) =', &
                                 math_transpose33(math_rotate_backward33(F_aim,params%rotation_BC))
   write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') ' deformation gradient aim =', &
                                 math_transpose33(F_aim)
   flush(6)
 endif

!--------------------------------------------------------------------------------------------------
! evaluate inertia
 dynamic: if (params%density > 0.0_pReal) then
   residual_F = ((F - reshape(F_lastInc,[9,grid(1),grid(2),grid(3)]))/params%timeinc - &
                 reshape(F_lastInc - F_lastInc2, [9,grid(1),grid(2),grid(3)])/params%timeincOld)/&
                ((params%timeinc + params%timeincOld)/2.0_pReal)
   residual_F = params%density*product(geomSize/grid)*residual_F
   field_real = 0.0_pReal
   field_real(1:grid(1),1:grid(2),1:grid(3),1:3,1:3) = reshape(residual_F,[grid(1),grid(2),grid(3),3,3],&
                                                               order=[4,5,1,2,3])                   ! field real has a different order
   call Utilities_FFTforward()
   call Utilities_inverseLaplace()
   inertiaField_fourier = field_fourier
 else dynamic
   inertiaField_fourier = cmplx(0.0_pReal,0.0_pReal,pReal)
 endif dynamic 

!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
 do i = 1, params%nActivePhaseFields
   if(params%phaseFieldData(i)%label == 'thermal') &
     crystallite_temperature(1,1_pInt:product(grid)) = &
       reshape(x_scal(9+i,1:grid(1),1:grid(2),1:grid(3)),[product(grid)])
 enddo
 
 call Utilities_constitutiveResponse(F_lastInc,F,params%temperature,params%timeinc, &
                                     residual_F,C_volAvg,C_minmaxAvg,P_av,ForwardData,params%rotation_BC)
 ForwardData = .false.
 
 do i = 1, params%nActivePhaseFields
   if(params%phaseFieldData(i)%label == 'fracture') &
     residual_F = residual_F * spread(x_scal(9+i,1:grid(1),1:grid(2),1:grid(3)),dim=1,ncopies=9)
 enddo
  
!--------------------------------------------------------------------------------------------------
! stress BC handling
 F_aim_lastIter = F_aim
 F_aim = F_aim - math_mul3333xx33(S, ((P_av - params%P_BC)))                                        ! S = 0.0 for no bc
 err_stress = maxval(abs(mask_stress * (P_av - params%P_BC)))                                       ! mask = 0.0 for no bc
 
!--------------------------------------------------------------------------------------------------
! updated deformation gradient using fix point algorithm of basic scheme
 field_real = 0.0_pReal
 field_real(1:grid(1),1:grid(2),1:grid(3),1:3,1:3) = reshape(residual_F,[grid(1),grid(2),grid(3),3,3],&
                                                             order=[4,5,1,2,3]) ! field real has a different order
 call Utilities_FFTforward()
 field_fourier = field_fourier + inertiaField_fourier
 err_divDummy = Utilities_divergenceRMS()
 call Utilities_fourierConvolution(math_rotate_backward33(F_aim_lastIter-F_aim,params%rotation_BC))
 call Utilities_FFTbackward()
 
!--------------------------------------------------------------------------------------------------
! constructing phase field residual 
 do i = 1, params%nActivePhaseFields
   select case (params%phaseFieldData(i)%label)
     case ('thermal')
       phaseField_real = 0.0_pReal
       phaseField_real(1:grid(1),1:grid(2),1:grid(3)) = &
          phaseField_lastInc(i,1:grid(1),1:grid(2),1:grid(3))
       call utilities_scalarFFTforward()
       call utilities_diffusion(params%phaseFieldData(i)%diffusion,params%timeinc)
       call utilities_scalarFFTbackward()
       f_scal(9+i,1:grid(1),1:grid(2),1:grid(3)) = &
           phaseField_real(1:grid(1),1:grid(2),1:grid(3))
       
       phaseFieldRHS(i,1:grid(1),1:grid(2),1:grid(3)) = &
           reshape(materialpoint_heat(1,1_pInt:product(grid)),[grid(1),grid(2),grid(3)])
       phaseField_real = 0.0_pReal
       phaseField_real(1:grid(1),1:grid(2),1:grid(3)) = &
         params%timeinc*params%phaseFieldData(i)%mobility* &
           (phaseFieldRHS_lastInc(i,1:grid(1),1:grid(2),1:grid(3)) + &
            phaseFieldRHS        (i,1:grid(1),1:grid(2),1:grid(3)))/2.0_pReal
       call utilities_scalarFFTforward()
       call utilities_diffusion(params%phaseFieldData(i)%diffusion,params%timeinc/2.0_pReal)
       call utilities_scalarFFTbackward()
       f_scal(9+i,1:grid(1),1:grid(2),1:grid(3)) = &
           x_scal(9+i,1:grid(1),1:grid(2),1:grid(3)) - &
           f_scal(9+i,1:grid(1),1:grid(2),1:grid(3)) - &
           phaseField_real(1:grid(1),1:grid(2),1:grid(3))
       err_phaseField(i) = maxval(abs(f_scal(9+i,1:grid(1),1:grid(2),1:grid(3)))) 
       phaseField_Avg(i) = sum(x_scal(9+i,1:grid(1),1:grid(2),1:grid(3)))*wgt
 
     case ('fracture')
       
       phaseField_real = 0.0_pReal
       phaseField_real(1:grid(1),1:grid(2),1:grid(3)) = &
          phaseField_lastInc(i,1:grid(1),1:grid(2),1:grid(3))
       call utilities_scalarFFTforward()
       call utilities_diffusion(2.0_pReal*maxval(geomSize/real(grid,pReal))* &
                                params%phaseFieldData(i)%diffusion,params%timeinc)
       call utilities_scalarFFTbackward()
       f_scal(9+i,1:grid(1),1:grid(2),1:grid(3)) = &
           phaseField_real(1:grid(1),1:grid(2),1:grid(3))
       
       phaseFieldRHS(i,1:grid(1),1:grid(2),1:grid(3)) = &
          - params%phaseFieldData(i)%mobility* &
            sum(residual_F* &
                (F-reshape(spread(spread(spread(math_I3,3,grid(1)),4,grid(2)),5,grid(3)),[9,grid(1),grid(2),grid(3)])),dim=1) &
          - params%phaseFieldData(i)%diffusion*(x_scal(9+i,1:grid(1),1:grid(2),1:grid(3)) - 1.0_pReal)/ &
            8.0_pReal/maxval(geomSize/real(grid,pReal))
       phaseField_real = 0.0_pReal
       phaseField_real(1:grid(1),1:grid(2),1:grid(3)) = &
         params%timeinc*params%phaseFieldData(i)%mobility* &
           (phaseFieldRHS_lastInc(i,1:grid(1),1:grid(2),1:grid(3)) + &
            phaseFieldRHS        (i,1:grid(1),1:grid(2),1:grid(3)))/2.0_pReal
       call utilities_scalarFFTforward()
       call utilities_diffusion(2.0_pReal*maxval(geomSize/real(grid,pReal))* &
                                params%phaseFieldData(i)%diffusion,params%timeinc/2.0_pReal)
       call utilities_scalarFFTbackward()
       f_scal(9+i,1:grid(1),1:grid(2),1:grid(3)) = &
           x_scal(9+i,1:grid(1),1:grid(2),1:grid(3)) - &
           f_scal(9+i,1:grid(1),1:grid(2),1:grid(3)) - &
           phaseField_real(1:grid(1),1:grid(2),1:grid(3))
       err_phaseField(i) = maxval(abs(f_scal(9+i,1:grid(1),1:grid(2),1:grid(3)))) 
       phaseField_Avg(i) = sum(x_scal(9+i,1:grid(1),1:grid(2),1:grid(3)))*wgt
     
   end select  
 enddo  
 
!--------------------------------------------------------------------------------------------------
! constructing residual
 residual_F = reshape(field_real(1:grid(1),1:grid(2),1:grid(3),1:3,1:3),&
                                      [9,grid(1),grid(2),grid(3)],order=[2,3,4,1])
                                      
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
   stressTol, &
   phaseField_err = 0.0_pReal 
 
 divTol    = max(maxval(abs(P_av))*err_div_tolRel,err_div_tolAbs)
 stressTol = max(maxval(abs(P_av))*err_stress_tolrel,err_stress_tolabs)
 err_divPrev = err_div; err_div = err_divDummy
 
 if (params%nActivePhaseFields .ne. 0_pInt) phaseField_err = maxval(err_phaseField/phaseField_Avg)
 converged: if ((totalIter >= itmin .and. &
                 all([ err_div/divTol, err_stress/stressTol] < 1.0_pReal) .and. &
                 phaseField_err < 1.0e-3_pReal) &
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
 write(6,'(a,f10.2,a,es8.2,a,es9.2,a)')   ' error phase field =  ', &
            maxval(err_phaseField/phaseField_Avg)/1.0e-3, ' (',maxval(err_phaseField/phaseField_Avg), ' Pa,  tol =',1.0e-3,')' 
 write(6,'(/,a)') ' ==========================================================================='
 flush(6) 
 
end subroutine BasicPETSc_converged


!--------------------------------------------------------------------------------------------------
!> @brief convergence check
!--------------------------------------------------------------------------------------------------
subroutine BasicPETSc_convergedKSP(ksp_local,PETScIter,fnorm,reason,dummy,ierr)
 use numerics, only: &
   itmax, &
   itmin, &
   err_div_tolRel, &
   err_div_tolAbs
 use FEsolving, only: &
   terminallyIll
 use DAMASK_spectral_Utilities, only: &
   wgt
 
 implicit none
 KSP :: ksp_local
 PetscInt :: PETScIter, SNESIter
 PetscReal :: &
   fnorm, &
   SNESfnorm, &
   estimatedErrDiv
 KSPConvergedReason :: reason
 PetscObject :: dummy
 PetscErrorCode :: ierr
 real(pReal) :: &
   divTol, &
   r_tol 

 call SNESGetIterationNumber(snes,SNESIter,ierr); CHKERRQ(ierr)
 call SNESGetFunctionNorm(snes,SNESfnorm,ierr); CHKERRQ(ierr)
 
 if (SNESIter == 0_pInt) then                                                              ! Eisenstat-Walker calculation of relative tolerance for inexact newton 
   r_tol = 0.3
 else
   r_tol = (err_div/err_divPrev)**1.618
 endif     
 
 divTol    = max(maxval(abs(P_av))*err_div_tolRel,err_div_tolAbs)
 estimatedErrDiv = fnorm*err_div/SNESfnorm                                                 ! Estimated error divergence
 
 converged: if ((PETScIter >= itmin .and. &
                 any([fnorm/snesFnorm/r_tol, &
                     estimatedErrDiv/divTol] < 1.0_pReal)) &
                .or. terminallyIll) then   
   reason = 1
 elseif (totalIter >= itmax) then converged
   reason = -1
 else converged
   reason = 0
 endif converged
 
end subroutine BasicPETSc_convergedKSP


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
 call PetscFinalize(ierr); CHKERRQ(ierr)
 call Utilities_destroy()

end subroutine BasicPETSc_destroy

end module DAMASK_spectral_SolverBasicPETSc
