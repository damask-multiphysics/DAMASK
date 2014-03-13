! Copyright 2011-13 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced Material Simulation Kit.
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
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Koen Janssens, Paul Scherrer Institut
!> @author Arun Prakash, Fraunhofer IWM
!> @brief interfaces DAMASK with Abaqus/Standard
!> @details put the included file abaqus_v6.env in either your home or model directory, 
!> it is a minimum Abaqus environment file  containing all changes necessary to use the 
!> DAMASK subroutine (see Abaqus documentation for more information on the use of abaqus_v6.env)
!--------------------------------------------------------------------------------------------------

#ifndef INT
#define INT 4
#endif

#ifndef FLOAT
#define FLOAT 8
#endif

#define Abaqus

#include "prec.f90"

module DAMASK_interface

implicit none
character(len=4), dimension(2),  parameter :: INPUTFILEEXTENSION = ['.pes','.inp']
character(len=4),                parameter :: LOGFILEEXTENSION   =  '.log'

contains

!--------------------------------------------------------------------------------------------------
!> @brief just reporting 
!--------------------------------------------------------------------------------------------------
subroutine DAMASK_interface_init()

 write(6,'(/,a)') ' <<<+-  DAMASK_abaqus init  -+>>>'
 write(6,'(a)')   ' $Id$'
#include "compilation_info.f90"  

end subroutine DAMASK_interface_init


!--------------------------------------------------------------------------------------------------
!> @brief using Abaqus/Standard function to get working directory name
!--------------------------------------------------------------------------------------------------
character(1024) function getSolverWorkingDirectoryName()

 implicit none
 integer :: lenOutDir

 getSolverWorkingDirectoryName=''
 call getoutdir(getSolverWorkingDirectoryName, lenOutDir)
 getSolverWorkingDirectoryName=trim(getSolverWorkingDirectoryName)//'/'
 
end function getSolverWorkingDirectoryName


!--------------------------------------------------------------------------------------------------
!> @brief using Abaqus/Standard function to get solver job name
!--------------------------------------------------------------------------------------------------
character(1024) function getSolverJobName()
 
 implicit none
 integer :: lenJobName

 getSolverJobName=''
 call getJobName(getSolverJobName, lenJobName)

end function getSolverJobName

end module DAMASK_interface

#include "IO.f90"
#include "libs.f90"
#include "numerics.f90"
#include "debug.f90"
#include "math.f90"
#include "FEsolving.f90"
#include "mesh.f90"
#include "material.f90"
#include "lattice.f90"
#include "constitutive_none.f90"
#include "constitutive_j2.f90"
#include "constitutive_phenopowerlaw.f90"
#include "constitutive_titanmod.f90"
#include "constitutive_dislotwin.f90"
#include "constitutive_nonlocal.f90"
#include "constitutive.f90"
#include "crystallite.f90"
#include "homogenization_none.f90"
#include "homogenization_isostrain.f90"
#include "homogenization_RGC.f90"
#include "homogenization.f90"
#include "CPFEM.f90"
 
subroutine UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,&
                RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,&
                TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,&
                NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,&
                DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
 
 use prec, only: &
   pReal, &
   pInt
 use numerics, only: &
   usePingPong
 use FEsolving, only: &
   cycleCounter, &
   theInc, &
   calcMode, &
   theTime, &
   theDelta, &
   lastIncConverged, &
   outdatedByNewInc, &
   outdatedFFN1, &
   terminallyIll, &
   symmetricSolver, &
   lastStep
 use math, only: &
   invnrmMandel
 use debug, only: &
   debug_info, &
   debug_reset, &
   debug_levelBasic, &
   debug_level, &
   debug_abaqus
 use mesh, only:  &
   mesh_unitlength, &
   mesh_FEasCP, &
   mesh_ipCoordinates
 use CPFEM, only: &
   CPFEM_general, &
   CPFEM_init_done, &
   CPFEM_initAll, &
   CPFEM_CALCRESULTS, &
   CPFEM_AGERESULTS, &
   CPFEM_COLLECT, &
   CPFEM_RESTOREJACOBIAN, &
   CPFEM_BACKUPJACOBIAN
 use homogenization, only: &
   materialpoint_sizeResults, &
   materialpoint_results

 implicit none
 integer(pInt),                       intent(in) :: &
   nDi, &                                                                                           !< Number of direct stress components at this point
   nShr, &                                                                                          !< Number of engineering shear stress components at this point
   nTens, &                                                                                         !< Size of the stress or strain component array (NDI + NSHR)
   nStatV, &                                                                                        !< Number of solution-dependent state variables
   nProps, &                                                                                        !< User-defined number of material constants
   noEl, &                                                                                          !< element number
   nPt,&                                                                                            !< integration point number
   kSlay, &                                                                                         !< layer number (shell elements etc.)
   kSpt, &                                                                                          !< section point within the current layer
   kStep, &                                                                                         !< step number
   kInc                                                                                             !< increment number
 character(len=80),                   intent(in) :: &
   cmname                                                                                           !< uses-specified material name, left justified
 real(pReal),                         intent(in) :: &
   DTIME, &
   TEMP, &
   DTEMP, &
   CELENT
 real(pReal), dimension(1),           intent(in) :: & 
   PREDEF, & 
   DPRED
 real(pReal), dimension(2),           intent(in) :: &
   TIME                                                                                             !< step time/total time at beginning of the current increment
 real(pReal), dimension(3),           intent(in) :: &
   COORDS
 real(pReal), dimension(nTens),       intent(in) :: &
   STRAN, &                                                                                         !< total strains at beginning of the increment
   DSTRAN                                                                                           !< strain increments
 real(pReal), dimension(nProps),      intent(in) :: &
   PROPS
 real(pReal), dimension(3,3),         intent(in) :: &
   DROT, &                                                                                          !< rotation increment matrix
   DFGRD0, &                                                                                        !< F at beginning of increment
   DFGRD1                                                                                           !< F at end of increment
 real(pReal),                         intent(inout) :: &                                                             
   PNEWDT, &                                                                                        !< ratio of suggested new time increment
   SSE, &                                                                                           !< specific elastic strain engergy
   SPD, &                                                                                           !< specific plastic dissipation
   SCD, &                                                                                           !< specific creep dissipation
   RPL, &                                                                                           !< volumetric heat generation per unit time at the end of the increment 
   DRPLDT                                                                                           !< varation of RPL with respect to the temperature
 real(pReal), dimension(nTens),       intent(inout) :: &
   STRESS                                                                                           !< stress tensor at the beginning of the increment, needs to be updated
 real(pReal), dimension(nStatV),      intent(inout) :: &
   STATEV                                                                                           !< solution-dependent state variables
 real(pReal), dimension(nTens),       intent(out) :: &
   DDSDDT, &
   DRPLDE
 real(pReal), dimension(nTens,nTens), intent(out) :: &
   DDSDDE                                                                                           !< Jacobian matrix of the constitutive model

 real(pReal) :: temperature                                                                         ! temp by Abaqus is intent(in)
 real(pReal), dimension(6) ::   stress_h
 real(pReal), dimension(6,6) :: ddsdde_h
 integer(pInt) :: computationMode, i, cp_en
 logical :: cutBack

 temperature = temp                                                                                 ! temp is intent(in)
 DDSDDT = 0.0_pReal
 DRPLDE = 0.0_pReal

 if (iand(debug_level(debug_abaqus),debug_levelBasic) /= 0 .and. noel == 1 .and. npt == 1) then
   !$OMP CRITICAL (write2out)
   write(6,*) 'el',noel,'ip',npt
   write(6,*) 'got kInc as',kInc
   write(6,*) 'got dStran',dStran
   flush(6)
   !$OMP END CRITICAL (write2out)
 endif

 if (.not. CPFEM_init_done) call CPFEM_initAll(temperature,noel,npt)

 computationMode = 0
 cp_en = mesh_FEasCP('elem',noel)
 if (time(2) > theTime .or. kInc /= theInc) then                                                    ! reached convergence
   terminallyIll = .false.
   cycleCounter = -1                                                                                ! first calc step increments this to cycle = 0
   if (kInc == 1) then                                                                              ! >> start of analysis << 
     lastIncConverged = .false.                                                                     ! no Jacobian backup
     outdatedByNewInc = .false.                                                                     ! no aging of state
     calcMode = .false.                                                                             ! pretend last step was collection
     !$OMP CRITICAL (write2out)
     write (6,'(i8,1x,i2,1x,a)') noel,npt,'<< UMAT >> start of analysis..!'
     flush(6)
     !$OMP END CRITICAL (write2out)
   else if (kInc - theInc > 1) then                                                                 ! >> restart of broken analysis <<
     lastIncConverged = .false.                                                                     ! no Jacobian backup
     outdatedByNewInc = .false.                                                                     ! no aging of state
     calcMode = .true.                                                                              ! pretend last step was calculation
     !$OMP CRITICAL (write2out)
     write (6,'(i8,1x,i2,1x,a)') noel,npt,'<< UMAT >> restart of analysis..!'
     flush(6)
     !$OMP END CRITICAL (write2out)
   else                                                                                             ! >> just the next inc << 
     lastIncConverged = .true.                                                                      ! request Jacobian backup
     outdatedByNewInc = .true.                                                                      ! request aging of state
     calcMode = .true.                                                                              ! assure last step was calculation
     !$OMP CRITICAL (write2out)
     write (6,'(i8,1x,i2,1x,a)') noel,npt,'<< UMAT >> new increment..!'
     flush(6)
     !$OMP END CRITICAL (write2out)
   endif
 else if ( dtime < theDelta ) then                                                                  ! >> cutBack <<
   lastIncConverged = .false.                                                                       ! no Jacobian backup
   outdatedByNewInc = .false.                                                                       ! no aging of state
   terminallyIll = .false.
   cycleCounter = -1                                                                                ! first calc step increments this to cycle = 0
   calcMode = .true.                                                                                ! pretend last step was calculation
   !$OMP CRITICAL (write2out)
   write(6,'(i8,1x,i2,1x,a)') noel,npt,'<< UMAT >> cutback detected..!'
   flush(6)
   !$OMP END CRITICAL (write2out)
 endif                                                                                              ! convergence treatment end


 if (usePingPong) then
   calcMode(npt,cp_en) = .not. calcMode(npt,cp_en)                                                  ! ping pong (calc <--> collect)
   if (calcMode(npt,cp_en)) then                                                                    ! now --- CALC ---                  
     computationMode = CPFEM_CALCRESULTS
     if ( lastStep /= kStep ) then                                                                  ! first after ping pong
       call debug_reset()                                                                           ! resets debugging
       outdatedFFN1 = .false.
       cycleCounter = cycleCounter + 1_pInt
     endif
     if(outdatedByNewInc) then
       computationMode = ior(computationMode,CPFEM_AGERESULTS)                                      ! calc and age results
       outdatedByNewInc = .false.                                                                   ! reset flag
     endif
   else                                                                                             ! now --- COLLECT ---
     computationMode = CPFEM_COLLECT                                                                ! plain collect
     if(lastStep /= kStep .and. .not. terminallyIll) &
       call debug_info()                                                                            ! first after ping pong reports (meaningful) debugging
     if (lastIncConverged) then
       computationMode = ior(computationMode,CPFEM_BACKUPJACOBIAN)                                  !  collect and backup Jacobian after convergence
       lastIncConverged = .false.                                                                   ! reset flag
     endif
     mesh_ipCoordinates(1:3,npt,cp_en) = mesh_unitlength * COORDS
   endif
 else                                                                                               ! --- PLAIN MODE ---
   computationMode = CPFEM_CALCRESULTS                                                              ! always calc
   if (lastStep /= kStep) then
     if (.not. terminallyIll) &
       call debug_info()                                                                            ! first reports (meaningful) debugging
     call debug_reset()                                                                             ! and resets debugging
     outdatedFFN1  = .false.
     cycleCounter  = cycleCounter + 1_pInt
   endif
   if (outdatedByNewInc) then
     computationMode = ior(computationMode,CPFEM_AGERESULTS)
     outdatedByNewInc = .false.                                                                     ! reset flag
   endif
   if (lastIncConverged) then
     computationMode = ior(computationMode,CPFEM_BACKUPJACOBIAN)                                    ! backup Jacobian after convergence
     lastIncConverged = .false.                                                                     ! reset flag
   endif
 endif
   
   
 theTime  = time(2)                                                                                 ! record current starting time
 theDelta = dtime                                                                                   ! record current time increment
 theInc   = kInc                                                                                    ! record current increment number
 lastStep = kStep                                                                                   ! record step number

 if (iand(debug_level(debug_abaqus),debug_levelBasic) /= 0) then
   !$OMP CRITICAL (write2out)
   write(6,'(a16,1x,i2,1x,a,i8,a,i8,1x,i5,a)') 'computationMode',computationMode,'(',cp_en,':',noel,npt,')'
   flush(6)
   !$OMP END CRITICAL (write2out)
 endif
   
 call CPFEM_general(computationMode,usePingPong,dfgrd0,dfgrd1,temperature,dtime,noel,npt,stress_h,ddsdde_h)

!     Mandel:              11, 22, 33, SQRT(2)*12, SQRT(2)*23, SQRT(2)*13
!     straight:            11, 22, 33, 12, 23, 13
!     ABAQUS explicit:     11, 22, 33, 12, 23, 13
!     ABAQUS implicit:     11, 22, 33, 12, 13, 23
!     ABAQUS implicit:     11, 22, 33, 12

 forall(i=1:ntens) ddsdde(1:ntens,i) = invnrmMandel(i)*ddsdde_h(1:ntens,i)*invnrmMandel(1:ntens)
 stress(1:ntens) = stress_h(1:ntens)*invnrmMandel(1:ntens)
 if(symmetricSolver) ddsdde(1:ntens,1:ntens) = 0.5_pReal*(ddsdde(1:ntens,1:ntens) + transpose(ddsdde(1:ntens,1:ntens)))
 if(ntens == 6) then
   stress_h = stress
   stress(5) = stress_h(6)
   stress(6) = stress_h(5)
   ddsdde_h = ddsdde
   ddsdde(:,5) = ddsdde_h(:,6)
   ddsdde(:,6) = ddsdde_h(:,5)
   ddsdde_h = ddsdde
   ddsdde(5,:) = ddsdde_h(6,:)
   ddsdde(6,:) = ddsdde_h(5,:)
 end if

 statev = materialpoint_results(1:min(nstatv,materialpoint_sizeResults),npt,mesh_FEasCP('elem', noel))

 if ( terminallyIll ) pnewdt = 0.5_pReal                                                            ! force cutback directly ?

end subroutine UMAT


!--------------------------------------------------------------------------------------------------
!> @brief calls the exit function of Abaqus/Standard
!--------------------------------------------------------------------------------------------------
subroutine quit(mpie_error)
 use prec, only: &
   pInt
 
 implicit none
 integer(pInt) mpie_error

 flush(6)
 call xit

end subroutine quit
