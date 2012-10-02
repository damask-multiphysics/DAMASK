! Copyright 2011 Max-Planck-Institut für Eisenforschung GmbH
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
!##############################################################
!* $Id$
!********************************************************************
! Material subroutine for Abaqus
!
! written by P. Eisenlohr,
!            F. Roters,
!            K. Janssens 2,
!            A. Prakash 3
!
!   MPI fuer Eisenforschung, Duesseldorf
! 2 PSI, Switzerland
! 3 Fraunhofer IWM, Freiburg
!
! REMARK: put the included file abaqus_v6.env in either your home
!         or model directory, it is a minimum Abaqus environment file
!         containing all changes necessary to use the MPIE subroutine
!        (see Abaqus documentation for more information on the use of abaqus_v6.env)
!
!********************************************************************

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
character(len=4),  parameter :: InputFileExtension = '.inp'
character(len=4),  parameter :: LogFileExtension = '.log'

contains

!--------------------
subroutine DAMASK_interface_init()
!--------------------
  write(6,*)
  write(6,*) '<<<+-  DAMASK_abaqus init  -+>>>'
  write(6,*) '$Id$'
  write(6,*)

end subroutine DAMASK_interface_init

!--------------------
function getSolverWorkingDirectoryName()
!--------------------
 use prec

 implicit none
 character(1024) getSolverWorkingDirectoryName
 integer(pInt) LENOUTDIR

 getSolverWorkingDirectoryName=''
 CALL GETOUTDIR( getSolverWorkingDirectoryName, LENOUTDIR )
 getSolverWorkingDirectoryName=trim(getSolverWorkingDirectoryName)//'/'
! write(6,*) 'getSolverWorkingDirectoryName', getSolverWorkingDirectoryName

end function getSolverWorkingDirectoryName


!--------------------
function getSolverJobName()
!--------------------
 use prec
 
 implicit none
 character(1024) getSolverJobName, JOBNAME
 integer(pInt) LENJOBNAME

 getSolverJobName=''
 CALL GETJOBNAME(getSolverJobName , LENJOBNAME )
! write(6,*) 'getSolverJobName', getSolverJobName

end function getSolverJobName

end module DAMASK_interface

#include "IO.f90"
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
#include "homogenization_isostrain.f90"
#include "homogenization_RGC.f90"
#include "homogenization.f90"
#include "CPFEM.f90"
 
subroutine UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,&
       RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,&
       TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,&
       NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,&
       DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
 
 use prec, only:      pReal, &
                      pInt
 use FEsolving, only: cycleCounter, &
                      theInc, &
                      calcMode, &
                      lastMode, &
                      theTime, &
                      theDelta, &
                      lastIncConverged, &
                      outdatedByNewInc, &
                      outdatedFFN1, &
                      terminallyIll, &
                      symmetricSolver
 use math, only:      invnrmMandel
 use debug, only:     debug_info, &
                      debug_reset, &
                      debug_levelBasic, &
                      debug_level, &
                      debug_abaqus
 use mesh, only:      mesh_FEasCP
 use CPFEM, only:     CPFEM_general,CPFEM_init_done, CPFEM_initAll
 use homogenization, only: materialpoint_sizeResults, materialpoint_results


 implicit none
 CHARACTER(80) CMNAME
 integer(pInt) ndi, nshr, ntens, nstatv, nprops, noel, npt,&
               kslay, kspt, kstep, kinc
 real(pReal) STRESS(NTENS),STATEV(NSTATV),&
       DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),&
       STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),&
       PROPS(NPROPS),COORDS(3),DROT(3,3),&
       DFGRD0(3,3),DFGRD1(3,3)
 real(pReal) SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP,&
        DTEMP, PNEWDT, CELENT
 real(pReal), dimension (3,3) :: pstress             ! not used, but needed for call of cpfem_general
 real(pReal), dimension (3,3,3,3) :: dPdF            ! not used, but needed for call of cpfem_general

! local variables
 real(pReal), dimension(6) ::   stress_h
 real(pReal), dimension(6,6) :: ddsdde_h
 integer(pInt) computationMode, i, cp_en
 logical :: cutBack

 if (iand(debug_level(debug_abaqus),debug_levelBasic) /= 0 .and. noel == 1 .and. npt == 1) then
   !$OMP CRITICAL (write2out)
     write(6,*) 'el',noel,'ip',npt
     write(6,*) 'got kinc as',kinc
     write(6,*) 'got dStran',dstran
     call flush(6)
   !$OMP END CRITICAL (write2out)
 endif

 if (.not. CPFEM_init_done) call CPFEM_initAll(temp,noel,npt)

 cp_en = mesh_FEasCP('elem',noel)
 if (time(2) > theTime .or. kinc /= theInc) then                         ! reached convergence
    terminallyIll = .false.
    cycleCounter = -1                                                    ! first calc step increments this to cycle = 0
    if (kinc == 1) then                                                  ! >> start of analysis <<
        lastIncConverged = .false.                                       ! no Jacobian backup
        outdatedByNewInc = .false.                                       ! no aging of state
        lastMode = .false.                                               ! pretend last step was collection
        calcMode = .false.                                               ! pretend last step was collection
        !$OMP CRITICAL (write2out)
        write (6,'(i8,1x,i2,1x,a)') noel,npt,'<< UMAT >> start of analysis..!'; call flush(6)
        !$OMP END CRITICAL (write2out)
    else if (kinc - theInc > 1) then                                     ! >> restart of broken analysis <<
        lastIncConverged = .false.                                       ! no Jacobian backup
        outdatedByNewInc = .false.                                       ! no aging of state
        lastMode = .true.                                                ! pretend last step was calculation
        calcMode = .true.                                                ! pretend last step was calculation
        !$OMP CRITICAL (write2out)
        write (6,'(i8,1x,i2,1x,a)') noel,npt,'<< UMAT >> restart of analysis..!'; call flush(6)
        !$OMP END CRITICAL (write2out)
    else                                                                 ! >> just the next inc <<
        lastIncConverged = .true.                                        ! request Jacobian backup
        outdatedByNewInc = .true.                                        ! request aging of state
        lastMode = .true.                                                ! assure last step was calculation
        calcMode = .true.                                                ! assure last step was calculation
        !$OMP CRITICAL (write2out)
        write (6,'(i8,1x,i2,1x,a)') noel,npt,'<< UMAT >> new increment..!'; call flush(6)
        !$OMP END CRITICAL (write2out)
    endif
    
 else if ( dtime < theDelta ) then                                     ! >> cutBack <<

    cutBack = .true.                                                    
    terminallyIll = .false.
    cycleCounter = -1                                                   ! first calc step increments this to cycle = 0
    calcMode = .true.                                                   ! pretend last step was calculation
    !$OMP CRITICAL (write2out)
    write(6,'(i8,1x,i2,1x,a)') noel,npt,'<< UMAT >> cutback detected..!'; call flush(6)
    !$OMP END CRITICAL (write2out)

 endif                                                                  ! convergence treatment end

 calcMode(npt,cp_en) = .not. calcMode(npt,cp_en)                        ! ping pong (calc <--> collect)

 if ( calcMode(npt,cp_en) ) then                                        ! now calc
    if ( lastMode .neqv. calcMode(npt,cp_en) ) then                         ! first after ping pong
        call debug_reset()                                              ! resets debugging
        outdatedFFN1 = .false.
        cycleCounter = cycleCounter + 1
    endif
    if ( outdatedByNewInc ) then
        outdatedByNewInc = .false.
        computationMode = 1                                                ! calc and age results
    else
        computationMode = 2                                                ! plain calc
    endif
 else                                                                  ! now collect
    if ( lastMode .neqv. calcMode(npt,cp_en) .and. &
         .not. terminallyIll) then
        call debug_info()                                              ! first after ping pong reports debugging
    endif
    if ( lastIncConverged ) then
        lastIncConverged = .false.
        computationMode = 4                                            ! collect and backup Jacobian after convergence
    elseif ( cutBack ) then
        cutBack = .false.
        computationMode = 5                                            ! collect and restore Jacobian after cutback
    else
        computationMode = 3                                            ! plain collect
    endif
 endif

 theTime  = time(2)                                                    ! record current starting time
 theDelta = dtime                                                      ! record current time increment
 theInc   = kinc                                                       ! record current increment number
 lastMode = calcMode(npt,cp_en)                                        ! record calculationMode

 if (iand(debug_level(debug_abaqus),debug_levelBasic) /= 0) then
   !$OMP CRITICAL (write2out)
     write(6,'(a16,1x,i2,1x,a,i8,a,i8,1x,i5,a)') 'computationMode',computationMode,'(',cp_en,':',noel,npt,')'; call flush(6)
   !$OMP END CRITICAL (write2out)
 endif
   
 call CPFEM_general(computationMode,COORDS,dfgrd0,dfgrd1,temp,dtime,noel,npt,stress_h,ddsdde_h, pstress, dPdF)

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

 if ( terminallyIll ) pnewdt = 0.5_pReal                               ! force cutback directly ?

 end subroutine UMAT

!********************************************************************
!     This subroutine replaces the corresponding Marc subroutine
!********************************************************************
 subroutine quit(mpie_error)

 use prec, only: pInt
 
 implicit none
 integer(pInt) mpie_error
 flush(6)
 call xit
 end subroutine quit
