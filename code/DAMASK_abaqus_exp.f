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
character(len=4), dimension(2),  parameter :: InputFileExtension = ['.pes','.inp']
character(len=4),                parameter :: LogFileExtension = '.log'

contains

!--------------------
subroutine DAMASK_interface_init()
!--------------------
  use IO, only: IO_timeStamp
  
  write(6,*)
  write(6,*) '<<<+-  DAMASK_abaqus init  -+>>>'
  write(6,*) '$Id$'
  write(6,'(a16,a)')   ' Current time : ',IO_timeStamp()
#include "compilation_info.f90"  
  write(6,*)

end subroutine DAMASK_interface_init

!--------------------
function getSolverWorkingDirectoryName()
!--------------------
 use prec, only: pInt

 implicit none
 character(1024) getSolverWorkingDirectoryName
 integer(pInt) LENOUTDIR

 getSolverWorkingDirectoryName=''
 CALL VGETOUTDIR( getSolverWorkingDirectoryName, LENOUTDIR )
 getSolverWorkingDirectoryName=trim(getSolverWorkingDirectoryName)//'/'
 
end function getSolverWorkingDirectoryName

!--------------------
function getSolverJobName()
!--------------------
 use prec, only: pInt
 
 implicit none
 character(1024) getSolverJobName, JOBNAME
 integer(pInt) LENJOBNAME

 getSolverJobName=''
 CALL VGETJOBNAME(getSolverJobName , LENJOBNAME )
 
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

subroutine vumat (jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, &
                  stepTime, totalTime, dt, cmname, coordMp, charLength, &
                  props, density, strainInc, relSpinInc, &
                  tempOld, stretchOld, defgradOld, fieldOld, &
                  stressOld, stateOld, enerInternOld, enerInelasOld, &
                  tempNew, stretchNew, defgradNew, fieldNew, &
                  stressNew, stateNew, enerInternNew, enerInelasNew )

 include 'vaba_param.inc'

 dimension jblock(*), props(nprops), density(*), coordMp(*), &
           charLength(*), strainInc(*), &
           relSpinInc(*), tempOld(*), &
           stretchOld(*), &
           defgradOld(*), &
           fieldOld(*), stressOld(*), &
           stateOld(*), enerInternOld(*), &
           enerInelasOld(*), tempNew(*), &
           stretchNew(*), &
           defgradNew(*), &
           fieldNew(*), &
           stressNew(*), stateNew(*), &
           enerInternNew(*), enerInelasNew(*)

 character(80) cmname


 call vumatXtrArg ( jblock(1), &
                    ndir, nshr, nstatev, nfieldv, nprops, lanneal, &
                    stepTime, totalTime, dt, cmname, coordMp, charLength, &
                    props, density, strainInc, relSpinInc, &
                    tempOld, stretchOld, defgradOld, fieldOld, &
                    stressOld, stateOld, enerInternOld, enerInelasOld, &
                    tempNew, stretchNew, defgradNew, fieldNew, &
                    stressNew, stateNew, enerInternNew, enerInelasNew, &
                    jblock(5), jblock(2))

 end subroutine vumat


 subroutine vumatXtrArg (nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, &
                         stepTime, totalTime, dt, cmname, coordMp, charLength, &
                         props, density, strainInc, relSpinInc, &
                         tempOld, stretchOld, defgradOld, fieldOld, &
                         stressOld, stateOld, enerInternOld, enerInelasOld, &
                         tempNew, stretchNew, defgradNew, fieldNew, &
                         stressNew, stateNew, enerInternNew, enerInelasNew, &
                         nElement, nMatPoint)

 use prec, only:      pReal, &
                      pInt
 use numerics, only:  numerics_unitlength
 use FEsolving, only: cycleCounter, &
                      theTime, &
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
 use mesh, only:      mesh_FEasCP, &
                      mesh_ipCoordinates
 use CPFEM, only:     CPFEM_general,CPFEM_init_done, CPFEM_initAll
 use homogenization, only: materialpoint_sizeResults, materialpoint_results

 include 'vaba_param.inc'      ! Abaqus exp initializes a first step in single prec. for this a two-step compilation is used.
                               ! symbolic linking switches between .._sp.inc and .._dp.inc for both consecutive compilations...


 dimension props(nprops), density(nblock), &
           strainInc(nblock,ndir+nshr), &
           relSpinInc(nblock,nshr), defgradOld(nblock,ndir+nshr+nshr), &
           stressOld(nblock,ndir+nshr), &
           stateOld(nblock,nstatev), enerInternOld(nblock), &
           enerInelasOld(nblock), tempNew(nblock), tempOld(nblock), &
           stretchNew(nblock,ndir+nshr), defgradNew(nblock,ndir+nshr+nshr), &
           stressNew(nblock,ndir+nshr), coordMp(nblock,3) 

 dimension enerInelasNew(nblock),stateNew(nblock,nstatev),enerInternNew(nblock)
 dimension nElement(nblock),nMatPoint(nblock)

 character(80) cmname
 real(pReal), dimension (3,3) :: pstress                ! not used, but needed for call of cpfem_general
 real(pReal), dimension (3,3,3,3) :: dPdF               ! not used, but needed for call of cpfem_general
! local variables
 real(pReal), dimension(3) :: coordinates
 real(pReal), dimension(3,3) :: defgrd0,defgrd1
 real(pReal), dimension(6) ::   stress
 real(pReal), dimension(6,6) :: ddsdde
 real(pReal) temp, timeInc
 integer(pInt) computationMode, n, i, cp_en
 logical :: cutBack

 do n = 1,nblock                                                       ! loop over vector of IPs

   temp    = tempOld(n)
   if ( .not. CPFEM_init_done ) then
     call CPFEM_initAll(temp,nElement(n),nMatPoint(n))
     outdatedByNewInc = .false.

     if (iand(debug_level(debug_abaqus),debug_levelBasic) /= 0) then
       !$OMP CRITICAL (write2out)
         write(6,'(i8,1x,i2,1x,a)') nElement(n),nMatPoint(n),'first call special case..!'; call flush(6)
       !$OMP END CRITICAL (write2out)
     endif

   else if (theTime < totalTime) then                                  ! reached convergence
     outdatedByNewInc = .true.

     if (iand(debug_level(debug_abaqus),debug_levelBasic) /= 0) then
       !$OMP CRITICAL (write2out)
         write (6,'(i8,1x,i2,1x,a)') nElement(n),nMatPoint(n),'lastIncConverged + outdated'; call flush(6)
       !$OMP END CRITICAL (write2out)
     endif

   endif

   outdatedFFN1 = .false.
   terminallyIll = .false.
   cycleCounter = 1
   if ( outdatedByNewInc ) then
     outdatedByNewInc = .false.
     call debug_info()                                             ! first after new inc reports debugging
     call debug_reset()                                            ! resets debugging
     computationMode = 8                                           ! calc and age results with implicit collection
   else
     computationMode = 9                                           ! plain calc with implicit collection
   endif

   theTime  = totalTime                                            ! record current starting time

   if (iand(debug_level(debug_abaqus),debug_levelBasic) /= 0) then
     !$OMP CRITICAL (write2out)
       write(6,'(a16,1x,i2,1x,a,i8,1x,i5,a)') 'computationMode',computationMode,'(',nElement(n),nMatPoint(n),')'; call flush(6)
     !$OMP END CRITICAL (write2out)
   endif
  
   defgrd0 = 0.0_pReal
   defgrd1 = 0.0_pReal
   timeInc = dt

  !     ABAQUS explicit:     deformation gradient as vector 11, 22, 33, 12, 23, 31, 21, 32, 13
  !     ABAQUS explicit:     deformation gradient as vector 11, 22, 33, 12, 21
  
   forall (i=1:ndir)
     defgrd0(i,i) = defgradOld(n,i)
     defgrd1(i,i) = defgradNew(n,i)
   end forall
   if (nshr == 1) then
     defgrd0(1,2) = defgradOld(n,4)
     defgrd1(1,2) = defgradNew(n,4)
     defgrd0(2,1) = defgradOld(n,5)
     defgrd1(2,1) = defgradNew(n,5)
   else
     defgrd0(1,2) = defgradOld(n,4)
     defgrd1(1,2) = defgradNew(n,4)
     defgrd0(1,3) = defgradOld(n,9)
     defgrd1(1,3) = defgradNew(n,9)
     defgrd0(2,1) = defgradOld(n,7)
     defgrd1(2,1) = defgradNew(n,7)
     defgrd0(2,3) = defgradOld(n,5)
     defgrd1(2,3) = defgradNew(n,5)
     defgrd0(3,1) = defgradOld(n,6)
     defgrd1(3,1) = defgradNew(n,6)
     defgrd0(3,2) = defgradOld(n,8)
     defgrd1(3,2) = defgradNew(n,8)
   endif

   cp_en = mesh_FEasCP('elem',nElement(n))
   mesh_ipCoordinates(1:3,n,cp_en) = numerics_unitlength * coordMp(n,1:3)

   call CPFEM_general(computationMode,defgrd0,defgrd1,temp,timeInc,cp_en,nMatPoint(n),stress,ddsdde, pstress, dPdF)
  
  !     Mandel:     11, 22, 33, SQRT(2)*12, SQRT(2)*23, SQRT(2)*13
  !     straight:   11, 22, 33, 12, 23, 13
  !     ABAQUS implicit:     11, 22, 33, 12, 13, 23
  !     ABAQUS explicit:     11, 22, 33, 12, 23, 13
  !     ABAQUS explicit:     11, 22, 33, 12

   stressNew(n,1:ndir+nshr) = stress(1:ndir+nshr)*invnrmMandel(1:ndir+nshr)
  
   stateNew(n,:) = materialpoint_results(1:min(nstatev,materialpoint_sizeResults),nMatPoint(n),mesh_FEasCP('elem', nElement(n)))
   tempNew(n) = temp
  
 enddo

 end subroutine vumatXtrArg

!********************************************************************
!     This subroutine replaces the corresponding Marc subroutine
!********************************************************************
 subroutine quit(mpie_error)

 use prec, only: pInt
 
 implicit none
 integer(pInt) mpie_error

 call xit
 flush(6)
 end subroutine quit
