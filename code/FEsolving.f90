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
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief triggering reading in of restart information when doing a restart
!> @todo Descriptions for public variables needed
!--------------------------------------------------------------------------------------------------
module FEsolving
 use prec, only: &
   pInt, &
   pReal
 
 implicit none
 private
 integer(pInt), public :: &                                                                         !< needs description
   cycleCounter =  0_pInt, &                                                                        !< needs description
   theInc       = -1_pInt, &                                                                        !< needs description
   restartInc   =  1_pInt, &                                                                        !< needs description
   lastLovl     =  0_pInt, &                                                                        !< lovl in previous call to marc hypela2
   lastStep     =  0_pInt                                                                           !< kstep in previous call to abaqus umat
   
 real(pReal), public :: &
   theTime      = 0.0_pReal, &                                                                      !< needs description
   theDelta     = 0.0_pReal                                                                         !< needs description

 logical, public :: & 
   outdatedFFN1      = .false., &                                                                   !< needs description
   symmetricSolver   = .false., &                                                                   !< use a symmetric solver (FEM)
   restartWrite      = .false., &                                                                   !< write current state to enable restart
   restartRead       = .false., &                                                                   !< restart information to continue calculation from saved state
   terminallyIll     = .false., &                                                                   !< at least one material point is terminally ill
   lastIncConverged  = .false., &                                                                   !< needs description
   outdatedByNewInc  = .false.                                                                      !< needs description

 integer(pInt), dimension(:,:), allocatable, public :: &
   FEsolving_execIP                                                                                 !< needs description
   
 integer(pInt), dimension(2), public :: &
   FEsolving_execElem                                                                               !< needs description
   
 character(len=1024), public :: &
   modelName                                                                                        !< needs description
   
 logical, dimension(:,:), allocatable, public :: &
   calcMode                                                                                         !< needs description

 public :: FE_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief determine whether a symmetric solver is used and whether restart is requested
!> @details restart information is found in input file in case of FEM solvers, in case of spectal
!> solver the information is provided by the interface module
!--------------------------------------------------------------------------------------------------
subroutine FE_init
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_level, &
   debug_FEsolving, &
   debug_levelBasic
 use IO, only: &
   IO_stringPos, &
   IO_stringValue, &
   IO_intValue, &
   IO_lc, &
#ifndef Spectral
   IO_open_inputFile, &
   IO_open_logFile, &
#endif
   IO_warning, &
   IO_timeStamp
 use DAMASK_interface
 
 implicit none
#ifndef Spectral
 integer(pInt), parameter :: &
   fileunit = 222_pInt, &
   maxNchunks = 6_pInt
 integer(pInt) :: j
 character(len=64)   :: tag
 character(len=1024) :: line
 integer(pInt), dimension(1_pInt+2_pInt*maxNchunks) :: positions
#endif

 write(6,'(/,a)') ' <<<+-  FEsolving init  -+>>>'
 write(6,'(a)')   ' $Id$'
 write(6,'(a16,a)')   ' Current time : ',IO_timeStamp()
#include "compilation_info.f90"

 modelName = getSolverJobName()
#ifdef Spectral
 restartInc = spectralRestartInc
 if(restartInc <= 0_pInt) then
   call IO_warning(warning_ID=34_pInt)
   restartInc = 1_pInt
 endif
 restartRead = restartInc > 1_pInt                                                                  ! only read in if "true" restart requested
#else
 call IO_open_inputFile(fileunit,modelName)
 rewind(fileunit)
 do
   read (fileunit,'(a1024)',END=100) line
   positions = IO_stringPos(line,maxNchunks)
   tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                               ! extract key
   select case(tag)
     case ('solver')
       read (fileunit,'(a1024)',END=100) line                                                       ! next line
       positions = IO_stringPos(line,maxNchunks)
       symmetricSolver = (IO_intValue(line,positions,2_pInt) /= 1_pInt)
     case ('restart')
       read (fileunit,'(a1024)',END=100) line                                                       ! next line
       positions = IO_stringPos(line,maxNchunks)
       restartWrite = iand(IO_intValue(line,positions,1_pInt),1_pInt) > 0_pInt
       restartRead  = iand(IO_intValue(line,positions,1_pInt),2_pInt) > 0_pInt
     case ('*restart')
       do j=2_pInt,positions(1)
         restartWrite = (IO_lc(IO_StringValue(line,positions,j)) == 'write') .or. restartWrite
         restartRead  = (IO_lc(IO_StringValue(line,positions,j)) == 'read')  .or. restartRead
       enddo
       if(restartWrite) then
         do j=2_pInt,positions(1)
           restartWrite = (IO_lc(IO_StringValue(line,positions,j)) /= 'frequency=0') .and. restartWrite
         enddo
       endif
   end select
 enddo
 100 close(fileunit)

 if (restartRead) then
#ifdef Marc4DAMASK
   call IO_open_logFile(fileunit)
   rewind(fileunit)
   do
     read (fileunit,'(a1024)',END=200) line
     positions = IO_stringPos(line,maxNchunks)
     if ( IO_lc(IO_stringValue(line,positions,1_pInt)) == 'restart' .and. &
          IO_lc(IO_stringValue(line,positions,2_pInt)) == 'file' .and. &
          IO_lc(IO_stringValue(line,positions,3_pInt)) == 'job' .and. &
          IO_lc(IO_stringValue(line,positions,4_pInt)) == 'id' ) &
        modelName = IO_StringValue(line,positions,6_pInt)
   enddo
#else
   call IO_open_inputFile(fileunit,modelName)
   rewind(fileunit)
   do
     read (fileunit,'(a1024)',END=200) line
     positions = IO_stringPos(line,maxNchunks)
     if ( IO_lc(IO_stringValue(line,positions,1_pInt))=='*heading') then
       read (fileunit,'(a1024)',END=200) line
       positions = IO_stringPos(line,maxNchunks)
       modelName = IO_StringValue(line,positions,1_pInt)
     endif
   enddo
#endif
 200 close(fileunit)
 endif

!--------------------------------------------------------------------------------------------------
! the following array are allocated by mesh.f90 and need to be deallocated in case of regridding
 if (allocated(calcMode)) deallocate(calcMode)
 if (allocated(FEsolving_execIP)) deallocate(FEsolving_execIP)
#endif
 if (iand(debug_level(debug_FEsolving),debug_levelBasic) /= 0_pInt) then
   write(6,*) 'restart writing:    ', restartWrite
   write(6,*) 'restart reading:    ', restartRead
   if (restartRead) write(6,'(a,/)') 'restart Job:        '//trim(modelName)
 endif

end subroutine FE_init

end module FEsolving
