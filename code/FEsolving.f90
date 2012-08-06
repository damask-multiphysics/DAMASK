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
!##############################################################
module FEsolving
!##############################################################
 use prec, only: &
   pInt, &
   pReal
 
 implicit none
 private
 integer(pInt), public :: &
   cycleCounter =  0_pInt, &
   theInc       = -1_pInt, &
   restartInc   =  1_pInt
   
 real(pReal), public :: &
   theTime      = 0.0_pReal, &
   theDelta     = 0.0_pReal
   
 logical, public :: & 
   outdatedFFN1      = .false., &
   symmetricSolver   = .false., &
   restartWrite      = .false., &
   restartRead       = .false., &
   terminallyIll     = .false., &
   parallelExecution = .true.,  & 
   lastMode          = .true.,  &
   lastIncConverged  = .false., &
   outdatedByNewInc  = .false., &
   cutBack           = .false.

 integer(pInt), dimension(:,:), allocatable, public :: &
   FEsolving_execIP
   
 integer(pInt), dimension(2), public :: &
   FEsolving_execElem
   
 character(len=1024), public :: &
   modelName
   
 logical, dimension(:,:), allocatable, public :: &
   calcMode

 public :: FE_init

contains

!***********************************************************
! determine whether a symmetric solver is used
! and whether restart is requested
!***********************************************************
subroutine FE_init
 
 use, intrinsic :: iso_fortran_env                                ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
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
   IO_warning

 use DAMASK_interface
 
 implicit none
 integer(pInt), parameter :: &
   fileunit = 222_pInt, &
   maxNchunks = 6_pInt

#ifndef Spectral
 integer(pInt) :: j
 character(len=64)   :: tag
 character(len=1024) :: line
 integer(pInt), dimension(1_pInt+2_pInt*maxNchunks) :: positions
#endif
!$OMP CRITICAL (write2out)
 write(6,*)
 write(6,*) '<<<+-  FEsolving init  -+>>>'
 write(6,*) '$Id$'
#include "compilation_info.f90"
!$OMP END CRITICAL (write2out)

 modelName = getSolverJobName()
#ifdef Spectral
 restartInc = spectralRestartInc
 if(restartInc <= 0_pInt) then
   call IO_warning(warning_ID=34_pInt)
   restartInc = 1_pInt
 endif
 restartRead = restartInc > 1_pInt                           ! only read in if "true" restart requested
#else
 call IO_open_inputFile(fileunit,modelName)
 rewind(fileunit)
 do
   read (fileunit,'(a1024)',END=100) line
   positions = IO_stringPos(line,maxNchunks)
   tag = IO_lc(IO_stringValue(line,positions,1_pInt))        ! extract key
   select case(tag)
     case ('solver')
       read (fileunit,'(a1024)',END=100) line  ! next line
       positions = IO_stringPos(line,maxNchunks)
       symmetricSolver = (IO_intValue(line,positions,2_pInt) /= 1_pInt)
     case ('restart')
       read (fileunit,'(a1024)',END=100) line  ! next line
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
#ifdef Marc
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
 ! the following array are allocated by mesh.f90 and need to be deallocated in case of regridding
 if (allocated(calcMode)) deallocate(calcMode)
 if (allocated(FEsolving_execIP)) deallocate(FEsolving_execIP)
#endif
!$OMP CRITICAL (write2out)
 if (iand(debug_level(debug_FEsolving),debug_levelBasic) /= 0_pInt) then
   write(6,*) 'restart writing:    ', restartWrite
   write(6,*) 'restart reading:    ', restartRead
   if (restartRead) write(6,*) 'restart Job:        ', trim(modelName)
   write(6,*)
 endif
!$OMP END CRITICAL (write2out)

end subroutine FE_init

end module FEsolving
