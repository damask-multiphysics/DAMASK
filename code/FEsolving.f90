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
   parallelExecution = .true., & 
   lastMode          = .true.

 integer(pInt), dimension(:,:), allocatable, public :: &
   FEsolving_execIP
   
 integer(pInt), dimension(2), public :: &
   FEsolving_execElem
   
 character(len=1024), public :: &
   FEmodelGeometry
   
 logical, dimension(:,:), allocatable, public :: &
   calcMode
      
 logical, private :: & 
   lastIncConverged  = .false., &
   outdatedByNewInc  = .false., &
   cutBack           = .false.
   
 public :: FE_init

contains

!***********************************************************
! determine whether a symmetric solver is used
! and whether restart is requested
!***********************************************************
subroutine FE_init
 
 use, intrinsic :: iso_fortran_env                                ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_what, &
   debug_FEsolving, &
   debug_levelBasic
 
 use IO, only: &
   IO_open_inputFile, &
   IO_stringPos, &
   IO_stringValue, &
   IO_intValue, &
   IO_lc, &
   IO_open_logFile, &
   IO_warning

 use DAMASK_interface
 
 implicit none
 integer(pInt), parameter :: &
   fileunit = 222_pInt, &
   maxNchunks = 6_pInt

 integer :: i, start = 0, length                                 ! is save for FE_init (only called once)
 integer(pInt) :: j
 integer(pInt), dimension(1_pInt+2_pInt*maxNchunks) :: positions
 character(len=64)   :: tag
 character(len=1024) :: line, &
                        commandLine
 
 FEmodelGeometry = getModelName()
 call IO_open_inputFile(fileunit,FEmodelGeometry)

 if (trim(FEsolver) == 'Spectral') then
   call get_command(commandLine)                                                 ! may contain uppercase
   do i=1,len(commandLine)
     if(64 < iachar(commandLine(i:i)) .and. iachar(commandLine(i:i)) < 91) &
       commandLine(i:i) = achar(iachar(commandLine(i:i))+32)                     ! make lowercase
   enddo
   if (index(commandLine,'-r ',.true.)>0) &                                      ! look for -r
     start = index(commandLine,'-r ',.true.) + 3                                 ! set to position after trailing space
   if (index(commandLine,'--restart ',.true.)>0) &                               ! look for --restart
     start = index(commandLine,'--restart ',.true.) + 10                         ! set to position after trailing space
   if(start /= 0) then                                                           ! found something
     length = verify(commandLine(start:len(commandLine)),'0123456789',.false.)   ! where is first non number after argument?
     read(commandLine(start:start+length),'(I12)') restartInc                    ! read argument
     restartRead  = restartInc > 0_pInt
     if(restartInc <= 0_pInt) then
       call IO_warning(warning_ID=34_pInt)
       restartInc = 1_pInt
     endif
   endif
 else
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
 endif
 100 close(fileunit)
 
 if (restartRead) then
   if(FEsolver == 'Marc') then
     call IO_open_logFile(fileunit)
     rewind(fileunit)
     do
       read (fileunit,'(a1024)',END=200) line
       positions = IO_stringPos(line,maxNchunks)
       if ( IO_lc(IO_stringValue(line,positions,1_pInt)) == 'restart' .and. &
            IO_lc(IO_stringValue(line,positions,2_pInt)) == 'file' .and. &
            IO_lc(IO_stringValue(line,positions,3_pInt)) == 'job' .and. &
            IO_lc(IO_stringValue(line,positions,4_pInt)) == 'id' ) &
          FEmodelGeometry = IO_StringValue(line,positions,6_pInt)
     enddo
   elseif (FEsolver == 'Abaqus') then
     call IO_open_inputFile(fileunit,FEmodelGeometry)
     rewind(fileunit)
     do
       read (fileunit,'(a1024)',END=200) line
       positions = IO_stringPos(line,maxNchunks)
       if ( IO_lc(IO_stringValue(line,positions,1_pInt))=='*heading') then
         read (fileunit,'(a1024)',END=200) line
         positions = IO_stringPos(line,maxNchunks)
         FEmodelGeometry = IO_StringValue(line,positions,1_pInt)
       endif
     enddo
   endif
 endif

200 close(fileunit)

!$OMP CRITICAL (write2out)
 write(6,*)
 write(6,*) '<<<+-  FEsolving init  -+>>>'
 write(6,*) '$Id$'
#include "compilation_info.f90"
 if (iand(debug_what(debug_FEsolving),debug_levelBasic) /= 0_pInt) then
   write(6,*) 'restart writing:    ', restartWrite
   write(6,*) 'restart reading:    ', restartRead
   if (restartRead) write(6,*) 'restart Job:        ', trim(FEmodelGeometry)
   write(6,*)
 endif
!$OMP END CRITICAL (write2out)

end subroutine FE_init

end module FEsolving
