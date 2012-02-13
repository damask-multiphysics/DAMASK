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
 MODULE FEsolving
!##############################################################
 use prec, only: pInt,pReal
 implicit none

 integer(pInt) :: cycleCounter = 0_pInt, theInc = -1_pInt, restartInc = 1_pInt
 real(pReal)   :: theTime = 0.0_pReal, theDelta = 0.0_pReal
 logical :: lastIncConverged = .false.,outdatedByNewInc = .false.,outdatedFFN1 = .false.,terminallyIll = .false.
 logical :: symmetricSolver = .false. 
 logical :: parallelExecution = .true. 
 logical :: restartWrite = .false.
 logical :: restartRead  = .false.
 logical :: lastMode = .true., cutBack = .false.
 logical, dimension(:,:), allocatable :: calcMode
 integer(pInt), dimension(:,:), allocatable :: FEsolving_execIP
 integer(pInt), dimension(2) :: FEsolving_execElem
 character(len=1024) FEmodelGeometry

 CONTAINS

!***********************************************************
! determine whether a symmetric solver is used
! and whether restart is requested
!***********************************************************
 subroutine FE_init()
 
 use, intrinsic :: iso_fortran_env  
 use prec, only: pInt
 use debug, only: debug_verbosity
 use DAMASK_interface
 use IO
 implicit none
 
 integer(pInt), parameter :: fileunit = 222
 integer(pInt), parameter :: maxNchunks = 6
 integer(pInt):: i, start = 0_pInt, length=0_pInt
integer(pInt), dimension(1_pInt+2_pInt*maxNchunks) :: positions
 character(len=64) tag
 character(len=1024) line, commandLine

 FEmodelGeometry = getModelName()
 call IO_open_inputFile(fileunit,FEmodelGeometry)
 if (trim(FEsolver) == 'Spectral') then
   call get_command(commandLine)                                                 ! may contain uppercase
   do i=1,len(commandLine)
     if(64 < iachar(commandLine(i:i)) .and. iachar(commandLine(i:i)) < 91) &
       commandLine(i:i) = achar(iachar(commandLine(i:i))+32)                     ! make lowercase
   enddo
   if (index(commandLine,'-r ',.true.)>0) &                                      ! look for -r
     start = index(commandLine,'-r ',.true.) + 3_pInt                            ! set to position after trailing space
   if (index(commandLine,'--restart ',.true.)>0) &                               ! look for --restart
     start = index(commandLine,'--restart ',.true.) + 10_pInt                    ! set to position after trailing space

   if(start /= 0_pInt) then                                                      ! found something
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
         do i=2,positions(1)
           restartWrite = (IO_lc(IO_StringValue(line,positions,i)) == 'write') .or. restartWrite
           restartRead  = (IO_lc(IO_StringValue(line,positions,i)) == 'read')  .or. restartRead
         enddo
         if(restartWrite) then
           do i=2,positions(1)
             restartWrite = (IO_lc(IO_StringValue(line,positions,i)) /= 'frequency=0') .and. restartWrite
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
 if (debug_verbosity > 0) then
   write(6,*) 'restart writing:    ', restartWrite
   write(6,*) 'restart reading:    ', restartRead
   if (restartRead) write(6,*) 'restart Job:        ', trim(FEmodelGeometry)
   write(6,*)
 endif
!$OMP END CRITICAL (write2out)

 end subroutine

 END MODULE FEsolving
