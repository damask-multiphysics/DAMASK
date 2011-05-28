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

 integer(pInt) :: cycleCounter = 0_pInt, theInc = -1_pInt
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
! determine wether a symmetric solver is used
! and whether restart is requested
!***********************************************************
 subroutine FE_init()
 
 use prec, only: pInt
 use debug, only: debug_verbosity
 use DAMASK_interface, only: getModelName, FEsolver
 use IO
 implicit none
 
 integer(pInt), parameter :: fileunit = 222
 integer(pInt), parameter :: maxNchunks = 6
 integer(pInt) i
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 character(len=64) tag
 character(len=1024) line

 FEmodelGeometry = getModelName()

 if (IO_open_inputFile(fileunit,FEmodelGeometry)) then
 
   rewind(fileunit)
   do
     read (fileunit,'(a1024)',END=100) line
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
     select case(tag)
       case ('solver')
         read (fileunit,'(a1024)',END=100) line  ! next line
         positions = IO_stringPos(line,maxNchunks)
         symmetricSolver = (IO_intValue(line,positions,2) /= 1_pInt)
       case ('restart')
         read (fileunit,'(a1024)',END=100) line  ! next line
         positions = IO_stringPos(line,maxNchunks)
         restartWrite = iand(IO_intValue(line,positions,1),1_pInt) > 0_pInt
         restartRead  = iand(IO_intValue(line,positions,1),2_pInt) > 0_pInt
       case ('*restart')
         do i=2,positions(1)
           restartWrite = IO_lc(IO_StringValue(line,positions,i)) == 'write'
!           restartRead  = IO_lc(IO_StringValue(line,positions,i)) == 'read'
         enddo
     end select
   enddo
 else
   call IO_error(101) ! cannot open input file
 endif

100 close(fileunit)
 
 if (restartRead) then
   if(FEsolver == 'Marc' .and. IO_open_logFile(fileunit)) then
     rewind(fileunit)
     do
       read (fileunit,'(a1024)',END=200) line
       positions = IO_stringPos(line,maxNchunks)
       if ( IO_lc(IO_stringValue(line,positions,1)) == 'restart' .and. &
            IO_lc(IO_stringValue(line,positions,2)) == 'file' .and. &
            IO_lc(IO_stringValue(line,positions,3)) == 'job' .and. &
            IO_lc(IO_stringValue(line,positions,4)) == 'id' ) &
          FEmodelGeometry = IO_StringValue(line,positions,6)
     enddo
   elseif (FEsolver == 'Abaqus' .and. IO_open_jobFile(fileunit, 'com')) then
     rewind(fileunit)
     do
       read (fileunit,'(a1024)',END=200) line
       positions = IO_stringPos(line,maxNchunks)
!       if ( IO_lc(IO_stringValue(line,positions,?)) == 'oldjob?') &
!          FEmodelGeometry = IO_StringValue(line,positions,?)       
     enddo
   else
     call IO_error(106) ! cannot open file for old job info
   endif
 endif

200 close(fileunit)

!$OMP CRITICAL (write2out)
 write(6,*)
 write(6,*) '<<<+-  FEsolving init  -+>>>'
 write(6,*) '$Id$'
 write(6,*)
 if (debug_verbosity > 0) then
   write(6,*) 'restart writing:    ', restartWrite
   write(6,*) 'restart reading:    ', restartRead
   if (restartRead) write(6,*) 'restart Job:        ', trim(FEmodelGeometry)
   write(6,*)
 endif
!$OMP END CRITICAL (write2out)

 end subroutine

 END MODULE FEsolving
