!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief triggering reading in of restart information when doing a restart
!> @todo Descriptions for public variables needed
!--------------------------------------------------------------------------------------------------
module FEsolving
 use prec
  
 implicit none
 private
 integer, public :: &
   restartInc   =  1                                                                                !< needs description

 logical, public :: & 
   symmetricSolver   = .false., &                                                                   !< use a symmetric FEM solver
   restartWrite      = .false., &                                                                   !< write current state to enable restart
   restartRead       = .false., &                                                                   !< restart information to continue calculation from saved state
   terminallyIll     = .false.                                                                      !< at least one material point is terminally ill

 integer, dimension(:,:), allocatable, public :: &
   FEsolving_execIP                                                                                 !< for ping-pong scheme always range to max IP, otherwise one specific IP
   
 integer, dimension(2), public :: &
   FEsolving_execElem                                                                               !< for ping-pong scheme always whole range, otherwise one specific element
   
 character(len=1024), public :: &
   modelName                                                                                        !< needs description
   
 logical, dimension(:,:), allocatable, public :: &
   calcMode                                                                                         !< do calculation or simply collect when using ping pong scheme

 public :: FE_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief determine whether a symmetric solver is used and whether restart is requested
!> @details restart information is found in input file in case of FEM solvers, in case of spectal
!> solver the information is provided by the interface module
!--------------------------------------------------------------------------------------------------
subroutine FE_init
 use debug, only: &
   debug_level, &
   debug_FEsolving, &
   debug_levelBasic
 use IO, only: &
   IO_stringPos, &
   IO_stringValue, &
   IO_intValue, &
   IO_lc, &
#if defined(Marc4DAMASK) || defined(Abaqus)
   IO_open_inputFile, &
   IO_open_logFile, &
#endif
   IO_warning
 use DAMASK_interface
 
#if defined(Marc4DAMASK) || defined(Abaqus)
 integer, parameter :: &
   FILEUNIT = 222
 integer :: j
 character(len=65536) :: tag, line
 integer, allocatable, dimension(:) :: chunkPos
#endif

 write(6,'(/,a)')   ' <<<+-  FEsolving init  -+>>>'

 modelName = getSolverJobName()

#if defined(Grid) || defined(FEM)
 restartInc = interface_RestartInc

 if(restartInc < 0) then
   call IO_warning(warning_ID=34)
   restartInc = 0
 endif
 restartRead = restartInc > 0                                                                       ! only read in if "true" restart requested
#else
 call IO_open_inputFile(FILEUNIT,modelName)
 rewind(FILEUNIT)
 do
   read (FILEUNIT,'(a1024)',END=100) line
   chunkPos = IO_stringPos(line)
   tag = IO_lc(IO_stringValue(line,chunkPos,1))                                                     ! extract key
   select case(tag)
     case ('solver')
       read (FILEUNIT,'(a1024)',END=100) line                                                       ! next line
       chunkPos = IO_stringPos(line)
       symmetricSolver = (IO_intValue(line,chunkPos,2) /= 1)
     case ('restart')
       read (FILEUNIT,'(a1024)',END=100) line                                                       ! next line
       chunkPos = IO_stringPos(line)
       restartWrite = iand(IO_intValue(line,chunkPos,1),1) > 0
       restartRead  = iand(IO_intValue(line,chunkPos,1),2) > 0
     case ('*restart')
       do j=2,chunkPos(1)
         restartWrite = (IO_lc(IO_StringValue(line,chunkPos,j)) == 'write') .or. restartWrite
         restartRead  = (IO_lc(IO_StringValue(line,chunkPos,j)) == 'read')  .or. restartRead
       enddo
       if(restartWrite) then
         do j=2,chunkPos(1)
           restartWrite = (IO_lc(IO_StringValue(line,chunkPos,j)) /= 'frequency=0') .and. restartWrite
         enddo
       endif
   end select
 enddo
 100 close(FILEUNIT)

 if (restartRead) then
#ifdef Marc4DAMASK
   call IO_open_logFile(FILEUNIT)
   rewind(FILEUNIT)
   do
     read (FILEUNIT,'(a1024)',END=200) line
     chunkPos = IO_stringPos(line)
     if (   IO_lc(IO_stringValue(line,chunkPos,1)) == 'restart' &
      .and. IO_lc(IO_stringValue(line,chunkPos,2)) == 'file'    &
      .and. IO_lc(IO_stringValue(line,chunkPos,3)) == 'job'     &
      .and. IO_lc(IO_stringValue(line,chunkPos,4)) == 'id' )    &
        modelName = IO_StringValue(line,chunkPos,6)
   enddo
#else                                                                                               ! QUESTION: is this meaningful for the spectral/FEM case?
   call IO_open_inputFile(FILEUNIT,modelName)
   rewind(FILEUNIT)
   do
     read (FILEUNIT,'(a1024)',END=200) line
     chunkPos = IO_stringPos(line)
     if (IO_lc(IO_stringValue(line,chunkPos,1))=='*heading') then
       read (FILEUNIT,'(a1024)',END=200) line
       chunkPos = IO_stringPos(line)
       modelName = IO_StringValue(line,chunkPos,1)
     endif
   enddo
#endif
 200 close(FILEUNIT)
 endif

#endif
 if (iand(debug_level(debug_FEsolving),debug_levelBasic) /= 0) then
   write(6,'(a21,l1)') ' restart writing:    ', restartWrite
   write(6,'(a21,l1)') ' restart reading:    ', restartRead
   if (restartRead) write(6,'(a,/)') ' restart Job:        '//trim(modelName)
 endif

end subroutine FE_init

end module FEsolving
