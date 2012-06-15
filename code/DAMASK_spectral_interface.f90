! Copyright 2011 Max-Planck-Institut für Eisenforschung GmbH
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
!* $Id$
!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Interfacing between the spectral solver and the material subroutines provided
!! by DAMASK
!--------------------------------------------------------------------------------------------------
module DAMASK_interface
 use prec, only: &
   pInt

 implicit none
 private
 logical,                       public  :: &
   appendToOutFile = .false.                                                                        !< Append to existing spectralOut file (in case of restart, not in case of regridding)
 integer(pInt),                 public  :: &
   spectralRestart = 1_pInt                                                                         !< Increment at which calculation starts
 character(len=1024),           public  :: &
   geometryFile = '', &                                                                             !< parameter given for geometry file
   loadCaseFile = ''                                                                                !< parameter given for load case file

 public  :: getSolverWorkingDirectoryName, &
            getSolverJobName, &
            DAMASK_interface_init
 private :: getGeometryFile, &
            getLoadCaseFile, &
            rectifyPath, &
            makeRelativePath, &
            getPathSep, &
            IO_stringValue, &
            IO_intValue, &
            IO_lc, &
            IO_stringPos

contains

!--------------------------------------------------------------------------------------------------
!> @brief initializes the solver by interpreting the command line arguments. Also writes
!! information on computation to screen
!--------------------------------------------------------------------------------------------------
subroutine DAMASK_interface_init(loadCaseParameterIn,geometryParameterIn)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)

 implicit none
 character(len=1024), optional, intent(in) :: &
   loadCaseParameterIn, &                                        
   geometryParameterIn    
 character(len=1024) :: &
   commandLine, &                                                                                   !< command line call as string
   geometryParameter, &
   loadCaseParameter, &
   hostName, &                                                                                      !< name of computer
   userName, &                                                                                         !< name of user calling the executable
   tag
 integer :: &
   i
 integer, parameter :: &
   maxNchunks = 7
 integer, dimension(1+ 2* maxNchunks) :: &
   positions
 integer, dimension(8) :: &
   dateAndTime                                                                                      ! type default integer
 logical :: &
   gotLoadCase = .false., &
   gotGeometry = .false., &
   gotRestart  = .false.
 
 write(6,'(a)') ''
 write(6,'(a)') '<<<+-  DAMASK_spectral_interface init  -+>>>'
 write(6,'(a)') '$Id$'
#include "compilation_info.f90"

 if ( present(loadcaseParameterIn) .and. present(geometryParameterIn)) then                         ! both mandatory parameters given in function call 
   geometryParameter = geometryParameterIn
   loadcaseParameter = loadcaseParameterIn
   commandLine='n/a'
   gotLoadCase = .true.
   gotGeometry = .true.
 else if ( .not.( present(loadcaseParameterIn) .and. present(geometryParameterIn))) then            ! none parameters given in function call, trying to get them from command line
   call get_command(commandLine)
   positions = IO_stringPos(commandLine,maxNchunks)
   do i = 1, maxNchunks
     tag = IO_lc(IO_stringValue(commandLine,positions,i))                                           ! extract key
     select case(tag)
       case ('-h','--help')
         write(6,'(a)')  '#############################################################'
         write(6,'(a)')  'DAMASK spectral:'
         write(6,'(a)')  'The spectral method boundary value problem solver for'
         write(6,'(a)')  'the Duesseldorf Advanced Material Simulation Kit'
         write(6,'(a)')  '#############################################################'
         write(6,'(a)')  'Valid command line switches:'
         write(6,'(a)')  '   --geom    (-g, --geometry)'
         write(6,'(a)')  '   --load    (-l, --loadcase)'
         write(6,'(a)')  '   --restart (-r, --rs)'
         write(6,'(a)')  '   --regrid  (--rg)'
         write(6,'(a)')  '   --help    (-h)'
         write(6,'(a)')  ' '
         write(6,'(a)')  'Mandatory Arguments:'
         write(6,'(a)')  '  --load PathToLoadFile/NameOfLoadFile.load'
         write(6,'(a)')  '       "PathToLoadFile" will be the working directory.'
         write(6,'(a)')  '       Make sure the file "material.config" exists in the working'
         write(6,'(a)')  '           directory'   
         write(6,'(a)')  '       For further configuration place "numerics.config"'
         write(6,'(a)')  '           and "numerics.config" in that directory.'
         write(6,'(a)')  ' '
         write(6,'(a)')  '  --geom PathToGeomFile/NameOfGeom.geom'
         write(6,'(a)')  ' '
         write(6,'(a)')  'Optional Argument:'
         write(6,'(a)')  '  --restart XX'
         write(6,'(a)')  '       Reads in total increment No. XX-1 and continous to'
         write(6,'(a)')  '           calculate total increment No. XX.'
         write(6,'(a)')  '       Appends to existing results file '
         write(6,'(a)')  '           "NameOfGeom_NameOfLoadFile.spectralOut".'
         write(6,'(a)')  '       Works only if the restart information for total increment'
         write(6,'(a)')  '            No. XX-1 is available in the working directory.'
         write(6,'(a)')  ' '
         write(6,'(a)')  '  --regrid XX'
         write(6,'(a)')  '       Reads in total increment No. XX-1 and continous to'
         write(6,'(a)')  '           calculate total increment No. XX.'
         write(6,'(a)')  '       Attention: Overwrites existing results file '
         write(6,'(a)')  '           "NameOfGeom_NameOfLoadFile.spectralOut".'
         write(6,'(a)')  '       Works only if the restart information for total increment'
         write(6,'(a)')  '            No. XX-1 is available in the working directory.'
         write(6,'(a)')  'Help:'
         write(6,'(a)')  '  --help'
         write(6,'(a)')  '       Prints this message and exits'
         write(6,'(a)')  ' '
         call quit(0_pInt)                                                                                ! normal Termination
       case ('-l', '--load', '--loadcase')
         if (gotLoadCase) then
           write(6,'(a)') 'Got 2nd time loadcase keyword'
           call quit(1_pInt)
         endif
         loadcaseParameter = IO_stringValue(commandLine,positions,i+1_pInt)
         gotLoadCase = .true.
       case ('-g', '--geom', '--geometry')
         if (gotGeometry) then
           write(6,'(a)') 'Got 2nd time geometry keyword'
           call quit(1_pInt)
          endif
         geometryParameter = IO_stringValue(commandLine,positions,i+1_pInt)
         gotGeometry = .true.
       case ('-r', '--rs', '--restart')
         spectralRestart = IO_IntValue(commandLine,positions,i+1_pInt)
         appendToOutFile = .true.
         if (gotRestart) then
           write(6,'(a)') 'Got 2nd time restart/regrid keyword'
           call quit(1_pInt)
         endif
       case ('--rg', '--regrid')
         spectralRestart = IO_IntValue(commandLine,positions,i+1_pInt)
         if (gotRestart) then
           write(6,'(a)') 'Got 2nd time restart/regrid keyword'
           call quit(1_pInt)
         endif
     end select
   enddo
 endif
 
 if (.not. (gotLoadCase .and. gotGeometry)) then
   write(6,'(a)') 'Please specify Geometry AND Load Case'
   call quit(1_pInt)
 endif

 geometryFile = getGeometryFile(geometryParameter)
 loadCaseFile = getLoadCaseFile(loadCaseParameter)

 call get_environment_variable('HOST',hostName)
 call get_environment_variable('USER',userName)
 call date_and_time(values = dateAndTime)

 write(6,'(a,2(i2.2,a),i4.4)') 'Date:                  ',dateAndTime(3),'/',&
                                                         dateAndTime(2),'/',&
                                                         dateAndTime(1) 
 write(6,'(a,2(i2.2,a),i2.2)') 'Time:                  ',dateAndTime(5),':',&
                                                         dateAndTime(6),':',&
                                                         dateAndTime(7)  
 write(6,'(a,a)')  'Host name:             ', trim(hostName)
 write(6,'(a,a)')  'User name:             ', trim(userName)
 write(6,'(a,a)')  'Path separator:        ', getPathSep()
 write(6,'(a,a)')  'Command line call:     ', trim(commandLine)
 write(6,'(a,a)')  'Geometry parameter:    ', trim(geometryParameter)
 write(6,'(a,a)')  'Loadcase parameter:    ', trim(loadcaseParameter)
 if (SpectralRestart > 1) write(6,'(a,i6.6)') &
                   'Restart at increment:  ', spectralRestart 
 write(6,'(a,l1)') 'Append to result file: ', appendToOutFile

end subroutine DAMASK_interface_init

!--------------------------------------------------------------------------------------------------
!> @brief  extract working directory from loadcase file possibly based on current working dir
!--------------------------------------------------------------------------------------------------
character(len=1024) function getSolverWorkingDirectoryName()

 implicit none
 character(len=1024) :: cwd
 character           :: pathSep

 pathSep = getPathSep()

 if (geometryFile(1:1) == pathSep) then                                                             ! absolute path given as command line argument
   getSolverWorkingDirectoryName = geometryFile(1:scan(geometryFile,pathSep,back=.true.))
 else
   call getcwd(cwd)                                                                                 ! relative path given as command line argument
   getSolverWorkingDirectoryName = trim(cwd)//pathSep//&
                                   geometryFile(1:scan(geometryFile,pathSep,back=.true.))
 endif

 getSolverWorkingDirectoryName = rectifyPath(getSolverWorkingDirectoryName)
 
end function getSolverWorkingDirectoryName


!--------------------------------------------------------------------------------------------------
!> @brief solver job name (no extension) as combination of geometry and load case name
!--------------------------------------------------------------------------------------------------
character(len=1024) function getSolverJobName()

 implicit none
 integer :: posExt,posSep
 character :: pathSep
 character(len=1024) :: tempString 

 pathSep = getPathSep()

 tempString = geometryFile
 posExt = scan(tempString,'.',back=.true.)
 posSep = scan(tempString,pathSep,back=.true.)

 getSolverJobName = tempString(posSep+1:posExt-1)

 tempString = loadCaseFile
 posExt = scan(tempString,'.',back=.true.)
 posSep = scan(tempString,pathSep,back=.true.)

 getSolverJobName = trim(getSolverJobName)//'_'//tempString(posSep+1:posExt-1)

end function getSolverJobName


!--------------------------------------------------------------------------------------------------
!> @brief basename of geometry file with extension from command line arguments
!--------------------------------------------------------------------------------------------------
character(len=1024) function getGeometryFile(geometryParameter)

 implicit none
 character(len=1024), intent(in) :: &
   geometryParameter
 character(len=1024) :: &
   cwd
 integer :: posExt, posSep
 character :: pathSep

 getGeometryFile = geometryParameter
 pathSep = getPathSep()
 posExt = scan(getGeometryFile,'.',back=.true.)
 posSep = scan(getGeometryFile,pathSep,back=.true.)

 if (posExt <= posSep) getGeometryFile = trim(getGeometryFile)//('.geom')                           ! no extension present
 if (scan(getGeometryFile,pathSep) /= 1) then                                                       ! relative path given as command line argument
   call getcwd(cwd)
   getGeometryFile = rectifyPath(trim(cwd)//pathSep//getGeometryFile)
 else
   getGeometryFile = rectifyPath(getGeometryFile)
 endif

 getGeometryFile = makeRelativePath(getSolverWorkingDirectoryName(), getGeometryFile)

end function getGeometryFile


!--------------------------------------------------------------------------------------------------
!> @brief relative path of loadcase from command line arguments
!--------------------------------------------------------------------------------------------------
character(len=1024) function getLoadCaseFile(loadCaseParameter)

 implicit none
 character(len=1024), intent(in) :: &
   loadCaseParameter
 character(len=1024) :: &
   cwd
 integer :: posExt, posSep
 character :: pathSep
 
 getLoadCaseFile = loadcaseParameter
 pathSep = getPathSep()
 posExt = scan(getLoadCaseFile,'.',back=.true.)
 posSep = scan(getLoadCaseFile,pathSep,back=.true.)

 if (posExt <= posSep) getLoadCaseFile = trim(getLoadCaseFile)//('.load')                           ! no extension present
 if (scan(getLoadCaseFile,pathSep) /= 1) then                                                       ! relative path given as command line argument
   call getcwd(cwd)
   getLoadCaseFile = rectifyPath(trim(cwd)//pathSep//getLoadCaseFile)
 else
   getLoadCaseFile = rectifyPath(getLoadCaseFile)
 endif

 getLoadCaseFile = makeRelativePath(getSolverWorkingDirectoryName(), getLoadCaseFile)

end function getLoadCaseFile


!--------------------------------------------------------------------------------------------------
!> @brief remove ../ and ./ from path
!--------------------------------------------------------------------------------------------------
function rectifyPath(path)

 implicit none
 character(len=*) :: path
 character(len=len_trim(path)) :: rectifyPath
 character :: pathSep
 integer :: i,j,k,l !no pInt

 pathSep = getPathSep()

 !remove ./ from path
 l = len_trim(path)
 rectifyPath = path
 do i = l,3,-1
    if ( rectifyPath(i-1:i) == '.'//pathSep .and. rectifyPath(i-2:i-2) /= '.' ) &
      rectifyPath(i-1:l) = rectifyPath(i+1:l)//'  '
 enddo

 !remove ../ and corresponding directory from rectifyPath
 l = len_trim(rectifyPath)
 i = index(rectifyPath(i:l),'..'//pathSep)
 j = 0
 do while (i > j)
    j = scan(rectifyPath(1:i-2),pathSep,back=.true.)
    rectifyPath(j+1:l) = rectifyPath(i+3:l)//repeat(' ',2+i-j)
    if (rectifyPath(j+1:j+1) == pathSep) then                                                       !search for '//' that appear in case of XXX/../../XXX
      k = len_trim(rectifyPath)
      rectifyPath(j+1:k-1) = rectifyPath(j+2:k)
      rectifyPath(k:k) = ' '
    endif
    i = j+index(rectifyPath(j+1:l),'..'//pathSep)
 enddo
 if(len_trim(rectifyPath) == 0) rectifyPath = pathSep

end function rectifyPath

 
!--------------------------------------------------------------------------------------------------
!> @brief relative path from absolute a to absolute b
!--------------------------------------------------------------------------------------------------
character(len=1024) function makeRelativePath(a,b)

 implicit none
 character (len=*) :: a,b
 character :: pathSep
 integer :: i,posLastCommonSlash,remainingSlashes !no pInt

 pathSep = getPathSep()
 posLastCommonSlash = 0
 remainingSlashes = 0

 do i = 1, min(1024,len_trim(a),len_trim(b))
   if (a(i:i) /= b(i:i)) exit
   if (a(i:i) == pathSep) posLastCommonSlash = i
 enddo
 do i = posLastCommonSlash+1,len_trim(a)
   if (a(i:i) == pathSep) remainingSlashes = remainingSlashes + 1
 enddo
 makeRelativePath = repeat('..'//pathSep,remainingSlashes)//b(posLastCommonSlash+1:len_trim(b))

end function makeRelativePath


!--------------------------------------------------------------------------------------------------
!> @brief counting / and \ in $PATH System variable the character occuring more often is assumed
!! to be the path separator
!--------------------------------------------------------------------------------------------------
character function getPathSep()

 implicit none
 character(len=2048) path
 integer(pInt) :: backslash = 0_pInt, slash = 0_pInt
 integer :: i

 call get_environment_variable('PATH',path)
 do i=1, len(trim(path))
   if (path(i:i)=='/') slash     =     slash + 1_pInt
   if (path(i:i)=='\') backslash = backslash + 1_pInt
 enddo

 if (backslash>slash) then
   getPathSep = '\'
 else
   getPathSep = '/'
 endif

end function getPathSep

!********************************************************************
! read string value at myPos from line
!********************************************************************
 pure function IO_stringValue(line,positions,myPos)
 
 implicit none

 integer(pInt),                                intent(in) :: positions(*), &
                                                             myPos

 character(len=1+positions(myPos*2+1)-positions(myPos*2)) :: IO_stringValue

 character(len=*),                             intent(in) :: line

 if (positions(1) < myPos) then
   IO_stringValue = ''
 else
   IO_stringValue = line(positions(myPos*2):positions(myPos*2+1))
 endif

end function IO_stringValue

!********************************************************************
! read int value at myPos from line
!********************************************************************
integer(pInt) pure function IO_intValue(line,positions,myPos)

 implicit none
 character(len=*), intent(in) :: line
 integer(pInt),    intent(in) :: positions(*), &
                                 myPos

 if (positions(1) < myPos) then
   IO_intValue = 0_pInt
 else
   read(UNIT=line(positions(myPos*2):positions(myPos*2+1)),ERR=100,FMT=*) IO_intValue
 endif
 return
100 IO_intValue = huge(1_pInt)

end function IO_intValue

!********************************************************************
! change character in line to lower case
!********************************************************************
pure function IO_lc(line)

 implicit none
 character(26), parameter :: lower = 'abcdefghijklmnopqrstuvwxyz'
 character(26), parameter :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
 character(len=*), intent(in) :: line
 character(len=len(line))     :: IO_lc
 
 integer                      :: i,n                      ! no pInt (len returns default integer)

 IO_lc = line
 do i=1,len(line)
   n = index(upper,IO_lc(i:i))
   if (n/=0) IO_lc(i:i) = lower(n:n)
 enddo

end function IO_lc

!********************************************************************
! locate at most N space-separated parts in line
! return array containing number of parts in line and
! the left/right positions of at most N to be used by IO_xxxVal
!********************************************************************
pure function IO_stringPos(line,N)

 implicit none
 integer(pInt),    intent(in) :: N
 integer(pInt)                :: IO_stringPos(1_pInt+N*2_pInt)
 
 character(len=*), intent(in) :: line
 
 character(len=*), parameter  :: sep=achar(44)//achar(32)//achar(9)//achar(10)//achar(13) ! comma and whitespaces

 integer                      :: left, right                      !no pInt (verify and scan return default integer)


 IO_stringPos = -1_pInt
 IO_stringPos(1) = 0_pInt
 right = 0
 
 do while (verify(line(right+1:),sep)>0)
   left  = right + verify(line(right+1:),sep)
   right = left + scan(line(left:),sep) - 2
   if ( line(left:left) == '#' ) then
     exit
   endif
   if ( IO_stringPos(1)<N ) then
     IO_stringPos(1_pInt+IO_stringPos(1)*2_pInt+1_pInt) = int(left, pInt)
     IO_stringPos(1_pInt+IO_stringPos(1)*2_pInt+2_pInt) = int(right, pInt)
   endif
   IO_stringPos(1) = IO_stringPos(1)+1_pInt
 enddo

end function IO_stringPos

end module
