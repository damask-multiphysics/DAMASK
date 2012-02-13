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

MODULE DAMASK_interface
 use prec, only: pInt, pReal
 implicit none

 character(len=64), parameter :: FEsolver = 'Spectral'
 character(len=5),  parameter :: InputFileExtension = '.geom'
 character(len=4),  parameter :: LogFileExtension = '.log'    !until now, we don't have a log file. But IO.f90 requires it
 character(len=1024) :: geometryParameter,loadcaseParameter
CONTAINS

!********************************************************************
! initialize interface module
!
!********************************************************************
subroutine DAMASK_interface_init()
 use, intrinsic :: iso_fortran_env                                          ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 implicit none

 character(len=1024) commandLine, hostName, userName
 integer(pInt):: i, start = 0_pInt, length=0_pInt
 integer, dimension(8) ::  date_and_time_values                    ! type default integer
 call get_command(commandLine)
 call DATE_AND_TIME(VALUES=date_and_time_values)
 do i=1,len(commandLine)                                           ! remove capitals
   if(64<iachar(commandLine(i:i)) .and. iachar(commandLine(i:i))<91) commandLine(i:i) =achar(iachar(commandLine(i:i))+32)
 enddo

 if(index(commandLine,' -h ',.true.)>0_pInt .or. index(commandLine,' --help ',.true.)>0_pInt) then   ! search for ' -h ' or '--help'
   write(6,*) '$Id$'
#include "compilation_info.f90"
   print '(a)',  '#############################################################'
   print '(a)',  'DAMASK spectral:'
   print '(a)',  'The spectral method boundary value problem solver for'
   print '(a)',  'the Duesseldorf Advanced Material Simulation Kit'
   print '(a)',  '#############################################################'
   print '(a)',  'Valid command line switches:'
   print '(a)',  '   --geom    (-g, --geometry)'
   print '(a)',  '   --load    (-l, --loadcase)'
   print '(a)',  '   --restart (-r)'
   print '(a)',  '   --help    (-h)'
   print '(a)',  ' '
   print '(a)',  'Mandatory Arguments:'
   print '(a)',  '  --load PathToLoadFile/NameOfLoadFile.load'
   print '(a)',  '       "PathToGeomFile" will be the working directory.'
   print '(a)',  '       Make sure the file "material.config" exists in the working'
   print '(a)',  '           directory'   
   print '(a)',  '       For further configuration place "numerics.config"'
   print '(a)',  '           and "numerics.config" in that directory.'
   print '(a)',  ' '
   print '(a)',  '  --geom PathToGeomFile/NameOfGeom.geom'
   print '(a)',  ' '
   print '(a)',  'Optional Argument:'
   print '(a)',  '  --restart XX'
   print '(a)',  '       Reads in total increment No. XX-1 and continous to'
   print '(a)',  '           calculate total increment No. XX.'
   print '(a)',  '       Attention: Overwrites existing results file '
   print '(a)',  '           "NameOfGeom_NameOfLoadFile_spectralOut".'
   print '(a)',  '       Works only if the restart information for total increment'
   print '(a)',  '            No. XX-1 is available in the working directory.'
   print '(a)',  'Help:'
   print '(a)',  '  --help'
   print '(a)',  '       Prints this message and exits'
   print '(a)',  ' '
   stop
 endif
 if (.not.(command_argument_count()==4 .or. command_argument_count()==6)) & ! check for correct number of given arguments (no --help)
   stop 'Wrong Nr. of Arguments. Run DAMASK_spectral.exe --help'            ! Could not find valid keyword (position 0 +3). Functions from IO.f90 are not available
 start = index(commandLine,'-g',.true.) + 3_pInt                            ! search for '-g' and jump to first char of geometry
 if (index(commandLine,'--geom',.true.)>0) then                             ! if '--geom' is found, use that (contains '-g')
   start = index(commandLine,'--geom',.true.) + 7_pInt
 endif               
 if (index(commandLine,'--geometry',.true.)>0) then                ! again, now searching for --geometry'
   start = index(commandLine,'--geometry',.true.) + 11_pInt
 endif
 if(start==3_pInt) stop 'No Geometry specified, terminating DAMASK'! Could not find valid keyword (position 0 +3). Functions from IO.f90 are not available
 length = index(commandLine(start:len(commandLine)),' ',.false.)

 call get_command(commandLine)                                     ! may contain capitals
 geometryParameter = ''                                            ! should be empty
 geometryParameter(1:length)=commandLine(start:start+length)
 
 do i=1,len(commandLine)                                           ! remove capitals
   if(64<iachar(commandLine(i:i)) .and. iachar(commandLine(i:i))<91) commandLine(i:i) =achar(iachar(commandLine(i:i))+32)
 enddo
 
 start = index(commandLine,'-l',.true.) + 3_pInt                   ! search for '-l' and jump forward iby 3 to given name
 if (index(commandLine,'--load',.true.)>0) then                    ! if '--load' is found, use that (contains '-l')
   start = index(commandLine,'--load',.true.) + 7_pInt
 endif               
 if (index(commandLine,'--loadcase',.true.)>0) then                ! again, now searching for --loadcase'
   start = index(commandLine,'--loadcase',.true.) + 11_pInt
 endif
 if(start==3_pInt) stop 'No Loadcase specified, terminating DAMASK'! Could not find valid keyword (position 0 +3). Functions from IO.f90 are not available
 length = index(commandLine(start:len(commandLine)),' ',.false.)
 
 call get_command(commandLine)                                     ! may contain capitals
 loadcaseParameter = ''                                            ! should be empty
 loadcaseParameter(1:length)=commandLine(start:start+length)
 
 do i=1,len(commandLine)                                           ! remove capitals
   if(64<iachar(commandLine(i:i)) .and. iachar(commandLine(i:i))<91) commandLine(i:i) =achar(iachar(commandLine(i:i))+32)
 enddo

 start = index(commandLine,'-r',.true.) + 3_pInt                   ! search for '-r' and jump forward iby 3 to given name
 if (index(commandLine,'--restart',.true.)>0) then                 ! if '--restart' is found, use that (contains '-l')
   start = index(commandLine,'--restart',.true.) + 7_pInt
 endif 
 length = index(commandLine(start:len(commandLine)),' ',.false.)

 call get_command(commandLine)                                     ! may contain capitals
 call GET_ENVIRONMENT_VARIABLE('HOST',hostName)
 call GET_ENVIRONMENT_VARIABLE('USER',userName)

 write(6,*)
 write(6,*) '<<<+-  DAMASK_spectral_interface init  -+>>>'
 write(6,*) '$Id$'
#include "compilation_info.f90"
 write(6,'(a,2(i2.2,a),i4.4)') ' Date:               ',date_and_time_values(3),'/',&
                                                       date_and_time_values(2),'/',&
                                                       date_and_time_values(1) 
 write(6,'(a,2(i2.2,a),i2.2)') ' Time:               ',date_and_time_values(5),':',&
                                                       date_and_time_values(6),':',&
                                                       date_and_time_values(7)  
 write(6,*) 'Host Name:          ', trim(hostName)
 write(6,*) 'User Name:          ', trim(userName)
 write(6,*) 'Command line call:  ', trim(commandLine)
 write(6,*) 'Geometry Parameter: ', trim(geometryParameter)
 write(6,*) 'Loadcase Parameter: ', trim(loadcaseParameter)
 if (start/=3_pInt) write(6,*) 'Restart Parameter:  ', trim(commandLine(start:start+length))

endsubroutine DAMASK_interface_init

!********************************************************************
! extract working directory from loadcase file
! possibly based on current working dir
!********************************************************************
function getSolverWorkingDirectoryName()

 use prec, only: pInt
 implicit none

 character(len=1024) cwd,getSolverWorkingDirectoryName
 character(len=*), parameter :: pathSep = achar(47) //achar(92)              !forwardslash, backwardslash

 if (scan(geometryParameter,pathSep) == 1) then                              ! absolute path given as command line argument
   getSolverWorkingDirectoryName = geometryParameter(1:scan(geometryParameter,pathSep,back=.true.))
 else
   call getcwd(cwd)
   getSolverWorkingDirectoryName = trim(cwd)//'/'//geometryParameter(1:scan(geometryParameter,pathSep,back=.true.))
 endif

 getSolverWorkingDirectoryName = rectifyPath(getSolverWorkingDirectoryName)
 
endfunction getSolverWorkingDirectoryName

!********************************************************************
! basename of geometry file from command line arguments
!
!********************************************************************
function getSolverJobName()

 implicit none

 character(1024) :: getSolverJobName

 getSolverJobName = trim(getModelName())//'_'//trim(getLoadCase())

endfunction getSolverJobName

!********************************************************************
! basename of geometry file from command line arguments
!
!********************************************************************
function getModelName()

 use prec, only: pInt

 implicit none

 character(1024) getModelName, cwd
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! forwardslash, backwardslash
 integer(pInt) :: posExt,posSep
 
 posExt = scan(geometryParameter,'.',back=.true.)
 posSep = scan(geometryParameter,pathSep,back=.true.)

 if (posExt <= posSep) posExt = len_trim(geometryParameter)+1_pInt     ! no extension present
 getModelName = geometryParameter(1:posExt-1_pInt)                   ! path to geometry file (excl. extension)

 if (scan(getModelName,pathSep) /= 1) then                ! relative path given as command line argument
   call getcwd(cwd)
   getModelName = rectifyPath(trim(cwd)//'/'//getModelName)
 else
   getModelName = rectifyPath(getModelName)
 endif

 getModelName = makeRelativePath(getSolverWorkingDirectoryName(),&
                                 getModelName)
                                 
endfunction getModelName

!********************************************************************
! name of load case file exluding extension
!
!********************************************************************
function getLoadCase()

 use prec, only: pInt

 implicit none

 character(1024) getLoadCase
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! forwardslash, backwardslash
 integer(pInt) posExt,posSep

 posExt = scan(loadcaseParameter,'.',back=.true.)
 posSep = scan(loadcaseParameter,pathSep,back=.true.)

 if (posExt <= posSep) posExt = len_trim(loadcaseParameter)+1_pInt            ! no extension present
 getLoadCase = loadcaseParameter(posSep+1_pInt:posExt-1_pInt)                 ! name of load case file exluding extension

endfunction getLoadCase


!********************************************************************
! relative path of loadcase from command line arguments
!
!********************************************************************
function getLoadcaseName()

 use prec, only: pInt

 implicit none

 character(len=1024) getLoadcaseName,cwd
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! forwardslash, backwardslash
 integer(pInt) posExt,posSep
 posExt = 0_pInt

 getLoadcaseName = loadcaseParameter
 posExt = scan(getLoadcaseName,'.',back=.true.)
 posSep = scan(getLoadcaseName,pathSep,back=.true.)

 if (posExt <= posSep) getLoadcaseName = trim(getLoadcaseName)//('.load')   ! no extension present
 if (scan(getLoadcaseName,pathSep) /= 1) then          ! relative path given as command line argument
   call getcwd(cwd)
   getLoadcaseName = rectifyPath(trim(cwd)//'/'//getLoadcaseName)
 else
   getLoadcaseName = rectifyPath(getLoadcaseName)
 endif

 getLoadcaseName = makeRelativePath(getSolverWorkingDirectoryName(),&
                                    getLoadcaseName)
endfunction getLoadcaseName


!********************************************************************
! remove ../ and ./ from path
!
!********************************************************************
function rectifyPath(path)

 use prec, only: pInt

 implicit none

 character(len=*) path
 character(len=len_trim(path)) rectifyPath
 integer(pInt) i,j,k,l

 !remove ./ from path
 l = len_trim(path)
 rectifyPath = path
 do i = l,3_pInt,-1_pInt
    if ( rectifyPath(i-1_pInt:i) == './' .and. rectifyPath(i-2_pInt:i-2_pInt) /= '.' ) &
      rectifyPath(i-1_pInt:l) = rectifyPath(i+1_pInt:l)//'  '
 enddo

 !remove ../ and corresponding directory from rectifyPath
 l = len_trim(rectifyPath)
 i = index(rectifyPath(i:l),'../')
 j = 0_pInt
 do while (i > j)
    j = scan(rectifyPath(1:i-2_pInt),'/',back=.true.)
    rectifyPath(j+1_pInt:l) = rectifyPath(i+3_pInt:l)//repeat(' ',2_pInt+i-j)
    if (rectifyPath(j+1_pInt:j+1_pInt) == '/') then !search for '//' that appear in case of XXX/../../XXX
      k = len_trim(rectifyPath)
      rectifyPath(j+1_pInt:k-1_pInt) = rectifyPath(j+2_pInt:k)
      rectifyPath(k:k) = ' '
    endif
    i = j+index(rectifyPath(j+1_pInt:l),'../')
 enddo
 if(len_trim(rectifyPath) == 0) rectifyPath = '/'

 end function rectifyPath



!********************************************************************
! relative path from absolute a to absolute b
!
!********************************************************************
function makeRelativePath(a,b)

 use prec, only: pInt

 implicit none

 character (len=*) :: a,b
 character (len=1024) :: makeRelativePath
 integer(pInt) i,posLastCommonSlash,remainingSlashes

 posLastCommonSlash = 0_pInt
 remainingSlashes = 0_pInt
 do i = 1_pInt,min(1024,len_trim(a),len_trim(b))
   if (a(i:i) /= b(i:i)) exit
   if (a(i:i) == '/') posLastCommonSlash = i
 enddo
 do i = posLastCommonSlash+1_pInt,len_trim(a)
   if (a(i:i) == '/') remainingSlashes = remainingSlashes + 1_pInt
 enddo
 makeRelativePath = repeat('../',remainingSlashes)//b(posLastCommonSlash+1_pInt:len_trim(b))

endfunction makeRelativePath

END MODULE
