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
 character(len=64), parameter :: FEsolver = 'Spectral'
 character(len=5),  parameter :: InputFileExtension = '.geom'
 character(len=4),  parameter :: LogFileExtension = '.log'    !until now, we don't have a log file. But IO.f90 requires it

CONTAINS

!********************************************************************
! initialize interface module
!
!********************************************************************
subroutine DAMASK_interface_init()

 write(6,*)
 write(6,*) '<<<+-  DAMASK_spectral_interface init  -+>>>'
 write(6,*) '$Id$'
 write(6,*)

endsubroutine DAMASK_interface_init

!********************************************************************
! extract working directory from loadcase file
! possibly based on current working dir
!********************************************************************
function getSolverWorkingDirectoryName()

 use prec, only: pInt
 implicit none

 character(len=1024) cwd,commandLine,outName,getSolverWorkingDirectoryName
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! forwardslash, backwardslash
 integer(pInt):: i, start, length
 
 call get_command(commandLine)
 do i=1,len(commandLine)                                           ! remove capitals
   if(64<iachar(commandLine(i:i)) .and. iachar(commandLine(i:i))<91) commandLine(i:i) =achar(iachar(commandLine(i:i))+32)
 enddo
 
 start = index(commandLine,'-g',.true.) + 3_pInt                   ! search for '-g' and jump to first char of geometry
 if (index(commandLine,'--geom',.true.)>0) then                    ! if '--geom' is found, use that (contains '-g')
   start = index(commandLine,'--geom',.true.) + 7_pInt
 endif               
 if (index(commandLine,'--geometry',.true.)>0) then                ! again, now searching for --geometry'
   start = index(commandLine,'--geometry',.true.) + 11_pInt
 endif
 if(start==3_pInt) stop 'No Geometry Specified, terminating DAMASK'! Could not find valid keyword functions from IO.f90 are not available
 length = index(commandLine(start:len(commandLine)),' ',.false.)

 call get_command(commandLine)                                     ! may contain capitals
 outName = ' '                                                     ! should be empty
 outName(1:length)=commandLine(start:start+length)
 
 if (scan(outName,pathSep) == 1) then                              ! absolute path given as command line argument
   getSolverWorkingDirectoryName = outName(1:scan(outName,pathSep,back=.true.))
 else
   call getcwd(cwd)
   getSolverWorkingDirectoryName = trim(cwd)//'/'//outName(1:scan(outName,pathSep,back=.true.))
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

 character(1024) getModelName, outName, cwd, commandLine
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! forwardslash, backwardslash
 integer(pInt) :: i,posExt,posSep,start,length
 
 call get_command(commandLine)
 do i=1,len(commandLine)                                           ! remove capitals
   if(64<iachar(commandLine(i:i)) .and. iachar(commandLine(i:i))<91) commandLine(i:i) =achar(iachar(commandLine(i:i))+32)
 enddo
 
 start = index(commandLine,'-g',.true.) + 3_pInt                   ! search for '-g' and jump to first char of geometry
 if (index(commandLine,'--geom',.true.)>0) then                    ! if '--geom' is found, use that (contains '-g')
   start = index(commandLine,'--geom',.true.) + 7_pInt
 endif               
 if (index(commandLine,'--geometry',.true.)>0) then                ! again, now searching for --geometry'
   start = index(commandLine,'--geometry',.true.) + 11_pInt
 endif
 if(start==3_pInt) stop 'No Geometry Specified, terminating DAMASK'! Could not find valid keyword functions from IO.f90 are not available
 length = index(commandLine(start:len(commandLine)),' ',.false.)
 
 call get_command(commandLine)                                     ! may contain capitals
 getModelName = ' '
 outName = ' '                                                     ! should be empty
 outName(1:length)=commandLine(start:start+length)
 posExt = scan(outName,'.',back=.true.)
 posSep = scan(outName,pathSep,back=.true.)

 if (posExt <= posSep) posExt = len_trim(outName)+1       ! no extension present
 getModelName = outName(1:posExt-1)                       ! path to geometry file (excl. extension)

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

 character(1024) getLoadCase, outName, commandLine
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! forwardslash, backwardslash
 integer(pInt) posExt,posSep,i,start,length

 call get_command(commandLine)
 do i=1,len(commandLine)                                           ! remove capitals
   if(64<iachar(commandLine(i:i)) .and. iachar(commandLine(i:i))<91) commandLine(i:i) =achar(iachar(commandLine(i:i))+32)
 enddo
 
 start = index(commandLine,'-l',.true.) + 3_pInt                   ! search for '-l' and jump forward to given name
 if (index(commandLine,'--load',.true.)>0) then                    ! if '--load' is found, use that (contains '-l')
   start = index(commandLine,'--load',.true.) + 7_pInt
 endif               
 if (index(commandLine,'--loadcase',.true.)>0) then                ! again, now searching for --loadcase'
   start = index(commandLine,'--loadcase',.true.) + 11_pInt
 endif
 if(start==3_pInt) stop 'No Loadcase  Specified, terminating DAMASK'! Could not find valid keyword functions from IO.f90 are not available
 length = index(commandLine(start:len(commandLine)),' ',.false.)
 
 call get_command(commandLine)                                     ! may contain capitals
 getLoadCase = ''
 outName = ' '                                                     ! should be empty
 outName(1:length)=commandLine(start:start+length)
 posExt = scan(outName,'.',back=.true.)
 posSep = scan(outName,pathSep,back=.true.)

 if (posExt <= posSep) posExt = len_trim(outName)+1                ! no extension present
 getLoadCase = outName(posSep+1:posExt-1)                          ! name of load case file exluding extension

endfunction getLoadCase


!********************************************************************
! relative path of loadcase from command line arguments
!
!********************************************************************
function getLoadcaseName()

 use prec, only: pInt

 implicit none

 character(len=1024) getLoadcaseName, cwd, commandLine
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! forwardslash, backwardslash
 integer(pInt) posExt,posSep,i,start,length
 posExt = 0

 call get_command(commandLine)
 do i=1,len(commandLine)                                           ! remove capitals
   if(64<iachar(commandLine(i:i)) .and. iachar(commandLine(i:i))<91) commandLine(i:i) =achar(iachar(commandLine(i:i))+32)
 enddo
 
 start = index(commandLine,'-l',.true.) + 3_pInt                   ! search for '-l' and jump forward to given name
 if (index(commandLine,'--load',.true.)>0) then                    ! if '--load' is found, use that (contains '-l')
   start = index(commandLine,'--load',.true.) + 7_pInt
 endif               
 if (index(commandLine,'--loadcase',.true.)>0) then                ! again, now searching for --loadcase'
   start = index(commandLine,'--loadcase',.true.) + 11_pInt
 endif
 if(start==3_pInt) stop 'No Loadcase  Specified, terminating DAMASK'! Could not find valid keyword functions from IO.f90 are not available
 length = index(commandLine(start:len(commandLine)),' ',.false.)

 call get_command(commandLine)                                     ! may contain capitals
 getLoadCaseName = ' '
 getLoadCaseName(1:length)=commandLine(start:start+length)

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
 do i = l,3,-1
    if ( rectifyPath(i-1:i) == './' .and. rectifyPath(i-2:i-2) /= '.' ) &
      rectifyPath(i-1:l) = rectifyPath(i+1:l)//'  '
 enddo

 !remove ../ and corresponding directory from rectifyPath
 l = len_trim(rectifyPath)
 i = index(rectifyPath(i:l),'../')
 j = 0_pInt
 do while (i > j)
    j = scan(rectifyPath(:i-2),'/',back=.true.)
    rectifyPath(j+1:l) = rectifyPath(i+3:l)//repeat(' ',2+i-j)
    i = j+index(rectifyPath(j+1:l),'../')
 enddo
 if(len_trim(rectifyPath) == 0) rectifyPath = '/'

 endfunction rectifyPath


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

 posLastCommonSlash = 0
 remainingSlashes = 0
 do i = 1,min(1024,len_trim(a),len_trim(b))
   if (a(i:i) /= b(i:i)) exit
   if (a(i:i) == '/') posLastCommonSlash = i
 enddo
 do i = posLastCommonSlash+1,len_trim(a)
   if (a(i:i) == '/') remainingSlashes = remainingSlashes + 1
 enddo
 makeRelativePath = repeat('../',remainingSlashes)//b(posLastCommonSlash+1:len_trim(b))

endfunction makeRelativePath

END MODULE
