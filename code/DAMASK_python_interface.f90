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

 implicit none
 private
 character(len=64),   parameter, public  :: FEsolver = 'Spectral'                                   !< Keyword for spectral solver
 character(len=5),    parameter, public  :: inputFileExtension = '.geom'                            !< File extension for geometry description
 character(len=4),    parameter, public  :: logFileExtension = '.log'                               !< Dummy variable as the spectral solver has no log
 character(len=1024),            private :: geometryParameter, &                                    !< Interpretated parameter given at command line
                                            loadcaseParameter                                       !< Interpretated parameter given at command line

 public  :: getSolverWorkingDirectoryName, & !< Interpretated parameter given at command line
            getSolverJobName, &
            getLoadCase, &
            getLoadCaseName, &
            getModelName, &
            DAMASK_interface_init
 private :: rectifyPath, &
            makeRelativePath, &
            getPathSep

contains

!--------------------------------------------------------------------------------------------------
!> @brief initializes the solver by interpreting the command line arguments. Also writes
!! information on computation on screen
!--------------------------------------------------------------------------------------------------

subroutine DAMASK_interface_init(loadcaseParameterIn,geometryParameterIn)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)

 implicit none
 character(len=1024)   :: hostName, &                                                               !< name of computer
                          userName                                                                  !< name of user calling the executable
 character(len=1024), intent(in)   :: loadcaseParameterIn                                        
 character(len=1024), intent(in)   :: geometryParameterIn                                                                    
 integer, dimension(8) :: dateAndTime                                                               ! type default integer
 
 geometryParameter = loadcaseParameterIn
 loadcaseParameter = geometryParameterIn
 
 call date_and_time(values = dateAndTime)

 call GET_ENVIRONMENT_VARIABLE('HOST',hostName)
 call GET_ENVIRONMENT_VARIABLE('USER',userName)

 write(6,*)
 write(6,*) '<<<+-  DAMASK_spectral_interface init  -+>>>'
 write(6,*) '$Id$'
#include "compilation_info.f90"
 write(6,'(a,2(i2.2,a),i4.4)') ' Date:               ',dateAndTime(3),'/',&
                                                       dateAndTime(2),'/',&
                                                       dateAndTime(1) 
 write(6,'(a,2(i2.2,a),i2.2)') ' Time:               ',dateAndTime(5),':',&
                                                       dateAndTime(6),':',&
                                                       dateAndTime(7)  
 write(6,*) 'Host Name:          ', trim(hostName)
 write(6,*) 'User Name:          ', trim(userName)
 write(6,*) 'Path Separator:     ', getPathSep()
 write(6,*) 'Geometry Parameter: ', trim(geometryParameter)
 write(6,*) 'Loadcase Parameter: ', trim(loadcaseParameter)

end subroutine DAMASK_interface_init

!--------------------------------------------------------------------------------------------------
!> @brief  extract working directory from loadcase file possibly based on current working dir
!--------------------------------------------------------------------------------------------------
 character(len=1024) function getSolverWorkingDirectoryName()

 implicit none
 character(len=1024) :: cwd
 character           :: pathSep

 pathSep = getPathSep()

 if (geometryParameter(1:1) == pathSep) then                                                        ! absolute path given as command line argument
   getSolverWorkingDirectoryName = geometryParameter(1:scan(geometryParameter,pathSep,back=.true.))
 else
   call getcwd(cwd)
   getSolverWorkingDirectoryName = trim(cwd)//pathSep//geometryParameter(1:scan(geometryParameter,pathSep,back=.true.))
 endif

 getSolverWorkingDirectoryName = rectifyPath(getSolverWorkingDirectoryName)
 
end function getSolverWorkingDirectoryName


!--------------------------------------------------------------------------------------------------
!> @brief  basename of geometry file from command line arguments
!--------------------------------------------------------------------------------------------------
character(len=1024) function getSolverJobName()

 implicit none
 getSolverJobName = trim(getModelName())//'_'//trim(getLoadCase())

end function getSolverJobName


!--------------------------------------------------------------------------------------------------
!> @brief  basename of geometry file from command line arguments
!--------------------------------------------------------------------------------------------------
character(len=1024) function getModelName()

 use prec, only: pInt

 implicit none
 character(len=1024) :: cwd
 integer :: posExt,posSep
 character :: pathSep

 pathSep = getPathSep()
 posExt = scan(geometryParameter,'.',back=.true.)
 posSep = scan(geometryParameter,pathSep,back=.true.)

 if (posExt <= posSep) posExt = len_trim(geometryParameter)+1                                       ! no extension present
 getModelName = geometryParameter(1:posExt-1_pInt)                                                  ! path to geometry file (excl. extension)

 if (scan(getModelName,pathSep) /= 1) then                                                          ! relative path given as command line argument
   call getcwd(cwd)
   getModelName = rectifyPath(trim(cwd)//'/'//getModelName)
 else
   getModelName = rectifyPath(getModelName)
 endif

 getModelName = makeRelativePath(getSolverWorkingDirectoryName(),&
                                 getModelName)
                                 
end function getModelName


!--------------------------------------------------------------------------------------------------
!> @brief  name of load case file exluding extension
!--------------------------------------------------------------------------------------------------
character(len=1024) function getLoadCase()

 implicit none
 integer :: posExt,posSep
 character :: pathSep

 pathSep = getPathSep()
 posExt = scan(loadcaseParameter,'.',back=.true.)
 posSep = scan(loadcaseParameter,pathSep,back=.true.)

 if (posExt <= posSep) posExt = len_trim(loadcaseParameter)+1                                       ! no extension present
 getLoadCase = loadcaseParameter(posSep+1:posExt-1)                                                 ! name of load case file exluding extension

end function getLoadCase


!--------------------------------------------------------------------------------------------------
!> @brief relative path of loadcase from command line arguments
!--------------------------------------------------------------------------------------------------
character(len=1024) function getLoadcaseName()

 implicit none
 character(len=1024) :: cwd
 integer :: posExt = 0, posSep
 character :: pathSep

 pathSep = getPathSep()
 getLoadcaseName = loadcaseParameter
 posExt = scan(getLoadcaseName,'.',back=.true.)
 posSep = scan(getLoadcaseName,pathSep,back=.true.)

 if (posExt <= posSep) getLoadcaseName = trim(getLoadcaseName)//('.load')                           ! no extension present
 if (scan(getLoadcaseName,pathSep) /= 1) then                                                       ! relative path given as command line argument
   call getcwd(cwd)
   getLoadcaseName = rectifyPath(trim(cwd)//pathSep//getLoadcaseName)
 else
   getLoadcaseName = rectifyPath(getLoadcaseName)
 endif

 getLoadcaseName = makeRelativePath(getSolverWorkingDirectoryName(),&
                                    getLoadcaseName)
end function getLoadcaseName


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

 use prec, only: pInt

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

end function

end module
