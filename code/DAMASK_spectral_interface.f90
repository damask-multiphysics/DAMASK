!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author   Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief    Interfacing between the spectral solver and the material subroutines provided
!!           by DAMASK
!> @details  Interfacing between the spectral solver and the material subroutines provided
!>           by DAMASK. Interpretating the command line arguments or, in case of called from f2py, 
!>           the arguments parsed to the init routine to get load case, geometry file, working 
!>           directory, etc.
!--------------------------------------------------------------------------------------------------
module DAMASK_interface
 use prec, only: &
   pInt

 implicit none
 private
#include <finclude/petscsys.h>

 logical,             public, protected :: appendToOutFile = .false.                                !< Append to existing spectralOut file (in case of restart, not in case of regridding)
 integer(pInt),       public, protected :: spectralRestartInc = 1_pInt                              !< Increment at which calculation starts
 character(len=1024), public, protected :: &
   geometryFile = '', &                                                                             !< parameter given for geometry file
   loadCaseFile = ''                                                                                !< parameter given for load case file
 character(len=1024), private           :: workingDirectory                                         !< accessed by getSolverWorkingDirectoryName for compatibility reasons

 public :: &
   getSolverWorkingDirectoryName, &
   getSolverJobName, &
   DAMASK_interface_init
 private :: &
   storeWorkingDirectory, &
   getGeometryFile, &
   getLoadCaseFile, &
   rectifyPath, &
   makeRelativePath, &
   getPathSep, &
   IIO_stringValue, &
   IIO_intValue, &
   IIO_lc, &
   IIO_stringPos

contains

!--------------------------------------------------------------------------------------------------
!> @brief initializes the solver by interpreting the command line arguments. Also writes
!! information on computation to screen
!--------------------------------------------------------------------------------------------------
subroutine DAMASK_interface_init(loadCaseParameterIn,geometryParameterIn)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)

 implicit none
 character(len=1024), optional, intent(in) :: &
   loadCaseParameterIn, &                                                                           !< if using the f2py variant, the -l argument of DAMASK_spectral.exe
   geometryParameterIn                                                                              !< if using the f2py variant, the -g argument of DAMASK_spectral.exe
 character(len=1024) :: &
   commandLine, &                                                                                   !< command line call as string
   loadCaseArg   ='', &                                                                             !< -l argument given to DAMASK_spectral.exe
   geometryArg   ='', &                                                                             !< -g argument given to DAMASK_spectral.exe
   workingDirArg ='', &                                                                             !< -w argument given to DAMASK_spectral.exe
   hostName, &                                                                                      !< name of machine on which DAMASK_spectral.exe is execute (might require export HOSTNAME)
   userName, &                                                                                      !< name of user calling DAMASK_spectral.exe
   tag
 integer :: &
   i, &
   worldrank = 0
 integer, parameter :: &
   maxNchunks = 128                                                                                 !< DAMASK_spectral + (l,g,w,r)*2 + h
 integer, dimension(1+ 2* maxNchunks) :: &
   positions
 integer, dimension(8) :: &
   dateAndTime                                                                                      ! type default integer
 PetscErrorCode :: ierr
 external :: &
   quit,&
   MPI_Comm_rank,&
   PETScInitialize, &
   MPI_abort

!--------------------------------------------------------------------------------------------------
! PETSc Init
 call PetscInitialize(PETSC_NULL_CHARACTER,ierr)                                                    ! according to PETSc manual, that should be the first line in the code
 CHKERRQ(ierr)                                                                                      ! this is a macro definition, it is case sensitive

 open(6, encoding='UTF-8')                                                                          ! modern fortran compilers (gfortran >4.4, ifort >11 support it)
 call MPI_Comm_rank(PETSC_COMM_WORLD,worldrank,ierr);CHKERRQ(ierr)
 mainProcess: if (worldrank == 0) then
   write(6,'(/,a)') ' <<<+-  DAMASK_spectral_interface init  -+>>>'
   write(6,'(a)')   ' $Id$'
#include "compilation_info.f90"
 endif mainProcess

 if ( present(loadcaseParameterIn) .and. present(geometryParameterIn)) then                         ! both mandatory parameters given in function call 
   geometryArg = geometryParameterIn
   loadcaseArg = loadcaseParameterIn
   commandLine = 'n/a'
 else if ( .not.( present(loadcaseParameterIn) .and. present(geometryParameterIn))) then            ! none parameters given in function call, trying to get them from command line
   call get_command(commandLine)
   positions = IIO_stringPos(commandLine,maxNchunks)
   do i = 1, positions(1)
     tag = IIO_lc(IIO_stringValue(commandLine,positions,i))                                         ! extract key
     select case(tag)
       case ('-h','--help')
         mainProcess2: if (worldrank == 0) then
           write(6,'(a)')  ' #######################################################################'
           write(6,'(a)')  ' DAMASK_spectral:'
           write(6,'(a)')  ' The spectral method boundary value problem solver for'
           write(6,'(a)')  ' the Düsseldorf Advanced Material Simulation Kit'
           write(6,'(a,/)')' #######################################################################'
           write(6,'(a,/)')' Valid command line switches:'
           write(6,'(a)')  '    --geom         (-g, --geometry)'
           write(6,'(a)')  '    --load         (-l, --loadcase)'
           write(6,'(a)')  '    --workingdir   (-w, --wd, --workingdirectory, -d, --directory)'
           write(6,'(a)')  '    --restart      (-r, --rs)'
           write(6,'(a)')  '    --regrid       (--rg)'
           write(6,'(a)')  '    --help         (-h)'
           write(6,'(/,a)')' -----------------------------------------------------------------------'
           write(6,'(a)')  ' Mandatory arguments:'
           write(6,'(/,a)')'   --geom PathToGeomFile/NameOfGeom.geom'
           write(6,'(a)')  '        Specifies the location of the geometry definition file,' 
           write(6,'(a)')  '            if no extension is given, .geom will be appended.' 
           write(6,'(a)')  '        "PathToGeomFile" will be the working directory if not specified'
           write(6,'(a)')  '            via --workingdir.'
           write(6,'(a)')  '        Make sure the file "material.config" exists in the working'
           write(6,'(a)')  '            directory.'   
           write(6,'(a)')  '        For further configuration place "numerics.config"'
           write(6,'(a)')'            and "numerics.config" in that directory.'
           write(6,'(/,a)')'   --load PathToLoadFile/NameOfLoadFile.load'
           write(6,'(a)')  '        Specifies the location of the load case definition file,' 
           write(6,'(a)')  '            if no extension is given, .load will be appended.' 
           write(6,'(/,a)')' -----------------------------------------------------------------------'
           write(6,'(a)')  ' Optional arguments:'
           write(6,'(/,a)')'   --workingdirectory PathToWorkingDirectory'
           write(6,'(a)')  '        Specifies the working directory and overwrites the default'
           write(6,'(a)')  '            "PathToGeomFile".'
           write(6,'(a)')  '        Make sure the file "material.config" exists in the working'
           write(6,'(a)')  '            directory.'   
           write(6,'(a)')  '        For further configuration place "numerics.config"'
           write(6,'(a)')'            and "numerics.config" in that directory.'
           write(6,'(/,a)')'   --restart XX'
           write(6,'(a)')  '        Reads in total increment No. XX-1 and continues to'
           write(6,'(a)')  '            calculate total increment No. XX.'
           write(6,'(a)')  '        Appends to existing results file '
           write(6,'(a)')  '            "NameOfGeom_NameOfLoadFile.spectralOut".'
           write(6,'(a)')  '        Works only if the restart information for total increment'
           write(6,'(a)')  '             No. XX-1 is available in the working directory.'
           write(6,'(/,a)')'   --regrid XX'
           write(6,'(a)')  '        Reads in total increment No. XX-1 and continues to'
           write(6,'(a)')  '            calculate total increment No. XX.'
           write(6,'(a)')  '        Attention: Overwrites existing results file '
           write(6,'(a)')  '            "NameOfGeom_NameOfLoadFile.spectralOut".'
           write(6,'(a)')  '        Works only if the restart information for total increment'
           write(6,'(a)')  '             No. XX-1 is available in the working directory.'
           write(6,'(/,a)')' -----------------------------------------------------------------------'
           write(6,'(a)')  ' Help:'
           write(6,'(/,a)')'   --help'
           write(6,'(a,/)')'        Prints this message and exits'
           call quit(0_pInt)                                                                              ! normal Termination
         endif mainProcess2
       case ('-l', '--load', '--loadcase')
         loadcaseArg = IIO_stringValue(commandLine,positions,i+1_pInt)
       case ('-g', '--geom', '--geometry')
         geometryArg = IIO_stringValue(commandLine,positions,i+1_pInt)
       case ('-w', '-d', '--wd', '--directory', '--workingdir', '--workingdirectory')
         workingDirArg = IIO_stringValue(commandLine,positions,i+1_pInt)
       case ('-r', '--rs', '--restart')
         spectralRestartInc = IIO_IntValue(commandLine,positions,i+1_pInt)
         appendToOutFile = .true.
       case ('--rg', '--regrid')
         spectralRestartInc = IIO_IntValue(commandLine,positions,i+1_pInt)
         appendToOutFile = .false.
     end select
   enddo
 endif
 
 if (len(trim(loadcaseArg)) == 0 .or. len(trim(geometryArg)) == 0) then
   write(6,'(a)') ' Please specify geometry AND load case (-h for help)'
   call quit(1_pInt)
 endif

 workingDirectory = storeWorkingDirectory(trim(workingDirArg),trim(geometryArg))
 geometryFile = getGeometryFile(geometryArg)
 loadCaseFile = getLoadCaseFile(loadCaseArg)

 call get_environment_variable('HOSTNAME',hostName)
 call get_environment_variable('USER',userName)
 call date_and_time(values = dateAndTime)
 mainProcess3: if (worldrank == 0) then
   write(6,'(a,2(i2.2,a),i4.4)') ' Date:                  ',dateAndTime(3),'/',&
                                                            dateAndTime(2),'/',&
                                                            dateAndTime(1) 
   write(6,'(a,2(i2.2,a),i2.2)') ' Time:                  ',dateAndTime(5),':',&
                                                            dateAndTime(6),':',&
                                                            dateAndTime(7)  
   write(6,'(a,a)')      ' Host name:             ', trim(hostName)
   write(6,'(a,a)')      ' User name:             ', trim(userName)
   write(6,'(a,a)')      ' Path separator:        ', getPathSep()
   write(6,'(a,a)')      ' Command line call:     ', trim(commandLine)
   if (len(trim(workingDirArg))>0) &
     write(6,'(a,a)')    ' Working dir argument:  ', trim(workingDirArg)
   write(6,'(a,a)')      ' Geometry argument:     ', trim(geometryArg)
   write(6,'(a,a)')      ' Loadcase argument:     ', trim(loadcaseArg)
   write(6,'(a,a)')      ' Working directory:     ', trim(getSolverWorkingDirectoryName())
   write(6,'(a,a)')      ' Geometry file:         ', trim(geometryFile)
   write(6,'(a,a)')      ' Loadcase file:         ', trim(loadCaseFile)
   write(6,'(a,a)')      ' Solver job name:       ', trim(getSolverJobName())
   if (SpectralRestartInc > 1_pInt) &
     write(6,'(a,i6.6)') ' Restart at increment:  ', spectralRestartInc
   write(6,'(a,l1,/)')   ' Append to result file: ', appendToOutFile
 endif mainProcess3

end subroutine DAMASK_interface_init


!--------------------------------------------------------------------------------------------------
!> @brief extract working directory from given argument or from location of geometry file,
!!        possibly converting relative arguments to absolut path
!> @todo  change working directory with call chdir(storeWorkingDirectory)?
!--------------------------------------------------------------------------------------------------
character(len=1024) function storeWorkingDirectory(workingDirectoryArg,geometryArg)
#ifdef __INTEL_COMPILER
 use IFPORT
#endif

 implicit none
 character(len=*),  intent(in) :: workingDirectoryArg                                               !< working directory argument
 character(len=*),  intent(in) :: geometryArg                                                       !< geometry argument
 character(len=1024)           :: cwd
 character                     :: pathSep
 logical                       :: dirExists
 external                      :: quit
 integer                       :: error

 pathSep = getPathSep()
 if (len(workingDirectoryArg)>0) then                                                               ! got working directory as input
   if (workingDirectoryArg(1:1) == pathSep) then                                                    ! absolute path given as command line argument
     storeWorkingDirectory = workingDirectoryArg
   else
     error = getcwd(cwd)                                                                            ! relative path given as command line argument
     storeWorkingDirectory = trim(cwd)//pathSep//workingDirectoryArg
   endif
   if (storeWorkingDirectory(len(trim(storeWorkingDirectory)):len(trim(storeWorkingDirectory))) &   ! if path seperator is not given, append it
      /= pathSep) storeWorkingDirectory = trim(storeWorkingDirectory)//pathSep
#ifdef __INTEL_COMPILER
   inquire(directory = trim(storeWorkingDirectory)//'.', exist=dirExists)
#else
   inquire(file = trim(storeWorkingDirectory), exist=dirExists)
#endif
   if(.not. dirExists) then                                                                         ! check if the directory exists
     write(6,'(a20,a,a16)') ' working directory "',trim(storeWorkingDirectory),'" does not exist'
     call quit(1_pInt)
   endif
 else                                                                                               ! using path to geometry file as working dir
   if (geometryArg(1:1) == pathSep) then                                                            ! absolute path given as command line argument
     storeWorkingDirectory = geometryArg(1:scan(geometryArg,pathSep,back=.true.))
   else
     error = getcwd(cwd)                                                                            ! relative path given as command line argument
     storeWorkingDirectory = trim(cwd)//pathSep//&
                              geometryArg(1:scan(geometryArg,pathSep,back=.true.))
   endif
 endif
 storeWorkingDirectory = rectifyPath(storeWorkingDirectory)
 
end function storeWorkingDirectory


!--------------------------------------------------------------------------------------------------
!> @brief simply returns the private string workingDir 
!--------------------------------------------------------------------------------------------------
character(len=1024) function getSolverWorkingDirectoryName()

 implicit none
 getSolverWorkingDirectoryName = workingDirectory

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
#ifdef __INTEL_COMPILER
 use IFPORT
#endif

 implicit none
 character(len=1024), intent(in) :: &
   geometryParameter
 character(len=1024) :: &
   cwd
 integer :: posExt, posSep
 character :: pathSep
 integer :: error

 getGeometryFile = geometryParameter
 pathSep = getPathSep()
 posExt = scan(getGeometryFile,'.',back=.true.)
 posSep = scan(getGeometryFile,pathSep,back=.true.)

 if (posExt <= posSep) getGeometryFile = trim(getGeometryFile)//('.geom')                           ! no extension present
 if (scan(getGeometryFile,pathSep) /= 1) then                                                       ! relative path given as command line argument
   error = getcwd(cwd)
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
#ifdef __INTEL_COMPILER
 use IFPORT
#endif

 implicit none
 character(len=1024), intent(in) :: &
   loadCaseParameter
 character(len=1024) :: &
   cwd
 integer :: posExt, posSep, error
 character :: pathSep

 getLoadCaseFile = loadcaseParameter
 pathSep = getPathSep()
 posExt = scan(getLoadCaseFile,'.',back=.true.)
 posSep = scan(getLoadCaseFile,pathSep,back=.true.)

 if (posExt <= posSep) getLoadCaseFile = trim(getLoadCaseFile)//('.load')                           ! no extension present
 if (scan(getLoadCaseFile,pathSep) /= 1) then                                                       ! relative path given as command line argument
   error = getcwd(cwd)
   getLoadCaseFile = rectifyPath(trim(cwd)//pathSep//getLoadCaseFile)
 else
   getLoadCaseFile = rectifyPath(getLoadCaseFile)
 endif

 getLoadCaseFile = makeRelativePath(getSolverWorkingDirectoryName(), getLoadCaseFile)

end function getLoadCaseFile


!--------------------------------------------------------------------------------------------------
!> @brief remove ../ and /./ from path
!--------------------------------------------------------------------------------------------------
function rectifyPath(path)

 implicit none
 character(len=*) :: path
 character(len=len_trim(path)) :: rectifyPath
 character :: pathSep
 integer :: i,j,k,l                                                                                 ! no pInt

 pathSep = getPathSep()

!--------------------------------------------------------------------------------------------------
! remove /./ from path
 l = len_trim(path)
 rectifyPath = path
 do i = l,3,-1
   if (rectifyPath(i-2:i) == pathSep//'.'//pathSep) &
     rectifyPath(i-1:l) = rectifyPath(i+1:l)//'  '
 enddo

!--------------------------------------------------------------------------------------------------
! remove ../ and corresponding directory from rectifyPath
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
! to be the path separator
!--------------------------------------------------------------------------------------------------
character function getPathSep()

 implicit none
 character(len=2048) :: &
   path
 integer(pInt) :: &
   backslash = 0_pInt, &
   slash = 0_pInt
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


!--------------------------------------------------------------------------------------------------
!> @brief taken from IO, check IO_stringValue for documentation 
!--------------------------------------------------------------------------------------------------
pure function IIO_stringValue(line,positions,myPos)
 
 implicit none
 integer(pInt),                                intent(in) :: positions(*), &
                                                             myPos
 character(len=1+positions(myPos*2+1)-positions(myPos*2)) :: IIO_stringValue
 character(len=*),                             intent(in) :: line

 if (positions(1) < myPos) then
   IIO_stringValue = ''
 else
   IIO_stringValue = line(positions(myPos*2):positions(myPos*2+1))
 endif

end function IIO_stringValue


!--------------------------------------------------------------------------------------------------
!> @brief taken from IO, check IO_stringValue for documentation 
!--------------------------------------------------------------------------------------------------
integer(pInt) pure function IIO_intValue(line,positions,myPos)

 implicit none
 character(len=*), intent(in) :: line
 integer(pInt),    intent(in) :: positions(*), &
                                 myPos

 if (positions(1) < myPos) then
   IIO_intValue = 0_pInt
 else
   read(UNIT=line(positions(myPos*2):positions(myPos*2+1)),ERR=100,FMT=*) IIO_intValue
 endif
 return
100 IIO_intValue = huge(1_pInt)

end function IIO_intValue


!--------------------------------------------------------------------------------------------------
!> @brief taken from IO, check IO_lc for documentation 
!--------------------------------------------------------------------------------------------------
pure function IIO_lc(line)

 implicit none
 character(26), parameter :: lower = 'abcdefghijklmnopqrstuvwxyz'
 character(26), parameter :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
 character(len=*), intent(in) :: line
 character(len=len(line))     :: IIO_lc
 
 integer                      :: i,n                                                                ! no pInt (len returns default integer)

 IIO_lc = line
 do i=1,len(line)
   n = index(upper,IIO_lc(i:i))
   if (n/=0) IIO_lc(i:i) = lower(n:n)
 enddo

end function IIO_lc


!--------------------------------------------------------------------------------------------------
!> @brief taken from IO, check IO_stringPos for documentation 
!--------------------------------------------------------------------------------------------------
pure function IIO_stringPos(line,N)

 implicit none
 integer(pInt),    intent(in) :: N
 integer(pInt)                :: IIO_stringPos(1_pInt+N*2_pInt)
 
 character(len=*), intent(in) :: line
 
 character(len=*), parameter  :: sep=achar(44)//achar(32)//achar(9)//achar(10)//achar(13)           ! comma and whitespaces

 integer                      :: left, right                                                        ! no pInt (verify and scan return default integer)


 IIO_stringPos = -1_pInt
 IIO_stringPos(1) = 0_pInt
 right = 0
 
 do while (verify(line(right+1:),sep)>0)
   left  = right + verify(line(right+1:),sep)
   right = left + scan(line(left:),sep) - 2
   if ( line(left:left) == '#' ) then
     exit
   endif
   if ( IIO_stringPos(1)<N ) then
     IIO_stringPos(1_pInt+IIO_stringPos(1)*2_pInt+1_pInt) = int(left, pInt)
     IIO_stringPos(1_pInt+IIO_stringPos(1)*2_pInt+2_pInt) = int(right, pInt)
   endif
   IIO_stringPos(1) = IIO_stringPos(1)+1_pInt
 enddo

end function IIO_stringPos


end module
