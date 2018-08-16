!--------------------------------------------------------------------------------------------------
!> @author   Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief    Interfacing between the FEM solvers and the material subroutines provided
!!           by DAMASK
!> @details  Interfacing between the FEM solvers and the material subroutines provided
!>           by DAMASK. Interpretating the command line arguments to the init routine to 
!>           get load case, geometry file, working directory, etc.
!--------------------------------------------------------------------------------------------------
module DAMASK_interface
 use prec, only: &
   pInt

 implicit none
 private
 logical,             public, protected :: appendToOutFile = .false.                                !< Append to existing output file
 integer(pInt),       public, protected :: FEMRestartInc = 0_pInt                                   !< Increment at which calculation starts
 character(len=1024), public, protected :: &
   geometryFile = '', &                                                                             !< parameter given for geometry file
   loadCaseFile = ''                                                                                !< parameter given for load case file
 character(len=1024), private           :: workingDirectory

 public :: &
   getSolverJobName, &
   DAMASK_interface_init
 private :: &
   setWorkingDirectory, &
   getGeometryFile, &
   getLoadCaseFile, &
   rectifyPath, &
   makeRelativePath, &
   IIO_stringValue, &
   IIO_intValue, &
   IIO_stringPos
contains

!--------------------------------------------------------------------------------------------------
!> @brief initializes the solver by interpreting the command line arguments. Also writes
!! information on computation to screen
!--------------------------------------------------------------------------------------------------
subroutine DAMASK_interface_init()
 use, intrinsic :: &
   iso_fortran_env
#include <petsc/finclude/petscsys.h>
#if PETSC_VERSION_MAJOR!=3 || PETSC_VERSION_MINOR!=9
=================================================================================================== 
=========================   THIS VERSION OF DAMASK REQUIRES PETSc 3.9.x   ========================= 
===================================================================================================
#endif
 use PETScSys
 use system_routines, only: &
   getHostName

 implicit none
 character(len=1024) :: &
   commandLine, &                                                                                   !< command line call as string
   loadcaseArg   = '', &                                                                            !< -l argument given to DAMASK_FEM.exe
   geometryArg   = '', &                                                                            !< -g argument given to DAMASK_FEM.exe
   workingDirArg = '', &                                                                            !< -w argument given to DAMASK_FEM.exe
   hostName, &                                                                                      !< name of machine on which DAMASK_FEM.exe is execute (might require export HOSTNAME)
   userName, &                                                                                      !< name of user calling DAMASK_FEM.exe
   tag
 integer :: &
   i, &
#ifdef _OPENMP
   threadLevel, &
#endif
   worldrank = 0, &
   worldsize = 0
 integer, allocatable, dimension(:) :: &
   chunkPos
 integer, dimension(8) :: &
   dateAndTime                                                                                      ! type default integer
 PetscErrorCode :: ierr
 logical        :: error
 external :: &
   quit,&
   PETScErrorF, &                                                                                   ! is called in the CHKERRQ macro
   PETScInitialize

 open(6, encoding='UTF-8')                                                                          ! for special characters in output

!--------------------------------------------------------------------------------------------------
! PETSc Init
#ifdef _OPENMP
 ! If openMP is enabled, check if the MPI libary supports it and initialize accordingly.
 ! Otherwise, the first call to PETSc will do the initialization.
 call MPI_Init_Thread(MPI_THREAD_FUNNELED,threadLevel,ierr);CHKERRQ(ierr)
 if (threadLevel<MPI_THREAD_FUNNELED) then
   write(6,'(a)') ' MPI library does not support OpenMP'
   call quit(1_pInt)
 endif
#endif
 call PETScInitialize(PETSC_NULL_CHARACTER,ierr)                                                    ! according to PETSc manual, that should be the first line in the code
 CHKERRQ(ierr)                                                                                      ! this is a macro definition, it is case sensitive
 call MPI_Comm_rank(PETSC_COMM_WORLD,worldrank,ierr);CHKERRQ(ierr)
 call MPI_Comm_size(PETSC_COMM_WORLD,worldsize,ierr);CHKERRQ(ierr)
 mainProcess: if (worldrank == 0) then
   if (output_unit /= 6) then
     write(output_unit,'(a)') ' STDOUT != 6'
     call quit(1_pInt)
   endif
   if (error_unit /= 0) then
     write(output_unit,'(a)') ' STDERR != 0'
     call quit(1_pInt)
   endif
 else mainProcess
   close(6)                                                                                         ! disable output for non-master processes (open 6 to rank specific file for debug)
   open(6,file='/dev/null',status='replace')                                                        ! close(6) alone will leave some temp files in cwd
 endif mainProcess

 call date_and_time(values = dateAndTime)
 write(6,'(/,a)') ' <<<+-  DAMASK_FEM  -+>>>'
 write(6,'(a,/)') ' Roters et al., Computational Materials Science, 2018'
 write(6,'(/,a)')              ' Version: '//DAMASKVERSION
 write(6,'(a,2(i2.2,a),i4.4)') ' Date:    ',dateAndTime(3),'/',&
                                            dateAndTime(2),'/',&
                                            dateAndTime(1) 
 write(6,'(a,2(i2.2,a),i2.2)') ' Time:    ',dateAndTime(5),':',&
                                            dateAndTime(6),':',&
                                            dateAndTime(7)  
 write(6,'(/,a,i4.1)') ' MPI processes: ',worldsize
 write(6,'(/,a)') ' <<<+-  DAMASK_interface init  -+>>>'
#include "compilation_info.f90"
 
 call get_command(commandLine)
 chunkPos = IIO_stringPos(commandLine)
 do i = 2_pInt, chunkPos(1)
   select case(IIO_stringValue(commandLine,chunkPos,i))                                             ! extract key
     case ('-h','--help')
       write(6,'(a)')  ' #######################################################################'
       write(6,'(a)')  ' DAMASK_FEM:'
       write(6,'(a)')  ' FEM solvers for the Düsseldorf Advanced Material Simulation Kit'
       write(6,'(a,/)')' #######################################################################'
       write(6,'(a,/)')' Valid command line switches:'
       write(6,'(a)')  '    --geom         (-g, --geometry)'
       write(6,'(a)')  '    --load         (-l, --loadcase)'
       write(6,'(a)')  '    --workingdir   (-w, --wd, --workingdirectory, -d, --directory)'
       write(6,'(a)')  '    --restart      (-r, --rs)'
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
       write(6,'(a)')'            and "debug.config" in that directory.'
       write(6,'(/,a)')'   --restart XX'
       write(6,'(a)')  '        Reads in increment XX and continues with calculating'
       write(6,'(a)')  '            increment XX+1 based on this.'
       write(6,'(a)')  '        Appends to existing results file'
       write(6,'(a)')  '            "NameOfGeom_NameOfLoadFile.YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY".'
       write(6,'(a)')  '        Works only if the restart information for increment XX'
       write(6,'(a)')  '            is available in the working directory.'
       write(6,'(/,a)')' -----------------------------------------------------------------------'
       write(6,'(a)')  ' Help:'
       write(6,'(/,a)')'   --help'
       write(6,'(a,/)')'        Prints this message and exits'
       call quit(0_pInt)                                                                            ! normal Termination
     case ('-l', '--load', '--loadcase')
       if ( i < chunkPos(1)) loadcaseArg = trim(IIO_stringValue(commandLine,chunkPos,i+1_pInt))
     case ('-g', '--geom', '--geometry')
       if (i < chunkPos(1)) geometryArg = trim(IIO_stringValue(commandLine,chunkPos,i+1_pInt))
     case ('-w', '-d', '--wd', '--directory', '--workingdir', '--workingdirectory')
       if (i < chunkPos(1)) workingDirArg = trim(IIO_stringValue(commandLine,chunkPos,i+1_pInt))
     case ('-r', '--rs', '--restart')
       if (i < chunkPos(1)) then
         FEMRestartInc = IIO_IntValue(commandLine,chunkPos,i+1_pInt)
         appendToOutFile = .true.
       endif
   end select
 enddo

 if (len_trim(loadcaseArg) == 0 .or. len_trim(geometryArg) == 0) then
   write(6,'(a)') ' Please specify geometry AND load case (-h for help)'
   call quit(1_pInt)
 endif

 workingDirectory = trim(setWorkingDirectory(trim(workingDirArg)))
 geometryFile = getGeometryFile(geometryArg)
 loadCaseFile = getLoadCaseFile(loadCaseArg)

 call get_environment_variable('USER',userName)
 error = getHostName(hostName)
 write(6,'(a,a)')      ' Host name:              ', trim(hostName)
 write(6,'(a,a)')      ' User name:              ', trim(userName)
 write(6,'(a,a)')      ' Command line call:      ', trim(commandLine)
 if (len(trim(workingDirArg)) > 0) &
   write(6,'(a,a)')    ' Working dir argument:   ', trim(workingDirArg)
 write(6,'(a,a)')      ' Geometry argument:      ', trim(geometryArg)
 write(6,'(a,a)')      ' Loadcase argument:      ', trim(loadcaseArg)
 write(6,'(a,a)')      ' Working directory:      ', trim(workingDirectory)
 write(6,'(a,a)')      ' Geometry file:          ', trim(geometryFile)
 write(6,'(a,a)')      ' Loadcase file:          ', trim(loadCaseFile)
 write(6,'(a,a)')      ' Solver job name:        ', trim(getSolverJobName())
 if (SpectralRestartInc > 0_pInt) &
   write(6,'(a,i6.6)') ' Restart from increment: ', FEMRestartInc
 write(6,'(a,l1,/)')   ' Append to result file:  ', appendToOutFile

end subroutine DAMASK_interface_init


!--------------------------------------------------------------------------------------------------
!> @brief extract working directory from given argument or from location of geometry file,
!!        possibly converting relative arguments to absolut path
!--------------------------------------------------------------------------------------------------
character(len=1024) function setWorkingDirectory(workingDirectoryArg)
 use system_routines, only: &
   getCWD, &
   setCWD

 implicit none
 character(len=*),  intent(in) :: workingDirectoryArg                                               !< working directory argument
 logical                       :: error
 external                      :: quit

 wdGiven: if (len(workingDirectoryArg)>0) then
   absolutePath: if (workingDirectoryArg(1:1) == '/') then
     setWorkingDirectory = workingDirectoryArg
   else absolutePath
     error = getCWD(setWorkingDirectory)
     if (error) call quit(1_pInt)
     setWorkingDirectory = trim(setWorkingDirectory)//'/'//workingDirectoryArg
   endif absolutePath
 else wdGiven
   error = getCWD(setWorkingDirectory)                                                              ! relative path given as command line argument
   if (error) call quit(1_pInt)
 endif wdGiven

 setWorkingDirectory = trim(rectifyPath(setWorkingDirectory))

 error = setCWD(trim(setWorkingDirectory))
 if(error) then
   write(6,'(a20,a,a16)') ' working directory "',trim(setWorkingDirectory),'" does not exist'
   call quit(1_pInt)
 endif

end function setWorkingDirectory


!--------------------------------------------------------------------------------------------------
!> @brief solver job name (no extension) as combination of geometry and load case name
!--------------------------------------------------------------------------------------------------
character(len=1024) function getSolverJobName()

 implicit none
 integer :: posExt,posSep
 character(len=1024) :: tempString 


 tempString = geometryFile
 posExt = scan(tempString,'.',back=.true.)
 posSep = scan(tempString,'/',back=.true.)

 getSolverJobName = tempString(posSep+1:posExt-1)

 tempString = loadCaseFile
 posExt = scan(tempString,'.',back=.true.)
 posSep = scan(tempString,'/',back=.true.)

 getSolverJobName = trim(getSolverJobName)//'_'//tempString(posSep+1:posExt-1)

end function getSolverJobName


!--------------------------------------------------------------------------------------------------
!> @brief basename of geometry file with extension from command line arguments
!--------------------------------------------------------------------------------------------------
character(len=1024) function getGeometryFile(geometryParameter)

 implicit none
 character(len=1024), intent(in) :: &
   geometryParameter
 integer :: posExt, posSep
 external  :: quit

 getGeometryFile = trim(geometryParameter)
 posExt = scan(getGeometryFile,'.',back=.true.)
 posSep = scan(getGeometryFile,'/',back=.true.)

 if (posExt <= posSep) getGeometryFile = trim(getGeometryFile)//('.geom')
 if (scan(getGeometryFile,'/') /= 1) &
   getGeometryFile = trim(workingDirectory)//'/'//trim(getGeometryFile)

 getGeometryFile = makeRelativePath(workingDirectory, getGeometryFile)


end function getGeometryFile


!--------------------------------------------------------------------------------------------------
!> @brief relative path of loadcase from command line arguments
!--------------------------------------------------------------------------------------------------
character(len=1024) function getLoadCaseFile(loadCaseParameter)

 implicit none
 character(len=1024), intent(in) :: &
   loadCaseParameter
 integer :: posExt, posSep
 external  :: quit

 getLoadCaseFile = trim(loadCaseParameter)
 posExt = scan(getLoadCaseFile,'.',back=.true.)
 posSep = scan(getLoadCaseFile,'/',back=.true.)

 if (posExt <= posSep) getLoadCaseFile = trim(getLoadCaseFile)//('.load')
 if (scan(getLoadCaseFile,'/') /= 1) &
   getLoadCaseFile = trim(workingDirectory)//'/'//trim(getLoadCaseFile)

 getLoadCaseFile = makeRelativePath(workingDirectory, getLoadCaseFile)

end function getLoadCaseFile


!--------------------------------------------------------------------------------------------------
!> @brief remove ../, /./, and // from path.
!> @details works only if absolute path is given
!--------------------------------------------------------------------------------------------------
function rectifyPath(path)

 implicit none
 character(len=*) :: path
 character(len=len_trim(path)) :: rectifyPath
 integer :: i,j,k,l                                                                                 ! no pInt

!--------------------------------------------------------------------------------------------------
! remove /./ from path
 l = len_trim(path)
 rectifyPath = path
 do i = l,3,-1
   if (rectifyPath(i-2:i) == '/./') rectifyPath(i-1:l) = rectifyPath(i+1:l)//'  '
 enddo

!--------------------------------------------------------------------------------------------------
! remove // from path
 l = len_trim(path)
 rectifyPath = path
 do i = l,2,-1
   if (rectifyPath(i-1:i) == '//') rectifyPath(i-1:l) = rectifyPath(i:l)//' '
 enddo

!--------------------------------------------------------------------------------------------------
! remove ../ and corresponding directory from rectifyPath
 l = len_trim(rectifyPath)
 i = index(rectifyPath(i:l),'../')
 j = 0
 do while (i > j)
    j = scan(rectifyPath(1:i-2),'/',back=.true.)
    rectifyPath(j+1:l) = rectifyPath(i+3:l)//repeat(' ',2+i-j)
    if (rectifyPath(j+1:j+1) == '/') then                                                           !search for '//' that appear in case of XXX/../../XXX
      k = len_trim(rectifyPath)
      rectifyPath(j+1:k-1) = rectifyPath(j+2:k)
      rectifyPath(k:k) = ' '
    endif
    i = j+index(rectifyPath(j+1:l),'../')
 enddo
 if(len_trim(rectifyPath) == 0) rectifyPath = '/'

end function rectifyPath

 
!--------------------------------------------------------------------------------------------------
!> @brief relative path from absolute a to absolute b
!--------------------------------------------------------------------------------------------------
character(len=1024) function makeRelativePath(a,b)

 implicit none
 character (len=*), intent(in) :: a,b
 character (len=1024)          :: a_cleaned,b_cleaned
 integer :: i,posLastCommonSlash,remainingSlashes !no pInt

 posLastCommonSlash = 0
 remainingSlashes = 0
 a_cleaned = rectifyPath(trim(a)//'/')
 b_cleaned = rectifyPath(b)

 do i = 1, min(1024,len_trim(a_cleaned),len_trim(rectifyPath(b_cleaned)))
   if (a_cleaned(i:i) /= b_cleaned(i:i)) exit
   if (a_cleaned(i:i) == '/') posLastCommonSlash = i
 enddo
 do i = posLastCommonSlash+1,len_trim(a_cleaned)
   if (a_cleaned(i:i) == '/') remainingSlashes = remainingSlashes + 1
 enddo

 makeRelativePath = repeat('..'//'/',remainingSlashes)//b_cleaned(posLastCommonSlash+1:len_trim(b_cleaned))

end function makeRelativePath


!--------------------------------------------------------------------------------------------------
!> @brief taken from IO, check IO_stringValue for documentation 
!--------------------------------------------------------------------------------------------------
pure function IIO_stringValue(string,chunkPos,myChunk)
 
 implicit none
 integer(pInt),    dimension(:),               intent(in)   :: chunkPos                             !< positions of start and end of each tag/chunk in given string
 integer(pInt),                                intent(in)   :: myChunk                              !< position number of desired chunk
 character(len=chunkPos(myChunk*2+1)-chunkPos(myChunk*2)+1) :: IIO_stringValue
 character(len=*),                             intent(in)   :: string                               !< raw input with known start and end of each chunk

 IIO_stringValue = string(chunkPos(myChunk*2):chunkPos(myChunk*2+1))

end function IIO_stringValue


!--------------------------------------------------------------------------------------------------
!> @brief taken from IO, check IO_intValue for documentation 
!--------------------------------------------------------------------------------------------------
integer(pInt) pure function IIO_intValue(string,chunkPos,myChunk)                                           
                                                                                                    
 implicit none                                                                                      
 character(len=*),               intent(in) :: string                                               !< raw input with known start and end of each chunk
 integer(pInt),                  intent(in) :: myChunk                                              !< position number of desired sub string
 integer(pInt),   dimension(:),  intent(in) :: chunkPos                                             !< positions of start and end of each tag/chunk in given string


 valuePresent: if (myChunk > chunkPos(1) .or. myChunk < 1_pInt) then
   IIO_intValue = 0_pInt
 else valuePresent
   read(UNIT=string(chunkPos(myChunk*2):chunkPos(myChunk*2+1)),ERR=100,FMT=*) IIO_intValue
 endif valuePresent
 return
100 IIO_intValue = huge(1_pInt)

end function IIO_intValue


!--------------------------------------------------------------------------------------------------
!> @brief taken from IO, check IO_stringPos for documentation 
!--------------------------------------------------------------------------------------------------
pure function IIO_stringPos(string)

 implicit none
 integer(pInt), dimension(:), allocatable            :: IIO_stringPos
 character(len=*),                        intent(in) :: string                                      !< string in which chunks are searched for
 
 character(len=*), parameter  :: SEP=achar(44)//achar(32)//achar(9)//achar(10)//achar(13)           ! comma and whitespaces
 integer                      :: left, right                                                        ! no pInt (verify and scan return default integer)

 allocate(IIO_stringPos(1), source=0_pInt)
 right = 0
 
 do while (verify(string(right+1:),SEP)>0)
   left  = right + verify(string(right+1:),SEP)
   right = left + scan(string(left:),SEP) - 2
   if ( string(left:left) == '#' ) exit
   IIO_stringPos = [IIO_stringPos,int(left, pInt), int(right, pInt)]
   IIO_stringPos(1) = IIO_stringPos(1)+1_pInt
 enddo

end function IIO_stringPos

end module
