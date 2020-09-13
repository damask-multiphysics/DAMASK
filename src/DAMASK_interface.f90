!--------------------------------------------------------------------------------------------------
!> @author   Jaeyong Jung, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief    Interfacing between the PETSc-based solvers and the material subroutines provided
!!           by DAMASK
!> @details  Interfacing between the PETSc-based solvers and the material subroutines provided
!>           by DAMASK. Interpreting the command line arguments to get load case, geometry file,
!>           and working directory.
!--------------------------------------------------------------------------------------------------
#define PETSC_MAJOR 3
#define PETSC_MINOR_MIN 10
#define PETSC_MINOR_MAX 13

module DAMASK_interface
  use, intrinsic :: iso_fortran_env

  use PETScSys

  use prec
  use parallelization
  use system_routines

  implicit none
  private
  logical,          volatile,    public, protected :: &
    interface_SIGTERM, &                                                                            !< termination signal
    interface_SIGUSR1, &                                                                            !< 1. user-defined signal
    interface_SIGUSR2                                                                               !< 2. user-defined signal
  integer,                       public, protected :: &
    interface_restartInc = 0                                                                        !< Increment at which calculation starts
  character(len=:), allocatable, public, protected :: &
    interface_geomFile, &                                                                           !< parameter given for geometry file
    interface_loadFile                                                                              !< parameter given for load case file

  public :: &
    getSolverJobName, &
    DAMASK_interface_init, &
    interface_setSIGTERM, &
    interface_setSIGUSR1, &
    interface_setSIGUSR2

contains

!--------------------------------------------------------------------------------------------------
!> @brief initializes the solver by interpreting the command line arguments. Also writes
!! information on computation to screen
!--------------------------------------------------------------------------------------------------
subroutine DAMASK_interface_init
#include <petsc/finclude/petscsys.h>

#if PETSC_VERSION_MAJOR!=3 || PETSC_VERSION_MINOR<PETSC_MINOR_MIN || PETSC_VERSION_MINOR>PETSC_MINOR_MAX
===================================================================================================
--  WRONG PETSc VERSION --- WRONG PETSc VERSION --- WRONG PETSc VERSION ---  WRONG PETSc VERSION --
===================================================================================================
============   THIS VERSION OF DAMASK REQUIRES A DIFFERENT PETSc VERSION   ========================
===============   THIS VERSION OF DAMASK REQUIRES A DIFFERENT PETSc VERSION   =====================
==================   THIS VERSION OF DAMASK REQUIRES A DIFFERENT PETSc VERSION   ==================
===================================================================================================
--  WRONG PETSc VERSION --- WRONG PETSc VERSION --- WRONG PETSc VERSION ---  WRONG PETSc VERSION --
===================================================================================================
#endif

  character(len=pPathLen*3+pStringLen) :: &
    commandLine                                                                                     !< command line call as string
  character(len=pPathLen) :: &
    arg, &                                                                                          !< individual argument
    loadCaseArg   = '', &                                                                           !< -l argument given to the executable
    geometryArg   = '', &                                                                           !< -g argument given to the executable
    workingDirArg = ''                                                                              !< -w argument given to the executable
  character(len=pStringLen) :: &
    userName                                                                                        !< name of user calling the executable
  integer :: &
    stat, &
    i
  integer, dimension(8) :: &
    dateAndTime
  integer        :: err
  PetscErrorCode :: petsc_err
  external :: &
    quit

  write(6,'(/,a)') ' <<<+-  DAMASK_interface init  -+>>>'

  open(6, encoding='UTF-8')                                                                         ! for special characters in output

 ! http://patorjk.com/software/taag/#p=display&f=Lean&t=DAMASK%203
#ifdef DEBUG
  print*, achar(27)//'[31m'
  write(6,'(a,/)') ' debug version - debug version - debug version - debug version - debug version'
#else
  print*, achar(27)//'[94m'
#endif
  print*, '     _/_/_/      _/_/    _/      _/    _/_/      _/_/_/  _/    _/    _/_/_/'
  print*, '    _/    _/  _/    _/  _/_/  _/_/  _/    _/  _/        _/  _/            _/'
  print*, '   _/    _/  _/_/_/_/  _/  _/  _/  _/_/_/_/    _/_/    _/_/          _/_/'
  print*, '  _/    _/  _/    _/  _/      _/  _/    _/        _/  _/  _/            _/'
  print*, ' _/_/_/    _/    _/  _/      _/  _/    _/  _/_/_/    _/    _/    _/_/_/'
#ifdef DEBUG
  write(6,'(/,a)') ' debug version - debug version - debug version - debug version - debug version'
#endif
  print*, achar(27)//'[0m'

  write(6,'(a)') ' Roters et al., Computational Materials Science 158:420–478, 2019'
  write(6,'(a)')   ' https://doi.org/10.1016/j.commatsci.2018.04.030'

  write(6,'(/,a)') ' Version: '//DAMASKVERSION

  ! https://github.com/jeffhammond/HPCInfo/blob/master/docs/Preprocessor-Macros.md
#if defined(__PGI)
  write(6,'(/,a,i4.4,a,i8.8)')   ' Compiled with PGI fortran version :', __PGIC__,&
                                                                    '.', __PGIC_MINOR__
#else
  write(6,'(/,a)') ' Compiled with: '//compiler_version()
  write(6,'(a)')   ' Compiler options: '//compiler_options()
#endif

  write(6,'(/,a)') ' Compiled on: '//__DATE__//' at '//__TIME__

  call date_and_time(values = dateAndTime)
  write(6,'(/,a,2(i2.2,a),i4.4)') ' Date: ',dateAndTime(3),'/',dateAndTime(2),'/', dateAndTime(1)
  write(6,'(a,2(i2.2,a),i2.2)')   ' Time: ',dateAndTime(5),':', dateAndTime(6),':', dateAndTime(7)

  do i = 1, command_argument_count()
    call get_command_argument(i,arg,status=err)
    if (err /= 0) call quit(1)
    select case(trim(arg))                                                                          ! extract key
      case ('-h','--help')
        write(6,'(a)')  ' #######################################################################'
        write(6,'(a)')  ' DAMASK Command Line Interface:'
        write(6,'(a)')  ' For PETSc-based solvers for the Düsseldorf Advanced Material Simulation Kit'
        write(6,'(a,/)')' #######################################################################'
        write(6,'(a,/)')' Valid command line switches:'
        write(6,'(a)')  '    --geom         (-g, --geometry)'
        write(6,'(a)')  '    --load         (-l, --loadcase)'
        write(6,'(a)')  '    --workingdir   (-w, --wd, --workingdirectory)'
        write(6,'(a)')  '    --restart      (-r, --rs)'
        write(6,'(a)')  '    --help         (-h)'
        write(6,'(/,a)')' -----------------------------------------------------------------------'
        write(6,'(a)')  ' Mandatory arguments:'
        write(6,'(/,a)')'   --geom PathToGeomFile/NameOfGeom'
        write(6,'(a)')  '        Specifies the location of the geometry definition file.'
        write(6,'(/,a)')'   --load PathToLoadFile/NameOfLoadFile'
        write(6,'(a)')  '        Specifies the location of the load case definition file.'
        write(6,'(/,a)')' -----------------------------------------------------------------------'
        write(6,'(a)')  ' Optional arguments:'
        write(6,'(/,a)')'   --workingdirectory PathToWorkingDirectory'
        write(6,'(a)')  '        Specifies the working directory and overwrites the default ./'
        write(6,'(a)')  '        Make sure the file "material.config" exists in the working'
        write(6,'(a)')  '            directory.'
        write(6,'(a)')  '        For further configuration place "numerics.config"'
        write(6,'(a)')'            and "debug.config" in that directory.'
        write(6,'(/,a)')'   --restart N'
        write(6,'(a)')  '        Reads in increment N and continues with calculating'
        write(6,'(a)')  '            increment N+1 based on this.'
        write(6,'(a)')  '        Appends to existing results file'
        write(6,'(a)')  '            "NameOfGeom_NameOfLoadFile.hdf5".'
        write(6,'(a)')  '        Works only if the restart information for increment N'
        write(6,'(a)')  '            is available in the working directory.'
        write(6,'(/,a)')' -----------------------------------------------------------------------'
        write(6,'(a)')  ' Help:'
        write(6,'(/,a)')'   --help'
        write(6,'(a,/)')'        Prints this message and exits'
        call quit(0)                                                                                ! normal Termination
      case ('-l', '--load', '--loadcase')
        call get_command_argument(i+1,loadCaseArg,status=err)
      case ('-g', '--geom', '--geometry')
        call get_command_argument(i+1,geometryArg,status=err)
      case ('-w', '--wd', '--workingdir', '--workingdirectory')
        call get_command_argument(i+1,workingDirArg,status=err)
      case ('-r', '--rs', '--restart')
        call get_command_argument(i+1,arg,status=err)
        read(arg,*,iostat=stat) interface_restartInc
        if (interface_restartInc < 0 .or. stat /=0) then
          write(6,'(/,a)') ' ERROR: Could not parse restart increment: '//trim(arg)
          call quit(1)
        endif
    end select
    if (err /= 0) call quit(1)
  enddo

  if (len_trim(loadcaseArg) == 0 .or. len_trim(geometryArg) == 0) then
    write(6,'(/,a)') ' ERROR: Please specify geometry AND load case (-h for help)'
    call quit(1)
  endif

  if (len_trim(workingDirArg) > 0) call setWorkingDirectory(trim(workingDirArg))
  interface_geomFile = getGeometryFile(geometryArg)
  interface_loadFile = getLoadCaseFile(loadCaseArg)

  call get_command(commandLine)
  call get_environment_variable('USER',userName)
  ! ToDo: https://stackoverflow.com/questions/8953424/how-to-get-the-username-in-c-c-in-linux
  write(6,'(/,a,i4.1)') ' MPI processes: ',worldsize
  write(6,'(a,a)')      ' Host name: ', trim(getHostName())
  write(6,'(a,a)')      ' User name: ', trim(userName)

  write(6,'(/a,a)')     ' Command line call:      ', trim(commandLine)
  if (len_trim(workingDirArg) > 0) &
    write(6,'(a,a)')    ' Working dir argument:   ', trim(workingDirArg)
  write(6,'(a,a)')      ' Geometry argument:      ', trim(geometryArg)
  write(6,'(a,a)')      ' Load case argument:     ', trim(loadcaseArg)
  write(6,'(a,a)')      ' Working directory:      ', getCWD()
  write(6,'(a,a)')      ' Geometry file:          ', interface_geomFile
  write(6,'(a,a)')      ' Loadcase file:          ', interface_loadFile
  write(6,'(a,a)')      ' Solver job name:        ', getSolverJobName()
  if (interface_restartInc > 0) &
    write(6,'(a,i6.6)') ' Restart from increment: ', interface_restartInc

  !call signalterm_c(c_funloc(catchSIGTERM))
  call signalusr1_c(c_funloc(catchSIGUSR1))
  call signalusr2_c(c_funloc(catchSIGUSR2))
  call interface_setSIGTERM(.false.)
  call interface_setSIGUSR1(.false.)
  call interface_setSIGUSR2(.false.)

end subroutine DAMASK_interface_init


!--------------------------------------------------------------------------------------------------
!> @brief extract working directory from given argument or from location of geometry file,
!!        possibly converting relative arguments to absolut path
!--------------------------------------------------------------------------------------------------
subroutine setWorkingDirectory(workingDirectoryArg)

  character(len=*),  intent(in) :: workingDirectoryArg                                              !< working directory argument
  character(len=pPathLen)       :: workingDirectory
  logical                       :: error
  external                      :: quit

  absolutePath: if (workingDirectoryArg(1:1) == '/') then
    workingDirectory = workingDirectoryArg
  else absolutePath
    workingDirectory = getCWD()
    workingDirectory = trim(workingDirectory)//'/'//workingDirectoryArg
  endif absolutePath

  workingDirectory = trim(rectifyPath(workingDirectory))
  error = setCWD(trim(workingDirectory))
  if(error) then
    write(6,'(/,a)') ' ERROR: Invalid Working directory: '//trim(workingDirectory)
    call quit(1)
  endif

end subroutine setWorkingDirectory


!--------------------------------------------------------------------------------------------------
!> @brief solver job name (no extension) as combination of geometry and load case name
!--------------------------------------------------------------------------------------------------
function getSolverJobName()

  character(len=:), allocatable :: getSolverJobName
  integer :: posExt,posSep

  posExt = scan(interface_geomFile,'.',back=.true.)
  posSep = scan(interface_geomFile,'/',back=.true.)

  getSolverJobName = interface_geomFile(posSep+1:posExt-1)

  posExt = scan(interface_loadFile,'.',back=.true.)
  posSep = scan(interface_loadFile,'/',back=.true.)

  getSolverJobName = getSolverJobName//'_'//interface_loadFile(posSep+1:posExt-1)

end function getSolverJobName


!--------------------------------------------------------------------------------------------------
!> @brief basename of geometry file with extension from command line arguments
!--------------------------------------------------------------------------------------------------
function getGeometryFile(geometryParameter)

  character(len=:), allocatable :: getGeometryFile
  character(len=*),  intent(in) :: geometryParameter
  logical                       :: file_exists
  external                      :: quit

  getGeometryFile = trim(geometryParameter)
  if (scan(getGeometryFile,'/') /= 1) getGeometryFile = getCWD()//'/'//trim(getGeometryFile)
  getGeometryFile = trim(makeRelativePath(getCWD(), getGeometryFile))

  inquire(file=getGeometryFile, exist=file_exists)
  if (.not. file_exists) then
    write(6,'(/,a)') ' ERROR: Geometry file does not exists ('//trim(getGeometryFile)//')'
    call quit(1)
  endif

end function getGeometryFile


!--------------------------------------------------------------------------------------------------
!> @brief relative path of load case from command line arguments
!--------------------------------------------------------------------------------------------------
function getLoadCaseFile(loadCaseParameter)

  character(len=:), allocatable :: getLoadCaseFile
  character(len=*),  intent(in) :: loadCaseParameter
  logical                       :: file_exists
  external                      :: quit

  getLoadCaseFile = trim(loadCaseParameter)
  if (scan(getLoadCaseFile,'/') /= 1) getLoadCaseFile = getCWD()//'/'//trim(getLoadCaseFile)
  getLoadCaseFile = trim(makeRelativePath(getCWD(), getLoadCaseFile))

  inquire(file=getLoadCaseFile, exist=file_exists)
  if (.not. file_exists) then
    write(6,'(/,a)') ' ERROR: Load case file does not exists ('//trim(getLoadCaseFile)//')'
    call quit(1)
  endif

end function getLoadCaseFile


!--------------------------------------------------------------------------------------------------
!> @brief remove ../, /./, and // from path.
!> @details works only if absolute path is given
!--------------------------------------------------------------------------------------------------
function rectifyPath(path)

  character(len=*), intent(in)  :: path
  character(len=:), allocatable :: rectifyPath
  integer :: i,j,k,l

!--------------------------------------------------------------------------------------------------
! remove /./ from path
  rectifyPath = trim(path)
  l = len_trim(rectifyPath)
  do i = l,3,-1
    if (rectifyPath(i-2:i) == '/./') rectifyPath(i-1:l) = rectifyPath(i+1:l)//'  '
  enddo

!--------------------------------------------------------------------------------------------------
! remove // from path
  l = len_trim(rectifyPath)
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
     if (rectifyPath(j+1:j+1) == '/') then                                                          !search for '//' that appear in case of XXX/../../XXX
       k = len_trim(rectifyPath)
       rectifyPath(j+1:k-1) = rectifyPath(j+2:k)
       rectifyPath(k:k) = ' '
     endif
     i = j+index(rectifyPath(j+1:l),'../')
  enddo
  if(len_trim(rectifyPath) == 0) rectifyPath = '/'

  rectifyPath = trim(rectifyPath)

end function rectifyPath


!--------------------------------------------------------------------------------------------------
!> @brief relative path from absolute a to absolute b
!--------------------------------------------------------------------------------------------------
function makeRelativePath(a,b)

  character (len=*), intent(in) :: a,b
  character (len=pPathLen)      :: a_cleaned,b_cleaned
  character(len=:), allocatable :: makeRelativePath
  integer :: i,posLastCommonSlash,remainingSlashes

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
!> @brief Set global variable interface_SIGTERM to .true.
!> @details This function can be registered to catch signals send to the executable.
!--------------------------------------------------------------------------------------------------
subroutine catchSIGTERM(signal) bind(C)

  integer(C_INT), value :: signal
  interface_SIGTERM = .true.

  write(6,'(a,i2.2,a)') ' received signal ',signal, ', set SIGTERM=TRUE'

end subroutine catchSIGTERM


!--------------------------------------------------------------------------------------------------
!> @brief Set global variable interface_SIGTERM.
!--------------------------------------------------------------------------------------------------
subroutine interface_setSIGTERM(state)

  logical, intent(in) :: state
  interface_SIGTERM = state

end subroutine interface_setSIGTERM


!--------------------------------------------------------------------------------------------------
!> @brief Set global variable interface_SIGUSR1 to .true.
!> @details This function can be registered to catch signals send to the executable.
!--------------------------------------------------------------------------------------------------
subroutine catchSIGUSR1(signal) bind(C)

  integer(C_INT), value :: signal
  interface_SIGUSR1 = .true.

  write(6,'(a,i2.2,a)') ' received signal ',signal, ', set SIGUSR1=TRUE'

end subroutine catchSIGUSR1


!--------------------------------------------------------------------------------------------------
!> @brief Set global variable interface_SIGUSR.
!--------------------------------------------------------------------------------------------------
subroutine interface_setSIGUSR1(state)

  logical, intent(in) :: state
  interface_SIGUSR1 = state

end subroutine interface_setSIGUSR1


!--------------------------------------------------------------------------------------------------
!> @brief Set global variable interface_SIGUSR2 to .true.
!> @details This function can be registered to catch signals send to the executable.
!--------------------------------------------------------------------------------------------------
subroutine catchSIGUSR2(signal) bind(C)

  integer(C_INT), value :: signal
  interface_SIGUSR2 = .true.

  write(6,'(a,i2.2,a)') ' received signal ',signal, ', set SIGUSR2=TRUE'

end subroutine catchSIGUSR2


!--------------------------------------------------------------------------------------------------
!> @brief Set global variable interface_SIGUSR2.
!--------------------------------------------------------------------------------------------------
subroutine interface_setSIGUSR2(state)

  logical, intent(in) :: state
  interface_SIGUSR2 = state

end subroutine interface_setSIGUSR2


end module
