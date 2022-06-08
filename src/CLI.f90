!--------------------------------------------------------------------------------------------------
!> @author Jaeyong Jung, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Parse command line interface for PETSc-based solvers
!--------------------------------------------------------------------------------------------------
#define PETSC_MAJOR 3
#define PETSC_MINOR_MIN 12
#define PETSC_MINOR_MAX 17

module CLI
  use, intrinsic :: ISO_fortran_env

  use PETScSys

  use prec
  use parallelization
  use system_routines

  implicit none
  private
  integer,                       public, protected :: &
    CLI_restartInc = 0                                                                              !< Increment at which calculation starts
  character(len=:), allocatable, public, protected :: &
    CLI_geomFile, &                                                                                 !< parameter given for geometry file
    CLI_loadFile                                                                                    !< parameter given for load case file

  public :: &
    getSolverJobName, &
    CLI_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief initializes the solver by interpreting the command line arguments. Also writes
!! information on computation to screen
!--------------------------------------------------------------------------------------------------
subroutine CLI_init
#include <petsc/finclude/petscsys.h>

#if PETSC_VERSION_MAJOR!=3 || PETSC_VERSION_MINOR<PETSC_MINOR_MIN || PETSC_VERSION_MINOR>PETSC_MINOR_MAX
--  UNSUPPORTED PETSc VERSION --- UNSUPPORTED PETSc VERSION --- UNSUPPORTED PETSc VERSION ---
#endif

  character(len=pPathLen*3+pStringLen) :: &
    commandLine                                                                                     !< command line call as string
  character(len=pPathLen) :: &
    arg, &                                                                                          !< individual argument
    loadCaseArg   = '', &                                                                           !< -l argument given to the executable
    geometryArg   = '', &                                                                           !< -g argument given to the executable
    workingDirArg = ''                                                                              !< -w argument given to the executable
  integer :: &
    stat, &
    i
  integer, dimension(8) :: &
    dateAndTime
  integer        :: err
  external :: &
    quit


  print'(/,1x,a)', '<<<+-  CLI init  -+>>>'

 ! http://patorjk.com/software/taag/#p=display&f=Lean&t=DAMASK%203
#ifdef DEBUG
  print*, achar(27)//'[31m'
  print'(a,/)', ' debug version - debug version - debug version - debug version - debug version'
#else
  print*, achar(27)//'[94m'
#endif
  print*, '     _/_/_/      _/_/    _/      _/    _/_/      _/_/_/  _/    _/    _/_/_/'
  print*, '    _/    _/  _/    _/  _/_/  _/_/  _/    _/  _/        _/  _/            _/'
  print*, '   _/    _/  _/_/_/_/  _/  _/  _/  _/_/_/_/    _/_/    _/_/          _/_/'
  print*, '  _/    _/  _/    _/  _/      _/  _/    _/        _/  _/  _/            _/'
  print*, ' _/_/_/    _/    _/  _/      _/  _/    _/  _/_/_/    _/    _/    _/_/_/'
#if   defined(GRID)
  print*, ' Grid solver'
#elif defined(MESH)
  print*, ' Mesh solver'
#endif
#ifdef DEBUG
  print'(/,a)', ' debug version - debug version - debug version - debug version - debug version'
#endif
  print*, achar(27)//'[0m'

  print*, 'F. Roters et al., Computational Materials Science 158:420–478, 2019'
  print*, 'https://doi.org/10.1016/j.commatsci.2018.04.030'

  print'(/,a)', ' Version: '//DAMASKVERSION

  print'(/,a)', ' Compiled with: '//compiler_version()
  print'(a)',   ' Compiled on: '//CMAKE_SYSTEM
  print'(a)',   ' Compiler options: '//compiler_options()

  ! https://github.com/jeffhammond/HPCInfo/blob/master/docs/Preprocessor-Macros.md
  print'(/,a)', ' Compiled on: '//__DATE__//' at '//__TIME__

  print'(/,a,i0,a,i0,a,i0)', &
                ' PETSc version: ',PETSC_VERSION_MAJOR,'.',PETSC_VERSION_MINOR,'.',PETSC_VERSION_SUBMINOR

  call date_and_time(values = dateAndTime)
  print'(/,a,2(i2.2,a),i4.4)', ' Date: ',dateAndTime(3),'/',dateAndTime(2),'/', dateAndTime(1)
  print'(a,2(i2.2,a),i2.2)',   ' Time: ',dateAndTime(5),':', dateAndTime(6),':', dateAndTime(7)

  do i = 1, command_argument_count()
    call get_command_argument(i,arg,status=err)
    if (err /= 0) call quit(1)
    select case(trim(arg))                                                                          ! extract key
      case ('-h','--help')
        print'(/,a)',' #######################################################################'
        print'(a)',  ' DAMASK Command Line Interface:'
        print'(a)',  ' Düsseldorf Advanced Material Simulation Kit with PETSc-based solvers'
        print'(a,/)',' #######################################################################'
        print'(a,/)',' Valid command line switches:'
        print'(a)',  '    --geom         (-g, --geometry)'
        print'(a)',  '    --load         (-l, --loadcase)'
        print'(a)',  '    --workingdir   (-w, --wd, --workingdirectory)'
        print'(a)',  '    --restart      (-r, --rs)'
        print'(a)',  '    --help         (-h)'
        print'(/,a)',' -----------------------------------------------------------------------'
        print'(a)',  ' Mandatory arguments:'
        print'(/,a)','   --geom PathToGeomFile/NameOfGeom'
        print'(a)',  '        Specifies the location of the geometry definition file.'
        print'(/,a)','   --load PathToLoadFile/NameOfLoadFile'
        print'(a)',  '        Specifies the location of the load case definition file.'
        print'(/,a)',' -----------------------------------------------------------------------'
        print'(a)',  ' Optional arguments:'
        print'(/,a)','   --workingdirectory PathToWorkingDirectory'
        print'(a)',  '        Specifies the working directory and overwrites the default ./'
        print'(a)',  '        Make sure the file "material.yaml" exists in the working'
        print'(a)',  '            directory.'
        print'(a)',  '        For further configuration place "numerics.yaml"'
        print'(a)','            and "debug.yaml" in that directory.'
        print'(/,a)','   --restart N'
        print'(a)',  '        Reads in increment N and continues with calculating'
        print'(a)',  '            increment N+1 based on this.'
        print'(a)',  '        Appends to existing results file'
        print'(a)',  '            "NameOfGeom_NameOfLoadFile.hdf5".'
        print'(a)',  '        Works only if the restart information for increment N'
        print'(a)',  '            is available in the working directory.'
        print'(/,a)',' -----------------------------------------------------------------------'
        print'(a)',  ' Help:'
        print'(/,a)','   --help'
        print'(a,/)','        Prints this message and exits'
        call quit(0)                                                                                ! normal Termination
      case ('-l', '--load', '--loadcase')
        call get_command_argument(i+1,loadCaseArg,status=err)
      case ('-g', '--geom', '--geometry')
        call get_command_argument(i+1,geometryArg,status=err)
      case ('-w', '--wd', '--workingdir', '--workingdirectory')
        call get_command_argument(i+1,workingDirArg,status=err)
      case ('-r', '--rs', '--restart')
        call get_command_argument(i+1,arg,status=err)
        read(arg,*,iostat=stat) CLI_restartInc
        if (CLI_restartInc < 0 .or. stat /=0) then
          print'(/,a)', ' ERROR: Could not parse restart increment: '//trim(arg)
          call quit(1)
        end if
    end select
    if (err /= 0) call quit(1)
  end do

  if (len_trim(loadcaseArg) == 0 .or. len_trim(geometryArg) == 0) then
    print'(/,a)', ' ERROR: Please specify geometry AND load case (-h for help)'
    call quit(1)
  end if

  if (len_trim(workingDirArg) > 0) call setWorkingDirectory(trim(workingDirArg))
  CLI_geomFile = getGeometryFile(geometryArg)
  CLI_loadFile = getLoadCaseFile(loadCaseArg)

  call get_command(commandLine)
  print'(/,a)',      ' Host name: '//getHostName()
  print'(a)',        ' User name: '//getUserName()

  print'(/a)',       ' Command line call:      '//trim(commandLine)
  if (len_trim(workingDirArg) > 0) &
    print'(a)',      ' Working dir argument:   '//trim(workingDirArg)
  print'(a)',        ' Geometry argument:      '//trim(geometryArg)
  print'(a)',        ' Load case argument:     '//trim(loadcaseArg)
  print'(/,a)',      ' Working directory:      '//getCWD()
  print'(a)',        ' Geometry file:          '//CLI_geomFile
  print'(a)',        ' Load case file:         '//CLI_loadFile
  print'(a)',        ' Solver job name:        '//getSolverJobName()
  if (CLI_restartInc > 0) &
    print'(a,i6.6)', ' Restart from increment: ', CLI_restartInc

end subroutine CLI_init


!--------------------------------------------------------------------------------------------------
!> @brief extract working directory from given argument or from location of geometry file,
!!        possibly converting relative arguments to absolut path
!--------------------------------------------------------------------------------------------------
subroutine setWorkingDirectory(workingDirectoryArg)

  character(len=*), intent(in)  :: workingDirectoryArg                                              !< working directory argument
  character(len=:), allocatable :: workingDirectory
  logical                       :: error
  external                      :: quit

  absolutePath: if (workingDirectoryArg(1:1) == '/') then
    workingDirectory = workingDirectoryArg
  else absolutePath
    workingDirectory = getCWD()
    workingDirectory = trim(workingDirectory)//'/'//workingDirectoryArg
  end if absolutePath

  workingDirectory = trim(rectifyPath(workingDirectory))
  error = setCWD(trim(workingDirectory))
  if(error) then
    print*, 'ERROR: Invalid Working directory: '//trim(workingDirectory)
    call quit(1)
  end if

end subroutine setWorkingDirectory


!--------------------------------------------------------------------------------------------------
!> @brief solver job name (no extension) as combination of geometry and load case name
!--------------------------------------------------------------------------------------------------
function getSolverJobName()

  character(len=:), allocatable :: getSolverJobName
  integer :: posExt,posSep

  posExt = scan(CLI_geomFile,'.',back=.true.)
  posSep = scan(CLI_geomFile,'/',back=.true.)

  getSolverJobName = CLI_geomFile(posSep+1:posExt-1)

  posExt = scan(CLI_loadFile,'.',back=.true.)
  posSep = scan(CLI_loadFile,'/',back=.true.)

  getSolverJobName = getSolverJobName//'_'//CLI_loadFile(posSep+1:posExt-1)

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
    print*, 'ERROR: Geometry file does not exists: '//trim(getGeometryFile)
    call quit(1)
  end if

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
    print*, 'ERROR: Load case file does not exists: '//trim(getLoadCaseFile)
    call quit(1)
  end if

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
  end do

!--------------------------------------------------------------------------------------------------
! remove // from path
  l = len_trim(rectifyPath)
  do i = l,2,-1
    if (rectifyPath(i-1:i) == '//') rectifyPath(i-1:l) = rectifyPath(i:l)//' '
  end do

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
     end if
     i = j+index(rectifyPath(j+1:l),'../')
  end do
  if(len_trim(rectifyPath) == 0) rectifyPath = '/'

  rectifyPath = trim(rectifyPath)

end function rectifyPath


!--------------------------------------------------------------------------------------------------
!> @brief Determine relative path from absolute a to absolute b
!--------------------------------------------------------------------------------------------------
function makeRelativePath(a,b)

  character(len=*), intent(in)  :: a,b
  character(len=pPathLen)       :: a_cleaned,b_cleaned
  character(len=:), allocatable :: makeRelativePath
  integer :: i,posLastCommonSlash,remainingSlashes

  posLastCommonSlash = 0
  remainingSlashes = 0
  a_cleaned = rectifyPath(trim(a)//'/')
  b_cleaned = rectifyPath(b)

  do i = 1, min(len_trim(a_cleaned),len_trim(rectifyPath(b_cleaned)))
    if (a_cleaned(i:i) /= b_cleaned(i:i)) exit
    if (a_cleaned(i:i) == '/') posLastCommonSlash = i
  end do
  do i = posLastCommonSlash+1,len_trim(a_cleaned)
    if (a_cleaned(i:i) == '/') remainingSlashes = remainingSlashes + 1
  end do

  makeRelativePath = repeat('..'//'/',remainingSlashes)//b_cleaned(posLastCommonSlash+1:len_trim(b_cleaned))

end function makeRelativePath

end module CLI
