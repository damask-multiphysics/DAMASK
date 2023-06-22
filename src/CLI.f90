!--------------------------------------------------------------------------------------------------
!> @author Jaeyong Jung, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Parse command line interface for PETSc-based solvers
!--------------------------------------------------------------------------------------------------
#define PETSC_MINOR_MIN 12
#define PETSC_MINOR_MAX 19

module CLI
  use, intrinsic :: ISO_fortran_env

  use PETScSys

  use prec
  use parallelization
  use system_routines
  use IO

  implicit none(type,external)
  private
  integer,                       public, protected :: &
    CLI_restartInc = 0                                                                              !< Increment at which calculation starts
  character(len=:), allocatable, public, protected :: &
    CLI_geomFile, &                                                                                 !< parameter given for geometry file
    CLI_loadFile, &                                                                                 !< parameter given for load case file
    CLI_materialFile

  public :: &
    getSolverJobName, &
    CLI_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief initializes the solver by interpreting the command line arguments. Also writes
!! information on computation to screen
!--------------------------------------------------------------------------------------------------
subroutine CLI_init()
#include <petsc/finclude/petscsys.h>

#if PETSC_VERSION_MAJOR!=3 || PETSC_VERSION_MINOR<PETSC_MINOR_MIN || PETSC_VERSION_MINOR>PETSC_MINOR_MAX
--  UNSUPPORTED PETSc VERSION --- UNSUPPORTED PETSc VERSION --- UNSUPPORTED PETSc VERSION ---
#endif

  character(len=:), allocatable :: &
    commandLine, &                                                                                  !< command line call as string
    arg, &                                                                                          !< individual argument
    loadCaseArg, &                                                                                  !< -l argument given to the executable
    geometryArg, &                                                                                  !< -g argument given to the executable
    materialArg, &                                                                                  !< -m argument given to the executable
    workingDirArg                                                                                   !< -w argument given to the executable
  integer :: &
    stat, &
    i
  integer, dimension(8) :: &
    dateAndTime
  integer        :: err
  external :: &
    quit


  workingDirArg = getCWD()

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
    arg = getArg(i)
    select case(trim(arg))                                                                          ! extract key
      case ('-h','--help')
        print'(/,a)',' #######################################################################'
        print'(a)',  ' DAMASK Command Line Interface:'
        print'(a)',  ' Düsseldorf Advanced Material Simulation Kit with PETSc-based solvers'
        print'(a,/)',' #######################################################################'
        print'(a,/)',' Valid command line switches:'
        print'(a)',  '    --geom         (-g, --geometry)'
        print'(a)',  '    --load         (-l, --loadcase)'
        print'(a)',  '    --material     (-m, --materialconfig)'
        print'(a)',  '    --workingdir   (-w, --wd, --workingdirectory)'
        print'(a)',  '    --restart      (-r, --rs)'
        print'(a)',  '    --help         (-h)'
        print'(/,a)',' -----------------------------------------------------------------------'
        print'(a)',  ' Mandatory arguments:'
        print'(/,a)','   --geom PathToGeomFile/NameOfGeom'
        print'(a)',  '        Specifies the location of the geometry definition file.'
        print'(/,a)','   --load PathToLoadFile/NameOfLoadFile'
        print'(a)',  '        Specifies the location of the load case definition file.'
        print'(/,a)','   --material PathToMaterialConfigurationFile/NameOfMaterialConfigurationFile'
        print'(a)',  '        Specifies the location of the material configuration file.'
        print'(/,a)',' -----------------------------------------------------------------------'
        print'(a)',  ' Optional arguments:'
        print'(/,a)','   --workingdirectory PathToWorkingDirectory'
        print'(a)',  '        Specifies the working directory and overwrites the default ./'
        print'(a)',  '        Make sure the file "material.yaml" exists in the working'
        print'(a)',  '            directory.'
        print'(a)',  '        For further configuration place "numerics.yaml"'
        print'(a)','              in that directory.'
        print'(/,a)','   --restart N'
        print'(a)',  '        Reads in increment N and continues with calculating'
        print'(a)',  '            increment N+1, N+2, ... based on this.'
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
        loadCaseArg = getArg(i+1)
      case ('-g', '--geom', '--geometry')
        geometryArg = getArg(i+1)
      case ('-m', '--material', '--materialconfig')
        materialArg = getArg(i+1)
      case ('-w', '--wd', '--workingdir', '--workingdirectory')
        workingDirArg = getArg(i+1)
      case ('-r', '--rs', '--restart')
        arg = getArg(i+1)
        read(arg,*,iostat=stat) CLI_restartInc
        if (CLI_restartInc < 0 .or. stat /=0) then
          print'(/,a)', ' ERROR: Could not parse restart increment: '//trim(arg)
          call quit(1)
        end if
    end select
    if (err /= 0) call quit(1)
  end do

  if (.not. all([allocated(loadcaseArg),allocated(geometryArg),allocated(materialArg)])) then
    print'(/,a)', ' ERROR: Please specify geometry AND load case AND material configuration (-h for help)'
    call quit(1)
  end if

  call setWorkingDirectory(trim(workingDirArg))
  CLI_geomFile = getPathRelCWD(geometryArg,'geometry')
  CLI_loadFile = getPathRelCWD(loadCaseArg,'load case')
  CLI_materialFile = getPathRelCWD(materialArg,'material configuration')

  commandLine = getArg(0)
  print'(/,a)',      ' Host name: '//getHostName()
  print'(a)',        ' User name: '//getUserName()

  print'(/a/)',      ' Command line call:      '//trim(commandLine)
  print'(a)',        ' Working directory:      '//IO_glueDiffering(getCWD(),workingDirArg)
  print'(a)',        ' Geometry:               '//IO_glueDiffering(CLI_geomFile,geometryArg)
  print'(a)',        ' Load case:              '//IO_glueDiffering(CLI_loadFile,loadCaseArg)
  print'(a)',        ' Material config:        '//IO_glueDiffering(CLI_materialFile,materialArg)
  print'(a)',        ' Solver job name:        '//getSolverJobName()
  if (CLI_restartInc > 0) &
    print'(a,i6.6)', ' Restart from increment: ', CLI_restartInc

  contains

  !------------------------------------------------------------------------------------------------
  !> @brief Get argument from command line.
  !------------------------------------------------------------------------------------------------
  function getArg(n)

    integer, intent(in) :: n                                                                        !< number of the argument
    character(len=:), allocatable :: getArg

    integer :: l,err


    allocate(character(len=0)::getArg)
    call get_command_argument(n,getArg,length=l)
    deallocate(getArg)
    allocate(character(len=l)::getArg)
    call get_command_argument(n,getArg,status=err)
    if (err /= 0) call quit(1)

  end function getArg

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

  workingDirectory = trim(normpath(workingDirectory))
  error = setCWD(trim(workingDirectory))
  if (error) then
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
!> @brief Translate path as relative to CWD and check for existence.
!--------------------------------------------------------------------------------------------------
function getPathRelCWD(path,fileType)

  character(len=:), allocatable :: getPathRelCWD
  character(len=*),  intent(in) :: path
  character(len=*),  intent(in) :: fileType

  logical                       :: file_exists
  external                      :: quit


  getPathRelCWD = trim(path)
  if (scan(getPathRelCWD,'/') /= 1) getPathRelCWD = getCWD()//'/'//trim(getPathRelCWD)
  getPathRelCWD = trim(relpath(getPathRelCWD,getCWD()))

  inquire(file=getPathRelCWD, exist=file_exists)
  if (.not. file_exists) then
    print*, 'ERROR: '//fileType//' file does not exist: '//trim(getPathRelCWD)
    call quit(1)
  end if

end function getPathRelCWD


!--------------------------------------------------------------------------------------------------
!> @brief Remove ../, /./, and // from path.
!> @details Works only if absolute path is given.
!--------------------------------------------------------------------------------------------------
function normpath(path)

  character(len=*), intent(in)  :: path
  character(len=:), allocatable :: normpath

  integer :: i,j,k,l


!--------------------------------------------------------------------------------------------------
! remove /./ from path
  normpath = trim(path)
  l = len_trim(normpath)
  do i = l,3,-1
    if (normpath(i-2:i) == '/./') normpath(i-1:l) = normpath(i+1:l)//'  '
  end do

!--------------------------------------------------------------------------------------------------
! remove // from path
  l = len_trim(normpath)
  do i = l,2,-1
    if (normpath(i-1:i) == '//') normpath(i-1:l) = normpath(i:l)//' '
  end do

!--------------------------------------------------------------------------------------------------
! remove ../ and corresponding directory from path
  l = len_trim(normpath)
  i = index(normpath(i:l),'../')
  j = 0
  do while (i > j)
     j = scan(normpath(1:i-2),'/',back=.true.)
     normpath(j+1:l) = normpath(i+3:l)//repeat(' ',2+i-j)
     if (normpath(j+1:j+1) == '/') then                                                             !search for '//' that appear in case of XXX/../../XXX
       k = len_trim(normpath)
       normpath(j+1:k-1) = normpath(j+2:k)
       normpath(k:k) = ' '
     end if
     i = j+index(normpath(j+1:l),'../')
  end do
  if (len_trim(normpath) == 0) normpath = '/'

  normpath = trim(normpath)

end function normpath


!--------------------------------------------------------------------------------------------------
!> @brief Determine relative path.
!--------------------------------------------------------------------------------------------------
function relpath(path,start)

  character(len=*), intent(in)  :: start,path
  character(len=:), allocatable :: relpath

  character(len=:), allocatable :: start_cleaned,path_cleaned
  integer :: i,posLastCommonSlash,remainingSlashes


  posLastCommonSlash = 0
  remainingSlashes = 0
  start_cleaned = normpath(trim(start)//'/')
  path_cleaned = normpath(path)

  do i = 1, min(len_trim(start_cleaned),len_trim(path_cleaned))
    if (start_cleaned(i:i) /= path_cleaned(i:i)) exit
    if (start_cleaned(i:i) == '/') posLastCommonSlash = i
  end do
  do i = posLastCommonSlash+1,len_trim(start_cleaned)
    if (start_cleaned(i:i) == '/') remainingSlashes = remainingSlashes + 1
  end do

  relpath = repeat('..'//'/',remainingSlashes)//path_cleaned(posLastCommonSlash+1:len_trim(path_cleaned))

end function relpath

end module CLI
