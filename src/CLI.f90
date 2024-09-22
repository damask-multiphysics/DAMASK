!--------------------------------------------------------------------------------------------------
!> @author Jaeyong Jung, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Parse command line interface for PETSc-based solvers
!--------------------------------------------------------------------------------------------------
#define PETSC_MINOR_MIN 12
#define PETSC_MINOR_MAX 22

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
    CLI_restartInc = 0                                                                              !< increment at which calculation starts
  character(len=:), allocatable, public, protected :: &
    CLI_geomFile, &                                                                                 !< location of the geometry file
    CLI_loadFile, &                                                                                 !< location of the load case file
    CLI_materialFile, &                                                                             !< location of the material configuration file
    CLI_numericsFile, &                                                                             !< location of the numerics configuration file
    solverJobname

  public :: &
    getSolverJobname, &
    CLI_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief Initialize the solver by interpreting the command line arguments. Write
!! information on computation to screen
!--------------------------------------------------------------------------------------------------
subroutine CLI_init()
#include <petsc/finclude/petscsys.h>

#if PETSC_VERSION_MAJOR!=3 || PETSC_VERSION_MINOR<PETSC_MINOR_MIN || PETSC_VERSION_MINOR>PETSC_MINOR_MAX
--  UNSUPPORTED PETSc VERSION --- UNSUPPORTED PETSc VERSION --- UNSUPPORTED PETSc VERSION ---
#endif
#if    PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==18
#define PETSC_DOI '10.2172/1893326'
#elif  PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==19
#define PETSC_DOI '10.2172/1968587'
#elif  PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==20
#define PETSC_DOI '10.2172/2205494'
#elif  PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==21
#define PETSC_DOI '10.2172/2337606'
#endif
  character(len=:), allocatable :: &
    commandLine, &                                                                                  !< command line call as string
    flag, &                                                                                         !< individual flag
    val, &
    geomArg, &                                                                                      !< -g CLI argument
    loadArg, &                                                                                      !< -l CLI argument
    materialArg, &                                                                                  !< -m CLI argument
    numericsArg, &                                                                                  !< -n CLI argument
    workingDirArg                                                                                   !< -w CLI argument
  integer :: &
    i, s
  integer, dimension(8) :: &
    dateAndTime
#ifdef PETSC_DOI
  character(len=*), parameter :: PETSc_DOI = PETSC_DOI
#endif

  workingDirArg = getCWD()

  print'(/,1x,a)', '<<<+-  CLI init  -+>>>'

 ! http://patorjk.com/software/taag/#p=display&f=Lean&t=DAMASK%203
#ifdef DEBUG
  print'(a)', IO_color([255,0,0])
  print'(1x,a)', 'debug version - debug version - debug version - debug version - debug version'
#endif
  print'(a)', IO_color([67,128,208])
  print'(1x,a)', '    _/_/_/      _/_/    _/      _/    _/_/      _/_/_/  _/    _/    _/_/_/'
  print'(1x,a)', '   _/    _/  _/    _/  _/_/  _/_/  _/    _/  _/        _/  _/            _/'
  print'(1x,a)', '  _/    _/  _/_/_/_/  _/  _/  _/  _/_/_/_/    _/_/    _/_/          _/_/'
  print'(1x,a)', ' _/    _/  _/    _/  _/      _/  _/    _/        _/  _/  _/            _/'
  print'(1x,a)', '_/_/_/    _/    _/  _/      _/  _/    _/  _/_/_/    _/    _/    _/_/_/'
#if   defined(GRID)
  print'(a)', IO_color([123,207,68])
  print'(1x,a)', 'Grid solver'
#elif defined(MESH)
  print'(a)', IO_color([230,150,68])
  print'(1x,a)', 'Mesh solver'
#endif
#ifdef DEBUG
  print'(a)', IO_color([255,0,0])
  print'(1x,a)', 'debug version - debug version - debug version - debug version - debug version'
#endif
  print'(a)', IO_color()

  print'(1x,a)', 'F. Roters et al., Computational Materials Science 158:420–478, 2019'
  print'(1x,a)', 'https://doi.org/10.1016/j.commatsci.2018.04.030'//IO_EOL
#if PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=18
  print'(1x,a,i0,a,i0)', 'S. Balay et al., PETSc/TAO User Manual Revision ',PETSC_VERSION_MAJOR,'.',PETSC_VERSION_MINOR
#ifdef PETSC_DOI
  print'(1x,a)', 'https://doi.org/'//PETSc_DOI
#endif
#endif
  print'(/,1x,a)', 'Version: '//DAMASKVERSION

  print'(/,1x,a)', 'Compiled with: '//compiler_version()
  print'(1x,a)',   'Compiled on: '//CMAKE_SYSTEM
  print'(1x,a)',   'Compiler options: '//compiler_options()

  ! https://github.com/jeffhammond/HPCInfo/blob/master/docs/Preprocessor-Macros.md
  print'(/,1x,a)', 'Compiled on: '//__DATE__//' at '//__TIME__

  print'(/,1x,a,1x,i0,a,i0,a,i0)', &
                'PETSc version:',PETSC_VERSION_MAJOR,'.',PETSC_VERSION_MINOR,'.',PETSC_VERSION_SUBMINOR

  call date_and_time(values = dateAndTime)
  print'(/,1x,a,1x,2(i2.2,a),i4.4)', 'Date:',dateAndTime(3),'/',dateAndTime(2),'/',dateAndTime(1)
  print'(1x,a,1x,2(i2.2,a),i2.2)',   'Time:',dateAndTime(5),':',dateAndTime(6),':',dateAndTime(7)


  if (command_argument_count() == 0) call help()
  do i = 1, command_argument_count()
    flag = getArg(i)
    if (flag == '-h' .or. flag == '--help') call help()
  end do

  i = 1
  do while (i <= command_argument_count())
    flag = getArg(i)
    i = i + 1
    s = scan(flag,'=')
    if (s /= 0) then
      val = flag(s+1:)
      flag = flag(:s-1)
    else
      if (i > command_argument_count()) call IO_error(610,ext_msg=flag)
      val = getArg(i)
      i = i + 1
    end if

    select case(flag)
      case ('-g', '--geom', '--geometry')
        geomArg = val
      case ('-l', '--load', '--loadcase')
        loadArg = val
      case ('-m', '--material', '--materialconfig')
        materialArg = val
      case ('-n', '--numerics', '--numericsconfig')
        numericsArg = val
      case ('-j', '--job', '--jobname')
        solverJobname = val
      case ('-w', '--wd', '--workingdir', '--workingdirectory')
        workingDirArg = val
#if defined(GRID)
      case ('-r', '--rs', '--restart')
        CLI_restartInc = IO_strAsInt(val)
        if (CLI_restartInc < 0) call IO_error(611,ext_msg=val,label1='--restart')
#endif
      case default
        call IO_error(613,ext_msg=flag)
    end select
  end do

  if (.not. allocated(geomArg))     call IO_error(612,ext_msg='--geom')
  if (.not. allocated(loadArg))     call IO_error(612,ext_msg='--load')
  if (.not. allocated(materialArg)) call IO_error(612,ext_msg='--material')

  call setWorkingDirectory(trim(workingDirArg))
  CLI_geomFile = getPathRelCWD(geomArg,'geometry')
  CLI_loadFile = getPathRelCWD(loadArg,'load case')
  CLI_materialFile = getPathRelCWD(materialArg,'material configuration')
  if (allocated(numericsArg)) &
    CLI_numericsFile = getPathRelCWD(numericsArg,'numerics configuration')

  if (.not. allocated(solverJobname)) then
    solverJobname = jobname(CLI_geomFile,CLI_loadFile,CLI_materialFile,CLI_numericsFile)
  elseif (scan(solverJobname,'/') > 0) then
    call IO_error(612,ext_msg=solverJobname,label1='--jobname')
  endif

  commandLine = getArg(-1)

  print'(/,1x,a)',      'Host name: '//getHostName()
  print'(1x,a)',        'User name: '//getUserName()

  print'(/,1x,a,/)',    'Command line call:      '//trim(commandLine)
  print'(1x,a)',        'Working directory:      '//IO_glueDiffering(getCWD(),workingDirArg)
  print'(1x,a)',        'Geometry:               '//IO_glueDiffering(CLI_geomFile,geomArg)
  print'(1x,a)',        'Load case:              '//IO_glueDiffering(CLI_loadFile,loadArg)
  print'(1x,a)',        'Material config:        '//IO_glueDiffering(CLI_materialFile,materialArg)
  if (allocated(numericsArg)) &
    print'(1x,a)',      'Numerics config:        '//IO_glueDiffering(CLI_numericsFile,numericsArg)
  print'(1x,a)',        'Solver job name:        '//getSolverJobname()
  if (CLI_restartInc > 0) &
    print'(1x,a,i6.6)', 'Restart from increment: ', CLI_restartInc


end subroutine CLI_init


!--------------------------------------------------------------------------------------------------
!> @brief Get argument from command line.
!--------------------------------------------------------------------------------------------------
function getArg(n)

  integer, intent(in) :: n                                                                          !< number of the argument
  character(len=:), allocatable :: getArg

  integer :: l,err


  allocate(character(len=0)::getArg)
  if (n<0) then
    call get_command(getArg, length=l)
  else
    call get_command_argument(n,getArg,length=l)
  endif
  deallocate(getArg)
  allocate(character(len=l)::getArg)
  if (n<0) then
    call get_command(getArg, status=err)
  else
    call get_command_argument(n,getArg,status=err)
  endif
  if (err /= 0) error stop 'getting command arguments failed'

end function getArg


!--------------------------------------------------------------------------------------------------
!> @brief Extract working directory from given argument or from location of geometry file,
!!        possibly converting relative arguments to absolut path.
!--------------------------------------------------------------------------------------------------
subroutine setWorkingDirectory(workingDirectoryArg)

  character(len=*), intent(in)  :: workingDirectoryArg                                              !< working directory argument
  character(len=:), allocatable :: workingDirectory


  absolutePath: if (workingDirectoryArg(1:1) == '/') then
    workingDirectory = workingDirectoryArg
  else absolutePath
    workingDirectory = getCWD()
    workingDirectory = trim(workingDirectory)//'/'//workingDirectoryArg
  end if absolutePath

  workingDirectory = trim(normpath(workingDirectory))
  if (setCWD(trim(workingDirectory))) call IO_error(640,ext_msg=workingDirectory)

end subroutine setWorkingDirectory


!--------------------------------------------------------------------------------------------------
!> @brief Return solver job name (MSC.Marc compatible).
!--------------------------------------------------------------------------------------------------
function getSolverJobname()

  character(len=:), allocatable :: getSolverJobname


  getSolverJobname = solverJobname

end function getSolverJobname


!--------------------------------------------------------------------------------------------------
!> @brief Determine solver job name.
!--------------------------------------------------------------------------------------------------
function jobname(geomFile,LoadFile,materialsFile,numericsFile)

  character(len=:), allocatable :: jobname
  character(len=*), intent(in)  :: geomFile,loadFile,materialsFile
  character(len=:), allocatable, intent(in) :: numericsFile


  jobname = stem(geomFile)//'_'//stem(loadFile)//'_'//stem(materialsFile)
  if (allocated(numericsFile)) jobname = jobname//'_'//stem(numericsFile)

  contains

  function stem(fullname)

    character(len=:), allocatable :: stem
    character(len=*), intent(in)  :: fullname


    stem = fullname(scan(fullname,'/',back=.true.)+1:scan(fullname,'.',back=.true.)-1)

  end function stem

end function jobname


!--------------------------------------------------------------------------------------------------
!> @brief Translate path as relative to CWD and check for existence.
!--------------------------------------------------------------------------------------------------
function getPathRelCWD(path,fileType)

  character(len=:), allocatable :: getPathRelCWD
  character(len=*),  intent(in) :: path
  character(len=*),  intent(in) :: fileType

  logical                       :: file_exists


  getPathRelCWD = trim(path)
  if (scan(getPathRelCWD,'/') /= 1) getPathRelCWD = getCWD()//'/'//trim(getPathRelCWD)
  getPathRelCWD = trim(relpath(getPathRelCWD,getCWD()))

  inquire(file=getPathRelCWD, exist=file_exists)
  if (.not. file_exists) call IO_error(100,ext_msg=fileType//' "'//trim(getPathRelCWD)//'"')

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


!--------------------------------------------------------------------------------------------------
!> @brief Print usage instructions to STDOUT and terminate program.
!--------------------------------------------------------------------------------------------------
subroutine help()

  external :: quit


  print'(/,1x,a)','#######################################################################'
  print'(1x,a)',  'DAMASK Command Line Interface:'
  print'(1x,a)',  'Düsseldorf Advanced Material Simulation Kit with PETSc-based solvers'
  print'(1x,a,/)','#######################################################################'
  print'(1x,a,/)','Valid command line flags:'
  print'(1x,a)',  '   --geom         (-g, --geometry)'
  print'(1x,a)',  '   --load         (-l, --loadcase)'
  print'(1x,a)',  '   --material     (-m, --materialconfig)'
  print'(1x,a)',  '   --numerics     (-n, --numericsconfig)'
  print'(1x,a)',  '   --jobname      (-j, --job)'
  print'(1x,a)',  '   --workingdir   (-w, --wd, --workingdirectory)'
#if defined(GRID)
  print'(1x,a)',  '   --restart      (-r, --rs)'
#endif
  print'(1x,a)',  '   --help         (-h)'
  print'(/,1x,a)','-----------------------------------------------------------------------'
  print'(1x,a)',  'Mandatory flags:'
  print'(/,1x,a)','  --geom GEOMFILE'
  print'(1x,a)',  '       specify the file path of the geometry definition'
  print'(/,1x,a)','  --load LOADFILE'
  print'(1x,a)',  '       specify the file path of the load case definition'
  print'(/,1x,a)','  --material MATERIALFILE'
  print'(1x,a)',  '       specify the file path of the material configuration'
  print'(/,1x,a)','-----------------------------------------------------------------------'
  print'(1x,a)',  'Optional flags:'
  print'(/,1x,a)','  --numerics NUMERICSFILE'
  print'(1x,a)',  '       Specify the file path of the numerics configuration'
  print'(/,1x,a)','  --jobname JOBNAME'
  print'(1x,a)',  '       specify the job name.'
  print'(1x,a)',  '       Defaults to GEOM_LOAD_MATERIAL[_NUMERICS].'
  print'(/,1x,a)','  --workingdirectory WORKINGDIRECTORY'
  print'(1x,a)',  '       specify the base directory of relative paths.'
  print'(1x,a)',  '       Defaults to the current working directory'
#if defined(GRID)
  print'(/,1x,a)','  --restart N'
  print'(1x,a)',  '       read in increment N and continues with calculating'
  print'(1x,a)',  '           increment N+1, N+2, ... based on this'
  print'(1x,a)',  '       works only if the restart information for increment N'
  print'(1x,a)',  '           is available in JOBNAME_restart.hdf5'
  print'(1x,a)',  '       append to existing results file JOBNAME.hdf5'
#endif
  print'(/,1x,a)','-----------------------------------------------------------------------'
  print'(1x,a)',  'Help:'
  print'(/,1x,a)','  --help'
  print'(1x,a,/)','       Prints this message and exits'

  call quit(0)

end subroutine help

end module CLI
