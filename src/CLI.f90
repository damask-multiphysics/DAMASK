! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @author Jaeyong Jung, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Parse command line interface for PETSc-based solvers
!--------------------------------------------------------------------------------------------------
module CLI
  use, intrinsic :: ISO_Fortran_env
  use, intrinsic :: ISO_C_binding

  use PETScSys

  use prec
  use parallelization
  use OS
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
    CLI_jobName, &                                                                                  !< name of the job (will be used for DADF5 result file)
    CLI_jobID                                                                                       !< unique job ID (UUID)

#if (defined(BOOST) && !defined(OLD_STYLE_C_TO_FORTRAN_STRING))
   type :: tCLIBuffer
      character(len=:, kind=C_CHAR), pointer :: buf
   end type tCLIBuffer

   type :: tCLIArgs
      integer(C_INT)                              :: argc = 0
      type(C_PTR),      allocatable, dimension(:) :: argv
      type(tCLIBuffer), allocatable, dimension(:) :: cliBuffer
   contains
      procedure :: copyCommandLineArgs
   end type tCLIArgs

interface
  function C_CLI__new(argc, argv, worldrank) result(this) bind(C, name='CLI__new')
    use, intrinsic :: ISO_C_binding, only: C_INT, C_PTR
    integer(C_INT), intent(in) :: argc
    type(C_PTR),    intent(in) :: argv(*) ! MD I think this should be dimension(:), but wait for working ifx
    integer(C_INT), intent(in) :: worldrank
    type(C_PTR)                :: this
  end function C_CLI__new

  subroutine C_CLI_getParsedArgs(cli, geom, load, material, numerics, jobname, uuid, restart, stat) &
      bind(C, name='CLI_getParsedArgs')
    use ISO_C_binding
    type(C_PTR), value :: cli
    character(kind=C_CHAR,len=:), allocatable, intent(out) :: geom, load, material, numerics, jobname, uuid
    integer(C_INT), intent(out) :: restart, stat
  end subroutine C_CLI_getParsedArgs
end interface

#endif

public :: &
  CLI_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief Initialize the solver by interpreting the command line arguments and write
!!        information on computation to screen.
!--------------------------------------------------------------------------------------------------
subroutine CLI_init()
#include <petsc/finclude/petscsys.h>

#if PETSC_VERSION_MAJOR!=3 || PETSC_VERSION_MINOR<PETSC_MINOR_MIN || PETSC_VERSION_MINOR>PETSC_MINOR_MAX
--  UNSUPPORTED PETSc VERSION --- UNSUPPORTED PETSc VERSION --- UNSUPPORTED PETSc VERSION ---
#endif
#if   PETSC_VERSION_MINOR==19
#define PETSC_DOI '10.2172/1968587'
#elif PETSC_VERSION_MINOR==20
#define PETSC_DOI '10.2172/2205494'
#elif PETSC_VERSION_MINOR==21
#define PETSC_DOI '10.2172/2337606'
#elif PETSC_VERSION_MINOR==22
#define PETSC_DOI '10.2172/2476320'
#elif PETSC_VERSION_MINOR==23
#define PETSC_DOI '10.2172/2565610'
#elif PETSC_VERSION_MINOR==24
#define PETSC_DOI '10.2172/2998643'
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
#ifdef PETSC_DOI
  character(len=*), parameter :: PETSc_DOI = PETSC_DOI
#endif
#if (defined(BOOST) && !defined(OLD_STYLE_C_TO_FORTRAN_STRING))
  type(C_PTR) :: CLI_ = C_NULL_PTR
  type(tCLIArgs) :: cliArgs
  integer(C_INT) :: stat
#endif

  print'(/,1x,a)', '<<<+-  CLI init  -+>>>'

#if (defined(BOOST) && !defined(OLD_STYLE_C_TO_FORTRAN_STRING))
  print'(/,1x,a)', 'Using C++ parser'

  call cliArgs%copyCommandLineArgs()
  ! https://fortran-lang.discourse.group/t/c-interoperability-command-line-arguments/5773/7
  CLI_ = C_CLI__new(cliArgs%argc, cliArgs%argv, worldrank)
  call C_CLI_getParsedArgs(CLI_, CLI_geomFile, CLI_loadFile, CLI_materialFile, &
                           CLI_numericsFile, CLI_jobID, CLI_jobName, CLI_restartInc, stat)
  if (stat /= 0) error stop 'could not collect parsed args from CLI.cpp'
  call parallelization_bcast_str(CLI_jobID)
#else
  print'(/,1x,a)', 'Using Fortran parser'
  workingDirArg = OS_getCWD()
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

  print'(1x,a)',           'F. Roters et al., Computational Materials Science 158:420–478, 2019'
  print'(1x,a)',           'https://doi.org/10.1016/j.commatsci.2018.04.030'
  print'(/,1x,a,i0,a,i0)', 'S. Balay et al., PETSc/TAO User Manual Revision ',PETSC_VERSION_MAJOR,'.',PETSC_VERSION_MINOR
#ifdef PETSC_DOI
  print'(1x,a)',           'https://doi.org/'//PETSc_DOI
#endif
  print'(/,1x,a)', 'Version: '//DAMASK_VERSION

  ! https://github.com/jeffhammond/HPCInfo/blob/master/docs/Preprocessor-Macros.md
  print'(/,1x,a)', 'Compiled with: '//compiler_version()
  call printCompileOptions()
  print'(1x,a,1x,i0,a,i0,a,i0)', &
                   'PETSc version:',PETSC_VERSION_MAJOR,'.',PETSC_VERSION_MINOR,'.',PETSC_VERSION_SUBMINOR

  print'(/,1x,a)', 'Compiled at: '//__DATE__//' at '//__TIME__

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
        CLI_jobName = val
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

  if (.not. allocated(CLI_jobName)) then
    CLI_jobName = jobname(CLI_geomFile,CLI_loadFile,CLI_materialFile,CLI_numericsFile)
  elseif (scan(CLI_jobName,'/') > 0) then
    call IO_error(612,ext_msg=CLI_jobName,label1='--jobname')
  endif

  commandLine = getArg(-1)

  print'(/,1x,a)',      'Host name: '//OS_getHostName()
  print'(1x,a)',        'User name: '//OS_getUserName()

  print'(/,1x,a,/)',    'Command line call:  '//trim(commandLine)
  print'(1x,a)',        'Working directory:  '//IO_glueDiffering(OS_getCWD(),workingDirArg)
  print'(1x,a)',        'Geometry:           '//IO_glueDiffering(CLI_geomFile,geomArg)
  print'(1x,a)',        'Load case:          '//IO_glueDiffering(CLI_loadFile,loadArg)
  print'(1x,a)',        'Material config:    '//IO_glueDiffering(CLI_materialFile,materialArg)
  if (allocated(numericsArg)) &
    print'(1x,a)',      'Numerics config:    '//IO_glueDiffering(CLI_numericsFile,numericsArg)
  print'(1x,a)',        'Job name:           '//CLI_jobName
#if (defined(BOOST) && !defined(OLD_STYLE_C_TO_FORTRAN_STRING))
  print'(1x,a)',        'Job ID:             '//CLI_jobID
#endif
  if (CLI_restartInc > 0) &
    print'(1x,a,i0)',   'Restart increment:  ', CLI_restartInc
#endif

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
    workingDirectory = OS_getCWD()
    workingDirectory = trim(workingDirectory)//'/'//workingDirectoryArg
  end if absolutePath

  workingDirectory = trim(normpath(workingDirectory))
  if (OS_setCWD(trim(workingDirectory))) call IO_error(640,ext_msg=workingDirectory)

end subroutine setWorkingDirectory


!--------------------------------------------------------------------------------------------------
!> @brief Print fortran compiler options to stdout.
!--------------------------------------------------------------------------------------------------
subroutine printCompileOptions() bind(C, name="F_printCompileOptions")

  print'(1x,a)', 'Compiler options: '//compiler_options()
  print'(1x,a)', 'Compiled for: '//CMAKE_SYSTEM_NAME//' on '//CMAKE_SYSTEM_PROCESSOR

end subroutine printCompileOptions

#if (defined(BOOST) && !defined(OLD_STYLE_C_TO_FORTRAN_STRING))
!--------------------------------------------------------------------------------------------------
!> @brief Copy command line args to c strings for boost-processing,
!         allocate one extra element for the NULL termination.
!--------------------------------------------------------------------------------------------------
subroutine copyCommandLineArgs(this)

  class(tCLIArgs), intent(out) :: this

  integer :: n_args, i, len_arg


  n_args = command_argument_count() + 1
  this%argc = n_args
  allocate(this%argv(n_args))
  allocate(this%cliBuffer(n_args))

  do i = 0, n_args-1
    call get_command_argument(i, length=len_arg)
    allocate(character(len=len_arg+1, kind=C_CHAR) :: this%cliBuffer(i+1)%buf)
    if (len_arg > 0) &
      call get_command_argument(i, value=this%cliBuffer(i+1)%buf(1:len_arg))
    this%cliBuffer(i+1)%buf(len_arg+1:len_arg+1) = C_NULL_CHAR                                      !< need to use substring syntax because buf is defined as a fixed-length character
    this%argv(i+1) = c_loc(this%cliBuffer(i+1)%buf)
  end do

end subroutine copyCommandLineArgs
#endif

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

  pure function stem(fullname)
#ifndef __GFORTRAN__
    import, none
#endif
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
  if (scan(getPathRelCWD,'/') /= 1) getPathRelCWD = OS_getCWD()//'/'//trim(getPathRelCWD)
  getPathRelCWD = trim(relpath(getPathRelCWD,OS_getCWD()))

  inquire(file=getPathRelCWD, exist=file_exists)
  if (.not. file_exists) call IO_error(100_pI16, 'file does not exist',getPathRelCWD,emph=[2])

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
#if defined(GRID)
  print'(1x,a)',  '       Relative or absolute path to a VTK image data file (*.vti)'
  print'(1x,a)',  '       with mandatory "material" field variable.'
#elif defined(MESH)
  print'(1x,a)',  '       Relative or absolute path to a Gmsh file (*.msh)'
  print'(1x,a)',  '       with definitions of physical groups/tags for material IDs'
  print'(1x,a)',  '       and boundary conditions.'
#endif
  print'(/,1x,a)','  --load LOADFILE'
  print'(1x,a)',  '       Relative or absolute path to a load case definition'
  print'(1x,a)',  '       in YAML format.'
  print'(/,1x,a)','  --material MATERIALFILE'
  print'(1x,a)',  '       Relative or absolute path to a material configuration'
  print'(1x,a)',  '       in YAML format.'
  print'(/,1x,a)','-----------------------------------------------------------------------'
  print'(1x,a)',  'Optional flags:'
  print'(/,1x,a)','  --numerics NUMERICSFILE'
  print'(1x,a)',  '       Relative or absolute path to a numerics configuration'
  print'(1x,a)',  '       in YAML format.'
  print'(/,1x,a)','  --jobname JOBNAME'
  print'(1x,a)',  '       Job name, defaults to GEOM_LOAD_MATERIAL[_NUMERICS].'
  print'(/,1x,a)','  --workingdir WORKINGDIRECTORY'
  print'(1x,a)',  '       Working directory, defaults to current directory and'
  print'(1x,a)',  '       serves as base directory of relative paths.'
#if defined(GRID)
  print'(/,1x,a)','  --restart N'
  print'(1x,a)',  '       Restart simulation from given increment.'
  print'(1x,a)',  '       Read in increment N and, based on this, continue with'
  print'(1x,a)',  '       calculating increments N+1, N+2, ...'
  print'(1x,a)',  '       Requires restart information for increment N to be present in'
  print'(1x,a)',  '       JOBNAME_restart.hdf5 and will append subsequent results to'
  print'(1x,a)',  '       existing file JOBNAME.hdf5.'
#endif
  print'(/,1x,a)','-----------------------------------------------------------------------'
  print'(1x,a)',  'Help:'
  print'(/,1x,a)','  --help'
  print'(1x,a,/)','       Display help and exit.'

  call quit(0)

end subroutine help

end module CLI
