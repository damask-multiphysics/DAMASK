!--------------------------------------------------------------------------------------------------
!> @author   Jaeyong Jung, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author   Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief    Interfacing between the PETSc-based solvers and the material subroutines provided
!!           by DAMASK
!> @details  Interfacing between the PETSc-based solvers and the material subroutines provided
!>           by DAMASK. Interpretating the command line arguments to get load case, geometry file,
!>           and working directory.
!--------------------------------------------------------------------------------------------------
#define GCC_MIN 6
#define INTEL_MIN 1600
#define PETSC_MAJOR 3
#define PETSC_MINOR_MIN 10
#define PETSC_MINOR_MAX 11

module DAMASK_interface
  use, intrinsic :: iso_fortran_env
  use PETScSys

  use prec
  use system_routines
  
  implicit none
  private
  logical, public, protected :: &
    SIGTERM, &                                                                                      !< termination signal 
    SIGUSR1, &                                                                                      !< 1. user-defined signal
    SIGUSR2                                                                                         !< 2. user-defined signal
  integer, public, protected :: &
    interface_restartInc = 0                                                                        !< Increment at which calculation starts
  character(len=1024), public, protected :: &
    geometryFile = '', &                                                                            !< parameter given for geometry file
    loadCaseFile = ''                                                                               !< parameter given for load case file

  public :: &
    getSolverJobName, &
    DAMASK_interface_init, &
    setSIGTERM, &
    setSIGUSR1, &
    setSIGUSR2

contains

!--------------------------------------------------------------------------------------------------
!> @brief initializes the solver by interpreting the command line arguments. Also writes
!! information on computation to screen
!--------------------------------------------------------------------------------------------------
subroutine DAMASK_interface_init
#include <petsc/finclude/petscsys.h>
#if defined(__GFORTRAN__) &&  __GNUC__<GCC_MIN
===================================================================================================
   -----  WRONG COMPILER VERSION ----- WRONG COMPILER VERSION ----- WRONG COMPILER VERSION -----
===================================================================================================
===============   THIS VERSION OF DAMASK REQUIRES A NEWER gfortran VERSION   ======================
==================   THIS VERSION OF DAMASK REQUIRES A NEWER gfortran VERSION   ===================
=====================   THIS VERSION OF DAMASK REQUIRES A NEWER gfortran VERSION   ================
===================================================================================================
   -----  WRONG COMPILER VERSION ----- WRONG COMPILER VERSION ----- WRONG COMPILER VERSION -----
===================================================================================================
#endif

#if defined(__INTEL_COMPILER) && __INTEL_COMPILER<INTEL_MIN
===================================================================================================
   -----  WRONG COMPILER VERSION ----- WRONG COMPILER VERSION ----- WRONG COMPILER VERSION -----
===================================================================================================
=================   THIS VERSION OF DAMASK REQUIRES A NEWER ifort VERSION   =======================
====================   THIS VERSION OF DAMASK REQUIRES A NEWER ifort VERSION   ====================
=======================   THIS VERSION OF DAMASK REQUIRES A NEWER ifort VERSION   =================
===================================================================================================
   -----  WRONG COMPILER VERSION ----- WRONG COMPILER VERSION ----- WRONG COMPILER VERSION -----
===================================================================================================
#endif

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

  character(len=1024) :: &
    commandLine, &                                                                                  !< command line call as string
    arg, &                                                                                          !< individual argument
    loadCaseArg   = '', &                                                                           !< -l argument given to the executable
    geometryArg   = '', &                                                                           !< -g argument given to the executable
    workingDirArg = '', &                                                                           !< -w argument given to the executable
    userName                                                                                        !< name of user calling the executable
  integer :: &
    stat, &
    i, &
#ifdef _OPENMP
    threadLevel, &
#endif
    worldrank = 0, &
    worldsize = 0, &
    typeSize
  integer, dimension(8) :: &
    dateAndTime
  integer        :: mpi_err
  PetscErrorCode :: petsc_err
  external :: &
    quit

  open(6, encoding='UTF-8')                                                                         ! for special characters in output

!--------------------------------------------------------------------------------------------------
! PETSc Init
#ifdef _OPENMP
  ! If openMP is enabled, check if the MPI libary supports it and initialize accordingly.
  ! Otherwise, the first call to PETSc will do the initialization.
  call MPI_Init_Thread(MPI_THREAD_FUNNELED,threadLevel,mpi_err)
  if (mpi_err /= 0) call quit(1)
  if (threadLevel<MPI_THREAD_FUNNELED) then
    write(6,'(/,a)') ' ERROR: MPI library does not support OpenMP'
    call quit(1)
  endif
#endif
  call PETScInitializeNoArguments(petsc_err)                                                        ! according to PETSc manual, that should be the first line in the code
  CHKERRQ(petsc_err)                                                                                ! this is a macro definition, it is case sensitive

  call MPI_Comm_rank(PETSC_COMM_WORLD,worldrank,mpi_err)
  if (mpi_err /= 0) call quit(1)
  call MPI_Comm_size(PETSC_COMM_WORLD,worldsize,mpi_err)
  if (mpi_err /= 0) call quit(1)

  mainProcess: if (worldrank == 0) then
    if (output_unit /= 6) then
      write(output_unit,'(/,a)') ' ERROR: STDOUT != 6'
      call quit(1)
    endif
    if (error_unit /= 0) then
      write(output_unit,'(/,a)') ' ERROR: STDERR != 0'
      call quit(1)
    endif
  else mainProcess
    close(6)                                                                                        ! disable output for non-master processes (open 6 to rank specific file for debug)
    open(6,file='/dev/null',status='replace')                                                       ! close(6) alone will leave some temp files in cwd
  endif mainProcess

  write(6,'(/,a)') ' <<<+-  DAMASK_interface init  -+>>>'

  ! http://patorjk.com/software/taag/#p=display&f=Lean&t=DAMASK
  write(6,*) achar(27)//'[94m'
  write(6,*) '     _/_/_/      _/_/    _/      _/    _/_/      _/_/_/  _/    _/'
  write(6,*) '    _/    _/  _/    _/  _/_/  _/_/  _/    _/  _/        _/  _/'
  write(6,*) '   _/    _/  _/_/_/_/  _/  _/  _/  _/_/_/_/    _/_/    _/_/'
  write(6,*) '  _/    _/  _/    _/  _/      _/  _/    _/        _/  _/  _/'
  write(6,*) ' _/_/_/    _/    _/  _/      _/  _/    _/  _/_/_/    _/    _/'
  write(6,*) achar(27)//'[0m'

  write(6,'(/,a)') ' Roters et al., Computational Materials Science 158:420–478, 2019'
  write(6,'(a)')   ' https://doi.org/10.1016/j.commatsci.2018.04.030'

  write(6,'(/,a)') ' Version: '//DAMASKVERSION

 ! https://github.com/jeffhammond/HPCInfo/blob/master/docs/Preprocessor-Macros.md
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
  write(6,'(/,a)') ' Compiled with: '//compiler_version()
  write(6,'(a)')   ' Compiler options: '//compiler_options()
#elif defined(__INTEL_COMPILER)
  write(6,'(/,a,i4.4,a,i8.8)') ' Compiled with Intel fortran version :', __INTEL_COMPILER,&
                                                       ', build date :', __INTEL_COMPILER_BUILD_DATE
#elif defined(__PGI)
  write(6,'(a,i4.4,a,i8.8)')   ' Compiled with PGI fortran version :', __PGIC__,&
                                                                  '.', __PGIC_MINOR__
#endif

  write(6,'(/,a)') ' Compiled on: '//__DATE__//' at '//__TIME__

  call date_and_time(values = dateAndTime)
  write(6,'(/,a,2(i2.2,a),i4.4)') ' Date: ',dateAndTime(3),'/',dateAndTime(2),'/', dateAndTime(1)
  write(6,'(a,2(i2.2,a),i2.2)')   ' Time: ',dateAndTime(5),':', dateAndTime(6),':', dateAndTime(7)

  call MPI_Type_size(MPI_INTEGER,typeSize,mpi_err)
  if (mpi_err /= 0) call quit(1)
  if (typeSize*8 /= bit_size(0)) then
    write(6,'(a)') ' Mismatch between MPI and DAMASK integer' 
    call quit(1)
  endif

  call MPI_Type_size(MPI_DOUBLE,typeSize,mpi_err)
  if (mpi_err /= 0) call quit(1)
  if (typeSize*8 /= storage_size(0.0_pReal)) then
    write(6,'(a)') ' Mismatch between MPI and DAMASK real' 
    call quit(1)
  endif

  do i = 1, command_argument_count()
    call get_command_argument(i,arg)
    select case(trim(arg))                                                                          ! extract key
      case ('-h','--help')
        write(6,'(a)')  ' #######################################################################'
        write(6,'(a)')  ' DAMASK Command Line Interface:'
        write(6,'(a)')  ' For PETSc-based solvers for the Düsseldorf Advanced Material Simulation Kit'
        write(6,'(a,/)')' #######################################################################'
        write(6,'(a,/)')' Valid command line switches:'
        write(6,'(a)')  '    --geom         (-g, --geometry)'
        write(6,'(a)')  '    --load         (-l, --loadcase)'
        write(6,'(a)')  '    --workingdir   (-w, --wd, --workingdirectory, -d, --directory)'
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
        write(6,'(/,a)')'   --restart XX'
        write(6,'(a)')  '        Reads in increment XX and continues with calculating'
        write(6,'(a)')  '            increment XX+1 based on this.'
        write(6,'(a)')  '        Appends to existing results file'
        write(6,'(a)')  '            "NameOfGeom_NameOfLoadFile".'
        write(6,'(a)')  '        Works only if the restart information for increment XX'
        write(6,'(a)')  '            is available in the working directory.'
        write(6,'(/,a)')' -----------------------------------------------------------------------'
        write(6,'(a)')  ' Help:'
        write(6,'(/,a)')'   --help'
        write(6,'(a,/)')'        Prints this message and exits'
        call quit(0)                                                                                ! normal Termination
      case ('-l', '--load', '--loadcase')
        call get_command_argument(i+1,loadCaseArg)
      case ('-g', '--geom', '--geometry')
        call get_command_argument(i+1,geometryArg)
      case ('-w', '-d', '--wd', '--directory', '--workingdir', '--workingdirectory')
        call get_command_argument(i+1,workingDirArg)
      case ('-r', '--rs', '--restart')
        call get_command_argument(i+1,arg)
        read(arg,*,iostat=stat) interface_restartInc
        if (interface_restartInc < 0 .or. stat /=0) then
          write(6,'(/,a)') ' ERROR: Could not parse restart increment: '//trim(arg)
          call quit(1)
        endif
    end select
  enddo

  if (len_trim(loadcaseArg) == 0 .or. len_trim(geometryArg) == 0) then
    write(6,'(/,a)') ' ERROR: Please specify geometry AND load case (-h for help)'
    call quit(1)
  endif

  if (len_trim(workingDirArg) > 0) call setWorkingDirectory(trim(workingDirArg))
  geometryFile = getGeometryFile(geometryArg)
  loadCaseFile = getLoadCaseFile(loadCaseArg)

  call get_command(commandLine)
  call get_environment_variable('USER',userName)
  ! ToDo: https://stackoverflow.com/questions/8953424/how-to-get-the-username-in-c-c-in-linux
  write(6,'(/,a,i4.1)') ' MPI processes: ',worldsize
  write(6,'(a,a)')      ' Host name: ', trim(getHostName())
  write(6,'(a,a)')      ' User name: ', trim(userName)

  write(6,'(/a,a)')     ' Command line call:      ', trim(commandLine)
  if (len(trim(workingDirArg)) > 0) &
    write(6,'(a,a)')    ' Working dir argument:   ', trim(workingDirArg)
  write(6,'(a,a)')      ' Geometry argument:      ', trim(geometryArg)
  write(6,'(a,a)')      ' Load case argument:     ', trim(loadcaseArg)
  write(6,'(a,a)')      ' Working directory:      ', trim(getCWD())
  write(6,'(a,a)')      ' Geometry file:          ', trim(geometryFile)
  write(6,'(a,a)')      ' Loadcase file:          ', trim(loadCaseFile)
  write(6,'(a,a)')      ' Solver job name:        ', trim(getSolverJobName())
  if (interface_restartInc > 0) &
    write(6,'(a,i6.6)') ' Restart from increment: ', interface_restartInc

  !call signalterm_c(c_funloc(catchSIGTERM))
  call signalusr1_c(c_funloc(catchSIGUSR1))
  call signalusr2_c(c_funloc(catchSIGUSR2))
  call setSIGTERM(.false.)
  call setSIGUSR1(.false.)
  call setSIGUSR2(.false.)


end subroutine DAMASK_interface_init


!--------------------------------------------------------------------------------------------------
!> @brief extract working directory from given argument or from location of geometry file,
!!        possibly converting relative arguments to absolut path
!--------------------------------------------------------------------------------------------------
subroutine setWorkingDirectory(workingDirectoryArg)

  character(len=*),  intent(in) :: workingDirectoryArg                                              !< working directory argument
  character(len=1024)           :: workingDirectory                                                 !< working directory argument
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
    write(6,'(/,a)') ' ERROR: Working directory "'//trim(workingDirectory)//'" does not exist'
    call quit(1)
  endif

end subroutine setWorkingDirectory


!--------------------------------------------------------------------------------------------------
!> @brief solver job name (no extension) as combination of geometry and load case name
!--------------------------------------------------------------------------------------------------
character(len=1024) function getSolverJobName()

  integer :: posExt,posSep

  posExt = scan(geometryFile,'.',back=.true.)
  posSep = scan(geometryFile,'/',back=.true.)

  getSolverJobName = geometryFile(posSep+1:posExt-1)

  posExt = scan(loadCaseFile,'.',back=.true.)
  posSep = scan(loadCaseFile,'/',back=.true.)

  getSolverJobName = trim(getSolverJobName)//'_'//loadCaseFile(posSep+1:posExt-1)

end function getSolverJobName


!--------------------------------------------------------------------------------------------------
!> @brief basename of geometry file with extension from command line arguments
!--------------------------------------------------------------------------------------------------
character(len=1024) function getGeometryFile(geometryParameter)

  character(len=1024), intent(in) :: geometryParameter
  logical                         :: file_exists
  external                        :: quit

  getGeometryFile = trim(geometryParameter)
  if (scan(getGeometryFile,'/') /= 1) getGeometryFile = trim(getCWD())//'/'//trim(getGeometryFile)
  getGeometryFile = makeRelativePath(trim(getCWD()), getGeometryFile)

  inquire(file=trim(getGeometryFile), exist=file_exists)
  if (.not. file_exists) then
    write(6,'(/,a)') ' ERROR: Geometry file does not exists ('//trim(getGeometryFile)//')'
    call quit(1)
  endif

end function getGeometryFile


!--------------------------------------------------------------------------------------------------
!> @brief relative path of loadcase from command line arguments
!--------------------------------------------------------------------------------------------------
character(len=1024) function getLoadCaseFile(loadCaseParameter)

  character(len=1024), intent(in) :: loadCaseParameter
  logical                         :: file_exists
  external                        :: quit

  getLoadCaseFile = trim(loadCaseParameter)
  if (scan(getLoadCaseFile,'/') /= 1) getLoadCaseFile = trim(getCWD())//'/'//trim(getLoadCaseFile)
  getLoadCaseFile = makeRelativePath(trim(getCWD()), getLoadCaseFile)

  inquire(file=trim(getLoadCaseFile), exist=file_exists)
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

  character(len=*)    :: path
  character(len=1024) :: rectifyPath
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

end function rectifyPath


!--------------------------------------------------------------------------------------------------
!> @brief relative path from absolute a to absolute b
!--------------------------------------------------------------------------------------------------
character(len=1024) function makeRelativePath(a,b)

  character (len=*), intent(in) :: a,b
  character (len=1024)          :: a_cleaned,b_cleaned
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
!> @brief sets global variable SIGTERM to .true.
!--------------------------------------------------------------------------------------------------
subroutine catchSIGTERM(signal) bind(C)

  integer(C_INT), value :: signal
  SIGTERM = .true.

  write(6,'(a,i2.2,a)') ' received signal ',signal, ', set SIGTERM'

end subroutine catchSIGTERM


!--------------------------------------------------------------------------------------------------
!> @brief sets global variable SIGTERM
!--------------------------------------------------------------------------------------------------
subroutine setSIGTERM(state)

  logical, intent(in) :: state
  SIGTERM = state

end subroutine setSIGTERM


!--------------------------------------------------------------------------------------------------
!> @brief sets global variable SIGUSR1 to .true.
!--------------------------------------------------------------------------------------------------
subroutine catchSIGUSR1(signal) bind(C)

  integer(C_INT), value :: signal
  SIGUSR1 = .true.

  write(6,'(a,i2.2,a)') ' received signal ',signal, ', set SIGUSR1'

end subroutine catchSIGUSR1


!--------------------------------------------------------------------------------------------------
!> @brief sets global variable SIGUSR1
!--------------------------------------------------------------------------------------------------
subroutine setSIGUSR1(state)

  logical, intent(in) :: state
  SIGUSR1 = state

end subroutine setSIGUSR1


!--------------------------------------------------------------------------------------------------
!> @brief sets global variable SIGUSR2 to .true. if program receives SIGUSR2
!--------------------------------------------------------------------------------------------------
subroutine catchSIGUSR2(signal) bind(C)

  integer(C_INT), value :: signal
  SIGUSR2 = .true.

  write(6,'(a,i2.2,a)') ' received signal ',signal, ', set SIGUSR2'

end subroutine catchSIGUSR2


!--------------------------------------------------------------------------------------------------
!> @brief sets global variable SIGUSR2
!--------------------------------------------------------------------------------------------------
subroutine setSIGUSR2(state)

  logical, intent(in) :: state
  SIGUSR2 = state

end subroutine setSIGUSR2


end module
