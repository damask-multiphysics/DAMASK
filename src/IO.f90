! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  input/output functions
!--------------------------------------------------------------------------------------------------
module IO
  use, intrinsic :: ISO_C_binding
  use, intrinsic :: ISO_fortran_env, only: &
    IO_STDOUT => OUTPUT_UNIT, &
    IO_STDERR => ERROR_UNIT, &
    IO_STDIN => INPUT_UNIT

  use prec
  use constants
  use misc
#ifndef MARC_SOURCE
  use OS
#endif
implicit none(type,external)
  private


  interface
#ifndef MARC_SOURCE
    function isatty_stdout_C() bind(C)
      use, intrinsic :: ISO_C_binding, only: C_BOOL

      implicit none(type,external)
      logical(C_BOOL) :: isatty_stdout_C
    end function isatty_stdout_C

    function isatty_stderr_C() bind(C)
      use, intrinsic :: ISO_C_binding, only: C_BOOL

      implicit none(type,external)
      logical(C_BOOL) :: isatty_stderr_C
    end function isatty_stderr_C

    function isatty_stdin_C() bind(C)
      use, intrinsic :: ISO_C_binding, only: C_BOOL

      implicit none(type,external)
      logical(C_BOOL) :: isatty_stdin_C
    end function isatty_stdin_C
#else
    subroutine quit(stop_id)

      implicit none(type,external)
      integer, intent(in) :: stop_id
    end subroutine quit
#endif
  end interface

  ! For transition period
  interface IO_error
    module procedure IO_error_new
    module procedure IO_error_old
  end interface IO_error

  character, parameter, public :: &
    IO_ESC = achar(27), &                                                                           !< escape character
    IO_EOL = LF                                                                                     !< end of line character
  character(len=*), parameter, public :: &
    IO_FORMATRESET = IO_ESC//'[0m', &                                                               !< reset formatting
    IO_EMPH = IO_ESC//'[3m', &                                                                      !< emphasize (italics)
    IO_QUOTES  = "'"//'"' , &                                                                       !< quotes for strings
    IO_WHITESPACE = achar(44)//achar(32)//achar(9)//achar(10)//achar(13)                            !< whitespace characters

#ifndef MARC_SOURCE
    logical(C_BOOL), bind(C, name='IO_redirectedSTDOUT') :: IO_redirectedSTDOUT = .false.           !< STDOUT writes to file 'out.X' where X is the world rank
    logical(C_BOOL), bind(C, name='IO_redirectedSTDERR') :: IO_redirectedSTDERR = .false.           !< STDERR writes to file 'err.X' where X is the world rank
#endif
   logical :: IO_colored = .true.                                                                   !< status of colored output

  public :: &
    quit, &
    IO_init, &
    IO_selfTest, &
    IO_read, &
    IO_wrapLines, &
    IO_lc, &
    IO_glueDiffering, &
    IO_intAsStr, &
    IO_realAsStr, &
    IO_strAsInt, &
    IO_strAsReal, &
    IO_strAsBool, &
    IO_color, &
    IO_error, &
    IO_warning, &
    IO_STDOUT, &
    IO_STDERR, &
    tokenize

contains


!--------------------------------------------------------------------------------------------------
!> @brief Set options related to use of ANSI escape codes and do self test.
!--------------------------------------------------------------------------------------------------
subroutine IO_init()

  character(len=pSTRLEN) :: fname
  integer :: status


  print'(/,1x,a)', '<<<+-  IO init  -+>>>'; flush(IO_STDOUT)

#ifndef MARC_SOURCE
  ! redirection occurs in parallelization_init before any output is written
  inquire(unit=IO_STDOUT,name=fname)
  IO_redirectedSTDOUT = logical(fname(:4) == 'out.',C_BOOL)
  inquire(unit=IO_STDERR,name=fname)
  IO_redirectedSTDERR = logical(fname(:4) == 'err.',C_BOOL)
#endif
  call get_environment_variable('NO_COLOR',status=status)                                           !< https://no-color.org
  IO_colored = 0 /= status

  call IO_selfTest()

end subroutine IO_init


!--------------------------------------------------------------------------------------------------
!> @brief Read ASCII file.
!> @details Proper Unix style (LF line endings and LF at EOF) is ensured.
!--------------------------------------------------------------------------------------------------
function IO_read(fileName) result(fileContent)

  character(len=*),  intent(in) :: fileName
  character(len=:), allocatable :: fileContent

  integer ::  &
    fileUnit, &
    myStat
  integer(pI64) ::  &
    fileLength


  inquire(file = fileName, size=fileLength)
  open(newunit=fileUnit, file=fileName, access='stream',&
       status='old', position='rewind', action='read',iostat=myStat)
  if (myStat /= 0) call IO_error(100_pI16,'cannot open file',fileName,emph=[2])
  allocate(character(len=fileLength)::fileContent)
  if (fileLength==0) then
    close(fileUnit)
    return
  end if

  read(fileUnit,iostat=myStat) fileContent
  if (myStat /= 0) call IO_error(100_pI16,'cannot read from file',fileName,emph=[2])
  close(fileUnit)

  if (index(fileContent,CR//LF,kind=pI64) /= 0)     fileContent = CRLF2LF(fileContent)
  if (fileContent(fileLength:fileLength) /= IO_EOL) fileContent = fileContent//IO_EOL               ! ensure EOL@EOF

end function IO_read


!--------------------------------------------------------------------------------------------------
!> @brief Insert EOL at separator trying to keep line length below limit.
!--------------------------------------------------------------------------------------------------
function IO_wrapLines(str,separator,filler,length)

  character(len=*),    intent(in) :: str                                                            !< string to split
  character, optional, intent(in) :: separator                                                      !< line breaks are possible after this character, defaults to ','
  character(len=*), optional, intent(in) :: filler                                                  !< character(s) to insert after line break, defaults to none
  integer,   optional, intent(in) :: length                                                         !< (soft) line limit, defaults to 80
  character(len=:), allocatable :: IO_wrapLines

  integer, dimension(:), allocatable :: pos_sep, pos_split
  integer :: i,s,e


  i = index(str,misc_optional(separator,','))
  if (i == 0) then
    IO_wrapLines = str
  else
    pos_sep = [0]
    s = i
    do while (i /= 0 .and. s < len(str))
      pos_sep = [pos_sep,s]
      i = index(str(s+1:),misc_optional(separator,','))
      s = s + i
    end do
    pos_sep = [pos_sep,len(str)]

    pos_split = emptyIntArray
    s = 1
    e = 2
    IO_wrapLines = ''
    do while (e < size(pos_sep))
      if (pos_sep(e+1) - pos_sep(s) >= misc_optional(length,80)) then
        IO_wrapLines = IO_wrapLines//adjustl(str(pos_sep(s)+1:pos_sep(e)))//IO_EOL//misc_optional(filler,'')
        s = e
      end if
      e = e + 1
    end do
    IO_wrapLines = IO_wrapLines//adjustl(str(pos_sep(s)+1:))
  end if

end function IO_wrapLines


!--------------------------------------------------------------------------------------------------
!> @brief Convert characters in string to lower case.
!--------------------------------------------------------------------------------------------------
pure function IO_lc(str)

  character(len=*), intent(in) :: str                                                               !< string to convert
  character(len=len(str))   :: IO_lc

  integer :: i,n


  do i = 1,len(str)
    n = index(UPPER,str(i:i))
    if (n==0) then
      IO_lc(i:i) = str(i:i)
    else
      IO_lc(i:i) = LOWER(n:n)
    end if
  end do

end function IO_lc


!--------------------------------------------------------------------------------------------------
! @brief Return first (with glued on second if they differ).
!--------------------------------------------------------------------------------------------------
pure function IO_glueDiffering(first,second,glue)

  character(len=*),           intent(in)  :: first
  character(len=*),           intent(in)  :: second
  character(len=*), optional, intent(in)  :: glue
  character(len=:), allocatable :: IO_glueDiffering

  character(len=:), allocatable           :: glue_


  glue_ = misc_optional(glue,'<--')
  IO_glueDiffering = trim(first)
  if (trim(first) /= trim(second)) IO_glueDiffering = IO_glueDiffering//' '//trim(glue_)//' '//trim(second)

end function IO_glueDiffering


!--------------------------------------------------------------------------------------------------
!> @brief Return given int value as string.
!--------------------------------------------------------------------------------------------------
pure function IO_intAsStr(i)

  integer, intent(in)            :: i
  character(len=:), allocatable  :: IO_intAsStr


  allocate(character(len=merge(2,1,i<0) + floor(log10(real(abs(merge(1,i,i==0))))))::IO_intAsStr)
  write(IO_intAsStr,'(i0)') i

end function IO_intAsStr


!--------------------------------------------------------------------------------------------------
!> @brief Return given float value as string.
!--------------------------------------------------------------------------------------------------
pure function IO_realAsStr(f)

  real(pREAL), intent(in)        :: f
  character(len=:), allocatable  :: IO_realAsStr
  character(len=15)              :: tmp


  write(tmp,'(g15.7)') f
  tmp = adjustl(tmp)
  allocate(IO_realAsStr,source=tmp(:len_trim(tmp)))

end function IO_realAsStr


!--------------------------------------------------------------------------------------------------
!> @brief Return integer value from given string.
!--------------------------------------------------------------------------------------------------
integer function IO_strAsInt(str)

  character(len=*), intent(in) :: str                                                               !< string for conversion to int value

  integer :: readStatus


  read(str,*,iostat=readStatus) IO_strAsInt
  if (readStatus /= 0) call IO_error(111_pI16,'cannot represent',str,'as integer',emph=[2])

end function IO_strAsInt


!--------------------------------------------------------------------------------------------------
!> @brief Return real value from given string.
!--------------------------------------------------------------------------------------------------
real(pREAL) function IO_strAsReal(str)

  character(len=*), intent(in) :: str                                                               !< string for conversion to real value

  integer :: readStatus


  read(str,*,iostat=readStatus) IO_strAsReal
  if (readStatus /= 0) call IO_error(111_pI16,'cannot represent',str,'as real',emph=[2])

end function IO_strAsReal


!--------------------------------------------------------------------------------------------------
!> @brief Return logical value from given string.
!> @details: 'True' and 'true' are converted to .true.
!> @details: 'False' and 'false' are converted to .false.
!--------------------------------------------------------------------------------------------------
logical function IO_strAsBool(str)

  character(len=*), intent(in) :: str                                                               !< string for conversion to boolean


  if     (trim(adjustl(str)) == 'True' .or.  trim(adjustl(str)) == 'true') then
    IO_strAsBool = .true.
  elseif (trim(adjustl(str)) == 'False' .or. trim(adjustl(str)) == 'false') then
    IO_strAsBool = .false.
  else
    call IO_error(111_pI16,'cannot represent',str,'as boolean',emph=[2])
  end if

end function IO_strAsBool


!--------------------------------------------------------------------------------------------------
!> @brief Return string to set foreground and/or background color.
!> @details Only active if unit is a TTY. Does nothing for MSC.Marc. No color disables formatting.
!> @details https://stackoverflow.com/questions/4842424
!--------------------------------------------------------------------------------------------------
function IO_color(fg,bg,unit)

  character(len=:), allocatable :: IO_color
  integer, intent(in), dimension(3), optional :: &
    fg, &                                                                                           !< foreground color (8 bit RGB)
    bg                                                                                              !< background color (8 bit RGB)
  integer, intent(in), optional :: unit                                                             !< output unit (default STDOUT)


  IO_color = ''

  if (.not. IO_colored .or. .not. IO_isaTTY(misc_optional(unit,int(IO_STDOUT)))) return

  if (present(fg)) &
    IO_color = IO_color//IO_ESC//'[38;2;'//IO_intAsStr(fg(1))//';' &
                                         //IO_intAsStr(fg(2))//';' &
                                         //IO_intAsStr(fg(3))//'m'
  if (present(bg)) &
    IO_color = IO_color//IO_ESC//'[48;2;'//IO_intAsStr(bg(1))//';' &
                                         //IO_intAsStr(bg(2))//';' &
                                         //IO_intAsStr(bg(3))//'m'

  if (.not. present(fg) .and. .not. present(bg)) IO_color = IO_FORMATRESET

end function IO_color


!--------------------------------------------------------------------------------------------------
!> @brief Write error statements and terminate the run with exit #9xxx.
!> @details Should become "IO_error" after completed migration.
!--------------------------------------------------------------------------------------------------
subroutine IO_error_new(error_ID, &
                        info_1,info_2,info_3,info_4,info_5,info_6,info_7,info_8,info_9, &
                        emph)


  integer(pI16),      intent(in) :: error_ID        ! should go back to default integer after completed migration.
  class(*), optional, intent(in) :: info_1,info_2,info_3,info_4,info_5,info_6,info_7,info_8,info_9
  integer, dimension(:), optional, intent(in) :: emph                                               !< which info(s) to emphasize

  character(len=:), allocatable :: msg


  select case (error_ID)

!--------------------------------------------------------------------------------------------------
! file handling errors
    case (100)
      msg = 'file error'

!--------------------------------------------------------------------------------------------------
! file parsing errors
    case (110)
      msg = 'invalid chunk selected'
    case (111)
      msg = 'invalid string for conversion'
    case (114)
      msg = 'cannot decode base64 string:'

!--------------------------------------------------------------------------------------------------
! lattice error messages
    case (130)
      msg = 'invalid crystal parameters'
    case (138)
      msg = 'not enough interaction parameters given'

!--------------------------------------------------------------------------------------------------
! errors related to the parsing of material.yaml
    case (141)
      msg = 'number of chunks in string differs'
    case (142)
      msg = 'empty list'
    case (143)
      msg = 'key error'
    case (144)
      msg = 'negative number systems requested'
    case (145)
      msg = 'too many systems requested'
    case (146)
      msg = 'number of values does not match'
    case (147)
      msg = 'V_e needs to be symmetric'
    case (148)
      msg = 'Nconstituents mismatch between homogenization and material'

!--------------------------------------------------------------------------------------------------
! material error messages and related messages in geometry
    case (150)
      msg = 'index out of bounds'
    case (153)
      msg = 'sum of phase fractions differs from 1'
    case (155)
      msg = 'material index out of bounds'
    case (180)
      msg = 'missing/invalid material definition'
    case (190)
      msg = 'unknown element type:'
    case (191)
      msg = 'mesh contains more than one element type'

!--------------------------------------------------------------------------------------------------
! plasticity error messages
    case (200)
      msg = 'unknown type specified:'

    case (211)
      msg = 'material parameter out of bounds:'
    case (212)
      msg = 'nonlocal model not supported'

!--------------------------------------------------------------------------------------------------
! numerics error messages
    case (301)
      msg = 'numerics parameter out of bounds:'

!--------------------------------------------------------------------------------------------------
! math errors
    case (402)
      msg = 'invalid orientation specified'

!-------------------------------------------------------------------------------------------------
! homogenization errors
    case (500)
      msg = 'unknown homogenization specified'
    case (501)
      msg = 'homogenization description absent'

!--------------------------------------------------------------------------------------------------
! user errors
    case (600)
      msg = 'only one source entry allowed'
    case (601)
      msg = 'invalid option'
    case (603)
      msg = 'invalid data for table'
    case (610)
      msg = 'invalid command line options'
    case (640)
      msg = 'invalid working directory'


!------------------------------------------------------------------------------------------------
! errors related to YAML data
    case (701)
      msg = 'incorrect indent/Null value not allowed'
    case (702)
      msg = 'invalid use of flow YAML'
    case (703)
      msg = 'invalid YAML'
    case (704)
      msg = 'space expected after a colon for <key>: <value> pair'
    case (706)
      msg = 'type mismatch in YAML data node'
    case (708)
      msg = '"---" expected after YAML file header'
    case (709)
      msg = 'length mismatch'
    case (710)
      msg = 'closing quotation mark missing in string'

!-------------------------------------------------------------------------------------------------
! errors related to the mesh solver
    case (800)
      msg = 'invalid mesh'
    case (812)
      msg = 'invalid boundary conditions'

!-------------------------------------------------------------------------------------------------
! errors related to the grid solver
    case (830)
      msg = 'invalid load case'
    case (844)
      msg = 'invalid VTI file'
    case (894)
      msg = 'MPI error'

    case (950)
      msg = 'max number of cutbacks exceeded, terminating'

    case default
      error stop 'invalid error number'

  end select

  call panel('error',int(error_ID),msg, &
             info_1,info_2,info_3,info_4,info_5,info_6,info_7,info_8,info_9, &
             emph)
  call quit(9000+int(error_ID))

end subroutine IO_error_new


!--------------------------------------------------------------------------------------------------
!> @brief Write error statements and terminate the run with exit #9xxx.
!> @details Deprecated.
!--------------------------------------------------------------------------------------------------
subroutine IO_error_old(error_ID,ext_msg,label1,ID1,label2,ID2)

  integer,                    intent(in) :: error_ID
  character(len=*), optional, intent(in) :: ext_msg,label1,label2
  integer,          optional, intent(in) :: ID1,ID2

  character(len=:), allocatable :: msg_extra


  if (.not. present(label1) .and. present(ID1)) error stop 'missing label for value 1'
  if (.not. present(label2) .and. present(ID2)) error stop 'missing label for value 2'

  msg_extra = ''
  if (present(ext_msg)) msg_extra = msg_extra//ext_msg//IO_EOL
  if (present(label1)) then
    msg_extra = msg_extra//'at '//label1
    if (present(ID1)) msg_extra = msg_extra//' '//IO_intAsStr(ID1)
    msg_extra = msg_extra//IO_EOL
  end if
  if (present(label2)) then
    msg_extra = msg_extra//'at '//label2
    if (present(ID2)) msg_extra = msg_extra//' '//IO_intAsStr(ID2)
    msg_extra = msg_extra//IO_EOL
  end if

  call IO_error_new(int(error_ID,pI16),msg_extra,IO_EOL)

end subroutine IO_error_old


!--------------------------------------------------------------------------------------------------
!> @brief Write warning statements.
!--------------------------------------------------------------------------------------------------
subroutine IO_warning(warning_ID, &
                      info_1,info_2,info_3,info_4,info_5,info_6,info_7,info_8,info_9, &
                      emph)

  integer,                         intent(in) :: warning_ID
  class(*),              optional, intent(in) :: info_1,info_2,info_3,info_4,info_5,info_6,info_7,info_8,info_9
  integer, dimension(:), optional, intent(in) :: emph                                               !< which info(s) to emphasize

  character(len=:), allocatable :: msg


  select case (warning_ID)
    case (10)
      msg = 'deprecated keyword'
    case (207)
      msg = 'line truncated'
    case (600)
      msg = 'failed to converge'
    case (601)
      msg = 'unexpected stiffness'
    case (709)
      msg = 'read only the first document'

    case default
      error stop 'invalid warning number'
  end select

  call panel('warning',int(warning_ID),msg, &
             info_1,info_2,info_3,info_4,info_5,info_6,info_7,info_8,info_9, &
             emph)

end subroutine IO_warning


!--------------------------------------------------------------------------------------------------
!> @brief Test whether a file descriptor refers to a terminal.
!> @detail A terminal is neither a file nor a redirected STDOUT/STDERR/STDIN.
!>         This function cannot detect redirection when invoked via mpirun/mpiexec.
!--------------------------------------------------------------------------------------------------
logical function IO_isaTTY(unit)

  integer, intent(in) :: unit


  select case(unit)
#ifndef MARC_SOURCE
    case (IO_STDOUT)
      IO_isaTTY = .not. logical(IO_redirectedSTDOUT) .and. logical(isatty_stdout_C())
    case (IO_STDERR)
      IO_isaTTY = .not. logical(IO_redirectedSTDERR) .and. logical(isatty_stderr_C())
    case (IO_STDIN)
      IO_isaTTY = logical(isatty_stdin_C())
#endif
    case default
      IO_isaTTY = .false.
  end select

end function IO_isaTTY


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some IO functions.
!--------------------------------------------------------------------------------------------------
subroutine IO_selfTest()

  character(len=:),      allocatable :: str
  character(len=:), dimension(:), allocatable :: tokens


  if (dNeq(1.0_pREAL, IO_strAsReal('1.0')))          error stop 'IO_strAsReal'
  if (dNeq(1.0_pREAL, IO_strAsReal('1e0')))          error stop 'IO_strAsReal'
  if (dNeq(0.1_pREAL, IO_strAsReal('1e-1')))         error stop 'IO_strAsReal'
  if (dNeq(0.1_pREAL, IO_strAsReal('1.0e-1')))       error stop 'IO_strAsReal'
  if (dNeq(0.1_pREAL, IO_strAsReal('1.00e-1')))      error stop 'IO_strAsReal'
  if (dNeq(10._pREAL, IO_strAsReal(' 1.0e+1 ')))     error stop 'IO_strAsReal'

  if (3112019  /= IO_strAsInt( '3112019'))           error stop 'IO_strAsInt'
  if (3112019  /= IO_strAsInt(' 3112019'))           error stop 'IO_strAsInt'
  if (-3112019 /= IO_strAsInt('-3112019'))           error stop 'IO_strAsInt'
  if (3112019  /= IO_strAsInt('+3112019 '))          error stop 'IO_strAsInt'
  if (3112019  /= IO_strAsInt('03112019 '))          error stop 'IO_strAsInt'
  if (3112019  /= IO_strAsInt('+03112019'))          error stop 'IO_strAsInt'

  if (.not. IO_strAsBool(' true'))                   error stop 'IO_strAsBool'
  if (.not. IO_strAsBool(' True '))                  error stop 'IO_strAsBool'
  if (      IO_strAsBool(' false'))                  error stop 'IO_strAsBool'
  if (      IO_strAsBool('False'))                   error stop 'IO_strAsBool'

  if ('1234' /= IO_intAsStr(1234))                   error stop 'IO_intAsStr'
  if ('-12'  /= IO_intAsStr(-0012))                  error stop 'IO_intAsStr'

  if ('-0.1200000' /= IO_realAsStr(-0.12_pREAL))        error stop 'IO_realAsStr'
  if ('0.1234000E-31' /= IO_realAsStr(123.4e-34_pREAL)) error stop 'IO_realAsStr'

  if (CRLF2LF('') /= '')                             error stop 'CRLF2LF/0'
  if (CRLF2LF(LF)     /= LF)                         error stop 'CRLF2LF/1a'
  if (CRLF2LF(CR//LF) /= LF)                         error stop 'CRLF2LF/1b'
  if (CRLF2LF(' '//LF)     /= ' '//LF)               error stop 'CRLF2LF/2a'
  if (CRLF2LF(' '//CR//LF) /= ' '//LF)               error stop 'CRLF2LF/2b'
  if (CRLF2LF('A'//CR//LF//'B') /= 'A'//LF//'B')     error stop 'CRLF2LF/3'
  if (CRLF2LF('A'//CR//LF//'B'//CR//LF) /= &
              'A'//LF//'B'//LF)                      error stop 'CRLF2LF/4'
  if (CRLF2LF('A'//LF//CR//'B') /= 'A'//LF//CR//'B') error stop 'CRLF2LF/5'

  str='*(HiU!)3';if ('*(hiu!)3' /= IO_lc(str))       error stop 'IO_lc'

  if ('abc, def' /= IO_wrapLines('abc, def')) &
                                                     error stop 'IO_wrapLines/1'
  if ('abc,'//IO_EOL//'def' /= IO_wrapLines('abc,def',length=3)) &
                                                     error stop 'IO_wrapLines/2'
  if ('abc,'//IO_EOL//'def' /= IO_wrapLines('abc,def',length=5)) &
                                                     error stop 'IO_wrapLines/3'
  if ('abc, def' /= IO_wrapLines('abc, def',length=3,separator='.')) &
                                                     error stop 'IO_wrapLines/4'
  if ('abc.'//IO_EOL//'def' /= IO_wrapLines('abc. def',length=3,separator='.')) &
                                                     error stop 'IO_wrapLines/5'
  if ('abc,'//IO_EOL//'defg,'//IO_EOL//'hij' /= IO_wrapLines('abc,defg,hij',length=4)) &
                                                     error stop 'IO_wrapLines/6'
  if ('abc,'//IO_EOL//'xxdefg,'//IO_EOL//'xxhij' /= IO_wrapLines('abc,defg, hij',filler='xx',length=4)) &
                                                     error stop 'IO_wrapLines/7'

#if ((defined(__INTEL_COMPILER) && __INTEL_COMPILER_BUILD_DATE < 20240000) || !defined(__INTEL_COMPILER))
  call tokenize('','$',tokens)
  if (size(tokens) /= 0 .or. len(tokens) /=0) error stop 'tokenize empty'
  call tokenize('abcd','dcba',tokens)
  if (size(tokens) /= 0 .or. len(tokens) /=0) error stop 'tokenize only separators'

  tokens=['a']
  call test_tokenize('a','#',tokens)
  call test_tokenize('#a','#',tokens)
  call test_tokenize('a#','#',tokens)

  tokens=['aa']
  call test_tokenize('aa','#',tokens)
  call test_tokenize('$aa','$',tokens)
  call test_tokenize('aa$','$',tokens)

  tokens=['a','b']
  call test_tokenize('a$b','$',tokens)
  call test_tokenize('@a@$b@','$@',tokens)

  tokens=['aa','bb']
  call test_tokenize('aa$bb','$',tokens)
  call test_tokenize('aa$$bb','$',tokens)
  call test_tokenize('aa$bb$','$',tokens)

  tokens=['aa  ','bbb ','cccc']
  call test_tokenize('aa$bbb$cccc','$',tokens)
  call test_tokenize('$aa$bbb$cccc$','$',tokens)
  call tokenize('#aa@@bbb!!!cccc#','#@!',tokens)


  contains
  pure subroutine test_tokenize(input,delimiter,solution)
    character(len=*), intent(in) :: input, delimiter
    character(len=*), dimension(:), intent(in) :: solution

    character(len=:), dimension(:), allocatable :: tok
    integer :: i


    call tokenize(input,delimiter,tok)
    do i = 1,size(tok)
      !if (solution(i) /= tok(i)) error stop 'tokenize "'//solution(i)//'" vs. "'//tok(i)//'"'      ! requires 2018 standard
      if (solution(i) /= tok(i)) error stop 'tokenize'
    end do

  end subroutine test_tokenize
#endif

end subroutine IO_selfTest


#ifndef MARC_SOURCE
!--------------------------------------------------------------------------------------------------
!> @brief Stop execution and report status.
!> @details exits the program and reports current time and duration. Exit code 0 signals
!> everything is fine. Exit code 1 signals an error, message according to IO_error.
!--------------------------------------------------------------------------------------------------
subroutine quit(stop_id)
  use, intrinsic :: ISO_fortran_env, only: ERROR_UNIT, OUTPUT_UNIT
#include <petsc/finclude/petscsys.h>
  use PETScSys
#ifndef PETSC_HAVE_MPI_F90MODULE_VISIBILITY
  use MPI_f08
#endif
  use HDF5

#ifndef PETSC_HAVE_MPI_F90MODULE_VISIBILITY
  implicit none(type,external)
#else
  implicit none
#endif

  integer, intent(in) :: stop_id

  integer, dimension(8) :: date_time
  integer :: err_HDF5
  integer(MPI_INTEGER_KIND) :: err_MPI, worldsize
  PetscErrorCode :: err_PETSc


  call H5Open_f(err_HDF5)                                                                           ! prevents error if not opened yet
  if (err_HDF5 < 0) write(ERROR_UNIT,'(a,i0)') ' Error in H5Open_f ',err_HDF5
  call H5Close_f(err_HDF5)
  if (err_HDF5 < 0) write(ERROR_UNIT,'(a,i0)') ' Error in H5Close_f ',err_HDF5

  call PetscFinalize(err_PETSc)

  call date_and_time(values = date_time)
  write(OUTPUT_UNIT,'(/,a)') ' DAMASK terminated on:'
  print'(3x,a,1x,2(i2.2,a),i4.4)', 'Date:',date_time(3),'/',date_time(2),'/',date_time(1)
  print'(3x,a,1x,2(i2.2,a),i2.2)', 'Time:',date_time(5),':',date_time(6),':',date_time(7)

  if (stop_id == 0 .and. err_HDF5 == 0 .and. err_PETSC == 0) then
    call MPI_Finalize(err_MPI)
    if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI_Finalize error'
    stop 0                                                                                          ! normal termination
  else
    call MPI_Comm_size(MPI_COMM_WORLD,worldsize,err_MPI)
    if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI_Comm error'
    if (stop_id /= 0 .and. worldsize > 1) call MPI_Abort(MPI_COMM_WORLD,1,err_MPI)
    stop 1                                                                                          ! error (message from IO_error)
  endif

end subroutine quit

!--------------------------------------------------------------------------------------------------
!> @brief Print C string to Fortran stdout.
!--------------------------------------------------------------------------------------------------
subroutine IO_printCppString(C_STR) bind(C, name='F_IO_printCppString')

  character(kind=C_CHAR), intent(in), dimension(*) :: c_str


  write (IO_STDOUT, '(a)', advance='no') c_f_string(c_str)
  flush(IO_STDOUT)

end subroutine IO_printCppString
#endif


!--------------------------------------------------------------------------------------------------
!> @brief Convert Windows (CRLF) to Unix (LF) line endings.
!--------------------------------------------------------------------------------------------------
pure function CRLF2LF(str)

  character(len=*), intent(in)  :: str
  character(len=:), allocatable :: CRLF2LF

  integer(pI64) :: c,n


  allocate(character(len=len_trim(str,pI64))::CRLF2LF)
  if (len(CRLF2LF,pI64) == 0) return

  n = 0_pI64
  do c=1_pI64, len_trim(str,pI64)
    CRLF2LF(c-n:c-n) = str(c:c)
    if (c == len_trim(str,pI64)) exit
    if (str(c:c+1_pI64) == CR//LF) n = n + 1_pI64
  end do

  CRLF2LF = CRLF2LF(:c-n)

end function CRLF2LF

#if ((defined(__INTEL_COMPILER) && __INTEL_COMPILER_BUILD_DATE < 20240000) || !defined(__INTEL_COMPILER))
!--------------------------------------------------------------------------------------------------
!> @brief Fortran 2023 "tokenize" (first form, without optional argument).
!--------------------------------------------------------------------------------------------------
pure subroutine tokenize(string,set,tokens)

  character(len=*), intent(in) :: string, set
  character(len=:), dimension(:), allocatable, intent(out) :: tokens

  integer, allocatable, dimension(:,:) :: pos
  integer :: i, s, e


  allocate(pos(2,0))
  e = 0
  do while (e < verify(string,set,back=.true.))
    s = e + merge(verify(string(e+1:),set),1,scan(string(e+1:),set)/=0)
    e = s + merge(scan(string(s:),set)-2,len(string(s:))-1,scan(string(s:),set)/=0)
    pos = reshape([pos,[s,e]],[2,size(pos)/2+1])
  end do
  allocate(character(len=merge(maxval(pos(2,:)-pos(1,:))+1,0,size(pos)>0))::tokens(size(pos,2)))
  do i = 1, size(pos,2)
    tokens(i) = string(pos(1,i):pos(2,i))
  end do

end subroutine tokenize
#endif

!--------------------------------------------------------------------------------------------------
!> @brief Write statements to standard error.
!--------------------------------------------------------------------------------------------------
subroutine panel(paneltype,ID,msg, &
                 info_1,info_2,info_3,info_4,info_5,info_6,info_7,info_8,info_9, &
                 emph)

  character(len=*),           intent(in) :: paneltype, &                                            !< either 'error' or 'warning'
                                            msg                                                     !< general error/warning message
  integer,                    intent(in) :: ID                                                      !< error/warning ID
  class(*),         optional, intent(in) :: info_1,info_2,info_3,info_4,info_5,info_6,info_7,info_8,info_9 !< extra info
  integer, dimension(:), optional, intent(in) :: emph                                               !< which info(s) to emphasize

  integer, parameter :: panelwidth = 69
  character(len=*), parameter :: DIVIDER = repeat('─',panelwidth)
  character(len=pSTRLEN) :: formatString
  character(len=:), allocatable :: heading, msg_, info_extra
  character(len=:), dimension(:), allocatable :: info_split
  integer :: len_corrected, &                                                                       !< string length corrected for control characters
             i


  ! Needed to avoid output glitches observed with Gfortran.
  ! see https://fortran-lang.discourse.group/t/openmp-and-thread-safety-of-i-os-write-read/4567/19
  !$OMP CRITICAL (internal_IO)
  heading = paneltype//' '//IO_intAsStr(ID)

  select case (paneltype)

    case ('error')
       msg_ = IO_color([255,0,0],  unit=IO_STDERR)//trim(msg)//IO_color(unit=IO_STDERR)
    case ('warning')
       msg_ = IO_color([255,192,0],unit=IO_STDERR)//trim(msg)//IO_color(unit=IO_STDERR)
    case default
       error stop 'invalid panel type: '//trim(paneltype)

  end select

  info_extra = as_str(info_1,is_emph(1,emph)) &
            // as_str(info_2,is_emph(2,emph)) &
            // as_str(info_3,is_emph(3,emph)) &
            // as_str(info_4,is_emph(4,emph)) &
            // as_str(info_5,is_emph(5,emph)) &
            // as_str(info_6,is_emph(6,emph)) &
            // as_str(info_7,is_emph(7,emph)) &
            // as_str(info_8,is_emph(8,emph)) &
            // as_str(info_9,is_emph(9,emph))
  !$OMP END CRITICAL (internal_IO)

  !$OMP CRITICAL (output_to_screen)
  write(IO_STDERR,'(/,a)')                ' ┌'       //DIVIDER//        '┐'
  write(formatString,'(a,i2,a)') '(a,24x,a,',max(1,panelwidth-24-len_trim(heading)),'x,a)'
  write(IO_STDERR,formatString)           ' │',    trim(heading),       '│'
  write(IO_STDERR,'(a)')                  ' ├'       //DIVIDER//        '┤'
  write(formatString,'(a,i3.3,a,i3.3,a)') '(a,a',max(1,len_trim(msg_)),',',&
                                                 max(1,panelwidth+3-len_trim(msg)-4),'x,a)'
  write(IO_STDERR,formatString)           ' │ ',     trim(msg_),        '│'
  if (len_trim(info_extra) > 0) then
    call tokenize(info_extra,IO_EOL,info_split)
    do i = 1, size(info_split)
      info_extra = adjustl(info_split(i))
      if (len_trim(info_extra) == 0) then
        write(IO_STDERR,'(a)')            ' │'//repeat(' ',panelwidth)//'│'
      else
        len_corrected = len_trim(info_extra) - count([(info_extra(i:i)==IO_ESC,i=1,len_trim(info_extra))])*4
        write(formatString,'(a,i3.3,a,i3.3,a)') '(a,a',max(1,len_trim(info_extra)),',',&
                                                       max(1,panelwidth+3-len_corrected-4),'x,a)'
        write(IO_STDERR,formatString)     ' │ ',    trim(info_extra),   '│'
      end if
    end do
  endif
  write(IO_STDERR,'(a)')                  ' └'       //DIVIDER//        '┘'
  flush(IO_STDERR)
  !$OMP END CRITICAL (output_to_screen)

end subroutine panel


!-----------------------------------------------------------------------------------------------
!> @brief Convert to string with white space prefix and optional emphasis.
!-----------------------------------------------------------------------------------------------
function as_str(info,emph)

  character(len=:), allocatable :: as_str
  class(*), optional, intent(in) :: info                                                           !< info message
  logical, intent(in) :: emph                                                                      !< whether info should be emphasized


  if (present(info)) then
    select type(info)
      type is (character(*))
        as_str = info
      type is (integer)
        as_str = IO_intAsStr(info)
      type is (real(pREAL))
        as_str = IO_realAsStr(info)
      class default
        error stop 'cannot convert info argument to string'
    end select

    if (emph) then
      if (IO_colored .and. IO_isaTTY(IO_STDERR)) then
        as_str = IO_EMPH//as_str//IO_FORMATRESET
      else
        as_str = IO_QUOTES(2:2)//as_str//IO_QUOTES(2:2)
      end if
    end if
    as_str = ' '//as_str
  else
    as_str = ''
  end if

end function as_str

!-----------------------------------------------------------------------------------------------
!> @brief Determine whether info at given position has to be emphasized.
!-----------------------------------------------------------------------------------------------
pure logical function is_emph(idx,emph)

  integer, intent(in) :: idx                                                                        !< index of considered info
  integer, dimension(:), optional, intent(in) :: emph                                               !< which info(s) to emphasize


  if (present(emph)) then
    is_emph = any(emph == idx)
  else
    is_emph = .false.
  end if

end function is_emph


end module IO
