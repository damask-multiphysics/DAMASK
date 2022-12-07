!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  input/output functions
!--------------------------------------------------------------------------------------------------
module IO
  use, intrinsic :: ISO_fortran_env, only: &
    IO_STDOUT => OUTPUT_UNIT, &
    IO_STDERR => ERROR_UNIT

  use prec

  implicit none(type,external)
  private

  character(len=*), parameter, public :: &
    IO_WHITESPACE = achar(44)//achar(32)//achar(9)//achar(10)//achar(13), &                         !< whitespace characters
    IO_QUOTES  = "'"//'"'
  character, parameter, public :: &
    IO_EOL = new_line('DAMASK'), &                                                                  !< end of line character
    IO_COMMENT = '#'
  character, parameter :: &
    CR = achar(13), &
    LF = IO_EOL

  public :: &
    IO_init, &
    IO_read, &
    IO_readlines, &
    IO_isBlank, &
    IO_stringPos, &
    IO_stringValue, &
    IO_intValue, &
    IO_floatValue, &
    IO_lc, &
    IO_rmComment, &
    IO_intAsString, &
    IO_stringAsInt, &
    IO_stringAsFloat, &
    IO_stringAsBool, &
    IO_error, &
    IO_warning, &
    IO_STDOUT

contains


!--------------------------------------------------------------------------------------------------
!> @brief Do self test.
!--------------------------------------------------------------------------------------------------
subroutine IO_init()

  print'(/,1x,a)', '<<<+-  IO init  -+>>>'; flush(IO_STDOUT)

  call selfTest()

end subroutine IO_init


!--------------------------------------------------------------------------------------------------
!> @brief Read ASCII file and split at EOL.
!--------------------------------------------------------------------------------------------------
function IO_readlines(fileName) result(fileContent)

  character(len=*),          intent(in)                :: fileName
  character(len=pStringLen), dimension(:), allocatable :: fileContent                               !< file content, separated per lines

  character(len=pStringLen)                            :: line
  character(len=:),                        allocatable :: rawData
  integer ::  &
    startPos, endPos, &
    N_lines, &                                                                                      !< # lines in file
    l
  logical :: warned


  rawData = IO_read(fileName)

  N_lines = count([(rawData(l:l) == IO_EOL,l=1,len(rawData))])
  allocate(fileContent(N_lines))

!--------------------------------------------------------------------------------------------------
! split raw data at end of line
  warned = .false.
  startPos = 1
  l = 1
  do while (l <= N_lines)
    endPos = startPos + scan(rawData(startPos:),IO_EOL) - 2
    if (endPos - startPos > pStringLen-1) then
      line = rawData(startPos:startPos+pStringLen-1)
      if (.not. warned) then
        call IO_warning(207,trim(fileName),label1='line',ID1=l)
        warned = .true.
      end if
    else
      line = rawData(startPos:endpos)
    end if
    startPos = endPos + 2                                                                           ! jump to next line start

    fileContent(l) = trim(line)//''
    l = l + 1
  end do

end function IO_readlines


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
  if (myStat /= 0) call IO_error(100,trim(fileName))
  allocate(character(len=fileLength)::fileContent)
  if (fileLength==0) then
    close(fileUnit)
    return
  end if

  read(fileUnit,iostat=myStat) fileContent
  if (myStat /= 0) call IO_error(102,trim(fileName))
  close(fileUnit)

  if (scan(fileContent(:index(fileContent,LF)),CR//LF) /= 0) fileContent = CRLF2LF(fileContent)
  if (fileContent(fileLength:fileLength) /= IO_EOL)          fileContent = fileContent//IO_EOL      ! ensure EOL@EOF

end function IO_read


!--------------------------------------------------------------------------------------------------
!> @brief Identifiy strings without content.
!--------------------------------------------------------------------------------------------------
logical pure function IO_isBlank(string)

  character(len=*), intent(in) :: string                                                            !< string to check for content

  integer :: posNonBlank


  posNonBlank = verify(string,IO_WHITESPACE)
  IO_isBlank = posNonBlank == 0 .or. posNonBlank == scan(string,IO_COMMENT)

end function IO_isBlank


!--------------------------------------------------------------------------------------------------
!> @brief Locate all whitespace-separated chunks in given string and returns array containing
!! number them and the left/right position to be used by IO_xxxVal.
!! Array size is dynamically adjusted to number of chunks found in string
!! IMPORTANT: first element contains number of chunks!
!--------------------------------------------------------------------------------------------------
pure function IO_stringPos(string)

  character(len=*),                  intent(in) :: string                                           !< string in which chunk positions are searched for
  integer, dimension(:), allocatable            :: IO_stringPos

  integer                      :: left, right


  allocate(IO_stringPos(1), source=0)
  right = 0

  do while (verify(string(right+1:),IO_WHITESPACE)>0)
    left  = right + verify(string(right+1:),IO_WHITESPACE)
    right = left + scan(string(left:),IO_WHITESPACE) - 2
    if ( string(left:left) == IO_COMMENT) exit
    IO_stringPos = [IO_stringPos,left,right]
    IO_stringPos(1) = IO_stringPos(1)+1
    endOfString: if (right < left) then
      IO_stringPos(IO_stringPos(1)*2+1) = len_trim(string)
      exit
    end if endOfString
  end do

end function IO_stringPos


!--------------------------------------------------------------------------------------------------
!> @brief Read string value at myChunk from string.
!--------------------------------------------------------------------------------------------------
function IO_stringValue(string,chunkPos,myChunk)

  character(len=*),             intent(in) :: string                                                !< raw input with known start and end of each chunk
  integer,   dimension(:),      intent(in) :: chunkPos                                              !< positions of start and end of each tag/chunk in given string
  integer,                      intent(in) :: myChunk                                               !< position number of desired chunk
  character(len=:), allocatable            :: IO_stringValue

  validChunk: if (myChunk > chunkPos(1) .or. myChunk < 1) then
    IO_stringValue = ''
    call IO_error(110,'IO_stringValue: "'//trim(string)//'"',label1='chunk',ID1=myChunk)
  else validChunk
    IO_stringValue = string(chunkPos(myChunk*2):chunkPos(myChunk*2+1))
  end if validChunk

end function IO_stringValue


!--------------------------------------------------------------------------------------------------
!> @brief Read integer value at myChunk from string.
!--------------------------------------------------------------------------------------------------
integer function IO_intValue(string,chunkPos,myChunk)

  character(len=*),      intent(in) :: string                                                       !< raw input with known start and end of each chunk
  integer, dimension(:), intent(in) :: chunkPos                                                     !< positions of start and end of each tag/chunk in given string
  integer,               intent(in) :: myChunk                                                      !< position number of desired chunk

  IO_intValue = IO_stringAsInt(IO_stringValue(string,chunkPos,myChunk))

end function IO_intValue


!--------------------------------------------------------------------------------------------------
!> @brief Read float value at myChunk from string.
!--------------------------------------------------------------------------------------------------
real(pReal) function IO_floatValue(string,chunkPos,myChunk)

  character(len=*),        intent(in) :: string                                                     !< raw input with known start and end of each chunk
  integer,   dimension(:), intent(in) :: chunkPos                                                   !< positions of start and end of each tag/chunk in given string
  integer,                 intent(in) :: myChunk                                                    !< position number of desired chunk

  IO_floatValue = IO_stringAsFloat(IO_stringValue(string,chunkPos,myChunk))

end function IO_floatValue


!--------------------------------------------------------------------------------------------------
!> @brief Convert characters in string to lower case.
!--------------------------------------------------------------------------------------------------
pure function IO_lc(string)

  character(len=*), intent(in) :: string                                                            !< string to convert
  character(len=len(string))   :: IO_lc

  character(len=*),          parameter :: LOWER = 'abcdefghijklmnopqrstuvwxyz'
  character(len=len(LOWER)), parameter :: UPPER = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  integer :: i,n


  do i = 1,len(string)
    n = index(UPPER,string(i:i))
    if (n==0) then
      IO_lc(i:i) = string(i:i)
    else
      IO_lc(i:i) = LOWER(n:n)
    end if
  end do

end function IO_lc


!--------------------------------------------------------------------------------------------------
! @brief Remove comments (characters beyond '#') and trailing space.
! ToDo: Discuss name (the trim aspect is not clear)
!--------------------------------------------------------------------------------------------------
function IO_rmComment(line)

  character(len=*), intent(in)  :: line
  character(len=:), allocatable :: IO_rmComment
  integer :: split


  split = index(line,IO_COMMENT)

  if (split == 0) then
    IO_rmComment = trim(line)
  else
    IO_rmComment = trim(line(:split-1))
  end if

end function IO_rmComment


!--------------------------------------------------------------------------------------------------
!> @brief Return given int value as string.
!--------------------------------------------------------------------------------------------------
function IO_intAsString(i)

  integer, intent(in)            :: i

  character(len=:), allocatable  :: IO_intAsString

  allocate(character(len=merge(2,1,i<0) + floor(log10(real(abs(merge(1,i,i==0))))))::IO_intAsString)
  write(IO_intAsString,'(i0)') i

end function IO_intAsString


!--------------------------------------------------------------------------------------------------
!> @brief Return integer value from given string.
!--------------------------------------------------------------------------------------------------
integer function IO_stringAsInt(string)

  character(len=*), intent(in) :: string                                                            !< string for conversion to int value

  integer                      :: readStatus
  character(len=*), parameter  :: VALIDCHARS = '0123456789+- '


  valid: if (verify(string,VALIDCHARS) == 0) then
    read(string,*,iostat=readStatus) IO_stringAsInt
    if (readStatus /= 0) call IO_error(111,string)
  else valid
    IO_stringAsInt = 0
    call IO_error(111,string)
  end if valid

end function IO_stringAsInt


!--------------------------------------------------------------------------------------------------
!> @brief Return float value from given string.
!--------------------------------------------------------------------------------------------------
real(pReal) function IO_stringAsFloat(string)

  character(len=*), intent(in) :: string                                                            !< string for conversion to float value

  integer                      :: readStatus
  character(len=*), parameter  :: VALIDCHARS = '0123456789eE.+- '


  valid: if (verify(string,VALIDCHARS) == 0) then
    read(string,*,iostat=readStatus) IO_stringAsFloat
    if (readStatus /= 0) call IO_error(112,string)
  else valid
    IO_stringAsFloat = 0.0_pReal
    call IO_error(112,string)
  end if valid

end function IO_stringAsFloat


!--------------------------------------------------------------------------------------------------
!> @brief Return logical value from given string.
!--------------------------------------------------------------------------------------------------
logical function IO_stringAsBool(string)

  character(len=*), intent(in) :: string                                                            !< string for conversion to int value


  if     (trim(adjustl(string)) == 'True' .or.  trim(adjustl(string)) == 'true') then
    IO_stringAsBool = .true.
  elseif (trim(adjustl(string)) == 'False' .or. trim(adjustl(string)) == 'false') then
    IO_stringAsBool = .false.
  else
    IO_stringAsBool = .false.
    call IO_error(113,string)
  end if

end function IO_stringAsBool



!--------------------------------------------------------------------------------------------------
!> @brief Write error statements and terminate the run with exit #9xxx.
!--------------------------------------------------------------------------------------------------
subroutine IO_error(error_ID,ext_msg,label1,ID1,label2,ID2)

  integer,                    intent(in) :: error_ID
  character(len=*), optional, intent(in) :: ext_msg,label1,label2
  integer,          optional, intent(in) :: ID1,ID2

  external                      :: quit
  character(len=:), allocatable :: msg

  select case (error_ID)

!--------------------------------------------------------------------------------------------------
! file handling errors
    case (100)
      msg = 'could not open file:'
    case (102)
      msg = 'could not read file:'

!--------------------------------------------------------------------------------------------------
! file parsing errors
    case (110)
      msg = 'invalid chunk selected'
    case (111)
      msg = 'invalid character for int:'
    case (112)
      msg = 'invalid character for float:'
    case (113)
      msg = 'invalid character for logical:'
    case (114)
      msg = 'cannot decode base64 string:'

!--------------------------------------------------------------------------------------------------
! lattice error messages
    case (130)
      msg = 'unknown lattice structure encountered'
    case (131)
      msg = 'hex lattice structure with invalid c/a ratio'
    case (132)
      msg = 'trans_lattice_structure not possible'
    case (134)
      msg = 'negative lattice parameter'
    case (135)
      msg = 'zero entry on stiffness diagonal'
    case (137)
      msg = 'not defined for lattice structure'
    case (138)
      msg = 'not enough interaction parameters given'

!--------------------------------------------------------------------------------------------------
! errors related to the parsing of material.yaml
    case (140)
      msg = 'key not found'
    case (141)
      msg = 'number of chunks in string differs'
    case (142)
      msg = 'empty list'
    case (143)
      msg = 'no value found for key'
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
      msg = 'unknown elasticity specified:'
    case (201)
      msg = 'unknown plasticity specified:'

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

!--------------------------------------------------------------------------------------------------
! user errors
    case (602)
      msg = 'invalid selection for debug'
    case (603)
      msg = 'invalid data for table'

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
    case (705)
      msg = 'unsupported feature'
    case (706)
      msg = 'type mismatch in YAML data node'
    case (707)
      msg = 'abrupt end of file'
    case (708)
      msg = '"---" expected after YAML file header'
    case (709)
      msg = 'length mismatch'
    case (710)
      msg = 'closing quotation mark missing in string'
    case (711)
      msg = 'incorrect type'

!-------------------------------------------------------------------------------------------------
! errors related to the mesh solver
    case (821)
      msg = 'order not supported'

!-------------------------------------------------------------------------------------------------
! errors related to the grid solver
    case (831)
      msg = 'mask consistency violated in grid load case'
    case (833)
      msg = 'non-positive ratio for geometric progression'
    case (834)
      msg = 'negative time increment in grid load case'
    case (835)
      msg = 'non-positive increments in grid load case'
    case (836)
      msg = 'non-positive result frequency in grid load case'
    case (837)
      msg = 'incomplete loadcase'
    case (838)
      msg = 'mixed boundary conditions allow rotation'
    case (839)
      msg = 'non-positive restart frequency in grid load case'
    case (844)
      msg = 'invalid VTI file'
    case (891)
      msg = 'unknown solver type selected'
    case (892)
      msg = 'unknown filter type selected'
    case (894)
      msg = 'MPI error'

    case (950)
      msg = 'max number of cutbacks exceeded, terminating'

    case default
      error stop 'invalid error number'

  end select

  call panel('error',error_ID,msg, &
                 ext_msg=ext_msg, &
                 label1=label1,ID1=ID1, &
                 label2=label2,ID2=ID2)
  call quit(9000+error_ID)

end subroutine IO_error


!--------------------------------------------------------------------------------------------------
!> @brief Write warning statements.
!--------------------------------------------------------------------------------------------------
subroutine IO_warning(warning_ID,ext_msg,label1,ID1,label2,ID2)

  integer,                    intent(in) :: warning_ID
  character(len=*), optional, intent(in) :: ext_msg,label1,label2
  integer,          optional, intent(in) :: ID1,ID2

  character(len=:), allocatable :: msg


  select case (warning_ID)
    case (47)
      msg = 'invalid parameter for FFTW'
    case (207)
      msg = 'line truncated'
    case (600)
      msg = 'crystallite responds elastically'
    case (601)
      msg = 'stiffness close to zero'
    case (709)
      msg = 'read only the first document'

    case default
      error stop 'invalid warning number'
  end select

  call panel('warning',warning_ID,msg, &
                 ext_msg=ext_msg, &
                 label1=label1,ID1=ID1, &
                 label2=label2,ID2=ID2)

end subroutine IO_warning


!--------------------------------------------------------------------------------------------------
!> @brief Convert Windows (CRLF) to Unix (LF) line endings.
!--------------------------------------------------------------------------------------------------
pure function CRLF2LF(string)

  character(len=*), intent(in)  :: string
  character(len=:), allocatable :: CRLF2LF

  integer :: c,n


  allocate(character(len=len_trim(string))::CRLF2LF)
  if (len(CRLF2LF) == 0) return

  n = 0
  do c=1, len_trim(string)
    CRLF2LF(c-n:c-n) = string(c:c)
    if (c == len_trim(string)) exit
    if (string(c:c+1) == CR//LF) n = n + 1
  end do

  CRLF2LF = CRLF2LF(:c-n)

end function CRLF2LF


!--------------------------------------------------------------------------------------------------
!> @brief Write statements to standard error.
!--------------------------------------------------------------------------------------------------
subroutine panel(paneltype,ID,msg,ext_msg,label1,ID1,label2,ID2)

  character(len=*),           intent(in) :: paneltype,msg
  character(len=*), optional, intent(in) :: ext_msg,label1,label2
  integer,                    intent(in) :: ID
  integer,          optional, intent(in) :: ID1,ID2

  character(len=pStringLen)              :: formatString
  integer, parameter                     :: panelwidth = 69
  character(len=*), parameter            :: DIVIDER = repeat('─',panelwidth)


  if (.not. present(label1) .and.       present(ID1)) error stop 'missing label for value 1'
  if (.not. present(label2) .and.       present(ID2)) error stop 'missing label for value 2'
  if (      present(label1) .and. .not. present(ID1)) error stop 'missing value for label 1'
  if (      present(label2) .and. .not. present(ID2)) error stop 'missing value for label 2'

  !$OMP CRITICAL (write2out)
  write(IO_STDERR,'(/,a)')                ' ┌'//DIVIDER//'┐'
  write(formatString,'(a,i2,a)') '(a,24x,a,',max(1,panelwidth-24-len_trim(paneltype)),'x,a)'
  write(IO_STDERR,formatString)          ' │',trim(paneltype),                                      '│'
  write(formatString,'(a,i2,a)') '(a,24x,i3,',max(1,panelwidth-24-3),'x,a)'
  write(IO_STDERR,formatString)          ' │',ID,                                                   '│'
  write(IO_STDERR,'(a)')                  ' ├'//DIVIDER//'┤'
  write(formatString,'(a,i3.3,a,i3.3,a)') '(1x,a4,a',max(1,len_trim(msg)),',',&
                                                     max(1,panelwidth+3-len_trim(msg)-4),'x,a)'
  write(IO_STDERR,formatString)            '│ ',trim(msg),                                          '│'
  if (present(ext_msg)) then
    write(formatString,'(a,i3.3,a,i3.3,a)') '(1x,a4,a',max(1,len_trim(ext_msg)),',',&
                                                       max(1,panelwidth+3-len_trim(ext_msg)-4),'x,a)'
    write(IO_STDERR,formatString)          '│ ',trim(ext_msg),                                      '│'
  end if
  if (present(label1)) then
    write(formatString,'(a,i3.3,a,i3.3,a)') '(1x,a7,a',max(1,len_trim(label1)),',i9,',&
                                                       max(1,panelwidth+3-len_trim(label1)-9-7),'x,a)'
    write(IO_STDERR,formatString)          '│ at ',trim(label1),ID1,                                '│'
  end if
  if (present(label2)) then
    write(formatString,'(a,i3.3,a,i3.3,a)') '(1x,a7,a',max(1,len_trim(label2)),',i9,',&
                                                       max(1,panelwidth+3-len_trim(label2)-9-7),'x,a)'
    write(IO_STDERR,formatString)          '│ at ',trim(label2),ID2,                                '│'
  end if
  write(formatString,'(a,i2.2,a)') '(a,',max(1,panelwidth),'x,a)'
  write(IO_STDERR,formatString)          ' │',                                                     '│'
  write(IO_STDERR,'(a)')                  ' └'//DIVIDER//'┘'
  flush(IO_STDERR)
  !$OMP END CRITICAL (write2out)

end subroutine panel


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some IO functions.
!--------------------------------------------------------------------------------------------------
subroutine selfTest()

  integer, dimension(:), allocatable :: chunkPos
  character(len=:),      allocatable :: str,out


  if (dNeq(1.0_pReal, IO_stringAsFloat('1.0')))      error stop 'IO_stringAsFloat'
  if (dNeq(1.0_pReal, IO_stringAsFloat('1e0')))      error stop 'IO_stringAsFloat'
  if (dNeq(0.1_pReal, IO_stringAsFloat('1e-1')))     error stop 'IO_stringAsFloat'
  if (dNeq(0.1_pReal, IO_stringAsFloat('1.0e-1')))   error stop 'IO_stringAsFloat'
  if (dNeq(0.1_pReal, IO_stringAsFloat('1.00e-1')))  error stop 'IO_stringAsFloat'
  if (dNeq(10._pReal, IO_stringAsFloat(' 1.0e+1 '))) error stop 'IO_stringAsFloat'

  if (3112019  /= IO_stringAsInt( '3112019'))        error stop 'IO_stringAsInt'
  if (3112019  /= IO_stringAsInt(' 3112019'))        error stop 'IO_stringAsInt'
  if (-3112019 /= IO_stringAsInt('-3112019'))        error stop 'IO_stringAsInt'
  if (3112019  /= IO_stringAsInt('+3112019 '))       error stop 'IO_stringAsInt'
  if (3112019  /= IO_stringAsInt('03112019 '))       error stop 'IO_stringAsInt'
  if (3112019  /= IO_stringAsInt('+03112019'))       error stop 'IO_stringAsInt'

  if (.not. IO_stringAsBool(' true'))                error stop 'IO_stringAsBool'
  if (.not. IO_stringAsBool(' True '))               error stop 'IO_stringAsBool'
  if (      IO_stringAsBool(' false'))               error stop 'IO_stringAsBool'
  if (      IO_stringAsBool('False'))                error stop 'IO_stringAsBool'

  if ('1234' /= IO_intAsString(1234))                error stop 'IO_intAsString'
  if ('-12'  /= IO_intAsString(-0012))               error stop 'IO_intAsString'

  if (any([1,1,1]     /= IO_stringPos('a')))         error stop 'IO_stringPos'
  if (any([2,2,3,5,5] /= IO_stringPos(' aa b')))     error stop 'IO_stringPos'

  str = ' 1.0 xxx'
  chunkPos = IO_stringPos(str)
  if (dNeq(1.0_pReal,IO_floatValue(str,chunkPos,1))) error stop 'IO_floatValue'

  str = 'M 3112019 F'
  chunkPos = IO_stringPos(str)
  if (3112019 /= IO_intValue(str,chunkPos,2))        error stop 'IO_intValue'

  if (CRLF2LF('') /= '')                             error stop 'CRLF2LF/0'
  if (CRLF2LF(LF)     /= LF)                         error stop 'CRLF2LF/1a'
  if (CRLF2LF(CR//LF) /= LF)                         error stop 'CRLF2LF/1b'
  if (CRLF2LF(' '//LF)     /= ' '//LF)               error stop 'CRLF2LF/2a'
  if (CRLF2LF(' '//CR//LF) /= ' '//LF)               error stop 'CRLF2LF/2b'
  if (CRLF2LF('A'//CR//LF//'B') /= 'A'//LF//'B')     error stop 'CRLF2LF/3'
  if (CRLF2LF('A'//CR//LF//'B'//CR//LF) /= &
              'A'//LF//'B'//LF)                      error stop 'CRLF2LF/4'

  str='  ';        if (.not. IO_isBlank(str))        error stop 'IO_isBlank/1'
  str='  #isBlank';if (.not. IO_isBlank(str))        error stop 'IO_isBlank/2'
  str='  i#s';     if (      IO_isBlank(str))        error stop 'IO_isBlank/3'

  str='*(HiU!)3';if ('*(hiu!)3' /= IO_lc(str))       error stop 'IO_lc'

  str='#';out=IO_rmComment(str)
  if (out /= ''   .or. len(out) /= 0)                error stop 'IO_rmComment/1'
  str=' #';out=IO_rmComment(str)
  if (out /= ''   .or. len(out) /= 0)                error stop 'IO_rmComment/2'
  str=' # ';out=IO_rmComment(str)
  if (out /= ''   .or. len(out) /= 0)                error stop 'IO_rmComment/3'
  str=' # a';out=IO_rmComment(str)
  if (out /= ''   .or. len(out) /= 0)                error stop 'IO_rmComment/4'
  str=' a#';out=IO_rmComment(str)
  if (out /= ' a' .or. len(out) /= 2)                error stop 'IO_rmComment/5'
  str=' ab #';out=IO_rmComment(str)
  if (out /= ' ab'.or. len(out) /= 3)                error stop 'IO_rmComment/6'

end subroutine selfTest

end module IO
