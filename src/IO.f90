!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  input/output functions
!--------------------------------------------------------------------------------------------------
module IO
  use prec

  implicit none
  private

  character(len=*), parameter, public :: &
    IO_WHITESPACE = achar(44)//achar(32)//achar(9)//achar(10)//achar(13)                            !< whitespace characters
  character, parameter, public :: &
    IO_EOL = new_line('DAMASK'), &                                                                  !< end of line character
    IO_COMMENT = '#'
  character(len=*), parameter, private :: &
    IO_DIVIDER = '───────────────────'//&
                 '───────────────────'//&
                 '───────────────────'//&
                 '────────────'

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
    IO_stringAsInt, &
    IO_stringAsFloat, &
    IO_stringAsBool, &
    IO_error, &
    IO_warning

contains


!--------------------------------------------------------------------------------------------------
!> @brief do self test
!--------------------------------------------------------------------------------------------------
subroutine IO_init

  write(6,'(/,a)') ' <<<+-  IO init  -+>>>'; flush(6)

  call selfTest

end subroutine IO_init


!--------------------------------------------------------------------------------------------------
!> @brief read ASCII file and split at EOL
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

!--------------------------------------------------------------------------------------------------
! count lines to allocate string array
  N_lines = 0
  do l=1, len(rawData)
    if (rawData(l:l) == IO_EOL) N_lines = N_lines+1
  enddo
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
        call IO_warning(207,ext_msg=trim(fileName),el=l)
        warned = .true.
      endif
    else
      line = rawData(startPos:endpos)
    endif
    startPos = endPos + 2                                                                           ! jump to next line start

    fileContent(l) = trim(line)//''
    l = l + 1
  enddo

end function IO_readlines


!--------------------------------------------------------------------------------------------------
!> @brief read ASCII file into a string
!> @details ensures that the string ends with a new line (expected UNIX behavior)
!--------------------------------------------------------------------------------------------------
function IO_read(fileName) result(fileContent)

  character(len=*),  intent(in) :: fileName
  character(len=:), allocatable :: fileContent
  integer ::  &
    fileLength, &
    fileUnit, &
    myStat

  inquire(file = fileName, size=fileLength)
  open(newunit=fileUnit, file=fileName, access='stream',&
       status='old', position='rewind', action='read',iostat=myStat)
  if(myStat /= 0) call IO_error(100,ext_msg=trim(fileName))
  allocate(character(len=fileLength)::fileContent)
  if(fileLength==0) then
    close(fileUnit)
    return
  endif

  read(fileUnit,iostat=myStat) fileContent
  if(myStat /= 0) call IO_error(102,ext_msg=trim(fileName))
  close(fileUnit)

  if(fileContent(fileLength:fileLength) /= IO_EOL) fileContent = fileContent//IO_EOL                ! ensure EOL@EOF

end function IO_read


!--------------------------------------------------------------------------------------------------
!> @brief identifies strings without content
!--------------------------------------------------------------------------------------------------
logical pure function IO_isBlank(string)

  character(len=*), intent(in) :: string                                                            !< string to check for content

  integer :: posNonBlank

  posNonBlank = verify(string,IO_WHITESPACE)
  IO_isBlank = posNonBlank == 0 .or. posNonBlank == scan(string,IO_COMMENT)

end function IO_isBlank


!--------------------------------------------------------------------------------------------------
!> @brief locates all whitespace-separated chunks in given string and returns array containing
!! number them and the left/right position to be used by IO_xxxVal
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
    endif endOfString
  enddo

end function IO_stringPos


!--------------------------------------------------------------------------------------------------
!> @brief reads string value at myChunk from string
!--------------------------------------------------------------------------------------------------
function IO_stringValue(string,chunkPos,myChunk)

  character(len=*),             intent(in) :: string                                                !< raw input with known start and end of each chunk
  integer,   dimension(:),      intent(in) :: chunkPos                                              !< positions of start and end of each tag/chunk in given string
  integer,                      intent(in) :: myChunk                                               !< position number of desired chunk
  character(len=:), allocatable            :: IO_stringValue

  validChunk: if (myChunk > chunkPos(1) .or. myChunk < 1) then
    IO_stringValue = ''
    call IO_error(110,el=myChunk,ext_msg='IO_stringValue: "'//trim(string)//'"')
  else validChunk
    IO_stringValue = string(chunkPos(myChunk*2):chunkPos(myChunk*2+1))
  endif validChunk

end function IO_stringValue


!--------------------------------------------------------------------------------------------------
!> @brief reads integer value at myChunk from string
!--------------------------------------------------------------------------------------------------
integer function IO_intValue(string,chunkPos,myChunk)

  character(len=*),      intent(in) :: string                                                       !< raw input with known start and end of each chunk
  integer, dimension(:), intent(in) :: chunkPos                                                     !< positions of start and end of each tag/chunk in given string
  integer,               intent(in) :: myChunk                                                      !< position number of desired chunk

  IO_intValue = IO_stringAsInt(IO_stringValue(string,chunkPos,myChunk))

end function IO_intValue


!--------------------------------------------------------------------------------------------------
!> @brief reads float value at myChunk from string
!--------------------------------------------------------------------------------------------------
real(pReal) function IO_floatValue(string,chunkPos,myChunk)

  character(len=*),        intent(in) :: string                                                     !< raw input with known start and end of each chunk
  integer,   dimension(:), intent(in) :: chunkPos                                                   !< positions of start and end of each tag/chunk in given string
  integer,                 intent(in) :: myChunk                                                    !< position number of desired chunk

  IO_floatValue = IO_stringAsFloat(IO_stringValue(string,chunkPos,myChunk))

end function IO_floatValue


!--------------------------------------------------------------------------------------------------
!> @brief changes characters in string to lower case
!--------------------------------------------------------------------------------------------------
pure function IO_lc(string)

  character(len=*), intent(in) :: string                                                            !< string to convert
  character(len=len(string))   :: IO_lc

  character(len=*),          parameter :: LOWER = 'abcdefghijklmnopqrstuvwxyz'
  character(len=len(LOWER)), parameter :: UPPER = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  integer :: i,n

  do i=1,len(string)
    n = index(UPPER,string(i:i))
    if(n/=0) then
      IO_lc(i:i) = LOWER(n:n)
    else
      IO_lc(i:i) = string(i:i)
    endif
  enddo

end function IO_lc


!--------------------------------------------------------------------------------------------------
! @brief Remove comments (characters beyond '#') and trailing space
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
  endif

end function IO_rmComment


!--------------------------------------------------------------------------------------------------
!> @brief return verified integer value in given string
!--------------------------------------------------------------------------------------------------
integer function IO_stringAsInt(string)

  character(len=*), intent(in) :: string                                                            !< string for conversion to int value

  integer                      :: readStatus
  character(len=*), parameter  :: VALIDCHARS = '0123456789+- '

  valid: if (verify(string,VALIDCHARS) == 0) then
    read(string,*,iostat=readStatus) IO_stringAsInt
    if (readStatus /= 0) call IO_error(111,ext_msg=string)
  else valid
    IO_stringAsInt = 0
    call IO_error(111,ext_msg=string)
  endif valid

end function IO_stringAsInt


!--------------------------------------------------------------------------------------------------
!> @brief return verified float value in given string
!--------------------------------------------------------------------------------------------------
real(pReal) function IO_stringAsFloat(string)

  character(len=*), intent(in) :: string                                                            !< string for conversion to float value

  integer                      :: readStatus
  character(len=*), parameter  :: VALIDCHARS = '0123456789eE.+- '

  valid: if (verify(string,VALIDCHARS) == 0) then
    read(string,*,iostat=readStatus) IO_stringAsFloat
    if (readStatus /= 0) call IO_error(112,ext_msg=string)
  else valid
    IO_stringAsFloat = 0.0_pReal
    call IO_error(112,ext_msg=string)
  endif valid

end function IO_stringAsFloat


!--------------------------------------------------------------------------------------------------
!> @brief return verified logical value in given string
!--------------------------------------------------------------------------------------------------
logical function IO_stringAsBool(string)

  character(len=*), intent(in) :: string                                                            !< string for conversion to int value

  if     (trim(adjustl(string)) == 'True' .or.  trim(adjustl(string)) == 'true') then
    IO_stringAsBool = .true.
  elseif (trim(adjustl(string)) == 'False' .or. trim(adjustl(string)) == 'false') then
    IO_stringAsBool = .false.
  else
    IO_stringAsBool = .false.
    call IO_error(113,ext_msg=string)
  endif

end function IO_stringAsBool


!--------------------------------------------------------------------------------------------------
!> @brief write error statements to standard out and terminate the Marc/spectral run with exit #9xxx
!--------------------------------------------------------------------------------------------------
subroutine IO_error(error_ID,el,ip,g,instance,ext_msg)

  integer,                    intent(in) :: error_ID
  integer,          optional, intent(in) :: el,ip,g,instance
  character(len=*), optional, intent(in) :: ext_msg

  external                      :: quit
  character(len=:), allocatable :: msg
  character(len=pStringLen)     :: formatString

  select case (error_ID)

!--------------------------------------------------------------------------------------------------
! internal errors
    case (0)
      msg = 'internal check failed:'

!--------------------------------------------------------------------------------------------------
! file handling errors
    case (100)
      msg = 'could not open file:'
    case (101)
      msg = 'write error for file:'
    case (102)
      msg = 'could not read file:'
    case (103)
      msg = 'could not assemble input files'
    case (106)
      msg = 'working directory does not exist:'

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
    case (133)
      msg = 'transformed hex lattice structure with invalid c/a ratio'
    case (134)
      msg = 'negative lattice parameter'
    case (135)
      msg = 'zero entry on stiffness diagonal'
    case (136)
      msg = 'zero entry on stiffness diagonal for transformed phase'
    case (137)
      msg = 'not defined for lattice structure'
    case (138)
      msg = 'not enough interaction parameters given'

!--------------------------------------------------------------------------------------------------
! errors related to the parsing of material.config
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
    case (148)
      msg = 'Nconstituents mismatch between homogenization and microstructure'

!--------------------------------------------------------------------------------------------------
! material error messages and related messages in mesh
    case (150)
      msg = 'index out of bounds'
    case (151)
      msg = 'microstructure has no constituents'
    case (153)
      msg = 'sum of phase fractions differs from 1'
    case (155)
      msg = 'microstructure index out of bounds'
    case (180)
      msg = 'missing/invalid microstructure definition via State Variable 2'
    case (190)
      msg = 'unknown element type:'
    case (191)
      msg = 'mesh consists of more than one element type'

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
    case (400)
      msg = 'matrix inversion error'
    case (401)
      msg = 'error in Eigenvalue calculation'
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

!------------------------------------------------------------------------------------------------
! errors related to YAML data
    case (701)
      msg = 'Incorrect indent/Null value not allowed'
    case (702)
      msg = 'Invalid use of flow yaml'
    case (704)
      msg = 'Space expected after a colon for <key>: <value> pair'
    case (705)
      msg = 'Unsupported feature'
    case (706)
      msg = 'Access by incorrect node type'

!-------------------------------------------------------------------------------------------------
! errors related to the grid solver
    case (809)
      msg = 'initializing FFTW'
    case (810)
      msg = 'FFTW plan creation'
    case (831)
      msg = 'mask consistency violated in grid load case'
    case (832)
      msg = 'ill-defined L (line partly defined) in grid load case'
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
    case (841)
      msg = 'missing header length info in grid mesh'
    case (842)
      msg = 'incomplete information in grid mesh header'
    case (843)
      msg = 'microstructure count mismatch'
    case (844)
      msg = 'invalid VTR file'
    case (846)
      msg = 'rotation for load case rotation ill-defined (R:RT != I)'
    case (891)
      msg = 'unknown solver type selected'
    case (892)
      msg = 'unknown filter type selected'
    case (894)
      msg = 'MPI error'


!-------------------------------------------------------------------------------------------------
! general error messages
    case default
      msg = 'unknown error number...'

  end select

  !$OMP CRITICAL (write2out)
  write(0,'(/,a)')                ' ┌'//IO_DIVIDER//'┐'
  write(0,'(a,24x,a,40x,a)')      ' │','error',                                             '│'
  write(0,'(a,24x,i3,42x,a)')     ' │',error_ID,                                            '│'
  write(0,'(a)')                  ' ├'//IO_DIVIDER//'┤'
  write(formatString,'(a,i6.6,a,i6.6,a)') '(1x,a4,a',max(1,len_trim(msg)),',',&
                                                     max(1,72-len_trim(msg)-4),'x,a)'
  write(0,formatString)            '│ ',trim(msg),                                          '│'
  if (present(ext_msg)) then
    write(formatString,'(a,i6.6,a,i6.6,a)') '(1x,a4,a',max(1,len_trim(ext_msg)),',',&
                                                       max(1,72-len_trim(ext_msg)-4),'x,a)'
    write(0,formatString)          '│ ',trim(ext_msg),                                      '│'
  endif
  if (present(el)) &
    write(0,'(a19,1x,i9,44x,a3)') ' │ at element    ',el,                                   '│'
  if (present(ip)) &
    write(0,'(a19,1x,i9,44x,a3)') ' │ at IP         ',ip,                                   '│'
  if (present(g)) &
    write(0,'(a19,1x,i9,44x,a3)') ' │ at constituent',g,                                    '│'
  if (present(instance)) &
    write(0,'(a19,1x,i9,44x,a3)') ' │ at instance   ',instance,                             '│'
  write(0,'(a,69x,a)')            ' │',                                                     '│'
  write(0,'(a)')                  ' └'//IO_DIVIDER//'┘'
  flush(0)
  call quit(9000+error_ID)
  !$OMP END CRITICAL (write2out)

end subroutine IO_error


!--------------------------------------------------------------------------------------------------
!> @brief writes warning statement to standard out
!--------------------------------------------------------------------------------------------------
subroutine IO_warning(warning_ID,el,ip,g,ext_msg)

  integer,                    intent(in) :: warning_ID
  integer,          optional, intent(in) :: el,ip,g
  character(len=*), optional, intent(in) :: ext_msg

  character(len=:), allocatable :: msg
  character(len=pStringLen)     :: formatString

  select case (warning_ID)
    case (1)
      msg = 'unknown key'
    case (34)
      msg = 'invalid restart increment given'
    case (35)
      msg = 'could not get $DAMASK_NUM_THREADS'
    case (40)
      msg = 'found spectral solver parameter'
    case (42)
      msg = 'parameter has no effect'
    case (43)
      msg = 'main diagonal of C66 close to zero'
    case (47)
      msg = 'no valid parameter for FFTW, using FFTW_PATIENT'
    case (50)
      msg = 'not all available slip system families are defined'
    case (51)
      msg = 'not all available twin system families are defined'
    case (52)
      msg = 'not all available parameters are defined'
    case (53)
      msg = 'not all available transformation system families are defined'
    case (101)
      msg = 'crystallite debugging off'
    case (201)
      msg = 'position not found when parsing line'
    case (207)
      msg = 'line truncated'
    case (600)
      msg = 'crystallite responds elastically'
    case (601)
      msg = 'stiffness close to zero'
    case (650)
      msg = 'polar decomposition failed'
    case (700)
      msg = 'unknown crystal symmetry'
    case (850)
      msg = 'max number of cut back exceeded, terminating'
    case default
      msg = 'unknown warning number'
    end select

  !$OMP CRITICAL (write2out)
  write(6,'(/,a)')                ' ┌'//IO_DIVIDER//'┐'
  write(6,'(a,24x,a,38x,a)')      ' │','warning',                                           '│'
  write(6,'(a,24x,i3,42x,a)')     ' │',warning_ID,                                          '│'
  write(6,'(a)')                  ' ├'//IO_DIVIDER//'┤'
  write(formatString,'(a,i6.6,a,i6.6,a)') '(1x,a4,a',max(1,len_trim(msg)),',',&
                                                     max(1,72-len_trim(msg)-4),'x,a)'
  write(6,formatString)            '│ ',trim(msg),                                          '│'
  if (present(ext_msg)) then
    write(formatString,'(a,i6.6,a,i6.6,a)') '(1x,a4,a',max(1,len_trim(ext_msg)),',',&
                                                       max(1,72-len_trim(ext_msg)-4),'x,a)'
    write(6,formatString)          '│ ',trim(ext_msg),                                      '│'
  endif
  if (present(el)) &
    write(6,'(a19,1x,i9,44x,a3)') ' │ at element    ',el,                                   '│'
  if (present(ip)) &
    write(6,'(a19,1x,i9,44x,a3)') ' │ at IP         ',ip,                                   '│'
  if (present(g)) &
    write(6,'(a19,1x,i9,44x,a3)') ' │ at constituent',g,                                    '│'
  write(6,'(a,69x,a)')            ' │',                                                     '│'
  write(6,'(a)')                  ' └'//IO_DIVIDER//'┘'
  flush(6)
  !$OMP END CRITICAL (write2out)

end subroutine IO_warning


!--------------------------------------------------------------------------------------------------
!> @brief check correctness of some IO functions
!--------------------------------------------------------------------------------------------------
subroutine selfTest

  integer, dimension(:), allocatable :: chunkPos
  character(len=:),      allocatable :: str

  if(dNeq(1.0_pReal, IO_stringAsFloat('1.0')))      call IO_error(0,ext_msg='IO_stringAsFloat')
  if(dNeq(1.0_pReal, IO_stringAsFloat('1e0')))      call IO_error(0,ext_msg='IO_stringAsFloat')
  if(dNeq(0.1_pReal, IO_stringAsFloat('1e-1')))     call IO_error(0,ext_msg='IO_stringAsFloat')

  if(3112019  /= IO_stringAsInt( '3112019'))        call IO_error(0,ext_msg='IO_stringAsInt')
  if(3112019  /= IO_stringAsInt(' 3112019'))        call IO_error(0,ext_msg='IO_stringAsInt')
  if(-3112019 /= IO_stringAsInt('-3112019'))        call IO_error(0,ext_msg='IO_stringAsInt')
  if(3112019  /= IO_stringAsInt('+3112019 '))       call IO_error(0,ext_msg='IO_stringAsInt')

  if(.not. IO_stringAsBool(' true'))                call IO_error(0,ext_msg='IO_stringAsBool')
  if(.not. IO_stringAsBool(' True '))               call IO_error(0,ext_msg='IO_stringAsBool')
  if(      IO_stringAsBool(' false'))               call IO_error(0,ext_msg='IO_stringAsBool')
  if(      IO_stringAsBool('False'))                call IO_error(0,ext_msg='IO_stringAsBool')

  if(any([1,1,1]     /= IO_stringPos('a')))         call IO_error(0,ext_msg='IO_stringPos')
  if(any([2,2,3,5,5] /= IO_stringPos(' aa b')))     call IO_error(0,ext_msg='IO_stringPos')

  str=' 1.0 xxx'
  chunkPos = IO_stringPos(str)
  if(dNeq(1.0_pReal,IO_floatValue(str,chunkPos,1))) call IO_error(0,ext_msg='IO_floatValue')

  str='M 3112019 F'
  chunkPos = IO_stringPos(str)
  if(3112019 /= IO_intValue(str,chunkPos,2))        call IO_error(0,ext_msg='IO_intValue')

  if(.not. IO_isBlank('  '))                        call IO_error(0,ext_msg='IO_isBlank/1')
  if(.not. IO_isBlank('  #isBlank'))                call IO_error(0,ext_msg='IO_isBlank/2')
  if(      IO_isBlank('  i#s'))                     call IO_error(0,ext_msg='IO_isBlank/3')

  str = IO_rmComment('#')
  if (str /= ''   .or. len(str) /= 0)               call IO_error(0,ext_msg='IO_rmComment/1')
  str = IO_rmComment(' #')
  if (str /= ''   .or. len(str) /= 0)               call IO_error(0,ext_msg='IO_rmComment/2')
  str = IO_rmComment(' # ')
  if (str /= ''   .or. len(str) /= 0)               call IO_error(0,ext_msg='IO_rmComment/3')
  str = IO_rmComment(' # a')
  if (str /= ''   .or. len(str) /= 0)               call IO_error(0,ext_msg='IO_rmComment/4')
  str = IO_rmComment(' # a')
  if (str /= ''   .or. len(str) /= 0)               call IO_error(0,ext_msg='IO_rmComment/5')
  str = IO_rmComment(' a#')
  if (str /= ' a' .or. len(str) /= 2)               call IO_error(0,ext_msg='IO_rmComment/6')
  str = IO_rmComment(' ab #')
  if (str /= ' ab'.or. len(str) /= 3)               call IO_error(0,ext_msg='IO_rmComment/7')

end subroutine selfTest

end module IO
