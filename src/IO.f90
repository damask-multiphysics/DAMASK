!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  input/output functions, partly depending on chosen solver
!--------------------------------------------------------------------------------------------------
module IO
  use prec
  use DAMASK_interface
  
  implicit none
  private
  character(len=*), parameter, public :: &
    IO_EOF = '#EOF#'                                                                                !< end of file string
  character, parameter, public :: &
    IO_EOL = new_line(' ')                                                                          !< end of line str
  character(len=*), parameter, private :: &
    IO_DIVIDER = '───────────────────'//&
                 '───────────────────'//&
                 '───────────────────'//&
                 '────────────'
  public :: &
    IO_init, &
    IO_read_ASCII, &
    IO_open_file, &                                                                                 ! deprecated, use IO_read_ASCII
    IO_open_jobFile_binary, &
    IO_isBlank, &
    IO_getTag, &
    IO_stringPos, &
    IO_stringValue, &
    IO_floatValue, &
    IO_intValue, &
    IO_lc, &
    IO_error, &
    IO_warning
#if defined(Marc4DAMASK)
  public :: &
    IO_open_inputFile
#endif

contains


!--------------------------------------------------------------------------------------------------
!> @brief does nothing.
! ToDo: needed?
!--------------------------------------------------------------------------------------------------
subroutine IO_init
 
  write(6,'(/,a)') ' <<<+-  IO init  -+>>>'; flush(6)
 
end subroutine IO_init


!--------------------------------------------------------------------------------------------------
!> @brief reads an entire ASCII file into an array
!--------------------------------------------------------------------------------------------------
function IO_read_ASCII(fileName) result(fileContent)

  character(len=*),          intent(in)                :: fileName

  character(len=pStringLen), dimension(:), allocatable :: fileContent                               !< file content, separated per lines
  character(len=pStringLen)                            :: line
  character(len=:),                        allocatable :: rawData
  integer ::  &
    fileLength, &
    fileUnit, &
    startPos, endPos, &
    myTotalLines, &                                                                                 !< # lines read from file
    l, &
    myStat
  logical :: warned
  
!--------------------------------------------------------------------------------------------------
! read data as stream
  inquire(file = fileName, size=fileLength)
  if (fileLength == 0) then
    allocate(fileContent(0))
    return
  endif
  open(newunit=fileUnit, file=fileName, access='stream',&
       status='old', position='rewind', action='read',iostat=myStat)
  if(myStat /= 0) call IO_error(100,ext_msg=trim(fileName))
  allocate(character(len=fileLength)::rawData)
  read(fileUnit) rawData
  close(fileUnit)

!--------------------------------------------------------------------------------------------------
! count lines to allocate string array
  myTotalLines = 1
  do l=1, len(rawData)
    if (rawData(l:l) == IO_EOL) myTotalLines = myTotalLines+1
  enddo
  allocate(fileContent(myTotalLines))

!--------------------------------------------------------------------------------------------------
! split raw data at end of line
  warned = .false.
  startPos = 1
  l = 1
  do while (l <= myTotalLines)
    endPos = merge(startPos + scan(rawData(startPos:),IO_EOL) - 2,len(rawData),l /= myTotalLines)
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

    fileContent(l) = line
    l = l + 1
  enddo

end function IO_read_ASCII


!--------------------------------------------------------------------------------------------------
!> @brief   opens existing file for reading to given unit. Path to file is relative to working
!!          directory
!--------------------------------------------------------------------------------------------------
subroutine IO_open_file(fileUnit,path)
 
  integer,            intent(in) :: fileUnit                                                        !< file unit
  character(len=*),   intent(in) :: path                                                            !< relative path from working directory
 
  integer                        :: myStat
 
  open(fileUnit,status='old',iostat=myStat,file=path,action='read',position='rewind')
  if (myStat /= 0) call IO_error(100,el=myStat,ext_msg=path)
 
end subroutine IO_open_file


!--------------------------------------------------------------------------------------------------
!> @brief opens an existing file for reading or a new file for writing. Name is the job name
!> @details replaces an existing file when writing
!--------------------------------------------------------------------------------------------------
integer function IO_open_jobFile_binary(extension,mode)

  character(len=*), intent(in)           :: extension
  character,        intent(in), optional :: mode
 
  if (present(mode)) then
    IO_open_jobFile_binary = IO_open_binary(trim(getSolverJobName())//'.'//trim(extension),mode)
  else
    IO_open_jobFile_binary = IO_open_binary(trim(getSolverJobName())//'.'//trim(extension))
  endif

end function IO_open_jobFile_binary


!--------------------------------------------------------------------------------------------------
!> @brief opens an existing file for reading or a new file for writing.
!> @details replaces an existing file when writing
!--------------------------------------------------------------------------------------------------
integer function IO_open_binary(fileName,mode)

  character(len=*), intent(in)           :: fileName
  character,        intent(in), optional :: mode
 
  character :: m
  integer   :: ierr 

  if (present(mode)) then
    m = mode
  else
    m = 'r'
  endif

  if    (m == 'w') then
    open(newunit=IO_open_binary, file=trim(fileName),&
         status='replace',access='stream',action='write',iostat=ierr)
    if (ierr /= 0) call IO_error(100,ext_msg='could not open file (w): '//trim(fileName))
  elseif(m == 'r') then
    open(newunit=IO_open_binary, file=trim(fileName),&
         status='old',    access='stream',action='read', iostat=ierr)
    if (ierr /= 0) call IO_error(100,ext_msg='could not open file (r): '//trim(fileName))
  else
    call IO_error(100,ext_msg='unknown access mode: '//m)
  endif

end function IO_open_binary


!--------------------------------------------------------------------------------------------------
!> @brief identifies strings without content
!--------------------------------------------------------------------------------------------------
logical pure function IO_isBlank(string)

  character(len=*), intent(in) :: string                                                            !< string to check for content

  character(len=*),  parameter :: blankChar = achar(32)//achar(9)//achar(10)//achar(13)             ! whitespaces
  character(len=*),  parameter :: comment = achar(35)                                               ! comment id '#'

  integer :: posNonBlank, posComment

  posNonBlank = verify(string,blankChar)
  posComment  = scan(string,comment)
  IO_isBlank = posNonBlank == 0 .or. posNonBlank == posComment

end function IO_isBlank


!--------------------------------------------------------------------------------------------------
!> @brief get tagged content of string
!--------------------------------------------------------------------------------------------------
pure function IO_getTag(string,openChar,closeChar)

  character(len=*), intent(in)  :: string                                                           !< string to check for tag
  character(len=len_trim(string)) :: IO_getTag
 
  character, intent(in)  :: openChar, &                                                             !< indicates beginning of tag
                            closeChar                                                               !< indicates end of tag
 
  character(len=*), parameter   :: SEP=achar(32)//achar(9)//achar(10)//achar(13)                    ! whitespaces
  integer :: left,right
 
  IO_getTag = ''
 
 
  if (openChar /= closeChar) then
    left  = scan(string,openChar)
    right = scan(string,closeChar)
  else
    left  = scan(string,openChar)
    right = left + merge(scan(string(left+1:),openChar),0,len(string) > left)
  endif
 
  if (left == verify(string,SEP) .and. right > left) &                                              ! openChar is first and closeChar occurs
    IO_getTag = string(left+1:right-1)

end function IO_getTag


!--------------------------------------------------------------------------------------------------
!> @brief locates all space-separated chunks in given string and returns array containing number
!! them and the left/right position to be used by IO_xxxVal
!! Array size is dynamically adjusted to number of chunks found in string
!! IMPORTANT: first element contains number of chunks!
!--------------------------------------------------------------------------------------------------
pure function IO_stringPos(string)

 integer, dimension(:), allocatable            :: IO_stringPos
 character(len=*),                  intent(in) :: string                                            !< string in which chunk positions are searched for

 character(len=*), parameter  :: SEP=achar(44)//achar(32)//achar(9)//achar(10)//achar(13)           ! comma and whitespaces
 integer                      :: left, right

 allocate(IO_stringPos(1), source=0)
 right = 0

 do while (verify(string(right+1:),SEP)>0)
   left  = right + verify(string(right+1:),SEP)
   right = left + scan(string(left:),SEP) - 2
   if ( string(left:left) == '#' ) exit
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
function IO_stringValue(string,chunkPos,myChunk,silent)

  integer,   dimension(:),                intent(in) :: chunkPos                                    !< positions of start and end of each tag/chunk in given string
  integer,                                intent(in) :: myChunk                                     !< position number of desired chunk
  character(len=*),                       intent(in) :: string                                      !< raw input with known start and end of each chunk
  character(len=:), allocatable                      :: IO_stringValue

  logical,                       optional,intent(in) :: silent                                      !< switch to trigger verbosity
  character(len=*),  parameter                       :: MYNAME = 'IO_stringValue: '

  logical                                            :: warn

  if (present(silent)) then
    warn = .not. silent
  else
    warn = .false.
  endif

  IO_stringValue = ''
  valuePresent: if (myChunk > chunkPos(1) .or. myChunk < 1) then
    if (warn) call IO_warning(201,el=myChunk,ext_msg=MYNAME//trim(string))
  else valuePresent
    IO_stringValue = string(chunkPos(myChunk*2):chunkPos(myChunk*2+1))
  endif valuePresent

end function IO_stringValue


!--------------------------------------------------------------------------------------------------
!> @brief reads float value at myChunk from string
!--------------------------------------------------------------------------------------------------
real(pReal) function IO_floatValue(string,chunkPos,myChunk)

  integer,   dimension(:),        intent(in) :: chunkPos                                            !< positions of start and end of each tag/chunk in given string
  integer,                        intent(in) :: myChunk                                             !< position number of desired chunk
  character(len=*),               intent(in) :: string                                              !< raw input with known start and end of each chunk
  character(len=*),               parameter  :: MYNAME = 'IO_floatValue: '
  character(len=*),               parameter  :: VALIDCHARACTERS = '0123456789eEdD.+-'

  IO_floatValue = 0.0_pReal

  valuePresent: if (myChunk > chunkPos(1) .or. myChunk < 1) then
    call IO_warning(201,el=myChunk,ext_msg=MYNAME//trim(string))
  else valuePresent
    IO_floatValue = verifyFloatValue(trim(adjustl(string(chunkPos(myChunk*2):chunkPos(myChunk*2+1)))),&
                                     VALIDCHARACTERS,MYNAME)
  endif valuePresent

end function IO_floatValue


!--------------------------------------------------------------------------------------------------
!> @brief reads integer value at myChunk from string
!--------------------------------------------------------------------------------------------------
integer function IO_intValue(string,chunkPos,myChunk)

  character(len=*),      intent(in) :: string                                                       !< raw input with known start and end of each chunk
  integer,               intent(in) :: myChunk                                                      !< position number of desired chunk
  integer, dimension(:), intent(in) :: chunkPos                                                     !< positions of start and end of each tag/chunk in given string
  character(len=*),      parameter  :: MYNAME = 'IO_intValue: '
  character(len=*),      parameter  :: VALIDCHARACTERS = '0123456789+-'

  IO_intValue = 0

  valuePresent: if (myChunk > chunkPos(1) .or. myChunk < 1) then
    call IO_warning(201,el=myChunk,ext_msg=MYNAME//trim(string))
  else valuePresent
    IO_intValue = verifyIntValue(trim(adjustl(string(chunkPos(myChunk*2):chunkPos(myChunk*2+1)))),&
                                    VALIDCHARACTERS,MYNAME)
  endif valuePresent

end function IO_intValue


!--------------------------------------------------------------------------------------------------
!> @brief changes characters in string to lower case
!--------------------------------------------------------------------------------------------------
pure function IO_lc(string)

  character(len=*), intent(in) :: string                                                            !< string to convert
  character(len=len(string))   :: IO_lc

  character(26), parameter :: LOWER = 'abcdefghijklmnopqrstuvwxyz'
  character(26), parameter :: UPPER = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  integer                  :: i,n

  do i=1,len(string)
    IO_lc(i:i) = string(i:i)
    n = index(UPPER,IO_lc(i:i))
    if (n/=0) IO_lc(i:i) = LOWER(n:n)
  enddo

end function IO_lc


!--------------------------------------------------------------------------------------------------
!> @brief write error statements to standard out and terminate the Marc/spectral run with exit #9xxx
!--------------------------------------------------------------------------------------------------
subroutine IO_error(error_ID,el,ip,g,instance,ext_msg)

  integer,                    intent(in) :: error_ID
  integer,          optional, intent(in) :: el,ip,g,instance
  character(len=*), optional, intent(in) :: ext_msg

  external                               :: quit
  character(len=pStringLen)              :: msg
  character(len=pStringLen)              :: formatString

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
    case (147)
      msg = 'not supported anymore'

!--------------------------------------------------------------------------------------------------
! material error messages and related messages in mesh
    case (150)
      msg = 'index out of bounds'
    case (151)
      msg = 'microstructure has no constituents'
    case (153)
      msg = 'sum of phase fractions differs from 1'
    case (154)
      msg = 'homogenization index out of bounds'
    case (155)
      msg = 'microstructure index out of bounds'
    case (157)
      msg = 'invalid texture transformation specified'
    case (160)
      msg = 'no entries in config part'
    case (161)
      msg = 'config part found twice'
    case (165)
      msg = 'homogenization configuration'
    case (170)
      msg = 'no homogenization specified via State Variable 2'
    case (180)
      msg = 'no microstructure specified via State Variable 3'
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

    case (210)
      msg = 'unknown material parameter:'
    case (211)
      msg = 'material parameter out of bounds:'

!--------------------------------------------------------------------------------------------------
! numerics error messages
    case (300)
      msg = 'unknown numerics parameter:'
    case (301)
      msg = 'numerics parameter out of bounds:'

!--------------------------------------------------------------------------------------------------
! math errors
    case (400)
      msg = 'matrix inversion error'
    case (401)
      msg = 'math_check failed'
    case (402)
      msg = 'invalid orientation specified'

!-------------------------------------------------------------------------------------------------
! homogenization errors
    case (500)
      msg = 'unknown homogenization specified'

!--------------------------------------------------------------------------------------------------
! user errors
    case (600)
      msg = 'Ping-Pong not possible when using non-DAMASK elements'
    case (601)
      msg = 'Ping-Pong needed when using non-local plasticity'
    case (602)
      msg = 'invalid selection for debug'

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
  write(formatString,'(a,i6.6,a,i6.6,a)') '(1x,a4,a',max(1,len(trim(msg))),',',&
                                                     max(1,72-len(trim(msg))-4),'x,a)'
  write(0,formatString)            '│ ',trim(msg),                                          '│'
  if (present(ext_msg)) then
    write(formatString,'(a,i6.6,a,i6.6,a)') '(1x,a4,a',max(1,len(trim(ext_msg))),',',&
                                                       max(1,72-len(trim(ext_msg))-4),'x,a)'
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
 
  character(len=pStringLen)              :: msg
  character(len=pStringLen)              :: formatString
 
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
    case (202)
      msg = 'invalid character in string chunk'
    case (203)
      msg = 'interpretation of string chunk failed'
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
  write(formatString,'(a,i6.6,a,i6.6,a)') '(1x,a4,a',max(1,len(trim(msg))),',',&
                                                     max(1,72-len(trim(msg))-4),'x,a)'
  write(6,formatString)            '│ ',trim(msg),                                          '│'
  if (present(ext_msg)) then
    write(formatString,'(a,i6.6,a,i6.6,a)') '(1x,a4,a',max(1,len(trim(ext_msg))),',',&
                                                       max(1,72-len(trim(ext_msg))-4),'x,a)'
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
! internal helper functions

!--------------------------------------------------------------------------------------------------
!> @brief returns verified integer value in given string
!--------------------------------------------------------------------------------------------------
integer function verifyIntValue(string,validChars,myName)
 
  character(len=*), intent(in) :: string, &                                                         !< string for conversion to int value. Must not contain spaces!
                                  validChars, &                                                     !< valid characters in string
                                  myName                                                            !< name of caller function (for debugging)
  integer                      :: readStatus, invalidWhere
 
  verifyIntValue = 0
 
  invalidWhere = verify(string,validChars)
  if (invalidWhere == 0) then
    read(string,*,iostat=readStatus) verifyIntValue                                                ! no offending chars found
    if (readStatus /= 0) &                                                                          ! error during string to integer conversion
      call IO_warning(203,ext_msg=myName//'"'//string//'"')
  else
    call IO_warning(202,ext_msg=myName//'"'//string//'"')                                           ! complain about offending characters
    read(string(1:invalidWhere-1),*,iostat=readStatus) verifyIntValue                               ! interpret remaining string
    if (readStatus /= 0) &                                                                          ! error during string to integer conversion
      call IO_warning(203,ext_msg=myName//'"'//string(1:invalidWhere-1)//'"')
  endif
 
end function verifyIntValue


!--------------------------------------------------------------------------------------------------
!> @brief returns verified float value in given string
!--------------------------------------------------------------------------------------------------
real(pReal) function verifyFloatValue(string,validChars,myName)
 
  character(len=*), intent(in) :: string, &                                                         !< string for conversion to int value. Must not contain spaces!
                                  validChars, &                                                     !< valid characters in string
                                  myName                                                            !< name of caller function (for debugging)
 
  integer                      :: readStatus, invalidWhere
 
  verifyFloatValue = 0.0_pReal
 
  invalidWhere = verify(string,validChars)
  if (invalidWhere == 0) then
    read(string,*,iostat=readStatus) verifyFloatValue                                               ! no offending chars found
    if (readStatus /= 0) &                                                                          ! error during string to float conversion
      call IO_warning(203,ext_msg=myName//'"'//string//'"')
  else
    call IO_warning(202,ext_msg=myName//'"'//string//'"')                                           ! complain about offending characters
    read(string(1:invalidWhere-1),*,iostat=readStatus) verifyFloatValue                             ! interpret remaining string
    if (readStatus /= 0) &                                                                          ! error during string to float conversion
      call IO_warning(203,ext_msg=myName//'"'//string(1:invalidWhere-1)//'"')
  endif
 
end function verifyFloatValue

end module IO
