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
  character(len=5), parameter, public :: &
    IO_EOF = '#EOF#'                                                                                !< end of file string
  character(len=207), parameter, private :: &
    IO_DIVIDER = '───────────────────'//&
                 '───────────────────'//&
                 '───────────────────'//&
                 '────────────'
  public :: &
    IO_init, &
    IO_read_ASCII, &
    IO_open_file, &
    IO_open_jobFile_binary, &
    IO_write_jobFile, &
    IO_isBlank, &
    IO_getTag, &
    IO_stringPos, &
    IO_stringValue, &
    IO_floatValue, &
    IO_intValue, &
    IO_lc, &
    IO_error, &
    IO_warning, &
    IO_intOut
#if defined(Marc4DAMASK) || defined(Abaqus)
  public :: &
    IO_open_inputFile, &
    IO_open_logFile, &
    IO_countContinuousIntValues, &
    IO_continuousIntValues, &
#if defined(Abaqus)
    IO_extractValue, &
    IO_countDataLines
#elif defined(Marc4DAMASK)
    IO_fixedNoEFloatValue, &
    IO_fixedIntValue, &
    IO_countNumericalDataLines
#endif
#endif

 private :: &
   IO_verifyFloatValue, &
   IO_verifyIntValue

contains


!--------------------------------------------------------------------------------------------------
!> @brief does nothing.
! ToDo: needed?
!--------------------------------------------------------------------------------------------------
subroutine IO_init
 
  write(6,'(/,a)')   ' <<<+-  IO init  -+>>>'
 
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
    if (rawData(l:l) == new_line('')) myTotalLines = myTotalLines+1
  enddo
  allocate(fileContent(myTotalLines))

!--------------------------------------------------------------------------------------------------
! split raw data at end of line
  warned = .false.
  startPos = 1
  l = 1
  do while (l <= myTotalLines)
    endPos = merge(startPos + scan(rawData(startPos:),new_line('')) - 2,len(rawData),l /= myTotalLines)
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


#if defined(Marc4DAMASK) || defined(Abaqus)
!--------------------------------------------------------------------------------------------------
!> @brief opens FEM input file for reading located in current working directory to given unit
!--------------------------------------------------------------------------------------------------
subroutine IO_open_inputFile(fileUnit)

  integer,          intent(in) :: fileUnit                                                          !< file unit
 
  integer                  :: myStat
  character(len=1024)      :: path
#if defined(Abaqus)
  integer                  :: fileType

  fileType = 1                                                                                      ! assume .pes
  path = trim(getSolverJobName())//inputFileExtension(fileType)                                     ! attempt .pes, if it exists: it should be used
  open(fileUnit+1,status='old',iostat=myStat,file=path,action='read',position='rewind')
  if(myStat /= 0) then                                                                              ! if .pes does not work / exist; use conventional extension, i.e.".inp"
     fileType = 2
     path = trim(getSolverJobName())//inputFileExtension(fileType)
     open(fileUnit+1,status='old',iostat=myStat,file=path,action='read',position='rewind')
  endif
  if (myStat /= 0) call IO_error(100,el=myStat,ext_msg=path)
 
  path = trim(getSolverJobName())//inputFileExtension(fileType)//'_assembly'
  open(fileUnit,iostat=myStat,file=path)
  if (myStat /= 0) call IO_error(100,el=myStat,ext_msg=path)
     if (.not.abaqus_assembleInputFile(fileUnit,fileUnit+1)) call IO_error(103)                     ! strip comments and concatenate any "include"s
  close(fileUnit+1)
  
  contains
  
  !--------------------------------------------------------------------------------------------------
  !> @brief create a new input file for abaqus simulations by removing all comment lines and
  !> including "include"s
  !--------------------------------------------------------------------------------------------------
  recursive function abaqus_assembleInputFile(unit1,unit2) result(createSuccess)
  
    integer, intent(in)                :: unit1, &
                                          unit2
   
   
    integer, allocatable, dimension(:) :: chunkPos
    character(len=65536)               :: line,fname
    logical                            :: createSuccess,fexist
   
   
    do
      read(unit2,'(A65536)',END=220) line
      chunkPos = IO_stringPos(line)
   
      if (IO_lc(IO_StringValue(line,chunkPos,1))=='*include') then
        fname = trim(line(9+scan(line(9:),'='):))
        inquire(file=fname, exist=fexist)
        if (.not.(fexist)) then
          write(6,*)'ERROR: file does not exist error in abaqus_assembleInputFile'
          write(6,*)'filename: ', trim(fname)
          createSuccess = .false.
          return
        endif
        open(unit2+1,err=200,status='old',file=fname)
        if (abaqus_assembleInputFile(unit1,unit2+1)) then
          createSuccess=.true.
          close(unit2+1)
        else
          createSuccess=.false.
          return
        endif
      else if (line(1:2) /= '**' .OR. line(1:8)=='**damask') then
        write(unit1,'(A)') trim(line)
      endif
    enddo
  
220 createSuccess = .true.
    return
  
200 createSuccess =.false.
  
  end function abaqus_assembleInputFile
#elif defined(Marc4DAMASK)
  path = trim(getSolverJobName())//inputFileExtension
  open(fileUnit,status='old',iostat=myStat,file=path)
  if (myStat /= 0) call IO_error(100,el=myStat,ext_msg=path)
#endif

end subroutine IO_open_inputFile


!--------------------------------------------------------------------------------------------------
!> @brief opens existing FEM log file for reading to given unit. File is named after solver job
!!        name and located in current working directory
!--------------------------------------------------------------------------------------------------
subroutine IO_open_logFile(fileUnit)

  integer,      intent(in) :: fileUnit                                                              !< file unit

  integer                  :: myStat
  character(len=1024)      :: path

  path = trim(getSolverJobName())//LogFileExtension
  open(fileUnit,status='old',iostat=myStat,file=path,action='read',position='rewind')
  if (myStat /= 0) call IO_error(100,el=myStat,ext_msg=path)

end subroutine IO_open_logFile
#endif


!--------------------------------------------------------------------------------------------------
!> @brief opens ASCII file to given unit for writing. File is named after solver job name plus
!!        given extension and located in current working directory
!--------------------------------------------------------------------------------------------------
subroutine IO_write_jobFile(fileUnit,ext)

  integer,          intent(in) :: fileUnit                                                          !< file unit
  character(len=*), intent(in) :: ext                                                               !< extension of file

  integer                      :: myStat
  character(len=1024)          :: path

  path = trim(getSolverJobName())//'.'//ext
  open(fileUnit,status='replace',iostat=myStat,file=path)
  if (myStat /= 0) call IO_error(100,el=myStat,ext_msg=path)

end subroutine IO_write_jobFile


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
  character(len=16), parameter                       :: MYNAME = 'IO_stringValue: '

  logical                                            :: warn

  if (present(silent)) then
    warn = silent
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
real(pReal) function IO_floatValue (string,chunkPos,myChunk)

  integer,   dimension(:),        intent(in) :: chunkPos                                            !< positions of start and end of each tag/chunk in given string
  integer,                        intent(in) :: myChunk                                             !< position number of desired chunk
  character(len=*),               intent(in) :: string                                              !< raw input with known start and end of each chunk
  character(len=15),              parameter  :: MYNAME = 'IO_floatValue: '
  character(len=17),              parameter  :: VALIDCHARACTERS = '0123456789eEdD.+-'

  IO_floatValue = 0.0_pReal

  valuePresent: if (myChunk > chunkPos(1) .or. myChunk < 1) then
    call IO_warning(201,el=myChunk,ext_msg=MYNAME//trim(string))
  else  valuePresent
    IO_floatValue = &
                IO_verifyFloatValue(trim(adjustl(string(chunkPos(myChunk*2):chunkPos(myChunk*2+1)))),&
                                        VALIDCHARACTERS,MYNAME)
  endif  valuePresent

end function IO_floatValue


!--------------------------------------------------------------------------------------------------
!> @brief reads integer value at myChunk from string
!--------------------------------------------------------------------------------------------------
integer function IO_intValue(string,chunkPos,myChunk)

  character(len=*),      intent(in) :: string                                                       !< raw input with known start and end of each chunk
  integer,               intent(in) :: myChunk                                                      !< position number of desired chunk
  integer, dimension(:), intent(in) :: chunkPos                                                     !< positions of start and end of each tag/chunk in given string
  character(len=13),     parameter  :: MYNAME = 'IO_intValue: '
  character(len=12),     parameter  :: VALIDCHARACTERS = '0123456789+-'

  IO_intValue = 0

  valuePresent: if (myChunk > chunkPos(1) .or. myChunk < 1) then
    call IO_warning(201,el=myChunk,ext_msg=MYNAME//trim(string))
  else valuePresent
    IO_intValue = IO_verifyIntValue(trim(adjustl(string(chunkPos(myChunk*2):chunkPos(myChunk*2+1)))),&
                                    VALIDCHARACTERS,MYNAME)
  endif valuePresent

end function IO_intValue


#ifdef Marc4DAMASK
!--------------------------------------------------------------------------------------------------
!> @brief reads float x.y+z value at myChunk from format string
!--------------------------------------------------------------------------------------------------
real(pReal) function IO_fixedNoEFloatValue (string,ends,myChunk)

  character(len=*),      intent(in) :: string                                                       !< raw input with known ends of each chunk
  integer,               intent(in) :: myChunk                                                      !< position number of desired chunk
  integer, dimension(:), intent(in) :: ends                                                         !< positions of end of each tag/chunk in given string
  character(len=22),     parameter  :: MYNAME = 'IO_fixedNoEFloatValue '
  character(len=13),     parameter  :: VALIDBASE = '0123456789.+-'
  character(len=12),     parameter  :: VALIDEXP  = '0123456789+-'
 
  real(pReal)   :: base
  integer :: expon
  integer       :: pos_exp
 
  pos_exp = scan(string(ends(myChunk)+1:ends(myChunk+1)),'+-',back=.true.)
  hasExponent: if (pos_exp > 1) then
    base  = IO_verifyFloatValue(trim(adjustl(string(ends(myChunk)+1:ends(myChunk)+pos_exp-1))),&
                                VALIDBASE,MYNAME//'(base): ')
    expon = IO_verifyIntValue(trim(adjustl(string(ends(myChunk)+pos_exp:ends(myChunk+1)))),&
                                VALIDEXP,MYNAME//'(exp): ')
  else hasExponent
    base  = IO_verifyFloatValue(trim(adjustl(string(ends(myChunk)+1:ends(myChunk+1)))),&
                                VALIDBASE,MYNAME//'(base): ')
    expon = 0
  endif hasExponent
  IO_fixedNoEFloatValue = base*10.0_pReal**real(expon,pReal)

end function IO_fixedNoEFloatValue


!--------------------------------------------------------------------------------------------------
!> @brief reads integer value at myChunk from fixed format string
!--------------------------------------------------------------------------------------------------
integer function IO_fixedIntValue(string,ends,myChunk)

  character(len=*),      intent(in) :: string                                                       !< raw input with known ends of each chunk
  integer,               intent(in) :: myChunk                                                      !< position number of desired chunk
  integer, dimension(:), intent(in) :: ends                                                         !< positions of end of each tag/chunk in given string
  character(len=20),     parameter  :: MYNAME = 'IO_fixedIntValue: '
  character(len=12),     parameter  :: VALIDCHARACTERS = '0123456789+-'
 
  IO_fixedIntValue = IO_verifyIntValue(trim(adjustl(string(ends(myChunk)+1:ends(myChunk+1)))),&
                                       VALIDCHARACTERS,MYNAME)

end function IO_fixedIntValue
#endif


!--------------------------------------------------------------------------------------------------
!> @brief changes characters in string to lower case
!--------------------------------------------------------------------------------------------------
pure function IO_lc(string)

  character(len=*), intent(in) :: string                                                            !< string to convert
  character(len=len(string))   :: IO_lc

  character(26), parameter :: LOWER = 'abcdefghijklmnopqrstuvwxyz'
  character(26), parameter :: UPPER = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  integer                      :: i,n

  IO_lc = string
  do i=1,len(string)
    n = index(UPPER,IO_lc(i:i))
    if (n/=0) IO_lc(i:i) = LOWER(n:n)
  enddo

end function IO_lc


!--------------------------------------------------------------------------------------------------
!> @brief returns format string for integer values without leading zeros
!--------------------------------------------------------------------------------------------------
pure function IO_intOut(intToPrint)

  integer, intent(in) :: intToPrint
  character(len=41)   :: IO_intOut
  integer             :: N_digits
  character(len=19)   :: width                                                                      ! maximum digits for 64 bit integer
  character(len=20)   :: min_width                                                                  ! longer for negative values

  N_digits =  1 + int(log10(real(max(abs(intToPrint),1))))
  write(width, '(I19.19)') N_digits
  write(min_width, '(I20.20)') N_digits + merge(1,0,intToPrint < 0)
  IO_intOut = 'I'//trim(min_width)//'.'//trim(width)

end function IO_intOut


!--------------------------------------------------------------------------------------------------
!> @brief write error statements to standard out and terminate the Marc/spectral run with exit #9xxx
!> in ABAQUS either time step is reduced or execution terminated
!--------------------------------------------------------------------------------------------------
subroutine IO_error(error_ID,el,ip,g,instance,ext_msg)

  integer,                    intent(in) :: error_ID
  integer,          optional, intent(in) :: el,ip,g,instance
  character(len=*), optional, intent(in) :: ext_msg

  external                               :: quit
  character(len=1024)                    :: msg
  character(len=1024)                    :: formatString

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
    case (104)
      msg = '{input} recursion limit reached'
    case (105)
      msg = 'unknown output:'
    case (106)
      msg = 'working directory does not exist:'
    case (107)
      msg = 'line length exceeds limit of 256'

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
! DAMASK_marc errors
    case (700)
      msg = 'invalid materialpoint result requested'

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
! error messages related to parsing of Abaqus input file
    case (900)
      msg = 'improper definition of nodes in input file (Nnodes < 2)'
    case (901)
      msg = 'no elements defined in input file (Nelems = 0)'
    case (902)
      msg = 'no element sets defined in input file (No *Elset exists)'
    case (903)
      msg = 'no materials defined in input file (Look into section assigments)'
    case (904)
      msg = 'no elements could be assigned for Elset: '
    case (905)
      msg = 'error in mesh_abaqus_map_materials'
    case (906)
      msg = 'error in mesh_abaqus_count_cpElements'
    case (907)
      msg = 'size of mesh_mapFEtoCPelem in mesh_abaqus_map_elements'
    case (908)
      msg = 'size of mesh_mapFEtoCPnode in mesh_abaqus_map_nodes'
    case (909)
      msg = 'size of mesh_node in mesh_abaqus_build_nodes not equal to mesh_Nnodes'


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
 
  character(len=1024)                    :: msg
  character(len=1024)                    :: formatString
 
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


#if defined(Abaqus) || defined(Marc4DAMASK)


!--------------------------------------------------------------------------------------------------
!> @brief reads a line from a text file.
!--------------------------------------------------------------------------------------------------
function IO_read(fileUnit) result(line)
 
  integer, intent(in) :: fileUnit                                                                   !< file unit
 
  character(len=pStringLen) :: line
 
 
  read(fileUnit,'(a256)',END=100) line
 
100 end function IO_read


#ifdef Abaqus
!--------------------------------------------------------------------------------------------------
!> @brief extracts string value from key=value pair and check whether key matches
!--------------------------------------------------------------------------------------------------
character(len=300) pure function IO_extractValue(pair,key)

  character(len=*), intent(in) :: pair, &                                                           !< key=value pair
                                  key                                                               !< key to be expected
 
  character(len=*), parameter  :: SEP = achar(61)                                                   ! '='
 
  integer                      :: myChunk                                                           !< position number of desired chunk
 
  IO_extractValue = ''
 
  myChunk = scan(pair,SEP)
  if (myChunk > 0 .and. pair(:myChunk-1) == key) IO_extractValue = pair(myChunk+1:)                 ! extract value if key matches

end function IO_extractValue


!--------------------------------------------------------------------------------------------------
!> @brief count lines containig data up to next *keyword
!--------------------------------------------------------------------------------------------------
integer function IO_countDataLines(fileUnit)

  integer, intent(in)                :: fileUnit                                                    !< file handle
 
 
  integer, allocatable, dimension(:) :: chunkPos
  character(len=65536)               :: line, &
                                        tmp
 
  IO_countDataLines = 0
  line = ''
 
  do while (trim(line) /= IO_EOF)
    line = IO_read(fileUnit)
    chunkPos = IO_stringPos(line)
    tmp = IO_lc(IO_stringValue(line,chunkPos,1))
    if (tmp(1:1) == '*' .and. tmp(2:2) /= '*') then                                                 ! found keyword
      exit
    else
      if (tmp(2:2) /= '*') IO_countDataLines = IO_countDataLines + 1
    endif
  enddo
  backspace(fileUnit)

end function IO_countDataLines
#endif


#ifdef Marc4DAMASK
!--------------------------------------------------------------------------------------------------
!> @brief count lines containig data up to next *keyword
!--------------------------------------------------------------------------------------------------
integer function IO_countNumericalDataLines(fileUnit)

  integer, intent(in)                :: fileUnit                                                    !< file handle
 
 
  integer, allocatable, dimension(:) :: chunkPos
  character(len=65536)               :: line, &
                                        tmp
 
  IO_countNumericalDataLines = 0
  line = ''
 
  do while (trim(line) /= IO_EOF)
    line = IO_read(fileUnit)
    chunkPos = IO_stringPos(line)
    tmp = IO_lc(IO_stringValue(line,chunkPos,1))
    if (verify(trim(tmp),'0123456789') == 0) then                                                   ! numerical values
      IO_countNumericalDataLines = IO_countNumericalDataLines + 1
    else
      exit
    endif
  enddo
  backspace(fileUnit)

end function IO_countNumericalDataLines
#endif


!--------------------------------------------------------------------------------------------------
!> @brief count items in consecutive lines depending on lines
!> @details Marc:      ints concatenated by "c" as last char or range of values a "to" b
!> Abaqus:    triplet of start,stop,inc
!--------------------------------------------------------------------------------------------------
integer function IO_countContinuousIntValues(fileUnit)

 integer, intent(in) :: fileUnit

#ifdef Abaqus
 integer                            :: l,c
#endif
 integer, allocatable, dimension(:) :: chunkPos
 character(len=65536)               :: line

 IO_countContinuousIntValues = 0
 line = ''

#if defined(Marc4DAMASK)
 do while (trim(line) /= IO_EOF)
   line = IO_read(fileUnit)
   chunkPos = IO_stringPos(line)
   if (chunkPos(1) < 1) then                                                                        ! empty line
     exit
   elseif (IO_lc(IO_stringValue(line,chunkPos,2)) == 'to' ) then                                    ! found range indicator
     IO_countContinuousIntValues = 1 + abs( IO_intValue(line,chunkPos,3) &
                                           -IO_intValue(line,chunkPos,1))
     exit                                                                                           ! only one single range indicator allowed                              
   else
     IO_countContinuousIntValues = IO_countContinuousIntValues+chunkPos(1)-1                        ! add line's count when assuming 'c'
     if ( IO_lc(IO_stringValue(line,chunkPos,chunkPos(1))) /= 'c' ) then                            ! line finished, read last value
       IO_countContinuousIntValues = IO_countContinuousIntValues+1
       exit                                                                                         ! data ended
     endif
   endif
 enddo
#elif defined(Abaqus)
 c = IO_countDataLines(fileUnit)
 do l = 1,c
   backspace(fileUnit)
 enddo

 l = 1
 do while (trim(line) /= IO_EOF .and. l <= c)                                                       ! ToDo: is this correct?
   l = l + 1
   line = IO_read(fileUnit)
   chunkPos = IO_stringPos(line)
   IO_countContinuousIntValues = IO_countContinuousIntValues + 1 + &                                ! assuming range generation
                            (IO_intValue(line,chunkPos,2)-IO_intValue(line,chunkPos,1))/&
                                                     max(1,IO_intValue(line,chunkPos,3))
 enddo
#endif

end function IO_countContinuousIntValues


!--------------------------------------------------------------------------------------------------
!> @brief return integer list corresponding to items in consecutive lines.
!! First integer in array is counter
!> @details Marc:      ints concatenated by "c" as last char, range of a "to" b, or named set
!! Abaqus:    triplet of start,stop,inc or named set
!--------------------------------------------------------------------------------------------------
function IO_continuousIntValues(fileUnit,maxN,lookupName,lookupMap,lookupMaxN)

 integer,                           intent(in) :: maxN
 integer,           dimension(1+maxN)          :: IO_continuousIntValues

 integer,                           intent(in) :: fileUnit, &
                                                  lookupMaxN
 integer,           dimension(:,:), intent(in) :: lookupMap
 character(len=64), dimension(:),   intent(in) :: lookupName
 integer :: i,first,last
#ifdef Abaqus
 integer :: j,l,c
#endif
 integer, allocatable, dimension(:) :: chunkPos
 character(len=65536) line
 logical rangeGeneration

 IO_continuousIntValues = 0
 rangeGeneration = .false.

#if defined(Marc4DAMASK)
 do
   read(fileUnit,'(A65536)',end=100) line
   chunkPos = IO_stringPos(line)
   if (chunkPos(1) < 1) then                                                                        ! empty line
     exit
   elseif (verify(IO_stringValue(line,chunkPos,1),'0123456789') > 0) then                           ! a non-int, i.e. set name
     do i = 1, lookupMaxN                                                                           ! loop over known set names
       if (IO_stringValue(line,chunkPos,1) == lookupName(i)) then                                   ! found matching name
         IO_continuousIntValues = lookupMap(:,i)                                                    ! return resp. entity list
         exit
       endif
     enddo
     exit
   else if (chunkPos(1) > 2 .and. IO_lc(IO_stringValue(line,chunkPos,2)) == 'to' ) then             ! found range indicator
     first = IO_intValue(line,chunkPos,1)
     last  = IO_intValue(line,chunkPos,3)
     do i = first, last, sign(1,last-first)
       IO_continuousIntValues(1) = IO_continuousIntValues(1) + 1
       IO_continuousIntValues(1+IO_continuousIntValues(1)) = i
     enddo
     exit
   else
     do i = 1,chunkPos(1)-1                                                                         ! interpret up to second to last value
       IO_continuousIntValues(1) = IO_continuousIntValues(1) + 1
       IO_continuousIntValues(1+IO_continuousIntValues(1)) = IO_intValue(line,chunkPos,i)
     enddo
     if ( IO_lc(IO_stringValue(line,chunkPos,chunkPos(1))) /= 'c' ) then                            ! line finished, read last value
       IO_continuousIntValues(1) = IO_continuousIntValues(1) + 1
       IO_continuousIntValues(1+IO_continuousIntValues(1)) = IO_intValue(line,chunkPos,chunkPos(1))
       exit
     endif
   endif
 enddo
#elif defined(Abaqus)
 c = IO_countDataLines(fileUnit)
 do l = 1,c
   backspace(fileUnit)
 enddo

!--------------------------------------------------------------------------------------------------
! check if the element values in the elset are auto generated
 backspace(fileUnit)
 read(fileUnit,'(A65536)',end=100) line
 chunkPos = IO_stringPos(line)
 do i = 1,chunkPos(1)
   if (IO_lc(IO_stringValue(line,chunkPos,i)) == 'generate') rangeGeneration = .true.
 enddo

 do l = 1,c
   read(fileUnit,'(A65536)',end=100) line
   chunkPos = IO_stringPos(line)
   if (verify(IO_stringValue(line,chunkPos,1),'0123456789') > 0) then                               ! a non-int, i.e. set names follow on this line
     do i = 1,chunkPos(1)                                                                           ! loop over set names in line
       do j = 1,lookupMaxN                                                                          ! look through known set names
         if (IO_stringValue(line,chunkPos,i) == lookupName(j)) then                                 ! found matching name
           first = 2 + IO_continuousIntValues(1)                                                    ! where to start appending data
           last  = first + lookupMap(1,j) - 1                                                       ! up to where to append data
           IO_continuousIntValues(first:last) = lookupMap(2:1+lookupMap(1,j),j)                     ! add resp. entity list
           IO_continuousIntValues(1) = IO_continuousIntValues(1) + lookupMap(1,j)                   ! count them
         endif
       enddo
     enddo
   else if (rangeGeneration) then                                                                   ! range generation
     do i = IO_intValue(line,chunkPos,1),&
            IO_intValue(line,chunkPos,2),&
            max(1,IO_intValue(line,chunkPos,3))
       IO_continuousIntValues(1) = IO_continuousIntValues(1) + 1
       IO_continuousIntValues(1+IO_continuousIntValues(1)) = i
     enddo
   else                                                                                             ! read individual elem nums
     do i = 1,chunkPos(1)
       IO_continuousIntValues(1) = IO_continuousIntValues(1) + 1
       IO_continuousIntValues(1+IO_continuousIntValues(1)) = IO_intValue(line,chunkPos,i)
     enddo
   endif
 enddo
#endif

100 end function IO_continuousIntValues
#endif

!--------------------------------------------------------------------------------------------------
! internal helper functions

!--------------------------------------------------------------------------------------------------
!> @brief returns verified integer value in given string
!--------------------------------------------------------------------------------------------------
integer function IO_verifyIntValue (string,validChars,myName)
 
  character(len=*), intent(in) :: string, &                                                         !< string for conversion to int value. Must not contain spaces!
                                  validChars, &                                                     !< valid characters in string
                                  myName                                                            !< name of caller function (for debugging)
  integer                      :: readStatus, invalidWhere
 
  IO_verifyIntValue = 0
 
  invalidWhere = verify(string,validChars)
  if (invalidWhere == 0) then
    read(UNIT=string,iostat=readStatus,FMT=*) IO_verifyIntValue                                     ! no offending chars found
    if (readStatus /= 0) &                                                                          ! error during string to integer conversion
      call IO_warning(203,ext_msg=myName//'"'//string//'"')
  else
    call IO_warning(202,ext_msg=myName//'"'//string//'"')                                           ! complain about offending characters
    read(UNIT=string(1:invalidWhere-1),iostat=readStatus,FMT=*) IO_verifyIntValue                   ! interpret remaining string
    if (readStatus /= 0) &                                                                          ! error during string to integer conversion
      call IO_warning(203,ext_msg=myName//'"'//string(1:invalidWhere-1)//'"')
  endif
 
end function IO_verifyIntValue


!--------------------------------------------------------------------------------------------------
!> @brief returns verified float value in given string
!--------------------------------------------------------------------------------------------------
real(pReal) function IO_verifyFloatValue (string,validChars,myName)
 
  character(len=*), intent(in) :: string, &                                                         !< string for conversion to int value. Must not contain spaces!
                                  validChars, &                                                     !< valid characters in string
                                  myName                                                            !< name of caller function (for debugging)
 
  integer                      :: readStatus, invalidWhere
 
  IO_verifyFloatValue = 0.0_pReal
 
  invalidWhere = verify(string,validChars)
  if (invalidWhere == 0) then
    read(UNIT=string,iostat=readStatus,FMT=*) IO_verifyFloatValue                                   ! no offending chars found
    if (readStatus /= 0) &                                                                          ! error during string to float conversion
      call IO_warning(203,ext_msg=myName//'"'//string//'"')
  else
    call IO_warning(202,ext_msg=myName//'"'//string//'"')                                           ! complain about offending characters
    read(UNIT=string(1:invalidWhere-1),iostat=readStatus,FMT=*) IO_verifyFloatValue                 ! interpret remaining string
    if (readStatus /= 0) &                                                                          ! error during string to float conversion
      call IO_warning(203,ext_msg=myName//'"'//string(1:invalidWhere-1)//'"')
  endif
 
end function IO_verifyFloatValue

end module IO
