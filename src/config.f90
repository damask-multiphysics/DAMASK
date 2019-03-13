!-------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Reads in the material configuration from file
!> @details Reads the material configuration file, where solverJobName.materialConfig takes
!! precedence over material.config. Stores the raw strings and the positions of delimiters for the
!! parts 'homogenization', 'crystallite', 'phase', 'texture', and 'microstucture'
!--------------------------------------------------------------------------------------------------
module config
 use prec, only: &
   pReal

 implicit none
 private 
  type, private :: tPartitionedString
   character(len=:),            allocatable :: val
   integer, dimension(:), allocatable :: pos
 end type tPartitionedString
 
 type, private :: tPartitionedStringList
   type(tPartitionedString)               :: string
   type(tPartitionedStringList),  pointer :: next => null()
   contains
     procedure :: add            => add
     procedure :: show           => show
     procedure :: free           => free

! currently, a finalize is needed for all shapes of tPartitionedStringList.
! with Fortran 2015, we can define one recursive elemental function
! https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/543326
     final     :: finalize, &
                  finalizeArray 

     procedure :: keyExists      => keyExists
     procedure :: countKeys      => countKeys

     procedure :: getFloat       => getFloat
     procedure :: getInt         => getInt
     procedure :: getString      => getString

     procedure :: getFloats      => getFloats
     procedure :: getInts        => getInts
     procedure :: getStrings     => getStrings


 end type tPartitionedStringList

 type(tPartitionedStringList), public, protected, allocatable, dimension(:) :: &
   config_phase, &
   config_microstructure, &
   config_homogenization, &
   config_texture, &
   config_crystallite
   
  type(tPartitionedStringList), public, protected :: &
   config_numerics, &
   config_debug
 
 character(len=64), dimension(:), allocatable, public, protected :: &
   phase_name, &                                                                                    !< name of each phase
   homogenization_name, &                                                                           !< name of each homogenization
   crystallite_name, &                                                                              !< name of each crystallite setting
   microstructure_name, &                                                                           !< name of each microstructure
   texture_name                                                                                     !< name of each texture


! ToDo: Remove, use size(config_phase) etc
 integer, public, protected :: &
   material_Nphase, &                                                                               !< number of phases
   material_Nhomogenization                                                                         !< number of homogenizations

 public :: &
   config_init, &
   config_deallocate

contains

!--------------------------------------------------------------------------------------------------
!> @brief reads material.config and stores its content per part
!--------------------------------------------------------------------------------------------------
subroutine config_init
 use prec, only: &
   pStringLen
 use DAMASK_interface, only: &
   getSolverJobName
 use IO, only: &
   IO_read_ASCII, &
   IO_error, &
   IO_lc, &
   IO_getTag
 use debug, only: &
   debug_level, &
   debug_material, &
   debug_levelBasic

 implicit none
 integer :: myDebug,i

 character(len=pStringLen) :: &
  line, &
  part
 character(len=pStringLen), dimension(:), allocatable :: fileContent
 logical :: fileExists

 write(6,'(/,a)') ' <<<+-  config init  -+>>>'

 myDebug = debug_level(debug_material)

 inquire(file=trim(getSolverJobName())//'.materialConfig',exist=fileExists)
 if(fileExists) then
   write(6,'(/,a)') ' reading '//trim(getSolverJobName())//'.materialConfig'; flush(6)
   fileContent = read_materialConfig(trim(getSolverJobName())//'.materialConfig')
 else
   inquire(file='material.config',exist=fileExists)
   if(.not. fileExists) call IO_error(100,ext_msg='material.config')
   write(6,'(/,a)') ' reading material.config'; flush(6)
   fileContent = read_materialConfig('material.config')
 endif

 do i = 1, size(fileContent)
   line = trim(fileContent(i))
   part = IO_lc(IO_getTag(line,'<','>'))
   select case (trim(part))
    
     case (trim('phase'))
       call parse_materialConfig(phase_name,config_phase,line,fileContent(i+1:))
       if (iand(myDebug,debug_levelBasic) /= 0) write(6,'(a)') ' Phase          parsed'; flush(6)
    
     case (trim('microstructure'))
       call parse_materialConfig(microstructure_name,config_microstructure,line,fileContent(i+1:))
       if (iand(myDebug,debug_levelBasic) /= 0) write(6,'(a)') ' Microstructure parsed'; flush(6)
    
     case (trim('crystallite'))
       call parse_materialConfig(crystallite_name,config_crystallite,line,fileContent(i+1:))
       if (iand(myDebug,debug_levelBasic) /= 0) write(6,'(a)') ' Crystallite    parsed'; flush(6)
    
     case (trim('homogenization'))
       call parse_materialConfig(homogenization_name,config_homogenization,line,fileContent(i+1:))
       if (iand(myDebug,debug_levelBasic) /= 0) write(6,'(a)') ' Homogenization parsed'; flush(6)
    
     case (trim('texture'))
       call parse_materialConfig(texture_name,config_texture,line,fileContent(i+1:))
       if (iand(myDebug,debug_levelBasic) /= 0) write(6,'(a)') ' Texture        parsed'; flush(6)

   end select

 enddo
 
 material_Nhomogenization = size(config_homogenization)
 material_Nphase          = size(config_phase)
 
 if (material_Nhomogenization    < 1) call IO_error(160,ext_msg='<homogenization>')
 if (size(config_microstructure) < 1) call IO_error(160,ext_msg='<microstructure>')
 if (size(config_crystallite)    < 1) call IO_error(160,ext_msg='<crystallite>')
 if (material_Nphase             < 1) call IO_error(160,ext_msg='<phase>')
 if (size(config_texture)        < 1) call IO_error(160,ext_msg='<texture>')
 
 
  inquire(file='numerics.config', exist=fileExists)
  if (fileExists) then
    write(6,'(/,a)') ' reading numerics.config'; flush(6)
    fileContent = IO_read_ASCII('numerics.config')
    call parse_debugAndNumericsConfig(config_numerics,fileContent)
  endif
  
  inquire(file='debug.config', exist=fileExists)
  if (fileExists) then
    write(6,'(/,a)') ' reading debug.config'; flush(6)
    fileContent = IO_read_ASCII('debug.config')
    call parse_debugAndNumericsConfig(config_debug,fileContent)
  endif

contains


!--------------------------------------------------------------------------------------------------
!> @brief reads material.config
!!        Recursion is triggered by "{path/to/inputfile}" in a line
!--------------------------------------------------------------------------------------------------
recursive function read_materialConfig(fileName,cnt) result(fileContent)
  use IO, only: &
    IO_warning

  implicit none
  character(len=*),          intent(in)                :: fileName
  integer,                   intent(in), optional      :: cnt                                       !< recursion counter
  character(len=pStringLen), dimension(:), allocatable :: fileContent                               !< file content, separated per lines
  character(len=pStringLen), dimension(:), allocatable :: includedContent
  character(len=pStringLen)                            :: line
  character(len=pStringLen), parameter                 :: dummy = 'https://damask.mpie.de'          !< to fill up remaining array
  character(len=:),                 allocatable :: rawData
  integer ::  &
    fileLength, &
    fileUnit, &
    startPos, endPos, &
    myTotalLines, &                                                                                 !< # lines read from file without include statements
    l,i, &
    myStat
  logical :: warned
  
  if (present(cnt)) then
    if (cnt>10) call IO_error(106,ext_msg=trim(fileName))
  endif

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
! split raw data at end of line and handle includes
  warned = .false.
  startPos = 1
  l = 1
  do while (l <= myTotalLines)
    endPos = merge(startPos + scan(rawData(startPos:),new_line('')) - 2,len(rawData),l /= myTotalLines)
    if (endPos - startPos > pStringLen -1) then
      line = rawData(startPos:startPos+pStringLen-1)
      if (.not. warned) then
        call IO_warning(207,ext_msg=trim(fileName),el=l)
        warned = .true.
      endif
    else
      line = rawData(startPos:endpos)
    endif
    startPos = endPos + 2                                                                           ! jump to next line start

    recursion: if (scan(trim(adjustl(line)),'{') == 1 .and. scan(trim(line),'}') > 2) then
      includedContent = read_materialConfig(trim(line(scan(line,'{')+1:scan(line,'}')-1)), &
                        merge(cnt,1,present(cnt)))                                                  ! to track recursion depth
      fileContent     = [ fileContent(1:l-1), includedContent, [(dummy,i=1,myTotalLines-l)] ]       ! add content and grow array
      myTotalLines    = myTotalLines - 1 + size(includedContent)
      l               = l            - 1 + size(includedContent)
    else recursion
      fileContent(l) = line
      l = l + 1
    endif recursion

  enddo

end function read_materialConfig


!--------------------------------------------------------------------------------------------------
!> @brief parses the material.config file
!--------------------------------------------------------------------------------------------------
subroutine parse_materialConfig(sectionNames,part,line, &
                                fileContent)
 implicit none
 character(len=64),            allocatable, dimension(:), intent(out)   :: sectionNames
 type(tPartitionedStringList), allocatable, dimension(:), intent(inout) :: part
 character(len=pStringLen),                               intent(inout) :: line
 character(len=pStringLen),                 dimension(:), intent(in)    :: fileContent

 integer,                allocatable, dimension(:)                :: partPosition             ! position of [] tags + last line in section
 integer        :: i, j
 logical              :: echo

 echo = .false. 

 if (allocated(part)) call IO_error(161,ext_msg=trim(line))
 allocate(partPosition(0))
 
 do i = 1, size(fileContent)
   line = trim(fileContent(i))
   if (IO_getTag(line,'<','>') /= '') exit
   nextSection: if (IO_getTag(line,'[',']') /= '') then
     partPosition = [partPosition, i]
     cycle
   endif nextSection
   if (size(partPosition) < 1) &
     echo = (trim(IO_getTag(line,'/','/')) == 'echo') .or. echo
 enddo

 allocate(sectionNames(size(partPosition)))
 allocate(part(size(partPosition)))

 partPosition = [partPosition, i]                                                                   ! needed when actually storing content

 do i = 1, size(partPosition) -1
   sectionNames(i) = trim(adjustl(IO_getTag(fileContent(partPosition(i)),'[',']')))
   do j = partPosition(i) + 1,  partPosition(i+1) -1
     call part(i)%add(trim(adjustl(fileContent(j))))
   enddo
   if (echo) then
     write(6,*) 'section',i, '"'//trim(sectionNames(i))//'"'
     call part(i)%show()
   endif
 enddo

end subroutine parse_materialConfig


!--------------------------------------------------------------------------------------------------
!> @brief parses the material.config file
!--------------------------------------------------------------------------------------------------
subroutine parse_debugAndNumericsConfig(config_list, &
                                        fileContent)
 implicit none
 type(tPartitionedStringList),              intent(out) :: config_list
 character(len=pStringLen),   dimension(:), intent(in)  :: fileContent
 integer :: i

 do i = 1, size(fileContent)
   call config_list%add(trim(adjustl(fileContent(i))))
 enddo

end subroutine parse_debugAndNumericsConfig

end subroutine config_init


!--------------------------------------------------------------------------------------------------
!> @brief deallocates the linked lists that store the content of the configuration files
!--------------------------------------------------------------------------------------------------
subroutine config_deallocate(what)
 use IO, only: &
   IO_error

 implicit none
 character(len=*), intent(in) :: what

 select case(trim(what))

   case('material.config/phase')
     deallocate(config_phase)

   case('material.config/microstructure')
     deallocate(config_microstructure)

   case('material.config/crystallite')
     deallocate(config_crystallite)

   case('material.config/homogenization')
     deallocate(config_homogenization)

   case('material.config/texture')
     deallocate(config_texture)
     
   case('debug.config')
     call config_debug%free
     
   case('numerics.config')
     call config_numerics%free
     
   case default
     call IO_error(0,ext_msg='config_deallocate')

 end select

end subroutine config_deallocate


!##################################################################################################
! The folowing functions are part of the tPartitionedStringList object
!##################################################################################################



!--------------------------------------------------------------------------------------------------
!> @brief add element
!> @details Adds a string together with the start/end position of chunks in this string. The new 
!! element is added at the end of the list. Empty strings are not added. All strings are converted
!! to lower case. The data is not stored in the new element but in the current.
!--------------------------------------------------------------------------------------------------
subroutine add(this,string)
  use IO, only: &
    IO_isBlank, &
    IO_lc, &
    IO_stringPos

  implicit none
  class(tPartitionedStringList),  target, intent(in) :: this
  character(len=*),                       intent(in) :: string
  type(tPartitionedStringList),   pointer            :: new, temp

  if (IO_isBlank(string)) return

  allocate(new)
  temp => this
  do while (associated(temp%next))
    temp => temp%next
  enddo
  temp%string%val = IO_lc       (trim(string))
  temp%string%pos = IO_stringPos(trim(string))
  temp%next => new

end subroutine add


!--------------------------------------------------------------------------------------------------
!> @brief prints all elements
!> @details Strings are printed in order of insertion (FIFO)
!--------------------------------------------------------------------------------------------------
subroutine show(this)

 implicit none
 class(tPartitionedStringList), target, intent(in) :: this
 type(tPartitionedStringList),  pointer            :: item

 item => this
 do while (associated(item%next))
   write(6,'(a)') ' '//trim(item%string%val)
   item => item%next
 enddo

end subroutine show


!--------------------------------------------------------------------------------------------------
!> @brief empties list and frees associated memory
!> @details explicit interface to reset list. Triggers final statement (and following chain reaction)
!--------------------------------------------------------------------------------------------------
subroutine free(this)

  implicit none
  class(tPartitionedStringList),  intent(inout) :: this

  if(associated(this%next)) deallocate(this%next)

end subroutine free


!--------------------------------------------------------------------------------------------------
!> @brief empties list and frees associated memory
!> @details called when variable goes out of scope. Triggers chain reaction for list
!--------------------------------------------------------------------------------------------------
recursive subroutine finalize(this)

  implicit none
  type(tPartitionedStringList),  intent(inout) :: this

  if(associated(this%next)) deallocate(this%next)

end subroutine finalize


!--------------------------------------------------------------------------------------------------
!> @brief cleans entire array of linke lists
!> @details called when variable goes out of scope and deallocates the list at each array entry
!--------------------------------------------------------------------------------------------------
subroutine finalizeArray(this)

  implicit none
  integer :: i
  type(tPartitionedStringList),  intent(inout), dimension(:) :: this
  type(tPartitionedStringList),  pointer :: temp ! bug in Gfortran?

  do i=1, size(this)
    if (associated(this(i)%next)) then
      temp => this(i)%next
      !deallocate(this(i)) !internal compiler error: in gfc_build_final_call, at fortran/trans.c:975
      deallocate(temp)
    endif
  enddo

end subroutine finalizeArray


!--------------------------------------------------------------------------------------------------
!> @brief reports wether a given key (string value at first position) exists in the list
!--------------------------------------------------------------------------------------------------
logical function keyExists(this,key)
 use IO, only: &
   IO_stringValue

 implicit none
 class(tPartitionedStringList), target, intent(in) :: this
 character(len=*),                      intent(in) :: key
 type(tPartitionedStringList),  pointer            :: item

 keyExists = .false.

 item => this
 do while (associated(item%next) .and. .not. keyExists)
   keyExists = trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)
   item => item%next
 enddo

end function keyExists


!--------------------------------------------------------------------------------------------------
!> @brief count number of key appearances
!> @details traverses list and counts each occurrence of specified key
!--------------------------------------------------------------------------------------------------
integer function countKeys(this,key)
 use IO, only: &
   IO_stringValue

 implicit none

 class(tPartitionedStringList), target, intent(in) :: this
 character(len=*),                      intent(in) :: key
 type(tPartitionedStringList),  pointer            :: item

 countKeys = 0

 item => this
 do while (associated(item%next))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) &
     countKeys = countKeys + 1
   item => item%next
 enddo

end function countKeys


!--------------------------------------------------------------------------------------------------
!> @brief gets float value of for a given key from a linked list
!> @details gets the last value if the key occurs more than once. If key is not found exits with 
!! error unless default is given
!--------------------------------------------------------------------------------------------------
real(pReal) function getFloat(this,key,defaultVal)
 use IO, only : &
   IO_error, &
   IO_stringValue, &
   IO_FloatValue

 implicit none
 class(tPartitionedStringList), target, intent(in)           :: this
 character(len=*),                      intent(in)           :: key
 real(pReal),                           intent(in), optional :: defaultVal
 type(tPartitionedStringList), pointer                       :: item
 logical                                                     :: found

 found = present(defaultVal)
 if (found) getFloat = defaultVal
 
 item => this
 do while (associated(item%next))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) then
     found = .true.
     if (item%string%pos(1) < 2) call IO_error(143,ext_msg=key)
     getFloat = IO_FloatValue(item%string%val,item%string%pos,2)
   endif
   item => item%next
 enddo

 if (.not. found) call IO_error(140,ext_msg=key)

end function getFloat


!--------------------------------------------------------------------------------------------------
!> @brief gets integer value of for a given key from a linked list
!> @details gets the last value if the key occurs more than once. If key is not found exits with 
!! error unless default is given
!--------------------------------------------------------------------------------------------------
integer function getInt(this,key,defaultVal)
 use IO, only: &
   IO_error, &
   IO_stringValue, &
   IO_IntValue

 implicit none
 class(tPartitionedStringList), target, intent(in)           :: this
 character(len=*),                      intent(in)           :: key
 integer,                         intent(in), optional :: defaultVal
 type(tPartitionedStringList), pointer                       :: item
 logical                                                     :: found

 found = present(defaultVal)
 if (found) getInt = defaultVal
 
 item => this
 do while (associated(item%next))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) then
     found = .true.
     if (item%string%pos(1) < 2) call IO_error(143,ext_msg=key)
     getInt = IO_IntValue(item%string%val,item%string%pos,2)
   endif
   item => item%next
 enddo

 if (.not. found) call IO_error(140,ext_msg=key)

end function getInt


!--------------------------------------------------------------------------------------------------
!> @brief gets string value of for a given key from a linked list
!> @details gets the last value if the key occurs more than once. If key is not found exits with 
!! error unless default is given. If raw is true, the the complete string is returned, otherwise 
!! the individual chunks are returned
!--------------------------------------------------------------------------------------------------
character(len=65536) function getString(this,key,defaultVal,raw)
 use IO, only: &
   IO_error, &
   IO_stringValue

 implicit none
 class(tPartitionedStringList), target, intent(in)           :: this
 character(len=*),                      intent(in)           :: key
 character(len=65536),                  intent(in), optional :: defaultVal
 logical,                               intent(in), optional :: raw
 type(tPartitionedStringList),  pointer                      :: item
 logical                                                     :: found, &
                                                                whole
 if (present(raw)) then
   whole = raw
 else
   whole = .false.
 endif

 found = present(defaultVal)
 if (found) then
   getString = trim(defaultVal)
   if (len_trim(getString) /= len_trim(defaultVal)) call IO_error(0,ext_msg='getString')
 endif

 item => this
 do while (associated(item%next))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) then
     found = .true.
     if (item%string%pos(1) < 2) call IO_error(143,ext_msg=key)

     if (whole) then
       getString = trim(item%string%val(item%string%pos(4):))                                       ! raw string starting a second chunk
     else
       getString = IO_StringValue(item%string%val,item%string%pos,2)
     endif
   endif
   item => item%next
 enddo

 if (.not. found) call IO_error(140,ext_msg=key)

end function getString


!--------------------------------------------------------------------------------------------------
!> @brief gets array of float values of for a given key from a linked list
!> @details for cumulative keys, "()", values from all occurrences are return. Otherwise only all
!! values from the last occurrence. If key is not found exits with error unless default is given.
!--------------------------------------------------------------------------------------------------
function getFloats(this,key,defaultVal,requiredSize)
 use IO, only: &
   IO_error, &
   IO_stringValue, &
   IO_FloatValue

 implicit none
 real(pReal),     dimension(:), allocatable          :: getFloats
 class(tPartitionedStringList), target, intent(in)   :: this
 character(len=*),              intent(in)           :: key
 real(pReal),   dimension(:),   intent(in), optional :: defaultVal
 integer,                 intent(in), optional :: requiredSize
 type(tPartitionedStringList),  pointer              :: item
 integer                                       :: i
 logical                                             :: found, &
                                                        cumulative

 cumulative = (key(1:1) == '(' .and. key(len_trim(key):len_trim(key)) == ')')
 found = .false.

 allocate(getFloats(0))

 item => this
 do while (associated(item%next))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) then
     found = .true.
     if (.not. cumulative) getFloats = [real(pReal)::]
     if (item%string%pos(1) < 2) call IO_error(143,ext_msg=key)
     do i = 2, item%string%pos(1)
       getFloats = [getFloats,IO_FloatValue(item%string%val,item%string%pos,i)]
     enddo
   endif
   item => item%next
 enddo

 if (.not. found) then
   if (present(defaultVal)) then; getFloats = defaultVal; else; call IO_error(140,ext_msg=key); endif
 endif
 if (present(requiredSize)) then
   if(requiredSize /= size(getFloats)) call IO_error(146,ext_msg=key)
 endif

end function getFloats


!--------------------------------------------------------------------------------------------------
!> @brief gets array of integer values of for a given key from a linked list
!> @details for cumulative keys, "()", values from all occurrences are return. Otherwise only all
!! values from the last occurrence. If key is not found exits with error unless default is given.
!--------------------------------------------------------------------------------------------------
function getInts(this,key,defaultVal,requiredSize)
 use IO, only: &
   IO_error, &
   IO_stringValue, &
   IO_IntValue

 implicit none
 integer, dimension(:), allocatable            :: getInts
 class(tPartitionedStringList), target, intent(in)   :: this
 character(len=*),              intent(in)           :: key
 integer, dimension(:),   intent(in), optional :: defaultVal
 integer,                 intent(in), optional :: requiredSize
 type(tPartitionedStringList),  pointer              :: item
 integer                                       :: i
 logical                                             :: found, &
                                                        cumulative

 cumulative = (key(1:1) == '(' .and. key(len_trim(key):len_trim(key)) == ')')
 found = .false.

 allocate(getInts(0))

 item => this
 do while (associated(item%next))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) then
     found = .true.
     if (.not. cumulative) getInts = [integer::]
     if (item%string%pos(1) < 2) call IO_error(143,ext_msg=key)
     do i = 2, item%string%pos(1)
       getInts = [getInts,IO_IntValue(item%string%val,item%string%pos,i)]
     enddo
   endif
   item => item%next
 enddo

 if (.not. found) then
   if (present(defaultVal)) then; getInts = defaultVal; else; call IO_error(140,ext_msg=key); endif
 endif
 if (present(requiredSize)) then
   if(requiredSize /= size(getInts)) call IO_error(146,ext_msg=key)
 endif

end function getInts


!--------------------------------------------------------------------------------------------------
!> @brief gets array of string values of for a given key from a linked list
!> @details for cumulative keys, "()", values from all occurrences are return. Otherwise only all
!! values from the last occurrence. If key is not found exits with error unless default is given.
!! If raw is true, the the complete string is returned, otherwise the individual chunks are returned
!--------------------------------------------------------------------------------------------------
function getStrings(this,key,defaultVal,raw)
 use IO, only: &
   IO_error, &
   IO_StringValue

 implicit none
 character(len=65536),dimension(:), allocatable           :: getStrings
 class(tPartitionedStringList), target, intent(in)        :: this
 character(len=*),                   intent(in)           :: key
 character(len=65536),dimension(:),  intent(in), optional :: defaultVal
 logical,                            intent(in), optional :: raw
 type(tPartitionedStringList), pointer                    :: item
 character(len=65536)                                     :: str
 integer                                            :: i
 logical                                                  :: found, &
                                                             whole, &
                                                             cumulative

 cumulative = (key(1:1) == '(' .and. key(len_trim(key):len_trim(key)) == ')')
 if (present(raw)) then
   whole = raw
 else
   whole = .false.
 endif
 found = .false.

 item => this
 do while (associated(item%next))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) then
     found = .true.
     if (allocated(getStrings) .and. .not. cumulative) deallocate(getStrings)
     if (item%string%pos(1) < 2) call IO_error(143,ext_msg=key)
     
     notAllocated: if (.not. allocated(getStrings)) then
       if (whole) then
         str = item%string%val(item%string%pos(4):)
         getStrings = [str]
       else
         str = IO_StringValue(item%string%val,item%string%pos,2)
         allocate(getStrings(1),source=str)
         do i=3,item%string%pos(1)
           str = IO_StringValue(item%string%val,item%string%pos,i)
           getStrings = [getStrings,str]
         enddo
       endif
     else notAllocated
       if (whole) then
         str = item%string%val(item%string%pos(4):)
         getStrings = [getStrings,str]
       else
         do i=2,item%string%pos(1)
           str = IO_StringValue(item%string%val,item%string%pos,i)
           getStrings = [getStrings,str]
         enddo
       endif
     endif notAllocated
   endif
   item => item%next
 enddo

 if (.not. found) then
   if (present(defaultVal)) then; getStrings = defaultVal; else; call IO_error(140,ext_msg=key); endif
 endif

end function getStrings


end module config
