!--------------------------------------------------------------------------------------------------
!> @author   Martin Dieh, Max-Planck-Institut für Eisenforschung GmbH
!> @brief    Chained list to store string together with position of delimiters
!--------------------------------------------------------------------------------------------------
module linked_list
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 type, private :: tPartitionedString
   character(len=:),            allocatable :: val
   integer(pInt), dimension(:), allocatable :: pos
 end type tPartitionedString
 
 type, public :: tPartitionedStringList
   type(tPartitionedString)               :: string
   type(tPartitionedStringList),  pointer :: next => null()
   type(tPartitionedStringList),  pointer :: prev => null()
   contains
     procedure :: add            => add
     procedure :: show           => show

     procedure :: keyExists      => exist
     procedure :: countKeys      => count
     procedure :: getStringsRaw  => strings

     procedure :: getRaw         => getRaw
     procedure :: getRaws        => getRaws

     procedure :: getFloat       => getFloat
     procedure :: getFloatArray  => getFloatArray

     procedure :: getInt         => getInt
     procedure :: getIntArray    => getIntArray

     procedure :: getString      => getString
     procedure :: getStrings     => getStrings

 end type tPartitionedStringList

 type(tPartitionedStringList), public :: emptyList

contains

!--------------------------------------------------------------------------------------------------
!> @brief add element
!> @details Adds a string together with the start/end position of chunks in this string. The new 
!! element is added at the end of the list. Empty strings are not added. All strings are converted
!! to lower case
!--------------------------------------------------------------------------------------------------
subroutine add(this,string)
  use IO, only: &
    IO_isBlank, &
    IO_lc, &
    IO_stringPos

  implicit none
  class(tPartitionedStringList),  target, intent(in) :: this
  character(len=*),                       intent(in) :: string
  type(tPartitionedStringList),   pointer            :: new, item

  if (IO_isBlank(string)) return

  allocate(new)
  new%string%val = IO_lc       (trim(string))
  new%string%pos = IO_stringPos(trim(string))

  item => this
  do while (associated(item%next))
    item => item%next
  enddo
  item%next => new

end subroutine add


!--------------------------------------------------------------------------------------------------
!> @brief prints all elements
!> @details Strings are printed in order of insertion (FIFO)
!--------------------------------------------------------------------------------------------------
subroutine show(this)

 implicit none
 class(tPartitionedStringList) :: this
 type(tPartitionedStringList), pointer :: item

 item => this%next
 do while (associated(item))
   write(6,'(a)') trim(item%string%val)
   item => item%next
 end do

end subroutine show


!--------------------------------------------------------------------------------------------------
!> @brief deallocates all elements of a given list
!> @details Strings are printed in order of insertion (FIFO)
!--------------------------------------------------------------------------------------------------
!    subroutine free_all()
!      implicit none
!                 
!      type(node), pointer :: item
!         
!      do        
!        item => first
!         
!        if (associated(item) .eqv. .FALSE.) exit
!          
!        first => first%next
!        deallocate(item)
!      end do                     
!    end subroutine free_all


!--------------------------------------------------------------------------------------------------
!> @brief reports wether a given key (string value at first position) exists in the list
!--------------------------------------------------------------------------------------------------
logical function exist(this,key)
 use IO, only: &
   IO_stringValue

 implicit none
 class(tPartitionedStringList), intent(in) :: this
 character(len=*), intent(in)              :: key
 type(tPartitionedStringList), pointer     :: item

 exist = .false.

 item => this%next
 do while (associated(item) .and. .not. exist)
   exist = trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)
   item => item%next
 end do

end function exist


!--------------------------------------------------------------------------------------------------
!> @brief count number of key appearances
!> @details traverses list and counts each occurrence of specified key
!--------------------------------------------------------------------------------------------------
integer(pInt) function count(this,key)
 use IO, only: &
   IO_stringValue

 implicit none

 class(tPartitionedStringList), intent(in) :: this
 character(len=*), intent(in)              :: key
 type(tPartitionedStringList), pointer     :: item
 integer(pInt) :: i

 count = 0_pInt

 item => this%next
 do while (associated(item))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) &
     count = count + 1_pInt
   item => item%next
 end do

end function count


!--------------------------------------------------------------------------------------------------
!> @brief returns all strings in the list
!> @details returns raw string without start/end position of chunks
!--------------------------------------------------------------------------------------------------
function strings(this)
 use IO, only: &
   IO_error, &
   IO_stringValue

 implicit none
 class(tPartitionedStringList),      intent(in)  :: this
 character(len=65536), dimension(:), allocatable :: strings
 character(len=65536)                            :: string
 type(tPartitionedStringList),  pointer          :: item

 item => this%next
 do while (associated(item))
   string = item%string%val
   GfortranBug86033: if (.not. allocated(strings)) then
     allocate(strings(1),source=string)
   else GfortranBug86033
     strings = [strings,string]
   endif GfortranBug86033
   item => item%next
 end do

 if (size(strings) < 0_pInt) call IO_error(142_pInt)                           ! better to check for "allocated"?

end function strings


!--------------------------------------------------------------------------------------------------
!> @brief gets first string that matches given key (i.e. first chunk)
!> @details returns raw string and start/end position of chunks in this string
!--------------------------------------------------------------------------------------------------
subroutine getRaw(this,key,string,stringPos)
 use IO, only : &
   IO_error, &
   IO_stringValue

 implicit none
 class(tPartitionedStringList),            intent(in)  :: this
 character(len=*),                         intent(in)  :: key
 character(len=*),                         intent(out) :: string
 integer(pInt), dimension(:), allocatable, intent(out) :: stringPos
 type(tPartitionedStringList), pointer                 :: item
 logical                                               :: found
 
 found = .false.

 item => this%next
 do while (associated(item) .and. .not. found)
   found = trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)
   if (found) then
     stringPos = item%string%pos
     string    = item%string%val
   endif
   item => item%next
 end do

 if (.not. found) call IO_error(140_pInt,ext_msg=key)

end subroutine getRaw


!--------------------------------------------------------------------------------------------------
!> @brief gets all strings that matches given key (i.e. first chunk)
!> @details returns raw strings and start/end positions of chunks in these strings.
! Will fail if number of positions in strings differs.
!--------------------------------------------------------------------------------------------------
subroutine getRaws(this,key,string,stringPos)
 use IO, only: &
   IO_error, &
   IO_stringValue

 implicit none
 class(tPartitionedStringList),                     intent(in)  :: this
 character(len=*),                                  intent(in)  :: key
 character(len=65536), dimension(:),   allocatable, intent(out) :: string
 integer(pInt),        dimension(:,:), allocatable, intent(out) :: stringPos
 
 character(len=65536)                      :: string_tmp
 integer(pInt)                             :: posSize
 integer(pInt),  dimension(:), allocatable :: stringPosFlat
 type(tPartitionedStringList), pointer     :: item

 posSize = -1_pInt
 item => this%next
 do 
   if (.not. associated(item)) then
     if (posSize < 0_pInt) call IO_error(140_pInt,ext_msg=key)
     stringPos = reshape(stringPosFlat,[posSize,size(string)])
     exit
   endif
   foundKey: if (trim(IO_stringValue(item%string%val,item%string%pos,1))==trim(key)) then
     if (posSize < 0_pInt) then
       posSize = size(item%string%pos)
       stringPosFlat = item%string%pos
       allocate(string(1))
       string(1) = item%string%val
     else
       if (size(item%string%pos) /= posSize) &
         call IO_error(141_pInt,ext_msg=trim(item%string%val),el=posSize)
       stringPosFlat = [stringPosFlat,item%string%pos]
       string_tmp = item%string%val
       string = [string,string_tmp]
     endif 
   endif foundKey
   item => item%next
 end do

end subroutine getRaws


!--------------------------------------------------------------------------------------------------
!> @brief gets float value of first string that matches given key (i.e. first chunk)
!> @details gets one float value. If key is not found exits with error unless default is given
!--------------------------------------------------------------------------------------------------
real(pReal) function getFloat(this,key,defaultVal)
 use IO, only : &
   IO_error, &
   IO_stringValue, &
   IO_FloatValue

 implicit none
 class(tPartitionedStringList), intent(in)           :: this
 character(len=*),              intent(in)           :: key
 real(pReal),                   intent(in), optional :: defaultVal
 type(tPartitionedStringList),  pointer              :: item
 logical                                             :: found

 if (present(defaultVal)) getFloat = defaultVal
 found = present(defaultVal)
 
 item => this%next
 do while (associated(item))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) then
     found = .true.
     if (item%string%pos(1) < 2_pInt) call IO_error(143_pInt,ext_msg=key)
     getFloat = IO_FloatValue(item%string%val,item%string%pos,2)
   endif
   item => item%next
 end do

 if (.not. found) call IO_error(140_pInt,ext_msg=key)

end function getFloat


!--------------------------------------------------------------------------------------------------
!> @brief gets integer value for given key
!> @details gets one integer value. If key is not found exits with error unless default is given
!--------------------------------------------------------------------------------------------------
integer(pInt) function getInt(this,key,defaultVal)
 use IO, only: &
   IO_error, &
   IO_stringValue, &
   IO_IntValue

 implicit none
 class(tPartitionedStringList), intent(in)           :: this
 character(len=*),              intent(in)           :: key
 integer(pInt),                 intent(in), optional :: defaultVal
 type(tPartitionedStringList),  pointer              :: item
 logical                                             :: found

 if (present(defaultVal)) getInt = defaultVal
 found = present(defaultVal)
 
 item => this%next
 do while (associated(item))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) then
     found = .true.
     if (item%string%pos(1) < 2_pInt) call IO_error(143_pInt,ext_msg=key)
     getInt = IO_IntValue(item%string%val,item%string%pos,2)
   endif
   item => item%next
 end do

 if (.not. found) call IO_error(140_pInt,ext_msg=key)

end function getInt


!--------------------------------------------------------------------------------------------------
!> @brief gets string value for given key
!> @details if key is not found exits with error unless default is given
!--------------------------------------------------------------------------------------------------
character(len=65536) function getString(this,key,defaultVal,raw)
 use IO, only: &
   IO_error, &
   IO_stringValue

 implicit none
 class(tPartitionedStringList), intent(in)           :: this
 character(len=*),              intent(in)           :: key
 character(len=65536),          intent(in), optional :: defaultVal
 logical,                       intent(in), optional :: raw
 type(tPartitionedStringList),  pointer              :: item
 logical                                             :: found, &
                                                        split

 if (present(defaultVal)) getString = defaultVal
 split = merge(raw,.true.,present(raw))
 found = present(defaultVal)

 item => this%next
 do while (associated(item))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) then
     found = .true.
     if (split) then
       if (item%string%pos(1) < 2_pInt) call IO_error(143_pInt,ext_msg=key)
       getString = IO_StringValue(item%string%val,item%string%pos,2)
     else
       getString = trim(item%string%val(item%string%pos(4):))                                  ! raw string starting a second chunk
     endif
   endif
   item => item%next
 end do

 if (.not. found) call IO_error(140_pInt,ext_msg=key)

end function getString


!--------------------------------------------------------------------------------------------------
!> @brief ...
!> @details ...
!--------------------------------------------------------------------------------------------------
function getStrings(this,key)
  use IO

  implicit none
  character(len=64),dimension(:), allocatable :: getStrings
  class(tPartitionedStringList),   intent(in) :: this
  character(len=*),                intent(in) :: key
  type(tPartitionedStringList), pointer       :: item
  character(len=64)                           :: str
  integer(pInt)                               :: i


  item => this%next
  do 
    if (.not. associated(item)) then
      if (.not. allocated(getStrings)) allocate(getStrings(0),source=str)
      exit
    endif
    if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) then
      if (item%string%pos(1) < 2) print*, "NOT WORKING"
      str = IO_StringValue(item%string%val,item%string%pos,2)

 GfortranBug86033: if (.not. allocated(getStrings)) then
   allocate(getStrings(1),source=str)
 else GfortranBug86033
   getStrings  = [getStrings,str]
 endif GfortranBug86033
    endif
    item => item%next
  end do
end function


!--------------------------------------------------------------------------------------------------
!> @brief gets array of int values for given key
!> @details if key is not found exits with error unless default is given
!--------------------------------------------------------------------------------------------------
function getIntArray(this,key,defaultVal)
 use IO, only: &
   IO_error, &
   IO_stringValue, &
   IO_IntValue

 implicit none
 integer(pInt), dimension(:), allocatable            :: getIntArray
 class(tPartitionedStringList), intent(in)           :: this
 character(len=*),              intent(in)           :: key
 integer(pInt), dimension(:),   intent(in), optional :: defaultVal
 type(tPartitionedStringList),  pointer              :: item
 integer(pInt)                                       :: i
 logical                                             :: found, &
                                                        cumulative

 cumulative = (key(1:1) == '(' .and. key(len_trim(key):len_trim(key)) == ')')
 found = .false.

 if (present(defaultVal)) then
   getIntArray = defaultVal
 else
   allocate(getIntArray(0))
 endif

 item => this%next
 do while (associated(item) .and. (.not. found .or. cumulative))
   found = trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)
   if (found) then
     if (.not. cumulative) then
       deallocate(getIntArray) ! use here rhs allocation with empty list
       allocate(getIntArray(0))
     endif
     if (item%string%pos(1) < 2_pInt) call IO_error(143_pInt,ext_msg=key)
     do i = 2_pInt, item%string%pos(1)
       getIntArray = [getIntArray,IO_IntValue(item%string%val,item%string%pos,i)]
     enddo
   endif
   item => item%next
 end do

 if (.not. found .and. .not. present(defaultVal)) call IO_error(140_pInt,ext_msg=key)

end function getIntArray



!--------------------------------------------------------------------------------------------------
!> @brief gets array of float values for given key
!> @details if key is not found exits with error unless default is given
!--------------------------------------------------------------------------------------------------
function getFloatArray(this,key,defaultVal)
 use IO, only: &
   IO_error, &
   IO_stringValue, &
   IO_FloatValue

 implicit none
 real(pReal), dimension(:), allocatable              :: getFloatArray
 class(tPartitionedStringList), intent(in)           :: this
 character(len=*),              intent(in)           :: key
 real(pReal), dimension(:),     intent(in), optional :: defaultVal
 type(tPartitionedStringList),  pointer              :: item
 integer(pInt)                                       :: i
 logical                                             :: found

 found = .false.

 if (present(defaultVal)) then
   getFloatArray = defaultVal
 else
   allocate(getFloatArray(0))
 endif

 item => this%next
 do while (associated(item) .and. .not. found)
   found = trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)
   if (found) then
     if (item%string%pos(1) < 2_pInt) call IO_error(143_pInt,ext_msg=key)
     do i = 2_pInt, item%string%pos(1)
       getFloatArray = [getFloatArray,IO_FloatValue(item%string%val,item%string%pos,i)]
     enddo
   endif
   item => item%next
 end do

 if (.not. found .and. .not. present(defaultVal)) call IO_error(140_pInt,ext_msg=key)

end function getFloatArray



end module linked_list
