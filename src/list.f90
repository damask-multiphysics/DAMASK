!--------------------------------------------------------------------------------------------------
!> @author   Martin Dieh, Max-Planck-Institut für Eisenforschung GmbH
!> @brief    Chained list to store string together with position of delimiters
!--------------------------------------------------------------------------------------------------
module chained_list
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 type, private :: tPartitionedString
   character(len=:),      allocatable :: val
   integer(pInt), dimension(:), allocatable :: pos
 end type tPartitionedString
 
 type, public :: tPartitionedStringList
   type(tPartitionedString)    :: string
   type(tPartitionedStringList),  pointer :: next => null()
   type(tPartitionedStringList),  pointer :: prev => null()
   contains
     procedure :: add            => add
     procedure :: show           => show

     procedure :: keyExists      => keyExists
     procedure :: countKeys      => countKeyAppearances
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
  type(tPartitionedStringList),   pointer            :: new, tmp

  if (IO_isBlank(string)) return

  allocate(new)
  new%string%val=IO_lc(trim(string))
  new%string%pos=IO_stringPos(trim(string))

  tmp => this
  do while (associated(tmp%next))
    tmp => tmp%next
  enddo
  tmp%next => new

end subroutine add


!--------------------------------------------------------------------------------------------------
!> @brief prints all elements
!> @details Strings are printed in order of insertion (FIFO)
!--------------------------------------------------------------------------------------------------
subroutine show(this)

 implicit none
 class(tPartitionedStringList) :: this
 type(tPartitionedStringList),  pointer  :: tmp

 tmp => this%next
 do 
   if (.not. associated(tmp)) exit
   write(6,'(a)') trim(tmp%string%val)
   tmp => tmp%next
 end do

end subroutine show


!--------------------------------------------------------------------------------------------------
!> @brief deallocates all elements of a given list
!> @details Strings are printed in order of insertion (FIFO)
!--------------------------------------------------------------------------------------------------
!    subroutine free_all()
!      implicit none
!                 
!      type(node), pointer :: tmp
!         
!      do        
!      tmp => first
!         
!        if (associated(tmp) .eqv. .FALSE.) exit
!          
!        first => first%next
!        deallocate(tmp)
!      end do                     
!    end subroutine free_all


!--------------------------------------------------------------------------------------------------
!> @brief reports wether a given key (string value at first position) exists in the list
!--------------------------------------------------------------------------------------------------
logical function keyExists(this,key)
 use IO, only: &
   IO_stringValue

 implicit none
 class(tPartitionedStringList), intent(in) :: this
 character(len=*), intent(in)              :: key
 type(tPartitionedStringList), pointer     :: tmp

 keyExists = .false.

 tmp => this%next
 do 
   if (.not. associated(tmp)) exit
   if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
     keyExists = .true.
     exit
   endif
   tmp => tmp%next
 end do

end function keyExists


!--------------------------------------------------------------------------------------------------
!> @brief prints all elements
!> @details Strings are printed in order of insertion (FIFO)
!--------------------------------------------------------------------------------------------------
integer(pInt) function countKeyAppearances(this,key)
 use IO, only: &
   IO_stringValue

 implicit none

 class(tPartitionedStringList), intent(in) :: this
 character(len=*), intent(in)              :: key
 type(tPartitionedStringList), pointer     :: tmp
 integer(pInt) :: i

 countKeyAppearances = 0_pInt

 tmp => this%next
 do 
   if (.not. associated(tmp)) exit
   if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
     countKeyAppearances = countKeyAppearances + 1_pInt
   endif
   tmp => tmp%next
 end do

end function countKeyAppearances


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
 character(len=65536)                            :: stringTmp
 type(tPartitionedStringList),  pointer          :: tmp

 tmp => this%next
 do 
   if (.not. associated(tmp)) then
     if(size(strings) < 0_pInt) call IO_error(142_pInt)
     exit
   endif
     stringTmp = tmp%string%val
     GfortranBug86033: if (.not. allocated(strings)) then
       allocate(strings(1),source=stringTmp)
     else GfortranBug86033
       strings = [strings,stringTmp]
     endif GfortranBug86033
   tmp => tmp%next
 end do
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
 integer(pInt), dimension(:), allocatable, intent(out) :: stringPos
 character(len=*),                         intent(out) :: string
 type(tPartitionedStringList),  pointer                :: tmp

 tmp => this%next
 do 
   if (.not. associated(tmp)) call IO_error(140_pInt,ext_msg=key)
   foundKey: if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
     stringPos = tmp%string%pos
     string    = tmp%string%val
     exit
   endif foundKey
   tmp => tmp%next
 end do
end subroutine getRaw


!--------------------------------------------------------------------------------------------------
!> @brief gets all strings that matches given key (i.e. first chunk)
!> @details returns raw strings and start/end positions of chunks in these strings. Will fail if
! number of positions in strings differs
!--------------------------------------------------------------------------------------------------
subroutine getRaws(this,key,string,stringPos)
 use IO, only: &
   IO_error, &
   IO_stringValue

 implicit none
 class(tPartitionedStringList),                     intent(in)  :: this
 character(len=*),                                  intent(in)  :: key
 integer(pInt),        dimension(:,:), allocatable, intent(out) :: stringPos
 character(len=65536), dimension(:),   allocatable, intent(out) :: string
 
 character(len=65536)                      :: stringTmp
 integer(pInt)                             :: posSize
 integer(pInt),  dimension(:), allocatable :: stringPosFlat
 type(tPartitionedStringList), pointer     :: tmp

 posSize = -1_pInt
 tmp => this%next
 do 
   if (.not. associated(tmp)) then
     if(posSize < 0_pInt) call IO_error(140_pInt,ext_msg=key)
     stringPos = reshape(stringPosFlat,[posSize,size(string)])
     exit
   endif
   foundKey: if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
     if (posSize < 0_pInt) then
       posSize = size(tmp%string%pos)
       stringPosFlat = tmp%string%pos
       allocate(string(1))
       string(1) = tmp%string%val
     else
       if (size(tmp%string%pos) /= posSize) &
         call IO_error(141_pInt,ext_msg=trim(tmp%string%val),el=posSize)
       stringPosFlat = [stringPosFlat,tmp%string%pos]
       stringTmp = tmp%string%val
       string = [string,stringTmp]
     endif 
   endif foundKey
   tmp => tmp%next
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
 type(tPartitionedStringList),  pointer              :: tmp

 tmp => this%next
 do 
   endOfList: if (.not. associated(tmp)) then
     if(present(defaultVal)) then
       getFloat = defaultVal
       exit
     else
       call IO_error(140_pInt,ext_msg=key)
     endif
   endif endOfList
   foundKey: if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
     if (tmp%string%pos(1) < 2_pInt) call IO_error(143_pInt,ext_msg=key)
     getFloat = IO_FloatValue(tmp%string%val,tmp%string%pos,2)
     exit
   endif foundKey
   tmp => tmp%next
 end do

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
 type(tPartitionedStringList),  pointer              :: tmp

 tmp => this%next
 do 
   endOfList: if (.not. associated(tmp)) then
     if(present(defaultVal)) then
       getInt = defaultVal
       exit
     else
       call IO_error(140_pInt,ext_msg=key)
     endif
   endif endOfList
   foundKey: if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
     if (tmp%string%pos(1) < 2_pInt) call IO_error(143_pInt,ext_msg=key)
     getInt = IO_IntValue(tmp%string%val,tmp%string%pos,2)
     exit
   endif foundKey
   tmp => tmp%next
 end do

end function getInt


!--------------------------------------------------------------------------------------------------
!> @brief gets string value for given key
!> @details if key is not found exits with error unless default is given
!--------------------------------------------------------------------------------------------------
character(len=65536) function getString(this,key,defaultVal)
 use IO, only: &
   IO_error, &
   IO_stringValue

 implicit none
 class(tPartitionedStringList), intent(in)           :: this
 character(len=*),              intent(in)           :: key
 character(len=65536),          intent(in), optional :: defaultVal
 type(tPartitionedStringList),  pointer              :: tmp

 tmp => this%next
 do 
   endOfList: if (.not. associated(tmp)) then
     if(present(defaultVal)) then
       getString = defaultVal
       exit
     else
       call IO_error(140_pInt,ext_msg=key)
     endif
   endif endOfList
   foundKey: if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
     if (tmp%string%pos(1) < 2_pInt) call IO_error(143_pInt,ext_msg=key)
     getString = IO_StringValue(tmp%string%val,tmp%string%pos,2)
     exit
   endif foundKey
   tmp => tmp%next
 end do

end function getString


function getStrings(this,key)
  use IO

  implicit none
  character(len=64),dimension(:),allocatable :: getStrings
  character(len=64) :: str

  class(tPartitionedStringList), intent(in) :: this
  character(len=*), intent(in) :: key
  type(tPartitionedStringList), pointer :: tmp
  integer(pInt) :: i


  tmp => this%next
  do 
    if (.not. associated(tmp)) then
      if (.not. allocated(getStrings)) allocate(getStrings(0),source=str)
      exit
    endif
    if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
      if (tmp%string%pos(1) < 2) print*, "NOT WORKKING"
      str = IO_StringValue(tmp%string%val,tmp%string%pos,2)

 GfortranBug86033: if (.not. allocated(getStrings)) then
   allocate(getStrings(1),source=str)
 else GfortranBug86033
   getStrings  = [getStrings,str]
 endif GfortranBug86033
    endif
    tmp => tmp%next
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
 integer(pInt),dimension(:),    intent(in), optional :: defaultVal
 type(tPartitionedStringList),  pointer              :: tmp
 integer(pInt) :: i

 allocate(getIntArray(0))

 tmp => this%next
 do 
   endOfList: if (.not. associated(tmp)) then
     if(present(defaultVal)) then
       getIntArray = defaultVal
       exit
     else
       call IO_error(140_pInt,ext_msg=key)
     endif
   endif endOfList
   foundKey: if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
     if (tmp%string%pos(1) < 2_pInt) call IO_error(143_pInt,ext_msg=key)
     do i = 2_pInt, tmp%string%pos(1)
       getIntArray = [getIntArray,IO_IntValue(tmp%string%val,tmp%string%pos,i)]
     enddo
     exit
   endif foundKey
   tmp => tmp%next
 end do
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
 real(pReal),dimension(:),      intent(in), optional :: defaultVal
 type(tPartitionedStringList),  pointer              :: tmp
 integer(pInt) :: i

 allocate(getFloatArray(0))

 tmp => this%next
 do 
   endOfList: if (.not. associated(tmp)) then
     if(present(defaultVal)) then
       getFloatArray = defaultVal
       exit
     else
       call IO_error(140_pInt,ext_msg=key)
     endif
   endif endOfList
   foundKey: if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
     if (tmp%string%pos(1) < 2_pInt) call IO_error(143_pInt,ext_msg=key)
     do i = 2_pInt, tmp%string%pos(1)
       getFloatArray = [getFloatArray,IO_FloatValue(tmp%string%val,tmp%string%pos,i)]
     enddo
     exit
   endif foundKey
   tmp => tmp%next
 end do
end function getFloatArray



end module chained_list
