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
   contains
     procedure :: add            => add
     procedure :: getRaw         => getRaw
     procedure :: getRaws         => getRaws

     procedure :: getFloat       => getFloat
     procedure :: getFloatArray  => getFloatArray

     procedure :: getInt         => getInt
     procedure :: getIntArray    => getIntArray

     procedure :: getStrings     => getStrings
     procedure :: keyExists      => keyExists

 end type tPartitionedStringList
 
contains

!--------------------------------------------------------------------------------------------------
!> @brief add element
!> @details adds raw string and start/end position of chunks in this string
!--------------------------------------------------------------------------------------------------
subroutine add(this,string,stringPos)
  implicit none
  class(tPartitionedStringList) :: this
  type(tPartitionedStringList), pointer :: &
    new, &
    tmp
  character(len=*), intent(in) :: string
  integer(pInt), dimension(:), intent(in) :: stringPos

  allocate(new)
  new%string%val=string
  new%string%pos=stringPos

  if (.not. associated(this%next)) then
    this%next => new
  else
    tmp => this%next
    this%next => new
    this%next%next => tmp
  end if

end subroutine add


!--------------------------------------------------------------------------------------------------
!> @brief gets raw data
!> @details returns raw string and start/end position of chunks in this string
!--------------------------------------------------------------------------------------------------
subroutine getRaw(this,key,string,stringPos)
 use IO, only : &
   IO_error, &
   IO_stringValue

 implicit none
 class(tPartitionedStringList), intent(in)           :: this
 character(len=*),              intent(in)           :: key
 integer(pInt), dimension(:),allocatable,   intent(out)          :: stringPos
 character(len=*),              intent(out)          :: string
 type(tPartitionedStringList),  pointer              :: tmp

 tmp => this%next
 do 
   if (.not. associated(tmp)) call IO_error(1_pInt,ext_msg=key)
   foundKey: if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
     stringPos = tmp%string%pos
     string    = tmp%string%val
     exit
   endif foundKey
   tmp => tmp%next
 end do
end subroutine getRaw


!--------------------------------------------------------------------------------------------------
!> @brief gets raw data
!> @details returns raw string and start/end position of chunks in this string
!--------------------------------------------------------------------------------------------------
subroutine getRaws(this,key,string,stringPos)
 use IO, only: &
   IO_error, &
   IO_stringValue

 implicit none
 class(tPartitionedStringList), intent(in)           :: this
 character(len=*),              intent(in)           :: key
 integer(pInt), dimension(:,:),allocatable, intent(out)          :: stringPos
 character(len=256), dimension(:),allocatable,   intent(out)          :: string
 character(len=256)          :: stringTmp
 integer(pInt)  :: posSize
 integer(pInt), dimension(:),allocatable        :: stringPosFlat
 type(tPartitionedStringList),  pointer              :: tmp

 posSize = -1_pInt
 tmp => this%next
 do 
   if (.not. associated(tmp)) then
     if(posSize < 0_pInt) call IO_error(1_pInt,ext_msg=key)
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
       if (size(tmp%string%pos) /= posSize) call IO_error(1_pInt,ext_msg=key)
       stringPosFlat = [stringPosFlat,tmp%string%pos]
       stringTmp = tmp%string%val
       string = [string,stringTmp]
     endif 
   endif foundKey
   tmp => tmp%next
 end do
end subroutine getRaws


!--------------------------------------------------------------------------------------------------
!> @brief gets float value for given key
!> @details if key is not found exits with error unless default is given
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
       call IO_error(1_pInt,ext_msg=key)
     endif
   endif endOfList
   foundKey: if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
     if (tmp%string%pos(1) < 2_pInt) call IO_error(1_pInt,ext_msg=key)
     getFloat = IO_FloatValue(tmp%string%val,tmp%string%pos,2)
     exit
   endif foundKey
   tmp => tmp%next
 end do
end function getFloat


!--------------------------------------------------------------------------------------------------
!> @brief gets float value for given key
!> @details if key is not found exits with error unless default is given
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
       call IO_error(1_pInt,ext_msg=key)
     endif
   endif endOfList
   foundKey: if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
     if (tmp%string%pos(1) < 2_pInt) call IO_error(1_pInt,ext_msg=key)
     getInt = IO_IntValue(tmp%string%val,tmp%string%pos,2)
     exit
   endif foundKey
   tmp => tmp%next
 end do
end function getInt


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
       call IO_error(1_pInt,ext_msg=key)
     endif
   endif endOfList
   foundKey: if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
     if (tmp%string%pos(1) < 2_pInt) call IO_error(1_pInt,ext_msg=key)
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
       call IO_error(1_pInt,ext_msg=key)
     endif
   endif endOfList
   foundKey: if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
     if (tmp%string%pos(1) < 2_pInt) call IO_error(1_pInt,ext_msg=key)
     do i = 2_pInt, tmp%string%pos(1)
       getFloatArray = [getFloatArray,IO_FloatValue(tmp%string%val,tmp%string%pos,i)]
     enddo
     exit
   endif foundKey
   tmp => tmp%next
 end do
end function getFloatArray

! reports wether a key exists at least once
    function keyExists(this,key)
      use IO

      implicit none
      logical :: keyExists

      class(tPartitionedStringList), intent(in) :: this
      character(len=*), intent(in) :: key
      type(tPartitionedStringList), pointer :: tmp

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
    end function


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
        if (.not. associated(tmp)) exit
        if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
          if (tmp%string%pos(1) < 2) print*, "NOT WORKKING"
          str = IO_StringValue(tmp%string%val,tmp%string%pos,2)
          if (.not. allocated(getStrings)) then
            getStrings = [str]
          else
            getStrings = [getStrings,str]
          endif
        endif
        tmp => tmp%next
      end do
    end function

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

end module chained_list
