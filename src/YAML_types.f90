!--------------------------------------------------------------------------------------------------
!> @brief yaml_types
!> @details module describes the various functions to store and get the yaml data.
!! tNode is the fundamental derived data type. It can be of tScalar, &
!! tList or tDict.
!! Every 'value' in a key: value pair is of tNode and is a pointer.
!! If 'value' is of tScalar, it can either be a string, real, integer or logical, &
!! functions exist to convert this scalar type to its respective primitive data type.
!--------------------------------------------------------------------------------------------------

module YAML_types

  use IO
  use prec

  implicit none

  private

  public :: &
    tNode, &
    tScalar, &
    tDict, &
    tList, &
    YAML_types_init

  type, abstract :: tNode
    integer :: length = 0
    contains
    procedure(asFormattedString), deferred :: asFormattedString
    procedure :: &
      asScalar     => tNode_asScalar
    procedure :: &
      asList       => tNode_asList
    procedure :: &
      asDict       => tNode_asDict
   procedure :: &
      tNode_get_byIndex           => tNode_get_byIndex
    procedure :: &
      tNode_get_byIndex_asFloat   => tNode_get_byIndex_asFloat
    procedure :: &
      tNode_get_byIndex_asFloats  => tNode_get_byIndex_asFloats
    procedure :: &
      tNode_get_byIndex_asInt     => tNode_get_byIndex_asInt
    procedure :: &
      tNode_get_byIndex_asInts    => tNode_get_byIndex_asInts
    procedure :: &
      tNode_get_byIndex_asBool    => tNode_get_byIndex_asBool
    procedure :: &
      tNode_get_byIndex_asBools   => tNode_get_byIndex_asBools
    procedure :: &
      tNode_get_byIndex_asString  => tNode_get_byIndex_asString
    procedure :: &
      tNode_get_byIndex_asStrings => tNode_get_byIndex_asStrings
    procedure :: &
      tNode_get_byKey             => tNode_get_byKey
    procedure :: &
      tNode_get_byKey_asFloat     => tNode_get_byKey_asFloat
    procedure :: &
      tNode_get_byKey_asFloats    => tNode_get_byKey_asFloats
    procedure :: &
      tNode_get_byKey_asInt       => tNode_get_byKey_asInt
    procedure :: &
      tNode_get_byKey_asInts      => tNode_get_byKey_asInts
    procedure :: &
      tNode_get_byKey_asBool      => tNode_get_byKey_asBool
    procedure :: &
      tNode_get_byKey_asBools     => tNode_get_byKey_asBools
    procedure :: &
      tNode_get_byKey_asString    => tNode_get_byKey_asString
    procedure :: &
      tNode_get_byKey_asStrings   => tNode_get_byKey_asStrings
    procedure :: &  
      getIndex                    => tNode_get_byKey_asIndex

    generic :: &
      get           => tNode_get_byIndex, &
                       tNode_get_byKey
    generic :: &
      get_asFloat   => tNode_get_byIndex_asFloat, &
                       tNode_get_byKey_asFloat
    generic :: &
      get_asFloats  => tNode_get_byIndex_asFloats, &
                       tNode_get_byKey_asFloats
    generic :: &
      get_asInt     => tNode_get_byIndex_asInt, &
                       tNode_get_byKey_asInt
    generic :: &
      get_asInts    => tNode_get_byIndex_asInts, &
                       tNode_get_byKey_asInts
    generic :: &
      get_asBool    => tNode_get_byIndex_asBool, &
                       tNode_get_byKey_asBool
    generic :: &
      get_asBools   => tNode_get_byIndex_asBools, &
                       tNode_get_byKey_asBools
    generic :: &
      get_asString  => tNode_get_byIndex_asString, &
                       tNode_get_byKey_asString
    generic :: &
      get_asStrings => tNode_get_byIndex_asStrings, &
                       tNode_get_byKey_asStrings
  end type tNode


  type, extends(tNode) :: tScalar

    character(len=:), allocatable, private :: value

    contains
    procedure :: asFormattedString => tScalar_asFormattedString
    procedure :: &
      asFloat   => tScalar_asFloat
    procedure :: &
      asInt     => tScalar_asInt
    procedure :: &
      asBool    => tScalar_asBool
    procedure :: &
      asString  => tScalar_asString
  end type tScalar

  type, extends(tNode) :: tList

    class(tItem), pointer  :: first => null()

    contains
    procedure :: asFormattedString => tList_asFormattedString
    procedure :: append            => tList_append
    procedure :: &
      asFloats  => tList_asFloats
    procedure :: &
      asInts    => tList_asInts
    procedure :: &
      asBools   => tList_asBools
    procedure :: &
      asStrings => tList_asStrings
    final :: tList_finalize
  end type tList

  type, extends(tList) :: tDict
    contains
    procedure :: asFormattedString => tDict_asFormattedString
    procedure :: set               => tDict_set
  end type tDict


  type :: tItem
    character(len=:), allocatable :: key
    class(tNode),     pointer     :: node => null()
    class(tItem),     pointer     :: next => null()
    
    contains
    final :: tItem_finalize
  end type tItem

  abstract interface

    recursive subroutine asFormattedString(self,indent)
      import tNode
      class(tNode), intent(in), target   :: self
      integer,      intent(in), optional :: indent
    end subroutine asFormattedString

  end interface

  interface tScalar
    module procedure tScalar_init__
  end interface tScalar

  interface assignment (=)
    module procedure tScalar_assign__
  end interface assignment (=)

contains

!--------------------------------------------------------------------------------------------------
!> @brief do sanity checks
!--------------------------------------------------------------------------------------------------
subroutine YAML_types_init

  write(6,'(/,a)') ' <<<+-  YAML_types init  -+>>>'

  call unitTest

end subroutine YAML_types_init


!--------------------------------------------------------------------------------------------------
!> @brief check correctness of some type bound procedures
!--------------------------------------------------------------------------------------------------
subroutine unitTest

  class(tNode), pointer  :: s1,s2
  allocate(tScalar::s1)
  allocate(tScalar::s2)
  select type(s1)
    class is(tScalar) 
      s1 = '1'
      if(s1%asInt() /= 1)              call IO_error(0,ext_msg='tScalar_asInt')
      if(dNeq(s1%asFloat(),1.0_pReal)) call IO_error(0,ext_msg='tScalar_asFloat')
      s1 = 'True'
      if(.not. s1%asBool())            call IO_error(0,ext_msg='tScalar_asBool')
      if(s1%asString() /= 'True')      call IO_error(0,ext_msg='tScalar_asString')
  end select

  block
    class(tNode), pointer :: l1, l2, n
    select type(s1)
      class is(tScalar)
        s1 = '2'
    endselect
   
    select type(s2)
      class is(tScalar)
        s2 = '3'
    endselect

    allocate(tList::l1)
    select type(l1)
      class is(tList)
        call l1%append(s1)
        call l1%append(s2)
        n => l1
        if(any(l1%asInts() /= [2,3]))                            call IO_error(0,ext_msg='tList_asInts')
        if(any(dNeq(l1%asFloats(),[2.0_pReal,3.0_pReal])))       call IO_error(0,ext_msg='tList_asFloats')
        if(n%get_asInt(1) /= 2)                                  call IO_error(0,ext_msg='byIndex_asInt')
        if(dNeq(n%get_asFloat(2),3.0_pReal))                     call IO_error(0,ext_msg='byIndex_asFloat')
    endselect
    
    allocate(tList::l2)
    select type(l2)
      class is(tList)
        call l2%append(l1)
        if(any(l2%get_asInts(1) /= [2,3]))                       call IO_error(0,ext_msg='byIndex_asInts')
        if(any(dNeq(l2%get_asFloats(1),[2.0_pReal,3.0_pReal])))  call IO_error(0,ext_msg='byIndex_asFloats')
        n => l2
    end select
    deallocate(n)
   end block

   block
     type(tList),  target  :: l1
     type(tScalar),pointer :: s3,s4
     class(tNode), pointer :: n
    
     allocate(tScalar::s1)
     allocate(tScalar::s2)
     s3 => s1%asScalar()
     s4 => s2%asScalar()
     s3 = 'True'
     s4 = 'False'
    
     call l1%append(s1)
     call l1%append(s2)
     n => l1
     
     if(any(l1%asBools() .neqv. [.true., .false.]))               call IO_error(0,ext_msg='tList_asBools')
     if(any(l1%asStrings() /=   ['True ','False']))               call IO_error(0,ext_msg='tList_asStrings')
     if(n%get_asBool(2))                                          call IO_error(0,ext_msg='byIndex_asBool')
     if(n%get_asString(1) /= 'True')                              call IO_error(0,ext_msg='byIndex_asString')
   end block

end subroutine unitTest


!---------------------------------------------------------------------------------------------------
!> @brief init from string
!---------------------------------------------------------------------------------------------------
type(tScalar) pure function tScalar_init__(value)

  character(len=*), intent(in) ::  value

  tScalar_init__%value =value

end function tScalar_init__


!---------------------------------------------------------------------------------------------------
!> @brief set value from string
!---------------------------------------------------------------------------------------------------
elemental pure subroutine tScalar_assign__(self,value)

  type(tScalar),    intent(out) :: self
  character(len=*), intent(in)  :: value

  self%value = value

end subroutine tScalar_assign__


!--------------------------------------------------------------------------------------------------
!> @brief Type guard, guarantee scalar
!--------------------------------------------------------------------------------------------------
function tNode_asScalar(self) result(scalar)

  class(tNode),   intent(in), target :: self
  class(tScalar), pointer            :: scalar

  select type(self)
    class is(tScalar)
      scalar => self
    class default
      call IO_error(0)
  end select

end function tNode_asScalar


!--------------------------------------------------------------------------------------------------
!> @brief Type guard, guarantee list
!--------------------------------------------------------------------------------------------------
function tNode_asList(self) result(list)

  class(tNode), intent(in), target :: self
  class(tList), pointer            :: list

  select type(self)
    class is(tList)
      list => self
    class default
      call IO_error(0)
  end select

end function tNode_asList


!--------------------------------------------------------------------------------------------------
!> @brief Type guard, guarantee dict
!--------------------------------------------------------------------------------------------------
function tNode_asDict(self) result(dict)

  class(tNode), intent(in), target :: self
  class(tDict), pointer            :: dict

  select type(self)
    class is(tDict)
      dict => self
    class default
      call IO_error(0)
  end select

end function tNode_asDict


!--------------------------------------------------------------------------------------------------
!> @brief Access by index
!--------------------------------------------------------------------------------------------------
function tNode_get_byIndex(self,i) result(node)

  class(tNode), intent(in), target :: self
  integer,      intent(in)         :: i
  class(tNode),  pointer :: node

  class(tList),  pointer :: self_
  class(tItem),  pointer :: item
  integer :: j

  self_ => self%asList()
  if(i < 1 .or. i > self_%length) call IO_error(0)

  j = 1
  item => self_%first
  do while(j<i)
    item => item%next
    j = j + 1
  enddo
  node => item%node

end function tNode_get_byIndex


!--------------------------------------------------------------------------------------------------
!> @brief Access by index and convert to float
!--------------------------------------------------------------------------------------------------
function tNode_get_byIndex_asFloat(self,i) result(nodeAsFloat)

  class(tNode), intent(in), target :: self
  integer,      intent(in)         :: i
  real(pReal) :: nodeAsFloat

  class(tNode),  pointer :: node
  type(tScalar), pointer :: scalar

  node   => self%get(i)
  scalar => node%asScalar()
  nodeAsFloat = scalar%asFloat()

end function tNode_get_byIndex_asFloat


!--------------------------------------------------------------------------------------------------
!> @brief Access by index and convert to int
!--------------------------------------------------------------------------------------------------
function tNode_get_byIndex_asInt(self,i) result(nodeAsInt)

  class(tNode), intent(in), target :: self
  integer,      intent(in)         :: i
  integer :: nodeAsInt

  class(tNode),  pointer :: node
  type(tScalar), pointer :: scalar

  node   => self%get(i)
  scalar => node%asScalar()
  nodeAsInt = scalar%asInt()

end function tNode_get_byIndex_asInt


!--------------------------------------------------------------------------------------------------
!> @brief Access by index and convert to bool
!--------------------------------------------------------------------------------------------------
function tNode_get_byIndex_asBool(self,i) result(nodeAsBool)

  class(tNode), intent(in), target :: self
  integer,      intent(in)         :: i
  logical :: nodeAsBool

  class(tNode),  pointer :: node
  type(tScalar), pointer :: scalar

  node   => self%get(i)
  scalar => node%asScalar()
  nodeAsBool = scalar%asBool()

end function tNode_get_byIndex_asBool


!--------------------------------------------------------------------------------------------------
!> @brief Access by index and convert to string
!--------------------------------------------------------------------------------------------------
function tNode_get_byIndex_asString(self,i) result(nodeAsString)

  class(tNode), intent(in), target :: self
  integer,      intent(in)         :: i
  character(len=:), allocatable :: nodeAsString

  class(tNode),  pointer :: node
  type(tScalar), pointer :: scalar

  node   => self%get(i)
  scalar => node%asScalar()
  nodeAsString = scalar%asString()

end function tNode_get_byIndex_asString


!--------------------------------------------------------------------------------------------------
!> @brief Access by index and convert to float array
!--------------------------------------------------------------------------------------------------
function tNode_get_byIndex_asFloats(self,i) result(nodeAsFloats)

  class(tNode), intent(in), target :: self
  integer,      intent(in)         :: i
  real(pReal), dimension(:), allocatable :: nodeAsFloats

  class(tNode), pointer :: node
  class(tList),  pointer :: list

  node => self%get(i)
  list => node%asList()
  nodeAsFloats = list%asFloats()

end function tNode_get_byIndex_asFloats


!--------------------------------------------------------------------------------------------------
!> @brief Access by index and convert to int array
!--------------------------------------------------------------------------------------------------
function tNode_get_byIndex_asInts(self,i) result(nodeAsInts)

  class(tNode), intent(in), target :: self
  integer,      intent(in)         :: i
  integer, dimension(:), allocatable :: nodeAsInts

  class(tNode), pointer :: node
  class(tList), pointer :: list

  node => self%get(i)
  list => node%asList()
  nodeAsInts = list%asInts()

end function tNode_get_byIndex_asInts


!--------------------------------------------------------------------------------------------------
!> @brief Access by index and convert to bool array
!--------------------------------------------------------------------------------------------------
function tNode_get_byIndex_asBools(self,i) result(nodeAsBools)

  class(tNode), intent(in), target :: self
  integer,      intent(in)         :: i
  logical, dimension(:), allocatable :: nodeAsBools

  class(tNode), pointer :: node
  class(tList), pointer :: list

  node => self%get(i)
  list => node%asList()
  nodeAsBools = list%asBools()

end function tNode_get_byIndex_asBools


!--------------------------------------------------------------------------------------------------
!> @brief Access by index and convert to string array
!--------------------------------------------------------------------------------------------------
function tNode_get_byIndex_asStrings(self,i) result(nodeAsStrings)

  class(tNode), intent(in), target :: self
  integer,      intent(in)         :: i
  character(len=:), allocatable, dimension(:) :: nodeAsStrings

  class(tNode), pointer :: node
  type(tList),  pointer :: list

  node => self%get(i)
  list => node%asList()
  nodeAsStrings = list%asStrings()

end function tNode_get_byIndex_asStrings


!--------------------------------------------------------------------------------------------------
!> @brief Access by index
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey(self,k) result(node)

  class(tNode),     intent(in), target :: self
  character(len=*), intent(in)         :: k
  class(tNode),  pointer :: node

  type(tDict),  pointer :: self_
  type(tItem), pointer :: item
  integer :: j

  self_ => self%asDict()

  j = 1
  item => self_%first
  do while(j <= self_%length)
    if (item%key == k) exit
    item => item%next
    j = j + 1
  enddo
  if (.not. item%key == k) call IO_error(0)
  node => item%node

end function tNode_get_byKey


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to float
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asFloat(self,k) result(nodeAsFloat)

  class(tNode),     intent(in), target :: self
  character(len=*), intent(in)         :: k
  real(pReal) :: nodeAsFloat

  class(tNode),  pointer :: node
  type(tScalar), pointer :: scalar

  node   => self%get(k)
  scalar => node%asScalar()
  nodeAsFloat = scalar%asFloat()

end function tNode_get_byKey_asFloat


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to int
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asInt(self,k) result(nodeAsInt)

  class(tNode),     intent(in), target :: self
  character(len=*), intent(in)         :: k
  integer :: nodeAsInt

  class(tNode),  pointer :: node
  type(tScalar), pointer :: scalar

  node   => self%get(k)
  scalar => node%asScalar()
  nodeAsInt = scalar%asInt()

end function tNode_get_byKey_asInt


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to bool
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asBool(self,k) result(nodeAsBool)

  class(tNode),     intent(in), target :: self
  character(len=*), intent(in)         :: k
  logical :: nodeAsBool

  class(tNode),  pointer :: node
  type(tScalar), pointer :: scalar

  node   => self%get(k)
  scalar => node%asScalar()
  nodeAsBool = scalar%asBool()

end function tNode_get_byKey_asBool


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to string
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asString(self,k) result(nodeAsString)

  class(tNode),     intent(in), target :: self
  character(len=*), intent(in)         :: k
  character(len=:), allocatable :: nodeAsString

  class(tNode),  pointer :: node
  type(tScalar), pointer :: scalar

  node   => self%get(k)
  scalar => node%asScalar()
  nodeAsString = scalar%asString()

end function tNode_get_byKey_asString


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to float array
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asFloats(self,k) result(nodeAsFloats)

  class(tNode),     intent(in), target :: self
  character(len=*), intent(in)         :: k
  real(pReal), dimension(:), allocatable :: nodeAsFloats

  class(tNode), pointer :: node
  type(tList),  pointer :: list

  node   => self%get(k)
  list => node%asList()
  nodeAsFloats = list%asFloats()

end function tNode_get_byKey_asFloats


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to int array
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asInts(self,k) result(nodeAsInts)

  class(tNode),     intent(in), target :: self
  character(len=*), intent(in)         :: k
  integer, dimension(:), allocatable :: nodeAsInts

  class(tNode), pointer :: node
  type(tList),  pointer :: list

  node   => self%get(k)
  list => node%asList()
  nodeAsInts = list%asInts()

end function tNode_get_byKey_asInts


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to bool array
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asBools(self,k) result(nodeAsBools)

  class(tNode),     intent(in), target :: self
  character(len=*), intent(in)         :: k
  logical, dimension(:), allocatable :: nodeAsBools

  class(tNode), pointer :: node
  type(tList),  pointer :: list

  node   => self%get(k)
  list => node%asList()
  nodeAsBools = list%asBools()

end function tNode_get_byKey_asBools


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to string array
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asStrings(self,k) result(nodeAsStrings)

  class(tNode),     intent(in), target :: self
  character(len=*), intent(in)         :: k
  character(len=:), allocatable, dimension(:) :: nodeAsStrings

  class(tNode), pointer :: node
  type(tList),  pointer :: list

  node   => self%get(k)
  list => node%asList()
  nodeAsStrings = list%asStrings()

end function tNode_get_byKey_asStrings


!-------------------------------------------------------------------------------------------------------
!> @brief Returns the index of a key in a dictionary
!-------------------------------------------------------------------------------------------------------
function tNode_get_byKey_asIndex(self,key)  result(keyIndex)

  class(tNode),     intent(in), target  :: self
  character(len=*), intent(in)          :: key

  integer              :: keyIndex
  integer              :: i
  type(tDict), pointer :: dict
  type(tItem), pointer :: item

  dict => self%asDict()
  item => dict%first
  do i = 1, dict%length
    if(key == item%key) then
      keyIndex = i
      exit
    else
      item => item%next
    endif
  enddo

end function tNode_get_byKey_asIndex


!--------------------------------------------------------------------------------------------------
!> @brief Prints scalar as string
!--------------------------------------------------------------------------------------------------
recursive subroutine tScalar_asFormattedString(self,indent)

  class (tScalar), intent(in), target    :: self
  integer,         intent(in), optional  :: indent

  integer :: indent_

  if(present(indent)) then
    indent_ = indent
  else
    indent_ = 0
  endif

  write (6,'(a)') trim(self%value)

end subroutine tScalar_asFormattedString


!--------------------------------------------------------------------------------------------------
!> @brief Prints list as string (YAML block style)
!--------------------------------------------------------------------------------------------------
recursive subroutine tList_asFormattedString(self,indent)

  class (tList),intent(in),target      :: self
  integer,      intent(in),optional    :: indent

  type (tItem), pointer  :: item
  integer :: i, indent_

  if(present(indent)) then
    indent_ = indent
  else
    indent_ = 0
  endif

  item => self%first
  do i = 1, self%length
    if( i /= 1) write (6,'(a)',advance='NO') repeat(' ',indent_)
    write (6,'(a)',advance='NO') '- '
    call item%node%asFormattedString(indent_+2)
    item => item%next
  end do

end subroutine tList_asFormattedString


!--------------------------------------------------------------------------------------------------
!> @brief Prints dictionary as string (YAML block style)
!--------------------------------------------------------------------------------------------------
recursive subroutine tDict_asFormattedString(self,indent)

  class (tDict),intent(in),target      :: self
  integer,      intent(in),optional    :: indent
  
  type (tItem),pointer  :: item
  integer :: i, indent_

  if(present(indent)) then
    indent_ = indent
  else
    indent_ = 0
  endif

  item => self%first
  do i = 1, self%length
    if( i /= 1) write (6,'(a)',advance='NO') repeat(' ',indent_)
    select type (node_ => item%node)
      class is (tScalar)
        write (6,'(a)',advance='NO') trim(item%key)//': '
        call node_%asFormattedString(indent_+len_trim(item%key)+2)
      class default
        write (6,'(a)') trim(item%key)//':'
        write (6,'(a)',advance='NO') repeat(' ',indent_+2)
        call node_%asFormattedString(indent_+2)
    end select
    item => item%next
  end do

end subroutine tDict_asFormattedString


!--------------------------------------------------------------------------------------------------
!> @brief Convert to float
!--------------------------------------------------------------------------------------------------
function tScalar_asFloat(self)

  class(tScalar), intent(in), target :: self
  real(pReal) :: tScalar_asFloat

  tScalar_asFloat = IO_stringAsFloat(self%value)

end function tScalar_asFloat


!--------------------------------------------------------------------------------------------------
!> @brief Convert to int
!--------------------------------------------------------------------------------------------------
function tScalar_asInt(self)

  class(tScalar), intent(in), target :: self
  integer :: tScalar_asInt

  tScalar_asInt = IO_stringAsInt(self%value)

end function tScalar_asInt


!--------------------------------------------------------------------------------------------------
!> @brief Convert to bool
!--------------------------------------------------------------------------------------------------
function tScalar_asBool(self)

  class(tScalar), intent(in), target :: self
  logical :: tScalar_asBool

  tScalar_asBool = IO_stringAsBool(self%value)

end function tScalar_asBool


!--------------------------------------------------------------------------------------------------
!> @brief Convert to string
!--------------------------------------------------------------------------------------------------
function tScalar_asString(self)

  class(tScalar), intent(in), target :: self
  character(len=:), allocatable :: tScalar_asString

  tScalar_asString = self%value

end function tScalar_asString


!--------------------------------------------------------------------------------------------------
!> @brief Convert to float array
!--------------------------------------------------------------------------------------------------
function tList_asFloats(self)

  class(tList), intent(in), target :: self
  real(pReal), dimension(:), allocatable :: tList_asFloats

  integer :: i
  type(tItem),   pointer :: item
  type(tScalar), pointer :: scalar

  allocate(tList_asFloats(self%length))
  item => self%first
  do i = 1, self%length
    scalar => item%node%asScalar()
    tList_asFloats(i) = scalar%asFloat()
    item => item%next
  enddo

end function tList_asFloats


!--------------------------------------------------------------------------------------------------
!> @brief Convert to int array
!--------------------------------------------------------------------------------------------------
function tList_asInts(self)

  class(tList), intent(in), target :: self
  integer, dimension(:), allocatable :: tList_asInts

  integer :: i
  type(tItem),   pointer :: item
  type(tScalar), pointer :: scalar

  allocate(tList_asInts(self%length))
  item => self%first
  do i = 1, self%length
    scalar => item%node%asScalar()
    tList_asInts(i) = scalar%asInt()
    item => item%next
  enddo

end function tList_asInts


!--------------------------------------------------------------------------------------------------
!> @brief Convert to bool array
!--------------------------------------------------------------------------------------------------
function tList_asBools(self)

  class(tList), intent(in), target :: self
  logical, dimension(:), allocatable :: tList_asBools

  integer :: i
  type(tItem),   pointer :: item
  type(tScalar), pointer :: scalar

  allocate(tList_asBools(self%length))
  item => self%first
  do i = 1, self%length
    scalar => item%node%asScalar()
    tList_asBools(i) = scalar%asBool()
    item => item%next
  enddo

end function tList_asBools


!--------------------------------------------------------------------------------------------------
!> @brief Convert to string array
!--------------------------------------------------------------------------------------------------
function tList_asStrings(self)

  class(tList), intent(in), target :: self
  character(len=:), allocatable, dimension(:) :: tList_asStrings

  integer :: i,len_max
  type(tItem),   pointer :: item
  type(tScalar), pointer :: scalar

  len_max = 0
  allocate(character(len=pStringLen) :: tList_asStrings(self%length))
  item => self%first
  do i = 1, self%length
    scalar => item%node%asScalar()
    tList_asStrings(i) = scalar%asString()
    len_max = max(len_max, len_trim(tList_asStrings(i)))
    item => item%next
  enddo

  !ToDo: trim to len_max

end function tList_asStrings


!--------------------------------------------------------------------------------------------------
!> @brief Append element
!--------------------------------------------------------------------------------------------------
subroutine tList_append(self,node)

  class(tList), intent(inout)         :: self
  class(tNode), intent(in), target    :: node

  type(tItem), pointer :: item

  if (.not. associated(self%first)) then
    allocate(self%first)
    item => self%first
  else
    item => self%first
    do while (associated(item%next))
      item => item%next
    end do
    allocate(item%next)
    item => item%next
  end if

  item%node => node
  self%length = self%length + 1

end subroutine tList_append


!--------------------------------------------------------------------------------------------------
!> @brief Set the value of a key (either replace or add new)
!--------------------------------------------------------------------------------------------------
subroutine tDict_set(self,key,node)

  class (tDict),    intent(inout)         :: self
  character(len=*), intent(in)            :: key
  class(tNode),     intent(in), target    :: node

  type(tItem), pointer :: item

  if (.not. associated(self%first)) then
    allocate(self%first)
    item => self%first
    self%length = 1
  else
    item => self%first
    searchExisting: do while (associated(item%next))
      if (item%key == key) exit
      item => item%next
    end do searchExisting
    if (.not. item%key == key) then
      allocate(item%next)
      item => item%next
      self%length = self%length + 1
    end if
  end if

  item%key = key
  item%node => node

end subroutine tDict_set


!--------------------------------------------------------------------------------------------------
!> @brief empties lists and dicts and free associated memory
!> @details called when variable goes out of scope.
!--------------------------------------------------------------------------------------------------
recursive subroutine tList_finalize(self)

  type (tList),intent(inout) :: self

  deallocate(self%first)

end subroutine tList_finalize


!--------------------------------------------------------------------------------------------------
!> @brief empties nodes and frees associated memory
!--------------------------------------------------------------------------------------------------
recursive subroutine tItem_finalize(self)

  type(tItem),intent(inout) :: self
  
  deallocate(self%node)
  if(associated(self%next)) deallocate(self%next)

end subroutine tItem_finalize

end module YAML_types
