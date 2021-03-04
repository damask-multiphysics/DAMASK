!--------------------------------------------------------------------------------------------------
!> @author Sharan Roongta, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Data types to create a scalar, a list, and a dictionary/hash
!> @details module describes the various functions to store and get the yaml data.
!! A node is the base class for scalar, list and dictionary, list items and dictionary entries point
!! to a node.
!--------------------------------------------------------------------------------------------------

module YAML_types
  use IO
  use prec

  implicit none
  private

  type, abstract, public :: tNode
    integer :: length = 0
    contains
    procedure(asFormattedString), deferred :: asFormattedString
    procedure :: &
      asScalar     => tNode_asScalar
    procedure :: &
      isScalar     => tNode_isScalar
    procedure :: &
      asList       => tNode_asList
    procedure :: &
      isList       => tNode_isList
    procedure :: &
      asDict       => tNode_asDict
    procedure :: &
      isDict       => tNode_isDict
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
    procedure :: &
      getKey                      => tNode_getKey_byIndex
    procedure :: &
      contains                    => tNode_contains

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


  type, extends(tNode), public :: tScalar

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

  type, extends(tNode), public :: tList

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

  type, extends(tList), public :: tDict
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

  type(tDict), target, public :: &
    emptyDict
  type(tList), target, public :: &
    emptyList

  abstract interface

    recursive function asFormattedString(self,indent)
      import tNode
      character(len=:), allocatable      :: asFormattedString
      class(tNode), intent(in), target   :: self
      integer,      intent(in), optional :: indent
    end function asFormattedString

  end interface

  interface tScalar
    module procedure tScalar_init__
  end interface tScalar

  interface assignment (=)
    module procedure tScalar_assign__
  end interface assignment (=)

  public :: &
    YAML_types_init, &
    output_asStrings, &                                        !ToDo: Hack for GNU. Remove later
    assignment(=)

contains

!--------------------------------------------------------------------------------------------------
!> @brief Do sanity checks.
!--------------------------------------------------------------------------------------------------
subroutine YAML_types_init

  print'(/,a)', ' <<<+-  YAML_types init  -+>>>'

  call selfTest

end subroutine YAML_types_init


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some type bound procedures.
!--------------------------------------------------------------------------------------------------
subroutine selfTest

  class(tNode), pointer  :: s1,s2
  allocate(tScalar::s1)
  allocate(tScalar::s2)
  select type(s1)
    class is(tScalar)
      s1 = '1'
      if (s1%asInt() /= 1)              error stop 'tScalar_asInt'
      if (dNeq(s1%asFloat(),1.0_pReal)) error stop 'tScalar_asFloat'
      s1 = 'true'
      if (.not. s1%asBool())            error stop 'tScalar_asBool'
      if (s1%asString() /= 'true')      error stop 'tScalar_asString'
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
        if (any(l1%asInts() /= [2,3]))                      error stop 'tList_asInts'
        if (any(dNeq(l1%asFloats(),[2.0_pReal,3.0_pReal]))) error stop 'tList_asFloats'
        if (n%get_asInt(1) /= 2)                            error stop 'byIndex_asInt'
        if (dNeq(n%get_asFloat(2),3.0_pReal))               error stop 'byIndex_asFloat'
    endselect

    allocate(tList::l2)
    select type(l2)
      class is(tList)
        call l2%append(l1)
        if (any(l2%get_asInts(1) /= [2,3]))                      error stop 'byIndex_asInts'
        if (any(dNeq(l2%get_asFloats(1),[2.0_pReal,3.0_pReal]))) error stop 'byIndex_asFloats'
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
     s3 = 'true'
     s4 = 'False'

     call l1%append(s1)
     call l1%append(s2)
     n => l1

     if (any(l1%asBools() .neqv. [.true., .false.])) error stop 'tList_asBools'
     if (any(l1%asStrings() /=   ['true ','False'])) error stop 'tList_asStrings'
     if (n%get_asBool(2))                            error stop 'byIndex_asBool'
     if (n%get_asString(1) /= 'true')                error stop 'byIndex_asString'
  end block

end subroutine selfTest


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
      call IO_error(706,ext_msg='Expected "scalar"')
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
      call IO_error(706,ext_msg='Expected "list"')
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
      call IO_error(706,ext_msg='Expected "dict"')
  end select

end function tNode_asDict


!--------------------------------------------------------------------------------------------------
!> @brief Checks if node is a scalar
!--------------------------------------------------------------------------------------------------
function tNode_isScalar(self) result(scalar)

  class(tNode),   intent(in), target :: self
  logical                            :: scalar

  scalar = .false.
  select type(self)
    class is(tScalar)
      scalar = .true.
  end select

end function tNode_isScalar


!--------------------------------------------------------------------------------------------------
!> @brief Checks if node is a list
!--------------------------------------------------------------------------------------------------
function tNode_isList(self) result(list)

  class(tNode), intent(in), target :: self
  logical                          :: list

  list = .false.
  select type(self)
    class is(tList)
      list = .true.
  end select

end function tNode_isList


!--------------------------------------------------------------------------------------------------
!> @brief Checks if node is a dict
!--------------------------------------------------------------------------------------------------
function tNode_isDict(self) result(dict)

  class(tNode), intent(in), target :: self
  logical                          :: dict

  dict = .false.
  select type(self)
    class is(tDict)
      dict = .true.
  end select

end function tNode_isDict


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
  if (i < 1 .or. i > self_%length) call IO_error(150,ext_msg='tNode_get_byIndex')

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
!> @brief Returns the key in a dictionary as a string
!--------------------------------------------------------------------------------------------------
function tNode_getKey_byIndex(self,i)  result(key)

  class(tNode),     intent(in), target  :: self
  integer,          intent(in)          :: i

  character(len=:), allocatable         :: key
  integer              :: j
  type(tDict), pointer :: dict
  type(tItem), pointer :: item

  dict => self%asDict()
  item => dict%first
  do j = 1, dict%length
    if (j == i) then
      key = item%key
      exit
    else
      item => item%next
    endif
  enddo

end function tNode_getKey_byIndex


!-------------------------------------------------------------------------------------------------
!> @brief Checks if a given key/item is present in the dict/list
!-------------------------------------------------------------------------------------------------
function tNode_contains(self,k)  result(exists)

  class(tNode),     intent(in), target  :: self
  character(len=*), intent(in)          :: k

  logical                               :: exists
  integer   :: j
  type(tList), pointer :: list
  type(tDict), pointer :: dict

  exists = .false.
  if (self%isDict()) then
    dict => self%asDict()
    do j=1, dict%length
      if (dict%getKey(j) == k) then
        exists = .true.
        return
      endif
    enddo
  elseif (self%isList()) then
    list => self%asList()
    do j=1, list%length
      if (list%get_asString(j) == k) then
        exists = .true.
        return
      endif
    enddo
  else
     call IO_error(706,ext_msg='Expected "list" or "dict"')
  endif

end function tNode_contains


!--------------------------------------------------------------------------------------------------
!> @brief Access by key
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey(self,k,defaultVal) result(node)

  class(tNode),     intent(in), target         :: self
  character(len=*), intent(in)                 :: k
  class(tNode),     intent(in),optional,target :: defaultVal
  class(tNode),     pointer :: node

  type(tDict), pointer :: self_
  type(tItem), pointer :: item
  integer :: j
  logical :: found

  found = present(defaultVal)
  if (found) node => defaultVal

  self_ => self%asDict()

  j = 1
  item => self_%first
  do while(j <= self_%length)
    if (item%key == k) then
      found = .true.
      exit
    endif
    item => item%next
    j = j + 1
  enddo

  if (.not. found) then
    call IO_error(143,ext_msg=k)
  else
    if (associated(item)) node => item%node
  endif

end function tNode_get_byKey


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to float
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asFloat(self,k,defaultVal) result(nodeAsFloat)

  class(tNode),     intent(in), target  :: self
  character(len=*), intent(in)          :: k
  real(pReal),      intent(in),optional :: defaultVal
  real(pReal) :: nodeAsFloat

  class(tNode),  pointer :: node
  type(tScalar), pointer :: scalar

  if (self%contains(k)) then
    node => self%get(k)
    scalar => node%asScalar()
    nodeAsFloat = scalar%asFloat()
  elseif (present(defaultVal)) then
    nodeAsFloat = defaultVal
  else
    call IO_error(143,ext_msg=k)
  endif

end function tNode_get_byKey_asFloat


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to int
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asInt(self,k,defaultVal) result(nodeAsInt)

  class(tNode),     intent(in), target  :: self
  character(len=*), intent(in)          :: k
  integer,          intent(in),optional :: defaultVal
  integer :: nodeAsInt

  class(tNode),  pointer :: node
  type(tScalar), pointer :: scalar

  if (self%contains(k)) then
    node => self%get(k)
    scalar => node%asScalar()
    nodeAsInt = scalar%asInt()
  elseif (present(defaultVal)) then
    nodeAsInt = defaultVal
  else
    call IO_error(143,ext_msg=k)
  endif

end function tNode_get_byKey_asInt


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to bool
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asBool(self,k,defaultVal) result(nodeAsBool)

  class(tNode),     intent(in), target  :: self
  character(len=*), intent(in)          :: k
  logical,          intent(in),optional :: defaultVal
  logical :: nodeAsBool

  class(tNode),  pointer :: node
  type(tScalar), pointer :: scalar

  if (self%contains(k)) then
    node => self%get(k)
    scalar => node%asScalar()
    nodeAsBool = scalar%asBool()
  elseif (present(defaultVal)) then
    nodeAsBool = defaultVal
  else
    call IO_error(143,ext_msg=k)
  endif

end function tNode_get_byKey_asBool


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to string
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asString(self,k,defaultVal) result(nodeAsString)

  class(tNode),     intent(in), target  :: self
  character(len=*), intent(in)          :: k
  character(len=*), intent(in),optional :: defaultVal
  character(len=:), allocatable :: nodeAsString

  class(tNode),  pointer :: node
  type(tScalar), pointer :: scalar

  if (self%contains(k)) then
    node => self%get(k)
    scalar => node%asScalar()
    nodeAsString = scalar%asString()
  elseif (present(defaultVal)) then
    nodeAsString = defaultVal
  else
    call IO_error(143,ext_msg=k)
  endif

end function tNode_get_byKey_asString


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to float array
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asFloats(self,k,defaultVal,requiredSize) result(nodeAsFloats)

  class(tNode),     intent(in), target                 :: self
  character(len=*), intent(in)                         :: k
  real(pReal),      intent(in), dimension(:), optional :: defaultVal
  integer,          intent(in),               optional :: requiredSize

  real(pReal), dimension(:), allocatable :: nodeAsFloats

  class(tNode), pointer :: node
  type(tList),  pointer :: list

  if (self%contains(k)) then
    node => self%get(k)
    list => node%asList()
    nodeAsFloats = list%asFloats()
  elseif (present(defaultVal)) then
    nodeAsFloats = defaultVal
  else
    call IO_error(143,ext_msg=k)
  endif

  if (present(requiredSize)) then
    if (requiredSize /= size(nodeAsFloats)) call IO_error(146,ext_msg=k)
  endif

end function tNode_get_byKey_asFloats


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to int array
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asInts(self,k,defaultVal,requiredSize) result(nodeAsInts)

  class(tNode),          intent(in), target   :: self
  character(len=*),      intent(in)           :: k
  integer, dimension(:), intent(in), optional :: defaultVal
  integer,               intent(in), optional :: requiredSize
  integer, dimension(:), allocatable :: nodeAsInts

  class(tNode), pointer :: node
  type(tList),  pointer :: list

  if (self%contains(k)) then
    node => self%get(k)
    list => node%asList()
    nodeAsInts = list%asInts()
  elseif (present(defaultVal)) then
    nodeAsInts = defaultVal
  else
    call IO_error(143,ext_msg=k)
  endif

  if (present(requiredSize)) then
    if (requiredSize /= size(nodeAsInts)) call IO_error(146,ext_msg=k)
  endif

end function tNode_get_byKey_asInts


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to bool array
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asBools(self,k,defaultVal) result(nodeAsBools)

  class(tNode),          intent(in), target   :: self
  character(len=*),      intent(in)           :: k
  logical, dimension(:), intent(in), optional :: defaultVal
  logical, dimension(:), allocatable          :: nodeAsBools

  class(tNode), pointer :: node
  type(tList),  pointer :: list

  if (self%contains(k)) then
    node => self%get(k)
    list => node%asList()
    nodeAsBools = list%asBools()
  elseif (present(defaultVal)) then
    nodeAsBools = defaultVal
  else
    call IO_error(143,ext_msg=k)
  endif

end function tNode_get_byKey_asBools


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to string array
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asStrings(self,k,defaultVal) result(nodeAsStrings)

  class(tNode),     intent(in), target                 :: self
  character(len=*), intent(in)                         :: k
  character(len=*), intent(in), dimension(:), optional :: defaultVal
  character(len=:), allocatable, dimension(:)          :: nodeAsStrings

  class(tNode), pointer :: node
  type(tList),  pointer :: list

  if (self%contains(k)) then
    node => self%get(k)
    list => node%asList()
    nodeAsStrings = list%asStrings()
  elseif (present(defaultVal)) then
    nodeAsStrings = defaultVal
  else
    call IO_error(143,ext_msg=k)
  endif

end function tNode_get_byKey_asStrings


!--------------------------------------------------------------------------------------------------
!> @brief Returns string output array (hack for GNU)
!--------------------------------------------------------------------------------------------------
function output_asStrings(self)  result(output)                   !ToDo: SR: Remove whenever GNU works

  class(tNode), pointer,intent(in)   ::  self
  character(len=pStringLen), allocatable, dimension(:) :: output

  class(tNode), pointer :: output_list
  integer :: o

  output_list   => self%get('output',defaultVal=emptyList)
  allocate(output(output_list%length))
  do o = 1, output_list%length
    output(o) = output_list%get_asString(o)
  enddo


end function output_asStrings


!--------------------------------------------------------------------------------------------------
!> @brief Returns the index of a key in a dictionary
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asIndex(self,key)  result(keyIndex)

  class(tNode),     intent(in), target  :: self
  character(len=*), intent(in)          :: key

  integer              :: keyIndex
  integer              :: i
  type(tDict), pointer :: dict
  type(tItem), pointer :: item

  dict => self%asDict()
  item => dict%first
  keyIndex = -1
  do i = 1, dict%length
    if (key == item%key) then
      keyIndex = i
      exit
    else
      item => item%next
    endif
  enddo

  if (keyIndex == -1) call IO_error(140,ext_msg=key)


end function tNode_get_byKey_asIndex


!--------------------------------------------------------------------------------------------------
!> @brief Scalar as string (YAML block style)
!--------------------------------------------------------------------------------------------------
recursive function tScalar_asFormattedString(self,indent)

  character(len=:), allocatable          :: tScalar_asFormattedString
  class (tScalar), intent(in), target    :: self
  integer,         intent(in), optional  :: indent

  tScalar_asFormattedString = trim(self%value)//IO_EOL

end function tScalar_asFormattedString


!--------------------------------------------------------------------------------------------------
!> @brief List as string (YAML block style)
!--------------------------------------------------------------------------------------------------
recursive function tList_asFormattedString(self,indent) result(str)

  class (tList),intent(in),target      :: self
  integer,      intent(in),optional    :: indent

  type (tItem), pointer  :: item
  character(len=:), allocatable :: str
  integer :: i, indent_

  str = ''
  if (present(indent)) then
    indent_ = indent
  else
    indent_ = 0
  endif

  item => self%first
  do i = 1, self%length
    if (i /= 1) str = str//repeat(' ',indent_)
    str = str//'- '//item%node%asFormattedString(indent_+2)
    item => item%next
  end do

end function tList_asFormattedString


!--------------------------------------------------------------------------------------------------
!> @brief Dictionary as string (YAML block style)
!--------------------------------------------------------------------------------------------------
recursive function tDict_asFormattedString(self,indent) result(str)

  class (tDict),intent(in),target      :: self
  integer,      intent(in),optional    :: indent

  type (tItem),pointer  :: item
  character(len=:), allocatable :: str
  integer :: i, indent_

  str = ''
  if (present(indent)) then
    indent_ = indent
  else
    indent_ = 0
  endif

  item => self%first
  do i = 1, self%length
    if (i /= 1) str = str//repeat(' ',indent_)
    select type(node_1 =>item%node)
      class is(tScalar)
        str = str//trim(item%key)//': '//item%node%asFormattedString(indent_+len_trim(item%key)+2)
      class default
        str = str//trim(item%key)//':'//IO_EOL//repeat(' ',indent_+2)//item%node%asFormattedString(indent_+2)
    endselect
    item => item%next
  end do

end function tDict_asFormattedString


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
  item => self%first
  do i = 1, self%length
    scalar => item%node%asScalar()
    len_max = max(len_max, len_trim(scalar%asString()))
    item => item%next
  enddo

  allocate(character(len=len_max) :: tList_asStrings(self%length))
  item => self%first
  do i = 1, self%length
    scalar => item%node%asScalar()
    tList_asStrings(i) = scalar%asString()
    item => item%next
  enddo

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
  if (associated(self%next)) deallocate(self%next)

end subroutine tItem_finalize

end module YAML_types
