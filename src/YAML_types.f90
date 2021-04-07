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
      asList       => tNode_asList
    procedure :: &
      asDict       => tNode_asDict
    procedure :: &
      tNode_get_byIndex            => tNode_get_byIndex
    procedure :: &
      tNode_get_byIndex_asFloat    => tNode_get_byIndex_asFloat
    procedure :: &
      tNode_get_byIndex_as1dFloat  => tNode_get_byIndex_as1dFloat
    procedure :: &
      tNode_get_byIndex_asInt      => tNode_get_byIndex_asInt
    procedure :: &
      tNode_get_byIndex_as1dInt    => tNode_get_byIndex_as1dInt
    procedure :: &
      tNode_get_byIndex_asBool     => tNode_get_byIndex_asBool
    procedure :: &
      tNode_get_byIndex_as1dBool   => tNode_get_byIndex_as1dBool
    procedure :: &
      tNode_get_byIndex_asString   => tNode_get_byIndex_asString
    procedure :: &
      tNode_get_byIndex_as1dString => tNode_get_byIndex_as1dString
    procedure :: &
      tNode_get_byKey              => tNode_get_byKey
    procedure :: &
      tNode_get_byKey_asFloat      => tNode_get_byKey_asFloat
    procedure :: &
      tNode_get_byKey_as1dFloat    => tNode_get_byKey_as1dFloat
    procedure :: &
      tNode_get_byKey_asInt        => tNode_get_byKey_asInt
    procedure :: &
      tNode_get_byKey_as1dInt      => tNode_get_byKey_as1dInt
    procedure :: &
      tNode_get_byKey_asBool       => tNode_get_byKey_asBool
    procedure :: &
      tNode_get_byKey_as1dBool     => tNode_get_byKey_as1dBool
    procedure :: &
      tNode_get_byKey_asString     => tNode_get_byKey_asString
    procedure :: &
      tNode_get_byKey_as1dString   => tNode_get_byKey_as1dString
    procedure :: &
      getKey                       => tNode_get_byIndex_asKey
    procedure :: &
      getIndex                     => tNode_get_byKey_asIndex
    procedure :: &
      contains                     => tNode_contains
    procedure :: &
      get_as2dFloat                => tNode_get_byKey_as2dFloat

    generic :: &
      get            => tNode_get_byIndex, &
                        tNode_get_byKey
    generic :: &
      get_asFloat    => tNode_get_byIndex_asFloat, &
                        tNode_get_byKey_asFloat
    generic :: &
      get_as1dFloat  => tNode_get_byIndex_as1dFloat, &
                        tNode_get_byKey_as1dFloat
    generic :: &
      get_asInt      => tNode_get_byIndex_asInt, &
                        tNode_get_byKey_asInt
    generic :: &
      get_as1dInt    => tNode_get_byIndex_as1dInt, &
                        tNode_get_byKey_as1dInt
    generic :: &
      get_asBool     => tNode_get_byIndex_asBool, &
                        tNode_get_byKey_asBool
    generic :: &
      get_as1dBool   => tNode_get_byIndex_as1dBool, &
                        tNode_get_byKey_as1dBool
    generic :: &
      get_asString   => tNode_get_byIndex_asString, &
                        tNode_get_byKey_asString
    generic :: &
      get_as1dString => tNode_get_byIndex_as1dString, &
                        tNode_get_byKey_as1dString
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
      as1dFloat  => tList_as1dFloat
    procedure :: &
      as2dFloat  => tList_as2dFloat
    procedure :: &
      as1dInt    => tList_as1dInt
    procedure :: &
      as1dBool   => tList_as1dBool
    procedure :: &
      as1dString => tList_as1dString
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
    output_as1dString, &                                       !ToDo: Hack for GNU. Remove later
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

  class(tNode), pointer  :: s1,s2,s3,s4
  allocate(tScalar::s1)
  allocate(tScalar::s2)
  allocate(tScalar::s3)
  allocate(tScalar::s4)
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
    class(tNode), pointer :: l1, l2, l3, n
    real(pReal), allocatable, dimension(:,:) :: x

    select type(s1)
      class is(tScalar)
        s1 = '2'
    endselect

    select type(s2)
      class is(tScalar)
        s2 = '3'
    endselect

    select type(s3)
      class is(tScalar)
        s3 = '4'
    endselect

    select type(s4)
      class is(tScalar)
        s4 = '5'
    endselect


    allocate(tList::l1)
    select type(l1)
      class is(tList)
        call l1%append(s1)
        call l1%append(s2)
        n => l1
        if (any(l1%as1dInt() /= [2,3]))                      error stop 'tList_as1dInt'
        if (any(dNeq(l1%as1dFloat(),[2.0_pReal,3.0_pReal]))) error stop 'tList_as1dFloat'
        if (n%get_asInt(1) /= 2)                             error stop 'byIndex_asInt'
        if (dNeq(n%get_asFloat(2),3.0_pReal))                error stop 'byIndex_asFloat'
    endselect

    allocate(tList::l3)
    select type(l3)
      class is(tList)
        call l3%append(s3)
        call l3%append(s4)
    endselect

    allocate(tList::l2)
    select type(l2)
      class is(tList)
        call l2%append(l1)
        if(any(l2%get_as1dInt(1) /= [2,3]))                      error stop 'byIndex_as1dInt'
        if(any(dNeq(l2%get_as1dFloat(1),[2.0_pReal,3.0_pReal]))) error stop 'byIndex_as1dFloat'
        call l2%append(l3)
        x = l2%as2dFloat()
        if(dNeq(x(2,1),4.0_pReal))                               error stop 'byKey_as2dFloat'
        if(any(dNeq(pack(l2%as2dFloat(),.true.),&
               [2.0_pReal,4.0_pReal,3.0_pReal,5.0_pReal])))      error stop 'byKey_as2dFloat'
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

     if (any(l1%as1dBool() .neqv. [.true., .false.])) error stop 'tList_as1dBool'
     if (any(l1%as1dString() /=   ['true ','False'])) error stop 'tList_as1dString'
     if (n%get_asBool(2))                             error stop 'byIndex_asBool'
     if (n%get_asString(1) /= 'true')                 error stop 'byIndex_asString'
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

  select type(self)
    class is(tList)
      self_ => self%asList()
    class default
      call IO_error(706,ext_msg='Expected list')
  endselect

  item => self_%first

  if (i < 1 .or. i > self_%length) call IO_error(150,ext_msg='tNode_get_byIndex')

  do j = 2,i
    item => item%next
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
  select type(node)
    class is(tScalar)
      scalar => node%asScalar()
      nodeAsFloat = scalar%asFloat()
    class default
      call IO_error(706,ext_msg='Expected scalar float')
  end select

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
  select type(node)
    class is(tScalar)
      scalar => node%asScalar()
      nodeAsInt = scalar%asInt()
    class default
      call IO_error(706,ext_msg='Expected scalar integer')
  end select

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
  select type(node)
    class is(tScalar)
      scalar => node%asScalar()
      nodeAsBool = scalar%asBool()
    class default
      call IO_error(706,ext_msg='Expected scalar Boolean')
  endselect

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
  select type(node)
    class is(tScalar)
      scalar => node%asScalar()
      nodeAsString = scalar%asString()
    class default
      call IO_error(706,ext_msg='Expected scalar string')
  endselect

end function tNode_get_byIndex_asString


!--------------------------------------------------------------------------------------------------
!> @brief Access by index and convert to float array (1D)
!--------------------------------------------------------------------------------------------------
function tNode_get_byIndex_as1dFloat(self,i) result(nodeAs1dFloat)

  class(tNode), intent(in), target :: self
  integer,      intent(in)         :: i
  real(pReal), dimension(:), allocatable :: nodeAs1dFloat

  class(tNode), pointer  :: node
  class(tList),  pointer :: list

  node => self%get(i)
  select type(node)
    class is(tList)
      list => node%asList()
      nodeAs1dFloat = list%as1dFloat()
    class default
      call IO_error(706,ext_msg='Expected list of floats')
  endselect

end function tNode_get_byIndex_as1dFloat


!--------------------------------------------------------------------------------------------------
!> @brief Access by index and convert to int array (1D)
!--------------------------------------------------------------------------------------------------
function tNode_get_byIndex_as1dInt(self,i) result(nodeAs1dInt)

  class(tNode), intent(in), target :: self
  integer,      intent(in)         :: i
  integer, dimension(:), allocatable :: nodeAs1dInt

  class(tNode), pointer :: node
  class(tList), pointer :: list

  node => self%get(i)
  select type(node)
    class is(tList)
      list => node%asList()
      nodeAs1dInt = list%as1dInt()
    class default
      call IO_error(706,ext_msg='Expected list of integers')
  endselect

end function tNode_get_byIndex_as1dInt


!--------------------------------------------------------------------------------------------------
!> @brief Access by index and convert to bool array (1D)
!--------------------------------------------------------------------------------------------------
function tNode_get_byIndex_as1dBool(self,i) result(nodeAs1dBool)

  class(tNode), intent(in), target :: self
  integer,      intent(in)         :: i
  logical, dimension(:), allocatable :: nodeAs1dBool

  class(tNode), pointer :: node
  class(tList), pointer :: list

  node => self%get(i)
  select type(node)
    class is(tList)
      list => node%asList()
      nodeAs1dBool = list%as1dBool()
    class default
      call IO_error(706,ext_msg='Expected list of Booleans')
  endselect

end function tNode_get_byIndex_as1dBool


!--------------------------------------------------------------------------------------------------
!> @brief Access by index and convert to string array (1D)
!--------------------------------------------------------------------------------------------------
function tNode_get_byIndex_as1dString(self,i) result(nodeAs1dString)

  class(tNode), intent(in), target :: self
  integer,      intent(in)         :: i
  character(len=:), allocatable, dimension(:) :: nodeAs1dString

  class(tNode), pointer :: node
  type(tList),  pointer :: list

  node => self%get(i)
  select type(node)
    class is(tList)
      list => node%asList()
      nodeAs1dString = list%as1dString()
    class default
      call IO_error(706,ext_msg='Expected list of strings')
  endselect

end function tNode_get_byIndex_as1dString


!--------------------------------------------------------------------------------------------------
!> @brief Returns the key in a dictionary as a string
!--------------------------------------------------------------------------------------------------
function tNode_get_byIndex_asKey(self,i)  result(key)

  class(tNode),     intent(in), target  :: self
  integer,          intent(in)          :: i

  character(len=:), allocatable         :: key
  integer              :: j
  type(tDict), pointer :: dict
  type(tItem), pointer :: item

  select type(self)
    class is(tDict)
      dict => self%asDict()
      item => dict%first
      do j = 1, min(i,dict%length)-1
        item => item%next
      enddo
    class default
      call IO_error(706,ext_msg='Expected dict')
  endselect

  key = item%key

end function tNode_get_byIndex_asKey


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
  select type(self)
    class is(tDict)
      dict => self%asDict()
      do j=1, dict%length
        if (dict%getKey(j) == k) then
          exists = .true.
          return
        endif
      enddo
    class is(tList)
      list => self%asList()
      do j=1, list%length
        if (list%get_asString(j) == k) then
          exists = .true.
          return
        endif
      enddo
    class default
      call IO_error(706,ext_msg='Expected list or dict')
  endselect

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

  select type(self)
    class is(tDict)
      self_ => self%asDict()
    class default
      call IO_error(706,ext_msg='Expected dict for key '//k)
  endselect

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
    select type(node)
      class is(tScalar)
        scalar => node%asScalar()
        nodeAsFloat = scalar%asFloat()
      class default
        call IO_error(706,ext_msg='Expected scalar float for key '//k)
    endselect
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
    select type(node)
      class is(tScalar)
        scalar => node%asScalar()
        nodeAsInt = scalar%asInt()
      class default
        call IO_error(706,ext_msg='Expected scalar integer for key '//k)
    endselect
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
    select type(node)
      class is(tScalar)
        scalar => node%asScalar()
        nodeAsBool = scalar%asBool()
      class default
        call IO_error(706,ext_msg='Expected scalar Boolean for key '//k)
    endselect
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
    select type(node)
      class is(tScalar)
        scalar => node%asScalar()
        nodeAsString = scalar%asString()
      class default
        call IO_error(706,ext_msg='Expected scalar string for key '//k)
    endselect
  elseif (present(defaultVal)) then
    nodeAsString = defaultVal
  else
    call IO_error(143,ext_msg=k)
  endif

end function tNode_get_byKey_asString


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to float array (1D)
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_as1dFloat(self,k,defaultVal,requiredSize) result(nodeAs1dFloat)

  class(tNode),     intent(in), target                 :: self
  character(len=*), intent(in)                         :: k
  real(pReal),      intent(in), dimension(:), optional :: defaultVal
  integer,          intent(in),               optional :: requiredSize

  real(pReal), dimension(:), allocatable :: nodeAs1dFloat

  class(tNode), pointer :: node
  type(tList),  pointer :: list

  if (self%contains(k)) then
    node => self%get(k)
    select type(self)
      class is(tList)
        list => node%asList()
        nodeAs1dFloat = list%as1dFloat()
      class default
        call IO_error(706,ext_msg='Expected 1D float array for key '//k)
    endselect
  elseif (present(defaultVal)) then
    nodeAs1dFloat = defaultVal
  else
    call IO_error(143,ext_msg=k)
  endif

  if (present(requiredSize)) then
    if (requiredSize /= size(nodeAs1dFloat)) call IO_error(146,ext_msg=k)
  endif

end function tNode_get_byKey_as1dFloat


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to float array (2D)
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_as2dFloat(self,k,defaultVal) result(nodeAs2dFloat)

  class(tNode),     intent(in), target                   :: self
  character(len=*), intent(in)                           :: k
  real(pReal),      intent(in), dimension(:,:), optional :: defaultVal

  real(pReal), dimension(:,:), allocatable :: nodeAs2dFloat

  class(tNode), pointer :: node
  type(tList),  pointer :: rows

  if(self%contains(k)) then
    node => self%get(k)
    select type(node)
      class is(tList)
        rows => node%asList()
        nodeAs2dFloat = rows%as2dFloat()
      class default
        call IO_error(706,ext_msg='Expected 2D float array for key '//k)
    endselect
  elseif(present(defaultVal)) then
    nodeAs2dFloat = defaultVal
  else
    call IO_error(143,ext_msg=k)
  endif

end function tNode_get_byKey_as2dFloat


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to int array (1D)
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_as1dInt(self,k,defaultVal,requiredSize) result(nodeAs1dInt)

  class(tNode),          intent(in), target   :: self
  character(len=*),      intent(in)           :: k
  integer, dimension(:), intent(in), optional :: defaultVal
  integer,               intent(in), optional :: requiredSize
  integer, dimension(:), allocatable :: nodeAs1dInt

  class(tNode), pointer :: node
  type(tList),  pointer :: list

  if (self%contains(k)) then
    node => self%get(k)
    select type(node)
      class is(tList)
        list => node%asList()
        nodeAs1dInt = list%as1dInt()
      class default
        call IO_error(706,ext_msg='Expected 1D integer array for key '//k)
    endselect
  elseif (present(defaultVal)) then
    nodeAs1dInt = defaultVal
  else
    call IO_error(143,ext_msg=k)
  endif

  if (present(requiredSize)) then
    if (requiredSize /= size(nodeAs1dInt)) call IO_error(146,ext_msg=k)
  endif

end function tNode_get_byKey_as1dInt


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to bool array (1D)
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_as1dBool(self,k,defaultVal) result(nodeAs1dBool)

  class(tNode),          intent(in), target   :: self
  character(len=*),      intent(in)           :: k
  logical, dimension(:), intent(in), optional :: defaultVal
  logical, dimension(:), allocatable          :: nodeAs1dBool

  class(tNode), pointer :: node
  type(tList),  pointer :: list

  if (self%contains(k)) then
    node => self%get(k)
    select type(node)
      class is(tList)
        list => node%asList()
        nodeAs1dBool = list%as1dBool()
      class default
        call IO_error(706,ext_msg='Expected 1D Boolean array for key '//k)
    endselect
  elseif (present(defaultVal)) then
    nodeAs1dBool = defaultVal
  else
    call IO_error(143,ext_msg=k)
  endif

end function tNode_get_byKey_as1dBool


!--------------------------------------------------------------------------------------------------
!> @brief Access by key and convert to string array (1D)
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_as1dString(self,k,defaultVal) result(nodeAs1dString)

  class(tNode),     intent(in), target                 :: self
  character(len=*), intent(in)                         :: k
  character(len=*), intent(in), dimension(:), optional :: defaultVal
  character(len=:), allocatable, dimension(:)          :: nodeAs1dString

  class(tNode), pointer :: node
  type(tList),  pointer :: list

  if (self%contains(k)) then
    node => self%get(k)
    select type(node)
      class is(tList)
        list => node%asList()
        nodeAs1dString = list%as1dString()
      class default
        call IO_error(706,ext_msg='Expected 1D string array for key '//k)
    endselect
  elseif (present(defaultVal)) then
    nodeAs1dString = defaultVal
  else
    call IO_error(143,ext_msg=k)
  endif

end function tNode_get_byKey_as1dString


!--------------------------------------------------------------------------------------------------
!> @brief Returns string output array (1D) (hack for GNU)
!--------------------------------------------------------------------------------------------------
function output_as1dString(self)  result(output)                   !ToDo: SR: Remove whenever GNU works

  class(tNode), pointer,intent(in)   ::  self
  character(len=pStringLen), allocatable, dimension(:) :: output

  class(tNode), pointer :: output_list
  integer :: o

  output_list   => self%get('output',defaultVal=emptyList)
  allocate(output(output_list%length))
  do o = 1, output_list%length
    output(o) = output_list%get_asString(o)
  enddo

end function output_as1dString


!--------------------------------------------------------------------------------------------------
!> @brief Returns the index of a key in a dictionary
!--------------------------------------------------------------------------------------------------
function tNode_get_byKey_asIndex(self,key)  result(keyIndex)

  class(tNode),     intent(in), target  :: self
  character(len=*), intent(in)          :: key

  integer              :: keyIndex
  type(tDict), pointer :: dict
  type(tItem), pointer :: item

  dict => self%asDict()
  item => dict%first
  keyIndex = 1
  do while (associated(item%next) .and. item%key /= key)
    item => item%next
    keyIndex = keyIndex+1
  enddo

  if (item%key /= key) call IO_error(140,ext_msg=key)

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
!> @brief Convert to float array (1D)
!--------------------------------------------------------------------------------------------------
function tList_as1dFloat(self)

  class(tList), intent(in), target :: self
  real(pReal), dimension(:), allocatable :: tList_as1dFloat

  integer :: i
  type(tItem),   pointer :: item
  type(tScalar), pointer :: scalar

  allocate(tList_as1dFloat(self%length))
  item => self%first
  do i = 1, self%length
    scalar => item%node%asScalar()
    tList_as1dFloat(i) = scalar%asFloat()
    item => item%next
  enddo

end function tList_as1dFloat


!--------------------------------------------------------------------------------------------------
!> @brief Convert to float array (2D)
!--------------------------------------------------------------------------------------------------
function tList_as2dFloat(self)

  class(tList), intent(in), target :: self
  real(pReal), dimension(:,:), allocatable :: tList_as2dFloat

  integer :: i
  class(tNode), pointer :: row
  type(tList),  pointer :: row_data

  row => self%get(1)
  row_data => row%asList()
  allocate(tList_as2dFloat(self%length,row_data%length))

  do i=1,self%length
    row => self%get(i)
    row_data => row%asList()
    if(row_data%length /= size(tList_as2dFloat,2)) call IO_error(709,ext_msg='Varying number of columns')
    tList_as2dFloat(i,:) = self%get_as1dFloat(i)
  enddo

end function tList_as2dFloat


!--------------------------------------------------------------------------------------------------
!> @brief Convert to int array (1D)
!--------------------------------------------------------------------------------------------------
function tList_as1dInt(self)

  class(tList), intent(in), target :: self
  integer, dimension(:), allocatable :: tList_as1dInt

  integer :: i
  type(tItem),   pointer :: item
  type(tScalar), pointer :: scalar

  allocate(tList_as1dInt(self%length))
  item => self%first
  do i = 1, self%length
    scalar => item%node%asScalar()
    tList_as1dInt(i) = scalar%asInt()
    item => item%next
  enddo

end function tList_as1dInt


!--------------------------------------------------------------------------------------------------
!> @brief Convert to bool array (1D)
!--------------------------------------------------------------------------------------------------
function tList_as1dBool(self)

  class(tList), intent(in), target :: self
  logical, dimension(:), allocatable :: tList_as1dBool

  integer :: i
  type(tItem),   pointer :: item
  type(tScalar), pointer :: scalar

  allocate(tList_as1dBool(self%length))
  item => self%first
  do i = 1, self%length
    scalar => item%node%asScalar()
    tList_as1dBool(i) = scalar%asBool()
    item => item%next
  enddo

end function tList_as1dBool


!--------------------------------------------------------------------------------------------------
!> @brief Convert to string array (1D)
!--------------------------------------------------------------------------------------------------
function tList_as1dString(self)

  class(tList), intent(in), target :: self
  character(len=:), allocatable, dimension(:) :: tList_as1dString

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

  allocate(character(len=len_max) :: tList_as1dString(self%length))
  item => self%first
  do i = 1, self%length
    scalar => item%node%asScalar()
    tList_as1dString(i) = scalar%asString()
    item => item%next
  enddo

end function tList_as1dString


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
    if (item%key /= key) then
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
