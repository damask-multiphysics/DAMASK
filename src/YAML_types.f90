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

  implicit none(type,external)
  private

  type, abstract, public :: tNode
    integer :: &
      length = 0
    contains
    procedure(asFormattedString), deferred :: &
      asFormattedString
    procedure :: &
      asScalar => tNode_asScalar, &
      asList   => tNode_asList, &
      asDict   => tNode_asDict
  end type tNode

  type, extends(tNode), public :: tScalar
    character(len=:), allocatable, private :: &
      value
    contains
    procedure :: &
      asFormattedString => tScalar_asFormattedString, &
      asFloat   => tScalar_asFloat, &
      asInt     => tScalar_asInt, &
      asBool    => tScalar_asBool, &
      asString  => tScalar_asString
  end type tScalar

  type, extends(tNode), public :: tList
    class(tItem), pointer :: &
      first => NULL(), &
      last => NULL()
    contains
    procedure :: &
      asFormattedString => tList_asFormattedString, &
      append            => tList_append, &
      as1dFloat         => tList_as1dFloat, &
      as2dFloat         => tList_as2dFloat, &
      as1dInt           => tList_as1dInt, &
      as1dBool          => tList_as1dBool, &
      as1dString        => tList_as1dString, &
      contains          => tList_contains, &
      tList_get, &
      tList_get_scalar, &
      tList_get_list, &
      tList_get_dict, &
      tList_get_asFloat, &
      tList_get_as1dFloat, &
      tList_get_asInt, &
      tList_get_as1dInt, &
      tList_get_asBool, &
      tList_get_as1dBool, &
      tList_get_asString, &
      tList_get_as1dString
    generic :: get            => tList_get
    generic :: get_scalar     => tList_get_scalar
    generic :: get_list       => tList_get_list
    generic :: get_dict       => tList_get_dict
    generic :: get_asFloat    => tList_get_asFloat
    generic :: get_as1dFloat  => tList_get_as1dFloat
    generic :: get_asInt      => tList_get_asInt
    generic :: get_as1dInt    => tList_get_as1dInt
    generic :: get_asBool     => tList_get_asBool
    generic :: get_as1dBool   => tList_get_as1dBool
    generic :: get_asString   => tList_get_asString
    generic :: get_as1dString => tList_get_as1dString
    final :: tList_finalize
  end type tList

  type, extends(tList), public :: tDict
    contains
    procedure ::  &
      asFormattedString      => tDict_asFormattedString, &
      set                    => tDict_set, &
      index                  => tDict_index, &
      key                    => tDict_key, &
      keys                   => tDict_keys, &
      contains               => tDict_contains, &
      tDict_get, &
      tDict_get_scalar, &
      tDict_get_list, &
      tDict_get_dict, &
      tDict_get_asFloat, &
      tDict_get_as1dFloat, &
      tDict_get_as2dFloat, &
      tDict_get_asInt, &
      tDict_get_as1dInt, &
      tDict_get_asBool, &
      tDict_get_as1dBool, &
      tDict_get_asString, &
      tDict_get_as1dString
    generic :: get            => tDict_get
    generic :: get_scalar     => tDict_get_scalar
    generic :: get_list       => tDict_get_list
    generic :: get_dict       => tDict_get_dict
    generic :: get_asFloat    => tDict_get_asFloat
    generic :: get_as1dFloat  => tDict_get_as1dFloat
    generic :: get_as2dFloat  => tDict_get_as2dFloat
    generic :: get_asInt      => tDict_get_asInt
    generic :: get_as1dInt    => tDict_get_as1dInt
    generic :: get_asBool     => tDict_get_asBool
    generic :: get_as1dBool   => tDict_get_as1dBool
    generic :: get_asString   => tDict_get_asString
    generic :: get_as1dString => tDict_get_as1dString
   end type tDict


  type, public :: tItem
    character(len=:), allocatable :: key
    class(tNode),     pointer     :: node => NULL()
    class(tItem),     pointer     :: next => NULL()
    contains
    final :: tItem_finalize
  end type tItem

  type(tDict), target, public :: &
    emptyDict
  type(tList), target, public :: &
    emptyList

  abstract interface

    recursive function asFormattedString(self)
      import tNode
      character(len=:), allocatable      :: asFormattedString
      class(tNode), intent(in), target   :: self
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
#ifdef __GFORTRAN__
    output_as1dString, &                                       !ToDo: Hack for GNU. Remove later
#endif
    assignment(=)

contains

!--------------------------------------------------------------------------------------------------
!> @brief Do sanity checks.
!--------------------------------------------------------------------------------------------------
subroutine YAML_types_init

  print'(/,1x,a)', '<<<+-  YAML_types init  -+>>>'

  call selfTest()

end subroutine YAML_types_init


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some type bound procedures.
!--------------------------------------------------------------------------------------------------
subroutine selfTest()

  scalar: block
    type(tScalar), target :: s
    type(tScalar), pointer :: s_pointer


    s_pointer => s%asScalar()
    s = '1'
    if (s%asInt() /= 1)                   error stop 'tScalar_asInt'
    if (s_pointer%asInt() /= 1)           error stop 'tScalar_asInt(pointer)'
    if (dNeq(s%asFloat(),1.0_pReal))      error stop 'tScalar_asFloat'
    s = 'true'
    if (.not. s%asBool())                 error stop 'tScalar_asBool'
    if (.not. s_pointer%asBool())         error stop 'tScalar_asBool(pointer)'
    if (s%asString() /= 'true')           error stop 'tScalar_asString'
    if (s%asFormattedString() /= 'true')  error stop 'tScalar_asFormattedString'


  end block scalar

  list: block
    type(tList),   pointer :: l, l_pointer
    type(tScalar), pointer :: s1,s2


    allocate(s1)
    allocate(s2)
    s1 = '1'
    s2 = '2'
    allocate(l)
    l_pointer => l%asList()
    if (l%contains('1'))                                 error stop 'empty tList_contains'
    if (l_pointer%contains('1'))                         error stop 'empty tList_contains(pointer)'
    call l%append(s1)
    call l%append(s2)
    if (l%length /= 2)                                   error stop 'tList%len'
    if (dNeq(l%get_asFloat(1),1.0_pReal))                error stop 'tList_get_asFloat'
    if (l%get_asInt(1) /= 1)                             error stop 'tList_get_asInt'
    if (l%get_asString(2) /= '2')                        error stop 'tList_get_asString'
    if (any(l%as1dInt() /= [1,2]))                       error stop 'tList_as1dInt'
    if (any(dNeq(l%as1dFloat(),real([1.0,2.0],pReal))))  error stop 'tList_as1dFloat'
    s1 = 'true'
    s2 = 'false'
    if (any(l%as1dBool() .neqv. [.true.,.false.]))       error stop 'tList_as1dBool'
    if (any(l%as1dString() /= ['true ','false']))        error stop 'tList_as1dString'
    if (l%asFormattedString() /= '[true, false]')        error stop 'tList_asFormattedString'
    if (   .not. l%contains('true') &
      .or. .not. l%contains('false'))                    error stop 'tList_contains'

  end block list

  dict: block
    type(tDict),   pointer :: d, d_pointer
    type(tList),   pointer :: l
    type(tScalar), pointer :: s1,s2,s3,s4


    allocate(s1)
    allocate(s2)
    s1 = '1'
    s2 = '2'
    allocate(l)
    call l%append(s1)
    call l%append(s2)

    allocate(s3)
    allocate(s4)
    s3 = '3'
    s4 = '4'
    allocate(d)
    d_pointer => d%asDict()
    if (d%contains('one-two'))                           error stop 'empty tDict_contains'
    if (d_pointer%contains('one-two'))                   error stop 'empty tDict_contains(pointer)'
    if (d%get_asInt('one-two',defaultVal=-1) /= -1)      error stop 'empty tDict_get'
    call d%set('one-two',l)
    call d%set('three',s3)
    call d%set('four',s4)
    if (d%asFormattedString() /= '{one-two: [1, 2], three: 3, four: 4}') &
                                                         error stop 'tDict_asFormattedString'
    if (d%get_asInt('three') /= 3)                       error stop 'tDict_get_asInt'
    if (dNeq(d%get_asFloat('three'),3.0_pReal))          error stop 'tDict_get_asFloat'
    if (d%get_asString('three') /= '3')                  error stop 'tDict_get_asString'
    if (any(d%get_as1dInt('one-two') /= [1,2]))          error stop 'tDict_get_as1dInt'
    call d%set('one-two',s4)
    if (d%asFormattedString() /= '{one-two: 4, three: 3, four: 4}') &
                                                         error stop 'tDict_set overwrite'
    if (   .not. d%contains('one-two') &
      .or. .not. d%contains('three') &
      .or. .not. d%contains('four') &
      )                                                  error stop 'tDict_contains'

  end block dict

end subroutine selfTest


!---------------------------------------------------------------------------------------------------
!> @brief Init from string.
!---------------------------------------------------------------------------------------------------
type(tScalar) pure function tScalar_init__(value)

  character(len=*), intent(in) ::  value


  tScalar_init__%value = value

end function tScalar_init__


!---------------------------------------------------------------------------------------------------
!> @brief Set value from string.
!---------------------------------------------------------------------------------------------------
elemental pure subroutine tScalar_assign__(self,value)

  type(tScalar),    intent(out) :: self
  character(len=*), intent(in)  :: value


  self%value = value

end subroutine tScalar_assign__


!--------------------------------------------------------------------------------------------------
!> @brief Format as string (YAML flow style).
!--------------------------------------------------------------------------------------------------
recursive function tScalar_asFormattedString(self) result(str)

  class (tScalar), intent(in), target    :: self
  character(len=:), allocatable          :: str


  str = trim(self%value)

end function tScalar_asFormattedString


!--------------------------------------------------------------------------------------------------
!> @brief Type guard, guarantee scalar.
!--------------------------------------------------------------------------------------------------
function tNode_asScalar(self) result(scalar)

  class(tNode),   intent(in), target :: self
  class(tScalar), pointer            :: scalar


  select type(self)
    class is(tScalar)
      scalar => self
    class default
      nullify(scalar)
      call IO_error(706,'"'//trim(self%asFormattedString())//'" is not a scalar')
  end select

end function tNode_asScalar


!--------------------------------------------------------------------------------------------------
!> @brief Type guard, guarantee list.
!--------------------------------------------------------------------------------------------------
function tNode_asList(self) result(list)

  class(tNode), intent(in), target :: self
  class(tList), pointer            :: list


  select type(self)
    class is(tList)
      list => self
    class default
      nullify(list)
      call IO_error(706,'"'//trim(self%asFormattedString())//'" is not a list')
  end select

end function tNode_asList


!--------------------------------------------------------------------------------------------------
!> @brief Type guard, guarantee dict.
!--------------------------------------------------------------------------------------------------
function tNode_asDict(self) result(dict)

  class(tNode), intent(in), target :: self
  class(tDict), pointer            :: dict


  select type(self)
    class is(tDict)
      dict => self
    class default
      nullify(dict)
      call IO_error(706,'"'//trim(self%asFormattedString())//'" is not a dict')
  end select

end function tNode_asDict


!--------------------------------------------------------------------------------------------------
!> @brief Convert to float.
!--------------------------------------------------------------------------------------------------
function tScalar_asFloat(self)

  class(tScalar), intent(in), target :: self
  real(pReal) :: tScalar_asFloat


  tScalar_asFloat = IO_stringAsFloat(self%value)

end function tScalar_asFloat


!--------------------------------------------------------------------------------------------------
!> @brief Convert to int.
!--------------------------------------------------------------------------------------------------
function tScalar_asInt(self)

  class(tScalar), intent(in), target :: self
  integer :: tScalar_asInt


  tScalar_asInt = IO_stringAsInt(self%value)

end function tScalar_asInt


!--------------------------------------------------------------------------------------------------
!> @brief Convert to bool.
!--------------------------------------------------------------------------------------------------
function tScalar_asBool(self)

  class(tScalar), intent(in), target :: self
  logical :: tScalar_asBool


  tScalar_asBool = IO_stringAsBool(self%value)

end function tScalar_asBool


!--------------------------------------------------------------------------------------------------
!> @brief Convert to string.
!--------------------------------------------------------------------------------------------------
function tScalar_asString(self)

  class(tScalar), intent(in), target :: self
  character(len=:), allocatable :: tScalar_asString


  tScalar_asString = self%value

end function tScalar_asString


!--------------------------------------------------------------------------------------------------
!> @brief Format as string (YAML flow style).
!--------------------------------------------------------------------------------------------------
recursive function tList_asFormattedString(self) result(str)

  class(tList),intent(in),target      :: self

  type(tItem), pointer  :: item
  character(len=:), allocatable :: str
  integer :: i

  str = '['
  item => self%first
  do i = 2, self%length
    str = str//item%node%asFormattedString()//', '
    item => item%next
  end do
  str = str//item%node%asFormattedString()//']'

end function tList_asFormattedString


!--------------------------------------------------------------------------------------------------
!> @brief Append element.
!--------------------------------------------------------------------------------------------------
subroutine tList_append(self,node)

  class(tList), intent(inout)         :: self
  class(tNode), intent(in), target    :: node

  type(tItem), pointer :: item


  if (.not. associated(self%first)) then
    allocate(item)
    self%first => item
    self%last => item
  else
    allocate(self%last%next)
    item => self%last%next
    self%last => item
  end if

  item%node => node
  self%length = self%length + 1

end subroutine tList_append


!--------------------------------------------------------------------------------------------------
!> @brief Convert to float array (1D).
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
  end do

end function tList_as1dFloat


!--------------------------------------------------------------------------------------------------
!> @brief Convert to float array (2D).
!--------------------------------------------------------------------------------------------------
function tList_as2dFloat(self)

  class(tList), intent(in), target :: self
  real(pReal), dimension(:,:), allocatable :: tList_as2dFloat

  integer :: i
  type(tList), pointer :: row_data


  row_data => self%get_list(1)
  allocate(tList_as2dFloat(self%length,row_data%length))

  do i = 1, self%length
    row_data => self%get_list(i)
    if (row_data%length /= size(tList_as2dFloat,2)) call IO_error(709,ext_msg='inconsistent column count in tList_as2dFloat')
    tList_as2dFloat(i,:) = self%get_as1dFloat(i)
  end do

end function tList_as2dFloat


!--------------------------------------------------------------------------------------------------
!> @brief Convert to int array (1D).
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
  end do

end function tList_as1dInt


!--------------------------------------------------------------------------------------------------
!> @brief Convert to bool array (1D).
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
  end do

end function tList_as1dBool


!--------------------------------------------------------------------------------------------------
!> @brief Convert to string array (1D).
!--------------------------------------------------------------------------------------------------
function tList_as1dString(self)

  class(tList), intent(in), target :: self
#ifdef __GFORTRAN__
  character(len=pStringLen), allocatable, dimension(:) :: tList_as1dString
#else
  character(len=:), allocatable, dimension(:) :: tList_as1dString
#endif

  integer :: j
  type(tItem),   pointer :: item
  type(tScalar), pointer :: scalar


#ifdef __GFORTRAN__
  allocate(tList_as1dString(self%length))
#else
  integer :: len_max
  len_max = 0
  item => self%first
  do j = 1, self%length
    scalar => item%node%asScalar()
    len_max = max(len_max, len_trim(scalar%asString()))
    item => item%next
  end do

  allocate(character(len=len_max) :: tList_as1dString(self%length))
#endif
  item => self%first
  do j = 1, self%length
    scalar => item%node%asScalar()
    tList_as1dString(j) = scalar%asString()
    item => item%next
  end do

end function tList_as1dString


!-------------------------------------------------------------------------------------------------
!> @brief Check for existence of (string) value.
!-------------------------------------------------------------------------------------------------
function tList_contains(self,k)  result(exists)

  class(tList),     intent(in), target  :: self
  character(len=*), intent(in)          :: k
  logical                               :: exists

  integer :: j
  type(tItem),   pointer :: item
  type(tScalar), pointer :: scalar


  item => self%first
  exists = .false.
  j = 1
  do while (j <= self%length .and. .not. exists)
    scalar => item%node%asScalar()
    exists = scalar%value == k
    item => item%next
    j = j + 1
  end do

end function tList_contains


!--------------------------------------------------------------------------------------------------
!> @brief Get by index.
!--------------------------------------------------------------------------------------------------
function tList_get(self,i) result(node)

  class(tList), intent(in), target :: self
  integer,      intent(in)         :: i
  class(tNode), pointer :: node

  class(tItem), pointer :: item
  integer :: j


  if (i < 1 .or. i > self%length) call IO_error(150,ext_msg='tList_get @ '//IO_intAsString(i) &
                                                                  //' of '//IO_intAsString(self%length) )
  item => self%first
  do j = 2, i
    item => item%next
  end do
  node => item%node

end function tList_get


!--------------------------------------------------------------------------------------------------
!> @brief Get scalar by index.
!--------------------------------------------------------------------------------------------------
function tList_get_scalar(self,i) result(nodeAsScalar)

  class(tList), intent(in) :: self
  integer,      intent(in) :: i
  type(tScalar), pointer :: nodeAsScalar

  class(tNode),  pointer :: node


  node => self%get(i)
  nodeAsScalar => node%asScalar()

end function tList_get_scalar


!--------------------------------------------------------------------------------------------------
!> @brief Get list by index.
!--------------------------------------------------------------------------------------------------
function tList_get_list(self,i) result(nodeAsList)

  class(tList), intent(in) :: self
  integer,      intent(in) :: i
  type(tList),  pointer :: nodeAsList

  class(tNode),  pointer :: node


  node => self%get(i)
  nodeAsList => node%asList()

end function tList_get_list


!--------------------------------------------------------------------------------------------------
!> @brief Get dict by index.
!--------------------------------------------------------------------------------------------------
function tList_get_dict(self,i) result(nodeAsDict)

  class(tList), intent(in) :: self
  integer,      intent(in) :: i
  type(tDict),  pointer :: nodeAsDict

  class(tNode),  pointer :: node


  node => self%get(i)
  nodeAsDict => node%asDict()

end function tList_get_dict


!--------------------------------------------------------------------------------------------------
!> @brief Get scalar by index and convert to float.
!--------------------------------------------------------------------------------------------------
function tList_get_asFloat(self,i) result(nodeAsFloat)

  class(tList), intent(in) :: self
  integer,      intent(in) :: i
  real(pReal) :: nodeAsFloat

  class(tScalar),  pointer :: scalar


  scalar => self%get_scalar(i)
  nodeAsFloat = scalar%asFloat()

end function tList_get_asFloat


!--------------------------------------------------------------------------------------------------
!> @brief Get list by index and convert to float array (1D).
!--------------------------------------------------------------------------------------------------
function tList_get_as1dFloat(self,i) result(nodeAs1dFloat)

  class(tList), intent(in) :: self
  integer,      intent(in) :: i
  real(pReal), dimension(:), allocatable :: nodeAs1dFloat

  class(tList),  pointer :: list


  list => self%get_list(i)
  nodeAs1dFloat = list%as1dFloat()

end function tList_get_as1dFloat


!--------------------------------------------------------------------------------------------------
!> @brief Get scalar by index and convert to int.
!--------------------------------------------------------------------------------------------------
function tList_get_asInt(self,i) result(nodeAsInt)

  class(tList), intent(in) :: self
  integer,      intent(in) :: i
  integer :: nodeAsInt

  class(tScalar),  pointer :: scalar


  scalar => self%get_scalar(i)
  nodeAsInt = scalar%asInt()

end function tList_get_asInt


!--------------------------------------------------------------------------------------------------
!> @brief Get list by index and convert to int array (1D).
!--------------------------------------------------------------------------------------------------
function tList_get_as1dInt(self,i) result(nodeAs1dInt)

  class(tList), intent(in) :: self
  integer,      intent(in) :: i
  integer, dimension(:), allocatable :: nodeAs1dInt

  class(tList),  pointer :: list


  list => self%get_list(i)
  nodeAs1dInt = list%as1dInt()

end function tList_get_as1dInt


!--------------------------------------------------------------------------------------------------
!> @brief Get scalar by index and convert to bool
!--------------------------------------------------------------------------------------------------
function tList_get_asBool(self,i) result(nodeAsBool)

  class(tList), intent(in) :: self
  integer,      intent(in) :: i
  logical :: nodeAsBool

  class(tScalar),  pointer :: scalar


  scalar => self%get_scalar(i)
  nodeAsBool = scalar%asBool()

end function tList_get_asBool


!--------------------------------------------------------------------------------------------------
!> @brief Get list by index and convert to bool array (1D).
!--------------------------------------------------------------------------------------------------
function tList_get_as1dBool(self,i) result(nodeAs1dBool)

  class(tList), intent(in) :: self
  integer,      intent(in) :: i
  logical, dimension(:), allocatable :: nodeAs1dBool

  class(tList),  pointer :: list


  list => self%get_list(i)
  nodeAs1dBool = list%as1dBool()

end function tList_get_as1dBool


!--------------------------------------------------------------------------------------------------
!> @brief Get scalar by index and convert to string.
!--------------------------------------------------------------------------------------------------
function tList_get_asString(self,i) result(nodeAsString)

  class(tList), intent(in) :: self
  integer,      intent(in) :: i
  character(len=:), allocatable :: nodeAsString

  class(tScalar),  pointer :: scalar


  scalar => self%get_scalar(i)
  nodeAsString = scalar%asString()

end function tList_get_asString


!--------------------------------------------------------------------------------------------------
!> @brief Get list by index and convert to string array (1D).
!--------------------------------------------------------------------------------------------------
function tList_get_as1dString(self,i) result(nodeAs1dString)

  class(tList), intent(in) :: self
  integer,      intent(in) :: i
  character(len=:), allocatable, dimension(:) :: nodeAs1dString

  type(tList),  pointer :: list


  list => self%get_list(i)
  nodeAs1dString = list%as1dString()

end function tList_get_as1dString


!--------------------------------------------------------------------------------------------------
!> @brief Free associated memory.
!--------------------------------------------------------------------------------------------------
recursive subroutine tList_finalize(self)

  type (tList),intent(inout) :: self

  deallocate(self%first)

end subroutine tList_finalize


!--------------------------------------------------------------------------------------------------
!> @brief Format as string (YAML flow style).
!--------------------------------------------------------------------------------------------------
recursive function tDict_asFormattedString(self) result(str)

  class(tDict),intent(in),target      :: self

  type(tItem),pointer  :: item
  character(len=:), allocatable :: str
  integer :: i


  str = '{'
  item => self%first
  do i = 2, self%length
    str = str//trim(item%key)//': '//item%node%asFormattedString()//', '
    item => item%next
  end do
  str = str//trim(item%key)//': '//item%node%asFormattedString()//'}'

end function tDict_asFormattedString


!--------------------------------------------------------------------------------------------------
!> @brief Set value (either replace or add new).
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
    searchExisting: do while (associated(item%next) .and. item%key /= key)
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
!> @brief Return the index of a key.
!--------------------------------------------------------------------------------------------------
function tDict_index(self,key)  result(keyIndex)

  class(tDict),     intent(in), target  :: self
  character(len=*), intent(in)          :: key

  integer              :: keyIndex
  type(tItem), pointer :: item


  item => self%first
  keyIndex = 1
  do while (associated(item%next) .and. item%key /= key)
    item => item%next
    keyIndex = keyIndex+1
  end do

  if (item%key /= key) call IO_error(140,ext_msg=key)

end function tDict_index


!--------------------------------------------------------------------------------------------------
!> @brief Get key of given index.
!--------------------------------------------------------------------------------------------------
function tDict_key(self,i)  result(key)

  class(tDict),     intent(in), target  :: self
  integer,          intent(in)          :: i

  character(len=:), allocatable         :: key
  integer              :: j
  type(tItem), pointer :: item


  if (i < 1 .or. i > self%length) call IO_error(150,ext_msg='tDict_key @ '//IO_intAsString(i) &
                                                                  //' of '//IO_intAsString(self%length) )
  item => self%first
  do j = 2, i
    item => item%next
  end do

  key = item%key

end function tDict_key


!--------------------------------------------------------------------------------------------------
!> @brief Get all keys.
!--------------------------------------------------------------------------------------------------
function tDict_keys(self) result(keys)

  class(tDict), intent(in) :: self
  character(len=:), dimension(:), allocatable :: keys

  character(len=pStringLen), dimension(:), allocatable :: temp
  integer :: j, l


  allocate(temp(self%length))
  l = 0
  do j = 1, self%length
    temp(j) = self%key(j)
    l = max(len_trim(temp(j)),l)
  end do

  allocate(character(l)::keys(self%length))
  do j = 1, self%length
    keys(j) = trim(temp(j))
  end do

end function tDict_keys


!-------------------------------------------------------------------------------------------------
!> @brief Check whether a given key is present.
!-------------------------------------------------------------------------------------------------
function tDict_contains(self,k)  result(exists)

  class(tDict),     intent(in), target  :: self
  character(len=*), intent(in)          :: k
  logical                               :: exists

  integer   :: j


  exists = .false.
  j = 1
  do while(j <= self%length .and. .not. exists)
    exists = self%key(j) == k
    j = j + 1
  end do

end function tDict_contains


!--------------------------------------------------------------------------------------------------
!> @brief Get by key.
!--------------------------------------------------------------------------------------------------
function tDict_get(self,k,defaultVal) result(node)

  class(tDict),     intent(in), target         :: self
  character(len=*), intent(in)                 :: k
  class(tNode),     intent(in),optional,target :: defaultVal
  class(tNode),     pointer :: node

  type(tItem), pointer :: item
  integer :: j

  item => self%first

  do j=1, self%length
    if (item%key == k) then
      node => item%node
      return
    end if
    item => item%next
  end do

  if (present(defaultVal)) then
    node => defaultVal
  else
    call IO_error(143,ext_msg=k)
  end if

end function tDict_get


!--------------------------------------------------------------------------------------------------
!> @brief Get scalar by key.
!--------------------------------------------------------------------------------------------------
function tDict_get_scalar(self,k,defaultVal) result(nodeAsScalar)

  class(tDict),     intent(in) :: self
  character(len=*), intent(in) :: k
  type(tScalar),    intent(in), optional, target :: defaultVal
  type(tScalar), pointer :: nodeAsScalar

  class(tNode),  pointer :: node


  node => self%get(k,defaultVal)
  nodeAsScalar => node%asScalar()

end function tDict_get_scalar


!--------------------------------------------------------------------------------------------------
!> @brief Get list by key.
!--------------------------------------------------------------------------------------------------
function tDict_get_list(self,k,defaultVal) result(nodeAsList)

  class(tDict),     intent(in) :: self
  character(len=*), intent(in) :: k
  type(tList),      intent(in), optional, target :: defaultVal
  type(tList),  pointer :: nodeAsList

  class(tNode), pointer :: node


  node => self%get(k,defaultVal)
  nodeAsList => node%asList()

end function tDict_get_list


!--------------------------------------------------------------------------------------------------
!> @brief Get dict by key.
!--------------------------------------------------------------------------------------------------
function tDict_get_dict(self,k,defaultVal) result(nodeAsDict)

  class(tDict),     intent(in) :: self
  character(len=*), intent(in) :: k
  type(tDict),      intent(in), optional, target :: defaultVal
  type(tDict),  pointer :: nodeAsDict

  class(tNode), pointer :: node


  node => self%get(k,defaultVal)
  nodeAsDict => node%asDict()

end function tDict_get_dict


!--------------------------------------------------------------------------------------------------
!> @brief Get scalar by key and convert to float.
!--------------------------------------------------------------------------------------------------
function tDict_get_asFloat(self,k,defaultVal) result(nodeAsFloat)

  class(tDict),     intent(in) :: self
  character(len=*), intent(in) :: k
  real(pReal),      intent(in), optional :: defaultVal
  real(pReal) :: nodeAsFloat

  type(tScalar), pointer :: scalar


  if (self%contains(k)) then
    scalar => self%get_scalar(k)
    nodeAsFloat = scalar%asFloat()
  elseif (present(defaultVal)) then
    nodeAsFloat = defaultVal
  else
    call IO_error(143,ext_msg=k)
  end if

end function tDict_get_asFloat


!--------------------------------------------------------------------------------------------------
!> @brief Get list by key and convert to float array (1D).
!--------------------------------------------------------------------------------------------------
function tDict_get_as1dFloat(self,k,defaultVal,requiredSize) result(nodeAs1dFloat)

  class(tDict),     intent(in) :: self
  character(len=*), intent(in) :: k
  real(pReal),      intent(in), dimension(:), optional :: defaultVal
  integer,          intent(in),               optional :: requiredSize
  real(pReal), dimension(:), allocatable :: nodeAs1dFloat

  type(tList), pointer :: list


  if (self%contains(k)) then
    list => self%get_list(k)
    nodeAs1dFloat = list%as1dFloat()
  elseif (present(defaultVal)) then
    nodeAs1dFloat = defaultVal
  else
    call IO_error(143,ext_msg=k)
  end if

  if (present(requiredSize)) then
    if (requiredSize /= size(nodeAs1dFloat)) call IO_error(146,ext_msg=k)
  end if

end function tDict_get_as1dFloat


!--------------------------------------------------------------------------------------------------
!> @brief Get list of lists by key and convert to float array (2D).
!--------------------------------------------------------------------------------------------------
function tDict_get_as2dFloat(self,k,defaultVal,requiredShape) result(nodeAs2dFloat)

  class(tDict),     intent(in) :: self
  character(len=*), intent(in) :: k
  real(pReal),      intent(in), dimension(:,:), optional :: defaultVal
  integer,          intent(in), dimension(2),   optional :: requiredShape
  real(pReal), dimension(:,:), allocatable :: nodeAs2dFloat

  type(tList), pointer :: list


  if (self%contains(k)) then
    list => self%get_list(k)
    nodeAs2dFloat = list%as2dFloat()
  elseif (present(defaultVal)) then
    nodeAs2dFloat = defaultVal
  else
    call IO_error(143,ext_msg=k)
  end if

  if (present(requiredShape)) then
    if (any(requiredShape /= shape(nodeAs2dFloat))) call IO_error(146,ext_msg=k)
  end if

end function tDict_get_as2dFloat


!--------------------------------------------------------------------------------------------------
!> @brief Get scalar by key and convert to int.
!--------------------------------------------------------------------------------------------------
function tDict_get_asInt(self,k,defaultVal) result(nodeAsInt)

  class(tDict),     intent(in) :: self
  character(len=*), intent(in) :: k
  integer,          intent(in), optional :: defaultVal
  integer :: nodeAsInt

  type(tScalar), pointer :: scalar


  if (self%contains(k)) then
    scalar => self%get_scalar(k)
    nodeAsInt = scalar%asInt()
  elseif (present(defaultVal)) then
    nodeAsInt = defaultVal
  else
    call IO_error(143,ext_msg=k)
  end if

end function tDict_get_asInt


!--------------------------------------------------------------------------------------------------
!> @brief Get list by key and convert to int array (1D).
!--------------------------------------------------------------------------------------------------
function tDict_get_as1dInt(self,k,defaultVal,requiredSize) result(nodeAs1dInt)

  class(tDict),          intent(in) :: self
  character(len=*),      intent(in) :: k
  integer, dimension(:), intent(in), optional :: defaultVal
  integer,               intent(in), optional :: requiredSize
  integer, dimension(:), allocatable :: nodeAs1dInt

  type(tList), pointer :: list


  if (self%contains(k)) then
    list => self%get_list(k)
    nodeAs1dInt = list%as1dInt()
  elseif (present(defaultVal)) then
    nodeAs1dInt = defaultVal
  else
    call IO_error(143,ext_msg=k)
  end if

  if (present(requiredSize)) then
    if (requiredSize /= size(nodeAs1dInt)) call IO_error(146,ext_msg=k)
  end if

end function tDict_get_as1dInt


!--------------------------------------------------------------------------------------------------
!> @brief Get scalar by key and convert to bool.
!--------------------------------------------------------------------------------------------------
function tDict_get_asBool(self,k,defaultVal) result(nodeAsBool)

  class(tDict),     intent(in) :: self
  character(len=*), intent(in) :: k
  logical,          intent(in), optional :: defaultVal
  logical :: nodeAsBool

  type(tScalar), pointer :: scalar


  if (self%contains(k)) then
    scalar => self%get_scalar(k)
    nodeAsBool = scalar%asBool()
  elseif (present(defaultVal)) then
    nodeAsBool = defaultVal
  else
    call IO_error(143,ext_msg=k)
  end if

end function tDict_get_asBool


!--------------------------------------------------------------------------------------------------
!> @brief Get list by key and convert to bool array (1D).
!--------------------------------------------------------------------------------------------------
function tDict_get_as1dBool(self,k,defaultVal) result(nodeAs1dBool)

  class(tDict),          intent(in) :: self
  character(len=*),      intent(in) :: k
  logical, dimension(:), intent(in), optional :: defaultVal
  logical, dimension(:), allocatable          :: nodeAs1dBool

  type(tList), pointer :: list


  if (self%contains(k)) then
    list => self%get_list(k)
    nodeAs1dBool = list%as1dBool()
  elseif (present(defaultVal)) then
    nodeAs1dBool = defaultVal
  else
    call IO_error(143,ext_msg=k)
  end if

end function tDict_get_as1dBool


!--------------------------------------------------------------------------------------------------
!> @brief Get scalar by key and convert to string.
!--------------------------------------------------------------------------------------------------
function tDict_get_asString(self,k,defaultVal) result(nodeAsString)

  class(tDict),     intent(in) :: self
  character(len=*), intent(in) :: k
  character(len=*), intent(in), optional :: defaultVal
  character(len=:), allocatable :: nodeAsString

  type(tScalar), pointer :: scalar


  if (self%contains(k)) then
    scalar => self%get_scalar(k)
    nodeAsString = scalar%asString()
  elseif (present(defaultVal)) then
    nodeAsString = defaultVal
  else
    call IO_error(143,ext_msg=k)
  end if

end function tDict_get_asString


!--------------------------------------------------------------------------------------------------
!> @brief Get list by key and convert to string array (1D).
!--------------------------------------------------------------------------------------------------
function tDict_get_as1dString(self,k,defaultVal) result(nodeAs1dString)

  class(tDict),     intent(in) :: self
  character(len=*), intent(in) :: k
  character(len=*), intent(in), dimension(:), optional :: defaultVal
  character(len=:), allocatable, dimension(:)          :: nodeAs1dString

  type(tList), pointer :: list


  if (self%contains(k)) then
    list => self%get_list(k)
    nodeAs1dString = list%as1dString()
  elseif (present(defaultVal)) then
    nodeAs1dString = defaultVal
  else
    call IO_error(143,ext_msg=k)
  end if

end function tDict_get_as1dString


#ifdef __GFORTRAN__
!--------------------------------------------------------------------------------------------------
!> @brief Returns string output array (1D) (hack for GNU).
!--------------------------------------------------------------------------------------------------
function output_as1dString(self)  result(output)

  class(tDict), pointer,intent(in)   ::  self
  character(len=pStringLen), allocatable, dimension(:) :: output

  type(tList), pointer :: output_list
  integer :: o

  output_list => self%get_list('output',defaultVal=emptyList)
  allocate(output(output_list%length))
  do o = 1, output_list%length
    output(o) = output_list%get_asString(o)
  end do

end function output_as1dString
#endif


!--------------------------------------------------------------------------------------------------
!> @brief Free associated memory.
!--------------------------------------------------------------------------------------------------
recursive subroutine tItem_finalize(self)

  type(tItem),intent(inout) :: self

  deallocate(self%node)
  if (associated(self%next)) deallocate(self%next)

end subroutine tItem_finalize

end module YAML_types
