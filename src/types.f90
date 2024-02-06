!--------------------------------------------------------------------------------------------------
!> @author Sharan Roongta, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Data types to create a scalar, a list, and a dictionary/hash
!> @details module describes the various functions to store and get the yaml data.
!! A node is the base class for scalar, list and dictionary, list items and dictionary entries point
!! to a node.
!--------------------------------------------------------------------------------------------------

module types
  use IO
  use prec
  use misc

  implicit none(type,external)
  private

  type, abstract, public :: tNode
    integer :: &
      length = 0
    contains
    procedure(asFormattedStr), deferred :: &
      asFormattedStr
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
      asFormattedStr => tScalar_asFormattedStr, &
      asReal => tScalar_asReal, &
      asInt  => tScalar_asInt, &
      asBool => tScalar_asBool, &
      asStr  => tScalar_asStr
  end type tScalar

  type, extends(tNode), public :: tList
    class(tItem), pointer :: &
      first => NULL(), &
      last => NULL()
    contains
    procedure :: &
      asFormattedStr => tList_asFormattedStr, &
      append            => tList_append, &
      as1dReal          => tList_as1dReal, &
      as2dReal          => tList_as2dReal, &
      as1dInt           => tList_as1dInt, &
      as1dBool          => tList_as1dBool, &
      as1dStr           => tList_as1dStr, &
      contains          => tList_contains, &
      tList_get, &
      tList_get_scalar, &
      tList_get_list, &
      tList_get_dict, &
      tList_get_asReal, &
      tList_get_as1dReal, &
      tList_get_asInt, &
      tList_get_as1dInt, &
      tList_get_asBool, &
      tList_get_as1dBool, &
      tList_get_asStr, &
      tList_get_as1dStr
    generic :: get          => tList_get
    generic :: get_scalar   => tList_get_scalar
    generic :: get_list     => tList_get_list
    generic :: get_dict     => tList_get_dict
    generic :: get_asReal   => tList_get_asReal
    generic :: get_as1dReal => tList_get_as1dReal
    generic :: get_asInt    => tList_get_asInt
    generic :: get_as1dInt  => tList_get_as1dInt
    generic :: get_asBool   => tList_get_asBool
    generic :: get_as1dBool => tList_get_as1dBool
    generic :: get_asStr    => tList_get_asStr
    generic :: get_as1dStr  => tList_get_as1dStr
    final :: tList_finalize
  end type tList

  type, extends(tList), public :: tDict
    contains
    procedure ::  &
      asFormattedStr => tDict_asFormattedStr, &
      set            => tDict_set, &
      index          => tDict_index, &
      key            => tDict_key, &
      keys           => tDict_keys, &
      contains       => tDict_contains, &
      tDict_get, &
      tDict_get_scalar, &
      tDict_get_list, &
      tDict_get_dict, &
      tDict_get_asReal, &
      tDict_get_as1dReal_sized, &
      tDict_get_as1dReal_chunked, &
      tDict_get_as2dReal, &
      tDict_get_asInt, &
      tDict_get_as1dInt, &
      tDict_get_asBool, &
      tDict_get_as1dBool, &
      tDict_get_asStr, &
      tDict_get_as1dStr
    generic :: get          => tDict_get
    generic :: get_scalar   => tDict_get_scalar
    generic :: get_list     => tDict_get_list
    generic :: get_dict     => tDict_get_dict
    generic :: get_asReal   => tDict_get_asReal
    generic :: get_as1dReal => tDict_get_as1dReal_sized
    generic :: get_as1dReal => tDict_get_as1dReal_chunked
    generic :: get_as2dReal => tDict_get_as2dReal
    generic :: get_asInt    => tDict_get_asInt
    generic :: get_as1dInt  => tDict_get_as1dInt
    generic :: get_asBool   => tDict_get_asBool
    generic :: get_as1dBool => tDict_get_as1dBool
    generic :: get_asStr    => tDict_get_asStr
    generic :: get_as1dStr  => tDict_get_as1dStr
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

    recursive function asFormattedStr(self)
      import tNode
      character(len=:), allocatable      :: asFormattedStr
      class(tNode), intent(in), target   :: self
    end function asFormattedStr

  end interface

  interface tScalar
    module procedure tScalar_init__
  end interface tScalar

  interface assignment (=)
    module procedure tScalar_assign__
  end interface assignment (=)

  public :: &
    types_init, &
    types_selfTest, &
#ifdef __GFORTRAN__
    output_as1dStr, &                                       !ToDo: Hack for GNU. Remove later
#endif
    assignment(=)

contains

!--------------------------------------------------------------------------------------------------
!> @brief Do sanity checks.
!--------------------------------------------------------------------------------------------------
subroutine types_init

  print'(/,1x,a)', '<<<+-  types init  -+>>>'

  call types_selfTest()

end subroutine types_init


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some type bound procedures.
!--------------------------------------------------------------------------------------------------
subroutine types_selfTest()

  scalar: block
    type(tScalar), target :: s
    type(tScalar), pointer :: s_pointer


    s_pointer => s%asScalar()
    s = '1'
    if (s%asInt() /= 1)                error stop 'tScalar_asInt'
    if (s_pointer%asInt() /= 1)        error stop 'tScalar_asInt(pointer)'
    if (dNeq(s%asReal(),1.0_pREAL))    error stop 'tScalar_asReal'
    s = 'true'
    if (.not. s%asBool())              error stop 'tScalar_asBool'
    if (.not. s_pointer%asBool())      error stop 'tScalar_asBool(pointer)'
    if (s%asStr() /= 'true')           error stop 'tScalar_asStr'
    if (s%asFormattedStr() /= 'true')  error stop 'tScalar_asFormattedStr'


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
    if (l%contains('1'))                               error stop 'empty tList_contains'
    if (l_pointer%contains('1'))                       error stop 'empty tList_contains(pointer)'
    call l%append(s1)
    call l%append(s2)
    if (l%length /= 2)                                 error stop 'tList%len'
    if (dNeq(l%get_asReal(1),1.0_pREAL))               error stop 'tList_get_asReal'
    if (l%get_asInt(1) /= 1)                           error stop 'tList_get_asInt'
    if (l%get_asStr(2) /= '2')                         error stop 'tList_get_asStr'
    if (any(l%as1dInt() /= [1,2]))                     error stop 'tList_as1dInt'
    if (any(dNeq(l%as1dReal(),real([1.0,2.0],pREAL)))) error stop 'tList_as1dReal'
    s1 = 'true'
    s2 = 'false'
    if (any(l%as1dBool() .neqv. [.true.,.false.]))     error stop 'tList_as1dBool'
    if (any(l%as1dStr() /= ['true ','false']))         error stop 'tList_as1dStr'
    if (l%asFormattedStr() /= '[true, false]')         error stop 'tList_asFormattedStr'
    if (   .not. l%contains('true') &
      .or. .not. l%contains('false'))                  error stop 'tList_contains'

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
    if (d%contains('one-two'))                        error stop 'empty tDict_contains'
    if (d_pointer%contains('one-two'))                error stop 'empty tDict_contains(pointer)'
    if (d%get_asInt('one-two',defaultVal=-1) /= -1)   error stop 'empty tDict_get'
    call d%set('one-two',l)
    call d%set('three',s3)
    call d%set('four',s4)
    if (d%asFormattedStr() /= '{one-two: [1, 2], three: 3, four: 4}') &
                                                      error stop 'tDict_asFormattedStr'
    if (d%get_asInt('three') /= 3)                    error stop 'tDict_get_asInt'
    if (dNeq(d%get_asReal('three'),3.0_pREAL))        error stop 'tDict_get_asReal'
    if (any(d%get_as1dReal('one-two') /= real([1,2],pReal))) &
                                                      error stop 'tDict_get_as1dReal'
    if (any(d%get_as1dReal('three',requiredSize=3) /= real([3,3,3],pReal))) &
                                                      error stop 'tDict_get_as1dReal/size'
    if (d%get_asStr('three') /= '3')                  error stop 'tDict_get_asStr'
    if (any(d%get_as1dInt('one-two') /= [1,2]))       error stop 'tDict_get_as1dInt'
    if (any(d%get_as1dInt('three',requiredSize=3) /= [3,3,3])) &
                                                      error stop 'tDict_get_as1dInt/size'
    call d%set('one-two',s4)
    if (d%asFormattedStr() /= '{one-two: 4, three: 3, four: 4}') &
                                                      error stop 'tDict_set overwrite'
    if (   .not. d%contains('one-two') &
      .or. .not. d%contains('three') &
      .or. .not. d%contains('four') &
      )                                               error stop 'tDict_contains'


  end block dict

#ifdef __GFORTRAN__
  dict_get_as1dReal_chunked: block
    type(tDict),   pointer :: d
    type(tList),   pointer :: l_outer, l_inner
    type(tScalar), pointer :: s1,s2,s3
    real(pREAL), dimension(5) :: a


    allocate(s1)
    allocate(s2)
    allocate(s3)
    s1 = '1.'
    s2 = '2'
    s3 = '3.0'

    allocate(l_inner)
    call l_inner%append(s1)
    call l_inner%append(s2)

    allocate(l_outer)
    call l_outer%append(l_inner)
    call l_outer%append(s3)

    allocate(d)
    call d%set('list',l_outer)
    call d%set('scalar',s1)

    a = d%get_as1dReal('list',requiredChunks=[2,3])
    if (any(dNeq(a,real([1.0,2.0,3.0,3.0,3.0],pReal)))) &
      error stop 'dict_get_as1dReal_shape list'

    if (any(dNeq(d%get_as1dReal('non-existing',a,[2,3]),a))) &
      error stop 'dict_get_as1dReal_shape default individual'

    a = real([42.0, 42.0, 5.0, 5.0, 5.0],pREAL)
    if (any(dNeq(d%get_as1dReal('non-existing',[42._pREAL, 5._pREAL],[2,3]),a))) &
      error stop 'dict_get_as1dReal_shape default group'

    if (any(dNeq(d%get_as1dReal('scalar',requiredChunks=[3,5,2]),misc_ones(10)))) &
      error stop 'dict_get_as1dReal_shape scalar'

  end block dict_get_as1dReal_chunked
#endif

end subroutine types_selfTest


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
recursive function tScalar_asFormattedStr(self) result(str)

  class (tScalar), intent(in), target    :: self
  character(len=:), allocatable          :: str


  str = trim(self%value)

end function tScalar_asFormattedStr


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
      call IO_error(706,'"'//trim(self%asFormattedStr())//'" is not a scalar')
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
      call IO_error(706,'"'//trim(self%asFormattedStr())//'" is not a list')
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
      call IO_error(706,'"'//trim(self%asFormattedStr())//'" is not a dict')
  end select

end function tNode_asDict


!--------------------------------------------------------------------------------------------------
!> @brief Convert to real.
!--------------------------------------------------------------------------------------------------
function tScalar_asReal(self)

  class(tScalar), intent(in), target :: self
  real(pREAL) :: tScalar_asReal


  tScalar_asReal = IO_strAsReal(self%value)

end function tScalar_asReal


!--------------------------------------------------------------------------------------------------
!> @brief Convert to int.
!--------------------------------------------------------------------------------------------------
function tScalar_asInt(self)

  class(tScalar), intent(in), target :: self
  integer :: tScalar_asInt


  tScalar_asInt = IO_strAsInt(self%value)

end function tScalar_asInt


!--------------------------------------------------------------------------------------------------
!> @brief Convert to bool.
!--------------------------------------------------------------------------------------------------
function tScalar_asBool(self)

  class(tScalar), intent(in), target :: self
  logical :: tScalar_asBool


  tScalar_asBool = IO_strAsBool(self%value)

end function tScalar_asBool


!--------------------------------------------------------------------------------------------------
!> @brief Convert to string.
!--------------------------------------------------------------------------------------------------
function tScalar_asStr(self)

  class(tScalar), intent(in), target :: self
  character(len=:), allocatable :: tScalar_asStr


  tScalar_asStr = self%value

end function tScalar_asStr


!--------------------------------------------------------------------------------------------------
!> @brief Format as string (YAML flow style).
!--------------------------------------------------------------------------------------------------
recursive function tList_asFormattedStr(self) result(str)

  class(tList),intent(in),target      :: self

  type(tItem), pointer  :: item
  character(len=:), allocatable :: str
  integer :: i

  str = '['
  item => self%first
  do i = 2, self%length
    str = str//item%node%asFormattedStr()//', '
    item => item%next
  end do
  str = str//item%node%asFormattedStr()//']'

end function tList_asFormattedStr


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
!> @brief Convert to real array (1D).
!--------------------------------------------------------------------------------------------------
function tList_as1dReal(self)

  class(tList), intent(in), target :: self
  real(pREAL), dimension(:), allocatable :: tList_as1dReal

  integer :: i
  type(tItem),   pointer :: item
  type(tScalar), pointer :: scalar


  allocate(tList_as1dReal(self%length))
  item => self%first
  do i = 1, self%length
    scalar => item%node%asScalar()
    tList_as1dReal(i) = scalar%asReal()
    item => item%next
  end do

end function tList_as1dReal


!--------------------------------------------------------------------------------------------------
!> @brief Convert to real array (2D).
!--------------------------------------------------------------------------------------------------
function tList_as2dReal(self)

  class(tList), intent(in), target :: self
  real(pREAL), dimension(:,:), allocatable :: tList_as2dReal

  integer :: i
  type(tList), pointer :: row_data


  row_data => self%get_list(1)
  allocate(tList_as2dReal(self%length,row_data%length))

  do i = 1, self%length
    row_data => self%get_list(i)
    if (row_data%length /= size(tList_as2dReal,2)) call IO_error(709,ext_msg='inconsistent column count in tList_as2dReal')
    tList_as2dReal(i,:) = self%get_as1dReal(i)
  end do

end function tList_as2dReal


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
function tList_as1dStr(self)

  class(tList), intent(in), target :: self
#ifdef __GFORTRAN__
  character(len=pSTRLEN), allocatable, dimension(:) :: tList_as1dStr
#else
  character(len=:), allocatable, dimension(:) :: tList_as1dStr
#endif

  integer :: j
  type(tItem),   pointer :: item
  type(tScalar), pointer :: scalar


#ifdef __GFORTRAN__
  allocate(tList_as1dStr(self%length))
#else
  integer :: len_max
  len_max = 0
  item => self%first
  do j = 1, self%length
    scalar => item%node%asScalar()
    len_max = max(len_max, len_trim(scalar%asStr()))
    item => item%next
  end do

  allocate(character(len=len_max) :: tList_as1dStr(self%length))
#endif
  item => self%first
  do j = 1, self%length
    scalar => item%node%asScalar()
    tList_as1dStr(j) = scalar%asStr()
    item => item%next
  end do

end function tList_as1dStr


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


  if (i < 1 .or. i > self%length) call IO_error(150,ext_msg='tList_get @ '//IO_intAsStr(i) &
                                                                  //' of '//IO_intAsStr(self%length) )
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
!> @brief Get scalar by index and convert to real.
!--------------------------------------------------------------------------------------------------
function tList_get_asReal(self,i) result(nodeAsReal)

  class(tList), intent(in) :: self
  integer,      intent(in) :: i
  real(pREAL) :: nodeAsReal

  class(tScalar),  pointer :: scalar


  scalar => self%get_scalar(i)
  nodeAsReal = scalar%asReal()

end function tList_get_asReal


!--------------------------------------------------------------------------------------------------
!> @brief Get list by index and convert to real array (1D).
!--------------------------------------------------------------------------------------------------
function tList_get_as1dReal(self,i) result(nodeAs1dReal)

  class(tList), intent(in) :: self
  integer,      intent(in) :: i
  real(pREAL), dimension(:), allocatable :: nodeAs1dReal

  class(tList),  pointer :: list


  list => self%get_list(i)
  nodeAs1dReal = list%as1dReal()

end function tList_get_as1dReal


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
function tList_get_asStr(self,i) result(nodeAsStr)

  class(tList), intent(in) :: self
  integer,      intent(in) :: i
  character(len=:), allocatable :: nodeAsStr

  class(tScalar),  pointer :: scalar


  scalar => self%get_scalar(i)
  nodeAsStr = scalar%asStr()

end function tList_get_asStr


!--------------------------------------------------------------------------------------------------
!> @brief Get list by index and convert to string array (1D).
!--------------------------------------------------------------------------------------------------
function tList_get_as1dStr(self,i) result(nodeAs1dStr)

  class(tList), intent(in) :: self
  integer,      intent(in) :: i
  character(len=:), allocatable, dimension(:) :: nodeAs1dStr

  type(tList),  pointer :: list


  list => self%get_list(i)
  nodeAs1dStr = list%as1dStr()

end function tList_get_as1dStr


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
recursive function tDict_asFormattedStr(self) result(str)

  class(tDict),intent(in),target      :: self

  type(tItem),pointer  :: item
  character(len=:), allocatable :: str
  integer :: i


  str = '{'
  item => self%first
  do i = 2, self%length
    str = str//trim(item%key)//': '//item%node%asFormattedStr()//', '
    item => item%next
  end do
  str = str//trim(item%key)//': '//item%node%asFormattedStr()//'}'

end function tDict_asFormattedStr


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


  if (i < 1 .or. i > self%length) call IO_error(150,ext_msg='tDict_key @ '//IO_intAsStr(i) &
                                                                  //' of '//IO_intAsStr(self%length) )
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

  character(len=pSTRLEN), dimension(:), allocatable :: temp
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
!> @brief Get scalar by key and convert to real.
!--------------------------------------------------------------------------------------------------
function tDict_get_asReal(self,k,defaultVal) result(nodeAsReal)

  class(tDict),     intent(in) :: self
  character(len=*), intent(in) :: k
  real(pREAL),      intent(in), optional :: defaultVal
  real(pREAL) :: nodeAsReal

  type(tScalar), pointer :: scalar


  if (self%contains(k)) then
    scalar => self%get_scalar(k)
    nodeAsReal = scalar%asReal()
  elseif (present(defaultVal)) then
    nodeAsReal = defaultVal
  else
    call IO_error(143,ext_msg=k)
  end if

end function tDict_get_asReal


!--------------------------------------------------------------------------------------------------
!> @brief Get list by key and convert to real array (1D).
!> @details If a size is required, scalars are valid input and are broadcasted to the required size.
!--------------------------------------------------------------------------------------------------
function tDict_get_as1dReal_sized(self,k,defaultVal,requiredSize) result(nodeAs1dReal)

  class(tDict),     intent(in) :: self
  character(len=*), intent(in) :: k
  real(pREAL),      intent(in), dimension(:), optional :: defaultVal
  integer,          intent(in),               optional :: requiredSize
  real(pREAL), dimension(:), allocatable :: nodeAs1dReal

  class(tNode), pointer :: content


  if (self%contains(k)) then
    content => self%get(k)
    select type(content)
      class is(tScalar)
        if (present(requiredSize)) then
          allocate(nodeAs1dReal(requiredSize),source = content%asReal())
        else
          call IO_error(706,'"'//trim(content%asFormattedStr())//'" is not a list of reals')
        end if
      class is(tList)
        nodeAs1dReal = content%as1dReal()
    end select
  elseif (present(defaultVal)) then
    nodeAs1dReal = defaultVal
  else
    call IO_error(143,ext_msg=k)
  end if

  if (present(requiredSize)) then
    if (requiredSize /= size(nodeAs1dReal)) &
      call IO_error(146,ext_msg=k, &
                    label1='actual',ID1=size(nodeAs1dReal), &
                    label2='required',ID2=requiredSize)
  end if

end function tDict_get_as1dReal_sized


!--------------------------------------------------------------------------------------------------
!> @brief Get entry by key and convert to real array (1D).
!> @details Values will be broadcasted. A List content can be composed from mixture of scalar
!> or list entries. [2., [1., 3.]] with required chunks [3, 2] gives [2., 2., 2., 1., 3.].
!--------------------------------------------------------------------------------------------------
function tDict_get_as1dReal_chunked(self,k,defaultVal,requiredChunks) result(nodeAs1dReal)

  class(tDict),     intent(in) :: self
  character(len=*), intent(in) :: k
  real(pREAL),      intent(in), dimension(:), optional :: defaultVal
  integer,          intent(in), dimension(:)           :: requiredChunks
  real(pREAL),                  dimension(sum(requiredChunks)) :: nodeAs1dReal

  type(tList), pointer :: list_outer, list_inner
  class(tNode), pointer :: node_outer, node_inner
  integer :: i


  if (self%contains(k)) then
    node_outer => self%get(k)
    select type(node_outer)
      class is(tScalar)
        nodeAs1dReal = node_outer%asReal()
      class is(tList)
        list_outer => self%get_list(k)
        if (list_outer%length /= size(requiredChunks)) &
          call IO_error(709,'list "'//list_outer%asFormattedStr()//'" is not of length '//IO_intAsStr(size(requiredChunks)))
        do i = 1, size(requiredChunks)
          node_inner => list_outer%get(i)
          select type(node_inner)
            class is(tScalar)
              nodeAs1dReal(sum(requiredChunks(:i-1))+1:sum(requiredChunks(:i))) = node_inner%asReal()
            class is(tList)
              list_inner => node_inner%asList()
              if (size(list_inner%as1dReal()) /= requiredChunks(i)) &
                  call IO_error(709,'entry "'//k//'" is not of length '//IO_intAsStr(requiredChunks(i)),&
                                'position',i)
              nodeAs1dReal(sum(requiredChunks(:i-1))+1:sum(requiredChunks(:i))) = list_inner%as1dReal()
            class default
              call IO_error(706,'entry "'//k//'" is neither scalar nor list','position',i)
          end select
        end do
    end select
  elseif (present(defaultVal)) then
    if (size(defaultVal) == size(nodeAs1dReal)) then
      nodeAs1dReal = defaultVal
    elseif (size(defaultVal) == size(requiredChunks)) then
      do i = 1, size(requiredChunks)
        nodeAs1dReal(sum(requiredChunks(:i-1))+1:sum(requiredChunks(:i))) = defaultVal(i)
      end do
    else
      call IO_error(709,'default values not of required shape')
    end if
  else
    call IO_error(143,ext_msg=k)
  end if

end function tDict_get_as1dReal_chunked


!--------------------------------------------------------------------------------------------------
!> @brief Get list of lists by key and convert to real array (2D).
!--------------------------------------------------------------------------------------------------
function tDict_get_as2dReal(self,k,defaultVal,requiredShape) result(nodeAs2dReal)

  class(tDict),     intent(in) :: self
  character(len=*), intent(in) :: k
  real(pREAL),      intent(in), dimension(:,:), optional :: defaultVal
  integer,          intent(in), dimension(2),   optional :: requiredShape
  real(pREAL), dimension(:,:), allocatable :: nodeAs2dReal

  type(tList), pointer :: list


  if (self%contains(k)) then
    list => self%get_list(k)
    nodeAs2dReal = list%as2dReal()
  elseif (present(defaultVal)) then
    nodeAs2dReal = defaultVal
  else
    call IO_error(143,ext_msg=k)
  end if

  if (present(requiredShape)) then
    if (any(requiredShape /= shape(nodeAs2dReal))) call IO_error(146,ext_msg=k)
  end if

end function tDict_get_as2dReal


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
!> @details If a size is required, scalars are valid input and are broadcasted to the required size.
!--------------------------------------------------------------------------------------------------
function tDict_get_as1dInt(self,k,defaultVal,requiredSize) result(nodeAs1dInt)

  class(tDict),          intent(in) :: self
  character(len=*),      intent(in) :: k
  integer, dimension(:), intent(in), optional :: defaultVal
  integer,               intent(in), optional :: requiredSize
  integer, dimension(:), allocatable :: nodeAs1dInt

  class(tNode), pointer :: content


  if (self%contains(k)) then
    content => self%get(k)
    select type(content)
      class is(tScalar)
        if (present(requiredSize)) then
          allocate(nodeAs1dInt(requiredSize),source = content%asInt())
        else
          call IO_error(706,'"'//trim(content%asFormattedStr())//'" is not a list of integers')
        end if
      class is(tList)
        nodeAs1dInt = content%as1dInt()
    end select
  elseif (present(defaultVal)) then
    nodeAs1dInt = defaultVal
  else
    call IO_error(143,ext_msg=k)
  end if

  if (present(requiredSize)) then
    if (requiredSize /= size(nodeAs1dInt)) &
      call IO_error(146,ext_msg=k, &
                    label1='actual',ID1=size(nodeAs1dInt), &
                    label2='required',ID2=requiredSize)
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
function tDict_get_asStr(self,k,defaultVal) result(nodeAsStr)

  class(tDict),     intent(in) :: self
  character(len=*), intent(in) :: k
  character(len=*), intent(in), optional :: defaultVal
  character(len=:), allocatable :: nodeAsStr

  type(tScalar), pointer :: scalar


  if (self%contains(k)) then
    scalar => self%get_scalar(k)
    nodeAsStr = scalar%asStr()
  elseif (present(defaultVal)) then
    nodeAsStr = defaultVal
  else
    call IO_error(143,ext_msg=k)
  end if

end function tDict_get_asStr


!--------------------------------------------------------------------------------------------------
!> @brief Get list by key and convert to string array (1D).
!--------------------------------------------------------------------------------------------------
function tDict_get_as1dStr(self,k,defaultVal) result(nodeAs1dStr)

  class(tDict),     intent(in) :: self
  character(len=*), intent(in) :: k
  character(len=*), intent(in), dimension(:), optional :: defaultVal
  character(len=:), allocatable, dimension(:)          :: nodeAs1dStr

  type(tList), pointer :: list


  if (self%contains(k)) then
    list => self%get_list(k)
    nodeAs1dStr = list%as1dStr()
  elseif (present(defaultVal)) then
    nodeAs1dStr = defaultVal
  else
    call IO_error(143,ext_msg=k)
  end if

end function tDict_get_as1dStr


#ifdef __GFORTRAN__
!--------------------------------------------------------------------------------------------------
!> @brief Returns string output array (1D) (hack for GNU).
!--------------------------------------------------------------------------------------------------
function output_as1dStr(self)  result(output)

  class(tDict), pointer,intent(in)   ::  self
  character(len=pSTRLEN), allocatable, dimension(:) :: output

  type(tList), pointer :: output_list
  integer :: o

  output_list => self%get_list('output',defaultVal=emptyList)
  allocate(output(output_list%length))
  do o = 1, output_list%length
    output(o) = output_list%get_asStr(o)
  end do

end function output_as1dStr
#endif


!--------------------------------------------------------------------------------------------------
!> @brief Free associated memory.
!--------------------------------------------------------------------------------------------------
recursive subroutine tItem_finalize(self)

  type(tItem),intent(inout) :: self

  deallocate(self%node)
  if (associated(self%next)) deallocate(self%next)

end subroutine tItem_finalize

end module types
