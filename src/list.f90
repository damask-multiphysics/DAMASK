module chained_list
  use prec
  implicit none
  
  type tPartitionedString
    character(len=:),      allocatable :: val
    integer(pInt), dimension(:), allocatable :: pos
  end type
 
  type, public :: tPartitionedStringList
    type(tPartitionedString)    :: string
    type(tPartitionedStringList),  pointer :: next => null()
    contains 
      procedure :: add            => add
      procedure :: getFloat       => getFloat
      procedure :: getFloatArray  => getFloatArray
      procedure :: getStrings     => getStrings
      procedure :: keyExists     => keyExists
  end type tPartitionedStringList
 

  contains
    subroutine add(self,string,stringPos)
      implicit none
      class(tPartitionedStringList) :: self
      type(tPartitionedStringList), pointer :: new,tmp
      character(len=*), intent(in) :: string
      integer(pInt), dimension(:), intent(in) :: stringPos

      allocate(new)

      new%string%val=string
      new%string%pos=stringPos

      if (.not. associated(self%next)) then
        self%next => new
      else
        tmp => self%next
        self%next => new
        self%next%next => tmp
      end if
  
    end subroutine add


! gets float value, if key is not found exits with error unless default is given
    function getFloat(self,key,default)
      use IO

      implicit none
      real(pReal) :: getFloat

      class(tPartitionedStringList), intent(in) :: self
      character(len=*), intent(in) :: key
      real(pReal), intent(in), optional :: default
      type(tPartitionedStringList), pointer :: tmp

      tmp => self%next
      do 
        if (.not. associated(tmp)) then
          if(present(default)) then
            getFloat = default
            exit
          else
            call IO_error(1_pInt,ext_msg=key)
          endif
        endif
        if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
          if (tmp%string%pos(1) > 2) call IO_error(1_pInt,ext_msg=key)
          getFloat = IO_FloatValue(tmp%string%val,tmp%string%pos,2)
          exit
        endif
        tmp => tmp%next
      end do
    end function

! reports wether a key exists at least once
    function keyExists(self,key)
      use IO

      implicit none
      logical :: keyExists

      class(tPartitionedStringList), intent(in) :: self
      character(len=*), intent(in) :: key
      type(tPartitionedStringList), pointer :: tmp

      keyExists = .false.

      tmp => self%next
      do 
        if (.not. associated(tmp)) exit
        if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
          keyExists = .true.
          exit
        endif
        tmp => tmp%next
      end do
    end function

    function getFloatArray(self,key)
      use IO

      implicit none
      real(pReal),dimension(:),allocatable :: getFloatArray

      class(tPartitionedStringList), intent(in) :: self
      character(len=*), intent(in) :: key
      type(tPartitionedStringList), pointer :: tmp
      integer(pInt) :: i

      allocate(getFloatArray(0))

      tmp => self%next
      do 
        if (.not. associated(tmp)) exit
        if (trim(IO_stringValue(tmp%string%val,tmp%string%pos,1))==trim(key)) then
          do i = 2_pInt, tmp%string%pos(1)
            getFloatArray = [getFloatArray,IO_FloatValue(tmp%string%val,tmp%string%pos,i)]
          enddo
          exit
        endif
        tmp => tmp%next
      end do
    end function


    function getStrings(self,key)
      use IO

      implicit none
      character(len=64),dimension(:),allocatable :: getStrings
      character(len=64) :: str

      class(tPartitionedStringList), intent(in) :: self
      character(len=*), intent(in) :: key
      type(tPartitionedStringList), pointer :: tmp
      integer(pInt) :: i

      tmp => self%next
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
