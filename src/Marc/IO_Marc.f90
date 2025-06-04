!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Basic IO functions for parsing MSC.Marc input file(s)
!--------------------------------------------------------------------------------------------------
module IO_Marc

  use prec
  use IO

  public :: &
    strPos, &
    strValue, &
    realValue, &
    continuousIntValues, &
    containsRange

contains

!--------------------------------------------------------------------------------------------------
!> @brief Locate all whitespace-separated chunks in the given string and
!! return array containing their number and each left/right position
!! for later use in IO_xxxVal.
!! Array size is dynamically adjusted to number of chunks found in string
!! IMPORTANT: first element contains number of chunks!
!--------------------------------------------------------------------------------------------------
pure function strPos(str)

  character(len=*),                  intent(in) :: str                                              !< string in which chunk positions are searched for
  integer, dimension(:), allocatable            :: strPos

  integer :: left, right


  allocate(strPos(1), source=0)
  right = 0

  do while (verify(str(right+1:),IO_WHITESPACE)>0)
    left  = right + verify(str(right+1:),IO_WHITESPACE)
    right = left + scan(str(left:),IO_WHITESPACE) - 2
    strPos = [strPos,left,right]
    strPos(1) = strPos(1)+1
    endOfStr: if (right < left) then
      strPos(strPos(1)*2+1) = len_trim(str)
      exit
    end if endOfStr
  end do

end function strPos

!--------------------------------------------------------------------------------------------------
!> @brief Read string value at myChunk from string.
!--------------------------------------------------------------------------------------------------
function strValue(str,chunkPos,myChunk)

  character(len=*),             intent(in) :: str                                                   !< raw input with known start and end of each chunk
  integer,   dimension(:),      intent(in) :: chunkPos                                              !< positions of start and end of each tag/chunk in given string
  integer,                      intent(in) :: myChunk                                               !< position number of desired chunk
  character(len=:), allocatable            :: strValue


  validChunk: if (myChunk > chunkPos(1) .or. myChunk < 1) then
    strValue = ''
    call IO_error(110,'strValue: "'//trim(str)//'"',label1='chunk',ID1=myChunk)
  else validChunk
    strValue = str(chunkPos(myChunk*2):chunkPos(myChunk*2+1))
  end if validChunk

end function strValue

!--------------------------------------------------------------------------------------------------
!> @brief Read integer value at myChunk from string.
!--------------------------------------------------------------------------------------------------
integer function intValue(str,chunkPos,myChunk)

  character(len=*),      intent(in) :: str                                                          !< raw input with known start and end of each chunk
  integer, dimension(:), intent(in) :: chunkPos                                                     !< positions of start and end of each tag/chunk in given string
  integer,               intent(in) :: myChunk                                                      !< position number of desired chunk


  intValue = IO_strAsInt(strValue(str,chunkPos,myChunk))

end function intValue


!--------------------------------------------------------------------------------------------------
!> @brief Read real value at myChunk from string.
!--------------------------------------------------------------------------------------------------
real(pREAL) function realValue(str,chunkPos,myChunk)

  character(len=*),        intent(in) :: str                                                        !< raw input with known start and end of each chunk
  integer,   dimension(:), intent(in) :: chunkPos                                                   !< positions of start and end of each tag/chunk in given string
  integer,                 intent(in) :: myChunk                                                    !< position number of desired chunk


  realValue = IO_strAsReal(strValue(str,chunkPos,myChunk))

end function realValue


!--------------------------------------------------------------------------------------------------
!> @brief return integer list corresponding to items in consecutive lines.
!! First integer in array is counter
!> @details ints concatenated by "c" as last char, range of a "to" b, or named set
!--------------------------------------------------------------------------------------------------
function continuousIntValues(fileContent,maxN,lookupName,lookupMap,lookupMaxN)

  character(len=*), dimension(:),   intent(in) :: fileContent                                       !< file content, separated per lines
  integer,                          intent(in) :: maxN
  integer,                          intent(in) :: lookupMaxN
  integer,          dimension(:,:), intent(in) :: lookupMap
  character(len=*), dimension(:),   intent(in) :: lookupName

  integer,           dimension(1+maxN)         :: continuousIntValues

  integer :: l,i,first,last
  integer, allocatable, dimension(:) :: chunkPos
  logical :: rangeGeneration

  continuousIntValues = 0
  rangeGeneration = .false.

  do l = 1, size(fileContent)
    chunkPos = strPos(fileContent(l))
    if (chunkPos(1) < 1) then                                                                       ! empty line
      exit
    elseif (verify(strValue(fileContent(l),chunkPos,1),'0123456789') > 0) then                   ! a non-int, i.e. set name
      do i = 1, lookupMaxN                                                                          ! loop over known set names
        if (strValue(fileContent(l),chunkPos,1) == lookupName(i)) then                           ! found matching name
          continuousIntValues = lookupMap(:,i)                                                      ! return resp. entity list
          exit
        end if
      end do
      exit
    elseif (containsRange(fileContent(l),chunkPos)) then
      first = intValue(fileContent(l),chunkPos,1)
      last  = intValue(fileContent(l),chunkPos,3)
      do i = first, last, sign(1,last-first)
        continuousIntValues(1) = continuousIntValues(1) + 1
        continuousIntValues(1+continuousIntValues(1)) = i
      end do
      exit
    else
      do i = 1,chunkPos(1)-1                                                                        ! interpret up to second to last value
        continuousIntValues(1) = continuousIntValues(1) + 1
        continuousIntValues(1+continuousIntValues(1)) = intValue(fileContent(l),chunkPos,i)
      end do
      if ( IO_lc(strValue(fileContent(l),chunkPos,chunkPos(1))) /= 'c' ) then                    ! line finished, read last value
        continuousIntValues(1) = continuousIntValues(1) + 1
        continuousIntValues(1+continuousIntValues(1)) = intValue(fileContent(l),chunkPos,chunkPos(1))
        exit
      end if
    end if
  end do

end function continuousIntValues


!--------------------------------------------------------------------------------------------------
!> @brief return whether a line contains a range ('X to Y')
!--------------------------------------------------------------------------------------------------
logical function containsRange(str,chunkPos)

  character(len=*),      intent(in) :: str
  integer, dimension(:), intent(in) :: chunkPos                                                     !< positions of start and end of each tag/chunk in given string


  containsRange = .False.
  if (chunkPos(1) == 3) then
    if (IO_lc(strValue(str,chunkPos,2)) == 'to') containsRange = .True.
  end if

end function containsRange

end module IO_Marc
