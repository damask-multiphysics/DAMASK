!----------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Sharan Roongta, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Parser for YAML files
!> @details module converts a YAML input file to an equivalent YAML flow style which is then parsed.
!----------------------------------------------------------------------------------------------------
module YAML_parse
  use prec
  use IO
  use YAML_types

  implicit none
  private

  public :: &
    YAML_parse_init, &
    YAML_parse_file

contains

!--------------------------------------------------------------------------------------------------
!> @brief Do sanity checks.
!--------------------------------------------------------------------------------------------------
subroutine YAML_parse_init

  call selfTest

end subroutine YAML_parse_init


!--------------------------------------------------------------------------------------------------
!> @brief Parse a YAML file into a a structure of nodes.
!--------------------------------------------------------------------------------------------------
function YAML_parse_file(fname) result(node)

  character(len=*), intent(in) :: fname
  class (tNode), pointer :: node

  node => parse_flow(to_flow(IO_read(fname)))

end function YAML_parse_file


!--------------------------------------------------------------------------------------------------
!> @brief reads the flow style string and stores it in the form of dictionaries, lists and scalars.
!> @details A node type pointer can either point to a dictionary, list or scalar type entities.
!--------------------------------------------------------------------------------------------------
recursive function parse_flow(YAML_flow) result(node)

  character(len=*), intent(in)    :: YAML_flow                                                      !< YAML file in flow style
  class (tNode), pointer          :: node

  class (tNode),    pointer       :: &
    myVal
  character(len=:), allocatable   :: &
    flow_string, &
    key
  integer :: &
    e, &                                                                                            ! end position of dictionary or list
    s, &                                                                                            ! start position of dictionary or list
    d                                                                                               ! position of key: value separator (':')

  flow_string = trim(adjustl(YAML_flow(:)))
  if (len_trim(flow_string) == 0) then
    node => emptyDict
    return
  elseif (flow_string(1:1) == '{') then                                                             ! start of a dictionary
    e = 1
    allocate(tDict::node)
    do while (e < len_trim(flow_string))
      s = e
      d = s + scan(flow_string(s+1:),':')
      e = d + find_end(flow_string(d+1:),'}')

      key = trim(adjustl(flow_string(s+1:d-1)))
      myVal => parse_flow(flow_string(d+1:e-1))                                                     ! parse items (recursively)

      select type (node)
        class is (tDict)
          call node%set(key,myVal)
      end select
    end do
  elseif (flow_string(1:1) == '[') then                                                             ! start of a list
    e = 1
    allocate(tList::node)
    do while (e < len_trim(flow_string))
      s = e
      e = s + find_end(flow_string(s+1:),']')
      myVal => parse_flow(flow_string(s+1:e-1))                                                     ! parse items (recursively)

      select type (node)
        class is (tList)
          call node%append(myVal)
      end select
    end do
  else                                                                                              ! scalar value
    allocate(tScalar::node)
      select type (node)
        class is (tScalar)
          node = trim(adjustl(flow_string))
      end select
  endif

end function parse_flow


!--------------------------------------------------------------------------------------------------
!> @brief finds location of chunk end: ',' or '}' or  ']'
!> @details leaves nested lists ( '[...]' and dicts '{...}') intact
!--------------------------------------------------------------------------------------------------
integer function find_end(str,e_char)

  character(len=*), intent(in) :: str                                                               !< chunk of YAML flow string
  character,        intent(in) :: e_char                                                            !< end of list/dict  ( '}' or ']')

  integer                      :: N_sq, &                                                           !< number of open square brackets
                                  N_cu, &                                                           !< number of open curly brackets
                                  i

  N_sq = 0
  N_cu = 0
  do i = 1, len_trim(str)
    if (N_sq==0 .and. N_cu==0 .and. scan(str(i:i),e_char//',') == 1) exit
    N_sq = N_sq + merge(1,0,str(i:i) == '[')
    N_cu = N_cu + merge(1,0,str(i:i) == '{')
    N_sq = N_sq - merge(1,0,str(i:i) == ']')
    N_cu = N_cu - merge(1,0,str(i:i) == '}')
  enddo
  find_end = i

end function find_end


!--------------------------------------------------------------------------------------------------
! @brief Returns Indentation.
! @details It determines the indentation level for a given block/line.
! In cases for nested lists, an offset is added to determine the indent of the item block (skip
! leading dashes)
!--------------------------------------------------------------------------------------------------
integer function indentDepth(line,offset)

  character(len=*), intent(in) :: line
  integer, optional,intent(in) :: offset

  indentDepth = verify(line,IO_WHITESPACE) -1
  if(present(offset)) indentDepth = indentDepth + offset

end function indentDepth


!--------------------------------------------------------------------------------------------------
! @brief check whether a string is in flow style, i.e. starts with '{' or '['
!--------------------------------------------------------------------------------------------------
logical function isFlow(line)

  character(len=*), intent(in) :: line

  isFlow = index(adjustl(line),'[') == 1 .or. index(adjustl(line),'{') == 1

end function isFlow


!--------------------------------------------------------------------------------------------------
! @brief check whether a string is a scalar item, i.e. starts without any special symbols
!--------------------------------------------------------------------------------------------------
logical function isScalar(line)

  character(len=*), intent(in) :: line

  isScalar = (.not.isKeyValue(line) .and. .not.isKey(line) .and. .not.isListItem(line) &
                                                            .and. .not.isFlow(line))

end function isScalar


!--------------------------------------------------------------------------------------------------
! @brief check whether a string is a list item, i.e. starts with '-'
!--------------------------------------------------------------------------------------------------
logical function isListItem(line)

  character(len=*), intent(in) :: line

  isListItem = .false.
  if(len_trim(adjustl(line))> 2 .and. index(trim(adjustl(line)), '-') == 1) then
    isListItem = scan(trim(adjustl(line)),' ') == 2
  else
    isListItem = trim(adjustl(line)) == '-'
  endif

end function isListItem


!--------------------------------------------------------------------------------------------------
! @brief check whether a string contains a key value pair of the for '<key>: <value>'
!--------------------------------------------------------------------------------------------------
logical function isKeyValue(line)

  character(len=*), intent(in) :: line
  isKeyValue = .false.

  if( .not. isKey(line) .and. index(IO_rmComment(line),':') > 0 .and. .not. isFlow(line)) then
    if(index(IO_rmComment(line),': ') > 0) then
      isKeyValue = .true.
    else
      call IO_error(704,ext_msg=line)
    endif
  endif

end function isKeyValue


!--------------------------------------------------------------------------------------------------
! @brief check whether a string contains a key without a value, i.e. it ends with ':'
! ToDo: check whether this is safe for trailing spaces followed by a new line character
!--------------------------------------------------------------------------------------------------
logical function isKey(line)

  character(len=*), intent(in) :: line

  if(len(IO_rmComment(line)) == 0) then
    isKey = .false.
  else
    isKey = index(IO_rmComment(line),':',back=.false.) == len(IO_rmComment(line)) .and. &
            index(IO_rmComment(line),':',back=.true.)  == len(IO_rmComment(line)) .and. &
            .not. isFlow(line)
  endif

end function isKey


!--------------------------------------------------------------------------------------------------
! @brief check whether a string is a list in flow style
!--------------------------------------------------------------------------------------------------
logical function isFlowList(line)

  character(len=*), intent(in) :: line

  isFlowList = index(adjustl(line),'[') == 1

end function isFlowList


!--------------------------------------------------------------------------------------------------
! @brief skip empty lines
! @details update start position in the block by skipping empty lines if present.
!--------------------------------------------------------------------------------------------------
subroutine skip_empty_lines(blck,s_blck)

  character(len=*), intent(in)     :: blck
  integer,          intent(inout)  :: s_blck

  logical :: empty

  empty = .true.
  do while(empty .and. len_trim(blck(s_blck:)) /= 0)
    empty = len_trim(IO_rmComment(blck(s_blck:s_blck + index(blck(s_blck:),IO_EOL) - 2))) == 0
    if(empty) s_blck = s_blck + index(blck(s_blck:),IO_EOL)
  enddo

end subroutine skip_empty_lines
 

!--------------------------------------------------------------------------------------------------
! @brief skip file header
! @details update start position in the block by skipping file header if present.
!--------------------------------------------------------------------------------------------------
subroutine skip_file_header(blck,s_blck)

  character(len=*), intent(in)     :: blck
  integer,          intent(inout)  :: s_blck

  character(len=:), allocatable    :: line

  line = IO_rmComment(blck(s_blck:s_blck + index(blck(s_blck:),IO_EOL) - 2))
  if(index(adjustl(line),'%YAML') == 1) then
    s_blck = s_blck + index(blck(s_blck:),IO_EOL)
    call skip_empty_lines(blck,s_blck)
    if(trim(IO_rmComment(blck(s_blck:s_blck + index(blck(s_blck:),IO_EOL) - 2))) == '---') then
      s_blck = s_blck + index(blck(s_blck:),IO_EOL)
    else
      call IO_error(708,ext_msg = line)
    endif
  endif
 
end subroutine skip_file_header


!--------------------------------------------------------------------------------------------------
!> @brief check if a line in flow YAML starts and ends in the same line
!--------------------------------------------------------------------------------------------------
logical function flow_is_closed(str,e_char)

  character(len=*), intent(in) :: str
  character,        intent(in) :: e_char                                                            !< end of list/dict  ( '}' or ']')
  integer                      :: N_sq, &                                                           !< number of open square brackets
                                  N_cu, &                                                           !< number of open curly brackets
                                  i
  character(len=:), allocatable:: line

  flow_is_closed = .false.
  N_sq = 0
  N_cu = 0
  if(e_char == ']') line = str(index(str(:),'[')+1:)
  if(e_char == '}') line = str(index(str(:),'{')+1:)

  do i = 1, len_trim(line)
    flow_is_closed = (N_sq==0 .and. N_cu==0 .and. scan(line(i:i),e_char) == 1)
    N_sq = N_sq + merge(1,0,line(i:i) == '[')
    N_cu = N_cu + merge(1,0,line(i:i) == '{')
    N_sq = N_sq - merge(1,0,line(i:i) == ']')
    N_cu = N_cu - merge(1,0,line(i:i) == '}')
  enddo

end function flow_is_closed


!--------------------------------------------------------------------------------------------------
!> @brief return the flow YAML line without line break
!--------------------------------------------------------------------------------------------------
subroutine remove_line_break(blck,s_blck,e_char,flow_line)

  character(len=*), intent(in)               :: blck                                                !< YAML in mixed style
  integer,          intent(inout)            :: s_blck
  character,        intent(in)               :: e_char                                              !< end of list/dict  ( '}' or ']')
  character(len=:), allocatable, intent(out) :: flow_line
  logical :: line_end

  line_end =.false.
  flow_line = ''

  do while(.not.line_end)
    flow_line = flow_line//IO_rmComment(blck(s_blck:s_blck + index(blck(s_blck:),IO_EOL) - 2))//' '
    line_end  = flow_is_closed(flow_line,e_char)
    s_blck    = s_blck + index(blck(s_blck:),IO_EOL)
  enddo

end subroutine remove_line_break


!--------------------------------------------------------------------------------------------------
! @brief reads a line of YAML block which is already in flow style
! @details Dicts should be enlcosed within '{}' for it to be consistent with DAMASK YAML parser
!--------------------------------------------------------------------------------------------------
recursive subroutine line_isFlow(flow,s_flow,line)

  character(len=*), intent(inout) :: flow                                                           !< YAML in flow style only
  integer,          intent(inout) :: s_flow                                                         !< start position in flow
  character(len=*), intent(in)    :: line

  integer :: &
    s, &
    list_chunk, &
    dict_chunk

  if(index(adjustl(line),'[') == 1) then
    s = index(line,'[')
    flow(s_flow:s_flow) = '['
    s_flow = s_flow +1
    do while(s < len_trim(line))
      list_chunk = s + find_end(line(s+1:),']')
      if(iskeyValue(line(s+1:list_chunk-1))) then
        flow(s_flow:s_flow) = '{'
        s_flow = s_flow +1
        call keyValue_toFlow(flow,s_flow,line(s+1:list_chunk-1))
        flow(s_flow:s_flow) = '}'
        s_flow = s_flow +1
      elseif(isFlow(line(s+1:list_chunk-1))) then
        call line_isFlow(flow,s_flow,line(s+1:list_chunk-1))
      else
        call line_toFlow(flow,s_flow,line(s+1:list_chunk-1))
      endif
      flow(s_flow:s_flow+1) = ', '
      s_flow = s_flow +2
      s = s + find_end(line(s+1:),']')
    enddo
    s_flow = s_flow - 1
    if (flow(s_flow-1:s_flow-1) == ',') s_flow = s_flow - 1
    flow(s_flow:s_flow) = ']'
    s_flow = s_flow+1

  elseif(index(adjustl(line),'{') == 1) then
    s = index(line,'{')
    flow(s_flow:s_flow) = '{'
    s_flow = s_flow +1
    do while(s < len_trim(line))
      dict_chunk = s + find_end(line(s+1:),'}')
      if( .not. iskeyValue(line(s+1:dict_chunk-1))) call IO_error(705,ext_msg=line)
      call keyValue_toFlow(flow,s_flow,line(s+1:dict_chunk-1))
      flow(s_flow:s_flow+1) = ', '
      s_flow = s_flow +2
      s = s + find_end(line(s+1:),'}')
    enddo
    s_flow = s_flow -1
    if(flow(s_flow-1:s_flow-1) == ',') s_flow = s_flow -1
    flow(s_flow:s_flow) = '}'
    s_flow = s_flow +1
  else
    call line_toFlow(flow,s_flow,line)
  endif

end subroutine line_isFlow


!-------------------------------------------------------------------------------------------------
! @brief reads a line of YAML block of type <key>: <value> and places it in the YAML flow style structure
! @details Makes sure that the <value> is consistent with the input required in DAMASK YAML parser
!-------------------------------------------------------------------------------------------------
recursive subroutine keyValue_toFlow(flow,s_flow,line)

  character(len=*), intent(inout) :: flow                                                           !< YAML in flow style only
  integer,          intent(inout) :: s_flow                                                         !< start position in flow
  character(len=*), intent(in)    :: line

  character(len=:), allocatable   :: line_asStandard                                                ! standard form of <key>: <value>
  integer :: &
    d_flow, &
    col_pos, &
    offset_value

  col_pos = index(line,':')
  if(isFlow(line(col_pos+1:))) then
    d_flow = len_trim(adjustl(line(:col_pos)))
    flow(s_flow:s_flow+d_flow+1) = trim(adjustl(line(:col_pos)))//' '
    s_flow = s_flow + d_flow+1
    call line_isFlow(flow,s_flow,line(col_pos+1:))
  else
    offset_value = indentDepth(line(col_pos+2:))
    line_asStandard = line(:col_pos+1)//line(col_pos+2+offset_value:)
    call line_toFlow(flow,s_flow,line_asStandard)
  endif

end subroutine keyValue_toFlow


!-------------------------------------------------------------------------------------------------
! @brief reads a line of YAML block and places it in the YAML flow style structure
!-------------------------------------------------------------------------------------------------
subroutine line_toFlow(flow,s_flow,line)

  character(len=*), intent(inout) :: flow                                                           !< YAML in flow style only
  integer,          intent(inout) :: s_flow                                                         !< start position in flow
  character(len=*), intent(in)    :: line

  integer :: &
    d_flow

  d_flow = len_trim(adjustl(line))
  flow(s_flow:s_flow+d_flow) = trim(adjustl(line))
  s_flow = s_flow + d_flow

end subroutine line_toFlow


!-------------------------------------------------------------------------------------------------
! @brief convert a yaml list in block style to a yaml list in flow style
! @details enters the function when encountered with the list indicator '- '
! reads each scalar list item and separates each other with a ','
! If list item is non scalar, it stores the offset for that list item block
! Increase in the indentation level or when list item is not scalar -> 'decide' function is called.
! decrease in indentation level indicates the end of an indentation block
!-------------------------------------------------------------------------------------------------
recursive subroutine lst(blck,flow,s_blck,s_flow,offset)

  character(len=*), intent(in)    :: blck                                                           !< YAML in mixed style
  character(len=*), intent(inout) :: flow                                                           !< YAML in flow style only
  integer,          intent(inout) :: s_blck, &                                                      !< start position in blck
                                     s_flow, &                                                      !< start position in flow
                                     offset                                                         !< stores leading '- ' in nested lists
  character(len=:), allocatable :: line,flow_line
  integer :: e_blck,indent

  indent = indentDepth(blck(s_blck:),offset)
  do while (s_blck <= len_trim(blck))
    e_blck = s_blck + index(blck(s_blck:),IO_EOL) - 2
    line = IO_rmComment(blck(s_blck:e_blck))
    if(trim(line) == '---' .or. trim(line) == '...') then
      exit
    elseif (len_trim(line) == 0) then
      s_blck = e_blck + 2                                                                           ! forward to next line
      cycle
    elseif(indentDepth(line,offset) > indent) then
      call decide(blck,flow,s_blck,s_flow,offset)
      offset = 0
      flow(s_flow:s_flow+1) = ', '
      s_flow = s_flow + 2
    elseif(indentDepth(line,offset) < indent .or. .not. isListItem(line)) then
      offset = 0
      exit                                                                                          ! job done (lower level)
    else
      if(trim(adjustl(line)) == '-') then                                                           ! list item in next line
        s_blck = e_blck + 2
        call skip_empty_lines(blck,s_blck)
        e_blck = s_blck + index(blck(s_blck:),IO_EOL) - 2
        line = IO_rmComment(blck(s_blck:e_blck))
        if(trim(line) == '---') call IO_error(707,ext_msg=line)
        if(indentDepth(line) < indent .or. indentDepth(line) == indent) &
          call IO_error(701,ext_msg=line)

        if(isScalar(line)) then
          call line_toFlow(flow,s_flow,line)
          s_blck = e_blck +2
          offset = 0
        elseif(isFlow(line)) then
          if(isFlowList(line)) then
            call remove_line_break(blck,s_blck,']',flow_line)
          else
            call remove_line_break(blck,s_blck,'}',flow_line)
          endif
          call line_isFlow(flow,s_flow,flow_line)
          offset = 0
        endif
      else                                                                                          ! list item in the same line
        line = line(indentDepth(line)+3:)
        if(isScalar(line)) then
          call line_toFlow(flow,s_flow,line)
          s_blck = e_blck +2
          offset = 0
        elseif(isFlow(line)) then
          s_blck = s_blck + index(blck(s_blck:),'-')
          if(isFlowList(line)) then
            call remove_line_break(blck,s_blck,']',flow_line)
          else
            call remove_line_break(blck,s_blck,'}',flow_line)
          endif
          call line_isFlow(flow,s_flow,flow_line)
          offset = 0
        else                                                                                        ! non scalar list item
          offset = offset + indentDepth(blck(s_blck:))+1                                            ! offset in spaces to be ignored
          s_blck = s_blck + index(blck(s_blck:e_blck),'-')                                          ! s_blck after '-' symbol
         endif
      end if
    end if

    if(isScalar(line) .or. isFlow(line)) then
      flow(s_flow:s_flow+1) = ', '
      s_flow = s_flow + 2
    endif

  end do

  s_flow = s_flow - 1
  if (flow(s_flow-1:s_flow-1) == ',') s_flow = s_flow - 1

end subroutine lst


!--------------------------------------------------------------------------------------------------
! @brief convert a yaml dict in block style to a yaml dict in flow style
! @details enters the function when encountered with the dictionary indicator ':'
! parses each line in the block and compares indentation of a line with the preceding line
! upon increase in indentation level -> 'decide' function decides if the line is a list or dict
! decrease in indentation indicates the end of an indentation block
!--------------------------------------------------------------------------------------------------
recursive subroutine dct(blck,flow,s_blck,s_flow,offset)

  character(len=*), intent(in)    :: blck                                                           !< YAML in mixed style
  character(len=*), intent(inout) :: flow                                                           !< YAML in flow style only
  integer,          intent(inout) :: s_blck, &                                                      !< start position in blck
                                     s_flow, &                                                      !< start position in flow
                                     offset

  character(len=:), allocatable :: line,flow_line
  integer :: e_blck,indent,col_pos
  logical :: previous_isKey

  previous_isKey = .false.

  indent = indentDepth(blck(s_blck:),offset)

  do while (s_blck <= len_trim(blck))
    e_blck = s_blck + index(blck(s_blck:),IO_EOL) - 2
    line = IO_rmComment(blck(s_blck:e_blck))
    if(trim(line) == '---' .or. trim(line) == '...') then
      exit
    elseif (len_trim(line) == 0) then
      s_blck = e_blck + 2                                                                           ! forward to next line
      cycle
    elseif(indentDepth(line,offset) < indent) then
      if(isScalar(line) .or. isFlow(line) .and. previous_isKey) &
        call IO_error(701,ext_msg=line)
      offset = 0
      exit                                                                                          ! job done (lower level)
    elseif(indentDepth(line,offset) > indent .or. isListItem(line)) then
      offset = 0
      call decide(blck,flow,s_blck,s_flow,offset)
    else
      if(isScalar(line)) call IO_error(701,ext_msg=line)
      if(isFlow(line))   call IO_error(702,ext_msg=line)

      line = line(indentDepth(line)+1:)
      if(previous_isKey) then
        flow(s_flow-1:s_flow) = ', '
        s_flow = s_flow + 1
      endif

      if(isKeyValue(line)) then
        col_pos = index(line,':')
        if(isFlow(line(col_pos+1:))) then
          if(isFlowList(line(col_pos+1:))) then
            call remove_line_break(blck,s_blck,']',flow_line)
          else
            call remove_line_break(blck,s_blck,'}',flow_line)
          endif
          call keyValue_toFlow(flow,s_flow,flow_line)
        else
          call keyValue_toFlow(flow,s_flow,line)
          s_blck = e_blck + 2
        endif
      else
        call line_toFlow(flow,s_flow,line)
        s_blck = e_blck + 2
      endif
    end if

    if(isScalar(line) .or. isKeyValue(line)) then
      flow(s_flow:s_flow) = ','
      s_flow = s_flow + 1
      previous_isKey = .false.
    else
      previous_isKey = .true.
    endif

    flow(s_flow:s_flow) = ' '
    s_flow = s_flow + 1
    offset = 0
  end do

  s_flow = s_flow - 1
  if (flow(s_flow-1:s_flow-1) == ',') s_flow = s_flow - 1

end subroutine dct


!--------------------------------------------------------------------------------------------------
! @brief decide whether next block is list or dict
!--------------------------------------------------------------------------------------------------
recursive subroutine decide(blck,flow,s_blck,s_flow,offset)

  character(len=*), intent(in)    :: blck                                                           !< YAML in mixed style
  character(len=*), intent(inout) :: flow                                                           !< YAML in flow style only
  integer,          intent(inout) :: s_blck, &                                                      !< start position in blck
                                     s_flow, &                                                      !< start position in flow
                                     offset
  integer :: e_blck
  character(len=:), allocatable :: line,flow_line

  if(s_blck <= len(blck)) then
    call skip_empty_lines(blck,s_blck)
    e_blck = s_blck + index(blck(s_blck:),IO_EOL) - 2
    line = IO_rmComment(blck(s_blck:e_blck))
    if(trim(line) == '---' .or. trim(line) == '...') then
      continue                                                                                      ! end parsing at this point but not stop the simulation
    elseif(len_trim(line) == 0) then
      s_blck = e_blck +2
      call decide(blck,flow,s_blck,s_flow,offset)
    elseif    (isListItem(line)) then
      flow(s_flow:s_flow) = '['
      s_flow = s_flow + 1
      call lst(blck,flow,s_blck,s_flow,offset)
      flow(s_flow:s_flow) = ']'
      s_flow = s_flow + 1
    elseif(isKey(line) .or. isKeyValue(line)) then
      flow(s_flow:s_flow) = '{'
      s_flow = s_flow + 1
      call dct(blck,flow,s_blck,s_flow,offset)
      flow(s_flow:s_flow) = '}'
      s_flow = s_flow + 1
    elseif(isFlow(line)) then
      if(isFlowList(line)) then
        call remove_line_break(blck,s_blck,']',flow_line)
      else
        call remove_line_break(blck,s_blck,'}',flow_line)
      endif
      call line_isFlow(flow,s_flow,line)
    else
      line = line(indentDepth(line)+1:)
      call line_toFlow(flow,s_flow,line)
      s_blck = e_blck +2
    endif
  endif

end subroutine


!--------------------------------------------------------------------------------------------------
! @brief convert all block style YAML parts to flow style
!--------------------------------------------------------------------------------------------------
function to_flow(blck)

  character(len=:), allocatable :: to_flow
  character(len=*), intent(in)  :: blck                                                             !< YAML mixed style

  character(len=:), allocatable :: line
  integer                       :: s_blck, &                                                        !< start position in blck
                                   s_flow, &                                                        !< start position in flow
                                   offset, &                                                        !< counts leading '- ' in nested lists
                                   end_line
 
  allocate(character(len=len(blck)*2)::to_flow)
  s_flow = 1
  s_blck = 1
  offset = 0

  if(len_trim(blck) /= 0) then
    call skip_empty_lines(blck,s_blck)
    call skip_file_header(blck,s_blck)
    line = IO_rmComment(blck(s_blck:s_blck + index(blck(s_blck:),IO_EOL) - 2))
    if(trim(line) == '---') s_blck = s_blck + index(blck(s_blck:),IO_EOL)
    call decide(blck,to_flow,s_blck,s_flow,offset)
  endif
  line = IO_rmComment(blck(s_blck:s_blck+index(blck(s_blck:),IO_EOL)-2))
  if(trim(line)== '---') call IO_warning(709,ext_msg=line)
  to_flow = trim(to_flow(:s_flow-1))
  end_line = index(to_flow,IO_EOL)
  if(end_line > 0) to_flow = to_flow(:end_line-1)

end function to_flow


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some YAML functions.
!--------------------------------------------------------------------------------------------------
subroutine selfTest

  if (indentDepth(' a') /= 1)     error stop 'indentDepth'
  if (indentDepth('a')  /= 0)     error stop 'indentDepth'
  if (indentDepth('x ') /= 0)     error stop 'indentDepth'

  if (      isFlow(' a'))         error stop 'isFLow'
  if (.not. isFlow('{'))          error stop 'isFlow'
  if (.not. isFlow(' ['))         error stop 'isFlow'

  if (      isListItem(' a'))     error stop 'isListItem'
  if (      isListItem(' -b'))    error stop 'isListItem'
  if (.not. isListItem('- a '))   error stop 'isListItem'
  if (.not. isListItem('- -a '))  error stop 'isListItem'

  if (      isKeyValue(' a'))     error stop 'isKeyValue'
  if (      isKeyValue(' a: '))   error stop 'isKeyValue'
  if (.not. isKeyValue(' a: b'))  error stop 'isKeyValue'

  if (      isKey(' a'))          error stop 'isKey'
  if (      isKey('{a:b}'))       error stop 'isKey'
  if (      isKey(' a:b'))        error stop 'isKey'
  if (.not. isKey(' a: '))        error stop 'isKey'
  if (.not. isKey(' a:'))         error stop 'isKey'
  if (.not. isKey(' a: #'))       error stop 'isKey'

  if(       isScalar('a:  '))     error stop 'isScalar'
  if(       isScalar('a: b'))     error stop 'isScalar'
  if(       isScalar('{a:b}'))    error stop 'isScalar'
  if(       isScalar('- a:'))     error stop 'isScalar'
  if(.not.  isScalar('   a'))     error stop 'isScalar'

  basic_list: block
  character(len=*), parameter :: block_list = &
    " - Casablanca"//IO_EOL//&
    " - North by Northwest"//IO_EOL
  character(len=*), parameter :: block_list_newline = &
    " -"//IO_EOL//&
    "   Casablanca"//IO_EOL//&
    " -"//IO_EOL//&
    "   North by Northwest"//IO_EOL
  character(len=*), parameter :: flow_list = &
    "[Casablanca, North by Northwest]"

  if (.not. to_flow(block_list)         == flow_list) error stop 'to_flow'
  if (.not. to_flow(block_list_newline) == flow_list) error stop 'to_flow'
  end block basic_list

  basic_dict: block
  character(len=*), parameter :: block_dict = &
    " aa: Casablanca"//IO_EOL//&
    " bb: North by Northwest"//IO_EOL
  character(len=*), parameter :: block_dict_newline = &
    " aa:"//IO_EOL//&
    "   Casablanca"//IO_EOL//&
    " bb:"//IO_EOL//&
    "   North by Northwest"//IO_EOL
  character(len=*), parameter :: flow_dict = &
    "{aa: Casablanca, bb: North by Northwest}"

  if (.not. to_flow(block_dict)         == flow_dict) error stop 'to_flow'
  if (.not. to_flow(block_dict_newline) == flow_dict) error stop 'to_flow'
  end block basic_dict

  only_flow: block
  character(len=*), parameter :: flow_dict = &
    " {a: [b,c: {d: e}, f: g, e]}"//IO_EOL
  character(len=*), parameter :: flow_list = &
    " [a,b: c,    d,e: {f: g}]"//IO_EOL
  character(len=*), parameter :: flow_1 = &
    "{a: [b, {c: {d: e}}, {f: g}, e]}"
  character(len=*), parameter :: flow_2 = &
    "[a, {b: c}, d, {e: {f: g}}]"

  if (.not. to_flow(flow_dict)        == flow_1) error stop 'to_flow'
  if (.not. to_flow(flow_list)        == flow_2) error stop 'to_flow'
  end block only_flow

  basic_flow: block
  character(len=*), parameter :: flow_braces = &
    " source: [{param: 1}, {param: 2}, {param: 3}, {param: 4}]"//IO_EOL
  character(len=*), parameter :: flow_mixed_braces = &
    " source: [param: 1, {param: 2}, param: 3, {param: 4}]"//IO_EOL
  character(len=*), parameter :: flow = &
    "{source: [{param: 1}, {param: 2}, {param: 3}, {param: 4}]}"

  if (.not. to_flow(flow_braces)        == flow) error stop 'to_flow'
  if (.not. to_flow(flow_mixed_braces)  == flow) error stop 'to_flow'
  end block basic_flow

  multi_line_flow1: block
  character(len=*), parameter :: flow_multi = &
    "%YAML 1.1"//IO_EOL//&
    "---"//IO_EOL//&
    "a: [b,"//IO_EOL//&
    "c: "//IO_EOL//&
    "d, e]"//IO_EOL

  character(len=*), parameter :: flow = &
    "{a: [b, {c: d}, e]}"
  
  if( .not. to_flow(flow_multi)        == flow) error stop 'to_flow'
  end block multi_line_flow1

  multi_line_flow2: block
  character(len=*), parameter :: flow_multi = &
    "%YAML 1.1"//IO_EOL//&
    "---"//IO_EOL//&
    "-"//IO_EOL//&
    " a: {b:"//IO_EOL//&
    "[c,"//IO_EOL//&
    "d"//IO_EOL//&
    "e, f]}"//IO_EOL
 
  character(len=*), parameter :: flow = &
    "[{a: {b: [c, d e, f]}}]"

  if( .not. to_flow(flow_multi)        == flow) error stop 'to_flow'
  end block multi_line_flow2

  basic_mixed: block
  character(len=*), parameter :: block_flow = &
    "%YAML 1.1"//IO_EOL//&
    " "//IO_EOL//&
    " "//IO_EOL//&
    "---"//IO_EOL//&
    " aa:"//IO_EOL//&
    " - "//IO_EOL//&
    " "//IO_EOL//&
    " "//IO_EOL//&
    "                 param_1: [a:                   b, c, {d: {e: [f: g, h]}}]"//IO_EOL//&
    " - c: d"//IO_EOL//&
    " bb:"//IO_EOL//&
    " "//IO_EOL//&
    "  - "//IO_EOL//&
    "   {param_1: [{a: b}, c, {d: {e: [{f: g}, h]}}]}"//IO_EOL//&
    "..."//IO_EOL
  character(len=*), parameter :: mixed_flow = &
    "{aa: [{param_1: [{a: b}, c, {d: {e: [{f: g}, h]}}]}, {c: d}], bb: [{param_1: [{a: b}, c, {d: {e: [{f: g}, h]}}]}]}"

  if(.not. to_flow(block_flow) == mixed_flow)    error stop 'to_flow'
  end block basic_mixed

end subroutine selfTest

end module YAML_parse
