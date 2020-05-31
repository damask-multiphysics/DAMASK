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
 
  public :: YAML_init
  public :: parse_flow,to_flow

contains

!--------------------------------------------------------------------------------------------------
!> @brief do sanity checks
!--------------------------------------------------------------------------------------------------
subroutine YAML_init

  call selfTest

end subroutine YAML_init


!--------------------------------------------------------------------------------------------------
!> @brief reads the flow style string and stores it in the form of dictionaries, lists and scalars.
!> @details  A node type pointer can either point to a dictionary, list or scalar type entities.
!--------------------------------------------------------------------------------------------------
recursive function parse_flow(flow_string) result(node)

  character(len=*), intent(inout) :: flow_string
  class (tNode), pointer          :: node

  class (tNode),    pointer       :: myVal
  character(len=pStringLen)       :: key

  integer                         :: e, &                                                           !> end position of dictionary or list
                                     s, &                                                           !> start position of dictionary or list
                                     d                                                              !> position of key: value separator (':')

  flow_string = trim(adjustl(flow_string(:)))
  if (flow_string(1:1) == '{') then                                                                 ! start of a dictionary
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

  character(len=*), intent(in) :: str
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

  isListItem = index(adjustl(line),'-') == 1

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
    isKey = IO_rmComment(line(len(IO_rmComment(line)):len(IO_rmComment(line)))) == ':' & 
                                                                .and. .not. isFlow(line)
  endif

end function isKey


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
  character(len=pStringLen) :: line
  integer :: e_blck,indent

  indent = indentDepth(blck(s_blck:),offset)
  do while (s_blck <= len_trim(blck))
    e_blck = s_blck + index(blck(s_blck:),IO_EOL) - 2
    line = IO_rmComment(blck(s_blck:e_blck))
    if (len_trim(line) == 0) then
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
        e_blck = e_blck + index(blck(e_blck+2:),IO_EOL)
        line = IO_rmComment(blck(s_blck:e_blck))
        if(indentDepth(line) < indent .or. indentDepth(line) == indent) &
          call IO_error(701,ext_msg=line)

        if(isScalar(line)) then
          call line_toFlow(flow,s_flow,line)
          s_blck = e_blck +2
          offset = 0
        elseif(isFlow(line)) then
          call line_isFlow(flow,s_flow,line)
          s_blck = e_blck +2
          offset = 0
        endif
      else                                                                                          ! list item in the same line
        if(line(indentDepth(line)+2:indentDepth(line)+2) /= ' ') &
          call IO_error(703,ext_msg=line)
        line = line(indentDepth(line)+3:)
        if(isScalar(line)) then
          call line_toFlow(flow,s_flow,line)
          s_blck = e_blck +2
          offset = 0
        elseif(isFlow(line)) then
          call line_isFlow(flow,s_flow,line)
          s_blck = e_blck +2
          offset = 0
        else                                                                                        ! non scalar list item
          offset = offset + indentDepth(blck(s_blck:))+1                                            ! offset in spaces to be ignored
          s_blck = s_blck + index(blck(s_blck:e_blck),'-')                                          ! s_blck after '-' symbol
         endif
      end if
    end if

    if(isScalar(line) .or. isFlow(line)) then
      flow(s_flow:s_flow+1) = ', '
      s_flow = s_flow +2 
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

  character(len=pStringLen) :: line
  integer :: e_blck,indent
  logical :: previous_isKey

  previous_isKey = .false.

  indent = indentDepth(blck(s_blck:),offset)

  do while (s_blck <= len_trim(blck))
    e_blck = s_blck + index(blck(s_blck:),IO_EOL) - 2
    line = IO_rmComment(blck(s_blck:e_blck))
    if (len_trim(line) == 0) then
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
        call keyValue_toFlow(flow,s_flow,line)
      else
        call line_toFlow(flow,s_flow,line)
      endif
      
      s_blck = e_blck +2
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
  character(len=pStringLen) :: line

  if(s_blck <= len(blck)) then
    e_blck = s_blck + index(blck(s_blck:),IO_EOL) - 2
    line = IO_rmComment(blck(s_blck:e_blck))

    ! exit here if '---' is found
    if    (isListItem(line)) then
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
      call line_isFlow(flow,s_flow,line)
      s_blck = e_blck +2
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
  integer                       :: s_blck, &                                                        !< start position in blck
                                   s_flow, &                                                        !< start position in flow
                                   offset, &                                                        !< counts leading '- ' in nested lists
                                   end_line
  if(isFlow(blck)) then                                                                       
    to_flow = trim(adjustl(blck))
  else
    allocate(character(len=len(blck)*2)::to_flow)
    ! move forward here (skip empty lines) and remove '----' if found
    s_flow = 1
    s_blck = 1
    offset = 0
    call decide(blck,to_flow,s_blck,s_flow,offset)
    to_flow = trim(to_flow(:s_flow-1))
  endif
    end_line = index(to_flow,new_line(''))
    if(end_line > 0) to_flow = to_flow(:end_line-1) 

end function to_flow


!--------------------------------------------------------------------------------------------------
subroutine selfTest()

  if (indentDepth(' a') /= 1)     call IO_error(0,ext_msg='indentDepth')
  if (indentDepth('a')  /= 0)     call IO_error(0,ext_msg='indentDepth')
  if (indentDepth('x ') /= 0)     call IO_error(0,ext_msg='indentDepth')

  if (      isFlow(' a'))         call IO_error(0,ext_msg='isFLow')
  if (.not. isFlow('{'))          call IO_error(0,ext_msg='isFlow')
  if (.not. isFlow(' ['))         call IO_error(0,ext_msg='isFlow')

  if (      isListItem(' a'))     call IO_error(0,ext_msg='isListItem')
  if (.not. isListItem('- a '))   call IO_error(0,ext_msg='isListItem')
  if (.not. isListItem(' -b'))    call IO_error(0,ext_msg='isListItem')

  if (      isKeyValue(' a'))     call IO_error(0,ext_msg='isKeyValue')
  if (      isKeyValue(' a: '))   call IO_error(0,ext_msg='isKeyValue')
  if (.not. isKeyValue(' a: b'))  call IO_error(0,ext_msg='isKeyValue')

  if (      isKey(' a'))          call IO_error(0,ext_msg='isKey')
  if (      isKey('{a:b}'))       call IO_error(0,ext_msg='isKey') 
  if (      isKey(' a:b'))        call IO_error(0,ext_msg='isKey')
  if (.not. isKey(' a: '))        call IO_error(0,ext_msg='isKey')
  if (.not. isKey(' a:'))         call IO_error(0,ext_msg='isKey')
  if (.not. isKey(' a: #'))       call IO_error(0,ext_msg='isKey')

  if(       isScalar('a:  '))     call IO_error(0,ext_msg='isScalar')
  if(       isScalar('a: b'))     call IO_error(0,ext_msg='isScalar')
  if(       isScalar('{a:b}'))    call IO_error(0,ext_msg='isScalar')
  if(       isScalar('- a:'))     call IO_error(0,ext_msg='isScalar')
  if(.not.  isScalar('   a'))     call IO_error(0,ext_msg='isScalar')

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

  if (.not. to_flow(block_list)         == flow_list)     call IO_error(0,ext_msg='to_flow')
  if (.not. to_flow(block_list_newline) == flow_list)     call IO_error(0,ext_msg='to_flow')
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

  if (.not. to_flow(block_dict)         == flow_dict)     call IO_error(0,ext_msg='to_flow')
  if (.not. to_flow(block_dict_newline) == flow_dict)     call IO_error(0,ext_msg='to_flow')
  end block basic_dict
  
  basic_flow: block
  character(len=*), parameter :: flow_braces = &
    " source: [{param: 1}, {param: 2}, {param: 3}, {param: 4}]"//IO_EOL
  character(len=*), parameter :: flow_mixed_braces = &
    " source: [param: 1, {param: 2}, param: 3, {param: 4}]"//IO_EOL
  character(len=*), parameter :: flow = &
    "{source: [{param: 1}, {param: 2}, {param: 3}, {param: 4}]}"
 
  if (.not. to_flow(flow_braces)        == flow)           call IO_error(0,ext_msg='to_flow')
  if (.not. to_flow(flow_mixed_braces)  == flow)           call IO_error(0,ext_msg='to_flow')
  end block basic_flow

  basic_mixed: block
  character(len=*), parameter :: block_flow = &
    " aa:"//IO_EOL//&
    " - "//IO_EOL//&
    "  param_1: [a:                   b, c, {d: {e: [f: g, h]}}]"//IO_EOL//&
    " - c: d"//IO_EOL//&
    " bb:"//IO_EOL//&
    "  - {param_1: [{a: b}, c, {d: {e: [{f: g}, h]}}]}"//IO_EOL
  character(len=*), parameter :: mixed_flow = &
    "{aa: [{param_1: [{a: b}, c, {d: {e: [{f: g}, h]}}]}, {c: d}], bb: [{param_1: [{a: b}, c, {d: {e: [{f: g}, h]}}]}]}"
 
  if(.not. to_flow(block_flow)           == mixed_flow)     call IO_error(0,ext_msg='to_flow')
  end block basic_mixed

end subroutine selfTest

end module YAML_parse
